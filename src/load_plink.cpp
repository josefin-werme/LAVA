#include "load_plink.h"



list load_plink(string prefix, int total_snps, int total_indiv, integers snp_index, list parameters) {return load(prefix, total_snps, total_indiv, snp_index, parameters, Genotypes);}
list compute_ld(string prefix, int total_snps, int total_indiv, integers snp_index, list parameters) {return load(prefix, total_snps, total_indiv, snp_index, parameters, Correlations);}

list load(string prefix, int total_snps, int total_indiv, integers snp_index, list parameters, OutputMode mode) {
	unsigned int mac = (unsigned int) max(Utils::list_value(parameters, "mac", 1.0), 1.0);
	double maf = min(max(Utils::list_value(parameters, "maf", 0.0), 0.0), 0.5);
	double missing = min(max(Utils::list_value(parameters, "missing", 0.05), 0.0), 0.25);
	
	BinaryPlink plink_input(prefix, total_snps > 0 ? total_snps : 0, total_indiv > 0 ? total_indiv : 0);
	plink_input.set_thresholds(mac, maf, missing).read_data(snp_index);
	
	integers include = plink_input.get_included(); writable::list output({"include"_nm = include});
	if (include.size() > 0) {
		if (mode == Genotypes) output.push_back("data"_nm = plink_input.load_genotypes());
		else if (mode == Correlations) output.push_back("ld"_nm = plink_input.compute_ld());
	} 
	return output;
}



doubles_matrix<> compute_info(string prefix, int total_snps, int total_indiv, integers snp_index) {
	BinaryPlink plink_input(prefix, total_snps > 0 ? total_snps : 0, total_indiv > 0 ? total_indiv : 0);
	return plink_input.compute_info(snp_index);
}





BinaryPlink::BinaryPlink(const string& file_prefix, unsigned int total_snps, unsigned int total_indiv) : prefix(file_prefix), no_snps(total_snps), no_individuals(total_indiv) {this->prep_bed();}


// data is stored in blocks of four genotypes in one byte (three bytes offset at start of file), ie. two bits per genotype
// - bit coding: hom1 (00), miss (01), het (10), hom2 (11) (NB: this is when read in stored order)
// - individuals are stored in reverse order within each block; if not a multiple of four, final block for a SNP is padded with 00's (in front, ie. dummy individuals)
// function sets up indexing array to map byte code of block of four genotypes onto corresponding sequence of four unsigned char values
// - coding is: hom1 (0), het (1), hom2 (2), miss (3); order of individuals is reversed to match natural order in .fam file
void BinaryPlink::prep_bed() {
	string filename = Utils::check_file(this->prefix + ".bed");
	
	this->line_size = ceil(this->no_individuals/4.0);
	this->bed_file.open(filename.c_str(), ios::in|ios::binary|ios::ate);
	unsigned long long bed_size = this->bed_file.tellg(), exp_bed_size = (unsigned long long) this->line_size * this->no_snps + 3; 
	
	char buffer[3]; this->bed_file.seekg(0, ios::beg); this->bed_file.read(buffer, 3);
	if (this->bed_file.fail() || ((unsigned short) buffer[0] != 108 || (unsigned short) buffer[1] != 27)) error("file is not a valid .bed file");
	if ((unsigned short) buffer[2] != 1) {if (buffer[2] == 0) error("file is in individual-major format"); else error("file-format specifier is not valid");}
	if (bed_size != exp_bed_size) error("size of .bed file is not consistent with number of SNPs and individuals in .bim and .fam files");

	memset(this->geno_index, 0, 256*4);
	for (int i = 0; i < 256; i++) {
		unsigned char* target = reinterpret_cast<unsigned char*>(this->geno_index + i);
		for (unsigned char j = 0; j < 4; j++) {
			unsigned char value = (unsigned char) i >> (2*j) & 0b0000'0011;
			target[j] = (value != 0b0000'0001) ? ( (value >= 0b0000'00010) ? value - 1 : value ) : 3; 
		}
	}
}


bool BinaryPlink::check_snp(unsigned int* counts) {
	if (counts[Missing] / (double) this->no_individuals <= this->missing_thresh) {
		unsigned int nobs = this->no_individuals - counts[Missing], mac = counts[HetGen] + 2*counts[Hom2];
		if (mac > nobs) mac = 2*nobs - mac;
		return mac >= this->mac_thresh && mac / (2.0 * nobs) >= this->maf_thresh;
	} else return false;
}


BinaryPlink& BinaryPlink::read_data(integers snp_index) {
	if (snp_index.size() == 0) error("SNP index is empty");
	
	Array<unsigned int> step_index = Utils::process_index(snp_index, this->no_snps, true);
	unsigned int no_used = step_index.size(); Array<char> raw_buffer(this->line_size); 
	this->read_buffer.set_size(4*this->line_size, no_used); this->count_buffer.set_size(4, no_used).zero_fill();
	this->snp_include.set_size(no_used); unsigned int* include = this->snp_include.data();
	
	unsigned long long file_offset = 3; unsigned int total_missing = 0;
	for (unsigned int i_snp = 0; i_snp < no_used; i_snp++) {
		file_offset += ((unsigned long long) this->line_size) * step_index[i_snp];
		this->bed_file.seekg(file_offset).read(raw_buffer.data(), this->line_size);
		
		char* read = raw_buffer.data(); unsigned int* counts = this->count_buffer.column(i_snp);
		for (unsigned char *write = this->read_buffer.column(i_snp), *write_end = this->read_buffer.column(i_snp+1); write < write_end; write += 4) memcpy(write, this->geno_index + (unsigned char) *read++, 4);
		for (unsigned char *value = this->read_buffer.column(i_snp), *value_end = value + this->no_individuals; value < value_end; value++) counts[*value]++;
		
		total_missing += counts[Missing];
		if (this->check_snp(counts)) *include++ = i_snp;
	}
	this->snp_include.set_size(include - this->snp_include.data());
	
	this->missing_buffer.set_size(total_missing); this->missing_index.set_size(no_used+1);
	unsigned int* curr = this->missing_buffer.data(); this->missing_index[0] = curr;
	for (unsigned int i_snp = 0; i_snp < no_used; i_snp++) {
		if (count_buffer(Missing, i_snp) > 0) {
			unsigned char *read_start = this->read_buffer.column(i_snp), *read_end = read_start + this->no_individuals;
			for (unsigned char *read = read_start; read < read_end; read++) {
				if (*read == Missing) {*read = 0; *curr++ = read - read_start;}
			}
		}
		this->missing_index[i_snp + 1] = curr;
	}

	return *this;
}



unsigned int BinaryPlink::check_input() {	
	if (this->read_buffer.is_empty()) error("no data has been read");
	if (this->snp_include.size() == 0) error("no SNPs remaining after filtering");
	return this->snp_include.size();
}


integers_matrix<> BinaryPlink::load_genotypes() {
	int no_used = this->check_input();
	
	writable::integers_matrix<> data(this->no_individuals, no_used);
	for (int i_out = 0; i_out < no_used; i_out++) {
		unsigned int i_in = this->snp_include[i_out];
		auto column = data[i_out]; unsigned char* read = this->read_buffer.column(i_in);
		for (int i_row = 0; i_row < (int) this->no_individuals; i_row++) column[i_row] = *read++;
		if (this->count_buffer(Missing,i_in) > 0) for (unsigned int *curr = this->missing_index[i_in], *end = this->missing_index[i_in+1]; curr < end; curr++) column[*curr] = NA_INTEGER;	
	}

	return data;
}


doubles_matrix<> BinaryPlink::compute_ld() {
	int no_used = this->check_input();
	
	writable::doubles_matrix<> ld(no_used, no_used);
	for (int i_out = 0; i_out < no_used; i_out++) {
		unsigned int i_in = this->snp_include[i_out];
		
		unsigned int *counts = this->count_buffer.column(i_in); Array<uint8_t> missing(this->no_individuals, true); uint8_t* miss_ptr = missing.data();
		for (unsigned int *curr = this->missing_index[i_in], *end = this->missing_index[i_in+1]; curr < end; curr++) miss_ptr[*curr] = 1;
		unsigned char* data = this->read_buffer.column(i_in);
		
		ld(i_out,i_out) = 1;
		for (int i_out_other = i_out + 1; i_out_other < no_used; i_out_other++) {
			unsigned int i_in_other = this->snp_include[i_out_other], *counts_other = this->count_buffer.column(i_in_other);
			unsigned char* data_other = this->read_buffer.column(i_in_other);
			
			unsigned int no_obs = this->no_individuals - counts[Missing] - counts_other[Missing];
			if (counts[Missing] > 0 && counts_other[Missing] > 0) {for (unsigned int *curr = this->missing_index[i_in_other], *end = this->missing_index[i_in_other+1]; curr < end; curr++) if (miss_ptr[*curr]) no_obs++;}

			double corr = 0;
			if (no_obs >= 2) {
				unsigned int s_snp = counts[HetGen] + 2*counts[Hom2], ssq_snp = s_snp + 2*counts[Hom2];
				unsigned int s_other = counts_other[HetGen] + 2*counts_other[Hom2], ssq_other = s_other + 2*counts_other[Hom2];
				unsigned int s_prod = this->product_sum(data, data_other);
		
				if (counts_other[Missing] > 0) this->update_sums(data, this->missing_index[i_in_other], counts_other[Missing], s_snp, ssq_snp);
				if (counts[Missing] > 0) this->update_sums(data_other, this->missing_index[i_in], counts[Missing], s_other, ssq_other);	
		
				double mean_snp = s_snp / double(no_obs), var_snp = (ssq_snp - no_obs*mean_snp*mean_snp) / (no_obs - 1);
				double mean_other = s_other / double(no_obs), var_other = (ssq_other - no_obs*mean_other*mean_other) / (no_obs - 1);
				
				if (var_snp > 0 && var_other > 0) {
					double cov = (s_prod - no_obs*mean_snp*mean_other) / (no_obs - 1);
					corr = min(max(cov / sqrt(var_snp*var_other), -1.0), 1.0);
				}
			}
			ld(i_out, i_out_other) = corr; ld(i_out_other, i_out) = corr;
		}
	}
	
	return ld;
}


// due to reverse storage in PLINK files, bit coding / count index is hom1 (00 / 0), miss (01 / 1), het (10 / 2), hom2 (11 / 3)
doubles_matrix<> BinaryPlink::compute_info(integers snp_index) {
	if (snp_index.size() == 0) error("SNP index is empty");
	
	Array<unsigned int> step_index = Utils::process_index(snp_index, this->no_snps, true); unsigned int no_used = step_index.size(); 
	Array<uint8_t> raw_buffer(this->line_size); uint8_t *buffer_start = raw_buffer.data(), *buffer_end = raw_buffer.end();
	
	writable::doubles_matrix<> info(no_used, 3); //columns: NOBS, MISS, FREQ
		
	unsigned long long file_offset = 3; 
	for (unsigned int i_snp = 0; i_snp < no_used; i_snp++) {
		file_offset += ((unsigned long long) this->line_size) * step_index[i_snp];
		this->bed_file.seekg(file_offset).read(reinterpret_cast<char*>(buffer_start), this->line_size);
		
		unsigned int counts[4] = {0,0,0,0};
		for (uint8_t* data = buffer_start; data < buffer_end; data++) while (*data) {counts[*data & 0b0000'0011]++; *data >>= 2;}
		counts[0] = this->no_individuals - (counts[1]	+ counts[2] + counts[3]);
		
		double nobs = this->no_individuals - counts[1];
		info(i_snp, 0) = nobs;
		info(i_snp, 1) = 1 - nobs / this->no_individuals;
		info(i_snp, 2) = (counts[2] + 2*counts[3]) / (2 * nobs);
	}

	return info;
}




unsigned int BinaryPlink::product_sum(unsigned char* data1, unsigned char* data2) {
	unsigned int sum = 0;	unsigned char* end = data1 + this->no_individuals;
	while (data1 < end) sum += *data1++ * *data2++;
	return sum;
}

void BinaryPlink::update_sums(unsigned char* data, unsigned int* miss_ptr, unsigned int miss_count, unsigned int& sum, unsigned int& sq_sum) {
	for (unsigned int* miss_end = miss_ptr + miss_count; miss_ptr < miss_end; miss_ptr++) {unsigned char value = data[*miss_ptr]; sum -= value; sq_sum -= value*value;}
}








