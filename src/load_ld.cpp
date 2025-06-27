#include "load_ld.h"



doubles_matrix<> load_ld(string prefix, int chromosome, int total_snps, integers snp_index) {
	BinaryLD ld_input(prefix, Utils::check_chromosome(chromosome), total_snps > 0 ? total_snps : 0);
	return ld_input.read_ld(snp_index);
}



BinaryLD::BinaryLD(const string& file_prefix, unsigned int chromosome, unsigned int total_snps) : no_snps(total_snps) {
	this->chromosome = Utils::check_chromosome(chromosome);
	this->prefix = file_prefix + "_chr" + Utils::to_string(this->chromosome);
	
	this->prep_bcor(); this->load_index();
}

void BinaryLD::prep_bcor() {
	string filename = Utils::check_file(this->prefix + ".bcor");
	
	this->bcor_file.open(filename.c_str(), ios::in|ios::binary|ios::ate);
	unsigned long long file_size = this->bcor_file.tellg();
	if (file_size < 12) error("header is incomplete");

	unsigned int check_value = this->header_value(CheckValue); if (check_value != 21775) error("check value " + Utils::to_string(check_value) + " not recognized"); 
	unsigned int header_size = this->header_value(HeaderSize); if (file_size < 4*header_size) error("header is incomplete");		
	
	unsigned int snp_count = this->header_value(NoSNPs); 
	if (snp_count != this->no_snps) error("number of SNPs in input is inconsistent");

	unsigned int value_size = this->header_value(ValueSize), exp_size = 0; this->storage_type = static_cast<StorageType>(this->header_value(ValueType));
	if (this->storage_type == Float) exp_size = sizeof(float);
	else if (this->storage_type == Integer8) exp_size = sizeof(uint8_t);
	else if (this->storage_type == Integer9 || this->storage_type == Integer10 || this->storage_type == Integer12) exp_size = sizeof(uint8_t);	
	else if (this->storage_type == Integer16) exp_size = sizeof(uint16_t);
	else error("incompatible value format");
	if (value_size != exp_size) error("incorrect value size");

	this->read_offset = header_size * 4; this->read_size = value_size;
	unsigned int format_id = this->header_value(FormatID);
	if (format_id == Indexed || format_id == Truncated) {
		this->storage_format = static_cast<StorageFormat>(format_id);

		this->index_offset[0] = this->format_param_long(0);
		if (file_size != index_offset[0] + 4*this->no_snps) error("incorrect file size");
		
		if (this->storage_format == Truncated) {
			this->index_offset[1] = this->format_param_long(2);
			if (index_offset[0] - index_offset[1] != 4*this->no_snps) error("invalid index offsets");
		}
	} else if (format_id == EigenDecomposed) {
		if (this->storage_type != Float) error("invalid storage type for eigen-decomposed LD format"); 
		this->storage_format = EigenDecomposed;
		
		this->index_offset[0] = this->format_param_long(0); this->index_offset[1] = this->format_param_long(2);
		this->no_blocks = this->format_param(4); this->no_entries = this->format_param(5);
		if (file_size != index_offset[1] + 5 * 4 * this->no_entries) error("incorrect file size"); // five values in index: id1, id2, no_values, no_bytes, as_integer/storage_mode (compound)
		if (index_offset[1] - index_offset[0] != 2 * 4 *this->no_blocks) error("invalid index offsets"); // two values in index: no_snps, no_pcs
		if (this->no_blocks == 0 || this->no_entries < this->no_blocks) error("file indices are empty or inconsistent");
	} else error("storage format code " + Utils::to_string(format_id) + " not recognized");
}		


void BinaryLD::load_index() {
	if (this->storage_format == Indexed || this->storage_format == Truncated) {
		this->no_corrs.set_size(this->no_snps); this->get_values(index_offset[0], this->no_corrs.data(), this->no_snps);
	
		Array<uint32_t> line_sizes; uint32_t* size_ptr;
		if (this->storage_format == Truncated) {
			line_sizes.set_size(this->no_snps);
			size_ptr = this->get_values(this->index_offset[1], line_sizes.data(), this->no_snps);
		} else size_ptr = this->no_corrs.data();
		
		this->value_offset.set_size(this->no_snps+1); this->value_offset[0] = 0;
		for (unsigned int i = 1; i <= this->no_snps; i++) this->value_offset[i] = value_offset[i-1] + size_ptr[i-1];	
	} else if (this->storage_format == EigenDecomposed) {
		Buffer<uint32_t> block_info(2, this->no_blocks); this->get_values(index_offset[0], block_info.data(), block_info.size());
		this->block_sizes.reserve(this->no_blocks); 
		this->snp_map.set_size(this->no_snps); unsigned int *snp_curr = this->snp_map.data(), id_offset = 0;
		for (unsigned int i = 0; i < this->no_blocks; i++) {
			unsigned int nsnps = block_info(0,i), npcs = block_info(1,i), *snp_end = snp_curr + nsnps; 
			if (snp_end > this->snp_map.end()) error("inconsistent number of SNPs");
			
			this->block_sizes.push_back(BlockInfo(nsnps, npcs, id_offset)); 
			while (snp_curr < snp_end) *snp_curr++ = i;

			id_offset += nsnps;
		}
		
		unsigned long long total_bytes = 0; Buffer<uint32_t> data_blocks(5, this->no_entries); 
		this->get_values(index_offset[1], data_blocks.data(), data_blocks.size());
		this->block_index.resize(this->no_entries); this->product_info.reserve(this->no_entries);
		this->value_offset.set_size(this->no_entries+1); this->value_offset[0] = 0;
		for (unsigned int i = 0; i < this->no_entries; i++) {
			unsigned int id1 = data_blocks(0,i), id2 = data_blocks(1,i), no_values = data_blocks(2,i), no_bytes = data_blocks(3,i);
			vector<unsigned int>& target = this->block_index[id1];
			if (id2 < id1 || target.size() != (id2 - id1)) error("invalid ordering of data blocks");
			target.push_back(i); 
			
			uint32_t tail = data_blocks(4,i), mode = tail & 0b1111'1111, as_integer = (tail >> 16);
			if (!(as_integer == 0 || as_integer == 1) || mode < Raw || mode > Full) error("invalid block product parameter settings");

			ProductInfo info(no_values, no_bytes, static_cast<ProductMode>(mode), as_integer);
			if (id1 == id2 && (info.as_integer || info.mode != Raw)) error("invalid block product parameters for main block entry");
			
			this->product_info.push_back(info);
			this->value_offset[i+1] = value_offset[i] + no_bytes;
			total_bytes += no_bytes;
		}
		if (total_bytes != (this->index_offset[0] - this->read_offset)) error("file size is inconsistent with block index");
	}
}


BinaryLD::BlockData* BinaryLD::load_block(unsigned int id, vector<unsigned int>& snps) {
	static_assert(sizeof(float) == 4, "code assumes float values have a size of four bytes");	
	BlockData* block = new BlockData(id, this->block_sizes[id], snps); unsigned int entry_offset = this->block_index[id][0];
	ProductInfo& info = this->product_info[entry_offset];

	unsigned long long offset = this->value_offset[entry_offset], no_bytes = this->value_offset[entry_offset + 1] - offset;
	this->read_buffer.set_size(block->no_pcs, block->no_snps + 1);
	if (info.no_values != this->read_buffer.size() || 4*info.no_values != no_bytes) error("inconsistent primary block size");

	this->bcor_file.seekg(offset + this->read_offset);
	this->bcor_file.read(reinterpret_cast<char*>(this->read_buffer.data()), no_bytes);

	block->eigenvectors.set_size(block->no_pcs, block->used.size());
	if (block->used.size() < block->no_snps) {
		for (unsigned int i = 0; i < block->used.size(); i++) {
			unsigned int col = block->used[i];
			memcpy(block->eigenvectors.column(i), this->read_buffer.column(col), block->no_pcs * sizeof(float));
		}
	} else memcpy(block->eigenvectors.data(), this->read_buffer.data(), no_bytes - 4*block->no_pcs);

	block->eigenvalues.set_size(block->no_pcs); 
	memcpy(block->eigenvalues.data(), this->read_buffer.column(block->no_snps), block->no_pcs * sizeof(float));

	return block;
}

Buffer<float>& BinaryLD::load_block_product(unsigned int id1, unsigned int id2) {
	if (id2 <= id1) error("invalid block ID pair");

	unsigned int block_offset = id2 - id1, entry_offset = this->block_index[id1][block_offset]; 
	ProductInfo& info = this->product_info[entry_offset];
	if (block_offset < this->block_index[id1].size()) {
		unsigned long long offset = this->value_offset[entry_offset], no_bytes = this->value_offset[entry_offset + 1] - offset;
		static Array<char> load_buffer; load_buffer.set_size(no_bytes);
		this->bcor_file.seekg(offset + this->read_offset);
		this->bcor_file.read(load_buffer.data(), no_bytes);
		
		unsigned int nrows = this->block_sizes[id1].no_pcs, ncols = this->block_sizes[id2].no_pcs;
		if (info.mode == MQ) {nrows = this->block_sizes[id2].no_pcs; ncols = this->block_sizes[id1].no_snps;} // stored as transpose 
		else if (info.mode == Raw) {nrows = this->block_sizes[id1].no_snps; ncols = this->block_sizes[id2].no_snps;}
		else if (info.mode == QtM) ncols = this->block_sizes[id2].no_snps;
		
		if (info.no_values != nrows * ncols) error("inconsistent block product size");
		this->read_buffer.set_size(nrows, ncols); float* write = this->read_buffer.data();
		
		if (info.as_integer) {
			uint16_t zero_value = numeric_limits<uint16_t>::max(), max_value = zero_value - 1;
			float *range = reinterpret_cast<float*>(load_buffer.data()), min = range[0], max = range[1], scale = (max - min) / max_value;
			
			for (uint16_t *read = reinterpret_cast<uint16_t*>(load_buffer.data() + 8), *end = reinterpret_cast<uint16_t*>(load_buffer.end()); read < end; read++) {
				if (*read == zero_value) for (int count = *(++read); count > 0; count--) *write++ = 0;
				else *write++ = *read * scale + min;
			}
		} else {
			for (float *read = reinterpret_cast<float*>(load_buffer.data()), *end = reinterpret_cast<float*>(load_buffer.end()); read < end; read++) {
				if (*read == 0)	for (int count = *reinterpret_cast<const uint32_t*>(++read); count > 0; count--) *write++ = 0;
				else *write++ = *read;
			}
		}
	} else this->read_buffer.set_size(0);
	
	return this->read_buffer;
}


doubles_matrix<> BinaryLD::read_decomposed(integers snp_index) {
	Array<unsigned int> used_index = Utils::process_index(snp_index, this->no_snps, false); unsigned int no_used = used_index.size();
	writable::doubles_matrix<> ld(no_used,no_used);
	
	if (no_used > 0) {
		for (int i = 0; i < (int) no_used; i++) {auto column = ld[i]; for (int j = 0; j < (int) no_used; j++) column[j] = 0; column[i] = 1;}

		vector<BlockData*> blocks; vector<unsigned int> curr_ids; unsigned int curr_block = this->snp_map[used_index[0]];
		for (unsigned int i = 0; i <= no_used; i++) {
			if (i == no_used || this->snp_map[used_index[i]] > curr_block) {
				blocks.push_back(this->load_block(curr_block, curr_ids));
				if (i == no_used) break;	
				curr_ids.clear(); curr_block = this->snp_map[used_index[i]];
			}
			curr_ids.push_back(used_index[i]);
		}

		Array<double> prod; unsigned int col_offset = 0;
		for (unsigned int b1 = 0; b1 < blocks.size(); b1++) {
			BlockData& curr = *(blocks[b1]); prod.set_size(curr.no_pcs);

			for (unsigned int s1 = 0; s1 < curr.used.size(); s1++) {
				unsigned int snp_id1 = col_offset + s1;
				for (unsigned int i = 0; i < curr.no_pcs; i++) prod[i] = (double) curr.eigenvectors(i,s1) * (double) curr.eigenvalues[i];

				for (unsigned int s2 = s1 + 1; s2 < curr.used.size(); s2++) {
					unsigned int snp_id2 = col_offset + s2; double *read1 = prod.data(), value = 0;
					for (float *read2 = curr.eigenvectors.column(s2), *end = curr.eigenvectors.column(s2+1); read2 < end; read2++) value += *read1++ * (double) *read2;
					ld(snp_id1, snp_id2) = value;
				}
			}

			unsigned int row_offset = 0;
			for (unsigned int b2 = 0; b2 < b1; b2++) {
				BlockData& prev = *(blocks[b2]);
				Buffer<float>& snp_prod = this->load_block_product(prev.id, curr.id);

				if (!snp_prod.is_empty()) {
					ProductMode mode = this->product_info[this->block_index[prev.id][curr.id - prev.id]].mode;
				
					for (unsigned int s1 = 0; s1 < prev.used.size(); s1++) {
						unsigned int snp_id1 = row_offset + s1, u1 = prev.used[s1];

						if (mode == Full) {
							for (unsigned int v = 0; v < curr.no_pcs; v++) {
								double& value = prod[v]; value = 0;
								for (float *read1 = prev.eigenvectors.column(s1), *read2 = snp_prod.column(v), *end = snp_prod.column(v+1); read2 < end; read1++, read2++) value += (double) *read1 * (double) *read2;
							}
						
							for (unsigned int s2 = 0; s2 < curr.used.size(); s2++) {
								double *read1 = prod.data(), value = 0;
								for (float *read2 = curr.eigenvectors.column(s2), *end = curr.eigenvectors.column(s2+1); read2 < end; read2++) value += *read1++ * (double) *read2;
								ld(snp_id1, col_offset + s2) = value;
							}	
						} else if (mode == QtM || mode == MQ) { // MQ is stored transposed, as t(MQ)
							for (unsigned int s2 = 0; s2 < curr.used.size(); s2++) {
								unsigned int c1 = (mode == QtM) ? s1 : u1, c2 = (mode == QtM) ? curr.used[s2] : s2; double value = 0;
								Buffer<float> &R1 = (mode == QtM) ? prev.eigenvectors : snp_prod, &R2 = (mode == QtM) ? snp_prod : curr.eigenvectors;
								for (float *read1 = R1.column(c1), *read2 = R2.column(c2), *end = R2.column(c2+1); read2 < end; read1++, read2++) value += (double) *read1 * (double) *read2;
								ld(snp_id1, col_offset + s2) = value;
							}
						} else for (unsigned int s2 = 0; s2 < curr.used.size(); s2++) {unsigned int u2 = curr.used[s2]; ld(snp_id1, col_offset + s2) = snp_prod(u1,u2);}
						
					}
				}
				row_offset += prev.used.size();
			}
			col_offset += curr.used.size();
		}

		for (unsigned int c = 0; c < no_used; c++) {
			for (unsigned int r = c+1; r < no_used; r++) ld(r,c) = (double) ld(c,r);
		}
		
		for (unsigned int i = 0; i < blocks.size(); i++) delete blocks[i];				
	} else if (no_used == 1) ld(0,0) = 1;

	return(ld);
}



doubles_matrix<> BinaryLD::read_indexed(integers snp_index) {
	Array<unsigned int> step_index = Utils::process_index(snp_index, this->no_snps, true); int no_used = step_index.size();
	writable::doubles_matrix<> ld(no_used,no_used);
	
	if (no_used > 1) {
		for (int i_col = 0, col_id = 0; i_col < no_used; i_col++) {col_id += step_index[i_col];
			auto column = ld[i_col]; column[i_col] = 1;
			for (int i_row = 0; i_row < i_col; i_row++) column[i_row] = (double) ld(i_col,i_row);
			
			Buffer<float>& corrs = this->read_line(col_id); unsigned int* step = step_index.data() + i_col + 1;
			float *read = corrs.data() - 1, *read_end = corrs.end(); // shifting read back by one since buffer starts at first correlation rather than current SNP itself
			for (int i_row = i_col + 1; i_row < no_used; i_row++) {
				read += *step++;
				if (read >= read_end) while (i_row < no_used) column[i_row++] = 0;
				else column[i_row] = *read;
			}
		}
	} else if (no_used == 1) ld(0,0) = 1;
	
	return(ld);
}


doubles_matrix<> BinaryLD::read_ld(integers snp_index) {
	if (this->storage_format == EigenDecomposed) return(this->read_decomposed(snp_index));
	else return(this->read_indexed(snp_index));
}



// NB: will skip trailing zeroes
Buffer<float>& BinaryLD::read_line(unsigned int snp_id) {
	if (snp_id >= this->no_snps) error("SNP ID out of bounds");
	
	if (this->no_corrs[snp_id] > 0) {
		unsigned long long offset = this->value_offset[snp_id], read_count = this->value_offset[snp_id+1] - offset; unsigned int no_values = this->no_corrs[snp_id];
		if (this->storage_format != Truncated && read_count != no_values) error("inconsistent read count");
		
		if (this->storage_type == Float) return this->read_float(offset, read_count, no_values);
		else if (this->storage_type == Integer8) return this->read_integer<uint8_t>(offset, read_count, no_values);		
		else if (this->storage_type == Integer16) return this->read_integer<uint16_t>(offset, read_count, no_values);
		else if (this->storage_type == Integer9) return this->read_integer<uint8_t,1>(offset, read_count, no_values);
		else if (this->storage_type == Integer10) return this->read_integer<uint8_t,2>(offset, read_count, no_values);		
		else if (this->storage_type == Integer12) return this->read_integer<uint8_t,4>(offset, read_count, no_values);		
		else error("storage type is unknown");
	} else this->read_buffer.set_size(0);
	
	return this->read_buffer;
}
	

Buffer<float>& BinaryLD::read_float(unsigned long long offset, unsigned long long read_count, unsigned int no_values) {
	static_assert(sizeof(float) == 4, "code assumes float values have a size of four bytes");	
	this->read_buffer.set_size(read_count); this->bcor_file.seekg(offset * this->read_size + this->read_offset);
	this->bcor_file.read(reinterpret_cast<char*>(this->read_buffer.data()), read_count * this->read_size);
	
	if (this->storage_format == Truncated) {
		static Buffer<float> process_buffer; process_buffer.set_size(no_values);	
		float* write = process_buffer.data();	int zero_count = 0;
		for (float *read = this->read_buffer.data(), *read_end = this->read_buffer.end(); read < read_end; read++) {
			if (*read != 0) {
				while (zero_count-- > 0) *write++ = 0;
				*write++ = *read;
			} else zero_count = *reinterpret_cast<uint32_t*>(++read);
		}
		
		process_buffer.set_size(write - process_buffer.data(), 1, true);
		if ((process_buffer.size() + (zero_count > 0 ? zero_count : 0)) != no_values) stop("inconsistent number of loaded values");
		return process_buffer;
	} else return this->read_buffer;
}


BinaryLD::BlockData::BlockData(unsigned int id, BlockInfo info, vector<unsigned int>& snps) : id(id), no_snps(info.no_snps), no_pcs(info.no_pcs) {
	this->used.set_size(snps.size());

	for (unsigned int i = 0; i < snps.size(); i++) {
		if (snps[i] < info.id_offset || snps[i] >= info.id_offset + this->no_snps) error("SNP index in block data out of bounds");
		this->used[i] = snps[i] - info.id_offset;
	}
}



