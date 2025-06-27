#ifndef LOAD_PLINK_H
#define LOAD_PLINK_H


#include "global.h"
#include "utils.h"

#include <fstream>




[[cpp11::register]]
list load_plink(std::string prefix, int total_snps, int total_indiv, integers snp_index, list parameters);

[[cpp11::register]]
list compute_ld(std::string prefix, int total_snps, int total_indiv, integers snp_index, list parameters); 


[[cpp11::register]]
doubles_matrix<> compute_info(std::string prefix, int total_snps, int total_indiv, integers snp_index); 



enum OutputMode {Genotypes = 10, Correlations = 20};
list load(std::string prefix, int total_snps, int total_indiv, integers snp_index, list parameters, OutputMode mode); 


class BinaryPlink {
	enum ValueType {Hom1 = 0, HetGen = 1, Hom2 = 2, Missing = 3};

	ifstream bed_file; string prefix; unsigned int no_snps, no_individuals;
	uint32_t geno_index[256]; unsigned int line_size;
	Buffer<unsigned char> read_buffer; Buffer<unsigned int> count_buffer;
	Array<unsigned int> missing_buffer; Array<unsigned int*> missing_index; 
	Array<unsigned int> snp_include; unsigned int mac_thresh = 1; double maf_thresh = 0, missing_thresh = 0.05;

	void prep_bed();

	bool check_snp(unsigned int* counts);
	unsigned int check_input();

	unsigned int product_sum(unsigned char* data1, unsigned char* data2);
	void update_sums(unsigned char* data, unsigned int* miss_ptr, unsigned int miss_count, unsigned int& sum, unsigned int& sq_sum);
		
public:
	BinaryPlink(const string& file_prefix, unsigned int total_snps, unsigned int total_indiv);
	
	void clear() {this->read_buffer.clear(); this->count_buffer.clear(); this->missing_buffer.clear(); this->missing_index.clear();}
	BinaryPlink& read_data(integers snp_index);
	
	BinaryPlink& set_thresholds(unsigned int mac, double maf, double missing) {this->mac_thresh = mac; this->maf_thresh = maf; this->missing_thresh = missing; return *this;}
	integers get_included() {writable::integers index(this->snp_include.size()); for (int i = 0; i < index.size(); i++) index[i] = this->snp_include[i]; return index;}

	doubles_matrix<> compute_info(integers snp_index);
		
	integers_matrix<> load_genotypes();
	doubles_matrix<> compute_ld();
};






#endif /* LOAD_PLINK_H */
