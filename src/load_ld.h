#ifndef LOAD_LD_H
#define LOAD_LD_H


#include "global.h"
#include "utils.h"

#include <fstream>
#include <vector>





[[cpp11::register]]
doubles_matrix<> load_ld(std::string prefix, int chromosome, int total_snps, integers snp_index); 



class BinaryLD {
	enum HeaderField {CheckValue = 0, HeaderSize = 1, NoSNPs = 2, ValueType = 3, ValueSize = 4, FormatID = 5, FormatParamStart = 6};
	
	ifstream bcor_file; string prefix; unsigned int chromosome, no_snps;
	enum StorageFormat{UnknownFormat = 0, Indexed = 11, Truncated = 12, EigenDecomposed = 14} storage_format = UnknownFormat; 
	enum StorageType {UnknownType = 0, Float = 3, Integer8 = 4, Integer16 = 5, Integer9 = 6, Integer10 = 7, Integer12 = 8} storage_type = UnknownType; 
	Buffer<float> read_buffer; unsigned int read_offset, read_size; unsigned long long index_offset[2] = {0,0}; // file offsets for indices stored at end of file
	Array<unsigned long long> value_offset; // file offsets for aggregate data entries (lines of correlations or data blocks); in bytes if EigenDecomposed, in units of size read_size otherwise
	Array<uint32_t> no_corrs; // variables specific to Indexed/Truncated format
	
	// variables specific to EigenDecomposed format
	enum ProductMode {Raw = 1, QtM = 2, MQ = 3, Full = 4}; 
	struct ProductInfo {unsigned int no_values, no_bytes; ProductMode mode; bool as_integer; ProductInfo(unsigned int no_values, unsigned int no_bytes, ProductMode mode, bool as_integer) : no_values(no_values), no_bytes(no_bytes), mode(mode), as_integer(as_integer) {}};
	struct BlockInfo {unsigned int no_snps, no_pcs, id_offset; BlockInfo(unsigned int no_snps, unsigned int no_pcs, unsigned int id_offset) : no_snps(no_snps), no_pcs(no_pcs), id_offset(id_offset) {}}; struct BlockData;
	unsigned int no_blocks, no_entries; vector<BlockInfo> block_sizes; vector<ProductInfo> product_info;
	vector<vector<unsigned int>> block_index; // map of block IDs onto value_offset (main block, then cross blocks)
	Array<unsigned int> snp_map; // map of SNPs to block IDs
	
	
	template<typename T> T* get_values(unsigned long long offset, T* target, unsigned int total=1) {streampos current = this->bcor_file.tellg(); this->bcor_file.seekg(offset).read(reinterpret_cast<char*>(target), sizeof(T)*total).seekg(current); return target;}
	uint32_t header_value(HeaderField field) {uint32_t out; this->get_values(4*field, &out); return out;}
	uint32_t format_param(unsigned int param_id) {uint32_t out; this->get_values(4*(FormatParamStart+param_id), &out); return out;}
	uint64_t format_param_long(unsigned int param_id) {uint64_t out; this->get_values(4*(FormatParamStart+param_id), &out); return out;}
	float format_param_float(unsigned int param_id) {float out; this->get_values(4*(FormatParamStart+param_id), &out); return out;}

	void prep_bcor();
	void load_index();
	
	Buffer<float>& read_line(unsigned int snp_id);
	Buffer<float>& read_float(unsigned long long offset, unsigned long long read_count, unsigned int no_values);
	template<typename T> Buffer<float>& read_integer(unsigned long long offset, unsigned long long read_count, unsigned int no_values);
	template<typename T, unsigned int ADD> Buffer<float>& read_integer(unsigned long long offset, unsigned long long read_count, unsigned int no_values);
	BlockData* load_block(unsigned int id, vector<unsigned int>& snps);
	Buffer<float>& load_block_product(unsigned int id1, unsigned int id2);

	doubles_matrix<> read_indexed(integers snp_index);
	doubles_matrix<> read_decomposed(integers snp_index);
public:
	BinaryLD(const string& file_prefix, unsigned int chromosome, unsigned int total_snps);
	
	doubles_matrix<> read_ld(integers snp_index);
};

struct BinaryLD::BlockData {
	unsigned int id, no_snps, no_pcs; Array<unsigned int> used;
	Buffer<float> eigenvectors; Array<float> eigenvalues;
	
	BlockData(unsigned int id, BlockInfo info, vector<unsigned int>& snps);
};


#include "load_ld.tpl.h"



#endif /* LOAD_LD_H */
