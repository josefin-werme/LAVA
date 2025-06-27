#ifndef UTILS_H
#define UTILS_H


#include "global.h"


// retains capacity, does not reallocate when resizing (except by call to clear())
template<typename T>
class Array {
	T* content = 0;	unsigned long capacity = 0, elements = 0;

	Array<T>& operator=(const Array<T>& other) = delete;
	Array(const Array<T>&) = delete; 
	

public:
	Array(unsigned long init_size=0, bool zero=false) {this->set_size(init_size); if (zero) this->zero_fill();}
	Array(Array<T>&& other) : content(other.content), capacity(other.capacity), elements(other.elements) {other.content = 0;}
	~Array() {this->clear();}

	void clear() {delete[] this->content; this->content = 0; this->elements = this->capacity = 0;}
	
	//if keep_data=true, keeps (remaining) data after resizing, even if reallocating
	Array<T>& set_size(unsigned long e, bool keep_data=false) {
		if (e > 0) {
			if (e > this->capacity) {
				T* new_content = new T[e]; 
				if (keep_data) memcpy(new_content, this->content, this->elements * sizeof(T));
				delete[] this->content; this->content = new_content;
				this->elements = this->capacity = e;
			} else this->elements = e;
		} else this->clear();
		return *this;
	}

	Array<T>& zero_fill(bool used_only=true) {
		if (!this->is_empty()) Utils::set_zero(this->content, used_only ? this->elements : this->capacity);
		return *this;
	}

	unsigned long size() const {return !this->is_empty() ? this->elements : 0;}
	unsigned long available() const {return this->capacity;}
	bool is_empty() const {return this->content == 0 || this->elements == 0;}

	T& operator() (unsigned long e) {return this->content[e];}
	const T& operator() (unsigned long e) const {return this->content[e];}	
	T& operator[] (unsigned long e) {return this->content[e];}
	const T& operator[] (unsigned long e) const {return this->content[e];}	

	T* data() {return this->content;}
	const T* data() const {return this->content;}
	T* end() {return this->content + this->size();}
	const T* end() const {return this->content + this->size();}
};



// retains capacity, does not reallocate when resizing (except by call to clear())
template<typename T>
class Buffer {
	T* content = 0; unsigned long capacity = 0, rows = 0, cols = 0;
	
	Buffer<T>& operator=(const Buffer<T>& other) = delete;
	Buffer(const Buffer<T>&) = delete; 
	
	void resize(unsigned long r, unsigned long c, bool keep_data) {
		T* new_content = new T[r*c];
		if (keep_data && !this->is_empty()) memcpy(new_content, this->content, min(this->cols*this->rows, c*r) * sizeof(T));
		delete[] this->content; this->content = new_content;
		this->rows = r; this->cols = c;	this->capacity = c*r;
	}
	
public:
	Buffer(unsigned long init_rows=0, unsigned long init_cols=1, bool zero=false) {this->set_size(init_rows, init_cols); if (zero) this->zero_fill();} 
	Buffer(Buffer<T>&& other) : content(other.content), capacity(other.capacity), rows(other.elements), cols(other.cols) {other.content = 0;}
	~Buffer() {this->clear();}
	
	void clear() {delete[] this->content; this->content = 0; this->rows = this->cols = this->capacity = 0;}
	
	Buffer<T>& set_size(unsigned long r, unsigned long c=1, bool keep_data=false) {
		if (r*c > 0) {
			if (r*c > this->capacity) this->resize(r, c, keep_data);
			else {this->rows = r; this->cols = c;}
		} else this->clear();
		return *this;
	}

	Buffer<T>& zero_fill(bool used_only=true) {
		if (!this->is_empty()) Utils::set_zero(this->content, used_only ? this->rows * this->cols : this->capacity);
		return *this;
	}
	
	unsigned long nrow() const {return !this->is_empty() ? this->rows : 0;}
	unsigned long ncol() const {return !this->is_empty() ? this->cols : 0;}
	unsigned long size() const {return !this->is_empty() ? this->rows*this->cols : 0;}
	unsigned long available(bool columns=false) const {return columns ? (this->rows > 0 ? this->capacity / this->rows : 0) : this->capacity;}
	bool is_empty() const {return this->content == 0 || this->rows*this->cols == 0;}
	
	T& operator()(const unsigned long r, const unsigned long c=0) {return this->content[r+c*this->rows];}	
	const T& operator()(const unsigned long r, const unsigned long c=0) const {return this->content[r+c*this->rows];}	
	
	T* column(const unsigned long c) {return this->content+c*this->rows;}
	const T* column(const unsigned long c) const {return this->content+c*this->rows;}
	T* data() {return this->content;}
	const T* data() const {return this->content;}
	T* end() {return this->content + this->size();}
	const T* end() const {return this->content + this->size();}
};



namespace Utils {
	//validate index of (numeric) SNP IDs, and convert into an array of step sizes if as_steps=true
	inline Array<unsigned int> process_index(integers snp_index, int max_snps, bool as_steps=true) {	
		Array<unsigned int> out_index(snp_index.size()); int no_used = snp_index.size(), previous = -1; 
		for (int i = 0; i < no_used; i++) {
			if (snp_index[i] < 0 || snp_index[i] >= max_snps) error("SNP ID out of bounds");
			if (snp_index[i] <= previous) error(snp_index[i] == previous ? "duplicate SNP ID" : "SNP index is not ordered");
			previous = snp_index[i];
			
			if (as_steps) out_index[i] = (i > 0) ? snp_index[i] - snp_index[i-1] : snp_index[i];
			else out_index[i] = snp_index[i];
		}
		return out_index;
	}

	inline cpp11::doubles list_values(list input, const string& name) {if (!Rf_isNull(input[name])) return doubles(input[name]); else return(doubles());}
	inline double list_value(list input, const string& name, double default_value=0) {
		doubles values = list_values(input, name);
		return (values.size() > 0 ? values[0] : default_value);
	}
		

};



#endif /* UTILS_H */
