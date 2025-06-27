#ifndef GLOBAL_H
#define GLOBAL_H

#include <cpp11.hpp>

#include <sys/stat.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace cpp11;


class ErrorMsg : public exception {
	string message;
public:
	ErrorMsg(const string& msg) : message(msg) {}
	
	virtual const char* what() const throw() {return this->message.c_str();}
};
	

static void error(const string& msg) {throw ErrorMsg(msg);}



namespace Utils {
	inline unsigned int check_chromosome(int chr) {if (chr <= 0 || chr > 24) error("invalid chromosome specified"); return chr;}

	inline bool is_dir(const string& basename) {struct stat status; return stat(basename.c_str(), &status) == 0 && S_ISDIR(status.st_mode);}
	inline bool is_file(const string& filename) {struct stat status; int code = stat(filename.c_str(), &status); return code == 0 && S_ISREG(status.st_mode);}
	inline const string& check_file(const string& file) {if (is_dir(file)) error("file '" + file + "' is a directory"); if (!is_file(file)) error("file '" + file + "' not found"); return file;}

	template<typename T> string to_string(T num) {ostringstream ostr; return static_cast<ostringstream*>( &(ostr << num) )->str();}
	
	
	template<typename T> bool check_zero_bits() {T zero = 0, bits; memset(&bits, 0, sizeof(T)); return bits == zero;}
	template<typename T> void set_zero(T* target, unsigned int size) {
		static bool zero_loop = check_zero_bits<T>();		
		if (zero_loop) for (T *end = target + size; target < end; target++) *target = 0;
		else memset(target, 0, size*sizeof(T));
	}
};




#endif /* GLOBAL_H */
