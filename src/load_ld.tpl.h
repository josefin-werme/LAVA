

template<typename T, unsigned int ADD> 
Buffer<float>& BinaryLD::read_integer(unsigned long long offset, unsigned long long read_count, unsigned int no_values) {
	static_assert(numeric_limits<T>::is_integer && !numeric_limits<T>::is_signed, "type must be an unsigned integer type");
	static_assert(ADD < sizeof(T) * 8 && (ADD == 0 || (ADD & (ADD - 1)) == 0), "number of added bits must be zero or a power of two, and no more than half of the bit size of the used type"); 

	static Array<T> conversion_buffer; conversion_buffer.set_size(read_count); T *input = conversion_buffer.data(), *input_end = conversion_buffer.end(); 
	this->bcor_file.seekg(offset * this->read_size + this->read_offset).read(reinterpret_cast<char*>(input), read_count * this->read_size);

	const unsigned int range_scale = pow(2, ADD), add_count = (ADD > 0) ? 8 * sizeof(T) / ADD : 1; 
	const T range_increment = numeric_limits<T>::max(), zero_value = range_increment, bitmask = range_scale - 1;
	float value_scale = (range_scale * range_increment - 1) / 2.0;
	
	this->read_buffer.set_size(no_values); float* write = read_buffer.data();
	if (this->storage_format == Truncated) {
		int zero_count = 0, no_added = add_count; T add_value; 
		for (; input < input_end; input++) {
			if (*input != zero_value) {
				while (zero_count-- > 0) *write++ = 0;
				unsigned int value = *input;
				if (ADD > 0) {
					if (no_added >= (int) add_count) {no_added = 0; add_value = *(++input);}
					value += ((add_value >> ADD * no_added++) & bitmask) * range_increment;
				}
				*write++ = value / value_scale - 1;
			} else zero_count = (zero_count > 0 ? zero_count : 0) + *(++input);
		}
		
		read_buffer.set_size(write - read_buffer.data(), 1, true);
		if ((read_buffer.size() + (zero_count > 0 ? zero_count : 0)) != no_values) stop("inconsistent number of loaded values");
	} else while (input < input_end) *write++ = *input++ / value_scale - 1;	
	
	return this->read_buffer;
}



