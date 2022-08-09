#pragma once

// libraries

#include "container_traits.h"

#include <fstream>

// utils

using serialization_traits_container_size_t = uint32_t;


//===================================//
//===standart serialization_traits===//
//===================================//


// standart serialization_traits
template<typename T, typename enable_if_type = void>
struct serialization_traits {

	static void serialize(std::ostream& os, const T& value) {
		value.serialize(os);
	}

	static T deserialize(std::istream& is) {
		T result;
		result.deserialize(is);
		return result;
	}
};


//==================================//
//========byte serialization========//
//==================================//


template<typename T>
void byte_serialization(std::ostream& os, const T& value) {

	os.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T>
T byte_deserialization(std::istream& is) {
	T result;

	is.read(reinterpret_cast<char*>(&result), sizeof(T));

	return result;
}


// macro for bytes serialization in structures
#define serialization_traits_byte(type)\
void serialize(std::ostream& os) const {\
	byte_serialization(os, *this);\
}\
void deserialize(std::istream& is) {\
	*this = byte_deserialization<type>(is);\
}


//===================================//
//====others serialization_traits====//
//===================================//


// serialization_traits for arithmetic: char, int, long long, double ...
template<typename T>
struct serialization_traits<T, typename std::enable_if<std::is_arithmetic<T>::value>::type> {

	static void serialize(std::ostream& os, const T& x) {
		byte_serialization(os, x);
	}

	static T deserialize(std::istream& is) {
		return byte_deserialization<T>(is);
	}
};


// serialization_traits for container
template<typename container_t>
struct serialization_traits<container_t, typename std::enable_if<is_container<container_t>::value>::type> {

	static void serialize(std::ostream& os, const container_t& v) {
		byte_serialization<serialization_traits_container_size_t>(os, v.size());

		for (auto& it : v) {
			serialization_traits<container_t::value_type>::serialize(os, it);
		}
	}

	static container_t deserialize(std::istream& is) {

		serialization_traits_container_size_t len;
		len = byte_deserialization<decltype(len)>(is);

		std::vector<container_t::value_type> tmp(len);

		for (size_t i = 0; i < len; i++) {
			tmp[i] = serialization_traits<container_t::value_type>::deserialize(is);
		}

		return container_t(tmp.begin(), tmp.end());
	}
};


//==================================//
//===============file===============//
//==================================//


template<typename T>
void save_to_file(const std::string& file_name, const T& data) {
	std::ofstream os(file_name, std::ios::binary);
	serialization_traits<T>::serialize(os, data);
}

template<typename T>
T read_from_file(const std::string& file_name) {
	std::ifstream is(file_name, std::ios::binary);
	return serialization_traits<T>::deserialize(is);
}

