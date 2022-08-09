#include "serializer.h"


#include <iostream>
#include <vector>


int main() {

	{
		std::string s = "hello world";

		save_to_file("in.txt", s);

		s.clear();

		auto s_new = read_from_file<std::string>("in.txt");

		std::cout << s_new << "\n";
	}
	std::cout << "\n";

	{
		struct dot {
			double x, y;

			serialization_traits_byte(dot)
		};

		dot p;
		p.x = 8123.123;
		p.y = 12.091;

		save_to_file("in.txt", p);

		dot g = read_from_file<dot>("in.txt");

		std::cout << g.x << " " << g.y << "\n";
	}
	std::cout << "\n";

	{
		std::vector<double> dbls = { 12.213, 912e1, 812.1, INFINITY };

		save_to_file("in.txt", dbls);

		dbls.clear();

		auto dbls_new = read_from_file<std::vector<double>>("in.txt");

		for (auto& it : dbls_new) {
			std::cout << it << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	{
		struct my_vector {
			std::vector<std::string> Strs;

			int g;

			void serialize(std::ostream& os) const {

				serialization_traits<std::vector<std::string>>::serialize(os, Strs);

				serialization_traits<int>::serialize(os, g);
			}
			void deserialize(std::istream& is) {
				
				Strs = serialization_traits<std::vector<std::string>>::deserialize(is);

				g = serialization_traits<int>::deserialize(is);
			}

		};

		std::vector<my_vector> V = {
			{{"hi", "me", "I"}, 10},

			{{"hello", "world", "!"}, 8712},

			{{"aqk", "all", "home"}, -129},
		};

		save_to_file("in.txt", V);

		V.clear();

		auto nV = read_from_file<std::vector<my_vector>>("in.txt");

		for (auto& it : nV) {
			for (auto& jt : it.Strs) {
				std::cout << jt << " ";
			}
			std::cout << "\n" << it.g << "\n";
		}
	}
	std::cout << "\n";

	{
		std::vector<std::vector<long long>> A = {
			{1, -8, 17},
			{9, 12, 5},
			{-4, 3, 2},
		};

		save_to_file("in.txt", A);

		A.clear();

		auto B = read_from_file<std::vector<std::vector<long long>>>("in.txt");

		for (auto& row : B) {
			for (auto& col : row) {
				std::cout << col << " ";
			}
			std::cout << "\n";
		}
	}
	std::cout << "\n";

    return 0;
}
