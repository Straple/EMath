// tester
// Developed by Mob


#include <bits/stdc++.h>
using namespace std;

#include<random>
mt19937 rnd(time(NULL));

// random value between [l, r]
int randint(int l, int r) {
	return rnd() % (r - l + 1) + l;
}

// строит тест
void build_test() {
	ofstream cout("test.txt");

	// вывести входные данные
}

vector<string> input(ifstream& in) {
	vector<string> A;
	string s;
	while (in >> s) {
		A.push_back(s);
	}
	return A;
}

bool check_for_WA() {
	ifstream in1("result1.txt");
	ifstream in2("result2.txt");

	auto A1 = input(in1);
	auto A2 = input(in2);

	return A1 != A2;
}

#include "solve1.cpp"
#include "solve2.cpp"

int main() {

	for (int cnt_test = 0; ; cnt_test++) {
		
		build_test();

		solve1();
		solve2();

		// проверить, что они равны
		if (check_for_WA()) {
			cout << "WA\n";
			break;
		}

		cout << cnt_test << "\n";
	}
	
	return 0;
}
