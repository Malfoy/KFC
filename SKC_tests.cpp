#include <iostream>

#include "Kmers.hpp"
#include "SuperKmerCount.hpp"

using namespace std;


uint64_t str2num(const string& str) {
	uint64_t res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}

int main(int argc, char** argv) {
	// Init skc
	SKC skc(str2num("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC"), 1);
	cout << skc << endl;
	// Add kmer on right
	skc.add_kmer({2, str2num("AAAAAAAAAAAAAAAAAAAAAAAAAAAAACG"), str2num("CGTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")});
	cout << skc << endl;
	// Add kmer already present
	skc.add_kmer({1, str2num("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC"), str2num("GTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")});
	cout << skc << endl;
	// Compact left
	skc.add_kmer({0, str2num("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), str2num("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTA")});
	cout << skc << endl;
	// Same on reverse complement strand
	skc.add_kmer({24, str2num("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTA"), str2num("TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")});
	cout << skc << endl;
	// Insert after on reverse complement
	skc.add_kmer({3, str2num("CCGTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), str2num("AAAAAAAAAAAAAAAAAAAAAAAAAAAACGG")});
	cout << skc << endl;
}
