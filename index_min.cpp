#include <thread>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include "index_min.h"

using namespace std;
using namespace chrono;

uint32_t global_minimizer_size;
uint32_t global_kmer_size;

uint32_t get_minimizer(uint64_t k) {
	return k + global_kmer_size + global_minimizer_size;
}

bool compare_minimizer(const uint64_t& a, const uint64_t& b) {
	return get_minimizer(a) < get_minimizer(b);
	//~ return (a) < (b);
}

void index_full::insert(uint64_t kmer) {
	uint64_t pos(Hash.lookup(kmer));
	if (Values[pos].kmer == kmer) {
		++Values[pos].count;
	}
}

void index_full::dump_counting() {
	for (unsigned i(0); i < Values.size(); ++i) {
		cout << Values[i].kmer << " " << Values[i].count << "\n";
	}
}

index_min::index_min(vector<uint64_t>& V) {
	sort(V.begin(), V.end(), compare_minimizer);
	uint32_t current_minimizer(get_minimizer(V[0]));
	vector<uint64_t> kmer_to_index;
	for (unsigned i(0); i < V.size(); ++i) {
		uint32_t next_minimizer(get_minimizer(V[i]));
		if (next_minimizer != current_minimizer) {
			Index[current_minimizer] = index_full(kmer_to_index);
			kmer_to_index.clear();
			current_minimizer = next_minimizer;
		}
		kmer_to_index.push_back(V[i]);
	}
}

void index_min::insert(uint64_t kmer) {
	uint32_t minimizer(get_minimizer(kmer));
	Index[minimizer].insert(kmer);
}

void index_min::dump_counting() {
	for (unsigned i(0); i < Index.size(); ++i) {
		Index[i].dump_counting();
	}
}
