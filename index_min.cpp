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
	} else {
		weak_kmer_buffer.push_back(kmer);
	}
}

uint64_t str2num(const string& str) {
	uint64_t res(0);
	for (unsigned i(0); i < str.size(); i++) {
		res <<= 2;
		switch (str[i]) {
			case 'A': res += 0; break;
			case 'C': res += 1; break;
			case 'G': res += 2; break;
			default: res += 3; break;
		}
	}
	return res;
}

void index_full::insert_seq(const string& read) {
	for (unsigned i(0); i + kmer_size < read.size(); ++i) {
		uint64_t seq(str2num(read.substr(i, kmer_size)));
		insert(seq);
	}
}

void index_full::dump_counting() {
	for (unsigned i(0); i < Values.size(); ++i) {
		cout << Values[i].kmer << " " << Values[i].count << "\n";
	}
}

void index_full::clear(bool force = false) {
	if (weak_kmer_buffer.size() > 10 * 1000 or force) {
		for (unsigned i(0); i < weak_kmer_buffer.size(); ++i) {
			dump_weak << weak_kmer_buffer[i];
		}
		weak_kmer_buffer.clear();
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
