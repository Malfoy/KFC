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
	return static_cast<uint32_t>(k + global_kmer_size + global_minimizer_size); // WTF?
}

bool compare_minimizer(const uint64_t& a, const uint64_t& b) {
	return get_minimizer(a) < get_minimizer(b);
	//~ return (a) < (b);
}

void index_full::insert(uint64_t kmer) {
	uint64_t pos(Hash.lookup(kmer));
	if (pos == ULLONG_MAX) {
		weak_kmer_buffer.push_back(kmer);
	} else {
		if (Values[pos].kmer == kmer) {
			++Values[pos].count;
			//~ cout<<Values[pos].count<<endl;
		} else {
			weak_kmer_buffer.push_back(kmer);
		}
	}
}

uint64_t str2num(const string& str) {
	uint64_t res(0);
	for (unsigned i(0); i < str.size(); i++) {
		res <<= 2;
		switch (str[i]) {
			case 'A': res += 0; break;
			case 'C': res += 1; break;
			case 'T': res += 2; break;
			default: res += 3; break;
		}
	}
	return res;
}

#define hash_letter(letter) ((letter >> 1) & 3)

uint64_t rcb(uint64_t in, unsigned n) {
	//~ assume(n <= 32, "n=%u > 32", n);
	// Complement, swap byte order
	uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
	res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

	// Realign to the right
	res >>= 64 - 2 * n;
	return res;
}

void index_full::insert_seq(const string& seq) {
	uint64_t hash = 0;
	uint64_t canon_hash = 0;
	for (unsigned i = 0; i < 31; i++)
		hash = hash << 2 | hash_letter(seq[i]);

	for (unsigned idx = 31; idx < seq.size() /**/; idx++) {
		hash = hash << 2 | hash_letter(seq[idx]);
		canon_hash = min(hash, rcb(hash, 31));
		//~ print_kmer(hash);
		//~ print_kmer(rcb(hash,31));
		//~ print_kmer(rcb(rcb(hash,31),31))	;
		//~ cin.get();
		insert(canon_hash);
	}
}

//~ void index_full::insert_seq(const string&  read){
//~ for (unsigned i(0);i+kmer_size<read.size();++i){
//~ uint64_t seq(str2num(read.substr(i,kmer_size)));
//~ insert(seq);
//~ }
//~ }

void index_full::print_kmer(uint64_t num, std::ostream& stream) {
	uint64_t anc(1);
	anc <<= (2 * (kmer_size));
	for (unsigned i(0); i < kmer_size; ++i) {
		auto nuc = num / anc;
		if (nuc >= 4) {
			stream << nuc << endl;
			stream << "WTF" << endl;
		}
		stream << "ACTG"[num];
		anc >>= 2;
	}
}

void index_full::dump_counting(std::ostream& stream) {
	for (unsigned i(0); i < Values.size(); ++i) {
		if (Values[i].count != 0) {
			print_kmer(Values[i].kmer, stream);
			stream << " " << to_string(Values[i].count) << "\n";
		}
	}
}

void index_full::clear(bool force) {
	//~ return;
	if (weak_kmer_buffer.size() > weak_kmer_buffer_waterline or force) {
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

void index_min::dump_counting(std::ostream& stream) {
	for (unsigned i(0); i < Index.size(); ++i) {
		Index[i].dump_counting(stream);
	}
}
