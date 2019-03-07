#ifndef __INDEX8MIN
#define __INDEX8MIN

#include <stdio.h>
#include <fstream>
#include <vector>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <algorithm>
#include "BooPHF.h"

typedef boomphf::SingleHashFunctor<uint64_t> hasher;
typedef boomphf::mphf<uint64_t, hasher> MPHF;

// FIXME: maybe it would be better to use a structure of arrays instead of an array of packed struct:
struct __attribute__((packed)) value {
	uint64_t kmer;
	uint8_t count;
};

class index_full {
  public:
	MPHF Hash;
	std::vector<value> Values;
	std::vector<uint64_t> weak_kmer_buffer;
	std::ofstream dump_weak;
	static const uint32_t kmer_size = 31;
	static const size_t weak_kmer_buffer_waterline = 1 << 14;

	index_full(const std::vector<uint64_t>& V) {
		weak_kmer_buffer.reserve(weak_kmer_buffer_waterline);
		dump_weak.open("weak_kmers", std::ofstream::out | std::ofstream::binary);
		Hash = boomphf::mphf<uint64_t, hasher>(V.size(), V, 4, 5, false, true, 0.003f);
		Values.resize(V.size());
		for (auto& kmer : V)
			Values[Hash.lookup(kmer)].kmer = kmer;
	}

	void insert(uint64_t);
	void dump_counting(std::ostream& stream = std::cout);
	void insert_seq(const std::string& read);
	void clear(bool = false);
	void print_kmer(uint64_t num, std::ostream& steam);
};

uint64_t rcb(uint64_t min, unsigned n);

class index_min {
  public:
	unsigned kmer_size;
	unsigned minimizer_size;

	std::vector<index_full> Index;

	index_min(std::vector<uint64_t>& V);
	void insert(uint64_t kmer);
	void dump_counting(std::ostream& stream = std::cout);
};

#endif
