
#include <cstdint>

#ifndef KMERS_H
#define KMERS_H

extern uint64_t k;
extern uint64_t k_mask;
extern uint64_t minimizer_size;
extern uint64_t min_mask;

// ----- Usefull binary kmer functions -----

void print_kmer(__uint128_t num,uint64_t n);

// ----- Kmer class -----

class kmer_full {
public:
	uint8_t minimizer_idx;
	uint64_t kmer_s;
	uint64_t kmer_rc;

	kmer_full(uint8_t minimizer_idx, uint64_t value, uint64_t reverse_comp_value);
};

#endif
