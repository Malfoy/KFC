#include <emmintrin.h>
#include <tmmintrin.h>
#include <cstdint>
#include <string>

#ifndef KMERS_H
#define KMERS_H

extern uint64_t k;
extern uint64_t k_mask;
extern uint64_t minimizer_size;
extern uint64_t min_mask;

// ----- Usefull binary kmer functions -----

void print_kmer(__uint128_t num,uint64_t n);
std::string kmer2str(__uint128_t num, uint k);
uint64_t str2num(const std::string& str);

// RC functions
uint64_t rcb(const uint64_t& in);
__uint128_t rcb(const __uint128_t& in);
uint64_t rcbc(uint64_t in, uint64_t n);

// Hash functions (TODO: move them)
__m128i mm_bitshift_right(__m128i x, unsigned count);
__m128i mm_bitshift_left(__m128i x, unsigned count);
uint64_t hash64shift(uint64_t key);

/** Return the canonical minimizer for a uint64 sequence.
	*/
uint64_t get_minimizer(uint64_t seq, uint64_t& position);

// ----- Kmer class -----

class kmer_full {
public:
	uint8_t minimizer_idx;
	uint64_t kmer_s;
	uint64_t kmer_rc;

	kmer_full(uint8_t minimizer_idx, uint64_t value, uint64_t reverse_comp_value);
	/** Return the minimizer regarding the minimizer_idx property
		* Warning: The minimizer can be non canonical !
		*/
	uint64_t get_minimizer();
};

#endif
