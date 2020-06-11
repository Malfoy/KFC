#include <emmintrin.h>
#include <tmmintrin.h>
#include <cstdint>
#include <string>
#include "robin_hood.h"




using namespace std;



#ifndef KMERS_H
#define KMERS_H
extern robin_hood::unordered_flat_map<string, uint64_t> real_count;
extern uint64_t counting_errors;
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
// Return the canonical minimizer for a uint64 sequence.
int64_t get_minimizer(uint64_t seq, int8_t& position);
uint64_t reversebits(uint64_t b);
string getCanonical(const string& str);




// ----- Kmer class -----
class kmer_full {
public:
	int8_t minimizer_idx;
	uint64_t kmer_s;
	uint64_t prefix;
	uint64_t suffix;
	kmer_full(int8_t minimizer_idx, uint64_t value);
	uint64_t get_compacted();
	/** Return the minimizer regarding the minimizer_idx property
		* Warning: The minimizer should be canon
		*/
	uint8_t get_minimizer_idx() const;
	// uint64_t get_minimizer() const;
	bool contains_multi_minimizer() const;
};

#endif
