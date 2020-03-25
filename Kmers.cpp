#include <iostream>

#include "Kmers.hpp"
#include "pow2.hpp"


using namespace std;

uint64_t k = 31;
uint64_t k_mask = (((uint64_t)1) << (2*k)) - 1;
uint64_t minimizer_size = 15;
uint64_t min_mask = (((uint64_t)1) << (2*minimizer_size)) - 1;


// ----- Useful binary kmer functions -----

void print_kmer(__uint128_t num,uint64_t n){
	__uint128_t anc((__uint128_t)1<<(2*(n-1)));
	for(uint64_t i(0);i<n and anc!=0;++i){
		uint64_t nuc=num/anc;
		num=num%anc;
		if(nuc==2){
			cout<<"T";
		}
		if(nuc==3){
			cout<<"G";
		}
		if(nuc==1){
			cout<<"C";
		}
		if(nuc==0){
			cout<<"A";
		}
		if (nuc>=4){
			cout<<nuc<<endl;
			cout<<"WTF"<<endl;
		}
		anc>>=2;
	}
	cout<<endl;
}

string kmer2str(__uint128_t num, uint k = 31) {
	string res;
	Pow2<__uint128_t> anc(2 * (k - 1));
	for (uint64_t i(0); i < k; ++i) {
		__uint128_t nuc = num / anc;
		num             = num % anc;
		if (nuc == 3) {
			res += "G";
		}
		if (nuc == 2) {
			res += "T";
		}
		if (nuc == 1) {
			res += "C";
		}
		if (nuc == 0) {
			res += "A";
		}
		if (nuc >= 4) {
			cout << "WTF" << endl;
		}
		anc >>= 2;
	}
	return res;
}


uint64_t str2num(const string& str) {
	uint64_t res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}


__uint128_t rcb(const __uint128_t& in) {

	assume(k <= 64, "k=%u > 64", k);

	union kmer_u {
		__uint128_t k;
		__m128i m128i;
		uint64_t u64[2];
		uint8_t u8[16];
	};

	kmer_u res = {.k = in};

	// static_assert(sizeof(res) == sizeof(__uint128_t), "kmer sizeof mismatch");

	// Swap byte order

	kmer_u shuffidxs = {.u8 = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};

	res.m128i = _mm_shuffle_epi8(res.m128i, shuffidxs.m128i);

	// Swap nuc order in bytes

	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;

	const uint64_t c2 = 0x3333333333333333;

	for (uint64_t& x : res.u64) {

		x = ((x & c1) << 4) | ((x & (c1 << 4)) >> 4); // swap 2-nuc order in bytes

		x = ((x & c2) << 2) | ((x & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

		x ^= 0xaaaaaaaaaaaaaaaa; // Complement;
	}

	// Realign to the right

	res.m128i = mm_bitshift_right(res.m128i, 128 - 2 * k);

	return res.k;
}

uint64_t rcb(const uint64_t& in) {
	assume(k <= 32, "k=%u > 32", k);
	// Complement, swap byte order
	uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
	res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc
	// Realign to the right
	res >>= 64 - 2 * k;
	return res;
}

uint64_t rcbc(uint64_t in, uint64_t n) {
	assume(n <= 32, "n=%u > 32", n);
	// Complement, swap byte order
	uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
	res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

	// Realign to the right
	res >>= 64 - 2 * n;
	return res;
}

uint64_t canonize(uint64_t x, uint64_t n) {
	return min(x, rcbc(x, n));
}

// It's quite complex to bitshift mmx register without an immediate (constant) count
// See: https://stackoverflow.com/questions/34478328/the-best-way-to-shift-a-m128i
__m128i mm_bitshift_left(__m128i x, unsigned count) {
	assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_slli_si128(x, 8);
	if (count >= 64)                           //TODO: bench: Might be faster to skip this fast-path branch
		return _mm_slli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
	// else
	carry = _mm_srli_epi64(carry, 64 - count);

	x = _mm_slli_epi64(x, count);
	return _mm_or_si128(x, carry);
}

__m128i mm_bitshift_right(__m128i x, unsigned count) {
	assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_srli_si128(x, 8);
	if (count >= 64)
		return _mm_srli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
	// else
	carry = _mm_slli_epi64(carry, 64 - count);

	x = _mm_srli_epi64(x, count);
	return _mm_or_si128(x, carry);
}

uint64_t hash64shift(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

Pow2<uint64_t> minimizer_number(2 * minimizer_size);
/** Get the minimizer from a sequence and modify the position parameter.
	*/
uint64_t get_minimizer(uint64_t seq, uint64_t& min_position) {
	// print_kmer(seq,31);

	// Init with the first possible minimizer
	uint64_t mini, mmer;
	uint64_t fwd_mini = seq % minimizer_number;
	mini = mmer = canonize(fwd_mini, minimizer_size);

	uint64_t hash_mini = hash64shift(mmer);
	min_position = 0;
	// Search in all possible position (from 1) the minimizer
	uint64_t i(1), i_rc(k-minimizer_size);
	for (; i <= k - minimizer_size; i++, i_rc--) {
		seq >>= 2;
		fwd_mini = seq % minimizer_number;
		mmer = canonize(fwd_mini, minimizer_size);
		bool current_reversed = mmer != fwd_mini;
		uint64_t hash = (hash64shift(mmer));
		if (hash_mini > hash) {
			min_position = i;
			mini = mmer;
			hash_mini = hash;
		}
		else if (current_reversed and (hash_mini == hash) and i_rc < current_reversed) {
			min_position = i;
		}
	}
	return ((uint64_t)mini);
}


// ----- Kmer class -----

kmer_full::kmer_full(uint8_t minimizer_idx, uint64_t value, uint64_t reverse_comp_value) {
	this->minimizer_idx = minimizer_idx;
	this->kmer_s = value;
	this->kmer_rc = reverse_comp_value;
}

uint64_t kmer_full::get_minimizer() const {
	return (this->kmer_s >> (2*this->minimizer_idx)) & min_mask;
}
