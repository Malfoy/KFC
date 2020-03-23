#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <fstream>
#include <gatbl/common.hpp>
#include <gatbl/kmer.hpp>
#include <iostream>
#include <math.h>
#include <mutex>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <tmmintrin.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include "robin_hood.h"

using namespace std;

bool check = false;
robin_hood::unordered_map<string, uint64_t> real_count;

uint64_t k = (31);
uint64_t minimizer_size(4);
uint64_t subminimizer_size(minimizer_size-1);
uint64_t nb_kmer_read(0);
uint64_t line_count(0);

vector<omp_lock_t> MutexWall(4096);

// Represents the cardinality of a pow2 sized set. Allows div/mod arithmetic operations on indexes.
template <typename T>
struct Pow2 {
	Pow2(uint_fast8_t bits)
	    : _bits(bits) {
		assume(bits < 8 * sizeof(T), "Pow2(%u > %u)", unsigned(bits), unsigned(8 * sizeof(T)));
	}

	uint_fast8_t bits() const { return _bits; }
	T value() const { return T(1) << _bits; }
	explicit operator T() const { return value(); }
	T max() const { return value() - T(1); }
	Pow2& operator=(const Pow2&) = default;

	Pow2()
	    : _bits(0) {}
	Pow2(const Pow2&) = default;

	friend T operator*(const T& x, const Pow2& y) { return x << y._bits; }
	friend T& operator*=(T& x, const Pow2& y) { return x <<= y._bits; }
	friend T operator/(const T& x, const Pow2& y) { return x >> y._bits; }
	friend T& operator/=(T& x, const Pow2& y) { return x >>= y._bits; }
	friend T operator%(const T& x, const Pow2& y) { return x & y.max(); }
	friend T& operator%=(T& x, const Pow2& y) { return x &= y.max(); }
	Pow2& operator>>=(uint_fast8_t d) {
		_bits -= d;
		return *this;
	}
	Pow2& operator<<=(uint_fast8_t d) {
		_bits += d;
		return *this;
	}
	friend bool operator<(const T& x, const Pow2& y) { return x < y.value(); }
	friend bool operator<=(const T& x, const Pow2& y) { return x < y.value(); }
	friend T operator+(const T& x, const Pow2& y) { return x + y.value(); }
	friend T& operator+=(T& x, const Pow2& y) { return x += y.value(); }
	friend T operator-(const T& x, const Pow2& y) { return x - y.value(); }
	friend T& operator-=(T& x, const Pow2& y) { return x -= y.value(); }

	private:
	uint_fast8_t _bits;
};

Pow2<uint64_t> minimizer_number(2 * minimizer_size);
Pow2<uint64_t> bucket_number(2 * subminimizer_size);

Pow2<uint64_t> offsetUpdateAnchor(2 * k);
Pow2<uint64_t> offsetUpdateAnchorMin(2 * minimizer_size);

struct SKC {
	__uint128_t sk;
	uint8_t size;
	uint8_t counts[15] = {0};
	uint8_t minimizer_idx;
	SKC(uint64_t kmer, uint8_t mini_idx) {
		sk        = kmer;
		size      = 1;
		counts[0] = 1;
		minimizer_idx = mini_idx;
	};
};

struct kmer_full {
	uint8_t minimizer_idx;
	uint64_t kmer_s;
	uint64_t kmer_rc;
};

//START LOW LEVEL FUNCTIONS

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

string getLineFasta(ifstream* in) {
	string line, result;
	getline(*in, line);
	char c = static_cast<char>(in->peek());
	while (c != '>' and c != EOF) {
		getline(*in, line);
		result += line;
		c = static_cast<char>(in->peek());
	}
	return result;
}

uint32_t revhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x2c1b3c6d;
	x = ((x >> 16) ^ x) * 0x297a2d39;
	x = ((x >> 16) ^ x);
	return x;
}

uint32_t unrevhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ((x >> 16) ^ x) * 0x64ea2d65;
	x = ((x >> 16) ^ x);
	return x;
}

uint64_t revhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x);
	return x;
}

uint64_t unrevhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
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

inline __uint128_t rcb(const __uint128_t& in) {

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

inline uint64_t rcb(const uint64_t& in) {
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

uint64_t canonize(uint64_t x, uint64_t n) {
	return min(x, rcbc(x, n));
}


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

uint64_t get_minimizer_pos(uint64_t seq, uint64_t& position) {
	print_kmer(seq,31);
	uint64_t mini, mmer;
	mmer = seq % minimizer_number;
	mini = mmer        = canonize(mmer, minimizer_size);
	uint64_t hash_mini = hash64shift(mmer);
	position           = 0;
	for (uint64_t i(1); i <= k - minimizer_size; i++) {
		seq >>= 2;
		mmer          = seq % minimizer_number;
		mmer          = canonize(mmer, minimizer_size);
		uint64_t hash = (hash64shift(mmer));
		if (hash_mini > hash) {
			position  = i;
			mini      = mmer;
			hash_mini = hash;
		}
	}
	cout<<"get min pos"<<endl;
	print_kmer(mini,minimizer_size);

	cout<<position<<endl;
	// return revhash((uint64_t)mini)%minimizer_number;
	return ((uint64_t)mini);
}

inline uint64_t nuc2int(char c) {
	return (c / 2) % 4;
}

inline uint64_t nuc2intrc(char c) {
	return ((c / 2) % 4) ^ 2;
}

inline void updateK(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchor;
}

inline void add_nuc_superkmer(SKC& min, char nuc) {
	min.sk <<= 2;
	min.sk += nuc2int(nuc);
	min.counts[min.size++] = 0;
}

inline void updateM(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchorMin;
}

inline void updateRCK(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * k - 2));
}

inline void updateRCM(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * minimizer_size - 2));
}

char revCompChar(char c) {
	switch (c) {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
	}
	return 'A';
}

string revComp(const string& s) {
	string rc(s.size(), 0);
	for (int i((int)s.length() - 1); i >= 0; i--) {
		rc[s.size() - 1 - i] = revCompChar(s[i]);
	}
	return rc;
}

string getCanonical(const string& str) {
	return (min(str, revComp(str)));
}

uint64_t str2num(const string& str) {
	uint64_t res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}

string intToString(uint64_t n) {
	if (n < 1000) {
		return to_string(n);
	}
	string end(to_string(n % 1000));
	if (end.size() == 3) {
		return intToString(n / 1000) + "," + end;
	}
	if (end.size() == 2) {
		return intToString(n / 1000) + ",0" + end;
	}
	return intToString(n / 1000) + ",00" + end;
}

//START HIGH LEVEL FUNCTIONS

void dump_count(const SKC& skc) {
	string skmer(kmer2str(skc.sk, skc.size + 30));
	for (uint64_t i(0); i < skc.size; ++i) {
		cout << getCanonical(skmer.substr(i, k)) << " " << (uint64_t)skc.counts[skc.size - 1 - i] << endl;
		if (check) {
			if (real_count[getCanonical(skmer.substr(i, k))] != (uint64_t)skc.counts[skc.size - 1 - i]) {
				cout << i << " " << skc.size - 1 - i << endl;
				cerr << "fail " << (skmer.substr(i, k)) << " " << revComp(skmer.substr(i, k)) << " real: " << real_count[getCanonical(skmer.substr(i, k))] << " output:"
				     << (uint64_t)skc.counts[skc.size - 1 - i] << endl;
				cin.get();
			}
		}
	}
}

void dump_counting(const vector<vector<SKC> >& buckets) {
	for (uint64_t i(0); i < buckets.size(); ++i) {
		for (uint64_t ii(0); ii < buckets[i].size(); ++ii) {
			dump_count(buckets[i][ii]);
		}
	}
	if (check) {
		cerr << "The results were checked" << endl;
	}
}

void dump_stats(const vector<vector<SKC> >& buckets) {
	uint64_t total_super_kmers(0);
	uint64_t total_kmers(0);
	uint64_t non_null_buckets(0);
	uint64_t null_buckets(0);
	uint64_t largest_bucket(0);

	//FOREACH BUCKETS
	for (uint64_t i(0); i < buckets.size(); ++i) {
		if (buckets[i].size() != 0) {
			largest_bucket=max(largest_bucket,(uint64_t)buckets[i].size());
			non_null_buckets++;
			total_super_kmers += buckets[i].size();
			for (uint64_t j(0); j < buckets[i].size(); ++j) {
				total_kmers += buckets[i][j].size;
			}
		} else {
			null_buckets++;
		}
	}
	cout << endl;
	cout << "Empty buckets:	" << intToString(null_buckets) << endl;
	cout << "Useful buckets:	" << intToString(non_null_buckets) << endl;
	cout << "#Superkmer:	" << intToString(total_super_kmers) << endl;
	cout << "#kmer:	" << intToString(total_kmers) << endl;
	cout << "super_kmer per useful buckets:	" << intToString(total_super_kmers / non_null_buckets) << endl;
	cout << "kmer per useful buckets:	" << intToString(total_kmers / non_null_buckets) << endl;
	cout << "kmer per super_kmer:	" << intToString(total_kmers / total_super_kmers) << endl;
	cout<<"Largest_bucket:	"<<intToString(largest_bucket)<<endl;
}





int64_t kmer_in_super_kmer(const SKC& super_kmer, const kmer_full& kmer) {
	__uint128_t superkmer = super_kmer.sk;
	uint64_t skc_size     = super_kmer.size;
	for (uint64_t i = 0; i < skc_size; ++i) {
		if (kmer.kmer_s == (uint64_t)(superkmer & (((uint64_t)1 << 62) - 1)) or kmer.kmer_rc == (uint64_t)(superkmer & (((uint64_t)1 << 62) - 1))) {
			return i;
		}
		superkmer >>= 2;
	}
	return -1;
}

uint64_t mask = ((uint64_t)1 << 62) - 1;
/** Return true if the kmer is present in the superkmer.
  * The technic used is based on the knoledge of the minimizer idx.
  * The superkmer is aligned to the kmer using the minimizer as anchor.
  * Then a comparison is done using a xor and the result is returned
  * @return the position of the kmer in SKC, -1 if not found.
  */
int64_t kmer_in_super_kmer_short(const SKC& super_kmer,const kmer_full& kmer){
	int64_t start_idx = (int64_t)super_kmer.minimizer_idx - (int64_t)kmer.minimizer_idx;
	uint64_t sub_sk = (super_kmer.sk >> (2 * start_idx)) & mask;
	cout<<"GO KMER IN SUPERKMER"<<endl;
	if (sub_sk == kmer.kmer_s or sub_sk == kmer.kmer_rc) {
		cout<<"SUCCESS!!!!!!!!!!!!!!!!"<<endl;

		cout<<"idx_super_kmer"<<(int64_t)super_kmer.minimizer_idx<<endl;
		cout<<"idx kmer"<<(int64_t)kmer.minimizer_idx<<endl;
		cout<<"start_idx"<<start_idx<<endl;
		print_kmer(kmer.kmer_s,31);
		// print_kmer(kmer.kmer_s,31);
		print_kmer(super_kmer.sk,33	);
		cout<<endl;
		return start_idx;
	}
	cout<<"FAIL!"<<endl;
	cout<<"idx_super_kmer"<<(int64_t)super_kmer.minimizer_idx<<endl;
	cout<<"idx kmer"<<(int64_t)kmer.minimizer_idx<<endl;
	cout<<"start_idx"<<start_idx<<endl;
	print_kmer(kmer.kmer_s,31);
	// print_kmer(kmer.kmer_s,31);
	print_kmer(super_kmer.sk,33	);
		cout<<endl;
	// print_kmer(kmer.kmer_s,31);
	// If one of the two is zero, return the position
	return -1;
}

vector<bool> kmers_in_super_kmer(vector<SKC>& super_kmers, const vector<kmer_full>& kmers) {
	uint64_t size_bucket = (super_kmers.size());
	uint64_t buffer_size = kmers.size();
	vector<bool> res(buffer_size, false);
	for (uint64_t b = (0); b < size_bucket; b++) {
		__uint128_t superkmer = super_kmers[b].sk;
		uint64_t sk_size      = super_kmers[b].size;
		for (uint64_t i = 0; i < sk_size; ++i) {
			for (uint64_t j = 0; j < buffer_size; ++j) {
				if (kmers[j].kmer_s == (uint64_t)(superkmer & (((uint64_t)1 << 62) - 1)) or kmers[j].kmer_rc == (uint64_t)(superkmer & (((uint64_t)1 << 62) - 1))) {
					++super_kmers[b].counts[i];
					res[j] = true;
				}
			}
			superkmer >>= 2;
		}
	}
	return res;
}

int compact(SKC& super_kmer, kmer_full kmer) {
	if (super_kmer.size >= 13) {
		return -1;
	}
	uint64_t superkmer = super_kmer.sk;
	uint64_t end_sk    = superkmer;
	end_sk &= (((uint64_t)1 << 60) - 1);
	uint64_t beg_k = kmer.kmer_s >> 2;
	if (end_sk == beg_k) {
		super_kmer.sk <<= 2;
		super_kmer.sk += (((kmer.kmer_s % 4)));
		super_kmer.counts[super_kmer.size++] = 1;
		super_kmer.minimizer_idx++;
		return super_kmer.size;
	}
	// beg_k = kmer.kmer_rc >> 2;
	// if (end_sk == beg_k) {
	// 	// super_kmer.sk <<= 2;//TODO IDONOTKNOW
	// 	super_kmer.sk += (((kmer.kmer_rc % 4)));
	// 	super_kmer.counts[super_kmer.size++] = 1;
	// 	return super_kmer.size;
	// }
	return -1;
}


void insert_kmers2(vector<kmer_full>& kmers, vector<SKC>& skc) {
	uint64_t size_sk(kmers.size());
	uint64_t size_skc(skc.size());
	bool init(false);
	//FOREACH KMER
	for (uint64_t ik = 0; ik < size_sk; ++ik) {
		kmer_full kmer = kmers[ik];
		bool placed(false);
		// // FOREACh SUPERKMER
		for (uint64_t i = (0); i < size_skc and (not placed); i++) {
			//IS IT HERE?
			// int64_t pos = (kmer_in_super_kmer_short(skc[i], {kmer}));
			int64_t pos = (kmer_in_super_kmer(skc[i], {kmer}));
			if (pos >= 0) {
				skc[i].counts[pos]++;
				placed = true;
			}
		}
		//TRY TO COMPACT IT TO THE LAST ONE
		if (init and not placed) {
			int64_t pos = (compact(skc[size_skc - 1], kmer));
			if (pos > 0) {
				placed = true;
			}
		}
		//MAKE IT A NEW ONE
		if (not placed) {
			skc.push_back(SKC(kmer.kmer_s, kmer.minimizer_idx));
			init = true;
			size_skc++;
		}
	}
	kmers.clear();
}

void insert_kmers(vector<kmer_full>& kmers, vector<SKC>& skc) {
	uint64_t size_sk(kmers.size());
	uint64_t size_skc(skc.size());
	uint64_t elementfound(0);
	vector<bool> placed(size_sk, false);
	// // FOREACh SUPERKMER
	for (uint64_t i = 0; i < size_skc; i++) {
		auto localsk = skc[i];
		//FOREACH KMER
		for (uint64_t ik = 0; ik < size_sk; ++ik) {
			if (not placed[ik]) {
				kmer_full kmer = kmers[ik];
				//IS IT HERE?
				int64_t pos = (kmer_in_super_kmer_short(localsk, kmer));
				// int64_t pos = (kmer_in_super_kmer(localsk, kmer));
				if (pos >= 0) {
					++skc[i].counts[pos];
					placed[ik] = true;
					if (++elementfound == size_sk) {
						kmers.clear();
						return;
					}
				}
			}
		}
	}
	bool init(false); //TODO  WHY DO WE NEED THIS MECANISM???
	//FOREACH KMER
	for (uint64_t ik = 0; ik < size_sk; ++ik) {
		if (not placed[ik]) {
			kmer_full kmer = kmers[ik];
			//TRY TO COMPACT IT TO THE LAST ONE
			if (init) {
				int64_t pos = (compact(skc[size_skc - 1], kmer));
				if (pos > 0) {
					if (++elementfound == size_sk) {
						kmers.clear();
						return;
					}
				} else {
					//MAKE IT A NEW ONE
					skc.push_back(SKC(kmer.kmer_s, kmer.minimizer_idx));
					++size_skc;
					if (++elementfound == size_sk) {
						kmers.clear();
						return;
					}
					init = true;
				}
			} else {
				//MAKE IT A NEW ONE
				skc.push_back(SKC(kmer.kmer_s, kmer.minimizer_idx));
				++size_skc;
				if (++elementfound == size_sk) {
					kmers.clear();
					return;
				}
				init = true;
			}
		} else {
			init = false;
		}
	}
	kmers.clear();
}



void count_line(const string& line, vector<vector<SKC> >& buckets) {
	if (line.size() < k) {
		return;
	}
	vector<kmer_full> kmers;
	uint64_t seq     = (str2num(line.substr(0, k))), rcSeq(rcb(seq));
	uint64_t min_seq = (str2num(line.substr(k - minimizer_size, minimizer_size))), min_rcseq(rcbc(min_seq, minimizer_size)), min_canon(min(min_seq, min_rcseq));
	uint64_t position_min;
	uint64_t position_minimizer_in_kmer;
	uint64_t minimizer = get_minimizer_pos(seq, position_min);
	position_minimizer_in_kmer=position_min;
	uint64_t hash_mini = hash64shift(minimizer);
	kmers.push_back({position_min, seq, rcSeq});
	if (check) {
		real_count[getCanonical(line.substr(0, k))]++;
	}
	uint64_t line_size = line.size();
	for (uint64_t i = 0; i + k < line_size; ++i) {
		if (check) {
			real_count[getCanonical(line.substr(i + 1, k))]++;
		}
		updateK(seq, line[i + k]);
		updateRCK(rcSeq, line[i + k]);
		updateM(min_seq, line[i + k]);
		updateRCM(min_rcseq, line[i + k]);
		min_canon = (min(min_seq, min_rcseq));

		//THE NEW mmer is a MINIMIZER
		uint64_t new_hash = (hash64shift(min_canon));
		if (new_hash < hash_mini) {
			uint64_t bucketindice = revhash(hash_mini) % bucket_number;
			omp_set_lock(&MutexWall[bucketindice % 4096]);
			insert_kmers(kmers, buckets[bucketindice]);
			omp_unset_lock(&MutexWall[bucketindice % 4096]);
			minimizer    = (min_canon);
			hash_mini    = new_hash;
			position_min = i + k - minimizer_size + 1;
			position_minimizer_in_kmer=0;
			cout<<"new min go 0"<<endl;
		} else {
			//the previous MINIMIZER is outdated
			if (i >= position_min) {
				uint64_t bucketindice = revhash(hash_mini) % bucket_number;
				omp_set_lock(&MutexWall[bucketindice % 4096]);
				insert_kmers(kmers, buckets[bucketindice]);
				omp_unset_lock(&MutexWall[bucketindice % 4096]);
				minimizer = get_minimizer_pos(seq, position_min);
				cout<<"new min get min pos"<<endl;
				cout<<position_minimizer_in_kmer<<endl;
				position_minimizer_in_kmer=position_min;
				hash_mini = hash64shift(minimizer);
				position_min += i + 1;
			}
		}
		// if(position_min+kmers.size()>31){
		// 	cout<<"WTF"<<endl;
		// }
		position_minimizer_in_kmer++;
		cout<<"add"<<position_minimizer_in_kmer<<endl;
		kmers.push_back({position_minimizer_in_kmer, seq, rcSeq});
	}
	uint64_t bucketindice = revhash(hash_mini) % bucket_number;
	omp_set_lock(&MutexWall[bucketindice % 4096]);
	insert_kmers(kmers, buckets[bucketindice]);
	omp_unset_lock(&MutexWall[bucketindice % 4096]);
}

void read_fasta_file(const string& filename, vector<vector<SKC> >& buckets) {
	for (uint64_t i(0); i < 4096; ++i) {
		omp_init_lock(&MutexWall[i]);
	}
	ifstream in(filename);
	uint8_t nb_core(4);
	if(check){
		nb_core=1;
	}
#pragma omp parallel num_threads(nb_core)
	{
		string line;
		while (in.good()) {
#pragma omp critical
			{
				line = getLineFasta(&in);
			}
			count_line(line, buckets);
			line_count++;
			if (line_count % 100000 == 0) {
				cerr << "-" << flush;
			}
		}
	}
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "[fasta file]" << endl;
		cout << "0 DEFAULT, OUTPUT THE COUNT WITH NO CHECKING" << endl;
		cout << "1 PERFORMANCE MODE NO COUNT OUTPUT" << endl;
		cout << "2 DEBUG MODE COUNT AND CHECK WITH HASHTABLE" << endl;
		exit(0);
	}
	int mode(0);
	if (argc > 2) {
		mode = stoi(argv[2]);
	}

	if (mode > 1) {
		check = true;
	}
	vector<vector<SKC> > buckets(bucket_number.value());
	read_fasta_file(argv[1], buckets);
	cout << endl;
	if (mode % 2 == 0) {
		dump_counting(buckets);
	}
	dump_stats(buckets);
	exit(0);
}
