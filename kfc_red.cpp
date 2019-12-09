#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include <stdint.h>
#include <stdio.h>
#include <atomic>
#include <mutex>
#include <omp.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <tmmintrin.h>
#include "robin_hood.h"


using namespace std;


bool check=false;
robin_hood::unordered_map<string,uint64_t> real_count;



// Represents the cardinality of a pow2 sized set. Allows div/mod arithmetic operations on indexes.
template<typename T>
struct Pow2 {
	Pow2(uint_fast8_t bits) : _bits(bits) {
		// assume(bits < CHAR_BIT*sizeof(T), "Pow2(%u > %u)", unsigned(bits), unsigned(CHAR_BIT*sizeof(T)));
	}

	uint_fast8_t bits() const { return _bits; }
	T value() const { return T(1) << _bits; }
	explicit operator T() const { return value(); }
	T max() const { return value() - T(1); }
	Pow2& operator=(const Pow2&) = default;

	Pow2():_bits(0){}
	Pow2(const Pow2&) = default;

	friend T operator*(const T& x, const Pow2& y) { return x << y._bits; }
	friend T& operator*=(T& x, const Pow2& y) { return x <<= y._bits; }
	friend T operator/(const T& x, const Pow2& y) { return x >> y._bits; }
	friend T& operator/=(T& x, const Pow2& y) { return x >>= y._bits; }
	friend T operator%(const T& x, const Pow2& y) { return x & y.max(); }
	friend T& operator%=(T& x, const Pow2& y) { return x &= y.max(); }
	Pow2& operator>>=(uint_fast8_t d) { _bits -= d; return *this; }
	Pow2& operator<<=(uint_fast8_t d) { _bits += d; return *this; }
	friend bool operator<(const T& x, const Pow2& y) { return x < y.value(); }
	friend bool operator<=(const T& x, const Pow2& y) { return x < y.value(); }
	friend T operator+(const T& x, const Pow2& y) { return x + y.value(); }
	friend T& operator+=(T& x, const Pow2& y) { return x += y.value(); }
	friend T operator-(const T& x, const Pow2& y) { return x - y.value(); }
	friend T& operator-=(T& x, const Pow2& y) { return x -= y.value(); }
private:
	uint_fast8_t _bits;
};



string kmer2str(__uint128_t num,uint k=31){
	string res;
	Pow2<__uint128_t> anc(2*(k-1));
	for(uint64_t i(0);i<k;++i){
		__uint128_t nuc=num/anc;
		num=num%anc;
		if(nuc==3){
			res+="G";
		}
		if(nuc==2){
			res+="T";
		}
		if(nuc==1){
			res+="C";
		}
		if(nuc==0){
			res+="A";
		}
		if (nuc>=4){
			cout<<"WTF"<<endl;
		}
		anc>>=2;
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



struct SKC{
  __uint128_t sk;
  vector<uint8_t> counts;
};



uint64_t k=(31);



uint64_t minimizer_size(9);


uint64_t minimizer_number((uint64_t)1<<(2*minimizer_size));



uint32_t revhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x2c1b3c6d;
	x = ( ( x >> 16 ) ^ x ) * 0x297a2d39;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



uint32_t unrevhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ( ( x >> 16 ) ^ x ) * 0x64ea2d65;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



uint64_t revhash ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



uint64_t unrevhash ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x );
	return x;
}




// It's quite complex to bitshift mmx register without an immediate (constant) count
// See: https://stackoverflow.com/questions/34478328/the-best-way-to-shift-a-m128i
  __m128i mm_bitshift_left(__m128i x, unsigned count)
{
	// assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_slli_si128(x, 8);
	if (count >= 64) //TODO: bench: Might be faster to skip this fast-path branch
		return _mm_slli_epi64(carry, count-64);  // the non-carry part is all zero, so return early
	// else
	carry = _mm_srli_epi64(carry, 64-count);

	x = _mm_slli_epi64(x, count);
	return _mm_or_si128(x, carry);
}



  __m128i mm_bitshift_right(__m128i x, unsigned count)
{
	// assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_srli_si128(x, 8);
	if (count >= 64)
		return _mm_srli_epi64(carry, count-64);  // the non-carry part is all zero, so return early
	// else
	carry = _mm_slli_epi64(carry, 64-count);

	x = _mm_srli_epi64(x, count);
	return _mm_or_si128(x, carry);
}




uint64_t rcbc(uint64_t in, uint64_t n){
  //assume(n <= 32, "n=%u > 32", n);
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



inline __uint128_t rcb(const __uint128_t& in){

	// assume(k <= 64, "k=%u > 64", k);

	union kmer_u { __uint128_t k; __m128i m128i; uint64_t u64[2]; uint8_t u8[16];};

	kmer_u res = { .k = in };

	// static_assert(sizeof(res) == sizeof(__uint128_t), "kmer sizeof mismatch");

	// Swap byte order

	kmer_u shuffidxs = { .u8 = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0} };

	res.m128i = _mm_shuffle_epi8 (res.m128i, shuffidxs.m128i);

	// Swap nuc order in bytes

	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;

	const uint64_t c2 = 0x3333333333333333;

	for(uint64_t& x : res.u64) {

		x = ((x & c1) << 4) | ((x & (c1 << 4)) >> 4); // swap 2-nuc order in bytes

		x = ((x & c2) << 2) | ((x & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

		x ^= 0xaaaaaaaaaaaaaaaa; // Complement;

	}

	// Realign to the right

	res.m128i = mm_bitshift_right(res.m128i, 128 - 2*k);

	return res.k;
}



inline uint64_t  rcb(const uint64_t& in){
    // assume(k <= 32, "k=%u > 32", k);
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




uint64_t canonize(uint64_t x,uint64_t n){
	return min(x,rcbc(x,n));
}



uint64_t get_minimizer(uint64_t seq){
	uint64_t mini,mmer;
	mmer=seq%minimizer_number;
	mini=mmer=canonize(mmer,minimizer_size);
	uint64_t hash_mini=unrevhash(mmer);
	for(uint64_t i(1);i<=k-minimizer_size;i++){
		seq>>=2;
		mmer=seq%minimizer_number;
		mmer=canonize(mmer,minimizer_size);
		// uint64_t hash = ((mmer));
		uint64_t hash = (unrevhash(mmer));
		if(hash_mini>hash){
			mini=mmer;
			hash_mini=hash;
		}
	}
	return ((uint64_t)mini)%minimizer_number;
}



uint64_t get_minimizer_pos(uint64_t seq,uint64_t& position){
	uint64_t mini,mmer;
	mmer=seq%minimizer_number;
	mini=mmer=canonize(mmer,minimizer_size);
	uint64_t hash_mini=unrevhash(mmer);
	position=0;
	for(uint64_t i(1);i<=k-minimizer_size;i++){
		seq>>=2;
		mmer=seq%minimizer_number;
		mmer=canonize(mmer,minimizer_size);
		uint64_t hash = (unrevhash(mmer));
		// uint64_t hash = ((mmer));
		if(hash_mini>hash){
			// cout<<i<<endl;
			// cout<<"kmer"<<kmer2str(seq,k)<<endl;
			// cout<<"min"<<kmer2str(hash,minimizer_size)<<endl;
			position=i;
			mini=mmer;
			hash_mini=hash;
		}
	}
	return ((uint64_t)mini)%minimizer_number;
}



uint64_t nuc2int(char c){
	return (c/2)%4;
}



uint64_t nuc2intrc(char c){
	return ((c/2)%4)^2;
}



uint64_t offsetUpdateAnchor((uint64_t)1<<(2*k));
uint64_t offsetUpdateAnchorMin((uint64_t)1<<(2*(minimizer_size)));



inline void updateK(uint64_t& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateAnchor;
}


inline void updateM(uint64_t& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateAnchorMin;
}



inline void updateRCK(uint64_t& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-2));
}



inline void updateRCM(uint64_t& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*minimizer_size-2));
}



char revCompChar(char c) {
switch (c) {
	case 'A': return 'T';
	case 'C': return 'G';
	case 'G': return 'C';
}
return 'A';
}



string revComp(const string& s){
string rc(s.size(),0);
for (int i((int)s.length() - 1); i >= 0; i--){
	rc[s.size()-1-i] = revCompChar(s[i]);
}
return rc;
}



string getCanonical(const string& str){
 return (min(str,revComp(str)));
}


uint64_t str2num(const string& str){
uint64_t res(0);
for(uint64_t i(0);i<str.size();i++){
	res<<=2;
	res+=(str[i]/2)%4;
}
return res;
}



struct kmer_full{
	uint64_t kmer_s;
	uint64_t kmer_rc;
};



int64_t kmer_in_super_kmer(const SKC& super_kmer,const kmer_full& kmer){
	__uint128_t superkmer=super_kmer.sk;
	uint64_t skc_size=super_kmer.counts.size();
	for(uint64_t i=0;i<skc_size;++i){
		if(kmer.kmer_s==(uint64_t)(superkmer&((1<<62)-1)) or kmer.kmer_rc==(uint64_t)(superkmer&((1<<62)-1))) {
			return i;
		}
		superkmer>>=2;
	}
	return -1;
}



uint64_t compaction(0);



int compact(SKC& super_kmer, kmer_full kmer){
	__uint128_t superkmer=super_kmer.sk;
	uint64_t sizesk=super_kmer.counts.size();
	__uint128_t end_sk=superkmer>>(2*sizesk);
	__uint128_t beg_k=kmer.kmer_s&((1<<62)-1);
	if(end_sk==beg_k){
		super_kmer.sk+=( ( (__uint128_t)(kmer.kmer_s>>60)<<(2*(31+sizesk)) ));
		super_kmer.counts.push_back(1);
		compaction++;
		return sizesk+1;
	}
	beg_k=kmer.kmer_rc&((1<<62)-1);
	if(end_sk==beg_k){
		super_kmer.sk+=( ( (__uint128_t)(kmer.kmer_rc>>60)<<(2*(31+sizesk)) ));
		super_kmer.counts.push_back(1);
		compaction++;
		return sizesk+1;
	}
	return -1;
}



void dump_count(const SKC& skc){
	string skmer(kmer2str(skc.sk,skc.counts.size()+30));
	for(uint64_t i(0);i<skc.counts.size();++i){
			cout<<getCanonical(skmer.substr(i,k))<<" "<<(uint64_t)skc.counts[skc.counts.size()-1-i]<<endl;
			if(check){
				if(real_count[getCanonical(skmer.substr(i,k))] != (uint64_t)skc.counts[skc.counts.size()-1-i]) {
					cerr<<"fail "<<getCanonical(skmer.substr(i,k))<<" "<<real_count[getCanonical(skmer.substr(i,k))]<<" "
					<<(uint64_t)skc.counts[skc.counts.size()-1-i]<<endl;
				}
			}
	}
}



void insert_kmers( vector<kmer_full>& kmers, vector<SKC>& skc){
	//FOREACh KMER
	uint64_t size_sk(kmers.size());
	uint64_t size_skc(skc.size());
	for(uint64_t ik=0;ik<size_sk;++ik){
		kmer_full kmer=kmers[ik];
		bool placed(false);
		// // FOREACh SUPERKMER

		for(uint64_t i=(0); i < size_skc and (not placed);i++){
			//IS IT HERE?
		 	int64_t pos=(kmer_in_super_kmer(skc[i],kmer));
	    if(pos>=0){
	      skc[i].counts[pos]++;
				placed=true;
	    }
	  }
		//FOREACh SUPERKMER
	  for(uint64_t i=(0); i < size_skc  and (not placed);i++){
			//CAN WE MERGE IT ?
		 	int64_t pos=(compact(skc[i],kmer));
	    if(pos>0){
				placed=true;
	    }
	  }
		if(not placed){
			skc.push_back({kmer.kmer_s,{1}});
		}
	}
	kmers.clear();
}



void dump_counting(const vector<vector<SKC>>& buckets){
	for(uint64_t i(0);i< buckets.size();++i) {
		for(uint64_t ii(0);ii<buckets[i].size();++ii){
			dump_count(buckets[i][ii]);
		}
	}
	if(check){cerr<<"The results were checked"<<endl;}
}


void dump_counting_stats(const vector<vector<SKC>>& buckets){
	for(uint64_t i(0);i< buckets.size();++i) {
		if(buckets[i].size()!=0)
			cout<<buckets[i].size()<<" ";
	}
	cout<<endl;
}



uint64_t nb_kmer_read(0);



void count_line(const string& line, vector<vector<SKC>>& buckets){
  vector<kmer_full> kmers;
	uint64_t seq=(str2num(line.substr(0,k))),rcSeq(rcb(seq)),canon(min(seq,rcSeq));
	uint64_t min_seq=(str2num(line.substr(k-minimizer_size,minimizer_size))),min_rcseq(rcbc(min_seq,minimizer_size)),min_canon(min(min_seq,min_rcseq));
	uint64_t position_min;
	uint64_t minimizer=get_minimizer_pos(rcSeq,position_min);
	kmers.push_back({seq,rcSeq});
	// insert_kmers(kmers,buckets[minimizer]);
	if(check){
		real_count[getCanonical(line.substr(0,k))]++;
	}
	uint64_t line_size=line.size();
  for(uint64_t i=0;i+k<line_size;++i){
		if(check){
			real_count[getCanonical(line.substr(i+1,k))]++;
		}
		updateK(seq,line[i+k]);
		updateRCK(rcSeq,line[i+k]);
		canon=(min(seq,rcSeq));
		updateM(min_seq,line[i+k]);
		updateRCM(min_rcseq,line[i+k]);
		min_canon=(min(min_seq,min_rcseq));

		//THE NEW mmer is a MINIMIZER
		if(unrevhash(min_canon)<unrevhash(minimizer)){
			insert_kmers(kmers,buckets[minimizer]);
			minimizer=min_canon;
			position_min=i+k-minimizer_size+1;
		}else{
			//the previous MINIMIZER is outdated
			if(i>=position_min){
				insert_kmers(kmers,buckets[minimizer]);
				minimizer=get_minimizer_pos(rcSeq,position_min);
				position_min+=i+1;
			}
		}
		kmers.push_back({seq,rcSeq});
  }
	insert_kmers(kmers,buckets[minimizer]);
}




void read_fasta_file(const string& filename,vector<vector<SKC>>& buckets){
  ifstream in(filename);
  string line;
	uint64_t line_count;
  while(in.good()){
    line=getLineFasta(&in);
    count_line(line,buckets);
		line_count++;
		if(line_count%1000==0){
			cerr<<"-"<<flush;
		}
  }
	cerr<<compaction<<endl;
	cerr<<line_count<<endl;
}




int main(int argc, char ** argv){
	if(argc<2){
		cout<<"[fasta file]"<<endl;
		exit(0);
	}
	if(argc>2){
		check=true;
	}
	vector<vector<SKC>> buckets(minimizer_number);
  read_fasta_file(argv[1],buckets);
	dump_counting(buckets);
	// dump_counting_stats(buckets);
  return 0;
}
