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
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <tmmintrin.h>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <tmmintrin.h>
#include <math.h>
#include <algorithm>
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
  string sk;
  vector<uint8_t> counts;
};


struct SKC2{
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
	uint64_t mini,mmer,hash_mini(-1);
	mmer=seq%minimizer_number;
	mini=mmer=canonize(mmer,minimizer_size);
	for(uint64_t i(1);i<=k-minimizer_size;i++){
		seq>>=2;
		mmer=seq%minimizer_number;
		mmer=canonize(mmer,minimizer_size);
		uint64_t hash = (unrevhash(mmer));
		if(hash_mini>hash){
			mini=mmer;
			hash_mini=hash;
		}
	}
	return revhash((uint64_t)mini)%minimizer_number;
}



uint64_t nuc2int(char c){
	return (c/2)%4;
}



uint64_t nuc2intrc(char c){
	return ((c/2)%4)^2;
}



uint64_t offsetUpdateAnchor((uint64_t)1<<62);



inline void updateK(uint64_t& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateAnchor;
}



inline void updateRCK(uint64_t& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-2));
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



int64_t kmer_in_super_kmer(const string& super_kmer,const string& kmer){
	size_t found = super_kmer.find(kmer);
    if (found != string::npos){
			return found;
		}
		found = super_kmer.find(revComp(kmer));
	    if (found != string::npos){
				return found;
			}
		return -1;
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



int64_t kmer_in_super_kmer2(const SKC2& super_kmer,const uint64_t& kmer){
	__uint128_t superkmer=super_kmer.sk;
	uint64_t rc(rcb(kmer));
	for(uint64_t i=0;i<super_kmer.counts.size();++i){
		if(kmer==(uint64_t)(superkmer&((1<<62)-1)) or rc ==(uint64_t)(superkmer&((1<<62)-1))) {
			return i;
		}
		superkmer>>=2;
	}
	return -1;
}



int compact(string& super_kmer,const string& kmer){
	uint64_t ssk(super_kmer.size()),sk(kmer.size());
	if(sk==0){return -1;}
	if(ssk==0){
		super_kmer=kmer;
		return 0;
	}
	string rc_k(revComp(kmer)),end_sk(super_kmer.substr(ssk-k+1,k-1)), beg_k(kmer.substr(0,k-1));
	if(end_sk==beg_k){
		super_kmer=super_kmer+(kmer.substr(k-1));
		return super_kmer.size()-k+1;
	}
	string beg_rck(rc_k.substr(0,k-1));
	if(end_sk==beg_rck){
		super_kmer=super_kmer+(rc_k.substr(k-1));
		return super_kmer.size()-k+1;
	}
	return -1;
}


uint64_t compaction(0);


int compact2(SKC2& super_kmer, uint64_t kmer){
	// return -1;
	__uint128_t superkmer=super_kmer.sk;
	uint64_t sizesk=super_kmer.counts.size();
	// cout<<"super	"<<kmer2str(super_kmer.sk,sizesk+30)<<endl;
	__uint128_t end_sk=superkmer>>(2*sizesk);
	// cout<<"super	"<<kmer2str(end_sk,sizesk+30)<<endl;
	__uint128_t beg_k=kmer&((1<<62)-1);
	// cout<<"beg kmer	"<<kmer2str(beg_k,33)<<endl;
	// cout<<"end skmer	"<<kmer2str(end_sk,33)<<endl;
	if(end_sk==beg_k){
		// cout<<"compact1"<<endl;
		// cout<<"kmer	"<<kmer2str(kmer,31)<<endl;
		// cout<<"super	"<<kmer2str(super_kmer.sk,30+sizesk)<<endl;
			// cout<<"kmer	"<<kmer2str(kmer>>60,31)<<endl;
		super_kmer.sk+=( ( (__uint128_t)(kmer>>60)<<(2*(31+sizesk)) ));
		// cout<<"super	"<<kmer2str(super_kmer.sk,31+sizesk)<<endl;
		// cin.get();
		super_kmer.counts.push_back(1);
		compaction++;
		return sizesk+1;
	}
	kmer=rcb(kmer);
	beg_k=kmer&((1<<62)-1);
	if(end_sk==beg_k){
		// cout<<"compact2"<<endl;
		// cout<<"kmer	"<<kmer2str(kmer,31)<<endl;
		// cout<<"super	"<<kmer2str(super_kmer.sk,sizesk+30)<<endl;
		// super_kmer.sk<<=2;
		// cout<<"super	"<<kmer2str(super_kmer.sk,sizesk+31)<<endl;
		super_kmer.sk+=( ( (__uint128_t)(kmer>>60)<<(2*(31+sizesk)) ));
		// cout<<"super	"<<kmer2str(super_kmer.sk,sizesk+31)<<endl;
		super_kmer.counts.push_back(1);
		compaction++;
		return sizesk+1;
	}
	return -1;
}



void dump_count(const SKC& skc){
	for(uint64_t i(0);i<skc.counts.size();++i){
		cout<<">"<<(uint64_t)skc.counts[i]<<"\n"<<getCanonical(skc.sk.substr(i,k))<<"\n";
		// cout<<getCanonical(skc.sk.substr(i,k))<<" "<<(uint64_t)skc.counts[i]<<"\n";
	}
}


void dump_count2(const SKC2& skc){
	string skmer(kmer2str(skc.sk,skc.counts.size()+30));
	// cout<<skc.counts.size()<<endl;
	// cout<<skmer.size()<<endl;
	// reverse(skc.counts.begin(),skc.counts.end());
	for(uint64_t i(0);i<skc.counts.size();++i){
			// cout<<">"<<(uint64_t)skc.counts[i]<<"\n"<<getCanonical(skmer.substr(i,k))<<"\n";
			// cout<<getCanonical(skmer.substr(i,k))<<" "<<(uint64_t)skc.counts[i]<<endl;
			// if( skmer.substr(i,k)=="CTCATTCCCCAGTAGTACTGCAAGAGGTTCC") {
			// 	cout<<skmer<<endl;
			// 	for(uint64_t ii(0);ii<skc.counts.size();++ii){
			// 		cout<<(uint64_t)skc.counts[ii]<<endl;
			// 	}
			//
			// }
			if(check){
				if(real_count[getCanonical(skmer.substr(i,k))] != (uint64_t)skc.counts[skc.counts.size()-1-i]) {
					cerr<<"fail "<<getCanonical(skmer.substr(i,k))<<" "<<real_count[getCanonical(skmer.substr(i,k))]<<" "
					<<(uint64_t)skc.counts[skc.counts.size()-1-i]<<endl;
				}
			}

			// 	// cin.get();
			// }
	}
}



void insert_kmer(const string& str_kmer, vector<vector<SKC>>& skc){
  uint64_t min(get_minimizer(str2num(str_kmer)));
	for(uint64_t i(0); i < skc[min].size();i++){
	 	int64_t pos(kmer_in_super_kmer(skc[min][i].sk,str_kmer));
    if(pos>=0){
      skc[min][i].counts[pos]++;
			return;
    }
  }

  // for(uint64_t i(0); i < skc[min].size();i++){
	//  	int64_t pos(compact(skc[min][i].sk,str_kmer));
  //   if(pos>	0){
  //     if(skc[min][i].counts.size()<pos){
  //       skc[min][i].counts.resize(pos,0);
  //     }
  //     skc[min][i].counts[pos-1]++;
	// 		return;
  //   }
  // }
	skc[min].push_back({str_kmer,{1}});
}



void insert_kmer2(const uint64_t& kmer, vector<vector<SKC2>>& skc){
  uint64_t min=(get_minimizer(kmer));
	for(uint64_t i=(0); i < skc[min].size();i++){
		// cout<<"searchsk"<<endl;
	 	int64_t pos=(kmer_in_super_kmer2(skc[min][i],kmer));
    if(pos>=0){
      skc[min][i].counts[pos]++;
			return;
    }
  }

  for(uint64_t i=(0); i < skc[min].size();i++){
			// cout<<"search compact"<<endl;
	 	int64_t pos=(compact2(skc[min][i],kmer));
    if(pos>0){
			return;
    }
  }
	skc[min].push_back({kmer,{1}});
}



void dump_counting(vector<vector<SKC>>& buckets){
	for(uint64_t i(0);i< buckets.size();++i) {
		for(uint64_t ii(0);ii<buckets[i].size();++ii){
			dump_count(buckets[i][ii]);
		}
	}
}



void dump_counting2(vector<vector<SKC2>>& buckets){
	for(uint64_t i(0);i< buckets.size();++i) {
		for(uint64_t ii(0);ii<buckets[i].size();++ii){
			dump_count2(buckets[i][ii]);
		}
	}
	if(check){cout<<"The results were checked"<<endl;}
}



uint64_t nb_kmer_read(0);



void count_line(const string& line, vector<vector<SKC>>& buckets){
  string str_kmer;
  for(uint64_t i=0;i+k<=line.size();++i){
    str_kmer=getCanonical(line.substr(i,k));
    insert_kmer(str_kmer,buckets);
		nb_kmer_read++;
  }
}



void count_line2(const string& line, vector<vector<SKC2>>& buckets){
  uint64_t kmer;
	uint64_t seq=(str2num(line.substr(0,k))),rcSeq(rcb(seq)),canon(min(seq,rcSeq));
	insert_kmer2(canon,buckets);
	// cout<<kmer2str(canon)<<endl;;
	if(check){
		real_count[getCanonical(line.substr(0,k))]++;
	}
  for(uint64_t i=0;i+k<line.size();++i){
		if(check){
			real_count[getCanonical(line.substr(i+1,k))]++;
		}
		updateK(seq,line[i+k]);
		updateRCK(rcSeq,line[i+k]);
		canon=(min(seq,rcSeq));
		// cout<<kmer2str(canon)<<endl;;
		insert_kmer2(canon,buckets);
  }
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
}



void read_fasta_file2(const string& filename,vector<vector<SKC2>>& buckets){
  ifstream in(filename);
  string line;
	uint64_t line_count;
  while(in.good()){
    line=getLineFasta(&in);
    count_line2(line,buckets);
		line_count++;
		if(line_count%1000==0){
			cerr<<"-"<<flush;
		}
  }
	cout<<compaction<<endl;
}




int main(int argc, char ** argv){
	if(argc<2){
		cout<<"[fasta file]"<<endl;
		exit(0);
	}
	if(false){
		vector<vector<SKC>> buckets(minimizer_number);
	  read_fasta_file(argv[1],buckets);
		dump_counting(buckets);
	}else{
		vector<vector<SKC2>> buckets(minimizer_number);
	  read_fasta_file2(argv[1],buckets);
		dump_counting2(buckets);
	}
  return 0;
}
