#ifndef __INDEX8MIN
#define __INDEX8MIN



#include <stdio.h>
#include <fstream>
#include <vector>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <algorithm>
#include "BBHash/BooPHF.h"


using namespace std;



typedef boomphf::SingleHashFunctor<uint64_t>  hasher;
typedef boomphf::mphf<  uint64_t, hasher > MPHF;






struct value{
	uint64_t kmer;
	uint8_t count;
};





class index_full{
public:
	uint32_t kmer_size;
	MPHF Hash;
	vector<value> Values;
	vector<uint64_t> weak_kmer_buffer;
	ofstream dump_weak;

	index_full(const vector<uint64_t>& V){
		kmer_size=(31);
		Values.resize(V.size());
		dump_weak.open(("weak_kmers"));
		auto data_iterator = boomphf::range(static_cast<const uint64_t*>(&((V)[0])), static_cast<const uint64_t*>((&(V)[0])+V.size()));
		Hash= boomphf::mphf<uint64_t,hasher>(V.size(),data_iterator,4,5,false);
		for(uint i(0);i<V.size();++i){
			int64_t pos(Hash.lookup(V[i]));
			Values[pos].kmer=V[i];
		}
	}

	void insert(uint64_t);
	void dump_counting();
	void insert_seq(const string&  read);
	void clear(bool=false);
	void print_kmer(uint64_t num);



};


uint64_t rcb(uint64_t min,uint n);




class index_min{
public:
	uint kmer_size;
	uint minimizer_size;

	vector<index_full> Index;

	index_min(vector<uint64_t>& V);
	void insert(uint64_t kmer);
	void dump_counting();

};



#endif
