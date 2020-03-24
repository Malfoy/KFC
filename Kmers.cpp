
#include <iostream>
#include "Kmers.hpp"

using namespace std;

uint64_t k = 31;
uint64_t k_mask = (((uint64_t)1) << (2*k)) - 1;
uint64_t minimizer_size = 7;
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


// ----- Kmer class -----

kmer_full::kmer_full(uint8_t minimizer_idx, uint64_t value, uint64_t reverse_comp_value) {
	this->minimizer_idx = minimizer_idx;
	this->kmer_s = value;
	this->kmer_rc = reverse_comp_value;
}
