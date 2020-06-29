#include <iostream>
#include <cstdint>
#include "Kmers.hpp"



#ifndef SKC_H
#define SKC_H



class SKC {
public:
	uint64_t sk;//8B
	uint8_t size;//1B
	uint8_t minimizer_idx;//1B
	uint8_t counts[6] = {0};//22B
	//TOTAL 16

	/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int 64 used to represent the binary kmer.
	* @param mini_idx The minimizer position in the kmer.
	*/
	SKC(const uint64_t kmer, const uint8_t mini_idx);
	/** Look for the minimizer between the SKC and both of the fwd and rev kmer. Call the right compact function if needed.*/
	bool add_kmer_old(const kmer_full& kmer);
	bool add_kmer(const kmer_full& kmer);
	bool operator< (const  SKC& str);
	friend std::ostream & operator << (std::ostream& out, const SKC& skc);
	bool suffix_is_prefix(const kmer_full&);
	void print_count(const string& out,const string minimizer) const;
	bool  operator < (const  SKC& str) const {
		uint64_t s1((get_suffix()));
		uint64_t s2((str.get_suffix()));
		if(minimizer_idx>=str.minimizer_idx){
			s1>>=(2*(minimizer_idx-str.minimizer_idx));
		}else{
			s2>>=(2*(str.minimizer_idx-minimizer_idx));
		}
		if(s1==s2){
			return minimizer_idx < str.minimizer_idx;
		}
		return  s1<s2;
	}
	bool compact_right(const uint64_t);
	bool compact_left(const uint64_t);
	bool is_present(uint64_t, uint64_t);
	bool is_present_brutforce(kmer_full kmer, uint8_t & mini_k_idx);
	bool compact_right(const kmer_full& kmf);
	const bool is_present(kmer_full kmf);
	uint64_t get_prefix()  const;
	uint64_t get_suffix() const;
};



#endif
