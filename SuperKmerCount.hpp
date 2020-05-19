#include <iostream>
#include <cstdint>
#include "Kmers.hpp"



#ifndef SKC_H
#define SKC_H



class SKC {
public:
	uint64_t sk;
	uint8_t counts[12] = {0};
	uint8_t size;
	uint8_t minimizer_idx;
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

private:
	bool compact_right(const uint64_t);
	bool compact_left(const uint64_t);
	bool is_present(uint64_t, uint64_t);
	bool is_present_brutforce(kmer_full kmer, uint8_t & mini_k_idx);
	bool compact_right(const kmer_full kmf);
	bool is_present(kmer_full kmf);
	 uint64_t get_prefix()  const;
	 uint64_t get_suffix() const;
};

#endif
