#include "SuperKmerCount.hpp"
#include "Kmers.hpp"



uint64_t nb_superkmer(1);
uint64_t nb_kmer_read(0);



using namespace std;



/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int 64 used to represent the binary kmer.
	* @param mini_idx The minimizer position in the kmer.
	*/
SKC::SKC(const uint64_t kmer, const uint8_t mini_idx) {
	this->sk            = kmer;
	this->size          = 1;
	this->counts[0]     = 1;
	this->minimizer_idx = mini_idx;
};



//THE PREFIX IS AT RIGHT!
 uint64_t SKC::get_prefix()  const {
	return sk>>(2*minimizer_idx);
}




//THE SUFFIX IS AT LEFT!
 uint64_t SKC::get_suffix() const  {
	return sk%((uint64_t)1<<(2*minimizer_idx));
}



bool  SKC::operator< (const  SKC& str){
	return  (reversebits(get_suffix()) < reversebits(str.get_suffix()));
}



/** Used to compact a new nucleotide from a kmer on the right of the superkmer.
	* @param kmer The binary representation of the kmer to compact on the right. Must be on the same strand than the superkmer.
	* @return True if the kmer is inserted false otherwise.
	*/
bool SKC::compact_right(const uint64_t kmer_val) {
	uint64_t end_sk = this->sk;
	end_sk &= (((uint64_t)1 << (2 * (k - 1))) - 1);
	uint64_t begin_kmer = kmer_val >> 2;
	if (end_sk == begin_kmer) {
		this->sk <<= 2;
		this->sk += kmer_val % 4;
		this->size ++;
		this->minimizer_idx ++;
		this->counts[this->size-1]               = 1;
		return true;
	}
	return false;
}



bool SKC::compact_right(const kmer_full kmf) {
	uint64_t prefix(get_prefix());
	uint64_t suffix(get_suffix());
	if((kmf.suffix>>2)==suffix){
		if(kmf.prefix==(prefix%((uint64_t)1<<(2*(31-minimizer_size-kmf.minimizer_idx))))) {
			sk<<=2;
			sk += (kmf.suffix % 4);
			size++;
			minimizer_idx++;
			this->counts[this->size-1]               = 1;
			return true;
		}else{
		}
	}else{
	}
	return false;
}



/** Used to compact a new nucleotide from a kmer on the right of the superkmer.
	* @param kmer The binary representation of the kmer to compact on the right. Must be on the same strand than the superkmer.
	* @return True if the kmer is inserted false otherwise.
	*/
	//NOT USED
bool SKC::compact_left(const uint64_t kmer_val) {
	uint64_t begin_sk = this->sk >> (this->size * 2);
	begin_sk &= (((uint64_t)1 << (2 * (k - 1))) - 1);
	uint64_t end_kmer = kmer_val & (((uint64_t)1 << (2 * (k - 1))) - 1);
	if (begin_sk == end_kmer) {
		this->sk += ((__uint128_t)(kmer_val >> (2 * (k - 1)))) << (2 * (k + this->size - 1));
		//           Select 2 left bits         Shift to the beginning of the sk
		this->size += 1;
		this->counts[this->size - 1] = 1;
		return true;
	}
	return false;
}



/** Look for the presence of the kmer inside of the superkmer.
	* The function supposes that the minimizer in the kmer is in the same strand than in sk.
	* @param kmer_val The binary value for the kmer.
	* @param minimizer_idx The index of the minimizer in the kmer_val
	* @return True if the kmer is present inside of the sk.
	*/
bool SKC::is_present(uint64_t kmer_val, uint64_t kmer_minimizer_idx) {
	int64_t start_idx  = (int64_t)this->minimizer_idx - (int64_t)kmer_minimizer_idx;
	if(start_idx<0 or (start_idx>=this->size)){return false;}
	uint64_t aligned_sk = (this->sk >> (2 * start_idx)) & k_mask;
	return aligned_sk == kmer_val;
}



bool SKC::is_present(kmer_full kmf) {
	if(minimizer_idx>=kmf.minimizer_idx){
		if((int)kmf.minimizer_idx>=(int)minimizer_idx+1-size){
			if(kmf.prefix==get_prefix()%(1<<(2*(31-minimizer_size-kmf.minimizer_idx)))) {
				if(kmf.suffix==get_suffix()>>(2*(minimizer_idx -kmf.minimizer_idx))) {
					return true;
				}
			}
		}
	}
	return false;
}



bool SKC::is_present_brutforce(kmer_full kmer, uint8_t & mini_k_idx) {
	for (mini_k_idx=0 ; mini_k_idx<=k-minimizer_size ; mini_k_idx++) {
		if (this->is_present(kmer.kmer_s, mini_k_idx))
			return true;
		// if (this->is_present(kmer.kmer_rc, k-minimizer_size-mini_k_idx)){
		// 	mini_k_idx = k-minimizer_size-mini_k_idx;
		// 	return true;
		// }
	}
	return false;
}







bool SKC::add_kmer(const kmer_full& kmer) {
	bool present = this->is_present(kmer);
	if(present){
		this->counts[kmer.minimizer_idx - (this->minimizer_idx-size+1)] ++;
		return true;
	}else{
		// The kmer is not found in the skc, try to compact
		if (this->size<32-minimizer_size) {
			return this->compact_right(kmer);
		}else{
			return false;
		}
	}
	return false;
}



// --- Pretty printing functions ---
void _out_kmer(ostream& out, __uint128_t kmer, uint64_t size) {
	kmer &= (((__uint128_t)1) << (2 * size)) - 1;
	__uint128_t anc((__uint128_t)1 << (2 * (size - 1)));
	for (uint64_t i(0); i < size and anc != 0; ++i) {
		uint64_t nuc = kmer / anc;
		kmer         = kmer % anc;
		if (nuc == 2) {
			out << "T";
		}
		if (nuc == 3) {
			out << "G";
		}
		if (nuc == 1) {
			out << "C";
		}
		if (nuc == 0) {
			out << "A";
		}
		if (nuc >= 4) {
			out << nuc << endl;
			out << "WTF" << endl;
		}
		anc >>= 2;
	}
}



ostream& operator<<(ostream& out, const SKC& skc) {
	// Print the superkmer
	_out_kmer(out, skc.sk, k + skc.size - 1);
	out << endl;
	// Print the minimizer
	uint64_t left_spaces = k + skc.size - minimizer_size - skc.minimizer_idx - 1;
	for (uint64_t i = 0; i < left_spaces; i++)
		out << ' ';
	_out_kmer(out, skc.sk >> (2 * skc.minimizer_idx), minimizer_size);
	out << " (mini idx " << (uint64_t)skc.minimizer_idx << ")" << endl;
	// Print the counters (same direction than superkmer)
	for (uint64_t i = skc.size; i > 0; i--)
		out << static_cast<uint64_t>(skc.counts[i - 1]) << "\t";
	return out;
}
