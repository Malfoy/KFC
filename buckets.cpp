#include <algorithm>
#include "buckets.hpp"



using namespace std;



void  Bucket::add_kmers( vector<kmer_full>& kmers){
	if(not add_kmers_sorted(kmers)){
		add_kmers_buffer(kmers);
	}
	if(skml.size()-sorted_size>1000000){
		insert_buffer();
	}
	kmers.clear();
}



void Bucket::insert_buffer(){
	sort(skml.begin()+sorted_size , skml.end(),[ ]( const SKC& lhs, const SKC& rhs ){return lhs < rhs;});
	inplace_merge(skml.begin(), skml.begin()+sorted_size ,skml.end() ,[ ]( const SKC& lhs, const SKC& rhs ){return lhs < rhs;});
	sorted_size=skml.size();
}



void  Bucket::add_kmers_buffer( vector<kmer_full>& kmers){
	uint64_t inserted(0);
	uint64_t buffsize(skml.size());
	//HERE IF CHECK IF THE KMER ARE IN THE UNSORTED BUFFER
	//FOREACH SUPERKMER
	for (uint64_t i = sorted_size; i <buffsize; i++) {
		SKC& skc = skml[i];
		//FOREACH KMER
		for (uint64_t ik = 0; ik < kmers.size(); ++ik) {
			kmer_full& kmer = kmers[ik];
			if (kmer.minimizer_idx!=69) {
				if (skc.add_kmer(kmer)) {
					++inserted;
					kmers[ik].minimizer_idx=69;
				}
			}
		}
	}
	//HERE WE CREATE NEW SUPERKMERS
	if(inserted!=kmers.size()){
		//FOREACH KMER
		for (uint64_t ik = 0; ik < kmers.size(); ++ik) {
			kmer_full& kmer = kmers[ik];
			if (kmer.minimizer_idx!=69) {
				//WE TRY TO COMPACT IT TO THE LAST SUPERKMER
				if(not skml.empty()){
					if(skml[skml.size()-1].compact_right(kmer)){
						continue;
					}
				}
				bool isinserted=false;
				//FOREACH new SUPERKMER
				for (uint64_t i = buffsize; i < skml.size(); i++) {
					if (skml[i].add_kmer(kmer)) {
						isinserted=true;
						break;
					}else{
					}
				}
				if(not isinserted){
					skml.push_back(SKC(kmer.get_compacted(), kmer.get_minimizer_idx()));
				}
			}
		}
	}
}



bool compSKM(const SKC& s1, const SKC& s2){
	return s1<s2;
}



bool  Bucket::add_kmers_sorted( vector<kmer_full>& kmers){
	if(sorted_size==0){
		return false;
	}
	uint64_t inserted(0);
	//FOREACH KMER
	for (uint64_t ik = 0; ik < 1; ++ik) {
		kmer_full kmer = kmers[ik];
		if(kmer.minimizer_idx==69){continue;}
		SKC mockskm( kmer.kmer_s, kmer.minimizer_idx);
		int64_t low=lower_bound (skml.begin(), skml.begin()+sorted_size,mockskm,[ ]( const SKC& lhs, const SKC& rhs ){return lhs < rhs;}) - skml.begin();
		//FOREACH SUPERKMER
		while (low<(int)sorted_size and skml[low].suffix_is_prefix(kmer)) {
			//FOREACH KMER
			for (uint64_t iikk = ik; iikk < kmers.size(); ++iikk) {
				if(kmers[iikk].minimizer_idx==69){continue;}
				if (skml[low].add_kmer(kmers[iikk])) {
					kmers[iikk].minimizer_idx=69;
					inserted++;
					if(inserted==kmers.size()){
						return true;
					}
				}
			}
			low++;
		}
		// cout<<sorted_size<<" "<<(low-lowmem)<<" "<<sorted_size-(low-lowmem)<<" "<<(double)(low-lowmem)/sorted_size<<endl;
		// cout<<kmers.size()<<endl;
		// cout<<KMERCHECK-memKMER<<endl;

		// cin.get();
		//FOREACH KMER
		// for (uint64_t iikk = ik; iikk < kmers.size(); ++iikk) {
		// 	kmer_full& kmerl = kmers[iikk];
		// 	// if(kmerl.minimizer_idx==69){continue;}
		// 	int llow=low;
		// 	while (llow<(int)sorted_size and skml[llow].suffix_is_prefix(kmerl)) {
		// 		if (skml[llow].add_kmer(kmerl)) {
		// 			kmerl.minimizer_idx=69;
		// 			inserted++;
		// 			if(inserted==kmers.size()){
		// 				return true;
		// 			}
		// 			break;
		// 		}
		// 		llow++;
		// 	}
		// }
	}
	return false;
}



void  Bucket::print_kmers(string& result,const  string& mini){
	insert_buffer();
	for(uint64_t i(0);i<skml.size();++i){
		skml[i].print_count(result,mini);
	}
}



uint64_t Bucket::size()const{
	return skml.size();
}



uint64_t Bucket::number_kmer()const{
	uint64_t result(0);
	for(uint64_t i(0);i<skml.size();++i){
		result+=skml[i].size;
	}
	return result;
}
