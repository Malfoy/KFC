#include <thread>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include "index_min.h"


using namespace std;
using namespace chrono;


void index_full::insert(uint64_t kmer){
	uint64_t pos(Hash.lookup(kmer));
	if(Values[pos].kmer==kmer){
		++Values[pos].count;
	}
}

void index_full::dump_counting(){
	for(uint i(0);i<Values.size();++i){
		cout<<Values[i].kmer<<" "<<Values[i].count<<"\n";
	}
}
