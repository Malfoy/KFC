#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>
#include <iostream>

using namespace std;


struct SKC{
  __uint128_t sk;
  vector<uint8_t> counts;
};

struct kmer_full{
  uint64_t kmer_s;
  uint64_t kmer_rc;
};


int64_t kmer_in_super_kmer(const SKC& super_kmer,const kmer_full& kmer){
  __uint128_t superkmer=super_kmer.sk;
  uint64_t skc_size=super_kmer.counts.size();
  for(uint64_t i=0;i<skc_size;++i){
    if(kmer.kmer_s==(uint64_t)(superkmer&(((uint64_t)1<<62)-1)) or kmer.kmer_rc==(uint64_t)(superkmer&(((uint64_t)1<<62)-1))) {
      return i;
    }
    superkmer>>=2;
  }
  return -1;
}

void print_128(const __uint128_t val) {
  uint8_t * pointer = (uint8_t *)&val;
  for (int8_t i=15; i>=0 ; i--)
    cout << hex << (uint32_t)pointer[i] << "\t";
  cout << endl;

  cout << dec;
}

int64_t kmer_in_super_kmer2(const SKC& super_kmer,const kmer_full& kmer){
  __uint128_t sk = super_kmer.sk;
  __uint128_t sk_shifted = super_kmer.sk >> 48;

  uint64_t kmer_val[2];
  kmer_val[0] = kmer.kmer_s;
  kmer_val[1] = kmer.kmer_rc;

  for (uint64_t sk_shift=0 ; sk_shift<4 ; sk_shift++) {
    sk >>= 2;
    sk_shifted >>= 2;

    for (uint64_t idx=0 ; idx<2 ; idx++) {
      uint64_t km = kmer_val[idx];
      // Load and repeat the first and last byte of the kmer into 2 mm_128 registers
      const __m128i first_letter = _mm_set1_epi8((char)(km));
      const __m128i last_letter  = _mm_set1_epi8((char)(km >> 48));
      // Load the super kmer and the super kmer shifter into two 128 bits registers
      const __m128i sk_reg = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&sk));
      const __m128i sk_shift_reg  = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&sk_shifted));

      // Test byte equalities
      const __m128i eq_first = _mm_cmpeq_epi8(first_letter, sk_reg);
      const __m128i eq_last  = _mm_cmpeq_epi8(last_letter, sk_shift_reg);
      // Mask results
      uint32_t mask = _mm_movemask_epi8(_mm_and_si128(eq_first, eq_last));

      while (mask != 0) {
        // Get number of 0 on the less significant side
        // Same as position of the first 1 from right
        int pos = __builtin_ctz(mask);

        // Multiply position by 8 for byte shifting
        int match_shift = pos << 3;
        // Check the full sub-superker against the kmer
        if (km == ((sk >> match_shift) & 0x3FFFFFFFFFFFFFFF)) {
          return pos*4+sk_shift;
        }

        // Remove used bit
        mask = mask & (mask-1);
      }
    }
  }

  return -1;
}




int main (int argc, char ** argv){
  struct kmer_full kmer = {rand(), rand()};
  struct SKC skc;
  for (uint i=0 ; i<31 ; i++)
    skc.counts.push_back(0);
  skc.sk = 0;
  // print_128(skc.sk);
  skc.sk = rand();
  // print_128(skc.sk);
  skc.sk <<= 62;
  // print_128(skc.sk);
  skc.sk |= kmer.kmer_s;
  // print_128(skc.sk);
  skc.sk <<= 12;
  // print_128(skc.sk);
  // print_128((__uint128_t)kmer.kmer_s);

  auto pos = kmer_in_super_kmer2(skc, kmer);
  cout << pos << endl;

  return 0;
}