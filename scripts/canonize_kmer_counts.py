import argparse
import re
from Bio.Seq import Seq

def parse_arguments():
    parser = argparse.ArgumentParser(description='Parse a kmer list and transform all the kmers to keep the smallest between itself and the reverse complement')
    parser.add_argument('kmer_file', help='A file containing a text kmer list with counts')

    args = parser.parse_args()
    return args


def canonize(kmer):
    s = Seq(kmer)
    rev_comp = str(s.reverse_complement())

    return s if s <= rev_comp else rev_comp


if __name__ == "__main__":
    args = parse_arguments()
    with open(args.kmer_file) as fp:
        for line in fp:
            line = line.strip()
            kmer, count = re.split("\t| ", line)
            kmer = canonize(kmer)
            print(f"{kmer}\t{count}")
