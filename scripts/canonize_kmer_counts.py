import argparse
import sys
stdout = sys.stdout.buffer

def parse_arguments():
    parser = argparse.ArgumentParser(description='Parse a kmer list and transform all the kmers to keep the smallest between itself and the reverse complement')
    parser.add_argument('kmer_file', help='A file containing a text kmer list with counts')

    args = parser.parse_args()
    return args

complement_lkt = bytes.maketrans(b'ATGCatgc', b'TACGtacg')

def infer_separator(fin):
    with open(args.kmer_file, 'rb') as fp:
        for line in fp:
            line = line.strip()
            if b' ' in line:
                return b' '
            elif b'\t' in line:
                return b'\t'
            else:
                raise ValueError(f"Unknown line speparator {line:r}")

if __name__ == "__main__":
    args = parse_arguments()
    separator = infer_separator(args.kmer_file)
    with open(args.kmer_file, 'rb') as fp:
        for line in fp:
            kmer, _, count_and_newline = line.partition(separator)
            rev_comp = kmer.translate(complement_lkt)[::-1]
            stdout.write((kmer if kmer <= rev_comp else rev_comp) + b'\t' + count_and_newline)
