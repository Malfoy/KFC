
import argparse
import re
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description='Verify the coherence between a kmer count and a reference count.')
    parser.add_argument('kmer_count', help='The kmer count to verify (space or tab delimiter)')
    parser.add_argument('reference_count', help='The reference kmer count.')

    args = parser.parse_args()
    return args


def compare_files(verif, ref):
    # Init counters
    identical = 0
    diff_count = 0
    absent = 0
    invalide = 0

    v_line = verif.readline()
    r_line = ref.readline()

    # Compare files
    while v_line != '' and r_line != '':
        v_kmer, v_count = re.split("\t| ", v_line.strip())
        r_kmer, r_count = re.split("\t| ", r_line.strip())

        if v_kmer < r_kmer:
            invalide += 1
            print(f"+{v_kmer}")
            v_line = verif.readline()
        elif v_kmer > r_kmer:
            absent += 1
            print(f"-{r_kmer}")
            r_line = ref.readline()
        else:
            v_line = verif.readline()
            r_line = ref.readline()
            if v_count == r_count:
                identical += 1
            else:
                invalide += 1
                print(f"~{r_kmer} {v_count}->{r_count}")

    # Go to the end of the files
    while v_line != '':
        v_kmer, v_count = re.split("\t| ", v_line.strip())
        print(f"+{v_kmer}")
        invalide += 1
        v_line = verif.readline()
    while r_line != '':
        r_kmer, r_count = re.split("\t| ", r_line.strip())
        print(f"-{r_kmer}")
        absent += 1
        r_line = ref.readline()

    return identical, diff_count, absent, invalide


if __name__ == "__main__":
    args = parse_arguments()
    with open(args.kmer_count) as verif, open(args.reference_count) as ref:
        counts = compare_files(verif, ref)
        print(f"{counts[0]} valid kmers", file=sys.stderr)
        print(f"{counts[1]} kmers with wrong count", file=sys.stderr)
        print(f"{counts[2]} missed kmers", file=sys.stderr)
        print(f"{counts[3]} kmers that should not be present", file=sys.stderr)
