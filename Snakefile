
KMC_REPO="https://github.com/refresh-bio/KMC.git"
KMC_DIR="bench/softwares/kmc"
INPUT="data/ecoli_150bp_50kreads.fasta"
K=31
KMC_BLUE_MEM=67

def raw_files(wildcards):
    if "color" in config and (config["color"] == "red" or config["color"] == "blue"):
        return f"bench/diff_{config['color']}.txt"
    return ["bench/diff_red.txt", "bench/diff_blue.txt"]

rule all:
    input: raw_files


rule compare_kfcs_kmc:
    input:
        kfc="bench/kfc_{color}_kmers.txt",
        kmc="bench/kmc_kmers.txt"
    output:
        kfc_o="bench/diff_{color}.txt"
    shell:
        f"python3 scripts/eval_kmers.py {{input.kfc}} {{input.kmc}} > {{output.kfc_o}};"


rule kfc_red_exec:
    input:
        data=INPUT,
    output:
        "bench/kfc_red_kmers.txt"
    shell:
        "./bin/kfc_red {input.data} > {output};"
        "python3 scripts/canonize_kmer_counts.py {output} > tmp_red_kmers.csv;"
        "sort tmp_red_kmers.csv -o {output};"
        "rm tmp_red_kmers.csv;"


rule kfc_blue_exec:
    input:
        data=INPUT,
    output:
        "bench/kfc_blue_kmers.txt"
    shell:
        "./bin/kfc_blue {input.data} {KMC_BLUE_MEM} > {output};"
        "python3 scripts/canonize_kmer_counts.py {output} > tmp_blue_kmers.csv;"
        "sort tmp_blue_kmers.csv -o {output};"
        "rm tmp_blue_kmers.csv;"


rule kmc_exec:
    input:
        data=INPUT,
        bin=f"{KMC_DIR}/bin/kmc",
        update=f"{KMC_DIR}/update.lock"
    output:
        "bench/kmc_kmers.txt"
    shell:
        f"{{input.bin}} -fm -k{K} -ci0 -cs2048 -cx4294967295 {{input.data}} tmp_bin . ;"
        "{input.bin}_dump -ci0 -cx4294967295 tmp_bin {output} ;"
        "python3 scripts/canonize_kmer_counts.py {output} > tmp_kmc_kmers.csv;"
        "sort tmp_kmc_kmers.csv -o {output};"
        "rm tmp_bin* tmp_kmc_kmers.csv {input.update}"


rule kmc_update_compile:
    input:
        f"{KMC_DIR}/clone.lock"
    output:
        version=f"{KMC_DIR}/kmc_version_checked.lock",
        up=f"{KMC_DIR}/update.lock",
        bin=f"{KMC_DIR}/bin/kmc"
    shell:
        f"cd {KMC_DIR}/;"
        "git pull origin master;"
        "make -j;"
        "./bin/kmc | head -n 1 | cut -d ' ' -f5 > ./kmc_version_checked.lock;"
        "git rev-parse HEAD >> ./kmc_version_checked.lock;"
        "cd -;"
        "touch {output.up}"

rule kmc_clone:
    output:
        f"{KMC_DIR}/clone.lock"
    shell:
        f"if [ ! -f {KMC_DIR}/README.md ]; then git clone {KMC_REPO} {KMC_DIR} ; fi;"
        "touch {output}"
