
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
    run:
        import os
        failed = False
        for file in input:
            if os.stat(file).st_size != 0:
                print(f"Kmer difference with the reference is not empty for {file}", file=sys.stderr)
                failed = True
        if failed:
            exit(1)

rule compare_kfcs_kmc:
    input:
        kfc="bench/kfc_{color}_sortedkmers.txt",
        kmc="bench/kmc_sortedkmers.txt"
    output:
        kfc_o="bench/diff_{color}.txt"
    run:
        shell(f"python3 -O scripts/eval_kmers.py {{input.kfc}} {{input.kmc}} > {{output.kfc_o}}")

rule sort_connonize_counts:
    input:
        input="bench/{name}_kmers.txt"
    output:
        output="bench/{name}_sortedkmers.txt"
    shell:
        "python3 -O scripts/canonize_kmer_counts.py {input} | LC_ALL=C sort -o {output}"

rule kfc_red_exec:
    input:
        data=INPUT,
        bin="bin/kfc_red"
    output:
        "bench/kfc_red_kmers.txt"
    shell:
        "{input.bin} {input.data} > {output};"

rule kfc_blue_exec:
    input:
        data=INPUT,
        bin="bin/kfc_blue"
    output:
        "bench/kfc_blue_kmers.txt"
    shell:
        "{input.bin} {input.data} {KMC_BLUE_MEM} > {output};"

rule build_kfc:
    output:
        "bin/{exec}"
    shell:
        "make DEBUG=0 ASSERTS=1 {wildcards.exec}"

rule kmc_exec:
    input:
        data=INPUT,
        bin=f"{KMC_DIR}/bin/kmc",
        update=f"{KMC_DIR}/update.lock"
    output:
        "bench/kmc_kmers.txt"
    shell:
        f"{{input.bin}} -fm -k{K} -ci0 -cs2048 -cx4294967295 -m6 {{input.data}} tmp_bin . ;"
        "{input.bin}_dump -ci0 -cx4294967295 tmp_bin {output} ;"

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
        "make;"
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
