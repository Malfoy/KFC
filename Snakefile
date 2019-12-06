
KMC_REPO="https://github.com/refresh-bio/KMC.git"
KMC_DIR="bench/softwares/kmc"
INPUT="data/ecoli_150bp_50kreads.fasta"


rule kfc_red_build_exec:
    input:
        data=INPUT,
    output:
        "bench/kfc_red_kmers.txt"
    shell:
        "make kfc_red ;"
        "kfc_red {input.data} > {output}"


rule kmc_exec:
    input:
        data=INPUT,
        bin=f"{KMC_DIR}/bin/kmc",
        update=f"{KMC_DIR}/update.lock"
    output:
        "bench/kmc_kmers.txt"
    shell:
        "{input.bin} -fm {input.data} tmp_bin . ;"
        "{input.bin}_dump tmp_bin {output} ;"
        "rm tmp_bin* {input.update}"


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
