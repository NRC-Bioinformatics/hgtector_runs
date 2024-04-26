import pandas as pd
stems = glob_wildcards(config["bins"]+"/{stem}.faa").stem  #[0:1]
wdir = config["work_directory"]
print(config["bins"])
print(stems)

### get maps
dfm = pd.read_csv("resources/class.csv")
bs = ["nod_"+b for b in dfm.bin]
fids = [int(f.split()[1]) for f in dfm.family]
oids = [int(f.split()[1]) for f in dfm.order]


fmap = dict(zip( bs, fids))
omap = dict(zip( bs, oids))



rule target:
    input:
        # wdir + "/merged_genomes/chosen_species.txt"
        # #wdir + "/merged_genomes/genomes_length_info.txt",
        # #wdir + "/merged_genomes/genomes_info.txt"
        #wdir + "/analyses/most_abundant_groups_assignments.tsv"
        wdir + "/analyses/order_level_transfers.csv",
        # expand(
        #     #  wdir + "/renamed_genomes/{stem}.fasta",
        #     #wdir + "/kaiju/{stem}_info.txt",
        #     #wdir + "/kaiju/{stem}__length_info.txt",
        #     #wdir + "/analyses/{stem}_groups_assignments.txt",
        #     #wdir + "/analyses/{stem}_info.txt",
        #     wdir + "/analyses/{stem}_matches_to_community.tsv",
        #     #wdir + "/searches/{stem}_info.txt",
        #     stem = stems
        # )
rule analyse_kaiju_assign_taxonomy:
    input:
        linf = wdir + "/merged_genomes/genomes_length_info.txt",
        kaiju = wdir + "/merged_genomes/genomes_info.txt"
    output:
        wdir + "/merged_genomes/chosen_species.txt"
    notebook:
        "notebooks/choose_species.r.ipynb"



rule search:
    input: config["bins"] + "/{stem}.faa"
    output:
        od = directory(wdir + "/searches/{stem}"),
        of = wdir + "/searches/{stem}_info.txt"
    params:
        db = config["hgdb"],
        tmp = wdir + "/analyses/{stem}_tmp"
    threads: 96
    conda: "hgtector"
    shell:
        """
        mkdir -p {params.tmp}
        unbuffer hgtector search -i {input} -o {output.od} -m diamond -p {threads} --tmpdir {params.tmp} -d \
        {params.db}/diamond/db -t {params.db}/taxdump  --aln-method native  > {output.of}
        """

rule get_length:
    input: fa = wdir + "/merged_genomes/genomes.fasta"
    output:
        wdir + "/merged_genomes/genomes_length_info.txt"
    shell:
        """
        seqkit fx2tab -n -l {input} -o {output}
        """
rule merge:
    input: 
        expand("/mnt/workspace/analyses/omics/hgtector/nodbins/" + "{stem}.fasta",
        stem = stems)
    output:
        wdir + "/merged_genomes/genomes.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule search_kaiju:
    input: fa = wdir + "/merged_genomes/genomes.fasta"
    output:
        wdir + "/merged_genomes/genomes_info.txt"
    params:
        db = "/mnt/workspace/databases/kaiju/kaiju_db_refseq_nr_2023-06-17/",
    conda:
        "kaiju"
    threads: 96
    shell:
        """
        kaiju -t {params.db}/nodes.dmp -f {params.db}/kaiju_db_refseq_nr.fmi -i {input} -v -o {output} -z {threads} -E 100 
        """


def get_self_gr(wildcards):
    out = fmap[wildcards.stem]
    if wildcards.stem == "nod_bin.13":
        out = 1649453
    return(out)


rule analyse:
    input: directory(wdir + "/searches/{stem}")
    output:
        od = directory(wdir + "/analyses/{stem}"),
        of = wdir + "/analyses/{stem}_info.txt"
    params:
        db = config["hgdb"],
        tmp = wdir + "/analyses//{stem}_tmp",
        selftax = get_self_gr,
        closetax = lambda wildcards: omap[wildcards.stem]
    threads: 96
    conda: "hgtector"
    shell:
        """
        mkdir -p {params.tmp} 
        unbuffer hgtector analyze  --self-tax {params.selftax} --close-tax {params.closetax}  -i {input} -o {output.od} -t {params.db}/taxdump \
        --bandwidth grid   --identity 50 --coverage 50 > {output.of}
        """

rule extract_matches_to_comunity:
    input:
        od = directory(wdir + "/analyses/{stem}"),
        close = "resources/class.csv"
    output:
        wdir + "/analyses/{stem}_matches_to_community.tsv"
    notebook:
        "notebooks/close_to_community.r.ipynb"



rule analyse_close_group:
    input:
        search =  wdir + "/searches/{stem}",
        info = wdir + "/merged_genomes/genomes_length_info.txt",
    output:
         wdir + "/analyses/{stem}_groups_assignments.txt"
        #limit = wdir + "/analyses/{stem}_groups_limits.txt",
    threads: 40
    notebook:
        "notebooks/extract_close_group.r.ipynb"

rule merge_group_assignments:
    input:
            expand(
            wdir + "/analyses/{stem}_groups_assignments.txt",
            stem = stems)
    output:
        wdir + "/analyses/most_abundant_groups_assignments.tsv"
    notebook:
        "notebooks/merge_assignments.r.ipynb"

rule analyse_donnors:
    input:
        info = expand(
            wdir + "/analyses/{stem}_matches_to_community.tsv",
            stem = stems)
    output:
        plot = wdir + "/analyses/order_level_transfers.svg",
        csv = wdir + "/analyses/order_level_transfers.csv"
    notebook:
        "notebooks/order_level_analysis.r.ipynb"



# rule exctract_not_cloes:
#     input: directory(wdir + "/analyses/{stem}")
#     output:
#         wdir + "/analyses/{stem}_notclose.tsv"
#     notebook:
#         "notebooks/extract_not_close.r.ipynb"1122c