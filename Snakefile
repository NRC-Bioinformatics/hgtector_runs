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
        #wdir + "/analyses/order_level_transfers.csv",
        #wdir + "/merged_genomes/proteins.fasta",
        #wdir + "/merged_genomes/proteins.fasta.dmnd"
        #expand(wdir + "/searches/{stem}_vs_matches.csv", stem = stems),
        wdir + "/analyses/self_level_transfers.csv"

rule target2:
    input:
        expand(wdir + "/searches/{stem}_vs_matches.csv", stem = stems)

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

rule merge_proteins:
    input: 
        expand(config["bins"] + "/{stem}.faa",
        stem = stems)
    output:
        wdir + "/merged_genomes/proteins.fasta"
    shell:
        """
        cat {input} > {output}
        """

rule make_dimond_db4prot:
    input:
        wdir + "/merged_genomes/proteins.fasta"
    output:
        wdir + "/merged_genomes/proteins.fasta.dmnd"
    shell:
        """
        makedb --in {input} --db {input}

        """

rule searchVSallprots:
    input: config["bins"] + "/{stem}.faa"
    output:
        of = wdir + "/searches/{stem}_vsallbins.txt"
    params:
        db = wdir + "/merged_genomes/proteins.fasta",
        tmp = wdir + "/searches/{stem}_tmp_vsallprots"
    threads: 96
    conda: "hgtector"
    shell:
        """
        mkdir -p {params.tmp}
         diamond blastp -d {params.db} --threads {threads}  --query-cover 50 --max-target-seqs 100000000 \
         --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp --query {input} --out {output.of} --evalue 10 --tmpdir {params.tmp}
        """

rule searchVSnr:
    input: config["bins"] + "/{stem}.faa"
    output:
        of = wdir + "/searches/{stem}_vsnr.txt"
    params:
        db = config["hgdb"]+"/diamond/db.dmnd",
        tmp = wdir + "/searches/{stem}_tmp_vsnr"
    threads: 96
    conda: "hgtector"
    shell:
        """
        mkdir -p {params.tmp}
         diamond blastp -d {params.db} --threads {threads}  --query-cover 50 --max-target-seqs 100000000 \
         --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids --query {input} --out {output.of} --evalue 10 --tmpdir {params.tmp}
        """

rule wxtract_matches:
    input:
        vsbalt =  wdir + "/searches/{stem}_vsnr.txt",
        vsself =  wdir + "/searches/{stem}_vsallbins.txt",
    output:
        vsself =  wdir + "/searches/{stem}_vs_matches.csv"
    threads: 16
    notebook:
        "notebooks/analyse_sekf.r.ipynb"        

#1.0e-5 

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

rule analyse_donors_fine:
    input:
        info = expand(wdir + "/searches/{stem}_vs_matches.csv", stem = stems),
    output:
        plot = wdir + "/analyses/self_level_transfers.svg",
        csv = wdir + "/analyses/self_level_transfers.csv"
    notebook:
        "notebooks/self_level_analysis.r.ipynb"


# rule exctract_not_cloes:
#     input: directory(wdir + "/analyses/{stem}")
#     output:
#         wdir + "/analyses/{stem}_notclose.tsv"
#     notebook:
#         "notebooks/extract_not_close.r.ipynb"1122c