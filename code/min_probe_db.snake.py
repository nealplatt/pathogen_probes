import os
import pandas as pd
from snakemake.utils import min_version

###################################################################################

#NOTES:

#snakemake \
#    --use-conda \
#    --printshellcmds \
#    --cluster 'qsub -V -cwd -S /bin/bash -pe smp {threads} -o {log}.log -j y -q all.q' \
#    --jobs 100 \
#    --latency-wait 200 \
#    --rerun-incomplete \
#    --keep-going \
#    --snake /master/nplatt/pathogen_probes/code/min_probe_db.snake.py

#using
# - snakemake v6.12.3
# - using conda v4.11.0


###################################################################################
##### set minimum snakemake version #####
min_version("6.12.3")

#set main project dir and work from there
proj_dir   = "/master/nplatt/pathogen_probes/results/seq_analyses/min_probe_db"
genome_dir = "/master/nplatt/pathogen_probes/results/seq_analyses/min_probe_db/refseq_212_fas/v1/files"
envs_dir   = "/master/nplatt/pathogen_probes/env"
logs_dir   = "{}/logs".format(proj_dir)

probe_fas   = "/master/nplatt/pathogen_probes/results/seq_analyses/min_probe_db/final_probes_v1.0.fas"
genome_list = "/master/nplatt/pathogen_probes/results/seq_analyses/min_probe_db/genomes.list"

#get sample info
genomes_df  = pd.read_csv(genome_list, header=None)
genomes     = list(genomes_df[0])

localrules:
    all,

rule all:
    input:
        expand("{dir}/{id}_genomic.fna.gz",     dir = genome_dir,
                                                id  = genomes),
        expand("{dir}/target_regions/{id}.fas", dir = proj_dir,
                                                id  = genomes),
        proj_dir + "/target_loci.fas"

rule process_genome:
    input:
        fas_gz    = genome_dir + "/{id}_genomic.fna.gz"
    output:
        fas       = temp(proj_dir + "/{id}_genomic.fna"),
        fai       = temp(proj_dir + "/{id}_genomic.fna.fai"),
        bed       = temp(proj_dir + "/{id}_genomic.windows.bed"),
        split_fas = temp(proj_dir + "/{id}_genomic.windows.fas"),
        split_fai = temp(proj_dir + "/{id}_genomic.windows.fas.fai")
    threads:
        1
    log:
        logs_dir + "/process_genome#{id}"
    conda:
        envs_dir + "/min_probe_db.yml"
    shell:
        """
        zcat {input.fas_gz} >{output.fas}

        samtools faidx {output.fas}

        bedtools makewindows  -g {output.fai} -w 500000000 >{output.bed}

        bedtools getfasta \
            -fi {output.fas} \
            -fo {output.split_fas} \
            -bed {output.bed}

        samtools faidx {output.split_fas}
        """

rule map_probes:
    input:
        split_fas = proj_dir + "/{id}_genomic.windows.fas",
        probe_fas = probe_fas
    output:
        sam       = temp(proj_dir + "/{id}.sam")
    params:
        index_dir = proj_dir + "/scratch/{id}/"
    threads:
        24
    log:
        logs_dir + "/map_probes#{id}"
    conda:
        envs_dir + "/min_probe_db.yml"
    shell:
        """
        bbmap.sh \
            ref={input.split_fas} \
            path={params.index_dir} \
            in={input.probe_fas} \
            slow=t \
            minid=0.85 \
            threads={threads} \
            ambiguous=all \
            overwrite=t \
            secondary=t \
            sssr=0.8 \
            maxsites=10 \
            sam=1.4 \
            out={output.sam} \
            -Xmx200g 

        rm -r {params.index_dir}
        """

rule extract_tagets:
    input:
        split_fas = proj_dir + "/{id}_genomic.windows.fas",
        sam       = proj_dir + "/{id}.sam",
        split_fai = proj_dir + "/{id}_genomic.windows.fas.fai"
    output:
        bed       = temp(proj_dir + "/{id}.loci.bed"),
        fas       = proj_dir + "/target_regions/{id}.fas"
    threads:
        1
    log:
        logs_dir + "/extract_tagets#{id}"
    conda:
        envs_dir + "/min_probe_db.yml"
    shell:
        """
        #filter unmapped
        samtools view \
            -h \
            -Sb \
            -F 4 \
            {input.sam} \
            | samtools sort \
            | bedtools bamtobed \
            | bedtools slop \
                -b 1000 \
                -g {input.split_fai} \
                -i - \
            | bedtools merge \
                -d 0 \
                >{output.bed}

        bedtools getfasta \
            -fi {input.split_fas} \
            -fo {output.fas} \
            -bed {output.bed}

        sed -i 's/:/ /' {output.fas}
        """

rule gen_target_file:
    input:
        fas       = expand("{dir}/target_regions/{id}.fas", dir = proj_dir, 
                                                            id  = genomes)
    output:
        fas       = proj_dir + "/target_loci.fas"
    threads:
        1
    log:
        logs_dir + "/gen_target_file"
    conda:
        envs_dir + "/min_probe_db.yml"
    shell:
        """
        cat {input.fas} >{output.fas}
        rm -r scratch
        """
