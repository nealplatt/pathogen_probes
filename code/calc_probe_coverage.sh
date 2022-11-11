conda activate calc_probe_coverage

cd /master/nplatt/pathogen_probes/results/seq_analyses/coverage

#get reads
cp ../../seq_analyses/01_clean-fastq/*/split-adapter-quality-trimmed/c*.gz .

#get genomes
zcat ../../../data/genomes/trematoda/Schistosoma_mansoni_GCA_000237925.fas.gz            | cut -f1 -d" " >sm.fas
zcat ../../../data/genomes/apicomplexa/Plasmodium_falciparum_GCA_000002765.fas.gz        | cut -f1 -d" " >pf.fas
zcat ../../../data/genomes/apicomplexa/Plasmodium_vivax_GCF_000002415.fas.gz             | cut -f1 -d" " >pv.fas
zcat ../../../data/genomes/mycobacterium/Mycobacterium_tuberculosis_GCF_000195955.fas.gz | cut -f1 -d" " >mt.fas

samtools faidx sm.fas
samtools faidx pf.fas
samtools faidx pv.fas
samtools faidx mt.fas

#get probes
cp ../../../final_probes_v1.1.fas .
cat final_probes_v1.1.fas | paste - - | grep -i -e mycobacterium | sed 's/\t/\n/' >mycobacterium_probes.fas
cat final_probes_v1.1.fas | paste - - | grep -i -e trematod | sed 's/\t/\n/' >trematoda_probes.fas
cat final_probes_v1.1.fas | paste - - | grep -i -e apicomplexa | sed 's/\t/\n/' >plasmodium_probes.fas

#get probe locations###################################################################################

#for mt
conda activate bbmap
bbmap.sh \
    ref=mt.fas \
    in=mycobacterium_probes.fas \
    maxindel=10 \
    threads=48 \
    ambiguous=all \
    out=mycobacterium_probes.sam

conda activate samtools

samtools view -F4 -Sb mycobacterium_probes.sam | samtools sort - >mycobacterium_probes.bam
samtools index mycobacterium_probes.bam

conda activate bedtools
bedtools bamtobed -i mycobacterium_probes.bam >mycobacterium_probes.bed
bedtools merge -d 100 -i mycobacterium_probes.bed >mycobacterium_probes.merged.bed
#49 loci


#get probe locations for sm
bbmap.sh \
    ref=sm.fas \
    in=trematoda_probes.fas \
    maxindel=10 \
    threads=48 \
    ambiguous=all \
    out=trematoda_probes.sam

samtools view -F4 -Sb trematoda_probes.sam | samtools sort - >trematoda_probes.bam
samtools index trematoda_probes.bam

bedtools bamtobed -i trematoda_probes.bam >trematoda_probes.bed
bedtools merge -d 100 -i trematoda_probes.bed >trematoda_probes.merged.bed
#49 loci



#get probe locations for pf
bbmap.sh \
    ref=pf.fas \
    in=plasmodium_probes.fas \
    maxindel=10 \
    threads=48 \
    ambiguous=all \
    out=plasmodium_probes.sam

samtools view -F4 -Sb plasmodium_probes.sam | samtools sort - >plasmodium_probes.bam
samtools index plasmodium_probes.bam

bedtools bamtobed -i plasmodium_probes.bam >plasmodium_probes.bed
bedtools merge -d 100 -i plasmodium_probes.bed >plasmodium_probes.merged.bed
#49 loci

###################################################################################


#map reads to reference genomes ###################################################################################

for REF in mt pf sm; do
    for SAMPLE in c1e_3p_control_enrich c1p_control_enrich c1p_control_no_enrich; do
        echo ${REF} vs ${SAMPLE}
        bbmap.sh \
            ref=${REF}.fas \
            in=${SAMPLE}-READ1.fastq.gz \
            in2=${SAMPLE}-READ2.fastq.gz \
            maxindel=10 \
            minid=0.85 \
            threads=48 \
            ambiguous=all \
            out=${SAMPLE}_VS_${REF}.sam

        samtools view -F4 -Sb ${SAMPLE}_VS_${REF}.sam | samtools sort - >${SAMPLE}_VS_${REF}.bam
        samtools index ${SAMPLE}_VS_${REF}.bam
    done
done

# get coverages at target loci ###################################################################################

#cp the probes just to standardize things
cp mycobacterium_probes.merged.bed mt_probes.merged.bed 
cp plasmodium_probes.merged.bed  pf_probes.merged.bed
cp trematoda_probes.merged.bed sm_probes.merged.bed


for REF in mt pf sm; do
    for SAMPLE in c1e_3p_control_enrich c1p_control_enrich c1p_control_no_enrich; do
        mosdepth -b ${REF}_probes.merged.bed -t 4 -m ${SAMPLE}_VS_${REF} ${SAMPLE}_VS_${REF}.bam
    done

    gunzip *${REF}.regions.bed.gz

    paste \
        c1p_control_no_enrich_VS_${REF}.regions.bed \
        c1e_3p_control_enrich_VS_${REF}.regions.bed \
        c1p_control_enrich_VS_${REF}.regions.bed \
        | cut -f1,2,3,4,8,12 \
            >tmp.tsv

    echo -e "chrom\tstart\tstop\tc1p_control_no_enrich\tc1e_3p_control_enrich\tc1p_control_enrich\n$(cat tmp.tsv)" >${REF}_sample_loci_covs.tsv
    rm tmp.tsv
done


################################################################

#now get random coverages outside of probed regions (+- 500 bp)
bedtools slop -g mt.fas.fai -i mycobacterium_probes.merged.bed -b 500 >mt_excl_loci.bed
bedtools slop -g pf.fas.fai -i plasmodium_probes.merged.bed    -b 500 >pf_excl_loci.bed
bedtools slop -g sm.fas.fai -i trematoda_probes.merged.bed     -b 500 >sm_excl_loci.bed

#####
for REF in mt pf sm; do
    bedtools slop -g ${REF}.fas.fai -i ${REF}_probes.merged.bed -b 500 >${REF}_excl_loci.bed


    for SAMPLE in c1e_3p_control_enrich c1p_control_enrich c1p_control_no_enrich; do
        echo -e "chrom\tstart\tstop\tcov">${REF}_VS_${SAMPLE}_random_loci_covs.bed
        for i in $(seq 1 100); do
            echo $SAMPLE: $i
            bedtools shuffle -seed $RANDOM -i ${REF}_probes.merged.bed -g ${REF}.fas.fai -noOverlapping -excl ${REF}_excl_loci.bed >shuff.bed

            mosdepth -b shuff.bed -t 4 -m random ${SAMPLE}_VS_${REF}.bam
            zcat random.regions.bed.gz >>${REF}_VS_${SAMPLE}_random_loci_covs.bed
        done
    done

    paste \
        ${REF}_VS_c1p_control_no_enrich_random_loci_covs.bed \
        ${REF}_VS_c1e_3p_control_enrich_random_loci_covs.bed \
        ${REF}_VS_c1p_control_enrich_random_loci_covs.bed \
        | cut -f4,8,12 \
        >tmp.tsv

    echo -e "c1p_control_no_enrich\tc1e_3p_control_enrich\tc1p_control_enrich\n$(cat tmp.tsv | sed 1d)" >${REF}_random_loci_covs.tsv
    rm tmp.tsv
done

