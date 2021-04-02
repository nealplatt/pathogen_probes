#on my local desktop

cd ~/Desktop/patho_detect

#activate environment
conda activate patho_detect

#quality control and trimming
for SAMPLE in 0x-0c-01p-control_S8 \
    0x-0c-10p-control_S7 \
    15x-60c-01p-exp2-b2_S5 \
    15x-60c-10p-exp2-a1_S1 \
    15x-63c-01p-exp2-b1_S4 \
    15x-63c-10p-exp2-a2_S2 \
    15x-65c-01p-exp2-b3_S6 \
    15x-65c-10p-exp2-a3_S3; do

    echo $SAMPLE
    R1_READS=data/seq_data/20200904_113012/Fastq/"$SAMPLE"_L001_R1_001.fastq.gz 
    R2_READS=data/seq_data/20200904_113012/Fastq/"$SAMPLE"_L001_R2_001.fastq.gz 

    #------------------
    #check quality with fastqc
    if [ ! -d results/fastqc/20200904 ]; then
        mkdir results/fastqc/20200904
    fi

    fastqc \
        -o results/fastqc/20200904 \
        -t 2 \
        $R1_READS \
        $R2_READS
    #------------------

    #------------------
    #clean up reads with trimmomatic
    if [ ! -d results/filter_reads/20200904 ]; then
        mkdir results/filter_reads/20200904
    fi

    trimmomatic \
        PE \
        -threads 2 \
        -phred33 \
        $R1_READS \
        $R2_READS \
        results/filter_reads/20200904/"$SAMPLE"_filtered_R1_PE.fastq.gz \
        results/filter_reads/20200904/"$SAMPLE"_filtered_R1_SE.fastq.gz \
        results/filter_reads/20200904/"$SAMPLE"_filtered_R2_PE.fastq.gz \
        results/filter_reads/20200904/"$SAMPLE"_filtered_R2_SE.fastq.gz \
        LEADING:15 \
        TRAILING:15 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36 \
        >results/filter_reads/20200904/"$SAMPLE".trimmomatic.log 2>&1 
 
        #combine SE reads and clean up
        zcat results/filter_reads/20200904/"$SAMPLE"_filtered_R1_SE.fastq.gz \
            results/filter_reads/20200904/"$SAMPLE"_filtered_R2_SE.fastq.gz \
            | gzip >results/filter_reads/20200904/"$SAMPLE"_filtered_RX_SE.fastq.gz
    
       rm results/filter_reads/20200904/"$SAMPLE"_filtered_R1_SE.fastq.gz
       rm results/filter_reads/20200904/"$SAMPLE"_filtered_R2_SE.fastq.gz
    #------------------

done 


#------------------
##classify and quantify reads
if [ ! -d results/kraken2/20200904/  ]; then
    mkdir results/kraken2/20200904/ 
fi

##build kraken2 db
kraken2-build --standard --threads 2 --db results/kraken2/kraken2_standard_db

for SAMPLE in 0x-0c-01p-control_S8 \
    0x-0c-10p-control_S7 \
    15x-60c-01p-exp2-b2_S5 \
    15x-60c-10p-exp2-a1_S1 \
    15x-63c-01p-exp2-b1_S4 \
    15x-63c-10p-exp2-a2_S2 \
    15x-65c-01p-exp2-b3_S6 \
    15x-65c-10p-exp2-a3_S3; do

    echo $SAMPLE

    #quanify/classify
    kraken2 \
        --use-names \
        --threads 12 \
        --db results/kraken2/kraken2_standard_db \
        --report results/kraken2/20200904/"$SAMPLE"_standard_kraken_report.out \
        --classified-out results/kraken2/20200904/"$SAMPLE"_standard_classifed#.fq \
        --unclassified-out results/kraken2/20200904/"$SAMPLE"_standard_unclassifed#.fq \
        --gzip-compressed \
        --paired \
        results/filter_reads/20200904/"$SAMPLE"_filtered_R1_PE.fastq.gz \
        results/filter_reads/20200904/"$SAMPLE"_filtered_R2_PE.fastq.gz \
        >results/kraken2/20200904/"$SAMPLE"_standard_kraken_results.tbl
    #------------------

done



#####################################################################################
conda activate patho_detect

mkdir results/bbmap/20200904
cd results/bbmap/20200904


#download reference genomes
#GCF_005156105.1_ASM515610v1_genomic        mycBov_GCF_005156105.fa
#GCF_000195955.2_ASM19595v2_genomic         mycTub_GCF_000195955.fa
#GCF_000002765.5_GCA_000002765_genomic      plaFal_GCF_000002765.fa
#GCF_000002415.2_ASM241v2_genomic           plaViv_GCF_000002415.fa
#GCA_003958945.1_ASM395894v1_genomic        schBov_GCA_003958945.fa
#GCF_000237925.1_ASM23792v2_genomic         schMan_GCF_000237925.fa
#GCF_000001635.26_GRCm38.p6_genomic.fna     musMus_GCF_000001635.fa
#GCF_000001405.39_GRCh38.p13_genomic.fna    homSap_GCF_000001405.fa

#mv GCF_000001635.26_GRCm38.p6_genomic.fna musMus_GCF_000001635.fa
#mv GCF_000001405.39_GRCh38.p13_genomic.fna homSap_GCF_000001405.fa

#clean up the genome
for ASSEMBLY in musMus_GCF_000001635 \
    homSap_GCF_000001405 \
    mycBov_GCF_005156105 \
    mycTub_GCF_000195955 \
    plaFal_GCF_000002765 \
    plaViv_GCF_000002415 \
    schBov_GCA_003958945 \
    schMan_GCF_000237925 ; do

    cut -f1 -d" " ../$ASSEMBLY.fa >$ASSEMBLY.fas &
done

wait

#map reads to genomes
for SAMPLE in 15x-60c-10p-exp2-a1_S1 \
    0x-0c-01p-control_S8 \
    0x-0c-10p-control_S7 \
    15x-60c-01p-exp2-b2_S5 \
    15x-63c-01p-exp2-b1_S4 \
    15x-63c-10p-exp2-a2_S2 \
    15x-65c-01p-exp2-b3_S6 \
    15x-65c-10p-exp2-a3_S3 ; do

    echo "####################################################################"
    echo "......................$SAMPLE........................................"
    echo "####################################################################"
     bbsplit.sh \
        in1=../../filter_reads/20200904/"$SAMPLE"_filtered_R1_PE.fastq.gz \
        in2=../../filter_reads/20200904/"$SAMPLE"_filtered_R2_PE.fastq.gz \
        ref=musMus_GCF_000001635.fas,homSap_GCF_000001405.fas,mycBov_GCF_005156105.fas,mycTub_GCF_000195955.fas,plaFal_GCF_000002765.fas,plaViv_GCF_000002415.fas,schBov_GCA_003958945.fas,schMan_GCF_000237925.fas  \
        basename="$SAMPLE/$SAMPLE"_%.sam \
        refstats="$SAMPLE"_refstats.file \
        ambiguous=best \
        ambiguous2=all \
        threads=24
done

#map probes
awk '/^>/{print ">" ++i; next}{print}' < decon_probes.fas >renamed_decon_probes.fas

for ASSEMBLY in musMus_GCF_000001635 \
    homSap_GCF_000001405 \
    mycBov_GCF_005156105 \
    mycTub_GCF_000195955 \
    plaFal_GCF_000002765 \
    plaViv_GCF_000002415 \
    schBov_GCA_003958945 \
    schMan_GCF_000237925 ; do

    echo "####################################################################"
    echo "......................$ASSEMBLY........................................"
    echo "####################################################################"

    #bbmap.sh in=renamed_decon_probes.fas out=$ASSEMBLY/$ASSEMBLY.sam ref=$ASSEMBLY.fas threads=36

    samtools sort $ASSEMBLY/$ASSEMBLY.sam >$ASSEMBLY/$ASSEMBLY.bam
    samtools index $ASSEMBLY/$ASSEMBLY.bam 
    
    bedtools bamtobed -i $ASSEMBLY/$ASSEMBLY.bam | bedtools merge >$ASSEMBLY/$ASSEMBLY"_loci.bed"
done

#get non-target regions
for ASSEMBLY in musMus_GCF_000001635 \
    homSap_GCF_000001405 \
    mycBov_GCF_005156105 \
    mycTub_GCF_000195955 \
    plaFal_GCF_000002765 \
    plaViv_GCF_000002415 \
    schBov_GCA_003958945 \
    schMan_GCF_000237925 ; do

    #samtools faidx $ASSEMBLY.fas

    bedtools slop -g $ASSEMBLY.fas.fai -i $ASSEMBLY/$ASSEMBLY"_loci.bed" -b 1000  | bedtools complement -g $ASSEMBLY.fas.fai -i - >$ASSEMBLY/$ASSEMBLY"_nontarget".bed
done


#calculate coverage of target vs. non-target regions
echo "taxa,protocol,target_cov,nontarget_cov" >covs.csv
for ASSEMBLY in mycBov_GCF_005156105 \
    mycTub_GCF_000195955 \
    plaFal_GCF_000002765 \
    plaViv_GCF_000002415 \
    schBov_GCA_003958945 \
    schMan_GCF_000237925 \
    musMus_GCF_000001635 \
    homSap_GCF_000001405 ; do

    echo $ASSEMBLY
    for SAMPLE in 0x-0c-01p-control_S8 \
        0x-0c-10p-control_S7 \
        15x-60c-01p-exp2-b2_S5 \
        15x-60c-10p-exp2-a1_S1 \
        15x-63c-01p-exp2-b1_S4 \
        15x-63c-10p-exp2-a2_S2 \
        15x-65c-01p-exp2-b3_S6 \
        15x-65c-10p-exp2-a3_S3 ; do
    
        echo $SAMPLE
        #samtools sort $SAMPLE/$SAMPLE"_"$ASSEMBLY".sam" >$SAMPLE/$SAMPLE"_"$ASSEMBLY".bam"
        #samtools index $SAMPLE/$SAMPLE"_"$ASSEMBLY".bam"

        #target
        mosdepth -t 4 -b $ASSEMBLY/$ASSEMBLY"_loci.bed"      $SAMPLE/$ASSEMBLY"_target"    $SAMPLE/$SAMPLE"_"$ASSEMBLY".bam"
        TARGET_COV=$(zcat $SAMPLE/$ASSEMBLY"_target.regions.bed.gz" | awk '{SUM += $4} END {print SUM/NR}')

        #non-target        
        mosdepth -t 4 -b $ASSEMBLY/$ASSEMBLY"_nontarget.bed" $SAMPLE/$ASSEMBLY"_nontarget" $SAMPLE/$SAMPLE"_"$ASSEMBLY".bam"
        NONTARGET_COV=$(zcat $SAMPLE/$ASSEMBLY"_nontarget.regions.bed.gz" | awk '{SUM += $4} END {print SUM/NR}')

        echo "$ASSEMBLY,$SAMPLE,$TARGET_COV,$NONTARGET_COV" >>covs.csv
            
    done
done


echo "taxa,protocol,num_reads,num_target_reads" >num_reads.csv
for SAMPLE in 0x-0c-01p-control_S8 \
    0x-0c-10p-control_S7 \
    15x-60c-01p-exp2-b2_S5 \
    15x-60c-10p-exp2-a1_S1 \
    15x-63c-01p-exp2-b1_S4 \
    15x-63c-10p-exp2-a2_S2 \
    15x-65c-01p-exp2-b3_S6 \
    15x-65c-10p-exp2-a3_S3 ; do
    
    echo $SAMPLE
       
    NUM_READS=$(zcat ../../filter_reads/20200904/"$SAMPLE"_filtered_R*_PE.fastq.gz | wc -l)
    NUM_READS=$(echo $(( NUM_READS / 4 )))

    for ASSEMBLY in mycBov_GCF_005156105 \
        mycTub_GCF_000195955 \
        plaFal_GCF_000002765 \
        plaViv_GCF_000002415 \
        schBov_GCA_003958945 \
        schMan_GCF_000237925 \
        musMus_GCF_000001635 \
        homSap_GCF_000001405 ; do

        echo ..$ASSEMBLY

        NUM_SPECIES_READS=$(samtools view $SAMPLE/$SAMPLE"_"$ASSEMBLY".bam" | wc -l)

        echo "$ASSEMBLY,$SAMPLE,$NUM_READS,$NUM_SPECIES_READS" >>num_reads.csv
            
    done
done

#clean up results
