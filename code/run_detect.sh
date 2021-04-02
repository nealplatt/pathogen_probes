#on my local desktop

cd ~/Desktop/patho_detect

#activate environment
conda activate main


Pat-C1_S1 Pat-C2_S2 Pat-C3_S3 Pat-E1_S4 Pat-E2_S5 Pat-E3_S6 Undetermined_S0



#quality control and trimming
for SAMPLE in Pat-C1_S1 Pat-C2_S2 Pat-C3_S3 Pat-E1_S4 Pat-E2_S5 Pat-E3_S6 Undetermined_S0; do

    echo $SAMPLE
    R1_READS=data/seq_data/"$SAMPLE"_L001_R1_001.fastq.gz 
    R2_READS=data/seq_data/"$SAMPLE"_L001_R2_001.fastq.gz 

    #------------------
    #check quality with fastqc
    if [ ! -d results/fastqc ]; then
        mkdir results/fastqc
    fi

    fastqc \
        -o results/fastqc \
        -t 12 \
        $R1_READS \
        $R2_READS
    #------------------

    #------------------
    #clean up reads with trimmomatic
    if [ ! -d results/filter_reads ]; then
        mkdir results/filter_reads
    fi

    trimmomatic \
        PE \
        -threads 12 \
        -phred33 \
        $R1_READS \
        $R2_READS \
        results/filter_reads/"$SAMPLE"_filtered_R1_PE.fastq.gz \
        results/filter_reads/"$SAMPLE"_filtered_R1_SE.fastq.gz \
        results/filter_reads/"$SAMPLE"_filtered_R2_PE.fastq.gz \
        results/filter_reads/"$SAMPLE"_filtered_R2_SE.fastq.gz \
        LEADING:15 \
        TRAILING:15 \
        SLIDINGWINDOW:4:20 \
        MINLEN:36 \
        >results/filter_reads/"$SAMPLE".trimmomatic.log 2>&1 
 
        #combine SE reads and clean up
        zcat results/filter_reads/"$SAMPLE"_filtered_R1_SE.fastq.gz \
            results/filter_reads/"$SAMPLE"_filtered_R2_SE.fastq.gz \
            | gzip >results/filter_reads/"$SAMPLE"_filtered_RX_SE.fastq.gz
    
       rm results/filter_reads/"$SAMPLE"_filtered_R1_SE.fastq.gz
       rm results/filter_reads/"$SAMPLE"_filtered_R2_SE.fastq.gz
    #------------------

done 


#------------------
##classify and quantify reads
if [ ! -d results/kraken2  ]; then
    mkdir results/kraken2 
fi

##build kraken2 db
kraken2-build --standard --threads 2 --db results/kraken2/kraken2_standard_db

for SAMPLE in Pat-C1_S1 Pat-C2_S2 Pat-C3_S3 Pat-E1_S4 Pat-E2_S5 Pat-E3_S6 Undetermined_S0; do
    echo $SAMPLE

    #quanify/classify
    kraken2 \
        --use-names \
        --threads 12 \
        --db results/kraken2/kraken2_standard_db \
        --report results/kraken2/"$SAMPLE"_standard_kraken_report.out \
        --classified-out results/kraken2/"$SAMPLE"_standard_classifed#.fq \
        --unclassified-out results/kraken2/"$SAMPLE"_standard_unclassifed#.fq \
        --gzip-compressed \
        --paired \
        results/filter_reads/"$SAMPLE"_filtered_R1_PE.fastq.gz \
        results/filter_reads/"$SAMPLE"_filtered_R2_PE.fastq.gz \
        >results/kraken2/"$SAMPLE"_standard_kraken_results.tbl
    #------------------

done


