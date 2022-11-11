conda activate phyluce

TAXID=119060
TAXA="burkholderiaceae"
OUTGROUP_TAXA=""
OUTGROUP_TAXID=""
PROBE_GENUS="burkholderia"

TAXID=772
TAXA="bartonella"
OUTGROUP_TAXA=""
OUTGROUP_TAXID=""
PROBE_GENUS="bartonella"

TAXID=31245
TAXA="schistosomatidae"
OUTGROUP_TAXA=""
OUTGROUP_TAXID=""
PROBE_GENUS="TEST"

TAXID=1763
TAXA="mycobacterium"
OUTGROUP_TAXA=""
OUTGROUP_TAXID=""
PROBE_GENUS="mycobacterium"

# TAXID=5820
# TAXA="plasmodium"
# OUTGROUP_TAXA=""
# OUTGROUP_TAXID=""
# PROBE_GENUS="plasmodium"

N_THREADS=48
PROJ_DIR=/master/nplatt/pathogen_probes/results/seq_analyses
TAXA_DIR=${PROJ_DIR}/03_${TAXA}_phylogeny

SAMPLE_CSV=${PROJ_DIR}/sample_info.csv

mkdir -p $TAXA_DIR
cd $TAXA_DIR

mkdir ${TAXA_DIR}/00_confs ${TAXA_DIR}/00_logs

for SAMPLE in $(cut -f2 -d"," ${SAMPLE_CSV} | sed 1d); do
    #input files
    IN_FQ_R1=${PROJ_DIR}/01_clean-fastq/${SAMPLE}/split-adapter-quality-trimmed/${SAMPLE}-READ1.fastq.gz
    IN_FQ_R2=${PROJ_DIR}/01_clean-fastq/${SAMPLE}/split-adapter-quality-trimmed/${SAMPLE}-READ2.fastq.gz
    K_REPORT=${PROJ_DIR}/02_kraken/${SAMPLE}_kraken.report
    K_TABLE=${PROJ_DIR}/02_kraken/${SAMPLE}_kraken.tbl
    
    #output files
    OUT_FQ_R1=${TAXA_DIR}/01_reads/${SAMPLE}/${TAXA}_taxid-${TAXID}_${SAMPLE}_R1.fq
    OUT_FQ_R2=${TAXA_DIR}/01_reads/${SAMPLE}/${TAXA}_taxid-${TAXID}_${SAMPLE}_R2.fq
    
    mkdir -p 01_reads/${SAMPLE}

    #extract reads
    /master/nplatt/pathogen_probes/bin/KrakenTools/extract_kraken_reads.py \
        -k ${K_TABLE} \
        -1 ${IN_FQ_R1} \
        -2 ${IN_FQ_R2} \
        --report ${K_REPORT} \
        --output ${OUT_FQ_R1} \
        --output2 ${OUT_FQ_R2} \
        --fastq-output \
        --taxid ${TAXID} \
        --include-children

    gzip ${OUT_FQ_R1} ${OUT_FQ_R2}

done

#count the number of reads and make a new samples list
for SAMPLE in $(ls ${TAXA_DIR}/01_reads/); do

    FQ_R1=${TAXA_DIR}/01_reads/${SAMPLE}/${TAXA}_taxid-${TAXID}_${SAMPLE}_R1.fq.gz

    SAMPLE=$(basename $SAMPLE)
    FQ_LEN=$(zcat $FQ_R1 | head | wc -l | cut -f1 -d" ")
    
    if test $FQ_LEN -gt "2"; then
        echo $SAMPLE
    fi
done >${TAXA_DIR}/${TAXA}_samples.txt

#now proceed and make assembly conf file
#[samples]
#sample_id:read_dir
echo "[samples]">${TAXA_DIR}/00_confs/assembly.conf
for SAMPLE in $(cat ${TAXA_DIR}/${TAXA}_samples.txt); do

    IN_DIR="${TAXA_DIR}/01_reads/${SAMPLE}"
    echo "${SAMPLE}:${IN_DIR}"

done >>${TAXA_DIR}/00_confs/assembly.conf

phyluce_assembly_assemblo_spades \
    --conf ${TAXA_DIR}/00_confs/assembly.conf \
    --output ${TAXA_DIR}/02_spades_assemblies \
    --cores ${N_THREADS} \
    --memory 768 \
    --log-path 00_logs 

############################ COVERAGE FILTER
for FASTA in $(ls ${TAXA_DIR}/02_spades_assemblies/contigs/*.fasta); do
    python ~/pathogen_probes/code/filter_uce_assemblies.py \
        ${FASTA} \
        $TAXA_DIR/02_spades_assemblies/filtered_contigs \
        10 \
        100
done

#process to appropriate 2bit format (and rename)
genome_updater.sh \
    -d refseq,genbank \
    -T ${TAXID},${OUTGROUP_TAXID} \
    -f genomic.fna.gz \
    -A 1 \
    -o ${TAXA_DIR}/03_ncbi_genomes \
    -b v1 \
    -t ${N_THREADS} \
    -m \
    -c "reference genome","representative genome" \
    -p

#rename
mkdir -p ${TAXA_DIR}/03_ncbi_genomes/fas
for ACCESSION in $(cat ${TAXA_DIR}/03_ncbi_genomes/assembly_summary.txt | cut -f1 | sed 's/\.[0-9]//'); do

    #get genome file
    OLD_FAS=$(ls ${TAXA_DIR}/03_ncbi_genomes/v1/files/ | grep $ACCESSION)
    OLD_DIR=${TAXA_DIR}/03_ncbi_genomes/v1/files

    #new_name
    NEW_FAS=$(grep $ACCESSION ${TAXA_DIR}/03_ncbi_genomes/assembly_summary.txt \
              | awk -F'\t' '{print $8"_"$1}' \
              | sed 's/\.[0-9]$//' \
              | sed 's/ /_/g' \
              | sed 's/(/_/g' \
              | sed 's/)/_/g' \
              | sed 's/\./_/g' \
              | sed 's/:/_/g' \
              | sed 's/=/_/g' \
              | sed 's/\[\|\]//g').fa.gz

    NEW_DIR=${TAXA_DIR}/03_ncbi_genomes/fas

    echo ${NEW_DIR}/${NEW_FAS}

    cp ${OLD_DIR}/${OLD_FAS} ${NEW_DIR}/${NEW_FAS}
done

#get two bit info
for FAS in $(ls ${TAXA_DIR}/03_ncbi_genomes/fas/*.fa.gz); do

    echo $FAS

    NAME=$(basename $FAS .fa.gz)
    NAME=$(echo $NAME | sed 's/-/_/g')

    OUT_DIR="${TAXA_DIR}/03_ncbi_genomes/fas/${NAME}"

    mkdir $OUT_DIR

    TWOBIT=${NAME}.2bit
    faToTwoBit $FAS $OUT_DIR/$TWOBIT
    
    INFO=${NAME}.tab
    twoBitInfo $OUT_DIR/$TWOBIT $OUT_DIR/$INFO
    
    mv $FAS $OUT_DIR

done

#get target specific  probes
cat ~/pathogen_probes/final_probes_v1.0.fas | grep -A1 --no-group-separator -i ${PROBE_GENUS} >${TAXA_DIR}/${PROBE_GENUS}_probes.fas

#and clean up the names so they match the Faircloth format
sed -i "s/uce_[a-z,A-Z]*_/uce-/" ${PROBE_GENUS}_probes.fas
sed -i "s/uce_nematoda-clade[1-5]/uce/" ${PROBE_GENUS}_probes.fas

ulimit -n 4000

#get list of genomes
ls ${TAXA_DIR}/03_ncbi_genomes/fas/ >${TAXA_DIR}/00_confs/genomes.list

#if there are to many genomes it can cause a ulmimit error so we need to process them in batches.
#  Batches of 50 seem to work but using 25 to be safe
split -d --additional-suffix=.list --lines 25 ${TAXA_DIR}/00_confs/genomes.list ${TAXA_DIR}/00_confs/genomes_sub_

#process the first batch
phyluce_probe_run_multiple_lastzs_sqlite \
    --db ${TAXA_DIR}/04_ncbi_genomes_lastz/ncbi_genomes.sqlite \
    --output ${TAXA_DIR}/04_ncbi_genomes_lastz \
    --identity 0.85 \
    --scaffoldlist $(cat ${TAXA_DIR}/00_confs/genomes_sub_00.list) \
    --genome-base-path ${TAXA_DIR}/03_ncbi_genomes/fas \
    --probefile ${TAXA_DIR}/${PROBE_GENUS}_probes.fas \
    --cores $N_THREADS \
    --log-path ${TAXA_DIR}/00_logs/

#MODIFIED VERSION OF phyluce_probe_run_multiple_lastzs_sqlite2 that allows your to run append (by skipping deletion)
#if there are multiple batches run each and append to the db created in the prior step
for SUB_LIST in $(ls ${TAXA_DIR}/00_confs/genomes_sub_*.list | grep -v genomes_sub_00.list); do
    #echo $SUB_LIST
    phyluce_probe_run_multiple_lastzs_sqlite_append \
        --append \
        --db ${TAXA_DIR}/04_ncbi_genomes_lastz/ncbi_genomes.sqlite \
        --output ${TAXA_DIR}/04_ncbi_genomes_lastz \
        --identity 0.85 \
        --scaffoldlist $(cat $SUB_LIST) \
        --genome-base-path ${TAXA_DIR}/03_ncbi_genomes/fas \
        --probefile ${TAXA_DIR}/${PROBE_GENUS}_probes.fas \
        --cores $N_THREADS \
        --log-path ${TAXA_DIR}/00_logs/
done

echo "[scaffolds]" >${TAXA_DIR}/00_confs/genomes.conf
for TWOBIT in $(ls ${TAXA_DIR}/03_ncbi_genomes/fas/*/*2bit); do
    SAMPLE=$(basename ${TWOBIT} .2bit)
    echo $SAMPLE":"$TWOBIT
done >>${TAXA_DIR}/00_confs/genomes.conf

phyluce_probe_slice_sequence_from_genomes \
    --lastz ${TAXA_DIR}/04_ncbi_genomes_lastz \
    --conf ${TAXA_DIR}/00_confs/genomes.conf \
    --flank 1000 \
    --name-pattern "${PROBE_GENUS}_probes.fas_v_{}.lastz.clean" \
    --output ${TAXA_DIR}/05_ncbi_genomes_uce_fasta \
    --log-path ${TAXA_DIR}/00_logs/

mkdir ${TAXA_DIR}/06_all_loci
cp ${TAXA_DIR}/05_ncbi_genomes_uce_fasta/*.fasta ${TAXA_DIR}/06_all_loci
cp ${TAXA_DIR}/02_spades_assemblies/filtered_contigs/*.contigs.fasta ${TAXA_DIR}/06_all_loci
rename .contigs.fasta .fasta ${TAXA_DIR}/06_all_loci/*.contigs.fasta

#find probed regions in assemblies
phyluce_assembly_match_contigs_to_probes \
    --contigs ${TAXA_DIR}/06_all_loci \
    --probes ${TAXA_DIR}/${PROBE_GENUS}_probes.fas \
    --min-identity 90 \
    --output ${TAXA_DIR}/07_${TAXA}_uce_search_results \
    --log-path ${TAXA_DIR}/00_logs 

echo "[all]" >${TAXA_DIR}/00_confs/taxon-set.conf
for FAS in $(ls ${TAXA_DIR}/06_all_loci/*fasta); do
    NAME=$(basename ${FAS} .fasta)
    echo $NAME
done >>${TAXA_DIR}/00_confs/taxon-set.conf

mkdir -p ${TAXA_DIR}/08_taxon_sets/

# create the data matrix configuration file
phyluce_assembly_get_match_counts \
    --locus-db ${TAXA_DIR}/07_${TAXA}_uce_search_results/probe.matches.sqlite \
    --taxon-list-config ${TAXA_DIR}/00_confs/taxon-set.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output ${TAXA_DIR}/08_taxon_sets/all-taxa-incomplete.conf \
    --log-path ${TAXA_DIR}/00_logs 

# get FASTA data for taxa in our taxon set
phyluce_assembly_get_fastas_from_match_counts \
    --contigs ${TAXA_DIR}/06_all_loci  \
    --locus-db ${TAXA_DIR}/07_${TAXA}_uce_search_results/probe.matches.sqlite \
    --match-count-output ${TAXA_DIR}/08_taxon_sets/all-taxa-incomplete.conf \
    --output ${TAXA_DIR}/08_taxon_sets/all-taxa-incomplete.fasta \
    --incomplete-matrix ${TAXA_DIR}/08_taxon_sets/all-taxa-incomplete.incomplete \
    --log-path ${TAXA_DIR}/00_logs 

phyluce_assembly_explode_get_fastas_file \
    --input ${TAXA_DIR}/08_taxon_sets/all-taxa-incomplete.fasta \
    --output ${TAXA_DIR}/09_raw_uce_msa_taxon \
    --by-taxon 

phyluce_align_seqcap_align \
    --input ${TAXA_DIR}/08_taxon_sets/all-taxa-incomplete.fasta \
    --output ${TAXA_DIR}/10_uce_msa \
    --taxa $(ls ${TAXA_DIR}/09_raw_uce_msa_taxon | wc -l) \
    --aligner muscle \
    --cores $N_THREADS \
    --incomplete-matrix \
    --output-format fasta \
    --log-path ${TAXA_DIR}/00_logs 

#make sure the names are short enough
sed -i 's/\(>.\{47\}\).*\(.\{13\}\)/\1_\2/' ${TAXA_DIR}/10_uce_msa/*.fasta

phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments ${TAXA_DIR}/10_uce_msa \
    --output ${TAXA_DIR}/11_uce_msa_trimmed \
    --cores ${N_THREADS} \
    --log-path ${TAXA_DIR}/00_logs 

phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments ${TAXA_DIR}/10_uce_msa \
    --output ${TAXA_DIR}/11_uce_msa_trimmed \
    --cores ${N_THREADS} \
    --log-path ${TAXA_DIR}/00_logs 

phyluce_align_remove_locus_name_from_files \
    --alignments ${TAXA_DIR}/11_uce_msa_trimmed \
    --output ${TAXA_DIR}/12_uce_phylo_prep \
    --log-path ${TAXA_DIR}/00_logs 

mkdir 13_count_loci_per_taxa
phyluce_align_get_taxon_locus_counts_in_alignments \
    --alignments  ${TAXA_DIR}/12_uce_phylo_prep \
    --input-format nexus \
    --output ${TAXA_DIR}/13_count_loci_per_taxa/counts.csv \
    --log-path ${TAXA_DIR}/00_logs 


awk -F"," '{if ($2<3) print $1}' ${TAXA_DIR}/13_count_loci_per_taxa/counts.csv >${TAXA_DIR}/13_count_loci_per_taxa/samples_to_remove.txt

NUM_TO_REMOVE=$(wc -l ${TAXA_DIR}/13_count_loci_per_taxa/samples_to_remove.txt | cut -f1 -d" ")

if [ $NUM_TO_REMOVE -gt 1 ]; then
    phyluce_align_extract_taxa_from_alignments \
    --alignments  ${TAXA_DIR}/12_uce_phylo_prep \
    --input-format nexus \
    --output ${TAXA_DIR}/14_filter_taxa_by_n_loci \
    --output-format nexus \
    --exclude $(tr '\n' ' ' <${TAXA_DIR}/13_count_loci_per_taxa/samples_to_remove.txt) 

    else
        cp -r ${TAXA_DIR}/12_uce_phylo_prep ${TAXA_DIR}/14_filter_taxa_by_n_loci

fi

phyluce_align_get_taxon_locus_counts_in_alignments \
    --alignments  ${TAXA_DIR}/14_filter_taxa_by_n_loci \
    --input-format nexus \
    --output ${TAXA_DIR}/14_filter_taxa_by_n_loci/counts.csv \
    --log-path ${TAXA_DIR}/00_logs 

NUM_SAMPLES=$(grep -v "taxon,count" ${TAXA_DIR}/14_filter_taxa_by_n_loci/counts.csv | wc -l)

phyluce_align_get_only_loci_with_min_taxa \
    --alignments ${TAXA_DIR}/14_filter_taxa_by_n_loci/ \
    --taxa $NUM_SAMPLES \
    --percent 0.50 \
    --output ${TAXA_DIR}/15_uce_min_taxa_filt \
    --log-path ${TAXA_DIR}/00_logs 

phyluce_align_concatenate_alignments \
    --alignments ${TAXA_DIR}/15_uce_min_taxa_filt \
    --output ${TAXA_DIR}/16_concat_nexus \
    --phylip \
    --log-path ${TAXA_DIR}/00_logs

mkdir -p ${TAXA_DIR}/17_raxml

raxml-ng \
    --all \
    --prefix ${TAXA_DIR}/17_raxml/${TAXA} \
    --seed 12345 \
    --msa  ${TAXA_DIR}/16_concat_nexus/16_concat_nexus.phylip \
    --msa-format PHYLIP \
    --data-type DNA \
    --model GTR+G \
    --tree pars{100} \
    --bs-trees 1000 \
    --threads $N_THREADS

nw_ed  ${TAXA_DIR}/17_raxml/${TAXA}.raxml.support 'i & b<=50' o >${TAXA_DIR}/17_raxml/${TAXA}.raxml.support-BS50.tre

cp ${TAXA_DIR}/17_raxml/${TAXA}.raxml.support-BS50.tre  ${TAXA_DIR}/${TAXA}-BS50.tre