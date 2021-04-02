conda activate phyluce

cd ~/patho_detect/results

mkdir -p postseq_phyluce/raw-fastq

cd postseq_phyluce/raw-fastq

#get raw data
cp ~/patho_detect/data/19047-23/*.gz .


#process reads with illumina processor -> clean_fastq/

#make illuminaprocessor.conf file

illumiprocessor \
    --input raw-fastq/ \
    --output clean-fastq \
    --config illumiprocessor.conf \
    --cores 24


#assemble reads
phyluce_assembly_assemblo_spades \
    --conf assembly.conf \
    --output spades-assemblies \
    --cores 48 \
    --memory 768


#find UCE loci

#get probes
cp ../../decon_probes.fas .

#find UCE loci
phyluce_assembly_match_contigs_to_probes \
    --contigs spades-assemblies/contigs \
    --probes decon_probes.fas \
    --output uce-search-results







