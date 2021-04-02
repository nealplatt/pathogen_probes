SAMPLE=$1

#read in accession
ACCESSION=$(echo $SAMPLE | cut -f3,4 -d "_")

NEW_GENOME=$SAMPLE".fas"
GFF=$SAMPLE".gff"

chmod u+w *

###################################
# mask genome

#mask lowercased from NCBI
cat $NEW_GENOME | sed '/^>/! s/[atcg]/N/g' >ncbi_masked_genome.fas

#repeatmasker - mammals
~/bin/RepeatMasker/RepeatMasker \
    -pa 12 \
    -s \
    -species mammal \
    ncbi_masked_genome.fas

#repeatmasker - species
SPECIES=$(echo $SAMPLE | cut -f1,2 -d "_" | sed 's/_/ /')
GENUS=$(echo $SPECIES | cut -f1 -d"_")
MASKED_GENOME=$SAMPLE.fullymasked.fas

~/bin/RepeatMasker/RepeatMasker \
    -pa 12 \
    -s \
    -species "$SPECIES" \
    ncbi_masked_genome.fas.masked

SPECIES_RESULT=$?

if [ $SPECIES_RESULT != 0 ]; then
  ~/bin/RepeatMasker/RepeatMasker \
    -pa 12 \
    -s \
   -species "$GENUS" \
    ncbi_masked_genome.fas.masked
    
    GENUS_RESULT=$?

    if [ $GENUS_RESULT != 0 ]; then
      ln -s ncbi_masked_genome.fas.masked $MASKED_GENOME
    fi

else
    ln -s ncbi_masked_genome.fas.masked.masked $MASKED_GENOME
fi


###################################
# get 200 bp windows

awk '{ if ($3 == "CDS") print $0}' $GFF >cds.gff
bedtools makewindows -w 120 -b cds.gff >cds_120bp_windows.bed
bedtools getfasta -fi $MASKED_GENOME -fo cds_120bp_windows.fas -bed cds_120bp_windows.bed

###################################
# filter out probes with any repeats

python ~/pathogen_probes/code/remove_repetitive_targets.py cds_120bp_windows.fas #outputs repeat_cleaned_seqs.fas

###################################
# blast 200bp windows against genome

#make blast index
makeblastdb -in $NEW_GENOME -dbtype nucl

#blast
blastn \
    -query repeat_cleaned_seqs.fas \
    -db $NEW_GENOME \
    -out cds_120bp_windows_v_genome.blastn.out \
    -evalue 1e-10 \
    -num_threads 4 \
    -outfmt 6


#sort for single hits
cut -f1 cds_120bp_windows_v_genome.blastn.out \
    | sort \
    | uniq -c \
    | awk '{ if ($1 == 1) print $2}' \
    >single_copy_windows.txt

#make sure that the probed region is 120bp
sed 's/:/\t/' single_copy_windows.txt \
    | sed 's/-/\t/' \
    | awk '{ if ($3-$2 == 120) print $0}' \
    >single_copy_windows_120bp.bed

######################################
# make sure that sites aren't to close

# index genome
samtools faidx $NEW_GENOME

#find 5K random sites in the genome
bedtools random \
    -l 1 \
    -n 5000 \
    -seed 12345 \
    -g $NEW_GENOME.fai \
    | bedtools sort \
        -i - \
    >random_5k.bed

#sort the put. probe file and...
bedtools sort \
    -i single_copy_windows_120bp.bed \
    >sorted_single_copy_windows_120bp.bed

#and find the closest probe to those random regions
bedtools closest \
    -t first \
    -a random_5k.bed \
    -b sorted_single_copy_windows_120bp.bed \
    | cut -f1,2,3 \
    | shuf \
    >distributed_5k.bed

#reformat headers
python ~/pathogen_probes/code/dispersed_probes.py distributed_5k.bed $SAMPLE #outputs dispersed_tiled_and_reformatted_probes.bed

######################################
# get fasta of probes
bedtools \
    getfasta \
    -nameOnly \
    -fi $NEW_GENOME \
    -bed dispersed_tiled_and_reformatted_probes.bed \
    -fo "$SAMPLE"-tiled_probes.fas

