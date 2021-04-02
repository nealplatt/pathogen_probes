import sys

in_bed = sys.argv[1]
sample_info = sys.argv[2]

genus     = sample_info.split("_")[0]
species   = sample_info.split("_")[1]
accession = sample_info.split("_")[-2:]
accession = "_".join(accession)

clusters_seen = []
unique_beds_per_cluster = []


with open("dispersed_tiled_and_reformatted_probes.bed", 'w') as out_file:

    with open(in_bed, 'r') as in_file: 
        bed_entries = in_file.readlines() 
        
        cluster_id = 0
        for bed_entry in bed_entries:
            chrom, start, stop = bed_entry.rstrip().split("\t")

            #tile
            tile_1=[start,   int(start) + 80]
            tile_2=[int(stop)-80, stop]
            
            
            
            tile_1_header = "pp_{}_{}_{}_p1 | species=\"{} {}\"; accession={}; coords={}:{}-{}".format( genus,
                                                                                                        species, 
                                                                                                        cluster_id, 
                                                                                                        genus, 
                                                                                                        species,
                                                                                                        accession, 
                                                                                                        chrom, 
                                                                                                        tile_1[0], 
                                                                                                        tile_1[1] )

            tile_2_header = "pp_{}_{}_{}_p2 | species=\"{} {}\"; accession={}; coords={}:{}-{}".format( genus,
                                                                                                        species, 
                                                                                                        cluster_id, 
                                                                                                        genus, 
                                                                                                        species, 
                                                                                                        accession, 
                                                                                                        chrom, 
                                                                                                        tile_2[0], 
                                                                                                        tile_2[1] )

            outline_tile_1 = "\t".join([chrom, start, stop, tile_1_header])
            outline_tile_2 = "\t".join([chrom, start, stop, tile_2_header])

            out_file.write(outline_tile_1 + "\n" + outline_tile_2 + "\n")
            cluster_id = cluster_id + 1
