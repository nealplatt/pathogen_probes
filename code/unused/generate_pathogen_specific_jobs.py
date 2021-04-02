import os                                                                                                                                                                                                 
import shutil                                                                                                                                                                                             
import sys                                                                                                                                                                                                
                         
os.chdir('/master/nplatt/pathogen_probes')                                                                                                                                                                                                     
                                                                                                                                                                                 
#accessions_csv = sys.argv[1]                                                                                                                                                                             
#accessions_csv = "data/pathogen_specific_probe_accessions.csv"
accessions_csv = "data/tmp.csv"

targeted_parasites = []
with open(accessions_csv, 'r') as in_file:                                                                                                                                                                
    csv_entries = in_file.readlines()                                                                                                                                                                     
    for csv_entry in csv_entries:                                                                                                                                                                         
        if csv_entry != "genus,species,accession":
            targeted_parasites.append(csv_entry.rstrip())


parasites_with_gff = []
parasites_wo_gff = []

for parasite in targeted_parasites:                                                                                                                                                  
    genus, species, accession = parasite.split(",")                                                                                                                                              
                                                                                                                                                                                                              
    sample_name="{}_{}_{}".format(genus, species, ''.join(accession.split())[:-2])                                                                                                                
                                                                                                                                                                                                          
    o_dir = "results/pathogen_specific_probes/" + sample_name                                                                                                                                     
    if os.path.exists(o_dir):                                                                                                                                                                     
        shutil.rmtree(o_dir)                                                                                                                                                                      
                                                                                                                                                                                                          
    os.makedirs(o_dir)                                                                                                                                                                            
    os.chdir(o_dir)                                                                                                                                                                               
                                                                                                                                                                                                          
    #get the ftp site from the accession                                                                                                                                                          
    ftp=!(esearch -db assembly -query {accession} | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank)                                                                          
    ftp=ftp[0].replace("ftp://", "")                                                                                                                                                              
            
    #generate a fas and download                                                                                                                                                                  
    fas_ftp_link  = ftp + "/" + ftp.split("/")[-1] + "_genomic.fna.gz"                                                                                                                            
    fas_rsync_cmd ="rsync --copy-links --times --verbose rsync://{} .".format(fas_ftp_link, sample_name)                                                                                                                                                                                                                                                                   
    !{fas_rsync_cmd}                                                                                                                                                                              
                                                                                                                                                                                                          
    #generate a gff and download                                                                                                                                                                  
    gff_ftp_link  = ftp + "/" + ftp.split("/")[-1] + "_genomic.gff.gz"                                                                                                                            
    gff_rsync_cmd ="rsync --copy-links --times --verbose rsync://{} .".format(gff_ftp_link, sample_name)                                                                                          
    !{gff_rsync_cmd}                                  
            
    clean_fas_name = "_".join([genus, species, accession.split(".")[0]]) + ".fas.gz"
    orig_fas_name  = ftp.split("/")[-1] + "_genomic.fna.gz"
    os.rename(orig_fas_name, clean_fas_name)
                        
    clean_gff_name = "_".join([genus, species, accession.split(".")[0]]) + ".gff.gz"
    orig_gff_name  = ftp.split("/")[-1] + "_genomic.gff.gz"
    if not os.path.exists(orig_gff_name):
        print(parasite + " no gff")
        os.chdir('/master/nplatt/pathogen_probes')
        shutil.rmtree(o_dir) 
        parasites_wo_gff.append(parasite)
        continue
    else:
        os.rename(orig_gff_name, clean_gff_name)
        parasites_with_gff.append(parasite)
        !gunzip *.gz
        qsub = "qsub -V -cwd -S /bin/bash -q all.q -j y -N {} -o {}.log -pe smp 12".format(sample_name, sample_name)
        cmd = '"conda activate pathogen_probes; ~/pathogen_probes/code/single_genome_probes.sh {}"'.format(sample_name)

        !echo {cmd} | {qsub}
        os.chdir('/master/nplatt/pathogen_probes')


with open('results/pathogen_specific_probes/no_gffs.csv', 'w') as out_file:
    out_file.write("genus,species,accession\n")
    for parasite_entry in parasite_entries:                                                                                                                                                                         
        out_file.write(parasite_entry + "\n")


