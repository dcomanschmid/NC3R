########################################################################
# Agilent microarray probe sequence mapping to reference transcriptomes#
#	>> personal communication & use scripts from                       #
#                https://github.com/MWSchmid/microarray                #    
#                prepareMicroarrayProbes.py                            #
# Diana Coman Schmid                                                   #
# Eawag 2015                                                           #
# diana.comanschmid@eawag.ch                                           #
########################################################################

    # map probe sequences from the zebrafish Agilent microarray platform to the latest zebrafish transcriptome (Zv9)
    # build the bowtie index 

python ./MicroarrayCrossPlatform/prepareMicroarrayProbes.py BUILD Danio_rerio.Zv9.cdna.all.fa Danio_rerioZv9cdna_index

    # get probe ID and nucleotide sequence information
    # remove the first 9 lines (not needed) and select only the "ProbeUID" and "Sequence" columns
    
head Cy3_AA_High_Cy5_AA_Lind_Ctrl_253100010016_2_2.txt

sed '1,9d' Cy3_AA_High_Cy5_AA_Lind_Ctrl_253100010016_2_2.txt > zebrafish_seq_probe_tmp.txt

cut -f 10,13 zebrafish_seq_probe_tmp.txt > zebrafish_seq_probe.txt
awk '{ print $2 " " $1}' zebrafish_seq_probe.txt > zebrafish_probe_seq.txt

    # check and format (TAB delim.) the file fields 
    
python ./MicroarrayCrossPlatform/prepareMicroarrayProbes.py TABTOFASTA zebrafish_seq_probe.txt 1 2 1 zebrafish_seq_probe.fasta

    # align the probes to the reference transcriptome ("cDNA")
    
python ./MicroarrayCrossPlatform/prepareMicroarrayProbes.py ALIGN ./zebrafish/genome/Danio_rerioZv9cdna_index zebrafish_seq_probe.fasta unaligned.txt aligned.txt

    # extract the ID mappings
    
python ./MicroarrayCrossPlatform/prepareMicroarrayProbes.py EXTRACT aligned.txt probeNameToID.txt IDtoProbeName.txt
