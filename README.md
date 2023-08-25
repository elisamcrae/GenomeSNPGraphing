# Sars-CoV-2 Genome Analysis
2023 June-August internship

## Purpose
This project aimed to replicate NextStrain.org's covid database internally within bioMerieux. This allows the company to own the data and parse it accordingly. It also allows for more extensive analyses of current covid tests and the position of the primers in relation to where mutational SNPs occurred.

Each script has a summary at the top. The first script to run is main.py, after which different scripts can be run on the output of main.py based on desired result.

The final product of this internship project was an IVG which mapped all SNPs from GISAID sars-cov-2 genomes to the reference genome and current bioMerieux covid primers.

To split 15.5 million genomes into smaller files of 1 million each, run this script: gt splitfasta -numfiles 15 completeGenomeFile.fasta

Command to run nucmer and pipe it into main: NUCMER GISAID_REF.fasta {GISAID sequences} ; show-snps -q -T -H -l out.delta | python3 main.py GISAID_REF.fasta

Command to run any script that takes in the sparse matrices (placeholder of findNaturalDrift shown): nohup python3 ~/findNaturalDrift.py sequencesC.fasta..1.deltaSM sequencesC.fasta..2.deltaSM sequencesC.fasta..3.deltaSM sequencesC.fasta..4.deltaSM sequencesC.fasta..5.deltaSM sequencesC.fasta..6.deltaSM sequencesC.fasta..7.deltaSM sequencesC.fasta..8.deltaSM sequencesC.fasta..9.deltaSM sequencesC.fasta..10.deltaSM sequencesC.fasta..11.deltaSM sequencesC.fasta..12.deltaSM sequencesC.fasta..13.deltaSM sequencesC.fasta..14.deltaSM sequencesC.fasta..15.deltaSM &

## Files
Each .png file was made from this code from the sparseMatrix file and can be used for reference.

hundredThousandGenomes is a test file to start with via the main.py. This will create a sparseMatrix file (similar to the one already provided). Next, the sparseMatrix file of choice can be used as input for nearly every script. For the few scripts that do not take a sparse Matrix input file, this is commented at the head of the script.
