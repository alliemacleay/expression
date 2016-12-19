Expression
==========
Calculate RNA to DNA ratio to quantify expression and plot the values along with the fusion information found

Required files
--------------
Indexed BAM files
targets.gtf
exome locations - gencodeV24lift37.bed
Fusion information file

Some of the file names are hardcoded so you will have to go into mk_header.sh and expression.R to edit the hardcoded names.

Steps
-----
1. Copy all bams and their index files to a directory
2. ls bam_dir > bam_files.txt
3. sed -e 's/\n/ /g' bam_files.txt > bam_files.oneline.txt
4. bedtools multicov -abam \`cat bam_files.oneline.txt\` -b targets.gtf > coverage.txt
5. ./mk_header.sh
6. expression_allie.R


