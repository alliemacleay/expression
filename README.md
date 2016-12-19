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
1. Create a target file from the targets file and an exon reference.
  * `python multicov_bed_file_generator.py --exome-file gencodeV24lift37.bed --primer-bed targets.gtf --output-file expression.bed`
2. Copy all bams and their index files to a directory
3. Make a file with the bam paths
  * `ls bam_dir > bam_files.txt`
4. Put all bam paths on one line
  * `sed -e 's/\n/ /g' bam_files.txt > bam_files.oneline.txt`
5. Calculate coverage of all targets for all samples
  * `bedtools multicov -abam \`cat bam_files.oneline.txt\` -b expression.bed > coverage.txt`
6. Fix the coverage file so that it has a header
  * `./mk_header.sh`
7. Run analysis
  * `expression_allie.R`


