Expression
==========
Calculate RNA to DNA ratio to quantify expression and plot the values along with the fusion information found

Allison MacLeay, M.S.

Zongli Zheng, Ph.D.

Richard Schmidt, M.D., Ph.D.

![alt-text](https://github.com/alliemacleay/expression/blob/master/target_creation.png "Exon and intron targets")

Choose targets at the exon closest to the primer 3 bases before the end of the exon and 3 bases after into to intron to measure rna coverage and dna coverage.

Required files
--------------
* Indexed BAM files
* primer locations - targets.gtf
* exome locations - gencodeV24lift37.bed
* Fusion information file

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
5. Calculate coverage of all targets for all samples (single quote is a backtick)
  * `bedtools multicov -bams 'cat bam_files.oneline.txt' -bed expression.bed > coverage.txt`
6. Fix the coverage file so that it has a header
  * `./mk_header.sh`
7. Run analysis
  * `expression_allie.R`
  
Output
--------
RNA_to_DNA_.raw.full.pdf  - box plot with all genes and labeled fusions

RNA_to_DNA_\<gene\>.raw.pdf – box plot by gene with labeled fusions

RNA_to_DNA_\<gene\>_exp.pdf – scatter plot with log10 ratio of rna to dna expression with labeled fusions

RNA_to_DNA_normalized_by_housekeeping.full.pdf – box plot with all genes and labeled fusions normalized by MSMB, ERG, and TFEB

RNA_to_DNA_normalized_by_housekeeping.pdf – subset of genes (ALK, BRD4, and ROS1)

RNA_normalized_by_housekeeping.full.pdf – RNA expression normalized by MSMB, ERG, and TFEB

RNA_normalized_by_housekeeping_subset.pdf - subset of genes (ALK, BRD4, and ROS1) RNA expression


