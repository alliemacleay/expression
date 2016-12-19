BAMS=rna_bams_list_121416.txt
COVERAGE=rna_coverage_121416.txt
echo "chrom   assay   gene    start   stop    n1      strand  n2      gene_id " > head.txt 
sed -e 's/ /\t/g' $BAMS >> head.txt 
cat $COVERAGE >> head.txt 
mv head.txt $COVERAGE
