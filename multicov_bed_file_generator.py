'''
Ryan Schmidt
Script to generate a bed file that will serve as an input to bedtools multicov
This bed file should take a bed file of GSP2s and output a line in the bed file for a location adjacent to the primer
in the intron and another in the exon

INPUTS
1. gencode basic file from ucsc
2. fusion gsp2 bed file

OUTPUT
bed file

'''

import click
from collections import defaultdict


def closest_position(position_list, primer_end, strand):
    distance_list = sorted([int(x) for x in position_list])
    if strand == '+':
        for i, pos in enumerate(distance_list):
            if int(pos) - primer_end > 0:
                break
        next = distance_list[i + 1] if i < len(distance_list) - 1 else None
        if next is not None and next == pos:
            j = 2
            while next == pos and i + j < len(distance_list):
                #print 'adjusting for +'
                next = distance_list[i + j]
                j += 1
            next = None if next == pos else next

        return int(pos), next
    elif strand == '-':
        last_pos = 0
        for i, pos in enumerate(distance_list):
            if int(pos) - primer_end > 0:
                break
            last_pos = int(pos)
        next = int(distance_list[i - 2]) if i > 1 else None
        j = 3
        while next == last_pos and i + j >= 0:
            #print 'adjusting for -'
            next = distance_list[i - j]
            j -= 1
        next = None if next == last_pos else next

        return last_pos, next  # next was for Zongli's gdna value which was in the adjacent exon, but I removed this


def get_exome_dict(exome_file):
    # parse the gencode file and record

    exon_start = defaultdict(list) # dictionary with a list as the value
    exon_end = defaultdict(list) # dictionary with a list as the value
    exon_start_end_dict = defaultdict(list) # dictionary with a list as the value
    cds_start_end_list = []
    with open(exome_file, 'rb') as gencode_basic:
        for line in gencode_basic:
            if line.startswith('#'):
                pass
            else:
                line_list = line.rstrip('\n').split('\t')
                chr = line_list[2]
                strand = line_list[3]
                cds_start = int(line_list[6])
                cds_end = int(line_list[7])
                ex_starts = list(line_list[9].split(','))
                del ex_starts[-1]
                ex_ends = list(line_list[10].split(','))
                del ex_ends[-1]
                gene_name = line_list[12]

                exon_start[chr] += ex_starts
                exon_end[chr] += ex_ends
                exon_start_end_dict[chr] += ex_starts + ex_ends # populates the dictionary by chromosome with all exon start/end coordinates in an associated list

                cds_start_end_list.extend([cds_start, cds_end])
    return exon_start_end_dict, cds_start_end_list, exon_start, exon_end


def write_to_bed(bed_outfile, chrom, primer_name, strand, intron_coordinate, exon_coordinate, next_coordinate=None):
    """
    intron - Zongli's gcorr - close by intron
    exon   - Zongli's cdna - exon
    gdna   - Zongli's gdna - intron by adjacent exon ** Removed
    """
    bed_outfile.write(chrom + '\t' + str(intron_coordinate) + '\t' + str(intron_coordinate) + '\t' + primer_name + '\t' + 'intron' + '\t' + strand + '\n')
    bed_outfile.write(chrom + '\t' + str(exon_coordinate) + '\t' + str(exon_coordinate) + '\t' + primer_name + '\t' + 'exon' + '\t' + strand + '\n')
    if next_coordinate is not None:
        bed_outfile.write(chrom + '\t' + str(next_coordinate) + '\t' + str(next_coordinate) + '\t' + primer_name + '\t' + 'gdna' + '\t' + strand + '\n')


@click.group(invoke_without_command=True)
@click.option('--primer-bed', '-p', type=click.Path(exists=True), help='Path to input VarVetter input file (required)')
@click.option('--exome-file', '-e', required=True, help='Path to output VarVetter input file (required)')
@click.option('--output-file', '-o', help='Path to bam file')
def cli(primer_bed, exome_file, output_file):
    #gencode_basic = open('/Users/r/Documents/Residency/MGP_Fellowship/MGH/Projects/MET_exon_14/expression/EncodeGencodeBasicV24lift37_RefSeq_genes')
    #primer_bed = open('/Users/r/Documents/Residency/MGP_Fellowship/MGH/Projects/MET_exon_14/expression/FusionPlex_Solid_Tumor_Panel_V1.gtf')
    #bed_outfile = open('/Users/r/Documents/Residency/MGP_Fellowship/MGH/Projects/MET_exon_14/expression/fusion_expression_multicov.bed', 'w')
    primer_bed = open(primer_bed, 'rb')
    bed_outfile = open(output_file, 'wb')

    exon_start_end_dict, cds_start_end_list, exon_start, exon_end = get_exome_dict(exome_file)

    # loop through the primer bed and create a bed outfile line for the intronic and exonic position named based on the primer name

    i = 0
    append_chr = 'chr' if 'chr' in exon_start_end_dict.keys()[0] else ''
    dist = 3  # bases between boundary and positions for intron and exon
    flip = 0
    skip = 0
    for line in primer_bed:
        print "Reading Line " + str(i)
        line_list = line.split('\t')
        chr = line_list[0]
        start = int(line_list[3])
        end = int(line_list[4])
        strand = line_list[6]
        info = line_list[8]
        primer_name = info.split('"')[1]
        if strand == '+':
            primer_end = end
        elif strand == '-':
            primer_end = start
        chr = chr.replace('chr', '')
        chr = chr.replace('chr', '')
        closest_boundary, next_exon = closest_position(exon_start_end_dict['{}{}'.format(append_chr, chr)], primer_end, strand) # find coordinate closest to primer end

        bound_type = 'start' if str(closest_boundary) in exon_start['{}{}'.format(append_chr, chr)] else 'stop'

        if closest_boundary in cds_start_end_list: # will exclude primers that face into UTRs
            print "Excluding " + primer_name + ": UTR exon"
        else:
            if strand == '+' and primer_end < closest_boundary - dist:  # added the check to make sure we're not taking primer coverage
                intron_coordinate = closest_boundary + dist
                exon_coordinate = closest_boundary - dist
                if bound_type != 'stop':
                    flip += 1
                    tmp = intron_coordinate
                    intron_coordinate = exon_coordinate
                    exon_coordinate = tmp
                write_to_bed(bed_outfile, chr, primer_name, strand, intron_coordinate, exon_coordinate)
            elif strand == '-' and primer_end > closest_boundary + dist:
                intron_coordinate = closest_boundary - dist
                exon_coordinate = closest_boundary + dist
                if bound_type != 'start':
                    flip += 1
                    tmp = intron_coordinate
                    intron_coordinate = exon_coordinate
                    exon_coordinate = tmp
                write_to_bed(bed_outfile, chr, primer_name, strand, intron_coordinate, exon_coordinate)
            else:
                print primer_name + '\t' + str(closest_boundary)
                skip += 1
                print

        # example outfile line:  7    116411552       116411552       MET_5_exon013_exon      exon    -

        i += 1
    primer_bed.close()
    bed_outfile.close()
    print '{} flipped, {} skipped'.format(flip, skip)

if __name__ == '__main__':
    cli()