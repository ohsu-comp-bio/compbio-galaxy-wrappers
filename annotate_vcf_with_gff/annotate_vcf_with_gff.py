import argparse

def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('gff', help='')
    parser.add_argument('vcf', help='')
    parser.add_argument('gene_list', help='')
    parser.add_argument('--ordered_test', '-o')
    args = parser.parse_args()
    return args

def filter_genes(genelist_filepath, ot_input):

    gene_list = open(genelist_filepath, 'r')
    gene_filter = []

    match = False
    for line in gene_list:
        line_array = line.split('\t')
        ordered_test = line_array[1]

        if ot_input == ordered_test:
            match = True
            gene_filter.append(line_array[2])

    if match == False:
        raise NameError('Please enter proper ordered test name.')

    return gene_filter


def annotate_vcf(gff_filepath, vcf_filepath, genes=[]):
    gff = open(gff_filepath, 'r')
    vcf = open(vcf_filepath, 'r')
    new_vcf_filepath = vcf_filepath.partition('.vcf')[0] + '_MODIFIED' + vcf_filepath.partition('.vcf')[1]
    new_vcf = open(new_vcf_filepath, 'w')

    # Init dictionary of GENE:[chr, start pos, end pos] from GFF
    gene_dict = {}

    # loop over each of the records in the vcf
    for line in gff:
        if line.startswith('#'):
            continue
        else:
            # Split the line into workable entries
            line_array = line.split('\t')
            if line_array[2] == 'region':
                info = line_array[8].split(';')[3]
                current_chr = info.partition('chromosome=')[2]
            elif line_array[2] == 'gene':
                start = line_array[3]
                end = line_array[4]
                gene = line_array[8].split(';')[1].partition('Name=')[2]
                gene_dict[gene] = [current_chr, start, end]

    for line in vcf:
        # If it's a header line, just write it to new VCF
        if line.startswith('#'):
            new_vcf.write(line)
            continue
        else:
            # Get INFO entry
            line_array2 = line.split('\t')
            info_array = line_array2[7].split(';')

            # iterate through .gff genes
            for k, v in gene_dict.items():
                if v[0] == line_array2[0] and line_array2[1] >= v[1] and line_array2[1] <= v[2]:

                    if k not in genes and len(genes) != 0:
                        break
                    else:
                        info_array.append('Gene=' + k)
                        break

            # Filter out non-annotated lines by skipping it if gene info isn't available
            if info_array[-1][0:4] != 'Gene' and len(genes) != 0:
                continue

            first = True
            for entry in info_array:
                if first:
                    corrected_info = entry
                    first = False
                else:
                    corrected_info += ';' + entry
            line_array2[7] = corrected_info

            # write the updated line array to the new vcf
            first = True
            for entry in line_array2:
                if first:
                    new_vcf.write(entry)
                    first = False
                else:
                    new_vcf.write('\t' + entry)

    vcf.close()
    new_vcf.close()


def main():

    args = supply_args()

    if args.ordered_test != None:
        genes = filter_genes(args.gene_list, args.ordered_test)
        annotate_vcf(args.gff, args.vcf, genes)
    else:
        annotate_vcf(args.gff, args.vcf)

if __name__ == '__main__':
    main()