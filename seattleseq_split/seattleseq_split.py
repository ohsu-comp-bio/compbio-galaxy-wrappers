#!/usr/bin/env python

### Split a writeGenotypes carrierAll.txt file by multiallelic sites and write two new VCF's.
### These will then be processed back through SeattleSeq and prepped to be sent to CGD.
### ./seattleseq_split.py <carrierAll.txt> <main VCF> <alt VCF>


import sys

def writeVcfHeader(vcf):

    vcf.write("##fileformat=VCFv4.1\n")
    vcf.write("##FILTER=<ID=DPFilter,Description=\"DP < 10\">\n")
#    vcf.write("##FILTER=<ID=LowQual,Description=\"Low quality\">\n")
    vcf.write("##FILTER=<ID=QDFilter,Description=\"QD < 5\">\n")
    vcf.write("##FILTER=<ID=QUALFilter,Description=\"QUAL < 100\">\n")
    vcf.write("##FILTER=<ID=SnpCluster,Description=\"SNPs found in clusters\">\n")
    vcf.write("##reference=file:///Users/letaw/Local/resources/genomes/hg19/Homo_sapiens_assembly19.fasta\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

def main():

    """
    << carrierAll.txt header >>
    #chrom
    positionHg19
    type
    referenceBase
    alternateBase
    filterFlagGATK
    QUAL
    PG0003583-BLDGtype
    PG0003583-BLDDepth
    PG0003583-BLDQual
    geneList
    rsID
    createBuildfunctionGVSaminoAcidsproteinPositioncDNAPositiongenomeHGVSgranthamScorescorePhastConsconsScoreGERPscoreCADDpolyPhenSIFTESPEurAlleleCountsESPEurMinorPercentESPAfrAlleleCountsESPAfrMinorPercent1000GenomesEur%1000GenomesAfr%1000GenomesAsn%ExACExomesEur%ExACExomesAfr%ExACExomesAsn%local507clinicalAssnOMIMClinVarphenotypeHGMDreferenceHGMD
    """

    CHROMS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    handle_in = open(sys.argv[1], 'rU')
    handle_out = open(sys.argv[2], 'w')
    handle_out_alt = open(sys.argv[3], 'w')

    writeVcfHeader(handle_out)
    writeVcfHeader(handle_out_alt)

### 1    12783    .    G    A    115.85    PASS    AC=2;AF=1.00;AN=2;DP=6;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=26.68;MQ0=0;QD=19.31    GT:AD:DP:GQ:PL    1/1:0,6:6:17:144,17,0

    
    with handle_in as seattleseq:

        next(seattleseq)

        format = "GT:AD:DP:GQ:PL"

        for line in seattleseq:

            alt1 = None

            line = line.rstrip('\n').split('\t')
            
            ### These are the same whether the site is multiallelic or not.

            chrom = line[0]
            position = line[1]
            filter = line[5]
            qual = line[6]
            ref = line[3]
            gq = line[9]
            pl = "."
            af = "0.5"
            fs = "0.000"

            if line[11] != "0":
                dbsnp = "rs" + line[11]
            else:
                dbsnp = "."

            ### These need to be handled if we run in to a multiallelic site.
                
            if ',' not in line[4]:
                alt = line[4]
                ad = line[8].split('(')[1].split(')')[0]
                info_dp = line[8].split('(')[0]
                dp = str(int(line[8].split('(')[1].split(')')[0].split(',')[0]) + int(line[8].split('(')[1].split(')')[0].split(',')[1]))
            else:
                alt = line[4].split(',')[0]
                alt1 = line[4].split(',')[1]
                ad = line[8].split('(')[1].split(')')[0].split(',')[0] + ',' + line[8].split('(')[1].split(')')[0].split(',')[1]
                ad1 = line[8].split('(')[1].split(')')[0].split(',')[0] + ',' + line[8].split('(')[1].split(')')[0].split(',')[2]
                info_dp = str(int(line[8].split('(')[1].split(')')[0].split(',')[0]) + int(line[8].split('(')[1].split(')')[0].split(',')[1]) + int(line[8].split('(')[1].split(')')[0].split(',')[2]))
                dp = str(int(line[8].split('(')[1].split(')')[0].split(',')[0]) + int(line[8].split('(')[1].split(')')[0].split(',')[1]))
                dp1 = str(int(line[8].split('(')[1].split(')')[0].split(',')[0]) + int(line[8].split('(')[1].split(')')[0].split(',')[2]))
                qd1 = str("%.2f" % (float(qual)/float(dp)))
                info1 = "AF=" + af + ";DP=" + info_dp + ";FS=" + fs + ";QD=" + qd1
                sample1 = ':'.join(["0/1", ad1, dp1, gq, pl])

            qd = str("%.2f" % (float(qual)/float(dp)))
            sample = ':'.join(["0/1", ad, dp, gq, pl])
            info = "AF=" + af + ";DP=" + info_dp + ";FS=" + fs + ";QD=" + qd


            if chrom in CHROMS:
                to_write = "\t".join([chrom, position, dbsnp, ref, alt, qual, filter, info, format, sample])
                if alt1 != None:
                    to_write_alt = "\t".join([chrom, position, dbsnp, ref, alt1, qual, filter, info1, format, sample1])
                    handle_out_alt.write(to_write_alt + '\n')

                handle_out.write(to_write + '\n')

            else:
                # The chromosome must currently be either 1-22, X, or Y.  Should add support for MT soon.
                pass

    handle_out.close()
    handle_out_alt.close()

main()
