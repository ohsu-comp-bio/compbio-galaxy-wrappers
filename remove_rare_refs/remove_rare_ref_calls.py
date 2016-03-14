#!/usr/bin/env python

### Remove the rare reference calls from SeattleSeq tsv writeGenotypes output.
### These should not be sent to CGD, they are not called in the VCF.

import argparse

RARE_REFS = [("1","169519049"),("7","6026775"),("13","32929387"),("14","75513883")]
FAKE_DEPTH = "100(100,0)"

def compareCoords(chrom, pos):

    if (chrom, pos) in RARE_REFS:
        return True
    else:
        return False


def parseSeattleSeqTSV(handle, handle_out):
    
    """
    #chrompositionHg19typereferenceBasealternateBasefilterFlagGATKQUALSAMPLEGtypeSAMPLEDepthSAMPLEQualgeneListrsIDcreateBuildfunction\
    GVSaminoAcidsproteinPositioncDNAPositiongenomeHGVSgranthamScorescorePhastConsconsScoreGERPscoreCADDpolyPhenSIFTESPEurAlleleCounts\
    ESPEurMinorPercentESPAfrAlleleCountsESPAfrMinorPercent1000GenomesEur%1000GenomesAfr%1000GenomesAsn%ExACExomesEur%ExACExomesAfr%\
    ExACExomesAsn%local507clinicalAssnOMIMClinVarphenotypeHGMDreferenceHGMD
    167504231indelTATPASS774.19TA/T17(1,16)99SLC35D111335108120intronnoneNANANC_000001.10:g.67504232
    _67504232delANA0.457-100.0-100.0unknownunknownNA-100.0NA-100.0unknownunknownunknownunknownunknownunknownNAunknownnononenonenone

    These are tab-separated.
    """

    with handle as seattle_seq:
        for line in seattle_seq:

            sline = line.rstrip('\n').split('\t')
            chrom = sline[0]
            pos = sline[1]
            depth = sline[8]

            if compareCoords(chrom, pos) == False and depth != FAKE_DEPTH:
                handle_out.write(line)
            else:
                print("Coordinate " + chrom + ":" + pos + " found in input, removing.")
    
    handle_out.close()


def main():
    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(dest='input', help='Input SeattleSeq TSV.')
    parser.add_argument(dest='output', help='Output SeattleSeq TSV.')
    args = parser.parse_args()

    handle_in_tsv = open(args.input, 'rU')
    handle_out_tsv = open(args.output, 'w')

    parseSeattleSeqTSV(handle_in_tsv, handle_out_tsv)

if __name__ == "__main__":
    main()
