#!/usr/bin/env python
import argparse
import pysam

#############
# author: skadri
# This script will parse a bam file and remove reads that:
# 1. are not uniquely mapped to the reference

def main(inputbam,outputbam):

    bamfile = pysam.AlignmentFile(inputbam, "rb")
    outputfile = pysam.AlignmentFile(outputbam, "wb", template=bamfile)
    supplementary =0
    secondary =0
    unmapped =0
    output = 0
    mapq = 0

    for read in bamfile.fetch():
        if read.is_secondary == True:
            secondary += 1
        else:
            if read.is_supplementary == True:
                supplementary += 1
            else:
                if read.is_unmapped == True:
                    unmapped += 1
                else:
                    if read.mapping_quality==0:
                        mapq +=1
                    else:
                        output += 1
                        outputfile.write(read)

    print("secondary={} \nsupplementary={} \nunmapped={} \nmapping quality 0={}\noutput={}".format(secondary, supplementary, unmapped,mapq,output))

    # close the bam files
    bamfile.close()
    outputfile.close()

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input","-i",help="Input BAM file",type=str)
    parser.add_argument("--out","-o",help="Output BAM file",type=str)

    args = parser.parse_args()

    main(args.input,args.out)
