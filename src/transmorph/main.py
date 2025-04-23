import gtf_utils as ul
import argparse
import os
import sys

def main():
   
    parser = argparse.ArgumentParser(
        description='a tool that converts a BAM file with read alignments in genomic coordinates into transcriptomic coordinates using a given transcript annotation file.'
    )
    parser.add_argument('-i', '--input-bamfile', required=True,
                        help='Input bam file')
    parser.add_argument('-g', '--input-gtffile', required=True,
                        help='Input gtf file.')
    parser.add_argument('-o', '--output-filename', required=True,
                        help='Output filename for bamfile.')
    
    args = parser.parse_args()
    if os.path.isfile(args.input_bamfile):
            bam_file = ul.open_bam(args.input_bamfile)
    else:
            print(f"Error: {args.input_bamfile} is not a valid file.")
            sys.exit(1)
    if os.path.isfile(args.input_gtffile):
            print("Processing gtf file...")
            gtf_file = ul.load_gtf(args.input_gtffile)
    else:
            print(f"Error: {args.input_gtffile} is not a valid file.")
            sys.exit(1)
    print("processing reads...")
    ul.process_reads(bam_file,args.output_filename,gtf_file)
    bam_file.close()
    
if __name__ == "__main__":
    main()
