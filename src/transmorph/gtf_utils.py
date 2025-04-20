
import pandas as pd
import sys
import pysam

# read GTF file
def load_gtf(gtf_file):
    gtf = pd.read_csv(gtf_file,sep="\t",comment="#", header=None,
                      names=["seqname","source","feature",
                      "start","end","score","strand","frame","attribute"])
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"')
    return gtf

# extract transcripts and exons
def extract_transcripts_and_exons(gtf):
    transcripts=gtf[gtf["feature"]=='transcript']
    exons = gtf[gtf["feature"] == 'exon']
    return transcripts,exons

# open bam file
def open_bam(bam_file):
    return pysam.AlignmentFile(bam_file, "rb")

def build_transcript_header(bam, transcripts, exons):
    header = bam.header.to_dict()
    header["SQ"] = []
    transcripts = transcripts.drop_duplicates("transcript_id")
    # 用 groupby 先组织好
    exon_groups = exons.groupby("transcript_id")
    for _, row in transcripts.iterrows():
        transcript_id = row["transcript_id"]
        if transcript_id not in exon_groups.groups:
            continue  # this transcript doesn't have exon，skip
        exon_rows = exon_groups.get_group(transcript_id)
        transcript_len = (exon_rows["end"] - exon_rows["start"] + 1).sum()
        header["SQ"].append({
            "LN": int(transcript_len),
            "SN": str(transcript_id)
        })
    print("return head")
    return header

# process every reads
def process_reads(bam_file,output_bamname,gtf):
    transcripts,exons = extract_transcripts_and_exons(gtf)
    header_out = build_transcript_header(bam_file, transcripts,exons)
    out_bam= pysam.AlignmentFile(output_bamname, "wb", header=header_out)

    for read in bam_file:
        if read.is_unmapped:
            continue  # skip unmapped reads
        genome_to_transcript_coords(bam_file,read, exons,transcripts,out_bam)
    out_bam.close()
    print("process finished!")
    
def genome_to_transcript_coords(bam_file,read, exons, transcripts,out_bam):
    """
    change genome coordinates to transcript coordinates, if reads mapped to intron, skip this read
    """
    try:
        if not read or not hasattr(read, 'reference_start'):
            return None
        genome_start = read.reference_start + 1  # change to 1-based 
        genome_end = read.reference_end
        reference_id = read.reference_id
        chrom = bam_file.get_reference_name(reference_id)
        chrom_list=list(exons["seqname"])
        if chrom not in chrom_list:
            print(chrom)
            print("chr name in bam do not match with gtf file, skip")
            return None
        else:
            sub_transcript = transcripts[(transcripts["seqname"]==chrom)&(transcripts["start"] <= genome_start) & (transcripts["end"] >= genome_end)]
            sub_transcript_id=list(sub_transcript["transcript_id"])
        if sub_transcript_id is not None:
            exon_chrom = exons[(exons["transcript_id"].isin(sub_transcript_id)) & (exons["start"]<=genome_start)]
        else:
            return None
        print(exon_chrom.shape[0])
        for _, row in exon_chrom.iterrows():
            trans_id = row["transcript_id"]
           #find transcript start position
            transcript = transcripts[transcripts["transcript_id"] == trans_id]
            if transcript.empty:
                continue
            trans_start = transcript["start"].values[0]  # transcripts start position
            trans_end=transcript["end"].values[0]
            transcript_len=trans_end-trans_start+1
            exon_start = row['start']
            exon_end = row['end']
            # check if read mapped to exons 
            if genome_start >= exon_start and genome_end <= exon_end:
                coords_start = genome_start - trans_start
                read.reference_start = coords_start-1
                if coords_start - 1 >= transcript_len:
                    print(f"Warning: {trans_id} coords_start {coords_start} > length {transcript_len}")
                    continue
                trans_id = str(trans_id).strip()
                print(f"write {trans_id}")
                out_bam.write(read)

    except ValueError as e:
        print(f"Skipping read {trans_id} name: {read.query_name} due to error: {e}")

