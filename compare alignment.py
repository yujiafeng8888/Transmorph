import pandas as pd
import pysam

gtf_path = "/Users/yujiafeng/Desktop/project/chess3.0.1.primary.gtf"
isoform_list_output = "isoform_list.txt"
isoform_ids = set()

with open(gtf_path) as f:
    for line in f:
        if line.startswith("#"):
            continue  
        fields = line.strip().split("\t")
        if fields[2] == "transcript":  
            attributes = fields[8]
            attrs = {item.strip().split(' ')[0]: item.strip().split(' ')[1].replace('"', '') 
                     for item in attributes.strip(';').split(';') if item.strip()}
            transcript_id = attrs.get('transcript_id')
            if transcript_id:
                isoform_ids.add(transcript_id)

with open(isoform_list_output, "w") as f:
    for tid in sorted(isoform_ids):
        f.write(tid + "\n")

print(f"{len(isoform_ids)} isoforms are saved to {isoform_list_output}")



bam1_path = "/Users/yujiafeng/Desktop/project/bowtie.bam"
bam2_path = "transmorph.bam"
isoform_list_path = "isoform_list.txt"
output_path = "isoform_read_counts.csv"

# read isoform_list
with open(isoform_list_path) as f:
    isoform_list = [line.strip() for line in f if line.strip()]


read_counts = pd.DataFrame({
    'isoform_id': isoform_list,
    'bam1_count': [0] * len(isoform_list),
    'bam2_count': [0] * len(isoform_list)
}).set_index('isoform_id')

def count_reads(bam_path, column_name):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    for read in bamfile.fetch(until_eof=True):
        if not read.is_unmapped:  
            isoform = read.reference_name
            if isoform in read_counts.index:
                read_counts.at[isoform, column_name] += 1
    bamfile.close()

count_reads(bam1_path, "bam1_count")
count_reads(bam2_path, "bam2_count")
read_counts = read_counts[~((read_counts['bam1_count'] == 0) & (read_counts['bam2_count'] == 0))]
read_counts.to_csv(output_path)

print(f"finished! output path: {output_path}")
