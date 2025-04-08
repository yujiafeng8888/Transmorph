
import pandas as pd
import re
import sys
# 解析 GTF 文件
def load_gtf(gtf_file):
    gtf = pd.read_csv(gtf_file,sep="\t",comment="#", header=None,
                      names=["seqname","source","feature",
                      "start","end","score","strand","frame","attribute"])
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"')
    # gtf["transcript_id"]=gtf["attribute"].apply(extract_transcript_id)
    # gtf["gene_id"]=gtf["gene_id"].apply(extract_gene_id)
    gtf["gene_id"] = gtf["attribute"].str.extract(r'gene_id "([^"]+)"')
    return gtf

# def extract_transcript_id(attr):
#     match = re.search(r'\btranscript_id\s+"([^"]+)"', attr)
#     if match:
#         return match.group(1)
#     else:
#         print(f"do not have transcript_id{attr}")
#         return None
    
# def extract_gene_id(attr):
#     match=re.search(r'gene_id "([^"]+)"',attr)
#     return match.group(1)

# 提取转录本和外显子信息
def extract_transcripts_and_exons(gtf):
    # 提取转录本
    # transcripts = gtf[gtf["feature"] == 'transcript']
    # 提取外显子
    transcripts=gtf[gtf["feature"]=='transcript']
    exons = gtf[gtf["feature"] == 'exon']
    return transcripts,exons

import pysam

# 打开 BAM 文件
def open_bam(bam_file):
    return pysam.AlignmentFile(bam_file, "rb")

def build_transcript_header(bam, transcripts):
    header = bam.header.to_dict()
    # 清空原有 SQ，或者保留原始参考
    header["SQ"] = []

    for _, row in transcripts.iterrows():
        transcript_id = row["transcript_id"]
        transcript_len = row["end"] - row["start"] + 1
        header["SQ"].append({"SN": transcript_id, "LN": transcript_len})

    return header

# 读取每个 read
def process_reads(bam_file, out_bam,gtf):
    transcripts,exons = extract_transcripts_and_exons(gtf)
    header = build_transcript_header(bam_file, transcripts)
    output_bam_name = out_bam + ".bam"
    output_bam = pysam.AlignmentFile(output_bam_name, "wb", header=header)
    # for sq in output_bam.header["SQ"]:
    #     print(sq["SN"])
    for read in bam_file:
        if read.is_unmapped:
            continue  # 跳过未比对的 reads
        
        tx_coords = genome_to_transcript_coords(bam_file,read, exons,transcripts,output_bam)
        
        if tx_coords !=None:
            output_bam.write(tx_coords)
    output_bam.close()

def genome_to_transcript_coords(bam_file,read, exons, transcripts,output_bam):
    """
    将基因组坐标转换为转录本坐标，返回转录本坐标范围（如果在外显子内）。
    如果在内含子或不在外显子内，返回 None。
    """
    genome_start = read.reference_start + 1  # 转换为1-based坐标
    genome_end = read.reference_end
    reference_id = read.reference_id
    chrom = bam_file.get_reference_name(reference_id)
    # print(chrom)
    # print("start:",genome_start)
    # print("end:",genome_end)
    chrom_list=list(exons["seqname"])
    if chrom not in chrom_list:
        print("chr name in bam do not match with gtf file")
        sys.exit(1)
    else:
        sub_transcript = transcripts[(transcripts["seqname"]==chrom)&(transcripts["start"] <= genome_start) & (transcripts["end"] >= genome_end)]
        sub_transcript_id=list(sub_transcript["transcript_id"])
    # 提取与当前染色体匹配的外显子
    # print(sub_transcript_id)
    if sub_transcript_id is not None:
        exon_chrom = exons[(exons["transcript_id"].isin(sub_transcript_id)) & (exons["start"]<=genome_start)]
    # print(sub_transcript)
    else:
        return None
    # 遍历每个外显子
    for _, row in exon_chrom.iterrows():
        trans_id = row["transcript_id"]
        # 找到相应的转录本起始位置
        transcript = transcripts[transcripts["transcript_id"] == trans_id]
        # print(transcript)
        if transcript.empty:
            continue
        trans_start = transcript["start"].values[0]  # 转录本起始位置
        exon_start = row['start']
        exon_end = row['end']
        # 检查 read 是否完全位于外显子内
        if genome_start >= exon_start and genome_end <= exon_end:
            # 计算相对转录本的坐标
            print("find a read can be process")
            coords_start = genome_start - trans_start
            coords_end = genome_end - trans_start

            # 创建新的 read 对象或者直接修改 read 的坐标
            # 这里直接修改 read 对象，也可以返回新对象
            read.reference_start = coords_start
            # read.reference_end = coords_end
            read.reference_id = output_bam.get_tid(trans_id)  # 更新为转录本的 ID

            return read

    # 如果没有在外显子内，返回 None
        return None



# def cal_end(read):
#     cigar=read.cigarstring
#     start_pos=read.reference_start
#     aligned_length=0
#     for length,operation in read.cigartuples:
#         if operation in (0,2,7,8):
#             aligned_length+=length
#     end_pos=start_pos + aligned_length-1
#     return end_pos
