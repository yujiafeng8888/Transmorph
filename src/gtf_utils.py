
import pandas as pd
import sys
import pysam

# 解析 GTF 文件
def load_gtf(gtf_file):
    gtf = pd.read_csv(gtf_file,sep="\t",comment="#", header=None,
                      names=["seqname","source","feature",
                      "start","end","score","strand","frame","attribute"])
    gtf["transcript_id"] = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"')
    return gtf

# 提取转录本和外显子信息
def extract_transcripts_and_exons(gtf):
    transcripts=gtf[gtf["feature"]=='transcript']
    exons = gtf[gtf["feature"] == 'exon']
    return transcripts,exons

# 打开 BAM 文件
def open_bam(bam_file):
    return pysam.AlignmentFile(bam_file, "rb")

def build_transcript_header(bam, transcripts,exons):
    header = bam.header.to_dict()
    # 清空 SQ 部分
    header["SQ"] = []
    transcripts = transcripts.drop_duplicates("transcript_id")
    for _, row in transcripts.iterrows():
        transcript_id = row["transcript_id"]
        exon_rows = exons[exons["transcript_id"] == transcript_id]
        transcript_len = sum(exon_rows["end"] - exon_rows["start"] + 1)
        header["SQ"].append({"LN": int(transcript_len),"SN": str(transcript_id) })
    return header

# 读取每个 read
def process_reads(bam_file,output_bamname,gtf):
    transcripts,exons = extract_transcripts_and_exons(gtf)
    header_out = build_transcript_header(bam_file, transcripts,exons)
    # output_bam_name = output_bamname + ".bam"
    out_bam= pysam.AlignmentFile(output_bamname, "wb", header=header_out)
    seen = set()
    for read in bam_file.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        if read.query_name in seen:
            continue  # 跳过已经处理过的 read
        seen.add(read.query_name)
        genome_to_transcript_coords(bam_file, read, exons, transcripts, out_bam)

    for read in bam_file:
        if read.is_unmapped:
            continue  # 跳过未比对的 reads
        # print(out_bam.get_reference_name(0))
        # read.reference_id=875
        genome_to_transcript_coords(bam_file,read, exons,transcripts,out_bam)
    out_bam.close()
    print("process finished!")
    

def genome_to_transcript_coords(bam_file,read, exons, transcripts,out_bam):
    """
    将基因组坐标转换为转录本坐标，返回转录本坐标范围（如果在外显子内）。
    如果在内含子或不在外显子内，返回 None。
    """
    try:
        if not read or not hasattr(read, 'reference_start'):
            return None
        genome_start = read.reference_start + 1  # 转换为1-based坐标
        genome_end = read.reference_end
        reference_id = read.reference_id
        chrom = bam_file.get_reference_name(reference_id)
        chrom_list=list(exons["seqname"])
        if chrom not in chrom_list:
            print(chrom)
            print("chr name in bam do not match with gtf file, skip")
            return None
            # sys.exit(1)
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
            trans_end=transcript["end"].values[0]
            transcript_len=trans_end-trans_start+1
            exon_start = row['start']
            exon_end = row['end']
            # 检查 read 是否完全位于外显子内
            if genome_start >= exon_start and genome_end <= exon_end:
                # 计算相对转录本的坐标
                # print("find a read can be process")
                coords_start = genome_start - trans_start
                # 创建新的 read 对象或者直接修改 read 的坐标
                # 这里直接修改 read 对象，也可以返回新对象
                read.reference_start = coords_start-1
                # read.reference_end = coords_end
                if coords_start - 1 >= transcript_len:
                    print(f"Warning: {trans_id} coords_start {coords_start} > length {transcript_len}")
                    continue
                trans_id = str(trans_id).strip()
                if read.is_paired:
                    if read.is_read1:  # 如果是 read1，更新 read2
                        mate_read = bam_file.mate(read)
                        # 更新配对读取的坐标
                        mate_read.reference_start = coords_start - 1
                        out_bam.write(mate_read)
                out_bam.write(read)

    except ValueError as e:
        print(f"Skipping read {trans_id} name: {read.query_name} due to error: {e}")

