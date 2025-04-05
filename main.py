# main.py
import pysam
from gtf_utils import load_transcripts, genome_to_transcript_coords

# 加载注释文件
transcripts = load_transcripts("annotation.gtf")

# 打开输入输出 BAM 文件
in_bam = pysam.AlignmentFile("input_genomic.bam", "rb")
out_bam = pysam.AlignmentFile("output_transcriptomic.bam", "wb", template=in_bam)

for read in in_bam:
    if read.is_unmapped:
        continue

    # 获取所有可能匹配的 transcript 坐标
    tx_alignments = genome_to_transcript_coords(read, transcripts)

    if not tx_alignments:
        continue  # 没有匹配的 transcript，过滤掉

    for tx_read in tx_alignments:
        out_bam.write(tx_read)

in_bam.close()
out_bam.close()
