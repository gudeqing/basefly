import os
import sys
import re
from collections import Counter
import pandas as pd
from pysam import VariantFile, VariantHeader
from pysam import FastaFile
from pysam import AlignmentFile
from Bio import Align
import logging

__author__ = 'gdq'

""""
pip install -i https://mirrors.aliyun.com/pypi/simple/ biopython pysam
测试路径：/home/report/gdqtest
dragen镜像：/home/novaseq/singularity_images/trusight-oncology-500-ruo.img 
/home/Novaseq_Data/dragen/230604_A01765_0080_BHMNLGDSX5_tissue/Logs_Intermediates/VariantCaller/DCA0140456801T/DCA0140456801T.filtered.genome.vcf
/staging/analysis/230604_A01765_0080_BHMNLGDSX5_tissue_from_fastq
/mnt/TSO500/dragen/230604_A01765_0080_BHMNLGDSX5_tissue/Logs_Intermediates/VariantCaller/DCA0140456801T/
------------------------------------------------------------------------------------------------------------------------
突变表示形式一（非期望的表达形式，因在vcf中需要用【1个del突变】和【1个SNP突变】才能一起完整表示下面的突变）：
Ref: G GAAxxxTTAAGA GAAxxx G
Alt: G ------------ GAAxxx C  ---> GGAAxxxC

突变表示形式二（期望的表达形式，因在vcf中可以用一个复合突变表示）：
Ref: G GAAxxxTTAAGA GAAxxx G
Alt: G GAAxxx------ ------ C  ---> GGAAxxxC，del+snp的混合
或者：
Alt: G GAAxxxC----- ------ -  ---> GGAAxxxC，snp + del的混合

1. 可以发现，本质上相同的突变可以至少用2种形式表示，但表面看却有比较大的区别
2. 对于第一种形式，突变软件可能会报出2个突变,且2个突变的距离长度为len(GAAxxx),可见x越多，距离越远
3. 对于第二种形式，则只需要报一个复合突变

同理：上述突变反转看待，软件可能报出一个插入，同时间隔一个突变，我们也期望可以直接变成连续插入和突变，使用complex表示
Ref1：G ------------ GAAxxx C
Ref2：G GAAxxxC----- ------ -  (本质和Ref1序列相同）
ALT： G GAAxxxTTAAGA GAAxxx G  这组突变相对Ref1，需要标为2个不连续的突变，而相对Ref2则可以表示为连续的突变。

------------------------------------------------------------------------------------------------------------------------
问题：如何将第一种形式的突变转化为第二种最精简的表示形式
解答思路：
a. 最简单的情况，如果2个突变都是SNP或等长替换，且两个坐标差<=2则直接合并。
b. 一个deletion，一个为snp，两个坐标相差<=50, 则采用如下策略合并。
    1. 将vcf中2个可能需要合并的突变合并成一条突变，获得突变后的序列,命名为FinalSeq，如上述“GGAAxxxC"
    2. 将Finalseq和参考序列如上述”G GAAxxxTTAAGA GAAxxx G“进行比对
    3. 比对方法选择”全局比对“，使用NCBI全局比对工具得到的结果如下：
    -------------------------------------------
    NW Score	Identities	Gaps	Strand
    -------------------------------------------
    -18	4/14(29%)	9/14(64%)	Plus/Plus
    -------------------------------------------
    Query  1   GGAATTAAGAGAAG  14
               ||||          
    Sbjct  1   GGAAC           5
    --------------------------------------------
    4. 采用biopython包中的PairwiseAligner比对替代NCBI网页工具。
    如果algined部分只有一组索引，即去除match的部分，剩下的都是gap，则可以成功的合并2个突变为一条。
c. 一个insertion和snp的情况同上也可以处理
d. 合并阈值说明：
(1）对于snp，要求同时支持2个snp的read(read1和read2只计数一次）数量大于等于3, 3是凭经验选定
(2）对于50bp距离以内的2个突变（可以是snp+indel, indel+indel)，符合以下任一情况，建议合并：
    情况1：同时支持2个突变read数量大于等于3，而且合并后的突变的支持reads数量大于等于2，3和2都是凭经验选定，可根据需要调整
    情况2：仅仅要求支持合并后的突变的read数量大于等于20，20是根据经验选定
    说明：测试结果发现，通过pysam找到的支持合并后的突变的read数量往往比较少，推测是因为比对算法本身导致的，后续可以考虑针对这个进行优化，
    比如把可能的read序列和突变单独重新进行局部比对，看是否支持突变。

# 下面是一个典型的处于重复区域或微卫星区域的情况，也可以被识别
target="GGGAGGAGGAGGAGGAGGAGGAGGATGAG"
query ="GGGAGGAGGAGGAGGAGGAGGAT"

对于距离超过2且突变类型都是snp或等长替换的，均不考虑合并突变

## 新的问题，在微卫星序列附近的突变
chr21	42864143	.	AAC	A
chr21	42864145	.	C	CA
chr21	42864145	.	C	CAA
chr21	42864145	.	CA	C
上述突变实际是一种MSI现象，可见同一个位点可以有多种突变形式，当前程序不考虑处理此种特殊情况。
"""


def init_vcf(vcf_path, genome='hg19', chrom_name_is_numeric=False):
    vcf = VariantFile(vcf_path, 'w', header=VariantHeader())
    contig_info = [
        '##fileformat=VCFv4.2',
        f'##assembly={genome}',
        "##contig=<ID=chr1,length=249250621>",
        "##contig=<ID=chr2,length=243199373>",
        "##contig=<ID=chr3,length=198022430>",
        "##contig=<ID=chr4,length=191154276>",
        "##contig=<ID=chr5,length=180915260>",
        "##contig=<ID=chr6,length=171115067>",
        "##contig=<ID=chr7,length=159138663>",
        "##contig=<ID=chr8,length=146364022>",
        "##contig=<ID=chr9,length=141213431>",
        "##contig=<ID=chr10,length=135534747>",
        "##contig=<ID=chr11,length=135006516>",
        "##contig=<ID=chr12,length=133851895>",
        "##contig=<ID=chr13,length=115169878>",
        "##contig=<ID=chr14,length=107349540>",
        "##contig=<ID=chr15,length=102531392>",
        "##contig=<ID=chr16,length=90354753>",
        "##contig=<ID=chr17,length=81195210>",
        "##contig=<ID=chr18,length=78077248>",
        "##contig=<ID=chr19,length=59128983>",
        "##contig=<ID=chr20,length=63025520>",
        "##contig=<ID=chr21,length=48129895>",
        "##contig=<ID=chr22,length=51304566>",
        "##contig=<ID=chrMT,length=16569>",
        "##contig=<ID=chrX,length=155270560>",
        "##contig=<ID=chrY,length=59373566>",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=Gene,Number=1,Type=String,Description="Gene Symbol">',
        '##INFO=<ID=pHgvs,Number=1,Type=String,Description="P-Dot Notation">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">',
    ]
    if chrom_name_is_numeric:
        contig_info = [x.replace('=chr', '=') for x in contig_info]
    for line in contig_info:
        vcf.header.add_line(line)
    return vcf


def read_tso500_small_variant_tsv(infile):
    discarded_lines = []
    with open(infile) as fr:
        while True:
            line = next(fr)
            discarded_lines.append(line)
            if line.startswith("[Small Variants]"):
                break
        table = pd.read_csv(fr, header=0, skip_blank_lines=True, sep='\t', dtype=str)
    # skip_blank_lines不起作用，用dropna处理最后一个空行
    table = table.dropna(how='all')
    return table, discarded_lines


def small_variant_tsv_to_vcf(infile, out=None, sample_name=None):
    table, _ = read_tso500_small_variant_tsv(infile)
    out = out or infile[:-3] + 'SmallVariant.vcf'
    vcf = init_vcf(out, genome='hg19', chrom_name_is_numeric=False)
    sample_name = sample_name or os.path.basename(infile).split('_Combine')[0]
    vcf.header.add_sample(sample_name)
    for idx, row in table.iterrows():
        record = vcf.new_record()
        record.contig = row['Chromosome']
        record.pos = int(row['Genomic Position'])
        record.ref = row['Reference Call']
        record.alts = [row['Alternative Call']]
        record.qual = None
        record.filter.add('PASS')
        record.info.update(dict(
            Gene=str(row['Gene']),
            DP=int(row['Depth']),
            AF=float(row['Allele Frequency']),
            pHgvs=str(row['P-Dot Notation'])
        ))
        record.samples[sample_name]['DP'] = int(row['Depth'])
        record.samples[sample_name]['AF'] = float(row['Allele Frequency'])
        vcf.write(record)
    else:
        vcf.close()
    return out


def vep_annotation(infile, fasta=None, cache_dir=None, outfile=None, threads=4):
    cache_dir = cache_dir or "/home/report/ecloud/nas_002/biodb/vep/grch37"
    fasta = fasta or "/home/report/ecloud/nas_002/genome_and_annotation/hg19/hg19.fa"
    outfile = outfile or infile[:-3] + 'vep.vcf'
    mount_vols = [
        os.path.dirname(cache_dir),
        os.path.dirname(os.path.abspath(fasta)),
        os.path.dirname(os.path.abspath(infile)),
        os.path.dirname(os.path.abspath(infile))
    ]
    cmd = 'singularity exec -B {} '.format(','.join(set(mount_vols)))
    cmd += "/home/novaseq/singularity_images/ensembl-vep:104.3--pl5262h4a94de4_0 "
    cmd += "vep "
    cmd += f'-i {infile} '
    cmd += f'--fasta {fasta} '
    cmd += f'-o {outfile} '
    cmd += f'--dir_cache {cache_dir} '
    cmd += f'--fork {threads} '
    cmd += f'--vcf '
    cmd += f'--hgvs '
    cmd += f'--use_given_ref '
    cmd += f'--stats_file {outfile}.summary.html '
    # cmd += f'--compress_output bgzip '
    cmd += f'--force_overwrite '
    cmd += f'--offline '
    cmd += f'--merged '
    cmd += f'--variant_class '
    cmd += f'--terms SO '
    cmd += f'--sift b '
    cmd += f'--polyphen b '
    cmd += f'--symbol '
    cmd += f'--tsl '
    cmd += f'--no_escape '
    cmd += f'--canonical '
    cmd += f'--flag_pick '
    print(cmd)
    os.system(cmd)
    return outfile


def add_vep_annot_back_to_tsv(vcf_file, tsv_file, out=None):
    table, discarded_lines = read_tso500_small_variant_tsv(tsv_file)
    annot_dict = dict()
    with VariantFile(vcf_file) as f:
        # get csq format
        csq_format = f.header.info['CSQ'].description.split('Format: ')[1]
        for r in f:
            picked = None
            canonical = None
            for each in r.info['CSQ']:
                csq_dict = dict(zip(csq_format.split('|'), each.split('|')))
                # 找到被flag为1的记录作为报告
                if csq_dict['PICK'] == "1":
                    picked = csq_dict
                # 找到canonical转录本，且
                if csq_dict['CANONICAL'] == 'YES' and csq_dict['SOURCE'] == 'RefSeq':
                    canonical = csq_dict
                if picked and canonical:
                    break
            # print(picked, canonical)
            csq_dict = canonical or picked
            target_info = [
                csq_dict['SYMBOL'],  # Gene
                r.contig,  # Chromosome
                r.pos,  # Genomic Position
                r.ref,  # Reference Call
                r.alts[0],  # Alternative Call
                r.info['AF'],  # Allele Frequency
                r.info['DP'],  # Depth
                csq_dict['HGVSp'].replace('%3D', '='),  # P-Dot Notation
                csq_dict['HGVSc'].replace('%3D', '='),  # C-Dot Notation
                csq_dict['Consequence'],  # Consequence(s)
                csq_dict['EXON']  # Affected Exon(s)
            ]
            annot_dict[tuple(r.info['MergeFrom'].split('|'))] = target_info
    new = pd.DataFrame(annot_dict.values(), columns=table.columns).round(4)
    result = pd.concat([table, new])
    # result = result.infer_objects()
    result['Genomic Position'] = result['Genomic Position'].astype(int)
    result['Depth'] = result['Depth'].astype(int)
    out = out or os.path.basename(tsv_file)[:-3] + 'updated.tsv'
    with open(out, 'w') as fw:
        for line in discarded_lines:
            fw.write(line)
        result.to_csv(fw, sep='\t', index=False)
    return out


def detect_snp(bam, contig, start, min_bq=0, ignore_overlaps=False, fasta_file=None):
    cols = bam.pileup(
        contig, start, start + 1,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
        # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
        max_depth=300000,
        fastafile=fasta_file,
        # flag_filter=4,
    )
    # ref = fasta_file.fetch(contig, start, start+1)
    base_dict = dict()
    for col in cols:
        for base, read in zip(
                # 如不加add_indels参数，那么将无法知晓插入的碱基序列
                col.get_query_sequences(),
                col.get_query_names()
        ):
            base_dict.setdefault(base, 0)
            base_dict[base] += 1
    depth = sum(base_dict.values())
    sorted_base_freq = sorted(zip(base_dict.keys(), base_dict.values()), key=lambda x: x[1])
    alt, freq = sorted_base_freq[-1]
    return alt, freq, depth


def get_snp_support_reads(bam, contig, start, alt, min_bq=0, ignore_overlaps=False, fasta_file=None):
    cols = bam.pileup(
        contig, start, start+1,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
        # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
        max_depth=300000,
        fastafile=fasta_file,
        # flag_filter=4,
    )
    support_reads = set()
    for col in cols:
        for base, read in zip(
                # 如不加add_indels参数，那么将无法知晓插入的碱基序列
                col.get_query_sequences(),
                col.get_query_names()
        ):
            if base == alt:
                support_reads.add(read)
    return support_reads


def get_insert_support_reads(bam, contig, start, alt, min_bq=0, ignore_overlaps=False, fasta_file=None):
    cols = bam.pileup(
        contig, start, start+1,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
        # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
        max_depth=300000,
        fastafile=fasta_file,
        # flag_require=read_type,
        # flag_filter=4,
    )
    support_reads = set()
    expected_insert = alt.upper()
    for col in cols:
        for base, read in zip(
                # 如不加add_indels参数，那么将无法知晓插入的碱基序列
                col.get_query_sequences(add_indels=True),
                col.get_query_names()
        ):
            if '+' in base:
                # such as 'G+2AT', +号前面是参考序列
                insertion = re.sub(r'\+\d+', '', base).upper()
                # 暂时几乎不考虑insertion恰巧出现在read的两端且不包含完整的insertion
                mismatch_allowed = round(len(expected_insert)*0.15)
                mismatched_num = sum(x != y for x, y in zip(expected_insert, insertion))
                if mismatched_num <= mismatch_allowed:
                    support_reads.add(read)
    return support_reads


def get_del_support_reads(bam, contig, start, del_len, min_bq=0, ignore_overlaps=False, fasta_file=None):
    cols = bam.pileup(
        contig, start, start+1,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
        # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
        max_depth=300000,
        fastafile=fasta_file,
        # flag_filter=4,
    )
    support_reads = set()
    for col in cols:
        for pileup_read in col.pileups:
            if -pileup_read.indel == del_len:
                support_reads.add(pileup_read.alignment.query_name)
    return support_reads


def get_substitution_support_reads(bam, contig, start, alt, min_bq=0, ignore_overlaps=False, fasta_file=None):
    # 该函数可以被get_complex_support_reads替代
    cols = bam.pileup(
        contig, start, start+len(alt),
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
        # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
        max_depth=300000,
        fastafile=fasta_file,
        # flag_filter=4,
    )
    support_reads = list()
    for idx, col in enumerate(cols):
        for base, qual, read in zip(
                # 如不加add_indels参数，那么将无法知晓插入的碱基序列
                col.get_query_sequences(add_indels=True),
                col.get_query_qualities(),
                col.get_query_names()
        ):
            if base == alt[idx]:
                support_reads.append(read)
    # 如果一个read匹配到的次数等于alt的长度，则认为该read支持alt
    count = Counter(support_reads)
    alt_length = len(alt)
    support_reads = set(x for x, y in count.items() if y == alt_length)
    return support_reads


def get_complex_support_reads(bam, mut, min_bq=0, ignore_overlaps=False, fasta_file=None):
    """
    思路：提取vcf中突变起始位置的前一个碱基对应的比对信息，分析每条比对到该位置的read序列是否支持alt
    当突变是AAA->A这种del形式时，上面的思路不可行, 因为没有突变的read也将符合程序的判定，如下突变就是个典型案例
    12      49444957        .       ACTGGGGGGACAGGTGTGATTCCTCAGGTTGGGGGGACAAGCATGGCTCCTCAGGCACAGGAGACAGGTGCGGCTCCTCAGT      A
    当突变时ATC->G或ATC->GC这种complex形式时，上述思路可行。
    注意：上述思路并非严谨的寻找支持reads，还是以"ATC->G或ATC->GC"为例，容易想到支持”ATC->GC”或“AT->G"的read同样也支持"ATC->G"
    为增加严谨性：
    1. 假设alt后的3个碱基一定是和reference匹配的
    2. 对于complex突变，可以理解为是：删除ref，加入alt, 假设ref是ATC，alt是GT，ref前后是X，那么突变的read应该是XGTX
    3. 取X的长度为3
    """
    contig, start = mut.contig, mut.start
    ref, alt = mut.ref, mut.alts[0]
    # mut_type = get_mutation_type(mut)
    cols = bam.pileup(
        # 提取ref的前一个碱基对应的比对信息
        contig, start-1, start,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
        # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
        max_depth=300000,
        fastafile=fasta_file,
    )
    support_reads = set()
    # 提取删除ref后的3个参考碱基, 然后将alt序列与之拼接，得到的应该是突变后的真实read信息
    # 注意如果遇到串联重复序列，该方法将无效
    expected_seq_after = fasta_file.fetch(contig, start + len(ref), start + len(ref) + 3)
    expected_seq_before = fasta_file.fetch(contig, start - 3, start)
    #
    for col in cols:
        for pileup_read in col.pileups:
            query_seq = pileup_read.alignment.query_sequence
            query_pos = pileup_read.query_position
            if query_pos is not None:
                # 我们期望该位置下一位|mismatch|insertion|deletion|
                # 提取真实read信息： 期望是（3个和参考一致的碱基+可能的插入序列+3个和参考一致的碱基） | 6个和参考一致的碱基
                if query_pos >= 2:
                    back_extend = 3
                elif query_pos == 1:
                    back_extend = 2
                else:
                    back_extend = 1

                if query_pos + len(alt) + 4 <= len(query_seq) - 3:
                    forward_extend = 3
                elif query_pos + len(alt) + 4 == len(query_seq) - 2:
                    forward_extend = 2
                elif query_pos + len(alt) + 4 == len(query_seq) - 1:
                    forward_extend = 1
                else:
                    forward_extend = 0

                real_seq = query_seq[query_pos + 1 - back_extend:query_pos + 1 + len(alt) + 3]
                # print("xxx", contig, start, real_seq, alt, expected_seq)
                if real_seq == expected_seq_before.upper()[-back_extend:] + alt + expected_seq_after.upper()[:forward_extend]:
                    support_reads.add(pileup_read.alignment.query_name)
                    # print(query_seq, query_pos, expected_seq_before, ref+'>'+alt, expected_seq_after)
    return support_reads


def get_support_reads(mut, mut_type, bam, fasta):
    # 提取同时支持2个突变的read信息
    supports = None
    if mut_type == 'snp':
        supports = get_snp_support_reads(bam, mut.contig, mut.start, mut.alts[0], fasta_file=fasta)
    elif mut_type == 'substitution':
        supports = get_substitution_support_reads(bam, mut.contig, mut.start, mut.alts[0], fasta_file=fasta)
        # supports = get_complex_support_reads(bam, mut, fasta_file=fasta)
    elif mut_type == 'deletion':
        del_size = len(mut.ref) - len(mut.alts[0])
        supports = get_del_support_reads(bam, mut.contig, mut.start, del_size, fasta_file=fasta)
    elif mut_type == 'insertion':
        supports = get_insert_support_reads(bam, mut.contig, mut.start, mut.alts[0], fasta_file=fasta)
    elif mut_type == 'complex':
        supports = get_complex_support_reads(bam, mut, fasta_file=fasta)
    else:
        raise Exception(f"{mut_type} is unexpected")
    # 去掉区分read1和read2的字符串
    supports = set(x[:-2] for x in supports)
    return supports


def global_pair_seq_align(target_pos=55242465, target='GGAATTAAGAGAAG', query='GGAAC'):
    """
    EGFR突变的例子：
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    7	55242465	.	GGAATTAAGA	G	.	.	.
    7	55242478	.	G	C	.	.	.
    上述等同下面的突变
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    7	55242469	.	TTAAGAGAAG	C	.	.	.
    ---------------------------------
    Score = 7.5:
    target            0 GGAA-TTAAGAGAAG 14
                      0 ||||----------- 15
    query             0 GGAAC----------  5
    """
    aligner = Align.PairwiseAligner()
    # aligner.mode = 'global'
    # aligner.match_score = 2  # 2
    # aligner.mismatch_score = -1  # -3
    # aligner.open_gap_score = -0.5  # -5
    # aligner.extend_gap_score = -0.1  # -2
    # aligner.target_end_gap_score = 0.0
    # aligner.query_end_gap_score = 0.0
    alignments = aligner.align(target.upper(), query.upper())
    # for alignment in alignments:
    #     print("Score = %.1f:" % alignment.score)
    #     print(alignment)
    # 仅仅取第一个比对结果进行判断
    alignment = alignments[0]
    # print(alignment)
    target_aligned_idx, query_aligned_idx = alignment.aligned
    map_start_from_0 = (target_aligned_idx[0][0] == 0) and (query_aligned_idx[0][0] == 0)
    aligned_section_number = len(target_aligned_idx)
    second_align_is_1bp = False
    if aligned_section_number == 2:
        target_align_2_len = target_aligned_idx[1][1] - target_aligned_idx[1][0]
        query_align_2_len = query_aligned_idx[1][1] - query_aligned_idx[1][0]
        second_align_is_1bp = target_align_2_len + query_align_2_len <= 2

    if (aligned_section_number == 1 and map_start_from_0) or (aligned_section_number == 2 and map_start_from_0 and second_align_is_1bp):
        # alt和ref对比，从最左边开始计算，有且仅有一段完全匹配的区域
        ref = alignment.target[target_aligned_idx[0][1]:]
        alt = alignment.query[query_aligned_idx[0][1]:]
        new_pos = target_pos + target_aligned_idx[0][1]
        # print(target_aligned_idx)
        # print(query_aligned_idx)
        if alt == "" or ref == "":
            # print(alignment)
            # print("说明两个突变可以合并为一个缺失或插入的突变，此时需修改ref和alt")
            ref = alignment.target[target_aligned_idx[0][1] - 1:]
            alt = alignment.query[query_aligned_idx[0][1] - 1:]
            new_pos = target_pos + target_aligned_idx[0][1] - 1
        return new_pos, ref, alt, alignment
    else:
        # print(alignment)
        # print('比对结果存在超过2段以上的匹配，因此不建议合并突变')
        return None, None, None, alignment


def distance(mut_a, mut_b):
    """
    必须按顺序输入，要求mut_a的坐标小于mut_b的坐标
    """
    # 跳过alt为”."的记录
    if (mut_a.alts is None) or (mut_b.alts is None):
        return 10000
    if mut_a.contig != mut_b.contig:
        return 10000
    if mut_a.pos > mut_b.pos:
        print('第二个突变的起始坐标居然小于第一个突变的起始坐标:')
        print(mut_a.__str__())
        print(mut_b.__str__())
        return 10000
    if mut_a.pos == mut_b.pos:
        # 对于同一个位点的不同突变，跳过合并
        return 10000
    if (mut_b.start < mut_a.stop < mut_b.stop) and (get_mutation_type(mut_a) == get_mutation_type(mut_b) == 'deletion'):
        # 针对2个同时为del的突变进行处理
        # 两个突变区域首尾部分重叠
        # 如果重叠区域分别对应的突变序列也一致，则考虑合并
        overlap_len = mut_b.start - mut_a.stop
        if mut_a.alts[0][-overlap_len:] == mut_b.alts[0][:overlap_len]:
            return mut_b.start - mut_a.stop
        else:
            return 10000
    if mut_b.start < mut_a.stop:
        # 可能是MSI导致的
        print(f'忽略特殊突变：后一个突变{mut_b.pos} 在前一个突变的坐标范围内{mut_a.pos}-{mut_a.stop}')
        return 10000
    dis = mut_b.pos - mut_a.pos
    return dis


def get_mutation_type(mutation):
    if mutation.alts is None:
        return None
    if len(mutation.ref) == len(mutation.alts[0]) == 1:
        return "snp"
    if len(mutation.ref) == len(mutation.alts[0]) > 1:
        return "substitution"
    if len(mutation.ref) > len(mutation.alts[0]) and mutation.ref.startswith(mutation.alts[0]):
        return 'deletion'
    elif len(mutation.ref) < len(mutation.alts[0]) and mutation.alts[0].startswith(mutation.ref):
        return 'insertion'
    else:
        return 'complex'


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w+')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger


def find_complex_variant(vcf, out=None, genome=None, bam=None, filter_by_pass=False, logger=None):
    """
    1. 要求vcf是排序好的
    2. 检查相邻的SNP突变的距离是否在2个碱基以内，如果是，如果是，将尝试合并突变
    3. 检查非snp突变的距离是否在50bp以内，如果是，将尝试合并突变
    4. 合并突变时，以测序深度相对较低的记录为主
    :param vcf: input split and left-normalization vcf, 假设最后一个样本为肿瘤样本
    :param out: 输出vcf名称
    :param genome: 参考基因组的fasta文件
    :param bam: bam file
    :param logger: 输出日志
    :param filter_by_pass: 是否仅仅考虑filter==PASS的突变
    """
    fasta = FastaFile(genome) if genome else None
    bam = AlignmentFile(bam) if bam else None
    out = out or (vcf[:-3]+'AddComplex.vcf')
    out2 = out.replace(".AddComplex", "")[:-3] + 'OnlyComplex.vcf'
    if logger is None:
        logger = set_logger(out[:-3]+'log', 'find_complex_variants')

    def read_vcf(infile, keep_pass=filter_by_pass):
        # 只保留filter=PASS和alt不为空的行
        with VariantFile(infile) as temp:
            vcf_header = temp.header
            # if 'MergeFrom' not in header.info: 有些pysam版本中这个判断失效
            vcf_header.info.add('MergeFrom', number=1, type='String', description='Merged Variant Support Number | VariantPos | VariantPos')
            for r in temp:
                if keep_pass:
                    if list(r.filter)[0] == "PASS" and (r.alts is not None):
                        yield r
                else:
                    if r.alts is not None:
                        yield r

    with VariantFile(vcf) as fr:
        header = fr.header
        if 'MergeFrom' not in header.info:
            header.info.add('MergeFrom', number=1, type='String', description='Merged Variant Support Number | VariantPos | VariantPos')
        samples = list(header.samples)
        # 假设最后一个样本是目标样本
        target_idx = len(samples) - 1

    variants = read_vcf(vcf)

    with VariantFile(out, 'w', header=header) as fw, VariantFile(out2, 'w', header=header) as fw2:
        mut_a = next(variants, None)
        if mut_a is None:
            print('! No variant found')
            return
        mut_a_supports = None
        while True:
            merged = False
            # 保证下一次复用的信息是对的
            mut_b_supports = None
            mut_b = next(variants, None)
            if mut_b is None:
                break
            dis = distance(mut_a, mut_b)
            mut_a_type = get_mutation_type(mut_a)
            mut_b_type = get_mutation_type(mut_b)
            mut_a_is_s = mut_a_type in ['snp', 'substitution']
            mut_b_is_s = mut_b_type in ['snp', 'substitution']
            if 'AD' in mut_a.samples[target_idx]:
                a_alt_depth = mut_a.samples[target_idx]['AD'][1]
                b_alt_depth = mut_b.samples[target_idx]['AD'][1]
            elif 'DP' in mut_a.info and 'AF' in mut_a.info:
                a_alt_depth = round(mut_a.info['AF']*mut_a.info['DP'])
                b_alt_depth = round(mut_b.info['AF']*mut_b.info['DP'])
            else:
                a_alt_depth = None
                b_alt_depth = None
            # depth_ratio = a_alt_depth / b_alt_depth

            if dis <= 2 and mut_a_is_s and mut_b_is_s:
                logger.info(f">>>Try to merge mutations at {mut_a.contig}:{mut_a.pos}|{mut_b.pos}")
                logger.info(mut_a.__str__())
                logger.info(mut_b.__str__())
                if mut_a_supports is None:
                    mut_a_supports = get_support_reads(mut_a, mut_a_type, bam, fasta)
                mut_b_supports = get_support_reads(mut_b, mut_b_type, bam, fasta)
                intersection = mut_a_supports & mut_b_supports
                if len(intersection) >= 3:
                    # 这里仅仅合并距离在2bp以内的SNP或等长替换(substitution)
                    ref = fasta.fetch(mut_a.contig, mut_a.start, mut_b.stop)
                    in_low_complex_region = not ref.isupper()
                    ref = ref.upper()
                    alt = mut_a.alts[0] + fasta.fetch(mut_a.contig, mut_a.stop, mut_b.start) + mut_b.alts[0]
                    # 参考基因组的序列中可能存在小写
                    alt = alt.upper()
                    # record = mut_a.copy() if mut_a.info['DP'] <= mut_b.info['DP'] else mut_b.copy()
                    record = mut_a.copy() if a_alt_depth < b_alt_depth else mut_b.copy()
                    record.pos = mut_a.pos
                    record.ref = ref
                    record.alts = [alt]
                    merged = True

                    if in_low_complex_region:
                        logger.info("warn: this mutation is likely in repeating region")
                    logger.info(f"--Vcf reported supporting read number:{a_alt_depth} vs {b_alt_depth}")
                    logger.info(f"--The program saw supporting read number:{len(mut_a_supports)} vs {len(mut_b_supports)}")
                    # logger.info(f"Both supporting read number: {len(intersection)}, egg. {intersection.pop()}")
                    logger.info(f"--Both variant supporting read number:{len(intersection)}, they are {intersection}")

                    # 搜索支持merge突变的reads
                    mut_c_supports = get_support_reads(record, get_mutation_type(record), bam, fasta)
                    logger.info(f'--Merged variant supporting read number:{len(mut_c_supports)}, they are {mut_c_supports}')
                    if len(mut_c_supports) < 2 and len(alt) >= 3:
                        logger.info('支持合并后的突变的reads少于2个，说明2个snp中间也是snp，但原vcf中没有')
                        mid_alt, mid_req, mid_depth = detect_snp(bam, mut_a.contig, mut_a.start + 1, fasta_file=fasta)
                        record.alts = [alt[0] + mid_alt + alt[2]]
                        logger.info(f'现在根据bam将中间的snp报告出来:{alt[1]}>{mid_alt}, {mid_req}/{mid_depth}')
                        mut_c_supports = get_support_reads(record, get_mutation_type(record), bam, fasta)
                        logger.info(f'New Merged variant supporting read number:{len(mut_c_supports)}, they are {mut_c_supports}')

                    if 'MergeFrom' in mut_a.info:
                        logger.info('Super merging found!')
                        record.info['MergeFrom'] = mut_a.info['MergeFrom'] + '|' + str(mut_b.pos)
                    else:
                        record.info['MergeFrom'] = '|'.join([str(len(mut_c_supports)), str(mut_a.pos), str(mut_b.pos)])
                    fw2.write(record)
                    logger.info("Suggested Complex/Simple formation as follow:")
                    logger.info(record.__str__())
                    logger.info('------------------------------------------------------------------')
                    # 思考：如果这里不直接将合并结果写出，而是把record作为mut_b传递，是否就可以实现滚动合并2个以上的突变？
                    # merge得到的record会进入到下一轮作为mut_a与新读取的mut_b进一步比较是否可以合并
                    mut_b = record
                else:
                    logger.info(f"Give up merging for few supporting read {intersection}")
                    logger.info("--------------------------------------------------------------")

            elif dis <= 50 and not (mut_a_is_s and mut_b_is_s):
                # 对于距离超过2且突变类型都是snp或等长替换的，均不考虑合并突变
                logger.info(f">>>Try to merge mutations at {mut_a.contig}:{mut_a.pos}|{mut_b.pos}")
                logger.info(mut_a.__str__())
                logger.info(mut_b.__str__())
                ref = fasta.fetch(mut_a.contig, mut_a.start, mut_b.stop)
                in_low_complex_region = not ref.isupper()
                ref = ref.upper()
                if dis > 0:
                    alt = mut_a.alts[0] + fasta.fetch(mut_a.contig, mut_a.stop, mut_b.start) + mut_b.alts[0]
                    alt = alt.upper()
                else:
                    logger.info('The above two deletion variants are part overlapped')
                    alt = mut_a.alts[0]
                new_pos, new_ref, new_alt, alignment = global_pair_seq_align(mut_a.pos, ref, alt)
                logger.info(alignment.__str__())
                if new_pos is not None:
                    if mut_a_supports is None:
                        mut_a_supports = get_support_reads(mut_a, mut_a_type, bam, fasta)
                    mut_b_supports = get_support_reads(mut_b, mut_b_type, bam, fasta)
                    intersection = mut_a_supports & mut_b_supports
                    record = mut_a.copy() if a_alt_depth < b_alt_depth else mut_b.copy()
                    record.pos = new_pos
                    # print(record.ref, pair_align)
                    record.ref = new_ref
                    record.alts = [new_alt]
                    # 搜索支持merge突变的reads
                    mut_c_supports = get_support_reads(record, get_mutation_type(record), bam, fasta)
                    if 'MergeFrom' in mut_a.info:
                        logger.info('Super merging found!')
                        record.info['MergeFrom'] = mut_a.info['MergeFrom'] + '|' + str(mut_b.pos)
                    else:
                        record.info['MergeFrom'] = '|'.join([str(len(mut_c_supports)), str(mut_a.pos), str(mut_b.pos)])

                    if in_low_complex_region:
                        logger.info("warn: this mutation is likely in repeating region")
                    logger.info(f"--Vcf reported supporting read number:{a_alt_depth} vs {b_alt_depth}")
                    logger.info(f"--The program saw supporting read number:{len(mut_a_supports)} vs {len(mut_b_supports)}")
                    # logger.info(f"Both supporting read number: {len(intersection)}, egg. {intersection.pop()}")
                    logger.info(f"--Both variant supporting read number:{len(intersection)}, they are {intersection}")
                    logger.info(f'--Merged variant supporting read number:{len(mut_c_supports)}, they are {mut_c_supports}')
                    if (len(intersection) >= 3 and len(mut_c_supports) >= 2) or (len(mut_c_supports) >= 20):
                        merged = True
                        fw2.write(record)
                        logger.info("Suggested Complex/Simple formation as follow:")
                        logger.info(record.__str__())
                        logger.info("--------------------------------------------------------------")
                        # merge得到的record会进入到下一轮作为mut_a与新读取的mut_b进一步比较是否可以合并
                        mut_b = record
                    else:
                        logger.info(f"Give up merging for few supporting read {intersection} and {mut_c_supports}")
                        logger.info(f'Possible merged formation as follow:\n{record.__str__()}')
                        logger.info("--------------------------------------------------------------")
                else:
                    logger.info('Based on alignment result, give up merging!')
                    logger.info("--------------------------------------------------------------")

            if not merged:
                # 写入没有被合并的突变mut_a, 他可能是之前成功合并得到的record，也可能是原始record
                fw.write(mut_a)

            # 把mut_b替换为mut_a,用于下一轮比较
            mut_a = mut_b

            # 下面的操作是为了复用mut_b_supports
            if not merged:
                mut_a_supports = mut_b_supports
            else:
                # 发生了merge，则不能复用mut_b_supports
                mut_a_supports = None

        # 循环结束, 如果最后一轮如果没有merge发生，需要把最后读取的一个突变写入
        if not merged and (mut_a is not None):
            fw.write(mut_a)
        bam.close()

    return out, out2


def add_complex_variant_to_tso500_report(infile, out=None, genome=None, bam=None, cache_dir=None):
    base_name = os.path.basename(infile)
    outdir = os.path.dirname(out) if out else os.getcwd()
    vcf = small_variant_tsv_to_vcf(infile, out=os.path.join(outdir, base_name[:-3] + 'SmallVariant.vcf'))
    _, vcf = find_complex_variant(vcf, genome=genome, bam=bam, filter_by_pass=False)
    vep_vcf = vep_annotation(vcf, fasta=genome, cache_dir=cache_dir)
    add_vep_annot_back_to_tsv(vep_vcf, infile, out=out)


if __name__ == '__main__':
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', type=Path, required=True, help='vcf file annotated with vep. 假设最后一个样本是肿瘤样本')
    parser.add_argument('-genome', type=Path, required=True, help='path to indexed genome fasta')
    parser.add_argument('-bam', type=Path, required=True, help='bam file')
    parser.add_argument('--filter_by_pass', default=False, action='store_true', help='if filter record which is flagged as PASS')
    parser.add_argument('-out', required=True,  help='output vcf file')
    args = parser.parse_args()
    find_complex_variant(args.vcf, out=args.out, genome=args.genome, bam=args.bam, filter_by_pass=args.filter_by_pass)
