"""
测序错误过滤突变的思路：
基于置信区间，可以参考此链接理解统计原理：https://www.statisticshowto.com/binomial-confidence-interval/
假设事先知道测序过程中发生的碱基转换的概率，例如，知道A被错误测成T的概率为0.001，
针对一个突变位点A-T，其AF=0.008，是否可以将当前位点判定为假阳性呢?
假设还知道当前位点深度为1000，根据这1000 Reads，在95%的置信水平下，估计出来的测序错误概率与真实的错误概率偏差可以使用如下公式推算：
（下面的公式是根据正态分布来的，即用正太分布近似二型分布，理论依据是当N和p都比较大时，二型分布趋近正太分布。）
N = P*(1-P)*(Z/E)**2
=> E = Z/(N/P/(1-P))**0.5
即：E = 1.96/(N/(P*(1-P)))**0.5
99的置信水平：E = 2.58/(N/(P*(1-P)))**0.5
代入 N=1000, P=0.001得 E=0.00196，也就是，如果A>T是测序错误导致的，那么这个错误频率的置信区间为
[P-E, P+E] = [0, 0.00296], 由于现在测得的AF=0.008 已经远在置信区间水平以外，所以不应该过滤掉。

实际上，在python的scipy包中：
    from scipy.stats import poisson, binom, norm
    # binom.interval(confidence, n, p, loc=0)
    # 可以看到p的上限是3.0/1000 = 0.003，和基于正太分布的计算公式结果非常相近
    binom.interval(0.95, n=1000, p=0.001)
    Out: (0.0, 3.0)

    p-value = 1 - binom.cdf(k=8, n=1000, p=0.001)
    Out: 1.0936095204971963e-06

思考2：
当有对照样本时，且对照样本中该位点也存在突变，设其突变频率为ctrP, 如果使用上述方法判定该突变为假阳，
则可以使用该P作为真实测序错误率对测试样本进行上述类似的统计.

思考3：
当没有对照样本时，而只有阴性设计样本如NA12878时，我们可以通过统计阴性样本中碱基测序错误率作为上述过滤思路的输入

思考4：
当既没有对照样本，又没有阴性样本，我们假设肿瘤样本中call出来的低频突变中绝大部分为假阳性突变，那么我们也可以根据肿瘤样本粗略估计测序错误。
最后可以用肿瘤样本自己作为输入并实现过滤。

新思考，确定LOD，如何估计测序错误率：
设一次测序错误理论概率为P
由样本量估计公式 N = P*(1-P)*(Z/E)**2 推导得到N次测序估计出来的误差幅度为：E = Z/(N/P/(1-P))**0.5
同样，已知某个突变的理论值时，也可以应用上述公式。

设LOD为检测下限，假设某个突变的真实频率为L，为减少假阳性，在某次实验中，我们应该要求“AF最低估计值”大于“测序错误最大估计值”，即：
=> L > 测序错误的上限
=> L - Ei > P + E
=>  L - Z/(N/L/(1-L))**0.5 > P + Z/(N/P/(1-P))**0.5
=> (L-P) > Z/(N/P/(1-P))**0.5 + Z/(N/L/(1-L))**0.5
=> N > ( ( Z/(1/P/(1-P))**0.5 + Z/(1/L/(1-L))**0.5 ) / (L-P) )**2

99%置信区间对应的标准正态分布统计量Z=2.58
95%置信区间对应的标准正态分布统计量Z=1.96

P=0.001; L=0.01; Z=1.96
N = ( ( Z/(1/P/(1-P))**0.5 + Z/(1/L/(1-L))**0.5 ) / (L-P) )**2 = 815.21
即95%置信区间时，要保证LOD=0.01,测序深度需要大于816

P=0.001; L=0.01; Z=2.58
N = ( ( Z/(1/P/(1-P))**0.5 + Z/(1/L/(1-L))**0.5 ) / (L-P) )**2 = 1412.53
即选择99%置信区间时，要保证LOD=0.01,测序深度需要大于1413

新思考，下面的过滤是否合理？
1. 提供error_rate, depth, confidence，计算得到error_rate的估计上限
2. 根据error_rate, AF, confidence, 计算得到最低测序深度min_alt_depth
3. 根据error_rate_upper 和 min_alt_depth同时进行突变过滤

"""
import math
import os
import re
import json
import logging
import pandas as pd
import scipy.stats as stats
import statistics
from collections import Counter
import pysam
from pysam import FastaFile, VariantFile, AlignmentFile, VariantHeader, VariantRecord
from Bio import Align
from stat_3bases_error import estimate_context_seq_error
__author__ = 'gdq'
"""
要求vcf每一行只包含一个突变，这个可以通过bcftools norm 快速实现
变异类型的分类
https://www.ebi.ac.uk/training-beta/online/courses/human-genetic-variation-introduction/what-is-genetic-variation/types-of-genetic-variation/

考虑MSI的识别, 这样可以过滤MSI突变
"""


def set_logger(name='log.info', logger_id='x'):
    logger = logging.getLogger(logger_id)
    fh = logging.FileHandler(name, mode='w+')
    logger.addHandler(fh)
    logger.setLevel(logging.INFO)
    return logger


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


class ValidateMutationByBam(object):
    def __init__(self, bam_file, genome_file):
        self.bam = AlignmentFile(bam_file)
        self.genome = FastaFile(genome_file)

    def get_pileup_column_tags(self, pileup_column, tag_name='cD'):
        tags = []
        for pileup_read in pileup_column.pileups:
            tag = pileup_read.alignment.get_tag(tag_name)
            tags.append(tag)
        return tags

    def get_pileup_column_cigars(self, pileup_column):
        tags = []
        for pileup_read in pileup_column.pileups:
            tag = pileup_read.alignment.cigarstring
            tags.append(tag)
        return tags

    @staticmethod
    def qual_diff_test(alt_reads:set, qual_dict:dict, qual_dict2:dict=None):
        if qual_dict2 is None:
            alt_quals = [qual_dict[r] for r in alt_reads if r in qual_dict]
            ref_quals = [qual_dict[r] for r in (qual_dict.keys() - alt_reads)]
        else:
            # qual_dict2是旁边位置信息，因此需要取交集，进行碱基质量配对比较分析
            alt_reads = alt_reads & qual_dict.keys() & qual_dict2.keys()
            alt_reads = list(alt_reads)
            alt_quals = [qual_dict[r] for r in alt_reads]
            ref_quals = [qual_dict2[r] for r in alt_reads]

        alt_quals = [x for x in alt_quals if x is not None]
        ref_quals = [x for x in ref_quals if x is not None]
        if len(alt_quals) > 8 and len(ref_quals) > 8:
            # u, pvalue = stats.mannwhitneyu(alt_quals, ref_quals, alternative='less')
            u, pvalue = stats.ranksums(alt_quals, ref_quals)
            # u, pvalue = stats.ttest_ind(alt_quals, ref_quals)
        elif len(alt_quals) >= 3 and len(alt_quals) >= 3:
            if qual_dict2:
                # 配对检验
                s, pvalue = stats.ttest_rel(alt_quals, ref_quals)
            else:
                # 非配对检验
                s, pvalue = stats.ttest_ind(alt_quals, ref_quals)
        else:
            # 统计意义不足
            pvalue = 1.0
        # 对为nan的pvalue值进行转换
        if pvalue != pvalue:
            pvalue = 1.0
        return pvalue, alt_quals, ref_quals

    @staticmethod
    def consensus_seqs(seqs):
        # 对支持突变的reads的突变中心进行consensus
        if len(seqs) == 1:
            return seqs[0]
        else:
            # ''.join(Counter(bases).most_common(1)[0][0] for bases in zip(*seqs))
            # 考虑多个分支
            total = len(seqs)
            consensus = ['']
            for bases in zip(*seqs):
                two_bases = Counter(bases).most_common()
                base1, freq1 = two_bases[0]
                if len(two_bases) > 1:
                    base2, freq2 = two_bases[1]
                else:
                    base2, freq2 = '', 0
                tmp_consensus = []
                for each in consensus:
                    tmp_consensus.append(each + base1)
                    if freq2 / total >= 0.4:
                        tmp_consensus.append(each+base2)
                consensus = tmp_consensus
            return consensus

    @staticmethod
    def foxog_calc(ref, F1R2, F2R1):
        """
        C>A|G>T variants can be oxidation artifacts introduced during library preparation (Costello et al, 2013).
        These "OxoG" artifacts have a telltale read orientation, with the majority of ALT reads in the artifact orientation.
        """
        FoxoG = None
        if (F1R2 + F2R1) > 6:
            if ref == "C" or ref == "A":
                FoxoG = F2R1 / (F1R2 + F2R1)
            elif ref == "G" or ref == "T":
                FoxoG = F1R2 / (F1R2 + F2R1)
            FoxoG = round(FoxoG, 3)
        return FoxoG

    def get_snp_support_reads(self, contig, start, ref, alt, tag_names:tuple=tuple(), min_bq=13, logger=None):
        pileup_columns = self.bam.pileup(
            contig, start-1, start + 2,  # 提取左右碱基，方便后续突变碱基与左右碱基质量的比较
            stepper='samtools',
            truncate=True,
            min_base_quality=min_bq,
            ignore_orphans=False,
            ignore_overlaps=True,  # set 为True则意味着取质量高的base作为代表
            # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
            max_depth=300000,
            fastafile=self.genome,
            # flag_filter: The default is BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP.
        )
        support_reads = set()
        all_read_tag_dict = dict()
        support_read_seqs = []
        base_qual_dict = dict()
        left_base_qual_dict = dict()
        mid_base_qual_dict = dict()
        right_base_qual_dict = dict()
        query_pos_dict = dict()
        map_qual_dict = dict()
        # init value
        mean_alt_bq = None
        mean_left_bq = None
        mean_right_bq = None
        mean_alt_pos = None
        alt_pos_pstd = None
        ref_pos_std = None,
        mean_alt_mq = None,
        bq_pvalue = None
        left_bq_pvalue = None
        right_bq_pvalue = None
        mq_pvalue = None
        pos_pvalue = None
        mean_ref_bq = None
        mean_ref_mq = None
        mean_ref_pos = None
        max_alt_pos_diff = None
        # F=foward strand, R=reverse strand, 1:read1, 2:read2
        F2R1 = 0
        F1R2 = 0
        FoxoG = None
        # init value end
        for idx, col in enumerate(pileup_columns):
            base_quals = col.get_query_qualities()
            query_names = col.get_query_names()
            base_qual_dict = dict(zip(query_names, base_quals))
            if idx == 0:
                left_base_qual_dict = base_qual_dict
            elif idx == 2:
                right_base_qual_dict = base_qual_dict
            else:
                map_quals = col.get_mapping_qualities()
                map_qual_dict = dict(zip(query_names, map_quals))
                mid_base_qual_dict = base_qual_dict
                # 获取支持变异的reads/序列/所有read的相关tag信息
                for pileup_read in col.pileups:
                    alignment = pileup_read.alignment
                    read_name = alignment.query_name
                    query_seq = alignment.query_sequence
                    query_pos = alignment.query_position
                    query_len = alignment.query_length
                    if query_pos is not None:
                        query_pos_dict[read_name] = min(query_pos, query_len - query_pos)
                        if query_seq[query_pos].upper() == alt:
                            # 由于未知原因，这里找到的read_name不一定在前面的query_names中
                            support_reads.add(read_name)
                            # 提取围绕SNV位点前后5个碱基，用于检查支持突变的read之间的一致性
                            if (query_pos >= 5) and (query_len > query_pos+5):
                                support_read_seqs.append(query_seq[query_pos-5:query_pos+5])
                            # Counters for FoxoG calculations
                            if alignment.is_proper_pair:
                                if alignment.is_read1:
                                    if alignment.is_reverse:
                                        F2R1 += 1
                                    else:
                                        F1R2 += 1
                                elif alignment.is_read2:
                                    if alignment.is_reverse:
                                        F1R2 += 1
                                    else:
                                        F2R1 += 1
                    # 由于要获取下面的tag，才不得不对pileups进行循环
                    if read_name not in all_read_tag_dict:
                        target_tag_dict = dict()
                        for tag in tag_names:
                            target_tag_dict[tag] = pileup_read.alignment.get_tag(tag)
                        all_read_tag_dict[read_name] = target_tag_dict
        # 开始统计分析
        if support_reads:
            # try:
            # 检验支持alt的碱基质量和不支持alt的碱基质量差异
            bq_pvalue, alt_bqs, ref_bqs = self.qual_diff_test(support_reads, mid_base_qual_dict)
            # 检验支持alt的read质量和不支持alt的read的比对质量差异
            mq_pvalue, alt_mqs, ref_mqs = self.qual_diff_test(support_reads, map_qual_dict)
            # 检验支持alt的碱基位置和不支持alt的碱基位置的差异
            pos_pvalue, alt_pos, ref_pos = self.qual_diff_test(support_reads, query_pos_dict)
            mean_alt_bq = statistics.mean(alt_bqs) if alt_mqs else None
            mean_ref_bq = statistics.mean(ref_bqs) if ref_bqs else None
            mean_alt_mq = statistics.mean(alt_mqs) if alt_mqs else None
            mean_ref_mq = statistics.mean(ref_mqs) if ref_mqs else None
            mean_alt_pos = statistics.median(alt_pos) if alt_pos else None
            mean_ref_pos = statistics.median(ref_pos) if ref_pos else None
            alt_pos_pstd = statistics.pstdev(alt_pos) if alt_pos else None
            max_alt_pos_diff = max(alt_pos) - min(alt_pos)
            ref_pos_std = statistics.pstdev(ref_pos) if ref_pos else None
            FoxoG = self.foxog_calc(ref, F1R2, F2R1)

        if left_base_qual_dict and support_reads:
            left_bq_pvalue, _, left_bqs = self.qual_diff_test(support_reads, mid_base_qual_dict, left_base_qual_dict)
            if left_bqs:
                mean_left_bq = statistics.mean(left_bqs)
            else:
                print('Left base info lost', [query_pos_dict[x] for x in support_reads])
                print(contig, start, ref, alt, support_reads)

        if right_base_qual_dict and support_reads:
            right_bq_pvalue, _, right_bqs = self.qual_diff_test(support_reads, mid_base_qual_dict, right_base_qual_dict)
            if right_bqs:
                mean_right_bq = statistics.mean(right_bqs)
            else:
                print('Right base info lost', [query_pos_dict[x] for x in support_reads])
                print(contig, start, ref, alt, support_reads)
        # 判断snv的质量
        good_snp = True
        if bq_pvalue and left_bq_pvalue and right_bq_pvalue:
            reasonable_pos = mean_alt_pos >= 5 and max_alt_pos_diff > 1
            small_bq_diff = (bq_pvalue > 0.001) and (left_bq_pvalue > 0.05) and (right_bq_pvalue > 0.05)
            small_mq_diff = (mq_pvalue > 0.001) if mq_pvalue else True
            good_alt_bq = mean_alt_bq >= 30
            OxoG = FoxoG < 0.9 if FoxoG else False
            if not (reasonable_pos and small_bq_diff and small_mq_diff and good_alt_bq and OxoG):
                good_snp = False
        all_read_tag_dict['snp_is_good'] = good_snp
        if logger:
            info = dict(
                contig=contig,
                start=start+1,
                ref=ref,
                alt=alt,
                alt_dp=len(support_reads),
                total_dp=len(base_qual_dict),
                mean_left_bq=mean_left_bq,
                mean_alt_bq=mean_alt_bq,
                mean_right_bq=mean_right_bq,
                mean_ref_bq=mean_ref_bq,
                mean_alt_pos=mean_alt_pos,
                max_alt_pos_diff=max_alt_pos_diff,
                FoxoG=FoxoG,
                mean_ref_pos=mean_ref_pos,
                alt_pos_pstd=alt_pos_pstd,
                ref_pos_std=ref_pos_std,
                mean_alt_mq=mean_alt_mq,
                mean_ref_mq=mean_ref_mq,
                bq_pvalue=bq_pvalue,
                left_bq_pvalue=left_bq_pvalue,
                right_bq_pvalue=right_bq_pvalue,
                mq_pvalue=mq_pvalue,
                pos_pvalue=pos_pvalue,
                good_snp=good_snp,
            )
            logger.info(str(info))
        return support_reads, all_read_tag_dict, support_read_seqs

    def get_insert_support_reads(self, contig, start, alt, tag_names:tuple=tuple(), min_bq=13, logger=None):
        cols = self.bam.pileup(
            contig, start, start + 1,
            stepper='samtools',
            truncate=True,
            min_base_quality=min_bq,
            ignore_orphans=False,
            ignore_overlaps=True,  # set 为True则意味着取质量高的base作为代表
            # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
            max_depth=300000,
            fastafile=self.genome,
            # flag_require=read_type,
            # flag_filter=4,
        )
        # 对于处于重复区域的插入，可能有多种比对形式的记录，因此我们找的证据可能偏少
        support_reads = set()
        all_read_tag_dict = dict()
        expected_insert = alt.upper()
        alt_len = len(alt)
        support_read_seqs = []
        for col in cols:
            for base, read in zip(
                    # 如不加add_indels参数，那么将无法知晓插入的碱基序列
                    col.get_query_sequences(add_indels=True),
                    col.get_query_names(),
            ):
                if '+' in base:
                    # such as 'G+2AT', +号前面是参考序列
                    insertion = re.sub(r'\+\d+', '', base).upper()
                    # 暂时几乎不考虑insertion恰巧出现在read的两端且不包含完整的insertion
                    mismatch_allowed = round(len(expected_insert) * 0.15)
                    mismatched_num = sum(x != y for x, y in zip(expected_insert, insertion))
                    if mismatched_num <= mismatch_allowed:
                        support_reads.add(read)

            for pileup_read in col.pileups:
                target_tag_dict = dict()
                read_name = pileup_read.alignment.query_name
                query_seq = pileup_read.alignment.query_sequence
                query_pos = pileup_read.query_position
                # clip的base或许是插入的序列
                if alt_len > 2:
                    cigar = pileup_read.alignment.cigarstring
                    if cigar.endswith('S'):
                        # read右端clip的序列是insertion的起始
                        query_alignment_end = pileup_read.alignment.query_alignment_end
                        if expected_insert[1:].startswith(query_seq[query_alignment_end:].upper()):
                            support_reads.add(read_name)
                    if pileup_read.alignment.cigartuples[0][0] == 4:
                        # read左端clip掉的序列是insertion尾部
                        query_alignment_start = pileup_read.alignment.query_alignment_start
                        if expected_insert.endswith(query_seq[:query_alignment_start]):
                            support_reads.add(read_name)
                if read_name not in all_read_tag_dict:
                    for tag in tag_names:
                        target_tag_dict[tag] = pileup_read.alignment.get_tag(tag)
                    # 不关心query_name是否重复，以最后一个为准
                    all_read_tag_dict[read_name] = target_tag_dict
                # 提取支持突变的read的序列信息，围绕突变周围展开
                if (query_pos is not None) and (read_name in support_reads):
                    if (query_pos >= 4) and (len(query_seq) > query_pos+alt_len+4):
                            target = query_seq[query_pos-4: query_pos] + query_seq[query_pos+alt_len:query_pos+alt_len+4]
                            support_read_seqs.append(target)

        return support_reads, all_read_tag_dict, support_read_seqs

    def get_insert_support_reads_v2(self, contig, start, alt, tag_names:tuple=tuple(), min_bq=13, logger=None):
        cols = self.bam.pileup(
            contig, start, start + 1,
            stepper='samtools',
            truncate=True,
            min_base_quality=min_bq,
            ignore_orphans=False,
            ignore_overlaps=True,  # set 为True则意味着取质量高的base作为代表
            # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
            max_depth=300000,
            fastafile=self.genome,
            # flag_require=read_type,
            # flag_filter=4,
        )
        # 对于处于重复区域的插入，可能有多种比对形式的记录，因此我们找的证据可能偏少
        # 我们围绕插入序列，提取出50bp，根据已经获得的证据进行consensus，作为alt参考序列
        # 对插入位点前后50bp的区域进行read检查，如果read包含alt参考序列，则判定为突变证据
        support_reads = set()
        all_read_tag_dict = dict()
        expected_insert = alt.upper()
        alt_len = len(alt)
        alt_centered_seqs = []
        alt_flank_seqs = []
        all_read_names = set()
        for col in cols:
            for base, read in zip(
                    # 如不加add_indels参数，那么将无法知晓插入的碱基序列
                    col.get_query_sequences(add_indels=True),
                    col.get_query_names(),
            ):
                if '+' in base:
                    # such as 'G+2AT', +号前面是参考序列
                    insertion = re.sub(r'\+\d+', '', base).upper()
                    # 暂时几乎不考虑insertion恰巧出现在read的两端且不包含完整的insertion
                    mismatch_allowed = round(len(expected_insert) * 0.15)
                    mismatched_num = sum(x != y for x, y in zip(expected_insert, insertion))
                    if mismatched_num <= mismatch_allowed:
                        support_reads.add(read)

            for pileup_read in col.pileups:
                read_name = pileup_read.alignment.query_name
                query_seq = pileup_read.alignment.query_sequence.upper()
                query_pos = pileup_read.query_position
                all_read_names.add(read_name)
                # 进一步提取可能的证据，clip的base或许是插入的序列
                if alt_len >= 2:
                    cigar = pileup_read.alignment.cigarstring
                    if cigar.endswith('S'):
                        # read右端clip的序列是insertion的起始
                        query_alignment_end = pileup_read.alignment.query_alignment_end
                        if expected_insert[1:].startswith(query_seq[query_alignment_end:]):
                            support_reads.add(read_name)
                    if pileup_read.alignment.cigartuples[0][0] == 4:
                        # read左端clip掉的序列是insertion尾部
                        query_alignment_start = pileup_read.alignment.query_alignment_start
                        if expected_insert.endswith(query_seq[:query_alignment_start]):
                            support_reads.add(read_name)

                # 提取支持突变的read的序列信息，围绕突变周围展开
                if (query_pos is not None) and (read_name in support_reads):
                    if (query_pos >= 22) and (len(query_seq) > query_pos+alt_len+22):
                            alt_centered_seqs.append(query_seq[query_pos-22:query_pos+alt_len+22])
                    if (query_pos >= 5) and (len(query_seq) > query_pos+alt_len+5):
                            flank5 = query_seq[query_pos-5: query_pos] + query_seq[query_pos+alt_len:query_pos+alt_len+5]
                            alt_flank_seqs.append(flank5)

        # 进一步捞回其他表达形式的insertion证据
        if alt_centered_seqs:
            # 如果插入位点附近还有1个或多个杂合碱基，那么consensus序列将只能代表其中一个
            consenus_alt_seqs = self.consensus_seqs(alt_centered_seqs)
            for each in self.bam.fetch(contig, start-50, start+51):
                read_name = each.query_name
                if read_name in support_reads:
                    continue
                read_seq = each.query_sequence.upper()
                # 查看read是否包含目标序列
                for each in consenus_alt_seqs:
                    if each in read_seq:
                        support_reads.add(read_name)
                        break
                # 提取tag信息
                if (read_name in all_read_names) or (read_name in support_reads):
                    target_tag_dict = dict()
                    if read_name not in all_read_tag_dict:
                        for tag in tag_names:
                            target_tag_dict[tag] = each.get_tag(tag)
                        # 不关心query_name是否重复，以最后一个为准
                        all_read_tag_dict[read_name] = target_tag_dict
        # 将consensus_alt_seq和参考基因组序列进行比对，能找到

        return support_reads, all_read_tag_dict, alt_flank_seqs

    def get_del_support_reads(self, contig, start, del_len, tag_names:tuple=tuple(), min_bq=13, ignore_overlaps=True, logger=None):
        cols = self.bam.pileup(
            contig, start, start + 1,
            stepper='samtools',
            truncate=True,
            min_base_quality=min_bq,
            ignore_orphans=False,
            ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
            # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
            max_depth=300000,
            fastafile=self.genome,
            # flag_filter=4,
        )
        # 对于一个del，比对方式多种多样
        # ref: ACCTTTATGCTG
        # alt: ACCTTT---CTG (中间删除3个碱基）
        # alt: ACCTTTCTG    (中间是2个替换）
        # 不能alt的形式如何变幻，终究表示的是一个突变，也就是本质包含同一个突变的reads序列应该一致
        support_reads = set()
        all_read_tag_dict = dict()
        support_read_seqs = []
        for col in cols:
            for pileup_read in col.pileups:
                read_name = pileup_read.alignment.query_name
                query_seq = pileup_read.alignment.query_sequence
                query_pos = pileup_read.query_position_or_next
                if -pileup_read.indel == del_len:
                    # 由于indel的比对方式有多种，这种方法找到的del证据偏少，需要进行重新比对才能找到足够的证据，目前为未实现
                    support_reads.add(pileup_read.alignment.query_name)
                    if (query_pos >= 4) and (len(query_seq) > query_pos + 4):
                        support_read_seqs.append(query_seq[query_pos - 4:query_pos + 4])
                target_tag_dict = dict()
                if read_name not in all_read_tag_dict:
                    for tag in tag_names:
                        target_tag_dict[tag] = pileup_read.alignment.get_tag(tag)
                    all_read_tag_dict[read_name] = target_tag_dict
        return support_reads, all_read_tag_dict, support_read_seqs

    def get_substitution_support_reads(self, contig, start, alt, tag_names:tuple=tuple(), min_bq=13, ignore_overlaps=True, logger=None):
        # 该函数可以被get_complex_support_reads替代
        cols = self.bam.pileup(
            contig, start, start + 1,
            stepper='samtools',
            truncate=True,
            min_base_quality=min_bq,
            ignore_orphans=False,
            ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
            # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
            max_depth=300000,
            fastafile=self.genome,
            # flag_filter=4,
        )
        alt = alt.upper()
        support_reads = list()
        all_read_tag_dict = dict()
        for idx, col in enumerate(cols):
            # for base, read in zip(
            #         # 如不加add_indels参数，那么将无法知晓插入的碱基序列
            #         col.get_query_sequences(add_indels=True),
            #         col.get_query_names(),
            # ):
            #     if base.upper() == alt[idx]:
            #         support_reads.append(read)
            for pileup_read in col.pileups:
                query_seq = pileup_read.alignment.query_sequence
                query_pos = pileup_read.query_position
                read_name = pileup_read.alignment.query_name
                if query_pos is not None:
                    if query_seq[query_pos].upper() == alt[idx]:
                        support_reads.append(read_name)
                target_tag_dict = dict()
                if read_name not in all_read_tag_dict:
                    for tag in tag_names:
                        target_tag_dict[tag] = pileup_read.alignment.get_tag(tag)
                    # 不关心query_name是否重复，以最后一个为准
                    all_read_tag_dict[read_name] = target_tag_dict

        # 如果一个read匹配到的次数等于alt的长度，则认为该read支持alt
        count = Counter(support_reads)
        alt_length = len(alt)
        support_reads = set(x for x, y in count.items() if y == alt_length)
        return support_reads, all_read_tag_dict

    def get_complex_support_reads(self, contig, start, ref, alt, tag_names:tuple=tuple(), min_bq=13, ignore_overlaps=True, logger=None):
        """
        思路：分析每条能够比对到突变起始位置的read序列是否支持alt
        1. 假设alt前后的3个碱基一定是和reference匹配的
        2. 对于complex突变，可以理解为是：删除ref，加入alt, 假设ref=ATC，alt=GT，ref前后序列是XY，那么突变的read应该是XGTY
        3. 如果 X + alt + Y = read_alt_left + alt + read_alt_right
        注意：
            1. 该思路不适合重复区域的扩增或缺失情形，因为扩增或缺失后，突变前后，序列仅仅发生了长短发生变化
            2. 因为如果要用，还得扩大XY到重负区域以外才可以
        """
        cols = self.bam.pileup(
            # 提取ref的前一个碱基对应的比对信息
            contig, start, start + 1,
            stepper='samtools',
            truncate=True,
            min_base_quality=min_bq,
            ignore_orphans=False, # 不会忽略单端read
            ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
            # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
            max_depth=300000,
            fastafile=self.genome,
        )
        support_reads = set()
        all_read_tag_dict = dict()
        expected_left = self.genome.fetch(contig, start - 3, start)
        expected_right = self.genome.fetch(contig, start + len(ref), start + len(ref) + 3)
        alt_len = len(alt)
        support_read_seqs = []
        for col in cols:
            for pileup_read in col.pileups:
                query_seq = pileup_read.alignment.query_sequence
                query_len = pileup_read.alignment.query_length
                query_pos = pileup_read.query_position
                read_name = pileup_read.alignment.query_name
                target_tag_dict = dict()
                if read_name not in all_read_tag_dict:
                    for tag in tag_names:
                        target_tag_dict[tag] = pileup_read.alignment.get_tag(tag)
                    # 不关心query_name是否重复，以最后一个为准
                    all_read_tag_dict[read_name] = target_tag_dict
                if query_pos is not None:
                    # position of the read base at the pileup site, 0-based. None if is_del or is_refskip is set.
                    # 我们期望该位置下一位|mismatch|insertion|deletion|
                    # 提取真实read信息： 期望是（3个和参考一致的碱基+可能的插入序列+3个和参考一致的碱基） | 6个和参考一致的碱基
                    if query_pos >= 3:
                        # read的第4个或之后的碱基比对到当前位点
                        back_extend = 3
                    else:
                        # read的第一个碱基比对到当前位点
                        back_extend = query_pos

                    alt_end_to_read_end = query_len - query_pos - len(alt)
                    if alt_end_to_read_end >= 3:
                        forward_extend = 3
                    else:
                        forward_extend = alt_end_to_read_end

                    real_seq = query_seq[query_pos - back_extend:query_pos + len(alt) + 3]
                    expected_seq = expected_left[-back_extend:] + alt + expected_right[:forward_extend]
                    # print("xxx", contig, start, real_seq, alt, expected_seq)
                    if real_seq.lower() == expected_seq.lower():
                        support_reads.add(read_name)
                        # 提取支持突变的read的序列信息，围绕突变周围展开
                        if (query_pos >= 4) and (len(query_seq) > query_pos + alt_len + 4):
                            target = query_seq[query_pos - 4: query_pos] + query_seq[query_pos + alt_len:query_pos + alt_len + 4]
                            support_read_seqs.append(target)

                        # print(query_seq, query_pos, expected_seq_before, ref+'>'+alt, expected_seq_after)
        return support_reads, all_read_tag_dict, support_read_seqs

    def get_all_read_tag_dict(self, mut, min_bq=13, tag_names=('cD',), ignore_overlaps=True):
        pileup_columns = self.bam.pileup(
            mut.contig, mut.start, mut.start + len(mut.alts[0]),
            stepper='samtools',
            truncate=True,
            min_base_quality=min_bq,
            ignore_orphans=False,
            ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
            # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
            max_depth=300000,
            fastafile=self.genome,
            # flag_filter=4,
        )
        result = dict()
        for col in pileup_columns:
            for pileup_read in col.pileups:
                target_tag_dict = dict()
                read_name = pileup_read.alignment.query_name
                if read_name not in result:
                    for tag in tag_names:
                        target_tag_dict[tag] = pileup_read.alignment.get_tag(tag)
                    # 不关心query_name是否重复，以最后一个为准
                    result[read_name] = target_tag_dict
        return result

    def _get_mutation_type(self, mutation:VariantRecord):
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

    def get_mut_support_reads(self, mut: VariantRecord, tag_names:tuple=tuple(), logger=None):
        mut_type= self._get_mutation_type(mut)
        if mut_type == 'snp':
            supports, tag_dict, support_seqs = self.get_snp_support_reads(mut.contig, mut.start, mut.ref, mut.alts[0], tag_names=tag_names, logger=logger)
        elif mut_type == 'substitution':
            # supports = self.get_substitution_support_reads(mut.contig, mut.start, mut.alts[0])
            supports, tag_dict, support_seqs = self.get_complex_support_reads(mut.contig, mut.start, mut.ref, mut.alts[0], tag_names=tag_names, logger=logger)
        elif mut_type == 'deletion':
            del_size = len(mut.ref) - len(mut.alts[0])
            supports, tag_dict, support_seqs = self.get_del_support_reads(mut.contig, mut.start, del_size, tag_names=tag_names, logger=logger)
        elif mut_type == 'insertion':
            supports, tag_dict, support_seqs = self.get_insert_support_reads(mut.contig, mut.start, mut.alts[0], tag_names=tag_names, logger=logger)
        elif mut_type == 'complex':
            supports, tag_dict, support_seqs = self.get_complex_support_reads(mut.contig, mut.start, mut.ref, mut.alts[0], tag_names=tag_names, logger=logger)
        else:
            raise Exception(f"{mut_type} is unexpected")
        return supports, tag_dict, support_seqs

    def check_support_consistency(self, read_seq_lst, max_disagree=2):
        # 检查支持同一个突变的reads的一致性
        # 假设这些证据都是指向同一个突变的可靠证据，那么这些reads在突变周围的碱基一致性应该非常高
        # 当某个突变附近存在一个杂合的突变，而当前的突变是纯合突变时，那么理论上支持纯合突变的reads的包含了附近的杂合突变
        # 这些reads在附近杂合突变位置的一致性则只有50%
        dis_agree = 0
        total = len(read_seq_lst)
        if total == 1:
            return True
        for bases in zip(*read_seq_lst):
            # 确保碱基大小写一致
            bases = [x.upper() for x in bases]
            freq = Counter(bases).most_common(1)[0][1]
            if freq/total < 0.75:
                dis_agree += 1
            if dis_agree >= max_disagree:
                return False
        return dis_agree <= max_disagree

    def detect_snp(self, contig, start, min_bq=13, ignore_overlaps=False):
        cols = self.bam.pileup(
            contig, start, start + 1,
            stepper='samtools',
            truncate=True,
            min_base_quality=min_bq,
            ignore_orphans=False,
            ignore_overlaps=ignore_overlaps,  # set 为True则意味着取质量高的base作为代表
            # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
            max_depth=300000,
            fastafile=None,
            # flag_filter=4,
        )
        alt, freq, depth = '', 0, 0
        for col in cols:
            bases = col.get_query_sequences()
            bases = [x.upper() for x in bases]
            base_dict = Counter(bases)
            depth = col.get_num_aligned()
            alt, freq = base_dict.most_common(1)
        return alt, freq, depth


class FindComplexVariant(ValidateMutationByBam):
    def __init__(self, vcf_file, bam_file, genome_file, tumor_sample, out_prefix=None):
        super().__init__(bam_file, genome_file)
        self.vcf_path = vcf_file
        if out_prefix is None:
            if self.vcf_path.endswith('vcf.gz'):
                out_prefix = self.vcf_path[:-6]
            else:
                out_prefix = self.vcf_path[:-4]
        self.out_prefix = out_prefix
        self.logger = set_logger(f'{self.out_prefix}.log', 'find_complex_variants')
        self.with_added_complex_vcf = self.out_prefix + '.AddComplex.vcf'
        self.with_only_complex_vcf = self.out_prefix + '.OnlyComplex.vcf'

    def read_vcf(self, only_pass=True, new_info:dict=None):
        with VariantFile(self.vcf_path) as temp:
            vcf_header = temp.header
            if new_info:
                vcf_header.info.add(
                    id = new_info['name'],
                    number=new_info['number'],
                    type=new_info['type'],
                    description=new_info['desc']
                )
            for r in temp:
                if r.alts is not None:
                    if only_pass:
                        if list(r.filter)[0] == "PASS":
                            yield r
                    else:
                        if r.alts is not None:
                            yield r

    def find_complex_variant(self, only_pass=True):
        """
        1. 要求vcf是按照坐标排序好的，经过normalization
        2. 检查相邻的SNP突变的距离是否在2个碱基以内，如果是，将尝试合并突变
        3. 检查非snp突变的距离是否在50bp以内，如果是，将尝试合并突变
        """
        with VariantFile(self.vcf_path) as fr:
            header = fr.header
            if 'MergeFrom' not in header.info:
                header.info.add('MergeFrom', number=1, type='String', description='Merged Variant Support Number | VariantPos | VariantPos')
            samples = list(header.samples)
            # 假设最后一个样本是目标样本
            target_idx = len(samples) - 1

        variants = read_vcf(self.vcf_path)

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
                dis = self.get_mut_distance(mut_a, mut_b)
                mut_a_type = self._get_mutation_type(mut_a)
                mut_b_type = self._get_mutation_type(mut_b)
                mut_a_is_s = mut_a_type in ['snp', 'substitution']
                mut_b_is_s = mut_b_type in ['snp', 'substitution']
                a_alt_depth = self.get_depth(mut_a, self.tumor)
                b_alt_depth = self.get_depth(mut_b, self.tumor)
                # depth_ratio = a_alt_depth / b_alt_depth
                if dis <= 2 and mut_a_is_s and mut_b_is_s:
                    logger.info(f">>>Try to merge mutations at {mut_a.contig}:{mut_a.pos}|{mut_b.pos}")
                    logger.info(mut_a.__str__())
                    logger.info(mut_b.__str__())
                    if mut_a_supports is None:
                        mut_a_supports, *_ = self.get_mut_support_reads(mut_a)
                        mut_b_supports, *_ = self.get_mut_support_reads(mut_b)
                    intersection = mut_a_supports & mut_b_supports
                    if len(intersection) >= 3:
                        # 这里仅仅合并距离在2bp以内的SNP或等长替换(substitution)
                        ref = self.genome.fetch(mut_a.contig, mut_a.start, mut_b.stop)
                        in_low_complex_region = not ref.isupper()
                        ref = ref.upper()
                        alt = mut_a.alts[0] + self.genome.fetch(mut_a.contig, mut_a.stop, mut_b.start) + mut_b.alts[0]
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
                        logger.info(f"--Both variant supporting read number:{len(intersection)}, they are {intersection}")

                        # 搜索支持merge突变的reads
                        mut_c_supports, *_ = self.get_mut_support_reads(record)
                        logger.info(f'--Merged variant supporting read number:{len(mut_c_supports)}, they are {mut_c_supports}')
                        if len(mut_c_supports) < 2 and len(alt) >= 3:
                            logger.info('支持合并后的突变的reads少于2个，说明2个snp中间也是snp，但原vcf中没有')
                            mid_alt, mid_req, mid_depth = self.detect_snp(mut_a.contig, mut_a.start + 1)
                            record.alts = [alt[0] + mid_alt + alt[2]]
                            logger.info(f'现在根据bam将中间的snp报告出来:{alt[1]}>{mid_alt}, {mid_req}/{mid_depth}')
                            mut_c_supports = self.get_mut_support_reads(record)
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
                    ref = self.genome.fetch(mut_a.contig, mut_a.start, mut_b.stop)
                    in_low_complex_region = not ref.isupper()
                    ref = ref.upper()
                    if dis > 0:
                        alt = mut_a.alts[0] + self.genome.fetch(mut_a.contig, mut_a.stop, mut_b.start) + mut_b.alts[0]
                        alt = alt.upper()
                    else:
                        logger.info('The above two deletion variants are part overlapped')
                        alt = mut_a.alts[0]
                    new_pos, new_ref, new_alt, alignment = global_pair_seq_align(mut_a.pos, ref, alt)
                    logger.info(alignment.__str__())
                    if new_pos is not None:
                        if mut_a_supports is None:
                            mut_a_supports, *_ = self.get_mut_support_reads(mut_a)
                        mut_b_supports, *_ = self.get_mut_support_reads(mut_b)
                        intersection = mut_a_supports & mut_b_supports
                        record = mut_a.copy() if a_alt_depth < b_alt_depth else mut_b.copy()
                        record.pos = new_pos
                        # print(record.ref, pair_align)
                        record.ref = new_ref
                        record.alts = [new_alt]
                        # 搜索支持merge突变的reads
                        mut_c_supports, *_ = self.get_mut_support_reads(record)
                        if 'MergeFrom' in mut_a.info:
                            logger.info('Super merging found!')
                            record.info['MergeFrom'] = mut_a.info['MergeFrom'] + '|' + str(mut_b.pos)
                        else:
                            record.info['MergeFrom'] = '|'.join(
                                [str(len(mut_c_supports)), str(mut_a.pos), str(mut_b.pos)])

                        if in_low_complex_region:
                            logger.info("warn: this mutation is likely in repeating region")
                        logger.info(f"--Vcf reported supporting read number:{a_alt_depth} vs {b_alt_depth}")
                        logger.info(f"--The program saw supporting read number:{len(mut_a_supports)} vs {len(mut_b_supports)}")
                        logger.info(f"--Both variant supporting read number:{len(intersection)}, they are {intersection}")
                        logger.info(f'--Merged variant supporting read number:{len(mut_c_supports)}, they are {mut_c_supports}')
                        if (len(intersection) >= 3 and len(mut_c_supports) >= 1) or (len(mut_c_supports) >= 20):
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
            self.bam.close()

        return out, out2


class VcfFilter(ValidateMutationByBam):
    def __init__(self, vcf_path, bam_file, genome_file, tumor=None, normal=None, normal_vcf=None, gene_primer_used=False):
        super().__init__(bam_file, genome_file)
        self.vcf = VariantFile(vcf_path)
        self.vcf_path = vcf_path
        self.tumor = tumor
        self.normal = normal
        self.gene_primer_used = gene_primer_used
        samples = list(self.vcf.header.samples)
        self.use_depth_bias = False
        if len(samples) == 1:
            self.tumor = samples[0]
            self.use_depth_bias = False
        elif len(samples) == 2:
            if self.tumor is None:
                self.tumor = samples[1]
                self.normal = samples[0]
                print(f'assume tumor and normal samples are {self.tumor} and {self.normal}')
            else:
                self.normal = list(set(samples) - {tumor})[0]
        else:
            if tumor is None:
                print('More than 2 samples found in the vcf, we consider the last sample as target sample')
                self.tumor = samples[-1]

        if normal_vcf:
            # 如果单独提供normal vcf时做如下处理
            normal_af_dict = dict()
            ctrl_vcf = VariantFile(normal_vcf)
            self.normal = list(ctrl_vcf.header.samples)[0]
            for r in ctrl_vcf:
                # if list(r.filter)[0] != "PASS":
                #     continue
                af = self.get_af_value(r, self.normal)
                dp = self.get_depth(r, self.normal)
                # 对照样本包含的突变可能是germline突变
                # 对照样本也可能包含假阳性突变、造血干细胞亚克隆突变，但突变会比较低，一般会低于1%
                # 因此如果要用对照样本过滤假阳性突变和造血性突变，需要对测序深度有一定要求
                if dp >= 5 and af >= 0.25:
                    key = (r.contig, r.pos, r.ref, r.alts[0])
                    # 如果对照样本测序深度不够，germline的突变AF可能被低估，而肿瘤样本由于测序深度足够，可以检测到远高于对照样本的AF
                    # 此时如果对照样本测到的AF作为error_rate, 肿瘤样本有可能通过，所以直接把germline突变AF提高到1
                    normal_af_dict[key] = 1.0
                elif dp >= 300:
                    # 对于测序深度足够，但频率没有满足germline标准的突变，可以作为测序错误率对肿瘤样本进行过滤
                    key = (r.contig, r.pos, r.ref, r.alts[0])
                    normal_af_dict[key] = af
                else:
                    continue
            ctrl_vcf.close()
            self.normal_af = normal_af_dict
        else:
            self.normal_af = dict()

    def add_contig_header(self, ref_dict):
        contig_info = []
        with open(ref_dict) as f:
            _ = f.readline()
            for line in f:
                lst = line.strip().split()
                name = lst[1].split(':')[1]
                length = lst[2].split(':')[1]
                md5 = lst[3].split(':')[1]
                ref_file = lst[4].split(':')[2]
                contig_info.append((name, length))
        # update header
        self.vcf.header.add_line(f'##reference={ref_file}')
        for line in contig_info:
            self.vcf.header.contigs.add(line[0], length=line[1])

    def poll_error_binomial_conf(self, error_rate, depth, confidence=0.95):
        lower, upper = stats.binom.interval(confidence, n=depth, p=error_rate)
        lower, upper = lower/depth, upper/depth
        if lower < 1e-6:
            lower = 1e-6
        return lower, upper

    def poll_error_norm_conf(self, error_rate, depth, confidence=0.95):
        # 可以根据样本量估算公式反推已知测序错误率和测序深度的条件下，计算测序错误率的上限，作为检测下限
        z = stats.norm.interval(confidence)[1]
        e = z/(depth/(error_rate*(1-error_rate)))**0.5
        lower = 0 if (error_rate - e <= 0) else (error_rate - e)
        upper = error_rate + e
        return lower, upper

    def estimate_min_required_depth(self, error_rate, af, confidence):
        """
        设LOD为检测下限，假设某个突变的真实频率为L，为减少假阳性，在某次实验中，我们应该要求“AF最低估计值”大于“测序错误最大估计值”，即：
        即要求：E由样本量估算公式反推 N = P*(1-P)*(Z/E)**2  用正太分布近似二型分布
            => L - E‘ > P + E
            =>  L - Z/(N/L/(1-L))**0.5 > P + Z/(N/P/(1-P))**0.5
            => (L-P) > Z/(N/P/(1-P))**0.5 + Z/(N/L/(1-L))**0.5
            => N > ( ( Z/(1/P/(1-P))**0.5 + Z/(1/L/(1-L))**0.5 ) / (L-P) )**2
        """
        P = error_rate
        L = af
        Z = stats.norm.interval(confidence)[1]
        if L - P > 0.00001:
            N = ((Z / (1 / P / (1 - P)) ** 0.5 + Z / (1 / L / (1 - L)) ** 0.5) / (L - P)) ** 2
        else:
            N = 1000000
        return int(N)

    def get_alt_binomial_pvalue(self, alt_depth, depth, error_rate):
        # 假设测序错误服从二型分布，可以计算alt_depth全部来自错误的概率
        return 1 - stats.binom.cdf(k=alt_depth, n=depth, p=error_rate)

    def qual_to_error_rate(self, base_qual):
        # 碱基质量值还原为概率值
        # Phred = -10 * log(error_p)
        error_prob = 10**(base_qual*(-0.1))
        return error_prob

    def get_mut_distance(self, mut_a, mut_b):
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
        mut_a_type = self.get_mutation_type(mut_a)
        mut_b_type = self.get_mutation_type(mut_b)
        if (mut_b.start < mut_a.stop < mut_b.stop) and (mut_a_type == mut_b_type == 'Deletion'):
            # 针对2个同时为del的突变进行处理
            # 两个突变区域首尾部分重叠
            # 如果重叠区域分别对应的突变序列也一致，则考虑合并的可能
            overlap_len = mut_b.start - mut_a.stop
            if mut_a.alts[0][-overlap_len:] == mut_b.alts[0][:overlap_len]:
                return mut_b.start - mut_a.stop
            else:
                return 10000
        if mut_b.start < mut_a.stop:
            # 这种情况常见于MSI
            print(f'忽略特殊突变：后一个突变{mut_b.pos} 在前一个突变的坐标范围内{mut_a.pos}-{mut_a.stop}')
            return 10000
        dis = mut_b.pos - mut_a.pos
        return dis

    def get_af_value(self, record, sample):
        if 'HIAF' in record.info and 'AF' in record.info:
            # vardict style
            af = max([record.info['HIAF'], record.info['AF'][0]])
        elif 'AF' in record.samples[sample]:
            # mutect2中的输出结果符合
            af = record.samples[sample]['AF'][0]
        elif 'FREQ' in record.samples[sample]:
            # for varscan2 style
            af = record.samples[sample]['FREQ']
            if '%' in af:
                af = float(af.strip('%')) * 0.01
        else:
            raise Exception('No AF or FREQ field found !')
        return af

    def get_normal_af(self, record):
        af = None
        if self.normal_af:
            # 对照信息来源独立的vcf
            mut_id = (record.contig, record.pos, record.ref, record.alts[0])
            if mut_id in self.normal_af:
                af = self.normal_af[mut_id]
        else:
            # 对照信息来源于paired模式出来的vcf，即和肿瘤样本的信息在同一个vcf
            # 在vardict的vcf中，DP是一个数，为总的测序深度，AD是一个数组，第一个是ref的depth，剩下的是alt的depth
            # 在mutect2的vcf中，DP是一个数，而且是一个经过过滤后的测序深度，AD和vardict中的AD含义一样
            normal_dp = record.samples[self.normal]['DP']
            if type(normal_dp) == tuple:
                normal_dp = sum(normal_dp)
            normal_ad = record.samples[self.normal]['AD']
            if type(normal_ad) == tuple:
                normal_ad = normal_ad[1]
            if normal_ad <= 0:
                normal_af = 0
            else:
                normal_af = normal_ad / normal_dp
            if normal_dp >= 5 and normal_af >= 0.25:
                # P=0.001; L=0.5; Z=1.96时，N > 4.3
                af = 1.0
            elif normal_dp >= 15 and normal_af >= 0.07:
                # 1/15 = 0.0666
                af = 1.0
            elif normal_dp >= 300 or normal_ad >= 3:
                # 无法判定为germline的突变作为背景噪音，用于后续过滤
                af = normal_af
            else:
                # 其他 (depth < 5) or (depth > 300 and af < 0.25)
                # 认为统计意义不足，不能用作过滤依据
                af = None
        return af

    def get_depth(self, record, sample):
        dp = record.samples[sample]['DP']
        if type(dp) == list:
            dp = sum(dp)
        return dp

    def get_mutation_type(self, record):
        if record.alts is None:
            return None
        if len(record.ref) == len(record.alts[0]) == 1:
            return "SNV"
        if len(record.ref) == len(record.alts[0]) > 1:
            return "Complex"
        if len(record.ref) > len(record.alts[0]) and record.ref.startswith(record.alts[0]):
            return 'Deletion'
        elif len(record.ref) < len(record.alts[0]) and record.alts[0].startswith(record.ref):
            return 'Insertion'
        else:
            return 'Complex'

    def get_raw_error_rate(self, record, min_error_rate=1e-6, read_len=150):
        error_rate = min_error_rate
        format_info = record.samples[self.tumor]
        if 'QUAL' in format_info:
            error_rate = self.qual_to_error_rate(format_info['QUAL']) * 0.9
        elif 'QUAL' in record.info:
            error_rate = self.qual_to_error_rate(record.info['QUAL']) * 0.9
        elif 'MBQ' in record.info:
            # mutect2 style
            base_qual = record.info['MBQ']
            if type(base_qual) == tuple:
                base_qual = record.info['MBQ'][1]
            error_rate = self.qual_to_error_rate(base_qual) * 0.9
        elif 'NM' in record.info and type(record.info['NM']) == float:
            # vardict style, 'NM':"Mean mismatches in reads"
            # 根据平均错配数量估计测序错误率，考虑真实突变和比对错误的存在，这个错误率可能偏大
            # 保守考虑，将这个错误率缩小5倍
            error_rate = record.info['NM'] / read_len * 0.2
        # print(f'we use QUAL info to infer error rate for {record.contig}:{record.start}:{record.ref}>{record.alts[0]}:', error_rate)
        return min(error_rate, min_error_rate)

    def pass_seq_error(self, record, sample, error_rate:float=None, alpha=0.05, factor=1.0):
        if error_rate == 0:
            raise Exception("ErrorRate=0? "+record.__str__())
        dp = self.get_depth(record, sample)
        af = self.get_af_value(record, sample)
        # 估计error_rate的置信区间
        confidence = 1 - alpha
        # confidence越大，则 error_upper越大, af_lower越小，min_depth越大， 这意味着过滤条件越严格
        error_lower, error_upper = self.poll_error_binomial_conf(error_rate=error_rate, depth=dp, confidence=confidence)
        af_lower, af_upper = self.poll_error_binomial_conf(error_rate=af, depth=dp, confidence=confidence)
        # 根据二型分布估计突变完全来自背景噪音或测序错误的概率值
        pvalue = self.get_alt_binomial_pvalue(alt_depth=round(dp*af), depth=dp, error_rate=error_rate)
        # error rate不能为0，af不能为1，否则会报错
        if af <= 0.99:
            min_depth = self.estimate_min_required_depth(error_rate, af, confidence)
            min_depth = int(min_depth * 0.9)
        else:
            min_depth = 5
        # print(dp, r.qual, error_rate, lower, upper)
        # 测试发现pvalue<alpha时，af 不一定小于upper，说明这可能是两种过滤策略
        if (af > error_upper*factor) and (pvalue < alpha) and (dp > min_depth) and (af_lower >= error_rate):
            return True, error_rate, error_upper, af_lower, pvalue, min_depth
        else:
            return False, error_rate, error_upper, af_lower, pvalue, min_depth

    def pass_depth_bias(self, record, cutoff=0.05, info_field='', normal=None, tumor=None):
        """
        当对照样本和正常样本同一批次测序且测序深度相当时才建议使用
        tumor vs normal depth bias
                    tumour  normal
        ref_depth   2475    2269
        var_depth   28    14
        odds, pvalue = stats.fisher_exact([[2475, 28], [2269, 14]], alternative='less')
        :param record: pysam.variant.record
        :param cutoff: 0.05
        :param info_field:
        :param fmt_filed:
        :param sample:
        :return:
        """
        passed = True
        pvalue = None
        possible_info_fields = [info_field] + ['SPV']
        for field in possible_info_fields:
            if field in record.info:
                if type(record.info[field]) == float:
                    pvalue = record.info[field]
                    break

        def get_depth(sample, field='DP'):
            dp = record.samples[sample][field]
            if type(dp) != int:
                if field == 'AD':
                    dp = dp[1]
                else:
                    dp = sum(dp)
            return dp

        if (normal and tumor) and (pvalue is None):
            normal_dp = get_depth(normal, 'DP')
            normal_ad = get_depth(normal, 'AD')
            tumor_dp = get_depth(tumor, 'DP')
            tumor_ad = get_depth(tumor, 'AD')
            odds, pvalue = stats.fisher_exact(
                [[tumor_dp-tumor_ad, tumor_ad], [normal_dp-normal_ad, normal_ad]],
                alternative='less'
            )
        # judge
        if pvalue is not None:
            if pvalue >= cutoff:
                passed = False
        else:
            # print('tumor vs normal depth bias filter is not applied for no valid field is found !')
            pass
        return passed, pvalue

    def pass_strand_bias(self, record, info_field='', fmt_field='', cutoff=0.001, sample=None):
        passed = True
        pvalue = None
        possible_info_fields = [info_field] + ['SBF']
        for field in possible_info_fields:
            if field in record.info:
                if type(record.info[field]) == float:
                    pvalue = record.info[field]
                    break
        # find info in "format"
        possible_fmt_fields = [fmt_field] + ['SB', 'SBF','DP4']
        if (pvalue is None) and (sample is not None):
            ref_fwd, ref_bwd, alt_fwd, alt_bwd = 0, 0, 0, 0
            for field in possible_fmt_fields:
                if field in record.samples[sample]:
                    if field == 'SBF':
                        # vardict的配对模式输出该值
                        pvalue = record.samples[sample][field]
                        break
                    else:
                        if type(record.samples[sample][field]) == tuple:
                            ref_fwd, ref_bwd, alt_fwd, alt_bwd = record.samples[sample][field]
                            break
            if sum([ref_fwd, ref_bwd, alt_fwd, alt_bwd]) > 0:
                odds_ratio, pvalue = stats.fisher_exact([[ref_fwd, ref_bwd], [alt_fwd, alt_bwd]])

        # judge
        if pvalue is not None:
            if pvalue <= cutoff:
                passed = False
        else:
            print('strand bias filter is not applied for no valid field is found !')
        # done
        return passed, pvalue

    def pass_pstd(self, record, cutoff=1e-4, gene_primer_used=False):
        """
        vardict中PSTD表示突变在reads中位置的标准差，如果突变在reads中出现的位置一样，那么该值为0，假阳的可能性极高
        通常，alt reads序列完全一致会导致这种情况，这意味着支持该突变的uniq read只有1条
        对于QIA的数据预处理，由于会统一去除primer序列，这使得最终同一个primer捕获的插入片段的测序结果中read1的起始位置都一致
        对于测通的情况，read1和read2包含的序列信息还会完全一致
        所以针对上述过滤可能不大合适，因此改为仅在支持的reads数目<=2的情况使用该过滤条件。
        """
        passed = True
        pstd = 0
        if gene_primer_used:
            max_alt_depth = 2
        else:
            max_alt_depth = 1000
        if 'PSTD' in record.info and 'VD' in record.info:
            pstd = record.info['PSTD']
            support_reads = record.info['VD']
            if pstd <= cutoff and support_reads <= max_alt_depth:
                passed = False
        return passed, pstd

    def pass_population_af(self, record, cutoff=0.01):
        # 要求是VEP注释, 且使用了--max_af参数
        # Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD
        # csq_dict['SYMBOL'],
        # csq_dict['Gene'],
        # csq_dict['CANONICAL'],
        # csq_dict['HGVSc'].replace('%3D', '='),
        # csq_dict['HGVSp'].replace('%3D', '='),
        # csq_dict['VARIANT_CLASS'],
        # csq_dict['IMPACT'],
        # csq_dict['MAX_AF'],
        # csq_dict['MAX_AF_POPS']
        csq_format = self.vcf.header.info['CSQ'].description.split('Format: ')[1]
        pop_af = []
        for each in record.info['CSQ']:
            csq_dict = dict(zip(csq_format.split('|'), each.split('|')))
            if csq_dict['MAX_AF']:
                pop_af.append(float(csq_dict['MAX_AF']))
        if pop_af:
            max_pop_af = max(pop_af)
            passed = max_pop_af <= cutoff
        else:
            passed = True
            max_pop_af = None
        return passed, max_pop_af

    def pass_msi_filter(self, record):
        # vardict报出的MSI长度有时候偏短，不能作为唯一依据
        # 考虑聚合酶滑脱现象： 一个或几个重复单元的缺失或者插入，一般占比小于15%
        # 如果是经过UMI矫正的read，其实这种聚合酶滑脱导致的错误理论上应该可以矫正回来，因此该过滤需要小心对待
        # [SNV]      : AAAAGA -> AAAAAA;
        # [Insertion]: AAAAA -> AAAA[A]A;
        # [Deletion] : AAAAGAA -> AAAAAA
        passed = True
        af = self.get_af_value(record, self.tumor)
        # vardict的info中提供了突变侧翼碱基序列信息
        if 'LSEQ' in record.info and passed:
            left = record.info['LSEQ']
            right = record.info['RSEQ']
            mut_type = record.info['TYPE']
            if mut_type == 'SNV':
                # AAAAGA -> AAAAAA
                # 上述是非典型MSI错误，仅考虑单碱基串联重复
                alt = record.alts[0]
                alt_seq = left[:-5] + alt + right[:5]
                if re.search(alt + '{6,}', alt_seq):
                    # 找到不少于6次单碱基重复
                    if af < 0.05:
                        passed = False
            elif mut_type == 'Deletion':
                # CGG[CG]GGGGGA - > CGG[C]GGGGGA
                # 情形1：认为删除的碱基是重复单元, 那么考虑参考序列是否为重复串联单元
                alt = record.ref[1:]
                if len(alt) >= 3:
                    min_repeat = 4
                else:
                    min_repeat = 5
                if len(alt) * (min_repeat-1) <= 20:
                    ref_seq = left[:-len(alt)*(min_repeat-1)] + record.ref + right[:len(alt)*(min_repeat-1)]
                    if re.search(f'({alt})'+'{' + str(min_repeat) + ',}', ref_seq):
                        # 参考序列中可以找到以删除序列为基本单元的串联重复
                        if af < 0.05:
                            passed = False
                if passed:
                    # 再次检查，考虑删除的序列是重复单元的整数倍
                    if len(alt) > 1:
                        msi_len = int(record.info['MSILEN'])
                        if len(alt) % msi_len == 0:
                            alt = alt[:msi_len]
                            if len(alt) >= 3:
                                min_repeat = 4
                            else:
                                min_repeat = 5
                            if len(alt) * (min_repeat - 1) <= 20:
                                ref_seq = left[:-len(alt) * (min_repeat-1)] + record.ref + right[:len(alt) * (min_repeat-1)]
                                if re.search(f'({alt})' + '{' + str(min_repeat) + ',}', ref_seq):
                                    # 参考序列中可以找到以删除序列为基本单元的串联重复
                                    if af < 0.05:
                                        passed = False
                # 情形2：考虑最多删除2bp后，形成单碱基串联重复，如AAAAGGAA -> AAAAAA
                if passed:
                    alt = record.ref[0]
                    alt_seq = left[:-4] + alt + right[:4]
                    if len(record.ref) <= 3 and re.search(alt + '{5,}', alt_seq):
                        if af < 0.05:
                            passed = False
            elif mut_type == 'Insertion':
                # AAAAA -> AAAA[A]A; ATATATAT -> AT[AT]ATATAT
                # 考虑单碱基或多碱基串联重复, 重复单元的长度不超过4
                alt = record.alts[0][1:]
                if len(alt) > 1:
                    msi_len = int(record.info['MSILEN'])
                    if len(alt) % msi_len == 0:
                        alt = alt[:msi_len]
                if len(alt)*5 <= 20:
                    # vardict提供的侧翼长度信息为20
                    alt_seq = left[:-len(alt)*5] + record.alts[0] + right[:len(alt)*5]
                    if re.search(f'({alt})'+'{6,}', alt_seq):
                        if af < 0.05:
                            passed = False
            elif mut_type == 'Complex':
                # LSEQ=TTATATATATATATATATAT;RSEQ=TTTTTTTTTTCTGCAGCTGC；ATA>TTT
                if int(record.info['MSI']) >= 5:
                    # 假设alt是重复单元
                    alt = record.alts[0]
                    if len(alt) * 3 <= 20:
                        # vardict提供的侧翼长度信息为20
                        alt_seq = left[:-len(alt) * 3] + alt + right[:len(alt) * 3]
                        if re.search(f'({alt})' + '{4,}', alt_seq):
                            if af < 0.05:
                                passed = False

        return passed

    def format_txt_output(self, r) -> dict:
        csq_format = self.vcf.header.info['CSQ'].description.split('Format: ')[1]
        picked = None
        canonical = None
        for each in r.info['CSQ']:
            csq_dict = dict(zip(csq_format.split('|'), each.split('|')))
            # 找到被flag为1的记录作为报告
            if 'PICK' in csq_dict and csq_dict['PICK'] == "1":
                picked = csq_dict
            # 找到canonical转录本，且
            if csq_dict['CANONICAL'] == 'YES':
                if 'SOURCE' in csq_dict:
                    if csq_dict['SOURCE'] == 'RefSeq':
                        canonical = csq_dict
                else:
                    canonical = csq_dict
            if picked and canonical:
                break
        # 优先选择RefSeq经典转录本对应的注释，其次选择被PICK为1的注释
        csq_dict = canonical or picked
        # 转换坐标或格式
        mutation_type = self.get_mutation_type(r)
        if mutation_type == 'SNP':
            start_pos = r.pos
            end_pos = r.pos
            ref = r.ref
            alt = csq_dict['Allele']
        elif mutation_type == 'Insertion':
            start_pos = r.pos + 1
            end_pos = '-'
            ref = '-'
            alt = csq_dict['Allele']
        elif mutation_type == 'Deletion':
            start_pos = r.pos + 1
            ref = r.ref[1:]
            end_pos = start_pos + len(ref) - 1
            alt = csq_dict['Allele']
        else:
            # Complex
            start_pos = r.pos
            ref = r.ref
            end_pos = start_pos + len(ref) - 1
            alt = csq_dict['Allele']
        cHGVS = csq_dict['HGVSc'].replace('%3D', '=')
        if cHGVS:
            cHGVS = cHGVS.split(':')[1]
        pHGVS = csq_dict['HGVSp'].replace('%3D', '=')
        if pHGVS:
            pHGVS = pHGVS.split(':')[1]

        target_info = {
            'Sample': self.tumor,
            "Chr": r.contig,
            "Start": start_pos,
            "End": end_pos,
            "Ref": ref,
            "Alt": alt,
            "Gene": csq_dict['SYMBOL'],  # Gene
            "Type": mutation_type,
            "Transcript": csq_dict['Feature'],
            "cHGVS": cHGVS,  # C-Dot Notation
            "pHGVS": pHGVS,  # P-Dot Notation
            "VAF(%)": f'{self.get_af_value(r, self.tumor)*100:.2f}',  # Allele Frequency
            # additional information
            '=HYPERLINK("https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html","Consequence")': csq_dict['Consequence'],  # Consequence(s)
            "CLIN_SIG": csq_dict['CLIN_SIG'],
            "Exon": csq_dict['EXON'],  # Affected Exon(s)
            "Strand": csq_dict['STRAND'],
            "VariantClass": csq_dict['VARIANT_CLASS'],
            "Depth": self.get_depth(r, self.tumor),  # Depth
            "ExistingVariation": csq_dict['Existing_variation'],
            "MAX_AF": csq_dict['MAX_AF'],
            "MAX_AF_POPS": csq_dict['MAX_AF_POPS'],
            "gnomADe_EAS_AF": csq_dict['gnomADe_EAS_AF'],
            "CANONICAL": csq_dict['CANONICAL'],
            "HGVS_OFFSET": 0,
            'CtrlSample': self.normal,
        }
        if 'LSEQ' in r.info:
            target_info['LSEQ'] = r.info['LSEQ']  # from vardict output
            target_info['REF'] = ref
            target_info['RSEQ'] = r.info['RSEQ']

        if 'LOD' in r.info:
            target_info['LOD(error_rate,error_upper,af_lower,pvalue,min_depth)'] = r.info['LOD']
        if 'ConsInfo' in r.info:
            # consensus depth information
            target_info['(Alt_cD1, Alt_cD2+, cD1, cD2+, MeanError)'] = r.info['ConsInfo']
        # 根据注释结果判断是否报告，增加报告字段
        consequences = set(csq_dict['Consequence'].split('&'))
        if (consequences - {
            'intron_variant',
            '5_prime_UTR_variant',
            '3_prime_UTR_variant',
            'upstream_gene_variant',
            'downstream_gene_variant',
            'synonymous_variant',
            '5_prime_UTR_variant',
            '3_prime_UTR_variant',
            'splice_donor_5th_base_variant',
            'splice_region_variant',
            'splice_donor_region_variant',
            'splice_polypyrimidine_tract_variant',
            'incomplete_terminal_codon_variant',
            'start_retained_variant',
            'stop_retained_variant',
        }):
            target_info['SelectedByCsq'] = 'Yes'
        else:
            target_info['SelectedByCsq'] = 'No'

        # 还有突变坐标的潜在问题需要解决
        # https://hgvs.readthedocs.io/en/stable/examples/manuscript-example.html
        # VEP现在默认会给出符合HGVS右对齐的注释信息，但是不会更改突变的基因组上坐标信息
        # 因此正链上的基因的突变位置信息需要修改于hgvs保持一致，这是室间质评的要求。
        # ----解决问题方案一-----：
        # import hgvs.variantmapper
        # vm = hgvs.assemblymapper.AssemblyMapper(
        #     hdp, assembly_name='GRCh37', alt_aln_method='splign')
        # vm.c_to_g(var_c1) 输出如下:
        # SequenceVariant(ac=NC_000001.10, type = g, posedit = 150550916; G > A)
        # 从posedit中获取坐标信息，希望是正确的
        # https://hgvs.readthedocs.io/en/stable/installation.html
        # ---解决问题方案二---：其实可以根据VEP的注释结果直接提取，如下说明，VEP其实提供了HGVS_OFFSET值
        # http://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_shift_hgvs
        # Enable or disable 3' shifting of HGVS notations.
        # When enabled, this causes ambiguous insertions or deletions
        # (typically in repetetive sequence tracts) to be "shifted" to their most 3' possible coordinates
        # (relative to the transcript sequence and strand) before the HGVS notations are calculated;
        # the flag HGVS_OFFSET is set to the number of bases by which the variant has shifted,
        # relative to the input genomic coordinates.
        # If HGVS_OFFSET is equals to 0, no value will be added to HGVS_OFFSET column.
        # Disabling retains the original input coordinates of the variant. Default: 1 (shift)
        if 'HGVS_OFFSET' in csq_dict and csq_dict['HGVS_OFFSET']:
            target_info['HGVS_OFFSET'] = int(csq_dict['HGVS_OFFSET'])
            if mutation_type == 'Insertion':
                # 室间质评要求对于insertion，start是插入序列前一个碱基位置
                target_info['Start'] = start_pos + int(csq_dict['HGVS_OFFSET']) - 1
                if end_pos != '-':
                    target_info['End'] = start_pos + int(csq_dict['HGVS_OFFSET'])
            else:
                target_info['Start'] = start_pos + int(csq_dict['HGVS_OFFSET'])
                target_info['End'] = end_pos + int(csq_dict['HGVS_OFFSET'])

        return target_info

    def write_out_txt(self, vcf, out):
        result = []
        with VariantFile(vcf) as fr:
            for r in fr:
                result.append(self.format_txt_output(r))
        df = pd.DataFrame(result)
        df.to_csv(out, sep='\t', index=False)
        df.to_excel(out[:-4]+'.xlsx', index=False)

    def filtering(self, genome, bam=None, ref_dict=None, out_prefix=None, min_error_rate=None,
                  error_rate_file=None, min_af=0.0001, min_depth=5, alpha=0.05, disable_bg_model=False):
        # 先给vcf的header添加新字段定义才能往添加新的字段信息
        self.vcf.header.info.add(
            'LOD', number=5, type='Float',
            description='[error_rate, error_upper, af_lower, pvalue, min_depth]. '
                        'The first value is input error rate which will be used as theoretical frequency to '
                        f'calculate the second value. The second value is the upper bound of the {alpha} confidence interval of error rate. '
                        'The fourth value is the probability of alt observed from background noise'
        )
        self.vcf.header.info.add(
            'ConsInfo', number=5, type='Float', description='consensus status of reads'
        )
        self.vcf.header.filters.add('NoiseFromNormal', number=None, type=None, description='Noise from normal sample')
        self.vcf.header.filters.add('BackgroundNoise', number=None, type=None, description='Noise from background')
        self.vcf.header.filters.add('FromGermline', number=None, type=None, description='Considered as germline variant found in normal vcf')
        self.vcf.header.filters.add('HighPopFreq', number=None, type=None, description='Population frequency')
        # self.vcf.header.filters.add('LowAF', number=None, type=None, description=f'AF smaller than {min_af}')
        self.vcf.header.filters.add('MSIFilter', number=None, type=None, description=f'Likely PCR Slippage in MSI region')
        self.vcf.header.filters.add('LowUmiReadSupport', number=None, type=None, description=f'Low UMI read supporting')
        self.vcf.header.filters.add('BadSnpQual', number=None, type=None, description=f'Base qual of SNV is not good')
        self.vcf.header.filters.add('InconsistentSupports', number=None, type=None, description=f'Variant supporting reads are not concordant on other bases')
        if 'Bias' not in self.vcf.header.filters:
            self.vcf.header.filters.add('Bias', number=None, type=None, description="severe strand bias")
        if self.vcf.header.contigs.__len__() < 1:
            self.add_contig_header(ref_dict)

        if out_prefix is None:
            out_prefix = self.vcf_path.rsplit('.', 1)[0]
        out_vcf_name = out_prefix+'.final.vcf'
        vcf_out = VariantFile(out_vcf_name, "w", header=self.vcf.header.copy())
        vcf_discard = VariantFile(out_prefix+'.discarded.vcf', "w", header=self.vcf.header.copy())

        # 读取bam信息
        if bam:
            bamer = ValidateMutationByBam(bam_file=bam, genome_file=genome)
        else:
            bamer = None

        # 如果仅仅提供了一个全局的最低错误率信息
        key_left = key_right = 0
        if error_rate_file:
            seq_error_dict = json.load(open(error_rate_file))
            key_len = max(len(x) for x in seq_error_dict['A'].keys())//2
            key_left = key_right = key_len
        else:
            if min_error_rate:
                seq_error_dict = dict()
                for b in 'ATCGID':
                    for i in set('ATCG') - {b}:
                        seq_error_dict.setdefault(b, dict())[i] = {b: float(min_error_rate)}
                print(f'{seq_error_dict} is used as the background noise information')
            else:
                seq_error_dict = dict()

        # 读取参考基因组信息
        gn = FastaFile(genome)
        # print(seq_error_dict.keys())
        lod_list = []
        keep = 0
        total = 0
        filter_reasons = []
        log_file = open(out_prefix + '.filtering.log', 'w')
        vcf_loh = VariantFile(out_prefix+'.LOH.vcf', "w", header=self.vcf.header.copy())
        logger = set_logger(out_prefix+'.SNV.metrics.txt', 'var_metrics')
        for r in self.vcf:
            total += 1
            reasons = []
            if '<' in r.alts[0]:
                # 跳过vardict中输出的特殊突变
                print('skip special variant:', r.contig, r.pos, r.ref, list(r.alts), file=log_file)
                vcf_discard.write(r)
                continue

            if 'N' in r.alts[0]:
                print('skip "N" containing variant:', r.contig, r.pos, r.ref, list(r.alts), file=log_file)
                vcf_discard.write(r)
                continue

            if r.info['DP'] < min_depth:
                print(f'skip variant in shallow depth({r.info["DP"]}):', r.contig, r.pos, r.ref, list(r.alts), file=log_file)
                vcf_discard.write(r)
                continue

            if len(r.ref) > 51 or len(r.alts[0]) > 51:
                print(f'skip long mutation ({r.info["DP"]}):', r.contig, r.pos, r.ref, list(r.alts), file=log_file)
                vcf_discard.write(r)
                continue

            # 过滤VCF中原本没有被判定通过的突变
            if list(r.filter)[0] != "PASS":
                reasons = list(r.filter)

            # LOH
            af = self.get_af_value(r, self.tumor)
            if af <= 0:
                # print('discard AF=0 variant as following, maybe it called as LOH site', file=log_file)
                vcf_loh.write(r)
                continue

            # af cutoff
            if af < min_af:
                # reasons.append('LowAF')
                vcf_discard.write(r)
                continue

            # 4. 根据人群频率过滤
            judge4 = self.pass_population_af(r, cutoff=0.01)
            mutation_type = self.get_mutation_type(r)
            if not judge4[0]:
                # 如果人群频率没有通过，为加快速度，跳过后续过滤判定步骤
                reasons.append('HighPopFreq')
            else:
                # 1.根据测序错误率或germline突变频率过滤
                ctrl_af_as_error_rate = False
                # r.pos正好是1-based, r.start 是0-based
                error_rate = None
                if seq_error_dict:
                    if mutation_type == 'SNV':
                        key = gn.fetch(r.contig, r.start - key_left, r.start + 1 + key_right).upper()
                        # error_rate = seq_error_dict[key][r.alts[0]]
                        if key in seq_error_dict[r.alts[0]]:
                            error_rate = seq_error_dict[r.alts[0]][key][r.alts[0]]
                    elif mutation_type == 'Insertion':
                        key = gn.fetch(r.contig, r.start - key_left, r.start + 1 + key_right).upper()
                        if key in seq_error_dict['I']:
                            error_rate = seq_error_dict['I'][key]['I']
                            # if len(r.alts[0])-len(r.ref) >= 3:
                            #     error_rate = error_rate**2
                    elif mutation_type == 'Deletion':
                        # deletion
                        key = gn.fetch(r.contig, r.pos - key_left, r.pos + 1 + key_right).upper()
                        if 'D' in seq_error_dict['D']:
                            # error_rate = seq_error_dict[key]['']
                            if key in seq_error_dict['D']:
                                error_rate = seq_error_dict['D'][key]['D']
                                # if len(r.ref) - len(r.alts[0]) >= 3:
                                #     error_rate = error_rate ** 2
                        # print('del', key, error_rate, r.ref, list(r.alts))
                    elif mutation_type == 'Complex':
                        # 使用第一个碱基的信息，也即类似snp的方式处理
                        key = gn.fetch(r.contig, r.start - key_left, r.start + 1 + key_right).upper()
                        # print(r.pos, key, r.ref, r.alts)
                        if key in seq_error_dict[r.alts[0][0].upper()]:
                            error_rate = seq_error_dict[r.alts[0][0].upper()][key][r.alts[0][0].upper()]
                    else:
                        raise Exception(f'mutation type is unknown for {r.__str__()}')

                if self.normal:
                    # 当存在对照样本时，如果某个位点在对照样本也存在突变，且突变频率大于seq_error时，可以把对照样本中的突变频率作为测序错误率进行过滤
                    # 一开始认为这样可以过滤掉germline突变，如果对照样本测序深度足够，还可以过滤假阳性突变、克隆造血突变
                    # 但是如果对照样本测序深度不够，germline的突变AF可能被低估，而肿瘤样本由于测序深度足够，可以检测到远高于对照样本的AF
                    # 此时如果直接用对照样本的AF作为error_rate, 肿瘤样本有可能通过判定，所以最好的方式是直接扔掉被判定为germline的突变
                    normal_af = self.get_normal_af(r)
                    if normal_af and (normal_af > (error_rate or 0)):
                        ctrl_af_as_error_rate = True
                        error_rate = normal_af

                if (error_rate is None) or (error_rate < min_error_rate):
                    # 没有提取到错误率信息，尝试从record的字段信息提取error_rate的估计值
                    error_rate = self.get_raw_error_rate(r, read_len=150, min_error_rate=min_error_rate)

                # 1.seq error过滤或者germline突变过滤
                if error_rate > 0.999 or af > 0.999:
                    reasons.append('FromGermline')
                else:
                    if not disable_bg_model:
                        judge = self.pass_seq_error(r, self.tumor, error_rate, alpha=alpha)
                        if not judge[0]:
                            if ctrl_af_as_error_rate:
                                reasons.append('NoiseFromNormal')
                            else:
                                # print(key, error_rate, r.ref, list(r.alts))
                                reasons.append('BackgroundNoise')
                            pass
                        # error_rate, error_upper, af_lower, pvalue, min_depth
                        r.info['LOD'] = tuple(judge[1:])
                        lod_list.append(judge[2])

                # 2.根据strand bias进行过滤
                # 当使用链偏倚作为过滤指标时要注意：只有极端链偏倚的SNP才应被视为假阳性候选突变
                judge2 = self.pass_strand_bias(r, cutoff=1e-6, sample=self.tumor)
                if not judge2[0]:
                    if 'Bias' not in reasons:
                        reasons.append('Bias')
                        pass

                # 3. position std filtering, 对于使用gene_primer技术的文库，如果突变在reads中出现的位置基本不变,且支持的read<=2，需要过滤掉
                judge3 = self.pass_pstd(r, cutoff=0.00001, gene_primer_used=self.gene_primer_used)
                if not judge3[0]:
                    if 'pSTD' not in reasons:
                        reasons.append('pSTD')

                # 5. 根据MSI状况
                judge5 = self.pass_msi_filter(r)
                if not judge5:
                    reasons.append('MSIFilter')

                # 6. 获取突变支持信息
                if bamer and len(reasons) == 0:
                    # 该步骤比较耗时，为加快速度，仅仅对通过前面全部判定条件的突变进行分析
                    support_reads, tag_dict, support_seqs = bamer.get_mut_support_reads(r, tag_names=('cD', 'cE'), logger=logger)
                    if 'snp_is_good' in tag_dict:
                        snp_good = tag_dict.pop('snp_is_good')
                        if not snp_good:
                            reasons.append('BadSnpQual')
                    if support_reads:
                        fam1_support_num = sum(tag_dict[x]['cD'] == 1 for x in support_reads)
                        fam2_support_num = len(support_reads) - fam1_support_num
                        fam1_all_num = sum(tag_dict[x]['cD'] == 1 for x in tag_dict)
                        fam2_all_num = len(tag_dict) - fam1_all_num
                        consensus_error = [tag_dict[x]['cE'] for x in tag_dict]
                        # median_error = statistics.median(consensus_error)
                        mean_error = statistics.mean(consensus_error)
                        r.info['ConsInfo'] = (fam1_support_num, fam2_support_num, fam1_all_num, fam2_all_num, mean_error)
                        # 检查证据reads的一致性
                        consistency = bamer.check_support_consistency(support_seqs, max_disagree=2)
                        if not consistency:
                            reasons.append('InconsistentSupports')
                        if fam2_support_num == 0 and mutation_type == 'SNV':
                            if fam1_support_num < 5 or fam2_support_num < 2:
                                reasons.append('LowUmiReadSupport')
                        else:
                            # 针对较大的deletion或insertion,本程序找到的证据可能偏低，因为没有局部组装和重比对的功能
                            if mutation_type == 'SNV' and fam2_support_num < int(fam2_all_num * af * 0.9):
                                reasons.append('LowUmiReadSupport')
                        if not disable_bg_model:
                            # 根据umi统计得到的mean_error再过滤一次
                            if mean_error > 0:
                                judge = self.pass_seq_error(r, self.tumor, mean_error, alpha=alpha)
                                if not judge[0]:
                                    reasons.append('BackgroundNoise')

            # 更新filter的内容和输出结果
            if reasons:
                filter_reasons.append(tuple(reasons))
                for each in reasons:
                    r.filter.add(each)
                vcf_discard.write(r)
            else:
                keep += 1
                vcf_out.write(r)

        vcf_out.close()
        vcf_discard.close()
        vcf_loh.close()
        gn.close()
        bamer.bam.close()
        self.write_out_txt(out_vcf_name, out_prefix+'.final.txt')
        if lod_list:
            print('median LOD:', statistics.median(lod_list), file=log_file)
            print('min LOD:', min(lod_list), file=log_file)
            print('max LOD:', max(lod_list), file=log_file)
        print(f'discard {total - keep} variants while keep {keep} ones!', file=log_file)
        reason_couts = Counter(filter_reasons).most_common()
        for k, v in reason_couts:
            print(f'{v} variants are filtered out because of {k}', file=log_file)


def filter_vcf(vcf, genome, ref_dict=None, tumor_name=None, bam=None, bed=None, normal_vcf=None, alpha=0.01, min_af=0.0001,
              exclude_from=None, out_prefix=None, min_error_rate=None, error_rate_file=None, disable_bg_model=False, center_size:tuple=(1, 1)):
    # if bam and bed and (not error_rate_file):
    #     error_rate_file = estimate_context_seq_error(
    #         bed, bam, prefix=out_prefix, center_size=center_size,
    #         genome=genome, exclude_from=exclude_from
    #     )
    vcf = VcfFilter(vcf_path=vcf, bam_file=bam, genome_file=genome, tumor=tumor_name, normal_vcf=normal_vcf, gene_primer_used=False)
    vcf.filtering(
        genome=genome,
        bam=bam,
        ref_dict=ref_dict,
        out_prefix=out_prefix,
        min_error_rate=min_error_rate,
        error_rate_file=error_rate_file,
        alpha=alpha,
        min_af=min_af,
        disable_bg_model=disable_bg_model
    )


if __name__ == '__main__':
    # from xcmds import xcmds
    # xcmds.xcmds(locals())
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', type=Path, required=True, help='vcf file annotated with vep')
    parser.add_argument('-genome', type=Path, required=True, help='path to indexed genome fasta')
    parser.add_argument('-bam', type=Path, required=False, help='bam file which will be used to estimate background noise')
    parser.add_argument('-bed', type=Path,  required=False, help='target region file which will be used to estimate background noise')
    parser.add_argument('-center_size', type=int, nargs='+', default=(1, 1), help='extending size around ref base during background noise estimating')
    parser.add_argument('-exclude_from', type=Path, required=False, help='bed file containing known variant in input bam, these variants will be excluded during background noise estimating')
    parser.add_argument('-ref_dict', type=Path, required=False, help='path to genome dict file which will be used to add missed contig header in vcf')
    parser.add_argument('-tumor_name', required=False, help='tumor sample name in vcf')
    parser.add_argument('-normal_vcf', type=Path, required=False, help='normal sample vcf file')
    parser.add_argument('-error_rate_file', type=Path, required=False, help='Estimated background noise file, if not provided, bam file will be used to generate one')
    parser.add_argument('-min_error_rate', type=float, default=1e-6, help='global minimum error rate, if error rate cannot be aquired in other ways, this value will be used')
    parser.add_argument('-alpha', type=float, default=0.05, help='cutoff of pvalue from background noise model. The pvalue represents the probability of variants come from background noise')
    parser.add_argument('-min_af', type=float, default=0.0001, help='hard cutoff of AF')
    parser.add_argument('--disable_bg_model', default=False, action='store_true', help='disable background noise model filter')
    parser.add_argument('-out_prefix', help='output file prefix')
    args = parser.parse_args()
    filter_vcf(
        vcf=args.vcf,
        genome=args.genome,
        bam=args.bam,
        bed=args.bed,
        exclude_from=args.exclude_from,
        ref_dict=args.ref_dict,
        center_size=args.center_size,
        tumor_name=args.tumor_name,
        normal_vcf=args.normal_vcf,
        error_rate_file=args.error_rate_file,
        min_error_rate=args.min_error_rate,
        alpha=args.alpha,
        out_prefix=args.out_prefix,
        min_af=args.min_af,
        disable_bg_model=args.disable_bg_model
    )
