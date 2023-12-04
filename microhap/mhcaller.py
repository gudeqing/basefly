import os
import json
import statistics
import itertools
import json
import pandas as pd
from collections import Counter
from scipy import stats
from pysam import FastaFile, VariantFile, AlignmentFile, VariantHeader, VariantRecord

__author__ = 'gdq'


class MicroHapCaller(object):
    def __init__(self, bam_file, genome_file, micro_hap_file, out_prefix):
        self.bam = AlignmentFile(bam_file)
        self.genome = FastaFile(genome_file)
        self.markers = self.parse_micro_hap(micro_hap_file)
        self.run(out_prefix)

    @staticmethod
    def parse_micro_hap(def_file) -> dict:
        """
        :param def_file:
        --------------------------------------
        Marker  Offset  Chrom   OffsetHg38
        mh01CP-016      27      chr1    55559012
        mh01CP-016      65      chr1    55559050
        mh01CP-016      71      chr1    55559056
        mh01CP-010      38      chr1    85240117
        ----------------------------------------
        :return: {maker_name: {contig: chr1, snp_pos: [1, 2, 5]}}
        """
        result = dict()
        with open(def_file) as f:
            _header = f.readline().strip().split()
            for line in f:
                lst = line.strip().split()
                content = result.setdefault(lst[0], dict())
                content['contig'] = lst[2]
                # position should be 0-based
                content.setdefault('snp_pos', []).append(int(lst[3]))
        return result

    @staticmethod
    def qual_diff_test(alt_reads: set, ref_reads, qual_dict: dict, alternative='two-sided'):
        alt_quals = [qual_dict[r] for r in alt_reads if r in qual_dict]
        ref_quals = [qual_dict[r] for r in ref_reads if r in qual_dict]
        # 清洗
        alt_quals = [x for x in alt_quals if x is not None]
        ref_quals = [x for x in ref_quals if x is not None]
        # 比较
        # 优先进行mannwhitneyu检验，适合样本量较多的情况
        if len(alt_quals) >= 10 and len(ref_quals) >= 10:
            u, pvalue = stats.mannwhitneyu(alt_quals, ref_quals, alternative=alternative)
            # u, pvalue = stats.ranksums(alt_quals, ref_quals, alternative=alternative)
            # u, pvalue = stats.ttest_ind(alt_quals, ref_quals, alternative=alternative)
        elif len(alt_quals) >= 3 and len(ref_quals) >= 3:
            # 非配对独立t检验
            s, pvalue = stats.ttest_ind(alt_quals, ref_quals, alternative=alternative)
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
            return seqs
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
                FoxoG = round(F2R1 / (F1R2 + F2R1), 3)
            elif ref == "G" or ref == "T":
                FoxoG = round(F1R2 / (F1R2 + F2R1), 3)
        return FoxoG

    def marker_typing(self, marker, contig, positions, min_bq=10):
        result = dict()
        # result的key是position=(contig, start, ref_base)，value是相应position的碱基测序结果
        depths = []
        for pos in positions:
            # 期望这里的pos是0-based
            ref_base = self.genome.fetch(contig, pos, pos+1).upper()
            pileup_columns = self.bam.pileup(
                contig, pos, pos+1,
                stepper='samtools',
                truncate=True,
                min_base_quality=min_bq,
                ignore_orphans=False,
                # *  If called, mpileup will detect overlapping
                # *  read pairs and for each base pair set the base quality of the
                # *  lower-quality base to zero, thus effectively discarding it from
                # *  calling. If the two bases are identical, the quality of the other base
                # *  is increased to the sum of their qualities (capped at 200), otherwise
                # *  it is multiplied by 0.8.
                ignore_overlaps=False,  # 设置为True时，pileup 返回的很多碱基质量值为0，pysam的有issue提出这个问题，解决方案是ignore_overlaps=False
                # 这里的max_depth一定要设的足够大，否则有可能漏掉reads
                max_depth=300000,
                fastafile=self.genome,
                compute_baq=False,
                # flag_filter: The default is BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP.
            )
            support_reads = dict()
            # F=foward strand, R=reverse strand, 1:read1, 2:read2
            F2R1 = dict()
            F1R2 = dict()
            # 把指标放到metrics中
            metrics = dict()
            for idx, col in enumerate(pileup_columns):
                # 提取当前位点的碱基质量值信息，并构建字典，以read_name作为key
                base_quals = col.get_query_qualities()
                query_names = col.get_query_names()
                map_quals = col.get_mapping_qualities()
                predict_bases = col.get_query_sequences()
                # 找到出现频率最多的碱基, 在后续比较中作为参考
                main_base = Counter(predict_bases).most_common(1)[0][0].upper()
                # 存储目标位点的测序深度信息
                depth = col.get_num_aligned()
                depths.append(depth)
                # 注意，query_names有可能重复，因为ignore_overlaps=False
                base_qual_dict = dict(zip(query_names, base_quals))
                # 获取目标位点的信息，如比对质量信息和碱基质量信息
                map_qual_dict = dict(zip(query_names, map_quals))
                dist_to_read_end_dict = dict()
                # 获取支持变异的reads信息
                for pileup_read in col.pileups:
                    alignment = pileup_read.alignment
                    read_name = alignment.query_name
                    query_seq = alignment.query_sequence
                    query_pos = pileup_read.query_position
                    query_len = alignment.query_length
                    # 对于pair，只会记录其中一条记录
                    if read_name not in support_reads:
                        if query_pos is not None and pileup_read.indel == 0:
                            # 计算突变位点距离read两端的距离的最小值, 针对multiplexPCR得到的，这个值没有什么用处
                            dist_to_read_end_dict[read_name] = min(query_pos, query_len - query_pos)
                            # 提取目标位置检测到的具体碱基成分，其可能是和ref一致，也可能不一致
                            predict_base = query_seq[query_pos].upper()
                            # 注意，由于未知原因，这里找到的read_name不一定在前面的query_names中
                            # 存储支持突变的reads
                            support_reads.setdefault(predict_base, set()).add(read_name)
                            # 计算正负链的支持情况
                            if alignment.is_proper_pair:
                                if alignment.is_read1:
                                    if alignment.is_reverse:
                                        # 突变在read1,并且read1和参考基因组是反向互补的方向
                                        F2R1.setdefault(predict_base, 0)
                                        F2R1[predict_base] += 1
                                    else:
                                        # 突变在read1，并且read1和参考基因组是一致的方向
                                        F1R2.setdefault(predict_base, 0)
                                        F1R2[predict_base] += 1
                                elif alignment.is_read2:
                                    if alignment.is_reverse:
                                        # 突变在read2,并且read1和参考基因组是反向互补的方向
                                        F1R2.setdefault(predict_base, 0)
                                        F1R2[predict_base] += 1
                                    else:
                                        # 突变在read2，并且read1和参考基因组是一致的方向
                                        F2R1.setdefault(predict_base, 0)
                                        F2R1[predict_base] += 1

                # 开始统计分析; 提取主碱基作为参考
                for predict_base in support_reads.keys():
                    # 以主碱基即出现频次最多的碱基进行比较,看是否存在差异,如果存在差异,说明可能假阳性
                    test_reads = support_reads[predict_base]
                    ref_reads = support_reads[main_base]
                    # 后来根据实际情况检测,该比较可能没有大的参考意义,判断不了变异的可靠性
                    bq_pvalue, alt_bqs, ref_bqs = self.qual_diff_test(test_reads, ref_reads, base_qual_dict)
                    mq_pvalue, alt_mqs, ref_mqs = self.qual_diff_test(test_reads, ref_reads, map_qual_dict)
                    pos_pvalue, alt_pos, ref_pos = self.qual_diff_test(test_reads, ref_reads, dist_to_read_end_dict)
                    mean_alt_bq = statistics.mean(alt_bqs) if alt_mqs else None
                    mean_ref_bq = statistics.mean(ref_bqs) if alt_mqs else None
                    mean_alt_mq = statistics.mean(alt_mqs) if alt_mqs else None
                    mean_ref_mq = statistics.mean(ref_mqs) if alt_mqs else None
                    mean_alt_pos = statistics.median(alt_pos) if alt_pos else None
                    alt_pos_pstd = statistics.pstdev(alt_pos) if alt_pos else None
                    max_alt_pos_diff = max(alt_pos) - min(alt_pos)
                    # FoxoG = self.foxog_calc(ref_base, F1R2, F2R1)
                    # 收集指标
                    metric_info = metrics.setdefault(predict_base, dict())
                    # 实践表明,当前的碱基质量比较检验得到的pvalue无决定性的参考意义
                    # metric_info['bq_pvalue'] = bq_pvalue
                    metric_info['mean_bq'] = mean_alt_bq
                    metric_info['mean_ref_bq'] = mean_ref_bq
                    metric_info['mq_pvalue'] = mq_pvalue
                    metric_info['mean_mq'] = mean_alt_mq
                    metric_info['pos_pvalue'] = pos_pvalue
                    metric_info['pos_pstd'] = alt_pos_pstd
                    metric_info['mean_pos'] = mean_alt_pos
                    # 下面这个指标主要用于判断突变位点在read中的位置是否存在差异，如果完全一致，极可能是假的
                    metric_info['max_alt_pos_diff'] = max_alt_pos_diff
                    support_depth = len(support_reads[predict_base])
                    metric_info['support_depth'] = support_depth
                    metric_info['ref_base'] = ref_base
                    metric_info['main_base'] = main_base
                    metric_info['predict_base'] = predict_base
                    metric_info['F1R2'] = F1R2
                    metric_info['F2R1'] = F2R1

                    # 判断snv的质量
                    if mean_alt_bq >= 30 and support_depth >= 2:
                        metric_info['GoodQuality'] = True
                    else:
                        metric_info['GoodQuality'] = False

            # call完成
            call = result.setdefault((contig, pos), dict())
            for predict_base in support_reads:
                call[predict_base] = {
                    'reads': support_reads[predict_base],
                    'metrics': metrics[predict_base],
                }

        # 汇总变异结果
        typing_result = []
        possible_types = [list(call.keys()) for call in result.values()]
        for genotype in itertools.product(*possible_types):
            # 获取read交集的大小 = 支持一个基因型的read数量
            reads = [result[(contig, positions[idx])][base]['reads'] for idx, base in enumerate(genotype)]
            genotype_support_num = len(set.intersection(*reads))
            if genotype_support_num > 0:
                # 提取基因型对应的snp的metrics
                metric_lst = [result[(contig, positions[idx])][base]['metrics'] for idx, base in enumerate(genotype)]
                tmp_dict = dict(marker=marker)
                tmp_dict['haplotype'] = ','.join(genotype)
                tmp_dict['count'] = genotype_support_num
                tmp_dict['metrics'] = metric_lst
                tmp_dict['depths'] = depths
                typing_result.append(tmp_dict)
        return typing_result

    def predict_contributor_number(self):
        """
        可以推测出至少有多少个？
        :return:
        """

    def run(self, out_prefix='typing_result'):
        raw_result = []
        clean_result = []
        for marker_name, marker_info in self.markers.items():
            type_result = self.marker_typing(marker_name, marker_info['contig'], marker_info['snp_pos'])
            raw_result += type_result
            for each in type_result:
                # 要求每个SNP都是GoodQuality
                if all([x["GoodQuality"] for x in each['metrics']]):
                    clean_result.append(each)
        with open(out_prefix+'.raw.json', 'w') as f:
            json.dump(raw_result, f, indent=2)
        with open(out_prefix+'.clean.json', 'w') as f:
            json.dump(clean_result, f, indent=2)
        pd.DataFrame(clean_result).to_csv(out_prefix+'.csv', index=False)


def micro_hap_typing(bam_file, genome_file, micro_hap_file, out_prefix='typing_result'):
    MicroHapCaller(bam_file, genome_file, micro_hap_file, out_prefix=out_prefix)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['micro_hap_typing'])


