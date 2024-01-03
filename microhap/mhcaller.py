import os
import json
import statistics
import itertools
import json
import numpy as np
import pandas as pd
from collections import Counter
from scipy import stats
from pysam import FastaFile, VariantFile, AlignmentFile, VariantHeader, VariantRecord

__author__ = 'gdq'


class MicroHapCaller(object):
    def __init__(self, bam_file, genome_file, micro_hap_file, out_prefix, error_rate_file=None):
        self.bam = AlignmentFile(bam_file)
        self.genome = FastaFile(genome_file)
        self.markers = self.parse_micro_hap(micro_hap_file)
        if error_rate_file:
            if error_rate_file.endswith('json'):
                self.seq_error_dict = json.load(open(error_rate_file))
            else:
                error_dict = dict()
                base_mapping = dict(A='T', C='G', T='A', G="C")
                with open(error_rate_file) as f:
                    for line in f:
                        if line.startswith('#') or line.startswith('REF_BASE'):
                            continue
                        if line.strip():
                            lst = line.strip().split()
                            alt, rate = lst[2], lst[-1]
                            error_dict[alt] = float(rate)
                            r, a = alt.split('>')
                            error_dict[base_mapping[r]+'>'+base_mapping[a]] = float(rate)
                self.seq_error_dict = error_dict
                print(self.seq_error_dict)
        else:
            self.seq_error_dict = dict()

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
        if sum(x != y for x, y in zip(alt_quals, ref_quals)) == 0:
            # 两组数据一模一样
            return 1.0, alt_quals, ref_quals
        # 比较
        # 优先进行mannwhitneyu检验，适合样本量较多的情况
        if len(alt_quals) >= 10 and len(ref_quals) >= 10:
            # 如果alt_quals和ref_quals完全一样,有可能导致有些版本的stats.mannwhitneyu报错
            u, pvalue = stats.mannwhitneyu(alt_quals, ref_quals, alternative=alternative)
            # u, pvalue = stats.ranksums(alt_quals, ref_quals, alternative=alternative)
            # u, pvalue = stats.ttest_ind(alt_quals, ref_quals, alternative=alternative)
        elif len(alt_quals) >= 3 and len(ref_quals) >= 3:
            # 非配对独立t检验
            try:
                s, pvalue = stats.ttest_ind(alt_quals, ref_quals, alternative=alternative)
            except TypeError:
                s, pvalue = stats.ttest_ind(alt_quals, ref_quals)
        else:
            # 统计意义不足
            pvalue = 1.0
        # 对为nan的pvalue值进行转换
        if pvalue != pvalue:
            pvalue = 1.0
        return pvalue, alt_quals, ref_quals

    @staticmethod
    def pass_seq_error(dp, af, error_rate:float=None, alpha=0.05):
        if error_rate == 0:
            return True, error_rate, 1
        # 估计error_rate的置信区间
        confidence = 1 - alpha
        # confidence越大，则 error_upper越大, af_lower越小，min_depth越大， 这意味着过滤条件越严格
        try:
            lower, upper = stats.binom.interval(confidence, n=dp, p=error_rate)
        except:
            raise Exception(f'{dp}|{af}|{error_rate}')
        error_lower, error_upper = lower / dp, upper / dp
        # 假设测序错误服从二型分布，可以计算alt_depth全部来自错误的概率
        pvalue = 1 - stats.binom.cdf(k=round(dp*af), n=dp, p=error_rate)
        # error rate不能为0，af不能为1，否则会报错
        # print(dp, r.qual, error_rate, lower, upper)
        # 测试发现pvalue<alpha时，af 不一定小于upper，说明这可能是两种过滤策略
        if (af >= error_upper) and (pvalue <= alpha):
            return True, error_upper, pvalue
        else:
            return False, error_upper, pvalue

    def marker_typing(self, marker, contig, positions, min_bq=10, min_mean_bq=25, min_support_depth=2):
        """
        # 下面是一个值得思考的现象:
        marker	haplotype	count
        mh17CP-001	G,G,G	564 (true)
        mh17CP-001	G,G,A	15 (可以被error rate的过滤策略排除)
        mh17CP-001	G,A,G	26
        mh17CP-001	C,G,G	17
        mh17CP-001	C,A,G	666 (true)
        mh17CP-001	C,A,T	1 (可以被最低阈值2排除)
        mh17CP-001	C,T,G	1 (可以被最低阈值2排除)
        mh17CP-001	T,A,G	1 (可以被最低阈值2排除)
        mh17CP-001	A,A,G	1 (可以被最低阈值2排除)

        CAG可以是GGG中第2个碱基测序错误的结果, 也可能是等位基因型CAG中第1个碱基测序错误的结果
        CGG可以是GGG中第1个碱基测序错误的结果, 也可能是等位基因型CAG中第2个碱基测序错误的结果
        所以我们才会看到上述2个噪音值比较大

        再比如下面,也可以看到类似的高噪音值
        marker	haplotype	count
        mh17CP-002	G,T,T	324   [332, 324, 678]
        mh17CP-002	G,G,T	8     [332, 345, 678](高噪音值, 可以是GTT的中间碱基测序错误导致,也可以是等位基因型AGT的第一个碱基值测序错误导致)
        mh17CP-002	A,G,T	333   [335, 345, 678]
        mh17CP-002	A,G,C	1
        mh17CP-002	A,G,G	1
        mh17CP-002	C,G,T	1
        假设GTT, AGT, GGT都是真的, GTT和AGT是等位, GGT是mix进来的,大概率是纯合子,那么比例是324:4
        # 根据这个现象,我们需要发展更好的过滤策略: 见count_cutoff
        :param marker: marker名称
        :param contig: 染色体名称
        :param positions: marker的SNP坐标
        :param min_bq: 最小碱基质量值,pileup的输入
        :param min_mean_bq: 最小平均碱基质量
        :param min_support_depth: 支持alt的最小read数量
        :return:
        """
        result = dict()
        # result的key是position=(contig, start, ref_base)，value是相应position的碱基测序结果
        depths = []
        for pos in positions:
            # 期望这里的pos是0-based
            ref_base = self.genome.fetch(contig, pos-1, pos+2).upper()
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
                    if query_pos is not None and pileup_read.indel == 0:
                        # 对于pair，只会记录其中一条记录
                        # 提取目标位置检测到的具体碱基成分，其可能是和ref一致，也可能不一致
                        predict_base = query_seq[query_pos].upper()
                        support_reads.setdefault(predict_base, set())
                        F2R1.setdefault(predict_base, 0)
                        F1R2.setdefault(predict_base, 0)
                        if read_name not in support_reads[predict_base]:
                            # 存储支持突变的reads
                            support_reads[predict_base].add(read_name)
                            # 计算突变位点距离read两端的距离的最小值, 针对multiplexPCR得到的，这个值没有什么用处
                            dist_to_read_end_dict[read_name] = min(query_pos, query_len - query_pos)
                            # 计算正负链的支持情况, single end比对时将无记录, 对于
                            if alignment.is_read1:
                                if alignment.is_reverse:
                                    # 突变在read1,并且read1和参考基因组是反向互补的方向
                                    F2R1[predict_base] += 1
                                else:
                                    # 突变在read1，并且read1和参考基因组是一致的方向
                                    F1R2[predict_base] += 1
                            elif alignment.is_read2:
                                if alignment.is_reverse:
                                    # 突变在read2,并且read1和参考基因组是反向互补的方向
                                    F1R2[predict_base] += 1
                                else:
                                    # 突变在read2，并且read1和参考基因组是一致的方向
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
                    # mean_ref_bq = statistics.mean(ref_bqs) if alt_mqs else None
                    mean_alt_mq = statistics.mean(alt_mqs) if alt_mqs else None
                    # mean_ref_mq = statistics.mean(ref_mqs) if alt_mqs else None
                    mean_alt_pos = statistics.median(alt_pos) if alt_pos else None
                    alt_pos_pstd = statistics.pstdev(alt_pos) if alt_pos else None
                    max_alt_pos_diff = max(alt_pos) - min(alt_pos)
                    # FoxoG = self.foxog_calc(ref_base, F1R2, F2R1)
                    # 收集指标
                    metric_info = metrics.setdefault(predict_base, dict())
                    # 实践表明,当前的碱基质量比较检验得到的pvalue无决定性的参考意义
                    # metric_info['bq_pvalue'] = bq_pvalue
                    metric_info['mean_bq'] = mean_alt_bq
                    # metric_info['mean_ref_bq'] = mean_ref_bq
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
                    metric_info['F1R2'] = F1R2[predict_base]
                    metric_info['F2R1'] = F2R1[predict_base]
                    # metric_info['error_rate'] = 0
                    # if self.seq_error_dict:
                    #     if predict_base in self.seq_error_dict and ref_base in self.seq_error_dict[predict_base]:
                    #         metric_info['error_rate'] = self.seq_error_dict[predict_base][ref_base][predict_base]
                    metric_info['AF'] = support_depth/depth
                    # 判断snv的质量
                    metric_info['GoodQuality'] = False
                    if mean_alt_bq >= min_mean_bq and support_depth >= min_support_depth:
                        metric_info['GoodQuality'] = True
                        # judge, error_upper, pvalue = self.pass_seq_error(dp=depth, af=metric_info['AF'], error_rate=metric_info['error_rate'], alpha=0.05)
                        # metric_info['error_upper'] = error_upper
                        # metric_info['error_pvalue'] = pvalue
                        # if judge:
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
                tmp_dict['haplotype'] = '-'.join(genotype)
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

    def haplotype_convert_rate(self, current_hap, hap):
        hap1 = current_hap.split('-')
        hap2 = hap.split('-')
        error_rates = []
        for b1, b2 in zip(hap1, hap2):
            if b1 != b2:
                alt_mode = b2 + '>' + b1
                if alt_mode in self.seq_error_dict:
                    error_rates.append(self.seq_error_dict[alt_mode])
        if error_rates:
            multi_value = 1
            for each in error_rates:
                multi_value *= each
            return multi_value

    def run(self, out_prefix='typing_result'):
        result = []
        for marker_name, marker_info in self.markers.items():
            type_result = self.marker_typing(marker_name, marker_info['contig'], marker_info['snp_pos'])
            for each in type_result:
                # 要求每个SNP都是GoodQuality
                each["support_depths"] = [x["support_depth"] for x in each['metrics']]
                if all([x["GoodQuality"] for x in each['metrics']]):
                    each['Pass'] = True
                else:
                    each['Pass'] = False
                result.append(each)
        # 定制错误率相关cutoff,我们比较同一个marker的不同基因型的差异,根据碱基错误率计算一个阈值
        # 收集marker和haplotype的关系
        marker2type = dict()
        for marker_info in result:
            tmp_dict = marker2type.setdefault(marker_info['marker'], dict())
            tmp_dict[marker_info['haplotype']] = marker_info['count']
        # 开始统计计算
        for marker_info in result:
            count_ratio = marker_info['count'] / sum(marker2type[marker_info['marker']].values())
            marker_info['count_freq'] = round(count_ratio, 4)
            current_hap = marker_info['haplotype']
            cutoff = 0
            for hap, count in marker2type[marker_info['marker']].items():
                convert_rate = self.haplotype_convert_rate(current_hap, hap)
                if convert_rate:
                    lower, upper = stats.binom.interval(0.95, n=count+marker_info['count'], p=convert_rate)
                    cutoff += upper
            marker_info['count_cutoff'] = cutoff
            if marker_info['Pass']:
                # 进一步降低假阳性
                if marker_info['count'] <= marker_info['count_cutoff']:
                    marker_info['Pass'] = False

        # 汇总所有Pass信息,统计一个基因型被测错的指标, 然后用这个指标重新过滤
        total_good = 0
        total_bad = 0
        for type_result in result:
            if type_result['Pass']:
                total_good += type_result['count']
            else:
                total_bad += type_result['count']
        bad_ratio = total_bad / (total_good + total_bad)
        print('Bad typing ratio', bad_ratio)

        # save out file
        with open(out_prefix+'.json', 'w') as f:
            json.dump(result, f, indent=2)

        df = pd.DataFrame(result)
        df = df[["marker", "haplotype", "count", "count_cutoff", "count_freq", "Pass", "support_depths", "depths",  "metrics"]]
        df.to_csv(out_prefix+'.csv', index=False)


def micro_hap_caller(bam_file, genome_file, micro_hap_file, out_prefix='typing_result', error_rate_file=None):
    MicroHapCaller(bam_file, genome_file, micro_hap_file, out_prefix=out_prefix, error_rate_file=error_rate_file)


def call_diploid_profile(profile, freq_cutoff=0.15, out_json=None):
    """
    主要目的是确定一个人的marker的具体基因型信息, 不适用于存在DNA混合的样本分析
    :param profile: 来源于本工具的检测结果
    # 本工具的检测结果示例如下:
    marker	haplotype	count	count_cutoff	count_freq	Pass
    mh01CP-016	T-G-A	436	0	0.6469	TRUE
    mh01CP-016	T-A-A	230	0	0.3412	TRUE
    mh01CP-010	C-C-G	606	0	0.8978	TRUE
    mh01CP-010	T-C-A	63	0	0.0933	TRUE
    mh02CP-004	C-T-A	1584	0	0.9888	TRUE
    mh04CP-003	A-C-T	222	0	0.1336	TRUE
    mh04CP-003	A-C-C	841	0	0.506	TRUE
    mh04CP-003	G-C-C	576	0	0.3466	TRUE
    mh04CP-004	T-T-C	441	0	0.4168	TRUE
    mh04CP-004	C-C-T	591	0	0.5586	TRUE
    :param freq_cutoff: 如果一个基因型的频率低于该值则不参与统计
    :param out_json: 指定输出到json文件
    :return:
    """
    markers = dict()
    with open(profile) as f:
        header = f.readline()
        # {marker: allele: {sources: []}
        for line in f:
            marker, haplotype, count, count_cutoff, count_freq, flag = line.strip().split(',')[:6]
            count_freq = float(count_freq)
            # 这里我们假设是一个人结果,且marker都是双等位基因,因此预期频率应该都较高
            if count_freq >= freq_cutoff and flag == 'True':
                marker_info = markers.setdefault(marker, dict())
                marker_info.setdefault('Alleles', []).append(haplotype)
                marker_info.setdefault('count_freq', []).append(count_freq)
        # 检查marker是否是双等位的, 如果不是将提取top2
        for marker in markers:
            alleles = markers[marker]['Alleles']
            af = markers[marker]['count_freq']
            if len(alleles) >= 3:
                print(f'Marker {marker} has more than 2 genotypes, and we will use only top 2')
                allele_freq = zip(alleles, af)
                allele_freq = sorted(allele_freq, key=lambda x: x[1], reverse=True)
                markers[marker]['Alleles'] = allele_freq[0][:2]
                markers[marker]['count_freq'] = allele_freq[1][:2]
                print(f'The top 2 alleles are: {" | ".join(allele_freq[0][:2])}')
                if sum(allele_freq[0][:2]) < 0.8:
                    print(f'Warn: but the frequency of {marker}:{alleles} is too low (={af}))')
            else:
                if sum(af) < 0.8:
                    print(f'Warn: the frequency of {marker}:{alleles} is too low (={af}))')
    if out_json:
        with open(out_json, 'w') as f:
            json.dump(markers, f, indent=2)
    return markers


def get_marker_chemerism(donor_count: dict, recipient_count: dict):
    """
    在已知贡献者的marker基因型信息前提下,推算混合比例
    # 同为纯合子, 且不同
    Donor:      AA
    Recipient:  BB
    Percent recipient chimerism: B/(A+B)

    # 供者为纯合子, 受者为杂合子, 且和供者共享一个allele
    Donor:      AA
    Recipient:  AB
    Percent recipient chimerism: 2B/(A+B)

    # 供者为纯合子,受者为杂合子,但两者没有共享allele
    Donor:      AA
    Recipient:  CD
    Percent recipient chimerism: (C+D)/(A+C+D)

    # 供者为杂合子, 受者为纯合子, 且和供者分享一个allele
    Donor:      AB
    Recipient:  AA
    Percent recipient chimerism: (A-B)/(A+B)

    # 供者为杂合子, 受者为纯合子, 且和供者完全不同
    Donor:      AB
    Recipient:  CC
    Percent recipient chimerism: C/(A+B+C)

    # 供者为杂合子,受者也为杂合子,且共享一个allele
    Donor:      AB
    Recipient:  BC
    Percent recipient chimerism: 2C/(A+B+C)

    # 供者为杂合子, 受者也为杂合子,且和供者完全不同
    Donor:      AB
    Recipient:  CD
    Percent recipient chimerism: (C+D)/(A+B+C+D)

    最后归纳总结发现,本质只有3种计算方式方式
    :param donor_count: {'ATC': 5, 'AGC':3}
    :param recipient_count: {'ATC': 5, 'ATT':4}
    :return: Percent recipient chimerism or None
    """
    recipient_uniq = list(recipient_count.keys() - donor_count.keys())
    donor_uniq = list(donor_count.keys() - recipient_count.keys())
    if len(recipient_uniq) == len(donor_uniq) == 0:
        # print('两个来源的marker基因型一样, 因此信息无效')
        return None
    if len(recipient_count) == 1:
        recipient_is_homo = True
    else:
        recipient_is_homo = False
    total_count = sum(dict(donor_count, **recipient_count).values()) or 1
    # 针对不同情形计算
    comm = donor_count.keys() & recipient_count.keys()
    if len(comm) == 0:
        # 两者等位基因型完全不一样
        target_count = sum(recipient_count.values())
        return target_count / total_count, (target_count, total_count)
    else:
        # 两者等位基因型存在交集
        if not recipient_is_homo:
            target_count = 2 * recipient_count[recipient_uniq[0]]
            return target_count / total_count, (target_count, total_count)
        else:
            recipient = list(recipient_count.keys())[0]
            target_count = (recipient_count[recipient] - donor_count[donor_uniq[0]])
            return target_count / total_count, (target_count, total_count)


def chemerism(donor_profile, recipient_profile, test_profile, out_json='chimerism.json', adjust_count=False):
    # 提取有效marker, 相同基因的marker无法提供区分信息
    # 仅仅处理双等位基因的情况
    person2markers = dict()
    markers_lst = []
    sources = []
    for ref_profile in [donor_profile, recipient_profile]:
        source = os.path.basename(ref_profile).split('.')[0]
        sources.append(source)
        markers = call_diploid_profile(ref_profile)
        markers_lst.append(set(markers.keys()))
        person2markers[source] = markers

    # 检查是否使用了相同的marker
    comm_markers = set.intersection(*markers_lst)
    if not comm_markers:
        raise Exception('No common markers found between contributors')
    else:
        print(f'There are {len(comm_markers)} common markers found between contributors')

    # 检测有效marker
    informative_markers = []
    for marker in comm_markers:
        alleles1 = set(person2markers[sources[0]][marker]['Alleles'])
        alleles2 = set(person2markers[sources[1]][marker]['Alleles'])
        if (alleles1 - alleles2) != (alleles2 - alleles1):
            informative_markers.append([alleles1, alleles2])
    print(f'There are {len(informative_markers)} informative markers', informative_markers)

    # 分析比例
    detected = dict()
    detected_example = {
        "person1": {
            "marker1": {"A1": 3, "A2": 3},
            "marker2": {"A1": 3, "A2": 7},
        },
        'person2': {
            "marker1": {"A1": 3, "A2": 3},
            "marker2": {"A1": 2, "A2": 7},
        }
    }
    # 初始化detected
    for marker in comm_markers:
        for person in person2markers.keys():
            person_info = detected.setdefault(person, dict())
            person_marker_info = person_info.setdefault(marker, dict())
            for allele in person2markers[person][marker]['Alleles']:
                person_marker_info[allele] = 0
    # 更新detected
    with open(test_profile) as f:
        _header = f.readline()
        marker_found = set()
        # {marker: allele: {sources: []}
        not_assigned_hap = dict()
        for line in f:
            marker, haplotype, count, count_cutoff, count_freq, flag = line.strip().split(',')[:6]
            if flag != 'True':
                not_assigned_hap.setdefault(marker, dict())
                not_assigned_hap[marker].update({haplotype: [count, count_cutoff, count_freq, flag]})
                continue
            if marker not in comm_markers:
                print(f'skip marker {marker} that is not in ref contributors')
            else:
                marker_found.add(marker)
                hap_assigned = False
                for person in detected.keys():
                    if haplotype in detected[person][marker]:
                        # chimerism的精确性很大程度依赖这一步的分配是否正确了
                        # count比较小的,和其他主haplotype仅仅相差一个碱基的,往往很难区分是测序错误还是真实存在
                        hap_assigned = True
                        # detected[person][marker][haplotype] = int(count)
                        detected[person][marker][haplotype] = max(int(count) - float(count_cutoff), 0)
                if not hap_assigned:
                    not_assigned_hap.setdefault(marker, dict())
                    not_assigned_hap[marker].update({haplotype: [count, count_cutoff, count_freq, flag]})

    # 可以考虑尝试将没有分配的haplotype进行重新分配
    # 当然但这里的潜在假设还是: 没有分配的haplotype是人工产物
    # 分配策略: 只重新分配与目标haplotype相差一个碱基的haplotype
    # with open('not_assigned_haplotype.json', 'w') as f:
    #     json.dump(not_assigned_hap, f, indent=2)
    if adjust_count:
        for marker in not_assigned_hap:
            if not(marker in detected[sources[0]] and marker in detected[sources[1]]):
                continue
            donor_count = detected[sources[0]][marker]
            recipient_count = detected[sources[1]][marker]
            assign_to = dict()
            merged_count_dict = dict(donor_count, **recipient_count)
            for un_hap in not_assigned_hap[marker]:
                for target_hap in merged_count_dict.keys():
                    # 找到只有一个碱基差异的目标haplotype
                    if sum(x != y for x, y in zip(un_hap, target_hap)) == 1:
                        assign_to.setdefault(un_hap, []).append(target_hap)
            for un_hap, target_haps in assign_to.items():
                un_hap_count = int(not_assigned_hap[marker][un_hap][0])
                target_current_count = [merged_count_dict[x] for x in target_haps]
                for target in target_haps:
                    print(marker)
                    if target in donor_count:
                        ratio = merged_count_dict[target]/sum(target_current_count)
                        print(f'{un_hap}分配给donor{target}:', un_hap_count * ratio)
                        donor_count[target] += un_hap_count * ratio
                    if target in recipient_count:
                        ratio = merged_count_dict[target] / sum(target_current_count)
                        print(f'{un_hap}分配给recipient{target}:', un_hap_count * ratio)
                        recipient_count[target] += un_hap_count * ratio

    # 更新完成count信息后计算嵌合率
    # 统计是否存在预期marker丢失
    loss_markers = comm_markers - marker_found
    if loss_markers:
        print('Some marker are not found in testing sample:', loss_markers)
    percent_dict = dict()
    recipient_count_lst = []
    total_count_lst = []
    for marker in marker_found:
        # 这个过程没有考虑到有些marker的部分基因型的count为0, 从而误判基因型
        chimer_info = get_marker_chemerism(detected[sources[0]][marker], detected[sources[1]][marker])
        if chimer_info is not None:
            percent_dict[marker] = chimer_info[0]
            recipient_count_lst.append(chimer_info[1][0])
            total_count_lst.append(chimer_info[1][1])
    if percent_dict:
        mean_percent = statistics.mean(percent_dict.values())
        print(f"We found {len(percent_dict)} Informative Markers", list(percent_dict.keys()))
        print(f'Mean Recipient({sources[1]}) Chimerism: ', mean_percent)
        # 线性拟合的方法求得比例
        coeff, _ = np.polyfit(total_count_lst, recipient_count_lst, 1)
        if out_json:
            percent_dict = dict(sorted(zip(percent_dict.keys(), percent_dict.values()), key=lambda x: x[1]))
            percent_dict['mean_chimerism'] = mean_percent
            percent_dict['linear_chimerism'] = coeff
            # 加权平均,以depth作为加权
            weight = [x/sum(total_count_lst) for x in total_count_lst]
            percent_dict['depth_weighted_chimerism'] = sum(x * y for x, y in zip(weight, percent_dict.values()))
            with open(out_json, 'w') as f:
                json.dump(percent_dict, f, indent=2)
    else:
        print('No informative marker found!')
    return percent_dict


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

"""
测试:
cd /home/hxbio04/data/PRJNA508621/Result/analysis/SRR8284658-12-1v3/fullrefr
python /home/hxbio04/basefly/microhap/mhcaller.py micro_hap_caller -bam_file *-fullrefr.bam -genome_file /home/hxbio04/data/PRJNA807084/mh20panel_bymicrohabdb/MicroHapulator/microhapulator/data/hg38.fasta -micro_hap_file ../../../marker-definitions.tsv -out_prefix SRR8284658-12-1v3.call
python /home/hxbio04/basefly/microhap/mhcaller.py chemerism -donor_profile SRR8284710-1.call.csv -recipient_profile SRR8284713-2.call.csv -test_profile SRR8284658-12-1v3.call.csv -out 1v3.chimersim.json
"""
