"""
测序错误过滤突变的思路：
基于置信区间，可以参考此链接理解统计原理：https://www.statisticshowto.com/binomial-confidence-interval/
假设事先知道测序过程中发生的碱基转换的概率，例如，知道A被错误测成T的概率为0.001，
针对一个突变位点A-T，其AF=0.008，是否可以将当前位点判定为假阳性呢?
假设还知道当前位点深度为1000，根据这1000 Reads，在95%的置信水平下，估计出来的测序错误概率与真实的错误概率偏差可以使用如下公式推算：
（下面的公式是根据正态分布来的，即用正太分布近似二型分布，理论依据是当N和p都比较大时，二型分布趋近正太分布。）
N = p*(1-P)*(Z/E)**2
=> E = Z/(N/P/(1-P))**0.05
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

"""
import os
import json
import pandas as pd
import pysam
import scipy.stats as stats
import statistics
from collections import Counter
"""
要求vcf每一行只包含一个突变，这个可以通过bcftools norm 快速实现
变异类型的分类
https://www.ebi.ac.uk/training-beta/online/courses/human-genetic-variation-introduction/what-is-genetic-variation/types-of-genetic-variation/
"""


def get_seq_qual(contig, start, end, bam, min_bq=15):
    cols = bam.pileup(
        contig, start, end,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=False,
    )
    # 针对每一个位点，得到相应覆盖的碱基和碱基质量值
    # 如果相应位置出现删除或插入，则用D和I分别代替表示
    seq_qual = [[col.get_query_sequences(add_indels=True), col.get_query_qualities()] for col in cols]
    # A pattern
    # `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion
    # between this reference position and the next reference
    # position. The length of the insertion is given by the
    # integer in the pattern, followed by the inserted
    # sequence. Similarly, a pattern `-[0-9]+[ACGTNacgtn]+'
    # represents a deletion from the reference. The deleted bases
    # will be presented as `*' in the following lines.
    for index in range(len(seq_qual)):
        seqs = []
        for i in seq_qual[index][0]:
            if '+' in i:
                # 当前碱基后面有插入的碱基
                seqs.append('I')
            elif '*' in i:
                # 当前碱基被删除
                seqs.append('D')
            elif '-' in i:
                # 当前碱基的后一个碱基发生删除
                seqs.append(i.split('-')[0])
            else:
                seqs.append(i)
        seq_qual[index][0] = [x.upper() for x in seqs]
    return seq_qual


def get_center_seq(contig, pos, genome, sizes=(1, 1)):
    return genome.fetch(contig, pos-sizes[0], pos+sizes[1]+1)


def reverse_complement(seq):
    """
    :param seq: 输入序列
    :return: 返回reversed的序列
    """
    seq = seq.upper()
    complement = dict(zip(list("ATCG"), list("TAGC")))
    return ''.join(complement[base] if base in complement else base for base in seq[::-1])


def estimate_context_seq_error(bed, bam, prefix, center_size=(1, 1), exclude_from=None,
        genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta',):
    """
    对要目标区域的每一个位点进行统计，尝试得出测序错误条件下：碱基发生转换的条件概率
    分析思路：
    针对6种突变或测序错误情况'ref -> A/T/C/G/I/D' 进行6种统计，其中I:insertion D:deletion
    每一轮统计结果的意义：
    如关注的是ref突变为A，那么针对整个目标区域完成统计后，T/C/G突变为A的概率分别是多少.
    这里还可以假设前后碱基影响测序错误率，所以考虑把ref前后一个碱基纳入，如下统计示例：
    {
    "A": {
        "ACA": { "A": 0.003 }, # 以ACA为context的条件下，C被错误测为A的概率为0.003
        "AGC": { "A": 0.001 },
        "CTT": { "A": 0.0008 }，
        ...},
    ...}
    注:统计时会把测序DP>=10且AF>=10%的位点剔除，认为这些是真实突变，并不是测序错误导致
    :param bed: 要统计的目标区域
    :param bam: 比对结果文件
    :param prefix: 输出结果文件前缀
    :param center_size: 该值为统计时以当前位点为中心并向两边延申的长度，默认为1
    :param exclude_sites: 统计时需要排除的位点列表文件，可以是未经压缩的vcf或bed，通常为已知突变位点信息
    :param genome: 参考基因组序列文件路径, 用于获取参考序列
    :return: 默认输出两个文件。
        一个文件为*.each_site.txt，记录目标区域的每一个碱基的突变情况；
        一个文件为centered{center_size}_site.json'，记录位点中心序列对应的测序错误频率信息。
    """
    gn = pysam.FastaFile(genome)
    bam = pysam.AlignmentFile(bam)
    excludes = set()
    if exclude_from:
        with open(exclude_from) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                lst = line.strip().split()
                if exclude_from.endswith('vcf'):
                    chr_name, start, end = lst[0], int(lst[1]), int(lst[1]) + len(lst[3])-1
                else:
                    chr_name, start, end = lst[0], int(lst[1])+1, int(lst[2])
                for each in range(start, end+1):
                    # 统一转换为1-based
                    excludes.add((chr_name, str(each)))
    # 初始化变异字典
    alt_raw_dict = {x: dict() for x in ['A', 'T', 'C', 'G', 'I', 'D']}
    with open(bed) as bed_file, open(f'{prefix}.each_site.txt', 'w') as fw:
        header = ['contig', 'pos', 'ref', 'centered_ref', 'depth', 'alt_stat', 'base_qual_stat']
        fw.write('\t'.join(header)+'\n')
        for line in bed_file:
            if line.startswith('track') or line.startswith('#'):
                continue
            lst = line.strip().split()
            r, s, t = lst[:3]
            s, t = int(s), int(t)
            seq_quals = get_seq_qual(r, s, t, bam)

            for pos, seq_qual in zip(range(s, t), seq_quals):
                if (r, str(pos+1)) in excludes:
                    # print('skip', pos)
                    continue
                if len(seq_qual[0]) > 0:
                    # 统计当前位点A/T/C/G/I/D出现的频数
                    seq_counter = Counter(seq_qual[0])
                    qual_counter = Counter(seq_qual[1])
                    center_seq = get_center_seq(r, pos, gn, sizes=center_size).upper()
                    ref = center_seq[len(center_seq)//2]
                    # 统计alternative碱基的总数
                    alt_types = seq_counter.keys() - {ref}
                    alt_num = sum(seq_counter[x] for x in alt_types)
                    total = sum(seq_counter.values())
                    alt_freq = alt_num / (total)
                    if alt_freq > 0.05 and total > 20:
                        # 这里的判断可以较大程度的剔除真实突变的位点，目的是为了防止这些突变导致高估测序错误率
                        # 当然，对于测序比较浅的时候，这个过滤的作用不大
                        continue
                    for alt_type in alt_types:
                        alt_raw_dict[alt_type].setdefault(center_seq, Counter())
                        alt_raw_dict[alt_type][center_seq] += seq_counter

                    if alt_types:
                        info = [
                            r, pos, ref,
                            center_seq,
                            total,
                            dict(seq_counter.most_common()),
                            dict(qual_counter.most_common())
                        ]
                        fw.write('\t'.join(str(x) for x in info)+'\n')
    bam.close()
    # 查看碱基之间的转换频率
    print(dict(zip(alt_raw_dict.keys(), (len(v) for k, v in alt_raw_dict.items()))))

    alt_dict = dict()
    for alt_type, mdict in alt_raw_dict.items():
        # 合并：由于观察到A->C时，可能是pcr时把A错配成C, 也有可能pcr时把互补链的T错配为G
        result = dict()
        for base in ['A', 'T', 'C', 'G', 'I', 'D']:
            for center_seq in mdict:
                result.setdefault(center_seq, Counter())
                result[center_seq][base] = mdict[center_seq][base]
                r_key = reverse_complement(center_seq)
                r_base = reverse_complement(base)
                if r_key in mdict and r_base in mdict[r_key]:
                    result[center_seq][base] += mdict[r_key][r_base]
        # convert to freq
        freq_result = dict()
        for center_seq in sorted(result.keys()):
            v = result[center_seq]
            total = sum(v.values())
            # total是总的depth，如果总的depth不够多，则统计意义将不足
            if total < 1000:
                print(f'For context named "{center_seq}", the total depth is only {total}')
                continue
            v = dict(v.most_common())
            # 仅仅输出目标突变的概率
            freq = {x: v[x]/total for x in v if x==alt_type}
            freq_result[center_seq] = freq
        # 排序
        freq_result = dict(sorted(
            zip(freq_result.keys(), freq_result.values()),
            key=lambda x: sum(x[1].values()), reverse=True)
        )
        alt_dict[alt_type] = freq_result

    # print(alt_dict.keys())
    out_file = f'{prefix}.centered{center_size[0]}{center_size[1]}_site.json'
    with open(out_file, 'w') as fw:
        json.dump(alt_dict, fw, indent=4)
    return out_file


class VcfFilter(object):
    def __init__(self, vcf_path, tumor=None, normal_vcf=None, gene_primer_used=False):
        self.vcf = pysam.VariantFile(vcf_path)
        self.vcf_path = vcf_path
        self.tumor = tumor
        self.normal = None
        self.gene_primer_used = gene_primer_used
        samples = list(self.vcf.header.samples)
        self.use_depth_bias = False
        if len(samples) == 1:
            self.tumor = samples[0]
            self.use_depth_bias = False
        else:
            if tumor is None:
                self.tumor = samples[1]
                self.normal = samples[0]
                print(f'assume tumor and normal samples are {self.tumor} and {self.normal}')
            else:
                self.tumor = tumor
                self.normal = list(set(samples) - {tumor})[0]

        if normal_vcf:
            normal_af_dict = dict()
            ctrl_vcf = pysam.VariantFile(normal_vcf)
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
        # 可以根据样本量估算公式反推已知测序错误率和测序深度的条件下，计算测序错误率的上限，作为检测下限
        # e = z/(depth/(error_rate*(1-error_rate)))**0.5
        # lower = 0 if (error_rate - e <= 0) else (error_rate - e)
        lower, upper = stats.binom.interval(confidence=confidence, n=depth, p=error_rate)
        lower, upper = lower/depth, upper/depth
        if lower < 1e-6:
            lower = 1e-6
        return lower, upper

    def get_alt_binomial_pvalue(self, alt_depth, depth, error_rate):
        # 假设测序错误服从二型分布，可以计算alt_depth全部来自错误的概率
        return 1 - stats.binom.cdf(k=alt_depth, n=depth, p=error_rate)

    def qual_to_error_rate(self, base_qual):
        # 碱基质量值还原为概率值
        # Phred = -10 * log(error_p)
        error_prob = 10**(base_qual*(-0.1))
        return error_prob

    def get_af_value(self, record, sample):
        if 'HIAF' in record.info and 'AF' in record.info:
            # vardict style
            af = min([record.info['HIAF'], record.info['AF'][0]])
        elif 'AF' in record.samples[sample]:
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
            normal_dp = record.samples[self.normal]['DP']
            if type(normal_dp) == tuple:
                normal_dp = sum(normal_dp)
            normal_ad = record.samples[self.normal]['AD']
            if type(normal_ad) == tuple:
                normal_ad = normal_ad[1]
            normal_af = normal_ad / normal_dp
            if normal_dp >= 5 and normal_af >= 0.25:
                af = 1.0
            elif normal_af >= 500:
                af = normal_af
            else:
                pass
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
        if 'NM' in record.info and type(record.info['NM']) == float:
            # vardict style, 'NM':"Mean mismatches in reads"
            # 根据平均错配数量估计测序错误率，考虑真实突变和比对错误的存在，这个错误率肯定偏大
            # 保守考虑，将这个错误率缩小10倍
            error_rate = min([record.info['NM'] / read_len * 0.1, error_rate])
        elif 'QUAL' in record.info:
            error_rate = self.qual_to_error_rate(record.info['qual'])
        elif 'MBQ' in record.info:
            # mutect2 style
            error_rate = self.qual_to_error_rate(record.info['MBQ'])
        print(f'error rate for {record.start}:{record.ref}>{record.alts[0]}:', error_rate)
        return error_rate

    def pass_seq_error(self, record, sample, seq_error:float=None, alpha=0.05, read_len=150):
        # 置信水平99%对应Z值为2.58
        # 置信水平95%对应Z值为1.96
        dp = self.get_depth(record, sample)
        af = self.get_af_value(record, sample)
        if seq_error is None:
            # 尝试从record的字段信息提取error_rate的估计值
            error_rate = self.get_raw_error_rate(record, read_len=read_len)
        else:
            error_rate = seq_error
        # 估计error_rate的置信区间
        confidence = 1 - alpha
        lower, upper = self.poll_error_binomial_conf(error_rate=error_rate, depth=dp, confidence=confidence)
        # 根据二型分布估计突变完全来自背景噪音或测序错误的概率值
        pvalue = self.get_alt_binomial_pvalue(alt_depth=round(dp*af), depth=dp, error_rate=error_rate)
        # print(dp, r.qual, error_rate, lower, upper)
        # 理论上，下面两个不等式的判定应该是等效的
        if af >= upper or pvalue < alpha:
            return True, lower, upper, pvalue
        else:
            return False, lower, upper, pvalue

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

    def pass_strand_bias(self, record, info_field='', fmt_field='', cutoff=0.005, sample=None):
        passed = True
        pvalue = None
        possible_info_fields = [info_field] + ['SBF']
        for field in possible_info_fields:
            if field in record.info:
                if type(record.info[field]) == float:
                    pvalue = record.info[field]
                    break
        # find info in "format"
        possible_fmt_fields = [fmt_field] + ['SB', 'DP4']
        if (pvalue is None) and (sample is not None):
            ref_fwd, ref_bwd, alt_fwd, alt_bwd = 0, 0, 0, 0
            for field in possible_fmt_fields:
                if field in record.samples[sample]:
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
            "Consequence": csq_dict['Consequence'],  # Consequence(s)
            "CLIN_SIG": csq_dict['CLIN_SIG'],
            "Exon": csq_dict['EXON'],  # Affected Exon(s)
            "Strand": csq_dict['STRAND'],
            "VariantClass": csq_dict['VARIANT_CLASS'],
            "Depth": self.get_depth(r, self.tumor),  # Depth
            "ExistingVariation": csq_dict['Existing_variation'],
            "MAX_AF": csq_dict['MAX_AF'],
            "MAX_AF_POPS": csq_dict['MAX_AF_POPS'],
            "gnomADe_EAS_AF": csq_dict['gnomADe_EAS_AF'],
            "LSEQ": r.info['LSEQ'],  # from vardict output
            "RSEQ": r.info['RSEQ'],  # from vardict output
        }
        return target_info

    def write_out_txt(self, vcf, out):
        result = []
        with pysam.VariantFile(vcf) as fr:
            for r in fr:
                result.append(self.format_txt_output(r))
        df = pd.DataFrame(result)
        df.to_csv(out, sep='\t', index=False)
        df.to_excel(out[:-4]+'.xlsx', index=False)

    def filtering(self, genome, ref_dict=None, out_prefix=None, min_error_rate=1e-6, error_rate_file=None, alpha=0.05):
        # 先给vcf的header添加新字段定义才能往添加新的字段信息
        self.vcf.header.info.add(
            'LOD', number=3, type='Float',
            description='The first value is input error rate which will be used as theoretical frequency to '
                        f'calculate the second value. The second value is the upper bound of the {alpha} confidence interval of error rate. '
                        'The third value is the probability of alt observed from background noise'
        )
        self.vcf.header.filters.add('NoiseFromNormal', number=None, type=None, description='noise from normal sample')
        self.vcf.header.filters.add('BackgroundNoise', number=None, type=None, description='noise from background')
        self.vcf.header.filters.add('HighPopFreq', number=None, type=None, description='Population frequency')
        if 'Bias' not in self.vcf.header.filters:
            self.vcf.header.filters.add('Bias', number=None, type=None, description="severe strand bias")
        if self.vcf.header.contigs.__len__() < 1:
            self.add_contig_header(ref_dict)

        if out_prefix is None:
            out_prefix = self.vcf_path.rsplit('.', 1)[0]
        out_vcf_name = out_prefix+'.final.vcf'
        vcf_out = pysam.VariantFile(out_vcf_name, "w", header=self.vcf.header.copy())
        vcf_discard = pysam.VariantFile(out_prefix+'.discarded.vcf', "w", header=self.vcf.header.copy())

        # 如果仅仅提供了一个全局的最低错误率信息
        key_left = key_right = 0
        seq_error_dict = None
        if min_error_rate:
            seq_error_dict = dict()
            for b in 'ATCGID':
                for i in set('ATCG') - {b}:
                    seq_error_dict.setdefault(b, dict())[i] = {b: float(min_error_rate)}
        if error_rate_file:
            seq_error_dict = json.load(open(error_rate_file))
            key_len = max(len(x) for x in seq_error_dict['A'].keys())//2
            key_left = key_right = key_len

        if error_rate_file is None:
            print(f'{seq_error_dict} is used as the background noise information')

        # 读取参考基因组信息
        gn = pysam.FastaFile(genome)
        # print(seq_error_dict.keys())
        lod_list = []
        discard = 0
        total = 0
        filter_reasons = []
        for r in self.vcf:
            total += 1
            reasons = []
            if '<' in r.alts[0]:
                discard += 1
                # 跳过vardict中输出的特殊突变
                print('skip', r.contig, r.ref, list(r.alts))
                continue

            if 'N' in r.alts[0]:
                print('skip "N" containing variant', r.contig, r.ref, list(r.alts))
                continue

            # 过滤VCF中原本没有被判定通过的突变
            if list(r.filter)[0] != "PASS":
                reasons = list(r.filter)

            # 1.根据测序错误率或germline突变频率过滤
            ctrl_af_as_error_rate = False
            # r.pos正好是1-based, r.start 是0-based
            error_rate = min_error_rate
            mutation_type = self.get_mutation_type(r)
            if mutation_type == 'SNV':
                key = gn.fetch(r.contig, r.start - key_left, r.start + 1 + key_right).upper()
                # error_rate = seq_error_dict[key][r.alts[0]]
                if key in seq_error_dict[r.alts[0]]:
                    error_rate = seq_error_dict[r.alts[0]][key][r.alts[0]]
            elif mutation_type == 'Insertion':
                key = gn.fetch(r.contig, r.start - key_left, r.start + 1 + key_right).upper()
                if key in seq_error_dict['I']:
                    error_rate = seq_error_dict['I'][key]['I']
                    if len(r.alts[0])-len(r.ref) >= 3:
                        error_rate = error_rate**2
            elif mutation_type == 'Deletion':
                # deletion
                key = gn.fetch(r.contig, r.pos - key_left, r.pos + 1 + key_right).upper()
                if 'D' in seq_error_dict['D']:
                    # error_rate = seq_error_dict[key]['']
                    if key in seq_error_dict['D']:
                        error_rate = seq_error_dict['D'][key]['D']
                        if len(r.ref) - len(r.alts[0]) >= 3:
                            error_rate = error_rate ** 2
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
                # 这样既可以过滤掉germline突变
                # 如果对照样本测序深度足够，还可以过滤假阳性突变、克隆造血突变
                # 如果对照样本测序深度不够，germline的突变AF可能被低估，而肿瘤样本由于测序深度足够，可以检测到远高于对照样本的AF
                # 此时如果对照样本的AF作为error_rate, 肿瘤样本有可能通过，所以直接扔掉被判定为germline的突变
                normal_af = self.get_normal_af(r)
                if normal_af and (normal_af > error_rate):
                    ctrl_af_as_error_rate = True
                    error_rate = normal_af

            # 1.seq error过滤或者germline突变过滤
            judge = self.pass_seq_error(r, self.tumor, error_rate, alpha=alpha)
            if not judge[0]:
                if ctrl_af_as_error_rate:
                    reasons.append('NoiseFromNormal')
                else:
                    # print(key, error_rate, r.ref, list(r.alts))
                    reasons.append('BackgroundNoise')
                pass
            r.info['LOD'] = (round(error_rate, 5), round(judge[2], 5), round(judge[3], 7))
            lod_list.append(judge[2])

            # 2.根据strand bias进行过滤
            judge2 = self.pass_strand_bias(r, cutoff=0.003, sample=self.tumor)
            if not judge2[0]:
                if 'Bias' not in reasons:
                    reasons.append('Bias')
                    pass

            # 3. position std filtering, 如果突变在reads中出现的位置基本不变,且支持的read<=2，需要过滤掉
            judge3 = self.pass_pstd(r, cutoff=0.00001, gene_primer_used=self.gene_primer_used)
            if not judge3[0]:
                if 'pSTD' not in reasons:
                    reasons.append('pSTD')

            # 4. 根据人群频率过滤
            judge4 = self.pass_population_af(r, cutoff=0.01)
            if not judge4[0]:
                reasons.append('HighPopFreq')

            # 更新filter的内容和输出结果
            if reasons:
                filter_reasons.append(tuple(reasons))
                discard += 1
                for each in reasons:
                    r.filter.add(each)
                vcf_discard.write(r)
            else:
                vcf_out.write(r)

        vcf_out.close()
        vcf_discard.close()
        gn.close()
        self.write_out_txt(out_vcf_name, out_prefix+'.final.txt')
        print('median LOD:', statistics.median(lod_list))
        print('min LOD:', min(lod_list))
        print('max LOD:', max(lod_list))
        print(f'discard {discard} variants while keep {total-discard} ones!')
        reason_couts = Counter(filter_reasons).most_common()
        for k, v in reason_couts:
            print(f'{v} variants are filtered out because of {k}')


# 如果提供对照样本的vcf，利用对照样本的vcf信息进行过滤
# 输出时，按照室间质评的方式输出
# Chr Start End Ref Alt Gene Type Transcript cHGVS pHGVS VAF (%)

def filterVcf(vcf, genome, ref_dict=None, tumor_name=None, bam=None, bed=None, normal_vcf=None, alpha=0.05,
              exclude_from=None, out_prefix=None, min_error_rate=1e-6, error_rate_file=None, center_size:tuple=(1, 1)):
    if bam and bed:
        error_rate_file = estimate_context_seq_error(
            bed, bam, prefix=out_prefix, center_size=center_size,
            genome=genome, exclude_from=exclude_from
        )
    vcf = VcfFilter(vcf_path=vcf, tumor=tumor_name, normal_vcf=normal_vcf, gene_primer_used=False)
    vcf.filtering(
        genome=genome,
        ref_dict=ref_dict,
        out_prefix=out_prefix,
        min_error_rate=min_error_rate,
        error_rate_file=error_rate_file,
        alpha=alpha
    )


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())

