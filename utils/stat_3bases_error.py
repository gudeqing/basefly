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

# 重要的假设
当既没有对照样本，又没有阴性样本，我们假设肿瘤样本中call出来的低频突变中绝大部分为假阳性突变，那么我们也可以根据肿瘤样本粗略估计测序错误。
最后可以用肿瘤样本自己作为输入并实现过滤。

【确定检测下限和测错误率，反推所需测序深度】:
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
import pysam
import statistics
from collections import Counter
__author__ = 'gdq'


def get_seq_and_qual(contig, start, end, bam, min_bq=13, get_qual=True):
    cols = bam.pileup(
        contig, start, end,
        stepper='samtools',
        truncate=True,
        min_base_quality=min_bq,
        ignore_orphans=False,
        ignore_overlaps=True,
        max_depth=300000,
    )
    # 针对每一个位点，得到相应覆盖的碱基和碱基质量值
    # 如果相应位置出现删除或插入，则用D和I分别代替表示
    if get_qual:
        seq_qual = [[col.get_query_sequences(add_indels=True), col.get_query_qualities()] for col in cols]
    else:
        seq_qual = [[col.get_query_sequences(add_indels=True)] for col in cols]
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
                seqs.append(i.split('-')[0].upper())
            else:
                seqs.append(i.upper())
        seq_qual[index][0] = seqs
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


def estimate_context_seq_error(bed, bam, prefix, center_size=(1, 1), exclude_from=None, af_cutoff=0.01,
        genome='/nfs2/database/1_human_reference/hg19/ucsc.hg19.fasta'):
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
    :param bed: 要统计的目标区域
    :param bam: 比对结果文件
    :param prefix: 输出结果文件前缀
    :param center_size: 该值为统计时以当前位点为中心并向两边延申的长度，默认为1
    :param exclude_sites: 统计时需要排除的位点列表文件，可以是未经压缩的vcf或bed，通常为已知突变位点信息
    :param af_cutoff: 认为af低于这个频率的突变大部分是假的突变,通过累加这些位点的信息,估计出错误率
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
    # alt_raw_dict = {x: dict() for x in ['A', 'T', 'C', 'G', 'I', 'D']}
    alt_raw_dict = {}
    with open(bed) as bed_file, open(f'{prefix}.each_site.txt', 'w') as fw:
        # header = ['contig', 'pos', 'ref', 'centered_ref', 'depth', 'alt_stat', 'base_qual_stat']
        header = ['contig', 'pos', 'ref', 'centered_ref', 'depth', 'alt_stat']
        fw.write('\t'.join(header)+'\n')
        depth_lst = []
        for line in bed_file:
            if line.startswith('track') or line.startswith('#'):
                continue
            r, s, t = line.strip().split()[:3]
            s, t = int(s), int(t)
            # 返回每个位点的测序信息
            seq_quals = get_seq_and_qual(r, s, t, bam, get_qual=False)
            target_info_lst = []
            depths = []
            for pos, seq_qual in zip(range(s, t), seq_quals):
                seq_bases = seq_qual[0]
                if (r, str(pos+1)) in excludes:
                    # print('skip', pos)
                    continue
                if len(seq_bases) > 0:
                    # 统计当前位点A/T/C/G/I/D出现的频数
                    seq_counter = Counter(seq_bases)
                    # qual_counter = Counter(seq_qual[1])
                    center_seq = get_center_seq(r, pos, gn, sizes=center_size).upper()
                    ref = center_seq[center_size[0]]
                    # 统计
                    total_depth = len(seq_bases)
                    depths.append(total_depth)

                    # 根据碱基频率猜测真实测序模板
                    # 对于center_seq的中间碱基,我们应该要确认样本中真正的参考碱基是什么,我们可以用频率高低来判断
                    predict_base_freq = []
                    for each in seq_counter.keys():
                        base_freq = [each, seq_counter[each], seq_counter[each]/total_depth]
                        predict_base_freq.append(base_freq)
                    # 根据af倒序
                    sorted_base_freq = sorted(predict_base_freq, key=lambda x: x[2], reverse=True)
                    # 我们尽量要以germline信息去做参考,而不是输入的参考基因组作为参考,所以出现下面的center_seq1和center_seq2
                    center_seq1 = None
                    center_seq2 = None
                    if sorted_base_freq[0][2] >= 0.85:
                        # 认为频率超过85%的碱基才是真正的测序模板
                        center_seq = list(center_seq)
                        center_seq[center_size[0]] = sorted_base_freq[0][0]
                        center_seq = ''.join(center_seq)
                    elif len(sorted_base_freq) >= 2 and (sorted_base_freq[0][2] + sorted_base_freq[1][2]) >= 0.85:
                        # 认为频率top2的碱基都是测序模板, 这里我们最多也只考虑真实样本中含量最多的两位模板的情况
                        center_seq = list(center_seq)
                        center_seq[center_size[0]] = sorted_base_freq[0][0]
                        center_seq1 = ''.join(center_seq)
                        # 按照各自的占比去分配测序深度
                        top1_count = sorted_base_freq[0][1]
                        top2_count = sorted_base_freq[1][1]
                        center_seq1_depth = int(top1_count/(top1_count+top2_count) * total_depth) + 1
                        center_seq[center_size[0]] = sorted_base_freq[1][0]
                        center_seq2 = ''.join(center_seq)
                        center_seq2_depth = int(top2_count/(top1_count+top2_count) * total_depth) + 1

                    updated = False
                    for each in sorted_base_freq:
                        # 挑选可能错误的测序结果进行统计
                        if each[2] <= af_cutoff:
                            updated = True
                            # 只考虑低频突变
                            error_base = each[0]
                            error_base_depth = each[1]
                            # 例如: alt_raw_dict = {'CAT':{'A': 1, 'G':2, 'ref_depth':1222}}
                            if center_seq1:
                                alt_raw_dict.setdefault(center_seq1, Counter())
                                alt_raw_dict[center_seq1] += {error_base: error_base_depth}
                            else:
                                alt_raw_dict.setdefault(center_seq, Counter())
                                alt_raw_dict[center_seq] += {error_base: error_base_depth}
                            if center_seq2:
                                alt_raw_dict.setdefault(center_seq2, Counter())
                                alt_raw_dict[center_seq2] += {error_base: error_base_depth}
                    if updated:
                        # ref_depth, 理解为目标模板被测序了多少次, 下面是进行次数累加更新,夯实统计意义
                        if center_seq1:
                            if 'ref_depth' in alt_raw_dict[center_seq1]:
                                alt_raw_dict[center_seq1]['ref_depth'] += center_seq1_depth
                            else:
                                alt_raw_dict[center_seq1]['ref_depth'] = center_seq1_depth
                        else:
                            if 'ref_depth' in alt_raw_dict[center_seq]:
                                alt_raw_dict[center_seq]['ref_depth'] += total_depth
                            else:
                                alt_raw_dict[center_seq]['ref_depth'] = total_depth
                        if center_seq2:
                            if 'ref_depth' in alt_raw_dict[center_seq2]:
                                alt_raw_dict[center_seq2]['ref_depth'] += center_seq2_depth
                            else:
                                alt_raw_dict[center_seq2]['ref_depth'] = center_seq2_depth
                    # 下面按照常规方式去统计每个目标位点的测序信息,不会用于对测序错误率的估计
                    alt_types = seq_counter.keys() - {ref, 'N', 'n'}
                    if alt_types:
                        target_info_lst.append([
                            r, pos, ref, center_seq, total_depth,
                            dict(seq_counter.most_common()),
                            # dict(qual_counter.most_common())
                        ])

            # 区间信息统计完毕,现在输出该区间的测序信息
            for info in target_info_lst:
                fw.write('\t'.join(str(x) for x in info)+'\n')

            # 统计一下每个区间的测序深度信息
            if len(depths) == 0:
                depths = [1]
            median_depth = statistics.median(depths)
            depth_lst.append(median_depth)
    bam.close()
    gn.close()
    # 计算出区间测序深度的中位值的中位值
    median_depth = statistics.median(depth_lst)
    print('Median Depth:', median_depth)
    # print(alt_raw_dict.pop())

    # 转换成错误率信息
    freq_result = {}
    # 要求累加测序深度不能低于5000,否则认为统计意义不足,同时也不能低于中位值的5倍,这些都是经验阈值
    min_depth = max(median_depth * 5, 5000)
    for center_seq, alt_count_dict in alt_raw_dict.items():
        # 举例: center_seq='ACG', alt_count_dict={'A': 10, 'G':2, 'ref_depth':1000}
        freq = freq_result.setdefault(center_seq, dict())
        for error_base in alt_count_dict.keys():
            if error_base != 'ref_depth':
                if alt_count_dict['ref_depth'] >= min_depth:
                    freq[error_base] = dict(
                        error_rate=alt_count_dict[error_base] / alt_count_dict['ref_depth'],
                        depth=alt_count_dict['ref_depth']
                    )
                # 下面根据互补情况进行合并：例如观察到A->C时，可能是pcr时把A错配成C, 也有可能pcr时把互补链的T错配为G
                # r_key = reverse_complement(center_seq)
                # r_base = reverse_complement(error_base)
                # if r_key in alt_raw_dict:
                #     if r_base in alt_raw_dict[r_key]:
                #         # 累加相应alt的数量
                #         new_count = alt_count_dict[error_base] + alt_raw_dict[r_key][r_base]
                #         new_depth = alt_count_dict['ref_depth'] + alt_raw_dict[r_key]['ref_depth']
                #         # 更新具体error_base的信息
                #         freq[error_base] = dict(
                #             af=new_count/new_depth,
                #             depth=new_depth
                #         )

    # 转换为其他格式表述信息
    alt_dict = {}
    for center_seq, freq in freq_result.items():
        for each in ['A', 'T', 'C', 'G', 'I', 'D']:
            for error_base, info in freq.items():
                if error_base == each:
                    tmp_dict = alt_dict.setdefault(error_base, dict())
                    tmp_dict[center_seq] = {error_base: info['error_rate'], "Depth": info['depth']}
    # 按照error_rate排序
    alt_dict_sort = {}
    for alt_base, v_dict in alt_dict.items():
        sorted_v_dict = dict(sorted(
            zip(v_dict.keys(), v_dict.values()),
            key=lambda x: x[1][alt_base],
            reverse=True
        ))
        alt_dict_sort[alt_base] = sorted_v_dict

    # 保存结果
    out_file = f'{prefix}.centered{center_size[0]}{center_size[1]}_site.json'
    with open(out_file, 'w') as fw:
        json.dump(alt_dict_sort, fw, indent=4)
    return out_file


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('-genome', type=Path, required=True, help='path to indexed genome fasta')
    parser.add_argument('-bam', type=Path, required=True, help='bam file which will be used to estimate background noise')
    parser.add_argument('-bed', type=Path,  required=True, help='target region file which will be used to estimate background noise')
    parser.add_argument('-out_prefix', required=True, help='output file prefix')
    parser.add_argument('-center_size', type=int, nargs='+', default=(1, 1), help='extending size around ref base during background noise estimating')
    parser.add_argument('-exclude_from', type=Path, required=False, help='bed or vcf file containing known variant in input bam, these variants will be excluded during background noise estimating')
    parser.add_argument('-af_cutoff', type=float, default=0.01, help='Max af cutoff for error stats')
    args = parser.parse_args()
    estimate_context_seq_error(
        bed=args.bed,
        bam=args.bam,
        prefix=args.out_prefix,
        center_size=args.center_size,
        genome=args.genome,
        exclude_from=args.exclude_from,
        af_cutoff=args.af_cutoff,
    )
