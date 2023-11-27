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


def get_seq_and_qual(contig, start, end, bam, min_bq=13, get_qual=False):
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
        # header = ['contig', 'pos', 'ref', 'centered_ref', 'depth', 'alt_stat', 'base_qual_stat']
        header = ['contig', 'pos', 'ref', 'centered_ref', 'depth', 'alt_stat']
        fw.write('\t'.join(header)+'\n')
        depth_lst = []
        for line in bed_file:
            if line.startswith('track') or line.startswith('#'):
                continue
            r, s, t = line.strip().split()[:3]
            s, t = int(s), int(t)
            seq_quals = get_seq_and_qual(r, s, t, bam)
            target_info_lst = []
            depths = []
            for pos, seq_qual in zip(range(s, t), seq_quals):
                if (r, str(pos+1)) in excludes:
                    # print('skip', pos)
                    continue
                if len(seq_qual[0]) > 0:
                    # 统计当前位点A/T/C/G/I/D出现的频数
                    seq_counter = Counter(seq_qual[0])
                    # qual_counter = Counter(seq_qual[1])
                    center_seq = get_center_seq(r, pos, gn, sizes=center_size).upper()
                    ref = center_seq[len(center_seq)//2]
                    # 统计alternative碱基的总数
                    alt_types = seq_counter.keys() - {ref, 'N', 'n'}
                    total_depth = sum(seq_counter.values())
                    depths.append(total_depth)
                    # alt_freq = alt_num / (total)
                    for alt_type in alt_types:
                        alt_depth = seq_counter[alt_type]
                        alt_freq = alt_depth/total_depth
                        if alt_freq < 0.01:
                            # 只考虑低频突变
                            alt_raw_dict[alt_type].setdefault(center_seq, Counter())
                            # {'G': 'CAT':{'A': 10, 'G':2}}
                            alt_raw_dict[alt_type][center_seq] += seq_counter

                    if alt_types:
                        target_info_lst.append([
                            r, pos, ref,
                            center_seq,
                            total_depth,
                            dict(seq_counter.most_common()),
                            # dict(qual_counter.most_common())
                        ])
            # write out
            if len(depths) == 0:
                depths = [1]
            median_depth = statistics.median(depths)
            depth_lst.append(median_depth)
            for info in target_info_lst:
                fw.write('\t'.join(str(x) for x in info)+'\n')
    bam.close()
    gn.close()
    median_depth = statistics.median(depth_lst)
    print('Median Depth:', median_depth)

    # 查看字段结构大小
    print(dict(zip(alt_raw_dict.keys(), (len(v) for k, v in alt_raw_dict.items()))))
    alt_dict = dict()
    min_depth = max(median_depth * 5, 5000)
    for alt_type, mdict in alt_raw_dict.items():
        # 例如alt_type='G', mdict={'CAT':{'A': 10, 'G':2}}
        # 合并：由于观察到A->C时，可能是pcr时把A错配成C, 也有可能pcr时把互补链的T错配为G
        result = dict()
        for base in ['A', 'T', 'C', 'G', 'I', 'D']:
            for center_seq in mdict:
                result.setdefault(center_seq, Counter())
                if base not in mdict[center_seq]:
                    continue
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
            if total < min_depth:
                # total是总的depth，如果总的depth不够多，则统计意义将不足
                print(f'For context named "{center_seq}", the total deph < median_depth * 5, skip it for not enough statistic meaning')
                continue
            v = dict(v.most_common())
            # 仅仅输出目标突变的概率
            # freq = {x: v[x]/total for x in v if x==alt_type}
            freq = {alt_type: v[alt_type]/total}
            freq['Depth'] = total
            freq_result[center_seq] = freq
        # 排序
        freq_result = dict(sorted(
            zip(freq_result.keys(), freq_result.values()),
            key=lambda x: sum(x[1].values()), reverse=True)
        )
        if freq_result:
            alt_dict[alt_type] = freq_result

    # print(alt_dict.keys())
    out_file = f'{prefix}.centered{center_size[0]}{center_size[1]}_site.json'
    with open(out_file, 'w') as fw:
        json.dump(alt_dict, fw, indent=4)
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
    args = parser.parse_args()
    estimate_context_seq_error(
        bed=args.bed,
        bam=args.bam,
        prefix=args.out_prefix,
        center_size=args.center_size,
        genome=args.genome,
        exclude_from=args.exclude_from
    )
