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
import os
import re
import json
import pandas as pd
import pysam
import scipy.stats as stats
import statistics
from collections import Counter
from stat_3bases_error import estimate_context_seq_error
"""
要求vcf每一行只包含一个突变，这个可以通过bcftools norm 快速实现
变异类型的分类
https://www.ebi.ac.uk/training-beta/online/courses/human-genetic-variation-introduction/what-is-genetic-variation/types-of-genetic-variation/

考虑MSI的识别, 这样可以过滤MSI突变
"""


class VcfFilter(object):
    def __init__(self, vcf_path, tumor=None, normal=None, normal_vcf=None, gene_primer_used=False):
        self.vcf = pysam.VariantFile(vcf_path)
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
        if L - P > 0.0001:
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
            elif normal_dp >= 300:
                # 将测序深度超过300的，但又无法判定为germline的突变作为背景噪音，用于后续过滤
                af = normal_af
            else:
                # 其他 (depth < 5) or (5 <= depth < 300 and af < 0.25)
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
        print(f'we guess error rate for {record.start}:{record.ref}>{record.alts[0]}:', error_rate)
        return error_rate

    def pass_seq_error(self, record, sample, error_rate:float=None, alpha=0.05, factor=1.0):
        if error_rate == 0:
            raise Exception(record.__str__())
        dp = self.get_depth(record, sample)
        af = self.get_af_value(record, sample)
        # 估计error_rate的置信区间
        confidence = 1 - alpha
        # confidence越大，则 error_upper越大, af_lower越小，min_depth越大， 这意味着过滤条件越严格
        error_lower, error_upper = self.poll_error_binomial_conf(error_rate=error_rate, depth=dp, confidence=confidence)
        af_lower, af_upper = self.poll_error_binomial_conf(error_rate=af, depth=dp, confidence=confidence)
        # 根据二型分布估计突变完全来自背景噪音或测序错误的概率值
        pvalue = self.get_alt_binomial_pvalue(alt_depth=round(dp*af), depth=dp, error_rate=error_rate)
        min_depth = self.estimate_min_required_depth(error_rate, af, confidence)
        min_depth = int(min_depth) * 0.9
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
                if len(alt) > 1:
                    msi_len = int(record.info['MSILEN'])
                    if len(alt) % msi_len == 0:
                        alt = alt[:msi_len]
                if len(alt) * 4 <= 20:
                    ref_seq = left[:-len(alt)*4] + record.ref + right[:len(alt)*4]
                    if re.search(f'({alt})'+'{5,}', ref_seq):
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
            target_info['LOD(eError,UpperConf,Pvalue)'] = r.info['LOD']
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
        with pysam.VariantFile(vcf) as fr:
            for r in fr:
                result.append(self.format_txt_output(r))
        df = pd.DataFrame(result)
        df.to_csv(out, sep='\t', index=False)
        df.to_excel(out[:-4]+'.xlsx', index=False)

    def filtering(self, genome, ref_dict=None, out_prefix=None, min_error_rate=None, error_rate_file=None, min_af=0.001, min_depth=5, alpha=0.05):
        # 先给vcf的header添加新字段定义才能往添加新的字段信息
        self.vcf.header.info.add(
            'LOD', number=5, type='Float',
            description='[error_rate, error_upper, af_lower, pvalue, min_depth]. '
                        'The first value is input error rate which will be used as theoretical frequency to '
                        f'calculate the second value. The second value is the upper bound of the {alpha} confidence interval of error rate. '
                        'The fourth value is the probability of alt observed from background noise'
        )
        self.vcf.header.filters.add('NoiseFromNormal', number=None, type=None, description='noise from normal sample')
        self.vcf.header.filters.add('BackgroundNoise', number=None, type=None, description='noise from background')
        self.vcf.header.filters.add('FromGermline', number=None, type=None, description='considered as germline variant found in normal vcf')
        self.vcf.header.filters.add('HighPopFreq', number=None, type=None, description='Population frequency')
        self.vcf.header.filters.add('LowAF', number=None, type=None, description=f'AF smaller than {min_af}')
        self.vcf.header.filters.add('MSIFilter', number=None, type=None, description=f'likely PCR Slippage in MSI region')
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
        gn = pysam.FastaFile(genome)
        # print(seq_error_dict.keys())
        lod_list = []
        keep = 0
        total = 0
        filter_reasons = []
        log_file = open(out_prefix + '.filtering.log', 'w')
        vcf_loh = pysam.VariantFile(out_prefix+'.LOH.vcf', "w", header=self.vcf.header.copy())
        for r in self.vcf:
            total += 1
            reasons = []
            if '<' in r.alts[0]:
                # 跳过vardict中输出的特殊突变
                print('skip special variant:', r.contig, r.pos, r.ref, list(r.alts), file=log_file)
                continue

            if 'N' in r.alts[0]:
                print('skip "N" containing variant:', r.contig, r.pos, r.ref, list(r.alts), file=log_file)
                continue

            if r.info['DP'] < min_depth:
                print(f'skip variant in shallow depth({r.info["DP"]}):', r.contig, r.pos, r.ref, list(r.alts), file=log_file)
                continue

            # 过滤VCF中原本没有被判定通过的突变
            if list(r.filter)[0] != "PASS":
                reasons = list(r.filter)

            # af cutoff
            af = self.get_af_value(r, self.tumor)
            if af < min_af:
                reasons.append('LowAF')

            # LOH
            if af <= 0:
                # print('discard AF=0 variant as following, maybe it called as LOH site', file=log_file)
                vcf_loh.write(r)
                continue

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
                if normal_af and (normal_af > error_rate):
                    ctrl_af_as_error_rate = True
                    error_rate = normal_af

            if error_rate is None:
                # 没有提取到错误率信息，尝试从record的字段信息提取error_rate的估计值
                print(f'No error rate in seq_eror_file found for the following variant', r.__str__())
                error_rate = self.get_raw_error_rate(r, read_len=150, min_error_rate=1e-6)

            # 1.seq error过滤或者germline突变过滤
            if error_rate > 0.999:
                reasons.append('FromGermline')
            else:
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
            judge2 = self.pass_strand_bias(r, cutoff=0.001, sample=self.tumor)
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

            # 5. 根据MSI状况
            judge5 = self.pass_msi_filter(r)
            if not judge5:
                reasons.append('MSIFilter')

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
        self.write_out_txt(out_vcf_name, out_prefix+'.final.txt')
        print('median LOD:', statistics.median(lod_list), file=log_file)
        print('min LOD:', min(lod_list), file=log_file)
        print('max LOD:', max(lod_list), file=log_file)
        print(f'discard {total - keep} variants while keep {keep} ones!', file=log_file)
        reason_couts = Counter(filter_reasons).most_common()
        for k, v in reason_couts:
            print(f'{v} variants are filtered out because of {k}', file=log_file)


def filter_vcf(vcf, genome, ref_dict=None, tumor_name=None, bam=None, bed=None, normal_vcf=None, alpha=0.01, min_af=0.001,
              exclude_from=None, out_prefix=None, min_error_rate=None, error_rate_file=None, center_size:tuple=(1, 1)):
    if bam and bed and (not error_rate_file):
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
        alpha=alpha,
        min_af=min_af
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
    parser.add_argument('-min_af', type=float, default=0.001, help='hard cutoff of AF')
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
        out_prefix=args.out_prefix
    )
