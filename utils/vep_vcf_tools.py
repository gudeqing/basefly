import os
from collections import Counter
import statistics as sts
from pysam import VariantFile
try:
    import matplotlib
    matplotlib.use('agg')
    from matplotlib import pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    matplotlib.rcParams['hatch.linewidth'] = 0.1
    color_pool = plt.get_cmap('tab10').colors
except Exception as e:
    print(e)

import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame
import scipy.stats as stats

"""
%3A : (colon)
%3B ; (semicolon)
%3D = (equal sign)
%25 % (percent sign)
%2C , (comma)
%0D CR
%0A LF
%09 TAB
"""

def demo_one(target_vcfs:tuple):
    # 目标: 统计基因EGFR和KRAS在肺癌样本中出现somatic突变的频率
    # 获取肺癌样本cas_id
    # main_data_records = [MainData(), MainData(), ...]
    # vcfs = [VCF(), VCF(), ...]
    # lung_cancer_sample_case_id_list = [x.case_id for x in main_data_records if x.disease_type in ['Lung Adenocarcinoma', 'Lung Squamous Cell Carcinoma']]
    #
    # # 获取肺癌样本对应的somatic vcf
    # target_vcfs = [x for x in vcfs if (x.metadata.Tissue_type == "Tumor") and (x.metadata.Case_id in lung_cancer_sample_case_id_list)]

    # 搜索包含EGFR,KRAS突变的记录
    target_mutation_records = []
    for vcf in target_vcfs:
        sample = os.path.basename(vcf).split('.')[0]
        with VariantFile(vcf) as f:
            # get csq format
            csq_format = f.header.info['CSQ'].description.split('Format: ')[1]
            # parse line by line
            for r in f:
                if list(r.filter)[0] == "PASS":
                    for each in r.info['CSQ']:
                        record_dict = dict(zip(csq_format.split('|'), each.split('|')))
                        if record_dict['PICK'] == "1" and record_dict['SYMBOL'] in ["EGFR", "KRAS"]:
                            target_mutation_records.append([
                                sample,
                                record_dict['SYMBOL'],
                                record_dict['Gene'],
                                record_dict['CANONICAL'],
                                record_dict['HGVSc'].replace('%3D', '='),
                                record_dict['HGVSp'].replace('%3D', '='),
                                record_dict['VARIANT_CLASS'],
                                record_dict['IMPACT']
                            ])
                            break
    # 输出
    print(pd.DataFrame(target_mutation_records))  # 最后的输出表
    number_of_case_has_target_mutation = len(set([x[0] for x in target_mutation_records]))  # 计算case_id的uniq数量
    # freq = number_of_case_has_target_mutation / len(lung_cancer_sample_case_id_list)
    # print(freq)


def demo_two(target_vcfs: tuple):
    # 目标：统计在肺癌病人中，突变频率排名前10的是哪些基因, 每个基因对应的频次前5的具体突变是谁
    # 获取肺癌样本cas_id
    # main_data_records = []
    # vcfs = [VCF, VCF, ...]
    # lung_cancer_sample_case_id_list = [x.case_id for x in main_data_records if
    #                                    x.disease_type in ['Lung Adenocarcinoma', 'Lung Squamous Cell Carcinoma']]
    #
    # # 获取肺癌样本对应的somatic vcf
    # target_vcfs = [x for x in vcfs if
    #                (x.metadata.Tissue_type == "Tumor") and (x.metadata.Case_id in lung_cancer_sample_case_id_list)]

    # 统计基因的人群频率
    # mutation_records: {gene: {case_id: {HGVSc1, HGVSc2, ...}, gene2:...}
    CSQ_format = \
        """
        Allele|Consequence|IMPACT|SYMBOL|Gene|Feature
        _type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|
        FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|GENE_PHENO|NEAREST|SIFT|PolyPhen|HGVS_OFFSET|AFR_AF|
        AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_S
        AS_AF|MAX_AF|MAX_AF_POPS|FREQS|CLIN_SIG|SOMATIC|PHENO|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSeque
        nce|WildtypeProtein
        """
    mutation_records = dict()
    for vcf in target_vcfs:
        sample = os.path.basename(vcf).split('.')[0]
        with VariantFile(vcf) as f:
            # get csq format
            csq_format = f.header.info['CSQ'].description.split('Format: ')[1]
            # parse line by line
            for r in f:
                if list(r.filter)[0] == "PASS":
                    for each in r.info['CSQ']:
                        record_dict = dict(zip(csq_format.split('|'), each.split('|')))
                        if record_dict['PICK'] == "1":
                            tmp_dict = mutation_records.setdefault(record_dict['SYMBOL'], dict())
                            # case_id =  vcf.Case_id
                            case_id = sample
                            tmp_dict.setdefault(sample, set()).add(record_dict['HGVSc'].replace('%3D', '='))
                            break
    # statistic
    top10 = sorted([(k, len(v)) for k, v in mutation_records.items()], key=lambda x:x[1], reverse=True)[:10]
    table = []
    # total_case_number = len(lung_cancer_sample_case_id_list)
    total_case_number = len(target_vcfs)
    for gene, case_number in top10:
        hgvsc_lst = []
        _ = [hgvsc_lst.extend(list(x)) for x in mutation_records[gene].values()]
        # print(hgvsc_lst)
        top5_mutation = Counter(hgvsc_lst).most_common(5)
        print(top5_mutation)
        table.append([
            gene,  # Top10:Genes
            case_number,  # Case_number
            case_number/total_case_number,  # Frequency
            '|'.join(x[0] for x in top5_mutation),  # Top5_mutation
        ])
    # 输出
    df = pd.DataFrame(table)
    df.columns = ['Top10:Genes', 'Case_number', 'Frequency', 'Top5_mutation']
    print(df)


def demo_three():
    # 目标：统计在肠癌病人中，且肝转移灶数量大于5，突变频率排名前10的是哪些基因, 每个基因对应的频次前5的具体突变是谁
    # 实际和demo_two基本一样，只是在代码125行多加了一个条件：x.liver_leision >= 5
    # 获取肺癌样本cas_id
    main_data_records = []
    vcfs = [VCF, VCF, ...]
    cancer_sample_case_id_list = [x.case_id for x in main_data_records if x.disease_type in ['Rectum Adenocarcinoma', 'Colon Adenocarcinoma'] and x.liver_leision >= 5]

    # 获取肺癌样本对应的somatic vcf
    target_vcfs = [x for x in vcfs if
                   (x.metadata.Tissue_type == "Tumor") and (x.metadata.Case_id in cancer_sample_case_id_list)]

    # 统计基因的人群频率
    # mutation_records: {gene: {case_id_set: HGVSc_set}}
    mutation_records = dict()
    CSQ_format = \
        """
        Allele|Consequence|IMPACT|SYMBOL|Gene|Feature
        _type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|
        FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|GENE_PHENO|NEAREST|SIFT|PolyPhen|HGVS_OFFSET|AFR_AF|
        AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_S
        AS_AF|MAX_AF|MAX_AF_POPS|FREQS|CLIN_SIG|SOMATIC|PHENO|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSeque
        nce|WildtypeProtein
        """
    for vcf in target_vcfs:
        for line in vcf:
            if line.FILTER == "PASS":
                for each in line.INFO.CSQ.split(','):
                    record_dict = dict(zip(CSQ_format.split('|'), each.split('|')))
                    if record_dict['PICK'] == '1':
                        tmp_dict = mutation_records.setdefault(record_dict['SYMBOL'], dict())
                        tmp_dict.setdefault(vcf.Case_id, set()).add(record_dict['HGVSc'])
                        break
    # statistic
    top10 = sorted([(k, len(v)) for k, v in mutation_records.items()], key=lambda x:x[1], reverse=True)[:10]
    table = []
    for gene, case_number in top10:
        hgvsc_set = set()
        _ = [hgvsc_set.update(x) for x in mutation_records[gene].values()]
        top5_mutation = Counter(hgvsc_set).most_common(5)
        table.append([
            gene,  # Top10:Genes
            case_number,  # Case_number
            case_number/len(cancer_sample_case_id_list),  # Frequency
            top5_mutation,  # Top5_mutation
        ])
    # 输出
    print(table)


def demo_four():
    # 肠癌人群中，哪些病人接受过化疗，化疗的次数分布如何
    # 获取肺癌样本cas_id
    main_data_records = []

    cases_durations = []
    for x in main_data_records:
        condition = x.disease_type in ['Rectum Adenocarcinoma', 'Colon Adenocarcinoma']
        condition2 = x.chemotherapy == "是"
        if condition and condition2:
            duration_num = sum([int(x) for x in x.treatment_duration.split(';')])
            cases_durations.append([x.case_id, duration_num])
    print(cases_durations)


def get_top_mutated_genes(vcfs:tuple, top=20, tumor_index=None):
    gene_mutations = dict()
    mutation_af_dict = dict()
    mutation_hgvsp_dict = dict()
    gene_symbol_dict = dict()
    # {gene: {sample: set(muations)}}
    for vcf in vcfs:
        # sample = os.path.basename(vcf).split('.')[0]
        tumor_idx = guess_tumor_idx(vcf) if tumor_index is None else tumor_index
        with VariantFile(vcf) as f:
            sample = f.header.samples[tumor_idx]
            # get csq format
            csq_format = f.header.info['CSQ'].description.split('Format: ')[1]
            # parse line by line
            for r in f:
                if list(r.filter)[0] == "PASS":
                    alt_af_dict = dict(zip(format_alts(r), r.samples[tumor_idx]['AF']))
                    for each in r.info['CSQ']:
                        record_dict = dict(zip(csq_format.split('|'), each.split('|')))
                        if record_dict['PICK'] == "1":
                            if not record_dict['Gene']:
                                # print(f'Off gene mutation:{record_dict["HGVSc"]}')
                                continue
                            tmp_dict = gene_mutations.setdefault(record_dict['Gene'], dict())
                            mut_id = record_dict['HGVSc'].replace('%3D', '=')
                            mutation_hgvsp_dict[mut_id] = record_dict['HGVSp'].replace('%3D', '=')
                            try:
                                mutation_af_dict.setdefault(mut_id, []).append(alt_af_dict[record_dict['Allele']])
                            except:
                                print(alt_af_dict, r.alts)
                            gene_symbol_dict[record_dict['Gene']] = record_dict['SYMBOL']
                            tmp_dict.setdefault(sample, set()).add(mut_id)
                            break
    # statistic: 按照共有突变出现频次排序
    top_n = sorted([(k, len(v)) for k, v in gene_mutations.items()], key=lambda x: x[1], reverse=True)[:top]
    table = []
    # total_case_number = len(lung_cancer_sample_case_id_list)
    total_case_number = len(vcfs)
    for gene, case_number in top_n:
        hgvsc_lst = []
        _ = [hgvsc_lst.extend(list(x)) for x in gene_mutations[gene].values()]
        gene_symbol = gene_symbol_dict[gene]
        # print(hgvsc_lst)
        top5_mutation = Counter(hgvsc_lst).most_common(5)
        table.append([
            gene,  # Top10:Genes gene_id
            gene_symbol,  # symbol
            case_number,  # Case_number
            case_number / total_case_number,  # Frequency
            '|'.join(x[0] or "None" for x in top5_mutation),  # Top5_mutation
            '|'.join(f'{mutation_hgvsp_dict[x[0]] or "None"}' for x in top5_mutation),  # Top5_mutation hgvsp
            '|'.join(str(x[1]) for x in top5_mutation),  # Top5_mutation freq
            '|'.join(f'{sts.mean(mutation_af_dict[x[0]]):.3f}' for x in top5_mutation),  # Top5_mutation vaf
            '|'.join(gene_mutations[gene].keys())
        ])
    # 输出
    df = pd.DataFrame(table)
    df.columns = [f'Top{top}:Genes', 'Symbol', 'Case_number', 'Frequency', 'Top5_mut_hgvsc', 'Top5_mut_hgvsp',
                  'Top5_mut_freq', 'Top5_mut_meanAF', 'samples']
    df.to_csv(f'top{top}.mutated_genes.txt', sep='\t', index=False)
    print(df)


def filter_by_pick_flag(vcf, out=None):
    tumor_idx = guess_tumor_idx(vcf)
    with VariantFile(vcf) as fr:
        sample = fr.header.samples[tumor_idx]
        out = out if out else f'{sample}.pick1.vcf.gz'
        with VariantFile(out, 'w', header=fr.header) as fw:
            # get csq format
            csq_format = fr.header.info['CSQ'].description.split('Format: ')[1]
            # parse line by line
            for r in fr:
                if list(r.filter)[0] == "PASS":
                    for each in r.info['CSQ']:
                        record_dict = dict(zip(csq_format.split('|'), each.split('|')))
                        if record_dict['PICK'] == "1":
                            break
                    r.info['CSQ'] = [each]
                    fw.write(r)


def guess_tumor_idx(vcf_file, require_pass=True):
    tumor_is_first = 1
    tumor_is_second = 1
    with VariantFile(vcf_file, ignore_truncation=True) as fr:
        samples = list(fr.header.samples)
        try:
            tumor_sample = [(x.key, x.value) for x in fr.header.records if x.key == 'tumor_sample'][0][1]
            record_tumor_sample_idx = samples.index(tumor_sample)
        except Exception as e:
            tumor_sample = None
            record_tumor_sample_idx = None
            pass

        if len(samples) == 1:
            tumor_idx = 0
        else:
            formats = list(fr.header.formats)
            if 'AF' not in formats:
                raise Exception('No AF in <Format> info to detect tumor sample')
            for record in fr:
                if list(record.filter)[0] != "PASS" and require_pass:
                    continue
                try:
                    if record.samples[0]['AF'][0] > record.samples[1]['AF'][0]:
                        tumor_is_first += 1
                    else:
                        tumor_is_second += 1
                except Exception as e:
                    print(vcf_file, e)

            tumor_idx = tumor_is_second >= tumor_is_first
    print(f'we guess tumor sample is {int(tumor_idx) + 1}th {samples[tumor_idx]} base on ratio {tumor_is_first} vs {tumor_is_second}')
    if not(tumor_is_first/tumor_is_second >= 1.5 or tumor_is_second/tumor_is_first >= 1.5):
        print('The above guessing is not reliable since the ratio < 1.5')

    if tumor_sample is not None:
        if record_tumor_sample_idx != tumor_idx:
            print(f'Header info indicates tumor sample is {tumor_sample} while we guess tumor sample is {samples[tumor_idx]}')
            tumor_idx = record_tumor_sample_idx
        else:
            print('our guessing of tumor sample is consistent with the record of header')
    return tumor_idx


def format_alts(r):
    """
    把vcf格式的alt转换为txt格式的alt
    snp: A > T
    insertion: A > AGG; A > AAA
    deletion: AG > A; AGT > A
    substitution（same length between ref and alt): AT > TG
    indel: AC > TAA
    """
    alts = list(r.alts)
    for ind, each in enumerate(r.alts):
        if len(each) != len(r.ref):
            if len(r.ref) < len(each) and each.startswith(r.ref):
                # insertion
                alts[ind] = each.split(r.ref, 1)[1]
            elif len(r.ref) > len(each) and r.ref.startswith(each):
                # deletion
                alts[ind] = '-'
    return alts


def get_tsg(tsg_file, value_type='Ensembl'):
    """
    :param tsg_file: from https://bioinfo.uth.edu/TSGene/
    :return:
    """
    table = pd.read_table(tsg_file, index_col=0)
    targets = set()
    for x in table['Links']:
        if value_type + ':' in x:
            targets.add(x.split(value_type+':', 1)[1].split('|', 1)[0])
        else:
            print(f'Ignore tsg for not match value_type {value_type}', x)
    return targets


def merge_vcf_as_table(vcfs:tuple, out, min_af=0.001, min_alt_depth=2, min_depth=20, max_pop_freq=1e-6):
    """
    vcf "AD" style = [ref_depth, alt1_depth, alt2_depth]
    :param vcfs:
    :param out:
    :param min_af:
    :param min_alt_depth:
    :param min_depth:
    :param max_pop_freq:
    :return:
    """
    results = []
    samples = []
    for vcf_file in vcfs:
        tumor_idx = guess_tumor_idx(vcf_file)
        with VariantFile(vcf_file, ignore_truncation=True) as fr:
            # get csq format
            sample = fr.header.samples[tumor_idx]
            if sample in samples:
                raise Exception(f'we find duplicated sample {sample}')
            else:
                samples.append(sample)
            csq_header = fr.header.info['CSQ'].description.split('Format: ')[1]

            # parse line by line
            for r in fr:
                alt_depth_dict = dict(zip(format_alts(r), r.samples[tumor_idx]['AD'][1:]))
                # pass and depth filter
                dp = r.samples[tumor_idx]['DP']
                if list(r.filter)[0] != "PASS" or (dp < min_depth):
                    continue

                # pick 1 filter
                for each in r.info['CSQ']:
                    csq_dict = dict(zip(csq_header.split('|'), each.split('|')))
                    if csq_dict['PICK'] == "1":
                        break

                # ad and af filter
                ad = alt_depth_dict[csq_dict['Allele']]
                af = ad / dp
                if (af < min_af) or (ad < min_alt_depth):
                    continue

                # population af filter
                if 'MAX_AF' in csq_dict:
                    if csq_dict['MAX_AF'] and (float(csq_dict['MAX_AF']) > max_pop_freq):
                        continue
                # save record
                target_info = [
                    sample,
                    r.contig,
                    r.pos,
                    r.ref,
                    # '|'.join(r.alts),
                    csq_dict['Allele'],
                    round(af, 3),
                    ad,
                    dp,
                    csq_dict['SYMBOL'],
                    csq_dict['Gene'],
                    csq_dict['EXON'],
                    csq_dict['HGVSc'].replace('%3D', '='),
                    csq_dict['HGVSp'].replace('%3D', '='),
                    csq_dict['Consequence'],
                    csq_dict['VARIANT_CLASS'],
                    csq_dict['CANONICAL'],
                    csq_dict['Existing_variation'],
                    csq_dict['SOMATIC']
                ]
                if 'MAX_AF' in csq_dict:
                    target_info.append(csq_dict['MAX_AF'])
                else:
                    target_info.append('')
                results.append(target_info)
    header = ['sample', 'contig', 'position', 'ref', 'alt',
              'AF', 'AD', 'DP', 'SYMBOL', 'Gene',
              'EXON', 'HGVSc', 'HGVSp', 'Consequence',
              'VARIANT_CLASS', 'CANONICAL', 'Existing_variation', 'SOMATIC', 'MAX_POP_AF']
    df = pd.DataFrame(results, columns=header)
    df.to_csv(out, sep='\t', index=False)
    print(df.head())
    if len(set(df['sample'])) >= 2:
        # af boxplot
        ax = sns.boxplot(data=df, x='sample', y='AF')
        if len(set(df['sample'])) > 6:
            ax.xaxis.set_tick_params(labelsize=7, labelrotation=90)
        plt.savefig('AF.boxplot.pdf', bbox_inches="tight")
        plt.close()
        # af density plot
        sns.kdeplot(data=df, x='AF', hue='sample')
        plt.minorticks_on()
        plt.savefig('AF.density.pdf', bbox_inches="tight")
        plt.close()
        # variant_class stack bar
        target_cols = ['sample', 'VARIANT_CLASS', 'AF']
        df2 = df[target_cols].groupby(['sample', 'VARIANT_CLASS']).count()
        df3 = df2.reset_index().pivot(index='sample', columns='VARIANT_CLASS', values='AF')
        df3.plot.bar(stacked=True)
        plt.ylabel('Variant count')
        plt.tick_params(labelsize=7)
        plt.savefig('variant_class_count.boxplot.pdf', bbox_inches="tight")
        plt.close()
        # variant_class percent stack bar
        df4 = df3.div(df3.sum(axis=1), axis=0)
        df4.plot.bar(stacked=True)
        plt.tick_params(labelsize=7)
        plt.ylabel('Variant percent')
        plt.savefig('variant_class_percent.boxplot.pdf', bbox_inches="tight")
        plt.close()


def get_TMB(vcfs:tuple, out='TMB.txt', bed_size=59464418, tumor_index=None,
            min_af=0.05, min_alt_depth=3, min_depth=25, max_pop_freq=0.01,
            pick=True, tsg_file=None, include_synonymous=False, csq_tag='CSQ'):

    """
    https://zhuanlan.zhihu.com/p/378646738
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6510391/

    常见免疫治疗标志物：
    目前免疫治疗最为确定的 3 个标志物——TMB，MSI 和 PD-L1 表达常有不同：
    同时具备3个标志物的患者比例仅为 0.6%。
    MSI-H 和 TMB-H 似乎具有一定的相关性，83% 的 MSI-H 患者同时为 TMB-H，而只有 16% 的 TMB-H 是 MSI-H。
    TMB 和 PD-L1 则相互独立，同时 TMB-H 和 PD-L1 阳性者对 ICI 反应最佳。

    TMB的定义：
    TMB of a tumor sample is calculated by the number of non-synonymous somatic mutations (single nucleotide variants and small insertions/deletions) per mega-base in coding regions
    TMB 是指特定基因组编码区域内非同义体细胞突变的个数，通常用每兆碱基多少个突变表示（mut/Mb），在早期研究中也直接以突变数量表示。TMB 可以间接反映肿瘤产生新抗原的能力和程度，预测多种肿瘤的免疫治疗疗效。

    TMB的检测：
    TMB 的检测受样本质量、检测方法和分析方法等多种因素影响，临床应用前应充分了解 TMB 检测的条件。
    样本：肿瘤纯度要求 ≥ 20%，需要有正常对照为检测提供胚系变异信息。
    方法：全外显子测序（WES）是 TMB 检测的金标准。但是 WES 价格昂贵，检测时间长，需要新鲜标本，因而应用受限。靶向测序 panel 已经成为 WES 的有效替代，为准确性考虑，其覆盖范围应 ≥ 1.0 Mb，测序深度 ≥ 500×。
    分析：TMB 的中位值和分布范围在不同癌种中有所不同，因此，在各个癌种中分别确定界值十分重要。应使用相同的筛选策略，选择排序在 20% 以上的病例定义为 TMB-H，而前瞻性的临床疗效才是确定 TMB 界值的最佳标准。不同靶向测序 panel 的 TMB 不能通用。

    TMB的对于治疗的意义
    对于无标准治疗的晚期肿瘤患者，TMB-H 提供了免疫治疗获益的可能。 某些情况下 TMB 可预测 ICI 治疗反应，但结论并不一致，特别是 TMB 预测长期结局以及免疫联合治疗的疗效时应慎重。
    bTMB 和 tTMB 具有一致性，但临床应用还缺乏强证据。 使用 TMB 时，应结合瘤种、人口特征、基因特征和检测方法综合解读。联合使用 PD-L1 和 TMB 等多种生物标志物可能是筛选免疫获益人群的更好方法。

    胃癌高通量测序（NGS）临床应用中国专家共识（2023版）之 胃癌TMB-H的判定标准：
    TMB-H也是ICIs治疗获益肿瘤患者的分子标志物，与MSI-H具有显著相关性[31]。与MSI状态鉴定不同，TMB检测的“金标准”是WES[32]。
    胃癌患者中TMB值≤5 muts/Mb的约占66.2%，5～10 muts/Mb约11.3%，10～40 muts/Mb约19.5%，≥40 muts/Mb约3.0%[33]。
    一项泛癌种研究发现，TMB前20%患者接受ICIs治疗的OS更长；采用“五分位法”定义TMB-H（>8.8 muts/Mb），其中126例食管胃腺癌患者未观察到生存获益[34]。
    特瑞普利单抗临床试验发现，前20%的TMB-H（TMB≥12 muts/Mb）胃癌患者具有更优的客观缓解率和OS[35]。2022 V2版NCCN胃癌指南推荐使用帕博利珠单药治疗TMB截断值≥10 muts/Mb的患者[36]。
    通常检测的基因数目越多，NGS panel拟合的TMB值与WES的一致性越好，使用超过300个基因的大panel（覆盖率≥0.8 Mb，最好≥1.0 Mb）即可较好地评估TMB值[37]。
    选择合适的NGS方法（>300个基因的panel）判断TMB-H，有助于作为医生对胃癌患者采用ICIs治疗的参考（推荐等级：Ⅱ级）。

    TMB 的计算
    TMB was defined as the number of somatic mutations in the coding region per megabase, including single nucleotide variants (SNVs) and small INDELs (insertions and deletions, usually less than 20 bases).
    However, the means to determine reliable somatic coding mutations for TMB calculation was not trivial.
    In the F1CDx approach, synonymous mutations were included. However, stop-gain mutations in tumor suppressor genes and hotspot driver mutations were not included in order to reduce bias due to enrichment of cancer-related genes in the F1CDx panel.
    Here, we followed the F1CDx approach to calculate the TMB for our panel, defining the cutoff values as TMB-high (≥20 mutations/Mb), TMB-medium (<20 mutations/Mb ≥10 mutations/Mb) and TMB-low (<10 mutations/Mb).
    Others filtering parameters of somatic mutations identified by WES for TMB calculation were a mutated allele frequency greater than 5% and a sequence depth greater than 20X in tumor samples greater or 10X in normal samples.

    MAF文件中包含的突变类型，对于哪些突变应纳入计算目前似乎没有标准
    The MAF file contains 16 types of somatic mutations (Missense_Mutation, Silent, Nonsense_Mutation, Intron, 3'UTR, 5'UTR, Splice_Site, RNA,
    Frame_Shift_Ins, Frame_Shift_Del, In_Frame_Ins, Nonstop_Mutation, In_Frame_Del, 3'Flank, 5'Flank, Translation_Start_Site) flagged by variant calling software packages.
    Since TMB is defined as the mutations in the coding region in most studies, four types of mutations (Intron, RNA, 3'Flank, 5'Flank) outside the coding region are excluded.
    Considering the mutation types are different for the two FDA’s approval assays F1CDx (total point mutations in the coding region) and MSK-IMPACT (nonsynonymous mutations in the coding region), both of which are included in our study.

    根据参考文献，看看TMB的具体计算细节：
    一、根据文献，https://jitc.bmj.com/content/8/1/e000147（Establishing guidelines to harmonize tumor mutational burden (TMB):
    in silico assessment of variation in TMB quantification across diagnostic platforms: phase I of the Friends of Cancer Research TMB Harmonization Project）
    1. 分析区域选择：CCDS区域，计算区域大小时，记得去冗余
    2. 质控标准： DP>=25, AD>=3, AF>=5%; 如果一个样本中的50%突变没有通过QC，则样本去掉
    3. 保留以下6类型变异：
        Frame_Shift_Del
        Frame_Shift_Ins
        In_Frame_Del
        In_Frame_Ins
        Missense_Mutation
        Nonsense_Mutation

    二、肿瘤突变负荷检测国家参考品（v2.0)
    1. 分析区域选择：未详细说明
    2. 质控标准： DP>=25, AD>=3, AF>=5%; 如果一个样本中的50%突变没有通过QC，则样本去掉
    3. 对于WES，保留以下6类型变异：
        Frame_Shift_Del
        Frame_Shift_Ins
        In_Frame_Del
        In_Frame_Ins
        Missense_Mutation
        Splice_Site（注意，这里和上面的文章不一致）
    4. 肿瘤突变负荷检测及临床应用中国专家共识（2020 年版）: WES的平均测序深度需要大于100x

    其他文献：https://www.nature.com/articles/s41698-024-00504-1
    Here we show that the panel sizes beyond 1.04Mb and 389 genes are necessary for the basic discrete accuracy, as determined by over 40,000 synthetic panels.
    The somatic mutation detection should maintain a reciprocal gap of recall and precision less than 0.179 for reliable psTMB calculation results.
    The inclusion of synonymous, nonsense and hotspot mutations could enhance the accuracy of panel-based TMB assay.
    A 5% variant allele frequency cut-off is suitable for TMB assays using tumor samples with at least 20% tumor purity.
    该文章中的WES参考的是上述第一种计算方式：

    看看VEP给出的突变类型（对应的其实是consequence）https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html
    SO term	SO description	SO accession	Display term	IMPACT
    transcript_ablation	A feature ablation whereby the deleted region includes a transcript feature	SO:0001893	Transcript ablation	HIGH
    splice_acceptor_variant	A splice variant that changes the 2 base region at the 3' end of an intron	SO:0001574	Splice acceptor variant	HIGH
    splice_donor_variant	A splice variant that changes the 2 base region at the 5' end of an intron	SO:0001575	Splice donor variant	HIGH
    stop_gained	A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript	SO:0001587	Stop gained	HIGH
    frameshift_variant	A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three	SO:0001589	Frameshift variant	HIGH
    stop_lost	A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript	SO:0001578	Stop lost	HIGH
    start_lost	A codon variant that changes at least one base of the canonical start codon	SO:0002012	Start lost	HIGH
    transcript_amplification	A feature amplification of a region containing a transcript	SO:0001889	Transcript amplification	HIGH
    feature_elongation	A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence	SO:0001907	Feature elongation	HIGH
    feature_truncation	A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence	SO:0001906	Feature truncation	HIGH
    inframe_insertion	An inframe non synonymous variant that inserts bases into in the coding sequence	SO:0001821	Inframe insertion	MODERATE
    inframe_deletion	An inframe non synonymous variant that deletes bases from the coding sequence	SO:0001822	Inframe deletion	MODERATE
    missense_variant	A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved	SO:0001583	Missense variant	MODERATE
    protein_altering_variant	A sequence_variant which is predicted to change the protein encoded in the coding sequence	SO:0001818	Protein altering variant	MODERATE
    splice_donor_5th_base_variant	A sequence variant that causes a change at the 5th base pair after the start of the intron in the orientation of the transcript	SO:0001787	Splice donor 5th base variant	LOW
    splice_region_variant	A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron	SO:0001630	Splice region variant	LOW
    splice_donor_region_variant	A sequence variant that falls in the region between the 3rd and 6th base after splice junction (5' end of intron)	SO:0002170	Splice donor region variant	LOW
    splice_polypyrimidine_tract_variant	A sequence variant that falls in the polypyrimidine tract at 3' end of intron between 17 and 3 bases from the end (acceptor -3 to acceptor -17)	SO:0002169	Splice polypyrimidine tract variant	LOW
    incomplete_terminal_codon_variant	A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed	SO:0001626	Incomplete terminal codon variant	LOW
    start_retained_variant	A sequence variant where at least one base in the start codon is changed, but the start remains	SO:0002019	Start retained variant	LOW
    stop_retained_variant	A sequence variant where at least one base in the terminator codon is changed, but the terminator remains	SO:0001567	Stop retained variant	LOW
    synonymous_variant	A sequence variant where there is no resulting change to the encoded amino acid	SO:0001819	Synonymous variant	LOW
    coding_sequence_variant	A sequence variant that changes the coding sequence	SO:0001580	Coding sequence variant	MODIFIER
    mature_miRNA_variant	A transcript variant located with the sequence of the mature miRNA	SO:0001620	Mature miRNA variant	MODIFIER
    5_prime_UTR_variant	A UTR variant of the 5' UTR	SO:0001623	5 prime UTR variant	MODIFIER
    3_prime_UTR_variant	A UTR variant of the 3' UTR	SO:0001624	3 prime UTR variant	MODIFIER
    non_coding_transcript_exon_variant	A sequence variant that changes non-coding exon sequence in a non-coding transcript	SO:0001792	Non coding transcript exon variant	MODIFIER
    intron_variant	A transcript variant occurring within an intron	SO:0001627	Intron variant	MODIFIER
    NMD_transcript_variant	A variant in a transcript that is the target of NMD	SO:0001621	NMD transcript variant	MODIFIER
    non_coding_transcript_variant	A transcript variant of a non coding RNA gene	SO:0001619	Non coding transcript variant	MODIFIER
    coding_transcript_variant	A transcript variant of a protein coding gene	SO:0001968	Coding transcript variant	MODIFIER
    upstream_gene_variant	A sequence variant located 5' of a gene	SO:0001631	Upstream gene variant	MODIFIER
    downstream_gene_variant	A sequence variant located 3' of a gene	SO:0001632	Downstream gene variant	MODIFIER
    TFBS_ablation	A feature ablation whereby the deleted region includes a transcription factor binding site	SO:0001895	TFBS ablation	MODIFIER
    TFBS_amplification	A feature amplification of a region containing a transcription factor binding site	SO:0001892	TFBS amplification	MODIFIER
    TF_binding_site_variant	A sequence variant located within a transcription factor binding site	SO:0001782	TF binding site variant	MODIFIER
    regulatory_region_ablation	A feature ablation whereby the deleted region includes a regulatory region	SO:0001894	Regulatory region ablation	MODIFIER
    regulatory_region_amplification	A feature amplification of a region containing a regulatory region	SO:0001891	Regulatory region amplification	MODIFIER
    regulatory_region_variant	A sequence variant located within a regulatory region	SO:0001566	Regulatory region variant	MODIFIER
    intergenic_variant	A sequence variant located in the intergenic region, between genes	SO:0001628	Intergenic variant	MODIFIER
    sequence_variant	A sequence_variant is a non exact copy of a sequence_feature or genome exhibiting one or more sequence_alteration	SO:0001060	Sequence variant	MODIFIER

    根据理解，依据VEP的注释，可以纳入以下consequence的突变用于WES的TMB计算：
    frameshift_variant
    missense_variant
    stop_lost
    start_lost
    splice_acceptor_variant
    splice_donor_variant
    stop_gained
    inframe_insertion
    inframe_deletion

    :param vcfs: 一个或多个VCF文件
    :param out: 输出的csv文件名称
    :param bed_size: panel的大小,单位是bp，这个大小是否要考虑内含子区域在内？
    :param tumor_index: 从0开始计数，指示vcf中第几个样本是肿瘤样本，默认自己根据AF高低统计判断肿瘤样本
    :param min_af: 如果某个突变的频率AF小于该值，则不纳入统计范围
    :param min_alt_depth: 如果某个突变的支持reads数量小于该值，则不纳入统计范围
    :param min_depth: 如果某个突变所在位置总测序深度小于该值，则不纳入统计范围
    :param max_pop_freq: 如果某个突变的人群频率（vep给出的MAX_AF）大于该值，则不纳入统计范围
    :param pick: 是否只选择vep给出的pick flag为1的突变进行分析。这关系到突变的意义解读。
        VEP中“--pick” | “--flag_pick”的含义
        Pick one line or block of consequence data per variant, including transcript-specific columns.
        Consequences are chosen according to the criteria described here, and the order the criteria are applied may be customised with --pick_order.
        This is the best method to use if you are interested only in one consequence per variant
    :param tsg_file: 肿瘤抑制子基因文件https://bioinfo.uth.edu/TSGene/，默认选择Ensembl来源的基因作为候选。如果提供该文件，则针对这些基因进行判断，如果consequence为['stop_gained', 'start_lost']则不纳入TMB统计
    :param include_synonymous: 如果为False，则TMB统计时，不会纳入同义突变，根据VEP给出的consequence判断
    :param csq_tag: 指示VCF文件中，哪个字段为VEP的注释结果，默认为CSQ
    :return:
    """
    # sample = os.path.basename(vcf).split('.')[0]
    tsg = get_tsg(tsg_file, value_type='Ensembl') if tsg_file else set()
    count_lst = []
    tmb_lst = []
    sample_lst = []
    target_consequences = [
        "frameshift_variant",
        "missense_variant",
        "inframe_insertion",
        "inframe_deletion",
        "stop_lost",
        "start_lost",
        "stop_gained",
        # "splice_acceptor_variant",
        # "splice_donor_variant",
    ]

    synonymous_consequences = [
        "synonymous_variant",
        "start_retained_variant",
        "stop_retained_variant"
    ]

    for vcf in vcfs:
        count = 0
        tumor_idx = guess_tumor_idx(vcf) if tumor_index is None else tumor_index
        with VariantFile(vcf) as f:
            # get csq format
            csq_format = f.header.info['CSQ'].description.split('Format: ')[1]
            # parse line by line
            sample = f.header.samples[tumor_idx]
            for r in f:
                alt_depth_dict = dict(zip(format_alts(r), r.samples[tumor_idx]['AD'][1:]))
                dp = r.samples[tumor_idx]['DP']
                if list(r.filter)[0] != "PASS" or (dp < min_depth):
                    continue
                include = False
                for each in r.info[csq_tag]:
                    csq_dict = dict(zip(csq_format.split('|'), each.split('|')))
                    symbol = csq_dict['SYMBOL']
                    gene = csq_dict['Gene']

                    # 仅仅选择一条注释进行解读
                    if pick and csq_dict['PICK'] != "1":
                        continue

                    # population frequency filter
                    if 'MAX_AF' in csq_dict:
                        if csq_dict['MAX_AF'] and (float(csq_dict['MAX_AF']) >= max_pop_freq):
                            # print('max pop af', csq_dict['MAX_AF'], csq_dict['MAX_AF_POPS'])
                            continue

                    # ad and af
                    ad = alt_depth_dict[csq_dict['Allele']]
                    af = ad / dp
                    if (af < min_af) or (ad < min_alt_depth):
                        continue

                    # if ignore serious tumor suppressor gene mutation
                    if (symbol in tsg) or (gene in tsg):
                        if csq_dict['Consequence'] in ['stop_gained', 'start_lost']:
                            continue

                    # if include synonymous mutation
                    if not include_synonymous:
                        if csq_dict['Consequence'] not in target_consequences:
                            continue
                    else:
                        if csq_dict['Consequence'] not in target_consequences + synonymous_consequences:
                            continue

                    # pass all filter
                    include = True
                    print(csq_dict)
                    break

                if include:
                    count += 1
        tmb_value = count/bed_size*1e6
        print(f'TMB value of {sample}: {count}/{bed_size} = {tmb_value:.2f} mutations/Mb')
        count_lst.append(count)
        tmb_lst.append(tmb_value)
        sample_lst.append(sample)

    df = pd.DataFrame({'sample': sample_lst, 'variant_count': count_lst, 'TMB': tmb_lst})
    df.to_csv(out, index=False, sep='\t')
    return sample_lst, count_lst, tmb_lst


def merge_vcf_as_maf(vcfs:tuple, out, min_af=0.05, min_alt_depth=2, min_depth=15, max_pop_freq=0.01):
    """
    注意输出的Ref和alt格式可能不符合maf格式要求
    vcf "AD" style = [ref_depth, alt1_depth, alt2_depth]
    :param vcfs:
    :param out:
    :param min_af:
    :param min_alt_depth:
    :param min_depth:
    :param max_pop_freq:
    :return:
    """
    print('注意:本程序输出的Ref和alt格式可能不符合maf格式要求')
    results = []
    samples = []
    for vcf_file in vcfs:
        tumor_idx = guess_tumor_idx(vcf_file)
        with VariantFile(vcf_file, ignore_truncation=True) as fr:
            # get csq format
            sample = fr.header.samples[tumor_idx]
            normal_sample = fr.header.samples[1-tumor_idx]
            if sample in samples:
                raise Exception(f'we find duplicated sample {sample}')
            else:
                samples.append(sample)
            csq_header = fr.header.info['CSQ'].description.split('Format: ')[1]
            genome_file = [(x.key, x.value) for x in fr.header.records if x.key == 'reference'][0][1]
            # print(genome_file)
            genome_file = os.path.basename(genome_file)

            # parse line by line
            for r in fr:
                alt_depth_dict = dict(zip(format_alts(r), r.samples[tumor_idx]['AD'][1:]))
                # pass and depth filter
                dp = r.samples[tumor_idx]['DP']
                if list(r.filter)[0] != "PASS" or (dp < min_depth):
                    continue

                # pick 1 filter
                for each in r.info['CSQ']:
                    csq_dict = dict(zip(csq_header.split('|'), each.split('|')))
                    if csq_dict['PICK'] == "1":
                        break

                # ad and af filter
                ad = alt_depth_dict[csq_dict['Allele']]
                af = ad / dp
                if (af < min_af) or (ad < min_alt_depth):
                    continue

                # population af filter
                if 'MAX_AF' in csq_dict:
                    if csq_dict['MAX_AF'] and (float(csq_dict['MAX_AF']) > max_pop_freq):
                        continue
                # save record
                target_info = [
                    csq_dict['SYMBOL'],
                    csq_dict['Gene'],
                    'CenterUnknown',
                    genome_file,
                    r.contig,
                    r.start+1, # 坐标是否和maf格式要求一致？
                    r.stop,
                    csq_dict['STRAND'],
                    csq_dict['Consequence'],
                    csq_dict['VARIANT_CLASS'],
                    r.ref,  # 和csq_dict['Allele']的格式不对应？
                    'unknown',
                    csq_dict['Allele'],
                    csq_dict['Existing_variation'],
                    'unknown',
                    sample,
                    normal_sample,
                    # 后续为自定义
                    round(af, 3),
                    ad,
                    dp,
                    csq_dict['EXON'],
                    csq_dict['HGVSc'].replace('%3D', '='),
                    csq_dict['HGVSp'].replace('%3D', '='),
                    csq_dict['CANONICAL'],
                    csq_dict['SOMATIC']
                ]
                if 'MAX_AF' in csq_dict:
                    target_info.append(csq_dict['MAX_AF'])
                else:
                    target_info.append('')
                results.append(target_info)
    header = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome',
              'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type',
              'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
              'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
              'AF', 'AD', 'DP', 'EXON', 'HGVSc', 'HGVSp', 'CANONICAL', 'SOMATIC', 'MAX_POP_AF']
    df = pd.DataFrame(results, columns=header)
    df.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
