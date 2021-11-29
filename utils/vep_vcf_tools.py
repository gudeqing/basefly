import os
from collections import Counter
from pysam import VariantFile
import matplotlib
# matplotlib.use('agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame
import scipy.stats as stats
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['hatch.linewidth'] = 0.1
color_pool = plt.get_cmap('tab10').colors

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


def pick_flag(vcf, out):
    with VariantFile(vcf) as fr:
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


def guess_tumor_idx(vcf_file):
    tumor_is_first = 0
    tumor_is_second = 0

    with VariantFile(vcf_file, ignore_truncation=True) as fr:
        samples = list(fr.header.samples)
        if len(samples) == 1:
            tumor_idx = 0
        else:
            formats = list(fr.header.formats)
            if 'AF' not in formats:
                raise Exception('No AF in format info to detect tumor sample')
            for record in fr:
                try:
                    if record.samples[0]['AF'][0] > record.samples[1]['AF'][0]:
                        tumor_is_first += 1
                    else:
                        tumor_is_second += 1
                except Exception as e:
                    print(vcf_file, e)

            tumor_idx = tumor_is_second >= tumor_is_first
    print(f'we guess tumor sample is {samples[tumor_idx]} ')
    return tumor_idx


def format_alts(r):
    alts = list(r.alts)
    """
    snp: A > T
    insertion: A > AGG; A > AAA
    deletion: AG > A; AGT > A
    substitution（same length between ref and alt): AT > TG
    indel: AC > TAA
    """
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



def merge_vcf_as_table(vcfs:tuple, out, min_af=0.001, min_alt_depth=3, min_depth=20, max_pop_freq=1e-6):
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
        sns.boxplot(data=df, x='sample', y='AF')
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
        plt.savefig('variant_class_count.boxplot.pdf', bbox_inches="tight")
        plt.close()
        # variant_class percent stack bar
        df4 = df3.div(df3.sum(axis=1), axis=0)
        df4.plot.bar(stacked=True)
        plt.ylabel('Variant percent')
        plt.savefig('variant_class_percent.boxplot.pdf', bbox_inches="tight")
        plt.close()



def get_tmb(vcfs:tuple, bed_size=34, tumor_idx=None, min_af=0.05, min_alt_depth=3, min_depth=15, max_pop_freq=1e-3,
            pick=True, tsg_file=None, synonymous=False, tag='CSQ'):
    # sample = os.path.basename(vcf).split('.')[0]
    tsg = get_tsg(tsg_file, value_type='Ensembl') if tsg_file else set()
    count_lst = []
    tmb_lst = []
    sample_lst = []
    for vcf in vcfs:
        count = 0
        tumor_idx = guess_tumor_idx(vcf) if tumor_idx is None else tumor_idx
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
                for each in r.info[tag]:
                    csq_dict = dict(zip(csq_format.split('|'), each.split('|')))
                    symbol = csq_dict['SYMBOL']
                    gene = csq_dict['Gene']
                    if pick and csq_dict['PICK'] != "1":
                        continue
                    # ad and af
                    ad = alt_depth_dict[csq_dict['Allele']]
                    af = ad / dp
                    if (af < min_af) or (ad < min_alt_depth):
                        continue
                    # ignore serious tumor suppressor gene mutation
                    if symbol in tsg or gene in tsg:
                        if csq_dict['Consequence'] in ['stop_gained', 'start_lost']:
                            continue
                    # ignore synonymous mutation
                    if not synonymous:
                        if csq_dict['Consequence'] in ['synonymous_variant']:
                            continue
                    # population frequency filter
                    if 'MAX_AF' in csq_dict:
                        if csq_dict['MAX_AF'] and (float(csq_dict['MAX_AF']) > max_pop_freq):
                            print('max pop af', csq_dict['MAX_AF'], csq_dict['MAX_AF_POPS'])
                            continue
                    if not csq_dict['EXON']:
                        continue
                    # skip out of loop once hit
                    include = True
                    break

                if include:
                    count += 1
        tmb_value = count/bed_size
        print(f'TMB value: {count}/{bed_size} = {tmb_value:.2f} variant/M')
        count_lst.append(count)
        tmb_lst.append(tmb_value)
        sample_lst.append(sample)

    df = pd.DataFrame({'sample': sample_lst, 'variant_count': count_lst, 'TMB': tmb_lst})
    df.to_csv('TMB.txt', index=False, sep='\t')
    return sample_lst, count_lst, tmb_lst


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
