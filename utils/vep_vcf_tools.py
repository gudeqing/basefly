import os
from collections import Counter
from pysam import VariantFile
import pandas as pd
import gzip

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


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
