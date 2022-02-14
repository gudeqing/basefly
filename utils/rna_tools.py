import os
import csv
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from functools import reduce
from itertools import combinations
import math
from utils.gtfparser import GTF

"""
这里包含的工具是用于rnaseq分析的
"""

def filter_by_class_code(gtf, out, exclude_class_codes=('c', 's', 'p', 'r', '=')):
    """根据class_code过滤出目标新转录本"""
    gtf = GTF(gtf)
    filter = gtf.filter_by_attrs
    # 下面只能过滤出包含class_code的行
    trans_gtf = gtf.table[gtf.table['type'] == 'transcript']
    target_gtf = trans_gtf.loc[~trans_gtf['attrs'].apply(filter, args=('class_code', exclude_class_codes))]
    target_transcripts = {gtf.parse_col9(x)['transcript_id'] for x in target_gtf['attrs']}
    # print('hit', target_transcripts)
    # 根据符合条件的转录本id筛选出目标gtf内容
    target_gtf = gtf.table.loc[gtf.table['attrs'].apply(filter, args=('transcript_id', target_transcripts))]
    # quotechar='' 是为了避免多余的”出现导致gtf格式错误
    target_gtf.to_csv(out, sep='\t', header=False, index=False, quotechar='', quoting=csv.QUOTE_NONE)
    return target_gtf


def filter_by_control_sample(tumor_gtf, normal_gtf, out):
    """
    当前新抗原流程最后未使用该函数
    根据转录本的起始坐标判断两个转录本是否一致（虽然不严格，但绝大部分情况下应该可行）
    :param tumor_gtf:
    :param normal_gtf:
    :param out:
    :return:
    """
    tumor_gtf = GTF(tumor_gtf)
    normal_gtf = GTF(normal_gtf)
    tumor_gtf_table = tumor_gtf.table.set_index(['chr', 'start', 'end'])
    tumor_gtf_table = tumor_gtf_table[tumor_gtf_table['type']=='transcript']
    normal_gtf_table = normal_gtf.table.set_index(['chr', 'start', 'end'])
    normal_gtf_table = normal_gtf_table[normal_gtf_table['type'] == 'transcript']
    tumor_spec = set(tumor_gtf_table.index) - set(normal_gtf_table.index)
    target_transcripts = {GTF.parse_col9(x)['transcript_id'] for x in tumor_gtf_table.loc[tumor_spec]['attrs']}
    filter = tumor_gtf.filter_by_attrs
    target_gtf = tumor_gtf.table.loc[tumor_gtf.table['attrs'].apply(filter, args=('transcript_id', target_transcripts))]
    target_gtf.to_csv(out, sep='\t', header=False, index=False, quotechar='', quoting=csv.QUOTE_NONE)
    return target_gtf


def filter_gtf_for_neoantigen_prediction(tumor_gtf, normal_gtf, out, exclude_class_codes=('c', 's', 'p', 'r', '=')):
    tmp_out = os.path.join(os.path.dirname(out), 'tmp.tumor_specific.gtf')
    filter_by_control_sample(tumor_gtf, normal_gtf, out=tmp_out)
    filter_by_class_code(tmp_out, out, exclude_class_codes=exclude_class_codes)


def merge_to_interval(lst):
    """根据是否连续将区间列表尽可能合并为更大的区间列表"""
    init_lst = sorted(lst)
    s = init_lst[0]
    e = init_lst[0]
    interval = []
    for i in init_lst[1:]:
        if e + 1 == i:
            e += 1
        else:
            interval.append([s, e])
            s, e = i, i
    interval.append([s, e])
    return interval


def get_retained_intron_interval(new_gtf, ref_gtf):
    """
    通过比较新转录本和参考转录本，获得新转录本中的那些落在参考转录本的外显子区域以外的碱基的坐标区间
    因为内含子区域其实是一个相对性的概念，针对具体的转录本而言，注意:
    我这里将结果定义为内含子来源的区间retained_intron，是不准确的：因为其还可能来自基因间区；对于其他转录本而言，该区间可能又是外显子区间
    :param new_gtf: gffcompare注释后的gtf，最好已经经过过滤，从而减少计算量
    :param ref_gtf: 参考基因组的gtf
    :return: {transcript_1: [(0, 10), (17, 29)], transcript_2: Novel, ...}，对于全新的转录本即’cmp_ref'字段缺失
    """
    new_gtf_dict = GTF(new_gtf).to_transcript_dict()
    ref_gtf_dict = GTF(ref_gtf).to_transcript_dict()
    result = dict()
    for transcript, detail in new_gtf_dict.items():
        if 'cmp_ref' in detail:
            exon_pos = detail['exon_pos_lst']
            if detail['strand'] == '-':
                exon_pos = exon_pos[::-1]
            # 创建字典用于坐标转换
            genomic2trans_pos = {y:x for x, y in enumerate(exon_pos)}
            ref_transcript_detail = ref_gtf_dict[detail['cmp_ref']]
            ref_exon_pos = ref_transcript_detail['exon_pos_lst']
            # 找出新转录本外显子区域中没有落在参考转录本外显子区域的坐标
            pos_from_ref_intron = set(exon_pos) - set(ref_exon_pos)
            if pos_from_ref_intron:
                trans_position = [genomic2trans_pos[x] for x in pos_from_ref_intron]
                interval_from_intron = merge_to_interval(trans_position)
                result[transcript] = interval_from_intron
        else:
            print(f'novel transcript {transcript} has no cmp_ref information')
            result[transcript] = None
    # print(result)
    return result


def parse_transdecoder_pep(pep_file):
    """
    把pep file解析成字典，便于提取目标区域的肽段
    :param pep_file:
    :return:
    """
    trans_pep_dict = dict()
    with open(pep_file) as fr:
        for line in fr:
            if line.startswith('>'):
                # 开始整理当前序列表头信息
                header = line.strip()[1:].split()
                record_id = header[0]
                transcript, pos_exp = header[-1].split(':')
                pos_coord, strand = pos_exp.rsplit('(', 1)
                start, end = pos_coord.split('-')
                strand = strand[0]
                record = {
                    'type': header[-4].split(':')[1],
                    'length': int(header[-3].split(':')[1]),
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'transcript': transcript,
                    'pep': ''
                }
                trans_pep_dict[record_id] = record
            else:
                trans_pep_dict[record_id]['pep'] += line.strip()
    return trans_pep_dict


def get_overlap_pep_interval(a, b):
    """
    获取编码区域的交集
    :param a: 待查询区间
    :param b: 编码区间
    :return: 待查询区间中落在编码区间的区间
    """
    # 无交集
    if (a[0] > b[1]) or (a[1] < b[0]):
        return None

    # a 包含于b
    if a[0] >= b[0] and a[1] <= b[1]:
        return a

    # a 包含 b
    if b[0] >= a[0] and b[1] <= a[1]:
        # 因为只有b区间才是编码区，所以要返回b
        return b

    # a 和 b交叉，a在右边
    if (b[0] < a[0]) and (b[1] < a[1]) and (a[0] <= b[1]):
        return [a[0], b[1]]

    # a 和 b交叉，a在左边
    if (b[0] > a[0]) and (b[1] > a[1]) and (a[1] >= b[0]):
        return [b[0], a[1]]


def get_retained_intron_peptide(transdecoder_pep, out_prefix, intron_interval_dict:dict, extend_bases=21):
    """
    根据转录本中可能包含的内含子区间信息，提取对应的peptide信息
    :param transdecoder_pep: transdecoder的pep预测结果文件
    :param intron_interval_dict: 目标区间信息，以字典结果存储，其key必须和transdecoder_pep文件中的转录本名称一致
        such as {'MSTRG.1.1':[[570, 588],[590, 610],[650, 656], [940, 953]]}
    :param extend_bases: 以目标区间为中心向两边延申的碱基数目
    :param out_prefix: 输出文件前缀
    :return:
    """
    result = dict()
    trans_pep_dict = parse_transdecoder_pep(transdecoder_pep)
    novel_transcript_pep_dict = dict()
    for transcript, intervals in intron_interval_dict.items():
        hit_pep_ids = [x for x in trans_pep_dict.keys() if x.startswith(transcript+'.p')]
        if not hit_pep_ids:
            continue
        for pep_id in hit_pep_ids:
            pep_info = trans_pep_dict[pep_id]
            if intervals:
                for interval in intervals:
                    # 先判断interval是否和编码区有交集，并取交集区域
                    # 对交集区域前后进行延长7个氨基酸（21bases)（延长时注意是否会超出边界）
                    # 取出交集区域对应的peptides
                    code_region = [pep_info['start'], pep_info['end']]
                    overlap = get_overlap_pep_interval(interval, code_region)
                    if overlap:
                        if overlap[0] - extend_bases >= pep_info['start']:
                            overlap[0] = overlap[0] - extend_bases
                        else:
                            overlap[0] = pep_info['start']

                        if overlap[1] + extend_bases <= pep_info['end']:
                            overlap[1] = overlap[1] + extend_bases
                        else:
                            overlap[1] = pep_info['end']
                        # 提取overlap对应的peptide ATCATC --> NN 123456 -->  NN （4-1）/3
                        t_start = math.floor((overlap[0] - pep_info['start'])/3)
                        t_end = math.ceil((overlap[1] - pep_info['start'])/3)
                        new_pep_id = f'{pep_id}:Coding_POS={overlap[0]}-{overlap[1]}:IR_POS={interval[0]}-{interval[1]}'
                        result[new_pep_id] = pep_info['pep'][t_start:t_end].strip('*')
            else:
                # 对于全新的转录本(没有cmp_ref注释），则提取整个转录本的pep
                novel_transcript_pep_dict[pep_id] = pep_info['pep']
    if result:
        with open(out_prefix+'.intron_retained.faa', 'w') as f:
            for k, v in result.items():
                f.write(f'>{k}\n{v}\n')
    else:
        print('没有提取出任何疑似源于内含子区间的peptide')

    if novel_transcript_pep_dict:
        with open(out_prefix+'.novel_transcript.faa', 'w') as f:
            for k, v in novel_transcript_pep_dict.items():
                f.write(f'>{k}\n{v}\n')
    else:
        print('没有全新转录本的输出')

    return result, novel_transcript_pep_dict


def segment_peptides(pep_dict:dict, length=(8, 11)):
    """
    按照指定长度进行蛋白质切割，返回字典，key=(左侧5，中间，右侧5）, value=pep_id(from=提取坐标)
    :param pep_dict:
    :param length:
    :return:
    """
    result = dict()
    for k in range(length[0], length[1] + 1):
        for pep_id, sequence in pep_dict.items():
            for i in range(len(sequence) + 1 - k):
                n_flank = sequence[i-5:i] if i > 5 else ''
                c_flank = sequence[i+k:i+k+5]
                result.setdefault((n_flank, sequence[i:i+k], c_flank), set()).add(f'{pep_id}(from={i}-{i+k})')
    return result


def find_potential_intron_peptides(tumor_gtf, normal_gtf, ref_gtf, tumor_transdecoder_pep, normal_transdecoder_pep, out_prefix,
                                  mhc1_pep_len=(8, 11), mhc2_pep_len=(12, 18), ignore_novel_transcript=True, alleles:tuple=('HLA-A', 'HLA-B')):
    """
    :param tumor_gtf: 使用gffcompare注释过的gtf，是过滤后的gtf
    :param normal_gtf: 使用gffcompare注释过的gtf
    :param ref_gtf: 参考基因组的gtf
    :param tumor_transdecoder_pep: 基于tumor_gtf，使用transdecoder预测的pep文件
    :param normal_transdecoder_pep: 基于normal_gtf，使用transdecoder预测的pep文件
    :param out_prefix:
    :param pep_len: 肽段长度范围，闭区间
    :param ignore_novel_transcript: 是否忽略全新转录本
    :param alleles: HLA基因型，用于制作MHCflurry的输入
    :return:
    """
    # 获得肿瘤组织的结果
    retained_intron_interval_dict = get_retained_intron_interval(tumor_gtf, ref_gtf)
    t_intron_pep, t_novel_pep = get_retained_intron_peptide(tumor_transdecoder_pep, out_prefix=out_prefix+'.tumor',
                                                            intron_interval_dict=retained_intron_interval_dict)
    # t_intron_pep = {k:v['pep'] for k, v in parse_transdecoder_pep(tumor_transdecoder_pep).items()}
    # print('xxx', len(t_intron_pep))

    # 获取正常组织的结果
    # retained_intron_interval_dict = get_retained_intron_interval(normal_gtf, ref_gtf)
    # n_intron_pep, n_novel_pep = get_retained_intron_peptide(normal_transdecoder_pep, out_prefix=out_prefix+'.normal',
    #                                                         intron_interval_dict=retained_intron_interval_dict)
    # 对于正常组织，将所有新转录本对应的蛋白质都进行切割
    n_intron_pep = n_novel_pep = {k:v['pep'] for k, v in parse_transdecoder_pep(normal_transdecoder_pep).items()}

    # 获得肿瘤特异性的
    if not ignore_novel_transcript:
        t_intron_pep.update(t_novel_pep)
        n_intron_pep.update(n_novel_pep)

    for mhc_type in ['mhc1', 'mhc2']:
        print(f'>生成的{mhc_type}的候选neoEpitope序列......')
        len_range = mhc1_pep_len if mhc_type == 'mhc1' else mhc2_pep_len
        t_intron_pep_kmer_dict = segment_peptides(t_intron_pep, length=len_range)
        n_intron_pep_kmer_dict = segment_peptides(n_intron_pep, length=len_range)
        # t_uniq_pep_kmer_set = t_intron_pep_kmer_dict.keys() - n_intron_pep_kmer_dict.keys()
        n_intron_pep_kmers = {x[1] for x in n_intron_pep_kmer_dict.keys()}
        t_uniq_pep_kmer_set = {x for x in t_intron_pep_kmer_dict.keys() if x[1] not in n_intron_pep_kmers}
        print(f'过滤掉了{len(t_intron_pep_kmer_dict)-len(t_uniq_pep_kmer_set)}条在对照样本中也能检测到的肽段')
        print(f'最后找到{len(t_uniq_pep_kmer_set)}条肿瘤特异性肽段')
        t_uniq_pep_kmer_dict = {x: t_intron_pep_kmer_dict[x] for x in t_uniq_pep_kmer_set}

        # save
        if t_uniq_pep_kmer_dict:
            if mhc_type == 'mhc1':
                # 制备用于mhcflurry的输入
                with open(out_prefix+f'.{mhc_type}.uniq_intron_retained.pep_segments.csv', 'w') as f:
                    f.write('source_id,n_flank,peptide,c_flank,allele\n')
                    for v, k in t_uniq_pep_kmer_dict.items():
                        for hla_gene in alleles:
                            f.write(f'{"|".join(k)},{v[0]},{v[1]},{v[2]},{hla_gene}\n')

                # 输出切割前的肽段, 暂无用途
                with open(out_prefix+'.uniq_intron_retained.pep.faa', 'w') as f:
                    found = set()
                    for _, id_set in t_uniq_pep_kmer_dict.items():
                        for k in (id_set - found):
                            k = k.split('(')[0] # 提取原始ID
                            f.write(f'>{k}\n{t_intron_pep[k]}\n')
                            found.add(k)
            else:
                # 制备MixMHC2pred的输入
                with open(out_prefix+f'.{mhc_type}.uniq_intron_retained.pep_segments.faa', 'w') as f:
                    for v, k in t_uniq_pep_kmer_dict.items():
                        f.write(f'>{"|".join(k)}\n{v[1]}\n')
        else:
            print('没有提取出任何肿瘤样本特有的疑似源于内含子区间的peptide')


def check_and_convert_alleles_for_MixMHC2Pred(alleles:tuple):
    """
    为什么要考虑把DPA1和DPB1捆绑在一起考虑，或许可以参考https://pubmed.ncbi.nlm.nih.gov/22526601/
    :param alleles: list such as ['HLA-DRB1*01:01', 'HLA-DPA1*01:01']
    :return:
    """
    allele_mapping = {
        'HLA-DRB1*01:01': 'DRB1_01_01',
        'HLA-DRB1*01:02': 'DRB1_01_02',
        'HLA-DRB1*01:03': 'DRB1_01_03',
        'HLA-DRB1*03:01': 'DRB1_03_01',
        'HLA-DRB1*04:01': 'DRB1_04_01',
        'HLA-DRB1*04:04': 'DRB1_04_04',
        'HLA-DRB1*04:05': 'DRB1_04_05',
        'HLA-DRB1*04:08': 'DRB1_04_08',
        'HLA-DRB1*07:01': 'DRB1_07_01',
        'HLA-DRB1*08:01': 'DRB1_08_01',
        'HLA-DRB1*10:01': 'DRB1_10_01',
        'HLA-DRB1*11:01': 'DRB1_11_01',
        'HLA-DRB1*11:04': 'DRB1_11_04',
        'HLA-DRB1*12:01': 'DRB1_12_01',
        'HLA-DRB1*13:01': 'DRB1_13_01',
        'HLA-DRB1*13:03': 'DRB1_13_03',
        'HLA-DRB1*15:01': 'DRB1_15_01',
        'HLA-DRB1*16:01': 'DRB1_16_01',
        'HLA-DRB3*01:01': 'DRB3_01_01',
        'HLA-DRB3*02:02': 'DRB3_02_02',
        'HLA-DRB4*01:01': 'DRB4_01_01',
        'HLA-DRB4*01:03': 'DRB4_01_03',
        'HLA-DRB5*01:01': 'DRB5_01_01',
        'HLA-DRB5*02:02': 'DRB5_02_02',
        'HLA-DPA1*01:03 and HLA-DPA1*02:01 /DPB1*03:01': 'DPA1_01_03__DPA1_02_01__DPB1_03_01',
        'HLA-DPA1*01:03 and HLA-DPA1*02:01 /DPB1*11:01': 'DPA1_01_03__DPA1_02_01__DPB1_11_01',
        'HLA-DPA1*01:03 and HLA-DPA1*02:01 /DPB1*17:01': 'DPA1_01_03__DPA1_02_01__DPB1_17_01',
        'HLA-DPA1*01:03 and HLA-DPA1*02:02 /DPB1*05:01': 'DPA1_01_03__DPA1_02_02__DPB1_05_01',
        'HLA-DPA1*01:03 /HLA-DPB1*105:01 and HLA-DPB1*126:01': 'DPA1_01_03__DPB1_105_01__DPB1_126_01',
        'HLA-DPA1*01:03/DPB1*02:01': 'DPA1_01_03__DPB1_02_01',
        'HLA-DPA1*01:03/DPB1*03:01': 'DPA1_01_03__DPB1_03_01',
        'HLA-DPA1*01:03/DPB1*04:01': 'DPA1_01_03__DPB1_04_01',
        'HLA-DPA1*01:03/DPB1*06:01': 'DPA1_01_03__DPB1_06_01',
        'HLA-DPA1*01:03/DPB1*20:01': 'DPA1_01_03__DPB1_20_01',
        'HLA-DQA1*02:01/DQB1*02:02': 'DQA1_02_01__DQB1_02_02',
        'HLA-DQA1*03:01/DQB1*03:01': 'DQA1_03_01__DQB1_03_01',
        'HLA-DQA1*03:03/DQB1*03:01': 'DQA1_03_03__DQB1_03_01',
        'HLA-DQA1*05:05/DQB1*03:01': 'DQA1_05_05__DQB1_03_01'
    }
    mapping = dict()
    valid_alleles = set()
    for each in allele_mapping.values():
        official_names = ['HLA-' + x.replace('_', '*', 1).replace('_', ':') for x in each.split('__')]
        mapping[tuple(sorted(official_names))] = each
        valid_alleles.update(official_names)
    mapping = dict(sorted([[k, v] for k, v in mapping.items()], key=lambda x:-len(x[0])))

    valid = set(alleles) & valid_alleles
    if not valid:
        raise Exception('No valid alleles for MixMHC2Pred')
    print('These alleles are not available:', set(alleles) - valid_alleles)
    result = []
    for k, v in mapping.items():
        if not (set(k) - set(alleles)):
            result.append(v)
            break
    if not result:
        print('valid combination are', mapping)
        raise Exception('No valid alleles combination is valid for MixMHC2Pred')
    else:
        print('valid inputs will are', result)
    return result


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
