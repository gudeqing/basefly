import json
import os

import pandas as pd
import gzip
from collections import Counter
from itertools import combinations


"""
# microhapdb
https://microhapdb.readthedocs.io/en/latest/citations.html#published-marker-collections-and-allele-frequency-data

# 该数据库包含人群信息如下：
----------------------------------------------------------------------------
ID	Name	Source	翻译
MHDBP-936bc36f79	Asia	vanderGaag2018	亚洲（不对）
ChengduHan	Chengdu Han	Zou2022	成都汉族
CDX	Chinese Dai in Xishuangbanna, China	Byrska-Bishop2022	中国西双版纳傣族
DujiangyanTibetan	Dujiangyan Tibetan	Zou2022	西藏
EAS	East Asia	Byrska-Bishop2022	东亚
HainanHan	Hainan Han	Zou2022	海南汉族
HainanLi	Hainan Li	Zou2022	海南黎族
MHDBP-48c2cfb2aa	Han	Chen2019	汉族
CHB	Han Chinese in Beijing, China	Byrska-Bishop2022	北京汉族
MHDBP-63967b883e	Japanese	Hiroaki2015	日本
SA000010B	Japanese	Kidd2018	日本
JPT	Japanese in Tokyo, Japan	Byrska-Bishop2022	日本东京人
KHV	Kinh in Ho Chi Minh City, Vietnam	Byrska-Bishop2022	越南
SA000936S	Koreans	Kidd2018	韩国
SA004035M	Mongolian	Kidd2018	蒙古
MuliTibetan	Muli Tibetan	Zou2022	西藏
OrdosMongolian	Ordos Mongolian	Zou2022	蒙古
CHS	Southern Han Chinese	Byrska-Bishop2022	南方汉族
mMHseq-TWChinese	TWChinese	Gandotra2020	台湾
WuzhongHui	Wuzhong Hui	Zou2022	中国
XichangYi	Xichang Yi	Zou2022	中国
ZunyiGelao	Zunyi Gelao	Zou2022	中国
----------------------------------------------------------------------------

# 在千人基因组中:
Super Population Code	中文注释
EAS (东亚人群)	CHB：北京汉族；JPT：东京日本人；CHS：南方汉族；CDX：傣族；KHV：越南京族
EUR (欧洲人群)	CEU：犹他州白人；TSI：意大利托斯卡纳人；FIN：芬兰人；GBR：英国英格兰和苏格兰人；IBS：西班牙伊比利亚人
AFR (非洲人群)	YRI：尼日利亚约鲁巴人；LWK：肯尼亚卢希亚人；GWD：冈比亚西非人；MSL：塞拉利昂门德人；ESN：尼日利亚埃塞俄比亚人；ASW：美国非洲裔美国人；ACB：巴巴多斯非洲加勒比海人
AMR (拉丁美洲人群)	MXL：墨西哥裔美国人；PUR：波多黎各人；CLM：哥伦比亚哥伦比亚人；PEL：秘鲁利马人
SAS (南亚人群)	GIH：印度古吉拉特邦人；PJL：巴基斯坦旁遮普邦人；BEB：孟加拉国人；STU：斯里兰卡泰米尔人；ITU：英国印度泰卢固邦人

# 筛选策略思考：
1. 数据库中有不少群体没有对应的frequency信息，猜测是从文献中获得，文献中没有提供相应数据
2. 筛选人群频率比较高的位点，比如>5%
3. 避免频繁突变和频繁重组
4. 很明显，无数个等位基因可以最大程度区分任何两个2个体，即两两之间的等位基因都不同
5. 等位基因之间的出现频率比较接近时，在检查随机混合物时，区分度越好，比如2个等位基因，各自人群频率为50%最具备区分度
6. 有多个位点，但等位基因的分布频率相差较大，该如何对位点进行排序，筛选时可能需要考虑这一点
7. 经典人群遗传学给出的额答案是计算有效等位基因数量Ae: effective number of alleles 如，二等位基因Ae = 1/(p^2+q^2)
8. 如果我们过于强调高平均AE的选择，我们可能会忽视等位基因频率在地理区域之间的巨大差异，从而选择较少的祖先信息基因座。
9. 位于可以高度置信地测序的基因组区域中
10. 建议MRD和chimerism在同一个panel中实现？ctDNA一般有富集的实验过程，可能让嵌合率的结果失真？
11. 筛选多少个marker？https://www.sciencedirect.com/science/article/pii/S0009898122011792
    With guidelines suggesting to monitor at least three informative markers,
    we demonstrate that, for optimized assays, at least 40 biallelic markers need to be screened to achieve enough informative markers in over 99% of cases.
    根据这篇文献，为覆盖99%的情况，建议至少40个双等位基因marker，每次实验中，要求至少存在3个有效marker。
    通常，将markerhap视为一个整体的话，有不少是3等位以上的，因此可能需要的marker数量可能更少。
12. https://pubmed.ncbi.nlm.nih.gov/33582770/这篇文献中，使用了84对距离10bp以内的且连锁不平衡的SNP对，作为造血干细胞移植后检测的marker，检测下限：0.09%， 定量下限（相对背景噪音）0.39%
13. 基于之前microhap分析的实战经验，发现2个基因型如果仅仅相差一个碱基，而且C->T, G->A之间的转换，容易因为测序错误引来偏差或噪音
14. 我们可能需要尽量避免测序错误容易发生的位点，比如前后有3个以上的串联碱基重复的位点
15. 筛选标准，筛选40个以上的marker：
    0. 考虑捕获探针的长度一般为80-120nt
    1. 长度（第一个snp和最后一个snp的距离）：25-150nt
    2. 理论有效等位基因数量Ae: >= 2
    3. 为减少测序错误带来的定量偏差，每个snp前后不能出现超过3个及以上的碱基串联重复
    4. 每个microhap包含[2, 10]个变异位点
    6. 为减少测序错误带来的定量偏差，对于一个microhap，需要谨慎考虑基因型之间相差的碱基仅仅是{C>T} 或 {G>A}的情形
    7. 去除有重叠序列的MHs并选择Ae最大的MHs 


# 进一步了解microhapdb：
## 验证数据库中Ae的计算，以'mh01CP-007'为例，基于1KGP人群，计算其Ae
找人群频率：
    $ grep 'mh01CP-007' <(less frequency.csv.gz) |grep 1KGP
    mh01CP-007,1KGP,A|A|A,0.00115,5196,Byrska-Bishop2022
    mh01CP-007,1KGP,A|C|G,0.14088,5196,Byrska-Bishop2022
    mh01CP-007,1KGP,G|A|A,0.00038,5196,Byrska-Bishop2022
    mh01CP-007,1KGP,G|A|G,0.00192,5196,Byrska-Bishop2022
    mh01CP-007,1KGP,G|C|G,0.85566,5196,Byrska-Bishop2022
计算Ae：
    In [4]: 1/(0.00115**2+0.14088**2+0.00038**2+0.00192**2+0.85566**2)
    Out[4]: 1.3297759816974677
数据库中的Ae记录：
    $ grep 'mh01CP-007' marker-aes.csv|grep 1KGP
    mh01CP-007,1KGP,1.330
验证是否相等：
    1.3297759816974677 ~= 1.330

# 了解NGS的探针引物
在整个NGS实验流程中，会用到4种NGS引物，分别为接头引物、多重PCR引物、封阻引物和杂交捕获探针
1. 接头引物
接头引物的作用是将待测DNA与测序仪建立联系，其本质是一段碱基序列，长度15-80 nt，主要结构可以分为3部分：
a) P5/P7：是与Flow-cell上面寡核苷酸接头相同或互补的片段，能将待测DNA样本固定到测序芯片泳道上；
b) SP1/SP2：测序时测序引物结合部分，是测序启动的位点；
c) Index：用于区分不同样本，每个样本的Index都是唯一的；双端标签接头指的在两端均存在Index，能区分更多样本
d) DNA Insert：待测DNA序列，一般长度为200-300 bp，不会超过1 kb。

2. 封阻引物
捕获测序中阻止捕获探针与文库接头区的杂交，降低捕获探针的脱靶效应，提高捕获测序的特异性，
封阻引物长度一般根据靶DNA序列和接头区域序列相似度进行定制合成。

3. 杂交捕获探针
杂交捕获探针长度一般为80-120nt，它的工作原理是根据DNA碱基互补配对原则，其序列与靶序列互补且有生物素标记。
在捕获探针与NGS文库进行杂交后，探针会与靶序列结合，然后使用链霉素磁珠富集杂交后的DNA文库，
经过洗脱即可获得待测序列的文库。

4. 多重PCR（Multiplex PCR）
在多重PCR中，可以使用多个不同的引物，使多个目标序列在同一反应中被扩增。
多重PCR引物长度一般为15-120 nt，设计时需要注意引物之间的特异性和互补性，以避免引物间相互干扰和非特异性扩增。
多重PCR的优点包括可以在同一反应中扩增多个目标序列，减少反应时间和成本，同时还可以避免样本不足或失效的问题，
目前已经应用到基因分型、遗传育种、生育健康、肿瘤检测于NGS靶向富集等多个领域。

# microhap数据库的Ae信息
In [3]: ae = pd.read_csv('marker-aes.csv')
In [6]: ae['Population'].unique()
Out[6]:  发现包含AE信息的marker涉及人群只有如下这些
array(['1KGP', 'ACB', 'AFR', 'ASW', 'BEB', 'CDX', 'CEU', 'CHB', 'CHS',
       'CLM', 'EAS', 'ESN', 'EUR', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS',
       'ITU', 'JPT', 'KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR',
       'SAS', 'STU', 'TSI', 'YRI'], dtype=object)
@和我们感兴趣的人群取一个交集
In [7]: set(target_popultation_ids ) & set(ae_pops)
Out[7]: {'CDX', 'CHB', 'CHS', 'EAS', 'JPT', 'KHV'}
@也就是说，只有千人基因组计划的人群包含了Ae信息

# 其次想到的数据库是国内的女娲基因组资源, 我们依据这个资源可以对这些marker进行验算和重新排序，得到中国更准确的中国人群信息
# http://bigdata.ibp.ac.cn/NyuWa_variants/search.php

# NGS数据模拟,输入根据microhapdb工具生成的fasta文件和freq文件

"""
# 统计哪些目标群体没有AE信息可得
# hit_pops = {x.split('+')[1] for x in marker_pop_ae_value_dict.keys()}
# print('没有AE信息的群体有：', set(target_popultation_ids) - hit_pops)
# print('包含AE信息的人群有：', hit_pops)

# target_markers = {x.split('+')[0] for x in marker_pop_ae_value_dict.keys()}
# print(f'总共有{len(target_markers)}个marker，其在我们感兴趣的人群中至少有一个Ae值大于2')
#
# # 基于microhapdb的频率数据库，提取目标marker的基因型的人群频率信息
# freq_df = pd.read_csv('microhapdb/frequency.csv.gz', header=0)
# target_freq_df = freq_df.loc[[x in target_markers for x in freq_df['Marker']]]
# target_freq_df = target_freq_df.loc[[x in target_popultation_ids for x in target_freq_df['Population']]]
# target_freq_df.to_csv('target.freq.txt', index=False)


def screen(freq_file='microhapdb/frequency.csv.gz', marker_file='microhapdb/marker.csv',
           genome_file='/home/hxbio04/biosofts/MicroHapulator/microhapulator/data/hg38.fasta'):
    from pysam import FastaFile
    genome = FastaFile(genome_file)

    def genotypes_is_distinguishable(allele_lst):
        # 检查基因型两两之间的差异
        # 假设recipient基因型（A-T-C | A-T-T），而donar基因型（A-T-C | A-T-C)
        # 由于C容易测序错误为T，最后得到的定量误差会增加，因此我们尽量避免选用此类marker
        # 当前，我们选定比较松弛的限制即easy_mix_num>=2才考虑丢弃
        good = True
        easy_mix_num = 0
        for each, other in combinations(allele_lst, 2):
            diff_bases = [(x, y) for x, y in zip(each, other) if x != y]
            if len(diff_bases) == 1:
                if sorted(diff_bases[0]) in [['C', 'T'], ['A', 'G']]:
                    easy_mix_num += 1
        if easy_mix_num >= 2:
            good = False
        return good

    def near_very_short_tandem(row, repeat_len=4):
        # 串联重复区域容易测错，如AAAT有可能测成AAAA, 因此我们也希望能避免此类marker
        NumVars, Extent, Chrom, Start, End, Positions, Positions37, RSIDs, Source, *_ = row
        positions = [int(x) for x in Positions.split(';')]
        discard = False
        for pos in positions:
            up_3bases = genome.fetch(Chrom, pos - repeat_len - 1, pos - 1)
            down_3bases = genome.fetch(Chrom, pos, pos + repeat_len)
            # 如果参考序列有小写的，则认为处于低复杂度区域
            in_low_complex_region = any([x.islower() for x in up_3bases + down_3bases])
            if len(set(up_3bases)) == 1 or len(set(down_3bases)) == 1 or in_low_complex_region:
                discard = True
        return discard

    def get_snp_context(row, extend=10):
        NumVars, Extent, Chrom, Start, End, Positions, Positions37, RSIDs, Source, *_ = row
        positions = [int(x) for x in Positions.split(';')]
        contexts = []
        for pos in positions:
            up_3bases = genome.fetch(Chrom, pos - extend - 1, pos - 1)
            current_base = genome.fetch(Chrom, pos - 1, pos)
            down_3bases = genome.fetch(Chrom, pos, pos + extend)
            contexts.append(up_3bases + f'({current_base})' + down_3bases)
        return ";".join(contexts)

    # 实际上，根据频率表可以计算出每一个marker的AE信息, 经验算，marker-aes.csv这张表是可以根据这个频率计算表获得
    freq = pd.read_csv(freq_file, header=0)

    # 提取marker的在1KGP中的基因型信息
    marker_1KGP_info = freq.loc[freq['Population'] == '1KGP']
    distinguishable_info = marker_1KGP_info.groupby(['Marker', 'Population'])['Allele'].apply(genotypes_is_distinguishable).reset_index()
    marker_genotype_distinguishable = dict(zip(distinguishable_info['Marker'], distinguishable_info['Allele']))

    # 提取目标人群的信息，并计算有效等位基因数量Ae
    target_popultation_ids = [
        # 'JPT',  # 日本东京人，归为EAS
        # 'KHV',  # 越南, 归为EAS
        # 'CHB',  # 北京汉族, 归为EAS
        # 'CHS',  # 南方汉族, 归为EAS
        # 'CDX',  # 西双版纳傣族, 归为EAS
        'EAS',  # 东亚
        # 下面的都是人群相应marker的count信息一般不全，虽然都属东亚范畴，但无法和千人基因组的EAS合并计算，因此单独理出
        # 'MHDBP-936bc36f79',  # 通过阅读原文献，发现其实际是：研究者从欧洲和非洲人群的基因组参考数据中选择了高度变异的微单倍型位点，但确实在亚洲人群上也有测试实验
        'MHDBP-48c2cfb2aa',  # 汉族，至少可以确认在汉族人群上有测试
        'MHDBP-63967b883e',  # 日本
        'ChengduHan',  # 成都汉族？
        'DujiangyanTibetan',  # 西藏？
        'HainanHan',  # 海南汉族
        'HainanLi',  # 海南黎族
        'SA000010B',  # 日本
        'SA000936S',  # 韩国
        'SA004035M',  # 蒙古
        'MuliTibetan',  # 西藏
        'OrdosMongolian',  # 蒙古
        'mMHseq-TWChinese',  # 中国
        'WuzhongHui',  # 中国
        'XichangYi',  # 中国
        'ZunyiGelao'  # 中国
    ]
    freq = freq.loc[[x in target_popultation_ids for x in freq['Population']]]
    ae_df = freq.groupby(['Marker', 'Population', 'Source'])['Frequency'].apply(lambda freqs: 1/sum(x**2 for x in freqs))
    ae_df.name = 'Ae'
    ae_df = ae_df.sort_index()
    ae_df.round(3).to_csv('target_marker_Ae_byFrequency.txt', sep='\t')

    # 提取目标marker的信息
    target_markers = ae_df.reset_index()['Marker'].unique()
    marker_df = pd.read_csv(marker_file, header=0, index_col=0)
    # 有些marker居然不在marker信息表
    missed_markers = [x for x in target_markers if x not in marker_df.index]
    print(f'frequency中有些marker居然不在marker信息表, {missed_markers}')
    target_markers = [x for x in target_markers if x in marker_df.index]
    target_marker_df = marker_df.loc[list(target_markers)]

    # 计算Ae平均值, 即计算千人基因组中的EAS人群，然后加15个文献来源的其他人群信息
    mean_ae_df = ae_df.reset_index().groupby('Marker')['Ae'].mean().round(3)
    max_ae_df = ae_df.reset_index().groupby('Marker')['Ae'].max().round(3)
    min_ae_df = ae_df.reset_index().groupby('Marker')['Ae'].min().round(3)
    target_marker_df['meanAe'] = mean_ae_df.loc[target_marker_df.index]
    target_marker_df['maxAe'] = max_ae_df.loc[target_marker_df.index]
    target_marker_df['minAe'] = min_ae_df.loc[target_marker_df.index]
    target_marker_df = target_marker_df.sort_values(by='maxAe', ascending=False)

    # 过滤只包含一个snp的marker
    idx = target_marker_df['NumVars'] == 1
    print(f'过滤掉{sum(idx)}个只包含一个snp的marker, 他们是 {list(target_marker_df.loc[idx].index)}')
    # target_marker_df = target_marker_df.loc[~idx]
    # target_marker_df.to_csv('target_markers.csv')
    min_len = 25
    max_len = 300
    min_ae = 2
    min_NumVars = 2
    max_NumVars = 7
    filter_criterion = [min_len, max_len, min_NumVars, max_NumVars]

    init_marker_num = target_marker_df.shape[0]
    print('起始marker数量', init_marker_num)
    target_marker_df = target_marker_df.loc[target_marker_df['Extent'] <= max_len]
    target_marker_df = target_marker_df.loc[target_marker_df['Extent'] >= min_len]
    len_filtered_num = init_marker_num - target_marker_df.shape[0]
    print('进一步根据长度过滤掉的marker数量:', len_filtered_num)

    target_marker_df = target_marker_df.loc[target_marker_df['maxAe'] >= min_ae]
    ae_filtered_num = init_marker_num - len_filtered_num - target_marker_df.shape[0]
    print('进一步根据minAe过滤掉的marker数量:', ae_filtered_num)

    target_marker_df = target_marker_df.loc[target_marker_df['NumVars'] >= min_NumVars]
    target_marker_df = target_marker_df.loc[target_marker_df['NumVars'] <= max_NumVars]
    numvars_filtered_num = init_marker_num - len_filtered_num - ae_filtered_num - target_marker_df.shape[0]
    print('进一步根据SNP数量过滤掉的marker数量:', numvars_filtered_num)

    # 根据marker的序列特征过滤
    target_marker_df = target_marker_df.loc[~target_marker_df.apply(near_very_short_tandem, axis=1)]
    complex_filtered_num = init_marker_num - len_filtered_num - ae_filtered_num - numvars_filtered_num - target_marker_df.shape[0]
    print('进一步根据SNP位点前后特征（碱基串联重复）过滤掉的marker数量', complex_filtered_num)

    target_marker_df = target_marker_df.loc[[marker_genotype_distinguishable[x] for x in target_marker_df.index]]
    genotype_filtered_num = init_marker_num - len_filtered_num - ae_filtered_num - numvars_filtered_num - complex_filtered_num - target_marker_df.shape[0]
    print('进一步根据基因型之间是否容易因为测序错误易混淆而过滤掉的marker数量', genotype_filtered_num)

    # 增加每个marker的snp的前后碱基信息
    target_marker_df['SnpContext'] = target_marker_df.apply(get_snp_context, axis=1)

    # 总结
    print(f'Total {target_marker_df.shape[0]} markers left')
    print(target_marker_df.describe())
    target_marker_df.to_csv(f'target.markers.csv')
    print('marker 包含的snp数量分布统计：')
    print(target_marker_df['NumVars'].value_counts())
    print('marker 在染色体上的分布统计')
    print(target_marker_df['Chrom'].value_counts())
    print('marker 文献来源分布统计')
    print(target_marker_df['Source'].value_counts())
    print('相同染色体上marker之间的距离信息统计,仅显示距离最小的前5名, 单位为Mb')
    target_marker_df['center_pos'] = [int(x) for x in (target_marker_df['Start'] + target_marker_df['End'])/2]

    def compute_differences(lst, top=5):
        n = len(lst)
        if n <= 1:
            return []
        # 创建一个列表来存储差值
        differences = set()
        # 遍历列表中的每对元素
        for i in range(n):
            for j in range(i + 1, n):
                # 计算差值并添加到结果列表中
                dist = abs(lst[j] - lst[i])
                differences.add(round(dist*1e-6, 2))
        return sorted(differences)[:top]

    print(target_marker_df.groupby('Chrom')['center_pos'].apply(compute_differences))

    # 提取通过筛选的marker的人群Ae信息
    target_marker_freq = freq.loc[[x in target_marker_df.index for x in freq['Marker']]]
    target_marker_freq.to_csv('target.marker.freq.txt', sep='\t', index=False)

    # 生成bed格式的坐标文件，用于注释基因
    target_df = target_marker_df.reset_index()[['Chrom', 'Start', 'End', 'Name', 'meanAe', 'Extent', 'NumVars', 'Positions', 'RSIDs', 'Source']]
    target_df['Start'] = target_df['Start'] - 1
    target_df.columns = ['#Chrom', 'Start', 'End', 'Name', 'meanAe', 'Extent', 'NumVars', 'Positions', 'RSIDs', 'Source']
    target_df.to_csv('target.marker.zero-based.txt', sep='\t', index=False)


def prepare_simulation_reference(fasta_file, freq_file, out_prefix='SimDiploid', seed=11):
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    import random
    microhap_seq = dict()
    microhap_snp_dict = dict()
    for record in SeqIO.parse(fasta_file, format='fasta'):
        name, snp_pos = record.description.split()[-1].split('=')
        snp_pos = [int(x) for x in snp_pos.split(',')]
        ref_bases = [record.seq[x] for x in snp_pos]
        microhap_snp_dict[name] = [snp_pos, ref_bases]
        microhap_seq[name] = record.seq.__str__()

    microhap_genotypes = dict()
    with open(freq_file) as f:
        _header = f.readline()
        for line in f:
            lst = line.strip().split()
            microhap_genotypes.setdefault(lst[0], []).append(lst[1])

    random.seed(seed)
    chrom = dict()
    chrom2 = dict()
    marker_info = dict()
    for marker, genotypes in microhap_genotypes.items():
        marker_info.setdefault(marker, [])
        select_idx = random.randint(0, len(genotypes)-1)
        # print(genotypes, select_idx)
        genotype = genotypes[select_idx]
        positions, ref_bases = microhap_snp_dict[marker]
        marker_name = f'>{marker} genotype={genotype} snps={",".join((str(x) for x in positions))}'
        marker_seq_lst = list(microhap_seq[marker])
        for p, b in zip(positions, genotype.split('|')):
            marker_seq_lst[p] = b
        chrom[marker_name] = ''.join(marker_seq_lst)
        is_homo = random.randint(0, 1)
        if is_homo:
            chrom2[marker_name] = ''.join(marker_seq_lst)
            marker_info[marker].append([genotype])
        else:
            genotypes = sorted(set(genotypes) - {genotype})
            select_idx2 = random.randint(0, len(genotypes) - 1)
            genotype2 = genotypes[select_idx2]
            marker_info[marker].append([genotype, genotype2])
            marker_name = f'>{marker} genotype={genotype2} snps={",".join((str(x) for x in positions))}'
            for p, b in zip(positions, genotype2.split('|')):
                marker_seq_lst[p] = b
            chrom2[marker_name] = ''.join(marker_seq_lst)
    out_name = f'{out_prefix}.chrom1.fa'
    out_name2 = f'{out_prefix}.chrom2.fa'
    out_name3 = f'{out_prefix}.genotype.json'
    with open(out_name, 'w') as f1, open(out_name2, 'w') as f2, open(out_name3, 'w') as f3:
        for name, seq in chrom.items():
            f1.write(name+'\n')
            f1.write(seq+'\n')
        for name, seq in chrom2.items():
            f2.write(name+'\n')
            f2.write(seq+'\n')
        json.dump(marker_info, f3, indent=2)
    return out_name, out_name2


def simulate_data(chrom1=None, chrom2=None, fasta_file=None, freq_file=None, sample='SampleX', outdir='.', seed=11, depth=600, insert_size=350, insert_size_sd=50, simulator='/home/hxbio04/biosofts/ART/art_bin_VanillaIceCream/art_illumina'):
    if not chrom1:
        chrom1, chrom2 = prepare_simulation_reference(fasta_file, freq_file, out_prefix=os.path.join(outdir, sample), seed=seed)
    # simulator = 'docker run --rm -i vlr37/art_illumina:latest art_illumina -ss HS25'
    cmd1 = f'{simulator} --noALN --paired -l 150 -rs {seed} -m {insert_size} -s {insert_size_sd} -f {depth/2} -i {chrom1} --id {sample}c1 -o {sample}_1.'
    cmd2 = f'{simulator} --noALN --paired -l 150 -rs {seed} -m {insert_size} -s {insert_size_sd} -f {depth/2} -i {chrom2} --id {sample}c2 -o {sample}_2.'
    print(cmd1)
    print(cmd2)
    os.system(cmd1)
    os.system(cmd2)
    fq1 = f'{outdir}/{sample}.R1.fastq'
    fq2 = f'{outdir}/{sample}.R2.fastq'
    os.system(f'cat {sample}*.1.fq > {fq1}')
    os.system(f'cat {sample}*.2.fq > {fq2}')
    os.system(f'rm {sample}*.1.fq {sample}*.2.fq')
    return fq1, fq2


def batch_simulation(fasta_file, freq_file, depths=(300, 500, 700, 1000), insert_sizes=(200, 250, 300, 350), insert_size_sd=50, mix_ratios=(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.995, 0.999)):
    """
    供体和受体样本模拟及混合样本模拟
    :param fasta_file: 参考基因组或marker序列，由microhapdb的工具生成
    :param freq_file: marker的基因频率信息，有microhapdb的工具生成，模拟个人样本时，将从该基因型文件中随机挑选基因型进行模拟
    :param depths: 模拟的总测序深度，如果提供多个值，则按此模拟不同测序深度批次的数据
    :param insert_sizes: 平均插入片段长度, 如果提供多个值，则按此模拟不同批次的数据
    :param insert_size_sd: 插入片段长度的sd值，默认50
    :param mix_ratios: 指定要模拟的混合比例，如果depth*ratio < 1, 则放弃该比例的模拟
    :return:
    """
    donor_name = 'Donor'
    recipient_name = 'Recipient'
    fastq_dirs = []
    for depth in depths:
        for insert_size in insert_sizes:
            outdir = f'fastq_dp{depth}_ins{insert_size}'
            os.makedirs(outdir, exist_ok=True)
            fastq_dirs.append(outdir)
            # 供体和受体样本模拟, seed不一样体现区别
            simulate_data(fasta_file=fasta_file, freq_file=freq_file, sample=f'{donor_name}-{depth}v0.dp{depth}.ins{insert_size}', outdir=outdir, seed=11, depth=depth, insert_size=insert_size, insert_size_sd=insert_size_sd)
            simulate_data(fasta_file=fasta_file, freq_file=freq_file, sample=f'{recipient_name}-0v{depth}.dp{depth}.ins{insert_size}', outdir=outdir, seed=12, depth=depth, insert_size=insert_size, insert_size_sd=insert_size_sd)
            # 混合样本模拟
            for ratio in mix_ratios:
                recipient_depth = int(depth * ratio)
                donor_depth = int(depth * (1 - ratio))
                if recipient_depth >= 1 and donor_depth >= 1:
                    mix_name = 'Mix-DR-{d}v{v}'.format(d=donor_depth, v=recipient_depth)
                    donor_fq1, donor_fq2 = simulate_data(fasta_file=fasta_file, freq_file=freq_file, sample=donor_name+f'-tmp-{mix_name}', outdir=outdir, seed=11, depth=donor_depth, insert_size=insert_size, insert_size_sd=insert_size_sd)
                    recept_fq1, recept_fq2 = simulate_data(fasta_file=fasta_file, freq_file=freq_file, sample=recipient_name+f'-tmp-{mix_name}', outdir=outdir, seed=12, depth=recipient_depth, insert_size=insert_size, insert_size_sd=insert_size_sd)
                    os.system(f'cat {donor_fq1} {recept_fq1} > {outdir}/{mix_name}.dp{depth}.ins{insert_size}.R1.fastq')
                    os.system(f'cat {donor_fq2} {recept_fq2} > {outdir}/{mix_name}.dp{depth}.ins{insert_size}.R2.fastq')
                    os.system(f'rm {donor_fq1} {recept_fq1} {donor_fq2} {recept_fq2}')
                else:
                    print(f"for depth={depth} and mix ratio={ratio}, we cannot simulate by assuming depth < 1: donor_depth={donor_depth}, recipient_depth={recipient_depth}")
    print('simulation success')
    print("begin to analysis")
    result_dirs = []
    for each in fastq_dirs:
        depth = each.split("_")[1].replace("dp", '')
        d_name = f'{donor_name}-{depth}v0'
        r_name = f'{recipient_name}-0v{depth}'
        outdir = f'result_{each}'
        result_dirs.append(outdir)
        os.system(f'python ../../microhap.py -fastq {each} -r1 "(.*?).dp.*R1.fastq" -r2 "(.*?).dp.*R2.fastq" -microhaps ../panel/mypanel-defn.tsv -donor_name {d_name} -recipient_name {r_name} -skip GetSeqErrorMetrics -threads 5 --run --plot --docker -outdir {outdir}')
    # 合并分析结果
    print('merge analysis result')
    analysis_batch(result_dirs, out='merged_prediction_evaluation.xlsx')


def analysis_batch(result_dirs: list, out='merged_prediction_evaluation.xlsx'):
    # 分析同一测序深度，不同片段长度的情况
    # 分析同一片段长度，不同测序深度的情况
    # 每个marker的预测效果：统计每个marker得到的嵌合率和期望值的差异的绝对值的总和
    # exclude_markers = ["mh05KK-120.v1", "mh16HYP-36", "mh11WL-039", "mh05WL-049"]
    marker_effect = []
    mh_count_dfs = []
    mh_depth_dfs = []
    mh_chimerism_dfs = []
    for result_dir in result_dirs:
        result_name = os.path.basename(result_dir)
        depth = result_name.split("_")[2].replace("dp", '')
        insert = result_name.split("_")[3].replace("ins", '')
        a = pd.read_csv(os.path.join(result_dir, 'Report/all.chimerism.csv'), header=0, index_col=0)
        a.columns = [x+f'-dp{depth}-ins{insert}' for x in a.columns]
        mh_chimerism_dfs.append(a)
        # 计算每个marker的在一组模拟数据中的整体效果
        exp_diff = a.sub(a.loc['exp_chimerism'], axis='columns').abs().sum(axis=1)/a.shape[1]
        exp_diff.name = (depth, insert)
        marker_effect.append(exp_diff)
        b = pd.read_csv(os.path.join(result_dir, 'Report/haplotype_count.merged.csv'), header=0, index_col=0)
        b.columns = [x + f'-dp{depth}-ins{insert}' for x in b.columns]
        mh_count_dfs.append(b)
        c = pd.read_csv(os.path.join(result_dir, 'Report/marker_min_depth.merged.csv'), header=0, index_col=0)
        c.columns = [x + f'-dp{depth}-ins{insert}' for x in c.columns]
        mh_depth_dfs.append(c)

    result = pd.DataFrame(marker_effect)
    result.index.name = ('depth', 'insert_size')
    result.to_excel(out)
    chimerisms = pd.concat(mh_chimerism_dfs, axis=1)
    chimerisms.to_excel('all.chimerism.xlsx')
    counts = pd.concat(mh_count_dfs, axis=1)
    depths = pd.concat(mh_depth_dfs, axis=1)
    counts.to_excel('all.haplotype_count.merged.xlsx')
    depths.to_excel('all.marker_min_depth.merged.xlsx')


def y_snps(y_snp_file='YSNPs.txt', genome_file='/home/hxbio04/biosofts/MicroHapulator/microhapulator/data/hg38.fasta'):
    from pysam import FastaFile
    genome = FastaFile(genome_file)
    extend = 10
    with open(y_snp_file) as f, open('YSNPs.anno.txt', 'w') as fw:
        header = f.readline().strip().split('\t')
        header.append('SeqContext')
        fw.write('\t'.join(header)+'\n')
        for line in f:
            lst = line.strip().split('\t')
            YDNA_haplogroup, YSNP, Chrom, pos, *_ = lst
            Chrom = 'chrY'
            pos = int(pos)
            up_3bases = genome.fetch(Chrom, pos - extend - 1, pos - 1)
            current_base = genome.fetch(Chrom, pos - 1, pos)
            down_3bases = genome.fetch(Chrom, pos, pos + extend)
            context = up_3bases + f'({current_base})' + down_3bases
            lst[2] = Chrom
            lst.append(context)
            fw.write('\t'.join(lst)+'\n')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
