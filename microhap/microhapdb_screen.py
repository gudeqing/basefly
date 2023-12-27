import pandas as pd
import gzip
from collections import Counter


"""
https://microhapdb.readthedocs.io/en/latest/citations.html#published-marker-collections-and-allele-frequency-data
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

在千人基因组中:
    Super Population Code	中文注释
    EAS (东亚人群)	CHB：北京汉族；JPT：东京日本人；CHS：南方汉族；CDX：傣族；KHV：越南京族
    EUR (欧洲人群)	CEU：犹他州白人；TSI：意大利托斯卡纳人；FIN：芬兰人；GBR：英国英格兰和苏格兰人；IBS：西班牙伊比利亚人
    AFR (非洲人群)	YRI：尼日利亚约鲁巴人；LWK：肯尼亚卢希亚人；GWD：冈比亚西非人；MSL：塞拉利昂门德人；ESN：尼日利亚埃塞俄比亚人；ASW：美国非洲裔美国人；ACB：巴巴多斯非洲加勒比海人
    AMR (拉丁美洲人群)	MXL：墨西哥裔美国人；PUR：波多黎各人；CLM：哥伦比亚哥伦比亚人；PEL：秘鲁利马人
    SAS (南亚人群)	GIH：印度古吉拉特邦人；PJL：巴基斯坦旁遮普邦人；BEB：孟加拉国人；STU：斯里兰卡泰米尔人；ITU：英国印度泰卢固邦人
    
注意：
1. 这里有不少群体没有对应的frequency信息，可能是信息无法获取
2. 筛选人群频率比较高的位点，比如>5%
3. 避免频繁突变和频繁重组
4. 很明显，无数个等位基因可以最大程度区分任何两个2个体，即两两之间的等位基因都不同
5. 等位基因之间的出现频率比较接近时，在检查随机混合物时，区分度越好，比如2个等位基因，各自人群频率为50%最具备区分度
6. 有多个位点，但等位基因的分布频率相差较大，该如何对位点进行排序
7. 经典人群遗传学给出的额答案是计算有效等位基因数量Ae: effective number of alleles
   如，二等位基因Ae = 1/(p^2+q^2)
8. 如果我们过于强调高平均AE的选择，我们可能会忽视等位基因频率在地理区域之间的巨大差异，从而选择较少的祖先信息基因座。

了解microhapdb
验证数据库中Ae的计算，以'mh01CP-007'为例，基于1KGP人群，计算其Ae
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

"""
# 根据人群名称信息，筛选得到如下候选
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
    'XichangYi',   # 中国
    'ZunyiGelao'   # 中国
]


# marker_pop_ae_value_dict = dict()
# with open('microhapdb/marker-aes.csv') as f:
#     header = f.readline()
#     for line in f:
#         marker, pop, ae = line.strip().split(',')
#         ae = float(ae)
#         # 提取目标人群信息，筛选Ae >= 2的
#         if pop in target_popultation_ids:
#             if ae >= 2:
#                 marker_pop_ae_value_dict[marker+'+'+pop] = ae

"""
@全局检查数据库
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


# 实际上，根据频率表可以计算出每一个marker的AE信息
# 经验算，marker-aes.csv这张表是可以根据这个频率计算表获得，如下：
freq = pd.read_csv('microhapdb/frequency.csv.gz', header=0)
freq = freq.loc[[x in target_popultation_ids for x in freq['Population']]]
ae_df = freq.groupby(['Marker', 'Population', 'Source'])['Frequency'].apply(lambda freqs: 1/sum(x**2 for x in freqs))
ae_df.name = 'Ae'
ae_df = ae_df.sort_index()
ae_df.round(3).to_csv('target_marker_Ae_byFrequency.txt', sep='\t')


# 提取目标marker的信息
target_markers = ae_df.reset_index()['Marker'].unique()
marker_df = pd.read_csv('microhapdb/marker.csv', header=0, index_col=0)
# 有些marker居然不在marker信息表
missed_markers = [x for x in target_markers if x not in marker_df.index]
print(f'有些marker居然不在marker信息表, {missed_markers}')
target_markers = [x for x in target_markers if x in marker_df.index]
target_marker_df = marker_df.loc[list(target_markers)]

# 计算Ae平均值, 即计算千人基因组中的EAS人群，然后加15个其他人群的信息
mean_ae = ae_df.reset_index().groupby('Marker')['Ae'].mean().round(3)
target_marker_df['meanAe'] = mean_ae.loc[target_markers]
target_marker_df = target_marker_df.sort_values(by='meanAe', ascending=False)

# 过滤只包含一个snp的marker
idx = target_marker_df['NumVars'] == 1
print(f'过滤掉{sum(idx)}个只包含一个snp的marker, 他们是 {list(target_marker_df.loc[idx].index)}')
target_marker_df = target_marker_df.loc[~idx]
print(target_marker_df.describe())
target_marker_df.to_csv('target_markers.csv')

max_len = 100
min_len = 50
min_ae = 3
max_NumVars = 10
target_marker_df = target_marker_df.loc[target_marker_df['Extent'] <= max_len]
target_marker_df = target_marker_df.loc[target_marker_df['Extent'] >= min_len]
target_marker_df = target_marker_df.loc[target_marker_df['meanAe'] >= min_ae]
target_marker_df = target_marker_df.loc[target_marker_df['NumVars'] <= max_NumVars]
print(f'我们挑选: ({min_len} < extent < {max_len}) & (meanAe >={min_ae}) & (NumVars < {max_NumVars})的marker')
print(target_marker_df.describe())
target_marker_df.to_csv(f'target_markers.max{max_len}bp.minAe{min_ae}.csv')
print(target_marker_df['NumVars'].value_counts())

# 提取通过筛选的marker的人群Ae信息
target_marker_freq = freq.loc[[x in target_marker_df.index for x in freq['Marker']]]
target_marker_freq.to_csv('target.marker.freq.txt', sep='\t', index=False)

# 生成bed格式的坐标文件，用于注释基因
target_df = target_marker_df.reset_index()[['Chrom', 'Start', 'End', 'Name', 'meanAe', 'Extent', 'NumVars', 'Positions', 'RSIDs', 'Source']]
target_df['Start'] = target_df['Start'] - 1
target_df.to_csv('target.marker.zero-based.txt', sep='\t', index=False)

# 其次想到的数据库是国内的女娲基因组资源, 我们依据这个资源可以对这些marker进行验算和重新排序，得到中国更准确的中国人群信息
# http://bigdata.ibp.ac.cn/NyuWa_variants/search.php
# 通过上述网址可以逐一验证位点信息
"""
FLT3
KIT
RAS
PNPN11
JAK2
CBL
PML
RARA
RUNX1
RUNX1T1
CBFB-MYH11
MLL
CEBPA
NPM1
DNMT3a
TET2
IDH1
IDH2
ASXL1
WT1
"""
with open('allsample.txt') as f:
    lines = f.readlines()
    lines = [x for x in lines if x.strip()]
    for i in range(0, len(lines)+1, 6):
        samplename = lines[i].split()[1]
        accession = lines[i+5].split()[1]
        print(samplename + '\t' + accession)


