import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, ToWdlTask
from utils.get_fastq_info import get_fastq_info
__author__ = 'gdq'

"""
工作目录 /mnt/nas_101/genarsa/Pipeline_Test/gdq_0523
镜像：singularity shell -B $PWD /mnt/nas_002/gwas/dailsun/sig_image/plink:1.90b6.21--hec16e2b_2.sif
plink的参数特征：
PLINK 1.9 parses each command line as a collection of flags (each of which starts with two dashes1), plus parameters (which immediately follow a flag, and never start with a dash unless that dash is immediately followed by a digit) for those flags. 
1: Actually, that was a lie. With the exceptions of --1 and --23file, PLINK 1.9 allows you to use a single dash in front of each flag. In exchange for saving you some keystrokes, please do yourself a favor and avoid filenames that begin with a dash.

Plink文件格式
map文件格式:
    第一列：染色体编号，未知则为0
    第二列：SNP名称，和bed文件一一对应
    第三列：Position in centimorgans (safe to use dummy value of '0')
    第四列：SNP物理坐标

ped文件格式：
    第一列：Family ID，如果没有，用个体ID替代
    第二列：Individual ID，即个体ID编号
    第三列：Paternal ID，父亲编号
    第四列：Maternal ID，母亲编号
    第五列：Sex（1=male，2=female，0=unknown)
    第六列：Phenotype（0=unknown，1=unaffected，2=affected）
    第七列以后：SNP分型数据，可以是AT CG 或 11 12 或者 A T C G 或 1 1 2 2，每2列代表一个基因型
    
fam文件格式：
    第一列：Family ID，如果没有，用个体ID替代
    第二列：Individual ID，即个体ID编号
    第三列：Paternal ID，父亲编号
    第四列：Maternal ID，母亲编号
    第五列：Sex（1=male，2=female，0=unknown)
    第六列：Phenotype（0/-9 = 未知， 1 = 对照， 2 = 病例）

bim文件格式：
    第一列：染色体编号，0表未知
    第二列：变异标识符
    第三列：Position in centimorgans (safe to use dummy value of '0')
    第四列：碱基坐标
    第五列：ALT ("A1" in PLINK 1.x) allele code
    第六列：REF ("A2" in PLINK 1.x) allele code
    A few notes:
Yes, the ALT column comes before the REF column in a .bim file.
When .bed files are involved, the ALT and REF allele codes will sometimes be swapped, since that's PLINK 1.x's default behavior whenever the true REF allele is less common than the ALT allele in the current dataset. If that's a problem, you can use --ref-allele to swap them back.
It is safe to change a .bim file's extension to .pvar and use --pfile to load it.
Variants with negative bp coordinates are ignored by PLINK.
PLINK 1.9 and 2.0 permit the centimorgan column to be omitted. (However, omission is not recommended if the .bim file needs to be read by other software.)

过滤功能汇总：
功能	As summary statistic	As inclusion criteria
个体基因分型缺失率	--missing	--mind N
SNP基因分型缺失率	--missing	--geno N
等位基因频率	--freq	--maf N
哈迪-温伯格平衡	--hardy	--hwe N
孟德尔误差率	--mendel	--me N M

--geno和--maf命令执行优先级低于--mind，即当同时使用--mind和--geno时，会先剔除缺失率较高的个体，再计算SNP的检出率。

cd /mnt/nas_101/genarsa/gdq_0519/gwas/GWA_tutorial
"""


def plinK_missing():
    cmd = Command()
    cmd.meta.name = 'MissingDetection'
    cmd.meta.desc = "produces sample-based and variant-based missing data reports. " \
                    "If run with --within/--family, the variant-based report is stratified by cluster. 'gz'" \
                    " causes the output files to be gzipped."
    cmd.meta.version = '1.9'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'plink'
    cmd.args['missing'] = Argument(prefix='--missing ', default='', desc="produces sample-based and variant-based missing data reports. If run with --within/--family, the variant-based report is stratified by cluster. 'gz' causes the output files to be gzipped.")
    cmd.args['bed'] = Argument(prefix='--bed ', type='str', desc='bed file name')
    cmd.args['bim'] = Argument(prefix='--bim ', type='str', desc='bim file name')
    cmd.args['fam'] = Argument(prefix='--fam ', type='str', desc='fam file name')
    cmd.args['test_mishap'] = Argument(prefix='--test-mishap', type='bool', default=False, desc='tests whether genotype calls at the two adjacent variants can be used to predict missingness status of the current variant, writing results to plink.missing.hap. This can help one judge the safety of assuming missing calls are randomly distributed. Only autosomal diploid variants with at least 5 missing calls are included, and flanking haplotypes with frequency lower than the --maf threshold are ignored. (Nonfounders are no longer ignored.)')
    cmd.args['out_prefix'] = Argument(prefix='--out ', type='str', default='plink', desc="you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide.")
    # 定义输出
    cmd.outputs['log'] = Output(value="{out_prefix}.log", desc='log file')
    cmd.outputs['imiss'] = Output(value="{out_prefix}.imiss", desc='the proportion of missing SNPs per individual')
    cmd.outputs['lmiss'] = Output(value="{out_prefix}.lmiss", desc='the proportion of missing individuals per SNP')
    cmd.outputs['hh'] = Output(value="{out_prefix}.hh", desc='?')
    return cmd


def plink_filter_snp():
    """
    当同时使用--mind和--geno时，会先剔除缺失率较高的个体，再计算SNP的检出率。
    """
    cmd = Command()
    cmd.meta.name = 'FilterSNP'
    cmd.meta.desc = '针对具体的SNP，使用mind参数过滤掉人群频率低的,针对样本，过滤掉SNP检测缺失率较高的'
    cmd.meta.version = '1.9'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'plink'
    cmd.args['make_bed'] = Argument(prefix='--make-bed', type='fix', desc='creates a new PLINK 1 binary fileset, after applying sample/variant filters and other operations below')
    cmd.args['bed'] = Argument(prefix='--bed ', type='str', desc='bed file name')
    cmd.args['bim'] = Argument(prefix='--bim ', type='str', desc='bim file name')
    cmd.args['fam'] = Argument(prefix='--fam ', type='str', desc='fam file name')
    cmd.args['mind'] = Argument(prefix='--mind ', default='0.02', desc='maximum per-sample, filters out all sample with missing call rates exceeding the provided value (default 0.1) to be removed')
    cmd.args['geno'] = Argument(prefix='--geno ', level='optional', default='0.02', desc='maximum per-variant, filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.')
    cmd.args['oblig-missing '] = Argument(prefix='--oblig-missing', type='infile', default=[], array=True, desc='lets you specify blocks of missing genotype calls for --geno and --mind to ignore. The first file should be a text file with variant IDs in the first column and block names in the second, while the second file should be in .clst format. ')
    cmd.args['out_prefix'] = Argument(prefix='--out ', type='str', default='filterSnp', desc="you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide.")
    # 定义输出
    cmd.outputs['log'] = Output(value="{out_prefix}.log", desc='log file')
    cmd.outputs['bed'] = Output(value="{out_prefix}.bed", desc='bed file')
    cmd.outputs['bim'] = Output(value="{out_prefix}.bim", desc='bim file')
    cmd.outputs['fam'] = Output(value="{out_prefix}.fam", desc='fam file')
    cmd.outputs['hh'] = Output(value="{out_prefix}.hh", desc='?')
    return cmd


def plink_filter_person():
    cmd = Command()
    cmd.meta.name = 'FilterPerson'
    cmd.meta.desc = 'Delete individuals with missingness. 针对样本，过滤掉SNP检测缺失率较高的'
    cmd.meta.version = '1.9'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'plink'
    cmd.args['make_bed'] = Argument(prefix='--make-bed', type='fix', desc='creates a new PLINK 1 binary fileset, after applying sample/variant filters and other operations below')
    cmd.args['bed'] = Argument(prefix='--bed ', type='str', desc='bed file name')
    cmd.args['bim'] = Argument(prefix='--bim ', type='str', desc='bim file name')
    cmd.args['fam'] = Argument(prefix='--fam ', type='str', desc='fam file name')
    cmd.args['mind'] = Argument(prefix='--mind ', default='0.02', desc='maximum per-sample, filters out all sample with missing call rates exceeding the provided value (default 0.1) to be removed')
    cmd.args['oblig-missing '] = Argument(prefix='--oblig-missing', type='infile', default=[], array=True, desc='lets you specify blocks of missing genotype calls for --geno and --mind to ignore. The first file should be a text file with variant IDs in the first column and block names in the second, while the second file should be in .clst format. ')
    cmd.args['out_prefix'] = Argument(prefix='--out ', type='str', default='filterPerson', desc="you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide.")
    # 定义输出
    cmd.outputs['log'] = Output(value="{out_prefix}.log", desc='log file')
    cmd.outputs['bed'] = Output(value="{out_prefix}.bed", desc='bed file')
    cmd.outputs['bim'] = Output(value="{out_prefix}.bim", desc='bim file')
    cmd.outputs['fam'] = Output(value="{out_prefix}.fam", desc='fam file')
    cmd.outputs['hh'] = Output(value="{out_prefix}.hh", desc='?')
    return cmd


def plink_check_sex():
    cmd = Command()
    cmd.meta.name = 'CorrectSex'
    cmd.meta.desc = 'Delete individuals with missingness. 针对样本，过滤掉SNP检测缺失率较高的'
    cmd.meta.version = '1.9'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'plink'
    cmd.args['check_sex'] = Argument(prefix='--check-sex ', type='str', default='', desc='check-sex normally compares sex assignments in the input dataset with those imputed from X chromosome inbreeding coefficients, and writes a report to plink.sexcheck.')
    cmd.args['bed'] = Argument(prefix='--bed ', type='str', desc='bed file name')
    cmd.args['bim'] = Argument(prefix='--bim ', type='str', desc='bim file name')
    cmd.args['fam'] = Argument(prefix='--fam ', type='str', desc='fam file name')
    cmd.args['out_prefix'] = Argument(prefix='--out ', type='str', default='sex', desc="you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide.")
    # 定义输出
    cmd.outputs['log'] = Output(value="{out_prefix}.log", desc='log file')
    cmd.outputs['hh'] = Output(value="{out_prefix}.hh", desc='?')
    cmd.outputs['sexcheck'] = Output(value="{out_prefix}.sexcheck", desc='?')
    return cmd


def plink_extract():
    cmd = Command()
    cmd.meta.name = 'Extract'
    cmd.meta.desc = "extract normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis."
    cmd.meta.version = '1.9'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'plink'
    cmd.args['make_bed'] = Argument(prefix='--make-bed ', type='fix', desc='creates a new PLINK 1 binary fileset, after applying sample/variant filters and other operations below')
    cmd.args['extract'] = Argument(prefix='--extract ', type='infile', desc="normally accepts a text file with a list of variant IDs (usually one per line, but it's okay for them to just be separated by spaces), and removes all unlisted variants from the current analysis.")
    cmd.args['bed'] = Argument(prefix='--bed ', type='str', desc='bed file name')
    cmd.args['bim'] = Argument(prefix='--bim ', type='str', desc='bim file name')
    cmd.args['fam'] = Argument(prefix='--fam ', type='str', desc='fam file name')
    cmd.args['out_prefix'] = Argument(prefix='--out ', type='str', default='extracted', desc="you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide.")
    # 定义输出
    cmd.outputs['log'] = Output(value="{out_prefix}.log", desc='log file')
    cmd.outputs['bed'] = Output(value="{out_prefix}.bed", desc='bed file')
    cmd.outputs['bim'] = Output(value="{out_prefix}.bim", desc='bim file')
    cmd.outputs['fam'] = Output(value="{out_prefix}.fam", desc='fam file')
    cmd.outputs['hh'] = Output(value="{out_prefix}.hh", desc='?')
    return cmd


def plink_remove_low_maf_snp():
    cmd = Command()
    cmd.meta.name = 'FilterLowMaf'
    cmd.meta.desc = 'Remove SNPs with a low MAF frequency'
    cmd.meta.version = '1.9'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'plink'
    cmd.args['make_bed'] = Argument(prefix='--make-bed ', type='fix', desc='creates a new PLINK 1 binary fileset, after applying sample/variant filters and other operations below')
    cmd.args['bed'] = Argument(prefix='--bed ', type='str', desc='bed file name')
    cmd.args['bim'] = Argument(prefix='--bim ', type='str', desc='bim file name')
    cmd.args['fam'] = Argument(prefix='--fam ', type='str', desc='fam file name')
    cmd.args['maf'] = Argument(prefix='--maf ', default=0.05, desc='filters out all variants with minor allele frequency below the provided threshold')
    cmd.args['nonfounders'] = Argument(prefix='--nonfounders', type='bool', default=False, desc="By default, nonfounders are not counted by --freq[x] or --maf/--max-maf/--hwe. Use the --nonfounders flag to include them.")
    cmd.args['out_prefix'] = Argument(prefix='--out ', type='str', default='', desc="you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide.")
    # 定义输出
    cmd.outputs['log'] = Output(value="{out_prefix}.log", desc='log file')
    cmd.outputs['bed'] = Output(value="{out_prefix}.bed", desc='bed file')
    cmd.outputs['bim'] = Output(value="{out_prefix}.bim", desc='bim file')
    cmd.outputs['fam'] = Output(value="{out_prefix}.fam", desc='fam file')
    cmd.outputs['hh'] = Output(value="{out_prefix}.hh", desc='?')
    return cmd


def plink_hardy_weinberg():
    cmd = Command()
    cmd.meta.name = 'FilterHWE'
    cmd.meta.desc = 'Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).'
    cmd.meta.version = '1.9'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'plink'
    cmd.args['make_bed'] = Argument(prefix='--make-bed ', type='fix', desc='creates a new PLINK 1 binary fileset, after applying sample/variant filters and other operations below')
    cmd.args['bed'] = Argument(prefix='--bed ', type='str', desc='bed file name')
    cmd.args['bim'] = Argument(prefix='--bim ', type='str', desc='bim file name')
    cmd.args['fam'] = Argument(prefix='--fam ', type='str', desc='fam file name')
    cmd.args['hwe'] = Argument(prefix='--hwe ', default="1e-6", type='str', desc="默认只对control生效，如果希望对case组生效，需要加上'include-nonctrl'。filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold. We recommend setting a low threshold—serious genotyping errors often yield extreme p-values like 1e-50 which are detected by any reasonable configuration of this test, while genuine SNP-trait associations can be expected to deviate slightly from Hardy-Weinberg equilibrium (so it's dangerous to choose a threshold that filters out too many variants).")
    cmd.args['nonfounders'] = Argument(prefix='--nonfounders', type='bool', default=False, desc="By default, nonfounders are not counted by --freq[x] or --maf/--max-maf/--hwe. Use the --nonfounders flag to include them.")
    cmd.args['out_prefix'] = Argument(prefix='--out ', type='str', default='', desc="you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide.")
    # 定义输出
    cmd.outputs['log'] = Output(value="{out_prefix}.log", desc='log file')
    cmd.outputs['bed'] = Output(value="{out_prefix}.bed", desc='bed file')
    cmd.outputs['bim'] = Output(value="{out_prefix}.bim", desc='bim file')
    cmd.outputs['fam'] = Output(value="{out_prefix}.fam", desc='fam file')
    cmd.outputs['hh'] = Output(value="{out_prefix}.hh", desc='?')
    return cmd


def plink_exclude():
    cmd = Command()
    cmd.meta.name = 'FilterHWE'
    cmd.meta.desc = 'Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).'
    cmd.meta.version = '1.9'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'plink'
    cmd.args['bed'] = Argument(prefix='--bed ', type='str', desc='bed file name')
    cmd.args['bim'] = Argument(prefix='--bim ', type='str', desc='bim file name')
    cmd.args['fam'] = Argument(prefix='--fam ', type='str', desc='fam file name')
    cmd.args['exclude'] = Argument(prefix='--exclude ', type='infile', desc="high inversion regions (inversion.txt [High LD regions])")
    cmd.args['range'] = Argument(prefix='--range', type='bool', desc="With the 'range' modifier, the input file should be in set range format instead.")
    cmd.args['indep_pairwise'] = Argument(prefix='--indep-pairwise ', default='50 5 0.2', desc="stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously")
    cmd.args['out_prefix'] = Argument(prefix='--out ', type='str', default='', desc="you usually want to choose a different output file prefix for each run. --out causes 'plink' to be replaced with the prefix you provide.")
    cmd.outputs['log'] = Output(value="{out_prefix}.log", desc='log file')
    cmd.outputs['prune_in'] = Output(value="{out_prefix}.prune.in", desc='?')
    cmd.outputs['prune_out'] = Output(value="{out_prefix}.prune.out", desc='?')
    return cmd




