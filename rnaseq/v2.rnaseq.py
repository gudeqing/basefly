import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar
from basefly.commands import fastp, star, salmon, quantiseq, star_fusion, GTF2RefFlat, GTF2rRNAInterval, CollectRnaSeqMetrics, arcas_hla
from basefly.commands import stringtie, gffcompare, gffread, RNAmining, transdecoder_predict, TransDecoder_LongOrfs
from basefly.commands import mhcflurry_predict, makeblastdb, blastp, netMHCIIPan, raw2mgf_with_rawtools, comet
from basefly.commands import mixcr_rnaseq, pMTnet
from utils.tidy_tools import merge_metrics, merge_star_fusion, merge_arcasHLA_genetype, merge_salmon_quant, merge_mhcflurry_result
from utils.get_fastq_info import get_fastq_info
__author__ = 'gdq'

"""
to do list:
3. 引入输出数据库的概念，目前考虑使用mongodb
4. 定义输出结果的标准元数据

本流程功能
输入：
    bulk RNA-seq数据, 可以接收正常对照数据
过程：
    0. 使用fastp进行街头信息去除，使用star进行比对
    1. 使用salmon基于比对结果进行基因和转录本的表达定量
    2. 使用arcasHLA软件进行HLA定型,如果有配对的正常样本，将用正常样本的定型结果作为后续处理的输入
    3. 使用STAR-Fusion进行融合基因鉴定
    4. 使用Mixcr进行TCR分析，可获得CDR3序列
    5. 使用quantiseq基于基因定量结果进行免疫细胞浸润分析
    6. 内含子保留相关新抗原预测:
        1. 使用stringtie进行有参组装得到组装结果gtf文件
        2. 使用gffcompare比较组装结果gtf和参考基因组的gtf, 得到每条转录本的class_code信息, 使用--strict-match参数
        3. 根据class_code进行部分过滤
            a. 排除如下情况对应的转录本
                class_code = 'c' 其对应的是完全包含于参考转录本外显子区域的新转录本(因为我们关心的是内含子保留相关的新抗原, 所以考虑排除）
                class_code = 's' 其对应的转录本极有可能是比对错误导致的组装结果
                class_code = 'p' 其对应的可能是聚合酶跑过头了的产物，但和参考转录本没有任何交集
                class_code = 'r' 其对应的是有50%碱基在重复区域的转录本, 这要求gffcompare时提供soft masked reference.fa
                class_code = 'u' 全新转录本，来自基因间区（因为我们关心的是内含子保留相关的新抗原, 所以考虑排除）
            b. 如果class_code = '=', 则排除【相当于排除参考转录本，保留新转录本】，这要求gffcomapre时，使用--strict-match参数
            即match code '=' is only assigned when all exon boundaries match; code '~' is assigned for intron chain match or single-exon
        4. 根据gtf使用gffread提取肿瘤组织的新转录本序列
        5. 对新转录本编码能力进行预测:
            b. 预测编码与否：https://rnamining.integrativebioinformatics.me/about（该软件不能给出具体的编码预测，仅预测是否编码）
            c. 对于具有编码潜能的transcript，进一步使用transdecoder(https://github.com/TransDecoder/TransDecoder)进行pipetide预测
        6. 对于肿瘤样本，提取出内含子对应（gffcompare会尽量给组装出来的转录本一个最近似的参考转录本，这里的内含子是相对参考转录本而言，即没有落在参考转录本外显子区域）
            编码的肽段，且前后各延申7个碱基，这样可以得到肿瘤组织的新肽段集合Tumor_IR_Peptide(进行一定长度如[8-11]的切割，保留前后5个氨基酸的flank)
        7. 如果有正常对照数据，对正常组织也采用3，4，5的处理，最终得到正常组织对应的新转录本的所有肽段集合Normal_New_Peptide
        8. 得到初步新抗原预测的肽段集合: pre_new_peptide = Tumor_IR_New_Peptide - Normal_New_Peptide
        9. 最终的新抗原肽段集合：将pre_new_peptide和参考蛋白组进行比对（使用blastp-short)，过滤掉identity=100的peptide后得到最终的final_new_peptide用于后续新抗原预测
        10. 使用mhcflurry-predict进行MHC-I类新抗原预测
        11. 使用netMHCPanII进行MHC-II类新抗原预测
        13. 对预测结果进行注释，如添加表达量和基因等信息
        14. 使用pMTnet进行TCR和新抗原和HLA的互作预测，该步骤和10，11并行，作为两种新抗原的预测证据
    7. 如果有赛默飞的质谱数据【从MHC-肽段复合物进行免疫亲和捕获着手，然后提取其中的抗原肽并采用反向高效液相色谱-质谱技术（LC-MS/MS）对其进行分析，从而能鉴定各种细胞上的抗原肽）】，则以预测的新抗原肽段作为参考库，直接进行验证。

输出：
    基因和转录本的表达量
    融合基因
    HLA定型
    免疫浸润比例
    TCR clone鉴定（后续可以考虑用vdjtools做克隆多样性分析）
    内含子来源新抗原及质谱验证率

Mongo数据库:
    1. hg19/hg38基因信息表: 记录基因ID，基因名称，编码与否等功能
    2. 样本信息表:
        sample_id, r1, r2, patient_id, tissue_type, disease_type, primary_site, metastasis_site, ....
    3. 基因表达量表count/tpm
    3. 转录本表达量表count/tpm
    4. quantiseq免疫细胞比例表
    5. 融合基因检测表
    6. HLA定型结果表
    7. TCR结果表
    8. 新抗原结果表
    
    
"""


def filter_gtf_by_class_code():
    cmd = Command()
    cmd.meta.name = 'filterGTFByClassCode'
    cmd.meta.desc = 'filter gtf by class code'
    cmd.runtime.image = ''
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py filter_by_class_code'
    cmd.runtime.cpu = 2
    cmd.runtime.memory = 2 * 1024 ** 3
    cmd.args['gtf'] = Argument(prefix='-gtf ', type='infile', desc='gtf file')
    cmd.args['out_gtf'] = Argument(prefix='-out ', desc='out gtf file')
    cmd.args['exclude_class_codes'] = Argument(prefix='-exclude_class_codes ', array=True, default=('c', 's', 'p', 'r', '='))
    cmd.outputs['out_gtf'] = Output(value='{out_gtf}')
    return cmd


def filter_gtf_by_exp():
    cmd = Command()
    cmd.meta.name = 'filterGTFByExp'
    cmd.meta.desc = 'filter gtf by expression'
    cmd.runtime.image = ''
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py filter_gtf_by_exp'
    cmd.runtime.cpu = 2
    cmd.runtime.memory = 2 * 1024 ** 3
    cmd.args['gtf'] = Argument(prefix='-gtf ', type='infile', desc='gtf file')
    cmd.args['out_gtf'] = Argument(prefix='-out ', desc='out gtf file')
    cmd.args['min_tpm'] = Argument(prefix='-min_tpm ', type='float', default=1.0, desc='minimum tpm value')
    cmd.outputs['out_gtf'] = Output(value='{out_gtf}')
    # result of gffcompare
    cmd.outputs['tmap'] = Output(value='*.tmap')
    cmd.outputs['refmap'] = Output(value='*.refmap')
    return cmd


def find_potential_intron_peptides():
    cmd = Command()
    cmd.meta.name = 'findNovelPeptides'
    cmd.meta.desc = '提取内含子保留序列对应的潜在编码肽段，并根据对照样本进行过滤，保留肿瘤样本特有的'
    cmd.runtime.image = ''
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py find_potential_intron_peptides'
    cmd.runtime.cpu = 2
    cmd.runtime.memory = 2 * 1024 ** 3
    cmd.args['tumor_gtf'] = Argument(prefix='-tumor_gtf ', type='infile', desc='gtf file')
    # cmd.args['normal_gtf'] = Argument(prefix='-normal_gtf ', type='infile', desc='out gtf file')
    cmd.args['ref_gtf'] = Argument(prefix='-ref_gtf ', type='infile', desc='reference gtf')
    cmd.args['tumor_transdecoder_pep'] = Argument(prefix='-tumor_transdecoder_pep ', type='infile', desc='coding prediction result of transdecoder for tumor sample')
    cmd.args['normal_transdecoder_pep'] = Argument(prefix='-normal_transdecoder_pep ', type='infile', level='optional', desc='coding prediction result of transdecoder for normal sample')
    cmd.args['out_prefix'] = Argument(prefix='-out_prefix ', desc='prefix of output file')
    cmd.args['mhc1_pep_len'] = Argument(prefix='-mhc1_pep_len ', type='int', array=True, default=(8, 11), desc='肽段切割长度范围')
    cmd.args['ignore_novel_transcript'] = Argument(prefix='-ignore_novel_transcript', type='bool', default=False, desc='ignore peptides from brand new transcripts')
    cmd.args['alleles'] = Argument(prefix='-alleles ',  array=True, level='optional', desc='HLA gene list, will be used to format mhcflurry input csv file. 如果提供了alleles—file，该参数可以省略')
    cmd.args['alleles_file'] = Argument(prefix='-alleles_file ', type='infile', level='optional', desc='记录HLA定型结果的文件，比如arcasHLA的输出tsv文件或者merge_hisat_genotype的结果文件')
    cmd.args['sample_id'] = Argument(prefix='-sample_id ', level='optional', desc='sample_id用于从alleles_file中提取目标样本的基因信息，默认用alleles_file中的第一个样本')
    cmd.outputs['mhc1_csv'] = Output(value='{out_prefix}.mhc1.uniq_intron_retained.pep_segments.csv')
    cmd.outputs['mhc1_faa'] = Output(value='{out_prefix}.mhc1.uniq_intron_retained.pep_segments.faa')
    cmd.outputs['mhc2_faa'] = Output(value='{out_prefix}.mhc2.uniq_intron_retained.pep_segments.faa')
    cmd.outputs['mhc2_txt'] = Output(value='{out_prefix}.mhc2.uniq_intron_retained.pep_segments.txt')
    return cmd


def filter_pep_by_blast_id():
    cmd = Command()
    cmd.meta.name = 'filterPepByBlastProteome'
    cmd.meta.desc = '和参考蛋白组比对后，过滤掉能够100%比对上的petides，剩下的将用于新抗原预测'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py filter_pep_by_blast_id'
    cmd.runtime.cpu = 2
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.args['blast_result'] = Argument(prefix='-blast_result ', type='infile', desc='blastp result file')
    cmd.args['raw_files'] = Argument(prefix='-raw_files ', type='infile', array=True, desc='blasted fasta file/csv/txt file')
    cmd.args['out_prefix'] = Argument(prefix='-out_prefix ', desc='output file prefix')
    cmd.outputs['out_faa'] = Output(value='{out_prefix}.faa')
    cmd.outputs['out_csv'] = Output(value='{out_prefix}.csv')
    cmd.outputs['out_txt'] = Output(value='{out_prefix}.txt')
    return cmd


def annotate_mhcflurry_result():
    cmd = Command()
    cmd.meta.name = 'annotateMhcflurryResult'
    cmd.meta.desc = '根据gffcompare的结果文件tmap注释MHCfluryy的结果，主要补充表达量和基因信息'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py annotate_mhcflurry_result'
    cmd.runtime.cpu = 1
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.args['csv_file'] = Argument(prefix='-csv_file ', type='infile', desc='result file of mhcflurry')
    cmd.args['tmap'] = Argument(prefix='-tmap ', type='infile', desc='tmap file generatedby gffcompare')
    cmd.args['gtf'] = Argument(prefix='-gtf ', type='infile', desc='gtf file generated by gffcompare')
    cmd.args['comet_results'] = Argument(prefix='-comet_results ', level='optional', type='infile', desc='txt files generated by comet')
    cmd.args['out'] = Argument(prefix='-out ', desc='output file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def check_and_convert_alleles_for_netMHCIIpan4():
    cmd = Command()
    cmd.meta.name = 'convertAllelesForNetMHCIIpan4'
    cmd.meta.desc = '根据gffcompare的结果文件tmap注释netMHCPan的结果，主要补充表达量和基因信息'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py check_and_convert_alleles_for_netMHCIIpan4'
    cmd.runtime.cpu = 1
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.args['alleles_file'] = Argument(prefix='-alleles_file ', type='infile', level='optional', desc='记录HLA定型结果的文件，比如arcasHLA的输出tsv文件或者merge_hisat_genotype的结果文件')
    cmd.args['sample_id'] = Argument(prefix='-sample_id ', level='optional', desc='sample_id用于从alleles_file中提取目标样本的基因信息，默认用alleles_file中的第一个样本')
    cmd.args['alleles'] = Argument(prefix='-alleles ',  array=True, level='optional', desc='HLA gene list. 如果提供了alleles—file，该参数可以省略')
    cmd.args['support_list'] = Argument(prefix='-support_list ', type='infile', level='optional', desc='包含netMHCIIpan4所支持的HLA基因的文件')
    cmd.outputs['out'] = Output(value='valid_netMHCIIpan4_input_alleles.list')
    return cmd


def annotate_netMHCPan_result():
    cmd = Command()
    cmd.meta.name = 'annotateNetMHCPanResult'
    cmd.meta.desc = '根据gffcompare的结果文件tmap注释netMHCPan的结果，主要补充表达量和基因信息'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py annotate_netMHCpan_result'
    cmd.runtime.cpu = 1
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.args['net_file'] = Argument(prefix='-net_file ', type='infile', desc='result file of mhcflurry')
    cmd.args['pep2id_file'] = Argument(prefix='-pep2id_file ', type='infile', desc='记录peptide和相应蛋白id的文件，由find_potential_intron_peptides产生')
    cmd.args['tmap'] = Argument(prefix='-tmap ', type='infile', desc='tmap file generatedby gffcompare')
    cmd.args['gtf'] = Argument(prefix='-gtf ', type='infile', desc='gtf file generated by gffcompare')
    cmd.args['comet_results'] = Argument(prefix='-comet_results ', level='optional', type='infile', desc='txt files generated by comet')
    cmd.args['out'] = Argument(prefix='-out ', desc='output file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def prepare_pMTnet_input():
    cmd = Command()
    cmd.meta.name = 'PreparepMTnetInput'
    cmd.meta.desc = '根据gffcompare的结果文件tmap注释netMHCPan的结果，主要补充表达量和基因信息'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py prepare_pMTnet_input'
    cmd.args['mixcr_trb'] = Argument(prefix='-mixcr_trb ', type='infile', desc='input TRB file from mixcr')
    cmd.args['antigen_seq'] = Argument(prefix='-antigen_seq ', type='infile', desc='tab separated input peptide file with one column named "peptide"')
    cmd.args['hla_file'] = Argument(prefix='-hla_file ', type='infile', desc='hla genotyped file generated by arcasHLA')
    cmd.args['sample_id'] = Argument(prefix='-sample_id ', level='optional', desc='sample_id用于从alleles_file中提取目标样本的基因信息，默认用alleles_file中的第一个样本')
    cmd.args['out'] = Argument(prefix='-out ', default='pMTnet.input.txt', desc='output file name')
    cmd.args['alleles'] = Argument(prefix='-alleles ', array=True, default=('A', 'B', 'C', 'DRA', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DPB2'), desc='HLA gene to be used')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def prepare_quantiseq_input():
    cmd = Command()
    cmd.meta.name = 'PreparepQuantiseqInput'
    cmd.meta.desc = '根据gffcompare的结果文件tmap注释netMHCPan的结果，主要补充表达量和基因信息'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:2.0'
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py prepare_quantiseq_input'
    cmd.args['expr'] = Argument(prefix='-expr ', type='infile', desc='input salmon quan.gene.sf')
    cmd.args['id2name'] = Argument(prefix='-id2name ', type='infile', default='/opt/quantiseq/deconvolution/TIL10_signature_symbol_converstion.txt', desc='file with two column: gene_id gene_name_accepted_by_quantiseq')
    cmd.args['sample_name'] = Argument(prefix='-sample_name ', desc='sample name')
    cmd.outputs['out'] = Output(value='{sample_name}.quantiseq.input_expr.txt')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'RnaSeqPipeline'
    wf.meta.desc = 'This is a complex pipeline for rnaseq analysis'
    wf.init_argparser()
    # add workflow args
    wf.add_argument('-fastq_info', nargs='+', required=True, help='A list with elements from [fastq file path, fastq parent directory, fastq_info.txt, fastq_info.json]')
    wf.add_argument('-r1_name', default='(.*).R1.fastq', help="python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'")
    wf.add_argument('-r2_name', default='(.*).R2.fastq', help="python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'")
    wf.add_argument('-exclude_samples', default=tuple(), nargs='+', help='samples to exclude from analysis')
    wf.add_argument('-pair_info', required=False, help='tumor normal pair info, two-column txt file, first column is tumor sample name.')
    wf.add_argument('-strandness', default='none', help='strandness of the rna-seq library, value could be one of ["none", "rf", "fr"]')
    wf.add_argument('-ctat_genome_lib_build_dir', required=True, help='star-fusion database dir: ctat_genome_lib_build_dir from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz; In addition, you need to CreateSequenceDictionary for ref_genome.fa')
    wf.add_argument('-hla_db', required=False, help='arcasHLA database, please refer to https://github.com/RabadanLab/arcasHLA')
    wf.add_argument('-hla_file', required=False, help='hla genotype result file: sample_id as row name, gene name as header. At least one input of hla_db and hla_file must be provided!')
    wf.add_argument('-mhc1_genes', default=('A', 'B', 'C'), nargs='+', help='HLA alleles gene list, such as A,B,C')
    wf.add_argument('-mhc2_genes', default=('DRA', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DPB2'), nargs='+', help='HLA alleles gene list, such as DRA(B),DQA(B),DPA(B), 这里并不意味着预测软件能处理这里能提供的所有HLA基因，因此后续结果统计时，应该以实际能用的基因为主')
    wf.add_argument('-genome_pep', required=False, help='reference proteome fasta file, you may download this file from https://www.gencodegenes.org/. default will use pep file in ctat_genome_lib_build_dir, but you are recommended to use the latest pep from GENCODE to reduce false prediction of neoantigen')
    wf.add_argument('-ms_data', required=False, help='raw thermo MS data information, with two columns, first column is sample id, and the second column is raw data directory')
    wf.add_argument('-comet_params', required=False, help='parameter file for comet, please refer to https://comet-ms.sourceforge.net/parameters/parameters_202101/')
    wf.parse_args()

    # add top_vars
    wf.add_topvars(dict(
        fusionIndex=TopVar(value=wf.args.ctat_genome_lib_build_dir, type='indir'),
    ))
    if wf.args.hla_db:
        wf.topvars['hla_database'] = TopVar(value=wf.args.hla_db, type='indir')
    if wf.args.hla_file:
        wf.topvars['hla_file'] = TopVar(value=wf.args.hla_file, type='infile')
    if (wf.args.hla_db is None) and (wf.args.hla_file is None):
        raise Exception('At least one input of hla_db and hla_file must be provided!')
    if wf.args.ms_data:
        wf.topvars['ms_data'] = TopVar(value=wf.args.ms_data, type='infile')
        wf.topvars['comet_params'] = TopVar(value=wf.args.comet_params, type='infile')
    if wf.args.genome_pep:
        wf.topvars['genome_pep'] = TopVar(value=wf.args.genome_pep, type='infile')

    # 提取fastq信息
    fastq_info = get_fastq_info(fastq_info=wf.args.fastq_info, r1_name=wf.args.r1_name, r2_name=wf.args.r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found or matched!')

    # 提取配对信息
    pair_list = []
    sample_list = []
    if wf.args.pair_info:
        with open(wf.args.pair_info) as f:
            for line in f:
                if line.strip():
                    pairs = line.strip('\n').split('\t')[:2]
                    pair_list.append(pairs)
                    sample_list.extend(pairs)
    if not pair_list:
        # 如果没有配对信息，则所有样本均当作肿瘤样本进行分析
        pair_list = list(zip(fastq_info.keys(), [None]*len(fastq_info)))
        sample_list = list(fastq_info.keys())
    pair_dict = dict(pair_list)

    # 构建连特异性表示词汇字典
    picard_strandness = {'none': 'NONE', 'rf': 'SECOND_READ_TRANSCRIPTION_STRAND', 'fr': 'FIRST_READ_TRANSCRIPTION_STRAND'}

    # 检查ctat_genome_lib_build_dir，方便后面的输入
    expected_files = ['ref_genome.fa.star.idx', 'ref_genome.fa', 'ref_genome.dict', 'ref_annot.cdna.fa', 'ref_annot.gtf']
    your_list = os.listdir(wf.args.ctat_genome_lib_build_dir)
    missed_file_dir = set(expected_files) - set(your_list)
    if missed_file_dir:
        raise Exception(f'Expected file or directory not found in {wf.args.ctat_genome_lib_build_dir}:', missed_file_dir)
    star_index_dir = os.path.join(wf.args.ctat_genome_lib_build_dir, 'ref_genome.fa.star.idx')
    ref_genome_fa = os.path.join(wf.args.ctat_genome_lib_build_dir, 'ref_genome.fa')
    ref_genome_dict = os.path.join(wf.args.ctat_genome_lib_build_dir, 'ref_genome.dict')
    transcripts = os.path.join(wf.args.ctat_genome_lib_build_dir, 'ref_annot.cdna.fa')
    ref_genome_gtf = os.path.join(wf.args.ctat_genome_lib_build_dir, 'ref_annot.gtf')
    ref_genome_pep = os.path.join(wf.args.ctat_genome_lib_build_dir, 'ref_annot.pep')

    # GTF 转化 为refFlat
    gtf_ref_flat_task, args = wf.add_task(GTF2RefFlat(), name='GTF2RefFlat')
    args['gtf'].value = os.path.abspath(ref_genome_gtf)

    # GTF中提取rRNA的区间
    gtf_rrna_interval_task, args = wf.add_task(GTF2rRNAInterval(), name="GTF2rRNAInterval", depends=[gtf_ref_flat_task])
    args['gtf'].value = os.path.abspath(ref_genome_gtf)
    args['genome_dict'].value = os.path.abspath(ref_genome_dict)

    # 建参考基因组的蛋白库
    makedb_task, args = wf.add_task(makeblastdb(), tag='refProteome')
    args['input_file'].value = wf.topvars['genome_pep'] if 'genome_pep' in wf.topvars else os.path.abspath(ref_genome_pep)
    args['dbtype'].value = 'prot'
    args['out'].value = 'refProteome'

    # 解析质谱数据路径信息
    ms_data = dict()
    if wf.args.ms_data:
        # 解析质谱数据信息
        with open(wf.args.ms_data) as f:
            for line in f:
                lst = line.strip().split('\t')
                if lst[0] in ms_data:
                    raise Exception(f'sample {lst[0]} is duplicated!!')
                ms_data[lst[0]] = lst[1]

    # 逐一分析单样本
    for sample, reads in fastq_info.items():
        # 不对指定分析的样本进行分析，也不对没有包含在配对信息中的样本进行分析
        if (sample in wf.args.exclude_samples) or (sample not in sample_list):
            continue
        if len(reads) == 2:
            r1s, r2s = reads
        else:
            r1s = reads[0]
            r2s = [None] * len(r1s)
        # print(r1s, r2s)
        # fastq 处理
        fastp_tasks = []
        for ind, (r1, r2) in enumerate(zip(r1s, r2s)):
            fastp_task, args = wf.add_task(fastp(), name=f'fastp-{sample}-{ind}')
            args['read1'].value = r1
            args['out1'].value = f'{sample}.clean.R1.fq.gz'
            if r2 is not None:
                args['read2'].value = r2
                args['out2'].value = f'{sample}.clean.R2.fq.gz'
            args['html'].value = f'{sample}.fastp.html'
            args['json'].value = f'{sample}.fastp.json'
            fastp_tasks.append(fastp_task)

        # alignment
        fastp_task_ids = [x.task_id for x in fastp_tasks]
        star_task, args = wf.add_task(star(sample, sentieon=False), tag=sample, depends=fastp_task_ids)
        args['read1'].value = [x.outputs["out1"] for x in fastp_tasks]
        if len(reads) == 2:
            args['read2'].value = [x.outputs["out2"] for x in fastp_tasks]
        args['genomeDir'].value = os.path.abspath(star_index_dir)

        # salmon quant
        salmon_task, args = wf.add_task(salmon(), tag=sample, depends=[star_task.task_id])
        args['bam'].value = [star_task.outputs['transcript_bam']]
        args['transcripts'].value = os.path.abspath(transcripts)
        args['geneMap'].value = os.path.abspath(ref_genome_gtf)

        # prepare for quantiseq
        pre_quantiseq_task, args = wf.add_task(prepare_quantiseq_input(), tag=sample, depends=[salmon_task])
        args['expr'].value = salmon_task.outputs['gene_exp']
        args['sample_name'].value = sample

        # quantiseq
        quantiseq_task, args = wf.add_task(quantiseq(), tag=sample, depends=[pre_quantiseq_task])
        args['expr'].value = pre_quantiseq_task.outputs['out']
        if sample in pair_dict:
            args['is_tumor'].value = 'TRUE'
        else:
            args['is_tumor'].value = 'FALSE'
        args['prefix'].value = sample

        # fusion identification
        fusion_task, args = wf.add_task(star_fusion(), tag=sample, depends=[star_task.task_id])
        args['genomeLibDir'].value = wf.topvars['fusionIndex']
        args['chimeric_junction'].value = star_task.outputs['chimeric']
        # args['read1'].value = [x.outputs["out1"] for x in fastp_tasks]
        # args['read2'].value = [x.outputs["out2"] for x in fastp_tasks]

        # collectRNAseqMetrics
        metric_task, args = wf.add_task(CollectRnaSeqMetrics(sample), name='collectMetrics-'+sample, depends=[star_task.task_id, gtf_ref_flat_task, gtf_rrna_interval_task])
        args['bam'].value = star_task.outputs['bam']
        args['ref_flat'].value = gtf_ref_flat_task.outputs['out']
        args['ribosomal_interval'].value = gtf_rrna_interval_task.outputs['out_interval_list']
        args['strandness'].value = picard_strandness[wf.args.strandness]

        # HLA-typing，以bam文件输入
        if wf.args.hla_db:
            hla_task, args = wf.add_task(arcas_hla(), tag=sample, depends=[star_task.task_id])
            args['bam'].value = star_task.outputs['bam']
            args['database'].value = wf.topvars['hla_database']
            args['link'].value = True

        # TCR分析,直接以原始数据为输入
        mixcr_rnaseq_task, args = wf.add_task(mixcr_rnaseq(), tag=sample)
        args['out_prefix'].value = sample
        if len(r1s) > 1:
            args['read1'].value = '<(gzcat {})'.format(" ".join(r1s))
            if len(reads) == 2:
                args['read2'].value = '<(gzcat {})'.format(" ".join(r2s))
        else:
            args['read1'].value = r1s[0]
            if len(reads) == 2:
                args['read2'].value = r2s[0]
        args['out_prefix'].value = sample

        # stringtie组装转录本
        assemble_task, args = wf.add_task(stringtie(), tag=sample, depends=[star_task])
        args['bam'].value = [star_task.outputs['bam']]
        args['gene_model'].value = os.path.abspath(ref_genome_gtf)
        args['out_gtf'].value = sample + '.assembled.gtf'
        args['conservative'].value = True

        # filter gtf by exp
        filterbyexp_task, args = wf.add_task(filter_gtf_by_exp(), tag=sample, depends=[assemble_task])
        args['gtf'].value = assemble_task.outputs['out_gtf']
        args['out_gtf'].value = sample + '.filteredByExp.gtf'

        # gffcompare注释gtf
        annot_task, args = wf.add_task(gffcompare(), tag=sample, depends=[filterbyexp_task])
        args['gtfs'].value = [filterbyexp_task.outputs['out_gtf']]
        args['strict_match'].value = True
        args['ref'].value = os.path.abspath(ref_genome_gtf)
        args['genome'].value = os.path.abspath(ref_genome_fa)

        # 根据class_code过滤转录本，筛选出包含非编码区的转录本用于新抗原预测
        filter_task, args = wf.add_task(filter_gtf_by_class_code(), tag=sample, depends=[annot_task])
        args['gtf'].value = annot_task.outputs['annotated_gtf']
        args['out_gtf'].value = sample + '.filteredByClassCode.gtf'
        args['exclude_class_codes'].value = ('c', 's', 'p', 'r', '=', 'u')

        # 提取转录本序列
        get_seq_task, args = wf.add_task(gffread(), tag=sample, depends=[filter_task])
        args['input_gff'].value = filter_task.outputs['out_gtf']
        args['w'].value = sample + '.target_novel_transcript.fa'
        args['genome'].value = os.path.abspath(ref_genome_fa)

        # 使用RNAmining进行编码潜能预测
        coding_predict_task, args = wf.add_task(RNAmining(), tag=sample, depends=[get_seq_task])
        args['query'].value = get_seq_task.outputs['transcript_fa']

        # 转录本编码能力预测 transdecoder
        LongOrfs_task, args = wf.add_task(TransDecoder_LongOrfs(), tag=sample, depends=[coding_predict_task])
        args['t'].value = coding_predict_task.outputs['codings']

        # predict coding region
        decoder_task, args = wf.add_task(transdecoder_predict(), tag=sample, depends=[LongOrfs_task])
        args['t'].value = coding_predict_task.outputs['codings']
        args['output_dir'].value = LongOrfs_task.outputs['outdir']
        args['single_best_only'].value = True

    # ===============-以下进行内含子来源的新抗原分析-===================
    for tumor, normal in pair_list:
        # 结合肿瘤样本和对照样本的组装结果和蛋白编码预测结果进行特异性的肿瘤新肽段提取和筛选
        tumor_filter_task = wf.get_task_by_name('filterGTFByClassCode-'+tumor)
        tumor_decoder_task = wf.get_task_by_name('TransDecoderPredict-'+tumor)
        if normal in fastq_info:
            normal_filter_task = wf.get_task_by_name('filterGTFByClassCode-'+normal)
            normal_decoder_task = wf.get_task_by_name('TransDecoderPredict-'+normal)
            if 'hla_database' in wf.topvars:
                hla_task = wf.get_task_by_name('arcasHLA-'+normal)
            else:
                hla_task = None
        else:
            normal = None
            normal_filter_task = None
            normal_decoder_task = None
            if 'hla_database' in wf.topvars:
                hla_task = wf.get_task_by_name('arcasHLA-'+tumor)
            else:
                hla_task = None
        find_novel_peptide_task, args = wf.add_task(find_potential_intron_peptides(),
                                                    name='findNovelPeptides-' + tumor,
                                                    depends=[tumor_filter_task, tumor_decoder_task, hla_task])
        args['tumor_gtf'].value = tumor_filter_task.outputs['out_gtf']
        args['ref_gtf'].value = os.path.abspath(ref_genome_gtf)
        args['tumor_transdecoder_pep'].value = tumor_decoder_task.outputs['pep_file']
        if normal in fastq_info:
            find_novel_peptide_task.depends.append(normal_decoder_task)
            args['normal_transdecoder_pep'].value = normal_decoder_task.outputs['pep_file']
        args['out_prefix'].value = tumor
        args['alleles_file'].value = hla_task.outputs['hla_genotype_tsv'] if hla_task else wf.topvars['hla_file']
        args['sample_id'].value = normal or tumor

        # 先和最新的参考蛋白组（假设里面不应该包含肿瘤新抗原肽段）比对，再过滤，最后预测
        mhc1_blastp_task, args = wf.add_task(blastp(), tag=tumor + 'MHC1', depends=[makedb_task, find_novel_peptide_task])
        args['query'].value = find_novel_peptide_task.outputs['mhc1_faa']
        args['task'].value = 'blastp-short'
        args['db'].value = makedb_task.outputs['out']
        args['max_target_seqs'].value = 1
        args['num_threads'].value = 4
        args['ungapped'].value = True
        args['comp_based_stats'].value = '0'
        args['outfmt'].value = 6
        args['evalue'].value = 1000
        args['max_hsps'].value = 1
        args['qcov_hsp_perc'].value = 100
        args['out'].value = tumor + '.mhc1.blastp.txt'

        mhc2_blastp_task, args = wf.add_task(blastp(), tag=tumor + 'MHC2', depends=[makedb_task, find_novel_peptide_task])
        args['query'].value = find_novel_peptide_task.outputs['mhc2_faa']
        args['task'].value = 'blastp-short'
        args['db'].value = makedb_task.outputs['out']
        args['max_target_seqs'].value = 1
        args['num_threads'].value = 4
        args['ungapped'].value = True
        args['comp_based_stats'].value = '0'
        args['outfmt'].value = 6
        args['evalue'].value = 1000
        args['max_hsps'].value = 1
        args['qcov_hsp_perc'].value = 100
        args['out'].value = tumor + '.mhc2.blastp.txt'

        # 根据比对结果过滤
        filterMHC1_task, args = wf.add_task(filter_pep_by_blast_id(), tag=tumor + 'MHC1', depends=[mhc1_blastp_task, find_novel_peptide_task])
        args['blast_result'].value = mhc1_blastp_task.outputs['out']
        args['raw_files'].value = [find_novel_peptide_task.outputs['mhc1_csv'], find_novel_peptide_task.outputs['mhc1_faa']]
        args['out_prefix'].value = tumor + '.filtered.mhc1'

        filterMHC2_task, args = wf.add_task(filter_pep_by_blast_id(), tag=tumor + 'MHC2', depends=[mhc2_blastp_task, find_novel_peptide_task])
        args['blast_result'].value = mhc2_blastp_task.outputs['out']
        args['raw_files'].value = [find_novel_peptide_task.outputs['mhc2_txt'], find_novel_peptide_task.outputs['mhc2_faa']]
        args['out_prefix'].value = tumor + '.filtered.mhc2'

        # ------------------run comet------------------
        ms_search_task = None
        ms_search_task2 = None
        convert2mgf_task = None
        if ms_data and (tumor in ms_data):
            convert2mgf_task, args = wf.add_task(raw2mgf_with_rawtools(), tag=tumor)
            args['d'].value = ms_data[tumor]
            # 使用comet软件进行搜库，主要针对MHC-I类型的
            ms_search_task, args = wf.add_task(comet(), tag=tumor + 'MHC1', depends=[convert2mgf_task])
            if wf.args.genome_pep:
                ms_search_task.depends.append(filterMHC1_task)
                args['database'].value = filterMHC1_task.outputs['out_faa']
            else:
                ms_search_task.depends.append(find_novel_peptide_task)
                args['database'].value = find_novel_peptide_task.outputs['mhc1_faa']
            args['param_file'].value = wf.args.comet_params
            args['input_files'].value = convert2mgf_task.outputs['out_files']

            # 使用comet软件进行搜库，主要针对MHC-II类型的
            ms_search_task2, args = wf.add_task(comet(), tag=tumor + 'MHC2', depends=[convert2mgf_task])
            if wf.args.genome_pep:
                ms_search_task.depends.append(filterMHC2_task)
                args['database'].value = filterMHC2_task.outputs['out_faa']
            else:
                ms_search_task.depends.append(find_novel_peptide_task)
                args['database'].value = find_novel_peptide_task.outputs['mhc2_faa']
            args['param_file'].value = wf.topvars['comet_params']
            args['input_files'].value = convert2mgf_task.outputs['out_files']

        # --------------run mhcflurry and netMHCIIpan4--------------
        # mhcflurry prediction for MHC-1
        mhcflurry_task, args = wf.add_task(mhcflurry_predict(), tag=tumor, depends=[filterMHC1_task])
        args['input_csv'].value = filterMHC1_task.outputs['out_csv']
        args['out'].value = tumor + '.mhcflurry_prediction.csv'

        # check_and_convert_alleles_for_netMHCIIpan4
        convert_task, args = wf.add_task(check_and_convert_alleles_for_netMHCIIpan4(), tag=tumor, depends=[hla_task])
        args['alleles_file'].value = hla_task.outputs['hla_genotype_tsv'] if hla_task else wf.topvars['hla_file']
        args['sample_id'].value = normal or tumor

        # netMHCIIpan4 for MHC-2 prediction
        netmhcpanii_task, args = wf.add_task(netMHCIIPan(), tag=tumor, depends=[convert_task, filterMHC2_task])
        args['infile'].value = filterMHC2_task.outputs['out_txt']
        args['alleles_file'].value = convert_task.outputs['out']
        args['inptype'].value = '1'
        args['xls'].value = True
        args['xlsfile'].value = tumor + '.netMHCPanII.txt'
        args['stdout'].value = tumor + '.stdout.txt'
        # --------------end of run mhcflurry and netMHCIIpan4--------------

        # -------------------run pMTnet------------------
        # prepare_pMTnet_input for MHC1
        mixcr_task = wf.get_task_by_name('mixcrRNAseq-'+(tumor or normal))
        mhc1_pMTnet_input_task, args = wf.add_task(prepare_pMTnet_input(), tag=tumor + '-MHC1', depends=[mixcr_task, hla_task, filterMHC1_task])
        args['mixcr_trb'].value = mixcr_task.outputs['TRB']
        args['antigen_seq'].value = filterMHC1_task.outputs['out_faa']
        args['hla_file'].value = hla_task.outputs['hla_genotype_tsv'] if hla_task else wf.topvars['hla_file']
        args['sample_id'].value = normal or tumor

        # prepare_pMTnet_input for MHC2
        mhc2_pMTnet_input_task, args = wf.add_task(prepare_pMTnet_input(), tag=tumor + '-MHC2', depends=[mixcr_task, hla_task, filterMHC2_task])
        args['mixcr_trb'].value = mixcr_task.outputs['TRB']
        args['antigen_seq'].value = filterMHC2_task.outputs['out_faa']
        args['hla_file'].value = hla_task.outputs['hla_genotype_tsv'] if hla_task else wf.topvars['hla_file']
        args['sample_id'].value = normal or tumor

        # run pMTnet for MHC1
        mhc1_pmtnet_task, args = wf.add_task(pMTnet(), tag=tumor + '-MHC1', depends=[mhc1_pMTnet_input_task])
        args['input'].value = mhc1_pmtnet_task.outputs['out']

        # run pMTnet for MHC2
        mhc2_pmtnet_task, args = wf.add_task(pMTnet(), tag=tumor + '-MHC2', depends=[mhc2_pMTnet_input_task])
        args['input'].value = mhc2_pmtnet_task.outputs['out']
        # -----------------end of pMTnet------------------

        # annotate result of mhcflurry
        tumor_assemble_task = wf.get_task_by_name('filterGTFByExp-'+tumor)
        tumor_gffcompare_task = wf.get_task_by_name('gffcompare-'+tumor)
        annotateMhcflurry, args = wf.add_task(annotate_mhcflurry_result(), tag=tumor, depends=[mhcflurry_task, tumor_assemble_task, tumor_gffcompare_task])
        args['csv_file'].value = mhcflurry_task.outputs['out']
        args['tmap'].value = tumor_assemble_task.outputs['tmap']
        args['gtf'].value = tumor_gffcompare_task.outputs['annotated_gtf']
        if convert2mgf_task is not None:
            args['comet_results'].value = ms_search_task.outputs['out']
        args['out'].value = tumor + '.annotated.mhcflurry.csv'

        # annotate result of netMHCpanII
        if netmhcpanii_task is not None:
            annotNetMHCPan, args = wf.add_task(annotate_netMHCPan_result(), tag=tumor,
                                               depends=[netmhcpanii_task, find_novel_peptide_task,
                                                        tumor_assemble_task, tumor_gffcompare_task])
            args['net_file'].value = netmhcpanii_task.outputs['out']
            args['pep2id_file'].value = find_novel_peptide_task.outputs['mhc2_txt']
            args['tmap'].value = tumor_assemble_task.outputs['tmap']
            args['gtf'].value = tumor_gffcompare_task.outputs['annotated_gtf']
            if convert2mgf_task is not None:
                args['comet_results'].value = ms_search_task2.outputs['out']
            args['out'].value = tumor + '.annotated.netMHCPanII.txt'
    # final
    wf.run()

    # merger result
    if wf.success:
        print('Merging Results......')
        merge_metrics(wf.args.outdir, filter_ref='', outdir=os.path.join(wf.args.outdir, 'merged_qc'))
        merge_arcasHLA_genetype(wf.args.outdir, outdir=os.path.join(wf.args.outdir, 'merged_HLA'))
        merge_star_fusion(wf.args.outdir, outdir=os.path.join(wf.args.outdir, 'merged_starfusion'))
        merge_salmon_quant(wf.args.outdir, outdir=os.path.join(wf.args.outdir, 'merged_expression'))
        merge_mhcflurry_result(wf.args.outdir, outdir=os.path.join(wf.args.outdir, 'merged_mhcflurry'))
        # 还需要合并quantiseq、netMHCpanII、pmtnet
        print('...end...')


if __name__ == '__main__':
    pipeline()
