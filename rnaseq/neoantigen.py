"""
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
4. 如果有正常组织，将肿瘤组织的转录本组装结果 减去 正常组织对应的转录本组装结果，即两者gtf坐标完全一致时则排除。(下面的步骤8理论上已经包含该排除效果,所以未实施该步骤）
5. 转录本编码能力预测策略：
    0. 根据gtf使用gffread提取肿瘤组织的新转录本序列
    a. 预测编码与否：https://rnamining.integrativebioinformatics.me/about（该软件不能给出具体的编码预测，仅预测是否编码）
    b. 对于具有编码潜能的transcript，进一步使用transdecoder(https://github.com/TransDecoder/TransDecoder)进行pipetide预测
6. 提取出内含子对应（gffcompare会尽量给组装出来的转录本一个最近似的参考转录本，这里的内含子是相对参考转录本而言，即没有落在参考转录本外显子区域）
    编码的肽段，且前后各延申7个碱基，这样可以得到肿瘤组织的新肽段集合Tumor_New_Peptide(进行一定长度如[8-11]的切割，保留前后5个氨基酸的flank)
7. 如果有正常组织，对正常组织也采用3，5，6【但不仅仅是扣内含子区域对应的肽段，而是所有新转录本的肽段】的处理，得到正常组织对应的新肽段集合Normal_New_Peptide
8. 得到初步新抗原预测的肽段集合: pre_new_peptide = Tumor_New_Peptide - Normal_New_Peptide
9. 最终的新抗原肽段集合：将pre_new_peptide和参考蛋白组进行比对（使用diamond)，过滤掉identity=100的peptide后得到最终的final_new_peptide
10. 使用mhcflurry-predict-scan等进行MHC-I类新抗原预测，筛选出潜在的最终新抗原肽段

"""

import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, TmpVar
from utils.tidy_tools import get_4digits_hla_genetype
from utils.rna_tools import check_and_convert_alleles_for_MixMHC2Pred
from basefly.commands import stringtie, gffcompare, gffread, transdecoder_predict, \
    RNAmining, TransDecoder_LongOrfs, mhcflurry_predict, MixMHC2pred, makeblastdb, blastp
__author__ = 'gdq'


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


def find_potential_intron_peptides():
    cmd = Command()
    cmd.meta.name = 'findNovelPeptides'
    cmd.meta.desc = '提取内含子保留序列对应的潜在编码肽段，并根据对照样本进行过滤，保留肿瘤样本特有的'
    cmd.runtime.image = ''
    cmd.runtime.tool = f'python {script_path}/utils/rna_tools.py find_potential_intron_peptides'
    cmd.runtime.cpu = 2
    cmd.runtime.memory = 2 * 1024 ** 3
    cmd.args['tumor_gtf'] = Argument(prefix='-tumor_gtf ', type='infile', desc='gtf file')
    cmd.args['normal_gtf'] = Argument(prefix='-normal_gtf ', type='infile', desc='out gtf file')
    cmd.args['ref_gtf'] = Argument(prefix='-ref_gtf ', type='infile', desc='reference gtf')
    cmd.args['tumor_transdecoder_pep'] = Argument(prefix='-tumor_transdecoder_pep ', type='infile', desc='coding prediction result of transdecoder for tumor sample')
    cmd.args['normal_transdecoder_pep'] = Argument(prefix='-normal_transdecoder_pep ', type='infile', desc='coding prediction result of transdecoder for normal sample')
    cmd.args['out_prefix'] = Argument(prefix='-out_prefix ', desc='prefix of output file')
    cmd.args['mhc1_pep_len'] = Argument(prefix='-mhc1_pep_len ', type='int', array=True, default=(8, 11), desc='肽段切割长度范围')
    cmd.args['ignore_novel_transcript'] = Argument(prefix='-ignore_novel_transcript', type='bool', default=False, desc='ignore peptides from brand new transcripts')
    cmd.args['alleles'] = Argument(prefix='-alleles ', array=True, desc='HLA gene list, will be used to format mhcflurry input csv file')
    cmd.outputs['mhcflurry_csv'] = Output(value='{out_prefix}.mhc1.uniq_intron_retained.pep_segments.csv')
    cmd.outputs['MixMHC2pred_faa'] = Output(value='{out_prefix}.mhc2.uniq_intron_retained.pep_segments.faa')
    cmd.outputs['mhcflurry_faa'] = Output(value='{out_prefix}.mhc2.uniq_intron_retained.pep_segments.faa')
    return cmd


def parse_bam_list(bam_lst):
    result = dict()
    with open(bam_lst) as fr:
        for line in fr:
            if line.strip():
                tumor, normal, tumor_bam, normal_bam = line.strip().split()
                result[(tumor, normal)] = (tumor_bam, normal_bam)
    return result


def pipeline():
    wf = Workflow()
    wf.meta.name = 'NeoAntigenPrediction'
    wf.meta.desc = 'This is a pipeline for neo-antigen prediction based on RNA-seq analysis'
    wf.init_argparser()
    # add workflow args
    wf.add_argument('-bams', required=True, help='A file with four columns: [tumor_sample_name, normal_sample_name, tumor_bam_file_path, normal_bam_file_path]')
    wf.add_argument('-genome_fa', required=True, help='reference nucleotide fasta file')
    wf.add_argument('-gtf', required=True, help='genome annotation file in gtf format, gff format is not supported.')
    wf.add_argument('-mhcflurry_models', required=True, help='Directory containing models of mhcflurry')
    wf.add_argument('-hla_genotype', required=True, help='HLA genetype table with header consiting of HLA gene name, index is tumor sample name. Each element data is genetype detail.')
    wf.add_argument('-alleles', default=('DRA', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DPB2'),
                    nargs='+', help='HLA alleles gene list, such as A,B,C,DRA(B),DQA(B),DPA(B)')
    wf.add_argument('-genome_pep', required=False, help='reference proteome fasta file')
    wf.parse_args()

    # 建参考基因组的蛋白库
    if wf.args.genome_pep:
        makedb_task, args = wf.add_task(makeblastdb(), tag='refProteome')
        args['input_file'].value = wf.args.genome_pep
        args['dbtype'].value = 'prot'
        args['out'].value = 'refProteome'

    for sample_names, bams in parse_bam_list(wf.args.bams).items():
        # 组装获得gtf，并使用参考gtf进行注释
        tumor_sample, normal_sample = sample_names
        for ind, sample in enumerate(sample_names):
            assemble_task, args = wf.add_task(stringtie(), tag=sample)
            args['bam'].value = [bams[ind]]
            args['gene_model'].value = wf.args.gtf
            args['out_gtf'].value = sample + '.assembled.gtf'
            args['conservative'].value = True

            # 注释gtf
            annot_task, args = wf.add_task(gffcompare(), tag=sample, depends=[assemble_task])
            args['gtfs'].value = [assemble_task.outputs['out_gtf']]
            args['strict_match'].value = True
            args['ref'].value = wf.args.gtf
            args['genome'].value = wf.args.genome_fa

            # 根据class_code过滤
            filter_task, args = wf.add_task(filter_gtf_by_class_code(), tag=sample, depends=[annot_task])
            args['gtf'].value = annot_task.outputs['annotated_gtf']
            args['out_gtf'].value = sample + '.filteredByClassCode.gtf'
            args['exclude_class_codes'].value = ('c', 's', 'p', 'r', '=', 'u')
            if ind == 0:
                tumor_filter_task = filter_task
            else:
                normal_filter_task = filter_task

            # 提取转录本序列
            get_seq_task, args = wf.add_task(gffread(), tag=sample, depends=[filter_task])
            args['input_gff'].value = filter_task.outputs['out_gtf']
            args['w'].value = sample + '.target_novel_transcript.fa'
            args['genome'].value = wf.args.genome_fa

            # 使用RNAmining进行编码潜能预测
            coding_predict_task, args = wf.add_task(RNAmining(), tag=sample, depends=[get_seq_task])
            args['query'].value = get_seq_task.outputs['transcript_fa']

            # 转录本编码能力预测 transdecoder
            LongOrfs_task, args = wf.add_task(TransDecoder_LongOrfs(), tag=sample, depends=[coding_predict_task])
            args['t'].value = coding_predict_task.outputs['codings']
            # predict coding region
            decoder_task, args = wf.add_task(transdecoder_predict(), tag=sample, depends=[coding_predict_task, LongOrfs_task])
            args['t'].value = coding_predict_task.outputs['codings']
            args['output_dir'].value = LongOrfs_task.outputs['outdir']
            args['single_best_only'].value = True
            if ind == 0:
                tumor_decoder_task = decoder_task
            else:
                normal_decoder_task = decoder_task

        # 结合肿瘤样本和对照样本的组装结果和蛋白编码预测结果进行特异性的肿瘤新肽段提取
        find_novel_peptide_task, args = wf.add_task(find_potential_intron_peptides(),
                                                    name='findNovelPeptides-'+tumor_sample,
                                                    depends=[tumor_filter_task, normal_filter_task,
                                                             tumor_decoder_task, normal_decoder_task])
        args['tumor_gtf'].value = tumor_filter_task.outputs['out_gtf']
        args['normal_gtf'].value = normal_filter_task.outputs['out_gtf']
        args['ref_gtf'].value = wf.args.gtf
        args['tumor_transdecoder_pep'].value = tumor_decoder_task.outputs['pep_file']
        args['normal_transdecoder_pep'].value = normal_decoder_task.outputs['pep_file']
        args['out_prefix'].value = tumor_sample
        args['alleles'].value = ['HLA-A*02:01','HLA-A*03:01']

        # mhcflurry prediction for MHC-1
        mhcflurry_task, args = wf.add_task(mhcflurry_predict(), tag=tumor_sample, depends=[find_novel_peptide_task])
        args['input_csv'].value = find_novel_peptide_task.outputs['mhcflurry_csv']
        args['out'].value = tumor_sample + '.mhcflurry_prediction.csv'
        args['models'].value = wf.args.mhcflurry_models

        # MixMHC2Pred prediction for MHC-2
        mhc2pred_task, args = wf.add_task(MixMHC2pred(), tag=tumor_sample, depends=[find_novel_peptide_task])
        args['input'].value = find_novel_peptide_task.outputs['MixMHC2pred_faa']
        args['output'].value = tumor_sample + '.MixMHC2pred.txt'
        alleles = get_4digits_hla_genetype(wf.args.hla_genotype, normal_sample, alleles=wf.args.alleles)
        valid_inputs = check_and_convert_alleles_for_MixMHC2Pred(alleles)
        if valid_inputs:
            args['alleles'].value = valid_inputs
        else:
            print(f'skip MixMHC2Pred task for {tumor_sample}')
            wf.tasks.pop(mhc2pred_task.task_id)

        # blast against reference proteome
        if wf.args.genome_pep:
            mhc1_blastp_task, args = wf.add_task(blastp(), tag=tumor_sample+'MHC1', depends=[makedb_task, find_novel_peptide_task])
            args['query'].value = find_novel_peptide_task.outputs['mhcflurry_faa']
            args['task'].value = 'blastp-short'
            args['db'].value = makedb_task.outputs['out']
            args['max_target_seqs'].value = 1
            args['num_threads'].value = 4
            args['ungapped'].value = True
            args['comp_based_stats'].value = '0'
            args['outfmt'].value = 6
            args['out'].value = tumor_sample + '.mhc1.blastp.txt'

            mhc2_blastp_task, args = wf.add_task(blastp(), tag=tumor_sample+'MHC2', depends=[makedb_task, find_novel_peptide_task])
            args['query'].value = find_novel_peptide_task.outputs['MixMHC2pred_faa']
            args['task'].value = 'blastp-short'
            args['db'].value = makedb_task.outputs['out']
            args['max_target_seqs'].value = 1
            args['num_threads'].value = 4
            args['ungapped'].value = True
            args['comp_based_stats'].value = '0'
            args['outfmt'].value = 6
            args['out'].value = tumor_sample + '.mhc2.blastp.txt'

    wf.run()


if __name__ == '__main__':
    pipeline()
