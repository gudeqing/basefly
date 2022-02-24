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
5. 转录本编码能力预测策略:
    0. 根据gtf使用gffread提取肿瘤组织的新转录本序列
    a. 预测编码与否：https://rnamining.integrativebioinformatics.me/about（该软件不能给出具体的编码预测，仅预测是否编码）
    b. 对于具有编码潜能的transcript，进一步使用transdecoder(https://github.com/TransDecoder/TransDecoder)进行pipetide预测
6. 对于肿瘤样本，提取出内含子对应（gffcompare会尽量给组装出来的转录本一个最近似的参考转录本，这里的内含子是相对参考转录本而言，即没有落在参考转录本外显子区域）
    编码的肽段，且前后各延申7个碱基，这样可以得到肿瘤组织的新肽段集合Tumor_New_Peptide(进行一定长度如[8-11]的切割，保留前后5个氨基酸的flank)
7. 如果有正常组织，对正常组织也采用3，5，6【但不仅仅是扣内含子区域对应的肽段，而是所有新转录本的肽段】的处理，得到正常组织对应的新肽段集合Normal_New_Peptide
8. 得到初步新抗原预测的肽段集合: pre_new_peptide = Tumor_IR_New_Peptide - Normal_New_Peptide
9. 最终的新抗原肽段集合：将pre_new_peptide和参考蛋白组进行比对（使用blastp-short)，过滤掉identity=100的peptide后得到最终的final_new_peptide用于后续新抗原预测
10. 使用mhcflurry-predict进行MHC-I类新抗原预测
11. 使用MixMHC2Pred和netMHCPanII等进行MHC-II类新抗原预测
13. 对预测结果进行注释，如添加表达量和基因等信息
"""

import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, TmpVar
from utils.tidy_tools import get_4digits_hla_genetype
from utils.rna_tools import check_and_convert_alleles_for_MixMHC2Pred, check_and_convert_alleles_for_netMHCIIpan4
from basefly.commands import stringtie, gffcompare, gffread, transdecoder_predict, \
    RNAmining, TransDecoder_LongOrfs, mhcflurry_predict, MixMHC2pred, makeblastdb, blastp, netMHCIIPan
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
    cmd.args['min_TPM'] = Argument(prefix='-min_TPM ', type='float', default=1.0, desc='Minimum TPM value')
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
    # cmd.args['normal_gtf'] = Argument(prefix='-normal_gtf ', type='infile', desc='out gtf file')
    cmd.args['ref_gtf'] = Argument(prefix='-ref_gtf ', type='infile', desc='reference gtf')
    cmd.args['tumor_transdecoder_pep'] = Argument(prefix='-tumor_transdecoder_pep ', type='infile', desc='coding prediction result of transdecoder for tumor sample')
    cmd.args['normal_transdecoder_pep'] = Argument(prefix='-normal_transdecoder_pep ', type='infile', desc='coding prediction result of transdecoder for normal sample')
    cmd.args['out_prefix'] = Argument(prefix='-out_prefix ', desc='prefix of output file')
    cmd.args['mhc1_pep_len'] = Argument(prefix='-mhc1_pep_len ', type='int', array=True, default=(8, 11), desc='肽段切割长度范围')
    cmd.args['ignore_novel_transcript'] = Argument(prefix='-ignore_novel_transcript', type='bool', default=False, desc='ignore peptides from brand new transcripts')
    cmd.args['alleles'] = Argument(prefix='-alleles ', array=True, desc='HLA gene list, will be used to format mhcflurry input csv file')
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
    cmd.args['out'] = Argument(prefix='-out ', desc='output file')
    cmd.outputs['out'] = Output(value='{out}')
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
    cmd.args['out'] = Argument(prefix='-out ', desc='output file')
    cmd.outputs['out'] = Output(value='{out}')
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
    wf.add_argument('-mhc1_genes', default=('A', 'B', 'C'), nargs='+', help='HLA alleles gene list, such as A,B,C')
    wf.add_argument('-mhc2_genes', default=('DRA', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DPB2'), nargs='+', help='HLA alleles gene list, such as DRA(B),DQA(B),DPA(B), 这里并不意味着预测软件能处理这里能提供的所有HLA基因，因此统计时，应该以实际能用的基因为主')
    wf.add_argument('-genome_pep', required=False, help='reference proteome fasta file')
    wf.parse_args()

    if wf.args.genome_pep:
        # 建参考基因组的蛋白库
        makedb_task, args = wf.add_task(makeblastdb(), tag='refProteome')
        args['input_file'].value = wf.args.genome_pep
        args['dbtype'].value = 'prot'
        args['out'].value = 'refProteome'

    for sample_names, bams in parse_bam_list(wf.args.bams).items():
        # 组装获得gtf，并使用参考gtf进行注释
        tumor_sample, normal_sample = sample_names
        for ind, sample in enumerate(sample_names):
            assemble_task, args = wf.add_task(stringtie(), tag=sample)
            if ind == 0:
                tumor_assemble_task = assemble_task
            args['bam'].value = [bams[ind]]
            args['gene_model'].value = wf.args.gtf
            args['out_gtf'].value = sample + '.assembled.gtf'
            args['conservative'].value = True

            # 注释gtf
            annot_task, args = wf.add_task(gffcompare(), tag=sample, depends=[assemble_task])
            if ind == 0:
                tumor_gffcompare_task = annot_task
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
            decoder_task, args = wf.add_task(transdecoder_predict(), tag=sample, depends=[LongOrfs_task])
            args['t'].value = coding_predict_task.outputs['codings']
            args['output_dir'].value = LongOrfs_task.outputs['outdir']
            args['single_best_only'].value = True
            if ind == 0:
                tumor_decoder_task = decoder_task
            else:
                normal_decoder_task = decoder_task

        mhc1_alleles = get_4digits_hla_genetype(wf.args.hla_genotype, normal_sample, alleles=wf.args.mhc1_genes)
        mhc2_alleles = get_4digits_hla_genetype(wf.args.hla_genotype, normal_sample, alleles=wf.args.mhc2_genes)
        print('HLA Genes:', mhc1_alleles, mhc2_alleles)

        # 结合肿瘤样本和对照样本的组装结果和蛋白编码预测结果进行特异性的肿瘤新肽段提取和筛选
        find_novel_peptide_task, args = wf.add_task(find_potential_intron_peptides(),
                                                    name='findNovelPeptides-'+tumor_sample,
                                                    depends=[tumor_filter_task, normal_filter_task,
                                                             tumor_decoder_task, normal_decoder_task])
        args['tumor_gtf'].value = tumor_filter_task.outputs['out_gtf']
        # args['normal_gtf'].value = normal_filter_task.outputs['out_gtf']
        args['ref_gtf'].value = wf.args.gtf
        args['tumor_transdecoder_pep'].value = tumor_decoder_task.outputs['pep_file']
        args['normal_transdecoder_pep'].value = normal_decoder_task.outputs['pep_file']
        args['out_prefix'].value = tumor_sample
        args['alleles'].value = mhc1_alleles

        netmhcpanii_task = None
        if wf.args.genome_pep:
            # 先和参考蛋白组比对，再过滤，最后预测
            mhc1_blastp_task, args = wf.add_task(blastp(), tag=tumor_sample+'MHC1', depends=[makedb_task, find_novel_peptide_task])
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
            args['out'].value = tumor_sample + '.mhc1.blastp.txt'

            mhc2_blastp_task, args = wf.add_task(blastp(), tag=tumor_sample+'MHC2', depends=[makedb_task, find_novel_peptide_task])
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
            args['out'].value = tumor_sample + '.mhc2.blastp.txt'

            # 根据比对结果过滤
            filterMHC1_task, args = wf.add_task(filter_pep_by_blast_id(), tag=tumor_sample+'MHC1', depends=[mhc1_blastp_task, find_novel_peptide_task])
            args['blast_result'].value = mhc1_blastp_task.outputs['out']
            args['raw_files'].value = [find_novel_peptide_task.outputs['mhc1_csv'], find_novel_peptide_task.outputs['mhc1_faa']]
            args['out_prefix'].value = tumor_sample+'.filtered.mhc1'

            filterMHC2_task, args = wf.add_task(filter_pep_by_blast_id(), tag=tumor_sample+'MHC2', depends=[mhc2_blastp_task, find_novel_peptide_task])
            args['blast_result'].value = mhc2_blastp_task.outputs['out']
            args['raw_files'].value = [find_novel_peptide_task.outputs['mhc2_faa'], find_novel_peptide_task.outputs['mhc2_txt']]
            args['out_prefix'].value = tumor_sample + '.filtered.mhc2'

            # mhcflurry prediction for MHC-1
            mhcflurry_task, args = wf.add_task(mhcflurry_predict(), tag=tumor_sample, depends=[filterMHC1_task])
            args['input_csv'].value = filterMHC1_task.outputs['out_csv']
            args['out'].value = tumor_sample + '.mhcflurry_prediction.csv'
            args['models'].value = wf.args.mhcflurry_models

            valid_inputs = check_and_convert_alleles_for_MixMHC2Pred(mhc2_alleles)
            if valid_inputs:
                # MixMHC2Pred prediction for MHC-2
                mhc2pred_task, args = wf.add_task(MixMHC2pred(), tag=tumor_sample, depends=[filterMHC2_task])
                args['input'].value = filterMHC2_task.outputs['out_faa']
                args['output'].value = tumor_sample + '.MixMHC2pred.txt'
                args['alleles'].value = valid_inputs
            else:
                print(f'skip MixMHC2Pred task for {tumor_sample} since no valid HLA-gene combination available')

            # netMHCIIpan4 prediction
            valid_alleles = check_and_convert_alleles_for_netMHCIIpan4(mhc2_alleles)
            if valid_alleles:
                netmhcpanii_task, args = wf.add_task(netMHCIIPan(), tag=tumor_sample, depends=[filterMHC2_task])
                args['alleles'].value = valid_alleles
                args['infile'].value = filterMHC2_task.outputs['out_txt']
                args['inptype'].value = '1'
                args['xls'].value = True
                args['xlsfile'].value = tumor_sample + '.netMHCPanII.txt'
                args['stdout'].value = tumor_sample + '.stdout.txt'

        else:
            # mhcflurry prediction for MHC-1
            mhcflurry_task, args = wf.add_task(mhcflurry_predict(), tag=tumor_sample, depends=[find_novel_peptide_task])
            args['input_csv'].value = find_novel_peptide_task.outputs['mhcflurry_csv']
            args['out'].value = tumor_sample + '.mhcflurry_prediction.csv'
            args['models'].value = wf.args.mhcflurry_models

            # MixMHC2Pred prediction for MHC-2
            valid_alleles = check_and_convert_alleles_for_MixMHC2Pred(mhc2_alleles)
            if valid_alleles:
                mhc2pred_task, args = wf.add_task(MixMHC2pred(), tag=tumor_sample, depends=[find_novel_peptide_task])
                args['input'].value = find_novel_peptide_task.outputs['MixMHC2pred_faa']
                args['output'].value = tumor_sample + '.MixMHC2pred.txt'
                args['alleles'].value = valid_alleles
            else:
                print(f'skip MixMHC2Pred task for {tumor_sample} since no valid HLA-gene combination available')

            # netMHCIIpan4 prediction
            valid_alleles = check_and_convert_alleles_for_netMHCIIpan4(mhc2_alleles)
            if valid_alleles:
                netmhcpanii_task, args = wf.add_task(netMHCIIPan(), tag=tumor_sample, depends=[find_novel_peptide_task])
                args['alleles'].value = valid_alleles
                args['infile'].value = find_novel_peptide_task.outputs['netMHCPanII_txt']
                args['inptype'].value = '1'
                args['xls'].value = True
                args['xlsfile'].value = tumor_sample + '.netMHCPanII.txt'
                args['stdout'].value = tumor_sample + '.stdout.txt'

        # annotate result of mhcflurry
        annotateMhcflurry, args = wf.add_task(annotate_mhcflurry_result(), tag=tumor_sample, depends=[mhcflurry_task, tumor_assemble_task, tumor_gffcompare_task])
        args['csv_file'].value = mhcflurry_task.outputs['out']
        args['tmap'].value = tumor_assemble_task.outputs['tmap']
        args['gtf'].value = tumor_gffcompare_task.outputs['annotated_gtf']
        args['out'].value = tumor_sample + '.annotated.mhcflurry.csv'

        # annotate result of netMHCpanII
        if netmhcpanii_task is not None:
            annotNetMHCPan, args = wf.add_task(annotate_netMHCPan_result(), tag=tumor_sample,
                                               depends=[netmhcpanii_task, find_novel_peptide_task, tumor_assemble_task, tumor_gffcompare_task])
            args['net_file'].value = netmhcpanii_task.outputs['out']
            args['pep2id_file'].value = find_novel_peptide_task.outputs['mhc2_txt']
            args['tmap'].value = tumor_assemble_task.outputs['tmap']
            args['gtf'].value = tumor_gffcompare_task.outputs['annotated_gtf']
            args['out'].value = tumor_sample + '.annotated.netMHCPanII.txt'

    wf.run()


if __name__ == '__main__':
    pipeline()
