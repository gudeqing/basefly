import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from nestcmd.nestcmd import Workflow, Argument, Output, Command,  TopVar, TmpVar
from nestcmd.commands import add_exp_to_vcf, bam_read_count, add_read_count_to_vcf, pvacseq
from utils.tidy_tools import get_2digits_hla_genetype
import pandas as pd
__author__ = 'gdq'

"""
1. recommend to input vep annotated Vcf
2. add transcript expression info to Vcf
3. run pvacseq
当前测试发现bam_read_count, add_read_count_to_vcf两个步骤似乎不能按照预期工作，暂不建议使用
"""


def pipeline():
    wf = Workflow()
    wf.meta.name = 'neoantigen-pipeline'
    wf.meta.desc = 'neoantigen prediction pipeline for multiple samples using pvactools'
    wf.init_argparser()
    wf.add_argument('-gene_expr', required=True, help='gene expression table with header consisting of tumor sample name')
    wf.add_argument('-trans_expr', required=False, help='transcript expression table with header consisting of tumor sample name')
    wf.add_argument('-vcfs', required=True, help='somatic vcf list file with two columns, first column is tumor sample name, second column is vcf path')
    wf.add_argument('-pair_info', required=True, help='tumor-normal pair info without header, first column is tumor sample name')
    wf.add_argument('-hla_genotype', required=True, help='HLA genetype table with header consiting of HLA gene name, index is tumor sample name. Each element data is genetype detail.')
    wf.add_argument('-alleles', default=('A', 'B', 'C', 'DRA', 'DRB', 'DQA', 'DQB'), nargs='+', help='HLA alleles gene list, such as A,B,C,DRA(B),DQA(B),DPA(B)')
    wf.add_argument('-rna_bams', required=False, help='rnaseq bam list file with two columns, first column is normal sample name, second column is bam path. This bam will be used to add RNA coverage and RNA VAF using bam_readcount.')
    wf.add_argument('-ref_fasta', required=False, help='reference fasta file, which is need when rna_bams provided or fastq info provided')
    wf.parse_args()

    sample_has_expr = []
    with open(wf.args.gene_expr) as f:
        sample_has_expr = f.readline().strip().split()

    pair_list = []
    if wf.args.pair_info:
        with open(wf.args.pair_info) as f:
            for line in f:
                if line.strip():
                    pairs = line.strip('\n').split('\t')[:2]
                    pair_list.append(pairs)

    vcf_dict = dict()
    if wf.args.vcfs:
        with open(wf.args.vcfs) as f:
            for line in f:
                lst = line.strip().split('\t')[:2]
                vcf_dict[lst[0]] = lst[1]

    bam_dict = dict()
    if wf.args.rna_bams:
        with open(wf.args.rna_bams) as f:
            for line in f:
                lst = line.strip().split('\t')[:2]
                bam_dict[lst[0]] = lst[1]

    no_exp_samples = set(x[0] for x in pair_list) - set(sample_has_expr)
    if no_exp_samples:
        print('these samples are not in expression matrix', no_exp_samples)

    for tumor, normal in pair_list:
        tumor_exp_sample_name = tumor
        if tumor not in sample_has_expr:
            tmp_dict = {'P48_L': 'P48_L1', 'P31_L': 'P31_L1'}
            if tumor in tmp_dict:
                tumor_exp_sample_name = tmp_dict[tumor]
            else:
                print(f'{tumor} is not in expression matrix and the pair({tumor} vs {normal}) will be skipped for neoantigen analysis')
        # add gene exp to vcf
        add_exp_task, args = wf.add_task(add_exp_to_vcf(), name=f'addGeneExp-{tumor}')
        args['input-vcf'].value = os.path.abspath(vcf_dict[tumor])
        args['expression-file'].value = os.path.abspath(wf.args.gene_expr)
        args['sample-name'].value = tumor
        args['expression-column'].value = tumor_exp_sample_name
        args['exp-type'].value = 'gene'
        args['id-column'].value = 'Name'
        args['output-vcf'].value = f'{tumor}.somatic.gx.vcf'
        out_vcf = add_exp_task.outputs['output-vcf']
        if wf.args.trans_expr:
            add_exp_task, args = wf.add_task(add_exp_to_vcf(), name=f'addTransExp-{tumor}', depends=[add_exp_task.task_id])
            args['input-vcf'].value = out_vcf
            args['expression-file'].value = os.path.abspath(wf.args.trans_expr)
            args['sample-name'].value = tumor
            args['expression-column'].value = tumor_exp_sample_name
            args['id-column'].value = 'Name'
            args['exp-type'].value = 'transcript'
            args['output-vcf'].value = f'{tumor}.somatic.gx.tx.vcf'

        add_readcount_indel = None
        if bam_dict:
            readcount_task, args = wf.add_task(bam_read_count(), name=f'bamReadCount-{tumor}')
            args['vcf'].value = os.path.abspath(vcf_dict[tumor])
            args['bam'].value = os.path.abspath(bam_dict[tumor])
            args['sample'].value = tumor
            if not wf.args.ref_fasta:
                raise Exception('ref_fasta is not provided!')
            args['ref_fasta'].value = os.path.abspath(wf.args.ref_fasta)

            add_readcount_snv, args = wf.add_task(add_read_count_to_vcf(), name=f'addRnaSNV-{tumor}', depends=[readcount_task.task_id, add_exp_task.task_id])
            args['sample'].value = tumor
            args['variant_type'].value = 'snv'
            args['in_vcf'].value = add_exp_task.outputs['output-vcf']
            args['out_vcf'].value = f'{tumor}.somatic.withRNAsnv.vcf.gz'
            args['read_count_file'].value = readcount_task.outputs['snv_readcount']

            add_readcount_indel, args = wf.add_task(add_read_count_to_vcf(), name=f'addRnaIndel-{tumor}', depends=[readcount_task.task_id, add_readcount_snv.task_id])
            args['sample'].value = tumor
            args['variant_type'].value = 'indel'
            args['in_vcf'].value = add_readcount_snv.outputs['output-vcf']
            args['out_vcf'].value = f'{tumor}.somatic.final.vcf.gz'
            args['read_count_file'].value = readcount_task.outputs['indel_readcount']

        # pvacseq
        depend_task = add_readcount_indel or add_exp_task
        predict_task, args = wf.add_task(pvacseq(), name=f'pvacseq-{tumor}', depends=[depend_task.task_id])
        args['input_file'].value = depend_task.outputs['output-vcf']
        args['tumor-sample-name'].value = tumor
        args['normal-sample-name'].value = normal
        # use HLA gene-type from tumor or normal??
        args['allele'].value = ','.join(get_2digits_hla_genetype(wf.args.hla_genotype, normal, alleles=wf.args.alleles))

    wf.run()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pipeline'])
