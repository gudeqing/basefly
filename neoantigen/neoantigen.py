import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from nestcmd.nestcmd import Argument, Output, Command, Workflow, TopVar, TmpVar
import pandas as pd
__author__ = 'gdq'

"""
1. recommend to input vep annotated Vcf
2. add transcript expression info to Vcf
3. run pvacseq
4. aggregate
"""


def add_exp_to_vcf():
    cmd = Command()
    cmd.meta.name = 'vcf-expression-annotator'
    cmd.runtime.image = 'gudeqing/neoantigen:1.0'
    cmd.runtime.tool = 'vcf-expression-annotator'
    cmd.args['id-column'] = Argument(prefix='--id-column ', default='Name', desc='The column header in the expression_file for the column containing gene names/transcript ids.')
    cmd.args['expression-column'] = Argument(prefix='--expression-column ',  default='TPM', desc='The column header in the expression_file for the column containing expression data.')
    cmd.args['sample-name'] = Argument(prefix='--sample-name ',  desc='If the input_vcf contains multiple samples, the name of the sample to annotate')
    cmd.args['output-vcf'] = Argument(prefix='--output-vcf ', desc='Path to write the output VCF file.')
    cmd.args['ignore-ensembl-id-version'] = Argument(prefix='--ignore-ensembl-id-version', type='bool', default=True, desc='Assumes that the final period and number denotes the Ensembl ID version and ignores it ')
    cmd.args['input-vcf'] = Argument(prefix='', type='infile', desc='A VEP-annotated VCF file')
    cmd.args['expression-file'] = Argument(prefix='', type='infile', desc='A TSV file containing expression estimates')
    cmd.args['_x'] = Argument(value='custom', type='fix')
    cmd.args['exp-type'] = Argument(prefix='', default='gene', range=['gene', 'transcript'], desc='The type of expression data in the expression_file')
    cmd.outputs['output-vcf'] = Output(value='{output-vcf}', type='outfile')
    return cmd


def pvacseq():
    """
    usage: pvacseq run [-h] [-e1 CLASS_I_EPITOPE_LENGTH]
           [-e2 CLASS_II_EPITOPE_LENGTH]
           [--iedb-install-directory IEDB_INSTALL_DIRECTORY]
           [-b BINDING_THRESHOLD]
           [--percentile-threshold PERCENTILE_THRESHOLD]
           [--allele-specific-binding-thresholds] [-m {lowest,median}]
           [-r IEDB_RETRIES] [-k] [-t N_THREADS]
           [--net-chop-method {cterm,20s}] [--netmhc-stab]
           [--net-chop-threshold NET_CHOP_THRESHOLD]
           [--run-reference-proteome-similarity] [-a {sample_name}]
           [-s FASTA_SIZE] [--exclude-NAs]
           [-d DOWNSTREAM_SEQUENCE_LENGTH]
           [--normal-sample-name NORMAL_SAMPLE_NAME]
           [-p PHASED_PROXIMAL_VARIANTS_VCF] [-c MINIMUM_FOLD_CHANGE]
           [--normal-cov NORMAL_COV] [--tdna-cov TDNA_COV]
           [--trna-cov TRNA_COV] [--normal-vaf NORMAL_VAF]
           [--tdna-vaf TDNA_VAF] [--trna-vaf TRNA_VAF]
           [--expn-val EXPN_VAL]
           [--maximum-transcript-support-level {1,2,3,4,5}]
           [--pass-only]
           input_file sample_name allele
           {MHCflurry,MHCnuggetsI,MHCnuggetsII,NNalign,NetMHC,NetMHCIIpan,NetMHCcons,NetMHCpan,PickPocket,SMM,SMMPMBEC,SMMalign,all,all_class_i,all_class_ii}
           [{MHCflurry,MHCnuggetsI,MHCnuggetsII,NNalign,NetMHC,NetMHCIIpan,NetMHCcons,NetMHCpan,PickPocket,SMM,SMMPMBEC,SMMalign,all,all_class_i,all_class_ii} ...]
           output_dir
    :return:
    """
    cmd = Command()
    cmd.meta.name = 'pvacseq'
    cmd.meta.version = '2.0.4'
    cmd.runtime.image = 'griffithlab/pvactools:latest'
    # cmd.runtime.image = 'gudeqing/neoantigen:1.0'
    cmd.runtime.tool = 'pvacseq run'
    cmd.args['input_file'] = Argument(prefix='', type='infile', desc="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, Wildtype protein sequence, and Downstream protein sequence information.The VCF may be gzipped (requires tabix index)")
    cmd.args['tumor-sample-name'] = Argument(prefix='', desc='The name of the tumor sample being processed. When processing a multi-sample VCF the sample name must be a sample ID in the input VCF #CHROM header line.')
    cmd.args['allele'] = Argument(prefix='', desc="Name of the allele to use for epitope prediction. Multiple alleles can be specified using a comma-separated list.")
    cmd.args['algorithms'] = Argument(prefix='', array=True, default=['MHCflurry', 'MHCnuggetsII'], desc="The epitope prediction algorithms to use. Multiple prediction algorithms can be specified, separated by spaces. ['MHCflurry', 'MHCnuggetsI', 'MHCnuggetsII', 'NNalign', 'NetMHC', 'NetMHCIIpan', 'NetMHCcons', 'NetMHCpan', 'PickPocket', 'SMM', 'SMMPMBEC', 'SMMalign']")
    cmd.args['output_dir'] = Argument(prefix='', default='.', desc='The directory for writing all result files.')
    cmd.args['e1'] = Argument(prefix='-e1 ', type='int', array=True, delimiter=',', default=[8, 9, 10, 11], desc="CLASS_I_EPITOPE_LENGTH Length of MHC Class I subpeptides (neoepitopes) to predict. Multiple epitope lengths can be specified using a comma-separated list. Typical epitope lengths vary between 8-15.")
    cmd.args['e2'] = Argument(prefix='-e2 ', type='int', array=True, delimiter=',', default=[12, 13, 14, 15, 16, 17, 18], desc="CLASS_II_EPITOPE_LENGTH Length of MHC Class II subpeptides (neoepitopes) to predict. Multiple epitope lengths can be specified using a comma-separated list. Typical epitope lengths vary between 11-30. Required for Class II prediction algorithms.")
    cmd.args['iedb-install-directory'] = Argument(prefix='--iedb-install-directory ', default='/opt/iedb/', type='indir', level='optional', desc="Directory that contains the local installation of IEDB MHC I and/or MHC II")
    cmd.args['binding-threshold'] = Argument(prefix='--binding-threshold ', default=500, desc='Report only epitopes where the mutant allele has ic50 binding scores below this value')
    cmd.args['percentile-threshold'] = Argument(prefix='--percentile-threshold ', level='optional', desc="Report only epitopes where the mutant allele has a percentile rank below this value")
    cmd.args['allele-specific-binding-thresholds'] = Argument(prefix='--allele-specific-binding-thresholds', type='bool', default=False, desc="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `pvacseq allele_specific_cutoffs`. If an allele does not have a special threshold value, the `--binding-threshold` value will be used.")
    cmd.args['top-score-metric'] = Argument(prefix='--top-score-metric ', default='median', desc="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. lowest:Use the best MT Score and Corresponding Fold Change (i.e. the lowest MT ic50 binding score and corresponding fold change of all chosen prediction methods). median: Use the median MT Score and Median Fold Change (i.e. the median MT ic50 binding score and fold change of all chosen prediction methods).")
    cmd.args['keep-tmp-files'] = Argument(prefix='--keep-tmp-files', type='bool', default=False, desc='Keep intermediate output files')
    cmd.args['threads'] = Argument(prefix='--n-threads ', default=5, desc="Number of threads to use for parallelizing peptide-MHC binding prediction calls")
    cmd.args['netmhc-stab'] = Argument(prefix='--netmhc-stab', type='bool', default=False, desc='Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes.')
    cmd.args['run-reference-proteome-similarity'] = Argument(prefix='--run-reference-proteome-similarity', default=False, desc="Blast peptides against the reference proteome.")
    cmd.args['additional-report-columns'] = Argument(prefix='--additional-report-columns ', level='optional', default='sample_name', desc='Additional columns to output in the final report. If sample_name is chosen, this will add a column with the sample name in every row of the output. This can be useful if you later want to concatenate results from multiple individuals into a single file.')
    cmd.args['fasta-size'] = Argument(prefix='--fasta-size ', default=200, desc="Number of FASTA entries per IEDB request. For some resource-intensive prediction algorithms like Pickpocket and NetMHCpan it might be helpful to reduce this number. Needs to be an even number.")
    cmd.args['exclude-NAs'] = Argument(prefix='--exclude-NAs', default=False, desc="Exclude NA values from the filtered output.")
    cmd.args['downstream-sequence-length'] = Argument(prefix='--downstream-sequence-length ', type='str', default='1000', desc="DOWNSTREAM_SEQUENCE_LENGTH Cap to limit the downstream sequence length for frameshifts when creating the FASTA file. Use 'full' to include the full downstream sequence")
    cmd.args['normal-sample-name'] = Argument(prefix='--normal-sample-name ', level='optional', desc="In a multi-sample VCF, the name of the matched normal sample")
    cmd.args['phased-proximal-variants-vcf'] = Argument(prefix='--phased-proximal-variants-vcf ', type='infile', level='optional', desc='A VCF with phased proximal variant information. Must be gzipped and tabix indexed.')
    cmd.args['minimum-fold-change'] = Argument(prefix='--minimum-fold-change ', default=1.0, desc='Minimum fold change between mutant (MT) binding score and wild-type (WT) score (fold change = WT/MT). The default is 0, which filters no results, but 1 is often a sensible choice (requiring that binding is better to the MT than WT peptide). This fold change is sometimes referred to as a differential agretopicity index.')
    cmd.args['normal-cov'] = Argument(prefix='--normal-cov ', default=5, desc='Normal Coverage Cutoff')
    cmd.args['tdna-cov'] = Argument(prefix='--tdna-cov ', default=15, desc='Tumor DNA Coverage Cutoff.')
    cmd.args['trna-cov'] = Argument(prefix='--trna-cov ', default=15, desc='Tumor RNA Coverage Cutoff. Only sites above this read depth cutoff will be considered')
    cmd.args['tdna-vaf'] = Argument(prefix='--tdna-vaf ', default=0.05, desc='Tumor DNA VAF Cutoff. Only sites above this cutoff will be considered')
    cmd.args['trna-vaf'] = Argument(prefix='--trna-vaf ', default=0.0, desc='Tumor RNA VAF Cutoff. Only sites above this cutoff will be considered')
    cmd.args['expn-val'] = Argument(prefix='--expn-val ', default=1.0, desc='Gene and Transcript Expression cutoff. Only sites above this cutoff will be considered.')
    cmd.args['maximum-transcript-support-level'] = Argument(prefix='--maximum-transcript-support-level ', default=1, range=[1,2,3,4,5], desc="The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL). Keep all epitopes with a transcript support level <= to this cutoff")
    cmd.args['pass-only'] = Argument(prefix='--pass-only', type='bool', default=True, desc='Only process VCF entries with a PASS status.')
    cmd.outputs['outdir'] = Output(value='{output_dir}', type='outdir')
    cmd.outputs['all_epitopes'] = Output(value='{output_dir}/*.all_epitopes.tsv', type='outfile')
    cmd.outputs['aggregated_epitopes'] = Output(value='{output_dir}/*.aggregated.tsv', type='outfile')
    return cmd


def generate_aggregated_report():
    cmd = Command()
    cmd.meta.name = 'pvacseq'
    cmd.runtime.image = 'griffithlab/pvactools:latest'
    cmd.runtime.tool = 'pvacseq generate_aggregated_report'
    cmd.args['input_file'] = Argument(desc='A pVACseq .all_epitopes.tsv report file')
    cmd.args['output_file'] = Argument(desc='The file path to write the aggregated report tsv to')
    cmd.outputs['output_file'] = Output(value='{output_file}', type='outfile')
    return cmd


def get_2digits_hla_genetype(table, sample, alleles):
    """prepare for pvacseq tool"""
    df = pd.read_table(table, index_col=0)
    targets = list()
    for each in df.loc[sample]:
        a, b = each.split('|')[:2]
        if a.startswith(alleles):
            allele_2_digit = ':'.join(a.split(':')[:2])
            if allele_2_digit not in targets:
                if allele_2_digit.startswith(('A*', 'B*', 'C*', 'E*', 'F*', 'G*')):
                    allele_2_digit = 'HLA-' + allele_2_digit
                targets.append(allele_2_digit)
        if b.startswith(alleles):
            allele_2_digit = ':'.join(b.split(':')[:2])
            if allele_2_digit not in targets:
                if allele_2_digit.startswith(('A*', 'B*', 'C*', 'E*', 'F*', 'G*')):
                    allele_2_digit = 'HLA-' + allele_2_digit
                targets.append(allele_2_digit)
    return targets


def phasing_vcf():
    # hapcut2 to phase unzipped vcf with bam
    pass


def combine_germline_somatic_vcf():
    # bcftools
    pass


def pipeline():
    wf = Workflow()
    wf.meta.name = 'neoantigen-pipeline'
    wf.meta.desc = 'neoantigen prediction pipeline for multiple samples'
    wf.init_argparser()
    wf.add_argument('-gene_expr', required=True, help='gene expression table with header consisting of tumor sample name')
    wf.add_argument('-trans_expr', required=False, help='transcript expression table with header consisting of tumor sample name')
    wf.add_argument('-vcfs', required=True, help='somatic vcf list file with two columns, first column is tumor sample name, second column is vcf path')
    wf.add_argument('-pair_info', required=True, help='tumor-normal pair info without header, first column is tumor sample name')
    wf.add_argument('-hla_genotype', required=True, help='HLA genetype table with header consiting of HLA gene name, index is tumor sample name. Each element data is genetype detail.')
    wf.add_argument('-alleles', default=('A', 'B', 'C', 'DRA', 'DRB', 'DQA', 'DQB'), nargs='+', help='HLA alleles gene list, such as A,B,C,DRA(B),DQA(B),DPA(B)')
    wf.parse_args()

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
                lst = line.strip().split('\t')
                vcf_dict[lst[0]] = lst[1]

    for tumor, normal in pair_list:
        # add gene exp to vcf
        add_exp_task, args = wf.add_task(add_exp_to_vcf(), name=f'addGeneExp-{tumor}')
        args['input-vcf'].value = os.path.abspath(vcf_dict[tumor])
        args['expression-file'].value = os.path.abspath(wf.args.gene_expr)
        args['sample-name'].value = tumor
        args['expression-column'].value = tumor
        args['exp-type'].value = 'gene'
        args['id-column'].value = 'Name'
        args['output-vcf'].value = f'{tumor}.somatic.gx.vcf'
        out_vcf = add_exp_task.outputs['output-vcf']
        if wf.args.trans_expr:
            add_exp_task, args = wf.add_task(add_exp_to_vcf(), name=f'addTransExp-{tumor}', depends=[add_exp_task.task_id])
            args['input-vcf'].value = out_vcf
            args['expression-file'].value = os.path.abspath(wf.args.trans_expr)
            args['sample-name'].value = tumor
            args['expression-column'].value = tumor
            args['id-column'].value = 'Name'
            args['exp-type'].value = 'transcript'
            args['output-vcf'].value = f'{tumor}.somatic.gx.tx.vcf'

        # pvacseq
        predict_task, args = wf.add_task(pvacseq(), name=f'pvacseq-{tumor}', depends=[add_exp_task.task_id])
        args['input_file'].value = add_exp_task.outputs['output-vcf']
        args['tumor-sample-name'].value = tumor
        args['normal-sample-name'].value = normal
        # use HLA gene-type from tumor or normal??
        args['allele'].value = ','.join(get_2digits_hla_genetype(wf.args.hla_genotype, tumor, alleles=wf.args.alleles))

    wf.run()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pipeline'])





