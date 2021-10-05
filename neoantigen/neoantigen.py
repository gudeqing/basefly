import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from nestcmd.nestcmd import Argument, Output, Command, Workflow, TopVar, TmpVar
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
    cmd.runtime.image = 'griffithlab/pvactools:latest'
    cmd.runtime.tool = 'vcf-expression-annotator'
    cmd.args['id-column'] = Argument(prefix='--id-column ', default='transcript', desc='The column header in the expression_file for the column containing gene names/transcript ids.')
    cmd.args['expression-column'] = Argument(prefix='--expression-column ',  default='TPM', desc='The column header in the expression_file for the column containing expression data.')
    cmd.args['sample-name'] = Argument(prefix='--sample-name ',  desc='If the input_vcf contains multiple samples, the name of the sample to annotate')
    cmd.args['output-vcf'] = Argument(prefix='--output-vcf ', desc='Path to write the output VCF file.')
    cmd.args['ignore-ensembl-id-version'] = Argument(prefix='--ignore-ensembl-id-version', type='bool', default=False, desc='Assumes that the final period and number denotes the Ensembl ID version and ignores it ')
    cmd.args['input-vcf'] = Argument(prefix='', type='infile', desc='A VEP-annotated VCF file')
    cmd.args['expression-file'] = Argument(prefix='', type='infile', desc='A TSV file containing expression estimates')
    cmd.args['_x'] = Argument(value='custom', type='fix')
    cmd.args['exp-type'] = Argument(prefix='', default='transcript', range=['gene', 'transcript'], desc='The type of expression data in the expression_file')
    cmd.outputs['out_vcf'] = Output(value='{output-vcf}', type='outfile')
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
    cmd.runtime.image = 'griffithlab/pvactools:latest'
    cmd.runtime.tool = 'pvacseq'
    cmd.args['input_file'] = Argument(prefix='', type='infile', desc="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, Wildtype protein sequence, and Downstream protein sequence information.The VCF may be gzipped (requires tabix index)")
    cmd.args['tumor-sample-name'] = Argument(prefix='', desc='The name of the tumor sample being processed. When processing a multi-sample VCF the sample name must be a sample ID in the input VCF #CHROM header line.')
    cmd.args['allele'] = Argument(prefix='', desc="Name of the allele to use for epitope prediction. Multiple alleles can be specified using a comma-separated list.")
    cmd.args['algorithms'] = Argument(prefix='', array=True, default=['MHCflurry', 'MHCnuggetsI', 'MHCnuggetsII', 'NNalign', 'NetMHC', 'NetMHCIIpan', 'NetMHCcons', 'NetMHCpan', 'PickPocket', 'SMM', 'SMMPMBEC', 'SMMalign'], desc='The epitope prediction algorithms to use. Multiple prediction algorithms can be specified, separated by spaces.')
    cmd.args['output_dir'] = Argument(prefix='', default='.', desc='The directory for writing all result files.')
    cmd.args['e1'] = Argument(prefix='-e1 ', type='int', array=True, delimiter=',', default=[8, 9, 10, 11], desc="CLASS_I_EPITOPE_LENGTH Length of MHC Class I subpeptides (neoepitopes) to predict. Multiple epitope lengths can be specified using a comma-separated list. Typical epitope lengths vary between 8-15.")
    cmd.args['e2'] = Argument(prefix='-e2 ', type='int', array=True, delimiter=',', default=[12, 13, 14, 15, 16, 17, 18], desc="CLASS_II_EPITOPE_LENGTH Length of MHC Class II subpeptides (neoepitopes) to predict. Multiple epitope lengths can be specified using a comma-separated list. Typical epitope lengths vary between 11-30. Required for Class II prediction algorithms.")
    cmd.args['iedb-install-directory'] = Argument(prefix='--iedb-install-directory ', type='indir', level='optional', desc="Directory that contains the local installation of IEDB MHC I and/or MHC II")
    cmd.args['binding-threshold'] = Argument(prefix='--binding-threshold ', default=50, desc='Report only epitopes where the mutant allele has ic50 binding scores below this value')
    cmd.args['percentile-threshold'] = Argument(prefix='--percentile-threshold ', level='optional', desc="Report only epitopes where the mutant allele has a percentile rank below this value")
    cmd.args['allele-specific-binding-thresholds'] = Argument(prefix='--allele-specific-binding-thresholds', type='bool', default=False, desc="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `pvacseq allele_specific_cutoffs`. If an allele does not have a special threshold value, the `--binding-threshold` value will be used.")
    cmd.args['top-score-metric'] = Argument(prefix='--top-score-metric ', default='median', desc="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. lowest:Use the best MT Score and Corresponding Fold Change (i.e. the lowest MT ic50 binding score and corresponding fold change of all chosen prediction methods). median: Use the median MT Score and Median Fold Change (i.e. the median MT ic50 binding score and fold change of all chosen prediction methods).")
    cmd.args['keep-tmp-files'] = Argument(prefix='--keep-tmp-files', type='bool', default=False, desc='Keep intermediate output files')
    cmd.args['threads'] = Argument(prefix='--n-threads ', default=4, desc="Number of threads to use for parallelizing peptide-MHC binding prediction calls")
    cmd.args['netmhc-stab'] = Argument(prefix='--netmhc-stab', type='bool', default=False, desc='Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes.')
    cmd.args['run-reference-proteome-similarity'] = Argument(prefix='--run-reference-proteome-similarity', default=False, desc="Blast peptides against the reference proteome.")
    cmd.args['additional-report-columns'] = Argument(prefix='--additional-report-columns ', level='optional', desc='Additional columns to output in the final report. If sample_name is chosen, this will add a column with the sample name in every row of the output. This can be useful if you later want to concatenate results from multiple individuals into a single file.')
    cmd.args['fasta-size'] = Argument(prefix='--fasta-size ', default=200, desc="Number of FASTA entries per IEDB request. For some resource-intensive prediction algorithms like Pickpocket and NetMHCpan it might be helpful to reduce this number. Needs to be an even number.")
    cmd.args['exclude-NAs'] = Argument(prefix='--exclude-NAs', default=False, desc="Exclude NA values from the filtered output.")
    cmd.args['downstream-sequence-length'] = Argument(prefix='--downstream-sequence-length ', type='str', default='1000', desc="DOWNSTREAM_SEQUENCE_LENGTH Cap to limit the downstream sequence length for frameshifts when creating the FASTA file. Use 'full' to include the full downstream sequence")
    cmd.args['normal-sample-name'] = Argument(prefix='--normal-sample-name ', level='optional', desc="In a multi-sample VCF, the name of the matched normal sample")
    cmd.args['phased-proximal-variants-vcf'] = Argument(prefix='--phased-proximal-variants-vcf ', type='infile', level='optional', desc='A VCF with phased proximal variant information. Must be gzipped and tabix indexed.')
    cmd.args['minimum-fold-change'] = Argument(prefix='--minimum-fold-change ', default=1.0, desc='Minimum fold change between mutant (MT) binding score and wild-type (WT) score (fold change = WT/MT). The default is 0, which filters no results, but 1 is often a sensible choice (requiring that binding is better to the MT than WT peptide). This fold change is sometimes referred to as a differential agretopicity index.')
    cmd.args['normal-cov'] = Argument(prefix='--normal-cov ', default=5, desc='Normal Coverage Cutoff')
    cmd.args['tdna-cov'] = Argument(prefix='--tdna-cov ', default=10, desc='Tumor DNA Coverage Cutoff.')
    cmd.args['trna-cov'] = Argument(prefix='--trna-cov ', default=10, desc='Tumor RNA Coverage Cutoff. Only sites above this read depth cutoff will be considered')
    cmd.args['tdna-vaf'] = Argument(prefix='--tdna-vaf ', default=0.05, desc='Tumor DNA VAF Cutoff. Only sites above this cutoff will be considered')
    cmd.args['trna-vaf'] = Argument(prefix='--trna-vaf ', default=0.1, desc='Tumor RNA VAF Cutoff. Only sites above this cutoff will be considered')
    cmd.args['expn-val'] = Argument(prefix='--expn-val ', default=1.0, desc='Gene and Transcript Expression cutoff. Only sites above this cutoff will be considered.')
    cmd.args['maximum-transcript-support-level'] = Argument(prefix='--maximum-transcript-support-level ', default=1, range=[1,2,3,4,5], desc="The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL). Keep all epitopes with a transcript support level <= to this cutoff")
    cmd.args['pass-only'] = Argument(prefix='--pass-only', type='bool', default=True, desc='Only process VCF entries with a PASS status.')
    cmd.outputs['outdir'] = Output(value='{output_dir}', type='outdir')
    cmd.outputs['all_epitopes'] = Output(value='{output_dir}/*all_epitopes.tsv', type='outfile')
    return cmd


def generate_aggregated_report():
    cmd = Command()
    cmd.meta.name = 'pvacseq'
    cmd.runtime.image = 'griffithlab/pvactools:latest'
    cmd.runtime.tool = 'pvacseq'
    cmd.args['input_file'] = Argument(desc='A pVACseq .all_epitopes.tsv report file')
    cmd.args['output_file'] = Argument(desc='The file path to write the aggregated report tsv to')
    cmd.outputs['output_file'] = Output(value='{output_file}', type='outfile')
    return cmd


def pipeline(inputs=None, outdir='.', run=False, no_docker=False, workers=3, retry=2,
             no_monitor_resource=False, no_check_resource=False):
    wf = Workflow()
    wf.meta.name = 'neoantigen-pipeline'
    wf.meta.desc = 'neoantigen prediction pipeline for multiple samples'

    if inputs is None:
        print('this a simple test')
        input_detail = [['tumor', 'normal', 'HLA-A*02:01,HLA-B*35:01,DRB1*11:01', 'somatic.vcf', 'exp.txt']]
    else:
        input_detail = []
        with open(inputs) as f:
            for line in f:
                # tumor_sample normal_sample alleles somatic_vcf_file_path tumor_expression_txt_file_path
                input_detail.append(line.strip().split()[:5])

    for tumor, normal, alleles, vcf, expr in input_detail:
        # add exp to vcf
        add_exp_task, args = wf.add_task(add_exp_to_vcf(), name=f'addExp-{tumor}')
        args['input-vcf'].value = os.path.abspath(vcf)
        args['expression-file'].value = os.path.abspath(expr)
        args['sample-name'].value = tumor
        args['output-vcf'].value = f'{tumor}.somatic.tx.vcf'

        # pvacseq
        predict_task, args = wf.add_task(pvacseq(), name=f'pvacseq-{tumor}', depends=[add_exp_task.task_id])
        args['input_file'].value = add_exp_task.outputs['out_vcf']
        args['tumor-sample-name'].value = tumor
        args['normal-sample-name'].value = normal
        args['allele'].value = alleles

        # aggregate
        task, args = wf.add_task(generate_aggregated_report(), name=f'aggregate-{tumor}', depends=[predict_task.task_id])
        args['input_file'].value = predict_task.outputs['all_epitopes']
        args['output_file'].value = f'{tumor}.aggregated.report.txt'

    wf.to_nestcmd(outdir=outdir, run=run, no_docker=no_docker, threads=workers, retry=retry,
                  no_monitor_resource=no_monitor_resource, no_check_resource=no_check_resource)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pipeline'])





