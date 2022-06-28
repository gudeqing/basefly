import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, ToWdlTask
from utils.get_fastq_info import get_fastq_info
__author__ = 'gdq'


def fastp():
    cmd = Command()
    cmd.meta.name = 'fastp'
    cmd.runtime.image = 'gudeqing/fastp:latest'
    cmd.runtime.tool = 'fastp'
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', level='optional', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=7, desc='thread number')
    cmd.args['other_args'] = Argument(prefix='', default='', desc="other arguments you want to use, such as '-x val'")
    cmd.args['out1'] = Argument(prefix='-o ', type='str', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', level='optional', type='str', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='-h ', type='str', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', type='str', desc='json report file')
    # 定义输出
    cmd.outputs['out1'] = Output(value="{out1}")  # 这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd


def bismark_genome_preparation():
    cmd = Command()
    cmd.meta.name = 'GenomePreparation'
    cmd.meta.version = '0.23.1'
    cmd.runtime.image = 'chgyi/meth:latest'
    cmd.runtime.tool = 'bismark_genome_preparation'
    cmd.args['aligner'] = Argument(prefix='--path_to_aligner ', level='optional', desc='The full path to the Bowtie 2 or HISAT2 installation folder on your system (depending on which aligner/indexer you intend to use; please note that this is the folder and not any executable). Unless this path is specified, it is assumed that the aligner in question (Bowtie 2/HISAT2) is in the PATH')
    cmd.args['index_type'] = Argument(prefix='--', default='bowtie2', range=['bowtie2', 'hisat2'], desc='Bowtie2 is recommended for most bisulfite sequencing applications. For hisat2, only recommended for specialist applications such as RNA-methylation analyses or SLAM-seq type applications')
    cmd.args['threads'] = Argument(prefix='--parallel ', type='int', default=4, desc='Use several threads for each indexing process to speed up the genome preparation step.Remember that the indexing is run twice in parallel already (for the top and bottom strand separately), so e.g. "--parallel 4" will use 8 threads in total')
    cmd.args['genome_folder'] = Argument(prefix='', type='indir', desc='The path to the folder containing the genome to be bisulfite converted. The Bismark Genome Preparation expects one or more fastA files in the folder (with the file extension: .fa or .fasta (also ending in .gz)). Specifying this path is mandatory')
    cmd.args['slam'] = Argument(prefix='--slam', type='bool', default=False, desc='Instead of performing an in-silico bisulfite conversion, this mode transforms T to C (forward strand),or A to G (reverse strand).')
    cmd.outputs['outdir'] = Output(value='{genome_folder}', type='outdir')
    return cmd


def bismark_alignment():
    cmd = Command()
    cmd.meta.name = 'BismarkAlignment'
    cmd.meta.version = '0.23.1'
    cmd.runtime.image = 'chgyi/meth:latest'
    cmd.runtime.tool = 'bismark'
    cmd.args['output_dir'] = Argument(prefix='-o ', default='.', desc='Write all output files into this directory')
    cmd.args['temp_dir'] = Argument(prefix='--temp_dir ', default='.', desc='Write temporary files to this directory instead of into the same directory as the input files')
    cmd.args['prefix'] = Argument(prefix='--prefix ', level='optional', desc="Prefixes <prefix> to the output filenames. Trailing dots will be replaced by a single one. For example, '--prefix test' with 'file.fq' would result in the output file 'test.file.fq_bismark.sam' etc.")
    cmd.args['nucleotide_coverage'] = Argument(prefix='--nucleotide_coverage', type='bool', default=False, desc="Calculates the mono- and di-nucleotide sequence composition of covered positions in the analysed BAM file and compares it to the genomic average composition once alignments are complete by calling 'bam2nuc'.")
    cmd.args['genome_folder'] = Argument(prefix='--genome_folder ', type='indir', desc="The path to the folder containing the unmodified reference genome as well as the subfolders created by the Bismark_Genome_Preparation script (/Bisulfite_Genome/CT_conversion/ and /Bisulfite_Genome/GA_conversion/). Bismark expects one or more fastA files in this folder (file extension: .fa, .fa.gz or .fasta or .fasta.gz). The path can be relative or absolute.")
    cmd.args['read1'] = Argument(prefix='-1 ', type='infile', level='optional', array=True, delimiter=',', desc='Comma-separated list of files containing the first sequencing read')
    cmd.args['read2'] = Argument(prefix='-2 ', type='infile', level='optional', array=True, delimiter=',', desc='Comma-separated list of files containing the second sequencing read')
    cmd.args['singles'] = Argument(prefix='', type='infile', level='optional', array=True, delimiter=',', desc='A comma-separated list of files containing the reads to be aligned')
    cmd.outputs['bam'] = Output(value='*.bam')
    cmd.outputs['report'] = Output(value='*_report.txt*')
    cmd.outputs['nucleotide_stats'] = Output(value='*nucleotide_stats.txt*')
    return cmd


def bismark_deduplicate():
    desc = """
    This command will deduplicate the Bismark alignment BAM file and
     remove all reads but one which align to the the very same position and in the same orientation. 
    This step is recommended for whole-genome bisulfite samples, 
    but should not be used for reduced representation libraries such as RRBS, amplicon or target enrichment libraries.
    """
    cmd = Command()
    cmd.meta.desc = desc
    cmd.meta.name = 'BismarkDeduplicate'
    cmd.meta.version = '0.23.1'
    cmd.runtime.image = 'chgyi/meth:latest'
    cmd.runtime.tool = 'deduplicate_bismark'
    cmd.args['prefix'] = Argument(prefix='-o ', level='optional', desc="The basename of a desired output file. This basename is modified to end in '.deduplicated.bam', or '.multiple.deduplicated.bam' in '--multiple' mode, for consistency reasons.")
    cmd.args['multiple'] = Argument(prefix='--multiple', type='bool', default=False, desc='All specified input files are treated as one sample and concatenated together for deduplication.')
    cmd.args['outdir'] = Argument(prefix='-output_dir ', default='.', desc="Output directory, either relative or absolute. Output is written to the current directory if not specified explicitly")
    cmd.args['bam'] = Argument(prefix='', type='infile', array=True, desc='bam file generated by bismark alignment')
    cmd.outputs['bam'] = Output(value='*.deduplicated.bam')
    cmd.outputs['report'] = Output(value='*_report.txt*')
    return cmd


def bismark_methylation_extractor():
    cmd = Command()
    cmd.meta.name = 'MethylationExtractor'
    cmd.meta.version = '0.23.1'
    cmd.runtime.image = 'chgyi/meth:latest'
    cmd.runtime.tool = 'bismark_methylation_extractor'
    cmd.args['ignore'] = Argument(prefix='--ignore ', type='int', level='optional', desc="Ignore the first <int> bp from the 5' end of Read 1 (or single-end alignment; files) when processing the methylation call string. This can remove e.g. a restriction enzyme site at the start of each read or any other source of bias (such as PBAT-Seq data).")
    cmd.args['ignore_r2'] = Argument(prefix='--ignore_r2 ', type='int', level='optional', desc="Ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing results only. Since the first couple of bases in Read 2 of BS-Seq experiments show a severe bias towards non-methylation as a result of end-repairing sonicated fragments with unmethylated cytosines (see M-bias plot), it is recommended that the first couple of bp of Read 2 are removed before starting downstream analysis. Please see the section on M-bias plots in the Bismark User Guide for more details.")
    cmd.args['ignore_3prime'] = Argument(prefix='--ignore_3prime ', type='int', level='optional', desc="Ignore the last <int> bp from the 3' end of Read 1 (or single-end alignment files) when processing the methylation call string. This can remove unwanted biases from the end of reads.")
    cmd.args['ignore_3prime_r2'] = Argument(prefix='--ignore_3prime_r2 ', type='int', level='optional', desc="Ignore the last <int> bp from the 3' end of Read 2 of paired-end sequencing results only. This can remove unwanted biases from the end of reads.")
    cmd.args['comprehensive'] = Argument(prefix='--comprehensive', type='bool', default=False, desc="Specifying this option will merge all four possible strand-specific methylation info into context-dependent output files.")
    cmd.args['no_header'] = Argument(prefix='--no_header', type='bool', default=False, desc="Suppresses the Bismark version header line in all output files for more convenient batch processing.")
    cmd.args['gzip'] = Argument(prefix='--gzip', type='bool', default=True, desc="The methylation extractor files (CpG_OT_..., CpG_OB_... etc) will be written out in a GZIP compressed form to save disk space.")
    cmd.args['parallel'] = Argument(prefix='--parallel ', type='int', default=8, desc="Sets the number of cores to be used for the methylation extraction process. If system resources are plentiful this is a viable option to speed up the extraction process (we observed a near linear speed increase for up to 10 cores used). Please note that a typical process of extracting a BAM file and writing out '.gz' output streams will in fact use ~3 cores per value of --parallel <int> specified (1 for the methylation extractor itself, 1 for a Samtools stream, 1 for GZIP stream), so --parallel 10 is likely to use around 30 cores of system resources. This option has no bearing on the bismark2bedGraph or genome-wide cytosine report processes.")
    cmd.args['ucsc'] = Argument(prefix='--ucsc', type='bool', default=False, desc="Writes out an additional bedGraph file (ending in '_UCSC.bedGraph.gz') that is compatible with the UCSC genome browser. If alignments were carried out against an Ensembl version of the genome, this step will prefix chromosome names with 'chr', and rename the mitochondrial chromosome from 'MT' to 'chrM'. In addition, this wite out a tab delimited file ending in '.chromosome_sizes.txt' to enable tools such as bedGraphToBigWig.")
    cmd.args['yacht'] = Argument(prefix='--yacht', type='bool', default=False, desc="This option (Yet Another Context Hunting Tool) writes out additional information about the read a methylation call belongs to, and its output is meant to be fed into the NOMe_filtering script. This option writes out a single 'any_C_context' file that contains all methylation calls for a read consecutively. Its intended use is single-cell NOMe-Seq data, and thus this option works only in single-end mode (paired-end reads often suffer from chimaera problems...)")
    cmd.args['bedGraph'] = Argument(prefix='--bedGraph', type='bool', default=True, desc="After finishing the methylation extraction, the methylation output is written into a sorted bedGraph file that reports the position of a given cytosine and its methylation state (in %, see details below).")
    cmd.args['outdir'] = Argument(prefix='-o ', default='.', desc="Allows specification of a different output directory.If not specified explicitly, the output will be written to the current directory")
    cmd.args['bam'] = Argument(prefix='', type='infile', array=True, desc="A space-separated list of Bismark result files in SAM format from which methylation information is extracted for every cytosine in the reads.")
    cmd.outputs['outdir'] = Output(value='{outdir}', type='outdir')
    cmd.outputs['report'] = Output(value='*_report.txt*')
    cmd.outputs['mbias'] = Output(value='*.M-bias.txt*')
    return cmd


def bismark2report():
    cmd = Command()
    cmd.meta.name = 'bismark2report'
    cmd.meta.version = '0.23.1'
    cmd.runtime.image = 'chgyi/meth:latest'
    cmd.runtime.tool = 'bismark2report'
    cmd.args['outdir'] = Argument(prefix='--dir ', default='.', desc="Output directory. Output is written to the current directory if not specified explicitly.")
    cmd.args['alignment_report'] = Argument(prefix='--alignment_report ', type='infile', desc="Bismark alignment report file")
    cmd.args['dedup_report'] = Argument(prefix="--dedup_report ", level='optional', type='infile', desc='deduplication report file')
    cmd.args['splitting_report'] = Argument(prefix="--splitting_report ", level='optional', type='infile', desc='splitting report file (generated by the Bismark methylation extractor)')
    cmd.args['mbias_report'] = Argument(prefix="--mbias_report ", level='optional', type='infile', desc='a single M-bias report file (generated by the Bismark methylation extractor)')
    cmd.args['nucleotide_report'] = Argument(prefix="--nucleotide_report ", level='optional', type='infile', desc="a single nucleotide coverage report file (generated by Bismark with the option '--nucleotide_coverage')")
    cmd.outputs['report_html'] = Output(value='{outdir}/*_report.html')
    return cmd


def bismark2summary():
    cmd = Command()
    cmd.meta.name = 'bismark2summary'
    cmd.meta.version = '0.23.1'
    cmd.runtime.image = 'chgyi/meth:latest'
    cmd.runtime.tool = 'bismark2summary'
    cmd.args['basename'] = Argument(prefix='-o ', default="bismark_summary_report", desc="Basename of the output file (optional)")
    cmd.args['bams'] = Argument(prefix="", array=True, type='infile', desc='Bismark alignment files')
    cmd.outputs['summary_txt'] = Output(value='{basename}.txt')
    cmd.outputs['summary_html'] = Output(value='{basename}.html')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'Bismark'
    wf.meta.source = "https://github.com/FelixKrueger/Bismark"
    wf.meta.desc = """
    Bismark is a program to map bisulfite treated sequencing reads to a genome of interest and 
    perform methylation calls in a single step. The output can be easily imported into a genome viewer, 
    such as SeqMonk, and enables a researcher to analyse the methylation levels of their samples straight away. 
    It's main features are:
    * Bisulfite mapping and methylation calling in one single step
    * Supports single-end and paired-end read alignments
    * Supports ungapped, gapped or spliced alignments
    * Alignment seed length, number of mismatches etc. are adjustable
    * Output discriminates between cytosine methylation in CpG, CHG and CHH context
    """
    wf.meta.version = "0.23.1"

    # 定义流程输入参数
    wf.init_argparser()
    wf.add_argument('-fastq_info', nargs='+', required=True,
                    help='A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]')
    wf.add_argument('-r1_name', default='(.*).R1.fastq',
                    help="python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'")
    wf.add_argument('-r2_name', default='(.*).R2.fastq',
                    help="python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'")
    wf.add_argument('-exclude_samples', default=tuple(), nargs='+', help='samples to exclude from analysis')
    wf.add_argument('-genome_folder', help='The path to the folder containing the genome to be bisulfite converted.')
    wf.add_argument('--dedup', action='store_true', default=False, help='This step is recommended for whole-genome bisulfite samples, but should not be used for reduced representation libraries such as RRBS, amplicon or target enrichment libraries.')
    wf.parse_args()

    # 串联任务
    # 提取fastq信息
    fastq_info = get_fastq_info(fastq_info=wf.args.fastq_info, r1_name=wf.args.r1_name, r2_name=wf.args.r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found !')

    # 建索引
    index_dirs = os.listdir(wf.args.genome_folder)
    if 'Bisulfite_Genome' in index_dirs:
        print('Bisulfite_Genome already exist, we will skip bismark_genome_preparation')
        index_task_id = []
    else:
        index_task, args = wf.add_task(bismark_genome_preparation(), name='GenomePreparation')
        args['genome_folder'].value = wf.args.genome_folder
        index_task_id = [index_task.task_id]

    # fastq预处理和比对
    for sample, reads in fastq_info.items():
        if sample in wf.args.exclude_samples:
            continue
        if len(reads) == 2:
            r1s, r2s = reads
        else:
            r1s = reads[0]
            r2s = [None]*len(r1s)
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
        align_task, args = wf.add_task(bismark_alignment(), name='BismarkAlignment-' + sample, depends=fastp_task_ids + index_task_id)
        args['genome_folder'].value = wf.args.genome_folder
        if r2s[0] is not None:
            args['read1'].value = [x.outputs["out1"] for x in fastp_tasks]
            args['read2'].value = [x.outputs["out2"] for x in fastp_tasks]
        else:
            args['singles'].value = [x.outputs["out1"] for x in fastp_tasks]

        # deduplicate
        dedup_task = None
        if wf.args.dedup:
            dedup_task, args = wf.add_task(bismark_deduplicate(), name='BismarkDeduplicate-'+sample, depends=[align_task])
            args['bam'].value = [align_task.outputs['bam']]

        # methylation_task
        depend_tasks = [dedup_task if dedup_task is not None else align_task]
        methylation_task, args = wf.add_task(bismark_methylation_extractor(), name='MethylationExtractor-'+sample, depends=depend_tasks)
        args['bam'].value = [dedup_task.outputs['bam'] if dedup_task else align_task.outputs['bam']]

        # report
        if dedup_task:
            depends = [align_task, dedup_task, methylation_task]
            report_task, args = wf.add_task(bismark2report(), name='Bismark2report-' + sample, depends=depends)
            args['dedup_report'].value = dedup_task.outputs['report']
        else:
            depends = [align_task, methylation_task]
            report_task, args = wf.add_task(bismark2report(), name='Bismark2report-' + sample, depends=depends)
        args['alignment_report'].value = align_task.outputs['report']
        args['splitting_report'].value = methylation_task.outputs['report']
        args['mbias_report'].value = methylation_task.outputs['mbias']
        if align_task.cmd.args['nucleotide_coverage'].value:
            args['nucleotide_report'].value = align_task.outputs['nucleotide_stats']

    # summary
    depends = [k for k, v in wf.tasks.items() if v.name.startswith('BismarkAlignment')]
    summary_task, args = wf.add_task(bismark2summary(), name='bismark2summary', depends=depends)
    args['bams'].value = [wf.tasks[task_id].outputs['bam'] for task_id in depends]
    # finish
    wf.run()
    # to wdl
    # wf.to_wdl_tasks()


if __name__ == '__main__':
    pipeline()
