import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, ToWdlTask
from utils.get_fastq_info import get_fastq_info
__author__ = 'gdq'


def fastq_to_sam(sample):
    cmd = Command()
    cmd.meta.name = 'fastq2sam'
    cmd.meta.desc = 'convert fastq to sam'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.tool = 'java -jar picard.jar FastqToSam'
    cmd.args['read1'] = Argument(prefix='F1=', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='F2=', level='optional', type='infile', desc='read2 fastq file')
    cmd.args['out'] = Argument(prefix='O=', default=f'{sample}.unmapped.bam', desc='output sam file')
    cmd.args['read_group_name'] = Argument(prefix='--READ_GROUP_NAME ', default=sample, desc='read group name')
    cmd.args['sample_name'] = Argument(prefix='--SAMPLE_NAME ', default=sample, desc='sample name')
    cmd.args['library_name'] = Argument(prefix='--LIBRARY_NAME ', default=sample, desc='library name')
    cmd.args['platform'] = Argument(prefix='--PLATFORM ', default='illumina', desc='sequencing platform name')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def ubam_to_fastq_and_bwa_mem(sample):
    cmd = Command()
    cmd.meta.name = 'uBam2FastqBwaMem'
    cmd.meta.desc = 'ubam to fastq and then mapping'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5*1024**3
    cmd.runtime.cpu = 16
    cmd.args['_fix0'] = Argument(type='fix', value='java -jar picard.jar SamToFastq')
    cmd.args['ubam'] = Argument(prefix='I=', type='infile', desc='input ubam file')
    cmd.args['paired'] = Argument(prefix='INTERLEAVE=', default='true', range=['true', 'false'], desc='if input is paired fastq, set it be true, else set it be false')
    cmd.args['_fix1'] = Argument(type='fix', value='FASTQ=/dev/stdout NON_PF=true')
    cmd.args['_fix2'] = Argument(type='fix', value='| bwa mem -M -Y', desc='this is bwa base command')
    cmd.args['k'] = Argument(prefix='-K ', default=10000000)
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in mapping step')
    cmd.args['ref'] = Argument(prefix='', type='infile', desc='reference fasta file')
    cmd.args['_fix3'] = Argument(type='fix', value=f'/dev/stdin - 2> >(tee {sample}.bwa.stderr.log >&2) ')
    cmd.args['_fix4'] = Argument(type='fix', value=' | samtools view -1 - ', desc='input data to samtools view')
    cmd.args['out'] = Argument(prefix='> ',  value=f'{sample}.unmerged.bam', desc='output bam file')
    cmd.outputs['out'] = Output(value="{out}")
    cmd.outputs['out_log'] = Output(value=f'{sample}.bwa.stderr.log')
    return cmd


def merge_bam_alignment(sample):
    cmd = Command()
    cmd.meta.name = 'MergeBamAlignment'
    cmd.meta.desc = 'merge bam alignment'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.tool = 'java -jar picard.jar MergeBamAlignment'
    cmd.args['VALIDATION_STRINGENCY'] = Argument(prefix='--VALIDATION_STRINGENCY ', default='SILENT')
    cmd.args['EXPECTED_ORIENTATIONS'] = Argument(prefix='--EXPECTED_ORIENTATIONS ', level='optional')
    cmd.args['ATTRIBUTES_TO_RETAIN'] = Argument(prefix='--ATTRIBUTES_TO_RETAIN ', default='X0')
    cmd.args['ALIGNED_BAM'] = Argument(prefix='--ALIGNED_BAM ', type='infile', desc='SAM or BAM file')
    cmd.args['UNMAPPED_BAM'] = Argument(prefix='--UNMAPPED_BAM ', type='infile', desc='unmapped bam file')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.merged.unsorted.bam', desc='output bam file')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='--REFERENCE_SEQUENCE ', type='infile', desc='reference fasta file')
    cmd.args['SORT_ORDER'] = Argument(prefix='--SORT_ORDER ', default='"unsorted"')
    cmd.args['CLIP_ADAPTERS'] = Argument(prefix='--CLIP_ADAPTERS ', default='false', desc='Whether to clip adapters where identified.')
    cmd.args['MAX_RECORDS_IN_RAM'] = Argument(prefix='--MAX_RECORDS_IN_RAM ', default=2000000, desc='When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.')
    cmd.args['MAX_INSERTIONS_OR_DELETIONS'] = Argument(prefix='--MAX_INSERTIONS_OR_DELETIONS ', default=-1, desc='The maximum number of insertions or deletions permitted for an alignment to be included.')
    cmd.args['PRIMARY_ALIGNMENT_STRATEGY'] = Argument(prefix='--PRIMARY_ALIGNMENT_STRATEGY ', default='MostDistant')
    cmd.args['UNMAPPED_READ_STRATEGY'] = Argument(prefix='--UNMAPPED_READ_STRATEGY ', default='COPY_TO_TAG', desc='How to deal with alignment information in reads that are being unmapped (e.g. due to cross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS = true. Note that the DO_NOT_CHANGE strategy will actually reset the cigar and set the mapping quality on unmapped reads since otherwisethe result will be an invalid record. To force no change use the DO_NOT_CHANGE_INVALID strategy.')
    cmd.args['ALIGNER_PROPER_PAIR_FLAGS'] = Argument(prefix='--ALIGNER_PROPER_PAIR_FLAGS ', default='true', desc="Use the aligner's idea of what a proper pair is rather than computing in this program.")
    cmd.args['UNMAP_CONTAMINANT_READS'] = Argument(prefix='--UNMAP_CONTAMINANT_READS ', default='true', desc='Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial sample),and unmap + label those reads accordingly.')
    cmd.args['PROGRAM_RECORD_ID'] = Argument(prefix='--PROGRAM_RECORD_ID ', default='"bwamem"')
    cmd.args['PROGRAM_GROUP_VERSION'] = Argument(prefix='--PROGRAM_GROUP_VERSION ', default='0.7.17')
    cmd.args['PROGRAM_GROUP_COMMAND_LINE'] = Argument(prefix='--PROGRAM_GROUP_COMMAND_LINE ', default='')
    cmd.args['PROGRAM_GROUP_NAME'] = Argument(prefix='--PROGRAM_GROUP_NAME ', default='"bwamem"')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def picard_mark_duplicates(sample):
    cmd = Command()
    cmd.meta.name = 'MarkDuplicates'
    cmd.meta.desc = 'merge bam alignment'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.tool = 'java -jar picard.jar MarkDuplicates'
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', array=True, multi_times=True, desc='input bam file list')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.unsorted.dup_marked.bam', desc='output bam file')
    cmd.args['METRICS_FILE'] = Argument(prefix='--METRICS_FILE ', default=f'{sample}.dup_metrics.txt')
    cmd.args['VALIDATION_STRINGENCY'] = Argument(prefix='--VALIDATION_STRINGENCY ', default='SILENT')
    cmd.args['OPTICAL_DUPLICATE_PIXEL_DISTANCE'] = Argument(prefix='--OPTICAL_DUPLICATE_PIXEL_DISTANCE ', default=2500, desc='The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and models, users should experiment to find what works best.')
    cmd.args['ASSUME_SORT_ORDER'] = Argument(prefix='--ASSUME_SORT_ORDER ', default='"queryname')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def sort_and_fix_tags(sample):
    cmd = Command()
    cmd.meta.name = 'SortAndFixTags'
    cmd.meta.desc = 'sort bam and fix tags'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['_sort_sam'] = Argument(type='fix', value='java -jar picard.jar SortSam')
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', array=True, multi_times=True, desc='input bam file list')
    cmd.args['_OUTPUT'] = Argument(prefix='--OUTPUT ', default='/dev/stdout', desc='output bam file')
    cmd.args['SORT_ORDER'] = Argument(prefix='--SORT_ORDER ', default='"coordinate"')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='false')
    cmd.args['_fix_tags'] = Argument(type='fix', value='java -jar picard.jar SetNmMdAndUqTags --CREATE_INDEX true --INPUT /dev/stdin')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.sorted.dup_marked.bam', desc='output bam file')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='--REFERENCE_SEQUENCE ', type='infile', desc='reference fasta file')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    cmd.outputs['out_idx'] = Output(value='{OUTPUT}.bai')
    return cmd


def CreateSequenceGroupingTSV():
    cmd = Command()
    cmd.meta.name = 'SortAndFixTags'
    cmd.meta.desc = 'Generate sets of intervals for scatter-gathering over chromosomes'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5 * 1024 ** 3
    pass


def base_recalibrator(sample):
    cmd = Command()
    cmd.meta.name = 'BaseRecalibrator'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5*1024**3
    cmd.runtime.tool = 'gatk BaseRecalibrator'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['INPUT'] = Argument(prefix='-I ', type='infile', desc='input bam file')
    cmd.args['use-original-qualities'] = Argument(prefix='--use-original-qualities', type='bool', default=True)
    cmd.args['OUTPUT'] = Argument(prefix='-O ', default=f'{sample}.recal_table', desc='The output recalibration table file to create')
    cmd.args['known-sites'] = Argument(prefix='--known-sites ', type='infile', array=True, multi_times=True, desc='One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.')
    cmd.args['intervals'] = Argument(prefix='--intervals ', type='infile', array=True, multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def gather_bqsr_reports(sample):
    cmd = Command()
    cmd.meta.name = 'GatherBQSRReports'
    cmd.meta.desc = 'Gathers scattered BQSR recalibration reports into a single file'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 3*1024**3
    cmd.runtime.tool = 'gatk GatherBQSRReports'
    cmd.args['INPUT'] = Argument(prefix='-I ', type='infile', array=True, multi_times=True, desc='input bqsr file list')
    cmd.args['OUTPUT'] = Argument(prefix='-O ', default=f'{sample}.recal_table', desc='The output recalibration table file to create')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def apply_bqsr(sample):
    cmd = Command()
    cmd.meta.name = 'ApplyBQSR'
    cmd.meta.desc = 'Apply Base Quality Score Recalibration (BQSR) model'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk ApplyBQSR'
    cmd.args['INPUT'] = Argument(prefix='-I ', type='infile', desc='input bam file')
    cmd.args['OUTPUT'] = Argument(prefix='-O ', default=f'{sample}.recal_table', desc='The output recalibration table file to create')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['intervals'] = Argument(prefix='--intervals ', type='infile', array=True, multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['bqsr'] = Argument(prefix='-bqsr ', type='infile', desc='Input recalibration table for BQSR')
    cmd.args['static-quantized-quals'] = Argument(prefix='--static-quantized-quals ', type='int', array=True, multi_times=True, default=[10, 20, 30])
    cmd.args['add-output-sam-program-record'] = Argument(prefix='--add-output-sam-program-record', type='bool', default=True)
    cmd.args['use-original-qualities'] = Argument(prefix='--use-original-qualities', type='bool', default=True)
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def gather_bams(sample):
    cmd = Command()
    cmd.meta.name = 'GatherBamFiles'
    cmd.meta.desc = 'Concatenate efficiently BAM files that resulted from a scattered parallel analysis.'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk GatherBamFiles'
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', array=True, multi_times=True, desc='input bam file')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.gathered.bam', desc='output bam file')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='true')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    cmd.outputs['out_idx'] = Output(value='{OUTPUT}.bai')
    return cmd



