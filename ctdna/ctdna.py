import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, ToWdlTask
from utils.get_fastq_info import get_fastq_info
__author__ = 'gdq'


def creat_ref_dict():
    cmd = Command()
    cmd.meta.name = 'CreateSequenceDictionary'
    cmd.meta.desc = 'Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. '
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.args['copy_input_mode'] = Argument(prefix=f'cp -', default='L', range=['L', 'l', 's'], desc='indicate how to copy input fasta into work directory, "L": copy, "l": hard link, "s": softlink, do not use this if docker is used')
    cmd.args['ref_fasta'] = Argument(prefix='', type='infile', desc='reference fasta to index')
    cmd.args['_new_name'] = Argument(type='fix', value='ref_genome.fa && samtools faidx ref_genome.fa')
    cmd.args['_create_dict'] = Argument(type='fix', value='&& gatk CreateSequenceDictionary')
    cmd.args['_reference'] = Argument(prefix='-R ', type='fix', value='ref_genome.fa', desc='reference fasta')
    cmd.args['_ref_dict'] = Argument(prefix='-O ', type='fix', value='ref_genome.dict', desc="Output a dict file containing only the sequence dictionary. By default it will use the base name of the input reference with the .dict extension")
    cmd.args['alt_names'] = Argument(prefix='-ALT_NAMES ', type='infile', level='optional', desc="Optional file containing the alternative names for the contigs. Tools may use this information to consider different contig notations as identical (e.g: 'chr1' and '1'). The alternative names will be put into the appropriate @AN annotation for each contig. No header. First column is the original name, the second column is an alternative name. One contig may have more than one alternative name. Default value: null.")
    cmd.outputs['ref_genome'] = Output(value='ref_genome.fa')
    cmd.outputs['ref_dict'] = Output(value='{_ref_dict}')
    return cmd


def build_bwa_index():
    cmd = Command()
    cmd.meta.name = 'buildBwaIndex'
    cmd.meta.desc = 'bwa index and create sequence dictionary and fasta fai file'
    cmd.meta.version = '0.7.17'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 30 * 1024 ** 3
    cmd.args['copy_input_mode'] = Argument(prefix=f'cp -', default='s', range=['L', 'l', 's'], desc='indicate how to copy input fasta into work directory, "L": copy, "l": hard link, "s": softlink')
    cmd.args['ref_fasta'] = Argument(prefix='', type='infile', desc='reference fasta to index')
    cmd.args['_new_name'] = Argument(type='fix', value='ref_genome.fa && samtools faidx ref_genome.fa &&')
    cmd.args['_index_tool'] = Argument(type='fix', value='/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index')
    cmd.args['_ref_fasta'] = Argument(type='fix', value='ref_genome.fa', desc='reference fasta to be indexed')
    cmd.args['_create_dict'] = Argument(type='fix', value='&& gatk CreateSequenceDictionary')
    cmd.args['_reference'] = Argument(prefix='-R ', type='fix', value='ref_genome.fa', desc='reference fasta')
    cmd.args['ref_dict'] = Argument(prefix='-O ', default='ref_genome.dict', desc="Output a dict file containing only the sequence dictionary. By default it will use the base name of the input reference with the .dict extension")
    cmd.args['alt_names'] = Argument(prefix='-ALT_NAMES ', type='infile', level='optional', desc="Optional file containing the alternative names for the contigs. Tools may use this information to consider different contig notations as identical (e.g: 'chr1' and '1'). The alternative names will be put into the appropriate @AN annotation for each contig. No header. First column is the original name, the second column is an alternative name. One contig may have more than one alternative name. Default value: null.")
    cmd.outputs['ref_genome'] = Output(value='ref_genome.fa')
    cmd.outputs['ref_dict'] = Output(value='{ref_dict}')
    cmd.outputs['index_dir'] = Output(value='.')
    return cmd


def FastqToSam(sample):
    cmd = Command()
    cmd.meta.name = 'FastqToSam'
    cmd.meta.desc = 'convert fastq to sam'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'gatk FastqToSam'
    cmd.args['read1'] = Argument(prefix='-F1 ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-F2 ', level='optional', type='infile', desc='read2 fastq file')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.unmapped.bam', desc='output sam file')
    cmd.args['read_group_name'] = Argument(prefix='--READ_GROUP_NAME ', default=sample, desc='read group name')
    cmd.args['sample_name'] = Argument(prefix='--SAMPLE_NAME ', default=sample, desc='sample name')
    cmd.args['library_name'] = Argument(prefix='--LIBRARY_NAME ', default=sample, desc='library name')
    cmd.args['platform'] = Argument(prefix='--PLATFORM ', default='illumina', desc='sequencing platform name')
    cmd.args['tmpdir'] = Argument(prefix='--TMP_DIR ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def ExtractUmisFromBam():
    cmd = Command()
    cmd.meta.name = 'ExtractUmisFromBam'
    cmd.runtime.image = 'dceoy/fgbio:latest'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'fgbio ExtractUmisFromBam'
    cmd.args['input'] = Argument(prefix='-i ', type='infile', desc='Input BAM file')
    cmd.args['read-structure'] = Argument(prefix='-r ', array=True, default=['3M2S146T', '3M2S146T'], desc='The read structure, one per read in a template.')
    cmd.args['molecular-index-tags'] = Argument(prefix='-t ', array=True, default=['ZA', 'ZB'], desc='SAM tag(s) in which to store the molecular indices.')
    cmd.args['single-tag'] = Argument(prefix='-s ', default='RX', desc="Single tag into which to concatenate all molecular indices.")
    cmd.args['output'] = Argument(prefix='-o ', desc='Output BAM file')
    # 定义输出
    cmd.outputs['output'] = Output(value="{output}")
    return cmd


def MarkIlluminaAdapters():
    cmd = Command()
    cmd.meta.name = 'MarkIlluminaAdapters'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'gatk MarkIlluminaAdapters'
    cmd.args['input'] = Argument(prefix='--INPUT ', type='infile', desc='Input BAM file')
    cmd.args['metrics'] = Argument(prefix='--METRICS ', default='metrics.txt', desc='Histogram showing counts of bases_clipped in how many reads Required')
    cmd.args['output'] = Argument(prefix='--OUTPUT ', desc='Output BAM file')
    cmd.outputs['output'] = Output(value="{output}")
    cmd.outputs['metrics'] = Output(value="{metrics}")
    return cmd


def Bam2FastqBwaMem(sample):
    cmd = Command()
    cmd.meta.name = 'Bam2FastqBwaMem'
    cmd.meta.desc = 'bam to fastq and then mapping'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 15*1024**3
    cmd.runtime.cpu = 8
    cmd.args['_fix0'] = Argument(type='fix', value='gatk SamToFastq')
    cmd.args['input'] = Argument(prefix='-I ', type='infile', desc='input ubam file')
    cmd.args['CLIPPING_ATTRIBUTE'] = Argument(prefix='--CLIPPING_ATTRIBUTE ', level='optional', desc='The attribute that stores the position at which the SAM record should be clipped ')
    cmd.args['CLIPPING_ACTION'] = Argument(prefix='--CLIPPING_ACTION ', level='optional', desc="The action that should be taken with clipped reads: 'X' means the reads and qualities should be trimmed at the clipped position; 'N' means the bases should be changed to Ns in the clipped region; and any integer means that the base qualities should be set to that value in the clipped region. ")
    cmd.args['paired'] = Argument(prefix='--INTERLEAVE ', default='true', range=['true', 'false'], desc='if input is paired fastq, set it be true, else set it be false')
    cmd.args['_fix1'] = Argument(type='fix', value='--FASTQ /dev/stdout -NON_PF true')
    cmd.args['_fix2'] = Argument(type='fix', value='| /opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -M -Y -p -v 3', desc='this is bwa base command')
    cmd.args['k'] = Argument(prefix='-K ', default=10000000)
    cmd.args['t'] = Argument(prefix='-t ', default=cmd.runtime.cpu, desc='number of threads to use in mapping step')
    cmd.args['ref'] = Argument(prefix='', type='infile', desc='reference fasta file')
    cmd.args['_fix3'] = Argument(type='fix', value='/dev/stdin - ')
    cmd.args['_fix4'] = Argument(type='fix', value=' | samtools view -1 - ', desc='input data to samtools view')
    cmd.args['out'] = Argument(prefix='> ',  value=f'{sample}.unmerged.bam', desc='output bam file')
    cmd.outputs['out'] = Output(value="{out}")
    return cmd


def MergeBamAlignment(sample):
    cmd = Command()
    cmd.meta.name = 'MergeBamAlignment'
    cmd.meta.desc = 'merge bam alignment'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'gatk MergeBamAlignment'
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
    cmd.args['PROGRAM_GROUP_VERSION'] = Argument(prefix='--PROGRAM_GROUP_VERSION ', default='bwa-mem2-2.2.1_x64-linux')
    cmd.args['PROGRAM_GROUP_COMMAND_LINE'] = Argument(prefix='--PROGRAM_GROUP_COMMAND_LINE ', default='"bwa-mem2 mem -M -Y -p -v 3 -K 100000000 -t 4 ref.fa"')
    cmd.args['PROGRAM_GROUP_NAME'] = Argument(prefix='--PROGRAM_GROUP_NAME ', default='"bwamem"')
    cmd.args['tmpdir'] = Argument(prefix='--TMP_DIR ', default='.', desc='directory with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def MergeSamFiles():
    cmd = Command()
    cmd.meta.name = 'MergeSamFiles'
    cmd.meta.desc = 'Merges multiple SAM and/or BAM files into a single file.'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'gatk MergeSamFiles'
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', array=True, multi_times=True, desc='SAM or BAM input file')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', desc='SAM or BAM file to write merged result')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='true', range=['true', 'false'], desc='Whether to create a BAM index when writing a coordinate-sorted BAM file.')
    cmd.args['SORT_ORDER'] = Argument(prefix='--SORT_ORDER ', default='coordinate', range=['unsorted', 'queryname', 'coordinate', 'duplicate', 'unknown'], desc='Sort order of output file')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def GroupReadsByUmi(sample):
    cmd = Command()
    cmd.meta.name = 'GroupReadsByUmi'
    cmd.meta.source = 'https://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html'
    cmd.meta.desc = 'Groups reads together that appear to have come from the same original molecule. Reads are grouped by template, and then templates are sorted by the 5’ mapping positions of the reads from the template, used from earliest mapping position to latest. Reads that have the same end positions are then sub-grouped by UMI sequence.'
    cmd.runtime.image = 'dceoy/fgbio:latest'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'fgbio GroupReadsByUmi'
    cmd.args['input'] = Argument(prefix='-i ', type='infile', desc='the input BAM file.')
    cmd.args['output'] = Argument(prefix='-o ', default=f'{sample}.umi_grouped.bam', desc='The output BAM file.')
    cmd.args['family-size-histogram'] = Argument(prefix='-f ', default=f'{sample}.family.size.txt', desc='Optional output of tag family size counts.')
    cmd.args['raw-tag'] = Argument(prefix='-t ', default='RX', desc='The tag containing the raw UMI.')
    cmd.args['assign-tag'] = Argument(prefix='-T ', default='MI', desc='The output tag for UMI grouping.')
    cmd.args['min-map-q'] = Argument(prefix='-m ', default=1, desc='Minimum mapping quality for mapped reads.')
    cmd.args['strategy'] = Argument(prefix='-s ', default="paired", desc='The UMI assignment strategy. edit: reads are clustered into groups such that each read within a group has at least one other read in the group with <= edits differences and there are inter-group pairings with <= edits differences. Effective when there are small numbers of reads per UMI, but breaks down at very high coverage of UMIs. 3.adjacency: a version of the directed adjacency method described in umi_tools that allows for errors between UMIs but only when there is a count gradient.')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def CallDuplexConsensusReads():
    cmd = Command()
    cmd.meta.name = 'CallDuplexConsensusReads'
    cmd.meta.version = '2.1.0'
    cmd.meta.source = 'https://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html'
    cmd.meta.desc = 'Calls duplex consensus sequences from reads generated from the same double-stranded source molecule.'
    cmd.runtime.image = 'dceoy/fgbio:latest'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'fgbio CallDuplexConsensusReads'
    cmd.args['input'] = Argument(prefix='-i ', type='infile', desc='the input BAM file.')
    cmd.args['output'] = Argument(prefix='-o ', desc='The output BAM file.')
    cmd.args['error-rate-pre-umi'] = Argument(prefix='-1 ', default=45, desc='The Phred-scaled error rate for an error prior to the UMIs being integrated.')
    cmd.args['error-rate-post-umi'] = Argument(prefix='-2 ', default=30, desc='The Phred-scaled error rate for an error post the UMIs have been integrated.')
    cmd.args['min-input-base-quality'] = Argument(prefix='-m ', default=30, desc='Ignore bases in raw reads that have Q below this value.')
    cmd.args['threads'] = Argument(prefix='--threads ', default=cmd.runtime.cpu, desc='The number of threads to use while consensus calling')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def FilterConsensusReads():
    cmd = Command()
    cmd.meta.name = 'FilterConsensusReads'
    cmd.meta.source = 'https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html'
    cmd.meta.desc = 'Filters consensus reads generated by CallMolecularConsensusReads or CallDuplexConsensusReads.'
    cmd.runtime.image = 'dceoy/fgbio:latest'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'fgbio FilterConsensusReads'
    cmd.args['input'] = Argument(prefix='-i ', type='infile', desc='the input BAM file.')
    cmd.args['output'] = Argument(prefix='-o ', desc='The output BAM file.')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='Reference fasta file.')
    cmd.args['min-reads'] = Argument(prefix='-M ', default='2 1 1', desc='The minimum number of reads supporting a consensus base/read.')
    cmd.args['PhredScore'] = Argument(prefix='-N ', default=20, desc="Mask (make 'N') consensus bases with quality less than this threshold")
    cmd.args['require-single-strand-agreement'] = Argument(prefix='-s ', default='true', desc='Mask (make N) consensus bases where the AB and BA consensus reads disagree (for duplex-sequencing only).')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def ClipBam():
    cmd = Command()
    cmd.meta.name = 'ClipBam'
    cmd.meta.source = 'https://fulcrumgenomics.github.io/fgbio/tools/latest/ClipBam.html'
    cmd.meta.desc = 'Clips reads from the same template. Ensures that at least N bases are clipped from any end of the read (i.e. R1 5’ end, R1 3’ end, R2 5’ end, and R2 3’ end). Optionally clips reads from the same template to eliminate overlap between the reads.'
    cmd.runtime.image = 'dceoy/fgbio:latest'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'fgbio ClipBam'
    cmd.args['input'] = Argument(prefix='-i ', type='infile', desc='the input BAM file.')
    cmd.args['output'] = Argument(prefix='-o ', desc='The output BAM file.')
    cmd.args['ref'] = Argument(prefix='-r ', desc='Reference fasta file.')
    cmd.args['clipping-mode'] = Argument(prefix='-c ', default='Hard', range=['Hard', 'Soft', 'SoftWithMask'], desc='The type of clipping to perform.')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def FilterBam():
    cmd = Command()
    cmd.meta.name = 'FilterBam'
    cmd.meta.source = 'https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterBam.html'
    cmd.meta.desc = """
    Filters reads out of a BAM file. Removes reads that may not be useful in downstream processing or visualization. 
    By default will remove unmapped reads, reads with MAPQ=0, reads marked as secondary alignments, reads marked as duplicates, 
    and if a set of Intervals are provided, reads that do not overlap any of the intervals.
    """
    cmd.runtime.image = 'dceoy/fgbio:latest'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'fgbio FilterBam'
    cmd.args['input'] = Argument(prefix='-i ', type='infile', desc='The input BAM file.')
    cmd.args['output'] = Argument(prefix='-o ', desc='The output BAM file.')
    cmd.args['remove-duplicates'] = Argument(prefix='-D ', default='true', range=['true', 'false'], desc='If true remove all reads that are marked as duplicates')
    cmd.args['remove-unmapped-reads'] = Argument(prefix='-U ', default='true', desc='Remove all unmapped reads.')
    cmd.args['min-map-q'] = Argument(prefix='-M ', default=1, desc='Remove all mapped reads with MAPQ lower than this number.')
    cmd.args['remove-single-end-mappings'] = Argument(prefix='-P ', default='false', desc='Removes non-PE reads and any read whose mate pair is unmapped.')
    cmd.args['remove-secondary-alignments'] = Argument(prefix='-S ', default='true', desc='Remove all reads marked as secondary alignments.')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def SortAndIndexBam():
    cmd = Command()
    cmd.meta.name = 'SortAndIndexBam'
    cmd.meta.source = 'https://www.htslib.org/doc/samtools-sort.html'
    cmd.meta.desc = "sort and index bam using samtools"
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 12 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'samtools sort'
    cmd.args['output'] = Argument(prefix='-o ', desc='The output BAM file.')
    cmd.args['threads'] = Argument(prefix='--threads ', default=cmd.runtime.cpu, desc='Number of additional threads to use')
    cmd.args['output-fmt'] = Argument(prefix='--output-fmt ', default='BAM', desc='output format')
    cmd.args['input'] = Argument(prefix='', type='infile', desc='input bam file')
    cmd.args['_fix'] = Argument(type='fix', value=f'&& samtools index -@ {cmd.runtime.cpu} *.bam')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def Bamdst():
    cmd = Command()
    cmd.meta.name = 'SortAndIndexBam'
    cmd.meta.source = 'https://github.com/shiquan/bamdst'
    cmd.meta.desc = 'Bamdst is a lightweight tool to stat the depth coverage of target regions of bam file(s).'
    cmd.runtime.image = 'biocontainers/bamdst:1.0.9_cv1'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'bamdst'
    cmd.args['bed'] = Argument(prefix='-p ', type='infile', desc='probe bed file')
    cmd.args['outdir'] = Argument(prefix='-o ', default='.', desc='output directory')
    cmd.args['flank'] = Argument(prefix='-f ', default=200, desc='calculate the coverage of flank region')
    cmd.args['cutoffdepth'] = Argument(prefix='--cutoffdepth ', default=500, desc='the specified coverage')
    cmd.args['input'] = Argument(prefix='', type='infile', desc='input bam file')
    cmd.outputs['outdir'] = Output(value='{outdir}')
    return cmd


def VardictSingle():
    cmd = Command()
    cmd.meta.name = 'VardictSingle'
    cmd.meta.source = 'https://github.com/AstraZeneca-NGS/VarDictJava'
    cmd.meta.version = 'VarDict_v1.8.2'
    cmd.meta.desc = "VarDictJava is a variant discovery program written in Java and Perl."
    cmd.runtime.image = 'docker.io/truwl/vardict-java:latest'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'vardict-java'
    cmd.args['sample'] = Argument(prefix='-N ', desc='sample name')
    cmd.args['bam'] = Argument(prefix='-b ', type='infile', desc='The indexed BAM file')
    cmd.args['genome'] = Argument(prefix='-G ', type='infile', desc='The reference fasta. Should be indexed (.fai).')
    cmd.args['threads'] = Argument(prefix='-th ', default=8, desc='Threads count.')
    cmd.args['min-freq'] = Argument(prefix='-f ', default="0.00001", desc='The threshold for allele frequency')
    cmd.args['chromosome'] = Argument(prefix='-c ', default=1, desc='The column of chromosome')
    cmd.args['region_start'] = Argument(prefix='-S ', default=2, desc='The column of region start')
    cmd.args['region_end'] = Argument(prefix='-E ', default=3, desc='The column of region end')
    cmd.args['gene'] = Argument(prefix='-g ', default=4, desc='The column of gene name')
    cmd.args['fisher'] = Argument(prefix='--fisher', type='bool', default=True, desc='Fisher exact test.')
    cmd.args['realignment'] = Argument(prefix='-k', type='bool', default=True, desc='Indicate whether to perform local realignment')
    cmd.args['mfreq'] = Argument(prefix='-mfreq ', default=0.25, desc="The variant frequency threshold to determine variant as good in case of monomer MSI. Default: 0.25")
    cmd.args['nmfreq'] = Argument(prefix='-nmfreq ', default=0.1, desc='The variant frequency threshold to determine variant as good in case of non-monomer MSI')
    cmd.args['read_position_filter'] = Argument(prefix='-P ', default=3, desc="If the mean variants position is less that specified, it's considered false")
    cmd.args['min-reads'] = Argument(prefix='-r ', default=3, desc='The minimum # of variant reads')
    cmd.args['nosv'] = Argument(prefix='--nosv', type='bool', default=True,)
    cmd.args['UN'] = Argument(prefix='-UN', type='bool', default=True, desc='Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using first read only')
    cmd.args['bed'] = Argument(prefix='', type='infile', desc='region or bed file')
    cmd.args['_fix'] = Argument(type='fix', value='| var2vcf_valid.pl -E', desc='pipe to another script')
    cmd.args['min-freq2'] = Argument(prefix='-f ', default="0.00001", desc='The threshold for allele frequency')
    cmd.args['output'] = Argument(prefix='> ', desc='output vcf name')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def add_vcf_contig():
    cmd = Command()
    cmd.meta.name = 'AddVcfContig'
    cmd.meta.desc = 'Add contig info for vardict output vcf'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 2 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'python'
    cmd.args['script'] = Argument(prefix='', type='infile', value=f'{script_path}/utils/add_vcf_contig.py', desc='script path')
    cmd.args['vcf'] = Argument(prefix='-vcf ', type='infile', desc='path to vcf')
    cmd.args['ref_dict'] = Argument(prefix='-ref_dict ', level='optional', type='infile', desc='path to reference dict file')
    cmd.args['header_txt'] = Argument(prefix='-header_txt ', level='optional', type='infile', desc='path to file which contains header lines')
    cmd.args['out'] = Argument(prefix='-out ', desc='output file name')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def bcftools_norm():
    cmd = Command()
    cmd.meta.name = 'VcfLeftNorm'
    cmd.meta.desc = "Left-align and normalize indels; check if REF alleles match the reference; split multiallelic sites into multiple rows; recover multiallelics from multiple rows"
    cmd.runtime.image = 'dceoy/bcftools:latest'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = "bcftools norm"
    cmd.args['fasta-ref'] = Argument(prefix='-f ', type='infile', desc='reference fasta file')
    cmd.args['multiallelics'] = Argument(prefix='-m ', level='optional', desc='Split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]' )
    cmd.args['out'] = Argument(prefix='-o ', desc='Write output to a file [standard output]')
    cmd.args['output-type'] = Argument(prefix='--output-type ', default='v', desc="'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]")
    cmd.args['threads'] = Argument(prefix='--threads ', default=4, desc="Use multithreading with <int> worker threads")
    cmd.args['vcf'] = Argument(type='infile', desc='input vcf file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def vep(sample):
    cmd = Command()
    cmd.meta.name = 'VEP'
    cmd.runtime.image = 'ensemblorg/ensembl-vep:release_109.3'
    # 由于权限问题，更改docker_cmd_prefix
    cmd.runtime.docker_cmd_prefix = cmd.runtime.docker_cmd_prefix2
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'vep'
    cmd.args['input_file'] = Argument(prefix='-i ', type='infile', desc='input file')
    cmd.args['fasta'] = Argument(prefix='--fasta ', type='infile', desc="Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run VEP with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache).")
    cmd.args['output_file'] = Argument(prefix='-o ', default=f'{sample}.vep.vcf.gz', desc='output file')
    cmd.args['output_format'] = Argument(prefix='--', range={'vcf', 'json', 'tab'}, default='vcf', desc="If we choose to write output in VCF format. Consequences are added in the INFO field of the VCF file, using the key 'CSQ'. Data fields are encoded separated by '|'; the order of fields is written in the VCF header. Output fields in the 'CSQ' INFO field can be selected by using --fields.")
    cmd.args['compress_output'] = Argument(prefix='--compress_output ', default='bgzip', desc="Writes output compressed using either gzip or bgzip")
    cmd.args['force_overwrite'] = Argument(prefix="--force_overwrite ", type='bool', default=True, desc="Force overwriting of output file")
    cmd.args['fork'] = Argument(prefix='--fork ', type='int', default=cmd.runtime.cpu, desc='Use forking(multi-cpu/threads) to improve script runtime')
    cmd.args['species'] = Argument(prefix='--species ', default='homo_sapiens', desc='Species for your data. This can be the latin name e.g. homo_sapiens or any Ensembl alias e.g. mouse.')
    cmd.args['assembly_version'] = Argument(prefix='--assembly ', default='GRCh37', desc='Select the assembly version to use if more than one available.')
    # docker有可能报读写权限问题,可以考虑chmod -R a+rwx $HOME/vep_data，实际上需要用其他用户身份进入容器才有权限，比如用-u指定用户身份
    cmd.args['dir_cache'] = Argument(prefix='--dir_cache ', type='indir', desc='Specify the cache directory to use')
    cmd.args['dir_plugins'] = Argument(prefix='--dir_plugins ', type='indir', level='optional', desc='Specify the plugin directory to use')
    cmd.args['stats_file'] = Argument(prefix='--stats_file ', default=f'{sample}.vep.summary.html', desc='Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>.')
    cmd.args['cache'] = Argument(prefix='--cache ', type='bool', default=True, desc='Enables use of cache')
    cmd.args['offline'] = Argument(prefix='--offline ', type='bool', default=True, desc='Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation')
    cmd.args['merged'] = Argument(prefix='--merged ', type='bool', default=False, desc='Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.')
    cmd.args['refseq'] = Argument(prefix='--refseq ', type='bool', default=False, desc='Specify this option if you have installed the RefSeq cache in order for VEP to pick up the alternate cache directory. This cache contains transcript objects corresponding to RefSeq transcripts. Consequence output will be given relative to these transcripts in place of the default Ensembl transcripts')
    cmd.args['plugins'] = Argument(prefix='--plugin ', level='optional', multi_times=True, default=['Frameshift', 'Wildtype'], desc='Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times')
    cmd.args['variant_class'] = Argument(prefix='--variant_class ', type='bool', default=True, desc='Output the Sequence Ontology variant class.')
    cmd.args['sift'] = Argument(prefix='--sift ', default='b', range={'p', 's', 'b'}, desc="Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both.")
    cmd.args['polyphen'] = Argument(prefix='--polyphen ', default='b', range={'p', 's', 'b'}, desc="Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both.")
    cmd.args['nearest'] = Argument(prefix='--nearest ', default='transcript', range={'transcript', 'gene', 'symbol'}, desc='Retrieve the transcript or gene with the nearest protein-coding transcription start site (TSS) to each input variant. Use transcript to retrieve the transcript stable ID, gene to retrieve the gene stable ID, or symbol to retrieve the gene symbol. Note that the nearest TSS may not belong to a transcript that overlaps the input variant, and more than one may be reported in the case where two are equidistant from the input coordinates.')
    cmd.args['gene_phenotype'] = Argument(prefix='--gene_phenotype ', type='bool', default=True, desc='Indicates if the overlapped gene is associated with a phenotype, disease or trait.')
    cmd.args['regulatory'] = Argument(prefix='--regulatory ', type='bool', default=True, desc="Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature.")
    cmd.args['phased'] = Argument(prefix='--phased ', type='bool', default=True, desc="Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data.")
    cmd.args['numbers'] = Argument(prefix='--numbers ', type='bool', default=True, desc="Adds affected exon and intron numbering to to output. Format is Number/Total")
    cmd.args['hgvs'] = Argument(prefix='--hgvs ',  type='bool', default=True, desc="Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate.")
    cmd.args['transcript_version'] = Argument(prefix='--transcript_version ', type='bool', default=True, desc="Add version numbers to Ensembl transcript identifiers")
    cmd.args['symbol'] = Argument(prefix='--symbol ', type='bool', default=True, desc="Adds the gene symbol (e.g. HGNC) (where available) to the output.")
    cmd.args['tsl'] = Argument(prefix='--tsl ', type='bool', default=True, desc="Adds the transcript support level for this transcript to the output.")
    cmd.args['canonical'] = Argument(prefix='--canonical ', type='bool', default=True, desc="Adds a flag indicating if the transcript is the canonical transcript for the gene")
    cmd.args['biotype'] = Argument(prefix='--biotype ', type='bool', default=True, desc="Adds the biotype of the transcript or regulatory feature.")
    cmd.args['max_af'] = Argument(prefix='--max_af ', type='bool', default=True, desc="Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD")
    cmd.args['af_1kg'] = Argument(prefix='--af_1kg ', type='bool', default=True, desc="Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output.")
    cmd.args['af_gnomad'] = Argument(prefix='--af_gnomad ', type='bool', default=True, desc="Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included")
    cmd.args['af_esp'] = Argument(prefix='--af_esp ', type='bool', default=False, desc="Include allele frequency from NHLBI-ESP populations.")
    cmd.args['coding_only'] = Argument(prefix='--af_esp ', type='bool', default=False, desc="Only return consequences that fall in the coding regions of transcripts. Not used by default")
    cmd.args['pick'] = Argument(prefix='--pick', type='bool', default=False, desc="Pick one line or block of consequence data per variant, including transcript-specific columns. This is the best method to use if you are interested only in one consequence per variant")
    cmd.args['flag_pick'] = Argument(prefix='--flag_pick ', type='bool', default=True, desc="As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others.")
    cmd.args['filter_common'] = Argument(prefix='--filter_common ', type='bool', default=False, desc="Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters.")
    cmd.args['other_args'] = Argument(default='', desc='specify other arguments that you want to append to the command')
    cmd.args['_create_index'] = Argument(value='&& tabix *vcf.gz', type='fix')
    cmd.outputs['out_vcf'] = Output(value='{output_file}')
    cmd.outputs['out_vcf_idx'] = Output(value='{output_file}.tbi')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'ctDNA'
    wf.meta.source = "https://doi.org/10.1038/s41587-021-00857-z"
    wf.meta.desc = \
    """
    参考IDT的分析流程
    1. FastqToSam
    2. ExtractUmisFromBam
    3. MarkIlluminaAdapters
    4. SamToFastq and bwa-mem
    5. MergeBamAlignment
    6. GroupReadsByUmi
    7. CallDuplexConsensusReads
    8. bam2fastq and bwa-mem
    9. FilterConsensusReads
    10. CollectHsMetrics
    11. ClipBam(选择跳过，vardict call变异时可以通过参数指定不重复计数overlapped区域）
    12. filterbam:Supplementary aligned reads and not primary aligned reads were removed using samtools v1.5（可以用fgbio自带得filterbam）
    13. vardict
    14. (无)Low frequency mutations that were called in the Sample B replicates were removed in all other samples, and mutations flagged as p8 in the filter column were removed
    """
    wf.meta.version = "1.0"

    # 定义流程输入参数
    wf.init_argparser()
    # fastq 输入参数
    wf.add_argument('-fastq_info', nargs='+', required=True,
                    help='A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]')
    wf.add_argument('-r1_name', default='(.*)_S\d+_L\d+_R1_\d+.fastq.gz',
                    help="python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'")
    wf.add_argument('-r2_name', default='(.*)_S\d+_L\d+_R2_\d+.fastq.gz',
                    help="python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'")
    wf.add_argument('-exclude_samples', default=tuple(), nargs='+', help='samples to exclude from analysis')
    wf.add_argument('-bed', help="bed file for target region")
    # 参考数据库参数
    wf.add_argument('-ref', default='/home/hxbio04/dbs/hg19/hs37d5.fa', help='reference fasta file')
    wf.add_argument('-vep_cache', default='/home/hxbio04/dbs/vep', help='VEP cache directory')
    wf.add_argument('-vep_plugin', required=False, help='VEP plugin directory')

    # 收集参数
    wf.parse_args()
    top_vars = dict(
        ref=TopVar(value=os.path.abspath(wf.args.ref), type='infile'),
        bed=TopVar(value=os.path.abspath(wf.args.bed), type='infile'),
        vep_cache=TopVar(value=os.path.abspath(wf.args.vep_cache), type='indir'),
        vep_plugin=TopVar(value=os.path.abspath(wf.args.vep_plugin) if wf.args.vep_plugin else None, type='indir'),
    )
    wf.add_topvars(top_vars)

    # 提取fastq信息
    fastq_info = get_fastq_info(fastq_info=wf.args.fastq_info, r1_name=wf.args.r1_name, r2_name=wf.args.r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found !')

    # 建bwa索引, 也会同时建fai和dict文件
    make_index = False
    index_task = None
    if not os.path.exists(wf.topvars['ref'].value + '.0123'):
        make_index = True
        print('We will build bwa-mem2 index for provide reference')
        index_task, args = wf.add_task(build_bwa_index(), name='buildIndex')
        if not wf.args.docker:
            args['copy_input_mode'].value = 's'
        else:
            args['copy_input_mode'].value = 'L'
        args['ref_fasta'].value = wf.topvars['ref']

    # create dict or fai for reference if necessary
    dict_file = wf.topvars['ref'].value.rsplit('.', 1)[0] + '.dict'
    fai_file = wf.topvars['ref'].value + '.fai'
    wf.topvars['ref_dict'] = TopVar(value=dict_file, type='infile')
    if not os.path.exists(dict_file) or not os.path.exists(fai_file):
        print(f'dict or fai file for {wf.args.ref} does not exist, we will try to create one rightly')
        creat_dict_task = creat_ref_dict()
        creat_dict_task.args['ref_fasta'].value = wf.topvars['ref'].value
        tmp_wkdir = os.path.join(wf.args.outdir, "PreparedInputs")
        if wf.args.run:
            creat_dict_task.run_now(wkdir=tmp_wkdir, docker=wf.args.docker)
        wf.topvars['ref'].value = creat_dict_task.outputs['ref_genome'].value
        wf.topvars['ref_dict'].value = creat_dict_task.outputs['ref_dict'].value

    # 开始处理
    for sample, reads in fastq_info.items():
        # 跳过不需要分析的样本
        if sample in wf.args.exclude_samples:
            print(f'Skip {sample} for it is in excluded sample list')
            continue
        # fastq预处理和比对 考虑一个样本存在多对fastq的处理
        if sample in wf.args.exclude_samples:
            continue
        if len(reads) == 2:
            r1s, r2s = reads
        else:
            r1s = reads[0]
            r2s = [None]*len(r1s)

        merge_bam_tasks = []
        for ind, (r1, r2) in enumerate(zip(r1s, r2s)):
            uniq_tag = f'{sample}-{ind}' if len(r1s) > 1 else sample
            # fastq2sam
            fastq2sam_task, args = wf.add_task(FastqToSam(sample), tag=uniq_tag)
            args['read1'].value = r1
            args['read2'].value = r2
            args['platform'].value = 'Illumina'
            args['out'].value = f'{uniq_tag}.unmapped.bam'

            # ExtractUmisFromBam
            get_umi_task, args = wf.add_task(ExtractUmisFromBam(), tag=uniq_tag, depends=[fastq2sam_task])
            args['input'].value = fastq2sam_task.outputs['out']
            args['read-structure'].value = ['2S3M2S144T', '2S3M2S144T']
            args['output'].value = f'{uniq_tag}.umi.ubam'

            # 去除接头MarkIlluminaAdapters
            markadapter_task, args = wf.add_task(MarkIlluminaAdapters(), tag=uniq_tag, depends=[get_umi_task])
            args['input'].value = get_umi_task.outputs['output']
            args['metrics'].value = uniq_tag + '.markadapters.metircs.txt'
            args['output'].value = uniq_tag + '.markadapters.ubam'

            # bwa alignment
            bwa_task, args = wf.add_task(Bam2FastqBwaMem(uniq_tag), tag=uniq_tag, depends=[index_task, markadapter_task])
            args['input'].value = markadapter_task.outputs['output']
            args['CLIPPING_ATTRIBUTE'].value = 'XT'
            args['CLIPPING_ACTION'].value = 'X'
            args['ref'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']

            # merge mapped bam and unmapped bam
            merge_bam_task, args = wf.add_task(MergeBamAlignment(uniq_tag), tag=uniq_tag, depends=[bwa_task])
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
            args['ALIGNED_BAM'].value = bwa_task.outputs['out']
            args['UNMAPPED_BAM'].value = markadapter_task.outputs['output']
            merge_bam_tasks.append(merge_bam_task)

        if len(r1s) > 1:
            # 合并一个样本的多个fastq的比对结果
            merge_sam_task, args = wf.add_task(MergeSamFiles(), tag=sample, depends=merge_bam_tasks)
            args['INPUT'].value = [x.outputs['out'] for x in merge_bam_tasks]
            # GroupReadsByUmi 不强求bam是否是sorted
            args['SORT_ORDER'].value = 'unsorted'
            args['OUTPUT'].value = sample + '.merged.bam'
        else:
            merge_sam_task = merge_bam_task

        group_umi_task, args = wf.add_task(GroupReadsByUmi(sample), tag=sample, depends=[merge_sam_task])
        args['input'].value = merge_sam_task.outputs['out']

        consensus_task, args = wf.add_task(CallDuplexConsensusReads(), tag=sample, depends=[group_umi_task])
        args['input'].value = group_umi_task.outputs['output']
        args['output'].value = sample + '.consensus.bam'

        filter_consensus_task, args = wf.add_task(FilterConsensusReads(), tag=sample, depends=[consensus_task])
        args['input'].value = consensus_task.outputs['output']
        args['ref'].value = wf.topvars['ref']
        args['output'].value = sample + '.filtered_consensus.bam'

        map_task, args = wf.add_task(Bam2FastqBwaMem(sample), tag=sample+'-remap', depends=[filter_consensus_task])
        args['input'].value = filter_consensus_task.outputs['output']
        args['ref'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']

        filter_bam_task, args = wf.add_task(FilterBam(), tag=sample, depends=[map_task])
        args['input'].value = map_task.outputs['out']
        args['output'].value = sample + '.filtered.bam'

        bamdst_task, args = wf.add_task(Bamdst(), tag=sample, depends=[filter_bam_task])
        args['input'].value = filter_bam_task.outputs['output']
        args['bed'].value = wf.topvars['bed']

        sort_bam_task, args = wf.add_task(SortAndIndexBam(), tag=sample, depends=[filter_bam_task])
        args['input'].value = filter_bam_task.outputs['output']
        args['output'].value = sample + '.sorted.bam'

        vardict_task, args = wf.add_task(VardictSingle(), tag=sample, depends=[sort_bam_task])
        args['sample'].value = sample
        args['bam'].value = sort_bam_task.outputs['output']
        args['genome'].value = wf.topvars['ref']
        args['bed'].value = wf.topvars['bed']
        args['output'].value = sample + '.raw.vcf'

        add_contig_task, args = wf.add_task(add_vcf_contig(), tag=sample, depends=[vardict_task])
        args['vcf'].value = vardict_task.outputs['output']
        args['ref_dict'].value = wf.topvars['ref_dict']
        args['out'].value = sample + '.raw.vcf'

        vcf_norm_task, args = wf.add_task(bcftools_norm(), tag=sample, depends=[add_contig_task])
        args['fasta-ref'].value = wf.topvars['ref']
        args['multiallelics'].value = '+both'
        args['vcf'].value = add_contig_task.outputs['out']
        args['out'].value = sample + '.normed.raw.vcf'

        vep_task, args = wf.add_task(vep(sample), tag=sample, depends=[vcf_norm_task])
        args['input_file'].value = vcf_norm_task.outputs['out']
        args['fasta'].value = wf.topvars['ref']
        args['refseq'].value = True
        args['dir_cache'].value = top_vars['vep_cache']
        args['dir_plugins'].value = top_vars['vep_plugin']
        vep_task.outputs['out_vcf'].report = True
        vep_task.outputs['out_vcf_idx'].report = True
    # end
    wf.run()


if __name__ == '__main__':
    pipeline()

