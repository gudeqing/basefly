import os

import pandas as pd

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


def BedToIntervalList():
    cmd = Command()
    cmd.meta.name = 'BedToIntervalList'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk BedToIntervalList'
    cmd.args['bed'] = Argument(prefix='-I ', type='infile', desc='input bed file')
    cmd.args['ref_dict'] = Argument(prefix='-SD ', type='infile', desc='input sequence dictionary file')
    cmd.args['out'] = Argument(prefix='-O ',  desc='output interval file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def SplitIntervals(scatter_number):
    cmd = Command()
    cmd.meta.name = 'SplitIntervals'
    cmd.meta.desc = 'This tool takes in intervals via the standard arguments of IntervalArgumentCollection and splits them into interval files for scattering. The resulting files contain equal number of bases.'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk SplitIntervals'
    cmd.args['ref'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['intervals'] = Argument(prefix='-L ', level='optional', type='infile', multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['scatter'] = Argument(prefix='-scatter ', type='int', value=scatter_number, desc='number of output interval files to split into')
    cmd.args['mode'] = Argument(prefix='-mode ', default="INTERVAL_SUBDIVISION", range=['INTERVAL_SUBDIVISION', 'BALANCING_WITHOUT_INTERVAL_SUBDIVISION', 'BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW'], desc='How to divide intervals.')
    cmd.args['interval-merging-rule'] = Argument(prefix='--interval-merging-rule ', default='OVERLAPPING_ONLY')
    cmd.args['outdir'] = Argument(prefix='-O ', default='.', desc='The directory into which to write the scattered interval sub-directories.')
    for i in range(scatter_number):
        # 如果mode=BALANCING_WITHOUT_INTERVAL_SUBDIVISION，有可能生成的interval数量少于指定的scatter—number
        cmd.outputs[f'out{i}'] = Output(value='{outdir}/'+f'{i:04d}-scattered.interval_list', type='outfile')
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


def SamToFastq():
    cmd = Command()
    cmd.meta.name = 'SamToFastq'
    cmd.meta.desc = 'Use samtools to convert sam to Fastq. This tool enables user to copy tag to fastq header line'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'samtools fastq'
    cmd.args['read1'] = Argument(prefix='-1 ', desc='write paired reads flagged READ1 to FILE')
    cmd.args['read2'] = Argument(prefix='-2 ', desc='write paired reads flagged READ2 to FILE')
    cmd.args['tags'] = Argument(prefix='-T ', level='optional', array=True, delimiter=',', desc='copy arbitrary tags to the FASTQ header line')
    cmd.args['threads'] = Argument(prefix='--threads ', default=cmd.runtime.cpu, desc='Number of additional threads to use')
    cmd.args['bam'] = Argument(prefix='', desc='input bam file')
    cmd.outputs['read1'] = Output(value='{read1}')
    cmd.outputs['read2'] = Output(value='{read2}')
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
    cmd.args['k'] = Argument(prefix='-K ', default=100000000, desc='process INT input bases in each batch regardless of nThreads (for reproducibility)')
    cmd.args['t'] = Argument(prefix='-t ', default=cmd.runtime.cpu, desc='number of threads to use in mapping step')
    cmd.args['ref'] = Argument(prefix='', type='infile', desc='reference fasta file')
    cmd.args['_fix3'] = Argument(type='fix', value='/dev/stdin - ')
    cmd.args['_fix4'] = Argument(type='fix', value=f' | samtools view --threads {cmd.runtime.cpu} -1 - ', desc='input data to samtools view')
    cmd.args['out'] = Argument(prefix='> ',  value=f'{sample}.unmerged.bam', desc='output bam file')
    cmd.outputs['out'] = Output(value="{out}")
    return cmd


def bwa_mem(sample, platform):
    cmd = Command()
    cmd.meta.name = 'BwaMem2'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.tool = '/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -M -Y -v 3'
    cmd.runtime.memory = 15 * 1024 ** 3
    cmd.runtime.cpu = 8
    cmd.args['include_read_header'] = Argument(prefix='-C', type='bool', default=False, desc='Append FASTA/FASTQ comment to SAM output')
    cmd.args['readgroup'] = Argument(prefix='-R ', desc='read group info', value=f'"@RG\\tID:{sample}\\tSM:{sample}\\tPL:{platform}"')
    cmd.args['k'] = Argument(prefix='-K ', default=100000000, desc='process INT input bases in each batch regardless of nThreads (for reproducibility)')
    cmd.args['t'] = Argument(prefix='-t ', default=cmd.runtime.cpu, desc='number of threads to use in computation')
    cmd.args['ref'] = Argument(prefix='', type='infile', desc='reference fasta file')
    cmd.args['read1'] = Argument(prefix='', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='', level='optional', type='infile', desc='read2 fastq file')
    cmd.args['_fix'] = Argument(type='fix', value=f'| samtools view -O BAM --threads {cmd.runtime.cpu} - ')
    cmd.args['out'] = Argument(prefix='> ', value=f'{sample}.unmerged.bam', desc='output bam file')
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
    cmd.args['strategy'] = Argument(prefix='-s ', default="adjacency", desc='The UMI assignment strategy. edit: reads are clustered into groups such that each read within a group has at least one other read in the group with <= edits differences and there are inter-group pairings with <= edits differences. Effective when there are small numbers of reads per UMI, but breaks down at very high coverage of UMIs. 3.adjacency: a version of the directed adjacency method described in umi_tools that allows for errors between UMIs but only when there is a count gradient.')
    cmd.args['family-size-histogram'] = Argument(prefix='-f ', default=f'{sample}.family.size.txt', desc='Optional output of tag family size counts.')
    cmd.args['raw-tag'] = Argument(prefix='-t ', default='RX', desc='The tag containing the raw UMI.')
    cmd.args['assign-tag'] = Argument(prefix='-T ', default='MI', desc='The output tag for UMI grouping.')
    cmd.args['min-map-q'] = Argument(prefix='-m ', default=1, desc='Minimum mapping quality for mapped reads.')
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
    cmd.args['min-reads'] = Argument(prefix='--min-reads ', default=1, desc='The minimum number of reads to produce a consensus base')
    cmd.args['error-rate-pre-umi'] = Argument(prefix='-1 ', default=45, desc='The Phred-scaled error rate for an error prior to the UMIs being integrated.')
    cmd.args['error-rate-post-umi'] = Argument(prefix='-2 ', default=30, desc='The Phred-scaled error rate for an error post the UMIs have been integrated.')
    cmd.args['min-input-base-quality'] = Argument(prefix='-m ', default=30, desc='Ignore bases in raw reads that have Q below this value.')
    cmd.args['threads'] = Argument(prefix='--threads ', default=cmd.runtime.cpu, desc='The number of threads to use while consensus calling')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def CallMolecularConsensusReads():
    cmd = Command()
    cmd.meta.name = 'CallMolecularConsensusReads'
    cmd.meta.version = '2.1.0'
    cmd.meta.source = 'https://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html'
    cmd.meta.desc = 'Calls consensus sequences from reads with the same unique molecular tag.'
    cmd.runtime.image = 'dceoy/fgbio:latest'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'fgbio CallMolecularConsensusReads'
    cmd.args['input'] = Argument(prefix='-i ', type='infile', desc='the input BAM file.')
    cmd.args['output'] = Argument(prefix='-o ', desc='The output BAM file.')
    cmd.args['min-reads'] = Argument(prefix='--min-reads ', default=1, desc='The minimum number of reads to produce a consensus base')
    cmd.args['error-rate-pre-umi'] = Argument(prefix='-1 ', default=45, desc='The Phred-scaled error rate for an error prior to the UMIs being integrated.')
    cmd.args['error-rate-post-umi'] = Argument(prefix='-2 ', default=30, desc='The Phred-scaled error rate for an error post the UMIs have been integrated.')
    cmd.args['min-input-base-quality'] = Argument(prefix='-m ', default=20, desc='Ignore bases in raw reads that have Q below this value.')
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
    cmd.args['min-reads'] = Argument(prefix='-M ', default='1 1 1', desc='The minimum number of reads supporting a consensus base/read.')
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


def gencore():
    cmd = Command()
    cmd.meta.name = 'Gencore'
    cmd.meta.source = ''
    cmd.meta.desc = ''
    cmd.meta.version = '0.17.2'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 12 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'gencore'
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input sorted bam/sam file.')
    cmd.args['out'] = Argument(prefix='-o ', desc='output bam/sam file')
    cmd.args['ref'] = Argument(prefix='--ref ', type='infile', desc='reference fasta file name (should be an uncompressed .fa/.fasta file)')
    cmd.args['bed'] = Argument(prefix='--bed ', type='infile', desc='bed file to specify the capturing region')
    cmd.args['duplex_only'] = Argument(prefix='--duplex_only', type='bool', default=False, desc='only output duplex consensus sequences, which means single stranded consensus sequences will be discarded')
    cmd.args['no_duplex'] = Argument(prefix='--no_duplex', type='bool', default=False, desc="don't merge single stranded consensus sequences to duplex consensus sequences.")
    cmd.args['umi_prefix'] = Argument(prefix='--umi_prefix ', level='optional', desc='the prefix for UMI, if it has. None by default. Check the README for the defails of UMI formats. such as: "NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_GAGCATAC", prefix = "UMI", umi = "GAGCATAC"')
    cmd.args['supporting_reads'] = Argument(prefix='--supporting_reads ', default=1, desc='only output consensus reads/pairs that merged by >= <supporting_reads> reads/pairs. The valud should be 1~10')
    cmd.args['ratio_threshold'] = Argument(prefix='--ratio_threshold ', default=0.8, desc='if the ratio of the major base in a cluster is less than <ratio_threshold>, it will be further compared to the reference. The valud should be 0.5~1.0')
    cmd.args['score_threshold'] = Argument(prefix='--score_threshold ', default=6, desc='if the score of the major base in a cluster is less than <score_threshold>, it will be further compared to the reference.')
    cmd.args['umi_diff_threshold'] = Argument(prefix='--umi_diff_threshold ',  level='optional', default=1, desc='if two reads with identical mapping position have UMI difference <= <umi_diff_threshold>, then they will be merged to generate a consensus read.')
    cmd.args['duplex_diff_threshold'] = Argument(prefix='duplex_diff_threshold ', level='optional', default=2, desc='if the forward consensus and reverse consensus sequences have <= <duplex_diff_threshold> mismatches, then they will be merged to generatea duplex consensus sequence, otherwise will be discarded')
    cmd.args['coverage_sampling'] = Argument(prefix='--coverage_sampling ', default=10000, desc='the sampling rate for genome scale coverage statistics. Default 10000 means 1/10000.')
    cmd.args['json'] = Argument(prefix='--json ', default='gencore.json', desc='the json format report file name')
    cmd.args['html'] = Argument(prefix='--html ', default='gencore.html', desc='the html format report file name')
    cmd.outputs['json'] = Output(value='{json}')
    cmd.outputs['html'] = Output(value='{html}')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def Bamdst():
    cmd = Command()
    cmd.meta.name = 'Bamdst'
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
    cmd.args['mfreq'] = Argument(prefix='-mfreq ', default=0.25, desc="The variant frequency threshold to determine variant as good in case of monomer MSI. Default: 0.25")
    cmd.args['nmfreq'] = Argument(prefix='-nmfreq ', default=0.1, desc='The variant frequency threshold to determine variant as good in case of non-monomer MSI')
    cmd.args['read_position_filter'] = Argument(prefix='-P ', default=3, desc="If the mean variants position is less that specified, it's considered false")
    cmd.args['min-reads'] = Argument(prefix='-r ', default=2, desc='The minimum # of variant reads')
    cmd.args['nosv'] = Argument(prefix='--nosv', type='bool', default=True,)
    cmd.args['UN'] = Argument(prefix='-UN', type='bool', default=True, desc='Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using first read only')
    cmd.args['bed'] = Argument(prefix='', type='infile', desc='region or bed file')
    cmd.args['_fix'] = Argument(type='fix', value='| var2vcf_valid.pl -A -E -p 2 -q 15 -d 3 -v 1 -f 0.00001 ', desc='pipe to another script')
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


def Mutect2(prefix):
    cmd = Command()
    cmd.meta.name = 'Mutect2'
    cmd.meta.desc = 'Call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations.'
    cmd.meta.source = 'https://github.com/broadinstitute/gatk'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 5*1024**3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'gatk Mutect2'
    cmd.args['ref'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['tumor_bam'] = Argument(prefix='-I ', type='infile', desc='tumor bam')
    cmd.args['normal_bam'] = Argument(prefix='-I ', type='infile', level='optional', desc='normal bam')
    cmd.args['tumor_name'] = Argument(prefix='-tumor ', desc='tumor sample name')
    cmd.args['normal_name'] = Argument(prefix='-normal ', level='optional', desc='normal sample name')
    cmd.args['germline-resource'] = Argument(prefix='--germline-resource ', type='infile', level='optional', desc='optional database of known germline variants (and its index) (see http://gnomad.broadinstitute.org/downloads)')
    cmd.args['pon'] = Argument(prefix='-pon ', type='infile', level='optional', desc='')
    cmd.args['intervals'] = Argument(prefix='-L ', type='infile', level='optional', multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['string_intervals'] = Argument(prefix='-L ', level='optional', desc='interval string such as "chrM"')
    cmd.args['alleles'] = Argument(prefix='--alleles ', type='infile', level='optional', desc='The set of alleles to force-call regardless of evidence')
    cmd.args['initial-tumor-lod'] = Argument(prefix='--initial-tumor-lod ', default=3.0, desc='Log 10 odds threshold to consider pileup active.')
    cmd.args['normal-lod'] = Argument(prefix='--normal-lod ', default=2.2, desc='Log 10 odds threshold for calling normal variant non-germline. ')
    cmd.args['tumor-lod-to-emit'] = Argument(prefix='--tumor-lod-to-emit ', default=3.0, desc='Log 10 odds threshold to emit variant to VCF')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.vcf.gz', desc='output vcf')
    cmd.args['bam-output'] = Argument(prefix='--bam-output ', level='optional', desc='output bam file')
    cmd.args['f1r2-tar-gz'] = Argument(prefix='--f1r2-tar-gz ', type='infile', default=f'{prefix}.f1r2.tar.gz', desc='If specified, collect F1R2 counts and output files into this tar.gz file')
    cmd.args['mitochondria'] = Argument(prefix='--mitochondira', type='bool', desc='if to turn on mitochondria mode. Specifically, the mode sets --initial-tumor-lod to 0, --tumor-lod-to-emit to 0, --af-of-alleles-not-in-resource to 4e-3, and the advanced parameter --pruning-lod-thres')
    cmd.args['tmpdir'] = Argument(prefix='--tmp-dir ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{out}')
    cmd.outputs['f1r2'] = Output(value='{f1r2-tar-gz}')
    cmd.outputs['stats'] = Output(value='{out}.stats')
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


def stat_context_seq_error():
    cmd = Command()
    cmd.meta.name = 'StatSeqError'
    cmd.meta.desc = '假设大部分低频突变是测序错误，基于bam文件估计测序错误率'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 2 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'python'
    cmd.args['script'] = Argument(prefix='', type='infile', value=f'{script_path}/utils/stat_3bases_error.py', desc='script path')
    cmd.args['bam'] = Argument(prefix='-bam ', type='infile', level='optional', desc='path to bam file which will be used to estimate background noise')
    cmd.args['bed'] = Argument(prefix='-bed ', type='infile', level='optional', desc='path to target region file which will be used to estimate background noise')
    cmd.args['genome'] = Argument(prefix='-genome ', type='infile', desc='path to indexed genome fasta')
    cmd.args['exclude_from'] = Argument(prefix='-exclude_from ', type='infile', level='optional', desc='bed or vcf file containing known variant in input bam, these variants will be excluded during background noise estimating')
    cmd.args['center_size'] = Argument(prefix='-center_size ', default=(1, 1), array=True, desc='extending size around ref base during background noise estimating')
    cmd.args['out_prefix'] = Argument(prefix='-out_prefix ', desc='output file prefix')
    cmd.outputs['stat_per_site'] = Output(value='{out_prefix}.each_site.txt', report=True)
    cmd.outputs['context_error_rate'] = Output(value='{out_prefix}.centered11_site.json', report=True)
    return cmd


def VcfFilter():
    cmd = Command()
    cmd.meta.name = 'VcfFilter'
    cmd.meta.desc = 'filtering vcf, 输出符合室间质评的格式结果'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.memory = 2 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'python'
    cmd.args['script'] = Argument(prefix='', type='infile', value=f'{script_path}/utils/vcf_filter.py', desc='script path')
    cmd.args['vcf'] = Argument(prefix='-vcf ', type='infile', desc='path to vcf annotated with vep')
    cmd.args['genome'] = Argument(prefix='-genome ', type='infile', desc='path to indexed genome fasta')
    cmd.args['ref_dict'] = Argument(prefix='-ref_dict ', type='infile', level='optional', desc='path to genome dict file which will be used to add contig header in vcf')
    cmd.args['bam'] = Argument(prefix='-bam ', type='infile', level='optional', desc='path to bam file which will be used to estimate background noise')
    cmd.args['bed'] = Argument(prefix='-bed ', type='infile', level='optional', desc='path to target region file which will be used to estimate background noise')
    cmd.args['exclude_from'] = Argument(prefix='-exclude_from ', type='infile', level='optional', desc='bed or vcf file containing known variant in input bam, these variants will be excluded during background noise estimating')
    cmd.args['center_size'] = Argument(prefix='-center_size ', default=(1, 1), array=True,  desc='extending size around ref base during background noise estimating')
    cmd.args['tumor_name'] = Argument(prefix='-tumor_name ', level='optional', desc='tumor sample name in vcf')
    cmd.args['normal_vcf'] = Argument(prefix='-normal_vcf ', type='infile', level='optional', desc='normal sample vcf file')
    cmd.args['error_rate_file'] = Argument(prefix='-error_rate_file ', type='infile', level='optional', desc='Estimated background noise file, if not provided, bam file will be used')
    cmd.args['min_error_rate'] = Argument(prefix='-min_error_rate ', default=1e-6, desc='global minimum error rate, if error rate cannot be aquired in other ways, this value will be used')
    cmd.args['alpha'] = Argument(prefix='-alpha ', default=0.01, desc='cutoff of pvalue from background noise model')
    cmd.args['out_prefix'] = Argument(prefix='-out_prefix ', desc='output file prefix')
    cmd.outputs['final_vcf'] = Output(value='{out_prefix}.final.vcf', report=True)
    cmd.outputs['final_txt'] = Output(value='{out_prefix}.final.txt', report=True)
    cmd.outputs['final_xls'] = Output(value='{out_prefix}.final.xlsx', report=True)
    cmd.outputs['discarded_vcf'] = Output(value='{out_prefix}.discarded.vcf', report=True)
    cmd.outputs['log'] = Output(value='{out_prefix}.filtering.log', report=True)
    return cmd


def MergeVcfs(sample):
    cmd = Command()
    cmd.meta.name = 'MergeVcfs'
    cmd.meta.desc = 'Combines multiple variant files into a single variant file.'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360056969852-MergeVcfs-Picard-'
    cmd.runtime.tool = 'gatk MergeVcfs'
    cmd.args['inputs'] = Argument(prefix='-I ', type='infile', multi_times=True, desc='input vcf list')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.vcf.gz', desc='The merged VCF or BCF file. File format is determined by file extension.')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def LearnReadOrientationModel(sample):
    cmd = Command()
    cmd.meta.name = 'LearnReadOrientationModel'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.meta.desc = 'Learn the prior probability of read orientation artifact from the output of CollectF1R2Counts of Mutect2'
    cmd.runtime.tool = 'gatk LearnReadOrientationModel'
    cmd.args['inputs'] = Argument(prefix='-I ', type='infile', multi_times=True, desc='One or more .tar.gz containing outputs of CollectF1R2Counts')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.artifact-priors.tar.gz', desc='tar.gz of artifact prior tables')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def FilterMutectCalls(sample):
    cmd = Command()
    cmd.meta.name = 'FilterMutectCalls'
    cmd.meta.desc = 'FilterMutectCalls applies filters to the raw output of Mutect2'
    cmd.runtime.tool = 'gatk FilterMutectCalls'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='A VCF file containing variants')
    cmd.args['ref'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.filtered.vcf.gz', desc='output vcf file')
    cmd.args['contamination-table'] = Argument(prefix='--contamination-table ', level='optional', type='infile')
    cmd.args['tumor-segmentation'] = Argument(prefix='--tumor-segmentation ', level='optional', type='infile')
    cmd.args['ob-priors'] = Argument(prefix='--ob-priors ', type='infile', level='optional')
    cmd.args['stats'] = Argument(prefix='-stats ', type='infile', level='optional')
    cmd.args['filtering-stats'] = Argument(prefix='--filtering-stats ', default=f'{sample}.filtering.stats', desc='output filtering stat file')
    cmd.outputs['out'] = Output(value='{out}')
    cmd.outputs['filtering-stats'] = Output(value='{filtering-stats}')
    return cmd


def FilterAlignmentArtifacts(sample):
    cmd = Command()
    cmd.meta.name = 'FilterAlignmentArtifacts'
    cmd.meta.desc = 'Alignment artifacts can occur whenever there is sufficient sequence similarity between two or more regions in the genome to confuse the alignment algorithm. This can occur when the aligner for whatever reason overestimate how uniquely a read maps, thereby assigning it too high of a mapping quality. It can also occur through no fault of the aligner due to gaps in the reference, which can also hide the true position to which a read should map. By using a good alignment algorithm (the GATK wrapper of BWA-MEM), giving it sensitive settings (which may have been impractically slow for the original bam alignment) and mapping to the best available reference we can avoid these pitfalls. The last point is especially important: one can (and should) use a BWA-MEM index image corresponding to the best reference, regardless of the reference to which the bam was aligned.'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4418051467035-FilterAlignmentArtifacts-EXPERIMENTAL-'
    cmd.runtime.tool = 'gatk FilterAlignmentArtifacts'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='A VCF file containing variants')
    cmd.args['ref'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='input bam file')
    cmd.args['bwa-mem-index-image'] = Argument(prefix='--bwa-mem-index-image ', type='infile', desc='BWA-mem index image')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.align_artifacts_filtered.vcf.gz', desc='output vcf file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def MergeMutectStats(sample):
    cmd = Command()
    cmd.meta.name = 'MergeMutectStats'
    cmd.runtime.image = 'gudeqing/gatk-bwamem2-gencore:1.0'
    cmd.runtime.tool = 'gatk MergeMutectStats'
    cmd.args['stats'] = Argument(prefix='-stats ', type='infile', multi_times=True)
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.vcf.stats', desc='output merged stat files')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'ctDNA'
    wf.meta.source = ""
    wf.meta.desc = \
    """
    部分参考IDT的分析流程https://doi.org/10.1038/s41587-021-00857-z
    1. FastqToSam
    2. ExtractUmisFromBam
    3. MarkIlluminaAdapters
    4. SamToFastq and bwa-mem
    5. MergeBamAlignment
    6. GroupReadsByUmi
    7. CallDuplexConsensusReads (如果不是双端UMI，则不建议使用该方法，而用CallMolecularConsensusReads）
    8. FilterConsensusReads
    9. bam2fastq and bwa-mem (使用samtools fastq命令，带入consensus相关tag，然后bwa-mem比对时，加上-C参数，相关tag带到bam文件中，方便后续的使用）
    10. CollectHsMetrics或者bamdst
    11. ClipBam(选择跳过，vardict call变异时可以通过参数指定不重复计数overlapped区域）
    12. filterbam:Supplementary aligned reads and not primary aligned reads were removed using samtools v1.5（可以用fgbio自带得filterbam）
    13. vardict
    14. (无)Low frequency mutations that were called in the Sample B replicates were removed in all other samples, and mutations flagged as p8 in the filter column were removed
    
    尝试评估哪个假阳性更少:
    (1) CallDuplexConsensusReads vs CallMolecularConsensusReads
    预期CallDuplexConsensusReads效果大于CallMolecularConsensusReads
    （2）mutect2 vs vardict:
    预期mutect2效果好于vardict，个人认为最大的理由是因为vardict近几年都停止更新，而mutect2还在积极的维护
    (3) 检查每一个突变支持的read的UMI状态
    (4) 是否需要find_complex_variant
    * 由于分析结果中，很多的UMI family size==1, 可以考虑进一步使用gencore进行简单意义上的去重，进一步减少假阳性
    * 过滤时，对支持变异的read的family size进行检查，如果family size 小于3，考虑如何用一个更合适的背景测序错误率：
        可以考虑利用fgbio提供的cE信息，使用最后的3base context的策略统计背景噪音过滤
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
    wf.add_argument('-pair', required=False, help='Optional. pair information file, no header, tab separated, first column is tumor while second one is normal. Normal sample named "None" means tumor-only.')
    wf.add_argument('-umi', required=True, help='A string describes the read structural. Such as “1S3M3S144T,1S3M3S144T” denotes UMIs locate at 2-4bp of read1 and read2')
    wf.add_argument('-scatter', default=10, help='scatter number used for interval splitting of mutect2 variant calling steps')

    # 参考数据库参数
    wf.add_argument('-ref', default='/home/hxbio04/dbs/hg19/hs37d5.fa', help='reference fasta file')
    wf.add_argument('-vep_cache', default='/home/hxbio04/dbs/vep', help='VEP cache directory')
    wf.add_argument('-vep_plugin', required=False, help='VEP plugin directory')
    wf.add_argument('-germline_vcf', default='/home/hxbio04/dbs/hg19/af-only-gnomad.raw.sites.b37.vcf.gz', help='for Mutect2 input, will be used for germline variant filtering')
    wf.add_argument('-bwaMemIndexImage', required=False, help='bwa-mem-index-mage for artifact alignment filtering. you may created it with tool BwaMemIndexImageCreator with only fasta as input')


    # 收集参数
    wf.parse_args()
    top_vars = dict(
        ref=TopVar(value=os.path.abspath(wf.args.ref), type='infile'),
        bed=TopVar(value=os.path.abspath(wf.args.bed), type='infile'),
        vep_cache=TopVar(value=os.path.abspath(wf.args.vep_cache), type='indir'),
        vep_plugin=TopVar(value=os.path.abspath(wf.args.vep_plugin) if wf.args.vep_plugin else None, type='indir'),
        pair_info=TopVar(value=os.path.abspath(wf.args.pair) if wf.args.pair else None),
        germline_vcf=TopVar(value=os.path.abspath(wf.args.germline_vcf), type='infile'),
        bwaMemIndexImage=TopVar(value=os.path.abspath(wf.args.bwaMemIndexImage) if wf.args.bwaMemIndexImage else None, type='infile'),
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

    # interval creation for bed file
    if not os.path.exists(wf.topvars['bed'].value + '.interval_list'):
        bed2intervals_task, args = wf.add_task(BedToIntervalList())
        args['bed'].value = wf.topvars['bed']
        args['ref_dict'].value = wf.topvars['ref_dict']
        args['out'].value = os.path.basename(wf.topvars['bed'].value) + '.interval_list'
    else:
        wf.topvars['bed_interval'] = TopVar(value=wf.topvars['bed'].value + '.interval_list')
        bed2intervals_task = None

    # split intervals
    split_task, args = wf.add_task(SplitIntervals(scatter_number=int(wf.args.scatter)), tag='ForCaller')
    args['ref'].value = wf.topvars['ref']
    args['intervals'].value = [wf.topvars['bed_interval'] if not bed2intervals_task else bed2intervals_task.outputs['out']]
    args['outdir'].value = '.'
    interval_files = [split_task.outputs[f'out{i}'] for i in range(int(wf.args.scatter))]

    # 开始处理
    bam_task_dict = dict()
    vep_task_dict = dict()
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
            # args['read-structure'].value = ['1S3M3S144T', '1S3M3S144T']
            args['read-structure'].value = wf.args.umi.split(',')
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
            # 为了在mergeBamAlignment时和markIlluminaAdapters的输入保持一致，建议下面的设置为‘N'
            args['CLIPPING_ACTION'].value = 'N'
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
        # args['strategy'].value = 'paired'

        # consensus_task, args = wf.add_task(CallDuplexConsensusReads(), tag=sample, depends=[group_umi_task])
        consensus_task, args = wf.add_task(CallMolecularConsensusReads(), tag=sample, depends=[group_umi_task])
        args['input'].value = group_umi_task.outputs['output']
        args['output'].value = sample + '.consensus.bam'

        filter_consensus_task, args = wf.add_task(FilterConsensusReads(), tag=sample, depends=[consensus_task])
        args['input'].value = consensus_task.outputs['output']
        args['ref'].value = wf.topvars['ref']
        args['output'].value = sample + '.filtered_consensus.bam'

        sam2fastq_task, args = wf.add_task(SamToFastq(), tag=sample, depends=[filter_consensus_task])
        args['bam'].value = filter_consensus_task.outputs['output']
        args['read1'].value = sample + '.consensus.R1.fq.gz'
        args['read2'].value = sample + '.consensus.R2.fq.gz'
        args['tags'].value = ['cD', 'cE', 'cM']

        map_task, args = wf.add_task(bwa_mem(sample, 'Illumina'), tag=sample, depends=[sam2fastq_task])
        args['read1'].value = sam2fastq_task.outputs['read1']
        args['read2'].value = sam2fastq_task.outputs['read2']
        args['ref'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
        args['include_read_header'].value = True

        filter_bam_task, args = wf.add_task(FilterBam(), tag=sample, depends=[map_task])
        args['input'].value = map_task.outputs['out']
        args['output'].value = sample + '.filtered.bam'

        sort_bam_task, args = wf.add_task(SortAndIndexBam(), tag=sample, depends=[filter_bam_task])
        args['input'].value = filter_bam_task.outputs['output']
        args['output'].value = sample + '.sorted.bam'
        bam_task_dict[sample] = sort_bam_task

        gencore_task, args = wf.add_task(gencore(), tag=sample, depends=[sort_bam_task])
        args['bam'].value = sort_bam_task.outputs['output']
        args['bed'].value = wf.topvars['bed']
        args['ref'].value = wf.topvars['ref']
        args['out'].value = sample + '.gencore.unsorted.bam'
        args['json'].value = sample + '.gencore.json'
        args['html'].value = sample + '.gencore.html'

        sort_bam_task, args = wf.add_task(SortAndIndexBam(), tag=sample+'-re', depends=[gencore_task])
        args['input'].value = gencore_task.outputs['out']
        args['output'].value = sample + '.sorted.bam'
        bam_task_dict[sample] = sort_bam_task

        bamdst_task, args = wf.add_task(Bamdst(), tag=sample, depends=[sort_bam_task])
        args['input'].value = sort_bam_task.outputs['output']
        args['bed'].value = wf.topvars['bed']

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
        vep_task_dict[sample] = vep_task

    # filtering vcf
    if wf.topvars['pair_info']:
        pairs = [x.strip().split('\t')[:2] for x in open(wf.topvars['pair_info'].value)]
    else:
        pairs = zip(bam_task_dict.keys(), ['None']*len(bam_task_dict))

    vardict_filter_task_ids = []
    mutect2_filter_task_ids = []
    for tumor, normal in pairs:
        tumor_bam_task = bam_task_dict[tumor]
        tumor_vardict_vep_task = vep_task_dict[tumor]
        if normal != 'None':
            normal_vep_task = vep_task_dict[normal]
            normal_vcf = normal_vep_task.outputs['out_vcf']
            normal_bam_task = bam_task_dict[normal]
        else:
            normal_vcf = None
            normal_bam_task = None

        # ----mutect2---------------------------------------------------
        mutect_tasks = []
        for ind, interval_file in enumerate(interval_files):
            mutect_task, args = wf.add_task(Mutect2(f'{tumor}-{ind}'), tag=f'{tumor}-{ind}', parent_wkdir='Mutect2')
            mutect_tasks.append(mutect_task)
            args['ref'].value = wf.topvars['ref']
            args['tumor_bam'].value = tumor_bam_task.outputs['output']
            args['tumor_name'].value = tumor
            if normal_bam_task:
                args['normal_bam'].value = normal_bam_task.outputs['output']
                args['normal_name'].value = normal
            args['bam-output'].value = f'{tumor}-{ind}' + '.haplotypes.bam'
            args['germline-resource'].value = wf.topvars['germline_vcf']
            args['intervals'].value = [interval_file]

        # LearnReadOrientationModel
        lrom_task, args = wf.add_task(LearnReadOrientationModel(tumor), tag=tumor, depends=mutect_tasks)
        args['inputs'].value = [x.outputs['f1r2'] for x in mutect_tasks]

        # merge vcf
        merge_vcf_task, args = wf.add_task(MergeVcfs(tumor), tag=tumor, depends=mutect_tasks)
        args['inputs'].value = [x.outputs['out'] for x in mutect_tasks]

        # normalize vcf
        norm_vcf_task, args = wf.add_task(bcftools_norm(), tag=tumor, depends=[merge_vcf_task])
        args['fasta-ref'].value = wf.topvars['ref']
        args['multiallelics'].value = '+both'
        args['vcf'].value = merge_vcf_task.outputs['out']
        args['out'].value = tumor + '.mutect2.normed.raw.vcf'

        # merge stats
        merge_stat_task, args = wf.add_task(MergeMutectStats(tumor), tag=tumor, depends=mutect_tasks)
        args['stats'].value = [x.outputs['stats'] for x in mutect_tasks]
        merge_stat_task.outputs['out'].report = True

        # filter
        filter_task, args = wf.add_task(FilterMutectCalls(tumor), tag=tumor, depends=[norm_vcf_task, merge_stat_task, lrom_task])
        args['vcf'].value = norm_vcf_task.outputs['out']
        args['ref'].value = wf.topvars['ref']
        args['ob-priors'].value = lrom_task.outputs['out']
        args['stats'].value = merge_stat_task.outputs['out']
        filter_task.outputs['out'].report = True

        # filter alignment artifact
        filter_align_task = None
        if wf.topvars['bwaMemIndexImage'].value is not None:
            filter_align_task, args = wf.add_task(FilterAlignmentArtifacts(tumor), tag=tumor, depends=[filter_task])
            args['vcf'].value = filter_task.outputs['out']
            args['ref'].value = wf.topvars['ref']
            args['bam'].value = tumor_bam_task.outputs['output']
            args['bwa-mem-index-image'].value = wf.topvars['bwaMemIndexImage']
            filter_align_task.outputs['out'].report = True

        # vep 注释
        depend_task = filter_align_task or filter_task
        vep_task, args = wf.add_task(vep(tumor), tag=tumor, depends=[depend_task])
        args['input_file'].value = depend_task.outputs['out']
        args['fasta'].value = wf.topvars['ref']
        args['refseq'].value = True
        args['dir_cache'].value = top_vars['vep_cache']
        args['dir_plugins'].value = top_vars['vep_plugin']
        vep_task.outputs['out_vcf'].report = True
        vep_task.outputs['out_vcf_idx'].report = True
        # ----end mutect2---------------------------------------------------
        # 自研脚本估计错误率
        error_stat_task, args = wf.add_task(stat_context_seq_error(), tag=tumor)
        args['bam'].value = tumor_bam_task.outputs['output']
        args['genome'].value = wf.topvars['ref']
        args['bed'].value = wf.topvars['bed']
        args['out_prefix'].value = tumor

        # varidct结果过滤和整理
        filter_task, args = wf.add_task(VcfFilter(), tag=tumor)
        args['vcf'].value = tumor_vardict_vep_task.outputs['out_vcf']
        args['error_rate_file'].value = error_stat_task.outputs['context_error_rate']
        args['genome'].value = wf.topvars['ref']
        args['bed'].value = wf.topvars['bed']
        args['normal_vcf'].value = normal_vcf
        args['out_prefix'].value = tumor
        vardict_filter_task_ids.append(filter_task.task_id)

        # mutect2结果过滤和整理
        filter_task, args = wf.add_task(VcfFilter(), tag=tumor)
        args['vcf'].value = vep_task.outputs['out_vcf']
        args['error_rate_file'].value = error_stat_task.outputs['context_error_rate']
        args['genome'].value = wf.topvars['ref']
        args['bed'].value = wf.topvars['bed']
        args['normal_vcf'].value = normal_vcf
        args['out_prefix'].value = tumor
        mutect2_filter_task_ids.append(filter_task.task_id)

    # end
    wf.run()
    # tidy
    if wf.success:
        xls_lst = []
        for tid in vardict_filter_task_ids:
            xls_file = wf.tasks[tid].outputs['final_xls'].value
            xls_file = os.path.join(wf.wkdir, wf.tasks[tid].name, xls_file)
            xls_lst.append(xls_file)
        out_file = os.path.join(wf.wkdir, 'Outputs', 'All.vardict.variants.xlsx')
        df = pd.concat([pd.read_excel(xls_file, sheet_name='Sheet1') for xls_file in xls_lst])
        df.to_excel(out_file, index=False)

        xls_lst = []
        for tid in mutect2_filter_task_ids:
            xls_file = wf.tasks[tid].outputs['final_xls'].value
            xls_file = os.path.join(wf.wkdir, wf.tasks[tid].name, xls_file)
            xls_lst.append(xls_file)
        out_file = os.path.join(wf.wkdir, 'Outputs', 'All.mutect2.variants.xlsx')
        df = pd.concat([pd.read_excel(xls_file, sheet_name='Sheet1') for xls_file in xls_lst])
        df.to_excel(out_file, index=False)


if __name__ == '__main__':
    pipeline()

