import os
import time

script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar
from utils.get_fastq_info import get_fastq_info
__author__ = 'gdq'

"""references:
## https://github.com/gatk-workflows/gatk4-data-processing/tree/2.1.1
## https://github.com/gatk-workflows/gatk4-germline-snps-indels
## https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl
## https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl
"""


def fastp():
    cmd = Command()
    cmd.meta.version = '0.23.2'
    cmd.meta.name = 'fastp'
    cmd.meta.function = 'fastq QC, adapter trimming'
    cmd.meta.desc = """
    fastp is a tool used in bioinformatics for the quality control and preprocessing of raw sequence data. 
    fastp is known for its speed and efficiency, and it can process data in parallel, making it suitable for large datasets.
    fastp provides several key functions:
    * It can filter out low-quality reads, which are sequences that have a high probability of containing errors. This is done based on quality scores that are assigned to each base in a read.
    * It can trim adapter sequences, which are artificial sequences added during the preparation of sequencing libraries and are not part of the actual sample's genome.
    * It provides comprehensive quality control reports, including information on sequence quality, GC content, sequence length distribution, and more.
    """
    cmd.meta.source = 'https://github.com/OpenGene/fastp'
    cmd.runtime.image = 'gudeqing/dnaseq:1.0'
    cmd.runtime.memory = 8 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'fastp'
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', level='optional', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=cmd.runtime.cpu, desc='thread number')
    cmd.args['other_args'] = Argument(prefix='', default='', desc="other arguments you want to use, such as '-x val'")
    cmd.args['out1'] = Argument(prefix='-o ', type='outstr', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', level='optional', type='outstr', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='-h ', type='outstr', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', type='outstr', desc='json report file')
    # 定义输出,这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out1'] = Output(value="{out1}")
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd


def creat_ref_dict():
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'CreateSequenceDictionary'
    cmd.meta.desc = 'Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. '
    cmd.runtime.image = 'gudeqing/dnaseq:1.0'
    cmd.runtime.memory = 5 * 1024 ** 3
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
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'BedToIntervalList'
    cmd.meta.function = 'A gatk tool to create interval file from bed file'
    cmd.meta.desc = 'A gatk tool to create interval file from bed file'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360036726611-BedToIntervalList-Picard-'
    cmd.runtime.image = 'gudeqing/dnaseq:1.0'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g BedToIntervalList'
    cmd.args['bed'] = Argument(prefix='-I ', type='infile', desc='input bed file')
    cmd.args['ref_dict'] = Argument(prefix='-SD ', type='infile', desc='input sequence dictionary file')
    cmd.args['out'] = Argument(prefix='-O ',  desc='output interval file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def build_bwa_index():
    cmd = Command()
    cmd.meta.name = 'buildBwaIndex'
    cmd.meta.function = 'build bwa-mem2 index'
    cmd.meta.desc = 'bwa index and create sequence dictionary and fasta fai file'
    cmd.meta.source = 'https://github.com/bwa-mem2/bwa-mem2'
    cmd.meta.version = '0.7.17'
    cmd.runtime.image = 'gudeqing/dnaseq:1.0'
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
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'FastqToSam'
    cmd.meta.function = 'Converts a FASTQ file to an unaligned BAM or SAM file'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5257906838811-FastqToSam-Picard-'
    cmd.meta.desc = 'Converts a FASTQ file to an unaligned BAM or SAM file. Output read records will contain the original base calls and quality scores will be translated depending on the base quality score encoding: FastqSanger, FastqSolexa and FastqIllumina.'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.tool = 'gatk --java-options -Xmx12g FastqToSam'
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


def uBam2FastqBwaMem(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'uBam2FastqBwaMem'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4409897234331-SamToFastq-Picard-'
    cmd.meta.function = 'ubam to fastq and then mapping'
    cmd.meta.desc = 'convert ubam to fastq and then mapping to genome'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 15*1024**3
    cmd.runtime.cpu = 4
    cmd.args['_fix0'] = Argument(type='fix', value='set -o pipefail; gatk SamToFastq')
    cmd.args['ubam'] = Argument(prefix='-I ', type='infile', desc='input ubam file')
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
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'MergeBamAlignment'
    cmd.meta.function = 'merge bam alignment'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5257880651931-MergeBamAlignment-Picard-'
    cmd.meta.desc = """
    A command-line tool for merging BAM/SAM alignment info from a third-party aligner with the data in an unmapped BAM file, 
    producing a third BAM file that has alignment data (from the aligner) and all the remaining data from the unmapped BAM. 
    Quick note: this is not a tool for taking multiple sam files and creating a bigger file by merging them. 
    """
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g MergeBamAlignment'
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
    cmd.args['PROGRAM_GROUP_COMMAND_LINE'] = Argument(prefix='--PROGRAM_GROUP_COMMAND_LINE ', default='"bwa-mem2 mem -M -Y -p -v 3 -K 10000000 -t 4 ref.fa"')
    cmd.args['PROGRAM_GROUP_NAME'] = Argument(prefix='--PROGRAM_GROUP_NAME ', default='"bwamem"')
    cmd.args['tmpdir'] = Argument(prefix='--TMP_DIR ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def MarkDuplicates(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'MarkDuplicates'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5358880192027-MarkDuplicates-Picard-'
    cmd.meta.function = 'locates and tags duplicate reads in a BAM or SAM file'
    cmd.meta.desc = """
    This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for additional notes on PCR duplication artifacts. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates.
    The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).
    The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, please see the following blog post for additional information.
    """
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g MarkDuplicates'
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', multi_times=True, desc='input bam file list')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.unsorted.dup_marked.bam', desc='output bam file')
    cmd.args['METRICS_FILE'] = Argument(prefix='--METRICS_FILE ', default=f'{sample}.dup_metrics.txt')
    cmd.args['VALIDATION_STRINGENCY'] = Argument(prefix='--VALIDATION_STRINGENCY ', default='SILENT')
    cmd.args['OPTICAL_DUPLICATE_PIXEL_DISTANCE'] = Argument(prefix='--OPTICAL_DUPLICATE_PIXEL_DISTANCE ', default=2500, desc='The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and models, users should experiment to find what works best.')
    cmd.args['ASSUME_SORT_ORDER'] = Argument(prefix='--ASSUME_SORT_ORDER ', default='"queryname"')
    cmd.args['tmpdir'] = Argument(prefix='--TMP_DIR ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def SortAndFixTags(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'SortAndFixTags'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5257874207131-SetNmMdAndUqTags-Picard-'
    cmd.meta.function = 'sort bam and fix tags'
    cmd.meta.desc = """
    This tool takes in a coordinate-sorted SAM or BAM and calculatesthe NM, MD, and UQ tags by comparing with the reference.
    This may be needed when MergeBamAlignment was run with SORT_ORDER other than 'coordinate' and thuscould not fix these tags then.
    The input must be coordinate sorted in order to run. If specified,the MD and NM tags can be ignored and only the UQ tag be set.
    """
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['_sort_sam'] = Argument(type='fix', value='set -o pipefail; gatk --java-options -Xmx15g SortSam')
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', desc='input bam file list')
    cmd.args['_OUTPUT'] = Argument(prefix='--OUTPUT ', default='/dev/stdout', desc='output bam file')
    cmd.args['SORT_ORDER'] = Argument(prefix='--SORT_ORDER ', default='"coordinate"')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='false')
    cmd.args['tmpdir'] = Argument(prefix='--TMP_DIR ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.args['_fix_tags'] = Argument(type='fix', value='| gatk --java-options -Xmx15g SetNmMdAndUqTags --CREATE_INDEX true --INPUT /dev/stdin')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.sorted.dup_marked.bam', desc='output bam file')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='--REFERENCE_SEQUENCE ', type='infile', desc='reference fasta file')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    cmd.outputs['out_idx'] = Output(value='{OUTPUT}.bai')
    return cmd


def CreateSequenceGroupingTSV(ref_dict):
    # print('create the Sequencing Groupings used for BQSR and PrintReads Scatter.')
    with open(ref_dict, "r") as ref_dict_file:
        sequence_tuple_list = []
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        # 获得最长contig的长度
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]

    # with open(f"sequence_grouping_with_unmapped.txt", "w") as tsv_file_with_unmapped:
    #     tsv_file_with_unmapped.write(tsv_string)
    #     tsv_file_with_unmapped.close()

    sequence_grouping = [x.split('\t') for x in tsv_string.strip().split('\n')]
    sequence_grouping_with_unmapped = sequence_grouping.copy() + [['unmapped']]
    return sequence_grouping, sequence_grouping_with_unmapped


def BaseRecalibrator(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'BaseRecalibrator'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5358896138011-BaseRecalibrator'
    cmd.meta.function = 'Generates recalibration table for Base Quality Score Recalibration (BQSR)'
    cmd.meta.desc = """"
    This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. Since there is a large amount of data one can then calculate an empirical probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a table (of the several covariate values, num observations, num mismatches, empirical quality score).
    """
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 5*1024**3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g BaseRecalibrator'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['INPUT'] = Argument(prefix='-I ', type='infile', desc='input bam file')
    cmd.args['use-original-qualities'] = Argument(prefix='--use-original-qualities', type='bool', default=True)
    cmd.args['OUTPUT'] = Argument(prefix='-O ', default=f'{sample}.recal_table', desc='The output recalibration table file to create')
    cmd.args['known-sites'] = Argument(prefix='--known-sites ', type='infile', multi_times=True, desc='One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.')
    cmd.args['intervals'] = Argument(prefix='--intervals ', type='str', multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def GatherBQSRReports(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'GatherBQSRReports'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4413072344731-GatherBQSRReports'
    cmd.meta.function = 'Gathers scattered BQSR recalibration reports into a single file'
    cmd.meta.desc = """
    Gathers scattered BQSR recalibration reports into a single file
    """
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 3*1024**3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g GatherBQSRReports'
    cmd.args['INPUT'] = Argument(prefix='-I ', type='infile', multi_times=True, desc='input bqsr file list')
    cmd.args['OUTPUT'] = Argument(prefix='-O ', default=f'{sample}.recal_table', desc='The output recalibration table file to create')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def ApplyBQSR(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'ApplyBQSR'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/9570337264923-ApplyBQSR'
    cmd.meta.function = 'Apply Base Quality Score Recalibration (BQSR) model'
    cmd.meta.desc = 'Apply Base Quality Score Recalibration (BQSR) model'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g ApplyBQSR'
    cmd.args['INPUT'] = Argument(prefix='-I ', type='infile', desc='input bam file')
    cmd.args['OUTPUT'] = Argument(prefix='-O ', default=f'{sample}.recalibrated.bam', desc='The output recalibration table file to create')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['intervals'] = Argument(prefix='-L ', type='str', multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['bqsr'] = Argument(prefix='-bqsr ', type='infile', desc='Input recalibration table for BQSR')
    cmd.args['static-quantized-quals'] = Argument(prefix='--static-quantized-quals ', type='int', multi_times=True, default=[10, 20, 30])
    cmd.args['add-output-sam-program-record'] = Argument(prefix='--add-output-sam-program-record', type='bool', default=True)
    cmd.args['use-original-qualities'] = Argument(prefix='--use-original-qualities', type='bool', default=True)
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def GatherBamFiles(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5358848717595-GatherBamFiles-Picard-'
    cmd.meta.name = 'GatherBamFiles'
    cmd.meta.function = 'Concatenate efficiently BAM files that resulted from a scattered parallel analysis.'
    cmd.meta.desc = """
    This tool performs a rapid "gather" or concatenation on BAM files. This is often needed in operations that have been run in parallel across genomics regions by scattering their execution across computing nodes and cores thus resulting in smaller BAM files.
    """
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g GatherBamFiles'
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', multi_times=True, desc='input bam file')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.recalibrated.bam', desc='output bam file')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', level='optional', desc='reference fasta file')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='true')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    cmd.outputs['out_idx'] = Output(value='{OUTPUT}.bai')
    return cmd


def SplitIntervals(scatter_number=10):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'SplitIntervals'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4418054374939-SplitIntervals'
    cmd.meta.function = 'Split intervals into sub-interval files.'
    cmd.meta.desc = 'This tool takes in intervals via the standard arguments of IntervalArgumentCollection and splits them into interval files for scattering. The resulting files contain equal number of bases.'
    cmd.runtime.image = 'gudeqing/dnaseq:1.0'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g SplitIntervals'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['intervals'] = Argument(prefix='-L ', level='optional', type='infile', multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['scatter'] = Argument(prefix='-scatter ', type='fix', value=scatter_number, desc='number of output interval files to split into')
    cmd.args['mode'] = Argument(prefix='-mode ', level='optional', desc='How to divide intervals.')
    cmd.args['interval-merging-rule'] = Argument(prefix='--interval-merging-rule ', default='OVERLAPPING_ONLY')
    cmd.args['outdir'] = Argument(prefix='-O ', default='intervals_folder', desc='The directory into which to write the scattered interval sub-directories.')
    for i in range(scatter_number):
        cmd.outputs[f'out{i}'] = Output(value='{outdir}/'+f'{i:04d}-scattered.interval_list', type='outfile')
    return cmd


def Mutect2(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'Mutect2'
    cmd.meta.function = 'Call somatic short mutations'
    cmd.meta.desc = 'Call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations.'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.memory = 8*1024**3
    cmd.runtime.cpu = 3
    cmd.runtime.tool = 'gatk --java-options -Xmx15g Mutect2'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['tumor_bam'] = Argument(prefix='-I ', type='infile', desc='tumor bam')
    cmd.args['normal_bam'] = Argument(prefix='-I ', type='infile', level='optional', desc='normal bam')
    cmd.args['tumor_name'] = Argument(prefix='-tumor ', desc='tumor sample name')
    cmd.args['normal_name'] = Argument(prefix='-normal ', level='optional', desc='normal sample name')
    cmd.args['germline-resource'] = Argument(prefix='--germline-resource ', type='infile', level='optional', desc='optional database of known germline variants (and its index) (see http://gnomad.broadinstitute.org/downloads)')
    cmd.args['pon'] = Argument(prefix='-pon ', type='infile', level='optional', desc='')
    cmd.args['intervals'] = Argument(prefix='-L ', type='infile', level='optional', multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['string_intervals'] = Argument(prefix='-L ', level='optional', desc='interval string such as "chrM"')
    cmd.args['alleles'] = Argument(prefix='--alleles ', type='infile', level='optional', desc='The set of alleles to force-call regardless of evidence')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.vcf.gz', desc='output vcf')
    cmd.args['bam_output'] = Argument(prefix='--bam-output ', level='optional', desc='output bam file')
    cmd.args['f1r2-tar-gz'] = Argument(prefix='--f1r2-tar-gz ', type='infile', default=f'{prefix}.f1r2.tar.gz', desc='If specified, collect F1R2 counts and output files into this tar.gz file')
    cmd.args['mitochondria'] = Argument(prefix='--mitochondira', type='bool', desc='if to turn on mitochondria mode. Specifically, the mode sets --initial-tumor-lod to 0, --tumor-lod-to-emit to 0, --af-of-alleles-not-in-resource to 4e-3, and the advanced parameter --pruning-lod-thres')
    cmd.args['tmpdir'] = Argument(prefix='--tmp-dir ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{out}')
    cmd.outputs['f1r2'] = Output(value='{f1r2-tar-gz}')
    cmd.outputs['stats'] = Output(value='{out}.stats')
    return cmd


def bcftools_norm():
    cmd = Command()
    cmd.meta.version = '1.13'
    cmd.meta.name = 'VcfLeftNorm'
    cmd.meta.source = 'https://samtools.github.io/bcftools/bcftools.html'
    cmd.meta.function = 'Left-align and normalize indels'
    cmd.meta.desc = "Left-align and normalize indels; check if REF alleles match the reference; split multiallelic sites into multiple rows; recover multiallelics from multiple rows"
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.tool = "bcftools norm"
    cmd.args['fasta-ref'] = Argument(prefix='-f ', type='infile', desc='reference fasta file')
    cmd.args['multiallelics'] = Argument(prefix='-m ', level='optional', desc='Split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]' )
    cmd.args['out'] = Argument(prefix='-o ', desc='Write output to a file [standard output]')
    cmd.args['output-type'] = Argument(prefix='--output-type ', default='v', desc="'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]")
    cmd.args['threads'] = Argument(prefix='--threads ', default=2, desc="Use multithreading with <int> worker threads [2]")
    cmd.args['vcf'] = Argument(type='infile', desc='input vcf file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def GetPileupSummaries(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'GetPileupSummaries'
    cmd.meta.function = 'Tabulates pileup metrics for inferring contamination'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5358860217115-GetPileupSummaries'
    cmd.meta.desc = 'Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.tool = 'gatk --java-options -Xmx15g GetPileupSummaries'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='BAM/SAM/CRAM file containing reads')
    cmd.args['interval-set-rule'] = Argument(prefix='--interval-set-rule ', default='INTERSECTION', desc='Set merging approach to use for combining interval inputs')
    cmd.args['intervals'] = Argument(prefix='-L ', type='infile', multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['variants_for_contamination'] = Argument(prefix='-V ', type='infile', desc='A VCF file containing variants and allele frequencies')
    cmd.args['out'] = Argument(prefix='-O ', type='infile', default=f'{prefix}.tumor-pileups.table', desc='The output table')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def MergeVcfs(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'MergeVcfs'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5358883370011-MergeVcfs-Picard-'
    cmd.meta.function = 'Combines multiple variant files into a single variant file.'
    cmd.meta.desc = 'Combines multiple variant files into a single variant file.'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360056969852-MergeVcfs-Picard-'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g MergeVcfs'
    cmd.args['inputs'] = Argument(prefix='-I ', type='infile', multi_times=True, desc='input vcf list')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.vcf.gz', desc='The merged VCF or BCF file. File format is determined by file extension.')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def SortBam():
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'SortBam'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5358884334747-SortSam-Picard-'
    cmd.meta.function = 'sorts the input SAM or BAM file'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.meta.desc = 'This tool sorts the input SAM or BAM file by coordinate, queryname (QNAME), or some other property of the SAM record. The SortOrder of a SAM/BAM file is found in the SAM file header tag @HD in the field labeled SO.'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360036366932-SortSam-Picard-'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g SortSam'
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='BAM/SAM/CRAM file containing reads')
    cmd.args['out'] = Argument(prefix='-O ', desc='Sorted BAM or SAM output file.')
    cmd.args['SORT_ORDER'] = Argument(prefix='--SORT_ORDER ', default='coordinate', desc='Sort order of output file.')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='true')
    cmd.args['tmpdir'] = Argument(prefix='--TMP_DIR ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def MergeMutectStats(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'MergeMutectStats'
    cmd.meta.function = 'Merge Mutect2 Stats file'
    cmd.meta.desc = 'Merge Mutect2 Stats file'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.tool = 'gatk --java-options -Xmx15g MergeMutectStats'
    cmd.args['stats'] = Argument(prefix='-stats ', type='infile', multi_times=True)
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.vcf.stats', desc='output merged stat files')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def GatherPileupSummaries(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'GatherPileupSummaries'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-'
    cmd.meta.function = 'Gather Pileup Summaries'
    cmd.meta.desc = 'Gather Pileup Summaries'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.tool = 'gatk --java-options -Xmx15g GatherPileupSummaries'
    cmd.args['sequence-dictionary'] = Argument(prefix='--sequence-dictionary ', type='infile')
    cmd.args['inputs'] = Argument(prefix='-I ', type='infile', multi_times=True)
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.pileup.table', desc='output merged pileup summary file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def LearnReadOrientationModel(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.function = 'Learn Read Orientation Model'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-'
    cmd.meta.name = 'LearnReadOrientationModel'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.meta.desc = 'Learn the prior probability of read orientation artifact from the output of CollectF1R2Counts of Mutect2'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g LearnReadOrientationModel'
    cmd.args['inputs'] = Argument(prefix='-I ', type='infile', multi_times=True, desc='One or more .tar.gz containing outputs of CollectF1R2Counts')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.artifact-priors.tar.gz', desc='tar.gz of artifact prior tables')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def CalculateContamination(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'CalculateContamination'
    cmd.meta.function = 'Calculate Contamination'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-'
    cmd.meta.desc = 'Calculates the fraction of reads coming from cross-sample contamination, given results from GetPileupSummaries. The resulting contamination table is used with FilterMutectCalls.'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g CalculateContamination'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.args['tumor_pileups'] = Argument(prefix='-I ', type='infile', desc='input pileup table')
    cmd.args['normal_pileups'] = Argument(prefix='-matched ', type='infile', desc='The matched normal input table')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.contamination.table', desc='output table')
    cmd.args['tumor-segmentation'] = Argument(prefix='--tumor-segmentation ', default=f'{prefix}.tumor_segmentation', desc='The output table containing segmentation of the tumor by minor allele fraction')
    cmd.outputs['out'] = Output(value='{out}')
    cmd.outputs['tumor-segmentation'] = Output(value='{tumor-segmentation}')
    return cmd


def FilterMutectCalls(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'FilterMutectCalls'
    cmd.meta.function = 'FilterMutect Calls'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-'
    cmd.meta.desc = 'FilterMutectCalls applies filters to the raw output of Mutect2'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g FilterMutectCalls'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='A VCF file containing variants')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
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
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'FilterAlignmentArtifacts'
    cmd.meta.function = 'Filter Alignment Artifacts'
    cmd.meta.desc = 'Alignment artifacts can occur whenever there is sufficient sequence similarity between two or more regions in the genome to confuse the alignment algorithm. This can occur when the aligner for whatever reason overestimate how uniquely a read maps, thereby assigning it too high of a mapping quality. It can also occur through no fault of the aligner due to gaps in the reference, which can also hide the true position to which a read should map. By using a good alignment algorithm (the GATK wrapper of BWA-MEM), giving it sensitive settings (which may have been impractically slow for the original bam alignment) and mapping to the best available reference we can avoid these pitfalls. The last point is especially important: one can (and should) use a BWA-MEM index image corresponding to the best reference, regardless of the reference to which the bam was aligned.'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4418051467035-FilterAlignmentArtifacts-EXPERIMENTAL-'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g FilterAlignmentArtifacts'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='A VCF file containing variants')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='input bam file')
    cmd.args['bwa-mem-index-image'] = Argument(prefix='--bwa-mem-index-image ', type='infile', desc='BWA-mem index image')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.align_artifacts_filtered.vcf.gz', desc='output vcf file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def vep(sample):
    cmd = Command()
    cmd.meta.version = 'vep104.3'
    cmd.meta.name = 'VEP'
    cmd.meta.function = 'vcf annotation'
    cmd.meta.desc = """
    VEP determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions.
    Simply input the coordinates of your variants and the nucleotide changes to find out the:
    * Genes and Transcripts affected by the variants
    * Location of the variants (e.g. upstream of a transcript, in coding sequence, in non-coding RNA, in regulatory regions)
    * Consequence of your variants on the protein sequence (e.g. stop gained, missense, stop lost, frameshift), see variant consequences
    * Known variants that match yours, and associated minor allele frequencies from the 1000 Genomes Project
    * SIFT and PolyPhen-2 scores for changes to protein sequence
    * more...
    """
    cmd.runtime.image = 'gudeqing/dnaseq:1.0'
    cmd.runtime.tool = 'vep'
    cmd.args['input_file'] = Argument(prefix='-i ', type='infile', desc='input file')
    cmd.args['fasta'] = Argument(prefix='--fasta ', type='infile', desc="Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run VEP with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache).")
    cmd.args['output_file'] = Argument(prefix='-o ', default=f'{sample}.vep.vcf.gz', desc='output file')
    cmd.args['output_format'] = Argument(prefix='--', range={'vcf', 'json', 'tab'}, default='vcf', desc="If we choose to write output in VCF format. Consequences are added in the INFO field of the VCF file, using the key 'CSQ'. Data fields are encoded separated by '|'; the order of fields is written in the VCF header. Output fields in the 'CSQ' INFO field can be selected by using --fields.")
    cmd.args['compress_output'] = Argument(prefix='--compress_output ', default='bgzip', desc="Writes output compressed using either gzip or bgzip")
    cmd.args['force_overwrite'] = Argument(prefix="--force_overwrite ", type='bool', default=True, desc="Force overwriting of output file")
    cmd.args['fork'] = Argument(prefix='--fork ', type='int', default=4, desc='Use forking(multi-cpu/threads) to improve script runtime')
    cmd.args['species'] = Argument(prefix='--species ', default='homo_sapiens', desc='Species for your data. This can be the latin name e.g. homo_sapiens or any Ensembl alias e.g. mouse.')
    cmd.args['assembly_version'] = Argument(prefix='--assembly ', default='GRCh38', desc='Select the assembly version to use if more than one available.')
    cmd.args['dir_cache'] = Argument(prefix='--dir_cache ', type='indir', desc='Specify the cache directory to use')
    cmd.args['dir_plugins'] = Argument(prefix='--dir_plugins ', type='indir', level='optional', desc='Specify the plugin directory to use')
    cmd.args['merged'] = Argument(prefix='--merged ', type='bool', default=False, desc='Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.')
    cmd.args['refseq'] = Argument(prefix='--refseq ', type='bool', default=False, desc='Specify this option if you have installed the RefSeq cache in order for VEP to pick up the alternate cache directory. This cache contains transcript objects corresponding to RefSeq transcripts. Consequence output will be given relative to these transcripts in place of the default Ensembl transcripts')
    cmd.args['stats_file'] = Argument(prefix='--stats_file ', default=f'{sample}.vep.summary.html', desc='Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>.')
    cmd.args['cache'] = Argument(prefix='--cache ', type='bool', default=True, desc='Enables use of cache')
    cmd.args['offline'] = Argument(prefix='--offline ', type='bool', default=True, desc='Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation')
    cmd.args['merged'] = Argument(prefix='--merged ', type='bool', default=False, desc='Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.')
    cmd.args['plugins'] = Argument(prefix='--plugin ', multi_times=True, level='optional', default=['Frameshift', 'Wildtype'], desc='Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times')
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
    cmd.args['coding_only'] = Argument(prefix='--coding_only ', type='bool', default=False, desc="Only return consequences that fall in the coding regions of transcripts. Not used by default")
    cmd.args['pick'] = Argument(prefix='--pick', type='bool', default=False, desc="Pick one line or block of consequence data per variant, including transcript-specific columns. This is the best method to use if you are interested only in one consequence per variant")
    cmd.args['flag_pick'] = Argument(prefix='--flag_pick ', type='bool', default=True, desc="As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others.")
    cmd.args['filter_common'] = Argument(prefix='--filter_common ', type='bool', default=False, desc="Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters.")
    cmd.args['other_args'] = Argument(default='', desc='specify other arguments that you want to append to the command')
    cmd.args['_create_index'] = Argument(value='&& tabix *vcf.gz', type='fix')
    cmd.outputs['out_vcf'] = Output(value='{output_file}')
    cmd.outputs['out_vcf_idx'] = Output(value='{output_file}.tbi')
    return cmd


def Haplotyper(sample):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'Haplotyper'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5257862382875-HaplotypeCaller'
    cmd.meta.function = 'Call germline SNPs and indels via local re-assembly of haplotypes'
    cmd.meta.desc = """
    The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.
    """
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.tool = 'gatk --java-options -Xmx15g HaplotypeCaller'
    cmd.args['intervals'] = Argument(prefix='-L ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='reccaled tumor and normal bam list')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['emit_mode'] = Argument(prefix='-ERC ', default='GVCF', range=['NONE', 'GVCF', 'BP_RESOLUTION'], desc='The reference confidence mode makes it possible to emit a per-bp or summarized confidence estimate for a site being strictly homozygous-reference.')
    cmd.args['ploidy'] = Argument(prefix='-ploidy ', type='int', default=2, desc='determines the ploidy number of the sample being processed. The default value is 2.')
    cmd.args['annotation-group'] = Argument(prefix='-G ', multi_times=True, default=['StandardAnnotation', 'StandardHCAnnotation', 'AS_StandardAnnotation'])
    cmd.args['gvcf-gq-bands'] = Argument(prefix='-GQB ', multi_times=True, default=[10, 20, 30, 40, 50, 60, 70, 80, 90], desc='Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)')
    cmd.args['tmpdir'] = Argument(prefix='--tmp-dir ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.args['out_vcf'] = Argument(prefix='-O ', default=f'{sample}.g.vcf.gz', desc='output vcf file')
    cmd.outputs['out'] = Output(value='{out_vcf}')
    cmd.outputs['out_idx'] = Output(value='{out_vcf}.tbi')
    return cmd


def GenomicsDBImport():
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5358869876891-GenomicsDBImport'
    cmd.meta.name = 'GenomicsDBImport'
    cmd.meta.function = 'Import single-sample GVCFs into GenomicsDB before joint genotyping.'
    cmd.meta.desc = 'The GATK4 Best Practice Workflow for SNP and Indel calling uses GenomicsDBImport to merge GVCFs from multiple samples.  In brief, GenomicsDB utilises a data storage system optimized for storing/querying sparse arrays. Genomics data is typically sparse in that each sample has few variants with respect to the entire reference genome. GenomicsDB contains specialized code for genomics applications, such as VCF parsing and INFO field annotation calculation.'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g GenomicsDBImport'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.args['gvcfs'] = Argument(prefix='-V ', type='infile', multi_times=True, desc='input gVcf files')
    cmd.args['genomicsdb-workspace-path'] = Argument(prefix='--genomicsdb-workspace-path ', default='genomicsdb', desc='Workspace for GenomicsDB. Can be a POSIX file system absolute or relative path or a HDFS/GCS URL. Use this argument when creating a new GenomicsDB workspace. Either this or genomicsdb-update-workspace-path must be specified. Must be an empty or non-existent directory.')
    cmd.args['batch-size'] = Argument(prefix='--batch-size ', default=50)
    cmd.args['interval'] = Argument(prefix='-L ', type='infile', level='optional')
    cmd.args['reader-threads'] = Argument(prefix='--reader-threads ', default=5)
    cmd.args['merge-input-intervals'] = Argument(prefix='--merge-input-intervals', type='bool', default=True)
    cmd.args['consolidate'] = Argument(prefix='--consolidate', type='bool', default=True)
    cmd.args['tmpdir'] = Argument(prefix='--tmp-dir ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{genomicsdb-workspace-path}')
    return cmd


def GenotypeGVCFs(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'GenotypeGVCFs'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5257900485659-GenotypeGVCFs'
    cmd.meta.function = 'Perform joint genotyping on one or more samples pre-called with HaplotypeCaller'
    cmd.meta.desc = """
    Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
    This tool is designed to perform joint genotyping on a single input, which may contain one or many samples. In any case, the input samples must possess genotype likelihoods produced by HaplotypeCaller with `-ERC GVCF` or `-ERC BP_RESOLUTION`.
    * Input:
    The GATK4 GenotypeGVCFs tool can take only one input track. Options are 1) a single single-sample GVCF 2) a single multi-sample GVCF created by CombineGVCFs or 3) a GenomicsDB workspace created by GenomicsDBImport. A sample-level GVCF is produced by HaplotypeCaller with the `-ERC GVCF` setting.
    * Output:
    A final VCF in which all samples have been jointly genotyped.
    """
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.runtime.tool = 'gatk --java-options -Xmx15g GenotypeGVCFs'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.genotyped.gvcf', desc='output vcf file name')
    cmd.args['dbsnp'] = Argument(prefix='-D ', type='infile', desc='input dbsnp vcf file')
    cmd.args['annotation-group'] = Argument(prefix='-G ', multi_times=True, default=['StandardAnnotation', 'AS_StandardAnnotation'])
    cmd.args['only-output-calls-starting-in-intervals'] = Argument(prefix='--only-output-calls-starting-in-intervals', type='bool', default=True)
    cmd.args['gendb'] = Argument(prefix='-V gendb://', type='indir')
    cmd.args['interval'] = Argument(prefix='-L ', type='infile', level='optional')
    cmd.args['merge-input-intervals'] = Argument(prefix='--merge-input-intervals', type='bool', default=True)
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def VariantFiltration(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'VariantFiltration'
    cmd.meta.function = 'Filter variant calls based on INFO and/or FORMAT annotations'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/5257863504923-VariantFiltration'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.meta.desc = 'This tool is designed for hard-filtering variant calls based on certain criteria. Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records will be preserved in the output unless their removal is requested in the command line.'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g VariantFiltration'
    cmd.args['filters'] = Argument(prefix='', default=['--filter-name ExcessHet --filter-expression "ExcessHet > 54.69"',], multi_times=True, desc='filtering expressions, such as "AB < 0.2"')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.hard_filtered.vcf.gz', desc='output vcf file')
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='input vcf file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def MakeSitesOnlyVcf(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'MakeSitesOnlyVcf'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360045799552-MakeSitesOnlyVcf-Picard-'
    cmd.meta.function = 'Creates a VCF that contains all the site-level information for all records in the input VCF but no genotype information.'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.meta.desc = 'This tool reads a VCF/VCF.gz/BCF and removes all genotype information from it while retaining all site level information, including annotations based on genotypes (e.g. AN, AF).'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g MakeSitesOnlyVcf'
    cmd.args['vcf'] = Argument(prefix='-I ', type='infile', desc='input vcf file')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.site_only.vcf.gz', desc='output vcf file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def IndelsVariantRecalibrator(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.name = 'IndelsRecalibrator'
    cmd.meta.function = 'Indels Recalibration'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4413056319131-VariantRecalibrator'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.meta.desc = 'Build a recalibration model to score variant quality for filtering purposes'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g VariantRecalibrator'
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='site only variant filtered input vcf')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.indel.recal', desc='The output recal file used by ApplyVQSR')
    cmd.args['out-tranches-file'] = Argument(prefix='--tranches-file ', default=f'{prefix}.indel.tranches', desc='The output tranches file used by ApplyVQSR')
    cmd.args['tranche'] = Argument(prefix='-tranche ', multi_times=True, default=["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"], desc='recalibration tranche values')
    cmd.args['use-annotation'] = Argument(prefix='-an ', multi_times=True, default=["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"], desc='The names of the annotations which should used for calculations')
    cmd.args['trust-all-polymorphic'] = Argument(prefix='--trust-all-polymorphic', type='bool', default=True)
    cmd.args['use-allele-specific-annotations'] = Argument(prefix='--use-allele-specific-annotations', type='bool', default=True)
    cmd.args['mode'] = Argument(prefix='-mode ', type='fix', default='INDEL', desc='Recalibration mode to employ')
    cmd.args['max-gaussians'] = Argument(prefix='--max-gaussians ', default=4, desc='Max number of Gaussians for the positive model')
    cmd.args['mills'] = Argument(prefix='-resource:mills,known=false,training=true,truth=true,prior=12 ', type='infile', desc='mills resource vcf')
    cmd.args['axiomPloly'] = Argument(prefix='-resource:axiomPoly,known=false,training=true,truth=false,prior=10 ', type='infile', desc='axiomPoly resource vcf')
    cmd.args['dbsnp'] = Argument(prefix='-resource:dbsnp,known=true,training=false,truth=false,prior=2 ', type='infile', desc='dbsnp resource vcf')
    cmd.outputs['out'] = Output(value='{out}')
    cmd.outputs['tranches'] = Output(value='{out-tranches-file}')
    return cmd


def SNPsVariantRecalibrator(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4413056319131-VariantRecalibrator'
    cmd.meta.name = 'SNPsRecalibrator'
    cmd.meta.function = 'SNPs Recalibration'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.meta.desc = 'Build a recalibration model to score variant quality for filtering purposes'
    cmd.runtime.tool = 'gatk --java-options -Xmx15g VariantRecalibrator'
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='site only variant filtered input vcf')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.snp.recal', desc='The output recal file used by ApplyVQSR')
    cmd.args['out-tranches-file'] = Argument(prefix='--tranches-file ', default=f'{prefix}.snp.tranches', desc='The output tranches file used by ApplyVQSR')
    cmd.args['tranche'] = Argument(prefix='-tranche ', multi_times=True, default=["100.0", "99.95", "99.9", "99.8", "99.7", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0"], desc='recalibration tranche values')
    cmd.args['use-annotation'] = Argument(prefix='-an ', multi_times=True, default=["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"], desc='The names of the annotations which should used for calculations')
    cmd.args['trust-all-polymorphic'] = Argument(prefix='--trust-all-polymorphic', type='bool', default=True)
    cmd.args['use-allele-specific-annotations'] = Argument(prefix='--use-allele-specific-annotations', type='bool', default=True)
    cmd.args['mode'] = Argument(prefix='-mode ', type='fix', default='SNP', desc='Recalibration mode to employ')
    cmd.args['max-gaussians'] = Argument(prefix='--max-gaussians ', default=6, desc='Max number of Gaussians for the positive model')
    cmd.args['output-model'] = Argument(prefix='--output-model ', level='optional', desc='If specified, the variant recalibrator will output the VQSR model to this file path.')
    cmd.args['hapmap'] = Argument(prefix='-resource:hapmap,known=false,training=true,truth=true,prior=15 ', type='infile', desc='hapmap resource vcf')
    cmd.args['omni'] = Argument(prefix='-resource:omni,known=false,training=true,truth=true,prior=12 ', type='infile', desc='omni resource vcf')
    cmd.args['1000G'] = Argument(prefix='-resource:1000G,known=false,training=true,truth=false,prior=10 ', type='infile', desc='one thousand genomes resource vcf')
    cmd.args['dbsnp'] = Argument(prefix='-resource:dbsnp,known=true,training=false,truth=false,prior=7 ', type='infile', desc='dbsnp resource vcf')
    cmd.outputs['out'] = Output(value='{out}')
    cmd.outputs['tranches'] = Output(value='{out-tranches-file}')
    cmd.outputs['model'] = Output(value='{output-model}')
    return cmd


def ApplyVQSR(prefix):
    cmd = Command()
    cmd.meta.version = 'GATK v4.2.6.1'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4414602412827-ApplyVQSR'
    cmd.meta.function = 'Apply a score cutoff to filter variants based on a recalibration table'
    cmd.meta.name = 'ApplyVQSR'
    cmd.runtime.image = "gudeqing/dnaseq:1.0"
    cmd.meta.desc = "This tool performs the second pass in a two-stage process called Variant Quality Score Recalibration (VQSR). Specifically, it applies filtering to the input variants based on the recalibration table produced in the first step by VariantRecalibrator and a target sensitivity value, which the tool matches internally to a VQSLOD score cutoff based on the model's estimated sensitivity to a set of true variants."
    cmd.runtime.tool = 'gatk --java-options -Xmx15g ApplyVQSR'
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='input vcf')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.filtered.vcf.gz', desc='The output vcf')
    cmd.args['recal-file'] = Argument(prefix='--recal-file ', type='infile', desc='The input recal file used by ApplyVQSR')
    cmd.args['use-allele-specific-annotations'] = Argument(prefix='--use-allele-specific-annotations ', type='bool', default=True)
    cmd.args['tranches-file'] = Argument(prefix='--tranches-file ', type='infile', desc='The input tranches file used by ApplyVQSR')
    cmd.args['truth-sensitivity-filter-level'] = Argument(prefix='--truth-sensitivity-filter-level ', default=99.0, desc='The truth sensitivity level at which to start filtering')
    cmd.args['create-output-variant-index'] = Argument(prefix='--create-output-variant-index ', default='true', range=['true', 'false'])
    cmd.args['mode'] = Argument(prefix='-mode ', default='SNP', range=['SNP', 'INDEL'], desc="Recalibration mode to employ: 1.) SNP for recalibrating only SNPs (emitting indels untouched in the output VCF); 2.) INDEL for indels; and 3.) BOTH for recalibrating both SNPs and indels simultaneously.")
    cmd.args['tmpdir'] = Argument(prefix='--tmp-dir ', default='.', desc='directorie with space available to be used by this program for temporary storage of working files')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.version = "0.1"
    wf.meta.name = 'GATK-DNASeq-Workflow'
    wf.meta.source = """
    ## https://github.com/gatk-workflows/gatk4-data-processing/tree/2.1.1
    ## https://github.com/gatk-workflows/gatk4-germline-snps-indels
    ## https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl
    ## https://github.com/broadinstitute/warp/blob/develop/tasks/broad/JointGenotypingTasks.wdl
    """
    wf.meta.desc = """
    主要功能：
    当前流程是DNASeq分析流程，参考博得研究所的GATK-Best-Practice流程改写而成，主要用于SNP/Indel突变检测分析。流程包含的主要步骤如下：
    * 使用fastp进行测序接头自动去除
    * 使用BWA进行比对分析
    * 使用GATK检测small SNP/Indel, 同时支持tumor-only和tumor-normal配对模式
    * 使用GATK检测germline突变，如果输入多个normal样本，则会直接进行joint-calling
    * 基于VEP进行突变注释
    
    使用示例：？
    
    主要输入说明（以hg38为例）:
    1. 下载aws提供的公开数据集: aws s3 ls --no-sign-request s3://broad-references/hg38/v0/
        期望得到如下文件:
        Homo_sapiens_assembly38.fasta, 作为"-ref"参数的值
        Homo_sapiens_assembly38.dict
        Homo_sapiens_assembly38.fasta.fai
        dbsnp_146.hg38.vcf.gz, 作为参数'-dbsnp'的值
        dbsnp_146.hg38.vcf.gz.tbi
        Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz, 作为参数'-axiomPoly'的值
        Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
        Mills_and_1000G_gold_standard.indels.hg38.vcf.gz, 作为参数'-mills'的输入
        Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
        hapmap_3.3.hg38.vcf.gz, 作为参数'-hapmap'的输入
        hapmap_3.3.hg38.vcf.gz.tbi
        1000G_omni2.5.hg38.vcf.gz, 作为参数'-omni'的输入
        1000G_omni2.5.hg38.vcf.gz.tbi
        1000G_phase1.snps.high_confidence.hg38.vcf.gz, 作为参数'-G1000'的输入
        1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
        M2PoN_4.0_WGS_for_public.vcf, 作为参数"-pon"的输入
        M2PoN_4.0_WGS_for_public.vcf.idx
        af-only-gnomad.hg38.vcf.gz, 作为参数'-germline_vcf'的输入
        af-only-gnomad.hg38.vcf.gz.tbi
        small_exac_common_3.hg38.vcf.gz, 作为参数'-contamination_vcf'的输入
        small_exac_common_3.hg38.vcf.gz.tbi
        wgs_calling_regions.hg38.interval_list, 如果是WGS测序,可以作为-intervals参数的输入
    2. 使用bwamem2进行索引构建, 如进入容器gudeqing/dnaseq:1.0的容器: 执行/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index your_fasta_file
    3. fastq文件, 可以提供一个fastq所在目录,然后通过参数'-r1_name'和'-r2_name'提供fastq的文件名的正则表达式, 我们依据正则表达式提取样本路径和样本名称
    4. 输入参数call_mode决定分析模式，支持somatic，germline和mix三种模式，mix模式是对肿瘤样本进行somatic突变分析，并基于对照样本进行germline突变分析
    5. 输入参数interval提供捕获区域文件，可以提供bed文件格式或interval文件格式
    """
    # 定义流程输入参数
    wf.init_argparser()
    wf.add_argument('-fastq_info', nargs='+', required=True, help='A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]')
    wf.add_argument('-r1_name', default='(.*).R1.fastq', help="python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'")
    wf.add_argument('-r2_name', default='(.*).R2.fastq', help="python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'")
    wf.add_argument('-exclude_samples', default=tuple(), nargs='+', help='samples to exclude from analysis')
    wf.add_argument('-call_mode', default='somatic', help='somatic: somatic variant calling for tumor samples; germline: germline variant calling for all input samples, this mode is suitable for WES and WGS; mix: somatic variant calling for tumor samples and germline variant calling for paired normal samples')
    wf.add_argument('-pair_info', required=False, help='File containing tumor normal pair info, two-column txt file, first column is tumor sample name. sample not in pair info will be skipped. In "germline" mode, this input is ignored')
    wf.add_argument('-tumor_sample_name', required=False, help='If you only input two paired samples, for convenient, you may use this argument instead of "pair_info" to specify tumor sample. For multiple samples, you should provide sample pairing info file by "pair_info')
    wf.add_argument('-scatter', default=10, help='Scatter number used for interval splitting of variant calling steps')
    wf.add_argument('-intervals', required=False, help="Interval file, support bed file or picard interval file.")
    wf.add_argument('-alleles', required=False, help='The set of alleles to force-call regardless of evidence')
    wf.add_argument('-pon', required=False, help='Panel of normal vcf file for germline variant filtering, this will be required for tumor only analysis')
    wf.add_argument('-ref_dir', default='/enigma/datasets/*/hg38/', help='The directory containing all reference files for gatk input. The following input files need to be in this file')
    wf.add_argument('--ref', default='Homo_sapiens_assembly38.fasta', help='reference fasta file. It should be in "ref_dir"')
    wf.add_argument('--dbsnp', default='dbsnp_146.hg38.vcf.gz', help='dbsnp vcf file')
    wf.add_argument('--axiomPoly', default='Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz', help='high confidence known indel vcf file, such as Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz. 1,249 individuals are included, drawn from the International HapMap Project and 1000 Genomes Project sample collection. https://www.thermofisher.cn/cn/zh/home/life-science/microarray-analysis/microarray-data-analysis/microarray-analysis-sample-data/axiom-exome-sample-data-set.html')
    wf.add_argument('--mills', default="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", help='high confidence known indel vcf file. such as Mills_and_1000G_gold_standard.indels.hg38.vcf')
    wf.add_argument('--hapmap', default="hapmap_3.3.hg38.vcf.gz", help='high confidence known snp vcf file. 来自国际人类单倍体型图计划, 这个数据集包含了大量家系数据，并且有非常严格的质控和严密的实验验证，因此它的准确性是目前公认最高的')
    wf.add_argument('--omni', default="1000G_omni2.5.hg38.vcf.gz", help='high confidence known snp vcf file. 这个数据源自Illumina的Omni基因型芯片，大概2.5百万个位点，它的验证结果常常作为基因型的金标准')
    wf.add_argument('--G1000', default="1000G_phase1.snps.high_confidence.hg38.vcf.gz", help='high confidence known snp vcf file. source from 1000 genomes project')
    wf.add_argument('--germline_vcf', default='af-only-gnomad.hg38.vcf.gz', help='germline vcf, will be used for germline variant filtering')
    wf.add_argument('--contamination_vcf', default='small_exac_common_3.hg38.vcf.gz', help='germline vcf such as small_exac_common_3_b37.vcf, will be used for contamination analysis')
    wf.add_argument('--bwaMemIndexImage', required=False, help='bwa-mem-index-mage for artifact alignment filtering. you may created it with tool BwaMemIndexImageCreator with only fasta as input')
    wf.add_argument('-vep_cache_dir', required=False, help='VEP cache directory')
    wf.add_argument('-vep_plugin_dir', required=False, help='VEP plugin directory')
    wf.parse_args()

    top_vars = dict(
        ref=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.ref), type='infile'),
        dbsnp=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.dbsnp), type='infile'),
        axiomPoly=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.axiomPoly), type='infile'),
        mills=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.mills), type='infile'),
        hapmap=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.hapmap), type='infile'),
        omni=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.omni), type='infile'),
        G1000=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.G1000), type='infile'),
        germline_vcf=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.germline_vcf), type='infile'),
        contamination_vcf=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.contamination_vcf), type='infile'),
        bwaMemIndexImage=TopVar(value=os.path.join(wf.args.ref_dir, wf.args.bwaMemIndexImage) if wf.args.bwaMemIndexImage else None, type='infile'),
        pon=TopVar(value=wf.args.pon, type='infile'),
        vep_cache_dir=TopVar(value=wf.args.vep_cache_dir, type='indir'),
        vep_plugin_dir=TopVar(value=wf.args.vep_plugin_dir, type='indir'),
        intervals=TopVar(value=wf.args.intervals, type='infile'),
        alleles=TopVar(value=wf.args.alleles, type='infile'),
        # hisatgenotype_db=TopVar(value=wf.args.hisatgenotype_db, type='infile')
    )
    wf.add_topvars(top_vars)

    # 提取fastq信息
    fastq_info = get_fastq_info(fastq_info=wf.args.fastq_info, r1_name=wf.args.r1_name, r2_name=wf.args.r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found !')
    # 检测参数
    if wf.args.call_mode not in ['somatic', 'germline', 'mix']:
        raise Exception("call mode is invalid. candidate are from ['somatic', 'germline', 'mix']")
    # 提取配对信息
    pair_list = []
    sample_list = []
    if wf.args.pair_info:
        if wf.args.call_mode == 'germline':
            print('warn: you provide pair information, but "germline" mode is enabled !')
            time.sleep(8)
        with open(wf.args.pair_info) as f:
            for line in f:
                if line.strip():
                    pairs = line.strip('\n').split('\t')[:2]
                    pair_list.append(pairs)
                    sample_list.extend(pairs)
    elif wf.args.tumor_sample_name:
        if len(fastq_info) == 2:
            sample_list = list(fastq_info.keys())
            ctrl_sample = set(sample_list) - {wf.args.tumor_sample_name}
            pair_list = [[wf.args.tumor_sample_name, list(ctrl_sample)[0]]]
        else:
            raise Exception('you used the argument "tumor_sample_name", thus we expect only two input samples!')
        if wf.args.call_mode == 'germline':
            print('warn: you used the argument "tumor_sample_name", but "germline" mode is enabled !')
            time.sleep(8)
    if not sample_list:
        sample_list = list(fastq_info.keys())
    if wf.args.call_mode in ['somatic', 'mix']:
        if wf.args.call_mode == 'somatic':
            print('Analysis mode: somatic variant calling')
        else:
            print('Analysis mode: somatic variant calling for tumor samples and germline variant calling for normal samples')
        if not pair_list:
            print('No pair information provided, we assume all samples are tumor sample')
            pair_list = [[x, 'none'] for x in sample_list]
    else:
        print('Analysis mode: germline variant calling for all input samples')
        print('Caution: This mode is suitable for WES and WGS or large panel. For small panel, error may happen in recalibration stage!')
        pair_list = [['none', x] for x in sample_list]

    # create dict or fai if necessary
    dict_file = wf.topvars['ref'].value.rsplit('.', 1)[0] + '.dict'
    fai_file = wf.topvars['ref'].value + '.fai'
    wf.topvars['ref_dict'] = TopVar(value=dict_file, type='infile')
    if not os.path.exists(dict_file) or not os.path.exists(fai_file):
        print(f'dict or fai file for {wf.args.ref} does not exist, we will try to create one')
        creat_dict_task = creat_ref_dict()
        creat_dict_task.args['ref_fasta'].value = wf.topvars['ref'].value
        tmp_wkdir = os.path.join(wf.args.outdir, "PreparedInputs")
        creat_dict_task.run_now(wkdir=tmp_wkdir, docker=wf.args.docker)
        wf.topvars['ref'].value = creat_dict_task.outputs['ref_genome'].value
        wf.topvars['ref_dict'].value = creat_dict_task.outputs['ref_dict'].value

    # CreateSequenceGroupingTSV
    seq_groups, seq_groups2 = CreateSequenceGroupingTSV(wf.topvars['ref_dict'].value)

    # create interval list if necessary
    input_interval = wf.topvars['intervals'].value
    if input_interval and input_interval.endswith('.bed'):
        bed2interval_task = BedToIntervalList()
        bed2interval_task.args['bed'].value = input_interval
        bed2interval_task.args['ref_dict'].value = wf.topvars['ref_dict'].value
        bed2interval_task.args['out'].value = os.path.basename(input_interval) + '.interval_list'
        bed2interval_task.run_now(wkdir=os.path.join(wf.args.outdir, "PreparedInputs"), docker=wf.args.docker)
        wf.topvars['intervals'].value = bed2interval_task.outputs['out'].value

    # 建bwa索引
    make_index = False
    index_task = None
    if not os.path.exists(wf.topvars['ref'].value + '.0123'):
        make_index = True
        index_task, args = wf.add_task(build_bwa_index(), name='buildIndex')
        if not wf.args.docker:
            args['copy_input_mode'].value = 's'
        args['ref_fasta'].value = wf.topvars['ref']

    # 比对分析
    recal_dict = dict()
    bam_dict = dict()
    for sample, reads in fastq_info.items():
        # 跳过不需要分析的样本
        if sample in wf.args.exclude_samples or (sample not in sample_list):
            print(f'Skip {sample} for it is in excluded sample list or not found in pair info')
            continue
        # fastq预处理和比对 考虑一个样本存在多对fastq的处理
        if len(reads) == 2:
            # pair-end
            r1s, r2s = reads
        else:
            # single-end
            r1s = reads[0]
            r2s = [None]*len(r1s)

        merge_bam_tasks = []
        for ind, (r1, r2) in enumerate(zip(r1s, r2s)):
            if len(r1s) > 1:
                sample_tag = f'{sample}-{ind}'
            else:
                sample_tag = sample
            fastp_task, args = wf.add_task(fastp(), tag=sample_tag)
            args['read1'].value = r1
            args['out1'].value = f'{sample_tag}.clean.R1.fq.gz'
            if r2 is not None:
                args['read2'].value = r2
                args['out2'].value = f'{sample_tag}.clean.R2.fq.gz'
            args['html'].value = f'{sample_tag}.fastp.html'
            args['json'].value = f'{sample_tag}.fastp.json'
            fastp_task.outputs['html'].report = True
            fastp_task.outputs['json'].report = True

            # fastq2sam
            fastq2sam_task, args = wf.add_task(FastqToSam(sample), tag=sample_tag, depends=[fastp_task])
            args['read1'].value = fastp_task.outputs['out1']
            args['read2'].value = fastp_task.outputs['out2']
            args['out'].value = f'{sample_tag}.unmapped.bam'

            # bwa alignment
            if make_index:
                depend_tasks = [index_task, fastq2sam_task]
            else:
                depend_tasks = [fastq2sam_task]
            bwa_task, args = wf.add_task(uBam2FastqBwaMem(sample_tag), tag=sample_tag, depends=depend_tasks)
            args['ubam'].value = fastq2sam_task.outputs['out']
            args['ref'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']

            # merge
            merge_bam_task, args = wf.add_task(MergeBamAlignment(sample_tag), tag=sample_tag, depends=[fastq2sam_task, bwa_task])
            merge_bam_tasks.append(merge_bam_task)
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
            args['ALIGNED_BAM'].value = bwa_task.outputs['out']
            args['UNMAPPED_BAM'].value = fastq2sam_task.outputs['out']

        # mark duplicates, 利用mark-duplicate自动把一个样本的多个fastq比对结果汇总到一起
        markdup_task, args = wf.add_task(MarkDuplicates(sample), tag=sample)
        args['INPUT'].value = [tsk.outputs['out'] for tsk in merge_bam_tasks]

        # sort and fix tag
        sort_task, args = wf.add_task(SortAndFixTags(sample), tag=sample, depends=[markdup_task])
        args['INPUT'].value = markdup_task.outputs['out']
        args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']

        # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
        bqsr_tasks = []
        for ind, each in enumerate(seq_groups):
            bsqr_task, args = wf.add_task(BaseRecalibrator(f'{sample}-g{ind}'), tag=f'{sample}-g{ind}', depends=[sort_task], parent_wkdir='BaseRecalibrator')
            bqsr_tasks.append(bsqr_task)
            args['INPUT'].value = sort_task.outputs['out']
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
            args['known-sites'].value = [top_vars['dbsnp'], top_vars['mills'], top_vars['hapmap']]
            args['intervals'].value = each

        #  Merge the recalibration reports resulting from by-interval recalibration
        merge_bsqr_task, args = wf.add_task(GatherBQSRReports(sample), tag=sample, depends=bqsr_tasks)
        args['INPUT'].value = [x.outputs['out'] for x in bqsr_tasks]
        recal_dict[sample] = merge_bsqr_task

        # apply bqsr
        apply_tasks = []
        for ind, each in enumerate(seq_groups2):
            apply_task, args = wf.add_task(ApplyBQSR(f'{sample}-g{ind}'), tag=f'{sample}-g{ind}', parent_wkdir='ApplyBQSR')
            apply_tasks.append(apply_task)
            args['INPUT'].value = sort_task.outputs['out']
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
            args['intervals'].value = each
            args['bqsr'].value = merge_bsqr_task.outputs['out']

        # Merge the recalibrated BAM files resulting from by-interval recalibration
        merge_bam_task, args = wf.add_task(GatherBamFiles(sample), tag=sample, depends=apply_tasks)
        args['INPUT'].value = [x.outputs['out'] for x in apply_tasks]
        args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
        bam_dict[sample] = merge_bam_task
        merge_bam_task.outputs['out'].report = True

    # split intervals
    split_task, args = wf.add_task(SplitIntervals(scatter_number=int(wf.args.scatter)), tag='ForCaller')
    args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
    args['intervals'].value = [wf.topvars['intervals']]
    args['outdir'].value = '.'
    interval_files = [split_task.outputs[f'out{i}'] for i in range(int(wf.args.scatter))]
    # print(interval_files)

    # 变异检测分析
    germline_merge_vcf_tasks = []
    already_exist = set()
    for tumor_sample, normal_sample in pair_list:
        if tumor_sample not in bam_dict and tumor_sample.lower() != 'none':
            print(f'Warning: skip tumor sample {tumor_sample} since it is not in target list: {list(bam_dict.keys())}')
            continue
        if normal_sample not in bam_dict and normal_sample.lower() != 'none':
            print(f'Warning: skip normal sample {normal_sample} since it is not in target list: {list(bam_dict.keys())}')
            continue

        # paired somatic variant calling
        if normal_sample.lower() != 'none' and tumor_sample.lower() != 'none':
            print(f'# performing paired somatic variant calling: {tumor_sample} vs {normal_sample}')
            mutect_tasks = []
            tumor_pileup_tasks = []
            normal_pileup_tasks = []
            for ind, interval_file in enumerate(interval_files):
                mutect_task, args = wf.add_task(Mutect2(f'{tumor_sample}-{ind}'), tag=f'{tumor_sample}-{ind}', parent_wkdir='Mutect2')
                mutect_tasks.append(mutect_task)
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
                args['tumor_bam'].value = bam_dict[tumor_sample].outputs['out']
                args['normal_bam'].value = bam_dict[normal_sample].outputs['out']
                args['tumor_name'].value = tumor_sample
                args['normal_name'].value = normal_sample
                args['germline-resource'].value = wf.topvars['germline_vcf']
                args['intervals'].value = [interval_file]
                args['pon'].value = wf.topvars['pon']
                args['alleles'].value = wf.topvars['alleles']

                # get pileup summary
                tumor_pileup_task, args = wf.add_task(GetPileupSummaries(f'{tumor_sample}-{ind}'), tag=f'{tumor_sample}-{ind}', parent_wkdir='GetPileupSummaries')
                tumor_pileup_tasks.append(tumor_pileup_task)
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
                args['intervals'].value = [interval_file]
                args['variants_for_contamination'].value = wf.topvars['contamination_vcf']
                args['bam'].value = bam_dict[tumor_sample].outputs['out']

                normal_pileup_task, args = wf.add_task(GetPileupSummaries(f'{normal_sample}-{ind}'), tag=f'{normal_sample}-{ind}', parent_wkdir='GetPileupSummaries')
                normal_pileup_tasks.append(normal_pileup_task)
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
                args['intervals'].value = [interval_file]
                args['variants_for_contamination'].value = wf.topvars['contamination_vcf']
                args['bam'].value = bam_dict[normal_sample].outputs['out']

            # LearnReadOrientationModel
            lrom_task, args = wf.add_task(LearnReadOrientationModel(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['inputs'].value = [x.outputs['f1r2'] for x in mutect_tasks]

            # merge vcf
            merge_vcf_task, args = wf.add_task(MergeVcfs(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['inputs'].value = [x.outputs['out'] for x in mutect_tasks]

            # merge stats
            merge_stat_task, args = wf.add_task(MergeMutectStats(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['stats'].value = [x.outputs['stats'] for x in mutect_tasks]

            # merge pileup summary and calculate contamination
            contaminate_task = None
            if wf.topvars['contamination_vcf'].value is not None:
                # for tumor sample
                merge_tumor_pileup_task, args = wf.add_task(GatherPileupSummaries(tumor_sample), tag=tumor_sample, depends=tumor_pileup_tasks)
                args['sequence-dictionary'].value = index_task.outputs['ref_dict'] if make_index else wf.topvars['ref'].value.rsplit('.', 1)[0] + '.dict'
                args['inputs'].value = [x.outputs['out'] for x in tumor_pileup_tasks]
                # for normal sample
                merge_normal_pileup_task, args = wf.add_task(GatherPileupSummaries(normal_sample), tag=normal_sample, depends=normal_pileup_tasks)
                args['sequence-dictionary'].value = index_task.outputs['ref_dict'] if make_index else wf.topvars['ref'].value.rsplit('.', 1)[0] + '.dict'
                args['inputs'].value = [x.outputs['out'] for x in normal_pileup_tasks]
                # calculate contamination
                contaminate_task, args = wf.add_task(CalculateContamination(tumor_sample), tag=tumor_sample)
                args['tumor_pileups'].value = merge_tumor_pileup_task.outputs['out']
                args['normal_pileups'].value = merge_normal_pileup_task.outputs['out']

            # filtering variant
            filter_task, args = wf.add_task(FilterMutectCalls(tumor_sample), tag=tumor_sample)
            args['vcf'].value = merge_vcf_task.outputs['out']
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
            if contaminate_task:
                args['contamination-table'].value = contaminate_task.outputs['out']
                args['tumor-segmentation'].value = contaminate_task.outputs['tumor-segmentation']
            args['ob-priors'].value = lrom_task.outputs['out']
            args['stats'].value = merge_stat_task.outputs['out']
            # filter_task.outputs['out'].report = True

            # filter alignment artifact
            filter_align_task = None
            if wf.topvars['bwaMemIndexImage'].value is not None:
                filter_align_task, args = wf.add_task(FilterAlignmentArtifacts(tumor_sample), tag=tumor_sample)
                args['vcf'].value = filter_task.outputs['out']
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
                args['bam'].value = bam_dict[tumor_sample].outputs['out']
                args['bwa-mem-index-image'].value = wf.topvars['bwaMemIndexImage']
                # filter_align_task.outputs['out'].report = True

            # normalize vcf需要调整到过滤之后进行，在vep注释之前完成即可
            depend_task = filter_align_task or filter_task
            norm_vcf_task, args = wf.add_task(bcftools_norm(), tag=tumor_sample, depends=[depend_task])
            args['vcf'].value = depend_task.outputs['out']
            args['fasta-ref'].value = wf.topvars['ref']
            args['out'].value = tumor_sample + '.somatic.vcf'
            norm_vcf_task.outputs['out'].report = True

            # VEP annotation
            if wf.args.vep_cache_dir and wf.args.vep_plugin_dir:
                vep_task, args = wf.add_task(vep(tumor_sample), tag=tumor_sample, depends=[norm_vcf_task])
                args['input_file'].value = norm_vcf_task.outputs['out']
                args['fasta'].value = wf.topvars['ref']
                args['dir_cache'].value = top_vars['vep_cache_dir']
                args['dir_plugins'].value = top_vars['vep_plugin_dir']
                vep_task.outputs['out_vcf'].report = True
                vep_task.outputs['out_vcf_idx'].report = True

        # tumor only analysis
        if normal_sample.lower() == 'none' and tumor_sample.lower() != 'none':
            print(f'# performing tumor-only mode analysis for sample: {tumor_sample}')
            mutect_tasks = []
            for ind, interval_file in enumerate(interval_files):
                mutect_task, args = wf.add_task(Mutect2(f'{tumor_sample}-{ind}'), tag=f'{tumor_sample}-{ind}')
                mutect_tasks.append(mutect_task)
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
                args['tumor_bam'].value = bam_dict[tumor_sample].outputs['out']
                args['tumor_name'].value = tumor_sample
                args['germline-resource'].value = wf.topvars['germline_vcf']
                args['intervals'].value = [interval_file]
                args['pon'].value = wf.topvars['pon']
                args['alleles'].value = wf.topvars['alleles']

            # LearnReadOrientationModel
            lrom_task, args = wf.add_task(LearnReadOrientationModel(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['inputs'].value = [x.outputs['f1r2'] for x in mutect_tasks]

            # merge vcf
            merge_vcf_task, args = wf.add_task(MergeVcfs(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['inputs'].value = [x.outputs['out'] for x in mutect_tasks]

            # merge stats
            merge_stat_task, args = wf.add_task(MergeMutectStats(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['stats'].value = [x.outputs['stats'] for x in mutect_tasks]

            # filter
            filter_task, args = wf.add_task(FilterMutectCalls(tumor_sample), tag=tumor_sample)
            args['vcf'].value = merge_vcf_task.outputs['out']
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
            args['ob-priors'].value = lrom_task.outputs['out']
            args['stats'].value = merge_stat_task.outputs['out']
            # filter_task.outputs['out'].report = True

            # filter alignment artifact
            filter_align_task = None
            if wf.topvars['bwaMemIndexImage'].value is not None:
                filter_align_task, args = wf.add_task(FilterAlignmentArtifacts(tumor_sample), tag=tumor_sample, depends=[filter_task])
                args['vcf'].value = filter_task.outputs['out']
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
                args['bam'].value = bam_dict[tumor_sample].outputs['out']
                args['bwa-mem-index-image'].value = wf.topvars['bwaMemIndexImage']

            # normalize vcf
            depend_task = filter_align_task or filter_task
            norm_vcf_task, args = wf.add_task(bcftools_norm(), tag=tumor_sample)
            args['vcf'].value = depend_task.outputs['out']
            args['fasta-ref'].value = wf.topvars['ref']
            args['out'].value = tumor_sample + '.somatic.vcf'
            norm_vcf_task.outputs['out'].report = True

            if wf.args.vep_cache_dir and wf.args.vep_plugin_dir:
                vep_task, args = wf.add_task(vep(tumor_sample), tag=tumor_sample)
                args['input_file'].value = norm_vcf_task.outputs['out']
                args['fasta'].value = wf.topvars['ref']
                args['dir_cache'].value = top_vars['vep_cache_dir']
                args['dir_plugins'].value = top_vars['vep_plugin_dir']
                vep_task.outputs['out_vcf'].report = True
                vep_task.outputs['out_vcf_idx'].report = True

        # germline variant calling
        if normal_sample.lower() != 'none' and (normal_sample not in already_exist) and (wf.args.call_mode in ['mix', 'germline']):
            already_exist.add(normal_sample)
            hap_tasks = []
            for ind, interval_file in enumerate(interval_files):
                hap_task, args = wf.add_task(Haplotyper(f'{normal_sample}-{ind}'), tag=f'{normal_sample}-{ind}', parent_wkdir='Haplotyper')
                hap_tasks.append(hap_task)
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
                args['bam'].value = bam_dict[normal_sample].outputs['out']
                args['intervals'].value = [interval_file]

            # merge gvcf
            merge_vcf_task, args = wf.add_task(MergeVcfs(normal_sample), tag=normal_sample, depends=hap_tasks)
            args['inputs'].value = [x.outputs['out'] for x in hap_tasks]
            args['out'].value = f'{normal_sample}.gvcf.gz'
            germline_merge_vcf_tasks.append(merge_vcf_task)
            merge_vcf_task.outputs['out'].report = True

    # perform joint calling based on germline variants
    interval_files = []
    if germline_merge_vcf_tasks:
        # split intervals, 立即运算，获得结果
        split_task = SplitIntervals(scatter_number=int(wf.args.scatter))
        split_task.task_id = None
        args = split_task.args
        args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
        args['intervals'].value = [wf.topvars['intervals']]
        args['mode'].value = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
        args['outdir'].value = 'SplitIntervals-ForJointCalling'
        split_task.run_now(wkdir=wf.args.outdir, wf_tasks=wf.tasks, docker=wf.args.docker)
        interval_files = [x.value for x in split_task.outputs.values()]
        # print(interval_files)

        # split intervals，当实际生成的interval_files数量和指定的数量不一致时就会导致错误，因此采用上面立即运行的方式
        # split_task, args = wf.add_task(SplitIntervals(scatter_number=int(wf.args.scatter)), tag='ForJoint')
        # args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
        # args['intervals'].value = [wf.topvars['intervals']]
        # args['mode'].value = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
        # args['outdir'].value = '.'
        # interval_files = [split_task.outputs[f'out{i}'] for i in range(10)]

    scattered_tasks = []
    for ind, interval_file in enumerate(interval_files):
        import_vcf_task, args = wf.add_task(GenomicsDBImport(), tag=f'{ind}', parent_wkdir='GenomicsDBImport')
        args['gvcfs'].value = [x.outputs['out'] for x in germline_merge_vcf_tasks]
        args['interval'].value = interval_file

        genotype_task, args = wf.add_task(GenotypeGVCFs(ind), tag=f'{ind}', parent_wkdir='GenotypeGVCFs')
        args['REFERENCE_SEQUENCE'].value = wf.topvars['ref']
        args['dbsnp'].value = wf.topvars['dbsnp']
        args['gendb'].value = import_vcf_task.outputs['out']
        args['interval'].value = interval_file

        hard_filter_task, args = wf.add_task(VariantFiltration(f'scatter{ind}'), tag=f'{ind}', parent_wkdir='VariantHardFiltration')
        args['vcf'].value = genotype_task.outputs['out']

        make_site_only_task, args = wf.add_task(MakeSitesOnlyVcf(f'scatter{ind}'), tag=f'{ind}', parent_wkdir='MakeSitesOnlyVcf')
        args['vcf'].value = hard_filter_task.outputs['out']
        scattered_tasks.append(make_site_only_task)

    if scattered_tasks:
        gather_vcf_task, args = wf.add_task(MergeVcfs('Joint.SiteOnly'), tag='Joint', depends=scattered_tasks)
        args['inputs'].value = [x.outputs['out'] for x in scattered_tasks]

        indel_recal_task, args = wf.add_task(IndelsVariantRecalibrator('Joint'), tag='Joint', depends=[gather_vcf_task])
        args['vcf'].value = gather_vcf_task.outputs['out']
        args['mills'].value = wf.topvars['mills']
        args['axiomPloly'].value = wf.topvars['axiomPoly']
        args['dbsnp'].value = wf.topvars['dbsnp']

        snp_recal_task, args = wf.add_task(SNPsVariantRecalibrator('Joint'), tag='Joint', depends=[gather_vcf_task])
        args['vcf'].value = gather_vcf_task.outputs['out']
        args['hapmap'].value = wf.topvars['hapmap']
        args['omni'].value = wf.topvars['omni']
        args['dbsnp'].value = wf.topvars['dbsnp']
        args['1000G'].value = wf.topvars['G1000']

        final_recal_tasks = []
        for idx, task in enumerate(scattered_tasks):
            apply_indel_task, args = wf.add_task(ApplyVQSR(f'joint.{idx}.indel'), tag=f'Indel{idx}', parent_wkdir='ApplyIndelRecalibration')
            args['vcf'].value = task.outputs['out']
            args['recal-file'].value = indel_recal_task.outputs['out']
            args['tranches-file'].value = indel_recal_task.outputs['tranches']
            args['mode'].value = 'INDEL'

            apply_snp_task, args = wf.add_task(ApplyVQSR(f'joint.{idx}.snp'), tag=f'SNP{idx}', parent_wkdir='ApplySnpRecalibration')
            args['vcf'].value = apply_indel_task.outputs['out']
            args['recal-file'].value = snp_recal_task.outputs['out']
            args['tranches-file'].value = snp_recal_task.outputs['tranches']
            args['mode'].value = 'SNP'
            final_recal_tasks.append(apply_snp_task)

        if final_recal_tasks:
            gather_final_vcf_task, args = wf.add_task(MergeVcfs('Joint.final'), tag='Final', depends=final_recal_tasks)
            args['inputs'].value = [x.outputs['out'] for x in final_recal_tasks]

            # normalize vcf
            norm_vcf_task, args = wf.add_task(bcftools_norm(), tag='JointGermline', depends=[gather_final_vcf_task])
            args['vcf'].value = gather_final_vcf_task.outputs['out']
            args['fasta-ref'].value = wf.topvars['ref']
            args['out'].value = 'Germline.Joint.LeftNormalized.vcf'
            norm_vcf_task.outputs['out'].report = True

            if wf.args.vep_cache_dir and wf.args.vep_plugin_dir:
                depend_task = norm_vcf_task
                vep_task, args = wf.add_task(vep('Joint'), tag='JointGermline', depends=[depend_task])
                args['input_file'].value = depend_task.outputs['out']
                args['fasta'].value = wf.topvars['ref']
                args['dir_cache'].value = top_vars['vep_cache_dir']
                args['dir_plugins'].value = top_vars['vep_plugin_dir']
                vep_task.outputs['out_vcf'].report = True
                vep_task.outputs['out_vcf_idx'].report = True

    wf.run()


if __name__ == '__main__':
    pipeline()

