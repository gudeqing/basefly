import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, Task
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
    # 定义输出,这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out1'] = Output(value="{out1}")
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd


def build_index():
    cmd = Command()
    cmd.meta.name = 'buildIndex'
    cmd.meta.desc = 'bwa index and create sequence dictionary and fasta fai file'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.tool = '/usr/gitc/bwa index'
    cmd.args['ref_fasta'] = Argument(prefix='ln -s ', type='infile', desc='reference fasta to index')
    cmd.args['_new_name'] = Argument(type='fix', value='ref_genome.fa && samtools faidx ref_genome.fa &&')
    cmd.args['_tool'] = Argument(type='fix', value='/usr/gitc/bwa index')
    cmd.args['block_size'] = Argument(prefix='-b ', default=10000000, desc='block size for the bwtsw algorithm (effective with -a bwtsw)')
    cmd.args['name64'] = Argument(prefix='-6', type='bool', default=False, desc='index files named as <in.fasta>.64.* instead of <in.fasta>.*')
    cmd.args['_ref_fasta'] = Argument(type='fix', value='ref_genome.fa', desc='reference fasta to be indexed')
    cmd.args['_create_dict'] = Argument(type='fix', value='| gatk CreateSequenceDictionary')
    cmd.args['_reference'] = Argument(prefix='R=', type='fix', value='ref_genome.fa', desc='reference fasta')
    cmd.args['ref_dict'] = Argument(prefix='O=', default='ref_genome.dict', desc="Output a dict file containing only the sequence dictionary. By default it will use the base name of the input reference with the .dict extension")
    cmd.args['alt_names'] = Argument(prefix='ALT_NAMES=', type='infile', level='optional', desc="Optional file containing the alternative names for the contigs. Tools may use this information to consider different contig notations as identical (e.g: 'chr1' and '1'). The alternative names will be put into the appropriate @AN annotation for each contig. No header. First column is the original name, the second column is an alternative name. One contig may have more than one alternative name. Default value: null.")
    cmd.outputs['ref_genome'] = Output(value='ref_genome.fa')
    cmd.outputs['ref_dict'] = Output(value='{ref_dict}')
    cmd.outputs['index_dir'] = Output(value='.')
    return cmd


def FastqToSam(sample):
    cmd = Command()
    cmd.meta.name = 'FastqToSam'
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


def uBam2FastqBwaMem(sample):
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


def MergeBamAlignment(sample):
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
    cmd.args['PROGRAM_GROUP_VERSION'] = Argument(prefix='--PROGRAM_GROUP_VERSION ', default='0.7.?')
    cmd.args['PROGRAM_GROUP_COMMAND_LINE'] = Argument(prefix='--PROGRAM_GROUP_COMMAND_LINE ', default='')
    cmd.args['PROGRAM_GROUP_NAME'] = Argument(prefix='--PROGRAM_GROUP_NAME ', default='"bwamem"')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def MarkDuplicates(sample):
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


def SortAndFixTags(sample):
    cmd = Command()
    cmd.meta.name = 'SortAndFixTags'
    cmd.meta.desc = 'sort bam and fix tags'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['_sort_sam'] = Argument(type='fix', value='java -jar picard.jar SortSam')
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', desc='input bam file list')
    cmd.args['_OUTPUT'] = Argument(prefix='--OUTPUT ', default='/dev/stdout', desc='output bam file')
    cmd.args['SORT_ORDER'] = Argument(prefix='--SORT_ORDER ', default='"coordinate"')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='false')
    cmd.args['_fix_tags'] = Argument(type='fix', value='java -jar picard.jar SetNmMdAndUqTags --CREATE_INDEX true --INPUT /dev/stdin')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.sorted.dup_marked.bam', desc='output bam file')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='--REFERENCE_SEQUENCE ', type='infile', desc='reference fasta file')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    cmd.outputs['out_idx'] = Output(value='{OUTPUT}.bai')
    return cmd


def CreateSequenceGroupingTSV(ref_dict):
    print('create the Sequencing Groupings used for BQSR and PrintReads Scatter.')
    with open(ref_dict, "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
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

    with open("sequence_grouping_with_unmapped.txt", "w") as tsv_file_with_unmapped:
        tsv_file_with_unmapped.write(tsv_string)
        tsv_file_with_unmapped.close()

    sequence_grouping = [x.split('\t') for x in tsv_string.strip().split('\n')]
    sequence_grouping_with_unmapped = sequence_grouping.copy() + [['unmapped']]
    return sequence_grouping, sequence_grouping_with_unmapped


def BaseRecalibrator(sample):
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
    cmd.args['intervals'] = Argument(prefix='--intervals ', type='str', array=True, multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def GatherBQSRReports(sample):
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


def ApplyBQSR(sample):
    cmd = Command()
    cmd.meta.name = 'ApplyBQSR'
    cmd.meta.desc = 'Apply Base Quality Score Recalibration (BQSR) model'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk ApplyBQSR'
    cmd.args['INPUT'] = Argument(prefix='-I ', type='infile', desc='input bam file')
    cmd.args['OUTPUT'] = Argument(prefix='-O ', default=f'{sample}.recalibrated.bam', desc='The output recalibration table file to create')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['intervals'] = Argument(prefix='-L ', type='str', array=True, multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['bqsr'] = Argument(prefix='-bqsr ', type='infile', desc='Input recalibration table for BQSR')
    cmd.args['static-quantized-quals'] = Argument(prefix='--static-quantized-quals ', type='int', array=True, multi_times=True, default=[10, 20, 30])
    cmd.args['add-output-sam-program-record'] = Argument(prefix='--add-output-sam-program-record', type='bool', default=True)
    cmd.args['use-original-qualities'] = Argument(prefix='--use-original-qualities', type='bool', default=True)
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    return cmd


def GatherBamFiles(sample):
    cmd = Command()
    cmd.meta.name = 'GatherBamFiles'
    cmd.meta.desc = 'Concatenate efficiently BAM files that resulted from a scattered parallel analysis.'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'gatk GatherBamFiles'
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', array=True, multi_times=True, desc='input bam file')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', default=f'{sample}.recalibrated.bam', desc='output bam file')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', level='optional', desc='reference fasta file')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='true')
    cmd.outputs['out'] = Output(value='{OUTPUT}')
    cmd.outputs['out_idx'] = Output(value='{OUTPUT}.bai')
    return cmd


def SplitIntervals():
    cmd = Command()
    cmd.meta.name = 'SplitIntervals'
    cmd.meta.desc = 'This tool takes in intervals via the standard arguments of IntervalArgumentCollection and splits them into interval files for scattering. The resulting files contain equal number of bases.'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 3 * 1024 ** 3
    cmd.runtime.tool = 'mkdir intervals_folder && gatk SplitIntervals'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['intervals'] = Argument(prefix='-L ', type='infile', array=True, multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['scatter'] = Argument(prefix='-scatter ', default=10, desc='number of output interval files to split into')
    cmd.args['outdir'] = Argument(prefix='-O ', default='intervals_folder', desc='The directory into which to write the scattered interval sub-directories.')
    cmd.outputs['out'] = Output(value='{outdir}', type='outdir')
    return cmd


def Mutect2(prefix):
    cmd = Command()
    cmd.meta.name = 'Mutect2'
    cmd.meta.desc = 'Call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations.'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2'
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5*1024**3
    cmd.runtime.cpu = 3
    cmd.runtime.tool = 'gatk Mutect2'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['tumor_bam'] = Argument(prefix='-I ', type='infile', desc='tumor bam')
    cmd.args['normal_bam'] = Argument(prefix='-I ', type='infile', level='optional', desc='normal bam')
    cmd.args['tumor_name'] = Argument(prefix='-tumor ', desc='tumor sample name')
    cmd.args['normal_name'] = Argument(prefix='-normal ', level='optional', desc='normal sample name')
    cmd.args['germline-resource'] = Argument(prefix='--germline-resource ', type='infile', level='optional', desc='optional database of known germline variants (and its index) (see http://gnomad.broadinstitute.org/downloads)')
    cmd.args['pon'] = Argument(prefix='-pon ', type='infile', level='optional', desc='')
    cmd.args['intervals'] = Argument(prefix='-L ', type='infile', level='optional', array=True, multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['string_intervals'] = Argument(prefix='-L ', level='optional', desc='interval string such as "chrM"')
    cmd.args['alleles'] = Argument(prefix='--alleles ', type='infile', level='optional', desc='The set of alleles to force-call regardless of evidence')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.vcf.gz', desc='output vcf')
    cmd.args['bam_output'] = Argument(prefix='--bam-output ', level='optional', desc='output bam file')
    cmd.args['f1r2-tar-gz'] = Argument(prefix='--f1r2-tar-gz ', type='infile', default=f'{prefix}.f1r2.tar.gz', desc='If specified, collect F1R2 counts and output files into this tar.gz file')
    cmd.args['mitochondria'] = Argument(prefix='--mitochondira', type='bool', desc='if to turn on mitochondria mode. Specifically, the mode sets --initial-tumor-lod to 0, --tumor-lod-to-emit to 0, --af-of-alleles-not-in-resource to 4e-3, and the advanced parameter --pruning-lod-thres')
    cmd.outputs['out'] = Output(value='{out}')
    cmd.outputs['f1r2'] = Output(value='{f1r2-tar-gz}')
    cmd.outputs['stats'] = Output(value='{out}.stats')
    return cmd


def GetPileupSummaries(prefix):
    cmd = Command()
    cmd.meta.name = 'GetPileupSummaries'
    cmd.meta.desc = 'Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360037593451-GetPileupSummaries'
    cmd.runtime.tool = 'gatk GetPileupSummaries'
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='BAM/SAM/CRAM file containing reads')
    cmd.args['interval-set-rule'] = Argument(prefix='--interval-set-rule ', default='INTERSECTION', desc='Set merging approach to use for combining interval inputs')
    cmd.args['intervals'] = Argument(prefix='-L ', type='infile', array=True, multi_times=True, desc='One or more genomic intervals over which to operate')
    cmd.args['variants_for_contamination'] = Argument(prefix='-V ', type='infile', desc='A VCF file containing variants and allele frequencies')
    cmd.args['out'] = Argument(prefix='-O ', type='infile', default=f'{prefix}.tumor-pileups.table', desc='The output table')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def MergeVcfs(sample):
    cmd = Command()
    cmd.meta.name = 'MergeVcfs'
    cmd.meta.desc = 'Combines multiple variant files into a single variant file.'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360056969852-MergeVcfs-Picard-'
    cmd.runtime.tool = 'gatk MergeVcfs'
    cmd.args['inputs'] = Argument(prefix='-I ', type='infile', array=True, multi_times=True, desc='input vcf list')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.vcf.gz', desc='The merged VCF or BCF file. File format is determined by file extension.')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def SortBam():
    cmd = Command()
    cmd.meta.name = 'SortBam'
    cmd.meta.desc = 'This tool sorts the input SAM or BAM file by coordinate, queryname (QNAME), or some other property of the SAM record. The SortOrder of a SAM/BAM file is found in the SAM file header tag @HD in the field labeled SO.'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/360036366932-SortSam-Picard-'
    cmd.runtime.tool = 'gatk SortSam'
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='BAM/SAM/CRAM file containing reads')
    cmd.args['out'] = Argument(prefix='-O ', desc='Sorted BAM or SAM output file.')
    cmd.args['SORT_ORDER'] = Argument(prefix='--SORT_ORDER ', default='coordinate', desc='Sort order of output file.')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='true')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def MergeMutectStats(sample):
    cmd = Command()
    cmd.meta.name = 'MergeMutectStats'
    cmd.runtime.tool = 'gatk MergeMutectStats'
    cmd.args['stats'] = Argument(prefix='-stats ', type='infile', multi_times=True, array=True)
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.vcf.stats', desc='output merged stat files')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def GatherPileupSummaries(sample):
    cmd = Command()
    cmd.meta.name = 'GatherPileupSummaries'
    cmd.runtime.tool = 'gatk GatherPileupSummaries'
    cmd.args['sequence-dictionary'] = Argument(prefix='--sequence-dictionary ', type='infile')
    cmd.args['inputs'] = Argument(prefix='-I ', type='infile', multi_times=True, array=True)
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.pileup.table', desc='output merged pileup summary file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def LearnReadOrientationModel(sample):
    cmd = Command()
    cmd.meta.name = 'LearnReadOrientationModel'
    cmd.meta.desc = 'Learn the prior probability of read orientation artifact from the output of CollectF1R2Counts of Mutect2'
    cmd.runtime.tool = 'gatk LearnReadOrientationModel'
    cmd.args['inputs'] = Argument(prefix='-I ', type='infile', multi_times=True, array=True, desc='One or more .tar.gz containing outputs of CollectF1R2Counts')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.artifact-priors.tar.gz', desc='tar.gz of artifact prior tables')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def CalculateContamination(prefix):
    cmd = Command()
    cmd.meta.name = 'CalculateContamination'
    cmd.meta.desc = 'Calculates the fraction of reads coming from cross-sample contamination, given results from GetPileupSummaries. The resulting contamination table is used with FilterMutectCalls.'
    cmd.runtime.tool = 'gatk CalculateContamination'
    cmd.args['tumor_pileups'] = Argument(prefix='-I ', type='infile', desc='input pileup table')
    cmd.args['normal_pileups'] = Argument(prefix='-matched ', type='infile', desc='The matched normal input table')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{prefix}.contamination.table', desc='output table')
    cmd.args['tumor-segmentation'] = Argument(prefix='--tumor-segmentation ', default=f'{prefix}.tumor_segmentation', desc='The output table containing segmentation of the tumor by minor allele fraction')
    cmd.outputs['out'] = Output(value='{out}')
    cmd.outputs['tumor-segmentation'] = Output(value='{tumor-segmentation}')
    return cmd


def FilterMutectCalls(sample):
    cmd = Command()
    cmd.meta.name = 'FilterMutectCalls'
    cmd.meta.desc = 'FilterMutectCalls applies filters to the raw output of Mutect2'
    cmd.runtime.tool = 'gatk FilterMutectCalls'
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
    cmd.meta.name = 'FilterAlignmentArtifacts'
    cmd.meta.desc = 'Alignment artifacts can occur whenever there is sufficient sequence similarity between two or more regions in the genome to confuse the alignment algorithm. This can occur when the aligner for whatever reason overestimate how uniquely a read maps, thereby assigning it too high of a mapping quality. It can also occur through no fault of the aligner due to gaps in the reference, which can also hide the true position to which a read should map. By using a good alignment algorithm (the GATK wrapper of BWA-MEM), giving it sensitive settings (which may have been impractically slow for the original bam alignment) and mapping to the best available reference we can avoid these pitfalls. The last point is especially important: one can (and should) use a BWA-MEM index image corresponding to the best reference, regardless of the reference to which the bam was aligned.'
    cmd.meta.source = 'https://gatk.broadinstitute.org/hc/en-us/articles/4418051467035-FilterAlignmentArtifacts-EXPERIMENTAL-'
    cmd.runtime.tool = 'gatk FilterAlignmentArtifacts'
    cmd.args['vcf'] = Argument(prefix='-V ', type='infile', desc='A VCF file containing variants')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='input bam file')
    cmd.args['bwa-mem-index-image'] = Argument(prefix='--bwa-mem-index-image ', type='infile', desc='BWA-mem index image')
    cmd.args['out'] = Argument(prefix='-O ', default=f'{sample}.align_artifacts_filtered.vcf.gz', desc='output vcf file')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def vep(sample):
    cmd = Command()
    cmd.meta.name = 'VEP'
    cmd.runtime.image = 'ensemblorg/ensembl-vep:2.0.3'
    cmd.runtime.tool = 'vep'
    cmd.args['input_file'] = Argument(prefix='-i ', type='infile', desc='input file')
    cmd.args['fasta'] = Argument(prefix='--fasta ', type='infile', desc="Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run VEP with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache).")
    cmd.args['output_file'] = Argument(prefix='-o ', default=f'{sample}.vep.vcf.gz', desc='output file')
    cmd.args['output_format'] = Argument(prefix='--', range={'vcf', 'json', 'tab'}, default='vcf', desc="If we choose to write output in VCF format. Consequences are added in the INFO field of the VCF file, using the key 'CSQ'. Data fields are encoded separated by '|'; the order of fields is written in the VCF header. Output fields in the 'CSQ' INFO field can be selected by using --fields.")
    cmd.args['compress_output'] = Argument(prefix='--compress_output ', default='bgzip', desc="Writes output compressed using either gzip or bgzip")
    cmd.args['force_overwrite'] = Argument(prefix="--force_overwrite ", type='bool', default=True, desc="Force overwriting of output file")
    cmd.args['fork'] = Argument(prefix='--fork ', type='int', default=7, desc='Use forking(multi-cpu/threads) to improve script runtime')
    cmd.args['species'] = Argument(prefix='--species ', default='homo_sapiens', desc='Species for your data. This can be the latin name e.g. homo_sapiens or any Ensembl alias e.g. mouse.')
    cmd.args['assembly_version'] = Argument(prefix='--assembly ', default='GRCh37', desc='Select the assembly version to use if more than one available.')
    cmd.args['dir_cache'] = Argument(prefix='--dir_cache ', type='indir', desc='Specify the cache directory to use')
    cmd.args['dir_plugins'] = Argument(prefix='--dir_plugins ', type='indir', desc='Specify the plugin directory to use')
    cmd.args['stats_file'] = Argument(prefix='--stats_file ', default=f'{sample}.vep.summary.html', desc='Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>.')
    cmd.args['cache'] = Argument(prefix='--cache ', type='bool', default=True, desc='Enables use of cache')
    cmd.args['offline'] = Argument(prefix='--offline ', type='bool', default=True, desc='Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation')
    cmd.args['merged'] = Argument(prefix='--merged ', type='bool', default=False, desc='Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.')
    cmd.args['plugins'] = Argument(prefix='--plugin ', multi_times=True, default=['Frameshift', 'Wildtype'], desc='Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times')
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


def Haplotyper(sample):
    cmd = Command()
    cmd.meta.name = 'Haplotyper'
    cmd.runtime.image = 'registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02'
    cmd.runtime.tool = 'gatk HaplotypeCaller'
    cmd.args['intervals'] = Argument(prefix='-L ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='reccaled tumor and normal bam list')
    cmd.args['REFERENCE_SEQUENCE'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['emit_mode'] = Argument(prefix='-ERC ', default='GVCF', range=['NONE', 'GVCF', 'BP_RESOLUTION'], desc='The reference confidence mode makes it possible to emit a per-bp or summarized confidence estimate for a site being strictly homozygous-reference.')
    cmd.args['ploidy'] = Argument(prefix='-ploidy ', type='int', default=2, desc='determines the ploidy number of the sample being processed. The default value is 2.')
    cmd.args['annotation-group'] = Argument(prefix='-G ', array=True, multi_times=True, default=['StandardAnnotation', 'StandardHCAnnotation', 'AS_StandardAnnotation'])
    cmd.args['gvcf-gq-bands'] = Argument(prefix='-GQB ', array=True, multi_times=True, default=[10, 20, 30, 40, 50, 60, 70, 80, 90], desc='Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)')
    cmd.args['out_vcf'] = Argument(prefix='', default=f'{sample}.g.vcf.gz', desc='output vcf file')
    cmd.outputs['out'] = Output(value='{out_vcf}')
    cmd.outputs['out_idx'] = Output(value='{out_vcf}.tbi')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'GATK-DNAseq-Workflow'
    wf.meta.source = ""
    wf.meta.desc = """
    当前流程是参考博得研究所最新的GATK-BEST-Practice流程构建的DNAseq突变检测分析流程改写而成
    流程包含的主要功能：
    * 使用fastp进行接头自动去除
    * 使用BWA进行比对分析
    * 使用GATK检测small SNP/Indel, 支持tumor-only和tumor-normal配对模式
    * 基于VEP进行突变注释
    * 使用hisatGenotype进行HLA基因定型分析
    """
    wf.meta.version = "1.0"
    # 定义流程输入参数
    wf.init_argparser()
    wf.add_argument('-fastq_info', nargs='+', required=True,
                    help='A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]')
    wf.add_argument('-r1_name', default='(.*).R1.fastq',
                    help="python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'")
    wf.add_argument('-r2_name', default='(.*).R2.fastq',
                    help="python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'")
    wf.add_argument('-exclude_samples', default=tuple(), nargs='+', help='samples to exclude from analysis')
    wf.add_argument('-pair_info', required=True, help='tumor normal pair info, two-column txt file, first column is tumor sample name. sample not in pair info will be skipped')
    wf.add_argument('-ref', required=True, help='reference fasta file, require bwa index being created')
    wf.add_argument('-dbsnp', required=True, help='dbsnp vcf file')
    wf.add_argument('-known_indels', required=True, help='high confidence known indel vcf file')
    wf.add_argument('-known_mills', required=True, help='high confidence known indel vcf file')
    wf.add_argument('-pon', required=False, help='panel of normal vcf file for germline variant filtering, this will be required for tumor only analysis')
    wf.add_argument('-germline_vcf', required=False, help='germline vcf, will be used for germline variant filtering and contamination analysis')
    wf.add_argument('-alleles', required=False, help='The set of alleles to force-call regardless of evidence')
    wf.add_argument('-contamination_vcf', required=False, help='germline vcf such as small_exac_common_3_b37.vcf, will be used for contamination analysis')
    wf.add_argument('-bwaMemIndexImage', required=False, help='bwa-mem-index-mage for artifact alignment filtering')
    wf.add_argument('-vep_cache_dir', required=False, help='VEP cache directory')
    wf.add_argument('-vep_plugin_dir', required=False, help='VEP plugin directory')
    wf.add_argument('-intervals', required=False,
                    help="interval file, support bed file or picard interval or vcf format.")
    wf.add_argument('-hisatgenotype_db', required=False, help='indicies dir of hisat-genotype for HLA typing')
    wf.parse_args()

    top_vars = dict(
        ref=TopVar(value=wf.args.ref, type='infile'),
        known_dbsnp=TopVar(value=wf.args.dbsnp, type='infile'),
        known_indels=TopVar(value=wf.args.known_indels, type='infile'),
        known_mills=TopVar(value=wf.args.known_mills, type='infile'),
        pon=TopVar(value=wf.args.pon, type='infile'),
        germline_vcf=TopVar(value=wf.args.germline_vcf, type='infile'),
        contamination_vcf=TopVar(value=wf.args.contamination_vcf, type='infile'),
        bwaMemIndexImage=TopVar(value=wf.args.bwaMemIndexImage, type='infile'),
        vep_cache_dir=TopVar(value=wf.args.vep_cache_dir, type='indir'),
        vep_plugin_dir=TopVar(value=wf.args.vep_plugin_dir, type='indir'),
        intervals=TopVar(value=wf.args.intervals, type='infile'),
        alleles=TopVar(value=wf.args.alleles, type='infile'),
        hisatgenotype_db=TopVar(value=wf.args.hisatgenotype_db, type='infile')
    )
    wf.add_topvars(top_vars)

    # 提取fastq信息
    fastq_info = get_fastq_info(fastq_info=wf.args.fastq_info, r1_name=wf.args.r1_name, r2_name=wf.args.r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found !')

    # 提取配对信息
    pair_list = []
    sample_list = []
    if wf.args.pair_info:
        with open(wf.args.pair_info) as f:
            for line in f:
                if line.strip():
                    pairs = line.strip('\n').split('\t')[:2]
                    pair_list.append(pairs)
                    sample_list.extend(pairs)

    # CreateSequenceGroupingTSV
    seq_groups, seq_groups2 = CreateSequenceGroupingTSV(wf.args.ref)

    # 建bwa索引
    make_index = False
    if not os.path.exists(wf.topvars['ref'].value + '.bwt'):
        make_index = True
        index_task, args = wf.add_task(build_index(), name='buildIndex')
        args['ref_fasta'].value = wf.topvars['ref']

    recal_dict = dict()
    bam_dict = dict()
    for sample, reads in fastq_info.items():
        # 跳过不需要分析的样本
        if sample in wf.args.exclude_samples or (sample not in sample_list):
            continue
        # fastq预处理和比对 考虑一个样本存在多对fastq的处理
        if sample in wf.args.exclude_samples:
            continue
        if len(reads) == 2:
            r1s, r2s = reads
        else:
            r1s = reads[0]
            r2s = [None]*len(r1s)

        for ind, (r1, r2) in enumerate(zip(r1s, r2s)):
            fastp_task, args = wf.add_task(fastp(), tag=f'{sample}-{ind}')
            args['read1'].value = r1
            args['out1'].value = f'{sample}-{ind}.clean.R1.fq.gz'
            if r2 is not None:
                args['read2'].value = r2
                args['out2'].value = f'{sample}-{ind}.clean.R2.fq.gz'
            args['html'].value = f'{sample}-{ind}.fastp.html'
            args['json'].value = f'{sample}-{ind}.fastp.json'

            # fastq2sam
            fastq2sam_task, args = wf.add_task(FastqToSam(sample), tag=f'{sample}-{ind}', depends=[fastp_task])
            args['read1'].value = fastp_task.outputs['out1']
            args['read2'].value = fastp_task.outputs['out2']
            args['out'].value = f'{sample}-{ind}.unmapped.bam'

            # bwa alignment
            if make_index:
                depend_tasks = [index_task, fastq2sam_task]
            else:
                depend_tasks = [fastq2sam_task]
            bwa_task, args = wf.add_task(uBam2FastqBwaMem(f'{sample}-{ind}'), tag=f'{sample}-{ind}', depends=depend_tasks)
            args['ubam'].value = fastq2sam_task.outputs['out']
            args['ref'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']

            # merge
            merge_bam_task, args = wf.add_task(MergeBamAlignment(f'{sample}-{ind}'), tag=f'{sample}-{ind}', depends=[fastq2sam_task, bwa_task])
            args['ALIGNED_BAM'].value = bwa_task.outputs['out']
            args['UNMAPPED_BAM'].value = fastq2sam_task.outputs['out']

        # mark duplicates, 利用mark-duplicate自动把一个样本的多个fastq比对结果汇总到一起
        depend_task_ids = [task_id for task_id in wf.tasks if wf.tasks[task_id].name.startswith(f'{merge_bam_task.cmd.meta.name}-{sample}')]
        markdup_task, args = wf.add_task(MarkDuplicates(sample), tag=sample, depends=depend_task_ids)
        args['INPUT'].value = [wf.tasks[task_id].outputs['out'] for task_id in depend_task_ids]

        # sort and fix tag
        sort_task, args = wf.add_task(SortAndFixTags(sample), tag=sample, depends=[markdup_task])
        args['INPUT'].value = markdup_task.outputs['out']
        args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if os.path.exists(wf.topvars['ref'].value + '.bwt') else index_task.outputs['ref_genome']

        # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
        bqsr_tasks = []
        for ind, each in enumerate(seq_groups):
            bsqr_task, args = wf.add_task(BaseRecalibrator(f'{sample}-{ind}'), tag=f'{sample}-{ind}', depends=[sort_task])
            bqsr_tasks.append(bsqr_task)
            args['INPUT'].value = sort_task.outputs['out']
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
            args['known_sites'].value = [top_vars['known_dbsnp'], top_vars['known_mills'], top_vars['known_indels']]
            args['intervals'].value = each

        #  Merge the recalibration reports resulting from by-interval recalibration
        merge_bsqr_task, args = wf.add_task(GatherBQSRReports(sample), tag=sample, depends=bqsr_tasks)
        args['INPUT'].value = [x.outputs['out'] for x in bqsr_tasks]
        recal_dict[sample] = merge_bsqr_task

        # apply bqsr
        apply_tasks = []
        for ind, each in enumerate(seq_groups2):
            apply_task, args = wf.add_task(ApplyBQSR(f'{sample}-{ind}'), tag=f'{sample}-{ind}', depends=[merge_bsqr_task])
            apply_tasks.append(apply_task)
            args['INPUT'].value = sort_task.outputs['out']
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
            args['intervals'].value = each
            args['bqsr'].value = merge_bsqr_task.outputs['out']

        # Merge the recalibrated BAM files resulting from by-interval recalibration
        merge_bam_task, args = wf.add_task(GatherBamFiles(sample), tag=sample, depends=apply_tasks)
        args['INPUT'].value = [x.outputs['out'] for x in apply_tasks]
        args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
        bam_dict[sample] = merge_bam_task

    # split intervals, 立即运算，获得结果
    split_task = SplitIntervals()
    split_task.args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'].value
    split_task.args['intervals'].value = [wf.topvars['intervals'].value]
    split_task.run_now(wkdir=wf.args.outdir)
    interval_files = os.listdir(os.path.join(wf.args.outdir, split_task.args['outdir'].value))

    # 变异分析
    for tumor_sample, normal_sample in pair_list:
        if tumor_sample not in bam_dict and tumor_sample.lower() != 'none':
            print(f'Warning: skip tumor sample {tumor_sample} since it is not in target list: {list(bam_dict.keys())}')
            continue
        if normal_sample not in bam_dict and normal_sample.lower() != 'none':
            print(f'Warning: skip normal sample {normal_sample} since it is not in target list: {list(bam_dict.keys())}')
            continue

        # somatic variant calling
        if normal_sample.lower() != 'none' and tumor_sample.lower() != 'none':
            mutect_tasks = []
            tumor_pileup_tasks = []
            normal_pileup_tasks = []
            for ind, interval_file in enumerate(interval_files):
                mutect_task, args = wf.add_task(Mutect2(f'{tumor_sample}-{ind}'), tag=f'{tumor_sample}-{ind}')
                mutect_tasks.append(mutect_task)
                mutect_task.depends = [bam_dict[normal_sample].task_id, bam_dict[tumor_sample].task_id]
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['tumor_bam'].value = bam_dict[tumor_sample].outputs['out']
                args['normal_bam'].value = bam_dict[normal_sample].outputs['out']
                args['tumor_name'].value = tumor_sample
                args['normal_name'].value = normal_sample
                args['germline-resource'].value = wf.topvars['germline_vcf']
                args['intervals'].value = interval_file
                args['pon'].value = wf.topvars['pon']
                args['alleles'].value = wf.topvars['alleles']

                # get pileup summary
                tumor_pileup_task, args = wf.add_task(GetPileupSummaries(f'{tumor_sample}-{ind}'), tag=f'{tumor_sample}-{ind}', depends=[bam_dict[tumor_sample]])
                tumor_pileup_tasks.append(tumor_pileup_task)
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['intervals'].value = interval_file
                args['variants_for_contamination'].value = wf.topvars['contamination_vcf']
                args['bam'].value = bam_dict[tumor_sample].outputs['out']

                normal_pileup_task, args = wf.add_task(GetPileupSummaries(f'{normal_sample}-{ind}'), tag=f'{normal_sample}-{ind}', depends=[bam_dict[normal_sample]])
                normal_pileup_tasks.append(normal_pileup_task)
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['intervals'].value = interval_file
                args['variants_for_contamination'].value = wf.topvars['contamination_vcf']
                args['bam'].value = bam_dict[normal_sample].outputs['out']

            # LearnReadOrientationModel
            lrom_task, args = wf.add_task(LearnReadOrientationModel(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['inputs'].value = [x.outputs['f1r2-tar-gz'] for x in mutect_tasks]

            # merge vcf
            merge_vcf_task, args = wf.add_task(MergeVcfs(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['inputs'].value = [x.outputs['out'] for x in mutect_tasks]

            # merge stats
            merge_stat_task, args = wf.add_task(MergeMutectStats(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['stats'].value = [x.outputs['stats'] for x in mutect_tasks]

            # merge pileup summary and calculate contamination
            if wf.topvars['contamination_vcf'].value is not None:
                # for tumor sample
                merge_tumor_pileup_task, args = wf.add_task(GatherPileupSummaries(tumor_sample), tag=tumor_sample, depends=tumor_pileup_tasks)
                args['sequence-dictionary'].value = index_task.outputs['ref_dict'] if make_index else wf.topvars['ref'].rsplit('.', 1)[0] + '.dict'
                args['inputs'].value = [x.outputs['out'] for x in tumor_pileup_tasks]
                # for normal sample
                merge_normal_pileup_task, args = wf.add_task(GatherPileupSummaries(normal_sample), tag=normal_sample, depends=normal_pileup_tasks)
                args['sequence-dictionary'].value = index_task.outputs['ref_dict'] if make_index else wf.topvars['ref'].rsplit('.', 1)[0] + '.dict'
                args['inputs'].value = [x.outputs['out'] for x in normal_pileup_tasks]
                # calculate contamination
                contaminate_task, args = wf.add_task(CalculateContamination(tumor_sample), tag=tumor_sample, depends=[merge_normal_pileup_task, merge_tumor_pileup_task])
                args['tumor_pileups'].value = merge_tumor_pileup_task.outputs['out']
                args['normal_pileups'].value = merge_normal_pileup_task.outputs['out']

            # filtering variant
            depend_tasks = [merge_vcf_task, merge_stat_task, lrom_task]
            if wf.topvars['contamination_vcf'].value is not None:
                depend_tasks.append(contaminate_task)
            filter_task, args = wf.add_task(FilterMutectCalls(tumor_sample), tag=tumor_sample, depends=depend_tasks)
            args['vcf'].value = merge_vcf_task.outputs['out']
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
            if wf.topvars['contamination_vcf'].value is not None:
                args['contamination-table'].value = contaminate_task.outputs['out']
                args['tumor-segmentation'].value = contaminate_task.outputs['tumor-segmentation']
            args['ob-priors'].value = lrom_task.outputs['out']
            args['stats'].value = merge_stat_task.outputs['out']

            # filter alignment artifact
            if wf.topvars['bwaMemIndexImage'].value is not None:
                filter_align_task, args = wf.add_task(FilterAlignmentArtifacts(tumor_sample), tag=tumor_sample, depends=[filter_task])
                args['vcf'].value = filter_task.outputs['out']
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['bam'].value = bam_dict[tumor_sample].outputs['out']
                args['bwa-mem-index-image'].value = wf.topvars['bwaMemIndexImage']

            # VEP annotation
            if wf.args.vep_cache_dir and wf.args.vep_plugin_dir:
                depend_task = filter_align_task if wf.topvars['bwaMemIndexImage'].value is not None else filter_task
                vep_task, args = wf.add_task(vep(tumor_sample), tag=tumor_sample, depends=[depend_task])
                args['input_file'].value = depend_task.outputs['out']
                args['fasta'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['dir_cache'].value = top_vars['vep_cache_dir']
                args['dir_plugins'].value = top_vars['vep_plugin_dir']

        # tumor only analysis
        if normal_sample.lower() == 'none' and tumor_sample.lower() != 'none':
            mutect_tasks = []
            for ind, interval_file in enumerate(interval_files):
                mutect_task, args = wf.add_task(Mutect2(f'{tumor_sample}-{ind}'), tag=f'{tumor_sample}-{ind}')
                mutect_tasks.append(mutect_task)
                mutect_task.depends = [bam_dict[normal_sample].task_id, bam_dict[tumor_sample].task_id]
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['tumor_bam'].value = bam_dict[tumor_sample].outputs['out']
                args['tumor_name'].value = tumor_sample
                args['germline-resource'].value = wf.topvars['germline_vcf']
                args['intervals'].value = interval_file
                args['pon'].value = wf.topvars['pon']
                args['alleles'].value = wf.topvars['alleles']

            # LearnReadOrientationModel
            lrom_task, args = wf.add_task(LearnReadOrientationModel(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['inputs'].value = [x.outputs['f1r2-tar-gz'] for x in mutect_tasks]

            # merge vcf
            merge_vcf_task, args = wf.add_task(MergeVcfs(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['inputs'].value = [x.outputs['out'] for x in mutect_tasks]

            # merge stats
            merge_stat_task, args = wf.add_task(MergeMutectStats(tumor_sample), tag=tumor_sample, depends=mutect_tasks)
            args['stats'].value = [x.outputs['stats'] for x in mutect_tasks]

            # filter
            filter_task, args = wf.add_task(FilterMutectCalls(tumor_sample), tag=tumor_sample, depends=[merge_vcf_task, merge_stat_task, lrom_task])
            args['vcf'].value = merge_vcf_task.outputs['out']
            args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
            args['ob-priors'].value = lrom_task.outputs['out']
            args['stats'].value = merge_stat_task.outputs['out']

            if wf.topvars['bwaMemIndexImage'].value is not None:
                filter_align_task, args = wf.add_task(FilterAlignmentArtifacts(tumor_sample), tag=tumor_sample, depends=[filter_task])
                args['vcf'].value = filter_task.outputs['out']
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['bam'].value = bam_dict[tumor_sample].outputs['out']
                args['bwa-mem-index-image'].value = wf.topvars['bwaMemIndexImage']

            if wf.args.vep_cache_dir and wf.args.vep_plugin_dir:
                depend_task = filter_align_task if wf.topvars['bwaMemIndexImage'].value is not None else filter_task
                vep_task, args = wf.add_task(vep(tumor_sample), tag=tumor_sample, depends=[depend_task])
                args['input_file'].value = depend_task.outputs['out']
                args['fasta'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['dir_cache'].value = top_vars['vep_cache_dir']
                args['dir_plugins'].value = top_vars['vep_plugin_dir']

        # germline variant calling
        if normal_sample.lower() != 'none':
            hap_tasks = []
            for ind, interval_file in enumerate(interval_files):
                hap_task, args = wf.add_task(Haplotyper(normal_sample), tag=tumor_sample, depends=[bam_dict[normal_sample].task_id, recal_dict[normal_sample].task_id])
                hap_tasks.append(hap_task)
                args['REFERENCE_SEQUENCE'].value = wf.topvars['ref'] if not make_index else index_task.outputs['ref_genome']
                args['bam'].value = bam_dict[normal_sample].outputs['out']
                args['intervals'].value = [interval_file]

            # merge vcf
            merge_vcf_task, args = wf.add_task(MergeVcfs(normal_sample), tag=normal_sample, depends=hap_tasks)
            args['inputs'].value = [x.outputs['out'] for x in hap_tasks]






