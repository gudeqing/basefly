import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, TmpVar
from utils.tidy_tools import merge_metrics, merge_star_fusion, merge_arcasHLA_genetype
from utils.get_fastq_info import get_fastq_info
from basefly.commands import RNASplitReadsAtJunction, recalibration, Haplotyper
__author__ = 'gdq'

"""
1. fastp
2. star | sentieon star
3. salmon-alignment-based
4. star-fusion
5. picard collectrnaseqmetrics QC: rRNA, mRNA percent
6. arcasHLA
to do list:
7. https://icbi-lab.github.io/immunedeconv/
"""


def fastp(sample):
    cmd = Command()
    cmd.meta.name = 'fastp'
    # cmd.runtime.image = 'gudeqing/fastp:0.21.0'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.tool = 'fastp'
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=5, desc='thread number')
    cmd.args['other_args'] = Argument(prefix='', default='', desc="other arguments you want to use, such as '-x val'")
    cmd.args['out1'] = Argument(prefix='-o ', value=TmpVar(value=f'{sample}.clean.R1.fq.gz', name='~{sample}.clean.R1.fq.gz'), type='str', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', value=TmpVar(value=f'{sample}.clean.R2.fq.gz', name='~{sample}.clean.R2.fq.gz'), type='str', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='-h ', value=TmpVar(value=f'{sample}.fastp.html', name='~{sample}.fastp.html'), type='str', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', value=TmpVar(value=f'{sample}.fastp.json', name='~{sample}.fastp.json') , type='str', desc='html report file')
    # 下面的outputs设置起初是为了能够生成wdl设置,
    cmd.outputs['out1'] = Output(value="{out1}", type='outfile')  # 这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd


def star(sample, platform='illumina', sentieon=False):
    """
    star alignment
    """
    cmd = Command()
    cmd.meta.name = 'star'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 25*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'sentieon STAR' if sentieon else 'STAR'
    cmd.args['threads'] = Argument(prefix='--runThreadN ', default=4, desc='threads to use')
    cmd.args['genomeDir'] = Argument(prefix='--genomeDir ', type='indir', desc='genome index directory')
    cmd.args['read1'] = Argument(prefix='--readFilesIn ', type='infile', array=True, desc='input read1 fastq file list', delimiter=',')
    cmd.args['read2'] = Argument(prefix=' ', type='infile', array=True, level='optional', desc='input read2 fastq file list', delimiter=',')
    cmd.args['outFileNamePrefix'] = Argument(prefix='--outFileNamePrefix ', value=TmpVar(name='~{sample}.', value=f'{sample}.'), desc='output file name prefix')
    cmd.args["outSAMattrRGline"] = Argument(prefix="--outSAMattrRGline ", value=TmpVar(name="ID:~{sample} SM:~{sample} PL:~{platform}", value=f"ID:{sample} SM:{sample} PL:{platform}"), desc="SAM/BAM read group line. The first word contains the read group identifier and must start with 'ID:")
    cmd.args['outSAMtype'] = Argument(prefix='--outSAMtype ', default="BAM SortedByCoordinate", desc='out bam type')
    cmd.args["outSAMunmapped"] = Argument(prefix="--outSAMunmapped ", default="Within", desc='if to include unmapped reads in sam')
    cmd.args["readFilesCommand"] = Argument(prefix="--readFilesCommand ", default="zcat", desc='method to unzip fastq')
    cmd.args["twopassMode"] = Argument(prefix="--twopassMode ", default="Basic", desc='')
    cmd.args["outFilterMultimapNmax"] = Argument(prefix="--outFilterMultimapNmax ", default=20, desc="max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped")
    cmd.args["alignSJoverhangMin"] = Argument(prefix="--alignSJoverhangMin ", default=8, desc="minimum overhang for unannotated junctions")
    cmd.args["alignSJDBoverhangMin"] = Argument(prefix="--alignSJDBoverhangMin ", default=10, desc="minimum overhang for annotated junctions")
    cmd.args["chimSegmentMin"] = Argument(prefix="--chimSegmentMin ", default=12, desc="parameter controls the minimum mapped length of the two segments that is allowed. For example, if you have 2x75 reads and used --chimSegmentMin 20, a chimeric alignment with 130b on one chromosome and 20b on the other will be output, while 135 + 15 won’t be")
    cmd.args["outFilterMismatchNoverLmax"] = Argument(prefix="--outFilterMismatchNoverLmax ", default=0.04, desc="alignment will be output only if its ratio of mismatches to mapped length is less than or equal to this value.")
    cmd.args["outFilterType"] = Argument(prefix="--outFilterType ", default="BySJout", desc="type of filtering")
    cmd.args["outSAMstrandField"] = Argument(prefix="--outSAMstrandField ", default="intronMotif", desc="include for potential use with StringTie for assembly")
    cmd.args["quantMode"] = Argument(prefix="--quantMode ", default="TranscriptomeSAM", desc='output bam file for downstream gene/transcript quantification')
    cmd.args["limitBAMsortRAM"] = Argument(prefix="--limitBAMsortRAM ", default=35000000000, desc="int>=0: maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value can only be used with –genomeLoad NoSharedMemory option.")
    cmd.args["limitIObufferSize"] = Argument(prefix="--limitIObufferSize ", default=200000000, desc="int>0: max available buffers size (bytes) for input/output, per thread")
    cmd.args["outSAMattrIHstart"] = Argument(prefix="--outSAMattrIHstart ", default=0, desc="start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.")
    cmd.args["alignMatesGapMax"] = Argument(prefix="--alignMatesGapMax ", default=100000, desc="maximum gap between two mates")
    cmd.args["alignIntronMax"] = Argument(prefix="--alignIntronMax ", default=100000, desc="maximum intron size")
    cmd.args["alignSJstitchMismatchNmax"] = Argument(prefix="--alignSJstitchMismatchNmax ", default="5 -1 5 5", desc="maximum number of mismatches for stitching of the splice junctions.(1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.")
    cmd.args["chimJunctionOverhangMin"] = Argument(prefix="--chimJunctionOverhangMin ", default=8, desc="minimum overhang for a chimeric junction")
    cmd.args["chimMultimapScoreRange"] = Argument(prefix="--chimMultimapScoreRange ", default=3, desc="the score range for multi-mapping chimeras below the best chimeric score. Only works with –chimMultimapNmax > 1")
    cmd.args["chimMultimapNmax"] = Argument(prefix="--chimMultimapNmax ", default=20, desc="maximum number of chimeric multi-alignments")
    cmd.args["chimNonchimScoreDropMin"] = Argument(prefix="--chimNonchimScoreDropMin ", default=10, desc="to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be smaller than this value")
    cmd.args["chimOutJunctionFormat"] = Argument(prefix="--chimOutJunctionFormat ", default=1, desc="if 1: add comment lines at the end of the file: command line and Nreads: total, unique, multi")
    cmd.args["peOverlapNbasesMin"] = Argument(prefix="--peOverlapNbasesMin ", default=12, desc="minimum number of overlap bases to trigger mates merging and realignment")
    cmd.args["peOverlapMMp"] = Argument(prefix="--peOverlapMMp ", default=0.1, desc="maximum proportion of mismatched bases in the overlap area")
    cmd.args["chimOutType"] = Argument(prefix="--chimOutType ", default="Junctions WithinBAM", desc="type of chimeric output")
    cmd.args["alignInsertionFlush"] = Argument(prefix="--alignInsertionFlush ", default="Right", desc="how to flush ambiguous insertion positions")
    cmd.args["chimScoreJunctionNonGTAG"] = Argument(prefix="--chimScoreJunctionNonGTAG ", default=-4, desc="penalty for a non-GT/AG chimeric junction")
    cmd.args["alignSplicedMateMapLminOverLmate"] = Argument(prefix="--alignSplicedMateMapLminOverLmate ", default=0, desc="alignSplicedMateMapLmin normalized to mate length")
    cmd.args["alignSplicedMateMapLmin"] = Argument(prefix="--alignSplicedMateMapLmin ", default=30, desc="minimum mapped length for a read mate that is spliced")
    cmd.args["quantTranscriptomeBan"] = Argument(prefix="--quantTranscriptomeBan ", default="Singleend", range=["IndelSoftclipSingleend", "Singleend"], desc= "prohibit various alignment type. Singleend is allowed for salmon.")
    # index, flagstats,idxstats
    cmd.args["_1"] = Argument(value=' && ', type='fix')
    cmd.args['_index'] = Argument(value=f"samtools index {sample}.Aligned.sortedByCoord.out.bam", type='fix')
    cmd.args["_2"] = Argument(value=' && ', type='fix')
    cmd.args['_flagstats'] = Argument(value=f'samtools flagstat {sample}.Aligned.sortedByCoord.out.bam > {sample}.align.flagstat.txt', type='fix')
    cmd.args["_3"] = Argument(value=' && ', type='fix')
    cmd.args['_idxstats'] = Argument(value=f"samtools idxstats {sample}.Aligned.sortedByCoord.out.bam > {sample}.align.idxstats.txt", type='fix')
    cmd.args["_4"] = Argument(value=' && ', type='fix')
    # 增加一行命令检查是否有chimeric输出，如果没有，则写一个空文件出来，这样不会导致后面的star-fusion的步骤失败
    cmd.args['_5'] = Argument(value=f"if [ ! -f {sample}.Chimeric.out.junction ]; then echo -e '# No output of chimeric\\n# Nreads 1\tNreadsUnique 1\tNreadsMulti 0' >  {sample}.Chimeric.out.junction ; fi", type='fix')
    cmd.outputs['bam'] = Output(value=f"{sample}.Aligned.sortedByCoord.out.bam")
    cmd.outputs['bam_bai'] = Output(value=f"{sample}.Aligned.sortedByCoord.out.bam.bai")
    cmd.outputs['transcript_bam'] = Output(value=f"{sample}.Aligned.toTranscriptome.out.bam")
    cmd.outputs['align_log'] = Output(value=f"{sample}.Log.final.out")
    cmd.outputs['chimeric'] = Output(value=f"{sample}.Chimeric.out.junction")
    cmd.outputs['sj'] = Output(value=f"{sample}.SJ.out.tab")
    cmd.outputs['flagstat'] = Output(value=f"{sample}.align.flagstat.txt")
    cmd.outputs['idxstats'] = Output(value=f"{sample}.align.idxstats.txt")
    return cmd


def salmon():
    """
    salmon quant in alignment-based mode
    ./bin/salmon quant -t transcripts.fa -l <LIBTYPE> -a aln.bam -o salmon_quant
    """
    cmd = Command()
    cmd.meta.name = 'salmon'
    cmd.meta.desc = 'gene/transcript expression quantification'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 2*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'salmon quant'
    cmd.args['threads'] = Argument(prefix='--threads ', default=4, desc='The number of threads that will be used for quasi-mapping, quantification, and bootstrapping / posterior sampling (if enabled).')
    cmd.args['libType'] = Argument(prefix='--libType ', default='A', desc="Salmon also has the ability to automatically infer (i.e. guess) the library type based on how the first few thousand reads map to the transcriptome. ")
    cmd.args['transcripts'] = Argument(prefix='-t ', type='infile', desc='transcript fasta file')
    cmd.args['geneMap'] = Argument(prefix='-g ', type='infile', desc="The transcript to gene mapping should be provided as either a GTF file, or a in a simple tab-delimited format where each line contains the name of a transcript and the gene to which it belongs separated by a tab.")
    cmd.args['gencode'] = Argument(prefix='--gencode ', type='bool', default=True, desc='expect the input transcript fasta to be in GENCODE format, and will split the transcript name at the first | character.')
    cmd.args['bam'] = Argument(prefix='-a ', type='infile', array=True, delimiter=' ',  desc='transcriptome bam files')
    cmd.args['outDir'] = Argument(prefix='-o ', type='str', default='.', desc='output directory')
    cmd.args['no-version-check'] = Argument(prefix='--no-version-check ', type='bool', default=True, desc='check version')
    cmd.args['gcBias'] = Argument(prefix='--gcBias ', type='bool', default=True, desc='perform gc Bias correction')
    cmd.args['seqBias'] = Argument(prefix='--seqBias ', type='bool', default=True, desc='perform seq Bias correction')
    cmd.args['posBias'] = Argument(prefix='--posBias ', type='bool', default=True, desc='This is meant to model non-uniform coverage biases that are sometimes present in RNA-seq data (e.g. 5’ or 3’ positional bias).')
    cmd.outputs['trans_exp'] = Output(path="{outDir}" + "/quant.sf")
    cmd.outputs['gene_exp'] = Output(path="{outDir}" + "/quant.genes.sf")
    cmd.outputs['outDir'] = Output(path="{outDir}", type='outdir')
    return cmd


def star_fusion():
    """
    star-fusion
    """
    cmd = Command()
    cmd.meta.name = 'star-fusion'
    cmd.meta.source = "https://github.com/STAR-Fusion/STAR-Fusion"
    cmd.meta.version = 'v1.10.0'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 10*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = ' STAR-Fusion'
    cmd.args['threads'] = Argument(prefix='--CPU ', default=4, desc='The number of threads')
    cmd.args['read1'] = Argument(prefix='--left_fq ', level='optional', type='infile', array=True, delimiter=',', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='--right_fq ', level='optional', type='infile',  array=True, delimiter=',', desc='read2 fastq file')
    cmd.args['chimeric_junction'] = Argument(prefix='-J ', type='infile', desc="generated file called 'Chimeric.out.junction' by STAR alignment")
    cmd.args['genomeLibDir'] = Argument(prefix='--genome_lib_dir ', type='indir', desc='ctat_genome_lib_build_dir which contains fusion database')
    cmd.args['outdir'] = Argument(prefix='--output_dir ', default='.', desc='output dir')
    cmd.args['tmpdir'] = Argument(prefix='--tmpdir ', default='.', desc='temp file dir')
    cmd.args['FusionInspector'] = Argument(prefix='--FusionInspector ', level='optional', range=["inspect", "validate"], desc="FusionInspector that provides a more in-depth view of the evidence supporting the predicted fusions.")
    cmd.args['examine_coding_effect'] = Argument(prefix='--examine_coding_effect', type='bool', desc="explore impact of fusions on coding sequences")
    cmd.args['denovo_reconstruct'] = Argument(prefix='--denovo_reconstruct', type='bool', desc="attempt to reconstruct fusion transcripts using Trinity de novo assembly (requires --FusionInspector)")
    cmd.outputs['fusion_predictions'] = Output(value="{outdir}/star-fusion.fusion_predictions.tsv")
    cmd.outputs['fusion_predictions_abridged'] = Output(value="{outdir}/star-fusion.fusion_predictions.abridged.tsv")
    cmd.outputs['fusion_inspector_validate_fusions_abridged'] = Output(value="{outdir}/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv")
    cmd.outputs['fusion_inspector_validate_web'] = Output(value="{outdir}/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv")
    cmd.outputs['fusion_inspector_inspect_web'] = Output(value="{outdir}/FusionInspector-validate/finspector.fusion_inspector_web.html")
    cmd.outputs['fusion_inspector_inspect_fusions_abridged'] = Output(value="{outdir}/FusionInspector-inspect/finspector.FusionInspector.fusions.abridged.tsv")
    return cmd


def collect_metrics(sample):
    """
    picard to collect rnaseq metrics
    """
    cmd = Command()
    cmd.meta.name = 'CollectRnaSeqMetrics'
    # cmd.runtime.image = 'broadinstitute/picard:latest'
    # cmd.runtime.tool = 'java -jar /usr/picard/picard.jar CollectRnaSeqMetrics'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 10*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'java -jar /usr/local/src/picard.jar CollectRnaSeqMetrics'
    cmd.args['bam'] = Argument(prefix='I=', type='infile', desc='input bam file')
    cmd.args['out_metrics'] = Argument(prefix='O=', value=f'{sample}.RnaSeqMetrics.txt', desc='output metric file')
    cmd.args['strandness'] = Argument(prefix='STRAND_SPECIFICITY=', default='NONE', desc='strand-specificity')
    cmd.args['ref_flat'] = Argument(prefix='REF_FLAT=', type='infile', desc='Gene annotations in refFlat form.  Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat')
    cmd.args['ribosomal_interval'] = Argument(prefix='RIBOSOMAL_INTERVALS=', type='infile', desc="Location of rRNA sequences in genome, in interval_list format.  If not specified no bases will be identified as being ribosomal.")
    # cmd.args['cov_pdf'] = Argument(prefix='CHART_OUTPUT=', value=f'{sample}.coverage.pdf', desc='coverage output file')
    cmd.outputs['metrics'] = Output(value="{out_metrics}")
    # cmd.outputs['cov_pdf'] = Output(value="{cov_pdf}")
    return cmd


def arcas_hla(threads=4):
    """
    set -e
    arcasHLA extract --unmapped -t ~{threads} -o ./ ~{bam}
    arcasHLA genotype --min_count 75 -t ~{threads} -o ./ *.1.fq.gz *.2.fq.gz
    arcasHLA merge
    """
    cmd = Command()
    cmd.meta.name = 'arcasHLA'
    cmd.meta.version = '0.2.5'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 5*1024**3
    cmd.runtime.cpu = 2
    # 软链接数据库
    cmd.args['database'] = Argument(prefix='rm -r /home/arcasHLA-master/dat && ln -s ', type='indir', level='optional', desc='database of arcas_software')
    cmd.args['link'] = Argument(prefix='/home/arcasHLA-master/dat &', type='bool', default=False)
    # run software
    cmd.args['_1'] = Argument(value=f'arcasHLA extract --temp ./ --unmapped -t {threads} -o .', type='fix')
    cmd.args['bam'] = Argument(value='', type='infile', desc='input bam file')
    cmd.args['_2'] = Argument(value=' && arcasHLA genotype', type='fix')
    cmd.args['_3'] = Argument(value=f'--min_count 75 -t {threads} -o ./ *.1.fq.gz *.2.fq.gz &&', type='fix')
    cmd.args['_4'] = Argument(value='arcasHLA merge', type='fix')
    cmd.outputs['hla_genotype'] = Output(value='*.genotype.json')
    return cmd


def quant_merge():
    cmd = Command()
    cmd.meta.name = 'quantMerge'
    cmd.meta.desc = 'Merge multiple quantification results into a single file'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 2*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'salmon quantmerge'
    # 下面的quants参数对应的是目录，所以type='indir'
    cmd.args['quants'] = Argument(prefix="--quants ", array=True, type='indir', desc='salmon quant dir list')
    cmd.args['names'] = Argument(prefix='--names ', array=True, level='optional', desc='sample name to assign for each input')
    cmd.args['column'] = Argument(prefix='--column ', default='TPM', range=['TPM', 'NumReads'], desc='indicate which column to merge')
    cmd.args['genes'] = Argument(prefix='--genes ', type='bool', default=False, desc='indicate if to merge gene data, default to merge transcript data')
    cmd.args['out'] = Argument(prefix='--output ', desc='merged result file')
    cmd.outputs['out'] = Output(path="{out}")
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'RnaSeqPipeline'
    wf.meta.desc = 'This is a pipeline for rnaseq analysis'
    wf.init_argparser()
    # add workflow args
    wf.add_argument('-star_index', required=True, help='star alignment index dir')
    wf.add_argument('--sentieon', default=False, action="store_true", help='indicate to use sentieon STAR tool')
    wf.add_argument('-fusion_index', required=True, help='star-fusion database dir: CTAT_resource_lib from https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz')
    wf.add_argument('-transcripts_fa', required=True, help='transcriptome fasta file')
    wf.add_argument('-genome_fa', required=False, help='genome fasta file, needed for variant calling')
    wf.add_argument('-gtf', required=True, help='genome annotation file')
    wf.add_argument('-ref_flat', required=True, help='gene model file, genome annotation file, download site: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz')
    wf.add_argument('-rRNA_interval', required=True, help='rRNA interval file generated by using picard BedToIntervalList command')
    wf.add_argument('-hla_db', required=False, help='arcasHLA database, please refer to https://github.com/RabadanLab/arcasHLA')
    wf.add_argument('-fastq_info', nargs='+', required=True, help='A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]')
    wf.add_argument('-r1_name', default='(.*).R1.fastq', help="python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'")
    wf.add_argument('-r2_name', default='(.*).R2.fastq', help="python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'")
    wf.add_argument('-exclude_samples', default=tuple(), nargs='+', help='samples to exclude from analysis')
    wf.add_argument('-dbsnp', required=False, help='dbsnp vcf file, only needed for variant calling')
    wf.add_argument('-known_indels', required=False, help='high confidence known indel vcf file, only needed for variant calling')
    wf.add_argument('-known_mills', required=False, help='high confidence known indel vcf file, only needed for variant calling')
    wf.parse_args()

    # add top_vars
    wf.add_topvars(dict(
        starIndex=TopVar(value=wf.args.star_index, type='indir'),
        fusionIndex=TopVar(value=wf.args.fusion_index, type='indir'),
        transcripts=TopVar(value=wf.args.transcripts_fa, type='infile'),
        gtf=TopVar(value=wf.args.gtf, type='infile'),
        ref_flat=TopVar(value=wf.args.ref_flat, type='infile'),
        rRNA_interval=TopVar(value=wf.args.rRNA_interval, type='infile'),
        genome_fa=TopVar(value=wf.args.genome_fa, type='infile'),
        known_dbsnp=TopVar(value=wf.args.dbsnp, type='infile'),
        known_indels=TopVar(value=wf.args.known_indels, type='infile'),
        known_mills=TopVar(value=wf.args.known_mills)
    ))
    if wf.args.hla_db:
        wf.topvars['hla_database'] = TopVar(value=wf.args.hla_db, type='indir')

    fastq_info = get_fastq_info(fastq_info=wf.args.fastq_info, r1_name=wf.args.r1_name, r2_name=wf.args.r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found !')

    for sample, (r1s, r2s) in fastq_info.items():
        if sample in wf.args.exclude_samples:
            continue
        # 一个样本可能有多个fastq
        fastp_tasks = []
        for ind, (r1, r2) in enumerate(zip(r1s, r2s)):
            # 向流程中添加task
            fastp_task, args = wf.add_task(fastp(sample), name=f'fastp-{sample}-{ind}')
            args['read1'].value = TmpVar(name='read1', value=r1, type='infile')
            args['read2'].value = TmpVar(name='read2', value=r2, type='infile')
            fastp_tasks.append(fastp_task)

        # star alignment
        fastp_task_ids = [x.task_id for x in fastp_tasks]
        star_task, args = wf.add_task(star(sample, sentieon=wf.args.sentieon), name='star-'+sample, depends=fastp_task_ids)
        args['read1'].value = [x.outputs["out1"] for x in fastp_tasks]
        args['read2'].value = [x.outputs["out2"] for x in fastp_tasks]
        args['genomeDir'].value = wf.topvars['starIndex']

        # salmon quant
        salmon_task, args = wf.add_task(salmon(), name='salmon-'+sample, depends=[star_task.task_id])
        args['bam'].value = [star_task.outputs['transcript_bam']]
        args['transcripts'].value = wf.topvars['transcripts']
        args['geneMap'].value = wf.topvars['gtf']

        # fusion identification
        fusion_task, args = wf.add_task(star_fusion(), name='starfusion-'+sample, depends=[star_task.task_id])
        args['genomeLibDir'].value = wf.topvars['fusionIndex']
        args['chimeric_junction'].value = star_task.outputs['chimeric']
        # args['read1'].value = [x.outputs["out1"] for x in fastp_tasks]
        # args['read2'].value = [x.outputs["out2"] for x in fastp_tasks]

        # collectRNAseqMetrics
        metric_task, args = wf.add_task(collect_metrics(sample), name='collectMetrics-'+sample, depends=[star_task.task_id])
        args['bam'].value = star_task.outputs['bam']
        args['ref_flat'].value = wf.topvars['ref_flat']
        args['ribosomal_interval'].value = wf.topvars['rRNA_interval']

        # HLA-typing
        hla_task, args = wf.add_task(arcas_hla(), name='hla-'+sample, depends=[star_task.task_id])
        args['bam'].value = star_task.outputs['bam']
        if wf.args.hla_db:
            args['database'].value = wf.topvars['hla_database']
            args['link'].value = True

        # split bam
        split_task, args = wf.add_task(RNASplitReadsAtJunction(sample), name=f'splitBam-{sample}', depends=[star_task.task_id])
        args['bam'].value = star_task.outputs['bam']
        args['ref'].value = wf.topvars['genome_fa']

        # recalibration
        recal_task, args = wf.add_task(recalibration(sample), name=f'recal-{sample}', depends=[split_task.task_id])
        args['bam'].value = split_task.outputs['out_bam']
        args['ref'].value = wf.topvars['genome_fa']
        args['database'].value = [wf.topvars['known_dbsnp'], wf.topvars['known_indels'], wf.topvars['known_mills']]

        # haplotyper
        hap_task, args = wf.add_task(Haplotyper(sample), name=f'haplotyper-{sample}', depends=[recal_task.task_id])
        args['ref'].value = wf.topvars['genome_fa']
        args['bam'].value = split_task.outputs['out_bam']
        args['recal_data'].value = recal_task.outputs['recal_data']
        args['trim_soft_clip'].value = True
        args['call_conf'].value = 20
        args['emit_conf'].value = 20

    # merge transcript/gene TPM/Count
    depends = [k for k, v in wf.tasks.items() if v.name.startswith('salmon')]
    for feature in ['gene', 'trans']:
        for data_type in ['TPM', 'NumReads']:
            name = f'merge-{feature}-{data_type}'
            task, args = wf.add_task(quant_merge(), name=name, depends=depends)
            task.cmd.meta.name = name
            args['column'].value = data_type
            args['genes'].value = True if feature == 'gene' else False
            args['quants'].value = [wf.tasks[task_id].outputs['outDir'] for task_id in depends]
            args['names'].value = [wf.tasks[task_id].name.split('-', 1)[1] for task_id in depends]
            args['out'].value = f'{feature}.{data_type}.txt'
    wf.run()

    # merger result
    if wf.success:
        print('Merging Results......')
        merge_metrics(wf.args.outdir, filter_ref='', outdir=os.path.join(wf.args.outdir, 'merge_qc'))
        merge_arcasHLA_genetype(wf.args.outdir, outdir=os.path.join(wf.args.outdir, 'merge_HLA'))
        merge_star_fusion(wf.args.outdir, outdir=os.path.join(wf.args.outdir, 'merge_starfusion'))
        print('...end...')


if __name__ == '__main__':
    pipeline()