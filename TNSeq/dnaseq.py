import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from nestcmd.nestcmd import Argument, Output, Command, Workflow, TopVar, TmpVar, get_fastq_info
__author__ = 'gdq'

"""
参考：https://github.com/Sentieon/sentieon-scripts/blob/master/example_pipelines/somatic/TNseq/tumor_normal.sh
参考：https://support.sentieon.com/manual/TNseq_usage/tnseq/#
本版本不刻意考虑是否可以正常生成wdl
"""


def fastp(sample):
    cmd = Command()
    cmd.meta.name = 'fastp'
    # cmd.runtime.image = 'gudeqing/fastp:0.21.0'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.tool = 'fastp'
    # 可以直接用访问属性的方式添加参数，这个得益于使用Munch对象而不是原生字典
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', desc='read2 fastq file')
    cmd.args['other_args'] = Argument(prefix='', default='', desc="other arguments you want to use, such as '-x val'")
    # 当然，可以直接用字典的方式添加参数
    cmd.args['out1'] = Argument(prefix='-o ', value=TmpVar(value=f'{sample}.clean.R1.fq.gz', name='~{sample}.clean.R1.fq.gz'), type='str', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', value=TmpVar(value=f'{sample}.clean.R2.fq.gz', name='~{sample}.clean.R2.fq.gz'), type='str', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='-h ', value=TmpVar(value=f'{sample}.qc.html', name='~{sample}.qc.html'), type='str', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', value=TmpVar(value=f'{sample}.qc.json', name='~{sample}.qc.json') , type='str', desc='html report file')
    # 下面的outputs设置起初是为了能够生成wdl设置,
    cmd.outputs['out1'] = Output(value="{out1}", type='outfile')  # 这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd


def bwa_mem(sample, platform):
    cmd = Command()
    cmd.meta.name = 'bwa_mem'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon bwa mem -M'
    cmd.args['readgroup'] = Argument(prefix='-R ', desc='read group info', value=f'"@RG\\tID:{sample}\\tSM:{sample}\\tPL:{platform}"')
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['k'] = Argument(prefix='-K ', default=10000000)
    cmd.args['ref'] = Argument(type='infile', desc='reference fasta file')
    cmd.args['read1'] = Argument(type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(type='infile', desc='read2 fastq file')
    cmd.args['_x2'] = Argument(type='fix', value=' | sentieon util sort')
    cmd.args['ref2'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['out'] = Argument(prefix='-o ', desc='output bam file', value=f'{sample}.sorted.bam')
    cmd.args['t2'] = Argument(prefix='-t ', default=16, desc='number of threads to use')
    cmd.args['_x3'] = Argument(type='fix', value='--sam2bam -i -')
    cmd.outputs['out'] = Output(path="{out}", type='outfile')
    return cmd


def get_metrics(sample):
    cmd = Command()
    cmd.meta.name = 'get_metrics'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile',  desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile',  desc='input bam file')
    cmd.args['mq_metrics'] = Argument(prefix='--algo MeanQualityByCycle ', value=f'{sample}.mq_metrics.txt', desc='metric file of MeanQualityByCycle')
    cmd.args['qd_metrics'] = Argument(prefix='--algo QualDistribution ', value=f'{sample}.qd_metrics.txt', desc='metric file of QualDistribution')
    cmd.args['gc_summary'] = Argument(prefix='--algo GCBias --summary ', value=f'{sample}.gc_summary.txt', desc='summary file of GCBias')
    cmd.args['gc_metrics'] = Argument(desc='metrics file of GCBias', value=f'{sample}.gc_metrics.txt')
    cmd.args['aln_metrics'] = Argument(prefix='--algo AlignmentStat ', value=f'{sample}.aln_metrics.txt', desc='aln_metrics file of AlignmentStat')
    cmd.args['insert_metrics'] = Argument(prefix='--algo InsertSizeMetricAlgo ', value=f'{sample}.insert_metrics.txt', desc='insert_metrics file of InsertSizeMetricAlgo')
    cmd.outputs['mq_metrics'] = Output(path='{mq_metrics}')
    cmd.outputs['qd_metrics'] = Output(path='{qd_metrics}')
    cmd.outputs['gc_summary'] = Output(path='{gc_summary}')
    cmd.outputs['gc_metrics'] = Output(path='{gc_metrics}')
    cmd.outputs['aln_metrics'] = Output(path='{aln_metrics}')
    cmd.outputs['insert_metrics'] = Output(path='{insert_metrics}')
    return cmd


def plot_metrics(sample, method='GCBias'):
    cmd = Command()
    cmd.meta.name = f'plot{method}'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon plot'
    cmd.args['method'] = Argument(desc='method of plot', default=method)
    cmd.args['out'] = Argument(prefix='-o ', desc='plot file', value=f'{sample}.{method}.pdf')
    cmd.args['i'] = Argument(type='infile', desc='input metrics file for plot')
    cmd.outputs['out'] = Output(path='{out}')
    return cmd


def locus_collector(sample):
    cmd = Command()
    cmd.meta.name = 'LocusCollector'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['score'] = Argument(prefix='--algo LocusCollector --fun score_info ', desc='output score file', value=f'{sample}.score.txt')
    cmd.outputs['score'] = Output(path='{score}')
    return cmd


def dedup(sample):
    cmd = Command()
    cmd.meta.name = 'DeDup'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['_x'] = Argument(type='fix', value='--algo Dedup')
    cmd.args['score'] = Argument(prefix='--score_info ', type='infile', desc='score info file')
    cmd.args['dedup_metrics'] = Argument(prefix='--metrics ', desc='output metrics info file', value=f'{sample}.dedup.metrics.txt')
    cmd.args['deduped_bam'] = Argument(desc='output metrics info file', value=f'{sample}.deduped.bam')
    cmd.outputs['dedup_metrics'] = Output(path='{dedup_metrics}')
    cmd.outputs['out_bam'] = Output(path='{deduped_bam}')
    return cmd


def coverage_metrics(sample):
    cmd = Command()
    cmd.meta.name = 'CoverageMetrics'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['coverage_metrics'] = Argument(prefix='--algo CoverageMetrics ', value=f'{sample}.cov.metrics.txt', desc='output coverage metrics file')
    cmd.outputs['coverage_metrics'] = Output(path="{coverage_metrics}")
    return cmd


def realign(sample):
    cmd = Command()
    cmd.meta.name = 'realign'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args._x = Argument(type='fix', value='--algo Realigner')
    cmd.args['database'] = Argument(prefix='-k ', type='infile', multi_times=True, desc='known indel vcf file')
    cmd.args['realigned_bam'] = Argument(desc='output realigned bam file', value=f'{sample}.realigned.bam')
    cmd.outputs['out_bam'] = Output(path='{realigned_bam}')
    return cmd


def recalibration(sample):
    cmd = Command()
    cmd.meta.name = 'recalibration'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['_x'] = Argument(type='fix', value='--algo QualCal')
    cmd.args['database'] = Argument(prefix='-k ', type='infile', multi_times=True, desc='known indel vcf file')
    cmd.args['recal_data'] = Argument(desc="output recal_data.table", value=f'{sample}.recal_data.table')
    cmd.outputs['recal_data'] = Output(path='{recal_data}')
    return cmd


def TNhaplotyper2(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'TNhaplotyper2'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    # basic inputs
    cmd.args['bams'] = Argument(prefix='-i ', type='infile', multi_times=True, desc='reccaled tumor and normal bam list')
    cmd.args['recal_datas'] = Argument(prefix='-q ', type='infile', multi_times=True, desc='tumor and normal recal data list')
    cmd.args['method'] = Argument(prefix='--algo ', value='TNhaplotyper2', type='fix')
    cmd.args['tumor_sample'] = Argument(prefix='--tumor_sample ', desc='tumor sample name', default=tumor_sample)
    cmd.args['normal_sample'] = Argument(prefix='--normal_sample ', desc='normal sample name', level='optional')
    # optional inputs
    cmd.args['germline_vcf'] = Argument(prefix='--germline_vcf ', type='infile', level='optional', desc='the location of the population germline resource')
    cmd.args['pon'] = Argument(prefix='--pon ', type='infile', level='optional', desc='the location and name of panel of normal VCF file')
    cmd.args['out_vcf'] = Argument(value=f'{tumor_sample}.TNhaplotyper2.vcf.gz', desc='output vcf file of TNhaplotyper2, this will be used later for filtering')
    # orientation
    cmd.args['orientation_sample'] = Argument(prefix='--algo OrientationBias --tumor_sample ', level='optional', desc='tumor sample name')
    cmd.args['orientation_data'] = Argument(level='optional', default=f'{tumor_sample}.orientation.data', desc='output orientation bias result file')
    # contamination, 如果无对照样本或者germline vcf，则无该项分析
    cmd.args['contamination_tumor'] = Argument(prefix="--algo ContaminationModel --tumor_sample ", level='optional', desc='tumor sample name', default=tumor_sample)
    cmd.args['contamination_normal'] = Argument(prefix="--normal_sample ", level='optional', desc='normal sample name')
    cmd.args['germline_vcf2'] = Argument(prefix='--vcf ', type='infile', level='optional', desc='the location of the population germline resource')
    cmd.args['tumor_segments'] = Argument(prefix='--tumor_segments ', level='optional', default=f'{tumor_sample}.contamination.segments', desc='output file name of the file containing the tumor segments information produced by ContaminationModel')
    cmd.args['contamination_data'] = Argument(level='optional', default=f'{tumor_sample}.contamination.data', desc='output file containing the contamination information produced by ContaminationModel')
    cmd.outputs['out_vcf'] = Output(path='{out_vcf}')
    cmd.outputs['orientation_data'] = Output(path='{orientation_data}')
    cmd.outputs['tumor_segments'] = Output(path='{tumor_segments}')
    cmd.outputs['contamination_data'] = Output(path='{contamination_data}')
    return cmd


def TNfilter(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'TNfilter'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['method'] = Argument(type='fix', value='--algo TNfilter')
    cmd.args['tumor_sample'] = Argument(prefix='--tumor_sample ', desc='tumor sample name', default=tumor_sample)
    cmd.args['normal_sample'] = Argument(prefix='--normal_sample ', level='optional', desc='normal sample name')
    cmd.args['tmp_vcf'] = Argument(prefix='-v ', type='infile', desc='vcf file from TNhaplotyper2')
    cmd.args['contamination'] = Argument(prefix='--contamination ', type='infile', level='optional', desc='file containing the contamination information produced by ContaminationModel')
    cmd.args['tumor_segments'] = Argument(prefix='--tumor_segments ', type='infile', level='optional', desc='file containing the tumor segments information produced by ContaminationModel')
    cmd.args['orientation_data'] = Argument(prefix='--orientation_priors ', type='infile', level='optional', desc='file containing the orientation bias information produced by OrientationBias')
    cmd.args['out_vcf'] = Argument(desc='final output vcf', value=f'{tumor_sample}.final.vcf.gz')
    cmd.outputs['out_vcf'] = Output(path='{out_vcf}')
    return cmd


def Haplotyper(normal_sample):
    cmd = Command()
    cmd.meta.name = 'Haplotyper'
    cmd.runtime.image = 'registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='reccaled tumor and normal bam list')
    cmd.args['recal_data'] = Argument(prefix='-q ', type='infile', desc='tumor and normal recal data list')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['method'] = Argument(type='fix', value='--algo Haplotyper')
    cmd.args['emit_mode'] = Argument(prefix='--emit_mode ', default='gvcf', desc='determines what calls will be emitted. possible values:variant,confident,all,gvcf')
    cmd.args['ploidy'] = Argument(prefix='--ploidy ', type='int', default=2, desc='determines the ploidy number of the sample being processed. The default value is 2.')
    cmd.args['out_vcf'] = Argument(value=f'{normal_sample}.g.vcf.gz', desc='output vcf file')
    cmd.outputs['out_vcf'] = Output(path='{out_vcf}')
    cmd.outputs['out_vcf_idx'] = Output(path='{out_vcf}.tbi')
    return cmd


def GVCFtyper(normal_sample):
    cmd = Command()
    cmd.meta.name = 'GVCFtyper'
    cmd.runtime.image = 'registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['method'] = Argument(type='fix', value='--algo GVCFtyper')
    cmd.args['in_gvcf'] = Argument(prefix='-v ', type='infile', multi_times=True,  desc='input gvcf file')
    cmd.args['known_dbsnp'] = Argument(prefix='-d ', type='infile', desc='dbsnp file')
    cmd.args['call_conf'] = Argument(prefix='--call_conf ', type='int', default=30, desc="determine the threshold of variant quality to emit a variant. Variants with quality less than CONFIDENCE will be not be added to the output VCF file.")
    cmd.args['genotype_model'] = Argument(prefix='--genotype_model ', range={"coalescent", "multinomial"}, default='multinomial', desc="determines which model to use for genotyping and QUAL calculation")
    cmd.args['out_vcf'] = Argument(value=f'{normal_sample}.vcf.gz', desc='output vcf file')
    cmd.outputs['out_vcf'] = Output(path='{out_vcf}')
    cmd.outputs['out_vcf_idx'] = Output(path='{out_vcf}.tbi')
    return cmd


def snpEff(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'snpEff'
    cmd.runtime.image = '?'
    cmd.runtime.tool = 'java -Xmx9g snpEff.jar ann'
    cmd.args['genome_version'] = Argument(default='hg19', desc='human genome version')
    cmd.args['data_dir'] = Argument(prefix='-dataDir ', type='indir',  desc='Override data_dir parameter from config file')
    cmd.args['cancer'] = Argument(prefix='-cancer ', type='bool', default=True, desc='Perform cancer comparisons (Somatic vs Germline)')
    cmd.args['cancerSamples'] = Argument(prefix='-cancerSamples ', level='optional', type='infile', desc="Two column TXT file defining 'oringinal derived' samples. If '-cancer' used and the file is missing, then the last sample will be assumed as tumor sample.")
    cmd.args['canon'] = Argument(prefix='-canon ', type='bool', default=False, desc='Only use canonical transcripts')
    cmd.args['interval'] = Argument(prefix='-interval ', type='infile', level='optional', multi_times=True, desc="Use a custom intervals in TXT/BED/BigBed/VCF/GFF file (you may use this option many times)")
    cmd.args['other_args'] = Argument(level='optional', desc="other arguments that you want to input for the program, such as '-motif'")
    cmd.args['in_vcf'] = Argument(desc='input variant file', type='infile')
    cmd.args['_x'] = Argument(type='fix', value='>')
    cmd.args['out_vcf'] = Argument(desc='output annotated file', value=f'{tumor_sample}.final.annot.vcf')
    cmd.outputs['out_vcf'] = Output(path='{out_vcf}')
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
    cmd.args['fork'] = Argument(prefix='--fork ', type='int', default=4, desc='Use forking(multi-cpu/threads) to improve script runtime')
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
    cmd.args['filter_common'] = Argument(prefix='--filter_common ', type='bool', default=True, desc="Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters.")
    cmd.args['other_args'] = Argument(default='', desc='specify other arguments that you want to append to the command')
    cmd.outputs['out_vcf'] = Output(path='{output_file}')
    cmd.outputs['out_vcf_idx'] = Output(path='{output_file}.tbi')
    return cmd


def CombineVariants(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'CombineVariants'
    cmd.meta.desc = "Combine variants"
    cmd.runtime.image = 'broadinstitute/gatk3:3.8-1'
    cmd.runtime.tool = 'java -Xmx10g -jar GenomeAnalysisTK.jar -T CombineVariants'
    cmd.args['ref'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['variant'] = Argument(prefix='--variant ', multi_times=True, type='infile', desc='variant vcf file array')
    cmd.args['out_vcf'] = Argument(prefix='-o ', value=f'{tumor_sample}.combined_germline.vcf')
    cmd.args['assumeIdenticalSamples'] = Argument(prefix='--assumeIdenticalSamples', type='bool', desc='If true, assume input VCFs have identical sample sets and disjoint calls. This option allows the user to perform a simple merge (concatenation) to combine the VCFs.')
    cmd.outputs['combined_vcf'] = Output(path='{out_vcf}')
    return cmd


def SortVcf(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'SortVcf'
    cmd.meta.desc = "sort vcf"
    cmd.runtime.image = 'broadinstitute/picard:latest'
    cmd.runtime.tool = 'java -jar /usr/picard/picard.jar SortVcf'
    cmd.args['in_vcf'] = Argument(prefix='I=', type='infile', desc='input vcf to sort')
    cmd.args['out_vcf'] = Argument(prefix='O=', value=f'{tumor_sample}.combined_germline.sorted.vcf', type='infile', desc='output sorted vcf')
    cmd.outputs['sorted_vcf'] = Output(path='{out_vcf}')
    return cmd


def ReadBackedPhasing(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'ReadBackedPhasing'
    cmd.meta.desc = "ReadBackedPhasing"
    cmd.runtime.image = 'broadinstitute/gatk3:3.8-1'
    cmd.runtime.tool = 'java -Xmx10g -jar GenomeAnalysisTK.jar -T ReadBackedPhasing'
    cmd.args['ref'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='tumor bam file')
    cmd.args['variant'] = Argument(prefix='--variant ', type='infile', desc='input vcf file')
    cmd.args['interval'] = Argument(prefix='-L ', type='infile', desc='input vcf file')
    cmd.args['out_vcf'] = Argument(prefix='-o ', value=f'{tumor_sample}.phased.vcf', desc='output vcf file')
    cmd.outputs['phased_vcf'] = Output(path='{out_vcf}')
    return cmd


def HLA_ABC_typer(sample):
    cmd = Command()
    cmd.meta.name = 'OptiType'
    cmd.meta.desc = "OptiType:4-digit HLA typer"
    cmd.runtime.image = 'fred2/optitype:1.3.1'
    cmd.runtime.tool = 'OptiTypePipeline.py'
    cmd.args['reads'] = Argument(prefix='--input ', type='infile', array=True, desc='fastq file(s) (fished or raw) or .bam files stored for re-use, generated by an earlier OptiType run.')
    cmd.args['is_dna'] = Argument(prefix='--dna', type='bool', default=True, desc='use with DNA sequencing data')
    cmd.args['is_rna'] = Argument(prefix='--rna', type='bool', default=False, desc='use with RNA sequencing data')
    cmd.args['enumerate'] = Argument(prefix='--enumerate ', type='int', default=1, desc='Number of enumerations. OptiType will output the optimal solution and the top N-1 suboptimal solutions in the results CSV.')
    cmd.args['outdir'] = Argument(prefix='--outdir ', default='.', desc='Specifies the out directory to which all files should be written.')
    cmd.args['prefix'] = Argument(prefix='--prefix ', value=sample, desc='prefix of output files')
    cmd.args['config'] = Argument(prefix='--config ', default='config.ini', desc='config.ini file')
    cmd.outputs['result_tsv'] = Output(path='{prefix}_result.tsv')
    cmd.outputs['result_pdf'] = Output(path='{prefix}_coverage_plot.pdf')
    return cmd


def pipeline(fastq_dirs=None, fastq_files=None, r1_name="(.*).r1.fq", r2_name="(.*).r2.fq",
             pair_info_txt=None, optiType=False, ref_fa='hg19.fa', dbsnp='dbsnp.vcf', known_indels='known.indels.vcf',
             known_mills='mills.vcf', vep_cache='vep_cache/', vep_plugins='vep_plugins/',
             intervals=None, threads=15, germline_vcf=None, pon=None, do_realign=False,
             outdir='Results', run=False, no_docker=False, workers=3, retry=1,
             no_monitor_resource=False, no_check_resource=False):
    if fastq_files is None and fastq_dirs is None:
        print('this is a simple test')
        fastq_files = ('normal.r1.fq', 'normal.r2.fq', 'tumor.r1.fq', 'tumor.r2.fq')
        pair_info_txt = None
        run = False

    top_vars = dict(
        thread_number=TopVar(value=threads, type='int'),
        ref=TopVar(value=ref_fa, type='infile'),
        known_dbsnp=TopVar(value=dbsnp, type='infile'),
        known_indel=TopVar(value=known_indels, type='infile'),
        known_mills=TopVar(value=known_mills),
        pon=TopVar(value=pon if pon else ''),
        germline_vcf=TopVar(value=germline_vcf if germline_vcf else ''),
        vep_cache_dir=TopVar(value=vep_cache, type='indir'),
        vep_plugin_dir=TopVar(value=vep_plugins, type='indir'),
        intervals=TopVar(value=intervals if intervals else '', type='infile')
    )
    wf = Workflow(top_vars=top_vars)
    wf.meta.name = 'DNAseqPipeline'
    wf.meta.desc = 'typical bioinformatics pipeline using sentieon TNSeq and VEP, including HLA-typing'

    fastq_info = get_fastq_info(fastq_dirs=fastq_dirs, fastq_files=fastq_files, r1_name=r1_name, r2_name=r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found !')

    recal_dict = dict()
    bam_dict = dict()
    # batch是分组信息，用于wdl的scatter判断
    for sample, (r1s, r2s) in fastq_info.items():
        read1 = r1s[0]  # 假设每个样本只有对应一对fastq文件，不存在1对多的情况
        read2 = r2s[0]

        # fastp
        fastp_task, args = wf.add_task(fastp(sample), name=f'fastp-{sample}')
        args['read1'].value = TmpVar(name='read1', value=read1, type='infile')
        args['read2'].value = TmpVar(name='read2', value=read2, type='infile')

        # optiType
        if optiType:
            task, args = wf.add_task(HLA_ABC_typer(sample), name=f'optiType-{sample}', depends=[fastp_task.task_id])
            args['reads'].value = [fastp_task.outputs['out1'], fastp_task.outputs['out2']]

        # mapping
        mapping, args = wf.add_task(bwa_mem(sample, platform='ILLUMINA'), name=f'bwaMem-{sample}', depends=[fastp_task.task_id])
        args['t'].value = top_vars['thread_number']
        args['ref'].value = top_vars['ref']
        args['ref2'].value = top_vars['ref']
        args['read1'].value = fastp_task.outputs['out1']
        args['read2'].value = fastp_task.outputs['out2']
        args['t2'].value = top_vars['thread_number']

        # get_metrics
        depend_task = mapping
        task, args = wf.add_task(get_metrics(sample), name=f'getMetrics-{sample}', depends=[mapping.task_id])
        args['t'].value = top_vars['thread_number']
        args['ref'].value = top_vars['ref']
        args['bam'].value = depend_task.outputs['out']
        get_metrics_task_id = task.task_id

        # plot
        depend_task = task
        task, args = wf.add_task(plot_metrics(sample, method='GCBias'), name=f'plotGCBias-{sample}', depends=[get_metrics_task_id])
        task.cmd.meta.name = 'plotGCBias'
        args['i'].value = depend_task.outputs['gc_metrics']

        task, args = wf.add_task(plot_metrics(sample, method='MeanQualityByCycle'), name=f'plotMeanQual-{sample}', depends=[get_metrics_task_id])
        task.cmd.meta.name = 'plotMeanQualityByCycle'
        args['i'].value = depend_task.outputs['mq_metrics']

        task, args = wf.add_task(plot_metrics(sample, method='QualDistribution'), name=f'plotQualDistr-{sample}', depends=[get_metrics_task_id])
        task.cmd.meta.name = 'plotQualDistribution'
        args['i'].value = depend_task.outputs['qd_metrics']

        task, args = wf.add_task(plot_metrics(sample, method='InsertSizeMetricAlgo'), name=f'plotInsertSize-{sample}', depends=[get_metrics_task_id])
        task.cmd.meta.name = 'plotInsertSize'
        args['i'].value = depend_task.outputs['insert_metrics']

        # locus
        locus_get, args = wf.add_task(locus_collector(sample), name=f'locusCollector-{sample}', depends=[mapping.task_id])
        args['t'].value = top_vars['thread_number']
        args['bam'].value = mapping.outputs['out']

        # dedup
        dedup_task, args = wf.add_task(dedup(sample), name=f'dedup-{sample}', depends=[mapping.task_id])
        args['t'].value = top_vars['thread_number']
        args['bam'].value = mapping.outputs['out']
        args['score'].value = locus_get.outputs['score']

        # coverage
        cov_task, args = wf.add_task(coverage_metrics(sample), name=f'covMetrics-{sample}', depends=[dedup_task.task_id])
        args['t'].value = top_vars['thread_number']
        if intervals:
            args['intervals'].value = [top_vars['intervals']]
        args['bam'].value = dedup_task.outputs['out_bam']
        args['ref'].value = top_vars['ref']

        # realign
        if do_realign:
            realign_task, args = wf.add_task(realign(sample), name=f'realign-{sample}', depends=[dedup_task.task_id])
            args['t'].value = top_vars['thread_number']
            args['bam'].value = dedup_task.outputs['out_bam']
            args['ref'].value = top_vars['ref']
            args['database'].value = [top_vars['known_indel'], top_vars['known_mills']]

        # recalibration
        depend_task = realign_task if do_realign else dedup_task
        recal_task, args = wf.add_task(recalibration(sample), name=f'recal-{sample}', depends=[depend_task.task_id])
        args['t'].value = top_vars['thread_number']
        if intervals:
            args['intervals'].value = [top_vars['intervals']]
        args['bam'].value = depend_task.outputs['out_bam']
        args['ref'].value = top_vars['ref']
        args['database'].value = [top_vars['known_dbsnp'], top_vars['known_indel'], top_vars['known_mills']]
        recal_dict[sample] = recal_task
        bam_dict[sample] = depend_task

    if pair_info_txt is None:
        pair_info_txt = 'example_pair_info.txt'
        with open(pair_info_txt, 'w') as f:
            f.write('tumor\tnormal')

    for line in open(pair_info_txt):
        tumor_sample, normal_sample = line.split('\t')[:2]
        if tumor_sample not in bam_dict:
            print(f'Warning: tumor_sample {tumor_sample} is not in target list: {list(bam_dict.keys())}')
            continue
        if tumor_sample not in bam_dict:
            print(f'Warning: tumor_sample {tumor_sample} is not in target list: {list(bam_dict.keys())}')
            continue

        # germline variant calling
        # haplotyper
        hap_task, args = wf.add_task(Haplotyper(normal_sample), name=f'haplotyper-{normal_sample}', depends=[bam_dict[normal_sample].task_id])
        args['ref'].value = top_vars['ref']
        args['bam'].value = bam_dict[normal_sample].outputs['out_bam']
        args['recal_data'].value = recal_dict[normal_sample].outputs['recal_data']
        if intervals:
            args['intervals'].value = [top_vars['intervals']]

        # gvcf-typer
        germline_task, args = wf.add_task(GVCFtyper(normal_sample), name=f'gvcfTyper-{normal_sample}', depends=[hap_task.task_id])
        args['ref'].value = top_vars['ref']
        args['known_dbsnp'].value = top_vars['known_dbsnp']
        args['in_gvcf'].value = [hap_task.outputs['out_vcf']]

        # pair call
        task, args = wf.add_task(TNhaplotyper2(tumor_sample=tumor_sample), name=f'TNhaplotyper2-{tumor_sample}')
        task.depends = [bam_dict[normal_sample].task_id, bam_dict[tumor_sample].task_id]
        task.depends += [recal_dict[normal_sample].task_id, recal_dict[tumor_sample].task_id]
        args['ref'].value = top_vars['ref']
        args['t'].value = top_vars['thread_number']
        if intervals:
            args['intervals'].value = [top_vars['intervals']]
        args['bams'].value = [t.outputs['out_bam'] for t in bam_dict.values()]
        args['recal_datas'].value = [t.outputs['recal_data'] for t in recal_dict.values()]
        args['normal_sample'].value = normal_sample
        # pon and germline
        if pon:
            args['pon'].value = top_vars['pon']
        if germline_vcf:
            args['germline_vcf'].value = top_vars['germline_vcf']
        # orientation
        args['orientation_sample'].value = tumor_sample
        args['orientation_data'].value = f'{tumor_sample}.orientation.data'
        # contamination
        if germline_vcf:
            args['germline_vcf2'].value = top_vars['germline_vcf']
            args['contamination_tumor'].value = tumor_sample
            args['contamination_normal'].value = normal_sample
            args['contamination_data'].value = f'{tumor_sample}.contamination.data'
            args['tumor_segments'].value = f'{tumor_sample}.contamination.segments'

        # filter
        depend_task = task
        filter_task, args = wf.add_task(TNfilter(tumor_sample), name=f'tnfilter-{tumor_sample}', depends=[depend_task.task_id])
        somatic_task = task
        args['ref'].value = top_vars['ref']
        args['normal_sample'].value = normal_sample
        args['tmp_vcf'].value = depend_task.outputs['out_vcf']
        if germline_vcf:
            args['contamination'].value = depend_task.outputs['contamination_data']
            args['tumor_segments'].value = depend_task.outputs['tumor_segments']
        args['orientation_data'].value = depend_task.outputs['orientation_data']

        # annotation of somatic variant with VEP
        if os.path.exists(vep_cache):
            vep_task, args = wf.add_task(vep(tumor_sample), name=f'vep-{tumor_sample}', depends=[filter_task.task_id])
            args['input_file'].value = filter_task.outputs['out_vcf']
            args['fasta'].value = top_vars['ref']
            args['dir_cache'].value = top_vars['vep_cache_dir']
            args['dir_plugins'].value = top_vars['vep_plugin_dir']

        # ---combine germline vcf and somatic vcf---
        # combine_task, args = wf.add_task(CombineVariants(tumor_sample), depends=[germline_task.task_id, somatic_task.task_id])
        # args['ref'].value = top_vars['ref']
        # args['variant'].value = [germline_task.outputs['out_vcf'], somatic_task.outputs['out_vcf']]
        # sort_task, args = wf.add_task(SortVcf(tumor_sample), depends=[combine_task.task_id])
        # args['in_vcf'].value = combine_task.outputs['combined_vcf']
        # phase_task, args = wf.add_task(ReadBackedPhasing(tumor_sample), depends=[combine_task.task_id])
        # args['ref'].value = top_vars['ref']
        # args['variant'].value = sort_task.outputs['sorted_vcf']
        # args['interval'].value = sort_task.outputs['sorted_vcf']
        # args['bam'].value = bam_dict[tumor_sample].outputs['out_bam']
        # --- vep annot phased vcf---
        # task, args = wf.add_task(vep(tumor_sample), depends=[depend_task.task_id])
        # task.cmd.meta.name = 'vep_phased'
        # args['input_file'].value = phase_task.outputs['phased_vcf']
        # args['output_file'].value = f'{tumor_sample}.phased.vep.vcf.gz'
        # args['fasta'].value = top_vars['ref']
        # args['dir_cache'].value = top_vars['vep_cache_dir']
        # args['dir_plugins'].value = top_vars['vep_plugin_dir']

    wf.to_nestcmd(outdir=outdir, run=run, no_docker=no_docker, threads=workers, retry=retry,
                  no_monitor_resource=no_monitor_resource, no_check_resource=no_check_resource)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['pipeline'])
