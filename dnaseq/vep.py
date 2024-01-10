import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys;
sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar

__author__ = 'gdq'


def vep():
    cmd = Command()
    cmd.meta.version = 'vep104.3'
    cmd.meta.name = 'VEP'
    cmd.meta.function = 'Variant Annotation'
    cmd.meta.source = 'https://grch37.ensembl.org/info/docs/tools/vep/index.html'
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
    cmd.runtime.image = 'ensemblorg/ensembl-vep:release_110.1'
    cmd.runtime.tool = 'vep'
    cmd.runtime.docker_cmd_prefix = cmd.runtime.docker_cmd_prefix2
    cmd.args['input_file'] = Argument(prefix='-i ', type='infile', desc='input file')
    cmd.args['fasta'] = Argument(prefix='--fasta ', type='infile', default='/home/hxbio04/hg19/genome.fa',
                                 desc="Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run VEP with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache).")
    cmd.args['output_file'] = Argument(prefix='-o ', type='outstr', default='vep.vcf.gz', desc='output file')
    cmd.args['output_format'] = Argument(prefix='--', type='outstr', range={'vcf', 'json', 'tab'}, default='vcf',
                                         desc="If we choose to write output in VCF format. Consequences are added in the INFO field of the VCF file, using the key 'CSQ'. Data fields are encoded separated by '|'; the order of fields is written in the VCF header. Output fields in the 'CSQ' INFO field can be selected by using --fields.")
    cmd.args['compress_output'] = Argument(prefix='--compress_output ', default='bgzip',
                                           desc="Writes output compressed using either gzip or bgzip")
    cmd.args['force_overwrite'] = Argument(prefix="--force_overwrite ", type='bool', default=True,
                                           desc="Force overwriting of output file")
    cmd.args['fork'] = Argument(prefix='--fork ', type='int', default=4,
                                desc='Use forking(multi-cpu/threads) to improve script runtime')
    cmd.args['species'] = Argument(prefix='--species ', default='homo_sapiens',
                                   desc='Species for your data. This can be the latin name e.g. homo_sapiens or any Ensembl alias e.g. mouse.')
    cmd.args['assembly_version'] = Argument(prefix='--assembly ', default='GRCh38',
                                            desc='Select the assembly version to use if more than one available.')
    cmd.args['dir_cache'] = Argument(prefix='--dir_cache ', type='indir', default='/home/hxbio04/dbs/vep', desc='Specify the cache directory to use')
    cmd.args['dir_plugins'] = Argument(prefix='--dir_plugins ', type='indir', level='optional',
                                       desc='Specify the plugin directory to use')
    cmd.args['merged'] = Argument(prefix='--merged ', type='bool', default=False,
                                  desc='Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.')
    cmd.args['refseq'] = Argument(prefix='--refseq ', type='bool', default=False,
                                  desc='Specify this option if you have installed the RefSeq cache in order for VEP to pick up the alternate cache directory. This cache contains transcript objects corresponding to RefSeq transcripts. Consequence output will be given relative to these transcripts in place of the default Ensembl transcripts')
    cmd.args['stats_file'] = Argument(prefix='--stats_file ', type='outstr', default='vep.summary.html',
                                      desc='Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>.')
    cmd.args['cache'] = Argument(prefix='--cache ', type='bool', default=True, desc='Enables use of cache')
    cmd.args['offline'] = Argument(prefix='--offline ', type='bool', default=True,
                                   desc='Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation')
    cmd.args['merged'] = Argument(prefix='--merged ', type='bool', default=False,
                                  desc='Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.')
    cmd.args['plugins'] = Argument(prefix='--plugin ', multi_times=True, level='optional',
                                   desc='Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times')
    cmd.args['variant_class'] = Argument(prefix='--variant_class ', type='bool', default=True,
                                         desc='Output the Sequence Ontology variant class.')
    cmd.args['sift'] = Argument(prefix='--sift ', default='b', range={'p', 's', 'b'},
                                desc="Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both.")
    cmd.args['polyphen'] = Argument(prefix='--polyphen ', default='b', range={'p', 's', 'b'},
                                    desc="Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both.")
    cmd.args['nearest'] = Argument(prefix='--nearest ', default='transcript', range={'transcript', 'gene', 'symbol'},
                                   desc='Retrieve the transcript or gene with the nearest protein-coding transcription start site (TSS) to each input variant. Use transcript to retrieve the transcript stable ID, gene to retrieve the gene stable ID, or symbol to retrieve the gene symbol. Note that the nearest TSS may not belong to a transcript that overlaps the input variant, and more than one may be reported in the case where two are equidistant from the input coordinates.')
    cmd.args['gene_phenotype'] = Argument(prefix='--gene_phenotype ', type='bool', default=True,
                                          desc='Indicates if the overlapped gene is associated with a phenotype, disease or trait.')
    cmd.args['regulatory'] = Argument(prefix='--regulatory ', type='bool', default=True,
                                      desc="Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature.")
    cmd.args['phased'] = Argument(prefix='--phased ', type='bool', default=True,
                                  desc="Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data.")
    cmd.args['numbers'] = Argument(prefix='--numbers ', type='bool', default=True,
                                   desc="Adds affected exon and intron numbering to to output. Format is Number/Total")
    cmd.args['hgvs'] = Argument(prefix='--hgvs ', type='bool', default=True,
                                desc="Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate.")
    cmd.args['transcript_version'] = Argument(prefix='--transcript_version ', type='bool', default=True,
                                              desc="Add version numbers to Ensembl transcript identifiers")
    cmd.args['symbol'] = Argument(prefix='--symbol ', type='bool', default=True,
                                  desc="Adds the gene symbol (e.g. HGNC) (where available) to the output.")
    cmd.args['tsl'] = Argument(prefix='--tsl ', type='bool', default=True,
                               desc="Adds the transcript support level for this transcript to the output.")
    cmd.args['canonical'] = Argument(prefix='--canonical ', type='bool', default=True,
                                     desc="Adds a flag indicating if the transcript is the canonical transcript for the gene")
    cmd.args['biotype'] = Argument(prefix='--biotype ', type='bool', default=True,
                                   desc="Adds the biotype of the transcript or regulatory feature.")
    cmd.args['max_af'] = Argument(prefix='--max_af ', type='bool', default=True,
                                  desc="Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD")
    cmd.args['af_1kg'] = Argument(prefix='--af_1kg ', type='bool', default=True,
                                  desc="Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output.")
    cmd.args['af_gnomad'] = Argument(prefix='--af_gnomad ', type='bool', default=True,
                                     desc="Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included")
    cmd.args['af_esp'] = Argument(prefix='--af_esp ', type='bool', default=False,
                                  desc="Include allele frequency from NHLBI-ESP populations.")
    cmd.args['coding_only'] = Argument(prefix='--coding_only ', type='bool', default=False,
                                       desc="Only return consequences that fall in the coding regions of transcripts. Not used by default")
    cmd.args['pick'] = Argument(prefix='--pick', type='bool', default=False,
                                desc="Pick one line or block of consequence data per variant, including transcript-specific columns. This is the best method to use if you are interested only in one consequence per variant")
    cmd.args['flag_pick'] = Argument(prefix='--flag_pick ', type='bool', default=True,
                                     desc="As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others.")
    cmd.args['filter_common'] = Argument(prefix='--filter_common ', type='bool', default=False,
                                         desc="Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01. May be modified using any of the following freq_* filters.")
    cmd.args['other_args'] = Argument(default='', desc='specify other arguments that you want to append to the command')
    cmd.args['_create_index'] = Argument(value='&& if [ -f *.vcf.gz ]; then tabix *vcf.gz; fi', type='fix')
    cmd.outputs['out_vcf'] = Output(value='{output_file}', report=True)
    cmd.outputs['out_vcf_idx'] = Output(value='{output_file}.tbi', report=True)
    cmd.outputs['out_stats'] = Output(value='{stats_file}', report=True)
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.version = "vep104.3"
    wf.meta.name = 'Variant-Annotation-Workflow'
    wf.meta.source = 'https://grch37.ensembl.org/info/docs/tools/vep/index.html'
    wf.meta.function = 'Variant Annotation'
    wf.meta.desc = """
    主要功能：
    当前流程是突变注释分析流程，使用的工具是VEP（Variant Effect Predictor），其是一款强大的基因变异注释软件，由Ensembl开发。它支持多种基因组版本，可以识别和注释各种类型的变异，如单核苷酸变异、插入/删除、结构变异等。
    使用VEP进行基因变异注释的步骤如下：
    * 准备数据：首先需要准备基因组变异数据，可以是VCF格式或 Variant Call Format (VCF) 格式。
    * 安装VEP：根据操作系统和所需的数据类型，从VEP官方网站下载并安装适合的软件版本。
    * 运行VEP：使用命令行或脚本方式运行VEP，指定输入数据文件和参考基因组文件。
    * 参数设置：根据需求选择所需的注释选项和参数，例如是否进行基因表达注释、甲基化注释等。参数说明可以参考页面：https://grch37.ensembl.org/info/docs/tools/vep/script/vep_options.html
    * 结果输出：VEP将注释结果输出到标准输出或指定的结果文件中，包括变异的类型、位置、功能影响等信息。
    通过使用VEP，研究人员可以更全面地了解变异的生物学意义，从而为疾病研究、药物研发等方面提供有力支持。

    使用示例:
    * 确认已经准备好输入文件或加载所需数据集（如流程文件数据集，样本数据集，参考文件数据集）
    * 运行命令: python /enigma/datasets/*/scripts/vep/vep.py -i /enigma/datasets/*/your_testdata_dir/*vcf -dir_cache /enigma/datasets/*/vep -assembly GRCh38 --refseq --plot --run --docker -outdir /enigma/local_storage/Result
    * 分析结果：主要分析结果将汇总于输出目录如Result/Report

    主要输入说明（以hg38为例）:
    当前支持的VEP版本为104.3，需要下载相应的注释数据库：
    hg37版：https://ftp.ensembl.org/pub/grch37/release-104/variation/indexed_vep_cache/
    如：wget http://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/homo_sapiens_refseq_vep_104_GRCh37.tar.gz
    hg38版：http://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/
    如：wget http://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/homo_sapiens_vep_104_GRCh38.tar.gz
    
    关键参数说明：
    1. 如果使用的数据库版本是refseq，需加参数--refseq，如果使用的是merged版本的，则需加参数--merged
    2. assembly版本，默认是 --assembly GRCh38
    3. 通过参数input_files指定输入1个或多个文件，需vcf格式
    4. 通过output_format参数指定输出文件格式：{'vcf', 'json', 'tab'}，默认是vcf
    """
    wf.init_argparser()
    wf.add_argument('-input_files', nargs='+', help='input vcf files. one or more files are supported, and separated by white space')
    wf.add_argument('-output_format', default='vcf',  help="Output format, one of ['vcf', 'json', 'tab']")
    wf.add_argument('-assembly', default='GRCh38', help='Select the assembly version to use if more than one available.')
    wf.add_argument('-dir_cache', default='/enigma/datasets/*/vep/', help='Specify the cache directory to use')
    wf.add_argument('--merged', action='store_true', default=False, help='Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.')
    wf.add_argument('--refseq', action='store_true', default=False, help='Specify this option if you have installed the RefSeq cache in order for VEP to pick up the alternate cache directory. This cache contains; transcript objects corresponding to RefSeq transcripts. Consequence output will be given relative to these transcripts in place of the default Ensembl transcripts ')
    wf.parse_args()
    for each in wf.args.input_files:
        base_name = os.path.basename(each).rsplit('.vcf', 1)[0]
        if wf.args.output_format == 'vcf':
            out_file = base_name + '.vep.vcf.gz'
        elif wf.args.output_format == 'json':
            out_file = base_name + '.json'
        else:
            out_file = base_name + '.txt'
        vep_task, args = wf.add_task(vep(), tag=base_name)
        args['input_file'].value = each
        out_stat = os.path.basename(each).rsplit('.vcf', 1)[0] + '.vep.html'
        args['output_file'].value = out_file
        args['stats_file'].value = out_stat
        args['output_format'].value = wf.args.output_format
        args['assembly_version'].value = wf.args.assembly
        args['dir_cache'].value = wf.args.dir_cache
        args['merged'].value = wf.args.merged
        args['refseq'].value = wf.args.refseq
    wf.run()


if __name__ == '__main__':
    pipeline()
