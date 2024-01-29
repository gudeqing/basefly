from basefly import Argument, Output, Command, Workflow, TopVar

"""
HLA-typing tool benchmark:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8045758/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8253103/:
 * POLYSOLVER, OpiType, and xHLA ranked as top three HLA class I 
 * HLA‐HD, SOAP‐HLA, and xHLA obtain the relatively high performance (accuracy 0.708–0.904) on HLA class II alleles DRB1, DPB1, DQA1, and DQB1 calling.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9263794/

https://www.frontiersin.org/articles/10.3389/fimmu.2022.987655/full
https://www.frontiersin.org/articles/10.3389/fimmu.2021.688183/full
OptiType, HLA-HD, PHLAT, and HLA*LA showed high correction scores.
The complementary scores between OptiType and HLA-HD, HLA-HD and PHLAT, HLA-HD and HLAscan, HLA-HD and HLA*LA were high. 
OptiType, HLA-HD, and PHLAT showed high complementary ratios. 
We additionally calculated the accuracy of the combination tools considering the complementary ratio of each tool and the accuracy of each tool,
and the accuracy was over 98% in all groups with six or more concordant calls.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9679531/
On average, across the five HLA loci examined, HLA*LA was found to have the highest typing accuracy. For the individual loci, HLA-A, -B and -C, Optitype’s typing accuracy was the highest and HLA*LA had the highest typing accuracy for HLA-DRB1 and -DQB1.

https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-023-09351-z
The highest MHC-I calling accuracies were found for Optitype (98.0%) and arcasHLA (99.4%) on WES and RNA sequencing data respectively, 
while for MHC-II HLA-HD was the most accurate tool for both data types (96.2% and 99.4% on WES and RNA data respectively).
HLA*LA and HLA-HD are the best performing MHC class II genotyping tools on WES data

https://zhuanlan.zhihu.com/p/474124081
我们首先使用 HLA-HD、Athlates 和 OptiType 对 HLA 的 I 型和 II 型共计 12 个标准品分析，具体结果见表 3。
HLA-HD 和 Athlates 比 OptiType 的 HLA I 型分型更加准确，两者仅在 A 基因的 24 个 allele 上错误分型 1 个位点，而 OptiType 则在 A 和 C 基因上分型错误了 2 个和 1 个位点。
而 HLA 的 II 型分型的结果来看 HLA-HD 和 Athlates 各有优劣，HLA-HD 在 DRB3/4/5 基因上分型上错误了 4 个位点，而 Athlates 则在 DPB1 基因上效果不佳，错误了 3 个位点的分型。
进一步分析 HLA-HD 分型不准确的原因，我们发现主要原因在于使用 IMGT 数据库版本的差异。以 HLA I 型的 C1-106 样本为例，A 基因的一个 Allele 位点跟标准品的 HLA-A*02:01 不一致，而是 HLA-A*02:642。
但当使用更早的 IMGT 数据库版本（3.32.0）分析时，HLA-HD 的分型的结果为 HLA-A*02:01:01,与标准品的分型相一致。
实际上 HLA-A*02:642 在 IMGT 3.37 版中属于 HLA-A*02:01:01G group，也就是说 HLA-A*02:642 分型是正确的，只不过由于 IMGT 数据库新增分型产生了命名差异。存在此现象的还包括 C2-104、C2-112 样品，具体分型结果和说明见表 4。

HLA的数据库
https://www.ebi.ac.uk/ipd/imgt/hla/download/
https://github.com/ANHIG/IMGTHLA/

"""


def hisat_genotype():
    """
    hisatgenotype --base hla --locus-list A,B,C,DRB1,DQA1 -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz
    hisatgenotype_toolkit parse-results --csv --in-dir hisatgenotype_out
    """
    cmd = Command()
    cmd.meta.version = "1.3.3"
    cmd.meta.source = "https://github.com/DaehwanKimLab/hisat-genotype?"
    cmd.meta.name = 'HisatGenotype'
    cmd.meta.function = "HLA-typing using hisat"
    cmd.meta.desc = """
    HISAT-genotype is a next-generation genomic analysis software platform capable of assembling and genotyping human genes and genomic regions.
    The software leverages HISAT2s graph FM index and graph alignemnt algorithm to align reads to a specially constructed graph genome.
    An Expectation-Maximization (EM) algorithm finds the maximum likelihood estimates for each gene allele and a guided de Bruijn graph is used to construct the allele sequences.
    # prepare indices for hisat-genotype
    #    wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/genotype_genome_20180128.tar.gz
    #    tar xvzf genotype_genome_20180128.tar.gz
    #    rm genotype_genome_20180128.tar.gz
    #
    #    #grch38
    #    wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
    #    tar xvzf grch38.tar.gz
    #    rm grch38.tar.gz
    #    hisat2-inspect grch38/genome > genome.fa
    #    samtools faidx genome.fa
    #
    #    #HISATgenotpye Database
    #    git clone https://github.com/DaehwanKimLab/hisatgenotype_db.git
    其他用于建索引的参考：https://hxbiogit.asia/DrGBL/hisat_genotype_qgp
    
    """
    cmd.runtime.image = 'gudeqing/dnaseq:1.0'
    cmd.runtime.tool = 'hisatgenotype'
    cmd.runtime.cpu = 4
    cmd.runtime.memory = 8*1024**3
    cmd.args['base'] = Argument(prefix='--base ', default='hla', desc='Base file name for index, variants, haplotypes, etc. (e.g. hla, rbg, codis)')
    cmd.args['locus'] = Argument(prefix='--locus-list ', level='optional', array=True, delimiter=',', desc='A comma-separated list of gene names (default: empty, all genes)')
    cmd.args['read1'] = Argument(prefix='-1 ', type='infile', desc='absolute read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-2 ', type='infile', desc='absolute read2 fastq file')
    cmd.args['single_end'] = Argument(prefix='-U ', type='infile', level='optional', desc='filename for single-end reads')
    cmd.args['_read_dir'] = Argument(prefix='--in-dir ', value='/', type='fix')
    cmd.args['threads'] = Argument(prefix='--threads ', default=cmd.runtime.cpu, desc='Number of threads')
    cmd.args['hisat_threads'] = Argument(prefix='--pp ', default=4, desc='Number of threads')
    cmd.args['index_dir'] = Argument(prefix='--index_dir ', level='optional', type='indir', desc="Set location to use for indicies. you may prepared it by git clone https://github.com/DaehwanKimLab/hisatgenotype_db.git")
    cmd.args['_outdir'] = Argument(prefix='--out-dir ', value='./', type='fix')
    cmd.args['_parse_result'] = Argument(value='&& hisatgenotype_toolkit parse-results --csv --in-dir .', type='fix')
    cmd.args['level'] = Argument(prefix='-t ', default=2, desc='Trim allele to specific field level (example : A*01:01:01:01 trim 2 A*01:01)')
    cmd.args['out'] = Argument(prefix='--output-file ', default='HLA-gene-type.txt', desc='output of csv file')
    cmd.outputs['out'] = Output(value='{out}', type='outfile')
    return cmd


if __name__ == '__main__':
    # Workflow().to_cwl_tool(cmd=hisat_genotype(), write_out=True)
    print("目前我下载不到references")
    hisat_genotype().run_on_terminal()

