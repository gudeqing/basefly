import os
import re
import json
from glob import glob
import pandas as pd
import numpy as np

"""
这里包含的工具是用于结果整理的
"""


def find_files(query_dir, names:tuple, followlinks=True):
    results = [[] for x in names]
    for root, dirs, files in os.walk(query_dir, followlinks=followlinks):
        for each in files:
            for ind, target in enumerate(names):
                if re.fullmatch(target, each):
                    results[ind].append(os.path.join(root, each))
    return results


def CollectRnaSeqMetrics(files:list):
    if type(files) == str:
        files = glob(files)
        if not files:
            raise ValueError('No target file matched!')
    data = list()
    for each in files:
        if not os.path.exists(each):
            pass
        # histogram_line = [x[0] for x in enumerate(open(each)) if x[1].startswith('## HISTOGRAM')][0]
        sample = os.path.basename(each).split('.', 1)[0]
        summary = pd.read_table(each, comment='#', header=0, nrows=1)
        summary.index = [sample]
        data.append(summary)
    data = pd.concat(data, axis=0).dropna(axis=1).round(4)
    # print(data.head())
    data = data.drop(['CORRECT_STRAND_READS', 'INCORRECT_STRAND_READS', 'IGNORED_READS', 'PCT_CORRECT_STRAND_READS'], axis=1)
    # data = data.transpose()
    # out_table = os.path.join(outdir, 'RnaSeqMetrics.xls')
    # data.to_csv(out_table, index=True, header=True, sep='\t')
    return data


def CollectStarAlignmentMetrics(log_files: list, outdir=os.getcwd()):
    results = list()
    for logfile in log_files:
        sample = os.path.basename(logfile).split('.', 1)[0]
        with open(logfile) as fr:
            _ = [fr.readline() for i in range(5)]
            result = dict(sample=sample)
            for line in fr:
                if '|' in line:
                    desc, value = line.split('|')
                    desc = desc.strip()
                    value = value.strip()
                    result[desc] = value
            uniq_map = float(result['Uniquely mapped reads %'][:-1])
            multi_map = float(result['% of reads mapped to multiple loci'][:-1])
            too_many_map = float(result['% of reads mapped to too many loci'][:-1])
            result['PCT_mapped_reads'] =  (uniq_map + multi_map + too_many_map)*0.01
            results.append(result)
    df = pd.DataFrame(results).set_index('sample')
    outfile = os.path.join(outdir, 'star_alignment_stat.csv')
    df.to_csv(outfile)
    return df


def merge_metrics(query_dir, filter_ref=None, outdir=os.getcwd(), failed_cutoff=1, name_dict=''):
    if not filter_ref:
        filter_ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rna_qc_cutoff.txt')
    if name_dict:
        name_dict = dict(x.strip().split()[:2] for x in open(name_dict))
    else:
        name_dict = dict()
    rna_metric, star_metric = find_files(query_dir, (".*.RnaSeqMetrics.txt", '.*.Log.final.out'))
    print(f'we found there are {len(rna_metric)} *.RnaSeqMetrics.txt files')
    print(f'we found there are {len(rna_metric)} *.Log.final.out files')
    # rna_metric = glob(os.path.join(query_dir, '*', '*.RnaSeqMetrics.txt'))
    rna_metric = sorted(rna_metric, key=lambda x: os.path.basename(x).split('.', 1)[0])
    # star_metric = glob(os.path.join(query_dir, '*', '*.Log.final.out'))
    star_metric = sorted(star_metric, key=lambda x: os.path.basename(x).split('.', 1)[0])
    metric_dfs = list()
    if rna_metric:
        metric_dfs.append(CollectRnaSeqMetrics(rna_metric))
    if star_metric:
        os.makedirs(outdir, exist_ok=True)
        metric_dfs.append(CollectStarAlignmentMetrics(star_metric, outdir=outdir))
    if not metric_dfs:
        print('Found no qc file to merge!')
        return
    os.makedirs(outdir, exist_ok=True)
    raw = pd.concat([x.transpose() for x in metric_dfs], sort=False)
    raw.columns = [name_dict[x] if x in name_dict else x for x in raw.columns]
    raw.index.name = 'Metrics'
    raw.to_csv(os.path.join(outdir, 'MergedMetrics.csv'))

    # filter and annotate metrics
    ref = pd.read_table(filter_ref, index_col=0, header=0)
    target_metrics = [x for x in ref.index if x in raw.index]
    ref = ref.loc[target_metrics]
    data_type_dict = dict(zip(ref.index, [eval(x) for x in ref['type']]))
    out_table = raw.loc[target_metrics]
    out_table = out_table[~out_table.index.duplicated()]

    def str2float(x):
        if type(x) == str:
            if x.endswith('%'):
                return float(x[:-1])/100
            else:
                return float(x)
        else:
            return x

    out_table = out_table.transpose().applymap(str2float).round(4).astype(data_type_dict).transpose()
    lower_limit = [-np.inf if x.lower() == "none" else float(x) for x in ref['lower_limit']]
    upper_limit = [np.inf if x.lower() == "none" else float(x) for x in ref['upper_limit']]
    pass_lower = out_table.apply(lambda x: x >= lower_limit, axis=0)
    pass_upper = out_table.apply(lambda x: x <= upper_limit, axis=0)
    pass_state = pass_lower & pass_upper
    failed_samples = pass_state.columns[pass_state.sum() < ref.shape[0]-failed_cutoff+1]
    out_log = open(os.path.join(outdir, 'reason.txt'), 'w')
    print('failed sample number: {}'.format(len(failed_samples)), file=out_log)
    for sample in failed_samples:
        reason = pass_state.index[pass_lower[sample] == False]
        reason2 = pass_state.index[pass_upper[sample] == False]
        all_reason = reason.append(reason2)
        print(sample, 'failed for following reasons:', file=out_log)
        for metric in all_reason:
            metric_value = out_table.loc[metric, sample]
            limit = ref.loc[metric, ['lower_limit', 'upper_limit']]
            print('  ',metric, '=', metric_value, 'out of range [{x[0]}, {x[1]}]'.format(x=list(limit)), file=out_log)
    out_log.close()

    # save result
    out_table.loc['Pass'] = 'Yes'
    out_table.loc['Pass', pass_state.sum() < ref.shape[0]-failed_cutoff+1] = 'No'
    out_table['describe'] = list(ref['Description']) + ['Pass QC filter']
    out_table.set_index('describe', inplace=True, append=True)
    out_name = os.path.join(outdir, 'qc.all.summary.csv')
    out_table.to_csv(out_name, encoding='utf_8_sig')
    out_table.to_excel(out_name[:-3]+'xlsx')

    # save result2
    target_metrics =[
        'Number of input reads',
        'PCT_mapped_reads',
        'Average input read length',
        'Uniquely mapped reads number',
        'Uniquely mapped reads %',
        'Number of reads mapped to multiple loci',
        '% of reads mapped to multiple loci',
        'Number of chimeric reads',
        '% of chimeric reads',
        'Duplication_rate',
        'PF_ALIGNED_BASES',
        'PCT_RIBOSOMAL_BASES',
        'PCT_CODING_BASES',
        'PCT_INTERGENIC_BASES',
        'PCT_INTRONIC_BASES',
        'PCT_UTR_BASES',
        'PCT_MRNA_BASES',
        'MEDIAN_3PRIME_BIAS',
        'MEDIAN_5PRIME_BIAS'
    ]
    out_name = os.path.join(outdir, 'qc.report.summary.csv')
    hit_metrics = [x for x in target_metrics if x in out_table.index]
    out_table.loc[hit_metrics, :].to_csv(out_name, encoding='utf_8_sig')
    out_table.loc[hit_metrics, :].to_excel(out_name[:-3]+'xlsx')


def merge_salmon_quant(query_dir, outdir='.'):
    gene_files = []
    trans_files = []
    for root, dirs, files in os.walk(query_dir):
        for each in files:
            if each == 'quant.genes.sf':
                gene_files.append(os.path.join(root, each))
            if each == 'quant.sf':
                trans_files.append(os.path.join(root, each))
    # gene_files = glob(os.path.join(query_dir, 'salmon-*', 'quant.genes.sf'))
    # trans_files = glob(os.path.join(query_dir, 'salmon-*', 'quant.sf'))
    if not gene_files:
        print('no target salmon result file found')
        return
    os.makedirs(outdir, exist_ok=True)
    samples = [x.rsplit('/', 2)[1].split('-', 1)[1] for x in gene_files]
    for level, target_files in zip(['gene', 'transcript'], [gene_files, trans_files]):
        for feature in ['TPM', 'NumReads']:
            data = pd.read_table(target_files[0], index_col=0).loc[:, [feature]]
            data.columns = [samples[0]]
            for sample, file in zip(samples[1:], target_files[1:]):
                tmp = pd.read_table(file, index_col=0).loc[:, [feature]]
                tmp.columns = [sample]
                data = data.join(tmp)
            out_name = os.path.join(outdir, f'{level}.{feature}.txt')
            data.round(3).to_csv(out_name, sep='\t')


def merge_arcasHLA_genetype(query_dir, outdir='.'):
    files = find_files(query_dir, ("genotypes.tsv", ))
    print(f'we found there are {len(files[0])} genotypes.tsv files')
    if not files[0]:
        return
    os.makedirs(outdir, exist_ok=True)
    tables = []
    for each in files[0]:
        tables.append(pd.read_table(each, index_col=0))
    data = pd.concat(tables)
    out_name = os.path.join(outdir, 'HLA.genotypes.txt')
    data.to_csv(out_name, sep='\t')


def merge_star_fusion(query_dir, outdir='.'):
    files = find_files(query_dir, ("star-fusion.fusion_predictions.abridged.tsv", ))
    files = files[0]
    print(f'we found there are {len(files)} star-fusion.fusion_predictions.abridged.tsv files')
    if not files:
        return
    os.makedirs(outdir, exist_ok=True)
    samples = [x.rsplit('/', 2)[1].split('-', 1)[1] for x in files]
    tables = []
    for sample, each in zip(samples, files):
        tmp = pd.read_table(each)
        tmp['Sample'] = sample
        tables.append(tmp.set_index('Sample'))
    data = pd.concat(tables)
    out_name = os.path.join(outdir, 'star-fusion.abridged.txt')
    data.to_csv(out_name, sep='\t')


def merge_hisat_genotype(query_dir, outdir='.', name_pattern='.*.HLA-gene-type.txt'):
    tsv_files = find_files(query_dir, names=(name_pattern, ))[0]
    tables = []
    for each in tsv_files:
        tmp = pd.read_table(each, index_col=0)
        sample = os.path.basename(each).split('.', 1)[0]
        tmp.index = [sample]
        tables.append(tmp)
    df = pd.concat(tables)
    df.to_csv(os.path.join(outdir, 'hisat_genotype.raw.txt'), sep='\t')
    # 进一步将等位基因提取成2列，即按照二倍体的方式区分
    simplify_hisat_genotype(df, outdir=outdir)


def simplify_hisat_genotype(df, outdir='.'):
    if type(df) == str:
        df = pd.read_table(df, index_col=0)
    def apply_diploid_simplify(cell):
        if cell is None or type(cell) == float:
            return 'missed|missed'
        genes = [x.split()[0] for x in cell.split(',')]
        # scores = list(map(float, re.findall(r'abundance: (\d+\.\d+)%', cell)))
        scores = list(map(float, re.findall(r'score: (\d+\.\d+)', cell)))
        # print(scores)
        if len(genes) == 1:
            if sum(scores) > 0.3:
                genes = genes*2
            else:
                genes = ['lowScore']*2
        else:
            # print(scores)
            if scores[0]/scores[1] > 3 and sum(scores[:2]) > 0.75:
                genes = genes*2
            elif sum(scores[:2]) > 0.5 and (genes[0].startswith('DPB1') or genes[0].startswith('DRB1')):
                pass
            elif sum(scores[:2]) > 0.6:
                pass
            elif scores[0] > 0.45:
                genes = [genes[0]]*2
            else:
                genes = ['lowScore']*2

        return '|'.join(sorted(genes[:2]))

    # df2 = df[[x for x in df.columns if x.startswith('EM:')]].applymap(apply_diploid_simplify_with_EM)
    df2 = df[[x for x in df.columns if x.startswith('Allele splitting:')]].applymap(apply_diploid_simplify)
    df2.columns = [x.split(': ')[1] for x in df2.columns]
    df2.to_csv(os.path.join(outdir, 'hisat_genotype.diploid.txt'), sep='\t')
    # stat distribution


def get_4digits_hla_genetype(table, sample, alleles=('A', 'B', 'C', 'DRA', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1', 'DPB2')):
    """prepare for pvacseq tool based on result of merge_hisat_genotype"""
    df = pd.read_table(table, index_col=0, header=0)
    targets = list()
    for each in df.loc[sample]:
        try:
            a, b = each.split('|')[:2]
        except Exception as e:
            raise Exception('???', sample, each, e)
        if a.startswith(alleles):
            allele_2_digit = ':'.join(a.split(':')[:2])
            if allele_2_digit not in targets:
                if allele_2_digit.startswith(('A*', 'B*', 'C*', 'E*', 'F*', 'G*')):
                    allele_2_digit = 'HLA-' + allele_2_digit
                targets.append(allele_2_digit)
        if b.startswith(alleles):
            allele_2_digit = ':'.join(b.split(':')[:2])
            if allele_2_digit not in targets:
                if allele_2_digit.startswith(('A*', 'B*', 'C*', 'E*', 'F*', 'G*')):
                    allele_2_digit = 'HLA-' + allele_2_digit
                targets.append(allele_2_digit)
    return sorted(set(targets))


def merge_epitopes(query_dir, outdir='.'):
    """
    针对pvacseq的结果进行合并
    """
    patterns = ('.*.filtered.tsv', )
    files_list = find_files(query_dir, patterns)
    os.makedirs(outdir, exist_ok=True)
    for pattern, files in zip(patterns, files_list):
        print(f'we found there are {len(files)} {pattern} files')
        if not files:
            return
        tables = []
        for each in files:
            if '/MHC_Class_I/' in each:
                sample = os.path.basename(each).split('.filter')[0]
                tmp = pd.read_table(each)
                tmp['Sample'] = sample
                tables.append(tmp.set_index('Sample'))
        data = pd.concat(tables)
        out_name = os.path.join(outdir, 'MHC_I.merged.epitopes.txt')
        data.to_csv(out_name, sep='\t')

        tables = []
        for each in files:
            if '/MHC_Class_II/' in each:
                sample = os.path.basename(each).split('.filter')[0]
                tmp = pd.read_table(each)
                tmp['Sample'] = sample
                tables.append(tmp.set_index('Sample'))
        data = pd.concat(tables)
        out_name = os.path.join(outdir, 'MHC_II.merged.epitopes.txt')
        data.to_csv(out_name, sep='\t')

        tables = []
        for each in files:
            if '/combined/' in each:
                sample = os.path.basename(each).split('.filter')[0]
                tmp = pd.read_table(each)
                tmp['Sample'] = sample
                tables.append(tmp.set_index('Sample'))
        data = pd.concat(tables)
        out_name = os.path.join(outdir, 'All.merged.epitopes.txt')
        data.to_csv(out_name, sep='\t')


def merge_mhcflurry_result(query_dir, outdir='.', min_exp=0.5):
    patterns = ('.*.annotated.mhcflurry.csv',)
    files_list = find_files(query_dir, patterns)
    os.makedirs(outdir, exist_ok=True)
    for pattern, files in zip(patterns, files_list):
        print(f'we found there are {len(files)} {pattern} files')
        if not files:
            return
        tables = []
        for each in files:
            sample = os.path.basename(each).split('.annotated')[0]
            tmp = pd.read_csv(each, header=0)
            tmp = tmp.loc[tmp['TPM'] >= min_exp]
            tmp['Sample'] = sample
            tables.append(tmp.set_index('Sample'))
        data = pd.concat(tables).round(3)
        out_name = os.path.join(outdir, 'MHC_I.merged.epitopes.csv')
        data.to_csv(out_name)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
