import os
import re
import json
from glob import glob
import pandas as pd
import numpy as np


def find_files(query_dir, names:tuple):
    results = [[] for x in names]
    for root, dirs, files in os.walk(query_dir):
        for each in files:
            for ind, target in enumerate(names):
                if re.fullmatch(target, each):
                    results[ind].append(os.path.join(root, each))
    return results


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

    def apply_diploid_split(cell):
        if cell is None or type(cell) == float:
            return 'unknown|unknown'
        genes = [x.split(' (')[0] for x in cell.split(',')]
        abundances = list(map(float, re.findall(r'abundance: (\d+\.\d+)%', cell)))
        if len(genes) == 1:
            if abundances[0] > 30:
                return genes[0] + '|' + genes[0]
            else:
                return 'unknown|unknown'
        else:
            # print(abundances, genes)
            if abundances[0]/abundances[1] > 3 and sum(abundances[:2]) > 60:
                return genes[0] + '|' + genes[0]
            elif sum(abundances[:2]) > 0.75:
                return '|'.join(sorted(genes[:2]))
            else:
                return 'unknown|unknown'

    df2 = df[[x for x in df.columns if x.startswith('EM:')]].applymap(apply_diploid_split)
    df2.to_csv(os.path.join(outdir, 'hisat_genotype.diploid.txt'), sep='\t')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
