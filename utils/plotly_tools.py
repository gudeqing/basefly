#! /data/users/dqgu/anaconda3/bin/python
import os
import re
from functools import partial
from collections import OrderedDict
from glob import glob
import pandas as pd
import numpy as np
import json
import textwrap
from plotly import tools
import plotly.graph_objs as go
from plotly.offline import plot as plt
import plotly.io as pio
plt = partial(plt, auto_open=False)


def get_color_pool(n):
    # https://plot.ly/ipython-notebooks/color-scales/
    import colorlover
    if n <= 8:
        if n <= 3:
            n = 3
        return colorlover.scales[str(n)]['qual']['Set2']
    if n <= 12:
        return colorlover.scales[str(n)]['qual']['Set3']

    import random
    random.seed(666)

    def get_random_color(pastel_factor=0.5):
        return [(x + pastel_factor) / (1.0 + pastel_factor) for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]]

    def color_distance(c1, c2):
        return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])

    def generate_new_color(existing_colors, pastel_factor=0.5):
        max_distance = None
        best_color = None
        for i in range(0, 100):
            color = get_random_color(pastel_factor=pastel_factor)
            # exclude some colors
            if np.absolute(np.array(color) - np.array([1, 1, 1])).sum() < 0.08:
                continue
            if not existing_colors:
                return color
            best_distance = min([color_distance(color, c) for c in existing_colors])
            if not max_distance or best_distance > max_distance:
                max_distance = best_distance
                best_color = color
        return best_color

    color_pool = []
    for i in range(0, n):
        color_pool.append(generate_new_color(color_pool, pastel_factor=0.9))
    color_pool = [(int(x * 255), int(y * 255), int(z * 255)) for x, y, z in color_pool]
    color_pool = sorted(color_pool, key=lambda x: (x[0], x[1], x[2]))
    return colorlover.to_rgb(color_pool)


def get_marker_pool(n):
    maker_pool = [
        0, 'circle', 100, 'circle-open', 200, 'circle-dot', 300,
        'circle-open-dot', 1, 'square', 101, 'square-open', 201,
        'square-dot', 301, 'square-open-dot', 2, 'diamond', 102,
        'diamond-open', 202, 'diamond-dot', 302,
        'diamond-open-dot', 3, 'cross', 103, 'cross-open', 203,
        'cross-dot', 303, 'cross-open-dot', 4, 'x', 104, 'x-open',
        204, 'x-dot', 304, 'x-open-dot', 5, 'triangle-up', 105,
        'triangle-up-open', 205, 'triangle-up-dot', 305,
        'triangle-up-open-dot', 6, 'triangle-down', 106,
        'triangle-down-open', 206, 'triangle-down-dot', 306,
        'triangle-down-open-dot', 7, 'triangle-left', 107,
        'triangle-left-open', 207, 'triangle-left-dot', 307,
        'triangle-left-open-dot', 8, 'triangle-right', 108,
        'triangle-right-open', 208, 'triangle-right-dot', 308,
        'triangle-right-open-dot', 9, 'triangle-ne', 109,
        'triangle-ne-open', 209, 'triangle-ne-dot', 309,
        'triangle-ne-open-dot', 10, 'triangle-se', 110,
        'triangle-se-open', 210, 'triangle-se-dot', 310,
        'triangle-se-open-dot', 11, 'triangle-sw', 111,
        'triangle-sw-open', 211, 'triangle-sw-dot', 311,
        'triangle-sw-open-dot', 12, 'triangle-nw', 112,
        'triangle-nw-open', 212, 'triangle-nw-dot', 312,
        'triangle-nw-open-dot', 13, 'pentagon', 113,
        'pentagon-open', 213, 'pentagon-dot', 313,
        'pentagon-open-dot', 14, 'hexagon', 114, 'hexagon-open',
        214, 'hexagon-dot', 314, 'hexagon-open-dot', 15,
        'hexagon2', 115, 'hexagon2-open', 215, 'hexagon2-dot',
        315, 'hexagon2-open-dot', 16, 'octagon', 116,
        'octagon-open', 216, 'octagon-dot', 316,
        'octagon-open-dot', 17, 'star', 117, 'star-open', 217,
        'star-dot', 317, 'star-open-dot', 18, 'hexagram', 118,
        'hexagram-open', 218, 'hexagram-dot', 318,
        'hexagram-open-dot', 19, 'star-triangle-up', 119,
        'star-triangle-up-open', 219, 'star-triangle-up-dot', 319,
        'star-triangle-up-open-dot', 20, 'star-triangle-down',
        120, 'star-triangle-down-open', 220,
        'star-triangle-down-dot', 320,
        'star-triangle-down-open-dot', 21, 'star-square', 121,
        'star-square-open', 221, 'star-square-dot', 321,
        'star-square-open-dot', 22, 'star-diamond', 122,
        'star-diamond-open', 222, 'star-diamond-dot', 322,
        'star-diamond-open-dot', 23, 'diamond-tall', 123,
        'diamond-tall-open', 223, 'diamond-tall-dot', 323,
        'diamond-tall-open-dot', 24, 'diamond-wide', 124,
        'diamond-wide-open', 224, 'diamond-wide-dot', 324,
        'diamond-wide-open-dot', 25, 'hourglass', 125,
        'hourglass-open', 26, 'bowtie', 126, 'bowtie-open', 27,
        'circle-cross', 127, 'circle-cross-open', 28, 'circle-x',
        128, 'circle-x-open', 29, 'square-cross', 129,
        'square-cross-open', 30, 'square-x', 130, 'square-x-open',
        31, 'diamond-cross', 131, 'diamond-cross-open', 32,
        'diamond-x', 132, 'diamond-x-open', 33, 'cross-thin', 133,
        'cross-thin-open', 34, 'x-thin', 134, 'x-thin-open', 35,
        'asterisk', 135, 'asterisk-open', 36, 'hash', 136,
        'hash-open', 236, 'hash-dot', 336, 'hash-open-dot', 37,
        'y-up', 137, 'y-up-open', 38, 'y-down', 138,
        'y-down-open', 39, 'y-left', 139, 'y-left-open', 40,
        'y-right', 140, 'y-right-open', 41, 'line-ew', 141,
        'line-ew-open', 42, 'line-ns', 142, 'line-ns-open', 43,
        'line-ne', 143, 'line-ne-open', 44, 'line-nw', 144,
        'line-nw-open'
    ]
    return sorted([x for x in maker_pool if type(x) == int])[: n]


def draw(fig: go.Figure, prefix='', outdir=os.getcwd(), formats=('html',),
         height:int=None, width:int=None, scale=3, desc=None):
    for format in formats:
        out_name = os.path.join(outdir, '{}.{}'.format(prefix, format))
        if format == 'html':
            plt(fig, filename=out_name)
        else:
            if format in ['svg', 'pdf']:
                scale = 1
            pio.write_image(fig, format=format, file=out_name, height=height, width=width, scale=scale)
        if desc:
            with open(out_name+'.describe.txt', 'w') as fw:
                fw.write(desc+'\n')


def gene_body_coverage(files:list, outdir=os.getcwd(), file_from='RSeQC', desc=None,
                       formats=('html',), height:int=None, width:int=None, scale=3):
    layout = go.Layout(title="geneBodyCoverage")
    all_fig = go.Figure(layout=layout)
    files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        fig = go.Figure(layout=go.Layout(title='{}'.format(sample)))
        if file_from == 'RSeQC':
            data = pd.read_table(each, header=None, index_col=0)
            normalize_y = data.iloc[1, :]/data.iloc[1, :].max()
        else:
            # assume data from Picard
            data = pd.read_table(each, header=None, skiprows=11)
            data = data.transpose()
            normalize_y = data.iloc[1, :]
        fig.add_scatter(x=data.iloc[0, :], y=normalize_y, name=sample)
        all_fig.add_scatter(x=data.iloc[0, :], y=normalize_y, name=sample)
        prefix = '{}.geneBodyCoverage'.format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)
    prefix = '{}.geneBodyCoverage'.format('samples')
    draw(all_fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def fragment_length(files:list, outdir=os.getcwd(), min_len=50, max_len=600, desc=None,
                    formats=('html',), height:int=None, width:int=None, scale=3):
    layout = go.Layout(
        title="Fragment length distribution",
        xaxis=dict(title='Fragment length'),
        yaxis=dict(title='Probability'),
    )
    all_fig = go.Figure(layout=layout)
    files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=0, index_col=0)
        data = data.loc[(data['frag_median'] >= min_len) & (data['frag_median'] <= max_len)]
        fig = go.Figure(layout=go.Layout(
            title='{}'.format(sample),
            xaxis=dict(title='Fragment length'),
            yaxis=dict(title='probability'),
        ))
        fig.add_histogram(x=data["frag_median"], histnorm='probability', name=sample)
        all_fig.add_histogram(x=data["frag_median"], histnorm='probability', name=sample)
        prefix = "{}.fragmentLengthDistribution".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)
    prefix = '{}.fragmentLengthDistribution'.format('samples')
    draw(all_fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def inner_distance(files:list, outdir, min_dist=-250, max_dist=250, desc=None,
                   formats=('html',), height:int=None, width:int=None, scale=3):
    """
    ??????1000000???????????????????????????????????????????????????
    readPairOverlap
    sameTranscript=No,dist=genomic
    sameTranscript=Yes,nonExonic=Yes,dist=genomic
    sameTranscript=Yes,sameExon=No,dist=mRNA
    sameTranscript=Yes,sameExon=Yes,dist=mRNA
    """
    groups = [
        'sameTranscript=No,dist=genomic',
        'sameTranscript=Yes,nonExonic=Yes,dist=genomic',
        'readPairOverlap',
        'sameTranscript=Yes,sameExon=No,dist=mRNA',
        'sameTranscript=Yes,sameExon=Yes,dist=mRNA',
    ]
    layout = go.Layout(
        title="Inner distance distribution",
        xaxis=dict(title='Inner Distance'),
        yaxis=dict(title='Probability'),
    )
    all_fig = go.Figure(layout=layout)
    files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=None, index_col=0)
        values = [data[2][data[2] == x].shape[0] for x in groups]
        norm_values = [x/sum(values) for x in values]
        percents = dict(zip(groups, norm_values))
        data = data[(data[1] >= min_dist) & (data[1] <= max_dist)]
        fig = go.Figure(layout=go.Layout(
                title='{}'.format(sample),
                xaxis=dict(title='Inner Distance'),
                yaxis=dict(title='Frequency'),
        ))
        for g in groups:
            target = data[1][data[2] == g]
            name = "{}({:.2%})".format(g, percents[g])
            fig.add_histogram(x=target, name=name)
        all_fig.add_histogram(x=data[1][data[2] != "sameTranscript=No,dist=genomic"], histnorm='probability', name=sample)
        prefix = "{}.InnerDistanceDistribution".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)
    prefix = "{}.InnerDistanceDistribution".format('samples')
    draw(all_fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def read_distribution(files:list, outdir, formats=('html',), height:int=None, width:int=None, scale=3, name_dict='', desc=None):
    if name_dict:
        name_dict = dict(x.strip().split()[:2] for x in open(name_dict))
    all_data = list()
    files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        if sample in name_dict:
            sample = name_dict[sample]
        data = pd.read_table(each, index_col=0, skiprows=4, skipfooter=1, sep='\s+', engine='python')
        data = data.loc[['CDS_Exons', "5'UTR_Exons", "3'UTR_Exons", "Introns", "TSS_up_10kb", "TES_down_10kb"], 'Tag_count']
        data.name = sample
        all_data.append(data)
        fig = go.Figure(layout=go.Layout(
            title='{}'.format(sample),
        ))
        fig.add_pie(labels=data.index, values=data)
        prefix = "{}.ReadDistribution".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)
    df = pd.concat(all_data, axis=1).T
    order = sorted(df.index)
    df = df.loc[order]
    # print(df.head())
    data = [go.Bar(x=df.index, y=df[x]/df.sum(axis=1), name=x) for x in df.columns]
    layout = go.Layout(
        title="Read distribution",
        # xaxis=dict(title='Sample'),
        barmode='stack'
    )
    fig = go.Figure(data=data, layout=layout)
    prefix = "{}.ReadDistribution".format('samples')
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def read_duplication(files:list, outdir=os.getcwd(), max_dup=500,
                     formats=('html',), height:int=None, width:int=None, scale=3, desc=None):
    traces = list()
    duplication_rates = list()
    ordered_files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    for each in ordered_files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=0, index_col=None)
        duplication_rate = 1 - data['UniqReadNumber'][0]/data['UniqReadNumber'].sum()
        dup_info = dict()
        dup_info['sample'] = sample
        dup_info['duplication_rate'] = duplication_rate
        duplication_rates.append(dup_info)
        data = data[data.iloc[:, 0] <= max_dup]
        trace = go.Scatter(x=data.iloc[:, 0], y=data.iloc[:, 1], name=sample, mode='markers')
        traces.append(trace)
        layout = go.Layout(
            title=sample+f' duplication rate: {duplication_rate:.2f}',
            xaxis=dict(title='Occurrence'),
            yaxis=dict(title='UniqReadNumber', type='log'),
        )
        fig = go.Figure([trace], layout=layout)
        prefix = "{}.ReadDuplication".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)

    layout = go.Layout(
        title="Read duplication",
        xaxis=dict(title='Occurrence'),
        yaxis=dict(title='UniqReadNumber', type='log'),
    )
    fig = go.Figure(traces, layout=layout)
    prefix = "{}.ReadDuplication".format('samples')
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)
    dup_rate = pd.DataFrame(duplication_rates).set_index('sample')
    return dup_rate

def exp_saturation(files:list, outdir=os.getcwd(), outlier_limit=5, exp_lower=0.5,
                   formats=('html',), height:int=None, width:int=None, scale=3, desc=None):
    all_fig = tools.make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Q1',
            'Q2',
            'Q3',
            'Q4',
        ),
        shared_xaxes=False
    )
    files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    color_pool = get_color_pool(len(files))
    for exp_file, sample_color in zip(files, color_pool):
        sample = os.path.basename(exp_file).split('.', 1)[0]
        data = pd.read_table(exp_file, header=0, index_col=0)
        # plot deviation
        data = data[data['100'] >= exp_lower]
        data = data.sort_values(by='100')
        describe = data['100'].describe()
        regions = [
            (describe['min'], describe['25%']),
            (describe['25%'], describe['50%']),
            (describe['50%'], describe['75%']),
            (describe['75%'], describe['max']),
        ]
        errors = data.apply(lambda column: column / data.loc[:, '100'], axis=0)
        errors = ((errors - 1)*100).abs()
        plot_data = list()
        for lower, upper in regions:
            tmp = errors[(data['100'] >= lower) & (data['100'] <= upper)]
            plot_data.append(tmp)

        fig = tools.make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Q1: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[0]),
                'Q2: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[1]),
                'Q3: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[2]),
                'Q4: {x[0]:.2f}-{x[1]:.2f}'.format(x=regions[3]),
            ),
            shared_xaxes=False
        )
        for ind, each in enumerate(plot_data):
            x = 1 if (ind+1) / 2 <= 1 else 2
            y = 1 if (ind+1) % 2 == 1 else 2
            for col in each.columns:
                upper_limit = each[col].describe()['50%']*outlier_limit
                y_data = each[col][each[col] <= upper_limit]
                # box = go.Box(y=y_data, name=col, showlegend=False, boxpoints=False)
                box = go.Box(y=y_data, name=col, showlegend=False)
                fig.append_trace(box, x, y)
            line = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines', name='median')
            fig.append_trace(line, x, y)
            if ind != 0:
                line2 = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines',
                                   legendgroup=sample, name=sample, line=dict(color=sample_color), showlegend=False)
            else:
                line2 = go.Scatter(x=each.columns, y=each.describe().loc['50%', :], mode='lines',
                                   legendgroup=sample, name=sample, line=dict(color=sample_color))
            all_fig.append_trace(line2, x, y)
        prefix = "{}.TPM.Saturation".format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)

    # plot all sample together
    prefix = "{}.TPM.Saturation".format('samples')
    draw(all_fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def exp_pca(exp_table, row_sum_cutoff=0.1, exp_cutoff=0.1, cv_cutoff=0.01,
            explained_ratio=0.95, outdir=os.getcwd(), group_dict=None,
            log_transform=True, do_scale=True, prefix='',
            formats=('html',), height:int=None, width:int=None, scale=3, desc=None):
    from sklearn import decomposition, preprocessing
    data = pd.read_csv(exp_table, header=0, index_col=0, sep=None, engine='python')
    data = data[data.sum(axis=1) >= row_sum_cutoff]
    pass_state = data.apply(lambda x: sum(x > exp_cutoff), axis=1)
    data = data[pass_state >= int(data.shape[1])/3]
    data = data[data.std(axis=1)/data.mean(axis=1) > cv_cutoff]
    data.to_csv(os.path.join(outdir, '{}pca.filtered.data.txt'.format(prefix)), header=True, index=True, sep='\t')
    if log_transform:
        data = np.log(data+1)
    if do_scale:
        data = data.apply(preprocessing.scale, axis=0)
    data = data.transpose()
    pca = decomposition.PCA()
    pca.fit(data)
    _ratio = list(enumerate(pca.explained_variance_ratio_, start=1))
    total_ratio, n_components = 0, 0
    for ind, each in _ratio:
        total_ratio += each
        if total_ratio >= explained_ratio:
            n_components = ind
            break
    if n_components <= 1:
        n_components = 2
    _ratio = _ratio[:n_components]
    result = pd.DataFrame(pca.transform(data), index=data.index)
    result = result.iloc[:, :n_components]
    result.index.name = 'sample'
    # result.columns = ['PC'+str(n)+'('+'{:.2f}%'.format(r*100)+')' for n, r in _ratio]
    result.columns = ['PC'+str(n) for n in range(1, result.shape[1] + 1)]
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    out_dir = os.path.join(outdir, '{}PCA.xls'.format(prefix))
    result.to_csv(out_dir, sep='\t', header=True, index=True)
    pc_ratio = {'PC' + str(n): r for n, r in _ratio}
    out_dir2 = os.path.join(outdir, 'Explained_variance_ratio.xls')
    with open(out_dir2, 'w') as f:
        for each in sorted(pc_ratio.keys()):
            f.write(str(each) + '\t' + str(pc_ratio[each]) + '\n')

    # color and maker pool
    if type(group_dict) == str:
        if os.path.exists(group_dict):
            with open(group_dict) as f:
                group_dict = OrderedDict(x.strip().split()[:2] for x in f if x.strip())
        else:
            raise FileExistsError('group_dict file not existed!')
    if group_dict is None or not group_dict:
        group_dict = OrderedDict()
    for sample in result.index:
        if sample not in group_dict:
            group_dict[sample] = sample
    group_dict = OrderedDict([(x, y) for x, y in group_dict.items() if x in result.index])
    groups = list()
    for group in group_dict.values():
        if group not in groups:
            groups.append(group)
    colors = get_color_pool(len(group_dict))
    makers = get_marker_pool(len(groups))
    sample_colors = dict(zip(group_dict.keys(), colors))
    group_makers = dict(zip(groups, makers))
    sample_makers = dict()
    for sample in result.index:
        group = group_dict[sample]
        sample_makers[sample] = group_makers[group]

    # plot
    traces = list()
    for sample in result.index:
        trace = go.Scatter(
            x=[result.loc[sample, 'PC1']],
            y=[result.loc[sample, 'PC2']],
            mode='markers',
            name=sample,
            marker=dict(
                color=sample_colors[sample],
                line=dict(color='rgba(217, 217, 217, 0.14)', width=0.5),
                opacity=0.8,
                size=20,
                symbol=sample_makers[sample],
            ),
        )
        traces.append(trace)
    layout = dict(
        showlegend=True,
        # width=800,
        # height=600,
        hovermode='closest',
        xaxis=dict(title='PC1({:.2%})'.format(pc_ratio['PC1'])),
        yaxis=dict(title='PC2({:.2%})'.format(pc_ratio['PC2']))
    )
    fig = go.Figure(traces, layout=layout)
    prefix = "{}PC1_PC2".format(prefix)
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def exp_density(exp_table, outdir=os.getcwd(), row_sum_cutoff=0.1, exp_cutoff=0.1, logbase=2,
                formats=('html',), height:int=None, width:int=None, scale=3, name_dict='', target_samples='', desc=None):
    if name_dict:
        name_dict = dict(x.strip().split()[:2] for x in open(name_dict))
    else:
        name_dict = dict()
    if target_samples:
        target_samples = [x.strip().split()[0] for x in open(target_samples)]

    def get_density(all_exp_pd):
        """
        sampling 800 density point for each columns of the input pandas DataFrame
        data row with log transformation failed will be removed
        :param all_exp_pd: pandas DataFrame
        :return: a list with dict as element
        """
        from scipy import stats
        records = list()
        target_columns = all_exp_pd.columns
        for sample in target_columns:
            exp = all_exp_pd[sample]
            exp = exp.replace([np.inf, -np.inf], np.nan).dropna()
            exp = exp[exp != 0]
            density_func = stats.gaussian_kde(exp)
            min_exp, max_exp = exp.min(), exp.max()
            x_data = np.linspace(min_exp, max_exp, num=800, endpoint=False)
            y_data = density_func(x_data)
            point_df = pd.DataFrame({'exp': x_data, 'density': y_data})
            records.append(point_df)
        return records

    data = pd.read_table(exp_table, header=0, index_col=0)
    data.columns = [name_dict[x] if x in name_dict else x for x in data.columns]
    if target_samples:
        data = data.loc[:, [x for x in target_samples if x in data.columns]]
    if target_samples or name_dict:
        data.to_csv(os.path.join(outdir, 'expression_matrix.xls'), header=True, index=True, sep='\t')
    data = data[data.sum(axis=1) >= row_sum_cutoff]
    pass_state = data.apply(lambda x: sum(x > exp_cutoff), axis=1)
    data = data[pass_state >= int(data.shape[1]) / 3]
    data.to_csv(os.path.join(outdir, 'density.filtered.data.txt'), header=True, index=True, sep='\t')
    if logbase == 2:
        data = np.log2(data)
    elif logbase == 10:
        data = np.log10(data)
    else:
        raise Exception('logbase should be 2 or 10'.format(logbase))
    traces = list()
    density_point_df_list = get_density(data)
    color_pool = get_color_pool(len(data.columns))
    for ind, (sample, color) in enumerate(zip(data.columns, color_pool)):
        # print(sample)
        trace = go.Scatter(
            x=density_point_df_list[ind]['exp'],
            y=density_point_df_list[ind]['density'],
            mode='lines',
            fill='tonexty',
            fillcolor=color,
            name=sample,
            opacity=0.7,
            line=dict(color=color)
        )
        traces.append(trace)

    layout = go.Layout(
        barmode='overlay',
        title='Expression Density',
        xaxis=dict(title='Log{}(Expression)'.format(logbase), zeroline=False),
        yaxis=dict(title='Density', )
    )
    fig = go.Figure(data=traces, layout=layout)
    prefix = "Expression.density"
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def target_region_depth_distribution(files:list, outdir=os.getcwd(), formats=('html',), height:int=None, width:int=None, scale=3, desc=None):
    files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    data_dict = [json.load(open(x), object_pairs_hook=OrderedDict) for x in files]
    sample_list = [os.path.basename(x).split('.', 1)[0] for x in files]
    for distr_dict, sample in zip(data_dict, sample_list):
        layout = go.Layout(
            title="Base depth distribution",
            # xaxis=dict(title='Depth'),
            yaxis=dict(title='Base number ratio', type='log'),
        )
        total = sum(distr_dict.values())
        norm_y = [x/total for x in distr_dict.values()]
        trace = go.Bar(x=list(distr_dict.keys()), y=norm_y, )
        fig = go.Figure(data=[trace], layout=layout)
        prefix = '{}.baseDepthDistribution'.format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def chromosome_read_distribution(files:list, outdir=os.getcwd(), top=30, formats=('html',), height:int=None, width:int=None, scale=3, desc=None):
    all_data = list()
    samples = list()
    files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=None, index_col=0).iloc[:-1, :]
        data = data.sort_values(by=2, ascending=False).loc[:, [2]].iloc[:top, :]
        total = data[2].sum()
        trace = go.Bar(x=data.index, y=data[2])
        layout = go.Layout(
            title='Read distribution on chromosomes/scaffolds',
            # xaxis=dict(title='Chromosomes/Scaffolds'),
            yaxis=dict(title='Mapped read number')
        )
        fig = go.Figure(data=[trace], layout=layout)
        prefix = '{}.ChromosomeReadDistribution'.format(sample)
        draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)
        all_data.append(data[2]/total)
        samples.append(sample)
    # overall
    data = pd.concat(all_data, axis=1, sort=False).fillna(0)
    data.columns = samples
    out_table = os.path.join(outdir, 'ChromosomeReadDistribution.txt')
    data.to_csv(out_table, index=True, header=True, sep='\t')
    data = data.transpose()
    color_pool = get_color_pool(data.shape[1])
    traces = [go.Bar(x=data.index, y=data[x], name=x, marker=dict(color=c)) for x, c in zip(data.columns, color_pool)]
    layout = go.Layout(
        title='Read distribution on chromosomes/scaffolds across samples',
        # xaxis=dict(title='Sample'),
        yaxis=dict(title='Mapped read number percent'),
        barmode='stack'
    )
    fig = go.Figure(data=traces, layout=layout)
    prefix = '{}.ChromosomeReadDistribution'.format('samples')
    draw(fig, prefix=prefix, outdir=outdir, formats=formats, height=height, width=width, scale=scale, desc=desc)


def CollectDupRateMetrics(files:list):
    duplication_rates = list()
    files = sorted(files, key=lambda x: os.path.basename(x).split('.', 1)[0])
    for each in files:
        sample = os.path.basename(each).split('.', 1)[0]
        data = pd.read_table(each, header=0, index_col=None)
        duplication_rate = 1 - data['UniqReadNumber'][0] / data['UniqReadNumber'].sum()
        dup_info = dict()
        dup_info['sample'] = sample
        dup_info['Duplication_rate'] = duplication_rate
        duplication_rates.append(dup_info)
    dup_rate = pd.DataFrame(duplication_rates).set_index('sample')
    return dup_rate


def CollectAlignmentSummaryMetrics(files:list):
    if type(files) == str:
        files = glob(files)
    data = list()
    for each in files:
        if not os.path.exists(each):
            pass
        sample = os.path.basename(each).split('.', 1)[0]
        summary = pd.read_table(each, comment='#', index_col=0, header=0)
        summary =summary.loc[["PAIR"], :]
        summary.index = [sample]
        data.append(summary)
    data = pd.concat(data, axis=0).dropna(axis=1).drop('PCT_ADAPTER', axis=1).round(4)
    # data = data.transpose()
    # out_table = os.path.join(outdir, 'AlignmentSummaryMetrics.xls')
    # data.to_csv(out_table, index=True, header=True, sep='\t')
    return data


def CollectInsertSizeMetrics(files:list):
    if type(files) == str:
        files = glob(files)
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
    # data = data.transpose()
    # out_table = os.path.join(outdir, 'InsertSizeMetrics.xls')
    # data.to_csv(out_table, index=True, header=True, sep='\t')
    return data


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
    data = data.drop(['CORRECT_STRAND_READS', 'INCORRECT_STRAND_READS', 'IGNORED_READS', 'PCT_CORRECT_STRAND_READS'], axis=1)
    # data = data.transpose()
    # out_table = os.path.join(outdir, 'RnaSeqMetrics.xls')
    # data.to_csv(out_table, index=True, header=True, sep='\t')
    return data


def CollectTargetedPcrMetrics(files:list):
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
    # data = data.transpose()
    # out_table = os.path.join(outdir, 'TargetedPcrMetrics.xls')
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


def MergeMetrics(project_outdir, filter_ref, outdir=os.getcwd(), formats=('html',),
                 height:int=None, width:int=None, scale=3, failed_cutoff=1, name_dict=''):
    if name_dict:
        name_dict = dict(x.strip().split()[:2] for x in open(name_dict))
    else:
        name_dict = dict()
    rna_metric = glob(os.path.join(project_outdir, 'CollectRnaSeqMetrics', '*.xls'))
    rna_metric = sorted(rna_metric, key=lambda x: os.path.basename(x).split('.', 1)[0])
    pcr_metric = glob(os.path.join(project_outdir, 'CollectTargetedPcrMetrics', '*.xls'))
    pcr_metric = sorted(pcr_metric, key=lambda x: os.path.basename(x).split('.', 1)[0])
    star_metric = glob(os.path.join(project_outdir, 'Align', '*', '*.Log.final.out'))
    star_metric = sorted(star_metric, key=lambda x: os.path.basename(x).split('.', 1)[0])
    dup_metric = glob(os.path.join(project_outdir, 'ReadDuplication', '*.pos.DupRate.xls'))
    dup_metric = sorted(dup_metric, key=lambda x: os.path.basename(x).split('.', 1)[0])
    metric_dfs = list()
    if rna_metric:
        metric_dfs.append(CollectRnaSeqMetrics(rna_metric))
    if pcr_metric:
        metric_dfs.append(CollectTargetedPcrMetrics(pcr_metric))
    if star_metric:
        metric_dfs.append(CollectStarAlignmentMetrics(star_metric))
    if dup_metric:
        metric_dfs.append(CollectDupRateMetrics(dup_metric))
    if not metric_dfs:
        print('Found no qc file to merge!')
        return
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
    # plot
    pct_cols = [x for x in out_table.index if '_PCT' in x or 'PCT_' in x or '%' in x]
    orders = ['ZERO_CVG_TARGETS_PCT'] + sorted(list(set(pct_cols)-{'ZERO_CVG_TARGETS_PCT'}))
    pct_data = out_table.loc[pct_cols, :].sort_values(by=orders, axis=1)
    pct_data.to_csv(os.path.join(outdir, 'pct_metrics.csv'), encoding='utf_8_sig')
    traces = list()
    colors = get_color_pool(pct_data.shape[0])
    for metric, color in zip(pct_data.index, colors):
        trace = go.Scatter(
            x=pct_data.columns,
            y=pct_data.loc[metric],
            mode='lines',
            name=metric,
            marker=dict(color=color)
        )
        traces.append(trace)
    fig = go.Figure(data=traces)
    draw(fig, prefix="qc_percent_metric", outdir=outdir, formats=formats, height=height, width=width, scale=scale)

    # plot after excluding failed samples
    passed_samples = list()
    for sample, score in sorted(zip(pass_state.columns, pass_state.sum()), key=lambda x:x[1], reverse=True):
        if score > ref.shape[0]-failed_cutoff:
            passed_samples.append(sample)
    if passed_samples:
        pct_cols = [x for x in out_table.index if '_PCT' in x or 'PCT_' in x or '%' in x]
        orders = ['ZERO_CVG_TARGETS_PCT'] + sorted(list(set(pct_cols) - {'ZERO_CVG_TARGETS_PCT'}))
        pct_data = out_table.loc[pct_cols, passed_samples].sort_values(by=orders, axis=1)
        pct_data.to_csv(os.path.join(outdir, 'filtered.pct_metrics.xls'), index=True, header=True, sep='\t')
        traces = list()
        colors = get_color_pool(pct_data.shape[0])
        for metric, color in zip(pct_data.index, colors):
            trace = go.Scatter(
                x=pct_data.columns,
                y=pct_data.loc[metric],
                mode='lines',
                name=metric,
                marker=dict(color=color)
            )
            traces.append(trace)
        fig = go.Figure(data=traces)
        draw(fig, prefix="filtered_qc_percent_metric", outdir=outdir,
             formats=formats, height=height, width=width, scale=scale)

    # save result
    out_table.loc['Pass'] = 'Yes'
    out_table.loc['Pass', pass_state.sum() < ref.shape[0]-failed_cutoff+1] = 'No'
    out_table['describe'] = list(ref['Description']) + ['Pass QC filter']
    out_table.set_index('describe', inplace=True, append=True)
    out_name = os.path.join(outdir, 'QC_Summary.csv')
    out_table.to_csv(out_name, encoding='utf_8_sig')
    # save result2
    target_metrics =[
        'Number of input reads',
        'PCT_mapped_reads',
        'Average input read length',
        'Uniquely mapped reads number',
        'Uniquely mapped reads %',
        'Number of reads mapped to multiple loci',
        '% of reads mapped to multiple loci',
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
    out_name = os.path.join(outdir, 'QC_ReportSummary.csv')
    out_table.loc[target_metrics, :].to_csv(out_name, encoding='utf_8_sig')


def diff_volcano(files: list, outdir='', formats=('html', ), gene_annot=None,
                 fc_cutoff=2.0, sig_type='corrected_pvalue', pvalue_cutoff=0.05, # ????????????,???????????????
                 limit=8, height:int=None, width:int=None, scale=3):
    """
    :param files: ??????['log2fc', 'pvalue/padjust', 'significant'(='yes' or 'no'), 'regulate'(='up' or 'down')]?????????????????????
    :param outdir:
    :param formats:
    :param limit: ????????????outlier q3+(q3-q1)*limit
    :param height:
    :param width:
    :param scale:
    :return:
    """
    gene_annot = dict(x.strip().split('\t')[:2] for x in open(gene_annot)) if gene_annot else dict()
    table_desc = f"#display_columns: gene_symbol;log2fc;pvalue;padjust;significant;regulate;genecards\n#search_columns: 0;1\n?????????????????????????????????seq_id:Ensembl gene id; log2fc?????????log2???????????????????????? pvalue??????????????????????????????????????????	padjust?????????BH?????????pvalue, ???????????????untested???????????????????????????count???????????????????????????????????????significant?????????????????????????????????|fc|>={fc_cutoff}???{sig_type}<={pvalue_cutoff}???regulate???up???????????????down???????????????gene_symbol?????????????????????gene_type:????????????????????????https://www.gencodegenes.org/pages/biotypes.html???genecards???gene???genecards??????????????????*_count: raw count???*_normalized?????????????????????????????????\n"
    volcano_desc = "??????????????????????????????????????????:??????????????????????????????log2????????????????????????????????????BH??????????????????pvalue??????-log10??????????????????????????????????????????/???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????"
    p_type = 'padjust' if sig_type != 'pvalue' else 'pvalue'
    def pvalStr2Float(x):
        try:
            return float(x)
        except:
            return 1
    table_dir = os.path.join(outdir, '1.DiffTable')
    volcano_dir = os.path.join(outdir, '2.Volcano')
    os.mkdir(table_dir)
    os.mkdir(volcano_dir)
    link_files = []
    for table in files:
        link = os.path.join(table_dir, os.path.basename(table))
        link_files.append(os.path.abspath(link))
        os.symlink(os.path.abspath(table), link)

    for table in link_files:
        with open(table + '.describe.txt', 'w') as f:
            f.write(table_desc)
        ctrl, test = re.fullmatch(r'(.*)_vs_(.*?)\..*.xls', os.path.basename(table)).groups()
        df = pd.read_table(table, index_col=0, header=0)
        df[p_type] = df[p_type].apply(pvalStr2Float)
        exp_columns = [x for x in df.columns if x.lower().endswith(('tpm', 'fpkm', 'normalized'))]
        target_columns = exp_columns + ['log2fc', p_type, 'significant', 'regulate']
        df = df.loc[:, target_columns]
        df = df[df['regulate'] != 'untested']
        if exp_columns:
            df = df[df.loc[:, exp_columns].sum(axis=1) > 1e-5]
        fc_desc = df['log2fc'].describe()
        fc_upper_limit = fc_desc['75%'] + (fc_desc['75%'] - fc_desc['25%'])*limit
        fc_lower_limit = fc_desc['25%'] - (fc_desc['75%'] - fc_desc['25%'])*limit
        df.loc[(df['log2fc'] >= fc_upper_limit), ['log2fc']] = fc_upper_limit
        df.loc[(df['log2fc'] <= fc_lower_limit), ['log2fc']] = fc_lower_limit
        df.loc[df[p_type] == 0, p_type] = df[df[p_type] > 0][p_type].min()
        df[p_type] = -np.log10(df[p_type])
        p_desc = df[p_type].describe()
        p_limit = p_desc['75%'] + (p_desc['75%'] - p_desc['25%'])*limit
        df.loc[(df[p_type] >= p_limit), p_type] = p_limit

        df['colors'] = 'grey'
        df.loc[((df['significant'] == 'yes') & (df['regulate'] == 'up')), 'colors'] = 'red'
        df.loc[((df['significant'] == 'yes') & (df['regulate'] == 'down')), 'colors'] = 'green'

        # draw grey
        tmp_df = df[df['colors'] == 'grey']
        trace = go.Scatter(
            x=tmp_df['log2fc'],
            y=tmp_df[p_type],
            mode='markers',
            marker=dict(
                color='grey',
                opacity=0.8,
                # symbol='circle',
            ),
            name='Not Significant: {}'.format(tmp_df.shape[0])
        )
        # draw green
        tmp_df = df[df['colors'] == 'green']
        trace2 = go.Scatter(
            x=tmp_df['log2fc'],
            y=tmp_df[p_type],
            mode='markers',
            text=[gene_annot[x] if x in gene_annot else x for x in tmp_df.index],
            hoverinfo='text+x+y',
            marker=dict(
                color='green',
                opacity=0.8,
                # symbol='circle',
            ),
            name = 'Down-regulated: {}'.format(tmp_df.shape[0])
        )
        # draw red
        tmp_df = df[df['colors'] == 'red']
        trace3 = go.Scatter(
            x=tmp_df['log2fc'],
            y=tmp_df[p_type],
            mode='markers',
            text=[gene_annot[x] if x in gene_annot else x for x in tmp_df.index],
            hoverinfo='text+x+y',
            marker=dict(
                color='red',
                opacity=0.8,
                # symbol='circle',
            ),
            name='Up-regulated: {}'.format(tmp_df.shape[0])
        )

        layout = dict(
            title='{} vs {}'.format(ctrl, test),
            autosize=True,
            # margin=dict(t=25, l=10, r=10, b=10),
            showlegend=True,
            xaxis=dict(title='log2FC'),
            yaxis=dict(title=f'-log10({sig_type})'),
            hovermode='closest',
        )

        fig = go.Figure(data=[trace, trace2, trace3], layout=layout)
        prefix = '{}_vs_{}.volcano'.format(ctrl, test)
        draw(fig, prefix=prefix, outdir=volcano_dir, formats=formats, height=height, width=width, scale=scale, desc=volcano_desc)


def go_enriched_term_bubble(files: list, top=20, outdir='', formats=('html', ), gene_annot=None, use_pval=False,
                            color_scale='Rainbow', limit=5, height:int=None, width:int=None, scale=3):
    gene_annot = dict(x.strip().split('\t')[:2] for x in open(gene_annot)) if gene_annot else dict()
    p_type = 'p_uncorrected' if use_pval else 'p_corrected'
    descs = [
        "????????????????????????????????????BP(biological process) ????????????20???????????????(??????pvalue??????)?????????????????????????????????????????????????????????????????????????????????-log10???????????? p-value?????????BH????????????????????????????????????????????????????????????????????????????????????ratio_in_study/ratio_in_pop, ???????????????GO????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????GO???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????GO????????????????????????\n",
        "????????????????????????????????????CC(celluar component) ????????????20???????????????(??????pvalue??????)?????????????????????????????????????????????????????????????????????????????????-log10???????????? p-value?????????BH????????????????????????????????????????????????????????????????????????????????????ratio_in_study/ratio_in_pop, ???????????????GO????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????GO???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????GO????????????????????????\n",
        "????????????????????????????????????MF(molecular function) ????????????20???????????????(??????pvalue??????)?????????????????????????????????????????????????????????????????????????????????-log10???????????? p-value?????????BH????????????????????????????????????????????????????????????????????????????????????ratio_in_study/ratio_in_pop, ???????????????GO????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????GO???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????GO????????????????????????\n"
    ]
    table_desc = "#display_columns: NS;name;ratio_in_study;ratio_in_pop;p_corrected\nGene Ontology?????????????????????????????????GO???GO ID???NS??????????????????????????????BP(biological process)???CC(celluar component)???MF(molecular function)?????????enrichment: e????????????(enriched)???p???????????????(purified)???name???????????????ratio_in_study??????????????????GO????????????????????????????????????????????????????????????????????????????????????ratio_in_pop??????????????????GO???????????????????????????????????????????????????????????????????????????????????????????????????????????????p_uncoorected: fisher???????????????pvalue; depth: ?????????????????????GO??????????????????????????????????????????????????????????????????????????????????????????????????????study_count: ???????????????????????????????????????p_corrected?????????BH???????????????pvalue???study_items?????????????????????????????????????????????down???up??????????????????????????????\n"
    # *.dag.BP.svg
    dag_desc = "??????{}????????????20???????????????(??????pvalue??????)?????????????????????????????????????????????DAG??????????????????????????????????????????????????????Lxx?????????shortest distance from root node??????Dxx?????????longest distance from root node?????????????????????????????????(level???????????????????????????????????????????????????pvalue???????????????????????????????????????legend??????, ???????????????????????????pvalue>0.05???\n"
    bp_dag = glob(os.path.join(os.path.dirname(files[0]), '*.dag.BP.svg'))
    for dag in bp_dag:
        with open(dag+'.describe.txt', 'w') as f:
            f.write(dag_desc.format('BP???biological process???'))
    cc_dag = glob(os.path.join(os.path.dirname(files[0]), '*.dag.CC.svg'))
    for dag in cc_dag:
        with open(dag + '.describe.txt', 'w') as f:
            f.write(dag_desc.format('CC???celluar component???'))
    mf_dag = glob(os.path.join(os.path.dirname(files[0]), '*.dag.MF.svg'))
    for dag in mf_dag:
        with open(dag + '.describe.txt', 'w') as f:
            f.write(dag_desc.format('MF???molecular function???'))

    for table in files:
        print(table)
        with open(table+'.describe.txt', 'w') as f:
            f.write(table_desc)
        df_big = pd.read_table(table, index_col=0, header=0)
        for each, desc in zip(['BP', 'CC', 'MF'], descs):
            prefix = os.path.basename(table)[:-4]
            df = df_big[df_big['NS']==each].sort_values(by=[p_type]).head(top)
            gene_number = [int(x.split('/')[0]) for x in df['ratio_in_study']]
            pval_color = -np.log10(df[p_type])
            pval_color_desc = pval_color.describe()
            pval_color[pval_color > pval_color_desc['75%']*limit] = pval_color_desc['75%']*limit
            score = [round(eval(x)/eval(y), 2) for x, y in zip(df['ratio_in_study'], df['ratio_in_pop'])]
            text_list = []
            for x, y, h, g, n, p in zip(gene_number, df.index, df['NS'], df['study_items'], df['name'], df[p_type]):
                text = '{}: {} <br>'.format(y, n)
                text += 'type: {} <br>'.format(h)
                text += '{}: {} <br>'.format(p_type, format(p, '.3e') if p < 0.001 else round(p, 4))
                genes = [x.split('|')[0] for x in g.split(';')][:100]
                genes = [gene_annot[x] if x in gene_annot else x for x in genes]
                gene_regulate = [x.split('|')[1] for x in g.split(';')]
                if gene_regulate.count('up') or gene_regulate.count('down'):
                    text += 'gene regulate: {} up while {} down <br>'.format(gene_regulate.count('up'),
                                                                             gene_regulate.count('down'))
                else:
                    text += 'gene_number: {} <br>'.format(len(genes))
                text += 'genes: {}'.format('<br>'.join(textwrap.wrap(';'.join(genes))))
                text_list.append(text)
            x_data = score
            y_data = ['<br>'.join(textwrap.wrap(x, width=80)) for x in df['name']]
            trace = go.Scatter(
                y=y_data,
                x=score,
                mode='markers',
                text=text_list,
                hoverinfo='text+x',
                marker=dict(
                    size=gene_number,
                    sizemode='area',
                    sizeref=2. * max(gene_number) / (30. ** 2),
                    sizemin=4,
                    color=pval_color,
                    showscale=True,
                    colorbar=dict(
                        title='-log10(P-value)',
                    ),
                    colorscale=color_scale
                )
            )
            layout = dict(
                autosize=True,
                margin=dict(l=350),
                title='Top {} enriched {} GO term'.format(top, each),
                xaxis=dict(title='Enrichment Ratio (ratio_in_study/ratio_in_pop)'),
                yaxis=dict(dtick=1, tickfont={'size': 8}),
                hovermode='closest',
            )
            links = []
            for x, y, link in zip(x_data, y_data, df.index):
                links.append(
                    dict(x=x, y=y,
                         text="""<a href="http://amigo.geneontology.org/amigo/term/{}">{}</a>""".format(link, "   "),
                         showarrow=False, xanchor='center', yanchor='middle', )
                )
            layout['annotations'] = links
            fig = go.Figure(data=[trace], layout=layout)
            prefix = '{}.{}.bubble'.format(prefix, each)
            draw(fig, prefix=prefix, outdir=os.path.dirname(table), formats=formats, height=height, width=width, scale=scale, desc=desc)

    # reorganise
    result = os.path.dirname(os.path.abspath(files[0]))
    table_dir = os.path.join(outdir, '1.EnrichTable')
    bp_dir = os.path.join(outdir, '2.BP')
    mf_dir = os.path.join(outdir, '3.MF')
    cc_dir = os.path.join(outdir, '4.CC')
    os.mkdir(table_dir)
    os.mkdir(bp_dir)
    os.mkdir(mf_dir)
    os.mkdir(cc_dir)
    for each in glob(os.path.join(result, '*.goea.xls*')):
        os.symlink(each, os.path.join(outdir, table_dir, os.path.basename(each)))
    for each in glob(os.path.join(result, '*.BP.*')):
        os.symlink(each, os.path.join(outdir, bp_dir, os.path.basename(each)))
    for each in glob(os.path.join(result, '*.MF.*')):
        os.symlink(each, os.path.join(outdir, mf_dir, os.path.basename(each)))
    for each in glob(os.path.join(result, '*.CC.*')):
        os.symlink(each, os.path.join(outdir, cc_dir, os.path.basename(each)))


def kegg_enriched_term_bubble(files: list, top=20, outdir='', formats=('html', ), fdr_cutoff=0.05, use_pval=False,
                              gene_annot=None, color_scale='Rainbow', limit=5, height:int=None, width:int=None, scale=3):
    gene_annot = dict(x.strip().split('\t')[:2] for x in open(gene_annot)) if gene_annot else dict()
    desc = "???????????????????????????????????????20????????????pvalue??????)?????????????????????????????????????????????????????????????????????????????????????????????-log10???????????? p-value?????????BH????????????????????????????????????????????????????????????????????????????????????ratio_in_study/ratio_in_pop, ???????????????pathway??????????????????????????????????????????????????????pathway????????????????????????????????????????????????????????????pathway???????????????????????????????????????????????????pathway????????????????????????????????????????????????????????????pathway???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????kegg pathway????????????????????????????????????kegg????????????????????????pathway??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n"
    p_type = 'P-Value' if use_pval else 'Corrected P-Value'
    table_desc = '#display_columns: ID;Ratio_in_study;Ratio_in_pop;Corrected P-Value;typeII;Hyperlink\nKEGG pathway????????????????????? Term???pathway????????????DataBase???????????????????????????ID???Pathway ID???Ratio_in_study??????????????????pathway???????????????????????????????????????????????????????????????pathway???????????????Ratio_in_pop??????????????????pathway??????????????????????????????????????????????????????????????????????????????????????????pathway???????????????P-Value: fisher???????????????pvalue; Corrected Pvalue?????????BH???????????????pvalue???Genes?????????????????????pathway????????????????????????down???up??????????????????????????????Hyperlink???KEGG pathway?????????????????????kegg????????????????????????pathway??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????typeII?????????????????????????????????typeI????????????????????????????????????????????????<a href="https://www.genome.jp/kegg/pathway.html">KEGG</a>.\n'
    for table in files:
        with open(table+'.describe.txt', 'w') as f:
            f.write(table_desc)
        prefix = os.path.basename(table)[:-4]
        # df0 = pd.read_table(table, index_col=0, header=0, nrows=top)
        df = pd.read_table(table, index_col=0, header=0)
        df0 = df.sort_values(by=[p_type]).iloc[:top, :]
        df = df0.loc[df0[p_type]<=fdr_cutoff, :]
        if df.shape[0] >= 1:
            print('{} out of top {} are significant'.format(df.shape[0], top))
        else:
            print("No significant term found in {}, still we plot the top {} for you".format(table, top))
            df = df0
        gene_number = [int(x.split('/')[0]) for x in df['Ratio_in_study']]
        pval_color = -np.log10(df[p_type])
        pval_color[pval_color == np.inf] = 10
        pval_color_desc = pval_color.describe()
        pval_color[pval_color > pval_color_desc['75%']*limit] = pval_color_desc['75%']*limit
        score = [round(eval(x)/eval(y), 2) for x, y in zip(df['Ratio_in_study'], df['Ratio_in_pop'])]
        text_list = []
        for x, y, h, k, g, p in zip(gene_number, df.ID, df['typeI'], df['typeII'], df['Genes'], df[p_type]):
            text = ''
            text += 'path: {} <br>'.format(y)
            text += '{}: {} <br>'.format(p_type, format(p, '.3e') if p < 0.001 else round(p,4))
            text += 'typeI: {} <br>'.format(h)
            text += 'typeII: {} <br>'.format(k)
            genes = [x.split('|')[0] for x in g.split(';')][:100]
            genes = [gene_annot[x] if x in gene_annot else x for x in genes]
            gene_regulate = [x.split('|')[2].lower() for x in g.split(';')]
            if gene_regulate.count('up') or gene_regulate.count('down'):
                text += 'gene regulate: {} up while {} down <br>'.format(gene_regulate.count('up'), gene_regulate.count('down'))
            else:
                text += "gene number: {} <br>".format(len(genes))
            text += 'genes: {}'.format('<br>'.join(textwrap.wrap(';'.join(genes))))
            text_list.append(text)
        y_data = ['<br>'.join(textwrap.wrap(x, width=80)) for x in df.index]
        x_data = score
        trace = go.Scatter(
            y=y_data,
            x=x_data,
            mode='markers',
            text=text_list,
            hoverinfo='text+x+y',
            marker=dict(
                size=gene_number,
                sizemode='area',
                sizeref=2. * max(gene_number) / (30. ** 2),
                sizemin=4,
                color=pval_color,
                showscale=True,
                colorbar=dict(
                    title='-log10(P-value)',
                ),
                colorscale=color_scale
            )
        )
        layout = dict(
            autosize=True,
            title='Top {} enriched KEGG pathway'.format(top),
            margin=dict(l=350),
            xaxis=dict(title='Enrichment Ratio (ratio_in_study/ratio_in_pop)'),
            yaxis=dict(dtick=1, tickfont={'size': 9}),
            hovermode='closest',
        )
        links = []
        for x, y, link in zip(x_data, y_data, df['Hyperlink']):
            links.append(dict(x=x, y=y,
                              text="""<a href="{}">{}</a>""".format(link, "   "),
                              showarrow=False,
                              xanchor='center', yanchor='middle',))
        layout['annotations'] = links
        fig = go.Figure(data=[trace], layout=layout)
        prefix = '{}.bubble'.format(prefix)
        draw(fig, prefix=prefix, outdir=os.path.dirname(table),
             formats=formats, height=height, width=width, scale=scale, desc=desc)

    # reorganise
    result = os.path.dirname(os.path.abspath(files[0]))
    table_dir = os.path.join(outdir, '1.EnrichTable')
    plot_dir = os.path.join(outdir, '2.EnrichPlot')
    os.mkdir(table_dir)
    os.mkdir(plot_dir)
    for each in glob(os.path.join(result, '*.xls*')):
        os.symlink(each, os.path.join(outdir, table_dir, os.path.basename(each)))
    for each in glob(os.path.join(result, '*.bubble.*')):
        os.symlink(each, os.path.join(outdir, plot_dir, os.path.basename(each)))


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), log=False)
