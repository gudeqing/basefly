import os
import json
import pandas as pd


def merge_umi_qc(
        fastp_json_files: list,
        bamdst_cov_report_files: list,
        bamdst_depth_files: list,
        bamdst_insertsize_files: list,
        umi_family_size_files: list,
        outdir='.',
        prefix=None,
):
    """
    用于umi数据分析流程的质控汇总：fastp的分析结果，bamdst的深度统计，umi统计, 最好仅用于汇总一个样本的分析结果
    :param fastp_json_files: [json, json2], 根据文件名称信息提取样本名信息， sample_name = file_name.split('.')[0]
    :param : bamdst_cov_report_files:  第一个UMI处理前的，第二个是UMI处理后的; 需要和bamdst_depth_files一一对应
    :param : bamdst_depth_files: 第一个UMI处理前的，第二个是UMI处理后的; 需要和bamdst_cov_report_files一一对应
    :param : bamdst_insertsize_files: 第一个UMI处理前的，第二个是UMI处理后的; 需要和bamdst_cov_report_files一一对应
    :param umi_family_size_files: [family_size, family_size2, ...], 根据文件名称信息提取样本名信息， sample_name = file_name.split('.')[0]
    :param outdir: 输出结果目录
    :param prefix: 输出文件名前缀, 默认用样本名称作为前缀
    :return:
    """
    print('Merging QC table......')
    result = dict()

    # 汇总fastp信息
    samples = []
    for stat_file in fastp_json_files:
        if not os.path.exists(stat_file):
            print(f'skip un-existing file {stat_file}')
            continue
        # 文件名用‘.‘分割后，取第一个值作为样本名称
        sample = os.path.basename(stat_file).split('.', 1)[0]
        samples.append(sample)
        json_dict = json.load(open(stat_file))
        target_info = dict()
        # target_info['sequencing'] = json_dict['summary']['sequencing']
        target_info['Read1_length'] = json_dict['summary']['before_filtering']['read1_mean_length']
        target_info['Read2_length'] = json_dict['summary']['before_filtering']['read2_mean_length']
        target_info['Total_raw_reads'] = json_dict['summary']['before_filtering']['total_reads']
        target_info['Total_raw_bases'] = json_dict['summary']['before_filtering']['total_bases']
        target_info['Q20_rate'] = json_dict['summary']['before_filtering']['q20_rate']
        target_info['Q30_rate'] = json_dict['summary']['before_filtering']['q30_rate']
        target_info['Q20_bases'] = json_dict['summary']['before_filtering']['q20_bases']
        target_info['Q30_bases'] = json_dict['summary']['before_filtering']['q30_bases']
        target_info['GC_content'] = json_dict['summary']['before_filtering']['gc_content']
        target_info['Duplication(sequence-based)'] = json_dict['duplication']['rate']
        target_info['Insert_size_peak'] = json_dict['insert_size']['peak']
        result[sample] = target_info

    # 汇总bamdst的信息
    for idx, (stat_file, depth_file, size_file) in enumerate(zip(bamdst_cov_report_files, bamdst_depth_files, bamdst_insertsize_files)):
        sample = samples[0] if len(samples) == 1 else samples[idx]
        if not os.path.exists(stat_file):
            print(f'skip un-existing file {stat_file}')
            continue
        target_info = dict()
        with open(stat_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                name, value = line.strip().split('\t')
                if '%' in value:
                    value = round(float(value.replace('%', ''))*0.01, 4)
                else:
                    value = float(value)
                target_info[name] = value
        # 计算均一性
        data = pd.read_table(depth_file, header=0, low_memory=False)
        # And the "rmdup depth" is calculated after remove duplicated reads and secondary alignment
        # reads and low map quality reads(mapQ < 20), this value is similar with the output depth of samtools depth.
        mean_depth = data['Rmdup depth'].mean()
        panel_size = data.shape[0]
        target_info['[Target] Coverage (>=0.2*MeanDepth)'] = sum(data['Rmdup depth'] >= mean_depth*0.2)/panel_size
        target_info['[Target] Coverage (>=0.5*MeanDepth)'] = sum(data['Rmdup depth'] >= mean_depth*0.5)/panel_size
        target_info['[Target] Coverage (>=200x)'] = sum(data['Rmdup depth'] >= 200)/panel_size
        target_info['[Target] Coverage (>=300x)'] = sum(data['Rmdup depth'] >= 300)/panel_size
        target_info['[Target] Coverage (>=500x)'] = sum(data['Rmdup depth'] >= 500)/panel_size
        target_info['[Target] Coverage (>=1000x)'] = sum(data['Rmdup depth'] >= 1000)/panel_size
        target_info['[Target] Coverage (>=2000x)'] = sum(data['Rmdup depth'] >= 2000)/panel_size
        target_info['[Target] Coverage (>=5000x)'] = sum(data['Rmdup depth'] >= 5000)/panel_size
        target_info['[Target] Coverage (>=10000x)'] = sum(data['Rmdup depth'] >= 10000)/panel_size
        target_info['[Target] Fold80BasePenalty'] = mean_depth/data['Rmdup depth'].quantile(0.2)
        # 计算insert size
        data = pd.read_table(size_file, header=None)
        peak_insert_size = data[data[1] == data[1].max()][0].mean()
        mean_insert_size = (data[0] * data[1]).sum() / data[1].sum()
        target_info['[Total] Peak insert size'] = peak_insert_size
        target_info['[Total] Mean insert size'] = mean_insert_size
        # 假设第一组数据是UMI处理前的统计结果，第二组数据是UMI consensus后的结果，加上标签以示区别
        if idx == 0:
            target_info = {'[Raw]' + key: value for key, value in target_info.items()}
        else:
            target_info = {'[Consensus]' + key: value for key, value in target_info.items()}
        # 样本名称沿用前面fastp文件名中提取的，因为bamdst文件名是固定的，不包含样本名信息
        if sample in result:
            result[sample].update(target_info)
        else:
            result[sample] = target_info

    # 汇总UMI信息
    for stat_file in umi_family_size_files:
        if os.path.exists(stat_file):
            # 文件名用‘.‘分割后，取第一个值作为样本名称
            sample = os.path.basename(stat_file).split('.', 1)[0]
            with open(stat_file) as f:
                target_info = dict()
                header = f.readline()
                for line in f:
                    size, count, fraction, ratio = line.strip().split()
                    target_info[f'UmiFamilySize={size}:count'] = int(count)
                    target_info[f'UmiFamilySize={size}:fraction'] = round(float(fraction), 4)
                    target_info[f'UmiFamilySize>={size}:fraction'] = round(float(ratio), 4)
                if sample in result:
                    result[sample].update(target_info)
                else:
                    result[sample] = target_info

    # 保存结果
    if not result:
        print('Nothing merged!')
        return
    table = pd.DataFrame(result)
    table = table.loc[:, sorted(result.keys())]
    table.index.name = 'Metrics'
    table.loc["Total_raw_reads(M)", :] = table.loc["Total_raw_reads"] * 1e-6
    table.loc["Total_raw_bases(G)", :] = table.loc["Total_raw_bases"] * 1e-9
    prefix = prefix or sample
    table.to_csv(os.path.join(outdir, f'{prefix}.QC.all.metrics.txt'), sep='\t')
    table.to_excel(os.path.join(outdir, f'{prefix}.QC.all.metrics.xlsx'))

    # 精简QC报告
    target_metrics = {
        "Total_raw_bases(G)": "TotalRawBasesGb",
        "Total_raw_reads(M)": "TotalRawReadsMb",
        "GC_content": "GC",  # report
        "Duplication(sequence-based)": "SeqDupRate",  # report
        # "Insert_size_peak": "InsertSize",  # report
        "[Consensus][Total] Peak insert size": "PeakInsertSize",
        "[Consensus][Total] Mean insert size": "MeanInsertSize",
        "Q20_rate": "Q20",
        "Q30_rate": "Q30",  # report
        "[Raw][Total] Fraction of Mapped Reads": "MappingRate",  # report
        "[Raw][Target] Average depth": "RawMeanDepth",
        "[Raw][Total] Fraction of MapQ reads in mapped reads": "MapQ20Rate",  # report
        "[Raw][Target] Fraction of Target Reads in all reads": "OnTargetReadRate",  # report
        "[Raw][Target] Fraction of Target Data in all data": "OnTargetBaseRate",  # report
        "[Raw][Target] Coverage (>=0.2*MeanDepth)": "Uniformity20",
        "[Raw][Target] Coverage (>=0.5*MeanDepth)": "Uniformity50",
        "[Raw][Target] Fold80BasePenalty": "Fold80BasePenalty",  # report
        "[Raw][flank] Fraction of flank Reads in all reads": "OnFlankRate",
        "[Raw][Target] Coverage (>=200x)": "RawDepthOver200x",
        "[Raw][Target] Coverage (>=500x)": "RawDepthOver500x",
        "[Raw][Target] Coverage (>=2000x)": "RawDepthOver2000x",
        "[Raw][Target] Coverage (>=5000x)": "RawDepthOver5000x",
        "[Consensus][Target] Average depth": "UmiMeanDepth",  # report
        "[Consensus][Total] Fraction of MapQ reads in mapped reads": "UmiMapQ20Rate",
        "[Consensus][Target] Fraction of Target Reads in all reads": "UmiOnTargetReadRate",
        "[Consensus][Target] Fraction of Target Data in all data": "UmiOnTargetBaseRate",
        "[Consensus][Target] Coverage (>=0.2*MeanDepth)": "UmiUniformity20",
        "[Consensus][Target] Coverage (>=0.5*MeanDepth)": "UmiUniformity50",  # report
        "[Consensus][flank] Fraction of flank Reads in all reads": "UmiOnFlankRate",
        "[Consensus][Target] Coverage (>=200x)": "UmiDepthOver200x",
        "[Consensus][Target] Coverage (>=500x)": "UmiDepthOver500x",
        "[Consensus][Target] Coverage (>=1000x)": "UmiDepthOver1000x",  # report
        "[Consensus][Target] Coverage (>=2000x)": "UmiDepthOver2000x",
        "[Consensus][Target] Coverage (>=5000x)": "UmiDepthOver5000x",
        "UmiFamilySize=1:count": "UmiFamilySize1Count",
        "UmiFamilySize=2:count": "UmiFamilySize2Count",
        "UmiFamilySize=3:count": "UmiFamilySize3Count",
        "UmiFamilySize=1:fraction": "UmiFamilySize1",
        "UmiFamilySize=2:fraction": "UmiFamilySize2",
        "UmiFamilySize=3:fraction": "UmiFamilySize3",
        "UmiFamilySize>=2:fraction": "UmiFamilySizeOver2",  # report
        "UmiFamilySize>=3:fraction": "UmiFamilySizeOver3",  # report
    }
    target_rows = [x for x in target_metrics.keys()]
    table = table.loc[target_rows]
    table.index = [x for x in target_metrics.values()]
    table.to_excel(os.path.join(outdir, f'{prefix}.QC.target.metrics.xlsx'))
    table = table.T.round(4)
    table.index.name = 'Sample'
    table.to_csv(os.path.join(outdir, f'{prefix}.QC.target.metrics.txt'), sep='\t')
    table.to_csv(os.path.join(outdir, f'{prefix}.QC.target.metrics.csv'))


if __name__ == '__main__':
    # from xcmds import xcmds
    # xcmds.xcmds(locals())
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fastp_json_files', required=False, nargs='+', help='fastp output json file')
    parser.add_argument('-bamdst_cov_report_files', required=False, nargs='+', help='bamdist output coverage report file')
    parser.add_argument('-bamdst_depth_files', required=False, nargs='+', help='bamdst output coverage report file')
    parser.add_argument('-bamdst_insertsize_files', required=False, nargs='+', help='bamdst output coverage report file')
    parser.add_argument('-umi_family_size_files', required=False, nargs='+', help='UMI family size stat file from tool GroupReadsByUmi')
    args = parser.parse_args()
    merge_umi_qc(
        fastp_json_files=args.fastp_json_files,
        bamdst_cov_report_files=args.bamdst_cov_report_files,
        bamdst_depth_files=args.bamdst_depth_files,
        bamdst_insertsize_files=args.bamdst_insertsize_files,
        umi_family_size_files=args.umi_family_size_files
    )
