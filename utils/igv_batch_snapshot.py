import pandas as pd


def generate_igv_snapshot_script(table_file, out_prefix, bam_file_pattern=None, img_dir=None, add_chr='chr', min_af=0.0005):
    if bam_file_pattern is None:
        bam_file_pattern = f'D:\\haixi\\2023-EQA-LungCaner\\{out_prefix}\\' + 'ABRA2-{sample}\\{sample}.realigned.bam'
    if img_dir is None:
        img_dir = f'D:\\haixi\\2023-EQA-LungCaner\\{out_prefix}\\IGV'
    table = pd.read_excel(table_file)
    contents = []
    sample_loaded = set()
    for idx, row in table.iterrows():
        sample = row['Sample']
        if row['Report'] != 'yes':
            continue
        if row['VAF(%)'] <= min_af:
            continue
        if row['SelectedByCsq'] != 'Yes':
            continue
        if row['PassBGN'] == False:
            continue

        if sample not in sample_loaded:
            sample_loaded.add(sample)
            contents.append('\n')
            contents.append('new')
            contents.append('genome hg19')
            contents.append('load {}'.format(bam_file_pattern.format(sample=sample)))
            contents.append('snapshotDirectory {}'.format(img_dir))
        if row['HGVS_OFFSET']:
            row['Start'] = row['Start'] - row['HGVS_OFFSET']
            if row['End'] != '-':
                row['End'] = row['End'] - row['HGVS_OFFSET']
        # if row['End'] != '-':
        #     contents.append(f"goto {add_chr}{row['Chr']}:{row['Start']}-{row['End']}")
        # else:
        contents.append(f"goto {add_chr}{row['Chr']}:{row['Start']}")
        contents.append('sort base')
        contents.append('expand')
        # contents.append('collapse')
        contents.append('colorBy READ_STRAND')
        contents.append(f'snapshot {sample}.{row["Type"]}.{add_chr}{row["Chr"]}-{row["Start"]}-{row["Gene"]}-{row["pHGVS"]}.png')
        contents.append('\n')

    with open(f'{out_prefix}.igv.snapshot.scripts.txt', 'w') as f:
        for each in contents:
            f.write(each+'\n')

    target_info = []
    for idx, row in table.iterrows():
        if row['Report'] != 'yes':
            continue

        region = None
        if row['pHGVS'] == row['pHGVS']:
            region = 'exonic'
        elif "UTR_variant" in row['Consequence']:
            region = 'UTR'
        elif "splice_" in row['Consequence']:
            region = 'splicing'
        else:
            region = 'intron'

        mut_type = None
        if row["Type"] == 'SNV':
            mut_type = 'Substitution'
        elif row['Type'] == 'Complex':
            if len(row["Ref"]) == len(row["Alt"]):
                mut_type = 'Substitution'
            else:
                mut_type = 'Indel'
        elif row['cHGVS'] and row['cHGVS'].endswith('dup'):
            mut_type = "Duplication"
        else:
            mut_type = row['Type'].capitalize()

        target_info.append({
            "备注": row["备注"],
            "Tissue": row["Tissue"] if 'blood' in out_prefix.lower() else row['Blood'],
            "Report": row["Report"],
            "Sample": row["Sample"],
            "Gene": row['Gene'] + '('+row['Transcript'] + ')',
            "Chr": row["Chr"],
            "Start": row["Start"],
            "Ref": row["Ref"],
            "Alt": row["Alt"],
            "cHGVS": row["cHGVS"],
            "pHGVS": row["pHGVS"],
            "Region": region,
            "Exon": row["Exon"],
            "Depth": row["Depth"],
            "VAF(%)": row["VAF(%)"],
            "Type": mut_type,
            "ClinVar": row["CLIN_SIG"],
            "Consequence": row["Consequence"],
            "AltDepth": row["AltDepth"],
        })

        df = pd.DataFrame(target_info)
        df.to_excel(f'{out_prefix}.结果回报表.xlsx', index=False)


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
