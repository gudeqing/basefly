import pandas as pd


def generate_igv_snapshot_script(table_file, bam_file_pattern=None, img_dir=None, add_chr='chr', min_af=0.0005):
    if bam_file_pattern is None:
        bam_file_pattern = 'D:\\haixi\\2023-EQA-LungCaner\Blood\\ABRA2-{sample}\\{sample}.realigned.bam'
    if img_dir is None:
        img_dir = 'D:\\haixi\\2023-EQA-LungCaner\\Blood\\IGV'
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

    with open('igv.snapshot.scripts.txt', 'w') as f:
        for each in contents:
            f.write(each+'\n')


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
