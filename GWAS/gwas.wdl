version 1.0

workflow gwas_wf{
    input{
        # 输入文件 file_WDL.tsv grouping.txt
        File inputsamplefile
#        File grouping
        Array[Array[File]] inputsamples = read_tsv(inputsamplefile)
#        File Race_file
        # python脚本--流程中需要的脚本
#        File plot_R
#        File anno_assoc
#        File sort_assoc
#        File pca_select
#        File pheno_make
#        File MDS_merge
        # 结果输出路径
        String outdir
        String queue
    }

    scatter (inputsample in inputsamples) {
        call genotype {
            input:
                sample = inputsample[0],
                VCF = inputsample[1],
                genotype_outdir = outdir,
                queue=queue
        }
    }

    output {
        Array[File] out_bed = genotype.bed
    }
}

# -----------------task defined here----------------
task genotype{
    input{
        File VCF
        Int threads = 12
        Int cpu = 6
        String genotype_outdir
        String sample
        String singularity = "/mnt/nas_002/gwas/dailsun/sig_image/plink:1.90b6.21--hec16e2b_2.sif"
        String queue
    }

    runtime {
        cpus: cpu
        queue: queue
    }

    command<<<
            mkdir -p ~{genotype_outdir}/~{sample}/genotype && \
            singularity exec -B /mnt,/data ~{singularity} \
            plink --vcf ~{VCF}  --out ~{genotype_outdir}/~{sample}/genotype/~{sample} --const-fid --threads ~{threads}
    >>>

    output{
        File bed = "~{genotype_outdir}/~{sample}/genotype/~{sample}.bed"
        File bim = "~{genotype_outdir}/~{sample}/genotype/~{sample}.bim"
        File fam = "~{genotype_outdir}/~{sample}/genotype/~{sample}.fam"
    }
}
