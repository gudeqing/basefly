version 1.0
# docker exec yjwiddler widdler.py run $PWD/wgs_v3.wdl $PWD/fake.input.json -S 192.168.3.5
# docker exec yjwiddler widdler.py run $PWD/scatter_pipeline.wdl $PWD/input.json -S 192.168.3.5
# java -jar /mnt/nas_002/jars/womtool.jar inputs test.wdl > input.json
# java -jar /mnt/nas_002/jars/cromwell.jar run test.wdl --inputs input.js
# 推荐运行方式：cromshell submit wgs_v0.wdl fake.input.json
# 运行日志地址：/mnt/hnas_1001/genarsa_project/wdl/cromwell-executions/wgs_wf/
# 流程报错  /data/genarsa_project/wdl/logs/yjcromwell.log
# * 对大部分流程最开始的输入文件如fasta、fastq等，都改用String声明，从而避免复制，这得益于singularity运行时提供了挂载目录的功能
# 参考snakemake版本流程/mnt/nas_101/hongxhe/Pipeline/GenarsA_Test/wgs_v1/
# String singularity = "/data/singularity_images_wgsgvcf/bbmap_38.93--he522d1c_0"
# String singularity = "/data/singularity_images_wgsgvcf/bwa_samtools_bcftools_v1.simg"
# String singularity = "/data/singularity_images_wgsgvcf/ensembl-vep:104.3--pl5262h4a94de4_0"
# String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
# String singularity = "/data/singularity_images_wgsgvcf/qualimap:2.2.2d--hdfd78af_2"
# String singularity = "/data/singularity_images_wgsgvcf/sambamba:0.8.1--h41abebc_0"
# String singularity = "/data/singularity_images_wgsgvcf/samtools:1.13--h8c37831_0"

# 注意： call的流程中尽量不要使用默认函数来处理变量，尽量把变量处理的部分放入到task中
# 注意： 对于task的input中尽量不要加工变量，可以改为在command中加工，这样录入数据库的时候可以少定义一些输入参数


workflow wgs_wf {
    input{
        # 结果输出目录
        String outdir = "/mnt/nas_101/genarsa/Pipeline_Test/gdq_0523"
        # queue是特殊字符串，用于提交任务所需，可由IT提供
        String queue = "zryh"
        # 样本信息：样本名称,同一个样本不同测序批次的编号,测序平台类型,read1路径，read2路径
        File inputsamplefile = "/mnt/nas_101/bioworkertest/gdq_0519/in_files/units_WDL.tsv"
        # 样本名称信息列表文件
        File inputsamplenamefile = "/mnt/nas_101/bioworkertest/gdq_0519/in_files/samples_WDL.tsv"
        # contig信息，用于并行分区call突变
        File contigfile = "/mnt/nas_101/bioworkertest/gdq_0519/in_files/contig.tsv"

        # 参考基因组
#        String ref_fa = "/data/reference_genome/hg38/hg38.fa"
#        String ref_dict = "/data/reference_genome/hg38/hg38.dict"
#        String ref_fai = "/data/reference_genome/hg38/hg38.fa.fai"
#        String ref_bwt = "/data/reference_genome/hg38/hg38.fa.bwt"
#        String ref_amb = "/data/reference_genome/hg38/hg38.fa.amb"
#        String ref_ann = "/data/reference_genome/hg38/hg38.fa.ann"
#        String ref_pac = "/data/reference_genome/hg38/hg38.fa.pac"
#        String ref_sa = "/data/reference_genome/hg38/hg38.fa.sa"

        # 匹配目标bam和bai文件的脚本路径
#        File select_array_script = "/mnt/nas_101/genarsa/Pipeline_Test/20230510/script/select_array.py"

        # 突变类型 Possible values: {NO_VARIATION, SNP, MNP, INDEL,SYMBOLIC, MIXED}
        Array[String] variant_types =  ["SNP", "INDEL"]

        # 对输入变量进一步加工
        Array[Array[String]] inputsamples = read_tsv(inputsamplefile)
        Array[Array[String]] inputsamplenames = read_tsv(inputsamplenamefile)
        Array[Array[String]] contigs = read_tsv(contigfile)

    }

    scatter (inputsample in inputsamples) {
        call fastqc {
            input:
                read1=inputsample[3],
                read2=inputsample[4],
                outdir=outdir,
                queue=queue
        }

        call bbduk {
            input:
                sample=inputsample[0],
                unit=inputsample[1],
                Fastq1=inputsample[3],
                Fastq2=inputsample[4],
                outdir=outdir,
                queue=queue
        }

        call map_reads {
            input:
                sample=inputsample[0],
                unit=inputsample[1],
                Fastq1=bbduk.R1,
                Fastq2=bbduk.R2,
                outdir=outdir,
                queue=queue,
        }

        call mark_duplicates {
            input:
                sample=inputsample[0],
                unit=inputsample[1],
                sort_bam=map_reads.bam,
                outdir=outdir,
                queue=queue
        }

        call samtools_stats {
            input:
                input_bam = mark_duplicates.bam,
                sample_name = inputsample[0],
                sample_unit = inputsample[1],
                outdir = outdir,
                queue=queue
        }
    }

    scatter (inputsamplename in inputsamplenames) {
        call merge_bam {
            input:
                dedup_bam_array=mark_duplicates.bam,
                dedup_bambi_array=mark_duplicates.bambi,
                sample=inputsamplename[0],
                outdir=outdir,
                queue=queue
        }
   }

    scatter (inputsample in inputsamples) {
       call qualimap {
           input:
               merge_bam_array=merge_bam.bam,
               sample=inputsample[0],
               unit=inputsample[1],
               outdir=outdir,
               queue=queue
       }
    }

    scatter (inputsamplename in inputsamplenames) {
       scatter (contig in contigs) {
            call call_variants {
                input:
                    merge_bam_array=merge_bam.bam,
                    sample=inputsamplename[0],
                    outdir=outdir,
                    contig=contig[0],
                    queue=queue,
            }
       }

       call merge_variants_sample_contig {
           input:
               contig_gvcf=call_variants.gvcf,
               sample=inputsamplename[0],
               outdir=outdir,
               queue=queue
       }
   }

    call get_interval_list {
        input:
            outdir = outdir,
            queue = queue
    }

    scatter (interval in get_interval_list.interval_list) {
        call combine_calls {
            input:
                gvcfs = merge_variants_sample_contig.gvcf,
                tbis = merge_variants_sample_contig.tbi,
                interval_list_file = interval,
                outdir = outdir,
                queue = queue
        }

        call genotype_variants {
            input:
                gvcf = combine_calls.gvcf,
                gvcftbi = combine_calls.gvcftbi,
                interval_list_file = interval,
                outdir = outdir,
                queue = queue,
        }
   }

    call merge_variants_interval {
        input:
            input_vcfs = genotype_variants.vcf,
            outdir = outdir,
            queue = queue
    }

    scatter (contig in contigs) {
        scatter (vartype in variant_types) {
            call select_calls {
                input:
                    queue = queue,
                    input_vcf = merge_variants_interval.output_vcf,
                    input_vcf_index = merge_variants_interval.output_vcf_index,
                    outdir = outdir,
                    vartype = vartype,
                    contig = contig[0],
            }

            call hard_filter_calls {
                input:
                    queue = queue,
                    input_vcf = select_calls.output_vcf,
                    input_vcf_index = select_calls.output_vcf_index,
                    vartype = vartype,
                    contig = contig[0],
                    outdir = outdir,
            }

            call sort_vcf {
                input:
                    queue = queue,
                    input_vcf = hard_filter_calls.output_vcf,
                    vartype = vartype,
                    contig = contig[0],
                    outdir = outdir
            }
        }
    }

    call merge_calls {
        input:
            queue = queue,
            input_vcfs = sort_vcf.output_vcf,
            outdir = outdir
    }

    call annotate_variants {
        input:
            queue = queue,
            input_file = merge_calls.output_vcf,
            input_file_idx = merge_calls.output_vcf_index,
            outdir = outdir
    }

    call multiqc {
        input:
            fastqc_result_dir = fastqc.result_dir,
            samtools_stats_dir = samtools_stats.result_dir,
            trimmed_result_dir = bbduk.result_dir,
            qualimap_result_dir = qualimap.result_dir,
            dedup_result_dir = mark_duplicates.result_dir,
            queue = queue,
            outdir = outdir
    }
}


# ---task defined below---

task bbduk {
    input {
        String Fastq1
        String Fastq2
        String sample
        String unit
        String queue
        String outdir
        Int cpus = 4
        String mem_mb = 10000
        File adapters = "/mnt/nas_101/genarsa/Pipeline_Test/20230420/adapters.fa"
        String xmx = "8g"
        Int threads = 4
        String singularity = "/data/singularity_images_wgsgvcf/bbmap_38.93--he522d1c_0"
    }

    runtime{
        cpus: cpus
        memory: mem_mb
        queue: queue
    }

    command<<<
            singularity exec -B /mnt,/data ~{singularity} bbduk.sh \
            -Xmx~{xmx} \
            in1=~{Fastq1} in2=~{Fastq2} \
            out1=~{outdir}/result/trimmed/~{sample}-~{unit}.1.fastq.gz \
            out2=~{outdir}/result/trimmed/~{sample}-~{unit}.2.fastq.gz \
            threads=~{threads} \
            ref=~{adapters} \
            k=25 mink=8 ktrim=r \
            ordered=true \
            qtrim=rl hdist=1 stats=~{outdir}/result/trimmed/~{sample}-~{unit}.stats.txt

    >>>

    output{
        File R1 = "~{outdir}/result/trimmed/~{sample}-~{unit}.1.fastq.gz"
        File R2 = "~{outdir}/result/trimmed/~{sample}-~{unit}.2.fastq.gz"
        File stats = "~{outdir}/result/trimmed/~{sample}-~{unit}.stats.txt"
        String result_dir = "~{outdir}/result/trimmed/"
    }

}


task map_reads {
     input{
        File Fastq1
        File Fastq2
        String sample
        String unit
        String queue
        String outdir
        Int cpus = 4
        Int mem_mb = 10000
        Int threads = 12
        String genome_fa = "/data/reference_genome/hg38/hg38.fa"
#        String ref_dict = "/data/reference_genome/hg38/hg38.dict"
#        String ref_fai = "/data/reference_genome/hg38/hg38.fa.fai"
#        String ref_bwt = "/data/reference_genome/hg38/hg38.fa.bwt"
#        String ref_amb = "/data/reference_genome/hg38/hg38.fa.amb"
#        String ref_ann = "/data/reference_genome/hg38/hg38.fa.ann"
#        String ref_pac = "/data/reference_genome/hg38/hg38.fa.pac"
#        String ref_sa = "/data/reference_genome/hg38/hg38.fa.sa"
        String singularity = "/data/singularity_images_wgsgvcf/bwa_samtools_bcftools_v1.simg"
     }


    runtime {
        cpus: cpus
        memory: mem_mb
        queue: queue
    }


     command<<<
        mkdir -p ~{outdir}/result/mapped && \
        singularity exec -B /mnt,/data ~{singularity} bwa mem \
            -M -Y \
            -R '@RG\tID:~{sample}.L1-B1\tSM:~{sample}\tLB:L1-B1\tPL:ILLUMINA' \
            -t ~{threads} \
            ~{genome_fa} ~{Fastq1} ~{Fastq2} > ~{outdir}/result/mapped/~{sample}-~{unit}.unsorted.bam

         singularity exec -B /mnt,/data ~{singularity} samtools sort \
            --threads ~{threads} -OBAM \
         ~{outdir}/result/mapped/~{sample}-~{unit}.unsorted.bam \
         -o ~{outdir}/result/mapped/~{sample}-~{unit}.sorted.bam

         rm ~{outdir}/result/mapped/~{sample}-~{unit}.unsorted.bam

     >>>

     output{
        File bam = "~{outdir}/result/mapped/~{sample}-~{unit}.sorted.bam"
     }
}


task mark_duplicates{
    input {
        File sort_bam
        String sample
        String unit
        String outdir
        String queue
        Int cpus = 8
        Int mem_mb = 10000
        Int threads = 8
        String singularity = "/data/singularity_images_wgsgvcf/sambamba:0.8.1--h41abebc_0"
    }

    runtime {
        cpus: cpus
        memory: mem_mb
        queue: queue
    }

    command<<<
        mkdir -p ~{outdir}/result/dedup && touch ~{outdir}/result/dedup/~{sample}-~{unit}.log.txt \
        && singularity exec -B /mnt,/data ~{singularity} sambamba \
        markdup -t ~{threads} ~{sort_bam} ~{outdir}/result/dedup/~{sample}-~{unit}.bam 2>&1 | awk 'NR>2' > ~{outdir}/result/dedup/~{sample}-~{unit}.log.txt
    >>>

    output{
        File bam = "~{outdir}/result/dedup/~{sample}-~{unit}.bam"
        File bambi = "~{outdir}/result/dedup/~{sample}-~{unit}.bam.bai"
        File logtxt = "~{outdir}/result/dedup/~{sample}-~{unit}.log.txt"
        String result_dir = "~{outdir}/result/dedup"
    }

}


task merge_bam{
    input {
        Array[File] dedup_bam_array
        Array[File] dedup_bambi_array
        String sample
        String outdir
        String queue
        File select_array_script = "/mnt/nas_101/genarsa/Pipeline_Test/20230510/script/select_array.py"
        Int threads = 8
        String singularity = "/data/singularity_images_wgsgvcf/samtools:1.13--h8c37831_0"
    }

    runtime{
        queue: queue
    }

    command <<<

        dedup_bam=`/data/genarsa_project/bioinfo_miniconda3/envs/snakemake_python3.9.0/bin/python ~{select_array_script} ~{sep="," dedup_bam_array} ~{sep="," dedup_bambi_array} ~{sample} bam` && \
        dedup_bambi=`/data/genarsa_project/bioinfo_miniconda3/envs/snakemake_python3.9.0/bin/python  ~{select_array_script} ~{sep="," dedup_bam_array} ~{sep="," dedup_bambi_array} ~{sample} bambai` && \
        n_bam=`echo ${dedup_bam}|awk '{print NF}'` && \
        if (( $n_bam == 1 ));then \
            if [ -f "~{outdir}/result/dedup/~{sample}.merged.bam" ];then
                rm ~{outdir}/result/dedup/~{sample}.merged.bam ~{outdir}/result/dedup/~{sample}.merged.bam.bai
            fi
            ln -s ${dedup_bam} ~{outdir}/result/dedup/~{sample}.merged.bam && ln -s ${dedup_bambi} ~{outdir}/result/dedup/~{sample}.merged.bam.bai;\
        else \
            singularity exec -B /mnt,/data ~{singularity} samtools merge -f -o ~{outdir}/result/dedup/~{sample}.merged.bam --threads ~{threads} --write-index ${dedup_bam}; \
        fi

    >>>

    output{
        File bam = "~{outdir}/result/dedup/~{sample}.merged.bam"
    }
}


task qualimap{
    input{
        Array[File] merge_bam_array
        String sample
        String unit
        String outdir
        String queue
        Int cpus = 4
        Int mem_mb = 16000
        File select_array_script = "/mnt/nas_101/genarsa/Pipeline_Test/20230510/script/select_array.py"
        Int threads = 4
        String mem = "15000M"
        String singularity = "/data/singularity_images_wgsgvcf/qualimap:2.2.2d--hdfd78af_2"
    }

    runtime {
        cpus: cpus
        memory: mem_mb
        queue: queue
    }

    command <<<

     merge_bam=`/data/genarsa_project/bioinfo_miniconda3/envs/snakemake_python3.9.0/bin/python ~{select_array_script} ~{sep="," merge_bam_array} - ~{sample} bam` && \
     unset DISPLAY; singularity exec -B /mnt,/data ~{singularity} qualimap --java-mem-size=~{mem} bamqc -outformat HTML -outdir ~{outdir}/result/qc/qualimap/~{sample}-~{unit} -bam ${merge_bam} -nt ~{threads}

    >>>

    output{
        File genome_result = "~{outdir}/result/qc/qualimap/~{sample}-~{unit}/genome_results.txt"
        String result_dir = "~{outdir}/result/qc/qualimap/"

    }

}


task call_variants{
    input {
           Array[File] merge_bam_array
           String sample
           String outdir
           String contig
           String queue
           Int mem_mb = 10000
           Int cpus = 4
           File select_array_script = "/mnt/nas_101/genarsa/Pipeline_Test/20230510/script/select_array.py"
           String genome_fa = "/data/reference_genome/hg38/hg38.fa"
#           String ref_dict = "/data/reference_genome/hg38/hg38.dict"
#           String ref_fai = "/data/reference_genome/hg38/hg38.fa.fai"
           String java_opts = "-Xmx10g"
           String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
    }

    runtime {
        cpus: cpus
        memory: mem_mb
        queue: queue
    }

    command <<<
      merge_bam=`/data/genarsa_project/bioinfo_miniconda3/envs/snakemake_python3.9.0/bin/python ~{select_array_script} ~{sep="," merge_bam_array} - ~{sample} bam` && \
      mkdir -p ~{outdir}/result/called && \
        singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' \
        HaplotypeCaller \
            --annotate-with-num-discovered-alleles true \
            -R ~{genome_fa} \
            -I ${merge_bam} \
            -ERC GVCF \
            -O ~{outdir}/result/called/~{sample}.~{contig}.g.vcf.gz \
            -L ~{contig}

    >>>


    output{
        File gvcf = "~{outdir}/result/called/~{sample}.~{contig}.g.vcf.gz"
        File tbi = "~{outdir}/result/called/~{sample}.~{contig}.g.vcf.gz.tbi"

    }


}


task merge_variants_sample_contig{
    input {
       Array[File] contig_gvcf
       String sample
       String outdir
       String queue
       Int cpus = 4
       Int mem_mb = 10000
       String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
       String java_opts = "-Xmx10g"
    }

    runtime {
      cpus: cpus
      memory: mem_mb
      queue: queue
    }

    command <<<

       singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' MergeVcfs \
            ~{sep=" " prefix('INPUT=',contig_gvcf)} \
            OUTPUT=~{outdir}/result/called/~{sample}.g.vcf.gz \

    >>>

    output{

        File gvcf = "~{outdir}/result/called/~{sample}.g.vcf.gz"
        File tbi = "~{outdir}/result/called/~{sample}.g.vcf.gz.tbi"

    }

}


task get_interval_list {
    input {
#        String interval_dir = "/data/reference_genome/hg38/interval-files-folder"
        String interval_dir = "/mnt/nas_101/bioworkertest/gdq_0519/in_files/intervals"
        String queue
        String outdir
    }

    runtime {
        queue: queue
    }

    command <<<
        ls ~{interval_dir}/*-scattered.interval_list
    >>>

    output {
        Array[String] interval_list = read_lines(stdout())
    }
}


task combine_calls{
    input {
        Array[File] gvcfs
        Array[File] tbis
        File interval_list_file
        String interval_idx = sub(basename(interval_list_file), "-scattered.interval_list", "")
        String outdir
        String queue
        Int mem_mb = 10000
        Int cpus = 7
        String genome_fa = "/data/reference_genome/hg38/hg38.fa"
        # String ref_dict = "/data/reference_genome/hg38/hg38.dict"
        # String ref_fai = "/data/reference_genome/hg38/hg38.fa.fai"
        String java_opts = "-Xmx10g"
        String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
    }

    runtime {
       cpus: cpus
       memory: mem_mb
       queue: queue
    }

    command <<<
      singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' CombineGVCFs \
            -R ~{genome_fa} \
            -O ~{outdir}/result/called/all.~{interval_idx}.g.vcf.gz\
            ~{sep=" " prefix('-V ',gvcfs)} \
            -L ~{interval_list_file}
    >>>

    output{
        File gvcf="~{outdir}/result/called/all.~{interval_idx}.g.vcf.gz"
        File gvcftbi="~{outdir}/result/called/all.~{interval_idx}.g.vcf.gz.tbi"
    }
}


task genotype_variants{
    input {
        File gvcf
        File gvcftbi
        File interval_list_file
        String interval_idx = sub(basename(interval_list_file), "-scattered.interval_list", "")
        String outdir
        String queue
        Int cpus = "4"
        String genome_fa = "/data/reference_genome/hg38/hg38.fa"
        # String ref_dict = "/data/reference_genome/hg38/hg38.dict"
        # String ref_fai = "/data/reference_genome/hg38/hg38.fa.fai"
        String dbsnp = "/data/reference_genome/hg38/gatk/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
        # String dbsnp_idx = "/data/reference_genome/hg38/gatk/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
        String java_opts = "-Xmx10g"
        String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
    }

    runtime{

        cpus: cpus
        queue: queue

    }

    command <<<

       mkdir -p ~{outdir}/result/genotyped && singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' GenotypeGVCFs \
            -V ~{gvcf} \
            -R ~{genome_fa} \
            -O ~{outdir}/result/genotyped/all.~{interval_idx}.vcf.gz \
            --dbsnp ~{dbsnp} \

    >>>

    output{
        File vcf = "~{outdir}/result/genotyped/all.~{interval_idx}.vcf.gz"
        File vcf_tbi = "~{outdir}/result/genotyped/all.~{interval_idx}.vcf.gz.tbi"
    }
}


task merge_variants_interval {
    input {
        Array[File] input_vcfs
        String outdir
#        String output_filename = "~{outdir}/result/genotyped/all.vcf.gz"
        Int cpus = 4
        Int mem_mb = 10000
        String queue
        String java_opts = "-Xmx10g"
        String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    command <<<
        singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' \
        MergeVcfs \
        --INPUT ~{sep=' --INPUT ' input_vcfs} \
        --OUTPUT ~{outdir}/result/genotyped/all.vcf.gz
    >>>

    output {
        File output_vcf = "~{outdir}/result/genotyped/all.vcf.gz"
        File output_vcf_index = "~{outdir}/result/genotyped/all.vcf.gz.tbi"
  }
}


task select_calls {
    input {
        File input_vcf
        File input_vcf_index
        String vartype
        String contig
        String outdir
#        String output_filename = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.vcf.gz"
        String genome_fa = "/data/reference_genome/hg38/hg38.fa"
#        String ref_dict = "/data/reference_genome/hg38/hg38.dict"
#        String ref_fai = "/data/reference_genome/hg38/hg38.fa.fai"
        Int cpus = 4
        Int mem_mb = 10000
        String queue
        String java_opts = "-Xmx10g"
        String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    command <<<
        mkdir -p ~{outdir}/result/filtered
        singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' \
        SelectVariants \
        -select-type ~{vartype} \
        -R ~{genome_fa} \
        -V ~{input_vcf} \
        -L ~{contig} \
        -O ~{outdir}/result/filtered/all.~{contig}.~{vartype}.vcf.gz
    >>>

    output {
        File output_vcf = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.vcf.gz"
        File output_vcf_index = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.vcf.gz.tbi"
  }

}


task hard_filter_calls {
    input {
        File input_vcf
        File input_vcf_index
        String vartype = "SNP"
        String contig
        # hard-filter: https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
        String indel_hard_filter = '--filter-name "snv-hard-filter" --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"'
        String snv_hard_filter = '--filter-name "indel-hard-filter" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"'
        # 下面的参数加工命令没有办法加入到command中，因为会报错
        String filter_str = if vartype == "SNP" then snv_hard_filter else indel_hard_filter
        String outdir
#        String output_filename = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.hardfiltered.vcf.gz"
        String genome_fa = "/data/reference_genome/hg38/hg38.fa"
#        String ref_dict = "/data/reference_genome/hg38/hg38.dict"
#        String ref_fai = "/data/reference_genome/hg38/hg38.fa.fai"
        Int cpus = 4
        Int mem_mb = 10000
        String queue
        String java_opts = "-Xmx10g"
        String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    command <<<
        singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' \
        VariantFiltration \
        ~{filter_str} \
        -R ~{genome_fa} \
        -V ~{input_vcf} \
        -O ~{outdir}/result/filtered/all.~{contig}.~{vartype}.hardfiltered.vcf.gz
    >>>

    output {
        File output_vcf = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.hardfiltered.vcf.gz"
        File output_vcf_index = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.hardfiltered.vcf.gz.tbi"
  }

}


task sort_vcf {
    input {
        File input_vcf
        String vartype
        String contig
        String outdir
#        String output_filename = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.hardfiltered.sort.vcf.gz"
        Int cpus = 4
        Int mem_mb = 10000
        String queue
        String java_opts = "-Xmx10g"
        String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    command <<<
        singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' \
        SortVcf \
        --MAX_RECORDS_IN_RAM 20000 \
        -I ~{input_vcf} \
        -O ~{outdir}/result/filtered/all.~{contig}.~{vartype}.hardfiltered.sort.vcf.gz
    >>>

    output {
        File output_vcf = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.hardfiltered.sort.vcf.gz"
        File output_vcf_index = "~{outdir}/result/filtered/all.~{contig}.~{vartype}.hardfiltered.sort.vcf.gz.tbi"
  }

}


task merge_calls {
    input {
        Array[Array[File]] input_vcfs
        String outdir
#        String output_filename = "~{outdir}/result/filtered/all.vcf.gz"
        Int cpus = 4
        Int mem_mb = 10000
        String queue
        String java_opts = "-Xmx10g"
        String singularity = "/data/singularity_images_wgsgvcf/gatk4:4.2.2.0--hdfd78af_1"
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    command <<<
        singularity exec -B /mnt,/data ~{singularity} gatk --java-options '~{java_opts}' \
        MergeVcfs \
        --INPUT ~{sep=' --INPUT ' flatten(input_vcfs)} \
        -O ~{outdir}/result/filtered/all.vcf.gz
    >>>

    output {
        File output_vcf = "~{outdir}/result/filtered/all.vcf.gz"
        File output_vcf_index = "~{outdir}/result/filtered/all.vcf.gz.tbi"
  }

}


task annotate_variants {
    input {
        File input_file
        File? input_file_idx
        String genome_fa = "/data/reference_genome/hg38/hg38.fa"
        String outdir
        String output_file = "~{outdir}/result/annotated/all.vcf.gz"
        String stats_file = "~{outdir}/result/annotated/all.stats.html"
        String output_format = "vcf"
        Boolean output_vcf = output_format == "vcf"
        String compress_output = "bgzip"
        Boolean force_overwrite = true
        Int fork = 4
        # 直接给字符串路径作为输入，因singularity可以将数据盘挂载到执行环境
        String? vep_plugin_dir
        String vep_cache_dir = "/data/reference_genome/biodb/vep/grch38/"
#        Array[String] plugin_names = ['Frameshift', 'Wildtype']
        Array[String] plugin_names = []
        String species = "homo_sapiens"
        String assembly_version = "GRCh38" # 参考现有流程设置
        Boolean cache = true
        Boolean offline = true
        Boolean merged = true # 参考现有流程设置
        Boolean variant_class = false # 参考现有流程设置
        String sift = "b"
        String polyphen = "b"
        String nearest = "transcript"
        Boolean gene_phenotype = false # 参考现有流程设置
        Boolean regulatory = false # 参考现有流程设置
        Boolean phased = false # 参考现有流程设置
        Boolean numbers = false # 参考现有流程设置
        Boolean hgvs = true
        Boolean transcript_version = false # 参考现有流程设置
        Boolean symbol = false # 参考现有流程设置
        Boolean tsl = false # 参考现有流程设置
        Boolean canonical = false # 参考现有流程设置
        Boolean biotype = false # 参考现有流程设置
        Boolean mane = false # 参考现有流程设置
        Boolean max_af = false # 参考现有流程设置
        Boolean af_1kg = false # 参考现有流程设置
        Boolean af_gnomad = false # 参考现有流程设置
        Boolean af_esp = false
        Boolean coding_only = false
        Boolean pick = false
        Boolean flag_pick = false # 参考现有流程设置
        Boolean filter_common = false # 参考现有流程设置
        # for runtime
        String singularity = "/data/singularity_images_wgsgvcf/ensembl-vep:104.3--pl5262h4a94de4_0"
        Int cpus = fork + 1
        Int mem_mb = 10000
        String queue
    }

    command <<<
        set -e
        mkdir -p ~{outdir}/result/annotated/

        singularity exec -B /mnt,/data ~{singularity} vep \
        ~{"-i " + input_file} \
        ~{"--fasta " + genome_fa} \
        ~{"-o " + output_file} \
        ~{"--" + output_format} \
        ~{"--compress_output " + compress_output} \
        ~{if force_overwrite then "--force_overwrite  " else ""} \
        ~{"--fork " + fork} \
        ~{"--species " + species} \
        ~{"--assembly " + assembly_version} \
        --dir_cache ~{vep_cache_dir} \
        ~{if defined(vep_plugin_dir) then "--dir_plugins " + vep_plugin_dir else ""} \
        ~{sep=" " if length(plugin_names) > 0 then prefix("--plugin ", plugin_names) else [""]} \
        ~{"--stats_file " + stats_file} \
        ~{if cache then "--cache  " else ""} \
        ~{if offline then "--offline  " else ""} \
        ~{if merged then "--merged  " else ""} \
        ~{if variant_class then "--variant_class  " else ""} \
        ~{"--sift " + sift} \
        ~{"--polyphen " + polyphen} \
        ~{"--nearest " + nearest} \
        ~{if gene_phenotype then "--gene_phenotype  " else ""} \
        ~{if regulatory then "--regulatory  " else ""} \
        ~{if phased then "--phased  " else ""} \
        ~{if numbers then "--numbers  " else ""} \
        ~{if hgvs then "--hgvs  " else ""} \
        ~{if transcript_version then "--transcript_version  " else ""} \
        ~{if symbol then "--symbol  " else ""} \
        ~{if tsl then "--tsl  " else ""} \
        ~{if canonical then "--canonical  " else ""} \
        ~{if biotype then "--biotype  " else ""} \
        ~{if mane then "--mane " else ""} \
        ~{if max_af then "--max_af  " else ""} \
        ~{if af_1kg then "--af_1kg  " else ""} \
        ~{if af_gnomad then "--af_gnomad  " else ""} \
        ~{if af_esp then "--af_esp  " else ""} \
        ~{if coding_only then "--af_esp  " else ""} \
        ~{if pick then "--pick " else ""} \
        ~{if flag_pick then "--flag_pick  " else ""} \
        ~{if filter_common then "--filter_common  " else ""}
        # create index
        ~{if output_vcf then "singularity exec -B /mnt,/data " + singularity + " tabix " + output_file else ""}
    >>>

    output {
        File out_vcf = "~{output_file}"
        File? out_vcf_idx = "~{output_file}.tbi"
        File stats_file = "~{stats_file}"
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    meta {
        name: "VEP"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        input_file: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input file"}
        genome_fa: {prefix: "--fasta ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run VEP with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache)."}
        output_file: {prefix: "-o ", type: "str", level: "required", default: "tumor.vep.vcf.gz", range: "None", array: "False", desc: "output file"}
        output_format: {prefix: "--", type: "str", level: "required", default: "vcf", range: "{'vcf', 'json', 'tab'}", array: "False", desc: "If we choose to write output in VCF format. Consequences are added in the INFO field of the VCF file, using the key 'CSQ'. Data fields are encoded separated by '|'; the order of fields is written in the VCF header. Output fields in the 'CSQ' INFO field can be selected by using --fields."}
        compress_output: {prefix: "--compress_output ", type: "str", level: "required", default: "bgzip", range: "None", array: "False", desc: "Writes output compressed using either gzip or bgzip"}
        force_overwrite: {prefix: "--force_overwrite ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Force overwriting of output file"}
        fork: {prefix: "--fork ", type: "int", level: "required", default: "4", range: "None", array: "False", desc: "Use forking(multi-cpu/threads) to improve script runtime"}
        species: {prefix: "--species ", type: "str", level: "required", default: "homo_sapiens", range: "None", array: "False", desc: "Species for your data. This can be the latin name e.g. homo_sapiens or any Ensembl alias e.g. mouse."}
        assembly_version: {prefix: "--assembly ", type: "str", level: "required", default: "GRCh37", range: "None", array: "False", desc: "Select the assembly version to use if more than one available."}
        vep_cache_dir: {prefix: "--dir_cache ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Specify the cache files of VEP to use"}
        vep_plugin_dir: {prefix: "--dir_plugins ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Specify the plugin files to use"}
        stats_file: {prefix: "--stats_file ", type: "str", level: "required", default: "tumor.vep.summary.html", range: "None", array: "False", desc: "Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>."}
        cache: {prefix: "--cache ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Enables use of cache"}
        offline: {prefix: "--offline ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation"}
        merged: {prefix: "--merged ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used."}
        plugin_names: {prefix: "--plugin ", type: "str", level: "required", default: "['Frameshift', 'Wildtype']", range: "None", array: "False", desc: "Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times"}
        variant_class: {prefix: "--variant_class ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Output the Sequence Ontology variant class."}
        sift: {prefix: "--sift ", type: "str", level: "required", default: "b", range: "{'b', 's', 'p'}", array: "False", desc: "Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both."}
        polyphen: {prefix: "--polyphen ", type: "str", level: "required", default: "b", range: "{'b', 's', 'p'}", array: "False", desc: "Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both."}
        nearest: {prefix: "--nearest ", type: "str", level: "required", default: "transcript", range: "{'gene', 'transcript', 'symbol'}", array: "False", desc: "Retrieve the transcript or gene with the nearest protein-coding transcription start site (TSS) to each input variant. Use transcript to retrieve the transcript stable ID, gene to retrieve the gene stable ID, or symbol to retrieve the gene symbol. Note that the nearest TSS may not belong to a transcript that overlaps the input variant, and more than one may be reported in the case where two are equidistant from the input coordinates."}
        gene_phenotype: {prefix: "--gene_phenotype ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Indicates if the overlapped gene is associated with a phenotype, disease or trait."}
        regulatory: {prefix: "--regulatory ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature."}
        phased: {prefix: "--phased ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data."}
        numbers: {prefix: "--numbers ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds affected exon and intron numbering to to output. Format is Number/Total"}
        hgvs: {prefix: "--hgvs ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate."}
        transcript_version: {prefix: "--transcript_version ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add version numbers to Ensembl transcript identifiers"}
        symbol: {prefix: "--symbol ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the gene symbol (e.g. HGNC) (where available) to the output."}
        tsl: {prefix: "--tsl ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the transcript support level for this transcript to the output."}
        canonical: {prefix: "--canonical ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds a flag indicating if the transcript is the canonical transcript for the gene."}
        biotype: {prefix: "--biotype ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the biotype of the transcript or regulatory feature."}
        mane: {prefix: "--mane ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds a flag indicating if the transcript is the MANE Select or MANE Plus Clinical transcript for the gene."}
        max_af: {prefix: "--max_af ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD."}
        af_1kg: {prefix: "--af_1kg ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output."}
        af_gnomad: {prefix: "--af_gnomad ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included"}
        af_esp: {prefix: "--af_esp ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Include allele frequency from NHLBI-ESP populations."}
        coding_only: {prefix: "--af_esp ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Only return consequences that fall in the coding regions of transcripts. Not used by default"}
        pick: {prefix: "--pick", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Pick one line or block of consequence data per variant, including transcript-specific columns. This is the best method to use if you are interested only in one consequence per variant"}
        flag_pick: {prefix: "--flag_pick ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others."}
        filter_common: {prefix: "--filter_common ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters."}
    }

}


task fastqc {
    input {
        String read1
        String read2
        Int cpus = 6
        Int mem_mb = 10000
        String queue
        String outdir
        String? other_args = ""
#        String result_path = "~{outdir}/result/qc/fastqc"
        String singularity = "/data/singularity_images_wgsgvcf/fastqc:0.11.9--hdfd78af_1"
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    command <<<
        set -e
        mkdir -p ~{outdir}/result/qc/fastqc
        singularity exec -B /mnt,/data ~{singularity} fastqc \
        --quiet \
        ~{other_args} \
        -t ~{cpus} \
        --outdir ~{outdir}/result/qc/fastqc \
        ~{read1} \
        ~{read2}
    >>>

    output {
        String result_dir = "~{outdir}/result/qc/fastqc"
    }
}


task samtools_stats {
    input {
        File input_bam
        Int cpus = 1
        Int mem_mb = 10000
        String queue
        String outdir
        String sample_name
        String sample_unit
#        String output_filename = "~{outdir}/result/qc/samtools_stats/~{sample_name}-~{sample_unit}.txt"
        String singularity = "/data/singularity_images_wgsgvcf/samtools:1.13--h8c37831_0"
        String? other_args = ""
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    command <<<
        set -e
        mkdir -p ~{outdir}/result/qc/samtools_stats/
        singularity exec -B /mnt,/data ~{singularity} \
        samtools stats ~{other_args} ~{input_bam} > ~{outdir}/result/qc/samtools_stats/~{sample_name}-~{sample_unit}.txt
    >>>

    output {
        File bam_stats_file = "~{outdir}/result/qc/samtools_stats/~{sample_name}-~{sample_unit}.txt"
        String result_dir = "~{outdir}/result/qc/samtools_stats/"
    }
}


task multiqc {
    input {
        String outdir
        Array[String] fastqc_result_dir
        Array[String] samtools_stats_dir
        Array[String] trimmed_result_dir
        Array[String] qualimap_result_dir
        Array[String] dedup_result_dir
        Int cpus = 1
        Int mem_mb = 10000
        String queue
        String singularity = "/data/singularity_images_wgsgvcf/multiqc:1.11--pyhdfd78af_0"
        String? other_args = ""
    }

    runtime {
        cpus: cpus
        queue: queue
        memory: mem_mb
    }

    command <<<
        set -e
        singularity exec -B /mnt,/data ~{singularity} \
        multiqc --force ~{other_args} \
        -o ~{outdir}/result/qc/ \
        -n multiqc \
        ~{fastqc_result_dir[0]} \
        ~{samtools_stats_dir[0]} \
        ~{trimmed_result_dir[0]} \
        ~{qualimap_result_dir[0]} \
        ~{dedup_result_dir[0]}
    >>>

    output {
        File report_html = "~{outdir}/result/qc/multiqc.html"
    }
}


