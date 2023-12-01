# 使用cutadapt去除primer及primer前面的未知序列，同时丢掉不包含primer的序列
for each in `ls *_1.fastq.gz`; do name=${each%_1*};
 cutadapt -g "file:../primer-F.fa" -G "file:../primer-R.fa" --overlap 10 --pair-adapters --length-tag length= --discard-untrimmed -o ${name}.R1.fq.gz -p ${name}.R2.fq.gz --info-file ../cleanFastqs/${name}.primer-cutInfo.txt -e 0.25 ${name}_1.fastq.gz ${name}_2.fastq.gz > ../cleanFastqs/${name}.trim.log; done
