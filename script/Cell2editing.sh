#!/bin/bash
#Author------Yan
#Fix------LiuPinBo
#Function------RNA editing was detected in different cell types

if [ $# == 0 ]; then
        echo "Parameter error"
        echo "Eleven Parameter needed!!! Type --help or help for details "
        exit 1;
elif [[ $1 == *help ]] ; then
        echo "Parameter"
        echo "The first parameter is the file path of bam. The second parameter is the barcode list. The third parameter is the pattern of barcode pattern in bam file. The fourth parameter is the name of cell type. The fifth parameter is the path of output. The sixth parameter is the reference genome. The seventh parameter is the SNP database. The eighth parameter is the simpleRepeat.merge.bed file. The ninth parameter is the alu.bed file. The tenth parameter is the file path of Split_bacorde_v2.pl. The eleventh parameter is the file path of duplicate.v2.py. The twelve parameter is the file path of red_ML.pl"
        exit 1;
elif [ $# == 12 ] ; then
        bamfile=$1
        barcodelist=$2
        barcodepattern=$3
        celltype=$4
        outpath=$5
        refseq=$6
        refvcf=$7
        Repeat=$8
        alu=$9
	split="${10}"
	duplicate="${11}"
	red_ML="${12}"
#       echo ${bamfile} ${barcodelist} ${barcodepattern} ${celltype}
        cd ${outpath} &&
	echo "${celltype}"
        mkdir -p ${celltype} &&
        cd ${celltype} &&
        echo ==========Cell2bam start at : `date` ========== &&
        echo "Spliting bam file------------------------------" &&
        echo "perl $split -f1 ${barcodelist} -f2 ${bamfile} -f3 $barcodepattern -out1 ${celltype}.bam &&"
	perl $split -f1 ${barcodelist} -f2 ${bamfile} -f3 $barcodepattern -out1 ${celltype}.bam &&
        mv sorted.bam ${celltype}.sort.bam &&
        mv sorted.bam.bai ${celltype}.sort.bam.bai &&
        echo "index bam file AND get chromosome---------------" &&
#       samtools index ${celltype}.sort.bam &&
        samtools idxstats ${celltype}.sort.bam|awk '{print $1}'|sed 's/*//g'|sed '/^[[:blank:]]*$/d' | head -25 >  chromosome.txt &&
        echo "Marker duplicate, it will take long time!!!-----" &&
	rm -rf run.sh merge.path 
	awk '{print $1}' chromosome.txt| while read i;do echo "/mnt/sda/project/yjc1/miniconda3/envs/htseq/bin/python $duplicate ${celltype}.sort.bam ${celltype}.marked_duplicates.${i}.bam $i" >>run.sh; echo "${celltype}.marked_duplicates.${i}.bam" >>merge.path;done &&
        #for j in $(cat chromosome.txt|awk '{print $1}'); do  echo python $duplicate ${celltype}.sort.bam ${celltype}.marked_duplicates.${j}.bam ${j} >> run.sh ; echo ${celltype}.marked_duplicates.${j}.bam >> merge.path ;done &&
        sh run.sh &&
        echo "merge bam----------------------------------------" &&
        samtools merge -b merge.path ${celltype}.marked_duplicates.bam &&
	export JAVA_HOME=/mnt/sda/project/yjc1/pipeline/02_scRNA_edit/script/jdk-17.0.17 &&
	export CLASSPATH=.:$JAVA_HOME/lib/ &&
	export PATH=.:$JAVA_HOME/bin:$PATH &&
        /mnt/sda/project/yjc1/miniconda3/envs/scRE/bin/picard SortSam -I ${celltype}.marked_duplicates.bam -O ${celltype}.marked_duplicates.sort.bam -SORT_ORDER coordinate &&
        # echo "RNA editing--------------------------------------" &&
        # $red_ML --rnabam ${celltype}.marked_duplicates.sort.bam --reference $refseq --dbsnp $refvcf --simpleRepeat $Repeat --alu $alu --outdir . &&
	echo ==========end at : `date` ========== 

fi

