#!/bin/bash
#QUANTIFICATION

#Preparing directory and files...
echo Preparing directories and files...

directory=$1
genomeversion=$2
libtypeFeatureCounts=$3
#Este script es igual al quantSTAR_FeatureCounts.sh, pero activamos la variable $libPE=-p para indicar al programa que los ficheros .bam son paired end
libPE=$4

#validacion del numero de argumentos
if [ $# -eq 0 ]; then
   echo "No se ha utilizado ningún parámetro"
   exit
elif [ $# -ne 4 ]; then
   echo "Not equal to 4 parameters"
   exit
else
   while [ $# -gt 0 ]; do
      echo $1
      shift
   done
fi
#NO MODIFICAR
gtffile=/home/LAB/lab_glb_2/omicos2/gtf/*.gtf
echo ${gtffile}
############## SubReads -> FeatureCounts
echo .| echo .
echo featureCounts...
 
#featureCounts
featureCounts $libPE -O -s $libtypeFeatureCounts -T 12 -t exon -g gene_id -a $gtffile \
-o ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.txt ${directory}/filesSTAR/BAMfiles/*.bam

NofSample= $(ls -1 ${directory}/filesSTAR/BAMfiles/*.bam | grep -v ^l | wc -l )

echo "Numero de archivos BAM = " $NofSample

cut -f 1,$(seq -s , 7 $((7+$(($NofSample - 1))))) ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.txt > ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.mat 







