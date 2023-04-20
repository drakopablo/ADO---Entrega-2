#!/bin/bash
#QUANTIFICATION

#Preparing directory and files...
echo Preparing directories and files...

directory=$1
genomeversion=$2
libtypeFeatureCounts=$3

#libPE=$4  Al estar trabajando con muestras single-end, desactivamos esta opcion
 
#validacion del numero de argumentos
if [ $# -eq 0 ]; then
   echo "No se ha utilizado ningún parámetro"
   exit
elif [ $# -ne 3 ]; then
   echo "Not equal to 3 parameters"
   exit
else
   while [ $# -gt 0 ]; do
      echo $1
      shift
   done
fi
#NO MODIFICAR
#Ruta al archivo gtf del genoma, Este archivo es donde se encuentran las anotaciones de los genes en Ensembl
gtffile=/home/LAB/lab_glb_2/omicos2/gtf/*.gtf
echo ${gtffile}
############## SubReads -> FeatureCounts
echo .| echo .
echo featureCounts...
 
#featureCounts
featureCounts -O -s $libtypeFeatureCounts -T 12 -t exon -g gene_id -a $gtffile \
-o ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.txt ${directory}/filesSTAR/BAMfiles/*.bam

NofSample=$(ls -1 ${directory}/filesSTAR/BAMfiles/*.bam | grep -v ^l | wc -l)
echo "Numero de archivos BAM = " $NofSample

#Modificamos la matriz de cuentas para quitar columnas no necesarias y quedarnos solo con las cuentas de los genes. De esta manera ya tenemos la matriz lista para trabajar con ella en R
cut -f 1,$(seq -s , 7 $((7+$(($NofSample - 1))))) ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.txt > ${directory}/filesSTAR/FeatureCounts/featurecountstotal_STAR.mat 







