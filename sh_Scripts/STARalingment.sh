#!/bin/bash

echo it has begun at $(date)

echo .| echo .
#Preparing directory and files...
echo Preparing directories and files...

##### MODIFICAR 
#argumentos de entrada
directory=$1
fastqfilesB1=$2
genomeversion=$3

#validacion del numero de argumentos, que no sea igual que 0 y que no sea menor que 3
if [ $# -eq 0 ]; then
   echo "No se ha utilizado ningún parámetro"
   exit
elif [ $# -ne 3 ]; then
   echo "Not equal to 3 parameters"
   exit
# bucle while para que imprima el nombre del argumento dado mientras el numero de argumentos sea mayor a 0
else
   while [ $# -gt 0 ]; do
      echo $1
      shift
   done
fi

#NO MODIFICAR

genomeSTAR=/home/LAB/lab_glb_2/omicos2/genoma #ruta genoma STAR 
#genomeigv=/home/lpa/igv/genomes/mm39.genome #ruta genoma igv 
gtffile=/home/LAB/lab_glb_2/omicos2/gtf #ruta para el archivo .gtf del genoma
######################################
echo $genomeSTAR
echo $gtffile
echo .| echo .

#Alineamos los archivos fastq resultantes del trim-galore.
#Alingment with STAR
echo Alingment with STAR? Yes, of course\!\!
cd $directory/TrimFiles
for sample in $directory/TrimFiles/*.fq ; do
  echo ${sample}
   describer=$(basename -s .fastq $sample)
# Ejecución del programa STAR
   	STAR --runThreadN 16 --genomeDir $genomeSTAR --genomeLoad LoadAndKeep --limitBAMsortRAM 200000000000\
      	--readFilesIn $sample \
      	--outFileNamePrefix ${describer}_ \
      	--outSAMtype BAM SortedByCoordinate \
     	--quantMode TranscriptomeSAM GeneCounts 
#movemos los archivos al directorio del análisis
	mv -v ./*.tab $directory/filesSTAR/BAMfiles/logSTAR/
	mv -v ./*.out $directory/filesSTAR/BAMfiles/logSTAR/
	mv -v ./*Aligned.sortedByCoord.out.bam $directory/filesSTAR/BAMfiles/
	mv -v ./*.toTranscriptome.out.bam $directory/filesSTAR/BAMfiles/toTranscriptome/
done

STAR --genomeDir $genomeSTAR --genomeLoad Remove
cd ~/omicos2
