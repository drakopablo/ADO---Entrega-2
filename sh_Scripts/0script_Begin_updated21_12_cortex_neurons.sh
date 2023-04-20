#!/bin/bash

#### RUN directly from SCRIPT FOLDER
fecha=$(date --rfc-3339=date ) #Ejecuta date y lo almacena en la variable fecha


logfile=LoganalysisRNA_${fecha}.log

echo it has begun at $(date)

echo .| echo .
#Preparing directory and files...
echo Preparing directories and files...

#######
# Abrir 
#directory=./Files_2019_05_STAR #$ ./SCRIPTS/script_Begin.sh 
directory1=~/omicos2/analisis_neurons_cortex

#$ ./script_Begin.sh #Inside of the folder

#crea 2 directorios, FastQCfiles y filesSTAR. Dentro de fileSTAR crea 3 directorios,BAMfiles,FeatureCounts,TDFfiles.
#Dentro de BamQCFiles crea BamQCFiles, logSTAR,toTranscriptome.
mkdir -vp ${directory1}/{FastQCfiles,trimFastQCfiles,TrimFiles,filesSTAR/{BAMfiles/{BamQCFiles,logSTAR,toTranscriptome},FeatureCounts,TDFfiles}}


#${directory}/{FastQCfiles,filesSTAR/{BAMfiles/logSTAR,HTseq,TDFfiles}}

#ruta donde se encuentra el raw data
##### MODIFICAR 
fastqfiles1=~/omicos2/raw/neurons/fastqfiles/cortex_fastq

#version del genoma
genomeversion=95
#tipo de libreria que usará FeatureCounts
libtypeFeatureCounts=0 #-s < intorstring > (isStrandSpecific)  0 (un-stranded), 1 (stranded) and 2 (reversely stranded)
libPE=-p
#https://www.tldp.org/HOWTO/Bash-Prompt-HOWTO/x700.html 
#Counting Files in the Current Directory
#busca los archivos que empiezen por l sitados en el directorio Fastq y luego cuenta el numero de archivos que hay
NofSample1=$(ls -1 $fastqfiles1 | grep -v ^l | wc -l)



########################################
echo .| echo .
echo FASTQC Quality Control 
#bucle sobre todos los archivos .fq.gz ordenados que no estén repetidos
find $fastqfiles1 -type f -name "*.fastq" -printf "%p\n">>salida_temporal2

while read sampleraw;do
nombre=$(basename -s .fastq $sampleraw)
fastqc=$(find $directory1/FastQCfiles -type f -iname "${nombre}_fastqc.html")
echo $nombre
 if [ ! -f "$fastqc" ]
 then
    echo $sampleraw
    fastqc -t 16 $sampleraw -o ${directory1}/FastQCfiles/
 else 
	echo "el fichero ya habia sido procesado\n"
 fi
done<salida_temporal2
echo "comienza el trim_galore"

ls $fastqfiles1 > salida_temporal3
echo $salida_temporal3

while read samples;do
trimmed=$(find ${directory1}/TrimFiles -type f -iname "${samples}*" -printf "%p\n"| head -n 1 )
if [ ! -f "$trimmed" ]
then
describer=$(basename $samples | cut -d '_' -f 1-2 )
	
		#base=$(basename -s fq.gz $sampleraw)
		#describer=$(echo ${base}|sed 's/[0-9]$//')
  	trim_galore --paired ${fastqfiles1}/${describer}_1.fastq ${fastqfiles1}/${describer}_2.fastq -o ${directory1}/TrimFiles
else
	echo "el fichero ${samples} ya había sido procesado por trim-galore"
fi 
done<salida_temporal3

echo "comienza el FASTQC"
##FastQc a los archivos trimeados
for trimfile in ${directory1}/TrimFiles/*.fq;do
nombre=$(basename -s .fq $trimfile)
fastq=$(find $directory1/trimFastQCfiles -type f -iname "${nombre}_fastqc.html")
 if [ ! -f "$fastq" ]
 then 

	 fastqc -t 12 $trimfile -o ${directory1}/trimFastQCfiles
	echo "el fastqc ha terminado para $trimfile"
 else
	echo "el fichero ya había sido precesado por FastQC"
 fi
done
#llamamos al programa STARalingment para que ejecute el alineamiento por STAR 
cd
cd omicos2
echo RNA-seq ANALYSIS has begun at $(date)
#time ./STARalingment_paired.sh $directory1 $fastqfiles1 $genomeversion |& tee -a $logfile
#directory=$1
#fastqfiles=$2tee
#genomeversion=$3
echo STAR ANALYSIS has finished at $(date)

#Perform strand-specific read counting (use '-s 2' if reversely stranded, '-s 1' if forwardly stranded): 
time ./quantSTAR_FeatureCounts_paired.sh $directory1 $genomeversion $libtypeFeatureCounts $libPE |& tee -a $logfile
#directory=$1
#genomeversion=$2
#libtypeFeatureCounts=$3
#NofSample=$4
echo Quantification with Featurecounts has finished at $(date)


echo \n

echo Multiqc Report

multiqc ${directory1}/FastQCfiles/ ${directory1}/TrimFiles ${directory1}/trimFastQCfiles ${directory1}/filesSTAR/BAMfiles/logSTAR/ ${directory1}/filesSTAR/FeatureCounts/ -o $directory1 

echo "This script is sent to the log file, \"$logfile\"." 

