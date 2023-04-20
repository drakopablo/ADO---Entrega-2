
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
#Ruta y nombre de la carpeta en la que se van a alojar los ficheros resultantes del análisis.
directory1=~/omicos2/analisis_astrocytes


#crea 4 directorios,trimFastQCfiles,TrimFiles, FastQCfiles y filesSTAR. Dentro de fileSTAR crea 3 directorios,BAMfiles,FeatureCounts,TDFfiles.
#Dento de BAMFiles crea BamQCFiles, logSTAR,toTranscriptome.
mkdir -vp ${directory1}/{FastQCfiles,trimFastQCfiles,TrimFiles,filesSTAR/{BAMfiles/{BamQCFiles,logSTAR,toTranscriptome},FeatureCounts,TDFfiles}}


#${directory}/{FastQCfiles,filesSTAR/{BAMfiles/logSTAR,HTseq,TDFfiles}}

#ruta donde se encuentra el raw data
##### MODIFICAR 
fastqfiles1=~/omicos2/raw/astrocytes/fastqfiles

#version del genom
genomeversion=95

#tipo de libreria que usará FeatureCounts
libtypeFeatureCounts=0 #-s < intorstring > (isStrandSpecific)  0 (un-stranded), 1 (stranded) and 2 (reversely stranded)
#Variable que indica a FeatureCounts que el archivo es paired-end 
libPE=-p



########################################
echo .| echo .
echo FASTQC Quality Control 
#bucle sobre todos los archivos .fastq ordenados que no estén repetidos

for sampleraw in $fastqfiles1/*.fastq;do
   echo $sampleraw
   fastqc -t 12 $sampleraw -o ${directory1}/FastQCfiles/
done
echo "Fastqc ha finalizado"
echo " "

#TRIM-GALORE
echo "comienza el trim_galore"
for sample in $fastqfiles1/*.fastq; do

  	trim_galore $sample -o ${directory1}/TrimFiles
	echo trim-galore ha terminado para ${sample} 
done

echo "comienza el FASTQC de los archivos trimeados"
#FastQc a los archivos trimeados
for trimfile in ${directory1}/TrimFiles/*.fq;do

 fastqc -t 12 $trimfile -o ${directory1}/trimFastQCfiles
	echo "el fastqc ha terminado para $trimfile"
done

#llamamos al programa ./STARalingment.sh para que ejecute el alineamiento por STAR 
echo RNA-seq ANALYSIS has begun at $(date)
time ./STARalingment.sh $directory1 ${directory1}/trimFastQCfiles $genomeversion |& tee -a $logfile
echo STAR ANALYSIS has finished at $(date)

#Ejecutamos FeautureCaounts llamando al script quantSTAR
time ./quantSTAR_FeatureCounts.sh $directory1 $genomeversion $libtypeFeatureCounts $libPE |& tee -a $logfile
#directory=$1
#genomeversion=$2
#libtypeFeatureCounts=$3
#NofSample=$4
echo Quantification with Featurecounts has finished at $(date)


echo \n
#Análisis con el Multqc de todos los archivos que hemos obtenido del análisis
echo Multiqc Report

multiqc ${directory1}/FastQCfiles/ ${directory1}/TrimFiles ${directory1}/trimFastQCfiles ${directory1}/filesSTAR/BAMfiles/logSTAR/ ${directory1}/filesSTAR/FeatureCounts/ -o $directory1 

#Hemos mandado la salida vista en consola a un archivo de .log
echo "This script is sent to the log file, \"$logfile\"." 

