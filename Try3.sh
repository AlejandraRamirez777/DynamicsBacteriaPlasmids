#!/bin/bash

#SBATCH --job-name=Bac 		#Nombre del job
#SBATCH -p short			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem=2000		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=12:00:00			#Tiempo máximo de corrida, Default=2 horas
#SBATCH --mail-user=ma.ramirez@uniandes.edu.co
#SBATCH --mail-type=ALL			
#SBATCH -o TEST_job.o%j			#Nombre de archivo de salida

module load anaconda/python2
python Main_Code.py
module unload anaconda/python2

