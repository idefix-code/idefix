#!/bin/bash
#SBATCH -A cin2231
#SBATCH --constraint=MI250
#SBATCH --job-name=idefix           # nom du job
#SBATCH --nodes=1                   # nombre de noeuds
#SBATCH --ntasks-per-node=8         # nombre de tache MPI par noeud  (=nb de GPU)
#SBATCH --gres=gpu:8                # nombre de GPU par noeud
#SBATCH --cpus-per-task=4           # nombre de coeurs CPU par tache (un quart du noeud ici)
#SBATCH --hint=nomultithread        # hyperthreading desactive
#SBATCH --time=00:20:00             # temps maximum d'execution demande (HH:MM:SS)
#SBATCH --output=idefix%j.out     # nom du fichier de sortie
#SBATCH --error=idefix-error%j.out      # nom du fichier d'erreur 

# nettoyage des modules charges en interactif et herites par defaut
module purge

# chargement des modules
module load cpe/25.09
module load craype-accel-amd-gfx90a craype-x86-trento
module load PrgEnv-amd
module load rocm/6.4.3
module load cray-python/3.11.7

export MPICH_GPU_SUPPORT_ENABLED=1

set -x

# execution du code
srun --ntasks-per-node=8 --cpus-per-task=4 --threads-per-core=1 --gpu-bind=closest -- ./test_reduce

