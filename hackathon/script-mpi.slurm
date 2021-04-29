#!/bin/bash
#SBATCH --job-name=idefix_Gaylor     # nom du job
#SBATCH --nodes=1                   # nombre de noeuds
#SBATCH --ntasks-per-node=4         # nombre de tache MPI par noeud  (=nb de GPU)
#SBATCH --gres=gpu:4                # nombre de GPU par noeud
#SBATCH --cpus-per-task=10          # nombre de coeurs CPU par tache (un quart du noeud ici)
# /!\ Attention, "multithread" fait reference à l'hyperthreading dans la terminologie Slurm
#SBATCH --hint=nomultithread        # hyperthreading desactive
#SBATCH --time=00:10:00             # temps maximum d'execution demande (HH:MM:SS)
#SBATCH --output=idefix%j.out     # nom du fichier de sortie
#SBATCH --error=idefix-error%j.out      # nom du fichier d'erreur (ici commun avec la sortie)
##SBATCH --signal=SIGUSR2@100      # send USR2 signal 5 minutes before the end of the job

# nettoyage des modules charges en interactif et herites par defaut
module purge

# chargement des modules

module load cuda/10.2
module load gcc/8.3.1
module load openmpi/4.0.5-cuda

# echo des commandes lancees
set -x
#export KOKKOS_PROFILE_LIBRARY=~/src/kokkos-tools/kp_space_time_stack.so

# execution du code
srun ./idefix -dec 2 1 2  --kokkos-num-devices=4
