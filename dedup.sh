#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=interactive               #REQUIRED: which partition to use
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB

/usr/bin/time -v ./lindley_deduper.py -f sorted_C1_SE_uniqAlign.sam -u /projects/bgmp/klindley/bioinfo/Bi624/dedup/Deduper/STL96.txt -o outfile.txt 