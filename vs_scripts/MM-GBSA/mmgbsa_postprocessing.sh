#!/bin/bash

mkdir ../mmgbsa_input
mkdir ../mmgbsa_maegz
mkdir ../mmgbsa_csv
mkdir ../mmgbsa_error

for i in $(seq 0 15)  # Change the range according to the total number of ligands
do 
    if test -s "lig_${i}-out.maegz"
    then
        rm -rf "lig_${i}_tmp"
        mv "lig_${i}-out.maegz" ../mmgbsa_maegz
        mv "lig_${i}-out.csv" ../mmgbsa_csv
    else
        echo "lig_${i} has no results!"
        mv "lig_${i}_tmp" ../mmgbsa_error
    fi
    mv "lig_${i}.maegz" ../mmgbsa_input
done

cd ..
mv lig mmgbsa_log
mv mmgbsa_input lig
cd mmgbsa_log

rm mmgbsa.sbatch
rm mmgbsa_postprocessing.sh
