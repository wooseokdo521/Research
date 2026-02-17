#!/bin/bash
set -e

#In terminal, chmod +x codon_partition_code.sh before executing ./

#Extract sequences (no name)
awk 'NR % 2 == 0' DanielTree_Nucleotide_TranslationAlignment_Standard_MAFFT_LINS_BLOSUM62_Gap1_Offset0.1.fasta > sequences.txt

#Extract 1st, 2nd, 3rd nucleotides of each codon
n=3 
for ((p=1; p<=n; p++)); do
    > "partition${p}.txt"
done
awk -v n=$n '{
    for(p=1; p<=n; p++) {
        part=""
        for(i=p; i<=length($0); i+=n) part = part substr($0,i,1)
        file = "partition" p ".txt"
        print part >> file
    }
}' sequences.txt

#Change txt to fasta
n=3
for p in $(seq 1 $n); do
    lines=$(wc -l < "partition${p}.txt")                    
    paste -d'\n' <(seq -f ">seq%g" $lines) "partition${p}.txt" > "partition${p}.fasta"
done

#Move fasta files and remove unnecessary files
mv partition1.fasta partition2.fasta partition3.fasta /home/xtremodevolab/Desktop/Daniel/modeltest/bin
rm sequences.txt partition1.txt partition2.txt partition3.txt

#Change directory to modeltest
cd modeltest/bin

#Run modeltest-ng
n=3  
mkdir -p results
cd results
for p in $(seq 1 $n); do
    ../modeltest-ng -i "../partition${p}.fasta" -d nt -p 4 -o "partitions${p}"
done

#Extract best model and create partition file
seq_length=$(awk '/^>/ {if(seqlen>0) {exit} next} {gsub(/[ \t\r\n]/,""); seqlen += length($0)} END{print seqlen}' /home/xtremodevolab/Desktop/Daniel/RAxML/DanielTree_Nucleotide_TranslationAlignment_Standard_MAFFT_LINS_BLOSUM62_Gap1_Offset0.1.fasta)

#Read models
i=1
for file in *.out
do
    model[$i]=$(
    grep -A2 "Best model according to AIC" "$file" | grep "Model:" | awk '{print $2; exit}')
    ((i++))
done

#Loop codonpositions 1 2 3
for positions in 1 2 3
do
    echo "${model[$positions]}, CodonPosition$positions = $positions-$seq_length/3"
done > DanielTree_nt_bestmodels_AIC.txt

#Move partition file
mv DanielTree_nt_bestmodels_AIC.txt /home/xtremodevolab/Desktop/Daniel/RAxML
cd ../../../RAxML

#Run RAxML
./raxml-ng --all --msa DanielTree_Nucleotide_TranslationAlignment_Standard_MAFFT_LINS_BLOSUM62_Gap1_Offset0.1.fasta --model DanielTree_nt_bestmodels_AIC.txt --tree rand{100},pars{100} --bs-metric fbp --bs-trees 1000 --workers 2 --threads 2 --prefix DanielTree 
