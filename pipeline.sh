bin=/data/groups/heewlee/MA/wazim/TimeSeries/FINAL/bin/ClonalTREE2
ref=/data/groups/heewlee/MA/wazim/TimeSeries/FINAL/util/ecoli.v3
data=/data/groups/heewlee/MA/wazim/TimeSeries/FINAL/data
out=/data/groups/heewlee/MA/wazim/TimeSeries/FINAL/out1
trimmed=/data/groups/heewlee/MA/wazim/TimeSeries/FINAL/trimmed
TRIMHOME=/data/groups/heewlee/MA/bin/Trimmomatic-0.33/
TRIMJAR=/data/groups/heewlee/MA/bin/Trimmomatic-0.33/trimmomatic-0.33.jar

declare -a arr1=("t1" "t2" "t3" "t4" "t5" "t6")

$1 = pop125

mkdir ${out}/${1}
mkdir ${trimmed}/${1}
for i in "${arr1[@]}"
do
    java -jar ${TRIMJAR} PE -threads 1 -phred33 -trimlog $trimmed/${1}/${1}_${i}.trimlog $data/${1}/${1}_${i}_1.fq $data/${1}/${1}_${i}_2.fq $trimmed/${1}/${1}_${i}_1.fq.gz $trimmed/${1}/${1}_${i}_1.fq.singleton.gz $trimmed/${1}/${1}_${i}_2.fq.gz $trimmed/${1}/${1}_${i}_2.fq.singleton.gz ILLUMINACLIP:${TRIMHOME}adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70 &> $trimmed/${1}/${1}_${i}.trimstat
    bwa mem -M -t 2 $ref/K12MG1655.fna $trimmed/${1}/${1}_${i}_1.fq.gz $trimmed/${1}/${1}_${i}_2.fq.gz > ${out}/${1}/${1}_PE
    samtools view -Sb -t $ref/K12MG1655_gList.txt ${out}/${1}/${1}_PE -o ${out}/${1}/${1}_tmp.bam
    samtools sort ${out}/${1}/${1}_tmp.bam ${out}/${1}/${1}_tmp.sorted
    samtools mpileup -f $ref/K12MG1655.fna ${out}/${1}/${1}_tmp.sorted.bam > ${out}/${1}/${1}_${i}.pileup
    rm -f ${out}/${1}/${1}_tmp.bam ${out}/${1}/${1}_tmp.sorted.bam ${out}/${1}/${1}*PE
    perl ${bin}/printBases.pl ${out}/${1}/${1}_${i}.pileup > ${out}/${1}/${1}_${i}.printBases
    python3 ${bin}/bialleleFreq.py ${out}/${1}/${1}_${i}.printBases ${i}
done
cat ${out}/${1}/${1}_*.freqs > ${out}/${1}/${1}.freqs
python3 ${bin}/filterFreqs.py ${out}/${1}/${1}.freqs ${out}/${1}/${1}

# python3 cluster.py 0.05 ${out}/${1}/${1}


