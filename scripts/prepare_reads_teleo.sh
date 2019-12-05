################################################### Prepare reads for downstream analysis #######################################################

#it needs a three column file: sample name - path to R1 -  path to R2
#raw read files are in .gz
names=(`cut -f 1 bioman_teleo_samples.txt`)
R1=(`cut -f 3 bioman_teleo_samples.txt`)
R2=(`cut -f 4 bioman_teleo_samples.txt`)

echo -n "" > fasta_file.txt
echo -n "" > group.txt
echo -n "" > bioman_teleo.fasta
echo "Sample	raw	cut	assembled	retained" > reads_stats.txt

for i in `seq 0 51`; do

	echo -n ${names[i]} >> reads_stats.txt
	echo -n "	"   >> reads_stats.txt
	gzip -dc < /share/raw/MiSeq/2017/171215_M03233_0078_000000000-BK563/${R1[i]} > ${names[i]}_R1.fastq
	gzip -dc < /share/raw/MiSeq/2017/171215_M03233_0078_000000000-BK563/${R2[i]} > ${names[i]}_R2.fastq
	r="$(grep -c '^@M0' ${names[i]}_R1.fastq)"
    echo -n $r >> reads_stats.txt
	echo "" > ${names[i]}.log

#Remove forward and reverse primers from read1 and read2

	cutadapt -g ACACCGCCCGTCACTCT...CATGGTAAGTGTACCGGAAG -G CTTCCGGTACACTTACCATG...AGAGTGACGGGCGGTGT --discard-untrimmed --minimum-length 30 -e 0.2 -o ${names[i]}_R1_cut.fastq -p ${names[i]}_R2_cut.fastq ${names[i]}_R1.fastq ${names[i]}_R2.fastq >> ${names[i]}.log
	echo -n "	"   >> reads_stats.txt
	r="$(grep -c '^@M0' ${names[i]}_R1_cut.fastq)"
    echo -n $r >> reads_stats.txt

#merge pairs

	pear -f ${names[i]}_R1_cut.fastq -r ${names[i]}_R2_cut.fastq -v 20 -m 100 -n 40 -o ${names[i]}  >> ${names[i]}.log
    echo -n "	"   >> reads_stats.txt
    r="$(grep -c '^@M0' ${names[i]}.assembled.fastq)"
    echo -n $r >> reads_stats.txt

#remove low quality reads

	trimmomatic SE -phred33  ${names[i]}.assembled.fastq ${names[i]}.assembled_retained.fastq AVGQUAL:25 >> ${names[i]}.log
    echo -n "	"   >> reads_stats.txt
    r="$(grep -c '^@M0' ${names[i]}.assembled_retained.fastq)"
    echo  $r >> reads_stats.txt

#convert into fasta 

	mothur "#fastq.info(fastq=${names[i]}.assembled_retained.fastq)"

	if (($i > 50)) ; then
		echo ${names[i]}.assembled_retained.fasta >> fasta_file.txt
		echo ${names[i]} >> group.txt
	else
		echo -n ${names[i]}.assembled_retained.fasta"-" >> fasta_file.txt
		echo -n ${names[i]}"-" >> group.txt
	fi

	cat ${names[i]}.assembled_retained.fasta >> bioman_teleo.fasta

done

file="fasta_file.txt"
fastas=$(cat "$file")
file="group.txt"
groups=$(cat "$file")

mothur "#make.group(fasta=$fastas, groups=$groups)"
mv mergegroups bioman_teleo.groups
rm -f mothur*logfile
