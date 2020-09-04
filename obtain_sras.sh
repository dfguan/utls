if [ $# -lt 1 ]
then
	echo "Require an assembly summary as input"
	exit 1
fi

asmfl=$1
for projid in `cut -f15 $asmfl | tail -n +2 | sort | uniq`
do
	esearch -db sra -query $projid  | efetch -format runinfo | cut -d ',' -f 1 | grep [ES]RR  | xargs fastq-dump -n1 --split-files --gzip 
done

