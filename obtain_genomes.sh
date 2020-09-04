if [ ! -f "assembly_summary_genbank.txt" ]
then
	wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt 
fi
if [ ! -f "assembly_summary_refseq.txt" ]
then
	wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
fi

if [ $# -lt 1 ]
then
	echo "Require an assembly summary as input"
	exit 1
fi

asmfl=$1
for ga in `cut -f4 $asmfl | tail -n +2`
do
	gaftppath=`grep -w ^$ga assembly_summary_refseq.txt | cut -f19`	
	if [ $? -eq 0 ]
	then
		echo $gaftppath
		wget -c $gaftppath.fna.gz
	else
		gaftppath=`grep -w ^$ga assembly_summary_genbank.txt | cut -f19`	
		if [ $? -eq 0 ]
		then
			echo "downloading $ga from genbank"
			wget -c $gaftppath.fna.gz
		fi
	fi
done

