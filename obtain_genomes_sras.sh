#!/bin/bash
# File              : obtain_genomes_sras.sh
# Author            : Dengfeng Guan <dfguan9@gmail.com>
# Date              : 23.09.2020
# Last Modified Date: 23.09.2020
# Last Modified By  : Dengfeng Guan <dfguan9@gmail.com>
#set -x
USAGE="
`basename $0` [<OPTIONS>] <ARGUMENT>

ARGUMENT
 <input>                The list of assembly information delimited by tab         
 
OPTIONS
 -a <accession_ID>      The field number of accession id                   [default: 4] 
 -s <species_name>      The field number of species name                   [default: 5]
 -p <sample_id>         The field number of sample id                      [default: 5]
 -g <asmsum.txt>        Directory of assembly_summary_{genbank,refseq}.txt [default: .] 
 -O <output_directory>  Directory of output assemblies [default: .] 
 -v                     Verbose mode
"

afn=4
sfn=5
pfn=16
asmdir="."
header=1
rml=0
outdir="."

while getopts "a:O:s:g:p:fh" OPT "$@"; do
    case $OPT in
        a) afn="$OPTARG"
			;;
		s) sfn="$OPTARG"
			;;
		p) pfn="$OPTARG"
			;;
		f) rml=1
			;;
        g) asmdir="$OPTARG"
			;;
        h) header=0
			;;
        O) outdir="$OPTARG"
			;;
        \?) echo "Invalid option: -$OPTARG" >&2 
            exit 1 
			;;
        :) echo "Option -$OPTARG requires an argument." >&2
           exit 1
        ;;
    esac
done

shift $((OPTIND-1))
[[ $# -eq 1 ]] || { printf "$USAGE" >&2 && exit 1; }

asmfl=$1
outd=`readlink -f $outdir`
if [ ! -f $asmdir/"assembly_summary_genbank.txt" ]
then
	echo "assembly_summary_genbank.txt is not found under "$asmdir", now downloading......" 
	wget -N -P  $asmdir https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt 
	if [ $? -ne 0 ]
	then
		echo "Can not retrieve genbank summary, please check your network connection" 		
		exit 1
	fi
fi
if [ ! -f $asmdir/"assembly_summary_refseq.txt" ]
then
	echo "assembly_summary_genbank.txt is not found under "$asmdir", now downloading......" 
	wget -N -P $asmdir https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
	if [ $? -ne 0 ]
	then
		echo "Can not retrieve refseq summary, please check your network connection" 		
		exit 1
	fi
fi

#if [ $header -eq 1 ]
#then
	#cut -f$afn,$sfn -d$'\t'  $asmfl | head -n -1 | sed 's/ /_/g' > $outd/aid_spn_list
#else
cut -f$afn,$sfn,$pfn -d$'\t'  $asmfl | grep ^GC | sed 's/ /_/g' > $outd/aid_spn_sid_list
#fi

while read -r ga spn samid
do 
	#echo $ga $spn $samid
	if [ -f "$outd"/."$spn".asm.done ]
	then
		echo "already downloaded genome for $spn, skipped"
	else
		echo "Start step 1, downloading the assembly for $spn"
		gaftppath=`grep -w ^$ga "$asmdir"/assembly_summary_refseq.txt | cut -d$'\t' -f20`	
		# use -z to check if the string is empty 
		[ -z "$gaftppath" ] && gaftppath=`grep -w ^$ga "$asmdir"/assembly_summary_genbank.txt | cut -d$'\t' -f20`	
		if [ -z "$gaftppath" ]
		then
			echo "Warning: we can not find accession $ga, it may have been deleted on NCBI server"		
		else
			mkdir -p $outd/$spn
			gan=`basename $gaftppath`
			echo Now start downloading $gaftppath
			wget -N -c -o $outd/$spn/"$spn".wget.log -P  $outd/$spn "$gaftppath"/"$gan"_genomic.fna.gz 
			if [ $? -eq 0 ]
			then 
				touch $outd/."$spn".asm.done
			else
				echo "Failed to download assembly for $spn"
			fi
		fi	
	fi
	if [ -f "$outd"/.$spn.sra.done ]
	then
		echo "already downloaded sras for $spn, skipped"
	else
		echo "Start step 2, downloading SRAs for $spn"
		outputd=$outd/$spn/SRAs
		mkdir -p $outputd
		# be careful with esearch who is reading from stdin
		esearch -db sra -query $samid < /dev/null | efetch -format runinfo | grep WGS | grep GENOMIC | cut -d ',' -f1 | grep [ES]RR  > $outputd/sralist  
		prefetch -C yes -X 1000000000 -O $outputd  --option-file $outputd/sralist > $outputd/prefetch.log.o 2>$outputd/prefetch.log.e
		if [ $? -eq 0 ]
		then
			touch "$outd"/."$spn".sra.done
		else
			echo "Failed to download SRAs for $spn" 
		fi		
	fi
	if [ -f "$outd"/.$spn.upload.done ]
	then
		echo "Already uploaded data for $spn, skip step 3"
	else
		echo "Start step 3, uploading data for $spn to Huawei Cloud OBS"
		obsutil cp -vmd5 -u -r -f $outd/$spn obs://nextomics-customer/WHWLZ-201906006A/genomic_diversity 	
		if [ $? -eq 0 ] 
		then
			touch $outd/.$spn.upload.done
			if [ $rml -eq 1 ] && [ ! -z $spn ] 
			then
				rm -rf $outd/$spn
			fi
		fi
	fi
done < $outd/aid_spn_sid_list


