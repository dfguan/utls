#!/bin/bash
# File              : obtain_sras.sh
# Author            : Dengfeng Guan <dfguan9@gmail.com>
# Date              : 22.09.2020
# Last Modified Date: 23.09.2020
# Last Modified By  : Dengfeng Guan <dfguan9@gmail.com>
# Description: obtain sras from a given assembly list containing project id
#set -x
USAGE="
`basename $0` [<OPTIONS>] <ARGUMENT>

ARGUMENT
 <input>                The list of assembly information delimited by tab         
 
OPTIONS
 -p <sample_id>        The field number of project id                     [default: 16] 
 -s <species_name>      The field number of species name                   [default: 5]
 -g <asmsum.txt>        Directory of assembly_summary_{genbank,refseq}.txt [default: .] 
 -O <output_directory>  Directory of output assemblies [default: .] 
 -v                     Verbose mode
"
sfn=5
pfn=16
asmdir="."
#header=1
outd="."
while getopts "p:O:s:g:" OPT "$@"; do
    case $OPT in
        p) pfn="$OPTARG"
			;;
		s) sfn="$OPTARG"
			;;
        g) asmdir="$OPTARG"
			;;
        O) outd="$OPTARG"
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

#if [ $header -eq 1 ]
#then
	#cut -f$pfn,$sfn -d$'\t'  $asmfl | head -n -1 | sort | uniq | sed 's/ /_/g' > $outd/aid_pid_list
#else
cut -f$pfn,$sfn -d$'\t'  $asmfl | grep SAMN | sort | uniq | sed 's/ /_/g' > $outd/spn_sid_list
#fi
which esearch > /dev/null
if [ $? -ne 0 ]
then
	echo "esearch is not found, please follow this guide to install https://www.ncbi.nlm.nih.gov/books/NBK179288/"
	exit 1
fi
which prefetch > /dev/null
if [ $? -ne 0 ]
then
	echo "fastq-dump is not found, please follow this guide to install https://ncbi.github.io/sra-tools/install_config.html"
	exit 1
fi

while read -r dirn samid
do 
	outputd=$outd/$dirn/SRAs
	mkdir -p $outputd 
	esearch -db sra -query $samid  | efetch -format runinfo | cut -d ',' -f1 | grep [ES]RR  > $outputd/sralist 
	prefetch -C yes -X 1000000000 -O $outputd  --option-file $outputd/sralist > $outputd/prefetch.log.o 2>$outputd/prefetch.log.e
done < $outd/spn_sid_list


#for projid in `cut -f15 $asmfl | tail -n +2 | sort | uniq`
#do
	#mkdir -p
#done

