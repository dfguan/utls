#!/bin/bash
# Author: dfguan9
# Contact: dfguan9@gmail.com
# Insitute: Zoology Institute, CAS
# Date: 2020-09-05
# Function: obtain sras from a given assembly list containing project id

USAGE="
`basename $0` [<OPTIONS>] <ARGUMENT>

ARGUMENT
 <input>                The list of assembly information delimited by tab         
 
OPTIONS
 -p <project_id>        The field number of project id                     [default: 15] 
 -s <species_name>      The field number of species name                   [default: 5]
 -g <asmsum.txt>        Directory of assembly_summary_{genbank,refseq}.txt [default: .] 
 -v                     Verbose mode
"
sfn=5
pfn=15
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
        #h) header=0
			#;;
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
cut -f$pfn,$sfn -d$'\t'  $asmfl | grep ^PRJ | sort | uniq | sed 's/ /_/g' > $outd/aid_pid_list
#fi
which esearch > /dev/null
if [ $? -ne 0 ]
then
	echo "esearch is not found, please follow this guide to install https://www.ncbi.nlm.nih.gov/books/NBK179288/"
	exit 1
fi
which fastq-dump > /dev/null
if [ $? -ne 0 ]
then
	echo "fastq-dump is not found, please follow this guide to install https://ncbi.github.io/sra-tools/install_config.html"
	exit 1
fi

while read -r projid dirn 
do 
	outputd=$outd/$dirn/SRAs
	mkdir -p $outputd 
	esearch -db sra -query $projid  | efetch -format runinfo | cut -d ',' -f1 | grep [ES]RR  | xargs -n1 fastq-dump  -O $outputd --split-files --gzip 
done < $outd/aid_pid_list


#for projid in `cut -f15 $asmfl | tail -n +2 | sort | uniq`
#do
	#mkdir -p
#done

