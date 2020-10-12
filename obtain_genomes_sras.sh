#!/bin/bash
# File              : obtain_genomes_sras.sh
# Author            : Dengfeng Guan <dfguan9@gmail.com>
# Date              : 23.09.2020
# Last Modified Date: 23.09.2020
# Last Modified By  : Dengfeng Guan <dfguan9@gmail.com>
set -x
USAGE="
`basename $0` [<OPTIONS>] <ARGUMENT>

ARGUMENT
 <input>                The list of assembly information delimited by tab         
 
OPTIONS
 -a <accession_ID>      The field number of accession id                   [default: 4] 
 -s <species_name>      The field number of species name                   [default: 5]
 -p <sample_id>         The field number of sample id                      [default: 16]
 -g <asmsum.txt>        Directory of assembly_summary_{genbank,refseq}.txt [default: .] 
 -f <step number>       Three bits to force to run the three steps         [default: 0]
 -S <step number>       Three bits to force to skip the steps              [default: 0]
 -O <output_directory>  Directory of output assemblies                     [default: .] 
 -r                     Remove directory for the species                   [default: FALSE]
 -v                     Verbose mode
"

afn=4
sfn=5
pfn=16
asmdir="."
header=1
fstp=0
sstp=0
rml=0
outdir="."

while getopts "a:s:p:f:S:g:O:rh" OPT "$@"; do
    case $OPT in
        a) afn="$OPTARG"
			;;
		s) sfn="$OPTARG"
			;;
		p) pfn="$OPTARG"
			;;
		f) fstp="$OPTARG"
			;;
		S) sstp="$OPTARG"
			;;
		r) rml=1
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
mkdir -p $outdir
outd=`readlink -f $outdir`

stp1=$((fstp & 0x1))
stp2=$((fstp & 0x2))
stp3=$((fstp & 0x4))

sstp1=$((sstp & 0x1))
sstp2=$((sstp & 0x2))
sstp3=$((sstp & 0x4))
if [ $stp1 -ne 0 ] && [ $sstp1 -ne 0 ]
then
	echo "Not allow to force and skip step 1 at the same time, exiting"
	exit 1
elif [ $stp2 -ne 0 ] && [ $sstp2 -ne 0 ]
then
	echo "Not allow to force and skip step 2 at the same time, exiting"
	exit 1
elif [ $stp3 -ne 0 ] && [ $sstp3 -ne 0 ]
then
	echo "Not allow to force and skip step 3 at the same time, exiting"
	exit 1
fi

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
	stp1_ind=0
	downl=0
	if [ $sstp1 -eq 0 ]
	then
		if [ -f "$outd"/."$spn".asm.done ] && [ $stp1 -eq 0 ]
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
					downl=1
				else
					echo "Failed to download assembly for $spn"
					stp1_ind=1
				fi
			fi	
		fi
	else
		echo "User chooses to skip step 1 of downloading the assembly"
	fi

	stp2_ind=0	
	if [ $sstp2 -eq 0 ]
	then
		if [ -f "$outd"/.$spn.sra.done ] && [ $stp2 -eq 0 ]
		then
			echo "already downloaded sras for $spn, skipped"
		else
			echo "Start step 2, downloading SRAs for $spn"
			outputd=$outd/$spn/SRAs
			mkdir -p $outputd
			# be careful with esearch who is reading from stdin
			esearch -db sra -query $samid < /dev/null | efetch -format runinfo | grep WGS | grep GENOMIC | grep ILLUMINA | cut -d ',' -f1 | grep [ES]RR  > $outputd/sralist  
			if [ -s $outputd/sralist ]
			then
				prefetch -C yes -X 1000000000 -O $outputd  --option-file $outputd/sralist > $outputd/prefetch.log.o 2>$outputd/prefetch.log.e
				if [ $? -eq 0 ]
				then
					touch "$outd"/."$spn".sra.done
					downl=1
				else
					echo "Failed to download SRAs for $spn" 
					stp2_ind=1
				fi	
			else
				echo "$outputd/sralist is empty, will do nothing"	
				touch "$outd"/."$spn".sra.done
			fi	
		fi
	else
		echo "User chooses to skip step 2 of downloading the sras"
	fi
	# unless user choose to skip step 3 otherwise if some data is downloaded, we need to upload it to ther server 	
	if [ $sstp3 -eq 0 ]
	then
		#if [ -f "$outd"/.$spn.upload.done ] && [ $stp3 -eq 0 ]
		#then
			#echo "Already uploaded data for $spn, skip step 3"
		#else
		if [ $downl -eq 1 ]
		then
			echo "Start step 3, uploading data for $spn to Huawei Cloud OBS"
			obsutil cp -vmd5 -u -r -f $outd/$spn obs://nextomics-customer/WHWLZ-201906006A/genomic_diversity 	
			if [ $? -eq 0 ] && [ $rml -eq 1 ] && [ ! -z $spn ] && [ $stp1_ind -eq 0 ] && [ $stp2_ind -eq 0 ]
			#then
				#touch $outd/.$spn.upload.done
			then
				rm -rf $outd/$spn
			fi
			#fi
		fi
		#fi
	else
		echo "User chooses to skip step 3 of uploading the datasets"
	fi
done < $outd/aid_spn_sid_list


