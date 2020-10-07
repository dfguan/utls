#!/bin/bash
# Author: dfguan9
# Date: 2020-08-14
# Function: extract an assembly summary with taxonomy id or species name

USAGE="
`basename $0` [<OPTIONS>] <ARGUMENT>

ARGUMENT
 <output>               output file        
 
OPTIONS
 -t <taxids>            obtain the assembly summary for the taxid given.
 -s <species_name>      obtain the assembly summary for the species given.
 -n <number_of_records> obtain n records [default: -1 (all records)].
 -v                     Verbose mode
"
nrec=-1
while getopts "t:s:n:" OPT "$@"; do
    case $OPT in
        t) taxid="$OPTARG"
		   value=$taxid
		   usetax=1
		   exturl="/assembly_descriptors/taxid/${taxid}"
			;;
		s) species="$OPTARG"
		   value="$species"
		   usetax=0
		   exturl="/assembly_descriptors/organism/${species}"
			;;
		n) nrec="$OPTARG"
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
out=$1

#baseurl="https://api.ncbi.nlm.nih.gov/datasets/v1alpha"
#sufurl="?returned_content=COMPLETE"
#url=${baseurl}${exturl}${sufurl}

#[[ -f $value.json ]] || curl -s -X GET $url -H  "accept: application/json" -o $value.json

if [ ! -f "rankedlineage.dmp" ]
then
	if [ ! -f "new_taxdump.tar.gz" ]
	then
		wget -c -N https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
	fi
	tar zvxf  new_taxdump.tar.gz rankedlineage.dmp
fi

tag=`echo $value | sed 's/ /_/g'` 
if [ -f .$tag.asminfo.done ]
then
	echo "already pull assembly information, skipped" >&2
else
	echo "pull assembly information" >&2
	if [ $usetax -eq 1 ]
	then
		python3 get_assembly_summary.py -t $value -n $nrec > asminfo.tmp 
	else 
		python3 get_assembly_summary.py -s "$value" -n $nrec > asminfo.tmp
	fi
	if [ $? -eq 0 ]
	then
		echo "assembly infor obtained successfully" >&2
		touch .$tag.asminfo.done
	else
		echo "An error occured, please contact Dengfeng" >&2
		exit  1
	fi
fi
#cut -f8 $value.asm.tsv | xargs -n1 -I{} grep -w ^{} rankedlineage.dmp | awk -F'|' '{print $4"\t"$5"\t"$6}' | awk '{print $1"\t"$2"\t"$3}' > $value.lin.tsv
#paste $value.lin.tsv $value.asm.tsv > $out

if [ -f .$tag.srainfo.done ]
then
	echo "already obtained sra information, skipped" >&2
	exit 0
else
	echo "pull sra information" >&2
	cut -f16 -d$'\t' asminfo.tmp | grep ^SAM > sample_id.tmp
	for samid in `cat sample_id.tmp`
	do
		echo "processing $samid" >&2
		if [ $samid == "NA" ]
		then
			echo -e "0\t0" 
		else
			esearch -db sra -query $samid < /dev/null | efetch -format runinfo | grep WGS | grep GENOMIC | grep ILLUMINA | awk -F, -v sam_id=$samid 'BEGIN{tb=0;cnt=0}{tb+=$5;cnt+=1}END{print cnt"\t"tb}'
			if [ $? -ne 0 ]
			then
				echo "Fail to get sra information for $samid, now exit" >&2
				exit 1
			fi
		fi
	done > srainfo.tmp
	touch .$tag.srainfo.done
	if [ -f .$tag.collection.done ]
	then
		echo "already collected assembly and sra inforation, skipped" >&2
		exit 0
	else
		paste -d$'\t' asminfo.tmp srainfo.tmp > $out
		if [ $? -eq 0 ]
		then
			echo "Finished collecting assembly and sra information successfully" >&2
			touch .$tag.collection.done 
			rm -f asminfo.tmp srainfo.tmp sample_id.tmp
			exit 0
		fi
	fi
fi
