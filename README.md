# utilities
- seq2sam.c: convert single-end or paired-end reads from fastq/a format into bam format, users can define sample id and read id. 

# Installation
```
make 
```

# Description
```
Usage: seq2sam [options] -1 <FASTA/Q>
Options:
         -o    STR      output file [stdout]
         -1    STR      read1 file
         -2    STR      read2 file
         -r    STR      read group header line such as '@RG\tID:foo\tSM:bar'
         -h             help
```



