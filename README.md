# Utilities
- seq2sam: convert single-end or paired-end reads into sam format
- 10x\_trim: trim 10x data
- eutls: some algorithms, now has KMP alg

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



