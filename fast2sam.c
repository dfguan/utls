/*
 * =====================================================================================
 *
 *       Filename:  fastq2sam.c
 *
 *    Description:  convert fastq to sam format
 *
 *        Version:  1.0
 *        Created:  12/01/2019 10:04:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	char *s;
	size_t l,n;
}str_t;

typedef struct {
	str_t qn, tn, ntn, seq, qual, cigar;
	int flag, pos, maq, tl, npos; 
}sam_t;

int cp_str(str_t *t, kstring_t *ks)
{
	if (ks->l) {
		if (t->n <= ks->l) {
			if (t->l) free(t->s);
			t->s = (char *)malloc(sizeof(char)*(ks->l+1));
			t->n = ks->l + 1;
		}
		strcpy(t->s, ks->s);
		t->l = ks->l;
		return 0;
	} else 
		return 1;
}

int rls_str(str_t *t)
{
	if (t) {
		if (t->s) free(t->s);	
	}
	return 0;

}


int release_sam(sam_t *sts, int n)
{
	int i;
	if (!sts) return 0;
	for ( i = 0 ; i < n; ++i) {
		rls_str(&sts[i].cigar); rls_str(&sts[i].ntn);
		rls_str(&sts[i].qn); rls_str(&sts[i].qual);
		rls_str(&sts[i].seq); rls_str(&sts[i].tn);
	}
	free(sts);
	return 0;
}

int write_sam(sam_t *st, kseq_t *seq, int flag)
{
	if (!cp_str(&st->qn, &seq->name)) {
		st->flag = flag;
		if (!cp_str(&st->seq, &seq->seq)) {
			cp_str(&st->qual, &seq->qual);
			return 0;
		} else 
			return 1;
	} else 
		return 1;
}


int dump_hdr(FILE *fout, char *rg_line)
{
	fprintf(fout, "@HD\tVN:1.6\n%s\n",rg_line);
	return 0;	
}

int dump_sam_core(FILE *fout, sam_t *sts, int num)
{
	int i;
	for ( i = 0; i < num; ++i) {
		char *qn = sts[i].qn.s, *tn = sts[i].tn.s, *cigar=sts[i].cigar.s,*ntn = sts[i].ntn.s, *seq = sts[i].seq.s, *qual = sts[i].qual.s;
		int flag = sts[i].flag, pos = sts[i].pos, maq = sts[i].maq, npos = sts[i].npos, tl = sts[i].tl; 
		fprintf(fout, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n", qn, flag, tn?tn:"*", pos, maq, cigar? cigar:"*", ntn?ntn:"*", npos, tl, seq, qual?qual:"*");		
	}	
	return 0;
}

int se_fast2sam(char *fn, char *rg_line, char *out) // 
{
	long n = 0, slen = 0, qlen = 0;

	FILE *fout = out && *out ? fopen(out, "w") : stdout;	
	dump_hdr(fout, rg_line);
		
	int sam_num_lim = 50000;	
	sam_t *sts = (sam_t *)calloc(sam_num_lim, sizeof(sam_t));
	int sam_num = 0;
	
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *seq = kseq_init(fp);
	
	while (kseq_read(seq) >= 0) {
		write_sam(&sts[sam_num], seq, 4);
		if (++sam_num > sam_num_lim) {
			dump_sam_core(fout, sts, sam_num_lim);
			sam_num = 0;
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	release_sam(sts, sam_num_lim);
	return 0;			
}


//paired end reads to sam format
int pe_fast2sam(char *fn1, char *fn2, char *rg_line, char *out)
{
	long n = 0;

	FILE *fout = out && *out ? fopen(out, "w") : stdout;	
	dump_hdr(fout, rg_line);
		
	int sam_num_lim = 50000;	
	sam_t *sts = (sam_t *)calloc(sam_num_lim, sizeof(sam_t));
	int sam_num = 0;
	
	gzFile fp1, fp2;
	fp1 = fn1 && strcmp(fn1, "-")? gzopen(fn1, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *seq1 = kseq_init(fp1);
	kseq_t *seq2 = kseq_init(fp2);	
	int is_bad_file = 0;
	while (1) {
		int rtn1 = kseq_read(seq1);
	   	int rtn2 = kseq_read(seq2);
		if (rtn1 >= 0 && rtn2 >=0) {
			//compare read name
			if (!strcmp(seq1->name.s, seq2->name.s)) {
				if (!write_sam(&sts[sam_num], seq1, 77) && !write_sam(&sts[sam_num+1], seq2,141)) 
					sam_num += 2;
			} else {
				is_bad_file = 1;	
				break;
			} 
		} else if (rtn1 * rtn2 < 0) { //a better way to check both numbers are less than zero?
			is_bad_file = 1;
			break;
		} else 
			break;
		if (++sam_num > sam_num_lim) {
			dump_sam_core(fout, sts, sam_num_lim);
			sam_num = 0;
		}
	}
	kseq_destroy(seq1);
	kseq_destroy(seq2);
	gzclose(fp1);
	gzclose(fp2);
	release_sam(sts, sam_num_lim);
	return 0;			


}


int main(int argc, char *argv[])
{
	char *rg_line = 0;
	char *out = 0;
	char *fn1 = 0, *fn2 = 0;	
	int c;
	while (~(c=getopt(argc, argv, "p:1:2:r:h"))) {
		switch (c) {
			case 'p': 
				out = optarg;
				break;
			case '1':
				fn1 = optarg;
				break;
			case '2':
				fn2 = optarg; 
				break;
			case 'r':
				rg_line = optarg; 
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
	/*help:	*/
				fprintf(stderr, "\nUsage: fast2sam [options] ...\n");
				fprintf(stderr, "Options:\n");	
				fprintf(stderr, "         -o    STR      output file [stdout]\n");	
				fprintf(stderr, "         -1    STR      read1 file\n");
				fprintf(stderr, "         -2    STR      read2 file\n");	
				fprintf(stderr, "         -r    STR      read group header line such as \'@RG\tID:foo\tSM:bar\'\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	
	if (rg_line) { fprintf(stderr, "read group can not be ignored, process exit"); return 1;}

	if (fn1 && fn2) {
		fprintf(stderr, "Program starts\n");	
		pe_fast2sam(fn1, fn2, rg_line, out);
		fprintf(stderr, "Program ends\n");	
		return 0;	
	} else if (fn1) {
		fprintf(stderr, "Program starts\n");	
		se_fast2sam(fn1, rg_line, out);
		fprintf(stderr, "Program ends\n");	
		return 0;	
	} else {
		fprintf(stderr, "require at least one sequencing file\n");
		return 1;
	}	

}


