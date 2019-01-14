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
#include <stdarg.h>
#include <stdint.h>
#include <getopt.h>

#include "color.h"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define PKG "seq2sam"

typedef struct {
	char *s;
	size_t l,n;
}str_t;

typedef struct {
	str_t qn, tn, ntn, seq, qual, cigar;
	int flag, pos, maq, tl, npos; 
}sam_t;

int prtinfo(int type, FILE *stream, const char *format, ...) //type: 3:error 2:warning 1:message otherwise ignore
{
	va_list args;
	va_start (args, format);
	if (type == 1) 
		fprintf(stream, COLOR_YELLOW);	
	 else if (type == 2) 
		fprintf(stream, COLOR_BLUE);	
	else if (type == 3)
		fprintf(stream, COLOR_RED);
	vfprintf (stream, format, args);
	if (type > 0) fprintf(stream, COLOR_RESET);	
	va_end (args);
	return 0;
}

int trim_slide(char *s) 
{
	if (s) {
		int sl = strlen(s); 
		if (sl > 1) {
			if (s[sl-2] == '/' && (s[sl-1] == '2' || s[sl-1] == '1')) 
				s[sl-2] = 0;
		}
	}
	return 0; //or maybe return sl - 2/ sl?
}




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

//from bwa_escape 
int seq2sam_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return 0;
}



int dump_hdr(FILE *fout, char *rg_line)
{
	if (rg_line && rg_line[0] == '@') {
		seq2sam_escape(rg_line);	
		fprintf(fout, "@HD\tVN:1.6\n%s\n",rg_line);
	} else 
		return 1;
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

int single2sam(char *fn, char *rg_line, char *out) // 
{
	long n = 0, slen = 0, qlen = 0;

	FILE *fout = out && *out ? fopen(out, "w") : stdout;	
	if (dump_hdr(fout, rg_line)) {
		prtinfo(3, stderr,"[E::%s] read group header line error\n", __func__);
		return 1;
	}
		
	int sam_num_lim = 50000;	
	sam_t *sts = (sam_t *)calloc(sam_num_lim, sizeof(sam_t));
	int sam_num = 0;
	
	gzFile fp;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *seq = kseq_init(fp);
	
	while (1) {
		if (kseq_read(seq) >= 0)  {
			trim_slide(seq->name.s);
			if (!write_sam(&sts[sam_num], seq, 4)) {
				sam_num += 1;
				if (sam_num >= sam_num_lim) {
					n += sam_num;
					dump_sam_core(fout, sts, sam_num);
					sam_num = 0;
				}
			}
		} else {
			n += sam_num;
			if (sam_num) dump_sam_core(fout, sts, sam_num);
			break;	
		}
	}
	prtinfo(1, stderr, "[M::%s] process %ld reads\n", __func__, n);
	kseq_destroy(seq);
	gzclose(fp);
	release_sam(sts, sam_num_lim);
	return 0;			
}

//paired end reads to sam format
int pair2sam(char *fn1, char *fn2, char *rg_line, char *out)
{
	long n = 0;

	FILE *fout = out && *out ? fopen(out, "w") : stdout;	
	if (dump_hdr(fout, rg_line)) {
		prtinfo(3, stderr,"[E::%s] read group header line error\n", __func__);
		return 1;
	}
	
	int sam_num_lim = 50000;	
	sam_t *sts = (sam_t *)calloc(sam_num_lim, sizeof(sam_t));
	int sam_num = 0;
	
	gzFile fp1, fp2;
	fp1 = fn1 && strcmp(fn1, "-")? gzopen(fn1, "r") : gzdopen(fileno(stdin), "r");
	fp2 = fn2 && strcmp(fn2, "-")? gzopen(fn2, "r") : gzdopen(fileno(stdin), "r");
	kseq_t *seq1 = kseq_init(fp1);
	kseq_t *seq2 = kseq_init(fp2);	
	int is_bad_file = 0;
	while (1) {
		int rtn1 = kseq_read(seq1);
	   	int rtn2 = kseq_read(seq2);
		if (rtn1 >= 0 && rtn2 >=0) {
			//compare read name
			trim_slide(seq1->name.s);
			trim_slide(seq2->name.s);
			/*fprintf(stderr, "%s", seq1->name.s);*/
			/*fprintf(stderr, "%s", seq2->name.s);*/
			if (!strcmp(seq1->name.s, seq2->name.s)) {
				if (!write_sam(&sts[sam_num], seq1, 77) && !write_sam(&sts[sam_num+1], seq2,141))  {
					sam_num += 2;
					if (sam_num >= sam_num_lim) {
						n += sam_num;
						dump_sam_core(fout, sts, sam_num);
						sam_num = 0;
					}
				}
			} else {
				is_bad_file = 1;	
				break;
			} 
		} else {
			if (sam_num) dump_sam_core(fout, sts, sam_num);
			n += sam_num;
			if (rtn1 * rtn2 < 0)  //a better way to check both numbers are less than zero?
				is_bad_file = 1;
				break;
		} 		
		
	}
	if (is_bad_file) prtinfo(2, stderr,"[W::%s] sequence file is corrupted\n", __func__);
	prtinfo(1, stderr, "[M::%s] process %ld read pairs\n", __func__, n >> 1);
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
				if (c != 'h') prtinfo(3, stderr, "[E::%s] undefined option %c\n", __func__, c);
	help:	
				fprintf(stderr, "\nUsage: %s [options] -1 <FASTA/Q>\n", PKG);
				fprintf(stderr, "Options:\n");	
				fprintf(stderr, "         -o    STR      output file [stdout]\n");	
				fprintf(stderr, "         -1    STR      read1 file\n");
				fprintf(stderr, "         -2    STR      read2 file\n");	
				fprintf(stderr, "         -r    STR      read group header line such as \'@RG\\tID:foo\\tSM:bar\'\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	
	if (fn1 && fn2) {
		prtinfo(1, stderr, "[M::%s] Program starts\n", __func__);	
		pair2sam(fn1, fn2, rg_line, out);
		prtinfo(1, stderr, "[M::%s] Program ends\n", __func__);	
		return 0;	
	} else if (fn1) {
		prtinfo(1, stderr, "[M::%s] Program starts\n", __func__);	
		single2sam(fn1, rg_line, out);
		prtinfo(1, stderr, "[M::%s] Program ends\n", __func__);	
		return 0;	
	} else {
		prtinfo(3, stderr, "[E::%s] require at least one sequencing file\n", __func__);
		goto help;
		return 1;
	}	

}


