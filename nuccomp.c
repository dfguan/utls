/*
 * =====================================================================================
 *
 *       Filename:  nuccomp.c
 *
 *    Description: calculate nucletide compositions 
 *
 *        Version:  1.0
 *        Created:  07/01/2021 14:38:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"

KSEQ_INIT2(, gzFile, gzread)

int comp(char *fa)
{	
	gzFile fp;
	kseq_t *seq;
	char toint[256] = {0};
	toint['a'] = 1, toint['c'] = 2, toint['g'] = 3, toint['t'] = 4, toint['n'] = 5;
	toint['A'] = 6, toint['C'] = 7, toint['G'] = 8, toint['T'] = 9, toint['N'] = 10;
	char toDNA[] = "acgtnACGTNO";
	/*long cnt[11] = {0};*/
	int n = 0, slen = 0, qlen = 0;
	fp = gzopen(fa, "r");
	seq = kseq_init(fp);

	int i;

	fprintf(stdout, "chr");
	for (i = 0; i < 11; ++i) 
		fprintf(stdout, "\t#%c", toDNA[i]);
	fprintf(stdout, "\n");
	while (kseq_read(seq) >= 0) {
		long tmpcnt[11] = {0};
		for (i = 0; i < seq->seq.l; ++i) {
			int z = toint[seq->seq.s[i]];
			if (z) ++tmpcnt[z-1];
			else ++tmpcnt[9];
		}
		++n, slen += seq->seq.l;
		fprintf(stdout, "%s", seq->name.s);
		for (i = 0; i < 11; ++i) 
			fprintf(stdout, "\t%ld", tmpcnt[i]);
		fprintf(stdout, "\n");
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;

}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "nuccomp <fasta/fasta.gz>\n");
		return 1;
	}
	comp(argv[1]);	
	return 0;
}
