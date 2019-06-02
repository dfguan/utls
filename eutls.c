/*
 * =====================================================================================
 *
 *       Filename:  eutls.c
 *
 *    Description:  realization of functions defined in the header file
 *
 *        Version:  1.0
 *        Created:  30/05/2019 16:39:40
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
#include <string.h>

#include "eutls.h"
//compute partial match table
int ecal_pmt(char *p, int *t, int lp)
{
	int i = 1;
	int c = 0;
	t[0] = -1;
	
	for (i = 1; i < lp; ++i, ++c) {
		if (p[i] == p[c]) 
			t[i] = t[c];
		else {
			t[i] = c;
			c = t[c];
			while (c >= 0 && p[i] != p[c]) 
			   c = t[c];	
		}
	}
	t[i] = c;
	return 0;
}


size_t ekmp(char *s, char *p)
{
	//construct Tables	
	size_t tmp_lp = strlen(p); //potential bug: p is too large
	if (tmp_lp > 0X7FFFFFFF) {
		fprintf(stderr, "[W::%s] pattern size is beyond 2G, will not calculate\n", __func__);
		return 0;
	}
	int lp = tmp_lp;
	int *t = (int *)malloc(sizeof(char) * (lp + 1));
	ecal_pmt(p, t, lp);


	int	k; 
	size_t j, occ; 
	
	for ( j = 0, k = 0, occ = 0; j < strlen(s);) {
		if (p[k] == s[j]) {
			++j, ++k;
			if (k == lp) ++occ, k = t[k];
		} else {
			k = t[k];
			if (k < 0) ++j, ++k;
		}
	}
	free(t);
	return occ;
}

#ifdef TEST

int main()
{
	fprintf(stderr, "KMP test: \n");
	char *s = "ABABDABACDABABCABAB";
	char *p = "ABAB";
	fprintf(stderr, "s: %s, p: %s, occ: %lu\n", s, p, ekmp(s, p));
	return 0;
}

#endif

