/*
 * =====================================================================================
 *
 *       Filename:  eutls.h
 *
 *    Description:  utlities header
 *
 *        Version:  1.0
 *        Created:  30/05/2019 16:34:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef _E_UTLS_H
#define _E_UTLS_H

#ifdef __cplusplus
extern "C" {
#endif
	//KMP algorithm
	//return occurence of a given pattern p in s	
	int e_kmp(char *s, char *p);
#ifdef __cplusplus
}
#endif
#endif



