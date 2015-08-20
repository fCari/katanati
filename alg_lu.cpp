/***********************************************************
	alg_lu.c -- LUÊ¬²ò
***********************************************************/
#include <stdlib.h>
#include <math.h>
#include "alg_matutil.h"  /* ¹ÔÎóÁàºî€ÎŸ®Æ»¶ñœž */

#include <string>

double alg_lu(int n, matrix a, int *ip)
{
	int i, j, k, ii, ik;
	double t, u, det;
	vector weight;
	weight = alg_new_vector(n);    /* weight[0..n-1] €Îµ­²±ÎÎ°è³ÎÊÝ */
	det = 0;                       /* ¹ÔÎóŒ° */
	for (k = 0; k < n; k++) {      /* ³Æ¹Ô€Ë€Ä€€€Æ */
		ip[k] = k;             /* ¹ÔžòŽ¹ŸðÊó€ÎœéŽüÃÍ */
		u = 0;                 /* €œ€Î¹Ô€ÎÀäÂÐÃÍºÇÂç€ÎÍ×ÁÇ€òµá€á€ë */
		for (j = 0; j < n; j++) {
			t = fabs(a[k][j]);  if (t > u) u = t;
		}
		if (u == 0) goto EXIT; /* 0 €Ê€é¹ÔÎó€ÏLUÊ¬²ò€Ç€­€Ê€€ */
		weight[k] = 1 / u;     /* ºÇÂçÀäÂÐÃÍ€ÎµÕ¿ô */
	}
	det = 1;                       /* ¹ÔÎóŒ°€ÎœéŽüÃÍ */
	for (k = 0; k < n; k++) {      /* ³Æ¹Ô€Ë€Ä€€€Æ */
		u = -1;
		for (i = k; i < n; i++) {  /* €è€ê²Œ€Î³Æ¹Ô€Ë€Ä€€€Æ */
			ii = ip[i]; /* œÅ€ß¡ßÀäÂÐÃÍ €¬ºÇÂç€Î¹Ô€òž«€Ä€±€ë */
			t = fabs(a[ii][k]) * weight[ii];
			if (t > u) {  u = t;  j = i;  }
		}
		ik = ip[j];
		if (j != k) {
			ip[j] = ip[k];  ip[k] = ik;  /* ¹ÔÈÖ¹æ€òžòŽ¹ */
			det = -det;  /* ¹Ô€òžòŽ¹€¹€ì€Ð¹ÔÎóŒ°€ÎÉä¹æ€¬ÊÑ€ï€ë */
		}
		u = a[ik][k];  det *= u;       /* ÂÐ³ÑÀ®Ê¬ */
		if (u == 0) goto EXIT;         /* 0 €Ê€é¹ÔÎó€ÏLUÊ¬²ò€Ç€­€Ê€€ */
		for (i = k + 1; i < n; i++) {  /* GaussŸÃµîË¡ */
			ii = ip[i];
			t = (a[ii][k] /= u);
			for (j = k + 1; j < n; j++)
				a[ii][j] -= t * a[ik][j];
		}
	}
EXIT:
	alg_free_vector(weight);   /* µ­²±ÎÎ°è€ò²òÊü */
	return det;                /* Ìá€êÃÍ€Ï¹ÔÎóŒ° */
}
void alg_solve(int n, matrix a, vector b, int *ip, vector x)
{
	int i, j, ii;
	double t;
	for (i = 0; i < n; i++) {       /* GaussŸÃµîË¡€Î»Ä€ê */
		ii = ip[i];  t = b[ii];
		for (j = 0; j < i; j++) t -= a[ii][j] * x[j];
		x[i] = t;
	}
	for (i = n - 1; i >= 0; i--) {  /* žåÂàÂåÆþ */
		t = x[i];  ii = ip[i];
		for (j = i + 1; j < n; j++) t -= a[ii][j] * x[j];
		x[i] = t / a[ii][i];
	}
}
double alg_matinv(int n, matrix a, matrix a_inv)
{
	int i, j, k, ii;
	double t, det;
	int *ip;   /* ¹ÔžòŽ¹€ÎŸðÊó */
	ip = (int *)malloc(sizeof(int) * n);
    
    if (ip == NULL) alg_error(std::string("alg_error").c_str());
	det = alg_lu(n, a, ip);
	if (det != 0)
		for (k = 0; k < n; k++) {
			for (i = 0; i < n; i++) {
				ii = ip[i];  t = (ii == k);
				for (j = 0; j < i; j++)
					t -= a[ii][j] * a_inv[j][k];
				a_inv[i][k] = t;
			}
			for (i = n - 1; i >= 0; i--) {
				t = a_inv[i][k];  ii = ip[i];
				for (j = i + 1; j < n; j++)
					t -= a[ii][j] * a_inv[j][k];
				a_inv[i][k] = t / a[ii][i];
			}
		}
	free(ip);
	return det;
}
double alg_gauss(int n, matrix a, vector b, vector x)
{
	double det;  /* ¹ÔÎóŒ° */
	int *ip;     /* ¹ÔžòŽ¹€ÎŸðÊó */
	ip = (int *)malloc(sizeof(int) * n);      /* µ­²±ÎÎ°è³ÎÊÝ */

	if (ip == NULL) alg_error("alg_error");
	det = alg_lu(n, a, ip);                   /* LUÊ¬²ò */
	if (det != 0) alg_solve(n, a, b, ip, x);  /* Ï¢Î©ÊýÄøŒ°€ò²ò€¯ */
	free(ip);                                 /* µ­²±ÎÎ°è€Î²òÊü */
	return det;                               /* Ìá€êÃÍ€Ï¹ÔÎóŒ° */
}
double alg_det(int n, matrix a)
{
	double det;  /* ¹ÔÎóŒ° */
	int *ip;     /* ¹ÔžòŽ¹€ÎŸðÊó */
	ip = (int *)malloc(sizeof(int) * n);       /* µ­²±ÎÎ°è³ÎÊÝ */
	if (ip == NULL) alg_error("alg_error");
	det = alg_lu(n, a, ip);                    /* LUÊ¬²ò */
	free(ip);                                  /* µ­²±ÎÎ°è€Î²òÊü */
	return det;                                /* ¹ÔÎóŒ° */
}
