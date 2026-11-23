// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "global.h"

re *merge_plists(re *p1, int no1, re *p2, int no2, int *merged_no, double p_tol) {

	int merged_final_no = 0;
	int no_sum = 0;
	re *merged = NULL;
	re *merged_cleaned = NULL;

	void insertion_sort_re (re *a, int num, char swics);
	re *clean_plist(re *mean, int rows, double identical_tolerance, int *plist_no);

	no_sum = no1 + no2 ;

	merged = (re *) malloc(no_sum*sizeof(re));
	memcpy(merged, p1, no1*sizeof(re));
	memcpy(merged+no1, p2, no2*sizeof(re));

	insertion_sort_re (merged, no_sum, 'x');	

	merged_cleaned = (re *) malloc(no_sum*sizeof(re));
	merged_cleaned = clean_plist(merged, no_sum, p_tol, &merged_final_no);	

	free(merged);

	*merged_no = merged_final_no;

	return(merged_cleaned);
	}

