// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "global.h"

void insertion_sort_int_re_ss (int *a, re **b, int num) {
	
	int lb, ub, i, j;
	int temp_3;
	re *temp_2;
	
	
	lb = 0;
	ub = num - 1;

	for (i = lb + 1; i <= ub; i++) {
		temp_3 = *(a+i);
		temp_2 = *(b+i);
		for (j = i-1; j >= lb && *(a+j) < temp_3; j--) {
			*(a+j+1) = *(a+j);
			*(b+j+1) = *(b+j);
			}
		*(a+j+1) = temp_3;
		*(b+j+1) = temp_2;
		}

	}


void insertion_sort_int_ss (int **a, int num) {
	
	int lb, ub, i, j;
	int temp_1;
	int temp_2;
	int temp_3;	
	
	lb = 0;
	ub = num - 1;

	for (i = lb + 1; i <= ub; i++) {
		temp_3 = *(*(a+2)+i);
		temp_2 = *(*(a+1)+i);
		temp_1 = *(*(a+0)+i);
		for (j = i-1; j >= lb && *(*(a+2)+j) < temp_3; j--) {
			*(*(a+2)+j+1) = *(*(a+2)+j);
			*(*(a+1)+j+1) = *(*(a+1)+j);
			*(*(a+0)+j+1) = *(*(a+0)+j);
			}
		*(*(a+2)+j+1) = temp_3;
		*(*(a+1)+j+1) = temp_2;
		*(*(a+0)+j+1) = temp_1;
		}

	}


void insertion_sort_wr (wr *a, int num) {
	
	int lb, ub, i, j;
	int temp_int;
	wr temp_wr;
	
	lb = 0;
	ub = num - 1;

	for (i = lb + 1; i <= ub; i++) {
		temp_int = (a+i)->counter;
		temp_wr = *(a+i);	
		for (j = i-1; j >= lb && (a+j)->counter < temp_int; j--) {
			(a+j+1)->counter = (a+j)->counter;
			*(a+j+1) = *(a+j);
			}
		(a+j+1)->counter = temp_int;
		*(a+j+1) = temp_wr;
		}
	}


void insertion_sort_re (re *a, int num, char swics) {
	
	int lb, ub, i, j;
	int temp_int;
	re temp_re;
	
	lb = 0;
	ub = num - 1;

	if ( swics != 'a' ) {
		for (i = lb + 1; i <= ub; i++) {
			temp_int = (a+i)->water_no;
			temp_re = *(a+i);	
			for (j = i-1; j >= lb && (a+j)->water_no < temp_int; j--) {
				(a+j+1)->water_no = (a+j)->water_no;
				*(a+j+1) = *(a+j);
				}
			(a+j+1)->water_no = temp_int;
			*(a+j+1) = temp_re;
			}
		}

	else {
		for (i = lb + 1; i <= ub; i++) {
			temp_int = (a+i)->avg_water_no;
			temp_re = *(a+i);	
			for (j = i-1; j >= lb && (a+j)->avg_water_no < temp_int; j--) {
				(a+j+1)->avg_water_no = (a+j)->avg_water_no;
				*(a+j+1) = *(a+j);
				}
			(a+j+1)->avg_water_no = temp_int;
			*(a+j+1) = temp_re;
			}
		}

	}

void insertion_sort_res_ser (ap *a, int num) {
	
	int lb, ub, i, j;
	int temp_int;
	ap temp_re;
	
	lb = 0;
	ub = num - 1;


	for (i = lb + 1; i <= ub; i++) {
		temp_int = (a+i)->res_ser;
		temp_re = *(a+i);	
		for (j = i-1; j >= lb && (a+j)->res_ser > temp_int; j--) {
			(a+j+1)->res_ser = (a+j)->res_ser;
			*(a+j+1) = *(a+j);
			}
		(a+j+1)->res_ser = temp_int;
		*(a+j+1) = temp_re;
		}

	}

void insertion_sort_hs_wat_res_ser (hs *a, int num) {
	
	int lb, ub, i, j;
	int temp_int;
	hs temp_re;
	
	lb = 0;
	ub = num - 1;

	for (i = lb + 1; i <= ub; i++) {
		temp_int = (a+i)->wat_res_ser;
		temp_re = *(a+i);	
		for (j = i-1; j >= lb && (a+j)->wat_res_ser > temp_int; j--) {
			(a+j+1)->wat_res_ser = (a+j)->wat_res_ser;
			*(a+j+1) = *(a+j);
			}
		(a+j+1)->wat_res_ser = temp_int;
		*(a+j+1) = temp_re;
		}
	}

void insertion_sort_hs_dist (hs *a, int num) {
	
	int lb, ub, i, j;
	double temp_double;
	hs temp_re;
	
	lb = 0;
	ub = num - 1;

	for (i = lb + 1; i <= ub; i++) {
		temp_double = (a+i)->dist;
		temp_re = *(a+i);	
		for (j = i-1; j >= lb && (a+j)->dist > temp_double; j--) {
			(a+j+1)->dist = (a+j)->dist;
			*(a+j+1) = *(a+j);
			}
		(a+j+1)->dist = temp_double;
		*(a+j+1) = temp_re;
		}
	}


void insertion_sort_hd_res_ser (hd *a, int num) {
	
	int lb, ub, i, j;
	int temp_int;
	hd temp_re;
	
	lb = 0;
	ub = num - 1;

	for (i = lb + 1; i <= ub; i++) {
		temp_int = (a+i)->res_ser;
		temp_re = *(a+i);	
		for (j = i-1; j >= lb && (a+j)->res_ser > temp_int; j--) {
			(a+j+1)->res_ser = (a+j)->res_ser;
			*(a+j+1) = *(a+j);
			}
		(a+j+1)->res_ser = temp_int;
		*(a+j+1) = temp_re;
		}
	}

void insertion_sort_hd_type (hd *a, int num) {
	
	int lb, ub, i, j;
	char temp_char;
	hd temp_re;
	
	lb = 0;
	ub = num - 1;

	for (i = lb + 1; i <= ub; i++) {
		temp_char = (a+i)->type;
		temp_re = *(a+i);	
		for (j = i-1; j >= lb && (a+j)->type > temp_char; j--) {
			(a+j+1)->type = (a+j)->type;
			*(a+j+1) = *(a+j);
			}
		(a+j+1)->type = temp_char;
		*(a+j+1) = temp_re;
		}
	}

