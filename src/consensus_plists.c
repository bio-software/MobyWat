// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "consensus_plists.h"

re *consensus_plists(re *p1, int no1, re *p2, int no2, int *cons_no, double p_tol) {

	int i, j, label = 0, label_2 = 0;
	int cons_final_no = 0;
	re *cons = NULL;
	re *cons_cleaned = NULL;	
	re *x = NULL;
	
	for(j = 0; j < no1; j++)
		(p1+j)->label_mind_dist = 0;		// unused set to 0

	label = 0;
	for(i = 0; i < no2; i++) {
		label_2 = 0; 		// to have only one consensus result
		for(j = 0; j < no1; j++) 
			if ( label_2 == 0 && (p1+j)->label_mind_dist == 0 && ( ( ((p1+j)->x_coord - (p2+i)->x_coord) * ((p1+j)->x_coord - (p2+i)->x_coord) + ((p1+j)->y_coord - (p2+i)->y_coord) * ((p1+j)->y_coord - (p2+i)->y_coord) + ((p1+j)->z_coord - (p2+i)->z_coord) * ((p1+j)->z_coord - (p2+i)->z_coord) ) < CONS_SCALE * CONS_SCALE * p_tol * p_tol ) ) {
				if ( label == 0 ) {
					cons = (re *) malloc(sizeof(re));
					*(cons+label) = *(p2+i);
					(cons+label)->avg_water_no = (double) ( (p1+j)->water_no + (p2+i)->water_no ) / 2;
					(cons+label)->x_coord = ( (p1+j)->x_coord + (p2+i)->x_coord ) / 2;
					(cons+label)->y_coord = ( (p1+j)->y_coord + (p2+i)->y_coord ) / 2;
					(cons+label)->z_coord = ( (p1+j)->z_coord + (p2+i)->z_coord ) / 2;
					label++;
					label_2++;
					(p1+j)->label_mind_dist = 1;
					}
				else {
					label++;
					x = realloc(cons, label*sizeof(re));
					if (x) cons = x;
					x = NULL;
					*(cons+label-1) = *(p2+i);
					(cons+label-1)->avg_water_no = (double) ( (p1+j)->water_no + (p2+i)->water_no ) / 2;
					(cons+label-1)->x_coord = ( (p1+j)->x_coord + (p2+i)->x_coord ) / 2;
					(cons+label-1)->y_coord = ( (p1+j)->y_coord + (p2+i)->y_coord ) / 2;
					(cons+label-1)->z_coord = ( (p1+j)->z_coord + (p2+i)->z_coord ) / 2;
					label_2++;
					(p1+j)->label_mind_dist = 1;
					}
				}
		}
			
	insertion_sort_re (cons, label, 'a');	

	cons_cleaned = (re *) malloc(label*sizeof(re));
	cons_cleaned = clean_plist(cons, label, p_tol, &cons_final_no);	

	free(cons);

	*cons_no = cons_final_no;

	return(cons_cleaned);
	}

