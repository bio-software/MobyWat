// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "global.h"
extern ts *g_tolerances;

mr **match_clusters(re **mean, int *clusters, int rows, re *ref_waters, int num_ref_waters) {

	int i, j, k;
	double dist;

	mr **match = NULL;

	double dist_calc (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch);

	// Allocate memory for match
	match = (mr **) malloc(rows * sizeof(mr *));
	for (j = 0; j < rows; j++)
		*(match+j) = (mr *) malloc(*(clusters+j) * sizeof (mr));

	for(i = 0; i < rows; i++)
		for(j = 0; j < *(clusters+i); j++) {
			(*(match+i)+j)->ref_water_ser_no = 0;
			(*(match+i)+j)->distance = 0;
			for(k = 0; k < num_ref_waters; k++) {
				dist = dist_calc( (*(mean+i)+j)->x_coord, (*(mean+i)+j)->y_coord, (*(mean+i)+j)->z_coord, (ref_waters+k)->x_coord, (ref_waters+k)->y_coord, (ref_waters+k)->z_coord, "O", "O", 'q' );
				if (dist < (*(match+i)+j)->distance || k == 0) {
					(*(match+i)+j)->distance = dist;
					(*(match+i)+j)->ref_water_ser_no = k+1;
					(*(match+i)+j)->ref_water_atom_ser_no = (ref_waters+k)->target_atom_ser_no;
					(*(match+i)+j)->ref_water_res_ser_no = (ref_waters+k)->target_res_ser;
					}	
			}
		}

	return (match);
	}

mr *match_plist(re *predict, int rows, re *ref_waters, int num_ref_waters, double b_fact_threshold, double s_value) {

	int i, k, selected_k = 0, counter = 0;
	double dist;

	mr *match = NULL;

	double dist_calc (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch);

	// Allocate memory for match
	match = (mr *) malloc(rows * sizeof(mr));

	for(i = 0; i < rows; i++) {
                (match+i)->indicator = 0;       // 0 if no match found for the prediction
		(match+i)->ref_water_ser_no = 0;
		(match+i)->distance = 0;
		(match+i)->ref_water_atom_ser_no = 0;
		(match+i)->ref_water_res_ser_no = 0;
		(match+i)->b_factor = 0;				
		}


	for(i = 0; i < num_ref_waters; i++)
                (ref_waters+i)->label_mind_dist = 0;			// 0 to mark usable ref waters (all refs are used only once for match)

	for(i = 0; i < rows; i++) {
		selected_k = 0;
		counter = 0;
		for(k = 0; k < num_ref_waters; k++)
			if ( (ref_waters+k)->label_mind_dist == 0 ) {
				dist = dist_calc( (predict+i)->x_coord, (predict+i)->y_coord, (predict+i)->z_coord, (ref_waters+k)->x_coord, (ref_waters+k)->y_coord, (ref_waters+k)->z_coord, "O", "O", 'q' );
				if ( (dist < (match+i)->distance || counter == 0) ) {
					(match+i)->distance = dist;
					(match+i)->ref_water_ser_no = k+1;
					(match+i)->ref_water_atom_ser_no = (ref_waters+k)->target_atom_ser_no;
					(match+i)->ref_water_res_ser_no = (ref_waters+k)->target_res_ser;
					(match+i)->b_factor = (ref_waters+k)->b_factor;				
					selected_k = k;
					counter++;
					}
				}
		if ( (ref_waters+selected_k)->b_factor < b_fact_threshold && (match+i)->distance < s_value && (match+i)->ref_water_ser_no != 0) {
		        (match+i)->indicator = selected_k+1;
			(ref_waters+selected_k)->label_mind_dist = 1;	// one match is allowed for all ref waters
			}  
	  
	      }

	return (match);
	}

