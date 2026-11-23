// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 10-04-2016

#include "global.h"
extern double g_m_tolerance;

re *select_if_waters (ap *pdb_target, ap *pdb_ligand, ap *pdb_waters, int pdb_num_at_target, int pdb_num_at_ligand, int pdb_num_at_waters, int *num_at_common, double dist_tol) {
		
	int dist_num_wat_tar, dist_num_wat_lig, num_common_wat = 0;

	pd *distances_wat_tar = NULL;
	pd *distances_wat_lig = NULL;
	re *common_wat_res = NULL;

	pd *pdb_ligand_dist_full (ap *pdb_target, ap *pdb_ligand, int pdb_num_at_target, int pdb_num_at_ligand, int *dist_num, double dist_tolerance);
	re *extract_common_wat_res (pd *distances_t,  pd *distances_l, int num_dist_t, int num_dist_l, int *num_common_wat, int frame);

	distances_wat_tar = pdb_ligand_dist_full (pdb_target, pdb_waters, pdb_num_at_target, pdb_num_at_waters, &dist_num_wat_tar, dist_tol);
	distances_wat_lig = pdb_ligand_dist_full (pdb_ligand, pdb_waters, pdb_num_at_ligand, pdb_num_at_waters, &dist_num_wat_lig, dist_tol);
	common_wat_res = extract_common_wat_res (distances_wat_tar, distances_wat_lig, dist_num_wat_tar, dist_num_wat_lig, &num_common_wat, 0);

	free(distances_wat_tar);
	free(distances_wat_lig);

	*num_at_common = num_common_wat;

	return (common_wat_res);
	}

re *select_sf_waters (ap *pdb_target, ap *pdb_waters, int pdb_num_at_target, int pdb_num_at_waters, int *num_sf, double dist_tol) {
		
	int dist_num_wat_tar;

	pd *distances_wat_tar = NULL;
	re *surface_wat_res = NULL;

	pd *pdb_ligand_dist_full (ap *pdb_target, ap *pdb_ligand, int pdb_num_at_target, int pdb_num_at_ligand, int *dist_num, double dist_tolerance);
	re *copy_pd_to_re (pd *distances_t, int num_dist_t, int frame);

	distances_wat_tar = pdb_ligand_dist_full (pdb_target, pdb_waters, pdb_num_at_target, pdb_num_at_waters, &dist_num_wat_tar, dist_tol);
	surface_wat_res = copy_pd_to_re (distances_wat_tar, dist_num_wat_tar, 0);

	free(distances_wat_tar);
	*num_sf = dist_num_wat_tar;

	return (surface_wat_res);
	}


// The following function selects ligand residues (used for waters as ligands here) close to a target and produces a distance list
pd *pdb_ligand_dist_full (ap *pdb_target, ap *pdb_ligand, int pdb_num_at_target, int pdb_num_at_ligand, int *dist_num, double dist_tolerance) {	

	int j, k;
	int counter;
	double temp_dist;
	double dist_tolerance_sqr;
	int *label = NULL;

	pd *distances = NULL;
	pd *x = NULL;

	label = (int *) malloc(pdb_num_at_ligand * sizeof(int));
	for (j = 0; j < pdb_num_at_ligand; j++) *(label+j) = 0;		// 0 if not yet found

	dist_tolerance_sqr = dist_tolerance * dist_tolerance;

	counter = 0;

	for (j = 0; j < pdb_num_at_ligand; j++) 
		if ( *(label+j) == 0 && (pdb_ligand+j)->pdb_type[0] != 'H' )
			for (k = 0; k < pdb_num_at_target; k++) 
				if ( (pdb_target+k)->pdb_type[0] != 'H' ) {
					temp_dist = ((pdb_ligand+j)->x_coord - (pdb_target+k)->x_coord)*((pdb_ligand+j)->x_coord - (pdb_target+k)->x_coord) + ((pdb_ligand+j)->y_coord - (pdb_target+k)->y_coord)*((pdb_ligand+j)->y_coord - (pdb_target+k)->y_coord) + ((pdb_ligand+j)->z_coord - (pdb_target+k)->z_coord)*((pdb_ligand+j)->z_coord - (pdb_target+k)->z_coord);
					if (temp_dist < dist_tolerance_sqr) {
						if ( counter == 0) distances = (pd *) malloc(sizeof(pd));
							else {
								x = realloc( distances, (counter+1) * sizeof(pd) );
								if (x) distances = x;
								x = NULL;
								}
						(distances+counter)->ligand_atom_ser_no = (pdb_ligand+j)->atom_ser;
						(distances+counter)->ligand_res_ser = (pdb_ligand+j)->res_ser;
						(distances+counter)->ligand_x_coord = (pdb_ligand+j)->x_coord;
						(distances+counter)->ligand_y_coord = (pdb_ligand+j)->y_coord;
						(distances+counter)->ligand_z_coord = (pdb_ligand+j)->z_coord;
						(distances+counter)->ligand_b_factor = (pdb_ligand+j)->b_factor;
//						(distances+counter)->target_dist = temp_dist;
						counter++;
						*(label+j) = 1;
						break;
						}	
					}
			
	free(label);
	*(int *)dist_num = counter;
	return distances;
	}

// From 2 distance lists this function selects waters common in both lists (IF waters)
re *extract_common_wat_res (pd *distances_t,  pd *distances_l, int num_dist_t, int num_dist_l, int *num_common_wat, int frame) {	

	int i, j;
	int counter;

	re *residues = NULL;
	re *x = NULL;

	counter = 0;
	for (i = 0; i < num_dist_l; i++)
		for (j = 0; j < num_dist_t; j++)
			if ( (distances_l+i)->ligand_res_ser == (distances_t+j)->ligand_res_ser && (distances_l+i)->ligand_atom_ser_no == (distances_t+j)->ligand_atom_ser_no ) {
				if ( counter == 0) residues = (re *) malloc(sizeof(re));
					else {
						x = realloc( residues, (counter+1) * sizeof(re) );
						if (x) residues = x;
						x = NULL;
						}
				(residues+counter)->ligand_ser_no = frame;		// It is useful to mark in B-factor column when printing out OW-s after clustering
				(residues+counter)->target_atom_ser_no = (distances_l+i)->ligand_atom_ser_no;
				(residues+counter)->target_res_ser = (distances_l+i)->ligand_res_ser;
				(residues+counter)->x_coord = (distances_l+i)->ligand_x_coord;
				(residues+counter)->y_coord = (distances_l+i)->ligand_y_coord;
				(residues+counter)->z_coord = (distances_l+i)->ligand_z_coord;
				(residues+counter)->b_factor = (distances_l+i)->ligand_b_factor;
				counter++;
				}

	*(int *)num_common_wat = counter;
	return residues;
	}

// This selects waters from one distance list (SF waters)
re *copy_pd_to_re (pd *distances_t, int num_dist_t, int frame) {	

	int counter = 0;
	re *residues = NULL;
	re *x = NULL;

	for (counter = 0; counter < num_dist_t; counter++) {
		if ( counter == 0) residues = (re *) malloc(sizeof(re));
			else {
				x = realloc( residues, (counter+1) * sizeof(re) );
				if (x) residues = x;
				x = NULL;
				}
		(residues+counter)->ligand_ser_no = frame;		// It is useful to mark in B-factor column when printing out OW-s after clustering
		(residues+counter)->target_atom_ser_no = (distances_t+counter)->ligand_atom_ser_no;
		(residues+counter)->target_res_ser = (distances_t+counter)->ligand_res_ser;
		(residues+counter)->x_coord = (distances_t+counter)->ligand_x_coord;
		(residues+counter)->y_coord = (distances_t+counter)->ligand_y_coord;
		(residues+counter)->z_coord = (distances_t+counter)->ligand_z_coord;
		(residues+counter)->b_factor = (distances_t+counter)->ligand_b_factor;
		}

	return residues;
	}

pl *copy_re_to_pl (re *residues, int num_residues) {	

	int counter = 0;
	pl *pool = NULL;
	pl *x = NULL;

	for (counter = 0; counter < num_residues; counter++) {
		if ( counter == 0) pool = (pl *) malloc(sizeof(pl));
			else {
				x = realloc( pool, (counter+1) * sizeof(pl) );
				if (x) pool = x;
				x = NULL;
				}
		(pool+counter)->aser = (residues+counter)->target_atom_ser_no;
		(pool+counter)->rser = (residues+counter)->target_res_ser;
		(pool+counter)->x_coord = (residues+counter)->x_coord;
		(pool+counter)->y_coord = (residues+counter)->y_coord;
		(pool+counter)->z_coord = (residues+counter)->z_coord;
		}

	return pool;
	}

pd *pdb_ligand_dist_full_min (ap *pdb_target, ap *pdb_ligand, int pdb_num_at_target, int pdb_num_at_ligand, int *dist_num, double dist_tolerance, char h_switch) {	

	int j, k;
	int counter;
	double temp_dist;
	pd *distances = NULL;

	double dist_calc_min (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch);

	counter = 0;

	for (j = 0; j < pdb_num_at_ligand; j++) 
		for (k = 0; k < pdb_num_at_target; k++) 
			if ( dist_calc_min((pdb_ligand+j)->x_coord, (pdb_ligand+j)->y_coord, (pdb_ligand+j)->z_coord, (pdb_target+k)->x_coord, (pdb_target+k)->y_coord, (pdb_target+k)->z_coord, (pdb_ligand+j)->pdb_type, (pdb_target+k)->pdb_type, h_switch) < dist_tolerance ) 
				counter++;

	*(int *)dist_num = counter;
	distances = (pd *) malloc(counter*sizeof(pd));

	counter = 0;

	for (j = 0; j < pdb_num_at_ligand; j++) 
		for (k = 0; k < pdb_num_at_target; k++) {
			temp_dist = dist_calc_min((pdb_ligand+j)->x_coord, (pdb_ligand+j)->y_coord, (pdb_ligand+j)->z_coord, (pdb_target+k)->x_coord, (pdb_target+k)->y_coord, (pdb_target+k)->z_coord, (pdb_ligand+j)->pdb_type, (pdb_target+k)->pdb_type, h_switch);
			if (temp_dist < dist_tolerance) {
				(distances+counter)->ligand_ser_no = (pdb_ligand+j)->model_ser;
				(distances+counter)->ligand_atom_ser_no = (pdb_ligand+j)->atom_ser;
				strcpy((distances+counter)->ligand_pdb_type,(pdb_ligand+j)->pdb_type);							
				strcpy((distances+counter)->ligand_res_type,(pdb_ligand+j)->res_type);							
				strcpy((distances+counter)->ligand_chain,(pdb_ligand+j)->chain);							
				(distances+counter)->ligand_res_ser = (pdb_ligand+j)->res_ser;
				(distances+counter)->ligand_b_factor = (pdb_ligand+j)->b_factor;
				(distances+counter)->target_atom_ser_no = (pdb_target+k)->atom_ser;
				strcpy((distances+counter)->target_pdb_type,(pdb_target+k)->pdb_type);							
				strcpy((distances+counter)->target_res_type,(pdb_target+k)->res_type);							
				strcpy((distances+counter)->target_chain,(pdb_target+k)->chain);							
				(distances+counter)->target_res_ser = (pdb_target+k)->res_ser;
				(distances+counter)->target_b_factor = (pdb_target+k)->b_factor;
				(distances+counter)->target_dist = temp_dist;
				counter++;
				}	
			}
	return distances;
	}

double dist_calc_min (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch) {

	double distance;	
	
	if ((pdb_type1[0]  != 'H' && pdb_type2[0]  != 'H') || h_switch == 'h') {
		distance = sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));
		return distance;
		}
		else return VERY_BIG_DISTANCE;	
	}

re *extract_common_wat_res_min (pd *distances_t,  pd *distances_l, int num_dist_t, int num_dist_l, int *num_common_wat, char min_dist_switch, double min_dist_tolerance) {	

	int i, j;
	int counter;
	int num_ligand_residues;
	int num_target_residues;

	re *residues = NULL;
	re *ligand_residues = NULL;
	re *target_residues = NULL;

	re *extract_lig_res (pd *distances, int num_dist, int *num_dist_res);

	min_dist_switch = 'x';

	ligand_residues = extract_lig_res(distances_l, num_dist_l, &num_ligand_residues);

	for (i = 0; i < num_ligand_residues; i++)
		(ligand_residues+i)->label_mind_dist = 0;

	// Marking waters with close contact to the ligand

	if (min_dist_switch == 'x') {			
		for (i = 0; i < num_ligand_residues; i++)
			for (j = 0; j < num_dist_l; j++)
				if ((distances_l+j)->ligand_res_ser == (ligand_residues+i)->target_res_ser && abs((distances_l+j)->ligand_atom_ser_no - (ligand_residues+i)->target_atom_ser_no) < 3 && (distances_l+j)->target_dist < min_dist_tolerance)
					(ligand_residues+i)->label_mind_dist = 1;
		}

	target_residues = extract_lig_res(distances_t, num_dist_t, &num_target_residues);

	for (i = 0; i < num_target_residues; i++)
		(target_residues+i)->label_mind_dist = 0;

// Marking waters with close contact to the protein
//	if (min_dist_switch == 'x') {
//		for (i = 0; i < num_target_residues; i++)
//			for (j = 0; j < num_dist_t; j++)
//				if ((distances_t+j)->ligand_res_ser == (target_residues+i)->target_res_ser && abs((distances_t+j)->ligand_atom_ser_no - (target_residues+i)->target_atom_ser_no) < 3 && (distances_t+j)->target_dist < min_dist_tolerance)
//					(target_residues+i)->label_mind_dist = 1;
//		}


//		for (i = 0; i < num_target_residues; i++)
//			printf("%d%d\n",(target_residues+i)->target_res_ser, (target_residues+i)->label_mind_dist);


	counter = 0;
	for (i = 0; i < num_ligand_residues; i++)
		for (j = 0; j < num_target_residues; j++)
			if ((ligand_residues+i)->label_mind_dist == 0 && (target_residues+j)->label_mind_dist == 0  && (ligand_residues+i)->target_res_ser == (target_residues+j)->target_res_ser  && abs((ligand_residues+i)->target_atom_ser_no - (target_residues+j)->target_atom_ser_no ) < 3 ) {
				counter++;
				}

	*(int *)num_common_wat = counter;
	residues = (re *) malloc(counter*sizeof(re));

	counter = 0;
	for (i = 0; i < num_ligand_residues; i++)
		for (j = 0; j < num_target_residues; j++)
			if ((ligand_residues+i)->label_mind_dist == 0 && (target_residues+j)->label_mind_dist == 0 && (ligand_residues+i)->target_res_ser == (target_residues+j)->target_res_ser && abs((ligand_residues+i)->target_atom_ser_no - (target_residues+j)->target_atom_ser_no) < 3 ) {
				(residues+counter)->ligand_ser_no = (ligand_residues+i)->ligand_ser_no;
				(residues+counter)->target_atom_ser_no = (ligand_residues+i)->target_atom_ser_no;
				strcpy((residues+counter)->target_res_type, (ligand_residues+i)->target_res_type);
				strcpy((residues+counter)->target_chain, (ligand_residues+i)->target_chain);
				(residues+counter)->target_res_ser = (ligand_residues+i)->target_res_ser;
				counter++;
				}

	free(ligand_residues);
	free(target_residues);
		
	return residues;
	}

re *extract_lig_res (pd *distances, int num_dist, int *num_dist_res) {	

	int i, j;
	int counter;
	re *residues = NULL;

	int *label;
	label = (int *) malloc(num_dist*sizeof(int));

	counter = 0;
	for (i = 0; i < num_dist; i++) *(label+i) = 0;	

	for (i = 0; i < num_dist; i++) 
		if (*(label+i) == 0) {
			counter++;
			for (j = 0; j < num_dist; j++) 
				if (strcmp((distances+i)->ligand_res_type,(distances+j)->ligand_res_type) == 0 && strcmp((distances+i)->ligand_chain,(distances+j)->ligand_chain) == 0 && (distances+i)->ligand_res_ser == (distances+j)->ligand_res_ser)
					*(label+j) = 1;
			}

	*(int *)num_dist_res = counter;
	residues = (re *) malloc(counter*sizeof(re));

	counter = 0;
	for (i = 0; i < num_dist; i++) *(label+i) = 0;	

	for (i = 0; i < num_dist; i++) 
		if (*(label+i) == 0) {
			(residues+counter)->ligand_ser_no = (distances+i)->ligand_ser_no;
			(residues+counter)->target_atom_ser_no = (distances+i)->ligand_atom_ser_no;
			strcpy((residues+counter)->target_res_type,(distances+i)->ligand_res_type);
			strcpy((residues+counter)->target_chain,(distances+i)->ligand_chain);
			(residues+counter)->target_res_ser = (distances+i)->ligand_res_ser;
			for (j = 0; j < num_dist; j++)
				if ((distances+i)->ligand_ser_no == (distances+j)->ligand_ser_no && strcmp((distances+i)->ligand_res_type,(distances+j)->ligand_res_type) == 0 && strcmp((distances+i)->ligand_chain,(distances+j)->ligand_chain) == 0 && (distances+i)->ligand_res_ser == (distances+j)->ligand_res_ser)
					*(label+j) = 1;
			counter++;
			}
	free(label);
	return residues;
	}

hs *edge_lister (pd *distances, int num_distances, int *list_length, char ptr_type) {
	
	int i, j, indicator, counter;
	hs *edge_list = NULL;
	hs *x = NULL;	

	counter = 0;

	for ( i = 0; i < num_distances; i++ ) {
		indicator = 0;
		if ( counter == 0 && (distances+i)->target_dist != 0 ) edge_list = (hs *) malloc(sizeof(hs));
		for ( j = 0; j < counter; j++ ) {
			if ( (distances+i)->ligand_res_ser == (edge_list+j)->wat_res_ser && (distances+i)->target_res_ser == (edge_list+j)->ptr_res_ser /* && strcmp ((distances+i)->target_chain, (edge_list+j)->ptr_chn) == 0 */ ) indicator++;	// removes duplicate connections to a w,l,t,b node
			if ( ptr_type == 'w' && (distances+i)->ligand_res_ser == (edge_list+j)->ptr_res_ser && (distances+i)->target_res_ser == (edge_list+j)->wat_res_ser ) indicator++;	// removes duplicate w-w edges --- NB. chain ID is not considered at waters, residue ser. number must be different
			}
		if ( ( indicator == 0 || counter == 0 ) && (distances+i)->target_dist != 0 ) { 			// removes self-connections
			counter++;							
			if ( counter != 1 ) { 
				x = realloc( edge_list, counter * sizeof(hs) );
				if (x) edge_list = x;
				x = NULL;
				}
			(edge_list+counter-1)->wat_atm_ser = (distances+i)->ligand_atom_ser_no;
			(edge_list+counter-1)->wat_res_ser = (distances+i)->ligand_res_ser;
			(edge_list+counter-1)->ptr_atm_ser = (distances+i)->target_atom_ser_no;
			(edge_list+counter-1)->ptr_res_ser = (distances+i)->target_res_ser;
			strcpy ((edge_list+counter-1)->wat_res_typ, (distances+i)->ligand_res_type);
			strcpy ((edge_list+counter-1)->ptr_atm_typ, (distances+i)->target_pdb_type);
			strcpy ((edge_list+counter-1)->ptr_res_typ, (distances+i)->target_res_type);
			strcpy ((edge_list+counter-1)->ptr_chn, (distances+i)->target_chain);
			strcpy ((edge_list+counter-1)->wat_chn, (distances+i)->ligand_chain);
			(edge_list+counter-1)->wat_b_factor = (distances+i)->ligand_b_factor;			// here mobility
			(edge_list+counter-1)->ptr_b_factor = (distances+i)->target_b_factor;			// here mobility
			(edge_list+counter-1)->ptr_type = ptr_type;						// gives ptr type
			(edge_list+counter-1)->dist = (distances+i)->target_dist;
			}
//			else continue;
		}
	*list_length = counter;
	return (edge_list);
	}

hs *edge_filter_C (hs *edge_list, int num_edges, int *filt_num_edges_C) {

	int i, counter;
	hs *edge_list_C = NULL;
	hs *x = NULL;

	counter = 0;	

	for ( i = 0; i < num_edges; i++ )		
		if ( (edge_list+i)->ptr_atm_typ[0] != 'C' ) {		// remove edges to nodes of Carbon atoms
			counter++;
			if ( counter == 1 ) edge_list_C = (hs *) malloc(sizeof(hs));
				else { 
					x = realloc( edge_list_C, counter * sizeof(hs) );
					if (x) edge_list_C = x;
					x = NULL;
					}
			*(edge_list_C+counter-1) = *(edge_list+i);
			}
	*filt_num_edges_C = counter;
	return (edge_list_C);
	}

hs *edge_filter_4 (hs *edge_list, int num_edges, int *filt_num_edges_4) {

	int i, j, k, counter, counter_2, main_counter, wat_as_partner, wat_as_partner_2, available;
	hs *edge_list_4 = NULL;
	hs *x = NULL;
	hs *temp_list = NULL;
	hs *y = NULL;

	void insertion_sort_hs_dist (hs *a, int num);

	main_counter = 0;	

	for ( i = 0; i < num_edges; i++ ) {
		if ( (edge_list+i)->used != 'y' ) {
			counter = 0;
			for ( j = 0; j < num_edges; j++ )
				if ( (edge_list+i)->wat_atm_ser == (edge_list+j)->wat_atm_ser && (edge_list+i)->wat_res_ser == (edge_list+j)->wat_res_ser ) {
					counter++;
					if ( counter == 1 ) temp_list = (hs *) malloc(sizeof(hs));
						else {
							y = realloc( temp_list, counter * sizeof(hs) );
							if (y) temp_list = y;
							y = NULL;
							}
					*(temp_list+counter-1) = *(edge_list+j);
					(edge_list+j)->used = 'y';
					(temp_list+counter-1)->used = 'y';
					}

			if ( counter > 0 ) insertion_sort_hs_dist (temp_list, counter);		// order by increasing distances
			wat_as_partner = 0;	
			if ( main_counter != 0 ) {			
				for ( j = 0; j < main_counter; j++ )		// count how many times a WAT node was present before as a partner node in w-w edges
					if ( (edge_list+i)->wat_atm_ser == (edge_list_4+j)->ptr_atm_ser && (edge_list+i)->wat_res_ser == (edge_list_4+j)->ptr_res_ser )				
							wat_as_partner++;
				for ( k = 0; k < counter; k++ )	{	
					wat_as_partner_2 = 0;
					for ( j = 0; j < main_counter; j++ )		
						if ( (temp_list+k)->ptr_atm_ser == (edge_list_4+j)->ptr_atm_ser && (temp_list+k)->ptr_res_ser == (edge_list_4+j)->ptr_res_ser )				
							wat_as_partner_2++;
					if ( wat_as_partner_2 == 4 ) 		// mark edges with partner WAT nodes already present in 4 edges as partner
						(temp_list+k)->used = 'n';
					}
				}

			available = 4 - wat_as_partner; 		// no more than 4 partner allowed for one WAT, available can be zero
			counter_2 = 0; 					

			if ( main_counter == 0 && counter > 0 ) {
				edge_list_4 = (hs *) malloc(sizeof(hs));
				*(edge_list_4) = *(temp_list);
				counter_2++;
				for ( k = 1; k < counter; k++ )	
					if ( counter_2 <= 4 && (temp_list+k)->used == 'y' ) {		// no more than 4 partner allowed for one WAT
						counter_2++;
						x = realloc( edge_list_4, counter_2 * sizeof(hs) );
						if (x) edge_list_4 = x;
						x = NULL;
						*(edge_list_4+counter_2-1) = *(temp_list+k);
						}
				}

			if ( main_counter > 0 && counter > 0 && available > 0 ) {	
				for ( k = 0; k < counter; k++ )	
					if ( counter_2 < available && (temp_list+k)->used == 'y' ) {		// no more than 4 partner allowed for one WAT
						counter_2++;
						x = realloc( edge_list_4, (main_counter+counter_2) * sizeof(hs) );
						if (x) edge_list_4 = x;
						x = NULL;
						*(edge_list_4+main_counter+counter_2-1) = *(temp_list+k);
						}
				}

			main_counter = main_counter + counter_2;
			if ( counter > 0 ) free (temp_list);
			}
		}
	*filt_num_edges_4 = main_counter;
	return (edge_list_4);
	}

hd *list_nodes (hs *edge_list, int num_edges, int *num_nodes) {

	int i, j, indicator_p, indicator_w, counter;
	hd *nodes = NULL;
	hd *x = NULL;

	counter = 0;	
	
	for ( i = 0; i < num_edges; i++ ) {
		if ( i == 0 ) nodes = (hd *) malloc(sizeof(hd));
		indicator_w = 0;
		indicator_p = 0;		
		for ( j = 0; j < counter; j++ ) {
			if ( (edge_list+i)->wat_res_ser == (nodes+j)->res_ser && strncmp ((edge_list+i)->wat_res_typ, (nodes+j)->res_typ, 3) == 0 ) indicator_w = 1;
			if ( (edge_list+i)->ptr_res_ser == (nodes+j)->res_ser && strncmp ((edge_list+i)->ptr_res_typ, (nodes+j)->res_typ, 3) == 0 && strcmp ((edge_list+i)->ptr_chn, (nodes+j)->chain) == 0 ) indicator_p = 1;
			}

		if ( indicator_w == 0 || i == 0 ) {
			counter++;
			if ( counter > 1 ) {
				x = realloc( nodes, counter * sizeof(hd) );
				if (x) nodes = x;
				x = NULL;
				}
			(nodes+counter-1)->atm_ser = (edge_list+i)->wat_atm_ser;
			(nodes+counter-1)->res_ser = (edge_list+i)->wat_res_ser;
			strcpy((nodes+counter-1)->res_typ, (edge_list+i)->wat_res_typ);
			strcpy((nodes+counter-1)->chain, (edge_list+i)->wat_chn);
			(nodes+counter-1)->type = 'w';
			(nodes+counter-1)->b_factor = (edge_list+i)->wat_b_factor;
			}

		if ( indicator_p == 0 ) {
			counter++;
			if ( counter > 1 ) {
				x = realloc( nodes, counter * sizeof(hd) );
				if (x) nodes = x;
				x = NULL;
				}
			(nodes+counter-1)->atm_ser = (edge_list+i)->ptr_atm_ser;
			(nodes+counter-1)->res_ser = (edge_list+i)->ptr_res_ser;
			strcpy((nodes+counter-1)->res_typ, (edge_list+i)->ptr_res_typ);
			strcpy((nodes+counter-1)->chain, (edge_list+i)->ptr_chn);
			(nodes+counter-1)->type = (edge_list+i)->ptr_type;
			(nodes+counter-1)->b_factor = (edge_list+i)->ptr_b_factor;		// used as mobility for w only
			}

		}
	*num_nodes = counter;
	return (nodes);
	}

void count_degree (hs *edge_list, int num_edges, hd *node_list, int num_nodes) {

	int i, j, indicator;

	for ( i = 0; i < num_nodes; i++ )
		(node_list+i)->degree = 0;

	// total degree
	for ( i = 0; i < num_nodes; i++ )
		for ( j = 0; j < num_edges; j++ ) {
			indicator = 0;
			if ( (edge_list+j)->wat_res_ser == (node_list+i)->res_ser && strncmp ((edge_list+j)->wat_res_typ, (node_list+i)->res_typ, 3) == 0 ) indicator = 1;
			if ( (edge_list+j)->ptr_res_ser == (node_list+i)->res_ser && strncmp ((edge_list+j)->ptr_res_typ, (node_list+i)->res_typ, 3) == 0 && strcmp ((edge_list+j)->ptr_chn, (node_list+i)->chain) == 0 ) indicator = 1;
			if (indicator == 1)
				(node_list+i)->degree++;
			}

	// solute degree
	for ( i = 0; i < num_nodes; i++ ) {
		(node_list+i)->sol_degree = 0;
		for ( j = 0; j < num_edges; j++ )
			if ( (edge_list+j)->wat_res_ser == (node_list+i)->res_ser && strncmp ((edge_list+j)->wat_res_typ, (node_list+i)->res_typ, 3) == 0 && ( (edge_list+j)->ptr_type == 't' || (edge_list+j)->ptr_type == 'l' ) )
				(node_list+i)->sol_degree++;
		}

	}

void set_static (hs *edge_list, int num_edges, hd *node_list, int num_nodes) {

	int i, j;

	for ( i = 0; i < num_nodes; i++ )	// mark static nodes
		if ( (node_list+i)->type == 'w' && ( (node_list+i)->degree == 4 || (node_list+i)->sol_degree >= 3 || (node_list+i)->b_factor <= g_m_tolerance ) ) (node_list+i)->node_dyn = 's';
			else if ( (node_list+i)->type == 'l' || (node_list+i)->type == 't' ) (node_list+i)->node_dyn = 's';	// l,t partners are static by def.
				else (node_list+i)->node_dyn = 'd';		// bulk is dynamic by def.

	for ( i = 0; i < num_edges; i++ )		// mark deft. dynamic edges
		(edge_list+i)->edge_dyn = 'd';

	for ( i = 0; i < num_edges; i++ )		// mark wat nodes
		for ( j = 0; j < num_nodes; j++ )
			if ( (edge_list+i)->wat_res_ser == (node_list+j)->res_ser && strncmp ((edge_list+i)->wat_res_typ, (node_list+j)->res_typ, 3) == 0 )
				(edge_list+i)->wat_dyn = (node_list+j)->node_dyn;
	
	for ( i = 0; i < num_edges; i++ )		// mark partners
		for ( j = 0; j < num_nodes; j++ )
			if ( (edge_list+i)->ptr_res_ser == (node_list+j)->res_ser && strncmp ((edge_list+i)->ptr_res_typ, (node_list+j)->res_typ, 3) == 0 && strcmp ((edge_list+i)->ptr_chn, (node_list+j)->chain) == 0 )
				(edge_list+i)->ptr_dyn = (node_list+j)->node_dyn;

	for ( i = 0; i < num_edges; i++ )		// mark edges if static
		if ( (edge_list+i)->wat_dyn == 's' && (edge_list+i)->ptr_dyn == 's' )			
			(edge_list+i)->edge_dyn = 's';

	}

