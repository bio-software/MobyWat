// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "global.h"

extern ir *g_input_back;
extern int g_ref_target_ca_no;
extern int g_ref_target_no;
extern int g_ref_ligand_no;
extern char g_m;

double calc_lowb_no(re *ref_waters, int num_ref_waters, double b_fact_threshold) {

	int i;
	double number_of_lowb_wats = 0;

	for(i = 0; i < num_ref_waters; i++) 
		if ( (ref_waters+i)->b_factor < b_fact_threshold ) number_of_lowb_wats++;

	return (number_of_lowb_wats);
	}


int calc_success(mr *ml, int rows) {

	int i, j, label = 0;
	int number_of_matches = 0;

	for(i = 0; i < rows; i++) {
		label = 0;			// Calculate each reference waters as matches only once
		for(j = 0; j < i; j++)
			if ( (ml+i)->indicator == (ml+j)->indicator ) label = 1;
		if ( (ml+i)->indicator != 0 && label == 0 ) number_of_matches++;
		}

	return (number_of_matches);
	}


double dist_calc (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch) {

	double distance;	
	
	if ((pdb_type1[0]  != 'H' && pdb_type2[0]  != 'H') || h_switch == 'h') {
		distance = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
		return distance;
		}
		else return VERY_BIG_DISTANCE;	
	}

double dist_sqr_calc (double x1, double y1, double z1, double x2, double y2, double z2) {

	return ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
	}

void mean_coords3(pl *a, int members, re *mean) {

	int k;

	mean->x_coord = 0;
	mean->y_coord = 0;
	mean->z_coord = 0;

	for(k = 0; k < members; k++) {
		mean->x_coord = mean->x_coord + (a+k)->x_coord;
		mean->y_coord = mean->y_coord + (a+k)->y_coord;
		mean->z_coord = mean->z_coord + (a+k)->z_coord;
		}

	mean->x_coord = mean->x_coord / members;
	mean->y_coord = mean->y_coord / members;
	mean->z_coord = mean->z_coord / members;

	}

void mean_coords2(re *a, int members, re *mean) {

	int k;

	mean->x_coord = 0;
	mean->y_coord = 0;
	mean->z_coord = 0;

	for(k = 0; k < members; k++) {
		mean->x_coord = mean->x_coord + (a+k)->x_coord;
		mean->y_coord = mean->y_coord + (a+k)->y_coord;
		mean->z_coord = mean->z_coord + (a+k)->z_coord;
		}

	mean->x_coord = mean->x_coord / members;
	mean->y_coord = mean->y_coord / members;
	mean->z_coord = mean->z_coord / members;

	}

void mean_coords(re ***a, int **members, int *clusters, int rows, re **mean) {

	int i, j, k;

	for(i = 0; i < rows; i++)
		for(j = 0; j < *(clusters+i); j++) {
			(*(mean+i)+j)->x_coord = 0;
			(*(mean+i)+j)->y_coord = 0;
			(*(mean+i)+j)->z_coord = 0;
			for(k = 0; k < *(*(members+i)+j); k++) {
				(*(mean+i)+j)->x_coord = (*(mean+i)+j)->x_coord + (*(*(a+i)+j)+k)->x_coord;
				(*(mean+i)+j)->y_coord = (*(mean+i)+j)->y_coord + (*(*(a+i)+j)+k)->y_coord;
				(*(mean+i)+j)->z_coord = (*(mean+i)+j)->z_coord + (*(*(a+i)+j)+k)->z_coord;
				}
			(*(mean+i)+j)->x_coord = (*(mean+i)+j)->x_coord / *(*(members+i)+j);
			(*(mean+i)+j)->y_coord = (*(mean+i)+j)->y_coord / *(*(members+i)+j);
			(*(mean+i)+j)->z_coord = (*(mean+i)+j)->z_coord / *(*(members+i)+j);
			(*(mean+i)+j)->target_atom_ser_no = (*(*(a+i)+j)+0)->target_atom_ser_no;
			(*(mean+i)+j)->target_res_ser = (*(*(a+i)+j)+0)->target_res_ser;
			}

	}

int min_dbl_reno(double *a, int n) {

	int i;
	int min_no = 0;
	double temp;
	temp = *a;

	for(i = 0; i < n; i++)
		if (temp > *(a+i)) {
			temp = *(a+i);
			min_no = i;
			}
	return min_no;
	}

// 1 for a target RMSD atom: named CA, or $MOBYWAT_RMS_ATOM (e.g. C4' for nucleic acids), but not an ion residue of the same name (calcium CA)
int is_rms_atom (char *pdb_type, char *res_type) {

	static char *name = NULL;
	char res [RES_TYPE_LENGTH+1] = "";

	if ( name == NULL ) {
		name = getenv("MOBYWAT_RMS_ATOM");
		if ( name == NULL || name[0] == '\0' ) name = "CA";
		}
	if ( strcmp(pdb_type, name) != 0 ) return 0;
	sscanf(res_type, "%3s", res);
	return ( strcmp(res, pdb_type) != 0 );
	}

double rmsd_calc (ap *reference, ap *actual, int atom_number, char a_type_sw) {

 	int i, j, num_ignore;
	int ref_atom_number = 0;
	int ca_counter_ref = 0;
	int ca_counter_act = 0;
	int *ca_serial = NULL;
	double summa, rmsd;
	num_ignore = 0;
	summa = 0;

	if (g_m == 'p') ref_atom_number = atom_number;		// else refers to analysis mode
		else ref_atom_number = g_ref_target_no;

	switch (a_type_sw) {
	
		case 'a':	for (i = 0; i < ref_atom_number; i++)		// Calphas - does not require equal number of atoms only eq. no. of CA-s
					if (is_rms_atom((reference+i)->pdb_type, (reference+i)->res_type)) ca_counter_ref++;
				for (i = 0; i < atom_number; i++)		
					if (is_rms_atom((actual+i)->pdb_type, (actual+i)->res_type)) ca_counter_act++;
				g_input_back->target_traj_ca_no = ca_counter_act;
				if ( ca_counter_ref != ca_counter_act ) {
					printf("%s","ERROR! Numbers of CA atoms are not equal in the two compared target molecules! RMSD calculation is impossible!\n");
					printf("%s%6d%s","Number of CA atoms in target 1 (ref.):",ca_counter_ref,"\n");
					printf("%s%6d%s","Number of CA atoms in target 2 (act.):",ca_counter_act,"\n");
					exit (1);
					}
				ca_serial = (int *) malloc(ca_counter_ref * sizeof(int));	
				j = 0;
				for (i = 0; i < ref_atom_number; i++)
					if (is_rms_atom((reference+i)->pdb_type, (reference+i)->res_type)) {
						*(ca_serial+j) = i;
						j++;
						}
				j = 0;
				for (i = 0; i < atom_number; i++) {		
					if (is_rms_atom((actual+i)->pdb_type, (actual+i)->res_type)) {
						summa = summa + ((reference+*(ca_serial+j))->x_coord - (actual+i)->x_coord)*((reference+*(ca_serial+j))->x_coord - (actual+i)->x_coord) + ((reference+*(ca_serial+j))->y_coord - (actual+i)->y_coord)*((reference+*(ca_serial+j))->y_coord - (actual+i)->y_coord) + ((reference+*(ca_serial+j))->z_coord - (actual+i)->z_coord)*((reference+*(ca_serial+j))->z_coord - (actual+i)->z_coord);
						j++;
						}
					else num_ignore++;
					}
				free(ca_serial);
				break; 

		case 'h':	for (i = 0; i < atom_number; i++) {		// Non-hydrogen - requires equal number of atoms, but not checked here, used only between trajectory frames
					if ((actual+i)->pdb_type[0] != 'H')
					summa = summa + ((reference+i)->x_coord - (actual+i)->x_coord)*((reference+i)->x_coord - (actual+i)->x_coord) + ((reference+i)->y_coord - (actual+i)->y_coord)*((reference+i)->y_coord - (actual+i)->y_coord) + ((reference+i)->z_coord - (actual+i)->z_coord)*((reference+i)->z_coord - (actual+i)->z_coord);
					else num_ignore++;
					}
				break; 

		case 'l':	for (i = 0; i < atom_number; i++) 		// All
					summa = summa + ((reference+i)->x_coord - (actual+i)->x_coord)*((reference+i)->x_coord - (actual+i)->x_coord) + ((reference+i)->y_coord - (actual+i)->y_coord)*((reference+i)->y_coord - (actual+i)->y_coord) + ((reference+i)->z_coord - (actual+i)->z_coord)*((reference+i)->z_coord - (actual+i)->z_coord);
				num_ignore = 0;
				break; 


	    	}

	rmsd = sqrt (summa / (atom_number - num_ignore));
	return rmsd;	

	}

int wat_res_comp(char sample[]) {
	static char *wat_res_name[] = {"WAT", "SOL", "H2O", "HOH"};
	int indicator = 0;
	int i;
	for (i = 0; i < 4; i++) if ( strncmp(wat_res_name[i],sample,3) == 0 ) indicator = 1;
	return (indicator); 
	}

int wat_res_comp_excl(char sample[], char ref[]) {
	static char *wat_res_name[] = {"WAT", "SOL", "H2O", "HOH"};
	int indicator = 0;
	int i;
	for (i = 0; i < 4; i++) if ( strncmp(wat_res_name[i],sample,3) == 0 ) indicator = 1;
	if ( strncmp(ref,sample,3) == 0 ) indicator--;
	return (indicator); 
	}

int ch_comp(char sample[], char *chains, int no_chains) {
	int indicator = 0;
	int i;
	for (i = 0; i < no_chains; i++) if ( sample[0] == chains[i] ) indicator = 1;
	return (indicator); 
	}

double median(double *a, int n) {

	double temp;
	int i,j;

	for(i = 0; i < n; i++)
		for(j = i+1; j < n; j++) {
			if(*(a+i) > *(a+j)) {
				temp = *(a+j);
				*(a+j) = *(a+i);
				*(a+i) = temp;
				}
		}
	if(n%2 == 0)
		return (*(a+(n/2-1)) + *(a+n/2))/2;
	else
		return *(a+n/2);
	}

double mad(double *a, int n) {

	double median_value;
	int i;
	double *adeviation = NULL;
	adeviation = (double *) malloc(n*sizeof(double));

	median_value = median (a, n);
	for(i = 0; i < n; i++)
		*(adeviation+i) = fabs(*(a+i) - median_value);	
	return (median(adeviation, n));

	free(adeviation);
	}

double mean(double *a, int n) {

	int i;
	double sum = 0;

	for(i = 0; i < n; i++)
		sum = sum + *(a+i);
	return (sum / n);
	}

double sdev(double *a, int n) {

	int i;
	double mean_value, sumsq = 0;
	double *sq_diff;	
	sq_diff = (double *) malloc(n*sizeof(double));

	mean_value = mean(a, n);
	for(i = 0; i < n; i++)
		*(sq_diff+i) = pow((*(a+i) - mean_value),2);
	for(i = 0; i < n; i++)
		sumsq = sumsq + *(sq_diff+i);
	return (sqrt(sumsq/(n-1)));

	free(sq_diff);
	}

int max_int(int *b, int n) {

	int i, temp;
	temp = *b;

	for(i = 0; i < n; i++)
		if ( temp < *(b+i) ) 
			temp = *(b+i);
	return temp;
	}

double min_dbl(double *a, int n) {

	int i;
	double temp;
	temp = *a;

	for(i = 0; i < n; i++)
		if (temp > *(a+i)) 
			temp = *(a+i);
	return temp;
	}

double max_dbl(double *a, int n) {

	int i;
	double temp;
	temp = *a;

	for(i = 0; i < n; i++)
		if (temp < *(a+i)) 
			temp = *(a+i);
	return temp;
	}

double sum(double *a, int n) {

	int i;
	double sum = 0;

	for(i = 0; i < n; i++)
		sum = sum + *(a+i);
	return (sum);
	}

void b_factorize_target (ap *pdb_target, pd *distances, int num_target_atoms, int num_distances) {

	int i, j;

	for (i = 0; i < num_target_atoms; i++) {
		(pdb_target+i)->b_factor = MIN_B_FACTOR;
		(pdb_target+i)->occ = MIN_B_FACTOR;
		}

	for (i = 0; i < num_distances; i++)
		for (j = 0; j < num_target_atoms; j++) 
			if (strcmp((distances+i)->target_chain, (pdb_target+j)->chain) == 0 && (distances+i)->target_res_ser == (pdb_target+j)->res_ser) {
				(pdb_target+j)->occ = MAX_OCC;
				if ((distances+i)->target_atom_ser_no == (pdb_target+j)->atom_ser) (pdb_target+j)->b_factor = MAX_B_FACTOR;
				}
	}

int max_int_int(int *b, int n) {

	int i, temp;
	temp = *b;

	for(i = 0; i < n; i++)
		if (temp < *(b+i)) 
			temp = *(b+i);
	return temp;
	}

