// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "global.h"

extern int g_ref_target_ca_no;
extern char g_temp [MAX_FILENAME_LENGTH];
extern char g_log_file_name [MAX_FILENAME_LENGTH];
extern char g_e;
extern int g_frame_range_min;
extern int g_frame_range_max;
extern ir *g_input_back;
extern char g_v;

void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );
int wat_res_comp(char sample[]);
int wat_res_comp_excl(char sample[], char ref[]);

ap *box_to_pdb (ap *box_pdb, char limits [MAX_FILENAME_LENGTH], int *atom_num_pdb, int atom_num_box, char label) {

	int j, lower, upper, counter, expected_atom_num;	
	int no_chains = 0;
	char *delims = ":-\n[]";
	char temp1 [MAX_FILENAME_LENGTH]="";
	char temp2 [MAX_FILENAME_LENGTH]="";
	char def_mode = 'd';		// default: atom serial numbers in an xxx-xxx fashion
	char limits_copy [MAX_FILENAME_LENGTH] = "";

	ap *molecule = NULL;
	ap *x = NULL;

	int wat_res_comp(char sample[]);
	int ch_comp(char sample[], char *chains, int no_chains);

	void str0 (char string [], int length);

	str0 (temp1, MAX_FILENAME_LENGTH);
	str0 (temp2, MAX_FILENAME_LENGTH);

	strcpy(limits_copy, limits);

	if ( strncmp(limits,AUTO_KEYW, 4) == 0 && label == 'w' ) def_mode = 'a';		// a: auto, o: only RES
		else if ( label == 'w' ) def_mode = 'o';

	if ( strncmp(limits,AUTO_KEYW, 4) == 0 && label == 'b' ) def_mode = 'w';		// w: only WAT (auto mode here), e: exclude RES
		else if ( label == 'b' ) def_mode = 'e';

	if ( limits [0] == '[' && ( label == 't' || label == 'l' ) ) {				// chains mode (for target and ligand)
		def_mode = 'c';
		strcpy(temp1, strtok(limits_copy, delims));
		no_chains = strlen(temp1);
		}

	if ( def_mode == 'd' ) {
		strcpy(temp1, strtok(limits_copy, delims));					// default: d 
		strcpy(temp2, strtok(NULL, delims));
		}
	
	if ( def_mode == 'c' && strlen(temp1) == 0 ) {						// check range in c
		printf("ERROR! Chain ID range is erroneous or missing.\n"); 
		exit (1);
		}

	if ( def_mode == 'd' && (strlen(temp1) == 0 || strlen(temp2) == 0) ) {			// check range in d
		printf("ERROR! Atom serial number range is erroneous or missing.\n"); 
		exit (1);
		}
	
	if ( def_mode == 'd' ) {								// default: d check
		lower = atoi(temp1);
		upper = atoi(temp2);
		if (upper <= lower) {
			printf("%s%s%s","ERROR! Upper atom number limit is smaller than/equal to the lower limit here: ",limits,"\n"); 
			exit (1);
			}
		expected_atom_num = upper - lower + 1;
		}

	counter = 0;
	j = 0;	

	while ( j < atom_num_box ) {
		if ( ( def_mode == 'd' && (box_pdb+j)->atom_ser >= lower && (box_pdb+j)->atom_ser <= upper ) || ( def_mode == 'a' && wat_res_comp((box_pdb+j)->res_type) == 1 ) || ( def_mode == 'e' && wat_res_comp_excl((box_pdb+j)->res_type,limits) == 1 ) || ( def_mode == 'o' && strncmp((box_pdb+j)->res_type,limits, 3) == 0 ) || ( def_mode == 'w' && strncmp((box_pdb+j)->res_type,WAT_KEYW, 3) == 0 ) || ( def_mode == 'c' && ch_comp((box_pdb+j)->chain, temp1, no_chains) == 1 && wat_res_comp((box_pdb+j)->res_type) == 0 ) ) {
			if ( counter == 0) molecule = (ap *) malloc(sizeof(ap));
				else {
					x = realloc( molecule, (counter+1) * sizeof(ap) );
					if (x) molecule = x;
					x = NULL;
					}
			(molecule+counter)->model_ser = (box_pdb+j)->model_ser;
			strcpy((molecule+counter)->pdb_token, (box_pdb+j)->pdb_token);			
			(molecule+counter)->pdb_token [PDB_TOKEN_LENGTH] = '\0';
			(molecule+counter)->atom_ser = (box_pdb+j)->atom_ser;
			strcpy((molecule+counter)->pdb_type, (box_pdb+j)->pdb_type);			
			(molecule+counter)->pdb_type [PDB_TYPE_LENGTH] = '\0';
			strcpy((molecule+counter)->pdb_alt_loc, (box_pdb+j)->pdb_alt_loc);			
			(molecule+counter)->pdb_alt_loc [PDB_ALT_LOC_LENGTH] = '\0';
			strcpy((molecule+counter)->res_type, (box_pdb+j)->res_type);
			(molecule+counter)->res_type [RES_TYPE_LENGTH] = '\0';
			strcpy((molecule+counter)->chain, (box_pdb+j)->chain);
			(molecule+counter)->chain [CHAIN_LENGTH] = '\0';
			(molecule+counter)->res_ser = (box_pdb+j)->res_ser;			
			strcpy((molecule+counter)->pdb_achar, (box_pdb+j)->pdb_achar);
			(molecule+counter)->pdb_achar [PDB_ACHAR_LENGTH] = '\0';
			(molecule+counter)->x_coord = (box_pdb+j)->x_coord;
			(molecule+counter)->y_coord = (box_pdb+j)->y_coord;
			(molecule+counter)->z_coord = (box_pdb+j)->z_coord;
			(molecule+counter)->occ = (box_pdb+j)->occ;
			(molecule+counter)->b_factor = (box_pdb+j)->b_factor;
			counter++;
			}
		j++;
		}

	*(int *)atom_num_pdb = counter;	

	if ( counter == 0 ) {
		printf ("%s%c%s%s%s","ERROR! Range \"-",label," ",limits,"\" defines zero atoms.\n");
		exit (0);
		}

	if ( def_mode == 'd' ) {
		if ( *(int *)atom_num_pdb != expected_atom_num ) {
			printf("WARNING! There are missing/extra atoms in the following atom number range.\n");
			printf("%s%d%s", "WARNING! Lower atom range limit:   ",lower,"\n"); 
			printf("%s%d%s", "WARNING! Upper atom range limit:   ",upper,"\n"); 
			printf("%s%d%s", "WARNING! Found atom number in pdb: ",*(int *)atom_num_pdb,"\n"); 
			printf("%s%d%s", "WARNING! Difference (found-range): ",*(int *)atom_num_pdb - expected_atom_num,"\n");
			}
		}

	return molecule;
	}


ap *split_ref_file (ap *box_pdb, char **limits, int no_limits, int *atom_num_total, int atom_num_box, char label) {

	int i;
	int atom_num_temp = 0;
	int *atom_num_part = NULL;
	int ca_counter_ref = 0;

	ap **output_parts = NULL;
	ap *full_output = NULL;

//	ap *merge_mols (ap *mol1, ap *mol2, int atom_num1, int atom_num2, int *atom_numo);

	output_parts = (ap **) malloc( no_limits * sizeof(ap *) );
	atom_num_part = (int *) malloc( no_limits * sizeof(int) );
	
	ap *box_to_pdb (ap *box_pdb, char limits [MAX_FILENAME_LENGTH], int *atom_num_pdb, int atom_num_box, char label);

	*atom_num_total = 0;

	for (i = 0; i < no_limits; i++)
		*(output_parts+i) = box_to_pdb(box_pdb, *(limits+i), atom_num_part+i, atom_num_box, label);

	for (i = 0; i < no_limits; i++)
		*atom_num_total = *atom_num_total + *(atom_num_part+i);

	full_output = (ap *) malloc( *atom_num_total * sizeof(ap) );

	for (i = 0; i < no_limits; i++) {
		memcpy (full_output+atom_num_temp, *(output_parts+i), *(atom_num_part+i)*sizeof(ap));
		atom_num_temp = atom_num_temp + *(atom_num_part+i);
		}

	if (label == 't') {
		for (i = 0; i < *atom_num_total; i++)		// count Calphas 
			if ((full_output+i)->pdb_type[0]  == 'C' && (full_output+i)->pdb_type[1]  == 'A') ca_counter_ref++;
		g_ref_target_ca_no = ca_counter_ref;
		}

	free(atom_num_part);
	
	for (i = 0; i < no_limits; i++) free( *(output_parts+i) );
	free(output_parts);

	return(full_output);

	}

ap *merge_mols (ap *mol1, ap *mol2, int atom_num1, int atom_num2, int *atom_numo) {
	
	ap *molo = NULL;
	*atom_numo = atom_num1 + atom_num2;

	molo = (ap *) malloc(*atom_numo*sizeof(ap));

	memcpy (molo, mol1, atom_num1*sizeof(ap));
	memcpy (molo+atom_num1, mol2, atom_num2*sizeof(ap));

	return molo;
	}


hs *merge_edges (hs *mol1, hs *mol2, int atom_num1, int atom_num2, int *atom_numo) {
	
	hs *molo = NULL;
	*atom_numo = atom_num1 + atom_num2;

	molo = (hs *) malloc(*atom_numo*sizeof(hs));

	memcpy (molo, mol1, atom_num1*sizeof(hs));
	memcpy (molo+atom_num1, mol2, atom_num2*sizeof(hs));

	return molo;
	}


