// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 03-12-2015

#include "print.h"

void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag ) {

	FILE *outfile;

	if ( flag == 'w') outfile = fopen(file_name,"w");
		else outfile = fopen(file_name,"a");

	fprintf(outfile,"%s",text);

	fclose(outfile);
	}

void print_all_clust_waters(re ***a, int **members, int *clusters, int rows, char file_name [MAX_FILENAME_LENGTH], char i_label) {

	FILE *outfile;
	int i, j, k;
	int counter_atom = 1;
	int counter_residue = 1;

	outfile=fopen(file_name,"w");

	if (i_label == 'i') fprintf(outfile,"%s",ALL_CLUST_WAT_FILE_I_HEADER);
		else fprintf(outfile,"%s",ALL_CLUST_WAT_FILE_S_HEADER);

	fprintf(outfile,"%s",ALL_CLUST_WAT_FILE_HEADER);

	for (i = 0; i < rows; i++) {
		for (j = 0; j < *(clusters+i); j++) {
			fprintf(outfile,"%s%6d%s%6d%s","REMARK   Water#: ",i+1,"   Cluster#: ",j+1,"\n");
			for (k = 0; k < *(*(members+i)+j); k++) {
				fprintf(outfile,"%-6s","HETATM"); 
				fprintf(outfile,"%5d",counter_atom); 
				fprintf(outfile,"%1s",""); 	// just space
				fprintf(outfile,"%-4s","OW"); 	// atom type
				fprintf(outfile,"%1s",""); 	// alt loc
				fprintf(outfile,"%3s","WAT"); 	// res type
				fprintf(outfile,"%1s",""); 	// just space
				fprintf(outfile,"%1s",""); 	// chain ID
				fprintf(outfile,"%4d",counter_residue); 	// res ser = cluster number
				fprintf(outfile,"%1s",""); 	// achar
				fprintf(outfile,"%3s",""); 	// just space
				fprintf(outfile,"%8.3f",(*(*(a+i)+j)+k)->x_coord); 
				fprintf(outfile,"%8.3f",(*(*(a+i)+j)+k)->y_coord); 
				fprintf(outfile,"%8.3f",(*(*(a+i)+j)+k)->z_coord); 
				fprintf(outfile,"%6d",(*(*(a+i)+j)+k)->target_atom_ser_no); 	// occ = original_atom_ser
				fprintf(outfile,"%6d",(*(*(a+i)+j)+k)->ligand_ser_no); 		// b-factor = real frame ser no
				fprintf(outfile,"%s","\n");
				counter_atom++;
				if ( counter_atom == 100000 ) counter_atom = 1;
		    		}
			counter_residue++;
			if ( counter_residue == 10000 ) counter_residue = 1;		
			}
		}

	fprintf(outfile,"%s","END\n");
	fclose(outfile);		
	}


void print_all_IF_waters(re **a, int frames, int *total_if, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min) {

	FILE *outfile;
	int i, j;
	int counter_atom = 1;
	int counter_residue = 0;

	outfile=fopen(file_name,"w");

	// Project-specific header
	switch (flag) {
		case 'i': fprintf(outfile,"%s",TRAJ_ALLIFACE_FILE_HEADER);
		break; 
		case 's': fprintf(outfile,"%s",TRAJ_ALLSFACE_FILE_HEADER);
		break; 
	    	}


	for (i = 0; i <= frames; i++) {
		fprintf(outfile,"%s%6d%s%6d%s","REMARK   Frame#: ",i+frame_min,"   Count of IF waters in frame: ",*(total_if+i),"\n");
		for (j = 0; j < *(total_if+i); j++) {
			fprintf(outfile,"%-6s","HETATM"); 
			fprintf(outfile,"%5d",counter_atom); 
			fprintf(outfile,"%1s",""); 	// just space
			fprintf(outfile,"%-4s","OW"); 	// atom type
			fprintf(outfile,"%1s",""); 	// alt loc
			fprintf(outfile,"%3s","WAT"); 	// res type
			fprintf(outfile,"%1s",""); 	// just space
			fprintf(outfile,"%1s",""); 	// chain ID
			fprintf(outfile,"%4.4d",counter_residue); 
			fprintf(outfile,"%1s",""); 	// achar
			fprintf(outfile,"%3s",""); 	// just space
			fprintf(outfile,"%8.3f",(*(a+i)+j)->x_coord); 
			fprintf(outfile,"%8.3f",(*(a+i)+j)->y_coord); 
			fprintf(outfile,"%8.3f",(*(a+i)+j)->z_coord); 
			fprintf(outfile,"%6d",(*(a+i)+j)->target_atom_ser_no); 	// occ = original atom ser. number
			fprintf(outfile,"%6d",(*(a+i)+j)->target_res_ser); 	// b-factor = original residue ser. number
			fprintf(outfile,"%s","\n");
			counter_atom++;
			if ( counter_atom == 100000 ) counter_atom = 1;
	    		}		
		counter_residue++;
		if ( counter_residue == 10000 ) counter_residue = 1;
		}
	fprintf(outfile,"%s","END\n");
	fclose(outfile);		
	}


void print_clist_waters(re **a, int **members, int *clusters, int **ol, int rows, char file_name [MAX_FILENAME_LENGTH], char flag) {

	FILE *outfile;
	int i, j;
	int counter_atom = 1;
	int counter_residue = 1;

	outfile=fopen(file_name,"w");

	// Project-specific header
	switch (flag) {
		case 'i': fprintf(outfile,"%s",CLIST_IF_FILE_HEADER);
		break; 
		case 's': fprintf(outfile,"%s",CLIST_SF_FILE_HEADER);
		break; 
	    	}

	fprintf(outfile,"%s",CLIST_WAT_FILE_HEADER);

	for (i = 0; i < rows; i++) {
		fprintf(outfile,"%s%6d%s%6d%s","REMARK   Water#: ",i+1,"   Total occupacy: ",*(*(ol+2)+i),"\n");
		for (j = 0; j < *(clusters+i); j++) {
			fprintf(outfile,"%-6s","HETATM"); 
			fprintf(outfile,"%5d",counter_atom); 
			fprintf(outfile,"%1s",""); 	// just space
			fprintf(outfile,"%-4s","OW"); 	// atom type
			fprintf(outfile,"%1s",""); 	// alt loc
			fprintf(outfile,"%3s","WAT"); 	// res type
			fprintf(outfile,"%1s",""); 	// just space
			fprintf(outfile,"%1s",""); 	// chain ID
			fprintf(outfile,"%4d",counter_residue); 	// res ser = cluster number
			fprintf(outfile,"%1s",""); 	// achar
			fprintf(outfile,"%3s",""); 	// just space
			fprintf(outfile,"%8.3f",(*(a+i)+j)->x_coord); 
			fprintf(outfile,"%8.3f",(*(a+i)+j)->y_coord); 
			fprintf(outfile,"%8.3f",(*(a+i)+j)->z_coord); 
			fprintf(outfile,"%6d",j+1); 	// occ = cluster serial number
			fprintf(outfile,"%6d",*(*(members+i)+j)); 	// b-factor = partial occupacy of cluster 
			fprintf(outfile,"%s","\n");
			counter_atom++;
			if ( counter_atom == 100000 ) counter_atom = 1;
			}
		counter_residue++;
		if ( counter_residue == 10000 ) counter_residue = 1;		
		}

	fprintf(outfile,"%s","END\n");
	fclose(outfile);		
	}


void print_if_waters(re *a, int num_atom, char file_name [MAX_FILENAME_LENGTH], char swt) {

	FILE *outfile;
	int i;

	outfile=fopen(file_name,"w");

	if ( swt == 'i') fprintf(outfile,"%s",REF_IF_WAT_FILE_HEADER);
		else fprintf(outfile,"%s",REF_SF_WAT_FILE_HEADER);

	for (i = 0; i < num_atom; i++) 
		if ( (a+i)->b_factor < g_b_factor ) {
			fprintf(outfile,"%-6s","ATOM"); 
			fprintf(outfile,"%5d",(a+i)->target_atom_ser_no); 
			fprintf(outfile,"%1s",""); 	// just space
			fprintf(outfile,"%-4s","OW"); 	// atom type
			fprintf(outfile,"%1s",""); 	// alt loc
			fprintf(outfile,"%3s","WAT"); 	// res type
			fprintf(outfile,"%1s",""); 	// just space
			fprintf(outfile,"%1s",""); 	// chain ID
			fprintf(outfile,"%4d",(a+i)->target_res_ser); 	// res ser 
			fprintf(outfile,"%1s",""); 	// achar
			fprintf(outfile,"%3s",""); 	// just space
			fprintf(outfile,"%8.3f",(a+i)->x_coord); 
			fprintf(outfile,"%8.3f",(a+i)->y_coord); 
			fprintf(outfile,"%8.3f",(a+i)->z_coord); 
			fprintf(outfile,"%6d",i+1); 	// occ 
			fprintf(outfile,"%6.2f",(a+i)->b_factor); 	// b-factor 
			fprintf(outfile,"%s","\n");
			}

	fprintf(outfile,"%s","END\n");
	fclose(outfile);		
	}


char *print_pdb_file (ap *pdb, int atom_num, char file_out [MAX_FILENAME_LENGTH]) {

	FILE *outfile;
	int j;
	char *line = NULL;
	char ter [5] = "TER\n";
	char end [5] = "END\n";
	char header [MAX_LINE_LENGTH];

	char *print_pdb_line (ap *pdb, int index, char pdbqt_ind, int pdbqt_rank);

	outfile=fopen(file_out,"w");

	sprintf(header,"%s","REMARK Input coordinates.\n");
	fputs(header,outfile);

	for (j = 0; j < atom_num; j++) {
				line = print_pdb_line(pdb, j, ' ',j);
				fputs(line,outfile);
				free(line);
				if (j == atom_num-1) fputs(end,outfile);
					else if (strcmp ((pdb+j)->chain,(pdb+j+1)->chain) != 0) fputs(ter,outfile);
				}
	fclose(outfile);
	return(file_out);
	}


char *print_pdb_line (ap *pdb, int index, char pdbqt_ind, int pdbqt_rank) {

	char *line = NULL;
	line = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));

	if ( pdbqt_ind == 't' ) {
		sprintf(line,"%-6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%10.3f%s%-4s%s",

		"ATOM", 
		(pdb+index)->atom_ser,
		" ", /* just space */
 		(pdb+index)->pdb_type,
		(pdb+index)->pdb_alt_loc,  
		(pdb+index)->res_type, 
		" ", /* just space */
		(pdb+index)->chain, 
		(pdb+index)->res_ser,		// pdbqt_rank was here
		(pdb+index)->pdb_achar,
		"   ", /* just space */
		(pdb+index)->x_coord,
		(pdb+index)->y_coord,
		(pdb+index)->z_coord,
		(pdb+index)->occ,		
		(pdb+index)->b_factor,
		(pdb+index)->mol2_charge,
		" ", /* just space */
		(pdb+index)->pdbqt_type,	
		"\n");
		}

	else {
		sprintf(line,"%-6s%5d%1s%4s%1s%3s%1s%1s%4d%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%s",

		"ATOM", 
		(pdb+index)->atom_ser,
		" ", /* just space */
 		(pdb+index)->pdb_type,
		(pdb+index)->pdb_alt_loc,  
		(pdb+index)->res_type, 
		" ", /* just space */
		(pdb+index)->chain, 
		(pdb+index)->res_ser,
		(pdb+index)->pdb_achar,
		"   ", /* just space */
		(pdb+index)->x_coord,
		(pdb+index)->y_coord,
		(pdb+index)->z_coord,
		(pdb+index)->occ,
		(pdb+index)->b_factor,
	
		"\n");
		}
	return(line);
	}

void print_plist_waters(re *a, int members, char file_name [MAX_FILENAME_LENGTH], int no_frames, char elitism_switch) {

	FILE *outfile;
	int i;
	int counter = 1;
	double occ_norm_max = 0;
	double occ_norm_min = 0;
	double occ_norm_range = 0;
	double *occ_norm = NULL;

	double max_dbl(double *a, int n);
	double min_dbl(double *a, int n);

	occ_norm = (double *) malloc(members*sizeof(double));

	for (i = 0; i < members; i++) 
		*(occ_norm+i) = (a+i)->avg_water_no/(double)(no_frames+1);

	occ_norm_max = max_dbl(occ_norm,members);
	occ_norm_min = min_dbl(occ_norm,members);

	occ_norm_range = occ_norm_max - occ_norm_min;

	outfile=fopen(file_name,"w");

	switch (elitism_switch) {
		case 'e': fprintf(outfile,"%s","REMARK Prediction list (ID-elitist)");
		break; 
		case 'a': fprintf(outfile,"%s","REMARK Prediction list (ID-all-inclusive)");
		break;
		case 'p': fprintf(outfile,"%s","REMARK Prediction list (Positional)");
		break;
		case 'm': fprintf(outfile,"%s","REMARK Prediction list (Merged)");
		break;
		case 'c': fprintf(outfile,"%s","REMARK Prediction list (Consensus)");
		break;
		case 'l': fprintf(outfile,"%s","REMARK Occupancy-cluster list (Positional)");
		break;
		}

	fprintf(outfile,"%s","\n");

	if ( elitism_switch != 'l' && occ_norm_range != 0) fprintf(outfile,"%s","REMARK B-factor column holds calculated mobility values.\n");
	if ( elitism_switch == 'l' || occ_norm_range == 0 ) fprintf(outfile,"%s","REMARK B-factor column holds normalized occupancy counts.\n");
	  fprintf(outfile,"%s%8d%s","REMARK First frame #             : ",g_frame_range_min_orig,"\n");
	  fprintf(outfile,"%s%8d%s","REMARK Last frame #              : ",g_frame_range_max_orig,"\n");
	fprintf(outfile,"%s%8.3f%s","REMARK Distance limit (A)        : ",g_d_value,"\n");
	if ( elitism_switch != 'l' )
	fprintf(outfile,"%s%8.3f%s","REMARK Prediction tolerance (A)  : ",g_p_tolerance,"\n");
	fprintf(outfile,"%s%8.3f%s","REMARK Clustering tolerance (A)  : ",g_c_tolerance,"\n");
	if ( elitism_switch == 'c' )
	fprintf(outfile,"%s%8.3f%s","REMARK Consensus tolerance (A)   : ",CONS_SCALE*g_p_tolerance,"\n");


	for (i = 0; i < members; i++) {
		fprintf(outfile,"%-6s","HETATM"); 
		fprintf(outfile,"%5d",counter); 
		fprintf(outfile,"%1s",""); 	// just space
		fprintf(outfile,"%-4s","OW"); 	// atom type
		fprintf(outfile,"%1s",""); 	// alt loc
		fprintf(outfile,"%3s","WAT"); 	// res type
		fprintf(outfile,"%1s",""); 	// just space
		fprintf(outfile,"%1s",""); 	// chain ID
		fprintf(outfile,"%4d",counter); 	// res ser = cluster number
		fprintf(outfile,"%1s",""); 	// achar
		fprintf(outfile,"%3s",""); 	// just space
		fprintf(outfile,"%8.3f",(a+i)->x_coord); 
		fprintf(outfile,"%8.3f",(a+i)->y_coord); 
		fprintf(outfile,"%8.3f",(a+i)->z_coord); 
//		fprintf(outfile,"%6s",""); 	// occ = blank
		if ( elitism_switch == 'l' ) fprintf(outfile,"%6s",""); 	// occ = blank
			else fprintf(outfile,"%6.2f", *(occ_norm+i)*100);
		if ( elitism_switch == 'l' )
			fprintf(outfile,"%6.2f",(100*(a+i)->water_no)/(double)(no_frames+1)); 	// for pos olist only
			else { if ( occ_norm_range != 0) fprintf(outfile,"%6.2f",(occ_norm_max-*(occ_norm+i))*100/occ_norm_range); 	// 100 b-factor = occ_norm_max
				if ( occ_norm_range == 0 )fprintf(outfile,"%6.2f", *(occ_norm+i)*100);
				}
		fprintf(outfile,"%s","\n");
		counter++;
		if ( counter == 10000 ) counter = 1;
		}

	free (occ_norm);

	fprintf(outfile,"%s","END\n");
	fclose(outfile);		
	}


void print_double_vector(double *a, int total, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min) {

	FILE *outfile;
	int i;

	outfile=fopen(file_name,"w");

	// Project-specific header
	switch (flag) {
		case 't': fprintf(outfile,"%s","List of C-alpha RMSD values of target molecules of each frame vs. initial target structure.\n");
			  break; 
		case 'l': fprintf(outfile,"%s","List of non-H atom RMSD values of ligand molecules of each frame vs. initial ligand structure.\n");
			  break; 
		case 'r': fprintf(outfile,"%s","List of C-alpha RMSD values of target molecules of each frame vs. reference target structure.\n"); 
			  break;
		case 'd': fprintf(outfile,"%s","List of non-H atom RMSD values of ligand molecules of each frame vs. reference ligand structure.\n"); 
			  break;
		case 's': fprintf(outfile,"%s","Success rate values (%)\n");
			  fprintf(outfile,"%s%8.3f%s","B-factor limit            : ",g_b_factor,"\n");
			  fprintf(outfile,"%s%8.3f%s","Distance limit (A)        : ",g_d_value,"\n");
			  fprintf(outfile,"%s%8.3f%s","Match tolerance (A)       : ",g_s_tolerance,"\n");
			  break;
		case 'o': fprintf(outfile,"%s","Occupancy values (trajectory %)\n");
			  fprintf(outfile,"%s%8.3f%s","B-factor limit            : ",g_b_factor,"\n");
			  fprintf(outfile,"%s%8.3f%s","Distance limit (A)        : ",g_d_value,"\n");
			  fprintf(outfile,"%s%8.3f%s","Match tolerance (A)       : ",g_s_tolerance,"\n");
			  break;
	    	}

	// General part, prints serial number and double value
	for (i = 0; i <= total; i++) {
		fprintf(outfile,"%8d%8.3f%s",frame_min+i,*(a+i),"\n");
		}

	fclose(outfile);
	}


void print_int_vector(int *a, int total, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min) {

	FILE *outfile;
	int i;

	outfile=fopen(file_name,"w");

	// Project-specific header
	switch (flag) {
		case 'i': fprintf(outfile,"%s",TRAJ_NOIFACE_FILE_HEADER);
		break; 
		case 's': fprintf(outfile,"%s",TRAJ_NOSFACE_FILE_HEADER);
		break; 
		case 'x': fprintf(outfile,"%s","REMARK Number of exchanging water molecules for each water molecules in the reference pool.\n");
		break; 
	    	}

	// General part, prints serial number and int value
	for (i = 0; i <= total; i++) {
		fprintf(outfile,"%8d%8d%s",frame_min+i,*(a+i),"\n");
		}

	fclose(outfile);
	}

void print_mlist(mr **ml, int **members, int *clusters, int **ol, int rows, double max_dist, char file_name [MAX_FILENAME_LENGTH]) {

	FILE *outfile;
	int i, j, w;

	outfile=fopen(file_name,"w");

	fprintf(outfile,"%s","Match list between Cluster list and Reference pool\n");
	fprintf(outfile,"%s%8.3f%s","Match tolerance (A): ",max_dist,"\n");

	PRINT_LINE(outfile,68,"=",w)

	fprintf(outfile,"%18s","Cluster#"); 
	fprintf(outfile,"%6s","Count");
	fprintf(outfile,"%3s"," | ");
	fprintf(outfile,"%6s","Ref.#");
	fprintf(outfile,"%10s","Ref.atom#");
	fprintf(outfile,"%10s", "Ref.res.#");
	fprintf(outfile,"%10s", "Distance");
	fprintf(outfile,"%3s"," | ");
	fprintf(outfile,"%s","M");
	fprintf(outfile,"%s","\n");
	
	for (i = 0; i < rows; i++) {
		PRINT_LINE(outfile,68,"-",w)
		fprintf(outfile,"%-10s%-6d%8d%s","Water#: ",i+1,*(*(ol+2)+i),"\n");
		PRINT_LINE(outfile,68,"-",w)
		for (j = 0; j < *(clusters+i); j++) {
			fprintf(outfile,"%12s",""); 
			fprintf(outfile,"%6d",j+1); 
			fprintf(outfile,"%6d",*(*(members+i)+j));
			fprintf(outfile,"%3s"," | ");
			fprintf(outfile,"%6d", (*(ml+i)+j)->ref_water_ser_no);
			fprintf(outfile,"%10d", (*(ml+i)+j)->ref_water_atom_ser_no);
			fprintf(outfile,"%10d", (*(ml+i)+j)->ref_water_res_ser_no);
			fprintf(outfile,"%10.3f", (*(ml+i)+j)->distance);
			fprintf(outfile,"%3s"," | ");
			if ( (*(ml+i)+j)->distance < max_dist )
				fprintf(outfile,"%s","x");
				else fprintf(outfile,"%s","-");
			fprintf(outfile,"%s","\n");
			}
		}

		PRINT_LINE(outfile,68,"=",w)

	fclose(outfile);		
	}


void print_mlist_plist(mr *ml, int **members, re *pl, int rows_pl, double max_dist, double max_b, char file_name [MAX_FILENAME_LENGTH], char elitism_switch) {

	FILE *outfile;
	int i, w;
	int match_counter = 0;

	outfile=fopen(file_name,"w");

	fprintf(outfile,"%s","Match List between ");
	if ( elitism_switch == 'e' )
		fprintf(outfile,"%s","ID-elitist ");
	else if ( elitism_switch == 'a' ) fprintf(outfile,"%s","ID-all-inclusive ");
	else if ( elitism_switch == 'p' ) fprintf(outfile,"%s","Positional ");
	else if ( elitism_switch == 'm' ) fprintf(outfile,"%s","Merged ");
	else if ( elitism_switch == 'c' ) fprintf(outfile,"%s","Consensus ");
	fprintf(outfile,"%s","Prediction List and Reference Pool\n");
	  fprintf(outfile,"%s%8d%s","First frame #             : ",g_frame_range_min_orig,"\n");
	  fprintf(outfile,"%s%8d%s","Last frame #              : ",g_frame_range_max_orig,"\n");
	fprintf(outfile,"%s%8.3f%s","Match tolerance (A)       : ",max_dist,"\n");
	fprintf(outfile,"%s%8.3f%s","B-factor limit            : ",max_b,"\n");
	fprintf(outfile,"%s%8.3f%s","Distance limit (A)        : ",g_d_value,"\n");
	fprintf(outfile,"%s%8.3f%s","Prediction tolerance (A)  : ",g_p_tolerance,"\n");
	fprintf(outfile,"%s%8.3f%s","Clustering tolerance (A)  : ",g_c_tolerance,"\n");
	if ( elitism_switch == 'c' )
	fprintf(outfile,"%s%8.3f%s","Consensus tolerance (A)   : ",CONS_SCALE*g_p_tolerance,"\n");

	PRINT_LINE(outfile,103,"=",w)

	fprintf(outfile,"%10s","Water(pl)#"); 
	if ( elitism_switch == 'm' || elitism_switch == 'c' ) {
		fprintf(outfile,"%15s","Source-pl-type");
		fprintf(outfile,"%15s","Source(pl)#"); 
		}
	else {
		fprintf(outfile,"%15s","Water(ol)#");
		if ( elitism_switch == 'a' || elitism_switch == 'e' ) { 
			fprintf(outfile,"%15s","Cluster(ol)#"); 
			}
		else fprintf(outfile,"%15s","");
		}
	fprintf(outfile,"%10s","Count");
	fprintf(outfile,"%3s"," | ");
	fprintf(outfile,"%6s","Ref.#");
	fprintf(outfile,"%10s","Ref.atom#");
	fprintf(outfile,"%10s", "Ref.res.#");
	fprintf(outfile,"%10s", "B-factor");
	fprintf(outfile,"%10s", "Distance");
	fprintf(outfile,"%3s"," | ");
	fprintf(outfile,"%s","M");
	fprintf(outfile,"%s","\n");

	PRINT_LINE(outfile,103,"-",w)
	
	for (i = 0; i < rows_pl; i++) {
			fprintf(outfile,"%10d",i+1); 
			if ( elitism_switch == 'm' || elitism_switch == 'c' ) {
				fprintf(outfile,"%14s","");
				fprintf(outfile,"%c",(pl+i)->plist_type);
				fprintf(outfile,"%15d",(pl+i)->plist_ser_no); 
				}
			else {
				fprintf(outfile,"%15d",(pl+i)->cluster_no); 
				if ( elitism_switch == 'a' || elitism_switch == 'e' ) fprintf(outfile,"%15d",(pl+i)->cluster2_no); 	
				else fprintf(outfile,"%15s","");
				}
			if ( elitism_switch == 'm' || elitism_switch == 'c' ) fprintf(outfile,"%10.2f",(pl+i)->avg_water_no);
			else fprintf(outfile,"%10d",(pl+i)->water_no);
			fprintf(outfile,"%3s"," | ");
			if ( (ml+i)->ref_water_ser_no != 0 ) {
				fprintf(outfile,"%6d", (ml+i)->ref_water_ser_no);
				fprintf(outfile,"%10d", (ml+i)->ref_water_atom_ser_no);
				fprintf(outfile,"%10d", (ml+i)->ref_water_res_ser_no);
				fprintf(outfile,"%10.2f", (ml+i)->b_factor);			
				fprintf(outfile,"%10.3f", (ml+i)->distance);
				fprintf(outfile,"%3s"," | ");
				}
			else {
				fprintf(outfile,"%6s", "-");
				fprintf(outfile,"%10s", "-");
				fprintf(outfile,"%10s", "-");
				fprintf(outfile,"%10s", "-");			
				fprintf(outfile,"%10s", "-");
				fprintf(outfile,"%3s"," | ");
				}
			if ( (ml+i)->indicator != 0 ) {
				fprintf(outfile,"%s","x");
				match_counter++;
				}
				else fprintf(outfile,"%s","-");
			fprintf(outfile,"%s","\n");
		}

	PRINT_LINE(outfile,103,"=",w)

	fprintf(outfile,"%-5s%10.3f\n","SR",100 * (double) match_counter / g_tolerances->lowb_final_ref_wats );
//	fprintf(outfile,"%-5s%10.3f\n","PR",100 * (double) match_counter / (double) (rows_pl+1) );

	PRINT_LINE(outfile,103,"=",w)

	fclose(outfile);		
	}


void print_olist(int **members, int *clusters, int **ol, int rows, char file_name [MAX_FILENAME_LENGTH], char flag) {

	FILE *outfile;
	int i, j, w;

	outfile=fopen(file_name,"w");

	// Project-specific header
	switch (flag) {
		case 'i': fprintf(outfile,"%s",OLIST_IF_FILE_HEADER);
		break; 
		case 's': fprintf(outfile,"%s",OLIST_SF_FILE_HEADER);
		break; 
	    	}

	PRINT_LINE(outfile,54,"=",w)

	fprintf(outfile,"%18s","Cluster#"); 
	fprintf(outfile,"%10s","Atom#");
	fprintf(outfile,"%10s", "Residue#");
	fprintf(outfile,"%16s","Count");

	fprintf(outfile,"%s","\n");
	
	for (i = 0; i < rows; i++) {
		PRINT_LINE(outfile,54,"-",w)
		fprintf(outfile,"%-10s%-6d%12d%10d%16d%s","Water#: ",i+1,*(*(ol+0)+i),*(*(ol+1)+i),*(*(ol+2)+i),"\n");
		PRINT_LINE(outfile,54,"-",w)
		for (j = 0; j < *(clusters+i); j++) {
			fprintf(outfile,"%12s",""); 
			fprintf(outfile,"%6d",j+1); 
			fprintf(outfile,"%20s",""); 
			fprintf(outfile,"%16d",*(*(members+i)+j));
			fprintf(outfile,"%s","\n");
			}
		}

		PRINT_LINE(outfile,54,"=",w)

	fclose(outfile);		
	}


void print_succ_matr(int ****a, int stol_actual, int low_b_ref_wat, char file_name [MAX_FILENAME_LENGTH]) {

	FILE *outfile;
	int j, k;

	outfile=fopen(file_name,"w");

	fprintf(outfile,"%20s%8.3f%s","REMARK Match tolerance (A): ", (g_tolerances->min_stol + stol_actual*g_tolerances->step_stol), "\n");

	fprintf(outfile,"%15s","REMARK Rows: Const. clustering tolerance | Columns: Const. prediction tolerance\n");

	fprintf(outfile,"%15s","REMARK Success rates (%)\n");

	fprintf(outfile,"%8s","       ");
	for (j = 0; j < g_tolerances->num_ptol; j++) {
		fprintf(outfile,"%6.2f",(g_tolerances->min_ptol + j*g_tolerances->step_ptol));
		fprintf(outfile,"%2s","  ");
		}
	fprintf(outfile,"%s","\n");

	for (j = 0; j < g_tolerances->num_ctol; j++) {
		fprintf(outfile,"%6.2f",(g_tolerances->min_ctol + j*g_tolerances->step_ctol));
		fprintf(outfile,"%2s","  ");
		for (k = 0; k < g_tolerances->num_ptol; k++) {
			fprintf( outfile,"%6.2f", 100 * ( (double) *(*(*(*(a+1)+stol_actual)+k)+j) / (double) low_b_ref_wat ) );
			fprintf(outfile,"%2s","  ");
			} 
		fprintf(outfile,"%s","\n");
			}

/*	fprintf(outfile,"%s","\n");

	fprintf(outfile,"%15s","REMARK Prediction rates (%)\n");

	fprintf(outfile,"%8s","       ");
	for (j = 0; j < g_tolerances->num_ptol; j++) {
		fprintf(outfile,"%6.2f",(g_tolerances->min_ptol + j*g_tolerances->step_ptol));
		fprintf(outfile,"%2s","  ");
		}
	fprintf(outfile,"%s","\n");

	for (j = 0; j < g_tolerances->num_ctol; j++) {
		fprintf(outfile,"%6.2f",(g_tolerances->min_ctol + j*g_tolerances->step_ctol));
		fprintf(outfile,"%2s","  ");
		for (k = 0; k < g_tolerances->num_ptol; k++) {
			fprintf( outfile,"%6.2f", 100 * ( (double) *(*(*(*(a+1)+stol_actual)+k)+j) / (double) *(*(*(*(a+0)+stol_actual)+k)+j) ) );
			fprintf(outfile,"%2s","  ");
			} 
		fprintf(outfile,"%s","\n");
			}
*/

	fclose(outfile);		
	}


void print_traj_cwat_comp(wr **a, int frames, int refs, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min, int *classes) {

	FILE *outfile;
	int i, j;

	int max_int(int *b, int n);

	outfile=fopen(file_name,"w");

	switch (flag) {
		case 'a': fprintf(outfile,"%s","REMARK Oxygen atom serial numbers of closest waters in candidate and reference pools.\n");
		break; 
		case 'r': fprintf(outfile,"%s","REMARK Residue serial numbers of closest waters in candidate and reference pools.\n");
		break; 
		case 'd': fprintf(outfile,"%s","REMARK Distance of oxygen atoms of closest waters in candidate and reference pools.\n");
		break; 
		case '1': fprintf(outfile,"%s","REMARK Oxygen atom serial numbers of exchanging waters.\n");
		break; 
		case '2': fprintf(outfile,"%s","REMARK Residue serial numbers of exchanging waters.\n");
		break; 
		case '3': fprintf(outfile,"%s","REMARK Distance of oxygen atoms between exchanging waters and reference pool waters.\n");
		break; 
		case '4': fprintf(outfile,"%s","REMARK Occurence of exchanging waters (trajectory %).\n");
		break; 
	    	}


	fprintf(outfile,"%15s","");

	if ( flag == 'a' || flag == 'r'|| flag == 'd' ) {
		for (i = 0; i <= frames; i++) fprintf(outfile,"%15d",frame_min+i);
		fprintf(outfile,"%s","\n");
		}

	if ( flag == 'a' || flag == 'r'|| flag == 'd' )
		for (i = 0; i < refs; i++) {
			fprintf(outfile,"%15d",i+1);
			for (j = 0; j <= frames; j++)
				switch (flag) {
					case 'a': fprintf(outfile,"%15d",(*(a+j)+i)->atom_ser); 
					break; 
					case 'r': fprintf(outfile,"%15d",(*(a+j)+i)->res_ser); 
					break; 
					case 'd': fprintf(outfile,"%15.3f",(*(a+j)+i)->dist); 
					break; 
		    			}	
			fprintf(outfile,"%s","\n");
			}


	if ( flag == '1' || flag == '2'|| flag == '3' || flag == '4' ) {
		for (i = 0; i < max_int(classes, refs); i++) fprintf(outfile,"%15d",i+1);
		fprintf(outfile,"%s","\n");
		}

	if ( flag == '1' || flag == '2'|| flag == '3' || flag == '4' )
		for (i = 0; i < refs; i++) {
			fprintf(outfile,"%15d",i+1);
			for (j = 0; j < *(classes+i); j++)
				switch (flag) {
					case '1': fprintf(outfile,"%15d",(*(a+i)+j)->atom_ser); 
					break; 
					case '2': fprintf(outfile,"%15d",(*(a+i)+j)->res_ser); 
					break; 
					case '3': fprintf(outfile,"%15.3f",(*(a+i)+j)->dist); 
					break; 
					case '4': fprintf(outfile,"%15.3f",100*((double)(*(a+i)+j)->counter)/(double)frames); 
					break; 
		    			}	
			fprintf(outfile,"%s","\n");
			}


	fclose(outfile);		
	}

void print_comm_res_wat_file (re *residues, char wat_fname [MAX_FILENAME_LENGTH], int num_dist_res, double dist, double min_dist_tol) {

	FILE *outfile;
	int i;
	char file_out [MAX_FILENAME_LENGTH];
	char *delims = ".";
	char *line = NULL;
	char header [MAX_LINE_LENGTH];
		
	char *print_comm_res_wat_line (re *residues, int index);
	void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );

	sprintf(file_out,"%s%s%s","O_",strtok(wat_fname, delims),"_comm.res\0");
	outfile=fopen(file_out,"w");

	sprintf(header,"%s%6.3f%s","REMARK List of water molecules below ",dist," A distance from both ligand and target.\n");
	fputs(header,outfile);

	sprintf(header,"%s%6.3f%s","REMARK Minimum distance criterion (to avoid close contacts): ", min_dist_tol, " A\n");
	fputs(header,outfile);

	for (i = 0; i < num_dist_res; i++) {
		line = print_comm_res_wat_line (residues, i);
		fputs(line,outfile);
		free(line);				
		}

	fclose(outfile);
	print_log_file ( g_log_file_name, file_out, 'a' ); NLL
	}

char *print_comm_res_wat_line(re *residues, int index) {

	char *line = NULL;
	line = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));

	sprintf(line,"%s%4d%s%3s%1s%1s%s%4d%s",
	"WATER# ",
	index+1,
	"     ",
	(residues+index)->target_res_type, 
	" ", 
	(residues+index)->target_chain, 
	" ", 
	(residues+index)->target_res_ser,
	"\n");
	
	return(line);
	}


void print_dist_file (pd *distances, char ligand_fname [MAX_FILENAME_LENGTH], char target_fname [MAX_FILENAME_LENGTH], int num_distances, char mode_switch1, char mode_switch2, double dist) {

	FILE *outfile;
	int i;
	char file_out [MAX_FILENAME_LENGTH];
	char *line = NULL;
	char header [MAX_LINE_LENGTH];
		
	char *print_dist_line (pd *distances, int index, char mode_switch1, char mode_switch2);

	sprintf(file_out,"%s%s%s%s%s%s%s","O_",g_root_name,"_",ligand_fname,"_",target_fname,".dst\0"); 

	outfile=fopen(file_out,"w");

	if (mode_switch1 == 't' && mode_switch2 == 'l') {
		sprintf(header,"%s%6.3f%s","REMARK List of ligand-target atomic pairs & distances below ",dist," A\n");
		fputs(header,outfile);
		}

	if (mode_switch1 == 't' && mode_switch2 == 'o') {
		sprintf(header,"%s%6.3f%s","REMARK List of docked model-target atomic pairs & distances below ",dist," A\n");
		fputs(header,outfile);
		}

	if (mode_switch1 == 't' && mode_switch2 == 'w') {
		sprintf(header,"%s%6.3f%s","REMARK List of water-target atomic pairs & distances below ",dist," A\n");
		fputs(header,outfile);
		}

	if (mode_switch1 == 'l' && mode_switch2 == 'w') {
		sprintf(header,"%s%6.3f%s","REMARK List of water-ligand atomic pairs & distances below ",dist," A\n");
		fputs(header,outfile);
		}

	if (mode_switch1 == 'w' && mode_switch2 == 'w') {
		sprintf(header,"%s%6.3f%s","REMARK List of water-water atomic pairs & distances below ",dist," A\n");
		fputs(header,outfile);
		}

	if (mode_switch1 == 'b' && mode_switch2 == 'w') {
		sprintf(header,"%s%6.3f%s","REMARK List of water-bulk atomic pairs & distances below ",dist," A\n");
		fputs(header,outfile);
		}

	for (i = 0; i < num_distances; i++) {
		line = print_dist_line(distances, i, mode_switch1, mode_switch2);
		fputs(line,outfile);
		free(line);
		}

	fclose(outfile);
	print_log_file ( g_log_file_name, file_out, 'a' ); NLL

	}

char *print_dist_line(pd *distance, int index, char mode_switch1, char mode_switch2) {

	char *line = NULL;

	char molecule_1 [15];
	char molecule_2 [15];
	char molecule_2v [15];

	char molecule_1_t [] = "   TARGET\0";
	char molecule_1_l [] = "   LIGAND\0";
	char molecule_1_w [] = "    WATER\0";
	char molecule_1_b [] = "     BULK\0";

	char molecule_2_l [] = "LIGAND#\0";
	char molecule_2_o [] = "DOCKED MODEL#\0";
	char molecule_2_w [] = "WATER#\0";

	char molecule_2v_l [] = "LIGAND\0";
	char molecule_2v_o [] = "   DOCKED MODEL\0";
	char molecule_2v_w [] = "WATER\0";


	line = (char *) malloc(MAX_LINE_LENGTH*sizeof(char));

	if (mode_switch1 == 't' && mode_switch2 == 'l') {
		strcpy(molecule_1, molecule_1_t);
		strcpy(molecule_2, molecule_2_l);
		strcpy(molecule_2v, molecule_2v_l);

		}

	if (mode_switch1 == 't' && mode_switch2 == 'o') {
		strcpy(molecule_1, molecule_1_t);
		strcpy(molecule_2, molecule_2_o);
		strcpy(molecule_2v, molecule_2v_o);
		}

	if (mode_switch1 == 't' && mode_switch2 == 'w') {
		strcpy(molecule_1, molecule_1_t);
		strcpy(molecule_2, molecule_2_w);
		strcpy(molecule_2v, molecule_2v_w);
		}
	
	if (mode_switch1 == 'l' && mode_switch2 == 'w') {
		strcpy(molecule_1, molecule_1_l);
		strcpy(molecule_2, molecule_2_w);
		strcpy(molecule_2v, molecule_2v_w);
		}

	if (mode_switch1 == 'w' && mode_switch2 == 'w') {
		strcpy(molecule_1, molecule_1_w);
		strcpy(molecule_2, molecule_2_w);
		strcpy(molecule_2v, molecule_2v_w);
		}

	if (mode_switch1 == 'b' && mode_switch2 == 'w') {
		strcpy(molecule_1, molecule_1_b);
		strcpy(molecule_2, molecule_2_w);
		strcpy(molecule_2v, molecule_2v_w);
		}

	if (mode_switch2 == 'o')

	sprintf(line,"%s%4d%s%5d%1s%4s%1s%3s%1s%1s%4d%1s%s%6d%1s%4s%1s%3s%1s%1s%4d%s%8.3f%s",
	molecule_2,
	(distance+index)->ligand_ser_no,
	molecule_2v, 
	(distance+index)->ligand_atom_ser_no,
	" ", 
 	(distance+index)->ligand_pdb_type,
	" ",  
	(distance+index)->ligand_res_type, 
	" ", 
	(distance+index)->ligand_chain, 
	(distance+index)->ligand_res_ser,
	" ",
	molecule_1,
	(distance+index)->target_atom_ser_no,
	" ", 
 	(distance+index)->target_pdb_type,
	" ",  
	(distance+index)->target_res_type, 
	" ", 
	(distance+index)->target_chain, 
	(distance+index)->target_res_ser,
	"    DISTANCE (A)",
	(distance+index)->target_dist,
	"\n");

	else

	sprintf(line,"%s%5d%1s%4s%1s%3s%1s%1s%4d%1s%s%6d%1s%4s%1s%3s%1s%1s%4d%s%8.3f%s",
	molecule_2v, 
	(distance+index)->ligand_atom_ser_no,
	" ", 
 	(distance+index)->ligand_pdb_type,
	" ",  
	(distance+index)->ligand_res_type, 
	" ", 
	(distance+index)->ligand_chain, 
	(distance+index)->ligand_res_ser,
	" ",
	molecule_1,
	(distance+index)->target_atom_ser_no,
	" ", 
 	(distance+index)->target_pdb_type,
	" ",  
	(distance+index)->target_res_type, 
	" ", 
	(distance+index)->target_chain, 
	(distance+index)->target_res_ser,
	"    DISTANCE (A)",
	(distance+index)->target_dist,
	"\n");

	return(line);

	}



void print_edges(hs *a, int num_edges, char info [MAX_FILENAME_LENGTH/2]) {

	FILE *outfile;
	char file_out [MAX_FILENAME_LENGTH];
	int i;

	sprintf(file_out,"%s%s%s%s","O_",g_root_name,info,".dgs\0"); 
	outfile=fopen(file_out,"w");

	fprintf(outfile,"%s","List of edges in the hydration network\n");

	for (i = 0; i < num_edges; i++) {
		fprintf(outfile,"%-7s","WATER: ");
		fprintf(outfile,"%-6d",(a+i)->wat_atm_ser);
		fprintf(outfile,"%-4s",(a+i)->wat_res_typ);
		fprintf(outfile,"%-5d",(a+i)->wat_res_ser);	
		fprintf(outfile,"%6.2f",(a+i)->wat_b_factor);
		fprintf(outfile,"%-10s"," PARTNER: ");
		fprintf(outfile,"%-6d",(a+i)->ptr_atm_ser);
		fprintf(outfile,"%c",' ');
		fprintf(outfile,"%2s",(a+i)->ptr_chn);
		if ( (a+i)->ptr_type == 't' || (a+i)->ptr_type == 'l' ) fprintf(outfile,"%c",'_'); else fprintf(outfile,"%c",' ');
		fprintf(outfile,"%-4s",(a+i)->ptr_res_typ);
		fprintf(outfile,"%-5d",(a+i)->ptr_res_ser);
		fprintf(outfile,"%6.2f",(a+i)->ptr_b_factor);
		fprintf(outfile,"%c",' ');
		fprintf(outfile,"%c",(a+i)->ptr_type); 
		fprintf(outfile,"%-13s"," DISTANCE: ");	
		fprintf(outfile,"%-5.3f",(a+i)->dist);
		fprintf(outfile,"%c",' ');
		fprintf(outfile,"%c",(a+i)->wat_dyn); 
		fprintf(outfile,"%c",' ');
		fprintf(outfile,"%c",(a+i)->ptr_dyn); 
		fprintf(outfile,"%c",' ');
		fprintf(outfile,"%c",(a+i)->edge_dyn); 
		fprintf(outfile,"%s","\n");
		}

	print_log_file ( g_log_file_name, file_out, 'a' ); NLL

	fclose(outfile);		
	}

void print_nodes(hd *a, int num_nodes, char info [MAX_FILENAME_LENGTH/2]) {

	FILE *outfile;
	char file_out [MAX_FILENAME_LENGTH];
	int i;

	sprintf(file_out,"%s%s%s%s","O_",g_root_name,info,".nds\0"); 
	outfile=fopen(file_out,"w");

	fprintf(outfile,"%s","List of nodes in the hydration network\n");

	for (i = 0; i < num_nodes; i++) {
//		fprintf(outfile,"%-6d",(a+i)->atm_ser);
		fprintf(outfile,"%c",(a+i)->chain[0]);
		if ( (a+i)->type == 't' || (a+i)->type == 'l' ) fprintf(outfile,"%c",'_'); else fprintf(outfile,"%c",' ');
		fprintf(outfile,"%-4s",(a+i)->res_typ);
		fprintf(outfile,"%-7d",(a+i)->res_ser);	
		fprintf(outfile,"%c",' ');
		fprintf(outfile,"%c",(a+i)->type); 
		fprintf(outfile,"%3d",(a+i)->degree);
		if ( (a+i)->type == 'w' ) fprintf(outfile,"%3d",(a+i)->sol_degree); else fprintf(outfile,"%3s","   "); 
		fprintf(outfile,"%8.3f",(a+i)->b_factor);
		fprintf(outfile,"%c",' ');
		fprintf(outfile,"%c",(a+i)->node_dyn); 
		fprintf(outfile,"%s","\n");
		}

	print_log_file ( g_log_file_name, file_out, 'a' ); NLL

	fclose(outfile);		
	}

void print_dot(hd *a, int num_nodes, hs *b, int num_edges, char info [MAX_FILENAME_LENGTH/2]) {

	FILE *outfile;
	char file_out [MAX_FILENAME_LENGTH];
	int i;

	sprintf(file_out,"%s%s%s%s","O_",g_root_name,info,".dot\0"); 
	outfile=fopen(file_out,"w");

	fprintf(outfile,"%s%s","graph ","waternet\n");
	fprintf(outfile,"%s","\n");
	fprintf(outfile,"%s","{\n");
	fprintf(outfile,"%s","\n");

	for (i = 0; i < num_nodes; i++) {
		if ( (a+i)->type == 't' || (a+i)->type == 'l' ) fprintf(outfile,"%c",(a+i)->chain[0]);
		if ( (a+i)->type == 't' || (a+i)->type == 'l' ) fprintf(outfile,"%c",'_'); 
		fprintf(outfile,"%3s",(a+i)->res_typ);
		fprintf(outfile,"%d",(a+i)->res_ser);	
		fprintf(outfile,"%s","[color=\"");
		switch ((a+i)->type) { case 'w': fprintf(outfile,"%s",WAT_ND_COL); break;
				       case 'b': fprintf(outfile,"%s",BLK_ND_COL); break;
				       case 'l': fprintf(outfile,"%s",LIG_ND_COL); break;
				       case 't': fprintf(outfile,"%s",TAR_ND_COL); break;
				       }
		fprintf(outfile,"%s","\",style=");
		fprintf(outfile,"%s",ND_STYLE);
		fprintf(outfile,"%s","]");
		fprintf(outfile,"%s","\n");
		}

	for (i = 0; i < num_edges; i++) {
		fprintf(outfile,"%3s",(b+i)->wat_res_typ);
		fprintf(outfile,"%d",(b+i)->wat_res_ser);	
		fprintf(outfile,"%s"," -- ");
		if ( (b+i)->ptr_type == 't' || (b+i)->ptr_type == 'l' ) fprintf(outfile,"%c",(b+i)->ptr_chn[0]);
		if ( (b+i)->ptr_type == 't' || (b+i)->ptr_type == 'l' ) fprintf(outfile,"%c",'_'); 
		fprintf(outfile,"%3s",(b+i)->ptr_res_typ);
		fprintf(outfile,"%d",(b+i)->ptr_res_ser);
		fprintf(outfile,"%s","[color=\"");
		switch ((b+i)->edge_dyn) { case 'd': fprintf(outfile,"%s",DYN_COL); break;
				          case 's': fprintf(outfile,"%s",STA_COL); break;
 				       }
		fprintf(outfile,"%s","\",style=");
		fprintf(outfile,"%s",ED_STYLE);
		fprintf(outfile,"%s","]");
		fprintf(outfile,"%s","\n");
		}
	fprintf(outfile,"%s","\n");
	fprintf(outfile,"%s","}\n");

	print_log_file ( g_log_file_name, file_out, 'a' ); NLL

	fclose(outfile);		
	}

void print_gdf(hd *a, int num_nodes, hs *b, int num_edges, char info [MAX_FILENAME_LENGTH/2]) {

	FILE *outfile;
	char file_out [MAX_FILENAME_LENGTH];
	int i;

	sprintf(file_out,"%s%s%s%s","O_",g_root_name,info,".gdf\0"); 
	outfile=fopen(file_out,"w");


	fprintf(outfile,"%s","nodedef>name VARCHAR,color VARCHAR,style VARCHAR,mobility DOUBLE\n");

	for (i = 0; i < num_nodes; i++) {
		if ( (a+i)->type == 't' || (a+i)->type == 'l' ) fprintf(outfile,"%c",(a+i)->chain[0]);
		if ( (a+i)->type == 't' || (a+i)->type == 'l' ) fprintf(outfile,"%c",'_'); 
		fprintf(outfile,"%3s",(a+i)->res_typ);
		fprintf(outfile,"%d",(a+i)->res_ser);	
		fprintf(outfile,"%s",",");
		switch ((a+i)->type) { case 'w': fprintf(outfile,"%s",WAT_ND_COL_RGB); break;
				       case 'b': fprintf(outfile,"%s",BLK_ND_COL_RGB); break;
				       case 'l': fprintf(outfile,"%s",LIG_ND_COL_RGB); break;
				       case 't': fprintf(outfile,"%s",TAR_ND_COL_RGB); break;
				       }
		fprintf(outfile,"%s",",filled,");
		fprintf(outfile,"%f",(a+i)->b_factor);
		fprintf(outfile,"%s","\n");
		}

	fprintf(outfile,"%s","edgedef>node1 VARCHAR,node2 VARCHAR,style VARCHAR,color VARCHAR\n");

	for (i = 0; i < num_edges; i++) {
		fprintf(outfile,"%3s",(b+i)->wat_res_typ);
		fprintf(outfile,"%d",(b+i)->wat_res_ser);	
		fprintf(outfile,"%s",",");
		if ( (b+i)->ptr_type == 't' || (b+i)->ptr_type == 'l' ) fprintf(outfile,"%c",(b+i)->ptr_chn[0]);
		if ( (b+i)->ptr_type == 't' || (b+i)->ptr_type == 'l' ) fprintf(outfile,"%c",'_'); 
		fprintf(outfile,"%3s",(b+i)->ptr_res_typ);
		fprintf(outfile,"%d",(b+i)->ptr_res_ser);
		fprintf(outfile,"%s",",solid,");
		switch ((b+i)->edge_dyn) { case 'd': fprintf(outfile,"%s",DYN_COL_RGB); break;
				          case 's': fprintf(outfile,"%s",STA_COL_RGB); break;
 				       }
		fprintf(outfile,"%s","\n");
		}

//	fprintf(outfile,"%s","\n");

	print_log_file ( g_log_file_name, file_out, 'a' ); NLL

	fclose(outfile);		
	}

