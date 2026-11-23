// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "analysis.h"

void analysis (re *reference_waters, int num_reference_waters, double no_lowb_ref_wats, int *total_if) {
	
	// types
	FILE *infile;
	char file_in [MAX_FILENAME_LENGTH];
	char out_file_name [MAX_FILENAME_LENGTH];
	int i, j, l, min_serial, actual_xchange_class;
	int match_counter = 0;
	int offset = 0;
	int frames = g_frame_range_max - g_frame_range_min;		// +1 gives real frame numbers
	int xchange_counter = 0;
	int occupancy_counter = 0;
	double *distances = NULL;
	double *rmsd_value = NULL;
	double *success_rate = NULL;
	int *no_xchange_wats = NULL;
	double *occupancy = NULL;
	pl *waters = NULL;
	wr **matrices = NULL;
	wr **mobility = NULL;
	wr *x = NULL;

	str0 (file_in, MAX_FILENAME_LENGTH);
	str0 (out_file_name, MAX_FILENAME_LENGTH);
	
	// Initialize other output storage
	matrices = (wr **) malloc((frames+1) * sizeof(wr *)); 
	for (i = 0; i < frames+1; i++) 
		*(matrices+i) = (wr *) malloc(num_reference_waters * sizeof(wr));
	mobility = (wr **) malloc (num_reference_waters * sizeof (wr*));

	if ( g_v == 'v' || g_v == 'd') rmsd_value = (double *) malloc((frames+1) * sizeof(double));
	success_rate = (double *) malloc((frames+1) * sizeof(double));
	no_xchange_wats = (int *) malloc(num_reference_waters * sizeof(int));
	occupancy = (double *) malloc(num_reference_waters * sizeof(double));

	// Open binary pool
	sprintf(file_in,"%s%s%s","O_",g_root_name,".plw\0");	
	infile = fopen(file_in, "rb");
	if ( infile == NULL ) {
		printf("%s%s%s","ERROR The ",file_in," file could not be opened.\n");
		exit (1);
		}
	frames++; 		// Now, it is real frame number	

	// Calculating target_ref-target_frame calpha rmsds & Selecting closest waters in each frame 
	for ( i = 0; i < frames; i++ ) {
		match_counter = 0;
		if ( *(total_if+i) != 0 ) {
			waters = (pl *) malloc(*(total_if+i) * sizeof(pl));
			offset = 0;
			for (j = 0; j < i; j++) offset = offset + *(total_if+j) * sizeof(pl);
			fseek(infile, offset, SEEK_SET);			
			fread(waters, sizeof(pl), *(total_if+i), infile);
			for (j = 0; j < num_reference_waters; j++) {
				distances = (double *) malloc(*(total_if+i) * sizeof(double));
				for (l = 0; l < *(total_if+i); l++)
					*(distances+l) = dist_calc((reference_waters+j)->x_coord, (reference_waters+j)->y_coord, (reference_waters+j)->z_coord, (waters+l)->x_coord, (waters+l)->y_coord, (waters+l)->z_coord, "O", "O", ' ');
				min_serial = min_dbl_reno(distances, *(total_if+i));
				(*(matrices+i)+j)->atom_ser = (waters+min_serial)->aser;
				(*(matrices+i)+j)->res_ser = (waters+min_serial)->rser;
				(*(matrices+i)+j)->dist = *(distances+min_serial);
				if ( *(distances+min_serial) < g_s_tolerance && match_int ((waters+min_serial)->aser, (waters+min_serial)->rser, *(matrices+i), j) == 0  && (reference_waters+j)->b_factor < g_b_factor ) 
					match_counter++;		// to count the same match only once and use only references below b-factor limit
				free(distances);
				distances = NULL;
				}
			free(waters);
			}
		*(success_rate+i) = 100 * (double) match_counter / no_lowb_ref_wats;
		}

	fclose (infile);

	// Evaluation of mobility 
	for (j = 0; j < num_reference_waters; j++) {
		xchange_counter = 0;
		for ( i = 0; i < frames; i++ ) {
			actual_xchange_class = 0;
			if ( i == 0 ) {
				*(mobility+j) = (wr *) malloc(sizeof(wr));
				(*(mobility+j)+0)->atom_ser = (*(matrices+i)+j)->atom_ser;
				(*(mobility+j)+0)->res_ser  = (*(matrices+i)+j)->res_ser;
				(*(mobility+j)+0)->dist     = (*(matrices+i)+j)->dist;
				(*(mobility+j)+0)->counter = 1;
				xchange_counter = 1;
				}
				else {
					actual_xchange_class = match_int_2 ( (*(matrices+i)+j)->atom_ser, (*(matrices+i)+j)->res_ser, *(mobility+j), xchange_counter );
					if ( actual_xchange_class != 0 )
						(*(mobility+j)+actual_xchange_class-1)->counter++;					
					}
			if ( i != 0 && actual_xchange_class == 0 ) {
				xchange_counter++;
				x = realloc( *(mobility+j), xchange_counter * sizeof(wr) );
				if (x) *(mobility+j) = x;
				x = NULL;
				(*(mobility+j)+xchange_counter-1)->atom_ser = (*(matrices+i)+j)->atom_ser;
				(*(mobility+j)+xchange_counter-1)->res_ser  = (*(matrices+i)+j)->res_ser;
				(*(mobility+j)+xchange_counter-1)->dist  = (*(matrices+i)+j)->dist;
				(*(mobility+j)+xchange_counter-1)->counter = 1;
				}
			}
		*(no_xchange_wats+j) = xchange_counter;
		}


	for (j = 0; j < num_reference_waters; j++)
		insertion_sort_wr (*(mobility+j), *(no_xchange_wats+j));

	// Evaluation of occupancy
	for (j = 0; j < num_reference_waters; j++) {
		occupancy_counter = 0;
		for ( i = 0; i < frames; i++ ) 
			if ( (*(matrices+i)+j)->dist < g_s_tolerance ) occupancy_counter++;
		*(occupancy+j) = 100 * ((double) occupancy_counter / (double) frames ) ;
		}


	// Print occupancy file
	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_occy.txt\0");
	print_double_vector(occupancy, num_reference_waters-1, out_file_name, 'o', 1);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	free(occupancy);

	// Print success rate file
	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_succ.txt\0");
	print_double_vector(success_rate, frames-1, out_file_name, 's', g_frame_range_min_orig);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	free(success_rate);

	// Print 3 result matrices into files
	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_aser.txt\0");
	print_traj_cwat_comp(matrices, frames-1, num_reference_waters, out_file_name, 'a', g_frame_range_min_orig, no_xchange_wats);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_rser.txt\0");
	print_traj_cwat_comp(matrices, frames-1, num_reference_waters, out_file_name, 'r', g_frame_range_min_orig, no_xchange_wats);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_dist.txt\0");
	print_traj_cwat_comp(matrices, frames-1, num_reference_waters, out_file_name, 'd', g_frame_range_min_orig, no_xchange_wats);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	free(matrices);

	//Print number of EX waters
	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_NoEX.txt\0");
	print_int_vector(no_xchange_wats, num_reference_waters-1, out_file_name, 'x', 1);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	// Print 4 mobility statistics matrices into files
	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_asEX.txt\0");
	print_traj_cwat_comp(mobility, frames-1, num_reference_waters, out_file_name, '1', g_frame_range_min_orig, no_xchange_wats);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_rsEX.txt\0");
	print_traj_cwat_comp(mobility, frames-1, num_reference_waters, out_file_name, '2', g_frame_range_min_orig, no_xchange_wats);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_diEX.txt\0");
	print_traj_cwat_comp(mobility, frames-1, num_reference_waters, out_file_name, '3', g_frame_range_min_orig, no_xchange_wats);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_ocEX.txt\0");
	print_traj_cwat_comp(mobility, frames-1, num_reference_waters, out_file_name, '4', g_frame_range_min_orig, no_xchange_wats);
	print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	free(mobility);
	free(no_xchange_wats);

	}


int match_int (int atom_no, int ser_no, wr *actual_wr, int actual_j) {
	int i = 0;
	int label = 0;
	for ( i = 0; i < (actual_j-1); i++ ) 
		if ( atom_no == (actual_wr+i)->atom_ser && ser_no == (actual_wr+i)->res_ser && (actual_wr+i)->dist < g_s_tolerance ) label = 1;
	return label;
	}

int match_int_2 (int atom_no, int ser_no, wr *actual_wr, int actual_upper) {
	int i = 0;
	int label = 0;
	for ( i = 0; i < actual_upper; i++ ) 
		if ( atom_no == (actual_wr+i)->atom_ser && ser_no == (actual_wr+i)->res_ser ) label = i+1;
	return label;
	}


