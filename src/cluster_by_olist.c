// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "cluster_by_olist.h"

re **cluster_by_olist_bin (int **ol, int no_rows_ol, int *total_ifw, int frames, double clust_tol, int *num_cl, int **num_cl_member, char box_name [MAX_FILENAME_LENGTH], char printer, char i_label) {

	FILE *infile;
	int i, j, k, l, m, counter_cluster, label = 0, label_main = 0;
	double distance;
	char file_in [MAX_FILENAME_LENGTH];
	char out_file_name [MAX_FILENAME_LENGTH];
	int offset = 0;
	int *w = NULL;
	int **label_mind_dist = NULL;
	pl *waters = NULL;
	re ***clusters = NULL;
	re *y = NULL;
	re **z = NULL;
	re **m_coords = NULL;

	frames++; 		// Now, it is real frame number	

	str0 (out_file_name, MAX_FILENAME_LENGTH);

	label_mind_dist = (int **) malloc(frames * sizeof(int *)); 
	for (i = 0; i < frames; i++)
		*(label_mind_dist+i) = (int *) malloc(*(total_ifw+i) * sizeof(int)); 

	if ( g_xtc_ind != 'w' ) sprintf(file_in,"%s%s%s","O_",box_name,".plw\0");	
		else sprintf(file_in,"%s%s",box_name,".plw\0");
	infile = fopen(file_in, "rb");
	if ( infile == NULL ) {
		printf("%s%s%s","ERROR The ",file_in," file could not be opened.\n");
		exit (1);
		}
	
	// Zeros for non-clustered waters
	for (i = 0; i < frames; i++) 
		for (j = 0; j < *(total_ifw+i); j++)
			*(*(label_mind_dist+i)+j) = 0;

	// Allocate memory for output for o_list rows
	clusters = (re ***) malloc(no_rows_ol * sizeof(re **)); 

	// Clustering 1: initial fill-up, first member of first clusters (malloc)
	for (i = 0; i < no_rows_ol; i++) {
		label = 0;
		*(clusters+i) = (re **) malloc(sizeof (re *));			// Memory for one cluster per o_list row
		*(*(clusters+i)+0) = (re *) malloc(sizeof (re));		// Memory for one=first re element per cluster
		*(*(num_cl_member+i)+0) = 0;					// No cluster members in the beginning.
		for (j = 0; j < frames; j++) 
			if ( *(total_ifw+j) != 0 ) {
				waters = (pl *) malloc(*(total_ifw+j) * sizeof(pl));
				offset = 0;
				for (m = 0; m < j; m++) offset = offset + *(total_ifw+m) * sizeof(pl);
				fseek(infile, offset, SEEK_SET);			
				fread(waters, sizeof(pl), *(total_ifw+j), infile);
				for (k = 0; k < *(total_ifw+j); k++)
					if (label == 0 && *(*(ol+0)+i) == (waters+k)->aser && *(*(ol+1)+i) == (waters+k)->rser ) {
						label = 1;					// Only the first
						*(*(label_mind_dist+j)+k) = 1;			// Mark as clustered (= 1)
						(*(*(clusters+i)+0)+0)->x_coord = (waters+k)->x_coord;	// Copy re's
						(*(*(clusters+i)+0)+0)->y_coord = (waters+k)->y_coord;	// Copy re's
						(*(*(clusters+i)+0)+0)->z_coord = (waters+k)->z_coord;	// Copy re's
						(*(*(clusters+i)+0)+0)->target_atom_ser_no = (waters+k)->aser;	// Copy re's
						(*(*(clusters+i)+0)+0)->target_res_ser = (waters+k)->rser;	// Copy re's
						(*(*(clusters+i)+0)+0)->ligand_ser_no = g_frame_range_min + j; // will go into B-factor of clist
						(*(*(num_cl_member+i)+0))++;			// Increase cluster member number
						}	
				free(waters);
				waters = NULL;
				}
		*(num_cl+i) = 1;
		}


	// Clustering 2: add aditional clusters & elements (realloc) & first elements of additional clusters (malloc)
	for (i = 0; i < no_rows_ol; i++) {
		counter_cluster = 1;
		for (j = 0; j < frames; j++) {
			if ( *(total_ifw+j) != 0 ) {
				waters = (pl *) malloc(*(total_ifw+j) * sizeof(pl));
				offset = 0;
				for (m = 0; m < j; m++) offset = offset + *(total_ifw+m) * sizeof(pl);
				fseek(infile, offset, SEEK_SET);			
				fread(waters, sizeof(pl), *(total_ifw+j), infile);
				for (k = 0; k < *(total_ifw+j); k++)							// i,j,k: check if re belongs to o_list row & if not already used in a cluster
					if ( *(*(label_mind_dist+j)+k) == 0 && *(*(ol+0)+i) == (waters+k)->aser && *(*(ol+1)+i) == (waters+k)->rser ) {
						label_main = 1;			// Suppose that it is not member (=1) of any clusters in the row
						for (l = 0; l < counter_cluster; l++) {					// l,m: check if member of existing cluster or new cluster should be opened
							label = 0;							// Suppose that it is member of cluster l
							for (m = 0; m < *(*(num_cl_member+i)+l) ; m++) {	// Checks for all member of cluster l
								distance = dist_calc( (*(*(clusters+i)+l)+m)->x_coord, (*(*(clusters+i)+l)+m)->y_coord, (*(*(clusters+i)+l)+m)->z_coord, (waters+k)->x_coord, (waters+k)->y_coord, (waters+k)->z_coord, "0", "0", 'q');
								if ( distance > clust_tol ) label = 1;
								}
							if (label == 0 && label_main != 0) {							// If it belongs to cluster l, then realloc cluster size + add new member
								(*(*(num_cl_member+i)+l))++;			// Increase cluster member number
								y = realloc( *(*(clusters+i)+l), (*(*(num_cl_member+i)+l)) * sizeof (re) );
								if (y) *(*(clusters+i)+l) = y;
								y = NULL;
								(*(*(clusters+i)+l)+(*(*(num_cl_member+i)+l))-1)->x_coord = (waters+k)->x_coord;	// Copy re's
								(*(*(clusters+i)+l)+(*(*(num_cl_member+i)+l))-1)->y_coord = (waters+k)->y_coord;	// Copy re's
								(*(*(clusters+i)+l)+(*(*(num_cl_member+i)+l))-1)->z_coord = (waters+k)->z_coord;	// Copy re's
								(*(*(clusters+i)+l)+(*(*(num_cl_member+i)+l))-1)->target_atom_ser_no = (waters+k)->aser;	// Copy re's
								(*(*(clusters+i)+l)+(*(*(num_cl_member+i)+l))-1)->target_res_ser = (waters+k)->rser;	// Copy re's
								(*(*(clusters+i)+l)+(*(*(num_cl_member+i)+l))-1)->ligand_ser_no = g_frame_range_min + j; // will go into B-factor of clist
								*(*(label_mind_dist+j)+k) = 1;
								label_main = 0;
								}
							}					
						if (label_main == 1) {
							counter_cluster++;	
							z = realloc( *(clusters+i), counter_cluster * sizeof (re *) );	// Memory for new cluster
							if (z) *(clusters+i) = z;
							z = NULL;
							w = realloc ( *(num_cl_member+i), counter_cluster * sizeof (int) );	// One more int* cluster-member-holder per row
							if (w) *(num_cl_member+i) = w;
							w = NULL;
							*(*(clusters+i)+counter_cluster-1) = (re *) malloc(sizeof (re));		// Memory for 1st element of new cluster
							*(*(num_cl_member+i)+counter_cluster-1) = 1;					// Set cluster member number to 1
							(*(*(clusters+i)+counter_cluster-1)+0)->x_coord = (waters+k)->x_coord;
							(*(*(clusters+i)+counter_cluster-1)+0)->y_coord = (waters+k)->y_coord;
							(*(*(clusters+i)+counter_cluster-1)+0)->z_coord = (waters+k)->z_coord;
							(*(*(clusters+i)+counter_cluster-1)+0)->target_atom_ser_no = (waters+k)->aser;
							(*(*(clusters+i)+counter_cluster-1)+0)->target_res_ser = (waters+k)->rser;
							(*(*(clusters+i)+counter_cluster-1)+0)->ligand_ser_no = g_frame_range_min + j; // will go into B-factor of clist
							*(*(label_mind_dist+j)+k) = 1;
							}		
						}
				free(waters);
				waters = NULL;
				}
			}
		*(num_cl+i) = counter_cluster;
		}


	fclose (infile);

	// Sort clusters (by descending member number)
	for (j = 0; j < no_rows_ol; j++)
		insertion_sort_int_re_ss (*(num_cl_member+j), *(clusters+j), *(num_cl+j));

	// Allocate memory for mean coords
	m_coords = (re **) malloc(no_rows_ol * sizeof(re *));
	for (j = 0; j < no_rows_ol; j++)
		*(m_coords+j) = (re *) malloc(*(num_cl+j) * sizeof (re));

	// Calculate mean coords
	mean_coords(clusters, num_cl_member, num_cl, no_rows_ol, m_coords);

	if (printer == 'd') {
		// Print all clustered IF waters
		if (i_label == 'i') sprintf(out_file_name,"%s%s%s","O_",box_name,"_cIFw.pdb\0");
			else sprintf(out_file_name,"%s%s%s","O_",box_name,"_cSFw.pdb\0");
		print_all_clust_waters(clusters, num_cl_member, num_cl, no_rows_ol, out_file_name, i_label);	
		print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
		}

	for (i = 0; i < frames; i++)
		free(*(label_mind_dist+i));
	free (label_mind_dist);

	for(i = 0; i < no_rows_ol; i++)
		for(j = 0; j < *(num_cl+i); j++)
			free(*(*(clusters+i)+j));
	for(i = 0; i < no_rows_ol; i++) free(*(clusters+i));
	free(clusters);

	return(m_coords);
	}



