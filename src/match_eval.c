// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "prediction.h"

int ****match_eval (double distance_tol, char eval_selector, int row_no_occ_list, int **occ_list, int *no_if_waters, int no_frames, char box_name [MAX_FILENAME_LENGTH], re *reference_waters, int num_reference_waters, double b_factor, double no_lowb_ref_wats) {

	int i, j, k, l;	
	int no_members_plist_all = 0;
	int no_members_plist_elit = 0;
	int no_members_plist_pos = 0;
	int no_members_plist_merged_a_e = 0;
	int no_members_plist_merged_all = 0;
	int no_members_plist_con_a_e = 0;
	int no_members_plist_con_all = 0;
	int row_no_occ_list_pos = 0;

	int *no_clusters = NULL;
	int **no_cl_members = NULL;

	int ****match_numbers = NULL;

	re **mean_coords = NULL;
	re *plist_all = NULL;
	re *plist_elit = NULL;
	re *plist_pos = NULL;
	re *plist_merged_a_e = NULL;	
	re *plist_merged_all = NULL;
	re *pos_occ_list = NULL; 
	re *plist_con_a_e = NULL;
	re *plist_con_all = NULL;

	mr *match_plist_result = NULL;

	int calc_success(mr *ml, int rows);

	// Reserve memory for output
	match_numbers = (int ****) malloc(2 * sizeof(int ***));		// 0 = all plist members, 1 = matched plist members
	for (k = 0; k < 2; k++)
		*(match_numbers+k) = (int ***) malloc((g_tolerances->num_stol + 1) * sizeof (int **));
	for (k = 0; k < 2; k++)
		for (i = 0; i <= g_tolerances->num_stol; i++)
			*(*(match_numbers+k)+i) = (int **) malloc((g_tolerances->num_ptol + 1) * sizeof (int *));			
	for (k = 0; k < 2; k++)
		for (i = 0; i <= g_tolerances->num_stol; i++)
			for (l = 0; l <= g_tolerances->num_ptol; l++)
				*(*(*(match_numbers+k)+i)+l) = (int *) malloc((g_tolerances->num_ctol + 1) * sizeof (int));		

	// Calculation of success rates stored in match_numbers
	for (j = 0; j <= g_tolerances->num_stol; j++)	
		for (k = 0; k <= g_tolerances->num_ptol; k++) 
			for (l = 0; l <= g_tolerances->num_ctol; l++) { 
	
				// Clustering
				if ( eval_selector == 'a' || eval_selector == 'e' || eval_selector == 'm' || eval_selector == 'c' ) {
					no_clusters = (int *) malloc(row_no_occ_list * sizeof(int));
					no_cl_members = (int **) malloc(row_no_occ_list * sizeof(int *));
					for (i = 0; i < row_no_occ_list; i++) *(no_cl_members+i) = (int *) malloc(sizeof (int));
					mean_coords = cluster_by_olist_bin(occ_list, row_no_occ_list, no_if_waters, no_frames, (g_tolerances->min_ctol + l*g_tolerances->step_ctol), no_clusters, no_cl_members, box_name, ' ', ' ');
					}
				if ( eval_selector == 'p' || eval_selector == 'm' || eval_selector == 'c' )
					pos_occ_list = make_occ_list_pos_bin(no_frames, no_if_waters, &row_no_occ_list_pos, (g_tolerances->min_ctol + l*g_tolerances->step_ctol), box_name);

				// Prediction list
				if ( eval_selector == 'a' || eval_selector == 'm' || eval_selector == 'c' )
					plist_all = make_plist(mean_coords, no_clusters, no_cl_members, row_no_occ_list, (g_tolerances->min_ptol + k*g_tolerances->step_ptol), 'a', &no_members_plist_all);
				if ( eval_selector == 'e' || eval_selector == 'm' || eval_selector == 'c' )
					plist_elit = make_plist(mean_coords, no_clusters, no_cl_members, row_no_occ_list, (g_tolerances->min_ptol + k*g_tolerances->step_ptol), 'e', &no_members_plist_elit);
				if ( eval_selector == 'p' || eval_selector == 'm' || eval_selector == 'c' )
					plist_pos = make_plist_pos(pos_occ_list, row_no_occ_list_pos, (g_tolerances->min_ptol + k*g_tolerances->step_ptol), &no_members_plist_pos);
				if ( eval_selector == 'm' ) {
					plist_merged_a_e = merge_plists(plist_elit, no_members_plist_elit, plist_all, no_members_plist_all, &no_members_plist_merged_a_e, (g_tolerances->min_ptol + k*g_tolerances->step_ptol));
					plist_merged_all = merge_plists(plist_merged_a_e, no_members_plist_merged_a_e, plist_pos, no_members_plist_pos, &no_members_plist_merged_all, (g_tolerances->min_ptol + k*g_tolerances->step_ptol));
					}
				if ( eval_selector == 'c' ) {
					plist_con_a_e = consensus_plists(plist_elit, no_members_plist_elit, plist_all, no_members_plist_all, &no_members_plist_con_a_e, (g_tolerances->min_ptol + k*g_tolerances->step_ptol));
					plist_con_all = consensus_plists(plist_con_a_e, no_members_plist_con_a_e, plist_pos, no_members_plist_pos, &no_members_plist_con_all, (g_tolerances->min_ptol + k*g_tolerances->step_ptol));
					}

				// Match list & success 
				if ( eval_selector == 'a' ) {
					match_plist_result = match_plist(plist_all, no_members_plist_all, reference_waters, num_reference_waters, b_factor, (g_tolerances->min_stol + j*g_tolerances->step_stol));		
					*(*(*(*(match_numbers+0)+j)+k)+l) = no_members_plist_all;
					*(*(*(*(match_numbers+1)+j)+k)+l) = calc_success(match_plist_result, no_members_plist_all);		
					}
				if ( eval_selector == 'e' ) {
					match_plist_result = match_plist(plist_elit, no_members_plist_elit, reference_waters, num_reference_waters, b_factor, (g_tolerances->min_stol + j*g_tolerances->step_stol));		
					*(*(*(*(match_numbers+0)+j)+k)+l) = no_members_plist_elit;
					*(*(*(*(match_numbers+1)+j)+k)+l) = calc_success(match_plist_result, no_members_plist_elit);		
					}
				if ( eval_selector == 'p' ) {
					match_plist_result = match_plist(plist_pos, no_members_plist_pos, reference_waters, num_reference_waters, b_factor, (g_tolerances->min_stol + j*g_tolerances->step_stol));
					*(*(*(*(match_numbers+0)+j)+k)+l) = no_members_plist_pos;
					*(*(*(*(match_numbers+1)+j)+k)+l) = calc_success(match_plist_result, no_members_plist_pos);		
					}
				if ( eval_selector == 'm' ) {
					match_plist_result = match_plist(plist_merged_all, no_members_plist_merged_all, reference_waters, num_reference_waters, b_factor, (g_tolerances->min_stol + j*g_tolerances->step_stol));
					*(*(*(*(match_numbers+0)+j)+k)+l) = no_members_plist_merged_all;
					*(*(*(*(match_numbers+1)+j)+k)+l) = calc_success(match_plist_result, no_members_plist_merged_all);								
					}
				if ( eval_selector == 'c' ) {
					match_plist_result = match_plist(plist_con_all, no_members_plist_con_all, reference_waters, num_reference_waters, b_factor, (g_tolerances->min_stol + j*g_tolerances->step_stol));
					*(*(*(*(match_numbers+0)+j)+k)+l) = no_members_plist_con_all;
					*(*(*(*(match_numbers+1)+j)+k)+l) = calc_success(match_plist_result, no_members_plist_con_all);								
					}


				// Free memory
				if ( eval_selector == 'a' || eval_selector == 'e' || eval_selector == 'm' || eval_selector == 'c' ) {
					for (i = 0; i < row_no_occ_list; i++) free(*(mean_coords+i));
					free(mean_coords);	
					free(no_clusters);
					for (i = 0; i < row_no_occ_list; i++) free(*(no_cl_members+i));
					free(no_cl_members);
					}
				if ( eval_selector == 'a' || eval_selector == 'm' || eval_selector == 'c' )
					free(plist_all);
				if ( eval_selector == 'e' || eval_selector == 'm' || eval_selector == 'c' )
					free(plist_elit);
				if ( eval_selector == 'p' || eval_selector == 'm' || eval_selector == 'c' ) 
					free(plist_pos);
				if ( eval_selector == 'm' ) {
					free(plist_merged_a_e);
					free(plist_merged_all);
					}
				if ( eval_selector == 'c' ) {
					free(plist_con_a_e);
					free(plist_con_all);
					}
				if ( eval_selector == 'p' || eval_selector == 'm' || eval_selector == 'c' )
					free(pos_occ_list);
				free(match_plist_result);
				}

	return (match_numbers);

	}

