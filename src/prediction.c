// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "prediction.h"

void prediction (int *no_is_waters, re *reference_waters, char root_name [MAX_FILENAME_LENGTH], double dist_tol, double s_tol, double p_tol, double c_tol, double b_factor, int num_reference_waters, double no_lowb_ref_wats, char reference, char selector) {

	char file_in [MAX_FILENAME_LENGTH];
	char out_file_name [MAX_FILENAME_LENGTH];
	int i;
	int row_no_occ_list = 0;
	int row_no_occ_list_pos = 0;	
	int ****success_all = NULL;		// FREE
	int ****success_elit = NULL;		// FREE
	int ****success_pos = NULL;		// FREE
	int ****success_merged = NULL;		// FREE
	int *no_clusters = NULL;
	int **occ_list = NULL;
	int **no_cl_members = NULL;
	int no_members_plist_elit = 0;
	int no_members_plist_all = 0;
	int no_members_plist_pos = 0;
	int no_members_plist_merged_a_e = 0;
	int no_members_plist_merged_all = 0;
	int no_frames = g_frame_range_max - g_frame_range_min;		// +1 gives real no_frames
	re *plist_elit = NULL;
	re *plist_all = NULL;
	re **mn_coords = NULL;
	re *pos_occ_list = NULL;
	re *plist_pos = NULL;
	re *plist_merged_a_e = NULL;
	re *plist_merged_all = NULL;
	mr ** match_list = NULL;		// FREE
	mr *match_plist_elit = NULL;		// FREE
	mr *match_plist_all = NULL;		// FREE
	mr *match_plist_pos = NULL;
	mr *match_plist_merged_all = NULL;

	time_t clock = time(NULL);

	str0 (file_in, MAX_FILENAME_LENGTH);
	str0 (out_file_name, MAX_FILENAME_LENGTH);

	// Occupancy list
	if ( g_c_scheme == 'a' || g_c_scheme == 'e' || g_c_scheme == 'm' ) { occ_list = make_occ_list_bin(no_frames, no_is_waters, &row_no_occ_list, root_name);
		printf("\nID occupancy list done:\n%s", ctime(&clock)); }
	if ( g_c_scheme == 'p' || g_c_scheme == 'm' ) { pos_occ_list = make_occ_list_pos_bin(no_frames, no_is_waters, &row_no_occ_list_pos, c_tol, root_name);
		printf("\nPOS occupancy list done:\n%s", ctime(&clock)); }

	// Time for distances
	clock = time(NULL); printf("\nAll occupancies done:\n%s", ctime(&clock));

	// Clustering (default c_tol value)
	if ( g_c_scheme == 'a' || g_c_scheme == 'e' || g_c_scheme == 'm' ) {
		no_clusters = (int *) malloc(row_no_occ_list * sizeof(int));
		no_cl_members = (int **) malloc(row_no_occ_list * sizeof(int *));
		for (i = 0; i < row_no_occ_list; i++) *(no_cl_members+i) = (int *) malloc(sizeof (int));

		if ( selector == 'i' )
			mn_coords = cluster_by_olist_bin (occ_list, row_no_occ_list, no_is_waters, no_frames, c_tol, no_clusters, no_cl_members, root_name, g_v, 'i');
			else // surface
				mn_coords = cluster_by_olist_bin (occ_list, row_no_occ_list, no_is_waters, no_frames, c_tol, no_clusters, no_cl_members, root_name, g_v, 's');
		}

	// Prediction (default p_tol value)
	if ( g_c_scheme == 'e' || g_c_scheme == 'm' ) plist_elit = make_plist(mn_coords, no_clusters, no_cl_members, row_no_occ_list, p_tol, 'e', &no_members_plist_elit);
	if ( g_c_scheme == 'a' || g_c_scheme == 'm' ) plist_all = make_plist(mn_coords, no_clusters, no_cl_members, row_no_occ_list, p_tol, 'a', &no_members_plist_all);
	if ( g_c_scheme == 'p' || g_c_scheme == 'm' ) plist_pos = make_plist_pos(pos_occ_list, row_no_occ_list_pos, p_tol, &no_members_plist_pos);
	if ( g_c_scheme == 'm' ) plist_merged_a_e = merge_plists(plist_elit, no_members_plist_elit, plist_all, no_members_plist_all, &no_members_plist_merged_a_e, p_tol);
	if ( g_c_scheme == 'm' ) plist_merged_all = merge_plists(plist_merged_a_e, no_members_plist_merged_a_e, plist_pos, no_members_plist_pos, &no_members_plist_merged_all, p_tol);

	// Match lists & success rates  (default & scaled *_tol values)
	if ( reference == 'r') {
		if ( g_c_scheme == 'a' || g_c_scheme == 'e' || g_c_scheme == 'm' ) match_list = match_clusters(mn_coords, no_clusters, row_no_occ_list, reference_waters, num_reference_waters);
		if ( g_c_scheme == 'e' || g_c_scheme == 'm' ) match_plist_elit = match_plist(plist_elit, no_members_plist_elit, reference_waters, num_reference_waters, b_factor, s_tol);		
		if ( g_c_scheme == 'a' || g_c_scheme == 'm' ) match_plist_all = match_plist(plist_all, no_members_plist_all, reference_waters, num_reference_waters, b_factor, s_tol);
		if ( g_c_scheme == 'p' || g_c_scheme == 'm' ) match_plist_pos = match_plist(plist_pos, no_members_plist_pos, reference_waters, num_reference_waters, b_factor, s_tol);		
		if ( g_c_scheme == 'm' ) match_plist_merged_all = match_plist(plist_merged_all, no_members_plist_merged_all, reference_waters, num_reference_waters, b_factor, s_tol);		

		if ( g_tolerances->num_stol != 0 ) {
			time_t clock = time(NULL); printf("\nMatrix calculations started:\n%s", ctime(&clock));
			if ( g_c_scheme == 'a' || g_c_scheme == 'm' ) success_all  = match_eval (dist_tol, 'a', row_no_occ_list, occ_list, no_is_waters, no_frames, root_name, reference_waters, num_reference_waters, b_factor, no_lowb_ref_wats);
			if ( g_c_scheme == 'e' || g_c_scheme == 'm' ) success_elit = match_eval (dist_tol, 'e', row_no_occ_list, occ_list, no_is_waters, no_frames, root_name, reference_waters, num_reference_waters, b_factor, no_lowb_ref_wats);
			if ( g_c_scheme == 'p' || g_c_scheme == 'm' ) success_pos = match_eval (dist_tol, 'p', row_no_occ_list, occ_list, no_is_waters, no_frames, root_name, reference_waters, num_reference_waters, b_factor, no_lowb_ref_wats);
			if ( g_c_scheme == 'm' ) success_merged = match_eval (dist_tol, 'm', row_no_occ_list, occ_list, no_is_waters, no_frames, root_name, reference_waters, num_reference_waters, b_factor, no_lowb_ref_wats);
			}
		}

	// Printing
	if ( g_v == 'v'  || g_v == 'd' ) {		
		if ( g_c_scheme == 'a' || g_c_scheme == 'e' || g_c_scheme == 'm' ) {
			// Print occupancy list
			str0 (out_file_name, MAX_FILENAME_LENGTH);
			sprintf(out_file_name,"%s%s%s","O_",root_name,"_ocID.lst\0");
			print_olist(no_cl_members, no_clusters, occ_list, row_no_occ_list, out_file_name, 'i');
			print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
		
			// Print cluster list
			str0 (out_file_name, MAX_FILENAME_LENGTH);
			sprintf(out_file_name,"%s%s%s","O_",root_name,"_clID.pdb\0");
			print_clist_waters(mn_coords, no_cl_members, no_clusters, occ_list, row_no_occ_list, out_file_name, 'i');
			print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
			}
		if ( g_c_scheme == 'p' || g_c_scheme == 'm' )		
			// Print occupancy-cluster list of positional clustering
			PRINT_PLIST("_clPO.pdb",pos_occ_list,row_no_occ_list_pos,'l')
		}

	// Print elitist plist
	if ( g_c_scheme == 'e' || g_c_scheme == 'm' ) {
		PRINT_PLIST("_prIDe.pdb",plist_elit,no_members_plist_elit,'e')
		PRINT_PLIST("_prIDe_top.pdb",plist_elit,(int) no_members_plist_elit*g_p_top*0.01,'e')
		}

	// Print all-inclusive plist
	if ( g_c_scheme == 'a' || g_c_scheme == 'm' ) {
		PRINT_PLIST("_prIDa.pdb",plist_all,no_members_plist_all,'a')
		PRINT_PLIST("_prIDa_top.pdb",plist_all,(int) no_members_plist_all*g_p_top*0.01,'e')
		}

	// Print positional plist
	if ( g_c_scheme == 'p' || g_c_scheme == 'm' ) {
		PRINT_PLIST("_prPOS.pdb",plist_pos,no_members_plist_pos,'p')
		PRINT_PLIST("_prPOS_top.pdb",plist_pos,(int) no_members_plist_pos*g_p_top*0.01,'e')
		}

	// Print merged plist
	if ( g_c_scheme == 'm' ) {
		PRINT_PLIST("_prMER.pdb",plist_merged_all,no_members_plist_merged_all,'m')
		PRINT_PLIST("_prMER_top.pdb",plist_merged_all,(int) no_members_plist_merged_all*g_p_top*0.01,'e')
		}

	if ( reference == 'r') {
		if ( ( g_v == 'v'  || g_v == 'd' ) && ( g_c_scheme == 'a' || g_c_scheme == 'e' || g_c_scheme == 'm' ) ) {
			// Print cluster match list
			str0 (out_file_name, MAX_FILENAME_LENGTH);
			sprintf(out_file_name,"%s%s%s","O_",root_name,"_mtIDc.lst\0");   
			print_mlist(match_list, no_cl_members, no_clusters, occ_list, row_no_occ_list, s_tol, out_file_name);
			print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
			}

		// Print plist match list (Elitist)
		if ( g_c_scheme == 'e' || g_c_scheme == 'm' ) PRINT_MLIST_PLIST("_mtIDe.lst",match_plist_elit,plist_elit,no_members_plist_elit,'e')

		// Print plist match list (All-inclusive)
		if ( g_c_scheme == 'a' || g_c_scheme == 'm' ) PRINT_MLIST_PLIST("_mtIDa.lst",match_plist_all,plist_all,no_members_plist_all,'a')

		// Print plist match list (Positional)
		if ( g_c_scheme == 'p' || g_c_scheme == 'm' ) PRINT_MLIST_PLIST("_mtPOS.lst",match_plist_pos,plist_pos,no_members_plist_pos,'p')

		// Print plist match list (Merged all)
		if ( g_c_scheme == 'm' ) PRINT_MLIST_PLIST("_mtMER.lst",match_plist_merged_all,plist_merged_all,no_members_plist_merged_all,'m')

		// Print success rate matrices
		if ( g_tolerances->num_stol != 0 )
			for (i = 0; i < g_tolerances->num_stol; i++) {
				if ( g_c_scheme == 'a' || g_c_scheme == 'm' ) PRINT_SUCC_MATR("_IDa.mat",success_all)
				if ( g_c_scheme == 'e' || g_c_scheme == 'm' ) PRINT_SUCC_MATR("_IDe.mat",success_elit)
				if ( g_c_scheme == 'p' || g_c_scheme == 'm' ) PRINT_SUCC_MATR("_POS.mat",success_pos)
				if ( g_c_scheme == 'm' ) PRINT_SUCC_MATR("_MER.mat",success_merged)
				}
		}
	
	if ( g_c_scheme == 'a' || g_c_scheme == 'e' || g_c_scheme == 'm' ) {
		for (i = 0; i < 3; i++) free(*(occ_list+i));
		free(occ_list);
		}
	if ( g_c_scheme == 'p' || g_c_scheme == 'm' )
		free(pos_occ_list);

	if ( g_c_scheme == 'a' || g_c_scheme == 'e' || g_c_scheme == 'm' ) {
		for (i = 0; i < row_no_occ_list; i++) free(*(mn_coords+i));
		free(mn_coords);
		}

	if ( g_c_scheme == 'a' || g_c_scheme == 'e' || g_c_scheme == 'm' ) {
		free(no_clusters);
		for (i = 0; i < row_no_occ_list; i++) free(*(no_cl_members+i));
		free(no_cl_members);
		}

	}

