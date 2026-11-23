// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 10-04-2016

#include "netdraw.h"

void netdraw (char *file_in, char t_range [MAX_FILENAME_LENGTH], char l_range [MAX_FILENAME_LENGTH], char w_range [MAX_FILENAME_LENGTH], double min_dist_tol, double max_dist_tol) {

	int num_at_box, num_at_target, num_at_ligand, num_at_waters, num_at_bulk, num_at_solute, dist_num_wat_lig, dist_num_wat_tar, dist_num_wat_wat, dist_num_wat_blk, num_at_box_orig;
	char w_range_new [MAX_FILENAME_LENGTH];
	int edges_num_wat_tar, edges_num_wat_tar_flt, edges_num_wat_lig, edges_num_wat_lig_flt, edges_num_wat_wat, edges_num_wat_blk, edges_num_1, edges_num_2, edges_num, edges_num_wo_bulk;
	int edges_num_flt4, edges_num_wo_bulk_flt4;
	int nodes_num, nodes_num_wo_bulk;
	ap *pdb_box = NULL;
	ap *pdb_solute = NULL;
	ap *pdb_target = NULL;		
	ap *pdb_ligand = NULL;
	ap *pdb_waters = NULL;
	ap *pdb_bulk = NULL;
	ap *pdb_box_orig = NULL;

	pd *distances_wat_tar = NULL;
	pd *distances_wat_lig = NULL;
	pd *distances_wat_wat = NULL;
	pd *distances_wat_blk = NULL;

	hd *nodes = NULL;
	hd *nodes_wo_bulk = NULL;

	hs *edges_wat_tar = NULL;
	hs *edges_wat_tar_flt = NULL;
	hs *edges_wat_lig = NULL;
	hs *edges_wat_lig_flt = NULL;
	hs *edges_wat_wat = NULL;
	hs *edges_wat_blk = NULL;
	hs *edges_1 = NULL;
	hs *edges_2 = NULL;
	hs *edges = NULL;
	hs *edges_flt4 = NULL;
	hs *edges_wo_bulk = NULL;
	hs *edges_wo_bulk_flt4 = NULL;

	pdb_box = read_in_pdb (file_in, &num_at_box, -1);

	printf ("\nNetworking...\n");

	if ( strncmp(w_range,AUTO_KEYW, 4) == 0 ) strncpy(w_range_new,WAT_KEYW, 3);
		else strcpy(w_range_new,w_range);

	pdb_target = box_to_pdb (pdb_box, t_range, &num_at_target, num_at_box, 't');
	if ( g_l == 'l' ) pdb_ligand = box_to_pdb (pdb_box, l_range, &num_at_ligand, num_at_box, 'l');
	pdb_waters = box_to_pdb (pdb_box, w_range_new, &num_at_waters, num_at_box, 'w');
	pdb_bulk   = box_to_pdb (pdb_box, w_range_new, &num_at_bulk, num_at_box, 'b');

	distances_wat_tar = pdb_ligand_dist_full_min (pdb_target, pdb_waters, num_at_target, num_at_waters, &dist_num_wat_tar, max_dist_tol, ' ');
	if ( g_l == 'l' ) distances_wat_lig = pdb_ligand_dist_full_min (pdb_ligand, pdb_waters, num_at_ligand, num_at_waters, &dist_num_wat_lig, max_dist_tol, ' ');
	distances_wat_wat = pdb_ligand_dist_full_min (pdb_waters, pdb_waters, num_at_waters, num_at_waters, &dist_num_wat_wat, max_dist_tol, ' ');
	distances_wat_blk = pdb_ligand_dist_full_min (pdb_bulk, pdb_waters, num_at_bulk, num_at_waters, &dist_num_wat_blk, max_dist_tol, ' ');

	print_dist_file (distances_wat_tar, "water", "target", dist_num_wat_tar, 't', 'w', max_dist_tol);
	if ( g_l == 'l' ) print_dist_file (distances_wat_lig, "water", "ligand", dist_num_wat_lig, 'l', 'w', max_dist_tol);
	print_dist_file (distances_wat_wat, "water", "water", dist_num_wat_wat, 'w', 'w', max_dist_tol);
	print_dist_file (distances_wat_blk, "water", "bulk", dist_num_wat_blk, 'b', 'w', max_dist_tol);


	edges_wat_tar = edge_lister (distances_wat_tar, dist_num_wat_tar, &edges_num_wat_tar, 't');
	edges_wat_tar_flt = edge_filter_C (edges_wat_tar, edges_num_wat_tar, &edges_num_wat_tar_flt);
	if ( g_l == 'l' ) edges_wat_lig = edge_lister (distances_wat_lig, dist_num_wat_lig, &edges_num_wat_lig, 'l');
	if ( g_l == 'l' ) edges_wat_lig_flt = edge_filter_C (edges_wat_lig, edges_num_wat_lig, &edges_num_wat_lig_flt);
	edges_wat_wat = edge_lister (distances_wat_wat, dist_num_wat_wat, &edges_num_wat_wat, 'w');
	edges_wat_blk = edge_lister (distances_wat_blk, dist_num_wat_blk, &edges_num_wat_blk, 'b');

	if ( g_l == 'l' ) edges_1 = merge_edges (edges_wat_tar_flt, edges_wat_lig_flt, edges_num_wat_tar_flt, edges_num_wat_lig_flt, &edges_num_1);
	edges_2 = merge_edges (edges_wat_wat, edges_wat_blk, edges_num_wat_wat, edges_num_wat_blk, &edges_num_2);
	if ( g_l == 'l' ) edges = merge_edges (edges_1, edges_2, edges_num_1, edges_num_2, &edges_num);
		else edges = merge_edges (edges_wat_tar_flt, edges_2, edges_num_wat_tar_flt, edges_num_2, &edges_num);
	if ( g_l == 'l' ) edges_wo_bulk = merge_edges (edges_1, edges_wat_wat, edges_num_1, edges_num_wat_wat, &edges_num_wo_bulk);
		else edges_wo_bulk = merge_edges (edges_wat_tar_flt, edges_wat_wat, edges_num_wat_tar_flt, edges_num_wat_wat, &edges_num_wo_bulk);
	
	// All
	insertion_sort_hs_wat_res_ser (edges, edges_num);
	edges_flt4 = edge_filter_4 (edges, edges_num, &edges_num_flt4);
	nodes = list_nodes (edges_flt4, edges_num_flt4, &nodes_num);	
	count_degree (edges_flt4, edges_num_flt4, nodes, nodes_num);
	set_static (edges_flt4, edges_num_flt4, nodes, nodes_num);
	insertion_sort_hd_res_ser (nodes, nodes_num);
	insertion_sort_hd_type (nodes, nodes_num);
	print_nodes(nodes, nodes_num, "_all");
	print_edges(edges_flt4, edges_num_flt4, "_all");
	print_dot(nodes, nodes_num, edges_flt4, edges_num_flt4, "_all");
	print_gdf(nodes, nodes_num, edges_flt4, edges_num_flt4, "_all");

	// Without bulk
	insertion_sort_hs_wat_res_ser (edges_wo_bulk, edges_num_wo_bulk);
//	print_edges(edges_wo_bulk, edges_num_wo_bulk, "_wo_bulk_flt4_b");
	edges_wo_bulk_flt4 = edge_filter_4 (edges_wo_bulk, edges_num_wo_bulk, &edges_num_wo_bulk_flt4);
	nodes_wo_bulk = list_nodes (edges_wo_bulk_flt4, edges_num_wo_bulk_flt4, &nodes_num_wo_bulk);
	count_degree (edges_wo_bulk_flt4, edges_num_wo_bulk_flt4, nodes_wo_bulk, nodes_num_wo_bulk);
	set_static (edges_wo_bulk_flt4, edges_num_wo_bulk_flt4, nodes_wo_bulk, nodes_num_wo_bulk);
	insertion_sort_hd_res_ser (nodes_wo_bulk, nodes_num_wo_bulk);
	insertion_sort_hd_type (nodes_wo_bulk, nodes_num_wo_bulk);
	print_edges(edges_wo_bulk_flt4, edges_num_wo_bulk_flt4, "_wo_bulk");
	print_nodes(nodes_wo_bulk, nodes_num_wo_bulk, "_wo_bulk");
	print_dot(nodes_wo_bulk, nodes_num_wo_bulk, edges_wo_bulk_flt4, edges_num_wo_bulk_flt4, "_wo_bulk");
	print_gdf(nodes_wo_bulk, nodes_num_wo_bulk, edges_wo_bulk_flt4, edges_num_wo_bulk_flt4, "_wo_bulk");
	
	if ( g_l == 'l' ) pdb_solute  = merge_mols (pdb_target, pdb_ligand, num_at_target, num_at_ligand, &num_at_solute);

	if ( g_l == 'l' ) pdb_box_orig = merge_mols (pdb_solute, pdb_waters, num_at_solute, num_at_waters, &num_at_box_orig);
		else pdb_box_orig = merge_mols (pdb_target, pdb_waters, num_at_target, num_at_waters, &num_at_box_orig);

	PRINT_PDB_FILES ("_org.pdb",pdb_box_orig,num_at_box_orig)


	// Logfile time
	NLL
	print_log_file ( g_log_file_name, "Atom counts\n", 'a' );
	if ( g_l == 'l' ) PAR2LOGFILE_IC("Ligand all                    ",num_at_ligand)
	PAR2LOGFILE_IC("Target all                    ",num_at_target)
	if ( g_l == 'l' ) PAR2LOGFILE_IC("Water -- in interface         ",num_at_waters)
	if ( g_l != 'l' ) PAR2LOGFILE_IC("Water -- on surface           ",num_at_waters)
	PAR2LOGFILE_IC("Bulk                          ",num_at_bulk)
	print_log_file ( g_log_file_name, "Distance counts\n", 'a' );
	PAR2LOGFILE_IC("Water - target                ",dist_num_wat_tar)
	if ( g_l == 'l' ) PAR2LOGFILE_IC("Water - ligand                ",dist_num_wat_lig)
	PAR2LOGFILE_IC("Water - water                 ",dist_num_wat_wat)
	PAR2LOGFILE_IC("Water - bulk                  ",dist_num_wat_blk)
	print_log_file ( g_log_file_name, "Edge counts before filtering\n", 'a' );
	PAR2LOGFILE_IC("Water - target                ",edges_num_wat_tar)
	if ( g_l == 'l' ) PAR2LOGFILE_IC("Water - ligand                ",edges_num_wat_lig)
	PAR2LOGFILE_IC("Water - water                 ",edges_num_wat_wat)
	PAR2LOGFILE_IC("Water - bulk                  ",edges_num_wat_blk)
	print_log_file ( g_log_file_name, "Edge counts after C-filtering\n", 'a' );
	PAR2LOGFILE_IC("Water - target                ",edges_num_wat_tar_flt)
	if ( g_l == 'l' ) PAR2LOGFILE_IC("Water - ligand                ",edges_num_wat_lig_flt)
	PAR2LOGFILE_IC("Water - all                   ",edges_num)
	PAR2LOGFILE_IC("Water - all    -- without bulk",edges_num_wo_bulk)
	print_log_file ( g_log_file_name, "Final edge counts after 4-filtering\n", 'a' );
	PAR2LOGFILE_IC("Water - all                   ",edges_num_flt4)
	PAR2LOGFILE_IC("Water - all    -- without bulk",edges_num_wo_bulk_flt4)
	print_log_file ( g_log_file_name, "Final node counts\n", 'a' );
	PAR2LOGFILE_IC("Linked to all                 ",nodes_num)
	PAR2LOGFILE_IC("Linked to all  -- without bulk",nodes_num_wo_bulk)
	}
