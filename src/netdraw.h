// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 10-04-2016

#include "global.h"

extern char g_root_name [MAX_FILENAME_LENGTH];
extern char g_temp [MAX_FILENAME_LENGTH];
extern char g_log_file_name [MAX_FILENAME_LENGTH];
extern char g_v;
extern char g_l;

void str0 (char string [], int length);
ap *read_in_pdb (char *file_name, int *atom_num, int mdl_serno);
ap *box_to_pdb (ap *box_pdb, char limits [MAX_FILENAME_LENGTH], int *atom_num_pdb, int atom_num_box, char label);
char *print_pdb_file (ap *pdb, int atom_num, char file_out [MAX_FILENAME_LENGTH]);
void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );
ap *merge_mols (ap *mol1, ap *mol2, int atom_num1, int atom_num2, int *atom_numo);
pd *pdb_ligand_dist_full_min (ap *pdb_target, ap *pdb_ligand, int pdb_num_at_target, int pdb_num_at_ligand, int *dist_num, double dist_tolerance, char h_switch);
// re *extract_common_wat_res_min (pd *distances_t,  pd *distances_l, int num_dist_t, int num_dist_l, int *num_common_wat, char min_dist_switch, double min_dist_tolerance);
// void print_comm_res_wat_file (re *residues, char wat_fname [MAX_FILENAME_LENGTH], int num_dist_res, double dist, double min_dist_tol);
void print_dist_file (pd *distances, char ligand_fname [MAX_FILENAME_LENGTH], char target_fname [MAX_FILENAME_LENGTH], int num_distances, char mode_switch1, char mode_switch2, double dist);
hs *edge_lister (pd *distances, int num_distances, int *list_length, char ptr_type);
hs *merge_edges (hs *mol1, hs *mol2, int atom_num1, int atom_num2, int *atom_numo);
hs *edge_filter_C (hs *edge_list, int num_edges, int *filt_num_edges_C);
// hd *list_nodes (hs *edge_list, int num_edges, int *num_nodes);
void insertion_sort_hs_wat_res_ser (hs *a, int num);
hs *edge_filter_4 (hs *edge_list, int num_edges, int *filt_num_edges_4);
void print_edges(hs *a, int num_edges, char info [MAX_FILENAME_LENGTH/2]);
hd *list_nodes (hs *edge_list, int num_edges, int *num_nodes);
void print_nodes(hd *a, int num_nodes, char info [MAX_FILENAME_LENGTH/2]);
void insertion_sort_hd_res_ser (hd *a, int num);
void insertion_sort_hd_type (hd *a, int num);
void print_dot(hd *a, int num_nodes, hs *b, int num_edges, char info [MAX_FILENAME_LENGTH/2]);
void count_degree (hs *edge_list, int num_edges, hd *node_list, int num_nodes);
void set_static (hs *edge_list, int num_edges, hd *node_list, int num_nodes);
void print_gdf(hd *a, int num_nodes, hs *b, int num_edges, char info [MAX_FILENAME_LENGTH/2]);
