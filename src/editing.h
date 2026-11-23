// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 14-06-2015

#include "global.h"

extern char g_root_name [MAX_FILENAME_LENGTH];
extern char g_temp [MAX_FILENAME_LENGTH];
extern char g_log_file_name [MAX_FILENAME_LENGTH];
extern char g_v;

void str0 (char string [], int length);
ap *read_in_pdb (char *file_name, int *atom_num, int mdl_serno);
ap *box_to_pdb (ap *box_pdb, char limits [MAX_FILENAME_LENGTH], int *atom_num_pdb, int atom_num_box, char label);
char *print_pdb_file (ap *pdb, int atom_num, char file_out [MAX_FILENAME_LENGTH]);
void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );
ap *merge_mols (ap *mol1, ap *mol2, int atom_num1, int atom_num2, int *atom_numo);
pd *pdb_ligand_dist_full_min (ap *pdb_target, ap *pdb_ligand, int pdb_num_at_target, int pdb_num_at_ligand, int *dist_num, double dist_tolerance, char h_switch);
re *extract_common_wat_res_min (pd *distances_t,  pd *distances_l, int num_dist_t, int num_dist_l, int *num_common_wat, char min_dist_switch, double min_dist_tolerance);
void print_comm_res_wat_file (re *residues, char wat_fname [MAX_FILENAME_LENGTH], int num_dist_res, double dist, double min_dist_tol);
void print_dist_file (pd *distances, char ligand_fname [MAX_FILENAME_LENGTH], char target_fname [MAX_FILENAME_LENGTH], int num_distances, char mode_switch1, char mode_switch2, double dist);


