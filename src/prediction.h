// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "global.h"

extern char g_log_file_name [MAX_FILENAME_LENGTH];
extern ts *g_tolerances;
extern int g_frame_range_min;
extern int g_frame_range_max;
extern ir *g_input_back;
extern char g_v;
extern char g_c_scheme;
extern double g_p_top;

pd *pdb_ligand_dist_full (ap *pdb_target, ap *pdb_ligand, int pdb_num_at_target, int pdb_num_at_ligand, int *dist_num, double dist_tolerance);
int min_dbl_reno(double *a, int n);
void print_double_vector(double *a, int total, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min);
void print_int_vector(int *a, int total, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min);
int **make_occ_list(re **a, int frames, int *total_if, int *no_rows_ol);
void print_olist(int **members, int *clusters, int **ol, int rows, char file_name [MAX_FILENAME_LENGTH], char flag);
re **cluster_by_olist (int **ol, int no_rows_ol, re **ifw, int *total_ifw, int frames, double clust_tol, int *num_cl, int **num_cl_member, char box_nm [MAX_FILENAME_LENGTH], char printer, char i_label);
void print_clist_waters(re **a, int **members, int *clusters, int **ol, int rows, char file_name [MAX_FILENAME_LENGTH], char flag);
mr **match_clusters(re **mean, int *clusters, int rows, re *ref_waters, int num_ref_waters);
void print_mlist(mr **ml, int **members, int *clusters, int **ol, int rows, double max_dist, char file_name [MAX_FILENAME_LENGTH]);
re *make_plist(re **mean, int *clusters, int **members, int rows, double identical_tolerance, char switch_elit, int *plist_no);
void print_plist_waters(re *a, int members, char file_name [MAX_FILENAME_LENGTH], int no_frames, char elitism_switch);
mr *match_plist(re *predict, int rows, re *ref_waters, int num_ref_waters, double b_fact_threshold, double s_value);
void print_mlist_plist(mr *ml, int **members, re *pl, int rows_pl, double max_dist, double max_b, char file_name [MAX_FILENAME_LENGTH], char elitism_switch);
void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );
int ****match_eval (double distance_tol, char eval_selector, int row_no_occ_list, int **occ_list, int *no_if_waters, int no_frames, char box_name [MAX_FILENAME_LENGTH], re *reference_waters, int num_reference_waters, double b_factor, double no_lowb_ref_wats);
void print_succ_matr(int ****a, int stol_actual, int low_b_ref_wat, char file_name [MAX_FILENAME_LENGTH]);
re *make_plist_pos(re *mean, int rows, double identical_tolerance, int *plist_no);
re *merge_plists(re *p1, int no1, re *p2, int no2, int *merged_no, double p_tol);
re *consensus_plists(re *p1, int no1, re *p2, int no2, int *cons_no, double p_tol);
void str0 (char string [], int length);
int **make_occ_list_bin(int frames, int *total_if, int *no_rows_ol, char box_name [MAX_FILENAME_LENGTH]);
re *make_occ_list_pos_bin(int frames, int *total_if, int *no_rows_ol, int c_tolerance, char box_name [MAX_FILENAME_LENGTH]);
re **cluster_by_olist_bin (int **ol, int no_rows_ol, int *total_ifw, int frames, double clust_tol, int *num_cl, int **num_cl_member, char box_name [MAX_FILENAME_LENGTH], char printer, char i_label);

