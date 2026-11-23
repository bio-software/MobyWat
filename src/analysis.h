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
extern double g_b_factor;
extern double g_s_tolerance;
extern int g_frame_range_max_orig;
extern int g_frame_range_min_orig;
extern char g_root_name [MAX_FILENAME_LENGTH];

double dist_calc (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch);
void print_traj_cwat_comp(wr **a, int frames, int refs, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min, int *classes);
int match_int (int atom_no, int ser_no, wr *actual_wr, int actual_j);
int match_int_2 (int atom_no, int ser_no, wr *actual_wr, int actual_upper);
void insertion_sort_wr (wr *a, int num);
void str0 (char string [], int length);
ap *read_in_pdb (char *file_name, int *atom_num, int mdl_serno);
ap *box_to_pdb (ap *box_pdb, char limits [MAX_FILENAME_LENGTH], int *atom_num_pdb, int atom_num_box, char label);
double rmsd_calc (ap *reference, ap *actual, int atom_number, char a_type_sw);
int min_dbl_reno(double *a, int n);
void print_double_vector(double *a, int total, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min);
void print_int_vector(int *a, int total, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min);
void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );

