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
extern int g_frame_range_min;
extern char g_xtc_ind;

double dist_calc (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch);
void insertion_sort_int_re_ss (int *a, re **b, int num);
void mean_coords(re ***a, int **members, int *clusters, int rows, re **mean);
void print_all_clust_waters(re ***a, int **members, int *clusters, int rows, char file_name [MAX_FILENAME_LENGTH], char i_label);
void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );
void str0 (char string [], int length);

