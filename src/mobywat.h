// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 03-12-2015

#include "global.h"

FILE *g_outfile;

time_t t0, t1;
time_t ti; 

char g_m = DEF_M;
char g_e = DEF_E;
char g_v = DEF_V_MODE;
char g_r = ' ';
char g_f = ' ';		
char g_l = ' ';		
char g_w = ' ';		
char g_t = ' ';		
char g_n = ' ';		// Used for specifying NMR=n or separate=s PDB file input types by user
char g_mdl_ind = ' ';
char g_xtc_ind = ' ';
char g_tpy = ' ';
char g_c_scheme = 'a';
char g_pli = ' ';

double g_c_tolerance = DEF_CLUST_TOLERANCE;
double g_d_value = DEF_DIST_TOLERANCE;
double g_d_min_value = DEF_MIN_DIST_TOLERANCE;
double g_p_tolerance = DEF_IDENTICAL_PREDICTION_TOLERANCE;
double g_p_top = DEF_P_TOP;
char g_f_name [MAX_FILENAME_LENGTH] = DEF_FRAME_FILENAME;
char g_ref_name [MAX_FILENAME_LENGTH] = DEF_MOBYREF_FILENAME;
char g_root_name [MAX_FILENAME_LENGTH] = "";
char g_l_range [MAX_FILENAME_LENGTH] = ""; 
char g_m_mode [MAX_FILENAME_LENGTH] = DEF_PROG_MODE;
char g_n_frames [MAX_FILENAME_LENGTH] = DEF_FRAME_RANGE;
double g_s_tolerance = DEF_SEP_TOLERANCE; // = Match tolerance
double g_m_tolerance = DEF_MOB_TOLERANCE; // = Mobility tolerance
char g_t_range [MAX_FILENAME_LENGTH] = "";
char g_w_range [MAX_FILENAME_LENGTH] = "";
char g_pli_name [MAX_FILENAME_LENGTH] = "";

tp *g_topology = NULL;
ap *g_pdb_ref = NULL;
ap *g_ref_ligand = NULL;
ap *g_ref_target = NULL;
ap *g_ref_waters = NULL;

int g_pdb_no = 0;
int g_pdb_ref_no = 0;
int g_ref_ligand_no = 0;
int g_ref_target_no = 0;
int g_ref_waters_no = 0;
int g_ref_target_ca_no = 0;
int g_mdl_count = 0;
int g_offset = 0;
int g_frame_max_written_into_plw = 0;

char g_temp [MAX_FILENAME_LENGTH];
char g_log_file_name [MAX_FILENAME_LENGTH];
char g_top_file_name [MAX_FILENAME_LENGTH] = DEF_TOPOL_FILENAME;

re *g_final_ref_waters = NULL;
int g_final_ref_waters_no = 0;
char *g_delims = "-";
char *g_delims2 = "_.";
char *g_delims3 = ".";

char ***g_ref_ranges = NULL;
int g_frame_range_max = 10;
int g_frame_range_min = 0;
int g_frame_range_max_orig = 0;
int g_frame_range_min_orig = 0;

double g_b_factor = DEF_MAX_B_FACTOR;
double g_lowb_final_ref_wats = 0;

ts *g_tolerances;
int *g_ref_count_no = NULL;
//int *g_e_wish = NULL;		// Not used!

ir *g_input_back;

int *g_no_is_waters = NULL;

pf g_pool_info;

void get_args(int argc, char *argv[], char *f_value, char *n_value, char *l_value, char *t_value, char *w_value, double *b_factor, char *e_mode, double *c_tolerance, double *d_tolerance, char *f_name, char *l_range, char *m_mode, char *n_frames, double *s_tolerance, char *t_range, char *w_range, double *p_tolerance, char *v_mode, double *d_min_tolerance, double *m_tolerance);
char find_token_in_name (char name[MAX_FILENAME_LENGTH], char token [MAX_FILENAME_LENGTH], char *delims);
ap *read_in_pdb (char *file_name, int *atom_num, int mdl_serno);
re *select_if_waters (ap *pdb_target, ap *pdb_ligand, ap *pdb_waters, int pdb_num_at_target, int pdb_num_at_ligand, int pdb_num_at_waters, int *num_at_common, double dist_tol);
re *select_sf_waters (ap *pdb_target, ap *pdb_waters, int pdb_num_at_target, int pdb_num_at_waters, int *num_sf, double dist_tol);
void print_if_waters(re *a, int num_atom, char file_name [MAX_FILENAME_LENGTH], char swt);
void prediction (int *no_is_waters, re *reference_waters, char box_name [MAX_FILENAME_LENGTH], double dist_tol, double s_tol, double p_tol, double c_tol, double b_factor, int num_reference_waters, double no_lowb_ref_wats, char reference, char selector);
void analysis (re *reference_waters, int num_reference_waters, double no_lowb_ref_wats, int *total_if);
char ***read_in_ref_ranges (char *file_name, ts *tolerance, int *ref_count_no);
char *print_pdb_file (ap *pdb, int atom_num, char file_out [MAX_FILENAME_LENGTH]);
void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );
double calc_lowb_no(re *ref_waters, int num_ref_waters, double b_fact_threshold);
ap *split_ref_file (ap *box_pdb, char **limits, int no_limits, int *atom_num_total, int atom_num_box, char label);
void str0 (char string [], int length);
int count_mdl (char *file_name);
void test_xtc();
tp *read_in_top (char *file_name, int *atom_num);
int *make_pool (char *file_in, tp *topology, int topol_anum, char t_range [MAX_FILENAME_LENGTH], char l_range [MAX_FILENAME_LENGTH], char w_range [MAX_FILENAME_LENGTH], double dist_tol, char selector, char mode);
int pdb2xtc (char *file_in, char mdl_ind);
void editing (char *file_in, char t_range [MAX_FILENAME_LENGTH], char l_range [MAX_FILENAME_LENGTH], char w_range [MAX_FILENAME_LENGTH], double min_dist_tol, double max_dist_tol);
void netdraw (char *file_in, char t_range [MAX_FILENAME_LENGTH], char l_range [MAX_FILENAME_LENGTH], char w_range [MAX_FILENAME_LENGTH], double min_dist_tol, double max_dist_tol);

