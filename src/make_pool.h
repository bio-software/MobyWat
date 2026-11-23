// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "global.h"
#include <limits.h>
#include <float.h>
#include "xdrfile.h"
#include "xdrfile_xtc.h"

static void _die(char *msg, int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die(msg) _die(msg,__LINE__,__FILE__)

static void _die_r(char *msg, int result,int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr,"result = %d\n",result);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die_r(msg,res) _die_r(msg,res,__LINE__,__FILE__)

extern char g_log_file_name [MAX_FILENAME_LENGTH];
extern int g_frame_range_min;
extern int g_frame_range_max;
extern ir *g_input_back;
extern char g_v;
extern char g_temp [MAX_FILENAME_LENGTH];
extern char g_root_name [MAX_FILENAME_LENGTH];
extern int g_ref_ligand_no;
extern int g_ref_target_no;
extern ap *g_ref_ligand;
extern ap *g_ref_target;
extern int g_frame_range_min_orig;
extern int g_frame_range_max_orig;
extern int g_frame_max_written_into_plw;

pd *pdb_ligand_dist_full (ap *pdb_target, ap *pdb_ligand, int pdb_num_at_target, int pdb_num_at_ligand, int *dist_num, double dist_tolerance);
re *extract_common_wat_res (pd *distances_t,  pd *distances_l, int num_dist_t, int num_dist_l, int *num_common_wat, int frame);
ap *read_in_pdb (char *file_name, int *atom_num, int mdl_serno);
ap *box_to_pdb (ap *box_pdb, char limits [MAX_FILENAME_LENGTH], int *atom_num_pdb, int atom_num_box, char label);
void print_log_file ( char file_name [MAX_FILENAME_LENGTH], char text [MAX_OUT_TEXT_LENGTH], char flag );
re *copy_pd_to_re (pd *distances_t, int num_dist_t, int frame);
void str0 (char string [], int length);
pl *copy_re_to_pl (re *residues, int num_residues);
double rmsd_calc (ap *reference, ap *actual, int atom_number, char a_type_sw);
void print_double_vector(double *a, int total, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min);
void print_int_vector(int *a, int total, char file_name [MAX_FILENAME_LENGTH], char flag, int frame_min);
char *print_pdb_file (ap *pdb, int atom_num, char file_out [MAX_FILENAME_LENGTH]);
void print_if_waters(re *a, int num_atom, char file_name [MAX_FILENAME_LENGTH], char swt);
int pdb2xtc (char *file_in, char mdl_ind);

