// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 14-06-2015

#include "editing.h"

void editing (char *file_in, char t_range [MAX_FILENAME_LENGTH], char l_range [MAX_FILENAME_LENGTH], char w_range [MAX_FILENAME_LENGTH], double min_dist_tol, double max_dist_tol) {

	int num_at_box, num_at_target, num_at_ligand, num_at_waters, num_at_solute, num_at_box_dtd, dist_num_wat_lig, dist_num_wat_tar, num_extr_wat, j, k, count=0, num_at_box_orig;
	ap *pdb_box = NULL;
	ap *pdb_box_dtd = NULL;
	ap *pdb_solute = NULL;
	ap *pdb_target = NULL;		
	ap *pdb_ligand = NULL;
	ap *pdb_waters = NULL;
	ap *pdb_extr_waters = NULL;
	ap *pdb_box_orig = NULL;
	ap *x = NULL;

	pd *distances_wat_tar = NULL;
	pd *distances_wat_lig = NULL;
	re *extr_waters = NULL;

	pdb_box = read_in_pdb (file_in, &num_at_box, -1);

	printf ("\nEditing...\n");

	pdb_target = box_to_pdb (pdb_box, t_range, &num_at_target, num_at_box, 't');
	pdb_ligand = box_to_pdb (pdb_box, l_range, &num_at_ligand, num_at_box, 'l');
	pdb_waters = box_to_pdb (pdb_box, w_range, &num_at_waters, num_at_box, 'w');

	pdb_solute  = merge_mols (pdb_target, pdb_ligand, num_at_target, num_at_ligand, &num_at_solute);

	distances_wat_lig = pdb_ligand_dist_full_min (pdb_ligand, pdb_waters, num_at_ligand, num_at_waters, &dist_num_wat_lig, max_dist_tol, ' ');
	distances_wat_tar = pdb_ligand_dist_full_min (pdb_target, pdb_waters, num_at_target, num_at_waters, &dist_num_wat_tar, max_dist_tol, ' ');
	extr_waters = extract_common_wat_res_min (distances_wat_tar, distances_wat_lig, dist_num_wat_tar, dist_num_wat_lig, &num_extr_wat, 'x', min_dist_tol);

	print_dist_file (distances_wat_lig, "water", "ligand", dist_num_wat_lig, 'l', 'w', max_dist_tol);
	print_dist_file (distances_wat_lig, "water", "target", dist_num_wat_lig, 't', 'w', max_dist_tol);
	print_comm_res_wat_file (extr_waters, g_root_name, num_extr_wat, max_dist_tol, min_dist_tol);

	// Create extracted water pdb
	pdb_extr_waters = (ap *) malloc(sizeof(ap));

	count = 0;
	for (j = 0; j < num_at_waters; j++)
		for (k = 0; k < num_extr_wat; k++)
			if (strcmp((extr_waters+k)->target_res_type,(pdb_waters+j)->res_type) == 0 && strcmp((extr_waters+k)->target_chain,(pdb_waters+j)->chain) == 0 && (extr_waters+k)->target_res_ser == (pdb_waters+j)->res_ser && (extr_waters+k)->target_atom_ser_no + 3 > (pdb_waters+j)->atom_ser && (extr_waters+k)->target_atom_ser_no - 3 < (pdb_waters+j)->atom_ser) {
				count++;
				x = realloc( pdb_extr_waters, count * sizeof(ap) );
				if (x) pdb_extr_waters = x;
				x = NULL;	
				*(pdb_extr_waters+count-1) = *(pdb_waters+j);
				}

	pdb_box_dtd = merge_mols (pdb_solute, pdb_extr_waters, num_at_solute, count, &num_at_box_dtd);
	PRINT_PDB_FILES ("_dtd.pdb",pdb_box_dtd,num_at_box_dtd)

	pdb_box_orig = merge_mols (pdb_solute, pdb_waters, num_at_solute, num_at_waters, &num_at_box_orig);
	PRINT_PDB_FILES ("_org.pdb",pdb_box_orig,num_at_box_orig)


	// Logfile time
	NLL
	print_log_file ( g_log_file_name, "List of atom counts\n", 'a' );
	PAR2LOGFILE_IC("Ligand all                    ",num_at_ligand)
	PAR2LOGFILE_IC("Target all                    ",num_at_target)
	PAR2LOGFILE_IC("Water all (before editing)    ",num_at_waters)
	PAR2LOGFILE_IC("Water all (after editing)     ",count)

	free(pdb_box_orig);	
	free(pdb_extr_waters);
	free(extr_waters);
	free(distances_wat_tar);
	free(distances_wat_lig);
	free(pdb_box);
	free(pdb_box_dtd);
	free(pdb_solute);
	free(pdb_target);
	free(pdb_ligand);
	free(pdb_waters);

	}
