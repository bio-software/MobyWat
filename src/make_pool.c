// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "make_pool.h"

int *make_pool (char *file_in, tp *topology, int topol_anum, char t_range [MAX_FILENAME_LENGTH], char l_range [MAX_FILENAME_LENGTH], char w_range [MAX_FILENAME_LENGTH], double dist_tol, char selector, char mode) {

	FILE *outfile;
	XDRFILE *xd;
	char temporary1 [MAX_FILENAME_LENGTH];
	char file_out [MAX_FILENAME_LENGTH];
	char out_file_name [MAX_FILENAME_LENGTH];
	int i, k, num_at_traj_target, num_at_traj_ligand, num_at_traj_waters, dist_num_wat_tar, dist_num_wat_lig, num_common_wat,result;
	int no_frames = g_frame_range_max - g_frame_range_min;		// +1 gives real no_frames
	int natoms = 0;
	int step = 0;
	float time = 0;
	float prec=1000;
	double *rmsd_value = NULL;
	double *rmsd_value_lig = NULL;
	int *no_is_waters = NULL;
	ap *pdb_box = NULL;
	ap *pdb_target = NULL;		// reference first frame
	ap *pdb_ligand = NULL;
	ap *pdb_traj_target = NULL;	// actual
	ap *pdb_traj_ligand = NULL;
	ap *pdb_traj_waters = NULL;
	pd *distances_wat_tar = NULL;
	pd *distances_wat_lig = NULL;
	matrix box;
	rvec *r = NULL;
	re *is_waters = NULL;
	pl *pl_waters = NULL;

	str0 (temporary1, MAX_FILENAME_LENGTH);
	str0 (file_out, MAX_FILENAME_LENGTH);
	sprintf(file_out,"%s%s%s","O_",g_root_name,".plw\0");

	no_is_waters = (int *) malloc((no_frames+1) * sizeof(int)); 
	rmsd_value = (double *) malloc((no_frames+1) * sizeof(double)); 	
	if ( selector == 'i' ) rmsd_value_lig = (double *) malloc((no_frames+1) * sizeof(double)); 	

	// Read natoms in xtc and compare with topology
	result = read_xtc_natoms(file_in,&natoms); 
	if (exdrOK != result)
		die_r("The xtc file could not be opened.",result);			// old: "read_xtc_natoms"
	r = calloc(natoms,sizeof(r[0]));
	if (NULL == r)
		die("Allocating memory for r");
	if ( natoms != topol_anum ) {
		printf("%s","\nERROR! Atom counts in topology and xtc files are not equal!\n");
		printf("%s%d\n","Atom count in topology file = ",topol_anum);
		printf("%s%d\n","Atom count in xtc file      = ",natoms);
		exit(0);
		}

	// Read xtc, xtc->ap, RMSD, pool->binary
	xd = xdrfile_open(file_in,"r");
	if (NULL == xd)
		die("Opening xdrfile for reading");
	
	k = 0;
	do {
		result = read_xtc(xd,natoms,&step,&time,box,r,&prec);	// printf("%8.3f\n",prec);
		if (exdrENDOFFILE != result) {
			if (exdrOK != result) die_r("read_xtc",result);
			if ( k >= g_frame_range_min && k <= g_frame_range_max ) {
				pdb_box = (ap *) malloc(natoms * sizeof(ap));				// xtc->ap
				for ( i = 0; i < natoms; i++) {
					(pdb_box+i)->atom_ser = (topology+i)->atom_ser;
					(pdb_box+i)->res_ser = (topology+i)->res_ser;
					(pdb_box+i)->x_coord = 10 * r[i][0];				// 10 is for conversion from nm to A
					(pdb_box+i)->y_coord = 10 * r[i][1];
					(pdb_box+i)->z_coord = 10 * r[i][2];
					strcpy((pdb_box+i)->pdb_type, (topology+i)->pdb_type);
					strcpy((pdb_box+i)->pdb_alt_loc, (topology+i)->pdb_alt_loc);
					strcpy((pdb_box+i)->res_type, (topology+i)->res_type);
					strcpy((pdb_box+i)->chain, (topology+i)->chain);
					if ( g_v == 'd' && k == g_frame_range_min ) {			
						strcpy((pdb_box+i)->pdb_achar, (topology+i)->pdb_achar);
						(pdb_box+i)->occ = (topology+i)->occ;
						(pdb_box+i)->b_factor = (topology+i)->b_factor;
						}
					}
				strcpy(temporary1, t_range);
				pdb_traj_target = box_to_pdb (pdb_box, temporary1, &num_at_traj_target, natoms, 't');
				g_input_back->target_traj_atom_no = num_at_traj_target;
				str0 (temporary1, MAX_FILENAME_LENGTH);
				if ( mode == 'p' && k == g_frame_range_min ) { 
					pdb_target = (ap *) malloc(num_at_traj_target * sizeof(ap));
					memcpy(pdb_target, pdb_traj_target, num_at_traj_target * sizeof(ap));
					}
				else if ( mode == 'a' && k == g_frame_range_min ) {
					pdb_target = (ap *) malloc(g_ref_target_no * sizeof(ap));
					memcpy(pdb_target, g_ref_target, g_ref_target_no * sizeof(ap));
					}
				strcpy(temporary1, l_range);
				if ( selector == 'i' ) {
					pdb_traj_ligand = box_to_pdb (pdb_box, temporary1, &num_at_traj_ligand, natoms, 'l');
					g_input_back->ligand_traj_atom_no = num_at_traj_ligand;
					str0 (temporary1, MAX_FILENAME_LENGTH);
					if ( mode == 'p' && k == g_frame_range_min ) { 
						pdb_ligand = (ap *) malloc(num_at_traj_ligand * sizeof(ap));
						memcpy(pdb_ligand, pdb_traj_ligand, num_at_traj_ligand * sizeof(ap));
						}
					else if ( mode == 'a' && k == g_frame_range_min ) {
						pdb_ligand = (ap *) malloc(g_ref_ligand_no * sizeof(ap));
						memcpy(pdb_ligand, g_ref_ligand, g_ref_ligand_no * sizeof(ap));
						}
					}
				strcpy(temporary1, w_range);
				pdb_traj_waters = box_to_pdb (pdb_box, temporary1, &num_at_traj_waters, natoms, 'w'); 
				g_input_back->waters_traj_atom_no = num_at_traj_waters;
				str0 (temporary1, MAX_FILENAME_LENGTH);
				free(pdb_box);
				pdb_box = NULL;

				if ( g_v == 'd' && k == g_frame_range_min ) {		
					PRINT_PDB_FILES ("_frft.pdb",pdb_traj_target,num_at_traj_target)	// Kommentelni kiíratáskor!
					if ( selector == 'i') PRINT_PDB_FILES ("_frfl.pdb",pdb_traj_ligand,num_at_traj_ligand)
					PRINT_PDB_FILES ("_frfw.pdb",pdb_traj_waters,num_at_traj_waters) 
					}

				*(rmsd_value+k-g_frame_range_min) = rmsd_calc(pdb_target,pdb_traj_target,num_at_traj_target, 'a');
				if ( selector == 'i' ) *(rmsd_value_lig+k-g_frame_range_min) = rmsd_calc(pdb_ligand,pdb_traj_ligand,num_at_traj_ligand, 'h');
	
				if ( selector == 'i' )
					distances_wat_lig = pdb_ligand_dist_full (pdb_traj_ligand, pdb_traj_waters, num_at_traj_ligand, num_at_traj_waters, &dist_num_wat_lig, dist_tol);
				distances_wat_tar = pdb_ligand_dist_full (pdb_traj_target, pdb_traj_waters, num_at_traj_target, num_at_traj_waters, &dist_num_wat_tar, dist_tol);
				if ( selector == 'i' ) {	// interface
					is_waters = extract_common_wat_res (distances_wat_tar, distances_wat_lig, dist_num_wat_tar, dist_num_wat_lig, &num_common_wat, k);
					*(no_is_waters+k-g_frame_range_min) = num_common_wat;
					}
				else {				// surface
					is_waters = copy_pd_to_re (distances_wat_tar, dist_num_wat_tar, k);
					*(no_is_waters+k-g_frame_range_min) = dist_num_wat_tar;
					}
		
				free(pdb_traj_target);			  
				pdb_traj_target = NULL;
				if ( selector == 'i' ) {
					free(pdb_traj_ligand);			  
					pdb_traj_ligand = NULL;
					}
				free(pdb_traj_waters);			  
				pdb_traj_waters = NULL;
				free(distances_wat_tar);
				distances_wat_tar = NULL;
				if ( selector == 'i' ) {
					free(distances_wat_lig);
					distances_wat_lig = NULL;
					}
		
				pl_waters = copy_re_to_pl (is_waters, *(no_is_waters+k-g_frame_range_min)); 

				if ( g_v == 'd' && k == g_frame_range_min ) {
					if ( selector != 'i' ) {
						sprintf(temporary1,"%s%s%s","O_",g_root_name,"_fSFw.pdb\0");
						print_if_waters(is_waters, *(no_is_waters+k-g_frame_range_min), temporary1, 's');
						print_log_file ( g_log_file_name, temporary1, 'a' ); NLL
						str0 (temporary1, MAX_FILENAME_LENGTH);
						}
					else {
						sprintf(temporary1,"%s%s%s","O_",g_root_name,"_fIFw.pdb\0");
						print_if_waters(is_waters, *(no_is_waters+k-g_frame_range_min), temporary1, 'i');
						print_log_file ( g_log_file_name, temporary1, 'a' ); NLL
						str0 (temporary1, MAX_FILENAME_LENGTH);
						}
					}

				free (is_waters);
				is_waters = NULL;		
	
				if ( k == g_frame_range_min ) outfile = fopen(file_out,"wb");			// Print pl_waters into file 
					else outfile = fopen(file_out,"ab");
				fwrite (pl_waters, sizeof(pl), *(no_is_waters+k-g_frame_range_min), outfile);	
				fclose(outfile);		
				free (pl_waters);
				pl_waters = NULL;		
	
				if ( k == g_frame_range_min ) printf("\n");
				printf("%-6s%-8d","Number of frames written into pool = ",k+1);
				g_frame_max_written_into_plw = k;

				if ( k < g_frame_range_max ) printf("\r");
				
			}
			else if ( ( g_v == 'v' || g_v == 'd' ) && k < g_frame_range_min ) {
				printf("%s%d%s","\nWARNING! Frame not used from xtc: ",k,"\n");
				}

		}		
		else if (exdrENDOFFILE == result && k < g_frame_range_max ) {
			printf("%s","\nERROR! Number of frames in xtc file is less than required by maximal frame number or frame range \"-n\"!\n");
			printf("%s%d\n","Final frame number read from xtc    = ",k-1);
			printf("%s%d\n","Maximal frame number (frame range)  = ",g_frame_range_max);
			exit (0);
			}
		k++;
	} while (result == exdrOK);
	xdrfile_close(xd);
	free(r);
	
	print_log_file ( g_log_file_name, file_out, 'a' ); NLL

	printf("\n");

	free(pdb_target);
	if ( selector == 'i' ) free(pdb_ligand);

	// Printing time

	if ( g_v == 'v' || g_v == 'd' ) {	// Print target and ligand rmsd file
		sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_rmst.txt\0");
		if ( mode == 'a' ) print_double_vector(rmsd_value, no_frames, out_file_name, 'r', g_frame_range_min_orig);		
		else print_double_vector(rmsd_value, no_frames, out_file_name, 't', g_frame_range_min_orig);
		print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
		str0 (out_file_name, MAX_FILENAME_LENGTH);
		if ( selector == 'i' )	{
			sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_rmsl.txt\0");
			if ( mode == 'a' ) print_double_vector(rmsd_value_lig, no_frames, out_file_name, 'd', g_frame_range_min_orig);
			else print_double_vector(rmsd_value_lig, no_frames, out_file_name, 'l', g_frame_range_min_orig);
			print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
			str0 (out_file_name, MAX_FILENAME_LENGTH);
			}
						// Print No. of IF/SF waters file
		if ( selector == 'i' )	{
			sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_NoIF.txt\0");
			print_int_vector(no_is_waters, no_frames, out_file_name, 'i', g_frame_range_min_orig);
			}
		else {
			sprintf(out_file_name,"%s%s%s","O_",g_root_name,"_NoSF.txt\0");
			print_int_vector(no_is_waters, no_frames, out_file_name, 's', g_frame_range_min_orig);
			}
		print_log_file ( g_log_file_name, out_file_name, 'a' ); NLL
		str0 (out_file_name, MAX_FILENAME_LENGTH);
		}

	free(rmsd_value);
	if ( selector == 'i' ) free(rmsd_value_lig);

	return (no_is_waters);
	}
