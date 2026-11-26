// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 03-12-2015

#include "mobywat.h"

int main (int argc, char *argv[]) {

	str0 (g_ref_name, MAX_FILENAME_LENGTH);
	str0 (g_top_file_name, MAX_FILENAME_LENGTH);
	str0 (g_pli_name, MAX_FILENAME_LENGTH);
	str0 (g_root_name, MAX_FILENAME_LENGTH);
	
	get_args(argc, argv, &g_f, &g_n, &g_l, &g_t, &g_w, &g_b_factor, &g_e, &g_c_tolerance, &g_d_value, g_f_name, g_l_range, g_m_mode, g_n_frames, &g_s_tolerance, g_t_range, g_w_range, &g_p_tolerance, &g_v, &g_d_min_value, &g_m_tolerance);

	t0 = time(NULL); printf("MobyWat started:\n%s", ctime(&t0));

	// Program mode
	if ( strcmp(DEF_PROG_MODE,g_m_mode) == 0 ) g_m = 'p';		// pred/anal/edit/netdraw modes
		else if ( strcmp(ALT1_PROG_MODE,g_m_mode) == 0 ) g_m = 'a';
			else if ( strcmp(ALT2_PROG_MODE,g_m_mode) == 0 ) g_m = 'e';
				else if ( strcmp(ALT3_PROG_MODE,g_m_mode) == 0 ) g_m = 'd';
					else {
						printf("%s%s%s%s%s%s%s%s%s","ERROR! Program mode must be either ",DEF_PROG_MODE," or ",ALT1_PROG_MODE," or ",ALT2_PROG_MODE," or ",ALT3_PROG_MODE,"!\n");
						exit (1);
						}

	// Input files, initialize log file
	if ( g_m != 'e' && g_m != 'd' ) {
		FILENAME2GTEMP(g_f_name)
		if ( find_token_in_name ( g_temp, "mdl", g_delims2) == 'y' ) g_mdl_ind = 'm';
		FILENAME2GTEMP(g_f_name)
		if ( find_token_in_name ( g_temp, "plw", g_delims3) == 'y' ) g_xtc_ind = 'w';  
		FILENAME2GTEMP(g_f_name)
		if ( find_token_in_name ( g_temp, "xtc", g_delims3) == 'y' ) g_xtc_ind = 'x'; 
		FILENAME2GTEMP(g_f_name)
		if ( find_token_in_name ( g_temp, "pdb", g_delims3) == 'y' ) g_xtc_ind = 'p'; 
		if ( g_xtc_ind != 'w' &&  g_xtc_ind != 'x' && g_xtc_ind != 'p' ) {
			printf("%s","ERROR! The trajectory file must be one of the following three types: *.xtc, *.pdb, or *.plw!\n");
			exit (1);
			}
		FILENAME2GTEMP(g_f_name)
		if ( g_xtc_ind == 'x' ) strcpy(g_root_name, strtok(g_temp, g_delims3));
		if ( g_xtc_ind == 'p' ) strcpy(g_root_name, strtok(g_temp, g_delims2));
		if ( g_xtc_ind == 'w' ) strcpy(g_root_name, strtok(g_temp, g_delims3));
		}
		else {
			FILENAME2GTEMP(g_f_name)
			if ( find_token_in_name ( g_temp, "pdb", g_delims3) == 'y' ) g_xtc_ind = 'p'; 
			if ( g_xtc_ind != 'p' ) {
				printf("%s","ERROR! The edited file must be the type of *.pdb!\n");
				exit (1);
				}
			FILENAME2GTEMP(g_f_name)
			strcpy(g_root_name, strtok(g_temp, g_delims3));
			}
	sprintf(g_log_file_name,"%s%s%s","O_",g_root_name,".log\0");
	print_log_file ( g_log_file_name, INTRO, 'w' ); 

	// Frames
	if ( g_n == 'n' ) {
		FILENAME2GTEMP(g_n_frames)
		g_frame_range_min = atoi(strtok(g_temp, g_delims));	
		g_frame_range_max = atoi(strtok(NULL, g_delims));	
		}

	// else determine n from mdl or xtc - to be inserted

	g_frame_range_min_orig = g_frame_range_min;
	g_frame_range_max_orig = g_frame_range_max;

	g_input_back = (ir*) malloc( sizeof(ir) );

	// Reading in pool info file
	if ( g_xtc_ind == 'w') {
		if ( g_pli != 'y' )
			sprintf(g_pli_name,"%s%s",g_root_name,".pli\0");
		g_outfile = fopen(g_pli_name,"rb");
		if ( g_outfile == NULL ) {
			printf("%s%s%s","ERROR The ",g_pli_name," file could not be opened.\n");
			exit (1);
			}
		fseek(g_outfile, g_offset, SEEK_SET);
		str0 (g_temp, MAX_FILENAME_LENGTH);
		fread (g_temp, sizeof(char), MAX_FILENAME_LENGTH, g_outfile);		// First plw name
		if ( strcmp(g_temp,g_f_name) != 0 )
			printf("%s%s\n","WARNING! Pool info file was created for a pool file with different name: ",g_temp);
		g_offset = g_offset + MAX_FILENAME_LENGTH * sizeof(char);
		fseek(g_outfile, g_offset, SEEK_SET);
		fread (&g_pool_info, sizeof(pf), 1, g_outfile);				// Second parameters
		g_no_is_waters = (int *) malloc((g_pool_info.nmax - g_pool_info.nmin + 1) * sizeof(int)); 
		g_offset = g_offset + sizeof(pf);
		fseek(g_outfile, g_offset, SEEK_SET);
		fread (g_no_is_waters, (g_pool_info.nmax - g_pool_info.nmin + 1) * sizeof(int), 1, g_outfile);	// Third g_no_is_waters 
		fclose(g_outfile);
		}

	// Set pool info parameters as globals
	if ( g_xtc_ind == 'w') {
		if ( g_pool_info.dmax != g_d_value ) {
			g_d_value = g_pool_info.dmax;
			printf("%s%8.3f\n","WARNING! Distance limit is re-set according to pool info. New value: ",g_d_value);
			}
		if ( g_pool_info.nmin != g_frame_range_min ) {
			g_frame_range_min = g_pool_info.nmin;
			g_frame_range_min_orig = g_frame_range_min;
			printf("%s%8d\n","WARNING! Frame range (first frame) is re-set according to pool info. New value: ",g_frame_range_min);
			}
		if ( g_pool_info.nmax != g_frame_range_max ) {
			g_frame_range_max = g_pool_info.nmax;
			g_frame_range_max_orig = g_frame_range_max;
			printf("%s%8d\n","WARNING! Frame range (last frame) is re-set according to pool info. New value: ",g_frame_range_max);
			}
		g_t = g_pool_info.t;
		g_l = g_pool_info.l;
		g_w = g_pool_info.w;
		}

	// Set g_e (evaluation type) & check frame range
	if ( g_t == 't' && g_l == 'l' && g_w == 'w' ) g_e = 'i';	// Interface evaluation mode
		else if ( g_t == 't' && g_w == 'w' ) g_e = 's';		// Surface evaluation mode
			else { printf("%s\n", "ERROR! Target or water range was not specified!"); TERMINATE }
	if (g_frame_range_max < g_frame_range_min) {
		printf("%s","ERROR! First frame has larger serial number than last frame!\n");
		exit (1);
		}
	if ( g_frame_range_min < 0 || g_frame_range_max < 0 ) {
		printf("%s","ERROR! Frame number cannot be negative!\n");
		exit (1);
		}

	// Reading in topology 1
	if ( g_xtc_ind == 'x') {
		if ( g_tpy != 'y' ) 
			sprintf(g_top_file_name,"%s%s",g_root_name,"_tpy.pdb\0");	
		g_topology = read_in_top (g_top_file_name, &g_pdb_no);
		if ( g_topology == NULL ) {
			printf("%s%s%s","ERROR! There is no topology file ",g_top_file_name,"\n");
			exit (1);
			}
		}

	// Reading in ref file & matrix tolerances	
	if ( g_r != 'r' ) 
		sprintf(g_ref_name,"%s%s",g_root_name,"_ref.pdb\0");
	g_pdb_ref = read_in_pdb (g_ref_name, &g_pdb_ref_no, -2);
	
	g_tolerances = (ts*) malloc( sizeof(ts) );
	g_ref_count_no = (int *) malloc( 3 * sizeof(int) ); 
	g_ref_ranges = read_in_ref_ranges (g_ref_name, g_tolerances, g_ref_count_no); 

	// Check if reference structure exists
	if ( g_pdb_ref != NULL && g_ref_ranges != NULL ) g_r = 'r';	//  r = reference structure exists
		else g_r = ' ';

	// Re-set program mode
	if ( g_r == ' ' && g_m == 'a') {				// Predicition mode only if no ref file
		PDB_CANNOT_OPEN(g_ref_name)				
		printf("%s%s%s","WARNING! ",PROGRAM," cannot run in Analysis mode without reference file!\n");
		printf("%s%s%s","WARNING! Program mode is changed to ",DEF_PROG_MODE,"\n");
		g_m = 'p';
		}

	// Print program & evaluation mode
	NLL
	if ( g_m == 'p' ) PAR2LOGFILE_CC("Program Mode","Prediction")
	if ( g_m == 'p' && g_r == 'r' ) PAR2LOGFILE_CC("Program Sub-mode","Validation")
	if ( g_m == 'a' ) PAR2LOGFILE_CC("Program Mode","Analysis")
	if ( g_m == 'e' ) PAR2LOGFILE_CC("Program Mode","Editing")
	if ( g_m == 'd' ) PAR2LOGFILE_CC("Program Mode","NetDraw")
	NLL
	if ( g_e == 'i' ) PAR2LOGFILE_CC("Hydration Type","Interface (IF)")
	if ( g_e == 's' ) PAR2LOGFILE_CC("Hydration Type","Surface (SF)")
	NLL
	print_log_file ( g_log_file_name, "Input File(s)\n", 'a' );
	if ( g_xtc_ind == 'x')  PAR2LOGFILE_CC("Frames from xdr binary    ",g_f_name)
	if ( g_xtc_ind == 'p' && g_mdl_ind != 'm' && g_m != 'e' && g_m != 'd' ) PAR2LOGFILE_CC("Separate frames, root name ",g_root_name)
	if ( g_xtc_ind == 'p' && g_mdl_ind == 'm' && g_m != 'e' && g_m != 'd' ) PAR2LOGFILE_CC("Multiple frames (models)   ",g_f_name)
	if ( g_m == 'e' || g_m == 'd' ) PAR2LOGFILE_CC("Single pdb file of complex ",g_f_name)
	if ( g_xtc_ind == 'w' ) PAR2LOGFILE_CC("Water pool file            ",g_f_name)
	if ( g_xtc_ind == 'w' ) PAR2LOGFILE_CC("Water pool information file      ",g_pli_name)
	if ( g_r == 'r' )       PAR2LOGFILE_CC("Reference structure        ",g_ref_name)
	if ( g_xtc_ind == 'x')  PAR2LOGFILE_CC("Topology structure         ",g_top_file_name)
	NLL
	print_log_file ( g_log_file_name, "Actual Parameters\n", 'a' );
	if ( g_r == 'r' ) PAR2LOGFILE_DC("B-factor limit, A^2             ",g_b_factor)
	if ( g_m == 'p' ) switch (g_c_scheme) { case 'a': PAR2LOGFILE_CC("Clustering algorithm","IDa") break;
						case 'e': PAR2LOGFILE_CC("Clustering algorithm","IDe") break;
						case 'p': PAR2LOGFILE_CC("Clustering algorithm","POS") break;
						case 'm': PAR2LOGFILE_CC("Clustering algorithm","IDa, IDe, POS, MER") break;
						}
	if ( g_m == 'p' ) PAR2LOGFILE_DC("Clustering tolerance, A   ",g_c_tolerance)
	PAR2LOGFILE_DC("Max. distance limit, A     ",g_d_value)
	if ( g_m == 'e' ) PAR2LOGFILE_DC("Min. distance limit, A     ",g_d_min_value)
	if ( g_r == 'r' ) PAR2LOGFILE_DC("Match tolerance, A        ",g_s_tolerance)
	if ( g_m == 'd' ) PAR2LOGFILE_DC("Mobility tolerance     ",g_m_tolerance)
	if ( g_m == 'p' ) PAR2LOGFILE_DC("Prediction tolerance, A   ",g_p_tolerance)
	if ( g_m == 'p' ) PAR2LOGFILE_DC("Top cut of short prediction list, %",g_p_top)
	NLL
	if ( g_xtc_ind == 'w' ) {
		print_log_file ( g_log_file_name, "Frame range\n", 'a' );
		PAR2LOGFILE_IC("First frame (serial number) ",g_frame_range_min)
		PAR2LOGFILE_IC("Last frame (serial number) ",g_frame_range_max)
		}
	else {
		print_log_file ( g_log_file_name, "Actual Ranges (from command line)\n", 'a' );
		if ( g_m != 'e' && g_m != 'd' ) PAR2LOGFILE_CC("Frames, x-y             ",g_n_frames)
		if ( g_e == 'i')
			PAR2LOGFILE_CC("Ligand, x-y/[xy...] ",g_l_range)
		PAR2LOGFILE_CC("Target, x-y/[xy...]   ",g_t_range)
		PAR2LOGFILE_CC("Waters, x-y/WAT/Auto  ",g_w_range)
		}
	NLL
	print_log_file ( g_log_file_name, "Output File(s)\n", 'a' );

	// Handling reference file if exists & select reference waters
	if ( g_r == 'r') {
		str0 (g_temp, MAX_FILENAME_LENGTH);
		if ( g_e == 'i') g_ref_ligand = split_ref_file (g_pdb_ref, *(g_ref_ranges+0), *(g_ref_count_no+0), &g_ref_ligand_no, g_pdb_ref_no, 'l');
		g_ref_target = split_ref_file (g_pdb_ref, *(g_ref_ranges+1), *(g_ref_count_no+1), &g_ref_target_no, g_pdb_ref_no, 't');
		g_ref_waters = split_ref_file (g_pdb_ref, *(g_ref_ranges+2), *(g_ref_count_no+2), &g_ref_waters_no, g_pdb_ref_no, 'w');
		if ( g_v == 'd' ) {		
			if ( g_e == 'i') PRINT_PDB_FILES ("_refl.pdb",g_ref_ligand,g_ref_ligand_no)
			PRINT_PDB_FILES ("_reft.pdb",g_ref_target,g_ref_target_no)
			PRINT_PDB_FILES ("_refw.pdb",g_ref_waters,g_ref_waters_no) 
			}
		if ( g_e == 'i') {
			g_final_ref_waters = select_if_waters(g_ref_target, g_ref_ligand, g_ref_waters, g_ref_target_no, g_ref_ligand_no, g_ref_waters_no, &g_final_ref_waters_no, g_d_value);
			sprintf(g_temp,"%s%s%s","O_",g_root_name,"_rIFw.pdb\0");
			print_if_waters(g_final_ref_waters, g_final_ref_waters_no, g_temp, 'i');
			print_log_file ( g_log_file_name, g_temp, 'a' ); NLL
			str0 (g_temp, MAX_FILENAME_LENGTH);
			g_lowb_final_ref_wats = calc_lowb_no(g_final_ref_waters, g_final_ref_waters_no, g_b_factor);
			}
		else {  g_final_ref_waters = select_sf_waters(g_ref_target, g_ref_waters, g_ref_target_no, g_ref_waters_no, &g_final_ref_waters_no, g_d_value); 
			sprintf(g_temp,"%s%s%s","O_",g_root_name,"_rSFw.pdb\0");
			print_if_waters(g_final_ref_waters, g_final_ref_waters_no, g_temp, 's');
			print_log_file ( g_log_file_name, g_temp, 'a' ); NLL
			str0 (g_temp, MAX_FILENAME_LENGTH);
			g_lowb_final_ref_wats = calc_lowb_no(g_final_ref_waters, g_final_ref_waters_no, g_b_factor);
			}
		}

	g_tolerances->lowb_final_ref_wats = g_lowb_final_ref_wats;		//this is temporary stuff, to be removed
	free(g_pdb_ref); 

	// Handling mdl and sep pdb files, pdb2xtc conversion, generating topology pdb
	if ( g_xtc_ind == 'p' && g_m != 'e' && g_m != 'd' ) {
		g_mdl_count = pdb2xtc (g_f_name, g_mdl_ind);
		if ( g_mdl_count == 0 ) {
			if ( g_mdl_ind == 'm' ) printf("%s%s%s","ERROR! There are no frames in multiple frames file ",g_f_name,"\n");
			else printf("%s%s%s","ERROR! There are no separate frame files with root name ",g_root_name,"\n");
			exit (1);
			}
		if ( g_mdl_count < (g_frame_range_max-g_frame_range_min+1) ) {
			printf("%s","ERROR! Number of available frames is smaller than that is required by frame range \"-n\")!\n");
			exit (1);
			}
		str0 (g_f_name, MAX_FILENAME_LENGTH);
		sprintf(g_f_name,"%s%s%s","O_",g_root_name,".xtc\0");
		print_log_file ( g_log_file_name, g_f_name, 'a' ); NLL
		g_frame_range_max = g_frame_range_max - g_frame_range_min;
		g_frame_range_min = 0;
		NL ti = time(NULL); printf("Conversion to xtc done:\n%s", ctime(&ti));
		}

	// Reading in topology 2
	if ( g_xtc_ind == 'p' && g_m != 'e' && g_m != 'd' ) {
		if ( g_tpy != 'y' ) 
			sprintf(g_top_file_name,"%s%s%s","O_",g_root_name,"_tpy.pdb\0");	
		g_topology = read_in_top (g_top_file_name, &g_pdb_no);
		if ( g_topology == NULL ) {
			printf("%s%s%s","ERROR! There is no topology file ",g_top_file_name,"\n");
			exit (1);
			}
		}

	// Make pool & write pool_info file
	if ( g_xtc_ind != 'w' && g_m != 'e' && g_m != 'd' ) {
		g_no_is_waters = make_pool (g_f_name, g_topology, g_pdb_no, g_t_range, g_l_range, g_w_range, g_d_value, g_e, g_m);
		if ( g_pli != 'y' )
			sprintf(g_pli_name,"%s%s%s","O_",g_root_name,".pli\0");
		g_outfile = fopen(g_pli_name,"wb");
		str0 (g_temp, MAX_FILENAME_LENGTH);
		sprintf(g_temp,"%s%s%s","O_",g_root_name,".plw\0");
		fwrite (g_temp, sizeof(char), MAX_FILENAME_LENGTH, g_outfile);		// First plw name
		g_pool_info = (pf) { .dmax=g_d_value, .nmin=0, .nmax=g_frame_max_written_into_plw, .t=g_t, .l=g_l, .w=g_w }; 
		fwrite (&g_pool_info, sizeof(pf), 1, g_outfile);			// Second parameters
		fwrite (g_no_is_waters, (g_frame_range_max - g_frame_range_min + 1) * sizeof(int), 1, g_outfile);	// Third g_no_is_waters 
		fclose(g_outfile);
		print_log_file ( g_log_file_name, g_pli_name, 'a' ); NLL
		NL ti = time(NULL); printf("Pool in binary file done:\n%s", ctime(&ti));
		}

	// exit if no water / no frame - to be inserted

	// Prediction, analysis, editing, netdraw
		switch (g_m) {
				case 'p': if ( g_e == 'i') prediction (g_no_is_waters, g_final_ref_waters, g_root_name, g_d_value, g_s_tolerance, g_p_tolerance, g_c_tolerance, g_b_factor, g_final_ref_waters_no, g_lowb_final_ref_wats, g_r, 'i');
						else prediction (g_no_is_waters, g_final_ref_waters, g_root_name, g_d_value, g_s_tolerance, g_p_tolerance, g_c_tolerance, g_b_factor, g_final_ref_waters_no, g_lowb_final_ref_wats, g_r, 's');
					  break;
				case 'a': analysis (g_final_ref_waters, g_final_ref_waters_no, g_lowb_final_ref_wats, g_no_is_waters);
					  break;
				case 'e': editing (g_f_name, g_t_range, g_l_range, g_w_range, g_d_min_value, g_d_value);
					  break;
				case 'd': netdraw (g_f_name, g_t_range, g_l_range, g_w_range, g_d_min_value, g_d_value);
					  break;
				default : printf("ERROR! Undefined program mode!\n"); HELP
					  exit (1);
				}

	// Close logfile	
	if ( g_xtc_ind != 'w' && g_m != 'e' && g_m != 'd' ) { NLL
		print_log_file ( g_log_file_name, "Last frame (atom counts are listed)\n", 'a' );
		if ( g_e == 'i' && g_m == 'p' )
			PAR2LOGFILE_IC("Ligand all              ",g_input_back->ligand_traj_atom_no)
		PAR2LOGFILE_IC("Target all               ",g_input_back->target_traj_atom_no)
		if ( g_m == 'p' || (g_m == 'a' && g_v != 's') ) PAR2LOGFILE_IC("Target CA's               ",g_input_back->target_traj_ca_no)
		PAR2LOGFILE_IC("Waters                   ",g_input_back->waters_traj_atom_no)
		}

	if ( g_r == 'r') {
		NLL
		print_log_file ( g_log_file_name, "Reference structure (atom counts are listed)\n", 'a' );
		if ( g_e == 'i')
			PAR2LOGFILE_IC("Ligand all                ",g_ref_ligand_no)
		PAR2LOGFILE_IC("Target all                ",g_ref_target_no)
		PAR2LOGFILE_IC("Target CA's               ",g_ref_target_ca_no)
		PAR2LOGFILE_IC("Waters all                ",g_ref_waters_no)
		if ( g_e == 'i') PAR2LOGFILE_IC("Selected interface waters ",g_final_ref_waters_no)
			else PAR2LOGFILE_IC("Selected surface waters    ",g_final_ref_waters_no)
		PAR2LOGFILE_IC("Waters below B-factor limit ",(int) g_lowb_final_ref_wats)
		
		if ( g_v != DEF_V_MODE && g_ref_target_no != g_input_back->target_traj_atom_no)
			print_log_file ( g_log_file_name, "WARNING! Target atom counts in the last frame and reference structure are not equal!\n", 'a' );
		if ( g_v != DEF_V_MODE && g_e == 'i' && g_m == 'p' && g_ref_ligand_no != g_input_back->ligand_traj_atom_no)
			print_log_file ( g_log_file_name, "WARNING! Ligand atom counts in the last frame and reference structure are not equal!\n", 'a' );
		}

	if ( g_m != 'e' && g_m != 'd' && ( g_xtc_ind == 'p'|| g_xtc_ind == 'x' ) ) { NLL
		print_log_file ( g_log_file_name, "New frame range in binary output file(s)\n", 'a' );
		PAR2LOGFILE_IC("First frame (serial number) ",0)
		if ( g_xtc_ind == 'p' ) PAR2LOGFILE_IC("Last frame (serial number) ",g_frame_range_max)
		if ( g_xtc_ind == 'x' ) PAR2LOGFILE_IC("Last frame (serial number) ",g_frame_max_written_into_plw)
		}

	free(g_ref_target);
	free(g_ref_waters);
	if ( g_e == 'i') free(g_ref_ligand);

	print_log_file ( g_log_file_name, EPILOGUE2, 'a' );
	NL t1 = time(NULL); printf("MobyWat finished:\n%s", ctime(&t1));
	NLL PAR2LOGFILE_LC("Wall clock time elapsed (s)",(long) (t1-t0))
	print_log_file ( g_log_file_name, EPILOGUE, 'a' );

	return 0;
	}




