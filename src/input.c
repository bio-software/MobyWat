// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 03-12-2015

#include "input.h"

void get_args(int argc, char *argv[], char *f_value, char *n_value, char *l_value, char *t_value, char *w_value, double *b_factor, char *e_mode, double *c_tolerance, double *d_tolerance, char *f_name, char *l_range, char *m_mode, char *n_frames, double *s_tolerance, char *t_range, char *w_range, double *p_tolerance, char *v_mode, double *d_min_tolerance, double *m_tolerance) {

	int i;
	char v_value [MAX_FILENAME_LENGTH] = "";
	char c_value [MAX_FILENAME_LENGTH] = "";
	char valami = ' ';

	*e_mode = ' ';	

	if (argc < 2 && strncmp (*argv+strlen(argv[0])-NAME_LENGTH, MOBYWAT_NAME, NAME_LENGTH) == 0) {
		
		printf ("%s", INTRO);
		printf("%s", "Usage\n");
		printf("%s", "$ mobywat -x <value_x> -y <value_y> ...\n\n");
		printf("%s", "Input files\n");
		printf("%-9s%-16.16s%s", "-f",DEF_FRAME_FILENAME,"Trajectory file with frames\n");
		printf("%-9s%-16.16s%s", "-pli",DEF_PLI_FILENAME,"Pool information file\n");
		printf("%-9s%-16.16s%s", "-r",DEF_MOBYREF_FILENAME,"Reference file\n");
		printf("%-9s%-16.16s%s", "-tpy",DEF_TOPOL_FILENAME,"Topology file\n\n");
		printf("%s", "Input ranges\n");
		printf("%-9s%-16.16s%s", "-l","x-y/[xy...]","Ligand range, atom serial numbers/chain IDs\n");
		printf("%-9s%-16.16s%s", "-t","x-y/[xy...]","Target range, atom serial numbers/chain IDs\n");
		printf("%-9s%-16.16s%s", "-w ","x-y/WAT/Auto","Waters range, atom serial numbers/residue name/automatic\n\n");
		printf("%s", "Parameters & defaults\n");
		printf("%-9s%-16.3f%s", "-bmax",DEF_MAX_B_FACTOR,"B-factor limit, A^2\n");
		printf("%-9s%-16.16s%s", "-cls","IDa", "Clustering algorithm, IDa/IDe/POS/MER\n");
		printf("%-9s%-16.3f%s", "-ctol",DEF_CLUST_TOLERANCE,"Clustering tolerance, A\n");
		printf("%-9s%-16.3f%s", "-dmax",DEF_DIST_TOLERANCE,"Max. distance limit, A\n");
		printf("%-9s%-16.3f%s", "-dmin",DEF_MIN_DIST_TOLERANCE,"Min. distance limit, A\n");
		printf("%-9s%-16.16s%s", "-m",DEF_PROG_MODE, "Program mode, Analysis/Prediction/Editing/NetDraw\n");
		printf("%-9s%-16.3f%s", "-mmax",DEF_MOB_TOLERANCE,"Mobility tolerance\n");
		printf("%-9s%-16.3f%s", "-mtol",DEF_SEP_TOLERANCE,"Match tolerance, A\n");
		printf("%-9s%-16s%s", "-n",DEF_FRAME_RANGE,"Frame range, x-y\n");
		printf("%-9s%-16.3f%s", "-ptol",DEF_IDENTICAL_PREDICTION_TOLERANCE,"Prediction tolerance, A\n");
		printf("%-9s%-16.3f%s", "-ptop",DEF_PRED_LIST_LENGTH,"Top cut of short prediction list, %\n");
		printf("%-9s%-16.16s%s", "-v","Silent", "Verbosity, Silent/Verbose/Diagnostic\n");
		printf("%s",EPILOGUE2); NL
		printf("%s","For detailed instructions, please visit the web site at http://www.mobywat.com !\n");		
		printf("%s",EPILOGUE); NL
		exit (1);
		}

	for (i = 1; i < argc; i++) 
		if (argv[i][0] == '-') 
	    		switch (argv[i][1]) {

				case 'b':	READ_SWITCH_NUM("max",3,atof,*b_factor,DEF_MAX_B_FACTOR,"%s%8.3f%s","WARNING! No B-factor limit was specified at switch -bmax! The default value is used: ")
						break;

				case 'c':	READ_SWITCH_NUM("tol",3,atof,*c_tolerance,DEF_CLUST_TOLERANCE,"%s%8.3f%s","WARNING! No clustering tolerance was specified at switch -ctol! The default value is used: ")
						if ( i+1 < argc && strncmp(*(argv+i)+2,"ls",2) == 0 && argv[i+1][0] != '-' ) {
							strcpy(c_value, argv[i+1]);
							if ( strncmp(c_value, "IDa", 3) == 0 ) g_c_scheme = 'a';
								else if ( strncmp(c_value, "IDe", 3) == 0 ) g_c_scheme = 'e';
									else if ( strncmp(c_value, "POS", 3) == 0 ) g_c_scheme = 'p';
										else if ( strncmp(c_value, "MER", 3) == 0 ) g_c_scheme = 'm';
											else g_c_scheme = 'a';						
							}
						if ( strncmp(*(argv+i)+2,"ls",2) == 0 && (i+1 == argc || argv[i+1][0] == '-') ){
							g_c_scheme = 'a';
							printf("%s", "WARNING! No clustering scheme was specified at switch -cls! IDa scheme is evaluated.\n");
							}
						break;

				case 'd':	READ_SWITCH_NUM("max",3,atof,*d_tolerance,DEF_DIST_TOLERANCE,"%s%8.3f%s","WARNING! Max. distance tolerance was not specified at switch -dmax! The default value is used: ")
						READ_SWITCH_NUM2("min",3,atof,*d_min_tolerance,DEF_MIN_DIST_TOLERANCE,"%s%8.3f%s","WARNING! Min. distance tolerance was not specified at switch -dmin! The default value is used: ")
						break;			

				case 'f':	READ_SWITCH_TEXT(*f_value,'f',f_name,"WARNING! No input trajectory file name was specified at switch -f! The default file name is used: ",DEF_FRAME_FILENAME,NL)
						break; 

				case 'l':	READ_SWITCH_TEXT(*l_value,'l',l_range,"ERROR!"," No ligand range was specified at switch -l!\n",TERMINATE)
						break; 

				case 'm':	READ_SWITCH_NUM("max",3,atof,*m_tolerance,DEF_MOB_TOLERANCE,"%s%8.3f%s","WARNING! Max. mobility tolerance was specified at switch -mmax! The default value is used: ")
						READ_SWITCH_NUM("tol",3,atof,*s_tolerance,DEF_SEP_TOLERANCE,"%s%8.3f%s","WARNING! No match tolerance was specified at switch -mtol! The default value is used: ")
						if ( strlen(argv[i]) <= 3 ) READ_SWITCH_TEXT(valami,'y',m_mode,"WARNING! No program mode was specified at switch -m! The default program mode is used: ",DEF_PROG_MODE,NL)
						break; 

				case 'n':	READ_SWITCH_TEXT(*n_value,'n',n_frames,"WARNING! No frame range was specified at switch -n! The default range is used: ",DEF_FRAME_RANGE,NL)
						break;

				case 'p':	READ_SWITCH_NUM("tol",3,atof,*p_tolerance,DEF_IDENTICAL_PREDICTION_TOLERANCE,"%s%8.3f%s","WARNING! No prediction tolerance was specified at switch -ptol! The default value is used: ")
						READ_SWITCH_NUM("top",3,atof,g_p_top,DEF_P_TOP,"%s%8.3f%s","WARNING! No top cut value was specified at switch -ptop! The default value is used: ")
						READ_SWITCH_TEXT_2("li",2,g_pli,'y',g_pli_name,"WARNING! No water pool information filename was specified at switch -pli!"," The filename will be created from the pool binary filename.",NL)
						break;

				case 'r':	READ_SWITCH_TEXT(g_r,'r',g_ref_name,"WARNING! No reference filename was specified at switch -r! The default filename is used: ",DEF_MOBYREF_FILENAME,NL)
						break; 

				case 't':	READ_SWITCH_TEXT_2("py",2,g_tpy,'y',g_top_file_name,"WARNING! No topology filename was specified at switch -tpy! The default filename is used: ",DEF_TOPOL_FILENAME,NL)
						if ( strlen(argv[i]) <= 3 ) READ_SWITCH_TEXT(*t_value,'t',t_range,"ERROR!"," No target range was specified at switch -t!\n",TERMINATE)
						break; 

				case 'v':	READ_SWITCH_TEXT(valami,'y',v_value,"WARNING! No verbosity was specified at switch -v! The default verbosity is used: ","Silent",NL)
						if ( strncmp(v_value, "Silent", 6) == 0 ) *v_mode = 's';
							else if ( strncmp(v_value, "Verbose", 7) == 0 ) *v_mode = 'v';
								else if ( strncmp(v_value, "Diagnostic", 10) == 0 ) *v_mode = 'd';
						break;

				case 'w':	READ_SWITCH_TEXT(*w_value,'w',w_range,"ERROR!"," No waters range was specified at switch -w!\n",TERMINATE)
						break; 

				default:	printf("ERROR! Unknown switch: \"%s\"\n", argv[i]); HELP
						exit (1);
	    			}

	if (*f_value != 'f') printf("%s%s%s", "WARNING! No input trajectory file name was specified! The default file name is used: ", DEF_FRAME_FILENAME,"\n"); 

	}

int count_mdl (char *file_name) {

	FILE *f_pdb;
	char line [MAX_LINE_LENGTH]="";
	char model_token []="MODEL";
	char endmdl_token []="ENDMDL";
	int model_no = 0;
	int endmdl_no = 0;

	f_pdb = fopen(file_name, "r");
	if (f_pdb != NULL) { 
		while (fgets(line, sizeof(line), f_pdb)) {
			if ( strncmp (model_token, line, 5) == 0  ) model_no++; 
			if ( strncmp (endmdl_token, line, 6) == 0  ) endmdl_no++; 
			}
	if ( model_no != endmdl_no ) {
		printf("ERROR! Counts of MODEL and ENDMDL tokens are not equal. Check your input for missing tokens!\n"); HELP 
		exit (1);
		}
	else { fclose(f_pdb); return model_no; }
	
	} else return 0;

	}


ap *read_in_pdb (char *file_name, int *atom_num, int mdl_serno) {
	
	FILE *f_pdb;
	char line [MAX_LINE_LENGTH]="";
	char atom_token []="ATOM  ";
	char hetatm_token []="HETATM";
	char model_token []="MODEL";
	char endmdl_token []="ENDMDL";
	char temp [NUM_LENGTH]="";
	char temp2 [NUM_LENGTH]="";	
	char *delims = " \0";
	int counter = 0;
	int line_indicator = 0;
	int temp_mdl_no = -10000;
	int i;
	ap *molecule = NULL;
	ap *x = NULL;

	str0 (temp, NUM_LENGTH);
	str0 (temp2, NUM_LENGTH);
	str0 (line, MAX_LINE_LENGTH);

	f_pdb = fopen(file_name, "r");
	if (f_pdb != NULL) { 
		while (fgets(line, sizeof(line), f_pdb)) {
			if ( mdl_serno > -1 ) {					// mdl_serno > -1 if mdl file to be read
				if ( strncmp (model_token, line, 5) == 0 ) {
					strncpy (temp, line+10, MDL_NUM_LENGTH);
					temp [MDL_NUM_LENGTH] = '\0';
 					temp_mdl_no = atoi(temp);
					str0 (temp, NUM_LENGTH);
					if (temp_mdl_no == mdl_serno) line_indicator = 1; else line_indicator = 0;
					}
				if ( line_indicator == 1 && strncmp (endmdl_token, line, 6) == 0 ) break;
				}
			if ( mdl_serno <= -1 ) line_indicator = 1;		// mdl_serno <= -1 if separate pdb file to be read
			if ( line_indicator == 1 && ( strncmp (atom_token, line, PDB_TOKEN_LENGTH) == 0 || strncmp (hetatm_token, line, PDB_TOKEN_LENGTH) == 0 ) )  {
				if ( counter == 0) molecule = (ap *) malloc(sizeof(ap));
					else {
						x = realloc( molecule, (counter+1) * sizeof(ap) );
						if (x) molecule = x;
						x = NULL;
						}
				(molecule+counter)->model_ser = DEF_MODEL_NO;
				strncpy ((molecule+counter)->pdb_token, line, PDB_TOKEN_LENGTH);	
				(molecule+counter)->pdb_token[PDB_TOKEN_LENGTH] = '\0';				
				strncpy (temp, line+6, PDB_ATOM_NUM_LENGTH);
				temp [PDB_ATOM_NUM_LENGTH] = '\0';
 				(molecule+counter)->atom_ser = atoi(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+12, PDB_TYPE_LENGTH);
				strncpy (temp2, strtok(temp, delims), PDB_TYPE_LENGTH) ;			// strtok involved to delete empty first columns and spaces
				if ( temp2[0] == '0' || temp2[0] == '1' || temp2[0] == '2' || temp2[0] == '3' || temp2[0] == '4' || temp2[0] == '5' || temp2[0] == '6' || temp2[0] == '7' || temp2[0] == '8' || temp2[0] == '9' )
					strncpy ((molecule+counter)->pdb_type, temp2+1, PDB_TYPE_LENGTH); // delete numerical first columns
					else strncpy ((molecule+counter)->pdb_type, temp2, PDB_TYPE_LENGTH) ;	
				str0 (temp, NUM_LENGTH);
				str0 (temp2, NUM_LENGTH);
				(molecule+counter)->pdb_type[PDB_TYPE_LENGTH] = '\0';
				(molecule+counter)->pdb_alt_loc [0] = line[16];
				(molecule+counter)->pdb_alt_loc [PDB_ALT_LOC_LENGTH] = '\0';				
				strncpy ((molecule+counter)->res_type, line+17, RES_TYPE_LENGTH);	
				(molecule+counter)->res_type[RES_TYPE_LENGTH] = '\0';
				(molecule+counter)->chain [0] = line[21];
				(molecule+counter)->chain [CHAIN_LENGTH] = '\0';
				strncpy (temp, line+22, PDB_RES_NUM_LENGTH);			
				temp [PDB_RES_NUM_LENGTH] = '\0';
				(molecule+counter)->res_ser = atoi(temp);
				str0 (temp, NUM_LENGTH);
				(molecule+counter)->pdb_achar [0] = line[26];
				(molecule+counter)->pdb_achar [PDB_ACHAR_LENGTH] = '\0';				
				strncpy (temp, line+30, PDB_COORD_LENGTH);
				temp [PDB_COORD_LENGTH] = '\0';					
 				(molecule+counter)->x_coord = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+38, PDB_COORD_LENGTH);
				temp [PDB_COORD_LENGTH] = '\0';					
 				(molecule+counter)->y_coord = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+46, PDB_COORD_LENGTH);
				temp [PDB_COORD_LENGTH] = '\0';					
 				(molecule+counter)->z_coord = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+54, PDB_OCC_B_LENGTH);
				temp [PDB_OCC_B_LENGTH] = '\0';					
 				(molecule+counter)->occ = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+60, PDB_OCC_B_LENGTH);
				temp [PDB_OCC_B_LENGTH] = '\0';					
 				(molecule+counter)->b_factor = atof(temp);
				str0 (temp, NUM_LENGTH);
				counter++;
				}
			str0 (line, MAX_LINE_LENGTH); 
			}	
		*(int *)atom_num = counter;


		if ( (g_v == 'v' || g_v == 'd') && (molecule+counter-1)->atom_ser != *(int *)atom_num ) {
			printf("%s%s\n","\nWARNING! Number of ATOM/HETATM entries does not match with maximum atom serial number in file ",file_name); 
			printf("%s%d%s","WARNING! Number of ATOM/HETATM entries:   ",*(int *)atom_num,"\n");	
			printf("%s%d%s","WARNING! Maximum atom serial number:      ",(molecule+counter-1)->atom_ser,"\n");			
			printf("WARNING! This may happen e.g. if a unique atom serial number is assigned to TER tokens and/or >= 100000 atoms are listed.\n");
			}	
				
		// Check for alternate locations

		counter = 0;
		for (i = 0; i < *(int *)atom_num; i++)
			if ( (molecule+i)->pdb_alt_loc [0] != ' ' ) {
				if ( counter == 0 ) printf("%s%s%s", "ERROR! Problems occured during processing file ",file_name," !\n"); 
				printf("%s%d%s%d%s", "ERROR! Alternate location # ",counter+1," was detected at atom # ",(molecule+i)->atom_ser," !\n");
				counter++;
				}
		if (counter > 0) {
			printf("%s","Select only one location per atom and re-run the program!\n"); exit (1);
			}				
		fclose(f_pdb);
		return (molecule); 

	} else return (NULL); 
	}

tp *read_in_top (char *file_name, int *atom_num) {
	
	FILE *f_pdb;
	char line [MAX_LINE_LENGTH]="";
	char atom_token []="ATOM  ";
	char hetatm_token []="HETATM";
	char temp [NUM_LENGTH]="";
	char temp2 [NUM_LENGTH]="";	
	char *delims = " \0";
	int counter = 0;
	int i;
	tp *molecule = NULL;
	tp *x = NULL;

	str0 (temp, NUM_LENGTH);
	str0 (temp2, NUM_LENGTH);
	str0 (line, MAX_LINE_LENGTH);

	f_pdb = fopen(file_name, "r");
	if (f_pdb != NULL) { 
		while (fgets(line, sizeof(line), f_pdb)) {
			if ( ( strncmp (atom_token, line, PDB_TOKEN_LENGTH) == 0 || strncmp (hetatm_token, line, PDB_TOKEN_LENGTH) == 0 ) )  {
				if ( counter == 0) molecule = (tp *) malloc(sizeof(tp));
					else {
						x = realloc( molecule, (counter+1) * sizeof(tp) );
						if (x) molecule = x;
						x = NULL;
						}
				strncpy (temp, line+6, PDB_ATOM_NUM_LENGTH);
				temp [PDB_ATOM_NUM_LENGTH] = '\0';
 				(molecule+counter)->atom_ser = atoi(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+12, PDB_TYPE_LENGTH);
				strncpy (temp2, strtok(temp, delims), PDB_TYPE_LENGTH) ;			// strtok involved to delete empty first columns and spaces
				if ( temp2[0] == '0' || temp2[0] == '1' || temp2[0] == '2' || temp2[0] == '3' || temp2[0] == '4' || temp2[0] == '5' || temp2[0] == '6' || temp2[0] == '7' || temp2[0] == '8' || temp2[0] == '9' )
					strncpy ((molecule+counter)->pdb_type, temp2+1, PDB_TYPE_LENGTH); // delete numerical first columns
					else strncpy ((molecule+counter)->pdb_type, temp2, PDB_TYPE_LENGTH) ;	
				str0 (temp, NUM_LENGTH);
				str0 (temp2, NUM_LENGTH);
				(molecule+counter)->pdb_type[PDB_TYPE_LENGTH] = '\0';
				(molecule+counter)->pdb_alt_loc [0] = line[16];
				(molecule+counter)->pdb_alt_loc [PDB_ALT_LOC_LENGTH] = '\0';				
				strncpy ((molecule+counter)->res_type, line+17, RES_TYPE_LENGTH);	
				(molecule+counter)->res_type[RES_TYPE_LENGTH] = '\0';
				(molecule+counter)->chain [0] = line[21];
				(molecule+counter)->chain [CHAIN_LENGTH] = '\0';
				strncpy (temp, line+22, PDB_RES_NUM_LENGTH);			
				temp [PDB_RES_NUM_LENGTH] = '\0';
				(molecule+counter)->res_ser = atoi(temp);
				str0 (temp, NUM_LENGTH);
				(molecule+counter)->pdb_achar [0] = line[26];
				(molecule+counter)->pdb_achar [PDB_ACHAR_LENGTH] = '\0';				
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+54, PDB_OCC_B_LENGTH);
				temp [PDB_OCC_B_LENGTH] = '\0';					
 				(molecule+counter)->occ = atof(temp);
				str0 (temp, NUM_LENGTH);
				strncpy (temp, line+60, PDB_OCC_B_LENGTH);
				temp [PDB_OCC_B_LENGTH] = '\0';					
 				(molecule+counter)->b_factor = atof(temp);
				str0 (temp, NUM_LENGTH);
				counter++;
				}
			str0 (line, MAX_LINE_LENGTH); 
			}	
		*(int *)atom_num = counter;


		if ( (g_v == 'v' || g_v == 'd') && (molecule+counter-1)->atom_ser != *(int *)atom_num ) {
			printf("%s%s\n","\nWARNING! Number of ATOM/HETATM entries does not match with maximum atom serial number in file ",file_name); 
			printf("%s%d%s","WARNING! Number of ATOM/HETATM entries:   ",*(int *)atom_num,"\n");	
			printf("%s%d%s","WARNING! Maximum atom serial number:      ",(molecule+counter-1)->atom_ser,"\n");			
			printf("WARNING! This may happen e.g. if a unique atom serial number is assigned to TER tokens and/or >= 100000 atoms are listed.\n");
			}	
				
		// Check for alternate locations

		counter = 0;
		for (i = 0; i < *(int *)atom_num; i++)
			if ( (molecule+i)->pdb_alt_loc [0] != ' ' ) {
				if ( counter == 0 ) printf("%s%s%s", "ERROR! Problems occured during processing file ",file_name," !\n"); 
				printf("%s%d%s%d%s", "ERROR! Alternate location # ",counter+1," was detected at atom # ",(molecule+i)->atom_ser," !\n");
				counter++;
				}
		if (counter > 0) {
			printf("%s","Select only one location per atom and re-run the program!\n"); exit (1);
			}				
		fclose(f_pdb);
		return (molecule); 

	} else return (NULL); 
	}

char ***read_in_ref_ranges (char *file_name, ts *tolerance, int *ref_count_no) {
	
	FILE *f_pdb;
	int lig_count = 0;	// 0
	int tar_count = 0;	// 1
	int wat_count = 0;	// 2

	char line [MAX_LINE_LENGTH]="";
	char temporary1 [MAX_FILENAME_LENGTH];
	char remark_token [] = "REMARK";
	char ligand_token [] = "mobywat_reference_ligand";
	char target_token [] = "mobywat_reference_target";
	char waters_token [] = "mobywat_reference_waters";
	char min_ctol_token [] = "mobywat_min_ctol";
	char min_ptol_token [] = "mobywat_min_ptol";
	char min_stol_token [] = "mobywat_min_stol";
	char num_ctol_token [] = "mobywat_num_ctol";
	char num_ptol_token [] = "mobywat_num_ptol";
	char num_stol_token [] = "mobywat_num_stol";
	char step_ctol_token [] = "mobywat_step_ctol";
	char step_ptol_token [] = "mobywat_step_ptol";
	char step_stol_token [] = "mobywat_step_stol";
	char *delims = " \n\0";

	char ***ranges = NULL;
	char **x = NULL;

	ranges = (char ***) malloc(3 * sizeof(char **));

	tolerance->min_ctol = DEF_MIN_CTOL;
	tolerance->min_ptol = DEF_MIN_PTOL;
	tolerance->min_stol = DEF_MIN_STOL;
	tolerance->num_ctol = DEF_NUM_CTOL;
	tolerance->num_ptol = DEF_NUM_PTOL;
	tolerance->num_stol = DEF_NUM_STOL;
	tolerance->step_ctol = DEF_STEP_CTOL;
	tolerance->step_ptol = DEF_STEP_PTOL;
	tolerance->step_stol = DEF_STEP_STOL; 

	str0 (line, MAX_LINE_LENGTH);
	str0 (temporary1, MAX_FILENAME_LENGTH);

	f_pdb = fopen(file_name, "r");
	if (f_pdb != NULL) {
		while (fgets(line, sizeof(line), f_pdb)) {
			if ( strncmp (remark_token, line, REMARK_TOKEN_LENGTH) == 0 )  {
				strtok(line, delims);
				str0 (temporary1, MAX_FILENAME_LENGTH);
				strcpy(temporary1, strtok(NULL, delims));
				if ( strcmp (ligand_token, temporary1) == 0 ) {				
					if ( lig_count == 0) {
						*(ranges+0) = (char **) malloc( sizeof(char *) ); 
						*(*(ranges+0)+lig_count) = (char *) malloc( MAX_FILENAME_LENGTH * sizeof (char) );
						}
					else {
						x = realloc( *(ranges+0), (lig_count+1) * sizeof(char *) );
						if (x) *(ranges+0) = x;
						x = NULL;
						*(*(ranges+0)+lig_count) = (char *) malloc( MAX_FILENAME_LENGTH * sizeof (char) );
						}
					strcpy( *(*(ranges+0)+lig_count), strtok(NULL, delims));
					lig_count++; 
					} 
				if ( strcmp (target_token, temporary1) == 0 ) {				
					if ( tar_count == 0) {
						*(ranges+1) = (char **) malloc( sizeof(char *) ); 
						*(*(ranges+1)+tar_count) = (char *) malloc( MAX_FILENAME_LENGTH * sizeof (char) );
						}
					else {
						x = realloc( *(ranges+1), (tar_count+1) * sizeof(char *) );
						if (x) *(ranges+1) = x;
						x = NULL;
						*(*(ranges+1)+tar_count) = (char *) malloc( MAX_FILENAME_LENGTH * sizeof (char) );
						}
					strcpy( *(*(ranges+1)+tar_count), strtok(NULL, delims));
					tar_count++;
					}
				if ( strcmp (waters_token, temporary1) == 0 ) {				
					if ( wat_count == 0) {
						*(ranges+2) = (char **) malloc( sizeof(char *) ); 
						*(*(ranges+2)+wat_count) = (char *) malloc( MAX_FILENAME_LENGTH * sizeof (char) );
						}
					else {
						x = realloc( *(ranges+2), (wat_count+1) * sizeof(char *) );
						if (x) *(ranges+2) = x;
						x = NULL;
						*(*(ranges+2)+wat_count) = (char *) malloc( MAX_FILENAME_LENGTH * sizeof (char) );
						}
					strcpy( *(*(ranges+2)+wat_count), strtok(NULL, delims));
					wat_count++;
					}  
				if ( strcmp (min_ctol_token, temporary1) == 0 )
					tolerance->min_ctol = atof( strtok(NULL, delims) ); 
				if ( strcmp (min_ptol_token, temporary1) == 0 )
					tolerance->min_ptol = atof( strtok(NULL, delims) ); 
				if ( strcmp (min_stol_token, temporary1) == 0 )
					tolerance->min_stol = atof( strtok(NULL, delims) ); 
				if ( strcmp (step_ctol_token, temporary1) == 0 )
					tolerance->step_ctol = atof( strtok(NULL, delims) ); 
				if ( strcmp (step_ptol_token, temporary1) == 0 )
					tolerance->step_ptol = atof( strtok(NULL, delims) ); 
				if ( strcmp (step_stol_token, temporary1) == 0 )
					tolerance->step_stol = atof( strtok(NULL, delims) ); 
				if ( strcmp (num_ctol_token, temporary1) == 0 )
					tolerance->num_ctol = atoi( strtok(NULL, delims) ); 
				if ( strcmp (num_ptol_token, temporary1) == 0 )
					tolerance->num_ptol = atoi( strtok(NULL, delims) ); 
				if ( strcmp (num_stol_token, temporary1) == 0 )
					tolerance->num_stol = atoi( strtok(NULL, delims) );
				} 
			str0 (line, MAX_LINE_LENGTH);
			}	
		*(ref_count_no+0) = lig_count;
		*(ref_count_no+1) = tar_count;
		*(ref_count_no+2) = wat_count;
		fclose(f_pdb);
		return (ranges);

	} else {
//		printf("%s%s%s","WARNING The ",file_name," file could not be opened.\n");
//		fclose(f_pdb);
		return (NULL);
		}
	 
	}

char find_token_in_name (char name[MAX_FILENAME_LENGTH], char token [MAX_FILENAME_LENGTH], char *delims) {

	char *temp;
	int token_length = strlen(token);
	char token_indicator = 'n';	// not found

	temp = strtok(name, delims);

	while ( temp != NULL ) {
		if ( strncmp(temp, token, token_length) == 0 ) token_indicator = 'y'; 	// found
		temp = strtok(NULL, delims);
		}

	return token_indicator;
	}

void str0 (char string [], int length) {
	int i;
	for (i = 0; i < length; i++)
		string[i] = '\0';
	} 
