// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "make_pool.h"

int pdb2xtc (char *file_in, char mdl_ind) {

	XDRFILE *xd;
	char file_out [MAX_FILENAME_LENGTH];
	char temporary [MAX_FILENAME_LENGTH];
	int i = g_frame_range_min;
	int j,k,result;
	int mdl_num = 0;
	int num_at_box = 0;
	int step = 0;
	float time = 0;
	float prec = 10000;	// 10 x prec of GROMACS xtc !
	matrix box;
	rvec *r = NULL;
	ap *pdb_box = NULL;

	str0 (file_out, MAX_FILENAME_LENGTH);
	str0 (temporary, MAX_FILENAME_LENGTH);
	sprintf(file_out,"%s%s%s","O_",g_root_name,".xtc\0");

	xd = xdrfile_open(file_out,"w");
	if (NULL == xd)
		die("Opening xdrfile for writing");

	while ( i >= g_frame_range_min && i <= g_frame_range_max ) {
		if ( mdl_ind == 'm' )
			pdb_box = read_in_pdb (file_in, &num_at_box, i);
		if ( mdl_ind != 'm' ) {
			sprintf(temporary,"%s%s%d%s",g_root_name,"_",i,".pdb\0");
			pdb_box = read_in_pdb (temporary, &num_at_box, -2);
			str0 (temporary, MAX_FILENAME_LENGTH);
			}	
		if ( pdb_box == NULL) { PDB_CANNOT_OPEN(file_in) exit (1); }
		mdl_num++;
		for(j=0; (j<DIM); j++)
			for(k=0; (k<DIM); k++)
				box[j][k] = (j+1)*3.7 + (k+1);		// just for fun
		r = calloc(num_at_box,sizeof(*r));
		if (NULL == r)
			die("Allocating memory for r");
		for(j=0; (j<num_at_box); j++) {
			r[j][0] = 0.1 * (pdb_box+j)->x_coord;		// 0.1 is for conversion from A to nm
			r[j][1] = 0.1 * (pdb_box+j)->y_coord;
			r[j][2] = 0.1 * (pdb_box+j)->z_coord;
			}
		result = write_xtc(xd,num_at_box,step+i-g_frame_range_min,time+i-g_frame_range_min,box,r,prec);
		if (0 != result)
			die_r("Writing xtc file",result);
		if ( i == g_frame_range_min ) PRINT_PDB_FILES ("_tpy.pdb",pdb_box,num_at_box)		// generating topology pdb
		free(r);
		free(pdb_box);
		pdb_box = NULL;
		i++;
		}

	xdrfile_close(xd);

	return(mdl_num);
	}
