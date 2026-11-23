// MobyWat   
// Calculation of hydration structures of molecular surfaces and interfaces
// 
// Programmed by C. Hetenyi  
// Web www.mobywat.com
// E-mail questions@mobywat.com
// 
// Last modified 05-04-2015

#include "global.h"

extern char g_xtc_ind;

int **make_occ_list_bin(int frames, int *total_if, int *no_rows_ol, char box_name [MAX_FILENAME_LENGTH]) {

	FILE *infile;
	int i, j, k, l, counter = 0, label = 0;
	int offset = 0;
	int **o_list = NULL;

	char file_in [MAX_FILENAME_LENGTH];

	pl *waters = NULL;
	pl *waters_realloc = NULL;
	pl *frame = NULL;
	
	void str0 (char string [], int length);
	void insertion_sort_int_ss (int **a, int num);

	str0 (file_in, MAX_FILENAME_LENGTH);

	frames++; 		// Now, it is real frame number	

	if ( g_xtc_ind != 'w' ) sprintf(file_in,"%s%s%s","O_",box_name,".plw\0");	
		else sprintf(file_in,"%s%s",box_name,".plw\0");

	infile = fopen(file_in, "rb");
	if ( infile == NULL ) {
		printf("%s%s%s","ERROR The ",file_in," file could not be opened itt.\n");
		exit (1);
		}
	
	// Fill up the initial part of waters with the first frame having non-zero IF water number	
	for (i = 0; i < frames; i++) 
		if ( *(total_if+i) != 0 && waters == NULL ) {
			waters = (pl *) malloc(*(total_if+i) * sizeof(pl));
			offset = 0;
			for (j = 0; j < i; j++) offset = offset + *(total_if+j) * sizeof(pl);
			fseek(infile, offset, SEEK_SET);			
			fread(waters, sizeof(pl), *(total_if+i), infile);
			counter = counter + *(total_if+i);
			break;
			}

	// Add additional elements not included previously in the waters			
	for (i = 0; i < frames; i++) {
		if ( *(total_if+i) != 0 ) {
			frame = (pl *) malloc(*(total_if+i) * sizeof(pl));
			offset = 0;
			for (j = 0; j < i; j++) offset = offset + *(total_if+j) * sizeof(pl);
			fseek(infile, offset, SEEK_SET);			
			fread(frame, sizeof(pl), *(total_if+i), infile);
			for (j = 0; j < *(total_if+i); j++) {
				label = 0;
				for (k = 0; k < counter; k++) 
					if ( (frame+j)->aser == (waters+k)->aser && (frame+j)->rser == (waters+k)->rser )
						label = 1;
				if (label == 0) {
					counter++;
					waters_realloc = realloc(waters, counter * sizeof(re));
					if (waters_realloc) waters = waters_realloc;			
					waters_realloc = NULL;
					(waters+counter-1)->aser = (frame+j)->aser;
					(waters+counter-1)->rser = (frame+j)->rser;
					}
				}
			free(frame);
			frame = NULL;
			}
	 	}		

	*no_rows_ol = counter;

	o_list = (int **) malloc(3 * sizeof(int *)); 
	for (i = 0; i < 3; i++) 
		*(o_list+i) = (int *) malloc(counter * sizeof(int));

	// Just copying the atom & res serial numbers from the struct * into the int **
	for (i = 0; i < counter; i++) {
		*(*(o_list+0)+i) = (waters+i)->aser;
		*(*(o_list+1)+i) = (waters+i)->rser;
		}
	free(waters);

	// Third column holds occupancies - here collected from trajectory's IF lists
	for (l = 0; l < counter; l++) {			
		*(*(o_list+2)+l) = 0;
		for (i = 0; i < frames; i++) {
			if ( *(total_if+i) != 0 ) {
				frame = (pl *) malloc(*(total_if+i) * sizeof(pl));
				offset = 0;
				for (j = 0; j < i; j++) offset = offset + *(total_if+j) * sizeof(pl);
				fseek(infile, offset, SEEK_SET);			
				fread(frame, sizeof(pl), *(total_if+i), infile);
				label = 0;
				for (j = 0; j < *(total_if+i); j++) 
					if ((frame+j)->aser == *(*(o_list+0)+l) && (frame+j)->rser == *(*(o_list+1)+l)) {
						label = 1;
						break;
						}
				if (label == 1) (*(*(o_list+2)+l))++;
				free(frame);
				frame = NULL;
				}
	  		}		
		}
	insertion_sort_int_ss(o_list, counter);
	fclose (infile);
	return(o_list);
	}

re *make_occ_list_pos_bin(int frames, int *total_if, int *no_rows_ol, int c_tolerance, char box_name [MAX_FILENAME_LENGTH]) {

	FILE *infile;
	int i, j, k;
	int count_cluster = 0;
	int label = 0;
	int offset = 0;
	int pool_member_counter = 0;

	char file_in [MAX_FILENAME_LENGTH];

	pl *frame = NULL;
	pl **waters = NULL;
	pl **x = NULL;
	pl *y = NULL;

	re *o_list = NULL;		// ->water_no is occupancy; x_coord,y..,z.. are average coordinates
	re *z = NULL;
	
	void str0 (char string [], int length);
	void insertion_sort_re (re *a, int num, char swics);
	void mean_coords3(pl *a, int members, re *mean);
	
	str0 (file_in, MAX_FILENAME_LENGTH);

	frames++; 		// Now, it is real frame number	

	if ( g_xtc_ind != 'w' ) sprintf(file_in,"%s%s%s","O_",box_name,".plw\0");	
		else sprintf(file_in,"%s%s",box_name,".plw\0");

	infile = fopen(file_in, "rb");
	if ( infile == NULL ) {
		printf("%s%s%s","ERROR The ",file_in," file could not be opened.\n");
		exit (1);
		}

	for (i = 0; i < frames; i++) {
		if ( *(total_if+i) != 0 ) {
			pool_member_counter = 0;
			frame = (pl *) malloc(*(total_if+i) * sizeof(pl));
			offset = 0;
			for (j = 0; j < i; j++) offset = offset + *(total_if+j) * sizeof(pl);
			fseek(infile, offset, SEEK_SET);			
			fread(frame, sizeof(pl), *(total_if+i), infile);
			if (count_cluster == 0)	{				// open first cluster
				waters = (pl **) malloc(sizeof(pl *));
				*(waters+0) = (pl *) malloc(sizeof(pl));	// first element of first cluster
				o_list = (re *) malloc(sizeof(re));
				**waters = *frame;
				o_list->x_coord = frame->x_coord;
				o_list->y_coord = frame->y_coord;
				o_list->z_coord = frame->z_coord;
				(o_list+0)->water_no = 1;			// here it is count no.
				count_cluster++;
				pool_member_counter = 1;
				}
			for (j = pool_member_counter; j < *(total_if+i); j++) {
				label = 0;
				for (k = 0; k < count_cluster; k++)	// add to existing cluster if it fits to but only to one cluster (break)
					if ( (((frame+j)->x_coord - (o_list+k)->x_coord) * ((frame+j)->x_coord - (o_list+k)->x_coord) + ((frame+j)->y_coord - (o_list+k)->y_coord) * ((frame+j)->y_coord - (o_list+k)->y_coord) +  ((frame+j)->z_coord - (o_list+k)->z_coord) * ((frame+j)->z_coord - (o_list+k)->z_coord)) < c_tolerance * c_tolerance ) {
						(o_list+k)->water_no++;
						y = realloc(*(waters+k), (o_list+k)->water_no * sizeof(pl));
						if (y) *(waters+k) = y;			
						y = NULL;
						*(*(waters+k)+(o_list+k)->water_no-1) = *(frame+j);
						mean_coords3(*(waters+k), (o_list+k)->water_no, o_list+k);
						label = 1;
						break;						
						}
				if ( label == 0 ) {		// open new cluster
					count_cluster++;
					x = realloc(waters, count_cluster*sizeof(pl *));
					if (x) waters = x;
					x = NULL;
					*(waters+count_cluster-1) = (pl *) malloc(sizeof(pl));	// first element
					z = realloc(o_list, count_cluster*sizeof(re));
					if (z) o_list = z;
					z = NULL;
					*(*(waters+count_cluster-1)+0) = *(frame+j);
					(o_list+count_cluster-1)->x_coord = (frame+j)->x_coord;
					(o_list+count_cluster-1)->y_coord = (frame+j)->y_coord;
					(o_list+count_cluster-1)->z_coord = (frame+j)->z_coord;
					(o_list+count_cluster-1)->water_no = 1;	
					}					
				}
			free(frame);
			}
		}

	fclose (infile);
	*no_rows_ol = count_cluster;
	insertion_sort_re (o_list, count_cluster, 'x');
	for (k = 0; k < count_cluster; k++)		// put olist row number into cluster_no 
		(o_list+k)->cluster_no = k+1;	
	for (k = 0; k < count_cluster; k++)
		free ( *(waters+k) );
	free(waters);
	return(o_list);
	}

re *make_plist(re **mean, int *clusters, int **members, int rows, double identical_tolerance, char switch_elit, int *plist_no) {

	int i, j, k, l, max_water_ser_no = 0, max_cluster_ser_no = 0, label, stopv = 0;
	double dist;

	re *pred = NULL;

	double dist_calc (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch);

	// Allocate memory for pred
	pred = (re *) malloc(rows * sizeof(re));

	for(i = 0; i < rows; i++)
		for(j = 0; j < *(clusters+i); j++) 
			(*(mean+i)+j)->label_mind_dist = 0;		// used-up indicator set to 0

	// Elitist version = 1st clusters first

	if ( switch_elit == 'e') {
		for(i = 0; i < rows; i++) {
			if ( i == 0 ) {				// Just find the most occupied
				max_water_ser_no = 0;
				for(j = i+1; j < rows; j++)
					if (  *(*(members+j)+0) > *(*(members+max_water_ser_no)+0)  && (*(mean+j)+0)->label_mind_dist == 0 )
						max_water_ser_no = j;
				} 
			else {					// Disqualify the close ones to previous round
				for(j = 0; j < rows; j++) { 
					dist = dist_calc( (*(mean+j)+0)->x_coord, (*(mean+j)+0)->y_coord, (*(mean+j)+0)->z_coord, (pred+i-1)->x_coord, (pred+i-1)->y_coord, (pred+i-1)->z_coord, "O", "O", 'q');
					if ( dist < identical_tolerance ) 
						(*(mean+j)+0)->label_mind_dist = 1;
					}

				label = 0;			// Identify and count available ones and select first one as starting max
				for(j = 0; j < rows; j++) 				
					if ( (*(mean+j)+0)->label_mind_dist == 0 ) {
						if ( label == 0) max_water_ser_no = j;
						label++;
						}
				if ( label == 0 ) break;

				for(j = 0; j < rows; j++) 	// Select the max count
					if (  *(*(members+j)+0) > *(*(members+max_water_ser_no)+0) && (*(mean+j)+0)->label_mind_dist == 0 ) 
						max_water_ser_no = j;

				}

			(*(mean+max_water_ser_no)+0)->label_mind_dist = 1;
			*(pred+i) = *(*(mean+max_water_ser_no)+0);		
			(pred+i)->water_no = *(*(members+max_water_ser_no)+0);	//	count		
			(pred+i)->avg_water_no = (double) (pred+i)->water_no;	//	just for printing
			(pred+i)->cluster_no = max_water_ser_no+1;			// olist row ser no
			(pred+i)->cluster2_no = 1;				// olist cluster no
			}

	

		stopv = i;		// Continue with the non-first clusters if not all rows were used up

		if ( stopv < rows ) {
			for(j = 0; j < rows; j++) 				// Disqualify the close ones to available plist members up to stopv-1
				for(k = 0; k < *(clusters+j); k++)
					for(l = 0; l < stopv-1; l++) { 
						dist = dist_calc( (*(mean+j)+k)->x_coord, (*(mean+j)+k)->y_coord, (*(mean+j)+k)->z_coord, (pred+l)->x_coord, (pred+l)->y_coord, (pred+l)->z_coord, "O", "O", 'q');
							if ( dist < identical_tolerance ) 
								(*(mean+j)+k)->label_mind_dist = 1;
						}
			for(i = stopv; i < rows; i++) {
			 				
				for(j = 0; j < rows; j++) 	// Disqualify the close ones to previous round
					for(k = 0; k < *(clusters+j); k++) { 
						dist = dist_calc( (*(mean+j)+k)->x_coord, (*(mean+j)+k)->y_coord, (*(mean+j)+k)->z_coord, (pred+i-1)->x_coord, (pred+i-1)->y_coord, (pred+i-1)->z_coord, "O", "O", 'q');
						if ( dist < identical_tolerance ) 
						(*(mean+j)+k)->label_mind_dist = 1;
						}

				label = 0;			// Identify and count available ones and select first one as starting max
				for(j = 0; j < rows; j++) 
					for(k = 0; k < *(clusters+j); k++)
						if ( (*(mean+j)+k)->label_mind_dist == 0 ) {
							if ( label == 0) {
								max_water_ser_no = j;
								max_cluster_ser_no = k;
								}
							label++;
							}
				if ( label == 0 ) break;

				for(j = 0; j < rows; j++) 	// Select the max count
					for(k = 0; k < *(clusters+j); k++)
						if (  *(*(members+j)+k) > *(*(members+max_water_ser_no)+max_cluster_ser_no) && (*(mean+j)+k)->label_mind_dist == 0 ) {
							max_water_ser_no = j;
							max_cluster_ser_no = k;
							}

				(*(mean+max_water_ser_no)+max_cluster_ser_no)->label_mind_dist = 1;
				*(pred+i) = *(*(mean+max_water_ser_no)+max_cluster_ser_no);	
				(pred+i)->water_no = *(*(members+max_water_ser_no)+max_cluster_ser_no);				
				(pred+i)->avg_water_no = (double) (pred+i)->water_no;	//	just for printing
				(pred+i)->cluster_no = max_water_ser_no+1;		
				(pred+i)->cluster2_no = max_cluster_ser_no+1;
				}

			}

		*plist_no = i;
//NL	printf("%10d%10d\n",rows,i); NL

		 

		}

	else if ( switch_elit == 'a') {					// All-inclusive version: no preference for 1st memebrs of cluster lists only occupancy counts
		for(i = 0; i < rows; i++) {
			if ( i == 0 ) {				// Just find the most occupied
				max_water_ser_no = 0;
				for(j = i+1; j < rows; j++)
					if (  *(*(members+j)+0) > *(*(members+max_water_ser_no)+0)  && (*(mean+j)+0)->label_mind_dist == 0 )
						max_water_ser_no = j;
				} 
			else {
				for(j = 0; j < rows; j++) 	// Disqualify the close ones to previous round
					for(k = 0; k < *(clusters+j); k++) { 
						dist = dist_calc( (*(mean+j)+k)->x_coord, (*(mean+j)+k)->y_coord, (*(mean+j)+k)->z_coord, (pred+i-1)->x_coord, (pred+i-1)->y_coord, (pred+i-1)->z_coord, "O", "O", 'q');
						if ( dist < identical_tolerance ) 
						(*(mean+j)+k)->label_mind_dist = 1;
						}

				label = 0;			// Identify and count available ones and select first one as starting max
				for(j = 0; j < rows; j++) 
					for(k = 0; k < *(clusters+j); k++)
						if ( (*(mean+j)+k)->label_mind_dist == 0 ) {
							if ( label == 0 ) {
								max_water_ser_no = j;
								max_cluster_ser_no = k;
								}
							label++;
							}
				if ( label == 0 ) break;

				for(j = 0; j < rows; j++) 	// Select the max count
					for(k = 0; k < *(clusters+j); k++)
						if (  *(*(members+j)+k) > *(*(members+max_water_ser_no)+max_cluster_ser_no) && (*(mean+j)+k)->label_mind_dist == 0 ) {
							max_water_ser_no = j;
							max_cluster_ser_no = k;
							}
				}

			(*(mean+max_water_ser_no)+max_cluster_ser_no)->label_mind_dist = 1;
			*(pred+i) = *(*(mean+max_water_ser_no)+max_cluster_ser_no);	
			(pred+i)->water_no = *(*(members+max_water_ser_no)+max_cluster_ser_no);			
			(pred+i)->avg_water_no = (double) (pred+i)->water_no;	//	just for printing
			(pred+i)->cluster_no = max_water_ser_no + 1;
			(pred+i)->cluster2_no = max_cluster_ser_no+1;		

			}
		*plist_no = i;
//NL	printf("%10d%10d\n",rows,i); NL
		}

	for(i = 0; i < *plist_no; i++) {
		(pred+i)->plist_type = switch_elit;
		(pred+i)->plist_ser_no = i+1;
		}

	return (pred);
	}


re *make_plist_pos(re *mean, int rows, double identical_tolerance, int *plist_no) {

	int i, j, max_water_ser_no = 0, label = 0;
	double dist;

	re *pred = NULL;

	double dist_calc (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch);

	// Allocate memory for pred
	pred = (re *) malloc(rows * sizeof(re));

	for(i = 0; i < rows; i++)
		(mean+i)->label_mind_dist = 0;		// used-up indicator set to 0


	// Excluding duplicates from final list (a'la all inclusive plist)
	for(i = 0; i < rows; i++) {
		if ( i == 0 ) {				// Just find the most occupied
			max_water_ser_no = 0;
			for(j = i+1; j < rows; j++)
				if ( (mean+j)->water_no > (mean+max_water_ser_no)->water_no && (mean+j)->label_mind_dist == 0 )
					max_water_ser_no = j;
			} 
		else {
			for(j = 0; j < rows; j++) 	// Disqualify the close ones to previous round
				{ 
				dist = dist_calc( (mean+j)->x_coord, (mean+j)->y_coord, (mean+j)->z_coord, (pred+i-1)->x_coord, (pred+i-1)->y_coord, (pred+i-1)->z_coord, "O", "O", 'q');
				if ( dist < identical_tolerance ) 
				(mean+j)->label_mind_dist = 1;
				}
			label = 0;			// Identify and count available ones and select the first one as starting max
			for(j = 0; j < rows; j++) 
				if ( (mean+j)->label_mind_dist == 0 ) {
					if ( label == 0 ) max_water_ser_no = j;
					label++;
					}
			if ( label == 0 ) break;	// If all mean (=olist) list members are already used up
			for(j = 0; j < rows; j++) 	// Select the max count
				if ( (mean+j)->water_no > (mean+max_water_ser_no)->water_no && (mean+j)->label_mind_dist == 0  ) {
					max_water_ser_no = j;
					}
			}

		(mean+max_water_ser_no)->label_mind_dist = 1;
		(mean+max_water_ser_no)->avg_water_no = (double) (mean+max_water_ser_no)->water_no;	//	just for printing
		*(pred+i) = *(mean+max_water_ser_no);	
		}

	*plist_no = i;

	for(i = 0; i < *plist_no; i++) {
		(pred+i)->plist_type = 'p';
		(pred+i)->plist_ser_no = i+1;
		}

	return (pred);
	}


re *clean_plist(re *mean, int rows, double identical_tolerance, int *plist_no) {

	int i, j, max_water_ser_no = 0, label = 0;
	double dist;

	re *pred = NULL;

	double dist_calc (double x1, double y1, double z1, double x2, double y2, double z2, char pdb_type1 [PDB_TYPE_LENGTH+1], char pdb_type2 [PDB_TYPE_LENGTH+1], char h_switch);

	// Allocate memory for pred
	pred = (re *) malloc(rows * sizeof(re));

	for(i = 0; i < rows; i++)
		(mean+i)->label_mind_dist = 0;		// used-up indicator set to 0


	// Excluding duplicates from final list (a'la all inclusive plist)
	for(i = 0; i < rows; i++) {
		if ( i == 0 ) {				// Just find the most occupied
			max_water_ser_no = 0;
			for(j = i+1; j < rows; j++)
				if ( (mean+j)->avg_water_no > (mean+max_water_ser_no)->avg_water_no && (mean+j)->label_mind_dist == 0 )
					max_water_ser_no = j;
			} 
		else {
			for(j = 0; j < rows; j++) 	// Disqualify the close ones to previous round
				{ 
				dist = dist_calc( (mean+j)->x_coord, (mean+j)->y_coord, (mean+j)->z_coord, (pred+i-1)->x_coord, (pred+i-1)->y_coord, (pred+i-1)->z_coord, "O", "O", 'q');
				if ( dist < identical_tolerance ) 
				(mean+j)->label_mind_dist = 1;
				}
			label = 0;			// Identify and count available ones and select the first one as starting max
			for(j = 0; j < rows; j++) 
				if ( (mean+j)->label_mind_dist == 0 ) {
					if ( label == 0 ) max_water_ser_no = j;
					label++;
					}
			if ( label == 0 ) break;	// If all mean (=olist) list members are already used up
			for(j = 0; j < rows; j++) 	// Select the max count
				if ( (mean+j)->avg_water_no > (mean+max_water_ser_no)->avg_water_no && (mean+j)->label_mind_dist == 0  ) {
					max_water_ser_no = j;
					}
			}

		(mean+max_water_ser_no)->label_mind_dist = 1;
//		(mean+max_water_ser_no)->avg_water_no = (double) (mean+max_water_ser_no)->water_no;	//	just for printing
		*(pred+i) = *(mean+max_water_ser_no);	
		}

	*plist_no = i;

	for(i = 0; i < *plist_no; i++) {
//		(pred+i)->plist_type = 'p';
		(pred+i)->plist_ser_no = i+1;
		}

	return (pred);
	}



