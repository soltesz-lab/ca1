VERBATIM
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
ENDVERBATIM

NEURON {
	SUFFIX nothing
}

VERBATIM
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
ENDVERBATIM

:* v1.fastconn(gmin, gmax, v3) takes post_gid vector v1 and pre_gid range min and max, makes conns (or at least tells which ones to make)
VERBATIM

static double get_x_pos (int gid, int ncell, int gmin, int BinNumX, int BinNumYZ, int binSizeX) { //1-gid, 2-numCells, 3-startGid, 4-X binNum, 5- binNum Y*Z, 6-binSize; Return: x position of cell
	double pos;
	int CellNum, tmp;
	CellNum=gid - gmin+1;
	tmp = floor((CellNum-1)/BinNumYZ);
	pos =  (tmp%BinNumX)*binSizeX+binSizeX/2.0;
	return pos;
}

static double get_y_pos (int gid, int ncell, int gmin, int BinNumY, int BinNumZ, int binSizeY) {// Arguments: gid, numCells, startGid, Y binNum, binNum Z, binSize; Return: y position of cell
	double pos;
	int CellNum, tmp;
	CellNum=gid - gmin+1;
	tmp = floor((CellNum-1)/BinNumZ);
	pos =  (tmp%BinNumY)*binSizeY+binSizeY/2.0;
	return pos;
}

static double get_z_pos (int gid, int ncell, int gmin, int BinNumZ, int binSizeZ, int ZHeight) {// Arguments: 1-gid, 2-numCells, 3-startGid, 4-Z binNum, 5-binSize, cell layer Zo; Return: z position of cell
	double pos;
	int CellNum;
	CellNum=gid - gmin+1;
	pos = ((CellNum-1)%BinNumZ)*binSizeZ+binSizeZ/2+ZHeight;
	return pos;
}


static int fastconn (void* vv) {
  int finalconn, ny, nz, prflag, num_pre, num_post, size, gmin, gmax, nconn, ncell,  maxd, steps, myi, postgmin;
  double *x, *y, *z, a, b, c;
  finalconn = vector_instance_px(vv, &x);
 ny = vector_arg_px(1, &y);
 nz = vector_arg_px(2, &z);
  gmin = y[0];
  gmax = y[1];
  nconn = y[2];
  ncell = y[3];
  num_post = y[4];
  maxd = y[5];
  steps = y[6];
  a = y[7];
  b = y[8];
  c = y[9];
	postgmin = y[24];

	myi=1;
  num_pre = gmax - gmin + 1;
 
  //printf("int = %d, param = %d, numpost = %d, y[2]= %2.2f  y[3]= %2.2f\n", sizeof(int), sizeof(z), sizeof(y), z[0]*1.0, z[1]*1.0);

  double prepos [num_pre][3];
  double postpos [num_post][3];
  
  double mt [steps], tu [steps], tsum, conndist;
  int k, dln [steps], dlnlen [steps];
  
  for (k=0; k<num_pre; k++) {
	prepos [k] [0] = get_x_pos(k+gmin, gmax - gmin +1, gmin, y[10], y[11]*y[12], y[13]);
	prepos [k] [1] = get_y_pos(k+gmin, gmax - gmin +1, gmin, y[11], y[12], y[14]);
	prepos [k] [2] = get_z_pos(k+gmin, gmax - gmin +1, gmin, y[12], y[15], y[16]);
  }

  for (k=0; k<num_post; k++) {
	postpos [k] [0] = get_x_pos(z[k], ncell, postgmin, y[17], y[18]*y[19], y[20]);
	postpos [k] [1] = get_y_pos(z[k], ncell, postgmin, y[18], y[19], y[21]);
	postpos [k] [2] = get_z_pos(z[k], ncell, postgmin, y[19], y[22], y[23]);
  }

  /* calculate the distribution of connections*/  
  tsum = 0.;
  for (k=0; k<steps; k++) {
	  mt[k] = maxd*1.0*(k+1)/(steps); /* mt[k] = distance step (in terms of max distance)*/
	  tu[k] = (1.0/a)*exp(-((mt[k]-b)*1.0/c)*((mt[k]-b)*1.0/c))*maxd;
		//printf("mt[k] = %2.2f\ntu[k] = (1.0/a)*exp(-((mt[k]-b)*1.0/c)*((mt[k]-b)*1.0/c))*maxd = %2.2f\n", mt[k], tu[k]);
		//printf("tu[k] = (1.0/%2.2f)*exp(-((%2.2f-%2.2f)/%2.2f)*((%2.2f-%2.2f)/%2.2f))*%d\n", a, mt[k], b, c, mt[k], b, c, maxd);
		//printf("mt[%d] = %2.2f, tu[%d] = %2.2f\n", k, mt[k], k, tu[k]);
	  tsum = tsum + tu[k];
  }
//printf("tsum = %2.2f\n", tsum);

  /* dln gives the number of connections for each distance bin, per cell*/
  for (k=0; k<steps; k++) {
	  dln[k] = round((tu[k]/tsum)*(nconn/ncell)); // , per cell
	  //dln[k] = round((tu[k]/tsum)*(nconn)); // overall 
	  //printf("dln[%d] = round(%2.2f / %2.2f * %d / %d) = %d for distance: %2.2f\n", k, tu[k], tsum, nconn, ncell, dln[k], mt[k]);
  }  

	int m, n, dl, i, q, goupto, rem, szr, szp [steps], rt, ct;
	double pl;
	k=0;
for (n=0; n<num_post; n++) { // for each post cell
//printf("for (n=0; n<num_post; n++) --> n = %d\n", n);
	double sortedpos [num_pre][steps];
	for (i=0; i< steps; i++) {
		szp [i]=0; // initialize the szp array to 0 (it holds a number per bin, telling which index of the array you are on for that bin)
	}
	double dist;
	for(m=0; m<num_pre; m++) { // for each pre cell
	// calculate the distance between the pre and post cells
		pl = sqrt((1.0*prepos[m][1] - postpos[n][1])*(prepos[m][1] - postpos[n][1])+(prepos[m][2] - postpos[n][2])*(prepos[m][2] - postpos[n][2])+(prepos[m][3] - postpos[n][3])*(prepos[m][3] - postpos[n][3]));
		for (i=0; i< steps; i++) {
			if (pl<= mt[i]) // if the distance is less than the max distance for that step
			{
				sortedpos [szp [i]] [i] = m;	// add this pre cell to this particular bin's column (the next row, which szp keeps track of)
				szp [i]++;
				break;
			}
		}
	}
	
	// now, this particular post cell has an array (sortedpos) where each
	// column contains a bunch of pre-cell gids whose distances fit within
	// that column's "step" or "distance bin"
	// There is also a dln array that gives the ideal # of connections
	// for each step
	
	rem=0;
	for (dl=0; dl<steps; dl++) {	// for each step
//printf("\tfor (dl=0; dl<steps; dl++) --> dl = %d\n", dl);

		//printf("n: %d, dl=%d, dln[dl]=%d\n", n, dl, dln[dl]);
		if (dln[dl]>0) { // if this particular step wants any connections
			/* Find all the possible connections for each distance level  */
			szr = szp [dl]; // Get the number of available connections for this step
			int r[szr]; // Define an array the length of the number of available connections
			for (i=0; i< szr; i++) { 
				r[i] =  sortedpos [i] [dl]; // Fill the array with the available connections (in terms of the pre cell)
			}

			/* this random routine allows a pre-cell to make multiple connections on the post cell*/
			int randi, tmp;
			for (i=0; i<szr-1; i++) {
				randi = rand() % szr;
				tmp = r [i];	// randomly reorganize the pre cells in the r array
				r[i] = r[randi];
				r[randi] = tmp;
			}

//printf("\t\tDesired conns: %d, Available conns: %d\n", dln[dl]+rem, szr);
			if (dln[dl]+rem>szr) {	// if the number of desired connections (ones wanted in this step, plus unmade ones from previous steps)
									// is more than the available amt, set the remainder to the excess ones that can't be made in this step
				rem=dln[dl]+rem-szr;
				goupto=szr;			// and set the number to make (based on the desired and available amts)
			} else {
				goupto=dln[dl]+rem; // set the number to make (based on the desired and available amts)
				rem=0;
			}
//printf("\t\t\tWe will make %d conns\n", goupto);

			for (q=0; q<goupto; q++) { // for each one to make, r[q] gives the precell index in the pre_pos array (this program assumes
										// the gid range is continuous from gmin to gmax arguments to this mechanism.
										// n is the post-cell here. 
				x [myi] = r[q]+gmin;
				x [myi+nconn] = z[n];
				x [myi+2*nconn] = dl+1;
				myi++;
				//printf ("between precell gid: %d and postcell gidvec idx %d, distance: %2.2f\n", r[q]+gmin, n, conndist);
				/*conn(k,1:2)=[r[w[q]] c[w[q]]];
				k++;*/
			}
		
			/* Choose nt of those connections at random  */
			/* Make each of the nt connections  */
		}
	}
}
	//size = vector_capacity(vv)*(gmax-gmin+1);
	x [0] = myi-1;
	//printf("myi = %d,  nconn = %d\n", myi, nconn);
  return finalconn; //conmat;
}

ENDVERBATIM

:* PROCEDURE install_fastconn()
PROCEDURE install_fastconn () {
  VERBATIM
  install_vector_method("fastconn", fastconn);
ENDVERBATIM
}
