COMMENT
This repeatconn mechanism written by marianne.case@uci.edu on April 19, 2011.
It decreases the time taken to decide which connections to make between
cell by an order of magnitude, compared to Rob & Viji's code (even the
version that I parallelized). How it works, in context with hoc:

for each cell type (as postsynaptic cell type) - hoc
> for each cell type (as presynaptic cell type) - hoc
> > make a vector of gids of all cells owned by this processor
	of the postsynaptic type - hoc
> > send vector, presynaptic cell type axonal distribution,
	number of connections desired to NMODL - hoc
> > for each postsynaptic cell of postsynaptic cell type
	owned by this processor - NMODL
> > > calculate the distance between every presynaptic cell
	  of presynaptic cell type and this postsynaptic cell - NMODL
> > > choose the desired number of connections at various distances - NMODL
	  (as defined by the presynaptic axonal distribution
	  and total # of desired connections) - NMODL
> > > randomly pick the specific connections, up to the desired
	  and/or available number for each distance - NMODL
> > > for each desired connection, add the gid
	  of the presynaptic cell to a vector - NMODL
> > return the filled vector to hoc - NMODL
> > make all the connections listed in the vector - hoc

The way to call this from hoc:
- first, install the vector method with this argument:
	install_repeatconn()
	
- then, create a vector with the parameters for the connections
	(parameters described below)
- then, create a vector with the gids of the post synaptic cells
- and finally, create an empty placeholder vector with 3x the elements
	as the number of connections desired
- finally, populate the placeholder vector with the connection
	information using the following command:
conns2make.repeatconn(params, postcellgids)

params vector has 26 elements:
0 - start gid of the presynaptic cell type
1 - end gid of the presynaptic cell type (assume continuous range)
2 - total number of connections desired from this presynaptic cell type
	to this postsynaptic cell type
3 - total number of cells of the postsynaptic type
4 - number of cells of postsynaptic type with gids owned by this host
5 - total distance over which distance distribution applies (um)
6 - resolution of the fit to the distribution, in # steps to take
7 - distribution equation coefficient a (for this presynaptic cell type)
8 - distribution equation coefficient b (for this presynaptic cell type)
9 - distribution equation coefficient c (for this presynaptic cell type)
10 - number of position bins along the X axis for presynaptic cell type
11 - number of position bins along the Y axis for presynaptic cell type
12 - number of position bins along the Z axis for presynaptic cell type
13 - length (um) of bins along the X axis for presynaptic cell type
14 - length (um) of bins along the Y axis for presynaptic cell type
15 - length (um) of bins along the Z axis for presynaptic cell type
16 - height in the Z direction, of the layer in which the
	 presynaptic cell type exists
17 - number of position bins along the X axis for postsynaptic cell type
18 - number of position bins along the Y axis for postsynaptic cell type
19 - number of position bins along the Z axis for postsynaptic cell type
20 - length (um) of bins along the X axis for presynaptic cell type
21 - length (um) of bins along the Y axis for presynaptic cell type
22 - length (um) of bins along the Z axis for presynaptic cell type
23 - height in the Z direction, of the layer in which the
	 postsynaptic cell type exists
24 - start gid of the postsynaptic cell type
25 - high index at which to start the random number stream (for each gid)
	 should be set such that there are no overlaps per each gid, so err
	 on the side of caution when setting the next high index.
ENDCOMMENT

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
//void srand(unsigned seed);
ENDVERBATIM

VERBATIM

static double get_x_pos (int gid, int gmin, int BinNumX, int BinNumYZ, int binSizeX) {
	double pos;
	int CellNum, tmp;
	CellNum=gid - gmin+1;
	tmp = floor((CellNum-1)/BinNumYZ);
	pos =  (tmp%BinNumX)*binSizeX+binSizeX/2.0;
	return pos;
}

static double get_y_pos (int gid, int gmin, int BinNumY, int BinNumZ, int binSizeY) {
	double pos;
	int CellNum, tmp;
	CellNum=gid - gmin+1;
	tmp = floor((CellNum-1)/BinNumZ);
	pos =  (tmp%BinNumY)*binSizeY+binSizeY/2.0;
	return pos;
}

static double get_z_pos (int gid, int gmin, int BinNumZ, int binSizeZ, int ZHeight) {
	double pos;
	int CellNum;
	CellNum=gid - gmin+1;
	pos = ((CellNum-1)%BinNumZ)*binSizeZ+binSizeZ/2+ZHeight;
	return pos;
}

static int repeatconn (void* vv) {
  int repeatfinal, ny, nz, num_pre, num_post, gmin, gmax, maxd, steps, myflaggy, myi, postgmin, stepover;
  double *x, *y, *z, a, b, c, nconv, ncell;

	/* Get hoc vectors into c arrays */
	repeatfinal = vector_instance_px(vv, &x); // x is an array corresponding
											// to the placeholder vector
											// of connections to make

	ny = vector_arg_px(1, &y); // y is an array of parameters
	nz = vector_arg_px(2, &z); // z is an array of the postsynaptic gids
	
	/* Load the parameters from the param array */
	gmin = y[0];	// presynaptic start gid
	gmax = y[1];	// presynaptic end gid
	num_pre = gmax - gmin + 1;	// number of presynaptic cells
	
	nconv = y[2];	// total number of desired connections
	ncell = y[3];	// total number of postsynaptic cells
	num_post = y[4];	// number of postsynaptic cells owned by this host
	maxd = y[5];	// total distance over which distribution fits
	steps = y[6];	// resolution of the distribution fit (in steps)
	a = y[7];		// distribution fit coefficient a
	b = y[8];		// distribution fit coefficient b
	c = y[9];		// distribution fit coefficient c
	postgmin = y[24];	// postsynaptic start gid
	stepover = y[26];	// buffer size for number of conns for results vector

	myi=2;	// myi will give the next index into finalconn
			// 0 is reserved for # conns to make
			// 1 is reserved for the last high index used by nrnRan4int

	/* Get positions of the presynaptic and postsynaptic cells*/
	double prepos [num_pre][3];
	double postpos [num_post][3];
	int cell;

	for (cell=0; cell<num_pre; cell++) {
		prepos [cell] [0] = get_x_pos(cell+gmin, gmin, y[10], y[11]*y[12], y[13]);
		prepos [cell] [1] = get_y_pos(cell+gmin, gmin, y[11], y[12], y[14]);
		prepos [cell] [2] = get_z_pos(cell+gmin, gmin, y[12], y[15], y[16]);
	}

	for (cell=0; cell<num_post; cell++) {
		postpos [cell] [0] = get_x_pos(z[cell], postgmin, y[17], y[18]*y[19], y[20]);
		postpos [cell] [1] = get_y_pos(z[cell], postgmin, y[18], y[19], y[21]);
		postpos [cell] [2] = get_z_pos(z[cell], postgmin, y[19], y[22], y[23]);
	}

	/* calculate the distribution of desired connections*/   
	double mt [steps], tu [steps], tsum, conndist, mytmp;
	int step, dln [steps], fdln [steps];

	tsum = 0.0;
	mytmp = 0.0;
	int maxi; maxi=0;
	for (step=0; step<steps; step++) {
		mt[step] = maxd*1.0*(step+1)/(steps); /* mt[step] = distance step (in terms of max distance)*/
		tu[step] = (1.0/a)*exp(-((mt[step]-b)*1.0/c)*((mt[step]-b)*1.0/c))*maxd;
		if (tu[step]>tu[maxi]) {
			maxi=step;
		}
		tsum = tsum + tu[step];
	}

	if (tu[maxi]/tsum*nconv < 0.5) { //tsum) nconv=nconn*1.0/ncell
		for (step=0; step<steps; step++) {
			fdln[step] = round((2.0*tu[step]/tsum)*(nconv));// the number of desired
															// connections for each
															// distance bin step, per cell
			//printf("A. tu[%d]=%f, tsum=%f, nconv=%f\n", step, tu[step], tsum, nconv);
		}
	} else {
		for (step=0; step<steps; step++) {
			fdln[step] = round((tu[step]/tsum)*(nconv));// the number of desired
															// connections for each
															// distance bin step, per cell
			//printf("B. tu[%d]=%f, tsum=%f, nconv=%f\n", step, tu[step], tsum, nconv);
		}
	}
	
	/*for (step=0; step<steps; step++) {
		printf("fdln[%d]=%d\n", step, fdln[step]);
	}*/

	/* for each postsynaptic cell, find the possible connections and
	 * make the desired number of connections where possible */   
	int m, n, i, q, goupto, rem, extra, szr, szp [steps];
	double pl;
	u_int32_t idx1, idx2, maxidx1; // high and low index (seeds) for MCell_Ran4
	maxidx1 = y[25];

	for (n=0; n<num_post; n++) { // for each post cell
		int myx = (int)z[n]; // get the gid of the current postsynaptic cell in int form
		idx1 = y[25]; 	// reset the high index for the next postsynaptic
						// cell. It should be set to a value that is 
						// certainly higher than what would have been
						// used during previous calls for this low index/gid
						// We accomplish that by setting it equal to the
						// sum of all previous numconns/ncell (for this
						// post cell type, with the previous pre cell types)?
		idx2 = myx;		// set the low index equal to the gid
		
		/*if (myx==0 || myx==100 || myx==200 || myx==300 || myx==400) {
			printf("INDEX: gid=%d, idx1=%d\n", idx2, idx1);
		}*/

		double sortedpos [num_pre][steps];
		for (step=0; step< steps; step++) {
			szp [step]=0; 	// initialize the szp array to 0
							// (it holds a number per bin, telling which
							// index of the array you are on for that bin)
							// when filling the array with available
							// connections for each bin
			dln[step] = fdln[step];		
		}
		
		double dist;
		for(m=0; m<num_pre; m++) { // for each pre cell
			// calculate the distance between the pre and post cells
			pl = sqrt((1.0*prepos[m][0] - postpos[n][0])*(prepos[m][0] - postpos[n][0])+(prepos[m][1] - postpos[n][1])*(prepos[m][1] - postpos[n][1])+(prepos[m][2] - postpos[n][2])*(prepos[m][2] - postpos[n][2]));
			for (step=0; step< steps; step++) {
				/*if (ncell==2 && num_pre==3) {
					printf("distance=%f step=%d stepmax=%f gmin=%d postg=%d\n", pl, step, mt[step], gmin, postgmin);
				}*/

				if (pl<= mt[step]) // if the distance is less than the max distance for that step
				{
					sortedpos [szp [step]] [step] = m;	// add this pre cell to this particular bin's column (the next row, which szp keeps track of)
					szp [step]++;
					break;
				}
			}
		}

		/*if (ncell==2 && num_pre==3) {
			for (step=0; step< steps; step++) {
				printf("step: %d  szp: %d\n", step, szp [step]);
			}
		}*/

		// now, this particular post cell has an array (sortedpos) where each
		// column contains a bunch of pre-cell gids whose distances fit within
		// that column's "step" or "distance bin"
		// There is also a dln array that gives the ideal # of connections
		// for each step
			
		rem=0;extra=0;
		for (step=0; step<steps; step++) {	// for each step except the last one
			szr = szp [step]; // Get the number of available connections for this step
			if (szr < 1) { //Only if there are 0 available in that step, try it in a different step
				rem=dln[step]+rem-szr;
				// check the next level for extras
				if (step<steps-1) {
					if (szp [step+1] > dln[step+1]) {
						if (szp [step+1] - dln[step+1]>rem) {
							extra = rem;
						} else {
							extra = szp [step+1] - dln[step+1];
						}
						dln[step+1] = dln[step+1] + extra;
						dln[step] = dln[step] - extra;
						rem = rem - extra;
					}
				}
				if (rem>0 && step>0) { // if that still doesn't satisfy all the remainder
					if (szp [step-1] > dln[step-1]) {
						if (szp [step-1] - dln[step-1]>rem) {
							extra = rem;
						} else {
							extra = szp [step-1] - dln[step-1];
						}
						dln[step-1] = dln[step-1] + extra;
						dln[step] = dln[step] - extra;
						rem = rem - extra;				
					}
				}
			}
		}

		/*if (ncell==2 && num_pre==3) {
			for (step=0; step< steps; step++) {
				printf("step: %d  szp: %d  dln: %d\n", step, szp [step], dln[step]);
			}
		}*/
	
		rem=0;
		for (step=0; step<steps; step++) {	// for each step
			szr = szp [step]; // Get the number of available unique pre-cells for this step
			if (dln[step]>0 && szr>0) { // if this particular step wants any connections
				/* Find all the possible connections for each distance level  */
				
				/*if (ncell==2 && num_pre==3) {
					printf("precells=%d postcells=%d step=%d szr=%d\n", num_pre, ncell, step, szr);
				}*/
				
				int r[szr]; // Define an array the length of the number of available unique pre-cells
				int rout[dln[step]]; // Define an array the length of the number of desired connections
				for (i=0; i< szr; i++) { 
					r[i] =  sortedpos [i] [step]; // Fill the array with the available unique pre-cells
				}

				/* this random routine allows a pre-cell to make multiple connections on the post cell and makes the total number of desired connections*/
				u_int32_t randi;
				for (i=0; i<dln[step]; i++) {
					randi =  nrnRan4int(&idx1, idx2) % (u_int32_t)szr; // limit to the range of indices in the r array
					rout[i] = r[randi];
				}

				for (q=0; q<dln[step]; q++) { 	// for each one to make, r[q] gives the precell index in the pre_pos array (this program assumes
											// the gid range is continuous from gmin to gmax arguments to this mechanism.
											// n is the post-cell here. 
					x [myi] = (rout[q]+gmin)*1.0;				// presynaptic gid	
					x [myi+1*stepover] = (z[n])*1.0;	// postsynaptic gid
					x [myi+2*stepover] = (step+1)*1.0;	// distance step
					myi++;
				}
			} 
			//if (num_pre>10000) {
					   // if (z[0]==21504) {
					//printf("step=%d, gid=%f, myi=%d\n", step, z[n], myi);
				//}
			//}
		}
		if (idx1>maxidx1) { maxidx1=idx1;}
	}
	x [0] = myi-2;	// fill the first element of the array (vector)
					// with the total number of connections to make,
					// which may be less than the desired number (and
					// hence the size of the array)
	x [1] = (double)maxidx1;
	return repeatfinal;
}
ENDVERBATIM

: This PROCEDURE install_repeatconn() should be called from hoc
: to make the repeatconn procedure available there

PROCEDURE install_repeatconn () {
	VERBATIM
	install_vector_method("repeatconn", repeatconn);
	ENDVERBATIM
}
