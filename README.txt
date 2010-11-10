README for <Model-3.0.hoc>

INSTRUCTIONS FOR RUNNING THE PROGRAM
>SETUP:
1. Ensure neuron is installed on the head node and has its directory added to PATH
2. Ensure that mpich2 is installed on the head node and has its directory added to PATH
3. Compile the mod files into the directory from which the model program will be run
	a. To compile all mod files in the directory, simply enter:
		$ nrnivmodl
	b. To compile specific mod files, enter them after the command:
		$ nrnivmodl ccanl nca tca LcaMig CaBK gskch hyperde3 bgka ichan2 Gfluct2 Exp2Gid

>RUNNING:
1. Edit the parameter input files (listed below) as necessary
2. Turn on mpd:
		$ mpdboot	
3. Call the program with the number of available hosts where the '2' is:
		$ mpiexec -n 2 nrngui -nobanner -nogui -mpi Model-2.0.hoc	
4. When finished using mpi, shut down mpd:
		$ mpdallexit

VERSION HISTORY
Model-2.7
	Delay,rands,adj layer height-MC	2010-01-28
Model-2.6
	Multiple synapses allowed by MC	2009-12-03
Model-2.5
	3D model, adjustable bins by MC	2009-11-12
Model-2.4
	2D model, adjustable bins by MC	2009-11-12
Model-2.3
	more print options by MC		2009-10-15
Model-2.2
	Option-print ConMat/Pos by MC	2009-09-28
Model-2.3
	more print options by MC		2009-10-15
Model-2.2
	Option-print ConMat/Pos by MC	2009-09-28
Model-2.1
	Print outs ASAP by Marianne		2009-09-08
Model-2.0
	Restructured by Marianne Case	2009-09-02
	Parallelized by Michael Hines	2009-08-03
50knet
	Inherited from Rob Morgan		2009-04-20

PURPOSE OF CODE:
This code creates a scalable network of cells and runs a simulation.
It models the dentate gyrus using 4 cell types and can include
characteristics of an epileptic network (sprouting and cell death).

INPUT FILES:
Library files:
	nrngui.hoc
	netparmpi.hoc
	ranstream.hoc
	CellCategoryInfo.hoc

Per cell type files:
	celltemplatename.hoc
	dist_celltemplatename.hoc

Per cell type x cell type connection files:
	celltemplatename.celltemplatename
	
Parameter files:
	parameters.hoc
	cells2include.hoc

OUTPUT FILES:
spikeraster.dat (see VI.5.a)
connections.dat (see VI.5.b)
numcons.dat (see VI.5.b)
position.dat (see VI.5.c)
celltype.dat (see VI.5.d)
runtimes.dat (see VI.5.e)

CODE OUTLINE:
I.  LOAD LIBRARIES & PARAMETERS
	1. Load the main nrngui.hoc
	2. Load the template for the parallelnetmanager class (used to parallelize the code)
	3. Load the template for the randomstream class (used to generate random number streams)
	4. Load the template for CellCategoryInfo class, which generates 1 object per cell type
		to store celltype specific data and enable the removal of all hard-coded cell type refs
		
II. SET MODEL SIZE, CELL DEFINITIONS
	1. Load some parameters from file
	2. Set more static parameters directly in this code
	3. Load celltype specific data
		a. For each cell type specified in the cells2include.hoc file
			i.  Read in the cell name, # cells, layer specifier, stim input specifier
			ii. Compute the start and end of the gid # range
			iii.Create an object of the CellCategoryInfo class and store the data from i & ii
			iv. Load the hoc file containing the celltype class template
		b. Load the hoc file containing the perforant path stimulation template
	4. Kill off a % of cells (specified by sclerosis factor) that specified as being in the hilar layer
	5. Recalculate the gid ranges for each cell type, now that the number of cells for some types has decreased
	6. Calculate the total number of cells including and excluding the perforant path stimulation cell(s)
	
III.SET UP PARALLEL CAPABILITY
	1. Set up a ParallelNetManager object
	2. Create a ParallelContext
	3. Call the round robin command, which distributes all cells among all processors in a round
	4. Define an iterator that can iterate over all the cells in a given range that are owned by
		the host that called the iterator
		
IV. CREATE, UNIQUELY ID, AND POSITION CELLS
	1. For each cell type defined in cells2include.hoc:
		a. For each host in the cluster:
			i.  Check that the gid is owned by the host (and it should be, because we are using the iterator)
			ii. Create a cell of that cell type and reference it in a list contained by that cellType object
			iii.Add the cell to the 'cells' list (this is something proposed by NEURON developers,
				but we don't use explicitly use it in the code)
			iv. Create an empty connection for the cell to use as a spike detector
			v.  Associate the cell with its gid and with the spike generation location (empty connection)
			vi. Calculate and store the cell's position using an algorithm based on gid, cell type,
				# of cells of that type, # available 'bins'

V.	CONNECT THE CELLS AND CONNECT THE PERFORANT PATH TO SOME CELLS
	1. For each cell type x cell type combination, load in the connection properties
		(probability, weight, delay # synapses to choose from)
		a. If the probability of connection is not 0, for a given cell type x cell type combination:
		b. Iterate through each (potential) post-synaptic cell of that type that exists on the host executing this code
			i.  Iterate through each (potential) pre-synaptic cell of that type (regardless of where it lives)
			ii. Algorithmically obtain the positions of the pre- and post- cells
			iii.Calculate the distance between the cells and obtain a probability factor based on the distance and
				distribution of axon length for the pre-synaptic cell
			iv. Multiply the distribution number by the probability of those two types of cells connecting
				and by a connection factor that was artificially added to ensure the proper amount of connectivity
				between cells in a network of a given size
			v.  Pick a random number between 0 and 1 and check if it is less than this product.
			vi. If the above statement is true, make a connection between the cells with the given weight and delay
				and add the connection to a list 'nclist' (this is something proposed by NEURON developers)
	2. For each perforant path stimulator cell:
		a. Connect the cell with the middle 10% of granule cells in the network
		b. Connect the cell with the middle 10% of basket cells in the network
		c. Make 10 connections to the middle 10% of mossy cells in the network (note that with small network sizes,
			a given cell may receive more than 1 connection)
			
VI.	INITIALIZE AND RUN NETWORK, OUTPUT RESULT FILES
	1. Initialize the network using parameters specified in part I.
	2. Run a low resolution 'pre-simulation' to allow cells to 'settle' and all components to reach steady state
	3. Set the program to records all spikes of all cells on the host executing this code
	4. Run the simulation for the time specified in part I, at the resolution specified in part I.
	5. Output various result files:
		a. A spike raster file giving spike times # gids of spiking cells
		b. A connection file that gives pre- and post- synaptic cell gids and synapse types
		c. A position file that gives the gid and x, y, and z coordinates of each cell
		d. A cell type file that gives cell name and gid range for each cell type
		e. A runtimes file that gives the real time taken by each code section in seconds
		