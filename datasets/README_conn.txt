Each file in the connections folder defines properties for a connection
	between a particular type of presynaptic cell and a particular type
	of postsynaptic cell.

FILE NAME:
presynapticcell.postsynapticcell

FILE CONTENTS:
Line 1: weight of the synapse (in terms of conductance, nS)

Line 2: delay associated with synapse (in terms of time from spike in
		presynaptic cell to onset of postsynaptic response, ms)

Line 3: probability of connection (from 0 - 1, the chance of any
		particular presynaptic cell of this type being connected to any
		particular postsynaptic cell of this type, a function of the
		divergence from the presynaptic cell and convergence onto the
		postsynaptic cell)
		
Line 4:	number of connections desired. Either line 3 or line 4 is used
		to determine the connections; never both. Line 3 is used for the
		random one and line 4 is used for the randfast one.
