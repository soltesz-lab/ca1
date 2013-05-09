TITLE L-type calcium channel
 
COMMENT
L-Type Ca2+ channel
From: Aradi and Holmes, 1999
Updates:
20100910-MJCASE-documentation in progress
ENDCOMMENT

VERBATIM
#include <stdlib.h> /* 	Include this library so that the following
						(innocuous) warning does not appear:
						 In function '_thread_cleanup':
						 warning: incompatible implicit declaration of 
						          built-in function 'free'  */
ENDVERBATIM
 
UNITS {
	(mA) =(milliamp)
	(mV) =(millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
 
NEURON {
	SUFFIX ch_CavLaradi				: The name of the mechanism
	:USEION nca READ enca WRITE inca VALENCE 2 
	: note that CavT additionally uses ion 'ca' and reads cai, cao
	USEION ca READ eca WRITE ica VALENCE 2 
	RANGE  g
	RANGE gmax
	RANGE cinf, ctau, dinf, dtau, inca
	RANGE myi
	THREADSAFE
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}


PARAMETER {
	v (mV) 					: membrane potential
      celsius (degC) : temperature - set in hoc; default is 6.3
	gmax (mho/cm2)		: conductance flux - defined in CavT but not here
}
 
STATE {
	e		
}
 
ASSIGNED {			: assigned (where?)
	dt (ms) 				: simulation time step

	:inca (mA/cm2)	: current flux
	ica (mA/cm2)	: current flux
	g (mho/cm2)	: conductance flux
	:enca (mV)		: reversal potential
	eca (mV)		: reversal potential

	einf
	etau (ms)
	eexp      
	myi (mA/cm2)
}

BREAKPOINT {
	SOLVE states : what is the method? let's specify one
    g = gmax*e*e
	:inca = g*(v-enca)
	ica = g*(v-eca)
	myi = ica
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	e = einf
}

? states : verbatim blocks are not thread safe (perhaps related, this mechanism cannot be used with cvode)
PROCEDURE states() {	:Computes state variables m, h, and n 
        trates(v)	:      at the current v and dt.
	e = e + eexp*(einf-e)
        :VERBATIM				
        :return 0;
        :ENDVERBATIM
}
 
LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
       :q10 = 3^((celsius - 6.3)/10)
       q10 = 3^((celsius - 34)/10)
                :"e" LCa activation system
        alpha = 15.69*vtrap(81.5-v,10)
	beta = 0.29*exp(-v/10.86)
	sum = alpha+beta        
	etau = 1/sum      einf = alpha/sum
 }

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
        TABLE  einf, eexp, etau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	: not consistently executed from here if usetable_hh == 1
		: so don't expect the tau values to be tracking along with
		: the inf values in hoc

	       tinc = -dt * q10
	eexp = 1 - exp(tinc/etau)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON

