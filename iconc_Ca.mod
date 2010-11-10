TITLE intracellular calcium accumulation

COMMENT
intracellular Ca2+ accumulation
From: 
Notes:
	calcium accumulation into a volume of area*depth next to the
	membrane with a decay (time constant tau) to resting level
	given by the global calcium variable cai0_ca_ion
Updates:
	20100910-MJCASE-documented
ENDCOMMENT

VERBATIM
#include <stdlib.h> /* 	Include this library so that the following
						(innocuous) warning does not appear:
						 In function '_thread_cleanup':
						 warning: incompatible implicit declaration of 
						          built-in function 'free'  */
ENDVERBATIM

NEURON {
	SUFFIX iconc_Ca
USEION nca READ ncai, inca, enca WRITE enca, ncai VALENCE 2
USEION lca READ lcai, ilca, elca WRITE elca, lcai VALENCE 2
USEION tca READ tcai, itca, etca WRITE etca, tcai VALENCE 2
RANGE caiinf, catau, cai, ncai, lcai,tcai, eca, elca, enca, etca
    THREADSAFE
}

UNITS {
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (milli/liter)
	(mA) = (milliamp)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}

INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

PARAMETER {
      celsius (degC) : temperature - set in hoc; default is 6.3
	depth = 200 (nm)	: assume volume = area*depth
	catau = 9 (ms)
	caiinf = 50.e-6 (mM)	: takes precedence over cai0_ca_ion
			: Do not forget to initialize in hoc if different
			: from this default.
	cao = 2 (mM)
	ica (mA/cm2)
	inca (mA/cm2)
	ilca (mA/cm2)
	itca (mA/cm2)
	cai= 50.e-6 (mM)
}

ASSIGNED {
	enca (mV)
	elca (mV)
	etca (mV)
	eca (mV)
}

STATE {
	ncai (mM)
	lcai (mM)
	tcai (mM)
}

: verbatim blocks are not thread safe (perhaps related, this mechanism cannot be used with cvode)
INITIAL {
	VERBATIM	/* what is the point of this? */	
	ncai = _ion_ncai;
	lcai = _ion_lcai;
	tcai = _ion_tcai; 
	ENDVERBATIM
	ncai=caiinf/3
	lcai=caiinf/3
	tcai=caiinf/3
	cai = caiinf	
	eca = ktf() * log(cao/caiinf)	
	enca = eca
	elca = eca
	etca = eca
}


BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
	cai = ncai+lcai+tcai	
	eca = ktf() * log(cao/cai)	
	enca = eca
	elca = eca
	etca = eca
}

DERIVATIVE integrate {
ncai' = -(inca)/depth/FARADAY * (1e7) + (caiinf/3 - ncai)/catau
lcai' = -(ilca)/depth/FARADAY * (1e7) + (caiinf/3 - lcai)/catau
tcai' = -(itca)/depth/FARADAY * (1e7) + (caiinf/3 - tcai)/catau
}

FUNCTION ktf() (mV) {
	ktf = (1000)*R*(celsius +273.15)/(2*FARADAY)
} 
