TITLE ch_HCN.mod  
 
COMMENT
Chen K, Aradi I, Thon N, Eghbal-Ahmadi M, Baram TZ, Soltesz I: Persistently
modified
h-channels after complex febrile seizures convert the seizure-induced
enhancement of
inhibition to hyperexcitability. Nature Medicine, 7(3) pp. 331-337, 2001.
(modeling by Ildiko Aradi, iaradi@uci.edu)
distal dendritic Ih channel kinetics for both HT and Control anlimals
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
	(uF) = (microfarad)
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mM) = (millimolar)
	(um) = (micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 
	SUFFIX ch_HCN 
	USEION h READ eh WRITE ih VALENCE 1
	RANGE ghbar, gh, ih
	RANGE hinf
	RANGE fast_tau, slow_tau
	RANGE myi
}
 
 
PARAMETER {
	ghbar  (mho/cm2)
	eh (mV)
}
 
STATE {
	h
}
 
ASSIGNED {
	v (mV) 
	celsius (degC)
	dt (ms)    
  
	gh (mho/cm2)
 	ih (mA/cm2)
	
	hinf 
 	fast_tau (ms)
 	slow_tau (ms) 
 	myi (mA/cm2)
} 

BREAKPOINT {
	SOLVE states METHOD cnexp
		
	gh = ghbar*h*h
	ih = gh*(v - eh)
	myi = ih
}
 
UNITSOFF
 
INITIAL { : called from hoc to calculate hinf at resting potential
	trates(v)
	h = hinf
}

DERIVATIVE states {	:computes h at current v and dt 
	trates(v)
	h' = (hinf-h)/slow_tau + (hinf-h)/fast_tau
}
 
LOCAL q10
PROCEDURE trates(v) {  :Computes rate and other constants at current v.
	TABLE hinf, fast_tau, slow_tau	
	DEPEND celsius 
	FROM -120 TO 100 WITH 220
                           
    q10 = 3^((celsius - 6.3)/10)
       
	hinf =  1 / (1 + exp( (v+91)/10 ))

	:"hyf" FAST CONTROL Hype activation system
	fast_tau = (14.9 + 14.1 / (1+exp(-(v+95.2)/0.5)))/q10

	:"hys" SLOW CONTROL Hype activation system
	slow_tau = (80 + 172.7 / (1+exp(-(v+59.3)/-0.83)))/q10

}

UNITSON
