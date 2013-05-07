TITLE Voltage-gated sodium channel
 
COMMENT
Voltage-gated Na+ channel
From: 
Notes: none
Updates:
	20100916 - documented and cleaned - marianne.case@uci.edu
ENDCOMMENT

COMMENT
VERBATIM
#include <stdlib.h> 
/* 	Include this library so that the following (innocuous) warning does not appear:
		In function '_thread_cleanup':
		warning: incompatible implicit declaration of built-in function 'free'  */
ENDVERBATIM
ENDCOMMENT

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
	SUFFIX ch_Navngf 
	USEION na READ ena WRITE ina VALENCE 1
	RANGE g, gmax, minf, mtau, hinf, htau, ina
	RANGE myi, offset1, offset2, offset3, offset4, slope1, slope2, slope3, slope4
:	THREADSAFE
}

PARAMETER {
:	ena  (mV)
	gmax = 7.49968  (mho/cm2)   
	offset1 = 19 (mV)
	offset2 = 19 (mV)
	offset3 = 0.5816 (mV)
	offset4 = 0.35371 (mV)
	slope1 = 0.34133 (1)
	slope2 = 0.28483 (1)
	slope3 = 0.29648 (1)
	slope4 = 3.0931  (1)
}

STATE {
	m h
}

ASSIGNED {
:	offset1 (mV)
:	offset2 (mV)
:	offset3 (mV)
:	offset4 (mV)
:	slope1 (1)
:	slope2 (1)
:	slope3 (1)
:	slope4 (1)
	ena  (mV)
	v (mV) 
	celsius (degC) : temperature - set in hoc; default is 6.3
	dt (ms) 

	g (mho/cm2)
	ina (mA/cm2)
	minf
	hinf
	mtau (ms)
	htau (ms)
	mexp
	hexp 
	myi (mA/cm2)
}

BREAKPOINT {
	SOLVE states
	g = gmax*m*m*m*h  
	ina = g*(v - ena)
	myi = ina
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {	:Computes state variables m, h, and n 
	trates(v)			:      at the current v and dt.
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
}
 
LOCAL q10	: declare outside a block so available to whole mechanism
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
:	LOCAL  alpha, beta, sum	: only available to block; must be first line in block
	LOCAL  alpha, beta, sum, tinc	: only available to block; must be first line in block

	q10 = 3^((celsius - 6.3)/10)
	:q10 = 3^((celsius - 34)/10)

	:"m" sodium activation system - act and inact cross at -40
	alpha = -1*slope1*vtrap((v+43-offset1),-5) : -0.3*vtrap((v+60-27),-5) 
	beta = slope2*vtrap((v+15-offset2),5) : 0.3*vtrap((v+60-45),5)
	sum = alpha+beta        
	mtau = 1/sum 
	minf = alpha/sum
	
	:"h" sodium inactivation system
	alpha = slope3/exp((v+65-offset3)/20)
	beta = slope4/(1+exp((v+12.5-offset4)/-10))
	sum = alpha+beta
	htau = 1/sum 
	hinf = alpha/sum 	

	tinc = -dt * q10

	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
}
 
PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc	: only available to block; must be first line in block
	TABLE minf, mexp, hinf, hexp, mtau, htau
:	DEPEND dt, celsius
	DEPEND dt, celsius, slope1, slope2, slope3, slope4, offset1, offset2, offset3, offset4
  FROM -100 TO 100 WITH 200

	rates(v)	: not consistently executed from here if usetable_hh == 1
				: so don't expect the tau values to be tracking along with
				: the inf values in hoc

:	tinc = -dt * q10

:	mexp = 1 - exp(tinc/mtau)
:	hexp = 1 - exp(tinc/htau)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{  
		vtrap = x/(exp(x/y) - 1)
	}
}
 
UNITSON

