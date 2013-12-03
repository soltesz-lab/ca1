TITLE KDRS
: Slow K-DR current for hippocampal interneurons from Lien et al (2002)
: M.Migliore Jan. 2003

NEURON {
	SUFFIX IkdrsLAW
	USEION k READ ek WRITE ik
	RANGE  gmax, ik
	GLOBAL minf, mtau, hinf
	RANGE myi
}

PARAMETER {
	gmax = 0.0002   	(mho/cm2)	
								
	celsius		(degC)
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	a0m=0.015
	vhalfm=-25	(mV)
	zetam=0.15
	gmm=0.5
	htau=1000	(ms)
	mmin=7		(ms)
	q10=3
	f=0.93
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	minf 		mtau (ms)	 	
	hinf
	myi (mA/cm2)
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gmax*m*h*(v - ek)
	myi = ik
} 

INITIAL {
	trates(v)
	m=minf  
	h=hinf  
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(v(mV)) {  
	LOCAL qt
        qt=q10^((celsius-23(degC))/10(degC))
        minf = (1/(1 + exp(-(v+41.9(mV))/23.1(mV))))^4
	mtau = betm(v)/(qt*a0m*1(/ms)*(1+alpm(v)))
	if (mtau<mmin) {mtau=mmin}
        hinf = f*(1/(1 + exp((v+52.2(mV))/15.2(mV))))+(1-f)
}

FUNCTION alpm(v(mV)) {
  alpm = exp(zetam*(v-vhalfm)*1(/mV)) 
}

FUNCTION betm(v(mV)) {
  betm = exp(zetam*gmm*(v-vhalfm)*1(/mV)) 
}
