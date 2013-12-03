TITLE KDRF
: Fast K-DR current for hippocampal interneurons from Lien et al (2002)
: M.Migliore Jan. 2003

NEURON {
	SUFFIX IkdrfLAW
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
	a0m=0.036
	vhalfm=-33	(mV)
	zetam=0.1
	gmm=0.7
	htau=1000	(ms)
	q10=3
	f=0.92
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
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
        minf = (1/(1 + exp(-(v+36.2(mV))/16.1(mV))))^4
	mtau = betm(v)/(qt*a0m*1(/ms)*(1+alpm(v)))

        hinf = f*(1/(1 + exp((v+40.6(mV))/7.8(mV))))+(1-f)
}

FUNCTION alpm(v(mV)) {
  alpm = exp(zetam*(v-vhalfm)*1(/mV)) 
}

FUNCTION betm(v(mV)) {
  betm = exp(zetam*gmm*(v-vhalfm)*1(/mV)) 
}
