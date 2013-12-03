TITLE KA
: K-A current for hippocampal interneurons from Lien et al (2002)
: M.Migliore Jan. 2003

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
	SUFFIX IkaLAW
	USEION k READ ek WRITE ik
	RANGE  gmax, ik
	GLOBAL minf, hinf, htau, mtau
	RANGE myi
}

PARAMETER {
	gmax = 0.0002   (mho/cm2)
	celsius		(degC)
	ek		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	a0h=0.17
	vhalfh=-105	(mV)
	q10=3
	hmin=5		(ms)
}

ASSIGNED {
	ik 	(mA/cm2)
	minf 		
	mtau 	(ms)
	hinf	 	
	htau 	(ms)
	myi (mA/cm2)
	
}
 
INITIAL {
        trates(v)
        m=minf
        h=hinf
}

STATE {m h}

BREAKPOINT {
        SOLVE state METHOD cnexp
	ik = gmax*m*h*(v - ek)
	myi = ik
} 

DERIVATIVE state {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

PROCEDURE trates(v(mV)) {  
	LOCAL qt
        qt=q10^((celsius-23(degC))/10(degC))
        minf = (1/(1 + exp(-(v+41.4(mV))/26.6(mV))))^4
	mtau=0.5(ms)/qt
        hinf = 1/(1 + exp((v+78.5(mV))/6(mV)))
	htau = a0h*1(ms/mV)*(v-vhalfh)/qt
	if (htau<hmin/qt) {htau=hmin/qt}
}

