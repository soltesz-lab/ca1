:Migliore file Modify by Maciej Lazarewicz (mailto:mlazarew@gmu.edu) May/16/2001

TITLE T-calcium channel
: T-type calcium channel

NEURON {
	SUFFIX catLAW
	USEION ca READ cai,cao WRITE ica
        RANGE gmax, ica
	GLOBAL m_tau, h_tau, minf, hinf
	RANGE myi
}

UNITS {
	(mA) 	= 	(milliamp)
	(mV) 	= 	(millivolt)
	FARADAY =  	(faraday)  (kilocoulombs)
	R 	= 	(k-mole) (joule/degC)
	KTOMV 	= .0853 (mV/degC)	
}

PARAMETER {
	v (mV)
	celsius 	(degC)
	gmax	=.003 	(mho/cm2)
	cai 		(mM)
	cao 		(mM)
}

STATE {	m h }

ASSIGNED {
	ica (mA/cm2)
        gcat (mho/cm2)
	minf
	hinf
	m_tau
	h_tau
	myi (mA/cm2)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gcat = gmax*m*m*h
	ica  = gcat*ghk(v,cai,cao)
	myi = ica
}

INITIAL {
	rate(v)
	m = minf
	h = hinf
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (DegC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

DERIVATIVE state {
	rate(v)
	m' = (minf-m)/m_tau
	h' = (hinf-h)/h_tau
}

PROCEDURE rate (v (mV)) {
	LOCAL a,b
	TABLE hinf, minf, m_tau, h_tau FROM -150 TO 150 WITH 200
	a = 1.e-6*exp(-v/16.26)
	b = 1/(exp((-v+29.79)/10.)+1.)
	hinf = a/(a+b)
	a = 0.2*(-1.0*v+19.26)/(exp((-1.0*v+19.26)/10.0)-1.0)
	b = 0.009*exp(-v/22.03)
	minf = a/(a+b)
	a = 0.2*(-1.0*v+19.26)/(exp((-1.0*v+19.26)/10.0)-1.0)
	b = 0.009*exp(-v/22.03)
	m_tau = 1/(a+b)
	a = 1.e-5*exp(-v/16.26)
	b = 1/(exp((-v+29.79)/10.)+1.)
	h_tau = 1/(a+b)
}








