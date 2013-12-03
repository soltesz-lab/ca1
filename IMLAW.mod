UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX IMLAW
        USEION k READ ek WRITE ik
        RANGE gmax,ik, vhalf1, vhalf2, k1, k2, c1, c2
        GLOBAL minf, mtau
		RANGE myi
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gmax = 0.001 (mho/cm2)		
	ek
	mmin = 7	(ms)
	vhalf1 = -63
	vhalf2 = -63
	k1 = 15
	k2 = 15
	c1 = 0.003
	c2 = 0.003
}
 
STATE {
        m
}
 
ASSIGNED {
	ik
	minf 
	mtau	(ms)
	myi (mA/cm2)
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        ik = gmax*m*(v - ek)
		myi = ik
}
 
INITIAL {
	rates(v)
	m = minf
}

DERIVATIVE state { :Computes state variable h at current v and dt.
	rates(v)
	m' = (minf - m)/mtau
}

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
:        TABLE rinf, rexp, tau_r DEPEND dt, p FROM -200 TO 100 WITH 300
	LOCAL a, b
	minf = 1/(1 + exp(-(v+27(mV))/7(mV)))
	a = c1/exp(-(v-vhalf1)/k1)
	b = c2/exp((v-vhalf2)/k2)
	mtau = 1/(a+b)
	if (mtau<mmin) {mtau = mmin}
}
 
UNITSON

