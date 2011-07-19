TITLE HH channel that includes both a sodium and a delayed rectifier channel 
: and accounts for sodium conductance attenuation
: Bartlett Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)
: Terrence Brannon-added attenuation 
: Yiota Poirazi-modified Kdr and Na threshold and time constants to make it more stable
: Yiota Poirazi-modified threshold for soma/axon spike initiation (threshold about -57 mV),
: USC Los Angeles 2000, poirazi@LNC.usc.edu
: This file is used only in soma and axon sections


NEURON {
	SUFFIX ch_leakhh
	NONSPECIFIC_CURRENT il
	RANGE gl, el
	RANGE ar2, vhalfs
	RANGE inf, fac, tau
	RANGE taus
	RANGE W
	GLOBAL taumin
        RANGE myi
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {                     :parameters that can be entered when function is called in cell-setup
 a0r = 0.0003 (ms)
        b0r = 0.0003 (ms)
        zetar = 12    
	zetas = 12   
        gmr = 0.2   
	ar2 = 1.0               :initialized parameter for location-dependent
                                :Na-conductance attenuation, "s", (ar=1 -> zero attenuation)
	taumin = 3   (ms)       :min activation time for "s" attenuation system
        vvs  = 2     (mV)       :slope for "s" attenuation system
        vhalfr = -60 (mV)       :half potential for "s" attenuation system
	W = 0.016    (/mV)      :this 1/61.5 mV
	gl = 0       (mho/cm2)
	el = -70.0   (mV)       :steady state 
	celsius = 34 (degC)
	v            (mV)
        dt
}

ASSIGNED {			:parameters needed to solve DE
	il (mA/cm2)
	myi (mA/cm2)
}

BREAKPOINT {
	il = gl*(v - el)               :leak current
	myi = il
}

INITIAL {			:initialize the following parameter using states()
	il = gl*(v - el)
	myi = il
}

