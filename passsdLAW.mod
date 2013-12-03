TITLE passive membrane channel
                                                                                
UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
}
                                                                                
NEURON {
        SUFFIX passsdLAW
        NONSPECIFIC_CURRENT i
        RANGE gmax, erev, i
		RANGE myi
}
                                                                                
PARAMETER {
	v		(mV)
        gmax = .001        (mho/cm2)
        erev = -70      (mV)
}
                                                                                
ASSIGNED { i    (mA/cm2)
myi (mA/cm2)
}
                                                                                
BREAKPOINT {
        i = gmax*(v - erev)
		myi = i
}

