TITLE potassium leak channel
                                                                                
UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
}
                                                                                
NEURON {
        SUFFIX KleakaxonLAW
        USEION k READ ek WRITE ik
        RANGE gmax
		RANGE myi
}
                                                                                
PARAMETER {
	v		(mV)
        gmax = .001        (mho/cm2)
        ek = -70      (mV)
}
                                                                                
ASSIGNED { ik    (mA/cm2)
myi (mA/cm2)
}
                                                                                
BREAKPOINT {
        ik = gmax*(v - ek)
		myi = ik
}

