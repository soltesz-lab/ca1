COMMENT

inhsyn

ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS inhsyn
	POINTER vpre
	RANGE on, tau1, tau2, gsbar, es, is, vprethresh, deadtime
	NONSPECIFIC_CURRENT is
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau1 (ms)
	tau2 (ms)
	gsbar 	(umho)
	es	(mV)
	v	(mV)
	vprethresh (mV)
	on (ms) : initialize to  < -deadtime
}

PARAMETER {
	deadtime (ms)
}

ASSIGNED {
	is (nA)
	gs (umho)
	vpre (mV)
}

INITIAL {
	on = -1e-10
}

BREAKPOINT {
	SOLVE getonset
	gs = gsbar * (exp(-(t-on)/tau2)-exp(-(t-on)/tau1))
	is = gs*(v-es)
}


PROCEDURE getonset() {
	:will crash if user hasn't set vpre with the connect statement */

	if (vpre > vprethresh && t > on + deadtime) {
		on = t
	}

	VERBATIM
	return 0;
	ENDVERBATIM
}
