COMMENT

ppsyn

ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS ppsyn
	RANGE on, tau1, tau2, gsbar, es, is, deadtime
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
	on (ms) 
}

PARAMETER {
	deadtime (ms)
}

ASSIGNED {
	is (nA)
	gs (umho)
	}

BREAKPOINT {
	SOLVE getonset
gs=gsbar * on
	is = gs*(v-es)
}


PROCEDURE getonset() {
	if ( t < deadtime) {
		on = 0
	}
	if ( t> deadtime) {
	on = (exp(-(t-deadtime)/tau2)-exp(-(t-deadtime)/tau1))	}
				}
