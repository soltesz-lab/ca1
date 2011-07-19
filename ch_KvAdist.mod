TITLE K-A channel from Klee Ficker and Heinemann
: modified to account for Dax A Current ----------
: M.Migliore Jun 1997
: modified by Poirazi on 10/2/00 according to Hoffman_etal97 
: to account for I_A distal (>100microns)
: (n) activation, (l) inactivation


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        dt (ms)
	v (mV)
        ek (mV)              : must be explicitely def. in hoc
	gmax=0.018 (mho/cm2)
        vhalfn= -1   (mV)
        vhalfl=-56   (mV)
}

NEURON {
	SUFFIX ch_KvAdist
	USEION k READ ek WRITE ik
        RANGE gmax,g
        GLOBAL ninf,linf,taul,taun
        RANGE myi
}

STATE {
	n
        l
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        linf
        taul
        taun
        g
	myi (mA/cm2)
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
	g = gmax*n^4*l
	ik = g*(v-ek)
	myi = ik
}

BREAKPOINT {
	SOLVE states
	g = gmax*n^4*l
	ik = g*(v-ek)
	myi = ik
}

FUNCTION alpn(v(mV)) {
  alpn = -0.01*(v+34.4)/(exp((v+34.4)/-21)-1)
}


FUNCTION betn(v(mV)) {
  betn = 0.01*(v+34.4)/(exp((v+34.4)/21)-1)
}

FUNCTION alpl(v(mV)) {
  alpl = -0.01*(v+58)/(exp((v+58)/8.2)-1)
}

FUNCTION betl(v(mV)) {
  betl = 0.01*(v+58)/(exp((v+58)/-8.2)-1)
}

LOCAL facn,facl
PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v)
        n = n + facn*(ninf - n)
        l = l + facl*(linf - l)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,b
        a = alpn(v)
        b = betn(v)
        ninf = a/(a + b)
        taun = 0.2
        facn = (1 - exp(-dt/taun))
        a = alpl(v)
        b = betl(v)
        linf = a/(a + b)
        
        if (v > -20) {
	   taul = 5 + 2.6*(v+20)/10
        } else {
	   taul = 5
        }
        facl = (1 - exp(-dt/taul))
}
