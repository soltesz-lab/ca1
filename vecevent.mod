:  Vector stream of events

NEURON {
	THREADSAFE
	ARTIFICIAL_CELL VecStim
	BBCOREPOINTER ptr
        RANGE xpos, ypos, zpos, gid, randi   
    }
    
PARAMETER {
    xpos = 0
    ypos = 0
    zpos = 0
    randi = 0
    gid = 0
}

ASSIGNED {
	index
	etime (ms)
	ptr
}


INITIAL {
	index = 0
	element()
	if (index > 0) {
		net_send(etime - t, 1)
	}
}

NET_RECEIVE (w) {
	if (flag == 1) {
		net_event(t)
		element()
		if (index > 0) {
			net_send(etime - t, 1)
		}
	}
}

FUNCTION is_art() {
	is_art=1
}

PROCEDURE position(a, b, c) { 
	xpos = a
	ypos = b
	zpos = c
}

DESTRUCTOR {
VERBATIM
	void* vv = (void*)(_p_ptr);  
        if (vv) {
		hoc_obj_unref(*vector_pobj(vv));
	}
ENDVERBATIM
}

VERBATIM
#include <stdint.h>
#if NRNBBCORE
#include "coreneuron/nrniv/ivocvect.h"
#endif
ENDVERBATIM

PROCEDURE element() {
VERBATIM	
  { void* vv; int i, size; double* px;
	i = (int)index;
	if (i >= 0) {
		vv = (void*)(_p_ptr);
		if (vv) {
			size = vector_capacity(vv);
			px = vector_vec(vv);
			if (i < size) {
				etime = px[i];
				index += 1.;
			}else{
				index = -1.;
			}
		}else{
			index = -1.;
		}
	}
  }
ENDVERBATIM
}

PROCEDURE play() {
VERBATIM
#if !NRNBBCORE
	void** pv;
	void* ptmp = NULL;
	if (ifarg(1)) {
		ptmp = vector_arg(1);
		hoc_obj_ref(*vector_pobj(ptmp));
	}
	pv = (void**)(&_p_ptr);
	if (*pv) {
		hoc_obj_unref(*vector_pobj(*pv));
	}
	*pv = ptmp;
#endif
ENDVERBATIM
}

VERBATIM
#if !NRNBBCORE
static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
        uint32_t dsize = 0;
        if (_p_ptr) {
          dsize = (uint32_t)vector_capacity(_p_ptr);
        }
        if (iArray) {
                uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
                void* vec = _p_ptr;
                ia[0] = dsize;

                double *da = dArray + *doffset;
                double *dv;
                if(dsize) {
                  dv = vector_vec(vec);
                }
                int iInt;
                for (iInt = 0; iInt < dsize; ++iInt) {
                  da[iInt] = dv[iInt];
                }
        }
        *ioffset += 1;
        *doffset += dsize;
}
#endif

static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
        assert(!_p_ptr);
        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
        int dsize = ia[0];
        *ioffset += 1;

        double *da = dArray + *doffset;
        _p_ptr = vector_new1(dsize);  /* works for dsize=0 */
        double *dv = vector_vec(_p_ptr);
        int iInt;
        for (iInt = 0; iInt < dsize; ++iInt)
        {
          dv[iInt] = da[iInt];
        }
        *doffset += dsize;
}
ENDVERBATIM

