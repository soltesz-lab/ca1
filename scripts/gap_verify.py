'''
Gap junction consistency tests.

0) Every HalfGap.vgap receives information from some source.
1) After a transfer, every HalfGap.vgap value has a valid encoded source rank
2) All the HalfGap.vgap that have the same source, have the same value.
3) For all the  HalfGap.vgap value of encoded source sid and rank, that
     source sid actually exists on the source rank.

The tests are prepared by calling
  target_var(halfgap_instance, sid)
  source_var(x, sid)
for every call of the corresponding pc.target_var and pc.source_var calls.
Note that x is location in the currently accessed section. If the real
gaps are setup using Python, it may be more convenient to call the alternative
  src_var(seg, sid)

The tests are carried out by calling
  setup_transfer_verify()
just before stdinit() (after setup_transfer and set_maxstep).

Note, this presumes the name, HalfGap.vgap . If the POINT_PROCESS is not
HalfGap or the relevant range variable is not vgap, change appropriately
below.


'''

from neuron import h

pc = h.ParallelContext()
rank = int(pc.id())
nhost = int(pc.nhost())

sid2src = {}  # sid key specifies only one source (globally but only tested locally)
src2sids = {} # possibly several sids associated with same src
sid2tars = {} # possibly several targets associated with same sid
tar2sid = {}  # target key has only one sid

def target_var(g, id):
  sid = int(id)
  if g in tar2sid:
    print ('%d %s with sid %d reregistered with sid %d'%(rank,g.hname(),tar2sid[g],id))
    raise RuntimeError
  tar2sid[g] = sid
  if sid not in sid2tars:
    sid2tars[sid] = []
  sid2tars[sid].append(g)

def source_var(x, id):
  src_var(h.cas()(x), id)

def src_var(seg, id):
  sid = int(id)
  if sid in sid2src:
    s1 = name(seg)
    s2 = name(sid2src[sid])
    print ('%d source var sid %d for %s already in use for %s'%(rank,id,s1,s2))
    raise RuntimeError
  sid2src[sid] = seg
  if seg not in src2sids:
    src2sids[seg] = []
  src2sids[seg].append(sid)

def name(seg):
  return '%s(%g)' % (seg.sec.name(), seg.x)

def rank2frac(r):
  return float(r)/(2.**18)

def frac2rank(f):
  x = abs(f)
  frac = x - float(int(x))
  xi = frac*(2.**18)
  if xi != float(int(xi)):
    print ('%d frac2rank f=%g x=%g xi=%g'%(rank, f, x, xi))
    raise RunTimeError
  return int(xi)

def decode_vgap_val(v):
  r = frac2rank(v)
  sid = int(v)
  return (r, sid)

def pr(data, label):
  pc.barrier()
  for i in range(nhost):
    if i == rank:
      print ('%d %s'%(rank,label))
      print (data)
      pc.barrier()

def setup_transfer_verify():
  if rank == 0 : print ('setup_transfer_verify')
  h.finitialize(-65)
  frac = rank2frac(rank)
  for sec in h.allsec():
      sec.v = -1 - frac
  for seg in src2sids:
    seg.v = src2sids[seg][0] + frac
  for tar in tar2sid:
    tar.vgap = -2.0

  # do the parallel transfer
  h.finitialize()
  # note that FInitializeHandlers may destroy v but the transfer takes place
  # before any callback except type 3. So all the targets should have proper
  # value of firstsrcsid.rank

  # test 0. No transferred values are negative
  for tar in tar2sid:
    if tar.vgap < 0.0:
      print ('%d %s sid=%d vgap=%g'%(rank, tar.hname(), tar2sid[tar], tar.vgap))
      raise RuntimeError
  # test 1. All the transferred values make sense in terms of rank
  for tar in tar2sid:
    if frac2rank(tar.vgap) >= nhost:
      print ('%d %s has invalid rank code %d with value %g'%(rank, tar.hname(), frac2rank(tar.vgap), tar.vgap))
      raise RuntimeError

  # test 2. All the targets for a given sid should have the same value
  for sid in sid2tars:
    x0 = sid2tars[sid][0].vgap
    for tar in sid2tars[sid]:
      if tar.vgap != x0:
        print ('%d %s %g != %g %s'%(rank, sid2tars[sid][0].hname(),x0, tar.vgap, tar.hname()))
        raise RuntimeError

  # test 3. Send sid and vgap's sid back to rank it came from and verify on
  # the source rank that those sid's exist and have the same source
  # Because of test 2, only need to test first target of an sid
  data = [None]*nhost
  for sid, segs in sid2tars.items():
    r,ssid = decode_vgap_val(segs[0].vgap)
    if data[r] == None:
      data[r] = []
    data[r].append((sid, ssid))

  #pr(data, 'source')
  data = pc.py_alltoall(data)
  #pr(data, 'destination')

  for r,x in enumerate(data):
    if x:
      for pair in x:
        for sid in pair:
          if sid not in sid2src:
            print ('%d target sid %d from %d not associated with source here'%(rank,sid,r))
            raise RuntimeError
