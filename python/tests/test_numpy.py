import gc
import numpy as np
from zspace import core
from zspace.toolset import slicer

p = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]], dtype=np.float64)
m = core.zObjMesh()
f = core.zFnMesh(m)
f.create_numpy(p[:, ::-1], np.array([4]), np.array([0,1,2,3]))
a,c,i = f.get_mesh_data_numpy()
np.testing.assert_allclose(a,p[:, ::-1])
np.testing.assert_array_equal(c,[4])
np.testing.assert_array_equal(i,[0,1,2,3])
a[:] = 99
assert not np.all(f.get_positions_numpy()==99)
f.set_positions_numpy(p)
for args in [(p,np.array([3]),np.array([0,1,2,3])),(p,np.array([4]),np.array([0,1,2,9])),(p,np.array([4.]),np.array([0,1,2,3]))]:
    try: f.create_numpy(*args)
    except ValueError: pass
    else: raise AssertionError('invalid mesh accepted')
np.testing.assert_allclose(f.get_positions_numpy(),p)
g=core.zObjGraph(); fn=core.zFnGraph(g)
fn.create_numpy(p,np.array([0,1,1,2,2,3]))
q,e=fn.get_graph_data_numpy()
np.testing.assert_allclose(q,p)
np.testing.assert_array_equal(e,[0,1,1,2,2,3])
path=slicer.PrintPath(); path.graph=g
path.normals=[core.zVector(0,0,1) for _ in p]
path.widths=[2.]*4; path.heights=[1.]*4
out=path.to_numpy()
del path,fn,g,f,m
gc.collect()
np.testing.assert_allclose(out["points"],p)
np.testing.assert_allclose(out["widths"],2.)
print("NumPy batch smoke checks passed")
