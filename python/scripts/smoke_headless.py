"""Run with PYTHONPATH pointing at build_headlessXXX/package."""
import sys
from pathlib import Path
from zspace import core
from zspace.toolset import slicer

v = core.zVector
mesh = core.zObjMesh()
fn = core.zFnMesh(mesh)
fn.create([v(0,0,0),v(100,0,0),v(100,100,0),v(0,100,0),
           v(0,0,50),v(100,0,50),v(100,100,50),v(0,100,50)],
          [4]*6, [0,3,2,1,4,5,6,7,0,1,5,4,1,2,6,5,2,3,7,6,3,0,4,7])
assert fn.numVertices() == 8 and fn.numPolygons() == 6
engine = slicer.zTs3DP()
engine.setInputMesh(mesh, [1,2,3,3,3,3], [-1,-1,0,0,0,0])
engine.setSlotAnchor(v(0,0,0))
settings = engine.getSettings()
settings.BoundaryWidth = [5.0]
settings.LayerSpacing = 10.0
settings.FieldCellSize = 2.0
settings.DoubleBoundary = False
assert engine.computeSlicer(False), engine.getStatus()
assert engine.getIsoMeshes()
print("PASS geometry and topology slicing", sys.version.split()[0])
ps = engine.getPrintSynSetting()
ps.ResampleLength = 5.0
ps.SeamOffset = 2.0
ps.SeamSize = 1.0
if engine.computePrintPath(False) and engine.computePostProcess(False):
    assert engine.getResultantPrintPath().widths
    print("PASS topology PrintSyn")
else:
    print("ALGORITHM_DIAGNOSTIC", engine.getStatus())
# Harmonic heat distances require intermediate side vertices.
points = [v(x,y,z) for z in (0,25,50) for x,y in ((0,0),(100,0),(100,100),(0,100))]
faces = [0,3,2,1,8,9,10,11]
for ring in (0,4):
    for i in range(4):
        j = (i+1)%4
        faces.extend([ring+i,ring+j,ring+j+4,ring+i+4])
fn.create(points, [4]*10, faces)
engine.setInputMesh(mesh, [1,2]+[3]*8, [-1,-1]+[0]*8)
engine.setSlicerMode(slicer.SlicerMode.Harmonic)
assert engine.computeSlicer(False), engine.getStatus()
print("PASS harmonic field and slicing")

if len(sys.argv) > 2 and not slicer.HAS_USD:
    print("SKIP USD checks: built without USD")
if len(sys.argv) > 2 and slicer.HAS_USD:
    source, target = sys.argv[1:3]
    inp = slicer.read_input_usd(source)
    graph = core.zObjGraph()
    core.zFnGraph(graph).create([v(0,0,0),v(1,0,0),v(2,1,0)], [0,1,1,2])
    path = slicer.PrintPath()
    path.graph = graph
    path.normals = [v(0,0,1)]*3
    path.widths = [2.0]*3
    path.heights = [1.0]*3
    slicer.export_printpath_usd(path, source, target, 1)
    slicer.export_printpath_usd(path, target, target, 1)
    reread = slicer.read_input_usd(target)
    assert reread.face_tags == inp.face_tags
    text = Path(target).read_text()
    assert text.count('def Scope "TDP_gcode"') == 1
    assert 'def BasisCurves "printpath"' not in text
    assert 'custom float[] printHeight' in text
    assert 'primvars:printHeight' not in text
    print("PASS USD read/export/overwrite")
