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
if engine.supportsSlicerMode(slicer.SlicerMode.Harmonic):
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

else:
    print("SKIP harmonic: built without TetGen")

if len(sys.argv) > 2:
    if not core.USDdoc.supported():
        print("SKIP USD checks: Core built without USD")
    else:
        from zspace.io.slicing_usd import read_input, write_output
        source, target = sys.argv[1:3]
        inp = read_input(source)
        args = ([0,0,0,1,0,0,2,1,0], [0,0,1]*3, [2.0]*3, [1.0]*3, 1)
        write_output(source, target, *args)
        write_output(target, target, *args)
        assert read_input(target) == inp
        text = Path(target).read_text()
        assert text.count('def Scope "TDP_gcode"') == 1
        assert 'custom float[] printHeight' in text
        assert 'primvars:printHeight' not in text
        print("PASS USD read/export/overwrite")
