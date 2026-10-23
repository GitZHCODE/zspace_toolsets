"""Slicer USD IO delegates exclusively to the native Toolsets implementation."""
from zspace import core
from zspace.toolset import slicer

_MAPPING = {
 'boundaryWidths':'BoundaryWidth','bracingWidth':'BracingWidth',
 'layerSpacing':'LayerSpacing','fieldCellSize':'FieldCellSize','overlapWidth':'OverlapWidth',
 'fieldScaleUp':'FieldScaleUp','slotOffset':'SlotOffset','bracingSlot':'BracingSlot',
 'doubleBoundary':'DoubleBoundary','traceStepFraction':'_HTraceStepFraction','traceMaxSteps':'_HTraceMaxSteps',
 'equaliseIterations':'_HEqualiseIterations','arapIterations':'_HArapIterations','unrollStart':'_HUnrollStart'}
_PRINT = {'resampleLength':'ResampleLength','featureAngle':'FeatureAngleTol','seamOffset':'SeamOffset','seamSize':'SeamSize'}

def read_input(path):
    native=slicer.importUSD(str(path))
    positions,counts,connects=core.zFnMesh(native.Block).get_mesh_data_numpy()
    vertices,edges=core.zFnGraph(native.Bracing).get_graph_data_numpy()
    settings={key:getattr(native.Settings,name) for key,name in _MAPPING.items()}
    settings.update({key:getattr(native.PrintSettings,name) for key,name in _PRINT.items()})
    settings['bracingCurveWidths']=native.BracingPrintWidths
    anchor=native.StartPoint
    return dict(_native_input=native, mesh=dict(vertexPositions=positions.ravel().tolist(),faceCounts=counts.tolist(),faceConnects=connects.tolist()),
                bracing=dict(vertexPositions=vertices.ravel().tolist(),edgeConnects=edges.ravel().tolist()),
                faceTags=native.FaceTags,faceWidthIds=native.FaceWidthIds,anchor=[anchor.x,anchor.y,anchor.z],settings=settings)

def write_output(input_path,output_path,result):
    slicer.exportUSD(result['_native_path'],str(input_path),str(output_path),result['layercount'],solver=result.get('_solver'))
