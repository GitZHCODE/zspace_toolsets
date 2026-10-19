"""Select three meshes: block (faces), bracing (edges only), slotpoint (one vertex).
Block materials: bottom/top/side_N. Bracing: INT EDGE attribute widthid, default 0.
All width values live in settings. USD loose edges use custom edgeVertexIndices.
World coordinates and widths use Blender units. No pxr required.
"""
import bpy
import json
import math
import re
from pathlib import Path

OUTPUT_PATH = "//ts3dp_input.usda"
# Edit this plaintext preset. Names match ts3dp_settings.json exactly.
SETTINGS_JSON = r'''
{
  "slicer": {
    "BoundaryWidth": [0.06, 0.03, 0.015],
    "BracingWidth": [0.04, 0.03, 0.02],
    "OverlapWidth": 0.0,
    "DoubleBoundary": true,
    "BoundarySlotOffset": 1.0,
    "BracingSlotEven": 0.45,
    "BracingSlotOdd": 0.55,
    "FieldScaleUp": 1.1,
    "FieldCellSize": 0.005,
    "LayerSpacing": 0.01,
    "HTraceStepFraction": 0.2,
    "HTraceMaxSteps": 100000,
    "HEqualiseIterations": 2,
    "HArapIterations": 10,
    "HUnrollStart": 0
  },
  "printSyn": {
    "ResampleLength": 0.01,
    "FeatureAngleTol": 120.0,
    "SeamOffset": 0.04,
    "SeamSize": 0.02
  }
}
'''
def number(v):
    v = float(v)
    if not math.isfinite(v):
        raise ValueError("Non-finite coordinate or setting")
    return format(v, ".17g")


def array(values):
    return "[" + ", ".join(str(v) for v in values) + "]"


def points(values):
    return array("(" + ", ".join(number(x) for x in p) + ")" for p in values)


def dictionary(data, indent=8):
    lines = []
    pad = " " * indent
    for key, value in data.items():
        if not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", key):
            raise ValueError("Invalid settings key: " + key)
        if isinstance(value, dict):
            lines += [f"{pad}dictionary {key} = {{", dictionary(value, indent+4), pad+"}"]
        else:
            if isinstance(value, bool):
                kind, text = "bool", str(value).lower()
            elif isinstance(value, int):
                kind, text = "int", str(value)
            elif isinstance(value, float):
                kind, text = "double", number(value)
            elif isinstance(value, list):
                kind, text = "double[]", array(number(v) for v in value)
            elif isinstance(value, str):
                kind, text = "string", json.dumps(value, ensure_ascii=True)
            else:
                raise ValueError("Unsupported setting: " + key)
            lines.append(f"{pad}{kind} {key} = {text}")
    return "\n".join(lines)


def roles():
    selected = list(bpy.context.selected_objects)
    if bpy.context.mode != "OBJECT" or len(selected) != 3:
        raise ValueError("Select exactly three objects in Object mode")
    result = {}
    for obj in selected:
        if obj.type != "MESH":
            raise ValueError("All three selected objects must be Mesh")
        if obj.name in ("block", "bracing", "slotpoint"):
            role = obj.name
        elif not obj.data.polygons and obj.data.edges:
            role = "bracing"
        elif len(obj.data.vertices) == 1:
            role = "slotpoint"
        elif obj.type == "MESH" and len(obj.data.polygons):
            role = "block"
        else:
            raise ValueError(f"Cannot identify {obj.name}; bracing must be an edge-only mesh")
        if role in result:
            raise ValueError("Multiple selected objects for " + role)
        result[role] = obj
    if set(result) != {"block", "bracing", "slotpoint"}:
        raise ValueError("Need block mesh, edge-only bracing mesh and single-vertex slotpoint mesh")
    return result


def block_text(obj, widths):
    if obj.type != "MESH":
        raise ValueError("block must be a mesh")
    evaluated = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh = evaluated.to_mesh()
    try:
        if not mesh.polygons:
            raise ValueError("block has no faces")
        tags, ids, counts, indices = [], [], [], []
        reverse = evaluated.matrix_world.to_3x3().determinant() < 0
        for face in mesh.polygons:
            material = mesh.materials[face.material_index] if face.material_index < len(mesh.materials) else None
            name = material.name if material else ""
            if name in ("bottom", "top"):
                tag, wid = (1 if name == "bottom" else 2), -1
            else:
                match = re.fullmatch(r"side_(\d+)", name)
                if not match:
                    raise ValueError(f"Face {face.index}: material '{name}' must be bottom/top/side_N")
                tag, wid = 3, int(match[1])
                if wid >= len(widths):
                    raise ValueError(f"Material {name} exceeds BoundaryWidth preset")
            tags.append(tag); ids.append(wid)
            vertices = list(face.vertices)
            counts.append(len(vertices))
            indices.extend(reversed(vertices) if reverse else vertices)
        return '\n'.join([
            'def Mesh "block" {',
            '    uniform token subdivisionScheme = "none"',
            '    point3f[] points = ' + points(evaluated.matrix_world @ v.co for v in mesh.vertices),
            '    int[] faceVertexCounts = ' + array(counts),
            '    int[] faceVertexIndices = ' + array(indices),
            '    custom int[] FaceTag = ' + array(tags),
            '    custom int[] widthid = ' + array(ids), '}'])
    finally:
        evaluated.to_mesh_clear()


def bracing_text(obj, widths):
    if any(m.show_viewport for m in obj.modifiers):
        raise ValueError("Apply bracing modifiers before export to preserve edge widthid")
    mesh = obj.data
    if mesh.polygons or not mesh.edges:
        raise ValueError("bracing must contain edges and no faces")
    attr = mesh.attributes.get("widthid")
    if attr is not None and (attr.domain != "EDGE" or attr.data_type != "INT"):
        raise ValueError("bracing widthid must be an Integer attribute on the Edge domain")
    ids = [attr.data[e.index].value if attr else 0 for e in mesh.edges]
    if any(i < 0 or i >= len(widths) for i in ids):
        raise ValueError("Edge widthid is outside settings BracingWidth")
    return '\n'.join([
        'def Mesh "bracing" {',
        '    uniform token subdivisionScheme = "none"',
        '    point3f[] points = ' + points(obj.matrix_world @ v.co for v in mesh.vertices),
        '    int[] faceVertexCounts = []',
        '    int[] faceVertexIndices = []',
        '    custom int[] edgeVertexIndices = ' + array(v for e in mesh.edges for v in e.vertices),
        '    custom int[] widthid = ' + array(ids), '}'])


def export_selected():
    settings = json.loads(SETTINGS_JSON)
    widths = settings["slicer"]["BoundaryWidth"]
    bracing_widths = settings["slicer"]["BracingWidth"]
    for table in (widths, bracing_widths):
        if not isinstance(table, list) or not table or any(not math.isfinite(float(w)) or float(w) <= 0 for w in table):
            raise ValueError("Width presets must be nonempty positive finite arrays")
    objects = roles()
    slot = objects["slotpoint"]
    if len(slot.data.vertices) != 1 or slot.data.edges or slot.data.polygons or slot.modifiers:
        raise ValueError("slotpoint must be a single-vertex mesh without edges, faces or modifiers")
    anchor = slot.matrix_world @ slot.data.vertices[0].co
    text = '\n'.join([
        '#usda 1.0', '(', '    upAxis = "Z"',
        '    metersPerUnit = ' + number(bpy.context.scene.unit_settings.scale_length),
        '    customLayerData = {', dictionary(settings), '    }', ')',
        block_text(objects['block'], widths), bracing_text(objects['bracing'], bracing_widths),
        'def Mesh "slotpoint" {', '    point3f[] points = ' + points([anchor]),
        '    int[] faceVertexCounts = []', '    int[] faceVertexIndices = []', '}', ''])
    path = Path(bpy.path.abspath(OUTPUT_PATH)).expanduser()
    if path.suffix.lower() != '.usda':
        raise ValueError("OUTPUT_PATH must end in .usda")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding='utf-8')
    print(f"[TS3DP] Exported {path}")


if __name__ == '__main__':
    export_selected()
