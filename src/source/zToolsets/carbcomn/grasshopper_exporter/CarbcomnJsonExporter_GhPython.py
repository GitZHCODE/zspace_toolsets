"""Grasshopper GhPython component script for Carbcomn block JSON export.

Inputs expected on the GhPython component:
    M                 Rhino.Geometry.Mesh
    folder            Output folder path
    block_id          Integer block id; writes blockMesh_<block_id>.json
    medial_start      Mesh vertex id at medial start
    medial_end        Mesh vertex id at medial end
    feature_strides   List[int], ordered around the regular wall section
    start_corner_vid  Mesh vertex id used as Carbcomn's StartCornerVID
    write             Boolean trigger

Outputs:
    path              Output JSON path
    message           Status / validation message

Notes:
    Carbcomn currently assumes non-planar regular wall blocks. This exporter
    still writes legacy compatibility keys IsPlanar and IsCorner as false.
"""

import json
import os
import Rhino.Geometry as rg


def _vertex_normal(mesh, vid):
    if mesh.Normals.Count == mesh.Vertices.Count:
        n = mesh.Normals[vid]
        return [float(n.X), float(n.Y), float(n.Z)]
    return [0.0, 0.0, 1.0]


def _face_normal(mesh, fid):
    if mesh.FaceNormals.Count == mesh.Faces.Count:
        n = mesh.FaceNormals[fid]
        return [float(n.X), float(n.Y), float(n.Z)]
    return [0.0, 0.0, 1.0]


def _face_vertices(face):
    if face.IsTriangle:
        return [face.A, face.B, face.C]
    return [face.A, face.B, face.C, face.D]


def _build_zspace_halfedge_json(mesh):
    """Build zFnMesh::fromJSON compatible half-edge JSON from a Rhino mesh."""
    mesh = mesh.DuplicateMesh()
    mesh.Faces.ConvertQuadsToTriangles()
    mesh.Normals.ComputeNormals()
    mesh.FaceNormals.ComputeFaceNormals()
    mesh.Compact()

    directed = {}
    faces = []
    face_he_ids = []
    he_records = []
    he_attrs = []

    def ensure_pair(a, b):
        key = (min(a, b), max(a, b))
        if key in directed:
            return directed[(a, b)]

        he0 = len(he_records)
        he1 = he0 + 1
        lo, hi = key
        directed[(lo, hi)] = he0
        directed[(hi, lo)] = he1
        he_records.append([-1, -1, hi, -1])
        he_records.append([-1, -1, lo, -1])
        he_attrs.append([0.0, 0.0, 0.0])
        he_attrs.append([0.0, 0.0, 0.0])
        return directed[(a, b)]

    for fid, face in enumerate(mesh.Faces):
        verts = _face_vertices(face)
        loop_hes = []
        for i, a in enumerate(verts):
            b = verts[(i + 1) % len(verts)]
            he = ensure_pair(a, b)
            if he_records[he][3] != -1:
                raise Exception("Non-manifold or duplicate oriented edge at face {0}, edge {1}->{2}".format(fid, a, b))
            he_records[he][2] = b
            he_records[he][3] = fid
            loop_hes.append(he)

        for i, he in enumerate(loop_hes):
            he_records[he][0] = loop_hes[(i - 1) % len(loop_hes)]
            he_records[he][1] = loop_hes[(i + 1) % len(loop_hes)]

        face_he_ids.append(loop_hes[0])
        faces.append(verts)

    vertex_he = [-1] * mesh.Vertices.Count
    for he_id, rec in enumerate(he_records):
        end_v = rec[2]
        if end_v >= 0:
            start_v = he_records[he_id ^ 1][2]
            if start_v >= 0 and vertex_he[start_v] == -1:
                vertex_he[start_v] = he_id

    vertex_attrs = []
    for i, p in enumerate(mesh.Vertices):
        n = _vertex_normal(mesh, i)
        vertex_attrs.append([float(p.X), float(p.Y), float(p.Z), n[0], n[1], n[2], 0.0, 0.0, 0.0])

    face_attrs = []
    for i in range(mesh.Faces.Count):
        n = _face_normal(mesh, i)
        face_attrs.append([n[0], n[1], n[2], 0.5, 0.5, 0.5])

    return {
        "Vertices": vertex_he,
        "Halfedges": he_records,
        "Faces": face_he_ids,
        "VertexAttributes": vertex_attrs,
        "HalfedgeAttributes": he_attrs,
        "FaceAttributes": face_attrs,
    }


def _validate_inputs():
    if M is None:
        return "Input M is required."
    if not isinstance(M, rg.Mesh):
        return "Input M must be a Rhino.Geometry.Mesh."
    if M.Vertices.Count == 0 or M.Faces.Count == 0:
        return "Input mesh has no vertices or faces."
    ids = [medial_start, medial_end, start_corner_vid]
    for vid in ids:
        if vid is None or int(vid) < 0 or int(vid) >= M.Vertices.Count:
            return "Invalid vertex id: {0}. Mesh vertex count is {1}.".format(vid, M.Vertices.Count)
    if feature_strides is None or len(feature_strides) == 0:
        return "feature_strides must be a non-empty list of positive integers."
    if len(feature_strides) % 2 == 0:
        return "feature_strides should have an odd count so Carbcomn can find the middle bridging stride."
    for stride in feature_strides:
        if int(stride) <= 0:
            return "All feature_strides must be positive integers."
    return None


path = None
message = "Set write=True to export."

if write:
    error = _validate_inputs()
    if error:
        message = error
    else:
        data = _build_zspace_halfedge_json(M)
        data["IsPlanar"] = False
        data["IsCorner"] = False
        data["MedialStartEnd"] = [int(medial_start), int(medial_end)]
        data["FeaturedNumStrides"] = [int(x) for x in feature_strides]
        data["StartCornerVID"] = int(start_corner_vid)

        if not os.path.isdir(folder):
            os.makedirs(folder)

        path = os.path.join(folder, "blockMesh_{0}.json".format(int(block_id)))
        with open(path, "w") as fp:
            json.dump(data, fp, indent=2)

        message = "Wrote {0}".format(path)
