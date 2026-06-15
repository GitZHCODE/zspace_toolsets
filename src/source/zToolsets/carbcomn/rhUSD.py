import rhinoscriptsyntax as rs
import Rhino
import scriptcontext as sc

# --------------------------------------------------
# Sanitize USD names
# --------------------------------------------------

def sanitize(name):
    for c in " :/\\.-":
        name = name.replace(c, "_")
    return name

# --------------------------------------------------
# Select objects (Prompts user if nothing is selected)
# --------------------------------------------------

obj_ids = rs.GetObjects("Select geometry to export to USDA", preselect=True)

if not obj_ids:
    print("No objects selected. Export canceled.")

if obj_ids:
    # --------------------------------------------------
    # Output file (Brings up the browse dialog)
    # --------------------------------------------------
    
    filepath = rs.SaveFileName(
        title="Save USDA File",
        filter="USD ASCII (*.usda)|*.usda||",
        extension="usda"
    )

    if filepath:
        # --------------------------------------------------
        # Collect by layer
        # --------------------------------------------------
        
        layers = {}

        for obj_id in obj_ids:

            rh_obj = sc.doc.Objects.Find(obj_id)
            if rh_obj is None:
                continue

            layer_name = sc.doc.Layers[rh_obj.Attributes.LayerIndex].Name

            if layer_name not in layers:
                layers[layer_name] = {
                    "polylines": [],
                    "points": []
                }

            geo = rh_obj.Geometry

            # ---------------- Point ----------------

            if isinstance(geo, Rhino.Geometry.Point):
                layers[layer_name]["points"].append(geo.Location)
                continue

            # ---------------- PolylineCurve ----------------

            if isinstance(geo, Rhino.Geometry.PolylineCurve):
                ok, pl = geo.TryGetPolyline()
                if ok:
                    layers[layer_name]["polylines"].append(pl)
                continue

            # ---------------- General Curve ----------------

            if isinstance(geo, Rhino.Geometry.Curve):
                ok, pl = geo.TryGetPolyline()
                if ok:
                    layers[layer_name]["polylines"].append(pl)

        # --------------------------------------------------
        # Write USDA
        # --------------------------------------------------

        with open(filepath, "w") as f:

            f.write("#usda 1.0\n\n")
            f.write('def Xform "World"\n')
            f.write("{\n")

            for layer_name, data in layers.items():

                lname = sanitize(layer_name)

                f.write('    def Xform "{}"\n'.format(lname))
                f.write("    {\n")

                # ---------- Polylines ----------

                for i, pl in enumerate(data["polylines"]):

                    pts = list(pl)
                    closed = pl.IsClosed

                    # Remove duplicate last point if closed
                    if closed and len(pts) > 1:
                        if pts[0].DistanceTo(pts[-1]) < 1e-8:
                            pts = pts[:-1]

                    wrap = "periodic" if closed else "nonperiodic"

                    f.write('        def BasisCurves "Polyline_{}"\n'.format(i))
                    f.write("        {\n")
                    f.write('            uniform token type = "linear"\n')
                    f.write('            uniform token wrap = "{}"\n'.format(wrap))
                    f.write(
                        "            int[] curveVertexCounts = [{}]\n".format(len(pts))
                    )
                    f.write("            point3f[] points = [\n")

                    for j, p in enumerate(pts):
                        comma = "," if j < len(pts) - 1 else ""
                        f.write(
                            "                ({:.6f}, {:.6f}, {:.6f}){}\n".format(
                                p.X, p.Y, p.Z, comma
                            )
                        )

                    f.write("            ]\n")
                    
                    # ---- ADDED WIDTH FOR RENDERING ----
                    f.write('            float[] widths = [0.01]\n')
                    f.write('            uniform token interpolation = "constant"\n')
                    
                    f.write("        }\n\n")

                # ---------- Points ----------

                pts = data["points"]

                if pts:
                    f.write('        def Points "Points"\n')
                    f.write("        {\n")
                    f.write("            point3f[] points = [\n")

                    for j, p in enumerate(pts):
                        comma = "," if j < len(pts) - 1 else ""
                        f.write(
                            "                ({:.6f}, {:.6f}, {:.6f}){}\n".format(
                                p.X, p.Y, p.Z, comma
                            )
                        )

                    f.write("            ]\n")
                    
                    # ---- OPTIONAL: ADDED WIDTH FOR POINTS ----
                    f.write('            float[] widths = [0.01]\n')
                    
                    f.write("        }\n")

                f.write("    }\n\n")

            f.write("}\n")

        print("Export complete:")
        print(filepath)