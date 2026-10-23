"""Run with PYTHONPATH=<build>/package and a current-format fixture argument."""
import pathlib
import sys
import tempfile
import unittest
from zspace.core import USDdoc
from zspace.io.slicing_usd import read_input, write_output

FIXTURE = pathlib.Path(sys.argv.pop(1)).resolve()


class USDAdapterTest(unittest.TestCase):
    def test_relative_reference_survives_edit(self):
        with tempfile.TemporaryDirectory() as folder:
            root = pathlib.Path(folder)
            (root/'mesh.usda').write_text('#usda 1.0\ndef Mesh "mesh" {\n point3f[] points = [(1,2,3)]\n}\n')
            (root/'root.usda').write_text('#usda 1.0\ndef "ref" (references = @mesh.usda@</mesh>) {}\n')
            doc = USDdoc(); doc.open(str(root/'root.usda'))
            try:
                self.assertEqual(doc.readPoints('/ref'), [[1, 2, 3]])
                doc.applyEdits({'customLayerData': {'test': {'type': 'int', 'value': 1}}})
                self.assertEqual(doc.readPoints('/ref'), [[1, 2, 3]])
            finally: doc.close()

    def test_roundtrip_and_overwrite(self):
        source = read_input(FIXTURE)
        with tempfile.TemporaryDirectory() as folder:
            output = pathlib.Path(folder) / 'output.usda'
            args = ([0, 0, 0, 10, 0, 0], [0, 0, 1, 0, 0, 1], [2, 2], [1, 1], 5)
            write_output(FIXTURE, output, *args)
            self.assertEqual(read_input(output), source)
            write_output(output, output, *args)
            doc = USDdoc(); doc.open(str(output))
            try:
                self.assertEqual(doc.readPrim('/TDP_gcode')['type'], 'Scope')
                self.assertEqual(doc.readAttribute('/TDP_gcode', 'printWidths')['value'], [2, 2])
                self.assertEqual(sum(p['path'] == '/TDP_gcode' for p in doc.readMetadata()['prims']), 1)
                self.assertIn('double printvolume = 20', doc.text())
                for kind, value in [('int', .5), ('float', 1e100)]:
                    with self.assertRaises(RuntimeError):
                        doc.applyEdits({'remove': ['/AG_block'], 'prims': [{'path': '/bad', 'type': 'Scope',
                            'attributes': {'value': {'type': kind, 'value': value}}}]})
                    self.assertEqual(doc.readPrim('/AG_block')['type'], 'Mesh')
            finally: doc.close()

    def test_typed_arrays_and_atomic_failure(self):
        doc = USDdoc(); doc.open(str(FIXTURE))
        try:
            values = {'ids': {'type': 'int[]', 'value': [-3, 4]},
                      'unsigned': {'type': 'uint[]', 'value': [4294967295]},
                      'weights': {'type': 'double[]', 'value': [0.125, 1.25]},
                      'empty': {'type': 'float[]', 'value': []},
                      'labels': {'type': 'string[]', 'value': ['a', 'b']}}
            doc.applyEdits({'prims': [{'path': '/typed', 'type': 'Scope', 'attributes': values}],
                            'customLayerData': {'blockid': {'type': 'int', 'value': 7}}})
            for name, value in values.items():
                actual = doc.readAttribute('/typed', name)
                self.assertEqual(actual['type'], value['type'])
                self.assertEqual(actual['value'], value['value'])
            self.assertEqual(doc.readCustomLayerData()['blockid']['value'], 7)
            before = doc.text()
            with self.assertRaises(RuntimeError):
                doc.applyEdits({'remove': ['/typed'], 'prims': [{'path': '/bad', 'type': 'Scope',
                    'attributes': {'points': {'type': 'point3f[]', 'value': [1, 2]}}}]})
            self.assertEqual(doc.text(), before)
            for old in ('metadata', 'prim', 'attribute', 'points', 'edit'):
                self.assertFalse(hasattr(doc, old))
        finally: doc.close()

    def test_units_and_optional_defaults(self):
        text = FIXTURE.read_text(encoding='utf-8')
        with tempfile.TemporaryDirectory() as folder:
            path = pathlib.Path(folder) / 'metres.usda'
            path.write_text(text.replace('metersPerUnit = 0.001', 'metersPerUnit = 1'), encoding='utf-8')
            mm, metres = read_input(FIXTURE), read_input(path)
            self.assertEqual(metres['mesh']['vertexPositions'], [v*1000 for v in mm['mesh']['vertexPositions']])
            self.assertEqual(metres['settings']['slotOffset'], mm['settings']['slotOffset'])
            path.write_text('\n'.join(line for line in text.splitlines() if 'slicer:SlotOffset' not in line), encoding='utf-8')
            self.assertEqual(read_input(path)['settings']['slotOffset'], 1)


if __name__ == '__main__': unittest.main()
