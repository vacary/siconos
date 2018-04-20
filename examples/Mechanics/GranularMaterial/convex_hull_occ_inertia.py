
__all__ = ['test_occ_inertia']

import numpy


from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeSolid
from OCC.gp import gp_Pnt, gp_Ax2, gp_Dir


from siconos.io.mechanics_hdf5 import compute_inertia_and_center_of_mass

class assert_isdone(object):
    '''
    raises an assertion error when IsDone() returns false, with the error
    specified in error_statement
    '''
    def __init__(self, to_check, error_statement):
        self.to_check = to_check
        self.error_statement = error_statement

    def __enter__(self, ):
        if self.to_check.IsDone():
            pass
        else:
            raise AssertionError(self.error_statement)

    def __exit__(self, assertion_type, value, traceback):
        pass

def make_polygon(args, closed=False):

    poly = BRepBuilderAPI_MakePolygon()
    for pt in args:
        # support nested lists
        if isinstance(pt, list) or isinstance(pt, tuple):
            for i in pt:
                poly.Add(i)
        else:
            poly.Add(pt)
    if closed:
        poly.Close()
    poly.Build()

    with assert_isdone(poly, 'failed to produce wire'):
        result = poly.Wire()
        return result


from OCC.TopAbs import *
from OCC.TopoDS import topods, TopoDS_Shape
class ShapeToTopology(object):
    '''
    looks up the topology type and returns the corresponding topological entity
    '''
    def __init__(self):
        self.topoTypes = {TopAbs_VERTEX:      topods.Vertex,
                          TopAbs_EDGE:        topods.Edge,
                          TopAbs_FACE:        topods.Face,
                          TopAbs_WIRE:        topods.Wire,
                          TopAbs_SHELL:       topods.Shell,
                          TopAbs_SOLID:       topods.Solid,
                          TopAbs_COMPOUND:    topods.Compound,
                          TopAbs_COMPSOLID:   topods.CompSolid,
                          }

    def __call__(self, shape):
        if isinstance(shape, TopoDS_Shape):
            return self.topoTypes[shape.ShapeType()](shape)
        else:
            raise AttributeError('shape has not method `ShapeType`')

    def __getitem__(self, item):
        return self(item)



def sew_shapes(shapes, tolerance=0.001):
    sew = BRepBuilderAPI_Sewing(tolerance)
    for shp in shapes:
        if isinstance(shp, list):
            for i in shp:
                sew.Add(i)
        else:
            sew.Add(shp)
    sew.Perform()
    print("n degenerated shapes", sew.NbDegeneratedShapes())
    print("n deleted faces:", sew.NbDeletedFaces())
    print("n free edges", sew.NbFreeEdges())
    print("n multiple edges:", sew.NbMultipleEdges())
    result = ShapeToTopology()(sew.SewedShape())
    return result

def make_face(*args):
    face = BRepBuilderAPI_MakeFace(*args)
    with assert_isdone(face, 'failed to produce face'):
        result = face.Face()
        face.Delete()
        return result

def make_solid(*args):
    sld = BRepBuilderAPI_MakeSolid(*args)
    with assert_isdone(sld, 'failed to produce solid'):
        result = sld.Solid()
        sld.Delete()
        return result

def occ_make_face(a,b,c):

    p1 = gp_Pnt( a[0], a[1], a[2] )
    p2 = gp_Pnt( b[0], b[1], b[2] )
    p3 = gp_Pnt( c[0], c[1], c[2] )

    poly = make_polygon([p1, p2, p3], closed=True)
    print('poly', poly)

    face= make_face(poly)

    return face



def test_occ_inertia(cname, coordinates, density):

    from pyhull.convex_hull import ConvexHull
    hull = ConvexHull(coordinates)
    faces =[]
    for vertices in hull.vertices:
        a = coordinates[vertices[0]]
        b = coordinates[vertices[1]]
        c = coordinates[vertices[2]]
        faces.append(occ_make_face(a,b,c))

    solid = make_solid(sew_shapes(faces))
    print('solid', solid)
    # from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs
    # step_writer = STEPControl_Writer()
    # step_writer.Transfer(solid, STEPControl_AsIs)
    # step_writer.Write(cname+'.step')

    from OCC.GProp import GProp_GProps
    from OCC.BRepGProp import brepgprop_VolumeProperties
    system = GProp_GProps()

    brepgprop_VolumeProperties(solid, system)

    mass=  system.Mass()
    gp_mat= system.MatrixOfInertia()
    inertia_matrix = numpy.zeros((3,3))
    for i in range(0,3):
        for j in range(0,3):
            inertia_matrix[i,j]=  gp_mat.Value(i+1,j+1)
    print('mass', mass)
    print('inertia_matrix', inertia_matrix)
    system.Add(system, density)
    mass=  system.Mass()
    gp_mat= system.MatrixOfInertia()
    inertia_matrix = numpy.zeros((3,3))
    for i in range(0,3):
        for j in range(0,3):
            inertia_matrix[i,j]=  gp_mat.Value(i+1,j+1)
    print('mass', mass)
    print('inertia_matrix', inertia_matrix)
