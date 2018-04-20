#!/usr/bin/env python

__all__ = ['create_rocas', 'una_roca', 'una_roca_diag']

import numpy, random, math
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull
import sys

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
=======
>>>>>>> a5185cd3106a8e282dd23509f5b96d10feb3eba7


def rotation_matrix_to_quaternion(A):


    N=numpy.zeros((4,4))
    #on-diagonal elements
    N[0][0] =  A[0][0]+A[1][1]+A[2][2]
    N[1][1] =  A[0][0]-A[1][1]-A[2][2]
    N[2][2] = -A[0][0]+A[1][1]-A[2][2]
    N[3][3] = -A[0][0]-A[1][1]+A[2][2]

    #off-diagonal elements
    N[0][1] = N[1][0] = A[2][1]-A[1][2]
    N[0][2] = N[2][0] = A[0][2]-A[2][0]
    N[0][3] = N[3][0] = A[1][0]-A[0][1]

    N[1][2] = N[2][1] = A[1][0]+A[0][1]
    N[1][3] = N[3][1] = A[0][2]+A[2][0]
    N[2][3] = N[3][2] = A[2][1]+A[1][2]

    l,v=numpy.linalg.eig(N)
    idx = l.argsort()[::-1]
    l = l[idx]
    v = v[:,idx]
    print('quaternion', v[:,0])

    return v[:,0]


def una_roca(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]

    def vert(id1, id2, a, b, c):
        return (a*math.cos(rd[id1])*math.cos(rd[id2]),
                b*math.sin(rd[id1])*math.cos(rd[id2]),
                c*math.sin(rd[id2]))

    vertices = [ vert( 0,  1,   1,  1,  1),
                 vert( 2,  3,   1, -1,  1),
                 vert( 4,  5,  -1,  1,  1),
                 vert( 6,  7,  -1, -1,  1),
                 vert( 8,  9,   1,  1, -1),
                 vert(10, 11,   1, -1, -1),
                 vert(12, 13,  -1,  1, -1),
                 vert(14, 15,  -1, -1, -1) ]

    scale = roca_size / max(numpy.array(vertices).max(axis=0)
                            - numpy.array(vertices).min(axis=0))

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.1*roca_size)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())

    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    # print('inertia:', inertia*density)


    io.add_object(name,
                 [Contactor(cname)],
                 translation=trans,
                 #velocity=veloci,
                 mass=volume*density,
                 time_of_birth=tob,
                 inertia=inertia*density)
def una_roca_cube_shape(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]

    def vert(id1, id2, a, b, c):
        return (a*math.cos(rd[id1])*math.cos(rd[id2]),
                b*math.sin(rd[id1])*math.cos(rd[id2]),
                c*math.sin(rd[id2]))

    vertices = [ vert( 0,  1,   1,  1,  1),
                 vert( 2,  3,   1, -1,  1),
                 vert( 4,  5,  -1,  1,  1),
                 vert( 6,  7,  -1, -1,  1),
                 vert( 8,  9,   1,  1, -1),
                 vert(10, 11,   1, -1, -1),
                 vert(12, 13,  -1,  1, -1),
                 vert(14, 15,  -1, -1, -1) ]

    scale = roca_size / max(numpy.array(vertices).max(axis=0)
                            - numpy.array(vertices).min(axis=0))

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # Definition of a polyhedron as a convex shape
    #io.add_convex_shape(cname, vertices, insideMargin=0.1*roca_size)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())

    test_occ_inertia(cname, vertices, density)


    print('roca_size', roca_size)
    print('geometric inertia:', inertia)
    print('volume:', volume)
    print('mass:', volume*density)
    print('inertia:', inertia*density)



    vertices = numpy.array([[  1,  1,  1],
                            [  1, -1,  1],
                            [ -1,  1,  1],
                            [ -1, -1,  1],
                            [  1,  1, -1],
                            [  1, -1, -1],
                            [ -1,  1, -1],
                            [ -1, -1, -1]]) * roca_size


    # ch = ConvexHull(vertices)
    # cm = ch.centroid()

    # # Definition of a polyhedron as a convex shape
    # #io.add_convex_shape(cname, vertices, insideMargin=0.1*roca_size)

    # # computation of inertia and volume
    # inertia,volume=ch.inertia(ch.centroid())

    # print('roca_size', roca_size)
    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    # print('inertia:', inertia*density*.001)


    io.add_convex_shape(cname, vertices, insideMargin=0.1*roca_size)


    io.add_object(name,
                 [Contactor(cname)],
                  translation=trans,
                  #velocity=veloci,
                  mass=volume*density,
                  time_of_birth=tob,
                  inertia=inertia*density)


def una_roca_diagonal(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]

    def vert(id1, id2, a, b, c):
        return (a*math.cos(rd[id1])*math.cos(rd[id2]),
                b*math.sin(rd[id1])*math.cos(rd[id2]),
                c*math.sin(rd[id2]))

    vertices = [ vert( 0,  1,   1,  1,  1),
                 vert( 2,  3,   1, -1,  1),
                 vert( 4,  5,  -1,  1,  1),
                 vert( 6,  7,  -1, -1,  1),
                 vert( 8,  9,   1,  1, -1),
                 vert(10, 11,   1, -1, -1),
                 vert(12, 13,  -1,  1, -1),
                 vert(14, 15,  -1, -1, -1) ]

    scale = roca_size / max(numpy.array(vertices).max(axis=0)
                            - numpy.array(vertices).min(axis=0))

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.1*roca_size)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())

    print('geometric inertia:', inertia)
    print('volume:', volume)
    print('mass:', volume*density)
    print('inertia:', inertia*density)

    l,v =  numpy.linalg.eig(inertia)
    print('eigen value', l)
    print('eigen vector', v)
    cond= numpy.min(l)/numpy.max(l)
    print('cond', cond)

    inertia_diag=numpy.zeros((3,3))
    inertia_diag[0,0] = inertia[0,0]
    inertia_diag[1,1] = inertia[1,1]
    inertia_diag[2,2] = inertia[1,1]

    inertia= inertia_diag
    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    print('inertia:', inertia*density)

    io.add_object(name,
                 [Contactor(cname)],
                 translation=trans,
                 #velocity=veloci,
                 mass=volume*density,
                 time_of_birth=tob,
                 inertia=inertia*density)

def una_roca_diagonal_uniform(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]

    def vert(id1, id2, a, b, c):
        return (a*math.cos(rd[id1])*math.cos(rd[id2]),
                b*math.sin(rd[id1])*math.cos(rd[id2]),
                c*math.sin(rd[id2]))

    vertices = [ vert( 0,  1,   1,  1,  1),
                 vert( 2,  3,   1, -1,  1),
                 vert( 4,  5,  -1,  1,  1),
                 vert( 6,  7,  -1, -1,  1),
                 vert( 8,  9,   1,  1, -1),
                 vert(10, 11,   1, -1, -1),
                 vert(12, 13,  -1,  1, -1),
                 vert(14, 15,  -1, -1, -1) ]

    scale = roca_size / max(numpy.array(vertices).max(axis=0)
                            - numpy.array(vertices).min(axis=0))

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.1*roca_size)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())

    print('geometric inertia:', inertia)
    print('volume:', volume)
    print('mass:', volume*density)
    print('inertia:', inertia*density)

    l,v =  numpy.linalg.eig(inertia)
    idx = l.argsort()[::-1]
    l = l[idx]
    v = v[:,idx]
    print('eigen value', l)
    print('eigen vector', v)
    cond= numpy.min(l)/numpy.max(l)
    print('cond', cond)

    inertia_diag=numpy.zeros((3,3))
    inertia_diag[0,0] = inertia[0,0]
    inertia_diag[1,1] = inertia[0,0]
    inertia_diag[2,2] = inertia[0,0]

    inertia= inertia_diag
    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    print('inertia:', inertia*density)


    io.add_object(name,
                 [Contactor(cname)],
                 translation=trans,
                 #velocity=veloci,
                 mass=volume*density,
                 time_of_birth=tob,
                 inertia=inertia*density)

def una_roca_rotated(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]

    def vert(id1, id2, a, b, c):
        return (a*math.cos(rd[id1])*math.cos(rd[id2]),
                b*math.sin(rd[id1])*math.cos(rd[id2]),
                c*math.sin(rd[id2]))

    vertices = [ vert( 0,  1,   1,  1,  1),
                 vert( 2,  3,   1, -1,  1),
                 vert( 4,  5,  -1,  1,  1),
                 vert( 6,  7,  -1, -1,  1),
                 vert( 8,  9,   1,  1, -1),
                 vert(10, 11,   1, -1, -1),
                 vert(12, 13,  -1,  1, -1),
                 vert(14, 15,  -1, -1, -1) ]

    scale = roca_size / max(numpy.array(vertices).max(axis=0)
                            - numpy.array(vertices).min(axis=0))

    print('scale', scale)
    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale

    ch = ConvexHull(vertices)
    cm = ch.centroid()


    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())


    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.1*roca_size)


    io.add_object(name,
                 [Contactor(cname)],
                 translation=trans,
                 #velocity=veloci,
                 mass=volume*density,
                 time_of_birth=tob,
                 inertia=inertia*density)

    l,v = numpy.linalg.eig(inertia)
    idx = l.argsort()[::-1]
    l = l[idx]
    v = v[:,idx]
    print('eigen values', l)
    print('eigen vectors', v)

    determinant = numpy.linalg.det(v)
    print('determinant', determinant, abs(determinant) -1.0)

    if abs(abs(determinant) -1.0) < 1e10:
        print('abs(determinant)+determinant',abs(determinant)-determinant)

        if numpy.sign(determinant) == 1:
            print('rotation matrix')
        else :
            print('perform a permutation')
            v[:,[0,1]] = v[:,[1,0]]
            #c[:, [0, 4]] = c[:, [4, 0]]
            print(v)
            determinant = numpy.linalg.det(v)
            print('determinant', determinant)
    else:
        print('not a normal matrix')
        exit()

    print('R^-1=R^T', numpy.linalg.norm(v.transpose()-numpy.linalg.inv(v)))

    print('vertices', vertices)

    vertices = numpy.dot(v.transpose(),vertices.transpose()).transpose()
    #vertices = numpy.dot(v,vertices.transpose()).transpose()
    q=rotation_matrix_to_quaternion(v)


    print('rotated vertices', vertices)
    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:])

    ch = ConvexHull(vertices)
    cm = ch.centroid()


    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())

    l,v = numpy.linalg.eig(inertia)
    print('eigen values', l)
    print('eigen vectors', v)
    print('determinant', numpy.linalg.det(v))
    print('R^-1=R^T', numpy.linalg.norm(v.transpose()-numpy.linalg.inv(v)))

    print('center of mass', cm)
    print('geometric inertia:', inertia)
    print('volume:', volume)
    print('mass:', volume*density)
    print('inertia:', inertia*density)

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname+'rot', vertices, insideMargin=0.1*roca_size)


    io.add_object(name+'rot',
                 [Contactor(cname+'rot')],
                  translation=trans,
                  orientation=q,
                 #velocity=veloci,
                 mass=volume*density,
                 time_of_birth=tob,
                 inertia=inertia*density)


def un_cubo(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
    # Definition of a cube as a convex shape

    vertices = numpy.array([[  1,  1,  1],
                            [  1, -1,  1],
                            [ -1,  1,  1],
                            [ -1, -1,  1],
                            [  1,  1, -1],
                            [  1, -1, -1],
                            [ -1,  1, -1],
                            [ -1, -1, -1]]) * roca_size

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = numpy.array(vertices)[:] - cm[:]

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())




    io.add_object(name,
                  [Contactor(cname)],
                  translation=trans,
                  #velocity=veloci,
                  mass=volume*density,
                  time_of_birth=tob,
                  inertia=inertia*density)


def create_rocas(io, n_layer=5, n_row=5, n_col=5, x_shift=3.0,
                 roca_size=0.05, top=0, rate=0.01, density=1,
                 distribution = ('uniform', 0.1), inertia_type=None):

    N = n_layer*n_row*n_col

    dist, rng = distribution

    if dist == 'uniform':
        sizes = numpy.random.uniform(low = roca_size - rng/2,
                                     high = roca_size + rng/2,
                                     size = N)
    elif dist == 'double':
        sizes = numpy.hstack(
            (numpy.random.normal(scale = rng*0.2,
                                 loc   = roca_size - rng/2,
                                 size  = N/2),
             numpy.random.normal(scale = rng*0.2,
                                 loc   = roca_size + rng/2,
                                 size  = N/2)))
        numpy.random.shuffle(sizes)
        # Gives a rock size distribution with two sizes of rocks, with
        # the mean between both at roca_size
    elif dist == 'exp':
        # In a normal distribution, 10- and 90-percentiles are loc +/- rng*1.28.
        # Let's arrange the 10- and 90-percentiles of the distribution
        # to be in the desired range.
        sizes = numpy.random.exponential(1, N)
        bottom = numpy.percentile(sizes, 10)
        top = numpy.percentile(sizes, 90)
        scale = (rng*1.28) / (top - bottom)
        sizes = (sizes - bottom)*scale + roca_size - rng/2*1.28

    k=0
    print('Creation of the rocks')
    for n in range(n_layer):
        for i in range(n_row):
            for j in range(n_col):
                # initial translation
                if (k%100 == 0):
                    print('.', end='', flush=True)
                trans = [(i-n_row/2.0)*x_shift*roca_size,
                         (j-n_col/2.0)*x_shift*roca_size,
                         top]
                name = 'rock'+str(n)+'_'+str(i)+'_'+str(j)
                cname = 'RockCS'+str(n)+'_'+str(i)+'_'+str(j)

                delay= random.random()*rate
                # print('delay', delay)
                # delay=0.0
                if inertia_type == 'diagonal':
                    una_roca_diagonal(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + delay)
                elif inertia_type == 'diagonal_uniform':
                    una_roca_diagonal_uniform(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + delay)
                elif inertia_type == 'rotated':
                    una_roca_rotated(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + delay)
                elif inertia_type == 'cube_shape':
                    una_roca_cube_shape(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + delay)
                else:
                    una_roca(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + delay)

                k += 1
