#!/usr/bin/env python

__all__ = ['create_rocas', 'una_roca', 'una_roca_diag']

import numpy, random, math
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull
import sys



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

def una_roca_diagional_uniform(io, name, cname, roca_size=0.05, density=1, trans=None, tob=None):
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
                if inertia_type == 'diagonal':
                    una_roca_diagonal(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + random.random()*rate)
                elif inertia_type == 'diagonal_uniform':
                    una_roca_diagonal_uniform(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + random.random()*rate)
                elif inertia_type == 'rotated':
                    una_roca_rotated(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + random.random()*rate)
                elif inertia_type == 'cube_shape':
                    una_roca_cube_shape(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + random.random()*rate)
                else:
                    una_roca(io, name, cname, sizes[k], density, trans,
                             tob = n*rate + random.random()*rate)
                k += 1
