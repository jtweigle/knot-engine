"""--------------------------------------------------+
| This program writes a Wavefront .OBJ file for a    |
| 3D Lissajous knot.                                 |
|                                                    |
| https://en.wikipedia.org/wiki/Lissajous_knot       |
|                                                    |
| For a minimal example, just run:                   |
|    $ python3 ./knot-engine.py                      |
|                                                    |
+--------------------------------------------------"""

import argparse # reason: parsing arguments
import sys      # reason: writing files, exiting
import math     # reason: trig functions

"""--------------------------------------------+
|      PARAMETERS FROM THE COMMAND LINE        |
| (what follows will have huge line widths...) |
+--------------------------------------------"""
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", default="knot.obj", help="specify the filename")
parser.add_argument("-a", "--a", type=int, default=3, help="inner coefficient for x (default 3)")
parser.add_argument("-b", "--b", type=int, default=4, help="inner coefficient for y (default 4)")
parser.add_argument("-c", "--c", type=int, default=7, help="inner coefficient for z (default 7)")
parser.add_argument("-wd", "--wrap-density", type=int, default=20, help="number of points around the tube that wraps the curve (default 20)")
parser.add_argument("-wr", "--wrap-radius", type=int, default=5, help="radius of the tube wrapping around the curve (default 5)")        
parser.add_argument("-cd", "--curve-density", type=int, default=200, help="number of points along the curve (default 200)")
parser.add_argument("-xs", "--x-scale", type=float, default=100.0, help="scale in the x direction (default 100)")
parser.add_argument("-ys", "--y-scale", type=float, default=100.0, help="scale in the y direction (default 100)")
parser.add_argument("-zs", "--z-scale", type=float, default=100.0, help="scale in the z direction (default 100)")
args = parser.parse_args()

"""---------------------------+
| ADDITIONAL GLOBAL VARIABLES |
+---------------------------"""
# I decided to move the phases here, since some phase inputs yielded
# really horrible results with my wrap function.
args.y_phase = 1
args.x_phase = 0.5
args.z_phase = 0.7
new_vertex = 1 # vertex keys are natural numbers
Faces = []     # list of lists of vertex keys

# Let's throw out some absurd inputs:
def bound(variable, lower, upper):
    return max(lower, min(upper, variable))
def nonzero(variable):
    return 1 if variable == 0 else variable
args.a = bound(args.a, -20, 20)
args.b = bound(args.b, -20, 20)
args.c = bound(args.c, -20, 20)
args.wrap_density = bound(args.wrap_density, 2, args.wrap_density)
args.wrap_radius = nonzero(bound(args.wrap_radius, -200, 200))
args.curve_density = bound(args.curve_density, 3, args.curve_density)

def main():
    global args
    filename = args.output
    file_object = open(filename, "w+")
    file_object.write("# This knot was brought to you by James Weigle's knot engine.\n\n")
    file_object.write("# Vertices:\n")
    draw_lissajous_curve(file_object)
    file_object.close()

"""--------------------------------------------+
| draw_lissajous_curve                         |
+----------------------------------------------+
| Draws a lissajous curve with the given       |
| parameters to file_object, which it assumes  |
| is open, and which it doesn't close.         |
+--------------------------------------------"""
def draw_lissajous_curve(file_object):
    global args
    curve_density = args.curve_density
    wrap_density  = args.wrap_density
    wrap_radius   = args.wrap_radius
    a = args.a
    b = args.b
    c = args.c
    x_scale = args.x_scale
    y_scale = args.y_scale
    z_scale = args.z_scale
    x_phase = args.x_phase
    y_phase = args.y_phase
    z_phase = args.z_phase
    
    step_length = math.pi / curve_density * 2.0
    t = 0.0
    vertex_tuples = []

    for i in range(curve_density + 1):
        t += step_length
        vertex_tuples.append(get_lissajous_vertex_with_frame(t))

    first_vertex = vertex_tuples[0]
    current_circle = make_new_circle(file_object, first_vertex,
                                     wrap_radius, wrap_density)

    for current_vertex in vertex_tuples[1:]:
        current_circle = wrap_segment(file_object,
                                      current_circle, current_vertex,
                     wrap_radius, wrap_density)

    write_faces(file_object)

"""------------------------------------------------------+
| get_lissajous_vertex_with_frame                        |
+--------------------------------------------------------+ 
| The frame in question is Tangent, Normal, Binormal:    |
| three orthonormal vectors from the curve at point t,   |
| which we can put together to make the rotation matrix  |
| that will orient each circle properly.                 | 
|                                                        |
| The function returns a tuple of the form               |
|   (vertex, tangent, normal, binormal)                  |
| where each entry is a 3D vector.                       |
+------------------------------------------------------"""
def get_lissajous_vertex_with_frame(t):
    global args
    a = args.a
    b = args.b
    c = args.c
    x_scale = args.x_scale
    y_scale = args.y_scale
    z_scale = args.z_scale
    x_phase = args.x_phase
    y_phase = args.y_phase
    z_phase = args.z_phase
 
    # These are used more than once, so we compute them once
    x_inside = a * t + math.pi * x_phase
    y_inside = b * t + math.pi * y_phase
    z_inside = c * t + math.pi * z_phase
    
    # point along the curve
    x = x_scale * math.sin(x_inside)
    y = y_scale * math.sin(y_inside)
    z = z_scale * math.sin(z_inside)
    vertex = (x,y,z)

    # tangent to the curve at that point
    # is just the first derivative
    x_tangent = a * x_scale * math.cos(x_inside)
    y_tangent = b * y_scale * math.cos(y_inside)
    z_tangent = c * z_scale * math.cos(z_inside)
    tangent   = normalize((x_tangent, y_tangent, z_tangent))

    # this is the derivative of the NORMALIZED tangent,
    # which is complex enough to compute that I pulled it out
    # as its own function
    normal = get_normal(t, x_inside, y_inside, z_inside,
                        (x_tangent, y_tangent, z_tangent))

    # binormal is the cross product of tangent and normal
    binormal = normalize(cross_product(tangent, normal))
    return (vertex, tangent, normal, binormal)

"""-----------------------------------------------+
| get_normal                                      |
+-------------------------------------------------+
| Implements the giant, hairy function, explained |
| below, to return the normal vector at the given |
| point t along the Lissajous curve.              |
+-----------------------------------------------"""
def get_normal(t, x_inside, y_inside, z_inside,
               pre_normalized_tangent):
    global args
    a = args.a
    b = args.b
    c = args.c

    x_scale = args.x_scale
    y_scale = args.y_scale
    z_scale = args.z_scale

    x_phase = args.x_phase
    y_phase = args.y_phase
    z_phase = args.z_phase

    """
    All right, this is going to look obscenely complicated,
    and that's because I DON'T KNOW VECTOR CALCULUS, haha.

    I got my results from this friendly
    youtube tutorial on the Frenet-Serret frame:
    https://www.youtube.com/watch?v=S3-m1cRt80k

    What we're doing here is taking the derivative
    of the NORMALIZED tangent vector---which means we can't
    simply take the second derivative of the lissajous function.
    """

    # First we get the *NORM* of the *UN-NORMALIZED* tangent:
    norm = get_norm(pre_normalized_tangent)

    # Now we *DO* have to take the derivative of each tangent component
    def derivative_of_tangent(coefficient, scale, inside):
        return -1 * (coefficient ** 2) * scale * math.sin(inside)

    derivative_of_x_term = derivative_of_tangent(a, x_scale, x_inside)
    derivative_of_y_term = derivative_of_tangent(b, y_scale, y_inside)
    derivative_of_z_term = derivative_of_tangent(c, z_scale, z_inside)

    # The first term in the derivative is defined as follows:
    first_term_x = derivative_of_x_term / norm
    first_term_y = derivative_of_y_term / norm
    first_term_z = derivative_of_z_term / norm

    # The second term is a gigantic, hairy, fearsome quotient,
    # but fortunately all three have the same common factor:

    x_term = pre_normalized_tangent[0]
    y_term = pre_normalized_tangent[1]
    z_term = pre_normalized_tangent[2]

    def common_factor_numerator_term(derivative_term, term):
        return 2 * derivative_term * term

    common_factor_x_term = common_factor_numerator_term(derivative_of_x_term,
                                                        x_term)
    common_factor_y_term = common_factor_numerator_term(derivative_of_y_term,
                                                        y_term)
    common_factor_z_term = common_factor_numerator_term(derivative_of_z_term,
                                                        z_term)
    common_factor_numerator = sum([common_factor_x_term,
                                   common_factor_y_term,
                                   common_factor_z_term])
    
    common_factor_denominator = 2 * (norm ** 3)
    common_factor = common_factor_numerator / common_factor_denominator

    final_x_term = first_term_x - x_term * common_factor
    final_y_term = first_term_y - y_term * common_factor
    final_z_term = first_term_z - z_term * common_factor
    
    return normalize((final_x_term, final_y_term, final_z_term))

"""------------------------------------------------------+
| cross-product                                          |
+--------------------------------------------------------+
| Returns the cross product of the two given 3d vectors, |
| as a float.                                            |
+------------------------------------------------------"""
def cross_product(v0, v1):
    s_0 = v0[1] * v1[2] - v0[2] * v1[1]
    s_1 = v0[2] * v1[0] - v0[0] * v1[2]
    s_2 = v0[0] * v1[1] - v0[1] * v1[0]
    return (s_0, s_1, s_2)

"""------------------------------------------------------+
| get_norm                                               |
+--------------------------------------------------------+
| Returns the norm of the given vector, as a float.      |
|                                                        |
| I pulled this out because it comes in handy when we're |
| computing the Frenet-Serret frame.                     |
+------------------------------------------------------"""
def get_norm(vector):
    norm_squared = sum([(vector[i] ** 2) for i in range(3)])
    norm = math.sqrt(norm_squared)
    return norm

"""------------------------------------------------------+
| normalize                                              |
+--------------------------------------------------------+
| Returns the normalized version of the given 3d vector. |
|                                                        |
| Technically returns it as a list of length 3,          |
| instead of a tuple, but in Python the difference       | 
| isn't very important.                                  |
|                                                        |
| If you'd like to sneer at me, do so now.               |
+------------------------------------------------------"""
def normalize(vector):
    norm = get_norm(vector)
    new_vector = [vector[i] / norm for i in range(3)]
    return new_vector

"""------------------------------------------------------+
| rotate_vertex                                          |
+--------------------------------------------------------+
| rotates the given vertex (vector, ugh! I'm using these |
| interchangeably...) by the matrix given by the basis   |
| (tangent, normal, binormal).                           |
+------------------------------------------------------"""
def rotate_vertex(vertex, tangent, normal, binormal):
    matrix = transpose_3d_square_matrix([tangent, normal, binormal])
    return matrix_vector_product_3d(matrix, vertex)

def transpose_3d_square_matrix(matrix):
    row_0 = [matrix[i][0] for i in range(3)]
    row_1 = [matrix[i][1] for i in range(3)]
    row_2 = [matrix[i][2] for i in range(3)]
    return [row_0, row_1, row_2]

"""-------------------------------------------------+
| matrix_vector_product_3d                          |
+---------------------------------------------------+
| In the general case, maybe it would be better to  |
| use numpy matrices with the built-in function.    |
|                                                   | 
| On the other hand, it was fun to write these      |
| linear algebra operations from scratch.           |
|                                                   |
| Plus, the constraints on the input parameters     |
| mean we'll never encounter a use case where       |
| we're doing an exponential number of these guys,  |
| so time complexity isn't life-or-death.           | 
|                                                   |
| Returns the matrix-vector product as a tuple.     |
+-------------------------------------------------"""
def matrix_vector_product_3d(matrix, vector):
    p0 = sum([(matrix[0])[i] * vector[i] for i in range(3)])
    p1 = sum([(matrix[1])[i] * vector[i] for i in range(3)])
    p2 = sum([(matrix[2])[i] * vector[i] for i in range(3)])
    return (p0, p1, p2)

"""------------------------------------------------+
| translate_vertex                                 |
+--------------------------------------------------+
| Moves the given vertex by the given coordinates, |
| returning the new vertex as a tuple.             |
+------------------------------------------------"""
def translate_vertex(vertex, coordinates):
    return (vertex[0] + coordinates[0],
            vertex[1] + coordinates[1],
            vertex[2] + coordinates[2])

"""-----------------------------------------------+
| wrap_segment                                    |
+-------------------------------------------------+
| This function makes a tube between the PREVIOUS |
| CIRCLE, and a circle around the NEXT VERTEX.    |
| Note that it creates the new circle, and NOT    |
| the previous circle.                            |
|                                                 |
| In the course of doing this, updates Faces      |
| and writes the vertices to file_object.         |
|                                                 |
| Returns the next circle, as a list of           |
| vertex keys, which are integers.                |
+-----------------------------------------------"""
def wrap_segment(file_object,
                 previous_circle,
                 current_vertex,
                 wrap_radius,
                 wrap_density):
    if len(previous_circle) != wrap_density:
        raise ValueError("previous_circle has different wrap_density "+
                         "than the current circle")
        sys.exit(1)
    
    # updates new_vertex and file_object:
    new_circle = make_new_circle(file_object, current_vertex,
                                 wrap_radius, wrap_density)

    pairs_of_vertices = list(zip(previous_circle, new_circle))
    # updates Faces:
    make_cyclic_ribbon(pairs_of_vertices) 
    return new_circle

"""-----------------------------------------------+
| make_new_circle                                 |
+-------------------------------------------------+
| This function                                   |
| current_vertex is assumed to be in the form     |
| (center_vertex, tangent, normal, binormal),     |
| as computed in get_lissajous_vertex_with_frame. |
|                                                 |
| You should know that the frame in question is   |
| a Frenet-Serret frame: basically, imagine an    |
| orthonormal basis moving on a roller coaster    |
| around the curve.                               |
|                                                 |
| The tangent, normal, and binormal become        |
| our little unit vectors, and that allows us     |
| to just put them together into the rotation     |
| matrix that moves the circle to the proper      |
| orientation at that moment on the curve.        |
| That's the advantage of this approach over      |
| computing the angle and trying to get the       |
| right orientation through three rotation        |
| matrices.                                       |
|                                                 |
| Note that this calling this function writes     |
| to file_object and updates new_vertex.          |
|                                                 |
| Returns the list of vertex keys                 |
| in the new circle.                              |
+-----------------------------------------------"""
def make_new_circle(file_object, current_vertex,
                    wrap_radius, wrap_density):
    wrap_radius  = max(wrap_radius, 1.0)
    wrap_density = max(wrap_density, 2.0)

    center_vertex = current_vertex[0]
    tangent       = current_vertex[1] 
    normal        = current_vertex[2] 
    binormal      = current_vertex[3] 

    # a list of tuples
    circle_coords_list = []
    vertex_key_list    = []
    step_length = 2.0 * math.pi / wrap_density
    t = 0.0
    for i in range(wrap_density):
        t += step_length
         
        x_circle = 0.0 
        y_circle = wrap_radius * math.cos(t)
        z_circle = wrap_radius * math.sin(t)

        unrotated_vertex = (x_circle, y_circle, z_circle)

        rotated_vertex = rotate_vertex(unrotated_vertex,
                                       tangent,
                                       normal,
                                       binormal)
        new_vertex = translate_vertex(rotated_vertex, center_vertex)
        circle_coords_list.append(new_vertex)
        
    for (x,y,z) in circle_coords_list:
        # This updates new_vertex:
        vertex_key_list.append(write_new_vertex(file_object, (x, y, z)))

    return vertex_key_list

"""---------------------------------------------+
| make_cyclic_ribbon                            |
+-----------------------------------------------+
| This updates Faces with what I'm calling a    |
| cyclic ribbon, meaning a tube that connects   |
| the given vertices.                           |
|                                               |
| The input is assumed to be a list of pairs    |
| of vertices (v1, v2), where v1 is next to v2. |
+---------------------------------------- ----"""
def make_cyclic_ribbon(pairs_of_vertices):
    rotated_once = rotate_list(pairs_of_vertices,1)
    pairs_of_pairs = list(zip(pairs_of_vertices, rotated_once))
    # This order was worked out by trial-and-error:
    # the goal is of course to connect them in such a way
    # that we end up with a rectangle, not a bow-tie
    for ((v00, v01), (v10, v11)) in pairs_of_pairs:
        Faces.append([v00, v10, v11, v01])
   
"""----------------------------------------------+
| write_new_vertex                               |
+------------------------------------------------+
| This writes the given vector to file_object as |
| a new vertex and returns the vertex key, which |
| is taken from the global variable new_vertex.  |
+----------------------------------------------"""
def write_new_vertex(file_object, vector):
    global new_vertex

    x = vector[0]
    y = vector[1]
    z = vector[2]
    
    x = round(x, 3)
    y = round(y, 3)
    z = round(z, 3)
    
    file_object.write(f"v {x} {y} {z}\n")

    current_vertex = new_vertex
    new_vertex += 1

    return current_vertex

"""-----------------------------------------------------+
| write_faces                                           |
+-------------------------------------------------------+
| The list of faces is stored in the global list Faces. |
| This writes them to file_object.                      |
+-----------------------------------------------------"""
def write_faces(file_object):        
    file_object.write("# Faces\n");
    for f in Faces:
        file_object.write("f");
        for v in f:
            file_object.write(f" {v}")
        file_object.write("\n") 

"""------------------------------------------------------+
| rotate_list                                            |
+--------------------------------------------------------+
| Not a rotation in the sense of a rotation matrix!      |
| For example, rotate_list([1,2,3,4], 2) is [3,4,1,2].   |
+------------------------------------------------------"""
def rotate_list(given_list, n):
    return given_list[n:] + given_list[:n]
 
# Oh yeah, one more thing:
if __name__ == "__main__":
    main()
