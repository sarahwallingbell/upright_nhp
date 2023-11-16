import math
import numpy as np
from scipy.spatial.distance import euclidean


def _euclidean_distance(node1, node2):
    return euclidean(node1, node2)

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False
    
def find_translation(from_here, to_here):
    
    dx = to_here[0] - from_here[0]
    dy = to_here[1] - from_here[1]
    
    return dx, dy

def find_farthest(xs, ys, point):
    dist = 0
    opp = []
    for ix, x in enumerate(xs):
        y = ys[ix]
        this_dist = _euclidean_distance(point, [x, y])
        if this_dist > dist:
            dist = this_dist
            opp = [x, y]
            
    return opp

def _do_rotation(pts, angle):
    
    rot_x = []
    rot_y = []
    
    for p in pts:
#         print p
        rot = list(_rotate([0,0], p, angle))
        rot_x.append(rot[0])
        rot_y.append(rot[1])
        
    return rot_x, rot_y

def _unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def _angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = _unit_vector(v1)
    v2_u = _unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def _rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


def determine_mirror(lint, plx, ply, lx, ly):
    odx, ody = find_translation(lint, [0, 0])

    tp_x = plx + odx
    tp_y = ply + ody

    tm_x =  lx + odx
    tm_y = ly + ody

    opp = find_farthest(tp_x, tp_y, [0,0])

    
    opp.append(0)
    opp = tuple(opp)
    ##ccw angle
    ptheta= _angle_between((1, 0, 0), opp)#(tp_x[0], tp_y[0], 0))

    rotp_x, rotp_y = _do_rotation([[tp_x[0],tp_y[0]], [tp_x[-1], tp_y[-1]]], ptheta)
    rotm_x, rotm_y = _do_rotation([[tm_x[0],tm_y[0]], [tm_x[-1], tm_y[-1]]], ptheta)

    opp = find_farthest(rotp_x, rotp_y, [0,0])
    opp.append(0)
    opp = tuple(opp)
    ctheta= _angle_between((1, 0, 0), opp)#(tp_x[0], tp_y[0], 0))
    
    if ctheta != 0.0:
        
        rotp_x, rotp_y = _do_rotation([[tp_x[0],tp_y[0]], [tp_x[-1], tp_y[-1]]], -ptheta)
        rotm_x, rotm_y = _do_rotation([[tm_x[0],tm_y[0]], [tm_x[-1], tm_y[-1]]], -ptheta)

    mopp = find_farthest(rotm_x, rotm_y, [0,0])
    popp = find_farthest(rotp_x, rotp_y, [0,0])

    
    if mopp[1] < 0: #below the line
        if mopp[0] > popp[0]:
            return True
    else: #below line
        if mopp[0] < popp [0]:
            return True
        
    return False
