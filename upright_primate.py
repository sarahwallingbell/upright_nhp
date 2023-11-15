import numpy as np
import pandas as pd
import math
import os
import allensdk.core.swc as swc
# import psycopg2
from scipy.spatial.distance import euclidean
from lims_utils import query #, get_swc_from_lims

# CONNECTION_STRING = 'host=limsdb2 dbname=lims2 user=limsreader password=limsro'
NODE_ID = 'id'
NODE_TYPE = 'type'
NODE_X = 'x'
NODE_Y = 'y'
NODE_Z = 'z'
NODE_R = 'radius'
NODE_PN = 'parent'
SWC_COLUMNS = [NODE_ID, NODE_TYPE, NODE_X, NODE_Y, NODE_Z, NODE_R, NODE_PN]

# def query(sql, args):
#     conn = psycopg2.connect(CONNECTION_STRING)
#     cur = conn.cursor()

#     cur.execute(sql, args)
#     results = cur.fetchall()

#     cur.close()
#     conn.close()

#     return results

def mouse_layer_edges(avg_layer_depths):
    layer_edges = {
        "1": (0., avg_layer_depths["2/3"]),
        "2/3": (avg_layer_depths["2/3"], avg_layer_depths["4"]),
        "4": (avg_layer_depths["4"], avg_layer_depths["5"]),
        "5": (avg_layer_depths["5"], avg_layer_depths["6a"]),
        "6a": (avg_layer_depths["6a"], avg_layer_depths["6b"]),
        "6b": (avg_layer_depths["6b"], avg_layer_depths["wm"]),
    }

    return layer_edges

def query_lims_for_layers(specimen_id):
    sql = """
    SELECT sp.id as specimen_id, sp.name AS specimen, sp.cell_depth,
        imt.name AS image_type, agl.name AS drawing_layer, polygon.id AS polygon_id,
        bp.biospecimen_id,
        polygon.path, layer.mag, polygon.display_attributes, sc.resolution, struct.acronym
    FROM specimens sp JOIN specimens spp ON spp.id=sp.parent_id
    JOIN image_series iser ON iser.specimen_id=spp.id AND iser.type = 'FocalPlaneImageSeries' AND iser.is_stack = 'f'
    JOIN sub_images si ON si.image_series_id=iser.id
    JOIN avg_graphic_objects layer ON layer.sub_image_id=si.id
    JOIN avg_graphic_objects polygon ON polygon.parent_id=layer.id
    LEFT JOIN biospecimen_polygons bp ON polygon.id = bp.polygon_id
    JOIN images im ON im.id=si.image_id
    JOIN image_types imt ON imt.id=im.image_type_id
    JOIN scans sc ON sc.slide_id=im.slide_id
    LEFT JOIN structures struct ON struct.id = polygon.cortex_layer_id
    JOIN avg_group_labels agl ON layer.group_label_id=agl.id
    WHERE sp.id = %s
    ORDER BY 1, 4, 5, 6
    """

    results = query(sql, (specimen_id, ))
    df = pd.DataFrame(results, columns=["specimen_id", "specimen_name", "cell_depth",
                                         "img_type", "draw_type", "poly_id", "biospecimen_id", "poly_coords",
                                         "mag", "dispattr", "res", "layer_acronym"])

    # keep only draw types we are interested in
    used_draw_types = ["Pia", "White Matter", "Soma", "Cortical Layers"]
    df = df.loc[df["draw_type"].isin(used_draw_types), :].drop_duplicates(subset=["biospecimen_id", "poly_coords"])

    # keep soma, pia, and wm for specimen only
#     mask_out = df["draw_type"].isin(["Pia", "White Matter", "Soma"]).values & (df["biospecimen_id"].values != specimen_id)
#     df = df.loc[~mask_out, :]

    return df

def get_coords(df, layer_labels, layer_debug_flag=False,
        allow_partial=False):
    LAYER_ACRONYM_TRANSLATION = {
        "Layer1": "1",
        "Layer2/3": "2/3",
        "Layer2": "2",
        "Layer3": "3",
        "Layer4": "4",
        "Layer5": "5",
        "Layer6": "6",
        "Layer6a": "6a",
        "Layer6b": "6b",
    }
    base_types = ["Pia", "White Matter", "Soma"]
    path_types = ["Cortical Layers", "Pia", "White Matter"]

    path_strings = []
    for pt in base_types:
        try:
            coords_str = check_rows(df.loc[df["draw_type"] == pt, ], pt, df["specimen_id"].values[0])
            path_strings.append(coords_str)
        except ValueError:
            print("Missing {:s}".format(pt))
            return None, None, None, None

    pia_coords_str, wm_coords_str, soma_coords_str = path_strings

    res = df.loc[df["draw_type"] == "Pia", "res"].values[0]

    pia_x, pia_y = convert_coords_str(pia_coords_str)
    pia_coords = {"x": pia_x * res, "y": pia_y * res}
    wm_x, wm_y = convert_coords_str(wm_coords_str)
    wm_coords = {"x": wm_x * res, "y": wm_y * res}
    soma_x, soma_y = convert_coords_str(soma_coords_str)
    soma_coords = {"x": soma_x * res, "y": soma_y * res}

    layer_coords = {}
    for i, row in df.iterrows():
        if row["draw_type"] == "Cortical Layers":
            if row["layer_acronym"] not in LAYER_ACRONYM_TRANSLATION:
                continue

            layer_name = LAYER_ACRONYM_TRANSLATION[row["layer_acronym"]]
            coords_str = row["poly_coords"]
            res = row["res"]
            x, y = convert_coords_str(coords_str)
            layer_coords[layer_name] = {"x": x * res, "y": y * res}
    
    if len(layer_coords.keys()) < len(layer_labels) and not allow_partial:
        if layer_debug_flag:
            print("Not enough layers found (found {:d}). ID = {:d}".format(len(layer_coords.keys()), df["specimen_id"].values[0]))
        layer_coords = None

    return soma_coords, pia_coords, wm_coords, layer_coords

def check_rows(rows, label, sp_id):
    if len(rows) == 0:
        print("No drawing available")
        raise ValueError("No drawing available")
    elif len(rows) > 1:
        matched_rows = rows.loc[rows["biospecimen_id"] == sp_id, :]
        if len(matched_rows) > 1:
            print("Multiple drawings associated with same biospecimen_id")
            raise ValueError("Multiple drawings associated with same biospecimen_id")
        elif len(matched_rows) == 0:
            print("No matching rows")
            raise ValueError("No matching rows")
        path_string = matched_rows["poly_coords"].values[0]
    else:
        path_string = rows["poly_coords"].values[0]

    if len(path_string) == 0:
        print("Path string is empty")
        raise ValueError("Path string is empty")

    return path_string

def convert_coords_str(coords_str):
    vals = coords_str.split(',')
    x = np.array(vals[0::2], dtype=float)
    y = np.array(vals[1::2], dtype=float)
    return x, y

# def get_swc_from_lims(specimen_id):
#     conn = psycopg2.connect(CONNECTION_STRING)
#     cur = conn.cursor()

#     SQL = "SELECT f.filename, f.storage_directory FROM \
#      neuron_reconstructions n JOIN well_known_files f ON n.id = f.attachable_id \
#      AND n.specimen_id = %s AND n.manual AND NOT n.superseded AND f.well_known_file_type_id = 303941301"
#     cur.execute(SQL, (specimen_id,))
#     result = cur.fetchone()

#     if result is None:
#         raise("No SWC file found for specimen ID {}".format(specimen_id))

#     swc_filename = result[0]
#     swc_path = result[1] + result[0]
# #     print "SWC file: " + swc_path

#     cur.close()
#     conn.close()
#     return swc_filename, swc_path

def upright_angle(layer_coords, soma_coords, pia_coords, wm_coords):
    soma_x, soma_y = soma_coords["x"], soma_coords["y"]
    avg_x = soma_x.mean()
    avg_y = soma_y.mean()
    soma_point = np.array([avg_x, avg_y])
    pia_proj = project_to_polyline(pia_coords, soma_point)
    wm_proj = project_to_polyline(wm_coords, soma_point)

    # Implied that we are getting the angle between vector connecting the projections and [0, 1]
    return np.arctan2(wm_proj[0] - pia_proj[0], wm_proj[1] - pia_proj[1]), euclidean(soma_point, pia_proj)

def project_to_polyline_str(coords_str, target_point):
    x, y = convert_coords_str(coords_str)
    points = zip(x, y)
    dists_projs = [dist_proj_point_lineseg(target_point, np.array(q1), np.array(q2))
                   for q1, q2 in zip(points[:-1], points[1:])]
    min_idx = np.argmin(np.array([d[0] for d in dists_projs]))
    return dists_projs[min_idx][1]


def project_to_polyline(coords, target_point):
    x, y = coords["x"], coords["y"]
    points = list(map(np.array, zip(x, y)))
    dists_projs = [dist_proj_point_lineseg(target_point, q1, q2)
                   for q1, q2 in zip(points[:-1], points[1:])]
    min_idx = np.argmin(np.array([d[0] for d in dists_projs]))

    # check if the closes point is the endpoint of the whole polyline
    # - if so, extend past the edge
    if np.allclose(dists_projs[min_idx][0], points[0]) or np.allclose(dists_projs[min_idx][0], points[-1]):
        return dist_proj_point_lineseg(target_point, points[min_idx],
                                       points[min_idx + 1], clamp_to_segment=False)[1]
    else:
        return dists_projs[min_idx][1]

def dist_proj_point_lineseg(p, q1, q2, clamp_to_segment=True):
    # based on c code from http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    l2 = euclidean(q1, q2) ** 2
    if l2 == 0:
        return euclidean(p, q1), q1 # q1 == q2 case
    if clamp_to_segment:
        t = max(0, min(1, np.dot(p - q1, q2 - q1) / l2))
    else:
        t = np.dot(p - q1, q2 - q1) / l2
    proj = q1 + t * (q2 - q1)
    return euclidean(p, proj), proj


def euclidean_distance(node1, node2):
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
        this_dist = euclidean_distance(point, [x, y])
        if this_dist > dist:
            dist = this_dist
            opp = [x, y]
            
    return opp

def do_rotation(pts, angle):
    
    rot_x = []
    rot_y = []
    
    for p in pts:
        rot = list(rotate([0,0], p, angle))
        rot_x.append(rot[0])
        rot_y.append(rot[1])
        
    return rot_x, rot_y

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def rotate(origin, point, angle):
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
    flipX = False
    # flipY = False

    odx, ody = find_translation(lint, [0, 0])

    tp_x = plx + odx
    tp_y = ply + ody

    tm_x =  lx + odx
    tm_y = ly + ody

    opp = find_farthest(tp_x, tp_y, [0,0])

    opp.append(0)
    opp = tuple(opp)
    ##ccw angle
    ptheta= angle_between((1, 0, 0), opp)#(tp_x[0], tp_y[0], 0))

    rotp_x, rotp_y = do_rotation([[tp_x[0],tp_y[0]], [tp_x[-1], tp_y[-1]]], ptheta)
    rotm_x, rotm_y = do_rotation([[tm_x[0],tm_y[0]], [tm_x[-1], tm_y[-1]]], ptheta)

    opp = find_farthest(rotp_x, rotp_y, [0,0])
    opp.append(0)
    opp = tuple(opp)
    ctheta= angle_between((1, 0, 0), opp)#(tp_x[0], tp_y[0], 0))
    
    if ctheta != 0.0:
        rotp_x, rotp_y = do_rotation([[tp_x[0],tp_y[0]], [tp_x[-1], tp_y[-1]]], -ptheta)
        rotm_x, rotm_y = do_rotation([[tm_x[0],tm_y[0]], [tm_x[-1], tm_y[-1]]], -ptheta)

    mopp = find_farthest(rotm_x, rotm_y, [0,0])
    popp = find_farthest(rotp_x, rotp_y, [0,0])

    if mopp[1] < 0: #below the line
        if mopp[0] > popp[0]:
            # return True
            flipX = True
    else: #below line
        if mopp[0] < popp [0]:
            # return True
            flipX = True
    # return False
    return flipX

def dict_to_Morphology(neuron_dict):
    """
    Takes a neuron dictionary and converts to a Morphology
    """
    nodes = []
    for vals in neuron_dict.values():
        new_node = swc.Compartment({
            NODE_ID: vals['id'],
            NODE_TYPE: vals['type'],
            NODE_X: vals['x'],
            NODE_Y: vals['y'],
            NODE_Z: vals['z'],
            NODE_R: vals['radius'],
            NODE_PN: vals['parent']}
        )
        nodes.append(new_node)
    return swc.Morphology(compartment_list=nodes)#, strict_validation=False)

def morph_to_dict(morph):
    """
    Takes a neuron Morphology object (loaded from allensdk) and converts to a python dict 
    """
    nodes = {node['id'] : node for node in morph.compartment_list}

    return nodes


def to_dict(swc_file):
    nodes = {}
    with open(swc_file, "r") as f:
        for line in f:
            if line.lstrip().startswith('#'):
                continue
            toks = line.split()
            node_dict = {
                'id' : int(toks[0]),
                'type' : int(toks[1]),
                'x' : float(toks[2]),
                'y' : float(toks[3]),
                'z' : float(toks[4]),
                'radius' : float(toks[5]),
                'parent' : int(toks[6].rstrip())
            }
            nodes[int(toks[0])] = node_dict
    return nodes



#used in subcortical viewers
def upright_subcortical_morpho(specimen_id, swc_path, upright_swc_dir, print_on=False):
    try:
        ldf = query_lims_for_layers(specimen_id)
        layer_list = ["1", "2/3",  "4", "5", "6a", "6b"]
        soma_coords, pia_coords, wm_coords, layer_coords = get_coords(ldf, layer_list)
    except:
        if print_on: print("Warning: No layers found for specimen id {}".format(specimen_id))
        return None

    #ensure we have soma drawing for this cell 
    if soma_coords is None:
        if print_on: print("ERROR: No soma drawing for", specimen_id)
        return None

    #filename path for saving uprighted swc
    upright_swc_path = os.path.join(upright_swc_dir, '{}_upright.swc'.format(specimen_id))

    #load swc
    nrn=swc.read_swc(swc_path)

    try:
        #Edit WM coords
        row = ldf[ldf.draw_type == 'Pia']
        res = row.res.values[0]  
        pcoords = row.poly_coords.values[0]
        plx, ply = convert_coords_str(pcoords)
        plx = plx * res
        ply = ply * res
        L1 = line([plx[0],ply[0]], [plx[-1], ply[-1]])

        #Add White Matter
        row = ldf[ldf.draw_type == 'White Matter']
        res = row.res.values[0]  
        wcoords = row.poly_coords.values[0]
        lx, ly = convert_coords_str(wcoords)
        lx = lx * res
        ly = ly * res
        L2 = line([lx[0],ly[0]], [lx[-1], ly[-1]])

        lint = list(intersection(L1, L2))

        opp = find_farthest(lx, ly, lint)

        dx, dy = find_translation(lint, opp)
        new_lx = np.asarray(plx + dx)
        new_ly = np.asarray(ply + dy)

        wm_coords['x'] = pia_coords['x']
        wm_coords['y'] = pia_coords['y']

        pia_coords['x'] = new_lx
        pia_coords['y'] = new_ly


        theta, offset = upright_angle(layer_coords, soma_coords, pia_coords, wm_coords)
        theta += np.pi

        #upright
        soma_node = nrn.compartment_list_by_type(1)[0]
        aff = [1., 0., 0., 0., 1., 0., 0., 0., 1., -soma_node["x"], -soma_node["y"], -soma_node["z"]]
        nrn.apply_affine(aff)
        aff = [np.cos(theta), -np.sin(theta), 0., np.sin(theta), np.cos(theta), 0., 0., 0., 1., 0., -offset, 0.]
        nrn.apply_affine(aff)
        mdict = morph_to_dict(nrn) #to dict 

        flip = determine_mirror(lint, plx, ply, lx, ly)
        if flip:
            if print_on: print("\tflipping morph {}".format(specimen_id))                
            t = pd.DataFrame.from_dict(mdict).T
            t.x = t.x * -1
            tdict = t.to_dict(orient = 'index')
            tmorph = dict_to_Morphology(tdict)
            tmorph.save(upright_swc_path)
        else:
            mmorph = dict_to_Morphology(mdict)
            mmorph.save(upright_swc_path)
            
        if print_on: print("\tsaved uprighted swc: {}".format(upright_swc_path))

        return upright_swc_path
    
    except: return None