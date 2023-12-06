import numpy as np
from scipy.spatial.distance import euclidean


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
            coords_str = _check_rows(df.loc[df["draw_type"] == pt, ], pt, df["specimen_id"].values[0])
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

def _check_rows(rows, label, sp_id):
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

def upright_angle(layer_coords, soma_coords, pia_coords, wm_coords):
    soma_x, soma_y = soma_coords["x"], soma_coords["y"]
    avg_x = soma_x.mean()
    avg_y = soma_y.mean()
    soma_point = np.array([avg_x, avg_y])
    pia_proj = _project_to_polyline(pia_coords, soma_point)
    wm_proj = _project_to_polyline(wm_coords, soma_point)

    # Implied that we are getting the angle between vector connecting the projections and [0, 1]
    return np.arctan2(wm_proj[0] - pia_proj[0], wm_proj[1] - pia_proj[1]), euclidean(soma_point, pia_proj)

# def project_to_polyline_str(coords_str, target_point):
#     x, y = convert_coords_str(coords_str)
#     points = zip(x, y)
#     dists_projs = [_dist_proj_point_lineseg(target_point, np.array(q1), np.array(q2))
#                    for q1, q2 in zip(points[:-1], points[1:])]
#     min_idx = np.argmin(np.array([d[0] for d in dists_projs]))
#     return dists_projs[min_idx][1]


def _project_to_polyline(coords, target_point):
    x, y = coords["x"], coords["y"]
    points = list(map(np.array, zip(x, y)))
    dists_projs = [_dist_proj_point_lineseg(target_point, q1, q2)
                   for q1, q2 in zip(points[:-1], points[1:])]
    min_idx = np.argmin(np.array([d[0] for d in dists_projs]))

    # check if the closes point is the endpoint of the whole polyline
    # - if so, extend past the edge
    if np.allclose(dists_projs[min_idx][0], points[0]) or np.allclose(dists_projs[min_idx][0], points[-1]):
        return _dist_proj_point_lineseg(target_point, points[min_idx],
                                       points[min_idx + 1], clamp_to_segment=False)[1]
    else:
        return dists_projs[min_idx][1]

def _dist_proj_point_lineseg(p, q1, q2, clamp_to_segment=True):
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