import numpy as np
import pandas as pd
import allensdk.core.swc as swc
from lims_utils import get_swc_from_lims
from query import query_lims_for_layers
from morph import to_dict, dict_to_Morphology
from fiducials import get_coords, convert_coords_str, upright_angle
from geometry import line, intersection, find_translation, find_farthest, determine_mirror



def upright_nrn(specimen_id, oout, uout, error_dict={}):
    """ 
    Upright one cell 

    :param specimen_id: a cell specimen id
    :param oout: path to save original swc
    :param uout: path to save uprighted swc
    :param error_dict: dictionary to append message to if there's an error 
    :return: error_dict, unchanged if no error has occured 
    """

    try:
        print(specimen_id)

        ldf = query_lims_for_layers(specimen_id)

        try: soma_coords, pia_coords, wm_coords, layer_coords = get_coords(ldf, ["1", "2/3",  "4", "5", "6a", "6b"])
        except:
            print("ERROR: Couldn't load layers for {}".format(specimen_id))
            error_dict[specimen_id] = "Could not load layers"
            return error_dict

        try: swc_filename, swc_path = get_swc_from_lims(specimen_id)
        except TypeError:
            print("ERROR: Could not get swc from lims for ", specimen_id)
            error_dict[specimen_id] = "Could not get swc from lims"
            return error_dict
        
        swc_path = swc_path.replace('\\', '/')
        swc_path = swc_path.replace('/', '//', 1)
        morph = swc.read_swc(swc_path)
        morph.write(oout)

        if soma_coords is None:
            print("ERROR: No soma drawing for", specimen_id)
            error_dict[specimen_id] = "No soma drawing"
            return error_dict
        if pia_coords is None:
            print("ERROR: No 'pia' drawing for", specimen_id)
            error_dict[specimen_id] = "No 'pia' drawing"
            return error_dict
        if wm_coords is None:
            print("ERROR: No 'white matter' drawing for", specimen_id)
            error_dict[specimen_id] = "No 'white matter' drawing"
            return error_dict

        #Edit 'Pia' coords
        row = ldf[ldf.draw_type == 'Pia']
        res = row.res.values[0]  
        pcoords = row.poly_coords.values[0]
        plx, ply = convert_coords_str(pcoords)
        plx = plx * res
        ply = ply * res
        L1 = line([plx[0],ply[0]], [plx[-1], ply[-1]])

        #Add 'White Matter' coords
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
        soma_node = morph.compartment_list_by_type(1)[0]
        aff = [1., 0., 0., 0., 1., 0., 0., 0., 1., -soma_node["x"], -soma_node["y"], -soma_node["z"]]
        morph.apply_affine(aff)
        aff = [np.cos(theta), -np.sin(theta), 0., np.sin(theta), np.cos(theta), 0., 0., 0., 1., 0., -offset, 0.]
        morph.apply_affine(aff)

        print("\tsaving uprighted morph {}".format(specimen_id))
        morph.save(uout)

        flip = determine_mirror(lint, plx, ply, lx, ly)

        if flip:
            print("\tflipping morph {}".format(specimen_id))
            mdict = to_dict(uout)
            t = pd.DataFrame.from_dict(mdict).T
            t.x = t.x * -1
            tdict = t.to_dict(orient = 'index')
            tmorph = dict_to_Morphology(tdict)
            tmorph.save(uout)
    
        return error_dict
    
    except:
        print("ERROR: Unknown issue with cell {}", specimen_id)
        error_dict[specimen_id] = "Unknown issue with this cell"
        return error_dict
