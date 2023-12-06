import numpy as np
from morph_utils.query import get_swc_from_lims
from neuron_morphology.swc_io import morphology_from_swc, morphology_to_swc
from neuron_morphology.transforms.affine_transform import AffineTransform
from morph_utils.modifications import normalize_position
from upright_nhp.query import query_lims_for_layers
from upright_nhp.fiducials import get_coords, convert_coords_str, upright_angle
from upright_nhp.geometry import line, intersection, find_translation, find_farthest, determine_mirror

def upright_nrn(specimen_id, morph=None, oout=None, uout=None, error_dict={}, print_info=False):
    """ 
    Upright one cell - save original to oout and upright to uout and return error_dict 
    Final upright orientation has dorsal on top and medial on the right. 

    :param specimen_id: a cell specimen id
    :param morph: a Morphology object to upright, otherwise looks for swc from lims
    :param oout: path to save original swc
    :param uout: path to save uprighted swc
    :param error_dict: dictionary to append message to if there's an error 
    :param print_info: whether or not to print info
    :return: morph, the uprighted morphology object or None if an error occured
    :return: error_dict, unchanged if no error has occured 
    """

    try:
        if print_info: print(specimen_id)

        ldf = query_lims_for_layers(specimen_id)

        try: soma_coords, pia_coords, wm_coords, layer_coords = get_coords(ldf, ["1", "2/3",  "4", "5", "6a", "6b"])
        except:
            if print_info: print("ERROR: Couldn't load layers for {}".format(specimen_id))
            error_dict[specimen_id] = "Could not load layers"
            return None, error_dict

        if morph is None:
            try: _, swc_path = get_swc_from_lims(specimen_id)
            except TypeError:
                if print_info: print("ERROR: Could not get swc from lims for ", specimen_id)
                error_dict[specimen_id] = "Could not get swc from lims"
                return None, error_dict
            
            swc_path = swc_path.replace('\\', '/')
            swc_path = swc_path.replace('/', '//', 1)
            morph = morphology_from_swc(swc_path)
        if oout: morphology_to_swc(morph, oout)

        if soma_coords is None:
            if print_info: print("ERROR: No soma drawing for", specimen_id)
            error_dict[specimen_id] = "No soma drawing"
            return None, error_dict
        if pia_coords is None:
            if print_info: print("ERROR: No 'pia' drawing for", specimen_id)
            error_dict[specimen_id] = "No 'pia' drawing"
            return None, error_dict
        if wm_coords is None:
            if print_info: print("ERROR: No 'white matter' drawing for", specimen_id)
            error_dict[specimen_id] = "No 'white matter' drawing"
            return None, error_dict

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

        #upright transform 
        theta, offset = upright_angle(layer_coords, soma_coords, pia_coords, wm_coords)
        theta += np.pi
        morph = normalize_position(morph)
        aff = [np.cos(theta), -np.sin(theta), 0., np.sin(theta), np.cos(theta), 0., 0., 0., 1., 0., -offset, 0.]
        upright_transform = AffineTransform.from_list(aff)
        morph = upright_transform.transform_morphology(morph)

        hflip = determine_mirror(lint, plx, ply, lx, ly)
        if hflip:
            # horizontal flip (so medial is right side)
            aff = [-1.,0.,0.,   0.,1.,0.,   0.,0.,1.,   0.,0.,0.]
            hflip_transform = AffineTransform.from_list(aff)
            morph = hflip_transform.transform_morphology(morph)

        # vertical flip (so dorsal is on top)
        aff = [1.,0.,0.,   0.,-1.,0.,   0.,0.,1.,   0.,0.,0.]
        vflip_transform = AffineTransform.from_list(aff)
        morph = vflip_transform.transform_morphology(morph)

        #save 
        if uout: morphology_to_swc(morph, uout)

        return morph, error_dict
    
    except:
        if print_info: print("ERROR: Unknown issue with cell {}", specimen_id)
        error_dict[specimen_id] = "Unknown issue with this cell"
        return None, error_dict