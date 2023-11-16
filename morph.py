from allensdk.core.swc import Compartment, Morphology, NODE_ID, NODE_TYPE, NODE_X, NODE_Y, NODE_Z, NODE_R, NODE_PN


def dict_to_Morphology(neuron_dict):
    """
    Takes a neuron dictionary and converts to a Morphology
    """
    nodes = []
    for vals in neuron_dict.values():
        new_node = Compartment({
            NODE_ID: vals['id'],
            NODE_TYPE: vals['type'],
            NODE_X: vals['x'],
            NODE_Y: vals['y'],
            NODE_Z: vals['z'],
            NODE_R: vals['radius'],
            NODE_PN: vals['parent']}
        )
        nodes.append(new_node)
    return Morphology(compartment_list=nodes)#, strict_validation=False)

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