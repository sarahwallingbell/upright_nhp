import pandas as pd
import os
import allensdk.internal.core.lims_utilities as lu
from functools import partial

def default_query_engine():
    """Get Postgres query engine with environmental variable parameters"""

    return partial(
        lu.query,
        # host=os.getenv("LIMS_HOST"),
        # port=5432,
        # database=os.getenv("LIMS_DBNAME"),
        # user=os.getenv("LIMS_USER"),
        # password=os.getenv("LIMS_PASSWORD")
    )

def query_lims_for_layers(specimen_id, query_engine=None):
       
    if query_engine is None: query_engine = default_query_engine()

    sql = """
    SELECT sp.id AS specimen_id, 
        sp.name AS specimen_name, 
        sp.cell_depth,
        imt.name AS img_type, 
        agl.name AS draw_type, 
        polygon.id AS poly_id,
        bp.biospecimen_id, 
        polygon.path AS poly_coords, 
        layer.mag, 
        polygon.display_attributes AS dispattr, 
        sc.resolution AS res, 
        struct.acronym AS layer_acronym
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
    WHERE sp.id = {}
    ORDER BY 1, 4, 5, 6
    """.format(specimen_id)

    results = query_engine(sql) 
    df = pd.DataFrame(results)

    # keep only draw types we are interested in
    used_draw_types = ["Pia", "White Matter", "Soma", "Cortical Layers"]
    df = df.loc[df["draw_type"].isin(used_draw_types), :].drop_duplicates(subset=["biospecimen_id", "poly_coords"])

    return df