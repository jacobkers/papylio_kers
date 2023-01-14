import numpy as np
from shapely.geometry import Polygon, MultiPoint, LineString, Point

def overlap_vertices(vertices_A, vertices_B):
    polygon_A = Polygon(vertices_A)
    polygon_B = Polygon(vertices_B)
    if polygon_A.intersects(polygon_B):
        polygon_intersection = polygon_A.intersection(polygon_B)
    # return np.array(polygon_overlap.exterior.coords.xy).T[:-1]

        return np.array(polygon_intersection.boundary.coords)[:-1]
    else:
        return np.empty((0,2))

def area(vertices):
    if len(vertices) < 3:
        return 0
    else:
        return Polygon(vertices).area

def crop_coordinates_indices(coordinates, vertices):
    # return pth.Path(vertices).contains_points(coordinates)
    cropped_coordinates = crop_coordinates(coordinates, vertices)
    return np.array([c in cropped_coordinates for c in coordinates])

def crop_coordinates(coordinates, vertices):
    if len(vertices) > 0 and len(np.atleast_2d(coordinates)) > 0:
        # return np.atleast_2d(Polygon(vertices).intersection(MultiPoint(coordinates)))
        pointset_intersected = Polygon(vertices).intersection(MultiPoint(coordinates))
        if isinstance(pointset_intersected, MultiPoint):
            return np.atleast_2d(LineString(pointset_intersected.geoms).coords)
        elif isinstance(pointset_intersected, Point) and len(pointset_intersected.coords) > 0:
            return np.array(pointset_intersected.coords)
        else:
            return np.empty((0, 2))
        # Hopefully shapely 2.0 will be more consistent
    else:
        return np.empty((0, 2)) # np.atleast_2d([])

    # bounds.sort(axis=0)
    # selection = (coordinates[:, 0] > bounds[0, 0]) & (coordinates[:, 0] < bounds[1, 0]) & \
    #             (coordinates[:, 1] > bounds[0, 1]) & (coordinates[:, 1] < bounds[1, 1])
    # return coordinates[selection]


def determine_vertices(point_set, margin=0):
    return np.array(MultiPoint(point_set).convex_hull.buffer(margin, join_style=1).boundary.coords)[:-1]

def vertices_with_margin(vertices, margin):
    return np.array(Polygon(vertices).buffer(margin).exterior.coords)