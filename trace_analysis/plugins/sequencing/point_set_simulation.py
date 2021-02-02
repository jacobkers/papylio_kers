import numpy as np

def simulate_point_set(number_of_points, bounds = [[0,1],[0,1]]):
    bounds = np.array(bounds)
    dimensions = len(bounds)
    unit_coordinates = np.random.rand(number_of_points, dimensions)
    coordinates = unit_coordinates * (bounds[:,1]-bounds[:,0]) + bounds[:,0]

    return coordinates
