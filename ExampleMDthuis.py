import os
import numpy as np
from matplotlib import pyplot as plt

from trace_analysis import Experiment
from trace_analysis import InteractivePlot

# Define path to data, replace by your own directory
mainPath = r'D:\pathToDataFolder'
mainPath = r'D:\ivoseverins\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\20141017 - Holliday junction - Copy'
mainPath = r'C:\MargreetDATA\TUDelft\python\data' 
# Initialize an experiment
exp = Experiment(mainPath)
plt.close('all')
 
# Print files in experiment
print(exp.files)

## Perform mapping
mapping_file_index = -1
mapping_file = exp.files[mapping_file_index]

#instead of mapping_file.perform_mapping():
configuration = mapping_file.experiment.configuration['mapping']
configuration['peak_finding']={'method': 'AKAZE', 'threshold': 'None' , 'minimum_area': 5,'maximum_area': 50}

if 0: #linear mapping
    A=mapping_file
    A.perform_mapping(configuration,  transformation_type = 'linear')
    
    # different ways of getting the same plot
    At=A.mapping.transform_coordinates(A.mapping.source) # somehow uses non-linear, if you did B.perform mapping nonlinear before 
    if 0:
        plt.figure(), 
        plt.scatter(A.mapping.source[:,0],A.mapping.source[:,1])
        plt.scatter(A.mapping.destination[:,0],A.mapping.destination[:,1])
    
        plt.scatter(At[:,0],At[:,1])
    elif 1:
        plt.figure()
        from trace_analysis.mapping.icp import scatter_coordinates, show_point_connections, nearest_neighbor_pair
        scatter_coordinates([A.mapping.source,A.mapping.destination,At])
        Adist, Asi, Adi =           nearest_neighbor_pair(A.mapping.destination[:, :2], At[:, :2])
        show_point_connections(A.mapping.destination[Asi],At[Adi])
    elif 0:
        A.mapping.show_mapping_transformation() #does not show distances!!

if 1:
    B=mapping_file
    B.perform_mapping(configuration,  transformation_type = 'nonlinear')
    Bt=B.mapping.transform_coordinates(B.mapping.source) # somehow uses non-linear 
    plt.figure()
    from trace_analysis.mapping.icp import scatter_coordinates, show_point_connections, nearest_neighbor_pair
    scatter_coordinates([B.mapping.source,B.mapping.destination,Bt])
    Bdist, Bsi, Bdi = \
                    nearest_neighbor_pair(B.mapping.destination[:, :2], Bt[:, :2])
    show_point_connections(B.mapping.destination[Bsi],Bt[Bdi])