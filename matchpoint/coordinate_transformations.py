import numpy as np

def translate(displacement):
    T = np.array([[1, 0, displacement[0]], [0, 1, displacement[1]], [0, 0, 1]])
    return T


def rotate(angle, origin=np.array([0, 0, 1])):
    angle = np.array(angle)
    # angle = np.radians(angle)
    R = np.array([[np.cos(angle), np.sin(angle), 0], [-np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    return translate(origin) @ R @ translate(-origin)


def magnify(magnification, origin=np.array([0, 0, 1])):
    magnification = np.array(magnification)
    if magnification.size == 1:
        magnification = np.append(magnification, magnification)
    M = np.diag(np.append(magnification, 1))
    return translate(origin) @ M @ translate(-origin)


def reflect(axis=0):
    if axis == 0:
        R = np.diag([1, -1, 1])
    elif axis == 1:
        R = np.diag([-1, 1, 1])
    return R


def transform(pointSet, transformationMatrix=None, returnTransformationMatrix = False, **kwargs):
    if len(pointSet) == 0: return pointSet
    pointSet = np.append(pointSet, np.ones((pointSet.shape[0], 1)), axis=1)
    transformations = {
        'translation': translate,
        'rotation': rotate,
        'magnification': magnify,
        'reflection': reflect,

        't': translate,
        'r': rotate,
        'm': magnify
    }

    if transformationMatrix is None:
        transformationMatrix = np.identity(3)

    for key, value in kwargs.items():
        transformationMatrix = transformations.get(key)(value) @ transformationMatrix
        # print("%s == %s" %(key, value))

    if returnTransformationMatrix:
        return (transformationMatrix @ pointSet.T)[0:2, :].T, transformationMatrix
    else:
        return (transformationMatrix @ pointSet.T)[0:2, :].T


def scale(transformation_matrices):
    return np.sqrt(np.sum(transformation_matrices ** 2, axis=1))[:,:2]

def rotation(transformation_matrices):
    return np.atleast_2d(np.arctan2(transformation_matrices[:, 1, 0], transformation_matrices[:, 0, 0])).T

def shear(transformation_matrices):
    raise Warning('Not tested yet')
    beta = np.atleast_2d(np.arctan2(- transformation_matrices[:, 0, 1], transformation_matrices[:, 1, 1])).T
    return beta - rotation(transformation_matrices)

def translation(transformation_matrices):
    # raise Warning('Not tested yet')
    return transformation_matrices[:, 0:2, 2]


parameter_function_dict = {'translation': translation, 'rotation': rotation, 'scale': scale, 'shear': shear}

def parameters_from_transformation_matrices(transformation_matrices, parameters=['translation','rotation','scale','shear']):
    parameter_values = []
    for parameter in parameters:
        parameter_values.append(parameter_function_dict[parameter](transformation_matrices))
    return parameter_values

