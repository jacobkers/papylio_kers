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


def transform(pointSet, transformationMatrix=None, **kwargs):
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

    return (transformationMatrix @ pointSet.T)[0:2, :].T