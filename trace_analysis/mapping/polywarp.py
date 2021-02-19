import numpy as np
from skimage.transform._geometric import GeometricTransform

# If you can think of a better name ... PolynomialTransform2 possibly
class PolywarpTransform(GeometricTransform):
    # TODO: Update docstrings
    """2D polynomial transformation.

    Has the following form::

        X = sum[j=0:order]( sum[i=0:j]( a_ji * x**(j - i) * y**i ))
        Y = sum[j=0:order]( sum[i=0:j]( b_ji * x**(j - i) * y**i ))

    Parameters
    ----------
    params : (2, N) array, optional
        Polynomial coefficients where `N * 2 = (order + 1) * (order + 2)`. So,
        a_ji is defined in `params[0, :]` and b_ji in `params[1, :]`.

    Attributes
    ----------
    params : (2, N) array
        Polynomial coefficients where `N * 2 = (order + 1) * (order + 2)`. So,
        a_ji is defined in `params[0, :]` and b_ji in `params[1, :]`.

    """

    def __init__(self, params=None):
        if params is None:
            # default to transformation which preserves original coordinates
            params = (np.array([[0, 0], [1, 0]]), np.array([[0, 1], [0, 0]]))

        shapes = np.vstack([np.array(p.shape) for p in params])
        if not np.all((shapes-shapes[0,0])==0):
            raise ValueError("invalid shape of transformation parameters")
        self.params = params

        # TODO: Make order a class value, so it can be set on initialisation.

    def estimate(self, src, dst, order=3):
        """Estimate the transformation from a set of corresponding points.

        You can determine the over-, well- and under-determined parameters
        with the total least-squares method.

        Number of source and destination coordinates must match.

        The transformation is defined as::

            X = sum[j=0:order]( sum[i=0:j]( a_ji * x**(j - i) * y**i ))
            Y = sum[j=0:order]( sum[i=0:j]( b_ji * x**(j - i) * y**i ))

        These equations can be transformed to the following form::

            0 = sum[j=0:order]( sum[i=0:j]( a_ji * x**(j - i) * y**i )) - X
            0 = sum[j=0:order]( sum[i=0:j]( b_ji * x**(j - i) * y**i )) - Y

        which exist for each set of corresponding points, so we have a set of
        N * 2 equations. The coefficients appear linearly so we can write
        A x = 0, where::

            A   = [[1 x y x**2 x*y y**2 ... 0 ...             0 -X]
                   [0 ...                 0 1 x y x**2 x*y y**2 -Y]
                    ...
                    ...
                  ]
            x.T = [a00 a10 a11 a20 a21 a22 ... ann
                   b00 b10 b11 b20 b21 b22 ... bnn c3]

        In case of total least-squares the solution of this homogeneous system
        of equations is the right singular vector of A which corresponds to the
        smallest singular value normed by the coefficient c3.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.
        order : int, optional
            Polynomial order (number of coefficients is order + 1).

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """
        # xs = src[:, 0]
        # ys = src[:, 1]
        # xd = dst[:, 0]
        # yd = dst[:, 1]
        # rows = src.shape[0]
        #
        # # number of unknown polynomial coefficients
        # order = safe_as_int(order)
        # u = (order + 1) * (order + 2)
        #
        # A = np.zeros((rows * 2, u + 1))
        # pidx = 0
        # for j in range(order + 1):
        #     for i in range(j + 1):
        #         A[:rows, pidx] = xs ** (j - i) * ys ** i
        #         A[rows:, pidx + u // 2] = xs ** (j - i) * ys ** i
        #         pidx += 1
        #
        # A[:rows, -1] = xd
        # A[rows:, -1] = yd
        #
        # _, _, V = np.linalg.svd(A)
        #
        # # solution is right singular vector that corresponds to smallest
        # # singular value
        # params = - V[-1, :-1] / V[-1, -1]
        #
        # self.params = params.reshape((2, u // 2))
        #
        # return True
        self.params = polywarp(dst, src, degree=order)

        return True

    def __call__(self, coords):
        """Apply forward transformation.

        Parameters
        ----------
        coords : (N, 2) array
            source coordinates

        Returns
        -------
        coords : (N, 2) array
            Transformed coordinates.

        """
        # x = coords[:, 0]
        # y = coords[:, 1]
        # u = len(self.params.ravel())
        # # number of coefficients -> u = (order + 1) * (order + 2)
        # order = int((- 3 + math.sqrt(9 - 4 * (2 - u))) / 2)
        # dst = np.zeros(coords.shape)
        #
        # pidx = 0
        # for j in range(order + 1):
        #     for i in range(j + 1):
        #         dst[:, 0] += self.params[0, pidx] * x ** (j - i) * y ** i
        #         dst[:, 1] += self.params[1, pidx] * x ** (j - i) * y ** i
        #         pidx += 1

        dst = polywarp_apply(self.params[0], self.params[1], coords)

        return dst

    def inverse(self, coords):
        raise Exception(
            'There is no explicit way to do the inverse polynomial '
            'transformation. Instead, estimate the inverse transformation '
            'parameters by exchanging source and destination coordinates,'
            'then apply the forward transformation.')



# Original code source: Trey Wenger - August 2015
# Implementation of IDL's polar.pro
# Shamelessly copied, well tested against IDL procedure

def polywarp(xy_out, xy_in, degree=3):
    """
    TODO: Update docstring and function style
    originally polywarp(Xout,Yout,Xin,Yin,degree=3)
    Fit a function of the form
    Xout = sum over i and j from 0 to degree of: kx[i,j] * Xin^j * Yin^i
    Yout = sum over i and j from 0 to degree of: ky[i,j] * Xin^j * Yin^i
    Return kx, ky
    len(xo) must be greater than or equal to (degree+1)^2
    """
    x_out = xy_out[:, 0]
    y_out = xy_out[:, 1]
    x_in = xy_in[:, 0]
    y_in = xy_in[:, 1]

    if len(x_in) != len(y_in) or len(x_in) != len(x_out) or len(x_in) != len(y_out):
        print("Error: length of xo, yo, xi, and yi must be the same")
        return

    if len(x_in) < (degree + 1.) ** 2.:
        # print ("Error: length of arrays must be greater than (degree+1)^2")
        # return
        new_degree = int(np.floor(np.sqrt(len(x_in)) - 1))
        print(f'Too few datapoints for calculation with degree {degree}, reduced degree to {new_degree}')
        degree = new_degree

    # ensure numpy arrays
    x_in = np.array(x_in)
    y_in = np.array(y_in)
    x_out = np.array(x_out)
    y_out = np.array(y_out)

    # set up some useful variables
    degree2 = (degree + 1) ** 2
    x = np.array([x_out, y_out])
    u = np.array([x_in, y_in])
    ut = np.zeros([len(x_in), degree2])
    u2i = np.zeros(degree + 1)

    for i in range(len(x_in)):
        u2i[0] = 1.
        zz = u[1, i]
        for j in range(1, degree + 1):
            u2i[j] = u2i[j - 1] * zz

        # print ("u2i",u2i)

        ut[i, 0:degree + 1] = u2i

        for j in range(1, degree + 1):
            ut[i, j * (degree + 1):j * (degree + 1) + degree + 1] = u2i * u[0, i] ** j
    # print ("ut",ut)

    uu = ut.T
    #  print( "uu",uu)
    kk = np.dot(np.linalg.inv(np.dot(uu, ut)), uu).T
    # print( "kk",kk)
    # print( "x[0,:]",x[0,:])
    kx = np.dot(kk.T, x[0, :]).reshape(degree + 1, degree + 1)
    # print ("kx",kx)
    ky = np.dot(kk.T, x[1, :]).reshape(degree + 1, degree + 1)
    # print ("ky",ky)

    return kx, ky


def polywarp_apply(P, Q, pts1):
    # TODO: Add docstring and function style
    deg = len(P) - 1
    dst = np.ones(np.shape(pts1))
    dst[:, 0] = [
        np.sum([P[ii, jj] * pts1[kk, 0] ** ii * pts1[kk, 1] ** jj for ii in range(deg + 1) for jj in range(deg + 1)])
        for kk in range(len(pts1))]
    dst[:, 1] = [
        np.sum([Q[ii, jj] * pts1[kk, 0] ** ii * pts1[kk, 1] ** jj for ii in range(deg + 1) for jj in range(deg + 1)])
        for kk in range(len(pts1))]
    return (dst)


def translate(displacement, degree=3):
    kx = np.zeros((degree + 1, degree + 1))
    ky = np.zeros((degree + 1, degree + 1))
    kx[0, 0] = displacement[0]
    kx[1, 0] = 1
    ky[0, 0] = displacement[1]
    ky[0, 1] = 1
    return kx, ky


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    Npoints = 40

    np.random.seed(32)
    coordinates = np.random.rand(Npoints, 2) * 1000

    plt.figure()
    plt.scatter(coordinates[:, 0], coordinates[:, 1], c='green')
    plt.scatter(coordinates[:, 0] + 100, coordinates[:, 1] + 200, c='orange')
    Pt, Qt = translate([100, 200])
    new_coordinates = polywarp_apply(Pt, Qt, coordinates)
    plt.scatter(new_coordinates[:, 0], new_coordinates[:, 1], facecolors='none', edgecolors='r')

