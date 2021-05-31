import warnings
import numpy as np

# from skimage.transform._geometric import PolynomialTransform
#
# class PolynomialTransform(PolynomialTransform):
#     def estimate(self, src, dst, order=4):
#         src = src.astype(float)
#         dst = dst.astype(float)
#         return super().estimate(src, dst, order)
#
#     def __call__(self, coords):
#         coords = coords.astype(float)
#         return super().__call__(coords)

class PolynomialTransform:
    def estimate(self, source, destination, order=4):
        source = source.astype(float)
        A = np.vstack([source[:,0] ** i * source[:, 1] ** j for i in range(order + 1) for j in range(order + 1)]).T
        B = destination.astype(float)

        # Based on numpy.polyfit
        scale = np.sqrt(A * A).sum(axis=0)
        coeff, r, rank, s = np.linalg.lstsq(A / scale, B, rcond=None)
        coeff = (coeff.T / scale).T

        # warn on rank reduction, which indicates an ill conditioned matrix
        if rank != len(coeff):
            msg = "Estimation may be poorly conditioned"
            warnings.warn(msg, np.RankWarning, stacklevel=4)

        self.params = coeff
        print(rank)
        return True

    def __call__(self, coords):
        order = int(np.sqrt(len(self.params)))-1
        coords = coords.astype(float)
        A = np.vstack([coords[:, 0] ** i * coords[:, 1] ** j for i in range(order + 1) for j in range(order + 1)]).T
        return A@self.params
