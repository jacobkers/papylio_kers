from skimage.transform._geometric import PolynomialTransform

class PolynomialTransform(PolynomialTransform):
    def estimate(self, src, dst, order=4):
        src = src.astype(float)
        dst = dst.astype(float)
        return super().estimate(src, dst, order)

    def __call__(self, coords):
        coords = coords.astype(float)
        return super().__call__(coords)