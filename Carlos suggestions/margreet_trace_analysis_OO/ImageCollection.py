
from load_file import read_one_page_pma, read_one_page_tif
from image_adapt.rolling_ball import rollingball

from image_adapt.find_threshold import remove_background
from image_adapt.find_threshold import get_threshold
import numpy as np
from cached_property import cached_property
from Mapping import Mapping
from Image import Image
import analyze_label
import cv2

class ImageCollection(object):
    def __init__(self, tetra_fn, image_fn):
        self.image_fn = image_fn
        self.mapping = Mapping(tetra_fn)
        self.set_background_and_transformation()
        (self.background,
         self.threshold,
         self.pts_number,
         self.dstG,
         self.ptsG) = self.set_background_and_transformation()

    @cached_property
    def read_one_page(self):
        if '.pma' in self.image_fn:
            return read_one_page_pma
        elif '.tif' in self.image_fn:
            return read_one_page_tif

    def set_background_and_transformation(self):
        """
        sum 20 image, find spots& background. Then loop over all image, do background subtract+ extract traces
        :return:
        """
        _, hdim, vdim, n_images = self.read_one_page(self.image_fn, pageNb=0)
        im_array = np.dstack([(self.read_one_page(self.image_fn, pageNb=ii)[0]).astype(float) for ii in range(20)])
        im_mean20 = np.mean(im_array, axis=2).astype(int)
        bg = rollingball(im_mean20)[1]
        im_mean20_correct = im_mean20 - bg
        im_mean20_correct[im_mean20_correct < 0] = 0
        threshold = get_threshold(im_mean20_correct)

        pts_number, label_size, ptsG = analyze_label.analyze(im_mean20_correct[:, 0:int(vdim / 2)])
        dstG = cv2.perspectiveTransform(ptsG.reshape(-1, 1, 2),
                                        np.linalg.inv(self.mapping.transform_matrix))
        dstG = dstG.reshape(-1, 2)
        dstG = np.array([[ii[0] + 256, ii[1]] for ii in dstG])
        return bg, threshold, pts_number, dstG, ptsG

    def subtract_background(self, im):
        im_correct = im - self.background
        im_correct[im_correct < 0] = 0
        return remove_background(im_correct, self.threshold)

    def get_image(self, idx):
        img, hdim, vdim, n_images = self.read_one_page(self.image_fn, pageNb=idx)
        img = self.subtract_background(img)
        return Image(img, vdim, self.mapping.transform_matrix, self.ptsG, self.dstG, self.pts_number)
