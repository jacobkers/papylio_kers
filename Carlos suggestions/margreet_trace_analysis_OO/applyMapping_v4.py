
from ImageCollection import ImageCollection

tetra_fn='/home/carlos/PycharmProjects/margreet_code_review/rough_ave.tif'
image_fn ='/home/carlos/PycharmProjects/margreet_code_review/hel4.pma'

imc = ImageCollection(tetra_fn, image_fn)
img = imc.get_image(1)