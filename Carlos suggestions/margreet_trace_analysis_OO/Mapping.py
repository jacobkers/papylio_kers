import autopick.pick_spots_akaze_manual
import autopick.pick_spots_akaze_final

class Mapping(object):
    def __init__(self, tetra_fn):
        self.tetra_fn = tetra_fn
        self.gui_start()

    @property
    def tetra_fn(self):
        return self._tetra_fn

    @tetra_fn.setter
    def tetra_fn(self, fn):
        self.transform_matrix = autopick.pick_spots_akaze_final.mapping(fn,
                                                                        show=1,
                                                                        bg=None,
                                                                        tol=0, f=10000)[0]
        self._tetra_fn = fn

    def gui_start(self):
        x = input('Are you satisfied with the mapping (yes/no)?')
        if x[0] != 'y':  # do manual mapping
            self.transform_matrix = autopick.pick_spots_akaze_manual.mapping(self.tetra_fn,
                                                                             show=1,
                                                                             bg=None,
                                                                             tol=0, f=10000)[0]