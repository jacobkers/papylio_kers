import numpy as np #scientific computing with Python
import matplotlib.pyplot as plt #Provides a MATLAB-like plotting framework

class Molecule:

    def __init__(self, file):
        self.file = file
        self.index = None
        self.coordinates = None
        self.intensity = None

        self.isSelected = False

        #self.steps = None  #Defined in other classes as: pd.DataFrame(columns=['frame', 'trace', 'state', 'method','thres'])
        #self.kon_boolean = None  # 3x3 matrix that is indicates whether the kon will be calculated from the beginning, in-between molecules or for the end only

    def I(self, emission, Ioff=0):
        return self.intensity[emission, :] - Ioff - self.file.background[emission]

    def E(self, Imin=0, alpha=0):  # alpha correction is not implemented yet, this is just a reminder
        red = np.copy(self.I(1))
        green = self.I(0)
        np.putmask(red, red<Imin, 0)  # the mask makes all elements of acceptor that are below the Imin zero, for E caclulation
        return (red - alpha*green) / (green + red - alpha*green)

    def plot(self, figure = None):
        if not figure: figure = plt.gcf()
        figure.clf()
        if len(self.file.experiment.pairs) > 0:
            axis_I = figure.add_subplot(211)
        else:
            axis_I = figure.gca()

        axis_I.set_xlabel('Time (s)')
        axis_I.set_ylabel('Intensity (a.u.)')
        axis_I.set_ylim(0,1000)
        for i, colour in enumerate(self.file.experiment.colours):
            axis_I.plot(self.I(i), colour)

        if len(self.file.experiment.pairs) > 0:
            axis_E = figure.add_subplot(212, sharex = axis_I)
            axis_E.set_xlabel('Time (s)')
            axis_E.set_ylabel('FRET (-)')
            axis_E.set_ylim(0,1)
            for i, pair in enumerate(self.file.experiment.pairs):
                axis_E.plot(self.E())

        #plt.show()
#MD190104: why not add a subplot with FRET here as well, to match with display Matlab?

#    @property
#    def find_steps(self):
#        return stepfinder