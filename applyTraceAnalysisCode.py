# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 10:50:59 2018

@author: ivoseverins
"""

from traceAnalysisCode import *
import os
from pick_spots_akaze import mapping

#mainPath = r'D:\ivoseverins\SURFdrive\Promotie\Code\Python\traceAnalysis\twoColourExampleData\HJ A'
mainPath=r'H:\projects\research practicum\single molecule fluorescence\Matlab\HJA-data from Ivo '
os.chdir(os.path.dirname(os.path.abspath(__file__)))

#mainPath = './twoColourExampleData/HJ A'

T,fig=mapping()


exp = Experiment(mainPath)


exp.histogram()






#files = os.listdir(e.mainPath)
#filenames = list(set([re.search('hel[0-9]*',file).group() for file in files]))
#
#for filename in filenames: e.addFile(filename,r'.') 



