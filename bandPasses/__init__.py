__DESCRIPTION__="Handles bandpasses in a number of formats"
from bandPassToRimo import *

__MY_ABSOLUTE_PATH__ = '/'.join(__file__.split('/')[:-1])
DiodeWeights=diodeWeights(__MY_ABSOLUTE_PATH__+'/LFI_bandpasses_QUCS_20091109/weights_march15_2010.txt')

import rimo_fits
BPMeasuredByDiode=rimo_fits.LFI_BandPass_Fits(__MY_ABSOLUTE_PATH__ +'/LFI_bandpasses_measured_by_diode.fits')
