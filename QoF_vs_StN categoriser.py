#Part of the Bound (Realistic stars) data pipeline

#Produces any relevant QoF vs StN plots, which will form a key part of the false postive detection plan

import numpy as np
import os
import matplotlib.pyplot as plt
import math





def main():
	foldername = raw_input('Please enter the name of the folder in which contains the data: ') 
	os.chdir('{}'.format(foldername))
	DATAarray = np.load('QoFdataarray.npy')
	for star in DATAarray:
		log_StN = np.log10(star['Signal-to-Noise'])
		log_QoF = np.log10(star['QualityofFit'])
		False_Positive = False
		Planet_Candidate = False
		unsure_2 = False
		unsure_1 = False
		if log_QoF < -1.2:
			False_Positive = True
		elif log_QoF < (log_StN - 0.7):
			unsure_1 = True
		if log_QoF < ((2*log_StN)-1) and False_Positive == False:
			unsure_2 = True
		if False_Positive == False and unsure_2 == False and unsure_1 == False:
			Planet_Candidate = True
		dtypes = [('KeplerID', int), ('planetnumber', int), ('QualityofFit', float), ('Signal-to-Noise', float),\
			('False Positive', bool), ('Unsure y=2x-1', bool), ('Unsure y=x-0.7', bool), ('PlanetCandidate', bool)]
		values = [(star['KeplerID'], star['planetnumber'], star['QualityofFit'], star['Signal-to-Noise'], False_Positive, unsure_2, unsure_1, Planet_Candidate)]
		if not 'outarray' in locals():
			outarray = np.array(values, dtype=dtypes)
		else:
			array = np.array(values, dtype=dtypes)
			outarray = np.concatenate((outarray, array), axis=0)
	print outarray
	np.save('QoF and StN Classification', outarray)
	return

main()