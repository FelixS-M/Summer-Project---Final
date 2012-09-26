#Part of the Bound (Realistic stars) data pipeline

#Produces any relevant QoF vs StN plots, which will form a key part of the false postive detection plan

import numpy as np
import os
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import math

def listarraymaker():
	#First load the lists which tell us which category/array to place each planet in
	Main_FP_List = np.loadtxt('Main FP List.txt', dtype = {'names': ('KeplerID', 'planetnumber'), 'formats': (int, int)}, skiprows = 1)
	Slightly_FP_List = np.loadtxt('Slightly FP List.txt', dtype = {'names': ('KeplerID', 'planetnumber'), 'formats': (int, int)}, skiprows = 1)
	CT_FP_List = np.loadtxt('Cannot Tell FP List.txt', dtype = {'names': ('KeplerID', 'planetnumber'), 'formats': (int, int)}, skiprows = 1)
	Both_Flat_List = np.loadtxt('Both Flat.txt', dtype = {'names': ('KeplerID', 'planetnumber'), 'formats': (int, int)}, skiprows = 1)
	#Load the array which contains the QoF and StN data
	DATAarray = np.load('QoFdataarray.npy')
	#Data types for the arrays which are going to be generated
	dtypes = [('KeplerID', int), ('planetnumber', int), ('QualityofFit', float), ('sigmazero', float), ('Signal-to-Noise', float), ('BoundBetter', bool)]
	for star in DATAarray: #loop over all the stars
		star['QualityofFit'] = np.log10(star['QualityofFit'])
		star['Signal-to-Noise'] = np.log10(star['Signal-to-Noise'])
		values = [(star[0], star[1], star[2], star[3], star[4], star[5])] #values to be stored in the arrays
		lists = [Main_FP_List, Slightly_FP_List, CT_FP_List, Both_Flat_List] #A list of all the different lists
		names = ['Main_FP_', 'Slightly_FP_', 'CT_FP_', 'Both_Flat_'] #Names for the files which will contain the list arrays
		breakout = False #an indicator which is used to tell when the below loop should be broken
		for i in range(0, len(lists), 1): #loop over all the lists
			if breakout == True: break #When indicated, break
			for line in lists[i]: #now loop over the entries in the list
				if star['KeplerID'] == line['KeplerID'] and star['planetnumber'] == line['planetnumber']: #When the keplerID and planet number match, add the result to the current list_array
					savedataarray = np.array(values, dtype=dtypes)
					if os.path.exists('{}dataarray.npy'.format(names[i])) == False:
						outdataarray = np.copy(savedataarray)
					else:
						outdataarray = np.load('{}dataarray.npy'.format(names[i]))
						outdataarray = np.concatenate((outdataarray, savedataarray), axis=0)
					sortedoutdataarray = np.sort(outdataarray, order='KeplerID')
					np.save('{}dataarray'.format(names[i]), sortedoutdataarray) #save the sorted structured array
					breakout = True #Now that the data has been added to an list_array, the loop outside this loop can be broken!
					break
		if breakout == False: #If the star has not been added to any of the array lists, it should be added to the remaining candidates list
			savedataarray1 = np.array(values, dtype=dtypes)
			if os.path.exists('remainingdataarray.npy') == False:
				outdataarray1 = np.copy(savedataarray1)
			else:
				outdataarray1 = np.load('remainingdataarray.npy')
				outdataarray1 = np.concatenate((outdataarray1, savedataarray1), axis=0)
			sortedoutdataarray1 = np.sort(outdataarray1, order='KeplerID')
			np.save('remainingdataarray', sortedoutdataarray1) #save the sorted structured array
	#Now load the complete lists, which can then be returned to the main() function
	Main_FP_Array = np.load('Main_FP_dataarray.npy')
	Slightly_FP_Array = np.load('Slightly_FP_dataarray.npy')
	CT_FP_Array = np.load('CT_FP_dataarray.npy')
	Both_Flat_Array = np.load('Both_Flat_dataarray.npy')
	Remain_Array = np.load('remainingdataarray.npy')
	print len(Remain_Array['KeplerID']), (len(DATAarray['KeplerID']) - len(Main_FP_Array['KeplerID']) - len(Slightly_FP_Array['KeplerID']) - len(CT_FP_Array['KeplerID']) - len(Both_Flat_Array['KeplerID'])) #testing, confirming that there are no double entries
	return Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array


def largebasicplotfull(Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array):
	#Used to define the currently arbitary/placed by eye - should be easy enough to adjust the code in the future when proper dividing lines are chosen - data just need to be in an array
	xz = np.arange(-0.5, 3.1, 0.1)
	yz = np.arange(-1.2, 2.4, 0.1)
	xz2 = np.arange(-0.1, 2.1, 0.1)
	yz2 = np.arange(-1.2, 3.2, 0.2)
	xp = np.linspace(-4,3,10)
	yp = (np.copy(xp) * 0) - 1.2
	plt.clf() #clear all figures
	plt.figure(figsize=(23.5, 13.0)) #create a figure for this plot
	fig, ax = plt.subplots(1) #create a subplot within the figure
	f1 = ax.plot( Remain_Array['Signal-to-Noise'], Remain_Array['QualityofFit'], ',', label = ['Planetary Candidates'], color = 'green') #plot both sets of data
	f3 = ax.plot( Slightly_FP_Array['Signal-to-Noise'], Slightly_FP_Array['QualityofFit'], ',', label = ['Possible FP Selection'], color = 'pink') #plot both sets of data
	f4 = ax.plot( CT_FP_Array['Signal-to-Noise'], CT_FP_Array['QualityofFit'], ',', label = ['Hard To/Cannot Tell'], color = 'purple') #plot both sets of data
	f5 = ax.plot( Both_Flat_Array['Signal-to-Noise'], Both_Flat_Array['QualityofFit'], ',', label = ['Both LC\'s Flat'], color = 'orange') #plot both sets of data
	f2 = ax.plot( Main_FP_Array['Signal-to-Noise'], Main_FP_Array['QualityofFit'], ',', label = ['Main FP Selection'], color = 'red') #plot both sets of data
	ax.plot(xp,yp,'-', color = 'grey')
	ax.plot(xz,yz,'-', color = 'grey')
	ax.plot(xz2,yz2,'-', color = 'grey')
	leg = ax.legend([f1,f2,f3,f4,f5], ['Planetary Candidates', 'Main FP Selection', 'Possible FP Selection', 'Hard To/Cannot Tell', 'Both LC\'S Flat'], fancybox=True, loc = 'upper left') #plot the legend, in the lower righthand corner
	leg.get_frame().set_alpha(0.5) #make the legend slightly transparent
	plt.xlabel('log10(Signal-to-Noise Ratio)')
	plt.ylabel('log10(Quality of Fit)')
	plt.title('Signal-to-Noise Ratio vs Quality of Fit')
	plt.savefig('Large StN vs QoF [full].pdf') #save the plot
	return

def largebasicplotimpactpar(Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array):
	impactparlist = np.load('impactpardiscardarray.npy')
	fullarray = np.concatenate((Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array), axis=0)
	for star in fullarray:
		for entry in impactparlist:
			if star['KeplerID'] == entry['KeplerID'] and star['planetnumber'] == entry['planetnumber']:
				if not 'ximpar' in locals():
					ximpar = np.array([star['Signal-to-Noise']])
					yimpar = np.array([star['QualityofFit']])
				else:
					ximpar = np.append(ximpar, [star['Signal-to-Noise']])
					yimpar = np.append(yimpar, [star['QualityofFit']])
	#Used to define the currently arbitary/placed by eye - should be easy enough to adjust the code in the future when proper dividing lines are chosen - data just need to be in an array
	xz = np.arange(-0.5, 3.1, 0.1)
	yz = np.arange(-1.2, 2.4, 0.1)
	xz2 = np.arange(-0.1, 2.1, 0.1)
	yz2 = np.arange(-1.2, 3.2, 0.2)
	xp = np.linspace(-4,3,10)
	yp = (np.copy(xp) * 0) - 1.2
	plt.clf() #clear all figures
	plt.figure(figsize=(23.5, 13.0)) #create a figure for this plot
	fig, ax = plt.subplots(1) #create a subplot within the figure
	f1 = ax.plot( Remain_Array['Signal-to-Noise'], Remain_Array['QualityofFit'], ',', label = ['Planetary Candidates'], color = 'green') #plot both sets of data
	f3 = ax.plot( Slightly_FP_Array['Signal-to-Noise'], Slightly_FP_Array['QualityofFit'], ',', label = ['Possible FP Selection'], color = 'pink') #plot both sets of data
	f4 = ax.plot( CT_FP_Array['Signal-to-Noise'], CT_FP_Array['QualityofFit'], ',', label = ['Hard To/Cannot Tell'], color = 'purple') #plot both sets of data
	f5 = ax.plot( Both_Flat_Array['Signal-to-Noise'], Both_Flat_Array['QualityofFit'], ',', label = ['Both LC\'s Flat'], color = 'orange') #plot both sets of data
	f2 = ax.plot( Main_FP_Array['Signal-to-Noise'], Main_FP_Array['QualityofFit'], ',', label = ['Main FP Selection'], color = 'red') #plot both sets of data
	f6 = ax.plot( ximpar, yimpar, ',', label = ['Impact Parameter Exclude'], color='cyan')
	ax.plot(xp,yp,'-', color = 'grey')
	ax.plot(xz,yz,'-', color = 'grey')
	ax.plot(xz2,yz2,'-', color = 'grey')
	leg = ax.legend([f1,f2,f3,f4,f5,f6], ['Planetary Candidates', 'Main FP Selection', 'Possible FP Selection', 'Hard To/Cannot Tell', 'Both LC\'S Flat', 'Impact Parameter Exclude'], fancybox=True, loc = 'upper left') #plot the legend, in the lower righthand corner
	leg.get_frame().set_alpha(0.5) #make the legend slightly transparent
	plt.xlabel('log10(Signal-to-Noise Ratio)')
	plt.ylabel('log10(Quality of Fit)')
	plt.title('Signal-to-Noise Ratio vs Quality of Fit')
	plt.savefig('Large StN vs QoF [impactpar].pdf') #save the plot
	return

def multiplotfull(Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array):
	fullarray = np.concatenate((Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array), axis=0)
	for i in range(0, len(fullarray['KeplerID']), 1):
		if np.isfinite(fullarray['QualityofFit'][i]) == False:
			fullarray['QualityofFit'][i] = np.nan
		elif np.isfinite(fullarray['Signal-to-Noise'][i]) == False:
			fullarray['Signal-to-Noise'][i] = np.nan
	arrays = [Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array]
	colours = ['red', 'pink', 'purple', 'orange', 'yellow', 'green']
	names = ['Main FP Selection', 'Possible FP Selection', 'Hard To_Cannot Tell', 'Both LC\'S Flat', 'Planetary Candidates']
	for i in range(0, len(arrays), 1):
		array = arrays[i]
		plt.clf() #clear all figures
		plt.figure(figsize=(23.5, 13.0)) #create a figure for this plot
		fig, ax = plt.subplots(1) #create a subplot within the figure
		ax.plot( array['Signal-to-Noise'], array['QualityofFit'], '.', label = [names[i]], color = colours[i]) #plot both sets of data
		plt.xlabel('log10(Signal-to-Noise Ratio)')
		plt.ylabel('log10(Quality of Fit)')
		plt.title('Signal-to-Noise Ratio vs Quality of Fit ({})'.format(names[i]))
		plt.axis([-4, 3, -3, 4])
		#plt.axis([np.nanmin(fullarray['Signal-to-Noise']), np.nanmax(fullarray['Signal-to-Noise']), np.nanmin(fullarray['QualityofFit']), np.nanmax(fullarray['QualityofFit'])])
		plt.savefig('Large StN vs QoF ({}) [full].pdf'.format(names[i])) #save the plot
	return

def largebasicplot_fitted(Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array):
	z = np.polyfit(Main_FP_Array['Signal-to-Noise'], Main_FP_Array['QualityofFit'],1)
	print np.average(Main_FP_Array['QualityofFit'])
	print z
	p = np.poly1d(z)
	xp = np.linspace(0,10,100)
	yp = p(xp)
	print xp
	print yp
	plt.clf() #clear all figures
	plt.figure(figsize=(23.5, 13.0)) #create a figure for this plot
	fig, ax = plt.subplots(1) #create a subplot within the figure
	f1 = ax.plot( Remain_Array['Signal-to-Noise'], Remain_Array['QualityofFit'], ',', label = ['Planetary Candidates'], color = 'green') #plot both sets of data
	f3 = ax.plot( Slightly_FP_Array['Signal-to-Noise'], Slightly_FP_Array['QualityofFit'], ',', label = ['Possible FP Selection'], color = 'pink') #plot both sets of data
	f4 = ax.plot( CT_FP_Array['Signal-to-Noise'], CT_FP_Array['QualityofFit'], ',', label = ['Hard To/Cannot Tell'], color = 'purple') #plot both sets of data
	f5 = ax.plot( Both_Flat_Array['Signal-to-Noise'], Both_Flat_Array['QualityofFit'], ',', label = ['Both LC\'s Flat'], color = 'orange') #plot both sets of data
	f2 = ax.plot( Main_FP_Array['Signal-to-Noise'], Main_FP_Array['QualityofFit'], ',', label = ['Main FP Selection'], color = 'red') #plot both sets of data
	ax.plot(xp,yp,'-', color='grey')
	leg = ax.legend([f1,f2,f3,f4,f5], ['Planetary Candidates', 'Main FP Selection', 'Possible FP Selection', 'Hard To/Cannot Tell', 'Both LC\'S Flat'], fancybox=True, loc = 'upper right') #plot the legend, in the lower righthand corner
	leg.get_frame().set_alpha(0.5) #make the legend slightly transparent
	plt.xlabel('log10(Signal-to-Noise Ratio)')
	plt.ylabel('log10(Quality of Fit)')
	plt.title('Signal-to-Noise Ratio vs Quality of Fit')
	plt.savefig('Large StN vs QoF [fitted].pdf') #save the plot
	return

def main():
	foldername = raw_input('Please enter the name of the folder in which contains the data: ') 
	os.chdir('{}'.format(foldername))
	Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array = listarraymaker()
	largebasicplotfull(Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array)
	largebasicplotimpactpar(Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array)
	#largebasicplot_fitted(Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array)
	multiplotfull(Main_FP_Array, Slightly_FP_Array, CT_FP_Array, Both_Flat_Array, Remain_Array)
	return

main()