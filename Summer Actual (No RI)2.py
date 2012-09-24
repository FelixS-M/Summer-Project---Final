import pyfits
import math
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import os
#np.set_printoptions(threshold='nan')

def batchmaster():
	filelistlist = dataindexer('output.txt') #create the INDEX
	numberstars = len(filelistlist) #calculate the number of stars
	print numberstars #print the number of stars availible
	lowerlimit = int(raw_input('Please enter the lower limit (previous upperlimit): ')) #ask the user for the lower limit
	upperlimit = int(raw_input('Please enter the upperlimit limit: ')) #ask the user for the upper limit
	Qmax =  int(raw_input('Please enter Qmax (maximum quarter of data to use): ')) #ask the user for the highest quarter of data to use
	foldername = raw_input('Please enter the name of the folder in which the data should be saved: ') #ask the user for the folder to store the data in
	if not os.path.isdir("./" + foldername + "/"): #check to see if the given folder exists
		os.mkdir("./" + foldername + "/") #if not, it creates the folder
	batchloop(filelistlist, lowerlimit, upperlimit, Qmax, foldername)
	os.system('cp KeplerLDC.txt Test\ 5\ July')
	return

def dataindexer(inputname): #data must be in the correct order - keplerID then DATE - i.e. files are in chronological order
	f = open(inputname, 'r')
	i=0
	for line in f: #read linee from the file
		#remove the new line characters
		filename = line.rstrip()
		#now extract the keplerID
		keplername = filename.split('-')
		#now need to remove the short cadence files - so split at the underscore
		fileend = filename.split('_')
		if fileend[1] == 'slc.fits': #if the file is short cadence, ignore it
			pass
		else: 
			if i==0: #for the first entry that is not short cadence
				filestore = [filename]  #create an list to store this stars set of files
				counter = 0
			else: #for all other entries
				test = filestore[0].split('-') #extract the keplerID of the file kepler********
				if keplername[0] == test[0]: #if it is the same as the first entry in the list - add it to the list 
					filestore.append(filename)
				elif keplername[0] != test[0]: #otherwise when all files with the same keplerid are in a list
					#now need to store the previous stars data list
					if counter == 0: #if this is the first star/list 
						filelist = [filestore] #create a list to store the lists of filenames for the star
					elif counter > 0: #for the next stars, append the list to the list of lists 
						filelist.append(filestore)
					#ahhhhhhhh - yes i did say list to the list of lists
					del filestore #reset the temporary list - ready for the next star
					#now we create a new list for the new stars filenames
					filestore = [filename] #this will always red box but it is not an error
					counter = 1 #just makes things neater for the loop
			i = 1 #just makes things easier for the loop
	return filelist

def batchloop(filelistlist, lowerlimit, upperlimit, Qmax, foldername):
	dtypestar = [('KeplerID', int), ('Teff', int), ('Logg', float), ('Metalicity', float), ('StellarRadius', float), ('gamma1', float), ('gamma2', float)] #data types for the stellar data array
	for i in range(lowerlimit,upperlimit, 1): # loop over the chosen stars
		if os.path.exists('{0}/Stellardataarray.npy'.format(foldername)) == True: #If the saved planetary data array exists, load it
			stellardata = np.load('{0}/Stellardataarray.npy'.format(foldername))
		#print filelistlist[i]
		Tempoutput, STAR = Qmulticalc(filelistlist[i], Qmax, foldername) #create/manipulate the lightcurves
		print 'The number of stars completed is {0} of {1}'.format((i+1 - lowerlimit), (upperlimit - lowerlimit)) #keep track of progress
		#now need to save the stellar data
		values = [(STAR[0], STAR[1], STAR[2], STAR[4], STAR[5], 0.0, 0.0)] #values to be stored in the stellar data array
		if not 'stellardata' in locals(): #if the stellar data array does not exist (i.e. in the case of lowerlimit = 0), create it
			stellardata = np.array(values, dtype = dtypestar)
		else: #otherwise append to the stellar data list
			temparray = np.array(values, dtype = dtypestar)
			stellardata = np.concatenate((stellardata, temparray), axis=0)
		np.save('{0}/Stellardataarray'.format(foldername), stellardata) #save the stellar data array as a numpy binary file
	#print stellardata #testing
	#Store the array in a format that is easily readable by python
	print('\a') #Make a noise when finished + notification 
	return

def Qmulticalc(files, Qmax, foldername):
	#print files
	#load the first quarter of data
	hdulistQ1 = pyfits.open(files[0])
	Qfinale = reducetotable(hdulistQ1)
	#Creates a new array to hold the uncorrected data
	Qraw =  np.copy(Qfinale)
	#Get key stellar properties
	STAR = keydata(hdulistQ1)
	hdulistQ1.close() #close the file when we are done with it
	for i in range(0,len(files),1):
		#This loops such that each new quarter is data is appended to the original data using
		#the fit (which is also calculated)
		if i==0:
			continue #since the data has already been loaded
		#For more than one quarter of data
		elif i > 0:
			Qraw = anotherquarter(STAR, Qraw, files[i], Qmax, foldername )
			#simply join the data together
	transitmaster(STAR, Qraw, foldername)
	plt.clf()
	#basicplot(Qraw, 10, '{0}/kepler-{1}-spec.pdf'.format(foldername, STAR[0]))
	return Qraw, STAR

def reducetotable( hdulist ):
	#reduces the data to a table
	tbdata = hdulist[1].data
	time = tbdata.field(0) #column 0
	flux = tbdata.field(7) #column 7
	#makes the fluxes more manageable in scale
	fluxr = flux/1000
	tab = np.column_stack((time,fluxr))
	# Next line removes any lines in the table than contain nan
	tabcom = tab[~np.isnan(tab).any(1)]
	#changes the table from columns to rows
	finaltable = tabcom.transpose()
	return finaltable

def keydata( hdulist ):
	#Reads Key stellar data
	KEPLERID = hdulist[0].header['KEPLERID']
	TEFF = hdulist[0].header['TEFF']
	LOGG = hdulist[0].header['LOGG']
	KP = hdulist[0].header['KEPMAG']
	M_H = hdulist[0].header['FEH']
	STELLARRADIUS = hdulist[0].header['RADIUS']
	print KEPLERID
	try: #Set Teff to 0 if it is not found
		int(TEFF)
	except AttributeError:
		TEFF = 0
	try: #Set Logg to 0 if it is not found
		int(LOGG)
	except AttributeError:
		LOGG = 0
	try: #set the kepler magnitude to 0 if it is not found
		int(KP)
	except AttributeError:
		KP = 0
	try: #set the metalicity to 0 if it is not found
		int(M_H)
	except AttributeError:
		M_H = 0
	Metalicity = 10**M_H
	try: #set the Stellar Radius to 0 if it is not found
		int(STELLARRADIUS)
	except AttributeError:
		STELLARRADIUS = 0
	# Now define a tuple of Key stellar data to return
	STAR = KEPLERID, TEFF, LOGG, KP, Metalicity, STELLARRADIUS
	#print STAR #testing
	return STAR

def anotherquarter(STAR, Qstart, file, Qmax, foldername):
	#Load in the file to be joined to the given data
	hdulistQadd = pyfits.open(file)
	Quarter = hdulistQadd[0].header['QUARTER'] #Used for the limited quarters of data requirement
	if Quarter > Qmax:
		Qjoined = np.copy(Qstart)
	else:
		#reduce the data, that is to be added, into a table
		Qadd = reducetotable(hdulistQadd)
		#then attach the data to the lightcurve
		Qjoined = combinequartersaverage(Qstart, Qadd)
	hdulistQadd.close()
	return Qjoined

def	combinequartersaverage(Data1, Data2): #Not sure which to use - depends on new data
	preaverage = averagefitback(Data1)
	postaverage = averagefitfor(Data2)
	#calculates the correction
	correction = preaverage - postaverage
	#creates a new array to hold the corrected data and applies the correction
	Data2corrected =  np.copy(Data2)
	Data2corrected[1] = Data2[1] + correction
	#combines the orignal data with the new quarter
	combinedtable = np.hstack((Data1, Data2corrected))
	return combinedtable
	
def averagefitback( DATA ): #Part of combinequartersaverage
	#Take the average of the last 80 data points of the given data set
	averageflux = np.average(DATA[1][-1:-120:-1])
	return averageflux

def averagefitfor( DATA ): #Part of combinequartersaverage
	#Take the average of the first 80 data points of the given data set
	averageflux = np.average(DATA[1][0:120:1]) 
	return averageflux

def basicplot( Dataa , num , filename):
	# Plots the data ---> slow 
	#only plots one in every num data points 
	order = Dataa.size / 2
	datarange = range(1,order,num)
	Data = np.take(Dataa, datarange, 1)
	plt.clf() #make sure we have a clean slate
	plt.plot( [Data[0]], [Data[1]], 'b+')
	plt.xlabel( 'Kepler Barycentric Julian Date ', fontsize=12, color='blue' )
	plt.ylabel( 'Flux', fontsize=12, color='blue' )
	plt.title( 'Lightcurve', color='blue' )
	plt.grid(True)
	plt.savefig(filename)
	return

def advancedplot( time, flux , num , filename): #advanced becasue it does not assume that the data has a particular format
	# Plots the data ---> slow 
	#only plots one in every num data points 
	#By advanced i mean that the input is more flexable - should add an input for title, xlabel and ylabel.
	datarange = range(1,len(time),num)
	timea = np.take(time, datarange, 0)
	fluxa = np.take(flux, datarange, 0)
	plt.clf() #make sure we have a clean slate
	plt.plot( timea, fluxa, 'b+')
	plt.xlabel( 'Kepler Barycentric Julian Date ', fontsize=12, color='blue' )
	plt.ylabel( 'Flux', fontsize=12, color='blue' )
	plt.title( 'Lightcurve', color='blue' )
	plt.grid(True)
	plt.savefig(filename)
	return

def transitmaster(STAR, data, foldername):
	dataarray = basicarraymaker(STAR, data) #make the structured array which will contain the lightcurve and its expanded data
	tdata, noplanets = transitdata(STAR[0]) #extract the transit data from the kepler mast data
	if noplanets == 0:
		planetfound = 'N' #if no planets are no - there are no transits to remove
	else:
		planetfound = 'Y' #not currently used for anything - probably pointless
		translength = transitlength(tdata,noplanets) #calculate the transit length using the data from kepler mast
		saveplandata(tdata, translength, noplanets, STAR[0], foldername)
		for i in range(0,noplanets,1) : #loop over all the planets found around a star
			plandataarray = folddataandtrim(dataarray, tdata, translength[i], i+1, STAR, foldername) #fold the lightcuve over a particular planets orbit
			plandataarray2 = labeltransits(plandataarray, tdata, translength, i+1, STAR) #labels the data points and numbers the transits
			intranarray = reducedarray(plandataarray2, 'Intran') #creates an array which only contains data points that fall within the transit
			neartranarray = reducedarray(plandataarray2, 'onedayoftran') #creates an array for the data points within 1 day of transit
			fulltranarray = np.concatenate((intranarray, neartranarray), axis=0) #now merges these arrays to create the array that will be normalised/plotted
			normfactors = normcalc(neartranarray) #calculate the normalisation factors for the transits
			fulltranlist = listmaker(fulltranarray) #create a list in which each entry is an array of data points near and in the transit
			ACTUALtranarray = normaliser(normfactors, fulltranlist) #normalise the transits and rejoin the data so that all transits can be plotted at once
			#advancedplot(ACTUALtranarray['Foldtime'], ACTUALtranarray['Flux'], 1, '{0}/kepler-{2}-planet{1}-tran.pdf'.format(foldername, i, STAR[0])) #A nice plot which can be used to confirm that the transit has been probably found
			np.save('{0}/kepler-{1}-planet-{2}-trandata'.format(foldername, STAR[0], i+1), ACTUALtranarray) #save the transit data as a numpy binary file
	return

def basicarraymaker(STAR, data):
	global dtypesone #becasue I will use these dtypes multiple times to recreate and move the array
	dtypesone = [('Time', float), ('Flux', float), ('Foldtime', float), ('Intran', bool), ('onedayoftran', bool), ('trannumber', int), ('datapointid', int)] #contains the data types of each coloumn in the structured array
	for i in range(0, len(data[0]), 1): #now loop over all the data points in the lightcuve
		values = [ (data[0][i], data[1][i], data[0][i], False, False, 0, i) ] #a list containing the data to be stored in the array, at this point 'foldtime' is just a copy of 'time'
		if not 'dataarray' in locals(): #if the output array has not been created, create it
			dataarray = np.array(values, dtype=dtypesone) #
		else: #otherwise
			temparray = np.array(values, dtype=dtypesone) #create a temporary array
			dataarray = np.concatenate((dataarray, temparray), axis=0) #and append the temporary array to the output array
	return dataarray #return the output array

def transitdata( KeplerID ):
	#loads the data into an array, only loads relevant columns!
	# 0=keplerID, 1=KOI, 2=Transit Duration, 5=time of first transit (JD), 7=Orbital period (days),
	# 9=(a/R_star), 11=(r/R_star), 13=b (impact parameter), 15=Radius of planet (Earth R),
	# 16=a (semi major axis) ---> for more info see http://archive.stsci.edu/kepler/planet_candidates.html
	keplerdat = np.loadtxt('kepler2000++.dat', skiprows=1, usecols=(0,1,2,5,7,9,11,13,15,16))
	#defines kp and keplerdat as the same thing, makes lines shorter!
	kp = keplerdat
	# A count which is used to keep track of the number of exoplanets around a star
	count = 0
	#produces the array for the 'for loop' to loop over - current data set contains 2300 entries
	linerange = range(2300)
	#The for loop that searches for entries with an KeplerID matching the KeplerID of interest
	for i in linerange:
		#KeplerId should be an integer but numpy arrays contain floats
		testid = int(keplerdat[i][0])
		#checks each entries keplerid
		if KeplerID==testid:
			if count==0:
				#Creates an arrary when the first entry for the keplerID is found
				retdata = np.array([[kp[i][0],kp[i][1],kp[i][2],kp[i][3],kp[i][4],kp[i][5],kp[i][6],kp[i][7],kp[i][8],kp[i][9]]])
			elif count>0:
				#appends additional entries, for additional planets, to the above array
				retdata = np.append(retdata, [[kp[i][0],kp[i][1],kp[i][2],kp[i][3],kp[i][4],kp[i][5],kp[i][6],kp[i][7],kp[i][8],kp[i][9]]], axis=0)
			#increase the count when an exoplanet is identified around the star
			count +=1
	#Sets retdata to 0 if no entry is found - since we will always be using an old data list!
	if count==0:
		retdata = 0
	# returns an array containing the relevant planets and the number of planets orbiting the star
	return retdata, count

def transitlength( plandat, plancnt ): # calculates the duration of the planetary transit's 
	#define the array for the for loop to operate over
	planrange = range(0,plancnt,1)
	for i in planrange:
		#calculate the stellar radius
		rstar = (plandat[i][9]*149.598E9/plandat[i][5])
		#convert the radius of the planet into meters
		rplanet = (plandat[i][8]*6.378E6)
		#calculate the impact parameter
		b = plandat[i][7] * rstar
		#calculate arcsin[sqrt{(rstar+rplanet)^2-b^2} / a]
		parameter = np.arcsin(math.fabs((math.sqrt(math.fabs(((rstar+rplanet)**2)-(b**2))) / (plandat[i][9]*149.598E9))))
		#transit duration will be in days
		tdur = (plandat[i][4]/ math.pi)*parameter
		if i==0: #create an array for the first entry
			transdat = [tdur]
		elif i>0: #append additional planets to the data array (1D-Array)
			transdat.append(tdur)
	return transdat #return the transit durations

def saveplandata(tdata, translength, noplanets, KeplerID, foldername): #this function will save the planetary data
	dtypesplan = [('KeplerID', int), ('KOI', float), ('Transitduration', float), ('planetnumber', int), ('Tfirsttran', float), ('Orbitalperiod', float),\
		('a/R_star', float), ('R_plan/R_star', float), ('Impactpar', float), ('R_planet', float), ('semimajoraxis', float)] #data types for the planetary array
	for i in range(0, len(tdata), 1): #loop over all planets for which data exists
		values = [(tdata[i][0], tdata[i][1], (tdata[i][2] / 24), i+1, tdata[i][3], tdata[i][4], tdata[i][5], tdata[i][6], tdata[i][7], tdata[i][8], tdata[i][9])] #values to be stored in the planetary data array
		if not 'planarray' in locals(): #if the planetary data array has not been created, create it
			planarray = np.array(values, dtype = dtypesplan)
		else: #otherwise append further planets to the planetary data array
			temparray = np.array(values, dtype = dtypesplan)
			planarray = np.concatenate((planarray, temparray), axis=0)
	#print planarray #testing
	if os.path.exists('{0}/plandataarray.npy'.format(foldername)) == False: #If the saved planetary data array does not yet exist, create an array to contain this data
		plandataarray = np.copy(planarray)
	else: #otherwise load the saved planetary data array and append the new planetary data to it
		plandataarray = np.load('{0}/plandataarray.npy'.format(foldername))
		plandataarray = np.concatenate((plandataarray, planarray), axis=0)
	#print plandataarray #Testing
	np.save('{0}/plandataarray'.format(foldername, KeplerID), plandataarray) #save the updated/new planetary data array as a numpy binary file
	return


def folddataandtrim( spectra, tdata, tlength, planetnumber, STAR, foldername ):
	#loop over each of the data points
	spect = np.copy(spectra)
	for j in range(0, len(spect['Time']), 1):
		for i in range(0,5000,1):
			#if the data point is outside the orbital period and has a time greater than t0, remove 
			#1 orbital period.
			if spect['Foldtime'][j] > (tdata[planetnumber-1][3] + 67 + (tdata[planetnumber-1][4]/2)): #note: the Julian dates of the two sets of data differ by 67 days
				spect['Foldtime'][j] = spect['Foldtime'][j] - (tdata[planetnumber-1][4])
			#if the data point is outside the orbital period and has a time smaller than t0, add 
			#1 orbital period.
			elif spect['Foldtime'][j] < (tdata[planetnumber-1][3] + 67 - (tdata[planetnumber-1][4]/2)): #note: the Julian dates of the two sets of data differ by 67 days
				spect['Foldtime'][j] = spect['Foldtime'][j] + (tdata[planetnumber-1][4])
			else:
				break #end the loop when all data points are within the folded period
	#centre the fold around 0
	spect['Foldtime'] = spect['Foldtime'] - (tdata[planetnumber-1][3] + 67)
	return spect


def labeltransits(spect, tdata, tranlength, planetnumber, STAR):
	#this loop will label each data point wrt the tranists
	T_dur = tranlength[planetnumber-1]
	if (2*T_dur) > 0.5:
		Sideperiod = 2*T_dur
	else:
		Sideperiod = 0.5
	for j in range(0, len(spect['Foldtime']), 1): #loop over the lightcurve 
		if spect['Foldtime'][j] < tranlength[planetnumber-1]/2 and spect['Foldtime'][j] > (-tranlength[planetnumber-1]/2):
			spect['Intran'][j] = True #if the data point falls within the transit, mark the 'In Transit' label as True
		elif spect['Foldtime'][j] > tranlength[planetnumber-1]/2 and (spect['Foldtime'][j] - Sideperiod) < (tranlength[planetnumber-1]/2):
			spect['onedayoftran'][j] = True #if the data points is within 1 day of the end of the transit, mark the 'One day of tran' label as true
		elif spect['Foldtime'][j] < (-tranlength[planetnumber-1]/2) and (spect['Foldtime'][j] + Sideperiod) > (-tranlength[planetnumber-1]/2):
			spect['onedayoftran'][j] = True #if the data points is within 1 day of the end of the transit, mark the 'One day of tran' label as true (negative values)
		else:
			pass # otherwise do nothing
	count = 1
	#this loop will number each transit for easy analysis, it does this by searching for changes in the labels produced above
	for i in range(0, len(spect['Foldtime']), 1): 
		if spect['onedayoftran'][i-1] == False and spect['onedayoftran'][i] == True:
			spect['trannumber'][i] = count #if the start of the transit is detected (i.e. the start of the 1 day surroundings), start labeling the points (actually covered by the below if statement)
		elif spect['onedayoftran'][i] == True or spect['Intran'][i] == True:
			spect['trannumber'][i] = count #in the transit, label the data point with the transit number
		elif spect['onedayoftran'][i-1] == True and spect['onedayoftran'][i] == False  and spect['Intran'][i] == False and i != 0:
			count += 1 #increase the count by one at the end of an extended transit
		else:
			pass
	#print spect #testing
	return spect

def reducedarray(spect, colname):
	for i in range(0, len(spect['Foldtime']), 1): #loop over the main array
		if spect[colname][i] == True: #if the data point is part of the required group
			values = [ (spect['Time'][i], spect['Flux'][i], spect['Foldtime'][i], spect['Intran'][i],\
				 spect['onedayoftran'][i], spect['trannumber'][i], spect['datapointid'][i]) ] #store the values associated it 
			if not 'returnarray' in locals():
				returnarray = np.array(values, dtype=dtypesone) #in a new array
			else:
				temparray = np.array(values, dtype=dtypesone) #or by appending to the new array (i.e. this will create an array of data points that belong in the group represented by 'colname')
				returnarray = np.concatenate((returnarray, temparray), axis=0)
	return returnarray

def normcalc(neartranarray):
	tranarraylist = listmaker(neartranarray)
	#now we need to calculate the normalisation factors and create an array to store them 
	for i in tranarraylist: #loop through the list of transit arrays
		transitnumber = i['trannumber'][0] #just makes it simpler to know which transit we are dealing with
		fit = np.polyfit((i['Foldtime']), i['Flux'], 1) #caluclate the first order (straight line) fit
		if not 'returnarray' in locals(): #check to see if the return array has been created
			returnarray=np.array([(transitnumber, fit[0], fit[1])]) #if not, create the array
		else: returnarray = np.append(returnarray, [(transitnumber, fit[0], fit[1])], axis = 0) #if so, append to the array
	#print returnarray #testing
	return returnarray

def listmaker(array): #a simple function which makes a list in which each entry is the data for a particular transit
	#print array
	for i in range(1, (max(array['trannumber'])+1), 1): #loop over all the transits
		for j in range(0, len(array['trannumber']), 1): #loop over all the data points
			if array['trannumber'][j] == i: #if the data point is in the required transit
				values = [ (array['Time'][j], array['Flux'][j], array['Foldtime'][j], array['Intran'][j],\
					array['onedayoftran'][j], array['trannumber'][j], array['datapointid'][j]) ] #store the values associated with it 
				if not 'storearray' in locals(): #either in a new array (in the first case)
					storearray = np.array(values, dtype=dtypesone)
				else: #or by appending to the above array
					temparray = np.array(values, dtype=dtypesone)
					storearray = np.concatenate((storearray, temparray), axis=0)		
		if not 'arraylist' in locals(): #now it checks if the list which will contain the arrays has been created, and if not, it creates the list
			arraylist = [storearray]
		else: #if the list has been created, the array is appended to the list
			arraylist.append(storearray)
		del(storearray) #deletes the current array, so that the new array can be created by using a simple existance check
	return arraylist

def normaliser(normfactors, tranarraylist):
	for i in range(0, int(normfactors[-1][0])): #loop over the transits (remember the first column of the normfactors is the transit number)
		thistranarray = tranarraylist[i] #extract the correct element from the list
		thistranarrayreturn = np.copy(thistranarray) #make a copy of the array - this fixes a really weird bug with numpy
		for j in range(0, len(thistranarray['Time']), 1): #now loop over the data points
			thistranarrayreturn['Flux'][j] = thistranarray['Flux'][j] / ((normfactors[i][1]*(thistranarray['Foldtime'][j])) + normfactors[i][2]) #and normalise the data points to the fitted line
		#print thistranarray #Testing
		if not 'fulltranarray' in locals(): #now we want to merge all the transits data into a new array, which will be created and appended to in the normal way
			fulltranarray = np.copy(thistranarrayreturn)
		else:
			fulltranarray = np.concatenate((fulltranarray, thistranarrayreturn), axis=0)
	#print fulltranarray #Testing
	return fulltranarray


batchmaster()