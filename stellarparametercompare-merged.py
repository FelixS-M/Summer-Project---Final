#Part of the Bound (Realistic stars) data pipeline

#A rather large piece of code which produces a number of lists which can be used to decide which stars are false positives
#1) Produces lists for Mstar, Rstar, Rhostar, Teff, Logg and Metalicity which can be used to compare the pre and post fit values - the lists are ordered by the Fractional diffrence between the fit and kepler parameters
#2) Produces a list of chi^2reduced values which is ordered by the difference between the actual chisquaredreduced value and a perfect fit value of one
#3) Produces a list of orbital inclinations, ordered by the inclination in degrees
#4) Using the impact parameters as a metric, this produces two lists of stars, those with a normalised impact parameters of greater than 0.85, which are excluded from the sample as probable false positives
#	and those with an impact parameter of less then 0.85, which are fine. It also produces a list which contains entries from both. All Lists are ordered by normalised impact parameter.

import numpy as np
import os
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import math


def TLMRMrhocompareloop(): #1
	STARDAT = np.load('Stellardataarray(merged).npy')
	#PlANDAT = np.load('plandataarray(fitted).npy')
	for i in ['Mcalc', 'Rcalc', 'Densitycalc', 'Teff', 'Logg', 'Metalicity']: #loop over all the lists to be made	
		for star in STARDAT: #loop over all the stars/planets
			Mcalc_Kepler, Rcalc_Kepler, Densitycalc_Kepler = RMrhocalc(star['Teff_kep'], star['Logg_kep'], star['Metalicity_kep']) #Calculate Mass, Radius and Density
			dtypes = [('KeplerID', int), ('planetnumber', float), ('{}_fit'.format(i), float), ('{}_kep'.format(i), float), ('FracDifference', float), ('Ratio', float)] #Data types for the final array
			#These if statements will ensure that the correct variables are loaded for the array/comparison list
			if i == 'Mcalc':
				fittedparameter = star['Mstar(SI)']
				keplerparameter = Mcalc_Kepler
			elif i == 'Rcalc':
				fittedparameter = star['Rstar(SI)']
				keplerparameter = Rcalc_Kepler
			elif i == 'Densitycalc':
				fittedparameter = star['Rhostar(SI)']
				keplerparameter = Densitycalc_Kepler
			elif i == 'Teff':
				fittedparameter = star['Teff_fit']
				keplerparameter = star['Teff_kep']
			elif i == 'Logg':
				fittedparameter = star['Logg_fit']
				keplerparameter = star['Logg_kep']
			elif i == 'Metalicity':
				fittedparameter = star['Metalicity_fit']
				keplerparameter = star['Metalicity_kep']	
			#End of if statments
			ratio = fittedparameter / keplerparameter #calculate the ratio of the fitted and original kepler header pararmeters
			fracdifference = np.abs(1 - ratio) #calculate the absolute value of the fractional difference between fitted and kepler parameters 
			values = [(star['KeplerID'], star['planetnumber'], fittedparameter, keplerparameter, fracdifference, ratio)] #values to be stored in the sturctured output array
			dataarray = np.array(values, dtype=dtypes) #create a temporary array for this stars data
			if os.path.exists('{}dataarray(comparison).npy'.format(i)) == False: #If the saved stellar parameter data array does not yet exist, create an array to contain this data
				outdataarray = np.copy(dataarray)
			else: #otherwise load the saved stellar parameter data array and append the new stellar parameter data to it
				outdataarray = np.load('{}dataarray(comparison).npy'.format(i))
				outdataarray = np.concatenate((outdataarray, dataarray), axis=0)
			sortedarray = np.sort(outdataarray, order = 'FracDifference') #Sort the array by the fractional difference between the fitted and kepler values
			np.save('{}dataarray(comparison)'.format(i), sortedarray) #save the sorted structured array
		print i #Testing
		print sortedarray #Testing
	return

def RMrhocalc(Teff, Logg, Metalicity): #DONE - add comments 
	Mcalc = calculateM(Teff, Logg, Metalicity) #Calculate the mass of the star, in solar masses
	Rcalc = calculateR(Teff, Logg, Metalicity) #Calculate the radius of the star, in solar radii
	Density, Density_SI = densitycalc(Rcalc, Mcalc) #calculate the stellar density
	Mcalc_SI = Mcalc * 1.989E30 #convert the mass to SI units
	Rcalc_SI = Rcalc * 6.955E8 #convert the radius to SI units
	return Mcalc_SI, Rcalc_SI, Density_SI

def calculateM(Teff, Logg, Metalicity): #DONE
	if Teff == 0: #This means thats Teff was not found in the header of the file
		return 0
	else:
		avalues = [1.5689, 1.3787, 0.4243, 1.139, -0.1425, 0.01969, 0.1010] #coefficients for the equation 
		#errors = [0.058, 0.029, 0.029, 0.24, 0.011, 0.0019, 0.014]
		X = math.log10(Teff) - 4.1 #makes it simpler to write the equation
		logM = avalues[0] + avalues[1]*X + avalues[2]*(X**2) + avalues[3]*(X**3) + avalues[4]*(Logg**2) + avalues[5]*(Logg**3) + avalues[6]*Metalicity
		Mcalc = 10**logM #remove the log
		#Merror = errorcalc(Teff, Logg, Metalicity, avalues, errors, Mcalc)
	return Mcalc

def calculateR(Teff, Logg, Metalicity): #DONE
	if Teff == 0: #This means thats Teff was not found in the header of the file
		return 0
	else:
		bvalues = [2.4427, 0.6679, 0.1771, 0.705, -0.21415, 0.02306, 0.04173] #coefficients for the equation 
		#errors = [0.038, 0.016, 0.027, 0.13, 0.0075, 0.0013, 0.0082]
		X = math.log10(Teff) - 4.1 #makes it simpler to write the equation
		logR = bvalues[0] + bvalues[1]*X + bvalues[2]*(X**2) + bvalues[3]*(X**3) + bvalues[4]*(Logg**2) + bvalues[5]*(Logg**3) + bvalues[6]*Metalicity
		Rcalc = 10**logR #remove the log
		#Rerror = errorcalc(Teff, Logg, Metalicity, bvalues, errors, Rcalc)
	return Rcalc

def densitycalc(Radius, Mass): #DONE
	if Radius == 0: #If no radius is calculated - no density can be calculated 
		return 0, 0
	else: #otherwise calculate the stellar density 
		density_SI = (3 * Mass * 1.989E30) / (4 * math.pi * ((Radius * 6.955E8)**3))
		density = (3 * Mass) / (4 * math.pi * (Radius**3))
	return density, density_SI

def chilistmaker(): #Create an array/list which is ordered by the size of the reduced chi^2 values #2)
	STARdata = np.load('Stellardataarray(merged).npy') #Load the stellar data
	for star in STARdata: #loop over all the stars / planets
		dtypes = [('KeplerID', int), ('planetnumber', float), ('chi2reduced', float), ('difference', float)] #data types for the chi data array 
		values = [(star['KeplerID'], star['planetnumber'], star['chi2reduced'], np.abs(star['chi2reduced'] - 1))] #values for the chi data array
		dataarray = np.array(values, dtype=dtypes) #Create a temporary chi data array for the current planet/star
		if os.path.exists('chi2reddataarray(comparison).npy') == False: #If the saved chi data array does not yet exist, create an array to contain this data
			outdataarray = np.copy(dataarray)
		else: #otherwise load the saved chi data array and append the new chi data to it
			outdataarray = np.load('chi2reddataarray(comparison).npy')
			outdataarray = np.concatenate((outdataarray, dataarray), axis=0)
		sortedarray = np.sort(outdataarray, order = 'difference') #sort the array by the difference between the chi^2 reduced value and the ideal sitution of a ch^2 reduced of 1
		np.save('chi2reddataarray(comparison)', sortedarray)
	print 'chi2reduced(intran)' #Testing
	print sortedarray #Tesing
	for star in STARdata: #loop over all the stars / planets
		dtypes = [('KeplerID', int), ('planetnumber', float), ('chi2reduced(full)', float), ('difference', float)] #data types for the chi data array 
		values2 = [(star['KeplerID'], star['planetnumber'], star['chi2reducedfull'], np.abs(star['chi2reducedfull'] - 1))] #values for the chi data array
		dataarray2 = np.array(values2, dtype=dtypes) #Create a temporary chi data array for the current planet/star
		if os.path.exists('chi2reddataarray(comparison)[full].npy') == False: #If the saved chi data array does not yet exist, create an array to contain this data
			outdataarray2 = np.copy(dataarray2)
		else: #otherwise load the saved chi data array and append the new chi data to it
			outdataarray2 = np.load('chi2reddataarray(comparison)[full].npy')
			outdataarray2 = np.concatenate((outdataarray2, dataarray2), axis=0)
		sortedarray2 = np.sort(outdataarray2, order = 'difference') #sort the array by the difference between the chi^2 reduced value and the ideal sitution of a ch^2 reduced of 1
		np.save('chi2reddataarray(comparison)[full]', sortedarray2)
	print 'chi2reduced(full)' #Testing
	print sortedarray2 #Tesing
	return

def chifitjudge():
	STARdata = np.load('Stellardataarray(merged).npy') #load the stellar dataarray
	for star in STARdata: #loop over all the planets/stars
		chi2cap = (1 + math.sqrt(2.0 / star['DOF'])) #calculate the maximum chi^2reduced value for a good fit
		if star['chi2reduced'] < chi2cap: #if the chi^2 reduced value is below this value
			withincap = True #good fit 
		else: #otherwise
			withincap = False
		dtypes = [('KeplerID', int), ('planetnumber', float), ('chi2reduced', float), ('chi2cap', float), ('Withincap', bool)] #data types for the output array
		values = [(star['KeplerID'], star['planetnumber'], star['chi2reduced'], chi2cap, withincap)] #data to be stored in the output array
		dataarray = np.array(values, dtype=dtypes) #create an array for this particular star/planet
		if os.path.exists('chi2caparray.npy') == False: #if the saved output array has not been created, create it
			outdataarray = np.copy(dataarray)
		else: #otherwise load the saved array and append the new result to it
			outdataarray = np.load('chi2caparray.npy')
			outdataarray = np.concatenate((outdataarray, dataarray), axis=0)
		sortedarray = np.sort(outdataarray, order = 'Withincap') #sort the array that is about to be saved by the indicatior 'Withincap'
		np.save('chi2caparray', sortedarray) #save the array
	print 'chi2cap(intran)' #Testing
	print sortedarray #Testing
	#This loop is the same as above, but instead of only using the intransit data, it uses the full (full as in inclulde in output pdf's) data. Code is essentially the same
	for star in STARdata:
		chi2cap2 = 1 + math.sqrt(2.0 / star['DOFfull'])
		if star['chi2reducedfull'] < chi2cap2:
			withincap2 = True
		else:
			withincap2 = False
		dtypes = [('KeplerID', int), ('planetnumber', float), ('chi2reduced(full)', float), ('chi2capfull', float), ('Withincap', bool)]
		values2 = [(star['KeplerID'], star['planetnumber'], star['chi2reducedfull'], chi2cap2, withincap2)]
		dataarray2 = np.array(values2, dtype=dtypes)
		if os.path.exists('chi2caparray(full).npy') == False:
			outdataarray2 = np.copy(dataarray2)
		else:
			outdataarray2 = np.load('chi2caparray(full).npy')
			outdataarray2 = np.concatenate((outdataarray2, dataarray2), axis=0)
		sortedarray2 = np.sort(outdataarray2, order = 'Withincap')
		np.save('chi2caparray(full)', sortedarray2)
	print 'chi2cap(full)' #testing
	print sortedarray2 #testing
	return

def chifit1point1judge():
	#This funcition is essentially (almost exactly) the same as 'chifitjudge', except the upper limit for a good fit is always set to 1.1
	STARdata = np.load('Stellardataarray(merged).npy')
	for star in STARdata:
		chi2cap = 1.1
		if star['chi2reduced'] < chi2cap:
			withincap = True
		else:
			withincap = False
		dtypes = [('KeplerID', int), ('planetnumber', float), ('chi2reduced', float), ('chi2cap', float), ('Withincap', bool)]
		values = [(star['KeplerID'], star['planetnumber'], star['chi2reduced'], chi2cap, withincap)]
		dataarray = np.array(values, dtype=dtypes)
		if os.path.exists('chi2[1point1]array.npy') == False:
			outdataarray = np.copy(dataarray)
		else:
			outdataarray = np.load('chi2[1point1]array.npy')
			outdataarray = np.concatenate((outdataarray, dataarray), axis=0)
		sortedarray = np.sort(outdataarray, order = 'Withincap')
		np.save('chi2[1point1]array', sortedarray)
	print 'chi2[1point1](intran)'
	print sortedarray
	for star in STARdata:
		chi2cap2 = 1.1
		if star['chi2reducedfull'] < chi2cap2:
			withincap2 = True
		else:
			withincap2 = False
		dtypes = [('KeplerID', int), ('planetnumber', float), ('chi2reduced(full)', float), ('chi2capfull', float), ('Withincap', bool)]
		values2 = [(star['KeplerID'], star['planetnumber'], star['chi2reducedfull'], chi2cap2, withincap2)]
		dataarray2 = np.array(values2, dtype=dtypes)
		if os.path.exists('chi2[1point1]array(full).npy') == False:
			outdataarray2 = np.copy(dataarray2)
		else:
			outdataarray2 = np.load('chi2[1point1]array(full).npy')
			outdataarray2 = np.concatenate((outdataarray2, dataarray2), axis=0)
		sortedarray2 = np.sort(outdataarray2, order = 'Withincap')
		np.save('chi2[1point1]array(full)', sortedarray2)
	print 'chi2[1point1](full)'
	print sortedarray2
	return


def orbitalinclinationlist(): #Create an array / list which is ordered by orbital inclination #3
	PLANdata = np.load('plandataarray(merged).npy') #load the planetary data array
	for planet in PLANdata: #loop over all the planets
		inclination_rad = planet['Inclination'] #extract the inclination (in radians)
		inclination_rad_error = planet['Inclination_error'] #extract the error in the inclination (in radians)
		inclination_degrees = math.degrees(inclination_rad) #convert the inclination to degrees
		dtypes = [('KeplerID', int), ('planetnumber', float), ('Inclination_rad', float), ('Inclination_rad_error', float), ('Inclination_degrees', float)] #data types for the inclination array
		values = [(planet['KeplerID'], planet['planetnumber'], inclination_rad, inclination_rad_error, inclination_degrees)] #values to be stored in the inclination array
		dataarray = np.array(values, dtype=dtypes) #create an entry which will be added to the inclination array for the current planet
		if os.path.exists('inclinationdataarray.npy') == False: #If the saved inclination data array does not yet exist, create an array to contain this data
			outdataarray = np.copy(dataarray)
		else: #otherwise load the saved inclination data array and append the new inclination data to it
			outdataarray = np.load('inclinationdataarray.npy')
			outdataarray = np.concatenate((outdataarray, dataarray), axis=0)
		sortedarray = np.sort(outdataarray, order = 'Inclination_degrees') #sort the inclination array by inclination
		np.save('inclinationdataarray', sortedarray) #and save it
	print 'Inclination' #Testing
	print sortedarray #Testing
	return

def impactparameterlist(): #4)
	PLANdata = np.load('plandataarray(merged).npy')
	for planet in PLANdata: #loop over all the planets
		impactparameter = planet['Impactpar'] #extract the impactparameter
		if impactparameter > 0.85: #if the impact parameter is greater than 0.85, we want to exclude it from our sample
			exclude = True
		else:
			exclude = False
		dtypes = [('KeplerID', int), ('planetnumber', float), ('Impactparameter', float), ('Exclude', bool)] #data types for the output array
		values = [(planet['KeplerID'], planet['planetnumber'], impactparameter, exclude)] #values to be stored in the output array
		dataarray = np.array(values, dtype=dtypes) #create an array for the current planet
		if os.path.exists('impactpardataarray.npy') == False: #if there is not already a saved impactparameter data array, create one
			outdataarray = np.copy(dataarray)
		else: #otherwise load the saved array and append the current planets data to it
			outdataarray = np.load('impactpardataarray.npy')
			outdataarray = np.concatenate((outdataarray, dataarray), axis = 0)
		sortedarray2 = np.sort(outdataarray, order = 'Impactparameter') #sort the array by impactparameter
		np.save('impactpardataarray', sortedarray2) #and save it
		#self explanetary - produces 2 arrays depending upon if Exculde is true or false 
		if dataarray['Exclude'] == True:
			if not 'excludearray' in locals():
				excludearray = np.copy(dataarray)
			else:
				addarray = np.copy(dataarray)
				excludearray = np.concatenate((excludearray, addarray), axis=0)
		else:
			if not 'keeparray' in locals():
				keeparray = np.copy(dataarray)
			else:
				addarray = np.copy(dataarray)
				keeparray = np.concatenate((keeparray, addarray), axis=0)
	sortedkeeparray = np.sort(keeparray, order = 'Impactparameter')
	sortedexcludearray = np.sort(excludearray, order = 'Impactparameter')
	print 'Impactparameter' #testing
	print sortedarray2 #testing
	#now need to produce the keep and discard arrays
	np.save('impactparkeeparray', sortedkeeparray) #save the keep array
	np.save('impactpardiscardarray', sortedexcludearray) #save the exclude/discard array
	print 'Impactparameter keep array' #Testing
	print sortedkeeparray #Testing
	print 'Impactparameter discard array' #Testing
	print sortedexcludearray #Testing
	return

def main():
	foldername = raw_input('Please enter the name of the folder in which contains the data: ') 
	os.chdir('{}'.format(foldername)) 
	TLMRMrhocompareloop()
	orbitalinclinationlist()
	impactparameterlist()
	chilistmaker()
	chifitjudge()
	chifit1point1judge()
	return

main()