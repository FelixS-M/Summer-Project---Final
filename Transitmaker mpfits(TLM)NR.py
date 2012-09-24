#Part of the unbound (free parameters) fit data pipeline


#This is the Second Stage in the data processing pipline: It fits an analtical lightcurve to each transit
#This script stars by load the stellar and planetary parameters and using them to calculate an inital fit: Each fit requires Teff, Logg, Metalicity, Planetary/Stellar radius ratio, impact parameter and the transit duration 
#  - it is these values that are iterated over inorder to find the best fit.
# Becasue of the way that this fit works, all the parameters are free and so it can produce the best possible fit (given the limitations of mpfits) for each transit. I.e it does not worry about ensuring a realsitic host star
#This script uses mpfits for the iteration, which means that errors in each fitted parameters are returned. It also means that the quality of each fit is judged by its chi^2 value
#This script also saves the fitted and original stellar parameters in new structured arrays along with a seperate file which contains the analytical (fitted) lightcurve 
#NOTE: The way that the z values are calculate is very simple, but it is geometrically and physically accurate.

#NOTE: This only differs from it's partner scipt in that it will NOT reprocess anystars which already have saved analytical transit lightcuves - It is used when dealing with large samples
#The code can be interupted (command-c) and them resumed with very little performance loss - perfect for very large samples

import shutil
import numpy as np
import transit
from mpfit import mpfit
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import os
import warnings
from scipy.interpolate import griddata
warnings.simplefilter('ignore', RuntimeWarning) #The analytical lightcurve generator makes runtime errors but still works - so might as well suppress the errors



def transittest(): # was used to test the analytical transit generator - now not used
	z = np.linspace(0,1.4,200)
	gammaslist = [[0., 0.], [1., 0.], [2., -1.]]
	plt.figure()
	for gammas in gammaslist:
		f = transit.occultquad(z, 0.1, gammas)
		plt.plot(z, f)
	plt.savefig('Test.pdf')
	return

def zandtcalculator(tvalues, impactparameter, transitduration, R_p_over_R_s):
	#first calculate the minimum z value - turns out impact parameter is already a fraction of stellar radius
	zmin = impactparameter
	#print zmin #tetsing
	#and then calculate the max z value, which occurs as transitduration / 2
	zmax = 1 + R_p_over_R_s
	tmax = transitduration / 2
	#now caluclate the relationship between t and z
	z_t = (zmax - zmin) / tmax
	#now calculate the corresponding z values
	zvalues = np.copy(tvalues) #make a copy of the tvalues array to store the z values in
	for i in range(0, len(tvalues), 1): #now calculate the zvalue associated with each t, remembering that z cannot be negative
		zvalues[i] = np.abs(tvalues[i] * z_t) + zmin
	return zvalues

def stellarextractor(KeplerID, Stellardataarray): #a simple function which returns Teff, Logg and Metalicity associated with a particular keplerID
	for line in Stellardataarray: #loop through the array which contains the values
		if line['KeplerID'] == KeplerID: #when the entry for the given keplerId is found 
			Teff  = line['Teff'] #extract Teff, 
			Logg = line['Logg']
			Metalicity = line['Metalicity']
			Rheader = line['StellarRadius']
			break #end the loop
	return Teff, Logg, Metalicity, Rheader #and return them

def transitplot(t, f, KeplerID, planetnumber):
	plt.clf() #clear all figures
	plt.figure() #create a figure for this plot
	plt.plot(t, f) #plot the data
	plt.xlabel('Time From Mid-Transit (Days)')
	plt.ylabel('Normalised Flux')
	plt.title('Analytical Transit Lightcuves')
	plt.savefig('kepler-{0}-planet-{1}-analyticaltransitlc.pdf'.format(KeplerID, planetnumber))
	return

def comtransitplot(t, f, kt, kf, KeplerID, planetnumber, KOI):
	std = np.std(kf)
	plt.clf() #clear all figures
	plt.figure() #create a figure for this plot
	f = plt.plot( kt, kf, 'b+', t, f, 'r-', label = ['Kepler', 'Analytical']) #plot both sets of data
	plt.legend([f], ['Kepler', 'Analytical'] ) #adds a legend to the plot
	plt.xlabel('Time From Mid-Transit (Days)')
	plt.ylabel('Normalised Flux')
	plt.ylim([(1-6*std),(1+6*std)])
	plt.title('Kepler and Analytical Transit Lightcurves')
	plt.savefig('kepler-{0}-planet-{1}-comptranlc.pdf'.format(KeplerID, planetnumber))
	return

def transitcalcloop():
	plandataarray = np.load('plandataarray.npy') #load the planetary data array
	stellardataarray = np.load('Stellardataarray.npy') #load the stellar data array (No longer need to run LDCcalc beforehand
	for line in plandataarray: #now loop through the planetary data array, in which line represents a single planet
		print 'For Kepler {}, planet {}:'.format(line['KeplerID'], line['planetnumber'])
		if os.path.exists('kepler-{0}-planet-{1}-analyticaltransitlc.npy'.format(line['KeplerID'], line['planetnumber'])) == False:
			Teff, Logg, Metalicity, Rheader = stellarextractor(line['KeplerID'], stellardataarray) #extract the gamma values
			Keplerlc = np.load('kepler-{0}-planet-{1}-trandata.npy'.format(line['KeplerID'], line['planetnumber'])) #extract the kepler transit data
			tvalues =  np.sort(Keplerlc['Foldtime']) #create an array of times associated with each data point and sort them 
			Keplertransit = np.sort(Keplerlc, order = 'Foldtime') #sort the kepler data in the same way
			sigma = sigmacalc(Keplertransit)
			impactpar, T_dur, R_plan, Teff_fit, Logg_fit, Metalicity_fit,  bestfitarray, errors, covar = zandtfit(line, Keplerlc, Teff, Logg, Metalicity, sigma) #calculate the t and z values
			gamma1, gamma2, = gammainterpolator(Teff_fit, Logg_fit, np.log10(Metalicity_fit))
			zvalues = zandtcalculator(tvalues, impactpar, T_dur, R_plan)
			f = transit.occultquad(zvalues, R_plan, [gamma1, gamma2]) #calculate/find the analytical transit lightcurve
			#transitplot(tvalues, f, line['KeplerID'], line['planetnumber']) #plot the analytical transit lightcurve
			comtransitplot(tvalues, f, Keplertransit['Foldtime'], Keplertransit['Flux'], line['KeplerID'], line['planetnumber'], line['KOI']) #plot both the kepler and anaytical transit data on the same graph
			analyticaltransitlc = np.vstack((tvalues, f))
			chi, chireduced = chicalc(Keplertransit, analyticaltransitlc, sigma, T_dur)
			#print analyticaltransitlc
			print 'The new impactparameter is {}+-{}, the new duration is {}+-{} and the new radius is {}+-{}'. format(impactpar, errors[0], T_dur, errors[1], R_plan, errors[2])
			np.save('kepler-{0}-planet-{1}-analyticaltransitlc'.format(line['KeplerID'], line['planetnumber']), analyticaltransitlc) #save the analytical transit
			np.save('kepler-{0}-planet-{1}-covariance'.format(line['KeplerID'], line['planetnumber']), covar) #save the covariance array
			saveplandata(line, impactpar, T_dur, R_plan, line['planetnumber'], line['KeplerID'], errors)
			savefittedstellardata(line['KeplerID'], Teff, Teff_fit, errors[3], Logg,  Logg_fit, errors[4], Metalicity, Metalicity_fit, errors[5], Rheader, gamma1, gamma2, line['planetnumber'])
		else: print 'Already Done'
	return

def saveplandata(tdata, impactpar, translength, R_plan, plannumber, KeplerID, errors): #this function will save the planetary data
	dtypesplan = [('KeplerID', int), ('KOI', float), ('Transitduration', float), ('Transitsuration_error', float), ('planetnumber', int), ('Tfirsttran', float), ('Orbitalperiod', float),\
		('a/R_star', float), ('R_plan/R_star', float), ('R_plan_error', float), ('Impactpar', float), ('Impactpar_error', float), ('semimajoraxis', float)] #data types for the planetary array
	values = [(tdata['KeplerID'], tdata['KOI'], translength, errors[1], plannumber, tdata['Tfirsttran'], tdata['Orbitalperiod'], tdata['a/R_star'], R_plan, errors[2], impactpar, errors[0], tdata['semimajoraxis'])] #values to be stored in the planetary data array
	planarray = np.array(values, dtype = dtypesplan)
	#print planarray #testing
	if os.path.exists('plandataarray(fitted).npy') == False: #If the saved planetary data array does not yet exist, create an array to contain this data
		plandataarray = np.copy(planarray)
	else: #otherwise load the saved planetary data array and append the new planetary data to it
		plandataarray = np.load('plandataarray(fitted).npy')
		plandataarray = np.concatenate((plandataarray, planarray), axis=0)
	#print plandataarray #Testing
	np.save('plandataarray(fitted)', plandataarray) #save the updated/new planetary data array as a numpy binary file
	return

def savefittedstellardata(KeplerID, Teff_kep, Teff_fit, Teff_error, Logg_kep, Logg_fit, Logg_error, Metalicity_kep, Metalicity_fit, Metalicity_error, Rheader, gamma1, gamma2, planetnumber):
	dtypestar = [('KeplerID', int), ('Teff_kep', float), ('Teff_fit', float), ('Teff_error', float), ('Logg_kep', float), ('Logg_fit', float), ('Logg_error', float),\
		('Metalicity_kep', float), ('Metalicity_fit', float), ('Metalicity_error', float), ('Rheader', float), ('gamma1', float), ('gamma2', float), ('planetnumber', float)] #data types for the stellar data array
	values = [(KeplerID, Teff_kep, Teff_fit, Teff_error, Logg_kep, Logg_fit, Logg_error, Metalicity_kep, Metalicity_fit, Metalicity_error, Rheader, gamma1, gamma2, planetnumber)]
	STARarray = np.array(values, dtype = dtypestar)
	#print planarray #testing
	if os.path.exists('Stellardataarray(fitted).npy') == False: #If the saved planetary data array does not yet exist, create an array to contain this data
		STARdataarray = np.copy(STARarray)
	else: #otherwise load the saved planetary data array and append the new planetary data to it
		STARdataarray = np.load('Stellardataarray(fitted).npy')
		STARdataarray = np.concatenate((STARdataarray, STARarray), axis=0)
	#print plandataarray #Testing
	np.save('Stellardataarray(fitted)', STARdataarray) #save the updated/new planetary data array as a numpy binary file
	return

def zandtfit(plandata, Keplerlc, Teff, Logg, Metalicity, sigma):
	p0 = [plandata['Impactpar'], plandata['Transitduration'], (plandata['R_plan/R_star']/3), Teff, Logg, np.log10(Metalicity)] #inital guess (from kepler mast)
	#print p0
	fa = {'Keplerlc':Keplerlc, 'sigma':sigma} #additional variables to be past to the function
	parinfo = [{'value':p0[0], 'fixed':False, 'limited':[True, True], 'limits':[0,2]}, {'value':p0[1], 'fixed':False}, {'value':p0[2], 'fixed':False},\
		{'value':p0[3], 'fixed':False, 'limited':[True,True], 'limits':[0,40000]}, {'value':p0[4], 'fixed':False, 'limited':[True,True], 'limits':[0,5.0]},\
		{'value':p0[5], 'fixed':False, 'limited':[True,True], 'limits':[-5,1]} ]
	m = mpfit(minimisefunction, p0, functkw = fa, quiet = 0, maxiter = 25, parinfo = parinfo, ftol=0.0001) #run the mpfit for the best fit
	#print m #testing
	bestfitarray = m.params
	errors = m.perror #these need to be properly stored
	covar = m.covar #these need to be properly stored
	if errors == None:
		errors = [0,0,0,0,0,0,0]
	#print bestfitarray
	return np.abs(bestfitarray[0]), bestfitarray[1], bestfitarray[2], bestfitarray[3], bestfitarray[4], (10**bestfitarray[5]), bestfitarray, errors, covar

def minimisefunction(x, Keplerlc=None, fjac = None, sigma = None):
	impactparameter = x[0] #extract the impact parameter
	Transitduration = x[1] #extract the transit duration 
	R_plan = x[2] #extract the planetary radius 
	Teff = x[3]
	Logg = x[4]
	Metalicity = x[5]
	gamma1, gamma2 = gammainterpolator(Teff, Logg, Metalicity)
	tvalues =  np.sort(Keplerlc['Foldtime']) #create an array of times associated with each data point and sort them 
	Keplertransit = np.sort(Keplerlc, order = 'Foldtime') #sort the kepler data in the same way
	zvalues = zandtcalculator(tvalues, impactparameter, Transitduration, R_plan) #calculate the t and z values
	f = transit.occultquad(zvalues, R_plan, [gamma1, gamma2]) #calculate/find the analytical transit lightcurve
	#transitplot(tvalues, f, line['KeplerID'], line['planetnumber']) #plot the analytical transit lightcurve
	tempkeplervalues = np.copy(Keplertransit) #create a copy of the kepler data which will be used to caculate the standard deviation
	for i in range(0, len(Keplertransit['Flux']), 1): #loop over the kepler data array, normalising it to the analytical lightcuve 
		if not 'returnval' in locals():
			returnval = np.array((Keplertransit['Flux'][i] - f[i]) / sigma) #calculate the first value to return
		else:
			temparray = np.array((Keplertransit['Flux'][i] - f[i]) / sigma) #calculate the other values to return
			returnval = np.append(returnval, temparray) #and add them to the return array
	sigma = np.std(tempkeplervalues['Flux']) #calculate the standard deivation of the normalised lightcurve
	status = 0 #this will not error, so it exists purely as a formality 
	return [status, returnval] #return a list, which is apparently what is required

def sigmacalc(Keplerlc): #cauclates the deviation of all data points - essentailly an error value which can be used to calculate chi^2
	for i in range(0,len(Keplerlc), 1): #loops over all the data points iorder to create an array which only contains the near transit data
		if Keplerlc['onedayoftran'][i] == True:
			if not 'sigmacalcarray' in locals():
				sigmacalcarray = Keplerlc['Flux'][i]
			else:
				sigmacalcarray = np.append(sigmacalcarray, Keplerlc['Flux'][i])
	sigma = np.std(sigmacalcarray) #calculates the standard deviation of the near transit points 
	#print 'Sd is {}'.format(sigma) #testing 
	return sigma

def chicalc(Keplerlc, analyticallc, sigma, T_dur): #caculates the chi^2 and reduced chi^2 values whilst only using those data points that fall within a transit
	count = 0
	for i in range(0, len(Keplerlc['Foldtime']), 1): #loops over the kepler lc (and analytical lc) inorder to find each points chi^2 value
		if np.abs(Keplerlc['Foldtime'][i]) <= T_dur: #if the data point is part of the required group (i.e. falls within the transit)
			count += 1
			if not 'chi2' in locals(): #for the first point, create a chi^2 variable
				chi2 = (((Keplerlc['Flux'][i] - analyticallc[1][i])**2) / sigma**2)
			else: #otherwise add the next points chi^2 value to the total chi^2 value
				chi2 = chi2 + (((Keplerlc['Flux'][i] - analyticallc[1][i])**2) / sigma**2)
		else:
			continue
	if count == 0:
		return np.nan(), np.nan()
	print 'The chi^2 value is {} (intran)'.format(chi2) #Testing 
	chi2reduced = chi2 / (count + 3) #now calculate the reduced chi^2 value for 4 independent varaibles
	print 'The reduced chi^2 value is {} (intran)'.format(chi2reduced) #testing
	return chi2, chi2reduced #return both chi^2 values 

def gammainterpolator(Teff, Logg, M_H):
	#LDCarray = np.loadtxt('KeplerLDC.txt', skiprows=9, usecols = (0,1,2,4,5)) #loads the full array (but only with the quadratic LDC)
	#print LDCarray #Testing
	if not 'points'	in globals():
		global points
		global gammalist
		points = np.loadtxt('KeplerLDC.txt', skiprows=9, usecols=(0,1,2)) #load array of data points (Teff, Logg, M/H)
		gammalist = np.loadtxt('KeplerLDC.txt', skiprows=9, usecols=(4,5,5)) #load array of points to be interpolate (third point can be discarded but is required to prevent error)
	point = griddata(points, gammalist, (Teff, Logg, M_H), method = 'linear') #linearly interpolate to find the LDC for the given data point
	returnarray = np.array([point[0], point[1]]) #create an array to store the data points
	if np.isnan(returnarray[0]) == True: #if the 'linear' interpolation fails, use the nearest data point 
		point = griddata(points, gammalist, (Teff, Logg, M_H), method = 'nearest') #find the nearest data point (nearest gammas)
		returnarray = np.array([point[0], point[1]]) #create an array to store the data points
	#print returnarray #Testing
	return returnarray

def main():
	foldername = raw_input('Please enter the name of the folder in which contains the data: ')
	if os.path.exists('{}/KeplerLDC.txt') == False:
		currentdirectory = os.getcwd() #Get the current directory
		dst = '{}/{}/'.format(currentdirectory,foldername) #the destination for the merged data
		src = '{}/KeplerLDC.txt'.format(currentdirectory) #the destination for the merged data
		shutil.copy(src,dst)
	os.chdir('{}'.format(foldername)) #changes the current directory to the specified directory
	transitcalcloop()
	return

main()