#part of the unbound (ideal stars) pipeline

#a version of 'Transitmaker mpfit(Idealstars-SII)NR.py' which iterates any stars contained in the list 'flats.txt' far more slowly and carefully than normal - tends to fix the majority of the flat lightcurves  

import shutil
import numpy as np
import transit
import math
from mpfit import mpfit
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import os
import warnings
from scipy.interpolate import griddata
warnings.simplefilter('ignore', RuntimeWarning) #The analytical lightcurve generator makes runtime errors but still works - so might as well suppress the errors
np.set_printoptions(threshold='nan')



def zandtcalculator(tvalues, Stellardensity_SI, inclination, Orbitalperiod_SI): #DONE - add comments
	#Orbitalperiod_SI = Orbitalperiod_SI / 2 #FOR SOME REASON THIS FIXES A LOAD OF PROBLEMS
	Gravconst =  6.67300E-11
	tvalues2 = tvalues * 86400 #convert the times from days to seconds
	phasevalues = tvalues2 / (Orbitalperiod_SI) #convert the times to orbital phase
	zvalues = np.copy(tvalues) #create an array to store the zvalues
	M_over_R_cubed = Stellardensity_SI  * (4 * math.pi / 3) #note, this is not actually the stellar density
	for i in range(0, len(phasevalues), 1): #now calculate the zvalue associated with each t, using the formulae david gave me (i.e. the one i was originally going to use)
		zvalues[i] = (((Gravconst * M_over_R_cubed * (Orbitalperiod_SI**2)) / (4*math.pi*math.pi))**(1.0/3)) * math.sqrt(((math.sin(phasevalues[i] * math.pi * 2))**2) + ((math.cos(inclination)*math.cos(phasevalues[i] * math.pi * 2))**2))
	#print (((Gravconst * M_over_R_cubed * (Orbitalperiod_SI**2)) / (4*math.pi*math.pi))**(1.0/3)) * math.cos(inclination)
	return zvalues

def stellarextractor(KeplerID, Stellardataarray): #a simple function which returns Teff, Logg and Metalicity associated with a particular keplerID
	for line in Stellardataarray: #loop through the array which contains the values
		if line['KeplerID'] == KeplerID: #when the entry for the given keplerId is found 
			Teff  = line['Teff'] 
			Logg = line['Logg']
			Metalicity = line['Metalicity']
			Rheader = line['StellarRadius']
			break #end the loop
	if Teff == 0: #if no stellar data is available, assume initial conditions of a sol-like star
		Teff = 5778
	if Logg == 0:
		Logg = 1.45
	return Teff, Logg, Metalicity, Rheader #and return them

def transitplot(t, f, KeplerID, planetnumber): #DONE
	plt.clf() #clear all figures
	plt.figure() #create a figure for this plot
	plt.plot(t, f) #plot the data
	plt.xlabel('Time From Mid-Transit (Days)')
	plt.ylabel('Normalised Flux')
	plt.title('Analytical Transit Lightcuves')
	plt.savefig('kepler-{0}-planet-{1}-analyticaltransitlc.pdf'.format(KeplerID, planetnumber))
	return

def comtransitplot(t, f, kt, kf, KeplerID, planetnumber, KOI): #DONE
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

def transitcalcloop(): #Just need to add the save functions
	plandataarray = np.load('plandataarray.npy') #load the planetary data array
	stellardataarray = np.load('Stellardataarray.npy') #load the stellar data array (No longer need to run LDCcalc beforehand
	flats = np.loadtxt('flats.txt')
	for line in plandataarray: #now loop through the planetary data array, in which line represents a single planet
		print 'For Kepler {}, planet {}:'.format(line['KeplerID'], line['planetnumber'])
		if line['KeplerID'] in flats and os.path.exists('kepler-{0}-planet-{1}-analyticaltransitlc.npy'.format(line['KeplerID'], line['planetnumber'])) == False:
			Teff, Logg, Metalicity, Rheader = stellarextractor(line['KeplerID'], stellardataarray) #extract the gamma values
			Inclination = math.acos(line['Impactpar'] / line['a/R_star']) #calculate the inclination from the normalised impactpar
			#print Inclination
			Planetaryradius = line['R_planet'] #extract the planetary radius
			Planetaryradius_SI = Planetaryradius * 6.3675E6 #and convert it from earth radii to meters
			Orbitalperiod = line['Orbitalperiod'] #Extract the orbital period 
			Orbitalperiod_SI = Orbitalperiod * 86400 #and convert it from days to seconds
			Keplerlc = np.load('kepler-{0}-planet-{1}-trandata.npy'.format(line['KeplerID'], line['planetnumber'])) #extract the kepler transit data
			tvalues =  np.sort(Keplerlc['Foldtime']) #create an array of times associated with each data point and sort them 
			Keplertransit = np.sort(Keplerlc, order = 'Foldtime') #sort the kepler data in the same way
			sigma = sigmacalc(Keplertransit) #Calculate the standard deviation of the surrounding data points
			Inclination_fit, Planetaryradius_SI_fit, Teff_fit, Logg_fit, Metalicity_fit, Mstar_SI, Rstar_SI, Stellardensity_SI, semimajoraxis_SI, T_dur_SI, errors, covar = \
				zandtfit(Inclination, Planetaryradius_SI, Keplerlc, Teff, Logg, Metalicity, sigma, Orbitalperiod_SI, line['KeplerID'], line['planetnumber']) #calculate the t and z values
			gamma1, gamma2, = gammainterpolator(Teff_fit, Logg_fit, np.log10(Metalicity_fit)) #calculate the LDC
			zvalues = zandtcalculator(tvalues, Stellardensity_SI, Inclination_fit, Orbitalperiod_SI) #calculate the zvalues
			#print zvalues #testing
			#print tvalues
			#print 'The minimum zvalue is:'
			#print np.min(zvalues)
			norm_impactparameter = (semimajoraxis_SI * math.cos(Inclination_fit) / Rstar_SI)
			print 'The transit duration is:'
			print T_dur_SI / 86400
			Rplanratio = Planetaryradius_SI_fit / Rstar_SI #Calculate the ratio of planetary radius to stellar radius
			#print 'The Planetary Ratio is'
			print Logg, Logg_fit
			f = transit.occultquad(zvalues, Rplanratio, [gamma1, gamma2]) #calculate/find the analytical transit lightcurve
			##transitplot(tvalues, f, line['KeplerID'], line['planetnumber']) #plot the analytical transit lightcurve
			comtransitplot(tvalues, f, Keplertransit['Foldtime'], Keplertransit['Flux'], line['KeplerID'], line['planetnumber'], line['KOI']) #plot both the kepler and analytical transit data on the same graph
			analyticaltransitlc = np.vstack((tvalues, f)) #make a single array to contain the analytical transit lightcurve
			chi, chireduced, DOF, chifull, chireducedfull, DOFfull = chicalc(Keplertransit, analyticaltransitlc, sigma, T_dur_SI) #calculate the chisquared values
			#print analyticaltransitlc
			#print 'The new impactparameter is {}+-{}, the new duration is {}+-{} and the new radius is {}+-{}'. format(impactpar, errors[0], T_dur, errors[1], R_plan, errors[2])
			print 'Saving'
			np.save('kepler-{0}-planet-{1}-analyticaltransitlc'.format(line['KeplerID'], line['planetnumber']), analyticaltransitlc) #save the analytical transit
			np.save('kepler-{0}-planet-{1}-covariance'.format(line['KeplerID'], line['planetnumber']), covar) #save the covariance array
			saveplandata(line, Inclination_fit, T_dur_SI, Planetaryradius_SI_fit, semimajoraxis_SI, line['planetnumber'], line['KeplerID'], errors, chireduced, norm_impactparameter, DOF, chireducedfull, DOFfull)
			savefittedstellardata(line['KeplerID'], line['planetnumber'], Mstar_SI, Rstar_SI, Stellardensity_SI, Teff, Teff_fit, errors[2], Logg, Logg_fit, errors[3], Metalicity, Metalicity_fit, errors[4], Rheader, gamma1, gamma2, chireduced, DOF, chireducedfull, DOFfull)
		else: print 'Already Done'
	return

def saveplandata(tdata, Inclination, translength, Rplan_SI, semimajoraxis_SI, plannumber, KeplerID, errors, chireduced, norm_impactparameter, DOF, chireducedfull, DOFfull): #this function will save the planetary data
	#print tdata, Inclination, translength, Rplan_SI, semimajoraxis_SI, plannumber, KeplerID, errors
	dtypesplan = [('KeplerID', int), ('KOI', float),  ('planetnumber', float), ('Transitduration', float), ('Tfirsttran', float), ('Orbitalperiod', float),\
		('Rplan(SI)', float), ('R_plan_error', float), ('Inclination', float), ('Inclination_error', float), ('semimajoraxis', float), ('chi2reduced', float), ('DOF', int), ('chi2reducedfull', float), ('DOFfull', float), ('Impactpar', float)] #data types for the planetary array
	values = [(tdata['KeplerID'], tdata['KOI'], plannumber, translength, tdata['Tfirsttran'], tdata['Orbitalperiod'], Rplan_SI, errors[1], Inclination, errors[0], semimajoraxis_SI, chireduced, DOF, chireducedfull, DOFfull, norm_impactparameter)] #values to be stored in the planetary data array
	#print values
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

def savefittedstellardata(KeplerID, planetnumber, Mstar_SI, Rstar_SI, Stellardensity_SI, Teff_kep, Teff_fit, Teff_error, Logg_kep, Logg_fit, Logg_error, Metalicity_kep, Metalicity_fit, Metalicity_error, Rheader, gamma1, gamma2, chireduced, DOF, chireducedfull, DOFfull):
	#print KeplerID, planetnumber, Mstar_SI, Rstar_SI, Stellardensity_SI, Teff_kep, Teff_fit, Teff_error, Logg_kep, Logg_fit, Logg_error, Metalicity_kep, Metalicity_fit, Metalicity_error, Rheader, gamma1, gamma2
	dtypestar = [('KeplerID', int), ('planetnumber', float), ('Mstar(SI)', float), ('Rstar(SI)', float), ('Rhostar(SI)', float), ('Teff_kep', float), ('Teff_fit', float), ('Teff_error', float), ('Logg_kep', float), ('Logg_fit', float), ('Logg_error', float),\
		('Metalicity_kep', float), ('Metalicity_fit', float), ('Metalicity_error', float), ('Rheader', float), ('gamma1', float), ('gamma2', float), ('chi2reduced', float), ('DOF', int), ('chi2reducedfull', float), ('DOFfull', float),] #data types for the stellar data array
	values = [(KeplerID, planetnumber, Mstar_SI, Rstar_SI, Stellardensity_SI, Teff_kep, Teff_fit, Teff_error, Logg_kep, Logg_fit, Logg_error, Metalicity_kep, Metalicity_fit, Metalicity_error, Rheader, gamma1, gamma2, chireduced, DOF, chireducedfull, DOFfull)]
	STARarray = np.array(values, dtype = dtypestar)
	#print planarray #testing
	if os.path.exists('Stellardataarray(fitted).npy') == False: #If the saved stellar data array does not yet exist, create an array to contain this data
		STARdataarray = np.copy(STARarray)
	else: #otherwise load the saved stellar data array and append the new stellar data to it
		STARdataarray = np.load('Stellardataarray(fitted).npy')
		STARdataarray = np.concatenate((STARdataarray, STARarray), axis=0)
	#print plandataarray #Testing
	np.save('Stellardataarray(fitted)', STARdataarray) #save the updated/new stellar data array as a numpy binary file
	return

def zandtfit(Inclination, Planetaryradius_SI, Keplerlc, Teff, Logg, Metalicity, sigma, Orbitalperiod_SI, KeplerID, Planetnumberer):
	Planetaryradius_Earth = (Planetaryradius_SI) /  6.3675E6 #convert the planetary radius to earth radii
	p0 = [math.cos(1.570795), Planetaryradius_Earth, Teff, Logg, np.log10(Metalicity)] #inital guess (from kepler mast and headers)
	fa = {'Keplerlc':Keplerlc, 'sigma':sigma, 'Orbitalperiod_SI':Orbitalperiod_SI, 'KeplerID':KeplerID, 'Planetnumberer':Planetnumberer} #additional variables to be past to the function
	parinfo = [{'value':p0[0], 'fixed':0, 'limited':[1,1], 'limits':[-1,1], 'step':0.01}, {'value':p0[1], 'fixed':0, 'limited':[1,0], 'limits':[0.01,100], 'step':0},\
		{'value':p0[2], 'fixed':0, 'limited':[1,1], 'limits':[0,40000], 'step':0}, {'value':p0[3], 'fixed':0, 'limited':[1,1], 'limits':[0,5.0], 'step':0},\
		{'value':p0[4], 'fixed':0, 'limited':[1,1], 'limits':[-5,1], 'step':0} ]
	m = mpfit(minimisefunction, p0, functkw = fa, quiet = 0, maxiter = 35, parinfo = parinfo, ftol=0.0001) #run the mpfit for the best fit
	#print m #testing
	bestfitarray = m.params #extract the results 
	errors = m.perror #these need to be properly stored
	covar = m.covar #these need to be properly stored
	if errors == None: #i.e. if mpfits is a pain and decides to not return the errors or covariance properly
		errors = np.zeros(5) #create an empty errors array
	if covar == None: #if no covariance array is produced, create an array of zeros so that there are no problems later
		covar = np.zeros((5,5))
	#print bestfitarray
	Inclination_fit = math.acos(bestfitarray[0]) #Extract the fit inclination
	Planetaryradius_SI_fit =   bestfitarray[1] * 6.3675E6 #Extract the fit planetary radius and convert it to SI units
	Teff_fit = bestfitarray[2] #Extract the fit Teff
	Logg_fit = bestfitarray[3] #Extract the fit logg of the surface gravity
	Metalicity_fit = (10**bestfitarray[4]) #Extract the metalicity and remove the logg
	Mstar_SI, Rstar_SI, Stellardensity_SI = RMrhocalc(Teff_fit, Logg_fit, Metalicity_fit) #Calculate the stellar mass, radius and density
	semimajoraxis_SI = semimajoraxiscalc(Orbitalperiod_SI, Mstar_SI) #calculate the semi-major axis
	T_dur_SI = Transitdurationcalc(Rstar_SI, Planetaryradius_SI_fit, semimajoraxis_SI, Inclination_fit, Orbitalperiod_SI) #calculate the transit duration
	return Inclination_fit, Planetaryradius_SI_fit, Teff_fit, Logg_fit, Metalicity_fit, Mstar_SI, Rstar_SI, Stellardensity_SI, semimajoraxis_SI, T_dur_SI, errors, covar

def minimisefunction(x, Keplerlc=None, fjac = None, sigma = None, Orbitalperiod_SI = None, KeplerID = None, Planetnumberer = None): #DONE
	Inclination = math.acos(x[0]) #extract the impact parameter
	Planetaryradius = x[1] * 6.3675E6 #extract the transit duration 
	Teff = x[2] #Extract the effective surface temperature
	Logg = x[3] #Extract the log of the surface gravity
	Metalicity = x[4] #Extract the metalicity (logged)
	gamma1, gamma2 = gammainterpolator(Teff, Logg, Metalicity) #calculate the gamma values
	tvalues =  np.sort(Keplerlc['Foldtime']) #create an array of times associated with each data point and sort them 
	Keplertransit = np.sort(Keplerlc, order = 'Foldtime') #sort the kepler data in the same way
	Stellarmass_SI, StellarRadius_SI, Stellardensity_SI = RMrhocalc(Teff, Logg, (10**Metalicity))
	Planetaryradiusratio  = Planetaryradius / StellarRadius_SI #Planetary radius in terms of the stellar radius
	#print Teff, Logg, Metalicity, Planetaryradius, StellarRadius_SI, Planetaryradiusratio
	zvalues = zandtcalculator(tvalues, Stellardensity_SI, Inclination, Orbitalperiod_SI) #calculate the t and z values
	#print np.min(zvalues)
	f = transit.occultquad(zvalues, Planetaryradiusratio, [gamma1, gamma2]) #calculate/find the analytical transit lightcurve
	#transitplot(tvalues, f, line['KeplerID'], line['planetnumber']) #plot the analytical transit lightcurve
	#tempkeplervalues = np.copy(Keplertransit) #create a copy of the kepler data which will be used to calculate the standard deviation
	for i in range(0, len(Keplertransit['Flux']), 1): #loop over the kepler data array, normalising it to the analytical lightcurve 
		if not 'returnval' in locals():
			returnval = np.array((Keplertransit['Flux'][i] - f[i]) / sigma) #calculate the first value to return
		else:
			temparray = np.array((Keplertransit['Flux'][i] - f[i]) / sigma) #calculate the other values to return
			returnval = np.append(returnval, temparray) #and add them to the return array
	#sigma = np.std(tempkeplervalues['Flux']) #calculate the standard deviation of the normalised lightcurve
	status = 0 #this will not error, so it exists purely as a formality 
	if np.min(zvalues) > (1+Planetaryradiusratio):
		return [status, returnval*100]
	return [status, returnval] #return a list, which is apparently what is required

def semimajoraxiscalc(Orbitalperiod_SI, Mstar_SI): #DONE - Add Comments
	#Orbitalperiod_SI = Orbitalperiod_SI / 2
	Gravconst = 6.67300E-11 #Makes things easier in the equation
	semimajoraxis = ((Gravconst * Mstar_SI * (Orbitalperiod_SI**2)) / (4 * math.pi * math.pi))**(1.0/3) #calculate the semimajoraxis
	return semimajoraxis

def Transitdurationcalc(Rstar_SI, Rplan_SI, semimajoraxis_SI, inclination, Orbitalperiod_SI): #Done - Add Comments
	#Orbitalperiod_SI = Orbitalperiod_SI / 2
	try: #need to check that an error does not occur when a transit is not detected
		a = (Orbitalperiod_SI / math.pi) * math.asin(math.sqrt(np.abs((( (Rstar_SI + Rplan_SI) / semimajoraxis_SI)**2) - (math.cos(inclination)**2))))		
		del a #because we don't actually need it to be stored in memory
	except ValueError:
		return 0 #If the transit is not detected, return 0
	T_dur_SI = (Orbitalperiod_SI / math.pi) * math.asin(math.sqrt(np.abs((( (Rstar_SI + Rplan_SI) / semimajoraxis_SI)**2) - (math.cos(inclination)**2)))) #calculate the transit duration
	return T_dur_SI

def sigmacalc(Keplerlc): #cauclates the deviation of all data points - essentially an error value which can be used to calculate chi^2 #Still Should be ok
	for i in range(0,len(Keplerlc), 1): #loops over all the data points inorder to create an array which only contains the near transit data
		if Keplerlc['onedayoftran'][i] == True:
			if not 'sigmacalcarray' in locals():
				sigmacalcarray = Keplerlc['Flux'][i]
			else:
				sigmacalcarray = np.append(sigmacalcarray, Keplerlc['Flux'][i])
	sigma = np.std(sigmacalcarray) #calculates the standard deviation of the near transit points 
	#print 'Sd is {}'.format(sigma) #testing 
	return sigma

def chicalc(Keplerlc, analyticallc, sigma, T_dur_SI): #calculates the chi^2 and reduced chi^2 values whilst only using those data points that fall within a transit #Still should be ok
	T_dur = T_dur_SI / 86400 #need transit duration to be in days
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
	if count == 0: #if there are no intransit data points, then there are no chi values
		chi2 = 0
	print 'The chi^2 value is {} (intran)'.format(chi2) #Testing 
	chi2reduced = chi2 / (count + 5) #now calculate the reduced chi^2 value for 4 independent variables
	print 'The reduced chi^2 value is {} (intran)'.format(chi2reduced) #testing
	countfull = 0
	for i in range(0, len(Keplerlc['Foldtime']), 1): #loops over the kepler lc (and analytical lc) inorder to find each points chi^2 value
		countfull += 1
		if not 'chi2full' in locals(): #for the first point, create a chi^2 variable
			chi2full = (((Keplerlc['Flux'][i] - analyticallc[1][i])**2) / sigma**2)
		else: #otherwise add the next points chi^2 value to the total chi^2 value
			chi2full = chi2full + (((Keplerlc['Flux'][i] - analyticallc[1][i])**2) / sigma**2)
	if countfull == 0: #if there are no intransit data points, then there are no chi values
		chi2full = 0
	print 'The chi^2 value is {} (full)'.format(chi2full) #Testing 
	chi2reducedfull = chi2full / (countfull + 5) #now calculate the reduced chi^2 value for 4 independent variables
	print 'The reduced chi^2 value is {} (full)'.format(chi2reducedfull) #testing
	return chi2, chi2reduced, (count+5), chi2full, chi2reducedfull, (countfull+5) #return both chi^2 values 

def gammainterpolator(Teff, Logg, M_H): #Still should be ok - YUP
	#LDCarray = np.loadtxt('KeplerLDC.txt', skiprows=9, usecols = (0,1,2,4,5)) #loads the full array (but only with the quadratic LDC)
	#print LDCarray #Testing
	if not 'points'	in globals(): #load the arrays, since IO operations are time consuming 
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

def RMrhocalc(Teff, Logg, Metalicity): #DONE - add comments 
	Mcalc = calculateM(Teff, Logg, Metalicity) #Calculate the mass of the star, in solar masses
	Rcalc = calculateR(Teff, Logg, Metalicity) #Calculate the radius of the star, in solar radii
	Density, Density_SI = densitycalc(Rcalc, Mcalc) #calculate the stellar density
	Mcalc_SI = Mcalc * 1.989E30 #convert the mass to SI units
	Rcalc_SI = Rcalc * 6.955E8 #convert the radius to SI units
	return Mcalc_SI, Rcalc_SI, Density_SI

def calculateM(Teff, Logg, Metalicity): #DONE
	if Teff == 0: #This means that Teff was not found in the header of the file
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
	if Teff == 0: #This means that Teff was not found in the header of the file
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

def main():
	foldername = raw_input('Please enter the name of the folder in which contains the data: ')
	if os.path.exists('{}/KeplerLDC.txt') == False:
		currentdirectory = os.getcwd() #Get the current directory
		dst = '{}/{}/'.format(currentdirectory,foldername) #the destination
		src = '{}/KeplerLDC.txt'.format(currentdirectory) #the destination for the kepler LDC file
		shutil.copy(src,dst)
		src_flats = src = '{}/flats.txt'.format(currentdirectory) #the destination for the flats document
		shutil.copy(src_flats, dst)
	os.chdir('{}'.format(foldername)) #changes the current directory to the specified directory
	transitcalcloop()
	return

main()