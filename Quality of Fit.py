#Part of the Bound (Realistic Stars) data pipeline

#Simple script which recalculates (using the fitted parameters) the analytical lightcurves in such a way that a quality of fit can be found,
#which will then be used to tell if the transit is a false positive

import shutil
import numpy as np
import os
import matplotlib
matplotlib.use('PDF')
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import math
import transit
import warnings
warnings.simplefilter('ignore', RuntimeWarning) #The analytical lightcurve generator makes runtime errors but still works - so might as well suppress the errors
np.set_printoptions(threshold='nan')

def alccompare(foldernameold, foldernamenew, foldernamesave):
	STARdata = np.load('{}/Stellardataarray(merged).npy'.format(foldernamenew)) #Load the bound, fit, stellar data
	oldSTARdata = np.load('{}/Stellardataarray(fitted).npy'.format(foldernameold)) #Load the unbound, fit, stellar data
	newPLANdata = np.load('{}/plandataarray(merged).npy'.format(foldernamenew)) #load the bound, fit, planetary data
	oldPLANdata = np.load('{}/plandataarray(fitted).npy'.format(foldernameold)) #load the unbound, fit, planetary data
	for star in STARdata: #loop of all the stars/planets in the bound fit folder
		Teff_fit_bound, Logg_fit_bound, Metalicity_fit_bound, Teff_fit_unbound, Logg_fit_unbound, Metalicity_fit_unbound, Teff_kep, Logg_kep, Metalicity_kep, \
			Inclination_rad_bound, Planetaryradius_SI_bound, Tdur_SI_bound, Tdur_days_unbound, R_plan_over_R_star_unbound, Impactpar_norm_unbound,\
			Orbtialperiod_days = dataextractor(star, oldSTARdata, newPLANdata, oldPLANdata) #Extract Parameters for both the bound and unbound transits
		gamma1_bound, gamma2_bound = gammainterpolator(Teff_fit_bound, Logg_fit_bound, np.log10(Metalicity_fit_bound)) #Calculate limb darkening coefficients for the bound fit
		gamma1_unbound, gamma2_unbound = gammainterpolator(Teff_fit_unbound, Logg_fit_unbound, np.log10(Metalicity_fit_unbound)) #calculate limb darkening coefficients for the unbound fit
		Mcalc_SI_bound, Rcalc_SI_bound, Density_SI_bound = RMrhocalc(Teff_fit_bound, Logg_fit_bound, Metalicity_fit_bound) #Calculate M, R and rho for the bound fit
		unboundalc = np.load('{}/kepler-{}-planet-{}-analyticaltransitlc.npy'.format(foldernameold, star['KeplerID'], int(star['planetnumber']))) #load unbound analytical lightcurve
		unboundmin = np.min(unboundalc[1]) #find the minimum
		unbound_tdepth = 1 - unboundmin #and then find the transit depth
		boundalc = np.load('{}/kepler-{}-planet-{}-analyticaltransitlc.npy'.format(foldernamenew, star['KeplerID'], int(star['planetnumber']))) #load the bound analytical lightcurve
		boundmin = np.min(boundalc[1]) #find the minimum
		bound_tdepth = 1 - boundmin #and then find the transit depth
		if boundmin == 1.0 or unboundmin == 1.0: #this means that the lighturve is flat - i.e. the fit has failed in one, or both, of these cases
			flat_indicator = True
		else:
			flat_indicator = False
		print 'Lightcurve is flat: {}'.format(flat_indicator) #not really that useful
		#now, since I forgot to save the chi^2 values for the unbound fit, I need to calculate them: - not doing this at the moment
		keplerlc = np.load('{}/kepler-{}-planet-{}-trandata.npy'.format(foldernamenew, star['KeplerID'], int(star['planetnumber']))) #load the actual kepler lightcurve - used for the comparative transit plots
		keplerlc_unbound = np.load('{}/kepler-{}-planet-{}-trandata.npy'.format(foldernameold, star['KeplerID'], int(star['planetnumber']))) #load the, unbound, actual kepler lightcurve - only used to calculate the unbound chi^2 value
		#Now need to find out which fit gives the better chi2reduced value
		if (Tdur_SI_bound/86400) <= Tdur_days_unbound or (Tdur_SI_bound/86400) > (3*Tdur_days_unbound): #Need to use the longer transit duration, unless the bound transit is flat, in  which case we should use the other transit duration - once this is done, we should make a list containing the relevent data
			Trandata = [Tdur_days_unbound, unbound_tdepth]
		else:
			Trandata = [(Tdur_SI_bound / 86400), bound_tdepth] #Transit duration should be in days
		Time_datapoints = np.linspace((-Trandata[0]/2.0), (Trandata[0]/2.0), 200) #make an array of 200 data points which fall within the best fit transit array
		f_unbound = unbound_transit_calc(Impactpar_norm_unbound, Tdur_days_unbound, R_plan_over_R_star_unbound, gamma1_unbound, gamma2_unbound, Time_datapoints) #calculate the unbound transit
		f_bound = bound_transit_calc(Inclination_rad_bound, Planetaryradius_SI_bound, Orbtialperiod_days, Rcalc_SI_bound, Density_SI_bound, gamma1_bound, gamma2_bound, Time_datapoints) #calculate the bound transit
		comtransitplot(Time_datapoints, f_bound, f_unbound, keplerlc, star['KeplerID'], star['planetnumber'], foldernamesave) #Currently used to confirm that everything is working 
		#now need to normalise the transits - first normalise the transit durations
		Time_datapoints_norm = (1/Trandata[0]) * Time_datapoints
		# Now need to normalise the transit depth
		f_unbound_norm = (f_unbound - unboundmin) / unbound_tdepth
		f_bound_norm = (f_bound - boundmin) / bound_tdepth
		#Now to remove any nan values - set them to one since they only tend to occur at the very edges of transits
		for entry in [f_unbound, f_bound, f_unbound_norm, f_bound_norm]:
			for i in range(0, len(entry), 1):
				if np.isnan(entry[i]) == True:
					entry[i] = 1
		comtransitplot2(Time_datapoints_norm, f_bound_norm, f_unbound_norm, star['KeplerID'], star['planetnumber'], foldernamesave) #Currently used to confirm that everything is working - plot 2
		#Now to calculate the quality of the fit
		QoF = qualityoffit(Time_datapoints_norm, f_bound_norm, f_unbound_norm)
		#now going to invert the quality of fit, so that a high quality of fit indicates a high quality fit
		if QoF == 0: #if the quality of fit is 0, manually change the value to infinity inorder to prevent an error
			QoF == np.inf
		else:
			QoF = 1/QoF
		print 'The Quality of fit is: {}'.format(QoF) #Testing 
		#Now going to calculate the chi^2 values:
		Keplertransit = np.sort(keplerlc, order = 'Foldtime') #used for the bound fit
		Keplertransit_unbound = np.sort(keplerlc_unbound, order = 'Foldtime') #used for the unbound fit
		sigma, sigma_count = sigmacalc(Keplertransit, Trandata[0]) #calculate the out of transit standard deviation
		print 'Unbound'
		chi2_unbound, chi2reduced_unbound, count_unbound = chicalc(Keplertransit_unbound, unboundalc, sigma, Trandata[0]) #No longer in use
		# - bound fit:
		print 'Bound'
		chi2_bound, chi2reduced_bound, count_bound = chicalc(Keplertransit, boundalc, sigma, Trandata[0]) #No longer in use
		#now need to find out with analytical lightcurve is a better fit and then use that lightcurves data to calculate the Signal-to-Noise Ratio
		if chi2reduced_bound < chi2reduced_unbound:
			boundbetter = True
			depth_StN = bound_tdepth
			No_datapoints = count_bound
		else:
			boundbetter = False
			depth_StN = unbound_tdepth
			No_datapoints = count_unbound
		#Now going to calculate the Signal-to-Noise ratio
		StN_ratio = (depth_StN * (math.sqrt(No_datapoints)) / (sigma * math.sqrt(sigma_count)))
		print 'The Signal-to-Noise ratio is: {}'.format(StN_ratio)
		dtypes = [('KeplerID', int), ('planetnumber', int), ('QualityofFit', float), ('sigmazero', float), ('Signal-to-Noise', float), ('BoundBetter', bool)]
		values = [(star['KeplerID'], star['planetnumber'], QoF, sigma, StN_ratio, boundbetter)]
		dataarray = np.array(values, dtype=dtypes) #create a temporary array for this planets data
		if os.path.exists('{}/QoFdataarray.npy'.format(foldernamesave)) == False: #If the saved QoF data array does not yet exist, create an array to contain this data
			outdataarray = np.copy(dataarray)
		else: #otherwise load the saved QoF array and append the new QoF data to it
			outdataarray = np.load('{}/QoFdataarray.npy'.format(foldernamesave))
			outdataarray = np.concatenate((outdataarray, dataarray), axis=0)
		sortedarray = np.sort(outdataarray, order = 'QualityofFit') #Sort the array by the QoF value
		np.save('{}/QoFdataarray'.format(foldernamesave), sortedarray) #save the sorted structured array
		print 'Kepler{} complete:'.format(star['KeplerID'])
	Qof_StN_basic(sortedarray, foldernamesave)
	return

def dataextractor(star, oldSTARdata, newplandataarray, oldplandataarray):
	#New is bound, old is unbound!
	for line in oldSTARdata: #loop over the unbound stellar data
		if star['KeplerID'] == line['KeplerID'] and star['planetnumber'] == line['planetnumber']: #until the matching entry is found
			#And then extract the required parameters
			Teff_fit_new = star['Teff_fit']
			Logg_fit_new = star['Logg_fit']
			Metalicity_fit_new = star['Metalicity_fit']
			Teff_kep = star['Teff_kep']
			Logg_kep = star['Logg_kep']
			Metalicity_kep = star['Metalicity_kep']
			Teff_fit_old = line['Teff_fit']
			Logg_fit_old = line['Logg_fit']
			Metalicity_fit_old = line['Metalicity_fit']
			break
	for entry in newplandataarray: #loop over the bound planetary data array
		if star['KeplerID'] == entry['KeplerID'] and star['planetnumber'] == entry['planetnumber']: #until the matching entry is found
			Planetaryradius_SI_new = entry['Rplan(SI)']
			Inclination_new = entry['Inclination']
			Tdur_SI_new = entry['Transitduration']
			Orbtialperiod = entry['Orbitalperiod']
			break
	for planet in oldplandataarray: #loop ober the unbound planetary data array
		if star['KeplerID'] == planet['KeplerID'] and star['planetnumber'] == planet['planetnumber']: #until the matching entry is found
			Tdur_days_old = planet['Transitduration']
			R_plan_over_R_star_old = planet['R_plan/R_star']
			Impactpar_old = planet['Impactpar']
			break
	return 	Teff_fit_new, Logg_fit_new, Metalicity_fit_new, Teff_fit_old, Logg_fit_old, Metalicity_fit_old, Teff_kep, Logg_kep, Metalicity_kep, Inclination_new, Planetaryradius_SI_new, Tdur_SI_new, Tdur_days_old, R_plan_over_R_star_old, Impactpar_old, Orbtialperiod

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
		returnarray = [point[0], point[1]] #create an array to store the data points
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

def unbound_transit_calc(Impactpar, Transitduration_days, R_plan_over_R_star, gamma1, gamma2, Time_datapoints):
	zvalues = zandtcalculator_unbound(Time_datapoints, Impactpar, Transitduration_days, R_plan_over_R_star) #calculate the zvalues
	f = transit.occultquad(zvalues, R_plan_over_R_star, [gamma1, gamma2]) #and then use these to calculate the analytical transit
	return f

def zandtcalculator_unbound(tvalues, impactparameter, transitduration, R_p_over_R_s):
	#first calculate the minimum z value - turns out impact parameter is already a fraction of stellar radius
	zmin = impactparameter
	#print zmin #testing
	#and then calculate the max z value, which occurs as transitduration / 2
	zmax = 1 + R_p_over_R_s
	tmax = transitduration / 2
	#now calculate the relationship between t and z
	z_t = (zmax - zmin) / tmax
	#now calculate the corresponding z values
	zvalues = np.copy(tvalues) #make a copy of the tvalues array to store the z values in
	for i in range(0, len(tvalues), 1): #now calculate the zvalue associated with each t, remembering that z cannot be negative
		zvalues[i] = np.abs(tvalues[i] * z_t) + zmin
	return zvalues

def bound_transit_calc(Inclination_rad, Planetaryradius_SI, Orbtialperiod_days, Rstar_SI, RHOstar_SI, gamma1, gamma2, Time_datapoints):
	Orbitalperiod_SI = Orbtialperiod_days * 86400 #convert the orbital period to SI units
	zvalues = zandtcalculator_bound(Time_datapoints, RHOstar_SI, Inclination_rad, Orbitalperiod_SI) #calculate the zvalues
	Planetaryradiusratio = Planetaryradius_SI / Rstar_SI #calculate the ratio of R_plan/R_star
	f = transit.occultquad(zvalues, Planetaryradiusratio, [gamma1, gamma2]) #and then calculate the analytical transit
	return f

def zandtcalculator_bound(tvalues, Stellardensity_SI, inclination, Orbitalperiod_SI): #DONE - add comments
	#Orbitalperiod_SI = Orbitalperiod_SI / 2 #FOR SOME REASON THIS FIXES A LOAD OF PROBLEMS
	Gravconst =  6.67300E-11
	tvalues2 = tvalues * 86400 #convert the times from days to seconds
	phasevalues = tvalues2 / (Orbitalperiod_SI) #convert the times to orbital phase
	zvalues = np.copy(tvalues) #create an array to store the zvalues
	M_over_R_cubed = Stellardensity_SI  * (4 * math.pi / 3) #note, this is not actually the stellar density
	for i in range(0, len(phasevalues), 1): #now calculate the zvalue associated with each t, using the formulae David gave me (i.e. the one i was originally going to use)
		zvalues[i] = (((Gravconst * M_over_R_cubed * (Orbitalperiod_SI**2)) / (4*math.pi*math.pi))**(1.0/3)) * math.sqrt(((math.sin(phasevalues[i] * math.pi * 2))**2) + ((math.cos(inclination)*math.cos(phasevalues[i] * math.pi * 2))**2))
	#print (((Gravconst * M_over_R_cubed * (Orbitalperiod_SI**2)) / (4*math.pi*math.pi))**(1.0/3)) * math.cos(inclination)
	return zvalues

def qualityoffit(Time_datapoints_norm, f_bound_norm, f_unbound_norm): #used to calculate the Quality of Fit value, which is just the sum of the square of the difference between the two analytical lightcurves
	for i in range(0, len(Time_datapoints_norm), 1):
		if not 'F' in locals():
			F = (f_bound_norm[i] - f_unbound_norm[i])**2
		else:
			F += (f_bound_norm[i] - f_unbound_norm[i])**2
	return F

def sigmacalc(Keplerlc, Tdur): #calculates the deviation of all data points - essentially an error value which can be used to calculate chi^2 - No longer in use
	count = 0
	for i in range(0,len(Keplerlc), 1): #loops over all the data points inorder to create an array which only contains the near transit data
		if np.abs(Keplerlc['Foldtime'][i]) >= Tdur:
			count +=1
			if not 'sigmacalcarray' in locals():
				sigmacalcarray = Keplerlc['Flux'][i]
			else:
				sigmacalcarray = np.append(sigmacalcarray, Keplerlc['Flux'][i])
	if not 'sigmacalcarray' in locals():
		print 'Recursive Sigma'
		sigma, count = sigmacalc(Keplerlc, (Tdur/1.5))
	else:
		sigma = np.std(sigmacalcarray) #calculates the standard deviation of the near transit points 
	#print 'Sd is {}'.format(sigma) #testing 
	return sigma, count

def chicalc(Keplerlc, analyticallc, sigma, T_dur): #calculates the chi^2 and reduced chi^2 values whilst only using those data points that fall within a transit - No longer in use
	print len(Keplerlc['Foldtime'])
	print len(analyticallc[1])
	count = 0
	for i in range(0, len(analyticallc[1]), 1): #loops over the kepler lc (and analytical lc) inorder to find each points chi^2 value
		if np.abs(Keplerlc['Foldtime'][i]) <= T_dur or 1==1: #if the data point is part of the required group (i.e. falls within the transit)
			count += 1
			if not 'chi2' in locals(): #for the first point, create a chi^2 variable
				chi2 = (((Keplerlc['Flux'][i] - analyticallc[1][i])**2) / sigma**2)
			else: #otherwise add the next points chi^2 value to the total chi^2 value
				chi2 = chi2 + (((Keplerlc['Flux'][i] - analyticallc[1][i])**2) / sigma**2)
		else:
			continue
	if count == 0:
		return np.nan, np.nan, 0
	print 'The chi^2 value is {}'.format(chi2) #Testing 
	chi2reduced = chi2 / (count + 3) #now calculate the reduced chi^2 value for 4 independent variables
	print 'The reduced chi^2 value is {}'.format(chi2reduced) #testing
	print 'The number of data points is {}'.format(count)
	return chi2, chi2reduced, count #return both chi^2 values 

def main():
	foldernameold = raw_input('Please enter the name of the folder which contains the unbound analytical lightcurves: ') 
	foldernamenew = raw_input('Please enter the name of the folder which contains the bound analytical lightcurves: ')
	foldernamesave = raw_input('Please enter the name of the folder in which the output should be saved: ')
	if not os.path.isdir("./" + foldernamesave + "/"): #check to see if the given folder exists
		os.mkdir("./" + foldernamesave + "/") #if not, it creates the folder
	alccompare(foldernameold, foldernamenew, foldernamesave)
	currentdirectory = os.getcwd() #Get the current directory
	print 'Note: Actual processing is now complete. Now moving files around'
	dst = '{}/{}/'.format(currentdirectory,foldernamesave) #the destination for the moved data
	src_STARdata= '{}/{}/Stellardataarray(merged).npy'.format(currentdirectory, foldernamenew)
	src_PLANdata = '{}/{}/plandataarray(merged).npy'.format(currentdirectory, foldernamenew)
	shutil.copy(src_STARdata, dst) #Copy the stellar data array
	shutil.copy(src_PLANdata, dst) #Copy the planetary data array
	decide =  int(raw_input('Do you want to move the impact parameter data array from the merged data folder (YES=1, NO=0): '))
	if decide == 1:
		src_impact = '{}/{}/impactpardiscardarray.npy'.format(currentdirectory, foldernamenew)
		shutil.copy(src_impact, dst)
	return

def comtransitplot(Time_datapoints, f_bound, f_unbound, keplerlc, KeplerID, planetnumber, foldernamesave): #Produces a plain plot of the two analytical lightcurves and the original kepler data
	std = np.std(keplerlc['Flux'])
	plt.clf() #clear all figures
	plt.figure() #create a figure for this plot
	f1 = plt.plot( keplerlc['Foldtime'], keplerlc['Flux'], '+', label = ['Kepler'], color = 'grey') #plot both sets of data
	f2 = plt.plot( Time_datapoints, f_unbound, '-', label = ['Analytical(free)'], color = 'red') #plot both sets of data
	f3 = plt.plot( Time_datapoints, f_bound, '-', label = ['Analytical(bound)'], color = 'green') #plot both sets of data
	leg = plt.legend([f1,f2,f3], ['Kepler', 'Analytical(free)', 'Analytical(bound)'], fancybox=True, loc='lower right' ) #adds a legend to the plot
	leg.get_frame().set_alpha(0.5) #set the legand to be slightly transparent
	plt.xlabel('Time From Mid-Transit (Days)')
	plt.ylabel('Normalised Flux')
	plt.ylim([(1-6*std),(1+6*std)]) #Y limits based upon the standard deviation of the data points (prevents the plotting of extreme outliers)
	plt.title('Kepler and Analytical Transit Lightcurves')
	plt.savefig('{}/kepler-{}-planet-{}-intran-Analyticalcomptranlc.pdf'.format(foldernamesave, KeplerID, planetnumber)) #save the plot
	return

def comtransitplot2(Time_datapoints, f_bound, f_unbound, KeplerID, planetnumber, foldernamesave): #Produces a plain plot of the two analytical lightcurves and the original kepler data
	plt.clf() #clear all figures
	plt.figure() #create a figure for this plot
	f1 = plt.plot( Time_datapoints, f_unbound, '-', label = ['Analytical(free)'], color = 'red') #plot both sets of data
	f2 = plt.plot( Time_datapoints, f_bound, '-', label = ['Analytical(bound)'], color = 'green') #plot both sets of data
	leg = plt.legend([f1,f2], ['Analytical(free)', 'Analytical(bound)'], fancybox=True, loc='lower right' ) #adds a legend to the plot
	leg.get_frame().set_alpha(0.5) #set the legend to be slightly transparent
	plt.xlabel('Time From Mid-Transit (Days)')
	plt.ylabel('Normalised Flux')
	plt.title('Kepler and Analytical Transit Lightcurves')
	plt.savefig('{}/kepler-{}-planet-{}-intran-Analytical-only-comptranlc.pdf'.format(foldernamesave, KeplerID, planetnumber)) #save the plot
	return

def Qof_StN_basic(QoF_array, foldernamesave):
	plt.clf()
	plt.figure()
	plt.plot(QoF_array['Signal-to-Noise'], QoF_array['QualityofFit'], '.')
	plt.xlabel('Signal-to-Noise Ratio')
	plt.ylabel('Quality of Fit')
	plt.title('StN Ratio vs. QoF')
	plt.savefig('{}/StN_vs_QoF.pdf'.format(foldernamesave))
	return

main()