#! /usr/bin/env python
#A tool to manually perform the fits that GANDALF is supposed to do.
#Measure linear offset around OIII/Hbeta and subtract it
#Measure redshift - initial guess is imported
#Measure dispersion - initial guess is 7.5 (this is angstroms, might have to convert to velocity)
#Measure OIII/Hbeta flux
#Measure linear offset around NII/Halpha
#Measure NII/Halpha flux



from os.path import expanduser
import pyfits
import pylab as plt
import numpy as np
import sys
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import math
import os
import shutil
import gaussfitter
import random
from PyAstronomy import pyasl

plt.rcParams['legend.scatterpoints'] = 1
#Direct input 
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 18,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          'legend.fontsize': 14,
          }
plt.rcParams.update(params) 

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
def norm(x, mean, sd):
  norm = []
  for i in range(x.size):
    norm += [1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x[i] - mean)**2/(2*sd**2))]
  return np.array(norm)


def linear(x, *p):
	m, b = p
	return b + m*x
	
def func(x, a, b):
  	return a*x + b

HOME = expanduser("~")
galaxy = sys.argv[1]
mode = sys.argv[2]
sncut = sys.argv[3]
monte_iterations = int(sys.argv[4])
if len(sys.argv) > 5:
	color_excess = sys.argv[5]
SOF_DIR = HOME+'/Astro/reduced/'+galaxy+'pt1sof/'

#mode = "voronoi"
#mode = "one"


#Make sure the noisy values are not being saved to text files or being printed as plots.



data_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+galaxy+"all.fits"
if mode == "voronoi":
	bins_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/voronoi_2d_binning_output_emission.txt"
	voronoi_2d_bins_file = '/Users/jimmy/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/voronoi_2d_bins_emission.txt'
	gandalf_folder = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/gandalf/"
	folder_name = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/gandalf/manual_fits/"
	monte_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/monte_carlo_results.txt"
if mode == "one":
	bins_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/one_bin_output_emission.txt"
	voronoi_2d_bins_file = '/Users/jimmy/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/one_bin_bins_emission.txt'
	gandalf_folder = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/gandalf_one/"
	folder_name = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/gandalf_one/manual_fits/"
	monte_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/monte_carlo_results_one.txt"
voronoi_2d_bins_table=np.loadtxt(voronoi_2d_bins_file, dtype=float)
if os.path.isdir(folder_name):
	shutil.rmtree(folder_name)
if not os.path.isdir(gandalf_folder):
	os.mkdir(gandalf_folder) 
if not os.path.isdir(folder_name):
	os.mkdir(folder_name) 
bins_table = np.loadtxt(bins_file, dtype=float)
cube_hdu = pyfits.open(data_file)
img = cube_hdu[0].data
img_header = cube_hdu[0].header
cdelt = float(img_header["CDELT3"])
crval = float(img_header["CRVAL3"])
crpix = float(img_header["CRPIX3"])
wavelength_angstroms = crval + (np.arange(img.shape[0])*cdelt)
chosen_bin = 1
bins = bins_table[:,4]
x_pix = bins_table[:,2]
y_pix = bins_table[:,3]
if mode == "voronoi":
	xbin = voronoi_2d_bins_table[:,0]
	ybin = voronoi_2d_bins_table[:,1]
if mode == "one":
	xbin = voronoi_2d_bins_table[0]
	ybin = voronoi_2d_bins_table[1]
environmental_variables=np.loadtxt(HOME+"/Astro/reduced/"+galaxy+"pt1sof/"+galaxy+"_env.sh", delimiter="=", skiprows=2, dtype={'names': ('dummy', 'value'), 'formats': ('S10', 'f4')})
temp_variable = environmental_variables[1]
redshift = temp_variable[1]
print("Input redshift: "+str(redshift))


emitted_wavelengths = np.array([4861,4959, 5007])
wavelengths = emitted_wavelengths*(1.0+redshift)
print(wavelengths)
names = ['Hbeta 4861','OIII 4959','OIII 5007']
metallicities = np.zeros(int(max(bins+1)))
bpt_x = np.zeros(int(max(bins+1)))
bpt_y = np.zeros(int(max(bins+1)))
OIII_5007_flux_array = np.zeros(int(max(bins+1)))
h_beta_flux_array = np.zeros(int(max(bins+1)))
NII_6584_flux_array = np.zeros(int(max(bins+1)))
h_alpha_flux_array = np.zeros(int(max(bins+1)))
OIII_5007_wavelength_array = np.zeros(int(max(bins+1)))
OIII_5007_amplitude_array = np.zeros(int(max(bins+1)))
OIII_5007_dispersion_array = np.zeros(int(max(bins+1)))
h_beta_wavelength_array = np.zeros(int(max(bins+1)))
h_beta_amplitude_array = np.zeros(int(max(bins+1)))
h_beta_dispersion_array = np.zeros(int(max(bins+1)))
NII_6584_wavelength_array = np.zeros(int(max(bins+1)))
NII_6584_amplitude_array = np.zeros(int(max(bins+1)))
NII_6584_dispersion_array = np.zeros(int(max(bins+1)))
h_alpha_wavelength_array = np.zeros(int(max(bins+1)))
h_alpha_amplitude_array = np.zeros(int(max(bins+1)))
h_alpha_dispersion_array = np.zeros(int(max(bins+1)))

#noise_out = open('/Users/jimmy/Astro/noise.txt', "a")
#fits_out = open('/Users/jimmy/Astro/trouble3_fits.txt', "w")
#mask_out = open('/Users/jimmy/Astro/mask.txt','w')

if monte_iterations > 0:
	monte_out = open(monte_file, "a")
	num_lines = sum(1 for line in open(monte_file))
	if num_lines < 1:
		monte_out.write("#\tbinNum\txpix\typix\th_beta_flux\th_beta_wavelength\th_beta_amplitude\th_beta_dispersion\tOIII_5007_flux\tOIII_5007_wavelength\tOIII_5007_amplitude\tOIII_5007_dispersion\th_alpha_flux\th_alpha_wavelength\th_alpha_amplitude\th_alpha_dispersion\tNII_6584_flux\tNII_6584_wavelength\tNII_6584_amplitude\tNII_6584_dispersion\tLogOH_PP04_O3N2")
	print(num_lines)
	

#First fit like normal, then...
#For each spaxel:
#Measure the random noise
#Add the random noise to the whole spectrum
#Make a new measurement
#Append that new measurement to a text file with the bin number, x, and y coords, measured fluxes, amplitude, peak, dispersion, etc.

for iteration in np.arange(monte_iterations+1):
	print("Performing Monte Carlo iteration: "+str(iteration))
	if iteration == 0:
		print("Not adding random noise to this run.")
		print("Fitting to bin: "+str(chosen_bin))
	max_in_lines = 0.0
	for chosen_bin in range(int(max(bins)+1)):
	#for chosen_bin in range(5):
		fig = plt.figure(figsize=(15, 5))
		img_spectrum = img[:,0,0]*0
		for index, bin_num in enumerate(bins):
			if bin_num == chosen_bin:
				img_spectrum = img_spectrum+img[:,y_pix[index],x_pix[index]]
				spaxel_x = x_pix[index]+1
				spaxel_y = y_pix[index]+1
		spectrum = img_spectrum
		if max(spectrum) < 0.00001:
			continue
		#Correction for reddening
		#E_B_minus_V = 1.97*np.log10(balmer_dec/2.86)
		if len(sys.argv) > 5:
			RsubV = 3.1
			fluxUnred = pyasl.unred(wavelength_angstroms, spectrum, ebv=float(color_excess), R_V=RsubV)
			spectrum = fluxUnred
		noise_level = np.std(spectrum[310:390])*2
		#noise_out.write(str(noise_level)+'\n')
		if iteration != 0:
			#print("Measured noise level: "+str(noise_level))
			noise_array = [random.normalvariate(0.0, noise_level) for _ in xrange(len(spectrum))]
			spectrum = spectrum+noise_array

		
		fit_window_size = 500.0
		mask_window_size = 40.0
		first_cut = np.all([wavelength_angstroms>(wavelengths[1]-fit_window_size), wavelength_angstroms<(wavelengths[1]+fit_window_size)], axis=0)
		fit_cut = first_cut
		max_mask = np.all([wavelength_angstroms<3000, wavelength_angstroms>7000], axis=0)
		for line in wavelengths:
			line_mask = np.any([wavelength_angstroms<(line-mask_window_size), wavelength_angstroms>(line+mask_window_size)], axis=0)
			fit_cut = np.all([fit_cut,line_mask], axis=0)
			line_cut = np.all([wavelength_angstroms>(line-mask_window_size), wavelength_angstroms<(line+mask_window_size)], axis=0)
			max_mask = np.any([max_mask,line_cut], axis=0)
		try:
			lin_co, linear_fit = curve_fit(func, wavelength_angstroms[fit_cut], spectrum[fit_cut])
			linear_fit_result = func(wavelength_angstroms, *lin_co)
			#fits_out.write(str(chosen_bin)+'\t'+str(lin_co)+'\t')
		except:
			linear_fit_result = spectrum*0.0
		
		cont_subtracted_img = spectrum-linear_fit_result
		
		if iteration == -1:
			plot = fig.add_subplot(2,2,1)
			plt.plot(wavelength_angstroms[first_cut], spectrum[first_cut], color="blue", label="Data")
			try:
				plt.plot(wavelength_angstroms[first_cut], linear_fit_result[first_cut], color="red", label="Linear Continuum Fit")
			except:
				print(linear_fit_result)
			line_chooser = np.invert(max_mask)
			max_in_lines = max(spectrum[max_mask])
			min_in_lines = min(spectrum[max_mask])
			plt.ylim(-0.1, 0.4)
			if max(spectrum) > 10:
				plt.ylim(min_in_lines-max_in_lines*0.1, max_in_lines*1.3)
			plt.scatter(wavelength_angstroms[fit_cut], (wavelength_angstroms[fit_cut]*0.0)+(max(linear_fit_result)), color="lime", s=5, label="Continuum Fit Window")
			plt.legend(loc="upper left")
			#mask_out.write(str(wavelength_angstroms[fit_cut]))
			
			plot = fig.add_subplot(2,2,2)
			
			plt.plot(wavelength_angstroms[first_cut], cont_subtracted_img[first_cut], color="blue", label="Data")
			
		redshift_array = np.array([0.0, 0.0, 0.0])
		dispersion_array = np.array([0.0, 0.0, 0.0])
		intensity_array = np.array([0.0, 0.0, 0.0])
		peak_array = np.array([0.0, 0.0, 0.0])
		
		for index, line in enumerate(wavelengths):
			p0 = [0.25, line, 7.5]
			if max(spectrum) > 10:
				output = gaussfitter.multigaussfit(wavelength_angstroms, cont_subtracted_img, ngauss=1, err=None, params=[0.25, line, 7.5],
	    	    	fixed=[False,False,False], limitedmin=[True,False,True],
		    	    limitedmax=[False,False,True], minpars=[0,0,0], maxpars=[0,0,10.0],
					quiet=True, shh=True, veryverbose=False)
				coeff = output[0]
				gauss_line = gauss(wavelength_angstroms, *coeff)
			else:
				gauss_line = np.arange(len(spectrum))
				coeff=[0,line,7.5]
	
			redshift_array[index] = (coeff[1]-emitted_wavelengths[index])/emitted_wavelengths[index]
			dispersion_array[index] = coeff[2]
			intensity_array[index] = coeff[0]
			peak_array[index] = coeff[1]
			#print(redshift_array)
			#stop()
		
		
		redshift = np.median(redshift_array)
		if iteration == 0:
			print("Best fit redshift: "+str(redshift))	
		wavelengths = emitted_wavelengths*(1.0+redshift)
		median_dispersion = np.median(dispersion_array)
		if (median_dispersion < 6.5) or (median_dispersion > 8.5):
			median_dispersion = 7.5
			#print(" ")
			#print("Possible bad fiber: "+str(chosen_bin))
			#print(str(x_pix[chosen_bin])+" "+str(y_pix[chosen_bin]))
			#print(" ")
		dispersion_slop = (median_dispersion*0.1)
		wavelength_slop = 0.1
		
		residual = cont_subtracted_img
		
		colors = ['orange','purple','green']
		labels = ['H beta Fit','OIII 4959 Fit','OIII 5007 Fit']
		#print("peak_array[0]: "+str(peak_array[0]))
		#print("peak_array[1]: "+str(peak_array[1]))
		#print("peak_array[2]: "+str(peak_array[2]))
		#print("median_dispersion: "+str(median_dispersion))
		output = gaussfitter.multigaussfit(wavelength_angstroms, cont_subtracted_img, ngauss=3, err=None, params=[intensity_array[0], peak_array[0], median_dispersion, intensity_array[1], peak_array[1], median_dispersion, intensity_array[2], peak_array[2], median_dispersion],
		    	    fixed=[False,True,True], limitedmin=[True,True,True],
	    		    limitedmax=[False,True,True], minpars=[0,peak_array[0]-wavelength_slop,median_dispersion-dispersion_slop,0,peak_array[1]-wavelength_slop,median_dispersion-dispersion_slop,0,peak_array[2]-wavelength_slop,median_dispersion-dispersion_slop], maxpars=[0,peak_array[0]+wavelength_slop,median_dispersion+dispersion_slop,0,peak_array[1]+wavelength_slop,median_dispersion+dispersion_slop,0,peak_array[2]+wavelength_slop,median_dispersion+dispersion_slop], #minpars=[0,peak_array[0]*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop,0,peak_array[1]*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop,0,peak_array[2]*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop], ,
					quiet=True, shh=True, veryverbose=False)
		tripple_out = output[0]
		#fits_out.write(str(tripple_out)+'\t')
		h_beta_flux=tripple_out[0]*tripple_out[2]*math.sqrt(2*math.pi) #0 is height, 2 is width
		h_beta_wavelength = tripple_out[1]
		h_beta_amplitude = tripple_out[0]
		h_beta_dispersion = tripple_out[2]
		OIII_5007_flux=tripple_out[6]*tripple_out[8]*math.sqrt(2*math.pi) #0 is height, 2 is width
		OIII_5007_wavelength = tripple_out[7]
		OIII_5007_amplitude = tripple_out[6]
		OIII_5007_dispersion = tripple_out[8]
		if iteration == 0:
			print("H beta flux: "+str(h_beta_flux))
			print("OIII 5007 flux: "+str(OIII_5007_flux))
		if iteration == -1:
			hand_fit = gauss(wavelength_angstroms, *tripple_out[0:3])
			if max(hand_fit[first_cut] > 0):
				plt.plot(wavelength_angstroms[first_cut],hand_fit[first_cut], color='orange', label=r"H$\beta$ Fit")
			residual = residual-hand_fit
			hand_fit = gauss(wavelength_angstroms, *tripple_out[3:6])
			if max(hand_fit[first_cut] > 0):
				plt.plot(wavelength_angstroms[first_cut],hand_fit[first_cut], color='purple', label="OIII 4959 Fit")
			residual = residual-hand_fit
			hand_fit = gauss(wavelength_angstroms, *tripple_out[6:9])
			residual = residual-hand_fit
			plt.plot(wavelength_angstroms[first_cut], residual[first_cut], color="red", label="Residual")
			#plt.text(4450, ((min_in_lines-max(linear_fit_result)-max_in_lines*0.1)+(max_in_lines*1.3))/2, "Max residual: "+str(round(max(residual[max_mask]/noise_level),3)))#+" noise level: "+str(round(noise_level,3)))
			#plt.text(5200, max_in_lines/1.1, "in units of std dev")
			#plt.scatter(wavelength_angstroms[max_mask], (wavelength_angstroms[max_mask]*0.0)+(max_in_lines*0.15), color="blue", s=5, label="Lines Only")
			if max(hand_fit[first_cut] > 0):
				plt.plot(wavelength_angstroms[first_cut],hand_fit[first_cut], color='green', label="OIII 5007 Fit")
			plt.ylim(-0.1, 0.4)
			if max(spectrum) > 10:
				plt.ylim(min_in_lines-max(linear_fit_result)-max_in_lines*0.1, max_in_lines*1.3)
			plt.legend(loc="upper left")
		if iteration == 0:
			plot = fig.add_subplot(1,2,1)
			
		h_alpha_line = 6563*(1+redshift)
		
		params = np.genfromtxt(SOF_DIR+'circle.txt', delimiter='\t', names=True, dtype=None)
		
		#plot_cut = np.all([wavelength_angstroms>(h_alpha_line-1100), wavelength_angstroms<(h_alpha_line+500)], axis=0)
		first_cut = np.all([wavelength_angstroms>(5500), wavelength_angstroms<(6700)], axis=0)
		plot_cut = first_cut
		fit_cut = first_cut
		##mask = np.any([wavelength_angstroms<(6150), wavelength_angstroms>(6450)], axis=0)
		mask = np.any([wavelength_angstroms<(params['t1']), wavelength_angstroms>(params['t2'])], axis=0)
		fit_cut = np.all([fit_cut,mask], axis=0)
		#line_mask = np.any([wavelength_angstroms<(5860), wavelength_angstroms>(5930)], axis=0)
		line_mask = np.any([wavelength_angstroms<(5800), wavelength_angstroms>(5930)], axis=0)
		fit_cut = np.all([fit_cut,line_mask], axis=0)
		line_mask = np.any([wavelength_angstroms<(h_alpha_line-params['t3']), wavelength_angstroms>(h_alpha_line+params['t4'])], axis=0)
		##line_mask = np.any([wavelength_angstroms<(h_alpha_line-50), wavelength_angstroms>(h_alpha_line+50)], axis=0)
		fit_cut = np.all([fit_cut,line_mask], axis=0)
		line_mask = np.any([wavelength_angstroms<(5577-50), wavelength_angstroms>(5577+50)], axis=0)
		fit_cut = np.all([fit_cut,line_mask], axis=0)

		try:
			lin_co, linear_fit = curve_fit(func, wavelength_angstroms[fit_cut], spectrum[fit_cut])
			linear_fit_result = func(wavelength_angstroms, *lin_co)
			#fits_out.write(str(lin_co)+'\t')
		except:
			linear_fit_result = spectrum*0.0
		cont_subtracted_img = spectrum-linear_fit_result
		wavelength_slop = 2.0
		if iteration == 0:
			plt.plot(wavelength_angstroms[plot_cut], spectrum[plot_cut], color="blue", label="Data")
			try:
				plt.plot(wavelength_angstroms[plot_cut], linear_fit_result[plot_cut], color="red", label="Linear Continuum Fit")
			except:
				print(linear_fit_result)
			max_mask = np.all([wavelength_angstroms>(h_alpha_line-(mask_window_size/1)), wavelength_angstroms<(h_alpha_line+(mask_window_size/1))], axis=0)
			max_in_lines = max(spectrum[max_mask])
			h_alpha_assumed_amplitude = max(cont_subtracted_img[max_mask]) 
			print('Max in lines check: '+str(max_in_lines))
			min_in_lines = min(spectrum[max_mask])
			plt.ylim(-0.1, 0.4)
			if max(spectrum) > 10:
				plt.ylim(min_in_lines-max_in_lines*0.1, max_in_lines*1.3)
			plt.scatter(wavelength_angstroms[fit_cut], (wavelength_angstroms[fit_cut]*0.0)+(max(linear_fit_result)), color="lime", s=5, label="Continuum Fit Window")
			plt.scatter([h_alpha_line,h_alpha_line],[-10.0,-10.0], color='white', label='Spaxel: '+str(spaxel_x)+' '+str(spaxel_y))
			#plt.legend(loc="upper left")
			#mask_out.write(str(wavelength_angstroms[fit_cut])+'\n')
			plt.legend(loc="upper left")
			plt.xlabel(r'Wavelength \AA\ ')
			plt.ylabel(r'Flux [$10^{-16}$ erg cm$^{-2}$ s$^{-1}$]', labelpad=20)
			
			plot = fig.add_subplot(1,2,2)
			residual = cont_subtracted_img
			
			
		try:
			#if 1.0 > 0.0:
			#based on the amplitude in a window around H_alpha, and assuming the dispersion is 7.5 (or using the median_dispersion measurement) I can fix the flux ratio to be that of NII and H alpha for the integrated spectrum
			#per spaxels highest max is about 0.437941
			#integrated spaxels highest max is about 6.04156
			if h_alpha_assumed_amplitude < 1.0:
				#h_alpha_assumed_amplitude = max_in_lines
				#print('PASSED THE CHECK')
				h_alpha_assumed_flux = h_alpha_assumed_amplitude*median_dispersion*math.sqrt(2*math.pi) #0 is height, 2 is width
				#pull in h_alpha and NII ratio from one bin results and use those now to fix NII to the same assumed ratio
				one_bin_table = np.genfromtxt(HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/gandalf_table.txt',dtype=None)
				Ha_line = one_bin_table[6]
				stacked_Ha = Ha_line[1]
				NII_line = one_bin_table[5]
				stacked_NII = NII_line[1]
				assumed_ratio = stacked_Ha/stacked_NII
				NII_6584_assumed_flux = h_alpha_assumed_flux/assumed_ratio
				NII_6584_assumed_amplitude = NII_6584_assumed_flux/(median_dispersion*math.sqrt(2*math.pi))
				NII_6549_assumed_flux = NII_6584_assumed_flux/2.95
				NII_6549_assumed_amplitude = NII_6549_assumed_flux/(median_dispersion*math.sqrt(2*math.pi))
				
				#print('assumed_ratio: '+str(assumed_ratio))
				#print('h_alpha_assumed_amplitude: '+str(h_alpha_assumed_amplitude))
				#print('NII_6584_assumed_amplitude: '+str(NII_6584_assumed_amplitude))
				#print('h_alpha_assumed_flux: '+str(h_alpha_assumed_flux))
				#print('NII_6584_assumed_flux: '+str(NII_6584_assumed_flux))
				#print('NII_6549_assumed_flux: '+str(NII_6549_assumed_flux))
				#print(6563.0*(1.0+redshift))
				output = gaussfitter.multigaussfit(wavelength_angstroms, cont_subtracted_img, ngauss=3, err=None, params=[NII_6549_assumed_amplitude,6549.0*(1.0+redshift),median_dispersion,h_alpha_assumed_amplitude,6563.0*(1.0+redshift),median_dispersion,NII_6584_assumed_amplitude,6584.0*(1.0+redshift),median_dispersion],
	    			    fixed=[False,True,True], limitedmin=[True,False,False],
	    		    	limitedmax=[True,False,False], minpars=[NII_6549_assumed_amplitude*0.5,6549.0*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop,h_alpha_assumed_amplitude*0.90,6563.0*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop,NII_6584_assumed_amplitude*0.5,6584.0*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop], maxpars=[NII_6549_assumed_amplitude*1.0,6549.0*(1.0+redshift)+wavelength_slop,median_dispersion+dispersion_slop,h_alpha_assumed_amplitude*1.1,6563.0*(1.0+redshift)+wavelength_slop,median_dispersion+dispersion_slop,NII_6584_assumed_amplitude*1.0,6584.0*(1.0+redshift)+wavelength_slop,median_dispersion+dispersion_slop],
						quiet=True, shh=True, veryverbose=False)
				###output = gaussfitter.multigaussfit(wavelength_angstroms, cont_subtracted_img, ngauss=1, err=None, params=[h_alpha_assumed_amplitude,6563.0*(1.0+redshift),median_dispersion], fixed=[False,True,True], limitedmin=[True,False,False], limitedmax=[True,False,False], minpars=[h_alpha_assumed_amplitude*0.9,6563.0*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop], maxpars=[h_alpha_assumed_amplitude*1.5,6563.0*(1.0+redshift)+wavelength_slop,median_dispersion+dispersion_slop])
				tripple_out = output[0]
				###tripple_out = np.append(median_dispersion,tripple_out)
				###tripple_out = np.append(6549.0*(1.0+redshift),tripple_out)
				###tripple_out = np.append(0.0,tripple_out)
				###tripple_out = np.append(tripple_out,0.0)
				###tripple_out = np.append(tripple_out,6584.0*(1.0+redshift))
				###tripple_out = np.append(tripple_out,median_dispersion)
			else:
				output = gaussfitter.multigaussfit(wavelength_angstroms, cont_subtracted_img, ngauss=3, err=None, params=[0.5,6549.0*(1.0+redshift),median_dispersion,0.5,6563.0*(1.0+redshift),median_dispersion,0.5,6584.0*(1.0+redshift),median_dispersion],
	    			    fixed=[False,False,False], limitedmin=[True,True,True],
	    		    	limitedmax=[False,True,True], minpars=[0,6549.0*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop,0,6563.0*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop,0,6584.0*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop], maxpars=[0,6549.0*(1.0+redshift)+wavelength_slop,median_dispersion+dispersion_slop,0,6563.0*(1.0+redshift)+wavelength_slop,median_dispersion+dispersion_slop,0,6584.0*(1.0+redshift)+wavelength_slop,median_dispersion+dispersion_slop],
						quiet=True, shh=True, veryverbose=False)
				tripple_out = output[0]
			###output = gaussfitter.multigaussfit(wavelength_angstroms, cont_subtracted_img, ngauss=1, err=None, params=[0.5,6563.0*(1.0+redshift),median_dispersion], fixed=[False], limitedmin=[True], limitedmax=[False], minpars=[0,6563.0*(1.0+redshift)-wavelength_slop,median_dispersion-dispersion_slop], maxpars=[10000,6563.0*(1.0+redshift)+wavelength_slop,median_dispersion+dispersion_slop])
			#tripple_out = output[0]
			#print('output: '+str(tripple_out))
			#print(tripple_out)
			#fits_out.write(str(tripple_out)+'\n')
			fit = output[1]
		except:
			#else:
			tripple_out = [0.0,6549.0*(1.0+redshift),median_dispersion,0.0,6563.0*(1.0+redshift),median_dispersion,0.0,6584.0*(1.0+redshift),median_dispersion]
			fit = cont_subtracted_img*0.0
		
		if iteration == 0:
			#print(cont_subtracted_img[plot_cut])
			plt.plot(wavelength_angstroms[plot_cut], cont_subtracted_img[plot_cut], color="blue", label="Data")
			plt.plot(wavelength_angstroms[plot_cut], cont_subtracted_img[plot_cut]-fit[plot_cut], color="red", label="Residual")

			hand_fit = gauss(wavelength_angstroms, *tripple_out[0:3])
			if max(hand_fit[plot_cut] > 0):
				plt.plot(wavelength_angstroms[plot_cut],hand_fit[plot_cut], color='orange', label="NII 6549 Fit")
			residual = residual-hand_fit
			hand_fit = gauss(wavelength_angstroms, *tripple_out[3:6])
			if max(hand_fit[plot_cut] > 0):
				plt.plot(wavelength_angstroms[plot_cut],hand_fit[plot_cut], color='purple', label=r"H$\alpha$ Fit")
			max_alpha = max(hand_fit[plot_cut])
			residual = residual-hand_fit
			hand_fit = gauss(wavelength_angstroms, *tripple_out[6:9])
			residual = residual-hand_fit
			if max(hand_fit[plot_cut] > 0):
				plt.plot(wavelength_angstroms[plot_cut],hand_fit[plot_cut], color='green', label="NII 6584 Fit")
			#plt.scatter(wavelength_angstroms[max_mask], (wavelength_angstroms[max_mask]*0.0)+(max_in_lines*0.15), color="blue", s=5, label="Lines Only")
			#if max(hand_fit[plot_cut] > 0):
			#	plt.text(5500, ((min_in_lines-max(linear_fit_result)-max_in_lines*0.1)+(max_in_lines*1.3))/2, "Max residual: "+str(round(max(residual[max_mask]/noise_level),2)))#+" noise level: "+str(round(noise_level,2)))
			#plt.text(5500, max_in_lines/2.2, "in units of std dev")
			#print('maxhand_fit[plot_cut] check: '+str(max(hand_fit[plot_cut])))
			##if (max_alpha > 0):
				##plt.ylim(min_in_lines-max(linear_fit_result)-max_in_lines*0.1, max_in_lines*1.3)
			plt.ylim(-0.1, 0.4)
			if max(spectrum > 10.0):
				plt.ylim(min_in_lines-max(linear_fit_result)-max_in_lines*0.1, max_in_lines*1.3)
			plt.legend(loc="upper left")
			plt.xlabel(r'Wavelength \AA\ ')
			plt.subplots_adjust(bottom=0.14)
			plt.subplots_adjust(left=0.10)
			plt.subplots_adjust(right=0.98)
			plt.subplots_adjust(top=0.96)
		
			plt.savefig(folder_name+str(chosen_bin)+".pdf")
		
		############Try using the gandalf inst_resolution command...
		
		
		NII_6549_flux = tripple_out[0]*tripple_out[2]*math.sqrt(2*math.pi) #0 is height, 2 is width
		h_alpha_flux = tripple_out[3]*tripple_out[5]*math.sqrt(2*math.pi) #0 is height, 2 is width
		h_alpha_wavelength = tripple_out[4]
		h_alpha_amplitude = tripple_out[3]
		h_alpha_dispersion = tripple_out[5]
		NII_6584_flux = tripple_out[6]*tripple_out[8]*math.sqrt(2*math.pi) #0 is height, 2 is width
		NII_6584_wavelength = tripple_out[7]
		NII_6584_amplitude = tripple_out[6]
		NII_6584_dispersion = tripple_out[8]
		
		
		#print("NII 1 flux: "+str(NII_6549_flux))
		if iteration == 0:
			print("H alpha flux: "+str(h_alpha_flux))
			print("NII 2 flux: "+str(NII_6584_flux))
		
		#metallicity = log[(5007/hbeta)/(6584/halpha)]
		try:
			O3N2 = math.log10((OIII_5007_flux/h_beta_flux)/(NII_6584_flux/h_alpha_flux))
			LogOH_PP04_O3N2 = 8.73-(0.32*O3N2)
			N2 = math.log10(NII_6584_flux/h_alpha_flux)
			LogOH_PP04_N2 = 9.37+(2.03*N2)+(1.26*N2**2)+(0.32*N2**3)
			LogOH_D02 = 9.12+(0.73*N2)
			metallicities[chosen_bin] = LogOH_PP04_N2
			bpt_x[chosen_bin] = math.log10(NII_6584_flux/h_alpha_flux)
			bpt_y[chosen_bin] = math.log10(OIII_5007_flux/h_beta_flux)
		except:
			O3N2 = 0.0
			LogOH_PP04_O3N2 = 0.0
			LogOH_PP04_N2 = 0.0
			LogOH_D02 = 0.0
			metallicities[chosen_bin] = LogOH_PP04_O3N2
			bpt_x[chosen_bin] = 0.0
			bpt_y[chosen_bin] = 0.0
		
		OIII_5007_flux_array[chosen_bin] = OIII_5007_flux
		OIII_5007_wavelength_array[chosen_bin] = OIII_5007_wavelength
		OIII_5007_amplitude_array[chosen_bin] = OIII_5007_amplitude
		OIII_5007_dispersion_array[chosen_bin] = OIII_5007_dispersion
		h_beta_flux_array[chosen_bin] = h_beta_flux
		h_beta_wavelength_array[chosen_bin] = h_beta_wavelength
		h_beta_amplitude_array[chosen_bin] = h_beta_amplitude
		h_beta_dispersion_array[chosen_bin] = h_beta_dispersion
		NII_6584_flux_array[chosen_bin] = NII_6584_flux
		NII_6584_wavelength_array[chosen_bin] = NII_6584_wavelength
		NII_6584_amplitude_array[chosen_bin] = NII_6584_amplitude
		NII_6584_dispersion_array[chosen_bin] = NII_6584_dispersion
		h_alpha_flux_array[chosen_bin] = h_alpha_flux
		h_alpha_wavelength_array[chosen_bin] = h_alpha_wavelength
		h_alpha_amplitude_array[chosen_bin] = h_alpha_amplitude
		h_alpha_dispersion_array[chosen_bin] = h_alpha_dispersion
		
		
		if iteration == 0:
			print("metallicity: "+str(LogOH_PP04_O3N2))
			print("bpt x: "+str(bpt_x[chosen_bin]))
			print("bpt y: "+str(bpt_y[chosen_bin]))
		
		plt.close()
		
		if iteration != 0:
			monte_out.write('\n'+str(bins[chosen_bin])+'\t'+str(x_pix[chosen_bin])+'\t'+str(y_pix[chosen_bin])+'\t'+str(h_beta_flux)+'\t'+str(h_beta_wavelength)+'\t'+str(h_beta_amplitude)+'\t'+str(h_beta_dispersion)+'\t'+str(OIII_5007_flux)+'\t'+str(OIII_5007_wavelength)+'\t'+str(OIII_5007_amplitude)+'\t'+str(OIII_5007_dispersion)+'\t'+str(h_alpha_flux)+'\t'+str(h_alpha_wavelength)+'\t'+str(h_alpha_amplitude)+'\t'+str(h_alpha_dispersion)+'\t'+str(NII_6584_flux)+'\t'+str(NII_6584_wavelength)+'\t'+str(NII_6584_amplitude)+'\t'+str(NII_6584_dispersion)+'\t'+str(LogOH_PP04_O3N2))
			
#plt.clf()


	if iteration == 0:
		if mode == "voronoi":
			gandalf_file = '/Users/jimmy/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/gandalf.fits'
			pink_gandalf_file = '/Users/jimmy/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/pink_gandalf.fits'
		if mode == "one":
			gandalf_file = '/Users/jimmy/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/gandalf_one.fits'
			pink_gandalf_file = '/Users/jimmy/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/pink_gandalf_one.fits'
		if os.path.exists(gandalf_file):
			hdulist = pyfits.open(gandalf_file)
			hdulist[0].data = hdulist[0].data*0.0
			hdulist[1].data = hdulist[1].data*0.0
			hdulist[2].data = hdulist[2].data*0.0
			hdulist[3].data = hdulist[3].data*0.0
			hdulist[4].data = hdulist[4].data*0.0
			hdulist[5].data = hdulist[5].data*0.0
			hdulist[6].data = hdulist[6].data*0.0
			hdulist[7].data = hdulist[7].data*0.0
			#kinematic_measurements = hdulist[4].data
			kinematic_measurements = np.zeros([max(bins)+1,61]) #could auto generate the 61 value.
			if mode == "voronoi":
				kinematic_measurements[:,31] = OIII_5007_flux_array
				kinematic_measurements[:,33] = OIII_5007_wavelength_array
				kinematic_measurements[:,32] = OIII_5007_amplitude_array
				kinematic_measurements[:,34] = OIII_5007_dispersion_array
				kinematic_measurements[:,21] = h_beta_flux_array
				kinematic_measurements[:,23] = h_beta_wavelength_array
				kinematic_measurements[:,22] = h_beta_amplitude_array
				kinematic_measurements[:,24] = h_beta_dispersion_array
				kinematic_measurements[:,46] = NII_6584_flux_array
				kinematic_measurements[:,48] = NII_6584_wavelength_array
				kinematic_measurements[:,47] = NII_6584_amplitude_array
				kinematic_measurements[:,49] = NII_6584_dispersion_array
				kinematic_measurements[:,41] = h_alpha_flux_array
				kinematic_measurements[:,43] = h_alpha_wavelength_array
				kinematic_measurements[:,42] = h_alpha_amplitude_array
				kinematic_measurements[:,44] = h_alpha_dispersion_array
			
			hdulist[4].data = kinematic_measurements
			hdulist.writeto(pink_gandalf_file, clobber=True)
		if not os.path.exists(gandalf_file):
			kinematic_measurements = np.zeros([max(bins)+1,61]) #could auto generate the 61 value.
			if mode == "voronoi":
				kinematic_measurements[:,31] = OIII_5007_flux_array
				kinematic_measurements[:,33] = OIII_5007_wavelength_array
				kinematic_measurements[:,32] = OIII_5007_amplitude_array
				kinematic_measurements[:,34] = OIII_5007_dispersion_array
				kinematic_measurements[:,21] = h_beta_flux_array
				kinematic_measurements[:,23] = h_beta_wavelength_array
				kinematic_measurements[:,22] = h_beta_amplitude_array
				kinematic_measurements[:,24] = h_beta_dispersion_array
				kinematic_measurements[:,46] = NII_6584_flux_array
				kinematic_measurements[:,48] = NII_6584_wavelength_array
				kinematic_measurements[:,47] = NII_6584_amplitude_array
				kinematic_measurements[:,49] = NII_6584_dispersion_array
				kinematic_measurements[:,41] = h_alpha_flux_array
				kinematic_measurements[:,43] = h_alpha_wavelength_array
				kinematic_measurements[:,42] = h_alpha_amplitude_array
				kinematic_measurements[:,44] = h_alpha_dispersion_array
			pyfits.writeto(pink_gandalf_file, kinematic_measurements, img_header, clobber=True)
			pyfits.append(pink_gandalf_file, kinematic_measurements, img_header)
			pyfits.append(pink_gandalf_file, kinematic_measurements, img_header)
			pyfits.append(pink_gandalf_file, kinematic_measurements, img_header)
			pyfits.append(pink_gandalf_file, kinematic_measurements, img_header)
			pyfits.append(pink_gandalf_file, kinematic_measurements, img_header)
			pyfits.append(pink_gandalf_file, kinematic_measurements, img_header)
		if mode == "one":
			gandalf_table = open('/Users/jimmy/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/gandalf_table.txt', "w")
			gandalf_table.write("PP04_O3N2\t"+str(LogOH_PP04_O3N2)+"\n")
			gandalf_table.write("PP04_N2\t"+str(LogOH_PP04_N2)+"\n")
			gandalf_table.write("D02\t"+str(LogOH_D02)+"\n")
			gandalf_table.write("OIII_5007_flux\t"+str(OIII_5007_flux_array[0])+"\n")
			gandalf_table.write("h_beta_flux\t"+str(h_beta_flux_array[0])+"\n")
			gandalf_table.write("NII_6584_flux\t"+str(NII_6584_flux_array[0])+"\n")
			gandalf_table.write("h_alpha_flux\t"+str(h_alpha_flux_array[0])+"\n")
			gandalf_table.close()
			#np.savetxt(, np.column_stack([LogOH_PP04_O3N2, LogOH_PP04_O3N2]), fmt='%10.6f %10.6f')


