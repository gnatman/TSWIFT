#! /usr/bin/env python
#A tool to create a fake galaxy.

#Import a frame with no data in it
#Create gaussians of the right ratio for all the lines in the jimmy_lines file
#	What is the right ratio?  Maybe take AGC221000 SDSS results as a start
#Overlay those on a flat-ish field in a circle
#Save out to a file
#Bonus points: Vary S/N ratio with radius
#Need to make a gandalf.fits file to get the full output.

import pyfits
from os.path import expanduser
import numpy as np
import math
import pylab as plt
import sys

HOME = expanduser("~")

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


galaxy = 'AGC220261'
data_file = HOME+"/Astro/reduced/"+galaxy+"pro/temp.fits"
selector = float(sys.argv[1])
#Center the fake galaxy around x=28,y=26
#read in data file
cube_hdu = pyfits.open(data_file)
img = cube_hdu[0].data
segmentation_img = np.zeros_like(img[0,:,:])
var_field = cube_hdu[1].data
img_header = cube_hdu[0].header
cdelt = float(img_header["CDELT3"])
crval = float(img_header["CRVAL3"])
crpix = float(img_header["CRPIX3"])
wavelength_angstroms = crval + (np.arange(img.shape[0])*cdelt)
#define the spaxels which will hold the galaxy
x_gal = [18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38]
y_gal = [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]
x_center = np.median(x_gal)
y_center = np.median(y_gal)
#x_gal = [23]
#y_gal = [31]

#fig=plt.figure()

x_array = np.empty(0)
y_array = np.empty(0)
h_beta_flux_array = np.empty(0)
h_beta_wavelength_array = np.empty(0)
h_beta_amplitude_array = np.empty(0)
oiii_4959_flux_array = np.empty(0)
oiii_4959_wavelength_array = np.empty(0)
oiii_4959_amplitude_array = np.empty(0)
oiii_5007_flux_array = np.empty(0)
oiii_5007_wavelength_array = np.empty(0)
oiii_5007_amplitude_array = np.empty(0)
nii_6549_flux_array = np.empty(0)
nii_6549_wavelength_array = np.empty(0)
nii_6549_amplitude_array = np.empty(0)
h_alpha_flux_array = np.empty(0)
h_alpha_wavelength_array = np.empty(0)
h_alpha_amplitude_array = np.empty(0)
nii_6584_flux_array = np.empty(0)
nii_6584_wavelength_array = np.empty(0)
nii_6584_amplitude_array = np.empty(0)
dispersion_array = np.empty(0)

#start by taking the median value of all SDSS results
#table=np.loadtxt("/Users/jimmy/Astro/sdss_speclines.csv", delimiter=",", skiprows=1, dtype={'names': ('oiii_5007_flux', 'h_beta_flux', 'h_alpha_flux', 'nii_6584_flux', 'lgm_tot_p50', 'lgm_fib_p50'), 'formats': ('f4', 'f4', 'f4', 'f4', 'f4', 'f4')})
alfalfa_sdss_catalog_file = HOME+"/Astro/Catalogs/alfalfa_sdss_matched.txt"
alfalfa_sdss_catalog = np.genfromtxt(alfalfa_sdss_catalog_file, dtype=None, names=['AGCNr','dr7objid','oh_p50','lgm_tot_p50','sfr_tot_p50','h_beta_flux','oiii_5007_flux','h_alpha_flux','nii_6584_flux', 'PP04_Metallicity', 'hi_mass','z','petroR90_r','petroR50_r','petroMag_u','petroMag_g','petroMag_r','petroMag_i','petroMag_z'], delimiter=",", skiprows=1)
temp_h_beta_flux = np.array([])
temp_oiii_5007_flux = np.array([])
temp_h_alpha_flux = np.array([])
temp_nii_6584_flux = np.array([])
#desired_metallicity = 7.53
#for galaxy in alfalfa_sdss_catalog:
#	if (galaxy['oh_p50'] >= desired_metallicity) and (galaxy['oh_p50'] <= desired_metallicity+0.1):
#		temp_h_beta_flux = np.append(temp_h_beta_flux, galaxy['h_beta_flux'])
#		temp_h_alpha_flux = np.append(temp_h_alpha_flux, galaxy['h_alpha_flux'])
#		temp_oiii_5007_flux = np.append(temp_oiii_5007_flux, galaxy['oiii_5007_flux'])
#		temp_nii_6584_flux = np.append(temp_nii_6584_flux, galaxy['nii_6584_flux'])
#		print(galaxy['oh_p50'])
#h_beta = np.median(temp_h_beta_flux)
#h_alpha = np.median(temp_h_alpha_flux)
#nii_6584 = np.median(temp_nii_6584_flux)
#oiii_5007 = np.median(temp_oiii_5007_flux)
h_beta = 1.0
#line ratios taken form Berg et al. 2012
#if desired_metallicity == 7.2:
#	h_alpha = 2.81 #7.2
#	nii_6584 = 0.016 #7.2
#	oiii_5007 = 1.89 #7.2
#if desired_metallicity == 7.53:
#	h_alpha = 2.83 #7.53
#	nii_6584 = 0.03 #7.53
#	oiii_5007 = 2.89 #7.53
Berg_names_all = ['UGC 521 A', 'UGC695 E', 'UGC1056 A', 'UGC 1056 B', 'UGC 1176 A', 'NGC 784 A','NGC 784 B', 'UGC 2716 A', 'KKH 037 A', 'NGC 2537 A', 'NGC 2537 B', 'UGC 4278 B', 'UGC4278 A', 'NGC 2552 A', 'UGC 4393 B', 'UGC4393 C', 'CGCG 035-007 A', 'UGC 5139 A', 'IC 559 A', 'UGC 5272 A', 'UGC 5340 A', 'UGC 5423 A', 'UGC 5423 B', 'UGC 5672 A', 'UGC 5692 A', 'UGC 5797 A', 'UGC 5923 A', 'NGC 3741 A', 'NGC 3738 A', 'NGC 3738 B', 'UGC 6817 A', 'UGC 6900 A', 'NGC 4163 A', 'CGCG 269-049 C', 'CGCG 269-049 A', 'UGC 7577 A', 'NGC 4449 C', 'NGC 449 B', 'NGC 449 A', 'UGC 7605 A ', 'UGC 7639 A ', 'NGC 4656 A', 'UGC 8201 A', 'UGC 8245 A', 'UGC 8508 A', 'UGC 8638 A', 'UGC 8638 B', 'UGC 8837 A', 'NGC 5477 A', 'UGC 9405 A', 'UGC 10818 A', 'KKH 098 A']
Berg_OIII = np.array([3.64, 1.62, 2.35, 3.27, 3.61, 4.13, 3.32, 4.26, 0.53, 2.14, 1.78, 1.91, 2.53, 2.91, 3.28, 2.59, 2.03, 3.66, 2.80, 4.94, 1.89, 3.49, 3.71, 2.51, 1.70, 5.53, 2.48, 2.84, 2.96, 3.11, 2.89, 0.60, 0.49, 1.53, 2.51, 5.25, 2.33, 2.54, 3.46, 2.33, 1.33, 6.71, 2.94, 1.25, 3.25, 4.17, 4.15, 1.26, 4.64, 1.44, 2.22, 1.91])
Berg_Hb = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
Berg_Ha = np.array([2.77, 2.86, 2.80, 2.82, 2.82, 2.89, 2.81, 2.86, 2.80, 2.87, 2.79, 2.80, 2.86, 2.86, 2.86, 2.87, 2.84, 2.83, 2.86, 2.83, 2.81, 2.86, 2.86, 2.87, 2.79, 2.84, 2.78, 2.83, 2.83, 2.82, 2.83, 2.83, 2.75, 2.79, 2.83, 2.79, 2.84, 2.86, 2.83, 2.83, 2.83, 2.86, 2.81, 2.83, 2.83, 2.82, 2.81, 2.83, 2.84, 2.79, 2.79, 2.79])
Berg_NII = np.array([0.05, 0.125, 0.14, 0.12, 0.122, 0.109, 0.080, 0.092, 0.14, 0.48, 0.41, 0.060, 0.054, 0.25, 0.30, 0.31, 0.21, 0.07, 0.15, 0.035, 0.016, 0.120, 0.082, 0.23, 0.38, 0.09, 0.24, 0.051, 0.19, 0.20, 0.03, 0.52, 0.13, 0.0, 0.033, 0.09, 0.23, 0.205, 0.163, 0.066, 0.21, 0.023, 0.04, 0.12, 0.05, 0.07, 0.06, 0.202, 0.047, 0.40, 0.25, 0.08])
Berg_direct_metallicity = np.array([7.67, 7.69, 7.89, 7.98, 7.97, 8.03, 7.92, 7.97, 0.0, 8.24, 7.98, 7.69, 7.69, 8.15, 8.06, 7.99, 0.0, 7.92, 8.07, 7.87, 7.20, 7.77, 7.82, 0.0, 0.0, 7.96, 7.79, 7.68, 8.04, 8.03, 7.53, 0.0, 7.56, 0.0, 7.47, 7.97, 8.13, 8.16, 8.32, 7.66, 0.0, 8.09, 7.80, 0.0, 7.76, 7.97, 7.91, 7.82, 7.95, 0.0, 0.0, 0.0])
h_beta = Berg_Hb[selector]
h_alpha = Berg_Ha[selector]
nii_6584 = Berg_NII[selector]
oiii_5007 = Berg_OIII[selector]
desired_metallicity = Berg_direct_metallicity[selector]+(selector*0.000001)
N2 = np.log10(nii_6584/h_alpha)
O3N2 = np.log10((oiii_5007/h_beta)/(nii_6584/h_alpha))
text_file = open(HOME+"/Astro/reduced/AGC666pro/input_PP04_O3N2.txt", "w")
text_file.write(str(8.73-(0.32*O3N2)))
text_file.close()
text_file = open(HOME+"/Astro/reduced/AGC666pro/input_PP04_N2.txt", "w")
text_file.write(str(9.37 + (2.03*N2) + (1.26*N2**2) + (0.32*N2**3)))
text_file.close()
text_file = open(HOME+"/Astro/reduced/AGC666pro/input_D02_N2.txt", "w")
text_file.write(str(9.12 + 0.73*N2))
text_file.close()
text_file = open(HOME+"/Astro/reduced/AGC666pro/input_direct.txt", "w")
text_file.write(str(Berg_direct_metallicity[selector]))
text_file.close()
text_file = open(HOME+"/Astro/reduced/AGC666pro/h_alpha.txt", "w")
text_file.write(str(Berg_Ha[selector]))
text_file.close()
text_file = open(HOME+"/Astro/reduced/AGC666pro/h_beta.txt", "w")
text_file.write(str(Berg_Hb[selector]))
text_file.close()
#print(8.73-(0.32*O3N2))
#print(9.37 + (2.03*N2) + (1.26*N2**2) + (0.32*N2**3))
#print(9.12 + 0.73*N2)
fake_names = ['h_beta','oiii_4959','oiii_5007','nii_6549','h_alpha','nii_6584']
#h_beta = 6.53905219121 #9.11551681
#h_alpha = 25.9792229036
#nii_6584 = h_alpha*(10**((desired_metallicity-9.12)/0.73))
nii_6549 = nii_6584*0.33898305
#oiii_5007 =  ((nii_6584*h_beta)/h_alpha)*10**((8.73-desired_metallicity)/0.32)
oiii_4959 = oiii_5007*0.35
#fake_fluxes = [np.median(table['h_beta_flux']),np.median(table['oiii_5007_flux'])*0.35,np.median(table['oiii_5007_flux']),np.median(table['nii_6584_flux'])*0.33898305,np.median(table['h_alpha_flux']),np.median(table['nii_6584_flux'])]
#fake_fluxes = [301.682082853, 550.903491024*0.35, 550.903491024, 82.8573602808*0.33898305, 882.315032596, 82.8573602808] #NGC221000 integrated flux values, should probably divide by the number of spaxels...There are 82 bins
#fake_fluxes = [6.53905219121, 14.0100884928*0.35, 14.0100884928, 2.31209712408*0.33898305, 25.9792229036, 2.31209712408] #Values for the brightest spaxel in AGC221000
fake_fluxes = np.array([h_beta,oiii_4959,oiii_5007,nii_6549,h_alpha,nii_6584])*3.0

original_wavlengths = [4861, 4959, 5007, 6549, 6563, 6584]
fake_redshift = 0.004
fake_dispersion = 7.5
flux_threshhold = 5.0
flux_threshhold = 0.0
amplitude_threshhold = flux_threshhold/(fake_dispersion*math.sqrt(2*math.pi))

voronoi_2d_binning_emission = open(HOME+"/Astro/reduced/AGC666pro/voronoi_2d_binning_emission.txt", "w")
voronoi_2d_binning_output_emission = open(HOME+"/Astro/reduced/AGC666pro/voronoi_2d_binning_output_emission.txt", "w") 
voronoi_2d_bins_emission = open(HOME+"/Astro/reduced/AGC666pro/voronoi_2d_bins_emission.txt", "w")
one_bin_bins_emission = open(HOME+"/Astro/reduced/AGC666pro/one_bin_bins_emission.txt", "w") 
one_bin_output_emission = open(HOME+"/Astro/reduced/AGC666pro/one_bin_output_emission.txt", "w") 

voronoi_2d_binning_emission.write('#	x_arc, y_arc, x_pix, y_pix, signal, noise\n')
voronoi_2d_binning_output_emission.write('#          X"          Y"          Xpix          Ypix          BIN_NUM\n')
voronoi_2d_bins_emission.write('#          Xbin          Ybin          S/N     Num Pixels\n')
one_bin_bins_emission.write('#         Xbin          Ybin          S/N     Num Pixels		Total Noise\n')
x_med = np.median(x_gal)
y_med = np.median(y_gal)
counter=0
for x_coord in x_gal:
	for y_coord in y_gal:
		x_array = np.append(x_array, x_coord)
		y_array = np.append(y_array, y_coord)
		dispersion_array = np.append(dispersion_array, fake_dispersion)
		radius = math.sqrt(((x_coord-x_center)**2)+((y_coord-y_center)**2))
		for index,flux in enumerate(fake_fluxes):
			flux = flux/(math.e**(radius/2))
			wavelength = original_wavlengths[index]*(1+fake_redshift)
			amplitude = flux/(fake_dispersion*math.sqrt(2*math.pi))
			if amplitude > amplitude_threshhold/h_alpha:
				if index == 0:
					h_beta_flux_array = np.append(h_beta_flux_array, flux)
					h_beta_wavelength_array = np.append(h_beta_wavelength_array, wavelength)
					h_beta_amplitude_array = np.append(h_beta_amplitude_array, amplitude)
					params = [amplitude, wavelength, fake_dispersion]
					emission_line = gauss(wavelength_angstroms, *params)
					img[:,y_coord,x_coord] = img[:,y_coord,x_coord]+((emission_line))
			if index ==1:
				oiii_4959_flux_array = np.append(oiii_4959_flux_array, flux)
				oiii_4959_wavelength_array = np.append(oiii_4959_wavelength_array, wavelength)
				oiii_4959_amplitude_array = np.append(oiii_4959_amplitude_array, amplitude)
				params = [amplitude, wavelength, fake_dispersion]
				emission_line = gauss(wavelength_angstroms, *params)
				img[:,y_coord,x_coord] = img[:,y_coord,x_coord]+((emission_line))
			if index == 2:
				oiii_5007_flux_array = np.append(oiii_5007_flux_array, flux)
				oiii_5007_wavelength_array = np.append(oiii_5007_wavelength_array, wavelength)
				oiii_5007_amplitude_array = np.append(oiii_5007_amplitude_array, amplitude)
				params = [amplitude, wavelength, fake_dispersion]
				emission_line = gauss(wavelength_angstroms, *params)
				img[:,y_coord,x_coord] = img[:,y_coord,x_coord]+((emission_line))
			if index == 3:
				nii_6549_flux_array = np.append(nii_6549_flux_array, flux)
				nii_6549_wavelength_array = np.append(nii_6549_wavelength_array, wavelength)
				nii_6549_amplitude_array = np.append(nii_6549_amplitude_array, amplitude)
				params = [amplitude, wavelength, fake_dispersion]
				emission_line = gauss(wavelength_angstroms, *params)
				img[:,y_coord,x_coord] = img[:,y_coord,x_coord]+((emission_line))
			if amplitude > amplitude_threshhold:
				if index == 4:
					h_alpha_flux_array = np.append(h_alpha_flux_array, flux)
					h_alpha_wavelength_array = np.append(h_alpha_wavelength_array, wavelength)
					h_alpha_amplitude_array = np.append(h_alpha_amplitude_array, amplitude)
					params = [amplitude, wavelength, fake_dispersion]
					emission_line = gauss(wavelength_angstroms, *params)
					img[:,y_coord,x_coord] = img[:,y_coord,x_coord]+((emission_line))
					segmentation_img[y_coord,x_coord] = 1.
					voronoi_2d_binning_emission.write(str((x_coord-x_med-0.33)*0.66)+'\t'+str((y_coord-y_med)*0.66)+'\t'+str(x_coord)+'\t'+str(y_coord)+'\t5.0\t1.0\n')
					voronoi_2d_binning_output_emission.write(str((x_coord-x_med-0.33)*0.66)+'\t'+str((y_coord-y_med)*0.66)+'\t'+str(x_coord)+'\t'+str(y_coord)+'\t5.0\t'+str(counter)+'\n')
					counter=counter+1
					voronoi_2d_bins_emission.write(str((x_coord-x_med-0.33)*0.66)+'\t'+str((y_coord-y_med)*0.66)+'\t5.0\t1\n')
					one_bin_output_emission.write(str((x_coord-x_med-0.33)*0.66)+'\t'+str((y_coord-y_med)*0.66)+'\t'+str(x_coord)+'\t'+str(y_coord)+'\t0.0\n')
			if index == 5:
				nii_6584_flux_array = np.append(nii_6584_flux_array, flux)
				nii_6584_wavelength_array = np.append(nii_6584_wavelength_array, wavelength)
				nii_6584_amplitude_array = np.append(nii_6584_amplitude_array, amplitude)
				params = [amplitude, wavelength, fake_dispersion]
				emission_line = gauss(wavelength_angstroms, *params)
				img[:,y_coord,x_coord] = img[:,y_coord,x_coord]+((emission_line))
			#if amplitude > amplitude_threshhold:
			#params = [amplitude, wavelength, fake_dispersion]
			#emission_line = gauss(wavelength_angstroms, *params)
			#img[:,y_coord,x_coord] = img[:,y_coord,x_coord]+((emission_line))

one_bin_bins_emission.write('0\t0\t5'+str(counter)+'\t')
voronoi_2d_binning_emission.close()
voronoi_2d_binning_output_emission.close()
voronoi_2d_bins_emission.close()
one_bin_bins_emission.close()
one_bin_output_emission.close()

text_file = open(HOME+"/Astro/reduced/AGC666pro/h_alpha2.txt", "w")
text_file.write(str(np.sum(h_alpha_flux_array)))
text_file.close()
text_file = open(HOME+"/Astro/reduced/AGC666pro/h_beta2.txt", "w")
text_file.write(str(np.sum(h_beta_flux_array)))
text_file.close()

#cube_hdu[0].data = img
output_filename = HOME+"/Astro/reduced/AGC666pro/temp"+str(desired_metallicity)+".fits"
pyfits.writeto(output_filename, img, img_header, clobber=True)
pyfits.append(output_filename, var_field, img_header)
segmentation_filename = HOME+"/Astro/reduced/AGC666pro/segmentation.fits"
pyfits.writeto(segmentation_filename, segmentation_img, clobber=True)
#np.savetxt(HOME+"/Astro/reduced/AGC666pro/input_fake"+str(desired_metallicity)+".txt", np.column_stack([x_array,y_array,oiii_5007_flux_array,oiii_5007_wavelength_array,oiii_5007_amplitude_array,h_beta_flux_array,h_beta_wavelength_array,h_beta_amplitude_array,h_alpha_flux_array,h_alpha_wavelength_array,h_alpha_amplitude_array,nii_6584_flux_array,nii_6584_wavelength_array,nii_6584_amplitude_array,dispersion_array]), fmt='%8i %8i %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f', header="	x	y	oiii_5007_flux	oiii_5007_wavelength	oiii_5007_amplitude	h_beta_flux	h_beta_wavelength	h_beta_amplitude	nii_6584_flux	nii_6584_wavelength	nii_6584_amplitude	h_alpha_flux	h_alpha_wavelength	h_alpha_amplitude	dispersion")
text_file = open(HOME+"/Astro/reduced/AGC666pro/fake_metallicity.sh", "w")
text_file.write('export fake_metallicity='+str(desired_metallicity))
text_file.close()