#! /usr/bin/env python
#A tool to create any plot I will need for my 2nd paper.  One script to rule them all.

import numpy as np
from os.path import expanduser
import pyfits
import pylab as plt
from sauron_colormap import sauron
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import scipy.optimize as optimization
from scipy.odr import *
import math
import scipy.stats as stats
import random
import sys
plt.rcParams['legend.scatterpoints'] = 1
#Direct input 
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

dog = 'n' #to convert certain values into dog years

def linear_fit(x, a, b):
    return a + b*x
def second_order_fit(x,a,b,c):
	return a + b*x + c*x**2
def quad_func(p, x):
     a, b, c = p
     return a + b*x + c*x**2
def fourth_order_fit(x,a,b,c,d,e):
	return a + b*x + c*x**2 + d*x**3 + e*x**4
#def second_order_2d(m, s, c, ma, mb, sa, sb)
#	return c+ma*m+mb*m**2+sa*s+sb*s**2
def odr_moustakas(p, x):
	mto, asm, gamma = p
	return asm-np.log10(1+((10**(mto-x)))**gamma)
def moustakas(x, mto, asm, gamma):
		return asm-np.log10(1+((10**(mto-x)))**gamma)
def second_order_2d(m, a, b, c, d, e, f):
		return a+(b*m)+(c*SFR)+(d*m**2)+(e*m*SFR)+(f*SFR**2)

calibration='N2'
#calibration='O3N2'



#This is the pixel plotting routine that was devised by Cappellari, I'm using it here because I want to be able to adapt it to various scales with vmin and vmax
def display_pixels(x, y, val, colormap, pixelsize=None, angle=None, **kwargs):
    """
    Display vectors of square pixels at coordinates (x,y) coloured with "val".
    An optional rotation angle can be applied to the pixels.
    
    """     
    pixelsize=0.66

    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)
    nx = np.round((xmax-xmin)/pixelsize) + 1
    ny = np.round((ymax-ymin)/pixelsize) + 1
    img = np.empty((nx,ny)) + np.nan
    j = np.round((x-xmin)/pixelsize).astype(np.int)
    k = np.round((y-ymin)/pixelsize).astype(np.int)
    img[j,k] = val
    

    ax = plt.gca()       
    f = ax.imshow(np.rot90(img), interpolation='nearest', cmap=colormap, vmin=colorbar_min, vmax=colorbar_max, alpha=1.0, 
                      extent=[xmin-pixelsize/2, xmax+pixelsize/2,
                              ymin-pixelsize/2, ymax+pixelsize/2], **kwargs)

    ax.minorticks_on()
    ax.tick_params(length=10, width=1, which='major')
    ax.tick_params(length=5, width=1, which='minor')
    #print(ax) 
    
    
    return f

def fitLine(x, y, alpha=0.05, newx=[], plotFlag=1):
    ''' Fit a curve to the data using a least squares 1st order polynomial fit '''
    # Summary data
    n = len(x)			   # number of samples     
    
    Sxx = np.sum(x**2) - np.sum(x)**2/n
#    Syy = np.sum(y**2) - np.sum(y)**2/n    # not needed here
    Sxy = np.sum(x*y) - np.sum(x)*np.sum(y)/n    
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    
    # Linefit
    b = Sxy/Sxx
    a = mean_y - b*mean_x
    
    # Residuals
    fit = lambda xx: a + b*xx    
    residuals = y - fit(x)
    
    var_res = np.sum(residuals**2)/(n-2)
    sd_res = np.sqrt(var_res)
    
    # Confidence intervals
    se_b = sd_res/np.sqrt(Sxx)
    se_a = sd_res*np.sqrt(np.sum(x**2)/(n*Sxx))
    
    df = n-2                            # degrees of freedom
    tval = stats.t.isf(alpha/2., df) 	# appropriate t value
    
    ci_a = a + tval*se_a*np.array([-1,1])
    ci_b = b + tval*se_b*np.array([-1,1])

    # create series of new test x-values to predict for
    npts = 100
    px = np.linspace(np.min(x),np.max(x),num=npts)
    
    se_fit     = lambda x: sd_res * np.sqrt(  1./n + (x-mean_x)**2/Sxx)
    se_predict = lambda x: sd_res * np.sqrt(1+1./n + (x-mean_x)**2/Sxx)
    
    print 'Summary: a={0:5.4f}+/-{1:5.4f}, b={2:5.4f}+/-{3:5.4f}'.format(a,tval*se_a,b,tval*se_b)
    print 'Confidence intervals: ci_a=({0:5.4f} - {1:5.4f}), ci_b=({2:5.4f} - {3:5.4f})'.format(ci_a[0], ci_a[1], ci_b[0], ci_b[1])
    print 'Residuals: variance = {0:5.4f}, standard deviation = {1:5.4f}'.format(var_res, sd_res)
    print 'alpha = {0:.3f}, tval = {1:5.4f}, df={2:d}'.format(alpha, tval, df)
    
    # Return info
    ri = {'residuals': residuals, 
        'var_res': var_res,
        'sd_res': sd_res,
        'alpha': alpha,
        'tval': tval,
        'df': df}
    
    if plotFlag == 1:
        # Plot the data
        #plt.figure()
        
        #plt.plot(px, fit(px),'k', label='Regression line')
        plt.plot(px, fit(px), color=cmap(color_value), linestyle='--', linewidth=3.0)
        #plt.plot(x,y,'r.', label='Sample observations')
        
        x.sort()
        limit = (1-alpha)*100
        #plt.plot(x, fit(x)+tval*se_fit(x), 'r--', label='Confidence limit ({0:.1f}%)'.format(limit))
        #plt.plot(x, fit(x)-tval*se_fit(x), 'r--')
        #plt.fill_between(x, fit(x)+tval*se_fit(x), fit(x)-tval*se_fit(x), color="none", hatch="\\", edgecolor=cmap(color_value))
        plt.fill_between(x, fit(x)+sd_res, fit(x)-sd_res, color="none", hatch="\\", edgecolor=cmap(color_value))
        #plt.fill_between(x, fit(x)+, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(color_value))
        
        #plt.plot(x, fit(x)+tval*se_predict(x), 'c--', label='Prediction limit ({0:.1f}%)'.format(limit))
        #plt.plot(x, fit(x)-tval*se_predict(x), 'c--')
        #plt.fill_between(x, fit(x)+tval*se_predict(x), fit(x)-tval*se_predict(x), color="none", hatch="\\", edgecolor=cmap(color_value))
        

        ##plt.xlabel('X values')
        ##plt.ylabel('Y values')
        ##plt.title('Linear regression and confidence limits')
        
        # configure legend
        ##plt.legend(loc=0)
        ##leg = plt.gca().get_legend()
        ##ltext = leg.get_texts()
        ##plt.setp(ltext, fontsize=10)

        # show the plot
        ##plt.show()
        
    if newx != []:
        try:
            newx.size
        except AttributeError:
            newx = np.array([newx])
    
        print 'Example: x = {0}+/-{1} => se_fit = {2:5.4f}, se_predict = {3:6.5f}'\
        .format(newx[0], tval*se_predict(newx[0]), se_fit(newx[0]), se_predict(newx[0]))
        
        newy = (fit(newx), fit(newx)-se_predict(newx), fit(newx)+se_predict(newx))
        return (a,b,(ci_a, ci_b), ri, newy)
    else:
        return (a,b,(ci_a, ci_b), ri)

names = np.array(['AGC191702', 'AGC202218', 'AGC212838', 'AGC220755', 'AGC220837', 'AGC220860', 'AGC221000', 'AGC221004', 'AGC225852', 'AGC225882', 'AGC227897'])

HOME = expanduser("~")
PAPER_FOLDER = '/Astro/Tex/2nd Paper FMR Version/'
PAPER_FOLDER = '/Astro/Tex/3rd Paper/'
PAPER_FOLDER = '/Astro/Tex/Jimmy0901/'
#PAPER_FOLDER = '/Astro/Tex/talk delete/'


print("Available plots:")
print("1. Balmer Decrement")
print("2. BPT *needs adjustments*")
print("3. Gradients *needs fixing*")
print("4. Big H_alpha Plot")
print("5. Big 12+log(O/H) Plot")
print("6. Mass-Metallicity")
print("7. Coloful HI Mass-Mass-Metallicity")
print("8. 3D HI Mass Stellar Mass and Metallicity *Fails to Execute*")
print("9. HI Mass-SFR")
print('10. BPT Distance-Gradient Slope *using hard coded data*')
print('11. 3D HI Mass Stellar Mass and SFR *Fails to Execute*')
print('12. 3D HI Mass Metallicity and SFR *Fails to Execute*')
print('13. Metallicity and SFR')
print('14. Luminosity-Metallicity')
print('15. HI Mass-Stellar Mass')
print('16. HI Mass-Metallicity')
print('17. Density-SFR relation *Fails to Execute*')
print('18. Gas Fraction-Metallicity')
print('19. Big Table of Everything')
print('20. Colorful Mass-Metallicity-SFR')
print('21. 4D Plot of Everything *Fails to Execute*')
print('22. SFR Comparison *Fails to Execute*')
print('23. Mass Comparison')
print('24. Field of Views')
print('25. FMR HI Mass')
print('26. FMR SFR')
print('27. Colorful Gas Fraction SFR')
print('28. Colorful Luminosity Metallicity and SFR')
print('29. Colorful Luminosity Metallicity and HI Mass')
print('30. Display TSWIFT fit')
print('31. Display Fake Galaxy For Appendix')
print('32. Consistency between IFU and SDSS Observables')
print('33. Lowest Scatter in SFR FMR')
print('34. Lowest Scatter in HI-Mass FMR')
print('35. FMR Luminosity SFR')
print('36. FMR Luminosity HI Mass')
print('37. Lowest Scatter in SFR FMR Luminosity')
print('38. Lowest Scatter in SFR HI Mass Luminosity')
print('39. SFR in Segmentation Map vs AoN 3 cut')
if (calibration == 'O3N2'):
	print('Evil mode: Activated')
	
if len(sys.argv) > 1:
	plot_choice = int(sys.argv[1])
else:
	plot_choice = int(raw_input("Which plot would you like me to make? "))

	
if plot_choice == 0:
	print("No plot chosen.\n")
elif plot_choice == 1:
	plot_to_make = 'balmer'
elif plot_choice == 2:
	plot_to_make = 'bpt'
elif plot_choice == 3:
	plot_to_make = 'gradients'
elif plot_choice == 4:
	plot_to_make = 'h_alpha'
elif plot_choice == 5:
	plot_to_make = 'metallicity'
elif plot_choice == 6:
	plot_to_make = 'mass_metallicity'
elif plot_choice == 7:
	plot_to_make = 'colorful_hi_mass_metallicity'
elif plot_choice == 8:
	plot_to_make = '3d_hi_mass_metallicity'
elif plot_choice == 9:
	plot_to_make = 'hi_mass_sfr'
elif plot_choice == 10:
	plot_to_make = 'distance_slope'
elif plot_choice == 11:
	plot_to_make = '3d_hi_mass_sfr'
elif plot_choice == 12:
	plot_to_make = '3d_hi_metallicity_sfr'
elif plot_choice == 13:
	plot_to_make = 'metallicity_sfr'
elif plot_choice == 14:
	plot_to_make = 'luminosity_metallicity'
elif plot_choice == 15:
	plot_to_make = 'hi_mass_stellar_mass'
elif plot_choice == 16:
	plot_to_make = 'hi_mass_metallicity'
elif plot_choice == 17:
	plot_to_make = 'hi_mass_density_sfr'
elif plot_choice == 18:
	plot_to_make = 'gas_fraction_metallicity'
elif plot_choice == 19:
	plot_to_make = 'table_of_everything'
elif plot_choice == 20:
	plot_to_make = 'colorful_mass_metallicity_sfr'
elif plot_choice == 21:
	plot_to_make = '4d'
elif plot_choice == 22:
	plot_to_make = 'sfr_comparison'
elif plot_choice == 23:
	plot_to_make = 'mass_comparison'
elif plot_choice == 24:
	plot_to_make = 'field_of_views'
elif plot_choice == 25:
	plot_to_make = 'fmr_hi_mass'
elif plot_choice == 26:
	plot_to_make = 'fmr_sfr'
elif plot_choice == 27:
	plot_to_make = 'colorful_gas_fraction_sfr'
elif plot_choice == 28:
	plot_to_make = 'colorful_luminosity_metallicity_sfr'
elif plot_choice == 29:
	plot_to_make = 'colorful_hi_mass_luminosity_metallicity'
elif plot_choice == 30:
	plot_to_make = 'D02_line_fit_output'
elif plot_choice == 31:
	plot_to_make = 'fake_h_alpha_map'
elif plot_choice == 32:
	plot_to_make = 'consistency_check'
elif plot_choice == 33:
	plot_to_make = 'lowest_scatter_sfr'
elif plot_choice == 34:
	plot_to_make = 'lowest_scatter_hi_mass'
elif plot_choice == 35:
	plot_to_make = 'fmr_lzr_sfr'
elif plot_choice == 36:
	plot_to_make = 'fmr_lzr_hi_mass'
elif plot_choice == 37:
	plot_to_make = 'lowest_scatter_sfr_lzr'
elif plot_choice == 38:
	plot_to_make = 'lowest_scatter_hi_mass_lzr'
elif plot_choice == 39:
	plot_to_make = 'segmentation_v_aon'


#To choose between the gandalf fits or the new SWIFT fits, Shouldn't need this anymore, but might keep around just in case I need to go back to Gandalf
#if (plot_to_make == 'h_alpha') or (plot_to_make == 'balmer') or (plot_to_make == 'metallicity') or (plot_to_make == 'bpt') or (plot_to_make == 'gradients'):
#fit_mode = "gandalf"
fit_mode = "swift"

#If the plot requires a s/n cut, then this routine will ask for it
#if (plot_to_make == 'h_alpha') or (plot_to_make == 'balmer') or (plot_to_make == 'metallicity') or (plot_to_make == 'bpt') or (plot_to_make == 'gradients') or (plot_to_make == 'mass_metallicity'):
#sn_response = raw_input("Signal to Noise cut to use?");
sn_response = 5
sncut = "sn"+str(sn_response)
#sncut = 'snall'

#names = np.array(['AGC191702', 'AGC212838', 'AGC220755', 'AGC220837', 'AGC221004', 'AGC225882', 'AGC227897', 'AGC191702', 'AGC191702', 'AGC191702', 'AGC191702'])
#names = np.array(['AGC191702', 'AGC191702', 'AGC212838', 'AGC220755', 'AGC220837', 'AGC191702', 'AGC191702', 'AGC221004', 'AGC191702', 'AGC225882', 'AGC227897'])

#speed_of_light = 299792.458
#H0 = 73
#omega_r_0 = 8.4e-5
#omega_m_0 = 0.27
#omega_l_0 = 0.73

speed_of_light = 299792.458
speed_of_light_error = 0.0 #defined quantity
H0 = 67.0
H0_error = 1.4
#omega_r_0 = 9.23640e-5 #some guy on the internet
#omega_r_0 = 5.37539e-5 #Dustin
#omega_r_0_error = 2.21994e-6 #Dustin
omega_r_0 = 0.0
omega_r_0_error = 0.0
omega_m_0 = 0.32
omega_m_0_error = 0.02
omega_l_0 = 0.68
omega_l_0_error = 0.02

solar_line = 8.69-0.03
solar_x_axis = 5.80

#Parameters to putz with
mass_bins_lower_limit = np.linspace(5.00, 11.75, 28)
hi_mass_bins_lower_limit = np.linspace(6.75, 11.00, 6) #np.array([0.0, 7.00, 7.60, 8.20, 8.80, 9.40, 10.00, 10.60])
sfr_bins_lower_limit = np.linspace(-3.25, 1.0, 6) #np.array([-1000.0, -5.00, -4.00, -3.00, -2.00, -1.00, 0.00, 1.00])
if (dog == 'y'):
	sfr_bins_lower_limit = sfr_bins_lower_limit+np.log10(7)
#sfr_bins_lower_limit = np.array([-3.5, -1.7, -0.8, 0.1, 1.0])
M_B_bins_lower_limit = np.linspace(-22, -8.0, 28) #np.array([ -22, -21.5, -21, -20.5 -20, -19.5, -19, -18.5, -18, -17.5, -17, -16.5, -16, -15.5, -15, -14.5, -14, -13.5, -13, -12.5, -12, -11.5, -11, -10.5, -10, -9.5, -9.0, -8.5, -8.0]) #29
binning_cut = 20
if (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'fmr_lzr_sfr'):
	bounds = sfr_bins_lower_limit
if (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'lowest_scatter_hi_mass'):
	bounds = hi_mass_bins_lower_limit
sdss_alpha_value = 0.75
uncertainty_box_alpha = 0.25
fontsize=24
bottom_margin = 0.17
left_margin = 0.15
right_margin = 0.95
top_margin = 0.96
wspace_margin = 0.20
hspace_margin = 0.20
labelpad = 35

#If the plot could be shaped differently between horizontal slides or vertical paper, this will pick between the two
if (plot_to_make == 'h_alpha') or (plot_to_make == 'balmer') or (plot_to_make == 'metallicity') or (plot_to_make == 'gradients'):
	#plotmode = raw_input("Format for slides or paper?");
	plotmode = "paper"
	if plotmode == "paper":
		fig = plt.figure(figsize=(13.3, 10))
	if plotmode == "slides":
		fig = plt.figure(figsize=(20, 5))
if (plot_to_make == 'bpt'):
	plotmode = 'paper'
	fig = plt.figure(figsize=(12, 12))

#If the plot requires offsets for either the axis, this will load that up.
if (plot_to_make == 'h_alpha') or (plot_to_make == 'balmer') or (plot_to_make == 'metallicity') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr'):
	xmin = np.array([ -9, -15,   -4, -12,  -7, -15,  -10,  -10,  -10, -12, -12])
	xmax = np.array([ 14,   8,  19,  11,  16,   8,  13,  13,  13,  11,  11])
	ymin = np.array([-5, -12, -12,  -7,  -7, -17,  -9,  -9, -14, -12, -12])
	ymax = np.array([ 18, 11,  11,  16,  16,   6,  14,  14,   9,  11,  11])
if (plot_to_make == 'h_alpha'):
		#make gray go right, add to x
		#make gray go up, add to y
		x_nudge = [1, 3, -5,0,-3, 1, 0,0, 6,-5,0]
		y_nudge = [0,-3,-10,0,-3,-5,-1,3,-2,2,0]
#If the plot requires offsets for the SDSS images, this will load that up.
if (plot_to_make == '?'):
	xmin_mod = np.array([-14,  -8, -23.5, -10, -10, -14, -13, -12, -16,  -9, -10])
	xmax_mod = np.array([ 14,  20,   4.5,  10,  10,  14,  15,  16,  12,  19,  10])
	ymin_mod = np.array([-16, -12,   -10, -10, -10, -16, -11, -17, -12, -17, -10])
	ymax_mod = np.array([ 12,  16,    18,  10,  10,  12,  17,  11,  16,  11,  10])

#The bpt diagram could be done in several ways, this chooses it
if (plot_to_make == 'bpt'):
	plot_type = 'regular'
	#plot_type = 'radius'
	#plot_type = 'metallicity'

#We need the effective radius for gradient plotting
if (plot_to_make == 'gradients'):
	r_e = np.array([6.309886, 9.553676, 6.555505, 3.238886, 7.059399, 13.46125, 6.214602, 3.268023, 4.208799, 8.229763, 1.002498]) #PetroR50_r

#I need to get distance calculations for each galaxy
if (plot_to_make == 'metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == '3d_hi_mass_metallicity')  or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '4d') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
	distance = np.array([])
	Jimmy_metallicity = np.array([])
	Jimmy_metallicity_error = np.array([])

if (plot_to_make == 'hi_mass_sfr') or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == '4d') or (plot_to_make == 'sfr_comparison') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
	Jimmy_sfr = np.array([])
	Jimmy_sfr_error = np.array([])
	
if (plot_to_make == 'hi_mass_sfr') or (plot_to_make == 'hi_mass_stellar_mass'):
	huang_catalog_file = HOME+"/Astro/Catalogs/huang_alfalfa_catalog.tsv"
	huang_catalog = np.genfromtxt(huang_catalog_file, dtype=None, names=['AGCNr','Poor_fit','RAJ2000', 'DEJ2000', 'f_FUV', 'FUV', 'e_FUV', 'NUV', 'e_NUV', 'rmag', 'e_rmag', 'u-r', 'e_u-r', 'sFlag', 'Dist', 'e_Dist', 'logMHI', 'e_logMHI', 'logM', 'e_logM', 'logSFR', 'e_logSFR', 'SDSS'], delimiter='\t', skiprows=80)

if (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr'):
	Saintonge_names = ['AGC100244','AGC112521','AGC110482','AGC111946','AGC122206','AGC122226','AGC3755','AGC3775','AGC181491','AGC180167','AGC182462','AGC205072','AGC202042','AGC210023','AGC212064','AGC225882','AGC220242','AGC220478','AGC223205','AGC10282','AGC10351','AGC253923']
	#Saintonge_metallicity = np.array([8.17, 7.43, 7.82, 7.35, 8.44, 7.96, 8.16, 7.91, 8.31, 7.91, 7.97, 8.05, 7.90, 8.01, 7.95, 7.59, 7.63, 7.95, 8.24, 8.24, 8.59, 7.85])
	Saintonge_mass = np.array([8.78, 7.08, 7.63, 6.96, 8.87, 7.97, 8.17, 9.39, 9.08, 8.66, 8.27, 8.90, 7.62, 7.51, 9.03, 7.79, 7.93, 7.63, 7.25, 8.50, 9.10, 8.68])-0.3 #Convert Salpeter to Kroupa
	Saintonge_hi_mass = np.array([8.49, 6.90, 7.18, 6.94, 8.55, 7.35, 7.66, 9.16, 8.92, 9.06, 8.51, 8.24, 8.93, 7.61, 9.09, 7.82, 8.24, 7.99, 8.08, 8.50, 8.73, 8.07])
	Saintonge_B = np.array([8.54, 6.56, 7.39, 6.72, 8.32, 7.73, 7.93, 9.15, 9.13, 9.11, 8.07, 7.94, 8.65, 7.85, 9.23, 7.78, 8.00, 7.96, 7.72, 8.26, 8.80, 8.21])
	Saintonge_B = np.array([15.50, 18.52, 16.57, 18.22, 16.70, 15.73, 14.10, 15.53, 16.63, 15.34, 17.23, 18.97, 17.39, 15.77, 15.58, 17.15, 16.62, 16.80, 17.33, 16.45, 14.74, 18.13])
	Saintonge_D = np.array([19.0, 7.2, 7.2, 7.2, 21.2, 6.8, 4.3, 32.3, 61.2, 32.7, 23.7, 45.2, 50.0, 9.8, 43.5, 16.7, 16.7, 16.7, 16.7, 19.8, 17.3, 41.3])
	Saintonge_M_B = Saintonge_B + 5 - (5*np.log10(Saintonge_D*1e6))
	Saintonge_Ha = np.array([2.863, 2.784, 2.803, 2.862, 1.733, 2.935, 2.419, 2.856, 2.914, 1.699, 2.464, 2.376, 2.277, 1.656, 1.916, 2.900, 2.538, 2.819, 2.953, 2.832, 2.864, 2.796])
	Saintonge_OIII = np.array([3.089, 1.893, 3.612, 0.854, 3.315, 3.876, 6.761, 3.886, 2.210, 4.330, 2.352, 4.677, 4.465, 3.387, 3.193, 2.043, 2.645, 3.370, 2.604, 2.996, 3.474, 2.741])
	Saintonge_NII = np.array([0.120, 0.045, 0.044, 0.119, 0.184, 0.095, 0.125, 0.167, 0.334, 0.107, 0.086, 0.066, 0.064, 0.098, 0.196, 0.124, 0.0, 0.143, 0.162, 0.152, 0.135, 0.170])
	Saintonge_Hb = np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
	Saintonge_N2 = np.log10(Saintonge_NII/Saintonge_Ha)
	Saintonge_PP04_N2 = 9.37 + (2.03*Saintonge_N2) + (1.26*Saintonge_N2**2) + (0.32*Saintonge_N2**3)
	Saintonge_D02_N2 = 9.12 + 0.73*Saintonge_N2
	Saintonge_PP04_O3N2 = 8.73-(0.32*(np.log10((Saintonge_OIII/Saintonge_Hb)/(Saintonge_NII/Saintonge_Ha))))
	if (calibration == 'N2'):
		Saintonge_metallicity = Saintonge_D02_N2
	if (calibration == 'O3N2'):
		Saintonge_metallicity = Saintonge_PP04_O3N2
		
	
if (plot_to_make == 'mass_metallicity') or (plot_to_make == 'luminosity_metallicity'):
	#Berg_PP04 = [8.14, 8.34, 8.29, 8.18, 8.15, 8.15, 8.29, 8.19, 8.20, 8.19, 8.16, 8.26, 8.04, 8.09, 8.17, 8.07, 8.31, 8.24, 8.17, 8.10, 8.52, 8.05, 7.97, 8.19, 8.22, 7.96, 8.17, 8.14, 8.13, 8.38, 8.08]
	#Berg_names = ['UGC 521 A', 'UGC 695E', 'UGC 1056 A', 'UGC 1176 A', 'NGC 784 A', 'UGC 2716 A', 'NGC 2537 A', 'UGC 4278 A', 'NGC 2552 A', 'UGC 4393 B', 'UGC 5139 A', 'IC 559 A', 'UGC 5272 A', 'UGC 5340 A', 'UGC 5423 A', 'UGC 5797 A', 'UGC 5923 A', 'NGC 3738 A', 'NGC 3741 A', 'UGC 6817 A', 'NGC 4163 A', 'CGCG 269-049 A', 'UGC 7577 A', 'NGC 4449 A', 'UGC 7605 A', 'NGC 4656 A', 'UGC 8201 A', 'UGC 8508 A ', 'UGC 8638 A', 'UGC 8837 A', 'NGC 5477 A']
	#Berg_mass = [7.96, 8.08, 8.62, 8.48, 8.48, 8.13, 9.10, 8.50, 8.69, 9.43, 7.39, 7.86, 8.00, 7.97, 7.77, 7.75, 8.29, 8.50, 7.05, 6.97, 7.61, 5.90, 7.50, 9.25, 7.12, 9.04, 7.82, 7.00, 7.57, 8.41, 8.15]
	#Berg_M_B = [-15.16, -15.13, -15.09, -15.48, -16.50, -15.31, -17.14, -16.36, -16.72, -17.67, -14.42, -14.12, -14.98, -15.83, -13.77, -14.56, -14.70, -16.51, -13.18, -13.70, -13.65, -10.83, -14.12, -18.02, -13.49, -18.75, -15.17, -13.03, -13.77, -15.92, -15.22]
	#Berg_direct_metallicity = [7.67, 7.69, 7.89, 7.97, 8.03, 7.97, 8.24, 7.69, 8.15, 8.06, 7.92, 8.07, 7.87, 7.20, 7.77, 7.96, 7.79, 8.04, 7.68, 7.53, 7.56, 7.47, 7.97, 8.32, 7.66, 8.09, 7.80, 7.76, 7.97, 7.82, 7.95]
	#Berg_N2 = [7.87, 8.13, 8.15, 8.04, 8.02, 7.99, 8.27, 7.88, 8.12, 8.16, 7.94, 8.14, 7.73, 7.51, 8.02, 7.92, 8.26, 8.15, 7.87, 7.69, 8.18, 7.72, 7.94, 8.08, 7.94, 7.62, 7.85, 7.85, 7.93, 8.20, 7.81]
	#if (calibration == 'N2'):
	#	Berg_metallicity = Berg_N2
	#if (calibration == 'O3N2'):
	#	Berg_metallicity = Berg_PP04
	#print(len(Berg_metallicity))
	Berg_names = ['UGC 521 A', 'UGC695 E', 'UGC1056 A', 'UGC 1056 B', 'UGC 1176 A', 'NGC 784 A','NGC 784 B', 'UGC 2716 A', 'KKH 037 A', 'NGC 2537 A', 'NGC 2537 B', 'UGC 4278 B', 'UGC4278 A', 'NGC 2552 A', 'UGC 4393 B', 'UGC4393 C', 'CGCG 035-007 A', 'UGC 5139 A', 'IC 559 A', 'UGC 5272 A', 'UGC 5340 A', 'UGC 5423 A', 'UGC 5423 B', 'UGC 5672 A', 'UGC 5692 A', 'UGC 5797 A', 'UGC 5923 A', 'NGC 3741 A', 'NGC 3738 A', 'NGC 3738 B', 'UGC 6817 A', 'UGC 6900 A', 'NGC 4163 A', 'CGCG 269-049 C', 'CGCG 269-049 A', 'UGC 7577 A', 'NGC 4449 C', 'NGC 449 B', 'NGC 449 A', 'UGC 7605 A ', 'UGC 7639 A ', 'NGC 4656 A', 'UGC 8201 A', 'UGC 8245 A', 'UGC 8508 A', 'UGC 8638 A', 'UGC 8638 B', 'UGC 8837 A', 'NGC 5477 A', 'UGC 9405 A', 'UGC 10818 A', 'KKH 098 A']
	Berg_M_B = np.array([-15.16, -15.13, -15.09, -15.09, -15.48, -16.5, -16.5, -15.31, -11.98, -17.14, -17.14, -16.36, -16.36, -16.72, -17.67, -17.67, -13.38, -14.42, -14.12, -14.98, -15.83, -13.77, -13.77, -14.73, -14.68, -14.56, -14.70, -13.18, -16.51, -16.51, -13.70, -14.62, -13.65, -10.83, -10.83, -14.12, -18.02, -18.02, -18.02, -13.49, -15.55, -18.75, -15.17, -13.67, -13.03, -13.77, -13.77, -15.92, -15.22, -14.97, -18.59, -11.10])
	Berg_mass = np.array([7.96, 8.08, 8.62, 8.62, 8.48, 8.48, 8.48, 8.13, 7.01, 9.10, 9.10, 8.50, 8.50, 8.69, 9.43, 9.43, 7.69, 7.39, 7.86, 8.00, 7.97, 7.77, 7.77, 8.38, 8.16, 7.75, 8.29, 7.05, 8.50, 8.50, 6.97, 8.19, 7.61, 5.90, 5.90, 7.50, 9.25, 9.25, 9.25, 7.12, 8.25, 9.04, 7.82, 7.53, 7.00, 7.57, 7.57, 8.41, 8.15, 7.97, 9.45, 6.72])-0.3 #convert to a Kroupa IMF
	Berg_OIII = np.array([3.64, 1.62, 2.35, 3.27, 3.61, 4.13, 3.32, 4.26, 0.53, 2.14, 1.78, 1.91, 2.53, 2.91, 3.28, 2.59, 2.03, 3.66, 2.80, 4.94, 1.89, 3.49, 3.71, 2.51, 1.70, 5.53, 2.48, 2.84, 2.96, 3.11, 2.89, 0.60, 0.49, 1.53, 2.51, 5.25, 2.33, 2.54, 3.46, 2.33, 1.33, 6.71, 2.94, 1.25, 3.25, 4.17, 4.15, 1.26, 4.64, 1.44, 2.22, 1.91])
	Berg_Hb = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
	Berg_Ha = np.array([2.77, 2.86, 2.80, 2.82, 2.82, 2.89, 2.81, 2.86, 2.80, 2.87, 2.79, 2.80, 2.86, 2.86, 2.86, 2.87, 2.84, 2.83, 2.86, 2.83, 2.81, 2.86, 2.86, 2.87, 2.79, 2.84, 2.78, 2.83, 2.83, 2.82, 2.83, 2.83, 2.75, 2.79, 2.83, 2.79, 2.84, 2.86, 2.83, 2.83, 2.83, 2.86, 2.81, 2.83, 2.83, 2.82, 2.81, 2.83, 2.84, 2.79, 2.79, 2.79])
	Berg_NII = np.array([0.05, 0.125, 0.14, 0.12, 0.122, 0.109, 0.080, 0.092, 0.14, 0.48, 0.41, 0.060, 0.054, 0.25, 0.30, 0.31, 0.21, 0.07, 0.15, 0.035, 0.016, 0.120, 0.082, 0.23, 0.38, 0.09, 0.24, 0.051, 0.19, 0.20, 0.03, 0.52, 0.13, 0.0, 0.033, 0.09, 0.23, 0.205, 0.163, 0.066, 0.21, 0.023, 0.04, 0.12, 0.05, 0.07, 0.06, 0.202, 0.047, 0.40, 0.25, 0.08])
	Berg_direct_metallicity = np.array([7.67, 7.69, 7.89, 7.98, 7.97, 8.03, 7.92, 7.97, 0.0, 8.24, 7.98, 7.69, 7.69, 8.15, 8.06, 7.99, 0.0, 7.92, 8.07, 7.87, 7.20, 7.77, 7.82, 0.0, 0.0, 7.96, 7.79, 7.68, 8.04, 8.03, 7.53, 0.0, 7.56, 0.0, 7.47, 7.97, 8.13, 8.16, 8.32, 7.66, 0.0, 8.09, 7.80, 0.0, 7.76, 7.97, 7.91, 7.82, 7.95, 0.0, 0.0, 0.0])
	Berg_N2 = np.log10(Berg_NII/Berg_Ha)
	Berg_PP04_N2 = 9.37 + (2.03*Berg_N2) + (1.26*Berg_N2**2) + (0.32*Berg_N2**3)
	Berg_D02_N2 = 9.12 + 0.73*Berg_N2
	Berg_PP04_O3N2 = 8.73-(0.32*(np.log10((Berg_OIII/Berg_Hb)/(Berg_NII/Berg_Ha))))
	if (calibration == 'N2'):
		#Berg_D02_N2 = 9.12 + 0.73*(np.log10(Berg_NII/Berg_Ha))
		Berg_metallicity = Berg_D02_N2
	if (calibration == 'O3N2'):
		#Berg_PP04_O3N2 = 8.73-(0.32*(np.log10((Berg_OIII/Berg_Hb)/(Berg_NII/Berg_Ha))))
		Berg_metallicity = Berg_PP04_O3N2



if (plot_to_make == 'mass_metallicity') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr'):
	James_names = ['KJ 1.0','KJ 1.1','KJ 3.0','KJ 3.1','KJ 5.0','KJ 7.0','KJ 18.0','KJ 18.1','KJ 22.0','KJ 22,1','KJ 29.0','KJ 29.2','KJ 37.0','KJ 53.0','KJ 78.0']
	James_Hb = np.array([100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100])
	James_OIII = np.array([362, 477, 528, 207, 132, 341, 283, 188, 512, 458, 151, 182, 419, 338, 357]) #5007
	James_Ha = np.array([335, 265, 256, 346, 309, 292, 308, 371, 322, 313, 309, 322, 368, 280, 265])
	James_NII = np.array([8.29, 4.76, 6.00, 25.1, 7.01, 13, 8.53, 44.1, 8.65, 6.84, 6.42, 2.58, 12.5, 2.83, 4.68]) #6584
	James_LHa = [31.48, 98.20, 56.24, 6.70, 1.08, 0.99, 1.18, 0.39, 32.19, 20.03, 19.92, 31.34, 28.82, 0.47, 0.43] #x10^38 erg/s
	James_sfr = np.log10(np.array([0.0249, 0.0776, 0.0444, 0.0053, 0.0009, 0.0008, 0.0009, 0.0003, 0.0254, 0.0158, 0.0157, 0.0248, 0.0228, 0.0004, 0.0003]))
	James_mass = np.log10([1.1e9, 1.1e9, 8.0e8, 8.0e8, 2.7e8, 3.0e7, 2.7e8, 2.7e8, 2.1e8, 2.1e8, 6.1e8, 6.1e8, 4.9e8, 2.8e6, 3.3e6]) - 0.15 #convert Diet Salpeter to Kroupa IMF
	James_M_B = [-17.98, -17.98, -17.14, -17.14, -13.12, -13.28, -15.49, -15.49, -15.92, -15.92, -17.31, -17.31, -16.24, -11.17, -10.72]
	James_direct_metallicity = [7.74, 7.72, 7.99, 7.59, 7.67, 8.00, 7.67, 6.87, 7.98, 7.91, 7.58, 7.83, 7.83, 7.45, 7.62]
	James_PP04_O3N2 = 8.73-(0.32*(np.log10((James_OIII/James_Hb)/(James_NII/James_Ha))))
	James_N2 = np.log10(James_NII/James_Ha)
	James_PP04_N2 = 9.37 + (2.03*James_N2) + (1.26*James_N2**2) + (0.32*James_N2**3)
	James_D02_N2 = 9.12 + 0.73*James_N2
	if (calibration == 'N2'):
		James_metallicity = James_D02_N2
	if (calibration == 'O3N2'):
		James_metallicity = James_PP04_O3N2

if (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr'):
	SHIELD_names = ['AGC748778','AGC112521','AGC110482','AGC111946','AGC111977','AGC111164','AGC174585','AGC174605','AGC182595','AGC731457','AGC749237','AGC749241']
	SHIELD_mass = np.array([5.82, 6.93, 7.35, 6.83, 7.11, 6.64, 6.48, 6.94, 7.25, 0.0, 7.60, 0.0])-0.30#to adjust from Salpeter to Kroupa
	SHIELD_M_B = [-10.34, -10.59, -13.02, -11.87, -12.60, -11.16, -11.32, -12.22, -12.70, -13.73, -14.12, -10.25]
	SHIELD_sfr = np.array([-4.92, -4.10, -2.67, -3.10, -3.06, -3.40, -3.21, -3.19, -2.65, -2.52, -2.34, -5.00])
	SHIELD_hi_mass = np.array([6.72, 7.11, 7.31, 7.11, 6.85, 6.64, 6.94, 7.19, 6.92, 7.09, 7.49, 6.74])
	#SHIELD_mass = np.log10(np.array([3e6, 7e6, 55e6, 17e6, 37e6, 10e6, 9e6, 28e6, 50e6, 65e6, 53e6, 4e6])) #McQuinn:15 stellar mass
	#SHIELD_sfr = np.log10(np.array([0.32e-3, 0.75e-3, 5.5e-3, 1.7e-3, 3.8e-3, 1.0e-3, 0.90e-3, 2.8e-3, 5.1e-3, 6.6e-3, 5.4e-3, 0.36e-3])) #McQuinn:15 sfr
	#SHIELD_hi_mass = np.log10(np.array([4.5e6, 7.0e6, 19e6, 17e6, 7.1e6, 10e6, 7.9e6, 19e6, 8.1e6, 18e6, 57e6, 5.7e6])) #McQuinn:15 hi mass
	SHIELD_metallicity = [0.0, 7.33, 7.79, 7.86, 7.80, 7.59, 0.0, 0.0, 7.75, 8.0, 7.95, 0.0]
	SHIELD_Hb = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
	SHIELD_OIII = np.array([0.0, 1.701, 3.9485, 0.838, 3.375, 2.339, 0.0, 0.0, 4.159, 1.784, 2.23566667, 0.0])
	SHIELD_Ha = np.array([1.0, 2.860, 2.8225, 2.860, 2.860, 2.860, 1.0, 1.0, 2.786, 2.860, 2.860, 1.0])
	SHIELD_NII = np.array([999.0, 0.015, 0.040, 0.226, 0.176, 0.062, 999.0, 999.0, 0.064, 0.143, 0.14833333, 999.0])
	SHIELD_PP04_O3N2 = 8.73-(0.32*(np.log10((SHIELD_OIII/SHIELD_Hb)/(SHIELD_NII/SHIELD_Ha))))
	SHIELD_N2 = np.log10(SHIELD_NII/SHIELD_Ha)
	SHIELD_PP04_N2 = 9.37 + (2.03*SHIELD_N2) + (1.26*SHIELD_N2**2) + (0.32*SHIELD_N2**3)
	SHIELD_D02_N2 = 9.12 + 0.73*SHIELD_N2
	if (calibration == 'O3N2'):
		SHIELD_metallicity = SHIELD_PP04_O3N2
	if (calibration == 'N2'):
		SHIELD_metallicity = SHIELD_D02_N2

if (plot_to_make == 'mass_metallicity') or (plot_to_make == 'luminosity_metallicity'):
	Lee_names = ['DDO 154', 'DDO 165', 'DDO 53', 'GR 8', 'Ho I', 'Ho II', 'Ho IX', 'IC 1613', 'IC 2574', 'IC 5152', 'Leo A', 'M81 dwB', 'NGC 1569', 'NGC 1705', 'NGC 2366', 'NGC 3109', 'NGC 3738', 'NGC 4214', 'NGC 4449', 'NGC 5408', 'NGC 55', 'NGC 6682', 'Peg DIG ', 'Sextans A', 'Sextans B', 'UGC 6456', 'WLM']
	Lee_metallicity = np.array([7.67, 0, 7.62, 7.64, 7.70, 7.93, 0, 7.62, 8.15, 7.92, 7.35, 7.98, 8.19, 8.21, 7.91, 8.06, 8.23, 8.22, 8.31, 7.89, 8.05, 8.11, 7.93, 7.54, 7.53, 7.69, 7.83])
	Lee_M_B = [1.88-15.28, 2.46-17.62, 2.36-15.73, 2.52-14.72, 1.83-16.32, 1.82-18.53, 2.07-15.75, 1.42-15.92, 2.09-19.43, 2.78-18.28, 1.84-13.34, 1.96-16.50, 0.93-19.36, 2.62-18.44, 2.18-18.18, 1.11-17.61, 2.45-19.03, 2.55-19.53, 2.98-21.12, 2.29-18.79, 1.58-19.88, 2.98-18.14, 2.90-15.40, 0.85-14.85, 1.73-15.83, 1.83-15.87, 1.88-15.78]
	Lee_mass = [6.68, 7.73, 6.98, 6.62, 7.08, 7.96, 6.89, 6.82, 8.39, 8.11, 5.89, 7.19, 8.07, 8.13, 7.91, 7.41, 8.32, 8.55, 9.29, 8.19, 8.44, 8.10, 6.98, 6.24, 6.86, 6.90, 6.88]
	x_variable = Lee_metallicity
	#x_variable = 7.65
	a = 230.7820
	b = -75.79752
	c = 8.526986
	d = -0.3162894
	Lee_metallicity_PP04_O3N2_scale = a+(b*x_variable)+(c*x_variable**2)+(d*x_variable**3)
	a = 253.0031
	b = -87.03697
	c = 10.241880
	d = -0.3984731
	Lee_metallicity_D02_N2_scale = a+(b*x_variable)+(c*x_variable**2)+(d*x_variable**3)
	#print('XMP rescaled: '+str(Lee_metallicity_PP04_O3N2_scale))
	
if (plot_to_make == 'h_alpha') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == '3d_hi_mass_metallicity') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == '3d_hi_mass_sfr') or (plot_to_make == '3d_hi_metallicity_sfr')  or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '4d') or (plot_to_make == 'sfr_comparison') or (plot_to_make == 'mass_comparison') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
	z_array = np.array([])
	z_array_error = np.array([])
	E_B_minus_V_array = np.array([])
	
#Pull in the redshift of the galaxy which is needed for luminosity calculations for SFR estimates
if (plot_to_make == 'h_alpha') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == '3d_hi_mass_metallicity') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == '3d_hi_mass_sfr') or (plot_to_make == '3d_hi_metallicity_sfr')  or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '4d') or (plot_to_make == 'sfr_comparison') or (plot_to_make == 'mass_comparison') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
	for index, galaxy in enumerate(names):
		environmental_variables=np.loadtxt(HOME+"/Astro/reduced/"+galaxy+"pt1sof/"+galaxy+"_env.sh", delimiter="=", skiprows=2, dtype={'names': ('dummy', 'value'), 'formats': ('S10', 'f4')})
		temp_variable = environmental_variables[1]
		redshift = (temp_variable[1])
		#z_array = np.append(redshift, z_array)
		#z_array_error = (z_array*0.1)
	#This is all well and good, but it's too smart, and it turns out the stupid way works better or at least works just the same, but faster.
	alfalfa_sdss_catalog_file = HOME+"/Astro/Catalogs/alfalfa_sdss_averaged.csv"
	alfalfa_sdss_catalog = np.genfromtxt(alfalfa_sdss_catalog_file, dtype=None, names=['AGCNr','hi_mass','distance','oh_p50','lgm_tot_p50','sfr_tot_p50','h_beta_flux','h_beta_flux_err','oiii_5007_flux','oiii_5007_flux_err','h_alpha_flux','h_alpha_flux_err','nii_6584_flux','nii_6584_flux_err','petroR90_r','petroR90Err_r','petroR50_r','petroR50Err_r','petroMag_u','petroMagErr_u','petroMag_g','petroMagErr_g','petroMag_r','petroMagErr_r','petroMag_i','petroMagErr_i','petroMag_z','petroMagErr_z','z','zErr','fiberMag_r','fiberMagErr_r'], delimiter=",", skiprows=1)
	for index, galaxy in enumerate(names):
		#print(galaxy[3:])
		z_detected = 0
		for sdss_index, sdss_line in enumerate(alfalfa_sdss_catalog):
			if str(sdss_line['AGCNr']) == galaxy[3:]:
				#print(sdss_line['AGCNr'])
				z_array = np.append(sdss_line['z'], z_array)
				z_array_error = np.append(sdss_line['zErr'], z_array_error)
				#print(z_array_error)
				z_detected = 1
		if z_detected == 0:
			environmental_variables=np.loadtxt(HOME+"/Astro/reduced/"+galaxy+"pt1sof/"+galaxy+"_env.sh", delimiter="=", skiprows=2, dtype={'names': ('dummy', 'value'), 'formats': ('S10', 'f4')})
			temp_variable = environmental_variables[1]
			redshift = temp_variable[1]
			z_array = np.append(redshift, z_array)
			z_array_error = np.append(sdss_line['zErr'], z_array_error)
			#z_array_error = (z_array*0.1)

if (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == '3d_hi_mass_metallicity') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == '3d_hi_mass_sfr') or (plot_to_make == '3d_hi_metallicity_sfr') or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr')  or (plot_to_make == '4d') or (plot_to_make == 'mass_comparison') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'h_alpha') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
	#The corrected redshift corrects the velocity to the 3K CMB background, need to load b and l for each galaxy though
	#redshift correction v+vapex{sin(b)sin(bapex)+cos(b)cos(bapex)cos(l)cos(lapex)}
	vapex = 371.0
	bapex = 48.26
	lapex = 264.14
	l = np.array([224.858050, 233.270430, 251.072997, 289.191411, 291.796902, 293.501543, 298.791031, 298.55808, 310.088112, 261.1779, 301.9501]) #Galactic Longitude
	b = np.array([32.775573, 52.261359, 65.847993, 70.171420, 70.604647, 69.624125, 71.304957, 73.03325, 73.483619, 72.3239, 69.7147]) #Galactic Lattitude
	corrected_redshift = np.zeros_like(b)
	for z_index, z_variable in enumerate(z_array):
		corrected_redshift[z_index] = z_array[z_index] + vapex/speed_of_light*(math.sin(math.radians(b[z_index]))*math.sin(math.radians(bapex))+math.cos(math.radians(b[z_index]))*math.cos(math.radians(bapex))*math.cos(math.radians(l[z_index]-lapex)))
	redshift = corrected_redshift
	#print(z_array)
	#print(corrected_redshift)
	#z_array = corrected_redshift
	q0 = omega_r_0 + omega_m_0/2 - omega_l_0
	#print('q0: '+str(q0))
	q0_error = np.sqrt((omega_r_0_error**2)+((omega_m_0_error**2)/4)+(omega_l_0_error**2))
	#print('q0_error: '+str(q0_error))
	#Dp = (speed_of_light*z_array/H0)*(1-(((1+q0)/2)*z_array))
	#print(Dp)
	Dp = np.array([8.7, 19.6, 10.3, 16.4, 16.4, 16.4, 16.5, 16.7, 16.6, 23.6, 16.6])
	Dp_error = np.array(['AGC191702', 'AGC202218', 'AGC212838', 2.5, 2.5, 2.5, 'AGC221000', 1.7, 'AGC225852', 'AGC225882', 'AGC227897'])
	#print(Dp)
	#print(names)
	#print('Dp: '+str(Dp))
	#print('z_array_error: '+str(z_array_error))
	#Dp_error = np.sqrt(((((speed_of_light/H0)*(1-(((1+q0))*z_array)))**2)*z_array_error**2)+((((-speed_of_light*z_array/(H0*2))*(1-(((1+q0)/2)*z_array)))**2)*H0_error**2)+((((-speed_of_light*(z_array**2))/(2*H0))**2)*q0_error**2))
	#wrtz = ((speed_of_light/H0)*(1-(z_array*(1+q0))))
	#wrtH0 = (-((speed_of_light*z_array)/(H0**2))*(1-(((1+q0)/2)*z_array)))
	#wrtq0 = (-(speed_of_light*z_array**2)/(2*H0))
	#Dp_error = np.sqrt(((((speed_of_light/H0)*(1-(z_array*(1+q0))))**2)*z_array_error**2)+(((-((speed_of_light*z_array)/(H0**2))*(1-(((1+q0)/2)*z_array)))**2)*H0_error**2)+(((-(speed_of_light*z_array**2)/(2*H0))**2)*q0_error**2))
	#print('Dp_error: '+str(Dp_error))
	luminosity_distance = Dp#*(1+z_array)
	#luminosity_distance_error = np.sqrt(((Dp**2)*z_array_error**2)+(((1+z_array)**2)*Dp_error**2))
	luminosity_distance_error = np.array([2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 1.7, 2.5, 2.5, 2.5])
	#print('luminosity_distance_error: '+str(luminosity_distance_error))
	#petrosian magnitudes pulled form SDSS for each object in several bands
	#u = np.array([18.43, 16.18035, 18.81257, 17.75595, 15.97912, 17.25299, 16.18035, 17.20, 17.8402, 18.498, 19.10]) #SDSS catalog pulled
	#g = np.array([17.69, 15.34365, 18.40118, 16.59541, 14.76143, 16.73954, 15.34365, 16.67, 16.84656, 17.86655, 18.61]) #SDSS catalog pulled
	#r = np.array([17.39, 14.8, 18.27659, 16.17616, 14.24938, 16.4022, 14.8, 16.36, 16.46039, 17.82213, 18.44]) #SDSS catalog pulled
	#i = np.array([17.54, 14.84343, 18.31006, 16.04263, 13.93469, 16.21586, 14.84343, 16.42, 16.30921, 18.28569, 18.30]) #SDSS catalog pulled
	#u = np.array([26.2023, 22.8998, 24.3396, 23.857, 21.8131, 17.25299, 21.956, 23.2225, 23.4894, 25.1832, 25.3792]) #Source Extractor Fits missing light
	#g = np.array([26.0458, 22.0601, 23.8507, 22.8515, 20.5835, 16.73954, 21.033, 22.5684, 22.7595, 25.0957, 24.7202]) #Source Extractor Fits missing light
	#r = np.array([26.0049, 21.6824, 23.6969, 22.4939, 20.0879, 16.4022, 20.6258, 22.3486, 22.4769, 25.1647, 24.5838]) #Source Extractor Fits missing light
	#i = np.array([26.1788, 21.5526, 23.7786, 22.256, 19.8538, 16.21586, 20.4209, 22.1909, 22.3422, 25.3287, 24.6481]) #Source Extractor Fits missing light
	u = np.array([23.5818, 22.7641, 23.4786, 23.4888, 21.5112, 23.6755, 21.9577, 22.7118, 23.3406, 24.2919, 25.1304]) - 5.7453 #Source Extractor had the wrong zero point, probably from an older version of SDSS images, so instead of re-running Source Extractor on everything, I'll just put in a correction here.
	u_error = np.array([0.0625, 0.043, 0.0619, 0.0626, 0.0247, 0.0656, 0.0285, 0.0454, 0.0545, 0.0926, 0.1247])
	g = np.array([23.0497, 21.904, 22.8426, 22.5291, 20.3211, 22.5725, 21.0304, 21.9029, 22.6303, 23.7147, 24.4843]) - 5.7453
	g_error = np.array([0.046, 0.0272, 0.0421, 0.0362, 0.0131, 0.0368, 0.0181, 0.0271, 0.0378, 0.0624, 0.0887])
	r = np.array([22.792, 21.4984, 22.6522, 22.1278, 19.8047, 22.2631, 20.6143, 21.5129, 22.3208, 23.7029, 24.3246]) - 5.7453
	r_error = np.array([0.0419, 0.0228, 0.0397, 0.0307, 0.0104, 0.0324, 0.015, 0.0231, 0.0331, 0.0646, 0.0841])
	i = np.array([22.8993, 21.3417, 22.4427, 21.8419, 19.549, 22.1242, 20.4038, 21.3249, 22.151, 23.5616, 24.3239]) - 5.7453
	i_error = np.array([0.0442, 0.0213, 0.0367, 0.0269, 0.0093, 0.0305, 0.0136, 0.021, 0.0305, 0.0613, 0.0847])
	#B = g + 0.51 + 0.60*(g-r) #Windhorst, R. W., et al. 1991, ApJ, 380, 362 #https://www.astro.umd.edu/~ssm/ASTR620/mags.html
	#B = u - 0.8116*(u - g) + 0.1313 #Lupton (2005) http://classic.sdss.org/dr7/algorithms/sdssUBVRITransform.html
	B = g + 0.3130*(g - r) + 0.2271 #Lupton (2005)
	#print('B: '+str(B))
	B_error = np.sqrt(((1.313**2)*g_error**2)+(((-0.313)**2)*r_error**2))
	#print('B_error: '+str(B_error))
	huang_overlap = [3,6,7,8] #the galaxies which overlap with the huang subsample
	huang_r = np.array([16.01, 14.80, 15.53, 16.35])
	huang_g = np.array([16.01+1.32, 14.80+1.29, 15.53+1.13, 16.35+1.02])
	huang_D = np.array([16.4, 16.6, 16.7, 16.6])
	#huang_i = huang_r
	huang_i = i[huang_overlap]
	huang_M_i = huang_i + 5 - (5*np.log10(huang_D*1e6))
	#M = m + 5 - 5 log D
	#M = m - 5((log(D_L)) - 1)
	M_r = r + 5 - (5*np.log10(luminosity_distance*1e6))
	M_i = i + 5 - (5*np.log10(luminosity_distance*1e6))
	M_i_error = np.sqrt((i_error**2)+((((-5/luminosity_distance)*np.log10(np.e))**2)*luminosity_distance_error**2))
	#print('M_i_error: '+str(M_i_error))
	M_B = B + 5 - (5*np.log10(luminosity_distance*1e6))
	M_B_error = np.sqrt((B_error**2)+((((-5/luminosity_distance)*np.log10(np.e))**2)*luminosity_distance_error**2))
	#print('M_B_error: '+str(M_B_error))
	Jimmy_M_B = M_B
	Jimmy_M_B_error = M_B_error
	#Pull in petrosian u, g, r, i, and z and then calculate the luminosity distance from redshift
	bell_03_m_l = 0.222 + 0.864*(g-r)
	#print('bell_03_m_l: '+str(bell_03_m_l))
	#print('M_i: '+str(M_i))
	bell_03_mass = np.log10(bell_03_m_l)-((M_i-4.57)/2.51)
	#print('bell_03_mass: '+str(bell_03_mass))
	bell_03_mass_error = np.sqrt((((0.864*np.log(10)*bell_03_mass)**2)*g_error**2)+(((0.864*np.log(10)*bell_03_mass)**2)*r_error**2)+((((np.log(10)*bell_03_mass/2.51))**2)*M_i_error**2))
	bell_03_mass_error = np.log10(bell_03_mass_error)
	#print('bell_03_mass_error: '+str(bell_03_mass_error))
	#West 2010 masses
	west_10_log_m_l = -0.222+(0.864*(g-r)+np.log10(0.71))
	west_10_mass = west_10_log_m_l-((M_i-4.57)/2.51)
	#huang_west_10_log_m_l = -0.222+(0.864*(huang_g-huang_r)+np.log10(0.71))
	#huang_west_10_mass = huang_west_10_log_m_l-((huang_M_i-4.57)/2.51)
	
	HI_mass = np.array([7.74, 7.75, 7.60, 7.18, 7.41, 7.22, 7.46, 7.66, 7.67, 8.15, 7.43])
	sdss_dwarf_mass = [6.226313, -9.999990, -9.999990, 8.045958, 9.022079, -9.999990, 8.334774, 7.633332, 7.270247, -9.999990, -9.999990]
	FAST_mass = np.array([5.23, 6.30, 5.59, 6.44, 7.89, 4.35, 7.50, 6.49, 5.69, 5.51, 4.54])
	light_to_mass = np.array([ 6.27922169,  6.78568557,  6.67363125,  7.24777424,  8.20700453,  5.071811,  7.94262839,  7.35428572,  7.04113214,  6.59401507,  5.61404239])

	converted_HI_mass_high = (0.02*HI_mass)+9.52 #Catinella et al. 2010 mostly for high mass
	converted_HI_mass_low =  (0.712*HI_mass)+3.117 #low mass relation Shan Huang 2012b
	Jimmy_mass = west_10_mass
	#print('Jimmy_mass: '+str(Jimmy_mass))
	Jimmy_mass_error = np.sqrt(bell_03_mass_error**2+0.2**2)


if (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity')or (plot_to_make == '3d_hi_mass_metallicity') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == '3d_hi_mass_sfr') or (plot_to_make == '3d_hi_metallicity_sfr') or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr')  or (plot_to_make == '4d') or (plot_to_make == 'mass_comparison') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
	#alfalfa_sdss_catalog_file = HOME+"/Astro/Catalogs/alfalfa_sdss_matched.txt"
	alfalfa_sdss_catalog_file = HOME+"/Astro/Catalogs/alfalfa_sdss_averaged.csv"
	alfalfa_sdss_catalog = np.genfromtxt(alfalfa_sdss_catalog_file, dtype=None, names=['AGCNr','hi_mass','distance','oh_p50','lgm_tot_p50','sfr_tot_p50','h_beta_flux','h_beta_flux_err','oiii_5007_flux','oiii_5007_flux_err','h_alpha_flux','h_alpha_flux_err','nii_6584_flux','nii_6584_flux_err','petroR90_r','petroR90Err_r','petroR50_r','petroR50Err_r','petroMag_u','petroMagErr_u','petroMag_g','petroMagErr_g','petroMag_r','petroMagErr_r','petroMag_i','petroMagErr_i','petroMag_z','petroMagErr_z','z','zErr','fiberMag_r','fiberMagErr_r'], delimiter=",", skiprows=1)
	#alfalfa_sdss_photometry_catalog_file = HOME+"/Astro/Catalogs/ALFALFA-SDSS_photometry_and_redshift.csv" 
	#alfalfa_sdss_photometry_catalog = np.genfromtxt(alfalfa_sdss_photometry_catalog_file, dtype=None, names=['dr7objid'])
	
	#Jimmy_lgm_tot_p50 = np.zeros_like(Jimmy_mass)
	#for index_dummy_1, galaxy_name in enumerate(names):
	#	for index_dummy_2, sdss_galaxy in enumerate(alfalfa_sdss_catalog['AGCNr']):
	#		if (str(sdss_galaxy) == str(galaxy_name[3:])):
	#			print(sdss_galaxy)
	#			print(galaxy_name[3:])
	#			Jimmy_lgm_tot_p50[index_dummy_1] = (alfalfa_sdss_catalog['lgm_tot_p50'])[index_dummy_2]
	#print('Jimmy_lgm_tot_p50: '+str(Jimmy_lgm_tot_p50))
	#for index_dummy_3, lgm in enumerate(Jimmy_lgm_tot_p50):
	#	if lgm > 0:
	#		print('My Mass: '+str(lgm)+' sdss total mass: '+str(Jimmy_mass[index_dummy_3]))
	#		print('Percent difference between masses: '+str(round(100*abs(lgm-Jimmy_mass[index_dummy_3])/Jimmy_mass[index_dummy_3]),2))
	
	SDSS_B = alfalfa_sdss_catalog['petroMag_g'] + 0.3130*(alfalfa_sdss_catalog['petroMag_g'] - alfalfa_sdss_catalog['petroMag_r']) + 0.2271
	#SDSS_Dp = (speed_of_light*alfalfa_sdss_catalog['z']/H0)*(1-(((1+q0)/2)*alfalfa_sdss_catalog['z']))
	#SDSS_luminosity_distance = (1+alfalfa_sdss_catalog['z'])*SDSS_Dp
	SDSS_luminosity_distance = np.zeros_like(alfalfa_sdss_catalog['z'])
	for index_z, galaxy_z in enumerate(alfalfa_sdss_catalog['z']):
		#if index_z < 10:
			#print(galaxy_z)
		if galaxy_z > (6000/speed_of_light):
			SDSS_Dp = (speed_of_light*galaxy_z/H0)*(1-(((1+q0)/2)*galaxy_z))
			SDSS_luminosity_distance[index_z] = (1+galaxy_z)*SDSS_Dp
		else:
			SDSS_luminosity_distance[index_z] = (alfalfa_sdss_catalog['distance'])[index_z]
	
	#print('printing the first 10 values of the luminosity distances: '+str(SDSS_luminosity_distance[:10]))
	SDSS_M_B = SDSS_B + 5 - (5*np.log10(SDSS_luminosity_distance*1e6))
	
	SDSS_M_i = alfalfa_sdss_catalog['petroMag_i'] + 5 - (5*np.log10(SDSS_luminosity_distance*1e6))
	SDSS_bell_03_m_l = 0.222 + 0.864*(alfalfa_sdss_catalog['petroMag_g']-alfalfa_sdss_catalog['petroMag_r'])
	SDSS_bell_03_mass = np.log10(SDSS_bell_03_m_l)-((SDSS_M_i-4.57)/2.51)
	SDSS_west_10_log_m_l = -0.222+(0.864*(alfalfa_sdss_catalog['petroMag_g']-alfalfa_sdss_catalog['petroMag_r'])+np.log10(0.71))
	SDSS_west_10_mass = SDSS_west_10_log_m_l-((SDSS_M_i-4.57)/2.51)
	
	SDSS_O3N2 = np.log10((alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])/(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux']))
	SDSS_logOH_PP04_O3N2 = 8.73-(0.32*SDSS_O3N2)
	SDSS_N2 = np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])
	SDSS_logOH_PP04_N2 = 9.37 + (2.03*SDSS_N2) + (1.26*SDSS_N2**2) + (0.32*SDSS_N2**3)
	#SDSS_logOH_PP04_N2 = 8.90 + 0.57*SDSS_N2
	SDSS_logOH_D02_N2 = 9.12 + 0.73*SDSS_N2
	for index, galaxy in enumerate(SDSS_logOH_D02_N2):
		if galaxy == 9.12:
			SDSS_logOH_D02_N2[index] = 0
	#SDSS_mass = alfalfa_sdss_catalog['lgm_tot_p50']
	#SDSS_mass = SDSS_west_10_mass
	
	#h_alpha_flux units: 1e-17 erg/s/cm^2
	SDSS_cm_luminosity_distance = SDSS_luminosity_distance*3.08567758e24 
	SDSS_Ha_flux = alfalfa_sdss_catalog['h_alpha_flux']*1e-17#*angstroms
	SDSS_Luminosity = SDSS_Ha_flux*4*np.pi*SDSS_cm_luminosity_distance**2
	SDSS_balmer_dec = alfalfa_sdss_catalog['h_alpha_flux']/alfalfa_sdss_catalog['h_beta_flux']
	SDSS_E_B_minus_V = (1.97)*np.log10(SDSS_balmer_dec/2.86) #0.950 is from Reddy which assumes that the Cardelli et al. (1989) extinction curve applies to the HII regions #evaluates to 0.43429448190325182
	#print(SDSS_balmer_dec)
	SDSS_H_alpha_extinction = 3.33*SDSS_E_B_minus_V #Reddy, again using the Cardelli stuff.
	#SDSS_H_alpha_extinction = 5.91414*np.log10(SDSS_balmer_dec/2.86)
	SDSS_True_Luminosity = SDSS_Luminosity*10**(SDSS_H_alpha_extinction/2.5)*10**((alfalfa_sdss_catalog['fiberMag_r']-alfalfa_sdss_catalog['petroMag_r'])/2.5)
	aperature_correction = (10**((alfalfa_sdss_catalog['fiberMag_r']-alfalfa_sdss_catalog['petroMag_r'])/2.5))
	aperature_correction = (alfalfa_sdss_catalog['fiberMag_r']-alfalfa_sdss_catalog['petroMag_r'])
	#print('max(aperature_correction)'+str(max(aperature_correction)))
	#print('min(aperature_correction)'+str(min(aperature_correction)))
	#print('median(aperature_correction)'+str(np.median(aperature_correction)))
	#print('mean(aperature_correction)'+str(np.mean(aperature_correction)))
	#print('std(aperature_correction)'+str(np.std((aperature_correction)-np.mean(aperature_correction))))
	#Kennicutt 2012 logCX is 41.27 for H_alpha
	log_Cx = 41.27
	#SFR is in solar masses per year	
	SDSS_sfr = np.log10(SDSS_True_Luminosity)-log_Cx
	for sfr_index, sfr in enumerate(SDSS_sfr):
		if aperature_correction[sfr_index] > 5.0:
			#if (alfalfa_sdss_catalog['petroR90_r'])[sfr_index] > 6.0:
			SDSS_sfr[sfr_index] = np.inf
		if ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 191702) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 220755) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 220837) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 221000) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 221004) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 225852):
			print(str((alfalfa_sdss_catalog['AGCNr'])[sfr_index])+' aperature correction: '+str(aperature_correction[sfr_index]))
		if ( plot_to_make != 'consistency_check'):
			if ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 191702) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 220755) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 220837) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 221000) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 221004) or ((alfalfa_sdss_catalog['AGCNr'])[sfr_index] == 225852):
				SDSS_sfr[sfr_index] = np.inf
			
	#temp_variable = SDSS_sfr[aperature_correction<72]
	
	#sdss_table=np.loadtxt(HOME+"/Astro/sdss_speclines.csv", delimiter=",", skiprows=1, dtype={'names': ('oiii_5007_flux', 'h_beta_flux', 'h_alpha_flux', 'nii_6584_flux', 'lgm_tot_p50', 'lgm_fib_p50'), 'formats': ('f4', 'f4', 'f4', 'f4', 'f4', 'f4')})
	#O3N2 = np.log10((sdss_table['oiii_5007_flux']/sdss_table['h_beta_flux'])/(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux']))
	#sdss_sample_PP04_O3N2 = 8.73-(0.32*O3N2)
	
	#temp_metallicity = (alfalfa_sdss_catalog['PP04_Metallicity'])[np.logical_and(np.isfinite(alfalfa_sdss_catalog['PP04_Metallicity']),np.isfinite(alfalfa_sdss_catalog['lgm_tot_p50']))]
	#temp_mass = (alfalfa_sdss_catalog['lgm_tot_p50'])[np.logical_and(np.isfinite(alfalfa_sdss_catalog['PP04_Metallicity']),np.isfinite(alfalfa_sdss_catalog['lgm_tot_p50']))]
	#cleaned_metallicity = cleaned_metallicity[np.isfinite(cleaned_mass)]
	#cleaned_mass = cleaned_mass[np.isfinite(cleaned_mass)]
	#cleaned_mass = temp_mass[np.logical_and(temp_metallicity>0.0,temp_mass>0.0)]
	#cleaned_metallicity = temp_metallicity[np.logical_and(temp_metallicity>0.0,temp_mass>0.0)]
	#cleaned_mass = cleaned_mass[cleaned_mass>1.0]
	#cleaned_metallicity = cleaned_metallicity[cleaned_mass>1.0]
	#np.set_printoptions(threshold=np.inf)
	#print(cleaned_mass>0)
	
	print('SDSS Sample without any cuts: '+str(len(SDSS_sfr)))
	
	if (calibration == 'N2'):
		agn_cleaned_metallicity = SDSS_logOH_D02_N2[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	if (calibration == 'O3N2'):
		agn_cleaned_metallicity = SDSS_logOH_PP04_O3N2[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	agn_cleaned_mass = SDSS_west_10_mass[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	agn_cleaned_hi_mass = (alfalfa_sdss_catalog['hi_mass'])[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	agn_cleaned_M_B = SDSS_M_B[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	#agn_cleaned_sfr = (alfalfa_sdss_catalog['sfr_tot_p50'])[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]	
	agn_cleaned_sfr = SDSS_sfr[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]	
	agn_cleaned_redshift = (alfalfa_sdss_catalog['z'])[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]  #(alfalfa_sdss_catalog['z'])[alfalfa_sdss_catalog['z'] > 0.01]
	agn_cleaned_h_alpha = (alfalfa_sdss_catalog['h_alpha_flux'])[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	agn_cleaned_h_alpha_err = (alfalfa_sdss_catalog['h_alpha_flux_err'])[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	agn_cleaned_names = (alfalfa_sdss_catalog['AGCNr'])[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	agn_cleaned_mpa_jhu_mass = (alfalfa_sdss_catalog['lgm_tot_p50'])[np.logical_and(np.log10(alfalfa_sdss_catalog['oiii_5007_flux']/alfalfa_sdss_catalog['h_beta_flux'])<1.3+0.61/(np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])-0.05),np.log10(alfalfa_sdss_catalog['nii_6584_flux']/alfalfa_sdss_catalog['h_alpha_flux'])<0.05)]
	print('SDSS Sample after AGN cleaning: '+str(len(agn_cleaned_sfr)))
	
	#finite_metallicity = agn_cleaned_metallicity[np.isfinite(agn_cleaned_mass)]
	#finite_mass = agn_cleaned_mass[np.isfinite(agn_cleaned_mass)]
	#finite_hi_mass = agn_cleaned_hi_mass[np.isfinite(agn_cleaned_mass)]
	#finite_M_B = agn_cleaned_M_B[np.isfinite(agn_cleaned_mass)]
	#finite_sfr = agn_cleaned_sfr[np.isfinite(agn_cleaned_mass)]
	#finite_redshift = agn_cleaned_redshift[np.isfinite(agn_cleaned_mass)]
	#finite_h_alpha = agn_cleaned_h_alpha[np.isfinite(agn_cleaned_mass)]
	#finite_h_alpha_err = agn_cleaned_h_alpha_err[np.isfinite(agn_cleaned_mass)]
	#finite_names = agn_cleaned_names[np.isfinite(agn_cleaned_mass)]
	#print('SDSS Sample after removing infinite and nan masses: '+str(len(finite_sfr)))
	
	
	finite_metallicity = agn_cleaned_metallicity[np.isfinite(agn_cleaned_sfr)]
	finite_mass = agn_cleaned_mass[np.isfinite(agn_cleaned_sfr)]
	finite_hi_mass = agn_cleaned_hi_mass[np.isfinite(agn_cleaned_sfr)]
	finite_M_B = agn_cleaned_M_B[np.isfinite(agn_cleaned_sfr)]
	finite_sfr = agn_cleaned_sfr[np.isfinite(agn_cleaned_sfr)]
	finite_redshift = agn_cleaned_redshift[np.isfinite(agn_cleaned_sfr)]
	finite_h_alpha = agn_cleaned_h_alpha[np.isfinite(agn_cleaned_sfr)]
	finite_h_alpha_err = agn_cleaned_h_alpha_err[np.isfinite(agn_cleaned_sfr)]
	finite_names = agn_cleaned_names[np.isfinite(agn_cleaned_sfr)]
	finite_mpa_jhu_mass = agn_cleaned_mpa_jhu_mass[np.isfinite(agn_cleaned_sfr)]
	print('SDSS Sample after removing too high aperature corrections (infinite and nan sfrs): '+str(len(finite_sfr)))
	
	
	
	#hi_clean_metallicity = finite_metallicity[finite_hi_mass != 10.24]
	#hi_clean_mass = finite_mass[finite_hi_mass != 10.24]
	#hi_clean_hi_mass = finite_hi_mass[finite_hi_mass != 10.24]
	#hi_clean_M_B = finite_M_B[finite_hi_mass != 10.24]
	#hi_clean_sfr = finite_sfr[finite_hi_mass != 10.24]
	
	#print('SDSS Sample after removing that weird HI mass feature: '+str(len(hi_clean_sfr)))
	#print(finite_redshift)
	z_clean_metallicity = finite_metallicity[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_mass = finite_mass[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_hi_mass = finite_hi_mass[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_M_B = finite_M_B[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_sfr = finite_sfr[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_redshift = finite_redshift[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_h_alpha = finite_h_alpha[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_h_alpha_error = finite_h_alpha_err[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_names = finite_names[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	z_clean_mpa_jhu_mass = finite_mpa_jhu_mass[np.logical_and((finite_redshift < 0.05), (finite_redshift > (-6000/speed_of_light)))]
	print('SDSS Sample after selecting the redshift slice: '+str(len(z_clean_redshift)))
	
	#sn_cut_
	sn_cut_SDSS_metallicity = z_clean_metallicity[z_clean_h_alpha/z_clean_h_alpha_error > 3]
	sn_cut_SDSS_mass = z_clean_mass[z_clean_h_alpha/z_clean_h_alpha_error > 3]
	sn_cut_SDSS_hi_mass = z_clean_hi_mass[z_clean_h_alpha/z_clean_h_alpha_error > 3]
	sn_cut_SDSS_M_B = z_clean_M_B[z_clean_h_alpha/z_clean_h_alpha_error > 3]
	sn_cut_SDSS_sfr = z_clean_sfr[z_clean_h_alpha/z_clean_h_alpha_error > 3]
	sn_cut_SDSS_redshift = z_clean_redshift[z_clean_h_alpha/z_clean_h_alpha_error > 3]
	sn_cut_SDSS_names = z_clean_names[z_clean_h_alpha/z_clean_h_alpha_error > 3]
	sn_cut_SDSS_mpa_jhu_mass = z_clean_mpa_jhu_mass[z_clean_h_alpha/z_clean_h_alpha_error > 3]
	print('SDSS Sample after performing h_alpha S/N cut: '+str(len(sn_cut_SDSS_sfr)))
	
	
	#mass_missmatch_cut
	#plt.scatter(mpa_jhu_mass, SDSS_west_10_mass)
	#reject_counter = 0
	#for index, mpa_jhu_galaxy in enumerate(sn_cut_SDSS_mpa_jhu_mass):
	#	if (abs(mpa_jhu_galaxy-sn_cut_SDSS_mass[index]) < 1.0):
	#		plt.scatter(mpa_jhu_galaxy, sn_cut_SDSS_mass[index], color='blue')
	#	else:
	#		plt.scatter(mpa_jhu_galaxy, sn_cut_SDSS_mass[index], color='red')
	#		reject_counter = reject_counter+1
	#print('reject_counter: '+str(reject_counter))
	#plt.xlim(6.0,12.0)
	#plt.ylim(6.0,12.0)
	#plt.show()
	#stop()
	mass_missmatch_cut_SDSS_metallicity = sn_cut_SDSS_metallicity[abs(sn_cut_SDSS_mpa_jhu_mass-sn_cut_SDSS_mass) < 0.5]
	mass_missmatch_cut_SDSS_mass = sn_cut_SDSS_mass[abs(sn_cut_SDSS_mpa_jhu_mass-sn_cut_SDSS_mass) < 0.5]
	mass_missmatch_cut_SDSS_hi_mass = sn_cut_SDSS_hi_mass[abs(sn_cut_SDSS_mpa_jhu_mass-sn_cut_SDSS_mass) < 0.5]
	mass_missmatch_cut_SDSS_M_B = sn_cut_SDSS_M_B[abs(sn_cut_SDSS_mpa_jhu_mass-sn_cut_SDSS_mass) < 0.5]
	mass_missmatch_cut_SDSS_sfr = sn_cut_SDSS_sfr[abs(sn_cut_SDSS_mpa_jhu_mass-sn_cut_SDSS_mass) < 0.5]
	mass_missmatch_cut_SDSS_redshift = sn_cut_SDSS_redshift[abs(sn_cut_SDSS_mpa_jhu_mass-sn_cut_SDSS_mass) < 0.5]
	mass_missmatch_cut_SDSS_names = sn_cut_SDSS_names[abs(sn_cut_SDSS_mpa_jhu_mass-sn_cut_SDSS_mass) < 0.5]
	mass_missmatch_cut_mpa_jhu_mass = sn_cut_SDSS_mpa_jhu_mass[abs(sn_cut_SDSS_mpa_jhu_mass-sn_cut_SDSS_mass) < 0.5]
	print('SDSS Sample after performing mass miss-match cut: '+str(len(mass_missmatch_cut_mpa_jhu_mass)))
	
	#fig = plt.figure()
	#oiii = sdss_table['oiii_5007_flux']
	#nii = sdss_table['nii_6584_flux']
	#halpha = sdss_table['h_alpha_flux']
	#hbeta = sdss_table['h_beta_flux']
	#cleaned_oiii = oiii[np.logical_and(np.log10(sdss_table['oiii_5007_flux']/sdss_table['h_beta_flux'])<1.3+0.61/(np.log10(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux'])-0.05),np.log10(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux'])<0.05)]
	#cleaned_nii = nii[np.logical_and(np.log10(sdss_table['oiii_5007_flux']/sdss_table['h_beta_flux'])<1.3+0.61/(np.log10(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux'])-0.05),np.log10(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux'])<0.05)]
	#cleaned_halpha = halpha[np.logical_and(np.log10(sdss_table['oiii_5007_flux']/sdss_table['h_beta_flux'])<1.3+0.61/(np.log10(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux'])-0.05),np.log10(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux'])<0.05)]
	#cleaned_hbeta = hbeta[np.logical_and(np.log10(sdss_table['oiii_5007_flux']/sdss_table['h_beta_flux'])<1.3+0.61/(np.log10(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux'])-0.05),np.log10(sdss_table['nii_6584_flux']/sdss_table['h_alpha_flux'])<0.05)]
	#bpt_x = np.log10(nii/halpha)
	#bpt_y = np.log10(oiii/hbeta)
	#plt.scatter(bpt_x, bpt_y, color='red')
	#bpt_x = np.log10(cleaned_nii/cleaned_halpha)
	#bpt_y = np.log10(cleaned_oiii/cleaned_hbeta)
	#plt.scatter(bpt_x, bpt_y, color='blue')
	#plt.show()
	SDSS_mass = mass_missmatch_cut_SDSS_mass
	SDSS_metallicity = mass_missmatch_cut_SDSS_metallicity
	SDSS_hi_mass = mass_missmatch_cut_SDSS_hi_mass
	SDSS_M_B = mass_missmatch_cut_SDSS_M_B
	SDSS_sfr = mass_missmatch_cut_SDSS_sfr
	if (dog == 'y'):
		SDSS_sfr = SDSS_sfr+np.log10(7)
	SDSS_redshift = mass_missmatch_cut_SDSS_redshift
	SDSS_names = mass_missmatch_cut_SDSS_names
	mpa_jhu_mass = mass_missmatch_cut_mpa_jhu_mass
	#stop()
	
	#SDSS_mass = finite_mass 
	#SDSS_metallicity = finite_metallicity
	#SDSS_hi_mass = finite_hi_mass
	#SDSS_M_B = finite_M_B
	#SDSS_sfr = finite_sfr
	#SDSS_redshift = finite_redshift
	#np.set_printoptions(threshold='nan')
	#print(SDSS_mass)

#This is a general grid of plots with different switches turned on or off depending on which plots need them
if (plot_to_make == 'h_alpha') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'balmer') or (plot_to_make == 'metallicity') or (plot_to_make == 'bpt') or (plot_to_make == 'gradients') or (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == '3d_hi_mass_metallicity') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == '3d_hi_mass_sfr') or (plot_to_make == '3d_hi_metallicity_sfr')  or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass')  or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '4d') or (plot_to_make == 'sfr_comparison') or (plot_to_make == 'mass_comparison') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
	#loop through each galaxy as defined in the list of names above.
	reorder = [4, 3, 6, 7, 8, 1, 2, 10, 5, 0, 9]
	reorder = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	if (plot_to_make == 'balmer') or (plot_to_make == 'h_alpha') or (plot_to_make == 'metallicity') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr'):
		names = names[reorder]
		xmin = xmin[reorder]
		xmax = xmax[reorder]
		ymin = ymin[reorder]
		ymax = ymax[reorder]
	for index, galaxy in enumerate(names):
		if (plot_to_make == 'h_alpha') or (plot_to_make == 'balmer') or (plot_to_make == 'metallicity'):
			if plotmode == "paper":
				sp1 = fig.add_subplot(3,4,index+1)
				sp1.set_xticks([-10, -5, 0, 5, 10])
			if plotmode == "slides":
				sp1 = fig.add_subplot(3,8,(index*2)+1)
				sp1.set_xticklabels([])
				sp1.set_yticklabels([])
		if (plot_to_make == 'bpt'):
			if plot_type == 'radius':
				sp1 = fig.add_subplot(2,4,index+1)
			if plot_type == 'regular' or plot_type == 'metallicity':
				sp1 = fig.add_subplot(5,4,index+1)
		if (plot_to_make == 'gradients'):
			if plotmode == "paper":
				sp1 = fig.add_subplot(4,3,index+1)
			if plotmode == "slides":
				sp1 = fig.add_subplot(3,4,index+1)
		#2d bins file contains the x & y coordinates, as well as the s/n results
		voronoi_2d_bins_file = HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/voronoi_2d_bins_emission.txt'
		voronoi_2d_bins_table=np.loadtxt(voronoi_2d_bins_file, dtype=float)
		xbin = voronoi_2d_bins_table[:,0]
		ybin = voronoi_2d_bins_table[:,1]
		#binning output file contains the binning information necessary to read in the monte carlo stuff
		if (plot_to_make == 'bpt') or (plot_to_make == 'gradients'):
			voronoi_2d_binning_output_file = HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/voronoi_2d_binning_output_emission.txt'
			voronoi_2d_binning_output_table = np.loadtxt(voronoi_2d_binning_output_file, dtype=float)

		if fit_mode == "gandalf":
			gandalf_file = HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/gandalf.fits'
		if fit_mode == "swift":
			gandalf_file = HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/pink_gandalf.fits'
		#Open the fits file that contains the emission line fitting results
		hdulist = pyfits.open(gandalf_file)
		#In gandalf the emission line fit results are in the 4th extension
		kinematic_measurements = hdulist[4].data
		velocity = kinematic_measurements[:,43]
		flux_OIII2 = kinematic_measurements[:,31]
		flux_Hb = kinematic_measurements[:,21]
		flux_NII2 = kinematic_measurements[:,46]
		flux_Ha = kinematic_measurements[:,41]
		
		#Pull in the integrated results from all spaxels stacked into one spectrum
		one_bin_table = np.genfromtxt(HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/gandalf_table.txt',dtype=None)
		Hb_line = one_bin_table[4]
		stacked_Hb = Hb_line[1]
		OIII2_line = one_bin_table[3]
		stacked_OIII2 = OIII2_line[1]
		Ha_line = one_bin_table[6]
		stacked_Ha = Ha_line[1]
		#print(stacked_Ha)
		NII2_line = one_bin_table[5]
		stacked_NII2 = NII2_line[1]
		
		stacked_O3N2 = np.log10((stacked_OIII2/stacked_Hb)/(stacked_NII2/stacked_Ha))
		#print('Stacked O3N2 values: '+str(stacked_O3N2))
		
		#The balmer decrement is calculated for the balmer plots, and also to make corrections in the h_alpha plots
		if (plot_to_make == 'balmer') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'h_alpha') or (plot_to_make == 'hi_mass_sfr')  or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '4d') or (plot_to_make == 'sfr_comparison') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'consistency_check'):
			stacked_balmer_dec = stacked_Ha/stacked_Hb
			#print('stacked_balmer_dec:' +str(stacked_balmer_dec))
			E_B_minus_V = (1.97)*np.log10(stacked_balmer_dec/2.86)
			#print('E (B-V):'+str(E_B_minus_V))
			if (plot_to_make == 'table_of_everything'):
				E_B_minus_V_array = np.append(E_B_minus_V,E_B_minus_V_array)
		
		#Begin the heavy lifting of the SFR calculations and H alpha plotting
		if (plot_to_make == 'h_alpha') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == '3d_hi_mass_metallicity') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == '3d_hi_mass_sfr') or (plot_to_make == '3d_hi_metallicity_sfr')  or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass')  or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '4d') or (plot_to_make == 'sfr_comparison') or (plot_to_make == 'mass_comparison') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'consistency_check'):
			if (plot_to_make == 'mass_metallicity')  or (plot_to_make == 'luminosity_metallicity'):
				distance = np.append(distance,luminosity_distance[index])
			if (plot_to_make == 'hi_mass_density_sfr'):
				distance = np.append(distance,Dp)
			if (plot_to_make == 'h_alpha') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '4d') or (plot_to_make == 'sfr_comparison') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'consistency_check'):
				#print(luminosity_distance[index])
				#Used for testing the 'snall' masking technique.
				masked_one_bin_table = np.genfromtxt(HOME+'/Astro/reduced/'+galaxy+'pro/all/'+'snall'+'/gandalf_table.txt',dtype=None)
				masked_Hb_line = masked_one_bin_table[4]
				masked_stacked_Hb = masked_Hb_line[1]
				masked_OIII2_line = masked_one_bin_table[3]
				masked_stacked_OIII2 = masked_OIII2_line[1]
				masked_Ha_line = masked_one_bin_table[6]
				masked_stacked_Ha = masked_Ha_line[1]
				masked_NII2_line = masked_one_bin_table[5]
				masked_stacked_NII2 = masked_NII2_line[1]
				cm_luminosity_distance = luminosity_distance[index]*3.08567758e24
				cm_luminosity_distance_error = luminosity_distance_error[index]*3.08567758e24
				angstroms = 7.5*5.2
				Ha_flux = masked_stacked_Ha*1e-16#*angstroms
				#if (index == 1) or (index == 5) or (index == 6) or (index == 8):
				#	#stop()
				#	Ha_flux = 1e-14
				monte_one_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+'snall'+"/monte_carlo_results_one.txt"
				monte_one_table=np.loadtxt(monte_one_file, dtype=float)
				h_alpha_flux_error = np.nanstd(monte_one_table[:,11])
				Ha_flux_error = h_alpha_flux_error*1e-17
				h_beta_flux_error = np.nanstd(monte_one_table[:,3])
				Luminosity = Ha_flux*4*np.pi*cm_luminosity_distance**2
				H_alpha_extinction = 3.33*E_B_minus_V #Correct the h_alpha measurement values using the balmer decrement calculations above
				#print(H_alpha_extinction)
				#if H_alpha_extinction > 1000:
				#	H_alpha_extinction = 1.25
				True_Luminosity = Luminosity*10**(H_alpha_extinction/2.5)
				#print(cm_luminosity_distance) #This matches
				#print(Ha_flux)
				#print(True_Luminosity)
				#Kennicutt 2012 logCX is 41.27 for H_alpha
				log_Cx = 41.27
				#SFR is in solar masses per year	
				SFR = np.log10(True_Luminosity)-log_Cx
				if (dog == 'y'):
					SFR = SFR+np.log10(7)
				#if SFR > 5:
				#	print('H_alpha_extinction: '+str(H_alpha_extinction))
				#	print('True_Luminosity: '+str(True_Luminosity))
				#	print('SFR: '+str(SFR))
				#	stop()
				#print(SFR)
				#print('SFR: '+str(True_Luminosity/10**log_Cx))
				av_error = np.sqrt(((2.29374126/Ha_flux)*(Ha_flux_error))**2+((2.29374126*(stacked_Ha*1e-16/(stacked_Hb*1e-16)**2))*(h_beta_flux_error*1e-16))**2)
				#print('av_error: '+str(av_error))
				#print('Ha_flux_error: '+str(Ha_flux_error))
				#print('cm_luminosity_distance_error: '+str(cm_luminosity_distance_error))
				#SFR_error = (np.sqrt(((4*np.pi*cm_luminosity_distance**2*10**3.33*E_B_minus_V/2.5)*(Ha_flux_error))**2+((8*Ha_flux*np.pi*cm_luminosity_distance*10**(3.33*E_B_minus_V/2.5))*(cm_luminosity_distance_error))**2+((4*np.pi*Ha_flux*np.log(10)*0.4*cm_luminosity_distance**2*10**(3.33*E_B_minus_V/2.5))*(av_error))**2))
				SFR_error = 1e9*(np.sqrt(((4*np.pi*10**3.33*E_B_minus_V/2.5)*(Ha_flux_error))**2+((8*Ha_flux*np.pi*10**(3.33*E_B_minus_V/2.5))*(0.1))**2+((4*np.pi*Ha_flux*np.log(10)*0.4*10**(3.33*E_B_minus_V/2.5))*(av_error))**2))
				#print('sfr error: '+str(SFR_error))
				percent_error = SFR_error/(True_Luminosity/10**log_Cx)
				#print('percent error: '+str(percent_error))
				SFR_error = abs(percent_error*SFR)
				if (plot_to_make == 'hi_mass_sfr') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == '4d') or (plot_to_make == 'sfr_comparison') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr') or (plot_to_make == 'consistency_check'):
					Jimmy_sfr = np.append(Jimmy_sfr,SFR)
					Jimmy_sfr_error = np.append(Jimmy_sfr_error, SFR_error)
					#Jimmy_sfr = np.array([-3.22, -3.45, -2.60, -2.02, -2.50, -5.11, -2.37, -2.84, -1.33, -3.03, -4.46])
				if (plot_to_make == 'h_alpha'):
					SFR98 = np.log10((7.9e-42)*True_Luminosity) #units ergs/s
					if plotmode == "slides":
						sp1 = fig.add_subplot(3,8,(index*2)+2)
					colorbar_min = 37.0
					colorbar_max = 40.0
					colormap = 'winter_r'
					masked_voronoi_2d_bins_file = HOME+'/Astro/reduced/'+galaxy+'pro/all/'+'snall'+'/voronoi_2d_bins_emission.txt'
					masked_voronoi_2d_bins_table=np.loadtxt(masked_voronoi_2d_bins_file, dtype=float)
					masked_xbin = masked_voronoi_2d_bins_table[:,0]
					masked_ybin = masked_voronoi_2d_bins_table[:,1]
					cmap = plt.cm.winter_r
					cmap.set_under('0.75')
					#x_nudge = 10
					#y_nudge = 1
					fig2 = display_pixels(masked_xbin+(x_nudge[index]*0.66), masked_ybin+(y_nudge[index]*0.66), masked_xbin*0.0, cmap)
					#cmap = plt.cm.winter_r
					#cmap.set_under('black')
					fig2 = display_pixels(xbin, ybin, np.log10(flux_Ha*1e-16*angstroms*4*np.pi*cm_luminosity_distance**2), cmap)
			
		#Begin the heavy lifting of the balmer dec plotting
		if (plot_to_make == 'balmer'):
			balmer_dec = flux_Ha/flux_Hb
			colorbar_min = 0.86
			colorbar_max = 4.86
			colormap = 'Reds'
			fig2 = display_pixels(xbin, ybin, balmer_dec, colormap)
			plt.text(xmin[index]+1, ymin[index]+2,'E(B-V) = '+str(round(E_B_minus_V,3)), fontsize=14)
		
		if (plot_to_make == 'mass_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == '3d_hi_mass_metallicity') or (plot_to_make == 'metallicity_sfr') or (plot_to_make == 'luminosity_metallicity') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '4d') or (plot_to_make == 'metallicity') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'fmr_lzr_sfr')  or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
			#line_of_interest = one_bin_table[0] #0 is for PP04_O3N2, 1 is for PP04_N2, and 2 is for D02_N2
			#stacked_PP04_O3N2 = line_of_interest[0]	
			#line_of_interest = one_bin_table[1] #0 is for PP04_O3N2, 1 is for PP04_N2, and 2 is for D02_N2
			#stacked_PP04_N2 = line_of_interest[1]	
			#line_of_interest = one_bin_table[2] #0 is for PP04_O3N2, 1 is for PP04_N2, and 2 is for D02_N2
			#stacked_D02_N2 = line_of_interest[1]	
			#stacked_metallicity = np.append(stacked_metallicity, stacked_D02_N2)
			Hb_line = one_bin_table[4]
			stacked_Hb = Hb_line[1]
			OIII2_line = one_bin_table[3]
			stacked_OIII2 = OIII2_line[1]
			Ha_line = one_bin_table[6]
			stacked_Ha = Ha_line[1]
			NII2_line = one_bin_table[5]
			stacked_NII2 = NII2_line[1]
			Jimmy_N2 = np.log10(stacked_NII2/stacked_Ha)
			Jimmy_O3N2 = np.log10((stacked_OIII2/stacked_Hb)/(stacked_NII2/stacked_Ha))
			stacked_PP04_O3N2 = 8.73-(0.32*Jimmy_O3N2)
			stacked_PP04_N2 = 9.37 + (2.03*Jimmy_N2) + (1.26*Jimmy_N2**2) + (0.32*Jimmy_N2**3)
			stacked_PP04_N2_linear = 8.90 + 0.57*Jimmy_N2
			stacked_D02_N2 = 9.12 + 0.73*Jimmy_N2
			if (calibration == 'N2'):
				Jimmy_metallicity = np.append(Jimmy_metallicity, stacked_D02_N2)
			if (calibration == 'O3N2'):
				Jimmy_metallicity = np.append(Jimmy_metallicity, stacked_PP04_O3N2)
				
			monte_one_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/monte_carlo_results_one.txt"
			monte_one_table=np.loadtxt(monte_one_file, dtype=float)
			h_beta_flux_error = np.nanstd(monte_one_table[:,3])
			oiii_5007_flux_error = np.nanstd(monte_one_table[:,7])
			h_alpha_flux_error = np.nanstd(monte_one_table[:,11])
			nii_6584_flux_error = np.nanstd(monte_one_table[:,15])
			N2_error = np.sqrt(((1/(stacked_NII2*np.log(10)))**2)*((nii_6584_flux_error)**2)+((1/(stacked_Ha*np.log(10)))**2)*((h_alpha_flux_error)**2))
			D02_N2_error = 0.73*N2_error
			PP04_O3N2_error = 0.32*np.sqrt(((1/(stacked_OIII2*np.log(10)))**2)*((oiii_5007_flux_error)**2)+((1/(stacked_Hb*np.log(10)))**2)*((h_beta_flux_error)**2)+((1/(stacked_NII2*np.log(10)))**2)*((nii_6584_flux_error)**2)+((1/(stacked_Ha*np.log(10)))**2)*((h_alpha_flux_error)**2))
			if (calibration == 'N2'):
				Jimmy_metallicity_error = np.append(Jimmy_metallicity_error, D02_N2_error)
			if (calibration == 'O3N2'):
				Jimmy_metallicity_error = np.append(Jimmy_metallicity_error, PP04_O3N2_error)
			one_gandalf_file = HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/pink_gandalf_one.fits'
			#Open the fits file that contains the emission line fitting results
			one_hdulist = pyfits.open(one_gandalf_file)
			#In gandalf the emission line fit results are in the 4th extension
			one_kinematic_measurements = one_hdulist[4].data
			one_flux_Hb = one_kinematic_measurements[:,21]
			one_flux_OIII1 = one_kinematic_measurements[:,26]
			one_flux_OIII2 = one_kinematic_measurements[:,31]
			one_flux_NII1 = one_kinematic_measurements[:,36]
			one_flux_Ha = one_kinematic_measurements[:,41]
			one_flux_NII2 = one_kinematic_measurements[:,46]
			#print(one_flux_Ha)
			#print(galaxy+' & '+str(round(one_flux_Hb,1))+' $\pm$ '+str(round(h_beta_flux_error,1))+' & '+str(round(one_flux_OIII1,1))+' $\pm$ '+'str(round(oiii_4959_flux_error,1))'+' & '+str(round(stacked_OIII2,1))+' $\pm$ '+str(round(oiii_5007_flux_error,1))+' & '+str(round(one_flux_NII1,1))+' $\pm$ '+'str(round(nii_XXXX_error,1))'+' & '+str(round(stacked_Ha,1))+' $\pm$ '+str(round(h_alpha_flux_error,1))+' & '+str(round(stacked_NII2,1))+' $\pm$ '+str(round(nii_6584_flux_error,1)))
			print(galaxy+' & '+str(round(stacked_Hb,1))+' $\pm$ '+str(round(h_beta_flux_error,1))+' & '+str(round(stacked_OIII2,1))+' $\pm$ '+str(round(oiii_5007_flux_error,1))+' & '+str(round(stacked_Ha,1))+' $\pm$ '+str(round(h_alpha_flux_error,1))+' & '+str(round(stacked_NII2,1))+' $\pm$ '+str(round(nii_6584_flux_error,1))+' \\\\')
		
		#Begin the heavy lifting of the oxygen abundance plots
		if (plot_to_make == 'metallicity'):
			O3N2 = np.log10((flux_OIII2/flux_Hb)/(flux_NII2/flux_Ha))
			#print('Individual O3N2 values: '+str(O3N2))
			logOH_PP04_O3N2 = 8.73-(0.32*O3N2)
			N2 = np.log10(flux_NII2/flux_Ha)
			logOH_PP04_N2 = 9.37 + (2.03*N2) + (1.26*N2**2) + (0.32*N2**3) # FROM 
			logOH_PP04_N2_linear = 8.90 + 0.57*N2
			logOH_D02_N2 = 9.12 + 0.73*N2
			#logOH_PP04_O3N2 = 9.12+(0.73*N2)
			colorbar_min = 7.75
			colorbar_max = 8.75
			#colormap = 'winter' #This is blue to green
			colormap ='cool'
			if (calibration == 'N2'):
				fig2 = display_pixels(xbin, ybin, logOH_D02_N2, colormap)
			if (calibration == 'O3N2'):
				fig2 = display_pixels(xbin, ybin, logOH_PP04_O3N2, colormap)
			line_of_interest = one_bin_table[1] #0 is for PP04_O3N2, 1 is for PP04_N2, and 2 is for D02_N2
			#stacked_PP04_O3N2 = line_of_interest[1]
			plt.text(xmin[index]+1, ymin[index]+1.5,'12+log(O/H) = '+str(round(Jimmy_metallicity[index],3)), fontsize=14)
			#Temporary cheat
			temp_jimmy_mass = np.array([ 6.35926215,  7.08952765,  6.48030926,  6.88741696,  7.87090871,  6.70964029,  7.47002171,  7.08587284,  6.69903967,  5.81322981,  5.69995073])
			temp_jimmy_mass = temp_jimmy_mass[reorder]
			plt.text(xmin[index]+1, ymin[index]+3,r'log(M$_*$) = '+str(round(temp_jimmy_mass[index],2)), fontsize=14)
			
		
		#Begin the heavy lifting of the BPT plots
		if (plot_to_make == 'bpt'):
			if (plot_type == 'metallicity') or (plot_to_make == 'bpt'):
				O3N2 = np.log10((flux_OIII2/flux_Hb)/(flux_NII2/flux_Ha))
				N2 = np.log10((flux_NII2/flux_Ha))
				#print(N2)
				logOH_PP04_O3N2 = 8.73-(0.32*O3N2)
			bpt_y = np.log10(flux_OIII2/flux_Hb)
			bpt_x = np.log10(flux_NII2/flux_Ha)
			bpt_y_error = bpt_y*0.0 #easier to just copy the dimensions of these above then bother with append
			bpt_x_error = bpt_y*0.0
			monte_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/monte_carlo_results.txt"
			monte_table=np.loadtxt(monte_file, dtype=float)
			for sub_index,vbins_line in enumerate(voronoi_2d_binning_output_table):
				temp_monte_array = np.zeros([1,20])
				temp_monte_array[0] = monte_table[0]
				counter = 0
				for monte_line in monte_table:
					if [vbins_line[2],vbins_line[3]] == [monte_line[1],monte_line[2]]:
						if counter == 0:
							temp_monte_array[0] = monte_line
						if counter > 0:
							#print(monte_line)
							temp_monte_array = np.append(temp_monte_array, [monte_line], axis=0)
						counter = counter+1
				h_beta_flux_error = np.nanstd(temp_monte_array[:,3])
				oiii_5007_flux_error = np.nanstd(temp_monte_array[:,7])
				h_alpha_flux_error = np.nanstd(temp_monte_array[:,11])
				nii_6584_flux_error = np.nanstd(temp_monte_array[:,15])
				bpt_y_error[sub_index] = np.sqrt(((1/(flux_OIII2[sub_index]*np.log(10)))**2)*((oiii_5007_flux_error)**2)+((1/(flux_Hb[sub_index]*np.log(10)))**2)*((h_beta_flux_error)**2))
				bpt_x_error[sub_index] = np.sqrt(((1/(flux_NII2[sub_index]*np.log(10)))**2)*((nii_6584_flux_error)**2)+((1/(flux_Ha[sub_index]*np.log(10)))**2)*((h_alpha_flux_error)**2))
			print("Number of monte carlo iterations: "+str(len(temp_monte_array)))
			stacked_bpt_y = np.log10(stacked_OIII2/stacked_Hb)
			stacked_bpt_x = np.log10(stacked_NII2/stacked_Ha)
			print('stacked_bpt_x: '+str(stacked_bpt_x))
			print('stacked_bpt_y: '+str(stacked_bpt_y))
			if plot_type == 'regular':
				plt.scatter(bpt_x,bpt_y, color='blue', label="Dwarf Spaxels", alpha=0.5)
				plt.errorbar(bpt_x,bpt_y,yerr=bpt_y_error, xerr=bpt_x_error, linestyle="None", alpha=0.5)
				plt.scatter(stacked_bpt_x,stacked_bpt_y, color="red", label="Median Dwarf", s=75,zorder=300)
			#if index != 3:
			plt.text(-1.8, -0.8, "Starburst")
			plt.text(0.0, 1.2, "AGN")
			cutoff_x = np.log10(np.array(range(160))*0.01)
			cutoff_y = (0.61/(cutoff_x-0.47))+1.19
			plt.plot(cutoff_x,cutoff_y, color="black")
			new_cutoff_x = np.log10(np.array(range(100))*0.01)
			new_cutoff_y = (0.61/(new_cutoff_x-0.05))+1.3
			plt.plot(new_cutoff_x,new_cutoff_y, color="black")
			plt.ylim(-1.0,2.0)
			plt.xlim(-2.0,1.0)
			if index != 0:
				sp1.set_yticklabels([])
			if index == 4:
				sp1.set_yticklabels([0.0,-0.5,0.0,0.5,1.0,1.5,' '])
				plt.ylabel(r'               LOG([OIII]/H$\beta$)',fontsize=fontsize)
			if index == 0:
				sp1.set_yticklabels([' ',-0.5,0.0,0.5,1.0,1.5,2.0])
			if index < 4:
				sp1.set_xticklabels([])
			if index > 3:
				sp1.set_xticklabels([' ',-1.5,' ',-0.5,' ',0.5,' '])
			if index == 5:
				plt.xlabel(r'                   LOG([NII]/H$\alpha$)',fontsize=fontsize, labelpad=10)
			#if index == 3:
			#	plt.legend(loc='lower right')
			plt.tick_params(axis='y', labelsize=16)
			plt.tick_params(axis='x', labelsize=16)
			plt.text(-1.8, 1.6, galaxy)
			if (plot_type == 'radius') or (plot_type == 'metallicity'):
				plt.colorbar()
		
		#Begin the heavy lifting for the gradients plots
		if (plot_to_make == 'gradients'):
			O3N2 = np.log10((flux_OIII2/flux_Hb)/(flux_NII2/flux_Ha))
			logOH_PP04_O3N2 = np.array(8.73-(0.32*O3N2))
			logOH_PP04_O3N2_error = logOH_PP04_O3N2*0.0 #easier to just make this blank array and fill it in
			radius = np.sqrt(((xbin-np.median(xbin)))*(xbin-np.median(xbin))+((ybin-np.median(ybin))*(ybin-np.median(ybin))))
			monte_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/monte_carlo_results.txt"
			monte_table=np.loadtxt(monte_file, dtype=float)
			for sub_index,vbins_line in enumerate(voronoi_2d_binning_output_table):
				temp_monte_array = np.zeros([1,20])
				temp_monte_array[0] = monte_table[0]
				counter = 0
				for monte_line in monte_table:
					if [vbins_line[2],vbins_line[3]] == [monte_line[1],monte_line[2]]:
						if counter == 0:
							temp_monte_array[0] = monte_line
						if counter > 0:
							#print(monte_line)
							temp_monte_array = np.append(temp_monte_array, [monte_line], axis=0)
						counter = counter+1
				h_beta_flux_error = np.nanstd(temp_monte_array[:,3])
				oiii_5007_flux_error = np.nanstd(temp_monte_array[:,7])
				h_alpha_flux_error = np.nanstd(temp_monte_array[:,11])
				nii_6584_flux_error = np.nanstd(temp_monte_array[:,15])
				logOH_PP04_O3N2_error[sub_index] = 0.32*np.sqrt(((1/(flux_OIII2[sub_index]*np.log(10)))**2)*((oiii_5007_flux_error)**2)+((1/(flux_Hb[sub_index]*np.log(10)))**2)*((h_beta_flux_error)**2)+((1/(flux_NII2[sub_index]*np.log(10)))**2)*((nii_6584_flux_error)**2)+((1/(flux_Ha[sub_index]*np.log(10)))**2)*((h_alpha_flux_error)**2))
			print("Number of monte carlo iterations: "+str(len(temp_monte_array))+" in "+galaxy)
			plt.scatter(radius,logOH_PP04_O3N2)
			plt.errorbar(radius,logOH_PP04_O3N2,yerr=logOH_PP04_O3N2_error, linestyle="None")
			plt.ylim(7.5,9.5)
			plt.xlim(0.0,4.0)
			if plotmode == "paper":
				if index % 3 != 0:
					sp1.set_yticklabels([])
				if index < 8:
					sp1.set_xticklabels([])
				if index == 6:
					plt.ylabel('                          12+log(O/H)',fontsize=fontsize, labelpad=20)
				if index == 9:
					plt.xlabel('                                                   R/R$_e$',fontsize=fontsize, labelpad=20)
			if plotmode == "slides":
				sp1.set_yticklabels([])
				if (index == 0) or (index == 3):
					sp1.set_yticklabels([' ',8.0, 8.5, 9.0, 9.5])
				if index < 4:
					sp1.set_xticklabels([])
				if index > 3:
					sp1.set_xticklabels([' ',0.5,1.0,1.5,2.0,2.5,3.0,3.5,' '])
				if index == 4:
					plt.ylabel('                          12+log(O/H)',fontsize=fontsize, labelpad=20)
					sp1.set_yticklabels([7.5,8.0, 8.5, 9.0, ' '])
					sp1.set_xticklabels([' ',0.5,1.0,1.5,2.0,2.5,3.0,3.5,' '])
				if index == 5:
					plt.xlabel('                                                       R/R$_e$',fontsize=fontsize, labelpad=20)
				if index == 7:
					sp1.set_xticklabels([' ',0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0])
			
			sigma = logOH_PP04_O3N2_error #(logOH_PP04_O3N2*0.0)+1.0
			x0    = np.array([np.median(logOH_PP04_O3N2), 0.0])
			try:
				fit_results = optimization.curve_fit(linear_fit, radius, logOH_PP04_O3N2, x0, logOH_PP04_O3N2_error)
				fit_error = np.sqrt(np.diag(fit_results[1]))
				slope_error = fit_error[1]
				fit_results = fit_results[0]
			except:
				fit_results = np.array([8.2, 666.0])
				slope_error = 666
			x_axis = np.array(range(100))*0.1
			line = linear_fit(x_axis, *fit_results)
			plt.plot(x_axis,line, color='black', label="LSQ Fit")
			plt.text(0.1, 7.6, "Slope: "+str(round(fit_results[1],3))+" +\- "+str(round(slope_error,3)), fontsize=fontsize)
			plt.text(0.1, 9.2, galaxy, fontsize=fontsize)
			plt.tick_params(axis='y', labelsize=16)
			plt.tick_params(axis='x', labelsize=16)
		
			
			
		if (plot_to_make == 'h_alpha') or (plot_to_make == 'balmer') or (plot_to_make == 'metallicity'):
			plt.ylim(ymin[index],ymax[index])
			plt.xlim(xmin[index],xmax[index])
			plt.text(xmin[index]+1, ymax[index]-3,galaxy, fontsize=14)
		if (plot_to_make == 'h_alpha'):
			plt.text(xmin[index]+1, ymin[index]+1.5,r'log(SFR): '+str(round(SFR,2))+' M$_\odot$ yr$^{-1}$', fontsize=14) #+'$\pm$'+str(round(SFR_error,2))
			#plt.text(xmin[index]+1, ymin[index]+2,r'H$\alpha$ flux = '+str(round(stacked_Ha,3)), fontsize=14)
		if (plot_to_make == 'h_alpha') or (plot_to_make == 'metallicity') or (plot_to_make == 'balmer'):
			if plotmode == "paper":
				if index == 0:
					#sp1.set_yticklabels([-10,-5,0,5,10,15])
					sp1.set_yticklabels([-15,-10,-5,0,5,10])
					sp1.yaxis.set_ticks_position('left')
					sp1.xaxis.set_ticks_position('top')
					sp1.set_xticklabels([])
					plt.tick_params(axis='y', labelsize=20)
				if (index == 1) or (index == 2):
					sp1.set_yticks([])
					sp1.xaxis.set_ticks_position('top')
					plt.tick_params(axis='y', labelsize=20)
					sp1.set_xticklabels([])
				if index ==3:
					sp1.yaxis.set_ticks_position('right')
					sp1.xaxis.set_ticks_position('top')
					sp1.set_xticklabels([])
					sp1.set_yticklabels([])
				if index == 4:
					#sp1.set_yticklabels([-10,-5,0,5,10])
					sp1.set_yticklabels([-15,-10,-5,0,5])
					sp1.set_xticks([])
					sp1.yaxis.set_ticks_position('left')
					plt.tick_params(axis='y', labelsize=20)
					plt.ylabel('arcsec',fontsize=fontsize)
					sp1.set_xticklabels([])
				if (index == 5) or (index == 6):
					sp1.set_xticks([])
					sp1.set_yticks([])
				if index == 7:
					sp1.yaxis.set_ticks_position('right')
					sp1.xaxis.set_ticks_position('bottom')
					sp1.set_xticklabels([])
					sp1.set_yticklabels([])
				if index == 8:
					#sp1.set_yticklabels([-15,-10,-5,0,5])
					sp1.set_yticklabels([-15,-10,-5,0,5,10])
					#sp1.set_xticklabels([-15,-10,-5,0,5])
					sp1.set_xticklabels([-5,0,5,10,15])
					sp1.yaxis.set_ticks_position('left')
					sp1.xaxis.set_ticks_position('bottom')
					plt.tick_params(axis='y', labelsize=20)
					plt.tick_params(axis='x', labelsize=20)
				#else:
				#	sp1.set_yticklabels([])
				if index == 9:
					sp1.set_xticklabels([-10,-5,0,5,10])
					plt.tick_params(axis='x', labelsize=20)
					plt.xlabel('arcsec',fontsize=fontsize)
					sp1.xaxis.set_ticks_position('bottom')
					sp1.set_yticks([])
					plt.tick_params(axis='y', labelsize=20)
				if index == 10:
					sp1.set_xticklabels([-10,-5,0,5,10])
					plt.tick_params(axis='x', labelsize=20)
					sp1.xaxis.set_ticks_position('bottom')
					sp1.set_yticks([])
				#else:
					
			if plotmode == "slides":
				sp1.set_xticklabels([])
				sp1.set_yticklabels([])
			
	if (plot_to_make == 'h_alpha'):
		#plt.text(14, -2,r"H$\alpha$ flux", fontsize=fontsize)
		#plt.text(12.5, -6,r"$10^{-16}$ erg", fontsize=15)
		#plt.text(13, -8,r"cm$^{2}$ s $\AA$", fontsize=15)
		plt.text(12, -2,r"H$\alpha$ Luminosity", fontsize=22)
		plt.text(12, -6,r"[log(erg s$^{-1}$)]", fontsize=22)
	if (plot_to_make == 'balmer'):
		plt.text(14, -2,r"Balmer", fontsize=fontsize)
		plt.text(12.5, -6,r"Decrement", fontsize=fontsize)
	if (plot_to_make == 'metallicity'):
		plt.text(12, 0,"12+log(O/H)", fontsize=fontsize)
		if (calibration == 'N2'):
			plt.text(12, -2, "(N2)", fontsize=fontsize)
		if (calibration == 'O3N2'):
			plt.text(12, -2, "(O3N2)", fontsize=fontsize)
	if (plot_to_make == 'h_alpha') or (plot_to_make == 'balmer') or (plot_to_make == 'metallicity'):
		if plotmode == "paper":
			sp1 = fig.add_subplot(4,4,12)
		if plotmode == "slides":
			sp1 = fig.add_subplot(2,1,1)
		f = display_pixels([1], [1], [float("infinity")], colormap)
		plt.axis('off')
		if plotmode == "paper":
			cbaxes = fig.add_axes([0.925, 0.08, 0.02, 0.29])
		if plotmode == "slides":
			cbaxes = fig.add_axes([0.90, 0.35, 0.015, 0.6])
		cb = plt.colorbar(f, cax=cbaxes, extend='min')	
		cb.ax.tick_params(labelsize=fontsize-9)
		fig.subplots_adjust(left=0.08)
		fig.subplots_adjust(bottom=0.08)	
		fig.subplots_adjust(right=0.98)
		fig.subplots_adjust(top=0.98)
		fig.subplots_adjust(wspace=0.00)
		fig.subplots_adjust(hspace=0.00)
		if plotmode == "slides":	
			fig.subplots_adjust(left=0.01)
			fig.subplots_adjust(right=0.87)
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
		#plt.show()
	if (plot_to_make == 'gradients'):
		plt.legend(loc='upper right')
		fig.subplots_adjust(left=0.12)
		fig.subplots_adjust(bottom=0.11)
		fig.subplots_adjust(right=0.98)
		fig.subplots_adjust(top=0.98)
		fig.subplots_adjust(wspace=0.11)
		fig.subplots_adjust(hspace=0.15)
		if plotmode == "slides":
			fig.subplots_adjust(wspace=0.0)
			fig.subplots_adjust(hspace=0.0)
			fig.subplots_adjust(left=0.05)
			fig.subplots_adjust(bottom=0.18)
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
		#plt.show()

if (plot_to_make == 'bpt'):
	if plot_type == 'regular' or plot_type == 'metallicity':
		sp1 = fig.add_subplot(5,4,12)
		plt.ylim(-1.0,2.0)
		plt.xlim(-2.0,1.0)
		sp1.set_yticklabels([])
		sp1.set_xticklabels([])
		plt.scatter([0],[-0.5], color='blue', label="Dwarf Spaxels", alpha=0.5)
		plt.scatter([0],[-0.5], color="red", label="Median Dwarf", s=75)
		plt.legend(loc='lower right')
		table=np.loadtxt(HOME+"/Astro/sdss_speclines.csv", delimiter=",", skiprows=1, dtype={'names': ('oiii_5007_flux', 'h_beta_flux', 'h_alpha_flux', 'nii_6584_flux', 'lgm_tot_p50', 'lgm_fib_p50'), 'formats': ('f4', 'f4', 'f4', 'f4', 'f4', 'f4')})
		sdss_bpt_y = np.log10(table['oiii_5007_flux']/table['h_beta_flux'])
		sdss_bpt_x = np.log10(table['nii_6584_flux']/table['h_alpha_flux'])
		sdss_plot = fig.add_subplot(3,1,3)
		plt.ylim(-1.0,2.0)
		plt.xlim(-2.0,1.0)
		plt.plot(cutoff_x,cutoff_y, color="black")
		plt.plot(new_cutoff_x,new_cutoff_y, color="black")
		if plot_type == 'regular':
			plt.scatter(sdss_bpt_x,sdss_bpt_y, color="gray", label="ALFALFA/SDSS", alpha=sdss_alpha_value, rasterized=True)
		if plot_type == 'metallicity':
			sdss_O3N2 = np.log10((table['oiii_5007_flux']/table['h_beta_flux'])/(table['nii_6584_flux']/table['h_alpha_flux']))
			sdss_logOH_PP04_O3N2 = 8.73-(0.32*sdss_O3N2)
			plt.scatter(sdss_bpt_x,sdss_bpt_y, label="ALFALFA/SDSS", c=sdss_logOH_PP04_O3N2, rasterized=True)
			plt.colorbar()
		#plt.scatter(stacked_bpt_x,stacked_bpt_y, color="red", label="Dwarf Galaxy Medians", s=75)
		plt.legend(loc='upper right')
		plt.xlabel(r'LOG([NII]/H$\alpha$)',fontsize=fontsize, labelpad=10)
		plt.ylabel(r'LOG([OIII]/H$\beta$)',fontsize=fontsize)
		plt.tick_params(axis='y', labelsize=16)
		plt.tick_params(axis='x', labelsize=16)
	fig.subplots_adjust(left=0.13)
	fig.subplots_adjust(bottom=0.10)
	if (plot_type == 'radius'):
		fig.subplots_adjust(bottom=bottom_margin)
	fig.subplots_adjust(right=0.98)
	fig.subplots_adjust(top=0.98)
	fig.subplots_adjust(wspace=0.00)
	fig.subplots_adjust(hspace=0.00)
	if plot_type == 'metallicity':
		fig.subplots_adjust(hspace=0.10)
		fig.subplots_adjust(wspace=0.06)
		fig.subplots_adjust(left=0.10)
		fig.subplots_adjust(bottom=0.06)
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#plt.show()
	
if (plot_to_make == 'distance_slope'):
	from scipy.stats import spearmanr
	distance = [0.388, 0.164, 0.212, 0.141, 0.483, 0.415, 0.325, 0.715]
	distance_error = [0.186, 0.075, 0.101, 0.051, 0.148, 0.149, 0.018, 0.189]
	slope = [0.383, 0.257, 0.098, 0.088, -0.031, -0.093, -0.061, -0.139]
	slope_error = [0.048, 0.049, 0.102, 0.151, 0.032, 0.026, 0.016, 0.038]
	plt.scatter(slope, distance)
	plt.errorbar(slope, distance, xerr=slope_error, yerr=distance_error, linestyle="None")
	plt.xlabel("Gradient Slope")
	plt.ylabel("Distance from AGN division line")
	spearmint_test = spearmanr(slope,distance)
	print('Spearmen rank coefficient: '+str(spearmint_test))
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#plt.show()

	
if (plot_to_make == 'mass_metallicity'):
	fig, ax = plt.subplots()
	#logmstar = 7.0+(np.arange(400)*0.01)
	#PP04_line = 32.1488-(8.51258*logmstar)+(0.976384*logmstar*logmstar)-(0.0359763*logmstar*logmstar*logmstar)
	#Mannucci_y = 8.96+0.31*(Mannucci_x-10)-0.23*(Mannucci_x-10)**2-0.017*(Mannucci_x-10)**3+0.046*(Mannucci_x-10)**4
	#ax.plot(Mannucci_x, Mannucci_y, color='red', linewidth=5.0, label='Mannucci Fit')
	#fitting_SDSS_mass = SDSS_mass[SDSS_mass>0.75]
	#fitting_SDSS_metallicity = SDSS_metallicity[SDSS_mass>0.75]
	#fitting_SDSS_mass = fitting_SDSS_mass[fitting_SDSS_metallicity<90.1]
	#fitting_SDSS_metallicity = fitting_SDSS_metallicity[fitting_SDSS_metallicity<90.1]
	#fitting_SDSS_mass = fitting_SDSS_mass[fitting_SDSS_metallicity>0.6]
	#fitting_SDSS_metallicity = fitting_SDSS_metallicity[fitting_SDSS_metallicity>0.6]
	#final_fitting_SDSS_mass = fitting_SDSS_mass[fitting_SDSS_mass<110.0]
	#final_fitting_SDSS_metallicity = fitting_SDSS_metallicity[fitting_SDSS_mass<110.0]
	#fitting_SDSS_mass = final_fitting_SDSS_mass
	#fitting_SDSS_metallicity = final_fitting_SDSS_metallicity
	fitting_SDSS_mass = SDSS_mass[SDSS_metallicity>0.75]
	fitting_SDSS_metallicity = SDSS_metallicity[SDSS_metallicity>0.75]
	#fitting_SDSS_mass = np.append(fitting_SDSS_mass, Lee_mass)
	fitting_SDSS_mass = np.append(fitting_SDSS_mass, Berg_mass)
	fitting_SDSS_mass = np.append(fitting_SDSS_mass, James_mass)
	fitting_SDSS_mass = np.append(fitting_SDSS_mass, SHIELD_mass)
	##fitting_SDSS_mass = np.append(fitting_SDSS_mass, Saintonge_mass)
	fitting_SDSS_mass = np.append(fitting_SDSS_mass, Jimmy_mass)
	#fitting_SDSS_metallicity = np.append(fitting_SDSS_metallicity, Lee_metallicity)
	fitting_SDSS_metallicity = np.append(fitting_SDSS_metallicity, Berg_metallicity)
	fitting_SDSS_metallicity = np.append(fitting_SDSS_metallicity, James_metallicity)
	fitting_SDSS_metallicity = np.append(fitting_SDSS_metallicity, SHIELD_metallicity)
	##fitting_SDSS_metallicity = np.append(fitting_SDSS_metallicity, Saintonge_metallicity)
	fitting_SDSS_metallicity = np.append(fitting_SDSS_metallicity, Jimmy_metallicity)
	
	medians_array = np.array([])
	medians_x_axis = np.array([])
	running_residuals_array = np.array([])
	for bin_limit in mass_bins_lower_limit:
		#print(bin_limit)
		#print(bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))
		#np.logical_and
		#find median of the y-axis values for every galaxy in this bin
		in_bin = np.logical_and(fitting_SDSS_mass>bin_limit,fitting_SDSS_mass<(bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])))
		#print(fitting_SDSS_metallicity[in_bin])
		if len(fitting_SDSS_metallicity[in_bin]) > binning_cut:
			median = np.mean(fitting_SDSS_metallicity[in_bin])
			#median = np.sum(fitting_SDSS_metallicity[in_bin])/len(fitting_SDSS_metallicity[in_bin])
			##if bin_limit == 7.5:
			##	plt.scatter(fitting_SDSS_mass[in_bin],fitting_SDSS_metallicity[in_bin], color='blue')
			##if bin_limit == 7.25:
			##	plt.scatter(fitting_SDSS_mass[in_bin],fitting_SDSS_metallicity[in_bin], color='red')
			##if bin_limit == 7.0:
			##	plt.scatter(fitting_SDSS_mass[in_bin],fitting_SDSS_metallicity[in_bin], color='purple')
			if np.isfinite(median):
				#print(median)
				medians_array = np.append(medians_array,median)
				medians_x_axis = np.append(medians_x_axis, bin_limit+((mass_bins_lower_limit[1]-mass_bins_lower_limit[0])/2))
				print('mass bin: '+str(bin_limit)+' mean value: '+str(round(median,2))+' scatter in bin: '+str(round(np.std(fitting_SDSS_metallicity[in_bin]-median),2))+' number of galaxies in bin: '+str(len(fitting_SDSS_metallicity[in_bin])))
				running_residuals_array = np.append(running_residuals_array, fitting_SDSS_metallicity[in_bin]-median)
				
	#medians_array = np.array([8.51819964804, 8.5881260704, 8.68400184466, 8.74884150215, 8.7835023075, 8.8081325197, 8.82203300249, 8.83091902361, 8.8185881194, 8.81050094938, 8.47969102404])
	#medians_x_axis = np.array([9.125, 9.375, 9.625, 9.875, 10.125, 10.375, 10.625, 10.875, 11.125, 11.375, 11.875])
	#print('medians_array: '+str(medians_array))
	#print('medians_x_axis: '+str(medians_x_axis))
	#ax.scatter(medians_x_axis, medians_array, marker='.', s=400, color='blue', label='median points')
	#fit_results, fit_error = optimization.curve_fit(fourth_order_fit, fitting_SDSS_mass-10, fitting_SDSS_metallicity)
	#print('4th order fit attempt: '+str(fit_results))
	fit_results, fit_error = optimization.curve_fit(fourth_order_fit, medians_x_axis-10, medians_array, p0=[8.96,0.31,-0.23,-0.017,0.046])
	print('4th order fit results: '+str(fit_results))
	print('Predicted metallicities based on this fit: '+str(fourth_order_fit(Jimmy_mass-10,*fit_results)))
	print('1 sigma scatter of the means: '+str(np.std(medians_array-fourth_order_fit(medians_x_axis-10, *fit_results))))
	#fit_results = np.polyfit(medians_x_axis-10, medians_array, 4)
	#print('4th order fit attempt: '+str(fit_results))
	#fit_results = [8.76,0.31,-0.23,-0.017,0.046]
	#ax.scatter(fitting_SDSS_mass, fitting_SDSS_metallicity, color='green', alpha=0.5)
	n = len(fitting_SDSS_mass)			   # number of samples    
	#print('fit results: '+str(fit_results)) 
	polynomial_fit = fit_results[0]+fit_results[1]*(fitting_SDSS_mass-10)+fit_results[2]*(fitting_SDSS_mass-10)**2+fit_results[3]*(fitting_SDSS_mass-10)**3+fit_results[4]*(fitting_SDSS_mass-10)**4
	residuals = fitting_SDSS_metallicity - polynomial_fit
	var_res = np.sum(residuals**2)/(n-2)
	sd_res = np.sqrt(var_res)
	#print('polynomial_fit: '+str(polynomial_fit))
	#print('fitting_SDSS_metallicity: '+str(fitting_SDSS_metallicity))
	#print('residuals: '+str(residuals))
	#print('max residuals: '+str(residuals[residuals<-1]))
	#print('max residuals mass: '+str(fitting_SDSS_mass[residuals<-1]))
	#print('max residuals metallicity: '+str(fitting_SDSS_metallicity[residuals<-1]))
	#print('max residuals fit metallicity: '+str(polynomial_fit[residuals<-1]))
	#print('var res: '+str(var_res))
	#print('n: '+str(n))
	print('Scatter AKA standard deviation: '+str(round(sd_res,10)))
	#print(running_residuals_array)
	print('Alternative Scatter AKA standard deviation: '+str(round(np.std(running_residuals_array),10)))
	#H, xedges, yedges = np.histogram2d(sdss_sample_PP04_O3N2, sdss_table['lgm_tot_p50'])
	#~before a statment is the equivalent of not, I should remember that.
	#H, xedges, yedges = np.histogram2d(cleaned_metallicity, cleaned_mass)
	ax.fill_between(range(20), 9.1, 12.0, facecolor='yellow', alpha=0.25)
	ax.fill_between(range(20), 7.0, 7.6, facecolor='yellow', alpha=0.25)
	#extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
	#levels = (1.0e3, 2.5e2, 1.0e2, 4.0e1)
	#cset = plt.contour(H, levels, origin='lower',colors='purple',linewidths=(3.0, 2.5, 2.0, 1.5),extent=extent, label='SDSS Contours')
	#cset.collections[0].set_label("SDSS z<0.03")
	h, x, y, p = ax.hist2d(SDSS_mass, SDSS_metallicity, bins = [20,20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[7.75, 11],[7.6, 9.1]])
	x = x[:-1]
	y = y[:-1]
	if (calibration == 'N2'):
		ax.scatter(SDSS_mass[SDSS_mass<7.75], SDSS_metallicity[SDSS_mass<7.75], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
	if (calibration == 'O3N2'):
		ax.scatter(SDSS_mass[SDSS_mass<7.75], SDSS_metallicity[SDSS_mass<7.75], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (O3N2)', rasterized=True)
	for x_index, x_limit in enumerate(x):
		for y_index, y_limit in enumerate(y):
			if np.isnan(h[x_index,y_index]):
				for sdss_index, temp_mass in enumerate(SDSS_mass):
					if (temp_mass > x_limit) and (temp_mass < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
						ax.scatter(temp_mass, SDSS_metallicity[sdss_index], color='gray', marker='o', alpha=sdss_alpha_value)
	
	#print(min(x))
	#print(max(x))
	#print(x)
	#stop()
	Mannucci_x = np.linspace(min(medians_x_axis), max(medians_x_axis), 100)
	#Jimmy_y = fit_results[0]+fit_results[1]*(Mannucci_x-10)+fit_results[2]*(Mannucci_x-10)**2+fit_results[3]*(Mannucci_x-10)**3+fit_results[4]*(Mannucci_x-10)**4
	Jimmy_y = fourth_order_fit(Mannucci_x-10, *fit_results) #O3N2 requires adding this little guy: -10
	ax.scatter(medians_x_axis, medians_array, color='orange', s=500, linewidth='5', marker="_", label="Bin Means")
	ax.plot(Mannucci_x, Jimmy_y, color='green', linewidth=4.0, label='Fit to Means')
	print(Mannucci_x)
	print(Jimmy_y) 
	if (calibration == 'N2'):
		#plt.scatter(SDSS_mass, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
		#ax.scatter([0,0], [0,0], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
		#ax.scatter(Lee_mass, Lee_metallicity, s=40, facecolors='none', edgecolors='green', label='Lee+06 (Literature)')
		ax.scatter(Berg_mass, Berg_metallicity, s=40, color='red', label='Berg+12 (N2)')
		ax.scatter(James_mass, James_metallicity, s=40, color='purple', label='James+14 (N2)')
		ax.scatter(SHIELD_mass, SHIELD_metallicity, s=40, color='blue', label='SHIELD (N2)')
		ax.scatter(Saintonge_mass,Saintonge_metallicity, s=40, facecolors='black', edgecolors='black', marker='D', label="Saintonge (N2)")
		ax.errorbar(Jimmy_mass,Jimmy_metallicity,yerr=Jimmy_metallicity_error, xerr=Jimmy_mass_error, linestyle="None", color='black')
		ax.scatter(Jimmy_mass,Jimmy_metallicity, s=400, color='black', marker='*', label="Jimmy+15 (N2)")
		##ax.scatter(Berg_mass, Berg_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
		##ax.scatter(James_mass, James_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
		##ax.scatter(SHIELD_mass, SHIELD_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
		##ax.scatter(Saintonge_mass,Saintonge_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
		##ax.scatter(Jimmy_mass,Jimmy_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
		##ax.scatter(medians_x_axis, medians_array, s=400, color='blue')
	if (calibration == 'O3N2'):
		#ax.scatter(SDSS_mass, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (O3N2)', rasterized=True)
		#ax.scatter([0,0], [0,0], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (O3N2)', rasterized=True)
		#ax.scatter(Lee_mass, Lee_metallicity, s=40, facecolors='none', edgecolors='green', label='Lee+06 (Literature)')
		ax.scatter(Berg_mass, Berg_metallicity, s=40, color='red', label='Berg+12 (O3N2)')
		ax.scatter(James_mass, James_metallicity, s=40, color='purple', label='James+14 (O3N2)')
		ax.scatter(SHIELD_mass, SHIELD_metallicity, s=40, color='blue', label='SHIELD (O3N2)')
		ax.scatter(Jimmy_mass,Jimmy_metallicity, s=400, color='black', marker='*', label="Jimmy+15 (O3N2)")
		ax.scatter(Saintonge_mass,Saintonge_metallicity, s=40, facecolors='black', edgecolors='black', marker='D', label="Saintonge (O3N2)")
		ax.errorbar(Jimmy_mass,Jimmy_metallicity,yerr=Jimmy_metallicity_error, xerr=Jimmy_mass_error, linestyle="None", color='black')
	#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
	plt.xlabel('log($M_{*}$)',fontsize=fontsize, labelpad=20)
	plt.ylabel('12+log(O/H)      ',fontsize=fontsize, labelpad=10)
	plt.ylim(7.5,9.5)
	plt.xlim(6.0,11.0)
	plt.minorticks_on()
	plt.text(solar_x_axis, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	#plt.legend(loc='upper left', ncol=2)
	artist = plt.Rectangle((0, 0), 1, 1, fc="gray", ec=None, alpha = sdss_alpha_value)
	handles, labels = ax.get_legend_handles_labels()
	handles = handles+[artist]
	if (calibration == 'N2'):
		labels = labels+[u'ALFALFA/SDSS (N2)']
	if (calibration == 'O3N2'):
		labels = labels+[u'ALFALFA/SDSS (O3N2)']
	ax.legend(handles, labels, loc='upper left', ncol=2)
	plt.subplots_adjust(bottom=bottom_margin)
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	#initial_guess = np.array([-1.492, 1.847, -0.08026])
	# Create a model for fitting.
	#poly_model = Model(second_order_fit)
	#quad_model = Model(quad_func)
	#moustakas_model = Model(odr_moustakas)
	#moustakas_guess = np.array([8.901, 8.798, 0.640])
	# Create a RealData object using our initiated data from above.
	#data = RealData(Jimmy_mass, Jimmy_metallicity, sx=Jimmy_mass_error, sy=Jimmy_metallicity_error)
	#SDSS_mass_error = SDSS_mass*0.0+0.5
	#SDSS_metallicity_error = SDSS_metallicity*0.0+0.5
	#data = RealData(SDSS_mass, SDSS_metallicity, sx=SDSS_mass_error, sy=SDSS_metallicity_error)
	# Set up ODR with the model and data.
	#odr = ODR(data, poly_model, beta0=[0.0])
	#odr = ODR(data, quad_model, beta0=initial_guess)
	#odr = ODR(data, moustakas_model, beta0=moustakas_guess)
	# Run the regression.
	#out = odr.run()
	# Use the in-built pprint method to give us results.
	#out.pprint()
	#x_fit = np.linspace(6, 11, 1000)
	#asm = 8.798
	#mto = 8.901
	#gamma = 0.64
	#y_fit = asm-np.log10(1+((10**(mto-x_fit)))**gamma)
	#y_fit = moustakas(x_fit, *moustakas_guess)
	#fit_results, fit_error = optimization.curve_fit(moustakas, SDSS_mass, SDSS_metallicity, p0 = moustakas_guess)
	#fit_error = np.sqrt(np.diag(fit_results[1]))
	#print('fit results: '+str(fit_results))
	#print('fit error: '+str(fit_error))
	#slope_error = fit_error[1]
	#fit_results = fit_results[0]
	#tremonti_relation = -1.492+1.847*(logmstar)-0.08026*(logmstar**2)
	#fit = np.polyfit(SDSS_mass, SDSS_metallicity, 2)
	#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
	#plt.plot(x_fit, y_fit, color='black', label='Fit')
	#y_fit = odr_moustakas(out.beta, x_fit)
	#plt.plot(x_fit, y_fit, color='red', label='ODR Fit')
	#y_fit = moustakas(x_fit, 9.298, 8.91, 0.64)
	#plt.plot(x_fit, y_fit, color='green', label='BS Fit')
	if (calibration == 'N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	if (calibration == 'O3N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
	#plt.show()
	
if (plot_to_make == 'colorful_mass_metallicity_sfr') or (plot_to_make == 'colorful_luminosity_metallicity_sfr'):
	super_residuals_array = np.array([])
	big_table_array = np.zeros((5,len(sfr_bins_lower_limit), len(mass_bins_lower_limit))) #min_mass, max_mass, metallicity, scatter, number_of_galaxies_in_bin
	#print(big_table_array)
	#stop()
	#Bootstrapping shit
	#high_bin_offset = np.array([])
	#low_bin_offset = np.array([])
	#high_bin_slope = np.array([])
	#low_bin_slope = np.array([])
	#for iteration in xrange(10):
	#	random_index = plt.randint(0,11, 11)
	#	mass = Jimmy_mass[random_index]
	#	metallicity = Jimmy_metallicity[random_index]
	#	sfr = Jimmy_sfr[random_index]
	#	
	#	high_sfr_mass = mass[sfr>sfr_bins_lower_limit[1]]
	#	low_sfr_mass = mass[sfr<sfr_bins_lower_limit[1]]
	#	high_sfr_metallicity = metallicity[sfr>sfr_bins_lower_limit[1]]
	#	low_sfr_metallicity = metallicity[sfr<sfr_bins_lower_limit[1]]
	#	
	#	#print(sfr_bins_lower_limit)
	#	#print(sfr)
	#	#print(high_sfr_mass)
	#	#print(low_sfr_mass)
	#	#print(len(high_sfr_mass))
	#	#print(len(low_sfr_mass))
	#	
	#	if (len(low_sfr_mass) > 2) and (len(high_sfr_mass) > 2):
	#		high_fit_results, high_fit_error = optimization.curve_fit(linear_fit, high_sfr_mass, high_sfr_metallicity)
	#		low_fit_results, low_fit_error = optimization.curve_fit(linear_fit, low_sfr_mass, low_sfr_metallicity)
	#	
	#		high_bin_offset = np.append(high_bin_offset, high_fit_results[0])
	#		low_bin_offset = np.append(low_bin_offset, low_fit_results[0])
	#		high_bin_slope = np.append(high_bin_slope, high_fit_results[1])
	#		low_bin_slope = np.append(low_bin_slope, low_fit_results[1])
	#print('median high offset: '+str(np.median(high_bin_offset))+'standard deviation: '+str(np.std(high_bin_offset)))
	#print('median low offset: '+str(np.median(low_bin_offset))+'standard deviation: '+str(np.std(low_bin_offset)))
	#print('median high slope: '+str(np.median(high_bin_slope))+'standard deviation: '+str(np.std(high_bin_slope)))
	#print('median low slope: '+str(np.median(low_bin_slope))+'standard deviation: '+str(np.std(low_bin_slope)))
	#high_sfr_bin_std_offset = np.std(high_bin_offset)
	#low_sfr_bin_std_offset = np.std(low_bin_offset)
	fig, ax = plt.subplots()
	# define the colormap
	cmap = plt.cm.jet
	#cmap = plt.cm.jet_r
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# force the first color entry to be grey
	#cmaplist[0] = (.5,.5,.5,1.0)
	#cmaplist[0] = (0.0,0.0,0.0,1.0)
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	#cmap.set_under('gray')
	# define the bins and normalize
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	#plt.fill_between(range(20), 9.1, 12.0, facecolor='yellow', alpha=0.25)
	#plt.fill_between(range(20), 7.0, 7.6, facecolor='yellow', alpha=0.25)
	if (plot_to_make == 'colorful_mass_metallicity_sfr'):
		h, x, y, p = ax.hist2d(SDSS_mass, SDSS_metallicity, bins = [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[7.75, 11],[7.6, 9.1]])
		x = x[:-1]
		y = y[:-1]
		ax.scatter(SDSS_mass[SDSS_mass<7.75], SDSS_metallicity[SDSS_mass<7.75], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		for x_index, x_limit in enumerate(x):
			for y_index, y_limit in enumerate(y):
				if np.isnan(h[x_index,y_index]):
					for sdss_index, temp_mass in enumerate(SDSS_mass):
						if (temp_mass > x_limit) and (temp_mass < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
							ax.scatter(temp_mass, SDSS_metallicity[sdss_index], color='gray', alpha=sdss_alpha_value)
		#ax.scatter([0], [0], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		#H, xedges,yedges = np.histogram2d(SDSS_metallicity,SDSS_mass,bins=(metallicity_bins_lower_limit,mass_bins_lower_limit))
		#cax = (plt.imshow(H, interpolation='nearest', origin='lower', aspect=1.0))
		#nbins=100
		#H, xedges, yedges = np.histogram2d(SDSS_mass,SDSS_metallicity,bins=nbins)
		#Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
		#plt.pcolormesh(xedges,yedges,Hmasked)
		#temp_mass_array = np.array([])
		#temp_metallicity_array = np.array([])
		#top_bin_array = np.array([])
		#bottom_bin_array = np.array([])
		#mass_plot_array = np.array([])
		#for lower_limit in mass_bins_lower_limit:
		#	for index, galaxy_mass in enumerate(SDSS_mass):
		#		if (galaxy_mass > lower_limit) and (galaxy_mass < lower_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])):
		#			temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
		#	print(np.median(temp_metallicity_array)) #find the median
		#	print(len(temp_metallicity_array)) #find the number of points in that bin
		#	if len(temp_metallicity_array) > 20:
		#		sorted_metallicity_array = np.sort(temp_metallicity_array)#sort them by value
		#		#let's say there's 20 items in the bin, and I want 64%, that means 12 points (rounded) stay
		#		#so I need to cut 8 points
		#		#I cut 4 from the top, and 4 from the bottom
		#		cut = round(len(temp_metallicity_array)*0.64/2)#cut the top x% and bottom x% out of the array
		#		print(cut)
		#		cut_metallicity_array = sorted_metallicity_array[cut:-cut]
		#		#top_bin_array = np.append(top_bin_array, cut_metallicity_array[0])#save the top and bottom value to an array
		#		#bottom_bin_array = np.append(bottom_bin_array, cut_metallicity_array[-1])
		#		top_bin_array = np.append(top_bin_array, np.percentile(sorted_metallicity_array, 90))
		#		bottom_bin_array = np.append(bottom_bin_array, np.percentile(sorted_metallicity_array, 10))
		#		mass_plot_array = np.append(mass_plot_array, lower_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])/2)
		#print(top_bin_array)
		#print(bottom_bin_array)
		#plt.fill_between(mass_plot_array, bottom_bin_array, top_bin_array, color='gray')
	if (plot_to_make == 'colorful_luminosity_metallicity_sfr'):
		h, x, y, p = ax.hist2d(SDSS_M_B, SDSS_metallicity, bins = [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[-20.96, -13.7],[7.6, 9.1]])
		x = x[:-1]
		y = y[:-1]
		ax.scatter(SDSS_M_B[SDSS_M_B>-13.7], SDSS_metallicity[SDSS_M_B>-13.7], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		for x_index, x_limit in enumerate(x):
			for y_index, y_limit in enumerate(y):
				if np.isnan(h[x_index,y_index]):
					for sdss_index, temp_M_B in enumerate(SDSS_M_B):
						if (temp_M_B > x_limit) and (temp_M_B < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
							ax.scatter(temp_M_B, SDSS_metallicity[sdss_index], color='gray', marker='o', alpha=sdss_alpha_value)
		#ax.scatter([0], [0], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	
	
	#create bins and find means within them.
	for sfr_bin_index, sfr_bin_limit in enumerate(sfr_bins_lower_limit):
		#print(hi_bin_limit)
		temp_sdss_metallicity_array = np.array([])
		temp_sdss_mass_array = np.array([])
		temp_sdss_M_B_array = np.array([])
		temp_sdss_sfr_array = np.array([])
		for index, galaxy_sfr in enumerate(SDSS_sfr):
			if (galaxy_sfr > sfr_bin_limit) and (galaxy_sfr < (sfr_bin_limit+(sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0]))):
				temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,SDSS_metallicity[index])
				temp_sdss_mass_array = np.append(temp_sdss_mass_array, SDSS_mass[index])
				temp_sdss_M_B_array = np.append(temp_sdss_M_B_array, SDSS_M_B[index])
				temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, galaxy_sfr)
		temp_SDSS_mass = temp_sdss_mass_array
		temp_SDSS_M_B = temp_sdss_M_B_array
		temp_SDSS_metallicity = temp_sdss_metallicity_array
		temp_SDSS_sfr = temp_sdss_sfr_array
		temp_residuals_array = np.array([])
		
		#moustakas_guess = np.array([8.901, 8.798, 0.640])
		#fit_results, fit_error = optimization.curve_fit(moustakas, temp_sdss_mass_array, temp_sdss_metallicity_array, p0 = moustakas_guess)
		#moustakas_model = Model(odr_moustakas)
		#data = RealData(temp_sdss_mass_array, temp_sdss_metallicity_array)
		#odr = ODR(data, moustakas_model, beta0=moustakas_guess)
		#out = odr.run()
		#out.pprint()
		print('for star formation rate bin starting at :'+str(sfr_bin_limit))
		if (plot_to_make == 'colorful_mass_metallicity_sfr'):
			temp_median_x_axis = np.array([])
			temp_median_y_axis = np.array([])
			for mass_bins_index, bin_limit in enumerate(mass_bins_lower_limit):
				temp_sdss_metallicity_array = np.array([])
				temp_sdss_mass_array = np.array([])
				temp_sdss_sfr_array = np.array([])
				for index, galaxy_mass in enumerate(temp_SDSS_mass):
					if (galaxy_mass > bin_limit) and (galaxy_mass < (bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))):
						temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
						temp_sdss_mass_array = np.append(temp_sdss_mass_array, temp_SDSS_mass[index])
						temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, temp_SDSS_sfr[index])
				if len(temp_sdss_metallicity_array) > binning_cut:
					temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))	
					print('galaxy mass: '+str(bin_limit)+' mean metallicity: '+str(round(np.mean(temp_sdss_metallicity_array),2))+' standard deviation: '+str(round(np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)),2)))
					#print(sfr_bin_index)
					#print(mass_bins_index)
					big_table_array[0,sfr_bin_index,mass_bins_index] = bin_limit
					big_table_array[1,sfr_bin_index,mass_bins_index] = bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])
					big_table_array[2,sfr_bin_index,mass_bins_index] = round(np.mean(temp_sdss_metallicity_array),2)
					big_table_array[3,sfr_bin_index,mass_bins_index] = round(np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)),2)
					big_table_array[4,sfr_bin_index,mass_bins_index] = len(temp_sdss_metallicity_array)
					#print(big_table_array[0,sfr_bin_index,mass_bins_index])
					temp_median_x_axis = np.append(temp_median_x_axis, bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))
					temp_residuals_array = temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)
					super_residuals_array = np.append(super_residuals_array,temp_residuals_array)
			color_value = (sfr_bin_limit-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0]) #+ 1/((sfr_bins_lower_limit[-1]-sfr_bins_lower_limit[0]))
			if color_value <= 1.0:
				ax.plot(temp_median_x_axis, temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
				std_dev = np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array))
				#print('for star formation rate bin starting at :'+str(sfr_bins_lower_limit[sfr_bin_index])+' stellar masses: '+str(temp_median_x_axis)+' metallicities: '+str(temp_median_y_axis)+' standard deviation: '+str(std_dev))
			#x_fit = np.linspace(6, 11, 1000)
			#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
			#y_fit = odr_moustakas(out.beta, x_fit)
			#plt.plot(x_fit, y_fit, color=cmap(color_value))
		#print('length of metallicity array: '+str(len(temp_SDSS_metallicity))+' scatter in this array: '+str(round(np.std(temp_residuals_array),2)))

		
		#stop()
		if (plot_to_make == 'colorful_luminosity_metallicity_sfr'):
			temp_median_x_axis = np.array([])
			temp_median_y_axis = np.array([])
			for mass_bins_index, bin_limit in enumerate(M_B_bins_lower_limit):
				temp_sdss_metallicity_array = np.array([])
				temp_sdss_M_B_array = np.array([])
				temp_sdss_sfr_array = np.array([])
				for index, galaxy_M_B in enumerate(temp_SDSS_M_B):
					if (galaxy_M_B > bin_limit) and (galaxy_M_B < bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0])):
						temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
						temp_sdss_M_B_array = np.append(temp_sdss_M_B_array, temp_SDSS_M_B[index])
						temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, temp_SDSS_sfr[index])
				if len(temp_sdss_metallicity_array) > binning_cut:
					temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
					print('luminosity: '+str(round(bin_limit,2))+' mean metallicity: '+str(round(np.mean(temp_sdss_metallicity_array),2))+' standard deviation: '+str(round(np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)),2)))
					big_table_array[0,sfr_bin_index,mass_bins_index] = bin_limit
					big_table_array[1,sfr_bin_index,mass_bins_index] = bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0])
					big_table_array[2,sfr_bin_index,mass_bins_index] = round(np.mean(temp_sdss_metallicity_array),2)
					big_table_array[3,sfr_bin_index,mass_bins_index] = round(np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)),2)
					big_table_array[4,sfr_bin_index,mass_bins_index] = len(temp_sdss_metallicity_array)
					temp_median_x_axis = np.append(temp_median_x_axis, bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0])/2)
					temp_residuals_array = temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)
					super_residuals_array = np.append(super_residuals_array,temp_residuals_array)
			color_value = (sfr_bin_limit-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
			if (color_value <= 1.0):
				ax.plot(temp_median_x_axis, temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
			#x_fit = np.linspace(6, 11, 1000)
			#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
			#y_fit = odr_moustakas(out.beta, x_fit)
			#plt.plot(x_fit, y_fit, color=cmap(color_value))
		print('length of metallicity array: '+str(len(temp_SDSS_metallicity))+' scatter in this array: '+str(round(np.std(temp_residuals_array),2)))
		print('running length of super residuals array: '+str(len(super_residuals_array)))
		
	
	if (plot_to_make == 'colorful_mass_metallicity_sfr'):
		temp_metallicity_array = np.array([])
		temp_mass_array = np.array([])
		temp_mass_array_error = np.array([])
		temp_metallicity_array_error = np.array([])
		temp_sfr_array = np.array([])
		#print(Jimmy_sfr)
		#print(Jimmy_metallicity)
		#stop()
		for index, galaxy_sfr in enumerate(Jimmy_sfr):
			if (galaxy_sfr < sfr_bins_lower_limit[2]) and (galaxy_sfr > sfr_bins_lower_limit[1]):
				temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
				temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
				temp_mass_array_error = np.append(temp_mass_array_error, Jimmy_mass_error[index])
				temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
				temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
		
		#for index, galaxy_sfr in enumerate(SHIELD_sfr):
		#	if (SHIELD_metallicity[index] < 11) and (SHIELD_sfr[index] <  sfr_bins_lower_limit[2]) and (SHIELD_sfr[index] >  sfr_bins_lower_limit[1]):
		#		temp_sfr_array = np.append(temp_sfr_array, SHIELD_sfr[index])
		#		temp_mass_array = np.append(temp_mass_array, SHIELD_mass[index])
		#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
		#for index, galaxy_sfr in enumerate(SDSS_sfr):
		#	if (galaxy_sfr < sfr_bins_lower_limit[2]) and (galaxy_sfr > sfr_bins_lower_limit[1]):
		#		temp_sfr_array = np.append(temp_sfr_array, SDSS_sfr[index])
		#		temp_mass_array = np.append(temp_mass_array, SDSS_mass[index])
		#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
		#moustakas_guess = np.array([8.901, 8.798, 0.640])
		#print(temp_metallicity_array)
		##fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array, temp_metallicity_array)
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
		##fit_results, fit_error = optimization.curve_fit(fourth_order_fit, temp_mass_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
		x_fit = np.linspace(min(temp_mass_array), max(temp_mass_array), 10) #len(temp_mass_array)
		#perr = np.sqrt(np.diag(fit_error))
		#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
		y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
		##y_fit = fourth_order_fit(x_fit, *fit_results)
		color_value = (sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
		ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
		ax.errorbar(temp_mass_array,temp_metallicity_array,yerr=temp_metallicity_array_error, xerr=temp_mass_array_error, linestyle="None", ecolor=cmap(color_value))
		#goodIndex = np.invert(np.logical_or(np.isnan(temp_mass_array), np.isnan(temp_metallicity_array)))
		##(a,b,(ci_a, ci_b), ri) = fitLine(temp_mass_array,temp_metallicity_array)
		#bootstrap: high_sfr_bin_std_offset
		#standard deviation: np.sqrt(fit_error[0][0])
		#other standard deviation: ri['sd_res']
		#confidence interval: ci_a[0]
		#print('ci_a'+str(ci_a)print(temp_metallicity_array)
		##residuals = temp_metallicity_array - fourth_order_fit(temp_mass_array, *fit_results) #(fit_results[0]+(fit_results[1]*temp_mass_array))
		residuals = temp_metallicity_array - linear_fit(temp_mass_array, *fit_results) #(fit_results[0]+(fit_results[1]*temp_mass_array))
		var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
		sd_res = np.sqrt(var_res)
		#print('fit_results for color: '+str(color_value)+':'+' slope: '+str(round(fit_results[1],2))+' intercept: '+str(round(fit_results[0],2))+' scatter: '+str(round(sd_res,2)))
		print(str(sfr_bins_lower_limit[1])+' to '+str(sfr_bins_lower_limit[2])+' & '+str(round(fit_results[1],2))+' $\pm$ '+str(round((fit_error[1])[1],2))+' & '+str(round(fit_results[0],2))+' $\pm$ '+str(round((fit_error[0])[0],2))+' & '+str(round(sd_res,2))+'\\\\')
		#print('Calculated standard deviation: '+str(sd_res))
		y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
		y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
		##y_fit_upper = fourth_order_fit(x_fit, fit_results[0]+sd_res, fit_results[1], fit_results[2], fit_results[3], fit_results[4])
		##y_fit_lower = fourth_order_fit(x_fit, fit_results[0]-sd_res, fit_results[1], fit_results[2], fit_results[3], fit_results[4])
		
		#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(color_value))
		ax.fill_between(x_fit, y_fit_upper, y_fit_lower,color=cmap(color_value), alpha=uncertainty_box_alpha)
		#print('ri: '+str(((ci_a[1]-fit_results[0])/2)))
		
		temp_metallicity_array = np.array([])
		temp_mass_array = np.array([])
		temp_mass_array_error = np.array([])
		temp_metallicity_array_error = np.array([])
		temp_sfr_array = np.array([])
		for index, galaxy_sfr in enumerate(Jimmy_sfr):
			if (galaxy_sfr < sfr_bins_lower_limit[1]) and (galaxy_sfr > sfr_bins_lower_limit[0]):
				temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
				temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
				temp_mass_array_error = np.append(temp_mass_array_error, Jimmy_mass_error[index])
				temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
				temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
		#for index, galaxy_sfr in enumerate(SHIELD_sfr):
		#	if (SHIELD_metallicity[index] < 11) and (SHIELD_sfr[index] <  sfr_bins_lower_limit[3]) and (SHIELD_sfr[index] >  sfr_bins_lower_limit[2]):
		#		temp_sfr_array = np.append(temp_sfr_array, SHIELD_sfr[index])
		#		temp_mass_array = np.append(temp_mass_array, SHIELD_mass[index])
		#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
		#for index, galaxy_sfr in enumerate(SDSS_sfr):
		#	if (galaxy_sfr < sfr_bins_lower_limit[1]) and (galaxy_sfr > sfr_bins_lower_limit[0]):
		#		temp_sfr_array = np.append(temp_sfr_array, SDSS_sfr[index])
		#		temp_mass_array = np.append(temp_mass_array, SDSS_mass[index])
		#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
		#moustakas_guess = np.array([8.901, 8.798, 0.640])
		#print(temp_metallicity_array)
		##fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array, temp_metallicity_array)
		try:
			fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
		except:
			fit_results = [0,0]
			fit_error = [[0,0],[0,0]]
		x_fit = np.linspace(min(temp_mass_array), max(temp_mass_array), 10)
		#perr = np.sqrt(np.diag(fit_error))
		#print(perr)
		#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*temp_mass_array))
		var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
		sd_res = np.sqrt(var_res)
		y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
		color_value = (sfr_bins_lower_limit[0]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
		#print('fit_results for color: '+str(color_value)+':'+' slope: '+str(round(fit_results[1],2))+' intercept: '+str(round(fit_results[0],2))+' scatter: '+str(round(sd_res,2)))
		print(str(sfr_bins_lower_limit[0])+' to '+str(sfr_bins_lower_limit[1])+' & '+str(round(fit_results[1],2))+' $\pm$ '+str(round((fit_error[1])[1],2))+' & '+str(round(fit_results[0],2))+' $\pm$ '+str(round((fit_error[0])[0],2))+' & '+str(round(sd_res,2))+'\\\\')
		ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
		ax.errorbar(temp_mass_array,temp_metallicity_array,yerr=temp_metallicity_array_error, xerr=temp_mass_array_error, linestyle="None", ecolor=cmap(color_value))
		##(a,b,(ci_a, ci_b), ri) = fitLine(temp_mass_array,temp_metallicity_array)
		y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
		y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
		#mpl.patches.Rectangle((8,8),1,1, facecolor='black', hatch='///', label='95% Confidence Interval')
		#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(color_value))
		ax.fill_between(x_fit, y_fit_upper, y_fit_lower,color=cmap(color_value), alpha=uncertainty_box_alpha)
		
	if (plot_to_make == 'colorful_luminosity_metallicity_sfr'):
		temp_metallicity_array = np.array([])
		temp_M_B_array = np.array([])
		temp_M_B_array_error = np.array([])
		temp_metallicity_array_error = np.array([])
		temp_sfr_array = np.array([])
		for index, galaxy_sfr in enumerate(Jimmy_sfr):
			if (galaxy_sfr < sfr_bins_lower_limit[1]) and (galaxy_sfr > sfr_bins_lower_limit[0]):
				temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
				temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
				temp_M_B_array_error = np.append(temp_M_B_array_error, Jimmy_M_B_error[index])
				temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
				temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
		#for index, galaxy_sfr in enumerate(SHIELD_sfr):
		#	if (SHIELD_metallicity[index] < 11) and (SHIELD_M_B[index] > 1.00) and (SHIELD_metallicity[index] > 8.0):
		#		temp_sfr_array = np.append(temp_sfr_array, SHIELD_sfr[index])
		#		temp_M_B_array = np.append(temp_M_B_array, SHIELD_M_B[index])
		#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
		#for index, galaxy_sfr in enumerate(SDSS_sfr):
		#	if (galaxy_sfr < sfr_bins_lower_limit[4]) and (galaxy_sfr > sfr_bins_lower_limit[3]):
		#		temp_sfr_array = np.append(temp_sfr_array, SDSS_sfr[index])
		#		temp_M_B_array = np.append(temp_M_B_array, SDSS_M_B[index])
		#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
		#moustakas_guess = np.array([8.901, 8.798, 0.640])
		#print(temp_metallicity_array)
		#fit_results, fit_error = optimization.curve_fit(moustakas, temp_M_B_array, temp_metallicity_array, p0 = moustakas_guess)
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
		x_fit = np.linspace(min(temp_M_B_array), max(temp_M_B_array), 10)
		#perr = np.sqrt(np.diag(fit_error))
		#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
		y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*temp_M_B_array))
		var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
		sd_res = np.sqrt(var_res)
		y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
		y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
		#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(0.17))
		color_value = (sfr_bins_lower_limit[0]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
		#print('fit_results for color: '+str(color_value)+':'+' slope: '+str(round(fit_results[1],2))+' intercept: '+str(round(fit_results[0],2))+' scatter: '+str(round(sd_res,2)))
		print(str(sfr_bins_lower_limit[0])+' to '+str(sfr_bins_lower_limit[1])+' & '+str(round(fit_results[1],2))+' $\pm$ '+str(round((fit_error[1])[1],2))+' & '+str(round(fit_results[0],2))+' $\pm$ '+str(round((fit_error[0])[0],2))+' & '+str(round(sd_res,2))+'\\\\')
		ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
		ax.fill_between(x_fit, y_fit_upper, y_fit_lower,color=cmap(color_value), alpha=uncertainty_box_alpha)
		ax.errorbar(temp_M_B_array,temp_metallicity_array,yerr=temp_metallicity_array_error, xerr=temp_M_B_array_error, linestyle="None", ecolor=cmap(color_value))
		##(a,b,(ci_a, ci_b), ri) = fitLine(temp_M_B_array,temp_metallicity_array)
	
		temp_metallicity_array = np.array([])
		temp_M_B_array = np.array([])
		temp_M_B_array_error = np.array([])
		temp_metallicity_array_error = np.array([])
		temp_sfr_array = np.array([])
		for index, galaxy_sfr in enumerate(Jimmy_sfr):
			if (galaxy_sfr < sfr_bins_lower_limit[2]) and (galaxy_sfr > sfr_bins_lower_limit[1]):
				temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
				temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
				temp_M_B_array_error = np.append(temp_M_B_array_error, Jimmy_M_B_error[index])
				temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
				temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
		#for index, galaxy_sfr in enumerate(SHIELD_sfr):
		#	if (SHIELD_metallicity[index] < 11) and (SHIELD_M_B[index] > 1.00) and (SHIELD_metallicity[index] > 8.0):
		#		temp_sfr_array = np.append(temp_sfr_array, SHIELD_sfr[index])
		#		temp_M_B_array = np.append(temp_M_B_array, SHIELD_M_B[index])
		#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
		#for index, galaxy_sfr in enumerate(SDSS_sfr):
		#	if (galaxy_sfr < sfr_bins_lower_limit[4]) and (galaxy_sfr > sfr_bins_lower_limit[3]):
		#		temp_sfr_array = np.append(temp_sfr_array, SDSS_sfr[index])
		#		temp_M_B_array = np.append(temp_M_B_array, SDSS_M_B[index])
		#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
		#moustakas_guess = np.array([8.901, 8.798, 0.640])
		#print(temp_metallicity_array)
		#fit_results, fit_error = optimization.curve_fit(moustakas, temp_M_B_array, temp_metallicity_array, p0 = moustakas_guess)
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
		x_fit = np.linspace(min(temp_M_B_array), max(temp_M_B_array), 10)
		#perr = np.sqrt(np.diag(fit_error))
		#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*temp_M_B_array))
		var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
		sd_res = np.sqrt(var_res)
		y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
		y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
		y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
		color_value = (sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
		#print('fit_results for color: '+str(color_value)+':'+' slope: '+str(round(fit_results[1],2))+' intercept: '+str(round(fit_results[0],2))+' scatter: '+str(round(sd_res,2)))
		print(str(sfr_bins_lower_limit[1])+' to '+str(sfr_bins_lower_limit[2])+' & '+str(round(fit_results[1],2))+' $\pm$ '+str(round((fit_error[1])[1],2))+' & '+str(round(fit_results[0],2))+' $\pm$ '+str(round((fit_error[0])[0],2))+' & '+str(round(sd_res,2))+'\\\\')		
		#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(color_value))
		ax.fill_between(x_fit, y_fit_upper, y_fit_lower,color=cmap(color_value), alpha=uncertainty_box_alpha)
		ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
		ax.errorbar(temp_M_B_array,temp_metallicity_array,yerr=temp_metallicity_array_error, xerr=temp_M_B_array_error, linestyle="None", ecolor=cmap(color_value))
		##(a,b,(ci_a, ci_b), ri) = fitLine(temp_M_B_array,temp_metallicity_array)
	
	print('length of super residuals array: '+str(len(super_residuals_array))+' scatter in this array: '+str(round(np.std(super_residuals_array),2)))
	#print(big_table_array[1])
	lows = big_table_array[0]
	highs = big_table_array[1]
	metals = big_table_array[2]
	scatters = big_table_array[3]
	n_gals = big_table_array[4]
	for limits in range(len(lows[sfr])):
		if max(lows[:,limits] != 0.0):
			#if (plot_to_make == 'colorful_mass_metallicity_sfr'):
			print(str(round(min(lows[:,limits]),2))+' to '+str(round(min(highs[:,limits]),2))+' & '+str(metals[0,limits])+' & '+str(metals[1,limits])+' & '+str(metals[2,limits])+' & '+str(metals[3,limits])+' & '+str(metals[4,limits])+'\\\\') #'Yes' if fruit == 'Apple' else 'No'
			#print(' & '+str(scatters[0,limits])+' & '+str(scatters[1,limits])+' & '+str(scatters[2,limits])+' & '+str(scatters[3,limits])+' & '+str(scatters[4,limits])+'\\\\')
			#print(' & '+str(n_gals[0,limits])+' & '+str(n_gals[1,limits])+' & '+str(n_gals[2,limits])+' & '+str(n_gals[3,limits])+' & '+str(n_gals[4,limits])+'\\\\')
			#for sfr in range(6):
			#	print(lows[sfr])
			#if (plot_to_make == 'colorful_luminosity_metallicity_sfr'):
			#	print(str(round(min(lows[:,limits]),2))+' to '+str(round(min(highs[:,limits]),2))+' & '+str(metals[0,limits])+' & '+str(metals[1,limits])+' & '+str(metals[2,limits])+' & '+str(metals[3,limits])+' & '+str(metals[4,limits])+'\\\\') #'Yes' if fruit == 'Apple' else 'No'
		
	#stop()
	if (plot_to_make == 'colorful_mass_metallicity_sfr'):
		ax.plot([0,0],[0,0], color='black', label=r'M$_*$ Bin Mean')
		#plt.scatter(SDSS_mass, SDSS_metallicity, c=SDSS_sfr-SDSS_mass, alpha=0.5, cmap=cmap, norm=norm, label="ALFALFA", rasterized=True)
		#Jimmy_metallicity[1] = 6.0
		#Jimmy_metallicity[5] = 6.0
		#Jimmy_metallicity[6] = 6.0
		#Jimmy_metallicity[8] = 6.0
		s = ax.scatter(Jimmy_mass, Jimmy_metallicity, c=Jimmy_sfr, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
		#print(names)
		#print(Jimmy_sfr)
		#plt.scatter(Jimmy_mass, Jimmy_metallicity-0.15, c=Jimmy_sfr, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
		#plt.scatter(James_mass, James_metallicity, c=James_sfr, s=40, cmap=cmap, norm=norm, marker='s', label="James+14")
		#plt.scatter(SHIELD_mass, SHIELD_metallicity, c=SHIELD_sfr, s=40, cmap=cmap, norm=norm, label="SHIELD", color='black')
		#plt.scatter(Jimmy_mass, Jimmy_metallicity, c=Jimmy_sfr-Jimmy_mass, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15")
		#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
	if (plot_to_make == 'colorful_luminosity_metallicity_sfr'):
		ax.plot([0,0],[0,0], color='black', label=r'Luminosity Bin Mean')
		#plt.scatter(SDSS_mass, SDSS_metallicity, c=SDSS_sfr-SDSS_mass, alpha=0.5, cmap=cmap, norm=norm, label="ALFALFA", rasterized=True)
		s = ax.scatter(Jimmy_M_B, Jimmy_metallicity, c=Jimmy_sfr, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
		#plt.scatter(James_mass, James_metallicity, c=James_sfr, s=40, cmap=cmap, norm=norm, marker='s', label="James+14")
		#plt.scatter(SHIELD_M_B, SHIELD_metallicity, c=SHIELD_sfr, s=40, cmap=cmap, norm=norm, label="SHIELD", color='black')
		#plt.scatter(Jimmy_mass, Jimmy_metallicity, c=Jimmy_sfr-Jimmy_mass, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15")
		#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
	
	ax.text(solar_x_axis, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
	ax.plot([0,0],[0,0], color='black', linestyle='--', label='Linear Fit')
	artist = plt.Rectangle((0, 0), 1, 1, fc="gray", ec=None, alpha = sdss_alpha_value)
	handles, labels = ax.get_legend_handles_labels()
	handles = handles+[artist]
	labels = labels+[u'ALFALFA/SDSS']
	ax.legend(handles, labels, loc='upper left', ncol=2)
	# sort both labels and handles by labels
	#labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
	#ax.legend(handles, labels)
	#logmstar = 7.0+(np.arange(400)*0.01)
	#PP04_line = 32.1488-(8.51258*logmstar)+(0.976384*logmstar*logmstar)-(0.0359763*logmstar*logmstar*logmstar)
	#plt.plot(logmstar,PP04_line, color="black", linewidth=3.0, label="PP04 SDSS Relation")
	plt.ylabel('12+log(O/H)',fontsize=fontsize, labelpad=10)
	plt.ylim(7.6, 9.3)
	if (plot_to_make == 'colorful_mass_metallicity_sfr'):
		plt.xlim(6.0, 11)
		plt.xlabel(r'log(M$_{*}$) [log(M$_{\odot}$)]',fontsize=fontsize, labelpad=20)
	if (plot_to_make == 'colorful_luminosity_metallicity_sfr'):
		plt.xlim(-10, -21)
		plt.xlabel(r'M$_B$',fontsize=fontsize, labelpad=20)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	plt.minorticks_on()
	#plt.legend(loc='upper left', ncol=2)
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin-0.05)
	plt.subplots_adjust(top=top_margin)
	#cbaxes = fig.add_axes([0.875, 0.175, 0.03, 0.775])
	#cbar = plt.colorbar(cax=cbaxes)
	#cbar.set_label(r'SFR (M$_\odot$ yr$^{-1}$', fontsize=fontsize, rotation=270)
	#ax2 = fig.add_axes([0.855, 0.175, 0.03, 0.775])
	#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, extend='min')
	cb = plt.colorbar(s)
	#plt.ticks.set_fontsize(20)
	cb.ax.tick_params(labelsize=fontsize-6)
	if (dog != 'y'):
		cb.set_label(r'SFR (M$_\odot$ yr$^{-1})$', fontsize=fontsize, rotation=270, labelpad=labelpad)
	if (dog == 'y'):
		cb.set_label(r'SFR (M$_\odot$ dyr$^{-1})$', fontsize=fontsize, rotation=270, labelpad=labelpad)
	if (calibration == 'N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	if (calibration == 'O3N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')

if (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity') or (plot_to_make == 'colorful_gas_fraction_sfr') or (plot_to_make == '3d_hi_mass_metallicity') or (plot_to_make == 'hi_mass_sfr') or (plot_to_make == 'hi_mass_stellar_mass') or (plot_to_make == 'hi_mass_metallicity') or (plot_to_make == 'hi_mass_density_sfr') or (plot_to_make == 'gas_fraction_metallicity') or (plot_to_make == 'table_of_everything') or (plot_to_make == '4d') or (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_sfr') or (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'lowest_scatter_hi_mass_lzr') or (plot_to_make == 'fmr_lzr_hi_mass') or (plot_to_make == 'consistency_check'):
	super_residuals_array = np.array([])
	big_table_array = np.zeros((5,len(sfr_bins_lower_limit), len(mass_bins_lower_limit))) #min_mass, max_mass, metallicity, scatter, number_of_galaxies_in_bin

	fig, ax = plt.subplots()
	#Jimmy_hi_mass = np.array([])
	#for alfalfa_index in range(len(alfalfa_sdss_catalog['AGCNr'])):
	#	for galaxy in names:
	#		if "AGC"+str((alfalfa_sdss_catalog['AGCNr'])[alfalfa_index]) == galaxy:
	#			Jimmy_hi_mass = np.append(Jimmy_hi_mass,(alfalfa_sdss_catalog['hi_mass'])[alfalfa_index])
	#for some reason the loop above is missing  7.22, so I'll do it manually
	Jimmy_hi_mass = np.array([7.74, 7.75, 7.60, 7.18, 7.41, 7.22, 7.46, 7.66, 7.67, 8.15, 7.43])
	#print(Jimmy_hi_mass)
	s21 = np.array([3.07, 0.62, 1.59, 0.24, 0.41, 0.26, 0.45, 0.69, 0.73, 1.08, 0.42])
	Jimmy_hi_mass = np.log10(2.356e5*Dp**2*s21)
	#print('after change')
	#print(Jimmy_hi_mass)
	#Jimmy_hi_mass_error = (Jimmy_hi_mass*0.0)+0.2
	s21_error = np.array([0.07, 0.04, 0.04, 0.04, 0.03, 0.05, 0.05, 0.05, 0.05, 0.04, 0.05])
	percent_error = (s21_error)/s21
	#print(Jimmy_hi_mass*percent_error)
	#Dp = np.array([8.7, 19.6, 10.3, 16.4, 16.4, 16.4, 16.5, 16.7, 16.6, 23.6, 16.6])
	Jimmy_hi_mass_error = Jimmy_hi_mass*percent_error #np.log10(2.356e5*Dp**2*s21_error) #assuming distance errors are negligable
	if (plot_to_make == 'colorful_hi_mass_metallicity') or (plot_to_make == 'colorful_hi_mass_luminosity_metallicity'):
		# define the colormap
		cmap = plt.cm.jet
		# extract all colors from the .jet map
		cmaplist = [cmap(i) for i in range(cmap.N)]
		# force the first color entry to be grey
		#cmaplist[0] = (.5,.5,.5,1.0)
		#cmaplist[0] = (0.0,0.0,0.0,1.0)
		# create the new map
		cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
		#cmap.set_under('gray')
		# define the bins and normalize
		norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
		#plt.fill_between(range(20), 9.1, 12.0, facecolor='yellow', alpha=0.25)
		#plt.fill_between(range(20), 7.0, 7.6, facecolor='yellow', alpha=0.25)
		#plt.scatter(SDSS_mass, SDSS_metallicity, c=SDSS_hi_mass, alpha=0.5, cmap=cmap, norm=norm, label="ALFALFA", rasterized=True)
		if (plot_to_make == 'colorful_hi_mass_metallicity'):
			h, x, y, p = ax.hist2d(SDSS_mass, SDSS_metallicity, bins = [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[7.75, 11],[7.6, 9.1]])
			x = x[:-1]
			y = y[:-1]
			ax.scatter(SDSS_mass[SDSS_mass<7.75], SDSS_metallicity[SDSS_mass<7.75], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
			for x_index, x_limit in enumerate(x):
				for y_index, y_limit in enumerate(y):
					if np.isnan(h[x_index,y_index]):
						for sdss_index, temp_mass in enumerate(SDSS_mass):
							if (temp_mass > x_limit) and (temp_mass < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
								ax.scatter(temp_mass, SDSS_metallicity[sdss_index], color='gray', marker='o', alpha=sdss_alpha_value)
			#plt.scatter(SDSS_mass, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
			#ax.scatter([0], [0], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		if (plot_to_make == 'colorful_hi_mass_luminosity_metallicity'):
			h, x, y, p = ax.hist2d(SDSS_M_B, SDSS_metallicity, bins = [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[-20.96, -13.7],[7.6, 9.1]])
			x = x[:-1]
			y = y[:-1]
			ax.scatter(SDSS_M_B[SDSS_M_B>-13.7], SDSS_metallicity[SDSS_M_B>-13.7], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
			for x_index, x_limit in enumerate(x):
				for y_index, y_limit in enumerate(y):
					if np.isnan(h[x_index,y_index]):
						for sdss_index, temp_M_B in enumerate(SDSS_M_B):
							if (temp_M_B > x_limit) and (temp_M_B < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
								ax.scatter(temp_M_B, SDSS_metallicity[sdss_index], color='gray', marker='o', alpha=sdss_alpha_value)
			#ax.scatter([0], [0], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		
		#create bins and find means within them.
		for hi_mass_bin_index, hi_bin_limit in enumerate(hi_mass_bins_lower_limit):
			#print(hi_bin_limit)
			temp_sdss_metallicity_array = np.array([])
			temp_sdss_mass_array = np.array([])
			temp_sdss_M_B_array = np.array([])
			temp_sdss_hi_mass_array = np.array([])
			for index, galaxy_hi_mass in enumerate(SDSS_hi_mass):
				if (galaxy_hi_mass > hi_bin_limit) and (galaxy_hi_mass < (hi_bin_limit+(hi_mass_bins_lower_limit[1]-hi_mass_bins_lower_limit[0]))):
					temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,SDSS_metallicity[index])
					temp_sdss_mass_array = np.append(temp_sdss_mass_array, SDSS_mass[index])
					temp_sdss_M_B_array = np.append(temp_sdss_M_B_array, SDSS_M_B[index])
					temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, galaxy_hi_mass)
			temp_SDSS_mass = temp_sdss_mass_array
			temp_SDSS_M_B = temp_sdss_M_B_array
			temp_SDSS_metallicity = temp_sdss_metallicity_array
			temp_SDSS_hi_mass = temp_sdss_hi_mass_array
			temp_residuals_array = np.array([])
			
			#moustakas_guess = np.array([8.901, 8.798, 0.640])
			#fit_results, fit_error = optimization.curve_fit(moustakas, temp_sdss_mass_array, temp_sdss_metallicity_array, p0 = moustakas_guess)
			
			print('for HI mass bin starting at :'+str(hi_bin_limit))
			if (plot_to_make == 'colorful_hi_mass_metallicity'):
				temp_median_x_axis = np.array([])
				temp_median_y_axis = np.array([])
				for mass_bins_index, bin_limit in enumerate(mass_bins_lower_limit):
					#print('bin_limit'+str(bin_limit))
					temp_sdss_metallicity_array = np.array([])
					temp_sdss_mass_array = np.array([])
					temp_sdss_hi_mass_array = np.array([])
					for index, galaxy_mass in enumerate(temp_SDSS_mass):
						if (galaxy_mass > bin_limit) and (galaxy_mass < (bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))):
							temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
							temp_sdss_mass_array = np.append(temp_sdss_mass_array, temp_SDSS_mass[index])
							temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, temp_SDSS_hi_mass[index])
					if len(temp_sdss_metallicity_array) > binning_cut:
						#print(temp_sdss_metallicity_array)
						temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
						print('galaxy mass: '+str(bin_limit)+' mean metallicity: '+str(round(np.mean(temp_sdss_metallicity_array),2))+' standard deviation: '+str(round(np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)),2)))	
						big_table_array[0,hi_mass_bin_index,mass_bins_index] = bin_limit
						big_table_array[1,hi_mass_bin_index,mass_bins_index] = bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])
						big_table_array[2,hi_mass_bin_index,mass_bins_index] = round(np.mean(temp_sdss_metallicity_array),2)
						big_table_array[3,hi_mass_bin_index,mass_bins_index] = round(np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)),2)
						big_table_array[4,hi_mass_bin_index,mass_bins_index] = len(temp_sdss_metallicity_array)
						temp_median_x_axis = np.append(temp_median_x_axis, bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])/2)
						temp_residuals_array = temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)
						super_residuals_array = np.append(super_residuals_array,temp_residuals_array)
					color_value = (hi_bin_limit-hi_mass_bins_lower_limit[0])/(hi_mass_bins_lower_limit[-2]-hi_mass_bins_lower_limit[0])
				ax.plot(temp_median_x_axis, temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
				#x_fit = np.linspace(6, 11, 1000)
				#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
				#plt.plot(x_fit, y_fit, color=cmap(color_value))
			#print('length of metallicity array: '+str(len(temp_SDSS_metallicity))+' scatter in this array: '+str(round(np.std(temp_residuals_array),2)))
			
			if (plot_to_make == 'colorful_hi_mass_luminosity_metallicity'):
				temp_median_x_axis = np.array([])
				temp_median_y_axis = np.array([])
				for mass_bins_index, bin_limit in enumerate(M_B_bins_lower_limit):
					#print('bin_limit'+str(bin_limit))
					temp_sdss_metallicity_array = np.array([])
					temp_sdss_M_B_array = np.array([])
					temp_sdss_hi_mass_array = np.array([])
					for index, galaxy_M_B in enumerate(temp_SDSS_M_B):
						if (galaxy_M_B > bin_limit) and (galaxy_M_B < (bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0]))):
							temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
							temp_sdss_M_B_array = np.append(temp_sdss_M_B_array, temp_SDSS_M_B[index])
							temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, temp_SDSS_hi_mass[index])
					if len(temp_sdss_metallicity_array) > binning_cut:
						#print(temp_sdss_metallicity_array)
						temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
						print('luminosity: '+str(round(bin_limit,2))+' mean metallicity: '+str(round(np.mean(temp_sdss_metallicity_array),2))+' standard deviation: '+str(round(np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)),2)))
						big_table_array[0,hi_mass_bin_index,mass_bins_index] = bin_limit
						big_table_array[1,hi_mass_bin_index,mass_bins_index] = bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])
						big_table_array[2,hi_mass_bin_index,mass_bins_index] = round(np.mean(temp_sdss_metallicity_array),2)
						big_table_array[3,hi_mass_bin_index,mass_bins_index] = round(np.std(temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)),2)
						big_table_array[4,hi_mass_bin_index,mass_bins_index] = len(temp_sdss_metallicity_array)
						temp_median_x_axis = np.append(temp_median_x_axis, bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0])/2)
						temp_residuals_array = temp_sdss_metallicity_array-np.mean(temp_sdss_metallicity_array)
						super_residuals_array = np.append(super_residuals_array,temp_residuals_array)
					color_value = (hi_bin_limit-hi_mass_bins_lower_limit[0])/(hi_mass_bins_lower_limit[-2]-hi_mass_bins_lower_limit[0])
				ax.plot(temp_median_x_axis, temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
				#x_fit = np.linspace(6, 11, 1000)
				#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
				#plt.plot(x_fit, y_fit, color=cmap(color_value))
			print('length of metallicity array: '+str(len(temp_SDSS_metallicity))+' scatter in this array: '+str(round(np.std(temp_residuals_array),2)))

			
		if (plot_to_make == 'colorful_hi_mass_metallicity'):
			temp_metallicity_array = np.array([])
			temp_mass_array = np.array([])
			temp_mass_array_error = np.array([])
			temp_metallicity_array_error = np.array([])
			temp_hi_mass_array = np.array([])
			for index, galaxy_hi_mass in enumerate(Jimmy_hi_mass):
				if (galaxy_hi_mass < hi_mass_bins_lower_limit[1]) and (galaxy_hi_mass > hi_mass_bins_lower_limit[0]):
					temp_hi_mass_array = np.append(temp_hi_mass_array, Jimmy_hi_mass[index])
					temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
					temp_mass_array_error = np.append(temp_mass_array_error, Jimmy_mass_error[index])
					temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
					temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
			#for index, galaxy_hi_mass in enumerate(Saintonge_hi_mass):
			#	if (Saintonge_hi_mass[index] < hi_mass_bins_lower_limit[1]) and (Saintonge_hi_mass[index] > hi_mass_bins_lower_limit[0]): #and (Saintonge_metallicity[index] > 8.0)
			#		temp_hi_mass_array = np.append(temp_hi_mass_array, Saintonge_hi_mass[index])	
			#		temp_mass_array = np.append(temp_mass_array, Saintonge_mass[index]	)
			#		temp_metallicity_array = np.append(temp_metallicity_array, Saintonge_metallicity[index])
			#for index, galaxy_hi_mass in enumerate(SHIELD_hi_mass):
			#	if (galaxy_hi_mass < hi_mass_bins_lower_limit[1]) and (galaxy_hi_mass > hi_mass_bins_lower_limit[0]) and (SHIELD_mass[index] > 0.0) and (SHIELD_metallicity[index] < 11):
			#		temp_hi_mass_array = np.append(temp_hi_mass_array, SHIELD_hi_mass[index])
			#		temp_mass_array = np.append(temp_mass_array, SHIELD_mass[index])
			#		print(SHIELD_metallicity[index])
			#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
			#print(temp_metallicity_array)
			#print(Jimmy_hi_mass)
			#moustakas_guess = np.array([8.901, 8.798, 0.640])
			#linear_guess = np.array([6.24683786, 0.26661621])
			#fit_results, fit_error = optimization.curve_fit(moustakas, temp_mass_array, temp_metallicity_array, p0 = moustakas_guess)
			#print(temp_mass_array)
			fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
			x_fit = np.linspace(min(temp_mass_array), max(temp_mass_array), 10)
			#perr = np.sqrt(np.diag(fit_error))
			#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
			y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
			residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*temp_mass_array))
			var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
			sd_res = np.sqrt(var_res)
			color_value = (hi_mass_bins_lower_limit[0]-hi_mass_bins_lower_limit[0])/(hi_mass_bins_lower_limit[-2]-hi_mass_bins_lower_limit[0])
			print('fit_results for color: '+str(color_value)+':'+' slope: '+str(round(fit_results[1],2))+' intercept: '+str(round(fit_results[0],2))+' scatter: '+str(round(sd_res,2)))
			#print(fit_error)
			#print(str(round((fit_error[0])[0],2)))
			print(str(hi_mass_bins_lower_limit[0])+' to '+str(hi_mass_bins_lower_limit[1])+' & '+str(round(fit_results[1],2))+' $\pm$ '+str(round((fit_error[1])[1],2))+' & '+str(round(fit_results[0],2))+' $\pm$ '+str(round((fit_error[0])[0],2))+' & '+str(round(sd_res,2))+'\\\\')
			##(a,b,(ci_a, ci_b), ri) = fitLine(temp_mass_array,temp_metallicity_array)
			ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
			ax.errorbar(temp_mass_array,temp_metallicity_array,yerr=temp_metallicity_array_error, xerr=temp_mass_array_error, linestyle="None", ecolor=cmap(color_value))
			y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
			y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
			#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(color_value))
			ax.fill_between(x_fit, y_fit_upper, y_fit_lower,color=cmap(color_value), alpha=uncertainty_box_alpha)
			#plt.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
			
			
			temp_metallicity_array = np.array([])
			temp_mass_array = np.array([])
			temp_mass_array_error = np.array([])
			temp_metallicity_array_error = np.array([])
			temp_hi_mass_array = np.array([])
			for index, galaxy_hi_mass in enumerate(Jimmy_hi_mass):
				if (galaxy_hi_mass <= hi_mass_bins_lower_limit[2]) and (galaxy_hi_mass >= hi_mass_bins_lower_limit[1]):
					temp_hi_mass_array = np.append(temp_hi_mass_array, Jimmy_hi_mass[index])
					temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
					temp_mass_array_error = np.append(temp_mass_array_error, Jimmy_mass_error[index])
					temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
					temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
			#for index, galaxy_hi_mass in enumerate(Saintonge_hi_mass):
			#	if (Saintonge_hi_mass[index] <= hi_mass_bins_lower_limit[2]) and (Saintonge_hi_mass[index] >= hi_mass_bins_lower_limit[1]) and (Saintonge_metallicity[index] > 0): #and (Saintonge_metallicity[index] > 8.0)
			#		temp_hi_mass_array = np.append(temp_hi_mass_array, Saintonge_hi_mass[index])	
			#		temp_mass_array = np.append(temp_mass_array, Saintonge_mass[index]	)
			#		temp_metallicity_array = np.append(temp_metallicity_array, Saintonge_metallicity[index])
			#for index, galaxy_hi_mass in enumerate(SHIELD_hi_mass):
			#	if (galaxy_hi_mass < -2.00) and (galaxy_hi_mass > -2.75) and (SHIELD_metallicity[index] < 11):
			#		temp_hi_mass_array = np.append(temp_hi_mass_array, SHIELD_hi_mass[index])
			#		temp_mass_array = np.append(temp_mass_array, SHIELD_mass[index])
			#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
			#moustakas_guess = np.array([8.901, 8.798, 0.640])
			#linear_guess = np.array([6.24683786, 0.26661621])
			#fit_results, fit_error = optimization.curve_fit(moustakas, temp_mass_array, temp_metallicity_array, p0 = moustakas_guess)
			#print(temp_mass_array)
			fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
			#perr = np.sqrt(np.diag(fit_error))
			#print(perr)
			#print(fit_results)
			#print(fit_error[1][1])
			x_fit = np.linspace(min(temp_mass_array), max(temp_mass_array), 10)
			#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
			y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
			residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*temp_mass_array))
			var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
			sd_res = np.sqrt(var_res)
			color_value = (hi_mass_bins_lower_limit[1]-hi_mass_bins_lower_limit[0])/(hi_mass_bins_lower_limit[-2]-hi_mass_bins_lower_limit[0])
			#print('fit_results for color: '+str(color_value)+':'+' slope: '+str(round(fit_results[1],2))+' intercept: '+str(round(fit_results[0],2))+' scatter: '+str(round(sd_res,2)))
			print(str(hi_mass_bins_lower_limit[1])+' to '+str(hi_mass_bins_lower_limit[2])+' & '+str(round(fit_results[1],2))+' $\pm$ '+str(round((fit_error[1])[1],2))+' & '+str(round(fit_results[0],2))+' $\pm$ '+str(round((fit_error[0])[0],2))+' & '+str(round(sd_res,2))+'\\\\')
			##(a,b,(ci_a, ci_b), ri) = fitLine(temp_mass_array,temp_metallicity_array)
			ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
			ax.errorbar(temp_mass_array,temp_metallicity_array,yerr=temp_metallicity_array_error, xerr=temp_mass_array_error, linestyle="None", ecolor=cmap(color_value))
			y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
			y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
			#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(color_value))
			ax.fill_between(x_fit, y_fit_upper, y_fit_lower,color=cmap(color_value), alpha=uncertainty_box_alpha)
			#plt.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)	
			
		if (plot_to_make == 'colorful_hi_mass_luminosity_metallicity'):
			temp_metallicity_array = np.array([])
			temp_M_B_array = np.array([])
			temp_M_B_array_error = np.array([])
			temp_metallicity_array_error = np.array([])
			temp_hi_mass_array = np.array([])
			for index, galaxy_hi_mass in enumerate(Jimmy_hi_mass):
				if (galaxy_hi_mass < hi_mass_bins_lower_limit[1]) and (galaxy_hi_mass > hi_mass_bins_lower_limit[0]):
					temp_hi_mass_array = np.append(temp_hi_mass_array, Jimmy_hi_mass[index])
					temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
					temp_M_B_array_error = np.append(temp_M_B_array_error, Jimmy_M_B_error[index])
					temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
					temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
			#for index, galaxy_hi_mass in enumerate(Saintonge_hi_mass):
			#	if (Saintonge_hi_mass[index] <= hi_mass_bins_lower_limit[1]) and (Saintonge_hi_mass[index] >= hi_mass_bins_lower_limit[0]) and (Saintonge_metallicity[index] > 1.0): #and (Saintonge_metallicity[index] > 8.0)
			#		temp_hi_mass_array = np.append(temp_hi_mass_array, Saintonge_hi_mass[index])	
			#		temp_M_B_array = np.append(temp_M_B_array, Saintonge_M_B[index]	)
			#		temp_metallicity_array = np.append(temp_metallicity_array, Saintonge_metallicity[index])
			#print(temp_metallicity_array)
			#for index, galaxy_hi_mass in enumerate(SHIELD_hi_mass):
			#	if (galaxy_hi_mass <= hi_mass_bins_lower_limit[1]) and (galaxy_hi_mass >= hi_mass_bins_lower_limit[0]):
			#		temp_hi_mass_array = np.append(temp_hi_mass_array, SHIELD_hi_mass[index])
			#		temp_M_B_array = np.append(temp_M_B_array, SHIELD_M_B[index])
			#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
			#moustakas_guess = np.array([8.901, 8.798, 0.640])
			#linear_guess = np.array([6.24683786, 0.26661621])
			#fit_results, fit_error = optimization.curve_fit(moustakas, temp_mass_array, temp_metallicity_array, p0 = moustakas_guess)
			fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
			#perr = np.sqrt(np.diag(fit_error))
			x_fit = np.linspace(min(temp_M_B_array), max(temp_M_B_array), 10)
			#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
			y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
			residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*temp_M_B_array))
			var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
			sd_res = np.sqrt(var_res)
			color_value = (sfr_bins_lower_limit[0]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
			#print('fit_results for color: '+str(color_value)+':'+' slope: '+str(round(fit_results[1],2))+' intercept: '+str(round(fit_results[0],2))+' scatter: '+str(round(sd_res,2)))
			
			print(str(hi_mass_bins_lower_limit[0])+' to '+str(hi_mass_bins_lower_limit[1])+' & '+str(round(fit_results[1],2))+' $\pm$ '+str(round((fit_error[1])[1],2))+' & '+str(round(fit_results[0],2))+' $\pm$ '+str(round((fit_error[0])[0],2))+' & '+str(round(sd_res,2))+'\\\\')
			##(a,b,(ci_a, ci_b), ri) = fitLine(temp_M_B_array,temp_metallicity_array)
			ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
			ax.errorbar(temp_M_B_array,temp_metallicity_array,yerr=temp_metallicity_array_error, xerr=temp_M_B_array_error, linestyle="None", ecolor=cmap(color_value))
			y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
			y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
			#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(0.0))
			ax.fill_between(x_fit, y_fit_upper, y_fit_lower,color=cmap(color_value), alpha=uncertainty_box_alpha)
			#plt.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
		
			
			temp_metallicity_array = np.array([])
			temp_M_B_array = np.array([])
			temp_M_B_array_error = np.array([])
			temp_metallicity_array_error = np.array([])
			temp_hi_mass_array = np.array([])
			for index, galaxy_hi_mass in enumerate(Jimmy_hi_mass):
				if (galaxy_hi_mass < hi_mass_bins_lower_limit[2]) and (galaxy_hi_mass > hi_mass_bins_lower_limit[1]):
					temp_hi_mass_array = np.append(temp_hi_mass_array, Jimmy_hi_mass[index])
					temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
					temp_M_B_array_error = np.append(temp_M_B_array_error, Jimmy_M_B_error[index])
					temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
					temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
			#for index, galaxy_hi_mass in enumerate(Saintonge_hi_mass):
			#	if (Saintonge_hi_mass[index] < hi_mass_bins_lower_limit[2]) and (Saintonge_hi_mass[index] > hi_mass_bins_lower_limit[1]) and (Saintonge_metallicity[index] > 1.0):
			#		temp_hi_mass_array = np.append(temp_hi_mass_array, Saintonge_hi_mass[index])	
			#		temp_M_B_array = np.append(temp_M_B_array, Saintonge_M_B[index]	)
			#		temp_metallicity_array = np.append(temp_metallicity_array, Saintonge_metallicity[index])
			#print(temp_metallicity_array)
			#for index, galaxy_hi_mass in enumerate(SDSS_hi_mass):
			#	if (galaxy_hi_mass < -2.00) and (galaxy_hi_mass > -2.75):
			#		temp_hi_mass_array = np.append(temp_hi_mass_array, SDSS_hi_mass[index])
			#		temp_M_B_array = np.append(temp_M_B_array, SDSS_M_B[index])
			#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
			#print(temp_metallicity_array)
			#moustakas_guess = np.array([8.901, 8.798, 0.640])
			#linear_guess = np.array([6.24683786, 0.26661621])
			#fit_results, fit_error = optimization.curve_fit(moustakas, temp_M_B_array, temp_metallicity_array, p0 = moustakas_guess)
			fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array, temp_metallicity_array, sigma=temp_metallicity_array_error)
			#perr = np.sqrt(np.diag(fit_error))
			x_fit = np.linspace(min(temp_M_B_array), max(temp_M_B_array), 10)
			#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
			y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
			residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*temp_M_B_array))
			var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
			sd_res = np.sqrt(var_res)
			color_value = (sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
			#print('fit_results for color: '+str(color_value)+':'+' slope: '+str(round(fit_results[1],2))+' intercept: '+str(round(fit_results[0],2))+' scatter: '+str(round(sd_res,2)))
			print(str(hi_mass_bins_lower_limit[1])+' to '+str(hi_mass_bins_lower_limit[2])+' & '+str(round(fit_results[1],2))+' $\pm$ '+str(round((fit_error[1])[1],2))+' & '+str(round(fit_results[0],2))+' $\pm$ '+str(round((fit_error[0])[0],2))+' & '+str(round(sd_res,2))+'\\\\')
			##(a,b,(ci_a, ci_b), ri) = fitLine(temp_M_B_array,temp_metallicity_array)
			ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
			ax.errorbar(temp_M_B_array,temp_metallicity_array,yerr=temp_metallicity_array_error, xerr=temp_M_B_array_error, linestyle="None", ecolor=cmap(color_value))
			y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
			y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
			#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, color="none",hatch="\\",edgecolor=cmap(0.17))
			ax.fill_between(x_fit, y_fit_upper, y_fit_lower,color=cmap(color_value), alpha=uncertainty_box_alpha)
			#plt.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
		
		print('length of super residuals array: '+str(len(super_residuals_array))+' scatter in this array: '+str(round(np.std(super_residuals_array),2)))
		print(big_table_array[1])
		lows = big_table_array[0]
		highs = big_table_array[1]
		metals = big_table_array[2]
		scatters = big_table_array[3]
		n_gals = big_table_array[4]
		for limits in range(len(lows[sfr])):
			if max(lows[:,limits] != 0.0):
				print(str(round(min(lows[:,limits]),2))+' to '+str(round(min(highs[:,limits]),2))+' & '+str(metals[0,limits])+' & '+str(metals[1,limits])+' & '+str(metals[2,limits])+' & '+str(metals[3,limits])+' & '+str(metals[4,limits])+'\\\\') #'Yes' if fruit == 'Apple' else 'No'
				#print(' & '+str(scatters[0,limits])+' & '+str(scatters[1,limits])+' & '+str(scatters[2,limits])+' & '+str(scatters[3,limits])+' & '+str(scatters[4,limits])+'\\\\')
				#print(' & '+str(n_gals[0,limits])+' & '+str(n_gals[1,limits])+' & '+str(n_gals[2,limits])+' & '+str(n_gals[3,limits])+' & '+str(n_gals[4,limits])+'\\\\')
				#for sfr in range(6):
				#	print(lows[sfr])
		#stop()
		ax.plot([0,0],[0,0], color='black', linestyle='--', label='Linear Fit')
		if (plot_to_make == 'colorful_hi_mass_metallicity'):
			ax.plot([0,0], [0,0], color='black', label=r'M$_*$ Bin Mean')
			s = ax.scatter(Jimmy_mass, Jimmy_metallicity, c=Jimmy_hi_mass, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
			#plt.scatter(Saintonge_mass, Saintonge_metallicity, facecolor='black', c=Saintonge_hi_mass, s=40, cmap=cmap, norm=norm, marker='D', label="Saintonge")
			#plt.scatter(SHIELD_mass, SHIELD_metallicity, c=SHIELD_hi_mass, s=40, cmap=cmap, norm=norm, label="SHIELD", color='black')
		if (plot_to_make == 'colorful_hi_mass_luminosity_metallicity'):
			ax.plot([0,0], [0,0], color='black', label=r'Luminosity Bin Mean')
			s = ax.scatter(Jimmy_M_B, Jimmy_metallicity, c=Jimmy_hi_mass, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
			#plt.scatter(Saintonge_M_B, Saintonge_metallicity, facecolor='black', c=Saintonge_hi_mass, s=40, cmap=cmap, norm=norm, marker='D', label="Saintonge")
			#plt.scatter(SHIELD_M_B, SHIELD_metallicity, c=SHIELD_hi_mass, s=40, cmap=cmap, norm=norm, label="SHIELD", color='black')
		#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
		plt.text(solar_x_axis, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
		logmstar = 7.0+(np.arange(400)*0.01)
		PP04_line = 32.1488-(8.51258*logmstar)+(0.976384*logmstar*logmstar)-(0.0359763*logmstar*logmstar*logmstar)
		#plt.plot(logmstar,PP04_line, color="black", linewidth=3.0, label="PP04 SDSS Relation")
		plt.ylabel('12+log(O/H)',fontsize=fontsize, labelpad=10)
		if (plot_to_make == 'colorful_hi_mass_metallicity'):
			plt.xlim(6.0, 11)		
			plt.xlabel(r'log(M$_{*}$) [log(M$_{\odot}$)]',fontsize=fontsize, labelpad=20)
		if (plot_to_make == 'colorful_hi_mass_luminosity_metallicity'):
			plt.xlim(-10.0, -21)
			plt.xlabel(r'M$_B$',fontsize=fontsize, labelpad=20)
		plt.ylim(7.6, 9.3)	
		plt.minorticks_on()
		plt.tick_params(axis='x', labelsize=20)
		plt.tick_params(axis='y', labelsize=20)
		artist = plt.Rectangle((0, 0), 1, 1, fc="gray", ec=None, alpha = sdss_alpha_value)
		handles, labels = ax.get_legend_handles_labels()
		handles = handles+[artist]
		labels = labels+[u'ALFALFA/SDSS']
		ax.legend(handles, labels, loc='upper left', ncol=2)
		#plt.legend(loc='upper left', ncol=2)
		plt.subplots_adjust(bottom=bottom_margin)	
		plt.subplots_adjust(left=left_margin)
		plt.subplots_adjust(right=right_margin-0.05)
		plt.subplots_adjust(top=top_margin)
		#cbaxes = fig.add_axes([0.875, 0.175, 0.03, 0.775])
		#cbar = plt.colorbar(cax=cbaxes)
		#cbar.set_label(r'log(HI Mass) [log(M$_{\odot}$)]', fontsize=fontsize, rotation=270)
		#ax2 = fig.add_axes([0.855, 0.175, 0.03, 0.775])
		#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
		#cb.set_label(r'log(HI Mass) [log(M$_{\odot}$)]', fontsize=fontsize, rotation=270)
		#cb = plt.colorbar()
		cb = plt.colorbar(s)
		cb.ax.tick_params(labelsize=fontsize-6)
		cb.set_label(r'log(HI Mass) [log(M$_{\odot}$)]', fontsize=fontsize, rotation=270, labelpad=labelpad)
		if (calibration == 'N2'):
			plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
		if (calibration == 'O3N2'):
			plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
		#plt.show()
	
	if (plot_to_make == '3d_hi_mass_metallicity'):
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(alfalfa_sdss_catalog['hi_mass'], alfalfa_sdss_catalog['lgm_tot_p50'], alfalfa_sdss_catalog['PP04_Metallicity'], zdir='z', color='gray', s=5, alpha=sdss_alpha_value)
		ax.scatter(Jimmy_hi_mass, Jimmy_mass, Jimmy_metallicity, zdir='z', color='black', s=40)
		ax.set_xlim(5.5, 12)
		ax.set_ylim(6.0, 11)
		ax.set_zlim(7.5, 9.5)
		ax.set_xlabel('HI Mass')
		ax.set_ylabel('Stellar Mass')
		ax.set_zlabel('Metallicity')
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
		#plt.show()
	if (plot_to_make == '4d'):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		p = ax.scatter(alfalfa_sdss_catalog['hi_mass'], alfalfa_sdss_catalog['lgm_tot_p50'], alfalfa_sdss_catalog['PP04_Metallicity'], zdir='z', c=alfalfa_sdss_catalog['sfr_tot_p50']-alfalfa_sdss_catalog['lgm_tot_p50'], s=5, alpha=sdss_alpha_value, cmap='jet_r', vmin=-11, vmax=-8)
		q = ax.scatter(Jimmy_hi_mass, Jimmy_mass, Jimmy_metallicity, zdir='z', c=jimmy_sfr-Jimmy_mass, s=40, cmap='jet_r', vmin=-11, vmax=-8)
		ax.set_xlim(5.5, 12)
		ax.set_ylim(6.0, 11)
		ax.set_zlim(7.5, 9.5)
		ax.set_xlabel('HI Mass')
		ax.set_ylabel('Stellar Mass')
		ax.set_zlabel('Metallicity')
		cbar = fig.colorbar(q)
		cbar.set_label(r'SSFR (M$_\odot$ yr$^{-1})/M_\odot$', fontsize=fontsize, rotation=270)
		#plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
		plt.show()
	if (plot_to_make == 'hi_mass_metallicity'):
		#plt.fill_between(range(20), 9.1, 12.0, facecolor='yellow', alpha=0.25)
		#plt.fill_between(range(20), 7.0, 7.6, facecolor='yellow', alpha=0.25)
		#plt.scatter(SDSS_hi_mass, SDSS_metallicity, alpha=0.2, color='gray', label="ALFALFA/SDSS", rasterized=True)
		cmap = plt.cm.jet
		# extract all colors from the .jet map
		cmaplist = [cmap(i) for i in range(cmap.N)]
		# force the first color entry to be grey
		#cmaplist[0] = (0.0,0.0,0.0,1.0)
		# create the new map
		cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
		#cmap.set_under('gray')
		# define the bins and normalize
		norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
		plt.scatter(SDSS_hi_mass, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label="ALFALFA/SDSS", rasterized=True)
		for mass_bin_limit in mass_bins_lower_limit:
			#print(hi_bin_limit)
			temp_sdss_metallicity_array = np.array([])
			temp_sdss_mass_array = np.array([])
			temp_sdss_hi_mass_array = np.array([])
			for index, galaxy_mass in enumerate(SDSS_mass):
				if (galaxy_mass > mass_bin_limit) and (galaxy_mass < (mass_bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))):
					temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,SDSS_metallicity[index])
					temp_sdss_mass_array = np.append(temp_sdss_mass_array, SDSS_mass[index])
					temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, SDSS_hi_mass[index])
			temp_SDSS_mass = temp_sdss_mass_array
			temp_SDSS_metallicity = temp_sdss_metallicity_array
			temp_SDSS_hi_mass = temp_sdss_hi_mass_array
			
			#moustakas_guess = np.array([8.901, 8.798, 0.640])
			#fit_results, fit_error = optimization.curve_fit(moustakas, temp_sdss_mass_array, temp_sdss_metallicity_array, p0 = moustakas_guess)
			
			temp_median_x_axis = np.array([])
			temp_median_y_axis = np.array([])
			for hi_bin_limit in hi_mass_bins_lower_limit:
				#print('bin_limit'+str(bin_limit))
				temp_sdss_metallicity_array = np.array([])
				temp_sdss_mass_array = np.array([])
				temp_sdss_hi_mass_array = np.array([])
				for index, galaxy_hi_mass in enumerate(temp_SDSS_hi_mass):
					if (galaxy_hi_mass > hi_bin_limit) and (galaxy_hi_mass < (hi_bin_limit+(hi_mass_bins_lower_limit[1]-hi_mass_bins_lower_limit[0]))):
						temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
						temp_sdss_mass_array = np.append(temp_sdss_mass_array, temp_SDSS_mass[index])
						temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, temp_SDSS_hi_mass[index])
				if len(temp_sdss_metallicity_array) > binning_cut:
					#print(temp_sdss_metallicity_array)
					temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
					temp_median_x_axis = np.append(temp_median_x_axis, hi_bin_limit+(hi_mass_bins_lower_limit[1]-hi_mass_bins_lower_limit[0])/2)
			#print('temp x-axis line: '+str(temp_median_x_axis))
			#print('temp y-axis line: '+str(temp_median_y_axis))
			color_value = (mass_bin_limit-mass_bins_lower_limit[0])/(mass_bins_lower_limit[-2]-mass_bins_lower_limit[0])
			#print(color_value)
			plt.plot(temp_median_x_axis, temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
			#x_fit = np.linspace(6, 11, 1000)
			#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
			#plt.plot(x_fit, y_fit, color=cmap(color_value))
		#plt.scatter(SHIELD_hi_mass, SHIELD_D02_N2, s=40, color='cyan', label='SHIELD')
		plt.scatter(SHIELD_hi_mass, SHIELD_D02_N2, s=40, c=SHIELD_mass, norm=norm, label='SHIELD')
		#plt.scatter(Saintonge_hi_mass, Saintonge_metallicity, s=40, facecolor='black', edgecolor='black', marker='D', label="Saintonge")
		plt.scatter(Saintonge_hi_mass, Saintonge_metallicity, s=40, facecolor='black', c=Saintonge_mass, norm=norm, marker='D', label="Saintonge")
		#plt.scatter(Jimmy_hi_mass, Jimmy_metallicity, s=400, color='black', marker='*', label="Jimmy+15")
		plt.scatter(Jimmy_hi_mass, Jimmy_metallicity, s=400, c=Jimmy_mass, norm=norm, marker='*', label="Jimmy+15", color='black')
		#plt.errorbar(Jimmy_hi_mass,Jimmy_metallicity,yerr=Jimmy_metallicity_error, xerr=Jimmy_hi_mass, linestyle="None", color='black')
		#plt.scatter(alfalfa_sdss_catalog['hi_mass'], alfalfa_sdss_catalog['PP04_Metallicity'], alpha=0.1, cmap='jet', c=alfalfa_sdss_catalog['lgm_tot_p50'], label="ALFALFA/SDSS", vmin=6, vmax=12, rasterized=True)
		#plt.scatter(Jimmy_hi_mass, Jimmy_metallicity, s=400, cmap='jet', c=Jimmy_mass, vmin=6, vmax=12, label="Jimmy+15")
		#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
		plt.text(6.3, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
		plt.ylabel('12+log(O/H)',fontsize=fontsize, labelpad=10)
		plt.xlabel(r'HI Mass [log(M$_\odot$)]',fontsize=fontsize, labelpad=20)
		plt.xlim(6.5, 11)
		plt.ylim(7.6, 9.3)
		plt.tick_params(axis='x', labelsize=20)
		plt.tick_params(axis='y', labelsize=20)	
		plt.minorticks_on()
		#cbaxes = fig.add_axes([0.875, 0.175, 0.03, 0.775])
		#cbar = plt.colorbar(cax=cbaxes)
		#cbar.set_label(r'log(Mass) [log(M$_{\odot}$)]', fontsize=fontsize, rotation=270)	
		plt.plot([0,0],[0,0], color='black', linewidth=3.0, label="SDSS Mean")
		plt.legend(loc='upper left', ncol=2)
		plt.subplots_adjust(bottom=bottom_margin)	
		plt.subplots_adjust(left=left_margin)
		plt.subplots_adjust(right=right_margin)
		plt.subplots_adjust(top=top_margin)
		#ax2 = fig.add_axes([0.855, 0.175, 0.03, 0.775])
		#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
		cb = plt.colorbar()
		cb.ax.tick_params(labelsize=fontsize-6)
		cb.set_label(r'log(M$_*$) [log(M$_{\odot}$)]', fontsize=fontsize, rotation=270, labelpad=labelpad)
		if (calibration == 'N2'):
			plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
		if (calibration == 'O3N2'):
			plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
		#plt.show()
	if (plot_to_make == 'gas_fraction_metallicity'):
		plt.fill_between(np.array(range(20))-10, 9.1, 12.0, facecolor='yellow', alpha=0.25)
		plt.fill_between(np.array(range(20))-10, 7.0, 7.6, facecolor='yellow', alpha=0.25)
		plt.scatter(SDSS_hi_mass-SDSS_mass, SDSS_metallicity, alpha=sdss_alpha_value, color='gray', label="ALFALFA/SDSS", rasterized=True)
		#plt.scatter(Lee_mass-Lee_hi_mass, Lee_metallicity, s=40, color='green', label='Lee+06 (Literature)')
		plt.scatter(SHIELD_hi_mass-SHIELD_mass, SHIELD_D02_N2, s=40, color='blue', label='SHIELD')
		#print(Jimmy_hi_mass)
		#print(Jimmy_mass)
		plt.scatter(Saintonge_hi_mass-Saintonge_mass, Saintonge_metallicity, s=40, facecolor='black', edgecolor='black', marker='D', label="Saintonge")
		plt.scatter(Jimmy_hi_mass-Jimmy_mass, Jimmy_metallicity, s=400, color='black', marker='*', label="Jimmy+15")
		gas_fraction_error = np.sqrt((((1/Jimmy_mass)**2)*Jimmy_hi_mass_error**2)+(((-Jimmy_hi_mass/Jimmy_mass**2)**2)*Jimmy_mass_error**2))
		#print('gas_fraction_error: '+str(gas_fraction_error))
		plt.errorbar(Jimmy_hi_mass-Jimmy_mass,Jimmy_metallicity,yerr=Jimmy_metallicity_error, xerr=gas_fraction_error, linestyle="None", color='black')
		#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
		plt.text(-3.25, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
		plt.ylabel('12+log(O/H)',fontsize=fontsize, labelpad=10)
		plt.xlabel(r'log(M$_{HI}$/M$_*$)',fontsize=fontsize, labelpad=20)
		plt.xlim(-3, 3)
		plt.ylim(7.6, 9.3)
		plt.minorticks_on()
		plt.tick_params(axis='x', labelsize=20)
		plt.tick_params(axis='y', labelsize=20)
		plt.legend(loc='upper left')
		plt.subplots_adjust(bottom=bottom_margin)	
		plt.subplots_adjust(left=left_margin)
		plt.subplots_adjust(right=right_margin)
		plt.subplots_adjust(top=top_margin)
		if (calibration == 'N2'):
			plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
		if (calibration == 'O3N2'):
			plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
		#plt.show()

if (plot_to_make == 'hi_mass_sfr'): 
	plt.scatter(SDSS_hi_mass, SDSS_sfr, color='gray', alpha=sdss_alpha_value, label=r"ALFALFA/SDSS (H$\alpha$)", rasterized=True)
	plt.scatter(huang_catalog['logMHI'], huang_catalog['logSFR'], color='red', alpha = 1.0, label='Huang+12 (SED)')
	plt.scatter(Jimmy_hi_mass, Jimmy_sfr, color='black', s=400, marker="*", label=r"Jimmy+15 (H$\alpha$)") #pull in the dwarf data
	plt.scatter(SHIELD_hi_mass, SHIELD_sfr, color='blue', s=40, label="SHIELD (CMD)")
	#x_axis = np.array(range(20))-6
	#print(x_axis)
	#y_axis = x_axis+12
	#print(y_axis)
	#plt.plot(y_axis, x_axis, color='black', label="1:1 Relation")
	plt.ylabel('SFR [log(M$_{\odot}$ yr$^{-1}$)]',fontsize=fontsize, labelpad=10)
	plt.xlabel(r'HI Mass [log(M$_{\odot}$)]',fontsize=fontsize, labelpad=20)
	plt.xlim(6.0, 11)
	plt.ylim(-6, 4)	
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)	
	plt.minorticks_on()
	plt.legend(loc='upper left')
	plt.subplots_adjust(bottom=bottom_margin)	
	#plt.show()
	#plt.savefig(HOME+'/Astro/Tex/2nd Paper/hi_mass_sfr.pdf')
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	
if (plot_to_make == '3d_hi_mass_sfr'):
	fig=plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(alfalfa_sdss_catalog['hi_mass'], alfalfa_sdss_catalog['lgm_tot_p50'], alfalfa_sdss_catalog['sfr_tot_p50'], zdir='z', color='gray', s=5, alpha=sdss_alpha_value)
	ax.scatter(Jimmy_hi_mass, Jimmy_mass, Jimmy_sfr, zdir='z', color='black', s=40)
	ax.set_xlim(5.5, 12)
	ax.set_ylim(5.5, 12)
	ax.set_zlim(-6, 3)
	ax.set_xlabel('HI Mass')
	ax.set_ylabel('Stellar Mass')
	ax.set_zlabel('Star Formation Rate')
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#plt.show()

if (plot_to_make == '3d_hi_metallicity_sfr'):
	fig=plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(alfalfa_sdss_catalog['hi_mass'], alfalfa_sdss_catalog['PP04_Metallicity'], alfalfa_sdss_catalog['sfr_tot_p50'], zdir='z', color='gray', s=5, alpha=sdss_alpha_value)
	ax.scatter(Jimmy_hi_mass, Jimmy_metallicity, Jimmy_sfr, zdir='z', color='black', s=40)
	ax.set_xlim(5.5, 12)
	ax.set_ylim(7.5, 9.5)
	ax.set_zlim(-6, 3)
	ax.set_xlabel('HI Mass')
	ax.set_ylabel('Metallicity')
	ax.set_zlabel('SFR')
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#plt.show()
	
if (plot_to_make == 'metallicity_sfr'):
	fig=plt.figure()	
	#plt.fill_between(np.array(range(20))-6, 9.1, 12.0, facecolor='yellow', alpha=0.25)
	#plt.fill_between(np.array(range(20))-6, 7.0, 7.6, facecolor='yellow', alpha=0.25)
	plt.scatter(SDSS_sfr, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label="ALFALFA/SDSS", rasterized=True)
	cmap = plt.cm.jet
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# force the first color entry to be grey
	#cmaplist[0] = (0.0,0.0,0.0,1.0)
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	#cmap.set_under('gray')
	# define the bins and normalize
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	for mass_bin_limit in mass_bins_lower_limit:
		#print(hi_bin_limit)
		temp_sdss_metallicity_array = np.array([])
		temp_sdss_mass_array = np.array([])
		temp_sdss_sfr_array = np.array([])
		for index, galaxy_mass in enumerate(SDSS_mass):
			if (galaxy_mass > mass_bin_limit) and (galaxy_mass < (mass_bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))):
				temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,SDSS_metallicity[index])
				temp_sdss_mass_array = np.append(temp_sdss_mass_array, SDSS_mass[index])
				temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, SDSS_sfr[index])
		temp_SDSS_mass = temp_sdss_mass_array
		temp_SDSS_metallicity = temp_sdss_metallicity_array
		temp_SDSS_sfr = temp_sdss_sfr_array
		
		#moustakas_guess = np.array([8.901, 8.798, 0.640])
		#fit_results, fit_error = optimization.curve_fit(moustakas, temp_sdss_mass_array, temp_sdss_metallicity_array, p0 = moustakas_guess)
		
		temp_median_x_axis = np.array([])
		temp_median_y_axis = np.array([])
		for sfr_bin_limit in sfr_bins_lower_limit:
			#print('bin_limit'+str(bin_limit))
			temp_sdss_metallicity_array = np.array([])
			temp_sdss_mass_array = np.array([])
			temp_sdss_sfr_array = np.array([])
			for index, galaxy_sfr in enumerate(temp_SDSS_sfr):
				if (galaxy_sfr > sfr_bin_limit) and (galaxy_sfr < (sfr_bin_limit+(sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0]))):
					temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
					temp_sdss_mass_array = np.append(temp_sdss_mass_array, temp_SDSS_mass[index])
					temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, temp_SDSS_sfr[index])
			if len(temp_sdss_metallicity_array) > binning_cut:
				#print(temp_sdss_metallicity_array)
				temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
				temp_median_x_axis = np.append(temp_median_x_axis, sfr_bin_limit+(sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])/2)
		#print('temp x-axis line: '+str(temp_median_x_axis))
		#print('temp y-axis line: '+str(temp_median_y_axis))
		color_value = (mass_bin_limit-mass_bins_lower_limit[0])/(mass_bins_lower_limit[-2]-mass_bins_lower_limit[0])
		#print(color_value)
		plt.plot(temp_median_x_axis, temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
		#x_fit = np.linspace(6, 11, 1000)
		#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
		#plt.plot(x_fit, y_fit, color=cmap(color_value))
	#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
	plt.text(-4.2, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
	plt.scatter(Jimmy_sfr, Jimmy_metallicity, c=Jimmy_mass, s=400, marker="*", norm=norm, label="Jimmy+15", color='black') #pull in the dwarf data
	plt.scatter(SHIELD_sfr, SHIELD_metallicity, c=SHIELD_mass, s=40, norm=norm, label="SHIELD", color='black')
	plt.xlabel('SFR [log(M$_{\odot}$ yr$^{-1}$)]',fontsize=fontsize, labelpad=10)
	plt.ylabel(r'12+log(O/H)',fontsize=fontsize, labelpad=20)
	plt.ylim(7.6, 9.3)
	plt.xlim(-4.0, 1.0)	
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)	
	plt.minorticks_on()
	plt.plot([0,0],[0,0], color='black', linewidth=3.0, label="SDSS Mean")
	plt.legend(loc='upper left', ncol=2)
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	#ax2 = fig.add_axes([0.855, 0.175, 0.03, 0.775])
	#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
	cb = plt.colorbar()
	cb.ax.tick_params(labelsize=fontsize-6)
	cb.set_label(r'log(M$_*$) [log(M$_{\odot}$)]', fontsize=fontsize, rotation=270, labelpad=labelpad)
	#plt.show()
	if (calibration == 'N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	if (calibration == 'O3N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
#Huang has stellar mass, hi mass, sfr, 

if (plot_to_make == 'luminosity_metallicity'):
	fig, ax = plt.subplots()
	#Mannucci_x = np.linspace(-20.96, -13.7, 100)
	#Mannucci_y = 8.96+0.31*(Mannucci_x-10)-0.23*(Mannucci_x-10)**2-0.017*(Mannucci_x-10)**3+0.046*(Mannucci_x-10)**4
	#plt.plot(Mannucci_x, Mannucci_y, color='red', linewidth=5.0, label='Mannucci Fit')

	fitting_SDSS_M_B = SDSS_M_B[SDSS_M_B>-20.96]
	fitting_SDSS_metallicity = SDSS_metallicity[SDSS_M_B>-20.96]
	fitting_SDSS_M_B = fitting_SDSS_M_B[fitting_SDSS_metallicity<9.1]
	fitting_SDSS_metallicity = fitting_SDSS_metallicity[fitting_SDSS_metallicity<9.1]
	fitting_SDSS_M_B = fitting_SDSS_M_B[fitting_SDSS_metallicity>7.6]
	fitting_SDSS_metallicity = fitting_SDSS_metallicity[fitting_SDSS_metallicity>7.6]
	final_fitting_SDSS_M_B = fitting_SDSS_M_B[fitting_SDSS_M_B<-13.5]
	final_fitting_SDSS_metallicity = fitting_SDSS_metallicity[fitting_SDSS_M_B<-13.5]
	fitting_SDSS_M_B = final_fitting_SDSS_M_B
	fitting_SDSS_metallicity = final_fitting_SDSS_metallicity
	medians_array = np.array([])
	medians_x_axis = np.array([])
	for bin_limit in M_B_bins_lower_limit:
		#np.logical_and
		#find median of the y-axis values for every galaxy in this bin
		in_bin = np.logical_and(fitting_SDSS_M_B>bin_limit,fitting_SDSS_M_B<bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0]))
		median = np.mean(fitting_SDSS_metallicity[in_bin])
		if np.isfinite(median):
			#print(median)
			medians_array = np.append(medians_array,median)
			medians_x_axis = np.append(medians_x_axis, bin_limit+((M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0])/2))
			print('M_B bin: '+str(round(bin_limit,2))+' mean value: '+str(round(median,2))+' scatter in bin: '+str(round(np.std(fitting_SDSS_metallicity[in_bin]-median),2))+' number of galaxies in bin: '+str(len(fitting_SDSS_metallicity[in_bin])))
	#medians_array = np.array([8.51819964804, 8.5881260704, 8.68400184466, 8.74884150215, 8.7835023075, 8.8081325197, 8.82203300249, 8.83091902361, 8.8185881194, 8.81050094938, 8.47969102404])
	#medians_x_axis = np.array([9.125, 9.375, 9.625, 9.875, 10.125, 10.375, 10.625, 10.875, 11.125, 11.375, 11.875])
	#print('medians_array: '+str(medians_array))
	#print('medians_x_axis: '+str(medians_x_axis))
	#plt.scatter(medians_x_axis, medians_array, marker='.', s=400, color='blue', label='median points')
	#fit_results, fit_error = optimization.curve_fit(fourth_order_fit, fitting_SDSS_M_B-10, fitting_SDSS_metallicity)
	#print('4th order fit attempt: '+str(fit_results))
	fit_results, fit_error = optimization.curve_fit(fourth_order_fit, medians_x_axis, medians_array, p0=[8.96,0.31,-0.23,-0.017,0.046])
	print('4th order fit attempt: '+str(fit_results))
	#fit_results = np.polyfit(medians_x_axis-10, medians_array, 4)
	#print('4th order fit attempt: '+str(fit_results))
	#fit_results = [8.76,0.31,-0.23,-0.017,0.046]
	#plt.scatter(fitting_SDSS_M_B, fitting_SDSS_metallicity, color='green', alpha=0.5)
	#Jimmy_y = fit_results[0]+fit_results[1]*(Mannucci_x-10)+fit_results[2]*(Mannucci_x-10)**2+fit_results[3]*(Mannucci_x-10)**3+fit_results[4]*(Mannucci_x-10)**4
	Mannucci_x = np.linspace(min(medians_x_axis), max(medians_x_axis), 100)
	Jimmy_y = fourth_order_fit(Mannucci_x, *fit_results)
	ax.scatter(medians_x_axis, medians_array, color='orange', s=500, linewidth='5', marker="_", label="Bin Means")
	ax.plot(Mannucci_x, Jimmy_y, color='green', linewidth=4.0, label='Fit to Means')
	n = len(fitting_SDSS_M_B)			   # number of samples     
	polynomial_fit = fit_results[0]+fit_results[1]*(fitting_SDSS_M_B-10)+fit_results[2]*(fitting_SDSS_M_B-10)**2+fit_results[3]*(fitting_SDSS_M_B-10)**3+fit_results[4]*(fitting_SDSS_M_B-10)**4
	residuals = fitting_SDSS_metallicity - polynomial_fit
	var_res = np.sum(residuals**2)/(n-2)
	sd_res = np.sqrt(var_res)
	#print('polynomial_fit: '+str(polynomial_fit))
	#print('fitting_SDSS_metallicity: '+str(fitting_SDSS_metallicity))
	#print('residuals: '+str(residuals))
	#print('max residuals: '+str(residuals[residuals<-1]))
	#print('max residuals M_B: '+str(fitting_SDSS_M_B[residuals<-1]))
	#print('max residuals metallicity: '+str(fitting_SDSS_metallicity[residuals<-1]))
	#print('max residuals fit metallicity: '+str(polynomial_fit[residuals<-1]))
	#print('var res: '+str(var_res))
	#print('n: '+str(n))
	print('Scatter AKA standard deviation: '+str(round(sd_res,2)))
	
	
	plt.minorticks_on()
	#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
	plt.text(-9.6, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
	ax.fill_between(np.array(range(24))-22, 9.1, 12.0, facecolor='yellow', alpha=0.25)
	ax.fill_between(np.array(range(24))-22, 7.0, 7.6, facecolor='yellow', alpha=0.25)
	h, x, y, p = ax.hist2d(SDSS_M_B, SDSS_metallicity, [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[-20.96, -13.7],[7.6, 9.1]])
	x = x[:-1]
	y = y[:-1]
	if (calibration == 'N2'):
		ax.scatter(SDSS_M_B[SDSS_M_B>-13.7], SDSS_metallicity[SDSS_M_B>-13.7], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (N2)', rasterized=True)
	if (calibration == 'O3N2'):
		ax.scatter(SDSS_M_B[SDSS_M_B>-13.7], SDSS_metallicity[SDSS_M_B>-13.7], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS (O3N2)', rasterized=True)
	for x_index, x_limit in enumerate(x):
		for y_index, y_limit in enumerate(y):
			if np.isnan(h[x_index,y_index]):
				for sdss_index, temp_M_B in enumerate(SDSS_M_B):
					if (temp_M_B > x_limit) and (temp_M_B < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
						ax.scatter(temp_M_B, SDSS_metallicity[sdss_index], color='gray', alpha=sdss_alpha_value)
	#ax.scatter(SDSS_M_B[SDSS_M_B>-13.7], SDSS_metallicity[SDSS_M_B>-13.7], color='gray', alpha=sdss_alpha_value, label="ALFALFA/SDSS (N2)", rasterized=True)
	#ax.scatter([0,0], [0,0], color='gray', alpha=sdss_alpha_value, label="ALFALFA/SDSS (N2)", rasterized=True)
	#plt.scatter(SDSS_M_B, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label="ALFALFA/SDSS (N2)", rasterized=True)
	if (calibration == 'N2'):
		ax.scatter(James_M_B, James_metallicity, s=40, color='purple', label="James+14 (N2)")
		ax.scatter(Berg_M_B, Berg_metallicity, s=40, color='red', label="Berg+12 (N2)")
		ax.scatter(SHIELD_M_B, SHIELD_metallicity, s=40, color='blue', label="SHIELD (N2)")
		#ax.scatter(Lee_M_B, Lee_metallicity, s=40, facecolors='none', edgecolors='green', label="Lee+06 (Literature)")
		ax.scatter(Jimmy_M_B,Jimmy_metallicity, s=400, color='black', marker='*', label="Jimmy+15 (N2)")
		ax.scatter(Saintonge_M_B,Saintonge_metallicity, s=40, facecolors='black', edgecolors='black', marker='D', label="Saintonge (N2)")
		ax.errorbar(Jimmy_M_B,Jimmy_metallicity,yerr=Jimmy_metallicity_error, xerr=Jimmy_M_B_error, linestyle="None", color='black')
	if (calibration == 'O3N2'):
		ax.scatter(James_M_B, James_metallicity, s=40, color='purple', label="James+14 (O3N2)")
		ax.scatter(Berg_M_B, Berg_metallicity, s=40, color='red', label="Berg+12 (O3N2)")
		ax.scatter(SHIELD_M_B, SHIELD_metallicity, s=40, color='blue', label="SHIELD (O3N2)")
		#ax.scatter(Lee_M_B, Lee_metallicity, s=40, facecolors='none', edgecolors='green', label="Lee+06 (Literature)")
		ax.scatter(Jimmy_M_B,Jimmy_metallicity, s=400, color='black', marker='*', label="Jimmy+15 (O3N2)")
		ax.scatter(Saintonge_M_B,Saintonge_metallicity, s=40, facecolors='black', edgecolors='black', marker='D', label="Saintonge (O3N2)")
		ax.errorbar(Jimmy_M_B,Jimmy_metallicity,yerr=Jimmy_metallicity_error, xerr=Jimmy_M_B_error, linestyle="None", color='black')
	#plt.plot(lum_x, metals, label="Skillman+89")
	plt.xlabel(r'M$_B$',fontsize=fontsize, labelpad=20)
	plt.ylabel('12+log(O/H)      ',fontsize=fontsize, labelpad=10)
	plt.ylim(7.5,9.5)
	plt.xlim(-10,-21.0)
	#plt.legend(loc="upper left", ncol=2)
	artist = plt.Rectangle((0, 0), 1, 1, fc="gray", ec=None, alpha = sdss_alpha_value)
	handles, labels = ax.get_legend_handles_labels()
	handles = handles+[artist]
	if (calibration == 'N2'):
		labels = labels+[u'ALFALFA/SDSS (N2)']
	if (calibration == 'O3N2'):
		labels = labels+[u'ALFALFA/SDSS (O3N2)']
	ax.legend(handles, labels, loc='upper left', ncol=2)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	if (calibration == 'N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	if (calibration == 'O3N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
	#plt.show()

if (plot_to_make == 'hi_mass_stellar_mass'):
	fig=plt.figure()
	#plt.scatter(alfalfa_sdss_catalog['hi_mass'], alfalfa_sdss_catalog['lgm_tot_p50'], color='gray', s=5, alpha=0.1)
	#plt.scatter(Jimmy_hi_mass, Jimmy_mass, color='black', s=40)
	#plt.xlim(5.5, 12)
	#plt.ylim(5.5, 12)
	#plt.xlabel('HI Mass')
	#plt.ylabel('Stellar Mass')
	stellar_mass_x = 5+np.array(range(200))*0.1
	hi_mass_y = stellar_mass_x
	plt.plot(stellar_mass_x, hi_mass_y, color='black', label='1:1 Relation')
	plt.scatter(SDSS_mass, SDSS_hi_mass, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	plt.scatter(huang_catalog['logM'], huang_catalog['logMHI'], color='blue', s=40, label='Huang+12')
	plt.scatter(SHIELD_mass, SHIELD_hi_mass, s=40, color='blue', label='SHIELD')
	plt.scatter(Saintonge_mass, Saintonge_hi_mass, facecolor='black', edgecolor='black', s=400, marker='*', label='Saintonge')
	plt.scatter(Jimmy_mass, Jimmy_hi_mass, color='black', s=400, marker='*', label='Jimmy+15')
	plt.errorbar(Jimmy_mass,Jimmy_hi_mass, xerr=Jimmy_mass_error, yerr=Jimmy_hi_mass_error, linestyle="None", color='black')
	plt.ylim(6.0, 11)
	plt.xlim(6.0, 11)
	plt.minorticks_on()
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	plt.ylabel('HI Mass', fontsize=fontsize, labelpad=20)
	plt.xlabel('Stellar Mass', fontsize=fontsize, labelpad=10)
	plt.legend(loc="upper left")
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#plt.show()

if (plot_to_make == 'hi_mass_density_sfr'):
	r90 = np.array([12.58742, 1.342661, 12.00269, 17.82006, 29.33742, 12.00269, 15.70734, 26.1804, 15.42202, 14.39525, 9.693284])
	r50 = np.array([6.080357, 0.6931014, 4.746222, 8.17783, 12.36714, 4.746222, 6.273112, 11.31596, 6.666181, 5.572609, 4.79457])
	radius_rad = r90*4.84813681E-6
	radius_Mpc = distance*np.tan(radius_rad)
	radius_Mpc = radius_Mpc*1
	volume = (4)*np.pi*(radius_Mpc**3)
	density = Jimmy_hi_mass/volume 
	sdss_radius_rad = alfalfa_sdss_catalog['petroR90_r']*4.84813681E-6
	sdss_redshift = np.array(alfalfa_sdss_catalog['z'])
	sdss_Dp = (speed_of_light*sdss_redshift/H0)*(1-(((1+q0)/2)*sdss_redshift))
	sdss_radius_Mpc = sdss_Dp*np.tan(sdss_radius_rad)
	sdss_volume = (4)*np.pi*sdss_radius_Mpc**3
	sdss_density = alfalfa_sdss_catalog['hi_mass']/sdss_volume
	#print(radius_Mpc)
	#print(volume)
	
	
	plt.scatter(np.log10(sdss_density), alfalfa_sdss_catalog['sfr_tot_p50'], color='gray', alpha=sdss_alpha_value, label="ALFALFA/SDSS", rasterized=True)
	plt.scatter(np.log10(density), Jimmy_sfr, color='black', label="Jimmy")
	plt.xlim(3, 16)
	plt.ylim(-6, 3)
	volume = (4)*np.pi*(radius_Mpc**2)
	sdss_volume = (4)*np.pi*sdss_radius_Mpc**2
	density = Jimmy_hi_mass/volume 
	sdss_density = alfalfa_sdss_catalog['hi_mass']/sdss_volume
	plt.scatter(np.log10(sdss_density), alfalfa_sdss_catalog['sfr_tot_p50']-np.log10(sdss_volume), color='orange', alpha=sdss_alpha_value, label="ALFALFA/SDSS", rasterized=True)
	plt.scatter(np.log10(density), Jimmy_sfr-np.log10(volume), color='green', label="Jimmy")
	plt.xlim(3, 16)
	plt.ylim(-6, 3)
	plt.ylabel('SFR')
	plt.xlabel('HI Mass/Volume')
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#plt.show()
	#print(radius_Mpc)
	#distance = redshift*1.302773E26
	#radius_m = distance*tan(radius_rad)
	
	#convert arcseconds above into pc or something
	#convert pc into a density, either surface or volume
	#r90 = np.array([12.58742, 1.342661, 12.00269, 17.82006, 29.33742, 12.00269, 15.70734, 26.1804, 15.42202, 14.39525, 9.693284])
	#r50 = np.array([6.080357, 0.6931014, 4.746222, 8.17783, 12.36714, 4.746222, 6.273112, 11.31596, 6.666181, 5.572609, 4.79457])
	#divide stellar mass by volume or surface density
	#plot density against SFR for first of all

if (plot_to_make == 'table_of_everything'):
	print('attempted distance error: '+str(percent_error*Dp*1.1))
	print('based on cz: '+str(z_array_error*speed_of_light/H0))
	#print('units check! '+str(z_array_error*speed_of_light))
	distance_error = np.sqrt((163**2)+(z_array_error*speed_of_light)**2)/H0
	print('convolved with systematic errors: '+str(distance_error))
	hi_ra = np.array(['09:08:36.4', '10:28:54.5', '11:34:54.9', '12:32:48.6', '12:36:33.5', '12:38:13.7', '12:46:03.0', '12:46:16.0', '12:59:41.8', '12:03:26.2', '12:50:04.0'])
	hi_dec = np.array(['+05:18:00', '+09:51:30', '+11:01:27', '+07:48:00', '+08:03:11', '+07:00:26', '+08:27:50', '+10:11:43', '+10:43:52', '+13:27:56', '+06:50:44'])
	oc_ra = np.array(['09:08:36.5', '10:28:55.8', '11:34:53.4', '12:32:47.0', '12:36:34.9', '12:38:15.5', '12:46:04.4', '12:46:15.3', '12:59:41.9', '12:03:26.3', '12:50:04.2'])
	oc_dec = np.array(['+05:17:32', '+09:51:47', '+11:01:10', '+07:47:58', '+08:03:17', '+06:59:40', '+08:28:34', '+10:12:20', '+10:43:40', '+13:27:34', '+06:50:51'])
	for index, galaxy in enumerate(names):
		#print(str(galaxy)+' & '+oc_ra[index]+' & '+oc_dec[index]+' & '+str(round(z_array[index],4))+' & '+str(Jimmy_hi_mass[index])+' & '+str(round(Jimmy_mass[index],2))+' & '+str(round(Jimmy_metallicity[index],2))+' & '+str(round(10**jimmy_sfr[index],3))+' & '+str(round(E_B_minus_V_array[index],3))+' \\\\')
		#print(Jimmy_sfr)
		#print(SFR)
		##print(str(galaxy)+' & '+oc_ra[index]+' & '+oc_dec[index]+' & '+str(int(round(speed_of_light*z_array[index],0)))+' $\pm$ '+str(int(round(speed_of_light*abs(z_array_error[index]),0)))+' & '+str(round(Jimmy_hi_mass[index],2))+' $\pm$ '+str(round(Jimmy_hi_mass_error[index],2))+' & '+str(round(Jimmy_mass[index],2))+' $\pm$ '+str(round(Jimmy_mass_error[index],2))+' & '+str(round(Jimmy_M_B[index],2))+' $\pm$ '+str(round(Jimmy_M_B_error[index],2))+' & '+str(round(Jimmy_metallicity[index],2))+' $\pm$ '+str(round(Jimmy_metallicity_error[index],2))+' & '+str(round(Jimmy_sfr[index],2))+' $\pm$ '+str(round(Jimmy_sfr_error[index],2))+' \\\\')
		print(str(galaxy)+' & '+oc_ra[index]+' & '+oc_dec[index]+' & '+str(Dp[index])+' $\pm$ '+str(round(distance_error[index],2))+' & '+str(round(Jimmy_hi_mass[index],2))+' $\pm$ '+str(round(Jimmy_hi_mass_error[index],2))+' & '+str(round(Jimmy_mass[index],2))+' $\pm$ '+str(round(Jimmy_mass_error[index],2))+' & '+str(round(Jimmy_M_B[index],2))+' $\pm$ '+str(round(Jimmy_M_B_error[index],2))+' & '+str(round(Jimmy_metallicity[index],2))+' $\pm$ '+str(round(Jimmy_metallicity_error[index],2))+' & '+str(round(Jimmy_sfr[index],2))+' $\pm$ '+str(round(Jimmy_sfr_error[index],2))+' \\\\')
	#for ra in oc_ra:
	#	oc_ra_degrees = int(ra[:2])*15
	#	oc_ra_degrees = oc_ra_degrees + int(ra[3:5])/7.5i
	#	oc_ra_degrees = oc_ra_degrees + float(ra[6:10])/3600
	#	print(oc_ra_degrees)
	#print(' ')
	#for dec in oc_dec:
	#	oc_dec_degrees = int(dec[1:3])*15
	#	oc_dec_degrees = oc_dec_degrees + int(dec[4:6])/7.5
	#	oc_dec_degrees = oc_dec_degrees + int(dec[7:9])/3600
	#	print(oc_dec_degrees)
	#print(str(oc_ra_degrees) + ' ' + str(oc_dec_degrees))
	#clean_hi_mass = (alfalfa_sdss_catalog['hi_mass'])[alfalfa_sdss_catalog['lgm_tot_p50'] > 0]
	#clean_stellar_mass = (alfalfa_sdss_catalog['lgm_tot_p50'])[alfalfa_sdss_catalog['lgm_tot_p50'] > 0]
	#clean_metallicity = (alfalfa_sdss_catalog['PP04_Metallicity'])[alfalfa_sdss_catalog['lgm_tot_p50'] > 0]
	#clean_sfr = (alfalfa_sdss_catalog['sfr_tot_p50'])[alfalfa_sdss_catalog['lgm_tot_p50'] > 0]
	#super_clean_hi_mass = clean_hi_mass[np.isfinite(clean_metallicity)]
	#super_clean_stellar_mass = clean_stellar_mass[np.isfinite(clean_metallicity)]
	#super_clean_metallicity = clean_metallicity[np.isfinite(clean_metallicity)]
	#super_clean_sfr = clean_sfr[np.isfinite(clean_metallicity)]
	#from matplotlib import cm
	#temp_colors = super_clean_sfr-super_clean_stellar_mass
	#color_zero_to_one = (temp_colors +11)/(3)
	#for color in temp_colors:
	#	color_zero_to_one = (color +11)/(3)
	#	print(color_zero_to_one)
	#	print(plt.cm.jet_r(color_zero_to_one))
	#output = open('/Users/jimmy/Downloads/json2.txt', "w")
	#output.write('{\n')
	#output.write('\t"hi_mass": [\n')
	#output.write('\t"red": [\n')
	#for color in color_zero_to_one:
	#	output.write('\t\t'+str(plt.cm.jet_r(color)[0])+',\n')
	#for item in super_clean_hi_mass:
	#	output.write('\t\t'+str(item)+',\n')
	#output.write('\t],\n')
	#output.write('\t"stellar_mass": [\n')
	#output.write('\t"green": [\n')
	#for color in color_zero_to_one:
	#	output.write('\t\t'+str(plt.cm.jet_r(color)[1])+',\n')
	#for item in super_clean_stellar_mass:
	#	output.write('\t\t'+str(item)+',\n')
	#output.write('\t],\n')
	#output.write('\t"metallicity": [\n')
	#output.write('\t"blue": [\n')
	#for color in color_zero_to_one:
	#	output.write('\t\t'+str(plt.cm.jet_r(color)[2])+',\n')
	#for item in super_clean_metallicity:
	#	output.write('\t\t'+str(item)+',\n')
	#output.write('\t],\n')
	#output.write('\t"star_formation_rate": [\n')
	#for item in super_clean_sfr:
	#	output.write('\t\t'+str(item)+',\n')
	#output.write('\t]\n')
	#output.write('}')
	

if (plot_to_make == 'sfr_comparison'):
	FAST_sfr = np.array([-3.22, -3.45, -2.60, -2.02, -2.50, -5.11, -2.37, -2.84, -1.33, -3.03, -4.46])
	H_alpha_sfr = Jimmy_sfr
	plt.scatter(H_alpha_sfr, FAST_sfr, color='blue')
	x_axis = np.array(range(100))*0.1-6
	y_axis = np.array(range(100))*0.1-6
	plt.plot(x_axis,y_axis, color='black')
	plt.xlim(-5.5,1.0)
	plt.ylim(-5.5,1.0)
	plt.xlabel(r'H$\alpha$ Derrived SFR',fontsize=fontsize, labelpad=20)
	plt.ylabel('FAST (SED fitting) Derrived SFR',fontsize=fontsize, labelpad=20)
	plt.tick_params(axis='y', labelsize=16)
	plt.tick_params(axis='x', labelsize=16)
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	#plt.show()
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')

if (plot_to_make == 'mass_comparison'):
	#light_to_mass = np.array(light_to_mass)
	#print(FAST_mass)
	#print(bell_03_mass)
	#print(west_10_mass)
	#print(converted_HI_mass_low)
	#print(light_to_mass)
	plt.scatter(light_to_mass, Jimmy_mass, color='blue')
	#print(light_to_mass)
	#print(FAST_mass)
	x_axis = np.array(range(100))*0.1
	y_axis = np.array(range(100))*0.1
	plt.plot(x_axis,y_axis, color='black')
	plt.xlim(3,8.0)
	plt.ylim(3,8.0)
	plt.xlabel(r'r Band Luminosity Mass (M/L)',fontsize=fontsize, labelpad=20)
	plt.ylabel('Bell 03/West 10 Derrived Mass',fontsize=fontsize, labelpad=20)
	plt.tick_params(axis='y', labelsize=16)
	plt.tick_params(axis='x', labelsize=16)
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	#plt.show()
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')

if (plot_to_make == 'field_of_views'):
	#reorder = np.array([4, 3, 6, 7, 8, 1, 2, 10, 5, 0, 9])
	reorder = np.array([0,1,2,3,4,5,6,7,8,9,10])
	from astropy import wcs
	from astropy.io import fits
	fig = plt.figure(figsize=(12, 9))
	#fig = plt.gcf()
	filenames = np.array(['shifted_frame-r-003015-4-0154.fits', 'shifted_frame-r-003031-6-0249.fits', 'shifted_frame-r-003804-1-0108.fits', 'shifted_frame-r-003842-4-0048.fits', 'shifted_frame-r-003841-5-0162.fits', 'shifted_frame-r-003842-2-0057.fits', 'shifted_frame-r-003841-6-0177.fits', 'shifted_frame-r-003903-4-0069.fits', 'shifted_frame-r-003031-6-0496.fits', 'shifted_frame-r-003804-6-0155.fits', 'shifted_frame-r-006005-2-0209.fits'])
	fiber_ra = np.array([137.15227, 157.23223, 0.0, 188.19489, 189.14573, 189.56452, 191.51849, 191.56374, 194.92518, 0.0, 0.0])
	fiber_dec = np.array([5.2908198, 9.8616, 0.0, 7.7990335, 8.05475, 6.99407, 8.47625, 10.205598, 10.726815, 0.0, 0.0])
	#petro_ra = np.array([137.15225, 157.232226, 173.721301, 188.19489134, 189.14572811, 189.56450, 191.51838, 191.56389593, 194.92517977, 180.85952, 192.51792])
	#petro_dec = np.array([5.29083, 9.86159997, 11.02118606, 7.79903354, 8.05474689, 6.99408, 8.47646, 10.20456917, 10.72681536, 13.45950, 6.84748])
	petro_rad = np.array([12.80, 20.55, 10.27, 17.78, 28.88, 6.65, 13.55, 25.76, 14.46, 27.93, 9.37])
	radio_coords = np.array(['090836.4+051800', '102854.5+095130', '113454.9+110127', '123248.6+074800', '123633.5+080311', '123813.7+070026', '124603.0+082750 ', '124616.0+101143', '125941.8+104352', '120326.2+132756', '125004.0+065044'])
	hi_ra = np.array(['09:08:36.4', '10:28:54.5', '11:34:54.9', '12:32:48.6', '12:36:33.5', '12:38:13.7', '12:46:03.0', '12:46:16.0', '12:59:41.8', '12:03:26.2', '12:50:04.0'])
	hi_dec = np.array(['+05:18:00', '+09:51:30', '+11:01:27', '+07:48:00', '+08:03:11', '+07:00:26', '+08:27:50', '+10:11:43', '+10:43:52', '+13:27:56', '+06:50:44'])
	oc_ra = np.array(['09:08:36.5', '10:28:55.8', '11:34:53.4', '12:32:47.0', '12:36:34.9', '12:38:15.5', '12:46:04.4', '12:46:15.3', '12:59:41.9', '12:03:26.3', '12:50:04.2'])
	oc_dec = np.array(['+05:17:32', '+09:51:47', '+11:01:10', '+07:47:58', '+08:03:17', '+06:59:40', '+08:28:34', '+10:12:20', '+10:43:40', '+13:27:34', '+06:50:51'])
	petro_ra = np.array([])
	petro_dec = np.array([])
	for coords_string in oc_ra:
		ra = (float(coords_string[0:2])*15)+float(coords_string[3:5])/4+float(coords_string[6:10])/240
		petro_ra = np.append(petro_ra, ra)
	for coords_string in oc_dec:
		dec = float(coords_string[1:3])+float(coords_string[4:6])/60+float(coords_string[7:9])/3600
		petro_dec = np.append(petro_dec, dec)
	radio_ra = np.array([])
	radio_dec = np.array([])
	for coords_string in hi_ra:
		ra = (float(coords_string[0:2])*15)+float(coords_string[3:5])/4+float(coords_string[6:10])/240
		radio_ra = np.append(radio_ra, ra)
	for coords_string in hi_dec:
		dec = float(coords_string[1:3])+float(coords_string[4:6])/60+float(coords_string[7:9])/3600
		radio_dec = np.append(radio_dec, dec)
	names = names[reorder]
	filenames = filenames[reorder]
	fiber_ra = fiber_ra[reorder]
	fiber_dec = fiber_dec[reorder]
	petro_ra = petro_ra[reorder]
	petro_dec = petro_dec[reorder]
	petro_rad = petro_rad[reorder]
	radio_coords = radio_coords[reorder]
	radio_ra = radio_ra[reorder]
	radio_dec = radio_dec[reorder]
	for index, galaxy in enumerate(names):
		filename = '/Astro/SDSS/'+galaxy+'/'+filenames[index]
		sp1 = fig.add_subplot(3,4,index+1)
		if index == 0:
			sp1.set_yticklabels([-32, ' ', -16, -8, 0, 8, 16, 24])
			plt.tick_params(axis='y', labelsize=20)
			#plt.ylabel(r'',fontsize=fontsize)
			sp1.set_xticklabels([])
		if index == 1:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([])
		if index == 2:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([])
		if index == 3:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([])
		if index == 4:
			sp1.set_yticklabels([-32, ' ', -16, -8, 0, 8, 16, 24])
			plt.tick_params(axis='y', labelsize=20)
			plt.ylabel(r'arcsec',fontsize=fontsize)
			sp1.set_xticklabels([])
		if index == 5:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([])
		if index == 6:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([])
		if index == 7:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([-32, -24, -16, -8, 0, 8, 16, ' '])
			plt.tick_params(axis='x', labelsize=20)
		if index == 8:
			sp1.set_yticklabels([-32, -24, -16, -8, 0, 8, 16, 24])
			plt.tick_params(axis='y', labelsize=20)
			sp1.set_xticklabels([-32, -24, -16, -8, 0, 8, 16, ' '])
			plt.tick_params(axis='x', labelsize=20)
		if index == 9:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([-32, -24, -16, -8, 0, 8, 16, ' '])
			plt.tick_params(axis='x', labelsize=20)
			plt.xlabel(r'arcsec',fontsize=fontsize)
		if index == 10:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([-32, -24, -16, -8, 0, 8, 16, 24])
			plt.tick_params(axis='x', labelsize=20)
		if index == 11:
			sp1.set_yticklabels([])
			sp1.set_xticklabels([])
		data, header = fits.getdata(HOME+filename, header=True)
		#hdu_list = pyfits.open(HOME+'/Astro/SDSS/AGC191702/shifted_frame-r-003015-4-0154.fits')
		#image_data = hdu_list[0].data
		w = wcs.WCS(header)
		x_pix, y_pix = w.wcs_world2pix(fiber_ra[index], fiber_dec[index], 1)
		plt.imshow(data, cmap='gray_r', vmin=0.0, vmax=0.1)
		#sdss fibers are 3 arcsec in diameter, that's 1.5" for the radius, the pixel resolution is 0.396127 arcsec/pixel	
		fiber_radius = 1.5/0.396127
		fiber_circle=plt.Circle((x_pix,y_pix),fiber_radius,color='red',linewidth=3.0, fill=False)
		petro_radius = petro_rad[index]/0.396127
		x_pix, y_pix = w.wcs_world2pix(petro_ra[index], petro_dec[index], 1)
		vimos_ifu_fov_side = 54/0.396127
		plt.xlim(x_pix-68,x_pix+68)
		plt.ylim(y_pix-68,y_pix+68)
		plt.text(x_pix-60, y_pix-60, galaxy, color='purple', fontsize=fontsize, fontweight='bold')
		petro_circle=plt.Circle((x_pix,y_pix),petro_radius, color='blue', linewidth=3.0, fill=False)
		radio_radius = 20/0.396127
		ifu_side = 24/0.396127
		#coords_string = radio_coords[index]
		#ra = (float(coords_string[0:2])*15)+float(coords_string[2:4])/4+float(coords_string[4:8])/240
		#dec = float(coords_string[9:11])+float(coords_string[11:13])/60+float(coords_string[13:15])/3600
		#9,8,10
		#RA: 137.0417
		#DEC: 9.1361
		#x_pix, y_pix = w.wcs_world2pix(radio_ra[index], radio_dec[index], 1)
		#radio_circle=plt.Circle((x_pix,y_pix),radio_radius, color='green', linewidth=3.0, fill=False)
		x_pix, y_pix = w.wcs_world2pix(petro_ra[index], petro_dec[index], 1)
		if galaxy == 'AGC191702':
			y_pix = y_pix+20
			x_pix = x_pix+5
		ifu_footprint = plt.Rectangle((x_pix-(ifu_side/2), y_pix-(ifu_side/2)), ifu_side, ifu_side, fc='None', edgecolor='green', linewidth=2)
		fig.gca().add_artist(fiber_circle)
		fig.gca().add_artist(petro_circle)
		#fig.gca().add_artist(radio_circle)
		fig.gca().add_artist(ifu_footprint)
		
		#plt.colorbar()
	sp1 = fig.add_subplot(3,4,12)
	sp1.set_yticklabels([])
	sp1.set_xticklabels([])
	#plt.plot([0,10],[0,10], color='black')
	plt.text(0.05, 0.55, "SDSS Fiber Radius", color='red', fontsize=16, fontweight='bold')
	plt.text(0.05, 0.45, "SDSS Petrosian Radius", color='blue', fontsize=16, fontweight='bold')
	plt.text(0.05, 0.35, "IFU Analysis Footprint", color='green', fontsize=15)
	plt.tick_params(axis='y', left='off', right='off')
	plt.tick_params(axis='x', bottom='off', top='off')
	plt.subplots_adjust(bottom=0.07)
	plt.subplots_adjust(left=0.07)
	plt.subplots_adjust(right=0.95)
	plt.subplots_adjust(top=0.95)
	plt.subplots_adjust(wspace=0.0)
	plt.subplots_adjust(hspace=0.0)
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#plt.show()

if (plot_to_make == 'fmr_hi_mass') or (plot_to_make == 'fmr_lzr_hi_mass'):
	fig, ax = plt.subplots()
	if (plot_to_make == 'fmr_hi_mass'):
		weight = 0.32
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		weight = -0.50
	# define the colormap
	cmap = plt.cm.jet
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# force the first color entry to be grey
	#cmaplist[0] = (.5,.5,.5,1.0)
	#cmaplist[0] = (0.0,0.0,0.0,1.0)
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	# define the bins and normalize
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	#ax.scatter(SDSS_mass, SDSS_metallicity, c=SDSS_hi_mass, alpha=0.5, cmap=cmap, norm=norm, label="ALFALFA", rasterized=True)
	#ax.scatter(SDSS_mass-weight*(SDSS_hi_mass-9.80), SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	if (plot_to_make == 'fmr_hi_mass'):
		h, x, y, p = ax.hist2d(SDSS_mass-weight*(SDSS_hi_mass-9.80), SDSS_metallicity, bins = [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[7.75, 11],[7.6, 9.1]])
		x = x[:-1]
		y = y[:-1]
		ax.scatter(SDSS_mass[SDSS_mass<7.75], SDSS_metallicity[SDSS_mass<7.75], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		for x_index, x_limit in enumerate(x):
			for y_index, y_limit in enumerate(y):
				if np.isnan(h[x_index,y_index]):
					for sdss_index, temp_mass in enumerate(SDSS_mass):
						if (temp_mass-weight*(SDSS_hi_mass[sdss_index]-9.80) > x_limit) and (temp_mass-weight*(SDSS_hi_mass[sdss_index]-9.80) < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
							ax.scatter(temp_mass-weight*(SDSS_hi_mass[sdss_index]-9.80), SDSS_metallicity[sdss_index], color='gray', marker='o', alpha=sdss_alpha_value)
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		h, x, y, p = ax.hist2d(SDSS_M_B-weight*(SDSS_hi_mass-9.80), SDSS_metallicity, bins = [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[-20.96, -13.7],[7.6, 9.1]])
		x = x[:-1]
		y = y[:-1]
		ax.scatter(SDSS_M_B[SDSS_M_B>-13.7], SDSS_metallicity[SDSS_M_B>-13.7], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		for x_index, x_limit in enumerate(x):
			for y_index, y_limit in enumerate(y):
				if np.isnan(h[x_index,y_index]):
					for sdss_index, temp_M_B in enumerate(SDSS_M_B):
						if (temp_M_B-weight*(SDSS_hi_mass[sdss_index]-9.80) > x_limit) and (temp_M_B-weight*(SDSS_hi_mass[sdss_index]-9.80) < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
							ax.scatter(temp_M_B-weight*(SDSS_hi_mass[sdss_index]-9.80), SDSS_metallicity[sdss_index], color='gray', marker='o', alpha=sdss_alpha_value)
	#plt.scatter(SDSS_mass, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	#ax.scatter([0], [0], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	
		
	#create bins and find means within them.
	for hi_mass_bin_index, hi_bin_limit in enumerate(hi_mass_bins_lower_limit):
		temp_sdss_metallicity_array = np.array([])
		temp_sdss_mass_array = np.array([])
		temp_sdss_M_B_array = np.array([])
		temp_sdss_hi_mass_array = np.array([])
		for index, galaxy_hi_mass in enumerate(SDSS_hi_mass):
			if (galaxy_hi_mass > hi_bin_limit) and (galaxy_hi_mass < (hi_bin_limit+(hi_mass_bins_lower_limit[1]-hi_mass_bins_lower_limit[0]))):
				temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,SDSS_metallicity[index])
				temp_sdss_mass_array = np.append(temp_sdss_mass_array, SDSS_mass[index])
				temp_sdss_M_B_array = np.append(temp_sdss_M_B_array, SDSS_M_B[index])
				temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, galaxy_hi_mass)
		temp_SDSS_mass = temp_sdss_mass_array
		temp_SDSS_M_B = temp_sdss_M_B_array
		temp_SDSS_metallicity = temp_sdss_metallicity_array
		temp_SDSS_hi_mass = temp_sdss_hi_mass_array
		
		temp_median_x_axis = np.array([])
		temp_median_y_axis = np.array([])
		temp_median_z_axis = np.array([])
		if (plot_to_make == 'fmr_hi_mass'):
			for bin_limit in mass_bins_lower_limit:
				#print('bin_limit'+str(bin_limit))
				temp_sdss_metallicity_array = np.array([])
				temp_sdss_mass_array = np.array([])
				temp_sdss_hi_mass_array = np.array([])
				for index, galaxy_mass in enumerate(temp_SDSS_mass):
					if (galaxy_mass > bin_limit) and (galaxy_mass < (bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))):
						temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
						temp_sdss_mass_array = np.append(temp_sdss_mass_array, temp_SDSS_mass[index])
						temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, temp_SDSS_hi_mass[index])
				if len(temp_sdss_metallicity_array) > binning_cut:
					#print(temp_sdss_metallicity_array)
					temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
					temp_median_x_axis = np.append(temp_median_x_axis, bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])/2)
					temp_median_z_axis = np.append(temp_median_z_axis, np.mean(temp_sdss_hi_mass_array))
			#print('temp x-axis line: '+str(temp_median_x_axis))
			#print('temp y-axis line: '+str(temp_median_y_axis))
		if (plot_to_make == 'fmr_lzr_hi_mass'):
			for bin_limit in M_B_bins_lower_limit:
				#print('bin_limit'+str(bin_limit))
				temp_sdss_metallicity_array = np.array([])
				temp_sdss_M_B_array = np.array([])
				temp_sdss_hi_mass_array = np.array([])
				for index, galaxy_M_B in enumerate(temp_SDSS_M_B):
					if (galaxy_M_B > bin_limit) and (galaxy_M_B < (bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0]))):
						temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
						temp_sdss_M_B_array = np.append(temp_sdss_M_B_array, temp_SDSS_M_B[index])
						temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, temp_SDSS_hi_mass[index])
				if len(temp_sdss_metallicity_array) > binning_cut:
					#print(temp_sdss_metallicity_array)
					temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
					temp_median_x_axis = np.append(temp_median_x_axis, bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0])/2)
					temp_median_z_axis = np.append(temp_median_z_axis, np.mean(temp_sdss_hi_mass_array))
			#print('temp x-axis line: '+str(temp_median_x_axis))
			#print('temp y-axis line: '+str(temp_median_y_axis))
		color_value = (hi_bin_limit-hi_mass_bins_lower_limit[0])/(hi_mass_bins_lower_limit[-2]-hi_mass_bins_lower_limit[0])
		ax.plot(temp_median_x_axis-weight*(temp_median_z_axis-9.80), temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
					

	temp_metallicity_array = np.array([])
	temp_mass_array = np.array([])
	temp_M_B_array = np.array([])
	temp_metallicity_array_error = np.array([])
	temp_hi_mass_array = np.array([])
	for index, galaxy_hi_mass in enumerate(Jimmy_hi_mass):
		if (galaxy_hi_mass < hi_mass_bins_lower_limit[1]) and (galaxy_hi_mass > hi_mass_bins_lower_limit[0]):
			temp_hi_mass_array = np.append(temp_hi_mass_array, Jimmy_hi_mass[index])
			temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
			temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
			temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
			temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
	#for index, galaxy_hi_mass in enumerate(Saintonge_hi_mass):
	#	if (Saintonge_hi_mass[index] < hi_mass_bins_lower_limit[2]) and (Saintonge_hi_mass[index] > hi_mass_bins_lower_limit[1]): #and (Saintonge_metallicity[index] > 8.0)
	#		temp_hi_mass_array = np.append(temp_hi_mass_array, Saintonge_hi_mass[index])	
	#		temp_mass_array = np.append(temp_mass_array, Saintonge_mass[index]	)
	#		temp_metallicity_array = np.append(temp_metallicity_array, Saintonge_metallicity[index])
	#for index, galaxy_hi_mass in enumerate(SDSS_hi_mass):
	#	if (galaxy_hi_mass < -2.00) and (galaxy_hi_mass > -2.75):
	#		temp_hi_mass_array = np.append(temp_hi_mass_array, SDSS_hi_mass[index])
	#		temp_mass_array = np.append(temp_mass_array, SDSS_mass[index])
	#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
	#print(temp_metallicity_array)
	#print(temp_hi_mass_array)
	#moustakas_guess = np.array([8.901, 8.798, 0.640])
	#fit_results, fit_error = optimization.curve_fit(moustakas, temp_mass_array, temp_metallicity_array, p0 = moustakas_guess)
	if (plot_to_make == 'fmr_hi_mass'):
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array-weight*(temp_hi_mass_array-9.80), temp_metallicity_array, p0=[4.0, 0.5])
		x_fit = np.linspace(min(temp_mass_array-weight*(temp_hi_mass_array-9.80)), max(temp_mass_array-weight*(temp_hi_mass_array-9.80)), 100)
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array-weight*(temp_hi_mass_array-9.80), temp_metallicity_array, p0=[4.0, 0.5])
		x_fit = np.linspace(min(temp_M_B_array-weight*(temp_hi_mass_array-9.80)), max(temp_M_B_array-weight*(temp_hi_mass_array-9.80)), 100)
	#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
	y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
	#y_fit = linear_fit(x_fit, 4, 0.5)
	color_value = (sfr_bins_lower_limit[0]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
	ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
	if (plot_to_make == 'fmr_hi_mass'):
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_mass_array-weight*(temp_hi_mass_array-9.80))))
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_M_B_array-weight*(temp_hi_mass_array-9.80))))
	var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
	sd_res = np.sqrt(var_res)
	print('standard deviation: '+str(sd_res))
	y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
	y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
	ax.fill_between(x_fit, y_fit_upper, y_fit_lower, facecolor=cmap(color_value), alpha=0.25)
	#ax.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
	if (plot_to_make == 'fmr_hi_mass'):
		ax.errorbar(temp_mass_array-weight*(temp_hi_mass_array-9.80),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		ax.errorbar(temp_M_B_array-weight*(temp_hi_mass_array-9.80),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
		
	temp_metallicity_array = np.array([])
	temp_mass_array = np.array([])
	temp_M_B_array = np.array([])
	temp_metallicity_array_error = np.array([])
	temp_hi_mass_array = np.array([])
	for index, galaxy_hi_mass in enumerate(Jimmy_hi_mass):
		if (galaxy_hi_mass < hi_mass_bins_lower_limit[2]) and (galaxy_hi_mass > hi_mass_bins_lower_limit[1]):
			temp_hi_mass_array = np.append(temp_hi_mass_array, Jimmy_hi_mass[index])
			temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
			temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
			temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])				
			temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
	#for index, galaxy_hi_mass in enumerate(Saintonge_hi_mass)
	#	if (Saintonge_hi_mass[index] < hi_mass_bins_lower_limit[2]) and (Saintonge_hi_mass[index] > hi_mass_bins_lower_limit[1]): #and (Saintonge_metallicity[index] > 8.0)
	#		temp_hi_mass_array = np.append(temp_hi_mass_array, Saintonge_hi_mass[index])	
	#		temp_mass_array = np.append(temp_mass_array, Saintonge_mass[index]	)
	#		temp_metallicity_array = np.append(temp_metallicity_array, Saintonge_metallicity[index])
	#for index, galaxy_hi_mass in enumerate(SDSS_hi_mass):
	#	if (galaxy_hi_mass < -2.00) and (galaxy_hi_mass > -2.75):
	#		temp_hi_mass_array = np.append(temp_hi_mass_array, SDSS_hi_mass[index])
	#		temp_mass_array = np.append(temp_mass_array, SDSS_mass[index])
	#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
	#print(temp_metallicity_array)
	#moustakas_guess = np.array([8.901, 8.798, 0.640])
	##linear_guess = np.array([6.24683786, 0.26661621])
	#fit_results, fit_error = optimization.curve_fit(moustakas, temp_mass_array, temp_metallicity_array, p0 = moustakas_guess)
	if (plot_to_make == 'fmr_hi_mass'):
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array-weight*(temp_hi_mass_array-9.80), temp_metallicity_array, p0=[4.0, 0.5])
		x_fit = np.linspace(min(temp_mass_array-weight*(temp_hi_mass_array-9.80)), max(temp_mass_array-weight*(temp_hi_mass_array-9.80)), 100)
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array-weight*(temp_hi_mass_array-9.80), temp_metallicity_array, p0=[4.0, 0.5])
		x_fit = np.linspace(min(temp_M_B_array-weight*(temp_hi_mass_array-9.80)), max(temp_M_B_array-weight*(temp_hi_mass_array-9.80)), 100)
	color_value = (sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
	#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
	y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
	ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
	if (plot_to_make == 'fmr_hi_mass'):
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_mass_array-weight*(temp_hi_mass_array-9.80))))
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_M_B_array-weight*(temp_hi_mass_array-9.80))))
		print('residuals: '+str(residuals))
	var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
	sd_res = np.sqrt(var_res)
	print('standard deviation: '+str(sd_res))
	y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
	y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
	ax.fill_between(x_fit, y_fit_upper, y_fit_lower, facecolor=cmap(color_value), alpha=0.25)
	#ax.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
	if (plot_to_make == 'fmr_hi_mass'):
		ax.errorbar(temp_mass_array-weight*(temp_hi_mass_array-9.80),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		ax.errorbar(temp_M_B_array-weight*(temp_hi_mass_array-9.80),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
		print('temp_M_B_array: '+str(temp_M_B_array))

	#x_fit = np.linspace(5.50, 9.50, 100)
	#y_fit = linear_fit(x_fit, 6.40820096,  0.22217241)
	#color_value = 0.0
	#plt.plot(x_fit-weight*(7.0625-9.8), y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
	#x_fit = np.linspace(6.50, 10.25, 100)
	#y_fit = linear_fit(x_fit, 6.72859729,  0.18665824)
	#color_value = 0.25
	#plt.plot(x_fit-weight*(7.6875-9.8), y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)


	#print(Jimmy_hi_mass)
	ax.plot([0,0], [0,0], color='black', label='SDSS Mean')
	if (plot_to_make == 'fmr_hi_mass'):
		s = ax.scatter(Jimmy_mass-weight*(Jimmy_hi_mass-9.80), Jimmy_metallicity, c=Jimmy_hi_mass, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		s = ax.scatter(Jimmy_M_B-weight*(Jimmy_hi_mass-9.80), Jimmy_metallicity, c=Jimmy_hi_mass, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
	#plt.scatter(Saintonge_mass-weight*(Saintonge_hi_mass-9.80), Saintonge_metallicity, facecolor='black', c=Saintonge_hi_mass, s=40, cmap=cmap, norm=norm, marker='D', label="Saintonge")	
	#plt.scatter(SHIELD_mass-weight*(SHIELD_mass-9.80), SHIELD_metallicity, c=SHIELD_hi_mass, s=40, cmap=cmap, norm=norm, label="SHIELD", color='black')
	#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
	if (plot_to_make == 'fmr_hi_mass'):
		ax.text(solar_x_axis-0.1, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
	#logmstar = 7.0+(np.arange(400)*0.01)
	#PP04_line = 32.1488-(8.51258*logmstar)+(0.976384*logmstar*logmstar)-(0.0359763*logmstar*logmstar*logmstar)
	#plt.plot(logmstar,PP04_line, color="black", linewidth=3.0, label="PP04 SDSS Relation")
	plt.ylabel('12+log(O/H)',fontsize=fontsize, labelpad=10)
	if (plot_to_make == 'fmr_hi_mass'):
		plt.xlabel(r'log(M$_{*}$)-'+str(weight)+'(log(HI Mass)-9.80)',fontsize=fontsize, labelpad=20)
		plt.xlim(6.0, 11.0)
	if (plot_to_make == 'fmr_lzr_hi_mass'):
		plt.xlabel(r'M$_{B}$+'+str(abs(weight))+'(log(HI Mass)-9.80)',fontsize=fontsize, labelpad=20)
		plt.xlim(-12, -22)
	plt.ylim(7.6, 9.3)	
	plt.minorticks_on()
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	#plt.legend(loc='upper left')
	artist = plt.Rectangle((0, 0), 1, 1, fc="gray", ec=None, alpha = sdss_alpha_value)
	handles, labels = ax.get_legend_handles_labels()
	handles = handles+[artist]
	labels = labels+[u'ALFALFA/SDSS']
	ax.legend(handles, labels, loc='upper left', ncol=2)
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	#ax2 = fig.add_axes([0.855, 0.175, 0.03, 0.775])
	#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
	#cb.ax.tick_params(labelsize=fontsize-6)
	#cb.set_label(r'log(HI Mass) [log(M$_{\odot}$)]', fontsize=fontsize, rotation=270)
	cb = plt.colorbar(s)
	#plt.ticks.set_fontsize(20)
	cb.ax.tick_params(labelsize=fontsize-6)
	cb.set_label(r'log(HI Mass) [log(M$_{\odot}$)]', fontsize=fontsize, rotation=270, labelpad=labelpad)
	if (calibration == 'N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	if (calibration == 'O3N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
	#plt.show()
	

if (plot_to_make == 'fmr_sfr') or (plot_to_make == 'fmr_lzr_sfr'):
	combine_lowest_bins = 'n'
	fig, ax = plt.subplots()
	if (plot_to_make == 'fmr_sfr'):
		weight = 0.28
	if (plot_to_make == 'fmr_lzr_sfr'):
		weight = -0.35
	#print(Jimmy_sfr)
	# define the colormap
	cmap = plt.cm.jet
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# force the first color entry to be grey
	#cmaplist[0] = (.5,.5,.5,1.0)
	if (combine_lowest_bins == 'y'):
		print(cmaplist[0])
		for cmap_index in range(70):
			#cmaplist[cmap_index] = (0.0,0.0,0.5,1.0)
			cmaplist[cmap_index] = (0.0, 0.48823529411764705, 1.0, 1.0)
		print(cmaplist[0])
		print(cmaplist[1])
		print(len(cmaplist))
	#print(cmaplist.shape)
	print('done printing cmaps')
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	cmap.set_under('k')
	#cmap.set_clim(-1.7, 1.0)
	# define the bins and normalize
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	#plt.scatter(SDSS_mass-weight*(SDSS_sfr), SDSS_metallicity, c=SDSS_sfr, alpha=0.5, cmap=cmap, norm=norm, label="ALFALFA", rasterized=True)
	#plt.scatter(SDSS_mass-weight*(SDSS_sfr), SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	#print(SDSS_mass-weight*(SDSS_sfr))
	#print(SDSS_sfr[np.isnan(SDSS_sfr)])
	#SDSS_sfr[np.isnan(SDSS_sfr)] = 0.0
	#print(SDSS_sfr[np.isnan(SDSS_sfr)])
	if (plot_to_make == 'fmr_sfr'):
		h, x, y, p = ax.hist2d(SDSS_mass-weight*(SDSS_sfr), SDSS_metallicity, bins = [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[7.75, 11],[7.6, 9.1]])
		x = x[:-1]
		y = y[:-1]
		ax.scatter(SDSS_mass[SDSS_mass<7.75]-weight*(SDSS_sfr[SDSS_mass<7.75]), SDSS_metallicity[SDSS_mass<7.75], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		for x_index, x_limit in enumerate(x):
			for y_index, y_limit in enumerate(y):
				if np.isnan(h[x_index,y_index]):
					for sdss_index, temp_mass in enumerate(SDSS_mass):
						if (temp_mass-weight*(SDSS_sfr[sdss_index]) > x_limit) and (temp_mass-weight*(SDSS_sfr[sdss_index]) < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
							ax.scatter(temp_mass-weight*(SDSS_sfr[sdss_index]), SDSS_metallicity[sdss_index], color='gray', marker='o', alpha=sdss_alpha_value)
		#plt.scatter(SDSS_mass, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	if (plot_to_make == 'fmr_lzr_sfr'):
		h, x, y, p = ax.hist2d(SDSS_M_B-weight*(SDSS_sfr), SDSS_metallicity, bins = [20, 20], cmap = plt.cm.Greys, cmin=5.0, vmin=-50.0, vmax=150.0, range=[[-20.96, -13.7],[7.6, 9.1]])
		x = x[:-1]
		y = y[:-1]
		ax.scatter(SDSS_M_B[SDSS_M_B>-13.7], SDSS_metallicity[SDSS_M_B>-13.7], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
		for x_index, x_limit in enumerate(x):
			for y_index, y_limit in enumerate(y):
				if np.isnan(h[x_index,y_index]):
					for sdss_index, temp_M_B in enumerate(SDSS_M_B):
						if (temp_M_B-weight*(SDSS_sfr[sdss_index]) > x_limit) and (temp_M_B-weight*(SDSS_sfr[sdss_index]) < x_limit+(x[1]-x[0])) and (SDSS_metallicity[sdss_index] > y_limit) and (SDSS_metallicity[sdss_index] < y_limit+(y[1]-y[0])):
							ax.scatter(temp_M_B-weight*(SDSS_sfr[sdss_index]), SDSS_metallicity[sdss_index], color='gray', marker='o', alpha=sdss_alpha_value)
		#plt.scatter(SDSS_mass, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	#ax.scatter([0], [0], color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	
		
	#create bins and find means within them.
	#sfr_bins_lower_limit = np.linspace(-3.5, 1.0, 6) #np.array([-1000.0, -5.00, -4.00, -3.00, -2.00, -1.00, 0.00, 1.00])
	#sfr_bins_lower_limit = np.array([-3.5, -1.7, -0.8, 0.1, 1.0])
	for sfr_bin_index, sfr_bin_limit in enumerate(sfr_bins_lower_limit):
		#print(sfr_bin_limit)
		temp_sdss_metallicity_array = np.array([])
		temp_sdss_mass_array = np.array([])
		temp_sdss_M_B_array = np.array([])
		temp_sdss_sfr_array = np.array([])
		for index, galaxy_sfr in enumerate(SDSS_sfr):
			if (galaxy_sfr > sfr_bin_limit) and (galaxy_sfr < (sfr_bin_limit+(sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0]))):
				temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,SDSS_metallicity[index])
				temp_sdss_mass_array = np.append(temp_sdss_mass_array, SDSS_mass[index])
				temp_sdss_M_B_array = np.append(temp_sdss_M_B_array, SDSS_M_B[index])
				temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, galaxy_sfr)
		temp_SDSS_mass = temp_sdss_mass_array
		temp_SDSS_M_B = temp_sdss_M_B_array
		temp_SDSS_metallicity = temp_sdss_metallicity_array
		temp_SDSS_sfr = temp_sdss_sfr_array
		
		temp_median_x_axis = np.array([])
		temp_median_y_axis = np.array([])
		temp_median_z_axis = np.array([])
		if (plot_to_make == 'fmr_sfr'):
			for bin_limit in mass_bins_lower_limit:
				#print('bin_limit'+str(bin_limit))
				temp_sdss_metallicity_array = np.array([])
				temp_sdss_mass_array = np.array([])
				temp_sdss_sfr_array = np.array([])
				for index, galaxy_mass in enumerate(temp_SDSS_mass):
					if (galaxy_mass > bin_limit) and (galaxy_mass < (bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))):
						temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
						temp_sdss_mass_array = np.append(temp_sdss_mass_array, temp_SDSS_mass[index])
						temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, temp_SDSS_sfr[index])
				if len(temp_sdss_metallicity_array) > binning_cut:
					temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
					temp_median_x_axis = np.append(temp_median_x_axis, bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])/2)
					temp_median_z_axis = np.append(temp_median_z_axis, np.mean(temp_sdss_sfr_array))
			color_value = (sfr_bin_limit-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
			ax.plot(temp_median_x_axis-weight*(temp_median_z_axis), temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
		if (plot_to_make == 'fmr_lzr_sfr'):
			for bin_limit in M_B_bins_lower_limit:
				#print('bin_limit'+str(bin_limit))
				temp_sdss_metallicity_array = np.array([])
				temp_sdss_M_B_array = np.array([])
				temp_sdss_sfr_array = np.array([])
				for index, galaxy_M_B in enumerate(temp_SDSS_M_B):
					if (galaxy_M_B > bin_limit) and (galaxy_M_B < (bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0]))):
						temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
						temp_sdss_M_B_array = np.append(temp_sdss_M_B_array, temp_SDSS_M_B[index])
						temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, temp_SDSS_sfr[index])
				if len(temp_sdss_metallicity_array) > binning_cut:
					temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
					temp_median_x_axis = np.append(temp_median_x_axis, bin_limit+(M_B_bins_lower_limit[1]-M_B_bins_lower_limit[0])/2)
					temp_median_z_axis = np.append(temp_median_z_axis, np.mean(temp_sdss_sfr_array))
			color_value = (sfr_bin_limit-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
			
			if (color_value <= 1.0):
				print('color_value: '+str(color_value))
				ax.plot(temp_median_x_axis-weight*(temp_median_z_axis), temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
	#sfr_bins_lower_limit = np.linspace(-3.5, 1.0, 6) #np.array([-1000.0, -5.00, -4.00, -3.00, -2.00, -1.00, 0.00, 1.00])
	#sfr_bins_lower_limit = np.array([-3.5, -1.7, -0.8, 0.1, 1.0])		
	temp_metallicity_array = np.array([])
	temp_mass_array = np.array([])
	temp_M_B_array = np.array([])
	temp_metallicity_array_error = np.array([])
	temp_sfr_array = np.array([])
	for index, galaxy_sfr in enumerate(Jimmy_sfr):
		if (galaxy_sfr < sfr_bins_lower_limit[2]) and (galaxy_sfr > sfr_bins_lower_limit[1]):
			temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
			temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
			temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
			temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
			temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
	if (plot_to_make == 'fmr_sfr'):
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array-weight*(temp_sfr_array), temp_metallicity_array, sigma=temp_metallicity_array_error)
		x_fit = np.linspace(min(temp_mass_array-weight*(temp_sfr_array)), max(temp_mass_array-weight*(temp_sfr_array)), 100)
	if (plot_to_make == 'fmr_lzr_sfr'):
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array-weight*(temp_sfr_array), temp_metallicity_array, sigma=temp_metallicity_array_error)
		x_fit = np.linspace(min(temp_M_B_array-weight*(temp_sfr_array)), max(temp_M_B_array-weight*(temp_sfr_array)), 100)
	#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
	y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
	if (plot_to_make == 'fmr_sfr'):
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_mass_array-weight*(temp_sfr_array))))
	if (plot_to_make == 'fmr_lzr_sfr'):
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_M_B_array-weight*(temp_sfr_array))))
	var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
	sd_res = np.sqrt(var_res)
	print('standard deviation: '+str(sd_res))
	#y_fit = linear_fit(x_fit, 4, 0.5)
	color_value = (sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
	#print('color_value: '+str(cmap(color_value)))
	##(a,b,(ci_a, ci_b), ri) = fitLine(temp_mass_array-weight*(temp_sfr_array-9.80),temp_metallicity_array)
	if (combine_lowest_bins != 'y'):
		ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
	y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
	y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
	if (combine_lowest_bins != 'y'):
		ax.fill_between(x_fit, y_fit_upper, y_fit_lower, facecolor=cmap(color_value), alpha=0.25)
	#ax.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
	if (plot_to_make == 'fmr_sfr'):
		if (combine_lowest_bins != 'y'):
			ax.errorbar(temp_mass_array-weight*(temp_sfr_array),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
	if (plot_to_make == 'fmr_lzr_sfr'):
		ax.errorbar(temp_M_B_array-weight*(temp_sfr_array),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
	
		
	temp_metallicity_array = np.array([])
	temp_mass_array = np.array([])
	temp_M_B_array = np.array([])
	temp_metallicity_array_error = np.array([])
	temp_sfr_array = np.array([])
	for index, galaxy_sfr in enumerate(Jimmy_sfr):
		if (galaxy_sfr < sfr_bins_lower_limit[1]) and (galaxy_sfr > sfr_bins_lower_limit[0]):
			temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
			temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
			temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
			temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])				
			temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
	#for index, galaxy_sfr in enumerate(Saintonge_sfr):
	#	if (Saintonge_sfr[index] < sfr_bins_lower_limit[2]) and (Saintonge_sfr[index] > sfr_bins_lower_limit[1]): #and (Saintonge_metallicity[index] > 8.0)
	#		temp_sfr_array = np.append(temp_sfr_array, Saintonge_sfr[index])	
	#		temp_mass_array = np.append(temp_mass_array, Saintonge_mass[index]	)
	#		temp_metallicity_array = np.append(temp_metallicity_array, Saintonge_metallicity[index])
	#for index, galaxy_sfr in enumerate(SDSS_sfr):
	#	if (galaxy_sfr < -2.00) and (galaxy_sfr > -2.75):
	#		temp_sfr_array = np.append(temp_sfr_array, SDSS_sfr[index])
	#		temp_mass_array = np.append(temp_mass_array, SDSS_mass[index])
	#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
	#print(temp_metallicity_array)
	#moustakas_guess = np.array([8.901, 8.798, 0.640])
	#linear_guess = np.array([6.24683786, 0.26661621])
	#fit_results, fit_error = optimization.curve_fit(moustakas, temp_mass_array, temp_metallicity_array, p0 = moustakas_guess)
	if (plot_to_make == 'fmr_sfr'):
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array-weight*(temp_sfr_array), temp_metallicity_array, sigma=temp_metallicity_array_error)
		x_fit = np.linspace(min(temp_mass_array-weight*(temp_sfr_array)), max(temp_mass_array-weight*(temp_sfr_array)), 10)
	if (plot_to_make == 'fmr_lzr_sfr'):
		fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array-weight*(temp_sfr_array), temp_metallicity_array, sigma=temp_metallicity_array_error)
		x_fit = np.linspace(min(temp_M_B_array-weight*(temp_sfr_array)), max(temp_M_B_array-weight*(temp_sfr_array)), 10)
	#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
	y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
	#print('x_fit: '+str(x_fit))
	#print('y_fit: '+str(y_fit))
	if (plot_to_make == 'fmr_sfr'):
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_mass_array-weight*(temp_sfr_array))))
	if (plot_to_make == 'fmr_lzr_sfr'):
		residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_M_B_array-weight*(temp_sfr_array))))
	var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
	sd_res = np.sqrt(var_res)
	print('standard deviation: '+str(sd_res))
	color_value = (sfr_bins_lower_limit[0]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
	##(a,b,(ci_a, ci_b), ri) = fitLine(temp_mass_array-weight*(temp_sfr_array-9.80),temp_metallicity_array)
	if (combine_lowest_bins != 'y'):
		ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
	y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
	y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
	if (combine_lowest_bins != 'y'):
		ax.fill_between(x_fit, y_fit_upper, y_fit_lower, facecolor=cmap(color_value), alpha=0.25)
	#ax.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
	if (plot_to_make == 'fmr_sfr'):
		if (combine_lowest_bins != 'y'):
			ax.errorbar(temp_mass_array-weight*(temp_sfr_array),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
	if (plot_to_make == 'fmr_lzr_sfr'):
		ax.errorbar(temp_M_B_array-weight*(temp_sfr_array),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
	
	if (combine_lowest_bins == 'y'):
		temp_metallicity_array = np.array([])
		temp_mass_array = np.array([])
		temp_M_B_array = np.array([])
		temp_metallicity_array_error = np.array([])
		temp_sfr_array = np.array([])
		for index, galaxy_sfr in enumerate(Jimmy_sfr):
			if (galaxy_sfr < sfr_bins_lower_limit[2]) and (galaxy_sfr > sfr_bins_lower_limit[0]):
				temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
				temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
				temp_M_B_array = np.append(temp_M_B_array, Jimmy_M_B[index])
				temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
				temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
		if (plot_to_make == 'fmr_sfr'):
			fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array-weight*(temp_sfr_array), temp_metallicity_array, sigma=temp_metallicity_array_error)
			x_fit = np.linspace(min(temp_mass_array-weight*(temp_sfr_array)), max(temp_mass_array-weight*(temp_sfr_array)), 100)
		if (plot_to_make == 'fmr_lzr_sfr'):
			fit_results, fit_error = optimization.curve_fit(linear_fit, temp_M_B_array-weight*(temp_sfr_array), temp_metallicity_array, sigma=temp_metallicity_array_error)
			x_fit = np.linspace(min(temp_M_B_array-weight*(temp_sfr_array)), max(temp_M_B_array-weight*(temp_sfr_array)), 100)
		#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
		y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
		if (plot_to_make == 'fmr_sfr'):
			residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_mass_array-weight*(temp_sfr_array))))
		if (plot_to_make == 'fmr_lzr_sfr'):
			residuals = temp_metallicity_array - (fit_results[0]+(fit_results[1]*(temp_M_B_array-weight*(temp_sfr_array))))
		var_res = np.sum(residuals**2)/(len(temp_metallicity_array)-2)
		sd_res = np.sqrt(var_res)
		print('standard deviation: '+str(sd_res))
		#y_fit = linear_fit(x_fit, 4, 0.5)
		color_value = (sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])
		#print('color_value: '+str(cmap(color_value)))
		##(a,b,(ci_a, ci_b), ri) = fitLine(temp_mass_array-weight*(temp_sfr_array-9.80),temp_metallicity_array)
		ax.plot(x_fit, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0)
		y_fit_upper = linear_fit(x_fit, fit_results[0]+sd_res, fit_results[1])
		y_fit_lower = linear_fit(x_fit, fit_results[0]-sd_res, fit_results[1])
		ax.fill_between(x_fit, y_fit_upper, y_fit_lower, facecolor=cmap(color_value), alpha=0.25)
		#ax.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
		if (plot_to_make == 'fmr_sfr'):
			ax.errorbar(temp_mass_array-weight*(temp_sfr_array),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
		if (plot_to_make == 'fmr_lzr_sfr'):
			ax.errorbar(temp_M_B_array-weight*(temp_sfr_array),temp_metallicity_array,yerr=temp_metallicity_array_error, linestyle="None", ecolor=cmap(color_value))
	
	
	#x_fit = np.linspace(5.50, 9.50, 100)
	#y_fit = linear_fit(x_fit,4.5,0.46)
	#color_value = 0.25
	#ax.plot(x_fit-weight*-3.5, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0, alpha=0.5)
	#x_fit = np.linspace(6.50, 10.25, 100)
	#y_fit = linear_fit(x_fit, 4.92,0.46)
	#color_value = 0.0
	#ax.plot(x_fit-weight*-2.5, y_fit, color=cmap(color_value), linestyle='--', linewidth=3.0, alpha=0.5)
	
	
	#print(Jimmy_sfr)
	ax.plot([0,0], [0,0], color='black', label=r'M$_*$ Bin Mean')
	if (plot_to_make == 'fmr_sfr'):
		s = ax.scatter(Jimmy_mass-weight*(Jimmy_sfr), Jimmy_metallicity, c=Jimmy_sfr, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
		ax.text(solar_x_axis+0.5, solar_line, r'$\odot$', fontsize=14, fontweight='bold')		
		plt.xlabel(r'log(M$_{*}$)-'+str(weight)+'log(SFR)',fontsize=fontsize, labelpad=20)
		plt.xlim(6.5, 11.5)
	if (plot_to_make == 'fmr_lzr_sfr'):
		s = ax.scatter(Jimmy_M_B-weight*(Jimmy_sfr), Jimmy_metallicity, c=Jimmy_sfr, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
		ax.text(solar_x_axis+0.5, solar_line, r'$\odot$', fontsize=14, fontweight='bold')		
		plt.xlabel(r'M$_{B}$+'+str(abs(weight))+'log(SFR)',fontsize=fontsize, labelpad=20)
		plt.xlim(-12, -22)
	#logmstar = 7.0+(np.arange(400)*0.01)
	#PP04_line = 32.1488-(8.51258*logmstar)+(0.976384*logmstar*logmstar)-(0.0359763*logmstar*logmstar*logmstar)
	#ax.plot(logmstar,PP04_line, color="black", linewidth=3.0, label="PP04 SDSS Relation")
	plt.ylabel('12+log(O/H)',fontsize=fontsize, labelpad=10)
	plt.ylim(7.6, 9.3)	
	plt.minorticks_on()
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	#plt.legend(loc='upper left')
	artist = plt.Rectangle((0, 0), 1, 1, fc="gray", ec=None, alpha = sdss_alpha_value)
	handles, labels = ax.get_legend_handles_labels()
	handles = handles+[artist]
	labels = labels+[u'ALFALFA/SDSS']
	ax.legend(handles, labels, loc='upper left', ncol=2)
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin-0.02)
	if (dog == 'y'):
		plt.subplots_adjust(right=right_margin-0.04)
	plt.subplots_adjust(top=top_margin)
	#ax2 = fig.add_axes([0.855, 0.175, 0.03, 0.775])
	#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
	#cb.ax.tick_params(labelsize=fontsize-6)
	#cb.set_label(r'SFR M$_\odot$ yr$^{-1}$', fontsize=fontsize, rotation=270)
	cb = plt.colorbar(s)
	#plt.ticks.set_fontsize(20)
	cb.ax.tick_params(labelsize=fontsize-6)
	if dog != 'y':
		cb.set_label(r'SFR (M$_\odot$ yr$^{-1})$', fontsize=fontsize, rotation=270, labelpad=labelpad)
	if dog =='y':
		cb.set_label(r'SFR (M$_\odot$ dyr$^{-1})$', fontsize=fontsize, rotation=270, labelpad=labelpad)
	if (calibration == 'N2'):
		if (combine_lowest_bins != 'y'):
			plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
		if (combine_lowest_bins == 'y'):
			plt.savefig(HOME+PAPER_FOLDER+'combined_bins_fmr_sfr.pdf')
	if (calibration == 'O3N2'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
	#plt.show()

if (plot_to_make == 'fake_h_alpha_map'):
	galaxy = 'AGC666'
	#fig = plt.figure()
	voronoi_2d_bins_file = HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/voronoi_2d_bins_emission.txt'
	voronoi_2d_bins_table=np.loadtxt(voronoi_2d_bins_file, dtype=float)
	xbin = voronoi_2d_bins_table[:,0]
	ybin = voronoi_2d_bins_table[:,1]
	
	gandalf_file = HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/pink_gandalf.fits'
	#Open the fits file that contains the emission line fitting results
	hdulist = pyfits.open(gandalf_file)
	#In gandalf the emission line fit results are in the 4th extension	
	kinematic_measurements = hdulist[4].data
	velocity = kinematic_measurements[:,43]
	flux_OIII2 = kinematic_measurements[:,31]
	flux_Hb = kinematic_measurements[:,21]
	flux_NII2 = kinematic_measurements[:,46]
	flux_Ha = kinematic_measurements[:,41]
	
	one_bin_table = np.genfromtxt(HOME+'/Astro/reduced/'+galaxy+'pro/all/'+sncut+'/gandalf_table.txt',dtype=None)	
	Hb_line = one_bin_table[4]
	stacked_Hb = Hb_line[1]
	OIII2_line = one_bin_table[3]
	stacked_OIII2 = OIII2_line[1]
	Ha_line = one_bin_table[6]
	stacked_Ha = Ha_line[1]
	NII2_line = one_bin_table[5]
	stacked_NII2 = NII2_line[1]
	Jimmy_N2 = np.log10(stacked_NII2/stacked_Ha)
	Jimmy_O3N2 = np.log10((stacked_OIII2/stacked_Hb)/(stacked_NII2/stacked_Ha))
	stacked_PP04_O3N2 = 8.73-(0.32*Jimmy_O3N2)
	stacked_PP04_N2 = 9.37 + (2.03*Jimmy_N2) + (1.26*Jimmy_N2**2) + (0.32*Jimmy_N2**3)
	stacked_PP04_N2_linear = 8.90 + 0.57*Jimmy_N2
	stacked_D02_N2 = 9.12 + 0.73*Jimmy_N2
	if (calibration == 'N2'):
		Jimmy_metallicity = stacked_D02_N2
	if (calibration == 'O3N2'):	
		Jimmy_metallicity = stacked_PP04_O3N2
	
	O3N2 = np.log10((flux_OIII2/flux_Hb)/(flux_NII2/flux_Ha))
	#print('Individual O3N2 values: '+str(O3N2))
	logOH_PP04_O3N2 = 8.73-(0.32*O3N2)
	N2 = np.log10(flux_NII2/flux_Ha)
	logOH_PP04_N2 = 9.37 + (2.03*N2) + (1.26*N2**2) + (0.32*N2**3) # FROM 
	logOH_PP04_N2_linear = 8.90 + 0.57*N2
	logOH_D02_N2 = 9.12 + 0.73*N2
	#logOH_PP04_O3N2 = 9.12+(0.73*N2)
	colorbar_min = 0.00
	colorbar_max = 5.00
	#colormap = 'winter' #This is blue to green
	colormap ='winter_r'
	#xbin - xbin+20
	#if (calibration == 'N2'):
	#	fig2 = display_pixels(xbin+10, ybin+10, logOH_D02_N2, colormap)
	#if (calibration == 'O3N2'):
	#	fig2 = display_pixels(xbin, ybin, logOH_PP04_O3N2, colormap)
	fig2 = display_pixels(xbin+10, ybin+13, flux_Ha, colormap)
	line_of_interest = one_bin_table[1] #0 is for PP04_O3N2, 1 is for PP04_N2, and 2 is for D02_N2
	#stacked_PP04_O3N2 = line_of_interest[1]
	xmin = -12
	xmax = 11
	ymin = -12
	ymax = 11
	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)
	#plt.clim(colorbar_min, colorbar_max)
	#plt.text(-10, -10,'12+log(O/H) = '+str(round(Jimmy_metallicity,2)), fontsize=fontsize)
	plt.text(-11, -8,r'Integrated H$\alpha$ = '+str(round(stacked_Ha,2))+' $10^{-16}$ erg ', fontsize=fontsize)
	plt.text(5, -10,r'cm$^{-2}$ s$^{-1}$', fontsize=fontsize)
	plt.xlabel('arcsec',fontsize=fontsize)
	plt.ylabel('arcsec',fontsize=fontsize)
	
	#ax.set_xticklabels([-5,0,5,10,15])
	plt.tick_params(axis='y', labelsize=20)
	plt.tick_params(axis='x', labelsize=20)
	#cbaxes = fig.add_axes([0.90, 0.35, 0.015, 0.6])
	#plt.colorbar(fig2, cax=cbaxes)
	cb = plt.colorbar(fig2)
	cb.ax.tick_params(labelsize=fontsize-6)
	cb.set_label(r'H$\alpha$ flux [$10^{-16}$ erg cm$^{-2}$ s$^{-1}$]', fontsize=fontsize, rotation=270, labelpad=labelpad)
	
	plt.subplots_adjust(bottom=0.14)
	plt.subplots_adjust(left=0.12)
	plt.subplots_adjust(right=0.92)
	plt.subplots_adjust(top=0.98)
	#plt.show()
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	
if (plot_to_make == 'colorful_gas_fraction_sfr'):
	fig=plt.figure()
	# define the colormap
	cmap = plt.cm.jet
	#cmap = plt.cm.jet_r
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# force the first color entry to be grey
	#cmaplist[0] = (.5,.5,.5,1.0)
	#cmaplist[0] = (0.0,0.0,0.0,1.0)
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	# define the bins and normalize
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	plt.fill_between(range(20), 9.1, 12.0, facecolor='yellow', alpha=0.25)
	plt.fill_between(range(20), 7.0, 7.6, facecolor='yellow', alpha=0.25)
	plt.scatter(SDSS_mass, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS', rasterized=True)
	
	
	#create bins and find means within them.
	for sfr_bin_index, sfr_bin_limit in enumerate(sfr_bins_lower_limit):
		#print(sfr_bin_limit)
		temp_sdss_metallicity_array = np.array([])
		temp_sdss_mass_array = np.array([])
		temp_sdss_hi_mass_array = np.array([])
		temp_sdss_sfr_array = np.array([])
		for index, galaxy_sfr in enumerate(SDSS_sfr):
			if (galaxy_sfr > sfr_bin_limit) and (galaxy_sfr < (sfr_bin_limit+(sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0]))):
				temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,SDSS_metallicity[index])
				temp_sdss_mass_array = np.append(temp_sdss_mass_array, SDSS_mass[index])
				temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, SDSS_hi_mass[index])
				temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, galaxy_sfr)
		temp_SDSS_mass = temp_sdss_mass_array
		temp_SDSS_hi_mass = temp_sdss_hi_mass_array
		temp_SDSS_metallicity = temp_sdss_metallicity_array
		temp_SDSS_sfr = temp_sdss_sfr_array
		
		#moustakas_guess = np.array([8.901, 8.798, 0.640])
		#fit_results, fit_error = optimization.curve_fit(moustakas, temp_sdss_mass_array, temp_sdss_metallicity_array, p0 = moustakas_guess)
		#moustakas_model = Model(odr_moustakas)
		#data = RealData(temp_sdss_mass_array, temp_sdss_metallicity_array)
		#odr = ODR(data, moustakas_model, beta0=moustakas_guess)
		#out = odr.run()
		#out.pprint()
		
		temp_median_x_axis = np.array([])
		temp_median_y_axis = np.array([])
		for bin_limit in mass_bins_lower_limit:
			#print('bin_limit'+str(bin_limit))
			temp_sdss_metallicity_array = np.array([])
			temp_sdss_mass_array = np.array([])
			temp_sdss_hi_mass_array = np.array([])
			temp_sdss_sfr_array = np.array([])
			for index, galaxy_mass in enumerate(temp_SDSS_mass):
				if (galaxy_mass > bin_limit) and (galaxy_mass < (bin_limit+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0]))):
					temp_sdss_metallicity_array = np.append(temp_sdss_metallicity_array,temp_SDSS_metallicity[index])
					temp_sdss_mass_array = np.append(temp_sdss_mass_array, temp_SDSS_mass[index])
					temp_sdss_hi_mass_array = np.append(temp_sdss_hi_mass_array, temp_SDSS_hi_mass[index])
					temp_sdss_sfr_array = np.append(temp_sdss_sfr_array, temp_SDSS_sfr[index])
			if len(temp_sdss_metallicity_array) > binning_cut:
				#print(temp_sdss_metallicity_array)
				temp_median_y_axis = np.append(temp_median_y_axis, np.mean(temp_sdss_metallicity_array))
				temp_median_x_axis = np.append(temp_median_x_axis, (bin_limit-temp_SDSS_hi_mass[index])+(mass_bins_lower_limit[1]-mass_bins_lower_limit[0])/2)
		#print('temp x-axis line: '+str(temp_median_x_axis))
		#print('temp y-axis line: '+str(temp_median_y_axis))
		color_value = (sfr_bin_limit-sfr_bins_lower_limit[0])/(sfr_bins_lower_limit[-2]-sfr_bins_lower_limit[0])+0.10
		#print(color_value)
		#print(len(temp_median_y_axis))
		plt.plot(temp_median_x_axis, temp_median_y_axis, color=cmap(color_value), linewidth=3.0)
		#x_fit = np.linspace(6, 11, 1000)
		#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
		#y_fit = odr_moustakas(out.beta, x_fit)
		#plt.plot(x_fit, y_fit, color=cmap(color_value))
	
	
	temp_metallicity_array = np.array([])
	temp_mass_array = np.array([])
	temp_hi_mass_array = np.array([])
	temp_metallicity_array_error = np.array([])
	temp_sfr_array = np.array([])
	#print(Jimmy_sfr)
	for index, galaxy_sfr in enumerate(Jimmy_sfr):
		if (galaxy_sfr < sfr_bins_lower_limit[3]) and (galaxy_sfr > sfr_bins_lower_limit[2]):
			temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
			temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
			temp_hi_mass_array = np.append(temp_hi_mass_array, Jimmy_hi_mass[index])
			temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
			temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
	#for index, galaxy_sfr in enumerate(SHIELD_sfr):
	#	if (SHIELD_metallicity[index] < 11) and (SHIELD_mass[index] > 1.00) and (SHIELD_sfr[index] < sfr_bins_lower_limit[2]) and (SHIELD_sfr[index] > sfr_bins_lower_limit[1]): #and (SHIELD_metallicity[index] > 8.0)
	#		temp_sfr_array = np.append(temp_sfr_array, SHIELD_sfr[index])
	#		temp_mass_array = np.append(temp_mass_array, SHIELD_mass[index])
	#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
	#for index, galaxy_sfr in enumerate(SDSS_sfr):
	#	if (galaxy_sfr < -2.00) and (galaxy_sfr > -2.75):
	#		temp_sfr_array = np.append(temp_sfr_array, SDSS_sfr[index])
	#		temp_mass_array = np.append(temp_mass_array, SDSS_mass[index])
	#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
	#print(temp_metallicity_array)
	#moustakas_guess = np.array([8.901, 8.798, 0.640])
	#linear_guess = np.array([6.24683786, 0.26661621])
	#fit_results, fit_error = optimization.curve_fit(moustakas, temp_mass_array, temp_metallicity_array, p0 = moustakas_guess)
	#fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array, temp_metallicity_array, p0 = linear_guess)
	#print(fit_results)
	#perr = np.sqrt(np.diag(fit_error))
	#print(fit_error[1][1])
	#x_fit = np.linspace(5.5, 9.0, 100)
	#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
	#y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
	#plt.plot(x_fit, y_fit, color=cmap(0.20), linestyle='--', linewidth=3.0)
	#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, facecolor=cmap(0.20), alpha=0.25)
	#plt.fill_between(x_fit, 7.0, 7.6, facecolor='yellow', alpha=0.25)
	
	
	temp_metallicity_array = np.array([])
	temp_mass_array = np.array([])
	temp_hi_mass_array = np.array([])
	temp_metallicity_array_error = np.array([])
	temp_sfr_array = np.array([])
	for index, galaxy_sfr in enumerate(Jimmy_sfr):
		if (galaxy_sfr < sfr_bins_lower_limit[4]) and (galaxy_sfr > sfr_bins_lower_limit[3]):
			temp_sfr_array = np.append(temp_sfr_array, Jimmy_sfr[index])
			temp_mass_array = np.append(temp_mass_array, Jimmy_mass[index])
			temp_hi_mass_array = np.append(temp_hi_mass_array, Jimmy_hi_mass[index])
			temp_metallicity_array_error = np.append(temp_metallicity_array_error, Jimmy_metallicity_error[index])
			temp_metallicity_array = np.append(temp_metallicity_array, Jimmy_metallicity[index])
	#for index, galaxy_sfr in enumerate(SHIELD_sfr):
	#	if (SHIELD_metallicity[index] < 11) and (SHIELD_mass[index] > 1.00) and (SHIELD_metallicity[index] > 8.0):
	#		temp_sfr_array = np.append(temp_sfr_array, SHIELD_sfr[index])
	#		temp_mass_array = np.append(temp_mass_array, SHIELD_mass[index])
	#		temp_metallicity_array = np.append(temp_metallicity_array, SHIELD_metallicity[index])
	#for index, galaxy_sfr in enumerate(SDSS_sfr):
	#	if (galaxy_sfr < sfr_bins_lower_limit[4]) and (galaxy_sfr > sfr_bins_lower_limit[3]):
	#		temp_sfr_array = np.append(temp_sfr_array, SDSS_sfr[index])
	#		temp_mass_array = np.append(temp_mass_array, SDSS_mass[index])
	#		temp_metallicity_array = np.append(temp_metallicity_array, SDSS_metallicity[index])
	#moustakas_guess = np.array([8.901, 8.798, 0.640])
	#print(temp_metallicity_array)
	#fit_results, fit_error = optimization.curve_fit(moustakas, temp_mass_array, temp_metallicity_array, p0 = moustakas_guess)
	#fit_results, fit_error = optimization.curve_fit(linear_fit, temp_mass_array, temp_metallicity_array)
	#x_fit = np.linspace(7.0, 9.5, 100)
	#y_fit = moustakas(x_fit, fit_results[0], fit_results[1], fit_results[2])
	#y_fit = linear_fit(x_fit, fit_results[0], fit_results[1])
	#plt.fill_between(x_fit, y_fit_upper, y_fit_lower, facecolor=cmap(0.35), alpha=0.25)
	#plt.plot(x_fit, y_fit, color=cmap(0.35), linestyle='--', linewidth=3.0)
	plt.plot([0,0],[0,0], color='black', linestyle='--', label='Jimmy+15 Linear Fit')
	
	
	plt.plot([0,0],[0,0], color='black', label='SDSS Mean')
	#plt.scatter(SDSS_mass, SDSS_metallicity, c=SDSS_sfr-SDSS_mass, alpha=0.5, cmap=cmap, norm=norm, label="ALFALFA", rasterized=True)
	#print(Jimmy_sfr)
	#print(len(Jimmy_sfr))
	#print(len(Jimmy_mass))
	plt.scatter(Jimmy_mass-Jimmy_hi_mass, Jimmy_metallicity, c=Jimmy_sfr, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15", color='black')
	#plt.scatter(James_mass, James_metallicity, c=James_sfr, s=40, cmap=cmap, norm=norm, marker='s', label="James+14")
	plt.scatter(SHIELD_mass-SHIELD_hi_mass, SHIELD_metallicity, c=SHIELD_sfr, s=40, cmap=cmap, norm=norm, label="SHIELD", color='black')
	#plt.scatter(Jimmy_mass, Jimmy_metallicity, c=Jimmy_sfr-Jimmy_mass, s=400, cmap=cmap, norm=norm, marker='*', label="Jimmy+15")
	#plt.axhline(y=solar_line, color="red", linewidth=2.0, linestyle='--', label="Solar")
	plt.text(solar_x_axis, solar_line, r'$\odot$', fontsize=14, fontweight='bold')
	logmstar = 7.0+(np.arange(400)*0.01)
	PP04_line = 32.1488-(8.51258*logmstar)+(0.976384*logmstar*logmstar)-(0.0359763*logmstar*logmstar*logmstar)
	#plt.plot(logmstar,PP04_line, color="black", linewidth=3.0, label="PP04 SDSS Relation")
	plt.ylabel('12+log(O/H)',fontsize=fontsize, labelpad=10)
	plt.xlabel(r'log(M$_{*}$) [log(M$_{\odot}$)]',fontsize=fontsize, labelpad=20)
	plt.xlim(-3.0, 3.0)
	plt.ylim(7.6, 9.3)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	plt.minorticks_on()
	plt.legend(loc='upper left')
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	#cbaxes = fig.add_axes([0.875, 0.175, 0.03, 0.775])
	#cbar = plt.colorbar(cax=cbaxes)
	#cbar.set_label(r'SFR (M$_\odot$ yr$^{-1}$', fontsize=fontsize, rotation=270)
	ax2 = fig.add_axes([0.855, 0.175, 0.03, 0.775])
	#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
	#cb.ax.tick_params(labelsize=fontsize-6)
	#cb.set_label(r'SFR (M$_\odot$ yr$^{-1})$', fontsize=fontsize, rotation=270, labelpad=labelpad)
	cb = plt.colorbar()
	#plt.ticks.set_fontsize(20)
	cb.ax.tick_params(labelsize=fontsize-6)
	cb.set_label(r'SFR (M$_\odot$ yr$^{-1})$', fontsize=fontsize, rotation=270, labelpad=labelpad)
	#if (calibration == 'N2'):
	#	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#if (calibration == 'O3N2'):
	#	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'_O3N2.pdf')
	plt.show()

def gauss(x, *p):
	A, mu, sigma = p
	return A*np.exp(-(x-mu)**2/(2.*sigma**2))
def func(x, a, b):
	return a*x + b
if (plot_to_make == 'D02_line_fit_output'):
	galaxy = 'AGC221000'
	sncut = 'sn3'
	
	data_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+galaxy+"all.fits"
	cube_hdu = pyfits.open(data_file)
	img = cube_hdu[0].data
	img_header = cube_hdu[0].header
	cdelt = float(img_header["CDELT3"])
	crval = float(img_header["CRVAL3"])
	crpix = float(img_header["CRPIX3"])
	wavelength_angstroms = crval + (np.arange(img.shape[0])*cdelt)
	
	bins_file = HOME+"/Astro/reduced/"+galaxy+"pro/all/"+sncut+"/voronoi_2d_binning_output_emission.txt"
	bins_table = np.loadtxt(bins_file, dtype=float)
	bins = bins_table[:,4]
	x_pix = bins_table[:,2]
	y_pix = bins_table[:,3]
	chosen_bin = 28 #8
	
	img_spectrum = img[:,0,0]*0
	for index, bin_num in enumerate(bins):
		if bin_num == chosen_bin:
			img_spectrum = img_spectrum+img[:,y_pix[index],x_pix[index]]
	linear_fit_result = func(wavelength_angstroms, -6.97825685e-05, 4.88871756e-01)
	cont_subtracted_img = img_spectrum
	cont_subtracted_img[:350] = img_spectrum[:350]-linear_fit_result[:350]
	linear_fit_result = func(wavelength_angstroms, -3.28811279e-05, 2.92187196e-01)
	cont_subtracted_img[350:] = img_spectrum[350:]-linear_fit_result[350:]

	fig = plt.figure(figsize=(10,8))
	ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
	plt.plot(wavelength_angstroms, cont_subtracted_img, color="black", linewidth=6, label="Continuum Subtracted Spectrum", zorder=1)
	#gaussian_fit_results = gauss(wavelength_angstroms, 3.92202216e-01, 4.88613202e+03, 6.58036956e+00)
	#plt.plot(wavelength_angstroms, gaussian_fit_results, color="orange", linewidth=3, label=r"H$\beta$ Fit")
	#gaussian_fit_results = gauss(wavelength_angstroms, 1.92224458e-01, 4.98441310e+03, 6.77983032e+00)
	#plt.plot(wavelength_angstroms, gaussian_fit_results, color="purple", linewidth=3, label=r"[OIII] 4959 $\AA$ Fit")
	#gaussian_fit_results = gauss(wavelength_angstroms, 6.41776704e-01, 5.03245968e+03, 7.03859758e+00)
	#plt.plot(wavelength_angstroms, gaussian_fit_results, color="green", linewidth=3, label=r"[OIII] 5007 $\AA$ Fit")
	gaussian_fit_results = gauss(wavelength_angstroms, 1.04579050e+00, 6.59660324e+03, 7.46077479e+00)
	plt.plot(wavelength_angstroms, gaussian_fit_results, color="green", linewidth=2, label=r"H$\alpha$ Fit", linestyle='-', zorder=2)
	plt.fill_between(wavelength_angstroms, gaussian_fit_results, 0, color='green', alpha=0.5, zorder=3)
	residual = cont_subtracted_img-gaussian_fit_results
	gaussian_fit_results = gauss(wavelength_angstroms, 9.78426689e-02, 6.61964788e+03, 6.48601406e+00)
	plt.plot(wavelength_angstroms, gaussian_fit_results, color="red", linewidth=2, label=r"[NII] 6583 \AA\ Fit", linestyle='-', zorder=4)
	plt.fill_between(wavelength_angstroms, gaussian_fit_results, 0, color='red', alpha=0.5, zorder=5)
	residual = residual-gaussian_fit_results
	gaussian_fit_results = gauss(wavelength_angstroms, 8.76925018e-03, 6.58456257e+03, 7.46077479e+00)
	plt.plot(wavelength_angstroms, gaussian_fit_results, color="blue", linewidth=2, label=r"[NII] 6549 \AA\ Fit", linestyle='-', zorder=6)
	plt.fill_between(wavelength_angstroms, gaussian_fit_results, 0, color='blue', alpha=0.5, zorder=7)
	residual = residual-gaussian_fit_results
	#plt.plot(wavelength_angstroms, wavelength_angstroms*0.0, linewidth=3, color="blue")
	
	#ax1.set_yscale('log')
	plt.xlim(6420,6700)
	plt.ylim(-0.2, 1.05)
	plt.xticks([])
	#plt.xlabel(r'Wavelength (\AA)',fontsize=fontsize, labelpad=20)
	
	plt.legend(loc="upper left", ncol=1, prop={'size':17})
	plt.tick_params(axis='both', which='major', labelsize=20)
	plt.subplots_adjust(bottom=0.13)
	plt.subplots_adjust(left=0.13)
	plt.subplots_adjust(right=0.98)
	plt.subplots_adjust(top=0.98)
	#plt.subplots_adjust(wspace=0.23)
	#plt.subplots_adjust(hspace=0.17)
	
	
	ax2 = plt.subplot2grid((4,1), (3, 0))

	plt.plot(wavelength_angstroms, residual, linewidth=2, label="Residual", color='purple')

	plt.xlim(6420,6700)
	plt.ylim(-0.2, 0.2)
	plt.yticks([-0.2, -0.1, 0.0, 0.1, 0.2])
	plt.xlabel(r'Wavelength (\AA)',fontsize=fontsize, labelpad=20)
	#plt.ylabel(r'Flux $10^{-16}$ erg cm$^{-2}$ s$^{-1}$',fontsize=fontsize, labelpad=10)
	plt.legend(loc="upper left", ncol=1, prop={'size':17})
	plt.tick_params(axis='both', which='major', labelsize=20)
	plt.subplots_adjust(bottom=0.13)
	plt.subplots_adjust(left=0.13)
	plt.subplots_adjust(right=0.98)
	plt.subplots_adjust(top=0.98)
	#plt.subplots_adjust(wspace=0.23)
	#plt.subplots_adjust(hspace=0.17)
	fig.text(0.04, 0.5, r'Flux $10^{-16}$ erg cm$^{-2}$ s$^{-1}$',fontsize=fontsize, ha='center', va='center', rotation='vertical')
	#plt.show()
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	

if (plot_to_make == 'consistency_check'):
	#For galaxies with ALFALFA/SDSS variables, cross check between them and the IFU derived variables.
	Comparison_SDSS_names = np.array([])
	Comparison_SDSS_mass = np.array([])
	Comparison_mpa_jhu_mass = np.array([])
	Comparison_SDSS_mass_error = np.array([])
	Comaprison_SDSS_hi_mass = np.array([])
	Comparison_SDSS_M_B = np.array([])
	Comparison_SDSS_sfr = np.array([])
	Comparison_SDSS_metallicity = np.array([])
	Comparison_Jimmy_names = np.array([])
	Comparison_Jimmy_mass = np.array([])
	Comparison_Jimmy_mass_error = np.array([])
	Comaprison_Jimmy_hi_mass = np.array([])
	Comparison_Jimmy_M_B = np.array([])
	Comparison_Jimmy_M_B_error = np.array([])
	Comparison_Jimmy_sfr = np.array([])
	Comparison_Jimmy_sfr_error = np.array([])
	Comparison_Jimmy_metallicity = np.array([])
	Comparison_Jimmy_metallicity_error = np.array([])
	print(len(SDSS_names))
	print(len(SDSS_mass))
	for Jimmy_index, galaxy in enumerate(names):
		AGCN = galaxy[3:]
		for index, SDSS_AGCN in enumerate(SDSS_names):
			if str(SDSS_AGCN) == AGCN:
				print(SDSS_AGCN)
				#print(SDSS_mass[index])
				#print(SDSS_M_B[index])
				#print(SDSS_sfr[index])
				#print(SDSS_metallicity[index])
				Comparison_SDSS_names = np.append(Comparison_SDSS_names, 'AGC'+str(SDSS_names[index]))
				Comparison_SDSS_mass = np.append(Comparison_SDSS_mass, SDSS_mass[index])
				Comparison_mpa_jhu_mass = np.append(Comparison_mpa_jhu_mass, mass_missmatch_cut_mpa_jhu_mass[index])
				#Comparison_SDSS_mass_error = np.append(Comparison_SDSS_mass_error, SDSS_mass_error[index])
				Comparison_SDSS_M_B = np.append(Comparison_SDSS_M_B, SDSS_M_B[index])
				Comparison_SDSS_sfr = np.append(Comparison_SDSS_sfr, SDSS_sfr[index])
				Comparison_SDSS_metallicity = np.append(Comparison_SDSS_metallicity, SDSS_metallicity[index])
				Comaprison_SDSS_hi_mass = np.append(Comaprison_SDSS_hi_mass, SDSS_hi_mass[index])
				Comparison_Jimmy_names = np.append(Comparison_Jimmy_names, names[Jimmy_index])
				Comparison_Jimmy_mass = np.append(Comparison_Jimmy_mass, Jimmy_mass[Jimmy_index])
				Comparison_Jimmy_mass_error = np.append(Comparison_Jimmy_mass_error, Jimmy_mass_error[Jimmy_index])
				Comparison_Jimmy_M_B = np.append(Comparison_Jimmy_M_B, Jimmy_M_B[Jimmy_index])
				Comparison_Jimmy_M_B_error = np.append(Comparison_Jimmy_M_B_error, Jimmy_M_B_error[Jimmy_index])
				Comparison_Jimmy_sfr = np.append(Comparison_Jimmy_sfr, Jimmy_sfr[Jimmy_index])
				Comparison_Jimmy_sfr_error = np.append(Comparison_Jimmy_sfr_error, Jimmy_sfr_error[Jimmy_index])
				Comparison_Jimmy_metallicity = np.append(Comparison_Jimmy_metallicity, Jimmy_metallicity[Jimmy_index])
				Comparison_Jimmy_metallicity_error = np.append(Comparison_Jimmy_metallicity_error, Jimmy_metallicity_error[Jimmy_index])
				Comaprison_Jimmy_hi_mass = np.append(Comaprison_Jimmy_hi_mass, Jimmy_hi_mass[Jimmy_index])
	for Jimmy_index, galaxy in enumerate(names):
		AGCN = '221004'
		for index, SDSS_AGCN in enumerate(alfalfa_sdss_catalog['AGCNr']):
			if str(SDSS_AGCN) == AGCN:
				AGC221004_SDSS_mass = 6.01014444041
				AGC221004_SDSS_M_B = -10.5966351257
				AGC221004_SDSS_sfr = -3.60094857725
				AGC221004_SDSS_metallicity = 8.42292231033
				AGC221004_Jimmy_mass = Jimmy_mass[Jimmy_index]
				AGC221004_Jimmy_mass_error = Jimmy_mass_error[Jimmy_index]
				AGC221004_Jimmy_M_B = Jimmy_M_B[Jimmy_index]
				AGC221004_Jimmy_M_B_error = Jimmy_M_B_error[Jimmy_index]
				AGC221004_Jimmy_sfr = Jimmy_sfr[Jimmy_index]
				AGC221004_Jimmy_sfr_error = Jimmy_sfr_error[Jimmy_index]
				AGC221004_Jimmy_metallicity = Jimmy_metallicity[Jimmy_index]
				AGC221004_Jimmy_metallicity_error = Jimmy_metallicity_error[Jimmy_index]
	#print('AGC221004_Jimmy_mass: '+str(AGC221004_Jimmy_mass))
	#print(Comaprison_Jimmy_hi_mass)
	#print(Jimmy_hi_mass)
	#Stellar Mass		
	#Luminosity
	#Metallicity
	#SFR
	#Print a big table with SDSS derived, IFU derived, and %difference columns for each of the 4 above.

	for printindex, galaxy in enumerate(Comparison_Jimmy_names):
		#print(Comparison_SDSS_names[printindex]+' & '+str(round(Comparison_SDSS_mass[printindex],2))+' & '+str(round(Comparison_Jimmy_mass[printindex],2))+' & '+str(round(100*abs(10**Comparison_SDSS_mass[printindex]-10**Comparison_Jimmy_mass[printindex])/10**Comparison_Jimmy_mass[printindex],2))+' & '+str(round(Comparison_SDSS_M_B[printindex],2))+' & '+str(round(Comparison_Jimmy_M_B[printindex],2))+' & '+str(round(100*abs(2.5**Comparison_SDSS_M_B[printindex]-2.5*Comparison_Jimmy_M_B[printindex])/2.5**Comparison_Jimmy_M_B[printindex],2))+' & '+str(round(Comparison_SDSS_metallicity[printindex],2))+' & '+str(round(Comparison_Jimmy_metallicity[printindex],2))+' & '+str(round(100*abs(10**Comparison_SDSS_metallicity[printindex]-10**Comparison_Jimmy_metallicity[printindex])/10**Comparison_Jimmy_metallicity[printindex],2))+' & '+str(round(Comparison_SDSS_sfr[printindex],2))+' & '+str(round(Comparison_Jimmy_sfr[printindex],2))+' & '+str(round(100*abs(10**Comparison_SDSS_sfr[printindex]-10**Comparison_Jimmy_sfr[printindex])/10**Comparison_Jimmy_sfr[printindex],2)))
		#print(Comparison_SDSS_names[printindex]+' & '+str(round(Comparison_SDSS_mass[printindex],2))+' & '+str(round(Comparison_Jimmy_mass[printindex],2))+' & '+str(round(abs(Comparison_SDSS_mass[printindex]-Comparison_Jimmy_mass[printindex]),2))+' & '+str(round(Comparison_SDSS_M_B[printindex],2))+' & '+str(round(Comparison_Jimmy_M_B[printindex],2))+' & '+str(round(abs(Comparison_SDSS_M_B[printindex]-Comparison_Jimmy_M_B[printindex]),2))+' & '+str(round(Comparison_SDSS_metallicity[printindex],2))+' & '+str(round(Comparison_Jimmy_metallicity[printindex],2))+' & '+str(round(abs(Comparison_SDSS_metallicity[printindex]-Comparison_Jimmy_metallicity[printindex]),2))+' & '+str(round(Comparison_SDSS_sfr[printindex],2))+' & '+str(round(Comparison_Jimmy_sfr[printindex],2))+' & '+str(round(abs(Comparison_SDSS_sfr[printindex]-Comparison_Jimmy_sfr[printindex]),2)))
		print(Comparison_SDSS_names[printindex]+' & '+str(round(Comparison_SDSS_mass[printindex],2))+' & '+str(round(Comparison_mpa_jhu_mass[printindex],2))+' & '+str(round(Comparison_Jimmy_mass[printindex],2))+' & '+str(round(Comparison_SDSS_M_B[printindex],2))+' & '+str(round(Comparison_Jimmy_M_B[printindex],2))+' & '+str(round(Comparison_SDSS_metallicity[printindex],2))+' & '+str(round(Comparison_Jimmy_metallicity[printindex],2))+' & '+str(round(Comparison_SDSS_sfr[printindex],2))+' & '+str(round(Comparison_Jimmy_sfr[printindex],2)))
	
	#plot things
	#plt.clf()
	fig = plt.figure(figsize=(11, 10))
	sp1 = fig.add_subplot(2,2,1)
	plt.scatter(Comparison_Jimmy_mass, Comparison_SDSS_mass, marker="*", s=400, color='blue', label='Overlapping')
	plt.errorbar(Comparison_Jimmy_mass, Comparison_SDSS_mass, xerr=Comparison_Jimmy_mass_error, linestyle="None", color='blue')
	plt.scatter(AGC221004_Jimmy_mass, AGC221004_SDSS_mass, marker="*", s=400, color='red', label='AGC221004')
	plt.errorbar(AGC221004_Jimmy_mass, AGC221004_SDSS_mass, xerr=AGC221004_Jimmy_mass_error, linestyle="None", color='red')
	one_to_one = np.linspace(-100, 100, 100)
	plt.plot(one_to_one,one_to_one, label='1:1 Relation')
	plt.legend(loc='upper left', prop={'size':18})
	plt.xlabel(r'log(M$_*$) (IFU)',fontsize=fontsize, labelpad=10)
	plt.ylabel(r'Petrosian Mass log(M$_*$) (SDSS)',fontsize=fontsize, labelpad=10)
	plt.xlim(5,10)
	plt.ylim(5,10)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	sp1 = fig.add_subplot(2,2,2)
	plt.scatter(Comparison_Jimmy_M_B, Comparison_SDSS_M_B, marker="*", s=400, color='blue')
	plt.errorbar(Comparison_Jimmy_M_B, Comparison_SDSS_M_B, xerr=Comparison_Jimmy_M_B_error, linestyle="None", color='blue')
	plt.scatter(AGC221004_Jimmy_M_B, AGC221004_SDSS_M_B, marker="*", s=400, color='red', label='AGC221004')
	plt.errorbar(AGC221004_Jimmy_M_B, AGC221004_SDSS_M_B, xerr=AGC221004_Jimmy_M_B_error, linestyle="None", color='red')
	one_to_one = np.linspace(-100, 100, 100)
	plt.plot(one_to_one,one_to_one, label='1:1 Relation')
	#plt.legend(loc='upper left', prop={'size':18})
	plt.xlabel(r'M$_B$ (IFU)',fontsize=fontsize, labelpad=10)
	plt.ylabel(r'Petrosian Mag M$_B$ (SDSS)',fontsize=fontsize, labelpad=10)
	plt.xlim(-17,-9)
	plt.ylim(-17,-9)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	sp1 = fig.add_subplot(2,2,3)
	plt.scatter(Comparison_Jimmy_metallicity, Comparison_SDSS_metallicity, marker="*", s=400, color='blue')
	plt.errorbar(Comparison_Jimmy_metallicity, Comparison_SDSS_metallicity, xerr=Comparison_Jimmy_metallicity_error, linestyle="None", color='blue')
	plt.scatter(AGC221004_Jimmy_metallicity, AGC221004_SDSS_metallicity, marker="*", s=400, color='red', label='AGC221004')
	plt.errorbar(AGC221004_Jimmy_metallicity, AGC221004_SDSS_metallicity, xerr=AGC221004_Jimmy_metallicity_error, linestyle="None", color='red')
	one_to_one = np.linspace(-100, 100, 100)
	plt.plot(one_to_one,one_to_one, label='1:1 Relation')
	#plt.legend(loc='upper left', prop={'size':18})
	plt.xlabel('IFU 12+log(O/H)',fontsize=fontsize, labelpad=10)
	plt.ylabel('SDSS Fiber 12+log(O/H)',fontsize=fontsize, labelpad=10)
	plt.xlim(7.5,9.0)
	plt.ylim(7.5,9.0)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	sp1 = fig.add_subplot(2,2,4)
	plt.scatter(Comparison_Jimmy_sfr, Comparison_SDSS_sfr, marker="*", s=400, color='blue')
	plt.errorbar(Comparison_Jimmy_sfr, Comparison_SDSS_sfr, xerr=Comparison_Jimmy_sfr_error, linestyle="None", color='blue')
	plt.scatter(AGC221004_Jimmy_sfr, AGC221004_SDSS_sfr, marker="*", s=400, color='red', label='AGC221004')
	plt.errorbar(AGC221004_Jimmy_sfr, AGC221004_SDSS_sfr, xerr=AGC221004_Jimmy_sfr_error, linestyle="None", color='red')
	one_to_one = np.linspace(-100, 100, 100)
	plt.plot(one_to_one,one_to_one, label='1:1 Relation')
	#plt.legend(loc='upper left', prop={'size':18})
	plt.xlabel('IFU log(SFR)',fontsize=fontsize, labelpad=10)
	plt.ylabel('SDSS Fiber log(SFR)',fontsize=fontsize, labelpad=10)
	plt.xlim(-4,-1)
	plt.ylim(-4,-1)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	
	plt.subplots_adjust(bottom=0.13)
	plt.subplots_adjust(left=0.10)
	plt.subplots_adjust(right=0.96)
	plt.subplots_adjust(top=0.98)
	plt.subplots_adjust(wspace=0.31)
	plt.subplots_adjust(hspace=0.31)
	#plt.show()
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
			
	
	#mass = np.array([8.25, 8.50, 8.75, 9.0, 9.25, 9.50, 9.75, 10.0])+0.125
	#metallicity = np.array([8.17, 8.28, 8.35, 8.45, 8.53, 8.64, 8.71, 8.80])
	#SFR = -1.25
	#mass = np.array([8.00, 8.25, 8.50, 8.75, 9.0])+0.125
	#metallicity = np.array([8.25, 8.29, 8.34, 8.41, 8.45])
	#SFR = -2.15
	#mass = np.array([8.75, 9.0, 9.25, 9.50, 9.75, 10.0, 10.25, 10.50, 10.75])+0.125
	#metallicity = np.array([8.24, 8.29, 8.43, 8.56, 8.67, 8.74, 8.79, 8.83, 8.83])
	#SFR = -0.35
	#print(mass-10)
	#print(metallicity)
	#fit_results, fit_error = optimization.curve_fit(second_order_2d, mass-10, metallicity, p0=[8.90,0.37,-0.14,-0.19,0.12,-0.054])
	#print(fit_results)
	#plt.plot(mass, metallicity)
	#x_array = np.linspace(8.0, 10.0, 10)
	#y_array = second_order_2d(x_array-10, *fit_results)
	#plt.plot(x_array, y_array)
	#SFR = -1.25
	#y_array = second_order_2d(x_array-10, *fit_results)
	#plt.plot(x_array, y_array)
	#fit_result, fit_error = optim(mass, SFR, )
	#SDSS_SFR_CAT = np.array([-2.817625, 9999, 9999, -3.403556, -3.195252, 9999, -1.349520, -1.795884, -2.040241, 9999, 9999])
	#IFU_SFR = np.array([-2.84073278, -1.92158181, -2.89500917, -2.86735559, -2.46470433, -2.18502131, -1.78855362, -2.59464928, -3.00786333, -2.43873808, -3.16262098])
	#plt.scatter(SDSS_SFR_CAT, IFU_SFR)
	#plt.plot(np.linspace(-10,10,10), np.linspace(-10,10,10))
	#plt.xlim(-4,1)
	#plt.ylim(-4,1)
	#plt.xlabel('SDSS Fiber Value')
	#plt.ylabel('IFU Value')
	#plt.show()
	#plt.scatter(10**SDSS_SFR_CAT, 10**IFU_SFR)
	#plt.plot(np.linspace(-10,10,10), np.linspace(-10,10,10))
	#plt.xlim(-0.01,0.1)
	#plt.ylim(-0.01,0.1)
	#plt.xlabel('SDSS Fiber Value')
	#plt.ylabel('IFU Value')
	#percent_error = 100*abs((10**IFU_SFR-10**SDSS_SFR_CAT)/10**IFU_SFR)
	#print('Names: '+str(names[np.isfinite(percent_error)]))
	#print('SDSS_SFR_CAT: '+str(SDSS_SFR_CAT[np.isfinite(percent_error)]))
	#print('IFU_SFR: '+str(IFU_SFR[np.isfinite(percent_error)]))
	#percent_error = percent_error[np.isfinite(percent_error)]
	#print('Percent error: '+str(percent_error))
	#plt.show()

if (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr'):
	if (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_sfr_lzr'):
		mode='sfr'
	if (plot_to_make == 'lowest_scatter_hi_mass') or (plot_to_make == 'lowest_scatter_hi_mass_lzr'):
		mode='hi_mass'
	if (plot_to_make == 'lowest_scatter_sfr') or (plot_to_make == 'lowest_scatter_hi_mass'):
		mode2='mass'
	if (plot_to_make == 'lowest_scatter_sfr_lzr') or (plot_to_make == 'lowest_scatter_hi_mass_lzr'):
		mode2='luminosity'
	if (mode2 == 'mass'):
		start_weight = -1.0
		end_weight = 2.0
	if (mode2 == 'luminosity'):
		start_weight = -2.0
		end_weight = 1.0
	number_of_bins = 100
	# define the colormap
	cmap = plt.cm.jet
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# force the first color entry to be grey
	#cmaplist[0] = (.5,.5,.5,1.0)
	#cmaplist[0] = (0.0,0.0,0.0,1.0)
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
	# define the bins and normalize
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	#hi_mass: Minimum scatter: 0.0136222998012 at a weight of : 0.33 (bothwell 0.35)
	#sfr: Minimum scatter: 0.0233810215296 at a weight of : 0.34 (mannucci 0.32)
	#plot_mode='screen'
	plot_mode='file'
	if (plot_mode == 'screen'):
		plt.ion()
	dispersions_array_original_mzr = np.array([])
	dispersions_array_shifting_mzr = np.array([])
	weight_array = np.array([])
	weight = 0.0
	min_scatter = 10
	SDSS_mass = np.append(SDSS_mass, Jimmy_mass)
	SDSS_M_B = np.append(SDSS_M_B, Jimmy_M_B)
	SDSS_metallicity = np.append(SDSS_metallicity, Jimmy_metallicity)
	SDSS_sfr = np.append(SDSS_sfr, Jimmy_sfr)
	SDSS_hi_mass = np.append(SDSS_hi_mass, Jimmy_hi_mass)
	if (mode=='sfr'):
		if (mode2 == 'mass'):
			SDSS_mu = SDSS_mass-(weight*SDSS_sfr)
			mu_bins_lower_limit = np.linspace(mass_bins_lower_limit[0]-(weight*np.mean(SDSS_sfr)), mass_bins_lower_limit[-1]-(weight*np.mean(SDSS_sfr)), 28)
		if (mode2 == 'luminosity'):
			SDSS_mu = SDSS_M_B-(weight*SDSS_sfr)
			mu_bins_lower_limit = np.linspace(M_B_bins_lower_limit[0]-(weight*np.mean(SDSS_sfr)), M_B_bins_lower_limit[-1]-(weight*np.mean(SDSS_sfr)), 28)

	if (mode=='hi_mass'):
		if (mode2 == 'mass'):
			SDSS_mu = SDSS_mass-(weight*(SDSS_hi_mass-9.80))
			mu_bins_lower_limit = np.linspace(mass_bins_lower_limit[0]-(weight*(np.mean(SDSS_hi_mass)-9.80)), mass_bins_lower_limit[-1]-(weight*(np.mean(SDSS_hi_mass)-9.80)), 28)
		if (mode2 == 'luminosity'):
			SDSS_mu = SDSS_M_B-(weight*(SDSS_hi_mass-9.80))
			mu_bins_lower_limit = np.linspace(M_B_bins_lower_limit[0]-(weight*(np.mean(SDSS_hi_mass)-9.80)), M_B_bins_lower_limit[-1]-(weight*(np.mean(SDSS_hi_mass)-9.80)), 28)		
	#mu_bins_lower_limit = np.linspace(0.0, 20.0, 28)
	long_residuals_array = np.array([])
	#if (mode=='sfr'):
		#sfr_bins_lower_limit = np.linspace(-3.5, 1.0, 7)
	#if (mode=='hi_mass'):
		#hi_mass_bins_lower_limit = np.linspace(6.75, 11, 7)
	#establish a zero point for the MZR
	medians = np.array([])
	mus = np.array([])
	for lower_limit in mu_bins_lower_limit:
		upper_limit = lower_limit+(mu_bins_lower_limit[1]-mu_bins_lower_limit[0])
		if (mode=='sfr'):
			for sfr_lower_limit in sfr_bins_lower_limit:
				sfr_upper_limit = sfr_lower_limit+(sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])
				mu_cut = np.logical_and(SDSS_mu>lower_limit,SDSS_mu<upper_limit)
				sfr_cut = np.logical_and(SDSS_sfr>sfr_lower_limit,SDSS_sfr<sfr_upper_limit)
				#print('lower_limit: '+str(lower_limit))
				#print('upper_limit: '+str(upper_limit))
				passed_cut = np.logical_and(mu_cut, sfr_cut)
				#print('len(SDSS_metallicity[mu_cut]): '+str(len(SDSS_metallicity[mu_cut])))
				if len(SDSS_metallicity[passed_cut]) > binning_cut:
					#print('SDSS mu: '+str(SDSS_mu))
					#median = np.median(SDSS_metallicity[passed_cut])
					median = np.mean(SDSS_metallicity[passed_cut])
					medians = np.append(medians, median)
					mus = np.append(mus, lower_limit+((mu_bins_lower_limit[1]-mu_bins_lower_limit[0])/2))
					long_residuals_array = np.append(long_residuals_array, SDSS_metallicity[passed_cut]-median)
		if (mode=='hi_mass'):
			for hi_lower_limit in hi_mass_bins_lower_limit:
				hi_upper_limit = hi_lower_limit+(hi_mass_bins_lower_limit[1]-hi_mass_bins_lower_limit[0])
				mu_cut = np.logical_and(SDSS_mu>lower_limit,SDSS_mu<upper_limit)
				hi_cut = np.logical_and(SDSS_hi_mass>hi_lower_limit,SDSS_hi_mass<hi_upper_limit)
				passed_cut = np.logical_and(mu_cut, hi_cut)
				if len(SDSS_metallicity[passed_cut]) > binning_cut:
					#median = np.median(SDSS_metallicity[passed_cut])
					median = np.mean(SDSS_metallicity[passed_cut])
					medians = np.append(medians, median)
					mus = np.append(mus, lower_limit+((mu_bins_lower_limit[1]-mu_bins_lower_limit[0])/2))
					long_residuals_array = np.append(long_residuals_array, SDSS_metallicity[passed_cut]-median)
	min_mus = min(mus)
	max_mus = max(mus)
	original_mus = mus
	fit_results_original_mzr, fit_error_original_mzr = optimization.curve_fit(fourth_order_fit, mus, medians)
	residuals = medians - fourth_order_fit(mus, *fit_results_original_mzr)
	#print(np.std(residuals))
	#print(np.std(long_residuals_array))
	
	long_residuals_array = np.array([])
	for weight_index, weight in enumerate(np.linspace(start_weight, end_weight, number_of_bins+1)):
		if (mode=='sfr'):
			plot_color = 'black'
			if (mode2 == 'mass'):
				SDSS_mu = SDSS_mass-(weight*SDSS_sfr)
			if (mode2 == 'luminosity'):
				SDSS_mu = SDSS_M_B-(weight*SDSS_sfr)
			#mu_bins_lower_limit = np.linspace(mass_bins_lower_limit[0]-(weight*np.median(SDSS_sfr)), mass_bins_lower_limit[-1]-(weight*np.median(SDSS_sfr)), 28)
			if (mode2 == 'mass'):
				mu_bins_lower_limit = np.linspace(mass_bins_lower_limit[0]-(weight*np.mean(SDSS_sfr)), mass_bins_lower_limit[-1]-(weight*np.mean(SDSS_sfr)), 28)
			if (mode2 == 'luminosity'):
				mu_bins_lower_limit = np.linspace(M_B_bins_lower_limit[0]-(weight*np.mean(SDSS_sfr)), M_B_bins_lower_limit[-1]-(weight*np.mean(SDSS_sfr)), 28)
		if (mode=='hi_mass'):
			plot_color = 'black'
			if (mode2 == 'mass'):
				SDSS_mu = SDSS_mass-(weight*(SDSS_hi_mass-9.80))
			if (mode2 == 'luminosity'):
				SDSS_mu = SDSS_M_B-(weight*(SDSS_hi_mass-9.80))
			if (mode2 == 'mass'):
				mu_bins_lower_limit = np.linspace(mass_bins_lower_limit[0]-(weight*(np.mean(SDSS_hi_mass)-9.80)), mass_bins_lower_limit[-1]-(weight*(np.mean(SDSS_hi_mass)-9.80)), 28)
			if (mode2 == 'luminosity'):
				mu_bins_lower_limit = np.linspace(M_B_bins_lower_limit[0]-(weight*(np.mean(SDSS_hi_mass)-9.80)), M_B_bins_lower_limit[-1]-(weight*(np.mean(SDSS_hi_mass)-9.80)), 28)
		medians = np.array([])
		mus = np.array([])
		sfrs = np.array([])
		hi_masses = np.array([])
		for lower_limit in mu_bins_lower_limit:
			upper_limit = lower_limit+(mu_bins_lower_limit[1]-mu_bins_lower_limit[0])
			if (mode=='sfr'):
				for sfr_lower_limit in sfr_bins_lower_limit:
					sfr_upper_limit = sfr_lower_limit+(sfr_bins_lower_limit[1]-sfr_bins_lower_limit[0])
					mu_cut = np.logical_and(SDSS_mu>lower_limit,SDSS_mu<upper_limit)
					sfr_cut = np.logical_and(SDSS_sfr>sfr_lower_limit,SDSS_sfr<sfr_upper_limit)
					passed_cut = np.logical_and(mu_cut, sfr_cut)
					if len(SDSS_metallicity[passed_cut]) > binning_cut:
						#median = np.median(SDSS_metallicity[passed_cut])
						median = np.mean(SDSS_metallicity[passed_cut])
						medians = np.append(medians, median)
						mus = np.append(mus, lower_limit+((mu_bins_lower_limit[1]-mu_bins_lower_limit[0])/2))
						sfrs = np.append(sfrs, sfr_lower_limit)
						long_residuals_array = np.append(long_residuals_array, SDSS_metallicity[passed_cut]-median)
			if (mode=='hi_mass'):
				for hi_lower_limit in hi_mass_bins_lower_limit:
					hi_upper_limit = hi_lower_limit+(hi_mass_bins_lower_limit[1]-hi_mass_bins_lower_limit[0])
					mu_cut = np.logical_and(SDSS_mu>lower_limit,SDSS_mu<upper_limit)
					hi_cut = np.logical_and(SDSS_hi_mass>hi_lower_limit,SDSS_hi_mass<hi_upper_limit)
					passed_cut = np.logical_and(mu_cut, hi_cut)
					if len(SDSS_metallicity[passed_cut]) > binning_cut:
						#median = np.median(SDSS_metallicity[passed_cut])
						median = np.mean(SDSS_metallicity[passed_cut])
						medians = np.append(medians, median)
						mus = np.append(mus, lower_limit+((mu_bins_lower_limit[1]-mu_bins_lower_limit[0])/2))
						hi_masses = np.append(hi_masses, hi_lower_limit)
						long_residuals_array = np.append(long_residuals_array, SDSS_metallicity[passed_cut]-median)
		#print(medians)
		fit_results, fit_error = optimization.curve_fit(fourth_order_fit, mus, medians)
		residuals = medians - fourth_order_fit(mus, *fit_results)
		print('Working on weight: '+str(weight)+' scatter is: '+str(np.std(residuals)))
		dispersions_array_shifting_mzr = np.append(dispersions_array_shifting_mzr,np.std(residuals))
		if np.std(residuals) < min_scatter:
			min_scatter = np.std(residuals)
			min_weight = weight

		residuals_original_mzr = medians - fourth_order_fit(mus, *fit_results_original_mzr)
		dispersions_array_original_mzr = np.append(dispersions_array_original_mzr,np.std(residuals_original_mzr))
		weight_array = np.append(weight_array, weight)
		plt.scatter(SDSS_mu, SDSS_metallicity, color='gray', alpha=sdss_alpha_value, label='ALFALFA/SDSS')
		x_fit = np.linspace(min(mus), max(mus), 11)
		y_fit = fourth_order_fit(x_fit, *fit_results)
		plt.plot(x_fit, y_fit, color='blue', label='Best Fit')
		
		x_fit = np.linspace(min_mus, max_mus, 11)
		y_fit = fourth_order_fit(x_fit, *fit_results_original_mzr)
		plt.plot(x_fit, y_fit, color='black', label='Original MZR Fit')
		if (mode=='sfr'):
			plt.scatter(mus, medians, c=sfrs, color='black', cmap=cmap, norm=norm, s=100, label='Mean Metallicities')
		if (mode=='hi_mass'):
			plt.scatter(mus, medians, c=hi_masses, color='black', cmap=cmap, norm=norm, s=100, label='Mean Metallicities')
		plt.colorbar()
		plt.errorbar(mus,medians,yerr=np.std(residuals), linestyle="None")
		#plt.text(np.median(SDSS_mu)-2, 9.0, 'Weight being used: '+str(weight))
		#plt.text(np.median(SDSS_mu)-2, 7.0, 'Residual (shifting MZR): '+str(np.std(residuals)))
		#plt.text(np.median(SDSS_mu)-2, 7.25, 'Residual (original MZR): '+str(np.std(residuals_original_mzr)))
		if (mode2 == 'mass'):
			plt.text(6.0, 7.0, r'1$\sigma$ scatter: '+str(round(np.std(residuals),2)),fontsize=fontsize)
		if (mode2 == 'luminosity'):
			plt.text(-14.5, 7.0, r'1$\sigma$ scatter: '+str(round(np.std(residuals),2)),fontsize=fontsize)
		plt.xlim(np.median(SDSS_mu)-3, np.median(SDSS_mu)+3)
		plt.ylim(6.5, 9.5)
		if (mode2 == 'mass'):
			plt.xlim(5.0, 14)
		if (mode2 == 'luminosity'):
			plt.xlim(-14, -22)
		
		plt.tick_params(axis='x', labelsize=20)
		plt.tick_params(axis='y', labelsize=20)	
		plt.minorticks_on()
		plt.legend(loc='lower right')
		plt.subplots_adjust(bottom=bottom_margin)	
		plt.subplots_adjust(left=left_margin)
		plt.subplots_adjust(right=right_margin)
		plt.subplots_adjust(top=top_margin)
		if (mode=='sfr'):
			if (mode2 == 'mass'):
				plt.xlabel('M$_*$-('+str(weight)+'*SFR)',fontsize=fontsize, labelpad=20)
			if (mode2 == 'luminosity'):
				plt.xlabel('M$_B$-('+str(weight)+'*SFR)',fontsize=fontsize, labelpad=20)
		if (mode=='hi_mass'):
			plt.xlabel('M$_*$-('+str(weight)+'*(M$_{HI}$-9.80)',fontsize=fontsize, labelpad=20) #SDSS_mass-(weight*(SDSS_hi_mass-9.80))
		plt.ylabel('12+log(O/H)',fontsize=fontsize, labelpad=20)
		#if (plot_mode == 'screen'):
		#	plt.draw()
		if (plot_mode=='file'):
			plt.savefig(HOME+PAPER_FOLDER+'/weights/'+str(weight_index)+'weight_'+str(weight)+'.png')
		plt.clf()
	#plt.plot(weight_array, dispersions_array_original_mzr, color='red', label='original mzr')
	#plt.plot(weight_array, dispersions_array_shifting_mzr, color='red', linestyle='--', label='shifting mzr')
	plt.plot(weight_array, dispersions_array_shifting_mzr, color=plot_color, label=r'1$\sigma$ dispersion', linewidth=2.0)
	plt.scatter(min_weight, min_scatter, facecolor='none', label=r'Jimmy+15 Minimum', marker='o', s=200, edgecolor='blue', linewidth=3)
	if (mode=='sfr'):
		if (mode2 == 'mass'):
			plt.axvline(x=0.32, color='green', linewidth=2.0, linestyle='--', label="Mannucci+10 Minimum")
			plt.axvline(x=0.28, color='red', linewidth=2.0, linestyle='--', label="Bothwell+13 Minimum")
	if (mode=='hi_mass'):
		if (mode2 == 'mass'):
			plt.axvline(x=0.35, color='red', linewidth=2.0, linestyle='--', label="Bothwell+13 Minimum")
	plt.legend(loc='upper left')
	plt.xlim(start_weight,end_weight)
	plt.ylim(0.0,0.3)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)	
	plt.minorticks_on()
	if (mode=='sfr'):
		if (mode2 == 'mass'):
			plt.xlabel(r'$\alpha$',fontsize=fontsize, labelpad=20)
		if (mode2 == 'luminosity'):
			plt.xlabel(r'$\zeta$',fontsize=fontsize, labelpad=20)
	if (mode=='hi_mass'):
		if (mode2 == 'mass'):
			plt.xlabel(r'$\beta$',fontsize=fontsize, labelpad=20)
		if (mode2 == 'luminosity'):
			plt.xlabel(r'$\gamma$',fontsize=fontsize, labelpad=20)
	plt.ylabel('Dispersion [dex]',fontsize=fontsize, labelpad=20)
	if (plot_mode == 'screen'):
		plt.draw()
		raw_input("Press enter to continue");
	if (plot_mode == 'file'):
		plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	print('Minimum scatter: '+str(min_scatter)+' at a weight of : '+str(min_weight))
	#dummy = input('Press enter to quit')
	
	#stop()
	#for weight in np.linspace(start_weight, 2.0, 101):
	#	SDSS_mu = SDSS_mass-(weight*SDSS_hi_mass-9.80)
	#	#SDSS_mu = SDSS_mass-(weight*SDSS_sfr)
	#	mu_bins_lower_limit = np.linspace(np.median(SDSS_mu)-5, np.median(SDSS_mu)+5, 28)
	#	temp_sdss_mu = np.array([])
	#	temp_sdss_metallicity = np.array([])
	#	#print(np.median(SDSS_mu))
	#	
	#	fit_results_original_mzr, fit_error_original_mzr = optimization.curve_fit(fourth_order_fit, SDSS_mass, SDSS_metallicity)
	#	
	#	
	#	#for bin_limit in mu_bins_lower_limit:
	#	#	counting_array = np.array([])
	#	#	#print('bin limit: '+str(bin_limit))
	#	#	for mu_index, mu in enumerate(SDSS_mu):
	#	#		#if (mu > bin_limit) and (mu < bin_limit+(mu_bins_lower_limit[1]-mu_bins_lower_limit[0])) and (SDSS_metallicity[mu_index] > 7.0) and (SDSS_metallicity[mu_index] < 9.0):
	#	#		if (mu > 0) and (mu < 30) and (SDSS_metallicity[mu_index] > 7.5) and (SDSS_metallicity[mu_index] < 9.0):
	#	#			#temp_sdss_mu = np.append(temp_sdss_mu, SDSS_mu[mu_index])
	#	#			#temp_sdss_metallicity = np.append(temp_sdss_metallicity, SDSS_metallicity[mu_index])
	#	#			counting_array = np.append(counting_array, 1)
	#	#	#print('number of objects in this bin: '+str(len(counting_array)))
	#	#	if len(counting_array) > binning_cut:
	#	#		for mu_index, mu in enumerate(SDSS_mu):
	#	#			if (mu > 0) and (mu < 30) and (SDSS_metallicity[mu_index] > 7.5) and (SDSS_metallicity[mu_index] < 9.0):
	#	#				temp_sdss_mu = np.append(temp_sdss_mu, SDSS_mu[mu_index])
	#	#				temp_sdss_metallicity = np.append(temp_sdss_metallicity, SDSS_metallicity[mu_index])	
	#	#fit_results, fit_error = optimization.curve_fit(fourth_order_fit, temp_sdss_mu-10, temp_sdss_metallicity)
	#	fit_results, fit_error = optimization.curve_fit(fourth_order_fit, SDSS_mu, SDSS_metallicity)
	#	#print('fit_results: '+str(fit_results))
	#	#print('length of the array: '+str(len(temp_sdss_mu)))
	#	x_fit = np.linspace(-10, 50, 101)
	#	y_fit = fourth_order_fit(x_fit, *fit_results)
	#	#residuals = temp_sdss_metallicity - fourth_order_fit(temp_sdss_mu-10, *fit_results)
	#	residuals = SDSS_metallicity - fourth_order_fit(SDSS_mu, *fit_results)
	#	#var_res = np.sum(residuals**2)/(len(SDSS_metallicity)-2)
	#	#sd_res = np.sqrt(var_res)
	#	print('For weight: '+str(weight)+' standard deviation (shifting mzr): '+str(round(np.std(residuals),3)))
	#	dispersions_array_shifting_mzr = np.append(dispersions_array_shifting_mzr,np.std(residuals))
	#	residuals_original_mzr = SDSS_metallicityg - fourth_order_fit(SDSS_mu, *fit_results_original_mzr)
	#	print('For weight: '+str(weight)+' standard deviation (original mzr): '+str(round(np.std(residuals_original_mzr),3)))
	#	dispersions_array_original_mzr = np.append(dispersions_array_original_mzr,np.std(residuals_original_mzr))
	#	weight_array = np.append(weight_array, weight)
	#	plt.plot(x_fit, y_fit, color='red')
	#	#plt.scatter(temp_sdss_mu, temp_sdss_metallicity)
	#	plt.scatter(SDSS_mu, SDSS_metallicity)
	#	plt.xlim(np.median(SDSS_mu)-3, np.median(SDSS_mu)+3)
	#	plt.ylim(6.5, 9.5)
	#	plt.draw()
	#	plt.clf()
	#plt.plot(weight_array, dispersions_array_original_mzr, label='original mzr')
	#plt.plot(weight_array, dispersions_array_shifting_mzr, label='shifting mzr')
	#plt.legend()
	#plt.xlim(-1,2)
	#plt.ylim(0.0,0.3)
	#plt.draw()
	#dummy = input('Press enter to quit')

if (plot_to_make == 'segmentation_v_aon'):
	#names_aon = ['AGC191702', 'AGC212838', 'AGC220755', 'AGC220837', 'AGC221004', 'AGC225882', 'AGC227897']
	#sfr_aon = [-2.82980442, -2.88150391, -2.82393546, -2.4017753, -2.57922886, -2.44758737, -3.15932346]
	sfr_aon = [-2.82980442, -1.88, -2.88150391, -2.82393546, -2.4017753, -2.15, -1.79, -2.57922886, -3.00, -2.44758737, -3.15932346]
 
	#names_segmentation = ['AGC191702', 'AGC212838', 'AGC220755', 'AGC220837', 'AGC221004', 'AGC225882', 'AGC227897']
	#sfr_segmentation = [-2.71789665, -2.72473418, -2.45454168, -2.12113095, -2.23941557, -2.37781596, -2.87604468]
	#sfr_segmentation = [-2.71789665, -2.21430525, -2.72473418, -2.45454168, -2.12113095, -2.45600654, -2.37813348, -2.23941557, -2.67061478, -2.37781596, -2.87604468]
	#names = np.array(['AGC191702', 'AGC202218', 'AGC212838', 'AGC220755', 'AGC220837', 'AGC220860', 'AGC221000', 'AGC221004', 'AGC225852', 'AGC225882', 'AGC227897'])
	sfr_segmentation = [-2.76412536, -1.71675803, -2.72589109, -2.51608047, -2.25476628, -2.14638395, -1.65415823, -2.22718003, -2.52897292, -2.31773698, -2.8848191]
	
	plt.plot([-3.5,-2], [-3.5,-2], color='black', label='1:1 Relation')
	plt.scatter(sfr_aon, sfr_segmentation, s=40, color='blue', label='Jimmy+15')
	plt.plot([-3.5, -2.6], [-2.6,-2.6], color='red', linestyle='--', label='SFR Bin Limits')
	plt.plot([-2.6, -2.6], [-2.6,-3.5], color='red', linestyle='--')
	plt.legend(loc='upper left')
	plt.xlabel(r'AoN 3 SFR',fontsize=fontsize, labelpad=labelpad-20)
	plt.ylabel(r'Segmentation Mask SFR',fontsize=fontsize, labelpad=labelpad-25)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)	
	plt.subplots_adjust(bottom=bottom_margin)	
	plt.subplots_adjust(left=left_margin)
	plt.subplots_adjust(right=right_margin)
	plt.subplots_adjust(top=top_margin)
	plt.savefig(HOME+PAPER_FOLDER+plot_to_make+'.pdf')
	#plt.show()