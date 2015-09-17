#! /usr/bin/env python
#A tool to automatically reduce VIRUS-P data for a given night.

#O3N2 = np.log10((flux_OIII2/flux_Hb)/(flux_NII2/flux_Ha))
#LogOH_PP04_O3N2 = 8.73-(0.32*O3N2)
#N2 = np.log10(flux_NII2/flux_Ha)
#LogOH_PP04_N2 = 9.37 + (2.03*N2) + (1.26*N2**2) + (0.32*N2**3) # FROM 
#LogOH_PP04_N2_linear = 8.90 + 0.57*N2
#LogOH_D02_N2 = 9.12 + 0.73*N2

#pyfits to read in fits files
import pyfits
#numpy is required by pyfits
import numpy as np
#sys to read in the arguments
import sys
#To expand the tilda and define the home directory, also to call other python scripts from within this one
import os
#Import shell utilities for file copy
import shutil
import pylab as plt
import matplotlib as mpl
import scipy.optimize as optimization



plt.rcParams['legend.scatterpoints'] = 1
#Direct input 
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 20,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          'legend.fontsize': 14,
          }
plt.rcParams.update(params) 

def linear_fit(x, a, b):
    return a + b*x

fontsize=24
labelpad = 35

mode='toprow'
mode='tworow'
#mode='bottomrow'

HOME = os.path.expanduser("~")

INPUT_DIR = HOME+'/Astro/reduced/AGC666pro/metallicity_sorted_folders'
sncut = 'sn5'

# define the colormap
cmap = plt.cm.jet
cmaplist = [cmap(i) for i in range(cmap.N)]
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
# define the bins and normalize
bounds = np.linspace(7.4, 8.7, 6)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

recovered_metallicity_PP04_O3N2 = np.array([])
recovered_metallicity_PP04_N2 = np.array([])
recovered_metallicity_D02_N2 = np.array([])
input_metallicity_direct = np.array([])
input_metallicity_PP04_O3N2 = np.array([])
input_metallicity_PP04_N2 = np.array([])
input_metallicity_D02_N2 = np.array([])
Jimmy_metallicity_error_D02_N2 = np.array([])
Jimmy_metallicity_error_PP04_N2 = np.array([])
Jimmy_metallicity_error_PP04_O3N2 = np.array([])

dirs = [d for d in os.listdir(INPUT_DIR) if os.path.isdir(os.path.join(INPUT_DIR, d))]

for folder in dirs:
	if (folder != 'sn3') and (folder != 'sn5'):
		print(folder)
		#input_metallicity = np.append(input_metallicity, float(folder))
		one_bin_table = np.genfromtxt(INPUT_DIR+'/'+folder+'/'+sncut+'/gandalf_table.txt',dtype=None)
		line_of_interest = one_bin_table[0] #0 is for PP04_O3N2, 1 is for PP04_N2, and 2 is for D02_N2
		stacked_PP04_O3N2 = line_of_interest[1]
		line_of_interest = one_bin_table[1] #0 is for PP04_O3N2, 1 is for PP04_N2, and 2 is for D02_N2
		stacked_PP04_N2 = line_of_interest[1]
		line_of_interest = one_bin_table[2] #0 is for PP04_O3N2, 1 is for PP04_N2, and 2 is for D02_N2
		stacked_D02_N2 = line_of_interest[1]
		recovered_metallicity_PP04_O3N2 = np.append(recovered_metallicity_PP04_O3N2, stacked_PP04_O3N2)
		recovered_metallicity_PP04_N2 = np.append(recovered_metallicity_PP04_N2, stacked_PP04_N2)
		recovered_metallicity_D02_N2 = np.append(recovered_metallicity_D02_N2, stacked_D02_N2)
		temp_variable = np.loadtxt(INPUT_DIR+'/'+folder+'/fake_metallicity.sh',dtype=str, delimiter='=')
		input_metallicity_direct = np.append(input_metallicity_direct, float(temp_variable[1]))
		temp_variable = np.genfromtxt(INPUT_DIR+'/'+folder+'/input_PP04_O3N2.txt',dtype=None)
		input_metallicity_PP04_O3N2 = np.append(input_metallicity_PP04_O3N2, temp_variable)
		temp_variable = np.genfromtxt(INPUT_DIR+'/'+folder+'/input_PP04_N2.txt',dtype=None)
		input_metallicity_PP04_N2 = np.append(input_metallicity_PP04_N2, temp_variable)
		temp_variable = np.genfromtxt(INPUT_DIR+'/'+folder+'/input_D02_N2.txt',dtype=None)
		input_metallicity_D02_N2 = np.append(input_metallicity_D02_N2, temp_variable)
		
		Hb_line = one_bin_table[4]
		stacked_Hb = Hb_line[1]
		OIII2_line = one_bin_table[3]
		stacked_OIII2 = OIII2_line[1]
		Ha_line = one_bin_table[6]
		stacked_Ha = Ha_line[1]
		NII2_line = one_bin_table[5]
		stacked_NII2 = NII2_line[1]
		monte_one_file = INPUT_DIR+'/'+folder+'/'+sncut+'/monte_carlo_results_one.txt'
		monte_one_table=np.loadtxt(monte_one_file, dtype=float)
		h_beta_flux_error = np.nanstd(monte_one_table[:,3])
		oiii_5007_flux_error = np.nanstd(monte_one_table[:,7])
		h_alpha_flux_error = np.nanstd(monte_one_table[:,11])
		nii_6584_flux_error = np.nanstd(monte_one_table[:,15])
		N2_error = np.sqrt(((1/(stacked_NII2*np.log(10)))**2)*((nii_6584_flux_error)**2)+((1/(stacked_Ha*np.log(10)))**2)*((h_alpha_flux_error)**2))
		O3N2_error = np.sqrt(((1/(stacked_OIII2*np.log(10)))**2)*((oiii_5007_flux_error)**2)+((1/(stacked_Hb*np.log(10)))**2)*((h_beta_flux_error)**2)+((1/(stacked_NII2*np.log(10)))**2)*((nii_6584_flux_error)**2)+((1/(stacked_Ha*np.log(10)))**2)*((h_alpha_flux_error)**2))
		PP04_O3N2_error = (0.32*O3N2_error)
		D02_N2_error = 0.73*N2_error
		N2 = np.log10(stacked_NII2/stacked_Ha)
		PP04_N2_error = ((2.03) + (2.52*N2) + (0.96*N2**2))*N2_error
		Jimmy_metallicity_error_D02_N2 = np.append(Jimmy_metallicity_error_D02_N2, D02_N2_error)
		Jimmy_metallicity_error_PP04_N2 = np.append(Jimmy_metallicity_error_PP04_N2, PP04_N2_error)
		Jimmy_metallicity_error_PP04_O3N2 = np.append(Jimmy_metallicity_error_PP04_O3N2, PP04_O3N2_error)

#print(input_metallicity_direct)
#print(input_metallicity)
#print(recovered_metallicity_PP04_O3N2)


if (mode == 'tworow'):
	fig = plt.figure(figsize=(15, 10))
if (mode == 'toprow') or (mode == 'bottomrow'):
	fig = plt.figure(figsize=(15, 5))
xaxis = np.array(range(12))
yaxis = xaxis
#input_metallicity = np.array([7.2, 7.53, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2])
#expected_metallicity = input_metallicity
#input_metallicity_direct = input_metallicity_direct/max(input_metallicity_direct)
#for index, metallicity in enumerate(input_metallicity_direct):
#	if metallicity < 6.0:
#		input_metallicity_direct[index] = 7.1

if (mode == 'toprow'):
	sp1 = fig.add_subplot(1,3,1)
if (mode == 'tworow'):
	sp1 = fig.add_subplot(2,3,1)
if (mode != 'bottomrow'):
	input_metallicity = np.array([7.92326459474, 7.95061983246, 8.11421493915, 8.0637913383, 8.10524648292, 8.13012704363, 8.22092459292, 8.29186651711, 8.36239383837, 8.42761414775, 8.45179249432, 8.48006641579, 8.57647048958, 8.63984433731, 8.71373555579, 8.77370799006, 8.80667664955, 8.84932641166])
	plt.errorbar(input_metallicity_PP04_O3N2,recovered_metallicity_PP04_O3N2, yerr=Jimmy_metallicity_error_PP04_O3N2, linestyle="None", color='black')
	plt.scatter(input_metallicity_PP04_O3N2, recovered_metallicity_PP04_O3N2, c=input_metallicity_direct, cmap=cmap, norm=norm, s=60, label='PP04 O3N2', color='black')
	plt.plot(xaxis,yaxis, color='black', label='1:1 Relation')
	plt.xlim(7.4, 8.7)
	plt.ylim(7.4, 8.7)
	plt.legend(loc='lower right')
	plt.xlabel(r'Input Metallicity', labelpad=10)
	plt.ylabel(r'Recovered Metallicity')
	sp1.set_xticklabels(['',7.6,7.8,8.0,8.2,8.4,8.6])

if (mode == 'toprow'):
	sp1 = fig.add_subplot(1,3,2)
if (mode == 'tworow'):
	sp1 = fig.add_subplot(2,3,2)
if (mode != 'bottomrow'):
	#input_metallicity = np.array([7.54282834236, 7.8106094805, 8.01979371317, 8.02291375108, 8.09743011926, 8.12812645902, 8.20229943774, 8.2532421454, 8.31807081267, 8.39412850475, 8.41686998541, 8.44395951335, 8.52572692168, 8.58194915382, 8.68036104051, 8.67678606001, 8.7179453089, 8.79109747193])
	plt.errorbar(input_metallicity_PP04_N2,recovered_metallicity_PP04_N2, yerr=Jimmy_metallicity_error_PP04_N2, linestyle="None", color='black')
	plt.scatter(input_metallicity_PP04_N2, recovered_metallicity_PP04_N2, c=input_metallicity_direct, cmap=cmap, norm=norm, s=60, label='PP04 N2', color='black')
	plt.plot(xaxis,yaxis, color='black', label='1:1 Relation')
	plt.xlim(7.4, 8.7)
	plt.ylim(7.4, 8.7)
	plt.legend(loc='lower right')
	plt.xlabel(r'Input Metallicity', labelpad=10)
	#plt.ylabel(r'Recovered Metallicity')
	sp1.set_yticklabels([])
	sp1.set_xticklabels(['',7.6,7.8,8.0,8.2,8.4,8.6])

if (mode == 'toprow'):
	sp1 = fig.add_subplot(1,3,3)
if (mode == 'tworow'):
	sp1 = fig.add_subplot(2,3,3)
if (mode != 'bottomrow'):
	#input_metallicity = np.array([7.48145197381, 7.67849441801, 7.92524234189, 7.9300364074, 8.05658803495, 8.11445252403, 8.25682642352, 8.34691927629, 8.4461484819, 8.54242430499, 8.56777371904, 8.59625369837, 8.67281953138, 8.71892102365, 8.7900725883, 8.78766607241, 8.81467124057, 8.8592626592])
	plt.errorbar(input_metallicity_D02_N2,recovered_metallicity_D02_N2,yerr=Jimmy_metallicity_error_D02_N2, linestyle="None", color='black')	
	plt.scatter(input_metallicity_D02_N2, recovered_metallicity_D02_N2, c=input_metallicity_direct, cmap=cmap, norm=norm, s=60, label='D02 N2', color='black')
	plt.plot(xaxis,yaxis, color='black', label='1:1 Relation')
	plt.xlim(7.4, 8.7)
	plt.ylim(7.4, 8.7)
	plt.legend(loc='lower right')
	plt.xlabel(r'Input Metallicity', labelpad=10)
	#plt.ylabel(r'Recovered Metallicity')
	sp1.set_yticklabels([])
	cbaxes = fig.add_axes([0.88, 0.14, 0.02, 0.82])
	cb = plt.colorbar(cax=cbaxes)
	cb.ax.tick_params(labelsize=fontsize-6)
	cb.set_label(r'Direct Metallicity', fontsize=fontsize, rotation=270, labelpad=labelpad)
	sp1.set_xticklabels(['',7.6,7.8,8.0,8.2,8.4,8.6])

	plt.subplots_adjust(bottom=0.14)	
	plt.subplots_adjust(left=0.06)
	plt.subplots_adjust(right=0.87)
	plt.subplots_adjust(top=0.95)
	plt.subplots_adjust(wspace=0.0)

if (mode == 'toprow'):
	plt.show()

direct_select = input_metallicity_direct>7.0
if (mode == 'tworow'):
	sp1 = fig.add_subplot(2,3,4)
if (mode == 'bottomrow'):
	#plt.clf()
	sp1 = fig.add_subplot(1,3,1)
input_metallicity = np.array([7.92326459474, 7.95061983246, 8.11421493915, 8.0637913383, 8.10524648292, 8.13012704363, 8.22092459292, 8.29186651711, 8.36239383837, 8.42761414775, 8.45179249432, 8.48006641579, 8.57647048958, 8.63984433731, 8.71373555579, 8.77370799006, 8.80667664955, 8.84932641166])
if (mode == 'tworow'):
	plt.errorbar(input_metallicity_direct,recovered_metallicity_PP04_O3N2, yerr=Jimmy_metallicity_error_PP04_O3N2, linestyle="None", color='black')
	plt.scatter(input_metallicity_direct, recovered_metallicity_PP04_O3N2, c=input_metallicity_direct, cmap=cmap, norm=norm, s=60, label='PP04 O3N2', color='black')
if (mode == 'bottomrow'):
	plt.errorbar(input_metallicity_direct,recovered_metallicity_PP04_O3N2, yerr=Jimmy_metallicity_error_PP04_O3N2, linestyle="None", color='black')
	plt.scatter(input_metallicity_direct, recovered_metallicity_PP04_O3N2, s=60, label='PP04 O3N2', color='black')
plt.plot(xaxis,yaxis, color='black', label='1:1 Relation')
fit_results, fit_error = optimization.curve_fit(linear_fit, input_metallicity_direct[direct_select], recovered_metallicity_PP04_O3N2[direct_select], sigma=Jimmy_metallicity_error_PP04_O3N2[direct_select])
x_fit = np.linspace(min(input_metallicity_direct[direct_select]), max(input_metallicity_direct), 10)
y_fit = linear_fit(x_fit, *fit_results)
n=len(input_metallicity_direct[direct_select])
residuals = recovered_metallicity_PP04_O3N2[direct_select] - linear_fit(input_metallicity_direct[direct_select], *fit_results)
var_res = np.sum(residuals**2)/(n-2)
sd_res = np.sqrt(var_res)
#if (mode == 'tworow'):
	#plt.text(7.6, 8.6, 'Standard Deviation: '+str(round(sd_res,2)))
if (mode == 'bottomrow') or (mode == 'tworow'):
	plt.text(7.5, 8.6, 'Best Fit Slope: '+str(round(fit_results[1],3))+' +\- '+str(round(fit_error[1][1],3)))
plt.plot(x_fit, y_fit, linestyle='--', label='Best Linear Fit', linewidth=2.0)
line_fit = np.polyfit(input_metallicity_direct, recovered_metallicity_PP04_O3N2, 1)
fit = (xaxis*line_fit[0])+line_fit[1]
plt.xlim(7.4, 8.7)
plt.ylim(7.4, 8.7)
plt.legend(loc='lower right')
plt.xlabel(r'Direct Metallicity')
plt.ylabel(r'Recovered Metallicity')
#sp1.set_xticklabels(['',7.6,7.8,8.0,8.2,8.4,8.6])

if (mode == 'tworow'):
	sp1 = fig.add_subplot(2,3,5)
if (mode == 'bottomrow'):
	sp1 = fig.add_subplot(1,3,2)
#input_metallicity = np.array([7.54282834236, 7.8106094805, 8.01979371317, 8.02291375108, 8.09743011926, 8.12812645902, 8.20229943774, 8.2532421454, 8.31807081267, 8.39412850475, 8.41686998541, 8.44395951335, 8.52572692168, 8.58194915382, 8.68036104051, 8.67678606001, 8.7179453089, 8.79109747193])
if (mode == 'tworow'):
	plt.errorbar(input_metallicity_direct,recovered_metallicity_PP04_N2, yerr=Jimmy_metallicity_error_PP04_N2, linestyle="None", color='black')
	plt.scatter(input_metallicity_direct, recovered_metallicity_PP04_N2, c=input_metallicity_direct, cmap=cmap, norm=norm, s=60, label='PP04 N2', color='black')
if (mode == 'bottomrow'):
	plt.errorbar(input_metallicity_direct,recovered_metallicity_PP04_N2, yerr=Jimmy_metallicity_error_PP04_N2, linestyle="None", color='black')
	plt.scatter(input_metallicity_direct, recovered_metallicity_PP04_N2, s=60, label='PP04 N2', color='black')
plt.plot(xaxis,yaxis, color='black', label='1:1 Relation')
line_fit = np.polyfit(input_metallicity_direct, recovered_metallicity_PP04_N2, 1)
fit = (xaxis*line_fit[0])+line_fit[1]
fit_results, fit_error = optimization.curve_fit(linear_fit, input_metallicity_direct[direct_select], recovered_metallicity_PP04_N2[direct_select], sigma=Jimmy_metallicity_error_PP04_N2[direct_select])
x_fit = np.linspace(min(input_metallicity_direct[direct_select]), max(input_metallicity_direct), 10)
y_fit = linear_fit(x_fit, *fit_results)
n=len(input_metallicity_direct[direct_select])
residuals = recovered_metallicity_PP04_N2[direct_select] - linear_fit(input_metallicity_direct[direct_select], *fit_results)
var_res = np.sum(residuals**2)/(n-2)
sd_res = np.sqrt(var_res)
#if (mode == 'tworow'):
#	plt.text(7.6, 8.6, 'Standard Deviation: '+str(round(sd_res,2)))
if (mode == 'bottomrow') or (mode == 'tworow'):
	plt.text(7.5, 8.6, 'Best Fit Slope: '+str(round(fit_results[1],3))+' +\- '+str(round(fit_error[1][1],3)))
plt.plot(x_fit, y_fit, linestyle='--', label='Best Linear Fit', linewidth=2.0)
line_fit = np.polyfit(input_metallicity_direct, recovered_metallicity_PP04_O3N2, 1)
fit = (xaxis*line_fit[0])+line_fit[1]
plt.xlim(7.4, 8.7)
plt.ylim(7.4, 8.7)
plt.legend(loc='lower right')
plt.xlabel(r'Direct Metallicity')
sp1.set_yticklabels([])
#sp1.set_xticklabels(['',7.6,7.8,8.0,8.2,8.4,8.6])
#plt.ylabel(r'Recovered Metallicity')

if (mode == 'tworow'):
	sp1 = fig.add_subplot(2,3,6)
if (mode == 'bottomrow'):
	sp1 = fig.add_subplot(1,3,3)
#input_metallicity = np.array([7.48145197381, 7.67849441801, 7.92524234189, 7.9300364074, 8.05658803495, 8.11445252403, 8.25682642352, 8.34691927629, 8.4461484819, 8.54242430499, 8.56777371904, 8.59625369837, 8.67281953138, 8.71892102365, 8.7900725883, 8.78766607241, 8.81467124057, 8.8592626592])
if (mode == 'tworow'):
	plt.errorbar(input_metallicity_direct,recovered_metallicity_D02_N2,yerr=Jimmy_metallicity_error_D02_N2, linestyle="None", color='black')
	plt.scatter(input_metallicity_direct, recovered_metallicity_D02_N2, c=input_metallicity_direct, cmap=cmap, norm=norm, s=60, label='D02 N2', color='black')
if (mode == 'bottomrow'):
	plt.errorbar(input_metallicity_direct,recovered_metallicity_D02_N2,yerr=Jimmy_metallicity_error_D02_N2, linestyle="None", color='black')
	plt.scatter(input_metallicity_direct, recovered_metallicity_D02_N2, s=60, label='D02 N2', color='black')
plt.plot(xaxis,yaxis, color='black', label='1:1 Relation')
#plt.plot(xaxis,yaxis+0.3, color='black', linestyle='--', label='New Relation')
line_fit = np.polyfit(input_metallicity_direct, recovered_metallicity_D02_N2, 1)
fit = (xaxis*line_fit[0])+line_fit[1]
fit_results, fit_error = optimization.curve_fit(linear_fit, input_metallicity_direct[direct_select], recovered_metallicity_D02_N2[direct_select], sigma=Jimmy_metallicity_error_PP04_O3N2[direct_select])
x_fit = np.linspace(min(input_metallicity_direct[direct_select]), max(input_metallicity_direct), 10)
y_fit = linear_fit(x_fit, *fit_results)
n=len(input_metallicity_direct[direct_select])
residuals = recovered_metallicity_D02_N2[direct_select] - linear_fit(input_metallicity_direct[direct_select], *fit_results)
var_res = np.sum(residuals**2)/(n-2)
sd_res = np.sqrt(var_res)
#if (mode == 'tworow'):
#	plt.text(7.6, 8.6, 'Standard Deviation: '+str(round(sd_res,2)))
if (mode == 'bottomrow') or (mode == 'tworow'):
	plt.text(7.5, 8.6, 'Best Fit Slope: '+str(round(fit_results[1],3))+' +\- '+str(round(fit_error[1][1],3)))
plt.plot(x_fit, y_fit, linestyle='--', label='Best Linear Fit', linewidth=2.0)
line_fit = np.polyfit(input_metallicity_direct, recovered_metallicity_PP04_O3N2, 1)
fit = (xaxis*line_fit[0])+line_fit[1]
plt.xlim(7.4, 8.7)
plt.ylim(7.4, 8.7)
plt.legend(loc='lower right')
plt.xlabel(r'Direct Metallicity')
sp1.set_yticklabels([])
#sp1.set_xticklabels(['',7.6,7.8,8.0,8.2,8.4,8.6])
#plt.ylabel(r'Recovered Metallicity')


for horizontal in np.linspace(8.2, 8.2, 1):
	#print(horizontal-recovered_metallicity_PP04_O3N2)
	n=len(recovered_metallicity_PP04_O3N2)
	residuals = recovered_metallicity_PP04_O3N2 - horizontal
	var_res = np.sum(residuals**2)/(n-2)
	sd_res = np.sqrt(var_res)
	standard_deviation_PP04_O3N2 = round(sd_res,2)
	n=len(recovered_metallicity_PP04_N2)
	residuals = recovered_metallicity_PP04_N2 - horizontal
	var_res = np.sum(residuals**2)/(n-2)
	sd_res = np.sqrt(var_res)
	standard_deviation_PP04_N2 = round(sd_res,2)
	n=len(recovered_metallicity_D02_N2)
	residuals = recovered_metallicity_D02_N2 - horizontal
	var_res = np.sum(residuals**2)/(n-2)
	sd_res = np.sqrt(var_res)
	standard_deviation_D02_N2 = round(sd_res,2)
	print('For a horizontal line at x='+str(horizontal)+' STD PP04 O3N2 = '+str(standard_deviation_PP04_O3N2)+' STD PP04 N2 = '+str(standard_deviation_PP04_N2)+' STD D02 N2 = '+str(standard_deviation_D02_N2))

if (mode == 'bottomrow'):
	plt.subplots_adjust(bottom=0.14)	
	plt.subplots_adjust(left=0.06)
	plt.subplots_adjust(right=0.98)
	plt.subplots_adjust(top=0.98)
	plt.subplots_adjust(wspace=0.0)
	
if (mode == 'tworow') or (mode == 'bottomrow'):
	plt.subplots_adjust(hspace=0.30)
	plt.show()


#marino_table_file = HOME+'/Astro/Catalogs/marino_table.tsv'
#marino_table = np.genfromtxt(marino_table_file, dtype=None, names=['Ref','ID','n_ID', 'O3N2', 'N2', 'metallicity', 'SimbadName', 'RA', 'DEC'], delimiter='\t', skiprows=49)
#plt.scatter(marino_table['N2'],marino_table['metallicity'], color='gray', label='Marino+13')
#x = (np.array(range(100))*0.1)-3
#marino_y = 8.743+(0.462*x)
#plt.plot(x,marino_y, color='green', label='Marino Fit')
#denicolo_y = 9.12+(0.73*x)
#plt.plot(x,denicolo_y, color='blue', label='Denicolo Fit')
#plt.xlim(-2.6, 0.3)
#plt.ylim(7.0, 9.2)
#plt.legend(loc='upper left')
#plt.show()