#!/bin/bash
#Data reduction script written by Jimmy.

#This script takes two inputs, the SDSS-C4 target name of the galaxy, which will be used to create and identify the proper project folder.
#2nd input is the mask to be applied, either all, main, comp, or 2ndcomp
#Called like: vimos_data_reduction.sh 1050 all
#Should only need to edit the directories for the following files, and then most scripts will work when called from here.
#Of course that assumes you're following my naming scheme.

date

#Astro_Dir is the top level directory containing raw and reduced data.
ASTRO_DIR=$HOME/Astro
#Reduced dir is the directory containing the project directories for each galaxy and supporting files
REDUCED_DIR=$ASTRO_DIR/reduced
if [ $1 == "BCG_1261v" ]; then
	REDUCED_DIR=$ASTRO_DIR/reduced/mcdonald
fi
#Pro Dir is the project directory for that particular galaxy.  This is where the outputs are actually sent.
PRO_DIR=$REDUCED_DIR/$1pro
#SOF dir contains the "set of frames" files used by the VIMOS pipeline and is unique to each directory.
SOF_DIR=$REDUCED_DIR/$1pt1sof
#Cal Dir contains the calbration files used by both VIMOS Pipeline and IDL
CAL_DIR=$REDUCED_DIR/cal
#Python directory is where the python scripts are stored
PY_DIR=$ASTRO_DIR/python
#Used to automatically pull the IDL scripts out of testing mode.
export not_testing=1
run_r_e="n"
#low_sn=15
#high_sn=20
#low_sn_half=3
#high_sn_half=10
low_sn=5
high_sn=10
low_sn_half=5
high_sn_half=10



#Take the input from the command line argument.
if [ $# -ne 2 ]; then
        echo "No target given."
        exit 1
fi
echo "Input target number: " $1
echo "Input object of interest: " $2


#If the project directory for this galaxy doesn't exist, we will make it.
if [ ! -d $PRO_DIR ]; then
    mkdir $PRO_DIR
fi
if [ ! -d $PRO_DIR/all ]; then
    mkdir $PRO_DIR/all
fi
if [ ! -d $PRO_DIR/main ]; then
    mkdir $PRO_DIR/main
fi
if [ ! -d $PRO_DIR/comp ]; then
    mkdir $PRO_DIR/comp
fi
if [ ! -d $PRO_DIR/2ndcomp ]; then
    mkdir $PRO_DIR/2ndcomp
fi
export targetsn=$low_sn
if [ ! -d $PRO_DIR/all/sn$targetsn ]; then
    mkdir $PRO_DIR/all/sn$targetsn
    mkdir $PRO_DIR/main/sn$targetsn
    mkdir $PRO_DIR/comp/sn$targetsn
    mkdir $PRO_DIR/2ndcomp/sn$targetsn
fi
export targetsn=$high_sn
if [ ! -d $PRO_DIR/all/sn$targetsn ]; then
    mkdir $PRO_DIR/all/sn$targetsn
    mkdir $PRO_DIR/main/sn$targetsn
    mkdir $PRO_DIR/comp/sn$targetsn
    mkdir $PRO_DIR/2ndcomp/sn$targetsn
fi

#If there is a SOF directory for this galaxy, then we will remove the sym links to whatever galaxy may have been previously used, and add new symlinks for our current directory.
#symlinks are used because it's the easiest way of using one script to do everything without having to refer to a specific project or sof directory every time.

if [ -d $SOF_DIR ]; then
	#if [ $(whoami) != 'jimmyerickson' ]; then
	#if [ $clustermode != 'jimmyerickson' ]; then
    	rm -rf $REDUCED_DIR/pro
    	rm -rf $REDUCED_DIR/sof
    	ln -s $REDUCED_DIR/$1pro $REDUCED_DIR/pro
    	ln -s $REDUCED_DIR/$1sof $REDUCED_DIR/sof
    #fi
else
    echo "Target sof directory not found in reduced data folder."
    exit 1
fi



#QUESTIONS TIME
#This determines which steps in the process will be performed.  The default for all steps is no, so unless a y is entered, that step will be skipped.  They are structured in this manner so that one can just hit enter a bunch of times to skip those steps.

export default="n"

#This sets all the properties for running on the supercomputing cluster where everything runs better automatically.
if [ $clustermode ]; then
	calibration="n"
	science="n"
	idl_spectra_to_cube="n"
	vimos_combine="n"
	skyline_adjust="n"
	mosaic="n"
	mask="n"
	snc="n"
	aon="n"
	bin="n"
	ppxf="n"
	gandalf="y"
	plot="y"
	monte="n"
	lambda="n"
	wiki="n"
	newwiki="n"
else
	read -p "Use esorex to Create Master Bias and Calibration Files?: " -e t1
	if [ -n "$t1" ]; then
		calibration="$t1"
	else
		calibration="$default"
	fi

	read -p "Extract Science Spectra using esorex?: " -e t1
	if [ -n "$t1" ]; then
		science="$t1"
	else
		science="$default"
	fi
	
	read -p "Use esorex to create data cubes?: " -e t1
	if [ -n "$t1" ]; then
		vimos_combine="$t1"
	else
		vimos_combine="$default"
	fi

	read -p "Use IDL to create data cubes?: " -e t1
	if [ -n "$t1" ]; then
		idl_spectra_to_cube="$t1"
	else
		idl_spectra_to_cube="$default"
	fi

	read -p "Use make mosaic to stack cubes?: " -e t1
	if [ -n "$t1" ]; then
		mosaic="$t1"
	else
		mosaic="$default"
	fi

	read -p "Perform mask?: " -e t1
	if [ -n "$t1" ]; then
		mask="$t1"
	else
		mask="$default"
	fi

	read -p "Perform S/N cut?: " -e t1
	if [ -n "$t1" ]; then
		snc="$t1"
	else
		snc="$default"
	fi
	
	read -p "Perform AoN cut?: " -e t1
	if [ -n "$t1" ]; then
		aon="$t1"
	else
		aon="$default"
	fi

	read -p "Bin the data?: " -e t1
	if [ -n "$t1" ]; then
		bin="$t1"
	else
		bin="$default"
	fi

	read -p "Perform PPXF?: " -e t1
	if [ -n "$t1" ]; then
		ppxf="$t1"
	else
		ppxf="$default"
	fi
	
	read -p "Perform TSWIFT?: " -e t1
	if [ -n "$t1" ]; then
		gandalf="$t1"
	else
		gandalf="$default"
	fi

	read -p "Perform Monte Carlo Simulation?: " -e t1
	if [ -n "$t1" ]; then
		monte="$t1"
	else
		monte="$default"
	fi
	
	read -p "Plot Data?: " -e t1
	if [ -n "$t1" ]; then
		plot="$t1"
	else
		plot="$default"
	fi

	read -p "Perform Lambda?: " -e t1
	if [ -n "$t1" ]; then
		lambda="$t1"
	else
		lambda="$default"
	fi
	
	read -p "Old Update the wiki?: " -e t1
	if [ -n "$t1" ]; then
		wiki="$t1"
	else
		wiki="$default"
	fi
	read -p "New Update the wiki?: " -e t1
	if [ -n "$t1" ]; then
		newwiki="$t1"
	else
		newwiki="$default"
	fi

fi

#Pull in environmental variables like redshift and sky fibers
. $SOF_DIR/$1_env.sh


#MASTER BIAS AND CALIBRATION FILES CREATION
#Using the VIMOS pipeline, create the master bias frames that will be used for calibration.
if [ $calibration == y ]; then
	for pointing in 1 2 3 4 5 180; do
		SOF_DIR=$REDUCED_DIR/$1pt${pointing}sof
		if [ -d $SOF_DIR ]; then
			for quadrant in 1 2 3 4; do
				echo " "
				echo "EXECUTING for VIMOS quadrant ${quadrant}: --> esorex vmbias vmbias${quadrant}.sof"
		
				#Not all command line switches have been explored.
			
				esorex --suppress-prefix=true  \
        			--output-dir=$PRO_DIR   \
        			--output-readonly=false \
    	    		--suppress-link=true    \
        			--log-level=off		 \
         			vmbias                  \
        			$SOF_DIR/vmbias${quadrant}.sof
 
				echo " "
				echo "Renaming product to $PRO_DIR/master_bias${quadrant}_pt${pointing}.fits"
				mv $PRO_DIR/master_bias.fits $PRO_DIR/master_bias${quadrant}_pt${pointing}.fits
				echo " "
			done
		fi
	done
	#VIMOS IFU CALIBRATION
	#Creates the calibration files using the master bias above, as well as the other inputs specified in the SOF.
	for pointing in 1 2 3 4 5 180; do
		SOF_DIR=$REDUCED_DIR/$1pt${pointing}sof
		if [ -d $SOF_DIR ]; then
			for quadrant in 1 2 3 4; do
				echo " "
				echo "EXECUTING for VIMOS quadrant ${quadrant}: --> esorex vmifucalib vmifucalib${quadrant}.sof"
  
  				#Not all switches have been explored.
  
				esorex --suppress-prefix=true  \
        			--output-dir=$PRO_DIR     \
        			--output-readonly=false \
        			--suppress-prefix=true  \
        			--suppress-link=true    \
        			--log-level=off		 \
        			vmifucalib              \
        			--LineIdent=FirstGuess       \
        			--CleanBadPixel=true   \
        			$SOF_DIR/vmifucalib${quadrant}.sof

				echo " "
				echo "Renaming product to $PRO_DIR/ifu_arc_spectrum_extracted${quadrant}_pt${pointing}.fits"
				mv $PRO_DIR/ifu_arc_spectrum_extracted.fits $PRO_DIR/ifu_arc_spectrum_extracted${quadrant}_pt${pointing}.fits
				echo " "
				echo "Renaming product to $PRO_DIR/ifu_flat_spectrum_extracted${quadrant}_pt${pointing}.fits"
				mv $PRO_DIR/ifu_flat_spectrum_extracted.fits $PRO_DIR/ifu_flat_spectrum_extracted${quadrant}_pt${pointing}.fits
				echo " "
				echo "Renaming product to $PRO_DIR/ifu_ids${quadrant}_pt${pointing}.fits"
				mv $PRO_DIR/ifu_ids.fits $PRO_DIR/ifu_ids${quadrant}_pt${pointing}.fits
				echo " "
				echo "Renaming product to $PRO_DIR/ifu_master_screen_flat${quadrant}_pt${pointing}.fits"
				mv $PRO_DIR/ifu_master_screen_flat.fits $PRO_DIR/ifu_master_screen_flat${quadrant}_pt${pointing}.fits
				echo " "
				echo "Renaming product to $PRO_DIR/ifu_trace${quadrant}_pt${pointing}.fits"
				mv $PRO_DIR/ifu_trace.fits $PRO_DIR/ifu_trace${quadrant}_pt${pointing}.fits
				echo " "
				echo "Renaming product to $PRO_DIR/ifu_transmission${quadrant}_pt${pointing}.fits"
				mv $PRO_DIR/ifu_transmission.fits $PRO_DIR/ifu_transmission${quadrant}_pt${pointing}.fits
				echo " "
			done
		fi
	done
	rm *.paf
	SOF_DIR=$REDUCED_DIR/$1pt1sof
fi



#VIMOS IFU SCIENCE EXTRACTION
#Create the field of view for each quadrant, and reduce the science frames using the calibration files above.
if [ $science == y ]; then
	for pointing in 1 2 3 4 5 180; do
		SOF_DIR=$REDUCED_DIR/$1pt${pointing}sof
		if [ -d $SOF_DIR ]; then
			for quadrant in 1 2 3 4; do
				for run in a b c d; do
					if [ -e "$SOF_DIR/vmifuscience1${run}.sof" ]; then
						echo " "
						echo "EXECUTING for VIMOS quadrant ${quadrant}: --> esorex vmifuscience  vmifuscience${quadrant}${run}.sof"

						#Not all switches have been explored.

						esorex --suppress-prefix=true  \
							--output-dir=$PRO_DIR   \
		        	 		--output-readonly=false \
         					--suppress-link=true    \
			         		--log-level=off		 \
    			     		vmifuscience            \
        	 				--CleanBadPixel=true   \
		         			--CalibrateFlux=false   \
		         			--UseSkylines=false	\
		         			--UseSkyIndividual=false	\
        		 			$SOF_DIR/vmifuscience${quadrant}${run}.sof

						echo " "
						echo "Renaming product to $PRO_DIR/ifu_fov${quadrant}${run}_pt${pointing}.fits"
						mv $PRO_DIR/ifu_fov.fits $PRO_DIR/ifu_fov${quadrant}${run}_pt${pointing}.fits
						echo " "
						echo "Renaming product to $PRO_DIR/ifu_science_reduced${quadrant}${run}_pt${pointing}.fits"
						mv $PRO_DIR/ifu_science_reduced.fits $PRO_DIR/ifu_science_reduced${quadrant}${run}_pt${pointing}.fits
						echo " "
					fi
				done
			done
		fi
	done
	rm -rf *.paf
fi

if [ $vimos_combine == y ]; then
	#VIMOS COMBINE CUBE
	#Use the Vimos Pipeline to turn the reduced data into the data cubes.  This is probably completely useless because the VIMOS has no procedure for subtracting skylines or identifying dead fibers.
	for pointing in 1 2 3 4 5 180; do
		SOF_DIR=$REDUCED_DIR/$1pt${pointing}sof
		if [ -d $SOF_DIR ]; then
			for quadrant in 1 2 3 4; do
				for run in a b c d; do
					if [ -e "$SOF_DIR/vmifuscience1${run}.sof" ]; then
						echo " "
						echo "EXECUTING: --> esorex vmifucombine vmifucombine${run}.sof"

						#Not all switches have been explored, and probably never will because the VIMOS pipeline is almost entirely, but not quite, completely useless.
  
						esorex --suppress-prefix=true  \
							--output-dir=$PRO_DIR     \
							--output-readonly=false \
							--suppress-prefix=true  \
							--suppress-link=true    \
							--log-level=off		 \
							vmifucombinecube        \
							$SOF_DIR/vmifucombinecube${run}.sof
    
						echo "Renaming product to $PRO_DIR/ifu_science_cube${run}_vimos_pt${pointing}.fits"
						mv $PRO_DIR/ifu_science_cube.fits $PRO_DIR/ifu_science_cube${run}_vimos_pt${pointing}.fits
					fi
				done
			done
		fi
	done
	SOF_DIR=$REDUCED_DIR/$1pt1sof
fi


#if [ $skyline_adjust == y ]; then
#    echo "Fixing the flat field image artifacts"
#    for pointing in 1 2 3 4 5 180; do
#	    	SOF_DIR=$REDUCED_DIR/$1pt${pointing}sof
#			if [ -d $SOF_DIR ]; then
#				for quadrant in 1 2 3 4; do
#					for run in a b c d; do
#						if [ -e "$SOF_DIR/vmifuscience1${run}.sof" ]; then
#							$PY_DIR/fiber_shifter.py $1 ${quadrant} ${run} pt${pointing} noplot
#						fi
#					done
#				done
#			fi
#		done
#fi


#IDL DATA QUADRANT CREATION
#Redo some VIMOS Pipeline stuff in IDL because it's needed for the rest of the IDL scripts.  This takes the reduced spectra from the VIMOS Pipeline and creates a data cube for that quadrant.
if [ $idl_spectra_to_cube == y ]; then
    echo "Using IDL to make the cube quadrants"
    #Take in just the pro dir and cal dir, inside the IDL script the files to be used are defined.
    export infile1=$PRO_DIR/
    export infile2=$CAL_DIR/
    for pointing in 1 2 3 4 5 180; do
    	SOF_DIR=$REDUCED_DIR/$1pt${pointing}sof
    	#if [ $1 == "BCG_1261v" ]; then
    	#	SOF_DIR=$REDUCED_DIR/mcdonald/$1sof
    	#fi
		if [ -d $SOF_DIR ]; then
			if [ $1 == "BCG_1261v" ]; then
    			export quadrant='virus'
    			echo " "
	        	echo "Using IDL to extract spectra from Quadrant ${quadrant} Part ${pointing}"
				for run in a b c d; do
					if [ -e "$PRO_DIR/vp_0089_oextr1.fits" ]; then
			            export fiber_mask=$SOF_DIR/bad_fibers.txt
idl << EOF
.comp spectra_to_cube_quadrant.pro
spectra_to_cube_quadrant,'${quadrant}','${run}','pt${pointing}'
EOF
					fi
				done
    		fi
    		if [ $1 != "BCG_1261v" ]; then
				for quadrant in 1 2 3 4; do
					echo " "
	        		echo "Using IDL to extract spectra from Quadrant ${quadrant} Part ${pointing}"
					for run in a b c d; do
						if [ -e "$SOF_DIR/vmifuscience1${run}.sof" ]; then
				            export fiber_mask=$SOF_DIR/bad_fibers_${quadrant}.txt
idl << EOF
.comp spectra_to_cube_quadrant.pro
spectra_to_cube_quadrant,'${quadrant}','${run}','pt${pointing}'
EOF
						fi
					done
				done
			fi
		fi
	done
	if [ $1 != "BCG_1261v" ]; then
	    for pointing in 1 2 3 4 5 180; do
	    	SOF_DIR=$REDUCED_DIR/$1pt${pointing}sof
			if [ -d $SOF_DIR ]; then
				for run in a b c d; do
					if [ -e "$SOF_DIR/vmifuscience1${run}.sof" ]; then
						echo " "
	        			echo "Combining quadrants into a cube for Run ${run} Pointing ${pointing}"
idl << EOF
.comp combine_cube_quadrants.pro
combine_cube_quadrants,'${run}','pt${pointing}'
EOF
					fi
				done
			fi
		done
		SOF_DIR=$REDUCED_DIR/$1pt1sof
	fi
fi




#MAKE THE FINAL STACKED CUBE USING THE IDL MOSAIC SCRIPT
if [ $mosaic == y ]; then
	#Finding the brightest point to attempt automatic stacking
	echo "Finding brightest pixels"
	
	export infile1=$PRO_DIR
	export outfile1=$SOF_DIR
	
idl << EOF
.comp brightest_point.pro
brightest_point
EOF
	
    #Run the Mosaic procedure
    echo "Mosiac IFU Stacking"
    
    if [ -e "$SOF_DIR/idl.ifudat" ]; then
    	export infile=$SOF_DIR/idl.ifudat
    else
	    export infile=$SOF_DIR/auto.ifudat
    fi
    
    idl $ASTRO_DIR/coms/mosaic.com

    echo "Moving test file to project directory"
    mv test.fits /$PRO_DIR/temp.fits
fi


#MASK TO OBSERVE ONLY THE BCG OR COMPANION
if [ $mask == y ]; then
    export infile1=$PRO_DIR/
    export final_mask=$SOF_DIR/$2_mask.txt
    
    final_mask.py $1 $2
#idl << EOF
#.comp final_mask.pro
#final_mask,'$1','$2'
#EOF
fi

#for targetsn in $low_sn; do
for targetsn in $low_sn $high_sn; do
	#PERFORM THE SIGNAL TO NOISE CUT
	if [ $snc == y ]; then
	    #Perform Signal to Noise Cut
	    echo "S/N Cut"
	    #only input is the test.fits file   
	    #if [ $1 == "BCG_1261v" ]; then
    	#		export quadrant='virus'
    	#fi
	    export infile=$PRO_DIR/$2/$1$2.fits
	    export outfile1=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning.txt
	    if [ "$run_r_e" = y ]; then
	    	export outfile2=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning.txt
	    fi
	    export galaxy=$1
	    if [ $targetsn == $low_sn ]; then
		    export sncut=$low_sn/2
		fi
		if [ $targetsn == $high_sn ]; then
			export sncut=$high_sn/2
		fi
    	# Using the com file allows you to stop and check the s/n cut, EOF requires no user intervention
		#/Applications/itt/idl/idl/bin/idl /Users/jimmy/Astro/coms/sncut.com
idl <<EOF
.comp signal_noise_cut.pro
signal_noise_cut
EOF
		#export outfile1=$PRO_DIR/$2/sn$targetsn/radial_2d_binning.txt
		#export sncut=1.0
#idl <<EOF
#.comp signal_noise_cut.pro
#signal_noise_cut
#EOF
	fi
	
	#PERFORM THE AMPLITUDE TO NOISE CUT
	if [ $aon == y ]; then
	    #Perform Signal to Noise Cut
	    echo "AoN Cut"
	    export infile=$PRO_DIR/$2/$1$2.fits
	    export outfile=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_emission.txt
	    export maskfile=$SOF_DIR/$2_mask.txt
	    export lines_file=$ASTRO_DIR/idl/gandalf/jimmy_lines_emission.txt
	    if [ "$run_r_e" = y ]; then
	    	export outfile2=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning_emission.txt
	    fi
	    export galaxy=$1
	    if [ $targetsn == $low_sn ]; then
		    export sncut=$low_sn_half
		fi
		if [ $targetsn == $high_sn ]; then
			export sncut=$high_sn_half
		fi
		aon.py $infile $2 $outfile $sncut $redshift $lines_file $maskfile
	fi

	if [ $bin == y ]; then
    	#Perform the Binning
    	echo "Voronoi Binning"
    	export infile=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning.txt
    	export outfile1=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_output.txt
    	export outfile2=$PRO_DIR/$2/sn$targetsn/voronoi_2d_bins.txt    
    	vbinning.py $infile $outfile1 $outfile2 $targetsn voronoi
#idl <<EOF
#.comp vbinning.pro
#vbinning
#EOF
		export infile=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_emission.txt
    	export outfile1=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_output_emission.txt
    	export outfile2=$PRO_DIR/$2/sn$targetsn/voronoi_2d_bins_emission.txt    
    	vbinning.py $infile $outfile1 $outfile2 $targetsn voronoi
#idl <<EOF
#.comp vbinning.pro
#vbinning
#EOF

		if [ "$run_r_e" = y ]; then
		    export infile=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning.txt
		    export outfile1=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning_output.txt
	    	export outfile2=$PRO_DIR/$2/sn$targetsn/r_e_2d_bins.txt    
    	vbinning.py $infile $outfile1 $outfile2 $targetsn voronoi
#idl <<EOF
#.comp vbinning.pro
#vbinning
#EOF
		fi

    	export infile=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning.txt
    	export outfile1=$PRO_DIR/$2/sn$targetsn/one_bin_output.txt
    	export outfile2=$PRO_DIR/$2/sn$targetsn/one_bin_bins.txt 
    	vbinning.py $infile $outfile1 $outfile2 $targetsn one
#idl <<EOF
#.comp one_bin.pro
#one_bin
#EOF
		export infile=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_emission.txt
    	export outfile1=$PRO_DIR/$2/sn$targetsn/one_bin_output_emission.txt
    	export outfile2=$PRO_DIR/$2/sn$targetsn/one_bin_bins_emission.txt 
    	vbinning.py $infile $outfile1 $outfile2 $targetsn one
#idl <<EOF
#.comp one_bin.pro
#one_bin
#EOF

		if [ "$run_r_e" = y ]; then
	    	export infile=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning.txt
    		export outfile1=$PRO_DIR/$2/sn$targetsn/r_e_one_bin_output.txt
    		export outfile2=$PRO_DIR/$2/sn$targetsn/r_e_one_bin_bins.txt
    		vbinning.py $infile $outfile1 $outfile2 $targetsn one
#idl <<EOF
#.comp one_bin.pro
#one_bin
#EOF
		fi
	fi

	if [ $ppxf == y ]; then
	    #Find Velocities
	    echo "PPXF Velocities"
	    export ppxfimagesdir=$PRO_DIR/$2/sn$targetsn/ppxf_fits/
	    export infile1=$PRO_DIR/$2/$1$2.fits
	    export infile2=$PRO_DIR/$2/sn$targetsn/one_bin_output.txt
	    export infile3=$ASTRO_DIR/MILES_library
	    export infile4=$PRO_DIR/$2/sn$targetsn/one_bin_bins.txt
	    export outfile=$PRO_DIR/$2/sn$targetsn/ppxf_one_bin_output
	    export lines_file=$ASTRO_DIR/idl/gandalf/jimmy_lines.txt
	    export start_range=300
	    export end_range=1600
	    if [ $1 == "BCG_1261v" ]; then
			export start_range=10
	    	export end_range=950
		fi
	    export template_list="/*"
	    export monte_iterations=0
	    rm -rf $PRO_DIR/$2/sn$targetsn/ppxf_fit_one
	    export ppxfimagesdir=$PRO_DIR/$2/sn$targetsn/ppxf_fit_one/
	    if [ $monte == y ]; then
		    export monte_iterations=100
	    fi
idl << EOF
.comp jimmy_ppxf.pro
jimmy_ppxf
EOF
	    
    	rm -rf $PRO_DIR/$2/sn$targetsn/ppxf_fits
    	export ppxfimagesdir=$PRO_DIR/$2/sn$targetsn/ppxf_fits/
	    export infile2=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_output.txt
	    export infile4=$PRO_DIR/$2/sn$targetsn/voronoi_2d_bins.txt
	    export outfile=$PRO_DIR/$2/sn$targetsn/ppxf_v_bin_output  
idl << EOF
.comp jimmy_ppxf.pro
jimmy_ppxf
EOF
    	rm -rf $PRO_DIR/$2/sn$targetsn/r_e_ppxf_fits
    	export ppxfimagesdir=$PRO_DIR/$2/sn$targetsn/r_e_ppxf_fits/
		if [ "$run_r_e" = y ]; then
	    	export infile2=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning_output.txt
    		export infile4=$PRO_DIR/$2/sn$targetsn/r_e_2d_bins.txt
    		export outfile=$PRO_DIR/$2/sn$targetsn/r_e_ppxf_v_bin_output
idl << EOF
.comp jimmy_ppxf.pro
jimmy_ppxf
EOF
	    	
    		rm -rf /$PRO_DIR/$2/sn$targetsn/r_e_ppxf_fits_one
    		export ppxfimagesdir=$PRO_DIR/$2/sn$targetsn/r_e_ppxf_fits_one/
	    	export infile2=$PRO_DIR/$2/sn$targetsn/r_e_one_bin_output.txt
    		export infile4=$PRO_DIR/$2/sn$targetsn/r_e_one_bin_bins.txt
    		export outfile=$PRO_DIR/$2/sn$targetsn/r_e_ppxf_one_bin_output
idl << EOF
.comp jimmy_ppxf.pro
jimmy_ppxf
EOF
    		
    	fi

		#export infile2=$PRO_DIR/$2/sn$targetsn/radial_2d_binning_output.txt
		#export infile4=$PRO_DIR/$2/sn$targetsn/radial_2d_bins.txt
		#export outfile=$PRO_DIR/$2/sn$targetsn/ppxf_rad_bin_output
		#idl << EOF
		#.comp jimmy_ppxf.pro
		#jimmy_ppxf
		#EOF
		#rm -rf /$PRO_DIR/$2/sn$targetsn/ppxf_fits_rad
		#if [ -d ppxf_fits ]; then
		#mv ppxf_fits /$PRO_DIR/$2/sn$targetsn/ppxf_fits_rad
		#fi
	fi
	
	if [ $gandalf == y ]; then
	    #Find Velocities
	    echo "TSWIFT"
	    
	    export gandalfimagesdir=$PRO_DIR/$2/sn$targetsn/gandalf_one/
	    export infile1=$PRO_DIR/$2/$1$2.fits
	    export infile2=$PRO_DIR/$2/sn$targetsn/one_bin_output_emission.txt
	    export infile3=$ASTRO_DIR/MILES_library
	    export infile4=$PRO_DIR/$2/sn$targetsn/one_bin_bins_emission.txt
	    export outfile=$PRO_DIR/$2/sn$targetsn/gandalf_one.fits
	    export template_list="/*"
	    export lines_file=$ASTRO_DIR/idl/gandalf/jimmy_lines_emission.txt
	    export monte_iterations=0
	    export mdegree=1
	    #rm -rf $gandalfimagesdir
	    if [ $monte == y ]; then
		    export monte_iterations=100
	    fi
#idl <<EOF
#.comp jimmy_gandalf.pro
#jimmy_gandalf
#EOF
		trouble_trouble_trouble.py $1 one sn$targetsn $monte_iterations

	    export infile2=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_output_emission.txt
	    export infile4=$PRO_DIR/$2/sn$targetsn/voronoi_2d_bins_emission.txt
	    export outfile=$PRO_DIR/$2/sn$targetsn/gandalf.fits
	    export gandalfimagesdir=$PRO_DIR/$2/sn$targetsn/gandalf/
	    #rm -rf $gandalfimagesdir
#idl <<EOF
#.comp jimmy_gandalf.pro
#jimmy_gandalf
#EOF
		trouble_trouble_trouble.py $1 voronoi sn$targetsn $monte_iterations

	fi

	if [ $plot == y ]; then
    	#Make Pretty Plots
    	echo "Plot Data"
    	export prodir=$PRO_DIR
    	export infile1=$PRO_DIR/$2/sn$targetsn/one_bin_bins.txt
    	export infile2=$PRO_DIR/$2/sn$targetsn/ppxf_one_bin_output
    	export infile3=$PRO_DIR/$2/sn$targetsn/one_bin_output.txt
idl << EOF
.comp display_data.pro
display_data,'one','$1'
EOF

		export infile1=$PRO_DIR/$2/sn$targetsn/voronoi_2d_bins.txt
    	export infile2=$PRO_DIR/$2/sn$targetsn/ppxf_v_bin_output
    	export infile3=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_output.txt
idl << EOF
.comp display_data.pro
display_data,'vbinned','$1'
EOF
		
	    export infile2=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_output_emission.txt
	    export infile4=$PRO_DIR/$2/sn$targetsn/voronoi_2d_bins_emission.txt
		export lines_file=$ASTRO_DIR/idl/gandalf/jimmy_lines_emission.txt
	    export outfile=$PRO_DIR/$2/sn$targetsn/pink_gandalf.fits
	    export outfile_one=$PRO_DIR/$2/sn$targetsn/pink_gandalf_one.fits
	    export gandalfimagesdir=$PRO_DIR/$2/sn$targetsn/gandalf/
	    export table_one=$PRO_DIR/$2/sn$targetsn/idl_gandalf_table.txt
idl << EOF
.comp super_display_gandalf.pro
super_display_gandalf,'$1'
EOF


		#export infile1=$PRO_DIR/$2/sn$targetsn/radial_2d_bins.txt
		#export infile2=$PRO_DIR/$2/sn$targetsn/ppxf_rad_bin_output
		#export infile3=$PRO_DIR/$2/sn$targetsn/radial_2d_binning_output.txt
#idl << EOF
#.comp display_data.pro
#display_data,'radial','$1'
#EOF

		mv $PRO_DIR/table_one.txt /$PRO_DIR/$2/sn$targetsn/table_one.txt
		mv $PRO_DIR/sigma_scale.eps /$PRO_DIR/$2/sn$targetsn/sigma_scale.eps
		mv $PRO_DIR/sigma.eps /$PRO_DIR/$2/sn$targetsn/sigma.eps
		mv $PRO_DIR/velocity_scale.eps /$PRO_DIR/$2/sn$targetsn/velocity_scale.eps
		mv $PRO_DIR/table.txt /$PRO_DIR/$2/sn$targetsn/table.txt
		#mv $PRO_DIR/sigma_rad.eps /$PRO_DIR/$2/sn$targetsn/sigma_rad.eps
		#mv $PRO_DIR/velocity_rad.eps /$PRO_DIR/$2/sn$targetsn/velocity_rad.eps
		#mv $PRO_DIR/table_rad.txt /$PRO_DIR/$2/sn$targetsn/table_rad.txt
        mv $PRO_DIR/velocity.eps /$PRO_DIR/$2/sn$targetsn/velocity.eps
        mv $PRO_DIR/signal_noise.eps /$PRO_DIR/$2/sn$targetsn/signal_noise.eps
		#mv $PRO_DIR/signal_noise_rad.eps /$PRO_DIR/$2/sn$targetsn/signal_noise_rad.eps

		if [ "$run_r_e" = y ]; then
			export infile1=$PRO_DIR/$2/sn$targetsn/r_e_one_bin_bins.txt
			export infile2=$PRO_DIR/$2/sn$targetsn/r_e_ppxf_one_bin_output
			export infile3=$PRO_DIR/$2/sn$targetsn/r_e_one_bin_output.txt
idl << EOF
.comp display_data.pro
display_data,'one','$1'
EOF

			export infile1=$PRO_DIR/$2/sn$targetsn/r_e_2d_bins.txt
			export infile2=$PRO_DIR/$2/sn$targetsn/r_e_ppxf_v_bin_output
			export infile3=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning_output.txt
idl << EOF
.comp display_data.pro
display_data,'vbinned','$1'
EOF

	        mv $PRO_DIR/table_one.txt /$PRO_DIR/$2/sn$targetsn/r_e_table_one.txt
    	    mv $PRO_DIR/sigma_scale.eps /$PRO_DIR/$2/sn$targetsn/r_e_sigma_scale.eps
    	    mv $PRO_DIR/sigma.eps /$PRO_DIR/$2/sn$targetsn/r_e_sigma.eps
        	mv $PRO_DIR/velocity_scale.eps /$PRO_DIR/$2/sn$targetsn/r_e_velocity_scale.eps
	        mv $PRO_DIR/table.txt /$PRO_DIR/$2/sn$targetsn/r_e_table.txt
    	    mv $PRO_DIR/velocity.eps /$PRO_DIR/$2/sn$targetsn/r_e_velocity.eps
        	mv $PRO_DIR/signal_noise.eps /$PRO_DIR/$2/sn$targetsn/r_e_signal_noise.eps
		fi
	fi

	if [ $lambda == y ]; then
		echo "Running Lambda calculation"
		export indir=$PRO_DIR/$2/sn$targetsn/
		export infile1=$PRO_DIR/$2/sn$targetsn/voronoi_2d_bins.txt
    	export infile2=$PRO_DIR/$2/sn$targetsn/ppxf_v_bin_output
    	export infile3=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning_output.txt
    	export infile4=$PRO_DIR/$2/sn$targetsn/voronoi_2d_binning.txt
    	export outfile1=$PRO_DIR/$2/sn$targetsn/lambda.txt
    	export outfile2=$PRO_DIR/$2/sn$targetsn/lambda_re.txt
		export fitsdir=$PRO_DIR/$2/
idl << EOF
.comp lambda.pro
lambda,'$1','$2'
EOF
		if [ "$run_r_e" = y ]; then
			export infile1=$PRO_DIR/$2/sn$targetsn/r_e_2d_bins.txt
    		export infile2=$PRO_DIR/$2/sn$targetsn/r_e_ppxf_v_bin_output
    		export infile3=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning_output.txt
	    	export infile4=$PRO_DIR/$2/sn$targetsn/r_e_2d_binning.txt
    		export outfile1=$PRO_DIR/$2/sn$targetsn/r_e_lambda.txt
    		export outfile2=$PRO_DIR/$2/sn$targetsn/r_e_lambda_re.txt
			export fitsdir=$PRO_DIR/$2/
idl << EOF
.comp lambda.pro
lambda,'$1','$2'
EOF
		fi
	fi
done

if [ $wiki == y ]; then
	echo "Updating the wiki."
	export infile1=$PRO_DIR/$2/sn$high_sn/lambda_re.txt
	export infile2=$PRO_DIR/$2/sn$high_sn/table_one.txt
	export infile3=$PRO_DIR/$2/sn$low_sn/lambda_re.txt
	export infile4=$PRO_DIR/$2/sn$low_sn/table_one.txt
	export outfile=$PRO_DIR/output.csv
	PASSWORD=$(gpg --decrypt $ASTRO_DIR/wiki.gpg)
idl << EOF
.comp wiki.pro
wiki,'$1','$2'
EOF

	curl -c $PRO_DIR/cookies.txt -d "lgname=Jimmy&lgpassword="$PASSWORD"&action=login&format=xml" https://galaxies.physics.tamu.edu/api.php -o $PRO_DIR/token.txt
	TOKEN=$(cat $PRO_DIR/token.txt | cut -d \" -f 6)

	curl -b $PRO_DIR/cookies.txt -d "action=login&lgname=Jimmy&lgpassword="$PASSWORD"&format=xml&lgtoken="$TOKEN https://galaxies.physics.tamu.edu/api.php
	curl -b $PRO_DIR/cookies.txt -d "action=login&lgname=Jimmy&lgpassword="$PASSWORD"&format=xml&lgtoken="$TOKEN https://galaxies.physics.tamu.edu/api.php
	curl -b $PRO_DIR/cookies.txt -d "action=query&prop=info|revisions&titles="$1"_"$2"&format=xml&intoken=edit" https://galaxies.physics.tamu.edu/api.php > $PRO_DIR/edit.txt

	EDITTOKEN=$(cat $PRO_DIR/edit.txt | awk -F \(edittoken\) '{print $2}')
	EDITTOKEN2=$(echo $EDITTOKEN | cut -c 3-34)
	TIMESTAMP=$(cat $PRO_DIR/edit.txt | cut -d \" -f 22)

	INPUT=$outfile
	OLDIFS=$IFS
	IFS=,
	[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
	while read BCG target a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28 a29 a30 a31 a32 a33
	do
		TEXT="%3d%3d%3dKinematic%20Maps%20S%2fN%20$low_sn%3d%3d%3d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fvelocity_sn$low_sn.jpg%3c%2fimg%3e%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fdispersion_sn$low_sn.jpg%3c%2fimg%3e%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fsn_sn$low_sn.jpg%3c%2fimg%3e%0d%3d%3d%3dKinematic%20Maps%20S%2fN%2010%3d%3d%3d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fvelocity_sn$high_sn.jpg%3c%2fimg%3e%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fdispersion_sn$high_sn.jpg%3c%2fimg%3e%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fsn_sn$high_sn.jpg%3c%2fimg%3e%0d%3d%3d%3dTable%20Data%3d%3d%3d%0d{{BCG_Table|$a1|$a2|$a3|$a4|$a5|$a6|$a7|$a8|$a9|$a10|$a11|$a12|$a13|$a14|$a15|$a16|$a17|$a18|$a19|$a20|$a21|$a22|$a23|$a24|$a25|$a26|$a27|$a28|$a29|$a30|$a31|$a32|$a33}}"
	done < $INPUT
	IFS=$OLDIFS

	curl -b $PRO_DIR/cookies.txt -d "format=xml&action=edit&title="$1"_"$2"&recreate&summary=Automatically%20Generated%20Page&text="$TEXT"&token="$EDITTOKEN2"%2B%5C" https://galaxies.physics.tamu.edu/api.php

	convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$low_sn/velocity_scale.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/velocity_sn$low_sn.jpg
	convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$low_sn/sigma_scale.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/dispersion_sn$low_sn.jpg
	convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$low_sn/signal_noise.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/sn_sn$low_sn.jpg
	convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$high_sn/velocity_scale.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/velocity_sn$high_sn.jpg
	convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$high_sn/sigma_scale.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/dispersion_sn$high_sn.jpg
	convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$high_sn/signal_noise.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/sn_sn$high_sn.jpg	

	if [ "$run_r_e" = y ]; then
		export infile1=$PRO_DIR/$2/sn$high_sn/r_e_lambda_re.txt
		export infile2=$PRO_DIR/$2/sn$high_sn/r_e_table_one.txt
		export infile3=$PRO_DIR/$2/sn$low_sn/r_e_lambda_re.txt
		export infile4=$PRO_DIR/$2/sn$low_sn/r_e_table_one.txt
		export outfile=$PRO_DIR/r_e_output.csv
idl << EOF
.comp wiki.pro
wiki,'$1','$2'
EOF

	
		INPUT=$outfile
		OLDIFS=$IFS
		IFS=,
		[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
		while read BCG target a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 a26 a27 a28 a29 a30 a31 a32 a33
		do
			TEXT="%3d%3d%3dKinematic%20Maps%20S%2fN%205%3d%3d%3d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fr_e_velocity_sn$low_sn.jpg%3c%2fimg%3e%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fr_e_dispersion_sn$low_sn.jpg%3c%2fimg%3e%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fr_e_sn_sn$low_sn.jpg%3c%2fimg%3e%0d%3d%3d%3dKinematic%20Maps%20S%2fN%2010%3d%3d%3d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fr_e_velocity_sn$high_sn.jpg%3c%2fimg%3e%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fr_e_dispersion_sn$high_sn.jpg%3c%2fimg%3e%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2f$BCG%2f$target%2fr_e_sn_sn$high_sn.jpg%3c%2fimg%3e%0d%3d%3d%3dTable%20Data%3d%3d%3d%0d{{BCG_Table|$a1|$a2|$a3|$a4|$a5|$a6|$a7|$a8|$a9|$a10|$a11|$a12|$a13|$a14|$a15|$a16|$a17|$a18|$a19|$a20|$a21|$a22|$a23|$a24|$a25|$a26|$a27|$a28|$a29|$a30|$a31|$a32|$a33}}"
		done < $INPUT
		IFS=$OLDIFS
		
		curl -b $PRO_DIR/cookies.txt -d "format=xml&action=edit&title="$1"_"$2"_r_e&recreate&summary=Automatically%20Generated%20Page&text="$TEXT"&token="$EDITTOKEN2"%2B%5C" https://galaxies.physics.tamu.edu/api.php
	
		convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$low_sn/r_e_velocity_scale.eps" -flatten $PRO_DIR/temp.jpg
		scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/r_e_velocity_sn$low_sn.jpg
		convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$low_sn/r_e_sigma_scale.eps" -flatten $PRO_DIR/temp.jpg
		scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/r_e_dispersion_sn$low_sn.jpg
		convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$low_sn/r_e_signal_noise.eps" -flatten $PRO_DIR/temp.jpg
		scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/r_e_sn_sn$low_sn.jpg
		convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$high_sn/r_e_velocity_scale.eps" -flatten $PRO_DIR/temp.jpg
		scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/r_e_velocity_sn$high_sn.jpg
		convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$high_sn/r_e_sigma_scale.eps" -flatten $PRO_DIR/temp.jpg
		scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/r_e_dispersion_sn$high_sn.jpg
		convert -density 300 "Astro/reduced/"$1"pro/"$2"/sn$high_sn/r_e_signal_noise.eps" -flatten $PRO_DIR/temp.jpg
		scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/$1/$2/r_e_sn_sn$high_sn.jpg	
		
		ssh jimmy@io.physics.tamu.edu chmod a+r -R /home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/

		rm $PRO_DIR/token.txt
		rm $PRO_DIR/edit.txt
		rm $PRO_DIR/output.csv
		rm $PRO_DIR/temp.jpg
		rm $PRO_DIR/cookies.txt
		rm $PRO_DIR/r_e_output.csv
	fi
fi

if [ $newwiki == y ]; then
	echo "Updating the wiki."
	PASSWORD=$(gpg --decrypt $ASTRO_DIR/wiki.gpg)

	curl -c $PRO_DIR/cookies.txt -d "lgname=Jimmy&lgpassword="$PASSWORD"&action=login&format=xml" https://galaxies.physics.tamu.edu/api.php -o $PRO_DIR/token.txt
	TOKEN=$(cat $PRO_DIR/token.txt | cut -d \" -f 6)

	curl -b $PRO_DIR/cookies.txt -d "action=login&lgname=Jimmy&lgpassword="$PASSWORD"&format=xml&lgtoken="$TOKEN https://galaxies.physics.tamu.edu/api.php
	curl -b $PRO_DIR/cookies.txt -d "action=login&lgname=Jimmy&lgpassword="$PASSWORD"&format=xml&lgtoken="$TOKEN https://galaxies.physics.tamu.edu/api.php
	curl -b $PRO_DIR/cookies.txt -d "action=query&prop=info|revisions&titles="$1"_"$2"&format=xml&intoken=edit" https://galaxies.physics.tamu.edu/api.php > $PRO_DIR/edit.txt

	EDITTOKEN=$(cat $PRO_DIR/edit.txt | awk -F \(edittoken\) '{print $2}')
	EDITTOKEN2=$(echo $EDITTOKEN | cut -c 3-34)
	TIMESTAMP=$(cat $PRO_DIR/edit.txt | cut -d \" -f 22)
	TEXT="%3d%3d%3d$1%20S%2fN%205%3d%3d%3d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_logOH_PP04_03N2.jpg%3c%2fimg%3e%0d%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_logOH_PP04_N2.jpg%3c%2fimg%3e%0d%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_logOH_D02.jpg%3c%2fimg%3e%0d%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_PP04_O3N2_gradient.jpg%3c%2fimg%3e%0d%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_PP04_N2_gradient.jpg%3c%2fimg%3e%0d%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_D02_gradient.jpg%3c%2fimg%3e%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_BPT_NII2.jpg%3c%2fimg%3e%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_BPT_SII.jpg%3c%2fimg%3e%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_velocity_Hb.jpg%3c%2fimg%3e%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_dispersion_Hb.jpg%3c%2fimg%3e%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_flux_Ha.jpg%3c%2fimg%3e%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_flux_NII1.jpg%3c%2fimg%3e%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_flux_NII2.jpg%3c%2fimg%3e%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_flux_SII1.jpg%3c%2fimg%3e%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_flux_SII2.jpg%3c%2fimg%3e%0d%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_flux_OII.jpg%3c%2fimg%3e%0d%0d%0d%3cimg%20size%3d360%3ehttps%3a%2f%2fgalaxies.physics.tamu.edu%2fimages%2fjimmy%2fmetallicity%2f$1pro_gandalf_all_0.jpg%3c%2fimg%3e%0d"
	curl -b $PRO_DIR/cookies.txt -d "format=xml&action=edit&title=$1&recreate&summary=Automatically%20Generated%20Page&text=$TEXT&token="$EDITTOKEN2"%2B%5C" https://galaxies.physics.tamu.edu/api.php
	

	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/special_plots/LogOH_PP04_O3N2.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_logOH_PP04_03N2.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/special_plots/LogOH_PP04_N2.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_logOH_PP04_N2.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/special_plots/LogOH_D02.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_logOH_D02.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/special_plots/PP04_O3N2_gradient.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_PP04_O3N2_gradient.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/special_plots/PP04_N2_gradient.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_PP04_N2_gradient.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/special_plots/D02_gradient.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_D02_gradient.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/special_plots/BPT_NII2.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_BPT_NII2.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/special_plots/BPT_SII.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_BPT_SII.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/velocity_Hb.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_velocity_Hb.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/dispersion_Hb.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_dispersion_Hb.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/flux_OII.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_flux_OII.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/flux_Ha.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_flux_Ha.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/flux_NII1.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_flux_NII1.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/flux_NII2.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_flux_NII2.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/flux_SII1.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_flux_SII1.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/flux_SII2.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_flux_SII2.jpg
	convert -density 300 "Astro/reduced/$1pro/$2/sn3/gandalf/gandalf_all/gandalf_all           0.eps" -flatten $PRO_DIR/temp.jpg
	scp $PRO_DIR/temp.jpg jimmy@io.physics.tamu.edu:/home/websites/galaxies.physics.tamu.edu/htdocs/images/jimmy/metallicity/$1pro_gandalf_all_0.jpg
fi

date
