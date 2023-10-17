#!/bin/bash

# Version 7.0.1.1 to go with Arvert 7.0.1 most recent (Toffset version)
# January 23, 2022
# October 15 2023 bugfix to handle  proper use of final -P -O code at final gmt call

# Written by pkz

# This script plots output from Arvert 7 using gmt; obviously gmt 5 or gmt 6 needs to be installed.
# If run manually this script needs to be in the same directory as all the output files.

# As coded this is for gmt 5 but the script will run under gmt 6.

# I've tried to provide some documentation but keep in mind that the gmt calls are often mixed with bash
# variables and calls to awk, head, tail, and bc. Also, the script uses associative arrays, which have a complex
# syntax. Finally, the data files this script uses were generated within Arvert, creating a symbiosis - to fully
# understand what is going on you will likely need to look at the Arvert output code to understand what is in
# each file.

# If you have lots of samples and minerals, this script is slow - be patient for a minute or three while it churns.
# This particularly applies to using the "full" plotting options, which displays the initial MC history and also
# all instances of pool results for minerals and age spectra.

# gmt setup
gmt gmtset PS_CHAR_ENCODING ISOLatin1+   # easier to get symbols if you use ISOLATIN1+
gmt gmtset FONT_ANNOT_PRIMARY 10p,Helvetica,black
gmt gmtset FONT_LABEL 12p,Helvetica,black
gmt gmtset MAP_GRID_PEN_PRIMARY 0.20p,100/100/100,2_1:0p

# some bash variables to control parts of gmt commands
target="arvert_results.ps"  #  output file name
start="-P -K"
control="-P -K -O"
ending="-P -O"

# we assume we are in PLOTTING directory and that the needed plot files are local in this directory

# get needed info on plotting limits
# note that some data are pre-adjusted in Arvert (given 0-based arrays and loops) to minimize need for bash math

head -n 1 limits.in > scratch.out
read -r -d '' modelduration maxtemp agelimit nmembers use_minerals nworst nbest poolsize offset gmtplot < scratch.out

if ((use_minerals < 2)); then  # we have spectra somewhere so we need to learn about them
	tail -n 1 limits.in > scratch.out
	read -ra nspectra < scratch.out
else  # having no spectra we have no info so we need to pad out nspectra[] for later logic
	for ((i = 0 ; i <= $nmembers ; i++)); do
		nspectra[$i]=0
	done
fi

if ((use_minerals > 0)); then
	read -ra nmineral < member_minbest.in
fi

# more variables to fill with derived info from limits.in
general=$(expr $poolsize - $nbest - $nworst)
fgeneral=$(expr $poolsize - $nbest - $nworst + 1)

read -ra olimits < Toffset_worst.in

####################################################
##### 1. CRS tT PLOT for reference sample
####################################################
project=-JX3.4i/2.7i
range=-R0/$modelduration/0/$maxtemp

# tT plot backdrop
gmt psxy backdrop.in -A -Xf4.4i -Yf8.5i $project -L -G255/255/255 $range -P -K > $target # backdrop for tT plot

if [ "$gmtplot" -eq "2" ];  # for full plot option
then
	# tT plot range constrained by minerals
	if (($use_minerals > 0)); then
		gmt psxy mconstraint.in -A -L -G210/230/255 $project $range $control >> $target # mineral constrained range 
	fi

	# tT plot range constrained by age spectra
	if (($use_minerals < 2)); then
		gmt psxy sconstraint.in -A -L -G210/255/230 $project $range $control >> $target # spectra constrained range
	fi

	# tT plot MC starting pool
	for ((i = 0 ; i < $poolsize ; i++)); do
		b=$((2+$i))
		awk "{print \$1,\$$b}" montetT.in | gmt psxy -W0.3p,225/225/225 -A $project $range $control >> $target
	done

	# tT plot CRS general members
	for ((i = 0 ; i < $general ; i++)); do
		b=$((2+$i))
		awk "{print \$1,\$$b}" pool_general.in | gmt psxy -A -W0.4p,150/150/240 $project $range $control >> $target 
	done
fi

# tT plot CRS worst
for ((i = 0 ; i < $nworst ; i++)); do
	b=$((2+$i))
	awk "{print \$1,\$$b}" pool_worst.in | gmt psxy -A -W0.6p,red $project $range $control >> $target 
done

# tT plot CRS best
for ((i = 0 ; i < $nbest ; i++)); do
	b=$((2+$i))
	awk "{print \$1,\$$b}" pool_best.in | gmt psxy -A -W0.6p,0/150/85 $project $range $control >> $target 
done

# tT plot CRS lower envelope
gmt psxy envelope_low.in -A -W2.0p,black $project $range $control >> $target

# tT plot CRS upper envelope
gmt psxy envelope_high.in -A -W2.0p,black $project $range $control >> $target

# tT plot frame and axes - using hardwired plot scalings
if (( (modelduration > 0) && (modelduration <= 40) ));
then
gmt psbasemap -Bf1a5g5:'Time (Ma)':/a100f50g100:'Temperature (\260C)':WSne $project $range $control >> $target
fi
if (( (modelduration > 40) && (modelduration < 100) ));
then
gmt psbasemap -Bf2a10g10:'Time (Ma)':/a100f50g100:'Temperature (\260C)':WSne $project $range $control >> $target
fi
if (( (modelduration > 100) && (modelduration < 250) ));
then
gmt psbasemap -Bf5a25g25:'Time (Ma)':/a100f50g100:'Temperature (\260C)':WSne $project $range $control >> $target
fi
if (( (modelduration > 250) && (modelduration < 500) ));
then
gmt psbasemap -Bf10a50g50:'Time (Ma)':/a100f50g100:'Temperature (\260C)':WSne $project $range $control >> $target
fi
if ((modelduration >= 500));
then
gmt psbasemap -Bf100a500g100:'Time (Ma)':/a100f50g100:'Temperature (\260C)':WSne $project $range $control >> $target
fi

####################################################
# END OF tT plot - reference
####################################################

########################################################################
##### 2. Toffset plot 
########################################################################
# reduce size of some things
gmt gmtset FONT_ANNOT_PRIMARY 8p
gmt gmtset FONT_LABEL 10p
gmt gmtset MAP_GRID_PEN_PRIMARY 0.15p,100/100/100,2_1:0p

range=-R0/100/${olimits[$((0))]}/${olimits[$((1))]}
project=-JX2.2i/1.0i

o_range=$(( ${olimits[$((1))]} - ${olimits[$((0))]} ))

if (( o_range <= 20 ));
then
gmt psbasemap  -Xr-3.15i -Yr1.7i  -Bf100::/f5a5g5:'Toffset (\260C)':WSne $project $range $control >> $target
fi

if (( (o_range > 20) && (o_range <= 50) ));
then
gmt psbasemap  -Xr-3.15i -Yr1.7i  -Bf100::/f5a10g10:'Toffset (\260C)':WSne $project $range $control >> $target
fi

if (( (o_range > 50) && (o_range <= 200) ));
then
gmt psbasemap  -Xr-3.15i -Yr1.7i  -Bf100::/f25a25g25:'Toffset (\260C)':WSne $project $range $control >> $target
fi

if (( o_range > 200 ));
then
gmt psbasemap  -Xr-3.15i -Yr1.7i  -Bf100::/f50a100g100:'Toffset (\260C)':WSne $project $range $control >> $target
fi

#Toffsets - worst members
head -n $((nworst + 1)) Toffset_worst.in | tail -n $nworst > scratch.out
for ((i = 0 ; i < $nworst; i++)); do
	a=$((2*$i + 1))
	b=$((2*$i + 2))
	awk "{print \$$a,\$$b}" scratch.out | gmt psxy -A  -Sc0.08i -W0.6p,red $project $range $control >> $target
done	

# Toffsets - best members
head -n $((nbest + 1)) Toffset_best.in | tail -n $nbest > scratch.out
for ((i = 0 ; i < $nworst; i++)); do
	a=$((2*$i + 1))
	b=$((2*$i + 2))
	awk "{print \$$a,\$$b}" scratch.out | gmt psxy -A  -Sc0.08i -W0.6p,0/150/85 $project $range $control >> $target
done	

# label age plot
gmt pstext Toffset_labels.in -F+f9,,black+jcB -Gwhite $project $range $control >> $target

########################################################################
##### 3. Age Plot for reference sample (spectrum and/or minerals)
########################################################################

range=-R0/100/0/$agelimit
project=-JX1.5i/1.5i

# backdrop for reference age plot
gmt psxy age_backdrop.in -A -Xr0i -Yr-1.7i -L -G255/255/255 $project $range $control >> $target

## spectrum first
if (( (use_minerals < 2) && (${nspectra[0]} > 0) ));  # this means we have MDD data for the first sample
then
	if [ "$gmtplot" -eq "2" ];  # for full plot option
	then
		# model spectra - general members
		for ((i = 0 ; i < $general ; i++)); do
			b=$((2+$i))
			awk "{print \$1,\$$b}" reference_spec_general.in | gmt psxy -A -W0.4p,150/150/240 $project $range $control >> $target
		done
	fi
	
	# model spectra - worst members
	for ((i = 0 ; i < $nworst; i++)); do
		b=$((2+$i))
		awk "{print \$1,\$$b}" reference_spec_worst.in | gmt psxy -A -W0.6p,red $project $range $control >> $target
	done

	# model spectra - best members
	for ((i = 0 ; i < $nbest; i++)); do
		b=$((2+$i))
		awk "{print \$1,\$$b}" reference_spec_best.in | gmt psxy -A -W0.6p,0/150/85 $project $range $control >> $target
	done
	#awk -v x="1" -v y=$(expr $nbest + 1) 'BEGIN{OFS=IFS="\t"} {print $x,$y}' spectrum_best.in | gmt psxy -A -W0.5p,red $project $range $control >> $target

	# goal spectrum 
	gmt psxy goal_spectrum.in -A $project $range $control >> $target
fi

## minerals next
if (( (use_minerals > 0) ));  # this means we have mineral-age data 
then
	if (( (${nmineral[0]} > 0) )); # this means we have mineral data for this specific sample
	then
		# measured (goal) ages
		gmt psxy reference_mingoal.in -A -Sc0.14i -W0.4p,blue -Gblue $project $range $control >> $target

		# CRS-predicted worst ages
		for ((i = 0; i < $nworst; i++)); do
			a=$((2*$i + 1))
			b=$((2*$i + 2))
			awk "{print \$$a,\$$b}" reference_minworst.in | gmt psxy -A -Sc0.09i -W0.6p,red $project $range $control >> $target
		done

		# CRS-predicted best ages
		for ((i = 0; i < $nbest; i++)); do
			a=$((2*$i + 1))
			b=$((2*$i + 2))
			awk "{print \$$a,\$$b}" reference_minbest.in | gmt psxy -A -Sc0.09i -W0.6p,0/150/85 $project $range $control >> $target
		done
	fi
fi

# label axes depending on whether we have any MDD data ± mineral ages, or just mineral ages
if (( (use_minerals < 2) && (${nspectra[0]} > 0) ));
then
	# axes for sample having 4039 MDD data

	# age plot frame and axes
	if (( (agelimit > 0) && (agelimit <= 30) ));
	then
	gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f1a5g5:'Age (Ma)':WSne $project $range $control >> $target
	fi
	if (( (agelimit > 30) && (agelimit < 100) ));
	then
	gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f2a10g10:'Age (Ma)':WSne $project $range $control >> $target
	fi
	if (( (agelimit > 100) && (agelimit < 250) ));
	then
	gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f5a25g25:'Age (Ma)':WSne $project $range $control >> $target
	fi
	if (( (agelimit > 250) && (agelimit < 500) ));
	then
	gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f10a50g50:'Age (Ma)':WSne $project $range $control >> $target
	fi
	if ((agelimit >= 500));
	then
	gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f100a500g500:'Age (Ma)':WSne $project $range $control >> $target
	fi
else
	# axes if we just have mineral age data
	if (( (agelimit > 0) && (agelimit <= 30) ));
	then
	gmt psbasemap -Bf100:'Mineral':/f1a5g5:'Age (Ma)':WSne $project $range $control >> $target
	fi
	if (( (agelimit > 30) && (agelimit < 100) ));
	then
	gmt psbasemap -Bf100:'Mineral':/f2a10g10:'Age (Ma)':WSne $project $range $control >> $target
	fi
	if (( (agelimit > 100) && (agelimit < 250) ));
	then
	gmt psbasemap -Bf100:'Mineral':/f5a25g25:'Age (Ma)':WSne $project $range $control >> $target
	fi
	if (( (agelimit > 250) && (agelimit < 500) ));
	then
	gmt psbasemap -Bf100:'Mineral':/f10a50g50:'Age (Ma)':WSne $project $range $control >> $target
	fi
	if ((agelimit >= 500));
	then
	gmt psbasemap -Bf100:'Mineral':/f100a500g500:'Age (Ma)':WSne $project $range $control >> $target
	fi
fi

# label age plot
if ((nmembers < 1));
then
	echo 7 90 "s1 (ref)" | gmt pstext -F+f9,,black+jLB  $project -R0/100/0/100 $ending >> $target
else
	echo 7 90 "s1 (ref)" | gmt pstext -F+f9,,black+jLB  $project -R0/100/0/100 $control >> $target
fi

####################################################
### DONE WITH REFERENCE PLOTS
####################################################

########################################################################
##### 4. Age plots for member samples (spectrum and/or minerals)
########################################################################

# reduce size of some things
gmt gmtset FONT_ANNOT_PRIMARY 8p
gmt gmtset FONT_LABEL 10p
gmt gmtset MAP_GRID_PEN_PRIMARY 0.15p,100/100/100,2_1:0p

# some coordinate variables needed for automating subplot location
xleft=1.25
ytop=8.10 
xoffset=2.5
yoffset=2.0

project=-JX1.5i/1.5i

# read some needed info about the numbers and nature of data in files
if ((  (use_minerals < 2) && (nmembers > 0) )); then
	read -ra specsteps < member_specbest.in
	read -ra gspeclines < mgoallines.in
fi

# Counters needed for parsing complex files.
# Member plots could have variable amounts of data
# depending on the number of samples, presence of an MDD spectrum, and/or variable number of minerals
# writing a variable number of output files would be too awkward and produce too many files, so I decided to use a single
# data file for each main component (model spectra, goal spectrum, mineral ages both predicted and measured). This
# required some info to be either prefixed to a file's first line, or sent as a separate file.
# Note that the goal age spectra are particularly trying because depending on how many steps are skipped.
# the complex formatting file can be quite variable and long.
# In general we handle these complex files using a combination of head and tail commands plus an accumulating counter
# to snip out parts of the files for placing in a scratch file to be plotted by gmt. We also have to account for member
# samples that gave no mineral and/or no MDD spectrum data
# Thus, the need for... counters for parsing complex files

lines_goal=0 # cumulative number of lines in goal spectrum member file
specsum=0 # cumulative number of age spectra steps in member file
minsum=0 # cumulative number of minerals in member file
member=0 # id of which sample is being parsed and plotted

# Member plots will be in a 3*3 grid that can hold the maximum number of nine.
# so, we have to step through rows and columns but we also have to plot
# only the actual members, which will usually be well fewer than 9! Hence the if statement.
for ((yrow = 0; yrow <= 2; yrow++)); do
	for ((xcol = 0; xcol <= 2; xcol++)); do
		member=$(((yrow + 0) * 3 + (xcol + 1)))
		if (($member <= $nmembers));
		then
			if ((use_minerals > 0)); then
				minsum=$((minsum + ${nmineral[$((member))]}))
			fi

			if ((use_minerals < 2)); then			
				specsum=$((specsum + ${specsteps[$((member))]}))
				lines_goal=$((lines_goal + ${gspeclines[$((member))]}))
			fi
			
			# plot backdrop for age plot
			xposn=$(bc <<< "scale=2; $xleft + $xcol*$xoffset")i
			yposn=$(bc <<< "scale=2; $ytop - ($yrow + 1)*$yoffset")i
			gmt psxy age_backdrop.in -A -Xf$xposn -Yf$yposn -L -G255/255/255 $project $range $control >> $target

			# plot member model spectrum - best members only for now
			if (( (use_minerals < 2) && (${nspectra[member]} > 0) )); # this sample has a spectrum we are using
			then
				if [ "$gmtplot" -eq "2" ];  # for full plot option
				then
					# nworst worst-fit age spectra for this member
					head -n $((specsum + 1)) member_specworst.in | tail -n ${specsteps[$((member))]} > scratch.out
					for ((i = 0 ; i < $nbest; i++)); do
						b=$((2+$i))
						awk "{print \$1,\$$b}" scratch.out | gmt psxy -A -W0.6p,red $project $range $control >> $target
					done
				fi
			
				# nbest best-fit age spectra for this member
				head -n $((specsum + 1)) member_specbest.in | tail -n ${specsteps[$((member))]} > scratch.out
				for ((i = 0 ; i < $nbest; i++)); do
					b=$((2+$i))
					awk "{print \$1,\$$b}" scratch.out | gmt psxy -A -W0.6p,0/150/85 $project $range $control >> $target
				done
			
				# goal spectrum for this member
				head -n $lines_goal member_specgoal.in | tail -n ${gspeclines[$((member))]} > scratch.out
				gmt psxy scratch.out -A $project $range $control >> $target
			fi
	
			if (( use_minerals > 0 ));
			then			
				for ((n = 0; n < ${nmineral[$((member))]}; n++)); do
					head -n $minsum member_mingoal.in | tail -n ${nmineral[$((member))]} > scratch.out
					gmt psxy scratch.out -A -Sc0.14i -W0.4p,blue -Gblue $project $range $control >> $target
				done 
				
				for ((n = 0; n < ${nmineral[$((member))]}; n++)); do
					head -n $((1 + member * nworst)) member_minworst.in | tail -n $nworst > scratch.out
					# predicted ages
					for ((i = 0; i < $nworst; i++)); do
						a=$((2*$i + 1))
						b=$((2*$i + 2))
						awk "{print \$$a,\$$b}" scratch.out | gmt psxy -A -Sc0.09i -W0.6p,red $project $range $control >> $target
					done
				done				
				
				for ((n = 0; n < ${nmineral[$((member))]}; n++)); do
					head -n $((1 + member * nbest)) member_minbest.in | tail -n $nbest > scratch.out
					# predicted ages
					for ((i = 0; i < $nbest; i++)); do
						a=$((2*$i + 1))
						b=$((2*$i + 2))
						awk "{print \$$a,\$$b}" scratch.out | gmt psxy -A -Sc0.09i -W0.6p,0/150/85 $project $range $control >> $target
					done
				done
			fi
	
			# age plot frame and axes
			if (( (use_minerals < 2) && (${nspectra[member]} > 0) ));
#			if ((use_minerals < 2));
			then
				if (((agelimit > 0) && (agelimit <= 30) ));
				then
				gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f1a5g5:'Age (Ma)':WSne $project $range $control >> $target
				fi
				if (((agelimit > 30) && (agelimit < 100) ));
				then
				gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f2a10g10:'Age (Ma)':WSne $project $range $control >> $target
				fi
				if (((agelimit > 100) && (agelimit < 250) ));
				then
				gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f5a25g25:'Age (Ma)':WSne $project $range $control >> $target
				fi
				if (((agelimit > 250) && (agelimit < 500) ));
				then
				gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f10a50g50:'Age (Ma)':WSne $project $range $control >> $target
				fi
				if ((agelimit >= 500));
				then
				gmt psbasemap -Bf5a20g10:'@+39@+Ar Loss (%)':/f100a500g500:'Age (Ma)':WSne $project $range $control >> $target
				fi
			else
				if (((agelimit > 0) && (agelimit <= 30) ));
				then
				gmt psbasemap -Bf100:'Mineral':/f1a5g5:'Age (Ma)':WSne $project $range $control >> $target
				fi
				if (((agelimit > 30) && (agelimit < 100) ));
				then
				gmt psbasemap -Bf100:'Mineral':/f2a10g10:'Age (Ma)':WSne $project $range $control >> $target
				fi
				if (((agelimit > 100) && (agelimit < 250) ));
				then
				gmt psbasemap -Bf100:'Mineral':/f5a25g25:'Age (Ma)':WSne $project $range $control >> $target
				fi
				if (((agelimit > 250) && (agelimit < 500) ));
				then
				gmt psbasemap -Bf100:'Mineral':/f10a50g50:'Age (Ma)':WSne $project $range $control >> $target
				fi
				if ((agelimit >= 500));
				then
				gmt psbasemap -Bf100:'Mineral':/f100a500g500:'Age (Ma)':WSne $project $range $control >> $target
				fi
			fi
			if (($member < $nmembers));
			then
				echo 7 90 "s"$((member+1)) | gmt pstext -F+f9,,black+jLB  $project -R0/100/0/100 $control >> $target
			else
				echo 7 90 "s"$((member+1)) | gmt pstext -F+f9,,black+jLB  $project -R0/100/0/100 $ending >> $target		
			fi
		fi # check for valid member
	done
done

# Could convert $target .ps file to pdf using Ghostscript if it's installed
# ps2pdf $target 

# Probably better to use gmt for conversion, like this:
gmt psconvert $target -Tf -Z

# rm $target   # need this if we use ps2pdf

open arvert_results.pdf