# Arvert 7

# Description

Arvert is a numerical model that inverts geochronological data for thermal history. The current version can use blends of MDD, U-Th/He with radiation damage, and general volume-diffusion data. Up to ten samples can be jointly inverted, with up to 50 mineral ages and one MDD age spectrum per sample. Arvert can contextually place the samples into any 1D profile using temperature offsets between samples. Arvert uses a combination of the controlled random search algorithm and a custom search algorithm to explore parameter space.

# Getting Started

### Nature of Arvert

Arvert is plain-vanilla C++ code that you operate from the command line. I now use it under macOS running natively on Apple Silicon. Output from Arvert is formatted to produce graphics using gmt 6 (Generic Mapping Tools). It *should* run under Linux and Windows but you will need to configure your environment to run the gmt script.

### Dependencies and requirements

1. **gmt5 or 6**. Arvert is far more useful if you can visualize its results. Having a complete installation of gmt 6 (or 5) is strongly recommended. Packages for a variety of operating systems can be found at the [gmt download page](https://www.generic-mapping-tools.org/download/). If you dislike gmt it should be possible to script other software to plot results, using either Arvert's default output files or getting under the hood and tweaking what Arvert writes for output.

2. **C++ compiler**

### Installing

You will need to compile and link the Arvert source code for your machine's hardware. 

* Note that the header files `arvert.h, nr3.h, ran_pz.h, and tchar.h` all need to be in the source directory along with the actual source files.

There's no need for a makefile, so in a console like Apple's Terminal.app, just navigate to the directory containing the source code and type:

    g++ ZRDAAM.cpp main.cpp averages.cpp fits.cpp heatsched.cpp mineral.cpp lablovera.cpp geolovera.cpp monteTt.cpp newhistoryCRS.cpp ran3.cpp selectsubset.cpp sort.cpp -o arvert7 -Ofast -march=native -w

Use g++ for the Gnu compiler, clang++ if you prefer just Apple's command-line tools, or the command for whatever C++ compiler you have. For macOS, convenient frequently updated installer packages for g++ can be found at [https://hpc.sourceforge.net](https://hpc.sourceforge.net). Or, you can use a package management system like Homebrew.

### Executing

On a *nix system, just enter ./arvert7 and reply to the prompts to get going. However...

It's easiest to place the compiled executable in a directory that contains all the data files for the sample set you are inverting. If all goes well, Arvert will generate a RESULTS-xxx output directory containing both text output files as well as plotting files; the plotting information and rendered plot will be in the subdirectory PLOTTING, with the rendered plot being named arvert_results,pdf. If you run successive modeling attempts, Arvert will generate new numbered output directories -- nothing will be overwritten.

## Help

You will need to consult the Arvert 7 manual to be able to operate the code effectively. Remember to look at both the tT results AND the fits to the measured data! And remember that the results of any inversion will be conditioned by both the explicit and implicit constraints you supply (intentionally or otherwise).

### Examples

This directory contains several sets of input data for various combinations of samples as well as output results. You should try to run one or more of these if you think you installed Arvert properly but your input files might be causing Arvert to throw errors or complaints.

## Author

Peter Zeitler  
EES Department  
Lehigh University  
Bethlehem, PA 18015 USA

[peter.zeitler@lehigh.edu](mailto:peter.zeitler@lehigh.edu?subject=Arvert7%20question)

### Other indirect and thus innocent contributors

The Lovera() routine for calculating <sup>40</sup>Ar/<sup>39</sup>Ar age spectra is a rewritten C++ version of Oscar Lovera's original fortran code. The RDAAM/ZrDAAM routine for apatite and zircon helium ages is Rich Ketcham's code (circa HeFTY version ~1.7). Neither Oscar or Rich are to blame for any errors I introduced.
