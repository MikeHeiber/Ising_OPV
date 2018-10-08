Examples
=========

## 1) Isotropic Phase-Separation

#### Defining the Parameters

For this example, we will be creating a simple phase-separated morphology with two pure phases using the example parameter file, [parameters_ex1.txt](./example1/parameters_ex1.txt), located in the Ising_OPV/examples/example1 directory.
Let's start by taking a look at the parameters. 

First, we set the General Parameters:
```C++
## General Parameters
50 //Length 
50 //Width 
50 //Height 
true //Enable_z_periodic_boundary 
0.5 //Mix_fraction 
0.4 //Interaction_energy1 
0.4 //Interaction_energy2 
200 //MC_steps 
```
In this section, we set up a 50 x 50 x 50 lattice with periodic boundaries in all directions. 
Periodic boundaries are always enabled in the x- and y- directions, and here we also enable them in the z-direction. 
We then define the composition of the blend by setting the mix volume fraction to 0.5, which is then used to initialize a random blend. 
To control the rate of the phase-separation process we then set the interaction energies to 0.4. 
Large interaction energies will cause slower domain growth, but interaction energies that are too small will not cause enough driving force for phase-separation. 
For simple cases, the interaction energies can should be set to the same value, but they can also be made different so that one could see how this affects the morphology formation. 
Finally, we set the number of Monte Carlo steps to 200 to control the number of site swapping iterations to perform.
Performing more iterations will allow the domains to grow larger and is the primary handle with which to control the final domain size.
However, the number of MC steps needed to reach a particular domain size will also depend on the interaction energies used.

Once the phase-separation process is done, there are several options for further modifying the morphology.
The first of them is smoothing, so let's take a look at the Smoothing Options:
```C++
## Smoothing Options
true //Enable_smoothing 
0.52 //Smoothing_threshold 
```
Without smoothing, the morphologies will often have rough domain interfaces and some island sites.
This is especially true when using smaller interaction energies.
To create idealized model morphologies, one can use the built-in smoothing algorithm.
The smoothing algorithm detects island sites and sites at rough interfaces by calculating the fraction of neighboring sites that are of the opposite type.
When the fraction of dissimilar neighbor sites is above the threshold, the site is changed to the opposite type.
Here we enable smoothing and with the default smoothing threshold of 0.52.

The second option for further modifying the morphology is rescaling, so let's take a look at the Rescale Morphology Options:
```C++
## Rescale Morphology Options
true //Enable_rescale 
2 //Rescale_factor 
false //Enable_shrink 
```
Using the General Parameters specified, we are creating initial morphologies on a 50 x 50 x 50 lattice.
An efficient method for creating morphologies with a larger domain size is to rescale the morphology instead of increasing the MC steps.
Here we enable rescaling with a rescale factor of 2, which will produce morphologies on a 100 x 100 x 100 lattice with domains that are twice as large.

The final option for further modifying the morphology is interfacial mixing, so let's take a look at the Interfacial Mixing Options:
```C++
## Interfacial Mixing Options
false //Enable_interfacial_mixing 
4.0 //Interface_width (nm) 
0.5 //Interface_conc 
```
The interfacial mixing options allow users to create morphologies with controlled interfacial mixing features.
When interfacial mixing is enabled, one can set the interfacial width and the composition of final mixed interfacial regions.
However, in this simple example, we have disabled the interfacial mixing feature.

Now that we have defined the main parameters that control the creation of the morphologies, we must then set the parameters to determine what types of structural analysis will be performed in the Analysis Options section.
```C++
## Analysis Options
false //Enable_analysis_only 
true //Enable_correlation_calc 
100000 //N_sampling_max 
false //Enable_mix_frac_method 
true //Enable_e_method 
false //Enable_extended_correlation_calc 
10 // Extended_correlation_cutoff_distance 
true //Enable_interfacial_distance_calc 
true //Enable_tortuosity_calc 
false //Enable_reduced_memory_tortuosity_calc 
true //Enable_depth_dependent_calc 
true //Enable_areal_maps_calc 
```
When importing morphologies, one can perform analysis only without modifying the imported morphology set, but here we are creating a new morphology set, so this option is disabled.
Here, we enable calculation of the normalized compositional radial autocorrelation data, which is used to determine the domain size.
The calculation is done by averaging the correlation data over a set of randomly sampled sites, and we can define the maximum number of sites to sample.
Here, we choose a maximum of 100,000 sites, so that only 100,000 out of the total 500,000 sites of each type are sampled for the calculation.
To determine the domain size, there are two methods to choose from, the mix fraction method and the 1/e method.
The mix fraction method calculates the domain size by finding where the normalized autocorrelation data first crosses zero.
With some morphologies, this method can fail because the correlation function does not actually clearly cross zero and instead just converges to zero.
As a faster alternative, the 1/e method is the preferred method and set as the default.
The 1/e method calculates the domain size as twice the correlation length, where the correlation length is defined as the distance at which the normalized autocorrelation data first decays to 1/e.
By default, the autocorrelation data will only be calculated out radially as far as needed to determine the domain size.
However, if one would like to output and separately analyze more of the autocorrelation data, one can enable the extended correlation calculation and define the autocorrelation cutoff distance.
Here, we just perform the standard autocorrelation calculation using the 1/e method to determine the domain size.

In addition to domain size calculation options, one can enable or disable several other structural characterization calculations.
One is calculation and output of the interfacial distance probability histograms.
The interfacial distance probability histogram gives data for how close any given site is to the interface.
Another is calculation and output of the end-to-end tortuosity and tortuosity probability histograms.
This tortuosity data helps characterize how convoluted the charge transport pathways are through the film in the z-direction.
In cases where users are generating morphologies on large lattices, the default pathfinding algorithm used by the tortuosity calculation may use up all of the available RAM.
When reaching RAM limitations of the hardware, one can enable a reduced memory algorithm that is significantly slower, but uses much less RAM as a tradeoff.
Another option is the calculation of several of depth dependent characteristics, including the blend composition, domain size, and interfacial volume fraction.
Finally, one can enable or disable the calculation of areal maps for several characteristics, including the composition and tortuosity.
Here, we calculate and output all of the structural analysis metrics.

In addition, let's take a look at some of the more detailed options for modifying the morphology generation process in the Other Options section:
```C++
## Other Options
false //Enable_checkerboard_start 
false //Enable_growth_pref 
3 //Growth_direction 
-0.1 //Additional_interaction 
```
Instead of starting with a random blend, which is the default setting, one can choose to start with a 3D checkerboard morphology instead.
One can also enable a directional preference for domain growth by modifying the interaction energy in one direction. 
Here, we are just doing a standard simple isotropic phase-separation simulation, so we disable these options.

Once all of the options for morphology generation and analysis are set, we then specify how to save the generated morphologies in the Export Morphology Options section:
```C++
true //Enable_export_compressed_files 
false //Enable_export_cross_section 
```
Here, we enable the export of the morphologies in a custom compressed format to save disk space.
This compressed format can be imported later imported into Ising_OPV to further modify the morphologies and can also be directly imported into the [Excimontec](https://github.com/MikeHeiber/Excimontec) kinetic Monte Carlo simulation software.
However, if this is disabled, the morphologies will be saved in a uncompressed format that will be easier for a human to interpret and load into other software packages.
Finally, there is also an option to generate data for a cross-sectional image through the middle of the morphology that will be saved in the uncompressed format.

The last section in the parameter file is the Import Morphology Options:
```C++
false //Enable_import_morphologies (true or false) (choose whether or not to import a previously generated set of morphology files)
false //Enable_import_tomogram (true or false) (choose whether or not to import a tomogram dataset)
TOMO_000001 //Tomogram_name (specify the name of the tomogram dataset) (name cannot contain spaces)
1.0 //Desired_unit_size (nm) (specify the desired lattice unit size to use when importing the tomogram dataset)
0.0 //Mixed_frac (specify the volume fraction of the mixed phase)
0.5 //Mixed_conc (specify the type1 volume fraction within the mixed phase)
16 //N_extracted_segments (specify the number of cuboid segments to extract from the tomogram) (must be 4, 9, 16, 25, 64, 49, etc.)
1 //N_variants (specify the number of random variants to create from each extracted cuboid segment)
```
This section allows users to import a morphology set previously created by Ising_OPV in order to modify it further or perform additional analysis.
In addition, this section is where users set the options for importing and analyzing 3D experimental tomography data.
However, in this simple example, we are only generating a simple morphology set and have disabled all of the import options.

#### Running the Simulation

With all of the parameters set and explained, we now go through executing the simulation. 
Set the working directory to the Ising_OPV/examples/example1 directory and then run the simulation with:
```
mpiexec -n 4 ../../Ising_OPV.exe parameters_ex1.txt
```
If MPI is properly installed and you were successful in building Ising_OPV, the program should immediately start generating command line output informing you of the status of the simulation as it progresses.
If you have problems make sure you re-read the build instructions in the README and report any problems you have in the [Issues](https://github.com/MikeHeiber/Ising_OPV/issues) section.
With the parameters set in the parameters_ex1.txt file, the simulation should only take a few minutes to complete on modern hardware.
Once completed, the simulation should have created a number of new files in the example1 directory.
There should be four morphology files with the compressed morphology data (morphology_0.txt, morphology_1.txt, morphology_2.txt, morphology_3.txt) as well at all of the analysis output files.
A nice quick way to see the characteristics of the morphologies generated is to look at the analysis_summary.txt file.
This file has CSV formatted average metrics for the entire morphology set as well as the individual metrics calculated for each of the individual morphologies.
For this simple example, the final morphologies should have a domain size of around 11.9 nm, an anisotropy value near 1.0, and a tortuosity of about 1.08.
Feel free to explore the other structural analysis data saved in the other files.

#### Concluding Remarks

Well, that's it!  You are now ready to play around with making your own parameter files and adjusting the parameter values to reach your desired morphology.
A good workflow is to play with the parameters and generate small morphology sets on your own machine or on a single interactive node until you finalize a series of parameter files for further testing.
Once you are ready to generate larger morphology sets, it is recommended to submit batch jobs where you will request the resources needed to perform the simulation. 
Remember that one morphology will be generated for each processor that you request.
An example batch script for the SLURM job scheduling system is provided with this package (slurm_script.sh). 
Similar batch scripts can also be written for TORQUE or other job schedulers.