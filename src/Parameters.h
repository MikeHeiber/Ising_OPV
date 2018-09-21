// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Utils.h"
#include <fstream>
#include <stdexcept>

//! \brief This class contains the parameters used by Ising_OPV and functions to import parameter files.
//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
//! \author Michael C. Heiber
//! \date 2014-2018
class Parameters {

public:
	//! \brief Default constructor that creates an empty Parameters object.
	//! \warning An empty Parameters object should not be used until initialized.
	Parameters();

	//bool checkParameters();

	//! \brief Imports the parameters from the input parameter filestream
	//! \param parameterfile is an open input filestream that points to the beginning of a correctly formatted parameter file.
	//! \return true if the parameter file is read and the parameters are successfully set.
	//! \return false if there are any errors during parameter import.
	bool importParameters(std::ifstream& parameterfile);

	// General
	int Length = 0; // x-direction size of the lattice
	int Width = 0; // y-direction size of the lattice
	int Height = 0; // z-direction size of the lattice
	bool Enable_periodic_z = true; // z-direction periodic boundary option
	double Mix_fraction = 0.0; // volume fraction of donor
	double Interaction_energy1 = 0.0; // energetic favorability for type1-type1 interactions over type1-type2 interactions
	double Interaction_energy2 = 0.0; // energetic favorability for type2-type2 interactions over type1-type2 interactions
	int MC_steps = 0; // number of MC steps to be executed (determines number of Ising swapping iterations)
	// Export Morphology Options
	bool Enable_export_compressed_files = true; // choose whether the output morphology data file is in compressed format or not
	bool Enable_export_cross_section = false; // choose whether to output data for an uncompressed cross section (x=0 plane) of the morphology
	// Smoothing Options
	bool Enable_smoothing = false; // choose whether to perform domain smoothing or not
	double Smoothing_threshold = 0.0;  // specify the degree of smoothing
	 // Rescale Options
	bool Enable_rescale = false; // choose whether to perform lattice rescaling or not
	int Rescale_factor = 0; // specify the rescale factor to be used (must be an integer greater than 1)
	bool Enable_shrink = false; // chose whether to shrink the lattice instead of expand it
	// Interfacial Mixing Options
	bool Enable_interfacial_mixing = false; // choose whether to perform interfacial mixing or not
	double Interface_width = 0.0;  // specify the interfacial width
	double Interface_conc = 0.0; // specify the mixing concentration in the interfacial region
	// Analysis Options
	bool Enable_analysis_only = false;  // choose whether to only perform analysis on an imported morphology or not
	bool Enable_correlation_calc = false; // choose whether to perform the domain size calculation using the pair-pair correlation method or not
	int N_sampling_max = 0; // specify the maximum number of sites to be sampled for calculating the correlation function
	bool Enable_mix_frac_method = false; // choose whether to use the mix fraction method for determining the domain size or not
	bool Enable_e_method = false; // choose whether to use the 1/e method for determining the domain size or not
	bool Enable_extended_correlation_calc = false; // choose whether to extend the correlation function calculation to the specified cutoff distance or not
	int Extended_correlation_cutoff_distance = 0; // specify the maximum cutoff distnace for the extended correlation function calculation
	bool Enable_interfacial_distance_calc = false; // choose whether to calculate the interfacial distance histograms or not
	bool Enable_tortuosity_calc = false; // choose whether to calculate the tortuosity end-to-end histograms and island volume fraction or not
	bool Enable_reduced_memory_tortuosity_calc = false; // choose whether to enable a tortuosity calculation method that takes longer, but uses less memory or not
	bool Enable_depth_dependent_calc = false; // choose whether to enable calculation and output of film depth dependent morphology characteristics or not
	bool Enable_areal_maps_calc = false; // choose whether to enable calculation and output of areal mappings of morphology characteristics or not
	// Other Options
	bool Enable_checkerboard_start = false; //choose whether to start from a alternating checkerboard-like configuration instead of a random blend (creates 0.5 mix fraction) or not
	bool Enable_growth_pref = false;
	int Growth_direction = 0;
	double Additional_interaction = 0.0;
	// Tomogram Import Options
	//! Specifies the final unit size (resolution) of the output morphology dataset.
	double Desired_unit_size = 0.0;
	//! Choose whether to enable the image brightness cutoff threshold-based interpretation of the tomography data.
	bool Enable_cutoff_analysis = false;
	//! Specifies the pixel brightness range to assign to a distinct mixed third phase.
	int Mixed_greyscale_width = 0;
	//! Specifies the volume fraction of the mixed third phase.
	double Mixed_conc = 0.0;
	//! Choose whether to enable a probability-based analysis of pixel brightness for interpreting the tomography data. 
	bool Enable_probability_analysis = false;
	//! Specifies the probability scaling exponent use by the probability-based pixel brightness analysis method.
	double Probability_scaling_exponent = 0.0;
	//! Specify the number of equal size cuboids to extract from the tomogram.
	int N_extracted_segments = 0;
	//! Specify the number of random variants to create from each extracted cuboid segment.
	int N_variants = 0;

protected:

private:


};

#endif // PARAMETERS_H
