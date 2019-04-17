// Copyright (c) 2014-2019 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Utils.h"
#include "Version.h"
#include <fstream>
#include <stdexcept>

namespace Ising_OPV {

	//! \brief This class contains the parameters used by Ising_OPV and functions to import parameter files.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2014-2019
	class Parameters {

	public:
		//! \brief Default constructor that creates an empty Parameters object.
		//! \warning An empty Parameters object should not be used until initialized.
		Parameters();

		//! \brief Checks the current values of all parameters to ensure that they are valid.
		//! \return true if all parameters are valid.
		//! \return false if any parameter value is invalid.
		bool checkParameters() const;

		//! \brief Imports the parameters from the input parameter filestream
		//! \param parameterfile is an open input filestream that points to the beginning of a correctly formatted parameter file.
		//! \return true if the parameter file is read and the parameters are successfully set.
		//! \return false if there are any errors during parameter import.
		bool importParameters(std::ifstream& parameterfile);

		// General
		//! x-direction size of the lattice
		int Length = 0;
		//! y-direction size of the lattice
		int Width = 0;
		//! z-direction size of the lattice
		int Height = 0; 
		//! Choose whether or not to enable z-direction periodic boundary conditions
		bool Enable_periodic_z = true; 
		//! volume fraction of type1 sites
		double Mix_fraction = 0.0; 
		//! energetic favorability for type1-type1 interactions over type1-type2 interactions
		double Interaction_energy1 = 0.0;
		//! energetic favorability for type2-type2 interactions over type1-type2 interactions
		double Interaction_energy2 = 0.0; 
		//! number of MC steps to be executed (determines number of Ising swapping iterations), which sets the duration of the phase separation process
		int MC_steps = 0; 
		// Smoothing Options
		//! choose whether or not to perform domain smoothing
		bool Enable_smoothing = false;
		//! cutoff threshold for the smoothing algorithm
		double Smoothing_threshold = 0.0;
		 // Rescale Options
		 //! choose whether or nmot to perform lattice rescaling
		bool Enable_rescale = false; 
		//! rescale factor to be used
		int Rescale_factor = 0;
		//! chose whether or not to shrink the lattice by 1/rescale_factor instead of expand it 
		bool Enable_shrink = false;
		// Interfacial Mixing Options
		//! choose whether or not to perform interfacial mixing
		bool Enable_interfacial_mixing = false;
		//! interfacial width of the mixed interface
		double Interface_width = 0.0;
		//! concentration of type1 sites in the interfacial mixed region
		double Interface_conc = 0.0;
		// Analysis Options
		//! choose whether or not to only perform analysis on an imported morphology (no modification of the morphology)
		bool Enable_analysis_only = false;
		//! choose whether or not to perform the domain size calculation using the pair-pair autocorrelation function method
		bool Enable_correlation_calc = false; 
		//! maximum number of sites to be sampled for calculating the autocorrelation function
		int N_sampling_max = 0;
		//! choose whether or not to use the mix fraction method for determining the domain size from the autocorrelation data
		bool Enable_mix_frac_method = false;
		//! choose whether or not to use the 1/e method for determining the domain size from the autocorrelation data
		bool Enable_e_method = false;
		//! choose whether or not to extend the autocorrelation function calculation to the specified cutoff distance
		bool Enable_extended_correlation_calc = false;
		//! cutoff distance for the extended autocorrelation function calculation
		int Extended_correlation_cutoff_distance = 0;
		//! choose whether or not to calculate the interfacial distance histograms
		bool Enable_interfacial_distance_calc = false;
		//! choose whether or not to calculate the end-to-end tortuosity histograms and island volume fraction
		bool Enable_tortuosity_calc = false;
		//! choose whether or not to perform the tortuosity calculation using an algorithm that takes longer, but uses less memory
		bool Enable_reduced_memory_tortuosity_calc = false;
		//! choose whether or not to calculate and output the film depth dependent morphology characteristics
		bool Enable_depth_dependent_calc = false;
		//! choose whether or not to calculate and output areal mappings of the morphology characteristics
		bool Enable_areal_maps_calc = false;
		// Other Options
		//! choose whether or not to start from a alternating checkerboard-like configuration instead of a random blend (creates 0.5 mix fraction)
		bool Enable_checkerboard_start = false;
		//! choose whether or not to implement a directional-dependent interaction energy that caues directional-dependent domain growth
		bool Enable_growth_pref = false;
		//! direction that has the modified interaction energy with 1=x,2=y,3=z
		int Growth_direction = 0;
		//! amount that the interaction energy is modified by in the specified direction
		double Additional_interaction = 0.0;
		// Export Morphology Options
		//! choose whether or not the output morphology data file is in compressed format
		bool Enable_export_compressed_files = true;
		//! choose whether or not to output uncompressed data for a cross-section (x=Length/2 plane) of the morphology
		bool Enable_export_cross_section = false;
		// Import Morphology Options
		//! choose whether or not the import a morphology set from the working directory
		bool Enable_import_morphologies = false;
		//! choose whether or not to import a tomogram dataset from the working directory
		bool Enable_import_tomogram = false;
		//! name of the tomogram dataset corresponding to the Tomogram_name.xml and Tomogram_name.raw files
		std::string Tomogram_name = "";
		//! desired unit size (resolution) of the output morphology dataset extracted from the tomogram
		double Desired_unit_size = 0.0;
		//! volume fraction of the mixed phase in the tomogram dataset
		double Mixed_frac = 0.0;
		//! volume fraction of type1 sites in the mixed phase of the tomogram dataset
		double Mixed_conc = 0.0;
		//! number of equal size cuboids to extract from the tomogram dataset
		int N_extracted_segments = 0;
		//! number of random variants to create from each extracted cuboid segment
		int N_variants = 0;

	protected:

	private:

	};
}

#endif // PARAMETERS_H
