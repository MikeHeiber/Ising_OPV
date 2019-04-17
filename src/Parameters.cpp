// Copyright (c) 2014-2019 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "Parameters.h"

using namespace std;

namespace Ising_OPV {

	Parameters::Parameters() {

	}

	bool Parameters::checkParameters() const {
		bool Error_found = false;
		// Check for valid lattice dimensions
		if (Length <= 0 || Width <= 0 || Height <= 0) {
			cout << "Parameter error!  The input Length, Width, and Height of the lattice must be greater than zero." << endl;
			Error_found = true;
		}
		// Check the input mix fraction
		if (Mix_fraction < 0 || Mix_fraction > 1) {
			cout << "Parameter error! The input Mix_fraction must be between 0 and 1." << endl;
			Error_found = true;
		}
		// Check the input interaction energies
		if (Interaction_energy1 < 0 || Interaction_energy2 < 0) {
			cout << "Parameter error! The input Interaction_energy1 and Interaction_energy2 parameters cannot be negative." << endl;
			Error_found = true;
		}
		// Check the input number of Monte Carlo steps
		if (MC_steps < 0) {
			cout << "Parameter error! The input MC_steps parameter cannot be negative." << endl;
			Error_found = true;
		}
		// Check the smoothing parameters
		if (Enable_smoothing && !(Smoothing_threshold > 0)) {
			cout << "Parameter error! When performing smoothing, the input Smoothing_threshold must be greater than zero." << endl;
			Error_found = true;
		}
		// Check the lattice resclae parameters
		if (Enable_rescale && !(Rescale_factor > 0)) {
			cout << "Parameter error! When rescaling the lattice, the input Rescale_factor must be greater than zero." << endl;
			Error_found = true;
		}
		if (Enable_rescale && Enable_shrink && (Length % Rescale_factor != 0 || Width % Rescale_factor != 0 || Height % Rescale_factor != 0)) {
			cout << "Parameter error! When shrinking the lattice, the input Rescale_factor must be an integer multiple of the Length, Width, and Height." << endl;
			Error_found = true;
		}
		// Check the interfacial mixing parameters
		if (Enable_interfacial_mixing && !(Interface_width > 0)) {
			cout << "Parameter error! When performing interfacial mixing, the input Interface_width must be greater than zero." << endl;
			Error_found = true;
		}
		if (Enable_interfacial_mixing && (!(Interface_conc > 0) || !(Interface_conc < 1))) {
			cout << "Parameter error! When performing interfacial mixing, the input Interface_conc must be greater than zero and less than 1." << endl;
			Error_found = true;
		}
		// Check the correlation calculation parameters
		if (Enable_correlation_calc && !(N_sampling_max > 0)) {
			cout << "Parameter error! When performing the correlation calculation, the mix fraction method and the 1/e method cannot both be enabled." << endl;
			Error_found = true;
		}
		if (Enable_correlation_calc && Enable_mix_frac_method && Enable_e_method) {
			cout << "Parameter error! When performing the correlation calculation, the mix fraction method and the 1/e method cannot both be enabled." << endl;
			Error_found = true;
		}
		if (Enable_correlation_calc && !Enable_mix_frac_method && !Enable_e_method) {
			cout << "Parameter error! When performing the correlation calculation, either the mix fraction method or the 1/e method must be enabled." << endl;
			Error_found = true;
		}
		if (Enable_correlation_calc && Enable_extended_correlation_calc && !(Extended_correlation_cutoff_distance > 0)) {
			cout << "Parameter error! When performing the extended correlation calculation, Extended_correlation_cutoff_distance must be greater than zero." << endl;
			Error_found = true;
		}
		// Check the growth preference parameters
		if (Enable_growth_pref && (Growth_direction < 1 || Growth_direction > 3)) {
			cout << "Parameter error! When performing phase separation with a directional growth preference, the input Growth_direction paramter must be 1, 2, or 3." << endl;
			Error_found = true;
		}
		if (Enable_growth_pref && !(Additional_interaction > 0) && !(Additional_interaction < 0)) {
			cout << "Parameter error! When performing phase separation with a directional growth preference, the input Additional_interaction parameter must not be zero." << endl;
			Error_found = true;
		}
		// Check tomogram import parameters
		if (Enable_import_morphologies && Enable_import_tomogram) {
			cout << "Parameter error! The import morphologies and import tomogram options cannot both be enabled." << endl;
			Error_found = true;
		}
		if (Enable_import_tomogram && !(Desired_unit_size > 0)) {
			cout << "Parameter error! When importing a tomogram dataset, the input Desired_unit_size must not be zero." << endl;
			Error_found = true;
		}
		if (Enable_import_tomogram && !(Mixed_frac < 1)) {
			cout << "Parameter error! When importing a tomogram dataset, the Mixed_frac must be graeter than equal to 0 and less than 1." << endl;
			Error_found = true;
		}
		if (Enable_import_tomogram && (!(Mixed_conc > 0) || !(Mixed_conc < 1))) {
			cout << "Parameter error! When importing a tomogram dataset, the Mixed_conc must be greater than zero and less than 1." << endl;
			Error_found = true;
		}
		//if (Enable_import_tomogram && !Enable_cutoff_analysis && !Enable_probability_analysis) {
		//	cout << "Parameter error! When importing a tomogram dataset, the cutoff analysis or the probability analysis option must be enabled." << endl;
		//	Error_found = true;
		//}
		//if (Enable_import_tomogram && Enable_cutoff_analysis && Enable_probability_analysis) {
		//	cout << "Parameter error! When importing a tomogram dataset, the cutoff analysis and the probability analysis options cannot both be enabled." << endl;
		//	Error_found = true;
		//}
		//if (Enable_import_tomogram && Enable_probability_analysis && Probability_scaling_exponent < 0) {
		//	cout << "Parameter error! When importing a tomogram dataset and using the probability analysis option, the Probability_scaling_exponent must not be negative." << endl;
		//	Error_found = true;
		//}
		if (Enable_import_tomogram && N_extracted_segments <= 0) {
			cout << "Parameter error! When importing a tomogram dataset, the input N_extracted segments must be greater than zero." << endl;
			Error_found = true;
		}
		if (Enable_import_tomogram && N_variants <= 0) {
			cout << "Parameter error! When importing a tomogram dataset, the input N_variants segments must be greater than zero." << endl;
			Error_found = true;
		}
		if (Enable_import_tomogram && N_extracted_segments != intpow((int)floor(sqrt(N_extracted_segments)), 2)) {
			cout << "Parameter error! When importing a tomogram dataset, the input value for N_extracted_segments must be 1, 4, 9, 16, 25, 36, 49, 64, 89, 100, 121, 144, 169, or 196 but " << N_extracted_segments << " was entered." << endl;
			Error_found = true;
		}
		// Check other parameter conflicts
		if (Enable_analysis_only && !Enable_import_morphologies && !Enable_import_tomogram) {
			cout << "Parameter error!  The 'analysis only' option can only be used when importing morphologies." << endl;
			Error_found = true;
		}
		if (Error_found) {
			return false;
		}
		return true;
	}

	bool Parameters::importParameters(ifstream& parameterfile) {
		Version min_version("4.0.0-rc.1");
		string line;
		// Parse header file line
		getline(parameterfile, line);
		line = line.substr(line.find('v') + 1);
		Version file_version;
		try {
			file_version = Version(line);
		}
		catch (invalid_argument exception) {
			cout << "Error! Unable to load parameter file with version " << line << ". Only parameter files formatted for Ising_OPV v4.0.0-rc.1 or greater are supported." << endl;
			return false;
		}
		if (file_version < min_version) {
			cout << "Error! Unable to load parameter file with v" << file_version << ". Only parameter files formatted for Ising_OPV v4.0.0-rc.1 or greater are supported." << endl;
			return false;
		}
		string var;
		vector<string> stringvars;
		// Read input file line by line.
		while (getline(parameterfile, line)) {
			// Skip lines designated as section breaks and section headers.
			if ((line.substr(0, 2)).compare("--") != 0 && (line.substr(0, 2)).compare("##") != 0) {
				// Strip off trailing comments from each line.
				var = line.substr(0, line.find("//"));
				var = removeWhitespace(var);
				// Add parameter value strings to a vector.
				stringvars.push_back(var);
			}
		}
		// Check that correct number of parameters have been imported
		if ((int)stringvars.size() != 42) {
			cout << "Error! Incorrect number of parameters were loaded from the parameter file." << endl;
			return false;
		}
		bool Error_found = false;
		int i = 0;
		// Convert strings into the correct data type and assign them to their corresponding parameter variable.
		// General Parameters
		Length = atoi(stringvars[i].c_str());
		i++;
		Width = atoi(stringvars[i].c_str());
		i++;
		Height = atoi(stringvars[i].c_str());
		i++;
		//enable_z_periodic_boundary
		try {
			Enable_periodic_z = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting z-direction periodic boundary conditions!" << endl;
			Error_found = true;
		}
		i++;
		Mix_fraction = atof(stringvars[i].c_str());
		i++;
		Interaction_energy1 = atof(stringvars[i].c_str());
		i++;
		Interaction_energy2 = atof(stringvars[i].c_str());
		i++;
		MC_steps = atoi(stringvars[i].c_str());
		i++;
		//enable_smoothing
		try {
			Enable_smoothing = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting morphology smoothing options" << endl;
			Error_found = true;
		}
		i++;
		Smoothing_threshold = atof(stringvars[i].c_str());
		i++;
		//enable_rescale
		try {
			Enable_rescale = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting morphology rescale options" << endl;
			Error_found = true;
		}
		i++;
		Rescale_factor = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_shrink = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting morphology shrink options" << endl;
			Error_found = true;
		}
		i++;
		//enable_interfacial_mixing
		try {
			Enable_interfacial_mixing = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting interfacial mixing options" << endl;
			Error_found = true;
		}
		i++;
		Interface_width = atof(stringvars[i].c_str());
		i++;
		Interface_conc = atof(stringvars[i].c_str());
		i++;
		//enable_analysis_only
		try {
			Enable_analysis_only = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting analysis options" << endl;
			Error_found = true;
		}
		i++;
		//enable_correlation_calc
		try {
			Enable_correlation_calc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting correlation calculation options" << endl;
			Error_found = true;
		}
		i++;
		N_sampling_max = atoi(stringvars[i].c_str());
		i++;
		//enable_mix_frac_method
		try {
			Enable_mix_frac_method = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting correlation calculation domain size determination method options" << endl;
			Error_found = true;
		}
		i++;
		//enable_e_method
		try {
			Enable_e_method = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting correlation calculation domain size determination method  options" << endl;
			Error_found = true;
		}
		i++;
		//enable_extended_correlation_calc
		try {
			Enable_extended_correlation_calc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting extended correlation calculation options" << endl;
			Error_found = true;
		}
		i++;
		Extended_correlation_cutoff_distance = atoi(stringvars[i].c_str());
		i++;
		//enable_interfacial_distance_calc
		try {
			Enable_interfacial_distance_calc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting interfacial distance calculation options" << endl;
			Error_found = true;
		}
		i++;
		//enable_tortuosity_calc
		try {
			Enable_tortuosity_calc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting tortuosity calculation options" << endl;
			Error_found = true;
		}
		i++;
		//enable_reduced_memory_tortuosity_calc
		try {
			Enable_reduced_memory_tortuosity_calc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting reduced memory tortuosity calculation options" << endl;
			Error_found = true;
		}
		i++;
		//enable_depth_dependent_cal
		try {
			Enable_depth_dependent_calc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting areal mapping calculation options" << endl;
			Error_found = true;
		}
		i++;
		//enable_areal_maps_cal
		try {
			Enable_areal_maps_calc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting depth dependent calculation options" << endl;
			Error_found = true;
		}
		i++;
		//enable_checkerboard_start
		try {
			Enable_checkerboard_start = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting checkerboard starting condition" << endl;
			Error_found = true;
		}
		i++;
		//enable_growth_pref
		try {
			Enable_growth_pref = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting growth preference conditions" << endl;
			Error_found = true;
		}
		i++;
		Growth_direction = atoi(stringvars[i].c_str());
		i++;
		Additional_interaction = atof(stringvars[i].c_str());
		i++;
		// Export Morphology Parameters
		//Enable_export_compressed_files
		try {
			Enable_export_compressed_files = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting export options" << endl;
			Error_found = true;
		}
		i++;
		//Enable_export_cross_section
		try {
			Enable_export_cross_section = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting export cross-section options" << endl;
			Error_found = true;
		}
		i++;
		// Import Morphology Options
		//Enable_import_morphologies
		try {
			Enable_import_morphologies = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting morphology import option" << endl;
			Error_found = true;
		}
		i++;
		// Tomogram Import Options
		//Enable_import_tomogram
		try {
			Enable_import_tomogram = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting tomogram import option" << endl;
			Error_found = true;
		}
		i++;
		Tomogram_name = stringvars[i];
		i++;
		Desired_unit_size = atof(stringvars[i].c_str());
		i++;
		//try {
		//	Enable_cutoff_analysis = str2bool(stringvars[i]);
		//}
		//catch (invalid_argument& exception) {
		//	cout << exception.what() << endl;
		//	cout << "Error setting tomogram import conditions" << endl;
		//	Error_found = true;
		//}
		//i++;
		Mixed_frac = atof(stringvars[i].c_str());
		i++;
		Mixed_conc = atof(stringvars[i].c_str());
		i++;
		//try {
		//	Enable_probability_analysis = str2bool(stringvars[i]);
		//}
		//catch (invalid_argument& exception) {
		//	cout << exception.what() << endl;
		//	cout << "Error setting tomogram import conditions" << endl;
		//	Error_found = true;
		//}
		//i++;
		//Probability_scaling_exponent = atof(stringvars[i].c_str());
		//i++;
		N_extracted_segments = atoi(stringvars[i].c_str());
		i++;
		N_variants = atoi(stringvars[i].c_str());
		i++;
		return !Error_found;
	}
}
