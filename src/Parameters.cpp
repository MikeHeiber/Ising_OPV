// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "Parameters.h"

using namespace std;
using namespace Utils;

Parameters::Parameters() {

}

//  This function imports the parameters from an input parameter text file into the Input_Parameters data structure.
bool Parameters::importParameters(ifstream& parameterfile) {
	string line;
	string var;
	size_t pos;
	vector<string> stringvars;
	// Read input file line by line.
	while (parameterfile.good()) {
		getline(parameterfile, line);
		// Skip lines designated as section breaks and section headers.
		if ((line.substr(0, 2)).compare("--") != 0 && (line.substr(0, 2)).compare("##") != 0) {
			// Strip off trailing comments from each line.
			pos = line.find("/", 0);
			var = line.substr(0, pos - 1);
			// Add parameter value strings to a vector.
			stringvars.push_back(var);
		}
	}
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
		return false;
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
	//enable_export_compressed_files
	try {
		Enable_export_compressed_files = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting export options" << endl;
		return false;
	}
	i++;
	//enable_export_cross_section
	try {
		Enable_export_cross_section = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting export cross-section options" << endl;
		return false;
	}
	i++;
	//enable_smoothing
	try {
		Enable_smoothing = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting morphology smoothing options" << endl;
		return false;
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
		return false;
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
		return false;
	}
	i++;
	//enable_interfacial_mixing
	try {
		Enable_interfacial_mixing = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting interfacial mixing options" << endl;
		return false;
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
		return false;
	}
	i++;
	//enable_correlation_calc
	try {
		Enable_correlation_calc = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting correlation calculation options" << endl;
		return false;
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
		return false;
	}
	i++;
	//enable_e_method
	try {
		Enable_e_method = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting correlation calculation domain size determination method  options" << endl;
		return false;
	}
	i++;
	//enable_extended_correlation_calc
	try {
		Enable_extended_correlation_calc = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting extended correlation calculation options" << endl;
		return false;
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
		return false;
	}
	i++;
	//enable_tortuosity_calc
	try {
		Enable_tortuosity_calc = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting tortuosity calculation options" << endl;
		return false;
	}
	i++;
	//enable_reduced_memory_tortuosity_calc
	try {
		Enable_reduced_memory_tortuosity_calc = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting reduced memory tortuosity calculation options" << endl;
		return false;
	}
	i++;
	//enable_depth_dependent_cal
	try {
		Enable_depth_dependent_calc = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting areal mapping calculation options" << endl;
		return false;
	}
	i++;
	//enable_areal_maps_cal
	try {
		Enable_areal_maps_calc = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting depth dependent calculation options" << endl;
		return false;
	}
	i++;
	//enable_checkerboard_start
	try {
		Enable_checkerboard_start = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting checkerboard starting condition" << endl;
		return false;
	}
	i++;
	//enable_growth_pref
	try {
		Enable_growth_pref = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting growth preference conditions" << endl;
		return false;
	}
	i++;
	Growth_direction = atoi(stringvars[i].c_str());
	i++;
	Additional_interaction = atof(stringvars[i].c_str());
	i++;
	// Tomogram Import Options
	Desired_unit_size = atof(stringvars[i].c_str());
	i++;
	try {
		Enable_cutoff_analysis = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting tomogram import conditions" << endl;
		return false;
	}
	i++;
	Mixed_greyscale_width = atoi(stringvars[i].c_str());
	i++;
	Mixed_conc = atof(stringvars[i].c_str());
	i++;
	try {
		Enable_probability_analysis = str2bool(stringvars[i]);
	}
	catch (invalid_argument& exception) {
		cout << exception.what() << endl;
		cout << "Error setting tomogram import conditions" << endl;
		return false;
	}
	i++;
	Probability_scaling_exponent = atof(stringvars[i].c_str());
	i++;
	N_extracted_segments = atoi(stringvars[i].c_str());
	i++;
	N_variants = atoi(stringvars[i].c_str());
	i++;
	return true;
}
