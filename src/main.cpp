// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "Morphology.h"
#include "Utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <mpi.h>

using namespace std;
using namespace Utils;

struct Input_Params {
	// General
	int Length; // x-direction size of the lattice
	int Width; // y-direction size of the lattice
	int Height; // z-direction size of the lattice
	bool Enable_periodic_z; // z-direction periodic boundary option
	double Mix_fraction; // volume fraction of donor
	double Interaction_energy1; // energetic favorability for type1-type1 interactions over type1-type2 interactions
	double Interaction_energy2; // energetic favorability for type2-type2 interactions over type1-type2 interactions
	int MC_steps; // number of MC steps to be executed (determines number of Ising swapping iterations)
	// Export Morphology Options
	bool Enable_export_compressed_files; // choose whether the output morphology data file is in compressed format or not
	bool Enable_export_cross_section; // choose whether to output data for an uncompressed cross section (x=0 plane) of the morphology
	// Smoothing Options
	bool Enable_smoothing; // choose whether to perform domain smoothing or not
	double Smoothing_threshold;  // specify the degree of smoothing
	// Rescale Options
	bool Enable_rescale; // choose whether to perform lattice rescaling or not
	int Rescale_factor; // specify the rescale factor to be used (must be an integer greater than 1)
	bool Enable_shrink; // chose whether to shrink the lattice instead of expand it
	// Interfacial Mixing Options
	bool Enable_interfacial_mixing; // choose whether to perform interfacial mixing or not
	double Interface_width;  // specify the interfacial width
	double Interface_conc; // specify the mixing concentration in the interfacial region
	// Analysis Options
	bool Enable_analysis_only;  // choose whether to only perform analysis on an imported morphology or not
	bool Enable_correlation_calc; // choose whether to perform the domain size calculation using the pair-pair correlation method or not
	int N_sampling_max; // specify the maximum number of sites to be sampled for calculating the correlation function
	bool Enable_mix_frac_method; // choose whether to use the mix fraction method for determining the domain size or not
	bool Enable_e_method; // choose whether to use the 1/e method for determining the domain size or not
	bool Enable_extended_correlation_calc; // choose whether to extend the correlation function calculation to the specified cutoff distance or not
	int Correlation_cutoff_distance; // specify the maximum cutoff distnace for the extended correlation function calculation
	bool Enable_interfacial_distance_calc; // choose whether to calculate the interfacial distance histograms or not
	bool Enable_tortuosity_calc; // choose whether to calculate the tortuosity histograms, end-to-end tortuosity, and island volume fraction or not
	bool Enable_reduced_memory_tortuosity_calc; // choose whether to enable a tortuosity calculation method that takes longer, but uses less memory or not
	bool Enable_depth_dependent_calc; // choose whether to enable calculation and output of film depth dependent morphology characteristics or not
	bool Enable_areal_maps_calc; // choose whether to enable calculation and output of areal mappings of morphology characteristics or not
	// Other Options
	bool Enable_checkerboard_start; //choose whether to start from a alternating checkerboard-like configuration instead of a random blend (creates 0.5 mix fraction) or not
	bool Enable_growth_pref;
	int Growth_direction;
	double Additional_interaction;
	// Tomogram Import Options
	double Desired_unit_size;
	bool Enable_cutoff_analysis;
	int Mixed_greyscale_width;
	double Mixed_conc;
	bool Enable_probability_analysis;
	double Probability_scaling_exponent;
	int N_extracted_segments;
	int N_variants;
};

bool importParameters(ifstream& parameterfile, Input_Params& params, CorrelationCalc_Params& correlation_params);

int main(int argc, char * argv[]) {
	// Input parameters
	Input_Params parameters;
	CorrelationCalc_Params correlation_params;
	// Internal parameters
	string version = "v4.0-beta.2";
	bool Enable_import_morphology = false;
	bool Enable_import_tomogram = false;
	double mix_ratio = 0;
	double domain_size1 = 0;
	double domain_size2 = 0;
	double domain_anisotropy1 = 0;
	double domain_anisotropy2 = 0;
	double iav_ratio = 0;
	double iv_fraction = 0;
	double island_fraction1 = 0;
	double island_fraction2 = 0;
	int N_steps = 0;
	double elapsedtime = 0;
	time_t start_time, end_time;
	int procid = 0;
	int nproc = 1;
	string filename;
	ifstream parameter_file;
	ifstream morphology_input_file;
	ofstream analysis_file;
	ofstream areal_composition_file;
	ofstream areal_tortuosity_file;
	ofstream correlation_avg_file;
	ofstream correlation_file;
	ofstream depthdata_avg_file;
	ofstream depthdata_file;
	ofstream interfacial_dist_hist_file;
	ofstream morphology_output_file;
	ofstream morphology_cross_section_file;
	ofstream tortuosity_hist_file;
	ofstream path_data1_file;
	ofstream path_data2_file;
	bool success;
	double *mix_ratios = NULL;
	double *domain_sizes1 = NULL;
	double *domain_sizes2 = NULL;
	double *domain_anisotropies1 = NULL;
	double *domain_anisotropies2 = NULL;
	double *iav_ratios = NULL;
	double *iv_fractions = NULL;
	double *island_fractions1 = NULL;
	double *island_fractions2 = NULL;
	double *times = NULL;
	int pathdata1_size = 0;
	int pathdata2_size = 0;
	int pathdata1_count = 0;
	int pathdata2_count = 0;
	float *pathdata1 = NULL;
	float *pathdata2 = NULL;
	float *pathdata1_all = NULL;
	float *pathdata2_all = NULL;
	int *pathdata1_sizes = NULL;
	int *pathdata2_sizes = NULL;
	int *pathdata1_displacement = NULL;
	int *pathdata2_displacement = NULL;
	vector<double> tortuosity_hist1_vect;
	vector<double> tortuosity_hist2_vect;
	vector<double> interfacial_dist_hist1_vect;
	vector<double> interfacial_dist_hist2_vect;
	vector<double> correlation1_vect;
	vector<double> correlation2_vect;
	vector<double> depth_comp1_vect;
	vector<double> depth_comp2_vect;
	vector<double> depth_iv1_vect;
	vector<double> depth_iv2_vect;
	vector<double> depth_size1_vect;
	vector<double> depth_size2_vect;
	string input_morphology, input_file_path, filename_prefix;
	string tomo_info_filename, tomo_data_filename;
	// Begin
	start_time = time(NULL);
	// Initialize parallel processing.
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	// Import parameters from text file.
	if (argc > 3) {
		string argument = string(argv[2]);
		if (argument.compare("-import") == 0) {
			cout << procid << ": Attempting to import Ising_OPV morphology data..." << endl;
			Enable_import_morphology = true;
		}
		else if (argument.compare("-importTomogram") == 0) {
			cout << procid << ": Attempting to import tomogram data..." << endl;
			Enable_import_tomogram = true;
			Enable_import_morphology = true;
		}
	}
	else if (argc == 2) {
		// Standard operation
	}
	else {
		cout << procid << ": Incorrect arguments. There needs to be at least one argument which is the filename of a parameter file. Program will exit now!" << endl;
		return 0;
	}
	parameter_file.open((argv[1]), ifstream::in);
	success = importParameters(parameter_file, parameters, correlation_params);
	parameter_file.close();
	if (success) {
		cout << procid << ": Parameter file found and loaded successfully." << endl;
	}
	else {
		cout << procid << ": Error importing variables from file!  Program will exit now." << endl;
		return 0;
	}
	// Parameter checks
	if (parameters.Enable_analysis_only && !Enable_import_morphology && !Enable_import_tomogram) {
		cout << procid << ": Error!  The 'analysis only' option can only be used when a morphology is imported from a file." << endl;
		return 0;
	}
	if (parameters.Enable_analysis_only) {
		cout << procid << ": Warning! Only morphology analysis will be performed." << endl;
	}
	if (Enable_import_tomogram && ((parameters.Enable_cutoff_analysis && parameters.Enable_probability_analysis) || (!parameters.Enable_cutoff_analysis && !parameters.Enable_probability_analysis))) {
		cout << procid << " Error! When importing a tomogram the cutoff analysis or the probability analysis option must be enabled, but not both." << endl;
		return 0;
	}
	// Wait until all processors have loaded the parameters.
	MPI_Barrier(MPI_COMM_WORLD);
	// Create morphology data structure.
	Morphology morph(parameters.Length, parameters.Width, parameters.Height, parameters.Enable_periodic_z, procid);
	// Import morphology if enabled.
	// Import tomogram file in binary format
	if (Enable_import_tomogram) {
		if (procid == 0) {
			if (parameters.N_extracted_segments*parameters.N_variants != nproc) {
				cout << procid << ": Error! The number of processors used must be equal to N_extracted_segments*N_variants. ";
				cout << parameters.N_extracted_segments*parameters.N_variants << " processors are needed but " << nproc << " were requested." << endl;
				return 0;
			}
			// Create initial lattice based on filename to load the tomogram data into
			Morphology morph_tomo(procid);
			// Collect tomogram import options
			TomogramImport_Params import_params;
			import_params.Desired_unit_size = parameters.Desired_unit_size;
			import_params.Enable_cutoff_analysis = parameters.Enable_cutoff_analysis;
			import_params.Mixed_greyscale_width = parameters.Mixed_greyscale_width;
			import_params.Mixed_conc = parameters.Mixed_conc;
			import_params.Enable_probability_analysis = parameters.Enable_probability_analysis;
			import_params.Probability_scaling_exponent = parameters.Probability_scaling_exponent;
			import_params.N_extracted_segments = parameters.N_extracted_segments;
			cout << procid << ": Loading and analyzing tomogram data." << endl;
			tomo_info_filename = string(argv[3]);
			tomo_data_filename = string(argv[4]);
			vector<Morphology> morphology_set = morph_tomo.importTomogramMorphologyFile(tomo_info_filename, tomo_data_filename, import_params);
			// Check that a set of morphologies has been produced
			if (morphology_set.size() == 0) {
				cout << procid << ": Error! Morphology set could not be generated from the input tomogram." << endl;
				return 0;
			}
			for (int i = 1; i < parameters.N_variants; i++) {
				vector<Morphology> morphology_set2 = morph_tomo.importTomogramMorphologyFile(tomo_info_filename, tomo_data_filename, import_params);
				morphology_set.insert(morphology_set.end(), morphology_set2.begin(), morphology_set2.end());
			}
			// Output morphology set to separate files
			for (int i = 0; i < (int)morphology_set.size(); i++) {
				cout << procid << ": Writing morphology " << i << " to an output file." << endl;
				filename = "morphology_" + to_string(i) + ".txt";
				morphology_output_file.open(filename);
				morphology_set[i].outputMorphologyFile(version, morphology_output_file, parameters.Enable_export_compressed_files);
				morphology_output_file.close();
			}
		}
		// All processors must wait until the root proc finishes with morphology set generation.
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (Enable_import_morphology) {
		size_t pos_path, pos_id;
		// Filename from tomogram set
		if (Enable_import_tomogram) {
			input_morphology = "morphology_#.txt";
		}
		// Get filename of imported morphology from command line arguments.
		else {
			input_morphology = (argv[3]);
		}
		// Separate filepath and differentiate between Windows (1) and Linux (2) filepaths.
		// Linux filepaths
		if (input_morphology.substr(0, 1).compare("/") == 0) {
			pos_path = input_morphology.find_last_of("/");
			input_file_path = input_morphology.substr(0, pos_path) + "/";
			pos_path++;
		}
		// Windows filepaths
		else if (input_morphology.substr(0, 2).compare("\\") == 0) {
			pos_path = input_morphology.find_last_of("\\");
			input_file_path = input_morphology.substr(0, pos_path) + "\\";
			pos_path += 2;
		}
		// Default to be used when the morphology data file is in the current working directory
		else {
			input_file_path = "";
			pos_path = 0;
		}
		// Separate filename prefix from ID number.
		pos_id = input_morphology.find_last_of("_");
		filename_prefix = input_morphology.substr(pos_path, pos_id);
		filename = filename_prefix + "_" + to_string(procid) + ".txt";
		cout << procid << ": Opening morphology file " << input_file_path << filename << endl;
		morphology_input_file.open(filename);
		if (morphology_input_file.is_open()) {
			cout << procid << ": Morphology file successfully opened!" << endl;
		}
		else {
			cout << procid << ": Opening morphology file failed! Program will exit now!" << endl;
			return 0;
		}

		cout << procid << ": Importing morphology from file..." << flush;
		// Import the morphology from the given data file.
		morph.importMorphologyFile(morphology_input_file);
		morphology_input_file.close();
		cout << procid << ": Morphology import complete!" << endl;
	}
	// Create new morphology if import is disabled.
	else {
		if (parameters.Enable_checkerboard_start) {
			cout << procid << ": Generating initial checkerboard morphology..." << endl;
			morph.createCheckerboardMorphology();
		}
		else {
			cout << procid << ": Generating initial random morphology..." << endl;
			vector<double> mix_vec(2, 0);
			mix_vec[0] = parameters.Mix_fraction;
			mix_vec[1] = 1 - parameters.Mix_fraction;
			morph.createRandomMorphology(mix_vec);
		}
	}
	// Determine if any phase separation is to be executed on the morphology.
	if (parameters.MC_steps > 0 && !parameters.Enable_analysis_only) {
		N_steps = parameters.MC_steps;
	}
	else {
		N_steps = 0;
	}
	// Execute phase separation through Ising swapping.
	if (N_steps > 0) {
		cout << procid << ": Executing site swapping for " << N_steps << " MC steps..." << endl;
		morph.executeIsingSwapping(N_steps, parameters.Interaction_energy1, parameters.Interaction_energy2, parameters.Enable_growth_pref, parameters.Growth_direction, parameters.Additional_interaction);
	}
	// Perform lattice rescaling and domain smoothing if enabled.
	if (parameters.Enable_rescale && !parameters.Enable_analysis_only) {
		if (parameters.Enable_shrink) {
			cout << procid << ": Initial blend ratio is " << morph.getMixFraction((char)1) << endl;
			if (parameters.Enable_smoothing) {
				cout << procid << ": Executing standard smoothing with a smoothing threshold of " << parameters.Smoothing_threshold << "..." << endl;
				morph.executeSmoothing(parameters.Smoothing_threshold, 1);
				cout << procid << ": Blend ratio after smoothing is " << morph.getMixFraction((char)1) << endl;
			}
			cout << procid << ": Shrinking lattice by a factor of " << parameters.Rescale_factor << " ..." << endl;
			morph.shrinkLattice(parameters.Rescale_factor);
			parameters.Length = parameters.Length / parameters.Rescale_factor;
			parameters.Width = parameters.Width / parameters.Rescale_factor;
			parameters.Height = parameters.Height / parameters.Rescale_factor;
			cout << procid << ": Blend ratio after shrinking lattice is " << morph.getMixFraction((char)1) << endl;
		}
		else {
			cout << procid << ": Expanding lattice by a factor of " << parameters.Rescale_factor << " ..." << endl;
			morph.stretchLattice(parameters.Rescale_factor);
			parameters.Length = parameters.Length*parameters.Rescale_factor;
			parameters.Width = parameters.Width*parameters.Rescale_factor;
			parameters.Height = parameters.Height*parameters.Rescale_factor;
		}
	}
	if (parameters.Enable_smoothing && !parameters.Enable_shrink && !parameters.Enable_analysis_only) {
		if (!parameters.Enable_rescale) {
			cout << procid << ": Executing standard smoothing with a smoothing threshold of " << parameters.Smoothing_threshold << "..." << endl;
			morph.executeSmoothing(parameters.Smoothing_threshold, 1);
		}
		else {
			cout << procid << ": Executing rescale factor dependent smoothing with a smoothing threshold of " << parameters.Smoothing_threshold << "..." << endl;
			morph.executeSmoothing(parameters.Smoothing_threshold, parameters.Rescale_factor);
		}
	}
	// Perform interfacial mixing if enabled.
	if (parameters.Enable_interfacial_mixing && !parameters.Enable_analysis_only) {
		cout << "Executing interfacial mixing..." << endl;
		morph.executeMixing(parameters.Interface_width, parameters.Interface_conc);
	}
	// Calculate domain size if enabled.
	if (parameters.Enable_correlation_calc) {
		morph.calculateCorrelationDistances(correlation_params);
		filename = "correlation_data_" + to_string(procid) + ".txt";
		correlation_file.open(filename);
		morph.outputCorrelationData(correlation_file);
		correlation_file.close();
		domain_size1 = morph.getDomainSize((char)1);
		domain_size2 = morph.getDomainSize((char)2);
		morph.calculateAnisotropies(parameters.N_sampling_max);
		domain_anisotropy1 = morph.getDomainAnisotropy((char)1);
		domain_anisotropy2 = morph.getDomainAnisotropy((char)2);
	}
	// Calculate interfacial distance histogram if enabled.
	if (parameters.Enable_interfacial_distance_calc) {
		cout << procid << ": Calculating the interfacial distance histogram..." << endl;
		morph.calculateInterfacialDistanceHistogram();
	}
	// Calculate interfacial area to volume ratio.
	iav_ratio = morph.calculateInterfacialAreaVolumeRatio();
	// Calculate interfacial volume to total volume ratio.
	iv_fraction = morph.calculateInterfacialVolumeFraction();
	// Get Final Mix ratio
	mix_ratio = morph.getMixFraction((char)1);
	// Calculate end-to-end tortuosity, tortuosity histogram, and island volume fraction.
	if (parameters.Enable_tortuosity_calc) {
		if (parameters.Enable_reduced_memory_tortuosity_calc) {
			cout << procid << ": Calculating tortuosity using the reduced memory method..." << endl;
		}
		else {
			cout << procid << ": Calculating tortuosity using the standard method..." << endl;
		}
		success = morph.calculateTortuosity((char)1, parameters.Enable_reduced_memory_tortuosity_calc);
		success = morph.calculateTortuosity((char)2, parameters.Enable_reduced_memory_tortuosity_calc);
		if (!success) {
			cout << procid << ": Error calculating tortuosity! Program will exit now." << endl;
			return 0;
		}
		if (parameters.Enable_areal_maps_calc) {
			cout << procid << " Creating areal tortuosity map." << endl;
			filename = "areal_tortuosity_map_" + to_string(procid) + ".txt";
			areal_tortuosity_file.open(filename);
			morph.outputTortuosityMaps(areal_tortuosity_file);
			areal_tortuosity_file.close();
		}
		// Calculate island volume ratio.
		island_fraction1 = (double)morph.getIslandVolumeFraction((char)1);
		island_fraction2 = (double)morph.getIslandVolumeFraction((char)2);
	}
	if (parameters.Enable_depth_dependent_calc) {
		cout << procid << ": Calculating the depth dependent composition and domain size..." << endl;
		morph.calculateDepthDependentData(correlation_params);
		filename = "depth_dependent_data_" + to_string(procid) + ".txt";
		depthdata_file.open(filename);
		morph.outputDepthDependentData(depthdata_file);
		depthdata_file.close();
	}
	if (parameters.Enable_areal_maps_calc) {
		cout << procid << " Creating areal composition map." << endl;
		filename = "areal_composition_map_" + to_string(procid) + ".txt";
		areal_composition_file.open(filename);
		morph.outputCompositionMaps(areal_composition_file);
		areal_composition_file.close();
	}
	// Save final morphology to a text file.
	if (!parameters.Enable_analysis_only || Enable_import_tomogram) {
		cout << procid << ": Writing morphology to file..." << endl;
		if (!Enable_import_morphology || Enable_import_tomogram) {
			filename = "morphology_" + to_string(procid) + ".txt";
		}
		else {
			filename = filename_prefix + "mod_" + to_string(procid) + ".txt";
		}
		morphology_output_file.open(filename);
		morph.outputMorphologyFile(version, morphology_output_file, parameters.Enable_export_compressed_files);
		morphology_output_file.close();
	}
	// Save the cross-section of the x=0 plane to a file if enabled.
	if (parameters.Enable_export_cross_section) {
		filename = "morphology_" + to_string(procid) + "_cross_section.txt";
		morphology_cross_section_file.open(filename);
		morph.outputMorphologyCrossSection(morphology_cross_section_file);
		morphology_cross_section_file.close();
	}
	// Morphology generation is now finished.
	end_time = time(NULL);
	elapsedtime = (double)difftime(end_time, start_time) / 60;
	// All processors must finish with morphology generation before analysis can begin.
	MPI_Barrier(MPI_COMM_WORLD);
	// Gather morphology data from each processor and calculate set statistics.
	if (procid == 0) {
		cout << "Collecting morphology analysis data from each processor..." << endl;
	}
	if (parameters.Enable_tortuosity_calc) {
		// Get path data for end-to-end tortuosity calculation.
		vector<float> data = morph.getTortuosityData((char)1);
		// Determine size of the data vector.
		pathdata1_size = (int)data.size();
		// Allocate an array to store the data.
		pathdata1 = (float *)malloc(sizeof(float)*pathdata1_size);
		// Put data from the vector into the array.
		for (int i = 0; i < pathdata1_size; i++) {
			pathdata1[i] = data[i];
		}
		data.clear();
		// Repeat for type 2 path data.
		data = morph.getTortuosityData((char)2);
		pathdata2_size = (int)data.size();
		pathdata2 = (float *)malloc(sizeof(float)*pathdata2_size);
		for (int i = 0; i < pathdata2_size; i++) {
			pathdata2[i] = data[i];
		}
		data.clear();
		// Create array on the root processor that will contain the size of the path data arrays from each processor.
		if (procid == 0) {
			pathdata1_sizes = (int *)malloc(sizeof(int)*nproc);
			pathdata2_sizes = (int *)malloc(sizeof(int)*nproc);
		}
		// Gather the size of the data arrays from all processors into the size array on the root processor.
		MPI_Gather(&pathdata1_size, 1, MPI_INT, pathdata1_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(&pathdata2_size, 1, MPI_INT, pathdata2_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// Calculate the average tortuosity histograms.
		tortuosity_hist1_vect = MPI_calculateVectorAvg(morph.getTortuosityHistogram((char)1));
		tortuosity_hist2_vect = MPI_calculateVectorAvg(morph.getTortuosityHistogram((char)2));
	}
	// Calculate the average interfacial distance histograms.
	if (parameters.Enable_interfacial_distance_calc) {
		interfacial_dist_hist1_vect = MPI_calculateVectorAvg(morph.getInterfacialDistanceHistogram((char)1));
		interfacial_dist_hist2_vect = MPI_calculateVectorAvg(morph.getInterfacialDistanceHistogram((char)2));
	}
	// Calculate the average pair-pair correlation functions.
	if (parameters.Enable_correlation_calc) {
		correlation1_vect = MPI_calculateVectorAvg(morph.getCorrelationData((char)1));
		correlation2_vect = MPI_calculateVectorAvg(morph.getCorrelationData((char)2));
	}
	// Calculate the average depth dependent characteristics
	if (parameters.Enable_depth_dependent_calc) {
		depth_comp1_vect = MPI_calculateVectorAvg(morph.getDepthCompositionData((char)1));
		depth_comp2_vect = MPI_calculateVectorAvg(morph.getDepthCompositionData((char)2));
		depth_iv1_vect = MPI_calculateVectorAvg(morph.getDepthIVData((char)1));
		depth_iv2_vect = MPI_calculateVectorAvg(morph.getDepthIVData((char)2));
		depth_size1_vect = MPI_calculateVectorAvg(morph.getDepthDomainSizeData((char)1));
		depth_size2_vect = MPI_calculateVectorAvg(morph.getDepthDomainSizeData((char)2));
	}
	// Prepare root processor for data collection.
	if (procid == 0) {
		// Create arrays to store property values gathered from each processor.
		mix_ratios = (double *)malloc(sizeof(double)*nproc);
		iav_ratios = (double *)malloc(sizeof(double)*nproc);
		iv_fractions = (double *)malloc(sizeof(double)*nproc);
		times = (double *)malloc(sizeof(double)*nproc);
		if (parameters.Enable_correlation_calc) {
			domain_sizes1 = (double *)malloc(sizeof(double)*nproc);
			domain_sizes2 = (double *)malloc(sizeof(double)*nproc);
			domain_anisotropies1 = (double *)malloc(sizeof(double)*nproc);
			domain_anisotropies2 = (double *)malloc(sizeof(double)*nproc);
		}
		if (parameters.Enable_tortuosity_calc) {
			island_fractions1 = (double *)malloc(sizeof(double)*nproc);
			island_fractions2 = (double *)malloc(sizeof(double)*nproc);
			// The final path data array will have a size that is the sum of the sizes of the arrays coming from each processor.
			for (int i = 0; i < nproc; i++) {
				pathdata1_count += pathdata1_sizes[i];
				pathdata2_count += pathdata2_sizes[i];
			}
			pathdata1_all = (float *)malloc(sizeof(float)*pathdata1_count);
			pathdata2_all = (float *)malloc(sizeof(float)*pathdata2_count);
			// The path data array from the processors may have different sizes, and so displacements must be determined to know where to insert them into the full path data array.
			pathdata1_displacement = (int *)malloc(sizeof(int)*nproc);
			pathdata2_displacement = (int *)malloc(sizeof(int)*nproc);
			pathdata1_displacement[0] = 0;
			pathdata2_displacement[0] = 0;
			for (int i = 1; i < nproc; i++) {
				pathdata1_displacement[i] = pathdata1_displacement[i - 1] + pathdata1_sizes[i - 1];
				pathdata2_displacement[i] = pathdata2_displacement[i - 1] + pathdata2_sizes[i - 1];
			}
		}
	}
	// Gather the properties from each processor into the previously created arrays on the root processor.
	MPI_Gather(&mix_ratio, 1, MPI_DOUBLE, mix_ratios, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&iav_ratio, 1, MPI_DOUBLE, iav_ratios, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&iv_fraction, 1, MPI_DOUBLE, iv_fractions, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&elapsedtime, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (parameters.Enable_correlation_calc) {
		MPI_Gather(&domain_size1, 1, MPI_DOUBLE, domain_sizes1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&domain_size2, 1, MPI_DOUBLE, domain_sizes2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&domain_anisotropy1, 1, MPI_DOUBLE, domain_anisotropies1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&domain_anisotropy2, 1, MPI_DOUBLE, domain_anisotropies2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	// Gather the path data arrays from each processor into the previously created full path data array on the root processor at the positions defined by the displacement array.
	if (parameters.Enable_tortuosity_calc) {
		MPI_Gatherv(pathdata1, pathdata1_size, MPI_FLOAT, pathdata1_all, pathdata1_sizes, pathdata1_displacement, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Gatherv(pathdata2, pathdata2_size, MPI_FLOAT, pathdata2_all, pathdata2_sizes, pathdata2_displacement, MPI_FLOAT, 0, MPI_COMM_WORLD);
		// Gather the island volume fraction property from each processor into the previously created array on the root processor.
		MPI_Gather(&island_fraction1, 1, MPI_DOUBLE, island_fractions1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gather(&island_fraction2, 1, MPI_DOUBLE, island_fractions2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	// Output the analysis results to text files.
	if (procid == 0) {
		cout << "Writing morphology analysis data to files..." << endl;
		// Output the average pair-pair correlation function.
		if (parameters.Enable_correlation_calc) {
			if (!Enable_import_morphology || Enable_import_tomogram) {
				correlation_avg_file.open("correlation_data_avg.txt");
			}
			else {
				if (parameters.Enable_analysis_only) {
					correlation_avg_file.open("correlation_data_avg_new.txt");
				}
				else {
					correlation_avg_file.open("correlation_data_avg_mod.txt");
				}
			}
			correlation_avg_file << "Distance (nm),Correlation1,Correlation2" << endl;
			for (int i = 0; i < (int)correlation1_vect.size(); i++) {
				correlation_avg_file << morph.getUnitSize()*(double)i*0.5 << "," << correlation1_vect[i] << "," << correlation2_vect[i] << endl;
			}
			correlation_avg_file.close();
		}
		// Output the average tortuosity histograms and the end-to-end path data.
		if (parameters.Enable_tortuosity_calc) {
			if (!Enable_import_morphology || Enable_import_tomogram) {
				tortuosity_hist_file.open("tortuosity_histograms.txt");
				path_data1_file.open("end-to-end_path_data1.txt");
				path_data2_file.open("end-to-end_path_data2.txt");
			}
			else {
				if (parameters.Enable_analysis_only) {
					tortuosity_hist_file.open("tortuosity_histograms_new.txt");
					path_data1_file.open("end-to-end_path_data1_new.txt");
					path_data2_file.open("end-to-end_path_data2_new.txt");
				}
				else {
					tortuosity_hist_file.open("tortuosity_histograms_mod.txt");
					path_data1_file.open("end-to-end_path_data1_mod.txt");
					path_data2_file.open("end-to-end_path_data2_mod.txt");
				}
			}
			tortuosity_hist_file << "Distance (nm),Tortuosity1,Tortuosity2" << endl;
			tortuosity_hist_file << 0 << "," << tortuosity_hist1_vect[0] << "," << tortuosity_hist2_vect[0] << endl;
			int hist_size = (tortuosity_hist1_vect.size() > tortuosity_hist2_vect.size()) ? (int)tortuosity_hist1_vect.size() : (int)tortuosity_hist2_vect.size();
			for (int i = 1; i < hist_size; i++) {
				tortuosity_hist_file << morph.getUnitSize()*(i + 49.0) / 50.0 << ",";
				if (i < tortuosity_hist1_vect.size()) {
					tortuosity_hist_file << tortuosity_hist1_vect[i] << ",";
				}
				else {
					tortuosity_hist_file << "0,";
				}
				if (i < tortuosity_hist2_vect.size()) {
					tortuosity_hist_file << tortuosity_hist2_vect[i] << ",";
				}
				else {
					tortuosity_hist_file << "0,";
				}
			}
			for (int i = 0; i < pathdata1_count; i++) {
				path_data1_file << pathdata1_all[i] << "\n";
			}
			for (int i = 0; i < pathdata2_count; i++) {
				path_data2_file << pathdata2_all[i] << "\n";
			}
			tortuosity_hist_file.close();
			path_data1_file.close();
			path_data2_file.close();
		}
		// Output the interfacial distance histograms.
		if (parameters.Enable_interfacial_distance_calc) {
			if (!Enable_import_morphology || Enable_import_tomogram) {
				interfacial_dist_hist_file.open("interfacial_distance_histograms.txt");
			}
			else {
				if (parameters.Enable_analysis_only) {
					interfacial_dist_hist_file.open("interfacial_distance_histograms_new.txt");
				}
				else {
					interfacial_dist_hist_file.open("interfacial_distance_histograms_mod.txt");
				}
			}
			interfacial_dist_hist_file << "Distance (a),Probability1,Probability2" << endl;
			for (int i = 0; i < (int)interfacial_dist_hist1_vect.size(); i++) {
				interfacial_dist_hist_file << i + 1 << "," << interfacial_dist_hist1_vect[i] << "," << interfacial_dist_hist2_vect[i] << "\n";
			}
			interfacial_dist_hist_file.close();
		}
		// Output the average depth dependent data.
		if (parameters.Enable_depth_dependent_calc) {
			if (!Enable_import_morphology || Enable_import_tomogram) {
				depthdata_avg_file.open("depth_dependent_data_avg.txt");
			}
			else {
				if (parameters.Enable_analysis_only) {
					depthdata_avg_file.open("depth_dependent_data_avg_new.txt");
				}
				else {
					depthdata_avg_file.open("depth_dependent_data_avg_mod.txt");
				}
			}
			depthdata_avg_file << "Z-Position,Type1_composition,Type2_composition,Type1_IV_fraction,Type2_IV_fraction,Type1_domain_size,Type2_domain_size" << endl;
			for (int i = 0; i < (int)depth_size1_vect.size(); i++) {
				depthdata_avg_file << i << "," << depth_comp1_vect[i] << "," << depth_comp2_vect[i] << "," << depth_iv1_vect[i] << "," << depth_iv2_vect[i] << "," << depth_size1_vect[i] << "," << depth_size2_vect[i] << endl;
			}
			depthdata_avg_file.close();
		}
		// Output the final morphology set analysis summary to a text file.
		if (!Enable_import_morphology || Enable_import_tomogram) {
			analysis_file.open("analysis_summary.txt");
		}
		else {
			if (parameters.Enable_analysis_only) {
				analysis_file.open("analysis_summary_new.txt");
			}
			else {
				analysis_file.open("analysis_summary_mod.txt");
			}
		}
		analysis_file << "Summary of results for this morphology set containing " << nproc << " morphologies created using Ising_OPV " << version << ":" << endl;
		analysis_file << "length,width,height,mix_ratio_avg,mix_ratio_stdev,domain1_size_avg,domain1_size_stdev,domain2_size_avg,domain2_size_stdev,";
		analysis_file << "domain1_anisotropy_avg,domain1_anisotropy_stdev,domain2_anisotropy_avg,domain2_anisotropy_stdev,";
		analysis_file << "interfacial_area_volume_ratio_avg,interfacial_area_volume_ratio_stdev,interfacial_volume_ratio_avg,interfacial_volume_ratio_stdev,";
		analysis_file << "tortuosity1_avg,tortuosity1_stdev,tortuosity2_avg,tortuosity2_stdev,island_volume_ratio1_avg,island_volume_ratio1_stdev,";
		analysis_file << "island_volume_ratio2_avg,island_volume_ratio2_stdev,calc_time_avg(min),calc_time_stdev(min)" << endl;
		analysis_file << morph.getLength() << "," << morph.getWidth() << "," << morph.getHeight() << ",";
		analysis_file << array_avg(mix_ratios, nproc) << "," << array_stdev(mix_ratios, nproc) << ",";
		if (parameters.Enable_correlation_calc) {
			analysis_file << array_avg(domain_sizes1, nproc) << "," << array_stdev(domain_sizes1, nproc) << "," << array_avg(domain_sizes2, nproc) << "," << array_stdev(domain_sizes2, nproc) << ",";
			analysis_file << array_avg(domain_anisotropies1, nproc) << "," << array_stdev(domain_anisotropies1, nproc) << "," << array_avg(domain_anisotropies2, nproc) << "," << array_stdev(domain_anisotropies2, nproc) << ",";
		}
		else {
			analysis_file << "-" << "," << "-" << "," << "-" << "," << "-" << ",";
			analysis_file << "-" << "," << "-" << "," << "-" << "," << "-" << ",";
		}
		analysis_file << array_avg(iav_ratios, nproc) << "," << array_stdev(iav_ratios, nproc) << "," << array_avg(iv_fractions, nproc) << "," << array_stdev(iv_fractions, nproc) << ",";
		if (parameters.Enable_tortuosity_calc) {
			analysis_file << array_avg(pathdata1_all, pathdata1_count) << "," << array_stdev(pathdata1_all, pathdata1_count) << ",";
			analysis_file << array_avg(pathdata2_all, pathdata2_count) << "," << array_stdev(pathdata2_all, pathdata2_count) << ",";
			analysis_file << array_avg(island_fractions1, nproc) << "," << array_stdev(island_fractions1, nproc) << ",";
			analysis_file << array_avg(island_fractions2, nproc) << "," << array_stdev(island_fractions2, nproc) << ",";
		}
		else {
			analysis_file << "-" << "," << "-" << ",";
			analysis_file << "-" << "," << "-" << ",";
			analysis_file << "-" << "," << "-" << ",";
			analysis_file << "-" << "," << "-" << ",";
		}
		analysis_file << array_avg(times, nproc) << "," << array_stdev(times, nproc) << endl;
		analysis_file << endl;
		analysis_file << "Detailed results for each of the morphologies in the set:" << endl;
		analysis_file << "id#,length,width,height,mix_ratio,domain1_size,domain2_size,";
		analysis_file << "domain1_anisotropy,domain2_anisotropy,";
		analysis_file << "interfacial_area_volume_ratio,interfacial_volume_ratio,";
		analysis_file << "tortuosity1,tortuosity2,island_volume_ratio1,island_volume_ratio2,calc_time(min)" << endl;
		double *tortuosities1 = (double *)malloc(sizeof(double)*nproc);
		double *tortuosities2 = (double *)malloc(sizeof(double)*nproc);
		for (int i = 0; i < nproc; i++) {
			analysis_file << i << "," << morph.getLength() << "," << morph.getWidth() << "," << morph.getHeight() << "," << mix_ratios[i] << ",";
			if (parameters.Enable_correlation_calc) {
				analysis_file << domain_sizes1[i] << "," << domain_sizes2[i] << ",";
				analysis_file << domain_anisotropies1[i] << "," << domain_anisotropies2[i] << ",";
			}
			else {
				analysis_file << "-" << "," << "-" << ",";
				analysis_file << "-" << "," << "-" << ",";
			}
			analysis_file << iav_ratios[i] << "," << iv_fractions[i] << ",";
			if (parameters.Enable_tortuosity_calc) {
				pathdata1 = (float *)malloc(sizeof(float)*pathdata1_sizes[i]);
				pathdata2 = (float *)malloc(sizeof(float)*pathdata2_sizes[i]);
				copy(pathdata1_all + pathdata1_displacement[i], pathdata1_all + pathdata1_displacement[i] + pathdata1_sizes[i], pathdata1);
				copy(pathdata2_all + pathdata2_displacement[i], pathdata2_all + pathdata2_displacement[i] + pathdata2_sizes[i], pathdata2);
				tortuosities1[i] = array_avg(pathdata1, pathdata1_sizes[i]);
				tortuosities2[i] = array_avg(pathdata2, pathdata2_sizes[i]);
				if (std::isnan(tortuosities1[i])) {
					tortuosities1[i] = -1;
				}
				if (std::isnan(tortuosities2[i])) {
					tortuosities2[i] = -1;
				}
				analysis_file << tortuosities1[i] << "," << tortuosities2[i] << ",";
				analysis_file << island_fractions1[i] << "," << island_fractions2[i] << ",";
			}
			else {
				analysis_file << "-" << "," << "-" << ",";
				analysis_file << "-" << "," << "-" << ",";
			}
			analysis_file << times[i] << endl;
		}
		analysis_file << endl;
		if (parameters.Enable_correlation_calc) {
			analysis_file << "Morphology number " << array_which_median(domain_sizes1, nproc) << " has the median domain1 size of " << array_median(domain_sizes1, nproc) << endl;
		}
		if (parameters.Enable_tortuosity_calc) {
			analysis_file << "Morphology number " << array_which_median(tortuosities1, nproc) << " has the median tortuosity1 of " << array_median(tortuosities1, nproc) << endl;
		}
		if (Enable_import_tomogram) {
			analysis_file << endl;
			analysis_file << "Morphologies imported from tomogram file: " << tomo_data_filename << endl;
		}
		analysis_file.close();
		cout << "Finished!" << endl;
	}
	MPI_Finalize();
	// The morphology object is automatically deconstructed upon return.
	return 0;
}

//  This function imports the parameters from an input parameter text file into the Input_Parameters data structure.
bool importParameters(ifstream& parameterfile, Input_Params& params, CorrelationCalc_Params& correlation_params) {
	string line;
	string var;
	size_t pos;
	bool error_status = false;
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
	params.Length = atoi(stringvars[i].c_str());
	i++;
	params.Width = atoi(stringvars[i].c_str());
	i++;
	params.Height = atoi(stringvars[i].c_str());
	i++;
	//enable_z_periodic_boundary
	params.Enable_periodic_z = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting z-direction periodic boundary conditions!" << endl;
		return false;
	}
	i++;
	params.Mix_fraction = atof(stringvars[i].c_str());
	i++;
	params.Interaction_energy1 = atof(stringvars[i].c_str());
	i++;
	params.Interaction_energy2 = atof(stringvars[i].c_str());
	i++;
	params.MC_steps = atoi(stringvars[i].c_str());
	i++;
	//enable_export_compressed_files
	params.Enable_export_compressed_files = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting export options" << endl;
		return false;
	}
	i++;
	//enable_export_cross_section
	params.Enable_export_cross_section = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting export cross-section options" << endl;
		return false;
	}
	i++;
	//enable_smoothing
	params.Enable_smoothing = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting morphology smoothing options" << endl;
		return false;
	}
	i++;
	params.Smoothing_threshold = atof(stringvars[i].c_str());
	i++;
	//enable_rescale
	params.Enable_rescale = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting morphology rescale options" << endl;
		return false;
	}
	i++;
	params.Rescale_factor = atoi(stringvars[i].c_str());
	i++;
	params.Enable_shrink = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting morphology shrink options" << endl;
		return false;
	}
	i++;
	//enable_interfacial_mixing
	params.Enable_interfacial_mixing = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting interfacial mixing options" << endl;
		return false;
	}
	i++;
	params.Interface_width = atof(stringvars[i].c_str());
	i++;
	params.Interface_conc = atof(stringvars[i].c_str());
	i++;
	//enable_analysis_only
	params.Enable_analysis_only = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting analysis options" << endl;
		return false;
	}
	i++;
	//enable_correlation_calc
	params.Enable_correlation_calc = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting correlation calculation options" << endl;
		return false;
	}
	i++;
	params.N_sampling_max = atoi(stringvars[i].c_str());
	i++;
	//enable_mix_frac_method
	params.Enable_mix_frac_method = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting correlation calculation domain size determination method options" << endl;
		return false;
	}
	i++;
	//enable_e_method
	params.Enable_e_method = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting correlation calculation domain size determination method  options" << endl;
		return false;
	}
	i++;
	//enable_extended_correlation_calc
	params.Enable_extended_correlation_calc = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting extended correlation calculation options" << endl;
		return false;
	}
	i++;
	params.Correlation_cutoff_distance = atoi(stringvars[i].c_str());
	i++;
	//enable_interfacial_distance_calc
	params.Enable_interfacial_distance_calc = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting interfacial distance calculation options" << endl;
		return false;
	}
	i++;
	//enable_tortuosity_calc
	params.Enable_tortuosity_calc = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting tortuosity calculation options" << endl;
		return false;
	}
	i++;
	//enable_reduced_memory_tortuosity_calc
	params.Enable_reduced_memory_tortuosity_calc = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting reduced memory tortuosity calculation options" << endl;
		return false;
	}
	i++;
	//enable_depth_dependent_cal
	params.Enable_depth_dependent_calc = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting areal mapping calculation options" << endl;
		return false;
	}
	i++;
	//enable_areal_maps_cal
	params.Enable_areal_maps_calc = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting depth dependent calculation options" << endl;
		return false;
	}
	i++;
	//enable_checkerboard_start
	params.Enable_checkerboard_start = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting checkerboard starting condition" << endl;
		return false;
	}
	i++;
	//enable_growth_pref
	params.Enable_growth_pref = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting growth preference conditions" << endl;
		return false;
	}
	i++;
	params.Growth_direction = atoi(stringvars[i].c_str());
	i++;
	params.Additional_interaction = atof(stringvars[i].c_str());
	i++;
	// Tomogram Import Options
	params.Desired_unit_size = atof(stringvars[i].c_str());
	i++;
	params.Enable_cutoff_analysis = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting tomogram import conditions" << endl;
		return false;
	}
	i++;
	params.Mixed_greyscale_width = atoi(stringvars[i].c_str());
	i++;
	params.Mixed_conc = atof(stringvars[i].c_str());
	i++;
	params.Enable_probability_analysis = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting tomogram import conditions" << endl;
		return false;
	}
	i++;
	params.Probability_scaling_exponent = atof(stringvars[i].c_str());
	i++;
	params.N_extracted_segments = atoi(stringvars[i].c_str());
	i++;
	params.N_variants = atoi(stringvars[i].c_str());
	i++;
	correlation_params.N_sampling_max = params.N_sampling_max;
	correlation_params.Enable_mix_frac_method = params.Enable_mix_frac_method;
	correlation_params.Enable_e_method = params.Enable_e_method;
	correlation_params.Enable_extended_correlation_calc = params.Enable_extended_correlation_calc;
	correlation_params.Correlation_cutoff_distance = params.Correlation_cutoff_distance;
	return true;
}
