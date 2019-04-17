// Copyright (c) 2014-2019 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "Morphology.h"
#include "Parameters.h"
#include "Utils.h"
#include "Version.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <mpi.h>

using namespace Ising_OPV;
using namespace std;

int main(int argc, char * argv[]) {
	// Input parameters
	Parameters parameters;
	// Internal parameters
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
	bool success;
	vector<double> mix_ratios;
	vector<double> domain_sizes1;
	vector<double> domain_sizes2;
	vector<double> domain_anisotropies1;
	vector<double> domain_anisotropies2;
	vector<double> iav_ratios;
	vector<double> iv_fractions;
	vector<double> island_fractions1;
	vector<double> island_fractions2;
	vector<double> times;
	vector<double> tortuosity_data1;
	vector<double> tortuosity_data2;
	vector<double> tortuosity_hist1_vect;
	vector<double> tortuosity_hist2_vect;
	vector<pair<double, double>> interfacial_dist_probhist1;
	vector<pair<double, double>> interfacial_dist_probhist2;
	vector<double> correlation1_vect;
	vector<double> correlation2_vect;
	vector<double> depth_comp1_vect;
	vector<double> depth_comp2_vect;
	vector<double> depth_iv_vect;
	vector<double> depth_size1_vect;
	vector<double> depth_size2_vect;
	// Begin
	start_time = time(NULL);
	// Initialize parallel processing.
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	// Import parameters from text file.
	if (argc != 2) {
		cout << procid << ": Incorrect arguments. There needs to be at least one argument which is the filename of a parameter file. Program will exit now!" << endl;
		return 0;
	}
	parameter_file.open(argv[1]);
	success = parameters.importParameters(parameter_file);
	parameter_file.close();
	if (success) {
		cout << procid << ": Parameter file found and loaded successfully." << endl;
	}
	else {
		cout << procid << ": Error importing variables from file!  Program will exit now." << endl;
		return 0;
	}
	// Check validity of loaded parameters
	success = parameters.checkParameters();
	if (!success) {
		cout << procid << ": Error! One or more imported parameters are invalid.  Program will exit now." << endl;
		return 0;
	}
	// Wait until all processors have loaded the parameters.
	MPI_Barrier(MPI_COMM_WORLD);
	// Create morphology data structure.
	Morphology morph(parameters, procid);
	// Import morphology if enabled.
	// Import tomogram file in binary format
	if (parameters.Enable_import_tomogram) {
		if (procid == 0) {
			if (parameters.N_extracted_segments*parameters.N_variants != nproc) {
				cout << ": Error! The number of processors used must be equal to N_extracted_segments*N_variants. Program will exit now.";
				cout << parameters.N_extracted_segments * parameters.N_variants << " processors are needed but " << nproc << " were requested." << endl;
				return 0;
			}
			// Collect tomogram import options
			cout << procid << ": Loading and analyzing tomogram data." << endl;
			vector<Morphology> morphology_set = morph.importTomogramMorphologyFile();
			// Check that a set of morphologies has been produced
			if (morphology_set.size() == 0) {
				cout << procid << ": Error! Morphology set could not be generated from the input tomogram. Program will exit now." << endl;
				return 0;
			}
			for (int i = 1; i < parameters.N_variants; i++) {
				vector<Morphology> morphology_set2 = morph.importTomogramMorphologyFile();
				morphology_set.insert(morphology_set.end(), morphology_set2.begin(), morphology_set2.end());
			}
			// Output morphology set to separate files
			for (int i = 0; i < (int)morphology_set.size(); i++) {
				cout << procid << ": Writing morphology " << i << " to an output file." << endl;
				filename = "morphology_" + to_string(i) + ".txt";
				morphology_output_file.open(filename);
				morphology_set[i].outputMorphologyFile(morphology_output_file, parameters.Enable_export_compressed_files);
				morphology_output_file.close();
			}
		}
		// All processors must wait until the root proc finishes with morphology set generation.
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (parameters.Enable_import_morphologies || parameters.Enable_import_tomogram) {
		filename = "morphology_" + to_string(procid) + ".txt";
		cout << procid << ": Opening morphology file " << filename << endl;
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
		success = morph.importMorphologyFile(morphology_input_file);
		if (!success) {
			cout << procid << ": Importing morphology file failed! Program will exit now!" << endl;
			return 0;
		}
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
		morph.calculateCorrelationDistances();
		filename = "correlation_data_" + to_string(procid) + ".txt";
		correlation_file.open(filename);
		morph.outputCorrelationData(correlation_file);
		correlation_file.close();
		domain_size1 = morph.getDomainSize((char)1);
		domain_size2 = morph.getDomainSize((char)2);
		morph.calculateAnisotropies();
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
		morph.calculateDepthDependentData();
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
	if (!parameters.Enable_analysis_only || parameters.Enable_import_tomogram) {
		cout << procid << ": Writing morphology to file..." << endl;
		filename = "morphology_" + to_string(procid) + ".txt";
		morphology_output_file.open(filename);
		morph.outputMorphologyFile(morphology_output_file, parameters.Enable_export_compressed_files);
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
		// Collect end-to-end tortuosity data from all procs onto proc 0
		tortuosity_data1 = MPI_gatherVectors(morph.getTortuosityData((char)1));
		tortuosity_data2 = MPI_gatherVectors(morph.getTortuosityData((char)2));
	}
	// Calculate the average interfacial distance histograms.
	if (parameters.Enable_interfacial_distance_calc) {
		auto hist1 = morph.getInterfacialDistanceHistogram((char)1);
		auto hist2 = morph.getInterfacialDistanceHistogram((char)2);
		interfacial_dist_probhist1 = MPI_calculateProbHistAvg(hist1);
		interfacial_dist_probhist2 = MPI_calculateProbHistAvg(hist2);
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
		depth_iv_vect = MPI_calculateVectorAvg(morph.getDepthIVData());
		depth_size1_vect = MPI_calculateVectorAvg(morph.getDepthDomainSizeData((char)1));
		depth_size2_vect = MPI_calculateVectorAvg(morph.getDepthDomainSizeData((char)2));
	}
	// Gather the properties from each processor into the previously created arrays on the root processor.
	mix_ratios = MPI_gatherValues(mix_ratio);
	iav_ratios = MPI_gatherValues(iav_ratio);
	iv_fractions = MPI_gatherValues(iv_fraction);
	times = MPI_gatherValues(elapsedtime);
	if (parameters.Enable_correlation_calc) {
		domain_sizes1 = MPI_gatherValues(domain_size1);
		domain_sizes2 = MPI_gatherValues(domain_size2);
		domain_anisotropies1 = MPI_gatherValues(domain_anisotropy1);
		domain_anisotropies2 = MPI_gatherValues(domain_anisotropy2);
	}
	if (parameters.Enable_tortuosity_calc) {
		// Gather the island volume fraction property from each processor into the previously created array on the root processor.
		island_fractions1 = MPI_gatherValues(island_fraction1);
		island_fractions2 = MPI_gatherValues(island_fraction2);
	}
	// Output the analysis results to text files.
	if (procid == 0) {
		cout << "Writing morphology analysis data to files..." << endl;
		// Output the average pair-pair correlation function.
		if (parameters.Enable_correlation_calc) {
			correlation_avg_file.open("correlation_data_avg.txt");
			correlation_avg_file << "Distance (nm),Correlation1,Correlation2" << endl;
			for (int i = 0; i < (int)correlation1_vect.size(); i++) {
				correlation_avg_file << morph.getUnitSize()*(double)i*0.5 << "," << correlation1_vect[i] << "," << correlation2_vect[i] << endl;
			}
			correlation_avg_file.close();
		}
		// Output the average tortuosity histograms and the end-to-end path data.
		if (parameters.Enable_tortuosity_calc) {
			tortuosity_hist_file.open("tortuosity_histograms.txt");
			double bin_size = 0.01;
			auto probhist1 = calculateProbabilityHist(tortuosity_data1, 1.0, bin_size);
			auto probhist2 = calculateProbabilityHist(tortuosity_data2, 1.0, bin_size);
			outputVectorToFile(probhist1, "tortuosity_data");
			tortuosity_hist_file << "Tortuosity,Probability1,Probability2" << endl;
			int hist_size = (probhist1.size() > probhist2.size()) ? (int)probhist1.size() : (int)probhist2.size();
			for (int i = 0; i < hist_size; i++) {
				// output bin value
				tortuosity_hist_file << bin_size * i + 1.0 + bin_size / 2.0 << ",";
				if (i < (int)probhist1.size()) {
					tortuosity_hist_file << probhist1[i].second << ",";
				}
				else {
					tortuosity_hist_file << "0,";
				}
				if (i < (int)probhist2.size()) {
					tortuosity_hist_file << probhist2[i].second << endl;;
				}
				else {
					tortuosity_hist_file << "0" << endl;
				}
			}
			tortuosity_hist_file.close();
		}
		// Output the interfacial distance histograms.
		if (parameters.Enable_interfacial_distance_calc) {
			interfacial_dist_hist_file.open("interfacial_distance_histograms.txt");
			interfacial_dist_hist_file << "Distance (a),Probability1,Probability2" << endl;
			int hist_size = (interfacial_dist_probhist1.size() > interfacial_dist_probhist2.size()) ? (int)interfacial_dist_probhist1.size() : (int)interfacial_dist_probhist2.size();
			for (int i = 0; i < hist_size; i++) {
				interfacial_dist_hist_file << i + 1 << ",";
				if (i < (int)interfacial_dist_probhist1.size()) {
					interfacial_dist_hist_file << interfacial_dist_probhist1[i].second << ",";
				}
				else {
					interfacial_dist_hist_file << "0,";
				}
				if (i < (int)interfacial_dist_probhist2.size()) {
					interfacial_dist_hist_file << interfacial_dist_probhist2[i].second << endl;
				}
				else {
					interfacial_dist_hist_file << "0" << endl;
				}
			}
			interfacial_dist_hist_file.close();
		}
		// Output the average depth dependent data.
		if (parameters.Enable_depth_dependent_calc) {
			depthdata_avg_file.open("depth_dependent_data_avg.txt");
			depthdata_avg_file << "Z-Position,Type1_composition,Type2_composition,Type1_domain_size,Type2_domain_size,IV_fraction" << endl;
			for (int i = 0; i < (int)depth_size1_vect.size(); i++) {
				depthdata_avg_file << i << "," << depth_comp1_vect[i] << "," << depth_comp2_vect[i] << "," << depth_size1_vect[i] << "," << depth_size2_vect[i] << "," << depth_iv_vect[i] << endl;
			}
			depthdata_avg_file.close();
		}
		// Output the final morphology set analysis summary to a text file.
		analysis_file.open("analysis_summary.txt");
		analysis_file << "Summary of results for this morphology set containing " << nproc << " morphologies created using Ising_OPV v" << Current_version.getVersionStr() << ":" << endl;
		analysis_file << "length,width,height,mix_ratio_avg,mix_ratio_stdev,domain1_size_avg,domain1_size_stdev,domain2_size_avg,domain2_size_stdev,";
		analysis_file << "domain1_anisotropy_avg,domain1_anisotropy_stdev,domain2_anisotropy_avg,domain2_anisotropy_stdev,";
		analysis_file << "interfacial_area_volume_ratio_avg,interfacial_area_volume_ratio_stdev,interfacial_volume_ratio_avg,interfacial_volume_ratio_stdev,";
		analysis_file << "tortuosity1_avg,tortuosity1_stdev,tortuosity2_avg,tortuosity2_stdev,island_volume_ratio1_avg,island_volume_ratio1_stdev,";
		analysis_file << "island_volume_ratio2_avg,island_volume_ratio2_stdev,calc_time_avg(min),calc_time_stdev(min)" << endl;
		analysis_file << morph.getLength() << "," << morph.getWidth() << "," << morph.getHeight() << ",";
		analysis_file << vector_avg(mix_ratios) << "," << vector_stdev(mix_ratios) << ",";
		if (parameters.Enable_correlation_calc) {
			analysis_file << vector_avg(domain_sizes1) << "," << vector_stdev(domain_sizes1) << "," << vector_avg(domain_sizes2) << "," << vector_stdev(domain_sizes2) << ",";
			analysis_file << vector_avg(domain_anisotropies1) << "," << vector_stdev(domain_anisotropies1) << "," << vector_avg(domain_anisotropies2) << "," << vector_stdev(domain_anisotropies2) << ",";
		}
		else {
			analysis_file << "-" << "," << "-" << "," << "-" << "," << "-" << ",";
			analysis_file << "-" << "," << "-" << "," << "-" << "," << "-" << ",";
		}
		analysis_file << vector_avg(iav_ratios) << "," << vector_stdev(iav_ratios) << "," << vector_avg(iv_fractions) << "," << vector_stdev(iv_fractions) << ",";
		if (parameters.Enable_tortuosity_calc) {
			analysis_file << vector_avg(tortuosity_data1) << "," << vector_stdev(tortuosity_data1) << ",";
			analysis_file << vector_avg(tortuosity_data2) << "," << vector_stdev(tortuosity_data2) << ",";
			analysis_file << vector_avg(island_fractions1) << "," << vector_stdev(island_fractions1) << ",";
			analysis_file << vector_avg(island_fractions2) << "," << vector_stdev(island_fractions2) << ",";
		}
		else {
			analysis_file << "-" << "," << "-" << ",";
			analysis_file << "-" << "," << "-" << ",";
			analysis_file << "-" << "," << "-" << ",";
			analysis_file << "-" << "," << "-" << ",";
		}
		analysis_file << vector_avg(times) << "," << vector_stdev(times) << endl;
		analysis_file << endl;
		analysis_file << "Detailed results for each of the morphologies in the set:" << endl;
		analysis_file << "id#,length,width,height,mix_ratio,domain1_size,domain2_size,";
		analysis_file << "domain1_anisotropy,domain2_anisotropy,";
		analysis_file << "interfacial_area_volume_ratio,interfacial_volume_ratio,";
		analysis_file << "tortuosity1,tortuosity2,island_volume_ratio1,island_volume_ratio2,calc_time(min)" << endl;
		vector<double> tortuosities1(nproc);
		vector<double> tortuosities2(nproc);
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
				tortuosities1[i] = vector_avg(morph.getTortuosityData((char)1));
				tortuosities2[i] = vector_avg(morph.getTortuosityData((char)2));
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
			int index = vector_which_median(domain_sizes1);
			analysis_file << "Morphology number " << index << " has the median domain1 size of " << domain_sizes1[index] << endl;
		}
		if (parameters.Enable_tortuosity_calc) {
			int index = vector_which_median(tortuosities1);
			analysis_file << "Morphology number " << index << " has the median tortuosity1 of " << tortuosities1[index] << endl;
		}
		if (parameters.Enable_import_tomogram) {
			analysis_file << endl;
			analysis_file << "Morphologies imported from tomogram dataset: " << parameters.Tomogram_name << endl;
		}
		analysis_file.close();
		cout << "Finished!" << endl;
	}
	MPI_Finalize();
	// The morphology object is automatically deconstructed upon return.
	return 0;
}
