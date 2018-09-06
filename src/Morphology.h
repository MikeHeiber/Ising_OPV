// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#ifndef MORPHOLOGY_H
#define MORPHOLOGY_H

#include "Lattice.h"
#include "Utils.h"
#include "tinyxml2/tinyxml2.h"
#include <algorithm>
#include <ctime>
#include <fstream>
#include <functional>
#include <numeric>
#include <sstream>

struct CorrelationCalc_Params {
	int N_sampling_max;
	bool Enable_mix_frac_method;
	bool Enable_e_method;
	bool Enable_extended_correlation_calc;
	int Correlation_cutoff_distance;
};

struct TomogramImport_Params {
	double Desired_unit_size;
	bool Enable_cutoff_analysis;
	int Mixed_greyscale_width;
	double Mixed_conc;
	bool Enable_probability_analysis;
	double Probability_scaling_exponent;
	int N_extracted_segments;
};

//! \brief This class contains a lattice representation of a materials blend with the ability to simulate phase separation and perform a variety of structural analyses.
//! \details The class makes use of the Lattice class to store the morphology data, and phase separation simulations are implemented using an Ising-based method.
//! Morphological descriptors such as thh domain size, interfacial area, tortuosity, and more can be calculated for a given morphology.
//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
//! \author Michael C. Heiber
//! \date 2014-2018
class Morphology {

	// Custom Data Structures
	struct Node {
		long int neighbor_indices[26];
		char neighbor_distances_sq[26];
		float distance_est;
		long int site_index;
	};

	struct NodeIteratorCompare {
		bool operator()(const std::vector<Node>::const_iterator& a, const std::vector<Node>::const_iterator& b) const {
			if (a == b) {
				return false;
			}
			else if (a->distance_est == b->distance_est) {
				return a < b;
			}
			else {
				return a->distance_est < b->distance_est;
			}
		}
	};

	// Data structure that stores counts of the number of neighbors that have the same site type as the main site
	// sum1 keeps track of the first-nearest neighbors
	// sum2 keeps track of the second-nearest neighbors
	// sum3 keeps track of the third-nearest neighbors
	struct NeighborCounts {
		char sum1;
		char sum2;
		char sum3;
		bool operator==(const NeighborCounts& a) const {
			return (sum1 == a.sum1 && sum2 == a.sum2);
		}
	};

	struct NeighborInfo {
		long int first_indices[6];
		long int second_indices[12];
		long int third_indices[8];
		char total1;
		char total2;
		char total3;
	};

public:
	// Functions

	//! /brief This is the simplest constructor that creates a Morphology object with the default member variables and an empty lattice.
	//! /param id is the input integer ID number that will be assigned to the Morphology object.
	Morphology(const int id);

	//! \brief This constructor creates a Morphology object with a 3D lattice with a size defined by the input dimensions (length, width, height).
	//! \details Two-dimensional periodic boundaries in the x- and y- directions are implemented by default, but periodic boundaries in the z-direction can also be enabled upon construction.
	//! Morphology objects are also tagged with an integer identification number.
	//! \param length is the x-dimension size of the lattice.
	//! \param width is the y-dimension size of the lattice.
	//! \param height is the z-dimension size of the lattice.
	//! \param enable_periodic_z is a boolean parameter that specifies whether or not periodic boundaries are enabled in the z-direction.
	//! \param id is the input integer ID number that will be assigned to the Morphology object.
	Morphology(const int length, const int width, const int height, const bool enable_periodic_z, const int id);

	//! \brief This constructor creates a Morphology object with a 3D lattice defined by the input Lattice object.
	//! \param input_lattice is the input Lattice object that will be used to create the Morphology.
	//! \param id is the input integer ID number that will be assigned to the Morphology object.
	Morphology(const Lattice& input_lattice, const int id);

	//! \brief This is the default virtual destructor.
	virtual ~Morphology();

	//! \brief Calculates the domain size anisotropy of each phase 
	//! \param N_sampling_max defines the maximum number of sites that will be sampled from the lattice when the lattice has more sites than N_sampling_max.
	void calculateAnisotropies(const int N_sampling_max);

	//! \brief Calculates the correlation length data and the domain size using the input parameter options.
	//! \param parameters is the input data structure that contains all of the parameters needed by the correlation function calculation algorithm.
	void calculateCorrelationDistances(const CorrelationCalc_Params& parameters);

	//! \brief Calculates the lattice depth dependent (z-direction) characteristics of the morphology.
	//! \details Calculates the depth dependent composition, interfacial volume fraction, and domain size.
	//! \param correlation_params is the input data structure that contains all of the parameters needed by the correlation function calculation algorithm.
	void calculateDepthDependentData(const CorrelationCalc_Params& correlation_params);

	//! \brief Calculates the interfacial area in units of lattice units squared.
	double calculateInterfacialAreaVolumeRatio() const;

	//! \brief This function calculates the interfacial distance histograms, which gives the fraction of sites at a specified distance from the interface.
	void calculateInterfacialDistance();

	//  This function calculates the number of sites that are adjacent to a site of the opposite type, which represents the interfacial volume in lattice units cubed.
	double calculateInterfacialVolumeFraction() const;

	//  This function calculates the fraction of each type sites in the lattice to the total number of sites.
	void calculateMixFractions();

	//  This function calculates the tortuosity histograms that characterize the morphology.
	//  For all type 1 sites, the shortest distance from the site along a path through other type 1 sites to the boundary at z=0 is calculated.
	//  For all type 2 sites, the shortest distance from the site along a path through other type 2 sites to the boundary at z=Height-1 is calculated.
	//  The resulting shortest path divided by the straight vertical path is the tortuosity of the pathway.
	//  The shortest paths are calculated using Dijkstra's algorithm
	bool calculateTortuosity(const char site_tye, const bool electrode_num, const bool enable_reduced_memory);

	void createCheckerboardMorphology();

	//  This function creates a randomly mixed morphology on the lattice.
	//  Sites are randomly assigned based on the mix_fractions.
	void createRandomMorphology(const std::vector<double>& mix_fractions);

	//  This function enables interactions between third-neighbor sites that are a distance of sqrt(3) lattice units apart.
	//  By default third-neighbor interactions are disabled, so this function must be called to enable this option.
	void enableThirdNeighborInteraction();

	//  This function implements num_MCsteps iterations of the Ising site swapping process.
	//  This function uses the bond formation algorithm to determine the energy change in the system that results from swapping two neighboring sites.
	//  The energy change is determined by the input parameters interaction_energy1 and interaction_energy2, which are in units of kT.
	//  These parameters describe the preference for like-like interactions over like-unlike interactions for each site type.
	//  Positive values of the interaction energies result in a driving force for phase separation.
	void executeIsingSwapping(const int num_MCsteps, const double interaction_energy1, const double interaction_energy2, const bool enable_growth_pref, const int growth_direction, const double additional_interaction); // bond formation algorithm

	//  This function implements interfacial mixing with a specified interfacial width and a specified mixing concentration in the interfacial region.
	//  Mixing is implemented by first determining the bounds on either side of the interface where mixing should occur
	//  Then random swapping of type 1 and type 2 sites within the bounds creates mixing in the interfacial region.
	void executeMixing(const double width, const double interfacial_conc);

	//  This function smoothens out rough domain interfaces and removes small islands and island sites.
	//  This is done by determining a roughness factor for each site that is given by the fraction of surrounding sites that are a different type.
	//  Sites with a roughness factor is greater than the specified smoothing_threshold are switched to the opposite type.
	//  A rescale dependent smoothing process is executed when the rescale factor is greater than 1.
	void executeSmoothing(const double smoothing_threshold, const int rescale_factor);

	//  This function returns a vector containing the pair-pair correlation function data for the specified site type.
	std::vector<double> getCorrelationData(const char site_type) const;

	std::vector<double> getDepthCompositionData(const char site_type) const;

	std::vector<double> getDepthDomainSizeData(const char site_type) const;

	std::vector<double> getDepthIVData(const char site_type) const;

	//  This function returns the domain size determined for the specified site type.
	//  This function will return zero if the calculateCorrelationDistance function has not been called.
	double getDomainSize(const char site_type) const;

	//  This function returns the domain anisotropy determined for the specified site type.
	//  This function will return zero if the calculateAnisotropy function has not been called.
	double getDomainAnisotropy(const char site_type) const;

	//  This function returns the height or z-direction size of the lattice.
	int getHeight() const;

	//  This function returns a vector containing the interfacial distance histogram data for the specified site type.
	std::vector<double> getInterfacialHistogram(const char site_type) const;

	//  This function returns the island volume for the specified site type.
	double getIslandVolumeFraction(const char site_type) const;

	//  This function returns the length or x-direction size of the lattice.
	int getLength() const;

	//  This function return the mix fraction of the morphology.
	double getMixFraction(const char site_type) const;

	//  This function returns a vector containing the end-to-end tortuosity data for the specified site type.
	std::vector<float> getTortuosityData(const char site_type) const;

	//  This function returns a vector containing the overall tortuosity histogram data for all sites with the specified site type.
	std::vector<double> getTortuosityHistogram(const char site_type) const;

	double getUnitSize() const;

	//  This function returns the width or y-direction size of the lattice.
	int getWidth() const;

	//! \brief imports a tomogram dataset using a combination of information from an xml metadata file and a set of import parameters
	//! \param info_filename is the name of the xml metadata file
	//! \param data_filename is the name of the raw morphology data file
	//! \param params is the TomogramImportParams data structure that contains the additional information needed to handle the tomogram data
	//! \returns a vector of Morphology objects that consists of a series of subsections of the original tomogram data
	std::vector<Morphology> importTomogramMorphologyFile(const std::string& info_filename, const std::string& data_filename, const TomogramImport_Params& params);

	//  This function imports the morphology text file given by the input file stream.
	//  It must be specified whether or not the input file is in the compressed format.
	bool importMorphologyFile(std::ifstream& infile);

	void outputCompositionMaps(std::ofstream& outfile) const;

	void outputCorrelationData(std::ofstream& outfile) const;

	void outputDepthDependentData(std::ofstream& outfilt) const;

	//  This function outputs the morphology data to a text file specified by the output file stream.
	//  The user can specify whether to use the compress text format or not.
	void outputMorphologyFile(std::string version, std::ofstream& outfile, const bool enable_export_compressed_files) const;

	//  This function outputs to a text file a cross-section of the morphology at x=0 plane.
	void outputMorphologyCrossSection(std::ofstream& outfile) const;

	void outputTortuosityMaps(std::ofstream& outfile) const;

	//  This function shrinks the existing lattice by a fraction of 1 over the integer value called rescale_factor.
	//  Each of the original lattice dimensions must be divisible by the rescale factor
	//  This original lattice is overwritten by the newly created smaller lattice
	void shrinkLattice(const int rescale_factor);

	//  This function stretches the existing lattice by a integer value called rescale_factor.
	//  This original lattice is overwritten by the newly created larger rescale_factor lattice
	void stretchLattice(const int rescale_factor);

protected:
private:
	// Member variables
	int ID = 0;
	std::vector<double> Mix_fractions; // Volume fraction of each component to total
	bool Enable_third_neighbor_interaction = false;
	Lattice lattice;
	std::vector<char> Site_types;
	std::vector<int> Site_type_counts;
	std::vector<std::vector<double>> Correlation_data;
	std::vector<std::vector<float>> Tortuosity_data;
	std::vector<std::vector<double>> InterfacialHistogram_data;
	std::vector<std::vector<double>> TortuosityHistogram_data;
	std::vector<std::vector<double>> Depth_composition_data;
	std::vector<std::vector<double>> Depth_iv_data;
	std::vector<std::vector<double>> Depth_domain_size_data;
	std::vector<bool> Domain_anisotropy_updated;
	std::vector<double> Domain_sizes;
	std::vector<double> Domain_anisotropies;
	std::vector<int> Island_volume;
	std::vector<long int> Interfacial_sites;
	std::vector<NeighborCounts> Neighbor_counts;
	std::vector<NeighborInfo> Neighbor_info;
	NeighborCounts Temp_counts1;
	NeighborCounts Temp_counts2;
	std::mt19937_64 gen;

	// Functions
	void addSiteType(const char site_type);

	//  This function calculates the additional change in energy of the system that would occur if the adjacent sites at (x1,y1,z1) and (x2,y2,z2) were to be swapped with
	//  a given preferential domain growth direction.  This additional energy is to be used to modify the total energy change from swapping the two sites.
	//  Sites must be adjacent to each other for calculation to be correct. (Works for adjacent sites across periodic boundaries)
	//  When non-periodic/hard z-boundaries are used, it is assumed that neither site type has a preferential interaction with the z-boundary
	//  The values for growth_direction are 1 for x-direction, 2 for y-direction, and 3 for z-direction adjustment.
	double calculateAdditionalEnergyChange(const long int site_index_main, const long int site_index_neighbor, const int growth_direction, const double additional_interaction) const;

	//  This function calculates the anisotropy of the domains based on the directionally-dependent pair-pair correlation functions
	//  The correlation function is calculated from each starting site out to the cutoff distance.
	//  The correlation length in each direction is defined as the distance at which the pair-pair correlation function first crosses the value equal to the mixing fraction
	//  If this cross-over point is not reach within the cutoff distance, the function generates an error message and returns -1.
	//  For large lattices, the correlation function does not need to be calculated starting from every site to collect enough statistics and instead a sampling of starting sites can be used.
	//  When the total number of sites is greater than N_sampling_max, N_sampling_max sites are randomly selected and saved for performing a correlation function calculation by sampling.
	//  When the total number of sites is less than N_sampling_max, all sites will be used as starting points for the correlation function calculation.
	bool calculateAnisotropy(const std::vector<long int>& correlation_sites, const char site_type, const int cutoff_distance);

	//  This function calculates the domain size of the morphology based on the pair-pair correlation function
	//  The correlation function is calculated from each starting site out to the cutoff distance.
	//  The domain size is defined as the distance at which the pair-pair correlation function first crosses the value equal to the mixing fraction
	//  If this cross-over point is not reach within the cutoff distance, the function returns false.
	//  When the extended calculation is enabled the correlation function must reach the next peak, otherwise the function returns false.
	//  For large lattices, the correlation function does not need to be calculated starting from every site to collect enough statistics and instead a sampling of starting sites can be used.
	//  When the total number of sites is greater than N_sampling_max, N_sampling_max sites are randomly selected and saved for performing a correlation function calculation by sampling.
	//  When the total number of sites is less than N_sampling_max, all sites will be used as starting points for the correlation function calculation.
	//  If the function returns false and the function is re-called with a larger cutoff_distance, the correlation function is not recalculated for close distances and only fills in the missing data for larger distances.
	double calculateCorrelationDistance(const std::vector<long int>& correlation_sites, std::vector<double>& correlation_data, const double mix_fraction, const int cutoff_distance, const CorrelationCalc_Params& params);

	//  This function calculates the fraction of nearby sites the site at (x,y,z) that are not the same type.
	//  The radius that determines which sites are included as nearby sites is determined by the rescale factor parameter.
	//  This function is designed to be used by the executeSmoothing function and implement rescale factor dependent smoothing.
	double calculateDissimilarFraction(const Coords& coords, const int rescale_factor) const;

	//  This function calculates the change in energy of the system that would occur if the adjacent sites at (x1,y1,z1) and (x2,y2,z2) were to be swapped
	//  Sites must be adjacent to each other for calculation to be correct. (Works for adjacent sites across periodic boundaries)
	//  When non-periodic/hard z-boundaries are used, it is assumed that neither site type has a preferential interaction with the z-boundary
	double calculateEnergyChangeSimple(const long int site_index1, const long int site_index2, const double interaction_energy1, const double interaction_energy2);

	//  Calculates the change in energy of the system that would occur if the adjacent sites at (x1,y1,z1) and (x2,y2,z2) were to be swapped
	//  Sites must be adjacent to each other for calculation to be correct. (Works for adjacent sites across periodic boundaries)
	//  When non-periodic/hard z-boundaries are used, it is assumed that neither site type has a preferential interaction with the z-boundary
	double calculateEnergyChange(const Coords& coords1, const Coords& coords2, const double interaction_energy1, const double interaction_energy2) const;

	Morphology::NeighborCounts calculateNeighborCounts(const Coords& coords) const;

	//  This function calculates the shortest pathways through the domains in the morphology using Dijkstra's algorithm.
	//  For all type 1 sites, the shortest distance from each site along a path through other type 1 sites to the boundary at z=0 is calculated.
	//  For all type 2 sites, the shortest distance from each site along a path through other type 2 sites to the boundary at z=Height-1 is calculated.
	bool calculatePathDistances(std::vector<float>& path_distances);

	//  This function calculates the shortest pathways through the domains in the morphology using Dijkstra's algorithm.
	//  For all type 1 sites, the shortest distance from each site along a path through other type 1 sites to the boundary at z=0 is calculated.
	//  For all type 2 sites, the shortest distance from each site along a path through other type 2 sites to the boundary at z=Height-1 is calculated.
	//  As opposed to the calculatePathDistances function, this function uses less memory at the expense of more calculation time.
	bool calculatePathDistances_ReducedMemory(std::vector<float>& path_distances);

	// This function writes the node data for the site at the given x, y, z coordinates to the specified input node variable.
	// Each node contains a vector with indices of all first- ,second-, and third-nearest neighbors (at most 26 neighbors).
	// Another vector stores the squared distance to each of the neighbors.
	// Each node also has an estimated distance from the destination and the corresponding site index.
	void createNode(Node& node, const Coords& coords);

	void getSiteSampling(std::vector<long int>& sites, const char site_type, const int N_sites);

	void getSiteSamplingZ(std::vector<long int>& sites, const char site_type, const int N_sites, const int z);

	int getSiteTypeIndex(const char site_type) const;

	//  This function initializes the neighbor_info and neighbor_counts vectors for the morphology.  The neighbor_info vector contains counts of the number of first, second, and
	//  third nearest-neighbors and three site index vectors, one for each type of neighbors, that point to each of the neighbors.  The neighbor_counts vector contains counts of the
	//  number of similar type first, second and third nearest-neighbors.
	void initializeNeighborInfo();

	//  This function determines whether the site at (x,y,z) is within the specified distance from the interface.
	//  If so, the function returns true and if not, the function returns false.
	bool isNearInterface(const Coords& coords, const double distance) const;

	double rand01();

	//  This function is called after two sites are swapped, and it updates the neighbor_counts vector, which stores the number of similar type neighbors that each site has.
	void updateNeighborCounts(const long int site_index1, const long int site_index2);
};

#endif // MORPHOLOGY_H
