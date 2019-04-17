// Copyright (c) 2014-2019 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#ifndef MORPHOLOGY_H
#define MORPHOLOGY_H

#include "Lattice.h"
#include "Parameters.h"
#include "Utils.h"
#include "Version.h"
#include "tinyxml2/tinyxml2.h"
#include <algorithm>
#include <array>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <functional>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace Ising_OPV {

	//! \brief This class contains a lattice representation of a materials blend with the ability to simulate phase separation and perform a variety of structural analyses.
	//! \details The class makes use of the Lattice class to store the morphology data, and phase separation simulations are implemented using an Ising-based method.
	//! Morphological descriptors such as thh domain size, interfacial area, tortuosity, and more can be calculated for a given morphology.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2014-2019
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
			std::array<long int, 6> first_indices;
			std::array<long int, 12> second_indices;
			std::array<long int, 8> third_indices;
			char total1;
			char total2;
			char total3;
		};

	public:
		// Functions

		//! \brief Default constructor that creates a Morphology object with the default parameters and an empty lattice.
		//! \warning An empty Morphology object should not be used without defining a lattice.
		Morphology();

		//! \brief This constructor creates a Morphology object with a 3D lattice with a size defined by the input dimensions (length, width, height).
		//! \details Two-dimensional periodic boundaries in the x- and y- directions are implemented by default, but periodic boundaries in the z-direction can also be enabled upon construction.
		//! Morphology objects are also tagged with an integer identification number.
		//! \param params is the Parameters object that contains all parameters needed by the Morphology class.
		//! \param id is the input integer ID number that will be assigned to the Morphology object.
		Morphology(const Parameters& params, const int id);

		//! \brief This constructor creates a Morphology object with a 3D lattice defined by the input Lattice object.
		//! \param input_lattice is the input Lattice object that will be used to create the Morphology.
		//! \param params is the Parameters object that contains all parameters needed by the Morphology class.
		//! \param id is the input integer ID number that will be assigned to the Morphology object.
		Morphology(const Lattice& input_lattice, const Parameters& params, const int id);

		//! \brief This is the default virtual destructor.
		virtual ~Morphology();

		//! \brief Calculates the domain size anisotropy of each phase 
		void calculateAnisotropies();

		//! \brief Calculates the correlation length data and the domain size using the input parameter options.
		void calculateCorrelationDistances();

		//! \brief Calculates the lattice depth dependent (z-direction) characteristics of the morphology.
		//! \details Calculates the depth dependent composition, interfacial volume fraction, and domain size.
		void calculateDepthDependentData();

		//! \brief Calculates the interfacial area in units of lattice units squared.
		double calculateInterfacialAreaVolumeRatio() const;

		//! \brief Calculates the interfacial distance histograms, which gives the fraction of sites at a specified distance from the interface.
		void calculateInterfacialDistanceHistogram();

		//! \brief Calculates the fraction of sites adjacent to an interface.
		double calculateInterfacialVolumeFraction() const;

		//! \brief Calculates the volume fraction of each type site in the lattice to the total number of sites.
		void calculateMixFractions();

		//! \brief Calculates the tortuosity histogram for the specified site type.
		//! \details For all type 1 sites, the shortest paths through other type 1 sites to the boundary at z=0 is calculated.
		//! For all type 2 sites, the shortest pathes through other type 2 sites to the boundary at z=Height-1 is calculated.
		//! The shortest paths are calculated using Dijkstra's algorithm.
		//! \param site_type specifies which site type to perform the tortuosity calculation on.
		//! \param enable_reduced_memory allows users to choose to use a slower algorithm that uses less RAM.
		bool calculateTortuosity(const char site_type, const bool enable_reduced_memory);

		//! \brief Creates a split bilayer morphology in the z-direction.
		void createBilayerMorphology();

		//! \brief Creates a 3D checkerboard morphology.
		void createCheckerboardMorphology();

		//! \brief Creates a randomly mixed morphology with the specified blend ratios.
		//! \param mix_fractions is a vector that specifies the blend ratio of each site type.
		void createRandomMorphology(const std::vector<double>& mix_fractions);

		// //! \brief Enables interactions between third-neighbor sites that are a distance of sqrt(3) lattice units apart.
		// //! \details By default third-neighbor interactions are disabled, so this function must be called to enable this option.
		//void enableThirdNeighborInteraction();

		//! \brief Executes the Ising site swapping processes with the specified parameters for a given mumber of iterations.
		//! \details This function uses the bond formation algorithm to determine the energy change in the system that results from swapping two neighboring sites.
		//! Positive values for the interaction energies result in a driving force for phase separation.
		//! \param num_MCsteps is the number of Monte Carlo steps to execute, where the number of steps is defined by the total number of iterations divided by the total lattice volume.
		//! \param interaction_energy1 defines the energetic difference between like-like and unlike-unlike interactions for type 1 sites in units of kT.
		//! \param interaction_energy2 defines the energetic difference between like-like and unlike-unlike interactions for type 2 sites in units of kT.
		//! \param enable_growth_pref is a boolean option that allows users to enable preferential interations in one of the pricipal lattice directions.
		//! \param growth_direction is an integer used when directional interactions are enabled and specifies the direction with a modified interaction energy, with 1 = x-direction, 2 = y-direction, and 3 = z-direction.
		//! \param additional_interaction is used when directional interactions are enabled and specifies the additional interaction energy with sites in the specified direction.
		void executeIsingSwapping(const int num_MCsteps, const double interaction_energy1, const double interaction_energy2, const bool enable_growth_pref, const int growth_direction, const double additional_interaction);

		//! \brief Executes interfacial mixing with a specified interfacial width and interfacial mixing concentration.
		//! \details Random swapping of type 1 and type 2 sites near the interface is done to produce an mixed interfaction region with a controlled width and blend ratio.
		//! \param interfacial_width specifies the desired width of the mixed interfacial region.
		//! \param interfacial_conc specified the desired blend ratio of mixed interfacial region.
		void executeMixing(const double interfacial_width, const double interfacial_conc);

		//! \brief Executes a smoothing algorithm that smooths out rough domain interfaces and removes small islands and island sites.
		//! \details Smoothing is done by determining a roughness factor for each site that is given by the fraction of surrounding sites that are a different type.
		//! Sites with a roughness factor is greater than the specified smoothing_threshold are switched to the opposite type.
		//! A rescale dependent smoothing process is executed when the rescale factor is greater than 2.
		//! \param smoothing_threshold is the numerical parameter used to adjust how aggressive the smoothing should be.
		//! \param rescale_factor specifies whether the smoothing algorithm should be adjust to account for prior lattice rescaling by giving the rescaling factor used. 
		void executeSmoothing(const double smoothing_threshold, const int rescale_factor);

		//! \brief Returns a vector containing the pair-pair autocorrelation function data for the specified site type.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return a copy of the data vector.
		std::vector<double> getCorrelationData(const char site_type) const;

		//! \brief Returns a vector containing the film depth dependent blend composition data in the z-direction.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return a copy of the data vector.
		std::vector<double> getDepthCompositionData(const char site_type) const;

		//! \brief Returns a vector containing the film depth dependent domain size data in the z-direction.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return a copy of the data vector.
		std::vector<double> getDepthDomainSizeData(const char site_type) const;

		//! \brief Returns a vector containing the film depth dependent interfacial volume fraction data in the z-direction.
		//! \return a copy of the data vector.
		std::vector<double> getDepthIVData() const;

		//! \brief Returns the domain size determined for the specified site type.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return -1 if the calculateCorrelationDistance function has not been called.
		//! \return the domain size calculated from the correlation function data.
		double getDomainSize(const char site_type) const;

		//! \brief Returns the domain anisotropy determined for the specified site type.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return zero if the calculateAnisotropy function has not been called.
		//! \return the domain anisotropy calculated from the direction dependent correlation function data.
		double getDomainAnisotropy(const char site_type) const;

		//! \brief Gets the height or z-direction size of the lattice.
		//! \return an integer representing the height or z-direction size of the lattice.
		int getHeight() const;

		//! \brief Gets the ID number of the Morphology object.
		//! \return the integer ID number.
		int getID() const;

		//! \brief Returns a vector containing the interfacial distance probability histogram data for the specified site type.
		//! \details The histogram bins are distances given by the vector index + 1 in lattice units.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return a copy of the data vector.
		std::vector<std::pair<double, int>> getInterfacialDistanceHistogram(const char site_type) const;

		//! \brief Returns the island volume for the specified site type.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return -1 if the calculateTortuosity function has not been called.
		//! \return the volume fraction of island sites not connected to the top or bottom surface (z-dimension).
		double getIslandVolumeFraction(const char site_type) const;

		//! \brief Returns the length or x-direction size of the lattice.
		//! \return an integer representing the length or x-direction size of the lattice.
		int getLength() const;

		//! \brief Returns the mix fraction for the specified site type.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return -1 if the calculateMixFractions function has not been called.
		//! \return the volumetric blend ratio of the specified site type.
		double getMixFraction(const char site_type) const;

		//! \brief Returns a vector containing the end-to-end tortuosity data for the specified site type.
		//! \param site_type specifies for which site type the data should be retrieved.
		//! \return a copy of the data vector.
		std::vector<double> getTortuosityData(const char site_type) const;

		//! \brief Returns the lattice unit size in units of nm.
		double getUnitSize() const;

		//! \brief Returns the width or y-direction size of the lattice.
		//! \return an integer representing the width or y-direction size of the lattice.
		int getWidth() const;

		//! \brief Imports the tomogram dataset specified in the parameter file.
		//! \returns a vector of Morphology objects that consists of a series of subsections of the original tomogram data.
		std::vector<Morphology> importTomogramMorphologyFile();

		//! \brief Imports the Ising_OPV morphology text file given by the specified input filestream.
		//! \param infile is the already open input filestream pointing to an Ising_OPV morphology file.
		//! \return false if there is an error during file import.
		//! \return true if morphology file import is sucessful.
		bool importMorphologyFile(std::ifstream& infile);

		//! \brief Outputs the areal composition map data to the specified output filestream.
		//! \param outfile is the already open output filestream.
		void outputCompositionMaps(std::ofstream& outfile) const;

		//! \brief Outputs the domain autocorrelation function data to the specified output filestream.
		//! \param outfile is the already open output filestream.
		void outputCorrelationData(std::ofstream& outfile) const;

		//! \brief Outputs the lattice depth dependent (z-direction) characteristics to the specified output filestream.
		//! \details Outputs the depth dependent composition, interfacial volume fraction, and domain size.
		//! \param outfile is the already open output filestream.
		void outputDepthDependentData(std::ofstream& outfile) const;

		//! \brief Outputs the morphology data to a file specified by the output filestream.
		//! \details The user can specify whether to use the compressed format or not.
		//! \param outfile is the already open output filestream.
		//! \param enable_export_compressed is a boolean option that allows users to choose whether the output morphology file uses the compressed format or not.
		void outputMorphologyFile(std::ofstream& outfile, const bool enable_export_compressed) const;

		//! \brief Outputs a cross-section of the morphology at the x=Length/2 plane to the specified output filestream.
		//! \param outfile is the already open output filestream.
		void outputMorphologyCrossSection(std::ofstream& outfile) const;

		//! \brief Outputs the areal end-to-end tortuosity map data to the specified output filestream.
		//! \param outfile is the already open output filestream.
		void outputTortuosityMaps(std::ofstream& outfile) const;

		//! \brief Sets the parameters of the Morphology class using the input Parameters object.
		//! \param params is the Parameters object that contains all parameters needed by the Morphology class.
		void setParameters(const Parameters& params);

		//! \brief Shrinks the existing lattice by a fraction of 1 over the integer rescale_factor value.
		//! \details Each of the original lattice dimensions must be divisible by the rescale factor.
		//! The original lattice is overwritten by the newly created smaller lattice.
		void shrinkLattice(const int rescale_factor);

		//! \brief Stretches the existing lattice by a integer rescale_factor value.
		//! \details The original lattice is overwritten by the newly created larger lattice.
		void stretchLattice(const int rescale_factor);

	protected:
	private:
		// Member variables
		int ID = 0;
		Parameters Params;
		Lattice lattice;
		std::vector<double> Mix_fractions; // Volume fraction of each component to total
		std::vector<char> Site_types;
		std::vector<int> Site_type_counts;
		std::vector<std::vector<double>> Correlation_data;
		std::vector<std::vector<double>> Tortuosity_data;
		std::vector<std::vector<std::pair<double, int>>> InterfacialHistogram_data;
		std::vector<std::vector<double>> Depth_composition_data;
		std::vector<std::vector<double>> Depth_domain_size_data;
		std::vector<double> Depth_iv_data;
		std::vector<bool> Domain_anisotropy_updated;
		std::vector<double> Domain_sizes;
		std::vector<double> Domain_anisotropies;
		std::vector<int> Island_volume;
		std::vector<long int> Interfacial_sites;
		std::vector<NeighborCounts> Neighbor_counts;
		std::vector<NeighborInfo> Neighbor_info;
		NeighborCounts Temp_counts1;
		NeighborCounts Temp_counts2;
		std::mt19937_64 gen = std::mt19937_64((int)time(0));

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
		double calculateCorrelationDistance(const std::vector<long int>& correlation_sites, std::vector<double>& correlation_data, const double mix_fraction, const int cutoff_distance);

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
		//double calculateEnergyChange(const Coords& coords1, const Coords& coords2, const double interaction_energy1, const double interaction_energy2) const;

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
}

#endif // MORPHOLOGY_H
