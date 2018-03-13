// Copyright (c) 2017 Michael C. Heiber
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

// Data structure that stores counts of the number of neighbors that have the same site type as the main site
// sum1 keeps track of the first-nearest neighbors
// sum2 keeps track of the second-nearest neighbors
// sum3 keeps track of the third-nearest neighbors
struct NeighborCounts{
    char sum1;
    char sum2;
    char sum3;
    bool operator==(const NeighborCounts& a) const{
        return (sum1 == a.sum1 && sum2 == a.sum2);
    }
};

struct NeighborInfo{
    long int first_indices[6];
    long int second_indices[12];
    long int third_indices[8];
    char total1;
    char total2;
    char total3;
};

struct CorrelationCalcParams {
	int N_sampling_max;
	bool Enable_mix_frac_method;
	bool Enable_e_method;
	bool Enable_extended_correlation_calc;
	int Correlation_cutoff_distance;
};

struct TomogramImportParams {
	double Desired_unit_size;
	bool Enable_cutoff_analysis;
	int Mixed_greyscale_width;
	double Mixed_conc;
	bool Enable_probability_analysis;
	double Probability_scaling_exponent;
	int N_extracted_segments;
};

class Morphology {
public:
    // functions
	Morphology(const int id);
    Morphology(const int length, const int width, const int height, const bool enable_z_periodic_boundary, const int id);
	Morphology(const Lattice& input_lattice, const int id);
    virtual ~Morphology();
    bool calculateAnisotropies(const int N_sampling_max);
	void calculateCorrelationDistances(const CorrelationCalcParams& parameters);
	void calculateDepthDependentData(const CorrelationCalcParams& correlation_params_input);
    double calculateInterfacialAreaVolumeRatio() const;
    bool calculateInterfacialDistance();
    double calculateInterfacialVolumeFraction() const;
	void calculateMixFractions();
    bool calculateTortuosity(const char site_tye, const bool electrode_num, const bool enable_reduced_memory);
    void createCheckerboardMorphology();
    void createRandomMorphology(const std::vector<double>& mix_fractions);
    void enableThirdNeighborInteraction();
    void executeIsingSwapping(const int num_MCsteps, const double interaction_energy1, const double interaction_energy2, const bool enable_growth_pref, const int growth_direction, const double additional_interaction); // bond formation algorithm
    void executeMixing(const double width, const double interfacial_conc);
    void executeSmoothing(const double smoothing_threshold, const int rescale_factor);
    std::vector<double> getCorrelationData(const char site_type) const;
	std::vector<double> getDepthCompositionData(const char site_type) const;
	std::vector<double> getDepthDomainSizeData(const char site_type) const;
	std::vector<double> getDepthIVData(const char site_type) const;
    double getDomainSize(const char site_type) const;
    double getDomainAnisotropy(const char site_type) const;
    int getHeight() const;
    std::vector<double> getInterfacialHistogram(const char site_type) const;
    double getIslandVolumeFraction(const char site_type) const;
    int getLength() const;
    double getMixFraction(const char site_type) const;
    std::vector<float> getTortuosityData(const char site_type) const;
    std::vector<double> getTortuosityHistogram(const char site_type) const;
    int getWidth() const;

	//! \brief imports a tomogram dataset using a combination of information from an xml metadata file and a set of import parameters
	//! \param info_filename is the name of the xml metadata file
	//! \param data_filename is the name of the raw morphology data file
	//! \param params is the TomogramImportParams data structure that contains the additional information needed to handle the tomogram data
	//! \returns a vector of Morphology objects that consists of a series of subsections of the original tomogram data
	std::vector<Morphology> importTomogramMorphologyFile(const std::string& info_filename, const std::string& data_filename, const TomogramImportParams& params);

    bool importMorphologyFile(std::ifstream& infile);
	void outputCompositionMaps(std::ofstream& outfile) const;
    void outputMorphologyFile(std::ofstream& outfile, const bool enable_export_compressed_files) const;
    void outputMorphologyCrossSection(std::ofstream& outfile) const;
	void outputTortuosityMaps(std::ofstream& outfile) const;
    void shrinkLattice(const int rescale_factor);
    void stretchLattice(const int rescale_factor);
protected:
private:
    // custom data structures
    struct Node{
        long int neighbor_indices[26];
        char neighbor_distances_sq[26];
        float distance_est;
        long int site_index;
    };
    struct NodeIteratorCompare{
        bool operator()(const std::vector<Node>::const_iterator& a, const std::vector<Node>::const_iterator& b) const{
            if(a==b){
                return false;
            }
            else if(a->distance_est==b->distance_est){
                return a < b;
            }
            else{
                return a->distance_est < b->distance_est;
            }
        }
    };
    // properties
    int ID;
    std::vector<double> Mix_fractions; // Volume fraction of each component to total
    bool Enable_third_neighbor_interaction;
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
	std::mt19937 gen;
    // functions
	void addSiteType(const char site_type);
    double calculateAdditionalEnergyChange(const long int site_index_main, const long int site_index_neighbor,const int growth_direction,const double additional_interaction) const;
    bool calculateAnisotropy(const std::vector<long int>& correlation_sites, const char site_type,const int cutoff_distance);
	double calculateCorrelationDistance(const std::vector<long int>& correlation_sites, std::vector<double>& correlation_data, const char site_type, const double mix_fraction, const int cutoff_distance, const CorrelationCalcParams& params);
    double calculateDissimilarFraction(const Coords& coords, const int rescale_factor) const;
    double calculateEnergyChangeSimple(const long int site_index1, const long int site_index2, const double interaction_energy1, const double interaction_energy2);
    double calculateEnergyChange(const Coords& coords1, const Coords& coords2,const double interaction_energy1,const double interaction_energy2) const;
    NeighborCounts calculateNeighborCounts(const Coords& coords) const;
    bool calculatePathDistances(std::vector<float>& path_distances);
    bool calculatePathDistances_ReducedMemory(std::vector<float>& path_distances);
    void createNode(Node& node,const Coords& coords);
    void getSiteSampling(std::vector<long int>& sites, const char site_type, const int N_sites);
	void getSiteSamplingZ(std::vector<long int>& sites, const char site_type, const int N_sites, const int z);
	int getSiteTypeIndex(const char site_type) const;
    void initializeNeighborInfo();
	bool isNearInterface(const Coords& coords, const double distance) const;
	double rand01();
    void updateNeighborCounts(const long int site_index1, const long int site_index2);
};

#endif // MORPHOLOGY_H
