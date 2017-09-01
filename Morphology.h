// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#ifndef MORPHOLOGY_H
#define MORPHOLOGY_H

#include "Lattice.h"
#include "Utils.h"
#include <ctime>
#include <fstream>
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

class Morphology {
public:
    // functions
    Morphology(const int length, const int width, const int height, const bool enable_z_periodic_boundary, const int procid);
    virtual ~Morphology();
    bool calculateAnisotropies(const int cutoff_distance, const int N_sampling_max);
    bool calculateCorrelationDistance(const int cutoff_distance, const bool enable_extended_calc, const int N_sampling_max);
    double calculateInterfacialArea() const;
    bool calculateInterfacialDistance();
    //bool calculateInterfacialDistanceOld();
    double calculateInterfacialVolume() const;
    bool calculateTortuosity(const bool enable_reduced_memory);
    void createCheckerboardMorphology();
    void createRandomMorphology(const double mix_fraction);
    void enableThirdNeighborInteraction();
    void executeIsingSwapping(const int num_MCsteps, const double interaction_energy1, const double interaction_energy2, const bool enable_growth_pref, const int growth_direction, const double additional_interaction); // bond formation algorithm
    void executeMixing(const double width, const double interfacial_conc);
    void executeSmoothing(const double smoothing_threshold, const int rescale_factor);
    std::vector<double> getCorrelationData(const char site_type) const;
    double getDomainSize(const char site_type) const;
    //double getDomainSpacing(char site_type) const;
    double getDomainAnisotropy(const char site_type) const;
    int getHeight() const;
    std::vector<double> getInterfacialHistogram(const char site_type) const;
    double getIslandVolume(const char site_type) const;
    int getLength() const;
    double getMixFraction() const;
    std::vector<float> getTortuosityData(const char site_type) const;
    std::vector<double> getTortuosityHistogram(const char site_type) const;
    int getWidth() const;
    bool importMorphologyFile(std::ifstream * input, const bool compressed_files);
    bool outputMorphologyFile(std::ofstream * output, const bool enable_export_compressed_files);
    bool outputMorphologyCrossSection(std::ofstream * output);
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
    double Mix_fraction; // Fraction of donor
    bool Enable_third_neighbor_interaction;
	Lattice lattice;
    std::vector<double> Correlation1;
	std::vector<double> Correlation2;
	std::vector<float> TortuosityData1;
	std::vector<float> TortuosityData2;
	std::vector<double> InterfacialHistogram1;
	std::vector<double> InterfacialHistogram2;
	std::vector<double> TortuosityHistogram1;
	std::vector<double> TortuosityHistogram2;
    bool Domain_size1_updated;
    bool Domain_size2_updated;
    bool Domain_anisotropy1_updated;
    bool Domain_anisotropy2_updated;
    double Domain_size1;
    double Domain_size2;
    double Domain_anisotropy1;
    double Domain_anisotropy2;
    int Island_volume1;
    int Island_volume2;
	std::vector<long int> interfacial_sites;
	std::vector<long int> Correlation_sites;
	std::vector<NeighborCounts> neighbor_counts;
	std::vector<NeighborInfo> neighbor_info;
    NeighborCounts temp_counts1;
    NeighborCounts temp_counts2;
	std::mt19937 gen;
    // functions
    double calculateAdditionalEnergyChange(const long int site_index_main, const long int site_index_neighbor,const int growth_direction,const double additional_interaction) const;
    bool calculateAnisotropy(const char site_type,const int cutoff_distance,const int N_sampling_max);
    double calculateDissimilarFraction(const Coords& coords, const int rescale_factor) const;
    double calculateEnergyChangeSimple(const long int site_index1, const long int site_index2, const double interaction_energy1, const double interaction_energy2);
    double calculateEnergyChange(const Coords& coords1, const Coords& coords2,const double interaction_energy1,const double interaction_energy2) const;
    void calculateMixFraction();
    NeighborCounts calculateNeighborCounts(const Coords& coords) const;
    bool calculatePathDistances(std::vector<float>& path_distances);
    bool calculatePathDistances_ReducedMemory(std::vector<float>& path_distances);
    void createNode(Node& node,const Coords& coords);
    void getSiteSampling(std::vector<long int>& sites, const int N_sites);
    void initializeNeighborInfo();
	bool isNearInterface(const Coords& coords, const double distance) const;
	double rand01();
    void updateNeighborCounts(const long int site_index1, const long int site_index2);
};

#endif // MORPHOLOGY_H
