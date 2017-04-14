// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.

#ifndef MORPHOLOGY_H
#define MORPHOLOGY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <boost/random.hpp>
#include <list>
#include <set>
#include <iterator>

using namespace std;

struct Coords{
    int x;
    int y;
    int z;
};

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
    Morphology(int length,int width,int height,bool enable_z_periodic_boundary,int procid);
    virtual ~Morphology();
    bool calculateAnisotropies(int cutoff_distance,int N_sampling_max);
    bool calculateCorrelationDistance(int cutoff_distance, bool enable_extended_calc,int N_sampling_max);
    double calculateInterfacialArea();
    bool calculateInterfacialDistance();
    bool calculateInterfacialDistanceOld();
    double calculateInterfacialVolume();
    bool calculateTortuosity(bool enable_reduced_memory);
    void createCheckerboardMorphology();
    void createRandomMorphology(double mix_fraction);
    void enableThirdNeighborInteraction();
    void executeIsingSwapping(int num_MCsteps,double interaction_energy1,double interaction_energy2,bool enable_growth_pref,int growth_direction,double additional_interaction); // bond formation algorithm
    void executeMixing(double width,double interfacial_conc);
    void executeSmoothing(double smoothing_threshold, int rescale_factor);
    vector<double> getCorrelationData(char site_type);
    double getDomainSize(char site_type);
    double getDomainSpacing(char site_type);
    double getDomainAnisotropy(char site_type);
    int getHeight();
    vector<double> getInterfacialHistogram(char site_type);
    double getIslandVolume(char site_type);
    int getLength();
    double getMixFraction();
    vector<float> getTortuosityData(char site_type);
    vector<double> getTortuosityHistogram(char site_type);
    int getWidth();
    bool importMorphologyFile(ifstream * input,bool compressed_files);
    bool outputMorphologyFile(ofstream * output,bool enable_export_compressed_files);
    bool outputMorphologyCrossSection(ofstream * output);
    void shrinkLattice(int rescale_factor);
    void stretchLattice(int rescale_factor);
protected:
private:
    // custom data structures
    struct Site {
        char type;
    };
    struct Node{
        long int neighbor_indices[26];
        char neighbor_distances_sq[26];
        float distance_est;
        long int site_index;
    };
    struct NodeIteratorCompare{
        bool operator()(const vector<Node>::const_iterator& a, const vector<Node>::const_iterator& b) const{
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
    int Length;
    int Width;
    int Height;
    double Mix_fraction; // Fraction of donor
    bool Enable_z_periodic_boundary;
    bool Enable_third_neighbor_interaction;
    vector<double> Correlation1;
    vector<double> Correlation2;
    vector<float> TortuosityData1;
    vector<float> TortuosityData2;
    vector<double> InterfacialHistogram1;
    vector<double> InterfacialHistogram2;
    vector<double> TortuosityHistogram1;
    vector<double> TortuosityHistogram2;
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
    vector<Site> lattice;
    vector<long int> interfacial_sites;
    vector<long int> Correlation_sites;
    vector<NeighborCounts> neighbor_counts;
    vector<NeighborInfo> neighbor_info;
    NeighborCounts temp_counts1;
    NeighborCounts temp_counts2;
    boost::mt19937 gen;
    // functions
    double calculateAdditionalEnergyChange(long int main_site_index,long int neighbor_site_index,int growth_direction,double additional_interaction);
    bool calculateAnisotropy(char site_type,int cutoff_distance,int N_sampling_max);
    double calculateDissimilarFraction(int x,int y,int z,int rescale_factor);
    inline int calculateDX(int x, int i);
    inline int calculateDY(int y, int j);
    inline int calculateDZ(int z, int k);
    double calculateEnergyChangeSimple(long int main_site_index,long int neighbor_site_index,double interaction_energy1,double interaction_energy2);
    double calculateEnergyChange(int x1,int y1,int z1,int x2,int y2,int z2,double interaction_energy1,double interaction_energy2);
    void calculateMixFraction();
    NeighborCounts calculateNeighborCounts(int x,int y,int z);
    bool calculatePathDistances(vector<float>& path_distances);
    bool calculatePathDistances_ReducedMemory(vector<float>& path_distances);
    void createNode(Node& node,int x,int y,int z);
    Coords getCoords(long int site_index);
    inline long int getShrinkSite(int x,int y,int z,int rescale_factor);
    inline long int getStretchSite(int x,int y,int z,int rescale_factor);
    inline long int getSite(int x,int y,int z);
    void getSiteSampling(vector<long int>& sites, int N_sites);
    void initializeNeighborInfo();
    inline double intpow(double base,int exponent);
    inline int round_int(double num);
    bool isNearInterface(int x,int y,int z,double distance);
    void updateNeighborCounts(long int site_index1,long int site_index2);
};

#endif // MORPHOLOGY_H
