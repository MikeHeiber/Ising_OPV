// Copyright (c) 2015 Michael C. Heiber
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

using namespace std;

struct Coords{
    int x;
    int y;
    int z;
};

class Morphology {
public:
    // functions
    Morphology(int length,int width,int height,bool enable_z_periodic_boundary,int procid,int seed);
    virtual ~Morphology();
    bool calculateCorrelationDistance(int cutoff_distance,int N_sampling_max);
    double calculateInterfacialArea();
    bool calculateInterfacialDistance();
    double calculateInterfacialVolume();
    bool calculateTortuosity();
    void createRandomMorphology(double mix_fraction);
    void enableThirdNeighborInteraction();
    void executeIsingSwapping(int num_MCsteps,double interaction_energy1,double interaction_energy2); // bond formation algorithm
    void executeMixing(double width,double interfacial_conc);
    void executeSmoothing(double smoothing_threshold, int rescale_factor);
    vector<double> getCorrelationData(char site_type);
    double getDomainSize(char site_type);
    double getDomainSpacing(char site_type);
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
    void rescaleLattice(int scale);
protected:
private:
    // custom data structures
    struct Site {
        char type;
    };
    struct Node{
        float distance_est;
        unsigned long neighbor_indices[26];
        char neighbor_distances_sq[26];
    };
    struct NodeIteratorCompare{
        bool operator()(const vector<Node>::iterator& a, const vector<Node>::iterator& b) const{
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
    double Domain_size1;
    double Domain_size2;
    int Island_volume1;
    int Island_volume2;
    vector<Site> lattice;
    vector<unsigned long> interfacial_sites;
    vector<unsigned long> Correlation_sites;
    vector<Node> Node_vector;
    boost::mt19937 gen;
    // functions
    double calculateEnergyChange(int x1,int y1,int z1,int x2,int y2,int z2,double interaction_energy1,double interaction_energy2);
    void calculateMixFraction();
    double calculateDissimilarFraction(int x,int y,int z,int rescale_factor);
    Coords getCoords(unsigned long site_index);
    inline unsigned long getRescaleSite(int x,int y,int z,int scale);
    inline unsigned long getSite(int x,int y,int z);
    inline double intpow(double base,int exponent);
    inline int round_int(double num);
    bool isNearInterface(int x,int y,int z,double distance);
};

#endif // MORPHOLOGY_H
