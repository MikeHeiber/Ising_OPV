#ifndef MORPHOLOGY_H
#define MORPHOLOGY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <boost/random.hpp>

using namespace std;

class Morphology {
public:
    // functions
    Morphology(int length,int width,int height,int procid);
    virtual ~Morphology();
    bool calculateCorrelationDistance(int cutoff_distance);
    double calculateInterfacialArea();
    bool calculateTortuosity();
    void createRandomMorphology(double mix_fraction,int seed);
    void executeIsingSwappingAlt(int num_MCsteps,double interaction_energy,bool enable_third_neighbor_calc); // site energy algorithm
    void executeIsingSwapping(int num_MCsteps,double interaction_energy,bool enable_third_neighbor_calc); // bond formation algorithm
    void executeSmoothing(double smoothing_threshold);
    double getDomainSize(int site_type);
    int getHeight();
    int getLength();
    int getWidth();
    double getMixFraction();
    double getPathAvg();
    double getPathStdev();
    vector<double> getTortuosityData();
    vector<double> getCorrelationData(int site_type);
    int getSiteType(int x,int y,int z);
    bool importMorphologyFile(ifstream * input);
    bool outputMorphologyFile(ofstream * output);
    void rescaleLattice(int scale);
protected:
private:
    struct Site {
        int type;
        double energy;
        double path_distance;
    };
    int ProcID;
    // properties
    int Length;
    int Width;
    int Height;
    double Mix_fraction; // Fraction of donor
    vector<double> Correlation1;
    vector<double> Correlation2;
    vector<double> TortuosityData;
    bool Domain_size1_updated;
    bool Domain_size2_updated;
    double Domain_size1;
    double Domain_size2;
    bool Energies_initialized;
    double Path_avg;
    double Path_stdev;
    std::vector<Site> lattice;
    std::vector<double> lattice_temp;
    std::vector<int> interfacial_sites;
    boost::mt19937 gen;
    // functions
    double calculateEnergyChange1(int x1,int y1,int z1,int x2,int y2,int z2,double interaction_energy,bool enable_third_neighbor_calc);
    double calculateEnergyChange2(int x1,int y1,int z1,int x2,int y2,int z2,double interaction_energy,bool enable_third_neighbor_calc);
    void calculateSiteEnergies(double interaction_energy,bool enable_third_neighbor_calc);
    void calculateMixFraction();
    double calculateSiteEnergy(int x,int y,int z,double interaction_energy,bool enable_third_neighbor_calc);
    void calculateTempEnergies(int x1,int y1,int z1,int di,int dj,int dk,double interaction_energy,bool enable_third_neighbor_calc);
    int countDissimilarNeighbors(int x,int y,int z);
    int countNeighbors(int z);
    void executeTempEnergies(int x,int y,int z);
    int getRescaleSite(int x,int y,int z,int scale);
    int getSite(int x,int y,int z);
    int getTempSite(int x,int y,int z);
    bool isNearInterface(int x,int y,int z,double distance);
};

#endif // MORPHOLOGY_H
