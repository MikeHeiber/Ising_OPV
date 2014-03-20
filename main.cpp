// Copyright (c) 2014 Michael C. Heiber, The University of Akron, and Julius-Maximilains Universtiy of Würzburg
// For more information, see the LICENSE file that accompanies this software.

#include <iostream>
#include <fstream>
#include <mpi.h>
#include "Morphology.h"

using namespace std;

struct Input_Parameters{
    // General
    int length;
    int width;
    int height;
    double mix_fraction; // volume fraction of donor
    double interaction_energy;
    double target_size;
    int MC_steps;
    // Import Morphology Options
    bool enable_import_morphology;
    string morphology_import_prefix;
    // Smoothing Options
    bool enable_smoothing;
    double smoothing_threshold;
    // Rescale Options
    bool enable_rescale;
    int rescale_factor;
    // Analysis Options
    bool enable_correlation_calc;
    bool enable_tortuosity_calc;
};

double array_avg(const double data[],int size);
double array_stdev(const double data[],int size);
bool importParameters(ifstream * parameterfile,Input_Parameters * parameters);

int main(int argc, char * argv[]){
    // Input parameters
    Input_Parameters parameters;
    // Internal parameters
    bool enable_third_neighbor = false;
    double domain_size1 = 0;
    double domain_size2 = 0;
    double ratio = 0;
    int N_steps = 0;
    double elapsedtime = 0;
    time_t start, end;
    int procid = 0;
    int nproc = 1;
    stringstream ss;
    ifstream parameter_file;
    ifstream morphology_input_file;
    ofstream analysis_file;
    ofstream morphology_output_file;
    bool success;
    int cutoff_distance;
    bool correlation_status;
    int pathdata_size = 0;
    int pathdata_count = 0;
    double *domain_sizes1;
    double *domain_sizes2;
    double *ratios;
    double *times;
    double *pathdata;
    double *pathdata_all;
    int *pathdata_sizes;
    int *displacement;
    // Begin
    start = time(NULL);
    // Import parameters from text file
    parameter_file.open((argv[1]),ifstream::in);
    success = importParameters(&parameter_file,&parameters);
    parameter_file.close();
    if(!success){
        cout << "Error importing variables from file.  Program will now exit." << endl;
        return 0;
    }
    // Initialize parallel processing
    MPI::Init(argc,argv);
    nproc = MPI::COMM_WORLD.Get_size();
    procid = MPI::COMM_WORLD.Get_rank();
    // Set filename of imported morphology
    if(parameters.enable_import_morphology){
        ss << parameters.morphology_import_prefix << "_" << procid << "_uncompressed.txt";
        morphology_input_file.open(ss.str().c_str());
        ss.str("");
    }
    // Estimate the number of MC steps (N_steps) or the target domain size (target_size) based on user input
    else if(!parameters.enable_import_morphology){
        if(parameters.target_size>0){
            if(!parameters.enable_smoothing){
                // N_steps for J = 0.4, without smoothing
                if(floor(10*parameters.interaction_energy+0.5)==4){
                    N_steps = -55.65+0.823*pow(parameters.target_size,3.523);
                }
                // N_steps for J = 0.6, without smoothing
                else if(floor(10*parameters.interaction_energy+0.5)==6){
                    N_steps = -199.7+1.335*pow(parameters.target_size,3.841);
                }
                // N_steps for J = 0.8, without smoothing
                else if(floor(10*parameters.interaction_energy+0.5)==8){
                    // Regime 1
                    if(parameters.target_size<8.3){
                        N_steps = -226.2+0.563*pow(parameters.target_size,4.92);
                    }
                    // Regime 2
                    else{
                        N_steps = 288.2+9.638e-8*pow(parameters.target_size,12.3);
                    }
                }
                // N_steps for J = 1.0, without smoothing
                else if(floor(10*parameters.interaction_energy+0.5)==10){
                    N_steps = 111.6+0.408*pow(parameters.target_size,5.666);
                }
            }
            else{
                // N_steps for J = 0.4, with smoothing (assumes smoothing threshold of 0.52)
                if(floor(10*parameters.interaction_energy+0.5)==4){
                    N_steps = -100.9+0.785*pow(parameters.target_size,3.54);
                }
                // N_steps for J = 0.6, with smoothing (assumes smoothing threshold of 0.52)
                else if(floor(10*parameters.interaction_energy+0.5)==6){
                    N_steps = -267.3+1.306*pow(parameters.target_size,3.846);
                }
                // N_steps for J = 0.8, with smoothing (assumes smoothing threshold of 0.52)
                else if(floor(10*parameters.interaction_energy+0.5)==8){
                    // Regime 1
                    if(parameters.target_size<8.3){
                        N_steps = -307.4+0.463*pow(parameters.target_size,5.001);
                    }
                    // Regime 2
                    else{
                        N_steps = 2392+3.348e-9*pow(parameters.target_size,13.78);
                    }
                }
                // N_steps for J = 1.0, with smoothing (assumes smoothing threshold of 0.52)
                else if(floor(10*parameters.interaction_energy+0.5)==10){
                    N_steps = 111.6+0.408*pow(parameters.target_size,5.666);
                }
            }
        }
        else if(parameters.target_size<0){
            N_steps = parameters.MC_steps;
            // target_size for J = 0.4
            if(floor(10*parameters.interaction_energy+0.5)==4){
                parameters.target_size = 1.956+0.472*pow(N_steps,0.359);
            }
            // target_size for J = 0.6
            else if(floor(10*parameters.interaction_energy+0.5)==6){
                parameters.target_size = 1.588+0.528*pow(N_steps,0.304);
            }
            // target_size for J = 0.8
            else if(floor(10*parameters.interaction_energy+0.5)==8){
                // Regime 1
                if(N_steps<18000){
                    parameters.target_size = 0.492+1.001*pow(N_steps,0.209);
                }
                // Regime 2
                else{
                    parameters.target_size = 1.236+2.762*pow(N_steps,0.095);
                }
            }
            // target_size for J = 1.0
            else if(floor(10*parameters.interaction_energy+0.5)==10){
                parameters.target_size = -49.11+46.48*pow(N_steps,0.0185);
            }
        }
        if(procid==0){
            cout << "N_steps = " << N_steps << endl;
        }
    }
    // Create morphology data structure
    Morphology morph(parameters.length,parameters.width,parameters.height,procid);
    // Import uncompressed morphologies
    if(parameters.enable_import_morphology){
        cout << procid << ": Importing morphology from file..." << endl;
        morph.importMorphologyFile(&morphology_input_file);
        morphology_input_file.close();
        // Estimate new target domain size (target_size) if additional MC steps will be performed after importing
        if(parameters.target_size<0 && parameters.MC_steps>0){
            N_steps = parameters.MC_steps;
            int previous_steps;
            // target_size for J = 0.4
            if(floor(10*parameters.interaction_energy+0.5)==4){
                previous_steps = -55.65+0.823*pow(morph.getDomainSize(1),3.523);
                parameters.target_size = 1.956+0.472*pow(previous_steps+N_steps,0.359);
            }
            // target_size for J = 0.6
            else if(floor(10*parameters.interaction_energy+0.5)==6){
                previous_steps = -199.7+1.335*pow(morph.getDomainSize(1),3.841);
                parameters.target_size = 1.588+0.528*pow(previous_steps+N_steps,0.304);
            }
            // target_size for J = 0.8
            else if(floor(10*parameters.interaction_energy+0.5)==8){
                // Regime 1
                if(morph.getDomainSize(1)<8.3){
                    previous_steps = -226.2+0.563*pow(morph.getDomainSize(1),4.92);
                    parameters.target_size = 0.492+1.001*pow(previous_steps+N_steps,0.209);
                }
                // Regime 2
                else{
                    previous_steps = 288.2+9.638e-8*pow(morph.getDomainSize(1),12.3);
                    parameters.target_size = 1.236+2.762*pow(previous_steps+N_steps,0.095);
                }
            }
            // target_size for J = 1.0
            else if(floor(10*parameters.interaction_energy+0.5)==10){
                previous_steps = 111.6+0.408*pow(morph.getDomainSize(1),5.666);
                parameters.target_size = -49.11+46.48*pow(previous_steps+N_steps,0.0185);
            }
        }
        else{
            parameters.target_size = morph.getDomainSize(1);
            N_steps = 0;
        }
        if(procid==0){
            cout << "N_steps = " << N_steps << endl;
        }
    }
    // If not importing a morphology, generate random blend
    else{
        cout << procid << ": Generating initial random morphology..." << endl;
        morph.createRandomMorphology(parameters.mix_fraction,time(0)+procid);
    }
    // Execute phase separation
    if(N_steps>0){
        cout << procid << ": Executing site swapping..." << endl;
        morph.executeIsingSwapping(N_steps,parameters.interaction_energy,enable_third_neighbor);
    }
    // Perform domain smoothing
    if(parameters.enable_smoothing){
        cout << procid << ": Executing smoothing..." << endl;
        morph.executeSmoothing(parameters.smoothing_threshold);
    }
    // Perform lattice rescaling
    if(parameters.enable_rescale){
        cout << procid << ": Rescaling lattice by a factor of " << parameters.rescale_factor << " ..." << endl;
        morph.rescaleLattice(parameters.rescale_factor);
        parameters.length = parameters.length*parameters.rescale_factor;
        parameters.width = parameters.width*parameters.rescale_factor;
        parameters.height = parameters.height*parameters.rescale_factor;
        parameters.target_size = parameters.target_size*parameters.rescale_factor;
        // Perform final smoothing after rescaling
        if(parameters.enable_smoothing){
            cout << procid << ": Executing final smoothing..." << endl;
            morph.executeSmoothing(parameters.smoothing_threshold);
        }
    }
    // Calculate domain sizes
    if(parameters.enable_correlation_calc){
        cutoff_distance = (int)floor(parameters.target_size+0.5)+1;
        cout << procid << ": Calculating the domain size..." << endl;
        correlation_status = morph.calculateCorrelationDistance(cutoff_distance);
        while(!correlation_status){
            if(cutoff_distance<parameters.length || cutoff_distance<parameters.width){
                cutoff_distance++;
                correlation_status = morph.calculateCorrelationDistance(cutoff_distance);
            }
            else{
                break;
            }
        }
        if(correlation_status){
            domain_size1 = morph.getDomainSize(1);
            domain_size2 = morph.getDomainSize(2);
        }
    }
    // Calculate interfacial area to volume ratio
    ratio = morph.calculateInterfacialArea()/(parameters.length*parameters.width*parameters.height);
    // Calculate tortuosity
    if(parameters.enable_tortuosity_calc){
        cout << procid << ": Calculating tortuosity..." << endl;
        success = morph.calculateTortuosity();
        vector<double> data = morph.getTortuosityData();
        pathdata_size = (int)data.size();
        pathdata = (double *)malloc(sizeof(double)*pathdata_size);
        for(int i=0;i<pathdata_size;i++){
            pathdata[i] = data[i];
        }
    }
    // Save final morphology to a text file in uncompressed form
    cout << procid << ": Writing morphology to file..." << endl;
    if(!parameters.enable_import_morphology){
        ss << "morphology_" << procid << "_uncompressed.txt";
        morphology_output_file.open(ss.str().c_str());
        ss.str("");
    }
    else{
        ss << "morphology-modified_" << procid << "_uncompressed.txt";
        morphology_output_file.open(ss.str().c_str());
        ss.str("");
    }
    morph.outputMorphologyFile(&morphology_output_file);
    morphology_output_file.close();
    end = time(NULL);
    elapsedtime = (double) difftime(end,start)/60;
    MPI_Barrier(MPI_COMM_WORLD);
    // Gather morphology data from each processor and calculate set statistics
    if(procid==0){
        cout << "Collecting morphology analysis data..." << endl;
        pathdata_sizes = NULL;
        pathdata_sizes = (int *)malloc(sizeof(int)*nproc);
    }
    if(parameters.enable_tortuosity_calc){
        MPI_Gather(&pathdata_size,1,MPI_INT,pathdata_sizes,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    if(procid==0){
        if(parameters.enable_tortuosity_calc){
            for(int i=0;i<nproc;i++){
                pathdata_count += pathdata_sizes[i];
            }
        }
        domain_sizes1 = NULL;
        domain_sizes2 = NULL;
        ratios = NULL;
        times = NULL;
        pathdata_all = NULL;
        displacement = NULL;
        domain_sizes1 = (double *)malloc(sizeof(double)*nproc);
        domain_sizes2 = (double *)malloc(sizeof(double)*nproc);
        ratios = (double *)malloc(sizeof(double)*nproc);
        times = (double *)malloc(sizeof(double)*nproc);
        if(parameters.enable_tortuosity_calc){
            pathdata_all = (double *)malloc(sizeof(double)*pathdata_count);
            displacement = (int *)malloc(sizeof(int)*nproc);
            displacement[0] = 0;
            for(int i=1;i<nproc;i++){
                displacement[i] = displacement[i-1] + pathdata_sizes[i-1];
            }
        }
    }
    if(parameters.enable_correlation_calc){
        MPI_Gather(&domain_size1,1,MPI_DOUBLE,domain_sizes1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(&domain_size2,1,MPI_DOUBLE,domain_sizes2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    MPI_Gather(&ratio,1,MPI_DOUBLE,ratios,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(&elapsedtime,1,MPI_DOUBLE,times,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(parameters.enable_tortuosity_calc){
        MPI_Gatherv(pathdata,pathdata_size,MPI_DOUBLE,pathdata_all,pathdata_sizes,displacement,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    // Output final morphology set analysis to a text file
    if(procid==0){
        if(!parameters.enable_import_morphology){
            ss << "analysis_summary.txt";
            analysis_file.open(ss.str().c_str());
            ss.str("");
        }
        else{
            ss << "analysis_summary_mod.txt";
            analysis_file.open(ss.str().c_str());
            ss.str("");
        }
        analysis_file << "N_morphologies,length,width,height,interaction_energy,MC_steps,domain1_size_avg,domain1_size_stdev,domain2_size_avg,domain2_size_stdev,";
        analysis_file << "interfacial_area_volume_ratio_avg,interfacial_area_volume_ratio_stdev,tortuosity_avg,tortuosity_stdev,calc_time_avg(min),calc_time_stdev(min)" << endl;
        analysis_file << nproc << "," << parameters.length << "," << parameters.width << "," << parameters.height << ",";
        analysis_file << parameters.interaction_energy << "," << N_steps << "," << array_avg(domain_sizes1,nproc) << "," << array_stdev(domain_sizes1,nproc) << ",";
        analysis_file << array_avg(domain_sizes2,nproc) << "," << array_stdev(domain_sizes2,nproc) << "," << array_avg(ratios,nproc) << "," << array_stdev(ratios,nproc) << ",";
        analysis_file << array_avg(pathdata_all,pathdata_count) << "," << array_stdev(pathdata_all,pathdata_count) << "," << array_avg(times,nproc) << "," << array_stdev(times,nproc) << endl;
        analysis_file.close();
    }
    MPI_Finalize();
    cout << procid << ": Finished!" << endl;
    return 0;
}

bool importParameters(ifstream * parameterfile,Input_Parameters * parameters){
    string line;
    string var;
    size_t pos;
    vector<string> stringvars;
    while((*parameterfile).good()){
        getline(*parameterfile,line);
        if((line.substr(0,2)).compare("--")!=0 && (line.substr(0,2)).compare("##")!=0){
            pos = line.find("/",0);
            var = line.substr(0,pos-1);
            stringvars.push_back(var);
        }
    }
    int i = 0;
    // General Parameters
    (*parameters).length = atoi(stringvars[i].c_str());
    i++;
    (*parameters).width = atoi(stringvars[i].c_str());
    i++;
    (*parameters).height = atoi(stringvars[i].c_str());
    i++;
    (*parameters).mix_fraction = atof(stringvars[i].c_str());
    i++;
    (*parameters).interaction_energy = atof(stringvars[i].c_str());
    i++;
    (*parameters).target_size = atof(stringvars[i].c_str());
    i++;
    (*parameters).MC_steps = atoi(stringvars[i].c_str());
    i++;
    //enable_import_morphology
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_import_morphology = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_import_morphology = false;
    }
    else{
        cout << "Error setting morphology import options" << endl;
        return false;
    }
    i++;
    (*parameters).morphology_import_prefix = stringvars[i];
    i++;
    //enable_smoothing
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_smoothing = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_smoothing = false;
    }
    else{
        cout << "Error setting morphology smoothing options" << endl;
        return false;
    }
    i++;
    (*parameters).smoothing_threshold = atof(stringvars[i].c_str());
    i++;
    //enable_rescale
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_rescale = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_rescale = false;
    }
    else{
        cout << "Error setting morphology rescale options" << endl;
        return false;
    }
    i++;
    (*parameters).rescale_factor = atoi(stringvars[i].c_str());
    i++;
    //enable_correlation_calc
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_correlation_calc = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_correlation_calc = false;
    }
    else{
        cout << "Error setting correlation calculation options" << endl;
        return false;
    }
    i++;
    //enable_tortuosity_calc
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_tortuosity_calc = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_tortuosity_calc = false;
    }
    else{
        cout << "Error setting tortuosity calculation options" << endl;
        return false;
    }
    i++;
    return true;
}

double array_avg(const double data[],int size){
    double sum = 0;
    for(int i=0;i<size;i++){
        sum += data[i];
    }
    return sum/size;
}

double array_stdev(const double data[],int size){
    double sum = 0;
    double avg = array_avg(data,size);
    for(int i=0;i<size;i++){
        sum += pow(data[i]-avg,2);
    }
    return sqrt(sum/(size-1));
}
