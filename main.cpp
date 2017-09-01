// Copyright (c) 2017 Michael C. Heiber
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

struct Input_Parameters{
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
    bool Enable_smoothing; // choose whether to perform domain smoothing
    double Smoothing_threshold;  // specify the degree of smoothing
    // Rescale Options
    bool Enable_rescale; // choose whether to perform lattice rescaling
    int Rescale_factor; // specify the rescale factor to be used (must be an integer greater than 1)
    bool Enable_shrink; // chose whether to shrink the lattice instead of expand it
    // Interfacial Mixing Options
    bool Enable_interfacial_mixing; // choose whether to perform interfacial mixing
    double Interface_width;  // specify the interfacial width
    double Interface_conc; // specify the mixing concentration in the interfacial region
    // Analysis Options
    bool Enable_analysis_only;  // choose whether to only perform analysis on an imported morphology
    bool Enable_correlation_calc; // choose whether to perform the domain size calculation using the pair-pair correlation method
    int N_sampling_max; // specify the maximum number of sites to be sampled for calculating the correlation function
    bool Enable_extended_correlation_calc; // choose whether to extend the correlation function calculation to the second correlation maximum
    bool Enable_interfacial_distance_calc; // choose whether to calculate the interfacial distance histograms
    bool Enable_tortuosity_calc; // choose whether to calculate the tortuosity histograms, end-to-end tortuosity, and island volume fraction
    bool Enable_reduced_memory_tortuosity_calc; //choose whether or not to enable a tortuosity calculation method that takes longer, but uses less memory
    // Other Options
    bool Enable_checkerboard_start; //choose whether or not to start from a alternating checkerboard-like configuration instead of a random blend (creates 0.5 mix fraction)
    bool Enable_growth_pref;
    int Growth_direction;
    double Additional_interaction;
};

double array_median(const double data[],int size);
int array_which_median(const double data[],int size);
bool importParameters(ifstream * parameterfile,Input_Parameters& params);

int main(int argc, char * argv[]){
    // Input parameters
    Input_Parameters parameters;
    // Internal parameters
    string version = "v4.0-alpha.1";
    bool Enable_import_morphology;
    double mix_ratio = 0;
    double domain_size1 = 0;
    double domain_size2 = 0;
    double domain_anisotropy1 = 0;
    double domain_anisotropy2 = 0;
    double iav_ratio = 0;
    double iv_ratio = 0;
    double island_ratio1 = 0;
    double island_ratio2 = 0;
    int N_steps = 0;
    double elapsedtime = 0;
    time_t start_time, end_time;
    int procid = 0;
    int nproc = 1;
    stringstream ss;
    ifstream parameter_file;
    ifstream morphology_input_file;
    ofstream analysis_file;
    ofstream morphology_output_file;
    ofstream morphology_cross_section_file;
    ofstream correlationfile_avg;
    ofstream correlationfile;
    ofstream interfacial_dist_hist_file;
    ofstream tortuosity_hist_file;
    ofstream path_data1_file;
    ofstream path_data2_file;
    bool success;
    int cutoff_distance;
    double *mix_ratios = NULL;
    double *domain_sizes1 = NULL;
    double *domain_sizes2 = NULL;
    double *domain_anisotropies1 = NULL;
    double *domain_anisotropies2 = NULL;
    double *iav_ratios = NULL;
    double *iv_ratios = NULL;
    double *island_ratios1 = NULL;
    double *island_ratios2 = NULL;
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
    string path = "";
    string input_morphology, input_file_path, compression_str, filename_prefix;
    bool is_file_compressed;
    // Begin
    start_time = time(NULL);
    // Initialize parallel processing.
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    // Import parameters from text file.
    if(argc == 2){
        Enable_import_morphology=false;
    }
    else if(argc == 3){
        Enable_import_morphology=true;
    }
    else{
        cout << procid << ": Incorrect arguments. There needs to be at least one argument which is the filename of a parameter file. Second optional argument needs to be a morphology file of the form 'filename_procid_compressed.txt' or 'filename_procid_uncompressed.txt'. Program will exit now!" << endl;
        return 0;
    }
    parameter_file.open((argv[1]),ifstream::in);
    success = importParameters(&parameter_file,parameters);
    parameter_file.close();
    if(success){
        cout << procid << ": Parameter file found and loaded successfully." << endl;
    }
    else{
        cout << procid << ": Error importing variables from file!  Program will exit now." << endl;
        return 0;
    }
    if(parameters.Enable_analysis_only && !Enable_import_morphology){
        cout << procid << ": Error!  The 'analysis only' option can only be used when a morphology is imported from a file." << endl;
        return 0;
    }
    if(parameters.Enable_analysis_only){
        cout << procid << ": Warning! Only morphology analysis will be performed." << endl;
    }
    // Create morphology data structure.
    Morphology morph(parameters.Length,parameters.Width,parameters.Height,parameters.Enable_periodic_z,procid);
    // Import morphology if enabled.
    if(Enable_import_morphology){
        size_t pos_comp, pos_path, pos_id;
        // Get filename of imported morphology from command line arguments.
        input_morphology = (argv[2]);
        // Separate filepath and differentiate between Windows (1) and Linux (2) filepaths.
        int os;
        // Linux filepaths
        if(input_morphology.substr(0,1).compare("/")==0){
            pos_path = input_morphology.find_last_of("/");
            input_file_path = input_morphology.substr(0,pos_path)+"/";
            os=1;
        }
        // Windows filepaths
        else if(input_morphology.substr(0,2).compare("\\")==0){
            pos_path = input_morphology.find_last_of("\\");
            input_file_path = input_morphology.substr(0,pos_path)+"\\";
            os=2;
        }
        // Default to be used when the morphology data file is in the current working directory
        else{
            input_file_path = "";
            pos_path = -1;
            os=1;
        }
        // Determine if the input file is in the compressed format or not.
        pos_comp = input_morphology.find("compressed.txt");
        if(pos_comp==string::npos){
            cout << procid << ": Error! Format of the file needs to be 'filename_#_compressed.txt' or 'filename_#_uncompressed.txt'. Program will exit now." << endl;
            return 0;
        }
        compression_str = input_morphology.substr(pos_comp-1);
        if(compression_str.substr(0,1).compare("_")==0){
            is_file_compressed = true;
            compression_str = input_morphology.substr(pos_comp-1);
            filename_prefix = input_morphology.substr((pos_path+os),pos_comp-1-(pos_path+os));
        }
        else if(compression_str.substr(0,1).compare("n")==0){
            is_file_compressed = false;
            compression_str = input_morphology.substr(pos_comp-3);
            filename_prefix = input_morphology.substr((pos_path+os),pos_comp-3-(pos_path+os));
        }
        else{
            cout << procid << ": Error! Format of the file needs to be 'filename_#_compressed.txt' or 'filename_#_uncompressed.txt'. Program will exit now." << endl;
            return 0;
        }
        // Separate filename prefix from ID number.
        pos_id = filename_prefix.find_last_of("_");
        if(pos_id==string::npos){
            cout << procid << ": Error! Format of the file needs to be 'filename_#_compressed.txt' or 'filename_#_uncompressed.txt'. Program will exit now." << endl;
            return 0;
        }
        else{
            filename_prefix =input_morphology.substr((pos_path+os),pos_id+1);
        }
        ss << input_file_path << filename_prefix << procid << compression_str;
        morphology_input_file.open(ss.str().c_str());
        if(morphology_input_file.is_open()) {
            cout << procid << ": Morphology file successfully opened!" << endl;
        }
        else{
            cout << procid << ": Opening morphology file failed! Program will exit now!" << endl;
            return 0;
        }
        ss.str("");
        if(!is_file_compressed){
            cout << procid << ": Importing morphology from uncompressed file..." << flush;
        }
        else{
            cout << procid << ": Importing morphology from compressed file..." << flush;
        }
        // Import the morphology from the given data file.
        morph.importMorphologyFile(&morphology_input_file,is_file_compressed);
        morphology_input_file.close();
        cout << procid << ": Morphology import complete!" << endl;
    }
    // Create new morphology if import is disabled.
    else{
        if(parameters.Enable_checkerboard_start){
            cout << procid << ": Generating initial checkerboard morphology..." << endl;
            morph.createCheckerboardMorphology();
        }
        else{
            cout << procid << ": Generating initial random morphology..." << endl;
            morph.createRandomMorphology(parameters.Mix_fraction);
        }
    }
    // Determine if any phase separation is to be executed on the morphology.
    if(parameters.MC_steps>0  && !parameters.Enable_analysis_only){
        N_steps = parameters.MC_steps;
    }
    else{
        N_steps = 0;
    }
    // Execute phase separation through Ising swapping.
    if(N_steps>0){
        cout << procid << ": Executing site swapping for " << N_steps << " MC steps..." << endl;
        morph.executeIsingSwapping(N_steps,parameters.Interaction_energy1,parameters.Interaction_energy2,parameters.Enable_growth_pref,parameters.Growth_direction,parameters.Additional_interaction);
    }
    // Perform lattice rescaling and domain smoothing if enabled.
    if(parameters.Enable_rescale && !parameters.Enable_analysis_only){
        if(parameters.Enable_shrink){
            cout << procid << ": Initial blend ratio is " << morph.getMixFraction() << endl;
            if(parameters.Enable_smoothing){
                cout << procid << ": Executing standard smoothing with a smoothing threshold of " << parameters.Smoothing_threshold << "..." << endl;
                morph.executeSmoothing(parameters.Smoothing_threshold,2);
                cout << procid << ": Blend ratio after smoothing is " << morph.getMixFraction() << endl;
            }
            cout << procid << ": Shrinking lattice by a factor of " << parameters.Rescale_factor << " ..." << endl;
            morph.shrinkLattice(parameters.Rescale_factor);
            parameters.Length = parameters.Length/parameters.Rescale_factor;
            parameters.Width = parameters.Width/parameters.Rescale_factor;
            parameters.Height = parameters.Height/parameters.Rescale_factor;
            cout << procid << ": Blend ratio after shrinking lattice is " << morph.getMixFraction() << endl;
        }
        else{
            cout << procid << ": Expanding lattice by a factor of " << parameters.Rescale_factor << " ..." << endl;
            morph.stretchLattice(parameters.Rescale_factor);
            parameters.Length = parameters.Length*parameters.Rescale_factor;
            parameters.Width = parameters.Width*parameters.Rescale_factor;
            parameters.Height = parameters.Height*parameters.Rescale_factor;
        }
    }
    if(parameters.Enable_smoothing && !parameters.Enable_shrink && !parameters.Enable_analysis_only){
        if(!parameters.Enable_rescale){
            cout << procid << ": Executing standard smoothing with a smoothing threshold of " << parameters.Smoothing_threshold << "..." << endl;
            morph.executeSmoothing(parameters.Smoothing_threshold,2);
        }
        else{
            cout << procid << ": Executing rescale factor dependent smoothing with a smoothing threshold of " << parameters.Smoothing_threshold << "..." << endl;
            morph.executeSmoothing(parameters.Smoothing_threshold,parameters.Rescale_factor);
        }
    }
    // Perform interfacial mixing if enabled.
    if(parameters.Enable_interfacial_mixing && !parameters.Enable_analysis_only){
        cout << "Executing interfacial mixing..." << endl;
        morph.executeMixing(parameters.Interface_width,parameters.Interface_conc);
    }
    // Calculate domain size if enabled.
    if(parameters.Enable_correlation_calc){
        cutoff_distance = 5;
        if(parameters.Enable_extended_correlation_calc){
            cout << procid << ": Calculating the domain size using the extended pair-pair correlation function..." << endl;
        }
        else{
            cout << procid << ": Calculating the domain size from the pair-pair correlation function..." << endl;
        }
        success = false;
        // The correlation function calculation is called with an increasing cutoff distance until successful.
        while(!success){
            if(2*cutoff_distance>morph.getLength() || 2*cutoff_distance>morph.getWidth() || 2*cutoff_distance>morph.getHeight()){
                break;
            }
            success = morph.calculateCorrelationDistance(cutoff_distance,parameters.Enable_extended_correlation_calc,parameters.N_sampling_max);
            cutoff_distance++;
        }
        if(success){
            domain_size1 = morph.getDomainSize((char)1);
            domain_size2 = morph.getDomainSize((char)2);
            ss << path << "correlation_data_" << procid << ".txt";
            correlationfile.open(ss.str().c_str());
            ss.str("");
            vector<double> correlation_data1 = morph.getCorrelationData((char)1);
            vector<double> correlation_data2 = morph.getCorrelationData((char)2);
            for(int i=0;i<correlation_data1.size();i++){
                correlationfile << 0.5*(double)i << "," << correlation_data1[i] << "," << correlation_data2[i] << endl;
            }
            correlationfile.close();
        }
        success = false;
        cout << procid << ": Calculating the domain anisotropy..." << endl;
        while(!success){
            if(2*cutoff_distance>morph.getLength() && 2*cutoff_distance>morph.getWidth() && 2*cutoff_distance>morph.getHeight()){
                break;
            }
            success = morph.calculateAnisotropies(cutoff_distance,parameters.N_sampling_max);
            cutoff_distance++;
        }
        if(success){
            domain_anisotropy1 = morph.getDomainAnisotropy((char)1);
            domain_anisotropy2 = morph.getDomainAnisotropy((char)2);
        }
        else{
            cout << procid << ": Warning! Could not calculate the domain anisotropy." << endl;
        }
    }
    // Calculate interfacial distance histogram if enabled.
    if(parameters.Enable_interfacial_distance_calc){
        cout << procid << ": Calculating the interfacial distance histogram..." << endl;
        morph.calculateInterfacialDistance();
    }
    // Calculate interfacial area to volume ratio.
    iav_ratio = morph.calculateInterfacialArea()/(morph.getLength()*morph.getWidth()*morph.getHeight());
    // Calculate interfacial volume to total volume ratio.
    iv_ratio = morph.calculateInterfacialVolume()/(morph.getLength()*morph.getWidth()*morph.getHeight());
    // Get Final Mix ratio
    mix_ratio = morph.getMixFraction();
    // Calculate end-to-end tortuosity, tortuosity histogram, and island volume fraction.
    if(parameters.Enable_tortuosity_calc){
        if(parameters.Enable_reduced_memory_tortuosity_calc){
            cout << procid << ": Calculating tortuosity using the reduced memory method..." << endl;
        }
        else{
            cout << procid << ": Calculating tortuosity using the standard method..." << endl;
        }
        success = morph.calculateTortuosity(parameters.Enable_reduced_memory_tortuosity_calc);
        if(!success){
            cout << procid << ": Error calculating tortuosity! Program will exit now." << endl;
            return 0;
        }
        // Calculate island volume ratio.
        island_ratio1 = (double)morph.getIslandVolume((char)1)/(morph.getLength()*morph.getWidth()*morph.getHeight());
        island_ratio2 = (double)morph.getIslandVolume((char)2)/(morph.getLength()*morph.getWidth()*morph.getHeight());
    }
    // Save final morphology to a text file.
    if(!parameters.Enable_analysis_only){
        cout << procid << ": Writing morphology to file..." << endl;
        if(!Enable_import_morphology){
            if(!parameters.Enable_export_compressed_files){
                ss << path << "morphology_" << procid << "_uncompressed.txt";
                morphology_output_file.open(ss.str().c_str());
                ss.str("");
                morphology_output_file << "Ising_OPV " << version << " - uncompressed format" << endl;
            }
            else{
                ss << path << "morphology_" << procid << "_compressed.txt";
                morphology_output_file.open(ss.str().c_str());
                ss.str("");
                morphology_output_file << "Ising_OPV " << version << " - compressed format" << endl;
            }
        }
        else{
            if(!parameters.Enable_export_compressed_files){
                ss << path << filename_prefix << "mod_" << procid << "_uncompressed.txt";
                morphology_output_file.open(ss.str().c_str());
                ss.str("");
                morphology_output_file << "Ising_OPV " << version << " - uncompressed format" << endl;
            }
            else{
                ss << path << filename_prefix << "mod_" << procid << "_compressed.txt";
                morphology_output_file.open(ss.str().c_str());
                ss.str("");
                morphology_output_file << "Ising_OPV " << version << " - compressed format" << endl;
            }
        }
        morph.outputMorphologyFile(&morphology_output_file,parameters.Enable_export_compressed_files);
        morphology_output_file.close();
    }
    // Save the cross-section of the x=0 plane to a file if enabled.
    if(parameters.Enable_export_cross_section){
        ss << path << "morphology_" << procid << "_cross_section.txt";
        morphology_cross_section_file.open(ss.str().c_str());
        ss.str("");
        morph.outputMorphologyCrossSection(&morphology_cross_section_file);
        morphology_cross_section_file.close();
    }
    // Morphology generation is now finished.
    end_time = time(NULL);
    elapsedtime = (double) difftime(end_time,start_time)/60;
    // All processors must finish with morphology generation before analysis can begin.
    MPI_Barrier(MPI_COMM_WORLD);
    // Gather morphology data from each processor and calculate set statistics.
    if(procid==0){
        cout << "Collecting morphology analysis data from each processor..." << endl;
    }
    if(parameters.Enable_tortuosity_calc){
        // Get path data for end-to-end tortuosity calculation.
        vector<float> data = morph.getTortuosityData((char)1);
        // Determine size of the data vector.
        pathdata1_size = (int)data.size();
        // Allocate an array to store the data.
        pathdata1 = (float *)malloc(sizeof(float)*pathdata1_size);
        // Put data from the vector into the array.
        for(int i=0;i<pathdata1_size;i++){
            pathdata1[i] = data[i];
        }
        data.clear();
        // Repeat for type 2 path data.
        data = morph.getTortuosityData((char)2);
        pathdata2_size = (int)data.size();
        pathdata2 = (float *)malloc(sizeof(float)*pathdata2_size);
        for(int i=0;i<pathdata2_size;i++){
            pathdata2[i] = data[i];
        }
        data.clear();
        // Create array on the root processor that will contain the size of the path data arrays from each processor.
        if(procid==0){
            pathdata1_sizes = (int *)malloc(sizeof(int)*nproc);
            pathdata2_sizes = (int *)malloc(sizeof(int)*nproc);
        }
        // Gather the size of the data arrays from all processors into the size array on the root processor.
        MPI_Gather(&pathdata1_size,1,MPI_INT,pathdata1_sizes,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Gather(&pathdata2_size,1,MPI_INT,pathdata2_sizes,1,MPI_INT,0,MPI_COMM_WORLD);
        // Calculate the average tortuosity histograms.
        tortuosity_hist1_vect = MPI_calculateVectorAvg(morph.getTortuosityHistogram((char)1));
        tortuosity_hist2_vect = MPI_calculateVectorAvg(morph.getTortuosityHistogram((char)2));
    }
    // Calculate the average interfacial distance histograms.
    if(parameters.Enable_interfacial_distance_calc){
        interfacial_dist_hist1_vect = MPI_calculateVectorAvg(morph.getInterfacialHistogram((char)1));
        interfacial_dist_hist2_vect = MPI_calculateVectorAvg(morph.getInterfacialHistogram((char)2));
    }
    // Calculate the average pair-pair correlation functions.
    if(parameters.Enable_correlation_calc){
        correlation1_vect = MPI_calculateVectorAvg(morph.getCorrelationData((char)1));
        correlation2_vect = MPI_calculateVectorAvg(morph.getCorrelationData((char)2));
    }
    // Prepare root processor for data collection.
    if(procid==0){
        // Create arrays to store property values gathered from each processor.
        mix_ratios = (double *)malloc(sizeof(double)*nproc);
        iav_ratios = (double *)malloc(sizeof(double)*nproc);
        iv_ratios = (double *)malloc(sizeof(double)*nproc);
        times = (double *)malloc(sizeof(double)*nproc);
        if(parameters.Enable_correlation_calc){
            domain_sizes1 = (double *)malloc(sizeof(double)*nproc);
            domain_sizes2 = (double *)malloc(sizeof(double)*nproc);
            domain_anisotropies1 = (double *)malloc(sizeof(double)*nproc);
            domain_anisotropies2 = (double *)malloc(sizeof(double)*nproc);
        }
        if(parameters.Enable_tortuosity_calc){
            island_ratios1 = (double *)malloc(sizeof(double)*nproc);
            island_ratios2 = (double *)malloc(sizeof(double)*nproc);
            // The final path data array will have a size that is the sum of the sizes of the arrays coming from each processor.
            for(int i=0;i<nproc;i++){
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
            for(int i=1;i<nproc;i++){
                pathdata1_displacement[i] = pathdata1_displacement[i-1] + pathdata1_sizes[i-1];
                pathdata2_displacement[i] = pathdata2_displacement[i-1] + pathdata2_sizes[i-1];
            }
        }
    }
    // Gather the properties from each processor into the previously created arrays on the root processor.
    MPI_Gather(&mix_ratio,1,MPI_DOUBLE,mix_ratios,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(&iav_ratio,1,MPI_DOUBLE,iav_ratios,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(&iv_ratio,1,MPI_DOUBLE,iv_ratios,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(&elapsedtime,1,MPI_DOUBLE,times,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(parameters.Enable_correlation_calc){
        MPI_Gather(&domain_size1,1,MPI_DOUBLE,domain_sizes1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(&domain_size2,1,MPI_DOUBLE,domain_sizes2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(&domain_anisotropy1,1,MPI_DOUBLE,domain_anisotropies1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(&domain_anisotropy2,1,MPI_DOUBLE,domain_anisotropies2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    // Gather the path data arrays from each processor into the previously created full path data array on the root processor at the positions defined by the displacement array.
    if(parameters.Enable_tortuosity_calc){
        MPI_Gatherv(pathdata1,pathdata1_size,MPI_FLOAT,pathdata1_all,pathdata1_sizes,pathdata1_displacement,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Gatherv(pathdata2,pathdata2_size,MPI_FLOAT,pathdata2_all,pathdata2_sizes,pathdata2_displacement,MPI_FLOAT,0,MPI_COMM_WORLD);
        // Gather the island volume fraction property from each processor into the previously created array on the root processor.
        MPI_Gather(&island_ratio1,1,MPI_DOUBLE,island_ratios1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(&island_ratio2,1,MPI_DOUBLE,island_ratios2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    // Output the analysis results to text files.
    if(procid==0){
        cout << "Writing morphology analysis data to files..." << endl;
        // Output the average pair-pair correlation function.
        if(parameters.Enable_correlation_calc){
            if(!Enable_import_morphology){
                ss << path << "correlation_data_avg.txt";
                correlationfile_avg.open(ss.str().c_str());
                ss.str("");
            }
            else{
                if(parameters.Enable_analysis_only){
                    ss << path << "correlation_data_avg_new.txt";
                    correlationfile_avg.open(ss.str().c_str());
                    ss.str("");
                }
                else{
                    ss << path << "correlation_data_avg_mod.txt";
                    correlationfile_avg.open(ss.str().c_str());
                    ss.str("");
                }
            }
            for(int i=0;i<(int)correlation1_vect.size();i++){
                correlationfile_avg << (double)i/2 << "," << correlation1_vect[i] << "," << correlation2_vect[i] << endl;
            }
            correlationfile_avg.close();
        }
        // Output the average tortuosity histograms and the end-to-end path data.
        if(parameters.Enable_tortuosity_calc){
            if(!Enable_import_morphology){
                ss << path << "tortuosity_histograms.txt";
                tortuosity_hist_file.open(ss.str().c_str());
                ss.str("");
                ss << path << "end-to-end_path_data1.txt";
                path_data1_file.open(ss.str().c_str());
                ss.str("");
                ss << path << "end-to-end_path_data2.txt";
                path_data2_file.open(ss.str().c_str());
                ss.str("");
            }
            else{
                if(parameters.Enable_analysis_only){
                    ss << path << "tortuosity_histograms_new.txt";
                    tortuosity_hist_file.open(ss.str().c_str());
                    ss.str("");
                    ss << path << "end-to-end_path_data1_new.txt";
                    path_data1_file.open(ss.str().c_str());
                    ss.str("");
                    ss << path << "end-to-end_path_data2_new.txt";
                    path_data2_file.open(ss.str().c_str());
                    ss.str("");
                }
                else{
                    ss << path << "tortuosity_histograms_mod.txt";
                    tortuosity_hist_file.open(ss.str().c_str());
                    ss.str("");
                    ss << path << "end-to-end_path_data1_mod.txt";
                    path_data1_file.open(ss.str().c_str());
                    ss.str("");
                    ss << path << "end-to-end_path_data2_mod.txt";
                    path_data2_file.open(ss.str().c_str());
                    ss.str("");
                }
            }
            tortuosity_hist_file << 0 << "," << tortuosity_hist1_vect[0] << "," << tortuosity_hist2_vect[0] << "\n";
            for(int i=1;i<(int)tortuosity_hist1_vect.size();i++){
                tortuosity_hist_file << (i+49.0)/50.0 << "," << tortuosity_hist1_vect[i] << "," << tortuosity_hist2_vect[i] << "\n";
            }
            for(int i=0;i<pathdata1_count;i++){
                path_data1_file << pathdata1_all[i] << "\n";
            }
            for(int i=0;i<pathdata2_count;i++){
                path_data2_file << pathdata2_all[i] << "\n";
            }
            tortuosity_hist_file.close();
            path_data1_file.close();
            path_data2_file.close();
        }
        // Output the interfacial distance histograms.
        if(parameters.Enable_interfacial_distance_calc){
            if(!Enable_import_morphology){
                ss << path << "interfacial_distance_histograms.txt";
                interfacial_dist_hist_file.open(ss.str().c_str());
                ss.str("");
            }
            else{
                if(parameters.Enable_analysis_only){
                    ss << path << "interfacial_distance_histograms_new.txt";
                    interfacial_dist_hist_file.open(ss.str().c_str());
                    ss.str("");
                }
                else{
                    ss << path << "interfacial_distance_histograms_mod.txt";
                    interfacial_dist_hist_file.open(ss.str().c_str());
                    ss.str("");
                }
            }
            for(int i=0;i<(int)interfacial_dist_hist1_vect.size();i++){
                interfacial_dist_hist_file << i+1 << "," << interfacial_dist_hist1_vect[i] << "," << interfacial_dist_hist2_vect[i] << "\n";
            }
            interfacial_dist_hist_file.close();
        }
        // Output the final morphology set analysis summary to a text file.
        if(!Enable_import_morphology){
            ss << path << "analysis_summary.txt";
            analysis_file.open(ss.str().c_str());
            ss.str("");
        }
        else{
            if(parameters.Enable_analysis_only){
                ss << path << "analysis_summary_new.txt";
                analysis_file.open(ss.str().c_str());
                ss.str("");
            }
            else{
                ss << path << "analysis_summary_mod.txt";
                analysis_file.open(ss.str().c_str());
                ss.str("");
            }
        }
        analysis_file << "Summary of results for this morphology set containing " << nproc <<" morphologies created using Ising_OPV " << version << ":" << endl;
        analysis_file << "length,width,height,mix_ratio_avg,mix_ratio_stdev,domain1_size_avg,domain1_size_stdev,domain2_size_avg,domain2_size_stdev,";
        analysis_file << "domain1_anisotropy_avg,domain1_anisotropy_stdev,domain2_anisotropy_avg,domain2_anisotropy_stdev,";
        analysis_file << "interfacial_area_volume_ratio_avg,interfacial_area_volume_ratio_stdev,interfacial_volume_ratio_avg,interfacial_volume_ratio_stdev,";
        analysis_file << "tortuosity1_avg,tortuosity1_stdev,tortuosity2_avg,tortuosity2_stdev,island_volume_ratio1_avg,island_volume_ratio1_stdev,";
        analysis_file << "island_volume_ratio2_avg,island_volume_ratio2_stdev,calc_time_avg(min),calc_time_stdev(min)" << endl;
        analysis_file << morph.getLength() << "," << morph.getWidth() << "," << morph.getHeight() << ",";
        analysis_file << array_avg(mix_ratios,nproc) << "," << array_stdev(mix_ratios,nproc) << ",";
        if(parameters.Enable_correlation_calc){
            analysis_file << array_avg(domain_sizes1,nproc) << "," << array_stdev(domain_sizes1,nproc) << "," << array_avg(domain_sizes2,nproc) << "," << array_stdev(domain_sizes2,nproc) << ",";
            analysis_file << array_avg(domain_anisotropies1,nproc) << "," << array_stdev(domain_anisotropies1,nproc) << "," << array_avg(domain_anisotropies2,nproc) << "," << array_stdev(domain_anisotropies2,nproc) << ",";
        }
        else{
            analysis_file << "-" << "," << "-" << "," << "-" << "," << "-" << ",";
            analysis_file << "-" << "," << "-" << "," << "-" << "," << "-" << ",";
        }
        analysis_file << array_avg(iav_ratios,nproc) << "," << array_stdev(iav_ratios,nproc) << "," << array_avg(iv_ratios,nproc) << "," << array_stdev(iv_ratios,nproc) << ",";
        if(parameters.Enable_tortuosity_calc) {
            analysis_file << array_avg(pathdata1_all,pathdata1_count) << "," << array_stdev(pathdata1_all,pathdata1_count) << ",";
            analysis_file << array_avg(pathdata2_all,pathdata2_count) << "," << array_stdev(pathdata2_all,pathdata2_count) << ",";
            analysis_file << array_avg(island_ratios1,nproc) << "," << array_stdev(island_ratios1,nproc) << ",";
            analysis_file << array_avg(island_ratios2,nproc) << "," << array_stdev(island_ratios2,nproc) << ",";
        }
        else{
            analysis_file << "-" << "," << "-" << ",";
            analysis_file << "-" << "," << "-" << ",";
            analysis_file << "-" << "," << "-" << ",";
            analysis_file << "-" << "," << "-" << ",";
        }
        analysis_file << array_avg(times,nproc) << "," << array_stdev(times,nproc) << endl;
        analysis_file << endl;
        analysis_file << "Detailed results for each of the morphologies in the set:" << endl;
        analysis_file << "id#,length,width,height,mix_ratio,domain1_size,domain2_size,";
        analysis_file << "domain1_anisotropy,domain2_anisotropy,";
        analysis_file << "interfacial_area_volume_ratio,interfacial_volume_ratio,";
        analysis_file << "tortuosity1,tortuosity2,island_volume_ratio1,island_volume_ratio2,calc_time(min)" << endl;
        double *tortuosities1 = (double *)malloc(sizeof(double)*nproc);
        double *tortuosities2 = (double *)malloc(sizeof(double)*nproc);
        for(int i=0;i<nproc;i++){
            analysis_file << i << "," << morph.getLength() << "," << morph.getWidth() << "," << morph.getHeight() << "," << mix_ratios[i] << ",";
            if(parameters.Enable_correlation_calc) {
                analysis_file << domain_sizes1[i] << "," << domain_sizes2[i] << ",";
                analysis_file << domain_anisotropies1[i] << "," << domain_anisotropies2[i] << ",";
            }
            else{
                analysis_file << "-" << "," << "-" << ",";
                analysis_file << "-" << "," << "-" << ",";
            }
            analysis_file << iav_ratios[i] << "," << iv_ratios[i] << ",";
            if(parameters.Enable_tortuosity_calc){
                pathdata1 = (float *)malloc(sizeof(float)*pathdata1_sizes[i]);
                pathdata2 = (float *)malloc(sizeof(float)*pathdata2_sizes[i]);
                copy(pathdata1_all+pathdata1_displacement[i],pathdata1_all+pathdata1_displacement[i]+pathdata1_sizes[i],pathdata1);
                copy(pathdata2_all+pathdata2_displacement[i],pathdata2_all+pathdata2_displacement[i]+pathdata2_sizes[i],pathdata2);
                tortuosities1[i] = array_avg(pathdata1,pathdata1_sizes[i]);
                tortuosities2[i] = array_avg(pathdata2,pathdata2_sizes[i]);
                if(isnan(tortuosities1[i])){
                    tortuosities1[i] = -1;
                }
                if(isnan(tortuosities2[i])){
                    tortuosities2[i] = -1;
                }
                analysis_file << tortuosities1[i] << "," << tortuosities2[i] << ",";
                analysis_file << island_ratios1[i] << "," << island_ratios2[i] << ",";
            }
            else{
                analysis_file << "-" << "," << "-" << ",";
                analysis_file << "-" << "," << "-" << ",";
            }
            analysis_file << times[i] << endl;
        }
        analysis_file << endl;
        analysis_file << "Morphology number " << array_which_median(domain_sizes1,nproc) << " has the median domain1 size of " << array_median(domain_sizes1,nproc) << endl;
        analysis_file << "Morphology number " << array_which_median(tortuosities1,nproc) << " has the median tortuosity1 of " << array_median(tortuosities1,nproc) << endl;
        analysis_file.close();
        cout << "Finished!" << endl;
    }
    MPI_Finalize();
    // The morphology object is automatically deconstructed upon return.
    return 0;
}

//  This function imports the parameters from an input parameter text file into the Input_Parameters data structure.
bool importParameters(ifstream * parameterfile,Input_Parameters& params){
    string line;
    string var;
    size_t pos;
	bool error_status = false;
    vector<string> stringvars;
    // Read input file line by line.
    while((*parameterfile).good()){
        getline(*parameterfile,line);
        // Skip lines designated as section breaks and section headers.
        if((line.substr(0,2)).compare("--")!=0 && (line.substr(0,2)).compare("##")!=0){
            // Strip off trailing comments from each line.
            pos = line.find("/",0);
            var = line.substr(0,pos-1);
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
    if(error_status){
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
    if(error_status){
        cout << "Error setting export options" << endl;
        return false;
    }
    i++;
    //enable_export_cross_section
    params.Enable_export_cross_section = importBooleanParam(stringvars[i], error_status);
	if(error_status){
        cout << "Error setting export cross-section options" << endl;
        return false;
    }
    i++;
    //enable_smoothing
    params.Enable_smoothing = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting morphology smoothing options" << endl;
        return false;
    }
    i++;
    params.Smoothing_threshold = atof(stringvars[i].c_str());
    i++;
    //enable_rescale
    params.Enable_rescale = importBooleanParam(stringvars[i], error_status);
    if(error_status){
        cout << "Error setting morphology rescale options" << endl;
        return false;
    }
    i++;
    params.Rescale_factor = atoi(stringvars[i].c_str());
    i++;
    params.Enable_shrink = importBooleanParam(stringvars[i], error_status);
	if(error_status){
        cout << "Error setting morphology shrink options" << endl;
        return false;
    }
    i++;
    //enable_interfacial_mixing
    params.Enable_interfacial_mixing = importBooleanParam(stringvars[i], error_status);
    if(error_status){
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
    if(error_status){
        cout << "Error setting analysis options" << endl;
        return false;
    }
    i++;
    //enable_correlation_calc
    params.Enable_correlation_calc = importBooleanParam(stringvars[i], error_status);
    if(error_status){
        cout << "Error setting correlation calculation options" << endl;
        return false;
    }
    i++;
    params.N_sampling_max = atoi(stringvars[i].c_str());
    i++;
    //enable_extended_correlation_calc
    params.Enable_extended_correlation_calc = importBooleanParam(stringvars[i], error_status);
    if(error_status){
        cout << "Error setting extended correlation calculation options" << endl;
        return false;
    }
    i++;
    //enable_interfacial_distance_calc
    params.Enable_interfacial_distance_calc = importBooleanParam(stringvars[i], error_status);
    if(error_status){
        cout << "Error setting interfacial distance calculation options" << endl;
        return false;
    }
    i++;
    //enable_tortuosity_calc
    params.Enable_tortuosity_calc = importBooleanParam(stringvars[i], error_status);
	if(error_status){
        cout << "Error setting tortuosity calculation options" << endl;
        return false;
    }
    i++;
    //enable_reduced_memory_tortuosity_calc
    params.Enable_reduced_memory_tortuosity_calc = importBooleanParam(stringvars[i], error_status);
    if(error_status){
        cout << "Error setting reduced memory tortuosity calculation options" << endl;
        return false;
    }
    i++;
    //enable_checkerboard_start
    params.Enable_checkerboard_start = importBooleanParam(stringvars[i], error_status);
    if(error_status){
        cout << "Error setting checkerboard starting condition" << endl;
        return false;
    }
    i++;
    //enable_growth_pref
    params.Enable_growth_pref = importBooleanParam(stringvars[i], error_status);
    if(error_status){
        cout << "Error setting growth preference conditions" << endl;
        return false;
    }
    i++;
    params.Growth_direction = atoi(stringvars[i].c_str());
    i++;
    params.Additional_interaction = atof(stringvars[i].c_str());
    i++;
    return true;
}

double array_median(const double data[],int array_size){
    vector<double> data_vect;
    data_vect.assign(array_size,0);
    for(int i=0;i<array_size;i++){
        data_vect[i] = data[i];
    }
    sort(data_vect.begin(),data_vect.end());
    return data_vect[array_size/2];
}

int array_which_median(const double data[],int array_size){
    vector<double> data_vect;
    data_vect.assign(array_size,0);
    for(int i=0;i<array_size;i++){
        data_vect[i] = data[i];
    }
    sort(data_vect.begin(),data_vect.end());
    double median = data_vect[array_size/2];
    for(int i=0;i<array_size;i++){
        if(data[i]==median){
            return i;
        }
    }
    cout << "Error! Median not found." << endl;
    return -1;
}
