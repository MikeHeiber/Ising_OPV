// Copyright (c) 2015 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.

#include <iostream>
#include <fstream>
#include <mpi.h>
#include "Morphology.h"

using namespace std;

struct Input_Parameters{
    // General
    int length; // x-direction size of the lattice
    int width; // y-direction size of the lattice
    int height; // z-direction size of the lattice
    bool enable_z_periodic_boundary; // z-direction periodic boundary option
    double mix_fraction; // volume fraction of donor
    double interaction_energy1; // energetic favorability for type1-type1 interactions over type1-type2 interactions
    double interaction_energy2; // energetic favorability for type2-type2 interactions over type1-type2 interactions
    int MC_steps; // number of MC steps to be executed (determines number of Ising swapping iterations)
    // Export Morphology Options
    bool enable_export_compressed_files; // choose whether the output morphology data file is in compressed format or not
    bool enable_export_cross_section; // choose whether to output data for an uncompressed cross section (x=0 plane) of the morphology
    // Smoothing Options
    bool enable_smoothing; // choose whether to perform domain smoothing
    double smoothing_threshold;  // specify the degree of smoothing
    // Rescale Options
    bool enable_rescale; // choose whether to perform lattice rescaling
    int rescale_factor; // specify the rescale factor to be used (must be an integer greater than 1)
    // Interfacial Mixing Options
    bool enable_interfacial_mixing; // choose whether to perform interfacial mixing
    double interface_width;  // specify the interfacial width
    double interface_conc; // specify the mixing concentration in the interfacial region
    // Analysis Options
    bool enable_analysis_only;  // choose whether to only perform analysis on an imported morphology
    bool enable_correlation_calc; // chose whether to perform the domain size calculation using the pair-pair correlation method
    int N_sampling_max; // specify the maximum number of sites to be sampled for calculating the correlation function
    bool enable_interfacial_distance_calc; // chose whether to calculate the interfacial distance histograms
    bool enable_tortuosity_calc; // chose whether to calculate the tortuosity histograms, end-to-end path distances, and island volume fraction
};

double array_avg(const double data[],int size);
double array_stdev(const double data[],int size);
float array_avg(const float data[],int size);
float array_stdev(const float data[],int size);
vector<double> calculateAverageVector(vector<double> input_vector,int procid,int nproc);
bool importParameters(ifstream * parameterfile,Input_Parameters * parameters);

int main(int argc, char * argv[]){
    // Input parameters
    Input_Parameters parameters;
    // Internal parameters
    bool enable_import_morphology;
    double mix_ratio = 0;
    double domain_size1 = 0;
    double domain_size2 = 0;
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
    ofstream interfacial_dist_hist_file;
    ofstream tortuosity_hist_file;
    ofstream path_data1_file;
    ofstream path_data2_file;
    bool success;
    int cutoff_distance;
    double *mix_ratios = NULL;
    double *domain_sizes1 = NULL;
    double *domain_sizes2 = NULL;
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
    string path = "";
    string input_morphology, input_file_path, compression_str, filename_prefix;
    bool is_file_compressed;
    // Begin
    start_time = time(NULL);
    // Initialize parallel processing.
    MPI::Init(argc,argv);
    nproc = MPI::COMM_WORLD.Get_size();
    procid = MPI::COMM_WORLD.Get_rank();
    // Import parameters from text file.
    if(argc == 2){
        cout << procid << ": Parameter file found." << endl;
        enable_import_morphology=false;
    }
    else if(argc == 3){
        cout << procid << ": Parameter file and morphology file found." << endl;
        enable_import_morphology=true;
    }
    else{
        cout << procid << ": Incorrect arguments. There needs to be at least one argument which is the filename of a parameter file. Second optional argument needs to be a morphology file of the form 'filename_procid_compressed.txt' or 'filename_procid_uncompressed.txt'. Program will exit now!" << endl;
        return 0;
    }
    parameter_file.open((argv[1]),ifstream::in);
    success = importParameters(&parameter_file,&parameters);
    parameter_file.close();
    if(!success){
        cout << procid << ": Error importing variables from file.  Program will now exit." << endl;
        return 0;
    }
    if(parameters.enable_analysis_only && !enable_import_morphology){
        cout << procid << ": Error!  The 'analysis only' option can only be used when a morphology is imported from a file." << endl;
        return 0;
    }
    if(parameters.enable_analysis_only){
        cout << procid << ": Warning! Only morphology analysis will be performed." << endl;
    }
    // Create morphology data structure.
    Morphology morph(parameters.length,parameters.width,parameters.height,parameters.enable_z_periodic_boundary,procid,time(0)+procid);
    // Import morphology if enabled.
    if(enable_import_morphology){
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
            cout << procid << ": Error! Format of the file needs to be 'filename_#_compressed.txt' or 'filename_#_uncompressed.txt'. Program will exit now!" << endl;
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
            cout << procid << ": Error! Format of the file needs to be 'filename_#_compressed.txt' or 'filename_#_uncompressed.txt'. Program will exit now!" << endl;
            return 0;
        }
        // Separate filename prefix from ID number.
        pos_id = filename_prefix.find_last_of("_");
        if(pos_id==string::npos){
            cout << procid << ": Error! Format of the file needs to be 'filename_#_compressed.txt' or 'filename_#_uncompressed.txt'. Program will exit now!" << endl;
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
        cout << "complete!" << endl;
    }
    // Create random morphology if import is disabled.
    else{
        cout << procid << ": Generating initial random morphology..." << endl;
        morph.createRandomMorphology(parameters.mix_fraction);
    }
    // Determine if any phase separation is to be executed on the morphology.
    if(parameters.MC_steps>0  && !parameters.enable_analysis_only){
        N_steps = parameters.MC_steps;
    }
    else{
        N_steps = 0;
    }
    // Output the number of MC steps to be executed.
    if(procid==0){
        cout << N_steps << " MC steps of site swapping will be executed." << endl;
    }
    // Execute phase separation through Ising swapping.
    if(N_steps>0){
        cout << procid << ": Executing site swapping..." << endl;
        morph.executeIsingSwapping(N_steps,parameters.interaction_energy1,parameters.interaction_energy2);
    }
    // Perform lattice rescale if enabled.
    if(parameters.enable_rescale && !parameters.enable_analysis_only){
        cout << procid << ": Rescaling lattice by a factor of " << parameters.rescale_factor << " ..." << endl;
        morph.rescaleLattice(parameters.rescale_factor);
        parameters.length = parameters.length*parameters.rescale_factor;
        parameters.width = parameters.width*parameters.rescale_factor;
        parameters.height = parameters.height*parameters.rescale_factor;
    }
    // Perform domain smoothing if enabled.
    if(parameters.enable_smoothing && !parameters.enable_analysis_only){
        if(!parameters.enable_rescale){
            cout << procid << ": Executing standard smoothing..." << endl;
            morph.executeSmoothing(parameters.smoothing_threshold,2);
        }
        else{
            cout << procid << ": Executing rescale factor dependent smoothing..." << endl;
            morph.executeSmoothing(parameters.smoothing_threshold,parameters.rescale_factor);
        }
    }
    // Perform interfacial mixing if enabled.
    if(parameters.enable_interfacial_mixing && !parameters.enable_analysis_only){
        cout << "Executing interfacial mixing..." << endl;
        morph.executeMixing(parameters.interface_width,parameters.interface_conc);
    }
    // Calculate domain size if enabled.
    if(parameters.enable_correlation_calc){
        cutoff_distance = 5;
        cout << procid << ": Calculating the domain size from the pair-pair correlation function..." << endl;
        success = false;
        // The correlation function calculation is called with an increasing cutoff distance until successful.
        while(!success){
            if(cutoff_distance<morph.getLength() || cutoff_distance<morph.getWidth()){
                success = morph.calculateCorrelationDistance(cutoff_distance,parameters.N_sampling_max);
            }
            else{
                break;
            }
            cutoff_distance++;
        }
        if(success){
            domain_size1 = morph.getDomainSize(1);
            domain_size2 = morph.getDomainSize(2);
        }
    }
    // Calculate interfacial distance histogram if enabled.
    if(parameters.enable_interfacial_distance_calc){
        cout << procid << ": Calculating the interfacial distance histogram..." << endl;
        morph.calculateInterfacialDistance();
    }
    // Calculate interfacial area to volume ratio.
    iav_ratio = morph.calculateInterfacialArea()/(morph.getLength()*morph.getWidth()*morph.getHeight());
    // Calculate interfacial volume to total volume ratio.
    iv_ratio = morph.calculateInterfacialVolume()/(morph.getLength()*morph.getWidth()*morph.getHeight());
    // Get Final Mix ratio
    mix_ratio = morph.getMixFraction();
    // Calculate end-to-end distances, tortuosity histogram, and island volume fraction.
    if(parameters.enable_tortuosity_calc){
        cout << procid << ": Calculating tortuosity..." << endl;
        success = morph.calculateTortuosity();
        // Calculate island volume ratio.
        island_ratio1 = (double)morph.getIslandVolume(1)/(morph.getLength()*morph.getWidth()*morph.getHeight());
        island_ratio2 = (double)morph.getIslandVolume(2)/(morph.getLength()*morph.getWidth()*morph.getHeight());
    }
    // Save final morphology to a text file.
    if(!parameters.enable_analysis_only){
        cout << procid << ": Writing morphology to file..." << endl;
        if(!enable_import_morphology){
            if(!parameters.enable_export_compressed_files){
                ss << path << "morphology_" << procid << "_uncompressed.txt";
                morphology_output_file.open(ss.str().c_str());
                ss.str("");
            }
            else{
                ss << path << "morphology_" << procid << "_compressed.txt";
                morphology_output_file.open(ss.str().c_str());
                ss.str("");
            }
        }
        else{
            if(!parameters.enable_export_compressed_files){
                ss << path << filename_prefix << "mod_" << procid << "_uncompressed.txt";
                morphology_output_file.open(ss.str().c_str());
                ss.str("");
            }
            else{
                ss << path << filename_prefix << "mod_" << procid << "_compressed.txt";
                morphology_output_file.open(ss.str().c_str());
                ss.str("");
            }
        }
        morph.outputMorphologyFile(&morphology_output_file,parameters.enable_export_compressed_files);
        morphology_output_file.close();
    }
    // Save the cross-section of the x=0 plane to a file if enabled.
    if(parameters.enable_export_cross_section){
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
    if(parameters.enable_tortuosity_calc){
        // Get path data for end-to-end path distance calculation.
        vector<float> data = morph.getTortuosityData(1);
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
        data = morph.getTortuosityData(2);
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
        tortuosity_hist1_vect = calculateAverageVector(morph.getTortuosityHistogram(1),procid,nproc);
        tortuosity_hist2_vect = calculateAverageVector(morph.getTortuosityHistogram(2),procid,nproc);
    }
    // Calculate the average interfacial distance histograms.
    if(parameters.enable_interfacial_distance_calc){
        interfacial_dist_hist1_vect = calculateAverageVector(morph.getInterfacialHistogram(1),procid,nproc);
        interfacial_dist_hist2_vect = calculateAverageVector(morph.getInterfacialHistogram(2),procid,nproc);
    }
    // Prepare root processor for data collection.
    if(procid==0){
        // Create arrays to store property values gathered from each processor.
        mix_ratios = (double *)malloc(sizeof(double)*nproc);
        iav_ratios = (double *)malloc(sizeof(double)*nproc);
        iv_ratios = (double *)malloc(sizeof(double)*nproc);
        times = (double *)malloc(sizeof(double)*nproc);
        if(parameters.enable_correlation_calc){
            domain_sizes1 = (double *)malloc(sizeof(double)*nproc);
            domain_sizes2 = (double *)malloc(sizeof(double)*nproc);
        }
        if(parameters.enable_tortuosity_calc){
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
    if(parameters.enable_correlation_calc){
        MPI_Gather(&domain_size1,1,MPI_DOUBLE,domain_sizes1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(&domain_size2,1,MPI_DOUBLE,domain_sizes2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    // Gather the path data arrays from each processor into the previously created full path data array on the root processor at the positions defined by the displacement array.
    if(parameters.enable_tortuosity_calc){
        MPI_Gatherv(pathdata1,pathdata1_size,MPI_FLOAT,pathdata1_all,pathdata1_sizes,pathdata1_displacement,MPI_FLOAT,0,MPI_COMM_WORLD);
        MPI_Gatherv(pathdata2,pathdata2_size,MPI_FLOAT,pathdata2_all,pathdata2_sizes,pathdata2_displacement,MPI_FLOAT,0,MPI_COMM_WORLD);
        // Gather the island volume fraction property from each processor into the previously created array on the root processor.
        MPI_Gather(&island_ratio1,1,MPI_DOUBLE,island_ratios1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Gather(&island_ratio2,1,MPI_DOUBLE,island_ratios2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    // Output the analysis results to text files.
    if(procid==0){
        cout << "Writing morphology analysis data to files..." << endl;
        // Output the average tortuosity histograms and the end-to-end path data.
        if(parameters.enable_tortuosity_calc){
            if(!enable_import_morphology){
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
                if(parameters.enable_analysis_only){
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
        if(parameters.enable_interfacial_distance_calc){
            if(!enable_import_morphology){
                ss << path << "interfacial_distance_histograms.txt";
                interfacial_dist_hist_file.open(ss.str().c_str());
                ss.str("");
            }
            else{
                if(parameters.enable_analysis_only){
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
        if(!enable_import_morphology){
            ss << path << "analysis_summary.txt";
            analysis_file.open(ss.str().c_str());
            ss.str("");
        }
        else{
            if(parameters.enable_analysis_only){
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
        analysis_file << "length,width,height,mix_ratio_avg,mix_ratio_stdev,domain1_size_avg,domain1_size_stdev,domain2_size_avg,domain2_size_stdev,";
        analysis_file << "interfacial_area_volume_ratio_avg,interfacial_area_volume_ratio_stdev,interfacial_volume_ratio_avg,interfacial_volume_ratio_stdev,";
        analysis_file << "tortuosity1_avg,tortuosity1_stdev,tortuosity2_avg,tortuosity2_stdev,island_volume_ratio1_avg,island_volume_ratio1_stdev,";
        analysis_file << "island_volume_ratio2_avg,island_volume_ratio2_stdev,calc_time_avg(min),calc_time_stdev(min)" << endl;
        analysis_file << morph.getLength() << "," << morph.getWidth() << "," << morph.getHeight() << ",";
        analysis_file << array_avg(mix_ratios,nproc) << "," << array_stdev(mix_ratios,nproc) << ",";
        if(parameters.enable_correlation_calc) {
            analysis_file << array_avg(domain_sizes1,nproc) << "," << array_stdev(domain_sizes1,nproc) << "," << array_avg(domain_sizes2,nproc) << "," << array_stdev(domain_sizes2,nproc) << ",";
        }
        else{
            analysis_file << "-" << "," << "-" << "," << "-" << "," << "-" << ",";
        }
        analysis_file << array_avg(iav_ratios,nproc) << "," << array_stdev(iav_ratios,nproc) << "," << array_avg(iv_ratios,nproc) << "," << array_stdev(iv_ratios,nproc) << ",";
        if(parameters.enable_tortuosity_calc) {
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
        analysis_file.close();
        cout << "Finished!" << endl;
    }
    MPI_Finalize();
    // The morphology object is automatically deconstructed upon return.
    return 0;
}

//  This function imports the parameters from an input parameter text file into the Input_Parameters data structure.
bool importParameters(ifstream * parameterfile,Input_Parameters * parameters){
    string line;
    string var;
    size_t pos;
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
    (*parameters).length = atoi(stringvars[i].c_str());
    i++;
    (*parameters).width = atoi(stringvars[i].c_str());
    i++;
    (*parameters).height = atoi(stringvars[i].c_str());
    i++;
    //enable_z_periodic_boundary
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_z_periodic_boundary = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_z_periodic_boundary = false;
    }
    else{
        cout << "Error setting z-boundary options" << endl;
        return false;
    }
    i++;
    (*parameters).mix_fraction = atof(stringvars[i].c_str());
    i++;
    (*parameters).interaction_energy1 = atof(stringvars[i].c_str());
    i++;
    (*parameters).interaction_energy2 = atof(stringvars[i].c_str());
    i++;
    (*parameters).MC_steps = atoi(stringvars[i].c_str());
    i++;
    //enable_export_compressed_files
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_export_compressed_files = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_export_compressed_files = false;
    }
    else{
        cout << "Error setting export options" << endl;
        return false;
    }
    i++;
    //enable_export_cross_section
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_export_cross_section = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_export_cross_section = false;
    }
    else{
        cout << "Error setting export cross-section options" << endl;
        return false;
    }
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
    //enable_interfacial_mixing
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_interfacial_mixing = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_interfacial_mixing = false;
    }
    else{
        cout << "Error setting interfacial mixing options" << endl;
        return false;
    }
    i++;
    (*parameters).interface_width = atof(stringvars[i].c_str());
    i++;
    (*parameters).interface_conc = atof(stringvars[i].c_str());
    i++;
    //enable_analysis_only
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_analysis_only = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_analysis_only = false;
    }
    else{
        cout << "Error setting analysis options" << endl;
        return false;
    }
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
    (*parameters).N_sampling_max = atoi(stringvars[i].c_str());
    i++;
    //enable_interfacial_distance_calc
    if(stringvars[i].compare("true")==0){
        (*parameters).enable_interfacial_distance_calc = true;
    }
    else if(stringvars[i].compare("false")==0){
        (*parameters).enable_interfacial_distance_calc = false;
    }
    else{
        cout << "Error setting interfacial distance calculation options" << endl;
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
    return true;
}

//  This function calculates the average value of an array of doubles.
double array_avg(const double data[],int array_size){
    double sum = 0;
    for(int i=0;i<array_size;i++){
        sum += data[i];
    }
    return sum/array_size;
}

//  This function calculates the standard deviation of an array of doubles.
double array_stdev(const double data[],int array_size){
    double sum = 0;
    double avg = array_avg(data,array_size);
    for(int i=0;i<array_size;i++){
        sum += pow(data[i]-avg,2);
    }
    return sqrt(sum/(array_size-1));
}

//  This function calculates the average value of an array of floats.
float array_avg(const float data[],int array_size){
    float sum = 0;
    for(int i=0;i<array_size;i++){
        sum += data[i];
    }
    return sum/array_size;
}

//  This function calculates the standard deviation of an array of floats.
float array_stdev(const float data[],int array_size){
    float sum = 0;
    float avg = array_avg(data,array_size);
    for(int i=0;i<array_size;i++){
        sum += pow(data[i]-avg,2);
    }
    return sqrt(sum/(array_size-1));
}

// This function is uses the MPI Gather function to calculate the average vector from a set of vectors where one comes from each processor.
// The vectors do not need to be the same size, but vectors are assumed to have zero values beyond their length.
vector<double> calculateAverageVector(vector<double> input_vector,int procid,int nproc){
    int data_size = 0;
    int data_count = 0;
    double *data = NULL;
    double *data_all = NULL;
    int *data_sizes = NULL;
    int *data_displacement = NULL;
    int max_data_size = 0;
    double average = 0;
    vector<double> output_vector;
    if(procid==0){
        data_sizes = (int *)malloc(sizeof(int)*nproc);
    }
    data_size = (int)input_vector.size();
    data = (double *)malloc(sizeof(double)*data_size);
    for(int i=0;i<(int)input_vector.size();i++){
        data[i] = input_vector[i];
    }
    MPI_Gather(&data_size,1,MPI_INT,data_sizes,1,MPI_INT,0,MPI_COMM_WORLD);
    if(procid==0){
        for(int i=0;i<nproc;i++){
            data_count += data_sizes[i];
        }
        data_all = (double *)malloc(sizeof(double)*data_count);
        data_displacement = (int *)malloc(sizeof(int)*nproc);
        data_displacement[0] = 0;
        for(int i=1;i<nproc;i++){
            data_displacement[i] = data_displacement[i-1] + data_sizes[i-1];
        }
    }
    MPI_Gatherv(data,data_size,MPI_DOUBLE,data_all,data_sizes,data_displacement,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(procid==0){
        for(int i=0;i<nproc;i++){
            if(data_sizes[i]>max_data_size){
                max_data_size =data_sizes[i];
            }
        }
        for(int i=0;i<max_data_size;i++){
            average = 0;
            for(int j=0;j<nproc;j++){
                if(i<data_sizes[j]){
                    average += data_all[data_displacement[j]+i];
                }
            }
            average = average/nproc;
            output_vector.push_back(average);
        }
    }
    return output_vector;
}
