// Copyright (c) 2015 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this package.

#include "Morphology.h"

//  This constructor creates a Morphology object including a 3D lattice with a size defined by the input dimensions (length, width, height).
//  Two-dimensional periodic boundaries in the x- and y- directions are implemented by default, but periodic boundaries in the z-direction can also be enabled upon construction.
//  Morphology objects are also tagged with an integer identification number.
//  Each morphology object has a random number generator that is seeded by the seed input parameter upon creation.
Morphology::Morphology(int length, int width, int height, bool enable_z_periodic_boundary, int id, int seed){
    ID = id;
    Length = length;
    Width = width;
    Height = height;
    Mix_fraction = 0;
    Enable_z_periodic_boundary = enable_z_periodic_boundary;
    Enable_third_neighbor_interaction = false;
    Site site;
    site.type = (char)0;
    lattice.assign(length*width*height,site);
    Domain_size1_updated = false;
    Domain_size2_updated = false;
    Domain_size1 = 0;
    Domain_size2 = 0;
    Island_volume1 = 0;
    Island_volume1 = 0;
    gen.seed(static_cast<boost::uint32_t>(seed));
}

//  Default deconstructor
Morphology::~Morphology(){
    //dtor
}

//  This function calculates the domain size of the morphology based on the pair-pair correlation function
//  The correlation function is calculated from each starting site out to the cutoff distance.
//  The domain size is defined as the distance at which the pair-pair correlation function first crosses the value equal to the mixing fraction
//  If this cross-over point is not reach within the cutoff distance, the function returns false.
//  For large lattices, the correlation function does not need to be calculated starting from every site to collect enough statistics and instead a sampling of starting sites can be used.
//  When the total number of sites is greater than N_sampling_max, N_sampling_max sites are randomly selected and saved for performing a correlation function calculation by sampling.
//  When the total number of sites is less than N_sampling_max, all sites will be used as starting points for the correlation function calculation.
//  If the function returns false and the function is re-called with a larger cutoff_distance, the correlation function is not recalculated for close distances and only fills in the missing data for larger distances.
bool Morphology::calculateCorrelationDistance(int cutoff_distance,int N_sampling_max){
    int x,y,z;
    int dx,dy,dz;
    vector<int> site_count,total_count;
    double distance;
    int bin;
    int N1 = 0;
    int N2 = 0;
    double d1,y1,y2,slope,intercept,diff1,diff2;
    bool success;
    Coords site_coords;
    if(cutoff_distance>Length || cutoff_distance>Width){
        cout << ID << ": Error, cutoff distance is greater than the lattice length and/or width." << endl;
        return false;
    }
    if(Domain_size1_updated && Domain_size2_updated){
        Domain_size1_updated = false;
        Domain_size2_updated = false;
    }
    // Resolution of correlation distance data is 0.5 lattice units
    int correlation_size_old = Correlation1.size();
    int correlation_size_new = 2*cutoff_distance+1;
    if(correlation_size_old>=correlation_size_new){
        cout << ID << ": Error, new cutoff distance is not greater than the previous cutoff distance and no new calculations have been performed." << endl;
        return false;
    }
    // Initialize vectors to store correlation function data
    for(int m=0;m<(correlation_size_new-correlation_size_old);m++){
        Correlation1.push_back(0);
        Correlation2.push_back(0);
    }
    if((Length*Width*Height)<N_sampling_max){
        cout << ID << ": Performing complete domain size calculation with a cutoff of " << cutoff_distance << "." << endl;
    }
    else{
        cout << ID << ": Performing sampling domain size calculation with a cutoff of " << cutoff_distance << "." << endl;
    }
    // Select sites for correlation function calculation.
    // Site indices for each selected site are stored in the Correlation_sites vector.
    if((int)Correlation_sites.size()==0){
        if((Length*Width*Height)<N_sampling_max){
            // All sites in the lattice are selected
            Correlation_sites.assign(Length*Width*Height,-1);
            for(int n=0;n<Length*Width*Height;n++){
                Correlation_sites[n] = n;
            }
        }
        else{
            // Only N_sampling_max sites are randomly selected
            Correlation_sites.assign(N_sampling_max,-1);
            boost::uniform_int<> dist_x(0,Length-1);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_x(gen,dist_x);
            boost::uniform_int<> dist_y(0,Width-1);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_y(gen,dist_y);
            boost::uniform_int<> dist_z(0,Height-1);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_z(gen,dist_z);
            unsigned long site_index;
            int j;
            for(int n=0;n<N_sampling_max;n++){
                success=false;
                while(!success){
                    // Randomly select a site in the lattice
                    x = rand_x();
                    y = rand_y();
                    z = rand_z();
                    site_index=getSite(x,y,z);
                    j=0;
                    success=true;
                    // Make sure the randomly selected site is not a duplicate of one previously selected
                    while(Correlation_sites[j]>0){
                        if(site_index==Correlation_sites[j]){
                            success=false;
                            break;
                        }
                        j++;
                    }
                }
                Correlation_sites[n] = site_index;
            }
        }
    }
    // Loop through all selected sites and determine the correlation function for each
    // The pair-pair correlation is determined based on the fraction of sites that are the same as the starting site and this function is calculated as a function of distance from the starting site.
    // Sites surrounding the start site are placed into bins based on their distance from the starting site.
    // Bins covering a distance range of half a lattice unit are used, ex: second bin is from 0.25a to 0.7499a, third bin is from 0.75a to 1.2499a, etc.
    // site_count vector stores the number of sites that are are the same type as the starting site for each bin
    // total_count vector stores the total number of sites in each bin
    site_count.assign(2*cutoff_distance+1,0);
    total_count.assign(2*cutoff_distance+1,0);
    for(int m=0;m<(int)Correlation_sites.size();m++){
        site_coords = getCoords(Correlation_sites[m]);
        x = site_coords.x;
        y = site_coords.y;
        z = site_coords.z;
        // Count the number of sites of each type
        if(lattice[getSite(x,y,z)].type==(char)1){
            if(Domain_size1_updated){
                continue;
            }
            N1++;
        }
        else{
            if(Domain_size2_updated){
                continue;
            }
            N2++;
        }
        fill(site_count.begin(),site_count.end(),0);
        fill(total_count.begin(),total_count.end(),0);
        for(int i=-cutoff_distance;i<=cutoff_distance;i++){
            for(int j=-cutoff_distance;j<=cutoff_distance;j++){
                for(int k=-cutoff_distance;k<=cutoff_distance;k++){
                    if(i==0 && j==0 && k==0){
                        continue;
                    }
                    // The distance between two sites is rounded to the nearest half a lattice unit
                    bin = round_int(2*sqrt(i*i+j*j+k*k));
                    // Calculation is skipped for bin values that have already been calculated during previous calls to the calculateCorrelationDistance function
                    if(bin<(correlation_size_old-1)){
                        continue;
                    }
                    distance = (double)bin/2;
                    if(distance>cutoff_distance){
                        continue;
                    }
                    if(!Enable_z_periodic_boundary){
                        if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                            continue;
                        }
                        dz=0;
                    }
                    else{
                        if(z+k<0){ // Check for periodic z boundary
                            dz=Height;
                        }
                        else if(z+k>=Height){
                            dz=-Height;
                        }
                        else{
                            dz=0;
                        }
                    }
                    if(x+i<0){ // Check for x boundary
                        dx=Length;
                    }
                    else if(x+i>=Length){
                        dx=-Length;
                    }
                    else{
                        dx=0;
                    }
                    if(y+j<0){ // Check for y boundary
                        dy=Width;
                    }
                    else if(y+j>=Width){
                        dy=-Width;
                    }
                    else{
                        dy=0;
                    }
                    if(lattice[getSite(x,y,z)].type==lattice[getSite(x+i+dx,y+j+dy,z+k+dz)].type){
                        site_count[bin]++;
                    }
                    total_count[bin]++;
                }
            }
        }
        //  Calculate the fraction of similar sites for each bin
        for(int n=0;n<=2*cutoff_distance;n++){
            if(n<(correlation_size_old-1)){
                continue;
            }
            if(lattice[getSite(x,y,z)].type==(char)1){
                if(total_count[n]>0){
                    Correlation1[n] += (double)site_count[n]/(double)total_count[n];
                }
                else{
                    Correlation1[n] += 1;
                }
            }
            else{
                if(total_count[n]>0){
                    Correlation2[n] += (double)site_count[n]/(double)total_count[n];
                }
                else{
                    Correlation2[n] += 1;
                }
            }
        }
    }
    // Overall correlation function is an average of contributions from each starting site of the corresponding type
    for(int n=0;n<=2*cutoff_distance;n++){
        if(n<(correlation_size_old-1)){
            continue;
        }
        Correlation1[n] = Correlation1[n]/N1;
        Correlation2[n] = Correlation2[n]/N2;
    }
    // Determine if calculation is successful
    if(!Domain_size1_updated){
        // Find the bounds of where the pair-pair correlation1 function first crosses over the Mix_fraction
        success = false;
        for(int n=2;n<(int)Correlation1.size();n++){
            if(Correlation1[n]<Mix_fraction){
                d1 = ((double)n-1)/2;
                y1 = Correlation1[n-1];
                y2 = Correlation1[n];
                // Use linear interpolation to determine the cross-over point
                slope = (y2-y1)/0.5;
                intercept = y1-slope*d1;
                Domain_size1 = (Mix_fraction-intercept)/slope;
                success = true;
                break;
            }
        }
        if(!success){
            cout << ID << ": Cutoff distance of " << cutoff_distance << " is too small to calculate size of domain type 1." << endl;
            Domain_size1 = -1;
        }
    }
    if(!Domain_size2_updated){
        // Find the bounds of where the pair-pair correlation1 function first crosses over the Mix_fraction
        success = false;
        for(int n=2;n<(int)Correlation2.size();n++){
            if(Correlation2[n]<(1-Mix_fraction)){
                d1 = ((double)n-1)/2;
                y1 = Correlation2[n-1];
                y2 = Correlation2[n];
                // Use linear interpolation to determine the cross-over point
                slope = (y2-y1)/0.5;
                intercept = y1-slope*d1;
                Domain_size2 = ((1-Mix_fraction)-intercept)/slope;
                success = true;
                break;
            }
        }
        if(!success){
            cout << ID << ": Cutoff distance of " << cutoff_distance << " is too small to calculate size of domain type 2." << endl;
            Domain_size2 = -1;
        }
    }
    if(Domain_size1<0 || Domain_size2<0){
        return false;
    }
    return true;
}

//  Calculates the change in energy of the system that would occur if the adjacent sites at (x1,y1,z1) and (x2,y2,z2) were to be swapped
//  Sites must be adjacent to each other for calculation to be correct. (Works for adjacent sites across periodic boundaries)
//  When non-periodic/hard z-boundaries are used, it is assumed that neither site type has a preferential interaction with the z-boundary
double Morphology::calculateEnergyChange(int x1,int y1,int z1, int x2, int y2, int z2,double interaction_energy1,double interaction_energy2){
    // Used with bond formation algorithm
    int dx,dy,dz;
    int site1_type,site2_type;
    int sum1_1_delta, sum2_1_delta, sum3_1_delta, sum1_2_delta, sum2_2_delta, sum3_2_delta;
    double sum_1_delta,sum_2_delta;
    static const double one_over_sqrt2 = 1/sqrt(2);
    static const double one_over_sqrt3 = 1/sqrt(3);
    int sum1_1i = 0;
    int sum2_1i = 0;
    int sum3_1i = 0;
    int sum1_2i = 0;
    int sum2_2i = 0;
    int sum3_2i = 0;
    int sum1_1f = 0;
    int sum2_1f = 0;
    int sum3_1f = 0;
    int sum1_2f = 0;
    int sum2_2f = 0;
    int sum3_2f = 0;
    // There are in total 6 first-nearest, 12 second-nearest, and 8 third-nearest neighbors
    int total1 = 6;
    int total2 = 12;
    int total3 = 8;
    // Calculate change around x1,y1,z1
    site1_type = lattice[getSite(x1,y1,z1)].type;
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
            for(int k=-1;k<=1;k++){
                if(i==0 && j==0 && k==0){
                    continue;
                }
                if(!Enable_z_periodic_boundary){
                    if(z1+k>=Height || z1+k<0 ){ // Check for z boundary
                        // Total site counts must be reduced if next to a hard boundary
                        switch(i*i+j*j+k*k){
                            case 1:
                                total1--;
                                break;
                            case 2:
                                total2--;
                                break;
                            case 3:
                                total3--;
                                break;
                            default:
                                break;
                        }
                        continue;
                    }
                    dz = 0;
                }
                else{
                    if(z1+k<0){ // Check for periodic z boundary
                        dz=Height;
                    }
                    else if(z1+k>=Height){
                        dz=-Height;
                    }
                    else{
                        dz = 0;
                    }
                }
                if(x1+i<0){ // Check for x boundary
                    dx=Length;
                }
                else if(x1+i>=Length){
                    dx=-Length;
                }
                else{
                    dx = 0;
                }
                if(y1+j<0){ // Check for y boundary
                    dy=Width;
                }
                else if(y1+j>=Width){
                    dy=-Width;
                }
                else{
                    dy = 0;
                }
                // Count the number of similar neighbors
                if(lattice[getSite(x1+i+dx,y1+j+dy,z1+k+dz)].type==site1_type){
                    switch(i*i+j*j+k*k){
                        case 1:
                            sum1_1i++;
                            break;
                        case 2:
                            sum2_1i++;
                            break;
                        case 3:
                            sum3_1i++;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
    }
    sum1_2f = total1-sum1_1i-1;
    sum2_2f = total2-sum2_1i;
    sum3_2f = total3-sum3_1i;
    // Calculate change around x2,y2,z2
    site2_type = lattice[getSite(x2,y2,z2)].type;
    // There are in total 6 first-nearest, 12 second-nearest, and 8 third-nearest neighbors
    total1 = 6;
    total2 = 12;
    total3 = 8;
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
            for(int k=-1;k<=1;k++){
                if(i==0 && j==0 && k==0){
                    continue;
                }
                if(!Enable_z_periodic_boundary){
                    if(z2+k>=Height || z2+k<0 ){ // Check for hard z boundary
                        // Total site counts must be reduced if next to a hard boundary
                        switch(i*i+j*j+k*k){
                            case 1:
                                total1--;
                                break;
                            case 2:
                                total2--;
                                break;
                            case 3:
                                total3--;
                                break;
                            default:
                                break;
                        }
                        continue;
                    }
                    dz = 0;
                }
                else{
                    if(z2+k<0){ // Check for periodic z boundary
                        dz=Height;
                    }
                    else if(z2+k>=Height){
                        dz=-Height;
                    }
                    else{
                        dz = 0;
                    }
                }
                if(x2+i<0){ // Check for x boundary
                    dx=Length;
                }
                else if(x2+i>=Length){
                    dx=-Length;
                }
                else{
                    dx = 0;
                }
                if(y2+j<0){ // Check for y boundary
                    dy=Width;
                }
                else if(y2+j>=Width){
                    dy=-Width;
                }
                else{
                    dy = 0;
                }
                // Count the number of similar neighbors
                if(lattice[getSite(x2+i+dx,y2+j+dy,z2+k+dz)].type==site2_type){
                    switch(i*i+j*j+k*k){
                        case 1:
                            sum1_2i++;
                            break;
                        case 2:
                            sum2_2i++;
                            break;
                        case 3:
                            sum3_2i++;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
    }
    sum1_1f = total1-sum1_2i-1;
    sum2_1f = total2-sum2_2i;
    sum3_1f = total3-sum3_2i;
    sum1_1_delta = sum1_1f-sum1_1i;
    sum2_1_delta = sum2_1f-sum2_1i;
    sum3_1_delta = sum3_1f-sum3_1i;
    sum1_2_delta = sum1_2f-sum1_2i;
    sum2_2_delta = sum2_2f-sum2_2i;
    sum3_2_delta = sum3_2f-sum3_2i;
    sum_1_delta = -sum1_1_delta - sum2_1_delta*one_over_sqrt2;
    sum_2_delta = -sum1_2_delta - sum2_2_delta*one_over_sqrt2;
    // By default interactions with the third-nearest neighbors are not included, but when enabled they are added here
    if(Enable_third_neighbor_interaction){
        sum_1_delta -= sum3_1_delta*one_over_sqrt3;
        sum_2_delta -= sum3_2_delta*one_over_sqrt3;
    }
    if(site1_type==(char)1){
        return interaction_energy1*sum_1_delta+interaction_energy2*sum_2_delta;
    }
    else{
        return interaction_energy2*sum_1_delta+interaction_energy1*sum_2_delta;
    }
}

//  This function calculates the number of site faces that are between dissimilar sites, resulting in the interfacial area in units of square lattice units
double Morphology::calculateInterfacialArea(){
    unsigned long site_count = 0;
    int dx,dy,dz;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                if(lattice[getSite(x,y,z)].type==(char)1){
                    for(int i=-1;i<=1;i++){
                        for(int j=-1;j<=1;j++){
                            for(int k=-1;k<=1;k++){
                                if(abs(i)+abs(j)+abs(k)>1){
                                    continue;
                                }
                                if(i==0 && j==0 && k==0){
                                    continue;
                                }
                                if(!Enable_z_periodic_boundary) {
                                    if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                                        continue;
                                    }
                                    dz = 0;
                                }
                                else{
                                    if(z+k<0){ // Check for periodic z boundary
                                        dz=Height;
                                    }
                                    else if(z+k>=Height){
                                        dz=-Height;
                                    }
                                    else{
                                        dz = 0;
                                    }
                                }
                                if(x+i<0){ // Check for x boundary
                                    dx=Length;
                                }
                                else if(x+i>=Length){
                                    dx=-Length;
                                }
                                else{
                                    dx = 0;
                                }
                                if(y+j<0){ // Check for y boundary
                                    dy=Width;
                                }
                                else if(y+j>=Width){
                                    dy=-Width;
                                }
                                else{
                                    dy = 0;
                                }
                                if(lattice[getSite(x+i+dx,y+j+dy,z+k+dz)].type==(char)2){
                                    site_count++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return (double)site_count;
}

//  This function calculates the interfacial distance histograms that characterize the morphology, which gives the fraction of sites at a certain distance from the interface.
//  This histogram is compiled by calculating the shortest distance from each each site to an interface.
bool Morphology::calculateInterfacialDistance(){
    int dx,dy,dz;
    int count1;
    int count2;
    int d_int;
    float d;
    float d_temp;
    // The shortest distance from each site to the interface is stored in the path_distances vector
    vector<float> path_distances;
    path_distances.assign(Length*Width*Height,0);
    // Calculate distances to the interface by expanding outward from the interface
    // A first scan over the lattice is done to identify the sites at the interface, at a distance of less than 2 lattice units from the interface
    // Subsequent scans expand outward 1 lattice unit at a time from these interfacial sites until no sites are left.
    int calc_count = 1;
    float d_current = 1.99;
    while(calc_count>0){
        calc_count = 0;
        for(int x=0;x<Length;x++){
            for(int y=0;y<Width;y++){
                for(int z=0;z<Height;z++){
                    // Only perform the calculation on sites with a yet unknown interfacial distances
                    if(path_distances[getSite(x,y,z)]<0.1){
                        d = -1;
                        // Look around at neighboring sites
                        for(int i=-1;i<=1;i++){
                            for(int j=-1;j<=1;j++){
                                for(int k=-1;k<=1;k++){
                                    if(i==0 && j==0 && k==0){
                                        continue;
                                    }
                                    if(!Enable_z_periodic_boundary){
                                        if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                                            continue;
                                        }
                                        dz = 0;
                                    }
                                    else{
                                        if(z+k<0){ // Check for periodic z boundary
                                            dz=Height;
                                        }
                                        else if(z+k>=Height){
                                            dz=-Height;
                                        }
                                        else{
                                            dz = 0;
                                        }
                                    }
                                    if(x+i<0){ // Check for x boundary
                                        dx=Length;
                                    }
                                    else if(x+i>=Length){
                                        dx=-Length;
                                    }
                                    else{
                                        dx = 0;
                                    }
                                    if(y+j<0){ // Check for y boundary
                                        dy=Width;
                                    }
                                    else if(y+j>=Width){
                                        dy=-Width;
                                    }
                                    else{
                                        dy = 0;
                                    }
                                    // Initial scan identifies interfacial sites
                                    if(d_current<2){
                                        if(lattice[getSite(x+i+dx,y+j+dy,z+k+dz)].type!=lattice[getSite(x,y,z)].type){
                                            d_temp = sqrt(i*i+j*j+k*k);
                                            if(d<0 || d_temp<d){
                                                d = d_temp;
                                            }
                                        }
                                    }
                                    // Subsequent scans identify sites of the same type
                                    // A temporary distance to the interface from the target site by way of the identified neighbor site is calculated
                                    // The temporary distance is only stored if it is shorter than the previous temporary distance, ensuring that the shortest interfacial distance is calculated
                                    else if(lattice[getSite(x+i+dx,y+j+dy,z+k+dz)].type==lattice[getSite(x,y,z)].type && path_distances[getSite(x+i+dx,y+j+dy,z+k+dz)]>0.1){
                                        d_temp = path_distances[getSite(x+i+dx,y+j+dy,z+k+dz)] + sqrt(i*i+j*j+k*k);
                                        if(d<0 || d_temp<d){
                                            d = d_temp;
                                        }
                                    }
                                }
                            }
                        }
                        // The temporary distance is only accepted if it less than the expansion limit (d_current).
                        if(d>0 && d<d_current){
                            path_distances[getSite(x,y,z)] = d;
                            calc_count++;
                        }
                    }
                }
            }
        }
        // Incrementing the expansion limit (d_current) one lattice unit at a time ensures that the shortest path to the interface is determined.
        d_current += 1;
    }
    // Construct interfacial distance histograms
    InterfacialHistogram1.clear();
    InterfacialHistogram2.clear();
    count1 = 0;
    count2 = 0;
    // The interfacial distance data in path_data is rounded to the nearest integer lattice unit
    // One bin for each integer lattice unit is used to create the histograms.
    for(int i=0;i<(int)path_distances.size();i++){
        d_int = round_int(path_distances[i]);
        if(lattice[i].type==(char)1){
            if(d_int>(int)InterfacialHistogram1.size()){
                InterfacialHistogram1.push_back(1);
            }
            else{
                InterfacialHistogram1[d_int-1] += 1;
            }
            count1++;
        }
        else{
            if(d_int>(int)InterfacialHistogram2.size()){
                InterfacialHistogram2.push_back(1);
            }
            else{
                InterfacialHistogram2[d_int-1] += 1;
            }
            count2++;
        }
    }
    for(int i=0;i<(int)InterfacialHistogram1.size();i++){
        InterfacialHistogram1[i] /= count1;
    }
    for(int i=0;i<(int)InterfacialHistogram2.size();i++){
        InterfacialHistogram2[i] /= count2;
    }
    return true;
}

//  This function calculates the number of sites that are adjacent to a site of the opposite type, which represents the interfacial volume in cubic lattice units.
double Morphology::calculateInterfacialVolume(){
    unsigned long site_count = 0;
    int dx,dy,dz;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                // For each site in the lattice, the neighboring sites are checked to see if there is one that is not the same type.
                for(int i=-1;i<=1;i++){
                    for(int j=-1;j<=1;j++){
                        for(int k=-1;k<=1;k++){
                            if(i==0 && j==0 && k==0){
                                continue;
                            }
                            if(!Enable_z_periodic_boundary) {
                                if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                                    continue;
                                }
                                dz = 0;
                            }
                            else{
                                if(z+k<0){ // Check for periodic z boundary
                                    dz=Height;
                                }
                                else if(z+k>=Height){
                                    dz=-Height;
                                }
                                else{
                                    dz = 0;
                                }
                            }
                            if(x+i<0){ // Check for x boundary
                                dx=Length;
                            }
                            else if(x+i>=Length){
                                dx=-Length;
                            }
                            else{
                                dx=0;
                            }
                            if(y+j<0){ // Check for y boundary
                                dy=Width;
                            }
                            else if(y+j>=Width){
                                dy=-Width;
                            }
                            else{
                                dy=0;
                            }
                            if(lattice[getSite(x,y,z)].type!=lattice[getSite(x+i+dx,y+j+dy,z+k+dz)].type){
                                site_count++;
                                i = 2;
                                j = 2;
                                k = 2;
                            }
                        }
                    }
                }
            }
        }
    }
    return (double)site_count;
}

//  This function calculates the fraction of type 1 sites in the lattice
void Morphology::calculateMixFraction(){
    //Calculate final Mix_fraction
    int count1 = 0;
    int count2 = 0;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                if(lattice[getSite(x,y,z)].type==(char)1){
                    count1++;
                }
                else{
                    count2++;
                }
            }
        }
    }
    Mix_fraction = (double)count1/(count1+count2);
}

//  This function calculates the tortuosity histograms that characterize the morphology.
//  For all type 1 sites, the shortest distance from the site along a path through other type 1 sites to the boundary at z=0 is calculated.
//  For all type 2 sites, the shortest distance from the site along a path through other type 2 sites to the boundary at z=Height-1 is calculated.
//  The resulting shortest path divided by the straight vertical path is the tortuosity of the pathway.
//  The shortest paths are calculated using Dijkstra's algorithm
bool Morphology::calculateTortuosity(){
    int dx,dy;
    int bin;
    int neighbor_count;
    int count1;
    int count2;
    int current_index;
    float d;
    float d_temp;
    // The shortest path for each site is stored in the path_distances vector.
    // The path distances are initialized to zero.
    vector<float> path_distances;
    path_distances.assign(Length*Width*Height,0);
    // Create and initialize a blank node.
    // Each node contains a vector with indices of all first- ,second-, and third-nearest neighbors (at most 26 neighbors).
    // Another vector stores the squared distance to each of the neighbors.
    // Each node also has an estimated distance from the destination.
    Node temp_node;
    for(int i=0;i<26;i++){
        temp_node.neighbor_indices[i] = 0;
        temp_node.neighbor_distances_sq[i] = 0;
    }
    temp_node.distance_est = -1;
    // Create a node vector that is the same size as the lattice and initialize with blank nodes.
    Node_vector.assign(Length*Width*Height,temp_node);
    // The neighbor_nodes set is sorted by the estimated distance of nodes in the set.
    // This set is used in Dijsktra's algorithm to keep a sorted list of all nodes that are neighboring nodes that have their path distances already determined.
    // Once the path distance for a particular node is fixed, all of its neighboring nodes that have not yet been fixed are added to the neighbor_nodes set.
    set<vector<Node>::iterator,NodeIteratorCompare> neighbor_nodes;
    vector<Node>::iterator current_it;
    set<vector<Node>::iterator>::iterator current_set_it;
    set<vector<Node>::iterator>::iterator set_it;
    // Determine node connectivity.
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                neighbor_count = 0;
                for(int i=-1;i<=1;i++){
                    for(int j=-1;j<=1;j++){
                        for(int k=-1;k<=1;k++){
                            if(i==0 && k==0 && j==0){
                                continue;
                            }
                            if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                                continue;
                            }
                            if(x+i<0){ // Check for x boundary
                                dx=Length;
                            }
                            else if(x+i>=Length){
                                dx=-Length;
                            }
                            else{
                                dx = 0;
                            }
                            if(y+j<0){ // Check for y boundary
                                dy=Width;
                            }
                            else if(y+j>=Width){
                                dy=-Width;
                            }
                            else{
                                dy = 0;
                            }
                            if(lattice[getSite(x,y,z)].type==lattice[getSite(x+i+dx,y+j+dy,z+k)].type){
                                Node_vector[getSite(x,y,z)].neighbor_indices[neighbor_count] = getSite(x+i+dx,y+j+dy,z+k);
                                Node_vector[getSite(x,y,z)].neighbor_distances_sq[neighbor_count] = (char)(i*i+j*j+k*k);
                                neighbor_count++;
                            }
                        }
                    }
                }
            }
        }
    }
    // Initialize the path distances of top and bottom surfaces of the lattice.
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            if(lattice[getSite(x,y,0)].type==(char)1){
                path_distances[getSite(x,y,0)] = 1;
            }
            if(lattice[getSite(x,y,Height-1)].type==(char)2){
                path_distances[getSite(x,y,Height-1)] = 1;
            }
        }
    }
    int z;
    // The pathfinding algorithm is performed for one type at a time.
    for(int type=1;type<=2;type++){
        // Use Dijkstra's algorithm to fill in the remaining path distance data.
        cout << ID << ": Executing Dijkstra's algorithm to calculate shortest paths through domain type " << type << ".\n";
        // Initialize the neighbor node set.
        for(int x=0;x<Length;x++){
            for(int y=0;y<Width;y++){
                if(type==1){
                    z = 1;
                }
                else{
                    z = Height-2;
                }
                if(lattice[getSite(x,y,z)].type==(char)type){
                    d = -1;
                    int i=0;
                    while(Node_vector[getSite(x,y,z)].neighbor_indices[i]>0 && i<26){
                        if(path_distances[Node_vector[getSite(x,y,z)].neighbor_indices[i]]>0){
                            d_temp = path_distances[Node_vector[getSite(x,y,z)].neighbor_indices[i]]+sqrt((int)(Node_vector[getSite(x,y,z)]).neighbor_distances_sq[i]);
                            if(d<0 || d_temp<d){
                                d = d_temp;
                            }
                        }
                        i++;
                    }
                    if(d>0){
                        Node_vector[getSite(x,y,z)].distance_est = d;
                        neighbor_nodes.insert(Node_vector.begin()+getSite(x,y,z));
                    }
                }
            }
        }
        int i;
        int neighbor_index;
        while(!neighbor_nodes.empty()){
            // The neighbor nodes set is sorted, so the first node has the shortest estimated path distance and is set to the current node.
            current_set_it = neighbor_nodes.begin();
            current_it = *current_set_it;
            current_index = current_it-Node_vector.begin();
            // Insert neighbors of the current node into the neighbor node set.
            i=0;
            while(Node_vector[current_index].neighbor_indices[i]>0 && i<26){
                neighbor_index = Node_vector[current_index].neighbor_indices[i];
                // Check that the target neighbor node has not already been finalized.
                if(path_distances[neighbor_index]>0){
                    i++;
                    continue;
                }
                // Calculate the estimated path distance.
                d = Node_vector[current_index].distance_est + sqrt((int)(Node_vector[current_index]).neighbor_distances_sq[i]);
                // Check if node is already in the neighbor node set.
                set_it = neighbor_nodes.find(Node_vector.begin()+neighbor_index);
                // If the node is not already in the list, update the distance estimate for the target neighbor node and insert the node into the set.
                if(set_it==neighbor_nodes.end()){
                    Node_vector[neighbor_index].distance_est = d;
                    neighbor_nodes.insert(Node_vector.begin()+neighbor_index);
                }
                // If it already is in the list, replace it only if the new path distance estimate is smaller.
                else if(d<(*set_it)->distance_est){
                    neighbor_nodes.erase(set_it);
                    Node_vector[neighbor_index].distance_est = d;
                    neighbor_nodes.insert(Node_vector.begin()+neighbor_index);
                }
                i++;
            }
            // Finalize the path distance of current node and remove the current node from the neighbor node set.
            path_distances[current_index] = Node_vector[current_index].distance_est;
            neighbor_nodes.erase(current_set_it);
        }
    }
    // Clear allocated memory for the neighbor nodes set
    set<vector<Node>::iterator,NodeIteratorCompare>().swap(neighbor_nodes);
    // Construct tortuosity histograms from the path data
    // Tortuosity values are rounded to the nearest 0.02, resulting in bins centered at 1, 1.02, 10.4, etc.
    TortuosityHistogram1.assign(1,0);
    TortuosityHistogram2.assign(1,0);
    count1 = 0;
    count2 = 0;
    for(int i=0;i<(int)path_distances.size();i++){
        if(lattice[i].type==(char)1){
            if(path_distances[i]>0){
                bin = round_int(50*(path_distances[i]/(getCoords(i).z+1))-49);
                if(bin>=(int)TortuosityHistogram1.size()){
                    TortuosityHistogram1.push_back(1);
                }
                else{
                    TortuosityHistogram1[bin] += 1;
                }
                count1++;
            }
            else{
                bin = 0;
                TortuosityHistogram1[bin] += 1;
                count1++;
            }
        }
        else if(lattice[i].type==(char)2){
            if(path_distances[i]>0){
                bin = round_int(50*(path_distances[i]/(Height-getCoords(i).z))-49);
                if(bin>=(int)TortuosityHistogram2.size()){
                    TortuosityHistogram2.push_back(1);
                }
                else{
                    TortuosityHistogram2[bin] += 1;
                }
                count2++;
            }
            else{
                bin = 0;
                TortuosityHistogram2[bin] += 1;
                count2++;
            }
        }
    }
    for(int i=0;i<(int)TortuosityHistogram1.size();i++){
        TortuosityHistogram1[i] /= count1;
    }
    for(int i=0;i<(int)TortuosityHistogram2.size();i++){
        TortuosityHistogram2[i] /= count2;
    }
    // In addition to the overall tortuosity histograms from all sites, the end-to-end tortuosity distribution is collected.
    // the end-to-end describes the distribution of tortuosities for the all pathways from the top surface to the bottom surface of the lattice.
    TortuosityData1.clear();
    TortuosityData2.clear();
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            if(lattice[getSite(x,y,Height-1)].type==(char)1 && path_distances[getSite(x,y,Height-1)]>0){
                TortuosityData1.push_back(path_distances[getSite(x,y,Height-1)]/Height);
            }
            if(lattice[getSite(x,y,0)].type==(char)2 && path_distances[getSite(x,y,0)]>0){
                TortuosityData2.push_back(path_distances[getSite(x,y,0)]/Height);
            }
        }
    }
    // Any sites which are not connected their respective surface, will have a zero path distance and are identified as part of island domains.
    // Calculate island volume fraction
    Island_volume1 = 0;
    Island_volume2 = 0;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                if(lattice[getSite(x,y,z)].type==(char)1 && path_distances[getSite(x,y,z)]<1){
                    Island_volume1++;
                }
                if(lattice[getSite(x,y,z)].type==(char)2 && path_distances[getSite(x,y,z)]<1){
                    Island_volume2++;
                }
            }
        }
    }
    return true;
}

//  This function calculates the fraction of nearby sites the site at (x,y,z) that are not the same type.
//  The radius that determines which sites are included as nearby sites is determined by the rescale factor parameter.
//  This function is designed to be used by the executeSmoothing function and implement rescale factor dependent smoothing.
double Morphology::calculateDissimilarFraction(int x,int y,int z,int rescale_factor){
    int site_count = 0;
    int count_dissimilar = 0;
    int dx,dy,dz;
    // When the rescale factor is 1, the radius is 1, and the radius increases for larger rescale factors.
    static int radius = (int)ceil((double)(rescale_factor+1)/2);
    static int cutoff_squared = (int)floor(((double)(rescale_factor+1)/2)*((double)(rescale_factor+1)/2));
    for(int i=-radius;i<=radius;i++){
        for(int j=-radius;j<=radius;j++){
            for(int k=-radius;k<=radius;k++){
                if(i*i+j*j+k*k>cutoff_squared){
                    continue;
                }
                if(!Enable_z_periodic_boundary){
                    if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                        continue;
                    }
                    dz = 0;
                }
                else{
                    if(z+k<0){ // Check for periodic z boundary
                        dz=Height;
                    }
                    else if(z+k>=Height){
                        dz=-Height;
                    }
                    else{
                        dz = 0;
                    }
                }
                if(x+i<0){ // Check for x boundary
                    dx=Length;
                }
                else if(x+i>=Length){
                    dx=-Length;
                }
                else{
                    dx = 0;
                }
                if(y+j<0){ // Check for y boundary
                    dy=Width;
                }
                else if(y+j>=Width){
                    dy=-Width;
                }
                else{
                    dy = 0;
                }
                if(lattice[getSite(x,y,z)].type!=lattice[getSite(x+i+dx,y+j+dy,z+k+dz)].type){
                    count_dissimilar++;
                }
                site_count++;
            }
        }
    }
    return (double)count_dissimilar/(site_count-1);
}

//  This function enables interactions between third-neighbor sites that are a distance of sqrt(3) lattice units apart.
//  By default third-neighbor interactions are disabled, so this function must be called to enable this option.
void Morphology::enableThirdNeighborInteraction(){
    Enable_third_neighbor_interaction = true;
}

//  This function creates a randomly mixed morphology on the lattice.
//  Sites are randomly assigned as type 1 or type 2 based on the mix_fraction parameter.
void Morphology::createRandomMorphology(double mix_fraction){
    double value;
    if(mix_fraction>=1){
        cout << ID << ": Error creating morphology: Mix fraction must be less than one." << endl;
        return;
    }
    Mix_fraction = mix_fraction;
    boost::uniform_01<boost::mt19937> rand01(gen);
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                value = rand01();
                if(value<mix_fraction){
                    lattice[getSite(x,y,z)].type = (char)1;
                }
                else{
                    lattice[getSite(x,y,z)].type = (char)2;
                }
            }
        }
    }
    // This function calculates the actual mix fraction and updates Mix_fraction
    calculateMixFraction();
}

//  This function implements num_MCsteps iterations of the Ising site swapping process.
//  This function uses the bond formation algorithm to determine the energy change in the system that results from swapping two neighboring sites.
//  The energy change is determined by the input parameters interaction_energy1 and interaction_energy2, which are in units of kT.
//  These parameters describe the preference for like-like interactions over like-unlike interactions for each site type.
//  Positive values of the interaction energies result in a driving force for phase separation.
void Morphology::executeIsingSwapping(int num_MCsteps,double interaction_energy1,double interaction_energy2){
    boost::uniform_int<> distx(0,Length-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randx(gen,distx);
    boost::uniform_int<> disty(0,Width-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randy(gen,disty);
    boost::uniform_int<> distz(0,Height-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randz(gen,distz);
    boost::uniform_int<> dist_neighbor(-1,1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_neighbor(gen,dist_neighbor);
    boost::uniform_01<boost::mt19937> rand01(gen);
    int loop_count = 0;
    // N counts the number of MC steps that have been executed
    int N = 0;
    int neighbor_count;
    int main_site_index;
    int neighbor_site;
    int x,y,z,dx,dy,dz;
    char temp;
    double energy_delta,probability;
    Coords neighbor_coords;
    vector<int> neighbors;
    neighbors.assign(6,0);
    int m=1;
    while(N<num_MCsteps){
        // Randomly choose a site in the lattice
        x = randx();
        y = randy();
        z = randz();
        // If site is not an interfacial site, choose again
        if(!isNearInterface(x,y,z,1.1)){
            continue;
        }
        main_site_index = getSite(x,y,z);
        // Randomly choose a nearest neighbor site of a different type
        // First find all nearest neighbor sites of differing type
        neighbor_count = 0;
        for(int i=-1;i<=1;i++){
            for(int j=-1;j<=1;j++){
                for(int k=-1;k<=1;k++){
                    if(i*i+j*j+k*k>1){ // Limit to nearest neighbors only
                        continue;
                    }
                    if(i==0 && j==0 && k==0){ // Cannot be same site
                        continue;
                    }
                    if(!Enable_z_periodic_boundary) {
                        if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                            continue;
                        }
                        dz = 0;
                    }
                    else{
                        if(z+k<0){ // Check for periodic z boundary
                            dz=Height;
                        }
                        else if(z+k>=Height){
                            dz=-Height;
                        }
                        else{
                            dz = 0;
                        }
                    }
                    if(x+i<0){ // Check for periodic x boundary
                        dx=Length;
                    }
                    else if(x+i>=Length){
                        dx=-Length;
                    }
                    else{
                        dx = 0;
                    }
                    if(y+j<0){ // Check for periodic y boundary
                        dy=Width;
                    }
                    else if(y+j>=Width){
                        dy=-Width;
                    }
                    else{
                        dy = 0;
                    }
                    if(lattice[main_site_index].type!=lattice[getSite(x+i+dx,y+j+dy,z+k+dz)].type){
                        // Store site index of differing neighbor site
                        neighbors[neighbor_count] = getSite(x+i+dx,y+j+dy,z+k+dz);
                        neighbor_count++;
                    }
                }
            }
        }
        // Randomly select one of the differing neighbor sites
        neighbor_site = neighbors[round_int(rand01()*(neighbor_count-1))];
        neighbor_coords = getCoords(neighbor_site);
        // Calculate energy change and swapping probability
        energy_delta = calculateEnergyChange(x,y,z,neighbor_coords.x,neighbor_coords.y,neighbor_coords.z,interaction_energy1,interaction_energy2);
        probability = exp(-energy_delta)/(1+exp(-energy_delta));
        if(rand01()<=probability){
            // Swap Sites
            temp = lattice[main_site_index].type;
            lattice[main_site_index].type = lattice[neighbor_site].type;
            lattice[neighbor_site].type = temp;
        }
        loop_count++;
        // One MC step has been completed when loop_count is equal to the number of sites in the lattice
        if(loop_count==Length*Width*Height){
            N++;
            loop_count = 0;
        }
        if(N==100*m){
            cout << ID << ": " << N << " MC steps completed." << endl;
            m++;
        }
    }
}

//  This function implements interfacial mixing with a specified interfacial width and a specified mixing concentration in the interfacial region.
//  Mixing is implemented by first determining the bounds on either side of the interface where mixing should occur
//  Then random swapping of type 1 and type 2 sites within the bounds creates mixing in the interfacial region.
void Morphology::executeMixing(double width,double interfacial_conc){
    vector<int> sites_maj;
    vector<int> sites_min;
    int site_maj;
    int site_min;
    int target;
    int site_count = 0;
    char majority_type;
    char minority_type;
    double minority_conc;
    // Based on the interfacial concentration, the majority and minority types are defined
    if(interfacial_conc<=0.5){
        majority_type = (char)2;
        minority_type = (char)1;
        minority_conc = interfacial_conc;
    }
    else{
        majority_type = (char)1;
        minority_type = (char)2;
        minority_conc = 1-interfacial_conc;
    }
    // The minority type sites are mixed into the majority type sites at the interface.
    // First the majority site reservoir is determined by finding all majority type sites within (1-minority_conc)*width distance from the interface and adding them to a list
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                if(lattice[getSite(x,y,z)].type==majority_type){
                    if(isNearInterface(x,y,z,(1-minority_conc)*width)){
                        sites_maj.push_back(getSite(x,y,z));
                        site_count++;
                    }
                }
            }
        }
    }
    // Then the minority site reservoir that will be used to be mix into the majority reservoir to yield the desired final interfacial concentration is determined
    // More minority sites than the number of swaps must be included in the minority reservoir so that the desired interfacial concentration is also reached on the minority side of the interface.
    target = (int)(site_count*minority_conc/(1-minority_conc));
    int N = 0;
    double range = 1;
    // The minority type sites adjacent to the interface are identified, counted, and added to a list.
    // If this does not yield enough minority sites for swapping, then the range is slowly incremented and minority sites farther from the interface are included until there are enough sites.
    while(N<target){
        sites_min.clear();
        N = 0;
        for(int x=0;x<Length;x++){
            for(int y=0;y<Width;y++){
                for(int z=0;z<Height;z++){
                    if(lattice[getSite(x,y,z)].type==minority_type){
                        if(isNearInterface(x,y,z,range)){
                            sites_min.push_back(getSite(x,y,z));
                            N++;
                        }
                    }
                }
            }
        }
        range += 0.1;
    }
    // Sites are randomly chosen from the minority and majority reservoirs and swapped until the desired interfacial concentration is reached.
    for(int i=0;i<site_count*minority_conc;i++){
        boost::uniform_int<> dist_maj(0,sites_maj.size()-1);
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_maj(gen, dist_maj);
        boost::uniform_int<> dist_min(0,sites_min.size()-1);
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_min(gen, dist_min);
        site_maj = rand_maj();
        site_min = rand_min();
        lattice[sites_maj[site_maj]].type = minority_type;
        lattice[sites_min[site_min]].type = majority_type;
        // Both site types are removed from the lists once they are swapped to prevent unswapping.
        sites_maj.erase(sites_maj.begin()+site_maj);
        sites_min.erase(sites_min.begin()+site_min);
    }
}

//  This function smoothens out rough domain interfaces and removes small islands and island sites.
//  This is done by determining a roughness factor for each site that is given by the fraction of surrounding sites that are a different type.
//  Sites with a roughness factor is greater than the specified smoothing_threshold are switched to the opposite type.
//  A rescale dependent smoothing process is executed when the rescale factor is greater than 1.
void Morphology::executeSmoothing(double smoothing_threshold, int rescale_factor){
    double roughness_factor;
    int dx,dy,dz;
    int site_count;
    static int radius=(int)ceil((double)(rescale_factor+1)/2);
    static int cutoff_squared = (int)floor(((double)(rescale_factor+1)/2)*((double)(rescale_factor+1)/2));
    // The boolean vector consider_smoothing keeps track of whether each site is near the interface and should be considered for smoothing.
    // Sites in the interior of the domains or at very smooth interfaces do not need to be continually reconsidered for smoothing.
    vector<bool> consider_smoothing;
    consider_smoothing.assign(Length*Width*Height,true);
    site_count=1;
    while(site_count>0){
        site_count = 0;
        for(int x=0;x<Length;x++){
            for(int y=0;y<Width;y++){
                for(int z=0;z<Height;z++){
                    if(!consider_smoothing[getSite(x,y,z)]){
                        continue;
                    }
                    // Calculate the roughness factor of the target site.
                    roughness_factor = calculateDissimilarFraction(x,y,z,rescale_factor);
                    // Swap the site's type if the roughness_factor is greater than the smoothing_threshold.
                    if(roughness_factor>smoothing_threshold){
                        if(lattice[getSite(x,y,z)].type==(char)1){
                            lattice[getSite(x,y,z)].type=(char)2;
                        }
                        else if(lattice[getSite(x,y,z)].type==(char)2){
                            lattice[getSite(x,y,z)].type=(char)1;
                        }
                        site_count++;
                        // When a site swaps types, all surrounding sites must be reconsidered for smoothing.
                        for(int i=-radius;i<=radius;i++){
                            for(int j=-radius;j<=radius;j++){
                                for(int k=-radius;k<=radius;k++){
                                    if(i*i+j*j+k*k>cutoff_squared){
                                        continue;
                                    }
                                    if(!Enable_z_periodic_boundary){
                                        if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                                            continue;
                                        }
                                        dz = 0;
                                    }
                                    else{
                                        if(z+k<0){ // Check for periodic z boundary
                                            dz=Height;
                                        }
                                        else if(z+k>=Height){
                                            dz=-Height;
                                        }
                                        else{
                                            dz = 0;
                                        }
                                    }
                                    if(x+i<0){ // Check for x boundary
                                        dx=Length;
                                    }
                                    else if(x+i>=Length){
                                        dx=-Length;
                                    }
                                    else{
                                        dx=0;
                                    }
                                    if(y+j<0){ // Check for y boundary
                                        dy=Width;
                                    }
                                    else if(y+j>=Width){
                                        dy=-Width;
                                    }
                                    else{
                                        dy=0;
                                    }
                                    consider_smoothing[getSite(x+i+dx,y+j+dy,z+k+dz)]=true;
                                }
                            }
                        }
                    }
                    // Sites with a low roughness_factor are not swapped and removed from reconsideration.
                    else{
                        consider_smoothing[getSite(x,y,z)]=false;
                    }
                }
            }
        }
    }
    // The smoothing process can change the mix fraction, so the final mix fraction is recalculated and the Mix_fraction property is updated.
    calculateMixFraction();
}

//  This function returns a Coords object containing the x,y,z coordinates of the site with the given site_index.
Coords Morphology::getCoords(unsigned long site_index){
    Coords coords;
    coords.x = site_index/(Width*Height);
    unsigned long remainder = site_index % (Width*Height);
    coords.y = remainder/Height;
    coords.z = remainder % Height;
    return coords;
}

//  This function returns a vector containing the pair-pair correlation function data for the specified site type.
vector<double> Morphology::getCorrelationData(char site_type){
    if(site_type==(char)1){
        if(Correlation1[0]<0){
            cout << ID << ": Error Retrieving Correlation Data: Correlation data has not been calculated." << endl;
        }
        return Correlation1;
    }
    else if(site_type==(char)2){
        if(Correlation2[0]<0){
            cout << ID << ": Error Retrieving Correlation Data: Correlation data has not been calculated." << endl;
        }
        return Correlation2;
    }
    else{
        vector<double> temp(1,0);
        cout << ID << ": Error Retrieving Correlation Data: Invalid site type." << endl;
        return temp;
    }
}

//  This function returns the domain size determined for the specified site type.
//  This function will return zero if the calculateCorrelationDistance function has not been called.
double Morphology::getDomainSize(char site_type){
    if(site_type==(char)1){
        return Domain_size1;
    }
    else if(site_type==(char)2){
        return Domain_size2;
    }
    else{
        cout << ID << ": Error getting domain size: Invalid site type." << endl;
        return 0;
    }
}

//  This function return the height or z-direction size of the lattice.
int Morphology::getHeight(){
    return Height;
}

//  This function returns a vector containing the interfacial distance histogram data for the specified site type.
vector<double> Morphology::getInterfacialHistogram(char site_type){
    if(site_type==(char)1){
        return InterfacialHistogram1;
    }
    else if(site_type==(char)2){
        return InterfacialHistogram2;
    }
    else{
        cout << "Error getting interfacial histogram data! Invalid site type." << endl;
        return vector<double>();
    }
}

//  This function returns the island volume for the specified site type.
double Morphology::getIslandVolume(char site_type){
    if(site_type==(char)1){
        return Island_volume1;
    }
    else if(site_type==(char)2){
        return Island_volume2;
    }
    else{
        cout << ID << ": Error getting island volume! Invalid site type." << endl;
        return 0;
    }
}

//  This function returns the length or x-direction size of the lattice.
int Morphology::getLength(){
    return Length;
}

//  This function return the mix fraction of the morphology.
double Morphology::getMixFraction(){
    return Mix_fraction;
}

//  This function returns the site index of the site located at (x,y,z) on a rescaled lattice.
unsigned long Morphology::getRescaleSite(int x, int y, int z,int rescale_factor){
    return x*rescale_factor*Width*rescale_factor*Height+y*rescale_factor*Height+z;
}

//  This function returns the site index of the site located at (x,y,z).
unsigned long Morphology::getSite(int x, int y, int z){
    return x*Width*Height+y*Height+z;
}

//  This function returns a vector containing the end-to-end tortuosity data for the specified site type.
vector<float> Morphology::getTortuosityData(char site_type){
    if(site_type==(char)1){
        return TortuosityData1;
    }
    else if(site_type==(char)2){
        return TortuosityData2;
    }
    else{
        cout << "Error getting tortuosity data! Invalid site type." << endl;
        return vector<float>();
    }
}

//  This function returns a vector containing the overall tortuosity histogram data for all sites with the specified site type.
vector<double> Morphology::getTortuosityHistogram(char site_type){
    if(site_type==(char)1){
        return TortuosityHistogram1;
    }
    else if(site_type==(char)2){
        return TortuosityHistogram2;
    }
    else{
        cout << "Error getting tortuosity histogram! Invalid site type." << endl;
        return vector<double>();
    }
}

//  This function returns the width or y-direction size of the lattice.
int Morphology::getWidth(){
    return Width;
}

//  This function imports the morphology text file given by the input file stream.
//  It must be specified whether or not the input file is in the compressed format.
bool Morphology::importMorphologyFile(ifstream * input,bool compressed_files){
    string var;
    int x = 0;
    int y = 0;
    int z = 0;
    char type = (char)0;
    int j;
    string line;
    Site site;
    getline(*input,line);
    Length = atoi(line.c_str());
    getline(*input,line);
    Width = atoi(line.c_str());
    getline(*input,line);
    Height = atoi(line.c_str());
    site.type = (char)0;
    lattice.assign(Length*Width*Height,site);
    getline(*input,line);
    Domain_size1 = atof(line.c_str());
    getline(*input,line);
    Domain_size2 = atof(line.c_str());
    getline(*input,line);
    Mix_fraction = atof(line.c_str());
    if(!compressed_files){
        while((*input).good()){
            getline(*input,line);
            stringstream linestream(line);
            j=1;
            while(linestream.good()) {
                getline(linestream,var,',');
                switch(j) {
                case 1:
                    x = atoi(var.c_str());
                    break;
                case 2:
                    y = atoi(var.c_str());
                    break;
                case 3:
                    z = atoi(var.c_str());
                    break;
                case 4:
                    type = (char)atoi(var.c_str());
                    break;
                default:
                    break;
                }
                j++;
            }
            lattice[getSite(x,y,z)].type = type;
        }
    }
    else{
        int site_count = 0;
        for(int x=0;x<Length;x++){
            for(int y=0;y<Width;y++){
                for(int z=0;z<Height;z++){
                    if(site_count==0){
                        getline(*input,line);
                        stringstream linestream(line);
                        if(!linestream.good()){
                            cout << "Error parsing file.  End of file reached before expected." << endl;
                            return 0;
                        }
                        type = (char)atoi((line.substr(0,1)).c_str());
                        site_count = atoi((line.substr(1,line.length()-1)).c_str());
                    }
                    lattice[getSite(x, y, z)].type = type;
                    site_count--;
                }
            }
        }
    }
    calculateMixFraction();
    return true;
}

//  This function is an efficient implementation of the pow function when the exponent is an integer.
double Morphology::intpow(double base,int exponent){
    for(int i=1;i<exponent;i++){
        base *= base;
    }
    return base;
}

//  This function determines whether the site at (x,y,z) is within the specified distance from the interface.
//  If so, the function returns true and if not, the function returns false.
bool Morphology::isNearInterface(int x,int y,int z,double distance){
    int d = (int)floor(distance);
    int dx,dy,dz;
    for(int i=-d;i<=d;i++){
        for(int j=-d;j<=d;j++){
            for(int k=-d;k<=d;k++){
                if(d==1){
                    if(abs(i)+abs(j)+abs(k)>1 || i+j+k==0){
                        continue;
                    }
                }
                else if(sqrt((double)(i*i+j*j+k*k))>distance){
                    continue;
                }
                if(!Enable_z_periodic_boundary) {
                    if(z+k<0 || z+k>=Height){ // Check for hard z boundary
                        continue;
                    }
                    dz = 0;
                }
                else{
                    if(z+k<0){ // Check for periodic z boundary
                        dz=Height;
                    }
                    else if(z+k>=Height){
                        dz=-Height;
                    }
                    else{
                        dz = 0;
                    }
                }
                if(x+i<0){ // Check for x boundary
                    dx=Length;
                }
                else if(x+i>=Length){
                    dx=-Length;
                }
                else{
                    dx = 0;
                }
                if(y+j<0){ // Check for y boundary
                    dy=Width;
                }
                else if(y+j>=Width){
                    dy=-Width;
                }
                else{
                    dy = 0;
                }
                if(lattice[getSite(x,y,z)].type!=lattice[getSite(x+i+dx,y+j+dy,z+k+dz)].type){
                    return true;
                }
            }
        }
    }
    return false;
}

//  This function outputs to a text file a cross-section of the morphology at x=0 plane.
bool Morphology::outputMorphologyCrossSection(ofstream * output){
    for(int x=0;x<1;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                *output << x << "," << y << "," << z << "," << lattice[getSite(x,y,z)].type << endl;
            }
        }
    }
    return true;
}

//  This function outputs the morphology data to a text file specified by the output file stream.
//  The user can specify whether to use the compress text format or not.
bool Morphology::outputMorphologyFile(ofstream * output,bool enable_export_compressed_files){
    *output << Length << endl;
    *output << Width << endl;
    *output << Height << endl;
    *output << Domain_size1 << endl;
    *output << Domain_size2 << endl;
    *output << Mix_fraction << endl;
    if(!enable_export_compressed_files){
        for(int x=0;x<Length;x++){
            for(int y=0;y<Width;y++){
                for(int z=0;z<Height;z++){
                    *output << x << "," << y << "," << z << "," << lattice[getSite(x,y,z)].type << endl;
                }
            }
        }
    }
    else{
        int one_count = 0;
        int two_count = 0;
        for(int x=0;x<Length;x++){
            for(int y=0;y<Width;y++){
                for(int z=0;z<Height;z++){
                    if(lattice[getSite(x,y,z)].type==(char)1){
                        if(two_count<1){
                            one_count++;
                        }
                        else{
                            *output << "2" << two_count << "\n";
                            two_count = 0;
                            one_count++;
                        }
                    }
                    else if(lattice[getSite(x,y,z)].type==(char)2){
                        if(one_count<1){
                            two_count++;
                        }
                        else{
                            *output << "1" << one_count << "\n";
                            one_count = 0;
                            two_count++;
                        }
                    }
                }
            }
        }
        if(one_count>0){
            *output << "1" << one_count << "\n";
        }
        else if(two_count>0){
            *output << "2" << two_count << "\n";
        }
    }
    return true;
}

//  This function rescales the existing lattice by a integer value called scale.
//  This original lattice is overwritten by the newly created larger rerescale_factord lattice
void Morphology::rescaleLattice(int rescale_factor){
    // Crate and initialize a site with an undefined type
    Site site;
    site.type = (char)0;
    // Construct a larger lattice consisting of undefined sites
    vector<Site> lattice_rescale(Length*rescale_factor*Width*rescale_factor*Height*rescale_factor,site);
    // Assign site types to the new lattice based on the existing lattice
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                for(int i=rescale_factor*x;i<rescale_factor*x+rescale_factor;i++){
                    for(int j=rescale_factor*y;j<rescale_factor*y+rescale_factor;j++){
                        for(int k=rescale_factor*z;k<rescale_factor*z+rescale_factor;k++){
                            lattice_rescale[getRescaleSite(i,j,k,rescale_factor)].type = lattice[getSite(x,y,z)].type;
                        }
                    }
                }
            }
        }
    }
    // Update the lattice and its property variables
    Length = Length*rescale_factor;
    Width = Width*rescale_factor;
    Height = Height*rescale_factor;
    lattice = lattice_rescale;
}

// This function rounds a double to an integer value.
int Morphology::round_int(double num) {
    return (num > 0.0) ? (num + 0.5) : (num - 0.5);
}
