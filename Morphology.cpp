// Copyright (c) 2014 Michael C. Heiber, The University of Akron, and Julius-Maximilains Universtiy of Würzburg
// For more information, see the LICENSE file that accompanies this software.

#include "Morphology.h"

Morphology::Morphology(int length, int width, int height,int procid){
    ProcID = procid;
    Length = length;
    Width = width;
    Height = height;
    Site site;
    site.type = 0;
    site.energy = 0;
    site.path_distance = 0;
    lattice.assign(length*width*height,site);
    lattice_temp.assign(125,0);
    Domain_size1_updated = false;
    Domain_size2_updated = false;
    Energies_initialized = false;
    Correlation1.assign(1,-1);
    Correlation2.assign(1,-1);
    Path_avg = 0;
    Path_stdev = 0;
}

Morphology::~Morphology(){
    //dtor
}

bool Morphology::calculateCorrelationDistance(int cutoff_distance){
    int dx,dy;
    vector<int> count1,count2,total1,total2;
    double distance;
    int bin;
    int N1 = 0;
    int N2 = 0;
    double d1,y1,d2,y2,slope,intercept;
    bool success;
    int count = 0;
    if(cutoff_distance>Length || cutoff_distance>Width){
        cout << ProcID << ": Error, cutoff distance is greater than the lattice length and/or width." << endl;
        return false;
    }
    if(Domain_size1_updated && Domain_size2_updated){
        Domain_size1_updated = false;
        Domain_size2_updated = false;
    }
    // Resolution of correlation distance data is 0.5 lattice units
    Correlation1.assign(2*cutoff_distance+1,0);
    Correlation2.assign(2*cutoff_distance+1,0);
    // Loop through all sites on the lattice and determine the correlation function for each
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                count++;
                if(100*count%(25*Length*Width*Height)==0){
                    cout << ProcID << ": Domain size calculation " << (int)100*(double)count/(Length*Width*Height) << "% complete." << endl;
                }
                if(Domain_size1_updated && lattice[getSite(x,y,z)].type==1){
                    continue;
                }
                if(Domain_size2_updated && lattice[getSite(x,y,z)].type==2){
                    continue;
                }
                if(lattice[getSite(x,y,z)].type==1){
                    count1.assign(2*cutoff_distance+1,0);
                    total1.assign(2*cutoff_distance+1,0);
                    N1++;
                }
                else{
                    count2.assign(2*cutoff_distance+1,0);
                    total2.assign(2*cutoff_distance+1,0);
                    N2++;
                }
                for(int i=-cutoff_distance;i<=cutoff_distance;i++){
                    for(int j=-cutoff_distance;j<=cutoff_distance;j++){
                        for(int k=-cutoff_distance;k<=cutoff_distance;k++){
                            // The distance between two sites is rounded to the nearest half a lattice unit
                            bin = (int)floor(2*sqrt(pow((double)i,2)+pow((double)j,2)+pow((double)k,2))+0.5);
                            distance = (double)bin/2;
                            if(distance>cutoff_distance || (int)distance==0){
                                continue;
                            }
                            if(z+k>=Height || z+k<0){ // Check for z boundary
                                continue;
                            }
                            dx = 0;
                            dy = 0;
                            if(x+i<0){ // Check for x boundary
                                dx=Length;
                            }
                            else if(x+i>=Length){
                                dx=-Length;
                            }
                            if(y+j<0){ // Check for y boundary
                                dy=Width;
                            }
                            else if(y+j>=Width){
                                dy=-Width;
                            }
                            if(lattice[getSite(x,y,z)].type==lattice[getSite(x+i+dx,y+j+dy,z+k)].type){
                                if(lattice[getSite(x,y,z)].type==1){
                                    count1[bin]++;
                                }
                                else{
                                    count2[bin]++;
                                }
                            }
                            if(lattice[getSite(x,y,z)].type==1){
                                total1[bin]++;
                            }
                            else{
                                total2[bin]++;
                            }
                        }
                    }
                }
                for(int n=0;n<=2*cutoff_distance;n++){
                    if(lattice[getSite(x,y,z)].type==1){
                        if(total1[n]>0){
                            Correlation1[n] += (double)count1[n]/(double)total1[n];
                        }
                        else{
                            Correlation1[n] += 1;
                        }
                    }
                    else{
                        if(total2[n]>0){
                            Correlation2[n] += (double)count2[n]/(double)total2[n];
                        }
                        else{
                            Correlation2[n] += 1;
                        }
                    }
                }
            }
        }
    }
    // Total Correlation function is an average of contributions from each site of the corresponding type
    for(int n=0;n<=2*cutoff_distance;n++){
        Correlation1[n] = Correlation1[n]/N1;
        Correlation2[n] = Correlation2[n]/N2;
    }
    if(!Domain_size1_updated){
        // Find the bounds of where the pair-pair correlation1 function first crosses over the Mix_fraction
        d1 = 0;
        y1 = 0;
        d2 = 0;
        y2 = 0;
        success = false;
        for(int n=2;n<=2*cutoff_distance;n++){
            if(Correlation1[n]<Mix_fraction){
                d1 = ((double)n-1)/2;
                y1 = Correlation1[n-1];
                d2 = (double)n/2;
                y2 = Correlation1[n];
                success = true;
                break;
            }
            if(n==2*cutoff_distance){
                cout << ProcID << ": Error calculating correlation distances: Cutoff distance of " << cutoff_distance << " is too small to calculate size of domain type 1," << endl;
                Domain_size1 = -1;
                return false;
            }
        }
        // Use linear interpolation to determine the crossover point
        if(success){
            slope = (y2-y1)/(d2-d1);
            intercept = y1-slope*d1;
            Domain_size1 = (Mix_fraction-intercept)/slope;
            Domain_size1_updated = true;
        }
    }
    if(!Domain_size2_updated){
        // Find the bounds of where the pair-pair correlation2 function first crosses over 1-Mix_fraction
        d1 = 0;
        y1 = 0;
        d2 = 0;
        y2 = 0;
        success = false;
        for(int n=2;n<=2*cutoff_distance;n++){
            if(Correlation2[n]<(1-Mix_fraction)){
                d1 = ((double)n-1)/2;
                y1 = Correlation2[n-1];
                d2 = (double)n/2;
                y2 = Correlation2[n];
                success = true;
                break;
            }
            if(n==2*cutoff_distance){
                cout << ProcID << ": Error calculating correlation distances: Cutoff distance of " << cutoff_distance << " is too small to calculate size of domain type 2," << endl;
                Domain_size2 = -1;
                return false;
            }
        }
        // Use linear interpolation to determine the crossover point
        if(success){
            slope = (y2-y1)/(d2-d1);
            intercept = y1-slope*d1;
            Domain_size2 = ((1-Mix_fraction)-intercept)/slope;
            Domain_size2_updated = true;
        }
    }
    return true;
}

double Morphology::calculateEnergyChange1(int x1,int y1,int z1,int x2,int y2,int z2,double interaction_energy,bool enable_third_neighbor_calc){
    // Used with site energy algorithm
    double delta_e = 0;
    int di,dj,dk,dx,dy;
    if(x1==x2){
        di = 0;
    }
    else if((x2>x1 && x2-x1==1)||(x2<x1 && x1-x2>1)){
        di = 1;
    }
    else{
        di = -1;
    }
    if(y1==y2){
        dj = 0;
    }
    else if((y2>y1 && y2-y1==1)||(y2<y1 && y1-y2>1)){
        dj = 1;
    }
    else{
        dj = -1;
    }
    if(z1==z2){
        dk = 0;
    }
    else if(z2>z1){
        dk = 1;
    }
    else{
        dk = -1;
    }
    calculateTempEnergies(x1,y1,z1,di,dj,dk,interaction_energy,enable_third_neighbor_calc);
    for(int i=-2;i<=2;i++){
        for(int j=-2;j<=2;j++){
            for(int k=-2;k<=2;k++){
                // Blank sites in the temp lattice are marked as -1 and are not counted
                if(lattice_temp[getTempSite(i+2,j+2,k+2)]<0){
                    continue;
                }
                if(z1+k>=Height || z1+k<0){ // Check for z boundary
                    continue;
                }
                dx = 0;
                dy = 0;
                if(x1+i<0){ // Check for x boundary
                    dx=Length;
                }
                else if(x1+i>=Length){
                    dx=-Length;
                }
                if(y1+j<0){ // Check for y boundary
                    dy=Width;
                }
                else if(y1+j>=Width){
                    dy=-Width;
                }
                delta_e += lattice_temp[getTempSite(i+2,j+2,k+2)] - lattice[getSite(x1+i+dx,y1+j+dy,z1+k)].energy;
            }
        }
    }
    return delta_e;
}

double Morphology::calculateEnergyChange2(int x1,int y1,int z1, int x2, int y2, int z2,double interaction_energy,bool enable_third_neighbor_calc){
    // Used with bond formation algorithm
    double energy1,energy2;
    int dx,dy;
    int site1_type,site2_type;
    int sum1_delta,sum2_delta,sum3_delta;
    int sum1 = 0;
    int sum2 = 0;
    int sum3 = 0;
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
                if(z1+k>=Height || z1+k<0 ){ // Check for z boundary
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
                // Don't count site x2,y2,z2
                if(x1+i+dx==x2 && y1+j+dy==y2 && z1+k==z2){
                    total1--;
                    continue;
                }
                // Count the number of dissimilar neighbors
                if(lattice[getSite(x1+i+dx,y1+j+dy,z1+k)].type==site1_type){
                    switch(i*i+j*j+k*k){
                        case 1:
                            sum1++;
                            break;
                        case 2:
                            sum2++;
                            break;
                        case 3:
                            sum3++;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
    }
    sum1_delta = (total1-sum1)-sum1;
    sum2_delta = (total2-sum2)-sum2;
    sum3_delta = (total3-sum3)-sum3;
    energy1 = -sum1_delta;
    energy1 -= sum2_delta/sqrt((double)2);
    if(enable_third_neighbor_calc){
        energy1 -= sum3_delta/sqrt((double)3);
    }
    // Calculate change around x2,y2,z2
    site2_type = lattice[getSite(x2,y2,z2)].type;
    sum1 = 0;
    sum2 = 0;
    sum3 = 0;
    total1 = 6;
    total2 = 12;
    total3 = 8;
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
            for(int k=-1;k<=1;k++){
                if(i==0 && j==0 && k==0){
                    continue;
                }
                if(z2+k>=Height || z2+k<0 ){ // Check for z boundary
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
                // Don't count site x1,y1,z1
                if(x2+i+dx==x1 && y2+j+dy==y1 && z2+k==z1){
                    total1--;
                    continue;
                }
                // Count the number of dissimilar neighbors
                if(lattice[getSite(x2+i+dx,y2+j+dy,z2+k)].type==site2_type){
                    switch(i*i+j*j+k*k){
                        case 1:
                            sum1++;
                            break;
                        case 2:
                            sum2++;
                            break;
                        case 3:
                            sum3++;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
    }
    sum1_delta = (total1-sum1)-sum1;
    sum2_delta = (total2-sum2)-sum2;
    sum3_delta = (total3-sum3)-sum3;
    energy2 = -sum1_delta;
    energy2 -= sum2_delta/sqrt((double)2);
    if(enable_third_neighbor_calc){
        energy2 -= sum3_delta/sqrt((double)3);
    }
    return interaction_energy*(energy1+energy2);
}

double Morphology::calculateInterfacialArea(){
    int count = 0;
    int dx,dy;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                if(lattice[getSite(x,y,z)].type==1){
                    for(int i=-1;i<=1;i++){
                        for(int j=-1;j<=1;j++){
                            for(int k=-1;k<=1;k++){
                                if(abs(i)+abs(j)+abs(k)>1){
                                    continue;
                                }
                                if(i==0 && j==0 && k==0){
                                    continue;
                                }
                                if(z+k>=Height || z+k<0){ // Check for z boundary
                                    continue;
                                }
                                dx = 0;
                                dy = 0;
                                if(x+i<0){ // Check for x boundary
                                    dx=Length;
                                }
                                else if(x+i>=Length){
                                    dx=-Length;
                                }
                                if(y+j<0){ // Check for y boundary
                                    dy=Width;
                                }
                                else if(y+j>=Width){
                                    dy=-Width;
                                }
                                if(lattice[getSite(x+i+dx,y+j+dy,z+k)].type==2){
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return (double)count;
}

void Morphology::calculateMixFraction(){
    //Calculate final Mix_fraction
    int count1 = 0;
    int count2 = 0;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                if(lattice[getSite(x,y,z)].type==1){
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

bool Morphology::calculateTortuosity(){
    // Only calculates the tortuosity of the donor phase
    int dx,dy;
    double d,d_temp;
    int count = 1;
    double d_current = 2;
    // clear previous tortuosity calculations
    TortuosityData.clear();
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                lattice[getSite(x,y,z)].path_distance = 0;
            }
        }
    }
    // initialize sites next to the bottom interface
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            if(lattice[getSite(x,y,0)].type==1){
                lattice[getSite(x,y,0)].path_distance = 1;
            }
        }
    }
    // begin outward expansion
    while(count>0){
        count = 0;
        for(int x=0;x<Length;x++){
            for(int y=0;y<Width;y++){
    	        for(int z=0;z<Height;z++){
                    if(lattice[getSite(x,y,z)].type==1 && lattice[getSite(x,y,z)].path_distance==0){
                        d = -1;
                        for(int i=-1;i<=1;i++){
                            for(int j=-1;j<=1;j++){
                                for(int k=-1;k<=1;k++){
                                    if(i==0 && j==0 && k==0){
                                        continue;
                                    }
                                    if(z+k>=Height || z+k<0){ // Check for z boundary
                                        continue;
                                    }
                                    dx = 0;
                                    dy = 0;
                                    if(x+i<0){ // Check for x boundary
                                        dx=Length;
                                    }
                                    else if(x+i>=Length){
                                        dx=-Length;
                                    }
                                    if(y+j<0){ // Check for y boundary
                                        dy=Width;
                                    }
                                    else if(y+j>=Width){
                                        dy=-Width;
                                    }
                                    if(lattice[getSite(x+i+dx,y+j+dy,z+k)].type==1 && lattice[getSite(x+i+dx,y+j+dy,z+k)].path_distance>0){
                                        d_temp = sqrt(pow((double)i,2)+pow((double)j,2)+pow((double)k,2)) + lattice[getSite(x+i+dx,y+j+dy,z+k)].path_distance;
                                        if(d_temp<d || d<0){
                                            d = d_temp;
                                        }
                                    }
                                }
                            }
                        }
                        if(d>0 && d<d_current+0.01){
                            lattice[getSite(x,y,z)].path_distance = d;
                            count++;
                        }
                    }
                }
            }
        }
        d_current += 1;
    }
    // Check for islands
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                if(lattice[getSite(x,y,z)].type==1 && lattice[getSite(x,y,z)].path_distance==0){
                    lattice[getSite(x,y,z)].path_distance = -1;
                }
            }
        }
    }
    double sum = 0;
    count = 0;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            if(lattice[getSite(x,y,Height-1)].path_distance>0){
                TortuosityData.push_back(lattice[getSite(x,y,Height-1)].path_distance/(Height-1));
                sum += lattice[getSite(x,y,Height-1)].path_distance/(Height-1);
                count++;
            }
        }
    }
    Path_avg = sum/count;
    sum = 0;
    count = 0;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            if(lattice[getSite(x,y,Height-1)].path_distance>0){
                sum += pow(lattice[getSite(x,y,Height-1)].path_distance/(Height-1)-Path_avg,2);
                count++;
            }
        }
    }
    Path_stdev = sqrt(sum/(count-1));
    return true;
}

void Morphology::calculateSiteEnergies(double interaction_energy,bool enable_third_neighbor_calc){
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                lattice[getSite(x,y,z)].energy = calculateSiteEnergy(x,y,z,interaction_energy,enable_third_neighbor_calc);
            }
        }
    }
}

double Morphology::calculateSiteEnergy(int x, int y, int z,double interaction_energy,bool enable_third_neighbor_calc){
    int dx,dy;
    int sum1 = 0;
    int sum2 = 0;
    int sum3 = 0;
    double energy;
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
            for(int k=-1;k<=1;k++){
                if(i==0 && j==0 && k==0){
                    continue;
                }
                if(z+k>=Height || z+k<0){ // Check for z boundary
                    continue;
                }
                dx = 0;
                dy = 0;
                if(x+i<0){ // Check for x boundary
                    dx=Length;
                }
                else if(x+i>=Length){
                    dx=-Length;
                }
                if(y+j<0){ // Check for y boundary
                    dy=Width;
                }
                else if(y+j>=Width){
                    dy=-Width;
                }
                if(lattice[getSite(x+i+dx,y+j+dy,z+k)].type!=lattice[getSite(x,y,z)].type){
                    if(abs(i)+abs(j)+abs(k)==1){
                        sum1++;
                    }
                    else if(abs(i)+abs(j)+abs(k)==2){
                        sum2++;
                    }
                    else if(abs(i)+abs(j)+abs(k)==3){
                        sum3++;
                    }
                }
            }
        }
    }
    energy = (interaction_energy/2)*sum1;
    energy += (interaction_energy/(2*sqrt((double)2)))*sum2;
    if(enable_third_neighbor_calc){
        energy += (interaction_energy/(2*sqrt((double)3)))*sum3;
    }
    return energy;
}

void Morphology::calculateTempEnergies(int x1, int y1, int z1, int di, int dj, int dk, double interaction_energy,bool enable_third_neighbor_calc){
    lattice_temp.assign(125,-1);
    int sum1,sum2,sum3,dx1,dy1,dx2,dy2;
    double energy;
    for(int x=-2;x<=2;x++){
        if((di>0 && x==-2) || (di<0 && x==2)){
            continue;
        }
        for(int y=-2;y<=2;y++){
            if((dj>0 && y==-2) || (dj<0 && y==2)){
                continue;
            }
            for(int z=-2;z<=2;z++){
                if((dk>0 && z==-2) || (dk<0 && z==2)){
                    continue;
                }
                if(z1+z>=Height || z1+z<0){ // Check for z boundary
                    lattice_temp[getTempSite(x+2,y+2,z+2)] = 0;
                    continue;
                }
                sum1 = 0;
                sum2 = 0;
                sum3 = 0;
                dx1 = 0;
                dy1 = 0;
                if(x1+x<0){ // Check for x boundary
                    dx1=Length;
                }
                else if(x1+x>=Length){
                    dx1=-Length;
                }
                if(y1+y<0){ // Check for y boundary
                    dy1=Width;
                }
                else if(y1+y>=Width){
                    dy1=-Width;
                }
                for(int i=-1;i<=1;i++){
                    for(int j=-1;j<=1;j++){
                        for(int k=-1;k<=1;k++){
                            if(z1+z+k>=Height || z1+z+k<0 || (i==0 && j==0 && k==0)){ // Check for z boundary
                                continue;
                            }
                            dx2 = 0;
                            dy2 = 0;
                            if(x1+x+dx1+i<0){ // Check for x boundary
                                dx2=Length;
                            }
                            else if(x1+x+dx1+i>=Length){
                                dx2=-Length;
                            }
                            if(y1+y+dy1+j<0){ // Check for y boundary
                                dy2=Width;
                            }
                            else if(y1+y+dy1+j>=Width){
                                dy2=-Width;
                            }
                            if(lattice[getSite(x1+x+dx1+i+dx2,y1+y+dy1+j+dy2,z1+z+k)].type!=lattice[getSite(x1+x+dx1,y1+y+dy1,z1+z)].type){
                                if((abs(i)+abs(j)+abs(k))==1){
                                    sum1++;
                                }
                                else if((abs(i)+abs(j)+abs(k))==2){
                                    sum2++;
                                }
                                else if((abs(i)+abs(j)+abs(k))==3){
                                    sum3++;
                                }
                            }
                        }
                    }
                }
                energy = (interaction_energy/2)*(double)sum1;
                energy += (interaction_energy/(2*sqrt((double)2)))*(double)sum2;
                if(enable_third_neighbor_calc){
                    energy += (interaction_energy/(2*sqrt((double)3)))*(double)sum3;
                }
                lattice_temp[getTempSite(x+2,y+2,z+2)] = energy;
            }
        }
    }
}

int Morphology::countDissimilarNeighbors(int x,int y,int z){
    int count = 0;
    int dx,dy;
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
            for(int k=-1;k<=1;k++){
                if(z+k>=Height || z+k<0){ // Check for z boundary
                    continue;
                }
                dx = 0;
                dy = 0;
                if(x+i<0){ // Check for x boundary
                    dx=Length;
                }
                else if(x+i>=Length){
                    dx=-Length;
                }
                if(y+j<0){ // Check for y boundary
                    dy=Width;
                }
                else if(y+j>=Width){
                    dy=-Width;
                }
                if(lattice[getSite(x,y,z)].type!=lattice[getSite(x+i+dx,y+j+dy,z+k)].type){
                    count++;
                }
            }
        }
    }
    return count;
}

int Morphology::countNeighbors(int z){
    int count = 0;
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
            for(int k=-1;k<=1;k++){
                if(z+k>=Height || z+k<0){ // Check for z boundary
                    continue;
                }
                count++;
            }
        }
    }
    return count;
}

void Morphology::createRandomMorphology(double mix_fraction, int seed){
    double value;
    if(mix_fraction>=1){
        cout << ProcID << ": Error creating morphology: Mix fraction must be less than one." << endl;
        return;
    }
    Energies_initialized = false;
    Mix_fraction = mix_fraction;
    gen.seed(static_cast<boost::uint32_t>(seed));
    boost::uniform_01<boost::mt19937> rand01(gen);
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                value = rand01();
                if(value<mix_fraction){
                    lattice[getSite(x,y,z)].type = 1;
                }
                else{
                    lattice[getSite(x,y,z)].type = 2;
                }
            }
        }
    }
    calculateMixFraction();
}

void Morphology::executeIsingSwappingAlt(int num_MCsteps,double interaction_energy,bool enable_third_neighbor_calc){
    // Site energy algorithm
    boost::uniform_int<> distx(0,Length-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randx(gen,distx);
    boost::uniform_int<> disty(0,Width-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randy(gen,disty);
    boost::uniform_int<> distz(0,Height-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randz(gen,distz);
    boost::uniform_int<> dist_neighbor(-1,1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_neighbor(gen,dist_neighbor);
    boost::uniform_01<boost::mt19937> rand01(gen);
    int count=0;
    int neighbor_count;
    int x,y,z,x2,y2,z2,dx,dy,temp;
    int remainder;
    int main_site;
    int neighbor_site;
    double energy_delta,probability;
    int m=1;
    vector<int> neighbors;
    neighbors.assign(6,0);
    if(!Energies_initialized){
        calculateSiteEnergies(interaction_energy,enable_third_neighbor_calc);
        Energies_initialized = true;
    }
    while(count<num_MCsteps*Length*Width*Height){
        // Randomly a site
        x = randx();
        y = randy();
        z = randz();
        // If the site is not an interfacial site, choose again
        if(!isNearInterface(x,y,z,1.1)){
            continue;
        }
        main_site = getSite(x,y,z);
        // Randomly choose a nearest neighbor site of a different type
        // First find all nearest neighbor sites of differing type
        neighbor_count = 0;
        for(int i=-1;i<=1;i++){
            for(int j=-1;j<=1;j++){
                for(int k=-1;k<=1;k++){
                    if(abs(i)+abs(j)+abs(k)>1){ // Limit to first nearest neighbors only
                        continue;
                    }
                    if(i==0 && j==0 && k==0){ // Cannot be same site
                        continue;
                    }
                    if(z+k<0 || z+k>=Height){ // Check hard z boundaries
                        continue;
                    }
                    dx = 0;
                    dy = 0;
                    if(x+i<0){ // Check for periodic x boundary
                        dx=Length;
                    }
                    else if(x+i>=Length){
                        dx=-Length;
                    }
                    if(y+j<0){ // Check for periodic y boundary
                        dy=Width;
                    }
                    else if(y+j>=Width){
                        dy=-Width;
                    }
                    if(lattice[main_site].type!=lattice[getSite(x+i+dx,y+j+dy,z+k)].type){
                        // Store site index of differing neighbor site
                        neighbors[neighbor_count] = getSite(x+i+dx,y+j+dy,z+k);
                        neighbor_count++;
                    }
                }
            }
        }
        // Randomly select one of the differing neighbor sites
        neighbor_site = neighbors[floor(rand01()*(neighbor_count-1)+0.5)];
        // Convert site index back into x,y,z coordinates
        y2 = neighbor_site/(Length*Height);
        remainder = neighbor_site % (Length*Height);
        x2 = remainder/Height;
        z2 = remainder % Height;
        // Temporarily Swap Sites
        temp = lattice[main_site].type;
        lattice[main_site].type = lattice[neighbor_site].type;
        lattice[neighbor_site].type = temp;
        // Calculate energy change and determine whether to accept the swap
        energy_delta = calculateEnergyChange1(x,y,z,x2,y2,z2,interaction_energy,enable_third_neighbor_calc);
        probability = exp(-energy_delta)/(1+exp(-energy_delta));
        if(rand01()>probability){
            // Reswap sites
            temp = lattice[main_site].type;
            lattice[main_site].type = lattice[neighbor_site].type;
            lattice[neighbor_site].type = temp;
        }
        else{
            // Put temp energies into lattice
            executeTempEnergies(x,y,z);
        }
        count++;
        if(count/(Length*Width*Height)==100*m){
            cout << ProcID << ": " << count/(Length*Width*Height) << " MC steps completed." << endl;
            m++;
        }
    }
}

void Morphology::executeIsingSwapping(int num_MCsteps,double interaction_energy,bool enable_third_neighbor_calc){
    // Bond formation algorithm
    boost::uniform_int<> distx(0,Length-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randx(gen,distx);
    boost::uniform_int<> disty(0,Width-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randy(gen,disty);
    boost::uniform_int<> distz(0,Height-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randz(gen,distz);
    boost::uniform_int<> dist_neighbor(-1,1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rand_neighbor(gen,dist_neighbor);
    boost::uniform_01<boost::mt19937> rand01(gen);
    int count = 0;
    int N = 0;
    int neighbor_count;
    int main_site;
    int neighbor_site;
    int remainder;
    int x,y,z,x2,y2,z2,dx,dy,temp;
    double energy_delta,probability;
    vector<int> neighbors;
    neighbors.assign(6,0);
    int m=1;
    while(N<num_MCsteps){
        // Randomly choose a site
        x = randx();
        y = randy();
        z = randz();
        // If site is not an interfacial site, choose again
        if(!isNearInterface(x,y,z,1.1)){
            continue;
        }
        main_site = getSite(x,y,z);
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
                    if(z+k<0 || z+k>=Height){ // Check hard z boundaries
                        continue;
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
                    if(lattice[main_site].type!=lattice[getSite(x+i+dx,y+j+dy,z+k)].type){
                        // Store site index of differing neighbor site
                        neighbors[neighbor_count] = getSite(x+i+dx,y+j+dy,z+k);
                        neighbor_count++;
                    }
                }
            }
        }
        // Randomly select one of the differing neighbor sites
        neighbor_site = neighbors[floor(rand01()*(neighbor_count-1)+0.5)];
        // Convert site index back into x,y,z coordinates
        y2 = neighbor_site/(Length*Height);
        remainder = neighbor_site % (Length*Height);
        x2 = remainder/Height;
        z2 = remainder % Height;
        // Calculate energy change and swapping probability
        energy_delta = calculateEnergyChange2(x,y,z,x2,y2,z2,interaction_energy,enable_third_neighbor_calc);
        probability = exp(-energy_delta)/(1+exp(-energy_delta));
        if(rand01()<=probability){
            // Swap Sites
            temp = lattice[main_site].type;
            lattice[main_site].type = lattice[neighbor_site].type;
            lattice[neighbor_site].type = temp;
        }
        count++;
        if(count==Length*Width*Height){
            N++;
            count = 0;
        }
        if(N==100*m){
            cout << ProcID << ": " << N << " MC steps completed." << endl;
            m++;
        }
    }
}

void Morphology::executeSmoothing(double smoothing_threshold){
    int dissimilar_neighbors;
    int num_neighbors;
    int count = 1;
    double threshold = 0.9;
    while(threshold>=smoothing_threshold){
        while(count>0){
            count = 0;
            for(int x=0;x<Length;x++){
                for(int y=0;y<Width;y++){
                    for(int z=0;z<Height;z++){
                        dissimilar_neighbors = countDissimilarNeighbors(x,y,z);
                        num_neighbors = countNeighbors(z);
                        if(((double)dissimilar_neighbors/(double)num_neighbors)>smoothing_threshold){
                            if(lattice[getSite(x,y,z)].type==1){
                                lattice[getSite(x,y,z)].type=2;
                            }
                            else if(lattice[getSite(x,y,z)].type==2){
                                lattice[getSite(x,y,z)].type=1;
                            }
                            count++;
                        }
                    }
                }
            }
        }
        threshold -= 0.01;
    }
    calculateMixFraction();
}

void Morphology::executeTempEnergies(int x, int y, int z){
    int dx,dy;
    for(int i=-2;i<=2;i++){
        for(int j=-2;j<=2;j++){
            for(int k=-2;k<=2;k++){
                if(z+k>=Height || z+k<0){ // Check for z boundary
                    continue;
                }
                dx = 0;
                dy = 0;
                if(x+i<0){ // Check for x boundary
                    dx=Length;
                }
                else if(x+i>=Length){
                    dx=-Length;
                }
                if(y+j<0){ // Check for y boundary
                    dy=Width;
                }
                else if(y+j>=Width){
                    dy=-Width;
                }
                if(lattice_temp[getTempSite(i+2,j+2,k+2)]>=0){
                    lattice[getSite(x+i+dx,y+j+dy,z+k)].energy = lattice_temp[getTempSite(i+2,j+2,k+2)];
                }
            }
        }
    }
}

vector<double> Morphology::getCorrelationData(int site_type){
    if(site_type==1){
        if(Correlation1[0]<0){
            cout << ProcID << ": Error Retrieving Correlation Data: Correlation data has not been calculated." << endl;
        }
        return Correlation1;
    }
    else if(site_type==2){
        if(Correlation2[0]<0){
            cout << ProcID << ": Error Retrieving Correlation Data: Correlation data has not been calculated." << endl;
        }
        return Correlation2;
    }
    else{
        vector<double> temp(1,0);
        cout << ProcID << ": Error Retrieving Correlation Data: Invalid site type." << endl;
        return temp;
    }
}

double Morphology::getDomainSize(int site_type){
    if(site_type==1){
        return Domain_size1;
    }
    else if(site_type==2){
        return Domain_size2;
    }
    else{
        cout << ProcID << ": Error getting domain size: Invalid site type." << endl;
        return 0;
    }
}

int Morphology::getHeight(){
    return Height;
}

int Morphology::getLength(){
    return Length;
}

double Morphology::getMixFraction(){
    return Mix_fraction;
}

double Morphology::getPathAvg(){
    return Path_avg;
}

double Morphology::getPathStdev(){
    return Path_stdev;
}

int Morphology::getRescaleSite(int x, int y, int z,int scale){
    return y*scale*Length*scale*Height+x*scale*Height+z;
}

int Morphology::getSite(int x, int y, int z){
    return y*Length*Height+x*Height+z;
}

int Morphology::getTempSite(int x, int y, int z){
    return y*5*5+x*5+z;
}

vector<double> Morphology::getTortuosityData(){
    return TortuosityData;
}

int Morphology::getSiteType(int x, int y, int z){
    return lattice[getSite(x,y,z)].type;
}

int Morphology::getWidth(){
    return Width;
}

bool Morphology::importMorphologyFile(ifstream * input){
    string var;
    int x = 0;
    int y = 0;
    int z = 0;
    int type = 0;
    int j;
    string line;
    Site site;
    getline(*input,line);
    Length = atoi(line.c_str());
    getline(*input,line);
    Width = atoi(line.c_str());
    getline(*input,line);
    Height = atoi(line.c_str());
    site.type = 0;
    site.energy = 0;
    site.path_distance = 0;
    lattice.assign(Length*Width*Height,site);
    getline(*input,line);
    Domain_size1 = atof(line.c_str());
    getline(*input,line);
    Domain_size2 = atof(line.c_str());
    getline(*input,line);
    Mix_fraction = atof(line.c_str());
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
                type = atoi(var.c_str());
                break;
            default:
                break;
            }
            j++;
        }
        site.type = type;
        site.energy = 0;
        lattice[getSite(x,y,z)] = site;
    }
    if(Domain_size1<1 || Domain_size2<1){
        calculateCorrelationDistance(15);
        cout << ProcID << ": Domain size 1 = " << Domain_size1 << endl;
        cout << ProcID << ": Domain size 2 = " << Domain_size2 << endl;
    }
    calculateMixFraction();
    return true;
}

bool Morphology::isNearInterface(int x,int y,int z,double distance){
    int d = (int)floor(distance);
    int dx,dy;
    for(int i=-d;i<=d;i++){
        for(int j=-d;j<=d;j++){
            for(int k=-d;k<=d;k++){
                if(d==1){
                    if(abs(i)+abs(j)+abs(k)>1 || i+j+k==0){
                        continue;
                    }
                }
                else if(sqrt(pow((double)i,2)+pow((double)j,2)+pow((double)k,2))>distance){
                    continue;
                }
                if(z+k>=Height || z+k<0){ // Check for z boundary
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
                if(lattice[getSite(x,y,z)].type!=lattice[getSite(x+i+dx,y+j+dy,z+k)].type){
                    return true;
                }
            }
        }
    }
    return false;
}

bool Morphology::outputMorphologyFile(ofstream * output){
    *output << Length << endl;
    *output << Width << endl;
    *output << Height << endl;
    *output << Domain_size1 << endl;
    *output << Domain_size2 << endl;
    *output << Mix_fraction << endl;
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                *output << x << "," << y << "," << z << "," << lattice[getSite(x,y,z)].type << endl;
            }
        }
    }
    return true;
}

void Morphology::rescaleLattice(int scale){
    Site site;
    site.energy = 0;
    site.path_distance = 0;
    vector<Site> lattice_rescale(Length*scale*Width*scale*Height*scale,site);
    for(int x=0;x<Length;x++){
        for(int y=0;y<Width;y++){
            for(int z=0;z<Height;z++){
                for(int i=scale*x;i<scale*x+scale;i++){
                    for(int j=scale*y;j<scale*y+scale;j++){
                        for(int k=scale*z;k<scale*z+scale;k++){
                            lattice_rescale[getRescaleSite(i,j,k,scale)].type = lattice[getSite(x,y,z)].type;
                        }
                    }
                }
            }
        }
    }
    Length = Length*scale;
    Width = Width*scale;
    Height = Height*scale;
    lattice = lattice_rescale;
}

