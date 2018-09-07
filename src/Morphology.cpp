// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "Morphology.h"

using namespace std;
using namespace Utils;
using namespace tinyxml2;

Morphology::Morphology(const int id) {
	ID = id;
	gen.seed((int)time(0)*(id + 1));
}

Morphology::Morphology(const int length, const int width, const int height, const bool enable_periodic_z, const int id) {
	ID = id;
	Lattice_Params params;
	params.Enable_periodic_x = true;
	params.Enable_periodic_y = true;
	params.Enable_periodic_z = enable_periodic_z;
	params.Length = length;
	params.Width = width;
	params.Height = height;
	params.Unit_size = 1.0;
	Enable_third_neighbor_interaction = false;
	lattice.init(params, &gen);
	gen.seed((int)time(0)*(id + 1));
}

Morphology::Morphology(const Lattice& input_lattice, const int id) {
	ID = id;
	lattice = input_lattice;
	Enable_third_neighbor_interaction = false;
	gen.seed((int)time(0)*(id + 1));
	for (int i = 0; i < (int)lattice.getNumSites(); i++) {
		bool type_found = false;
		for (int n = 0; n < (int)Site_types.size(); n++) {
			if (lattice.getSiteType(i) == Site_types[n]) {
				type_found = true;
				break;
			}
		}
		if (!type_found) {
			addSiteType(lattice.getSiteType(i));
		}
	}
}

Morphology::~Morphology() {
	//dtor
}

void Morphology::addSiteType(const char site_type) {
	// check to make sure site type has not already been added
	for (int n = 0; n < (int)Site_types.size(); n++) {
		if (Site_types[n] == site_type) {
			return;
		}
	}
	Site_types.push_back(site_type);
	Site_type_counts.push_back(0);
	Mix_fractions.push_back(-1);
	Correlation_data.resize(Correlation_data.size() + 1);
	Tortuosity_data.resize(Tortuosity_data.size() + 1);
	InterfacialHistogram_data.resize(InterfacialHistogram_data.size() + 1);
	TortuosityHistogram_data.resize(TortuosityHistogram_data.size() + 1);
	Domain_anisotropy_updated.push_back(false);
	Domain_sizes.push_back(-1);
	Domain_anisotropies.push_back(-1);
	Island_volume.push_back(-1);
}

double Morphology::calculateAdditionalEnergyChange(const long int site_index_main, const long int site_index_neighbor, const int growth_direction, const double additional_interaction) const {
	int x1, y1, z1, x2, y2, z2;
	int dx, dy, dz;
	int total_sites = 0;
	int count1_i = 0;
	int count2_i = 0;
	int count1_f = 0;
	int count2_f = 0;
	char site1_type, site2_type;
	Coords coords_main = lattice.getSiteCoords(site_index_main);
	x1 = coords_main.x;
	y1 = coords_main.y;
	z1 = coords_main.z;
	site1_type = lattice.getSiteType(coords_main);
	Coords coords_neighbor = lattice.getSiteCoords(site_index_neighbor);
	x2 = coords_neighbor.x;
	y2 = coords_neighbor.y;
	z2 = coords_neighbor.z;
	site2_type = lattice.getSiteType(coords_neighbor);
	switch (growth_direction) {
	case 1: // x-direction
		total_sites = 2;
		for (int i = -1; i <= 1; i += 2) {
			dx = lattice.calculateDX(x1, i);
			// Count the number of similar neighbors
			if (lattice.getSiteType(x1 + i + dx, y1, z1) == site1_type) {
				count1_i++;
			}
		}
		count1_f = total_sites - count1_i;
		for (int i = -1; i <= 1; i += 2) {
			dx = lattice.calculateDX(x2, i);
			// Count the number of similar neighbors
			if (lattice.getSiteType(x2 + i + dx, y2, z2) == site2_type) {
				count2_i++;
			}
		}
		count2_f = total_sites - count2_i;
		break;
	case 2: // y-direction
		total_sites = 2;
		for (int j = -1; j <= 1; j += 2) {
			dy = lattice.calculateDY(y1, j);
			// Count the number of similar neighbors
			if (lattice.getSiteType(x1, y1 + j + dy, z1) == site1_type) {
				count1_i++;
			}
		}
		count1_f = total_sites - count1_i;
		for (int j = -1; j <= 1; j += 2) {
			dy = lattice.calculateDY(y2, j);
			// Count the number of similar neighbors
			if (lattice.getSiteType(x2, y2 + j + dy, z2) == site2_type) {
				count2_i++;
			}
		}
		count2_f = total_sites - count2_i;
		break;
	case 3: // z-direction
		total_sites = 2;
		for (int k = -1; k <= 1; k += 2) {
			if (!lattice.isZPeriodic()) {
				if (z1 + k >= lattice.getHeight() || z1 + k < 0) { // Check for z boundary
					total_sites--;
					continue;
				}
			}
			dz = lattice.calculateDZ(z1, k);
			// Count the number of similar neighbors
			if (lattice.getSiteType(x1, y1, z1 + k + dz) == site1_type) {
				count1_i++;
			}
		}
		count1_f = total_sites - count1_i;
		for (int k = -1; k <= 1; k += 2) {
			if (!lattice.isZPeriodic()) {
				if (z2 + k >= lattice.getHeight() || z2 + k < 0) { // Check for z boundary
					total_sites--;
					continue;
				}
			}
			dz = lattice.calculateDZ(z2, k);
			// Count the number of similar neighbors
			if (lattice.getSiteType(x2, y2, z2 + k + dz) == site2_type) {
				count2_i++;
			}
		}
		count2_f = total_sites - count2_i;
		break;
	default:
		cout << "Error calculating the additional energy for the preferential growth direction!" << endl;
		break;
	}
	return -additional_interaction * ((count1_f - count1_i) + (count2_f - count2_i));
}

void Morphology::calculateAnisotropies(const int N_sampling_max) {
	cout << ID << ": Calculating the domain anisotropy..." << endl;
	// Select sites for correlation function calculation.
	// Site indices for each selected site are stored in the Correlation_sites vector.
	vector<vector<long int>> correlation_sites_data(Site_types.size());
	for (int n = 0; n < (int)Site_types.size(); n++) {
		if ((int)correlation_sites_data[n].size() == 0) {
			// Only N_sampling_max sites are randomly selected
			getSiteSampling(correlation_sites_data[n], Site_types[n], N_sampling_max);
		}
	}
	Domain_anisotropy_updated.assign(Site_types.size(), false);
	bool success = false;
	int cutoff_distance = 5;
	while (!success) {
		if (2 * cutoff_distance > lattice.getLength() && 2 * cutoff_distance > lattice.getWidth() && 2 * cutoff_distance > lattice.getHeight()) {
			success = false;
			break;
		}
		for (int i = 0; i < (int)Site_types.size(); i++) {
			// Only perform anisotropy calculations for site types that have not yet been updated and that have at least 100 site counts.
			if (!Domain_anisotropy_updated[i] && Site_type_counts[i] > 100) {
				cout << ID << ": Performing sampling anisotropy calculation with " << (int)correlation_sites_data[i].size() << " sites for site type " << (int)Site_types[i] << " with a cutoff of " << cutoff_distance << "..." << endl;
				Domain_anisotropy_updated[i] = calculateAnisotropy(correlation_sites_data[i], Site_types[i], cutoff_distance);
			}
		}
		// Check if the anisotropy has been successfully calculated for all sites types
		for (int i = 0; i < (int)Site_types.size(); i++) {
			// If the anisotropy was not updated for any of the site types, increment the cutoff_distance, break the for loop, and repeat the anisotropy calculations
			if (!Domain_anisotropy_updated[i] && Site_type_counts[i] > 100) {
				success = false;
				cutoff_distance++;
				break;
			}
			// If the i index reaches its max value here, then the anisotropy calculation was successful for all site types.
			if (i == (int)Site_types.size() - 1) {
				success = true;
			}
		}
	}
	if (!success) {
		cout << ID << ": Warning! Could not calculate the domain anisotropy." << endl;
	}
}

bool Morphology::calculateAnisotropy(const vector<long int>& correlation_sites, const char site_type, const int cutoff_distance) {
	int type_index = getSiteTypeIndex(site_type);
	int N_sites = 0;
	double correlation_length_x = 0;
	double correlation_length_y = 0;
	double correlation_length_z = 0;
	double d1, y1, y2, slope, intercept;
	Coords site_coords, coords_dest;
	vector<double> correlation_x(cutoff_distance + 1, 0.0);
	vector<double> correlation_y(cutoff_distance + 1, 0.0);
	vector<double> correlation_z(cutoff_distance + 1, 0.0);
	vector<int> site_count(cutoff_distance + 1, 0);
	vector<int> site_total(cutoff_distance + 1, 0);
	// Check that correlation sites vector is not empty
	if (!(correlation_sites.size() > 0)) {
		cout << ID << ": Error! Vector of site tags to be used in the anisotropy calculation is empty." << endl;
		return false;
	}
	for (int m = 0; m < (int)correlation_sites.size(); m++) {
		if (lattice.getSiteType(correlation_sites[m]) != site_type) {
			continue;
		}
		site_coords = lattice.getSiteCoords(correlation_sites[m]);
		// Calculate correlation length in the x-direction
		site_count.assign(cutoff_distance + 1, 0);
		for (int i = -cutoff_distance; i <= cutoff_distance; i++) {
			if (!lattice.checkMoveValidity(site_coords, i, 0, 0)) {
				continue;
			}
			lattice.calculateDestinationCoords(site_coords, i, 0, 0, coords_dest);
			if (lattice.getSiteType(site_coords) == lattice.getSiteType(coords_dest)) {
				site_count[abs(i)]++;
			}
		}
		for (int n = 1; n <= cutoff_distance; n++) {
			correlation_x[n] += (double)site_count[n] / 2;
		}
		// Calculate correlation length in the y-direction
		site_count.assign(cutoff_distance + 1, 0);
		for (int j = -cutoff_distance; j <= cutoff_distance; j++) {
			if (!lattice.checkMoveValidity(site_coords, 0, j, 0)) {
				continue;
			}
			lattice.calculateDestinationCoords(site_coords, 0, j, 0, coords_dest);
			if (lattice.getSiteType(site_coords) == lattice.getSiteType(coords_dest)) {
				site_count[abs(j)]++;
			}
		}
		for (int n = 1; n <= cutoff_distance; n++) {
			correlation_y[n] += (double)site_count[n] / 2;
		}
		// Calculate correlation length in the z-direction
		site_total.assign(cutoff_distance + 1, 0);
		site_count.assign(cutoff_distance + 1, 0);
		for (int k = -cutoff_distance; k <= cutoff_distance; k++) {
			if (!lattice.checkMoveValidity(site_coords, 0, 0, k)) {
				continue;
			}
			lattice.calculateDestinationCoords(site_coords, 0, 0, k, coords_dest);
			if (lattice.getSiteType(site_coords) == lattice.getSiteType(coords_dest)) {
				site_count[abs(k)]++;
			}
			site_total[abs(k)]++;
		}
		for (int n = 1; n <= cutoff_distance; n++) {
			if (site_total[n] > 0) {
				correlation_z[n] += (double)site_count[n] / site_total[n];
			}
		}
		N_sites++;
	}
	// Average correlation data over all starting sites and normalize
	double averaging = 1.0 / N_sites;
	double norm = 1.0 / (1.0 - Mix_fractions[type_index]);
	correlation_x[0] = 1.0;
	correlation_y[0] = 1.0;
	correlation_z[0] = 1.0;
	for (int n = 1; n <= cutoff_distance; n++) {
		correlation_x[n] *= averaging;
		correlation_x[n] -= Mix_fractions[type_index];
		correlation_x[n] *= norm;
		correlation_y[n] *= averaging;
		correlation_y[n] -= Mix_fractions[type_index];
		correlation_y[n] *= norm;
		correlation_z[n] *= averaging;
		correlation_z[n] -= Mix_fractions[type_index];
		correlation_z[n] *= norm;

	}
	// Find the bounds of where the pair-pair correlation functions reach 1/e
	bool success_x = false;
	bool success_y = false;
	bool success_z = false;
	for (int n = 2; n <= cutoff_distance; n++) {
		if (!success_x && correlation_x[n] < (1.0 / exp(1.0))) {
			d1 = n - 1;
			y1 = correlation_x[n - 1];
			y2 = correlation_x[n];
			// Use linear interpolation to determine the cross-over point
			slope = (y2 - y1);
			intercept = y1 - slope * d1;
			correlation_length_x = 2.0 * (1.0 / exp(1.0) - intercept) / slope;
			success_x = true;
		}
		if (!success_y && correlation_y[n] < (1.0 / exp(1.0))) {
			d1 = n - 1;
			y1 = correlation_y[n - 1];
			y2 = correlation_y[n];
			// Use linear interpolation to determine the cross-over point
			slope = (y2 - y1);
			intercept = y1 - slope * d1;
			correlation_length_y = 2.0 * (1.0 / exp(1.0) - intercept) / slope;
			success_y = true;
		}
		if (!success_z && correlation_z[n] < (1.0 / exp(1.0))) {
			d1 = n - 1;
			y1 = correlation_z[n - 1];
			y2 = correlation_z[n];
			// Use linear interpolation to determine the cross-over point
			slope = (y2 - y1);
			intercept = y1 - slope * d1;
			correlation_length_z = 2.0 * (1.0 / exp(1.0) - intercept) / slope;
			success_z = true;
		}
		else if (!success_z && cutoff_distance > lattice.getHeight()) {
			correlation_length_z = lattice.getHeight();
			success_z = true;
		}
		// If cross-over point is not reached, return the point where the first minimum is reached.
		if (!success_x && correlation_x[n] > correlation_x[n - 1]) {
			correlation_length_x = n - 1;
			success_x = true;
		}
		if (!success_y && correlation_y[n] > correlation_y[n - 1]) {
			correlation_length_y = n - 1;
			success_y = true;
		}
		if (!success_z && correlation_z[n] > correlation_z[n - 1]) {
			correlation_length_z = n - 1;
			success_z = true;
		}
	}
	if (!success_x || !success_y || !success_z) {
		cout << ID << ": Cutoff distance of " << cutoff_distance << " is too small to calculate anisotropy of domain type " << (int)site_type << "." << endl;
		Domain_anisotropies[type_index] = -1;
		return false;
	}
	if (4 * correlation_length_x > lattice.getLength() || 4 * correlation_length_y > lattice.getWidth()) {
		cout << "Warning.  Correlation length in x- or y-direction is greater than L/4." << endl;
		cout << "x-direction correlation length is " << correlation_length_x << "." << endl;
		cout << "y-direction correlation length is " << correlation_length_y << "." << endl;
	}
	Domain_anisotropies[type_index] = (2 * correlation_length_z) / (correlation_length_x + correlation_length_y);
	return true;
}

double Morphology::calculateCorrelationDistance(const vector<long int>& correlation_sites, vector<double>& correlation_data, const double mix_fraction, const int cutoff_distance, const CorrelationCalc_Params& params) {
	vector<int> site_count, total_count;
	double distance;
	int bin;
	double d1, y1, y2, slope, intercept;
	Coords site_coords, coords_dest;
	if (cutoff_distance > lattice.getLength() || cutoff_distance > lattice.getWidth()) {
		cout << ID << ": Error, cutoff distance is greater than the lattice length and/or width." << endl;
		return -1;
	}
	// Resolution of correlation distance data is 0.5 lattice units
	int correlation_size_old = (int)correlation_data.size();
	int correlation_size_new = 2 * cutoff_distance + 1;
	if (correlation_size_old >= correlation_size_new) {
		cout << ID << ": Error, new cutoff distance is not greater than the previous cutoff distance and no new calculations have been performed." << endl;
		return -1;
	}
	// Initialize vectors to store correlation function data
	for (int m = 0; m < (correlation_size_new - correlation_size_old); m++) {
		correlation_data.push_back(0);
	}
	// Loop through all selected sites and determine the correlation function for each
	// The pair-pair correlation is determined based on the fraction of sites that are the same as the starting site and this function is calculated as a function of distance from the starting site.
	// Sites surrounding the start site are placed into bins based on their distance from the starting site.
	// Bins covering a distance range of half a lattice unit are used, ex: second bin is from 0.25a to 0.7499a, third bin is from 0.75a to 1.2499a, etc.
	// site_count vector stores the number of sites that are are the same type as the starting site for each bin
	// total_count vector stores the total number of sites in each bin
	site_count.assign(2 * cutoff_distance + 1, 0);
	total_count.assign(2 * cutoff_distance + 1, 0);
	for (int m = 0; m < (int)correlation_sites.size(); m++) {
		site_coords = lattice.getSiteCoords(correlation_sites[m]);
		fill(site_count.begin(), site_count.end(), 0);
		fill(total_count.begin(), total_count.end(), 0);
		for (int i = -cutoff_distance; i <= cutoff_distance; i++) {
			for (int j = -cutoff_distance; j <= cutoff_distance; j++) {
				for (int k = -cutoff_distance; k <= cutoff_distance; k++) {
					// The distance between two sites is rounded to the nearest half a lattice unit
					bin = round_int(2.0 * sqrt(i*i + j * j + k * k));
					// Calculation is skipped for bin values that have already been calculated during previous calls to the calculateCorrelationDistance function
					if (bin < (correlation_size_old - 1)) {
						continue;
					}
					distance = (double)bin / 2.0;
					if (distance > cutoff_distance) {
						continue;
					}
					if (!lattice.checkMoveValidity(site_coords, i, j, k)) {
						continue;
					}
					lattice.calculateDestinationCoords(site_coords, i, j, k, coords_dest);
					if (lattice.getSiteType(site_coords) == lattice.getSiteType(coords_dest)) {
						site_count[bin]++;
					}
					total_count[bin]++;
				}
			}
		}
		//  Calculate the fraction of similar sites for each bin
		for (int n = 0; n < (int)correlation_data.size(); n++) {
			if (n < correlation_size_old) {
				continue;
			}
			if (total_count[n] > 0) {
				correlation_data[n] += (double)site_count[n] / (double)total_count[n];
			}
			else {
				correlation_data[n] += 1;
			}
		}
	}
	// Average overall starting sites and normalize the correlation data
	double averaging = 1.0 / correlation_sites.size();
	double norm = 1.0 / (1.0 - mix_fraction);
	for (int n = 0; n < (int)correlation_data.size(); n++) {
		if (n < correlation_size_old) {
			continue;
		}
		correlation_data[n] *= averaging;
		correlation_data[n] -= mix_fraction;
		correlation_data[n] *= norm;
	}
	// Find the bounds of where the pair-pair correlation function first crosses over the Mix_fraction
	if (params.Enable_mix_frac_method) {
		for (int n = 2; n < (int)correlation_data.size(); n++) {
			if (correlation_data[n] < 0) {
				d1 = (double)(n - 1) * 0.5;
				y1 = correlation_data[n - 1];
				y2 = correlation_data[n];
				// Use linear interpolation to determine the cross-over point
				slope = (y2 - y1) * 2.0;
				intercept = y1 - slope * d1;
				return -intercept / slope;
			}
			if (correlation_data[n] > correlation_data[n - 1]) {
				return (double)(n - 1) / 2.0;
			}
		}
	}
	// Find the bounds of where the pair-pair correlation function first reaches within 1/e of the Mix_fraction
	if (params.Enable_e_method) {
		for (int n = 2; n < (int)correlation_data.size(); n++) {
			if (correlation_data[n] < (1.0 / exp(1.0))) {
				d1 = (double)(n - 1) * 0.5;
				y1 = correlation_data[n - 1];
				y2 = correlation_data[n];
				// Use linear interpolation to determine the cross-over point
				slope = (y2 - y1) * 2.0;
				intercept = y1 - slope * d1;
				return 2.0 * (1.0 / exp(1.0) - intercept) / slope;
			}
		}
	}
	return -1;
}

void Morphology::calculateCorrelationDistances(const CorrelationCalc_Params& params) {
	if (params.Enable_extended_correlation_calc) {
		cout << ID << ": Calculating the domain size using the extended pair-pair correlation function using a cutoff radius of " << params.Correlation_cutoff_distance << "..." << endl;
	}
	if (params.Enable_mix_frac_method) {
		cout << ID << ": Calculating the domain size from the pair-pair correlation function using the mix fraction method..." << endl;
	}
	else if (params.Enable_e_method) {
		cout << ID << ": Calculating the domain size from the pair-pair correlation function using the 1/e method..." << endl;
	}
	vector<vector<long int>> correlation_sites_data(Site_types.size());
	// Select sites for correlation function calculation.
	// Site indices for each selected site are stored in the Correlation_sites vector.
	for (int n = 0; n < (int)Site_types.size(); n++) {
		if ((int)correlation_sites_data[n].size() == 0) {
			getSiteSampling(correlation_sites_data[n], Site_types[n], params.N_sampling_max);
		}
	}
	vector<bool> domain_size_updated(Site_types.size(), false);
	int cutoff_distance;
	double domain_size;
	for (int n = 0; n < (int)Site_types.size(); n++) {
		if (params.Enable_extended_correlation_calc) {
			cutoff_distance = params.Correlation_cutoff_distance;
		}
		else {
			cutoff_distance = 5;
		}
		domain_size = -1;
		// The correlation function calculation is called with an increasing cutoff distance until successful.
		while (!domain_size_updated[n]) {
			if (2 * cutoff_distance > lattice.getLength() || 2 * cutoff_distance > lattice.getWidth() || (lattice.isZPeriodic() && 2 * cutoff_distance > lattice.getHeight())) {
				cout << ID << ": Correlation calculation cutoff radius is now too large to continue accurately calculating the correlation function for site type " << (int)Site_types[n] << "." << endl;
				break;
			}
			if (Site_type_counts[n] > 100) {
				cout << ID << ": Performing sampling domain size calculation with " << (int)correlation_sites_data[n].size() << " sites for site type " << (int)Site_types[n] << " with a cutoff radius of " << cutoff_distance << "..." << endl;
				domain_size = calculateCorrelationDistance(correlation_sites_data[n], Correlation_data[n], Mix_fractions[n], cutoff_distance, params);
			}
			if (domain_size > 0) {
				domain_size_updated[n] = true;
				Domain_sizes[n] = domain_size;
			}
			else {
				cout << ID << ": Cutoff distance of " << cutoff_distance << " is too small to calculate the size of domain type " << (int)Site_types[n] << "." << endl;
				cutoff_distance++;
			}
		}
	}
}

void Morphology::calculateDepthDependentData(const CorrelationCalc_Params& correlation_calc_params) {
	CorrelationCalc_Params correlation_params = correlation_calc_params;
	correlation_params.Enable_extended_correlation_calc = false;
	if (correlation_params.Enable_mix_frac_method) {
		cout << ID << ": Calculating the depth dependent domain size from the pair-pair correlation function using the mix fraction method..." << endl;
	}
	else if (correlation_params.Enable_e_method) {
		cout << ID << ": Calculating the depth dependent domain size from the pair-pair correlation function using the 1/e method..." << endl;
	}
	vector<vector<long int>> correlation_sites_data(Site_types.size());
	vector<double> correlation_data;
	vector<bool> domain_size_updated(Site_types.size(), false);
	Depth_composition_data.assign(Site_types.size(), vector<double>(lattice.getHeight(), 0));
	Depth_iv_data.assign(Site_types.size(), vector<double>(lattice.getHeight(), 0));
	Depth_domain_size_data.assign(Site_types.size(), vector<double>(lattice.getHeight(), 0));
	vector<int> counts(Site_types.size(), 0);
	int site_count;
	double mix_fraction_local;
	Coords coords, coords_dest;
	int counts_total = 0;
	int cutoff_distance;
	double domain_size;
	for (int z = 0; z < lattice.getHeight(); z++) {
		//cout << ID << ": Calculating depth dependent data for z = " << z << endl;
		// Calculate depth dependent composition
		counts.assign(Site_types.size(), 0);
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				counts[getSiteTypeIndex(lattice.getSiteType(x, y, z))]++;
			}
		}
		for (int n = 0; n < (int)Site_types.size(); n++) {
			Depth_composition_data[n][z] = (double)counts[n] / (double)(lattice.getLength()*lattice.getWidth());
		}
		// Calculate depth dependent interfacial volume fraction
		counts.assign(Site_types.size(), 0);
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				coords.setXYZ(x, y, z);
				if (isNearInterface(coords, 1)) {
					counts[getSiteTypeIndex(lattice.getSiteType(x, y, z))]++;
				}
			}
		}
		for (int n = 0; n < (int)Site_types.size(); n++) {
			Depth_iv_data[n][z] = (double)counts[n] / (double)(lattice.getLength()*lattice.getWidth());
		}
		// Select sites for correlation function calculation.
		correlation_sites_data.assign(Site_types.size(), vector<long int>(0));
		for (int n = 0; n < (int)Site_types.size(); n++) {
			if ((int)correlation_sites_data[n].size() == 0) {
				getSiteSamplingZ(correlation_sites_data[n], Site_types[n], correlation_params.N_sampling_max, z);
			}
		}
		// Calculate depth dependent domain size
		domain_size_updated.assign(Site_types.size(), false);
		for (int n = 0; n < (int)Site_types.size(); n++) {
			cutoff_distance = 5;
			correlation_data.clear();
			domain_size = -1;
			// The correlation function calculation is called with an increasing cutoff distance until successful.
			while (!domain_size_updated[n]) {
				if (2 * cutoff_distance > lattice.getLength() || 2 * cutoff_distance > lattice.getWidth() || 2 * cutoff_distance > lattice.getHeight()) {
					//cout << ID << ": Correlation calculation cutoff radius is now too large to continue accurately calculating the correlation function for site type " << (int)Site_types[n] << "." << endl;
					break;
				}
				// Calculate local mix fraction
				site_count = 0;
				counts_total = 0;
				for (int x = 0; x < lattice.getLength(); x++) {
					for (int y = 0; y < lattice.getWidth(); y++) {
						for (int k = -cutoff_distance; k <= cutoff_distance; k++) {
							coords.setXYZ(x, y, z);
							if (!lattice.checkMoveValidity(coords, 0, 0, k)) {
								continue;
							}
							lattice.calculateDestinationCoords(coords, 0, 0, k, coords_dest);
							if (lattice.getSiteType(coords_dest) == Site_types[n]) {
								site_count++;
							}
							counts_total++;
						}
					}
				}
				mix_fraction_local = (double)site_count / (double)counts_total;
				if (Site_type_counts[n] > 10) {
					//cout << ID << ": Performing sampling domain size calculation with " << (int)correlation_sites_data[n].size() << " sites for site type " << (int)Site_types[n] << " with a cutoff radius of " << cutoff_distance << "..." << endl;
					domain_size = calculateCorrelationDistance(correlation_sites_data[n], correlation_data, mix_fraction_local, cutoff_distance, correlation_params);
				}
				if (domain_size > 0) {
					Depth_domain_size_data[n][z] = domain_size;
					domain_size_updated[n] = true;
				}
				else {
					cutoff_distance++;
				}
			}
		}
		// Calculate depth dependent anisotropy
	}
}

double Morphology::calculateDissimilarFraction(const Coords& coords, const int rescale_factor) const {
	int site_count = 0;
	int count_dissimilar = 0;
	Coords coords_dest;
	// When the rescale factor is 1, the radius is 1, and the radius increases for larger rescale factors.
	static int radius = (rescale_factor <= 2) ? (int)ceil(3.0 / 2.0) : (int)ceil((double)(rescale_factor + 1) / 2);
	static int cutoff_squared = (rescale_factor <= 2) ? (int)floor((3.0 / 2.0)*(3.0 / 2.0)) : (int)floor(((double)(rescale_factor + 1) / 2)*((double)(rescale_factor + 1) / 2));
	for (int i = -radius; i <= radius; i++) {
		for (int j = -radius; j <= radius; j++) {
			for (int k = -radius; k <= radius; k++) {
				if ((i*i + j * j + k * k) > cutoff_squared) {
					continue;
				}
				if (!lattice.checkMoveValidity(coords, i, j, k)) {
					continue;
				}
				lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
				if (lattice.getSiteType(coords) != lattice.getSiteType(coords_dest)) {
					count_dissimilar++;
				}
				site_count++;
			}
		}
	}
	return (double)count_dissimilar / (double)site_count;
}

double Morphology::calculateEnergyChangeSimple(const long int site_index1, const long int site_index2, const double interaction_energy1, const double interaction_energy2) {
	// Used with bond formation algorithm
	static const double one_over_sqrt2 = 1 / sqrt(2);
	static const double one_over_sqrt3 = 1 / sqrt(3);
	char sum1_1_delta, sum2_1_delta, sum3_1_delta, sum1_2_delta, sum2_2_delta, sum3_2_delta;
	double sum_1_delta, sum_2_delta;
	char site1_type = lattice.getSiteType(site_index1);
	// Calculate change around site 1
	char sum1_1i = Neighbor_counts[site_index1].sum1;
	char sum2_1i = Neighbor_counts[site_index1].sum2;
	char sum3_1i = Neighbor_counts[site_index1].sum3;
	char sum1_2f = Neighbor_info[site_index1].total1 - sum1_1i - 1;
	char sum2_2f = Neighbor_info[site_index1].total2 - sum2_1i;
	char sum3_2f = Neighbor_info[site_index1].total3 - sum3_1i;
	// Calculate change around site 2
	char sum1_2i = Neighbor_counts[site_index2].sum1;
	char sum2_2i = Neighbor_counts[site_index2].sum2;
	char sum3_2i = Neighbor_counts[site_index2].sum3;
	char sum1_1f = Neighbor_info[site_index2].total1 - sum1_2i - 1;
	char sum2_1f = Neighbor_info[site_index2].total2 - sum2_2i;
	char sum3_1f = Neighbor_info[site_index2].total3 - sum3_2i;
	// Save swapped state into temp_counts1 and temp_counts2
	Temp_counts1.sum1 = sum1_2f;
	Temp_counts1.sum2 = sum2_2f;
	Temp_counts1.sum3 = sum3_2f;
	Temp_counts2.sum1 = sum1_1f;
	Temp_counts2.sum2 = sum2_1f;
	Temp_counts2.sum3 = sum3_1f;
	// Calculate change
	sum1_1_delta = sum1_1f - sum1_1i;
	sum2_1_delta = sum2_1f - sum2_1i;
	sum3_1_delta = sum3_1f - sum3_1i;
	sum1_2_delta = sum1_2f - sum1_2i;
	sum2_2_delta = sum2_2f - sum2_2i;
	sum3_2_delta = sum3_2f - sum3_2i;
	sum_1_delta = -(double)sum1_1_delta - (double)sum2_1_delta*one_over_sqrt2;
	sum_2_delta = -(double)sum1_2_delta - (double)sum2_2_delta*one_over_sqrt2;
	// By default interactions with the third-nearest neighbors are not included, but when enabled they are added here
	if (Enable_third_neighbor_interaction) {
		sum_1_delta -= (double)sum3_1_delta*one_over_sqrt3;
		sum_2_delta -= (double)sum3_2_delta*one_over_sqrt3;
	}
	if (site1_type == 1) {
		return interaction_energy1 * sum_1_delta + interaction_energy2 * sum_2_delta;
	}
	else {
		return interaction_energy2 * sum_1_delta + interaction_energy1 * sum_2_delta;
	}
}

double Morphology::calculateEnergyChange(const Coords& coords1, const Coords& coords2, const double interaction_energy1, const double interaction_energy2) const {
	// Used with bond formation algorithm
	char site1_type, site2_type;
	int sum1_1_delta, sum2_1_delta, sum3_1_delta, sum1_2_delta, sum2_2_delta, sum3_2_delta;
	double sum_1_delta, sum_2_delta;
	static const double one_over_sqrt2 = 1 / sqrt(2);
	static const double one_over_sqrt3 = 1 / sqrt(3);
	Coords coords_dest;
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
	site1_type = lattice.getSiteType(coords1);
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			for (int k = -1; k <= 1; k++) {
				if (!lattice.checkMoveValidity(coords1, i, j, k)) {
					// Total site counts must be reduced if next to a hard boundary
					switch (i*i + j * j + k * k) {
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
				lattice.calculateDestinationCoords(coords1, i, j, k, coords_dest);
				// Count the number of similar neighbors
				if (lattice.getSiteType(coords_dest) == site1_type) {
					switch (i*i + j * j + k * k) {
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
	sum1_2f = total1 - sum1_1i - 1;
	sum2_2f = total2 - sum2_1i;
	sum3_2f = total3 - sum3_1i;
	// Calculate change around x2,y2,z2
	site2_type = lattice.getSiteType(coords2);
	// There are in total 6 first-nearest, 12 second-nearest, and 8 third-nearest neighbors
	total1 = 6;
	total2 = 12;
	total3 = 8;
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			for (int k = -1; k <= 1; k++) {
				if (!lattice.checkMoveValidity(coords2, i, j, k)) {
					switch (i*i + j * j + k * k) {
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
				lattice.calculateDestinationCoords(coords2, i, j, k, coords_dest);
				// Count the number of similar neighbors
				if (lattice.getSiteType(coords_dest) == site2_type) {
					switch (i*i + j * j + k * k) {
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
	sum1_1f = total1 - sum1_2i - 1;
	sum2_1f = total2 - sum2_2i;
	sum3_1f = total3 - sum3_2i;
	sum1_1_delta = sum1_1f - sum1_1i;
	sum2_1_delta = sum2_1f - sum2_1i;
	sum3_1_delta = sum3_1f - sum3_1i;
	sum1_2_delta = sum1_2f - sum1_2i;
	sum2_2_delta = sum2_2f - sum2_2i;
	sum3_2_delta = sum3_2f - sum3_2i;
	sum_1_delta = -(double)sum1_1_delta - (double)sum2_1_delta*one_over_sqrt2;
	sum_2_delta = -(double)sum1_2_delta - (double)sum2_2_delta*one_over_sqrt2;
	// By default interactions with the third-nearest neighbors are not included, but when enabled they are added here
	if (Enable_third_neighbor_interaction) {
		sum_1_delta -= (double)sum3_1_delta*one_over_sqrt3;
		sum_2_delta -= (double)sum3_2_delta*one_over_sqrt3;
	}
	if (site1_type == (char)1) {
		return interaction_energy1 * sum_1_delta + interaction_energy2 * sum_2_delta;
	}
	else {
		return interaction_energy2 * sum_1_delta + interaction_energy1 * sum_2_delta;
	}
}

double Morphology::calculateInterfacialAreaVolumeRatio() const {
	unsigned long site_count = 0;
	Coords coords, coords_dest;
	for (int m = 0; m < (int)Site_types.size() - 1; m++) {
		for (int n = m + 1; n < (int)Site_types.size(); n++) {
			for (int x = 0; x < lattice.getLength(); x++) {
				for (int y = 0; y < lattice.getWidth(); y++) {
					for (int z = 0; z < lattice.getHeight(); z++) {
						coords.setXYZ(x, y, z);
						if (lattice.getSiteType(coords) == Site_types[m]) {
							for (int i = -1; i <= 1; i++) {
								for (int j = -1; j <= 1; j++) {
									for (int k = -1; k <= 1; k++) {
										if (abs(i) + abs(j) + abs(k) > 1) {
											continue;
										}
										if (!lattice.checkMoveValidity(coords, i, j, k)) {
											continue;
										}
										lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
										if (lattice.getSiteType(coords_dest) == Site_types[n]) {
											site_count++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return (double)site_count / (double)lattice.getNumSites();
}

void Morphology::calculateInterfacialDistance() {
	Coords coords, coords_dest;
	int d_int;
	float d;
	float d_temp;
	// The shortest distance from each site to the interface is stored in the path_distances vector
	vector<float> path_distances;
	path_distances.assign(lattice.getNumSites(), 0);
	// Calculate distances to the interface by expanding outward from the interface
	// A first scan over the lattice is done to identify the sites at the interface, at a distance of less than 2 lattice units from the interface
	// Subsequent scans expand outward 1 lattice unit at a time from these interfacial sites until no sites are left.
	int calc_count = 1;
	float d_current = (float)1.99;
	while (calc_count > 0) {
		calc_count = 0;
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					coords.setXYZ(x, y, z);
					// Only perform the calculation on sites with a yet unknown interfacial distance
					if (path_distances[lattice.getSiteIndex(x, y, z)] < 0.1) {
						d = -1;
						// Look around at neighboring sites
						for (int i = -1; i <= 1; i++) {
							for (int j = -1; j <= 1; j++) {
								for (int k = -1; k <= 1; k++) {
									if (!lattice.checkMoveValidity(coords, i, j, k)) {
										continue;
									}
									lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
									// Initial scan identifies interfacial sites
									if (d_current < 2) {
										if (lattice.getSiteType(coords_dest) != lattice.getSiteType(x, y, z)) {
											d_temp = sqrt((float)(i*i + j * j + k * k));
											if (d < 0 || d_temp < d) {
												d = d_temp;
											}
										}
									}
									// Subsequent scans identify sites of the same type
									// A temporary distance to the interface from the target site by way of the identified neighbor site is calculated
									// The temporary distance is only stored if it is shorter than the previous temporary distance, ensuring that the shortest interfacial distance is calculated
									else if (lattice.getSiteType(coords_dest) == lattice.getSiteType(x, y, z) && path_distances[lattice.getSiteIndex(coords_dest)] > 0.1) {
										d_temp = path_distances[lattice.getSiteIndex(coords_dest)] + sqrt((float)(i*i + j * j + k * k));
										if (d < 0 || d_temp < d) {
											d = d_temp;
										}
									}
								}
							}
						}
						// The temporary distance is only accepted if it less than the expansion limit (d_current).
						if (d > 0 && d < d_current) {
							path_distances[lattice.getSiteIndex(x, y, z)] = d;
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
	for (int i = 0; i < (int)Site_types.size(); i++) {
		InterfacialHistogram_data[i].clear();
	}
	vector<int> counts((int)Site_types.size(), 0);
	// The interfacial distance data in path_data is rounded to the nearest integer lattice unit
	// One bin for each integer lattice unit is used to create the histograms.
	// Set site of histogram vectors
	for (int m = 0; m < (int)path_distances.size(); m++) {
		d_int = round_int(path_distances[m]);
		for (int n = 0; n < (int)Site_types.size(); n++) {
			if (lattice.getSiteType(m) == Site_types[n]) {
				if (d_int > (int)InterfacialHistogram_data[n].size()) {
					InterfacialHistogram_data[n].assign(d_int,0.0);
				}
			}
		}
	}
	// Collect counts in each histogram bin
	for (int m = 0; m < (int)path_distances.size(); m++) {
		d_int = round_int(path_distances[m]);
		for (int n = 0; n < (int)Site_types.size(); n++) {
			if (lattice.getSiteType(m) == Site_types[n]) {
				InterfacialHistogram_data[n][d_int - 1] += 1.0;
				counts[n]++;
			}
		}
	}
	for (int n = 0; n < (int)Site_types.size(); n++) {
		for (int i = 0; i < (int)InterfacialHistogram_data[n].size(); i++) {
			InterfacialHistogram_data[n][i] /= (double)counts[n];
		}
	}
}

double Morphology::calculateInterfacialVolumeFraction() const {
	unsigned long site_count = 0;
	Coords coords, coords_dest;
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				coords.setXYZ(x, y, z);
				// For each site in the lattice, the neighboring sites are checked to see if there is one that is not the same type.
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						for (int k = -1; k <= 1; k++) {
							if (!lattice.checkMoveValidity(coords, i, j, k)) {
								continue;
							}
							lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
							if (lattice.getSiteType(x, y, z) != lattice.getSiteType(coords_dest)) {
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
	return (double)site_count / (double)lattice.getNumSites();
}

void Morphology::calculateMixFractions() {
	//Calculate final Mix_fraction
	vector<int> counts((int)Site_types.size(), 0);
	int type_index;
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				type_index = getSiteTypeIndex(lattice.getSiteType(x, y, z));
				counts[type_index]++;
			}
		}
	}
	for (int i = 0; i < (int)Site_types.size(); i++) {
		Mix_fractions[i] = (double)counts[i] / (double)lattice.getNumSites();
	}
}

Morphology::NeighborCounts Morphology::calculateNeighborCounts(const Coords& coords) const {
	Coords coords_dest;
	NeighborCounts counts;
	counts.sum1 = 0;
	counts.sum2 = 0;
	counts.sum3 = 0;
	// Calculate similar neighbors around x,y,z
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			for (int k = -1; k <= 1; k++) {
				if (!lattice.checkMoveValidity(coords, i, j, k)) {
					continue;
				}
				lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
				// Count the number of similar neighbors
				if (lattice.getSiteType(coords_dest) == lattice.getSiteType(coords)) {
					switch (i*i + j * j + k * k) {
					case 1:
						counts.sum1++;
						break;
					case 2:
						counts.sum2++;
						break;
					case 3:
						counts.sum3++;
						break;
					default:
						break;
					}
				}
			}
		}
	}
	return counts;
}

bool Morphology::calculatePathDistances(vector<float>& path_distances) {
	int z;
	Coords coords;
	long int current_index;
	long int neighbor_index;
	float d;
	float d_temp;
	const static float sqrt_two = sqrt((float)2.0);
	const static float sqrt_three = sqrt((float)3.0);
	// Create and initialize a blank node.
	// Each node contains a vector with indices of all first- ,second-, and third-nearest neighbors (at most 26 neighbors).
	// Another vector stores the squared distance to each of the neighbors.
	// Each node also has an estimated distance from the destination.
	Node temp_node;
	// Create a node vector that is the same size as the lattice and initialize with blank nodes.
	vector<Node> Node_vector;
	Node_vector.assign(lattice.getNumSites(), temp_node);
	// The neighbor_nodes set is sorted by the estimated distance of nodes in the set.
	// This set is used in Dijsktra's algorithm to keep a sorted list of all nodes that are neighboring nodes that have their path distances already determined.
	// Once the path distance for a particular node is fixed, all of its neighboring nodes that have not yet been fixed are added to the neighbor_nodes set.
	set<vector<Node>::iterator, NodeIteratorCompare> neighbor_nodes;
	vector<Node>::iterator current_it;
	set<vector<Node>::iterator>::iterator current_set_it;
	set<vector<Node>::iterator>::iterator set_it;
	// Determine node connectivity.
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				coords.setXYZ(x, y, z);
				createNode(temp_node, coords);
				Node_vector[lattice.getSiteIndex(x, y, z)] = temp_node;
			}
		}
	}
	// Initialize the path distances of top and bottom surfaces of the lattice.
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			if (lattice.getSiteType(x, y, 0) == Site_types[0]) {
				path_distances[lattice.getSiteIndex(x, y, 0)] = 1;
			}
			if (lattice.getSiteType(x, y, lattice.getHeight() - 1) == Site_types[1]) {
				path_distances[lattice.getSiteIndex(x, y, lattice.getHeight() - 1)] = 1;
			}
		}
	}
	// The pathfinding algorithm is performed for one domain type at a time.
	for (int n = 0; n < 2; n++) {
		// Use Dijkstra's algorithm to fill in the remaining path distance data.
		cout << ID << ": Executing Dijkstra's algorithm to calculate shortest paths through domain type " << (int)Site_types[n] << ".\n";
		if (n == 0) {
			z = 1;
		}
		else {
			z = lattice.getHeight() - 2;
		}
		// Initialize the neighbor node set.
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				if (lattice.getSiteType(x, y, z) == Site_types[n]) {
					d = -1;
					for (int i = 0; i < 26; i++) {
						if (Node_vector[lattice.getSiteIndex(x, y, z)].neighbor_indices[i] < 0) {
							break;
						}
						if (path_distances[Node_vector[lattice.getSiteIndex(x, y, z)].neighbor_indices[i]] > 0) {
							d_temp = path_distances[Node_vector[lattice.getSiteIndex(x, y, z)].neighbor_indices[i]] + sqrt((float)(Node_vector[lattice.getSiteIndex(x, y, z)]).neighbor_distances_sq[i]);
							if (d < 0 || d_temp < d) {
								d = d_temp;
							}
						}
					}
					if (d > 0) {
						Node_vector[lattice.getSiteIndex(x, y, z)].distance_est = d;
						neighbor_nodes.insert(Node_vector.begin() + lattice.getSiteIndex(x, y, z));
					}
				}
			}
		}
		while (!neighbor_nodes.empty()) {
			// The neighbor nodes set is sorted, so the first node has the shortest estimated path distance and is set to the current node.
			current_set_it = neighbor_nodes.begin();
			current_it = *current_set_it;
			//current_index = (long)(current_it - Node_vector.begin());
			current_index = (long int)distance(Node_vector.begin(), current_it);
			// Insert neighbors of the current node into the neighbor node set.
			for (int i = 0; i < 26; i++) {
				if (Node_vector[current_index].neighbor_indices[i] < 0) {
					break;
				}
				neighbor_index = Node_vector[current_index].neighbor_indices[i];
				// Check that the target neighbor node has not already been finalized.
				if (path_distances[neighbor_index] > 0) {
					continue;
				}
				// Calculate the estimated path distance.
				//d = Node_vector[current_index].distance_est + sqrt((int)(Node_vector[current_index]).neighbor_distances_sq[i]);
				switch (Node_vector[current_index].neighbor_distances_sq[i]) {
				case 1:
					d = Node_vector[current_index].distance_est + 1;
					break;
				case 2:
					d = Node_vector[current_index].distance_est + sqrt_two;
					break;
				case 3:
					d = Node_vector[current_index].distance_est + sqrt_three;
					break;
				default:
					d = Node_vector[current_index].distance_est;
					break;
				}
				// Check if node is already in the neighbor node set.
				set_it = neighbor_nodes.find(Node_vector.begin() + neighbor_index);
				// If the node is not already in the list, update the distance estimate for the target neighbor node and insert the node into the set.
				if (set_it == neighbor_nodes.end()) {
					Node_vector[neighbor_index].distance_est = d;
					neighbor_nodes.insert(Node_vector.begin() + neighbor_index);
				}
				// If it already is in the list, replace it only if the new path distance estimate is smaller.
				else if (d < (*set_it)->distance_est) {
					neighbor_nodes.erase(set_it);
					Node_vector[neighbor_index].distance_est = d;
					neighbor_nodes.insert(Node_vector.begin() + neighbor_index);
				}
			}
			// Finalize the path distance of current node and remove the current node from the neighbor node set.
			path_distances[current_index] = Node_vector[current_index].distance_est;
			neighbor_nodes.erase(current_set_it);
		}
	}
	// Clear allocated memory for the neighbor nodes set
	set<vector<Node>::iterator, NodeIteratorCompare>().swap(neighbor_nodes);
	return true;
}

bool Morphology::calculatePathDistances_ReducedMemory(vector<float>& path_distances) {
	float d;
	float d_temp;
	Coords coords;
	const static float sqrt_two = sqrt((float)2.0);
	const static float sqrt_three = sqrt((float)3.0);
	// Create a temporary node to be used throughout the function.
	// Each node contains a vector with indices of all first- ,second-, and third-nearest neighbors (at most 26 neighbors).
	// Another vector stores the squared distance to each of the neighbors.
	// Each node also has an estimated distance from the destination and the site index.
	Node temp_node;
	// Create Node vector to store neighbor nodes
	// This vector is used in Dijsktra's algorithm to keep a list of all nodes that are neighboring sites that already have their path distances determined.
	// Once the path distance for a particular node is finalized, it is removed from the neighbor node vector and
	// all of its neighboring nodes that have not yet been finalized are added to the neighbor node vector.
	vector<Node> Node_vector;
	int Node_vector_count = 0;
	long int node_index = -1;
	long int current_index = -1;
	long int neighbor_index = -1;
	// Create a boolean vector that keeps track of whether or not nodes have already been added to the node vector
	vector<bool> added;
	added.assign(lattice.getNumSites(), false);
	int z;
	// The pathfinding algorithm is performed for one domain type at a time.
	for (int n = 0; n < 2; n++) {
		// Use Dijkstra's algorithm to fill in the remaining path distance data.
		cout << ID << ": Executing Dijkstra's algorithm to calculate shortest paths through domain type " << (int)Site_types[n] << "." << endl;
		// Clear Node vector
		Node_vector.clear();
		Node_vector.assign(lattice.getLength()*lattice.getWidth(), temp_node);
		// Initialize the path distances with known values
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				if (n == 0 && lattice.getSiteType(x, y, 0) == Site_types[n]) {
					path_distances[lattice.getSiteIndex(x, y, 0)] = 1;
					added[lattice.getSiteIndex(x, y, 0)] = true;
				}
				if (n == 1 && lattice.getSiteType(x, y, lattice.getHeight() - 1) == Site_types[n]) {
					path_distances[lattice.getSiteIndex(x, y, lattice.getHeight() - 1)] = 1;
					added[lattice.getSiteIndex(x, y, lattice.getHeight() - 1)] = true;
				}
			}
		}
		// Initialize the neighbor node vector.
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				if (n == 0) {
					z = 1;
				}
				else if (n == 1) {
					z = lattice.getHeight() - 2;
				}
				else {
					continue;
				}
				if (lattice.getSiteType(x, y, z) == Site_types[n]) {
					d = -1;
					coords.setXYZ(x, y, z);
					createNode(temp_node, coords);
					for (int i = 0; i < 26; i++) {
						neighbor_index = temp_node.neighbor_indices[i];
						if (neighbor_index < 0) {
							break;
						}
						if (path_distances[neighbor_index] > 0) {
							d_temp = path_distances[neighbor_index] + sqrt((float)temp_node.neighbor_distances_sq[i]);
							if (d < 0 || d_temp < d) {
								d = d_temp;
							}
						}
					}
					if (d > 0) {
						temp_node.distance_est = d;
						Node_vector[Node_vector_count] = temp_node;
						Node_vector_count++;
						added[temp_node.site_index] = true;
					}
				}
			}
		}
		// The pathfinding algorithm proceeds until there are no nodes left in the neighbor node vector.
		while (!Node_vector_count == 0) {
			// Identify the node with the shortest estimated path distance as the current node.
			d_temp = -1;
			for (int i = 0; i < Node_vector_count; i++) {
				if (d_temp < 0 || Node_vector[i].distance_est < d_temp) {
					current_index = i;
					d_temp = Node_vector[i].distance_est;
				}
			}
			// Insert any unfinalized neighbors of the current node into the neighbor node vector, and check if any already added nodes need to be updated.
			for (int i = 0; i < 26; i++) {
				neighbor_index = Node_vector[current_index].neighbor_indices[i];
				// Check if the target neighbor node is valid.
				if (neighbor_index < 0) {
					break;
				}
				// Check if the target neighbor node has been finalized.
				else if (!(path_distances[neighbor_index] > 0)) {
					// Calculate the estimated path distance to the target neighbor node.
					switch (Node_vector[current_index].neighbor_distances_sq[i]) {
					case 1:
						d = Node_vector[current_index].distance_est + 1;
						break;
					case 2:
						d = Node_vector[current_index].distance_est + sqrt_two;
						break;
					case 3:
						d = Node_vector[current_index].distance_est + sqrt_three;
						break;
					default:
						d = Node_vector[current_index].distance_est;
						break;
					}
					// Check if the target neighbor node has already been added to the Node vector.
					// If not, create the node, update the distance estimate, and add it to the Node vector.
					if (!added[neighbor_index]) {
						coords = lattice.getSiteCoords(neighbor_index);
						createNode(temp_node, coords);
						temp_node.distance_est = d;
						if (Node_vector_count < (int)Node_vector.size()) {
							Node_vector[Node_vector_count] = temp_node;
							Node_vector_count++;
						}
						else {
							Node_vector.push_back(temp_node);
							Node_vector_count++;
						}
						added[neighbor_index] = true;
					}
					// If it has already been added to the node vector, find it, and update the distance estimate only if the new path distance estimate is smaller.
					else {
						// Find the location of the target neighbor node in the node vector
						node_index = -1;
						for (int j = 0; j < Node_vector_count; j++) {
							if (Node_vector[j].site_index == neighbor_index) {
								node_index = j;
								break;
							}
						}
						if (node_index < 0) {
							cout << ID << ": Error! A node designated as added could not be found in the node vector." << endl;
							return false;
						}
						// Update the distance estimate of the neighbor node if a shorter path has been located.
						if (d < Node_vector[node_index].distance_est) {
							Node_vector[node_index].distance_est = d;
						}
					}
				}
			}
			// Finalize the path distance of current node and remove it from the neighbor node vector.
			path_distances[Node_vector[current_index].site_index] = Node_vector[current_index].distance_est;
			Node_vector[current_index] = Node_vector[Node_vector_count - 1];
			Node_vector_count--;
		}
	}
	return true;
}

bool Morphology::calculateTortuosity(const char site_type, const bool enable_reduced_memory) {
	int bin;
	bool success;
	bool electrode_num;
	if (site_type == (char)1) {
		electrode_num = false;
	}
	else if (site_type == (char)2) {
		electrode_num = true;
	}
	else {
		cout << ID << ": Error! Tortuosity can only be calculated for site types 1 and 2." << endl;
		return false;
	}
	// The shortest path for each site is stored in the path_distances vector.
	// The path distances are initialized to -1.
	vector<float> path_distances(lattice.getNumSites(), -1.0);
	int site_type_index = getSiteTypeIndex(site_type);
	// Two different path distance calculation implementations are available.
	// The reduced memory implementation uses less memory but takes more calculation time.
	// It is designed to be used when creating large lattices to prevent running out of system memory.
	if (enable_reduced_memory) {
		success = calculatePathDistances_ReducedMemory(path_distances);
	}
	else {
		success = calculatePathDistances(path_distances);
	}
	if (!success) {
		cout << ID << ": Error calculating path distances!" << endl;
		return false;
	}
	// Construct tortuosity histograms from the path data
	// Tortuosity values are rounded to the nearest 0.02, resulting in bins centered at 1, 1.02, 1.04, etc.
	TortuosityHistogram_data[site_type_index].assign(1, 0);
	int counts = 0;
	for (int i = 0; i < (int)path_distances.size(); i++) {
		if (lattice.getSiteType(i) == site_type && path_distances[i] > 0) {
			if (!electrode_num) {
				bin = round_int(50 * (path_distances[i] / (lattice.getSiteCoords(i).z + 1)) - 49);
			}
			else {
				bin = round_int(50 * (path_distances[i] / (lattice.getHeight() - lattice.getSiteCoords(i).z)) - 49);
			}
			while (bin >= (int)TortuosityHistogram_data[site_type_index].size()) {
				TortuosityHistogram_data[site_type_index].push_back(0.0);
			}
			TortuosityHistogram_data[site_type_index][bin] += 1.0;
			counts++;

		}
	}
	for (int i = 0; i < (int)TortuosityHistogram_data[site_type_index].size(); i++) {
		TortuosityHistogram_data[site_type_index][i] /= (double)counts;
	}
	// In addition to the overall tortuosity histograms from all sites, the end-to-end tortuosity distribution is collected.
	// the end-to-end describes the distribution of tortuosities for the all pathways from the top surface to the bottom surface of the lattice.
	int index;
	Tortuosity_data[site_type_index].assign(lattice.getLength()*lattice.getWidth(), -1);
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			index = lattice.getWidth()*x + y;
			if (!electrode_num) {
				if (lattice.getSiteType(x, y, lattice.getHeight() - 1) == site_type && path_distances[lattice.getSiteIndex(x, y, lattice.getHeight() - 1)] > 0) {
					Tortuosity_data[site_type_index][index] = path_distances[lattice.getSiteIndex(x, y, lattice.getHeight() - 1)] / (float)lattice.getHeight();
				}
			}
			else {
				if (lattice.getSiteType(x, y, 0) == site_type && path_distances[lattice.getSiteIndex(x, y, 0)] > 0) {
					Tortuosity_data[site_type_index][index] = path_distances[lattice.getSiteIndex(x, y, 0)] / (float)lattice.getHeight();
				}
			}
		}
	}
	// Any sites which are not connected their respective surface, will have a zero path distance and are identified as part of island domains.
	// Calculate island volume fraction
	Island_volume.assign((int)Site_types.size(), 0);
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				if (lattice.getSiteType(x, y, z) == Site_types[0] && path_distances[lattice.getSiteIndex(x, y, z)] < 1) {
					Island_volume[0]++;
				}
				if (lattice.getSiteType(x, y, z) == Site_types[1] && path_distances[lattice.getSiteIndex(x, y, z)] < 1) {
					Island_volume[1]++;
				}
			}
		}
	}
	return true;
}

void Morphology::enableThirdNeighborInteraction() {
	Enable_third_neighbor_interaction = true;
}

void Morphology::createCheckerboardMorphology() {
	addSiteType((char)1);
	addSiteType((char)2);
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				if ((x + y + z) % 2 == 0) {
					lattice.setSiteType(x, y, z, (char)1);
					Site_type_counts[0]++;
				}
				else {
					lattice.setSiteType(x, y, z, (char)2);
					Site_type_counts[1]++;
				}
			}
		}
	}
	// This function calculates the actual mix fraction and updates Mix_fraction vector
	calculateMixFractions();
}

void Morphology::createNode(Node& node, const Coords& coords) {
	Coords coords_dest;
	for (int i = 0; i < 26; i++) {
		node.neighbor_indices[i] = -1;
		node.neighbor_distances_sq[i] = 0;
	}
	node.site_index = lattice.getSiteIndex(coords);
	int neighbor_count = 0;
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			for (int k = -1; k <= 1; k++) {
				if (!lattice.checkMoveValidity(coords, i, j, k)) {
					continue;
				}
				if (coords.z + k < 0 || coords.z + k >= lattice.getHeight()) {
					continue;
				}
				lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
				if (lattice.getSiteType(coords) == lattice.getSiteType(coords_dest)) {
					node.neighbor_indices[neighbor_count] = lattice.getSiteIndex(coords_dest);
					node.neighbor_distances_sq[neighbor_count] = (char)(i*i + j * j + k * k);
					neighbor_count++;
				}
			}
		}
	}
}

void Morphology::createRandomMorphology(const vector<double>& mix_fractions) {
	for (int n = 0; n < (int)mix_fractions.size(); n++) {
		addSiteType((char)(n + 1));
	}
	if (mix_fractions.size() != Site_types.size()) {
		cout << ID << ": Error creating random morphology: size of mix_fractions vector must be equal to the size of the Site_types vector." << endl;
		return;
	}
	double sum = 0;
	for (int n = 0; n < (int)Site_types.size(); n++) {
		if (mix_fractions[n] < 0) {
			cout << ID << ": Error creating random morphology: All mix fractions must be greater than or equal to zero." << endl;
			return;
		}
		sum += mix_fractions[n];
	}
	if ((sum - 1.0) > 1e-6) {
		cout << ID << ": Error creating random morphology: Sum of all mix fractions must be equal to one." << endl;
		return;
	}
	vector<double> thresholds = mix_fractions;
	for (int n = 1; n < (int)Site_types.size(); n++) {
		thresholds[n] += thresholds[n - 1];
	}
	Mix_fractions = mix_fractions;
	double randn;
	bool success = false;
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				success = false;
				randn = rand01();
				for (int n = 0; n < (int)Site_types.size(); n++) {
					if (randn < thresholds[n]) {
						lattice.setSiteType(x, y, z, Site_types[n]);
						Site_type_counts[n]++;
						success = true;
						break;
					}
				}
				if (!success) {
					lattice.setSiteType(x, y, z, Site_types.back());
					Site_type_counts.back()++;
				}
			}
		}
	}
	// This function calculates the actual mix fraction and updates Mix_fraction
	calculateMixFractions();
}

void Morphology::executeIsingSwapping(const int num_MCsteps, const double interaction_energy1, const double interaction_energy2, const bool enable_growth_pref, const int growth_direction, const double additional_interaction) {
	int loop_count = 0;
	// N counts the number of MC steps that have been executed
	int N = 0;
	int neighbor_count;
	long int main_site_index;
	long int neighbor_site_index;
	char temp;
	double energy_delta, probability;
	Coords main_site_coords;
	vector<long int> neighbors;
	const long int *neighbor_pointer;
	neighbors.assign(6, 0);
	initializeNeighborInfo();
	// Begin site swapping
	int m = 1;
	while (N < num_MCsteps) {
		// Randomly choose a site in the lattice
		main_site_coords = lattice.generateRandomCoords();
		main_site_index = lattice.getSiteIndex(main_site_coords);
		// If site is not an interfacial site, start over again
		// If total number of first-nearest neighbors = number of first-nearest neighbors of the same type, then the site is not at an interface
		if (Neighbor_info[main_site_index].total1 == Neighbor_counts[main_site_index].sum1) {
			continue;
		}
		// Randomly choose a nearest neighbor site that has a different type
		// First find all nearest neighbor sites of differing type
		neighbor_count = 0;
		neighbor_pointer = Neighbor_info[main_site_index].first_indices;
		for (int i = 0; i < 6; i++) {
			if (*(neighbor_pointer + i) >= 0 && lattice.getSiteType(main_site_index) != lattice.getSiteType(*(neighbor_pointer + i))) {
				// Store site index of differing neighbor site
				neighbors[neighbor_count] = *(neighbor_pointer + i);
				neighbor_count++;
			}
		}
		// Randomly select one of the differing neighbor sites
		neighbor_site_index = neighbors[(int)floor(rand01()*neighbor_count)];
		// Calculate energy change and swapping probability
		energy_delta = calculateEnergyChangeSimple(main_site_index, neighbor_site_index, interaction_energy1, interaction_energy2);
		if (enable_growth_pref) {
			energy_delta += calculateAdditionalEnergyChange(main_site_index, neighbor_site_index, growth_direction, additional_interaction);
		}
		probability = exp(-energy_delta) / (1 + exp(-energy_delta));
		if (rand01() <= probability) {
			// Swap Sites
			temp = lattice.getSiteType(main_site_index);
			lattice.setSiteType(main_site_index, lattice.getSiteType(neighbor_site_index));
			lattice.setSiteType(neighbor_site_index, temp);
			// Update neighbor counts
			updateNeighborCounts(main_site_index, neighbor_site_index);
		}
		loop_count++;
		// One MC step has been completed when loop_count is equal to the number of sites in the lattice
		if (loop_count == lattice.getNumSites()) {
			N++;
			loop_count = 0;
		}
		if (N == 100 * m) {
			cout << ID << ": " << N << " MC steps completed." << endl;
			m++;
		}
	}
	vector<NeighborCounts>().swap(Neighbor_counts);
	vector<NeighborInfo>().swap(Neighbor_info);
}

void Morphology::executeMixing(const double interfacial_width, const double interfacial_conc) {
	vector<int> sites_maj;
	vector<int> sites_min;
	int site_maj;
	int site_min;
	int target;
	int site_count = 0;
	char majority_type;
	char minority_type;
	double minority_conc;
	Coords coords;
	// Based on the interfacial concentration, the majority and minority types are defined
	if (interfacial_conc <= 0.5) {
		majority_type = (char)2;
		minority_type = (char)1;
		minority_conc = interfacial_conc;
	}
	else {
		majority_type = (char)1;
		minority_type = (char)2;
		minority_conc = 1 - interfacial_conc;
	}
	// The minority type sites are mixed into the majority type sites at the interface.
	// First the majority site reservoir is determined by finding all majority type sites within (1-minority_conc)*width distance from the interface and adding them to a list
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				coords.setXYZ(x, y, z);
				if (lattice.getSiteType(x, y, z) == majority_type) {
					if (isNearInterface(coords, (1 - minority_conc)*interfacial_width)) {
						sites_maj.push_back(lattice.getSiteIndex(x, y, z));
						site_count++;
					}
				}
			}
		}
	}
	// Then the minority site reservoir that will be used to be mix into the majority reservoir to yield the desired final interfacial concentration is determined
	// More minority sites than the number of swaps must be included in the minority reservoir so that the desired interfacial concentration is also reached on the minority side of the interface.
	target = (int)(site_count*minority_conc / (1 - minority_conc));
	int N = 0;
	double range = 1;
	// The minority type sites adjacent to the interface are identified, counted, and added to a list.
	// If this does not yield enough minority sites for swapping, then the range is slowly incremented and minority sites farther from the interface are included until there are enough sites.
	while (N < target) {
		sites_min.clear();
		N = 0;
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					coords.setXYZ(x, y, z);
					if (lattice.getSiteType(x, y, z) == minority_type) {
						if (isNearInterface(coords, range)) {
							sites_min.push_back(lattice.getSiteIndex(x, y, z));
							N++;
						}
					}
				}
			}
		}
		range += 0.1;
	}
	// Sites are randomly chosen from the minority and majority reservoirs and swapped until the desired interfacial concentration is reached.
	for (int i = 0; i < site_count*minority_conc; i++) {
		uniform_int_distribution<int> dist_maj(0, (int)sites_maj.size() - 1);
		auto rand_maj = bind(dist_maj, ref(gen));
		uniform_int_distribution<int> dist_min(0, (int)sites_min.size() - 1);
		auto rand_min = bind(dist_min, ref(gen));
		site_maj = rand_maj();
		site_min = rand_min();
		lattice.setSiteType(sites_maj[site_maj], minority_type);
		lattice.setSiteType(sites_min[site_min], majority_type);
		// Both site types are removed from the lists once they are swapped to prevent unswapping.
		sites_maj.erase(sites_maj.begin() + site_maj);
		sites_min.erase(sites_min.begin() + site_min);
	}
}

void Morphology::executeSmoothing(const double smoothing_threshold, const int rescale_factor) {
	double roughness_factor;
	Coords coords, coords_dest;
	static int radius = (rescale_factor<=2) ? (int)ceil(3.0/2.0) : (int)ceil((double)(rescale_factor + 1) / 2);
	static int cutoff_squared = (rescale_factor<=2) ? (int)floor((3.0/2.0)*(3.0 / 2.0)) : (int)floor(((double)(rescale_factor + 1) / 2)*((double)(rescale_factor + 1) / 2));
	// The boolean vector consider_smoothing keeps track of whether each site is near the interface and should be considered for smoothing.
	// Sites in the interior of the domains or at very smooth interfaces do not need to be continually reconsidered for smoothing.
	vector<bool> consider_smoothing;
	consider_smoothing.assign(lattice.getNumSites(), true);
	int site_count = 1;
	while (site_count > 0) {
		site_count = 0;
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					if (!consider_smoothing[lattice.getSiteIndex(x, y, z)]) {
						continue;
					}
					coords.setXYZ(x, y, z);
					// Calculate the roughness factor of the target site.
					roughness_factor = calculateDissimilarFraction(coords, rescale_factor);
					// Swap the site's type if the roughness_factor is greater than the smoothing_threshold.
					if (roughness_factor > smoothing_threshold) {
						if (lattice.getSiteType(x, y, z) == (char)1) {
							lattice.setSiteType(x, y, z, (char)2);
						}
						else if (lattice.getSiteType(x, y, z) == (char)2) {
							lattice.setSiteType(x, y, z, (char)1);
						}
						site_count++;
						// When a site swaps types, all surrounding sites must be reconsidered for smoothing.
						for (int i = -radius; i <= radius; i++) {
							for (int j = -radius; j <= radius; j++) {
								for (int k = -radius; k <= radius; k++) {
									if (i*i + j * j + k * k > cutoff_squared) {
										continue;
									}
									if (!lattice.checkMoveValidity(coords, i, j, k)) {
										continue;
									}
									lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
									consider_smoothing[lattice.getSiteIndex(coords_dest)] = true;
								}
							}
						}
					}
					// Sites with a low roughness_factor are not swapped and removed from reconsideration.
					else {
						consider_smoothing[lattice.getSiteIndex(x, y, z)] = false;
					}
				}
			}
		}
	}
	// The smoothing process can change the mix fraction, so the final mix fraction is recalculated and the Mix_fraction property is updated.
	calculateMixFractions();
}

vector<double> Morphology::getCorrelationData(const char site_type) const {
	if (Correlation_data[getSiteTypeIndex(site_type)][0] < 0) {
		cout << ID << ": Error getting correlation data: Correlation data has not been calculated." << endl;
	}
	return Correlation_data[getSiteTypeIndex(site_type)];
}

vector<double> Morphology::getDepthCompositionData(const char site_type) const {
	return Depth_composition_data[getSiteTypeIndex(site_type)];
}


vector<double> Morphology::getDepthDomainSizeData(const char site_type) const {
	return Depth_domain_size_data[getSiteTypeIndex(site_type)];
}

vector<double> Morphology::getDepthIVData(const char site_type) const {
	return Depth_iv_data[getSiteTypeIndex(site_type)];
}

double Morphology::getDomainAnisotropy(const char site_type) const {
	return Domain_anisotropies[getSiteTypeIndex(site_type)];
}

double Morphology::getDomainSize(char site_type) const {
	return Domain_sizes[getSiteTypeIndex(site_type)];
}

int Morphology::getHeight() const {
	return lattice.getHeight();
}

vector<double> Morphology::getInterfacialHistogram(char site_type) const {
	return InterfacialHistogram_data[getSiteTypeIndex(site_type)];
}

double Morphology::getIslandVolumeFraction(char site_type) const {
	return (double)Island_volume[getSiteTypeIndex(site_type)] / (double)lattice.getNumSites();
}

int Morphology::getLength() const {
	return lattice.getLength();
}

double Morphology::getMixFraction(const char site_type) const {
	return Mix_fractions[getSiteTypeIndex(site_type)];
}

void Morphology::getSiteSampling(vector<long int>& site_indices, const char site_type, const int N_sites_max) {
	vector<long int> all_sites(Site_type_counts[getSiteTypeIndex(site_type)], 0);
	int m = 0;
	for (long n = 0; n < lattice.getNumSites(); n++) {
		if (lattice.getSiteType(n) == site_type && m < (int)all_sites.size()) {
			all_sites[m] = n;
			m++;
		}
	}
	shuffle(all_sites.begin(), all_sites.end(), gen);
	if (N_sites_max > (int)all_sites.size()) {
		site_indices.assign(all_sites.begin(), all_sites.end());
	}
	else {
		site_indices.assign(all_sites.begin(), all_sites.begin() + N_sites_max);
	}
}

void Morphology::getSiteSamplingZ(vector<long int>& site_indices, const char site_type, const int N_sites_max, const int z) {
	vector<long int> all_sites;
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			if (lattice.getSiteType(x, y, z) == site_type) {
				all_sites.push_back(lattice.getSiteIndex(x, y, z));
			}
		}
	}
	shuffle(all_sites.begin(), all_sites.end(), gen);
	if (N_sites_max > (int)all_sites.size()) {
		site_indices.assign(all_sites.begin(), all_sites.end());
	}
	else {
		site_indices.assign(all_sites.begin(), all_sites.begin() + N_sites_max);
	}
}

int Morphology::getSiteTypeIndex(const char site_type) const {
	for (int i = 0; i < (int)Site_types.size(); i++) {
		if (site_type == Site_types[i]) {
			return i;
		}
	}
	cout << ID << ": Error! Site type " << (int)site_type << " was not found in the Site_types vector." << endl;
	return -1;
}

vector<float> Morphology::getTortuosityData(char site_type) const {
	vector<float> output_data;
	for (int i = 0; i < (int)Tortuosity_data[getSiteTypeIndex(site_type)].size(); i++) {
		if (Tortuosity_data[getSiteTypeIndex(site_type)][i] > 0) {
			output_data.push_back(Tortuosity_data[getSiteTypeIndex(site_type)][i]);
		}
	}
	return output_data;
}

vector<double> Morphology::getTortuosityHistogram(char site_type) const {
	return TortuosityHistogram_data[getSiteTypeIndex(site_type)];
}

double Morphology::getUnitSize() const {
	return lattice.getUnitSize();
}

int Morphology::getWidth() const {
	return lattice.getWidth();
}

vector<Morphology> Morphology::importTomogramMorphologyFile(const string& info_filename, const string& data_filename, const TomogramImport_Params& import_params) {
	vector<Morphology> morphologies;
	// Input file checking
	if ((int)info_filename.length() > 4 && info_filename.substr(info_filename.length() - 4, 4).compare(".xml") != 0) {
		cout << ID << ": Error! Input tomogram info file " << info_filename << " is not an xml file." << endl;
		return morphologies;
	}
	if ((int)data_filename.length() > 4 && data_filename.substr(data_filename.length() - 4, 4).compare(".raw") != 0) {
		cout << ID << ": Error! Input tomogram data file " << data_filename << "  is not a .raw file." << endl;
		return morphologies;
	}
	// Parameter checking
	switch (import_params.N_extracted_segments) {
	case 1:
		break;
	case 4:
		break;
	case 9:
		break;
	case 16:
		break;
	case 25:
		break;
	case 36:
		break;
	case 49:
		break;
	case 64:
		break;
	case 81:
		break;
	case 100:
		break;
	case 121:
		break;
	case 144:
		break;
	case 169:
		break;
	case 196:
		break;
	default:
		cout << "Error! Invalid value for N_extracted_segments parameter. Value must be 1, 4, 9, 16, 25, 36, 49, 64, 89, or 100 but " << import_params.N_extracted_segments << " was entered." << endl;
		return morphologies;
	}
	// Parse xml info file
	XMLDocument xml_doc;
	Lattice_Params lattice_params;
	string data_format;
	xml_doc.LoadFile(info_filename.c_str());
	if (xml_doc.Error()) {
		xml_doc.PrintError();
		return morphologies;
	}
	// Analyze xml info file based on metadata format version
	string schema_version = xml_doc.FirstChildElement("tomogram_metadata")->Attribute("schema_version");
	if (schema_version.compare("1.0") == 0) {
		// Initialize lattice params based on xml data
		lattice_params.Length = atoi(xml_doc.FirstChildElement("tomogram_metadata")->FirstChildElement("data_info")->FirstChildElement("length")->GetText());
		lattice_params.Width = atoi(xml_doc.FirstChildElement("tomogram_metadata")->FirstChildElement("data_info")->FirstChildElement("width")->GetText());
		lattice_params.Height = atoi(xml_doc.FirstChildElement("tomogram_metadata")->FirstChildElement("data_info")->FirstChildElement("height")->GetText());
		lattice_params.Unit_size = atof(xml_doc.FirstChildElement("tomogram_metadata")->FirstChildElement("data_info")->FirstChildElement("pixel_size")->FirstChildElement("value")->GetText());
		// Save additional xml data needed
		data_format = xml_doc.FirstChildElement("tomogram_metadata")->FirstChildElement("data_info")->FirstChildElement("data_format")->GetText();
		// Initialize Morphology site types
		int site_type_index = 0;
		XMLElement* element_ptr = xml_doc.FirstChildElement("tomogram_metadata")->FirstChildElement("sample_info")->FirstChildElement("composition_info")->FirstChildElement("chemical_component");
		while (element_ptr != 0) {
			addSiteType((char)(site_type_index + 1));
			Mix_fractions[site_type_index] = atof(element_ptr->FirstChildElement("vol_frac")->GetText());
			element_ptr = element_ptr->NextSiblingElement("chemical_component");
			site_type_index++;
		}
	}
	else {
		cout << ID << ": Error! XML info file format version not recognized." << endl;
		return morphologies;
	}
	lattice_params.Enable_periodic_x = false;
	lattice_params.Enable_periodic_y = false;
	lattice_params.Enable_periodic_z = false;
	Lattice lattice_i;
	lattice_i.init(lattice_params, &gen);
	// Open and read data file
	ifstream data_file(data_filename, ifstream::in | ifstream::binary);
	if (data_file.is_open()) {
		data_file.seekg(0, data_file.end);
		streampos N_bytes;
		N_bytes = data_file.tellg();
		data_file.seekg(0, data_file.beg);
		vector<float> data_vec;
		if (data_format.compare("8 bit") == 0) {
			unsigned char* data_block;
			data_block = (unsigned char*)malloc(sizeof(char)*N_bytes);
			data_file.read((char*)data_block, N_bytes);
			data_vec.assign(N_bytes, 0);
			for (int i = 0; i < (int)N_bytes; i++) {
				data_vec[i] = (float)data_block[i];
			}
		}
		else if (data_format.compare("16 bit") == 0) {
			char16_t* data_block;
			data_block = (char16_t*)malloc(sizeof(char)*N_bytes);
			data_file.read((char*)data_block, N_bytes);
			data_vec.assign(N_bytes / 2, 0);
			for (int i = 0; i < (int)(N_bytes / 2); i++) {
				data_vec[i] = (float)data_block[i];
			}
		}
		else if (data_format.compare("32 bit") == 0) {
			int32_t* data_block;
			data_block = (int32_t*)malloc(sizeof(char)*N_bytes);
			data_file.read((char*)data_block, N_bytes);
			data_vec.assign(N_bytes / 4, 0);
			for (int i = 0; i < (int)(N_bytes / 4); i++) {
				data_vec[i] = (float)data_block[i];
			}
		}
		data_file.close();
		if (data_vec.size() != lattice_i.getNumSites()) {
			cout << "Error! The imported tomogram binary file does not contain the correct number of sites." << endl;
			cout << "The initial lattice has " << lattice_i.getNumSites() << " sites but the data file has " << data_vec.size() << " entries." << endl;
			return morphologies;
		}
		// Determine final lattice dimensions based on desired unit size
		lattice_params.Length = (int)floor(lattice_params.Length*(lattice_params.Unit_size / import_params.Desired_unit_size));
		lattice_params.Width = (int)floor(lattice_params.Width*(lattice_params.Unit_size / import_params.Desired_unit_size));
		lattice_params.Height = (int)floor(lattice_params.Height*(lattice_params.Unit_size / import_params.Desired_unit_size));
		lattice_params.Unit_size = import_params.Desired_unit_size;
		lattice.init(lattice_params, &gen);
		int index = 0;
		vector<float> data_vec_new(data_vec.size());
		for (int k = 0; k < lattice_i.getHeight(); k++) {
			for (int j = 0; j < lattice_i.getWidth(); j++) {
				for (int i = 0; i < lattice_i.getLength(); i++) {
					data_vec_new[lattice_i.getSiteIndex(i, j, k)] = data_vec[index];
					index++;
				}
			}
		}
		vector<float>().swap(data_vec);
		vector<float> data_vec_final(lattice.getNumSites());
		for (int i = 0; i < lattice.getLength(); i++) {
			for (int j = 0; j < lattice.getWidth(); j++) {
				for (int k = 0; k < lattice.getHeight(); k++) {
					double x = i * lattice.getUnitSize();
					double y = j * lattice.getUnitSize();
					double z = k * lattice.getUnitSize();
					static vector<float> weights(8, 0.0);
					static vector<float> vals(8, 0.0);
					int x1 = (int)floor(x / lattice_i.getUnitSize());
					int x2 = (int)ceil(x / lattice_i.getUnitSize());
					int y1 = (int)floor(y / lattice_i.getUnitSize());
					int y2 = (int)ceil(y / lattice_i.getUnitSize());
					int z1 = (int)floor(z / lattice_i.getUnitSize());
					int z2 = (int)ceil(z / lattice_i.getUnitSize());
					if ((x - x1 * lattice_i.getUnitSize()) < 1e-6 && (y - y1 * lattice_i.getUnitSize()) < 1e-6 && (z - z1 * lattice_i.getUnitSize()) < 1e-6) {
						data_vec_final[lattice.getSiteIndex(i, j, k)] = data_vec_new[lattice_i.getSiteIndex(x1, y1, z1)];
						continue;
					}
					if (x2 >= lattice_i.getLength()) {
						x2 = lattice_i.getLength() - 1;
					}
					if (y2 >= lattice_i.getWidth()) {
						y2 = lattice_i.getWidth() - 1;
					}
					if (z2 >= lattice_i.getHeight()) {
						z2 = lattice_i.getHeight() - 1;
					}
					vals[0] = data_vec_new[lattice_i.getSiteIndex(x1, y1, z1)];
					vals[1] = data_vec_new[lattice_i.getSiteIndex(x1, y2, z1)];
					vals[2] = data_vec_new[lattice_i.getSiteIndex(x2, y1, z1)];
					vals[3] = data_vec_new[lattice_i.getSiteIndex(x2, y2, z1)];
					vals[4] = data_vec_new[lattice_i.getSiteIndex(x1, y1, z2)];
					vals[5] = data_vec_new[lattice_i.getSiteIndex(x1, y2, z2)];
					vals[6] = data_vec_new[lattice_i.getSiteIndex(x2, y1, z2)];
					vals[7] = data_vec_new[lattice_i.getSiteIndex(x2, y2, z2)];
					weights[0] = (float)(1.0 / (intpow(x - x1 * lattice_i.getUnitSize(), 2) + intpow(y - y1 * lattice_i.getUnitSize(), 2) + intpow(z - z1 * lattice_i.getUnitSize(), 2)));
					weights[1] = (float)(1.0 / (intpow(x - x1 * lattice_i.getUnitSize(), 2) + intpow(y - y2 * lattice_i.getUnitSize(), 2) + intpow(z - z1 * lattice_i.getUnitSize(), 2)));
					weights[2] = (float)(1.0 / (intpow(x - x2 * lattice_i.getUnitSize(), 2) + intpow(y - y1 * lattice_i.getUnitSize(), 2) + intpow(z - z1 * lattice_i.getUnitSize(), 2)));
					weights[3] = (float)(1.0 / (intpow(x - x2 * lattice_i.getUnitSize(), 2) + intpow(y - y2 * lattice_i.getUnitSize(), 2) + intpow(z - z1 * lattice_i.getUnitSize(), 2)));
					weights[4] = (float)(1.0 / (intpow(x - x1 * lattice_i.getUnitSize(), 2) + intpow(y - y1 * lattice_i.getUnitSize(), 2) + intpow(z - z2 * lattice_i.getUnitSize(), 2)));
					weights[5] = (float)(1.0 / (intpow(x - x1 * lattice_i.getUnitSize(), 2) + intpow(y - y2 * lattice_i.getUnitSize(), 2) + intpow(z - z2 * lattice_i.getUnitSize(), 2)));
					weights[6] = (float)(1.0 / (intpow(x - x2 * lattice_i.getUnitSize(), 2) + intpow(y - y1 * lattice_i.getUnitSize(), 2) + intpow(z - z2 * lattice_i.getUnitSize(), 2)));
					weights[7] = (float)(1.0 / (intpow(x - x2 * lattice_i.getUnitSize(), 2) + intpow(y - y2 * lattice_i.getUnitSize(), 2) + intpow(z - z2 * lattice_i.getUnitSize(), 2)));
					std::transform(vals.begin(), vals.end(), weights.begin(), vals.begin(), multiplies<float>());
					data_vec_final[lattice.getSiteIndex(i, j, k)] = accumulate(vals.begin(), vals.end(), 0.0f) / accumulate(weights.begin(), weights.end(), 0.0f);
				}
			}
		}

		vector<float>().swap(data_vec_new);
		// Normalize test data vector to have range from 0 to 100
		//vector<double> test_vec(data_vec_final.size());
		//for (int i = 0; i < data_vec_final.size(); i++) {
		//	test_vec[i] = (double)data_vec_final[i];
		//}
		//double maximum = *max_element(test_vec.begin(), test_vec.end());
		//double minimum = *min_element(test_vec.begin(), test_vec.end());
		//std::transform(test_vec.begin(), test_vec.end(), test_vec.begin(), bind(minus<double>(), std::placeholders::_1, minimum));
		//std::transform(test_vec.begin(), test_vec.end(), test_vec.begin(), bind(divides<double>(), std::placeholders::_1, 0.01*maximum));
		//auto prob = calculateProbabilityHist(test_vec, 1.0);
		//Utils::outputVectorToFile(prob, "hist_data1.txt");
		//vector<double>().swap(test_vec);
		// 
		// Interpret data vector using probability analysis method
		if (import_params.Enable_probability_analysis) {
			float avg = (float)vector_avg(data_vec_final);
			float min_val = *min_element(data_vec_final.begin(), data_vec_final.end());
			float max_val = *max_element(data_vec_final.begin(), data_vec_final.end());
			std::transform(data_vec_final.begin(), data_vec_final.end(), data_vec_final.begin(), bind(minus<float>(), std::placeholders::_1, avg));
			min_val -= avg;
			max_val -= avg;
			for (int i = 0; i < (int)data_vec_final.size(); i++) {
				if (data_vec_final[i] > 0) {
					data_vec_final[i] = pow(abs(data_vec_final[i]), (float)import_params.Probability_scaling_exponent);
				}
				else {
					data_vec_final[i] = -pow(abs(data_vec_final[i]), (float)import_params.Probability_scaling_exponent);
				}
			}
			min_val = -pow(abs(min_val), (float)import_params.Probability_scaling_exponent);
			max_val = pow(max_val, (float)import_params.Probability_scaling_exponent);
			for (int i = 0; i < (int)data_vec_final.size(); i++) {
				data_vec_final[i] = (data_vec_final[i] - min_val) / (max_val - min_val);
			}
			for (int x = 0; x < lattice.getLength(); x++) {
				for (int y = 0; y < lattice.getWidth(); y++) {
					for (int z = 0; z < lattice.getHeight(); z++) {
						if (rand01() < data_vec_final[lattice.getSiteIndex(x, y, z)]) {
							lattice.setSiteType(lattice.getSiteIndex(x, y, z), (char)1);
							Site_type_counts[getSiteTypeIndex((char)1)]++;
						}
						else {
							lattice.setSiteType(lattice.getSiteIndex(x, y, z), (char)2);
							Site_type_counts[getSiteTypeIndex((char)2)]++;
						}
					}
				}
			}
		}
		if (import_params.Enable_cutoff_analysis) {
			vector<float> temp_vec = data_vec_final;
			sort(temp_vec.begin(), temp_vec.end());
			float cutoff = temp_vec[round_int(Mix_fractions[1] * (temp_vec.size() - 1))];
			vector<float>().swap(temp_vec);
			int mix_bin = import_params.Mixed_greyscale_width;
			//addSiteType((char)3);
			for (int x = 0; x < lattice.getLength(); x++) {
				for (int y = 0; y < lattice.getWidth(); y++) {
					for (int z = 0; z < lattice.getHeight(); z++) {
						if (data_vec_final[lattice.getSiteIndex(x, y, z)] > cutoff + (mix_bin / 2)) {
							lattice.setSiteType(lattice.getSiteIndex(x, y, z), (char)1);
							Site_type_counts[getSiteTypeIndex((char)1)]++;
						}
						else if (data_vec_final[lattice.getSiteIndex(x, y, z)] < cutoff - (mix_bin / 2)) {
							lattice.setSiteType(lattice.getSiteIndex(x, y, z), (char)2);
							Site_type_counts[getSiteTypeIndex((char)2)]++;
						}
						else {
							//lattice.setSiteType(lattice.getSiteIndex(x, y, z), (char)3);
							//Site_type_counts[getSiteTypeIndex((char)3)]++;
							if (rand01() > 0.5) {
								lattice.setSiteType(lattice.getSiteIndex(x, y, z), (char)1);
								Site_type_counts[getSiteTypeIndex((char)1)]++;
							}
							else {
								lattice.setSiteType(lattice.getSiteIndex(x, y, z), (char)2);
								Site_type_counts[getSiteTypeIndex((char)2)]++;
							}
						}
					}
				}
			}
		}
		// Generate morphology set
		cout << ID << ": Creating morphology set from tomogram data." << endl;
		if (import_params.N_extracted_segments > 1) {
			int dim = round_int(sqrt(import_params.N_extracted_segments));
			int new_length = lattice.getLength() / dim;
			// make dimension even
			if (new_length % 2 != 0) {
				new_length--;
			}
			int new_width = lattice.getWidth() / dim;
			if (new_width % 2 != 0) {
				new_width--;
			}
			int size = (new_length < new_width) ? new_length : new_width;
			int offset_x = (lattice.getLength() - size * dim) / 2;
			int offset_y = (lattice.getWidth() - size * dim) / 2;
			Lattice sublattice;
			for (int x = 0; x < dim; x++) {
				for (int y = 0; y < dim; y++) {
					sublattice = lattice.extractSublattice(x*size + offset_x, size, y*size + offset_y, size, 0, lattice.getHeight());
					morphologies.push_back(Morphology(sublattice, x * dim + y));
					morphologies[x * dim + y].calculateMixFractions();
				}
			}
		}
		else {
			morphologies.push_back(*this);
		}
	}
	else {
		cout << "Error! Tomogram binary file " << data_filename << " could not be opened." << endl;
	}
	return morphologies;
}

bool Morphology::importMorphologyFile(ifstream& infile) {
	string var;
	int x = 0;
	int y = 0;
	int z = 0;
	char type = (char)0;
	int j;
	string line;
	bool isV4 = false;
	Lattice_Params params;
	params.Unit_size = 1.0;
	// Check if file is in compressed format or not
	getline(infile, line);
	bool is_file_compressed = !(line.find("uncompressed") != string::npos);
	if (line.substr(0, 9).compare("Ising_OPV") == 0) {
		// check if morphology is from Ising_OPV v4
		if (line.substr(10, 2).compare("v4") == 0) {
			isV4 = true;
		}
		// skip first line
		// get next line
		getline(infile, line);
	}
	params.Length = atoi(line.c_str());
	getline(infile, line);
	params.Width = atoi(line.c_str());
	getline(infile, line);
	params.Height = atoi(line.c_str());
	if (isV4) {
		getline(infile, line);
		params.Enable_periodic_x = (bool)atoi(line.c_str());
		getline(infile, line);
		params.Enable_periodic_y = (bool)atoi(line.c_str());
		getline(infile, line);
		params.Enable_periodic_z = (bool)atoi(line.c_str());
	}
	else {
		params.Enable_periodic_x = true;
		params.Enable_periodic_y = true;
		// Assume z-direction periodic boundary is disabled for an imported morphology
		params.Enable_periodic_z = false;
	}
	lattice.init(params, &gen);
	if (isV4) {
		getline(infile, line);
		int num_types = atoi(line.c_str());
		for (int n = 0; n < num_types; n++) {
			addSiteType((char)(n + 1));
		}
		for (int n = 0; n < num_types; n++) {
			getline(infile, line);
			Domain_sizes[n] = atof(line.c_str());
		}
		for (int n = 0; n < num_types; n++) {
			getline(infile, line);
			Mix_fractions[n] = atof(line.c_str());
		}
	}
	else {
		addSiteType((char)1);
		addSiteType((char)2);
		getline(infile, line);
		Domain_sizes[0] = atof(line.c_str());
		getline(infile, line);
		Domain_sizes[1] = atof(line.c_str());
		getline(infile, line);
		Mix_fractions[0] = atof(line.c_str());
		Mix_fractions[1] = 1 - Mix_fractions[0];
	}
	if (!is_file_compressed) {
		while ((infile).good()) {
			getline(infile, line);
			stringstream linestream(line);
			j = 1;
			while (linestream.good()) {
				getline(linestream, var, ',');
				switch (j) {
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
			lattice.setSiteType(x, y, z, type);
		}
	}
	else {
		int site_count = 0;
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					if (site_count == 0) {
						getline(infile, line);
						stringstream linestream(line);
						if (!linestream.good()) {
							cout << "Error parsing file.  End of file reached before expected." << endl;
							return false;
						}
						type = (char)atoi((line.substr(0, 1)).c_str());
						site_count = atoi((line.substr(1, line.length() - 1)).c_str());
					}
					lattice.setSiteType(x, y, z, type);
					Site_type_counts[getSiteTypeIndex(type)]++;
					site_count--;
				}
			}
		}
	}
	calculateMixFractions();
	return true;
}

void Morphology::initializeNeighborInfo() {
	Coords coords, coords_dest;
	char sum1, sum2, sum3;
	char total1, total2, total3;
	int first_neighbor_count, second_neighbor_count, third_neighbor_count;
	char site_type;
	// Initialize neighbor counts (this data is used in the calculateEnergyChangeSimple function)
	NeighborCounts counts;
	Neighbor_counts.assign(lattice.getNumSites(), counts);
	NeighborInfo info;
	Neighbor_info.assign(lattice.getNumSites(), info);
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				coords.setXYZ(x, y, z);
				sum1 = 0;
				sum2 = 0;
				sum3 = 0;
				// There are in total 6 first-nearest, 12 second-nearest, and 8 third-nearest neighbors
				total1 = 6;
				total2 = 12;
				total3 = 8;
				// Keep track of which neighbors have been calculated
				first_neighbor_count = 0;
				second_neighbor_count = 0;
				third_neighbor_count = 0;
				// Calculate similar neighbors around x,y,z
				site_type = lattice.getSiteType(x, y, z);
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						for (int k = -1; k <= 1; k++) {
							if (i == 0 && j == 0 && k == 0) {
								continue;
							}
							if (!lattice.isZPeriodic()) {
								if (z + k >= lattice.getHeight() || z + k < 0) { // Check for z boundary
									// Total site counts must be reduced if next to a hard boundary
									switch (i*i + j * j + k * k) {
									case 1:
										total1--;
										Neighbor_info[lattice.getSiteIndex(x, y, z)].first_indices[first_neighbor_count] = -1;
										first_neighbor_count++;
										break;
									case 2:
										total2--;
										Neighbor_info[lattice.getSiteIndex(x, y, z)].second_indices[second_neighbor_count] = -1;
										second_neighbor_count++;
										break;
									case 3:
										total3--;
										Neighbor_info[lattice.getSiteIndex(x, y, z)].third_indices[third_neighbor_count] = -1;
										third_neighbor_count++;
										break;
									default:
										break;
									}
									continue;
								}
							}
							lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
							// Count the number of similar neighbors
							if (lattice.getSiteType(coords_dest) == site_type) {
								switch (i*i + j * j + k * k) {
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
							// Determine neighbor site indices
							switch (i*i + j * j + k * k) {
							case 1:
								Neighbor_info[lattice.getSiteIndex(x, y, z)].first_indices[first_neighbor_count] = lattice.getSiteIndex(coords_dest);
								first_neighbor_count++;
								break;
							case 2:
								Neighbor_info[lattice.getSiteIndex(x, y, z)].second_indices[second_neighbor_count] = lattice.getSiteIndex(coords_dest);
								second_neighbor_count++;
								break;
							case 3:
								Neighbor_info[lattice.getSiteIndex(x, y, z)].third_indices[third_neighbor_count] = lattice.getSiteIndex(coords_dest);
								third_neighbor_count++;
								break;
							default:
								break;
							}
						}
					}
				}
				Neighbor_counts[lattice.getSiteIndex(x, y, z)].sum1 = sum1;
				Neighbor_counts[lattice.getSiteIndex(x, y, z)].sum2 = sum2;
				Neighbor_counts[lattice.getSiteIndex(x, y, z)].sum3 = sum3;
				Neighbor_info[lattice.getSiteIndex(x, y, z)].total1 = total1;
				Neighbor_info[lattice.getSiteIndex(x, y, z)].total2 = total2;
				Neighbor_info[lattice.getSiteIndex(x, y, z)].total3 = total3;
				if (!(Neighbor_counts[lattice.getSiteIndex(x, y, z)] == calculateNeighborCounts(coords))) {
					cout << "Error initializing neighbor counts!" << endl;
				}
			}
		}
	}
}

bool Morphology::isNearInterface(const Coords& coords, const double distance) const {
	int d = (int)floor(distance);
	Coords coords_dest;
	for (int i = -d; i <= d; i++) {
		for (int j = -d; j <= d; j++) {
			for (int k = -d; k <= d; k++) {
				if (d == 1 && abs(i) + abs(j) + abs(k) > 1) {
					continue;
				}
				else if ((double)(i*i + j * j + k * k) > distance*distance) {
					continue;
				}
				if (!lattice.checkMoveValidity(coords, i, j, k)) {
					continue;
				}
				lattice.calculateDestinationCoords(coords, i, j, k, coords_dest);
				if (lattice.getSiteType(coords) != lattice.getSiteType(coords_dest)) {
					return true;
				}
			}
		}
	}
	return false;
}

void Morphology::outputCompositionMaps(ofstream& outfile) const {
	vector<int> counts(Site_types.size(), 0);
	//outfile << "X-Position,Y-Position";
	for (int n = 0; n < (int)Site_types.size(); n++) {
		if (n > 0) {
			outfile << ",";
		}
		outfile << "Composition" << (int)Site_types[n];
	}
	outfile << endl;
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			counts.assign(Site_types.size(), 0);
			for (int z = 0; z < lattice.getHeight(); z++) {
				counts[getSiteTypeIndex(lattice.getSiteType(x, y, z))]++;
			}
			//outfile << x << "," << y;
			for (int n = 0; n < (int)Site_types.size(); n++) {
				if (n > 0) {
					outfile << ",";
				}
				outfile << (double)counts[n] / (double)lattice.getHeight();
			}
			outfile << endl;
		}
	}
}

void Morphology::outputCorrelationData(ofstream& outfile) const {
	// Output Calculation Results
	int max_correlation_size = (int)Correlation_data[0].size();
	for (int n = 1; n < (int)Site_types.size(); n++) {
		if ((int)Correlation_data[n].size() > max_correlation_size) {
			max_correlation_size = (int)Correlation_data[n].size();
		}
	}
	outfile << "Distance (nm)";
	for (int n = 0; n < (int)Site_types.size(); n++) {
		outfile << ",Correlation" << (int)Site_types[n];
	}
	outfile << endl;
	for (int i = 0; i < max_correlation_size; i++) {
		if (i < (int)Correlation_data[0].size()) {
			outfile << lattice.getUnitSize()*0.5*(double)i << "," << Correlation_data[0][i];
		}
		else {
			outfile << lattice.getUnitSize()*0.5*(double)i << "," << NAN;
		}
		for (int n = 1; n < (int)Site_types.size(); n++) {
			if (i < (int)Correlation_data[n].size()) {
				outfile << "," << Correlation_data[n][i];
			}
			else {
				outfile << "," << NAN;
			}
		}
		outfile << endl;
	}
}

void Morphology::outputDepthDependentData(ofstream& outfile) const {
	outfile << "Z-Position";
	for (int n = 0; n < (int)Site_types.size(); n++) {
		outfile << ",Type" << (int)Site_types[n] << "_composition";
	}
	for (int n = 0; n < (int)Site_types.size(); n++) {
		outfile << ",Type" << (int)Site_types[n] << "_IV_fraction";
	}
	for (int n = 0; n < (int)Site_types.size(); n++) {
		outfile << ",Type" << (int)Site_types[n] << "_domain_size";
	}
	outfile << endl;
	for (int z = 0; z < lattice.getHeight(); z++) {
		outfile << z;
		for (int n = 0; n < (int)Site_types.size(); n++) {
			outfile << "," << Depth_composition_data[n][z];
		}
		for (int n = 0; n < (int)Site_types.size(); n++) {
			outfile << "," << Depth_iv_data[n][z];
		}
		for (int n = 0; n < (int)Site_types.size(); n++) {
			outfile << "," << Depth_domain_size_data[n][z];
		}
		outfile << endl;
	}
}

void Morphology::outputMorphologyCrossSection(ofstream& outfile) const {
	int x = lattice.getLength() / 2;
	//for (int x = 0; x < lattice.getLength(); x++) {
	for (int y = 0; y < lattice.getWidth(); y++) {
		//int z = lattice.getHeight() / 2;
		for (int z = 0; z < lattice.getHeight(); z++) {
			outfile << x << "," << y << "," << z << "," << (int)lattice.getSiteType(x, y, z) << endl;
		}
	}
	//}
}

void Morphology::outputMorphologyFile(string version, ofstream& outfile, bool enable_export_compressed) const {
	if (enable_export_compressed) {
		outfile << "Ising_OPV " << version << " - compressed format" << endl;
	}
	else {
		outfile << "Ising_OPV " << version << " - uncompressed format" << endl;
	}
	outfile << lattice.getLength() << endl;
	outfile << lattice.getWidth() << endl;
	outfile << lattice.getHeight() << endl;
	outfile << lattice.isXPeriodic() << endl;
	outfile << lattice.isYPeriodic() << endl;
	outfile << lattice.isZPeriodic() << endl;
	outfile << (int)Site_types.size() << endl;
	for (int n = 0; n < (int)Site_types.size(); n++) {
		outfile << Domain_sizes[n] << endl;
	}
	for (int n = 0; n < (int)Site_types.size(); n++) {
		outfile << Mix_fractions[n] << endl;
	}
	if (!enable_export_compressed) {
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					outfile << x << "," << y << "," << z << "," << (int)lattice.getSiteType(x, y, z) << endl;
				}
			}
		}
	}
	else {
		int count = 1;
		char previous_type = lattice.getSiteType(0, 0, 0);
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					if (x == 0 && y == 0 && z == 0) {
						continue;
					}
					if (lattice.getSiteType(x, y, z) == previous_type) {
						count++;
					}
					else {
						outfile << (int)previous_type << count << "\n";
						count = 1;
						previous_type = lattice.getSiteType(x, y, z);
					}
				}
			}
		}
		outfile << (int)previous_type << count << "\n";
	}
}

void Morphology::outputTortuosityMaps(ofstream& outfile) const {
	int index;
	//outfile << "X-Position,Y-Position";
	for (int n = 0; n < (int)Site_types.size(); n++) {
		if (n > 0) {
			outfile << ",";
		}
		outfile << "Tortuosity" << (int)Site_types[n];
	}
	outfile << endl;
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			index = lattice.getWidth()*x + y;
			//outfile << x << "," << y;
			for (int n = 0; n < Site_types.size(); n++) {
				if (n > 0) {
					outfile << ",";
				}
				outfile << Tortuosity_data[n][index];
			}
			outfile << endl;
		}
	}
}

void Morphology::shrinkLattice(int rescale_factor) {
	// Error handling
	if (rescale_factor == 0) {
		cout << "Error! Lattice cannot be shrunken by a rescale factor of zero." << endl;
		return;
	}
	if (lattice.getLength() % rescale_factor != 0 || lattice.getWidth() % rescale_factor != 0 || lattice.getHeight() % rescale_factor != 0) {
		cout << "Error! All lattice dimensions are not divisible by the rescale factor." << endl;
		return;
	}
	// Construct the smaller lattice 
	Lattice lattice_rescale = lattice;
	lattice_rescale.resize(lattice.getLength() / rescale_factor, lattice.getWidth() / rescale_factor, lattice.getHeight() / rescale_factor);
	// Assign site types to the new lattice based on the existing lattice
	int type1_count;
	bool alternate = true;
	for (int x = 0; x < lattice_rescale.getLength(); x++) {
		for (int y = 0; y < lattice_rescale.getWidth(); y++) {
			for (int z = 0; z < lattice_rescale.getHeight(); z++) {
				type1_count = 0;
				for (int i = rescale_factor * x; i < (rescale_factor*x + rescale_factor); i++) {
					for (int j = rescale_factor * y; j < (rescale_factor*y + rescale_factor); j++) {
						for (int k = rescale_factor * z; k < (rescale_factor*z + rescale_factor); k++) {
							if (lattice.getSiteType(i, j, k) == (char)1) {
								type1_count++;
							}
						}
					}
				}
				if (2 * type1_count > (rescale_factor*rescale_factor*rescale_factor)) {
					lattice_rescale.setSiteType(x, y, z, (char)1);
				}
				else if (2 * type1_count < (rescale_factor*rescale_factor*rescale_factor)) {
					lattice_rescale.setSiteType(x, y, z, (char)2);
				}
				else {
					if (alternate) {
						lattice_rescale.setSiteType(x, y, z, (char)1);
					}
					else {
						lattice_rescale.setSiteType(x, y, z, (char)2);
					}
					alternate = !alternate;
				}
			}
		}
	}
	// Update the lattice
	lattice = lattice_rescale;
	// The shrinking process can change the mix fraction, so the Mix_fraction property is updated.
	calculateMixFractions();
}

void Morphology::stretchLattice(int rescale_factor) {
	// Construct the larger lattice
	Lattice lattice_rescale = lattice;
	lattice_rescale.resize(lattice.getLength()*rescale_factor, lattice.getWidth()*rescale_factor, lattice.getHeight()*rescale_factor);
	// Assign site types to the new lattice based on the existing lattice
	for (int x = 0; x < lattice.getLength(); x++) {
		for (int y = 0; y < lattice.getWidth(); y++) {
			for (int z = 0; z < lattice.getHeight(); z++) {
				for (int i = rescale_factor * x; i < rescale_factor*x + rescale_factor; i++) {
					for (int j = rescale_factor * y; j < rescale_factor*y + rescale_factor; j++) {
						for (int k = rescale_factor * z; k < rescale_factor*z + rescale_factor; k++) {
							lattice_rescale.setSiteType(i, j, k, lattice.getSiteType(x, y, z));
						}
					}
				}
			}
		}
	}
	// Update the lattice
	lattice = lattice_rescale;
}

double Morphology::rand01() {
	return generate_canonical<double, std::numeric_limits<double>::digits>(gen);
}

void Morphology::updateNeighborCounts(long int site_index1, long int site_index2) {
	char site_type1 = lattice.getSiteType(site_index1);
	char site_type2 = lattice.getSiteType(site_index2);
	long int neighbor_index;
	Neighbor_counts[site_index1] = Temp_counts1;
	Neighbor_counts[site_index2] = Temp_counts2;
	for (int i = 0; i < 6; i++) {
		neighbor_index = Neighbor_info[site_index1].first_indices[i];
		if (neighbor_index >= 0 && neighbor_index != site_index2) {
			if (lattice.getSiteType(neighbor_index) == site_type1) {
				Neighbor_counts[neighbor_index].sum1++;
			}
			else {
				Neighbor_counts[neighbor_index].sum1--;
			}
		}
	}
	for (int i = 0; i < 6; i++) {
		neighbor_index = Neighbor_info[site_index2].first_indices[i];
		if (neighbor_index >= 0 && neighbor_index != site_index1) {
			if (lattice.getSiteType(neighbor_index) == site_type2) {
				Neighbor_counts[neighbor_index].sum1++;
			}
			else {
				Neighbor_counts[neighbor_index].sum1--;
			}
		}
	}
	for (int i = 0; i < 12; i++) {
		neighbor_index = Neighbor_info[site_index1].second_indices[i];
		if (neighbor_index >= 0) {
			if (lattice.getSiteType(neighbor_index) == site_type1) {
				Neighbor_counts[neighbor_index].sum2++;
			}
			else {
				Neighbor_counts[neighbor_index].sum2--;
			}
		}
	}
	for (int i = 0; i < 12; i++) {
		neighbor_index = Neighbor_info[site_index2].second_indices[i];
		if (neighbor_index >= 0) {
			if (lattice.getSiteType(neighbor_index) == site_type2) {
				Neighbor_counts[neighbor_index].sum2++;
			}
			else {
				Neighbor_counts[neighbor_index].sum2--;
			}
		}
	}
	if (Enable_third_neighbor_interaction) {
		for (int i = 0; i < 8; i++) {
			neighbor_index = Neighbor_info[site_index1].third_indices[i];
			if (neighbor_index >= 0 && neighbor_index != site_index2) {
				if (lattice.getSiteType(neighbor_index) == site_type1) {
					Neighbor_counts[neighbor_index].sum3++;
				}
				else {
					Neighbor_counts[neighbor_index].sum3--;
				}
			}
		}
		for (int i = 0; i < 8; i++) {
			neighbor_index = Neighbor_info[site_index2].third_indices[i];
			if (neighbor_index >= 0 && neighbor_index != site_index1) {
				if (lattice.getSiteType(neighbor_index) == site_type2) {
					Neighbor_counts[neighbor_index].sum3++;
				}
				else {
					Neighbor_counts[neighbor_index].sum3--;
				}
			}
		}
	}
}
