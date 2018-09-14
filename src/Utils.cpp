// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "Utils.h"

namespace Utils {

	using namespace std;

	std::vector<std::pair<double, double>> calculateCumulativeHist(const std::vector<std::pair<double, double>>& hist) {
		auto result = hist;
		for (int i = 1; i < (int)hist.size(); i++) {
			result[i].second = result[i - 1].second + hist[i].second;
		}
		return result;
	}

	std::vector<std::pair<double, int>> calculateHist(const std::vector<int>& data, int bin_size) {
		// Check for valid input data
		if ((int)data.size() == 0) {
			cout << "Error! Cannot calculate probability histogram because data vector is empty." << endl;
			vector<pair<double, int>> null_output = { { 0.0,0 } };
			return null_output;
		}
		// Determine the starting bin position
		int min_val = *min_element(data.begin(), data.end());
		int max_val = *max_element(data.begin(), data.end());
		// Determine number of bins
		int num_bins = (int)ceil((double)(max_val - min_val + 1) / (double)bin_size);
		// Calculate bins
		vector<pair<double, int>> hist(num_bins, make_pair(0.0, 0));
		for (int i = 0; i < num_bins; i++) {
			hist[i].first = min_val + 0.5*(bin_size - 1) + (bin_size - 1) * i;
		}
		// Calculate histogram
		int index;
		for (int i = 0; i < (int)data.size(); i++) {
			index = (data[i] - min_val) / bin_size;
			hist[index].second++;
		}
		return hist;
	}

	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<std::pair<double, int>> hist) {
		// Check for valid input data
		if ((int)hist.size() == 0) {
			cout << "Error! Cannot calculate probability histogram because the input histogram is empty." << endl;
			vector<pair<double, double>> null_output = { { 0.0,0.0 } };
			return null_output;
		}
		// Add up the total counts in the histogram
		int total_counts = 0;
		for (const auto item : hist) {
			total_counts += item.second;
		}
		// Normalized histogram to get probability histogram
		vector<pair<double, double>> prob_hist(hist.size(), make_pair(0.0, 0.0));
		for (int i = 0; i < (int)hist.size(); i++) {
			prob_hist[i].first = hist[i].first;
			prob_hist[i].second = (double)hist[i].second / (double)(total_counts);
		}
		return prob_hist;
	}

	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<int>& data, int bin_size) {
		// Check for valid input data
		if ((int)data.size() == 0) {
			cout << "Error! Cannot calculate probability histogram because data vector is empty." << endl;
			vector<pair<double, double>> null_output = { { 0.0,0.0 } };
			return null_output;
		}
		// Determine the starting bin position
		int min_val = *min_element(data.begin(), data.end());
		int max_val = *max_element(data.begin(), data.end());
		// Determine number of bins
		int num_bins = (int)ceil((double)(max_val - min_val) / (double)bin_size);
		// Calculate bins
		vector<pair<double, double>> hist(num_bins, make_pair(0.0, 0.0));
		for (int i = 0; i < num_bins; i++) {
			hist[i].first = min_val + 0.5*(bin_size - 1) + (bin_size - 1) * i;
		}
		// Calculate histogram
		vector<int> counts(num_bins, 0);
		int index;
		for (int i = 0; i < (int)data.size(); i++) {
			index = (data[i] - min_val) / bin_size;
			counts[index]++;
		}
		// total counts
		int total_counts = accumulate(counts.begin(), counts.end(), 0);
		// Normalized histogram to get probability
		for (int i = 0; i < num_bins; i++) {
			hist[i].second = (double)counts[i] / (double)(total_counts);
		}
		return hist;
	}

	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<double>& data, int num_bins) {
		// Check for valid input data
		if ((int)data.size() == 0) {
			cout << "Error! Cannot calculate probability histogram because data vector is empty." << endl;
			std::vector<std::pair<double, double>> null_output = { { 0.0,0.0 } };
			return null_output;
		}
		// Determine data range
		double min_val = *min_element(data.begin(), data.end());
		double max_val = *max_element(data.begin(), data.end());
		// Limit the number of bins to the number of data entries
		if (num_bins > (int)data.size()) {
			num_bins = (int)data.size();
		}
		// Extend the range a little bit to ensure all data fits in the range
		min_val -= 1e-12*abs(min_val);
		max_val += 1e-12*abs(max_val);
		// Determine bin size
		double bin_size = (max_val - min_val) / num_bins;
		return calculateProbabilityHist(data, bin_size, num_bins);
	}

	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<double>& data, double bin_size) {
		// Check for valid input data
		if ((int)data.size() == 0) {
			cout << "Error! Cannot calculate probability histogram because data vector is empty." << endl;
			std::vector<std::pair<double, double>> null_output = { { 0.0,0.0 } };
			return null_output;
		}
		// Determine data range
		double min_val = *min_element(data.begin(), data.end());
		double max_val = *max_element(data.begin(), data.end());
		// Extend the range a little bit to ensure all data fits in the range
		min_val -= 1e-12*abs(min_val);
		max_val += 1e-12*abs(max_val);
		// Determine number of bins
		int num_bins = (int)ceil((max_val - min_val) / bin_size);
		// Limit the number of bins to the number of data entries
		if (num_bins > (int)data.size()) {
			num_bins = (int)data.size();
			bin_size = (max_val - min_val) / (double)num_bins;
		}
		return calculateProbabilityHist(data, bin_size, num_bins);
	}

	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<double>& data, const double bin_size, const int num_bins) {
		// Check for valid input data
		if ((int)data.size() == 0) {
			cout << "Error! Cannot calculate probability histogram because data vector is empty." << endl;
			vector<pair<double, double>> null_output = { { 0.0,0.0 } };
			return null_output;
		}
		// Determine the starting bin position
		double min_val = *min_element(data.begin(), data.end());
		// Extend the range a little bit to ensure all data fits in the range
		min_val -= 1e-12*abs(min_val);
		// Calculate bin-centered x values
		vector<pair<double, double>> hist(num_bins, make_pair(0.0, 0.0));
		for (int i = 0; i < num_bins; i++) {
			hist[i].first = min_val + 0.5*bin_size + bin_size * i;
		}
		// Calculate histogram
		vector<int> counts(num_bins, 0);
		int index;
		for (int i = 0; i < (int)data.size(); i++) {
			index = (int)floor((data[i] - min_val) / bin_size);
			counts[index]++;
		}
		// total counts
		int total_counts = accumulate(counts.begin(), counts.end(), 0);
		// Normalized histogram to get probability
		for (int i = 0; i < num_bins; i++) {
			hist[i].second = (double)counts[i] / (double)(total_counts);
		}
		return hist;
	}

	bool importBooleanParam(const std::string& input, bool& error_status) {
		string str = removeWhitespace(input);
		if (str.compare("true") == 0) {
			return true;
		}
		else if (str.compare("false") == 0) {
			return false;
		}
		else {
			cout << "Error importing boolean parameter." << endl;
			error_status = true;
			return false;
		}
	}

	double integrateData(const std::vector<std::pair<double, double>>& data) {
		double area = 0;
		for (int i = 1; i < (int)data.size(); i++) {
			area += ((data[i - 1].second + data[i].second) / 2.0)*(data[i].first - data[i - 1].first);
		}
		return area;
	}

	double interpolateData(const std::vector<std::pair<double, double>>& data, const double x_val) {
		for (int i = 1; i < (int)data.size(); i++) {
			if (data[i - 1].first < x_val && data[i].first > x_val) {
				return data[i - 1].second + ((data[i].second - data[i - 1].second) / (data[i].first - data[i - 1].first))*(x_val - data[i - 1].first);
			}
			if (abs(data[i].first - x_val) < 1e-6) {
				return data[i].second;
			}
		}
		cout << "Warning! The input x-value lies outside the range of the input data set." << endl;
		return NAN;
	}

	std::vector<std::pair<double, double>> MPI_calculateProbHistAvg(const std::vector<std::pair<double, int>>& input_hist) {
		if ((int)input_hist.size() < 2) {
			throw invalid_argument("Unable to calculate the average probability histogram because the input histogram must have more than one bin.");
		}
		int procid;
		int nproc;
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		// Determine the smallest bin and largest bin
		//double min_bin = 0;
		//double max_bin = 0;
		//double *min_bins = NULL;
		//double *max_bins = NULL;
		//if (procid == 0) {
		//	min_bins = new double[nproc];
		//	max_bins = new double[nproc];
		//}
		//min_bin = input_hist[0].first;
		//max_bin = input_hist.back().first;
		//MPI_Gather(&min_bin, 1, MPI_DOUBLE, min_bins, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//MPI_Gather(&max_bin, 1, MPI_DOUBLE, max_bins, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//double smallest_bin = min_bins[0];
		//for (int i = 1; i < nproc; i++) {
		//	if (min_bins[i] < smallest_bin) {
		//		smallest_bin = min_bins[i];
		//	}
		//}
		//double largest_bin = max_bins[0];
		//for (int i = 1; i < nproc; i++) {
		//	if (max_bins[i] > largest_bin) {
		//		largest_bin = max_bins[i];
		//	}
		//}
		// Determine the bin size
		double bin_size = input_hist[1].first - input_hist[0].first;
		// Gather the histogram sizes
		int data_size = 0;
		int *data_sizes = NULL;
		if (procid == 0) {
			data_sizes = new int[nproc];
		}
		data_size = (int)input_hist.size();
		MPI_Gather(&data_size, 1, MPI_INT, data_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// Determine the largest hist size
		int max_data_size = 0;
		if (procid == 0) {
			for (int i = 0; i < nproc; i++) {
				if (data_sizes[i] > max_data_size) {
					max_data_size = data_sizes[i];
				}
			}
		}
		MPI_Bcast(&max_data_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// Separate out the counts data from the histograms
		vector<int> counts(input_hist.size());
		for (int i = 0; i < (int)input_hist.size(); i++) {
			counts[i] = input_hist[i].second;
		}
		// Add zeroes padding to the end of the counts vector if needed so that all counts vectors are the same size
		counts.insert(counts.end(), max_data_size - counts.size(), 0);
		// Add up the counts from all processors
		auto counts_sum = MPI_calculateVectorSum(counts);
		// Create output probablilty histogram
		vector<pair<double, double>> prob_hist;
		if (procid == 0) {
			int total_counts = accumulate(counts_sum.begin(), counts_sum.end(), 0);
			for (int i = 0; i < max_data_size; i++) {
				prob_hist.push_back(make_pair(input_hist[0].first + bin_size * i, (double)counts_sum[i] / (double)total_counts));
			}
		}
		delete[] data_sizes;
		return prob_hist;
	}

	std::vector<double> MPI_calculateVectorAvg(const std::vector<double>& input_vector) {
		int data_size = 0;
		int data_count = 0;
		double *data = NULL;
		double *data_all = NULL;
		int *data_sizes = NULL;
		int *data_displacement = NULL;
		int max_data_size = 0;
		double average = 0;
		int procid;
		int nproc;
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		vector<double> output_vector;
		if (procid == 0) {
			data_sizes = new int[nproc];
		}
		data_size = (int)input_vector.size();
		data = new double[data_size];
		for (int i = 0; i < data_size; i++) {
			data[i] = input_vector[i];
		}
		MPI_Gather(&data_size, 1, MPI_INT, data_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < nproc; i++) {
				data_count += data_sizes[i];
			}
			data_all = new double[data_count];
			data_displacement = new int[nproc];
			data_displacement[0] = 0;
			for (int i = 1; i < nproc; i++) {
				data_displacement[i] = data_displacement[i - 1] + data_sizes[i - 1];
			}
		}
		MPI_Gatherv(data, data_size, MPI_DOUBLE, data_all, data_sizes, data_displacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < nproc; i++) {
				if (data_sizes[i] > max_data_size) {
					max_data_size = data_sizes[i];
				}
			}
			for (int i = 0; i < max_data_size; i++) {
				average = 0;
				for (int j = 0; j < nproc; j++) {
					if (i < data_sizes[j]) {
						average += data_all[data_displacement[j] + i];
					}
				}
				average = average / nproc;
				output_vector.push_back(average);
			}
		}
		delete[] data;
		delete[] data_all;
		delete[] data_sizes;
		delete[] data_displacement;
		return output_vector;
	}

	std::vector<double> MPI_calculateVectorSum(const std::vector<double>& input_vector) {
		int data_size = 0;
		double *data = NULL;
		double *sum = NULL;
		vector<double> output_vector;
		int procid;
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		data_size = (int)input_vector.size();
		data = new double[data_size];
		sum = new double[data_size];
		for (int i = 0; i < (int)input_vector.size(); i++) {
			data[i] = input_vector[i];
		}
		MPI_Reduce(data, sum, data_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < data_size; i++) {
				output_vector.push_back(sum[i]);
			}
		}
		delete[] data;
		delete[] sum;
		return output_vector;
	}

	std::vector<int> MPI_calculateVectorSum(const std::vector<int>& input_vector) {
		int data_size = 0;
		int *data = NULL;
		int *sum = NULL;
		vector<int> output_vector;
		int procid;
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		data_size = (int)input_vector.size();
		data = new int[data_size];
		sum = new int[data_size];
		for (int i = 0; i < (int)input_vector.size(); i++) {
			data[i] = input_vector[i];
		}
		MPI_Reduce(data, sum, data_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < data_size; i++) {
				output_vector.push_back(sum[i]);
			}
		}
		delete[] data;
		delete[] sum;
		return output_vector;
	}

	std::vector<double> MPI_gatherValues(double input_val) {
		double *data = NULL;
		vector<double> output_vector;
		int procid;
		int nproc;
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		data = new double[nproc];
		MPI_Gather(&input_val, 1, MPI_DOUBLE, data, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < nproc; i++) {
				output_vector.push_back(data[i]);
			}
		}
		delete[] data;
		return output_vector;
	}

	std::vector<double> MPI_gatherVectors(const std::vector<double>& input_vector) {
		int data_size = 0;
		int data_count = 0;
		double *data = NULL;
		double *data_all = NULL;
		int *data_sizes = NULL;
		int *data_displacement = NULL;
		vector<double> output_vector;
		int procid;
		int nproc;
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		if (procid == 0) {
			data_sizes = new int[nproc];
		}
		data_size = (int)input_vector.size();
		data = new double[data_size];
		for (int i = 0; i < (int)input_vector.size(); i++) {
			data[i] = input_vector[i];
		}
		MPI_Gather(&data_size, 1, MPI_INT, data_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < nproc; i++) {
				data_count += data_sizes[i];
			}
			data_all = new double[data_count];
			data_displacement = new int[nproc];
			data_displacement[0] = 0;
			for (int i = 1; i < nproc; i++) {
				data_displacement[i] = data_displacement[i - 1] + data_sizes[i - 1];
			}
		}
		MPI_Gatherv(data, data_size, MPI_DOUBLE, data_all, data_sizes, data_displacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < data_count; i++) {
				output_vector.push_back(data_all[i]);
			}
		}
		delete[] data;
		delete[] data_all;
		delete[] data_sizes;
		delete[] data_displacement;
		return output_vector;
	}

	std::vector<int> MPI_gatherVectors(const std::vector<int>& input_vector) {
		int data_size = 0;
		int data_count = 0;
		int *data = NULL;
		int *data_all = NULL;
		int *data_sizes = NULL;
		int *data_displacement = NULL;
		vector<int> output_vector;
		int procid;
		int nproc;
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		if (procid == 0) {
			data_sizes = new int[nproc];
		}
		data_size = (int)input_vector.size();
		data = new int[data_size];
		for (int i = 0; i < (int)input_vector.size(); i++) {
			data[i] = input_vector[i];
		}
		MPI_Gather(&data_size, 1, MPI_INT, data_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < nproc; i++) {
				data_count += data_sizes[i];
			}
			data_all = new int[data_count];
			data_displacement = new int[nproc];
			data_displacement[0] = 0;
			for (int i = 1; i < nproc; i++) {
				data_displacement[i] = data_displacement[i - 1] + data_sizes[i - 1];
			}
		}
		MPI_Gatherv(data, data_size, MPI_INT, data_all, data_sizes, data_displacement, MPI_INT, 0, MPI_COMM_WORLD);
		if (procid == 0) {
			for (int i = 0; i < data_count; i++) {
				output_vector.push_back(data_all[i]);
			}
		}
		delete[] data;
		delete[] data_all;
		delete[] data_sizes;
		delete[] data_displacement;
		return output_vector;
	}

	std::string removeWhitespace(const std::string& str_input) {
		// Remove tab characters
		string str_out = str_input;
		str_out.erase(remove(str_out.begin(), str_out.end(), '\t'), str_out.end());
		str_out.erase(remove(str_out.begin(), str_out.end(), ' '), str_out.end());
		return str_out;
	}

	int round_int(const double num) {
		return (num > 0.0) ? (int)(num + 0.5) : (int)(num - 0.5);
	}

}
