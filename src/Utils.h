// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <vector>

//! \brief This simple struct contains Cartesian coordinates specified by integers x,y,z.
//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
//! \author Michael C. Heiber
//! \date 2014-2018
struct Coords {
	//! The x Cartesian coordinate.
	int x;
	//! The y Cartesian coordinate.
	int y;
	//! The z Cartesian coordinate.
	int z;

	//! \brief Sets the x,y,z coordinates using the input values xval, yval, and zval.
	//! \param xval is the input x value.
	//! \param yval is the input y value.
	//! \param zval is the input z value.
	void setXYZ(const int xval, const int yval, const int zval) {
		x = xval;
		y = yval;
		z = zval;
	}

	//! \brief Defines the conditions when two Coords structs are deemed to be equal.
	//! \param rhs is the right hand side Coords struct of the comparison operator.
	bool operator==(const Coords& rhs) const {
		return (x == rhs.x && y == rhs.y && z == rhs.z);
	}

	//! \brief Defines the conditions when two Coords structs are deemed to be unequal.
	//! \param rhs is the right hand side Coords struct  of the comparison operator.
	bool operator!=(const Coords& rhs) const {
		return (x != rhs.x || y != rhs.y || z != rhs.z);
	}
};

//! \brief This namespace provides useful constants and utility functions.
//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
//! \author Michael C. Heiber
//! \date 2014-2018
namespace Utils {

	// Scientific Constants
	static constexpr double K_b = 8.61733035e-5;
	static constexpr double Elementary_charge = 1.602176621e-19; // C
	static constexpr double Vacuum_permittivity = 8.854187818e-12; // C/Vm
	static constexpr double Plank = 6.626070040e-34;
	static constexpr double Light_speed = 2.99792458e8;
	static constexpr double Pi = 3.14159265359;
	static constexpr double Coulomb_constant = 8.987551787e9; // N m^2 C^-2

	//! \brief Calculates the cumulative histogram from the input probability histogram data.
	//! \param data is the input probability histogram which is a vector of x-y pairs consisting of bin-centered x values and probability y values.
	//! \returns A vector of x-y pairs consisting of bin-centered x values and cumnulative probability y values.
	std::vector<std::pair<double, double>> calculateCumulativeHist(const std::vector<std::pair<double, double>>& data);

	//! \brief Calculates the probability histogram for the input double data vector using the input number of bins.
	//! \details Linearly spaced bins are automatically created from the minimum value to the maximum value of the data set. 
	//! The function outputs bin-centered x values and probability y values in a x-y pair vector.
	//! \param data is the input data vector.
	//! \param num_bins is the desired number of bins.
	//! \returns A vector of x-y pairs consisting of bin-centered x values and probability y values.
	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<double>& data, int num_bins);

	//! \brief Calculates the probability histogram for the input integer data vector using the input bin size.
	//! \details Linearly spaced bins are automatically created from the minimum value to the maximum value of the data set with the specified bin size.
	//! The function outputs bin-centered x values and probability y values in a x-y pair vector.
	//! \param data is the input integer data vector.
	//! \param bin_size is the specified bin size
	//! \returns A vector of x-y pairs consisting of bin-centered x values and probability y values.
	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<int>& data, int bin_size);

	//! \brief Calculates the probability histogram for the input data vector using the input bin size.
	//! \details Linearly spaced bins are automatically created from the minimum value to the maximum value of the data set. 
	//! with the specified bin spacing. The function outputs bin-centered x values and probability y values in a x-y pair vector.
	//! \param data is the input data vector.
	//! \param bin_size is the desired bin size.
	//! \returns A vector of x-y pairs consisting of bin-centered x values and probability y values.
	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<double>& data, double bin_size);

	//! \brief Calculates the probability histogram for the input data vector using the input bin size and input number of bins.
	//! \details Linearly spaced bins are automatically created starting from the minimum value of the data set. The function
	//! outputs bin-centered x values and probability y values in a x-y pair vector.
	//! \param data is the input data vector.
	//! \param num_bins is the number of bins that will be created.
	//! \param bin_size is the input bin size.
	//! \returns A vector of x-y pairs consisting of bin-centered x values and probability y values.
	std::vector<std::pair<double, double>> calculateProbabilityHist(const std::vector<double>& data, const double bin_size, const int num_bins);

	//! \brief Extracts a boolean value from a string containing "true" or "false".
	//! \param input is the input string.
	//! \param error_status is an input boolean that is used to indicate an error with the import process.
	//! \return true if the input string is "true".
	//! \return flase if the input string is "flase".
	bool importBooleanParam(const std::string& input, bool& error_status);

	//! \brief Numerically integrates a vector of x-y data using the trapezoid rule.
	//! \warning The function assumes that the data is sorted by the x values.
	//! \param data is the data vector containing x-y data pairs.
	//! \return the numerically calculated area.
	double integrateData(const std::vector<std::pair<double, double>>& data);

	//! \brief Linearly interpolates an x-y data set to determine the interpolated y-value corresponding to an input x-value.
	//! \warning The function assumes that the data is sorted by the x values.
	//! \param data is the data vector containing x-y data pairs.
	//! \param x_val is the x-value that will be interpolated to.
	//! \return the interpolated y-value when the input x-value lies within the range of the input data.
	double interpolateData(const std::vector<std::pair<double, double>>& data, const double x_val);

	//! \brief Uses MPI to calculate the element-wise average vector from separate vectors coming from different processors.
	//! \details Each processor calls this function and sends an input vector.  Upon function return, processor 0 
	//! receives the average vector and all of the other processors receive an empty vector.
	//! \param input_vector is the input data from the processor calling the function.
	//! \return A vector that is the element-wise average of all input vectors from each processor, when called on processor 0.
	//! \return An empty vector when called on other processors.
	std::vector<double> MPI_calculateVectorAvg(const std::vector<double>& input_vector);

	//! \brief Uses MPI to calculate the element-wise sum vector from separate vectors coming from different processors.
	//! \details Each processor calls this function and sends an input vector.  Upon function return, processor 0 
	//! receives the sum vector and all of the other processors receive an empty vector.
	//! \param input_vector is the input data from the processor calling the function.
	//! \return A vector that is the element-wise sum of all input vectors from each processor, when called on processor 0.
	//! \return An empty vector when called on other processors.
	std::vector<double> MPI_calculateVectorSum(const std::vector<double>& input_vector);

	//! \brief Uses MPI to calculate the element-wise sum vector from separate vectors coming from different processors.
	//! \details Each processor calls this function and sends an input vector.  Upon function return, processor 0 
	//! receives the sum vector and all of the other processors receive an empty vector.
	//! \param input_vector is the input data from the processor calling the function.
	//! \return A vector that is the element-wise sum of all input vectors from each processor, when called on processor 0.
	//! \return An empty vector when called on other processors.
	std::vector<int> MPI_calculateVectorSum(const std::vector<int>& input_vector);

	//! \brief Uses MPI to gather values from separate processors to build one vector containing all of the data.
	//! \details Each processor calls this function and sends an input value.  Upon function return, processor 0
	//! receives the data vector and all of the other processors receive an empty vector.
	//! \param input_val is the input data value from the processor calling the function.
	//! \return A vector that is a collection of the values from each processor, when called on processor 0.
	//! \return An empty vector when called on other processors.
	std::vector<double> MPI_gatherValues(const double input_val);

	//! \brief Uses MPI to gather vectors from separate processors to build one big vector containing all of the data.
	//! \details Each processor calls this function and sends an input vector.  Upon function return, processor 0
	//! receives the large data vector and all of the other processors receive an empty vector.
	//! \param input_vector is the input data from the processor calling the function.
	//! \return A vector that is a concatenation of all input vectors from each processor, when called on processor 0.
	//! \return An empty vector when called on other processors.
	std::vector<double> MPI_gatherVectors(const std::vector<double>& input_vector);

	//! \brief Uses MPI to gather vectors from separate processors to build one big vector containing all of the data.
	//! \details Each processor calls this function and sends an input vector.  Upon function return, processor 0
	//! receives the large data vector and all of the other processors receive an empty vector.
	//! \param input_vector is the input data from the processor calling the function.
	//! \return A vector that is a concatenation of all input vectors from each processor, when called on processor 0.
	//! \return An empty vector when called on other processors.
	std::vector<int> MPI_gatherVectors(const std::vector<int>& input_vector);

	//! \brief Removes all spaces and ta characters in a string.
	//! \param str is the input string
	//! \returns a new string that will have the whitespace removed.
	std::string removeWhitespace(const std::string& str);

	//! \brief Rounds the input double value to the nearest integer value.
	//! \param num is the input double value.
	//! \returns an integer value that is the rounded value.
	int round_int(const double num);

	// Template functions

	//! \brief This template function calculates and returns the average value in double format when given an array of numerical datatypes.
	//! \param data is the array of numerical data.
	//! \param array_size is the size of the input data array.
	//! \return The average of the data set in double format.
	template<typename T>
	double array_avg(const T data[], const int array_size) {
		double sum = 0;
		for (int i = 0; i < array_size; i++) {
			sum += data[i];
		}
		return sum / array_size;
	}

	//! \brief This template function calculates and returns the median value when given an array of numerical datatypes.
	//! \param data is the array of numerical data.
	//! \param array_size is the size of the input data array.
	//! \return The median value of the data set.
	template<typename T>
	T array_median(const T data[], const int array_size) {
		std::vector<T> data_vect(array_size);
		for (int i = 0; i < array_size; i++) {
			data_vect[i] = data[i];
		}
		sort(data_vect.begin(), data_vect.end());
		return data_vect[array_size / 2];
	}

	//! \brief This template function calculates and returns the standard deviation in double format when given an array of numerical datatypes.
	//! \param data is the input array of numerical data.
	//! \param array_size is the size of the input data array.
	//! \return The standard deviation of the data set in double format.
	template<typename T>
	double array_stdev(const T data[], const int array_size) {
		double sum = 0;
		double avg = array_avg(data, array_size);
		for (int i = 0; i < array_size; i++) {
			sum += (data[i] - avg)*(data[i] - avg);
		}
		return sqrt(sum / (array_size - 1));
	}

	//! \brief This template function determines which element in the data array has the median value of the dataset.
	//! \param data is the input array of numerical data.
	//! \param array_size is the size of the input data array.
	//! \return the integer array index that indicates the median value.
	//! \return -1 if the median cannot be determined.
	template<typename T>
	int array_which_median(const T data[], const int array_size) {
		T median = array_median(data, array_size);
		for (int i = 0; i < array_size; i++) {
			if (data[i] == median) {
				return i;
			}
		}
		return -1;
	}

	//! \brief This template function calculates and returns the results of an integer power operation on a base numerical datatype.
	//! \param base is the base of the integer power operation and can be any numerical type.
	//! \param exponent is the integer exponent for the integer power operation.
	//! \return The result of the integer power operation in whichever datatype format was input with the base.
	template<typename base_type>
	base_type intpow(const base_type base, const int exponent) {
		base_type result;
		if (exponent == 0) {
			result = 1;
			return result;
		}
		result = base;
		for (int i = 1; i < abs(exponent); i++) {
			result *= base;
		}
		if (exponent < 0) {
			result = 1.0 / result;
		}
		return result;
	}

	//! \brief This template function outputs the input data vector to a file with the specified filename.
	//! \warning This function may overwrite existing files if not used carefully.
	//! \param vec is the input data vector.
	//! \param filename is the input file name.
	template<typename T>
	void outputVectorToFile(std::vector<T>& vec, std::string filename) {
		std::ofstream outfile;
		outfile.open(filename);
		for (int i = 0; i < (int)vec.size(); i++) {
			outfile << vec[i] << "\n";
		}
		outfile.close();
	}

	//! \brief This template function outputs the input data pair vector to a file with the specified filename.
	//! \warning This function may overwrite existing files if not used carefully.
	//! \param vec is the input data pair vector.
	//! \param filename is the input file name.
	template<typename T>
	void outputVectorToFile(std::vector<std::pair<T, T>>& vec, std::string filename) {
		std::ofstream outfile;
		outfile.open(filename);
		for (int i = 0; i < (int)vec.size(); i++) {
			outfile << vec[i].first << "," << vec[i].second << "\n";
		}
		outfile.close();
	}

	//! \brief This template function efficienctly removes the duplicate entries from an input vector.
	//! \details This algorithm allow efficient removal of duplicate vector objects when > or < comparison operators do not exist.
	//! \param vec is the input vector to be operated on.
	template<typename T>
	void removeDuplicates(std::vector<T> &vec) {
		if ((int)vec.size() > 2) {
			for (auto it1 = vec.begin(); it1 != vec.end() - 2; ++it1) {
				for (auto it2 = it1 + 2; it2 != vec.end(); ++it2) {
					if (*it2 == *it1) {
						*it2 = *(it2 - 1);
					}
				}
			}
		}
		if ((int)vec.size() > 1) {
			auto it3 = unique(vec.begin(), vec.end());
			vec.resize(distance(vec.begin(), it3));
		}
	}

	//! \brief This template function calculates and returns the average value in double format when given a vector of numerical datatypes.
	//! \param data is the vector of numerical data.
	//! \return The average of the data set in double format.
	template<typename T, typename A>
	double vector_avg(const std::vector<T, A>& data) {
		double sum = 0;
		for (auto it = data.begin(); it != data.end(); ++it) {
			sum += *it;
		}
		return sum / data.size();
	}

	//! \brief This template function calculates and returns the median value when given a vector of numerical datatypes.
	//! \param data is the vector of numerical data.
	//! \return The median value in the data set.
	template<typename T, typename A>
	double vector_median(std::vector<T, A> data) {
		if (data.size() % 2 == 0) {
			const auto median_it1 = data.begin() + data.size() / 2 - 1;
			const auto median_it2 = data.begin() + data.size() / 2;
			std::nth_element(data.begin(), median_it1, data.end());
			T val1 = *median_it1;
			std::nth_element(data.begin(), median_it2, data.end());
			T val2 = *median_it2;
			return (double)(val1 + val2) / 2.0;
		}
		else {
			const auto median_it = data.begin() + data.size() / 2;
			std::nth_element(data.begin(), median_it, data.end());
			return *median_it;
		}
	}

	//! \brief This template function calculates and returns the standard deviation in double format when given a vector of numerical datatypes.
	//! \param data is the vector of numerical data.
	//! \return The standard deviation of the data set in double format.
	template<typename T, typename A>
	double vector_stdev(const std::vector<T, A>& data) {
		double sum = 0;
		double avg = vector_avg(data);
		for (auto it = data.begin(); it != data.end(); ++it) {
			sum += (*it - avg)*(*it - avg);
		}
		return sqrt(sum / (data.size() - 1));
	}

	//! \brief This template function determines which element in the data vector is closest to the median value of the dataset.
	//! \details if multiple values are equally close, the first one in the vector is selected.
	//! \param data is the input vector of numerical data.
	//! \return the integer vector index that indicates the median value.
	template<typename T, typename A>
	int vector_which_median(const std::vector<T, A>& data) {
		double median = vector_median(data);
		auto diff = data;
		for (auto& item : diff) {
			item = fabs(item - median);
		}
		auto it = min_element(diff.begin(), diff.end());
		return distance(diff.begin(),it);
	}
}

#endif // UTILS_H
