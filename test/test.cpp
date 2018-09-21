// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "gtest/gtest.h"
#include "Lattice.h"
#include "Morphology.h"
#include "Parameters.h"
#include "Utils.h"

using namespace std;
using namespace Utils;

namespace UtilsTests {

	TEST(UtilsTests, CoordsTests) {
		Coords coords{ 1, 1, 1 };
		EXPECT_EQ(1, coords.x);
		EXPECT_EQ(1, coords.y);
		EXPECT_EQ(1, coords.z);
		Coords coords2{ 0, 5, 10 };
		EXPECT_TRUE(coords != coords2);
		coords2.setXYZ(1, 1, 1);
		EXPECT_TRUE(coords == coords2);
	}

	TEST(UtilsTests, CalculateHistTests) {
		// Check behavior for an empty data set
		vector<int> data;
		EXPECT_THROW(calculateHist(data, 1), invalid_argument);
		// Create sample data
		data = { 0,0,2,1,2,4,3,1,0,4,2 };
		// Calculate the histogram
		auto hist = calculateHist(data, 1);
		// Check hist size
		EXPECT_EQ(5, (int)hist.size());
		// Check hist bins
		EXPECT_DOUBLE_EQ(0.0, hist[0].first);
		EXPECT_DOUBLE_EQ(1.0, hist[1].first);
		EXPECT_DOUBLE_EQ(2.0, hist[2].first);
		EXPECT_DOUBLE_EQ(3.0, hist[3].first);
		EXPECT_DOUBLE_EQ(4.0, hist[4].first);
		// Check hist values
		EXPECT_EQ(3, hist[0].second);
		EXPECT_EQ(2, hist[1].second);
		EXPECT_EQ(3, hist[2].second);
		EXPECT_EQ(1, hist[3].second);
		EXPECT_EQ(2, hist[4].second);
		// Caclulate the probability hist using the hist
		auto prob = calculateProbabilityHist(hist);
		// Check behavior if the hist is empty
		hist.clear();
		EXPECT_THROW(calculateProbabilityHist(hist), invalid_argument);
		// Check prob size
		EXPECT_EQ(5, (int)prob.size());
		// Check prob values
		EXPECT_DOUBLE_EQ(3.0 / 11.0, prob[0].second);
		EXPECT_DOUBLE_EQ(2.0 / 11.0, prob[1].second);
		EXPECT_DOUBLE_EQ(3.0 / 11.0, prob[2].second);
		EXPECT_DOUBLE_EQ(1.0 / 11.0, prob[3].second);
		EXPECT_DOUBLE_EQ(2.0 / 11.0, prob[4].second);
		// Check that the prob hist sums to 1
		auto cum_hist = calculateCumulativeHist(prob);
		EXPECT_DOUBLE_EQ(1.0, cum_hist.back().second);
		// Calculate the prob hist directly from the data vector
		prob = calculateProbabilityHist(data, 1);
		// Check prob size
		EXPECT_EQ(5, (int)prob.size());
		// Check prob values
		EXPECT_DOUBLE_EQ(3.0 / 11.0, prob[0].second);
		EXPECT_DOUBLE_EQ(2.0 / 11.0, prob[1].second);
		EXPECT_DOUBLE_EQ(3.0 / 11.0, prob[2].second);
		EXPECT_DOUBLE_EQ(1.0 / 11.0, prob[3].second);
		EXPECT_DOUBLE_EQ(2.0 / 11.0, prob[4].second);
		// Check that the prob hist sums to 1
		cum_hist = calculateCumulativeHist(prob);
		EXPECT_DOUBLE_EQ(1.0, cum_hist.back().second);
		// Check behavior how larger bin size
		// Create sample data
		data = { 0,0,2,1,2,4,3,1,0,4,2,5 };
		// Calculate the histogram
		hist = calculateHist(data, 2);
		// Check hist size
		EXPECT_EQ(3, (int)hist.size());
		// Check hist bins
		EXPECT_DOUBLE_EQ(0.5, hist[0].first);
		EXPECT_DOUBLE_EQ(2.5, hist[1].first);
		EXPECT_DOUBLE_EQ(4.5, hist[2].first);
		// Check hist values
		EXPECT_EQ(5, hist[0].second);
		EXPECT_EQ(4, hist[1].second);
		EXPECT_EQ(3, hist[2].second);
		// Check behavior with invalid bin size
		EXPECT_THROW(calculateHist(data, 0), invalid_argument);
		EXPECT_THROW(calculateHist(data, -1), invalid_argument);
	}

	TEST(UtilsTests, CalculateProbabilityHistTests) {
		vector<int> int_data;
		// Test with empty int data set that exception is thrown
		EXPECT_THROW(calculateProbabilityHist(int_data, 5), invalid_argument);
		// Generate a set of data from a uniform real distribution
		mt19937_64 gen(std::random_device{}());
		uniform_real_distribution<> dist(0, 100);
		vector<double> data((int)2e7, 0.0);
		for (int i = 0; i < (int)data.size(); i++) {
			data[i] = dist(gen);
		}
		// Calculate histogram with 9 bins
		auto prob = calculateProbabilityHist(data, 10);
		// Check for the correct number of bins
		EXPECT_EQ(10, (int)prob.size());
		// Check several random values from the uniform probability hist
		uniform_int_distribution<> dist2(0, 9);
		EXPECT_NEAR(10.0 / 100.0, prob[dist2(gen)].second, 2e-4);
		EXPECT_NEAR(10.0 / 100.0, prob[dist2(gen)].second, 2e-4);
		EXPECT_NEAR(10.0 / 100.0, prob[dist2(gen)].second, 2e-4);
		// Check that the prob hist sums to 1
		auto cum_hist = calculateCumulativeHist(prob);
		EXPECT_DOUBLE_EQ(1.0, cum_hist.back().second);
		// Calculate histogram with a bin size of 10.0
		prob = calculateProbabilityHist(data, 10.0);
		// Check for the correct number of bins
		EXPECT_EQ(10, (int)prob.size());
		// Check several random values from the uniform probability hist
		EXPECT_NEAR(10.0 / 100.0, prob[dist2(gen)].second, 2e-4);
		EXPECT_NEAR(10.0 / 100.0, prob[dist2(gen)].second, 2e-4);
		EXPECT_NEAR(10.0 / 100.0, prob[dist2(gen)].second, 2e-4);
		// Check that the prob hist sums to 1
		cum_hist = calculateCumulativeHist(prob);
		EXPECT_DOUBLE_EQ(1.0, cum_hist.back().second);
		// Clear data vector
		data.clear();
		// Check that empty double data vectors throw an exception
		EXPECT_THROW(calculateProbabilityHist(data, 10.0), invalid_argument);
		EXPECT_THROW(calculateProbabilityHist(data, 5);, invalid_argument);
		EXPECT_THROW(calculateProbabilityHist(data, 1.0, 5), invalid_argument);
		// Check behavior on a test dataset
		data = { 0.0, 1.0, 2.0, 3.0, 4.0 };
		prob = calculateProbabilityHist(data, 10);
		EXPECT_EQ(5, (int)prob.size());
		prob = calculateProbabilityHist(data, 0.1);
		EXPECT_EQ(5, (int)prob.size());
	}

	TEST(UtilsTests, Str2boolTests) {
		EXPECT_TRUE(str2bool("true"));
		EXPECT_TRUE(str2bool(" true  "));
		EXPECT_FALSE(str2bool("false	"));
		EXPECT_FALSE(str2bool("   false"));
		// Check that invalid input strings throw an exception as designed
		EXPECT_THROW(str2bool("blah"), invalid_argument);
		EXPECT_THROW(str2bool("	blah  "), invalid_argument);
	}

	TEST(UtilsTests, IntegrateDataTests) {
		vector<pair<double, double>> data_vec = { { 0.0,0.0 },{ 1.0,1.0 },{ 2.0,2.0 },{ 3.0,3.0 } };
		auto area = integrateData(data_vec);
		EXPECT_DOUBLE_EQ(area, 9.0 / 2.0);
	}

	TEST(UtilsTests, InterpolateDataTests) {
		vector<pair<double, double>> data_vec(100);
		for (int i = 0; i < (int)data_vec.size(); i++) {
			data_vec[i].first = 0.1*i;
			data_vec[i].second = exp(-data_vec[i].first / 2.85);
		}
		EXPECT_NEAR(1 / exp(1), interpolateData(data_vec, 2.85), 1e-4);
		for (int i = 0; i < (int)data_vec.size(); i++) {
			data_vec[i].first = 0.2*i;
			data_vec[i].second = 2.5*data_vec[i].first - 5.0;
		}
		EXPECT_NEAR(3.25, interpolateData(data_vec, 3.3), 1e-4);
		EXPECT_DOUBLE_EQ(7.5, interpolateData(data_vec, 5));
		// Attempt to interpolate beyond the data range
		EXPECT_TRUE(std::isnan(interpolateData(data_vec, -1.0)));
		EXPECT_TRUE(std::isnan(interpolateData(data_vec, 21.0)));

	}

	TEST(UtilsTests, RemoveWhitespaceTests) {
		string str = " text          ";
		EXPECT_EQ(removeWhitespace(str), "text");
		str = " text	";
		EXPECT_EQ(removeWhitespace(str), "text");
		str = "			text ";
		EXPECT_EQ(removeWhitespace(str), "text");
		str = "t e	xt";
		EXPECT_EQ(removeWhitespace(str), "text");
		str = "text";
		EXPECT_EQ(removeWhitespace(str), "text");
	}

	TEST(UtilsTests, ArrayStatsTests) {
		// positive ints
		int int_data[10];
		for (int i = 0; i < 10; i++) {
			int_data[i] = i + 1;
		}
		EXPECT_DOUBLE_EQ(5.5, array_avg(int_data, 10));
		EXPECT_NEAR(3.02765035409749, array_stdev(int_data, 10), 1e-14);
		// negative ints
		for (int i = 0; i < 10; i++) {
			int_data[i] = -(i + 1);
		}
		EXPECT_DOUBLE_EQ(-5.5, array_avg(int_data, 10));
		EXPECT_NEAR(3.02765035409749, array_stdev(int_data, 10), 1e-14);
		// positive doubles
		double double_data[10];
		for (int i = 0; i < 10; i++) {
			double_data[i] = 0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(2.75, array_avg(double_data, 10));
		EXPECT_NEAR(1.51382517704875, array_stdev(double_data, 10), 1e-14);
		// negative doubles
		for (int i = 0; i < 10; i++) {
			double_data[i] = -0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(-2.75, array_avg(double_data, 10));
		EXPECT_NEAR(1.51382517704875, array_stdev(double_data, 10), 1e-14);
	}

	TEST(UtilsTests, IntPowTests) {
		EXPECT_DOUBLE_EQ(1.0, intpow(2.5, 0));
		EXPECT_DOUBLE_EQ(2.5, intpow(2.5, 1));
		EXPECT_DOUBLE_EQ(6.25, intpow(2.5, 2));
		EXPECT_DOUBLE_EQ(9536.7431640625, intpow(2.5, 10));
		EXPECT_DOUBLE_EQ(0.4, intpow(2.5, -1));
		EXPECT_DOUBLE_EQ(0.16, intpow(2.5, -2));
		EXPECT_DOUBLE_EQ(1.048576e-4, intpow(2.5, -10));
		EXPECT_DOUBLE_EQ(1.0, intpow(15.04564, 0));
		EXPECT_DOUBLE_EQ(1e-21, intpow(1e-7, 3));
	}

	TEST(UtilsTests, RemoveDuplicatesTests) {
		vector<int> vec{ 0, 1, 1, 2, 3, 1, 4, 2 };
		removeDuplicates(vec);
		EXPECT_EQ(2, vec[2]);
		EXPECT_EQ(5, (int)vec.size());
		vec = { 0, 1, 2, 1 };
		removeDuplicates(vec);
		EXPECT_EQ(2, vec[2]);
		EXPECT_EQ(3, (int)vec.size());
		vector<double> vec2{ 0.0, 1.0, 1.0, 2.0, 3.0, 1.0, 4.0, 2.0 };
		removeDuplicates(vec2);
		EXPECT_DOUBLE_EQ(2.0, vec2[2]);
		EXPECT_EQ(5, (int)vec2.size());
		Coords coords1{ 1,2,3 };
		Coords coords2{ 1,2,3 };
		Coords coords3{ 4,5,6 };
		vector<Coords> vec3{ coords1, coords1, coords2, coords3, coords1, coords3, coords2 };
		removeDuplicates(vec3);
		EXPECT_EQ(4, vec3[1].x);
		EXPECT_EQ(2, (int)vec3.size());
		vector<int> vec4 = {};
		removeDuplicates(vec4);
		EXPECT_EQ(0, (int)vec4.size());
		vector<int> vec5 = { 0 };
		removeDuplicates(vec5);
		EXPECT_EQ(1, (int)vec5.size());
		vector<int> vec6 = { 0,0 };
		removeDuplicates(vec6);
		EXPECT_EQ(1, (int)vec6.size());
	}

	TEST(UtilsTests, VectorStatsTests) {
		// positive ints
		vector<int> int_data;
		int_data.assign(10, 0);
		for (int i = 0; i < 10; i++) {
			int_data[i] = i + 1;
		}
		EXPECT_DOUBLE_EQ(5.5, vector_avg(int_data));
		EXPECT_NEAR(3.02765035409749, vector_stdev(int_data), 1e-14);
		EXPECT_DOUBLE_EQ(5.5, vector_median(int_data));
		EXPECT_EQ(4, vector_which_median(int_data));
		// negative ints
		for (int i = 0; i < 10; i++) {
			int_data[i] = -(i + 1);
		}
		EXPECT_DOUBLE_EQ(-5.5, vector_avg(int_data));
		EXPECT_NEAR(3.02765035409749, vector_stdev(int_data), 1e-14);
		EXPECT_DOUBLE_EQ(-5.5, vector_median(int_data));
		EXPECT_EQ(4, vector_which_median(int_data));
		// positive doubles
		vector<double> double_data;
		double_data.assign(10, 0);
		for (int i = 0; i < 10; i++) {
			double_data[i] = 0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(2.75, vector_avg(double_data));
		EXPECT_NEAR(1.51382517704875, vector_stdev(double_data), 1e-14);
		EXPECT_DOUBLE_EQ(2.75, vector_median(double_data));
		EXPECT_EQ(4, vector_which_median(double_data));
		// negative doubles
		double_data.assign(10, 0);
		for (int i = 0; i < 10; i++) {
			double_data[i] = -0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(-2.75, vector_avg(double_data));
		EXPECT_NEAR(1.51382517704875, vector_stdev(double_data), 1e-14);
		EXPECT_DOUBLE_EQ(-2.75, vector_median(double_data));
		EXPECT_EQ(4, vector_which_median(double_data));
		// Check simple odd size vector median
		vector<double> data = { 3.0, 0.0, 1.0, 5.0, 10.0 };
		EXPECT_DOUBLE_EQ(3.0, vector_median(data));
		EXPECT_EQ(0, vector_which_median(data));
	}
}

namespace LatticeTests {

	class LatticeTest : public ::testing::Test {
	protected:
		Lattice_Params params_lattice;
		Lattice lattice;

		void SetUp() {
			// Setup params
			params_lattice.Enable_periodic_x = true;
			params_lattice.Enable_periodic_y = true;
			params_lattice.Enable_periodic_z = true;
			params_lattice.Length = 50;
			params_lattice.Width = 50;
			params_lattice.Height = 50;
			params_lattice.Unit_size = 1.0;
			// Initialize Lattice object
			lattice.init(params_lattice);
		}
	};

	TEST_F(LatticeTest, InitializationTests) {
		EXPECT_EQ(50, lattice.getLength());
		EXPECT_EQ(50, lattice.getWidth());
		EXPECT_EQ(50, lattice.getHeight());
		EXPECT_DOUBLE_EQ(1.0, lattice.getUnitSize());
		EXPECT_EQ((long int)50 * 50 * 50, lattice.getNumSites());
		EXPECT_DOUBLE_EQ(125000e-21, lattice.getVolume());
	}

	TEST_F(LatticeTest, CalculateDestCoordsTests) {
		Coords coords_i{ 49, 49, 49 };
		Coords coords_f;
		lattice.calculateDestinationCoords(coords_i, 2, -2, 0, coords_f);
		EXPECT_EQ(1, coords_f.x);
		EXPECT_EQ(47, coords_f.y);
		EXPECT_EQ(49, coords_f.z);
		coords_i.setXYZ(0, 0, 49);
		lattice.calculateDestinationCoords(coords_i, 1, -2, 2, coords_f);
		EXPECT_EQ(1, coords_f.x);
		EXPECT_EQ(48, coords_f.y);
		EXPECT_EQ(1, coords_f.z);
		coords_i.setXYZ(1, 48, 1);
		lattice.calculateDestinationCoords(coords_i, -2, 2, -3, coords_f);
		EXPECT_EQ(49, coords_f.x);
		EXPECT_EQ(0, coords_f.y);
		EXPECT_EQ(48, coords_f.z);
	}

	TEST_F(LatticeTest, PeriodicCrossingTests) {
		Coords coords_i{ 49, 49, 49 };
		Coords coords_f{ 0, 49, 48 };
		EXPECT_EQ(lattice.calculateDX(coords_i, coords_f), -50);
		EXPECT_EQ(lattice.calculateDY(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDZ(coords_i, coords_f), 0);
		coords_i = { 0, 49, 48 };
		coords_f = { 49, 49, 49 };
		EXPECT_EQ(lattice.calculateDX(coords_i, coords_f), 50);
		EXPECT_EQ(lattice.calculateDY(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDZ(coords_i, coords_f), 0);
		coords_i.setXYZ(0, 0, 49);
		coords_f.setXYZ(1, 49, 0);
		EXPECT_EQ(lattice.calculateDX(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDY(coords_i, coords_f), 50);
		EXPECT_EQ(lattice.calculateDZ(coords_i, coords_f), -50);
		coords_i = { 1, 49, 0 };
		coords_f = { 0, 0, 49 };
		EXPECT_EQ(lattice.calculateDX(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDY(coords_i, coords_f), -50);
		EXPECT_EQ(lattice.calculateDZ(coords_i, coords_f), 50);
		coords_i.setXYZ(4, 5, 6);
		coords_f.setXYZ(3, 6, 5);
		EXPECT_EQ(lattice.calculateDX(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDY(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDZ(coords_i, coords_f), 0);
	}

	TEST_F(LatticeTest, CheckMoveValidityTests) {
		params_lattice.Enable_periodic_x = false;
		params_lattice.Enable_periodic_y = false;
		params_lattice.Enable_periodic_z = true;
		Lattice lattice2;
		lattice2.init(params_lattice);
		EXPECT_FALSE(lattice2.isXPeriodic());
		EXPECT_FALSE(lattice2.isYPeriodic());
		EXPECT_TRUE(lattice2.isZPeriodic());
		Coords coords1{ 49, 0, 49 };
		EXPECT_FALSE(lattice2.checkMoveValidity(coords1, 0, 0, 0));
		EXPECT_FALSE(lattice2.checkMoveValidity(coords1, 1, 0, 0));
		EXPECT_FALSE(lattice2.checkMoveValidity(coords1, 0, -1, 0));
		EXPECT_TRUE(lattice2.checkMoveValidity(coords1, 0, 0, 1));
		EXPECT_TRUE(lattice2.checkMoveValidity(coords1, -1, 1, -1));
		params_lattice.Enable_periodic_z = false;
		lattice2.init(params_lattice);
		EXPECT_FALSE(lattice2.isZPeriodic());
		EXPECT_FALSE(lattice2.checkMoveValidity(coords1, 0, 0, 1));
	}

	TEST_F(LatticeTest, RandomSiteGenTests) {
		Coords coords;
		int N_test = 40000000;
		vector<int> xcoords(N_test);
		vector<int> ycoords(N_test);
		vector<int> zcoords(N_test);
		for (int i = 0; i < N_test; i++) {
			coords = lattice.generateRandomCoords();
			EXPECT_TRUE(coords.x >= 0 && coords.x < lattice.getLength());
			EXPECT_TRUE(coords.y >= 0 && coords.y < lattice.getWidth());
			EXPECT_TRUE(coords.z >= 0 && coords.z < lattice.getHeight());
			xcoords[i] = coords.x;
			ycoords[i] = coords.y;
			zcoords[i] = coords.z;
		}
		EXPECT_NEAR(24.5, vector_avg(xcoords), 2e-2);
		EXPECT_NEAR(lattice.getLength() / sqrt(12.0), vector_stdev(xcoords), 1e-2);
		EXPECT_NEAR(24.5, vector_avg(ycoords), 2e-2);
		EXPECT_NEAR(lattice.getWidth() / sqrt(12.0), vector_stdev(ycoords), 1e-2);
		EXPECT_NEAR(24.5, vector_avg(zcoords), 2e-2);
		EXPECT_NEAR(lattice.getHeight() / sqrt(12.0), vector_stdev(zcoords), 1e-2);
	}

	TEST_F(LatticeTest, LatticeDistanceTests) {
		Coords coords_i{ 49, 49, 49 };
		Coords coords_f{ 1, 49, 47 };
		EXPECT_EQ(8, lattice.calculateLatticeDistanceSquared(coords_i, coords_f));
		coords_i.setXYZ(4, 5, 6);
		coords_f.setXYZ(3, 6, 5);
		EXPECT_EQ(3, lattice.calculateLatticeDistanceSquared(coords_i, coords_f));
		coords_i.setXYZ(0, 49, 1);
		coords_f.setXYZ(48, 1, 49);
		EXPECT_EQ(12, lattice.calculateLatticeDistanceSquared(coords_i, coords_f));
	}

	TEST_F(LatticeTest, GetSiteTests) {
		Coords coords1, coords2;
		int index;
		for (int i = 0; i < 100; i++) {
			coords1 = lattice.generateRandomCoords();
			index = lattice.getSiteIndex(coords1);
			coords2 = lattice.getSiteCoords(index);
			EXPECT_EQ(coords1.x, coords2.x);
			EXPECT_EQ(coords1.y, coords2.y);
			EXPECT_EQ(coords1.z, coords2.z);
		}
		// Check request for coords given a site index that not in the lattice
		EXPECT_THROW(lattice.getSiteCoords(-1), invalid_argument);
		EXPECT_THROW(lattice.getSiteCoords(lattice.getLength()*lattice.getWidth()*lattice.getHeight()), out_of_range);
		// Check behavior of getSiteType
		coords1.setXYZ(10, 10, 10);
		// Check default site type
		EXPECT_EQ((char)0, lattice.getSiteType(coords1));
		EXPECT_EQ((char)0, lattice.getSiteType(0));
		// Check behavior of invalid site index
		EXPECT_THROW(lattice.getSiteType(-1), invalid_argument);
		EXPECT_THROW(lattice.getSiteType(50 * 50 * 50), out_of_range);
		// Check behavior of invalid coords
		coords1.setXYZ(-1, -1, -1);
		EXPECT_THROW(lattice.getSiteType(coords1), invalid_argument);
		coords1.setXYZ(50, 50, 50);
		EXPECT_THROW(lattice.getSiteType(coords1), out_of_range);
	}

	TEST_F(LatticeTest, ExtractSublatticeTests) {
		// Extract a smaller lattice
		Lattice lattice_new = lattice.extractSublattice(0, 20, 0, 20, 0, 20);
		// Check new lattice dimensions
		EXPECT_EQ(20, lattice_new.getLength());
		EXPECT_EQ(20, lattice_new.getWidth());
		EXPECT_EQ(20, lattice_new.getHeight());
		// Extract a smaller lattice
		lattice_new = lattice.extractSublattice(9, 40, 19, 30, 0, 30);
		// Check new lattice dimensions
		EXPECT_EQ(40, lattice_new.getLength());
		EXPECT_EQ(30, lattice_new.getWidth());
		EXPECT_EQ(30, lattice_new.getHeight());
		// Check that invalid input throws the correct excpetions
		// Check negative inputs
		EXPECT_THROW(lattice.extractSublattice(-1, 20, 0, 20, 0, 20), invalid_argument);
		EXPECT_THROW(lattice.extractSublattice(0, 20, -1, 20, 0, 20), invalid_argument);
		EXPECT_THROW(lattice.extractSublattice(0, 20, 0, 20, -1, 20), invalid_argument);
		EXPECT_THROW(lattice.extractSublattice(0, -20, 0, 20, 0, 20), invalid_argument);
		EXPECT_THROW(lattice.extractSublattice(0, 20, 0, -20, 0, 20), invalid_argument);
		EXPECT_THROW(lattice.extractSublattice(0, 20, 0, 20, 0, -20), invalid_argument);
		// Check out of range inputs
		// sublattice size extends beyond lattice
		EXPECT_THROW(lattice.extractSublattice(20, 50, 20, 20, 20, 20), out_of_range);
		EXPECT_THROW(lattice.extractSublattice(20, 20, 20, 50, 20, 20), out_of_range);
		EXPECT_THROW(lattice.extractSublattice(20, 20, 20, 20, 20, 50), out_of_range);
		// starting coordinates start outside the lattice
		EXPECT_THROW(lattice.extractSublattice(50, 20, 20, 20, 20, 20), out_of_range);
		EXPECT_THROW(lattice.extractSublattice(20, 20, 50, 20, 20, 20), out_of_range);
		EXPECT_THROW(lattice.extractSublattice(20, 20, 20, 20, 50, 20), out_of_range);
	}

}

namespace MorphologyTests {

	TEST(MorphologyTests, ConstructorTests) {
		// Test the default constructor
		Morphology morph;
		// Check that the object has the default ID number
		EXPECT_EQ(0, morph.getID());
		// Check that the lattice has default size of 0,0,0.
		EXPECT_EQ(0, morph.getLength());
		// Test simple constructor
		morph = Morphology(0);
		// Check that ID number is set properly
		EXPECT_EQ(0, morph.getID());
		// Check that the lattice has default size of 0,0,0.
		EXPECT_EQ(0, morph.getLength());
		// Test the standard constructor
		Parameters params;
		params.Length = 50;
		params.Width = 60;
		params.Height = 70;
		params.Enable_periodic_z = false;
		morph = Morphology(params, 1);
		// Check that ID number is set properly
		EXPECT_EQ(1, morph.getID());
		// Check that the lattice has the correct size.
		EXPECT_EQ(50, morph.getLength());
		EXPECT_EQ(60, morph.getWidth());
		EXPECT_EQ(70, morph.getHeight());
		// Check the lattice input constructor
		Lattice_Params params_lattice;
		params_lattice.Enable_periodic_x = true;
		params_lattice.Enable_periodic_y = true;
		params_lattice.Enable_periodic_z = false;
		params_lattice.Length = 20;
		params_lattice.Width = 30;
		params_lattice.Height = 40;
		params_lattice.Unit_size = 1.5;
		Lattice lattice;
		// Initialize Lattice object
		lattice.init(params_lattice);
		morph = Morphology(lattice, params, 2);
		// Check that ID number is set properly
		EXPECT_EQ(2, morph.getID());
		// Check that the lattice has the correct size.
		EXPECT_EQ(20, morph.getLength());
		EXPECT_EQ(30, morph.getWidth());
		EXPECT_EQ(40, morph.getHeight());
		// Check that the lattice has the correct unit size
		EXPECT_DOUBLE_EQ(1.5, morph.getUnitSize());
		// Check that the params are set properly

	}

	TEST(MorphologyTests, RandomMorphologyTests) {
		// Setup default parameters
		Parameters params;
		params.Length = 50;
		params.Width = 50;
		params.Height = 50;
		params.Enable_periodic_z = true;
		params.N_sampling_max = 50000;
		params.Enable_e_method = true;
		params.Enable_mix_frac_method = false;
		params.Enable_extended_correlation_calc = false;
		params.Extended_correlation_cutoff_distance = 3;
		// Initialize Morphology object
		Morphology morph(params, 0);
		vector<double> mix_fractions;
		mix_fractions.assign(2, 0.5);
		morph.createRandomMorphology(mix_fractions);
		// Check that the mix fractions were implemented properly
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.001);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.001);
		// Check correlation data before it has been calculated
		auto data1 = morph.getCorrelationData((char)1);
		auto data2 = morph.getCorrelationData((char)2);
		EXPECT_DOUBLE_EQ(0.0, data1[0]);
		EXPECT_DOUBLE_EQ(0.0, data2[0]);
		// Check correlation data request for invalid site type
		EXPECT_THROW(morph.getCorrelationData((char)3), invalid_argument);
		// Calculate the inital domain size values
		morph.calculateCorrelationDistances();
		double domain_size1 = morph.getDomainSize((char)1);
		double domain_size2 = morph.getDomainSize((char)2);
		// Check domain size request for invalid site type
		EXPECT_THROW(morph.getDomainSize((char)3), invalid_argument);
		// Check that random blend domain size is less than 2
		EXPECT_LT(domain_size1, 2.0);
		EXPECT_LT(domain_size2, 2.0);
		// Calculate the initial domain anisotropy
		morph.calculateAnisotropies();
		// Check that the random blend is isotropic
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)1), 0.05);
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)2), 0.05);
		// Check anisotropy request for invalid site type
		EXPECT_THROW(morph.getDomainAnisotropy((char)3), invalid_argument);
		// Try creating random morphology with invalid mix_fractions vectors
		mix_fractions[0] = -0.5;
		mix_fractions[1] = -0.5;
		EXPECT_THROW(morph.createRandomMorphology(mix_fractions), invalid_argument);
		mix_fractions[0] = 0.5;
		mix_fractions[1] = 0.2;
		EXPECT_THROW(morph.createRandomMorphology(mix_fractions), invalid_argument);
		mix_fractions[0] = 0.5;
		mix_fractions[1] = 0.7;
		EXPECT_THROW(morph.createRandomMorphology(mix_fractions), invalid_argument);
	}

	TEST(MorphologyTests, AnisotropicPhaseSeparationTests) {
		// Setup default parameters
		Parameters params;
		params.Length = 40;
		params.Width = 40;
		params.Height = 40;
		params.Enable_periodic_z = false;
		params.N_sampling_max = 100000;
		params.Enable_e_method = true;
		params.Enable_mix_frac_method = false;
		params.Enable_extended_correlation_calc = false;
		params.Extended_correlation_cutoff_distance = 3;
		Morphology morph(params, 0);
		vector<double> mix_fractions;
		mix_fractions.assign(2, 0.5);
		morph.createRandomMorphology(mix_fractions);
		// Check that the mix fractions were implemented properly
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.001);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.001);
		// Calculate the inital domain size values
		morph.calculateCorrelationDistances();
		double domain_size1_i = morph.getDomainSize((char)1);
		double domain_size2_i = morph.getDomainSize((char)2);
		// Check that random blend domain size is less than 2
		EXPECT_LT(domain_size1_i, 2.0);
		EXPECT_LT(domain_size2_i, 2.0);
		// Calculate the initial domain anisotropy
		morph.calculateAnisotropies();
		// Check that the initial random blend is isotropic
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)1), 0.05);
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)2), 0.05);
		// Perform some anisotropic phase separation that creates aligned structures out-of-plane in the z-direction
		morph.executeIsingSwapping(220, 0.35, 0.35, true, 3, 0.05);
		// Calculate the final domain size values
		morph.calculateCorrelationDistances();
		double domain_size1_f = morph.getDomainSize((char)1);
		double domain_size2_f = morph.getDomainSize((char)2);
		// Check that the domain size has increased
		EXPECT_LT(domain_size1_i, domain_size1_f);
		EXPECT_LT(domain_size2_i, domain_size2_f);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(5.0, domain_size1_f, 0.5);
		EXPECT_NEAR(5.0, domain_size2_f, 0.5);
		// Calculate the final domain anisotropy
		morph.calculateAnisotropies();
		// Check the approximate magnitude of the anisotropy factor
		EXPECT_NEAR(1.4, morph.getDomainAnisotropy((char)1), 0.25);
		EXPECT_NEAR(1.4, morph.getDomainAnisotropy((char)2), 0.25);
		// Reset morphology to a random blend
		morph.createRandomMorphology(mix_fractions);
		// Calculate the inital domain size values
		morph.calculateCorrelationDistances();
		domain_size1_i = morph.getDomainSize((char)1);
		domain_size2_i = morph.getDomainSize((char)2);
		// Check that random blend domain size is less than 2
		EXPECT_LT(domain_size1_i, 2.0);
		EXPECT_LT(domain_size2_i, 2.0);
		// Perform some anisotropic phase separation that creates aligned structures in the x-y plane
		morph.executeIsingSwapping(180, 0.4, 0.4, true, 3, -0.05);
		// Calculate the final domain size values
		morph.calculateCorrelationDistances();
		domain_size1_f = morph.getDomainSize((char)1);
		domain_size2_f = morph.getDomainSize((char)2);
		// Check that the domain size has increased
		EXPECT_LT(domain_size1_i, domain_size1_f);
		EXPECT_LT(domain_size2_i, domain_size2_f);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(5.0, domain_size1_f, 0.5);
		EXPECT_NEAR(5.0, domain_size2_f, 0.5);
		// Calculate the final domain anisotropy
		morph.calculateAnisotropies();
		// Check the approximate magnitude of the anisotropy factor
		EXPECT_NEAR(0.84, morph.getDomainAnisotropy((char)1), 0.125);
		EXPECT_NEAR(0.84, morph.getDomainAnisotropy((char)2), 0.125);
		// Reset morphology to a random blend
		morph.createRandomMorphology(mix_fractions);
		// Perform some anisotropic phase separation that creates aligned structures in the x-direction
		morph.executeIsingSwapping(200, 0.35, 0.35, true, 1, 0.05);
		// Calculate the domain anisotropy
		morph.calculateAnisotropies();
		auto anisotropy1 = morph.getDomainAnisotropy((char)1);
		auto anisotropy2 = morph.getDomainAnisotropy((char)2);
		// Reset morphology to a random blend
		morph.createRandomMorphology(mix_fractions);
		// Perform some anisotropic phase separation that creates aligned structures in the y-direction
		morph.executeIsingSwapping(200, 0.35, 0.35, true, 2, 0.05);
		// Calculate the domain anisotropy
		morph.calculateAnisotropies();
		// Check that anisotropic phase separation in the x- or y-directions produces approximately the same anisotropy metric
		EXPECT_NEAR(anisotropy1, morph.getDomainAnisotropy((char)1), 0.175);
		EXPECT_NEAR(anisotropy2, morph.getDomainAnisotropy((char)2), 0.175);
		// Test anisotropic morphologies on lattices that are too narrow
		params.Length = 10;
		params.Width = 10;
		params.Height = 25;
		morph = Morphology(params, 0);
		morph.createRandomMorphology(mix_fractions);
		// Perform some anisotropic phase separation that creates aligned structures in the x-y plane
		morph.executeIsingSwapping(500, 0.4, 0.4, true, 3, -0.1);
		// Check calculation of anisotropy with narrow lattice
		morph.calculateAnisotropies();
		// Calculation should have an error and result in default value of -1
		EXPECT_DOUBLE_EQ(-1.0, morph.getDomainAnisotropy((char)1));
		EXPECT_DOUBLE_EQ(-1.0, morph.getDomainAnisotropy((char)2));
		// Test anisotropic morphologies on lattices that are too thin
		params.Length = 25;
		params.Width = 25;
		params.Height = 10;
		params.Enable_periodic_z = true;
		morph = Morphology(params, 0);
		morph.createRandomMorphology(mix_fractions);
		// Perform some anisotropic phase separation that creates aligned structures in the x-y plane
		morph.executeIsingSwapping(600, 0.35, 0.35, true, 3, 0.05);
		// Check calculation of anisotropy with thin lattice
		morph.calculateAnisotropies();
		// Calculation should have an error and result in default value of -1
		EXPECT_DOUBLE_EQ(-1.0, morph.getDomainAnisotropy((char)1));
		EXPECT_DOUBLE_EQ(-1.0, morph.getDomainAnisotropy((char)2));
	}

	TEST(MorphologyTests, MorphologyAnalysisTests) {
		Parameters params;
		params.Length = 40;
		params.Width = 40;
		params.Height = 40;
		params.Enable_periodic_z = false;
		Morphology morph(params, 0);
		morph.createBilayerMorphology();
		// Calculate interfacial volume and interfacial area to volume ratio
		double iv_fraction = morph.calculateInterfacialVolumeFraction();
		double iav_ratio = morph.calculateInterfacialAreaVolumeRatio();
		// Check interfacial metrics for the well-defined bilayer morphology
		EXPECT_DOUBLE_EQ((double)(40 * 40 * 2) / (double)(40 * 40 * 40), iv_fraction);
		EXPECT_DOUBLE_EQ((double)(40 * 40) / (double)(40 * 40 * 40), iav_ratio);
		// Create a checkerboard morphology
		morph.createCheckerboardMorphology();
		// Check that the mix fractions were implemented properly
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.01);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.01);
		// Check the interfacial volume fraction
		EXPECT_DOUBLE_EQ(1.0, morph.calculateInterfacialVolumeFraction());
	}

	class MorphologyTest : public ::testing::Test {
	protected:
		static Morphology* morph_start;
		Parameters params;

		static void SetUpTestCase() {
			// Setup default parameters
			Parameters params;
			params.Length = 50;
			params.Width = 50;
			params.Height = 50;
			params.Enable_periodic_z = true;
			params.N_sampling_max = 50000;
			params.Enable_e_method = true;
			params.Enable_mix_frac_method = false;
			params.Enable_extended_correlation_calc = false;
			params.Extended_correlation_cutoff_distance = 3;
			// Initialize Morphology object
			morph_start = new Morphology(params, 0);
			vector<double> mix_fractions;
			mix_fractions.assign(2, 0.5);
			morph_start->createRandomMorphology(mix_fractions);
			// Perform some phase separation
			morph_start->executeIsingSwapping(290, 0.4, 0.4, false, 0, 0.0);
			// Calculate the domain size values
			morph_start->calculateCorrelationDistances();
		}

		static void TearDownTestCase() {
			delete morph_start;
			morph_start = NULL;
		}

		virtual void SetUp() {
			// Setup default parameters
			params.Length = 50;
			params.Width = 50;
			params.Height = 50;
			params.Enable_periodic_z = true;
			params.N_sampling_max = 50000;
			params.Enable_e_method = true;
			params.Enable_mix_frac_method = false;
			params.Enable_extended_correlation_calc = false;
			params.Extended_correlation_cutoff_distance = 3;
		}

		virtual void TearDown() { }
	};

	Morphology* MorphologyTest::morph_start = NULL;

	TEST_F(MorphologyTest, FixtureSetupTests) {
		// Check that the Morphology's Lattice object was setup properly
		EXPECT_EQ(50, morph_start->getLength());
		EXPECT_EQ(50, morph_start->getWidth());
		EXPECT_EQ(50, morph_start->getHeight());
		// Check that the mix fractions were implemented properly
		EXPECT_NEAR(0.5, morph_start->getMixFraction((char)1), 0.001);
		EXPECT_NEAR(0.5, morph_start->getMixFraction((char)2), 0.001);
	}

	TEST_F(MorphologyTest, DomainSizeTests) {
		Morphology morph = *morph_start;
		double domain_size1_i = morph.getDomainSize((char)1);
		double domain_size2_i = morph.getDomainSize((char)2);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(6.0, domain_size1_i, 0.5);
		EXPECT_NEAR(6.0, domain_size2_i, 0.5);
		// Try domain size calculation using more sampling sites
		params.N_sampling_max = 100000;
		morph.setParameters(params);
		morph.calculateCorrelationDistances();
		double domain_size1_f = morph.getDomainSize((char)1);
		double domain_size2_f = morph.getDomainSize((char)2);
		// Check that the domain size is almost the same
		EXPECT_NEAR(domain_size1_i, domain_size1_f, 0.02);
		EXPECT_NEAR(domain_size2_i, domain_size2_f, 0.02);
		// Try the extended correlation calculation
		params.Enable_extended_correlation_calc = true;
		params.Extended_correlation_cutoff_distance = 5;
		morph.setParameters(params);
		morph.calculateCorrelationDistances();
		domain_size1_f = morph.getDomainSize((char)1);
		domain_size2_f = morph.getDomainSize((char)2);
		// Check the length of the correlation data
		auto data1 = morph.getCorrelationData((char)1);
		auto data2 = morph.getCorrelationData((char)2);
		EXPECT_EQ(11, data1.size());
		EXPECT_EQ(11, data2.size());
		// Check that the domain size is the same
		EXPECT_NEAR(domain_size1_i, domain_size1_f, 0.02);
		EXPECT_NEAR(domain_size2_i, domain_size2_f, 0.02);
		// Calculate domain size using the regular mix fraction method
		params.Enable_e_method = false;
		params.Enable_mix_frac_method = true;
		params.Enable_extended_correlation_calc = false;
		params.Extended_correlation_cutoff_distance = 3;
		morph.setParameters(params);
		morph.calculateCorrelationDistances();
		domain_size1_f = morph.getDomainSize((char)1);
		domain_size2_f = morph.getDomainSize((char)2);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(6.0, domain_size1_f, 0.5);
		EXPECT_NEAR(6.0, domain_size2_f, 0.5);
		// Check domain anisotropy
		morph.calculateAnisotropies();
		// Check that the phase separated blend is isotropic
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)1), 0.1);
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)2), 0.1);
		// Calculate the depth dependent characteristics
		morph.calculateDepthDependentData();
		data1 = morph.getDepthDomainSizeData((char)1);
		data2 = morph.getDepthDomainSizeData((char)2);
		// Check the average depth dependent domain size 
		EXPECT_NEAR(domain_size1_i, vector_avg(data1), 0.5);
		EXPECT_NEAR(domain_size2_i, vector_avg(data2), 0.5);
		// Gather depth dependent composition data
		data1 = morph.getDepthCompositionData((char)1);
		data2 = morph.getDepthCompositionData((char)2);
		// Check the approximate composition at each depth
		for (auto item : data1) {
			EXPECT_NEAR(0.5, item, 0.15);
		}
		for (auto item : data2) {
			EXPECT_NEAR(0.5, item, 0.15);
		}
		// Calculate the bulk interfacial volume fraction
		double iv_frac_i = morph.calculateInterfacialVolumeFraction();
		// Check depth dependent interfacial volume fraction compared to bulk value
		auto data = morph.getDepthIVData();
		for (auto item : data) {
			EXPECT_NEAR(iv_frac_i, item, 0.1);
		}
		string line;
		// Output correlation data
		ofstream outfile2("correlation_data.txt");
		morph.outputCorrelationData(outfile2);
		outfile2.close();
		// Check that output looks valid
		ifstream infile2("correlation_data.txt");
		getline(infile2, line);
		// Check first column names line
		EXPECT_EQ("Distance (nm),Correlation1,Correlation2", line);
		// Check first data line
		getline(infile2, line);
		EXPECT_EQ("0,1,1", line);
		infile2.close();
		// Output depth dependent data
		ofstream outfile3("depth_data.txt");
		morph.outputDepthDependentData(outfile3);
		outfile3.close();
		// Check that output looks valid
		ifstream infile3("depth_data.txt");
		// Check first column names line
		getline(infile3, line);
		EXPECT_EQ("Z-Position,Type1_composition,Type2_composition,Type1_domain_size,Type2_domain_size,IV_fraction", line);
	}

	TEST_F(MorphologyTest, TortuosityTests) {
		Morphology morph = *morph_start;
		// Calculate the tortuosity
		morph.calculateTortuosity((char)1, false);
		morph.calculateTortuosity((char)2, false);
		// Check tortuosity calculation of invalid site type
		EXPECT_FALSE(morph.calculateTortuosity((char)3, false));
		// Check the tortuosity
		double tortuosity1 = vector_avg(morph.getTortuosityData((char)1));
		double tortuosity2 = vector_avg(morph.getTortuosityData((char)2));
		EXPECT_NEAR(1.1, tortuosity1, 0.025);
		EXPECT_NEAR(1.1, tortuosity2, 0.025);
		// Calculate the tortuosity using the reduced memory option
		morph.calculateTortuosity((char)1, true);
		morph.calculateTortuosity((char)2, true);
		// Check that both tortuosity methods give the same answer
		auto data1 = morph.getTortuosityData((char)1);
		auto data2 = morph.getTortuosityData((char)2);
		EXPECT_DOUBLE_EQ(tortuosity1, vector_avg(data1));
		EXPECT_DOUBLE_EQ(tortuosity2, vector_avg(data2));
		// Output tortuosity maps
		ofstream outfile("tortuosity_maps.txt");
		morph.outputTortuosityMaps(outfile);
		outfile.close();
		// Check that output looks valid
		ifstream infile("tortuosity_maps.txt");
		string line;
		// Check first column names line
		getline(infile, line);
		EXPECT_EQ("X-Position,Y-Position,Tortuosity1,Tortuosity2", line);
		// Check first data line
		getline(infile, line);
		EXPECT_TRUE(line.find("0,0,") != string::npos);
	}

	TEST_F(MorphologyTest, SmoothingAndResizeTests) {
		Morphology morph = *morph_start;
		// Calculate the interfacial volume fraction and interfacial area to volume ratio of the starting morphology
		double iv_frac_i = morph.calculateInterfacialVolumeFraction();
		double iav_ratio_i = morph.calculateInterfacialAreaVolumeRatio();
		// Get the initial domain size
		double domain_size1_i = morph.getDomainSize((char)1);
		double domain_size2_i = morph.getDomainSize((char)2);
		// Apply smoothing
		morph.executeSmoothing(0.52, 1);
		// Check that the mix fractions are not severely changed
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.02);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.02);
		// Recalculate the domain size
		morph.calculateCorrelationDistances();
		double domain_size1_f = morph.getDomainSize((char)1);
		double domain_size2_f = morph.getDomainSize((char)2);
		// Check that the domain size has not significantly increased
		EXPECT_NEAR(domain_size1_i, domain_size1_f, 1.0);
		EXPECT_NEAR(domain_size2_i, domain_size2_f, 1.0);
		// Calculate the new interfacial volume fraction and interfacial area to volume ratio
		double iv_frac_f = morph.calculateInterfacialVolumeFraction();
		double iav_ratio_f = morph.calculateInterfacialAreaVolumeRatio();
		// Check that the interfacial volume ratio and interfacial area to volume ratio have been reduced
		EXPECT_LT(iv_frac_f, iv_frac_i);
		EXPECT_LT(iav_ratio_f, iav_ratio_i);
		// Calculate the final tortuosity
		morph.calculateTortuosity((char)1, false);
		morph.calculateTortuosity((char)2, false);
		// Check the tortuosity of an isotropic smoothed phase separated blend
		auto tortuosity1_avg = vector_avg(morph.getTortuosityData((char)1));
		auto tortuosity2_avg = vector_avg(morph.getTortuosityData((char)2));
		EXPECT_NEAR(1.085, tortuosity1_avg, 0.025);
		EXPECT_NEAR(1.085, tortuosity2_avg, 0.025);
		// Check the island volume fraction
		EXPECT_LT(morph.getIslandVolumeFraction((char)1), 0.001);
		EXPECT_LT(morph.getIslandVolumeFraction((char)2), 0.001);
		// Stretch the lattice
		morph.stretchLattice(2);
		// Check the new dimensions
		EXPECT_EQ(100, morph.getLength());
		EXPECT_EQ(100, morph.getWidth());
		EXPECT_EQ(100, morph.getHeight());
		// Calculate domain size of new morphology
		morph.calculateCorrelationDistances();
		// Check the approximate magnitude of the new domain size relative to the original
		EXPECT_NEAR(2 * domain_size1_f, morph.getDomainSize((char)1), 0.5);
		EXPECT_NEAR(2 * domain_size2_f, morph.getDomainSize((char)2), 0.5);
		// Shrink the lattice
		morph.shrinkLattice(2);
		// Check the new dimensions
		EXPECT_EQ(50, morph.getLength());
		EXPECT_EQ(50, morph.getWidth());
		EXPECT_EQ(50, morph.getHeight());
		// Calculate domain size of new morphology
		morph.calculateCorrelationDistances();
		// Check the approximate magnitude of the new domain size relative to the original
		EXPECT_NEAR(domain_size1_f, morph.getDomainSize((char)1), 0.5);
		EXPECT_NEAR(domain_size2_f, morph.getDomainSize((char)2), 0.5);
		// Check behavior when trying to stretch and shrink with invalid rescale factors
		EXPECT_THROW(morph.stretchLattice(0), invalid_argument);
		EXPECT_THROW(morph.stretchLattice(0), invalid_argument);
		EXPECT_THROW(morph.shrinkLattice(0), invalid_argument);
		EXPECT_THROW(morph.shrinkLattice(0), invalid_argument);
		EXPECT_THROW(morph.shrinkLattice(3), invalid_argument);
		EXPECT_THROW(morph.shrinkLattice(3), invalid_argument);
	}

	TEST_F(MorphologyTest, InterfacialTests) {
		Morphology morph = *morph_start;
		// Calculate the interfacial distance histogram
		morph.calculateInterfacialDistanceHistogram();
		auto hist1 = morph.getInterfacialDistanceHistogram((char)1);
		auto hist2 = morph.getInterfacialDistanceHistogram((char)2);
		auto cum_hist1 = calculateCumulativeHist(calculateProbabilityHist(hist1));
		auto cum_hist2 = calculateCumulativeHist(calculateProbabilityHist(hist2));
		// Check that the probability histograms add up to 1
		EXPECT_DOUBLE_EQ(1.0, cum_hist1.back().second);
		EXPECT_DOUBLE_EQ(1.0, cum_hist2.back().second);
		// Get initial domain sizes
		double domain_size1_i = morph.getDomainSize((char)1);
		double domain_size2_i = morph.getDomainSize((char)2);
		// Apply smoothing
		morph.executeSmoothing(0.52, 1);
		// Save smoothed morphology
		Morphology morph_smoothed = morph;
		// Calculate the interfacial volume fraction and interfacial area to volume ratio of the smoothed morphology
		double iv_frac_i = morph.calculateInterfacialVolumeFraction();
		double iav_ratio_i = morph.calculateInterfacialAreaVolumeRatio();
		// Apply interfacial mixing to the smoothed morphology
		morph.executeMixing(2.0, 0.5);
		// Calculate the new interfacial volume fraction and interfacial area to volume ratio
		double iv_frac_mix = morph.calculateInterfacialVolumeFraction();
		double iav_ratio_mix = morph.calculateInterfacialAreaVolumeRatio();
		// Check that the interfacial volume ratio and interfacial area to volume ratio have been increased
		EXPECT_LT(iv_frac_i, iv_frac_mix);
		EXPECT_LT(iav_ratio_i, iav_ratio_mix);
		// Calculate domain size
		morph.calculateCorrelationDistances();
		// Check that the domain size has not greatly decreased
		EXPECT_NEAR(domain_size1_i, morph.getDomainSize((char)1), 0.5);
		EXPECT_NEAR(domain_size2_i, morph.getDomainSize((char)2), 0.5);
		// Reset morph to smoothed morphology
		morph = morph_smoothed;
		// Apply interfacial mixing with higher interfacial mixing ratio to the smoothed morphology
		morph.executeMixing(2.0, 0.6);
		// Calculate the new interfacial volume fraction and interfacial area to volume ratio
		iv_frac_mix = morph.calculateInterfacialVolumeFraction();
		iav_ratio_mix = morph.calculateInterfacialAreaVolumeRatio();
		// Check that the interfacial volume ratio and interfacial area to volume ratio have been increased
		EXPECT_LT(iv_frac_i, iv_frac_mix);
		EXPECT_LT(iav_ratio_i, iav_ratio_mix);
		// Calculate domain size
		morph.calculateCorrelationDistances();
		// Check that the domain size has not greatly decreased
		EXPECT_NEAR(domain_size1_i, morph.getDomainSize((char)1), 0.5);
		EXPECT_NEAR(domain_size2_i, morph.getDomainSize((char)2), 0.5);
	}

	TEST_F(MorphologyTest, ExportImportTests) {
		//// Create a local copy of the Morphology object
		Morphology morph = *morph_start;
		// Get the mix fractions of the original morphology
		double mix_fraction1 = morph.getMixFraction((char)1);
		double mix_fraction2 = morph.getMixFraction((char)2);
		// Get the domain size of the initial morphology
		double domain_size1 = morph.getDomainSize((char)1);
		double domain_size2 = morph.getDomainSize((char)2);
		// output original morphology in compressed format
		ofstream outfile1("morphology_file1.txt");
		morph.outputMorphologyFile("v4.0", outfile1, true);
		outfile1.close();
		// Import compressed morphology file
		cout << "Importing compressed morphology file." << endl;
		ifstream infile1("morphology_file1.txt");
		EXPECT_TRUE(morph.importMorphologyFile(infile1));
		infile1.close();
		// Check dimensions of imported morphology
		EXPECT_EQ(50, morph.getLength());
		EXPECT_EQ(50, morph.getWidth());
		EXPECT_EQ(50, morph.getHeight());
		// Check that the mix fractions remained the same
		EXPECT_NEAR(mix_fraction1, morph.getMixFraction((char)1), 0.001);
		EXPECT_NEAR(mix_fraction2, morph.getMixFraction((char)2), 0.001);
		// Calculate domain size of imported morphology
		morph.calculateCorrelationDistances();
		// Check that domain size of original and imported morphologies are the same
		EXPECT_NEAR(domain_size1, morph.getDomainSize((char)1), 0.05);
		EXPECT_NEAR(domain_size2, morph.getDomainSize((char)2), 0.05);
		// output original morphology in uncompressed format
		ofstream outfile2("morphology_file2.txt");
		morph.outputMorphologyFile("v4.0", outfile2, false);
		outfile2.close();
		// Import uncompressed morphology file
		cout << "Importing uncompressed morphology file." << endl;
		ifstream infile2("morphology_file2.txt", ifstream::in);
		EXPECT_TRUE(morph.importMorphologyFile(infile2));
		infile2.close();
		// Check dimensions of imported morphology
		EXPECT_EQ(50, morph.getLength());
		EXPECT_EQ(50, morph.getWidth());
		EXPECT_EQ(50, morph.getHeight());
		// Check that the mix fractions remained the same
		EXPECT_NEAR(mix_fraction1, morph.getMixFraction((char)1), 0.001);
		EXPECT_NEAR(mix_fraction2, morph.getMixFraction((char)2), 0.001);
		// Calculate domain size of imported morphology
		morph.calculateCorrelationDistances();
		// Check that domain size of original and imported morphologies are the same
		EXPECT_NEAR(domain_size1, morph.getDomainSize((char)1), 0.05);
		EXPECT_NEAR(domain_size2, morph.getDomainSize((char)2), 0.05);
		// Test attempt to load a corrupted compressed morphology file with missing data
		// Load compressed file into string vector
		string line;
		vector<string> file_data;
		ifstream infile3("morphology_file1.txt");
		while (getline(infile3, line)) {
			file_data.push_back(line);
		}
		infile3.close();
		// Delete last 2 lines
		file_data.pop_back();
		file_data.pop_back();
		// Save data back to file
		ofstream outfile3("morphology_file1.txt");
		for (auto& item : file_data) {
			outfile3 << item << endl;
		}
		outfile3.close();
		// Try importing corrupted compressed file
		ifstream infile4("morphology_file1.txt");
		EXPECT_FALSE(morph.importMorphologyFile(infile4));
		infile4.close();
		// Test attempt to load a corrupted uncompressed morphology file with missing data
		// Load uncompressed file into string vector
		file_data.clear();
		ifstream infile5("morphology_file2.txt");
		while (getline(infile5, line)) {
			file_data.push_back(line);
		}
		infile5.close();
		// Delete last 2 lines
		file_data.pop_back();
		file_data.pop_back();
		// Save data back to file
		ofstream outfile4("morphology_file2.txt");
		for (auto& item : file_data) {
			outfile4 << item << endl;
		}
		outfile4.close();
		// Try importing corrupted compressed file
		ifstream infile6("morphology_file2.txt");
		EXPECT_FALSE(morph.importMorphologyFile(infile6));
		infile6.close();
	}

	TEST_F(MorphologyTest, OtherOutputTests) {
		Morphology morph = *morph_start;
		string line;
		// Output composition maps
		ofstream outfile1("composition_maps.txt");
		morph.outputCompositionMaps(outfile1);
		outfile1.close();
		// Check that output looks valid
		ifstream infile1("composition_maps.txt");
		// Check first column names line
		getline(infile1, line);
		EXPECT_EQ("X-Position,Y-Position,Composition1,Composition2", line);
		// Check first data line
		getline(infile1, line);
		EXPECT_TRUE(line.find("0,0,") != string::npos);
		infile1.close();
		// Output morphology cross-section
		ofstream outfile2("morphology_cross.txt");
		morph.outputMorphologyCrossSection(outfile2);
		outfile2.close();
		// Check that output looks valid
		ifstream infile2("morphology_cross.txt");
		// Check first column names line
		getline(infile2, line);
		EXPECT_EQ("X-Position,Y-Position,Z-Position,Site_type", line);
		// Check first data line
		getline(infile2, line);
		EXPECT_TRUE(line.find(to_string(morph.getLength() / 2) + ",0,0,") != string::npos);
		infile2.close();

	}
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	// Redirect cout to NULL to suppress command line output during the tests
	//cout.rdbuf(NULL);
	return RUN_ALL_TESTS();
}
