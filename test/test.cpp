// Copyright (c) 2014-2018 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "gtest/gtest.h"
#include "Morphology.h"
#include "Lattice.h"
#include "Utils.h"

using namespace std;
using namespace Utils;

namespace MorphologyTests {

	TEST(MorphologyTests, ConstructorTests) {
		// Test the default constructor
		Morphology morph(0);
		// Check that ID number is set properly
		EXPECT_EQ(0, morph.getID());
		// Check that the lattice has default size of 0,0,0.
		EXPECT_EQ(0, morph.getLength());

		// Test the standard constructor
		morph = Morphology(50, 60, 70, false, 1);
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
		mt19937_64 gen;
		Lattice lattice;
		// Initialize Lattice object
		lattice.init(params_lattice, &gen);
		morph = Morphology(lattice, 2);
		// Check that ID number is set properly
		EXPECT_EQ(2, morph.getID());
		// Check that the lattice has the correct size.
		EXPECT_EQ(20, morph.getLength());
		EXPECT_EQ(30, morph.getWidth());
		EXPECT_EQ(40, morph.getHeight());
		// Check that the lattice has the correct unit size
		EXPECT_DOUBLE_EQ(1.5, morph.getUnitSize());
	}

	TEST(MorphologyTests, IsotropicPhaseSeparationTests) {
		Morphology morph(50, 50, 50, true, 0);
		vector<double> mix_fractions;
		mix_fractions.assign(2, 0.5);
		morph.createRandomMorphology(mix_fractions);
		// Check that the mix fractions were implemented properly
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.01);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.01);
		// Calculate initial domain size
		CorrelationCalc_Params params;
		params.N_sampling_max = 100000;
		params.Enable_e_method = true;
		params.Enable_mix_frac_method = false;
		params.Enable_extended_correlation_calc = false;
		params.Correlation_cutoff_distance = 3;
		// Calculate the inital domain size values
		morph.calculateCorrelationDistances(params);
		double domain_size1_i = morph.getDomainSize((char)1);
		double domain_size2_i = morph.getDomainSize((char)2);
		// Check that random blend domain size is less than 2
		EXPECT_LT(domain_size1_i, 2.0);
		EXPECT_LT(domain_size2_i, 2.0);
		// Calculate the initial domain anisotropy
		morph.calculateAnisotropies(params.N_sampling_max);
		// Check that the random blend is isotropic
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)1), 0.05);
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)2), 0.05);
		// Perform some phase separation
		morph.executeIsingSwapping(200, 0.4, 0.4, false, 0, 0.0);
		// Calculate the final domain size values
		morph.calculateCorrelationDistances(params);
		double domain_size1_f = morph.getDomainSize((char)1);
		double domain_size2_f = morph.getDomainSize((char)2);
		// Check that the domain size has increased
		EXPECT_LT(domain_size1_i, domain_size1_f);
		EXPECT_LT(domain_size2_i, domain_size2_f);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(5.0, domain_size1_f, 0.5);
		EXPECT_NEAR(5.0, domain_size2_f, 0.5);
		// Check new domain anisotropy
		morph.calculateAnisotropies(params.N_sampling_max);
		// Check that the phase separated blend is still isotropic
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)1), 0.1);
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)2), 0.1);
		// Check the tortuosity of an isotropic phase separated blend
		morph.calculateTortuosity((char)1, false);
		morph.calculateTortuosity((char)2, false);
		auto data = morph.getTortuosityData((char)1);
		EXPECT_NEAR(1.1, vector_avg(data), 0.02);
		data = morph.getTortuosityData((char)2);
		EXPECT_NEAR(1.1, vector_avg(data), 0.02);
		// Apply smoothing and re-recalculate the domain size
		morph.executeSmoothing(0.52, 1);
		morph.calculateCorrelationDistances(params);
		domain_size1_f = morph.getDomainSize((char)1);
		domain_size2_f = morph.getDomainSize((char)2);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(6.0, domain_size1_f, 0.5);
		EXPECT_NEAR(6.0, domain_size2_f, 0.5);
		// Check that the mix fractions are not severely changed
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.02);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.02);
		// Calculate domain size using the mix fraction method
		params.Enable_e_method = false;
		params.Enable_mix_frac_method = true;
		morph.calculateCorrelationDistances(params);
		domain_size1_f = morph.getDomainSize((char)1);
		domain_size2_f = morph.getDomainSize((char)2);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(6.0, domain_size1_f, 0.5);
		EXPECT_NEAR(6.0, domain_size2_f, 0.5);
		// Calculate the interfacial distance histogram
		morph.calculateInterfacialDistanceHistogram();
		auto prob1 = morph.getInterfacialDistanceHistogram((char)1);
		auto prob2 = morph.getInterfacialDistanceHistogram((char)2);
		// Check that the probability histograms sum to 1
		EXPECT_NEAR(1.0, accumulate(prob1.begin(), prob1.end(), 0.0), 0.05);
		EXPECT_NEAR(1.0, accumulate(prob2.begin(), prob2.end(), 0.0), 0.05);
		// Check that the size of the histograms are related to the domain size
		EXPECT_NEAR(5.0, prob1.size(), 1);
		EXPECT_NEAR(5.0, prob2.size(), 1);
	}

	TEST(MorphologyTests, AnisotropicPhaseSeparationTests) {
		Morphology morph(40, 40, 40, true, 0);
		vector<double> mix_fractions;
		mix_fractions.assign(2, 0.5);
		morph.createRandomMorphology(mix_fractions);
		// Check that the mix fractions were implemented properly
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.01);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.01);
		// Calculate initial domain size
		CorrelationCalc_Params params;
		params.N_sampling_max = 100000;
		params.Enable_e_method = true;
		params.Enable_mix_frac_method = false;
		params.Enable_extended_correlation_calc = false;
		params.Correlation_cutoff_distance = 3;
		// Calculate the inital domain size values
		morph.calculateCorrelationDistances(params);
		double domain_size1_i = morph.getDomainSize((char)1);
		double domain_size2_i = morph.getDomainSize((char)2);
		// Check that random blend domain size is less than 2
		EXPECT_LT(domain_size1_i, 2.0);
		EXPECT_LT(domain_size2_i, 2.0);
		// Calculate the initial domain anisotropy
		morph.calculateAnisotropies(params.N_sampling_max);
		// Check that the initial random blend is isotropic
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)1), 0.05);
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)2), 0.05);
		// Perform some anisotropic phase separation that creates aligned structures out-of-plane in the z-direction
		morph.executeIsingSwapping(220, 0.35, 0.35, true, 3, 0.05);
		// Calculate the final domain size values
		morph.calculateCorrelationDistances(params);
		double domain_size1_f = morph.getDomainSize((char)1);
		double domain_size2_f = morph.getDomainSize((char)2);
		// Check that the domain size has increased
		EXPECT_LT(domain_size1_i, domain_size1_f);
		EXPECT_LT(domain_size2_i, domain_size2_f);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(5.0, domain_size1_f, 0.5);
		EXPECT_NEAR(5.0, domain_size2_f, 0.5);
		// Calculate the final domain anisotropy
		morph.calculateAnisotropies(params.N_sampling_max);
		// Check the approximate magnitude of the anisotropy factor
		EXPECT_NEAR(1.4, morph.getDomainAnisotropy((char)1), 0.15);
		EXPECT_NEAR(1.4, morph.getDomainAnisotropy((char)2), 0.15);
		// Reset morphology to a random blend
		morph.createRandomMorphology(mix_fractions);
		// Calculate the inital domain size values
		morph.calculateCorrelationDistances(params);
		domain_size1_i = morph.getDomainSize((char)1);
		domain_size2_i = morph.getDomainSize((char)2);
		// Check that random blend domain size is less than 2
		EXPECT_LT(domain_size1_i, 2.0);
		EXPECT_LT(domain_size2_i, 2.0);
		// Perform some anisotropic phase separation that creates aligned structures in the x-y plane
		morph.executeIsingSwapping(180, 0.4, 0.4, true, 3, -0.05);
		// Calculate the final domain size values
		morph.calculateCorrelationDistances(params);
		domain_size1_f = morph.getDomainSize((char)1);
		domain_size2_f = morph.getDomainSize((char)2);
		// Check that the domain size has increased
		EXPECT_LT(domain_size1_i, domain_size1_f);
		EXPECT_LT(domain_size2_i, domain_size2_f);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(5.0, domain_size1_f, 0.5);
		EXPECT_NEAR(5.0, domain_size2_f, 0.5);
		// Calculate the final domain anisotropy
		morph.calculateAnisotropies(params.N_sampling_max);
		// Check the approximate magnitude of the anisotropy factor
		EXPECT_NEAR(0.80, morph.getDomainAnisotropy((char)1), 0.1);
		EXPECT_NEAR(0.80, morph.getDomainAnisotropy((char)2), 0.1);

	}
}

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

	TEST(UtilsTests, CalculateProbabilityHistTests) {
		mt19937_64 gen(std::random_device{}());
		uniform_real_distribution<> dist(0, 100);
		vector<double> data((int)1e7);
		for (int i = 0; i < (int)data.size(); i++) {
			data[i] = dist(gen);
		}
		auto hist = calculateProbabilityHist(data, 10);
		uniform_int_distribution<> dist2(0, 9);
		EXPECT_EQ(10, (int)hist.size());
		EXPECT_NEAR(1.0 / 100.0, hist[dist2(gen)].second, 1e-4);
		EXPECT_NEAR(1.0 / 100.0, hist[dist2(gen)].second, 1e-4);
		EXPECT_NEAR(1.0 / 100.0, hist[dist2(gen)].second, 1e-4);
		hist = calculateProbabilityHist(data, 10.0);
		EXPECT_EQ(10, (int)hist.size());
		EXPECT_NEAR(1.0 / 100.0, hist[dist2(gen)].second, 1e-4);
		EXPECT_NEAR(1.0 / 100.0, hist[dist2(gen)].second, 1e-4);
		EXPECT_NEAR(1.0 / 100.0, hist[dist2(gen)].second, 1e-4);
		data.clear();
		hist = calculateProbabilityHist(data, 10.0);
		EXPECT_DOUBLE_EQ(0.0, hist[0].first);
		EXPECT_DOUBLE_EQ(0.0, hist[0].second);
		hist = calculateProbabilityHist(data, 5);
		EXPECT_DOUBLE_EQ(0.0, hist[0].first);
		EXPECT_DOUBLE_EQ(0.0, hist[0].second);
		hist = calculateProbabilityHist(data, 1.0, 5);
		EXPECT_DOUBLE_EQ(0.0, hist[0].first);
		EXPECT_DOUBLE_EQ(0.0, hist[0].second);
		data = { 0.0, 1.0, 2.0, 3.0, 4.0 };
		hist = calculateProbabilityHist(data, 10);
		EXPECT_EQ(5, (int)hist.size());
		hist = calculateProbabilityHist(data, 0.1);
		EXPECT_EQ(5, (int)hist.size());
	}

	TEST(UtilsTests, ImportBooleanTests) {
		bool error_status;
		EXPECT_TRUE(importBooleanParam("true", error_status));
		EXPECT_TRUE(importBooleanParam(" true  ", error_status));
		EXPECT_FALSE(importBooleanParam("false	", error_status));
		EXPECT_FALSE(importBooleanParam("   false", error_status));
		EXPECT_FALSE(importBooleanParam("blah", error_status));
		EXPECT_FALSE(importBooleanParam("	blah  ", error_status));
		EXPECT_TRUE(error_status);
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
		// negative ints
		for (int i = 0; i < 10; i++) {
			int_data[i] = -(i + 1);
		}
		EXPECT_DOUBLE_EQ(-5.5, vector_avg(int_data));
		EXPECT_NEAR(3.02765035409749, vector_stdev(int_data), 1e-14);
		// positive doubles
		vector<double> double_data;
		double_data.assign(10, 0);
		for (int i = 0; i < 10; i++) {
			double_data[i] = 0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(2.75, vector_avg(double_data));
		EXPECT_NEAR(1.51382517704875, vector_stdev(double_data), 1e-14);
		// negative doubles
		double_data.assign(10, 0);
		for (int i = 0; i < 10; i++) {
			double_data[i] = -0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(-2.75, vector_avg(double_data));
		EXPECT_NEAR(1.51382517704875, vector_stdev(double_data), 1e-14);
	}
}

namespace LatticeTests {

	class LatticeTest : public ::testing::Test {
	protected:
		mt19937_64 gen;
		Lattice_Params params_lattice;
		Lattice lattice;

		void SetUp() {
			gen.seed(std::random_device{}());
			// Setup params
			params_lattice.Enable_periodic_x = true;
			params_lattice.Enable_periodic_y = true;
			params_lattice.Enable_periodic_z = true;
			params_lattice.Length = 50;
			params_lattice.Width = 50;
			params_lattice.Height = 50;
			params_lattice.Unit_size = 1.0;
			// Initialize Lattice object
			lattice.init(params_lattice, &gen);
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
		lattice2.init(params_lattice, &gen);
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
		lattice2.init(params_lattice, &gen);
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

	TEST_F(LatticeTest, GetSiteCoordsTests) {
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
	}

}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	// Redirect cout to NULL to suppress command line output during the tests
	//cout.rdbuf(NULL);
	return RUN_ALL_TESTS();
}
