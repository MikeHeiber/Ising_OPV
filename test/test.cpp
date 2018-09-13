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
		Lattice lattice;
		// Initialize Lattice object
		lattice.init(params_lattice);
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

	TEST(MorphologyTests, RandomMorphologyTests) {
		// Initialize Morphology object
		Morphology morph(50, 50, 50, true, 0);
		vector<double> mix_fractions;
		mix_fractions.assign(2, 0.5);
		morph.createRandomMorphology(mix_fractions);
		// Check that the mix fractions were implemented properly
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.001);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.001);
		// Calculate the inital domain size values
		CorrelationCalc_Params params;
		params.N_sampling_max = 100000;
		params.Enable_e_method = true;
		params.Enable_mix_frac_method = false;
		params.Enable_extended_correlation_calc = false;
		params.Extended_correlation_cutoff_distance = 3;
		morph.calculateCorrelationDistances(params);
		double domain_size1 = morph.getDomainSize((char)1);
		double domain_size2 = morph.getDomainSize((char)2);
		// Check that random blend domain size is less than 2
		EXPECT_LT(domain_size1, 2.0);
		EXPECT_LT(domain_size2, 2.0);
		// Calculate the initial domain anisotropy
		morph.calculateAnisotropies(params.N_sampling_max);
		// Check that the random blend is isotropic
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)1), 0.05);
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)2), 0.05);
	}

	TEST(MorphologyTests, AnisotropicPhaseSeparationTests) {
		Morphology morph(40, 40, 40, false, 0);
		vector<double> mix_fractions;
		mix_fractions.assign(2, 0.5);
		morph.createRandomMorphology(mix_fractions);
		// Check that the mix fractions were implemented properly
		EXPECT_NEAR(0.5, morph.getMixFraction((char)1), 0.001);
		EXPECT_NEAR(0.5, morph.getMixFraction((char)2), 0.001);
		// Calculate initial domain size
		CorrelationCalc_Params params;
		params.N_sampling_max = 100000;
		params.Enable_e_method = true;
		params.Enable_mix_frac_method = false;
		params.Enable_extended_correlation_calc = false;
		params.Extended_correlation_cutoff_distance = 3;
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
		EXPECT_NEAR(1.4, morph.getDomainAnisotropy((char)1), 0.25);
		EXPECT_NEAR(1.4, morph.getDomainAnisotropy((char)2), 0.25);
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
		EXPECT_NEAR(0.83, morph.getDomainAnisotropy((char)1), 0.1);
		EXPECT_NEAR(0.83, morph.getDomainAnisotropy((char)2), 0.1);
		// Reset morphology to a random blend
		morph.createRandomMorphology(mix_fractions);
		// Perform some anisotropic phase separation that creates aligned structures in the x-direction
		morph.executeIsingSwapping(200, 0.35, 0.35, true, 1, 0.05);
		// Calculate the domain anisotropy
		morph.calculateAnisotropies(params.N_sampling_max);
		auto anisotropy1 = morph.getDomainAnisotropy((char)1);
		auto anisotropy2 = morph.getDomainAnisotropy((char)2);
		// Reset morphology to a random blend
		morph.createRandomMorphology(mix_fractions);
		// Perform some anisotropic phase separation that creates aligned structures in the y-direction
		morph.executeIsingSwapping(200, 0.35, 0.35, true, 2, 0.05);
		// Calculate the domain anisotropy
		morph.calculateAnisotropies(params.N_sampling_max);
		// Check that anisotropic phase separation in the x- or y-directions produces approximately the same anisotropy metric
		EXPECT_NEAR(anisotropy1, morph.getDomainAnisotropy((char)1), 0.1);
		EXPECT_NEAR(anisotropy2, morph.getDomainAnisotropy((char)2), 0.1);
	}

	TEST(MorphologyTests, MorphologyAnalysisTests) {
		Morphology morph(40, 40, 40, false, 0);
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
		CorrelationCalc_Params params;

		static void SetUpTestCase() {
			// Correlation function calculation parameters	
			CorrelationCalc_Params params1;
			params1.N_sampling_max = 100000;
			params1.Enable_e_method = true;
			params1.Enable_mix_frac_method = false;
			params1.Enable_extended_correlation_calc = false;
			params1.Extended_correlation_cutoff_distance = 3;
			// Initialize Morphology object
			morph_start = new Morphology(50, 50, 50, true, 0);
			vector<double> mix_fractions;
			mix_fractions.assign(2, 0.5);
			morph_start->createRandomMorphology(mix_fractions);
			// Perform some phase separation
			morph_start->executeIsingSwapping(290, 0.4, 0.4, false, 0, 0.0);
			// Calculate the domain size values
			morph_start->calculateCorrelationDistances(params1);
		}

		static void TearDownTestCase() {
			delete morph_start;
			morph_start = NULL;
		}

		virtual void SetUp() {
			// Correlation function calculation parameters	
			params.N_sampling_max = 100000;
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

	TEST_F(MorphologyTest, IsotropicAnalysisTests) {
		Morphology morph = *morph_start;
		double domain_size1_i = morph.getDomainSize((char)1);
		double domain_size2_i = morph.getDomainSize((char)2);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(6.0, domain_size1_i, 0.5);
		EXPECT_NEAR(6.0, domain_size2_i, 0.5);
		// Calculate domain size using the mix fraction method
		params.Enable_e_method = false;
		params.Enable_mix_frac_method = true;
		morph.calculateCorrelationDistances(params);
		double domain_size1_f = morph.getDomainSize((char)1);
		double domain_size2_f = morph.getDomainSize((char)2);
		// Check the approximate magnitude of the domain size
		EXPECT_NEAR(6.0, domain_size1_f, 0.5);
		EXPECT_NEAR(6.0, domain_size2_f, 0.5);
		// Check domain anisotropy
		params.Enable_e_method = true;
		params.Enable_mix_frac_method = false;
		morph.calculateAnisotropies(params.N_sampling_max);
		// Check that the phase separated blend is isotropic
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)1), 0.1);
		EXPECT_NEAR(1.0, morph.getDomainAnisotropy((char)2), 0.1);
		// Calculate the tortuosity
		morph.calculateTortuosity((char)1, false);
		morph.calculateTortuosity((char)2, false);
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
		// Calculate the depth dependent characteristics
		morph.calculateDepthDependentData(params);
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
		// Calculate the interfacial volume fraction and interfacial area to volume ratio
		double iv_frac_i = morph.calculateInterfacialVolumeFraction();
		double iav_ratio_i = morph.calculateInterfacialAreaVolumeRatio();
		// Check depth dependent interfacial volume fraction compared to bulk value
		auto data = morph.getDepthIVData();
		for (auto item : data) {
			EXPECT_NEAR(iv_frac_i, item, 0.075);
		}
		// Calculate the interfacial distance histogram
		morph.calculateInterfacialDistanceHistogram();
		auto prob1 = morph.getInterfacialDistanceHistogram((char)1);
		auto prob2 = morph.getInterfacialDistanceHistogram((char)2);
		// Check that the probability histograms sum to 1
		EXPECT_NEAR(1.0, accumulate(prob1.begin(), prob1.end(), 0.0), 0.05);
		EXPECT_NEAR(1.0, accumulate(prob2.begin(), prob2.end(), 0.0), 0.05);
		// Check that the size of the histograms are related to the domain size
		EXPECT_NEAR(6.0 / 2.0, prob1.size(), 1);
		EXPECT_NEAR(6.0 / 2.0, prob2.size(), 1);
	}

	TEST_F(MorphologyTest, SmoothingTests) {
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
		morph.calculateCorrelationDistances(params);
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
	}

	TEST_F(MorphologyTest, InterfacialMixingTests) {
		Morphology morph = *morph_start;
		// Get initial domain sizes
		double domain_size1_i = morph.getDomainSize((char)1);
		double domain_size2_i = morph.getDomainSize((char)2);
		// Apply smoothing
		morph.executeSmoothing(0.52, 1);
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
		morph.calculateCorrelationDistances(params);
		// Check that the domain size has not greatly decreased
		EXPECT_NEAR(domain_size1_i, morph.getDomainSize((char)1), 0.5);
		EXPECT_NEAR(domain_size2_i, morph.getDomainSize((char)2), 0.5);
	}

	TEST_F(MorphologyTest, ExportImportTests) {
		// Get the domain size of the initial morphology
		double domain_size1 = morph_start->getDomainSize((char)1);
		double domain_size2 = morph_start->getDomainSize((char)2);
		// output original morphology
		ofstream outfile("morphology_file.txt", ofstream::out | ofstream::trunc);
		morph_start->outputMorphologyFile("v4.0", outfile, true);
		outfile.close();
		// Create a local copy of the Morphology object
		Morphology morph = *morph_start;
		ifstream infile("morphology_file.txt", ifstream::in);
		morph.importMorphologyFile(infile);
		infile.close();
		// Calculate domain size of imported morphology
		morph.calculateCorrelationDistances(params);
		// Check that domain size of original and imported morphologies are the same
		EXPECT_NEAR(domain_size1, morph.getDomainSize((char)1), 0.001);
		EXPECT_NEAR(domain_size2, morph.getDomainSize((char)2), 0.001);
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
