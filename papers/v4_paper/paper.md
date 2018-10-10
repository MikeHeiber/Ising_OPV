---
title: 'Ising_OPV v4.0: Experimental Tomography Data Import, Interpretation, and Analysis'
tags:
  - organic photovoltaics
  - bulk heterojunction morphology
  - tomography
  - Ising model
  - phase separation
  - image analysis
authors:
 - name: Michael C. Heiber
   orcid: 0000-0002-1567-5663
   affiliation: "1"
affiliations:
 - name: Center for Hierarchical Materials Design (CHiMaD), Northwestern University, Evanston, Illinois 60208, USA
   index: 1
date: 6 October 2018
bibliography: paper.bib
---

# Summary

Understanding the impact of the complex meso-scale morphology is critical for the development of organic semiconductor materials and devices.  This is particularly important in organic photovoltaics (OPVs), where a blend of two or more components phase separates to form a bulk heterojunction (BHJ) structure.  To build better structure-property models for organic BHJ photovoltaics, the simple Ising-based morphology model has proven to be a highly useful tool when coupled with kinetic Monte Carlo (KMC) simulations.[@heiber2019chapter] ``Ising_OPV`` was originally designed as an efficient, open-source C++ tool that would enable researchers in the community to have easy access to this morphology model and allow them to create well-controlled morphologies on an HPC cluster for KMC simulations.[@heiber2014prapp] Demonstrating the utility of this tool, the ability to systematically control the domain size allowed a detailed investigation of the charge carrier recombination kinetics in OPVs.[@heiber2015prl, heiber2016prb] The tool can also create controlled interfacial mixing, which can be important for simulating the exciton dissociation dynamics and charge separation yield in OPVs.[@lyons2012ees, @heiber2013jpcc] In addition, the tool was later updated to add new features that allow further structural control and quantification of important morphological features, most importantly the domain tortuosity.[@heiber2017prapp] The tool has also been used as a testbed for developing more advanced 3D image analysis methods.[@aboulhassan2017cgf]

Building on this foundation, v4.0 adds an exciting new feature that allows users to import three-dimensional morphology data sets from experimental techniques such as electron tomography [@vanbavel2009nl, @pfannmoeller2013ees] or atom probe tomography [@proudian2018arXiv] and prepare experimentally-derived morphology sets for KMC simulations using ``Excimontec``.[@heiber2018excimontec] A pictorial representation of the workflow when importing experimental morphology data is shown below. In addition, this update includes a major code overhaul to create a well-organized and well-documented object-oriented software package that is more reliable, testable, and extensible. The code has been updated to use many C++11 features and now includes rigorous unit testing with ``googletest``, integration testing with ``TravisCI``, and API documentation generated using ``Doxygen``.  The source code for ``Ising_OPV v4.0`` is archived with Zenodo.[@heiber2018ising4]

-![Tomogram analysis workflow.](tomogram_analysis.png)

# Acknowledgments

This work was developed under the financial assistance award 70NANB14H012 from U.S. Department of Commerce, National Institute of Standards and Technology as part of the Center for Hierarchical Materials Design (CHiMaD).  Thank you to Dr. Dean M. DeLongchamp for providing access to NIST's Raritan computing cluster and Dr. Andrew A. Herzing for providing TEM tomography data, which was helpful with software development and testing.

# References
