// Copyright (c) 2014-2019 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#ifndef LATTICE_H
#define LATTICE_H

#include "Utils.h"
#include <ctime>
#include <functional>
#include <stdexcept>

namespace Ising_OPV {

	//! \brief This class contains the properties of a three-dimensional lattice and the functions needed to interact with it.
	//! \details The class makes use of the Lattice_Params struct to load the neccessary input parameters and the Coords struct
	//! to record the Cartesian coordinates of each lattice site.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2014-2019
	class Lattice {

		struct Site {
			char type = 0;
		};

	public:

		//! \brief This struct contains all of the main input parameters needed by the Lattice class.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2014-2019
		struct Lattice_Params {
			//! Determines whether the x-direction periodic boundaries will be enabled. 
			bool Enable_periodic_x = false;
			//! Determines whether the y-direction periodic boundaries will be enabled. 
			bool Enable_periodic_y = false;
			//! Determines whether the z-direction periodic boundaries will be enabled. 
			bool Enable_periodic_z = false;
			//! Defines the desired x-direction size of the lattice.
			int Length = 0;
			//! Defines the desired y-direction size of the lattice.
			int Width = 0;
			//! Defines the desired z-direction size of the lattice.
			int Height = 0;
			//! Defines the desired lattice unit size, which is used to convert lattice units into real space units.
			double Unit_size = 0.0; // nm
		};

		//! \brief Default constructor that creates an empty Lattice object.
		//! \warning An empty lattice object should not be used until initialized using the init function.
		Lattice();

		//! \brief Initializes the Lattice object using the provided Parameters_Lattice input parameter struct.
		//! \param params is a Lattice_Params struct that contains all of the required
		//! parameters to initialize the Lattice object. 
		void init(const Lattice_Params& params);

		//! \brief Calculates the destination coordinates when given the starting coordinates and the displacement vector (i,j,k).
		//! \details When the starting coordinates are near one or more of the lattice boundaries and periodic boundary conditions are enabled,
		//! the function detemines the destination coordinates across the periodic boundary and assigns the calculated Coords struct to the input
		//! coords_dest argument.
		//! \param coords_initial is the Coords struct tht designates the starting coordinates.
		//! \param i is the displacement in the x-direction.
		//! \param j is the displacement in the y-direction.
		//! \param k is the displacement in the z-direction.
		//! \param coords_dest is Coords struct that indicates the output destination coordinates.
		void calculateDestinationCoords(const Coords& coords_initial, const int i, const int j, const int k, Coords& coords_dest) const;

		//! \brief Calculates a coordinate adjustment factor if the x-direction periodic boundary is crossed.
		//! \param x is the starting x coordinate.
		//! \param i is the displacement in the x-direction.
		//! \return Length (the x-direction size of the lattice) if the x periodic boundary is crossed in the negative direction.
		//! \return -Length if the x periodic boundary is crossed in the positive direction.
		//! \return 0 if the x periodic boundary is not enabled or if the x periodic boundary is not crossed.
		int calculateDX(const int x, const int i) const;

		//! \brief Calculates a coordinate adjustment factor if the x-direction periodic boundary is crossed.
		//! \param coords_initial is the Coords struct that represents the starting coordinates.
		//! \param coords_dest is the Coords struct that represents the destination coordinates.
		//! \return Length (the x-direction size of the lattice) if the x periodic boundary is crossed in the negative direction.
		//! \return -Length if the x periodic boundary is crossed in the positive direction.
		//! \return 0 if the x periodic boundary is not enabled or if the x periodic boundary is not crossed.
		int calculateDX(const Coords& coords_initial, const Coords& coords_dest) const;

		//! \brief Calculates a coordinate adjustment factor if the y-direction periodic boundary is crossed.
		//! \param y is the starting y coordinate.
		//! \param j is the displacement in the y-direction.
		//! \return Width (the y-direction size of the lattice) if the y periodic boundary is crossed in the negative direction.
		//! \return -Width if the y periodic boundary is crossed in the positive direction.
		//! \return 0 if the y periodic boundary is not enabled or if the y periodic boundary is not crossed.
		int calculateDY(const int y, const int j) const;

		//! \brief Calculates a coordinate adjustment factor if the y-direction periodic boundary is crossed.
		//! \param coords_initial is the Coords struct that represents the starting coordinates.
		//! \param coords_dest is the Coords struct that represents the destination coordinates.
		//! \return Width (the y-direction size of the lattice) if the y periodic boundary is crossed in the negative direction.
		//! \return -Width if the y periodic boundary is crossed in the positive direction.
		//! \return 0 if the y periodic boundary is not enabled or if the y periodic boundary is not crossed.
		int calculateDY(const Coords& coords_initial, const Coords& coords_dest) const;

		//! \brief Calculates a coordinate adjustment factor if the z-direction periodic boundary is crossed.
		//! \param z is the starting z coordinate.
		//! \param k is the displacement in the z-direction.
		//! \return Height (the z-direction size of the lattice) if the z periodic boundary is crossed in the negative direction.
		//! \return -Height if the z periodic boundary is crossed in the positive direction.
		//! \return 0 if the z periodic boundary is not enabled or if the z periodic boundary is not crossed.
		int calculateDZ(const int z, const int k) const;

		//! \brief Calculates a coordinate adjustment factor if the z-direction periodic boundary is crossed.
		//! \param coords_initial is the Coords struct that represents the starting coordinates.
		//! \param coords_dest is the Coords struct that represents the destination coordinates.
		//! \return Height (the z-direction size of the lattice) if the z periodic boundary is crossed in the negative direction.
		//! \return -Height if the z periodic boundary is crossed in the positive direction.
		//! \return 0 if the z periodic boundary is not enabled or if the z periodic boundary is not crossed.
		int calculateDZ(const Coords& coords_initial, const Coords& coords_dest) const;

		//! \brief Calculates the shortest distance between a pair of coordinates in squared lattice units.
		//! \param coords_start is the Coords struct that represents the starting coordinates.
		//! \param coords_dest is the Coords struct that represents the destination coordinates.
		//! \return The distance between the two sets of coordinates in squared lattice units.
		int calculateLatticeDistanceSquared(const Coords& coords_start, const Coords& coords_dest) const;

		//! \brief Checks to see if a generic move operation from the designated initial coordinates to a destination
		//! position specified by the displacement vector (i,j,k) is possible.
		//! \details The main use of this function is used to check if a proposed move event crosses a non-periodic boundary.
		//! \param coords_initial is the Coords struct that represents the starting coordinates.
		//! \param i is the displacement in the x-direction.
		//! \param j is the displacement in the y-direction.
		//! \param k is the displacement in the z-direction.
		//! \return true if a move event is possible.
		//! \return false if a move event is not possible.
		bool checkMoveValidity(const Coords& coords_initial, const int i, const int j, const int k) const;

		//! \brief Extracts a sub-lattice from a larger lattice.
		//! \param x is the starting x-coordinate of sub-lattice.
		//! \param sublength is the x-dimension of the sub-lattice.
		//! \param y is the starting y-coordinate of sub-lattice.
		//! \param subwidth is the y-dimension of the sub-lattice.
		//! \param z is the starting z-coordinate of sub-lattice.
		//! \param subheight is the z-dimension of the sub-lattice.
		//! \return A new Lattice object that contains the data from the specified sub-lattice.
		Lattice extractSublattice(const int x, const int sublength, const int y, const int subwidth, const int z, const int subheight) const;

		//! \brief Generates the coordinates for a randomly selected site in the lattice.
		//! \return A Coords struct containing the coordinates of a randomly selected site from the lattice.
		Coords generateRandomCoords();

		//! \brief Generates a random x coordinate that lies within the x-dimension size of the lattice.
		//! \return
		//! A randomly selected x coordinate value from in the range from to 0 to Length-1.
		int generateRandomX();

		//! \brief Generates a random y coordinate that lies within the y-dimension size of the lattice.
		//! \return
		//! A randomly selected y coordinate value in the range from 0 to Width-1.
		int generateRandomY();

		//! \brief Generates a random z coordinate that lies within the z-dimension size of the lattice.
		//! \return
		//! A randomly selected z coordinate value in the range from 0 to Height-1.
		int generateRandomZ();

		//! \brief Gets the z-direction size of the lattice, the height.
		//! \return The Height property of the lattice, which is the z-direction size.
		int getHeight() const;

		//! \brief Gets the x-direction size of the lattice, the length.
		//! \return The Length property of the lattice, which is the x-direction size.
		int getLength() const;

		//! \brief Gets the number of sites contained in the lattice.
		//! \return The number of sites in the lattice.
		long int getNumSites() const;

		//! \brief Gets the coordinates of the specified site.
		//! \param site_index is the vector index of the input site
		//! \return a Coords object that contains the coordinates of the site specified by the site index.
		Coords getSiteCoords(long int site_index) const;

		//! \brief Gets the vector index for the site corresponding to the input coordinates.
		//! \param coords is the Coords struct that represents the input coordinates.
		//! \return The vector index for the sites vector that is associated with the site located at the input coordinates.
		long int getSiteIndex(const Coords& coords) const;

		//! \brief Gets the vector index for the site corresponding to the input coordinates.
		//! \param x is the x coordinate.
		//! \param y is the y coordinate.
		//! \param z is the z coordinate.
		//! \return The vector index for the sites vector that is associated with the site located at the input coordinates.
		long int getSiteIndex(const int x, const int y, const int z) const;

		// //! \brief Gets the vector iterator for the site corresponding to the input coordinates.
		// //! \param coords is the Coords struct that represents the input coordinates.
		// //! \return The vector iterator for the sites vector that is associated with the site located at the input coordinates.
		// std::vector<Lattice::Site>::iterator getSiteIt(const Coords& coords);

		//! \brief Gets the type of the site wtih the specified site index
		//! \param site_index is the site index.
		//! \return A char datatype indicator of the site type.
		char getSiteType(const long int site_index) const;

		//! \brief Gets the type of the site located at the specified input coordinates
		//! \param coords is the Coords struct that represents the input coordinates.
		//! \return A char datatype indicator of the site type.
		char getSiteType(const Coords& coords) const;

		//! \brief Gets the type of the site located at the specified input coordinates
		//! \param x is the x coordinate.
		//! \param y is the y coordinate.
		//! \param z is the z coordinate.
		//! \return A char datatype indicator of the site type.
		char getSiteType(const int x, const int y, const int z) const;

		//! \brief Gets the lattice unit size, which is used to convert lattice units into real space units.
		//! \return The unit size property of the lattice.
		double getUnitSize() const;

		//! \brief Gets the volume of the lattice in cm^-3.
		double getVolume() const;

		//! \brief Gets the y-direction size of the lattice, the width.
		//! \return The Width property of the lattice, which is the y-direction size.
		int getWidth() const;

		//! \brief Checks whether the x-direction periodic boundaries are enabled or not.
		//! \return true if periodic boundaries are enabled in the x-direction.
		//! \return false if periodic boundaries are disabled in the x-direction.
		bool isXPeriodic() const;

		//! \brief Checks whether the y-direction periodic boundaries are enabled or not.
		//! \return true if periodic boundaries are enabled in the y-direction.
		//! \return false if periodic boundaries are disabled in the y-direction.
		bool isYPeriodic() const;

		//! \brief Checks whether the z-direction periodic boundaries are enabled or not.
		//! \return true if periodic boundaries are enabled in the z-direction.
		//! \return false if periodic boundaries are disabled in the z-direction.
		bool isZPeriodic() const;

		//! \brief Resizes the lattice and clears the sites.
		//! \param length_new is the new x-dimension of the lattice.
		//! \param width_new is the new y-dimension of the lattice.
		//! \param height_new is the new z-dimension of the lattice.
		void resize(const int length_new, const int width_new, const int height_new);

		//! \brief Sets the type of the site located at the specified input coordinates
		//! \param site_index is the site index.
		//! \param site_type is the char datatype designation for the site type.
		void setSiteType(const long int site_index, const char site_type);

		//! \brief Sets the type of the site located at the specified input coordinates
		//! \param x is the x coordinate.
		//! \param y is the y coordinate.
		//! \param z is the z coordinate.
		//! \param site_type is the char datatype designation for the site type.
		void setSiteType(const int x, const int y, const int z, const char site_type);

	protected:

	private:
		bool Enable_periodic_x = true;
		bool Enable_periodic_y = true;
		bool Enable_periodic_z = true;
		int Length = 0; // nm
		int Width = 0; // nm
		int Height = 0; // nm
		double Unit_size = 0.0; // nm
		std::vector<Site> sites;
		std::mt19937_64 gen;
	};
}

#endif // LATTICE_H
