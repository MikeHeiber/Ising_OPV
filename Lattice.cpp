// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Ising_OPV project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Ising_OPV project can be found on Github at https://github.com/MikeHeiber/Ising_OPV

#include "Lattice.h"

using namespace std;

Lattice::Lattice() {

}

void Lattice::init(const Parameters_Lattice& params, mt19937* generator_ptr) {
	Enable_periodic_x = params.Enable_periodic_x;
	Enable_periodic_y = params.Enable_periodic_y;
	Enable_periodic_z = params.Enable_periodic_z;
	Length = params.Length;
	Width = params.Width;
	Height = params.Height;
	Unit_size = params.Unit_size;
	Site site;
	sites.assign(Length*Width*Height, site);
	gen_ptr = generator_ptr;
}

void Lattice::calculateDestinationCoords(const Coords& coords_initial, const int i, const int j, const int k, Coords& coords_dest) const{
	coords_dest.x = coords_initial.x + i + calculateDX(coords_initial.x, i);
	coords_dest.y = coords_initial.y + j + calculateDY(coords_initial.y, j);
	coords_dest.z = coords_initial.z + k + calculateDZ(coords_initial.z, k);
}

int Lattice::calculateDX(const int x, const int i) const{
	if (Enable_periodic_x && x + i < 0) {
		return Length;
	}
	else if (Enable_periodic_x && x + i >= Length) {
		return -Length;
	}
	else {
		return 0;
	}
}

int Lattice::calculateDX(const Coords& coords_initial, const Coords& coords_dest) const{
	if (Enable_periodic_x && 2 * (coords_dest.x - coords_initial.x)>Length) {
		return Length;
	}
	else if (Enable_periodic_x && 2 * (coords_dest.x - coords_initial.x)<-Length) {
		return -Length;
	}
	else {
		return 0;
	}
}

int Lattice::calculateDY(const int y, const int j) const{
	if (Enable_periodic_y && y + j<0) {
		return Width;
	}
	else if (Enable_periodic_y && y + j >= Width) {
		return -Width;
	}
	else {
		return 0;
	}
}

int Lattice::calculateDY(const Coords& coords_initial, const Coords& coords_dest) const{
	if (Enable_periodic_y && 2 * (coords_dest.y - coords_initial.y)>Width) {
		return Width;
	}
	else if (Enable_periodic_y && 2 * (coords_dest.y - coords_initial.y)<-Width) {
		return -Width;
	}
	else {
		return 0;
	}
}

int Lattice::calculateDZ(const int z, const int k) const{
	if (Enable_periodic_z && z + k<0) {
		return Height;
	}
	else if (Enable_periodic_z && z + k >= Height) {
		return -Height;
	}
	else {
		return 0;
	}
}

int Lattice::calculateDZ(const Coords& coords_initial, const Coords& coords_dest) const{
	if (Enable_periodic_z && 2 * (coords_dest.z - coords_initial.z) > Height) {
		return Height;
	}
	else if (Enable_periodic_z && 2 * (coords_dest.z - coords_initial.z) < -Height) {
		return -Height;
	}
	else {
		return 0;
	}
}

int Lattice::calculateLatticeDistanceSquared(const Coords& coords_start, const Coords& coords_dest) const{
	int absx = abs(coords_dest.x - coords_start.x);
	int absy = abs(coords_dest.y - coords_start.y);
	int absz = abs(coords_dest.z - coords_start.z);
	int dx, dy, dz;
	if (Enable_periodic_x && 2*absx > Length) {
		dx = -Length;
	}
	else {
		dx = 0;
	}
	if (Enable_periodic_y && 2*absy > Width) {
		dy = -Width;
	}
	else {
		dy = 0;
	}
	if (Enable_periodic_z && 2*absz > Height) {
		dz = -Height;
	}
	else {
		dz = 0;
	}
	return (absx + dx)*(absx + dx) + (absy + dy)*(absy + dy) + (absz + dz)*(absz + dz);
}

bool Lattice::checkMoveValidity(const Coords& coords_initial, const int i, const int j, const int k) const {
	if (i == 0 && j == 0 && k == 0) {
		return false;
	}
	if (!Enable_periodic_x && (coords_initial.x + i >= Length || coords_initial.x + i<0)) {
		return false;
	}
	if (!Enable_periodic_y && (coords_initial.y + j >= Width || coords_initial.y + j<0)) {
		return false;
	}
	if (!Enable_periodic_z && (coords_initial.z + k >= Height || coords_initial.z + k<0)) {
		return false;
	}
	return true;
}

Coords Lattice::generateRandomCoords() {
	Coords coords;
	coords.x = generateRandomX();
	coords.y = generateRandomY();
	coords.z = generateRandomZ();
	return coords;
}

int Lattice::generateRandomX() {
	uniform_int_distribution<int> distx(0, Length - 1);
	auto randx = bind(distx, ref(*gen_ptr));
	return randx();
}

int Lattice::generateRandomY() {
	uniform_int_distribution<int> disty(0, Width - 1);
	auto randy = bind(disty, ref(*gen_ptr));
	return randy();
}

int Lattice::generateRandomZ() {
	uniform_int_distribution<int> distz(0, Height - 1);
	auto randz = bind(distz, ref(*gen_ptr));
	return randz();
}

int Lattice::getHeight() const{
	return Height;
}

int Lattice::getLength() const{
	return Length;
}

long int Lattice::getNumSites() const{
	return (long int)sites.size();
}

Coords Lattice::getSiteCoords(long int site_index) const{
	Coords coords;
	coords.x = site_index / (Width*Height);
	int remainder = site_index % (Width*Height);
	coords.y = remainder / Height;
	coords.z = remainder % Height;
	return coords;
}

long int Lattice::getSiteIndex(const Coords& coords) const{
	return (long int)coords.x*(long int)Width*(long int)Height + (long int)coords.y*(long int)Height + (long int)coords.z;
}

long int Lattice::getSiteIndex(const int x, const int y, const int z) const {
	return (long int)x*(long int)Width*(long int)Height + (long int)y*(long int)Height + (long int)z;
}

vector<Site>::iterator Lattice::getSiteIt(const Coords& coords){
	auto site_it = sites.begin();
	advance(site_it, getSiteIndex(coords));
	return site_it;
}

char Lattice::getSiteType(const long int site_index) const {
	return sites[site_index].type;
}

char Lattice::getSiteType(const Coords& coords) const {
	return sites[getSiteIndex(coords)].type;
}

char Lattice::getSiteType(const int x, const int y, const int z) const {
	return sites[getSiteIndex(x,y,z)].type;
}

double Lattice::getUnitSize() const{
	return Unit_size;
}

int Lattice::getWidth() const{
	return Width;
}

bool Lattice::isXPeriodic() const{
	return Enable_periodic_x;
}

bool Lattice::isYPeriodic() const{
	return Enable_periodic_y;
}

bool Lattice::isZPeriodic() const{
	return Enable_periodic_z;
}

void Lattice::resize(const int length_new, const int width_new, const int height_new) {
	Length = length_new;
	Width = width_new;
	Height = height_new;
	Site site;
	sites.assign(length_new*width_new*height_new, site);
}

void Lattice::setSiteType(const long int site_index, const char site_type) {
	sites[site_index].type = site_type;
}

void Lattice::setSiteType(const int x, const int y, const int z, const char site_type) {
	sites[getSiteIndex(x, y, z)].type = site_type;
}