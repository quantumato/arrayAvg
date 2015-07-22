#ifndef HALONEIGHBOR_NB_HPP
#define HALONEIGHBOR_NB_HPP

#include "Halo.hpp"

class HaloNeighbor_nb : public Halo
{
	HaloNeighbor_nb(matrix<double>& nA, int r, int np, int nx, int ny);
	~HaloNeighbor_nb();

	bool Halo_Init(matrix<double>& A);
	bool Halo_Finalize();
};
#endif
