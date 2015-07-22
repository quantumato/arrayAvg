#ifndef HALONEIGHBOR_HPP
#define HALONEIGHBOR_HPP

#include "Halo.hpp"

//see Halo.hpp for private variables
class HaloNeighbor : public Halo
{
	private:
		double * elementsToSend;
		double * elementsToRecv;
	public:
		HaloNeighbor(matrix<double> A, int r, int np, int nx, int ny);
		~HaloNeighbor();

		bool Halo_Init(matrix<double>& A);
		bool HalO_Finalize();
};

#endif
