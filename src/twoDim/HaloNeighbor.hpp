#ifndef HALONEIGHBOR_HPP
#define HALONEIGHBOR_HPP

#include "Halo.hpp"

//see Halo.hpp for private variables
class HaloNeighbor : public Halo
{
	private:
		//change size for different number of dimensions
		int ndim[2];
		int period[2];
		MPI_Comm cart;
		MPI_Request request;
		MPI_Status status;
	public:
		HaloNeighbor(int r, int np, int nx, int ny);
		~HaloNeighbor() {};

		bool Halo_Init(matrix<double>& A);
		void HalO_Finalize();
};

#endif
