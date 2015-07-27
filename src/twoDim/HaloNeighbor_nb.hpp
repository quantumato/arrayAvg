#ifndef HALONEIGHBOR_NB_HPP
#define HALONEIGHBOR_NB_HPP

#include "Halo.hpp"

class HaloNeighbor_nb : public Halo
{
	private:
		MPI_Comm cart;
		MPI_Request* request;
		MPI_Status* status;

		int ndim[2];
		int period[2];

	public:
		HaloNeighbor_nb(int r, int np, int nx, int ny);
		~HaloNeighbor_nb() {};

		bool Halo_Init(matrix<double>& A);
		void Halo_Finalize();
};
#endif
