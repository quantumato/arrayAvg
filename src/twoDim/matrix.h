#ifndef MATRIX_H
#define MATRIX_H

/*Program Name: Matrix.h
Author: Kevin Ye
Date: 3/5/15

This is the matrix class definition and implementation. It holds a 2D matrix.
 */
template<typename Tbuf>
class matrix
{
	public:
		matrix() {buf = NULL;} //constructor
		~matrix() {if (buf) delete [] buf;} //destructor

		void assign(int new_Nrows, int new_Ncols) //allocates space for rows and columns
		{
			Nrows = new_Nrows;
			Ncols = new_Ncols;

			buf = new Tbuf[Nrows*Ncols];
		}

		Tbuf *data() const {return buf; } //access to the raw data structure

		int get_nrows() const {return Nrows; } //gets number of columns
		int get_ncols() const {return Ncols; } //gets number of rows
		
		//local values needed for Halo Exchange
		Tbuf * getEdge(int e) //NOTE: simplified concept
		{
			//follows order of MPI Neighborhood collectives for consistency
			//NOTE: left and right edges are not continuous segments of memory
			switch (e)
			{
				case 1: //left edge
					return buf+Ncols+1; //would not actually give left edge
				case 2: //right edge
					return buf+Ncols*2-1; //would not actually give a right edge
				case 3: //top edge
					return buf+1;
				case 4: //bottom edge
					return buf+(Ncols*Nrows-Ncols+1);
				default:
					std::cout << "invalid argument. Arg must be <= 4 or >= 1.\n";
					//should throw exception but I dont' know how to handle that.
					break;
			}
		}

		//external values to be received
		Tbuf * getExternal(int h)
		{
			switch (h)
			{
				case 1: //left edge
					return buf+Ncols; //would not actually give left edge
				case 2: //right edge
					return buf+Ncols*2-1; //would not actually give a right edge
				case 3: //top edge
					return buf+Ncols+1;
				case 4: //bottom edge
					return buf+(Ncols*Nrows-Ncols*2+1);
				default:
					std::cout << "invalid argument. Arg must be <= 4 or >= 1.\n";
					//should throw exception but I dont' know how to handle that.
					break;
			}
		}
		Tbuf * getEdgeArray()
		{
			//some code to return an array of references to the edges
		}
		Tbuf * getExternalArray()
		{
			//some code to return an array of references to the externals
		}


		Tbuf *operator[] (int i) {return &buf[i*Ncols]; } //array access operator

	private:
		int Nrows, Ncols;
		Tbuf *buf; //1D array/pointer
};

#endif //MATRIX_H
