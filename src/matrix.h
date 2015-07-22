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
		Tbuf * get_edge(int e) //NOTE: this does not work but is a simplication for demonstration purposes
		{
			//follows order of MPI Neighborhood collectives for consistency
			//NOTE: left and right edges are not continuous segments of memory
			switch (e)
			{
				case 1: //left edge
					e = buf+Ncols+1; //would not actually give left edge
					break;
				case 2: //right edge
					e = buf+Ncols-2; //would not actually give a right edge
					break;
				case 3; //top edge
					e = buf+1;
					break;
				case 4: //bottom edge
					e = buf+(Ncols*Nrows-Ncols);
					break;
				default:
					std::cout << "invalid argument. Arg must be <= 4 or >= 1.\n";
					//should throw exception but I dont' know how to handle that.
					break;
			}
		}

		
		Tbuf *operator[] (int i) {return &buf[i*Ncols]; } //array access operator

	private:
		int Nrows, Ncols;
		Tbuf *buf; //1D array/pointer
};


//ported to matrix.h to make changing array operation easier
/*void average_local(matrix<double>& arr1, matrix<double>& arr2)
{
	//assumes n > 1
	for(int i=2; i<arr1.get_nrows()-2; i++)
		for(int j=2; j<arr2.get_ncols()-2; j++)
			//arr2[i][j]=pow(arr1[i-1][j]*arr1[i+1][j]*arr1[i][j+1]*arr1[i][j-1], 0.25);
			arr2[i][j]=(arr1[i-1][j]+arr1[i+1][j]+arr1[i][j+1]+arr1[i][j-1])/4;
}

void initialize_matrix(matrix<double>& A)
{
	int nx = A.get_nrows();
	int ny = A.get_ncols();
	for(int i=0; i<nx; i++)
		for(int j = 0; j<ny; j++)
			A[i][j] = rand()%10;
}

void average_external(matrix<double>& arr1, matrix<double>& arr2)
{
	//assumes array minimum size of 3
	int nx = arr1.get_nrows();
	int ny = arr1.get_ncols();
	for(int i=1; i<nx-1; i++)
	{
		arr2[1][i]=(arr1[0][i]+arr1[2][i]+arr1[1][i-1]+arr1[1][i+1])/4;
		arr2[nx-2][i]=(arr1[nx-3][i]+arr1[nx-1][i]+arr1[nx-2][i-1]+arr1[nx-2][i+1])/4;
	}
	for(int i = 1; i<ny-1; i++)
	{
		arr2[i][1]=(arr1[i][0]+arr1[i][2]+arr1[i-1][1]+arr1[i+1][1])/4;
		arr2[i][nx-2]=(arr1[i][nx-3]+arr1[i][nx-1]+arr1[i+1][nx-2]+arr1[i-1][nx-2])/4;
	}

}

void printarray(matrix<double>& A)
{
	int nx = A.get_nrows();
	int ny = A.get_nrows();
	for(int i=0; i<nx; i++)
	{
		for(int j = 0; j<ny; j++)
			std::cout << A[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << '\n';
	return;
}*/
#endif //MATRIX_H
