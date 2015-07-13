#ifndef MATRIX_H
#define MATRIX_H
#include<iostream>

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
		
		Tbuf *operator[] (int i) {return &buf[i*Ncols]; } //array access operator

	private:
		int Nrows, Ncols;
		Tbuf *buf; //1D array/pointer
};

#endif //MATRIX_H
