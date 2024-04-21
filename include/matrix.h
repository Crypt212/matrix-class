#ifndef MATRIX_CRYPT_10_10
#define MATRIX_CRYPT_10_10

#include <vector>
#include <cstdint>

namespace m {

	class Matrix {
		public:
			using value_type		= float;
			using length_type		= uint32_t;
			using row_type			= std::vector<value_type>;
			using matrix_type		= std::vector<row_type>;
			enum row_op_type {swap, scale, add_multiple};

			// creates a zero square matrix with t_length rows and columns
			Matrix (length_type t_length = 1);							
			
			// creates a zero matrix with t_rows rows and t_columns columns
			Matrix (length_type t_rows, length_type t_columns);					
			
			// creates a copty of matrix t_matrix
			Matrix (const Matrix& t_matrix);								
			
			// creates a copy submatrix of matrix t_matrix for rows within t_i0 and t_i1, and columns with t_j0 and t_j1
			Matrix (const Matrix& t_matrix, length_type t_i0, length_type t_j0, length_type t_i1, length_type t_j1);
			
			// if square, turns matrix into identity matrix
		   	void identity ();									

			// fills the matrix with a t_constant value
			void fill (value_type t_constant);								
			
			// prints the matrix
			void print () const;								

			// prints a portion of the matrix for rows within t_i0 and t_i1, and columns with t_j0 and t_j1
			void print (length_type t_i0, length_type t_j0, length_type t_i1, length_type t_j1) const;		

			// resize the matrix (note: could remove an augmentation seperator if shrank below its index!)
			void resize (length_type t_new_rows, length_type t_new_columns);				

			// enters values in portion of the matrix for rows within t_i0 and t_i1, and columns with t_j0 and t_j1
			void enter (length_type t_i0, length_type t_j0, length_type t_i1, length_type t_j1);
			
			// enters values in the whole matrix
			void enter ();										

			/* does row operation to the matrix and if determinant exists, it changes it accordingly.
			for t_operations_type set to:
			 	* swap			: swaps rows at t_row0 and t_row1 indices and multiplies determinant by -1
			 	* scale			: multiplies row at t_row0 index by t_multiple and multiplies determinant by t_multiple
			 	* add_multiple	: adds row at t_row0 index multiplied by t_multiple to row at t_row1 index
			*/
			void rowOperation (row_op_type t_operation_type, length_type t_row0 = 0, length_type t_row1 = 0,
				   	value_type t_multiple = 1.0f);

			// sets the matrix to its echelon form, calculates the determinant and returns the matrix
			Matrix& echelon ();									

			// sets the matrix to its reduced echelon form, calculates the determinant and returns the matrix
			Matrix& reduced_echelon ();							
																	
			// sets the matrix to its transpose, and returns it
			Matrix& transpose ();								

			// sets the matrix to its inverse, and returns it
			Matrix& invert ();									

			// creates augment t_matrix to the end of current matrix if number of rows match in both
			void augment (const Matrix& t_matrix);						
																	
			/* reduces the augmented matrix to two seperate matrices at the t_index.. assigns 
			  the second to the parameter matrix &t_matrix if given */
			void unaugment (length_type t_index, Matrix* t_matrix = nullptr);	
																    
			// adds augmentation seperator at t_index
			void addSeperator (length_type t_index);					

			// removes augmentation seperator from t_index
			bool removeSeperator (length_type t_index);					

			// removes all augmentation seperators
			bool removeSeperators ();							

			// sets a cell at t_row and t_column indices to a t_value
			void setCell (length_type t_row, length_type t_column, value_type t_value);

			// calculates the determinant if exists
			void calculateDeterminant ();						

			// returns value of cell at t_row and t_column indices
			value_type getCell (length_type t_row, length_type t_column) const;

			// returns the number of rows
			length_type getRows () const;						

			// returns the number of columns
			length_type getColumns () const;					

			// returns the determinant if there is
			value_type getDeterminant () const;					

			// adds matrix t_matrix to this matrix and returns it
			Matrix& add (const Matrix& t_matrix);						

			// multiplies matrix t_matrix by this matrix and returns it
			Matrix& multiply (const Matrix& t_matrix);					

			// copy-assignment
			Matrix& operator= (const Matrix& t_matrix);					

			// returns copty of the matrix with values multiplied by -1
			Matrix operator- () const;							

			/* adds the two matrices into a new one and returns it (if same dimentions,
			  returns 1x1 zero matrix otherwise) */
			Matrix operator+ (const Matrix& t_matrix) const;				

			/* subtracts the two matrices into a new one and returns it (if same dimentions,
			  returns 1x1 zero matrix otherwise) */
			Matrix operator- (const Matrix& t_matrix) const;				
																	

			/* multiplies the two matrices into a new one and returns it (if columns
			  of first = rows of second, returns 1x1 zero matrix otherwise) */
			Matrix operator* (const Matrix& t_matrix) const;				
																	

			// multiplies the matrix by t_constant value into a new matrix and returns it
			Matrix operator* (const value_type t_constant) const;			

			// sets precision of values in the matrix to t_precision
			void static setPrecision(int t_precision);					

			
			~Matrix ();
		private:
			matrix_type m_matrix {};
			length_type m_rows {0}, m_columns {0};
			std::vector<length_type> m_aug_sep {};	// it stores indices of columns at which matrix is augmented
			value_type * m_determinant {nullptr};	// if nullptr, then there is no determinant (Eg. rows != columns)

			static length_type m_matrices_count;
			static int m_precision;

			void m_defaultConstruct ();				// creates 1x1 zero matrix
	};
}

#endif
