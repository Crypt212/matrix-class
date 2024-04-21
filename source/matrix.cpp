#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include "matrix.h"

using namespace m;

Matrix::length_type Matrix::m_matrices_count = 0;
int Matrix::m_precision = 3;

Matrix::Matrix (length_type t_length) : m_rows(t_length), m_columns(t_length) {
	m_matrices_count++;
	if (t_length == 0) throw ("Number of rows and columns must be positive!");
	m_matrix.resize(m_rows, row_type{});
	for (auto& row : m_matrix)
		row.resize(m_columns, value_type(0));
	m_determinant = new value_type;
	*m_determinant = value_type{0};
}

Matrix::Matrix (length_type  t_rows, length_type t_columns) : m_rows(t_rows), m_columns(t_columns) { 
	m_matrices_count++;
	if (t_rows == 0 || t_columns == 0) throw ("Number of rows and columns must be positive!");
	m_matrix.resize(m_rows, row_type{});
	for (auto& row : m_matrix)
		row.resize(m_columns, value_type(0));
	if (m_rows == m_columns) {
		m_determinant = new value_type;
		*m_determinant = value_type(0);
	}
}

Matrix::Matrix (const Matrix& t_matrix) : m_rows(t_matrix.m_rows), m_columns(t_matrix.m_columns) { 
	m_matrices_count++;
	m_matrix.resize(m_rows, row_type{});
	for (length_type i = 0; i < m_rows; i++) {
		m_matrix[i] = t_matrix.m_matrix[i];
	}
	if (t_matrix.m_determinant) {
		m_determinant = new value_type;
		*m_determinant = *t_matrix.m_determinant;
	}
	m_aug_sep = t_matrix.m_aug_sep;
}

Matrix::Matrix (const Matrix& t_matrix, length_type t_i0, length_type t_j0, length_type t_i1, length_type t_j1) \
		: m_rows(t_i1 - t_i0 + 1), m_columns(t_j1 - t_j0 + 1) {
		m_matrices_count++;
		if (t_i0 > t_i1 || t_j0 > t_j1)
			throw("Number of rows and columns of submatrix must be positive!");
		if (t_i1 >= t_matrix.m_rows || t_j1 >= t_matrix.m_columns)
			throw("Rows and Columns of submatrix must be contained in the main matrix!");
		m_matrix.resize(m_rows, row_type{});
		for (length_type i = t_i0; i <= t_i1; i++) {
			m_matrix[i - t_i0].resize(m_columns);
			for (length_type j = t_j0; j <= t_j1; j++)
				m_matrix[i - t_i0][j - t_j0] = t_matrix.m_matrix[i][j];
		}
		for (length_type i = 0; i < t_matrix.m_aug_sep.size(); i++) {
			if (t_matrix.m_aug_sep[i] < t_j0) continue;
			if (t_matrix.m_aug_sep[i] >= t_j1) break;
			m_aug_sep.push_back(t_matrix.m_aug_sep[i] - t_j0);
		}
		if (m_determinant) {
			delete m_determinant;
			m_determinant = nullptr;
		}
}

void Matrix::identity () {
	if (m_rows != m_columns) throw ("Not square matrix!");
	for (length_type i = 0; i < m_rows; i++) {
		m_matrix[i].assign(m_columns, value_type(0));
		m_matrix[i][i] = value_type(1);
	}
	if (!m_determinant) m_determinant = new value_type;
	*m_determinant = value_type{1};
}

void Matrix::fill (value_type t_constant) {
	for (length_type i = 0; i < m_rows; i++)
			m_matrix[i].assign(m_columns, t_constant);
	if (m_rows == m_columns) {
		if (!m_determinant) m_determinant = new value_type;
	   	*m_determinant = value_type{0};
	}
}

void Matrix::print () const {
	std::cout << '\n';
	length_type index { 0 };
	std::cout << std::scientific << std::setprecision(m_precision);
	for (length_type i = 0; i < m_rows; i++) {
		std::cout << "[ ";
		for (length_type j = 0; j < m_columns; j++) {
			if (m_matrix[i][j] >= 0) std::cout << ' ';
			std::cout << m_matrix[i][j] << ' ';
			if (!m_aug_sep.empty() && m_aug_sep[index] == j) {
				std::cout << "| ";
				index++;
			}
		}
		index = 0;
		std::cout << "]\n";
	}
	std::cout << '\n';
} 

void Matrix::print (length_type t_i0, length_type t_j0, length_type t_i1, length_type t_j1) const {
	if (t_i0 > t_i1 || t_j0 > t_j1)
		throw("Number of rows and columns must be positive!");
	if (t_i1 >= m_rows || t_j1 >= m_columns)
		throw("Rows and Columns of submatrix must be contained in the main matrix!");
	std::cout << '\n';
	length_type index { 0 };
	std::cout << std::scientific << std::setprecision(m_precision);
	for (length_type i = t_i0; i <= t_i1; i++) {
		std::cout << "[ ";
		for (length_type j = t_j0; j <= t_j1; j++) {
			if (m_matrix[i][j] >= 0) std::cout << ' ';
			std::cout << m_matrix[i][j] << ' ';
			if (!m_aug_sep.empty() && j != t_j1 && m_aug_sep[index] == j) {
				std::cout << "| ";
				index++;
			}
		}
		index = 0;
		std::cout << "]\n";
	}
	std::cout << '\n';
} 

void Matrix::resize (length_type t_new_rows, length_type t_new_columns) {
	if (t_new_rows == m_rows && t_new_columns == m_columns) return;

	// removing out-of-bounds speerators
	if (t_new_columns < m_columns) for (length_type i = 0; i < m_aug_sep.size(); i++) if (m_aug_sep[i] >= t_new_columns - 1) {
		 m_aug_sep.resize(i);
		 break;
	}

	m_rows = t_new_rows;
	m_columns = t_new_columns;
	m_matrix.resize(m_rows, row_type{});
	for (auto& row : m_matrix)
		row.resize(m_columns, value_type(0));
	if (m_determinant) {
		delete m_determinant;
		m_determinant = nullptr;
	}
}

void Matrix::enter (length_type t_i0, length_type t_j0, length_type t_i1, length_type t_j1) {
	if (t_i0 > t_i1 || t_j0 > t_j1)
		throw("Number of rows and columns of submatrix must be positive!");
	if (t_i1 >= m_rows || t_j1 >= m_columns)
		throw("Rows and columns of submatrix must be contained in the main matrix!");
	for (length_type i = t_i0; i <= t_i1; i++) {
		for (length_type j = t_j0; j <= t_j1; j++) {
			std::cout << '[' << i << "][" << j << "]: ";
			std::cin >> m_matrix[i][j];
			while (std::cin.fail()) {
				std::cin.clear();
				std::cin.ignore(std::numeric_limits <std::streamsize>::max(), '\n');
				std::cout << "Invalid input.\n";
				std::cout << '[' << i << "][" << j << "]: ";
				std::cin >> m_matrix[i][j];
			}
			std::cin.clear();
			std::cin.ignore(std::numeric_limits <std::streamsize>::max(), '\n');
		}
	}
	if (m_determinant) {
		delete m_determinant;
		m_determinant = nullptr;
	}
}

void Matrix::enter () {
	for (length_type i = 0; i < m_rows; i++) {
		for (length_type j = 0; j < m_columns; j++) {
			std::cout << '[' << i << "][" << j << "]: ";
			std::cin >> m_matrix[i][j];
			while (std::cin.fail()) {
				std::cin.clear();
				std::cin.ignore(std::numeric_limits <std::streamsize>::max(), '\n');
				std::cout << "Invalid input.\n";
				std::cout << '[' << i << "][" << j << "]: ";
				std::cin >> m_matrix[i][j];
			}
			std::cin.clear();
			std::cin.ignore(std::numeric_limits <std::streamsize>::max(), '\n');
		}
	}
	if (m_determinant) {
		delete m_determinant;
		m_determinant = nullptr;
	}
}

void Matrix::rowOperation (row_op_type t_operation_type, length_type t_row0, length_type t_row1, value_type t_multiple) {
	if (t_row0 >= m_rows || t_row1 >= m_rows) throw("Rows must be contained in the main matrix!");
	switch (t_operation_type) {
		case row_op_type::swap:
			if (t_row0 == t_row1) break;
			m_matrix[t_row0].swap(m_matrix[t_row1]);
			if (m_determinant) *m_determinant *= value_type(-1);
			break;
		case row_op_type::scale:
			for (length_type i = 0; i < m_columns; i++) {
				m_matrix[t_row0][i] *= t_multiple;
			}
			if (m_determinant) *m_determinant *= t_multiple;
			break;
		case row_op_type::add_multiple:
			for (length_type i = 0; i < m_columns; i++) {
				m_matrix[t_row1][i] += m_matrix[t_row0][i] * t_multiple;
			}
			break;
		default:
			break;
	}
}

Matrix& Matrix::echelon () {
	bool non_zero_found = false;
	length_type row_limiter = 0;

	// setup for determinant
	if (m_rows == m_columns) {
		if (!m_determinant) m_determinant = new value_type;
		*m_determinant = value_type(1);
	}

	for (length_type j = 0; j < m_columns && row_limiter < m_rows; j++) {
		row_limiter += (non_zero_found ? 1 : 0);
		non_zero_found = false;
		for (length_type i = row_limiter; i < m_rows; i++) {	// non_zero_found is set to true when first non zero element is met in current column
			if (m_matrix[i][j] == 0) continue;
			if (!non_zero_found) {
				non_zero_found = true;
				if (i != row_limiter) rowOperation(row_op_type::swap, i, row_limiter);
				continue;
			}
			rowOperation(row_op_type::add_multiple, row_limiter, i, value_type(-1) * m_matrix[i][j] / m_matrix[row_limiter][j]);
		}
		if (m_determinant) *m_determinant *= m_matrix[row_limiter][j];
	}
	return *this;
}

Matrix& Matrix::reduced_echelon () {
	echelon();
	length_type column_limiter {m_columns};
	for (length_type row = m_rows - 1; row + 1 != 0; row--) {
		for (length_type j = 0; j < column_limiter; j++) {
			if (m_matrix[row][j] != value_type(0)) {
				column_limiter = j;
				for (length_type i = 0; i < row; i++)
					rowOperation(row_op_type::add_multiple, row, i, value_type(-1) * m_matrix[i][column_limiter] / m_matrix[row][column_limiter]);
			}
		}
	}
	return *this;
}

Matrix& Matrix::transpose () { 
	length_type short_side = 0;
	length_type tall_side = 0;
	value_type temp_val;
	if (m_rows > m_columns) {
		short_side = m_columns;
		tall_side = m_rows;
		m_columns = m_rows;
		for (length_type i = 0; i < short_side; i++)
			m_matrix[i].resize(tall_side);
	} else if (m_columns > m_rows) {
		short_side = m_rows;
		tall_side = m_columns;
		m_rows = m_columns;
		m_matrix.insert(m_matrix.end(), tall_side - short_side, row_type(short_side));
	}
	for (length_type i = 0; i < short_side; i++)
		for (length_type j = i + 1; j < tall_side; j++) {
			temp_val = m_matrix[i][j];
			m_matrix[i][j] = m_matrix[j][i];
			m_matrix[j][i] = temp_val;
		}
	if (short_side == tall_side) return *this; // no need to resize; square matrix
	if (tall_side == m_rows) {
		m_rows = short_side;
		m_matrix.resize(m_rows);
	}
	else if (tall_side == m_columns) {
		m_columns = short_side;
		for (auto& it : m_matrix)
			it.resize(m_columns);
	}
	// determinant stays same if exists
	return *this; 
}

Matrix& Matrix::invert () {
	value_type inverted_determinant;
	if (m_rows != m_columns) throw ("Matrix has no inverse!");
	if (this->getDeterminant() == value_type(0)) throw ("Matrix with zero determinant!");
	inverted_determinant = this->getDeterminant();
	Matrix augmented = *this, inverse = Matrix(augmented.m_rows);
	inverse.identity();
	augmented.augment(inverse);
	augmented.reduced_echelon();

	for (length_type i = 0; i < m_rows; i++) { // normalizing to identity in first matrix
		inverted_determinant /= augmented.m_matrix[i][i];
		augmented.rowOperation(row_op_type::scale, i, 0, value_type(1) / augmented.m_matrix[i][i]);
	}
	augmented.unaugment(0, &inverse);
	if (!inverse.m_determinant)
		inverse.m_determinant = new value_type;
	*inverse.m_determinant = inverted_determinant;
	*this = inverse;
	return *this;
}

void Matrix::augment (const Matrix& t_matrix) {
	if (m_rows != t_matrix.m_rows) throw ("Number of rows doesn't match!");
	m_aug_sep.push_back(m_columns - 1);
	m_aug_sep.insert(m_aug_sep.end(), t_matrix.m_aug_sep.begin(), t_matrix.m_aug_sep.end());
	m_columns += t_matrix.m_columns;
	for (length_type i = 0; i < m_rows; i++)
		m_matrix[i].insert(m_matrix[i].end(), t_matrix.m_matrix[i].begin(), t_matrix.m_matrix[i].end());
	if (m_determinant) {
		delete m_determinant;
		m_determinant = nullptr;
	}
}

void Matrix::unaugment (length_type t_index, Matrix* t_matrix) {
	if (m_aug_sep.empty())
		throw ("Matrix is already unaugmented!");
	if (t_index >= m_aug_sep.size()) 
		throw ("Augmentation index is out of bound!");
	length_type start { m_aug_sep[t_index] + 1 }, end;
	if (m_aug_sep.size() - 1 == t_index) end = m_columns - 1;
	else end = m_aug_sep[t_index + 1];
	if (t_matrix) {
		*t_matrix = Matrix(*this, 0, start, m_rows - 1, end);
	}
	*this = Matrix(*this, 0, 0, m_rows - 1, start - 1);
	if (m_determinant) {
		delete m_determinant;
		m_determinant = nullptr;
	}
}

void Matrix::addSeperator (length_type t_index) {
	for (auto s = m_aug_sep.begin(); s != m_aug_sep.end(); s++) {
		if (*s == t_index) return;
		if (*s > t_index) { m_aug_sep.insert(s - 1, t_index); return; }
	}
	m_aug_sep.push_back(t_index);
}

bool Matrix::removeSeperator (length_type t_index) {
	if (m_aug_sep.empty()) return false;
	for (auto s = m_aug_sep.begin(); s != m_aug_sep.end(); s++) {
		if (*s == t_index) { m_aug_sep.erase(s); return true; }
		if (*s > t_index) break;
	}
	return false;
}

bool Matrix::removeSeperators () {
	if (m_aug_sep.empty()) return false;
	m_aug_sep.clear();
	return true;
}

void Matrix::setCell (length_type t_row, length_type t_column, value_type t_value) {
	if (t_row >= m_rows || t_column >= m_columns) throw ("Indices are out of bound!");
	m_matrix[t_row][t_column] = t_value;
	if (m_determinant) {
		delete m_determinant;
		m_determinant = nullptr;
	}
}

void Matrix::calculateDeterminant () {
	if (m_rows != m_columns) {
		if (m_determinant) {
			delete m_determinant;
			m_determinant = nullptr;
		}
		return;
	}

	if (!m_determinant) m_determinant = new value_type;

	if (m_rows == 1) { *m_determinant = m_matrix[0][0]; return; }
	if (m_rows == 2) { *m_determinant = (m_matrix[0][0] * m_matrix[1][1]) - (m_matrix[0][1] * m_matrix[1][0]); return; }

	// calculating determinant of nxn matrix for n >= 3
	Matrix echelon_formed = *this;
	echelon_formed.echelon();
	*m_determinant = echelon_formed.getDeterminant();
}

Matrix::value_type Matrix::getCell (length_type t_row, length_type t_column) const {
	try { if (t_row >= m_rows || t_column >= m_columns) throw ("Indices are out bound!"); }
	catch (const char *e) {
		std::cerr << "Matrix.getCell: " << e << '\n';
		std::cerr << "Returned zero.\n";
		return value_type{0};
	}
	return m_matrix[t_row][t_column];
}

Matrix::length_type Matrix::getRows () const {
	return m_rows;
}

Matrix::length_type Matrix::getColumns () const {
	return m_columns;
}

Matrix::value_type Matrix::getDeterminant () const {
	if (!m_determinant) {
		Matrix& ref = const_cast<Matrix&>(*this);
		ref.calculateDeterminant();
		if (!m_determinant) throw("Matrix has no determinant!");
	}
	return *m_determinant;
}

Matrix& Matrix::add (const Matrix& t_matrix) {
	*this = *this + t_matrix;
	return *this;
}

Matrix& Matrix::multiply (const Matrix& t_matrix) {
	*this = *this * t_matrix;
	return *this;
}

Matrix& Matrix::operator= (const Matrix& t_matrix) {
	if (this == &t_matrix) return *this;
	// clean up
	if (m_determinant) {
		delete m_determinant;
		m_determinant = nullptr;
	}
	if (!m_aug_sep.empty()) m_aug_sep.clear();
	
	// create
	m_rows = t_matrix.m_rows; m_columns = t_matrix.m_columns;
	m_matrix.resize(m_rows, row_type{});
	for (length_type i = 0; i < m_rows; i++)
		m_matrix[i] = t_matrix.m_matrix[i];
	if (t_matrix.m_determinant) {
		m_determinant = new value_type;
		*m_determinant = *t_matrix.m_determinant;
	}
	if (!t_matrix.m_aug_sep.empty()) m_aug_sep = t_matrix.m_aug_sep;
	return *this;
}

Matrix Matrix::operator- () const {
	Matrix res(m_rows, m_columns);
	for (length_type i = 0; i < m_rows; i++)
		for (length_type j = 0; j < m_columns; j++)
			res.m_matrix[i][j] = - m_matrix[i][j];
	return res;
}

Matrix Matrix::operator+ (const Matrix& t_matrix) const {
	if (m_rows != t_matrix.m_rows || m_columns != t_matrix.m_columns) throw ("Matrices couldn't be added!");
	Matrix res(m_rows, m_columns);
	for (length_type i = 0; i < m_rows; i++)
		for (length_type j = 0; j < m_columns; j++)
			res.m_matrix[i][j] = m_matrix[i][j] + t_matrix.m_matrix[i][j];
	if (m_determinant)
		delete res.m_determinant;
	return res;
}

Matrix Matrix::operator- (const Matrix& t_matrix) const {
	return *this + (-t_matrix);
}

Matrix Matrix::operator* (const Matrix& t_matrix) const {
	if (m_columns != t_matrix.m_rows) throw ("Columns of first matrix not equal to rows of second one!");
	length_type common = m_columns;
	Matrix res(m_rows, t_matrix.m_columns);
	for (length_type i = 0; i < res.m_rows; i++)
		for (length_type j = 0; j < res.m_columns; j++) {
			res.m_matrix[i][j] = 0;
			for (length_type k = 0; k < common; k++)
				res.m_matrix[i][j] += m_matrix[i][k] * t_matrix.m_matrix[k][j];
		}
	if (t_matrix.m_determinant && m_determinant) {
		res.m_determinant = new value_type;
		*res.m_determinant = *t_matrix.m_determinant * *m_determinant;
	}
	return res;
}

Matrix Matrix::operator* (const value_type t_constant) const {
	Matrix res(*this);
	for (auto& row_extern : res.m_matrix)
		for (auto& cell : row_extern)
			cell *= t_constant;
	if (res.m_determinant) {
		res.m_determinant = new value_type;
		*res.m_determinant = t_constant * (*m_determinant) * value_type(m_rows);
	}
	return res;
}

void Matrix::setPrecision (int t_precision) {
	m_precision = t_precision;
}

Matrix::~Matrix () {
	m_matrices_count--;
	if (m_determinant) {
		delete m_determinant;
		m_determinant = nullptr;
	}
}

void Matrix::m_defaultConstruct () {
	m_rows = m_columns = 1;
	m_matrix.resize(1, row_type{value_type{0}});
	if (!m_determinant) m_determinant = new value_type;
	*m_determinant = value_type{1};
	if (!m_aug_sep.empty()) m_aug_sep.clear();
}
