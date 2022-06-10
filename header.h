#include <iostream>
#include <string>
#include <complex>

using namespace std;

class Vector;
class Vector_complex;

class Matrix
{
    int mat_size;
    void init();
    int max_elem(int);
    void change_p(int,int);
    void change_row(int,int);
public:
    double** A;
    int* P;//матрица перестановок
    string n; //номер матрицы
    Matrix(int);
	Matrix(const char*);
	int msize(){return mat_size;}
	void print();
	void p_print();
	void write(string);
	void p_write(string);
	void p_read(const char*);
	Matrix matrix_l();
	Matrix matrix_u();
    void lu();
	void lu_up();
	bool is_matrix_symmetrical();
	void tranpose();
	Matrix chol();
	Matrix mult(const Matrix&);
	Vector mult_vect(const Vector&);
	Vector column(int);
	void add_column(const Vector&,int);
	double off_diag();
	complex <double> complex_l(int ind);
};

class Matrix_complex
{
    int mat_size;
    void init();
    int max_elem(int);
    void change_p(int,int);
    void change_row(int,int);
public:
    complex <double> **A;
    int* P;//матрица перестановок
    string n; //номер матрицы
    Matrix_complex(int);
	Matrix_complex(const char*);
	int msize(){return mat_size;}
	void print();
	void write(string);
	void p_write(string);
	void p_read(const char*);
	Matrix_complex matrix_l();
	Matrix_complex matrix_u();
    void lu();
	void lu_up();
	bool is_matrix_hermit();
	void tranpose();
	Matrix_complex chol();
	Matrix_complex mult(const Matrix_complex&);
	Vector_complex mult_vect(const Vector_complex& V);
	Vector_complex column(int);
	void add_column(const Vector_complex&,int);
	double off_diag();
	complex <double> complex_l(int ind);
};

class Vector
{
    int vect_size;
public:
    double* B;
    Vector(int n);
    Vector(string path);
    int vsize(){return vect_size;}
    void print();
    void write(string);
    Vector operator*=(double scalar);
    Vector operator/=(double scalar);
    Vector operator+ (const Vector& V) const;
    Vector& operator+=(const Vector& V);
    Vector operator-(const Vector& V) const;
    Vector& operator-=(const Vector& V);
    double scal(const Vector& V);
    double norm();
};

class Vector_complex
{
    int vect_size;
public:
    complex <double> *B;
    Vector_complex(int n);
    Vector_complex(const char* path);
    int vsize(){return vect_size;}
    void print();
    void write(string);
    Vector_complex operator*=(complex <double> scalar);
    Vector_complex operator/=(complex <double> scalar);
    Vector_complex operator+ (const Vector_complex& V) const;
    Vector_complex& operator+=(const Vector_complex& V);
    Vector_complex operator-(const Vector_complex& V) const;
    Vector_complex& operator-=(const Vector_complex& V);
    complex <double> scal(const Vector_complex& V);
    complex <double> norm();
    void conjugate();
};
