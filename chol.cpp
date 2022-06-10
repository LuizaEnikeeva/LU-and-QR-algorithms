#include <iostream>
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <complex>
#include "Header.h"

using namespace std;

bool if_file_exist(string path);

/* Разложение по метолу Холецкого*/
/*
    Создается матрица, с помощью ее методов строится разложение по методу Холецкого.
    Для этого метода матрица должна быть симметричная.
*/
void matrix_chol(const char* path) //для вещественной матрицы
{
    Matrix M(path);
    if (!M.is_matrix_symmetrical())
    {
        cerr<<"Matrix is not symmetrical"<<endl;
        return;
    }
    Matrix C=M.chol();
    string out_path="Cmat"+M.n+".m";
    C.write(out_path);
}

void matrix_complex_chol(const char* path) //для комплексной матрицы
{
    Matrix_complex M(path);
    if (!M.is_matrix_hermit())
    {
        cerr<<"Matrix is not hermit"<<endl;
        return;
    }
    Matrix_complex C=M.chol();
    string out_path="Cmat"+M.n+".m";
    C.write(out_path);
}

/* Решение СЛАУ с помощью метода Холецкого */
Vector gauss_forward(Matrix M, Vector V); //фуникции из lu

Vector gauss_backward(Matrix M, Vector V);

/*
    Решается также как и в разложении LU
*/
void matrix_chol_gauss(const char* path) //chol_gauss для вещественных чисел
{
    Matrix M(path);
    string path_vect="bvec"+M.n+".m";
    Vector V(path_vect);
    if (M.msize()!=V.vsize())
    {
        cerr<<"Sizes is not equal"<<endl;
        return;
    }
    string path_c="Cmat"+M.n+".m";
    int n=path_c.length();
    char path_char[n+1];
    strcpy(path_char,path_c.c_str());
    if (!if_file_exist(path_c)) matrix_chol(path);
    Matrix C(path_char);
    Vector Y=gauss_forward(C,V);
    C.tranpose();
    Vector X=gauss_backward(C,Y);
    string out_path="xvec"+M.n+".m";
    X.write(out_path);
}

Vector_complex gauss_complex_forward(Matrix_complex M, Vector_complex V); //фуникции из lu

Vector_complex gauss_complex_backward(Matrix_complex M, Vector_complex V);

void matrix_complex_chol_gauss(const char* path) //chol_gauss для комплексных чисел
{
    Matrix_complex M(path);
    string path_vect="bvec"+M.n+".m";
    int l=path_vect.length();
    char path_char_vect[l+1];
    strcpy(path_char_vect,path_vect.c_str());
    Vector_complex V(path_char_vect);
    if (M.msize()!=V.vsize())
    {
        cerr<<"Sizes are not equal"<<endl;
        return;
    }

    string path_c="Cmat"+M.n+".m";
    int n=path_c.length();
    char path_char[n+1];
    strcpy(path_char,path_c.c_str());
    if (!if_file_exist(path_c)) matrix_complex_chol(path);
    Matrix_complex C(path_char);
    Vector_complex Y=gauss_complex_forward(C,V);
    C.tranpose();
    Vector_complex X=gauss_complex_backward(C,Y);
    string out_path="xvec"+M.n+".m";
    X.write(out_path);
}

/*
    Сначала нужно определить с какими числами будем работать: вещественными или комплексными.
    У файлов с комплексными числами должен стоять индефикатор complex.
    В этих функций в зависимости от наличия этого индефикатора выбирается
    правильная функция для создания матриц(для вещественных или комплексных чисел)
*/
void make_chol(const char* path)
{
    char str[12];
    char c;
    char* com="complex";
    FILE* fp;
    fp = fopen(path, "r");
    while (c != '=') fscanf(fp, "%c", &c);
    if (fp == NULL) { cerr<< "There is no file" << endl; return;}
    fgets(str,12,fp);
    fclose(fp);
    if (strstr(str,com)==NULL) matrix_chol(path);
    else matrix_complex_chol(path);
}

void chol_gauss(const char* path)
{
    char str[12];
    char c;
    char* com="complex";
    FILE* fp;
    fp = fopen(path, "r");
    if (fp == NULL) { cerr<< "There is no file" << endl; return;}
    while (c != '=') fscanf(fp, "%c", &c);
    fgets(str,12,fp);
    fclose(fp);
    if (strstr(str,com)==NULL) matrix_chol_gauss(path);
    else  matrix_complex_chol_gauss(path);
}
