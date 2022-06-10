#include <iostream>
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include "header.h"

using namespace std;


bool if_file_exist(string path) //проверяет существует ли файл
{
  ifstream file(path);
  if (file) return true;
  else return false;
}

/*
    В этих функциях создается экземпляр класса матриц, с помощью его методов
    создаются матрицы L и U и записываются в файлы с определенным названием
*/

/*LU-разложение*/
void matrix_lu(const char* path) //lu для вещественых чисел
{
    Matrix M(path);
    M.lu();
    Matrix L=M.matrix_l();
    Matrix U=M.matrix_u();
    string path_out="Lmat"+M.n+".m";
    L.write(path_out);
    path_out="Umat"+M.n+".m";
    U.write(path_out);
}

void matrix_complex_lu(const char* path) //lu для комплексных чисел
{
    Matrix_complex M(path);
    M.lu();
    Matrix_complex L=M.matrix_l();
    Matrix_complex U=M.matrix_u();
    string path_out="Lmat"+M.n+".m";
    L.write(path_out);
    path_out="Umat"+M.n+".m";
    U.write(path_out);
}

/*LUP-разложение*/
void matrix_lu_up (const char* path) //lu_upgrade для вещественых чисел
{
    Matrix M(path);
    M.lu_up();
    Matrix L=M.matrix_l();
    Matrix U=M.matrix_u();
    string path_out="LmatUp"+M.n+".m";
    L.write(path_out);
    path_out="UmatUp"+M.n+".m";
    U.write(path_out);
    path_out="Pmat"+M.n+".m";
    M.p_write(path_out);
}

void matrix_complex_lu_up (const char* path) //lu_upgrade для комплексных чисел
{
    Matrix_complex M(path);
    M.lu_up();
    Matrix_complex L=M.matrix_l();
    Matrix_complex U=M.matrix_u();
    string path_out="LmatUp"+M.n+".m";
    L.write(path_out);
    path_out="UmatUp"+M.n+".m";
    U.write(path_out);
    path_out="Pmat"+M.n+".m";
    M.p_write(path_out);
}


/* решение СЛАУ*/

/*
    В этой функции создается экземпляр класса матрици.
    Открывается файл со столбцем свободных членов с таким же номер как и у матрицы.
    Создается экземпляр класса векторов.
    Создается разложение LUP и с его помошью решается СЛАУ методом Гаусса.
*/

Vector gauss_forward(Matrix M, Vector V)
{
    int n=V.vsize();
    Vector Y(n);
    for (int i=0;i<n;i++)
    {
        Y.B[i]=V.B[M.P[i]];
        for (int j=0;j<i;j++)
            Y.B[i]-=M.A[i][j]*Y.B[j];
        Y.B[i]/=M.A[i][i];
    }
    return Y;
}

Vector gauss_backward(Matrix M, Vector V)
{
    int n=V.vsize();
    Vector X(n);
    for (int i=n-1;i>=0;i--)
    {
        X.B[i]=V.B[i];
        for (int j=n-1;j>i;j--)
            X.B[i]-=M.A[i][j]*X.B[j];
        X.B[i]/=M.A[i][i];
    }
    return X;
}

Matrix check_matrix(string path1, const char* path2) //проверяет есть ли разложение, если нет, то создает разложение, и открывает матрицу
{
    int n=path1.length();
    char path_char[n+1];
    strcpy(path_char,path1.c_str()); //перенос string в char[]
    if (!if_file_exist(path1)) matrix_lu_up(path2);
    Matrix M(path_char);
    return M;
}


void matrix_lu_gauss(const char* path) //lu-gauss для вещественных чисел
{
    Matrix M(path);
    string path_vect="bvec"+M.n+".m";
    Vector V(path_vect);
    if (M.msize()!=V.vsize())
    {
        cerr<<"Sizes are not equal"<<endl;
        return;
    }
    string path_lu="LmatUp"+M.n+".m";
    Matrix L=check_matrix(path_lu,path);
    path_lu="UmatUp"+M.n+".m";
    Matrix U=check_matrix(path_lu,path);
    path_lu="Pmat"+M.n+".m";
    int n=path_lu.length();
    char path_char[n+1];
    strcpy(path_char,path_lu.c_str());
    L.p_read(path_char); //считывается матрица перестановок
    Vector Y=gauss_forward(L,V);
    Vector X=gauss_backward(U,Y);
    string out_path="xvec"+M.n+".m";
    X.write(out_path);
}

/* Для комплексных чисел*/
Vector_complex gauss_complex_forward(Matrix_complex M, Vector_complex V)
{
    int n=V.vsize();
    Vector_complex Y(n);
    for (int i=0;i<n;i++)
    {
        Y.B[i]=V.B[M.P[i]];
        for (int j=0;j<i;j++)
            Y.B[i]-=M.A[i][j]*Y.B[j];
        Y.B[i]/=M.A[i][i];
    }
    return Y;
}

Vector_complex gauss_complex_backward(Matrix_complex M, Vector_complex V)
{
    int n=V.vsize();
    Vector_complex X(n);
    for (int i=n-1;i>=0;i--)
    {
        X.B[i]=V.B[i];
        for (int j=n-1;j>i;j--)
            X.B[i]-=M.A[i][j]*X.B[j];
        X.B[i]/=M.A[i][i];
    }
    return X;
}

Matrix_complex check_complex_matrix(string path1, const char* path2) //проверяет есть разложение и открывает матрицу
{
    int n=path1.length();
    char path_char[n+1];
    strcpy(path_char,path1.c_str()); //перенос string в char[]
    if (!if_file_exist(path1)) matrix_complex_lu_up(path2);
    Matrix_complex M(path_char);
    return M;
}

void matrix_complex_lu_gauss(const char* path) //lu-gauss для вещественных чисел
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
    string path_lu="LmatUp"+M.n+".m";
    Matrix_complex L=check_complex_matrix(path_lu,path);
    path_lu="UmatUp"+M.n+".m";
    Matrix_complex U=check_complex_matrix(path_lu,path);
    path_lu="Pmat"+M.n+".m";
    int n=path_lu.length();
    char path_char[n+1];
    strcpy(path_char,path_lu.c_str());
    L.p_read(path_char);
    Vector_complex Y=gauss_complex_forward(L,V);
    Vector_complex X=gauss_complex_backward(U,Y);
    string out_path="xvec"+M.n+".m";
    X.write(out_path);
}

/*
    Сначала нужно определить с какими числами будем работать: вещественными или комплексными.
    У файлов с комплексными числами должен стоять индефикатор complex.
    В этих функций в зависимости от наличия этого индефикатора выбирается
    правильная функция для создания матриц(для вещественных или комплексных чисел)
*/
void make_lu(const char* path)
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
    if (strstr(str,com)==NULL) matrix_lu(path);
    else matrix_complex_lu(path);
}

void make_lu_upgrade(const char* path)
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
    if (strstr(str,com)==NULL) matrix_lu_up(path);
    else matrix_complex_lu_up(path);
}

void lu_gauss(const char* path)
{
    char str[12];
    char c;
    char* com="complex";
    FILE* fp;
    fp = fopen(path, "r");
    while (c != '=') fscanf(fp, "%c", &c);
    if (fp == NULL) { cerr<<"There is no file"<<endl; return;}
    fgets(str,12,fp);
    fclose(fp);
    if (strstr(str,com)==NULL) matrix_lu_gauss(path);
    else matrix_complex_lu_gauss(path);
}
