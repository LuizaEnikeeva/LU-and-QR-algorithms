#include <iostream>
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <complex>
#include "header.h"

using namespace std;

void Matrix_complex::init() //инициализация матрицы
{
    A = new complex <double>* [mat_size];
    for (int i = 0; i < mat_size; i++)
        A[i] = new complex <double> [mat_size];
    P = new int[mat_size];
    for (int i=0;i<mat_size;i++)
        P[i]=i;
}

void Matrix_complex::print()
{
    for (int i = 0; i < mat_size; i++)
    {
        for (int j = 0; j < mat_size; j++)
            cout << A[i][j]<<' ';
        cout<<endl;
    }
}

Matrix_complex::Matrix_complex(int n)
{
    mat_size=n;
    init();
    for (int i=0;i<mat_size;i++)
        for (int j=0;j<mat_size;j++)
            A[i][j]=0;
}

Matrix_complex::Matrix_complex(const char* path)
{
    FILE* fp;
    fp = fopen(path, "r");
    char c = ' ';
    double buf;
    vector <double> tmp;
    while (c != '[') fscanf(fp, "%c", &c); //чтение до начала матрицы
    while (c != ';') //чтение первой строки матрицы
    {
        fscanf(fp, "%lf%c", &buf, &c);
        tmp.push_back(buf);
    }
    mat_size = tmp.size(); //размер матрицы равен количеству чисел в первой строке
    init();
    for (int i = 0; i < mat_size; i++) //перенос первой строки в матрицу
        A[0][i] = tmp[i];
    while (c != ']')
    {
        for (int i = 1; i < mat_size; i++)
        {
            for (int j = 0; j < mat_size; j++)
            {
                fscanf(fp, "%lf%c", &buf, &c);
                A[i][j] = buf;
            }
            if (c != ';' && c != ']') //в строках разное количество чисел
            {
                cerr<<"Not a matrix"<<endl;
                break;
            }
        }
        if (c != ']') //матрица не квадратная
        {
            cerr<<"Matrix is not square"<<endl;
            break;
        }
    }
    while (c != '[') fscanf(fp, "%c", &c); //чтение до начала матрицы
    while (c != ']')
    {
        for (int i = 0; i < mat_size; i++)
        {
            for (int j = 0; j < mat_size; j++)
            {
                fscanf(fp, "%lf%c", &buf, &c);
                A[i][j] += buf*1i;
            }
            if (c != ';' && c != ']') //в строках разное количество чисел
            {
                cerr<<"Not a matrix"<<endl;
                break;
            }
        }
        if (c != ']') //матрица не квадратная
        {
            cerr<<"Matrix is not square"<<endl;
            break;
        }
    }
    fclose(fp);
    int l=strlen(path); //запись номера матрицы в строку
    for (int i = 0; i < l; i++)
    {
        if (isdigit(path[i]))
        {
            n+=path[i];
        }
    }
}

void Matrix_complex::p_write(string path) //запись матрицы в файл
{
    ofstream out;
    out.open(path);
    if (!out.is_open()) cerr<<"Can't open file"<<endl;
    out<<"P = ...\n [";
    int k;
    for (int i=0;i<mat_size;i++)
    {
        for(int j=0;j<mat_size-1;j++)
        {
            k=j==P[i];
            out<<k<<' ';
        }
        k=mat_size-1==P[i];
        if (i==mat_size-1) out<<k<<"];";
        else out<<k<<";\n  ";
    }
    out.close();
}

void Matrix_complex::write(string path) //запись матрицы в файл
{
    ofstream out;
    out.open(path);
    if (!out.is_open()) { cerr<<"Can't open file"<<endl; return;}
    out<<"A = complex([";
    for (int i=0;i<mat_size;i++)
    {
        for (int j=0;j<mat_size-1;j++)
            out<<A[i][j].real()<<' ';
        if (i==mat_size-1) out<<A[i][mat_size-1].real()<<"],";
        else out<<A[i][mat_size-1].real()<<";\n   ";
    }
    out<<"[";
    for (int i=0;i<mat_size;i++)
    {
        for (int j=0;j<mat_size-1;j++)
            out<<A[i][j].imag()<<' ';
        if (i==mat_size-1) out<<A[i][mat_size-1].imag()<<"]);";
        else out<<A[i][mat_size-1].imag()<<";\n   ";
    }
    out.close();
}

//LU-разложение

Matrix_complex Matrix_complex::matrix_l() //выделение матрицы L
{
    Matrix_complex L(mat_size);
    for (int i=0;i<mat_size;i++)
    {
        L.A[i][i]=1;
        for (int j=0;j<i;j++)
            L.A[i][j]=A[i][j];
    }
    return L;
}

Matrix_complex Matrix_complex::matrix_u() //выделение матрицы U
{
    Matrix_complex U(mat_size);
    for (int i=0;i<mat_size;i++)
        for (int j=i;j<mat_size;j++)
            U.A[i][j]=A[i][j];
    return U;
}

void Matrix_complex::lu() //создание матрицы L+U
{
    complex <double> zero;
    for (int k=0;k<mat_size-1;k++)
    {
        if (A[k][k]==zero) { cerr<<"Division by zero"<<endl; return;}
        for (int i=k+1;i<mat_size;i++)
        {
            A[i][k]/=A[k][k];
            for (int j=k+1;j<mat_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
        }
    }

}

//LUP-разложение

void Matrix_complex::p_read(const char* path)
{
    P = new int[mat_size];
    char c = ' ';
    int buf;
    FILE* fp;
    fp = fopen(path, "r");
    if (fp == NULL) cerr<< "There is no file" << endl;
    while (c != '[') fscanf(fp, "%c", &c); //чтение до начала матрицы
    for(int i=0;i<mat_size;i++)
    {
        for(int j=0;j<mat_size;j++)
        {
            fscanf(fp, "%d%c", &buf, &c);
            if(buf==1) P[i]=j;
        }
    }
    fclose(fp);
}

int Matrix_complex::max_elem(int ind) //поиск наибольшого элемента
{
    double max_elem=abs(A[ind][ind]);
    int max_ind=ind;
    for (int i=ind+1;i<mat_size;i++)
    {
        if (abs(A[i][ind])>max_elem)
        {

            max_elem=abs(A[i][ind]);
            max_ind=i;
        }
    }
    return max_ind;
}

void Matrix_complex::change_p(int ind_1,int ind_2) //изменение матрицы перестановок
{
    int tmp=P[ind_1];
    P[ind_1]=P[ind_2];
    P[ind_2]=tmp;
}

void Matrix_complex::change_row(int ind_1,int ind_2) //обмен рядов в матрице
{
    for (int i=0;i<mat_size;i++)
    {
        complex <double> tmp=A[ind_1][i];
        A[ind_1][i]=A[ind_2][i];
        A[ind_2][i]=tmp;
    }
}

void Matrix_complex::lu_up() //создание матрицы L+U (улучшенный)
{
    complex <double> zero;
    for (int k=0;k<mat_size-1;k++)
    {
        if (A[k][k]==zero) { cerr<<"Division by zero"<<endl; return;}
        int max_ind=max_elem(k);
        change_p(k,max_ind);
        change_row(k,max_ind);
        for (int i=k+1;i<mat_size;i++)
        {
            A[i][k]/=A[k][k];
            for (int j=k+1;j<mat_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
        }
    }
}

//Разложение Холецкого

bool Matrix_complex::is_matrix_hermit()
{
    for (int i=0;i<mat_size;i++)
        for (int j=i+1;j<mat_size;j++)
            if (A[i][j]!=conj(A[j][i])) return false;
    return true;
}

void Matrix_complex::tranpose() //транспонирование матрицы
{
    complex <double> tmp;
    for (int i=0;i<mat_size;i++)
        for (int j=i+1;j<mat_size;j++)
        {
            tmp=conj(A[i][j]);
            A[i][j]=conj(A[j][i]);
            A[j][i]=tmp;
        }
}

Matrix_complex Matrix_complex::chol() //рассчет матрицы С
{
    Matrix_complex C(mat_size);
    for (int i=0;i<mat_size;i++)
    {
        complex <double> c=A[i][i];
        for (int j=0;j<mat_size-1;j++)
            c-=C.A[i][j]*conj(C.A[i][j]);
        C.A[i][i]=sqrt(c);
        for (int k=i+1;k<mat_size;k++)
        {
            c=A[k][i];
            for (int j=0;j<mat_size-1;j++)
                c-=C.A[i][j]*C.A[k][j];
            C.A[k][i]=c/C.A[i][i];
        }
    }
    return C;
}

//QR-разложение

Matrix_complex Matrix_complex::mult(const Matrix_complex& M) //умножение на матрицу
{
    if(mat_size!=M.mat_size) cerr<<"Matrix's size is different"<<endl;
    Matrix_complex res(mat_size);
    for (int i=0;i<mat_size;i++)
        for (int j=0;j<mat_size;j++)
            for (int k=0;k<mat_size;k++)
                res.A[i][j]+=M.A[i][k]*A[k][j];
    return res;
}

Vector_complex Matrix_complex::mult_vect(const Vector_complex& V)//умножение на вектор
{
    Vector_complex res(mat_size);
    for (int i=0;i<mat_size;i++)
        for (int j=0;j<mat_size;j++)
            res.B[i]+=A[i][j]*V.B[j];
    return res;
}

Vector_complex Matrix_complex::column(int ind)//выделение столбца матрицы
{
    Vector_complex V(mat_size);
    for (int i=0;i<mat_size;i++)
        V.B[i]=A[i][ind];
    return V;
}

void Matrix_complex::add_column(const Vector_complex& V,int ind)//запись столбца в матрицу
{
    for (int i=0;i<mat_size;i++)
        A[i][ind]=V.B[i];
}

double Matrix_complex::off_diag() //проверка малости элементов
{
    double e=0;
    for(int i=0;i<mat_size;i++)
        for(int j=i+1;j<mat_size;j++)
            e+=abs(pow(A[j][i],2));
    e=sqrt(e);
    return e;
}

complex <double> Matrix_complex::complex_l(int ind) //вычисление комплексного собственного значения
{
    complex <double> two(0,2);
    complex <double> D;//дискриминант
    complex <double> b=(A[ind][ind]+A[ind+1][ind+1])/two;
    complex <double> c=A[ind][ind]*A[ind+1][ind+1]-A[ind][ind+1]*A[ind+1][ind];
    D=pow(b,2)-c;
    complex <double> x=b+sqrt(D);
    return x;
}
