#include <iostream>
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <errno.h>
#include "header.h"

using namespace std;


void Matrix::init() //инициализация матрицы
{
    A = new double*[mat_size];
    for (int i = 0; i < mat_size; i++)
        A[i] = new double[mat_size];
    P = new int[mat_size];
    for (int i=0;i<mat_size;i++)
        P[i]=i;
}

void Matrix::print() //вывод матрицы на консоль
{
    for (int i = 0; i < mat_size; i++)
    {
        for (int j = 0; j < mat_size; j++)
            cout << A[i][j] <<' ';
        cout<<endl;
    }
}

void Matrix::p_print() //вывод матрицы перестановок на консоль
{
    for (int i=0;i<mat_size;i++)
    cout<<P[i]<<' ';
    cout<<endl;
}

Matrix::Matrix(int n) //создание пустой матрицы
{
    mat_size=n;
    init();
    for (int i=0;i<mat_size;i++)
        for (int j=0;j<mat_size;j++)
            A[i][j]=0;
}

Matrix::Matrix(const char* path) //создание матрицы из файла
{
    vector <double> tmp; //временный вектор для первой строки, чтобы узнать размер матрицы
    char c = ' ';
    float buf;
    FILE* fp;
    fp = fopen(path, "r");
    if (fp == NULL)  cerr<<"There is no file"<<endl;
    while (c != '[') fscanf(fp, "%c", &c); //чтение до начала матрицы
    while (c != ';') //чтение первой строки матрицы
    {
        fscanf(fp, "%f%c", &buf, &c);
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
                fscanf(fp, "%f%c", &buf, &c);
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
    fclose(fp);
    int l=strlen(path); //запись номера матрицы в строку, чтобы потом правильно назвать выходную матрицу
    for (int i = 0; i < l; i++)
        if (isdigit(path[i])) n+=path[i];

}

void Matrix::p_write(string path) //запись P-матрицы в файл
{
    ofstream out;
    out.open(path);
    if (!out.is_open()) {cerr<<"Can't open file"<<endl; return;}
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

void Matrix::write(string path) //запись матрицы в файл
{
    ofstream out;
    out.open(path);
    if (!out.is_open()) {cerr<<"Can't open file"<<endl; return;}
    out<<"A = ...\n [";
    for (int i=0;i<mat_size;i++)
    {
        for (int j=0;j<mat_size-1;j++)
            out<<A[i][j]<<' ';
        if (i==mat_size-1) out<<A[i][mat_size-1]<<"];";
        else out<<A[i][mat_size-1]<<";\n  ";
    }
    out.close();
}

/*LU-разложение*/

void Matrix::lu() //создание матрицы L+U
{
    for (int k=0;k<mat_size-1;k++)
    {
        if (A[k][k]==0) { cerr<<"Division by zero"<<endl; return;}
        for (int i=k+1;i<mat_size;i++)
        {
            A[i][k]/=A[k][k];
            for (int j=k+1;j<mat_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
        }
    }
}

Matrix Matrix::matrix_l() //выделение матрицы L
{
    Matrix L(mat_size);
    for (int i=0;i<mat_size;i++)
    {
        L.A[i][i]=1;
        for (int j=0;j<i;j++)
            L.A[i][j]=A[i][j];
    }
    return L;
}

Matrix Matrix::matrix_u() //выделение матрицы U
{
    Matrix U(mat_size);
    for (int i=0;i<mat_size;i++)
        for (int j=i;j<mat_size;j++)
            U.A[i][j]=A[i][j];
    return U;
}



/*LUP-разложение*/

void Matrix::p_read(const char* path)
{
    P = new int[mat_size];
    char c = ' ';
    int buf;
    FILE* fp;
    fp = fopen(path, "r");
    if (fp == NULL) {cerr<<"There is no file"<<endl; return;}
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

int Matrix::max_elem(int ind) //поиск наибольшого элемента по модулю
{
    double max_elem=fabs(A[ind][ind]);
    int max_ind=ind;
    for (int i=ind+1;i<mat_size;i++)
    {
        if (fabs(A[i][ind])>max_elem)
        {
            max_elem=fabs(A[i][ind]);
            max_ind=i;
        }
    }
    return max_ind;
}

void Matrix::change_p(int ind_1,int ind_2) //изменение матрицы перестановок
{
    int tmp=P[ind_1];
    P[ind_1]=P[ind_2];
    P[ind_2]=tmp;
}

void Matrix::change_row(int ind_1,int ind_2) //обмен рядов в матрице
{
    for (int i=0;i<mat_size;i++)
    {
        float tmp=A[ind_1][i];
        A[ind_1][i]=A[ind_2][i];
        A[ind_2][i]=tmp;
    }
}

void Matrix::lu_up() //создание матрицы L+U (улучшенный)
{
    for (int k=0;k<mat_size-1;k++)
    {
        if (A[k][k]==0) { cerr<<"Division by zero"<<endl; return;}
        int max_ind=max_elem(k);
        if(k!=max_ind)
        {
            change_p(k,max_ind);
            change_row(k,max_ind);
        }
        for (int i=k+1;i<mat_size;i++)
        {
            A[i][k]/=A[k][k];
            for (int j=k+1;j<mat_size;j++)
                A[i][j]-=A[k][j]*A[i][k];
        }
    }
}

/*Разложение Холецкого*/

bool Matrix::is_matrix_symmetrical()
{
    for (int i=0;i<mat_size;i++)
        for (int j=i;j<mat_size;j++)
            if (A[i][j]!=A[j][i]) return false;
    return true;
}

void Matrix::tranpose() //транспонирование матрицы
{
    double tmp;
    for (int i=0;i<mat_size;i++)
        for (int j=i+1;j<mat_size;j++)
        {
            tmp=A[i][j];
            A[i][j]=A[j][i];
            A[j][i]=tmp;
        }
}

Matrix Matrix::chol() //рассчет матрицы С
{
    Matrix C(mat_size);
    for (int i=0;i<mat_size;i++)
    {
        float c=A[i][i];
        for (int j=0;j<mat_size-1;j++)
            c-=pow(C.A[i][j],2);
        if (c>0) C.A[i][i]=sqrtf(c);
        else cerr<<"Matrix is not positive define"<<endl;
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

/*QR-разложение*/

Matrix Matrix::mult(const Matrix& M) //умножение на матрицу
{
    if(mat_size!=M.mat_size) cerr<<"Matrix's size is different"<<endl;
    Matrix res(mat_size);
    for (int i=0;i<mat_size;i++)
        for (int j=0;j<mat_size;j++)
            for (int k=0;k<mat_size;k++)
                res.A[i][j]+=M.A[i][k]*A[k][j];
    return res;
}

Vector Matrix::mult_vect(const Vector& V) //умножение на вектор
{
    Vector res(mat_size);
    for (int i=0;i<mat_size;i++)
        for (int j=0;j<mat_size;j++)
            res.B[i]+=A[i][j]*V.B[j];
    return res;
}

Vector Matrix::column(int ind) //выделение столбца матрицы
{
    Vector V(mat_size);
    for (int i=0;i<mat_size;i++)
        V.B[i]=A[i][ind];
    return V;
}

void Matrix::add_column(const Vector& V,int ind) //запись столбца в матрицу
{
    for (int i=0;i<mat_size;i++)
        A[i][ind]=V.B[i];
}

double Matrix::off_diag() //проверка малости элементов
{
    double e=0;
    for(int i=0;i<mat_size;i++)
        for(int j=i+1;j<mat_size;j++)
            e+=pow(A[j][i],2);
    e=sqrt(e);
    return e;
}

complex <double> Matrix::complex_l(int ind) //вычисление комплексного собственного значения
{
    complex <double> D;//дискриминант
    complex <double> b=(A[ind][ind]+A[ind+1][ind+1])/2;
    complex <double> c=A[ind][ind]*A[ind+1][ind+1]-A[ind][ind+1]*A[ind+1][ind];
    D=pow(b,2)-c;
    complex <double> x=b+sqrt(D);
    return x;
}
