#include <iostream>
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include "header.h"

Vector_complex::Vector_complex(int n)
{
    vect_size=n;
    B=new complex <double>[vect_size];
    for (int i = 0; i < vect_size; i++)
        B[i] = 0;
}

Vector_complex::Vector_complex(const char* path)
{
    FILE* fp;
    fp = fopen(path, "r");
    char c = ' ';
    double buf;
    vector <double> tmp;
    while (c != '[') fscanf(fp, "%c", &c); //чтение до начала матрицы
    while (c != ']') //чтение первой строки матрицы
    {
        fscanf(fp, "%lf%c", &buf, &c);
        tmp.push_back(buf);
    }
    vect_size = tmp.size(); //размер матрицы равен количеству чисел в первой строке
    B=new complex <double> [vect_size];
    for (int i = 0; i < vect_size; i++) //перенос первой строки в матрицу
        B[i] = tmp[i];
    while (c != '[') fscanf(fp, "%c", &c); //чтение до начала матрицы
    for (int i = 0; i < vect_size; i++)
    {
        fscanf(fp, "%lf%c", &buf, &c);
        B[i]+= buf*1i;
    }
    fclose(fp);
}

void Vector_complex::print()
{
    for (int i=0;i<vect_size;i++)
        cout<<B[i]<<' ';
    cout<<endl;
}

void Vector_complex::write(string path) //запись вектора в файл
{
    ofstream out;
    out.open(path);
    if (!out.is_open()) { cerr<<"Can't open file"<<endl; return;}
    out<<"x = complex([";
   for (int i=0;i<vect_size;i++)
    {
        if (i==vect_size-1) out<<B[vect_size-1].real()<<"],";
        else out<<B[i].real()<<"; ";
    }
    out<<"[";
    for (int i=0;i<vect_size;i++)
    {
        if (i==vect_size-1) out<<B[vect_size-1].imag()<<"]);";
        else out<<B[i].imag()<<"; ";


    }
    out.close();
}

Vector_complex Vector_complex::operator*=(complex <double> scalar)
{
    for (int i=0;i<vect_size;i++)
        B[i]*=scalar;
    return *this;
}

Vector_complex Vector_complex::operator/=(complex <double> scalar)
{
    for (int i=0;i<vect_size;i++)
        B[i]/=scalar;
    return *this;
}

Vector_complex Vector_complex::operator+ (const Vector_complex& V) const
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return *this;}
    Vector_complex res(vect_size);
    for (int i=0;i<vect_size;i++)
        res.B[i]=B[i]+V.B[i];
    return res;
}

Vector_complex& Vector_complex::operator+=(const Vector_complex& V)
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return *this;}
    for (int i=0;i<vect_size;i++)
        B[i]+=V.B[i];
    return *this;
}

Vector_complex Vector_complex::operator- (const Vector_complex& V) const
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return *this;}
    Vector_complex res(vect_size);
    for (int i=0;i<vect_size;i++)
        res.B[i]=B[i]-V.B[i];
    return res;
}

Vector_complex& Vector_complex::operator-=(const Vector_complex& V)
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return *this;}
    for (int i=0;i<vect_size;i++)
        B[i]-=V.B[i];
    return *this;
}

complex <double> Vector_complex::scal(const Vector_complex& V) //скалярное произведение
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return 0;}
    complex <double> res=0;
    for (int i=0;i<vect_size;i++)
        res+=B[i]*V.B[i];
    return res;
}

complex <double> Vector_complex::norm()//норма вектора
{
    complex <double> norm=0;
    for (int i=0;i<vect_size;i++)
        norm+=abs(pow(B[i],2));
    norm=sqrt(norm);
    return norm;
}

void Vector_complex::conjugate()//сопряженный вектор
{
    for (int i=0;i<vect_size;i++)
        B[i]=conj(B[i]);
}
