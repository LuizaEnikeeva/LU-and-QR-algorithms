#include <iostream>
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include "header.h"

Vector::Vector(int n)
{
    vect_size=n;
    B=new double[vect_size];
    for (int i = 0; i < vect_size; i++)
        B[i] = 0;
}

Vector::Vector(string path)
{
    ifstream file;
    file.open(path);
    if (!file) cerr<< "There is no file" << endl;
    char c = ' ';
    double buf;
    vector <double> tmp;
    while (c != '[') file>>c;
    while (c != ']')
    {
        file>>buf;
        file>>c;
        tmp.push_back(buf);
    }
    vect_size = tmp.size();
    B=new double[vect_size];
    for (int i = 0; i < vect_size; i++)
        B[i] = tmp[i];
    file.close();
}

void Vector::print()
{
    for (int i=0;i<vect_size;i++)
        cout<<B[i]<<' ';
    cout<<endl;
}

void Vector::write(string path) //запись вектора в файл
{
    ofstream out;
    out.open(path);
    if (!out.is_open()) { cerr<<"Can't open file"<<endl; return;}
    out<<"x = [";
    for (int i=0;i<vect_size;i++)
    {
        if (i==vect_size-1) out<<B[i]<<"];";
        else out<<B[i]<<"; ";
    }
    out.close();
}

Vector Vector::operator*=(double scalar)
{
    for (int i=0;i<vect_size;i++)
        B[i]*=scalar;
    return *this;
}

Vector Vector::operator/=(double scalar)
{
    for (int i=0;i<vect_size;i++)
        B[i]/=scalar;
    return *this;
}

Vector Vector::operator+(const Vector& V) const
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return *this;}
    Vector res(vect_size);
    for (int i=0;i<vect_size;i++)
        res.B[i]=B[i]+V.B[i];
    return res;
}

Vector& Vector::operator+=(const Vector& V)
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return *this;}
    for (int i=0;i<vect_size;i++)
        B[i]+=V.B[i];
    return *this;
}

Vector Vector::operator- (const Vector& V) const
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return *this;}
    Vector res(vect_size);
    for (int i=0;i<vect_size;i++)
        res.B[i]=B[i]-V.B[i];
    return res;
}

Vector& Vector::operator-=(const Vector& V)
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return *this;}
    for (int i=0;i<vect_size;i++)
        B[i]-=V.B[i];
    return *this;
}

double Vector::scal(const Vector& V) //скалярное произведение
{
    if (vect_size!=V.vect_size) { cerr<<"Vector's size is different"<<endl; return 0;}
    double res=0;
    for (int i=0;i<vect_size;i++)
        res+=B[i]*V.B[i];
    return res;
}

double Vector::norm() //норма вектора
{
    double norm=0;
    for (int i=0;i<vect_size;i++)
        norm+=pow(B[i],2);
    norm=sqrt(norm);
    return norm;
}

