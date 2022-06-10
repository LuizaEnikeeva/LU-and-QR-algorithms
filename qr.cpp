#include <iostream>
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#define Max 10000
#define eps 0.00001
#include "header.h"


using namespace std;

bool if_file_exist(string path);

/* QR разложение */

/*
    Для того чтобы найти матрицу Q, сделаем ортогонализацию Грамма-Шмидта для заданной матрицы.
    Матрицу R находим из уравнения R=QtA, получаем матрицу, состоящую из скалярных произведений векторов этих матриц.
    В первом цикле берем i-й вектор из матрицы А, а во втором считаем проекции всех уже ортогональных векторов
    и получаем i-й вектор из матрицы Q, в том же цикли кладем скалярное значение в матрицу R.
*/
void qr(Matrix& M,Matrix& Q, Matrix& R)
{
    int n=M.msize();
    Vector a(n);
    Vector p1(n);
    Vector p2(n);
    Vector q(n);
    double scal;
    for (int i=0;i<n;i++)
    {
        a=M.column(i);
        p1=a;
        p2=a;
        for (int j=0;j<i;j++)
        {
            q=Q.column(j);
            scal=q.scal(p1);
            R.A[j][i]=scal;
            q*=scal;
            p2-=q;
        }
        q=p2;
        R.A[i][i]=p2.norm();
        q/=p2.norm();
        Q.add_column(q,i);
    }
}


void qr_complex(Matrix_complex& M,Matrix_complex& Q,Matrix_complex& R)
{
    int n=M.msize();
    Vector_complex a(n);
    Vector_complex p1(n);
    Vector_complex p2(n);
    Vector_complex q(n);
    complex <double> scal;
    for (int i=0;i<n;i++)
    {
        a=M.column(i);
        p1=a;
        p2=a;
        for (int j=0;j<i;j++)
        {
            q=Q.column(j);
            q.conjugate();
            scal=q.scal(p1);
            R.A[j][i]=scal;
            q.conjugate();
            q*=scal;
            p2-=q;
        }
        q=p2;
        R.A[i][i]=p2.norm();
        q/=p2.norm();
        Q.add_column(q,i);
    }
}

/* Вывод матриц в нужные файлы*/

void matrix_qr(const char* path)
{
    Matrix M(path);
    int n=M.msize();
    Matrix Q(n);
    Matrix R(n);
    qr(M,Q,R);
    string out_path="Qmat"+M.n+".m";
    Q.write(out_path);
    out_path="Rmat"+M.n+".m";
    R.write(out_path);
}

void matrix_complex_qr(const char* path)
{
    Matrix_complex M(path);
    int n=M.msize();
    Matrix_complex Q(n);
    Matrix_complex R(n);
    qr_complex(M,Q,R);
    string out_path="Qmat"+M.n+".m";
    Q.write(out_path);
    out_path="Rmat"+M.n+".m";
    R.write(out_path);
}

/* Метод Гаусса*/

/* обратный ход метода Гаусса */
Vector qr_backward(Matrix M, Vector V)
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

Vector_complex qr_backward(Matrix_complex M, Vector_complex V)
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

 /*
    Находим нужный файл с вектором свободных членов. У него должна быть такая же размерность, что и у матрици.
    Пытаемся найти соответствующие файлы с разложением QR, если их нет, то создаем.
    Делим уравнение QRx=b на лва Rx=y Qy=b, так как Q - ортоганальная, то y находим умножая транспонированную
    матрицу Q на вектор b. И делаем обратный ход метода Гаусса, чтобы найти вектор x.
    И записываем этот вектор в нужный файл.
 */
void matrix_qr_gauss(const char* path)
{
    Matrix M(path);
    string path_vect="bvec"+M.n+".m";
    Vector V(path_vect);
    if (M.msize()!=V.vsize())
    {
        cerr<<"Sizes is not equal"<<endl;
        return;
    }
    string path_q="Qmat"+M.n+".m";
    int n=path_q.length();
    char path_char_q[n+1];
    strcpy(path_char_q,path_q.c_str());
    string path_r="Rmat"+M.n+".m";
    n=path_r.length();
    char path_char_r[n+1];
    strcpy(path_char_r,path_r.c_str());
    if (!if_file_exist(path_q) || !if_file_exist(path_r)) matrix_qr(path);
    Matrix Q(path_char_q);
    Matrix R(path_char_r);
    Q.tranpose();
    Vector Y=Q.mult_vect(V);
    Vector X=qr_backward(R,Y);
    string out_path="xvec"+M.n+".m";
    X.write(out_path);
}


void matrix_complex_qr_gauss(const char* path)
{
    Matrix_complex M(path);
    string path_vect="bvec"+M.n+".m";
    int l=path_vect.length();
    char path_char_vect[l+1];
    strcpy(path_char_vect,path_vect.c_str());
    Vector_complex V(path_char_vect);
    if (M.msize()!=V.vsize())
    {
        cerr<<"Sizes is not equal"<<endl;
        return;
    }
    string path_q="Qmat"+M.n+".m";
    int n=path_q.length();
    char path_char_q[n+1];
    strcpy(path_char_q,path_q.c_str());
    string path_r="Rmat"+M.n+".m";
    n=path_r.length();
    char path_char_r[n+1];
    strcpy(path_char_r,path_r.c_str());
    if (!if_file_exist(path_q) || !if_file_exist(path_r)) matrix_complex_qr(path);
    Matrix_complex Q(path_char_q);
    Matrix_complex R(path_char_r);
    Q.tranpose();
    Vector_complex Y=Q.mult_vect(V);
    Vector_complex X=qr_backward(R,Y);
    string out_path="xvec"+M.n+".m";
    X.write(out_path);

}


/*
    Сначала нужно определить с какими числами будем работать: вещественными или комплексными.
    У файлов с комплексными числами должен стоять индефикатор complex.
    В этих функций в зависимости от наличия этого индефикатора выбирается
    правильная функция для создания матриц(для вещественных или комплексных чисел)
*/

void make_qr(const char* path)
{
    char str[12];
    char c;
    char* com="complex";
    FILE* fp;
    fp = fopen(path, "r");
    if (fp == NULL) cerr<< "There is no file" << endl;
    while (c != '=') fscanf(fp, "%c", &c);
    fgets(str,12,fp);
    fclose(fp);
    if (strstr(str,com)==NULL) matrix_qr(path);
    else matrix_complex_qr(path);
}

void qr_gauss(const char* path)
{
    char str[12];
    char c;
    char* com="complex";
    FILE* fp;
    fp = fopen(path, "r");
    if (fp == NULL) cerr<< "There is no file" << endl;
    while (c != '=') fscanf(fp, "%c", &c);
    fgets(str,12,fp);
    fclose(fp);
    if (strstr(str,com)==NULL) matrix_qr_gauss(path);
    else matrix_complex_qr_gauss(path);
}
