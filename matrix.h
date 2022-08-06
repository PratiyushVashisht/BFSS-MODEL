#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <cstdlib> //seed and random number
#include<random>
using namespace std;

//const int site = 1;
//const int N = 1;  //size of matrix
const int D = 9 ;  //dimension of fields at one site

const int LX = 1; // ?
const int LT = 10; //?

//const int NSCALAR = 9;
const int SITES = 12;
const int NCOLOR = 6;
const int RANK = (NCOLOR*NCOLOR - 1);
const double GAUGETOL = 0.00000000000001;
const int DEGREE = 15;
const double TEMP = 0.3 ;
const double LAMBDA = 1.0/(SITES*SITES*SITES*TEMP*TEMP*TEMP);

// defining class=>

// *****1. Complex 
class Complex
{
private:
    double re,im;
public:
    Complex();    //constructor declaration
    Complex(double,double);
    double real(void) const;    // functions decleration
    double imag(void) const;
    double norm(void);
    void print(void) const;
    
    friend ostream& operator<<(ostream&,Complex); // & bit wise operation
    friend istream& operator<<(istream&,Complex &);
// friend is a function that is not a member of a class but has access to the class's private and protected members.

friend inline Complex conjug(const Complex &o1){
    return(Complex(o1.real(),-o1.imag()));}


inline friend Complex operator +(const Complex &o1, const Complex &o2) {
  return(Complex(o1.real() + o2.real(), o1.imag() + o2.imag()));}

inline friend Complex operator -(const Complex &o1, const Complex &o2) {
  return(Complex(o1.real() - o2.real(), o1.imag() - o2.imag()));}

inline friend Complex operator *(const Complex &o1, const Complex &o2) {
  return(Complex(o1.real() * o2.real() - o1.imag() * o2.imag(),
                 o1.real() * o2.imag() + o1.imag() * o2.real()));}

inline friend Complex operator *(const Complex &o1, const double o2) {
  return(Complex(o1.real() * o2, o1.imag() * o2));}

inline friend Complex operator *(const double o1, const Complex &o2) {
  return(Complex(o2.real() * o1, o2.imag() * o1));}

inline friend Complex operator /(const Complex &o1, const Complex &o2){
  Complex dum; // just Rationalised form 
  double norm;
  norm= o2.real()*o2.real() + o2.imag()*o2.imag();
  dum = Complex((o1.real()*o2.real()+ o1.imag()*o2.imag())/norm,
  (o1.imag()*o2.real()- o1.real()*o2.imag())/norm) ;
  return (dum);
}

inline friend Complex operator /(const Complex &o1, const double o2){
  Complex dum;
  dum = Complex((o1.real()/o2) , o1.imag()/o2) ;
  return (dum);
}


inline friend Complex pow( Complex o1, const int o2){
  Complex c(1,0);
  for (int i=0; i<o2 ; i++)
  c=c*o1;
  return c;
}
};


// **2**Unitary Matrix U******* Traceless,Hermitian,Unitary

    // HErmitian
    
class Matrix{               // hermiitian , traceless?
    private:
       Complex X[SITES][D][NCOLOR][NCOLOR];
    public:
        Matrix();
        Matrix(int ,int,int );
        void print(void) const;
        Complex loc(const int) const ;
        Complex scalar(const int,const int) const ;
        Complex get(int ,int ,int, int) const;  // element of particular row,column
        Complex mult(int,int) const ;
        void set(int,int,int,int, const Complex );
        void pratiyush(int,int , Complex mat[NCOLOR][NCOLOR] );
        friend ostream& operator<<(ostream& out , Matrix m);
        Complex Trace(const Matrix );
};

Matrix operator *(const double & , Matrix &);
Matrix operator -(const Matrix & ,const  Matrix&)   ;
Matrix operator +(const Matrix & ,const  Matrix&)   ;     
Matrix operator *(const  Matrix & , const Matrix &);
Matrix operator /(const  Matrix & , const double);
Matrix pow(const Matrix &, const int);
//****2. U MATRIX
class Umatrix 
{
private:
		Complex mat[NCOLOR][NCOLOR];
public:
		Umatrix ();
		Umatrix(int,int,Matrix );
		Umatrix(Complex [NCOLOR][NCOLOR]);
		Complex get(int, int) const; //The gets() function reads characters from stdin and stores them in str until a newline character or end of file is found.
		void set(int, int, const Complex);
		void print(void);
		friend ostream& operator<<(ostream &, Umatrix);
		friend istream& operator>>(istream &, Umatrix &);
};

Umatrix operator +(const Umatrix &o1, const Umatrix &o2);//MAtrix Operation
Umatrix operator -(const Umatrix &o1, const Umatrix &o2);
Umatrix operator *(const Umatrix &, const Umatrix &);
Umatrix operator *(const Umatrix &, const Complex &); // Matrix & Complex number operation
Umatrix operator *(const Complex &, const Umatrix &);
Umatrix operator *(const Umatrix &, const double);
Umatrix operator *(const double, const Umatrix &);
Umatrix comm(const Umatrix &, const Umatrix &);
Umatrix exp(const Umatrix &u);
Umatrix Adj(const Umatrix &u);
Complex Tr(const Umatrix &);
Umatrix operator /(const Umatrix &o1, const Complex & );
Umatrix operator/ (const Umatrix &o1, const double);

//Umatrix real_gaussian_Umatrix(void);



#endif