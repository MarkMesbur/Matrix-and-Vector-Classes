//
//  main.cpp
//  SCP1
//
//  Created by Mark Mesbur on 20/10/2018.
//  Copyright © 2018 Mark Mesbur. All rights reserved.
//
#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include "timing.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> // to output data to file


// Class that represents a mathematical vector
class MVector
{
public:
    // constructors
    MVector() {}
    explicit MVector(int n) : v(n) {}
    MVector(int n, double x) : v(n, x) {}
    
    // access element (lvalue) 
    double &operator[](int index) { return v[index]; }
    
    // access element (rvalue) 
    double operator[](int index) const { return v[index]; }
    
    int size() const {return v.size();} // number of elements
    
    double LInfNorm() const // 3.2.2 task 1
    {
        double maxAbs = 0;
        std::size_t s = size();
        for (int i=0; i<s; i++)
        {
            maxAbs = std::max(std::abs(v[i]), maxAbs);
        }
        return maxAbs;
    }
    
    double LsqfNorm() const // 3.2.2 task 2
    {
        double norm = 0;
        std::size_t s = size();
        for (int i=0; i<s; i++)
        {
            norm = norm + v[i]*v[i]; //summing the squared values of the vector
        }
        norm = sqrt(norm); // finally taking the entire summed values root
        return norm;
    }
    
private:
    std::vector<double> v;
};

#endif

// Operator overload for "scalar * vector"
inline MVector operator*(const double& lhs,const MVector& rhs)
{
    MVector temp(rhs);
        for (int i=0;i<temp.size();i++) temp[i]*=lhs;
    return temp;
}

// Operator overload for "vector * scalar"
inline MVector operator*(const MVector& lhs, const double& rhs)
{
    MVector temp(rhs);
        for (int i=0;i<temp.size();i++) temp[i]*=rhs;
    return temp;
}

// Operator overload for "scalar / vector"
inline MVector operator/(const double& lhs,const MVector& rhs)
{
    MVector temp(rhs);
        for (int i=0;i<temp.size();i++) temp[i]/=lhs;
    return temp;
}

// Operator overload for "vector + vector"
inline MVector operator+(const MVector& lhs,const MVector& rhs)
{
    MVector temp(lhs);
        for (int i=0;i<temp.size();i++) temp[i]+=rhs[i];
    return temp;
}

// Operator overload for "vector - vector"
inline MVector operator-(const MVector& lhs,const MVector& rhs)
{
    MVector temp(lhs);
        for (int i=0;i<temp.size();i++) temp[i]-=rhs[i];
    return temp;
}

// Overload the << operator to output MVectors to screen or file
std::ostream& operator<<(std::ostream& os, const MVector& v)
{
    int n = v.size();
    os << "(";
    for (int i=0; i < n; i++) //we less than as the index as the index starts from 0 but size calculates from 1 to ...n
        if (i<n-1) {os <<  v[i] << "," ;} //use of if statement just to remove final comma in vector print out
            else {os <<  v[i];}
    os << ")";
    return os; }


// Operator overload for Dot product of "vector*vector"
double dot(const MVector& lhs, const MVector& rhs) // 3.2.2 task 3
{
    MVector temp(lhs);
    double d = 0;
            for (int i=0;i<temp.size();i++) d = d + temp[i]*rhs[i];
        return d;
}

#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H

#include <vector>

// Class that represents a mathematical matrix
class MMatrix
{
public:
    // constructors
    MMatrix() : nRows(0), nCols(0) {}
    MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) {}
    
    // set all matrix entries equal to a double
    MMatrix &operator=(double x)
    {
        for (int i = 0; i < nRows * nCols; i++) A[i] = x;
        return *this;
    }
    
    // access element, indexed by (row, column) [rvalue]
    double operator()(int i, int j) const
    {
        return A[j + i * nCols];
    }
    
    // access element, indexed by (row, column) [lvalue]
    double &operator()(int i, int j)
    {
        return A[j + i * nCols];
    }
    
    // size of matrix
    int Rows() const { return nRows; }
    int Cols() const { return nCols; }
    
private:
    unsigned int nRows, nCols;
    std::vector<double> A;
};

#endif
    MVector operator*(const MMatrix& A, const MVector& x) //3.2.3 task 4
    {
        MVector b(A.Rows());
        for (int i=0;i< A.Rows();i++)
            for ( int j=0;j< A.Cols();j++)
                b[i]= b[i] + A(i,j)*x[j]; //since the index of vectors starts at 0
        return b;
    };

// Overload the << operator to output MMatrixs to screen, 3.2.3  additional task
std::ostream& operator<<(std::ostream& os, const MMatrix& A)
{
    for (int i=0; i< A.Rows();i++)
        for ( int j=0; j< A.Cols();j++)
            if (j<A.Cols()-1) {     // necessary to use column size minus 1 since the index and size starts from 0 and 1 respectively
                os.width(5);
                os << A(i,j);
                
            }
            else if (j == A.Cols()-1)
            {
                os.width(5);
                os << A(i,j)<< std::endl;
                
            }
    return os; }

MMatrix createM(int i, int j){  // 3.2.4  task 1
    MMatrix A(i,j);
    int x = 0;
    for (int r = 0; r< A.Rows(); r++){
        for (int c=0; c< A.Cols(); c++){
        x = r-c;
        x = abs(x);
            if (r == c) {A(r,c)= 2;}
            else if (x==1) {A(r,c)=-1;};}}
    return A;
}

MMatrix createM2(int i, int j,double m){  // 3.2.4  task 1 where m is the additionally parameter in the task
    MMatrix A(i,j);
    int x = 0;
    for (int r = 0; r< A.Rows(); r++){
        for (int c=0; c< A.Cols(); c++){
            x = r-c;
            x = abs(x);
            if (abs(x)==1) {A(r,c)= -1*(r+2)*(r+2);} //we add two here since the matrix indexing starts rows and cols from 0
            else if (r==c) {A(r,c)=2*(r+2)*(r+2)+m;};}}// as before add two rather than one for required matrix
    return A;
}

//calculate the residues task 3.2.4 task 3
MVector res(MMatrix A,MVector b, MVector x0){
    MVector r(A.Rows());
    r = b - A*x0;
    return r;
};


MVector CGM(MMatrix A,MVector b, MVector x){
    int maxIterations = 1000;
    double tolerance = 1e-6;
    // ...initialise vectors here...
    
    MVector r = res(A, b, x),p = r, r1(A.Rows()),Ap(A.Rows());
    double alpha = 0,scalar = 0,beta = 0,count=0; //we rename
    double startTime = Timer(); // take a time measurement
    for (int iter=0; iter<maxIterations; iter++)
    {
        // ...calculate new values for x and r here...
        // check if solution is accurate enough
        Ap = A*p;  //used so we don't have to do two large matrix vector multiplications for scalar and r
        scalar = dot(p,Ap);            //denominator of alpha
        alpha = dot(r, r)/scalar;       // alpha_k
        x = x + alpha*p; // K+1 of x
        r1 = r;
        r = r - alpha*(Ap);           // neccesary for extra r variable since we need both to calculate beta
        beta = dot(r, r)/dot(r1, r1);
        count = iter;
        if (r.LsqfNorm() < tolerance) break;
        // ...calculate new conjugate vector p here...
        p = r + beta*p;
    }
    double endTime = Timer(); // take another time measurement
    std::cout << count << std::endl;
    std::cout << "took " << endTime-startTime << "s" <<std::endl;
    return x;
};

#ifndef MBANDEDMATRIX_H // the 'include guard'
#define MBANDEDMATRIX_H
#include <vector>
class MBandedMatrix
{
public:
    // constructors
    MBandedMatrix() : nRows(0), nCols(0) {}
    MBandedMatrix(int n, int m, int lband, int rband, double x = 0) :
    nRows(n), nCols(m), A(n * (lband + rband + 1), x), l(lband), r(rband) {}
    
    
    // access element [rvalue]
    double operator()(int i, int j) const
    {
        int index = (j-i+l) + i*(r+l+1);
        if (i<0 || i>=nRows || j<0 || j>=nCols || index<0 || index>=A.size())
        {
            std::cout << "Vector access error" << std::endl;
            throw;
        }
            //(j+i-l) + i*(r + l + 1)
        return A[index];        //
    }
    
    // access element [lvalue]
    double &operator()(int i, int j)
    {
        int index = (j-i+l) + i*(r+l+1);
        if (i<0 || i>=nRows || j<0 || j>=nCols || index<0 || index>=A.size())
        {
            std::cout << "Vector access error" << std::endl;
            throw;
        }
        return A[index]; //
    }
    
    
    // size of matrix
    int Rows() const { return nRows; }
    int Cols() const { return nCols; }
    
    int Bands() const { return r + l + 1; } // total number of bands
    int LBands() const { return l; } // number of left bands
    int RBands() const { return r; } // number of right bands
private:
    unsigned int nRows, nCols;
    std::vector<double> A;
    int l, r; // number of left/right diagonals
};

#endif

std::ostream& operator<<(std::ostream& output, const MBandedMatrix&
                         banded)
{
    int r = banded.Rows(), c = banded.Cols();
    for (int i = 0; i < banded.Rows(); i++)
    {
        // calculate position of lower and upper band
        int jmin = std::max(std::min(i-banded.LBands(), banded.Cols()),0)
        ;
        int jmax = std::min(i+banded.RBands()+1, banded.Cols());
        output << "(";
        for (int j=0; j<jmin; j++)
            output << 0 << "\t";
        for (int j=jmin; j<jmax; j++)
            output << banded(i,j) << "\t";
        for (int j=jmax; j<c; j++)
            output << 0 << "\t";
        output << ")\n";
    }
    return output;
    }

    
MVector operator*(const MBandedMatrix& A, const MVector& x){
    MVector b(A.Rows(),0.0);
        int r = A.Rows(), c = A.Cols();
        for (int i = 0; i < A.Rows(); i++)
        {
            // calculate position of lower and upper band
            int jmin = std::max(std::min(i-A.LBands(), A.Cols()),0);
            int jmax = std::min(i+A.RBands()+1, A.Cols());
            for (int j=jmin; j<jmax; j++)
                b[i] = b[i]+A(i,j)*x[j];
        }
            return b;
        }

    
MBandedMatrix createM3(int i, int j, int l, int r, double p = 0.0){  // 3.2.5  task 2
    MBandedMatrix A(i,j,l,r,0.0);
    for (int b = 0; b< A.Rows(); b++)
    {
        for (int c=0; c< A.Cols(); c++)
        {
            if (b == c) {A(b,c)= 2;}
            else if (b == c+1) {A(b,c)=-1;}
            else if (b == c-1) {A(b,c)=-1;};
        }
    }
    return A;
}

MVector res2(MBandedMatrix A,MVector b, MVector x0){
    MVector r(A.Rows());
    r = b - A*x0;
    return r;
};

MVector CGM2(MBandedMatrix  A,MVector b, MVector x){ // 3.2.5  task 5
    int maxIterations = 1000;
    double tolerance = 1e-6;
    // ...initialise vectors here...
    
    MVector r = res2(A, b, x),p = r, r1(A.Rows()),Ap(A.Rows());
    double alpha = 0,scalar = 0,beta = 0,count=0; //we rename
    double startTime = Timer(); // take a time measurement
    for (int iter=0; iter<maxIterations; iter++)
    {
        // ...calculate new values for x and r here...
        // check if solution is accurate enough
        Ap = A*p; //used so we don't have to do two large matrix vector multiplications for scalar and r
        scalar = dot(p,Ap);            //denominator of alpha
        alpha = dot(r, r)/scalar;       // alpha_k
        x = x + alpha*p; // K+1 of x
        r1 = r;
        r = r - alpha*(Ap);           // neccesary for extra r variable since we need both to calculate beta
        beta = dot(r, r)/dot(r1, r1);
        count = iter;
        if (r.LsqfNorm() < tolerance) break;
        // ...calculate new conjugate vector p here...
        p = r + beta*p;
    }
    double endTime = Timer(); // take another time measurement
    std::cout << count << std::endl;
    std::cout << "took " << endTime-startTime << "s" <<std::endl;
    return x;
};

MMatrix createM4(int n){  // 3.2.5  task 1
    int m = n*n;
    MMatrix A(m,m);  //we replace n^2 with m here and rename it in the function
    int x = 0;
    for (int b = 0; b< A.Rows(); b++)
    {
        for (int c=0; c< A.Cols(); c++)
        {
            x = b-c;
            x = abs(x);
            if (b == c) {A(b,c)= 4;}
            else if (x==n) {A(b,c)=-1;}
            else if (x == 1 and (b+c)%(2*n) != 2*n-1) {A(b,c)=-1;};
        }
    }
    return A;
}

MBandedMatrix createM5(int n){  // 3.2.5  task 3
    int m = n*n;
    MBandedMatrix A(m,m,n,n); //we replace n^2 with n here and rename it in the function
    int x = 0;
    for (int b = 0; b< A.Rows(); b++)
    {
        for (int c=0; c< A.Cols(); c++)
        {
            x = b-c;
            x = abs(x);
            if (b == c) {A(b,c)= 4;}
            else if (x==n) {A(b,c)=-1;}
            else if (x == 1 and (b+c)%(2*n) != 2*n-1) {A(b,c)=-1;};
        }
    }
    return A;
}

int main()
    {

    MVector u(3),w(3),x(3),v(3);
    v[0]=0.1;   v[1]=4.8;   v[2]=3.7;
    w[0]=3.1;   w[1]=8.6;   w[2]=3.6;
    x[0]=5.8;   x[1]=7.4;   x[2]=12.4;
    
    u = 4.7*v + 1.3*w - 6.7*x; //task 4
    
    // inset task 5 here

    std::cout << u[0] << " " << u[1] << " " << u[2] << " " << std::endl; //primitive print out
    std::cout << u << std::endl; // Test of << operator
    std::cout << u.LInfNorm() << std::endl; //it works :)
    
    
    u[0]=1.5;   u[1]=1.3;   u[2]=2.8; //redefining based on required tasks in project
    v[0]=6.5;   v[1]=2.7;   v[2]=2.9;
    w[0]=0.1;   w[1]=-7.2;   w[2]=3.4;
    
    std::cout << dot(u,u)/(dot(v,w)) << std::endl; // 3.2.2 task 4 works correctly
    std::cout << u.LsqfNorm() << ", " << v.LsqfNorm() << ", " << w.LsqfNorm() << std::endl;  //it works :)
    
    MMatrix c(4,3,0);
    c(0,0)=1; c(1,1)=1; c(2,2)=1;c(3,3)=1;
    std::cout << c(0,0) << ", " << c(1,1) << ", "  << c(2,2) << ", " << c(3,3) << std::endl;
    
    MMatrix A(4,3);
    for (int i=0;i<A.Rows();i++)
        for ( int j=0;j<A.Cols();j++)
            A(i,j)=(3*(i))+j;
    
    x[0]=0.5;x[1]=1.6;x[2]=3.2; //redefining x
    MVector b(A.Rows());
    b = A*x; // 3.2.2 task 5
    std::cout <<  b << std::endl;
    
     std::cout << A << std::endl; // Task 3.2.3 additional task print matrix to screen
    
    A = createM(5,5);               //Task 3.2.4 task 2 and redefine A for 3.2.4 task 5
    std::cout << A << std::endl;    // Task 3.2.4 task 2 print createM to screen
     MVector x0(A.Rows(),0);
    
    std::cout << res(A,b,x0) << std::endl;
    MVector bb(5,1.0/(6.0*6.0));
    std::cout  << A*bb << std::endl;
    
        std::cout  << "task 3.2.4 task 5 & 6"<< std::endl;
    MVector b2(10,1.0/(11.0*11.0)); // part of 3.2.4 task 5 & 6 defining b
    A = createM(10,10);
    MVector x2(A.Rows(),0);         // residues task 3.2.4 task 3
    std::cout  << CGM(A, b2, x2) << std::endl;
    std::cout  << "should return b" << A*CGM(A, b2, x2) << std::endl; //test for correctness should return b
        
    MVector b1(25,1.0/(26.0*26.0)); // part of 3.2.4 task 5 & 6 defining b
    A = createM(25,25);
    MVector x1(A.Rows(),0);
    std::cout  << CGM(A, b1, x1) << std::endl;//test
    
    MVector b3(100,1.0/(101.0*101.0)); // part of 3.2.4 task 5 & 6 defining b
    A = createM(100,100);
    MVector x3(A.Rows(),0);
    std::cout  << CGM(A, b3, x3) << std::endl;//test

    std::cout << createM2(5,5,0) << std::endl;// task 3.2.4 task 8
        
        std::cout  << "task 3.2.4 task 8"<< std::endl;
        // part of 3.2.4 task 8
        A = createM2(10,10,100);
        std::cout  << CGM(A, b2, x2) << std::endl;
        std::cout  << A*CGM(A, b2, x2) << std::endl; //test for correctness should return b
        // part of 3.2.4 task 8
        A = createM2(25,25,100);
        std::cout  << CGM(A, b1, x1) << std::endl;
        
        // part of 3.2.4 task 8
        A = createM2(100,100,100);
        std::cout  << CGM(A, b3, x3) << std::endl;
        
    
    MBandedMatrix B4(5,5,1,2,0);  //3.2.5 task 1
    B4(0,0)=1;  B4(0,1)=6; B4(0,2)=10;
    B4(1,0)=13; B4(1,1)=2; B4(1,2)=0; B4(1,3)=11;
    B4(2,1)=14; B4(2,2)=3; B4(2,3)=8; B4(2,4)=12;
    B4(3,2)=0;  B4(3,3)=4; B4(3,4)=9;
    B4(4,3)=16; B4(4,4)=5;
    
    std::cout << B4 << std::endl;//3.2.5 task 1 works
    MVector b4(B4.Rows(),1.0);
    b4 = B4*b4;
    std::cout << b4 << std::endl;
   
    
   MBandedMatrix A5(5,5,1,1);
    A5 = createM3(5,5,1,1);   // 3.2.5  task 3
    MVector b5(A5.Rows(),1.0/(6.0*6.0)),c5(A5.Rows(),0.0);    // 3.2.5  task 3 & 4
    MVector x5(A5.Rows(),0.0);
    std::cout << A5 << std::endl;// 3.2.5  task 2 3.2.5  task 2
    
    c5 = A5*b5;
    std::cout << c5 << std::endl; // 3.2.5  task 3
    std::cout << res2(A5,b5,x5) << std::endl; // 3.2.5  task 4
    
    
    
    
    MBandedMatrix A6 = createM3(25,25,1,1);
    MVector b6(25,1.0/(26.0*26.0)); // part of 3.2.4 task 5 & 6 defining b
    MVector x6(25,0);         // residues task 3.2.4 task 3
    MVector c6 = CGM2(A6, b6, x6);
    std::cout << c6 << std::endl;//test
    std::cout << " if equal to b " <<A6*c6 << std::endl;
   
        std::cout  << "task 3.2.5 task 8"<< std::endl;
        MVector b7(10,1.0/(11.0*11.0)); // part of 3.2.4 task 5 & 6 defining b
        MBandedMatrix A7 = createM3(10,10,1,1);
        MVector x7(A7.Rows(),0);         // residues task 3.2.4 task 3
        std::cout  << CGM2(A7, b7, x7) << std::endl;//test

        
        MVector b8(25,1.0/(26.0*26.0)); // part of 3.2.4 task 5 & 6 defining b
        MBandedMatrix A8 = createM3(25,25,1,1);
        MVector x8(A8.Rows(),0);         // residues task 3.2.4 task 3
        std::cout  <<  CGM2(A8, b8, x8) << std::endl;//test
        
        MVector b9(100,1.0/(101.0*101.0)); // part of 3.2.4 task 5 & 6 defining b
        MBandedMatrix A9 = createM3(100,100,1,1);
        MVector x9(A9.Rows(),0);         // residues task 3.2.4 task 3
        std::cout  << CGM2(A9, b9, x9) << std::endl;//test
        
        
         // 3.2.5 task 1 & 2
        MMatrix A10 = createM4(5);
        MVector x10(A10.Rows(),0);
        MVector b10(A10.Rows(),1.0/(5.0*5.0));
        std::cout << " contour 5 " <<  CGM(A10, b10, x10) << std::endl;//test


        
        MBandedMatrix A11 = createM5(5); //for task 8 section 3.3
        MVector x11(A11.Rows(),0);
        MVector b11(A11.Rows(),1.0/(5.0*5.0));
        std::cout << "bcontour 5" <<  CGM2(A11, b11, x11) << std::endl;//test
        std::cout << A11 << std::endl;
        
        MMatrix A14 = createM4(25); //for task 8 section 3.3
        MVector x14(A14.Rows(),0);
        MVector b14(A14.Rows(),1.0/(26.0*26.0));
        std::cout << " contour 25 " <<  CGM(A14, b14, x14) << std::endl;//test
        
        
        
        MBandedMatrix A15 = createM5(25); //for task 8 section 3.3
        MVector x15(A15.Rows(),0);
        MVector b15(A15.Rows(),1.0/(26.0*26.0));
        std::cout << " bcontour 25 " <<  CGM2(A15, b15, x15) << std::endl;//test
    
        
        //edit of CGM function to output iterations and time for banded matrix below
        
        // From page 29 of lecture 4
        // declare a variable of type std::ofstream named 'file'
        std::ofstream HWFile;
        // try to open the file for writing, specifying the filename
        HWFile.open("banded.txt");
        // if the file failed to open, return from main with nonzero
        if (!HWFile) return 1;
        // write some text to the file; syntax is like std::cout
        
        
        
        int N = 100;
        for (double n = 1.0; n<= N; n++) //use n as a double since we need it's square root later
        {
            MBandedMatrix A12 = createM5(n);
            MVector x12(A12.Rows(),0);
            MVector b12(A12.Rows(),1.0/((n+1)*(n+1)));
            
            int maxIterations = 1000;
            double tolerance = 1e-6;
            // ...initialise vectors here...
            
             MVector r = res2(A12, b12, x12),p = r, r1(A12.Rows()),Ap(A12.Rows());
             double alpha = 0,scalar = 0,beta = 0,count=0; //we rename
             double startTime = Timer(); // take a time measurement
             for (int iter=0; iter<maxIterations; iter++)
             {
                 // ...calculate new values for x and r here...
                 // check if solution is accurate enough
                 Ap = A12*p; //used so we don't have to do two large matrix vector multiplications for scalar and r
                 scalar = dot(p,Ap);            //denominator of alpha
                 alpha = dot(r, r)/scalar;       // alpha_k
                 x12 = x12 + alpha*p; // K+1 of x
                 r1 = r;
                 r = r - alpha*(Ap);           // neccesary for extra r variable since we need both to calculate beta
                 beta = dot(r, r)/dot(r1, r1);
                 count = iter;
                 if (r.LsqfNorm() < tolerance) break;
                 // ...calculate new conjugate vector p here...
                 p = r + beta*p;
             }
            double endTime = Timer(); // take another time measurement
            
            HWFile.width(14);
            HWFile << n;//
            HWFile.width(14);
            HWFile << count;//
            HWFile.width(14);
            HWFile << endTime-startTime << std::endl;//
        }
        
         //edit of CGM function to output iterations and time for matrix class below
        ///From page 29 of lecture 4
        // declare a variable of type std::ofstream named 'file'
        std::ofstream HWFile2;
        // try to open the file for writing, specifying the filename
        HWFile2.open("SCHW2.txt");
        // if the file failed to open, return from main with nonzero
        if (!HWFile2) return 1;
        // write some text to the file; syntax is like std::cout
        
        for (double n = 1.0; n<= N; n++) //use n as a double since we need it's square root later
        
        {
            MMatrix A12 = createM4(n);
            MVector x12(A12.Rows(),0);
            MVector b12(A12.Rows(),1.0/((n+1)*(n+1)));
            
            int maxIterations = 1000;
            double tolerance = 1e-6;
            // ...initialise vectors here...
            
            MVector r = res(A12, b12, x12),p = r, r1(A12.Rows(),0.0),Ap(A12.Rows());
            double alpha = 0.0,scalar = 0.0,beta = 0.0,count=0.0; //we rename
            double startTime = Timer(); // take a time measurement
            for (int iter=0; iter<maxIterations; iter++)
            {
                // ...calculate new values for x and r here...
                // check if solution is accurate enough
                Ap =  A12*p;
                scalar = dot(p,Ap);            //denominator of alpha
                alpha = dot(r, r)/scalar;       // alpha_k
                x12 = x12 + alpha*p; // K+1 of x
                r1 = r;
                r = r - alpha*(Ap);           // neccesary for extra r variable since we need both to calculate beta
                beta = dot(r, r)/dot(r1, r1);
                count = iter;
                if (r.LsqfNorm() < tolerance) break;
                // ...calculate new conjugate vector p here...
                p = r + beta*p;
            }
            double endTime = Timer(); // take another time measurement
            
            HWFile2.width(14);
            HWFile2 << n;//
            HWFile2.width(14);
            HWFile2 << count;//
            HWFile2.width(14);
            HWFile2 << endTime-startTime << std::endl;//
        }
        return 0;

            }
