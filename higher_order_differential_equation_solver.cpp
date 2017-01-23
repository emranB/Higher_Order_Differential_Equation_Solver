

/*********************** Question 1 source file ********************/

/*  File: testsolve.cpp 
    test the solve function for the matrix class */
    
#include "matrix.h"
#include "matrix.cpp"

int main(void)
{
    ifstream fin("testsolvein.txt");

    matrix a(2,2), b(2,4), x;
    
    fin >> a >> b;

    cout << "a = " << endl;
    cout << a << endl << endl;

    cout << "b = " << endl;
    cout << b << endl << endl;
            
    x = solve(a,b);

    cout << "x = " << endl;
    cout << x << endl << endl;
        
    cout << "check a * x should be b" << endl;
    cout << a * x << endl << endl;
    
    return 0;
}

matrix solve(matrix a, matrix b)
{
    int p, i, j, k, m, n;
    double temp,factor;

    n = a.rows();
    m = b.cols();
          
    if(a.rows() != b.rows()) {
        cerr << "size mismatch\n";
        exit(1);
    }
   
    if(a.rows() != a.cols()) {
        cerr << "matrix is not square\n";
        exit(1);
    }
        
   
    for(p=0; p < n; p++) {
      
        /* find the largest entry in column p below row p and interchange
           the rows so that the largest entry is in position p */
         
        k = p;
        for(i = p+1; i <a.rows(); i++) {
            if ( fabs(a(i,p)) > fabs( a(k,p) ) ) k = i;
        }
         
        if(k != p) {
            for(j=0; j< n; j++) {
                temp = a(p,j);
                a(p,j) = a(k,j);
                a(k,j) = temp;
            }
         
            for(j=0; j< m; j++)  {
                temp = b(p,j);
                b(p,j) = b(k,j);
                b(k,j) = temp;
            }
        }
      
        if(fabs(a(p,p)) < matrix::tiny) {
            cerr << "the matrix is not invertible\n";
            exit(1);
        }
      
        // divide through row p by a(p,p) in both a and b
        factor = a(p,p);
        for(j = p; j< n; j++) {
            a(p,j) = a(p,j)/factor;
        }
      
        for(j=0; j < m; j++) {
            b(p,j) = b(p,j)/factor;
        }
         
        for(i = 0; i < n; i++) {
            if(i == p) continue;
            factor = a(i,p);
            a(i,p) = 0.0;
   
            for(j = p; j < n; j++) {
                a(i,j) = a(i,j) - factor*a(p,j);
            }
         
            for(j = 0; j < m; j++) {
                b(i,j) = b(i,j) - factor*b(p,j);
            }
        }
    }    
    return b;
}


/*********************** Question 1 output file ********************/

a =
     1.00000     2.00000
     2.00000     3.00000


b =
     5.00000     4.00000     1.00000     2.00000
     2.00000     1.00000     3.00000     4.00000


x =
   -11.00000   -10.00000     3.00000     2.00000
     8.00000     7.00000    -1.00000     0.00000


check a * x should be b
     5.00000     4.00000     1.00000     2.00000
     2.00000     1.00000     3.00000     4.00000


Press any key to continue . . .

/*********************** Question 2 source file ********************/

/*  File: testinverse.cpp 
    test the inverse function for the matrix class */
    
#include "matrix.h"
#include "matrix.cpp"

int main(void)
{
    ifstream fin("testinversein.txt");

    matrix a(2,2), x;
    
    fin >> a;

    cout << "a = " << endl;
    cout << a << endl << endl;

    x = inverse(a);

    cout << "inverse = " << endl;
    cout << x << endl << endl;
            
    cout << "check a * x should be the identity matrix" << endl;
    cout << a * x << endl << endl;
    
    return 0;
}


matrix inverse(const matrix& a)
{
    int n = a.rows();
    matrix b = eye(n);
    matrix y;
   
    y = solve(a, b );

    return y;
}


/*********************** Question 2 output file ********************/

a =
     1.00000     2.00000
     2.00000     3.00000


inverse =
    -3.00000     2.00000
     2.00000    -1.00000


check a * x should be the identity matrix
     1.00000     0.00000
     0.00000     1.00000


Press any key to continue . . .

/*********************** Question 3 source file ********************/

/*  File: testleastsquares.cpp 
    test the least squares function for the matrix class */
    
#include "matrix.h"
#include "matrix.cpp"

int main(void)
{
    ifstream fin("testleastsquaresin.txt");

    matrix a(3,2), b(3,1), x;
    
    fin >> a >> b;

    cout << "a = " << endl;
    cout << a << endl << endl;
    cout << "b = " << endl;
    cout << b << endl << endl;
            
    x = leastsquares(a,b);

    cout << "x = " << endl;
    cout << x << endl << endl;
        
    cout << "check a * x should be as close to b as possible" << endl;
    cout << a * x << endl << endl;
    
    return 0;
}


matrix leastsquares(const matrix& a, const matrix& b)
{
    matrix x;
    matrix a1 = transpose(a) * a;
    matrix b1 = transpose(a) * b;
        
    x = solve(a1, b1);
    
    return x;
}

/*********************** Question 3 output file ********************/

a =
     1.00000     2.00000
     2.00000     3.00000
     1.00000     1.00000


b =
     5.00000
     7.00000
     1.00000


x =
    -2.66667
     4.00000


check a * x should be as close to b as possible
     5.33333
     6.66667
     1.33333


Press any key to continue . . .

/*********************** Question 4 source file ********************/

/*  File: testpseudoinverse.cpp 
    test the pseudoinversefunction for the matrix class */
    
#include "matrix.h"
#include "matrix.cpp"

int main(void)
{
    ifstream fin("testpseudoinversein.txt");

    matrix a(3,2), x;
    
    fin >> a;

    cout << "a = " << endl;
    cout << a << endl << endl;
            
    x = pseudoinverse(a);

    cout << "x = " << endl;
    cout << x << endl << endl;
        
    cout << "check x * a should be the identity matrix" << endl;
    cout << x * a << endl << endl;
    
    return 0;
}

matrix pseudoinverse(const matrix& a)
{
    matrix b;
    
    b = inverse(transpose(a) * a) * transpose(a);
    
    return b;
}

/*********************** Question 4 output file ********************/

a =
     1.00000     2.00000
     2.00000     3.00000
     1.00000     1.00000


x =
    -1.33333     0.33333     1.66667
     1.00000     0.00000    -1.00000


check x * a should be the identity matrix
     1.00000    -0.00000
     0.00000     1.00000


Press any key to continue . . .
