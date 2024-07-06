#include <iostream>
#include <string>
#include <new>
#include <fstream>
#include <cmath>

using namespace std;

//Function that applies gaussian elimination with partial pivoting.
void gaussian_elimination(double** a, double* b, int n) {
    for (int j = 0; j < n - 1; j++) {         // starting from the first column j=0, applying the partial pivoting part.
        int max_row_index = j;
        for (int i = j + 1; i < n; i++) {    //looks below the pivot position to find whether higher value is there.
            if (abs(a[i][j]) > abs(a[max_row_index][j]))
                max_row_index = i;                 // the index of the maximum value in the row is found.
        }
        if (max_row_index != j) {   //change the order of the rows if the index has changed.

            double* temp = a[j];
            a[j] = a[max_row_index];   //change the row order of A.
            a[max_row_index] = temp;

            double temp2 = b[j];
            b[j] = b[max_row_index];   //change the row order of b.
            b[max_row_index] = temp2;
        }
        for (int i = j + 1; i < n; i++) { // gaussian eliminiation part, the formula applied for partial pivoting.
            double scaling = a[i][j] / a[j][j]; 
            for (int k = j; k < n; k++) {

                a[i][k] = a[i][k] - scaling * a[j][k]; //makes entries zero below the pivot.
            }
            b[i] = b[i] - scaling * b[j]; // same procedure also apllied to the b vector.
        }
    }
}

// Function that applies the backward substitution part.
double* backward_subs(double** a, double* b, int n) {    //return type of the function is matched with type of the function.
    for (int j = n - 1; j >= 0; j--) {     // loop starts from the bottom and goes to the top for the backwards substitution.
        if (abs(a[j][j]) < pow(10.0, -10.0)) {           // checking for singularity, if detects gives an error and quit.
            cout << "Error! Given matrix A is singular";
            return 0;                      //quits the function if it's singular.
        }
        b[j] = b[j] / a[j][j];            //formula for the backward substitution is applied.
        for (int i = 0; i < j; i++) {
            b[i] = b[i] - b[j] * a[i][j];       // updating the b to acquire x.
        }
    }
    return b;
}

int main(int argc, char* argv[]) {// main function takes the required comment line arguments.

    // Control the input number.
    if (argc != 3) {  //argument count argc must be 3 since argv[0] will always contain the name of the program.
        cout << "Wrong number of inputs,please check your inputs.";
    }

    // Read the A.txt file first to get the column (or row since it's square) number "n".
    ifstream mfile;
    mfile.open(argv[1]);
    string line;

    int n = 0;
    // Determine the column(or row) number "n" by counting the lines of the input matrix A.
    if (mfile.is_open()) {
        while (getline(mfile, line)) {
            n = n + 1;
        }
    }
    else cout << "Unable to open the file";

    // Resetting the mfile.
    mfile.clear();
    mfile.seekg(0);

    // After reading the A.txt, now we have created the nxn A matrix dynamically with the values of A.txt.

    double** matrixA;
    matrixA = new double* [n];
    for (int i = 0; i < n; i++) {
        matrixA[i] = new double[n];
        for (int j = 0; j < n; j++) {
            mfile >> matrixA[i][j];
        }
    }
    mfile.close();  //closed the file.

    // Read the b.txt input file.
    mfile.open(argv[2]);

    //Reset and go to the beginning of the nfile.
    mfile.clear();
    mfile.seekg(0);

    // Created dynamically allocated b vector by reading b.txt earlier.
    double* vectorb;
    vectorb = new double[n];
    for (int i = 0; i < n; i++) {
        mfile >> vectorb[i];
    }
    mfile.close(); //closed the file.

    //Part where condition number of 2x2 matrix is calculated, with respect to inifinte norm and one norm.

    if (n == 2) {
        double det = matrixA[0][0] * matrixA[1][1] - matrixA[1][0] * matrixA[0][1]; //determinant value of A is calculated.
        if (det == 0)
            cout << "matrix A is singular, therefore condition number is infinite.\n";  //stated if it is singular, cond number is infinite.
        else {
            double Anorminf = abs(matrixA[0][0]) + abs(matrixA[0][1]);   //nonsingular, we can calculate A and A inverse norms and take the product.                                                               
            if (Anorminf < abs(matrixA[1][0]) + abs(matrixA[1][1]))     //infinity norm of A is calculated.(abolsute row sum)
                Anorminf = abs(matrixA[1][0]) + abs(matrixA[1][1]);

            double inversedet = 1 / det;                 // entries of A inverse are calculated.
            double inv00 = inversedet * matrixA[1][1];
            double inv11 = inversedet * matrixA[0][0];
            double inv01 = -inversedet * matrixA[0][1];
            double inv10 = -inversedet * matrixA[1][0];

            double invAnorminf = abs(inv00) + abs(inv01);   //infinity norm of A inverse is calculated.
            if (invAnorminf < abs(inv10) + abs(inv11))
                invAnorminf = abs(inv10) + abs(inv11);

            double cond_number_inf = Anorminf * invAnorminf;  //take the product of infinity norms of A and A inverse to calculate cond number.


            double Anorm_one = abs(matrixA[0][0]) + abs(matrixA[1][0]);
            if (Anorm_one < abs(matrixA[0][1]) + abs(matrixA[1][1]))     //one norm of A is calculated.(absolute column sum)
                Anorm_one = abs(matrixA[0][1]) + abs(matrixA[1][1]);

            double invAnorm_one = abs(inv00) + abs(inv10);   //one norm of A inverse is calculated.
            if (invAnorm_one < abs(inv01) + abs(inv11))
                invAnorm_one = abs(inv01) + abs(inv11);

            double cond_number_one = Anorm_one * invAnorm_one; //take the product of one norms of A and A inverse to calculate cond number.


            cout << "The condition number of A with infinity norm: " << cond_number_inf << endl; //prints out the condition numbers.
            cout << "The condition number of A with one norm: " << cond_number_one << endl;
            cout << endl;
        }
    }


    gaussian_elimination(matrixA, vectorb, n);     //function is called to apply elimination on the matrix A and vector b.
    if (backward_subs == 0)
        return 0;   //function returns 0 if it is singular,with this code we also return 0 in main function and quit the code without further process.
    
    double* x = backward_subs(matrixA, vectorb, n); //solution vector points out the return values of the backward substitution function.

    for (int i = 0; i < n; i++) { //to print the solution vector x.
        cout << x[i] << endl;    //dereference of pointer x gives the updated b values by the function.
    }

    ofstream solfile;            //writes the solution vector x in a x.txt file.
    solfile.open("x.txt");
    for (int i = 0; i < n; i++) {
        solfile << x[i] << endl;
    }
    solfile.close();


    for (int i = 0; i < n; i++) {  // deletes dynamically allocated A matrix.
        delete[] matrixA[i];
    }
    delete[] matrixA;
    delete[] vectorb;// deletes dynamically allocated b vector.
    return 0;
}