# gaussian elimination
  Gaussian Elimination Implementation with Partial Pivoting
This program is coded in order to implement gaussian elimination with partial pivoting and backward substitution to solve the fundamental equation Ax=b. A is nxn square matrix and b is nx1 vector. The obtained solution x is also nx1 vector and it is written in the created x.txt file by the program. If the matrix A is singular,c the program gives an error and does not write on x.txt. Also the condition numbers with respect to infinity and one norms are calculated
for the 2x2 matrices only.


In order to use command line arguments in Visual Studio 2022, you must go to
Project ->  project_name Properties -> Configuration Properties -> Debugging -> Command Arguments
then you can write there: A. txt b.txt
This allows A.txt and b.txt becoming the inputs of the program.
They should be in the same folder with the program. Also, x.txt is created in teh same folder with the program.
Then run the code with choosing Debug and clicking start without debugging.

Example of matrix A:                     Example of vector b:             Example of x that is written at x.txt:
9.03260929555454 -8.94646004638415                -4.61761147202888             -0.302706979924999
8.40664079673127 4.75716191033993                 -1.54328769982384              0.210516515086913



*****The Case of High Condition Numbers*****
• Example: Solve Ax=b for both vectors and observe how a small change in b affects the result x.
A:
1.000 1.000
1.000 1.001
b1:   b2:
2.000   2.000
2.000   2.001

The condition number of A with infinity norm: 4004
The condition number of A with one norm: 4004
*calculated by the program.

x1:       while   x2:       
   2                 1
   0                 1
*calculated by the program. 
This shows that even though change is relatively small in input b,
the output changed significantly. This results from the high condition number of A. It means that matrix A
is almost singular.

