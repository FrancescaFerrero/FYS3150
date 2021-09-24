# Project 2

Codes for problem 3, 4, 5, 6, 7. 

Codes for plotting in problem 6 and 7 are written in Phython, on Jupyter Notebook ("prob_6_plot.ipynb" , "prob_7_plot.ipynb"). 

Plots folder contains plots from problem 6 and 7.

----------------------------------------

__Problem 3__

Sets up a tridiagonal matrix (a,d,a) 6x6, solves linear algebra eigenvalue problem with both Armadillo function and analythical equations. Compare the results with those two methods. 

For macOS users: 
g++ main_3.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

For Linux users:
g++ main_3.cpp utils.cpp -o main.exe -larmadillo

Run: ./main.exe

----------------------------------------

__Problem 4__

Identifies the largest off-diagonal element of a matrix and tests for a particular 4x4 matrix. 

For macOS users: 
g++ main_4.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

For Linux users:
g++ main_4.cpp utils.cpp -o main.exe -larmadillo

Run: ./main.exe

-----------------------------------------

__Problem 5__

Implements Jacobi’s rotation algorithm for solving  linear algebra eigenvalue problem, compares the results with analytical ones. 

For macOS users: 
g++ main_5.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

For Linux users:
g++ main_5.cpp utils.cpp -o main.exe -larmadillo

Run: ./main.exe

-----------------------------------------

__Problem 6a and 6b__

Estimates how the number of required transformations scale with the matrix size 𝑁, for both tridiagonal (6a) and dense matrix  (6b).
Results are plotted together with expected trend using "prob_6_plot.ipynb".

For macOS users: 
g++ main_6a.cpp utils.cpp -std=c++11 -o main.exe -larmadillo
g++ main_6b.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

For Linux users:
g++ main_6b.cpp utils.cpp -o main.exe -larmadillo
g++ main_6b.cpp utils.cpp -o main.exe -larmadillo

Run: ./main.exe

-----------------------------------------

__Problem 7__

Saves in textfiles the x values and the three eigenvectors that correspond to the first three eigenvalues for n=10 and n=100 steps. 
Results are plotted using "prob_7_plot.ipynb".

For macOS users: 
g++ main_7.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

For Linux users:
g++ main_7.cpp utils.cpp -o main.exe -larmadillo

Run: ./main.exe n


where n is desired number of steps, in our case 10 or 100.








