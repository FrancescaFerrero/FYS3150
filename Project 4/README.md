# Project 4

Codes for problem 4, 5, 6, 7, 8, 9.
Each code is implemented for the Ising model and uses the Markov Chain Monte Carlo approach to sample spin configurations (Metropolis algorithm). All the physical quantities mentioned below (energy, magnetization etc...) are to be considered per spin.

Plots folder contains plots.


## main_5.cpp

Computes the average energy, average magnetization, specific heat capacity and susceptibility numerically and computes them analytically for a lattice size L=2 and a temperature 1 J/k<sub>B</sub>. 
Saves in text file the running averages of energy and magnetization. 
For L=2, T=1 J/k<sub>B</sub> compares the numerical values with the analytical ones and displays the output to the monitor. 
The user needs to give as input the dimension L of the matrix (LxL), the desired temperature T and a keyword that indicates wheter the initial state should be ordered or random.

Compile and linking: g++ main_5.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

Run: ./main.exe L T "keyword"

where L is the size of the lattice, T is the temperature value and if "keyword" is "ordered" it will give the ordered initial state, while every other word will give the random initial state.

## main_6.cpp

Prints samples of the energy on a file, only after the burn-in time. The initial condition is chosen to be random.
The user needs to give as input the dimension L of the matrix (LxL) and the desired temperature T. 

Compile and linking: g++ main_6.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

Run: ./main.exe L T

where L is the size of the lattice and T is the temperature value.

## main_7.cpp

After reaching the burn-in time, computes and stores in a file average energy, average magnetization, specific heat capacity and susceptibility as function of temperature. The outer loop over temperature is parallelized using OpenMP. 
The initial condition is chosen to be random.
The user needs to give as input the dimension L of the matrix (LxL).

Compile and linking: g++ main_7.cpp utils.cpp -std=c++11 -Xpreprocessor -fopenmp -o main.exe -lomp -larmadillo

Set the number of threads to use by setting the environment variable OMP_NUM_THREADS: e.g. export OMP_NUM_THREADS = 4 

Run: ./main.exe L

where L is the size of the lattice.

## prob_5_plot.ipynb

Takes data from output files of main_5.cpp and main_6.cpp.
Produces plots for running average of energy and magnetization for L=20 and L=2, T=1 J/k<sub>B</sub> and T=2.4 J/k<sub>B</sub> as a function of the number of Monte Carlo cycles, for both random and ordered initial state.

## prob_6_plot.ipynb

Takes data from output files of main_6.cpp.
Produces plots of the probability distribution of the energy per spin by creating normalized histograms of generated energy samples, for T=1 J/k<sub>B</sub> and T=2.4 J/k<sub>B</sub>. 
 
## prob_7_plot.ipynb

Takes data from output files of main_7.cpp.
Produces plots of the average enery, average magnetization, specific heat capacity and susceptibility as function of temperature, for different lattice sizes.
Performs linear regression to verify the scaling relation T<sub>c</sub>(L)-T<sub>c</sub>(L=inifinity)=aL<sup>-1</sup> where T<sub>c</sub>(L) is the critical temperature for a certain lattice size L and T<sub>c</sub>(L=inifinity) is then compared to the analytical result found by Lars Onsager.



