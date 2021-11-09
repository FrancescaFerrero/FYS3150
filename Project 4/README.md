# Project 4

Codes for problem 4, 5, 6 and 8.

Codes for the plots in problem 5 and 8 are written in Phython, on Jupyter Notebook ("prob_5_plot.ipynb", "prob_8_plot.ipynb").

Plots folder contains plots.


## Problem 4-5-6

Computes the average energy, average magnetization, specific heat capacity and susceptibility computed with Metropolis algorithm. Saves in text files the average energy and magnetization. Compares the numerical values with the analytical ones. Results of problem 5 and histograms from problem 6 are plotted using "prob_5_plot.ipynb".

Compile: g++ main_5.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

Run: ./main.exe L T "keyword"

where L is the size of lattice (integer number), T is the temperature value (double number) and if "keyword" is "ordered" it will give the ordered initial state, while every other word will give the random initial state 


## Problem 8

For a given interval of temperatures values, computes and saves in text files the average energy, average magnetization, specific heat capacity and susceptibility computed with Metropolis algorithm. Results against temperatures are plotted using "prob_8_plot.ipynb".

Compile: g++ main_8.cpp utils.cpp -std=c++11 -o main.exe -larmadillo 

Run: ./main.exe L "keyword"

where L is the size of lattice (integer number) and if "keyword" is "ordered" it will give the ordered initial state, while every other word will give the random initial state


