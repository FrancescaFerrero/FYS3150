# Project 5

Codes for problem 2, 3, 4, 5, 6, 7, 8, 9.

Codes for the plots in problem 7, 8 and 9 are written in Phython, on Jupyter Notebook ("prob_7_plot.ipynb", "prob_8_plot.ipynb" and "prob_9_plot.ipynb").

Plots folder contains plots.
Animations folder contain animations of the 2D probability.


## main_7.cpp

Computes the total probability at each time step, with and without the double-slit barrier. It saves it together with the time in a txt file.

Compile and linking: g++ main_7.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

Run: ./main.exe h dt T x_c sigma_x p_x y_c sigma_y p_y v_0
where h and dt are the spatial and temporal step sizes, T is the total time, x_c and y_c are the coordinates of the centre of the initial wave packet, sigma_x and sigma_y  are the initial widths of the wave packet in the x and y directions, p_x and p_y are the wave packet momenta, v_0 constant potential inside the barrier. 

## main_8.cpp

Computes and saves in a binary file the probability for every time step of the simulation and the wavefunction for three different time steps: $t=0.000$, $t=0.001$ and $t=0.002$.

Compile and linking: g++ main_8.cpp utils.cpp -std=c++11 -o main.exe -larmadillo

Run: ./main.exe h dt T x_c sigma_x p_x y_c sigma_y p_y v_0 num_slits
where h and dt are the spatial and temporal step sizes, T is the total time, x_c and y_c are the coordinates of the centre of the initial wave packet, sigma_x and sigma_y  are the initial widths of the wave packet in the x and y directions, p_x and p_y are the wave packet momenta, v_0 constant potential inside the barrier and num_slits is the number of slits considered. 

## prob_7_plot.ipynb

Take data from output of main_7.cpp. Computes the absolute error of the probability. Produces plots of the absolute error as function of time for both the cases with and without double-slit.

## prob_8_plot.ipynb

Take data from output of main_8.cpp and insert edges of those data. Produces colourmap plots for three different time steps: $t=0.000$, $t=0.001$ and $t=0.002$. For the same time steps, also produces colourmap plots of the real and imaginary part of the wavefunction. Creates animations of the simulation for the single, double and triple slit.

## prob_9_plot.ipynb

Take data from output of main_8.cpp. Produces plots of the one-dimensional conditioned probability as a function of the y coordinate for the single, double and triple slit.


