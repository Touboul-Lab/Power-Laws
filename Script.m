% Script
% Destexhe, Touboul - PRL Comment, 2020.
% (C) Touboul J. jonathan.touboul@gmail.com.

clear all
close all

N_iter=2;
fprintf(sprintf('Simulating Brunel model with %d iterations\n',N_iter))
BrunelNetworkSimulation
Analyze_Dependence_Cutoff_2d
fprintf('Done, press any key\n')
pause()

close all
fprintf(sprintf('Simulating Poisson surrogate with %d iterations\n',N_iter))
IndependentPoisson_CV_Spikes
Analyze_Dependence_Cutoff_2d

