% Script
% Destexhe, Touboul - eNeuro Manuscript, 2021.
% (C) Touboul J. jonathan.touboul@gmail.com.

clear all
close all

PLOT_LEVEL = 1; % Different choices of illustrations 
                    % - PLOT_LEVEL=1 reproduces the heatmap from the paper
                    % for n_bin>1, or the top graphs for n_bin=1;
                    % - PLOT_LEVEL=2 shows the scalings for each choice of parameter
                    % - PLOT_LEVEL=3 shows in addition different graphs
                    % including the Akaike test

N_iter=20;   % Number of simulated networks. 
n_bin=10;   % Number of bins in the heatmap. 
            % For n_bins=1, the program selects the particular choice of 
            % thresholds chosen in the paper.  

fprintf(sprintf('Simulating Brunel model with %d iterations\n',N_iter))
BrunelNetworkSimulation         % Performs simulation of the Brunel network
Analyze_Dependence_Cutoff_2d    % Avalanche analysis

n_bin=1;
Analyze_Dependence_Cutoff_2d

fprintf('Done, press any key\n')


pause()

close all
fprintf(sprintf('Simulating Poisson surrogate with %d iterations\n',N_iter))
IndependentPoisson_CV_Spikes    % Performs simulation of the stochastic surrogate
                                % (same parameters as those of the Brunel
                                % network)
n_bin=10;
Analyze_Dependence_Cutoff_2d    % Analysis 
n_bin=1;
Analyze_Dependence_Cutoff_2d    % Analysis 

