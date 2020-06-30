% Simulation of the surrogate Poisson Network
% (Destexhe, Touboul - PRL Comment, 2020).
% (C) Touboul J.
% jonathan.touboul@gmail.com.

Nn=1000;                % Number of "neurons"

Spikes={};

if ~exist('N_iter')
    N_iter=1;           % N_iter is the number of independent simulations
end

% Numerical parameters
dt=1e-2;
Tf=1000;
N=round(Tf/dt);

for mc=1:N_iter
    fprintf(sprintf('Iteration %d/%d\n',mc,N_iter))
    % Parameters of the Ornstein-Uhlenbeck rate vary randomly
    theta=2+1*rand();           % Timescale parameter
    sigma=0.7+0.8*rand();       % Diffusion coefficient
    
    mu_r=0e-4;
    
    
    % B=randn(N,1);
    
    X=zeros(1,N);       % Stores the rate.
    allspikes=zeros(1,N);     % Stores spikes.
    X0=1;               % Rate initial condition
    B=randn(N,1);       % White noise (pre-computed)

    for i=1:N
        if i==1
            X(i)=X0+sigma*sqrt(dt)*B(i);
        else
            X(i)=(1-theta*dt)*(X(i-1))+sigma*sqrt(dt)*B(i);
        end
       
        allspikes(i)=(poissrnd(Nn*max(X(i),0)));
    end
    Spikes{mc}=allspikes;
    
end