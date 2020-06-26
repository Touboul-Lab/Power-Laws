% Simulation of the Brunel JCNS 2000 network with electrode recording
% spikes (Destexhe, Touboul - PRL Comment, 2020).
% (C) Touboul J.
% jonathan.touboul@gmail.com.

Spikes={};
for Iteration=1:200
    
    Duration=5;     % Total duration
    dt=0.1e-3;      % time-step
    N=1000;         % Number of neurons
    
    NTimes=ceil(Duration/dt);   % Number of time-steps
    
    
    %%%%%%%%   PARAMETERS OF THE NETWORK  %%%%%%%%
    
    NE=ceil(0.8*N);      % # excitatory neurons
    NI=N-NE;             % # excitatory neurons
    
    
    
    % (p.185, 2nd column, Brunel 2000)
    V_t=20e-3; % voltage threshold
    V_r=10e-3; % voltage reset
    tau=20e-3; % voltage time constant
    t_rp=2e-3; % refractory period
    
    
    
    %%%%%% BIFURCATION PARAMETERS in Brunel 2000 %%%%%%%
    
    J=0.12e-3;                  % Strength of exc. Connections
    J=J*10000/N;                % Rescaling to make behaviors independent of N
    
    Connectivity=0.1;           % Connectivity coefficients, parameter epsilon in Brunel 2000.
    D=1.8e-3;                        % Transmission Delay
    CE=round(Connectivity*NE);          % Number of exc connections
    CI=round(Connectivity*NI);          % Number of inh connections
    C_ext=CE;
    
    
    nu_thresh=V_t/(J*CE*tau);   % Frequency needed for a neuron to reach the threshold.
    
    
    %%%%%%%%   SIMULATION PARAMETERS %%%%%%%%
    
    N_rp=round(t_rp/dt);        % Refractory (in dt)
    N_del=round(D/dt);          % delay (in dt)
    
    
    g=4+4*(rand());                         % Strength of inh/exc (inh connections -gJ)
    ratioextthresh=1.1+0.2*(rand()-0.5);    % nu_ext/nu_thresh
    nu_ext=ratioextthresh*nu_thresh;        % external Poisson input rate
    
    W=zeros(N,N);
    
    %%%%%%%%   Generation of the connectivity matrix %%%%%%%%
    
    for i=1:N
        ECells=randperm(NE);
        ECells=ECells(1:CE);
        W(i,ECells)=J;
        
        ICells=randperm(NI);
        ICells=NE+ICells(1:CI);
        W(i,ICells)=-1;
    end
    
    
    W(W<0)=-g*J;
    Rasterplot=zeros(N,NTimes);
    Rasterplot_pert=zeros(N,NTimes);
    tic();
    
    V=V_r*ones(N,1);            % Voltage vector
    
    LS=-(N_del+1)*rand(N,1);
    
    allspikes=zeros(1,NTimes);
    
    
    i=0;
    while i<NTimes
        i=i+1;
        if (1+mod(i-1,1e4))==1
            ExInput=J*poissrnd(nu_ext*C_ext*dt,N,1e4);
        end
        
        V=(1-dt/tau)*V+ExInput(:,1+mod(i-1,1e4))+W*Rasterplot(:,max(1,i-N_del));        % Voltage update
        
        spike=V>=V_t;               % Spiking neurons have a "1"
        V(LS>i-N_del)=V_r;          % Refractory period.
        
        
        LS(spike)=i;                % Time of last spike
        V(spike)=V_r;               % Reset membrane potential
        % Store spike times
        Rasterplot(:,i)=spike;      % Each row is (neuron number,spike time)
        allspikes(1,i)=sum(spike);
        
%         progressbar(i,NTimes);
        
    end
    toc
    Spikes{Iteration}=allspikes;
    
end


