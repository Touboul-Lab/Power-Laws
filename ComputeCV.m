function CV=ComputeCV(allspikes,deltat)
% Computing unbiased CV - For Destexhe Touboul PRL 2020 Commentary.
% (c) J. Touboul. jonathan.touboul@gmail.com

    N=length(allspikes);
    Nbins=floor(N/deltat);
    binnedSpikes = squeeze(sum(reshape(allspikes(1:Nbins*deltat),deltat,[]),1))/deltat;

%     CV=std(binnedSpikes)./mean(binnedSpikes); % Biased estimator. Used in
%     Fontenele et al PRL?

CV=(1+1/(4*length(binnedSpikes)))*std(binnedSpikes)/mean(binnedSpikes);
end