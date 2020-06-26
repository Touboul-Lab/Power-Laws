% Avalanche Analysis (Destexhe, Touboul - PRL Comment, 2020).
% (C) Touboul J. jonathan.touboul@gmail.com.

load Spikes;  % Load here either Brunel or Poisson list of spikes

PLOT_PL=1;      % 1: Plots Power-law
PLOT_Font=1;    % 1: Plots like in Fontenele et al. 


N_Spikes=20;%length(Spikes);    % Number of spike trains considered
n=10;                           % Number of bins in the matrix. 
[Seuil_s_vect,Seuil_t_vect]=meshgrid(linspace(10,40,n),linspace(10,40,n));

% Storage assignment %
R_All_Kept=zeros(n,n,N_Spikes);
A_All_Kept=zeros(n,n,N_Spikes);

L_PL_s=zeros(n,n,N_Spikes);
L_PL_t=zeros(n,n,N_Spikes);
L_LN_s=zeros(n,n,N_Spikes);
L_LN_t=zeros(n,n,N_Spikes);


for u=1:n
    for v=1:n
        Seuil_s=Seuil_s_vect(u,v);
        Seuil_t=Seuil_t_vect(u,v);
        Ratio_All=[];
        Ratio_All_xmin=[];
        A_All=[];
        
        for i=1:N_Spikes %length(Spikes)
    
            allspikes=Spikes{i};
            CV=ComputeCV(allspikes,50);

            FullSize=length(allspikes);
            Sizes=FullSize;
            binsize=1;
    
            for sss=Sizes
                Ssum=allspikes(1:sss);

                Vbin=cumsum(Ssum);  % Number of spikes per bin
                Vbin=Vbin(1:binsize:length(Ssum));
                Vbin=diff(Vbin);

                Silences =find(Vbin==0);
                t= diff(Silences)-1; % Duration of avalanches
                TT=t;
                CumulativeSize=cumsum(Vbin);
                s=diff(CumulativeSize(Silences));
                SS=s;
                s=s(SS>0);
                t=t(SS>0);
                tabsize=tabulate(s);
                tabdur=tabulate(t);
                meansize=zeros(1,size(tabdur,1));
                
                for dur_ind=1:size(tabdur,1)
                    meansize(dur_ind)=mean(s(t==tabdur(dur_ind,1)));
                end
                MM=meansize(~isnan(meansize));
                TT=tabdur(~isnan(meansize));
            end
            
            s=s(s<Seuil_s);
            t=t(t<Seuil_t);
            [alphas,xmin,D,L_PL_s(u,v,i)]=plfitNoXmin(s); % Fits from the smallest avalanche
            [alphat,xmin,D,L_PL_t(u,v,i)]=plfitNoXmin(t);
            
            phats=lognfit(s);
            phatt=lognfit(t);
            L_LN_s(u,v,i)=-lognlike(phats,s);
            L_LN_t(u,v,i)=-lognlike(phatt,t); 
            
            T_no_zeros=tabdur(tabdur(:,2)>0,:);
            S_no_zeros=tabsize(tabsize(:,2)>0,:);


            alphats=polyfit(log(TT),log(MM),1);
            alphats=alphats(1);

            Ratio=(alphat-1)/(alphas-1);
            
            R_All_Kept(u,v,i)=Ratio;
            A_All_Kept(u,v,i)=alphats;
            
            if PLOT_Font
                figure(2*u-1)
                hold on
                plot(CV,Ratio,'rs','MarkerSize',10)
                plot(CV,alphats,'b*','MarkerSize',10)
                title(sprintf('Threshold size:%d, Threshold Duration=%d',Seuil_s,Seuil_t))

                figure(2*u);
            %     hold on
            %     plot(alphas,alphat,'*');
            %     plot(1:0.01:3,1+(1.28)*(0:0.01:2),'linewidth',2)
            %     title(sprintf('Threshold size:%d, Threshold Duration=%d',Seuil_s,Seuil_t))
                hold on
            %     plot(CV,Ratio_xmin,'rs','MarkerSize',10)
                plot(CV,alphats,'b*','MarkerSize',10)
                title(sprintf('With xmin, Threshold size:%d, Threshold Duration=%d',Seuil_s,Seuil_t))
            end
            if PLOT_PL
                  if (u+v==2)
                        figure(100);
                        loglog(tabdur(:,1),tabdur(:,2));
                        hold on
                        figure(101);
                        loglog(tabsize(:,1),tabsize(:,2));
                        hold on
                  end
            end
        end
    end      
end

%%
% Computation of the p-value for Sethna's relationship. 

N_kept=N_Spikes;
figure;
imagesc(mean(R_All_Kept(:,:,1:N_kept),3))
figure
imagesc(mean(A_All_Kept(:,:,1:N_kept),3))

P_val=zeros(n,n);
for u=1:n
    for v=1:n
        [h,P_val(u,v)]=ttest2(R_All_Kept(u,v,1:N_kept),A_All_Kept(u,v,1:N_kept));
    end
end
figure;
P_val(P_val>=0.1)=0.2;
cmap=parula(10)
cmap(1,:)=[0.1 0.1 0.5];


imagesc(linspace(10,40,n),linspace(10,40,n),P_val);
colormap(cmap)
caxis([0.,0.1])
xlabel('Size Cutoff')
ylabel('Duration Cutoff') 

figure;
imagesc(linspace(10,40,n),linspace(10,40,n),mean(R_All_Kept(:,:,1:N_kept),3));
hold on
% colormap(cmap)
% caxis([0.,0.1])
xlabel('Size Cutoff')
ylabel('Duration Cutoff') 
SS=linspace(10,40,n);
[a,b]=find(P_val>0.05)
plot(SS(b)-0.1,SS(a),'*k')
[a,b]=find(P_val>0.01)
plot(SS(b)+0.1,SS(a),'*k')
colorbar()
%% Akaike test

figure;hist(-L_PL_s(:)+L_LN_s(:),25)
figure;hist(-L_PL_t(:)+L_LN_t(:),25)
