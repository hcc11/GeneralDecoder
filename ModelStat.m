% Compute model statistics for different attentional modulation (Fig 2D-F): 
% firing rates, average correlations, and eigenvalues from Factor analysis   
% run CollectSpk('multith',1) first to generat spike counts matrices

rng('shuffle');

data_folder='';
addpath(genpath( 'fa_Yu/'));

inE=[0  0  0  0   0];
inI=[0 0.1  0.2  .4 .5];
Jx=[1;0.4]*[20 20.5 21 22 23];
sigma_n=3.5;

datafname=@(pid,ID) strrep(sprintf('%sSpkCounts_%s_sigma_n%.03g_Jex%.03g_muI%.03g_E2_%d',...
        data_folder,task,sigma_n,Jx(1,pid),inI(pid),ID),'.','d'),
    
Np=length(inI);
Nfile=6;  % number of spike counts files (ID in datafname) 

N=50; % number of neurons to sample each time 
Ne=4e4;
Stim=1; % stimulus orientation index 

M=5; % number of dimensions 
numFolds=2; % number of cross-validation folds for factor analysis 
zDimList=M;

fnamesave=strrep(sprintf('%s3Layer_ModelStat_sigma_n%.03g_Jx_mu_sum',...
    data_folder,sigma_n),'.','d'),


FR_th=1;  % firing rate threshold 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsample=30;
rate=zeros(Nsample,Np);  % firing rate 
Corr=zeros(Nsample,Np);  % average correlation
Cov=zeros(Nsample,Np);   % average covariance 
Lambda=zeros(M,Nsample,Np);  % eigenvalues from factor analysis 

for ss=1:Nsample
    Nid=randsample((1:Ne)', N);  % neuron index 
    
    for pid=1:Np
        % response to Stim
        data=load(datafname1(pid,1));
        Nstim=size(data.X,2);
        th=zeros(Nstim*Nfile,1);
        for ID=1:Nfile
            data=load(datafname1(pid,ID),'th_id');
            th((1:Nstim)+(ID-1)*Nstim)=data.th_id;
        end
        ns=nnz(th==Stim);
        E2=zeros(N,ns);
        tmp1=0; 
        for ID=1:Nfile
            data=load(datafname2(pid,ID));
            tmp2=nnz(data.th_id==Stim);
            E2(:,tmp1+(1:tmp2))=data.X(Nid2,data.th_id==Stim);
            tmp1=tmp1+tmp2;
        end
        rate(ss,pid)=mean(E2(:))/data.Tw*1e3;

        %  Correlation 
        id2=mean(E2,2)>(data.Tw*1e-3*FR_th);
        COV=cov(E2(id2>.5,:)');
        R=corrcov(COV);
        U=triu(ones(size(COV)),1);
        Corr(ss,pid)=mean(R(U==1));
        Cov(ss,pid)=mean(COV(U==1));
        
        % Factor analysis
        
        dim = crossvalidate_fa(E2, 'zDimList', zDimList,'showPlots',false,'numFolds',numFolds);
        L=dim(1).estParams.L;
        LL=L*L';
        [V,D]=eig(LL);
        la=diag(D);
        la=sort(la,'descend');
        Lambda(:,ss,pid)=la(1:M);
        clear dim
    end
end

Tw=data.Tw;  
save(fnamesave,'rate','Corr','Cov','Lambda','Nfile','ns','N','Nsample', ...
    'Stim','datafname','FR_th','Tw')



