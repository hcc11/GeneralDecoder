% To compute the linear Fisher information (FI) for specific decoders (Fig. 3B).  
% Modified codes from 
% I. Kanitscheider, R. Coen-Cagli, A. Kohn, A. Pouget, 
% Measuring fisher information accurately in correlated neural populations. 
% PLoS computational biology 11, e1004218 (2015).

% outputs: 
% FI_BC: the bias-corrected Fisher information 
% FIVAL: the Fisher information on the Validation set using early
% stopping method 
% FITR: the Fisher information on the Training set using early
% stopping method  
% (FI_BC is typically between FITR and FIVAL, and converge as # of trials
% increase.  Can set NR=0 to not run early stopping method, which can be slow for large N. )
% Nid: sampled neuron index for computing FI  
% compute and save correlation vs distance when N=1600

% need to define: 
% datafname: data filename for spike count matrix (#neurons x #trials) 
% fnamesave: filename to save computed FI 
% data_folder: folder name to read spike count data and to save FI 


data_folder=''; % folder name to read spike count data and to save FI 

rng('shuffle');
AI = getenv('SLURM_ARRAY_TASK_ID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset); 

Types={'2th_sigma_n3d5_Jex20_muI0_E2','2th_sigma_n3d5_Jex23_muI0d5_E2'};

ds=0.01;
sigma_n=3.5;
Np=length(Types),
task='N'; 
% task='Ntr'; 
switch task
    case 'N'
        Nn=9; % # of N to compute
        Nrep=10; % # of repetition of neuron sampling   
        ntype=ceil(job_dex/Nrep/Nn); % index for parameter set in 'Types'
        Nrun=mod(job_dex-1,Nrep)+1; % index for the repetition number  
        ipN= mod(ceil(job_dex/Nrep)-1,Nn)+1; % index for N in 'N_range'
    case 'Ntr'
        Nn=10; % # of Ntr to compute 
        Nrep=10; % # of repetition of neuron sampling 
        ntype=ceil(job_dex/Nrep/Nn); % index for parameter set in 'Types'
        Nrun=mod(job_dex-1,Nrep)+1; % index for the repetition number 
        ipN=mod(ceil(job_dex/Nrep)-1,Nn)+1; % index for Ntr in 'Ntr_range'
end

Type=Types{ntype}, 

datafname=@(ID) sprintf('%sSpkCounts_%s_%d',data_folder,Type,ID); % data filename for spike count matrix (#neurons x #trials)

Nfile= 3;  
data=load(datafname(1));
Nstim=size(data.X,2);
ns=Nstim*Nfile;

%%%%%%%%%%%% select neurons with FR > 1 Hz%%%%%%%%%%%%%%
Tw=0.2; % duration of each trial (sec) 
FR_th=1; % threshold (Hz) 
Fm = (mean(data.X(:,data.th_id==1),2)+mean(data.X(:,data.th_id==2),2))/2; % average spike count per neuron
ind_FR=find(Fm>FR_th*Tw);
sprintf('number of neurons w/ rate larger than %d Hz: ',FR_th, nnz(ind_FR)) 
%%%%%%%%% all neurons %%%%%%%%%%%%%%%%%%
% ind_FR=(1:4e4)';

switch task
    case 'N'
        %%%%%%%%%%%% vary N %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_range=[50 100 200 400 800 1600 3200 6400  12800 25600 4e4];
        N=N_range(ipN);
        fnamesave=strrep(sprintf('%sFI_%s_N%d_%d',data_folder,Type,N,Nrun),'.','d'),
        Ntr=ns/2;
        if N==4e4
            Nid=ind_FR;
        else
            Nid=randsample(ind_FR, N);
        end
        if Nrun==1 && N<5e3
            NR=5;  % how many cross-validation splits for early stopping
        else
            NR=0;
        end
    case 'Ntr'
        %%%%%%%%%%%% vary Ntr %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N=6400;
        Ntr_range=round(exp(linspace(log(N*2),log(ns/2),10)));    
        Ntr=Ntr_range(ipN),
        fnamesave=strrep(sprintf('%sFI_%s_N%d_Ntr%d_%d',data_folder,Type,N,Ntr,Nrun),'.','d'),
        if ipN==1
            Nid=randsample(ind_FR, N);
        else
            load(strrep(sprintf('%sFI_%s_N%d_Ntr%d_%d',data_folder,Type,N,Ntr_range(1),Nrun),'.','d'),'Nid')
        end 
        NR=0; 
end


if exist([fnamesave '.mat'], 'file')
    delete([fnamesave '.mat'])
end

COV=FIdecoder(datafname,fnamesave,Nid,NR,Nfile,Ntr,ds);

if N==1600
    U=triu(ones(size(COV)),1);
    covm=mean(COV(U==1))
    R=corrcov(COV);
    corr=mean(R(U==1))
    
    Ne1=200;
    Ix2=(ceil(Nid/Ne1))/Ne1;
    Iy2=(mod((Nid-1),Ne1)+1)/Ne1;
    D = pdist2([Ix2,Iy2],[Ix2,Iy2],'euclidean');
    D=D(U==1);
    R=R(U==1);
    dmax=0.5;
    dd=0.025;
    daxis=0:dd:dmax;
    [n, ind] = histc(D,daxis);
    n=n(1:end-1); % discard last bin (d>dmax)
    Cd=zeros(length(n),1);
    cov_m=zeros(length(n),1);
    for k=1:length(n)
        Cd(k)=mean(R(ind==k));
    end
    daxis=daxis(1:end-1)+dd/2;
    save(fnamesave,'covm','corr','Cd','daxis','-append')
end
