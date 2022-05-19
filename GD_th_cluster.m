% Computes Fisher information for the general decoder (Fig. 3B) 

rng('shuffle');
AI = getenv('SLURM_ARRAY_TASK_ID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);

data_folder='';

Types={'multith_sigma_n3d5_Jex20_muI0_E2','multith_sigma_n3d5_Jex23_muI0d5_E2'};

Np=length(Types),
Nfile=6;  

sigma_n=3.5;

testp.theta0=0.02:.02:1;
Nth=length(testp.theta0);
z=@(theta) [cos(2*pi*theta), sin(2*pi*theta)];

task='N'; 
% task='Ntr'; 
%%
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

data=load(datafname(1));
Nstim=size(data.X,2);
ns=Nstim*Nfile;

%%%%%%%%% select from all neurons %%%%%%%%%%%%%%%%%%
ind_FR=(1:4e4)';

switch task
    case 'N'
        %%%%%%%%%%%% vary N %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        N_range=[50 100 200 400 800 1600 3200 6400  12800 25600];
        if strcmp(Type,'X')
            N_range=[50 100 200 400 800 1600 2500];
        end
        N=N_range(ipN);
        Ntr=ns; 
        if Nrun==1 && N<5e3
            NR=5;
        else
            NR=0;
        end
        fnamesave=strrep(sprintf('%sFI_%s_N%d_%d',data_folder,Type,N,Nrun),'.','d'),       
        Nid=randsample(ind_FR, N);
    case 'Ntr'
        %%%%%%%%%%%% vary Ntr %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N=3200;
        NR=0; 
        Ntr_range=round(exp(linspace(log(N*2),log(ns),10)));    
        Ntr=Ntr_range(ipN),
        fnamesave=strrep(sprintf('%sFI_%s_N%d_Ntr%d_%d',data_folder,Type,N,Ntr,Nrun),'.','d'),
        if ipN==1
            Nid=randsample(ind_FR, N);
        else
            load(strrep(sprintf('%sFI_%s_N%d_Ntr%d_%d',data_folder,Type,N,Ntr_range(1),Nrun),'.','d'),'Nid')
        end     
end

if exist([fnamesave '.mat'], 'file')
    delete([fnamesave '.mat'])
end

X=zeros(N,ns);
th=zeros(ns,1);

for ID=1:Nfile
    data=load(datafname(ID));
    th((1:Nstim)+(ID-1)*Nstim)=testp.theta0(data.th_id);
    X(:,(1:Nstim)+(ID-1)*Nstim)=data.X(Nid,:);
end

if Ntr<ns
    idx=randsample(1:ns,Ntr);
    th=th(idx);
    X=X(:,idx);
    ns=Ntr;
end

Fm=mean(X,2);
N=nnz(Fm),
Nid=Nid(Fm>0);

A=X*X';
b=X*z(th);
E=b'*inv(A)*b/(size(th,1)),
mse0=1-(E(1,1)+E(2,2)),
w_opt=A\b;
z_hat=w_opt'*X;
eth0=mod(angle(z_hat(1,:)+1i*z_hat(2,:)),2*pi)/(2*pi)-th';
eth0(eth0>0.5)=eth0(eth0>0.5)-1;
eth0(eth0<-0.5)=eth0(eth0<-0.5)+1;
FI0=1/var(eth0),
FI_BC=(ns-1-N-1)./(ns-1)*FI0, % w/ bias correction 
save(fnamesave,'Nid','Nfile','E','FI0','FI_BC','mse0','z_hat','w_opt')


% train a linear decoder w/ gradient descent & early stopping 
if NR>0
    FITR=NaN(NR,1);
    FIVAL=NaN(NR,1);
    mseTR=NaN(NR,1);
    mseVAL=NaN(NR,1);
    Iters=NaN(NR,1);
    w0=zeros(N,2);
    
    for kk=1:NR
        idx=randperm(ns);
        nsTR=round(ns/3);  % per theta
        nsTE=round(ns/3);
        nsVAL=ns-nsTR-nsTE;
        idxTR=idx(1:nsTR);
        idxTE=idx((nsTR+1):(nsTR+nsTE));
        idxVAL=idx((nsTR+nsTE+1):end);
        
        ATR=X(Fm>0,idxTR)*X(Fm>0,idxTR)';
        bTR=X(Fm>0,idxTR)*z(th(idxTR));
        
        ATE=X(Fm>0,idxTE)*X(Fm>0,idxTE)';
        bTE=X(Fm>0,idxTE)*z(th(idxTE));
        
        w=w0;
        
        maxiters=1e6;
        % gradient descent
        alpha=1/10/max(eig(ATR)),
        iters=0;
        dETEdt=-1;
        tic
        while(dETEdt<0 && iters < maxiters)  % dETEdt<0
            iters=iters+1;
            r = bTR - ATR * w;
            w = w + alpha*r; % update w
            
            dETEdw = ATE*w - bTE;
            dETEdt = r(:,1)'*dETEdw(:,1)+r(:,2)'*dETEdw(:,2);
            
        end
        if(iters==maxiters)
            fprintf('Max iters reached -- run\n'),
        else
            fprintf('GD iter=%d \n',iters),
            dETEdt
        end
        toc
        
        Iters(kk)=iters;
        
        z_hat_test=w'*X(Fm>0,idxVAL);
        ethVAL=mod(angle(z_hat_test(1,:)+1i*z_hat_test(2,:)),2*pi)/(2*pi)-th(idxVAL)';
        ethVAL(ethVAL>0.5)=ethVAL(ethVAL>0.5)-1;
        ethVAL(ethVAL<-0.5)=ethVAL(ethVAL<-0.5)+1;
        FIVAL(kk)=1/var(ethVAL);
        
        e=z_hat_test-z(th(idxVAL))';
        mseVAL(kk) = mean(e(1,:).^2+e(2,:).^2);
        
        z_hat_TR=w'*X(Fm>0,idxTR);
        ethTR=mod(angle(z_hat_TR(1,:)+1i*z_hat_TR(2,:)),2*pi)/(2*pi)-th(idxTR)';
        ethTR(ethTR>0.5)=ethTR(ethTR>0.5)-1;
        ethTR(ethTR<-0.5)=ethTR(ethTR<-0.5)+1;
        FITR(kk)=1/var(ethTR);
        
        e=z_hat_TR-z(th(idxTR))';
        mseTR(kk) = mean(e(1,:).^2+e(2,:).^2);
        save(fnamesave,'w','Nid','FITR','FIVAL','mseVAL','mseTR','Iters','Ntr','-append')
        
    end
    clear ATR ATE;
    mean(FITR),
    mean(FIVAL),
end



