% Compare performance of specific and general decoders with small number of neurons (Fig. 3D)

rng('shuffle');

testp.theta0=0.02:.02:1;
Nth=length(testp.theta0);
z=@(theta) [cos(2*pi*theta), sin(2*pi*theta)];

Nrun=200;  % number of sampling per orientation 

data_folder='';

inE=[0  0  0  0   0];
inI=[0 0.1  0.2  .4 .5];
Jx=[1;0.4]*[20 20.5 21 22 23];
Types={'multith_sigma_n3d5_Jex20_muI0_E2','multith_sigma_n3d5_Jex20d5_muI0d1_E2',...
    'multith_sigma_n3d5_Jex21_muI0d2_E2','multith_sigma_n3d5_Jex22_muI0d4_E2',...
    'multith_sigma_n3d5_Jex23_muI0d5_E2'};

Np= length(Types),
Nfile=6;
Dth=2; % delta theta

N=100; % number of neurons to sample

fnamesave=@(ss,run) strrep(sprintf('%sPc_Jex_muI_th_N%d_ss%d_%d',data_folder,N,ss,run),'.','d');

task='part 1';

switch task
    case 'part 1'
        %% part 1
        
        %%%%%%%%%%%%%%% run on cluster %%%%%%%%%%%%%%
        AI = getenv('SLURM_ARRAY_TASK_ID');
        job_dex = str2num(AI);
        seed_offset = randi(floor(intmax/10));
        rng(job_dex + seed_offset);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        run=mod(job_dex-1,Nrun)+1,
        ss=5*ceil(job_dex/Nrun) -4,  % theta ID
        
        fnamesave(ss,run),
        
        Corr=zeros(1,Np);
        Pc_General=zeros(1,Np);  % P(correct) for general decoder 
        Pc_Specific=zeros(1,Np);  % P(correct) for specific decoder 
        
        datafname=@(pid,ID) sprintf('%sSpkCounts_%s_%d',data_folder,Types{pid},ID); % data filename for spike count matrix 
        data=load(datafname(1,1));
        Nstim=size(data.X,2);
        ns=Nstim*Nfile;
        ind_FR=(1:4e4)';
        if ss==1
            Nid=randsample(ind_FR, N);
        else
            load(fnamesave(1,run),'Nid')
        end
        
        th_id=zeros(ns,1);
        
        for pid=1:Np
            pid
            X=zeros(N,ns);
            for ID=1:Nfile
                data=load(datafname(pid,ID));
                th_id((1:Nstim)+(ID-1)*Nstim)=data.th_id;
                X(:,(1:Nstim)+(ID-1)*Nstim)=data.X(Nid,:);
            end
            clear data
            
            th=testp.theta0(th_id)';
            th_id1=mod(ss-1,Nth)+1;  % first stim.
            th_id2=mod(ss+Dth-1,Nth)+1; % second stim.
            idx1=(th_id<(th_id1+.01)&th_id>(th_id1-.01));  % select trials for stim 1
            idx2=(th_id<(th_id2+.01)&th_id>(th_id2-.01));  % select trials for stim 2
            
            %%%%%%%%%%   general decoder  %%%%%%%%%%
            A=X*X';
            b=X*z(th);
            w_opt=A\b;
            z_hat=w_opt'*X;
            
            a1=z_hat(:,idx1);
            a2=z_hat(:,idx2);
            X0=[a1';a2'];
            Y=[-ones(size(a1,2),1);ones(size(a2,2),1)];
            Mdl = fitcsvm(X0,Y,'CrossVal','on','kfold',2);
            Pc_General(pid)= 1- kfoldLoss(Mdl);
            
            R=corrcoef(X(:,idx1)');
            U=triu(ones(size(R)),1);
            Corr(pid)=mean(R(U==1));
            
            %%%%%%% specific decoder  %%%%%%%%%%%%%%
            a1=X(:,idx1);
            a2=X(:,idx2);
            X0=[a1';a2'];
            Y=[-1*ones(size(a1,2),1); ones(size(a2,2),1)];
            
            Mdl = fitcsvm(X0,Y,'CrossVal','on','kfold',2);
            Pc_Specific(pid)= 1 - kfoldLoss(Mdl);
            
            save(fnamesave(ss,run),'Pc_Specific','Pc_General','Corr','N','Nth','Nid','Dth','ns')
        end
        
    case 'part 2'
        %% collect results
        % run after part 1 is finished,  collect results 

        fnamesave0=strrep(sprintf('%sPc_Jex_inI_th_N%d_sum',data_folder,N),'.','d'),
        Nth=10;
        Pc_General=zeros(Nrun,Nth,Np);
        Pc_Specific=zeros(Nrun,Nth,Np);
        Corr=zeros(Nrun,Nth,Np);
        th_range = 5*(1:Nth)-4;  
        for run=1:Nrun
            for ss=1:Nth 
                data=load(fnamesave(th_range(ss),run));
                
                Pc_General(run,ss,:)=data.Pc_General;
                Pc_Specific(run,ss,:)=data.Pc_Specific;
                Corr(run,ss,:)=data.Corr;

                res(run).Nid=data.Nid;
                res(run).ns=data.ns;
                res(run).Dth=data.Dth;
            end
        end
        
        save(fnamesave0,'Pc_General','Corr','Pc_Specific','res','Jx','inE','inI','Types')
        
end

