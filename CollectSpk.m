function CollectSpk(task,RunPart)
% RunPart==1; Collect spike counts generated from Sim_Ori_gabor_L2.m
% RunPart==2;  compute mean and variance of spike counts 
data_folder='';

AI = getenv('SLURM_ARRAY_TASK_ID');
job_dex = str2num(AI);

switch task
    case 'multith'  % 50 orientations
        inE=[0  0  0  0   0];
        inI=[0 0.1  0.2  .4 .5];
        Jx=[1;0.4]*[20 20.5 21 22 23];
        Nth = 50; 
        
    case '2th'  % 2 orientations
        inE=[0  0 ];
        inI=[0 .5];
        Jx=[1;0.4]*[20 23];
        Nth = 2; 
end

Pop='E2';  N=4e4;  % E1 ;  E2 ; I1; I2
%  Pop='I2',  N=1e4;
%  Pop='E1',  N=4e4;
%  Pop='I1',  N=1e4;


if RunPart==1
 
    switch Pop
        case {'E2','I2'}
            Np=length(inI);
            pid=mod(job_dex-1,Np)+1,
            num=ceil(job_dex/Np),
        case {'E1','I1'}
            pid=1,
            Np=1;
            num=job_dex,
    end
    
    sigma_n=3.5; 
    Ntshift=0; % shift 1 Tw after stimulus turned on
    
    datafname=@(ID) strrep(sprintf('%sRF2D3layer_%s_sigma_n%.03g_Jex%.03g_Jix%.03g_ID%.0f',...
        data_folder,task,sigma_n,Jx(1,pid),Jx(2,pid),ID),'.','d');
    datafname(1),
    
    fnamesave=strrep(sprintf('%sSpkCounts_%s_sigma_n%.03g_Jex%.03g_muI%.03g_%s_%d',...
        data_folder,task,sigma_n,Jx(1,pid),inI(pid),Pop,num),'.','d'),
    
    %%%%%%%%%%%%%% check files %%%%%%%%%%%%%%%%%%%%%%
    %             IDval=[];
    %             for k=1:Ntrial
    %                 ID=k+(num-1)*Ntrial;
    %                 if exist([datafname(ID) '.mat'], 'file')==0
    %                     sprintf('file %s, does not exist\n',datafname(ID))
    %                 else
    %                     listOfVariables = who('-file', datafname(ID));
    %                     if ismember('th_id', listOfVariables)==0||ismember('E2', listOfVariables)==0
    %                         sprintf('file ID %s, variables not complete\n',datafname(ID))
    %                     else
    %                         IDval=[IDval k];
    %                     end
    %                 end
    %             end
    %             if nnz(IDval)<Ntrial
    %                 error('files not complete')
    %             end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    load(datafname(1),'T','p_stim','Tw','th_id')
    Tw,
    p_stim.T_on,
    p_stim.T_off,
    Nt=floor(T/Tw);
    Nt_on=p_stim.T_on/Tw;
    Nt_off=p_stim.T_off/Tw;
    Nseg=Nt_on+Nt_off;
    idt=1:Nt;
    idt(mod(idt-1,Nseg)+1<=Nt_off)=0;
    
    if Ntshift
        Nstim=length(th_id)-2;
    else
        Nstim=length(th_id)-1;
    end
    Ntrial= 500;
    ns=Nstim*Ntrial;
    
    if Ntshift
        idt=idt(Nt_off+1:end-Nt_on);
    end
    idt=idt(idt>0);
    idt=idt+Ntshift;
    nnz(idt)
    
    X=zeros(N,ns,'int8');
    th_id=zeros(ns,1,'int8');

    for k= 1:Ntrial
        ID=k+(num-1)*Ntrial;
        data=load(datafname(ID));
        if Ntshift
            th_id((1:Nstim)+(k-1)*Nstim)=data.th_id(2:end-1);
        else
            th_id((1:Nstim)+(k-1)*Nstim)=data.th_id(2:end);
        end
        switch Pop
            case 'E2'
                sc=permute(sum(reshape(data.E2{1}(:,idt),N,Nt_on,[]),2),[1 3 2]);
            case 'I2'
                sc=permute(sum(reshape(data.I2{1}(:,idt),N,Nt_on,[]),2),[1 3 2]);
            case 'E1'
                sc=permute(sum(reshape(data.E1(:,idt),N,Nt_on,[]),2),[1 3 2]);
            case 'I1'
                sc=permute(sum(reshape(data.I1(:,idt),N,Nt_on,[]),2),[1 3 2]);
        end
        X(:,(1:Nstim)+(k-1)*Nstim)=sc(:,2:end);
    end
    size(sc)
    
    save(fnamesave,'X','th_id','Ntrial','Ntshift','Tw')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
elseif RunPart==2
    %% compute Fm, Var

    Np=length(inI);
    
    sigma_n=3.5;
    
    for pid=1:Np

        fnamesave=strrep(sprintf('%sSpkCounts_%s_sigma_n%.03g_Jex%.03g_muI%.03g_%s_fm',...
            data_folder,task,sigma_n,Jx(1,pid),inI(pid),Pop),'.','d'),
        datafname=@(ID) strrep(sprintf('%sSpkCounts_%s_sigma_n%.03g_Jex%.03g_muI%.03g_%s_%d',...
        data_folder,task,sigma_n,Jx(1,pid),inI(pid),Pop,ID),'.','d'),

        Nfile=1;
        
        data=load(datafname(1));
        Nstim=size(data.X,2);
        N=size(data.X,1);
        
        ns=Nstim*Nfile;
        Fm=zeros(N,Nth);
        Var=zeros(N,Nth);
        
        X=zeros(N,ns);
        th=zeros(ns,1);
        
        for ID=1:Nfile
            data=load(datafname(ID));
            th((1:Nstim)+(ID-1)*Nstim)=data.th_id;
            X(:,(1:Nstim)+(ID-1)*Nstim)=double(data.X);
        end
        
        for ss=1:Nth
            Fm(:,ss)=mean(X(:,th==ss),2);
            Var(:,ss)=var(X(:,th==ss),[],2);
        end
        
        save(fnamesave,'Fm','Var','datafname','ns')
    end
    
end
