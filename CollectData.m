function CollectData(task)
    % collect Fisher information data (Fig. 3B)
    % run after FIdecoder_cluster_L1.m or GD_th_cluster.m
    switch task
        case 'specific'
            % run after FIdecoder_cluster_L1.m
            %%%%%%% FI vs. N %%%%%%%%%%
            data_folder='';
            
            Types={'2th_sigma_n3d5_Jex20_muI0_E2','2th_sigma_n3d5_Jex23_muI0d5_E2'};
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            sigma_n=3.5;
            
            Np=length(Types);
            Nn=9; % # of N
            Nrun=10;
            NR=5;
            fname=@(Type,N,run) strrep(sprintf('%sFI_%s_N%d_%d',data_folder,Type,N,Nrun),'.','d'),
            
            N_range=[50 100 200 400 800 1600 3200  6400  12800];
            
            fnamesave=strrep(sprintf('%sFI_%s_sigma_n%.03g_Nsum',data_folder,task,sigma_n),'.','d'),
            
            FITR=NaN(Nn,NR,Nrun,Np);
            FIVAL=NaN(Nn,NR,Nrun,Np);
            FI_w=NaN(Nn,NR,Nrun,Np);
            FIBC=NaN(Nn,Nrun,Np);
            FInaive=NaN(Nn,Nrun,Np);
            for pid=1:Np
                Type=Types{pid},
                for run=1:Nrun
                    for ipN=1:Nn
                        N=N_range(ipN);
                        if exist([fname(Type,N,run) '.mat'], 'file')
                            data=load(fname(Type,N,run));
                            listOfVariables = who('-file', fname(Type,N,run));
                            if ismember('FITR0', listOfVariables)==0
                                % sprintf('file %s, FITR0 does not exist\n',fname(Type,N,run))
                            else
                                FITR(ipN,:,run,pid)=data.FITR0;
                                FIVAL(ipN,:,run,pid)=data.FIVAL0;
                            end
                            if ismember('FI_BC', listOfVariables)==0
                                sprintf('file %s, FI_BC does not exist\n',fname(Type,N,run))
                            else
                                FIBC(ipN,run,pid)=data.FI_BC;
                                FInaive(ipN,run,pid)=data.FInaive;
                            end
                        else
                            sprintf('file %s does not exist\n',fname(Type,N,run))
                        end
                    end
                    if exist([fname(Type,1600,run) '.mat'], 'file')
                        data=load(fname(Type,1600,run),'covm','corr','Cd','daxis','Nfile','Ntr');
                        res(pid).covm(run)=data.covm;
                        res(pid).corr(run)=data.corr;
                        res(pid).Cd(:,run)=data.Cd;
                        res(pid).Nfile=data.Nfile;
                        res(pid).Ntr=data.Ntr;
                        res(pid).daxis=data.daxis;
                    end
                end
            end
            
            save(fnamesave,'FITR','FIVAL','FIBC','N_range','res','Types')
            
        case 'general'
            % run after GD_th_cluster.m
            data_folder='';
            Types={'multith_sigma_n3d5_Jex20_muI0_E2','multith_sigma_n3d5_Jex23_muI0d5_E2'};
            
            Np=length(Types);
            Type=Types{1};
            
            sigma_n=3.5;
            Nn=9;
            Nrun=10;
            NR=5;
            
            N_range=[50 100 200 400 800 1600 3200  6400  12800];
            
            fname=@(Type,N,run) strrep(sprintf('%sFI_%s_N%d_%d',data_folder,Type,N,Nrun),'.','d'),
            
            fnamesave=strrep(sprintf('%sFI_%s_sigma_n%.03g_Nsum',data_folder,task,sigma_n),'.','d'),

            FITR=NaN(Nn,NR,Nrun,Np);
            FIVAL=NaN(Nn,NR,Nrun,Np);
            mseTR=NaN(Nn,NR,Nrun,Np);
            mseVAL=NaN(Nn,NR,Nrun,Np);
            mse0=NaN(Nn,Nrun,Np);
            FI0=NaN(Nn,Nrun,Np);
            FIBC=NaN(Nn,Nrun,Np);
            Iters=NaN(Nn,NR,Nrun,Np);
            
            for run=1:Nrun
                for pid=1:Np
                    for ipN=1:Nn
                        N=N_range(ipN);
                        if exist([fname(N,pid,run) '.mat'], 'file')
                            data=load(fname(N,pid,run));
                            listOfVariables = who('-file', fname(N,pid,run));
                            if ismember('FITR', listOfVariables)
                                FITR(ipN,:,run,pid)=data.FITR;
                                FIVAL(ipN,:,run,pid)=data.FIVAL;
                                mseTR(ipN,:,run,pid)=data.mseTR;
                                mseVAL(ipN,:,run,pid)=data.mseVAL;
                                Iters(ipN,:,run,pid)=data.Iters;
                            end
                            if ismember('E', listOfVariables)
                                mse0(ipN,run,pid)=1-(data.E(1,1)+data.E(2,2));
                            end
                            if ismember('FI0', listOfVariables)
                                FI0(ipN,run,pid)=data.FI0;
                            end
                            if ismember('FI_BC', listOfVariables)==0
                                sprintf('file %s, FI_BC does not exist\n',fname(Type,N,run))
                            else
                                FIBC(ipN,run,pid)=data.FI_BC;
                            end
                            
                        else
                            sprintf('%s does not exist',fname(N,pid,run))
                        end
                    end
                end
            end
            
            save(fnamesave,'FIBC','FI0','FITR','FIVAL','mseTR','mseVAL','mse0','Iters','N_range','Types')
            
            
    end
    
    
    
    
