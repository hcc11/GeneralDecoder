% simulations of the three-layer network for Att. & Unatt. states. 
% mex EIF1DRFfastslowSyn.c
% mex spktime2count.c

% run following to generate orientation preference map before simulation: 
% Nx=50;
% Kc=5*2*pi; % spatial frequency
% theta_map=ori_map(Nx,Kc);
% save([data_folder 'theta_map'],'theta_map');

% all updated parameters need to be saved in struct 'ParamChange' 

clear 
%%%%%%%%%%%%% this part is for running as a job array on cluster %%%%%%%%%%%
rng('shuffle');
AI = getenv('SLURM_ARRAY_TASK_ID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_folder='';

%%%%%% for general decoder, 50 orientations %%%%%%%%%%
%%% parameters for different attentional modulation levels 
inE_range=[0  0  0  0   0];
inI_range=[0 0.1  0.2  .4 .5];
Jx_range=[1;0.4]*[20 20.5 21 22 23]; 
testp.theta0=0.02:.02:1;  % theta is normalized between 0 and 1 
Nth=length(testp.theta0); 
task = 'multith'; 

%%%%% for Fisher information of the specific decoder, two orientations  %%%%%%%%
%%%% parameters for different attentional modulation levels 
% inE_range=[0  0 ];
% inI_range=[0 .5];
% Jx_range=[1;0.4]*[20 23]; 
% dtheta=0.01;
% testp.theta0=[.5-dtheta/2, .5+dtheta/2]; % range of orientation theta 
% Nth=length(testp.theta0);  
% task = '2th'; 

%%%%%%%%  job index %%%%%%%%%%%%
Ntrial=100; % per param set
Np=length(inI_range);  % number of attentional modulation levels
trial=mod(job_dex-1,Ntrial)+1+Ntrial*(ceil(job_dex/Ntrial/Np)-1), % trial number 
pid=mod(ceil(job_dex/Ntrial)-1,Np)+1,  % index for attentional modulation level 
% % % run pid =1 first.  
% % % The simulation runs multiple attentional states, with same spike trains
% % % (s1) from layer 2, saved in files with pid = 1.  
% % % Different attentional states corresponds to different param(2).Iapp and param(2).Jx.   
testp.Iapp=[inE_range(pid);inI_range(pid)];
Jx = Jx_range(:,pid); 


T=20000; % total simulation time (ms) 
sigma_n=3.5; % variance of input noise
tau_n=40; % time constant of input noise 

Prx=[.05; .05]; % feedforward connection probability 
dt=.01;  % time step size (ms) 

opt.save=1; % save data
if pid ==1
    opt.saveS1=1;  % save spike times from Layer 1 (default is 1, save to filename, or s1_fname if specified)
    opt.loadS1=0;
else
    opt.saveS1=0;
    opt.loadS1=1;
    s1_fname = strrep(sprintf('%sRF2D3layer_%s_sigma_n%.03g_Jex%.03g_Jix%.03g_ID%.0f',...
        data_folder,task,sigma_n,Jx_range(1,1),Jx_range(2,1),trial),'.','d'),
end
opt.saveS2=0; % don't save spike times from Layer 2 (default is 1, save to filename. For multiple parameters (size(Jx,2)>1), s2 is saved in struct spk )
opt.CompCorr=0; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 2 & 1  
opt.plotPopR=0; % plot population rate
opt.fixW=1;  % generate weight matrices with specified random seeds 
    Wseed1=randi(1e4);
    Wseed2=randi(1e4);
opt.saveW=0;  % save weight matrices in W_fname
opt.useWfile=1;  % use weight matrices saved in W_fname  
    W_fname=[data_folder 'weight_sigRX2_0d05.mat'];

if trial==1
    opt.saveParam=1;  % save parameters
else
    opt.saveParam=0; 
end

filename=strrep(sprintf('%sRF2D3layer_%s_sigma_n%.03g_Jex%.03g_Jix%.03g_ID%.0f',...
    data_folder,task,sigma_n,Jx(1),Jx(2),trial),'.','d'),

%%%%%%%% define firing rate of L4 neurons %%%%%%%%
p_stim.stim_type='OriMap_gabor_Tseg';
rX=.01; % mean input rate (kHz) 
Nx=50; 
Kc=5*2*pi; % spatial frequency 
% theta_map=ori_map(Nx,Kc); 
load([data_folder 'theta_map'],'theta_map');
dx=0.04;
x=-.48:dx:.48;
[X,Y]=meshgrid(x,x);
X=X(:);Y=Y(:); 
sigma=0.2;
lambda=.6;
Imag=@(theta) (exp(-(X.^2+Y.^2)/(2*sigma^2))*ones(size(theta)))...
    .*cos(2*pi/lambda*(X*cos(theta)+Y*(sin(theta))));  % Gabor image (size 625x1) 

Filter=@(theta) exp(-(X.^2+Y.^2)/(2*sigma^2))*ones(size(theta))...
    .*cos(2*pi/lambda*(X*cos(theta)+Y*(sin(theta))))...
    /(sum(Imag(0).^2));  % Gabor filters for each L4 neuron (size(Nx^2 x NI) 
F=Filter(theta_map(:)'*pi)'; % Gabor filters for each input neuron
fr=F*Imag(mean(testp.theta0)*pi);
F=F/mean(fr)*rX;

fr=F*Imag(testp.theta0*pi);  %  firing rate for each theta  (size Nx^2 x Nth) 
NI=numel(Imag(0));

p_stim.F=F; % firing rate (Nx by Nx) 
p_stim.theta_map=theta_map;
p_stim.NI=NI;
p_stim.rX=rX;
p_stim.fr=fr;
p_stim.sigma_n=sigma_n;
p_stim.sigma=sigma;
p_stim.lambda=lambda;
p_stim.T=T;
p_stim.tau_n=tau_n; 
p_stim.T_on=200; % stim. is ON for 200 ms, then OFF for 300 ms  
p_stim.T_off=300; 
p_stim.rX_off=.005; % firing rate during OFF intervals (kHz)
Nseg=ceil(T/(p_stim.T_on+p_stim.T_off)); % number of stim. presentations per trial 
if opt.loadS1
    load(s1_fname,'th_id')
    p_stim.th_id=th_id;
else
    p_stim.th_id=int8(randsample(Nth,Nseg,1));  % randomly select orientation id for each stim. presentation 
end
p_stim.theta0=testp.theta0;

clear F Filter theta_map X Y Imag;

ParamChange={'param(2).Iapp', testp.Iapp; 'filename', filename;...
    'dt', dt; 'T', T; 'Nc',Nc;'param(2).Jx',Jx; 'param(2).Prx',Prx; 'p_stim',p_stim};

ParamChange=[ParamChange;{'param(2).sigmaRX',.05*ones(2,1)}];    

if exist('s1_fname','var')
    ParamChange=[ParamChange;{'s1_fname',s1_fname}];
end

if opt.useWfile==1
    if exist(W_fname, 'file')==0  % if W_fname does not exist, save W in simulation 
        opt.saveW=1;
        opt.useWfile=0; 
    end
end

if opt.fixW && opt.useWfile==0
    ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
end

if opt.saveW==1||opt.useWfile==1
    ParamChange=[ParamChange;{'W_fname',W_fname}];
end
   

RF2D3layer(opt, ParamChange) % main simulation 

th_id=p_stim.th_id;
save(filename,'Jx','testp','th_id','W_fname','-append')
