%Main CEAISGM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main file of Cross Entropy Based Adaptive Importance Sampling Using    %
%% Gaussian Mixtures (CE-AIS-GM) Algorithm by Nolan Kurtz and Junho Song, %
%% Department of Civil and Environmental Engineering,                     %
%% Univerisity of Illinois at Urbana-Champaign, Urbana, IL 61820, USA.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;

%initializing parameters
if strcmp(ceaisgmopt.default,'yes')
    input.N_rvs 	= size(probdata.marg,1);
    input.N_dps 	= 6;
    input.N_multi 	= 10000;
    input.g_roh 	= 0.1;
    input.R_factor 	= 0.5;
    input.cov_level = 0.05;
elseif strcmp(ceaisgmopt.default,'no')
    input.N_rvs 	= size(probdata.marg,1);
    input.N_dps 	= ceaisgmopt.N_search;
    input.N_multi 	= ceaisgmopt.N_presam;
    input.g_roh 	= 0.1;
    input.R_factor 	= ceaisgmopt.R_factor;
    input.cov_level = ceaisgmopt.target_cov;
else
    disp('Incorrect CEAISGM inputs, terminate program.')
    return;
end

%Using Nataf Distribution Model
R_mod = mod_corr(probdata, probdata.correlation);
Lo = chol(R_mod, 'lower');
Lu = chol(R_mod, 'upper');

%--CE algorithm--
disp('Main CE-AIS-GM');
fprintf('*01.Initializing');
fprintf('.'),pause(0.25);
fprintf('.'),pause(0.25);
fprintf('.\n'),pause(0.4);
fprintf('\t- Number of RVs: %i\n\t- Number of search points: %i\n',...
    input.N_rvs, input.N_dps);
fprintf('\t- Multi-gussian algorithm presample size: %i\n\t- Convergence Level of IS: %i%%\n\n',...
    input.N_multi, 100*input.cov_level);

%single gaussian algorithm
fprintf('*02.Running single gaussian algorithm');
fprintf('.'),pause(0.25);
fprintf('.'),pause(0.25);
fprintf('.\n'),pause(0.25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:2
    output = A1_single_update(input,probdata,Lo,gfundata,system_data);
    temp_sf(i,1) = output.sf;
    temp(i).mu_SG = output.mu_cur;
    temp(i).R_factor = 1;
    RRR(i,1) = A2_radius(temp(i));
end
[aaa,bbb] = sort(RRR,'descend');
input.mu_SG = temp(bbb(2)).mu_SG;
input.sf = temp_sf(bbb(2));
clear temp temp_sf RRR output aaa bbb i;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute scaled radius & initial point for GM
fprintf('*03.Computing initial points for gaussian mixtures');
fprintf('.'),pause(0.25);
fprintf('.'),pause(0.25);
fprintf('.\n'),pause(0.4);
input.radius = A2_radius(input);
input.mu_init_GM = A3_GM_dps(input);
fprintf('\t- %i%% of radius obtained from step 02\n\t  is used to compute initial points.\n\n',...
    input.R_factor*100);

fprintf('*04.Running multi-gaussian algorithm');
%multi-guassian algorithm
fprintf('.'),pause(0.25);
fprintf('.'),pause(0.25);
fprintf('.\n'),pause(0.25);
[output] = A4_multi_update(input,probdata,Lo,gfundata,system_data);
input.mu_final = output.mu_cur;
input.sig_final = output.sig;
input.Pi_final = output.Pi_cur;
clear output;

%final important sampling
fprintf('*05.Final important sampling');
fprintf('.'),pause(0.25);
fprintf('.'),pause(0.25);
fprintf('.\n'),pause(0.25);
[output] = B1_final_IS(input,probdata,Lo,gfundata,system_data);

fprintf('\n=END OF ANALYSIS=\n\n');

clear input outxls OUT R_mod Lu

N_presam=ceaisgmopt.N_presam;
sensitivity=C1_variance_decompose(Lo,output,N_presam,probdata,gfundata,system_data);
