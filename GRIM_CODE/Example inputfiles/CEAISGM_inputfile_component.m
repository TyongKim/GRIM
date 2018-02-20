%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Inputfile for Cross Entropy Based Adaptive Importance Sampling Using  %
%   Gaussian Mixtures (CE-AIS-GM) Algorithm by Nolan Kurtz and Junho Song,%
%   Department of Civil and Environmental Engineering,                    %
%   Univerisity of Illinois at Urbana-Champaign, Urbana, IL 61820, USA.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
%Marginal Probability Data For Random Variables
%probdata.marg(1,:) =  [ type  mean  stdev   0 0 0 0 0 0]

probdata.marg(1,:) =  [ 1    0   1   0     0 0 0 0 0];
probdata.marg(2,:) =  [ 1    0   1   0     0 0 0 0 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    - FERUM distributions library -                      %
% Type: 1 = Normal distribution                                           %                                         
%       2 = Lognormal distribution                                        %                                         
%       3 = Gamma                                                         %
%       4 = Shifted Exponential marginal distribution                     %
%       5 = Shifted Rayleigh marginal distribution                        %
%       6 = Uniform distribution                                          %
%       7 = Beta                                                          %
%       8 = Chi-square                                                    %
%                                                                         %
%      11 = Type I Largest Value marginal distribution                    % 
%      12 = Type I Smallest Value marginal distribution                   % 
%      13 = Type II Largest Value marginal distribution                   %
%      14 = Type III Smallest Value marginal distribution                 %
%      15 = Gumbel (same as type I largest value)                         %
%      16 = Weibull marginal distribution (same as Type III Smallest      % 
%           Value marginal distribution with epsilon = 0 )                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Correlation Matrix of Random Variables
probdata.correlation =  eye(2);


probdata.parameter = distribution_parameter(probdata.marg);


%CE-AIS-GM Analysis Option Values
ceaisgmopt.default  = 'no';
ceaisgmopt.N_search = 2;
ceaisgmopt.N_presam = 10000;
ceaisgmopt.R_factor = 0.5;
ceaisgmopt.target_cov = 0.05;


%Limit State Function Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For simulation efficiecy in Matlab, please input randon variables     %
%   as column vectors, for example x(:,1), x(:,2),..., x(:,n), and use    %
%   (.*), (./), (.^) in stead.                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gfundata(1).expression = '5 - x(:,2) - 0.5.*(x(:,1) - 0.1).^2';



%System Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - System Type: series, parallel and general.                          %
%   The system type should be defined according to the corresponding      %
%   system. For general system, system information shall be provided      %
%   using the vector system_data.general, where [1 3 0 2 3 0 1 4]         %
%   repersents event system of [(E1,E3)U(E2,E3)U(E1,E4)]. Each cut-sets   %
%   shall be separate by a "0". If the input system is not general, input %
%   dummy values for the vector or jsut remove it.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system_data.type    = 'no';      %(series/parallel/general/no)
%system_data.general = [1 3 0 2 0 1];


