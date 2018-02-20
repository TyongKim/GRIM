%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Manual of Cross-Entropy-Based Adaptive Importance Sampling Using      %
%   Gaussian Mixtures (CE-AIS-GM) Algorithm by Nolan Kurtz and Junho Song,%
%   (2013) Department of Civil and Environmental Engineering,             %
%   Univerisity of Illinois at Urbana-Champaign, Urbana, IL 61820, USA.   %
%   CEAISGM version 1.0 April-01, 2013                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;

fprintf(' ****************************************** \n');
fprintf('|     Welcome to CEAISGM Version 1.0       |\n');
fprintf('| (Cross-Entropy-Based Adaptive Impontance |\n');
fprintf('| Sampling) This algorithm is developed by |\n');
fprintf('| Nolan Kurtz and Junho Song (2013). For   |\n');
fprintf('| more information, please visit:          |\n');
fprintf('| http://systemreliability.wordpress.com   |\n');
fprintf('|                                          |\n');
fprintf('| WARNING: DO NOT run FERUM with CEAISGM''s |\n');
fprintf('| input file. Error may occur in FERUM.    |\n');
fprintf(' ******************************************\n\n');

if exist('probdata')~=1
    fprintf('WARNING: No input datum is found, program terminates.\n\n');
    return;
end

fprintf('\t- 0. Exit\n\t- 1: CE-AIS-GM Analysis\n\n');

selection = input('PLEASE SELECT OPTION FROM THE ABOVE: ');

switch selection

	case 0 % ----EXIT----
        fprintf('\n=END=\n');
        
    case 1 % ----CEAISGM----
        fprintf('\n')
        disp('************************************************************')
        disp('CEAISGM analysis is running, please wait... (Ctrl+C breaks)')
        disp('Due to imperfection of the code, if the code does not get ')
        disp('to the 5th step, final important sampling, in 20 iterations,')
        disp('please break the code using Ctrl+C and restart it manually.')
        disp('Updated code will be available soon.')
        disp('************************************************************')
        disp(' ')
        
        A0_main_CEAISGM
        
    otherwise 
        disp(' ');
        disp('  You entered an invalid choice.');
end
