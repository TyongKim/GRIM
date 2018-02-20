function output = A1_single_update(input,probdata,Lo,gfundata,system_data)
%Dec-1,2012

%--INPUT--
N_rvs = input.N_rvs;
g_roh = input.g_roh;
clear input;

%--initial parameters--
N = 1000;
mu_init = zeros(1,N_rvs);
mu_cur = mu_init;
sig_init = eye(N_rvs);
sig_cur = sig_init;
mu_plot(1,:) = mu_init;

%generate initial sample
scale_factor = 1;
while 1
    u = mvnrnd(mu_cur,sig_init, N);
    H_x = (lsf_CEAISGM(u,probdata,Lo,N,gfundata,system_data) <= 0);
    if sum(H_x) < 0.01*N && scale_factor <= 2
        %fprintf('\tBad path, restarting...\n');
        %pause(0.05);
        scale_factor = scale_factor + 0.2;
        sig_cur = scale_factor*sig_cur;
    elseif scale_factor >= 2 || sum(H_x) > 0
        break;
    else
        break;
    end
end

%single guassian update using CE
t = 1; %iteration
while 1
    %fprintf('\tIteration of single update: %i\n',t);
    %pause(0.1);
    f_x = mvnpdf(u,mu_init,sig_init);
    h_x = mvnpdf(u,mu_cur,sig_cur);
    W_xwu = f_x./h_x;
    HW_xwu = H_x.*W_xwu;
    clear W_xwu f_x h_x;
    
    %compute updated mean
    mu_temp = ((HW_xwu)'*u)./sum(HW_xwu);
    
    %Version .02 updated protion
    if isnan(mu_temp)
        mu_temp = mu_init;
    end
    
    if mu_temp ~= mu_init
        %compute updated sigma
        sig_up_sum = zeros(N_rvs,N_rvs);
        for iii = 1:N
            sig_up_temp = HW_xwu(iii,1).*(u(iii,:)-mu_cur)'*(u(iii,:)-mu_cur);
            sig_up_sum = sig_up_sum + sig_up_temp;
            clear sig_up_temp
        end
        sig_temp = sig_up_sum./sum(HW_xwu);
        sig_temp = round(sig_temp.*1e8)./1e8;
        clear sig_up_sum u iii;
    else
        sig_temp = sig_init;
    end
        
    %%CHECK VALID COV MATRIX FOR UPDATING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    %postive definitive check
     check_PD = mu_temp*sig_temp*mu_temp';
    
    %symmety check
    sym_check_t = sum(sum(sig_temp' == sig_temp));
    sym_check = N_rvs^2;
    
    %check if stdev is too small
    for i_SD_check = 1:N_rvs
        temp_diag(i_SD_check) = sig_temp(i_SD_check,i_SD_check);
    end
    SD_check = min(temp_diag);
    clear i_SD_check temp_diag
            
    if (SD_check >= 0.5 && SD_check <= 15 && check_PD > 0 && sym_check_t == sym_check...
            && det(sig_temp)>0.2)
        %fprintf('\t\t- covariance matrix value okay\n');
        mu_cur = mu_temp;
        sig_cur = sig_temp;
    else
        if (SD_check < 0.5 | sig_temp == NaN(N_rvs))
            %fprintf('\t\t- STD too small\n');
            mu_cur = mu_temp;
            sig_cur = scale_factor*sig_init;
            H_x = 0;
        elseif (SD_check >15 )
            %fprintf('\t\t- STD too big\n');
            mu_cur = mu_temp;
            sig_cur = scale_factor*sig_init;
            H_x = 0;
        elseif (check_PD <= 0)
            %fprintf('\t\t- covariance matrix not positive difinitive\n')
            mu_cur = mu_temp;
            sig_cur = scale_factor*sig_init;
            H_x = 0;
        elseif (det(sig_temp)<0.2)
            %fprintf('\t\t- det(covariance matrix) too small\n')
            mu_cur = mu_temp;
            sig_cur = scale_factor*sig_init;
            H_x = 0;
        elseif (sym_check_t ~= sym_check)
            %fprintf('\t\t- covariance matrix not symmetric\n')
            mu_cur = mu_temp;
            sig_cur = scale_factor*sig_init;
            H_x = 0;
        end
    end
    %%end checking%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear mu_temp sig_temp;
    
    %store mean updating path to plot
    mu_plot(t+1,:) = mu_cur;
    
    %generate new sample w/ updated mean and sigma
    u = mvnrnd(mu_cur, sig_cur, N);
    H_x = lsf_CEAISGM(u,probdata,Lo,N,gfundata,system_data) <= 0;
       
    %check convergence w/ g_roh
    if (sum(H_x)/N) >= g_roh %&& t >= 3
        break;
    else
        %one more round of update
        t = t+1;
    end
end

%--OUTPUT--
%fprintf('\tSamples in failure domain, g-quantile = %2.2f%%\n\n', sum(H_x)*100/N);
output.mu_cur = mu_cur;
output.sf = scale_factor;