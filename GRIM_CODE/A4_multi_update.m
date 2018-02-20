function output = A4_multi_update(input,probdata,Lo,gfundata,system_data)
%update gaussian mixtures algorithm

%--INPUT--
mu_cur = input.mu_init_GM; mu_init_GM = input.mu_init_GM;
N_rvs = input.N_rvs;
N_dps = input.N_dps;
g_roh = input.g_roh;
N = input.N_multi;
scale_factor = input.sf;
clear input;

%--initial parameters--
mu_init = zeros(1,N_rvs);
sig_init = eye(N_rvs);

for j = 1:N_dps
    sig(j).cur = scale_factor.*eye(N_rvs);
end
Pi_cur = ones(N_dps,1)./N_dps;
clear j;
g_init = 0.1; %this is the initial g_quantile

%gemerate initial sample with given latent variable
while 1
    z = mnrnd(N,Pi_cur);         
    zz = cumsum(z);
    for i = 1:N_dps              
        if i ==1,
            idx = 1:zz(i);
        else idx = (zz(i-1)+1):zz(i);
        end
        u(idx,:) = mvnrnd(mu_cur(i,:), sig(i).cur, z(i));
    end
    clear z zz;
    zj =  find(mnrnd(1,Pi_cur));
   % u = mvnrnd(mu_cur(zj,:),sig(zj).cur,N);
    H_x = (lsf_CEAISGM(u,probdata,Lo,N,gfundata,system_data) <= 0);
    if sum(H_x) < g_init*N && scale_factor <= 2.5
        %fprintf('\tBad path, restarting using lager sigma...\n');
        %pause(0.05);
        scale_factor = scale_factor + 0.1;
        sig(zj).cur = scale_factor*sig(zj).cur;
        clear zj;
    elseif sum(H_x) < g_init*N && scale_factor > 2.5 
        %fprintf('\tBad path, restarting...\n');
        %pause(0.05);
        %g_init = g_init - 0.005;
        break
        clear zj;
    else
        %fprintf('\n')
        break;
    end
end

%multi-gussian update algorthim
t = 1; %iteration
H_x = (lsf_CEAISGM(u,probdata,Lo,N,gfundata,system_data) <= 0);
while 1
    %fprintf('\n\tIteration of multi update: %i\n',t);
    pause(0.1);
    f_x = mvnpdf(u,mu_init,sig_init);
    
    %compute h(x) & gamma_bot, since gamma_bot = h(x)
    h_x = 0;
    for j = 1:N_dps
        h_x = h_x + Pi_cur(j,1).*mvnpdf(u,mu_cur(j,:),sig(j).cur);        
    end
    
    W_xwu = f_x./h_x;
    HW_xwu = H_x.*W_xwu;
    gamma_bot = h_x;
    clear j W_xwu f_x h_x;
    
    %this loop is to update parameters for each design point using the same
    %sample, index function and weight function.
    for j = 1:N_dps
        
        gamma_top = Pi_cur(j).*mvnpdf(u,mu_cur(j,:),sig(j).cur);
        gamma_z = gamma_top./gamma_bot;
        clear gamma_top;
        
        %compute new mean
        mu_temp = sum(repmat(HW_xwu.*gamma_z,1,N_rvs).*u)./sum(HW_xwu.*gamma_z);
        
        %%CHECK UPDATED MEAN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if max(isnan(mu_temp))
            mu_temp = mu_cur(j,:);
        end
        %%end checking%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %compute new sigma
        sig_top_sum{j} = zeros(N_rvs,N_rvs);
        for i = 1:N
            sig_top_temp = HW_xwu(i,1)*gamma_z(i)*(u(i,:)-mu_cur(j,:))'*(u(i,:)-mu_cur(j,:));
            sig_top_sum{j} = sig_top_sum{j} + sig_top_temp;
            clear sig_up_temp
        end
        sig_temp{j} = sig_top_sum{j}./sum(HW_xwu.*gamma_z);
        %sig_temp = round(sig_temp.*1e10)./1e10;
        %clear sig_up_sum i;
        
        %compute new Pi
        Pi_cur(j) = (sum(HW_xwu.*gamma_z))./sum(HW_xwu);
        
        %install new mean
        mu_cur(j,:) = mu_temp;
        clear mu_temp;
                    
        sig(j).cur = sig_temp{j};
 
        %%end checking%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear mu_temp sig_temp;
    end
    
    if sum(Pi_cur) ~= 1
        Pi_cur = Pi_cur/sum(Pi_cur);
    end
    
    %check convergence w/ g_roh 
    roh_check = 0;
    for jjj = 1:N_dps
        temp = sum(lsf_CEAISGM(mvnrnd(mu_cur(jjj,:),sig(jjj).cur,N),...
            probdata,Lo,N,gfundata,system_data) <= 0).*Pi_cur(jjj);
        roh_check = roh_check + temp;
    end
    clear temp jjj;

    if (roh_check/N) > 0.85 || t >= 30
        break;
    else
        %one more round of update
        t = t+1;
        clear roh_check;
    end
    
    %generate new sample w/ updated means and sigmas
    z = mnrnd(N,Pi_cur);         
    zz = cumsum(z);
    for i = 1:N_dps              
        if i ==1
            idx = 1:zz(i);
        else idx = (zz(i-1)+1):zz(i);
        end
        u(idx,:) = mvnrnd(mu_cur(i,:), sig(i).cur, z(i));
    end

    H_x = lsf_CEAISGM(u,probdata,Lo,N,gfundata,system_data) <= 0;
end


%--OUTPUT--       
%fprintf('\tSamples in failure domain, g-quantile = %2.2f%%\n\n', roh_check*100/N);
output.mu_cur = mu_cur;
for jjj = 1:N_dps
    output.sig(jjj).cur = sig(jjj).cur;
end
output.Pi_cur = Pi_cur;