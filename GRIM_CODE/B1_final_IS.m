function output = B1_final_IS(input,probdata,Lo,gfundata,system_data)
%final inportant sampling using given parameters

%--INPUTS--
N_rvs = input.N_rvs;
N_dps = input.N_dps;
mu_cur = input.mu_final;
sig = input.sig_final;
Pi = input.Pi_final;
cov_level = input.cov_level;
clear input;

%--initializing parameters--
mu_init = zeros(1,N_rvs);
sig_init = eye(N_rvs);
N = 10000;
n = 0;
convergence_check = 1;

%start sampling
while convergence_check > cov_level || n <= 30
    n = n + 1;
    
    %show proccess bar
    if n == 1000
        %fprintf('\t-Proccessing...(use Ctrl+C to break)\n\t');
    end
    
    if n >= 1000 && n/200 == ceil(n/200)
        %fprintf('> ');
        if n/4000 == ceil(n/4000)
            %fprintf('\n\t');
        end
    end
    
    %choose a gaussian density to generate samples
    z = mnrnd(N,Pi);         
    zz = cumsum(z);
    for i = 1:N_dps              
        if i ==1,
            idx = 1:zz(i);
        else idx = (zz(i-1)+1):zz(i);
        end
        u(idx,:) = mvnrnd(mu_cur(i,:), sig(i).cur, z(i));
    end
    
    H_x = (lsf_CEAISGM(u,probdata,Lo,N,gfundata,system_data) <= 0);
    f_x = mvnpdf(u,mu_init,sig_init);
        
    %compute updated pdf, h(x)
    h_x = 0;
    for j = 1:N_dps;
        h_temp = Pi(j,1).*mvnpdf(u,mu_cur(j,:),sig(j).cur);
        h_x = h_x + h_temp;
        clear h_temp;
    end
    
    %comput the failure probability of each run
    Pfn(n,1) = sum(H_x.*f_x./h_x)/N;
    clear h_x f_x;
    
    %compute current Pf
    Pf_cur = sum(Pfn)/n;
    cov_Pf(n,1) = std(Pfn)/Pf_cur/sqrt(n);
    
    if n >= 30
        convergence_check = mean(cov_Pf(end-9:end));
    end
end


ttt = 0;
Failure_Prob = Pf_cur;
while Failure_Prob < 1
    ttt = ttt + 1; 
    Failure_Prob = Failure_Prob*10;
end
output.Pf = Pf_cur;
output.ci_95 = [(1-2*cov_Pf(n))*Pf_cur,(1+2*cov_Pf(n))*Pf_cur];
output.num = n;
output.u = u;
output.z = z;
output.Pi=Pi;
output.sig=sig;
output.mu_cur=mu_cur;

uuu = 0;
low_bound = output.ci_95(1);
while low_bound < 1
    uuu = uuu + 1;
    low_bound = low_bound*10;
end

vvv = 0;
up_bound = output.ci_95(2);
while up_bound < 1
    vvv = vvv + 1;
    up_bound = up_bound*10;
end



fprintf('\t- Failure probability \t :\t%1.2f x 10^-%i\n', Failure_Prob,ttt);
fprintf('\t- 95%% Confidence Interval:\t[%1.2f x 10^-%i, %1.2f x 10^-%i]\n',...
   low_bound,uuu,up_bound,vvv);
fprintf('\t- Number of samples \t :\t%i\n\n', n);