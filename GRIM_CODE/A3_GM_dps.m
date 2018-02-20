function mu_init_GM = A3_GM_dps(input)
%generate points uniformly on a hypershpere using 
%Latin Hypercube Sampling

%--INPUT--
N_rvs = input.N_rvs;
N_dps = input.N_dps;
R = input.radius;
clear input;

while 1
    if N_dps == 1
        dps_up = lhsnorm(zeros(1,N_rvs),eye(N_rvs),1);
    else
        dps_up = lhsnorm(zeros(1,N_rvs),eye(N_rvs),N_dps,1);
    end
    norm = sqrt(sum((dps_up.^2)'))';

    for iii = 2:N_rvs
        norm(:,iii) = norm(:,1);
    end
    
    if max(isnan(dps_up./norm))
        clear dps_up norm;
    else
        %--OUTPUT--
        mu_init_GM = R*dps_up./norm;
        return;
    end
end
