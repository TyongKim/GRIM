%% C1 sensitivity

% Based on the near optimal Gaussian mixture model, one can esitmate
% the relative importance of random variables.

function sensitivity=C1_variance_decompose(Lo,output,N_presam,probdata,gfundata,system_data)

z = mnrnd(N_presam,output.Pi);         
zz = cumsum(z);
for i = 1:length(output.sig)              
    if i ==1
        idx = 1:zz(i);
    else idx = (zz(i-1)+1):zz(i);
    end
    temp_u = mvnrnd(output.mu_cur(i,:), output.sig(i).cur, z(i));
    u(idx,:) = temp_u;
    temp_G = lsf_CEAISGM(temp_u,probdata,Lo,z(i),gfundata,system_data);
    G(idx,:) = temp_G;
    G_var(i) = var(temp_G);
    G_mean(i) = mean(temp_G);
    
end


for i = 1:length(output.sig)
    temp_decom_var=0;
    for j=1:length(output.sig)
        temp_decom_var = temp_decom_var + G_mean(i)*output.Pi(i)*G_mean(j)*output.Pi(j);
    end
    decom_var(i) = G_var(i)*output.Pi(i)+G_mean(i)^2*output.Pi(i)-temp_decom_var;
end

var_par_fac = decom_var/sum(decom_var);


for i=1:length(output.sig)
    criti_cov=output.sig(i).cur;
    
    [kuu,muu]=eig(criti_cov);
    muu_temp = diag(muu);
    kty = find(muu_temp==min(muu_temp));
    aa=kuu(:,kty);
    if sign(aa(1,1))==sign(output.mu_cur(i,1))
        alpha(:,i) = kuu(:,kty);
    else
        alpha(:,i) = -kuu(:,kty);
    end

    gamma(:,i) = (alpha(:,i)'*Lo^(-1))'/norm(alpha(:,i)'*Lo^(-1));
end
temp_weigth_gamma=[];
for i = 1:length(var_par_fac)
    temp_weigth_gamma(:,i) =[gamma(:,i).^2*var_par_fac(i); sum(gamma(:,i).^2*var_par_fac(i))];
end
temp_weigth_gamma(:,size(temp_weigth_gamma,2)+1)=sum(temp_weigth_gamma');

sensitivity.alpha = alpha;
sensitivity.gamma = gamma;
sensitivity.var_par_fac = var_par_fac;
sensitivity.all = temp_weigth_gamma;

