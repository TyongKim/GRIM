function G_x = lsf_CEAISGM(u,probdata,Lo,N,gfundata,system_data)


%transformation using Nataf Distribution Model
u = u';
for i = 1:N
    x(:,i) = u_to_x(u(:,i),probdata,Lo);
end
clear i
x = x';

if strcmp(system_data.type,'series') || strcmp(system_data.type,'parallel')...
        || strcmp(system_data.type,'general')
    N_lsf = size(gfundata,2);
    
    for i = 1:N_lsf
        G_temp(:,i) = eval(gfundata(i).expression);
    end
    
    G_x = system_type_CEAISGM(system_data,G_temp); 
        
else
    %G_x = 5 - x(:,2) - 0.5.*(x(:,1) - 0.1).^2;
    G_x = eval(gfundata(1).expression);
end