function G_out = system_type_CEAISGM(system_data,G_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function is used to determine the system type for simulation     %
%   problems. For series system, the combined system value of limit       % 
%   states is the minimum among all g(x), i.e. min[g1(x),g2(x),...,gn(x)].%
%   For paralle system, g(x) = max[g1(x),g2(x),...,gn(x)]. For parallel,  %
%   g(x) equals a mixture of maximum and minimum functions.               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G_in = G_in';
system_data.N_lsf = size(G_in,1);

if strcmp(system_data.type,'series')
    
    G_x = min(G_in(1,:),G_in(2,:));
    
    if system_data.N_lsf >= 3
        for i = 3:system_data.N_lsf
            G_x = min(G_x,G_in(i,:));
        end
    end
    G_out = G_x';
    
elseif strcmp(system_data.type,'parallel')
    
    G_x = max(G_in(1,:),G_in(2,:));
    
    if system_data.N_lsf >= 3
        for i = 3:system_data.N_lsf
            G_x = max(G_x,G_in(i,:));
        end
    end
    G_out = G_x';
    
elseif strcmp(system_data.type,'general')
    Eves = system_data.general;
    num_E = length(Eves);
    j = 1;
    k = 1;
    for i = 1:num_E
        if Eves(i) ~= 0
            parallel(j,k) = Eves(i);
            k = k+1;
        elseif Eves(i) == 0
            j = j+1;
            k = 1;
        end
    end  
    clear i j k;
    
    [m, n] = size(parallel);
    
    for j = 1:m
        if parallel(j,2) == 0
            G_temp(j,:) = G_in(parallel(j,1),:);
        else
            for k = 2:n
                if parallel(j,k) ~= 0
                    G_temp_temp = max(G_in(parallel(j,1),:),G_in(parallel(j,2),:));
    
                    if k >= 3
                        G_temp_temp = max(G_temp_temp,G_in(parallel(j,k),:));
                    end
                    G_temp(j,:) = G_temp_temp;

                elseif parallel(j,k) == 0
                    break;
                end
            end
        end
    end
    clear j k;
    
    G_x = min(G_temp(1,:),G_temp(2,:));
    
    if m >= 3
        for j = 3:m
            G_x = min(G_x,G_temp(j,:));
        end
    end
    G_out = G_x';
end

     