function R = A2_radius(input)
%to determine the scaled length of radius from the origin
%to the single gaussian mean

%--INPUT--
mu_cur = input.mu_SG;
R_factor = input.R_factor;
clear input;

%compute radius
R = R_factor*sqrt(sum(mu_cur.^2));