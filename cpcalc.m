
%7/12/2011 WWC file created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREDICTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%table look-up for C_L(Re_c_design,alpha) using interp 1
% consider table look-up for C_L(Re_c,alpha) using interp2 (optional)

function [C_p alpha_p C_L_p C_D_p F_p a_prime_p a_p s]...
    =cpcalc(BEM,r_R,theta_p,Re,s)
global lambda B U rho R mu r_hub

%the last section is ommited because of zero lift

ind=1;
for i=1:length(BEM(end-1).phi)
    if BEM(end-1).phi(i)<0
        ind=i+1;
    end
end
    
%s=3;    %number of sections to omit from the end
for i=1:length(r_R)-s   
%     i;
%     length(Re(ind:end,2))
%     length(BEM(i).C_L(ind:end))
%    Re(ind:end,2)-BEM(i).C_L(ind:end)
    alpha_p(i)=interp1(Re(ind:end,2)-BEM(i).C_L(ind:end)...
        ,Re(ind:end,1),0);
    C_L_p(i)=interp1(BEM(i).alpha,BEM(i).C_L,alpha_p(i));
    C_D_p(i)=interp1(Re(ind:end,1),Re(ind:end,3),alpha_p(i));
    F_p(i)=interp1(BEM(i).alpha,BEM(i).F,alpha_p(i));
    a_prime_p(i)=interp1(BEM(i).alpha,BEM(i).a_prime,alpha_p(i));
    a_p(i)=interp1(BEM(i).alpha,BEM(i).a,alpha_p(i));
end

%[alpha_p' C_L_p' C_D_p' F_p' a_prime_p' a_p']




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solve for Cp%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_p_integrand=F_p.*a_prime_p.*(1-a_p).*...
    (1-cot((alpha_p+[theta_p(1:length(theta_p)-s)]')*pi/180)./...
    (C_L_p./C_D_p)).*[r_R(1:length(r_R)-s).^3]';

C_p=8*lambda^2*trapz(r_R(1:length(r_R)-s),C_p_integrand);

end