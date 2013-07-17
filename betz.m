
%7/12/2011 WWC file created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%BETZ ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_R,c_R,theta_p,a,a_prime,BETZ]=betz(C_L_D,alpha_D)

global lambda B U rho R mu r_hub

n=20; %number of sections of the blade
theta_po=atan(2/3/lambda)*180/pi-alpha_D; %(degrees)
r_R=[r_hub/R:(1-r_hub/R)/(n-1):1]';
phi=2/3*atan(1/lambda./(r_R))*180/pi;
c_R=8*pi*r_R/B/C_L_D.*(1-cos(phi*pi/180));
theta_T=phi-theta_po-alpha_D;
theta_p=theta_T+theta_po;
a=1./(1+2*pi*r_R*4.*sin(phi*pi/180).^2/B./(c_R)./cos(phi*pi/180));
a_prime=(1-3*a)./(4*a-1);
U_rel=U*(1-a)./sin(phi*pi/180);
Re_c=rho*c_R*R.*U_rel/mu;

BETZ=[r_R phi c_R theta_T theta_p a a_prime U_rel Re_c];


end
