%Wind Turbine Iterative Design

%7/19/2011 file created

%function C_p=wtd(lambda)

clc
clear all
close all

U=5; %wind speed (m/s)
rho=1.229; %density of air (kg/m^3)
mu=1.73e-5; %viscosity of air (Pa-s)
R=0.5; %rotor radius (m)
r_hub=3*2.54e-2; %hub radius (m)
B=3; %number of blades
lambda=5; %design tip speed ratio
alpha=7.6; %design angle of attack (degrees)
C_L=1.040; %lift coefficient (1.1617) based on Reynolds Number of 1e5 and 
            %angle of attack of 14.4 degrees

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%BETZ ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%

theta_po=atan(2/3/lambda)*180/pi-alpha;

r=.025/2:(.5-.025/2)/39:.5;
phi=2/3*atan(1/lambda./(r/R))*180/pi;
c=8*pi*r/B/C_L.*(1-cos(phi*pi/180));
theta_T=phi-theta_po-alpha;
theta_p=theta_T+theta_po;
a=1./(1+2*pi*r/R*4.*sin(phi*pi/180).^2/B./(c/R)./cos(phi*pi/180));
a_prime=(1-3*a)./(4*a-1);
U_rel=U*(1-a)./sin(phi*pi/180);
Re_c=rho*c.*U_rel/mu;

betz=[(r/R)'  phi'  (c/R)'  theta_T' theta_p'  a'  a_prime'...
    c' U_rel']

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%BEM ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%

BEM(20)=struct('alpha',[],'phi',[],'F',[],'sigma_prime',[],'C_L',[],...
    'a_prime',[],'a',[],'Re_c',[]);
for i=1:length(r)
    BEM(i).alpha=[0:.2:20]';
    BEM(i).phi=(theta_p(i)+BEM(i).alpha)*pi/180;
    BEM(i).F=2/pi*acos(exp(-B/2*(1-r(i)/R)/(r(i)/R)./sin(BEM(i).phi)));
    BEM(i).sigma_prime=B*c(i)/(2*pi*r(i));
    BEM(i).C_L=4*BEM(i).F.*sin(BEM(i).phi).*(cos(BEM(i).phi)-lambda*r(i)...
        /R*sin(BEM(i).phi))./(BEM(i).sigma_prime*(sin(BEM(i).phi)+lambda...
        *r(i)/R*cos(BEM(i).phi)));
    BEM(i).a_prime=1./(4*BEM(i).F.*cos(BEM(i).phi)/BEM(i).sigma_prime./...
        BEM(i).C_L-1);
    BEM(i).a=BEM(i).a_prime*lambda*r(i)/R./tan(BEM(i).phi);
    BEM(i).Re_c=rho*c(i)*U*(1-BEM(i).a)./sin(BEM(i).phi)/mu;
    
end

s=18;
[BEM(s).alpha BEM(s).phi BEM(s).F BEM(s).C_L BEM(s).a_prime BEM(s).a ]

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREDICTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%table look-up for C_L(Re_c_design,alpha) using interp 1
% consider table look-up for C_L(Re_c,alpha) using interp2 (optional)
load polar.mat

%the last section is ommited because of zero lift
for i=1:length(r)-1
    BEM(i).alpha_root=interp1(Re65e3(:,2)-BEM(i).C_L,Re65e3(:,1),0);
    BEM(i).C_L_root=interp1(BEM(i).alpha,BEM(i).C_L,BEM(i).alpha_root);
    BEM(i).C_D_root=interp1(Re65e3(:,1),Re65e3(:,3),BEM(i).alpha_root);
    
    BEM(i).F_root=interp1(BEM(i).alpha,BEM(i).F,...
        BEM(i).alpha_root);
    BEM(i).a_prime_root=interp1(BEM(i).alpha,BEM(i).a_prime,...
        BEM(i).alpha_root);
    BEM(i).a_root=interp1(BEM(i).alpha,BEM(i).a,...
        BEM(i).alpha_root);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct predicted alpha and predicted C_L vs relative radius graph
for i=1:length(r)-1
    alpha_predicted(i)=BEM(i).alpha_root;
end

for i=1:length(r)-1
    C_L_predicted(i)=BEM(i).C_L_root;
end

for i=1:length(r)-1
    C_D_predicted(i)=BEM(i).C_D_root;
end

% figure
% plot(r(1:length(r)-1)/R,alpha_predicted)
% grid on
% title('Predicted alpha vs Relative Radius')
% xlabel('Relative Radius r/R')
% ylabel('Predicted Alpha (deg.)')
% 
% figure
% plot(r(1:length(r)-1)/R,C_L_predicted)
% grid on
% title('Predicted C_L vs Relative Radius')
% xlabel('Relative Radius r/R')
% ylabel('Predicted C_L')
% 
% figure
% plot(r(1:length(r)-1)/R,C_D_predicted)
% grid on
% title('Predicted C_D vs Relative Radius')
% xlabel('Relative Radius r/R')
% ylabel('Predicted C_D')

%sigma_prime tabulation

for i=1:length(r)-1
    sigma_prime_predicted(i)=BEM(i).sigma_prime;
end

%correction factor F tabulation
for i=1:length(r)-1
    F_predicted(i)=BEM(i).F_root;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate a' and a from predicted angles for each section
%NOTE: phi has been changed to derive its value from alpha_root and
%theta_p



for s=1:length(r)-1
    
    a_prime_predicted(s)=1./(4*BEM(s).F_root.*cos((BEM(s).alpha_root+...
        theta_p(s))*pi/180)/BEM(s).sigma_prime./BEM(s).C_L_root-1);
    a_predicted(s)=a_prime_predicted(s)*lambda*r(s)/R./...
        tan((BEM(s).alpha_root+theta_p(s))*pi/180);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solve for Cp%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r(1:length(r)-1)'/R alpha_predicted' C_L_predicted']

C_p_integrand=F_predicted.*a_prime_predicted.*(1-a_predicted).*...
    (1-cot((alpha_predicted+theta_p(1:length(theta_p)-1))*pi/180)./...
    (C_L_predicted./C_D_predicted)).*(r(1:length(r)-1)/R).^3;

C_p=8*lambda^2*trapz(r(1:length(r)-1)/R,C_p_integrand)








