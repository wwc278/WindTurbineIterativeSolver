
%7/12/2011 WWC file created

%script to use wtd to solve in class problem
clc
clear all
close all

global lambda B U rho R mu
U=5;            %wind speed (m/s)
rho=1.229;      %density of air (kg/m^3)
mu=1.73e-5;     %viscosity of air (Pa-s)
R=0.5;          %rotor radius (m)
r_hub=3*2.54e-2;%hub radius (m)
B=3;            %number of blades
alpha_D=7.6;    %design angle of attack (degrees)
C_L_D=1.04;     %lift coefficient (1.1617) based on Reynolds Number of 1e5 and 
                %angle of attack of 14.4 degrees
            
                
                
lambda=5;      %design tip speed ratio
alpha_start=0;  %angle of attack to start the BEM sweep (deg.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function calls%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_R,c_R,theta_p,a,a_prime]=betz(C_L_D,alpha_D);
[BEM Re]=bem(r_R,c_R,theta_p,alpha_start);
[C_p alpha_p C_L_p C_D_p F_p a_prime_p a_p]=cpcalc(BEM,r_R,theta_p,Re);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RotorPowerCoefficient=C_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
plot(r_R(1:length(r_R)-1),alpha_p)
grid on
title('Predicted alpha vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Predicted Alpha (deg.)')

figure
plot(r_R,a,'b',r_R,a_prime,'g')
hold on
plot(r_R(1:end-1),a_p,'r',r_R(1:end-1),a_prime_p,'c')
hold off
grid on
title('a and a'' vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Induction Factors')
Legend('a: Design','a'': Design','a: Predicted','a'': Predicted')



figure
subplot(2,2,1)
plot(r_R(1:length(r_R)-1),C_L_p)
grid on
title('Predicted C_L vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Predicted C_L')

subplot(2,2,2)
plot(r_R(1:length(r_R)-1),C_D_p)
grid on
title('Predicted C_D vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Predicted C_D')

subplot(2,2,3)
plot(r_R,c_R)
grid on
title('Relative Chord Length vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Relative Chord Length c/R')

subplot(2,2,4)
plot(r_R,theta_p)
grid on
title('Section Twist \theta_p vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Section Twist \theta_p (deg.)')

% for i=1:10
%     lambda=i;
%     if lambda==[1:5]
%         alpha_start=0;
%     else
%         alpha_start=5;
%     end
%     [r_R,c_R,theta_p,a,a_prime]=betz(C_L_D,alpha_D);
%     [BEM Re]=bem(r_R,c_R,theta_p,alpha_start);
%     [C_p alpha_p C_L_p C_D_p F_p a_prime_p a_p]=cpcalc(BEM,r_R,theta_p,Re);
%     
%     x(i)=i;
%     y(i)=C_p;
% end
% 
% figure
% plot(x,y)
