
%7/12/2011 WWC file created

%script to solve the NACA 0018 at Re1e5
clc
clear all
close all

global lambda B U rho R mu r_hub
U=7;            %wind speed (m/s)
rho=1.229;      %density of air (kg/m^3)
mu=1.73e-5;     %viscosity of air (Pa-s)
R=0.5;          %rotor radius (m)
r_hub=3*2.54e-2;%hub radius (m)
B=3;            %number of blades
alpha_D=10;   %design angle of attack (degrees) from design Re of 85e3
C_L_D=1.0479;   %lift coefficient from design Re of 85e3
            
                
                
lambda=9;       %design tip speed ratio
s=1;            %number of sections to ignore from the end
alpha_start=8;  %angle of attack to start the BEM sweep (deg.)
alpha_end=15; %angle of attack to end the BEM sweep (deg.)
alpha_step=.5;  %angle step between alpha_start and alpha_end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function calls%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_R,c_R,theta_p,a,a_prime,BETZ]=betz(C_L_D,alpha_D);
[BEM Re]=bem(r_R,c_R,theta_p,alpha_start,alpha_end,alpha_step);
[C_p alpha_p C_L_p C_D_p F_p a_prime_p a_p s]=cpcalc(BEM,r_R,theta_p,Re,s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RotorPowerCoefficient=C_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2)
plot(r_R(1:length(r_R)-s),alpha_p)
grid on
title('Predicted alpha vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Predicted Alpha (deg.)')

figure
plot(r_R,a,'b',r_R,a_prime,'g')
hold on
plot(r_R(1:end-s),a_p,'r',r_R(1:end-s),a_prime_p,'c')
hold off
grid on
title('a and a'' vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Induction Factors')
legend('a: Design','a'': Design','a: Predicted','a'': Predicted')



figure
subplot(2,2,1)
plot(r_R(1:length(r_R)-s),C_L_p)
grid on
title('Predicted C_L vs. Relative Radius')
xlabel('Relative Radius r/R')
ylabel('Predicted C_L')

subplot(2,2,2)
plot(r_R(1:length(r_R)-s),C_D_p)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find C_p(lambda)

a_s=[0 0 0 0 5 6 7 7 8 8];
a_e=15*ones(10);
sec=[1 1 1 1 2 2 2 2 1 1];

for i=1:10
    lambda=i;
    alpha_start=a_s(i);
    alpha_end=a_e(i);
    s=sec(i);
    
    [r_R,c_R,theta_p,a,a_prime]=betz(C_L_D,alpha_D);
    [BEM Re]=bem(r_R,c_R,theta_p,alpha_start,alpha_end,alpha_step);
    [C_p alpha_p C_L_p C_D_p F_p a_prime_p a_p s]=...
        cpcalc(BEM,r_R,theta_p,Re,s);
    
    x(i)=i;
    y(i)=C_p;
end

figure
plot(x,y)
grid on
title('C_p vs. \lambda')
xlabel('\lambda')
ylabel('C_p')
