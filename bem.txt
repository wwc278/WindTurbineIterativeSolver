
%7/12/2011 WWC file created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%BEM ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%

function [BEM Re]=bem(r_R,c_R,theta_p,alpha_start,alpha_end,alpha_step)

global lambda B U rho R mu r_hub
% BEM(20)=struct('alpha',[],'phi',[],'F',[],'sigma_prime',[],'C_L',[],...
%     'a_prime',[],'a',[],'Re_c',[]);
load polar.mat
Re=Re1e5;


figure
subplot(1,2,1)
plot(Re(:,1),Re(:,2))
ylim([0 1.6])
hold on
grid on
title('C_L vs alpha')
xlabel('alpha (deg.)')
ylabel('C_L')

for i=1:length(r_R)
    %determines number of angle sweeps per section
    BEM(i).alpha=[alpha_start:alpha_step:alpha_end]'; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BEM(i).phi=(theta_p(i)+BEM(i).alpha)*pi/180;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BEM(i).F=2/pi*acos(exp(-B/2*(1-r_R(i))/r_R(i)./sin(BEM(i).phi)));

    BEM(i).sigma_prime=B*c_R(i)/(2*pi*r_R(i)); %just one number per section
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BEM(i).C_L=4*BEM(i).F.*sin(BEM(i).phi).*(cos(BEM(i).phi)-lambda*...
        r_R(i)*sin(BEM(i).phi))./(BEM(i).sigma_prime*(sin(BEM(i).phi)+lambda...
        *r_R(i)*cos(BEM(i).phi)));

    BEM(i).a_prime=1./(4*BEM(i).F.*cos(BEM(i).phi)/BEM(i).sigma_prime./...
        BEM(i).C_L-1);
    BEM(i).a=BEM(i).a_prime*lambda*r_R(i)./tan(BEM(i).phi);
    BEM(i).Re_c=rho*c_R(i)*R*U*(1-BEM(i).a)./sin(BEM(i).phi)/mu;
    
    plot(BEM(i).alpha,BEM(i).C_L)
    hold on
    %pause
end

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%chop off alpha sweep from low alpha
if BEM(1).alpha(1)~=Re(1,1)
    for i=1:length(Re(:,1))
        if Re(i,1)==BEM(1).alpha(1)
            ind=i;
        end
    end
    Re=Re(ind:end,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%chop off alpha sweep from the high alpha
if BEM(1).alpha(end)~=Re(end,1)
    for i=1:length(Re(:,1))
        if Re(i,1)==BEM(1).alpha(end)
            ind=i;
        end
    end
    Re=Re(1:ind,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Curve Fit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for i=1:length(r_R)
%     f=fittype('a*x^5+b*x^4+c*x^3+d*x^2+e*x+f');
%     BEM(i).c=fit(BEM(i).alpha, BEM(i).C_L, f);
% end
% 
% f=fittype('h*x^7+g*x^6+a*x^5+b*x^4+c*x^3+d*x^2+e*x+f');
% xfoil_fit=fit(Re65e3(:,1),Re65e3(:,2),f);
% 
% figure
% plot(Re65e3(:,1),Re65e3(:,2))
% hold on
% plot(xfoil_fit)



end