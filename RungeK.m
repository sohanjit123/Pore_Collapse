clc,clear all

%%%Dynamic pore collapse
d=1e-9; %nanosecond
alpha0=1.2;

a0=550e-6;
rho=8.93e3; % density for SLG
C=385; %specific heat for SLG
Tm=1356; %Melting Temp

T0=293; %Ambient Temperature

G=30e9; % Shear Modulus
Y=30e6;% Yield Strength

T0=293; %Ambient Temperature


 Y1=30e6/(1-(T0/Tm)); %Yield Strength for SLG Michael Ashby 2013


  A=-1.626;
    B=4239;
 TO=265.7;


% Y1=400e6;
% 
%  B=4239;

%   B=5765;

% eta_m=2e-3; %change later

eta_m=10^(A+B/(Tm-273-TO));

const=(d/(rho*C));

Pf=1*1e9; 
 
tau=0.25e-6;

time=0:d:tau;
tsteps=length(time);
P=linspace(0,Pf,tsteps);

 alpha=zeros(length(P),1);

 a=zeros(length(P),1);
 b=zeros(length(P),1);
 x=zeros(length(P),1);


alpha(1)=alpha0;


a(1)=a0;
% a(2)=a(1);
b(1)=bee(a(1),a0,alpha0);
x(1)=0;

r0=linspace(a(1),b(1),100);

T=zeros(length(P),length(r0));

rdotbyr=zeros(length(P),length(r0));

T(1,:)=T0;

 
for k=1:length(time)


  %a loop%
  
  k1=x(k);
  k2=x(k);
  k3=x(k);
  k4=x(k);

%   a(k+1)=a(k)+(d/6)*(k1+2*k2+2*k3+k4);

a(k+1)=a(k)+d*x(k);

  alpha(k+1)= 1+((alpha0-1)*(a(k+1)^3)/(a0^3));

  %x loop

  r=zeros(length(r0),1);
% 
      for i=1:length(r0)
             r(i)=arr(r0(i),a0,a(k));
      end
% 
        b(k)=bee(a(k),a0,alpha0);
        t1=zeros(length(r0),1);
        t2=zeros(length(r0),1);

        %term 1
        
        for i=1:length(r)
            t1(i)=yield(T(k,i),Y1,Tm)/r(i);
%             t2(i)=visc(T(k,i),Tm,B,eta_m)*a(k)*a(k)/(r(i)^4);
            t2(i)=visc2(T(k,i),TO,A,B)*a(k)*a(k)/(r(i)^4);
        end

        term1=2*trapz(r,t1);
        term2=-12*trapz(r,t2);

        term3=rho/2*(1-((a(k)^4)/(b(k)^4)));


        l1=(((term1+term2*x(k)+term3*x(k)*x(k)-P(k))/(rho*(1-(a(k)/b(k)))))-(2*x(k)*x(k)))/a(k);
        l2=(((term1+term2*(x(k)+d/2*l1)+term3*(x(k)+d/2*l1)*(x(k)+d/2*l1)-P(k))/(rho*(1-(a(k)/b(k)))))-(2*(x(k)+d/2*l1)*(x(k)+d/2*l1)))/a(k);
        l3=(((term1+term2*(x(k)+d/2*l2)+term3*(x(k)+d/2*l2)*(x(k)+d/2*l2)-P(k))/(rho*(1-(a(k)/b(k)))))-(2*(x(k)+d/2*l2)*(x(k)+d/2*l2)))/a(k);
        l4=(((term1+term2*(x(k)+d*l3)+term3*(x(k)+d*l3)*(x(k)+d*l3)-P(k))/(rho*(1-(a(k)/b(k)))))-(2*(x(k)+d*l3)*(x(k)+d*l3)))/a(k);

        x(k+1)=x(k) + (d/6)*(l1+2*l2+2*l3+l4);
        

%Temperature loop%    
for i=1:length(r0)

% 
%      m1=((-2*yield(T(k,i),Y1,Tm)*x(k)*(a(k)^2)/(r(i)^3))+(12*visc(T(k,i),Tm,B,eta_m)*x(k)*x(k)*(a(k)^4)/(r(i)^6)))/(rho*C);
%      m2=((-2*yield(T(k,i),Y1,Tm)*(x(k)+d/2*m1)*(a(k)^2)/(r(i)^3))+(12*visc(T(k,i),Tm,B,eta_m)*(x(k)+d/2*m1)*(x(k)+d/2*m1)*(a(k)^4)/(r(i)^6)))/(rho*C);
%      m3=((-2*yield(T(k,i),Y1,Tm)*(x(k)+d/2*m2)*(a(k)^2)/(r(i)^3))+(12*visc(T(k,i),Tm,B,eta_m)*(x(k)+d/2*m2)*(x(k)+d/2*m2)*(a(k)^4)/(r(i)^6)))/(rho*C);
%      m4=((-2*yield(T(k,i),Y1,Tm)*(x(k)+d*m3)*(a(k)^2)/(r(i)^3))+(12*visc(T(k,i),Tm,B,eta_m)*(x(k)+d*m3)*(x(k)+d*m3)*(a(k)^4)/(r(i)^6)))/(rho*C);


     m1=((-2*yield(T(k,i),Y1,Tm)*x(k)*(a(k)^2)/(r(i)^3))+(12*visc2(T(k,i),TO,A,B)*x(k)*x(k)*(a(k)^4)/(r(i)^6)))/(rho*C);
     m2=((-2*yield(T(k,i),Y1,Tm)*(x(k)+d/2*m1)*(a(k)^2)/(r(i)^3))+(12*visc2(T(k,i),TO,A,B)*(x(k)+d/2*m1)*(x(k)+d/2*m1)*(a(k)^4)/(r(i)^6)))/(rho*C);
     m3=((-2*yield(T(k,i),Y1,Tm)*(x(k)+d/2*m2)*(a(k)^2)/(r(i)^3))+(12*visc2(T(k,i),TO,A,B)*(x(k)+d/2*m2)*(x(k)+d/2*m2)*(a(k)^4)/(r(i)^6)))/(rho*C);
     m4=((-2*yield(T(k,i),Y1,Tm)*(x(k)+d*m3)*(a(k)^2)/(r(i)^3))+(12*visc2(T(k,i),TO,A,B)*(x(k)+d*m3)*(x(k)+d*m3)*(a(k)^4)/(r(i)^6)))/(rho*C);


    T(k+1,i)=T(k,i)+(d/6)*(m1+2*m2+2*m3+m4);

end


if alpha(k+1)<1
    break;
    k
end


end








% 
plot(time(1:(k))/1e-9,alpha(1:(k)),'r*-','LineWidth',3)
xlabel('Time (nanoseconds)')
ylabel('Distension \alpha')
ax=gca;
ax.FontSize = 16.0;
ax.LineWidth = 1.0;

function eta=visc2(T,TO,A,B)
% eta= eta_m*exp(B*((1/T)-(1/Tm)));
eta=10^(A+B/(T-273-TO));
end


function r=arr(r0,a0,a)
r=(r0^3-a0^3+a^3)^(1/3);
end

function b = bee(a,a0,alpha0)
b = (a^3+(a0^3/(alpha0-1)))^(1/3);
end

function Y=yield(T,Y1,Tm)
Y=Y1*(1-(T/Tm));
end

function eta=visc(T,Tm,B,eta_m)
eta= eta_m*exp(B*((1/T)-(1/Tm)));
end