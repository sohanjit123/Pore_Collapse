clc,clear all

%%%Dynamic pore collapse
d=1e-9; %nanosecond
alpha0=2;

alpha_star=1.66;

R=550e-6;

% c=R*nthroot((2-alpha_star),3);
% a0=R*nthroot((alpha0-alpha_star+1),3);
% b0=R*nthroot(alpha0,3);

a0=20e-6;
rho=8.93e3; % density for SLG
C=385; %specific heat for SLG
Tm=1356; %Melting Temp

T0=293; %Ambient Temperature

G=30e9; % Shear Modulus
Y=30e6;% Yield Strength

Y1=30e6/(1-(T0/Tm)); %Yield Strength for SLG Michael Ashby 2013
% Y1=400e6;
% Viscocity Parameters------https://en.wikipedia.org/wiki/Soda%E2%80%93lime_glass-------%

%Napolitano

 A=-1.626;
    B=4239;
 TO=265.7;

%     B=5765;

eta_m=10^(A+B/(Tm-273-TO));

% eta_m=2e-3; %change later

const=(d/(rho*C));


%  Pcrit=(2*Y1/3)*log(alpha0/(alpha0-1)); 
 Pf=14.3*1e9; 
 
%  P=[0:Prate*0.01:1*Prate];
tau=0.25e-6;

time=0:d:tau;
tsteps=length(time);
   P=linspace(0,Pf,tsteps);

%    P=Pf*ones(length(time),1);

alpha1=2*G*alpha0/(2*G+Y);

P1=2/3*Y*log(alpha1/(alpha1-1));

%  alpha=zeros(length(P),1);
%  time=zeros(length(P),1);
%  a=zeros(length(P),1);
%  b=zeros(length(P),1);


alpha(1)=alpha0;
alpha(2)=alpha(1);
% time(1)=0;

a(1)=a0;
a(2)=a(1);
b(1)=bee(a(1),a0,alpha0);

r0=linspace(a(1),b(1),100);
% r0(length(r0))=b(1);

T=zeros(length(P),length(r0));

rdotbyr=zeros(length(P),length(r0));

 T(1,:)=T0;


 for k=1:(length(P)-1)

  adot=(a(k+1)-a(k))/d;

  r=zeros(length(r0),1);

  for i=1:length(r0)
         r(i)=arr(r0(i),a0,a(k));
         rdotbyr(k,i)=a(k)*a(k)*(adot)/(r(i)^3);
  end

 %size loop

b(k)=bee(a(k),a0,alpha0);
t1=zeros(length(r0),1);
t2=zeros(length(r0),1);

for i=1:length(r)
    t1(i)=yield(T(k,i),Y1,Tm)/r(i);
%       t2(i)=visc(T(k,i),Tm,B,eta_m)*a(k)*a(k)*adot/(r(i)^4);
      t2(i)=visc2(T(k,i),TO,A,B)*a(k)*a(k)*adot/(r(i)^4);
end


term1=2*trapz(r,t1);
term2=-12*trapz(r,t2);

term3=rho*(adot^2)/2*(1-((a(k)^4)/(b(k)^4)));

term4=term1+term2+term3-P(k);
term4=term4/rho/(1-(a(k)/b(k)));

term4=(term4-2*adot*adot)/a(k);
a(k+2)=term4*d*d+2*a(k+1)-a(k);

alpha(k+2)= 1+((alpha0-1)*(a(k+2)^3)/(a0^3));

%Temperature loop%    
for i=1:length(r0)

     first_term=-2*yield(T(k,i),Y1,Tm)*rdotbyr(k,i);
%    second_term=12*visc(T(k,i),Tm,B,eta_m)*rdotbyr(k,i)*rdotbyr(k,i);
     second_term=12*visc2(T(k,i),TO,A,B)*rdotbyr(i)*rdotbyr(i);
   
    T(k+1,i)=T(k,i)+const*(first_term+second_term);

end

if alpha(k+2)<1
    break;
    k
end

 end


% 
plot(time(1:(k+1))/1e-6,alpha(1:(k+1)),'b','LineWidth',3)
xlabel('Time (microseconds)')
ylabel('Distension \alpha')
ax=gca;
ax.FontSize = 16.0;
ax.LineWidth = 1.0;


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

function eta=visc2(T,TO,A,B)
% eta= eta_m*exp(B*((1/T)-(1/Tm)));
eta=10^(A+B/(T-273-TO));
end

function eta=visc3(T,eta_m,E,R)
% eta= eta_m*exp(B*((1/T)-(1/Tm)));
eta=eta_m*exp(E/R/T);
end