% Barrier geometric model. Linear shoreface. Cte sea-level rise.
clc; clear all;
%%Input physical parameters%%%%%%%%%%%%%%%
B=0.001; %Basement slope
Dt=10;% Toe depth (meters). Typically in the range 10-20m
We=400; %Equilibrium width (meters)
He=2;  %Equilibrium heigth (meters)
Ae=0.02; %Equilibrium shoreface slope
Qow_max=100; %Maximum overwash flux (m^2/year)
Vd_max=300; %Maximum deficit volume (m^2/year)
K=1000;  %Shoreface Flux constant (m^2/year)
a=0.002;
b=0;
%%%% Dynamic equilibrium%%%%%%%%%%%%
W_de=We-a/B*Vd_max/Qow_max;
H_de=He-a*Vd_max/Qow_max;
p=4*K*B/Dt*((H_de+Dt)/(2*H_de+Dt));
bb=p*Ae-a;
A_de=(bb+(bb^2+8*B*p*a)^0.5)/(2*p);
Db_de=Dt*(1-B/A_de)-B*W_de;
Qow_de=a/B*(H_de+Db_de)+a*W_de;
Qsf_de=K*(Ae-A_de);
%%% Initial conditions %%%%%%%%%%%%
A=1.0*Ae;W=1.0*We;H=1.0*He; %Barrier initially in equlibrium
xt=0;xs=Dt/A;xb=xs+W;xso=xs;Z=1.0*Dt;
Db=Z-B*xb;  
Vexp=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computational parameters%%%%%%%%%%%%%%%
Tmax=1000; %Runing time (years)
dt=0.01*1;    
t=0:dt:Tmax;
n=length(t);
Heigth=zeros(1,n); Width=zeros(1,n);OverwashFlux=zeros(1,n); VolDef=zeros(1,n);
VolDef_H=zeros(1,n);VolDef_B=zeros(1,n);XB=zeros(1,n);XT= zeros(1,n);XS=zeros(1,n);
SFslope=zeros(1,n);Heigth(1)=He; Width(1)=He;XB(1)=xb;XT(1)=xt;XS(1)=xs;
SFFlux=zeros(1,n);SFFlux_bis=zeros(1,n);SHrate=zeros(1,n);W1=zeros(1,n);
XTrate=zeros(1,n);Ztoe=zeros(1,n);QT=zeros(1,n);
XBrate=zeros(1,n);Hrate=zeros(1,n);XS_de=zeros(1,n);
Vexport=zeros(1,n);

Weq=We*ones(1,n);Aeq=Ae*ones(1,n);Heq=He*ones(1,n);
Wdeq=W_de*ones(1,n);Hdeq=H_de*ones(1,n);Adeq=A_de*ones(1,n);Dbde=Db_de*ones(1,n);
Qowde=Qow_de*ones(1,n);Qsfde=Qsf_de*ones(1,n);xdotde=a/B*ones(1,n);

Depth_b=zeros(1,n);
SFslope_bis=zeros(1,n);

%% Main code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
zdot=a+2*b*t(i); %Base-level rise rate (m/year)
Z=Z+zdot*dt;
Vd_H=max(0,(He-H)*W);
Vd_B=max(0,(We-W)*(H+Db));
Vd=Vd_H+Vd_B;
if Vd<Vd_max; Qow_H=Qow_max*Vd_H/Vd_max;Qow_B=Qow_max*Vd_B/Vd_max;else
Qow_H=Qow_max*Vd_H/Vd;Qow_B=Qow_max*Vd_B/Vd;end
Qow=Qow_H+Qow_B;
Qsf=K*(Ae-A);
Hdot=Qow_H/W-zdot;
xbdot=Qow_B/(H+Db);
xsdot=2*Qow/(Dt+2*H)-4*Qsf*(H+Dt)/(2*H+Dt)^2;
xtdot=2*Qsf*(1/(Dt+2*H)+1/Dt)+2*zdot/A;
H=H+Hdot*dt;if H<0;tdrown_H = t(i);break; end
xb=xb+xbdot*dt;
xs=xs+xsdot*dt;
xt=xt+xtdot*dt;
A=Dt/(xs-xt);
W=xb-xs;if W<0;tdrown_W = t(i);break; end
Db=Z-xb*B;
ht=Z-Dt-B*xt;
qt=xtdot*ht;
Vexp=Vexp+qt*dt;
% Variable storage %%%%%%%%%%%%%%%%%%%%
Heigth(i)=H; Width(i)=W; OverwashFlux(i)=Qow;SFFlux(i)=Qsf;
VolDef(i)=Vd;VolDef_H(i)=Vd_H;VolDef_B(i)=Vd_B; SFslope(i)=A;
XB(i)=xb;XT(i)=xt;XS(i)=xs;SHrate(i)=xsdot;W1(i)=W;XTrate(i)=xtdot;
XBrate(i)=xbdot;Hrate(i)=Hdot;Depth_b(i)=Db;Ztoe(i)=ht;
XS_de(i)=xso+zdot/B*t(i);Vexport(i)=Vexp;
QT(i)=qt;
end
hold on
figure(1)
subplot(3,2,1);
hold on
plot(t,XS,'k', 'linewidth',2)
plot(t,XS_de,'k', 'linewidth',1)
box on
xlabel('time (y)')
title('Shoreline position x_S (m)')

subplot(3,2,2);
hold on
plot(t,Heigth,'k', 'linewidth',2)
plot(t,Hdeq,'k', 'linewidth',1)
plot(t,Heq,'k--', 'linewidth',2)
box on
xlabel('time (y)')
title('Height (m)')

subplot(3,2,3);
hold on
plot(t,Weq,'k--', 'linewidth',2)
plot(t,Wdeq,'k', 'linewidth',1)
plot(t,Width,'k', 'linewidth',2)
xlabel('time (y)')
title('Width W (m)')
box on

subplot(3,2,4);
hold on
plot(t,Aeq,'k--', 'linewidth',2)
plot(t,Adeq,'k', 'linewidth',1)
plot(t,SFslope,'k', 'linewidth',2)
box on
xlabel('time (y)')
title('Shoreface slope \alpha (-)')

subplot(3,2,5);
hold on
plot(t,OverwashFlux,'k', 'linewidth',2)
plot(t,Qowde,'k', 'linewidth',1)
box on
xlabel('time (y)')
title('Overwash flux (Q_{ow}) (m^3/m/y)')

subplot(3,2,6);
hold on
plot(t,SFFlux,'k', 'linewidth',2)
% plot(t,QT,'k', 'linewidth',2)
% plot(t,Depth_b,'r', 'linewidth',2)
plot(t,Qsfde,'k', 'linewidth',1)
box on
xlabel('time (y)')
title('Shoreface flux (Q_{sf}) (m^3/m/y)')


