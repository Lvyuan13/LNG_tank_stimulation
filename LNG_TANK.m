clear;
clc;
%% set global varibles
global fluid1 fluid2 fluid3 
global R L dt V Tair Kwet Kdry s
fluid1='methane';
fluid2='ethane';
fluid3='nitrogen';
L=1.516;                                                  %length of the LNG storage tank
R=0.251;                                                  %radius of the LNG storage tank
V=pi*L*R^2;  
Kwet=0.03;                                               %heat trainsimision coefficent of liquid zone
Kdry=0.03;                                               %heat trainsimision coefficent of vapor zone
dt=0.5;      
Tair=298.15;
s=0.21;
%%set required constants for iterative caculation
%% set stimulation time and time limits for iteration
itmax=80;   
limit=22; %100s
imax=limit/dt;  
%% set initial values
TV0=153;
TL0=153;
VL0=0.1;
Vv0=V-VL0;
%initial value of pressure should be set as the saturaed pressure of liquid
x0=[0.974652  0.000692 0.024656];
p0=refpropm('P','T',TV0,'Q',0,fluid1,fluid2,fluid3,x0);
[x0,y0]=refpropm('X','T',TV0,'Q',0,fluid1,fluid2,fluid3,x0)
DV0=refpropm('D','T',TV0,'Q',1,fluid1,fluid2,fluid3,y0);
DL0=refpropm('D','T',TL0,'Q',0,fluid1,fluid2,fluid3,x0);
T0=TV0;
M0=VL0*DL0+Vv0*DV0;
e=VL0*DL0/(M0);
q=1-e;
D0=M0/V;
z0=x0*e+y0*q;
h0=refpropm('H','T',T0,'P',p0,fluid1,fluid2,fluid3,z0);
hv0=refpropm('H','T',TV0,'Q',1,fluid1,fluid2,fluid3,y0);
%% 
for i=1:imax
    time(i)=i*dt;
    TV(i)=TV0;
    TL(i)=TL0;
    P(i)=p0;
    dltaT(i)=TV0-TL0;
   xliq1(i)=x0(1);
   xliq2(i)=x0(2);
   xliq3(i)=x0(3);
   xvap1(i)=y0(1);
   xvap2(i)=y0(2);
   xvap3(i)=y0(3);
    M(i)=M0;
    Passume=p0;
    [p1,TL1,TV1,BOG0,x1,y1,VL1]=runstim(TL0,x0,TV0,y0,VL0,Passume);
    M1=M0-BOG0;
    p0=p1;
    TL0=TL1;
    TV0=TV1;
    x0=x1;
    y0=y1;
    VL0=VL1;
    M0=M1;
    BOG(i)=BOG0/dt;
end
%pause
%% results and figures
figure(1);
plot(time,TV);
hold on;
plot(time,TL);
title('change of T');
xlabel('t/s')
ylabel('T/K')
figure(2);
plot(time,P)
xlabel('t/s')
ylabel('kPa')
title('change of P');

figure(3);
plot(time,xvap1);
hold on;
plot(time,xvap2);
hold on;
plot(time,xvap3);
xlabel('t/s')
title('change of vapor composition')

figure(4);
plot(time,xliq1);
hold on;
plot(time,xliq2);
hold on;
plot(time,xliq3);
xlabel('t/s')
title('change of liquid composition')

figure(5);
plot(time,BOG);
xlabel('t/s')
ylabel('kg/s')
title('change of BOG');

figure(6);
plot(time,dltaT);
xlabel('t/s')
ylabel('K')
title('temperature differnce between two zones');

figure(7);
plot(time,M);
xlabel('t/s');
ylabel('kg');
title('remaining mass in the tank')