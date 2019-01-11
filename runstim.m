function [P1,TL1,TV1,BOG,x1,y1,VL1]=runstim(TL,x,TV,y,VL,Passume)
% with last timepoints prsent calculate next time point's state
% input last time's state:TL x(liquid composition) TV y(gas compostion)
% Passume a assumed prssure of next timepoint the the ouput of next time
% status can be determined using a iterative way

%% rules of naming
%  Subscript 1 means new status
%  no subscript means old status
%  Pasume1 Pasume2 
%  Dm/Dp represents derivative fuction x0=X-f(X)/h(X); f(x) is the origin
%  fuction h(x) is Dm/Dp

global fluid1 fluid2 fluid3 
global R L dt V Tair Kwet Kdry s

% f(x)=fstep(P)-BOG0
% complete last status's parameters and parameters set for iteration
%P=refpropm('P','T',TL,'Q',0,fluid1,fluid2,fluid3,x);
[P,DL,hl]=refpropm('PDH','T',TL,'Q',0,fluid1,fluid2,fluid3,x);%liquid phase
% just for test if Q>1, which is the basic assumption
%Q=refpropm('Q','T',TV,'P',P,fluid1,fluid2,fluid3,y)
% there is a risk for Viscv is invalid when fluid is two phase state so try
% it
try
%[Dv,Viscv]=refpropm('DV','T',TV,'P',P,fluid1,fluid2,fluid3,y);
[kv,Dv,beta,Viscv,Cpv,hv] = refpropm('LDBVCH','T',TV,'P',P,fluid1,fluid2,fluid3,y);%气相参数
catch
[Dv,hv]=refpropm('DH','T',TV,'P',P,fluid1,fluid2,fluid3,y);
Viscv=6.1033e-06;
kv=0.0136;
beta=0.0093;
Cpv=2.0412e+03;
end
ML=DL*VL;%liquid mass for last time point
VV=V-VL;
MV=Dv*VV; %gas mass of last time point

H=geo_RLVltoh0(R,L,VL);%VL->H(liquid phase height
Avl=2*L*sqrt(R^2-(R-H)^2);%areas for heat transfer between gas phase and outside
Awet=2*R*((L+R)*(pi-acos((H-R)/R))+(H-R)*sin(acos((H-R)/R)));
Adry=2*pi*R^2+2*pi*R*L-2*R*((L+R)*(pi-acos((H-R)/R))+(H-R)*sin(acos((H-R)/R)));
m0=flow_TPxP0toFR(P,s,Dv,Viscv);%Deflation ratio calculted from last time point's valve model
BOG0=m0*dt;                %mass of deflation during one time step
eps=1e-7;       %the value of eps can have a great influence
tolr=1e-7;
itmax=80;
% maxium times for iteration
Passume1=Passume;
%% begin main loop 
for i=1:itmax
   %First-order Newton iteration
    DmDP=(fstep(TL,x,TV,y,VL,Passume1)-fstep(TL,x,TV,y,VL,Passume1-eps))/eps;
    Passume2=Passume1-(fstep(TL,x,TV,y,VL,Passume1)-BOG0)/DmDP;
    deltaBOG=fstep(TL,x,TV,y,VL,Passume2)-BOG0;
    if (abs(deltaBOG)<=tolr)
        disp('converged             converged        converged')
        break
    else
        Passume1=Passume2;
        %disp(['step:  ',num2str(i),'         havent converged'])
    end    
end
%% 
P1=Passume2;

g=9.8;%Gravity acceleration
h=6.7296588*kv^(0.857144)*(g*Dv^2*beta*(TV-TL)*Cpv/Viscv)^0.142856; %Gas-liquid two-phaseflow heat transfer coefficient
QVL=Avl*h*(TV - TL)*dt; %Total amount of gas-liquid phase heat transfer from the previous moment to the next time
QLin=Kwet*Awet*(Tair-TL)*dt;
QVin=Kdry*Adry*(Tair-TV)*dt;
%Heat exchange between gas phase and liquid phase and the outside world from the previous moment to the next time

hl1=hl+(QLin+QVL)/ML;  %next timepoint's liuqid phase's enthapy, take liquid phase as a control volume
[x1,y0]=refpropm('X','P',P1,'H',hl1,fluid1,fluid2,fluid3,x);%BIG's (Boiled Inside Gas)composition is y0，updated liquid pahse's composition is x1
[Q,TL1,DL1,Dv0]=refpropm('QT+-','P',P1,'H',hl1,fluid1,fluid2,fluid3,x);%previous liquid phase convert to two phase 
hv0=refpropm('H','P',P1,'Q',1,fluid1,fluid2,fluid3,y0);%BIG's enthapy
BIG=Q*ML;% 
ML1=ML-BIG;%liquid mass for next time
VL1=ML1/DL1;
VV1=V-VL1;
hv1=(BIG*hv0+MV*hv+QVin-QVL)/(BIG+MV);%ethapy for next timepoint
y1=(MV*y+BIG*y0)/(MV+BIG);
[Dv1,TV1]=refpropm('DT','P',P1,'H',hv1,fluid1,fluid2,fluid3,y1);
M1=Dv1*VV1+DL1*VL1;
BOG=ML+MV-M1;
end