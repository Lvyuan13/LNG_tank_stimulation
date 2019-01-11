function BOG=fstep(TL,x,TV,y,VL,Passume)
%the five inputs determine the state
%describtion:with state values of last time-point plus a assumed pressure to caculate boiled off gass

global fluid1 fluid2 fluid3
global R L dt V Tair Kwet Kdry
%backward differnce method
[P,DL,hl]=refpropm('PDH','T',TL,'Q',0,fluid1,fluid2,fluid3,x);
%determine the pressure density enthapy of liquid zone,which is assume
%to montain saturated
H=geo_RLVltoh0(R,L,VL);
% with VL presented, calculate the liquid height of last timepoint
Avl=2*L*sqrt(R^2-(R-H)^2);
% area of Heat transfer between gas phase and outside envrionment
Awet=2*R*((L+R)*(pi-acos((H-R)/R))+(H-R)*sin(acos((H-R)/R)));
Adry=2*pi*R^2+2*pi*R*L-2*R*((L+R)*(pi-acos((H-R)/R))+(H-R)*sin(acos((H-R)/R)));
g=9.8;%Gravity acceleration
try
[kv,Dv,beta,Viscv,Cpv,hv] = refpropm('LDBVCH','T',TV,'P',P,fluid1,fluid2,fluid3,y);
catch
[Dv,hv]=refpropm('DH','T',TV,'P',P,fluid1,fluid2,fluid3,y);
Viscv=6.1033e-06;
kv=0.0136;
beta=0.0093;
Cpv=2.0412e+03;
end
%The physical quantities can be found from refprop dll, and the heat 
%transfer characteristics of gas phase can be obtained. 
%The purpose is to determine the heat transfer coefficient of gas-liquid two-phaseflow.
ML=DL*VL;%Last moment liquid mass
VV=V-VL;
MV=Dv*VV; %Last moment gas phase mass
h=6.7296588*kv^(0.857144)*(g*Dv^2*beta*(TV-TL)*Cpv/Viscv)^0.142856; %Heat Transfer Coefficient of Gas-liquid 
%Two-phase Flow
QVL=Avl*h*(TV - TL)*dt; %Total Heat Transfer between Gas and Liquid Phase from the Last Moment to the Next Moment
QLin=Kwet*Awet*(Tair-TL)*dt;
QVin=Kdry*Adry*(Tair-TV)*dt;
%Heat transfer of gas phaseVSenvironment and liquid phaseVSenvironment from the last moment to the next moment
hl1=hl+(QLin+QVL)/ML;  
%The enthalpy of the liquid at the next moment (after heat transfer process);takes the liquid part as a control Volume.

[x1,y0]=refpropm('X','P',Passume,'H',hl1,fluid1,fluid2,fluid3,x);
%The component of BIG is y0, and the component of new liquid phase is x1.
[Q,TL1,DL1,Dv0]=refpropm('QT+-','P',Passume,'H',hl1,fluid1,fluid2,fluid3,x);
%the new state of liquid phase's differentiation after heat transfer is obtained.
%%
%Q     %just test if 0<Q<1 
%% differnt conditions when Q differs from refprop64.dll：
%          q--vapor quality on a MOLAR basis [moles vapor/total moles]
%          q < 0 indicates subcooled (compressed) liquid
%          q = 0 indicates saturated liquid
%          q = 1 indicates saturated vapor
%          q > 1 indicates superheated vapor
%          q = 998 superheated vapor, but quality not defined (in most situations, t > Tc)
%          q = 999 indicates supercritical state (t > Tc) and (p > Pc)
if Q<=0         % provided no gas evaporate from liquid phase, indicating that liuqid is supcooled
hv1=(MV*hv+QVin-QVL)/MV;% in such condition the enthapy of updated gas phase
[q0,DL0,Dv1]=refpropm('Q+-','P',Passume,'H',hv1,fluid1,fluid2,fluid3,y);% DL0 condensed liuqid density Dv1 evaporated gas density
% differnt conditions when there is condensed liuqid from gas phase or not
%q0     %just for test q>1 is valid
   if q0>=1 %gas is supheated or saturated
      %no change of liuqid phase volume
      BOG=MV-Dv1*VV;
   else %there is condensed liuqid from gas phase
      e=1-q0;%condensing ratio of gas phase
      MV0=MV*q0;%new gasphase of gas after condensing
      ML0=MV*e;% mass of condensed liquid
      VV1=VV-ML0/DL0;
      BOG=MV-MV0-VV1*Dv1;
   end
else     % there is gas evaporated from liquid, valid
    %假设下压力的液体分化出气相，属于猜测范围内
hv0=refpropm('H','P',Passume,'Q',1,fluid1,fluid2,fluid3,y0);%BIG的焓
BIG=Q*ML;%注BIG已经是总的蒸发量，不用再乘以dt
ML1=ML-BIG;%推出下一时刻新的液相值
VL1=ML1/DL1;
VV1=V-VL1;
hv1=(BIG*hv0+MV*hv+QVin-QVL)/(BIG+MV);%经历了传热和进入BIG后的下一时刻气相焓
%求传热后的总质量
y1=(MV*y+BIG*y0)/(MV+BIG);
[Dv1,Tv1]=refpropm('DT','P',Passume,'H',hv1,fluid1,fluid2,fluid3,y1);
%QLin QVin QVL 的量级需要调整
%BOG是此段时间内的总量
M1=Dv1*VV1+DL1*VL1;
BOG=ML+MV-M1;
end
   
%% 进行求取气相

end
