function [P1,TL1,TV1,BOG,x1,y1,VL1]=runstim(TL,x,TV,y,VL,Passume)
global fluid1 fluid2 fluid3 
global R L dt V Tair Kwet Kdry s
%% 变量名规则
%  下标为1代表新状态
%  不带下标为旧状态
%  Pasume1 Pasume2 
%  Dm/Dp 代表微分 x0=X-f(X)/h(X); f(x)为函数 h(x)为导数
% f(x)=fstep(P)-BOG0
%% 完善前一个时间的各参数状态及设定迭代所用的系数
%P=refpropm('P','T',TL,'Q',0,fluid1,fluid2,fluid3,x);
[P,DL,hl]=refpropm('PDH','T',TL,'Q',0,fluid1,fluid2,fluid3,x);%液相参数
%仅仅为测试:Q Q>1 为符合我们假设的条件
%Q=refpropm('Q','T',TV,'P',P,fluid1,fluid2,fluid3,y)
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

ML=DL*VL;%上一时刻液相质量
VV=V-VL;
MV=Dv*VV; %上一时刻气相质量
%此步可预料的风险在于是否为气相，若果处于两相，Viscv必然报错
H=geo_RLVltoh0(R,L,VL);%由给定上一时刻体积，求出高度
Avl=2*L*sqrt(R^2-(R-H)^2);%气相部分与外界的接触传热面积
Awet=2*R*((L+R)*(pi-acos((H-R)/R))+(H-R)*sin(acos((H-R)/R)));
Adry=2*pi*R^2+2*pi*R*L-2*R*((L+R)*(pi-acos((H-R)/R))+(H-R)*sin(acos((H-R)/R)));
m0=flow_TPxP0toFR(P,s,Dv,Viscv);%阀门模型推出泄漏速率，以上一时刻的速率为泄漏率
BOG0=m0*dt;                %按照阀门模型计算出的单位个步长时间内的泄露量
eps=1e-7;       %eps 的取值似乎可以影响到计算，应该是由dll的计算精度限制
tolr=1e-7;
%该程序运行时为定值
itmax=80;
Passume1=Passume;
%% 开启循环
for i=1:itmax
   %一阶牛顿迭代法
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
%% 存取值 求取值
P1=Passume2;
%按照fstep中推出出的符合要求的Passume2

g=9.8;%重力加速度
h=6.7296588*kv^(0.857144)*(g*Dv^2*beta*(TV-TL)*Cpv/Viscv)^0.142856; %气液两相对流换热系数
QVL=Avl*h*(TV - TL)*dt; %从上一时刻到下一时刻时间段里气液相换热总量
QLin=Kwet*Awet*(Tair-TL)*dt;
QVin=Kdry*Adry*(Tair-TV)*dt;
%从上一时刻到下一时刻时间段里气相与液相与外界的换热量
hl1=hl+(QLin+QVL)/ML;  %下一时刻液体的焓（经历了传热）,将液体部分看作一个控制体
[x1,y0]=refpropm('X','P',P1,'H',hl1,fluid1,fluid2,fluid3,x);%BIG的组分为y0，新的液相的组分为x1
[Q,TL1,DL1,Dv0]=refpropm('QT+-','P',P1,'H',hl1,fluid1,fluid2,fluid3,x);%求出经历了传热后的液相分化成的新态
hv0=refpropm('H','P',P1,'Q',1,fluid1,fluid2,fluid3,y0);%BIG的焓
BIG=Q*ML;%注BIG已经是总的蒸发量，不用再乘以dt
ML1=ML-BIG;%推出下一时刻新的液相值
VL1=ML1/DL1;
VV1=V-VL1;
hv1=(BIG*hv0+MV*hv+QVin-QVL)/(BIG+MV);%经历了传热和进入BIG后的下一时刻气相焓
%求传热后的总质量
y1=(MV*y+BIG*y0)/(MV+BIG);
[Dv1,TV1]=refpropm('DT','P',P1,'H',hv1,fluid1,fluid2,fluid3,y1);
M1=Dv1*VV1+DL1*VL1;
BOG=ML+MV-M1;
end