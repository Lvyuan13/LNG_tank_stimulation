function [P1,TL1,TV1,BOG,x1,y1,VL1]=runstim(TL,x,TV,y,VL,Passume)
global fluid1 fluid2 fluid3 
global R L dt V Tair Kwet Kdry s
%% ����������
%  �±�Ϊ1������״̬
%  �����±�Ϊ��״̬
%  Pasume1 Pasume2 
%  Dm/Dp ����΢�� x0=X-f(X)/h(X); f(x)Ϊ���� h(x)Ϊ����
% f(x)=fstep(P)-BOG0
%% ����ǰһ��ʱ��ĸ�����״̬���趨�������õ�ϵ��
%P=refpropm('P','T',TL,'Q',0,fluid1,fluid2,fluid3,x);
[P,DL,hl]=refpropm('PDH','T',TL,'Q',0,fluid1,fluid2,fluid3,x);%Һ�����
%����Ϊ����:Q Q>1 Ϊ�������Ǽ��������
%Q=refpropm('Q','T',TV,'P',P,fluid1,fluid2,fluid3,y)
try
%[Dv,Viscv]=refpropm('DV','T',TV,'P',P,fluid1,fluid2,fluid3,y);
[kv,Dv,beta,Viscv,Cpv,hv] = refpropm('LDBVCH','T',TV,'P',P,fluid1,fluid2,fluid3,y);%�������
catch
[Dv,hv]=refpropm('DH','T',TV,'P',P,fluid1,fluid2,fluid3,y);
Viscv=6.1033e-06;
kv=0.0136;
beta=0.0093;
Cpv=2.0412e+03;
end

ML=DL*VL;%��һʱ��Һ������
VV=V-VL;
MV=Dv*VV; %��һʱ����������
%�˲���Ԥ�ϵķ��������Ƿ�Ϊ���࣬�����������࣬Viscv��Ȼ����
H=geo_RLVltoh0(R,L,VL);%�ɸ�����һʱ�����������߶�
Avl=2*L*sqrt(R^2-(R-H)^2);%���ಿ�������ĽӴ��������
Awet=2*R*((L+R)*(pi-acos((H-R)/R))+(H-R)*sin(acos((H-R)/R)));
Adry=2*pi*R^2+2*pi*R*L-2*R*((L+R)*(pi-acos((H-R)/R))+(H-R)*sin(acos((H-R)/R)));
m0=flow_TPxP0toFR(P,s,Dv,Viscv);%����ģ���Ƴ�й©���ʣ�����һʱ�̵�����Ϊй©��
BOG0=m0*dt;                %���շ���ģ�ͼ�����ĵ�λ������ʱ���ڵ�й¶��
eps=1e-7;       %eps ��ȡֵ�ƺ�����Ӱ�쵽���㣬Ӧ������dll�ļ��㾫������
tolr=1e-7;
%�ó�������ʱΪ��ֵ
itmax=80;
Passume1=Passume;
%% ����ѭ��
for i=1:itmax
   %һ��ţ�ٵ�����
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
%% ��ȡֵ ��ȡֵ
P1=Passume2;
%����fstep���Ƴ����ķ���Ҫ���Passume2

g=9.8;%�������ٶ�
h=6.7296588*kv^(0.857144)*(g*Dv^2*beta*(TV-TL)*Cpv/Viscv)^0.142856; %��Һ�����������ϵ��
QVL=Avl*h*(TV - TL)*dt; %����һʱ�̵���һʱ��ʱ�������Һ�໻������
QLin=Kwet*Awet*(Tair-TL)*dt;
QVin=Kdry*Adry*(Tair-TV)*dt;
%����һʱ�̵���һʱ��ʱ�����������Һ�������Ļ�����
hl1=hl+(QLin+QVL)/ML;  %��һʱ��Һ����ʣ������˴��ȣ�,��Һ�岿�ֿ���һ��������
[x1,y0]=refpropm('X','P',P1,'H',hl1,fluid1,fluid2,fluid3,x);%BIG�����Ϊy0���µ�Һ������Ϊx1
[Q,TL1,DL1,Dv0]=refpropm('QT+-','P',P1,'H',hl1,fluid1,fluid2,fluid3,x);%��������˴��Ⱥ��Һ��ֻ��ɵ���̬
hv0=refpropm('H','P',P1,'Q',1,fluid1,fluid2,fluid3,y0);%BIG����
BIG=Q*ML;%עBIG�Ѿ����ܵ��������������ٳ���dt
ML1=ML-BIG;%�Ƴ���һʱ���µ�Һ��ֵ
VL1=ML1/DL1;
VV1=V-VL1;
hv1=(BIG*hv0+MV*hv+QVin-QVL)/(BIG+MV);%�����˴��Ⱥͽ���BIG�����һʱ��������
%���Ⱥ��������
y1=(MV*y+BIG*y0)/(MV+BIG);
[Dv1,TV1]=refpropm('DT','P',P1,'H',hv1,fluid1,fluid2,fluid3,y1);
M1=Dv1*VV1+DL1*VL1;
BOG=ML+MV-M1;
end