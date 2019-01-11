function m0=flow_TPxP0toFR(p,s,Dv,Viscv)
theta=2*pi/3;       %cone angle of teh valve seat
d0=0.005;      %orifice diameter
D=0.01;         %ball diameter
x0=0;             %initial opening
hmax=0.85;          %maximum valve opening
or=1;
Cd=0.65;          %flow discharge coefficient
Recr=10;          %critical Reynolds number
Aleak=1e-12;        %closed valve leakage area
Amax=1/4*pi*d0^2;     %maximum valve open area
p0=101.325;       %kPa
dp=(p-p0)*1000;       %pressure differential

h=x0+s*or;     %valve opening
if h<=0
    Ah=Aleak;
elseif (0<h)&&(h<hmax)
    Ah=pi*cos(theta/2)*sin(theta/2)*h*(D+sin(theta/2)*h);
else
    Ah=Amax+Aleak;
end
DH=sqrt(4*Ah/pi);   % valve instantaneous hydraulic diameter
%Dv=ref_TPtoDv(temp,p,z); %fluid density
%Viscv=ref_TPtoViscv(temp,p,z);  %fluid kinematic viscosity
Pcr=Dv/(2*(Recr*Viscv/(Cd*DH))^2);    %minimum pressure for turbulent folw
m0=Cd*Ah*sqrt(2/Dv)*dp/(dp^2+Pcr^2)^0.25;
end


