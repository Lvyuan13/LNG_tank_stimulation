function h0=geo_RLVltoh0(R,L,VL)
%for horizontal cylinder
%R,L,Vl to height
% R=0.265;
% L=1.7;
% VL=0.25
a=0;
b=2*R;
i=0;
tolr=1e-10;
itmax=100;

fpa=VL-L*((pi-acos((a-R)/R))*R^2+(a-R)*R*sin(acos((a-R)/R))); 
fpb=VL-L*((pi-acos((b-R)/R))*R^2+(b-R)*R*sin(acos((b-R)/R)));
e=b-a;
while(e>tolr)
    c=(a+b)/2;
    fpa=VL-L*((pi-acos((a-R)/R))*R^2+(a-R)*R*sin(acos((a-R)/R))); 
    fpb=VL-L*((pi-acos((b-R)/R))*R^2+(b-R)*R*sin(acos((b-R)/R)));
    fpc=VL-L*((pi-acos((c-R)/R))*R^2+(c-R)*R*sin(acos((c-R)/R)));
    if fpa*fpc<0
        b=c;
    else
        a=c;
    end
e=b-a;
i=i+1;
if i>itmax
    disp('error')
    break
end
end
h0=c;

end

    