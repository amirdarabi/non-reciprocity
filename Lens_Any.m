% This fucntion different lenses
function [si,Rank]=Lens_Any
x=-1;
y=0;
F0=1;
BC=1;
omega=50000*2*pi;
size=400;
row=55;%130;
S=5;
%%
xp=zeros(1,size);
yp=zeros(1,size);
% % EXP 
for i=1:size
    xp(i)=(-.05+.8*(i-1)/(size-1));
    yp(i)=(.1*(i-1)/(size-1));%(.1*(i-1)/(size-1)-.05);
end
[Bp,Be,k,d,N,diaa,Rank,si]=GRIN_Leamy(F0,omega,x,BC,row,S);
w_sca=zeros(size,size);
[X,Y] = meshgrid(xp,yp);
Min_b=(S-1)/2;
for iii=1:N
    for jjj=-Min_b:Min_b
                rr=((Y-d(iii,2)).^2+(X-d(iii,1)).^2).^.5;
                theta1=atan2(Y-d(iii,2),X-d(iii,1));
                w_sca=Be(iii,jjj+Min_b+1)*besselk(jjj,k*rr).*exp(1j*jjj*theta1)+Bp(iii,jjj+Min_b+1)*besselh(jjj,1,k*rr).*exp(1j*jjj*theta1)+w_sca;
    end
end

w_inc=F0*exp(1j*k*(X-x));
% w_tot=w_sca;%+w_inc;
%% hole disp
for iiii=1:size
    for jjjj=1:size
        for iii=1:N
            kkll=floor((iii-1)/row)+1;
                rr=sqrt((yp(jjjj)-d(iii,2))^2+(xp(iiii)-d(iii,1))^2);
                if rr<=diaa(kkll)
                w_sca(jjjj,iiii)=NaN;
                end
        end
    end
    iiii
end
%%
T=0.2;
w_tot=abs(real(w_sca*exp(1j*omega*T)));
w=abs(w_sca);
figure
contourf(X,Y,w)
hold on
contourf(X,-Y,w)
figure
contour(X,Y,w)
hold on
contour(X,-Y,w)
figure
mesh(X,Y,abs(w))
hold on
mesh(X,-Y,abs(w))
for i=1:size
    xp(i)=(-.05+.8*(i-1)/(size-1));
    yp(i)=(.1*(i-1)/(size-1));%(.1*(i-1)/(size-1)-.05);
end
[Bp,Be,k,d,N,diaa,Rank,si]=GRIN_Erturk_Case(F0,omega,x,row,S);
w_sca=zeros(size,size);
[X,Y] = meshgrid(xp,yp);
Min_b=(S-1)/2;
for iii=1:N
    for jjj=-Min_b:Min_b
                rr=((Y-d(iii,2)).^2+(X-d(iii,1)).^2).^.5;
                theta1=atan2(Y-d(iii,2),X-d(iii,1));
                w_sca=Be(iii,jjj+Min_b+1)*besselk(jjj,k*rr).*exp(1j*jjj*theta1)+Bp(iii,jjj+Min_b+1)*besselh(jjj,1,k*rr).*exp(1j*jjj*theta1)+w_sca;
    end
end

w_inc=F0*exp(1j*k*(X-x));
% w_tot=w_sca;%+w_inc;
%% hole disp
for iiii=1:size
    for jjjj=1:size
        for iii=1:N
            kkll=floor((iii-1)/row)+1;
                rr=sqrt((yp(jjjj)-d(iii,2))^2+(xp(iiii)-d(iii,1))^2);
                if rr<=diaa(kkll)
                w_sca(jjjj,iiii)=NaN;
                end
        end
    end
    iiii
end
%%
T=0.2;
w_tot=abs(real(w_sca*exp(1j*omega*T)));
w=abs(w_sca);
figure
contourf(X,Y,w)
hold on
contourf(X,-Y,w)
figure
contour(X,Y,w)
hold on
contour(X,-Y,w)
figure
mesh(X,Y,abs(w))
hold on
mesh(X,-Y,abs(w))

end