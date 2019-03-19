%This fucntion does GRIN for 75 scatters with new method
function [Bp,Be,k,d,N,diaa,Rank,DD,si]=GRIN_Leamy(F0,omega,x,BC,row,S)
type=11;
% digits(2);
N=type*row;
h=0.003175;
Ll2=0.0055;
d=zeros(N,3);
rho=2800;
E=71.1e9;
nu=0.34;
Ll=0.008;
dia1=0.00545/2;
dia2=0.0054/2;
dia3=0.00535/2;
dia4=0.00525/2;
dia5=0.0045/2;
dia6=0.0039/2;
diaa(11)=dia6;
diaa(10)=dia5;
diaa(9)=dia4;
diaa(8)=dia3;
diaa(7)=dia2;
diaa(6)=dia1;
diaa(5)=dia2;
diaa(4)=dia3;
diaa(3)=dia4;
diaa(2)=dia5;
diaa(1)=dia6;
DD=E*h^3/12/( 1-nu^2);
Min_b=(S-1)/2;
k=(rho*h/DD)^.25*omega^.5;
% %%
% x,y location of scatters; 3 is the diameter
for i=1:row     
    d(i,1)=Ll*(i-1);
    d(i,2)=5*Ll2;
    d(i,3)=dia6;
    d(i+row,1)=Ll*(i-1);
    d(i+row,2)=4*Ll2;
    d(i+row,3)=dia5;
    d(i+2*row,1)=Ll*(i-1);
    d(i+2*row,2)=3*Ll2;
    d(i+2*row,3)=dia4;
    d(i+3*row,1)=Ll*(i-1);
    d(i+3*row,2)=2*Ll2;
    d(i+3*row,3)=dia3;
    d(i+4*row,1)=Ll*(i-1);
    d(i+4*row,2)=1*Ll2;
    d(i+4*row,3)=dia2;
    d(i+5*row,1)=Ll*(i-1);
    d(i+5*row,2)=0*Ll2;
    d(i+5*row,3)=dia1;
    d(i+6*row,1)=Ll*(i-1);
    d(i+6*row,2)=-1*Ll2;
    d(i+6*row,3)=dia2;
    d(i+7*row,1)=Ll*(i-1);
    d(i+7*row,2)=-2*Ll2;
    d(i+7*row,3)=dia3;
    d(i+8*row,1)=Ll*(i-1);
    d(i+8*row,2)=-3*Ll2;
    d(i+8*row,3)=dia4;
    d(i+9*row,1)=Ll*(i-1);
    d(i+9*row,2)=-4*Ll2;
    d(i+9*row,3)=dia5;
    d(i+10*row,1)=Ll*(i-1);
    d(i+10*row,2)=-5*Ll2;
    d(i+10*row,3)=dia6;
end
figure
xxx=d(:,1);
yyy=d(:,2);
zzz=d(:,3);
plot3(xxx,yyy,zzz)
% distance between source and each scatterer
for ij=1:N
        dis_s(ij,1)=abs(d(ij,1)-x);
end
%%
% Incident plane wave coefficients
for uii=-Min_b:Min_b
    for uiu=1:N
        Ap(uii+Min_b+1,uiu)=F0*exp(1j*k*dis_s(uiu))*(1j)^(uii);
        Ae(uii+Min_b+1,uiu)=0;
    end
end
%%
% T-matrix initialization
Tpp=zeros(S,S,type);
Tpe=zeros(S,S,type);
Tep=zeros(S,S,type);
Tee=zeros(S,S,type);
%%
% Element of each matrix
chi21=zeros(S,type);
chi22=zeros(S,1);
chi11=zeros(S,type);
chi12=zeros(S,1);
chi31=zeros(S,type);
chi32=zeros(S,1);
chi41=zeros(S,type);
chi42=zeros(S,1);
delta=zeros(S,1);
%%
for rfr=1:type
if BC==1% void inclusions
    for iii=-Min_b:Min_b
    z=k*diaa(rfr);
    chi11(iii+Min_b+1)=DD*(-(1-nu)*z*(iii*besselj(iii,z)/z-besselj(iii+1,z))+((1-nu)*iii^2-z^2)*besselj(iii,z));
    chi12(iii+Min_b+1)=DD*(-((1-nu)*iii^2+z^2)*z*(iii*besselj(iii,z)/z-besselj(iii+1,z))+(1-nu)*iii^2*besselj(iii,z));
    chi21(iii+Min_b+1)=DD*(-(1-nu)*z*(iii*besselh(iii,1,z)/z-besselh(iii+1,1,z))+((1-nu)*iii^2-z^2)*besselh(iii,1,z));
    chi22(iii+Min_b+1)=DD*(-((1-nu)*iii^2+z^2)*z*(iii*besselh(iii,1,z)/z-besselh(iii+1,1,z))+(1-nu)*iii^2*besselh(iii,1,z));
    chi31(iii+Min_b+1)=DD*(-(1-nu)*z*(iii*besseli(iii,z)/z+besseli(iii+1,z))+((1-nu)*iii^2+z^2)*besseli(iii,z));
    chi32(iii+Min_b+1)=DD*(-((1-nu)*iii^2-z^2)*z*(iii*besseli(iii,z)/z+besseli(iii+1,z))+(1-nu)*iii^2*besseli(iii,z));
    chi41(iii+Min_b+1)=DD*(-(1-nu)*z*(iii*besselk(iii,z)/z-besselk(iii+1,z))+((1-nu)*iii^2+z^2)*besselk(iii,z));
    chi42(iii+Min_b+1)=DD*(-((1-nu)*iii^2-z^2)*z*(iii*besselk(iii,z)/z-besselk(iii+1,z))+(1-nu)*iii^2*besselk(iii,z));
    delta(iii+Min_b+1)=chi21(iii+Min_b+1)*chi42(iii+Min_b+1)-chi41(iii+Min_b+1)*chi22(iii+Min_b+1);
    Tpp(iii+Min_b+1,iii+Min_b+1,rfr)=(chi12(iii+Min_b+1)*chi41(iii+Min_b+1)-chi11(iii+Min_b+1)*chi42(iii+Min_b+1))/delta(iii+Min_b+1);
    Tep(iii+Min_b+1,iii+Min_b+1,rfr)=(chi22(iii+Min_b+1)*chi11(iii+Min_b+1)-chi21(iii+Min_b+1)*chi12(iii+Min_b+1))/delta(iii+Min_b+1);
    Tpe(iii+Min_b+1,iii+Min_b+1,rfr)=(chi32(iii+Min_b+1)*chi41(iii+Min_b+1)-chi31(iii+Min_b+1)*chi42(iii+Min_b+1))/delta(iii+Min_b+1);
    Tee(iii+Min_b+1,iii+Min_b+1,rfr)=(chi22(iii+Min_b+1)*chi31(iii+Min_b+1)-chi21(iii+Min_b+1)*chi32(iii+Min_b+1))/delta(iii+Min_b+1);
    end
    elseif BC==2 %rigid inclusions
        for iii=-Min_b:Min_b
            z=k*diaa(rfr);
%             delta(iii+Min_b+1)=besselh(iii+1,1,k*diaa(rfr))*besselk(iii,k*diaa(rfr))-besselh(iii,1,k*diaa(rfr))*besselk(iii+1,k*diaa(rfr));
%             Tpp(iii+Min_b+1,iii+Min_b+1,rfr)=(besselj(iii,k*diaa(rfr))*besselk(iii+1,k*diaa(rfr))-besselj(iii+1,k*diaa(rfr))*besselk(iii,k*diaa(rfr)))/delta(iii+Min_b+1);
%             Tep(iii+Min_b+1,iii+Min_b+1,rfr)=2*1j/pi/k/diaa(rfr)/delta(iii+Min_b+1);%(besselj(iii,k*diaa(rfr))*(besselh(iii,1,k*diaa(rfr))/k/diaa(rfr)*iii-besselh(iii+1,1,k*diaa(rfr)))-(besselj(iii-1,k*diaa(rfr))-besselj(iii+1,k*diaa(rfr)))/2*besselh(iii,1,k*diaa(rfr)))/delta(iii+Min_b+1);
%             Tpe(iii+Min_b+1,iii+Min_b+1,rfr)=1/k/diaa(rfr)/delta(iii+Min_b+1);%(besseli(iii,k*diaa(rfr))/2*(besselk(iii-1,k*diaa(rfr))+besselk(iii+1,k*diaa(rfr)))+(besseli(iii-1,k*diaa(rfr))+besseli(iii+1,k*diaa(rfr)))/2*besselk(iii,k*diaa(rfr)))/delta(iii+Min_b+1);
%             Tee(iii+Min_b+1,iii+Min_b+1,rfr)=-(besseli(iii,k*diaa(rfr))*besselh(iii+1,1,k*diaa(rfr))+besseli(iii+1,k*diaa(rfr))*besselh(iii,1,k*diaa(rfr)))/delta(iii+Min_b+1);
            deltaa=(besselj(iii,z)+1j*bessely(iii,z))*diff(besselk(iii,z),z)-diff(besselj(iii,z)+1j*bessely(iii,z),z)*besselk(iii,z);
            first=(-besselj(iii,z)*diff(besselk(iii,z),z)+diff(besselj(iii,z),z)*besselk(iii,z))/deltaa;
            second=2*1j/pi/z/deltaa;%(besselj(iii,k*diaa(rfr))*(besselh(iii,1,k*diaa(rfr))/k/diaa(rfr)*iii-besselh(iii+1,1,k*diaa(rfr)))-(besselj(iii-1,k*diaa(rfr))-besselj(iii+1,k*diaa(rfr)))/2*besselh(iii,1,k*diaa(rfr)))/delta(iii+Min_b+1);
            third=1/z/deltaa;%(besseli(iii,k*diaa(rfr))/2*(besselk(iii-1,k*diaa(rfr))+besselk(iii+1,k*diaa(rfr)))+(besseli(iii-1,k*diaa(rfr))+besseli(iii+1,k*diaa(rfr)))/2*besselk(iii,k*diaa(rfr)))/delta(iii+Min_b+1);
            fourth=(besseli(iii,z)*diff(besselj(iii,z)+1j*bessely(iii,z),z)-diff(besseli(iii,z),z)*(besselj(iii,z)+1j*bessely(iii,z)))/deltaa;
            Tpp(iii+Min_b+1,iii+Min_b+1,rfr)=subs(first);
            Tep(iii+Min_b+1,iii+Min_b+1,rfr)=subs(second);
            Tpe(iii+Min_b+1,iii+Min_b+1,rfr)=subs(third);
            Tee(iii+Min_b+1,iii+Min_b+1,rfr)=subs(fourth);
            
        end
end
end

%%
% Coordinate transformations for each pair of scatterers
for iii=1:N
    for jj=1:N
        for kk=-Min_b:Min_b
            for ll=-Min_b:Min_b
                if iii==jj    
                Rp(ll+Min_b+1,kk+Min_b+1,iii,jj)=0;
                Re(ll+Min_b+1,kk+Min_b+1,iii,jj)=0;
                else    
                Rp(ll+Min_b+1,kk+Min_b+1,iii,jj)=exp(1j*(kk-ll)*(atan2(d(jj,2)-d(iii,2),d(jj,1)-d(iii,1))))*besselh(kk-ll,1,k*(sqrt((d(iii,1)-d(jj,1))^2+((d(iii,2)-d(jj,2))^2))));
                Re(ll+Min_b+1,kk+Min_b+1,iii,jj)=((-1)*ll)*exp(1j*(kk-ll)*(atan2(d(jj,2)-d(iii,2),d(jj,1)-d(iii,1))))*besselk(kk-ll,k*(sqrt((d(iii,1)-d(jj,1))^2+((d(iii,2)-d(jj,2))^2))));
                end
            end
        end
    end
end

% L-matrix, A-matrix as in paper
L=zeros(2*S*N,2*S*N);
AA=zeros(2*S*N,1);
for ii=1:2*S*N
     for jj=1:2*S*N
     kk=floor((ii-1)/(2*S))+1;
     ll=floor((jj-1)/(2*S))+1;
     kkll=floor((kk-1)/row)+1;
     if ii==jj
         L(ii,jj)=1;
     else
             kkk=ii-(2*S)*(kk-1);
             lll=jj-(2*S)*(ll-1);
             if kkk <= S && lll<=S
                 L(ii,jj)=-(Tpp(kkk,:,kkll)*Rp(:,lll,ll,kk));
             elseif kkk <=S && lll>S
                 L(ii,jj)=-(Tpe(kkk,:,kkll)*Re(:,lll-S,ll,kk));
             elseif kkk>S && lll<= S
                 L(ii,jj)=-(Tep(kkk-S,:,kkll)*Rp(:,lll,ll,kk));
             else
                 L(ii,jj)=-(Tee(kkk-S,:,kkll)*Re(:,lll-S,ll,kk));
             end
     end
     end
    if kkk <= S 
        AA(ii,1)=(Tpp(kkk,:,kkll)*Ap(:,kk)+Tpe(kkk,:,kkll)*Ae(:,kk));
    else
        AA(ii,1)=(Tep(kkk-S,:,kkll)*Ap(:,kk)+Tee(kkk-S,:,kkll)*Ae(:,kk));
    end
end
si=size(L)
BB=L\AA;
Rank=rank(L)
Be=zeros(N,S);
Bp=zeros(N,S);
for iiii=1:N
    for jjii=1:S
        Bp(iiii,jjii)=BB((iiii-1)*2*S+jjii);
        Be(iiii,jjii)=BB((iiii-1)*2*S+jjii+S);
    end
end
end
%%Multiple_scattering_Any(-300,0)
