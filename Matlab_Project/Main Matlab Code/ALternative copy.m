
Lc=0.116;                 
Wc=0.116; 
Nc=1; 
Ac=Lc*Wc*Nc; 
Asc=0.003375; 
Wgap1=0.01;
Wgap2=0.01;
Wgap3=0.01;
Wgap4=0.005;


cell='silicon';
FLdFlt='valvoline';
kfl=0.1384;
GG=[200 400 600 800 1000];
Ta=308.0;
Twtin=298.0;
mwt=[0.0079 0.008 0.009 0.01 0.011 0.012 0.013 0.014 0.015]*Ac;




Wflg=0.0007;
Wgc=0.0038;

Emgc=0.9; 
[Trgc,Relgc,Abgc]=select_filter_M('anti_reflective_thin_glass_3.8mm');

[Trwt,Relwt,Abwt]=select_filter_M('pure_water_10mm');

Emgfl=0.9; 
[Trflg,Relflg,Abflg]=select_filter_M('anti_reflective_thin_glass_0.7mm');

if strcmp(FLdFlt,'nanofluid_1')==1
    [Trflf,Relflf,Abflf]=select_filter_M('nanofluid_1');
end

if strcmp(FLdFlt,'nanofluid_2')==1
    [Trflf,Relflf,Abflf]=select_filter_M('nanofluid_2');
end

if strcmp(FLdFlt,'nanofluid_3')==1
    [Trflf,Relflf,Abflf]=select_filter_M('nanofluid_3');
end


if strcmp(FLdFlt,'valvoline')==1
    [Trflf,Relflf,Abflf]=select_filter_M('valvoline');
end

if strcmp(FLdFlt,'glycerol')==1
    [Trflf,Relflf,Abflf]=select_filter_M('glycerol');
end

if strcmp(FLdFlt,'pure_water_10mm')==1
    [Trflf,Relflf,Abflf]=select_filter_M('pure_water_10mm');
end

Abpv=0.88; 
Empv=0.75*ones(2002,1); 


AM=xlsread('data_M', 'ambient conditions', 'C2:C2003'); 
GAM=1000.37;

Vwind=0;

Tsky=0.0552*Ta^1.5;

cpwt=4180;
Tss=329;
Emss=0.35;

Twt=Twtin;
Nm=size(mwt);
NG=size(GG);


for j=5
G=GG(j);    
G0=G/GAM*AM;

[G1pass,G1rel,G1abs]=Gfilter(G0,Trgc,Relgc,Abgc);
[G2pass,G2rel,G2abs]=Gfilter(G1pass,Trwt,Relwt,Abwt);
[G3pass,G3rel,G3abs]=Gfilter(G2pass,Trflg,Relflg,Abflg);
[G4pass,G4rel,G4abs]=Gfilter(G3pass,Trflf,Relflf,Abflf);
[G5pass,G5rel,G5abs]=Gfilter(G4pass,Trflg,Relflg,Abflg);
[G6pass,G6rel,G6abs]=Gfilter(G5pass,Trwt,Relwt,Abwt);
[G7pass,G7rel,G7abs]=Gfilter(G6pass,Trgc,Relgc,Abgc);

[G7pass1,G7rel1,G7abs1]=Gfilter(G7pass*(1-Abpv),Trgc,Relgc,Abgc);
G7abs=G7abs+G7abs1;

[G6pass1,G6rel1,G6abs1]=Gfilter(G7rel+G7pass1,Trwt,Relwt,Abwt);
G6abs=G6abs+G6abs1;

[G5pass1,G5rel1,G5abs1]=Gfilter(G6rel+G6pass1,Trflg,Relflg,Abflg);
G5abs=G5abs+G5abs1;

[G4pass1,G4rel1,G4abs1]=Gfilter(G5rel+G5pass1,Trflf,Relflf,Abflf);
G4abs=G4abs+G4abs1;

[G3pass1,G3rel1,G3abs1]=Gfilter(G4rel+G4pass1,Trflf,Relflf,Abflf);
G3abs=G3abs+G3abs1;

[G2pass1,G2rel1,G2abs1]=Gfilter(G3rel+G3pass1,Trwt,Relwt,Abwt);
G2abs=G2abs+G2abs1;

[G1pass1,G1rel1,G1abs1]=Gfilter(G2rel+G2pass1,Trgc,Relgc,Abgc);
G1abs=G1abs+G1abs1;

Solarcell=solarcell(Asc,G7pass*Abpv,G,cell);%()

for i=1:Nm(2)

Asol=zeros(9,9); 
Bsol=zeros(9,1); 
T=zeros(1,9)+300;
X=zeros(9,1)+300;
Ttemp=9999;
N=1;
error=1e-4;

while (max(abs(X.'-Ttemp))>error && N<100) 
    Ttemp=T;

[Kair,vair,Bair,Prair]=AirProperties(0.5*T(7)+0.5*T(8)); 
hr_gc2ss=0.01*h_radiation(Emgc,Emss,T(1),Tss); 
hwind=2*4.5+2.9*Vwind;
hr_gc2pv=0.01*h_radiation(Empv(1),Emgc,T(8),T(7)); 
hc_gc2pv=1/(2/h_enclosure(T(7),T(8),Wgap4,0)+Wgap4/Kair);
hr_gc22ss=0.01*h_radiation(Emgc,Emss,T(7),Tss); 
hc_wt1=h_watergap(T(2),T(3),0.1*mwt(i),Wgap1,Wc,Lc);
hc_wt2=h_watergap(T(2),T(9),0.9*mwt(i),Wgap3,Wc,Lc); 
hc_pv2amb=1.42*(1/(0.0005/0.35+0.00005/0.85+0.0001/0.2+0.001/205+1/hwind));

hc_filter=1/(Wgap2/2/kfl);

hr_pv2sky=h_radiation(Empv(1),Emss,T(8),Tss);

Asol(1,1)=-Ac*hr_gc2ss-Ac*(hc_wt1)-Ac*hwind;
Asol(1,2)=Ac*(hc_wt1)*0.5;
Asol(1,3)=Ac*(hc_wt1)*0.5;

Bsol(1)=-Ac*hr_gc2ss*Tss-Ac*hwind*Ta-Ac*G1abs;


Asol(2,1)=Ac*hc_wt1;
Asol(2,2)=-Ac*hc_wt1*0.5-Ac*hc_wt1*0.5+0.1*mwt(i)*cpwt;
Asol(2,3)=-Ac*hc_wt1*0.5-Ac*hc_wt1*0.5-0.1*mwt(i)*cpwt;
Asol(2,4)=Ac*hc_wt1;

Bsol(2)=-Ac*G2abs;



    Asol(3,2)=0.5*Ac*(hc_wt1);
    Asol(3,3)=0.5*Ac*(hc_wt1);
    Asol(3,4)=-Ac*hc_wt1-Ac*hc_filter;
    Asol(3,5)=Ac*hc_filter;
    Bsol(3)=-Ac*G2abs;

Asol(4,4)=Ac*hc_filter;
Asol(4,5)=-Ac*hc_filter-Ac*hc_filter;
Asol(4,6)=Ac*hc_filter;

Bsol(4)=-Ac*G4abs;


    Asol(5,2)=0.5*Ac*(hc_wt2);
    Asol(5,5)=Ac*hc_filter;
    Asol(5,6)=-Ac*hc_wt2-Ac*hc_filter;
    Asol(5,9)=0.5*Ac*(hc_wt2);
    Bsol(5)=-Ac*G5abs;

Asol(6,2)=-0.5*Ac*hc_wt2-0.5*Ac*hc_wt2+0.9*mwt(i)*cpwt;
Asol(6,6)=Ac*hc_wt2;
Asol(6,7)=Ac*hc_wt2;
Asol(6,9)=-Ac*hc_wt2*0.5-Ac*hc_wt2*0.5-0.9*mwt(i)*cpwt;

Bsol(6)=-Ac*G6abs;


    Asol(7,2)=Ac*hc_wt2*0.5;
    Asol(7,7)=-Ac*hc_wt2-Ac*hr_gc22ss-Asc*(hr_gc2pv+hc_gc2pv);
    Asol(7,8)=Asc*(hr_gc2pv+hc_gc2pv);
    Asol(7,9)=Ac*hc_wt2*0.5;
    
    Bsol(7)=-Ac*G7abs-Ac*hr_gc22ss*Tss;


effm=Solarcell(1)*(1-Solarcell(6)*(T(8)-273.15-25));%module efficiency
Gpvabs=Solarcell(5)-effm*G;%waste heat in PV

Asol(8,7)=Asc*(hr_gc2pv+hc_gc2pv)-Asc*hr_gc22ss;
Asol(8,8)=-Asc*(hr_gc2pv+hc_gc2pv)-Asc*hc_pv2amb;

Bsol(8)=-Asc*Gpvabs-Asc*hc_pv2amb*Ta-Asc*hr_gc22ss*Tss;

Asol(9,2)=1;
Bsol(9)=Twt;

X=linsolve(Asol(1:9,1:9),Bsol(1:9,1:1));
T=X.';

N=N+1;%count the iteration times
 
%end
end
%efficiencies
%Tr(i)=(0.5*T(3)+0.5*T(4)-Ta)/G;
%effthfl(i)=mwt(i)*cpwt*(0.5*T(3)+0.5*T(9)-T(2))/G/Ac;
T_ave(i)=(0.1*mwt(i)*T(3)+0.9*mwt(i)*T(9))/mwt(i);
effthfl(i)=mwt(i)*cpwt*(T_ave(i)-T(2))/G/Ac;

%effthfl(i)=mwt(i)*cpwt*(T(9)-T(2))/G/Ac;

effel(i)=effm;

%temperatures
Tsspvt(i,1)=T(1);
Tsspvt(i,2)=T(2);
Tsspvt(i,3)=T(3);
Tsspvt(i,4)=T(4);
Tsspvt(i,5)=T(5);
Tsspvt(i,6)=T(6);
Tsspvt(i,7)=T(7);
Tsspvt(i,8)=T(8);
Tsspvt(i,9)=T(9);

%energy balance

% if strcmp(Distance,'far')==1
Energy(i,1)=hr_gc2ss*(T(1)-Tsky)/G;
Energy(i,2)=hwind*(T(1)-Ta)/G;
% Energy(i,3)=hr_pv2sky*(T(8)-Tsky)/G;%far covered
Energy(i,3)=hr_gc2pv*(T(8)-T(7))/G;%close covered
% Energy(i,4)=hwind*(T(8)-Ta)/G;%Far covered
Energy(i,4)=hc_gc2pv*(T(8)-T(7))/G;%close covered
Energy(i,5)=hc_pv2amb*(T(8)-Ta)/G;%%%%%%%%%%%%%%%%CHECKKKKKKKKK
Energy(i,8)=Gpvabs/G;
end



end


[Trgcf,Relgcf,Abgcf]=select_filter_M('fake');

[G1pass2,G1rel2,G1abs2]=Gfilter(G0,Trgcf,Relgcf,Abgcf);%layer1:cover glass,Specrum after the cover glass%%%%CHECKKKKKKKKKK
%Gpv=ones(2002,1);
Solarcell1=solarcell(Asc,G0*Abpv,G,cell);
%simulation initials and settings
Asol1=zeros(2,2); % 2 Equations' left-hand coefficients
Bsol1=zeros(2,1); % 2 equations' right-hand coefficients
T1=zeros(1,2)+300;% T1...T2:Please see the structure of the S-PVT
X1=zeros(2,1)+300;% solution of each iteration
Ttemp1=9999;
N1=1;
error1=1e-4;


 while (max(abs(X1.'-Ttemp1))>error&& N1<100)
     Ttemp1=T1;
%     
% % heat transfer coefficients calculation
hr_pv2ss=h_radiation(Empv,Emss,T1(1),Tss); % radiation from top channel to solar simulator
hwind=4*4.5+2.9*Vwind; %8.3+2.2*vwind(i-1); % 4.5+2.9*vwind, convection heat transfer due to wind


%hc_pvg2pv=1/(0.002/0.8+0.0005/0.35);% heat conduct from pv-glass to pv: through the layers of glass, eva

hc_pv2amb=1/(0.0005/0.35+0.00005/0.85+0.0001/0.2+0.001/205+1/hwind);% heat conduct from pv to abiemt: through the layers of eva, ted, glue, back plate and back convection
hc_pvg2pv=1/(0.00000004/0.8+0.000005/0.35);% heat conduct from pv-glass to pv: through the layers of glass, eva

Asol1(1,1)=-Asc*hr_pv2ss(2)-Asc*hwind-Asc*hc_pvg2pv;
Asol1(1,2)=Asc*hc_pvg2pv;

Bsol1(1)=-Asc*hr_pv2ss(2)*Tss-Asc*hwind*Ta-Asc*G1abs2;


%Equation 2: PV
%Ac*hc_pvg2pv*(Tgpv-Tpv)+Ac*hc_pv2amb*(Ta-Tpv)+Ac*Gpvabs1=0
effm1=Solarcell1(1)*(1-Solarcell1(6)*(T1(2)-273.15-25));%module efficiency
Gpvabs1=Solarcell1(5)-effm1*G;%waste heat in PV


Asol1(2,1)=Asc*hc_pvg2pv;
Asol1(2,2)=-Asc*hc_pvg2pv-Asc*hc_pv2amb;

Bsol1(2)=-Asc*Gpvabs1-Asc*hc_pv2amb*Ta;



%****LINEAR SOLVER****
X1=linsolve(Asol1(1:2,1:2),Bsol1(1:2,1:1));
T1=X1.';

N1=N1+1;%count the iteration times
    
end

effel1=effm1;

 N1=N1+1;%count the iteration times

Tcpvt(1)=T1(1);
Tcpvt(2)=T1(2);

 Energy1(1)=hr_pv2ss(2)*(T1(1)-Tsky)/G;
 Energy1(2)=hwind*(T1(1)-Ta)/G;
 Energy1(3)=hc_pv2amb*(T1(2)-Ta)/G;
 Energy1(4)=Gpvabs1/G;
