function hwtr=h_watergap(Tin,Tout,m,wgap,wc,Lc)


Tavg=(Tin+Tout)/2;
vw=exp(-3.7188+578.919/(-137.546+0.5*Tin+0.5*Tout))/1000;
kw=(9.28516*10^(-7)*Tavg^3 - 1.06167*10^(-2)*Tavg^2 + 7.76041*Tavg - 7.87144*10^2)/1000;
Prw=vw*4180/kw; 
Rew=(m/wgap/wc)*(2*wgap)/vw;

Lentr=0.05*Rew*Prw*(2*wgap);
h_entrance=(7.54+0.03*(2*wgap/Lentr)*Rew*Prw/(1+0.016*((2*wgap/Lentr)*Rew*Prw)^(2/3)))*kw/(2*wgap);
h_fulldev=7.54*kw/(2*wgap);

if Lentr<Lc
    hwtr=h_entrance*(Lentr/Lc)+h_fulldev*((Lc-Lentr)/Lc);
else
    hwtr=h_entrance;
end
    

end