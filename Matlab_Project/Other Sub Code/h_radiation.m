
function hr=h_radiation(Em1,Em2,T1,T2)
SBoltz=5.670367e-8; 

if Em2==0
    hr=Em1*SBoltz*(T1^2+T2^2)*(T1+T2); 
else
    hr=SBoltz*(T1^2+T2^2)*(T1+T2)./(1/Em1+1/Em2-1); 
end

end