function CellCal=solarcell(Ac,Gpv,G,cell)

Tpv0=30.79+273.15;
kB=1.38064852E-23; 
e=1.602176634E-19; 

if G~=0
    
    if strcmp(cell,'silicon')==1
        k1=0.03;
        b=1.25;
        n=2.98; 

    Ebg=1.11*e; 
    A1=1; 
    Tcoeff=0.0045;
    
    SR=xlsread('data', 'solar cells', 'I3:I2004'); 
 
    end
    
    
    if strcmp(cell,'CdTe')==1  
        k1=0.03;
        b=1.15; 
        n=0.98; 

    Ebg=1.44*e; 
    A1=1;
    Tcoeff=0.0025;
    
    SR=xlsread('data', 'solar cells', 'D3:D2004'); 

    end
    
 
    
    N=size(SR);
    J0=k1*Tpv0^(3/n)*exp(-Ebg/(b*kB*Tpv0));
    Jsc=0;
    lamadaG=xlsread('data', 'ambient conditions', 'B2:B2004');

        for i=1:N(1)
           Jsc=Jsc+Ac*Gpv(i)*SR(i)*(lamadaG(i+1)-lamadaG(i)); 
        end
    
    Voc=A1*kB*Tpv0/e*log(Jsc/J0+1);
    voc=Voc/(A1*kB*Tpv0/e);
    FF=(voc-log(voc+0.72))/(voc+1);
    Effm=Voc*Jsc*FF/(G*Ac);
    
     Gpv1=0;
    
        for j=1:N(1)
            Gpv1=Gpv1+Gpv(j)*(lamadaG(j+1)-lamadaG(j));
        end

  

    CellCal(1)=Effm;
    CellCal(2)=FF;
    CellCal(3)=Voc;
    CellCal(4)=Jsc;
    CellCal(5)=Gpv1;
    CellCal(6)=Tcoeff;

else 
    CellCal(1)=0;
    CellCal(2)=0;
    CellCal(3)=0;
    CellCal(4)=0;
    CellCal(5)=0;
    CellCal(6)=0;
end

end