clear; close all;

%0--k-1

Dse=4*10^-10; %cm2/s
Dpe=5*10^-9;
Dsm=4.4*10^-9;
Dpm=4.8*10^-9;
S0=6*10^-6;%mol/cm3
Vmax0=60*10^-6;%mol/(cm2*s)
L1=20*10^-7;%cm
L2=100*10^-7;
L=L1+L2;
T=12;%s
tse=T*L*L/Dse;
tpe=T*L*L/Dpe;
Vmaxs=Vmax0/S0*tse;
Vmaxp=Vmax0/S0*tpe;
Km0=33*10^-6;%mol/cm3
Km=Km0/S0;


X=L/L;

dt=1*10^-7;
h=0.01;
k=round(L1/L/h)+1;

N=X/h+1;
x=0:h:X;

M=T/dt;

K1=dt/h/h*Dsm/Dse;
K2=dt/h/h*Dpm/Dpe;

Se=zeros(N,1);
Pe=zeros(N,1);
Set=zeros(N,1);
Pet=zeros(N,1);
Se(:,1)=0;
Se(end,1)=S0/S0;
Pe(:,1)=0;


j=0;

while(j<=M);
    
    %1
    
    Pet(1)=0;
    
    %2--k-1
    for i=2:k-1;
    Set(i)=Se(i)+dt/h/h*(Se(i+1)-2*Se(i)+Se(i-1))-dt/T*Vmaxs*Se(i)/(Km+Se(i));
    Pet(i)=Pe(i)+dt/h/h*(Pe(i+1)-2*Pe(i)+Pe(i-1))+dt/T*Vmaxp*Se(i)/(Km+Se(i));
    end
    
    Set(1)=Set(2);
    
    %k-1--k+1
    for i=k;
    Set(i)=Se(i)+dt/h/h*(Dsm/Dse*Se(i+1)-(Dsm/Dse+1)*Se(i)+Se(i-1));
    Pet(i)=Pe(i)+dt/h/h*(Dpm/Dpe*Pe(i+1)-(Dpm/Dpe+1)*Pe(i)+Pe(i-1));
    end


   %k+1--N-1
    for i=k+1:N-1;
    Set(i)=Se(i)+Dsm/Dse*dt/h/h*(Se(i+1)-2*Se(i)+Se(i-1));
    Pet(i)=Pe(i)+Dpm/Dpe*dt/h/h*(Pe(i+1)-2*Pe(i)+Pe(i-1));
    end
    %N
    Set(N)=S0/S0;
    Pet(N)=0;  
    
   
    Se1=Se;
    Pe1=Pe;
    
    Se=Set;
    Pe=Pet;
       
    j=j+1;
    
    

    
end
    error1=max(abs(Se-Se1));
    error2=max(abs(Pe-Pe1));
tse=T*L*L/Dse;
tpe=T*L*L/Dpe;
x1=0:h*L*10^7:L*10^7;%nm
Se1=Se*S0*10^6;%mM
Pe1=Pe*S0*10^9;%uM




 subplot(1,2,1);plot(x1,Se1);ylabel('Se (mM)');
 subplot(1,2,2);plot(x1,Pe1);ylabel('Pe (uM)');
