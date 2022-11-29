clear, clc

%Measurements and temperature, sealer vessel
sealer_barrel = 0.249/2; 
sealer_outer = 0.5;                                                 %change
sealer_thick = 0.04; 
sealer_inner = sealer_outer - sealer_thick;
sealer_area = pi*(sealer_inner^2-sealer_barrel^2);
sealer_temp = 550;
sealer_temp_min = 325;
flow_pb = 50;                                                       %change

%Measurements and temperature, tube
L=1.7;                                                                %change
tube_inner = 0.01;                                                  %change
tube_thick = 0.001;                                                 %change
tube_outer = tube_inner+tube_thick;
tube_area = pi*tube_inner^2;
tube_inlet_temp = 140;
tube_max = 520;
kt = 20; %Heat conduction coefficient 
Nmax=sealer_outer*pi/tube_outer;
N=0.5*Nmax;                                                              %change

%Liquid parameters
density_pb = 11342; %Förenkling (temp- och tryckberoende)
density_lbe = 10440; %Förenkling (temp- och tryckberoende)
c_pb=128;
c_lbe=146.5; 
k = 35; %Thermal conductivity
velocity_pb = flow_pb/(density_pb*sealer_area);
flow_lbe = c_pb*flow_pb*(sealer_temp-sealer_temp_min)/(c_lbe*(tube_max-tube_inlet_temp));
velocity_lbe = flow_lbe/(density_lbe*N*tube_area);

%Heat transfer coefficients
alfa_pb = k/(density_pb*c_pb);
Pe_pb = tube_inner*2*velocity_pb/alfa_pb;
%Nu = 4.003 + 0.228*Pe^0.67;
Nu_pb = 6 + 0.006*Pe_pb;
h_pb = k*Nu_pb/(tube_inner*2); 

alfa_lbe = k/(density_lbe*c_lbe);
Pe_lbe = tube_inner*2*velocity_lbe/alfa_lbe;
%Nu = 4.003 + 0.228*Pe^0.67;
Nu_lbe = 6 + 0.006*Pe_lbe;
h_lbe = k*Nu_lbe/(tube_inner*2); 

%Overall heat transfer coefficient
thickness=sealer_thick+tube_thick;
U=1/(1/h_pb+1/h_lbe+thickness/kt);
%rin=sealer_inner;
%rout=sealer_outer+tube_thick; %Förenkling
%Reff = 1/(2*pi*rin*h) + log(rout/rin)/(2*pi*kt) + 1/(2*pi*rout*h); %Förenkling

count=100;
n=100;
dL=L/n;
T_lbe=linspace(tube_max,tube_inlet_temp, n);
T_pb=linspace(sealer_temp,sealer_temp_min, n);
dq = zeros(1, 100);

for ii=0:count
    for l=1:n-1
        MTD=T_pb(l)-T_lbe(l);
        dq(l)=N*MTD*U*2*tube_outer*dL;
        T_pb(l+1)=T_pb(l)-dq(l)/(flow_pb*c_pb);
    end
    for l=n:-1:2
        MTD=T_pb(round(l))-T_lbe(round(l));
        dq(l)=N*MTD*U*2*tube_outer*dL;
        T_lbe(l-1)=T_lbe(l)+dq(l)/(flow_lbe*c_lbe);
    end
end

plot(1:n, T_lbe)
hold on 
plot(1:n, T_pb)

Q0 = sum(dq)-dq(1)-dq(end)
Q1 = c_pb*flow_pb*(T_pb(1)-T_pb(end))
Q2 = c_lbe*flow_lbe*(T_lbe(1)-T_lbe(end))

rho_hot = 10981.7-1.1369*(520+273);
rho_cold = 10981.7-1.1369*(T_lbe(end)+273);
P = 10*20*(rho_cold-rho_hot)

