sealer_barrel = 0.249/2; %innerradie på sealer 
sealer_outer = 0.5; %ytterradie på sealer
sealer_thick = 0.04; %Tjocklek på på sealer
sealer_inner = sealer_outer - sealer_thick;
sealer_temp = 550;
flow = 175; %or 50

tube_height = 1;
tube_inner = 0.01; %Tubtjocklek, ändra
tube_thick = 0.001;
tube_outer = tube_inner+tube_thick;
tube_inlet_temp = 127;
tube_outlet_temp = 520;
kt = 20; %Heat conduction coefficient 

lead_density = 11342; %Förenkling (temp- och tryckberoende)
lead_bi_density = 10440; %Förenkling (temp- och tryckberoende)
c=128; %Heat capacity
k = 35; %Thermal conductivity

sealer_area = pi*(sealer_inner^2-sealer_barrel^2);
tube_area = pi*tube_inner^2;
velocity = flow/(lead_density*sealer_area);

alfa = k/(lead_density*c);
Pe = tube_inner*2*velocity/alfa;
%Nu = 4.003 + 0.228*Pe^0.67;
Nu = 6 + 0.006*Pe;
h = k*Nu/(tube_inner*2); %Heat transfer coefficent

rin=sealer_inner;
rout=sealer_outer+tube_thick; %Förenkling
Reff = 1/(2*pi*rin*h) + log(rout/rin)/(2*pi*kt) + 1/(2*pi*rout*h); %Förenkling

dL=0.01;
L=tube_height;

Q=0;
T_pb=550;
T_lbe=520;
for l=0:dL:L
    MTD = (T_pb-T_lbe);
    dq=MTD*dL/Reff;
    %dQ=dq*dL*tube_outer*2;
    dQ = dq;
    Q = Q + dQ;
    T_lbe = T_lbe - (520-127)*dL/L; %Linapp
    deltaT= dQ/(c*flow);
    T_pb = T_pb-deltaT;
end
T_lbe
T_pb

N=sealer_outer*pi/tube_outer
Q


