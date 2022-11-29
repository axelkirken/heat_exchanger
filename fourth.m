%sealer_outer, flow_pb, L, tube_inner, tube_thick, n
%% Temperaturprofiler (för olika flow)
[Q, T_pb, T_lbe] = energy(0.5, 175, 2, 0.01, 0.001, 1);
[Q2, T_pb2, T_lbe2] = energy(0.5, 50, 2, 0.01, 0.001, 1);
plot(1:n, T_lbe)
hold on 
plot(1:n, T_pb)
plot(1:n, T_pb2)
plot(1:n, T_lbe2)

%% Längd varierar
l=1:0.1:4;
q = zeros(1,length(l));

for i=1:length(l)
    [Q, T_pb, T_lbe] = energy(0.5, 175, l(i), 0.01, 0.002, 1);
    q(i) = Q;
end

plot(l,q)
xlabel("Tube length")
ylabel("Heat flow")


%% Tjocklek varierar
tube_thick=0.001:0.001:0.04;
q = zeros(1,length(tube_thick));

for i=1:length(tube_thick)
    [Q, T_pb, T_lbe] = energy(0.5, 175, 2, 0.02 , tube_thick(i), 1);
    q(i) = Q;
end

plot(tube_thick,q)
xlabel("Tube thickness")
ylabel("Energy")

%% Innerradie varierar
tube_inner = 0.005:0.001:0.05;
q = zeros(1,length(tube_inner));

for i=1:length(tube_inner)
    [Q, T_pb, T_lbe] = energy(0.5, 175, 2, tube_inner(i) , 0.002, 1);
    q(i) = Q;
end

plot(tube_inner,q)
xlabel("Tube inner radius")
ylabel("Energy")

%% Längd och innerradie varierar
l=1:0.1:4;
tube_inner = 0.005:0.001:0.05;
q = zeros(length(l),length(tube_inner));

for i=1:length(l)
    for j=1:length(tube_inner)
        [Q, T_pb, T_lbe] = energy(0.5, 50, l(i), tube_inner(j), 0.002, 1);
        q(i,j)=Q;
    end 
end

mesh(l,tube_inner,q')
xlabel("Tube length")
ylabel("Tube inner radius")
zlabel("Heat flow")
%% Innerradie och tubtjocklek
tube_thick = 0.001:0.0001:0.003;
tube_inner = 0.005:0.001:0.05;
q = zeros(length(tube_thick),length(tube_inner));

for i=1:length(tube_thick)
    for j=1:length(tube_inner)
        [Q, T_pb, T_lbe] = energy(0.5, 50, 2, tube_inner(j), tube_thick(i), 1);
        q(i,j)=Q;
    end 
end

mesh(tube_thick,tube_inner,q')
xlabel("Tube wall thickness")
ylabel("Tube inner radius")
zlabel("Energy")

%Dmitris comments:
%%Which equations, temperature profiles, how close we can get to the goal. 
%Add the script
%Poster: goals and tasks, motivation

%%
tube_thick = 0.001:0.0001:0.003;
tube_inner = 0.005:0.001:0.05;
q = zeros(length(tube_thick),length(tube_inner));

for i=1:length(tube_thick)
    for j=1:length(tube_inner)
        [Q, T_pb, T_lbe] = energy(0.5, 50, 2, tube_inner(j), tube_thick(i), 1);
        q(i,j)=Q;
    end 
end

mesh(tube_thick,tube_inner,q')
xlabel("Tube wall thickness")
ylabel("Tube inner radius")
zlabel("Energy")

%% Materialmängd vs. tjocklek
tube_thick=0.001:0.001:0.04;
m = zeros(1,length(tube_thick));

for i=1:length(tube_thick)
    [Q, T_pb, T_lbe, material] = energy(0.5, 175, 2, 0.02 , tube_thick(i), 1);
    m(i) = material;
end

plot(tube_thick,m)
xlabel("Thickness")
ylabel("Material")

