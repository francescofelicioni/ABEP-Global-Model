function initial_conditions=Input_abep_cluster
global lambda gamma_rec R L ratio V S A R_foro n_cusps L_cusp exp_ratio
global k_0 atomic_mass T0 Kb Power_eff Conc in_flow_vec massflow_in 

lambda = 0.3; % 0 -> classical diffusion, 1 anomalous diffusion [Liebermann]
gamma_rec = 0.1; % recombination coefficient
Power_eff = 1; % Power efficiency

ratio = 1; % ratio between outlet and discharge surfaces [#]
V = R^2*pi*L; 
S = pi*R^2;
A = 2*pi*R*L+pi*R^2*(1+(1-ratio));
R_foro  = R*sqrt(ratio); %outlet radius [m]
n_cusps = 2; 
L_cusp  = 2*pi*R;
exp_ratio = 10; % expanding ratio throat-to-exit area


% Flow of intake
in_flow_vec   = (massflow_in*Conc(1:7))./(atomic_mass(1:7)*V);       % [m^-3 s^-1] Flow of O,N2,O2,N,CO2,NO2,CO

% Initial neutral density [m^-3] 
n_neutral_vec = (massflow_in.*Conc(1:7))./(S*ratio)./((k_0.*atomic_mass(1:7)/Kb).^0.5 .* 1/T0^0.5 .* (2./(k_0+1)).^((k_0+1)./2./(k_0-1)))./Kb./T0;   % O N2 O2 N CO2 NO2 CO
n_neutral     = sum(n_neutral_vec);

nie_0  = n_neutral/4000;       %initial chamber plasma density [m^-3]
np_vec = n_neutral_vec./4000;  %positive ion initial density [m^-3]
nn_vec = n_neutral_vec./12000; %negative ion initial density [m^-3]
Te_0 = 1;                   %initial chamber plasma temperature [eV]

initial_conditions = [n_neutral_vec(2), np_vec(2), n_neutral_vec(4), np_vec(4), n_neutral_vec(1), np_vec(1), nn_vec(1), n_neutral_vec(3), np_vec(3), nn_vec(3), ones(1,3)*nn_vec(3), n_neutral_vec(6), np_vec(6), nn_vec(6), ones(1,8)*nn_vec(7), ...
                      n_neutral_vec(5), n_neutral_vec(7), ones(1,4)*nn_vec(7),np_vec(5), nn_vec(7), np_vec(7), ones(1,21)*nn_vec(7) ,nie_0, Te_0]; 
% initial conditions = N2  N2+  N  N+  O  O+  O-  O2  O2+  O2-  3*O2-  NO2
% NO2+  NO2-  8*CO-  CO2  CO  4*CO-  CO2+  CO-  CO+  21*CO-  

end