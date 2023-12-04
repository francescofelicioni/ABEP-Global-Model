%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             GLOBAL MODEL                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main                                                                    %
%                                                                         %
% Written the 1/4/2020 by N. Souhair, University of Bologna               %
% Modified the 1/3/2022 by S. Dalle Fabbriche, University of Bologna      %
% Based on F. Bosi, PhD Thesis University of Padova                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vars_out,dens_out]= main_abep(Pw,flag_sim_mode,Gas_name,varargin)

global g q Kb B amu emass hP c_light eps0   
global k_0 atomic_mass Conc Rates R_const
global R L T0 Tpos_ion Tneg_ion 
global flag_electronegative flag_diffusion flag_sim_mode flag_atm Gas_name
global emass q Kb ratio V A S   
global k_energ Te Tneg_ion  
global flag_electronegative flag_atm gamma_rec 
global Pw q V n_e Te ntot Power_eff atomic_mass      
global Kb emass T0 V S ratio massflow_in  n_neutral ntot
global Te flag_electronegative 



%% INPUT DATA
addpath('./gas_data/');
load(strcat(Gas_name,'.mat'));

flag_electronegative = true;    % true=electronegative gas mixture diffusion model; false= electropositive gas mixture diffusion model;
flag_diffusion = 2;             % 0) LJ Mixture averaged ion diffusivity, 1) LJ self diffusion, 2) Langevin approximation
flag_atm = 'mars';              % 'mars','earth' --> choose planet atmosphere

% Geometry and init of parameters
L = 2;                       % [m] Length 
R = 0.5;                     % [m] Chamber radius
B = 1500e-4;                    % magnetic field [T]
%Pw= 400;
h_data = gas{5,2}; 

switch flag_sim_mode
    case 'inflow'
        mdot=varargin{1};
        [T0,massflow_in,Conc,~]=data_from_height_or_inflow(flag_sim_mode,mdot);
    case 'height'
        h=varargin{1};
        [T0,massflow_in,Conc,rho_air]=data_from_height_or_inflow(flag_sim_mode,h,h_data,flag_atm);
end
 
Tpos_ion = T0;
Tneg_ion = 1.5*T0;    % higher temperature for negative ions due to no wall recombination 

% Constants
Kb = 1.380649e-23;        % boltzmann constant [J/K]
hP = 4.135667696e-15; % Plank constant [eV*s]
c_light = 299792458;  % speed of light [m/s]
q  = 1.60218e-19;         % electron charge [C]
switch flag_atm
    case 'earth'
        g  = 9.81;    % gravitational acceleration [m/s^2]
    case 'mars'
        g = 3.71;     % gravitational acceleration [m/s^2]
end
emass = 9.1094e-31;     % electron mass [kg]
amu   = 1.66054e-27;     % proton mass [kg]
eps0  = 8.8542e-12;   % [C^2/N/m^2] vacuum permittivity
m_N   = 14*amu;                 % [kg] Atomic mass N 
m_O   = 16*amu;                 % [kg] Atomic mass O
m_C   = 12*amu;                 % [kg] Atomic mass C

% Reaction Rate Coefficients
Te_interp       = linspace(0.01,50,500);
[Rates,R_const] = Reaction_rates(T0,Te_interp);

% Gas data
k_0  = [5/3 7/5 7/5 5/3 4/3 4/3 7/5];     % specific heat constant for O,N2,O2,N,CO2,NO2,CO
atomic_mass  = amu*table2array(gas{3,2}); % atomic mass [kg]

%% Input
initial_conditions = Input_abep_cluster;

%% SOLVING
time       = [0 5];
options    = odeset('RelTol',1e-6,'AbsTol',1e-7,'Maxorder',2,'NonNegative',[]);
ode_input  = {time, initial_conditions,options};
tic; [t,y] = ode15s(@Plasma_dyn, ode_input{:}); toc

%% RESULTS (@ t=t_end)

sum_ions  = y(:,2)+y(:,4)+y(:,6)+y(:,9)+y(:,12)+y(:,15)+y(:,18)+y(:,31)+y(:,32)+y(:,33)+y(:,34)+y(:,35)+y(:,36)+y(:,37)+y(:,38)+y(:,39);
sum_neutr = y(:,1)+y(:,3)+y(:,5)+y(:,8)+y(:,11)+y(:,14)+y(:,17)+y(:,25)+y(:,26)+y(:,27)+y(:,28)+y(:,29)+y(:,30);

switch flag_sim_mode
    case 'height'
        n_air_vec = rho_air ./ [m_O  2*m_N  2*m_O  m_N  m_C+2*m_O  m_N+2*m_O  m_C+m_O];
        n_air = sum(n_air_vec);
               
        fprintf('nintake = %.2e m^-3\n', n_air)
    case 'inflow'
        mdot_vec = mdot*Conc;
        n_intake_vec = mdot_vec./[m_O  2*m_N  2*m_O  m_N  m_C+2*m_O  m_N+2*m_O  m_C+m_O];
        n_intake = sum(n_intake_vec);

        fprintf('intake mass flow = %.2e \nintake number flow = %.2e \n',mdot,n_intake)
end    
ionization_ratio = sum_ions(end)./sum(y(end,1:54));         % [#] Ionization ratio (ions/(all but electrons)) 
pressure = 760/101325*1e3.*sum_neutr.*T0*Kb;           % [mTorr] Neutrals pressure 
Pressure = pressure(end);
%NEUTRALS OUTPUT
n_N2 = y(:,1); n_N = y(:,3); n_O = y(:,5); n_O2 = y(:,8); n_NO = y(:,11); n_NO2 = y(:,14); n_N2O = y(:,17); 
n_CO2 = y(:,25); n_CO = y(:,26); n_C2O = y(:,27); n_C = y(:,28); n_C2 = y(:,29); n_O3 = y(:,30);
%IONS OUTPUT
n_N2i = y(:,2); n_Ni = y(:,4); n_Oi = y(:,6); n_O2i = y(:,9); n_NOi = y(:,12); n_NO2i = y(:,15); n_N2Oi = y(:,18); 
n_CO2i = y(:,31); n_CO4i = y(:,32); n_COi = y(:,33); n_C2O2i = y(:,34); n_C2O3i = y(:,35); n_C2O4i = y(:,36); n_Ci = y(:,37);
n_C2i = y(:,38); n_O4i = y(:,39); 
%ELECTRONS OUTPUT
ne = y(:,end-1);
Te = y(:,end);
[Thrust, Isp, eta] = propulsion_model(y(end,1:end-1), Te(end), T0); % PROPULSION 
T_over_P = Thrust*1e3/(Pw*1e-3);
vars_out = [Thrust*1e3, Isp, eta*1e3, ne(end), Te(end),Pressure,ionization_ratio(end)]; % main function output
dens_out = [t ne(:,end) Te];
fprintf('Te = %.2f eV \nne = %.2e m^-3 \nT = %.0f K \nionization_ratio = %.2f \npressure = %.8f mTorr\n', Te(end), ne(end), T0, ionization_ratio, Pressure(end));
fprintf('Thrust = %.1f mN \nIsp = %.0f s \nefficiency = %.2f\n', 1e3*Thrust, Isp, eta*100);
fprintf('DENSITIES OF NEUTRALS\nN2: %.2e \nN: %.2e \nO: %.2e \nO2: %.2e \nNO: %.2e \nNO2: %.2e \nN2O: %.2e \nCO2: %.2e \nCO: %.2e \nC2O: %.2e \nC: %.2e \nC2: %.2e \nO3: %.2e \n', n_N2(end), ...
    n_N(end), n_O(end), n_O2(end), n_NO(end), n_NO2(end), n_N2O(end), n_CO2(end), n_CO(end), n_C2O(end), n_C(end), n_C2(end), n_O3(end));
fprintf(['DENSITIES OF IONS\nN2i: %.2e \nNi: %.2e \nOi: %.2e \nO2i: %.2e \nNOi: %.2e \nNO2i: %.2e \nN2Oi: %.2e \nCO2i: %.2e \nCO4i: %.2e \nCOi: %.2e \nC2O2i: %.2e \nC2O3i: %.2e \nC2O4i: %.2e \n' ...
    'nCi: %.2e \nnC2i: %.2e \nnO4i: %.2e \n'], n_N2i(end), n_Ni(end), n_Oi(end), n_O2i(end), n_NOi(end), n_NO2i(end), n_N2Oi(end), n_CO2i(end), n_CO4i(end), ...
    n_COi(end), n_C2O2i(end), n_C2O3i(end), n_C2O4i(end), n_Ci(end), n_C2i(end), n_O4i(end));
clear global
end











%% DATA FROM HEIGHT OR INFLOW

function [T0,massflow_in,Conc,rho_air]=data_from_height_or_inflow(flag_sim_mode,varargin) %value=altitude  name=flag_sim_mode (height or inflow)  h_data=stored in gas (gas{5,2})  flag_atm=earth or mars

%this function evaluates the T0 and the inflow or the orbit height depending on the input
global massflow_in
amu   = 1.66054e-27;               % [kg] Proton mass
m_N   = 14*amu;                 % [kg] Atomic mass N 
m_O   = 16*amu;                 % [kg] Atomic mass O
m_C   = 12*amu;                 % [kg] Atomic mass C
Sez_intake = 1;
R_intake   = sqrt(Sez_intake/pi); % suppose GOCE like spacecraft frontal area and intake covering it
Sez_intake = pi*R_intake^2;
intake_eff = 0.43; % Romano et al. 2021
switch flag_sim_mode
    case 'height'
        h = varargin{1};
        h_data = varargin{2};
        flag_atm = varargin{3};
        height=find(h_data(:,1)==h);
        switch flag_atm
            case 'earth'
                n_air  = h_data(height,2:8)*1e6;
                v_orb  = sqrt(3.986e14/((h+6371)*1e3));       % Orbital speed (Earth) [m/s]
                for i=2:8
                    Conc(i-1)=h_data(height,i)/sum(h_data(height,2:8)); % Concentrations
                end
                T0     = h_data(height,9);
                rho_air= sum(n_air.*[m_O m_N*2 m_O*2 m_N (m_C+m_O*2) (m_N+m_O*2) (m_C+m_O)],2);
                inflow = v_orb*rho_air*Sez_intake*intake_eff;
                massflow_in = inflow;
            case 'mars' 
                n_air  = h_data(height,2:8);
                v_orb  = sqrt(4.26213e13/((h+3389.5)*1e3));   % Orbital speed (Mars)  [m/s]
                for i=2:8
                    Conc(i-1)=h_data(height,i)/sum(h_data(height,2:8));
                end
                T0     = h_data(height,9);
                rho_air= sum(n_air.*[m_O m_N*2 m_O*2 m_N (m_C+m_O*2) (m_N+m_O*2) (m_C+m_O)],2);
                inflow = v_orb*rho_air*Sez_intake*intake_eff;
                massflow_in = inflow;
        end
    case 'inflow'
        mdot=varargin{1};
        Conc = [0 0 0 0 1 0 0];             %Concentrations of the species at the intake (O N2 O2 N CO2 NO2 CO)%
        inflow = mdot;
        massflow_in = inflow;
        rho_air = 0;
        T0 = 400;
    otherwise
        fprintf('You have to choose between height or inflow');
end
end

%% PLASMA_DYN
function vars = Plasma_dyn(t,y)
global Pw q V n_e Te ntot Power_eff atomic_mass Kb      
global n_neutral n_ion total_mass_vector n_vec m_vector_neutral_charged m_vector_neutral
global n_N2 n_N2i n_N n_Ni n_O n_Oi n_On n_O2 n_O2i n_O2n n_NO n_NOi n_NOn n_NO2 n_NO2i n_NO2n n_N2O n_N2Oi n_N2On n_O2v1 n_O2v2 n_O2v3 n_O2e1 n_O2e2
global n_CO2 n_CO n_C2O n_C n_C2 n_O3 n_CO2i n_CO4i n_COi n_C2O2i n_C2O3i n_C2O4i n_Ci n_C2i n_O4i n_CO3n n_CO4n n_O3n n_O4n n_CO2v1 n_CO2v2 n_CO2v3 n_CO2v4 n_CO2e1 n_CO2e2 n_COv1 n_COe1 n_COe2 n_COe3 n_COe4
%% UPDATE FIELDS
n_N2    =  y(1);          % N2  density     [particle/m^3]
n_N2i   =  y(2);          % N2+ density
n_N     =  y(3);          % N density
n_Ni    =  y(4);          % N+ density
n_O     =  y(5);          % O density
n_Oi    =  y(6);          % O+ density
n_On    =  y(7);          % O- density
n_O2    =  y(8);          % O2 density
n_O2i   =  y(9);          % O2i density
n_O2n   =  y(10);         % O2n density
n_NO    =  y(11);         % NO density
n_NOi   =  y(12);         % NOi density
n_NOn   =  y(13);         % NOn density
n_NO2   =  y(14);         % NO2 density
n_NO2i  =  y(15);         % NO2i density
n_NO2n  =  y(16);         % NO2n density
n_N2O   =  y(17);         % N2O density
n_N2Oi  =  y(18);         % N2Oi density
n_N2On  =  y(19);         % N2On density
n_O2v1  =  y(20);         % O2v1 density
n_O2v2  =  y(21);         % O2v2 density
n_O2v3  =  y(22);         % O2v3 density
n_O2e1  =  y(23);         % O2e1 density
n_O2e2  =  y(24);         % O2e2 density
n_CO2   =  y(25);         % CO2 density
n_CO    =  y(26);         % CO density
n_C2O   =  y(27);         % C20 density
n_C     =  y(28);         % C density
n_C2    =  y(29);         % C2 density
n_O3    =  y(30);         % O3 density
n_CO2i  =  y(31);         % CO2+ density
n_CO4i  =  y(32);         % CO4+ density
n_COi   =  y(33);         % CO+ density
n_C2O2i =  y(34);         % C2O2+ density
n_C2O3i =  y(35);         % C2O3+ density
n_C2O4i =  y(36);         % C2O4+ density
n_Ci    =  y(37);         % C+ density
n_C2i   =  y(38);         % C2+ density
n_O4i   =  y(39);         % O4+ density
n_CO3n  =  y(40);         % CO3- density
n_CO4n  =  y(41);         % CO4- density
n_O3n   =  y(42);         % O3- density
n_O4n   =  y(43);         % O4- density
n_CO2v1 =  y(44);         % CO2v1 density
n_CO2v2 =  y(45);         % CO2v2 density
n_CO2v3 =  y(46);         % CO2v3 density
n_CO2v4 =  y(47);         % CO2v4 density
n_CO2e1 =  y(48);         % CO2e1 density
n_CO2e2 =  y(49);         % CO2e2 density
n_COv1  =  y(50);         % COv1 density
n_COe1  =  y(51);         % COe1 density
n_COe2  =  y(52);         % COe2 density
n_COe3  =  y(53);         % COe3 density
n_COe4  =  y(54);         % COe4 density
n_e     =  y(55);         % Electron density
Te      =  y(56);         % Electron temperature

P_absorbed = (Power_eff*Pw/q/V)*(1-exp(-t/.0001)); % This is a temporal profile of absorbed energy (facilitate the convergence)

% Useful density and mass vectors
n_neutral = [n_N2,n_N,n_O,n_O2,n_NO,n_NO2,n_N2O,n_CO2,n_CO,n_C2O,n_C,n_C2,n_O3]'; 
n_ion     = [n_N2i,n_Ni,n_Oi,n_O2i,n_NOi,n_NO2i,n_N2Oi,n_CO2i,n_COi,n_CO4i,n_C2O2i,n_C2O3i,n_C2O4i,n_Ci,n_C2i,n_O4i]; 
ntot      = sum(y(1:54));
n_vec     = y(1:54);
total_mass_vector        = [atomic_mass(2)*ones(1,2) atomic_mass(4)*ones(1,2) atomic_mass(1)*ones(1,3) atomic_mass(3)*ones(1,3) (atomic_mass(1)+atomic_mass(4))*ones(1,3) atomic_mass(6)*ones(1,3) (atomic_mass(1)+atomic_mass(2))*ones(1,3) atomic_mass(3)*ones(1,5) ...
                            atomic_mass(5) atomic_mass(7) 2*atomic_mass(7)-atomic_mass(1) atomic_mass(7)-atomic_mass(1) 2*atomic_mass(7)-atomic_mass(3) atomic_mass(3)+atomic_mass(1) atomic_mass(5) atomic_mass(5)+atomic_mass(3) atomic_mass(7) 2*atomic_mass(7) ...
                            2*atomic_mass(7)+atomic_mass(1) 2*atomic_mass(5) atomic_mass(7)-atomic_mass(1) 2*atomic_mass(7)-atomic_mass(3) 2*atomic_mass(3) atomic_mass(7)+atomic_mass(3) atomic_mass(5)+atomic_mass(3) 3*atomic_mass(1) 2*atomic_mass(3) atomic_mass(5)*ones(1,6) atomic_mass(7)*ones(1,5)];
m_vector_neutral_charged = [atomic_mass(2) atomic_mass(4) atomic_mass(1) atomic_mass(3) atomic_mass(1)+atomic_mass(4) atomic_mass(6) atomic_mass(1)+atomic_mass(2) atomic_mass(5) atomic_mass(7) atomic_mass(5)+atomic_mass(3) 2*atomic_mass(7) ...
                            2*atomic_mass(7)+atomic_mass(1) 2*atomic_mass(5) atomic_mass(7)-atomic_mass(1) 2*atomic_mass(7)-atomic_mass(3) 2*atomic_mass(3)];
m_vector_neutral         = [atomic_mass(2) atomic_mass(4) atomic_mass(1) atomic_mass(3) atomic_mass(1)+atomic_mass(4) atomic_mass(6) atomic_mass(1)+atomic_mass(2) atomic_mass(5) atomic_mass(7) 2*atomic_mass(7)-atomic_mass(1) atomic_mass(7)-atomic_mass(1) 2*atomic_mass(7)-atomic_mass(3) 3*atomic_mass(1)];


%----------------------------------------------------------------------------------------------------------------------------------------------------
% UPDATE MODELS
gas_diffusion;

[dn_bulk_chem_reac_plasma, Power_chem_reac_loss_plasma] = plasma_decomposition_abep(Te,y(1:end-1)); % plasma chemical reaction in the bulk plasma
[dn_wall, P_lost_wall]  =  update_walls; % wall reaction
[dn_inlet, dn_outlet, P_inlet, P_exhauste] = update_flows_in_out; % inlet/outlet fluxes

%----------------------------------------------------------------------------------------------------------------------------------------------------
% POPULATION BALANCE
dn = dn_bulk_chem_reac_plasma' + dn_inlet - dn_outlet + dn_wall';

%----------------------------------------------------------------------------------------------------------------------------------------------------
% POWER BALANCE
Ptot = P_absorbed + P_inlet - Power_chem_reac_loss_plasma - P_lost_wall - P_exhauste;
dTe = 2/(3*n_e)*(Ptot-3/2*Te*dn(end)); % d/dt(3/2*q*ne*Te) = P => 3/2*q*ne*d(Te)/dt+3/2*q*Te*d(ne)/dt
% dexhausted=-dn_outlet+dn_wall';  %Add this variable to check particle balances
% dnin=dn_inlet;                   %Add this variable to check particle balances
vars=[dn dTe]'; 
end

%% UPDATE WALLS
function [dn_wall, P_lost_wall]  =  update_walls
global emass q Kb ratio V A S   
global k_energ Te Tneg_ion m_vector_neutral  
global flag_electronegative flag_atm gamma_rec 
global n_ion n_neutral total_mass_vector n_vec n_neg n_e ntot 
global ub_vec Aeff_vec v_th_vec Beta_vec hL_vec T0

% WALL RECOMBINATION: loss by wall recombination 
wall_ion  = ub_vec.*n_ion.*Aeff_vec./V;
v_th_vecn = sqrt(8*Kb*T0./(pi*m_vector_neutral));   % [m/s] Neutrals Thermal velocity
wall_N = 0.25*n_neutral(2).*v_th_vecn(2)*gamma_rec*A/V;  
wall_O = 0.25*n_neutral(3).*v_th_vecn(3)*gamma_rec*A/V;  

% Wall reactions
%(N2+ + e -> N2; N+ + e -> N; O+ + e -> O; O2+ + e -> O2; NO+ + e -> NO;  NO2+ + e -> NO2;  N2O+ + e -> N2O; N -> 1/2 N2; O -> 1/2 O2) 
%(CO2+ + e -> CO2; CO+ + e -> CO; C+ + e- > C; C2+ + e -> C2;
wall_rec = [wall_ion wall_N wall_O]';

switch flag_atm
    case 'earth'
        M_wall_rec = [ 1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0.5 0;      %N2
                      -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %N2+
                       0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0;      %N
                       0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %N+
                       0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1;      %O
                       0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O+
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %On
                       0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 0.5;     %O2
                       0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2n
                       0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0;      %NO
                       0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0;      %NOi
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %NOn
                       0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0;      %NO2
                       0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0;      %NO2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %NO2n
                       0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0;      %N2O
                       0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0;      %N2Oi
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %N2On
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2v1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2v2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2v3
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2e1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2e2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2O
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O3
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO4i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COi
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2O2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2O3i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2O4i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %Ci
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O4i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO3n
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO4n
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O3n
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O4n
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2v1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2v2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2v3
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2v4
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2e1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2e2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COv1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COe1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COe2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COe3
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COe4
                      -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  0  0  0  0  0  0];     %e

    case 'mars'
        %(N2+ + e -> N2; N+ + e -> N; O+ + e -> O; O2+ + e -> O2; NO+ + e -> NO;  NO2+ + e -> NO2;  N2O+ + e -> N2O; N -> 1/2 N2; O -> 1/2 O2) 
        %(CO2+ + e -> CO2; CO+ + e -> CO; C+ + e- > C; C2+ + e -> C2;
        M_wall_rec = [ 1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0.5 0;      %N2
                      -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %N2+
                       0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0;      %N
                       0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %N+           
                       0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1;      %O
                       0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O+
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %On
                       0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 0.5;     %O2
                       0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2n
                       0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0;      %NO
                       0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0;      %NOi
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %NOn
                       0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0;      %NO2
                       0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0;      %NO2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %NO2n
                       0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0;      %N2O
                       0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0;      %N2Oi
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %N2On
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2v1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2v2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2v3
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2e1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O2e2
                       0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0;      %CO2
                       0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0;      %CO
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2O
                       0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0;      %C
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0;      %C2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O3
                       0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0;      %CO2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO4i
                       0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0;      %COi
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2O2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2O3i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %C2O4i
                       0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;      %Ci
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0;      %C2i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O4i
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO3n
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO4n
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O3n
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %O4n
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2v1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2v2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2v3
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2v4
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2e1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %CO2e2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COv1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COe1
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COe2
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COe3
                       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;      %COe4
                      -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0 -1 -1  0  0  0];     %e  
end

% WALL FLUXES
dn_wall = M_wall_rec*wall_rec;

if flag_electronegative==0       %electropositive    
    % Energy lost per electron hitting the walls (with reflection and
    % secondary emission)  from Zhou's paper
    Mr = total_mass_vector*n_vec/ntot; 
    E_S = 50;                                  % [eV] (for boron nitride walls)           
    E_R = 20;                                  % [eV] (for boron nitride walls)                   
    gamma_0 = 0.4;                             % [#] (for boron nitride walls)   
    Ts = 2;                                    % [eV] Temperature of electrons emitted from the wall 
    gamma_s = (2*Te)/E_S;                      % [#] Secondary emission coefficient 
    gamma_r = gamma_0*((E_R^2)/(Te +E_R)^2);   % [#] Reflection coefficient     
    phi_wall = Te*log(sqrt((Mr)/(2*pi*emass))*(1-gamma_r)*(1-gamma_s));    % Potentials at the wall
    k_energ = 2*Te/(1-gamma_s) - gamma_s*((2*Ts)/(1-gamma_s)) + phi_wall;  % Energy lost per electron hitting the walls [eV]
    
elseif flag_electronegative==1   %electronegative   
    wall_ion = ub_vec.*n_ion.*(Aeff_vec+S*ratio*Beta_vec.*hL_vec)./V;      % electronegative model includes exhaust flow in wall losses power evaluation
    alpha=n_neg/n_e;
    gamma_minus=(Te*q/Kb)/Tneg_ion;
    f=@(Vpsym) Vpsym - (Te*(alpha*exp(-(Vpsym*(gamma_minus - 1))/Te) + 1))/(2*(alpha*gamma_minus*exp(-(Vpsym*(gamma_minus - 1))/Te) + 1));
    options = optimoptions('fsolve','Display','off');
    Vp=fsolve(f,2,options);
    alpha_s=alpha*exp(Vp/Te*(1-gamma_minus));    
    ub_tilde_vec  = ub_vec*((1+alpha_s)/(1+alpha_s*gamma_minus))^0.5;
    ub_tilde_mean = (n_ion*ub_tilde_vec')/sum(n_ion);
    vth_e = sqrt(8*Te*q/pi/emass);
    v_th_mean = (v_th_vec(3:end)*n_ion(3:end)')/sum(n_ion(3:end));
    Vs = abs(Te*log(4*ub_tilde_mean/vth_e*(1+alpha_s)/(1+alpha_s*(v_th_mean/vth_e)^2)));    
    k_energ=(2*Te+Vp+Vs);
end
P_lost_wall = k_energ*(sum(wall_ion)); % Power loss at wall and exhaust

end

%% UPDATE FLOWS IN OUT
function [dn_inlet, dn_outlet, P_inlet, P_exhauste] = update_flows_in_out
global Kb emass T0 V S ratio massflow_in 
global n_ion n_neutral total_mass_vector m_vector_neutral n_vec ntot
global Te flag_electronegative 
global Beta_vec hL_vec ub_vec in_flow_vec

%----------------------------------------------------------------------------------------------------------------------------------------------------
% OUTLET (chocked nozzle + magnetic nozzle)
% Electron loss 
ns   =  n_ion.*hL_vec;
% Ion loss because of magnetic mirror to deliver momentum
exhaust_ion = Beta_vec.*ns.*ub_vec*ratio*S/V;
exhaust_e = sum(exhaust_ion);
k_vec     = [7/5,5/3,5/3,7/5,7/5,4/3,4/3,4/3,7/5,4/3,5/3,7/5,4/3];

W    = (S*ratio) .* (n_neutral*Kb*T0).* (k_vec.*m_vector_neutral'./Kb).^0.5 .* 1/T0^0.5 .* (2./(k_vec+1)).^((k_vec+1)./2./(k_vec-1));  % [kg/s] Neutral gas mass flow (isoentropic nozzle mass flow)
       
exhaust_all = W./(m_vector_neutral'.*V);
dn_outlet =[exhaust_all(1), exhaust_ion(1), exhaust_all(2), exhaust_ion(2), exhaust_all(3), exhaust_ion(3), 0, exhaust_all(4), exhaust_ion(4), 0,...
            exhaust_all(5), exhaust_ion(5), 0, exhaust_all(6), exhaust_ion(6), 0, exhaust_all(7), exhaust_ion(7), zeros(1,6), exhaust_all(8), exhaust_all(9), ...
            exhaust_all(10), exhaust_all(11), exhaust_all(12), exhaust_all(13), exhaust_ion(8), exhaust_ion(10), exhaust_ion(9), exhaust_ion(11), ...
            exhaust_ion(12), exhaust_ion(13), exhaust_ion(14), exhaust_ion(15), exhaust_ion(16), zeros(1,15), exhaust_e];
        
if flag_electronegative==0 %electropositive
    % Mean energy lost per electron leaving the discharge 
    % [10.2.4 Liebermann]
    % Reduced mass [kg]
    Mr = total_mass_vector*n_vec/ntot;       
    % Energy lost at exhaust [eV]
    k_energ_ex = (2 - log(sqrt(2*pi*emass/Mr)))*Te;     
    P_exhauste = exhaust_e*k_energ_ex; % magnetic nozzle loss [W/m^3]  % Power loss at exhaust (magnetic nozzle)
elseif flag_electronegative==1 %electronegative   
    P_exhauste = 0;  %already counted in wall losses    
end

%----------------------------------------------------------------------------------------------------------------------------------------------------
% INLET
u_i = in_flow_vec./ntot./S; % inlet velocity NOTE: this is approximated u_i=u_e
dn_inlet = [in_flow_vec(2),0,in_flow_vec(4),0,in_flow_vec(1),0,0,in_flow_vec(3),zeros(1,5),in_flow_vec(6),zeros(1,10),in_flow_vec(5),in_flow_vec(7),zeros(1,29)];
P_inlet = 0.5*massflow_in*sum(u_i.^2); % negligible inlet flow power

end

%% PROPULSION MODEL
function [Thrust, Isp, eta] = propulsion_model(n,Te,T)
% GIVEN THE EXPANDING RATIO WE CALCULATE THRUST, EXIT T AND PRESSURE...

global Kb q g Beta_vec B hL_vec m_vector_neutral_charged m_vector_neutral
global R ratio Pw exp_ratio T0
options = optimset('Display','off');

% GEOMETRICAL PARAMETERS
A_th = ratio*R^2*pi; % throat area
R_th = sqrt(A_th/pi);
A_exit = A_th*exp_ratio;   

% Densities and masses
n_neutral = [n(1) n(3) n(5) n(8) n(11) n(14) n(17) n(25) n(26) n(27) n(28) n(29) n(30)];
n_charged = [n(2) n(4) n(6) n(9) n(12) n(15) n(18) n(31) n(32) n(33) n(34) n(35) n(36) n(37) n(38) n(39)];
n_tot     = sum(n_neutral); 
k_vec     = [7/5,5/3,5/3,7/5,7/5,4/3,4/3,4/3,7/5,4/3,5/3,7/5,4/3];

%----------------------------------------------------------------------------------------------------------------------------------------------------
% CHARGED PART CALCULATION
ub_vec = sqrt(q.*Te./m_vector_neutral_charged);             % [m/s] Bohm velocity  
F0_vec = q*Beta_vec.*(hL_vec*2.*n_charged)*Te*pi*R_th^2;    % [N] Upstream plasma thrust 
v_t_vec = sqrt(q*T./m_vector_neutral_charged);      % CONTROLLA Q O KB
omega_ci_vec = q*B./m_vector_neutral_charged;
r_ci_vec = v_t_vec./omega_ci_vec;

for i=1:length(n_charged)
    if r_ci_vec(i) >= R_th
        M_det_vec(i) = 1;         % Already detached, demagnetized
    elseif r_ci_vec(i) < R_th
        % Detachment Mach number
        M_det_vec(i) = fsolve(@(M_det_vec)(log(q*B^2*R_th^2/m_vector_neutral_charged(i)/T)+log(M_det_vec)-1/2*(M_det_vec^2-1)),1.1,options);  
    end
end

% [eq 39 "Helicon plasma thruster discharge model"]
S_ion_vec = (M_det_vec.^2+1)./2./M_det_vec.*F0_vec;   % Ion Thrust

% Ion flux (with hL=hL)
% (Beta*n = radially averaged plasma density at a given axial location)
ideal_flux_vec = hL_vec.*Beta_vec.*n_charged.*m_vector_neutral_charged.*ub_vec.*A_th;

%----------------------------------------------------------------------------------------------------------------------------------------------------
% TOTAL CHARGED PART
S_ion      = sum(S_ion_vec);
ideal_flux = sum(ideal_flux_vec);

%----------------------------------------------------------------------------------------------------------------------------------------------------
% NEUTRAL PART CALCULATION
p0 = Kb*T0*sum(n_neutral);    % Chamber total pressure
X  = n_neutral./n_tot; % molar fractions
p0_vec = p0*X;
Rgas = Kb./m_vector_neutral; % gas constant of the mixture
% Neutral gas mass flow (isoentropic nozzle mass flow) [kg/s]
W_vec = A_th*p0_vec.*(k_vec./Rgas./T.*(2./(k_vec+1)).^((k_vec+1)./(k_vec-1))).^0.5;

% Exit Mach number
for j=1:length(n_neutral)
    M_vec(j) = fzero(@(M)(A_exit/A_th-(1/M)*((1+(k_vec(j)-1)/2*M^2)/((k_vec(j)+1)/2))^((k_vec(j)+1)/2/(k_vec(j)-1))),2,options);
end

p_exit_vec = p0_vec.*((1+(k_vec-1)./2.*M_vec.^2)).^(k_vec./(1-k_vec));       % Exit pressure
T_exit_vec = T./(1+(k_vec-1)./2.*M_vec.^2);                                  % Exit temperature
c_ex_vec   = (k_vec.*Rgas.*T_exit_vec).^0.5;                                 % Sonic velocity at the exit section
u_vec      = c_ex_vec.*M_vec;                                                % Exit velocity
S_neutral_vec = W_vec.*u_vec + p_exit_vec*A_exit;                            % thrust neutrals

%----------------------------------------------------------------------------------------------------------------------------------------------------
% TOTAL NEUTRAL PART
S_neutral = sum(S_neutral_vec);
W = sum(W_vec);

%----------------------------------------------------------------------------------------------------------------------------------------------------
% PROPULSION FIGURES OF MERIT
flow_rate = W + ideal_flux;              % total flow rate
Isp = (S_neutral + S_ion)/flow_rate/g;   % Specific Impulse
Thrust = (S_neutral + S_ion);            % Thrust = S_ion + S_neutral [N];
eta = Isp*Thrust/2*g/Pw;                 % efficiency
end

%% GAS DIFFUSION
function gas_diffusion   %%CHECK IF IT IS A MIXTURE

%ASSUMPTION from Marmuse: negative ions flux to the walls can be neglected due to low thermal speed and higher mass with respect to electrons

global Kb T0 eps0 emass q B lambda Tpos_ion Tneg_ion Te
global R_foro n_cusps L_cusp Rates L R  
global n_e ntot m_vector_neutral_charged total_mass_vector n_vec n_neg
global flag_electronegative flag_diffusion 
global Beta_vec hL_vec Aeff_vec ub_vec v_th_vec v_th_e
global n_N2 n_N2i n_N n_Ni n_O n_Oi n_On n_O2 n_O2i n_O2n n_NO n_NOi n_NOn n_NO2 n_NO2i n_NO2n n_N2O n_N2Oi n_N2On
global n_CO2 n_CO n_C2O n_C n_C2 n_O3 n_CO2i n_CO4i n_COi n_C2O2i n_C2O3i n_C2O4i n_Ci n_C2i n_O4i n_CO4n n_O4n n_CO3n n_O3n


% Parameters for electrons
v_th_e = sqrt(8*Kb*Te/(pi*emass));       % [m/s] Thermal velocity
omega_c = q*B/emass;                    % [rad/s] Cyclotron frequency
r_le = v_th_e/omega_c;                  % [m] Larmor radius
n_M = sum([n_N2,n_N,n_O,n_O2,n_NO,n_NO2,n_N2O,n_CO,n_CO2,n_O3,n_C2O,n_C,n_C2]);  % total neutral number density (third body in a process)
n_vector_ne=[n_Ni*n_M, n_e*n_Ni, n_N, n_N2i, n_N2i*n_M, n_e*n_N2i, n_N2, n_Oi*n_M, n_e*n_Oi,n_O,...
                     n_O*n_O2, n_O*n_O2, n_O2i, n_e*n_O2i, n_O2i*n_M, n_O2, n_O2, n_O2, n_O2, n_O2*n_O2,...
                     n_O2*n_N2, n_NOi, n_NOi*n_M, n_e*n_NOi, n_NO*n_M, n_NO2i, n_NO2*n_M, n_NO2, n_N2Oi, n_N2O,...
                     n_N2, n_N2, n_N2, n_N2i, n_N, n_N, n_Ni, n_O, n_O, n_Oi, n_O2, n_O2*ones(1,5), n_CO2*ones(1,14),...
                     n_CO*ones(1,11), n_C,n_C,n_C2,n_C2,n_C2,n_O3,n_O3,n_O3,n_CO2i,n_CO2i,n_CO4i,...
                     n_COi,n_C2O2i,n_C2O3i,n_C2O4i,n_C2i,n_O2*n_M,n_O3*n_M,n_O*n_M,n_O2i*n_M,n_O4i];


Rates_vec= NaN(92,1);

for i=1:46
    Rates_vec(i)= Rates{i,1}(Te);
end
for j=1:46
    Rates_vec(46+j)= Rates{j+49,1}(Te);
end

ve = n_vector_ne*Rates_vec;

v_th_vec = sqrt(8*Kb*T0./(pi*m_vector_neutral_charged));   % [m/s] Ion Thermal velocity
v_r  = sqrt(8*Kb*T0./(pi*m_vector_neutral_charged./2));    % [m/s] Thermal velocity with reduced mass

%----------------------------------------------------------------------------------------------------------------------------------------------------
% DIFFUSION COEFFICIENTS and PARAMETERS
% LENNARD-JONES PARAMETERS
LJ_A=1.06036;                             % Coefficient [#]
LJ_B=0.15610;                             % Coefficient [#]
LJ_C=0.19300;                             % Coefficient [#]
LJ_D=0.47635;                             % Coefficient [#]
LJ_E=1.03587;                             % Coefficient [#]
LJ_F=1.52996;                             % Coefficient [#]
LJ_G=1.76474;                             % Coefficient [#]
LJ_H=3.89411;                             % Coefficient [#] 

% ORDER N2 N O O2 NO NO2 N2O CO2 CO CO4 C2O2 C2O3 C2O4 C C2 O4
% FOR NOW N2 N O2 O NO NO2 N2O CO2 CO CO4 C2O2 C C2 
% [m^2] First Lennard Jhones parameter [table B.2 - Bookmatter  Pulverized
% - Coal Combustion], NO2 from (FROM SVEHLA-TRANSPORT COEFFICIENTS FOR THE NASA LEWIS, 1995)
% ALL VALUES FROM https://ntrs.nasa.gov/api/citations/19630012982/downloads/19630012982.pdf
sigma_LJ         = [3.798e-10,3.298e-10,3.467e-10,3.050e-10,3.492e-10,3.765e-10,3.828e-10,3.941e-10,3.690e-10, 0, 0, 0, 0,3.385e-10,3.913e-10, 0];
epsilonOverKb_LJ = [     71.4,     71.4,    106.7,    106.7,    116.7,      210,    232.4,    195.2,     91.7, 0, 0, 0, 0,     30.6,     78.8, 0];   % [K] Second Lennard Jhones parameter 
psi = T0./epsilonOverKb_LJ;             % [#] Reduced temperature 
Omega_D = LJ_A./psi.^LJ_B+LJ_C./exp(LJ_D*psi)+LJ_E./exp(LJ_F*psi)+LJ_G./exp(LJ_H*psi);   % [#] Reduced collision integral 

Comb = nchoosek(1:length(sigma_LJ),2);  % Combination indices
n_vector          = [n_N2+n_N2i n_N+n_Ni n_O+n_Oi+n_On n_O2+n_O2i+n_O2n n_NO+n_NOi+n_NOn n_NO2+n_NO2i+n_NO2n n_N2O+n_N2Oi+n_N2On n_CO2+n_CO2i n_CO+n_COi n_CO4i+n_CO4n n_C2O2i n_C2O3i n_C2O4i n_C+n_Ci n_C2+n_C2i n_O4i+n_O4n];
alphap            = [1.710, 1.100, 0.802, 1.562, 1.698, 2.910, 2.998, 2.507, 4.753, 1.953, 3.866, 4.677, 5.376, 1.76, 2.4, 3.04]*1e-30;  %species polarizability [Å3] from https://cccbdb.nist.gov/pollistx.asp

% Molar fractions  [#]
x_vector = n_vector/ntot;

% Mass fractions  [#]
total_mass_density = total_mass_vector*n_vec;
y_vector = m_vector_neutral_charged.*n_vector./total_mass_density;
        
for i=1:length(Comb)
    psi1 = psi(Comb(i,1)); 
    psi2 = psi(Comb(i,2));
    psi_comb(i) = sqrt(psi1*psi2);       % [#] Reduced temperature
    Omega_D_comb(i) = LJ_A/psi_comb(i)^LJ_B+LJ_C/exp(LJ_D*psi_comb(i))+LJ_E/exp(LJ_F*psi_comb(i))+LJ_G/exp(LJ_H*psi_comb(i));
    
    m1 = m_vector_neutral_charged(Comb(i,1));  sigma1 = sigma_LJ(Comb(i,1));
    m2 = m_vector_neutral_charged(Comb(i,2));  sigma2 = sigma_LJ(Comb(i,2));
    M_comb(i) = 2*m1*m2/(m1+m2);         % Reduced mass
    sigma_comb(i) = sqrt(sigma1*sigma2);
    
    D_comb(i) = 3/8*sqrt(Kb*T0/(pi*M_comb(i)))/(Omega_D_comb(i)*sigma_comb(i)^2)/ntot;
end

% Lennard-Jones Diffusion
if flag_diffusion == 1
    Di = 3/8*sqrt(Kb*T0./(pi*m_vector_neutral_charged))./(sigma_LJ.^2.*Omega_D)./ntot; 
% Langevin approximation
elseif flag_diffusion == 2
    eps_r = m_vector_neutral_charged/4.*v_r.^2;    %eps_r=1/2*reduced_mass*vr^2 where Reduced mass = m_vector/2 
    QL  = q*(pi/2/eps0)^0.5.*(alphap./eps_r).^0.5;
    mui = q./(m_vector_neutral_charged.*ntot.*v_th_vec.*QL*2./3);
    Di  = T0*Kb/q*mui;
elseif flag_diffusion == 0
    % Mixture averaged ion diffusivity 
    for i=1:length(sigma_LJ)
        x_local        = x_vector;
        x_local(i)     = [];
        x_over_D_local = 0;
        [row,columns]  = find(Comb==i);
        row            = sort(row);
        D_local        = D_comb(row);
        
        for j=1:length(row)
            x_over_D_local = x_over_D_local+x_local(j)/D_local(j);            
        end
        x_over_D(i) = x_over_D_local;    
    end    
    Di = (1-y_vector)./x_over_D;   
end

mui = Di*q./(Kb*T0);                                              % [m^2/(V s)] Ion mobility 
lambdai = sqrt(8.0*m_vector_neutral_charged.*Di.*mui./(pi*q));    % [m] Ion mean free path 
vi = q./(m_vector_neutral_charged.*mui);                          % [1/s] Ion collision frequency
ub_vec = sqrt(q*Te./m_vector_neutral_charged);                    % [m/s] Bohm velocity
        
if flag_electronegative == 0
    
    % Mobility and diffusion
    De  = q*Te/(emass*ve);
    mue = q/(emass*ve);

    % Cyclotron frequencies and radii
    omega_C = q*B./m_vector_neutral_charged;
    r_C     = v_th_vec./omega_C;

    % Magnetic factor
    fi = 1./(1+omega_C.^2./vi.^2);
    fe = 1./(1+omega_c^2/ve^2);

    % Mobility and diffusion transversal
    DeT  = De.*fe;
    mueT = mue.*fe;
    DiT  = Di.*fi;
    muiT = mui.*fi;

    % Ambipolar diffusion along and transversal (ions and electrons)
    Da  = (mui.*De+mue*Di)./(mui+mue);                   
    DaT = (muiT.*DeT+mueT*DiT)./(muiT+mueT);            
    DlT = muiT*Te; 
    DRT = exp( (1-lambda).*log(DaT) + lambda.*log(DlT) ); 

    % J1 bessel function of 2.405                                          
    J1       = 5.1911e-1;                                                      
    hL_vec   = 0.86./(3+L/2./lambdai+(0.86*L*ub_vec./pi./Da).^2).^0.5;     
    hR_vec   = 0.8./(4+R./lambdai+(0.8*R*ub_vec./(2.405*J1*DRT)).^2).^0.5;     % Godyak expression
    hCusp    = 0.8./(4 + R./lambdai + (0.8*R*ub_vec./(2.405*J1*Da)).^2).^0.5;  % Godyak expression for the diffusion towards the cusp regions   
    Beta_vec = 1./(7*(1-hR_vec.^(1/6))).*(((1-hR_vec.^(1/6))-1).^7+1);               

    % Compute effective wall recombination area
    if B == 0              % Avoid singularity due to infinite Larmor radii
        A_cusp = 0;                                       % [m^2] Cusp area
    else
        A_cusp = n_cusps*L_cusp*4*sqrt(r_C*r_le);         % [m^2] Cusp area
    end
    Ah_vec   = hL_vec.*(2*Beta_vec*pi*R^2 - Beta_vec*pi*R_foro^2);  % [m^2] Area of diffusion parallel to magnetic field
    Ar_vec   = hR_vec.*(2*pi*R*L - A_cusp) + hCusp.*A_cusp;         % [m^2] Area of diffusion perpendicular to magnetic field
    Aeff_vec = Ah_vec+Ar_vec;

elseif flag_electronegative == 1
    %lambdai from diffusion process
    n_neg       = n_On+n_O2n+n_NOn+n_NO2n+n_N2On+n_CO3n+n_CO4n+n_O3n+n_O4n;
    alpha       = n_neg/n_e;      
    gamma_plus  = Tpos_ion/(Te*q/Kb);
    gamma_minus = (Te*q/Kb)/Tneg_ion;
    tau_vec     = lambdai./v_th_vec;
    omega_C     = q*B./m_vector_neutral_charged;
    fb_vec      = (1+(omega_C.*tau_vec).^2).^-1;

    hL_vec = 0.86*(3+L/2./lambdai+(1+alpha).^0.5.*gamma_plus/5*(L./lambdai).^2).^-0.5*((gamma_minus-1)/(gamma_minus*(1+alpha)^2)+1/gamma_minus)^0.5;    
    hR_vec = 0.80*fb_vec.*(4+R./lambdai+(1+alpha).^0.5*gamma_plus*(R./lambdai).^2).^-0.5*((gamma_minus-1)/(gamma_minus*(1+alpha)^2)+1/gamma_minus)^0.5;
    hCusp  = 0.80*(4+R./lambdai+(1+alpha)^0.5*gamma_plus*(R./lambdai).^2).^-0.5*((gamma_minus-1)/(gamma_minus*(1+alpha)^2)+1/gamma_minus)^0.5;      %obtained from a comparison with the theory used in classic diffusion
    Beta_vec = 1./(7*(1-hR_vec.^(1/6))).*(((1-hR_vec.^(1/6))-1).^7+1);
    r_C      = v_th_vec./omega_C;
    if B == 0              % Avoid singularity due to infinite Larmor radii
        A_cusp = 0;                                       % [m^2] Cusp area
    else
        A_cusp = n_cusps*L_cusp*4*sqrt(r_C*r_le);         % [m^2] Cusp area
    end    
    Ah_vec   = hL_vec.*(2*Beta_vec*pi*R^2 - Beta_vec*pi*R_foro^2);  % [m^2]   Area of diffusion parallel to magnetic field
    Ar_vec   = hR_vec.*(2*pi*R*L - A_cusp) + hCusp.*A_cusp;         % [m^2]   Area of diffusion perpendicular to magnetic field
    Aeff_vec = Ah_vec+Ar_vec;
end
end

