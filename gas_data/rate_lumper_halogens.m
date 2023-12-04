
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             LUMPED MODEL                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates lumped curves of rates of reactions given Te    %
%                                                                         %
% Written the 1/4/2020 by N. Souhair, University of Bologna               %
% Modified the 1/3/2022 by S. Dalle Fabbriche, University of Bologna      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rate_lumper_halogens(Te,Gas_name,flag_cross)
global Reaction Ei gi Te_interp hP c_light e m gas_name gas
gas_name = Gas_name;
Te_interp = Te;

%% load stuff
load(strcat(gas_name,'.mat'));
Gas_model=table2struct(gas{1,2});

Reaction = [Gas_model(:).Reaction];
Ei = [Gas_model(:).Ei];
gi = [Gas_model(:).gi];
hP = 4.135667696e-15; % Plank constant [eV*s]
c_light = 299792458; % speed of light [m/s]
e = 1.6e-19;  % electron charge [C]
m = 9.109e-31; % electron mass [kg]

%% compute K(Te) detailed
K_detailed = compute_K_det(flag_cross);

K_mol_elast  = K_detailed.mol_elast{1,1}; 
K_mol_ion    = K_detailed.mol_ion{1,1};   
K_diss_att   = K_detailed.diss_att{1,1};   
K_mol_diss   = K_detailed.mol_diss{1,1}; 
K_diss_ion   = K_detailed.diss_ion{1,1}; 
K_atom_elast = K_detailed.atom_elast{1,1};
K_atom_exc   = K_detailed.atom_exc{1,1};
K_atom_ion   = K_detailed.atom_ion{1,1};
K_ion_diss   = K_detailed.ion_diss{1,1};
K_detach     = K_detailed.detach{1,1};
K_mol_neutr_overK       = K_detailed.mol_neutr_overK{1,1};
K_atom_neutr_overK      = K_detailed.atom_neutr_overK{1,1};
K_charge_exchange_overK = K_detailed.charge_exchange_overK{1,1};

save('K_detailed.mat','Te_interp','K_detailed');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot_rates; % uncomment for plotting
end

%% Service fun

% compute detailed rate coefficients
function K_detailed = compute_K_det(flag_cross)
global Te_interp gas T0

%% prepare rates of reaction
Gas_rates=gas{2,2};
for i=1:size(Gas_rates,1)
    % convert rates stored from cell to struct of strings (easily accessible and faster)
    % NOTE: check always rate coeff units [m3/s]!!!
    K.([Gas_rates{i,1}]) = Gas_rates{i,2};     
end

% initialize matrix
K_matrix   = NaN(numel(Te_interp),size(Gas_rates,1));  %unique matrix containing in each column i rates corresponding to reaction i


%% compute K(Te) detailed
% ELECTRONIC COLLISION:
if flag_cross == true % flag_cross = 1 get rates from cross-sections
    flag_interp = true;
    
    for i=1:size(Gas_rates,1)-5
        K_matrix(:,i) = rate_solver(Te_interp,K.([Gas_rates{i,1}]),flag_interp);
    end
    
else % flag_cross = 0 get rates from linear interpolation of given rates
    
    for m=1:size(Gas_rates,1)-5
        rates = K.([Gas_rates{m,1}]); 
        K_matrix(:,m) = interp1(rates(:,1),rates(:,2),Te_interp,'linear','extrap');
    end
    
end

%rates evaluation for functions of T0
for i=11:13
    if numel(Gas_rates{i,2})==1
        K_matrix(:,i)=Gas_rates{i,2}*ones(size(Te_interp',1),1);
    else
        rates_inK = K.([Gas_rates{i,1}]);
        K_matrix(:,i) = (interp1(rates_inK(:,1),rates_inK(:,2),T0,'linear','extrap'))*ones(size(Te_interp',1),1);
    end
end
if flag_cross == true % flag_cross = 1 get rates from cross-sections
    flag_interp = true;
    
    for i=14:size(Gas_rates,1)
        K_matrix(:,i) = rate_solver(Te_interp,K.([Gas_rates{i,1}]),flag_interp);
    end
    
else % flag_cross = 0 get rates from linear interpolation of given rates
    
    for m=14:size(Gas_rates,1)
        rates = K.([Gas_rates{m,1}]); 
        K_matrix(:,m) = interp1(rates(:,1),rates(:,2),Te_interp,'linear','extrap');
    end
    
end
%     keyboard    
% store data and save
for i=1:size(Gas_rates,1)
    reaction_rate = K_matrix(:,i);
    K_detailed.([Gas_rates{i,1}]) = {reaction_rate};
    
end
end

function [dir_rate]=rate_solver(T,reaction,flag_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         REACTION RATES SOLVER                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This solver integrates the rate coefficients of reaction given their    %
% relative cross-sections. A maxwellian distribution function is assumed. %
%                                                                         %
% ----------------------------------------------------------------------- %
% Input:                                                                  %
%       - electron temperature range, e.g.:     linspace(0.1,10,100)      %
%       - departure level, e.g.:                'gs'                      %
%       - arrival   level, e.g.:                '1 s 5'                   %
%       - flag interpolation:                   1 if yes, 0 if no         %
%                                                                         %
% Output:                                                                 %
%       - direct rate coefficient                                         %
%       - inverse rate coefficient                                        %
%                                                                         %
% ----------------------------------------------------------------------- %
% Example of utilization:                                                 %
%       Te=linspace(.01,10,1000);                                         %
%       [rd,ri]=rate_solver(Te,'gs','1 s 5',1);                           %
%                                                                         %
% ----------------------------------------------------------------------- %
% Written the 1/11/2020 by E. Majorana                                    %
% Reviewed the 30/12/2020 by E. Majorana and N. Souhair                   %
% Alma Propulsion Laboratory - Università di Bologna                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global e m Cross_section_dat

%% Load cross-sections database and look for input levels
% keyboard

pos=find(reaction == [Cross_section_dat{:,1}]');

%% Pre-pro the cross-sections
% set-up integration range
if flag_interp == 1  % extending cross-section range by interpolation
    eps=0:0.01:max(Cross_section_dat{pos, 2}{:,1}); % energy's integration range [eV]   %%CHIEDI NABIL PERCHE' PARTE DA 0
    sigma=interp1(Cross_section_dat{pos, 2}{:,1},Cross_section_dat{pos, 2}{:,2},eps,'linear','extrap'); % cross-section [m^2]
    for i=1:length(eps)
        if sigma(i)<0
            sigma(i)=0;
        end
    end
elseif flag_interp == 0  % NOT interpolating cross section(just experimental points)
    eps=Cross_section_dat{pos, 2}{:,1}; % energy's integration range [eV]
    sigma=Cross_section_dat{pos, 2}{:,2}; % cross-section [m^2]
end

% Initializing vectors
dir_rate=NaN(1,numel(T));

%% Reaction rate's coefficients integration
for i=1:length(T)
    % Direct reaction rate
    const=(8*e/(m*pi))^(1/2)*1/(T(i))^(3/2);          
    f=sigma.*eps.*exp(-eps/T(i));               
    dir_rate(i)=const*trapz(eps,f);              % Integral derived according to [Bosi, Phd Thesis, Eq. 2.5] [m^3/s]
    
%     plot_EEDF(eps, T, i) % uncomment to check EEDF
end

%% Show rates
% figure();
% box on, hold on, grid minor
% plot(T,dir_rate,T,inv_rate)
% xlabel('Te [eV]');
% ylabel('Rate [m^3/s]');
% legend('Direct reaction rate','Inverse reaction rate','Location','best')
% axis tight
end

function plot_EEDF(eps, T, iter)
figure(1)
box on, hold on
for j=1:length(eps) 
    f1(j)=(exp(-eps(j)/T(iter)))*sqrt((eps(j)));  % distribuzione maxwelliana (eV^-1)   
    f2(j)=f1(j)/sqrt(eps(j));                     % distribuzione maxwelliana (eV^-3/2)   
end
plot(eps,f2,'DisplayName',['run ',num2str(iter)])
xlabel('eps [eV]')
ylabel('EEDF [eV^-3/2]')
set(gca, 'YScale', 'log')
ylim([1e-9 1e0])
xlim([0 160])
legend


% FUN FACT: origine dei termini foo e foobar usate come var metasintattiche
% L'origine del termine inglese è incerta. La forma foobar potrebbe essere nata
% durante la seconda guerra mondiale dal termine FUBAR, acronimo della frase 
% Fucked Up Beyond All Repair ("fottuto oltre ogni (possibile) riparazione"), 
% poi modificato in foobar che in inglese ha una pronuncia simile. Fonte: "wikipedia"
% PS: la quarantena prosegue bene
end



%% plot stuff
function plot_rates
global gas
load('K_detailed.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_detailed=struct2cell(K_detailed);
Gas_rates=gas{2,2};

% plot detailed 
for idx_i=1:(size(K_detailed,1))
   figure(idx_i);
   set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % ingrandisci a tutto schermo 
   hold on
   box on
   KD=K_detailed{idx_i,1}{1,1}(:,1);
   plot(Te_interp,KD);
   lgd=([Gas_rates{idx_i,1}]);
   set(0, 'DefaultLegendInterpreter', 'none')
   legend(lgd,'Location','best');
   title('REACTION RATE COEFFICIENTS','FontSize',14,'FontWeight','bold','Color','k')
   xlabel('Te [eV]','Color','k','FontWeight','normal');
   ylabel('K [m^3/s]','Color','k','FontWeight','normal');
   axis tight
   grid minor
end

end