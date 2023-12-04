function launch_sim_cluster_abep
% close all
clear all

%% SET-UP PARPOOL (requires parallel toolbox)
% n_w = 3; %<-------- N° of workers (threads/processors ased on your machine)
% if isempty(gcp('nocreate'))
%     parpool(n_w) % starts parallel. Remember to delete pool before starting new sessions i.e., delete(gcp); parpool(n_w)
% end
% addpath('./gas_data/');
clc
flag_sim_mode ='height';        % 'inflow' for a se mass inflow, 'height' for the atmospheric quantities

switch flag_sim_mode
    case 'height'
        %% INPUT OF THE SIMULATIONS
        Gas_name = 'Mars';             % input: write gas name 'Air','Mars'
        addpath('./gas_data/');
        load(strcat(Gas_name,'.mat'));
        
        
        
        
        %% Function of the altitude 
        for alt=12
        Pw = 1000; % choose deposited power [w] (IMPLEMENT PRODUCT P_ANETNNA EFF_ANTENNA)
        h_data = gas{5,2}; 
        mdot=0;
        H=h_data(alt,1)';
        Thrust=[]; Isp=[]; Eta=[]; Ne=[]; Te=[]; Pres=[]; Iz_ratio=[];
        for jdx = 1:length(H)
            h=H(jdx);
            [vars_out,dens_out] = main_abep(Pw,flag_sim_mode,Gas_name,h);
            thrust   = vars_out(1);
            ne = dens_out(2);
            Thrust   = [Thrust; thrust];
            Ne = [Ne; ne];
        end
        % rimoltiplicare il vettore .*ones(length(Pw),3)
        Output_vars = [Pw*ones(length(Thrust),1) h' Thrust];
        Output_dens = [Ne];
        save('air_data.mat','Output_vars')    % choose height
        
        end
   case 'inflow'
        
        Pw = 200;              % Power provided
        mdot = 20*7.43583e-4*44e-5;            % Massflow at the intake
        
        h = 0;              % Variables necessary to run
        alt = 0;            % but not defined in this case

        Thrust=[]; Isp=[]; Eta=[]; Ne=[]; Te=[]; Pres=[]; Iz_ratio=[];

        [vars_out,dens_out] = main_abep(Pw,flag_sim_mode,Gas_name,mdot);
        thrust   = vars_out(1);
        ne = dens_out(2);
        Thrust   = [Thrust; thrust];
        Ne = [Ne; ne];

end
end