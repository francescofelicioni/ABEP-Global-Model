# ABEP-Global-Model
DO NOT READ IN GITHUB

Global model comprising Earth's atmosphere and Mars' atmosphere gases

Launch_sim_cluster_abep:
Launches the main_abep script defining power and altitude or mass flow at the intake. Outputs all the interesting
variables such as thrust and electron density.

main_abep:
Main body of the Global Model. It solves the number density and light species power balances. At last it uses the Lafleur model to 
evaluate thrust, specific impulse and efficiency.

Reaction_rates:
Generates a cell array with all the reaction rates of each species either with tabulated data or integrating the
Maxwellian distribution.

Input_abep_cluster:
Requires all the inputs necessary to evaluate the plasma diffusion properties and the initial conditions of the differential equations.

plasma_decomposition_abep:
Evaluates the power and density variations due to volume reactions.

INPUTS:
SCRIPT                      CASE        LINE            INPUT
launch_sim_cluster_abep                 12              flag_sim_mode (inflow or height)
                            height      17              Gas_name (Mars or Air)
                            height      25              alt (line correspondant to the altitude in mars.mat or Air.mat)
                            height      26              Pw (power provided to the antenna)
                            inflow      47              Pw (power provided to the antenna)
                            inflow      48              mdot (mass inflow)
main_abep                               32              flag_electronegative (electronegativity)
                                        33              flag_diffusion (diffusion model)
                                        34              flag_atm (mars or earth)
                                        37              L (length of the sicharge chamber)
                                        38              R (Radius of the discharge chamber)
                                        39              B (magnetic field applied)
                                        189             Conc (Concentration of the species for the inflow case)
                                        193             T0 (Initial heavy temperature for the inflow case)
Input_abep_cluster                      5               lambda (type of diffusion)
                                        6               gamma_rec (recombination coefficient)
                                        7               Power_eff (antenna efficiency)
                                        9               ratio (outlet/discharge surfaces ratio)
                                        14              n_cusps (number of magnetic cusps)
                                        16              exp_ratio (throat-to-exit area ratio)
                                        29              Te_0 (initial electron temperature)
