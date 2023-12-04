function [Rates,R_const] = Reaction_rates(Tg,T)
% Call the functions that give vectors of Reaction Rates in respect to Te
% for every species 
% Input: Neutral gas temperature Tg

global Kb q Gas_name flag_atm
Rates   = cell(95,1);    %last 10 positions are taken by old reactions
load(strcat(Gas_name,'.mat'));

%% NEW MODEL RATES & reaction set
TeK=@(Te_eV) q/Kb*Te_eV;   %[K] Electron temperature
Eps=@(Te_eV) 3/2*Te_eV;    %[eV] Mean electron energy

%% REACTION SET
% 1)  e + Ni + M   ==> N + M
% 2)  e + e + Ni   ==> N + e
% 3)  e + N        ==> Ni + e + e
% 4)  e + N2i      ==> N + N
% 5)  e + N2i + M  ==> N2 + M
% 6)  e + e + N2i  ==> N2 + e
% 7)  e + N2       ==> N2i + e + e
% 8)  e + Oi + M   ==> O + M
% 9)  e + e + Oi   ==> O + e
% 10) e + O        ==> Oi + e + e
% 11) e + O + O2   ==> On + O2
% 12) e + O + O2   ==> O + O2n
% 13) e + O2i      ==> 2O
% 14) e + e + O2i  ==> O2 + e
% 15) e + O2i + M  ==> O2 + M
% 16) e + O2       ==> O + Oi + e + e   
% 17) e + O2       ==> O + O + e
% 18) e + O2       ==> On + O 
% 19) e + O2       ==> O2i + e + e
% 20) e + O2 + O2  ==> O2n + O2
% 21) e + O2 + N2  ==> O2n + N2
% 22) e + NOi      ==> N + O
% 23) e + NOi + M  ==> NO + M
% 24) e + e  + NOi ==> NO + e
% 25) e + NO + M   ==> NOn + e
% 26) e + NO2i     ==> NO + O
% 27) e + NO2 + M  ==> NO2n + M
% 28) e + NO2      ==> On + NO
% 29) e + N2Oi     ==> N2 + O
% 30) e + N2O      ==> On + N2

% 31) Ni + O       ==> Oi + N 
% 32) Ni + O + M   ==> NOi + M
% 33) Ni + On      ==> O  + N
% 34) Ni + N + M   ==> N2i + M
% 35) Ni + NO      ==> NOi + N
% 36) Ni + NO      ==> N2i + O
% 37) Ni + NO      ==> Oi + N2
% 38) Ni + NOn     ==> NO + N
% 39) Ni + O2      ==> NOi + O
% 40) Ni + O2      ==> Oi + NO
% 41) Ni + O2      ==> O2i + N
% 42) Ni + O2n     ==> O2 + N
% 43) Ni + N2O     ==> NOi + N2
% 44) Ni + N2On    ==> N2O + N
% 45) Ni + NO2     ==> NO2i + N
% 46) Ni + NO2     ==> NOi + NO
% 47) Ni + NO2n    ==> NO2 + N

% 48) N + Oi + M   ==> NOi + M
% 49) N + O + M    ==> NO + M
% 50) N + On       ==> NO + e
% 51) N + N + M    ==> N2 + M
% 52) N + N2i      ==> Ni + N2
% 53) N + NOi + M  ==> N2Oi + M
% 54) N + NO       ==> N2 + O
% 55) N + O2i      ==> NOi + O
% 56) N + O2       ==> NO + O
% 57) N + O2n      ==> NO2 + e
% 58) N + NO2      ==> N2O + O
% 59) N + NO2      ==> N2 + O + O
% 60) N + NO2      ==> NO + NO
% 61) N + NO2      ==> N2 + O2
% 62) N + NO2n     ==> N2 + O2 + e

% 63) Oi + O +  M  ==> O2i + M
% 64) Oi + On      ==> O + O
% 65) Oi + NO      ==> NOi + O
% 66) Oi + NO      ==> O2i + N
% 67) Oi + N2 + M  ==> NOi + N + M
% 68) Oi + NOn     ==> NO + O
% 69) Oi + O2      ==> O2i + O
% 70) Oi + O2n     ==> O2 + O
% 71) Oi + N2O     ==> N2Oi + O
% 72) Oi + N2O     ==> NOi + NO
% 73) Oi + N2O     ==> O2i + N2
% 74) Oi + N2On    ==> N2O + O
% 75) Oi + NO2     ==> NOi + O2
% 76) Oi + NO2     ==> NO2i + O
% 77) Oi + NO2n    ==> NO2 + O

% 78) O + On       ==> O2 + e
% 79) O + N2i      ==> Oi + N2
% 80) O + N2i      ==> NOi + N
% 81) O + NO + M   ==> NO2 + M
% 82) O + NOn      ==> On + NO
% 83) O + O2n      ==> On + O2
% 84) O + O + M    ==> O2 + M
% 85) O + NO2      ==> NO + O2

% 86) On + N2      ==> N2O + e
% 87) On + NO      ==> NO2 + e
% 88) On + NO + M  ==> NO2n + M
% 89) On + N2i     ==> O + N2
% 90) On + N2i     ==> O + N + N
% 91) On + NOi     ==> O + O + N
% 92) On + NOi     ==> O + NO
% 93) On + O2i     ==> O + O + O
% 94) On + O2i     ==> O + O2
% 95) On + N2Oi    ==> O + N2O
% 96) On + N2Oi    ==> O + O + N2
% 97) On + N2O     ==> NOn + NO
% 98) On + N2O     ==> N2On + O
% 99) On + NO2i    ==> O + NO2
% 100) On + NO2i   ==> O + N + O2
% 101) On + NO2    ==> NO2n + O

% 102) N2i + NO    ==> NOi + N2
% 103) N2i + NOn   ==> NO + N2
% 104) N2i + NOn   ==> NO + N + N
% 105) N2i + O2    ==> O2i + N2
% 106) N2i + O2n   ==> O2 + N2
% 107) N2i + O2n   ==> O2 + N + N
% 108) N2i + N2O   ==> N2Oi + N2
% 109) N2i + N2O   ==> NOi + N2 + N
% 110) N2i + N2On  ==> N2O + N2
% 111) N2i + N2On  ==> N2O + N + N
% 112) N2i + NO2   ==> NOi + N2O
% 113) N2i + NO2   ==> NO2i + N2
% 114) N2i + NO2n  ==> NO2 + N2
% 115) N2i + NO2n  ==> NO2 + N + N

% 116) N2 + O2i    ==> NOi + NO
% 117) N2 + O2n    ==> N2 + O2 + e

% 118) NOi + NOn   ==> NO + NO
% 119) NOi + NOn   ==> NO + N + O
% 120) NOi + O2n   ==> O2 + N + O
% 121) NOi + O2n   ==> O2 + NO
% 122) NOi + N2On  ==> N2O + NO
% 123) NOi + N2On  ==> N2O + N + O
% 124) NOi + NO2n  ==> NO2 + NO
% 125) NOi + NO2n  ==> NO2 + N + O

% 126) NO + NOn    ==> NO + NO +e
% 127) NO + O2i    ==> NOi + O2
% 128) NO + N2Oi   ==> NOi + N2O
% 129) NO + NO2i   ==> NOi + NO2
% 130) NO + NO2n   ==> NOn + NO2

% 131) NOn + M     ==> NO + M +e
% 132) NOn + O2i   ==> NO + O2
% 133) NOn + O2i   ==> NO + O + O
% 134) NOn + O2    ==> O2n + NO
% 135) NOn + NO2i  ==> NO + NO2

% 137) NOn + N2Oi  ==> NO + N2O
% 138) NOn + N2Oi  ==> NO + N2 + O
% 139) NOn + N2O   ==> NO + N2O + e
% 140) NOn + N2O   ==> NO2n + N2
% 141) NOn + NO2   ==> NO2n + NO

% 142) O2i + O2n   ==> O2 + O2
% 143) O2i + O2n   ==> O2 + O + O
% 144) O2i + N2On  ==> N2O + O2
% 145) O2i + N2On  ==> N2O + O + O
% 146) O2i + NO2   ==> NO2i + O2
% 147) O2i + NO2n  ==> NO2 + O2
% 148) O2i + NO2n  ==> NO2 + O + O

% 149) O2 + O2n    ==> O2 + O2 + e
% 150) O2 + N2Oi   ==> NOi + NO2
% 151) O2 + N2Oi   ==> O2i + N2O

% 152) O2n + N2Oi  ==> O2 + N2O
% 153) O2n + N2Oi  ==> O2 + N2 + O
% 154) O2n + NO2i  ==> O2 + NO2
% 155) O2n + NO2i  ==> O2 + N + O2
% 156) O2n + NO2   ==> NO2n + O2

% 157) N2Oi + N2O  ==> NOi + NO + N2
% 158) N2Oi + N2On ==> N2O + N2O
% 159) N2Oi + N2On ==> N2O + N2 + O
% 160) N2Oi + NO2  ==> NOi + N2 + O2
% 161) N2Oi + NO2  ==> NO2i + N2O
% 162) N2Oi + NO2n ==> NO2 + N2O
% 163) N2Oi + NO2n ==> NO2 + N2 + O

% 164) N2On + NO2i ==> N2O + NO2
% 165) N2On + NO2i ==> N2O + N + O2

% 166) NO2i + NO2n ==> NO2 + NO2
% 167) NO2i + NO2n ==> NO2 + N + O2


% 168) e + N2      ==> N + N + e
% 169) e + N2      ==> e + N2*     
% 170) e + N2      ==> e + N2       
% 171) e + N2i     ==> e + N2i

% 172) e + N       ==> e + N*
% 173) e + N       ==> e + N
% 174) e + Ni      ==> e + Ni

% 175) e + O       ==> e + O*
% 176) e + O       ==> e + O
% 177) e + Oi      ==> e + Oi
% 178) e + O2      ==> e + O2
% 179) e + O2      ==> e + O2v1
% 180) e + O2      ==> e + O2v2
% 181) e + O2      ==> e + O2v3
% 182) e + O2      ==> e + O2e1
% 183) e + O2      ==> e + O2e2
% 184) O2v1 + O2   ==> O2 + O2
% 185) O2v2 + O2   ==> O2 + O2
% 186) O2v3 + O2   ==> O2 + O2



% 187) e + CO2     ==> e + CO2
% 188) e + CO2     ==> e + e + CO2i
% 189) e + CO2     ==> e + e + COi + O
% 190) e + CO2     ==> e + e + Ci + O2
% 191) e + CO2     ==> e + e + Oi + CO
% 192) e + CO2     ==> e + e + O2i + C
% 193) e + CO2     ==> On + CO
% 194) e + CO2     ==> e + CO + O
% 195) e + CO2     ==> e + CO2v1
% 196) e + CO2     ==> e + CO2v2
% 197) e + CO2     ==> e + CO2v3
% 198) e + CO2     ==> e + CO2v4
% 199) e + CO2     ==> e + CO2e1
% 200) e + CO2     ==> e + CO2e2

% 201) e + CO      ==> e + CO
% 202) e + CO      ==> e + e + COi
% 203) e + CO      ==> e + e + Ci + O
% 204) e + CO      ==> e + e + C + Oi
% 205) e + CO      ==> On + C
% 206) e + CO      ==> e + C + O
% 207) e + CO      ==> e + COv1
% 208) e + CO      ==> e + COe1
% 209) e + CO      ==> e + COe2
% 210) e + CO      ==> e + COe3
% 211) e + CO      ==> e + COe4

% 212) e + C       ==> e + C
% 213) e + C       ==> e + e + Ci
% 214) e + C2      ==> e + C2
% 215) e + C2      ==> e + C + C
% 216) e + C2      ==> e + e + C2i
% 217) e + O3      ==> e + O3
% 218) e + O3      ==> On + O2
% 219) e + O3      ==> O + O2n
% 220) e + CO2i    ==> CO + O
% 221) e + CO2i    ==> C + O2
% 222) e + CO4i    ==> CO2 + O2
% 223) e + COi     ==> C + O
% 224) e + C2O2i   ==> CO + CO
% 225) e + C2O3i   ==> CO2 + CO
% 226) e + C2O4i   ==> CO2 + CO2
% 227) e + C2i     ==> C + C
% 228) e + O2 + M  ==> 02n + M
% 229) e + O3 + M  ==> 03n + M
% 230) e + O + M   ==> 0n + M
% 231) e + O2i + M ==> O2 + M
% 232) e + O4i     ==> O2 + O2

%% 233) CO2 + M     ==> CO + O + M       R = 3.91e-16*exp(-49430/Tg)
% 233) O + CO2     ==> CO + O2
% 234) C + CO2     ==> CO + CO
% 235) O + CO + M  ==> CO2 + M
% 236) O2 + CO    ==> CO2 + O 
% 237) O3 + CO    ==> CO2 + O2
% 238) C + CO + M  ==> C2O + M
% 239) O2 + C      ==> CO + O
% 240) O + C + M   ==> CO + M
% 241) O + C2O     ==> CO + CO
% 242) O2 + C2O    ==> CO2 + CO
% 243) O + O3      ==> O2 + O2
% 244) O3 + M      ==> O2 + O + M
% 245) O + O2 + M  ==> O3 + M

% 246) O2i + CO2 + M   ==> CO4i + M
% 247) Oi + CO2        ==> O2i + CO
% 248) Oi + CO2        ==> CO2i + O
% 249) Ci + CO2        ==> COi + COi
% 250) COi + CO2       ==> CO2i + CO
% 251) On + CO2 + M    ==> CO3n + M
% 252) O2n + CO2 + M   ==> CO4n + M
% 253) O3n + CO2       ==> O2 + CO3n
% 254) O4n + CO2       ==> CO4n + O2
% 255) CO2i + CO2 + M  ==> C2O4i + M
% 256) Oi + CO         ==> COi + O
% 257) On + CO         ==> CO2 + e
% 258) CO3n + CO       ==> 2CO2 + e
% 259) C2O3i + CO      ==> CO2 + C2O2i
% 260) C2O4i + CO      ==> CO2 + C2O3i
% 261) C2O3i + CO + M  ==> C2O2i + CO2 + M
% 261) C2O4i + CO + M  ==> C2O3i + CO2 + M

% 263) Ci + CO         ==> COi + C
% 264) COi + C         ==> CO + Ci
% 265) O2i + C         ==> COi + O
% 266) O2i + C         ==> Ci + O2
% 267) C2i + C         ==> C2 + Ci
% 268) O + CO2i        ==> O2i + CO
% 269) O + CO2i        ==> Oi + CO2
% 270) O2 + CO2i        ==> O2i + CO2
% 271) CO3n + CO2i     ==> 2CO2 + O
% 272) CO4n + CO2i     ==> 2CO2 + O2
% 273) O2n + CO2i      ==> CO + O2 + O
% 274) O + COi         ==> CO + Oi
% 275) O2 + COi        ==> O2i + CO
% 276) O2 + C2O2i      ==> CO + CO + O2i
% 277) C2O2i + M       ==> COi + CO + M
% 278) CO3n + C2O2i    ==> CO2 + 2CO + O
% 279) CO4n + C2O2i    ==> CO2 + 2CO + O2
% 280) O2n + C2O2i     ==> 2CO + O2
% 281) CO3n + C2O3i    ==> 2CO2 + CO + O
% 282) CO4n + C2O3i    ==> 2CO2 + CO + O2
% 283) O2n + C2O3i     ==> CO2 + CO + O2
% 284) C2O4i + M       ==> CO2i + CO2 + M
% 285) CO3n + C2O4i    ==> 3CO2 + O
% 286) CO4n + C2O4i    ==> 3CO2 + O2
% 287) O2n + C2O4i     ==> 2CO2 + O2
% 288) O2i + CO3n      ==> CO2 + O2 + O
% 289) O2i + CO4n      ==> CO2 + O2 + O2
% 290) O + CO3n        ==> CO2 + O2n
% 291) O + CO4n        ==> CO3n + O2
% 292) O + CO4n        ==> CO2 + O2 + On
% 293) O + CO4n        ==> CO2 + O3n
% 294) O3 + CO4n       ==> CO2 + O3n + O2
% 295) O2 + Ci         ==> CO + Oi
% 296) O2 + Ci         ==> COi + O
% 297) O2i + O2 + M    ==> O4i + M
% 298) O2n + O2 + M    ==> O4n + M
% 299) On + O2         ==> O3 + e
% 300) On + O2 + M     ==> O3n + M 
% 301) On + O3         ==> O3n + O
% 302) On + O3         ==> O2 + O2 + e
% 303) O2n + O3        ==> O3n + O2
% 304) O3n + O3        ==> O2 + O2 + O2 + e
% 305) Oi + O3         ==> O2i + O2
% 306) Oi + O + M      ==> O2i + M
% 307) O2n + O         ==> O3 + e
% 308) O3n + O         ==> O3 + On
% 309) O3n + O         ==> O2 + O2 + e
% 310) O3n + O         ==> O2n + O2
% 311) O4n + O         ==> O3n + O2
% 312) O4n + O         ==> On + O2 + O2
% 313) O4i + O         ==> O2i + O3
% 314) Oi + O2n + M    ==> O3 + M
% 315) O2i + O2n + M   ==> O2 + O2 + M
% 316) O2n + O2        ==> O2 + O2 + e
% 317) O2i + O3n       ==> O2 + O3
% 318) O2i + O3n       ==> O + O + O3
% 319) Oi + O3n        ==> O3 + O
% 320) O2 + O3n        ==> O2 + O3 + e
% 321) O3n + M         ==> O3 + e
% 322) Oi + On + M     ==> O2 + M
% 323) O2i + On + M    ==> O3 + M
% 324) M + On          ==> O + M + e
% 325) M + O4n         ==> O2n + O2 + M
% 326) O4i + M         ==> O2i + O2 + M
% 327) CO2v1 + CO2     ==> CO2 + CO2
% 328) CO2v1 + CO      ==> CO2 + CO
% 329) CO2v1 + O2      ==> CO2 + O2
% 330) CO2v2 + CO2     ==> CO2 + CO2
% 331) CO2v2 + CO      ==> CO2 + CO
% 332) CO2v2 + O2      ==> CO2 + O2
% 333) CO2v2 + CO2     ==> CO2v1 + CO2
% 334) CO2v2 + CO      ==> CO2v1 + CO
% 335) CO2v2 + O2      ==> CO2v1 + O2

% 336) CO2v3 + CO2     ==> CO2v2 + CO2
% 337) CO2v3 + CO      ==> CO2v2 + CO
% 338) CO2v3 + O2      ==> CO2v2 + O2
% 339) CO2v3 + CO2     ==> CO2v4 + CO2
% 340) CO2v3 + CO      ==> CO2v4 + CO
% 341) CO2v3 + O2      ==> CO2v4 + O2
% 342) CO2v3 + CO2     ==> CO2v1 + CO2v2
% 343) CO2v3 + CO2     ==> CO2v1 + CO2
% 344) CO2v3 + CO      ==> CO2v1 + CO
% 345) CO2v3 + O2      ==> CO2v1 + O2
% 346) CO2v4 + CO2     ==> CO2v2 + CO2
% 347) CO2v4 + CO      ==> CO2v2 + CO
% 348) CO2v4 + O2      ==> CO2v2 + O2
% 349) CO2v4 + O2      ==> CO2v1 + CO2
% 350) CO2v4 + CO      ==> CO2v1 + CO
% 351) CO2v4 + O2      ==> CO2v1 + O2
% 352) COv1 + CO2      ==> CO + CO2
% 353) COv1 + CO       ==> CO + CO
% 354) COv1 + O2       ==> CO + O2
% 355) Ov1 + CO        ==> O + CO2
% 356) O2v1 + CO       ==> O2 + CO
% 357) O2v2 + CO2      ==> O2 + CO2
% 358) O2v2 + CO       ==> O2 + CO2
% 359) O2v3 + CO2      ==> O2 + CO2
% 360) O2v3 + CO       ==> O2 + CO


%% NEW MODEL RATES
TeK=@(Te_eV) q/Kb*Te_eV;   %[K] Electron temperature
Eps=@(Te_eV) 3/2*Te_eV;    %[eV] Mean electron energy


%------------------------------------------------------------------
% REACTIONS 1-30
Rates{1,1} = @(Te_eV) 3.12e-35/(TeK(Te_eV)^1.5);
Rates{2,1} = @(Te_eV) 1e-31*(Tg/TeK(Te_eV))^4.5;
% R{3,1} = @(Te_eV) 1.45e-17*Eps(Te_eV)^2.58*exp(-8.54/Eps(Te_eV));    OLD MODEL
Rates{4,1} = @(Te_eV) 2.8e-13*(300/TeK(Te_eV))^0.5;
Rates{5,1} = @(Te_eV) 3.12e-35/(TeK(Te_eV)^1.5);
Rates{6,1} = @(Te_eV) 1e-31*(Tg/TeK(Te_eV))^4.5;
% R{7,1} EEDF-CALCULATION     N2 IONIZATION     OLD MODEL
Rates{8,1} = @(Te_eV) 3.12e-35/(TeK(Te_eV)^1.5);
Rates{9,1} = @(Te_eV) 1e-31*(Tg/TeK(Te_eV))^4.5;
% R{10,1} = @(Te_eV) 4.75e-15*Eps(Te_eV)^0.61*exp(-22.1/Eps(Te_eV));   OLD MODEL

Rates{11,1} = @(Te_eV) 1e-43;            
Rates{12,1} = @(Te_eV) 1e-43;                  
Rates{13,1} = @(Te_eV) 2e-13*(300/TeK(Te_eV));
Rates{14,1} = @(Te_eV) 1e-31*(Tg/TeK(Te_eV))^4.5;
Rates{15,1} = @(Te_eV) 3.12e-35/(TeK(Te_eV)^1.5);
% R{16,1} = EEDF-CALCULATION   O2 DISSOCIATIVE IONIZATION
Rates{17,1} = @(Te_eV) 2.03e-14*Eps(Te_eV)^(-0.1)*exp(-8.47/Eps(Te_eV));                          
% R{18,1} = EEDF-CALCULATION   O2 DISSOCIATIVE ATTACHMENT             
% R{19,1} = EEDF-CALCULATION   O2 IONIZATION                  
Rates{20,1} = @(Te_eV) 1.4e-41*Tg/TeK(Te_eV)*exp(-600/Tg)*exp(700*(TeK(Te_eV)-Tg)/(Tg*TeK(Te_eV)));

Rates{21,1} = @(Te_eV) 1.4e-43*(Tg/TeK(Te_eV))^2*exp(-70/Tg)*exp(1500*(TeK(Te_eV)-Tg)/(Tg*TeK(Te_eV)));
Rates{22,1} = @(Te_eV) 1.07e-11/(TeK(Te_eV)^0.85);
Rates{23,1} = @(Te_eV) 3.12e-35/(TeK(Te_eV)^1.5);
Rates{24,1} = @(Te_eV) 1e-31*(Tg/TeK(Te_eV))^4.5;
Rates{25,1} = @(Te_eV) 8e-43;
Rates{26,1} = @(Te_eV) 3.46e-12/(TeK(Te_eV)^0.5);
Rates{27,1} = @(Te_eV) 1.5e-42;
Rates{28,1} = @(Te_eV) 1e-17;
Rates{29,1} = @(Te_eV) 3.46e-12/(TeK(Te_eV)^0.5);
Rates{30,1} = @(Te_eV) 2e-16;


%--------------------------------------------------------------------------
% REACTIONS 31-167
R_const(1)  = 1e-18;                    
R_const(2)  = 1e-41;                    
R_const(3)  = 2e-13*(300/Tg)^0.5;      
R_const(4)  = 1e-41;                    
R_const(5)  = 4.72e-16; 
R_const(6)  = 8.33e-17; 
R_const(7)  = 1e-18; 
R_const(8)  = 2e-13*(300/Tg)^0.5; 
R_const(9)  = 2.7e-16;             
R_const(10) = 2.8e-17;              

R_const(11) = 3e-16;                 
R_const(12) = 2e-13*(300/Tg)^0.5; 
R_const(13) = 5.5e-16;      
R_const(14) = 2e-13*(300/Tg)^0.5; 
R_const(15) = 3e-16;      
R_const(16) = 5e-16;      
R_const(17) = 2e-13*(300/Tg)^0.5;     
R_const(18) = 1e-41;
R_const(19) = 6.3e-45*exp(140/Tg);   
R_const(20) = 2.6e-16;               

R_const(21) = 8.3e-46*exp(500/Tg);   
R_const(22) = 1e-18;              
R_const(23) = 1e-41*(300/Tg);
R_const(24) = 2.1e-17*exp(100/Tg);
R_const(25) = 1.5e-16;                
R_const(26) = 1.5e-17*exp(-3600/Tg);   
R_const(27) = 5e-16;   
R_const(28) = 5.8e-18*exp(220/Tg); 
R_const(29) = 9.1e-19; 
R_const(30) = 6e-19; 

R_const(31) = 7e-19; 
R_const(32) = 1e-18; 
R_const(33) = 1e-41;                  
R_const(34) = 2e-13*(300/Tg)^0.5;    
R_const(35) = 1e-18; 
R_const(36) = 3e-18; 
R_const(37) = 6e-41*(300/Tg)^2; 
R_const(38) = 2e-13*(300/Tg)^0.5;  
R_const(39) = 2.1e-17*(300/Tg)^0.5;    
R_const(40) = 2e-13*(300/Tg)^0.5;    

R_const(41) = 6.3e-16; 
R_const(42) = 2.3e-16; 
R_const(43) = 2e-17; 
R_const(44) = 2e-13*(300/Tg)^0.5;  
R_const(45) = 5e-16; 
R_const(46) = 1.6e-15; 
R_const(47) = 2e-13*(300/Tg)^0.5;  
R_const(48) = 1.4e-16;                
R_const(49) = 1e-17*(300/Tg)^0.5;     
R_const(50) = 1.4e-16;                

R_const(51) = 1e-43*(300/Tg)^1.6; 
R_const(52) = 3e-16; 
R_const(53) = 3.3e-16;          
R_const(54) = 3.2e-47*exp(900/Tg);   
R_const(55) = 6.5e-18*exp(120/Tg);
R_const(56) = 1e-18; 
R_const(57) = 2.6e-16; 
R_const(58) = 1e-41; 
R_const(59) = 2e-13*(300/Tg)^0.5;     
R_const(60) = 1e-13;    

R_const(61) = 1e-13; 
R_const(62) = 2e-13*(300/Tg)^0.5;  
R_const(63) = 1e-13;                  
R_const(64) = 2e-13*(300/Tg)^0.5;     
R_const(65) = 2e-13*(300/Tg)^0.5;  
R_const(66) = 1e-13;   
R_const(67) = 2e-16;  
R_const(68) = 2e-18; 
R_const(69) = 2e-13*(300/Tg)^0.5;
R_const(70) = 1e-13; 

R_const(71) = 1.2e-15; 
R_const(72) = 3.9e-16; 
R_const(73) = 2e-13*(300/Tg)^0.5;
R_const(74) = 1e-13;
R_const(75) = 5e-17;                  
R_const(76) = 2e-13*(300/Tg)^0.5;    
R_const(77) = 1e-13;  
R_const(78) = 6e-16;
R_const(79) = 4e-16;
R_const(80) = 2e-13*(300/Tg)^0.5;

R_const(81) = 1e-13;
R_const(82) = 5e-17;
R_const(83) = 3e-16;
R_const(84) = 2e-13*(300/Tg)^0.5;
R_const(85) = 1e-13;
R_const(86) = 1e-23;            
R_const(87) = 1.9e-18*(Tg/300)^0.5*exp(-4990/Tg);
R_const(88) = 2e-13*(300/Tg)^0.5;
R_const(89) = 1e-13;
R_const(90) = 1e-13;

R_const(91) = 2e-13*(300/Tg)^0.5;
R_const(92) = 2e-13*(300/Tg)^0.5;
R_const(93) = 1e-13;
R_const(94) = 2e-13*(300/Tg)^0.5;
R_const(95) = 1e-13;
R_const(96) = 5e-18;
R_const(97) = 4.6e-16;
R_const(98) = 2.3e-16;
R_const(99) = 2.75e-16;
R_const(100) = 2.75e-16;

R_const(101) = 2.4e-19;
R_const(102) = 2e-13*(300/Tg)^0.5;
R_const(103) = 1e-13;
R_const(104) = 5e-16;
R_const(105) = 2e-13*(300/Tg)^0.5;
R_const(106) = 1e-13;
R_const(107) = 2e-13*(300/Tg)^0.5;
R_const(108) = 1e-13;
R_const(109) = 5.1e-18;
R_const(110) = 2.8e-20;

R_const(111) = 3e-16;
R_const(112) = 2e-13*(300/Tg)^0.5;      
R_const(113) = 1e-13;                   
R_const(114) = 2e-13*(300/Tg)^0.5;    
R_const(115) = 1e-13; 
R_const(116) = 6.6e-16;
R_const(117) = 2e-13*(300/Tg)^0.5;
R_const(118) = 1e-13; 
R_const(119) = 2.7e-16*(Tg/300)^0.5*exp(-5590/Tg);
R_const(120) = 4.59e-17; 

R_const(121) = 2.24e-16; 
R_const(122) = 2e-13*(300/Tg)^0.5;      
R_const(123) = 1e-13;                   
R_const(124) = 2e-13*(300/Tg)^0.5;    
R_const(125) = 1e-13; 
R_const(126) = 7e-16; 
R_const(127) = 1.2e-17; 
R_const(128) = 2e-13*(300/Tg)^0.5;      
R_const(129) = 1e-13; 
R_const(130) = 4.29e-16; 

R_const(131) = 2.21e-16; 
R_const(132) = 2e-13*(300/Tg)^0.5;      
R_const(133) = 1e-13; 
R_const(134) = 2e-13*(300/Tg)^0.5;      
R_const(135) = 1e-13;
R_const(136) = 2e-13*(300/Tg)^0.5;      
R_const(137) = 1e-13;


%--------------------------------------------------------------------------
% REACTONS 168-183
% Rates{31,1} = EEDF-CALCULATION    
% Rates{32,1} = EEDF-CALCULATION    
% Rates{33,1} = EEDF-CALCULATION     
% Rates{34,1} = EEDF-CALCULATION   
% Rates{35,1} = EEDF-CALCULATION
% Rates{36,1} = EEDF-CALCULATION
% Rates{37,1} = EEDF-CALCULATION
% Rates{38,1} = EEDF-CALCULATION
% Rates{39,1} = EEDF-CALCULATION
% Rates{40,1} = EEDF-CALCULATION
% Rates{41,1} = EEDF-CALCULATION
% Rates{42,1} = EEDF-CALCULATION    O2 vibrational excitation v1
% Rates{43,1} = EEDF-CALCULATION    O2 vibrational excitation v2
% Rates{44,1} = EEDF-CALCULATION    O2 vibrational excitation ve
% Rates{45,1} = EEDF-CALCULATION    O2 electronic excitation e1
% Rates{46,1} = EEDF-CALCULATION    O2 electronic excitation e2

switch flag_atm
    case 'earth'
        %------------------------------------------------------------------
        % REACTIONS 184-232
        %  = 0  for case 'earth'

        % REACTIONS 233-360
        R_const(138:265)= 0;


    case 'mars'
        %------------------------------------------------------------------
        % REACTIONS 184-232
        Rates{47,1} = @(Te_eV) 2.52*10^(-23);
        Rates{48,1} = @(Te_eV) 2.52*10^(-23);
        Rates{49,1} = @(Te_eV) 2.52*10^(-23);
        % Rates{50,1} = EEDF-CALCULATION    CO2 momentum transfer                                   187% 
        % Rates{51,1} = EEDF-CALCULATION    CO2 ionization                                          188%
        % Rates{52,1} = EEDF-CALCULATION    CO2 dissociative ionization                             189%
        % Rates{53,1} = EEDF-CALCULATION    CO2 dissociative ionization                             190%
        % Rates{54,1} = EEDF-CALCULATION    CO2 dissociative ionization                             191%
        % Rates{55,1} = EEDF-CALCULATION    CO2 dissociative ionization                             192%
        % Rates{56,1} = EEDF-CALCULATION    CO2 dissociative attachment                             193%
        Rates{57,1} = @(Te_ev) 1.6e-17;                                                            %194%
        % Rates{58,1} = EEDF-CALCULATION    CO2 vibrational excitation v1                           195%
        % Rates{59,1} = EEDF-CALCULATION    CO2 vibrational excitation v2                           196%
        % Rates{60,1} = EEDF-CALCULATION    CO2 vibrational excitation v3                           197%
        % Rates{61,1} = EEDF-CALCULATION    CO2 vibrational excitation v4                           198%
        % Rates{62,1} = EEDF-CALCULATION    CO2 electronic excitation e1                            199%
        % Rates{63,1} = EEDF-CALCULATION    CO2 electronic excitation e2                            200%
        % Rates{64,1} = EEDF-CALCULATION    CO momentum transfer                                    201%
        % Rates{65,1} = EEDF-CALCULATION    CO ionization                                          %202%
        % Rates{66,1} = EEDF-CALCULATION    CO dissociative ionization                             %203%
        % Rates{67,1} = EEDF-CALCULATION    CO dissociative ionization                             %204%
        % Rates{68,1} = EEDF-CALCULATION    CO dissociative attachment                             %205%
        % Rates{69,1} = EEDF-CALCULATION    CO dissociation                                        %206%
        % Rates{70,1} = EEDF-CALCULATION    CO vibrational excitation v1                            207%
        % Rates{71,1} = EEDF-CALCULATION    CO electronic excitation e1                             208%
        % Rates{72,1} = EEDF-CALCULATION    CO electronic excitation e2                             209%
        % Rates{73,1} = EEDF-CALCULATION    CO electronic excitation e3                             210%
        % Rates{74,1} = EEDF-CALCULATION    CO electronic excitation e4                             211%
        % Rates{75,1} = EEDF-CALCULATION    C momentum transfer                                     212%
        % Rates{76,1} = EEDF-CALCULATION    C dissociative ionization                               213%
        % Rates{77,1} = EEDF-CALCULATION    C2 momentum transfer                                    214%
        % Rates{78,1} = EEDF-CALCULATION    C2 dissociation                                         215%
        % Rates{79,1} = EEDF-CALCULATION    C2 ionization                                           216%
        % Rates{80,1} = EEDF-CALCULATION    O3 momentum transfer                                    217%
        Rates{81,1} = @(Te_eV) 10^(-21);                                                           %218%
        Rates{82,1} = @(Te_eV) 10^(-23);                                                           %219%
        Rates{83,1} = @(Te_eV) 2*10^(-11)*Tg^(-1)*Te_eV^(-0.5);                                    %220%
        Rates{84,1} = @(Te_eV) 3.94*10^(-13)*Te_eV^(-0.4);                                         %221%
        Rates{85,1} = @(Te_eV) 1.61*10^(-13)*Te_eV^(-0.5);%
        Rates{86,1} = @(Te_eV) 3.68*10^(-14)*Te_eV^(-0.55);%
        Rates{87,1} = @(Te_eV) 4.0*10^(-13)*Te_eV^(-0.34);%
        Rates{88,1} = @(Te_eV) 5.4*10^(-14)*Te_eV^(-0.7);%
        Rates{89,1} = @(Te_eV) 2.0*10^(-11)*Tg^(-1)*Te_eV^(-0.5);%
        Rates{90,1} = @(Te_eV) 1.79*10^(-14)*Te_eV^(-0.5);%
        Rates{91,1} = @(Te_eV) 3*10^(-42);%
        Rates{92,1} = @(Te_eV) 5*10^(-43)*Te_eV^(-0.5);%
        Rates{93,1} = @(Te_eV) 10^(-43);%
        Rates{94,1} = @(Te_eV) 10^(-38);%
        Rates{95,1} = @(Te_eV) 2.25*10^(-13)*Te_eV^(-0.5);%

        %------------------------------------------------------------------
        % REACTIONS 233-360
        R_const(138)  = 2.8*10^(-17)*exp(-26500/Tg); 
        R_const(139)  = 1.0*10^(-21);
        R_const(140)  = 8.2*10^(-46)*exp(-1510/Tg);
        R_const(141)  = 4.2*10^(-18)*exp(-24000/Tg);
        R_const(142)  = 4.0*10^(-31);
        R_const(143)  = 6.5*10^(-44);
        R_const(144)  = 3.0*10^(-17);
        R_const(145)  = 2.14*10^(-41)*(Tg/300)^(-3.08)*exp(-2114/Tg);
        R_const(146)  = 5.0*10^(-17);
        R_const(147)  = 3.3*10^(-19);
        % 243
        R_const(148)  = 3.1*10^(-20)*Tg^(0.75)*exp(-1575/Tg); 
        R_const(149)  = 4.12*10^(-16)*exp(-11430/Tg);
        R_const(150)  = 6.11e-46*(Tg/300)^(-2.6);
        R_const(151)  = 2.3*10^(-41);
        R_const(152)  = 9.4*10^(-16);
        R_const(153)  = 4.5*10^(-16);
        R_const(154)  = 1.1*10^(-15);
        R_const(155)  = 1.0*10^(-15);
        R_const(156)  = 9.0*10^(-41);
        R_const(157)  = 1.0*10^(-41);
        % 253
        R_const(158)  = 5.5*10^(-16); 
        R_const(159)  = 4.8*10^(-16);
        R_const(160)  = 3.0*10^(-40);
        R_const(161)  = 4.9*10^(-18)*(Tg/300)^(0.5)*exp(-4580/Tg);
        R_const(162)  = 5.5*10^(-16);
        R_const(163)  = 5*10^(-19);
        R_const(164)  = 1.1*10^(-15);
        R_const(165)  = 9.0*10^(-16);
        R_const(166)  = 2.6*10^(-38);
        R_const(167)  = 4.2*10^(-38);
        % 263
        R_const(168)  = 5.0*10^(-19); 
        R_const(169)  = 1.1*10^(-16);
        R_const(170)  = 5.2*10^(-17);
        R_const(171)  = 5.2*10^(-17);
        R_const(172)  = 1.1*10^(-16);
        R_const(173)  = 1.64*10^(-16);
        R_const(174)  = 9.62*10^(-17);
        R_const(175)  = 5.3*10^(-17);
        R_const(176)  = 5*10^(-13);
        R_const(177)  = 5*10^(-13);
        % 273
        R_const(178)  = 6*10^(-13);  
        R_const(179)  = 1.4*10^(-16); 
        R_const(180)  = 1.2*10^(-16);
        R_const(181)  = 5.0*10^(-18);
        R_const(182)  = 1.0*10^(-18);
        R_const(183)  = 5.0*10^(-13);
        R_const(184)  = 5.0*10^(-13);
        R_const(185)  = 6.0*10^(-13);
        R_const(186)  = 5.0*10^(-13);
        R_const(187)  = 5.0*10^(-13);
        % 283
        R_const(188)  = 6.0*10^(-13);  
        R_const(189)  = 1.0*10^(-20);
        R_const(190)  = 5.0*10^(-13);
        R_const(191)  = 5.0*10^(-13);
        R_const(192)  = 6.0*10^(-13);
        R_const(193)  = 3*10^(-13);
        R_const(194)  = 3*10^(-13);
        R_const(195)  = 8*10^(-13);
        R_const(196)  = 1.1*10^(-16);
        R_const(197)  = 1.4*10^(-17);
        %293
        R_const(198)  = 1.4*10^(-16); 
        R_const(199)  = 1.3*10^(-16);
        R_const(200)  = 6.2*10^(-16);
        R_const(201)  = 3.8*10^(-16);
        R_const(202)  = 2.4*10^(-42);
        R_const(203)  = 3.5*10^(-43);
        R_const(204)  = 1*10^(-18);
        R_const(205)  = 3.0*10^(-42)*(Tg/300)^(-1);
        R_const(206)  = 8*10^(-16);
        R_const(207)  = 3.0*10^(-16);
        % 303
        R_const(208)  = 4.0*10^(-16); 
        R_const(209)  = 3.0*10^(-16);
        R_const(210)  = 1.0*10^(-16);
        R_const(211)  = 1.0*10^(-41);
        R_const(212)  = 3.3*10^(-16);
        R_const(213)  = 1.0*10^(-19);
        R_const(214)  = 1.0*10^(-19);
        R_const(215)  = 2.5*10^(-16);
        R_const(216)  = 4.0*10^(-16);
        R_const(217)  = 3.0*10^(-16);
        % 313
        R_const(218)  = 3.0*10^(-16);
        R_const(219)  = 2.0*10^(-37);
        R_const(220)  = 2.0*10^(-37);
        R_const(221)  = 2.18*10^(-24);
        R_const(222)  = 2.0*10^(-13);
        R_const(223)  = 1.0*10^(-13);
        R_const(224)  = 1.0*10^(-13);
        R_const(225)  = 2.3*10^(-17);
        R_const(226)  = 2.3*10^(-17);
        R_const(227)  = 2.0*10^(-37);
        % 323
        R_const(228)  = 2.0*10^(-37);
        R_const(229)  = 4*10^(-18);
        R_const(230)  = 4*10^(-18);
        R_const(231)  = 1.73*10^(-19);
        R_const(232)  = 1.07*10^(-20); 
        R_const(233)  = 7.28*10^(-21);
        R_const(234)  = 7.28*10^(-21);
        R_const(235)  = 9.00*10^(-24);
        R_const(236)  = 2.79*10^(-23);
        R_const(237) = 2.79*10^(-23);
        % 333
        R_const(238) = 2.90*10^(-20);
        R_const(239) = 2.03*10^(-20);
        R_const(240) = 2.03*10^(-20);
        R_const(241) = 7.72*10^(-22);
        R_const(242) = 2.32*10^(-22);
        R_const(243) = 3.09*10^(-22);
        R_const(244) = 6.05*10^(-21);
        R_const(245) = 1.81*10^(-21);
        R_const(246) = 2.42*10^(-21);
        R_const(247) = 2.42*10^(-21);
        %343
        R_const(248) = 1.70*10^(-24);
        R_const(249) = 5.10*10^(-25);
        R_const(250) = 6.80*10^(-25);
        R_const(251) = 4.33*10^(-20);
        R_const(252) = 3.03*10^(-20);
        R_const(253) = 3.03*10^(-20);
        R_const(254) = 9.08*10^(-24);
        R_const(255) = 6.18*10^(-21);
        R_const(256) = 6.18*10^(-21);
        R_const(257) = 1.34*10^(-29);
        %353
        R_const(258) = 1.34*10^(-29);
        R_const(259) = 4.78*10^(-30);
        R_const(260) = 7.55*10^(-29);
        R_const(261) = 2.52*10^(-29);
        R_const(262) = 7.55*10^(-29);
        R_const(263) = 2.52*10^(-29);
        R_const(264) = 7.55*10^(-29);
        R_const(265) = 2.52*10^(-29);
end



Rrate_data = gas{2,2}{1,2};
for p=1:length(Rates)
    check = Rates{p,1};
    if isempty(check) == 1
        Reaction_rate = Rrate_data{p,2};
        if isempty(Reaction_rate) == 1
            Csection_data = Cross_section_dat{p,2};
            Reaction_rate_maxwell  = rate_solver_Maxwell(Csection_data,T,0);
            Rates{p,1}    = griddedInterpolant(T,Reaction_rate_maxwell,'spline');
        elseif Reaction_rate == 0
            Rates{p,1}    = @(Te_eV) 0;
        else
            Rates{p,1}    = griddedInterpolant(Reaction_rate(:,1),Reaction_rate(:,2),'linear');
        end
    end
end

end

function [dir_rate]=rate_solver_Maxwell(C_cross,T,flag_interp)
% Function evaluating thr reaction rates assuming a maxwellian distribution
% and intregrating the cross-sections
global q emass

%% Load cross-sections database and look for input levels
% C_O2ion = readmatrix('O2_ionization.txt');

%% Pre-pro the cross-sections
% set-up integration range

if flag_interp == 1     % extending cross-section range by interpolation
     eps  = min(C_cross(:,1)):0.01:max(C_cross(:,1));     % energy's integration range [eV]
     sigma= interp1(C_cross(:,1),C_cross(:,2),eps,'linear','extrap'); % cross-section [m^2]
    for i=1:length(eps)
        if sigma(i)<0
            sigma(i)=0;
        end
    end
elseif flag_interp == 0  % NOT interpolating cross section(just experimental points)
    eps  = C_cross(:,1);   % energy's integration range [eV]
    sigma= C_cross(:,2);   % cross-section [m^2]
end

% Initializing vectors
dir_rate=NaN(1,numel(T));

%% Reaction rate's coefficients integration
for i=1:length(T)
  % Direct reaction rate
    const=(8*q/(emass*pi))^(1/2)*1/(T(i))^(3/2);
    fun=const*eps.*exp(-eps/T(i)); 
    f=sigma.*fun;               
    dir_rate(i)=trapz(eps,f);              % Integral derived according to [Bosi, Phd Thesis, Eq. 2.5] [m^3/s]                 
end
end

