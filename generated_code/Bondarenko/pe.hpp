#ifndef _CML_bondarenko_model_2004_apex_pe_
#define _CML_bondarenko_model_2004_apex_pe_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: bondarenko_model_2004_apex
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translate: 1844, pycml: 1844)
//! on Wed Jul 20 18:35:40 2011
//! 
//! <autogenerated>



// 
// Settable parameters and readable variables
// 
void __equation__pe_(Solveode *classe)
{
    // Inputs:
    // Time units: millisecond
    double var_environment__time= classe->time_new;
    double var_membrane__i_stim;
    double var_membrane__i_CaL;
    double var_membrane__i_pCa;
    double var_membrane__i_NaCa;
    double var_membrane__i_Cab;
    double var_membrane__i_Na;
    double var_membrane__i_Nab;
    double var_membrane__i_NaK;
    double var_membrane__i_Kto_f;
    double var_membrane__i_Kto_s;
    double var_membrane__i_K1;
    double var_membrane__i_Ks;
    double var_membrane__i_Kur;
    double var_membrane__i_Kss;
    double var_membrane__i_ClCa;
    double var_membrane__i_Kr;
    double var_membrane__V = classe->V_old_; //var_membrane
    // Units: millivolt; Initial value: -82.4202
    double var_calcium_concentration__Cai = classe->Cai_old_; //var_calcium_concentration
    // Units: micromolar; Initial value: 0.115001
    double var_calcium_concentration__Cass = classe->Cass_old_; //var_calcium_concentration
    // Units: micromolar; Initial value: 0.115001
    double var_calcium_concentration__CaJSR = classe->CaJSR_old_; //var_calcium_concentration
    // Units: micromolar; Initial value: 1299.5
    double var_calcium_concentration__CaNSR = classe->CaNSR_old_; //var_calcium_concentration
    // Units: micromolar; Initial value: 1299.5
    double var_calcium_fluxes__P_RyR = classe->P_RyR_old_; //var_calcium_fluxes
    // Units: dimensionless; Initial value: 0
    double var_calcium_buffering__LTRPN_Ca = classe->LTRPN_Ca_old_; //var_calcium_buffering
    // Units: micromolar; Initial value: 11.2684
    double var_calcium_buffering__HTRPN_Ca = classe->HTRPN_Ca_old_; //var_calcium_buffering
    // Units: micromolar; Initial value: 125.29
    double var_ryanodine_receptors__P_O1 = classe->P_O1_old_; //var_ryanodine_receptors
    // Units: dimensionless; Initial value: 0.149102e-4
    double var_ryanodine_receptors__P_O2 = classe->P_O2_old_; //var_ryanodine_receptors
    // Units: dimensionless; Initial value: 0.951726e-10
    double var_ryanodine_receptors__P_C2 = classe->P_C2_old_; //var_ryanodine_receptors
    // Units: dimensionless; Initial value: 0.16774e-3
    double var_L_type_calcium_current__O = classe->O_old_; //var_L_type_calcium_current
    // Units: dimensionless; Initial value: 0.930308e-18
    double var_L_type_calcium_current__C2 = classe->C2_old_; //var_L_type_calcium_current
    // Units: dimensionless; Initial value: 0.124216e-3
    double var_L_type_calcium_current__C3 = classe->C3_old_; //var_L_type_calcium_current
    // Units: dimensionless; Initial value: 0.578679e-8
    double var_L_type_calcium_current__C4 = classe->C4_old_; //var_L_type_calcium_current
    // Units: dimensionless; Initial value: 0.119816e-12
    double var_L_type_calcium_current__I1 = classe->I1_old_; //var_L_type_calcium_current
    // Units: dimensionless; Initial value: 0.497923e-18
    double var_L_type_calcium_current__I2 = classe->I2_old_; //var_L_type_calcium_current
    // Units: dimensionless; Initial value: 0.345847e-13
    double var_L_type_calcium_current__I3 = classe->I3_old_; //var_L_type_calcium_current
    // Units: dimensionless; Initial value: 0.185106e-13
    double var_sodium_concentration__Nai = classe->Nai_old_; //var_sodium_concentration
    // Units: micromolar; Initial value: 14237.1
    double var_fast_sodium_current__O_Na = classe->O_Na_old_; //var_fast_sodium_current
    // Units: dimensionless; Initial value: 0.713483e-6
    double var_fast_sodium_current__C_Na1 = classe->C_Na1_old_; //var_fast_sodium_current
    // Units: dimensionless; Initial value: 0.279132e-3
    double var_fast_sodium_current__C_Na2 = classe->C_Na2_old_; //var_fast_sodium_current
    // Units: dimensionless; Initial value: 0.020752
    double var_fast_sodium_current__I1_Na = classe->I1_Na_old_; //var_fast_sodium_current
    // Units: dimensionless; Initial value: 0.673345e-6
    double var_fast_sodium_current__I2_Na = classe->I2_Na_old_; //var_fast_sodium_current
    // Units: dimensionless; Initial value: 0.155787e-8
    double var_fast_sodium_current__IF_Na = classe->IF_Na_old_; //var_fast_sodium_current
    // Units: dimensionless; Initial value: 0.153176e-3
    double var_fast_sodium_current__IC_Na2 = classe->IC_Na2_old_; //var_fast_sodium_current
    // Units: dimensionless; Initial value: 0.0113879
    double var_fast_sodium_current__IC_Na3 = classe->IC_Na3_old_; //var_fast_sodium_current
    // Units: dimensionless; Initial value: 0.34278
    double var_potassium_concentration__Ki = classe->Ki_old_; //var_potassium_concentration
    // Units: micromolar; Initial value: 143720
    double var_fast_transient_outward_potassium_current__ato_f = classe->ato_f_old_; //var_fast_transient_outward_potassium_current
    // Units: dimensionless; Initial value: 0.265563e-2
    double var_fast_transient_outward_potassium_current__ito_f = classe->ito_f_old_; //var_fast_transient_outward_potassium_current
    // Units: dimensionless; Initial value: 0.999977
    double var_slow_transient_outward_potassium_current__ato_s = classe->ato_s_old_; //var_slow_transient_outward_potassium_current
    // Units: dimensionless; Initial value: 0.417069e-3
    double var_slow_transient_outward_potassium_current__ito_s = classe->ito_s_old_; //var_slow_transient_outward_potassium_current
    // Units: dimensionless; Initial value: 0.998543
    double var_slow_delayed_rectifier_potassium_current__nKs = classe->nKs_old_; //var_slow_delayed_rectifier_potassium_current
    // Units: dimensionless; Initial value: 0.262753e-3
    double var_ultra_rapidly_activating_delayed_rectifier_potassium_current__aur = classe->aur_old_; //var_ultra_rapidly_activating_delayed_rectifier_potassium_current
    // Units: dimensionless; Initial value: 0.417069e-3
    double var_ultra_rapidly_activating_delayed_rectifier_potassium_current__iur = classe->iur_old_; //var_ultra_rapidly_activating_delayed_rectifier_potassium_current
    // Units: dimensionless; Initial value: 0.998543
    double var_non_inactivating_steady_state_potassium_current__aKss = classe->aKss_old_; //var_non_inactivating_steady_state_potassium_current
    // Units: dimensionless; Initial value: 0.417069e-3
    double var_non_inactivating_steady_state_potassium_current__iKss = classe->iKss_old_; //var_non_inactivating_steady_state_potassium_current
    // Units: dimensionless; Initial value: 1
    double var_rapid_delayed_rectifier_potassium_current__O_K = classe->O_K_old_; //var_rapid_delayed_rectifier_potassium_current
    // Units: dimensionless; Initial value: 0.175298e-3
    double var_rapid_delayed_rectifier_potassium_current__C_K1 = classe->C_K1_old_; //var_rapid_delayed_rectifier_potassium_current
    // Units: dimensionless; Initial value: 0.992513e-3
    double var_rapid_delayed_rectifier_potassium_current__C_K2 = classe->C_K2_old_; //var_rapid_delayed_rectifier_potassium_current
    // Units: dimensionless; Initial value: 0.641229e-3
    double var_rapid_delayed_rectifier_potassium_current__I_K = classe->I_K_old_; //var_rapid_delayed_rectifier_potassium_current
    // Units: dimensionless; Initial value: 0.319129e-4
    
    
    // Mathematics
    var_membrane__i_stim = ((var_environment__time >= 20.0) && (var_environment__time <= 100000.0) && (((var_environment__time - 20.0) - (floor((var_environment__time - 20.0) * 0.0139997200056) * 71.43)) <= 1.5)) ?  -80.0 : 0.0;
    const double var_L_type_calcium_current__i_CaL = 0.1729 * var_L_type_calcium_current__O * (var_membrane__V - 63.0);
    var_membrane__i_CaL = var_L_type_calcium_current__i_CaL;
    const double var_calcium_pump_current__i_pCa = (1.0 * pow(var_calcium_concentration__Cai, 2.0)) / (0.25 + pow(var_calcium_concentration__Cai, 2.0));
    var_membrane__i_pCa = var_calcium_pump_current__i_pCa;
    const double var_sodium_calcium_exchange_current__i_NaCa = (2.69705854643e-17 / (1.0 + (0.1 * exp(( -0.65 * var_membrane__V * 96.5) * 0.000403620964396)))) * ((exp((0.35 * var_membrane__V * 96.5) * 0.000403620964396) * pow(var_sodium_concentration__Nai, 3.0) * 1800.0) - (exp(( -0.65 * var_membrane__V * 96.5) * 0.000403620964396) * 2.744e+15 * var_calcium_concentration__Cai));
    var_membrane__i_NaCa = var_sodium_calcium_exchange_current__i_NaCa;
    const double var_calcium_background_current__i_Cab = 0.000367 * (var_membrane__V - (12.8371606218 * log(1800.0 / var_calcium_concentration__Cai)));
    var_membrane__i_Cab = var_calcium_background_current__i_Cab;
    const double var_fast_sodium_current__E_Na = 25.6743212435 * log(126540.0 / ((0.9 * var_sodium_concentration__Nai) + (0.1 * var_potassium_concentration__Ki)));
    const double var_fast_sodium_current__i_Na = 13.0 * var_fast_sodium_current__O_Na * (var_membrane__V - var_fast_sodium_current__E_Na);
    var_membrane__i_Na = var_fast_sodium_current__i_Na;
    const double var_sodium_background_current__i_Nab = 0.0026 * (var_membrane__V - var_fast_sodium_current__E_Na);
    var_membrane__i_Nab = var_sodium_background_current__i_Nab;
    const double var_sodium_potassium_pump_current__i_NaK = (((0.88 * (1.0 / (1.0 + (0.1245 * exp(( -0.1 * var_membrane__V * 96.5) * 0.000403620964396)) + (0.0365 * 1.00091030495 * exp(((-var_membrane__V) * 96.5) * 0.000403620964396)))) * 1.0) / (1.0 + pow(21000.0 / var_sodium_concentration__Nai, 1.5))) * 5400.0) * 0.000144927536232;
    var_membrane__i_NaK = var_sodium_potassium_pump_current__i_NaK;
    const double var_fast_transient_outward_potassium_current__E_K = 25.6743212435 * log(5400.0 / var_potassium_concentration__Ki);
    const double var_fast_transient_outward_potassium_current__i_Kto_f = 0.4067 * pow(var_fast_transient_outward_potassium_current__ato_f, 3.0) * var_fast_transient_outward_potassium_current__ito_f * (var_membrane__V - var_fast_transient_outward_potassium_current__E_K);
    var_membrane__i_Kto_f = var_fast_transient_outward_potassium_current__i_Kto_f;
    const double var_slow_transient_outward_potassium_current__i_Kto_s = 0.0 * var_slow_transient_outward_potassium_current__ato_s * var_slow_transient_outward_potassium_current__ito_s * (var_membrane__V - var_fast_transient_outward_potassium_current__E_K);
    var_membrane__i_Kto_s = var_slow_transient_outward_potassium_current__i_Kto_s;
    const double var_time_independent_potassium_current__i_K1 = (0.282802139037 * (var_membrane__V - var_fast_transient_outward_potassium_current__E_K)) / (1.0 + exp(0.0896 * (var_membrane__V - var_fast_transient_outward_potassium_current__E_K)));
    var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
    const double var_slow_delayed_rectifier_potassium_current__i_Ks = 0.00575 * pow(var_slow_delayed_rectifier_potassium_current__nKs, 2.0) * (var_membrane__V - var_fast_transient_outward_potassium_current__E_K);
    var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
    const double var_ultra_rapidly_activating_delayed_rectifier_potassium_current__i_Kur = 0.16 * var_ultra_rapidly_activating_delayed_rectifier_potassium_current__aur * var_ultra_rapidly_activating_delayed_rectifier_potassium_current__iur * (var_membrane__V - var_fast_transient_outward_potassium_current__E_K);
    var_membrane__i_Kur = var_ultra_rapidly_activating_delayed_rectifier_potassium_current__i_Kur;
    const double var_non_inactivating_steady_state_potassium_current__i_Kss = 0.05 * var_non_inactivating_steady_state_potassium_current__aKss * var_non_inactivating_steady_state_potassium_current__iKss * (var_membrane__V - var_fast_transient_outward_potassium_current__E_K);
    var_membrane__i_Kss = var_non_inactivating_steady_state_potassium_current__i_Kss;
    var_membrane__i_ClCa = ((10.0 * (0.2 / (1.0 + exp((-(var_membrane__V - 46.7)) * 0.128205128205))) * var_calcium_concentration__Cai) / (var_calcium_concentration__Cai + 10.0)) * (var_membrane__V -  -40.0);
    const double var_rapid_delayed_rectifier_potassium_current__i_Kr = 0.078 * var_rapid_delayed_rectifier_potassium_current__O_K * (var_membrane__V - (25.6743212435 * log(8092.0 / ((0.98 * var_potassium_concentration__Ki) + (0.02 * var_sodium_concentration__Nai)))));
    var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
    const double var_calcium_fluxes__J_leak = 1.74e-05 * (var_calcium_concentration__CaNSR - var_calcium_concentration__Cai);
    const double var_calcium_fluxes__J_rel = 4.5 * (var_ryanodine_receptors__P_O1 + var_ryanodine_receptors__P_O2) * (var_calcium_concentration__CaJSR - var_calcium_concentration__Cass) * var_calcium_fluxes__P_RyR;
    const double var_calcium_fluxes__J_up = (0.45 * pow(var_calcium_concentration__Cai, 2.0)) / (0.25 + pow(var_calcium_concentration__Cai, 2.0));
    const double var_calcium_fluxes__J_tr = (var_calcium_concentration__CaNSR - var_calcium_concentration__CaJSR) * 0.05;
    const double var_calcium_fluxes__J_xfer = (var_calcium_concentration__Cass - var_calcium_concentration__Cai) * 0.125;
    const double var_L_type_calcium_current__alpha = (0.4 * exp((var_membrane__V + 12.0) * 0.1) * ((1.0 + (0.7 * exp((-pow(var_membrane__V + 40.0, 2.0)) * 0.1))) - (0.75 * exp((-pow(var_membrane__V + 20.0, 2.0)) * 0.0025)))) / (1.0 + (0.12 * exp((var_membrane__V + 12.0) * 0.1)));
    const double var_L_type_calcium_current__beta = 0.05 * exp((-(var_membrane__V + 12.0)) * 0.0769230769231);
    const double var_L_type_calcium_current__gamma = (0.23324 * var_calcium_concentration__Cass) / (20.0 + var_calcium_concentration__Cass);
    const double var_L_type_calcium_current__Kpcf = 13.0 * (1.0 - exp((-pow(var_membrane__V + 14.5, 2.0)) * 0.01));
    const double var_fast_sodium_current__C_Na3 = 1.0 - (var_fast_sodium_current__O_Na + var_fast_sodium_current__C_Na1 + var_fast_sodium_current__C_Na2 + var_fast_sodium_current__IF_Na + var_fast_sodium_current__I1_Na + var_fast_sodium_current__I2_Na + var_fast_sodium_current__IC_Na2 + var_fast_sodium_current__IC_Na3);
    const double var_fast_sodium_current__alpha_Na11 = 3.802 / ((0.1027 * exp((-(var_membrane__V + 2.5)) * 0.0588235294118)) + (0.2 * exp((-(var_membrane__V + 2.5)) * 0.00666666666667)));
    const double var_fast_sodium_current__beta_Na11 = 0.1917 * exp((-(var_membrane__V + 2.5)) * 0.0492610837438);
    const double var_fast_sodium_current__alpha_Na12 = 3.802 / ((0.1027 * exp((-(var_membrane__V + 2.5)) * 0.0666666666667)) + (0.23 * exp((-(var_membrane__V + 2.5)) * 0.00666666666667)));
    const double var_fast_sodium_current__beta_Na12 = 0.2 * exp((-(var_membrane__V - 2.5)) * 0.0492610837438);
    const double var_fast_sodium_current__alpha_Na13 = 3.802 / ((0.1027 * exp((-(var_membrane__V + 2.5)) * 0.0833333333333)) + (0.25 * exp((-(var_membrane__V + 2.5)) * 0.00666666666667)));
    const double var_fast_sodium_current__beta_Na13 = 0.22 * exp((-(var_membrane__V - 7.5)) * 0.0492610837438);
    const double var_fast_sodium_current__alpha_Na3 = 7e-07 * exp((-(var_membrane__V + 7.0)) * 0.12987012987);
    const double var_fast_sodium_current__beta_Na3 = 0.0084 + (2e-05 * (var_membrane__V + 7.0));
    const double var_fast_sodium_current__alpha_Na2 = 1.0 / ((0.188495 * exp((-(var_membrane__V + 7.0)) * 0.0602409638554)) + 0.393956);
    const double var_fast_sodium_current__beta_Na2 = (var_fast_sodium_current__alpha_Na13 * var_fast_sodium_current__alpha_Na2 * var_fast_sodium_current__alpha_Na3) / (var_fast_sodium_current__beta_Na13 * var_fast_sodium_current__beta_Na3);
    const double var_fast_sodium_current__alpha_Na4 = var_fast_sodium_current__alpha_Na2 * 0.001;
    const double var_fast_sodium_current__beta_Na4 = var_fast_sodium_current__alpha_Na3;
    const double var_fast_sodium_current__alpha_Na5 = var_fast_sodium_current__alpha_Na2 * 1.05263157895e-05;
    const double var_fast_sodium_current__beta_Na5 = var_fast_sodium_current__alpha_Na3 * 0.02;
    const double var_slow_transient_outward_potassium_current__ass = 1.0 / (1.0 + exp((-(var_membrane__V + 22.5)) * 0.12987012987));
    const double var_slow_transient_outward_potassium_current__iss = 1.0 / (1.0 + exp((var_membrane__V + 45.2) * 0.175438596491));
    const double var_rapid_delayed_rectifier_potassium_current__alpha_a1 = 0.013733 * exp(0.038198 * var_membrane__V);
    const double var_rapid_delayed_rectifier_potassium_current__beta_a1 = 6.89e-05 * exp( -0.04178 * var_membrane__V);
    const double var_rapid_delayed_rectifier_potassium_current__alpha_i = 0.090821 * exp(0.023391 * (var_membrane__V + 5.0));
    const double var_rapid_delayed_rectifier_potassium_current__beta_i = 0.006497 * exp( -0.03268 * (var_membrane__V + 5.0));
    
    const double d_dt_membrane__V = -(var_membrane__i_CaL + var_membrane__i_pCa + var_membrane__i_NaCa + var_membrane__i_Cab + var_membrane__i_Na + var_membrane__i_Nab + var_membrane__i_NaK + var_membrane__i_Kto_f + var_membrane__i_Kto_s + var_membrane__i_K1 + var_membrane__i_Ks + var_membrane__i_Kur + var_membrane__i_Kss + var_membrane__i_Kr + var_membrane__i_ClCa + var_membrane__i_stim);
    const double d_dt_calcium_concentration__Cai = pow(1.0 + (11.9 / pow(0.238 + var_calcium_concentration__Cai, 2.0)),  -1.0) * ((var_calcium_fluxes__J_leak + var_calcium_fluxes__J_xfer) - (var_calcium_fluxes__J_up + (((0.00237 * var_calcium_concentration__Cai * (140.0 - var_calcium_buffering__HTRPN_Ca)) + (0.0327 * var_calcium_concentration__Cai * (70.0 - var_calcium_buffering__LTRPN_Ca))) - ((3.2e-05 * var_calcium_buffering__HTRPN_Ca) + (0.0196 * var_calcium_buffering__LTRPN_Ca))) + ((((var_calcium_background_current__i_Cab + var_calcium_pump_current__i_pCa) - (2.0 * var_sodium_calcium_exchange_current__i_NaCa)) * 0.0001534 * 1.0) * 200.516530583)));
    const double d_dt_calcium_concentration__Cass = pow(1.0 + (11.9 / pow(0.238 + var_calcium_concentration__Cass, 2.0)),  -1.0) * (((var_calcium_fluxes__J_rel * 1.2e-07) * 673400673.401) - (((var_calcium_fluxes__J_xfer * 2.584e-05) * 673400673.401) + ((var_L_type_calcium_current__i_CaL * 0.0001534 * 1.0) * 3489122.66011)));
    const double d_dt_calcium_concentration__CaJSR = pow(1.0 + (12000000.0 / pow(800.0 + var_calcium_concentration__CaJSR, 2.0)),  -1.0) * (var_calcium_fluxes__J_tr - var_calcium_fluxes__J_rel);
    const double d_dt_calcium_concentration__CaNSR = (((var_calcium_fluxes__J_up - var_calcium_fluxes__J_leak) * 2.584e-05) * 476644.42326) - ((var_calcium_fluxes__J_tr * 1.2e-07) * 476644.42326);
    const double d_dt_calcium_fluxes__P_RyR = ( -0.04 * var_calcium_fluxes__P_RyR) - (((0.1 * var_L_type_calcium_current__i_CaL) * 0.142857142857) * exp((-pow(var_membrane__V - 5.0, 2.0)) * 0.00154320987654));
    const double d_dt_calcium_buffering__LTRPN_Ca = (0.0327 * var_calcium_concentration__Cai * (70.0 - var_calcium_buffering__LTRPN_Ca)) - (0.0196 * var_calcium_buffering__LTRPN_Ca);
    const double d_dt_calcium_buffering__HTRPN_Ca = (0.00237 * var_calcium_concentration__Cai * (140.0 - var_calcium_buffering__HTRPN_Ca)) - (3.2e-05 * var_calcium_buffering__HTRPN_Ca);
    const double d_dt_ryanodine_receptors__P_O1 = ((0.006075 * pow(var_calcium_concentration__Cass, 4.0) * (1.0 - (var_ryanodine_receptors__P_C2 + var_ryanodine_receptors__P_O1 + var_ryanodine_receptors__P_O2))) + (0.965 * var_ryanodine_receptors__P_O2) + (0.0008 * var_ryanodine_receptors__P_C2)) - ((0.07125 * var_ryanodine_receptors__P_O1) + (0.00405 * pow(var_calcium_concentration__Cass, 3.0) * var_ryanodine_receptors__P_O1) + (0.009 * var_ryanodine_receptors__P_O1));
    const double d_dt_ryanodine_receptors__P_O2 = (0.00405 * pow(var_calcium_concentration__Cass, 3.0) * var_ryanodine_receptors__P_O1) - (0.965 * var_ryanodine_receptors__P_O2);
    const double d_dt_ryanodine_receptors__P_C2 = (0.009 * var_ryanodine_receptors__P_O1) - (0.0008 * var_ryanodine_receptors__P_C2);
    const double d_dt_L_type_calcium_current__O = ((var_L_type_calcium_current__alpha * var_L_type_calcium_current__C4) + (0.0005 * var_L_type_calcium_current__I1) + (0.001 * ((var_L_type_calcium_current__alpha * var_L_type_calcium_current__I2) - (var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__O)))) - ((4.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__O) + (var_L_type_calcium_current__gamma * var_L_type_calcium_current__O));
    const double d_dt_L_type_calcium_current__C2 = ((4.0 * var_L_type_calcium_current__alpha * (1.0 - (var_L_type_calcium_current__O + var_L_type_calcium_current__C2 + var_L_type_calcium_current__C2 + var_L_type_calcium_current__C3 + var_L_type_calcium_current__C4 + var_L_type_calcium_current__I1 + var_L_type_calcium_current__I2 + var_L_type_calcium_current__I3))) + (2.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__C3)) - ((var_L_type_calcium_current__beta * var_L_type_calcium_current__C2) + (3.0 * var_L_type_calcium_current__alpha * var_L_type_calcium_current__C2));
    const double d_dt_L_type_calcium_current__C3 = ((3.0 * var_L_type_calcium_current__alpha * var_L_type_calcium_current__C2) + (3.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__C4)) - ((2.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__C3) + (2.0 * var_L_type_calcium_current__alpha * var_L_type_calcium_current__C3));
    const double d_dt_L_type_calcium_current__C4 = ((2.0 * var_L_type_calcium_current__alpha * var_L_type_calcium_current__C3) + (4.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__O) + (0.01 * ((4.0 * 0.0005 * var_L_type_calcium_current__beta * var_L_type_calcium_current__I1) - (var_L_type_calcium_current__alpha * var_L_type_calcium_current__gamma * var_L_type_calcium_current__C4))) + (0.002 * ((4.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__I2) - (var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__C4))) + (4.0 * var_L_type_calcium_current__beta * 0.0005 * var_L_type_calcium_current__I3)) - ((3.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__C4) + (var_L_type_calcium_current__alpha * var_L_type_calcium_current__C4) + (1.0 * var_L_type_calcium_current__gamma * var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__C4));
    const double d_dt_L_type_calcium_current__I1 = ((var_L_type_calcium_current__gamma * var_L_type_calcium_current__O) + (0.001 * ((var_L_type_calcium_current__alpha * var_L_type_calcium_current__I3) - (var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__I1))) + (0.01 * ((var_L_type_calcium_current__alpha * var_L_type_calcium_current__gamma * var_L_type_calcium_current__C4) - (4.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__I1)))) - (0.0005 * var_L_type_calcium_current__I1);
    const double d_dt_L_type_calcium_current__I2 = ((0.001 * ((var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__O) - (var_L_type_calcium_current__alpha * var_L_type_calcium_current__I2))) + (0.0005 * var_L_type_calcium_current__I3) + (0.002 * ((var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__C4) - (4.0 * var_L_type_calcium_current__beta * var_L_type_calcium_current__I2)))) - (var_L_type_calcium_current__gamma * var_L_type_calcium_current__I2);
    const double d_dt_L_type_calcium_current__I3 = ((0.001 * ((var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__I1) - (var_L_type_calcium_current__alpha * var_L_type_calcium_current__I3))) + (var_L_type_calcium_current__gamma * var_L_type_calcium_current__I2) + (1.0 * var_L_type_calcium_current__gamma * var_L_type_calcium_current__Kpcf * var_L_type_calcium_current__C4)) - ((4.0 * var_L_type_calcium_current__beta * 0.0005 * var_L_type_calcium_current__I3) + (0.0005 * var_L_type_calcium_current__I3));
    const double d_dt_sodium_concentration__Nai = ((-(var_fast_sodium_current__i_Na + var_sodium_background_current__i_Nab + (3.0 * var_sodium_potassium_pump_current__i_NaK) + (3.0 * var_sodium_calcium_exchange_current__i_NaCa))) * 0.0001534 * 1.0) * 401.033061166;
    const double d_dt_fast_sodium_current__C_Na2 = ((var_fast_sodium_current__alpha_Na11 * var_fast_sodium_current__C_Na3) + (var_fast_sodium_current__beta_Na12 * var_fast_sodium_current__C_Na1) + (var_fast_sodium_current__alpha_Na3 * var_fast_sodium_current__IC_Na2)) - ((var_fast_sodium_current__beta_Na11 * var_fast_sodium_current__C_Na2) + (var_fast_sodium_current__alpha_Na12 * var_fast_sodium_current__C_Na2) + (var_fast_sodium_current__beta_Na3 * var_fast_sodium_current__C_Na2));
    const double d_dt_fast_sodium_current__C_Na1 = ((var_fast_sodium_current__alpha_Na12 * var_fast_sodium_current__C_Na2) + (var_fast_sodium_current__beta_Na13 * var_fast_sodium_current__O_Na) + (var_fast_sodium_current__alpha_Na3 * var_fast_sodium_current__IF_Na)) - ((var_fast_sodium_current__beta_Na12 * var_fast_sodium_current__C_Na1) + (var_fast_sodium_current__alpha_Na13 * var_fast_sodium_current__C_Na1) + (var_fast_sodium_current__beta_Na3 * var_fast_sodium_current__C_Na1));
    const double d_dt_fast_sodium_current__O_Na = ((var_fast_sodium_current__alpha_Na13 * var_fast_sodium_current__C_Na1) + (var_fast_sodium_current__beta_Na2 * var_fast_sodium_current__IF_Na)) - ((var_fast_sodium_current__beta_Na13 * var_fast_sodium_current__O_Na) + (var_fast_sodium_current__alpha_Na2 * var_fast_sodium_current__O_Na));
    const double d_dt_fast_sodium_current__IF_Na = ((var_fast_sodium_current__alpha_Na2 * var_fast_sodium_current__O_Na) + (var_fast_sodium_current__beta_Na3 * var_fast_sodium_current__C_Na1) + (var_fast_sodium_current__beta_Na4 * var_fast_sodium_current__I1_Na) + (var_fast_sodium_current__alpha_Na12 * var_fast_sodium_current__IC_Na2)) - ((var_fast_sodium_current__beta_Na2 * var_fast_sodium_current__IF_Na) + (var_fast_sodium_current__alpha_Na3 * var_fast_sodium_current__IF_Na) + (var_fast_sodium_current__alpha_Na4 * var_fast_sodium_current__IF_Na) + (var_fast_sodium_current__beta_Na12 * var_fast_sodium_current__IF_Na));
    const double d_dt_fast_sodium_current__I1_Na = ((var_fast_sodium_current__alpha_Na4 * var_fast_sodium_current__IF_Na) + (var_fast_sodium_current__beta_Na5 * var_fast_sodium_current__I2_Na)) - ((var_fast_sodium_current__beta_Na4 * var_fast_sodium_current__I1_Na) + (var_fast_sodium_current__alpha_Na5 * var_fast_sodium_current__I1_Na));
    const double d_dt_fast_sodium_current__I2_Na = (var_fast_sodium_current__alpha_Na5 * var_fast_sodium_current__I1_Na) - (var_fast_sodium_current__beta_Na5 * var_fast_sodium_current__I2_Na);
    const double d_dt_fast_sodium_current__IC_Na2 = ((var_fast_sodium_current__alpha_Na11 * var_fast_sodium_current__IC_Na3) + (var_fast_sodium_current__beta_Na12 * var_fast_sodium_current__IF_Na) + (var_fast_sodium_current__beta_Na3 * var_fast_sodium_current__IC_Na2)) - ((var_fast_sodium_current__beta_Na11 * var_fast_sodium_current__IC_Na2) + (var_fast_sodium_current__alpha_Na12 * var_fast_sodium_current__IC_Na2) + (var_fast_sodium_current__alpha_Na3 * var_fast_sodium_current__IC_Na2));
    const double d_dt_fast_sodium_current__IC_Na3 = ((var_fast_sodium_current__beta_Na11 * var_fast_sodium_current__IC_Na2) + (var_fast_sodium_current__beta_Na3 * var_fast_sodium_current__C_Na3)) - ((var_fast_sodium_current__alpha_Na11 * var_fast_sodium_current__IC_Na3) + (var_fast_sodium_current__alpha_Na3 * var_fast_sodium_current__IC_Na3));
    const double d_dt_potassium_concentration__Ki = ((-((var_fast_transient_outward_potassium_current__i_Kto_f + var_slow_transient_outward_potassium_current__i_Kto_s + var_time_independent_potassium_current__i_K1 + var_slow_delayed_rectifier_potassium_current__i_Ks + var_non_inactivating_steady_state_potassium_current__i_Kss + var_ultra_rapidly_activating_delayed_rectifier_potassium_current__i_Kur + var_rapid_delayed_rectifier_potassium_current__i_Kr) - (2.0 * var_sodium_potassium_pump_current__i_NaK))) * 0.0001534 * 1.0) * 401.033061166;
    const double d_dt_fast_transient_outward_potassium_current__ato_f = ((0.18064 * exp(0.03577 * (var_membrane__V + 30.0))) * (1.0 - var_fast_transient_outward_potassium_current__ato_f)) - ((0.3956 * exp( -0.06237 * (var_membrane__V + 30.0))) * var_fast_transient_outward_potassium_current__ato_f);
    const double d_dt_fast_transient_outward_potassium_current__ito_f = (((0.000152 * exp((-(var_membrane__V + 13.5)) * 0.142857142857)) / ((0.0067083 * exp((-(var_membrane__V + 33.5)) * 0.142857142857)) + 1.0)) * (1.0 - var_fast_transient_outward_potassium_current__ito_f)) - (((0.00095 * exp((var_membrane__V + 33.5) * 0.142857142857)) / ((0.051335 * exp((var_membrane__V + 33.5) * 0.142857142857)) + 1.0)) * var_fast_transient_outward_potassium_current__ito_f);
    const double d_dt_slow_transient_outward_potassium_current__ato_s = (var_slow_transient_outward_potassium_current__ass - var_slow_transient_outward_potassium_current__ato_s) / ((0.493 * exp( -0.0629 * var_membrane__V)) + 2.058);
    const double d_dt_slow_transient_outward_potassium_current__ito_s = (var_slow_transient_outward_potassium_current__iss - var_slow_transient_outward_potassium_current__ito_s) / (270.0 + (1050.0 / (1.0 + exp((var_membrane__V + 45.2) * 0.175438596491))));
    const double d_dt_slow_delayed_rectifier_potassium_current__nKs = (((4.81333e-06 * (var_membrane__V + 26.5)) / (1.0 - exp( -0.128 * (var_membrane__V + 26.5)))) * (1.0 - var_slow_delayed_rectifier_potassium_current__nKs)) - ((9.53333e-05 * exp( -0.038 * (var_membrane__V + 26.5))) * var_slow_delayed_rectifier_potassium_current__nKs);
    const double d_dt_ultra_rapidly_activating_delayed_rectifier_potassium_current__aur = (var_slow_transient_outward_potassium_current__ass - var_ultra_rapidly_activating_delayed_rectifier_potassium_current__aur) / ((0.493 * exp( -0.0629 * var_membrane__V)) + 2.058);
    const double d_dt_ultra_rapidly_activating_delayed_rectifier_potassium_current__iur = (var_slow_transient_outward_potassium_current__iss - var_ultra_rapidly_activating_delayed_rectifier_potassium_current__iur) / (1200.0 - (170.0 / (1.0 + exp((var_membrane__V + 45.2) * 0.175438596491))));
    const double d_dt_non_inactivating_steady_state_potassium_current__aKss = (var_slow_transient_outward_potassium_current__ass - var_non_inactivating_steady_state_potassium_current__aKss) / ((39.3 * exp( -0.0862 * var_membrane__V)) + 13.17);
    const double d_dt_non_inactivating_steady_state_potassium_current__iKss = 0.0;
    const double d_dt_rapid_delayed_rectifier_potassium_current__C_K2 = ((0.023761 * var_rapid_delayed_rectifier_potassium_current__C_K1) + (var_rapid_delayed_rectifier_potassium_current__beta_a1 * var_rapid_delayed_rectifier_potassium_current__O_K)) - ((0.036778 * var_rapid_delayed_rectifier_potassium_current__C_K2) + (var_rapid_delayed_rectifier_potassium_current__alpha_a1 * var_rapid_delayed_rectifier_potassium_current__C_K2));
    const double d_dt_rapid_delayed_rectifier_potassium_current__C_K1 = (((0.022348 * exp(0.01176 * var_membrane__V)) * (1.0 - (var_rapid_delayed_rectifier_potassium_current__C_K1 + var_rapid_delayed_rectifier_potassium_current__C_K2 + var_rapid_delayed_rectifier_potassium_current__O_K + var_rapid_delayed_rectifier_potassium_current__I_K))) + (0.036778 * var_rapid_delayed_rectifier_potassium_current__C_K2)) - (((0.047002 * exp( -0.0631 * var_membrane__V)) * var_rapid_delayed_rectifier_potassium_current__C_K1) + (0.023761 * var_rapid_delayed_rectifier_potassium_current__C_K1));
    const double d_dt_rapid_delayed_rectifier_potassium_current__O_K = ((var_rapid_delayed_rectifier_potassium_current__alpha_a1 * var_rapid_delayed_rectifier_potassium_current__C_K2) + (var_rapid_delayed_rectifier_potassium_current__beta_i * var_rapid_delayed_rectifier_potassium_current__I_K)) - ((var_rapid_delayed_rectifier_potassium_current__beta_a1 * var_rapid_delayed_rectifier_potassium_current__O_K) + (var_rapid_delayed_rectifier_potassium_current__alpha_i * var_rapid_delayed_rectifier_potassium_current__O_K));
    const double d_dt_rapid_delayed_rectifier_potassium_current__I_K = (var_rapid_delayed_rectifier_potassium_current__alpha_i * var_rapid_delayed_rectifier_potassium_current__O_K) - (var_rapid_delayed_rectifier_potassium_current__beta_i * var_rapid_delayed_rectifier_potassium_current__I_K);
    
    classe->V_lado_direito_ = d_dt_membrane__V;
    classe->Cai_lado_direito_ = d_dt_calcium_concentration__Cai;
    classe->Cass_lado_direito_ = d_dt_calcium_concentration__Cass;
    classe->CaJSR_lado_direito_ = d_dt_calcium_concentration__CaJSR;
    classe->CaNSR_lado_direito_ = d_dt_calcium_concentration__CaNSR;
    classe->P_RyR_lado_direito_ = d_dt_calcium_fluxes__P_RyR;
    classe->LTRPN_Ca_lado_direito_ = d_dt_calcium_buffering__LTRPN_Ca;
    classe->HTRPN_Ca_lado_direito_ = d_dt_calcium_buffering__HTRPN_Ca;
    classe->P_O1_lado_direito_ = d_dt_ryanodine_receptors__P_O1;
    classe->P_O2_lado_direito_ = d_dt_ryanodine_receptors__P_O2;
    classe->P_C2_lado_direito_ = d_dt_ryanodine_receptors__P_C2;
    classe->O_lado_direito_ = d_dt_L_type_calcium_current__O;
    classe->C2_lado_direito_ = d_dt_L_type_calcium_current__C2;
    classe->C3_lado_direito_ = d_dt_L_type_calcium_current__C3;
    classe->C4_lado_direito_ = d_dt_L_type_calcium_current__C4;
    classe->I1_lado_direito_ = d_dt_L_type_calcium_current__I1;
    classe->I2_lado_direito_ = d_dt_L_type_calcium_current__I2;
    classe->I3_lado_direito_ = d_dt_L_type_calcium_current__I3;
    classe->Nai_lado_direito_ = d_dt_sodium_concentration__Nai;
    classe->O_Na_lado_direito_ = d_dt_fast_sodium_current__O_Na;
    classe->C_Na1_lado_direito_ = d_dt_fast_sodium_current__C_Na1;
    classe->C_Na2_lado_direito_ = d_dt_fast_sodium_current__C_Na2;
    classe->I1_Na_lado_direito_ = d_dt_fast_sodium_current__I1_Na;
    classe->I2_Na_lado_direito_ = d_dt_fast_sodium_current__I2_Na;
    classe->IF_Na_lado_direito_ = d_dt_fast_sodium_current__IF_Na;
    classe->IC_Na2_lado_direito_ = d_dt_fast_sodium_current__IC_Na2;
    classe->IC_Na3_lado_direito_ = d_dt_fast_sodium_current__IC_Na3;
    classe->Ki_lado_direito_ = d_dt_potassium_concentration__Ki;
    classe->ato_f_lado_direito_ = d_dt_fast_transient_outward_potassium_current__ato_f;
    classe->ito_f_lado_direito_ = d_dt_fast_transient_outward_potassium_current__ito_f;
    classe->ato_s_lado_direito_ = d_dt_slow_transient_outward_potassium_current__ato_s;
    classe->ito_s_lado_direito_ = d_dt_slow_transient_outward_potassium_current__ito_s;
    classe->nKs_lado_direito_ = d_dt_slow_delayed_rectifier_potassium_current__nKs;
    classe->aur_lado_direito_ = d_dt_ultra_rapidly_activating_delayed_rectifier_potassium_current__aur;
    classe->iur_lado_direito_ = d_dt_ultra_rapidly_activating_delayed_rectifier_potassium_current__iur;
    classe->aKss_lado_direito_ = d_dt_non_inactivating_steady_state_potassium_current__aKss;
    classe->iKss_lado_direito_ = d_dt_non_inactivating_steady_state_potassium_current__iKss;
    classe->O_K_lado_direito_ = d_dt_rapid_delayed_rectifier_potassium_current__O_K;
    classe->C_K1_lado_direito_ = d_dt_rapid_delayed_rectifier_potassium_current__C_K1;
    classe->C_K2_lado_direito_ = d_dt_rapid_delayed_rectifier_potassium_current__C_K2;
    classe->I_K_lado_direito_ = d_dt_rapid_delayed_rectifier_potassium_current__I_K;
};


#endif
