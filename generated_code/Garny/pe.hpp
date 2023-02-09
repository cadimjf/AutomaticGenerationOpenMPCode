#ifndef _CML_zhang_SAN_model_2000_all_pe_
#define _CML_zhang_SAN_model_2000_all_pe_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: zhang_SAN_model_2000_all
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translate: 1844, pycml: 1844)
//! on Sun Jun 26 21:57:55 2011
//! 
//! <autogenerated>



// 
// Settable parameters and readable variables
// 
void __equation__pe_(Solveode *classe)
{
    // Inputs:
    // Time units: second
    double var_environment__time= classe->time_new;
    double var_membrane__i_Na;
    double var_membrane__i_Ca_L;
    double var_membrane__i_Ca_T;
    double var_membrane__i_to;
    double var_membrane__i_sus;
    double var_membrane__i_K_r;
    double var_membrane__i_K_s;
    double var_membrane__i_f_Na;
    double var_membrane__i_f_K;
    double var_membrane__i_b_Na;
    double var_membrane__i_b_Ca;
    double var_membrane__i_b_K;
    double var_membrane__i_NaCa;
    double var_membrane__i_p;
    double var_membrane__i_Ca_p;
    double var_membrane__V = classe->V_old_; //var_membrane
    // Units: millivolt; Initial value: -39.013558536
    double var_sodium_current_m_gate__m = classe->m_old_; //var_sodium_current_m_gate
    // Units: dimensionless; Initial value: 0.092361701692
    double var_sodium_current_h_gate__h1 = classe->h1_old_; //var_sodium_current_h_gate
    // Units: dimensionless; Initial value: 0.015905380261
    double var_sodium_current_h_gate__h2 = classe->h2_old_; //var_sodium_current_h_gate
    // Units: dimensionless; Initial value: 0.01445216109
    double var_L_type_Ca_channel_d_gate__d_L = classe->d_L_old_; //var_L_type_Ca_channel_d_gate
    // Units: dimensionless; Initial value: 0.04804900895
    double var_L_type_Ca_channel_f_gate__f_L = classe->f_L_old_; //var_L_type_Ca_channel_f_gate
    // Units: dimensionless; Initial value: 0.48779845203
    double var_T_type_Ca_channel_d_gate__d_T = classe->d_T_old_; //var_T_type_Ca_channel_d_gate
    // Units: dimensionless; Initial value: 0.42074047435
    double var_T_type_Ca_channel_f_gate__f_T = classe->f_T_old_; //var_T_type_Ca_channel_f_gate
    // Units: dimensionless; Initial value: 0.038968420558
    double var_four_AP_sensitive_currents_q_gate__q = classe->q_old_; //var_four_AP_sensitive_currents_q_gate
    // Units: dimensionless; Initial value: 0.29760539675
    double var_four_AP_sensitive_currents_r_gate__r = classe->r_old_; //var_four_AP_sensitive_currents_r_gate
    // Units: dimensionless; Initial value: 0.064402950262
    double var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af = classe->P_af_old_; //var_rapid_delayed_rectifying_potassium_current_P_af_gate
    // Units: dimensionless; Initial value: 0.13034201158
    double var_rapid_delayed_rectifying_potassium_current_P_as_gate__P_as = classe->P_as_old_; //var_rapid_delayed_rectifying_potassium_current_P_as_gate
    // Units: dimensionless; Initial value: 0.46960956028
    double var_rapid_delayed_rectifying_potassium_current_P_i_gate__P_i = classe->P_i_old_; //var_rapid_delayed_rectifying_potassium_current_P_i_gate
    // Units: dimensionless; Initial value: 0.87993375273
    double var_slow_delayed_rectifying_potassium_current_xs_gate__xs = classe->xs_old_; //var_slow_delayed_rectifying_potassium_current_xs_gate
    // Units: dimensionless; Initial value: 0.082293827208
    double var_hyperpolarisation_activated_current_y_gate__y = classe->y_old_; //var_hyperpolarisation_activated_current_y_gate
    // Units: dimensionless; Initial value: 0.03889291759
    
    
    // Mathematics
    const double var_sodium_current_h_gate__F_Na = ((0.09518 * exp( -0.06306 * (var_membrane__V + 34.4))) / (1.0 + (1.662 * exp( -0.2251 * (var_membrane__V + 63.7))))) + 0.08693;
    var_membrane__i_Na = ((((0.0 * pow(var_sodium_current_m_gate__m, 3.0) * (((1.0 - var_sodium_current_h_gate__F_Na) * var_sodium_current_h_gate__h1) + (var_sodium_current_h_gate__F_Na * var_sodium_current_h_gate__h2)) * 140.0 * 9378954025.0) * 3.87996927064e-07) * (exp(((var_membrane__V - 76.1718707053) * 96845.0) * 3.87996927064e-07) - 1.0)) / (exp((var_membrane__V * 96845.0) * 3.87996927064e-07) - 1.0)) * var_membrane__V;
    var_membrane__i_Ca_L = 0.0057938 * ((var_L_type_Ca_channel_f_gate__f_L * var_L_type_Ca_channel_d_gate__d_L) + (0.006 / (1.0 + exp((-(var_membrane__V + 14.1)) * 0.166666666667)))) * (var_membrane__V - 46.4);
    var_membrane__i_Ca_T = 0.00427806 * var_T_type_Ca_channel_d_gate__d_T * var_T_type_Ca_channel_f_gate__f_T * (var_membrane__V - 45.0);
    var_membrane__i_to = 0.004905 * var_four_AP_sensitive_currents_q_gate__q * var_four_AP_sensitive_currents_r_gate__r * (var_membrane__V -  -86.6319293974);
    var_membrane__i_sus = 6.645504e-05 * var_four_AP_sensitive_currents_r_gate__r * (var_membrane__V -  -86.6319293974);
    var_membrane__i_K_r = 0.00079704 * ((0.6 * var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af) + (0.4 * var_rapid_delayed_rectifying_potassium_current_P_as_gate__P_as)) * var_rapid_delayed_rectifying_potassium_current_P_i_gate__P_i * (var_membrane__V -  -86.6319293974);
    var_membrane__i_K_s = 0.0003445 * pow(var_slow_delayed_rectifying_potassium_current_xs_gate__xs, 2.0) * (var_membrane__V -  -71.3653228522);
    var_membrane__i_f_Na = 0.0005465 * var_hyperpolarisation_activated_current_y_gate__y * (var_membrane__V - 76.1718707053);
    var_membrane__i_f_K = 0.0005465 * var_hyperpolarisation_activated_current_y_gate__y * (var_membrane__V -  -86.6319293974);
    var_membrane__i_b_Na = 5.81818e-05 * (var_membrane__V - 76.1718707053);
    var_membrane__i_b_Ca = 1.3236e-05 * (var_membrane__V - 131.780962407);
    var_membrane__i_b_K = 2.523636e-05 * (var_membrane__V -  -86.6319293974);
    var_membrane__i_NaCa = (2.7229e-06 * ((512.0 * 2.0 * exp(0.03743 * var_membrane__V * 0.5)) - (2744000.0 * 0.0001 * exp(0.03743 * var_membrane__V *  -0.5)))) * 0.885081073426;
    var_membrane__i_p = 0.0124181288853 / (1.5 + exp((-(var_membrane__V + 60.0)) * 0.025));
    var_membrane__i_Ca_p = 0.0;
    const double var_sodium_current_h_gate__h1_infinity = 1.0 / (1.0 + exp((var_membrane__V + 66.1) * 0.15625));
    const double var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af_infinity = 1.0 / (1.0 + exp((-(var_membrane__V + 14.2)) * 0.0943396226415));
    
    const double d_dt_membrane__V =  -50000.0 * (var_membrane__i_Na + var_membrane__i_Ca_L + var_membrane__i_Ca_T + var_membrane__i_to + var_membrane__i_sus + var_membrane__i_K_r + var_membrane__i_K_s + var_membrane__i_f_Na + var_membrane__i_f_K + var_membrane__i_b_Na + var_membrane__i_b_Ca + var_membrane__i_b_K + var_membrane__i_NaCa + var_membrane__i_p + var_membrane__i_Ca_p);
    const double d_dt_sodium_current_m_gate__m = (pow(1.0 / (1.0 + exp((-(var_membrane__V + 30.32)) * 0.18315018315)), 0.333333333333) - var_sodium_current_m_gate__m) / ((0.0006247 / ((0.8322166 * exp( -0.33566 * (var_membrane__V + 56.7062))) + (0.6274 * exp(0.0823 * (var_membrane__V + 65.0131))))) + 4.569e-05);
    const double d_dt_sodium_current_h_gate__h1 = (var_sodium_current_h_gate__h1_infinity - var_sodium_current_h_gate__h1) / (((3.717e-06 * exp( -0.2815 * (var_membrane__V + 17.11))) / (1.0 + (0.003732 * exp( -0.3426 * (var_membrane__V + 37.76))))) + 0.0005977);
    const double d_dt_sodium_current_h_gate__h2 = (var_sodium_current_h_gate__h1_infinity - var_sodium_current_h_gate__h2) / (((3.186e-08 * exp( -0.6219 * (var_membrane__V + 18.8))) / (1.0 + (7.189e-05 * exp( -0.6683 * (var_membrane__V + 34.07))))) + 0.003556);
    const double d_dt_L_type_Ca_channel_d_gate__d_L = ((1.0 / (1.0 + exp((-(var_membrane__V + 22.3 + 0.0)) * 0.166666666667))) - var_L_type_Ca_channel_d_gate__d_L) / (2.0 / (((( -28.39 * (var_membrane__V + 35.0)) / (exp((-(var_membrane__V + 35.0)) * 0.4) - 1.0)) - ((84.9 * var_membrane__V) / (exp( -0.208 * var_membrane__V) - 1.0))) + ((11.43 * (var_membrane__V - 5.0)) / (exp(0.4 * (var_membrane__V - 5.0)) - 1.0))));
    const double d_dt_L_type_Ca_channel_f_gate__f_L = ((1.0 / (1.0 + exp((var_membrane__V + 45.0) * 0.2))) - var_L_type_Ca_channel_f_gate__f_L) / (1.2 / (((3.75 * (var_membrane__V + 28.0)) / (exp((var_membrane__V + 28.0) * 0.25) - 1.0)) + (30.0 / (1.0 + exp((-(var_membrane__V + 28.0)) * 0.25)))));
    const double d_dt_T_type_Ca_channel_d_gate__d_T = ((1.0 / (1.0 + exp((-(var_membrane__V + 37.0)) * 0.147058823529))) - var_T_type_Ca_channel_d_gate__d_T) / (1.0 / ((1068.0 * exp((var_membrane__V + 26.3) * 0.0333333333333)) + (1068.0 * exp((-(var_membrane__V + 26.3)) * 0.0333333333333))));
    const double d_dt_T_type_Ca_channel_f_gate__f_T = ((1.0 / (1.0 + exp((var_membrane__V + 71.0) * 0.111111111111))) - var_T_type_Ca_channel_f_gate__f_T) / (1.0 / ((15.3 * exp((-(var_membrane__V + 71.0 + 0.0)) * 0.0120048019208)) + (15.0 * exp((var_membrane__V + 71.0) * 0.0650195058518))));
    const double d_dt_four_AP_sensitive_currents_q_gate__q = ((1.0 / (1.0 + exp((var_membrane__V + 59.37) * 0.0763358778626))) - var_four_AP_sensitive_currents_q_gate__q) / (0.000333333333333 * (30.31 + (195.5 / ((0.5686 * exp( -0.08161 * (var_membrane__V + 39.0 + 0.0))) + (0.7174 * exp(0.2719 * 1.0 * (var_membrane__V + 40.93 + 0.0)))))));
    const double d_dt_four_AP_sensitive_currents_r_gate__r = ((1.0 / (1.0 + exp((-(var_membrane__V - 10.93)) * 0.0507614213198))) - var_four_AP_sensitive_currents_r_gate__r) / (0.0025 * (1.191 + (7.838 / ((1.037 * exp(0.09012 * (var_membrane__V + 30.61))) + (0.369 * exp( -0.119 * (var_membrane__V + 23.84)))))));
    const double d_dt_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af = (var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af_infinity - var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af) / (1.0 / ((37.2 * exp((var_membrane__V - 9.0) * 0.062893081761)) + (0.96 * exp((-(var_membrane__V - 9.0)) * 0.0444444444444))));
    const double d_dt_rapid_delayed_rectifying_potassium_current_P_as_gate__P_as = (var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af_infinity - var_rapid_delayed_rectifying_potassium_current_P_as_gate__P_as) / (1.0 / ((4.2 * exp((var_membrane__V - 9.0) * 0.0588235294118)) + (0.15 * exp((-(var_membrane__V - 9.0)) * 0.0462962962963))));
    const double d_dt_rapid_delayed_rectifying_potassium_current_P_i_gate__P_i = ((1.0 / (1.0 + exp((var_membrane__V + 18.6) * 0.0990099009901))) - var_rapid_delayed_rectifying_potassium_current_P_i_gate__P_i) * 500.0;
    const double d_dt_slow_delayed_rectifying_potassium_current_xs_gate__xs = ((14.0 / (1.0 + exp((-(var_membrane__V - 40.0)) * 0.111111111111))) * (1.0 - var_slow_delayed_rectifying_potassium_current_xs_gate__xs)) - ((1.0 * exp((-var_membrane__V) * 0.0222222222222)) * var_slow_delayed_rectifying_potassium_current_xs_gate__xs);
    const double d_dt_hyperpolarisation_activated_current_y_gate__y = ((1.0 * exp((-(var_membrane__V + 78.91)) * 0.0375516334961)) * (1.0 - var_hyperpolarisation_activated_current_y_gate__y)) - ((1.0 * exp((var_membrane__V + 75.13) * 0.0470588235294)) * var_hyperpolarisation_activated_current_y_gate__y);
    
    classe->V_lado_direito_ = d_dt_membrane__V;
    classe->m_lado_direito_ = d_dt_sodium_current_m_gate__m;
    classe->h1_lado_direito_ = d_dt_sodium_current_h_gate__h1;
    classe->h2_lado_direito_ = d_dt_sodium_current_h_gate__h2;
    classe->d_L_lado_direito_ = d_dt_L_type_Ca_channel_d_gate__d_L;
    classe->f_L_lado_direito_ = d_dt_L_type_Ca_channel_f_gate__f_L;
    classe->d_T_lado_direito_ = d_dt_T_type_Ca_channel_d_gate__d_T;
    classe->f_T_lado_direito_ = d_dt_T_type_Ca_channel_f_gate__f_T;
    classe->q_lado_direito_ = d_dt_four_AP_sensitive_currents_q_gate__q;
    classe->r_lado_direito_ = d_dt_four_AP_sensitive_currents_r_gate__r;
    classe->P_af_lado_direito_ = d_dt_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af;
    classe->P_as_lado_direito_ = d_dt_rapid_delayed_rectifying_potassium_current_P_as_gate__P_as;
    classe->P_i_lado_direito_ = d_dt_rapid_delayed_rectifying_potassium_current_P_i_gate__P_i;
    classe->xs_lado_direito_ = d_dt_slow_delayed_rectifying_potassium_current_xs_gate__xs;
    classe->y_lado_direito_ = d_dt_hyperpolarisation_activated_current_y_gate__y;
};


#endif
