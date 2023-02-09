#ifndef _CML_zhang_SAN_model_2000_all_pe_lut_
#define _CML_zhang_SAN_model_2000_all_pe_lut_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: zhang_SAN_model_2000_all
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translate: 1844, pycml: 1844)
//! on Sun Jun 26 21:58:01 2011
//! 
//! <autogenerated>



class CML_zhang_SAN_model_2000_all_pe_lut_LookupTables
{
public:
    static CML_zhang_SAN_model_2000_all_pe_lut_LookupTables* Instance()
    {
        if (mpInstance == NULL)
        {
            mpInstance = new CML_zhang_SAN_model_2000_all_pe_lut_LookupTables;
        }
        return mpInstance;
    }
    
    // Methods to look up values from lookup tables
    // using linear interpolation
    inline double _lookup_0(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][0];
        double y2 = _lookup_table_0[i+1][0];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_1(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][1];
        double y2 = _lookup_table_0[i+1][1];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_2(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][2];
        double y2 = _lookup_table_0[i+1][2];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_3(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][3];
        double y2 = _lookup_table_0[i+1][3];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_4(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][4];
        double y2 = _lookup_table_0[i+1][4];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_5(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][5];
        double y2 = _lookup_table_0[i+1][5];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_6(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][6];
        double y2 = _lookup_table_0[i+1][6];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_7(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][7];
        double y2 = _lookup_table_0[i+1][7];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_8(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][8];
        double y2 = _lookup_table_0[i+1][8];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_9(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][9];
        double y2 = _lookup_table_0[i+1][9];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_10(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][10];
        double y2 = _lookup_table_0[i+1][10];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_11(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][11];
        double y2 = _lookup_table_0[i+1][11];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_12(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][12];
        double y2 = _lookup_table_0[i+1][12];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_13(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][13];
        double y2 = _lookup_table_0[i+1][13];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_14(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][14];
        double y2 = _lookup_table_0[i+1][14];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_15(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][15];
        double y2 = _lookup_table_0[i+1][15];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_16(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][16];
        double y2 = _lookup_table_0[i+1][16];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_17(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][17];
        double y2 = _lookup_table_0[i+1][17];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_18(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][18];
        double y2 = _lookup_table_0[i+1][18];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_19(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][19];
        double y2 = _lookup_table_0[i+1][19];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_20(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][20];
        double y2 = _lookup_table_0[i+1][20];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_21(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][21];
        double y2 = _lookup_table_0[i+1][21];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_22(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][22];
        double y2 = _lookup_table_0[i+1][22];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_23(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][23];
        double y2 = _lookup_table_0[i+1][23];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_24(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][24];
        double y2 = _lookup_table_0[i+1][24];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_25(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][25];
        double y2 = _lookup_table_0[i+1][25];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_26(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][26];
        double y2 = _lookup_table_0[i+1][26];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_27(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][27];
        double y2 = _lookup_table_0[i+1][27];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_28(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][28];
        double y2 = _lookup_table_0[i+1][28];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_29(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][29];
        double y2 = _lookup_table_0[i+1][29];
        return y1 + (y2-y1)*factor;
    }
    
    inline double _lookup_30(unsigned i, double factor)
    {
        double y1 = _lookup_table_0[i][30];
        double y2 = _lookup_table_0[i+1][30];
        return y1 + (y2-y1)*factor;
    }
    
    
protected:
    CML_zhang_SAN_model_2000_all_pe_lut_LookupTables(const CML_zhang_SAN_model_2000_all_pe_lut_LookupTables&);
    CML_zhang_SAN_model_2000_all_pe_lut_LookupTables& operator= (const CML_zhang_SAN_model_2000_all_pe_lut_LookupTables&);
    CML_zhang_SAN_model_2000_all_pe_lut_LookupTables()
    {
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][0] = exp(((var_membrane__V - 76.1718707053) * 96845.0) * 3.87996927064e-07) - 1.0;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][1] = exp((var_membrane__V * 96845.0) * 3.87996927064e-07) - 1.0;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][2] = pow(1.0 / (1.0 + exp((-(var_membrane__V + 30.32)) * 0.18315018315)), 0.333333333333);
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][3] = (0.0006247 / ((0.8322166 * exp( -0.33566 * (var_membrane__V + 56.7062))) + (0.6274 * exp(0.0823 * (var_membrane__V + 65.0131))))) + 4.569e-05;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][4] = ((0.09518 * exp( -0.06306 * (var_membrane__V + 34.4))) / (1.0 + (1.662 * exp( -0.2251 * (var_membrane__V + 63.7))))) + 0.08693;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][5] = ((3.717e-06 * exp( -0.2815 * (var_membrane__V + 17.11))) / (1.0 + (0.003732 * exp( -0.3426 * (var_membrane__V + 37.76))))) + 0.0005977;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][6] = ((3.186e-08 * exp( -0.6219 * (var_membrane__V + 18.8))) / (1.0 + (7.189e-05 * exp( -0.6683 * (var_membrane__V + 34.07))))) + 0.003556;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][7] = 1.0 / (1.0 + exp((var_membrane__V + 66.1) * 0.15625));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][8] = 0.006 / (1.0 + exp((-(var_membrane__V + 14.1)) * 0.166666666667));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][9] = 1.0 / (1.0 + exp((-(var_membrane__V + 22.3 + 0.0)) * 0.166666666667));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][10] = 2.0 / (((( -28.39 * (var_membrane__V + 35.0)) / (exp((-(var_membrane__V + 35.0)) * 0.4) - 1.0)) - ((84.9 * var_membrane__V) / (exp( -0.208 * var_membrane__V) - 1.0))) + ((11.43 * (var_membrane__V - 5.0)) / (exp(0.4 * (var_membrane__V - 5.0)) - 1.0)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][11] = 1.0 / (1.0 + exp((var_membrane__V + 45.0) * 0.2));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][12] = 1.2 / (((3.75 * (var_membrane__V + 28.0)) / (exp((var_membrane__V + 28.0) * 0.25) - 1.0)) + (30.0 / (1.0 + exp((-(var_membrane__V + 28.0)) * 0.25))));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][13] = 1.0 / (1.0 + exp((-(var_membrane__V + 37.0)) * 0.147058823529));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][14] = 1.0 / ((1068.0 * exp((var_membrane__V + 26.3) * 0.0333333333333)) + (1068.0 * exp((-(var_membrane__V + 26.3)) * 0.0333333333333)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][15] = 1.0 / (1.0 + exp((var_membrane__V + 71.0) * 0.111111111111));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][16] = 1.0 / ((15.3 * exp((-(var_membrane__V + 71.0 + 0.0)) * 0.0120048019208)) + (15.0 * exp((var_membrane__V + 71.0) * 0.0650195058518)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][17] = 1.0 / (1.0 + exp((var_membrane__V + 59.37) * 0.0763358778626));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][18] = 0.000333333333333 * (30.31 + (195.5 / ((0.5686 * exp( -0.08161 * (var_membrane__V + 39.0 + 0.0))) + (0.7174 * exp(0.2719 * 1.0 * (var_membrane__V + 40.93 + 0.0))))));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][19] = 1.0 / (1.0 + exp((-(var_membrane__V - 10.93)) * 0.0507614213198));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][20] = 0.0025 * (1.191 + (7.838 / ((1.037 * exp(0.09012 * (var_membrane__V + 30.61))) + (0.369 * exp( -0.119 * (var_membrane__V + 23.84))))));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][21] = 1.0 / ((37.2 * exp((var_membrane__V - 9.0) * 0.062893081761)) + (0.96 * exp((-(var_membrane__V - 9.0)) * 0.0444444444444)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][22] = 1.0 / (1.0 + exp((-(var_membrane__V + 14.2)) * 0.0943396226415));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][23] = 1.0 / ((4.2 * exp((var_membrane__V - 9.0) * 0.0588235294118)) + (0.15 * exp((-(var_membrane__V - 9.0)) * 0.0462962962963)));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][24] = 1.0 / (1.0 + exp((var_membrane__V + 18.6) * 0.0990099009901));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][25] = 14.0 / (1.0 + exp((-(var_membrane__V - 40.0)) * 0.111111111111));
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][26] = 1.0 * exp((-var_membrane__V) * 0.0222222222222);
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][27] = 1.0 * exp((-(var_membrane__V + 78.91)) * 0.0375516334961);
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][28] = 1.0 * exp((var_membrane__V + 75.13) * 0.0470588235294);
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][29] = (2.7229e-06 * ((512.0 * 2.0 * exp(0.03743 * var_membrane__V * 0.5)) - (2744000.0 * 0.0001 * exp(0.03743 * var_membrane__V *  -0.5)))) * 0.885081073426;
        }
        
        for (int i=0 ; i<16001; i++)
        {
            double var_membrane__V = -100.0001 + i*0.01;
            _lookup_table_0[i][30] = 0.0124181288853 / (1.5 + exp((-(var_membrane__V + 60.0)) * 0.025));
        }
        
    }
    
private:
    /** The single instance of the class */
    static CML_zhang_SAN_model_2000_all_pe_lut_LookupTables *mpInstance;

    // Lookup tables
    double _lookup_table_0[16001][31];
    
};

CML_zhang_SAN_model_2000_all_pe_lut_LookupTables* CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::mpInstance = NULL;

// 
// Settable parameters and readable variables
// 
void __equation__pe_lut_(Solveode *classe)
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
    
    
    // Lookup table indexing
    if (var_membrane__V>59.9999 || var_membrane__V<-100.0001)
    {
#define COVERAGE_IGNORE
        printf("V outside lookup table range");
#undef COVERAGE_IGNORE
    }
    
    double _offset_0 = var_membrane__V - -100.0001;
    double _offset_0_over_table_step = _offset_0 * 100.0;
    unsigned _table_index_0 = (unsigned)(_offset_0_over_table_step);
    double _factor_0 = _offset_0_over_table_step - _table_index_0;
    
    // Mathematics
    const double var_sodium_current_h_gate__F_Na = CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_4(_table_index_0, _factor_0);
    var_membrane__i_Na = ((((0.0 * pow(var_sodium_current_m_gate__m, 3.0) * (((1.0 - var_sodium_current_h_gate__F_Na) * var_sodium_current_h_gate__h1) + (var_sodium_current_h_gate__F_Na * var_sodium_current_h_gate__h2)) * 140.0 * 9378954025.0) * 3.87996927064e-07) * CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_0(_table_index_0, _factor_0)) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_1(_table_index_0, _factor_0)) * var_membrane__V;
    var_membrane__i_Ca_L = 0.0057938 * ((var_L_type_Ca_channel_f_gate__f_L * var_L_type_Ca_channel_d_gate__d_L) + CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_8(_table_index_0, _factor_0)) * (var_membrane__V - 46.4);
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
    var_membrane__i_NaCa = CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_29(_table_index_0, _factor_0);
    var_membrane__i_p = CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_30(_table_index_0, _factor_0);
    var_membrane__i_Ca_p = 0.0;
    const double var_sodium_current_h_gate__h1_infinity = CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_7(_table_index_0, _factor_0);
    const double var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af_infinity = CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_22(_table_index_0, _factor_0);
    
    const double d_dt_membrane__V =  -50000.0 * (var_membrane__i_Na + var_membrane__i_Ca_L + var_membrane__i_Ca_T + var_membrane__i_to + var_membrane__i_sus + var_membrane__i_K_r + var_membrane__i_K_s + var_membrane__i_f_Na + var_membrane__i_f_K + var_membrane__i_b_Na + var_membrane__i_b_Ca + var_membrane__i_b_K + var_membrane__i_NaCa + var_membrane__i_p + var_membrane__i_Ca_p);
    const double d_dt_sodium_current_m_gate__m = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_2(_table_index_0, _factor_0) - var_sodium_current_m_gate__m) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_3(_table_index_0, _factor_0);
    const double d_dt_sodium_current_h_gate__h1 = (var_sodium_current_h_gate__h1_infinity - var_sodium_current_h_gate__h1) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_5(_table_index_0, _factor_0);
    const double d_dt_sodium_current_h_gate__h2 = (var_sodium_current_h_gate__h1_infinity - var_sodium_current_h_gate__h2) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_6(_table_index_0, _factor_0);
    const double d_dt_L_type_Ca_channel_d_gate__d_L = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_9(_table_index_0, _factor_0) - var_L_type_Ca_channel_d_gate__d_L) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_10(_table_index_0, _factor_0);
    const double d_dt_L_type_Ca_channel_f_gate__f_L = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_11(_table_index_0, _factor_0) - var_L_type_Ca_channel_f_gate__f_L) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_12(_table_index_0, _factor_0);
    const double d_dt_T_type_Ca_channel_d_gate__d_T = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_13(_table_index_0, _factor_0) - var_T_type_Ca_channel_d_gate__d_T) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_14(_table_index_0, _factor_0);
    const double d_dt_T_type_Ca_channel_f_gate__f_T = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_15(_table_index_0, _factor_0) - var_T_type_Ca_channel_f_gate__f_T) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_16(_table_index_0, _factor_0);
    const double d_dt_four_AP_sensitive_currents_q_gate__q = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_17(_table_index_0, _factor_0) - var_four_AP_sensitive_currents_q_gate__q) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_18(_table_index_0, _factor_0);
    const double d_dt_four_AP_sensitive_currents_r_gate__r = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_19(_table_index_0, _factor_0) - var_four_AP_sensitive_currents_r_gate__r) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_20(_table_index_0, _factor_0);
    const double d_dt_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af = (var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af_infinity - var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_21(_table_index_0, _factor_0);
    const double d_dt_rapid_delayed_rectifying_potassium_current_P_as_gate__P_as = (var_rapid_delayed_rectifying_potassium_current_P_af_gate__P_af_infinity - var_rapid_delayed_rectifying_potassium_current_P_as_gate__P_as) / CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_23(_table_index_0, _factor_0);
    const double d_dt_rapid_delayed_rectifying_potassium_current_P_i_gate__P_i = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_24(_table_index_0, _factor_0) - var_rapid_delayed_rectifying_potassium_current_P_i_gate__P_i) * 500.0;
    const double d_dt_slow_delayed_rectifying_potassium_current_xs_gate__xs = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_25(_table_index_0, _factor_0) * (1.0 - var_slow_delayed_rectifying_potassium_current_xs_gate__xs)) - (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_26(_table_index_0, _factor_0) * var_slow_delayed_rectifying_potassium_current_xs_gate__xs);
    const double d_dt_hyperpolarisation_activated_current_y_gate__y = (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_27(_table_index_0, _factor_0) * (1.0 - var_hyperpolarisation_activated_current_y_gate__y)) - (CML_zhang_SAN_model_2000_all_pe_lut_LookupTables::Instance()->_lookup_28(_table_index_0, _factor_0) * var_hyperpolarisation_activated_current_y_gate__y);
    
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
