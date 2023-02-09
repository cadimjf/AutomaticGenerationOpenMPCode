#include "sundials_types.h"
#include "cvode.h"
#include "cvode_dense.h"
#include "nvector_serial.h"
#include "sundials_dense.h"
#include "cvode_diag.h"
#include "cvode_band.h"
#define ABSTOL 1.0E-06
#define RELTOL 1.0E-04
static int check_flag(void *flagvalue, char *funcname, int opt);
static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__);
#include <stdio.h>
#include <math.h>
#include "MCutil.hpp"
#define _H_FACTOR_ 2
#define RECURSION_LIMIT 10
#define _EULER_ 0
 #define _ADAP_DT_ 4
#define _RK2_ 3
#define _ADAP_DT2_ 5
#include <omp.h>
#define _LOWER_LIMIT_   0.2
#define _UPPER_LIMIT_   1.0
/// equationType
#define eqAgos 1
#define eqPyCmlSimples 2
#define eqPyCmlPE 3
#define eqPyCmlLUT 4
#define eqPyCmlPE_LUT 5
#define forestSize 16
#include <sys/time.h>
#include "greedy.cpp"
#define _INCREASE_DT_ 1.5
#define _DECREASE_DT_ 0.65
#define _DECREASE_DT_2_ 2
#define numEDO 41
#define numAux 72

class Solveode
{
public:
	int it_countx;
	int *tree_thread;
	double time; 	 // millisecond
	double stim_amplitude; 	 // picoA_per_picoF
	double stim_start; 	 // millisecond
	double stim_end; 	 // millisecond
	double stim_period; 	 // millisecond
	double stim_duration; 	 // millisecond
	double Acap; 	 // cm2
	double Cm; 	 // microF_per_cm2
	double Vmyo; 	 // microlitre
	double F; 	 // coulomb_per_millimole
	double VJSR; 	 // microlitre
	double Vss; 	 // microlitre
	double VNSR; 	 // microlitre
	double CMDN_tot; 	 // micromolar
	double Km_CMDN; 	 // micromolar
	double CSQN_tot; 	 // micromolar
	double Km_CSQN; 	 // micromolar
	double v1; 	 // per_millisecond
	double tau_tr; 	 // millisecond
	double tau_xfer; 	 // millisecond
	double v2; 	 // per_millisecond
	double v3; 	 // micromolar_per_millisecond
	double Km_up; 	 // micromolar
	double k_plus_htrpn; 	 // per_micromolar_millisecond
	double HTRPN_tot; 	 // micromolar
	double k_plus_ltrpn; 	 // per_micromolar_millisecond
	double LTRPN_tot; 	 // micromolar
	double k_minus_htrpn; 	 // per_millisecond
	double k_minus_ltrpn; 	 // per_millisecond
	double i_CaL_max; 	 // picoA_per_picoF
	double k_plus_a; 	 // micromolar4_per_millisecond
	double n; 	 // dimensionless
	double k_minus_b; 	 // per_millisecond
	double k_minus_c; 	 // per_millisecond
	double k_minus_a; 	 // per_millisecond
	double k_plus_b; 	 // micromolar3_per_millisecond
	double m; 	 // dimensionless
	double k_plus_c; 	 // per_millisecond
	double g_CaL; 	 // milliS_per_microF
	double E_CaL; 	 // millivolt
	double Kpcb; 	 // per_millisecond
	double Kpc_max; 	 // per_millisecond
	double Kpc_half; 	 // micromolar
	double i_pCa_max; 	 // picoA_per_picoF
	double Km_pCa; 	 // micromolar
	double k_NaCa; 	 // picoA_per_picoF
	double K_mNa; 	 // micromolar
	double Nao; 	 // micromolar
	double K_mCa; 	 // micromolar
	double Cao; 	 // micromolar
	double k_sat; 	 // dimensionless
	double eta; 	 // dimensionless
	double R; 	 // joule_per_mole_kelvin
	double T; 	 // kelvin
	double g_Cab; 	 // milliS_per_microF
	double g_Na; 	 // milliS_per_microF
	double Ko; 	 // micromolar
	double g_Nab; 	 // milliS_per_microF
	double g_Kto_f; 	 // milliS_per_microF
	double g_Kto_s; 	 // milliS_per_microF
	double g_Ks; 	 // milliS_per_microF
	double g_Kur; 	 // milliS_per_microF
	double g_Kss; 	 // milliS_per_microF
	double g_Kr; 	 // milliS_per_microF
	double kf; 	 // per_millisecond
	double kb; 	 // per_millisecond
	double i_NaK_max; 	 // picoA_per_picoF
	double Km_Nai; 	 // micromolar
	double Km_Ko; 	 // micromolar
	double g_ClCa; 	 // milliS_per_microF
	double Km_Cl; 	 // micromolar
	double E_Cl; 	 // millivolt
	double calc_i_stim; 	 // picoA_per_picoF
	double calc_Bi; 	 // dimensionless
	double calc_Bss; 	 // dimensionless
	double calc_BJSR; 	 // dimensionless
	double calc_J_rel; 	 // micromolar_per_millisecond
	double calc_J_tr; 	 // micromolar_per_millisecond
	double calc_J_xfer; 	 // micromolar_per_millisecond
	double calc_J_leak; 	 // micromolar_per_millisecond
	double calc_J_up; 	 // micromolar_per_millisecond
	double calc_J_trpn; 	 // micromolar_per_millisecond
	double calc_P_C1; 	 // dimensionless
	double calc_i_CaL; 	 // picoA_per_picoF
	double calc_C1; 	 // dimensionless
	double calc_alpha; 	 // per_millisecond
	double calc_beta; 	 // per_millisecond
	double calc_gamma; 	 // per_millisecond
	double calc_Kpcf; 	 // per_millisecond
	double calc_i_pCa; 	 // picoA_per_picoF
	double calc_i_NaCa; 	 // picoA_per_picoF
	double calc_i_Cab; 	 // picoA_per_picoF
	double calc_E_CaN; 	 // millivolt
	double calc_i_Na; 	 // picoA_per_picoF
	double calc_E_Na; 	 // millivolt
	double calc_C_Na3; 	 // dimensionless
	double calc_alpha_Na11; 	 // per_millisecond
	double calc_alpha_Na12; 	 // per_millisecond
	double calc_alpha_Na13; 	 // per_millisecond
	double calc_beta_Na11; 	 // per_millisecond
	double calc_beta_Na12; 	 // per_millisecond
	double calc_beta_Na13; 	 // per_millisecond
	double calc_alpha_Na3; 	 // per_millisecond
	double calc_beta_Na3; 	 // per_millisecond
	double calc_alpha_Na2; 	 // per_millisecond
	double calc_beta_Na2; 	 // per_millisecond
	double calc_alpha_Na4; 	 // per_millisecond
	double calc_beta_Na4; 	 // per_millisecond
	double calc_alpha_Na5; 	 // per_millisecond
	double calc_beta_Na5; 	 // per_millisecond
	double calc_i_Nab; 	 // picoA_per_picoF
	double calc_i_Kto_f; 	 // picoA_per_picoF
	double calc_E_K; 	 // millivolt
	double calc_alpha_a; 	 // per_millisecond
	double calc_beta_a; 	 // per_millisecond
	double calc_alpha_i; 	 // per_millisecond
	double calc_beta_i; 	 // per_millisecond
	double calc_i_Kto_s; 	 // picoA_per_picoF
	double calc_ass; 	 // dimensionless
	double calc_iss; 	 // dimensionless
	double calc_tau_ta_s; 	 // millisecond
	double calc_tau_ti_s; 	 // millisecond
	double calc_i_K1; 	 // picoA_per_picoF
	double calc_i_Ks; 	 // picoA_per_picoF
	double calc_alpha_n; 	 // per_millisecond
	double calc_beta_n; 	 // per_millisecond
	double calc_i_Kur; 	 // picoA_per_picoF
	double calc_tau_aur; 	 // millisecond
	double calc_tau_iur; 	 // millisecond
	double calc_i_Kss; 	 // picoA_per_picoF
	double calc_tau_Kss; 	 // millisecond
	double calc_i_Kr; 	 // picoA_per_picoF
	double calc_C_K0; 	 // dimensionless
	double calc_alpha_a0; 	 // per_millisecond
	double calc_beta_a0; 	 // per_millisecond
	double calc_alpha_a1; 	 // per_millisecond
	double calc_beta_a1; 	 // per_millisecond
	double calc_i_NaK; 	 // picoA_per_picoF
	double calc_f_NaK; 	 // dimensionless
	double calc_sigma; 	 // dimensionless
	double calc_i_ClCa; 	 // picoA_per_picoF
	double calc_O_ClCa; 	 // dimensionless
	double calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current; 	 // (null)
	double calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current; 	 // (null)
	double dtime, *time_vec__;
	double time_new;

	//functions variables
	double *V;
	double V_new_, V_old_, V_ini_, V_lado_direito_;
	double *Cai;
	double Cai_new_, Cai_old_, Cai_ini_, Cai_lado_direito_;
	double *Cass;
	double Cass_new_, Cass_old_, Cass_ini_, Cass_lado_direito_;
	double *CaJSR;
	double CaJSR_new_, CaJSR_old_, CaJSR_ini_, CaJSR_lado_direito_;
	double *CaNSR;
	double CaNSR_new_, CaNSR_old_, CaNSR_ini_, CaNSR_lado_direito_;
	double *P_RyR;
	double P_RyR_new_, P_RyR_old_, P_RyR_ini_, P_RyR_lado_direito_;
	double *LTRPN_Ca;
	double LTRPN_Ca_new_, LTRPN_Ca_old_, LTRPN_Ca_ini_, LTRPN_Ca_lado_direito_;
	double *HTRPN_Ca;
	double HTRPN_Ca_new_, HTRPN_Ca_old_, HTRPN_Ca_ini_, HTRPN_Ca_lado_direito_;
	double *P_O1;
	double P_O1_new_, P_O1_old_, P_O1_ini_, P_O1_lado_direito_;
	double *P_O2;
	double P_O2_new_, P_O2_old_, P_O2_ini_, P_O2_lado_direito_;
	double *P_C2;
	double P_C2_new_, P_C2_old_, P_C2_ini_, P_C2_lado_direito_;
	double *O;
	double O_new_, O_old_, O_ini_, O_lado_direito_;
	double *C2;
	double C2_new_, C2_old_, C2_ini_, C2_lado_direito_;
	double *C3;
	double C3_new_, C3_old_, C3_ini_, C3_lado_direito_;
	double *C4;
	double C4_new_, C4_old_, C4_ini_, C4_lado_direito_;
	double *I1;
	double I1_new_, I1_old_, I1_ini_, I1_lado_direito_;
	double *I2;
	double I2_new_, I2_old_, I2_ini_, I2_lado_direito_;
	double *I3;
	double I3_new_, I3_old_, I3_ini_, I3_lado_direito_;
	double *Nai;
	double Nai_new_, Nai_old_, Nai_ini_, Nai_lado_direito_;
	double *C_Na2;
	double C_Na2_new_, C_Na2_old_, C_Na2_ini_, C_Na2_lado_direito_;
	double *C_Na1;
	double C_Na1_new_, C_Na1_old_, C_Na1_ini_, C_Na1_lado_direito_;
	double *O_Na;
	double O_Na_new_, O_Na_old_, O_Na_ini_, O_Na_lado_direito_;
	double *IF_Na;
	double IF_Na_new_, IF_Na_old_, IF_Na_ini_, IF_Na_lado_direito_;
	double *I1_Na;
	double I1_Na_new_, I1_Na_old_, I1_Na_ini_, I1_Na_lado_direito_;
	double *I2_Na;
	double I2_Na_new_, I2_Na_old_, I2_Na_ini_, I2_Na_lado_direito_;
	double *IC_Na2;
	double IC_Na2_new_, IC_Na2_old_, IC_Na2_ini_, IC_Na2_lado_direito_;
	double *IC_Na3;
	double IC_Na3_new_, IC_Na3_old_, IC_Na3_ini_, IC_Na3_lado_direito_;
	double *Ki;
	double Ki_new_, Ki_old_, Ki_ini_, Ki_lado_direito_;
	double *ato_f;
	double ato_f_new_, ato_f_old_, ato_f_ini_, ato_f_lado_direito_;
	double *ito_f;
	double ito_f_new_, ito_f_old_, ito_f_ini_, ito_f_lado_direito_;
	double *ato_s;
	double ato_s_new_, ato_s_old_, ato_s_ini_, ato_s_lado_direito_;
	double *ito_s;
	double ito_s_new_, ito_s_old_, ito_s_ini_, ito_s_lado_direito_;
	double *nKs;
	double nKs_new_, nKs_old_, nKs_ini_, nKs_lado_direito_;
	double *aur;
	double aur_new_, aur_old_, aur_ini_, aur_lado_direito_;
	double *iur;
	double iur_new_, iur_old_, iur_ini_, iur_lado_direito_;
	double *aKss;
	double aKss_new_, aKss_old_, aKss_ini_, aKss_lado_direito_;
	double *iKss;
	double iKss_new_, iKss_old_, iKss_ini_, iKss_lado_direito_;
	double *C_K2;
	double C_K2_new_, C_K2_old_, C_K2_ini_, C_K2_lado_direito_;
	double *C_K1;
	double C_K1_new_, C_K1_old_, C_K1_ini_, C_K1_lado_direito_;
	double *O_K;
	double O_K_new_, O_K_old_, O_K_ini_, O_K_lado_direito_;
	double *I_K;
	double I_K_new_, I_K_old_, I_K_ini_, I_K_lado_direito_;

public:
	Solveode(double,double,int);
	~Solveode();
	double STEP_TOLERANCE_;
	void getRightHandSide();
	void getRightHandSideParallel();
	void setMaxStep(double);
	int getErrorCode(double, double);
	double solveToFileOMP(double finalTime = 0, double savingRate = 0, int method=1, char *fileName="output.dat", int nThreads=1, double tol=0.2);
	void explicitEulerStep();
	void rungekutta2_Parallel( );
	void rungekutta2_Single( );
	void adaptiveDt_Parallel();
	int adaptiveDt_Single();
	int adaptiveDt();
	double errorAux;
	int count1; int count0;
	int count3; int count4;
	int recursion_counter;
	double maxStep;
	// CVODE VARIABLES
	realtype reltol__, abstol__;
	void *cvode_mem_cvode__;
	N_Vector dependent_variable__;
	int flag__, flagr__;
	double *depvar__;
	int setVariables(int, double);
	int setParameters(int, double);
	int setFreeVariable(double);
	double getVariables(int);
	double getLadoDireito(int);
	void executeTree(void* tree, int treeSize);
	double getParameters(int);
	double getFreeVariable();
	Variables get_Parameters();
	Variables get_Variables();
	Variables get_FreeVariable();
	void setParametersFromFile(char*);
	void setVariablesFromFile(char*);
	void setFreeVariableFromFile(char*);
	double solve(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0, int numThreads=1);
	double jacobian(double final, double svRate);
	double* getSolution(int indVariable);
	double* getIndependentVar();
	double solveToFile(char *filename, char *fileaccess = "", int firstcall__ = 0, int num_iterations__ = 0, int num_results__ = 0);
public:
	inline double ifnumber_0();
	void reInitCVODE();
	void setCVODEMaxStep(double maxstep);
	void solveCVODE(int firstcall__ = 0, int steps__ = 0, int num_results__ =0, int method__=0,char *fileName="saida.out", int cv_method__=1);
	void solver(int method=0,double finalTime=0, double svRate=0, char* fileName="", int threads=1);
private:
	double errorAD_;
	double edos_euler_[numEDO];
	double edos_rk2_[numEDO];
	double edos_fn_[numEDO];
	double ladodir_fn_[numEDO];
	void save_step(FILE *fileptr, int method);
	void euler(double finalTime, FILE *file);
	void euler_OMP(double finalTime, FILE *fileptr, int nThreads);
	void rungeKutta2ndOrder(double finalTime, FILE *fileptr);
	void addt(double finalTime, FILE *fileptr);
	void addt2(double finalTime, FILE *fileptr);
	void addt_OMP(double finalTime, FILE *fileptr, int nThreads);
	void addt2_OMP(double finalTime, FILE *fileptr, int nThreads);
	double timeSaving;
	double previous_dt;
	double savingRate;
};

#define AGOS_NAN (0.0/0.0)
#define AGOS_INF (1.0/0.0)
#define __agos_xor(a,b) (!(a && b) && (a || b))
float __agos_factorial(int);
double _agos_max(double*,int);
double _agos_min(double*,int);
double _agos_round( double, int);
 typedef struct str_funcao{
	void (*funcao)(Solveode*);
}typ_funcao;

typedef struct str_forest{
typ_funcao* tree;
int treeSize;
}typ_forest;

void __tree1__( Solveode *__AGOS){
	__AGOS->calc_i_stim = (((__AGOS->time_new>=__AGOS->stim_start)&&(__AGOS->time_new<=__AGOS->stim_end)&&(((__AGOS->time_new-__AGOS->stim_start)-(floor(((__AGOS->time_new-__AGOS->stim_start)/__AGOS->stim_period))*__AGOS->stim_period))<=__AGOS->stim_duration)))
?(__AGOS->stim_amplitude)
:(0.0000000000e+00);
	__AGOS->calc_Bi = pow((1.0000000000e+00+((__AGOS->CMDN_tot*__AGOS->Km_CMDN)/pow((__AGOS->Km_CMDN+__AGOS->Cai_old_),2.0000000000e+00))),(-1.0000000000e+00));
	__AGOS->calc_Bss = pow((1.0000000000e+00+((__AGOS->CMDN_tot*__AGOS->Km_CMDN)/pow((__AGOS->Km_CMDN+__AGOS->Cass_old_),2.0000000000e+00))),(-1.0000000000e+00));
	__AGOS->calc_BJSR = pow((1.0000000000e+00+((__AGOS->CSQN_tot*__AGOS->Km_CSQN)/pow((__AGOS->Km_CSQN+__AGOS->CaJSR_old_),2.0000000000e+00))),(-1.0000000000e+00));
	__AGOS->calc_J_rel = (__AGOS->v1*(__AGOS->P_O1_old_+__AGOS->P_O2_old_)*(__AGOS->CaJSR_old_-__AGOS->Cass_old_)*__AGOS->P_RyR_old_);
	__AGOS->calc_J_tr = ((__AGOS->CaNSR_old_-__AGOS->CaJSR_old_)/__AGOS->tau_tr);
	__AGOS->calc_J_xfer = ((__AGOS->Cass_old_-__AGOS->Cai_old_)/__AGOS->tau_xfer);
	__AGOS->calc_J_leak = (__AGOS->v2*(__AGOS->CaNSR_old_-__AGOS->Cai_old_));
	__AGOS->calc_J_up = ((__AGOS->v3*pow(__AGOS->Cai_old_,2.0000000000e+00))/(pow(__AGOS->Km_up,2.0000000000e+00)+pow(__AGOS->Cai_old_,2.0000000000e+00)));
	__AGOS->calc_J_trpn = (((__AGOS->k_plus_htrpn*__AGOS->Cai_old_*(__AGOS->HTRPN_tot-__AGOS->HTRPN_Ca_old_))+(__AGOS->k_plus_ltrpn*__AGOS->Cai_old_*(__AGOS->LTRPN_tot-__AGOS->LTRPN_Ca_old_)))-((__AGOS->k_minus_htrpn*__AGOS->HTRPN_Ca_old_)+(__AGOS->k_minus_ltrpn*__AGOS->LTRPN_Ca_old_)));
	__AGOS->calc_i_CaL = (__AGOS->g_CaL*__AGOS->O_old_*(__AGOS->V_old_-__AGOS->E_CaL));
	__AGOS->calc_i_pCa = ((__AGOS->i_pCa_max*pow(__AGOS->Cai_old_,2.0000000000e+00))/(pow(__AGOS->Km_pCa,2.0000000000e+00)+pow(__AGOS->Cai_old_,2.0000000000e+00)));
	__AGOS->calc_i_NaCa = (((((((__AGOS->k_NaCa*1.0000000000e+00)/(pow(__AGOS->K_mNa,3.0000000000e+00)+pow(__AGOS->Nao,3.0000000000e+00)))*1.0000000000e+00)/(__AGOS->K_mCa+__AGOS->Cao))*1.0000000000e+00)/(1.0000000000e+00+(__AGOS->k_sat*exp((((__AGOS->eta-1.0000000000e+00)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T))))))*((exp(((__AGOS->eta*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))*pow(__AGOS->Nai_old_,3.0000000000e+00)*__AGOS->Cao)-(exp((((__AGOS->eta-1.0000000000e+00)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))*pow(__AGOS->Nao,3.0000000000e+00)*__AGOS->Cai_old_)));
	__AGOS->calc_E_CaN = (((__AGOS->R*__AGOS->T)/(2.0000000000e+00*__AGOS->F))*log((__AGOS->Cao/__AGOS->Cai_old_)));
	__AGOS->calc_E_Na = (((__AGOS->R*__AGOS->T)/__AGOS->F)*log((((9.0000000000e-01*__AGOS->Nao)+(1.0000000000e-01*__AGOS->Ko))/((9.0000000000e-01*__AGOS->Nai_old_)+(1.0000000000e-01*__AGOS->Ki_old_)))));
	__AGOS->calc_E_K = (((__AGOS->R*__AGOS->T)/__AGOS->F)*log((__AGOS->Ko/__AGOS->Ki_old_)));
	__AGOS->calc_i_Kr = (__AGOS->g_Kr*__AGOS->O_K_old_*(__AGOS->V_old_-(((__AGOS->R*__AGOS->T)/__AGOS->F)*log((((9.8000000000e-01*__AGOS->Ko)+(2.0000000000e-02*__AGOS->Nao))/((9.8000000000e-01*__AGOS->Ki_old_)+(2.0000000000e-02*__AGOS->Nai_old_)))))));
	__AGOS->calc_sigma = ((1.0000000000e+00/7.0000000000e+00)*(exp((__AGOS->Nao/6.7300000000e+04))-1.0000000000e+00));
	__AGOS->calc_O_ClCa = (2.0000000000e-01/(1.0000000000e+00+exp(((-(__AGOS->V_old_-4.6700000000e+01))/7.8000000000e+00))));
	__AGOS->calc_i_Nab = (__AGOS->g_Nab*(__AGOS->V_old_-__AGOS->calc_E_Na));
	__AGOS->calc_i_Kto_s = (__AGOS->g_Kto_s*__AGOS->ato_s_old_*__AGOS->ito_s_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_K1 = ((((2.9380000000e-01*__AGOS->Ko)/(__AGOS->Ko+2.1000000000e+02))*(__AGOS->V_old_-__AGOS->calc_E_K))/(1.0000000000e+00+exp((8.9600000000e-02*(__AGOS->V_old_-__AGOS->calc_E_K)))));
	__AGOS->calc_i_Ks = (__AGOS->g_Ks*pow(__AGOS->nKs_old_,2.0000000000e+00)*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_Kur = (__AGOS->g_Kur*__AGOS->aur_old_*__AGOS->iur_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_Kss = (__AGOS->g_Kss*__AGOS->aKss_old_*__AGOS->iKss_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_Cab = (__AGOS->g_Cab*(__AGOS->V_old_-__AGOS->calc_E_CaN));
	__AGOS->calc_i_Na = (__AGOS->g_Na*__AGOS->O_Na_old_*(__AGOS->V_old_-__AGOS->calc_E_Na));
	__AGOS->calc_i_Kto_f = (__AGOS->g_Kto_f*pow(__AGOS->ato_f_old_,3.0000000000e+00)*__AGOS->ito_f_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_f_NaK = (1.0000000000e+00/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T))))+(3.6500000000e-02*__AGOS->calc_sigma*exp((((-__AGOS->V_old_)*__AGOS->F)/(__AGOS->R*__AGOS->T))))));
	__AGOS->calc_i_ClCa = (((__AGOS->g_ClCa*__AGOS->calc_O_ClCa*__AGOS->Cai_old_)/(__AGOS->Cai_old_+__AGOS->Km_Cl))*(__AGOS->V_old_-__AGOS->E_Cl));
	__AGOS->calc_i_NaK = ((((__AGOS->i_NaK_max*__AGOS->calc_f_NaK*1.0000000000e+00)/(1.0000000000e+00+pow((__AGOS->Km_Nai/__AGOS->Nai_old_),1.5000000000e+00)))*__AGOS->Ko)/(__AGOS->Ko+__AGOS->Km_Ko));
	__AGOS->V_lado_direito_= (-(__AGOS->calc_i_CaL+__AGOS->calc_i_pCa+__AGOS->calc_i_NaCa+__AGOS->calc_i_Cab+__AGOS->calc_i_Na+__AGOS->calc_i_Nab+__AGOS->calc_i_NaK+__AGOS->calc_i_Kto_f+__AGOS->calc_i_Kto_s+__AGOS->calc_i_K1+__AGOS->calc_i_Ks+__AGOS->calc_i_Kur+__AGOS->calc_i_Kss+__AGOS->calc_i_Kr+__AGOS->calc_i_ClCa+__AGOS->calc_i_stim));
	__AGOS->Cai_lado_direito_= (__AGOS->calc_Bi*((__AGOS->calc_J_leak+__AGOS->calc_J_xfer)-(__AGOS->calc_J_up+__AGOS->calc_J_trpn+((((__AGOS->calc_i_Cab+__AGOS->calc_i_pCa)-(2.0000000000e+00*__AGOS->calc_i_NaCa))*__AGOS->Acap*__AGOS->Cm)/(2.0000000000e+00*__AGOS->Vmyo*__AGOS->F)))));
	__AGOS->Cass_lado_direito_= (__AGOS->calc_Bss*(((__AGOS->calc_J_rel*__AGOS->VJSR)/__AGOS->Vss)-(((__AGOS->calc_J_xfer*__AGOS->Vmyo)/__AGOS->Vss)+((__AGOS->calc_i_CaL*__AGOS->Acap*__AGOS->Cm)/(2.0000000000e+00*__AGOS->Vss*__AGOS->F)))));
	__AGOS->CaJSR_lado_direito_= (__AGOS->calc_BJSR*(__AGOS->calc_J_tr-__AGOS->calc_J_rel));
	__AGOS->CaNSR_lado_direito_= ((((__AGOS->calc_J_up-__AGOS->calc_J_leak)*__AGOS->Vmyo)/__AGOS->VNSR)-((__AGOS->calc_J_tr*__AGOS->VJSR)/__AGOS->VNSR));
	__AGOS->P_RyR_lado_direito_= (((-4.0000000000e-02)*__AGOS->P_RyR_old_)-(((1.0000000000e-01*__AGOS->calc_i_CaL)/__AGOS->i_CaL_max)*exp(((-pow((__AGOS->V_old_-5.0000000000e+00),2.0000000000e+00))/6.4800000000e+02))));
	__AGOS->Nai_lado_direito_= (((-(__AGOS->calc_i_Na+__AGOS->calc_i_Nab+(3.0000000000e+00*__AGOS->calc_i_NaK)+(3.0000000000e+00*__AGOS->calc_i_NaCa)))*__AGOS->Acap*__AGOS->Cm)/(__AGOS->Vmyo*__AGOS->F));
	__AGOS->Ki_lado_direito_= (((-((__AGOS->calc_i_Kto_f+__AGOS->calc_i_Kto_s+__AGOS->calc_i_K1+__AGOS->calc_i_Ks+__AGOS->calc_i_Kss+__AGOS->calc_i_Kur+__AGOS->calc_i_Kr)-(2.0000000000e+00*__AGOS->calc_i_NaK)))*__AGOS->Acap*__AGOS->Cm)/(__AGOS->Vmyo*__AGOS->F));
} //fim

void __tree2__( Solveode *__AGOS){
	__AGOS->calc_P_C1 = (1.0000000000e+00-(__AGOS->P_C2_old_+__AGOS->P_O1_old_+__AGOS->P_O2_old_));
	__AGOS->P_O1_lado_direito_= (((__AGOS->k_plus_a*pow(__AGOS->Cass_old_,__AGOS->n)*__AGOS->calc_P_C1)+(__AGOS->k_minus_b*__AGOS->P_O2_old_)+(__AGOS->k_minus_c*__AGOS->P_C2_old_))-((__AGOS->k_minus_a*__AGOS->P_O1_old_)+(__AGOS->k_plus_b*pow(__AGOS->Cass_old_,__AGOS->m)*__AGOS->P_O1_old_)+(__AGOS->k_plus_c*__AGOS->P_O1_old_)));
} //fim

void __tree3__( Solveode *__AGOS){
	__AGOS->calc_C1 = (1.0000000000e+00-(__AGOS->O_old_+__AGOS->C2_old_+__AGOS->C2_old_+__AGOS->C3_old_+__AGOS->C4_old_+__AGOS->I1_old_+__AGOS->I2_old_+__AGOS->I3_old_));
	__AGOS->calc_alpha = ((4.0000000000e-01*exp(((__AGOS->V_old_+1.2000000000e+01)/1.0000000000e+01))*((1.0000000000e+00+(7.0000000000e-01*exp(((-pow((__AGOS->V_old_+4.0000000000e+01),2.0000000000e+00))/1.0000000000e+01))))-(7.5000000000e-01*exp(((-pow((__AGOS->V_old_+2.0000000000e+01),2.0000000000e+00))/4.0000000000e+02)))))/(1.0000000000e+00+(1.2000000000e-01*exp(((__AGOS->V_old_+1.2000000000e+01)/1.0000000000e+01)))));
	__AGOS->calc_beta = (5.0000000000e-02*exp(((-(__AGOS->V_old_+1.2000000000e+01))/1.3000000000e+01)));
	__AGOS->calc_gamma = ((__AGOS->Kpc_max*__AGOS->Cass_old_)/(__AGOS->Kpc_half+__AGOS->Cass_old_));
	__AGOS->calc_Kpcf = (1.3000000000e+01*(1.0000000000e+00-exp(((-pow((__AGOS->V_old_+1.4500000000e+01),2.0000000000e+00))/1.0000000000e+02))));
	__AGOS->O_lado_direito_= (((__AGOS->calc_alpha*__AGOS->C4_old_)+(__AGOS->Kpcb*__AGOS->I1_old_)+(1.0000000000e-03*((__AGOS->calc_alpha*__AGOS->I2_old_)-(__AGOS->calc_Kpcf*__AGOS->O_old_))))-((4.0000000000e+00*__AGOS->calc_beta*__AGOS->O_old_)+(__AGOS->calc_gamma*__AGOS->O_old_)));
	__AGOS->C2_lado_direito_= (((4.0000000000e+00*__AGOS->calc_alpha*__AGOS->calc_C1)+(2.0000000000e+00*__AGOS->calc_beta*__AGOS->C3_old_))-((__AGOS->calc_beta*__AGOS->C2_old_)+(3.0000000000e+00*__AGOS->calc_alpha*__AGOS->C2_old_)));
	__AGOS->C3_lado_direito_= (((3.0000000000e+00*__AGOS->calc_alpha*__AGOS->C2_old_)+(3.0000000000e+00*__AGOS->calc_beta*__AGOS->C4_old_))-((2.0000000000e+00*__AGOS->calc_beta*__AGOS->C3_old_)+(2.0000000000e+00*__AGOS->calc_alpha*__AGOS->C3_old_)));
	__AGOS->C4_lado_direito_= (((2.0000000000e+00*__AGOS->calc_alpha*__AGOS->C3_old_)+(4.0000000000e+00*__AGOS->calc_beta*__AGOS->O_old_)+(1.0000000000e-02*((4.0000000000e+00*__AGOS->Kpcb*__AGOS->calc_beta*__AGOS->I1_old_)-(__AGOS->calc_alpha*__AGOS->calc_gamma*__AGOS->C4_old_)))+(2.0000000000e-03*((4.0000000000e+00*__AGOS->calc_beta*__AGOS->I2_old_)-(__AGOS->calc_Kpcf*__AGOS->C4_old_)))+(4.0000000000e+00*__AGOS->calc_beta*__AGOS->Kpcb*__AGOS->I3_old_))-((3.0000000000e+00*__AGOS->calc_beta*__AGOS->C4_old_)+(__AGOS->calc_alpha*__AGOS->C4_old_)+(1.0000000000e+00*__AGOS->calc_gamma*__AGOS->calc_Kpcf*__AGOS->C4_old_)));
	__AGOS->I1_lado_direito_= (((__AGOS->calc_gamma*__AGOS->O_old_)+(1.0000000000e-03*((__AGOS->calc_alpha*__AGOS->I3_old_)-(__AGOS->calc_Kpcf*__AGOS->I1_old_)))+(1.0000000000e-02*((__AGOS->calc_alpha*__AGOS->calc_gamma*__AGOS->C4_old_)-(4.0000000000e+00*__AGOS->calc_beta*__AGOS->calc_Kpcf*__AGOS->I1_old_))))-(__AGOS->Kpcb*__AGOS->I1_old_));
	__AGOS->I2_lado_direito_= (((1.0000000000e-03*((__AGOS->calc_Kpcf*__AGOS->O_old_)-(__AGOS->calc_alpha*__AGOS->I2_old_)))+(__AGOS->Kpcb*__AGOS->I3_old_)+(2.0000000000e-03*((__AGOS->calc_Kpcf*__AGOS->C4_old_)-(4.0000000000e+00*__AGOS->calc_beta*__AGOS->I2_old_))))-(__AGOS->calc_gamma*__AGOS->I2_old_));
	__AGOS->I3_lado_direito_= (((1.0000000000e-03*((__AGOS->calc_Kpcf*__AGOS->I1_old_)-(__AGOS->calc_alpha*__AGOS->I3_old_)))+(__AGOS->calc_gamma*__AGOS->I2_old_)+(1.0000000000e+00*__AGOS->calc_gamma*__AGOS->calc_Kpcf*__AGOS->C4_old_))-((4.0000000000e+00*__AGOS->calc_beta*__AGOS->Kpcb*__AGOS->I3_old_)+(__AGOS->Kpcb*__AGOS->I3_old_)));
} //fim

void __tree4__( Solveode *__AGOS){
	__AGOS->calc_C_Na3 = (1.0000000000e+00-(__AGOS->O_Na_old_+__AGOS->C_Na1_old_+__AGOS->C_Na2_old_+__AGOS->IF_Na_old_+__AGOS->I1_Na_old_+__AGOS->I2_Na_old_+__AGOS->IC_Na2_old_+__AGOS->IC_Na3_old_));
	__AGOS->calc_alpha_Na11 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__AGOS->V_old_+2.5000000000e+00))/1.7000000000e+01)))+(2.0000000000e-01*exp(((-(__AGOS->V_old_+2.5000000000e+00))/1.5000000000e+02)))));
	__AGOS->calc_alpha_Na12 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__AGOS->V_old_+2.5000000000e+00))/1.5000000000e+01)))+(2.3000000000e-01*exp(((-(__AGOS->V_old_+2.5000000000e+00))/1.5000000000e+02)))));
	__AGOS->calc_alpha_Na13 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__AGOS->V_old_+2.5000000000e+00))/1.2000000000e+01)))+(2.5000000000e-01*exp(((-(__AGOS->V_old_+2.5000000000e+00))/1.5000000000e+02)))));
	__AGOS->calc_beta_Na11 = (1.9170000000e-01*exp(((-(__AGOS->V_old_+2.5000000000e+00))/2.0300000000e+01)));
	__AGOS->calc_beta_Na12 = (2.0000000000e-01*exp(((-(__AGOS->V_old_-2.5000000000e+00))/2.0300000000e+01)));
	__AGOS->calc_beta_Na13 = (2.2000000000e-01*exp(((-(__AGOS->V_old_-7.5000000000e+00))/2.0300000000e+01)));
	__AGOS->calc_alpha_Na3 = (7.0000000000e-07*exp(((-(__AGOS->V_old_+7.0000000000e+00))/7.7000000000e+00)));
	__AGOS->calc_beta_Na3 = (8.4000000000e-03+(2.0000000000e-05*(__AGOS->V_old_+7.0000000000e+00)));
	__AGOS->calc_alpha_Na2 = (1.0000000000e+00/((1.8849500000e-01*exp(((-(__AGOS->V_old_+7.0000000000e+00))/1.6600000000e+01)))+3.9395600000e-01));
	__AGOS->calc_beta_Na2 = ((__AGOS->calc_alpha_Na13*__AGOS->calc_alpha_Na2*__AGOS->calc_alpha_Na3)/(__AGOS->calc_beta_Na13*__AGOS->calc_beta_Na3));
	__AGOS->calc_alpha_Na4 = (__AGOS->calc_alpha_Na2/1.0000000000e+03);
	__AGOS->calc_beta_Na4 = __AGOS->calc_alpha_Na3;
	__AGOS->calc_alpha_Na5 = (__AGOS->calc_alpha_Na2/9.5000000000e+04);
	__AGOS->calc_beta_Na5 = (__AGOS->calc_alpha_Na3/5.0000000000e+01);
	__AGOS->C_Na2_lado_direito_= (((__AGOS->calc_alpha_Na11*__AGOS->calc_C_Na3)+(__AGOS->calc_beta_Na12*__AGOS->C_Na1_old_)+(__AGOS->calc_alpha_Na3*__AGOS->IC_Na2_old_))-((__AGOS->calc_beta_Na11*__AGOS->C_Na2_old_)+(__AGOS->calc_alpha_Na12*__AGOS->C_Na2_old_)+(__AGOS->calc_beta_Na3*__AGOS->C_Na2_old_)));
	__AGOS->C_Na1_lado_direito_= (((__AGOS->calc_alpha_Na12*__AGOS->C_Na2_old_)+(__AGOS->calc_beta_Na13*__AGOS->O_Na_old_)+(__AGOS->calc_alpha_Na3*__AGOS->IF_Na_old_))-((__AGOS->calc_beta_Na12*__AGOS->C_Na1_old_)+(__AGOS->calc_alpha_Na13*__AGOS->C_Na1_old_)+(__AGOS->calc_beta_Na3*__AGOS->C_Na1_old_)));
	__AGOS->O_Na_lado_direito_= (((__AGOS->calc_alpha_Na13*__AGOS->C_Na1_old_)+(__AGOS->calc_beta_Na2*__AGOS->IF_Na_old_))-((__AGOS->calc_beta_Na13*__AGOS->O_Na_old_)+(__AGOS->calc_alpha_Na2*__AGOS->O_Na_old_)));
	__AGOS->IF_Na_lado_direito_= (((__AGOS->calc_alpha_Na2*__AGOS->O_Na_old_)+(__AGOS->calc_beta_Na3*__AGOS->C_Na1_old_)+(__AGOS->calc_beta_Na4*__AGOS->I1_Na_old_)+(__AGOS->calc_alpha_Na12*__AGOS->IC_Na2_old_))-((__AGOS->calc_beta_Na2*__AGOS->IF_Na_old_)+(__AGOS->calc_alpha_Na3*__AGOS->IF_Na_old_)+(__AGOS->calc_alpha_Na4*__AGOS->IF_Na_old_)+(__AGOS->calc_beta_Na12*__AGOS->IF_Na_old_)));
	__AGOS->I1_Na_lado_direito_= (((__AGOS->calc_alpha_Na4*__AGOS->IF_Na_old_)+(__AGOS->calc_beta_Na5*__AGOS->I2_Na_old_))-((__AGOS->calc_beta_Na4*__AGOS->I1_Na_old_)+(__AGOS->calc_alpha_Na5*__AGOS->I1_Na_old_)));
	__AGOS->I2_Na_lado_direito_= ((__AGOS->calc_alpha_Na5*__AGOS->I1_Na_old_)-(__AGOS->calc_beta_Na5*__AGOS->I2_Na_old_));
	__AGOS->IC_Na2_lado_direito_= (((__AGOS->calc_alpha_Na11*__AGOS->IC_Na3_old_)+(__AGOS->calc_beta_Na12*__AGOS->IF_Na_old_)+(__AGOS->calc_beta_Na3*__AGOS->IC_Na2_old_))-((__AGOS->calc_beta_Na11*__AGOS->IC_Na2_old_)+(__AGOS->calc_alpha_Na12*__AGOS->IC_Na2_old_)+(__AGOS->calc_alpha_Na3*__AGOS->IC_Na2_old_)));
	__AGOS->IC_Na3_lado_direito_= (((__AGOS->calc_beta_Na11*__AGOS->IC_Na2_old_)+(__AGOS->calc_beta_Na3*__AGOS->calc_C_Na3))-((__AGOS->calc_alpha_Na11*__AGOS->IC_Na3_old_)+(__AGOS->calc_alpha_Na3*__AGOS->IC_Na3_old_)));
} //fim

void __tree5__( Solveode *__AGOS){
	__AGOS->calc_alpha_a = (1.8064000000e-01*exp((3.5770000000e-02*(__AGOS->V_old_+3.0000000000e+01))));
	__AGOS->calc_beta_a = (3.9560000000e-01*exp(((-6.2370000000e-02)*(__AGOS->V_old_+3.0000000000e+01))));
	__AGOS->ato_f_lado_direito_= ((__AGOS->calc_alpha_a*(1.0000000000e+00-__AGOS->ato_f_old_))-(__AGOS->calc_beta_a*__AGOS->ato_f_old_));
} //fim

void __tree6__( Solveode *__AGOS){
	__AGOS->calc_alpha_i = ((1.5200000000e-04*exp(((-(__AGOS->V_old_+1.3500000000e+01))/7.0000000000e+00)))/((6.7083000000e-03*exp(((-(__AGOS->V_old_+3.3500000000e+01))/7.0000000000e+00)))+1.0000000000e+00));
	__AGOS->calc_beta_i = ((9.5000000000e-04*exp(((__AGOS->V_old_+3.3500000000e+01)/7.0000000000e+00)))/((5.1335000000e-02*exp(((__AGOS->V_old_+3.3500000000e+01)/7.0000000000e+00)))+1.0000000000e+00));
	__AGOS->ito_f_lado_direito_= ((__AGOS->calc_alpha_i*(1.0000000000e+00-__AGOS->ito_f_old_))-(__AGOS->calc_beta_i*__AGOS->ito_f_old_));
} //fim

void __tree7__( Solveode *__AGOS){
	__AGOS->calc_ass = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+2.2500000000e+01))/7.7000000000e+00))));
	__AGOS->calc_tau_ta_s = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__AGOS->V_old_)))+2.0580000000e+00);
	__AGOS->calc_tau_aur = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__AGOS->V_old_)))+2.0580000000e+00);
	__AGOS->calc_tau_Kss = ((3.9300000000e+01*exp(((-8.6200000000e-02)*__AGOS->V_old_)))+1.3170000000e+01);
	__AGOS->ato_s_lado_direito_= ((__AGOS->calc_ass-__AGOS->ato_s_old_)/__AGOS->calc_tau_ta_s);
	__AGOS->aur_lado_direito_= ((__AGOS->calc_ass-__AGOS->aur_old_)/__AGOS->calc_tau_aur);
	__AGOS->aKss_lado_direito_= ((__AGOS->calc_ass-__AGOS->aKss_old_)/__AGOS->calc_tau_Kss);
} //fim

void __tree8__( Solveode *__AGOS){
	__AGOS->calc_iss = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+4.5200000000e+01)/5.7000000000e+00))));
	__AGOS->calc_tau_ti_s = (2.7000000000e+02+(1.0500000000e+03/(1.0000000000e+00+exp(((__AGOS->V_old_+4.5200000000e+01)/5.7000000000e+00)))));
	__AGOS->calc_tau_iur = (1.2000000000e+03-(1.7000000000e+02/(1.0000000000e+00+exp(((__AGOS->V_old_+4.5200000000e+01)/5.7000000000e+00)))));
	__AGOS->ito_s_lado_direito_= ((__AGOS->calc_iss-__AGOS->ito_s_old_)/__AGOS->calc_tau_ti_s);
	__AGOS->iur_lado_direito_= ((__AGOS->calc_iss-__AGOS->iur_old_)/__AGOS->calc_tau_iur);
} //fim

void __tree9__( Solveode *__AGOS){
	__AGOS->calc_alpha_n = ((4.8133300000e-06*(__AGOS->V_old_+2.6500000000e+01))/(1.0000000000e+00-exp(((-1.2800000000e-01)*(__AGOS->V_old_+2.6500000000e+01)))));
	__AGOS->calc_beta_n = (9.5333300000e-05*exp(((-3.8000000000e-02)*(__AGOS->V_old_+2.6500000000e+01))));
	__AGOS->nKs_lado_direito_= ((__AGOS->calc_alpha_n*(1.0000000000e+00-__AGOS->nKs_old_))-(__AGOS->calc_beta_n*__AGOS->nKs_old_));
} //fim

void __tree10__( Solveode *__AGOS){
	__AGOS->calc_C_K0 = (1.0000000000e+00-(__AGOS->C_K1_old_+__AGOS->C_K2_old_+__AGOS->O_K_old_+__AGOS->I_K_old_));
	__AGOS->calc_alpha_a0 = (2.2348000000e-02*exp((1.1760000000e-02*__AGOS->V_old_)));
	__AGOS->calc_beta_a0 = (4.7002000000e-02*exp(((-6.3100000000e-02)*__AGOS->V_old_)));
	__AGOS->C_K1_lado_direito_= (((__AGOS->calc_alpha_a0*__AGOS->calc_C_K0)+(__AGOS->kb*__AGOS->C_K2_old_))-((__AGOS->calc_beta_a0*__AGOS->C_K1_old_)+(__AGOS->kf*__AGOS->C_K1_old_)));
} //fim

void __tree11__( Solveode *__AGOS){
	__AGOS->calc_alpha_a1 = (1.3733000000e-02*exp((3.8198000000e-02*__AGOS->V_old_)));
	__AGOS->calc_beta_a1 = (6.8900000000e-05*exp(((-4.1780000000e-02)*__AGOS->V_old_)));
	__AGOS->calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (9.0821000000e-02*exp((2.3391000000e-02*(__AGOS->V_old_+5.0000000000e+00))));
	__AGOS->calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (6.4970000000e-03*exp(((-3.2680000000e-02)*(__AGOS->V_old_+5.0000000000e+00))));
	__AGOS->C_K2_lado_direito_= (((__AGOS->kf*__AGOS->C_K1_old_)+(__AGOS->calc_beta_a1*__AGOS->O_K_old_))-((__AGOS->kb*__AGOS->C_K2_old_)+(__AGOS->calc_alpha_a1*__AGOS->C_K2_old_)));
	__AGOS->O_K_lado_direito_= (((__AGOS->calc_alpha_a1*__AGOS->C_K2_old_)+(__AGOS->calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__AGOS->I_K_old_))-((__AGOS->calc_beta_a1*__AGOS->O_K_old_)+(__AGOS->calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__AGOS->O_K_old_)));
	__AGOS->I_K_lado_direito_= ((__AGOS->calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__AGOS->O_K_old_)-(__AGOS->calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__AGOS->I_K_old_));
} //fim

void __tree12__( Solveode *__AGOS){
	__AGOS->LTRPN_Ca_lado_direito_= ((__AGOS->k_plus_ltrpn*__AGOS->Cai_old_*(__AGOS->LTRPN_tot-__AGOS->LTRPN_Ca_old_))-(__AGOS->k_minus_ltrpn*__AGOS->LTRPN_Ca_old_));
} //fim

void __tree13__( Solveode *__AGOS){
	__AGOS->HTRPN_Ca_lado_direito_= ((__AGOS->k_plus_htrpn*__AGOS->Cai_old_*(__AGOS->HTRPN_tot-__AGOS->HTRPN_Ca_old_))-(__AGOS->k_minus_htrpn*__AGOS->HTRPN_Ca_old_));
} //fim

void __tree14__( Solveode *__AGOS){
	__AGOS->P_O2_lado_direito_= ((__AGOS->k_plus_b*pow(__AGOS->Cass_old_,__AGOS->m)*__AGOS->P_O1_old_)-(__AGOS->k_minus_b*__AGOS->P_O2_old_));
} //fim

void __tree15__( Solveode *__AGOS){
	__AGOS->P_C2_lado_direito_= ((__AGOS->k_plus_c*__AGOS->P_O1_old_)-(__AGOS->k_minus_c*__AGOS->P_C2_old_));
} //fim

void __tree16__( Solveode *__AGOS){
	__AGOS->iKss_lado_direito_= 0.0000000000e+00;
} //fim

void __AGOS_EQUATIONS__( Solveode *__AGOS){
		const double time_new = __AGOS->time_new;
		const double time = 0.0000000000e+00;
		const double stim_amplitude = -8.0000000000e+01;
		const double stim_start = 2.0000000000e+01;
		const double stim_end = 1.0000000000e+05;
		const double stim_period = 7.1430000000e+01;
		const double stim_duration = 1.5000000000e+00;
		const double Acap = 1.5340000000e-04;
		const double Cm = 1.0000000000e+00;
		const double Vmyo = 2.5840000000e-05;
		const double F = 9.6500000000e+01;
		const double VJSR = 1.2000000000e-07;
		const double Vss = 1.4850000000e-09;
		const double VNSR = 2.0980000000e-06;
		const double CMDN_tot = 5.0000000000e+01;
		const double Km_CMDN = 2.3800000000e-01;
		const double CSQN_tot = 1.5000000000e+04;
		const double Km_CSQN = 8.0000000000e+02;
		const double v1 = 4.5000000000e+00;
		const double tau_tr = 2.0000000000e+01;
		const double tau_xfer = 8.0000000000e+00;
		const double v2 = 1.7400000000e-05;
		const double v3 = 4.5000000000e-01;
		const double Km_up = 5.0000000000e-01;
		const double k_plus_htrpn = 2.3700000000e-03;
		const double HTRPN_tot = 1.4000000000e+02;
		const double k_plus_ltrpn = 3.2700000000e-02;
		const double LTRPN_tot = 7.0000000000e+01;
		const double k_minus_htrpn = 3.2000000000e-05;
		const double k_minus_ltrpn = 1.9600000000e-02;
		const double i_CaL_max = 7.0000000000e+00;
		const double k_plus_a = 6.0750000000e-03;
		const double n = 4.0000000000e+00;
		const double k_minus_b = 9.6500000000e-01;
		const double k_minus_c = 8.0000000000e-04;
		const double k_minus_a = 7.1250000000e-02;
		const double k_plus_b = 4.0500000000e-03;
		const double m = 3.0000000000e+00;
		const double k_plus_c = 9.0000000000e-03;
		const double g_CaL = 1.7290000000e-01;
		const double E_CaL = 6.3000000000e+01;
		const double Kpcb = 5.0000000000e-04;
		const double Kpc_max = 2.3324000000e-01;
		const double Kpc_half = 2.0000000000e+01;
		const double i_pCa_max = 1.0000000000e+00;
		const double Km_pCa = 5.0000000000e-01;
		const double k_NaCa = 2.9280000000e+02;
		const double K_mNa = 8.7500000000e+04;
		const double Nao = 1.4000000000e+05;
		const double K_mCa = 1.3800000000e+03;
		const double Cao = 1.8000000000e+03;
		const double k_sat = 1.0000000000e-01;
		const double eta = 3.5000000000e-01;
		const double R = 8.3140000000e+00;
		const double T = 2.9800000000e+02;
		const double g_Cab = 3.6700000000e-04;
		const double g_Na = 1.3000000000e+01;
		const double Ko = 5.4000000000e+03;
		const double g_Nab = 2.6000000000e-03;
		const double g_Kto_f = 4.0670000000e-01;
		const double g_Kto_s = 0.0000000000e+00;
		const double g_Ks = 5.7500000000e-03;
		const double g_Kur = 1.6000000000e-01;
		const double g_Kss = 5.0000000000e-02;
		const double g_Kr = 7.8000000000e-02;
		const double kf = 2.3761000000e-02;
		const double kb = 3.6778000000e-02;
		const double i_NaK_max = 8.8000000000e-01;
		const double Km_Nai = 2.1000000000e+04;
		const double Km_Ko = 1.5000000000e+03;
		const double g_ClCa = 1.0000000000e+01;
		const double Km_Cl = 1.0000000000e+01;
		const double E_Cl = -4.0000000000e+01;
		const double V_old_= __AGOS->V_old_;
		const double Cai_old_= __AGOS->Cai_old_;
		const double Cass_old_= __AGOS->Cass_old_;
		const double CaJSR_old_= __AGOS->CaJSR_old_;
		const double CaNSR_old_= __AGOS->CaNSR_old_;
		const double P_RyR_old_= __AGOS->P_RyR_old_;
		const double LTRPN_Ca_old_= __AGOS->LTRPN_Ca_old_;
		const double HTRPN_Ca_old_= __AGOS->HTRPN_Ca_old_;
		const double P_O1_old_= __AGOS->P_O1_old_;
		const double P_O2_old_= __AGOS->P_O2_old_;
		const double P_C2_old_= __AGOS->P_C2_old_;
		const double O_old_= __AGOS->O_old_;
		const double C2_old_= __AGOS->C2_old_;
		const double C3_old_= __AGOS->C3_old_;
		const double C4_old_= __AGOS->C4_old_;
		const double I1_old_= __AGOS->I1_old_;
		const double I2_old_= __AGOS->I2_old_;
		const double I3_old_= __AGOS->I3_old_;
		const double Nai_old_= __AGOS->Nai_old_;
		const double C_Na2_old_= __AGOS->C_Na2_old_;
		const double C_Na1_old_= __AGOS->C_Na1_old_;
		const double O_Na_old_= __AGOS->O_Na_old_;
		const double IF_Na_old_= __AGOS->IF_Na_old_;
		const double I1_Na_old_= __AGOS->I1_Na_old_;
		const double I2_Na_old_= __AGOS->I2_Na_old_;
		const double IC_Na2_old_= __AGOS->IC_Na2_old_;
		const double IC_Na3_old_= __AGOS->IC_Na3_old_;
		const double Ki_old_= __AGOS->Ki_old_;
		const double ato_f_old_= __AGOS->ato_f_old_;
		const double ito_f_old_= __AGOS->ito_f_old_;
		const double ato_s_old_= __AGOS->ato_s_old_;
		const double ito_s_old_= __AGOS->ito_s_old_;
		const double nKs_old_= __AGOS->nKs_old_;
		const double aur_old_= __AGOS->aur_old_;
		const double iur_old_= __AGOS->iur_old_;
		const double aKss_old_= __AGOS->aKss_old_;
		const double iKss_old_= __AGOS->iKss_old_;
		const double C_K2_old_= __AGOS->C_K2_old_;
		const double C_K1_old_= __AGOS->C_K1_old_;
		const double O_K_old_= __AGOS->O_K_old_;
		const double I_K_old_= __AGOS->I_K_old_;
	const double calc_i_stim = (((time_new>=stim_start)&&(time_new<=stim_end)&&(((time_new-stim_start)-(floor(((time_new-stim_start)/stim_period))*stim_period))<=stim_duration)))
?(stim_amplitude)
:(0.0000000000e+00);
	const double calc_Bi = pow((1.0000000000e+00+((CMDN_tot*Km_CMDN)/pow((Km_CMDN+Cai_old_),2.0000000000e+00))),(-1.0000000000e+00));
	const double calc_Bss = pow((1.0000000000e+00+((CMDN_tot*Km_CMDN)/pow((Km_CMDN+Cass_old_),2.0000000000e+00))),(-1.0000000000e+00));
	const double calc_BJSR = pow((1.0000000000e+00+((CSQN_tot*Km_CSQN)/pow((Km_CSQN+CaJSR_old_),2.0000000000e+00))),(-1.0000000000e+00));
	const double calc_J_rel = (v1*(P_O1_old_+P_O2_old_)*(CaJSR_old_-Cass_old_)*P_RyR_old_);
	const double calc_J_tr = ((CaNSR_old_-CaJSR_old_)/tau_tr);
	const double calc_J_xfer = ((Cass_old_-Cai_old_)/tau_xfer);
	const double calc_J_leak = (v2*(CaNSR_old_-Cai_old_));
	const double calc_J_up = ((v3*pow(Cai_old_,2.0000000000e+00))/(pow(Km_up,2.0000000000e+00)+pow(Cai_old_,2.0000000000e+00)));
	const double calc_J_trpn = (((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))+(k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_)))-((k_minus_htrpn*HTRPN_Ca_old_)+(k_minus_ltrpn*LTRPN_Ca_old_)));
	const double calc_i_CaL = (g_CaL*O_old_*(V_old_-E_CaL));
	const double calc_i_pCa = ((i_pCa_max*pow(Cai_old_,2.0000000000e+00))/(pow(Km_pCa,2.0000000000e+00)+pow(Cai_old_,2.0000000000e+00)));
	const double calc_i_NaCa = (((((((k_NaCa*1.0000000000e+00)/(pow(K_mNa,3.0000000000e+00)+pow(Nao,3.0000000000e+00)))*1.0000000000e+00)/(K_mCa+Cao))*1.0000000000e+00)/(1.0000000000e+00+(k_sat*exp((((eta-1.0000000000e+00)*V_old_*F)/(R*T))))))*((exp(((eta*V_old_*F)/(R*T)))*pow(Nai_old_,3.0000000000e+00)*Cao)-(exp((((eta-1.0000000000e+00)*V_old_*F)/(R*T)))*pow(Nao,3.0000000000e+00)*Cai_old_)));
	const double calc_E_CaN = (((R*T)/(2.0000000000e+00*F))*log((Cao/Cai_old_)));
	const double calc_E_Na = (((R*T)/F)*log((((9.0000000000e-01*Nao)+(1.0000000000e-01*Ko))/((9.0000000000e-01*Nai_old_)+(1.0000000000e-01*Ki_old_)))));
	const double calc_E_K = (((R*T)/F)*log((Ko/Ki_old_)));
	const double calc_i_Kr = (g_Kr*O_K_old_*(V_old_-(((R*T)/F)*log((((9.8000000000e-01*Ko)+(2.0000000000e-02*Nao))/((9.8000000000e-01*Ki_old_)+(2.0000000000e-02*Nai_old_)))))));
	const double calc_sigma = ((1.0000000000e+00/7.0000000000e+00)*(exp((Nao/6.7300000000e+04))-1.0000000000e+00));
	const double calc_O_ClCa = (2.0000000000e-01/(1.0000000000e+00+exp(((-(V_old_-4.6700000000e+01))/7.8000000000e+00))));
	const double calc_i_Nab = (g_Nab*(V_old_-calc_E_Na));
	const double calc_i_Kto_s = (g_Kto_s*ato_s_old_*ito_s_old_*(V_old_-calc_E_K));
	const double calc_i_K1 = ((((2.9380000000e-01*Ko)/(Ko+2.1000000000e+02))*(V_old_-calc_E_K))/(1.0000000000e+00+exp((8.9600000000e-02*(V_old_-calc_E_K)))));
	const double calc_i_Ks = (g_Ks*pow(nKs_old_,2.0000000000e+00)*(V_old_-calc_E_K));
	const double calc_i_Kur = (g_Kur*aur_old_*iur_old_*(V_old_-calc_E_K));
	const double calc_i_Kss = (g_Kss*aKss_old_*iKss_old_*(V_old_-calc_E_K));
	const double calc_i_Cab = (g_Cab*(V_old_-calc_E_CaN));
	const double calc_i_Na = (g_Na*O_Na_old_*(V_old_-calc_E_Na));
	const double calc_i_Kto_f = (g_Kto_f*pow(ato_f_old_,3.0000000000e+00)*ito_f_old_*(V_old_-calc_E_K));
	const double calc_f_NaK = (1.0000000000e+00/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*V_old_*F)/(R*T))))+(3.6500000000e-02*calc_sigma*exp((((-V_old_)*F)/(R*T))))));
	const double calc_i_ClCa = (((g_ClCa*calc_O_ClCa*Cai_old_)/(Cai_old_+Km_Cl))*(V_old_-E_Cl));
	const double calc_i_NaK = ((((i_NaK_max*calc_f_NaK*1.0000000000e+00)/(1.0000000000e+00+pow((Km_Nai/Nai_old_),1.5000000000e+00)))*Ko)/(Ko+Km_Ko));
	__AGOS->V_lado_direito_= (-(calc_i_CaL+calc_i_pCa+calc_i_NaCa+calc_i_Cab+calc_i_Na+calc_i_Nab+calc_i_NaK+calc_i_Kto_f+calc_i_Kto_s+calc_i_K1+calc_i_Ks+calc_i_Kur+calc_i_Kss+calc_i_Kr+calc_i_ClCa+calc_i_stim));
	__AGOS->Cai_lado_direito_= (calc_Bi*((calc_J_leak+calc_J_xfer)-(calc_J_up+calc_J_trpn+((((calc_i_Cab+calc_i_pCa)-(2.0000000000e+00*calc_i_NaCa))*Acap*Cm)/(2.0000000000e+00*Vmyo*F)))));
	__AGOS->Cass_lado_direito_= (calc_Bss*(((calc_J_rel*VJSR)/Vss)-(((calc_J_xfer*Vmyo)/Vss)+((calc_i_CaL*Acap*Cm)/(2.0000000000e+00*Vss*F)))));
	__AGOS->CaJSR_lado_direito_= (calc_BJSR*(calc_J_tr-calc_J_rel));
	__AGOS->CaNSR_lado_direito_= ((((calc_J_up-calc_J_leak)*Vmyo)/VNSR)-((calc_J_tr*VJSR)/VNSR));
	__AGOS->P_RyR_lado_direito_= (((-4.0000000000e-02)*P_RyR_old_)-(((1.0000000000e-01*calc_i_CaL)/i_CaL_max)*exp(((-pow((V_old_-5.0000000000e+00),2.0000000000e+00))/6.4800000000e+02))));
	__AGOS->Nai_lado_direito_= (((-(calc_i_Na+calc_i_Nab+(3.0000000000e+00*calc_i_NaK)+(3.0000000000e+00*calc_i_NaCa)))*Acap*Cm)/(Vmyo*F));
	__AGOS->Ki_lado_direito_= (((-((calc_i_Kto_f+calc_i_Kto_s+calc_i_K1+calc_i_Ks+calc_i_Kss+calc_i_Kur+calc_i_Kr)-(2.0000000000e+00*calc_i_NaK)))*Acap*Cm)/(Vmyo*F));
	const double calc_P_C1 = (1.0000000000e+00-(P_C2_old_+P_O1_old_+P_O2_old_));
	__AGOS->P_O1_lado_direito_= (((k_plus_a*pow(Cass_old_,n)*calc_P_C1)+(k_minus_b*P_O2_old_)+(k_minus_c*P_C2_old_))-((k_minus_a*P_O1_old_)+(k_plus_b*pow(Cass_old_,m)*P_O1_old_)+(k_plus_c*P_O1_old_)));
	const double calc_C1 = (1.0000000000e+00-(O_old_+C2_old_+C2_old_+C3_old_+C4_old_+I1_old_+I2_old_+I3_old_));
	const double calc_alpha = ((4.0000000000e-01*exp(((V_old_+1.2000000000e+01)/1.0000000000e+01))*((1.0000000000e+00+(7.0000000000e-01*exp(((-pow((V_old_+4.0000000000e+01),2.0000000000e+00))/1.0000000000e+01))))-(7.5000000000e-01*exp(((-pow((V_old_+2.0000000000e+01),2.0000000000e+00))/4.0000000000e+02)))))/(1.0000000000e+00+(1.2000000000e-01*exp(((V_old_+1.2000000000e+01)/1.0000000000e+01)))));
	const double calc_beta = (5.0000000000e-02*exp(((-(V_old_+1.2000000000e+01))/1.3000000000e+01)));
	const double calc_gamma = ((Kpc_max*Cass_old_)/(Kpc_half+Cass_old_));
	const double calc_Kpcf = (1.3000000000e+01*(1.0000000000e+00-exp(((-pow((V_old_+1.4500000000e+01),2.0000000000e+00))/1.0000000000e+02))));
	__AGOS->O_lado_direito_= (((calc_alpha*C4_old_)+(Kpcb*I1_old_)+(1.0000000000e-03*((calc_alpha*I2_old_)-(calc_Kpcf*O_old_))))-((4.0000000000e+00*calc_beta*O_old_)+(calc_gamma*O_old_)));
	__AGOS->C2_lado_direito_= (((4.0000000000e+00*calc_alpha*calc_C1)+(2.0000000000e+00*calc_beta*C3_old_))-((calc_beta*C2_old_)+(3.0000000000e+00*calc_alpha*C2_old_)));
	__AGOS->C3_lado_direito_= (((3.0000000000e+00*calc_alpha*C2_old_)+(3.0000000000e+00*calc_beta*C4_old_))-((2.0000000000e+00*calc_beta*C3_old_)+(2.0000000000e+00*calc_alpha*C3_old_)));
	__AGOS->C4_lado_direito_= (((2.0000000000e+00*calc_alpha*C3_old_)+(4.0000000000e+00*calc_beta*O_old_)+(1.0000000000e-02*((4.0000000000e+00*Kpcb*calc_beta*I1_old_)-(calc_alpha*calc_gamma*C4_old_)))+(2.0000000000e-03*((4.0000000000e+00*calc_beta*I2_old_)-(calc_Kpcf*C4_old_)))+(4.0000000000e+00*calc_beta*Kpcb*I3_old_))-((3.0000000000e+00*calc_beta*C4_old_)+(calc_alpha*C4_old_)+(1.0000000000e+00*calc_gamma*calc_Kpcf*C4_old_)));
	__AGOS->I1_lado_direito_= (((calc_gamma*O_old_)+(1.0000000000e-03*((calc_alpha*I3_old_)-(calc_Kpcf*I1_old_)))+(1.0000000000e-02*((calc_alpha*calc_gamma*C4_old_)-(4.0000000000e+00*calc_beta*calc_Kpcf*I1_old_))))-(Kpcb*I1_old_));
	__AGOS->I2_lado_direito_= (((1.0000000000e-03*((calc_Kpcf*O_old_)-(calc_alpha*I2_old_)))+(Kpcb*I3_old_)+(2.0000000000e-03*((calc_Kpcf*C4_old_)-(4.0000000000e+00*calc_beta*I2_old_))))-(calc_gamma*I2_old_));
	__AGOS->I3_lado_direito_= (((1.0000000000e-03*((calc_Kpcf*I1_old_)-(calc_alpha*I3_old_)))+(calc_gamma*I2_old_)+(1.0000000000e+00*calc_gamma*calc_Kpcf*C4_old_))-((4.0000000000e+00*calc_beta*Kpcb*I3_old_)+(Kpcb*I3_old_)));
	const double calc_C_Na3 = (1.0000000000e+00-(O_Na_old_+C_Na1_old_+C_Na2_old_+IF_Na_old_+I1_Na_old_+I2_Na_old_+IC_Na2_old_+IC_Na3_old_));
	const double calc_alpha_Na11 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(V_old_+2.5000000000e+00))/1.7000000000e+01)))+(2.0000000000e-01*exp(((-(V_old_+2.5000000000e+00))/1.5000000000e+02)))));
	const double calc_alpha_Na12 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(V_old_+2.5000000000e+00))/1.5000000000e+01)))+(2.3000000000e-01*exp(((-(V_old_+2.5000000000e+00))/1.5000000000e+02)))));
	const double calc_alpha_Na13 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(V_old_+2.5000000000e+00))/1.2000000000e+01)))+(2.5000000000e-01*exp(((-(V_old_+2.5000000000e+00))/1.5000000000e+02)))));
	const double calc_beta_Na11 = (1.9170000000e-01*exp(((-(V_old_+2.5000000000e+00))/2.0300000000e+01)));
	const double calc_beta_Na12 = (2.0000000000e-01*exp(((-(V_old_-2.5000000000e+00))/2.0300000000e+01)));
	const double calc_beta_Na13 = (2.2000000000e-01*exp(((-(V_old_-7.5000000000e+00))/2.0300000000e+01)));
	const double calc_alpha_Na3 = (7.0000000000e-07*exp(((-(V_old_+7.0000000000e+00))/7.7000000000e+00)));
	const double calc_beta_Na3 = (8.4000000000e-03+(2.0000000000e-05*(V_old_+7.0000000000e+00)));
	const double calc_alpha_Na2 = (1.0000000000e+00/((1.8849500000e-01*exp(((-(V_old_+7.0000000000e+00))/1.6600000000e+01)))+3.9395600000e-01));
	const double calc_beta_Na2 = ((calc_alpha_Na13*calc_alpha_Na2*calc_alpha_Na3)/(calc_beta_Na13*calc_beta_Na3));
	const double calc_alpha_Na4 = (calc_alpha_Na2/1.0000000000e+03);
	const double calc_beta_Na4 = calc_alpha_Na3;
	const double calc_alpha_Na5 = (calc_alpha_Na2/9.5000000000e+04);
	const double calc_beta_Na5 = (calc_alpha_Na3/5.0000000000e+01);
	__AGOS->C_Na2_lado_direito_= (((calc_alpha_Na11*calc_C_Na3)+(calc_beta_Na12*C_Na1_old_)+(calc_alpha_Na3*IC_Na2_old_))-((calc_beta_Na11*C_Na2_old_)+(calc_alpha_Na12*C_Na2_old_)+(calc_beta_Na3*C_Na2_old_)));
	__AGOS->C_Na1_lado_direito_= (((calc_alpha_Na12*C_Na2_old_)+(calc_beta_Na13*O_Na_old_)+(calc_alpha_Na3*IF_Na_old_))-((calc_beta_Na12*C_Na1_old_)+(calc_alpha_Na13*C_Na1_old_)+(calc_beta_Na3*C_Na1_old_)));
	__AGOS->O_Na_lado_direito_= (((calc_alpha_Na13*C_Na1_old_)+(calc_beta_Na2*IF_Na_old_))-((calc_beta_Na13*O_Na_old_)+(calc_alpha_Na2*O_Na_old_)));
	__AGOS->IF_Na_lado_direito_= (((calc_alpha_Na2*O_Na_old_)+(calc_beta_Na3*C_Na1_old_)+(calc_beta_Na4*I1_Na_old_)+(calc_alpha_Na12*IC_Na2_old_))-((calc_beta_Na2*IF_Na_old_)+(calc_alpha_Na3*IF_Na_old_)+(calc_alpha_Na4*IF_Na_old_)+(calc_beta_Na12*IF_Na_old_)));
	__AGOS->I1_Na_lado_direito_= (((calc_alpha_Na4*IF_Na_old_)+(calc_beta_Na5*I2_Na_old_))-((calc_beta_Na4*I1_Na_old_)+(calc_alpha_Na5*I1_Na_old_)));
	__AGOS->I2_Na_lado_direito_= ((calc_alpha_Na5*I1_Na_old_)-(calc_beta_Na5*I2_Na_old_));
	__AGOS->IC_Na2_lado_direito_= (((calc_alpha_Na11*IC_Na3_old_)+(calc_beta_Na12*IF_Na_old_)+(calc_beta_Na3*IC_Na2_old_))-((calc_beta_Na11*IC_Na2_old_)+(calc_alpha_Na12*IC_Na2_old_)+(calc_alpha_Na3*IC_Na2_old_)));
	__AGOS->IC_Na3_lado_direito_= (((calc_beta_Na11*IC_Na2_old_)+(calc_beta_Na3*calc_C_Na3))-((calc_alpha_Na11*IC_Na3_old_)+(calc_alpha_Na3*IC_Na3_old_)));
	const double calc_alpha_a = (1.8064000000e-01*exp((3.5770000000e-02*(V_old_+3.0000000000e+01))));
	const double calc_beta_a = (3.9560000000e-01*exp(((-6.2370000000e-02)*(V_old_+3.0000000000e+01))));
	__AGOS->ato_f_lado_direito_= ((calc_alpha_a*(1.0000000000e+00-ato_f_old_))-(calc_beta_a*ato_f_old_));
	const double calc_alpha_i = ((1.5200000000e-04*exp(((-(V_old_+1.3500000000e+01))/7.0000000000e+00)))/((6.7083000000e-03*exp(((-(V_old_+3.3500000000e+01))/7.0000000000e+00)))+1.0000000000e+00));
	const double calc_beta_i = ((9.5000000000e-04*exp(((V_old_+3.3500000000e+01)/7.0000000000e+00)))/((5.1335000000e-02*exp(((V_old_+3.3500000000e+01)/7.0000000000e+00)))+1.0000000000e+00));
	__AGOS->ito_f_lado_direito_= ((calc_alpha_i*(1.0000000000e+00-ito_f_old_))-(calc_beta_i*ito_f_old_));
	const double calc_ass = (1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+2.2500000000e+01))/7.7000000000e+00))));
	const double calc_tau_ta_s = ((4.9300000000e-01*exp(((-6.2900000000e-02)*V_old_)))+2.0580000000e+00);
	const double calc_tau_aur = ((4.9300000000e-01*exp(((-6.2900000000e-02)*V_old_)))+2.0580000000e+00);
	const double calc_tau_Kss = ((3.9300000000e+01*exp(((-8.6200000000e-02)*V_old_)))+1.3170000000e+01);
	__AGOS->ato_s_lado_direito_= ((calc_ass-ato_s_old_)/calc_tau_ta_s);
	__AGOS->aur_lado_direito_= ((calc_ass-aur_old_)/calc_tau_aur);
	__AGOS->aKss_lado_direito_= ((calc_ass-aKss_old_)/calc_tau_Kss);
	const double calc_iss = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+4.5200000000e+01)/5.7000000000e+00))));
	const double calc_tau_ti_s = (2.7000000000e+02+(1.0500000000e+03/(1.0000000000e+00+exp(((V_old_+4.5200000000e+01)/5.7000000000e+00)))));
	const double calc_tau_iur = (1.2000000000e+03-(1.7000000000e+02/(1.0000000000e+00+exp(((V_old_+4.5200000000e+01)/5.7000000000e+00)))));
	__AGOS->ito_s_lado_direito_= ((calc_iss-ito_s_old_)/calc_tau_ti_s);
	__AGOS->iur_lado_direito_= ((calc_iss-iur_old_)/calc_tau_iur);
	const double calc_alpha_n = ((4.8133300000e-06*(V_old_+2.6500000000e+01))/(1.0000000000e+00-exp(((-1.2800000000e-01)*(V_old_+2.6500000000e+01)))));
	const double calc_beta_n = (9.5333300000e-05*exp(((-3.8000000000e-02)*(V_old_+2.6500000000e+01))));
	__AGOS->nKs_lado_direito_= ((calc_alpha_n*(1.0000000000e+00-nKs_old_))-(calc_beta_n*nKs_old_));
	const double calc_C_K0 = (1.0000000000e+00-(C_K1_old_+C_K2_old_+O_K_old_+I_K_old_));
	const double calc_alpha_a0 = (2.2348000000e-02*exp((1.1760000000e-02*V_old_)));
	const double calc_beta_a0 = (4.7002000000e-02*exp(((-6.3100000000e-02)*V_old_)));
	__AGOS->C_K1_lado_direito_= (((calc_alpha_a0*calc_C_K0)+(kb*C_K2_old_))-((calc_beta_a0*C_K1_old_)+(kf*C_K1_old_)));
	const double calc_alpha_a1 = (1.3733000000e-02*exp((3.8198000000e-02*V_old_)));
	const double calc_beta_a1 = (6.8900000000e-05*exp(((-4.1780000000e-02)*V_old_)));
	const double calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (9.0821000000e-02*exp((2.3391000000e-02*(V_old_+5.0000000000e+00))));
	const double calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (6.4970000000e-03*exp(((-3.2680000000e-02)*(V_old_+5.0000000000e+00))));
	__AGOS->C_K2_lado_direito_= (((kf*C_K1_old_)+(calc_beta_a1*O_K_old_))-((kb*C_K2_old_)+(calc_alpha_a1*C_K2_old_)));
	__AGOS->O_K_lado_direito_= (((calc_alpha_a1*C_K2_old_)+(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*I_K_old_))-((calc_beta_a1*O_K_old_)+(calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*O_K_old_)));
	__AGOS->I_K_lado_direito_= ((calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*O_K_old_)-(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*I_K_old_));
	__AGOS->LTRPN_Ca_lado_direito_= ((k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_))-(k_minus_ltrpn*LTRPN_Ca_old_));
	__AGOS->HTRPN_Ca_lado_direito_= ((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))-(k_minus_htrpn*HTRPN_Ca_old_));
	__AGOS->P_O2_lado_direito_= ((k_plus_b*pow(Cass_old_,m)*P_O1_old_)-(k_minus_b*P_O2_old_));
	__AGOS->P_C2_lado_direito_= ((k_plus_c*P_O1_old_)-(k_minus_c*P_C2_old_));
	__AGOS->iKss_lado_direito_= 0.0000000000e+00;
} //fim
typedef struct str__rightHandSideFunction{
	void (*function)(Solveode*);
}typ_rightHandSideFunction;
typ_rightHandSideFunction rightHandSideFunction;
typ_rightHandSideFunction forest[16];
	#include "simple.hpp"
	#include "pe.hpp"
	#include "lut.hpp"
	#include "pe_lut.hpp"

	//Constructor, initializes all variables with 0.0 or default CellML initial values
	Solveode::Solveode(double abs, double rel, int equationType)
	{
		if(equationType==eqAgos){
			rightHandSideFunction.function = __AGOS_EQUATIONS__;
		}else if(equationType==eqPyCmlSimples){
			rightHandSideFunction.function = __equation__;
		}else if(equationType==eqPyCmlPE){
			rightHandSideFunction.function = __equation__pe_;
		}else if(equationType==eqPyCmlLUT){
			rightHandSideFunction.function = __equation__lut_;
		}else if(equationType==eqPyCmlPE_LUT){
			rightHandSideFunction.function = __equation__pe_lut_;
		}
		forest[0].function = __tree1__;
		forest[1].function = __tree2__;
		forest[2].function = __tree3__;
		forest[3].function = __tree4__;
		forest[4].function = __tree5__;
		forest[5].function = __tree6__;
		forest[6].function = __tree7__;
		forest[7].function = __tree8__;
		forest[8].function = __tree9__;
		forest[9].function = __tree10__;
		forest[10].function = __tree11__;
		forest[11].function = __tree12__;
		forest[12].function = __tree13__;
		forest[13].function = __tree14__;
		forest[14].function = __tree15__;
		forest[15].function = __tree16__;
		time = 0.00000000e+00;
		stim_amplitude = -8.00000000e+01;
		stim_start = 2.00000000e+01;
		stim_end = 1.00000000e+05;
		stim_period = 7.14300000e+01;
		stim_duration = 1.50000000e+00;
		Acap = 1.53400000e-04;
		Cm = 1.00000000e+00;
		Vmyo = 2.58400000e-05;
		F = 9.65000000e+01;
		VJSR = 1.20000000e-07;
		Vss = 1.48500000e-09;
		VNSR = 2.09800000e-06;
		CMDN_tot = 5.00000000e+01;
		Km_CMDN = 2.38000000e-01;
		CSQN_tot = 1.50000000e+04;
		Km_CSQN = 8.00000000e+02;
		v1 = 4.50000000e+00;
		tau_tr = 2.00000000e+01;
		tau_xfer = 8.00000000e+00;
		v2 = 1.74000000e-05;
		v3 = 4.50000000e-01;
		Km_up = 5.00000000e-01;
		k_plus_htrpn = 2.37000000e-03;
		HTRPN_tot = 1.40000000e+02;
		k_plus_ltrpn = 3.27000000e-02;
		LTRPN_tot = 7.00000000e+01;
		k_minus_htrpn = 3.20000000e-05;
		k_minus_ltrpn = 1.96000000e-02;
		i_CaL_max = 7.00000000e+00;
		k_plus_a = 6.07500000e-03;
		n = 4.00000000e+00;
		k_minus_b = 9.65000000e-01;
		k_minus_c = 8.00000000e-04;
		k_minus_a = 7.12500000e-02;
		k_plus_b = 4.05000000e-03;
		m = 3.00000000e+00;
		k_plus_c = 9.00000000e-03;
		g_CaL = 1.72900000e-01;
		E_CaL = 6.30000000e+01;
		Kpcb = 5.00000000e-04;
		Kpc_max = 2.33240000e-01;
		Kpc_half = 2.00000000e+01;
		i_pCa_max = 1.00000000e+00;
		Km_pCa = 5.00000000e-01;
		k_NaCa = 2.92800000e+02;
		K_mNa = 8.75000000e+04;
		Nao = 1.40000000e+05;
		K_mCa = 1.38000000e+03;
		Cao = 1.80000000e+03;
		k_sat = 1.00000000e-01;
		eta = 3.50000000e-01;
		R = 8.31400000e+00;
		T = 2.98000000e+02;
		g_Cab = 3.67000000e-04;
		g_Na = 1.30000000e+01;
		Ko = 5.40000000e+03;
		g_Nab = 2.60000000e-03;
		g_Kto_f = 4.06700000e-01;
		g_Kto_s = 0.00000000e+00;
		g_Ks = 5.75000000e-03;
		g_Kur = 1.60000000e-01;
		g_Kss = 5.00000000e-02;
		g_Kr = 7.80000000e-02;
		kf = 2.37610000e-02;
		kb = 3.67780000e-02;
		i_NaK_max = 8.80000000e-01;
		Km_Nai = 2.10000000e+04;
		Km_Ko = 1.50000000e+03;
		g_ClCa = 1.00000000e+01;
		Km_Cl = 1.00000000e+01;
		E_Cl = -4.00000000e+01;
		dtime = 0.0; time_vec__ = NULL;
		V = NULL;
		V_ini_ = -8.24202000e+01;
		Cai = NULL;
		Cai_ini_ = 1.15001000e-01;
		Cass = NULL;
		Cass_ini_ = 1.15001000e-01;
		CaJSR = NULL;
		CaJSR_ini_ = 1.29950000e+03;
		CaNSR = NULL;
		CaNSR_ini_ = 1.29950000e+03;
		P_RyR = NULL;
		P_RyR_ini_ = 0.00000000e+00;
		LTRPN_Ca = NULL;
		LTRPN_Ca_ini_ = 1.12684000e+01;
		HTRPN_Ca = NULL;
		HTRPN_Ca_ini_ = 1.25290000e+02;
		P_O1 = NULL;
		P_O1_ini_ = 1.49102000e-05;
		P_O2 = NULL;
		P_O2_ini_ = 9.51726000e-11;
		P_C2 = NULL;
		P_C2_ini_ = 1.67740000e-04;
		O = NULL;
		O_ini_ = 9.30308000e-19;
		C2 = NULL;
		C2_ini_ = 1.24216000e-04;
		C3 = NULL;
		C3_ini_ = 5.78679000e-09;
		C4 = NULL;
		C4_ini_ = 1.19816000e-13;
		I1 = NULL;
		I1_ini_ = 4.97923000e-19;
		I2 = NULL;
		I2_ini_ = 3.45847000e-14;
		I3 = NULL;
		I3_ini_ = 1.85106000e-14;
		Nai = NULL;
		Nai_ini_ = 1.42371000e+04;
		C_Na2 = NULL;
		C_Na2_ini_ = 2.07520000e-02;
		C_Na1 = NULL;
		C_Na1_ini_ = 2.79132000e-04;
		O_Na = NULL;
		O_Na_ini_ = 7.13483000e-07;
		IF_Na = NULL;
		IF_Na_ini_ = 1.53176000e-04;
		I1_Na = NULL;
		I1_Na_ini_ = 6.73345000e-07;
		I2_Na = NULL;
		I2_Na_ini_ = 1.55787000e-09;
		IC_Na2 = NULL;
		IC_Na2_ini_ = 1.13879000e-02;
		IC_Na3 = NULL;
		IC_Na3_ini_ = 3.42780000e-01;
		Ki = NULL;
		Ki_ini_ = 1.43720000e+05;
		ato_f = NULL;
		ato_f_ini_ = 2.65563000e-03;
		ito_f = NULL;
		ito_f_ini_ = 9.99977000e-01;
		ato_s = NULL;
		ato_s_ini_ = 4.17069000e-04;
		ito_s = NULL;
		ito_s_ini_ = 9.98543000e-01;
		nKs = NULL;
		nKs_ini_ = 2.62753000e-04;
		aur = NULL;
		aur_ini_ = 4.17069000e-04;
		iur = NULL;
		iur_ini_ = 9.98543000e-01;
		aKss = NULL;
		aKss_ini_ = 4.17069000e-04;
		iKss = NULL;
		iKss_ini_ = 1.00000000e+00;
		C_K2 = NULL;
		C_K2_ini_ = 6.41229000e-04;
		C_K1 = NULL;
		C_K1_ini_ = 9.92513000e-04;
		O_K = NULL;
		O_K_ini_ = 1.75298000e-04;
		I_K = NULL;
		I_K_ini_ = 3.19129000e-05;
		abstol__ = abs;
		reltol__ = rel;
		it_countx = 0;
	}
	Solveode::~Solveode()
	{
		if(V != NULL) free(V);
		if(Cai != NULL) free(Cai);
		if(Cass != NULL) free(Cass);
		if(CaJSR != NULL) free(CaJSR);
		if(CaNSR != NULL) free(CaNSR);
		if(P_RyR != NULL) free(P_RyR);
		if(LTRPN_Ca != NULL) free(LTRPN_Ca);
		if(HTRPN_Ca != NULL) free(HTRPN_Ca);
		if(P_O1 != NULL) free(P_O1);
		if(P_O2 != NULL) free(P_O2);
		if(P_C2 != NULL) free(P_C2);
		if(O != NULL) free(O);
		if(C2 != NULL) free(C2);
		if(C3 != NULL) free(C3);
		if(C4 != NULL) free(C4);
		if(I1 != NULL) free(I1);
		if(I2 != NULL) free(I2);
		if(I3 != NULL) free(I3);
		if(Nai != NULL) free(Nai);
		if(C_Na2 != NULL) free(C_Na2);
		if(C_Na1 != NULL) free(C_Na1);
		if(O_Na != NULL) free(O_Na);
		if(IF_Na != NULL) free(IF_Na);
		if(I1_Na != NULL) free(I1_Na);
		if(I2_Na != NULL) free(I2_Na);
		if(IC_Na2 != NULL) free(IC_Na2);
		if(IC_Na3 != NULL) free(IC_Na3);
		if(Ki != NULL) free(Ki);
		if(ato_f != NULL) free(ato_f);
		if(ito_f != NULL) free(ito_f);
		if(ato_s != NULL) free(ato_s);
		if(ito_s != NULL) free(ito_s);
		if(nKs != NULL) free(nKs);
		if(aur != NULL) free(aur);
		if(iur != NULL) free(iur);
		if(aKss != NULL) free(aKss);
		if(iKss != NULL) free(iKss);
		if(C_K2 != NULL) free(C_K2);
		if(C_K1 != NULL) free(C_K1);
		if(O_K != NULL) free(O_K);
		if(I_K != NULL) free(I_K);
	}

	int Solveode::setVariables(int indVariable, double value_new)
	{
		switch(indVariable)
		{
		case 0:		V_ini_ = V_old_= value_new;    break;
		case 1:		Cai_ini_ = Cai_old_= value_new;    break;
		case 2:		Cass_ini_ = Cass_old_= value_new;    break;
		case 3:		CaJSR_ini_ = CaJSR_old_= value_new;    break;
		case 4:		CaNSR_ini_ = CaNSR_old_= value_new;    break;
		case 5:		P_RyR_ini_ = P_RyR_old_= value_new;    break;
		case 6:		LTRPN_Ca_ini_ = LTRPN_Ca_old_= value_new;    break;
		case 7:		HTRPN_Ca_ini_ = HTRPN_Ca_old_= value_new;    break;
		case 8:		P_O1_ini_ = P_O1_old_= value_new;    break;
		case 9:		P_O2_ini_ = P_O2_old_= value_new;    break;
		case 10:		P_C2_ini_ = P_C2_old_= value_new;    break;
		case 11:		O_ini_ = O_old_= value_new;    break;
		case 12:		C2_ini_ = C2_old_= value_new;    break;
		case 13:		C3_ini_ = C3_old_= value_new;    break;
		case 14:		C4_ini_ = C4_old_= value_new;    break;
		case 15:		I1_ini_ = I1_old_= value_new;    break;
		case 16:		I2_ini_ = I2_old_= value_new;    break;
		case 17:		I3_ini_ = I3_old_= value_new;    break;
		case 18:		Nai_ini_ = Nai_old_= value_new;    break;
		case 19:		C_Na2_ini_ = C_Na2_old_= value_new;    break;
		case 20:		C_Na1_ini_ = C_Na1_old_= value_new;    break;
		case 21:		O_Na_ini_ = O_Na_old_= value_new;    break;
		case 22:		IF_Na_ini_ = IF_Na_old_= value_new;    break;
		case 23:		I1_Na_ini_ = I1_Na_old_= value_new;    break;
		case 24:		I2_Na_ini_ = I2_Na_old_= value_new;    break;
		case 25:		IC_Na2_ini_ = IC_Na2_old_= value_new;    break;
		case 26:		IC_Na3_ini_ = IC_Na3_old_= value_new;    break;
		case 27:		Ki_ini_ = Ki_old_= value_new;    break;
		case 28:		ato_f_ini_ = ato_f_old_= value_new;    break;
		case 29:		ito_f_ini_ = ito_f_old_= value_new;    break;
		case 30:		ato_s_ini_ = ato_s_old_= value_new;    break;
		case 31:		ito_s_ini_ = ito_s_old_= value_new;    break;
		case 32:		nKs_ini_ = nKs_old_= value_new;    break;
		case 33:		aur_ini_ = aur_old_= value_new;    break;
		case 34:		iur_ini_ = iur_old_= value_new;    break;
		case 35:		aKss_ini_ = aKss_old_= value_new;    break;
		case 36:		iKss_ini_ = iKss_old_= value_new;    break;
		case 37:		C_K2_ini_ = C_K2_old_= value_new;    break;
		case 38:		C_K1_ini_ = C_K1_old_= value_new;    break;
		case 39:		O_K_ini_ = O_K_old_= value_new;    break;
		case 40:		I_K_ini_ = I_K_old_= value_new;    break;
		default:	return 1;    break;
		}
		return 0;
	}

	int Solveode::setParameters(int indVariable, double value_new)
	{
		switch(indVariable)
		{
		case 0:		time = value_new;   break;
		case 1:		stim_amplitude = value_new;   break;
		case 2:		stim_start = value_new;   break;
		case 3:		stim_end = value_new;   break;
		case 4:		stim_period = value_new;   break;
		case 5:		stim_duration = value_new;   break;
		case 6:		Acap = value_new;   break;
		case 7:		Cm = value_new;   break;
		case 8:		Vmyo = value_new;   break;
		case 9:		F = value_new;   break;
		case 10:		VJSR = value_new;   break;
		case 11:		Vss = value_new;   break;
		case 12:		VNSR = value_new;   break;
		case 13:		CMDN_tot = value_new;   break;
		case 14:		Km_CMDN = value_new;   break;
		case 15:		CSQN_tot = value_new;   break;
		case 16:		Km_CSQN = value_new;   break;
		case 17:		v1 = value_new;   break;
		case 18:		tau_tr = value_new;   break;
		case 19:		tau_xfer = value_new;   break;
		case 20:		v2 = value_new;   break;
		case 21:		v3 = value_new;   break;
		case 22:		Km_up = value_new;   break;
		case 23:		k_plus_htrpn = value_new;   break;
		case 24:		HTRPN_tot = value_new;   break;
		case 25:		k_plus_ltrpn = value_new;   break;
		case 26:		LTRPN_tot = value_new;   break;
		case 27:		k_minus_htrpn = value_new;   break;
		case 28:		k_minus_ltrpn = value_new;   break;
		case 29:		i_CaL_max = value_new;   break;
		case 30:		k_plus_a = value_new;   break;
		case 31:		n = value_new;   break;
		case 32:		k_minus_b = value_new;   break;
		case 33:		k_minus_c = value_new;   break;
		case 34:		k_minus_a = value_new;   break;
		case 35:		k_plus_b = value_new;   break;
		case 36:		m = value_new;   break;
		case 37:		k_plus_c = value_new;   break;
		case 38:		g_CaL = value_new;   break;
		case 39:		E_CaL = value_new;   break;
		case 40:		Kpcb = value_new;   break;
		case 41:		Kpc_max = value_new;   break;
		case 42:		Kpc_half = value_new;   break;
		case 43:		i_pCa_max = value_new;   break;
		case 44:		Km_pCa = value_new;   break;
		case 45:		k_NaCa = value_new;   break;
		case 46:		K_mNa = value_new;   break;
		case 47:		Nao = value_new;   break;
		case 48:		K_mCa = value_new;   break;
		case 49:		Cao = value_new;   break;
		case 50:		k_sat = value_new;   break;
		case 51:		eta = value_new;   break;
		case 52:		R = value_new;   break;
		case 53:		T = value_new;   break;
		case 54:		g_Cab = value_new;   break;
		case 55:		g_Na = value_new;   break;
		case 56:		Ko = value_new;   break;
		case 57:		g_Nab = value_new;   break;
		case 58:		g_Kto_f = value_new;   break;
		case 59:		g_Kto_s = value_new;   break;
		case 60:		g_Ks = value_new;   break;
		case 61:		g_Kur = value_new;   break;
		case 62:		g_Kss = value_new;   break;
		case 63:		g_Kr = value_new;   break;
		case 64:		kf = value_new;   break;
		case 65:		kb = value_new;   break;
		case 66:		i_NaK_max = value_new;   break;
		case 67:		Km_Nai = value_new;   break;
		case 68:		Km_Ko = value_new;   break;
		case 69:		g_ClCa = value_new;   break;
		case 70:		Km_Cl = value_new;   break;
		case 71:		E_Cl = value_new;   break;
		default:	return 1;    break;
		}
		return 0;
	}

	int Solveode::setFreeVariable(double value_new)
	{
		dtime = value_new;
	}

	void Solveode::setMaxStep(double maxSt)
	{
		this->maxStep=maxSt;
	}
	void Solveode::solver(int method, double finalTime, double svRate, char* fileName, int threads)
	{
		this->savingRate = svRate;
		if(finalTime <= 0){
			fprintf(stderr,"ERROR - solveToFile - negative finalTime is not allowed\n");
			exit(1);
		}
		if(savingRate < 0){
			fprintf(stderr,"ERROR - solveToFile - negative saving rate is not allowed\n");
			exit(1);
		}
		FILE *fileptr;
		if(savingRate!=0.0)
			fileptr = fopen(fileName, "w+");
// 		system("ps aux | grep ./main");
		if(threads==1)
		{
			if(method==_EULER_){
				this->euler(finalTime, fileptr);
			}else if(method==_RK2_){
				this->rungeKutta2ndOrder(finalTime, fileptr);
			}else if(method == _ADAP_DT_){
				this->addt(finalTime, fileptr);
			}else if(method == _ADAP_DT2_){
				this->addt2(finalTime, fileptr);
			}else{
				printf("Invalid option!\n");
			}
		}else//OMP
		{
			tree_thread = (int*)malloc(sizeof(double)*forestSize);
			int *_jobs = (int*)malloc(sizeof(double)*forestSize);
			for(int i=0; i< forestSize; i++){
				timeval t1, t2;
				gettimeofday(&t1, NULL);
				for (int j=0; j< 100; j++)
					forest[i].function(this); //invoca cada uma das arvores da floresta
				//stop timer
				gettimeofday(&t2, NULL);
				//compute and print the elapsed time in millisec
				double elapsedTime = (t2.tv_sec - t1.tv_sec);      // sec to ms
				//elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
				_jobs[i] = (int)(t2.tv_usec - t1.tv_usec) + 1;
			}
			greedyPartitionProblem(_jobs, forestSize, threads, tree_thread);
			if(method==_EULER_){
				this->euler_OMP(finalTime, fileptr, threads);
			}else if(method == _ADAP_DT_){
				this->addt_OMP(finalTime, fileptr, threads);
			}else if(method == _ADAP_DT2_){
				this->addt2_OMP(finalTime, fileptr, threads);
			}else{
				printf("Invalid option!\n");
			}
		}
// 		system("ps aux | grep ./main");
		if(savingRate!=0.0)
			fclose(fileptr);
	}
	void Solveode::save_step(FILE *file, int method){
		double time_aux = this->time_new;
		double diff =  _agos_round(this->time_new - timeSaving, 10);
		if(diff==0){
			fprintf(file,  "%.8f %.8f %.8f\n", this->time_new, this->V_new_,previous_dt);
			this->timeSaving += savingRate;
		}else if(diff > 0){
			//salva resultados que passaram do tempo de salvar
			double *edos_new_aux_;
			double *edos_old_aux_;
			edos_old_aux_ = (double*)malloc(sizeof(double)*numEDO);
			edos_new_aux_ = (double*)malloc(sizeof(double)*numEDO);
			for(int i = 0; i < numEDO; i++){
				edos_old_aux_[i] = this->getVariables(i);
			}
			double _righthandside_[numEDO];
			double told = this->time_new - this->previous_dt;
			//encontra-se o dt adequado pra salvar na hora certa
			double dtTemp =   this->timeSaving - (told);
			this->time_new = told;
			//verifica se a iteração acima alcançou o tempo real
			double euler_ladodireito[numEDO], rk2_res[numEDO], old_aux[numEDO];
			while(time_aux >= timeSaving){
				this->time_new += dtTemp;
				//calcula mais uma iteração
				rightHandSideFunction.function(this);
				//euler
				for(int i = 0; i < numEDO; i++){
					old_aux[i] = this->getVariables(i);
					edos_new_aux_[i] = this->getVariables(i) + this->getLadoDireito(i)*dtTemp;
					euler_ladodireito[i] = this->getLadoDireito(i);
					this->setVariables(i, edos_new_aux_[i]);
				}
				//se RK2
				if(method==_RK2_)
				{
					this->time_new += dtTemp;
					rightHandSideFunction.function(this);
					//rk2
					for(int i = 0; i < numEDO; i++){
						rk2_res[i] = old_aux[i] +(this->getLadoDireito(i) + euler_ladodireito[i])*(dtTemp/2);
					}		
					this->time_new -= dtTemp;
				}
				if(method==_RK2_ ){
					for(int i = 0; i < numEDO; i++){
						this->setVariables(i, rk2_res[i] );
					}   
					fprintf(file,"%.8e %.8e %.8e\n", time_new, rk2_res[0],previous_dt);
				}else{
					for(int i = 0; i < numEDO; i++){
						this->setVariables(i, edos_new_aux_[i] );
					}   
					fprintf(file,"%.8e %.8e %.8e\n", time_new, edos_new_aux_[0],previous_dt);
				}
				this->timeSaving += savingRate;
				//coloca o dtime igual ao savingRate para bater exatamente com a proxima iteração de salvar
				dtTemp = savingRate;
			}//fimwhile
			//volta com os valores old originais
			for(int i = 0; i < numEDO; i++){
				this->setVariables(i, edos_old_aux_[i] );
			}
			this->time_new = time_aux;
		}//fim else if
	}
	void Solveode::euler(double finalTime, FILE *fileptr){
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->Cai_new_ = this->Cai_old_ = this->Cai_ini_;
		this->Cass_new_ = this->Cass_old_ = this->Cass_ini_;
		this->CaJSR_new_ = this->CaJSR_old_ = this->CaJSR_ini_;
		this->CaNSR_new_ = this->CaNSR_old_ = this->CaNSR_ini_;
		this->P_RyR_new_ = this->P_RyR_old_ = this->P_RyR_ini_;
		this->LTRPN_Ca_new_ = this->LTRPN_Ca_old_ = this->LTRPN_Ca_ini_;
		this->HTRPN_Ca_new_ = this->HTRPN_Ca_old_ = this->HTRPN_Ca_ini_;
		this->P_O1_new_ = this->P_O1_old_ = this->P_O1_ini_;
		this->P_O2_new_ = this->P_O2_old_ = this->P_O2_ini_;
		this->P_C2_new_ = this->P_C2_old_ = this->P_C2_ini_;
		this->O_new_ = this->O_old_ = this->O_ini_;
		this->C2_new_ = this->C2_old_ = this->C2_ini_;
		this->C3_new_ = this->C3_old_ = this->C3_ini_;
		this->C4_new_ = this->C4_old_ = this->C4_ini_;
		this->I1_new_ = this->I1_old_ = this->I1_ini_;
		this->I2_new_ = this->I2_old_ = this->I2_ini_;
		this->I3_new_ = this->I3_old_ = this->I3_ini_;
		this->Nai_new_ = this->Nai_old_ = this->Nai_ini_;
		this->C_Na2_new_ = this->C_Na2_old_ = this->C_Na2_ini_;
		this->C_Na1_new_ = this->C_Na1_old_ = this->C_Na1_ini_;
		this->O_Na_new_ = this->O_Na_old_ = this->O_Na_ini_;
		this->IF_Na_new_ = this->IF_Na_old_ = this->IF_Na_ini_;
		this->I1_Na_new_ = this->I1_Na_old_ = this->I1_Na_ini_;
		this->I2_Na_new_ = this->I2_Na_old_ = this->I2_Na_ini_;
		this->IC_Na2_new_ = this->IC_Na2_old_ = this->IC_Na2_ini_;
		this->IC_Na3_new_ = this->IC_Na3_old_ = this->IC_Na3_ini_;
		this->Ki_new_ = this->Ki_old_ = this->Ki_ini_;
		this->ato_f_new_ = this->ato_f_old_ = this->ato_f_ini_;
		this->ito_f_new_ = this->ito_f_old_ = this->ito_f_ini_;
		this->ato_s_new_ = this->ato_s_old_ = this->ato_s_ini_;
		this->ito_s_new_ = this->ito_s_old_ = this->ito_s_ini_;
		this->nKs_new_ = this->nKs_old_ = this->nKs_ini_;
		this->aur_new_ = this->aur_old_ = this->aur_ini_;
		this->iur_new_ = this->iur_old_ = this->iur_ini_;
		this->aKss_new_ = this->aKss_old_ = this->aKss_ini_;
		this->iKss_new_ = this->iKss_old_ = this->iKss_ini_;
		this->C_K2_new_ = this->C_K2_old_ = this->C_K2_ini_;
		this->C_K1_new_ = this->C_K1_old_ = this->C_K1_ini_;
		this->O_K_new_ = this->O_K_old_ = this->O_K_ini_;
		this->I_K_new_ = this->I_K_old_ = this->I_K_ini_;
		this->time_new = this->time;
		this->timeSaving = this->time;
// 		if(savingRate!=0.0)
// 			this->save_step(fileptr, _EULER_);//save the initial conditions
// 		int contador=0;
		while(this->time_new<=finalTime){
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			this->V_new_ = this->dtime*(this->V_lado_direito_) + this->V_old_;
			this->Cai_new_ = this->dtime*(this->Cai_lado_direito_) + this->Cai_old_;
			this->Cass_new_ = this->dtime*(this->Cass_lado_direito_) + this->Cass_old_;
			this->CaJSR_new_ = this->dtime*(this->CaJSR_lado_direito_) + this->CaJSR_old_;
			this->CaNSR_new_ = this->dtime*(this->CaNSR_lado_direito_) + this->CaNSR_old_;
			this->P_RyR_new_ = this->dtime*(this->P_RyR_lado_direito_) + this->P_RyR_old_;
			this->LTRPN_Ca_new_ = this->dtime*(this->LTRPN_Ca_lado_direito_) + this->LTRPN_Ca_old_;
			this->HTRPN_Ca_new_ = this->dtime*(this->HTRPN_Ca_lado_direito_) + this->HTRPN_Ca_old_;
			this->P_O1_new_ = this->dtime*(this->P_O1_lado_direito_) + this->P_O1_old_;
			this->P_O2_new_ = this->dtime*(this->P_O2_lado_direito_) + this->P_O2_old_;
			this->P_C2_new_ = this->dtime*(this->P_C2_lado_direito_) + this->P_C2_old_;
			this->O_new_ = this->dtime*(this->O_lado_direito_) + this->O_old_;
			this->C2_new_ = this->dtime*(this->C2_lado_direito_) + this->C2_old_;
			this->C3_new_ = this->dtime*(this->C3_lado_direito_) + this->C3_old_;
			this->C4_new_ = this->dtime*(this->C4_lado_direito_) + this->C4_old_;
			this->I1_new_ = this->dtime*(this->I1_lado_direito_) + this->I1_old_;
			this->I2_new_ = this->dtime*(this->I2_lado_direito_) + this->I2_old_;
			this->I3_new_ = this->dtime*(this->I3_lado_direito_) + this->I3_old_;
			this->Nai_new_ = this->dtime*(this->Nai_lado_direito_) + this->Nai_old_;
			this->C_Na2_new_ = this->dtime*(this->C_Na2_lado_direito_) + this->C_Na2_old_;
			this->C_Na1_new_ = this->dtime*(this->C_Na1_lado_direito_) + this->C_Na1_old_;
			this->O_Na_new_ = this->dtime*(this->O_Na_lado_direito_) + this->O_Na_old_;
			this->IF_Na_new_ = this->dtime*(this->IF_Na_lado_direito_) + this->IF_Na_old_;
			this->I1_Na_new_ = this->dtime*(this->I1_Na_lado_direito_) + this->I1_Na_old_;
			this->I2_Na_new_ = this->dtime*(this->I2_Na_lado_direito_) + this->I2_Na_old_;
			this->IC_Na2_new_ = this->dtime*(this->IC_Na2_lado_direito_) + this->IC_Na2_old_;
			this->IC_Na3_new_ = this->dtime*(this->IC_Na3_lado_direito_) + this->IC_Na3_old_;
			this->Ki_new_ = this->dtime*(this->Ki_lado_direito_) + this->Ki_old_;
			this->ato_f_new_ = this->dtime*(this->ato_f_lado_direito_) + this->ato_f_old_;
			this->ito_f_new_ = this->dtime*(this->ito_f_lado_direito_) + this->ito_f_old_;
			this->ato_s_new_ = this->dtime*(this->ato_s_lado_direito_) + this->ato_s_old_;
			this->ito_s_new_ = this->dtime*(this->ito_s_lado_direito_) + this->ito_s_old_;
			this->nKs_new_ = this->dtime*(this->nKs_lado_direito_) + this->nKs_old_;
			this->aur_new_ = this->dtime*(this->aur_lado_direito_) + this->aur_old_;
			this->iur_new_ = this->dtime*(this->iur_lado_direito_) + this->iur_old_;
			this->aKss_new_ = this->dtime*(this->aKss_lado_direito_) + this->aKss_old_;
			this->iKss_new_ = this->dtime*(this->iKss_lado_direito_) + this->iKss_old_;
			this->C_K2_new_ = this->dtime*(this->C_K2_lado_direito_) + this->C_K2_old_;
			this->C_K1_new_ = this->dtime*(this->C_K1_lado_direito_) + this->C_K1_old_;
			this->O_K_new_ = this->dtime*(this->O_K_lado_direito_) + this->O_K_old_;
			this->I_K_new_ = this->dtime*(this->I_K_lado_direito_) + this->I_K_old_;
			//save results on a file
// 			contador++;
// 			if(savingRate!=0.0)
// 			{
// 				this->save_step(fileptr, _EULER_);
// 			}
		    this->V_old_ = this->V_new_;
		    this->Cai_old_ = this->Cai_new_;
		    this->Cass_old_ = this->Cass_new_;
		    this->CaJSR_old_ = this->CaJSR_new_;
		    this->CaNSR_old_ = this->CaNSR_new_;
		    this->P_RyR_old_ = this->P_RyR_new_;
		    this->LTRPN_Ca_old_ = this->LTRPN_Ca_new_;
		    this->HTRPN_Ca_old_ = this->HTRPN_Ca_new_;
		    this->P_O1_old_ = this->P_O1_new_;
		    this->P_O2_old_ = this->P_O2_new_;
		    this->P_C2_old_ = this->P_C2_new_;
		    this->O_old_ = this->O_new_;
		    this->C2_old_ = this->C2_new_;
		    this->C3_old_ = this->C3_new_;
		    this->C4_old_ = this->C4_new_;
		    this->I1_old_ = this->I1_new_;
		    this->I2_old_ = this->I2_new_;
		    this->I3_old_ = this->I3_new_;
		    this->Nai_old_ = this->Nai_new_;
		    this->C_Na2_old_ = this->C_Na2_new_;
		    this->C_Na1_old_ = this->C_Na1_new_;
		    this->O_Na_old_ = this->O_Na_new_;
		    this->IF_Na_old_ = this->IF_Na_new_;
		    this->I1_Na_old_ = this->I1_Na_new_;
		    this->I2_Na_old_ = this->I2_Na_new_;
		    this->IC_Na2_old_ = this->IC_Na2_new_;
		    this->IC_Na3_old_ = this->IC_Na3_new_;
		    this->Ki_old_ = this->Ki_new_;
		    this->ato_f_old_ = this->ato_f_new_;
		    this->ito_f_old_ = this->ito_f_new_;
		    this->ato_s_old_ = this->ato_s_new_;
		    this->ito_s_old_ = this->ito_s_new_;
		    this->nKs_old_ = this->nKs_new_;
		    this->aur_old_ = this->aur_new_;
		    this->iur_old_ = this->iur_new_;
		    this->aKss_old_ = this->aKss_new_;
		    this->iKss_old_ = this->iKss_new_;
		    this->C_K2_old_ = this->C_K2_new_;
		    this->C_K1_old_ = this->C_K1_new_;
		    this->O_K_old_ = this->O_K_new_;
		    this->I_K_old_ = this->I_K_new_;
		}
// 		printf("%d iterações\n", contador);
	}
	void Solveode::rungeKutta2ndOrder(double finalTime, FILE *fileptr){
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->Cai_new_ = this->Cai_old_ = this->Cai_ini_;
		this->Cass_new_ = this->Cass_old_ = this->Cass_ini_;
		this->CaJSR_new_ = this->CaJSR_old_ = this->CaJSR_ini_;
		this->CaNSR_new_ = this->CaNSR_old_ = this->CaNSR_ini_;
		this->P_RyR_new_ = this->P_RyR_old_ = this->P_RyR_ini_;
		this->LTRPN_Ca_new_ = this->LTRPN_Ca_old_ = this->LTRPN_Ca_ini_;
		this->HTRPN_Ca_new_ = this->HTRPN_Ca_old_ = this->HTRPN_Ca_ini_;
		this->P_O1_new_ = this->P_O1_old_ = this->P_O1_ini_;
		this->P_O2_new_ = this->P_O2_old_ = this->P_O2_ini_;
		this->P_C2_new_ = this->P_C2_old_ = this->P_C2_ini_;
		this->O_new_ = this->O_old_ = this->O_ini_;
		this->C2_new_ = this->C2_old_ = this->C2_ini_;
		this->C3_new_ = this->C3_old_ = this->C3_ini_;
		this->C4_new_ = this->C4_old_ = this->C4_ini_;
		this->I1_new_ = this->I1_old_ = this->I1_ini_;
		this->I2_new_ = this->I2_old_ = this->I2_ini_;
		this->I3_new_ = this->I3_old_ = this->I3_ini_;
		this->Nai_new_ = this->Nai_old_ = this->Nai_ini_;
		this->C_Na2_new_ = this->C_Na2_old_ = this->C_Na2_ini_;
		this->C_Na1_new_ = this->C_Na1_old_ = this->C_Na1_ini_;
		this->O_Na_new_ = this->O_Na_old_ = this->O_Na_ini_;
		this->IF_Na_new_ = this->IF_Na_old_ = this->IF_Na_ini_;
		this->I1_Na_new_ = this->I1_Na_old_ = this->I1_Na_ini_;
		this->I2_Na_new_ = this->I2_Na_old_ = this->I2_Na_ini_;
		this->IC_Na2_new_ = this->IC_Na2_old_ = this->IC_Na2_ini_;
		this->IC_Na3_new_ = this->IC_Na3_old_ = this->IC_Na3_ini_;
		this->Ki_new_ = this->Ki_old_ = this->Ki_ini_;
		this->ato_f_new_ = this->ato_f_old_ = this->ato_f_ini_;
		this->ito_f_new_ = this->ito_f_old_ = this->ito_f_ini_;
		this->ato_s_new_ = this->ato_s_old_ = this->ato_s_ini_;
		this->ito_s_new_ = this->ito_s_old_ = this->ito_s_ini_;
		this->nKs_new_ = this->nKs_old_ = this->nKs_ini_;
		this->aur_new_ = this->aur_old_ = this->aur_ini_;
		this->iur_new_ = this->iur_old_ = this->iur_ini_;
		this->aKss_new_ = this->aKss_old_ = this->aKss_ini_;
		this->iKss_new_ = this->iKss_old_ = this->iKss_ini_;
		this->C_K2_new_ = this->C_K2_old_ = this->C_K2_ini_;
		this->C_K1_new_ = this->C_K1_old_ = this->C_K1_ini_;
		this->O_K_new_ = this->O_K_old_ = this->O_K_ini_;
		this->I_K_new_ = this->I_K_old_ = this->I_K_ini_;
		this->time_new = this->time;
		this->timeSaving = this->time;
		if(savingRate!=0.0)
		this->save_step(fileptr, _RK2_);//save the initial conditions
		double edos_old_aux_[numEDO];
		double edos_rightside_aux_[numEDO];
		while(this->time_new<=finalTime){
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			this->V_new_ = this->dtime*(this->V_lado_direito_) + this->V_old_;
			this->Cai_new_ = this->dtime*(this->Cai_lado_direito_) + this->Cai_old_;
			this->Cass_new_ = this->dtime*(this->Cass_lado_direito_) + this->Cass_old_;
			this->CaJSR_new_ = this->dtime*(this->CaJSR_lado_direito_) + this->CaJSR_old_;
			this->CaNSR_new_ = this->dtime*(this->CaNSR_lado_direito_) + this->CaNSR_old_;
			this->P_RyR_new_ = this->dtime*(this->P_RyR_lado_direito_) + this->P_RyR_old_;
			this->LTRPN_Ca_new_ = this->dtime*(this->LTRPN_Ca_lado_direito_) + this->LTRPN_Ca_old_;
			this->HTRPN_Ca_new_ = this->dtime*(this->HTRPN_Ca_lado_direito_) + this->HTRPN_Ca_old_;
			this->P_O1_new_ = this->dtime*(this->P_O1_lado_direito_) + this->P_O1_old_;
			this->P_O2_new_ = this->dtime*(this->P_O2_lado_direito_) + this->P_O2_old_;
			this->P_C2_new_ = this->dtime*(this->P_C2_lado_direito_) + this->P_C2_old_;
			this->O_new_ = this->dtime*(this->O_lado_direito_) + this->O_old_;
			this->C2_new_ = this->dtime*(this->C2_lado_direito_) + this->C2_old_;
			this->C3_new_ = this->dtime*(this->C3_lado_direito_) + this->C3_old_;
			this->C4_new_ = this->dtime*(this->C4_lado_direito_) + this->C4_old_;
			this->I1_new_ = this->dtime*(this->I1_lado_direito_) + this->I1_old_;
			this->I2_new_ = this->dtime*(this->I2_lado_direito_) + this->I2_old_;
			this->I3_new_ = this->dtime*(this->I3_lado_direito_) + this->I3_old_;
			this->Nai_new_ = this->dtime*(this->Nai_lado_direito_) + this->Nai_old_;
			this->C_Na2_new_ = this->dtime*(this->C_Na2_lado_direito_) + this->C_Na2_old_;
			this->C_Na1_new_ = this->dtime*(this->C_Na1_lado_direito_) + this->C_Na1_old_;
			this->O_Na_new_ = this->dtime*(this->O_Na_lado_direito_) + this->O_Na_old_;
			this->IF_Na_new_ = this->dtime*(this->IF_Na_lado_direito_) + this->IF_Na_old_;
			this->I1_Na_new_ = this->dtime*(this->I1_Na_lado_direito_) + this->I1_Na_old_;
			this->I2_Na_new_ = this->dtime*(this->I2_Na_lado_direito_) + this->I2_Na_old_;
			this->IC_Na2_new_ = this->dtime*(this->IC_Na2_lado_direito_) + this->IC_Na2_old_;
			this->IC_Na3_new_ = this->dtime*(this->IC_Na3_lado_direito_) + this->IC_Na3_old_;
			this->Ki_new_ = this->dtime*(this->Ki_lado_direito_) + this->Ki_old_;
			this->ato_f_new_ = this->dtime*(this->ato_f_lado_direito_) + this->ato_f_old_;
			this->ito_f_new_ = this->dtime*(this->ito_f_lado_direito_) + this->ito_f_old_;
			this->ato_s_new_ = this->dtime*(this->ato_s_lado_direito_) + this->ato_s_old_;
			this->ito_s_new_ = this->dtime*(this->ito_s_lado_direito_) + this->ito_s_old_;
			this->nKs_new_ = this->dtime*(this->nKs_lado_direito_) + this->nKs_old_;
			this->aur_new_ = this->dtime*(this->aur_lado_direito_) + this->aur_old_;
			this->iur_new_ = this->dtime*(this->iur_lado_direito_) + this->iur_old_;
			this->aKss_new_ = this->dtime*(this->aKss_lado_direito_) + this->aKss_old_;
			this->iKss_new_ = this->dtime*(this->iKss_lado_direito_) + this->iKss_old_;
			this->C_K2_new_ = this->dtime*(this->C_K2_lado_direito_) + this->C_K2_old_;
			this->C_K1_new_ = this->dtime*(this->C_K1_lado_direito_) + this->C_K1_old_;
			this->O_K_new_ = this->dtime*(this->O_K_lado_direito_) + this->O_K_old_;
			this->I_K_new_ = this->dtime*(this->I_K_lado_direito_) + this->I_K_old_;
			//stores the old variables in a vector
			for(int i=0;i<numEDO;i++){
				edos_old_aux_[i] = this->getVariables(i);
				edos_rightside_aux_[i] = this->getLadoDireito(i);
			}
			//steps one iteration ahead;
			this->V_old_ = this->V_new_;
			this->Cai_old_ = this->Cai_new_;
			this->Cass_old_ = this->Cass_new_;
			this->CaJSR_old_ = this->CaJSR_new_;
			this->CaNSR_old_ = this->CaNSR_new_;
			this->P_RyR_old_ = this->P_RyR_new_;
			this->LTRPN_Ca_old_ = this->LTRPN_Ca_new_;
			this->HTRPN_Ca_old_ = this->HTRPN_Ca_new_;
			this->P_O1_old_ = this->P_O1_new_;
			this->P_O2_old_ = this->P_O2_new_;
			this->P_C2_old_ = this->P_C2_new_;
			this->O_old_ = this->O_new_;
			this->C2_old_ = this->C2_new_;
			this->C3_old_ = this->C3_new_;
			this->C4_old_ = this->C4_new_;
			this->I1_old_ = this->I1_new_;
			this->I2_old_ = this->I2_new_;
			this->I3_old_ = this->I3_new_;
			this->Nai_old_ = this->Nai_new_;
			this->C_Na2_old_ = this->C_Na2_new_;
			this->C_Na1_old_ = this->C_Na1_new_;
			this->O_Na_old_ = this->O_Na_new_;
			this->IF_Na_old_ = this->IF_Na_new_;
			this->I1_Na_old_ = this->I1_Na_new_;
			this->I2_Na_old_ = this->I2_Na_new_;
			this->IC_Na2_old_ = this->IC_Na2_new_;
			this->IC_Na3_old_ = this->IC_Na3_new_;
			this->Ki_old_ = this->Ki_new_;
			this->ato_f_old_ = this->ato_f_new_;
			this->ito_f_old_ = this->ito_f_new_;
			this->ato_s_old_ = this->ato_s_new_;
			this->ito_s_old_ = this->ito_s_new_;
			this->nKs_old_ = this->nKs_new_;
			this->aur_old_ = this->aur_new_;
			this->iur_old_ = this->iur_new_;
			this->aKss_old_ = this->aKss_new_;
			this->iKss_old_ = this->iKss_new_;
			this->C_K2_old_ = this->C_K2_new_;
			this->C_K1_old_ = this->C_K1_new_;
			this->O_K_old_ = this->O_K_new_;
			this->I_K_old_ = this->I_K_new_;
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			//computes the runge kutta second order method
			this->V_new_ = ( this->V_lado_direito_ + edos_rightside_aux_[0] ) * this->dtime/2 + edos_old_aux_[0];
			this->Cai_new_ = ( this->Cai_lado_direito_ + edos_rightside_aux_[1] ) * this->dtime/2 + edos_old_aux_[1];
			this->Cass_new_ = ( this->Cass_lado_direito_ + edos_rightside_aux_[2] ) * this->dtime/2 + edos_old_aux_[2];
			this->CaJSR_new_ = ( this->CaJSR_lado_direito_ + edos_rightside_aux_[3] ) * this->dtime/2 + edos_old_aux_[3];
			this->CaNSR_new_ = ( this->CaNSR_lado_direito_ + edos_rightside_aux_[4] ) * this->dtime/2 + edos_old_aux_[4];
			this->P_RyR_new_ = ( this->P_RyR_lado_direito_ + edos_rightside_aux_[5] ) * this->dtime/2 + edos_old_aux_[5];
			this->LTRPN_Ca_new_ = ( this->LTRPN_Ca_lado_direito_ + edos_rightside_aux_[6] ) * this->dtime/2 + edos_old_aux_[6];
			this->HTRPN_Ca_new_ = ( this->HTRPN_Ca_lado_direito_ + edos_rightside_aux_[7] ) * this->dtime/2 + edos_old_aux_[7];
			this->P_O1_new_ = ( this->P_O1_lado_direito_ + edos_rightside_aux_[8] ) * this->dtime/2 + edos_old_aux_[8];
			this->P_O2_new_ = ( this->P_O2_lado_direito_ + edos_rightside_aux_[9] ) * this->dtime/2 + edos_old_aux_[9];
			this->P_C2_new_ = ( this->P_C2_lado_direito_ + edos_rightside_aux_[10] ) * this->dtime/2 + edos_old_aux_[10];
			this->O_new_ = ( this->O_lado_direito_ + edos_rightside_aux_[11] ) * this->dtime/2 + edos_old_aux_[11];
			this->C2_new_ = ( this->C2_lado_direito_ + edos_rightside_aux_[12] ) * this->dtime/2 + edos_old_aux_[12];
			this->C3_new_ = ( this->C3_lado_direito_ + edos_rightside_aux_[13] ) * this->dtime/2 + edos_old_aux_[13];
			this->C4_new_ = ( this->C4_lado_direito_ + edos_rightside_aux_[14] ) * this->dtime/2 + edos_old_aux_[14];
			this->I1_new_ = ( this->I1_lado_direito_ + edos_rightside_aux_[15] ) * this->dtime/2 + edos_old_aux_[15];
			this->I2_new_ = ( this->I2_lado_direito_ + edos_rightside_aux_[16] ) * this->dtime/2 + edos_old_aux_[16];
			this->I3_new_ = ( this->I3_lado_direito_ + edos_rightside_aux_[17] ) * this->dtime/2 + edos_old_aux_[17];
			this->Nai_new_ = ( this->Nai_lado_direito_ + edos_rightside_aux_[18] ) * this->dtime/2 + edos_old_aux_[18];
			this->C_Na2_new_ = ( this->C_Na2_lado_direito_ + edos_rightside_aux_[19] ) * this->dtime/2 + edos_old_aux_[19];
			this->C_Na1_new_ = ( this->C_Na1_lado_direito_ + edos_rightside_aux_[20] ) * this->dtime/2 + edos_old_aux_[20];
			this->O_Na_new_ = ( this->O_Na_lado_direito_ + edos_rightside_aux_[21] ) * this->dtime/2 + edos_old_aux_[21];
			this->IF_Na_new_ = ( this->IF_Na_lado_direito_ + edos_rightside_aux_[22] ) * this->dtime/2 + edos_old_aux_[22];
			this->I1_Na_new_ = ( this->I1_Na_lado_direito_ + edos_rightside_aux_[23] ) * this->dtime/2 + edos_old_aux_[23];
			this->I2_Na_new_ = ( this->I2_Na_lado_direito_ + edos_rightside_aux_[24] ) * this->dtime/2 + edos_old_aux_[24];
			this->IC_Na2_new_ = ( this->IC_Na2_lado_direito_ + edos_rightside_aux_[25] ) * this->dtime/2 + edos_old_aux_[25];
			this->IC_Na3_new_ = ( this->IC_Na3_lado_direito_ + edos_rightside_aux_[26] ) * this->dtime/2 + edos_old_aux_[26];
			this->Ki_new_ = ( this->Ki_lado_direito_ + edos_rightside_aux_[27] ) * this->dtime/2 + edos_old_aux_[27];
			this->ato_f_new_ = ( this->ato_f_lado_direito_ + edos_rightside_aux_[28] ) * this->dtime/2 + edos_old_aux_[28];
			this->ito_f_new_ = ( this->ito_f_lado_direito_ + edos_rightside_aux_[29] ) * this->dtime/2 + edos_old_aux_[29];
			this->ato_s_new_ = ( this->ato_s_lado_direito_ + edos_rightside_aux_[30] ) * this->dtime/2 + edos_old_aux_[30];
			this->ito_s_new_ = ( this->ito_s_lado_direito_ + edos_rightside_aux_[31] ) * this->dtime/2 + edos_old_aux_[31];
			this->nKs_new_ = ( this->nKs_lado_direito_ + edos_rightside_aux_[32] ) * this->dtime/2 + edos_old_aux_[32];
			this->aur_new_ = ( this->aur_lado_direito_ + edos_rightside_aux_[33] ) * this->dtime/2 + edos_old_aux_[33];
			this->iur_new_ = ( this->iur_lado_direito_ + edos_rightside_aux_[34] ) * this->dtime/2 + edos_old_aux_[34];
			this->aKss_new_ = ( this->aKss_lado_direito_ + edos_rightside_aux_[35] ) * this->dtime/2 + edos_old_aux_[35];
			this->iKss_new_ = ( this->iKss_lado_direito_ + edos_rightside_aux_[36] ) * this->dtime/2 + edos_old_aux_[36];
			this->C_K2_new_ = ( this->C_K2_lado_direito_ + edos_rightside_aux_[37] ) * this->dtime/2 + edos_old_aux_[37];
			this->C_K1_new_ = ( this->C_K1_lado_direito_ + edos_rightside_aux_[38] ) * this->dtime/2 + edos_old_aux_[38];
			this->O_K_new_ = ( this->O_K_lado_direito_ + edos_rightside_aux_[39] ) * this->dtime/2 + edos_old_aux_[39];
			this->I_K_new_ = ( this->I_K_lado_direito_ + edos_rightside_aux_[40] ) * this->dtime/2 + edos_old_aux_[40];
			this->time_new	-= this->dtime;//step back, to save in the right time
			for(int i=0;i<numEDO;i++){
				this->setVariables(i, edos_old_aux_[i]);
			}
			//save results on a file
			if(savingRate!=0.0)
			{
				this->save_step(fileptr,_RK2_);
			}
			this->V_old_ = this->V_new_;
			this->Cai_old_ = this->Cai_new_;
			this->Cass_old_ = this->Cass_new_;
			this->CaJSR_old_ = this->CaJSR_new_;
			this->CaNSR_old_ = this->CaNSR_new_;
			this->P_RyR_old_ = this->P_RyR_new_;
			this->LTRPN_Ca_old_ = this->LTRPN_Ca_new_;
			this->HTRPN_Ca_old_ = this->HTRPN_Ca_new_;
			this->P_O1_old_ = this->P_O1_new_;
			this->P_O2_old_ = this->P_O2_new_;
			this->P_C2_old_ = this->P_C2_new_;
			this->O_old_ = this->O_new_;
			this->C2_old_ = this->C2_new_;
			this->C3_old_ = this->C3_new_;
			this->C4_old_ = this->C4_new_;
			this->I1_old_ = this->I1_new_;
			this->I2_old_ = this->I2_new_;
			this->I3_old_ = this->I3_new_;
			this->Nai_old_ = this->Nai_new_;
			this->C_Na2_old_ = this->C_Na2_new_;
			this->C_Na1_old_ = this->C_Na1_new_;
			this->O_Na_old_ = this->O_Na_new_;
			this->IF_Na_old_ = this->IF_Na_new_;
			this->I1_Na_old_ = this->I1_Na_new_;
			this->I2_Na_old_ = this->I2_Na_new_;
			this->IC_Na2_old_ = this->IC_Na2_new_;
			this->IC_Na3_old_ = this->IC_Na3_new_;
			this->Ki_old_ = this->Ki_new_;
			this->ato_f_old_ = this->ato_f_new_;
			this->ito_f_old_ = this->ito_f_new_;
			this->ato_s_old_ = this->ato_s_new_;
			this->ito_s_old_ = this->ito_s_new_;
			this->nKs_old_ = this->nKs_new_;
			this->aur_old_ = this->aur_new_;
			this->iur_old_ = this->iur_new_;
			this->aKss_old_ = this->aKss_new_;
			this->iKss_old_ = this->iKss_new_;
			this->C_K2_old_ = this->C_K2_new_;
			this->C_K1_old_ = this->C_K1_new_;
			this->O_K_old_ = this->O_K_new_;
			this->I_K_old_ = this->I_K_new_;
		}
	}
	void Solveode::addt(double finalTime, FILE *fileptr){
		double _tolerances_[numEDO];
		double _aux_tol=0.0;
		double maxDt = this->dtime, minDt = this->dtime;
		int desc =0;
		//initializes the variables
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->Cai_new_ = this->Cai_old_ = this->Cai_ini_;
		this->Cass_new_ = this->Cass_old_ = this->Cass_ini_;
		this->CaJSR_new_ = this->CaJSR_old_ = this->CaJSR_ini_;
		this->CaNSR_new_ = this->CaNSR_old_ = this->CaNSR_ini_;
		this->P_RyR_new_ = this->P_RyR_old_ = this->P_RyR_ini_;
		this->LTRPN_Ca_new_ = this->LTRPN_Ca_old_ = this->LTRPN_Ca_ini_;
		this->HTRPN_Ca_new_ = this->HTRPN_Ca_old_ = this->HTRPN_Ca_ini_;
		this->P_O1_new_ = this->P_O1_old_ = this->P_O1_ini_;
		this->P_O2_new_ = this->P_O2_old_ = this->P_O2_ini_;
		this->P_C2_new_ = this->P_C2_old_ = this->P_C2_ini_;
		this->O_new_ = this->O_old_ = this->O_ini_;
		this->C2_new_ = this->C2_old_ = this->C2_ini_;
		this->C3_new_ = this->C3_old_ = this->C3_ini_;
		this->C4_new_ = this->C4_old_ = this->C4_ini_;
		this->I1_new_ = this->I1_old_ = this->I1_ini_;
		this->I2_new_ = this->I2_old_ = this->I2_ini_;
		this->I3_new_ = this->I3_old_ = this->I3_ini_;
		this->Nai_new_ = this->Nai_old_ = this->Nai_ini_;
		this->C_Na2_new_ = this->C_Na2_old_ = this->C_Na2_ini_;
		this->C_Na1_new_ = this->C_Na1_old_ = this->C_Na1_ini_;
		this->O_Na_new_ = this->O_Na_old_ = this->O_Na_ini_;
		this->IF_Na_new_ = this->IF_Na_old_ = this->IF_Na_ini_;
		this->I1_Na_new_ = this->I1_Na_old_ = this->I1_Na_ini_;
		this->I2_Na_new_ = this->I2_Na_old_ = this->I2_Na_ini_;
		this->IC_Na2_new_ = this->IC_Na2_old_ = this->IC_Na2_ini_;
		this->IC_Na3_new_ = this->IC_Na3_old_ = this->IC_Na3_ini_;
		this->Ki_new_ = this->Ki_old_ = this->Ki_ini_;
		this->ato_f_new_ = this->ato_f_old_ = this->ato_f_ini_;
		this->ito_f_new_ = this->ito_f_old_ = this->ito_f_ini_;
		this->ato_s_new_ = this->ato_s_old_ = this->ato_s_ini_;
		this->ito_s_new_ = this->ito_s_old_ = this->ito_s_ini_;
		this->nKs_new_ = this->nKs_old_ = this->nKs_ini_;
		this->aur_new_ = this->aur_old_ = this->aur_ini_;
		this->iur_new_ = this->iur_old_ = this->iur_ini_;
		this->aKss_new_ = this->aKss_old_ = this->aKss_ini_;
		this->iKss_new_ = this->iKss_old_ = this->iKss_ini_;
		this->C_K2_new_ = this->C_K2_old_ = this->C_K2_ini_;
		this->C_K1_new_ = this->C_K1_old_ = this->C_K1_ini_;
		this->O_K_new_ = this->O_K_old_ = this->O_K_ini_;
		this->I_K_new_ = this->I_K_old_ = this->I_K_ini_;
		this->time_new = this->time;
		this->timeSaving = this->time;
		if(savingRate!=0)
			save_step(fileptr, _ADAP_DT_);//save the initial conditions
		double edos_old_aux_[numEDO];
		double edos_new_rk2_[numEDO],edos_new_euler_[numEDO];
		double *_k1__ = (double*)malloc(sizeof(double)*numEDO);
		double *_k2__ = (double*)malloc(sizeof(double)*numEDO);
		double *_k_aux__;
		int iMaiorErro=0;
		int _cont_=0, cont_inv=0;
		double _soma_=0.0;
		this->time_new += this->dtime;
		rightHandSideFunction.function(this);
		for(int i=0;i<numEDO;i++){
			_k1__[i] = this->getLadoDireito(i);
		}
		const double __tiny_ = pow(this->abstol__, 2.0);
		while(this->time_new<=finalTime){
			for(int i=0; i<numEDO; i++){
				//stores the old variables in a vector
				edos_old_aux_[i] = this->getVariables(i);
				//computes euler method
				edos_new_euler_[i] = _k1__[i] * this->dtime + edos_old_aux_[i];
				//steps ahead to compute the rk2 method
				this->setVariables(i, edos_new_euler_[i]);
			}
			//steps time ahead
			this->time_new += this->dtime;
			//computes the right-hand side one step ahead
			rightHandSideFunction.function(this);
			//restore the original old value
			this->time_new -= this->dtime;//step back
			double greatestError=0.0, auxError=0.0;
			for(int i=0;i<numEDO;i++){
				//stores the new evaluation 
				_k2__[i] = this->getLadoDireito(i);
				_aux_tol = fabs(edos_new_euler_[i])*reltol__;
				_tolerances_[i] = (abstol__ > _aux_tol )?abstol__:_aux_tol;
				//finds the greatest error between  the steps
				auxError = fabs(( (this->dtime/2)*(_k1__[i] - _k2__[i])) / _tolerances_[i]);
				if(auxError > greatestError){
					greatestError = auxError;
					iMaiorErro = i;
				}
			}
			///adapt the time step
			greatestError+=__tiny_;
			int flag = this->getErrorCode(greatestError, 1.0);
			this->previous_dt = this->dtime;
			//it doesn't accept the solution and cut h in a half
			if(flag==-1){
			    cont_inv++;
				//restore the old values to do it again
				for(int i=0; i<numEDO; i++)
					this->setVariables(i, edos_old_aux_[i]);
				//throw the results away and compute again
				this->dtime = this->dtime/_DECREASE_DT_2_;//cut time step in a half
				desc++;
			}else{//it accepts the solutions
				if(savingRate!=0){
					//restore the previous value of old variables
					for(int i=0; i<numEDO; i++)
						this->setVariables(i, edos_old_aux_[i]);
						this->V_new_ = edos_new_euler_[0];
						this->Cai_new_ = edos_new_euler_[1];
						this->Cass_new_ = edos_new_euler_[2];
						this->CaJSR_new_ = edos_new_euler_[3];
						this->CaNSR_new_ = edos_new_euler_[4];
						this->P_RyR_new_ = edos_new_euler_[5];
						this->LTRPN_Ca_new_ = edos_new_euler_[6];
						this->HTRPN_Ca_new_ = edos_new_euler_[7];
						this->P_O1_new_ = edos_new_euler_[8];
						this->P_O2_new_ = edos_new_euler_[9];
						this->P_C2_new_ = edos_new_euler_[10];
						this->O_new_ = edos_new_euler_[11];
						this->C2_new_ = edos_new_euler_[12];
						this->C3_new_ = edos_new_euler_[13];
						this->C4_new_ = edos_new_euler_[14];
						this->I1_new_ = edos_new_euler_[15];
						this->I2_new_ = edos_new_euler_[16];
						this->I3_new_ = edos_new_euler_[17];
						this->Nai_new_ = edos_new_euler_[18];
						this->C_Na2_new_ = edos_new_euler_[19];
						this->C_Na1_new_ = edos_new_euler_[20];
						this->O_Na_new_ = edos_new_euler_[21];
						this->IF_Na_new_ = edos_new_euler_[22];
						this->I1_Na_new_ = edos_new_euler_[23];
						this->I2_Na_new_ = edos_new_euler_[24];
						this->IC_Na2_new_ = edos_new_euler_[25];
						this->IC_Na3_new_ = edos_new_euler_[26];
						this->Ki_new_ = edos_new_euler_[27];
						this->ato_f_new_ = edos_new_euler_[28];
						this->ito_f_new_ = edos_new_euler_[29];
						this->ato_s_new_ = edos_new_euler_[30];
						this->ito_s_new_ = edos_new_euler_[31];
						this->nKs_new_ = edos_new_euler_[32];
						this->aur_new_ = edos_new_euler_[33];
						this->iur_new_ = edos_new_euler_[34];
						this->aKss_new_ = edos_new_euler_[35];
						this->iKss_new_ = edos_new_euler_[36];
						this->C_K2_new_ = edos_new_euler_[37];
						this->C_K1_new_ = edos_new_euler_[38];
						this->O_K_new_ = edos_new_euler_[39];
						this->I_K_new_ = edos_new_euler_[40];
					save_step(fileptr, _ADAP_DT_);
					if(this->dtime>maxDt) maxDt = this->dtime;
					if(this->dtime<minDt) minDt = this->dtime;
					_soma_+=this->dtime;
					_cont_++;
				}
				if(flag==3){
					this->dtime = this->dtime*_DECREASE_DT_;
				}else if(flag==4){
					this->dtime = this->dtime*_INCREASE_DT_;
				}else if(flag==0){
					//it just doesnt do anything
				}else{
					printf("flag: %d\n", flag);
				}
				if(this->dtime > maxStep && maxStep!=0){
					this->dtime = maxStep;
				}else if(this->dtime==0){
					printf("Error: Time step is zero.\n");
					return;
				}
				//change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração
				_k_aux__	= _k2__;
				_k2__		= _k1__;
				_k1__		= _k_aux__;
				//sums the old dtime - the variable dtime is alreaady updated
				this->time_new += this->previous_dt;
				//it steps the method ahead, with euler solution
				for(int i=0;i<numEDO;i++){
					this->setVariables(i, edos_new_euler_[i]);
				}
			}
		}
		if(savingRate!=0){
		    printf("Dt max: %e dt min %e, %e %d %d\n", maxDt, minDt, _soma_/_cont_, _cont_, cont_inv);
		}
		printf("%d",desc);
	}
	void Solveode::addt2(double finalTime, FILE *fileptr){
		const double _beta_safety_ = 0.8;
		double _tolerances_[numEDO];
		double _aux_tol=0.0;
		double maxDt = this->dtime, minDt = this->dtime;
		int desc=0;
		//initializes the variables
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->Cai_new_ = this->Cai_old_ = this->Cai_ini_;
		this->Cass_new_ = this->Cass_old_ = this->Cass_ini_;
		this->CaJSR_new_ = this->CaJSR_old_ = this->CaJSR_ini_;
		this->CaNSR_new_ = this->CaNSR_old_ = this->CaNSR_ini_;
		this->P_RyR_new_ = this->P_RyR_old_ = this->P_RyR_ini_;
		this->LTRPN_Ca_new_ = this->LTRPN_Ca_old_ = this->LTRPN_Ca_ini_;
		this->HTRPN_Ca_new_ = this->HTRPN_Ca_old_ = this->HTRPN_Ca_ini_;
		this->P_O1_new_ = this->P_O1_old_ = this->P_O1_ini_;
		this->P_O2_new_ = this->P_O2_old_ = this->P_O2_ini_;
		this->P_C2_new_ = this->P_C2_old_ = this->P_C2_ini_;
		this->O_new_ = this->O_old_ = this->O_ini_;
		this->C2_new_ = this->C2_old_ = this->C2_ini_;
		this->C3_new_ = this->C3_old_ = this->C3_ini_;
		this->C4_new_ = this->C4_old_ = this->C4_ini_;
		this->I1_new_ = this->I1_old_ = this->I1_ini_;
		this->I2_new_ = this->I2_old_ = this->I2_ini_;
		this->I3_new_ = this->I3_old_ = this->I3_ini_;
		this->Nai_new_ = this->Nai_old_ = this->Nai_ini_;
		this->C_Na2_new_ = this->C_Na2_old_ = this->C_Na2_ini_;
		this->C_Na1_new_ = this->C_Na1_old_ = this->C_Na1_ini_;
		this->O_Na_new_ = this->O_Na_old_ = this->O_Na_ini_;
		this->IF_Na_new_ = this->IF_Na_old_ = this->IF_Na_ini_;
		this->I1_Na_new_ = this->I1_Na_old_ = this->I1_Na_ini_;
		this->I2_Na_new_ = this->I2_Na_old_ = this->I2_Na_ini_;
		this->IC_Na2_new_ = this->IC_Na2_old_ = this->IC_Na2_ini_;
		this->IC_Na3_new_ = this->IC_Na3_old_ = this->IC_Na3_ini_;
		this->Ki_new_ = this->Ki_old_ = this->Ki_ini_;
		this->ato_f_new_ = this->ato_f_old_ = this->ato_f_ini_;
		this->ito_f_new_ = this->ito_f_old_ = this->ito_f_ini_;
		this->ato_s_new_ = this->ato_s_old_ = this->ato_s_ini_;
		this->ito_s_new_ = this->ito_s_old_ = this->ito_s_ini_;
		this->nKs_new_ = this->nKs_old_ = this->nKs_ini_;
		this->aur_new_ = this->aur_old_ = this->aur_ini_;
		this->iur_new_ = this->iur_old_ = this->iur_ini_;
		this->aKss_new_ = this->aKss_old_ = this->aKss_ini_;
		this->iKss_new_ = this->iKss_old_ = this->iKss_ini_;
		this->C_K2_new_ = this->C_K2_old_ = this->C_K2_ini_;
		this->C_K1_new_ = this->C_K1_old_ = this->C_K1_ini_;
		this->O_K_new_ = this->O_K_old_ = this->O_K_ini_;
		this->I_K_new_ = this->I_K_old_ = this->I_K_ini_;
		this->time_new = this->time;
		this->timeSaving = this->time;
		if(savingRate!=0)
			save_step(fileptr, _ADAP_DT_);//save the initial conditions
		double edos_old_aux_[numEDO];
		double edos_new_rk2_[numEDO],edos_new_euler_[numEDO];
		double *_k1__ = (double*)malloc(sizeof(double)*numEDO);
		double *_k2__ = (double*)malloc(sizeof(double)*numEDO);
		double *_k_aux__;
		int iMaiorErro=0;
		int _cont_=0,cont_inv=0;
		double _soma_=0.0;
		this->time_new += this->dtime;
		rightHandSideFunction.function(this);
		for(int i=0;i<numEDO;i++){
			_k1__[i] = this->getLadoDireito(i);
		}
		const double __tiny_ = pow(abstol__, 2.0);
		while(this->time_new<=finalTime){
			for(int i=0; i<numEDO; i++){
				//stores the old variables in a vector
				edos_old_aux_[i] = this->getVariables(i);
				//computes euler method
				edos_new_euler_[i] = _k1__[i] * this->dtime + edos_old_aux_[i];
				//steps ahead to compute the rk2 method
				this->setVariables(i, edos_new_euler_[i]);
			}
			//steps time ahead
			this->time_new += this->dtime;
			//computes the right-hand side one step ahead
			rightHandSideFunction.function(this);
			//restore the original old value
			this->time_new -= this->dtime;//step back
			double greatestError=0.0, auxError=0.0;
			for(int i=0;i<numEDO;i++){
				//stores the new evaluation 
				_k2__[i] = this->getLadoDireito(i);
				_aux_tol = fabs(edos_new_euler_[i])*reltol__;
				_tolerances_[i] = (abstol__ > _aux_tol )?abstol__:_aux_tol;
				//finds the greatest error between  the steps
				auxError = fabs(( (this->dtime/2)*(_k1__[i] - _k2__[i])) / _tolerances_[i]);
				if(auxError > greatestError){
					greatestError = auxError;
					iMaiorErro = i;
				}
			}
			///adapt the time step
			greatestError += __tiny_;
			this->previous_dt = this->dtime;
			///adapt the time step
			this->dtime = _beta_safety_ * this->dtime * sqrt(1.0/greatestError);
			//it doesn't accept the solution
			if(greatestError>=1.0){
			    cont_inv++;
				//restore the old values to do it again
				for(int i=0; i<numEDO; i++)
					this->setVariables(i, edos_old_aux_[i]);
				//throw the results away and compute again
				}else{//it accepts the solutions
					if(savingRate!=0.0){
					  //restore the previous value of old variables
					  for(int i=0; i<numEDO; i++)
						  this->setVariables(i, edos_old_aux_[i]);
					  this->V_new_ = edos_new_euler_[0];
					  this->Cai_new_ = edos_new_euler_[1];
					  this->Cass_new_ = edos_new_euler_[2];
					  this->CaJSR_new_ = edos_new_euler_[3];
					  this->CaNSR_new_ = edos_new_euler_[4];
					  this->P_RyR_new_ = edos_new_euler_[5];
					  this->LTRPN_Ca_new_ = edos_new_euler_[6];
					  this->HTRPN_Ca_new_ = edos_new_euler_[7];
					  this->P_O1_new_ = edos_new_euler_[8];
					  this->P_O2_new_ = edos_new_euler_[9];
					  this->P_C2_new_ = edos_new_euler_[10];
					  this->O_new_ = edos_new_euler_[11];
					  this->C2_new_ = edos_new_euler_[12];
					  this->C3_new_ = edos_new_euler_[13];
					  this->C4_new_ = edos_new_euler_[14];
					  this->I1_new_ = edos_new_euler_[15];
					  this->I2_new_ = edos_new_euler_[16];
					  this->I3_new_ = edos_new_euler_[17];
					  this->Nai_new_ = edos_new_euler_[18];
					  this->C_Na2_new_ = edos_new_euler_[19];
					  this->C_Na1_new_ = edos_new_euler_[20];
					  this->O_Na_new_ = edos_new_euler_[21];
					  this->IF_Na_new_ = edos_new_euler_[22];
					  this->I1_Na_new_ = edos_new_euler_[23];
					  this->I2_Na_new_ = edos_new_euler_[24];
					  this->IC_Na2_new_ = edos_new_euler_[25];
					  this->IC_Na3_new_ = edos_new_euler_[26];
					  this->Ki_new_ = edos_new_euler_[27];
					  this->ato_f_new_ = edos_new_euler_[28];
					  this->ito_f_new_ = edos_new_euler_[29];
					  this->ato_s_new_ = edos_new_euler_[30];
					  this->ito_s_new_ = edos_new_euler_[31];
					  this->nKs_new_ = edos_new_euler_[32];
					  this->aur_new_ = edos_new_euler_[33];
					  this->iur_new_ = edos_new_euler_[34];
					  this->aKss_new_ = edos_new_euler_[35];
					  this->iKss_new_ = edos_new_euler_[36];
					  this->C_K2_new_ = edos_new_euler_[37];
					  this->C_K1_new_ = edos_new_euler_[38];
					  this->O_K_new_ = edos_new_euler_[39];
					  this->I_K_new_ = edos_new_euler_[40];
					  save_step(fileptr, _ADAP_DT_);
					  if(this->dtime>maxDt) maxDt = this->dtime;
					  if(this->dtime<minDt) minDt = this->dtime;
					_soma_+=this->dtime;
					_cont_++;
				}
				
				if(this->dtime > maxStep && maxStep!=0){
				    
					this->dtime = maxStep;
				}else if(this->dtime==0){
					printf("Error: Time step is zero.\n");
					return;
				}
				//change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração
				_k_aux__	= _k2__;
				_k2__	= _k1__;
				_k1__	= _k_aux__;
				//sums the old dtime - the variable dtime is alreaady updated
				this->time_new += this->previous_dt;
				//it steps the method ahead, with euler solution
				for(int i=0;i<numEDO;i++){
					this->setVariables(i, edos_new_euler_[i]);
				}
			}
		}
// 		if(savingRate!=0){
		    printf("Dt max: %e dt min %e, %e %d %d\n", maxDt, minDt, _soma_/_cont_ , _cont_,cont_inv);
// 		}
	}
	void Solveode::euler_OMP(double finalTime, FILE *fileptr, int nThreads){
		omp_set_num_threads(nThreads);
		double *__NEW_ = (double*)malloc(sizeof(double)*numEDO);
		double *__OLD_ = (double*)malloc(sizeof(double)*numEDO);
// 		int tempoespera;
		double *temp;
		this->time_new = this->time;
		this->timeSaving = this->time;
		this->previous_dt = this->dtime;
				
		
// 		tempoespera=0;
		
		#pragma omp parallel firstprivate(__OLD_, __NEW_, temp)
// 		shared(tempoespera)
		{
		    
// 			tempoespera=0;
		    
			//private parameters
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -8.0000000000e+01,  _prvt_stim_start = 2.0000000000e+01,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 7.1430000000e+01,  _prvt_stim_duration = 1.5000000000e+00,  _prvt_Acap = 1.5340000000e-04,  _prvt_Cm = 1.0000000000e+00,  _prvt_Vmyo = 2.5840000000e-05,  _prvt_F = 9.6500000000e+01,  _prvt_VJSR = 1.2000000000e-07,  _prvt_Vss = 1.4850000000e-09,  _prvt_VNSR = 2.0980000000e-06,  _prvt_CMDN_tot = 5.0000000000e+01,  _prvt_Km_CMDN = 2.3800000000e-01,  _prvt_CSQN_tot = 1.5000000000e+04,  _prvt_Km_CSQN = 8.0000000000e+02,  _prvt_v1 = 4.5000000000e+00,  _prvt_tau_tr = 2.0000000000e+01,  _prvt_tau_xfer = 8.0000000000e+00,  _prvt_v2 = 1.7400000000e-05,  _prvt_v3 = 4.5000000000e-01,  _prvt_Km_up = 5.0000000000e-01,  _prvt_k_plus_htrpn = 2.3700000000e-03,  _prvt_HTRPN_tot = 1.4000000000e+02,  _prvt_k_plus_ltrpn = 3.2700000000e-02,  _prvt_LTRPN_tot = 7.0000000000e+01,  _prvt_k_minus_htrpn = 3.2000000000e-05,  _prvt_k_minus_ltrpn = 1.9600000000e-02,  _prvt_i_CaL_max = 7.0000000000e+00,  _prvt_k_plus_a = 6.0750000000e-03,  _prvt_n = 4.0000000000e+00,  _prvt_k_minus_b = 9.6500000000e-01,  _prvt_k_minus_c = 8.0000000000e-04,  _prvt_k_minus_a = 7.1250000000e-02,  _prvt_k_plus_b = 4.0500000000e-03,  _prvt_m = 3.0000000000e+00,  _prvt_k_plus_c = 9.0000000000e-03,  _prvt_g_CaL = 1.7290000000e-01,  _prvt_E_CaL = 6.3000000000e+01,  _prvt_Kpcb = 5.0000000000e-04,  _prvt_Kpc_max = 2.3324000000e-01,  _prvt_Kpc_half = 2.0000000000e+01,  _prvt_i_pCa_max = 1.0000000000e+00,  _prvt_Km_pCa = 5.0000000000e-01,  _prvt_k_NaCa = 2.9280000000e+02,  _prvt_K_mNa = 8.7500000000e+04,  _prvt_Nao = 1.4000000000e+05,  _prvt_K_mCa = 1.3800000000e+03,  _prvt_Cao = 1.8000000000e+03,  _prvt_k_sat = 1.0000000000e-01,  _prvt_eta = 3.5000000000e-01,  _prvt_R = 8.3140000000e+00,  _prvt_T = 2.9800000000e+02,  _prvt_g_Cab = 3.6700000000e-04,  _prvt_g_Na = 1.3000000000e+01,  _prvt_Ko = 5.4000000000e+03,  _prvt_g_Nab = 2.6000000000e-03,  _prvt_g_Kto_f = 4.0670000000e-01,  _prvt_g_Kto_s = 0.0000000000e+00,  _prvt_g_Ks = 5.7500000000e-03,  _prvt_g_Kur = 1.6000000000e-01,  _prvt_g_Kss = 5.0000000000e-02,  _prvt_g_Kr = 7.8000000000e-02,  _prvt_kf = 2.3761000000e-02,  _prvt_kb = 3.6778000000e-02,  _prvt_i_NaK_max = 8.8000000000e-01,  _prvt_Km_Nai = 2.1000000000e+04,  _prvt_Km_Ko = 1.5000000000e+03,  _prvt_g_ClCa = 1.0000000000e+01,  _prvt_Km_Cl = 1.0000000000e+01,  _prvt_E_Cl = -4.0000000000e+01, 
			//private aux variables
			 _prvt_calc_i_stim=0,  _prvt_calc_Bi=0,  _prvt_calc_Bss=0,  _prvt_calc_BJSR=0,  _prvt_calc_J_rel=0,  _prvt_calc_J_tr=0,  _prvt_calc_J_xfer=0,  _prvt_calc_J_leak=0,  _prvt_calc_J_up=0,  _prvt_calc_J_trpn=0,  _prvt_calc_P_C1=0,  _prvt_calc_i_CaL=0,  _prvt_calc_C1=0,  _prvt_calc_alpha=0,  _prvt_calc_beta=0,  _prvt_calc_gamma=0,  _prvt_calc_Kpcf=0,  _prvt_calc_i_pCa=0,  _prvt_calc_i_NaCa=0,  _prvt_calc_i_Cab=0,  _prvt_calc_E_CaN=0,  _prvt_calc_i_Na=0,  _prvt_calc_E_Na=0,  _prvt_calc_C_Na3=0,  _prvt_calc_alpha_Na11=0,  _prvt_calc_alpha_Na12=0,  _prvt_calc_alpha_Na13=0,  _prvt_calc_beta_Na11=0,  _prvt_calc_beta_Na12=0,  _prvt_calc_beta_Na13=0,  _prvt_calc_alpha_Na3=0,  _prvt_calc_beta_Na3=0,  _prvt_calc_alpha_Na2=0,  _prvt_calc_beta_Na2=0,  _prvt_calc_alpha_Na4=0,  _prvt_calc_beta_Na4=0,  _prvt_calc_alpha_Na5=0,  _prvt_calc_beta_Na5=0,  _prvt_calc_i_Nab=0,  _prvt_calc_i_Kto_f=0,  _prvt_calc_E_K=0,  _prvt_calc_alpha_a=0,  _prvt_calc_beta_a=0,  _prvt_calc_alpha_i=0,  _prvt_calc_beta_i=0,  _prvt_calc_i_Kto_s=0,  _prvt_calc_ass=0,  _prvt_calc_iss=0,  _prvt_calc_tau_ta_s=0,  _prvt_calc_tau_ti_s=0,  _prvt_calc_i_K1=0,  _prvt_calc_i_Ks=0,  _prvt_calc_alpha_n=0,  _prvt_calc_beta_n=0,  _prvt_calc_i_Kur=0,  _prvt_calc_tau_aur=0,  _prvt_calc_tau_iur=0,  _prvt_calc_i_Kss=0,  _prvt_calc_tau_Kss=0,  _prvt_calc_i_Kr=0,  _prvt_calc_C_K0=0,  _prvt_calc_alpha_a0=0,  _prvt_calc_beta_a0=0,  _prvt_calc_alpha_a1=0,  _prvt_calc_beta_a1=0,  _prvt_calc_i_NaK=0,  _prvt_calc_f_NaK=0,  _prvt_calc_sigma=0,  _prvt_calc_i_ClCa=0,  _prvt_calc_O_ClCa=0,  _prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current=0,  _prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current=0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_Cai_lado_direito_,  _prvt_Cass_lado_direito_,  _prvt_CaJSR_lado_direito_,  _prvt_CaNSR_lado_direito_,  _prvt_P_RyR_lado_direito_,  _prvt_LTRPN_Ca_lado_direito_,  _prvt_HTRPN_Ca_lado_direito_,  _prvt_P_O1_lado_direito_,  _prvt_P_O2_lado_direito_,  _prvt_P_C2_lado_direito_,  _prvt_O_lado_direito_,  _prvt_C2_lado_direito_,  _prvt_C3_lado_direito_,  _prvt_C4_lado_direito_,  _prvt_I1_lado_direito_,  _prvt_I2_lado_direito_,  _prvt_I3_lado_direito_,  _prvt_Nai_lado_direito_,  _prvt_C_Na2_lado_direito_,  _prvt_C_Na1_lado_direito_,  _prvt_O_Na_lado_direito_,  _prvt_IF_Na_lado_direito_,  _prvt_I1_Na_lado_direito_,  _prvt_I2_Na_lado_direito_,  _prvt_IC_Na2_lado_direito_,  _prvt_IC_Na3_lado_direito_,  _prvt_Ki_lado_direito_,  _prvt_ato_f_lado_direito_,  _prvt_ito_f_lado_direito_,  _prvt_ato_s_lado_direito_,  _prvt_ito_s_lado_direito_,  _prvt_nKs_lado_direito_,  _prvt_aur_lado_direito_,  _prvt_iur_lado_direito_,  _prvt_aKss_lado_direito_,  _prvt_iKss_lado_direito_,  _prvt_C_K2_lado_direito_,  _prvt_C_K1_lado_direito_,  _prvt_O_K_lado_direito_,  _prvt_I_K_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate;
			__NEW_[0] = __OLD_[0] = -8.2420200000e+01;
			__NEW_[1] = __OLD_[1] = 1.1500100000e-01;
			__NEW_[2] = __OLD_[2] = 1.1500100000e-01;
			__NEW_[3] = __OLD_[3] = 1.2995000000e+03;
			__NEW_[4] = __OLD_[4] = 1.2995000000e+03;
			__NEW_[5] = __OLD_[5] = 0.0000000000e+00;
			__NEW_[6] = __OLD_[6] = 1.1268400000e+01;
			__NEW_[7] = __OLD_[7] = 1.2529000000e+02;
			__NEW_[8] = __OLD_[8] = 1.4910200000e-05;
			__NEW_[9] = __OLD_[9] = 9.5172600000e-11;
			__NEW_[10] = __OLD_[10] = 1.6774000000e-04;
			__NEW_[11] = __OLD_[11] = 9.3030800000e-19;
			__NEW_[12] = __OLD_[12] = 1.2421600000e-04;
			__NEW_[13] = __OLD_[13] = 5.7867900000e-09;
			__NEW_[14] = __OLD_[14] = 1.1981600000e-13;
			__NEW_[15] = __OLD_[15] = 4.9792300000e-19;
			__NEW_[16] = __OLD_[16] = 3.4584700000e-14;
			__NEW_[17] = __OLD_[17] = 1.8510600000e-14;
			__NEW_[18] = __OLD_[18] = 1.4237100000e+04;
			__NEW_[19] = __OLD_[19] = 2.0752000000e-02;
			__NEW_[20] = __OLD_[20] = 2.7913200000e-04;
			__NEW_[21] = __OLD_[21] = 7.1348300000e-07;
			__NEW_[22] = __OLD_[22] = 1.5317600000e-04;
			__NEW_[23] = __OLD_[23] = 6.7334500000e-07;
			__NEW_[24] = __OLD_[24] = 1.5578700000e-09;
			__NEW_[25] = __OLD_[25] = 1.1387900000e-02;
			__NEW_[26] = __OLD_[26] = 3.4278000000e-01;
			__NEW_[27] = __OLD_[27] = 1.4372000000e+05;
			__NEW_[28] = __OLD_[28] = 2.6556300000e-03;
			__NEW_[29] = __OLD_[29] = 9.9997700000e-01;
			__NEW_[30] = __OLD_[30] = 4.1706900000e-04;
			__NEW_[31] = __OLD_[31] = 9.9854300000e-01;
			__NEW_[32] = __OLD_[32] = 2.6275300000e-04;
			__NEW_[33] = __OLD_[33] = 4.1706900000e-04;
			__NEW_[34] = __OLD_[34] = 9.9854300000e-01;
			__NEW_[35] = __OLD_[35] = 4.1706900000e-04;
			__NEW_[36] = __OLD_[36] = 1.0000000000e+00;
			__NEW_[37] = __OLD_[37] = 6.4122900000e-04;
			__NEW_[38] = __OLD_[38] = 9.9251300000e-04;
			__NEW_[39] = __OLD_[39] = 1.7529800000e-04;
			__NEW_[40] = __OLD_[40] = 3.1912900000e-05;
			int *_prvt_tree_thread = tree_thread;
			while(_prvt_time_new<=_prvt_finalTime){
				_prvt_time_new += _prvt_dtime;
				if(omp_get_thread_num()==tree_thread[0])
				{
					_prvt_calc_i_stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
					_prvt_calc_Bi = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[1]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_Bss = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[2]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_BJSR = pow((1.0000000000e+00+((_prvt_CSQN_tot*_prvt_Km_CSQN)/pow((_prvt_Km_CSQN+__OLD_[3]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_J_rel = (_prvt_v1*(__OLD_[8]+__OLD_[9])*(__OLD_[3]-__OLD_[2])*__OLD_[5]);
					_prvt_calc_J_tr = ((__OLD_[4]-__OLD_[3])/_prvt_tau_tr);
					_prvt_calc_J_xfer = ((__OLD_[2]-__OLD_[1])/_prvt_tau_xfer);
					_prvt_calc_J_leak = (_prvt_v2*(__OLD_[4]-__OLD_[1]));
					_prvt_calc_J_up = ((_prvt_v3*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_up,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
					_prvt_calc_J_trpn = (((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))+(_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6])))-((_prvt_k_minus_htrpn*__OLD_[7])+(_prvt_k_minus_ltrpn*__OLD_[6])));
					_prvt_calc_i_CaL = (_prvt_g_CaL*__OLD_[11]*(__OLD_[0]-_prvt_E_CaL));
					_prvt_calc_i_pCa = ((_prvt_i_pCa_max*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_pCa,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
					_prvt_calc_i_NaCa = (((((((_prvt_k_NaCa*1.0000000000e+00)/(pow(_prvt_K_mNa,3.0000000000e+00)+pow(_prvt_Nao,3.0000000000e+00)))*1.0000000000e+00)/(_prvt_K_mCa+_prvt_Cao))*1.0000000000e+00)/(1.0000000000e+00+(_prvt_k_sat*exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))))*((exp(((_prvt_eta*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[18],3.0000000000e+00)*_prvt_Cao)-(exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Nao,3.0000000000e+00)*__OLD_[1])));
					_prvt_calc_E_CaN = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Cao/__OLD_[1])));
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((((9.0000000000e-01*_prvt_Nao)+(1.0000000000e-01*_prvt_Ko))/((9.0000000000e-01*__OLD_[18])+(1.0000000000e-01*__OLD_[27])))));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ko/__OLD_[27])));
					_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[39]*(__OLD_[0]-(((_prvt_R*_prvt_T)/_prvt_F)*log((((9.8000000000e-01*_prvt_Ko)+(2.0000000000e-02*_prvt_Nao))/((9.8000000000e-01*__OLD_[27])+(2.0000000000e-02*__OLD_[18])))))));
					_prvt_calc_sigma = ((1.0000000000e+00/7.0000000000e+00)*(exp((_prvt_Nao/6.7300000000e+04))-1.0000000000e+00));
					_prvt_calc_O_ClCa = (2.0000000000e-01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.6700000000e+01))/7.8000000000e+00))));
					_prvt_calc_i_Nab = (_prvt_g_Nab*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Kto_s = (_prvt_g_Kto_s*__OLD_[30]*__OLD_[31]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K1 = ((((2.9380000000e-01*_prvt_Ko)/(_prvt_Ko+2.1000000000e+02))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp((8.9600000000e-02*(__OLD_[0]-_prvt_calc_E_K)))));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[32],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Kur = (_prvt_g_Kur*__OLD_[33]*__OLD_[34]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Kss = (_prvt_g_Kss*__OLD_[35]*__OLD_[36]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Cab = (_prvt_g_Cab*(__OLD_[0]-_prvt_calc_E_CaN));
					_prvt_calc_i_Na = (_prvt_g_Na*__OLD_[21]*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Kto_f = (_prvt_g_Kto_f*pow(__OLD_[28],3.0000000000e+00)*__OLD_[29]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_f_NaK = (1.0000000000e+00/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.6500000000e-02*_prvt_calc_sigma*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_ClCa = (((_prvt_g_ClCa*_prvt_calc_O_ClCa*__OLD_[1])/(__OLD_[1]+_prvt_Km_Cl))*(__OLD_[0]-_prvt_E_Cl));
					_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_calc_f_NaK*1.0000000000e+00)/(1.0000000000e+00+pow((_prvt_Km_Nai/__OLD_[18]),1.5000000000e+00)))*_prvt_Ko)/(_prvt_Ko+_prvt_Km_Ko));
					_prvt_V_lado_direito_= (-(_prvt_calc_i_CaL+_prvt_calc_i_pCa+_prvt_calc_i_NaCa+_prvt_calc_i_Cab+_prvt_calc_i_Na+_prvt_calc_i_Nab+_prvt_calc_i_NaK+_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kur+_prvt_calc_i_Kss+_prvt_calc_i_Kr+_prvt_calc_i_ClCa+_prvt_calc_i_stim));
					__NEW_[0]= _prvt_V_lado_direito_ * _prvt_dtime + __OLD_[0];
					_prvt_Cai_lado_direito_= (_prvt_calc_Bi*((_prvt_calc_J_leak+_prvt_calc_J_xfer)-(_prvt_calc_J_up+_prvt_calc_J_trpn+((((_prvt_calc_i_Cab+_prvt_calc_i_pCa)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vmyo*_prvt_F)))));
					__NEW_[1]= _prvt_Cai_lado_direito_ * _prvt_dtime + __OLD_[1];
					_prvt_Cass_lado_direito_= (_prvt_calc_Bss*(((_prvt_calc_J_rel*_prvt_VJSR)/_prvt_Vss)-(((_prvt_calc_J_xfer*_prvt_Vmyo)/_prvt_Vss)+((_prvt_calc_i_CaL*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vss*_prvt_F)))));
					__NEW_[2]= _prvt_Cass_lado_direito_ * _prvt_dtime + __OLD_[2];
					_prvt_CaJSR_lado_direito_= (_prvt_calc_BJSR*(_prvt_calc_J_tr-_prvt_calc_J_rel));
					__NEW_[3]= _prvt_CaJSR_lado_direito_ * _prvt_dtime + __OLD_[3];
					_prvt_CaNSR_lado_direito_= ((((_prvt_calc_J_up-_prvt_calc_J_leak)*_prvt_Vmyo)/_prvt_VNSR)-((_prvt_calc_J_tr*_prvt_VJSR)/_prvt_VNSR));
					__NEW_[4]= _prvt_CaNSR_lado_direito_ * _prvt_dtime + __OLD_[4];
					_prvt_P_RyR_lado_direito_= (((-4.0000000000e-02)*__OLD_[5])-(((1.0000000000e-01*_prvt_calc_i_CaL)/_prvt_i_CaL_max)*exp(((-pow((__OLD_[0]-5.0000000000e+00),2.0000000000e+00))/6.4800000000e+02))));
					__NEW_[5]= _prvt_P_RyR_lado_direito_ * _prvt_dtime + __OLD_[5];
					_prvt_Nai_lado_direito_= (((-(_prvt_calc_i_Na+_prvt_calc_i_Nab+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
					__NEW_[18]= _prvt_Nai_lado_direito_ * _prvt_dtime + __OLD_[18];
					_prvt_Ki_lado_direito_= (((-((_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kss+_prvt_calc_i_Kur+_prvt_calc_i_Kr)-(2.0000000000e+00*_prvt_calc_i_NaK)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
					__NEW_[27]= _prvt_Ki_lado_direito_ * _prvt_dtime + __OLD_[27];
				}
				if(omp_get_thread_num()==tree_thread[1])
				{
					_prvt_calc_P_C1 = (1.0000000000e+00-(__OLD_[10]+__OLD_[8]+__OLD_[9]));
					_prvt_P_O1_lado_direito_= (((_prvt_k_plus_a*pow(__OLD_[2],_prvt_n)*_prvt_calc_P_C1)+(_prvt_k_minus_b*__OLD_[9])+(_prvt_k_minus_c*__OLD_[10]))-((_prvt_k_minus_a*__OLD_[8])+(_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])+(_prvt_k_plus_c*__OLD_[8])));
					__NEW_[8]= _prvt_P_O1_lado_direito_ * _prvt_dtime + __OLD_[8];
				}
				if(omp_get_thread_num()==tree_thread[2])
				{
					_prvt_calc_C1 = (1.0000000000e+00-(__OLD_[11]+__OLD_[12]+__OLD_[12]+__OLD_[13]+__OLD_[14]+__OLD_[15]+__OLD_[16]+__OLD_[17]));
					_prvt_calc_alpha = ((4.0000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01))*((1.0000000000e+00+(7.0000000000e-01*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.0000000000e+01))))-(7.5000000000e-01*exp(((-pow((__OLD_[0]+2.0000000000e+01),2.0000000000e+00))/4.0000000000e+02)))))/(1.0000000000e+00+(1.2000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01)))));
					_prvt_calc_beta = (5.0000000000e-02*exp(((-(__OLD_[0]+1.2000000000e+01))/1.3000000000e+01)));
					_prvt_calc_gamma = ((_prvt_Kpc_max*__OLD_[2])/(_prvt_Kpc_half+__OLD_[2]));
					_prvt_calc_Kpcf = (1.3000000000e+01*(1.0000000000e+00-exp(((-pow((__OLD_[0]+1.4500000000e+01),2.0000000000e+00))/1.0000000000e+02))));
					_prvt_O_lado_direito_= (((_prvt_calc_alpha*__OLD_[14])+(_prvt_Kpcb*__OLD_[15])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[11]))))-((4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(_prvt_calc_gamma*__OLD_[11])));
					__NEW_[11]= _prvt_O_lado_direito_ * _prvt_dtime + __OLD_[11];
					_prvt_C2_lado_direito_= (((4.0000000000e+00*_prvt_calc_alpha*_prvt_calc_C1)+(2.0000000000e+00*_prvt_calc_beta*__OLD_[13]))-((_prvt_calc_beta*__OLD_[12])+(3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])));
					__NEW_[12]= _prvt_C2_lado_direito_ * _prvt_dtime + __OLD_[12];
					_prvt_C3_lado_direito_= (((3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])+(3.0000000000e+00*_prvt_calc_beta*__OLD_[14]))-((2.0000000000e+00*_prvt_calc_beta*__OLD_[13])+(2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])));
					__NEW_[13]= _prvt_C3_lado_direito_ * _prvt_dtime + __OLD_[13];
					_prvt_C4_lado_direito_= (((2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])+(4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(1.0000000000e-02*((4.0000000000e+00*_prvt_Kpcb*_prvt_calc_beta*__OLD_[15])-(_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])))+(2.0000000000e-03*((4.0000000000e+00*_prvt_calc_beta*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[14])))+(4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17]))-((3.0000000000e+00*_prvt_calc_beta*__OLD_[14])+(_prvt_calc_alpha*__OLD_[14])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14])));
					__NEW_[14]= _prvt_C4_lado_direito_ * _prvt_dtime + __OLD_[14];
					_prvt_I1_lado_direito_= (((_prvt_calc_gamma*__OLD_[11])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[17])-(_prvt_calc_Kpcf*__OLD_[15])))+(1.0000000000e-02*((_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*_prvt_calc_Kpcf*__OLD_[15]))))-(_prvt_Kpcb*__OLD_[15]));
					__NEW_[15]= _prvt_I1_lado_direito_ * _prvt_dtime + __OLD_[15];
					_prvt_I2_lado_direito_= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[11])-(_prvt_calc_alpha*__OLD_[16])))+(_prvt_Kpcb*__OLD_[17])+(2.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*__OLD_[16]))))-(_prvt_calc_gamma*__OLD_[16]));
					__NEW_[16]= _prvt_I2_lado_direito_ * _prvt_dtime + __OLD_[16];
					_prvt_I3_lado_direito_= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[15])-(_prvt_calc_alpha*__OLD_[17])))+(_prvt_calc_gamma*__OLD_[16])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14]))-((4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17])+(_prvt_Kpcb*__OLD_[17])));
					__NEW_[17]= _prvt_I3_lado_direito_ * _prvt_dtime + __OLD_[17];
				}
				if(omp_get_thread_num()==tree_thread[3])
				{
					_prvt_calc_C_Na3 = (1.0000000000e+00-(__OLD_[21]+__OLD_[20]+__OLD_[19]+__OLD_[22]+__OLD_[23]+__OLD_[24]+__OLD_[25]+__OLD_[26]));
					_prvt_calc_alpha_Na11 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.7000000000e+01)))+(2.0000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_alpha_Na12 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+01)))+(2.3000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_alpha_Na13 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.2000000000e+01)))+(2.5000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_beta_Na11 = (1.9170000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_beta_Na12 = (2.0000000000e-01*exp(((-(__OLD_[0]-2.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_beta_Na13 = (2.2000000000e-01*exp(((-(__OLD_[0]-7.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_alpha_Na3 = (7.0000000000e-07*exp(((-(__OLD_[0]+7.0000000000e+00))/7.7000000000e+00)));
					_prvt_calc_beta_Na3 = (8.4000000000e-03+(2.0000000000e-05*(__OLD_[0]+7.0000000000e+00)));
					_prvt_calc_alpha_Na2 = (1.0000000000e+00/((1.8849500000e-01*exp(((-(__OLD_[0]+7.0000000000e+00))/1.6600000000e+01)))+3.9395600000e-01));
					_prvt_calc_beta_Na2 = ((_prvt_calc_alpha_Na13*_prvt_calc_alpha_Na2*_prvt_calc_alpha_Na3)/(_prvt_calc_beta_Na13*_prvt_calc_beta_Na3));
					_prvt_calc_alpha_Na4 = (_prvt_calc_alpha_Na2/1.0000000000e+03);
					_prvt_calc_beta_Na4 = _prvt_calc_alpha_Na3;
					_prvt_calc_alpha_Na5 = (_prvt_calc_alpha_Na2/9.5000000000e+04);
					_prvt_calc_beta_Na5 = (_prvt_calc_alpha_Na3/5.0000000000e+01);
					_prvt_C_Na2_lado_direito_= (((_prvt_calc_alpha_Na11*_prvt_calc_C_Na3)+(_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[19])+(_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na3*__OLD_[19])));
					__NEW_[19]= _prvt_C_Na2_lado_direito_ * _prvt_dtime + __OLD_[19];
					_prvt_C_Na1_lado_direito_= (((_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na3*__OLD_[22]))-((_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na3*__OLD_[20])));
					__NEW_[20]= _prvt_C_Na1_lado_direito_ * _prvt_dtime + __OLD_[20];
					_prvt_O_Na_lado_direito_= (((_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na2*__OLD_[22]))-((_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na2*__OLD_[21])));
					__NEW_[21]= _prvt_O_Na_lado_direito_ * _prvt_dtime + __OLD_[21];
					_prvt_IF_Na_lado_direito_= (((_prvt_calc_alpha_Na2*__OLD_[21])+(_prvt_calc_beta_Na3*__OLD_[20])+(_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na12*__OLD_[25]))-((_prvt_calc_beta_Na2*__OLD_[22])+(_prvt_calc_alpha_Na3*__OLD_[22])+(_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na12*__OLD_[22])));
					__NEW_[22]= _prvt_IF_Na_lado_direito_ * _prvt_dtime + __OLD_[22];
					_prvt_I1_Na_lado_direito_= (((_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na5*__OLD_[24]))-((_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na5*__OLD_[23])));
					__NEW_[23]= _prvt_I1_Na_lado_direito_ * _prvt_dtime + __OLD_[23];
					_prvt_I2_Na_lado_direito_= ((_prvt_calc_alpha_Na5*__OLD_[23])-(_prvt_calc_beta_Na5*__OLD_[24]));
					__NEW_[24]= _prvt_I2_Na_lado_direito_ * _prvt_dtime + __OLD_[24];
					_prvt_IC_Na2_lado_direito_= (((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_beta_Na12*__OLD_[22])+(_prvt_calc_beta_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_alpha_Na12*__OLD_[25])+(_prvt_calc_alpha_Na3*__OLD_[25])));
					__NEW_[25]= _prvt_IC_Na2_lado_direito_ * _prvt_dtime + __OLD_[25];
					_prvt_IC_Na3_lado_direito_= (((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_beta_Na3*_prvt_calc_C_Na3))-((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_alpha_Na3*__OLD_[26])));
					__NEW_[26]= _prvt_IC_Na3_lado_direito_ * _prvt_dtime + __OLD_[26];
				}
				if(omp_get_thread_num()==tree_thread[4])
				{
					_prvt_calc_alpha_a = (1.8064000000e-01*exp((3.5770000000e-02*(__OLD_[0]+3.0000000000e+01))));
					_prvt_calc_beta_a = (3.9560000000e-01*exp(((-6.2370000000e-02)*(__OLD_[0]+3.0000000000e+01))));
					_prvt_ato_f_lado_direito_= ((_prvt_calc_alpha_a*(1.0000000000e+00-__OLD_[28]))-(_prvt_calc_beta_a*__OLD_[28]));
					__NEW_[28]= _prvt_ato_f_lado_direito_ * _prvt_dtime + __OLD_[28];
				}
				if(omp_get_thread_num()==tree_thread[5])
				{
					_prvt_calc_alpha_i = ((1.5200000000e-04*exp(((-(__OLD_[0]+1.3500000000e+01))/7.0000000000e+00)))/((6.7083000000e-03*exp(((-(__OLD_[0]+3.3500000000e+01))/7.0000000000e+00)))+1.0000000000e+00));
					_prvt_calc_beta_i = ((9.5000000000e-04*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))/((5.1335000000e-02*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))+1.0000000000e+00));
					_prvt_ito_f_lado_direito_= ((_prvt_calc_alpha_i*(1.0000000000e+00-__OLD_[29]))-(_prvt_calc_beta_i*__OLD_[29]));
					__NEW_[29]= _prvt_ito_f_lado_direito_ * _prvt_dtime + __OLD_[29];
				}
				if(omp_get_thread_num()==tree_thread[6])
				{
					_prvt_calc_ass = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2500000000e+01))/7.7000000000e+00))));
					_prvt_calc_tau_ta_s = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
					_prvt_calc_tau_aur = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
					_prvt_calc_tau_Kss = ((3.9300000000e+01*exp(((-8.6200000000e-02)*__OLD_[0])))+1.3170000000e+01);
					_prvt_ato_s_lado_direito_= ((_prvt_calc_ass-__OLD_[30])/_prvt_calc_tau_ta_s);
					__NEW_[30]= _prvt_ato_s_lado_direito_ * _prvt_dtime + __OLD_[30];
					_prvt_aur_lado_direito_= ((_prvt_calc_ass-__OLD_[33])/_prvt_calc_tau_aur);
					__NEW_[33]= _prvt_aur_lado_direito_ * _prvt_dtime + __OLD_[33];
					_prvt_aKss_lado_direito_= ((_prvt_calc_ass-__OLD_[35])/_prvt_calc_tau_Kss);
					__NEW_[35]= _prvt_aKss_lado_direito_ * _prvt_dtime + __OLD_[35];
				}
				if(omp_get_thread_num()==tree_thread[7])
				{
					_prvt_calc_iss = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00))));
					_prvt_calc_tau_ti_s = (2.7000000000e+02+(1.0500000000e+03/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
					_prvt_calc_tau_iur = (1.2000000000e+03-(1.7000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
					_prvt_ito_s_lado_direito_= ((_prvt_calc_iss-__OLD_[31])/_prvt_calc_tau_ti_s);
					__NEW_[31]= _prvt_ito_s_lado_direito_ * _prvt_dtime + __OLD_[31];
					_prvt_iur_lado_direito_= ((_prvt_calc_iss-__OLD_[34])/_prvt_calc_tau_iur);
					__NEW_[34]= _prvt_iur_lado_direito_ * _prvt_dtime + __OLD_[34];
				}
				if(omp_get_thread_num()==tree_thread[8])
				{
					_prvt_calc_alpha_n = ((4.8133300000e-06*(__OLD_[0]+2.6500000000e+01))/(1.0000000000e+00-exp(((-1.2800000000e-01)*(__OLD_[0]+2.6500000000e+01)))));
					_prvt_calc_beta_n = (9.5333300000e-05*exp(((-3.8000000000e-02)*(__OLD_[0]+2.6500000000e+01))));
					_prvt_nKs_lado_direito_= ((_prvt_calc_alpha_n*(1.0000000000e+00-__OLD_[32]))-(_prvt_calc_beta_n*__OLD_[32]));
					__NEW_[32]= _prvt_nKs_lado_direito_ * _prvt_dtime + __OLD_[32];
				}
				if(omp_get_thread_num()==tree_thread[9])
				{
					_prvt_calc_C_K0 = (1.0000000000e+00-(__OLD_[38]+__OLD_[37]+__OLD_[39]+__OLD_[40]));
					_prvt_calc_alpha_a0 = (2.2348000000e-02*exp((1.1760000000e-02*__OLD_[0])));
					_prvt_calc_beta_a0 = (4.7002000000e-02*exp(((-6.3100000000e-02)*__OLD_[0])));
					_prvt_C_K1_lado_direito_= (((_prvt_calc_alpha_a0*_prvt_calc_C_K0)+(_prvt_kb*__OLD_[37]))-((_prvt_calc_beta_a0*__OLD_[38])+(_prvt_kf*__OLD_[38])));
					__NEW_[38]= _prvt_C_K1_lado_direito_ * _prvt_dtime + __OLD_[38];
				}
				if(omp_get_thread_num()==tree_thread[10])
				{
					_prvt_calc_alpha_a1 = (1.3733000000e-02*exp((3.8198000000e-02*__OLD_[0])));
					_prvt_calc_beta_a1 = (6.8900000000e-05*exp(((-4.1780000000e-02)*__OLD_[0])));
					_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (9.0821000000e-02*exp((2.3391000000e-02*(__OLD_[0]+5.0000000000e+00))));
					_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (6.4970000000e-03*exp(((-3.2680000000e-02)*(__OLD_[0]+5.0000000000e+00))));
					_prvt_C_K2_lado_direito_= (((_prvt_kf*__OLD_[38])+(_prvt_calc_beta_a1*__OLD_[39]))-((_prvt_kb*__OLD_[37])+(_prvt_calc_alpha_a1*__OLD_[37])));
					__NEW_[37]= _prvt_C_K2_lado_direito_ * _prvt_dtime + __OLD_[37];
					_prvt_O_K_lado_direito_= (((_prvt_calc_alpha_a1*__OLD_[37])+(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]))-((_prvt_calc_beta_a1*__OLD_[39])+(_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])));
					__NEW_[39]= _prvt_O_K_lado_direito_ * _prvt_dtime + __OLD_[39];
					_prvt_I_K_lado_direito_= ((_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])-(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]));
					__NEW_[40]= _prvt_I_K_lado_direito_ * _prvt_dtime + __OLD_[40];
				}
				if(omp_get_thread_num()==tree_thread[11])
				{
					_prvt_LTRPN_Ca_lado_direito_= ((_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6]))-(_prvt_k_minus_ltrpn*__OLD_[6]));
					__NEW_[6]= _prvt_LTRPN_Ca_lado_direito_ * _prvt_dtime + __OLD_[6];
				}
				if(omp_get_thread_num()==tree_thread[12])
				{
					_prvt_HTRPN_Ca_lado_direito_= ((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))-(_prvt_k_minus_htrpn*__OLD_[7]));
					__NEW_[7]= _prvt_HTRPN_Ca_lado_direito_ * _prvt_dtime + __OLD_[7];
				}
				if(omp_get_thread_num()==tree_thread[13])
				{
					_prvt_P_O2_lado_direito_= ((_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])-(_prvt_k_minus_b*__OLD_[9]));
					__NEW_[9]= _prvt_P_O2_lado_direito_ * _prvt_dtime + __OLD_[9];
				}
				if(omp_get_thread_num()==tree_thread[14])
				{
					_prvt_P_C2_lado_direito_= ((_prvt_k_plus_c*__OLD_[8])-(_prvt_k_minus_c*__OLD_[10]));
					__NEW_[10]= _prvt_P_C2_lado_direito_ * _prvt_dtime + __OLD_[10];
				}
				if(omp_get_thread_num()==tree_thread[15])
				{
					_prvt_iKss_lado_direito_= 0.0000000000e+00;
					__NEW_[36]= _prvt_iKss_lado_direito_ * _prvt_dtime + __OLD_[36];
				}
				//synchronizing all threads
// 				timeval t1, t2;
// 				gettimeofday(&t1, NULL);
				#pragma omp barrier
// 				gettimeofday(&t2, NULL);
// 				#pragma omp critical
// 				{
// 				    int aux = (int)((int)t2.tv_usec - (int)t1.tv_usec);
// 				    tempoespera += (aux<0)?0:aux;
// 				    
// 				}
				
			/*	if(_prvt_savingRate!=0){
					#pragma omp single
					{
						this->V_old_ = __OLD_[0];
						this->V_new_ = __NEW_[0];
						this->Cai_old_ = __OLD_[1];
						this->Cai_new_ = __NEW_[1];
						this->Cass_old_ = __OLD_[2];
						this->Cass_new_ = __NEW_[2];
						this->CaJSR_old_ = __OLD_[3];
						this->CaJSR_new_ = __NEW_[3];
						this->CaNSR_old_ = __OLD_[4];
						this->CaNSR_new_ = __NEW_[4];
						this->P_RyR_old_ = __OLD_[5];
						this->P_RyR_new_ = __NEW_[5];
						this->LTRPN_Ca_old_ = __OLD_[6];
						this->LTRPN_Ca_new_ = __NEW_[6];
						this->HTRPN_Ca_old_ = __OLD_[7];
						this->HTRPN_Ca_new_ = __NEW_[7];
						this->P_O1_old_ = __OLD_[8];
						this->P_O1_new_ = __NEW_[8];
						this->P_O2_old_ = __OLD_[9];
						this->P_O2_new_ = __NEW_[9];
						this->P_C2_old_ = __OLD_[10];
						this->P_C2_new_ = __NEW_[10];
						this->O_old_ = __OLD_[11];
						this->O_new_ = __NEW_[11];
						this->C2_old_ = __OLD_[12];
						this->C2_new_ = __NEW_[12];
						this->C3_old_ = __OLD_[13];
						this->C3_new_ = __NEW_[13];
						this->C4_old_ = __OLD_[14];
						this->C4_new_ = __NEW_[14];
						this->I1_old_ = __OLD_[15];
						this->I1_new_ = __NEW_[15];
						this->I2_old_ = __OLD_[16];
						this->I2_new_ = __NEW_[16];
						this->I3_old_ = __OLD_[17];
						this->I3_new_ = __NEW_[17];
						this->Nai_old_ = __OLD_[18];
						this->Nai_new_ = __NEW_[18];
						this->C_Na2_old_ = __OLD_[19];
						this->C_Na2_new_ = __NEW_[19];
						this->C_Na1_old_ = __OLD_[20];
						this->C_Na1_new_ = __NEW_[20];
						this->O_Na_old_ = __OLD_[21];
						this->O_Na_new_ = __NEW_[21];
						this->IF_Na_old_ = __OLD_[22];
						this->IF_Na_new_ = __NEW_[22];
						this->I1_Na_old_ = __OLD_[23];
						this->I1_Na_new_ = __NEW_[23];
						this->I2_Na_old_ = __OLD_[24];
						this->I2_Na_new_ = __NEW_[24];
						this->IC_Na2_old_ = __OLD_[25];
						this->IC_Na2_new_ = __NEW_[25];
						this->IC_Na3_old_ = __OLD_[26];
						this->IC_Na3_new_ = __NEW_[26];
						this->Ki_old_ = __OLD_[27];
						this->Ki_new_ = __NEW_[27];
						this->ato_f_old_ = __OLD_[28];
						this->ato_f_new_ = __NEW_[28];
						this->ito_f_old_ = __OLD_[29];
						this->ito_f_new_ = __NEW_[29];
						this->ato_s_old_ = __OLD_[30];
						this->ato_s_new_ = __NEW_[30];
						this->ito_s_old_ = __OLD_[31];
						this->ito_s_new_ = __NEW_[31];
						this->nKs_old_ = __OLD_[32];
						this->nKs_new_ = __NEW_[32];
						this->aur_old_ = __OLD_[33];
						this->aur_new_ = __NEW_[33];
						this->iur_old_ = __OLD_[34];
						this->iur_new_ = __NEW_[34];
						this->aKss_old_ = __OLD_[35];
						this->aKss_new_ = __NEW_[35];
						this->iKss_old_ = __OLD_[36];
						this->iKss_new_ = __NEW_[36];
						this->C_K2_old_ = __OLD_[37];
						this->C_K2_new_ = __NEW_[37];
						this->C_K1_old_ = __OLD_[38];
						this->C_K1_new_ = __NEW_[38];
						this->O_K_old_ = __OLD_[39];
						this->O_K_new_ = __NEW_[39];
						this->I_K_old_ = __OLD_[40];
						this->I_K_new_ = __NEW_[40];
						this->time_new = _prvt_time_new;
						save_step(fileptr, _EULER_);
					}
				}
				*/
				temp = __OLD_;
				__OLD_ = __NEW_;
				__NEW_= temp;
			}
			
// 			printf("%.1f\n", (float)tempoespera/(float)(nThreads));
		}
		
	}
	void Solveode::addt_OMP(double finalTime, FILE *fileptr, int nThreads){
		omp_set_num_threads(nThreads);
		
		//int *tempoespera = (int*)malloc(sizeof(int)*nThreads);
// 		int tempoespera =0;
		double *__NEW_ = (double*)malloc(sizeof(double)*numEDO);
		double *__OLD_ = (double*)malloc(sizeof(double)*numEDO);
		double *__OLD_AUX_ = (double*)malloc(sizeof(double)*numEDO);
		double *__TOL_ = (double*)malloc(sizeof(double)*numEDO);
		double *__ERROR_ = (double*)malloc(sizeof(double)*numEDO);
		double *__K1_  = (double*)malloc(sizeof(double)*numEDO);
		double *__K2_  = (double*)malloc(sizeof(double)*numEDO);
		double *__TEMP_;
		this->time_new = this->time;
		this->timeSaving = this->time;
		this->previous_dt = this->dtime;
		
		
		    
		#pragma omp parallel firstprivate(__NEW_, __OLD_, __TEMP_, __TOL_, __K1_,__K2_, __ERROR_,__OLD_AUX_)
// 		shared(tempoespera)
		{
// 		    tempoespera =0;
			int *_prvt_tree_thread = tree_thread;
			double _prvt_rel_tol_=reltol__, _prvt_abs_tol_ =abstol__, _prvt_aux_tol;
			const double __tiny_ = pow(_prvt_abs_tol_, 2.0);
			//private parameters
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -8.0000000000e+01,  _prvt_stim_start = 2.0000000000e+01,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 7.1430000000e+01,  _prvt_stim_duration = 1.5000000000e+00,  _prvt_Acap = 1.5340000000e-04,  _prvt_Cm = 1.0000000000e+00,  _prvt_Vmyo = 2.5840000000e-05,  _prvt_F = 9.6500000000e+01,  _prvt_VJSR = 1.2000000000e-07,  _prvt_Vss = 1.4850000000e-09,  _prvt_VNSR = 2.0980000000e-06,  _prvt_CMDN_tot = 5.0000000000e+01,  _prvt_Km_CMDN = 2.3800000000e-01,  _prvt_CSQN_tot = 1.5000000000e+04,  _prvt_Km_CSQN = 8.0000000000e+02,  _prvt_v1 = 4.5000000000e+00,  _prvt_tau_tr = 2.0000000000e+01,  _prvt_tau_xfer = 8.0000000000e+00,  _prvt_v2 = 1.7400000000e-05,  _prvt_v3 = 4.5000000000e-01,  _prvt_Km_up = 5.0000000000e-01,  _prvt_k_plus_htrpn = 2.3700000000e-03,  _prvt_HTRPN_tot = 1.4000000000e+02,  _prvt_k_plus_ltrpn = 3.2700000000e-02,  _prvt_LTRPN_tot = 7.0000000000e+01,  _prvt_k_minus_htrpn = 3.2000000000e-05,  _prvt_k_minus_ltrpn = 1.9600000000e-02,  _prvt_i_CaL_max = 7.0000000000e+00,  _prvt_k_plus_a = 6.0750000000e-03,  _prvt_n = 4.0000000000e+00,  _prvt_k_minus_b = 9.6500000000e-01,  _prvt_k_minus_c = 8.0000000000e-04,  _prvt_k_minus_a = 7.1250000000e-02,  _prvt_k_plus_b = 4.0500000000e-03,  _prvt_m = 3.0000000000e+00,  _prvt_k_plus_c = 9.0000000000e-03,  _prvt_g_CaL = 1.7290000000e-01,  _prvt_E_CaL = 6.3000000000e+01,  _prvt_Kpcb = 5.0000000000e-04,  _prvt_Kpc_max = 2.3324000000e-01,  _prvt_Kpc_half = 2.0000000000e+01,  _prvt_i_pCa_max = 1.0000000000e+00,  _prvt_Km_pCa = 5.0000000000e-01,  _prvt_k_NaCa = 2.9280000000e+02,  _prvt_K_mNa = 8.7500000000e+04,  _prvt_Nao = 1.4000000000e+05,  _prvt_K_mCa = 1.3800000000e+03,  _prvt_Cao = 1.8000000000e+03,  _prvt_k_sat = 1.0000000000e-01,  _prvt_eta = 3.5000000000e-01,  _prvt_R = 8.3140000000e+00,  _prvt_T = 2.9800000000e+02,  _prvt_g_Cab = 3.6700000000e-04,  _prvt_g_Na = 1.3000000000e+01,  _prvt_Ko = 5.4000000000e+03,  _prvt_g_Nab = 2.6000000000e-03,  _prvt_g_Kto_f = 4.0670000000e-01,  _prvt_g_Kto_s = 0.0000000000e+00,  _prvt_g_Ks = 5.7500000000e-03,  _prvt_g_Kur = 1.6000000000e-01,  _prvt_g_Kss = 5.0000000000e-02,  _prvt_g_Kr = 7.8000000000e-02,  _prvt_kf = 2.3761000000e-02,  _prvt_kb = 3.6778000000e-02,  _prvt_i_NaK_max = 8.8000000000e-01,  _prvt_Km_Nai = 2.1000000000e+04,  _prvt_Km_Ko = 1.5000000000e+03,  _prvt_g_ClCa = 1.0000000000e+01,  _prvt_Km_Cl = 1.0000000000e+01,  _prvt_E_Cl = -4.0000000000e+01, 
			//private aux variables
			 _prvt_calc_i_stim=0.0,  _prvt_calc_Bi=0.0,  _prvt_calc_Bss=0.0,  _prvt_calc_BJSR=0.0,  _prvt_calc_J_rel=0.0,  _prvt_calc_J_tr=0.0,  _prvt_calc_J_xfer=0.0,  _prvt_calc_J_leak=0.0,  _prvt_calc_J_up=0.0,  _prvt_calc_J_trpn=0.0,  _prvt_calc_P_C1=0.0,  _prvt_calc_i_CaL=0.0,  _prvt_calc_C1=0.0,  _prvt_calc_alpha=0.0,  _prvt_calc_beta=0.0,  _prvt_calc_gamma=0.0,  _prvt_calc_Kpcf=0.0,  _prvt_calc_i_pCa=0.0,  _prvt_calc_i_NaCa=0.0,  _prvt_calc_i_Cab=0.0,  _prvt_calc_E_CaN=0.0,  _prvt_calc_i_Na=0.0,  _prvt_calc_E_Na=0.0,  _prvt_calc_C_Na3=0.0,  _prvt_calc_alpha_Na11=0.0,  _prvt_calc_alpha_Na12=0.0,  _prvt_calc_alpha_Na13=0.0,  _prvt_calc_beta_Na11=0.0,  _prvt_calc_beta_Na12=0.0,  _prvt_calc_beta_Na13=0.0,  _prvt_calc_alpha_Na3=0.0,  _prvt_calc_beta_Na3=0.0,  _prvt_calc_alpha_Na2=0.0,  _prvt_calc_beta_Na2=0.0,  _prvt_calc_alpha_Na4=0.0,  _prvt_calc_beta_Na4=0.0,  _prvt_calc_alpha_Na5=0.0,  _prvt_calc_beta_Na5=0.0,  _prvt_calc_i_Nab=0.0,  _prvt_calc_i_Kto_f=0.0,  _prvt_calc_E_K=0.0,  _prvt_calc_alpha_a=0.0,  _prvt_calc_beta_a=0.0,  _prvt_calc_alpha_i=0.0,  _prvt_calc_beta_i=0.0,  _prvt_calc_i_Kto_s=0.0,  _prvt_calc_ass=0.0,  _prvt_calc_iss=0.0,  _prvt_calc_tau_ta_s=0.0,  _prvt_calc_tau_ti_s=0.0,  _prvt_calc_i_K1=0.0,  _prvt_calc_i_Ks=0.0,  _prvt_calc_alpha_n=0.0,  _prvt_calc_beta_n=0.0,  _prvt_calc_i_Kur=0.0,  _prvt_calc_tau_aur=0.0,  _prvt_calc_tau_iur=0.0,  _prvt_calc_i_Kss=0.0,  _prvt_calc_tau_Kss=0.0,  _prvt_calc_i_Kr=0.0,  _prvt_calc_C_K0=0.0,  _prvt_calc_alpha_a0=0.0,  _prvt_calc_beta_a0=0.0,  _prvt_calc_alpha_a1=0.0,  _prvt_calc_beta_a1=0.0,  _prvt_calc_i_NaK=0.0,  _prvt_calc_f_NaK=0.0,  _prvt_calc_sigma=0.0,  _prvt_calc_i_ClCa=0.0,  _prvt_calc_O_ClCa=0.0,  _prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current=0.0,  _prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current=0.0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_Cai_lado_direito_,  _prvt_Cass_lado_direito_,  _prvt_CaJSR_lado_direito_,  _prvt_CaNSR_lado_direito_,  _prvt_P_RyR_lado_direito_,  _prvt_LTRPN_Ca_lado_direito_,  _prvt_HTRPN_Ca_lado_direito_,  _prvt_P_O1_lado_direito_,  _prvt_P_O2_lado_direito_,  _prvt_P_C2_lado_direito_,  _prvt_O_lado_direito_,  _prvt_C2_lado_direito_,  _prvt_C3_lado_direito_,  _prvt_C4_lado_direito_,  _prvt_I1_lado_direito_,  _prvt_I2_lado_direito_,  _prvt_I3_lado_direito_,  _prvt_Nai_lado_direito_,  _prvt_C_Na2_lado_direito_,  _prvt_C_Na1_lado_direito_,  _prvt_O_Na_lado_direito_,  _prvt_IF_Na_lado_direito_,  _prvt_I1_Na_lado_direito_,  _prvt_I2_Na_lado_direito_,  _prvt_IC_Na2_lado_direito_,  _prvt_IC_Na3_lado_direito_,  _prvt_Ki_lado_direito_,  _prvt_ato_f_lado_direito_,  _prvt_ito_f_lado_direito_,  _prvt_ato_s_lado_direito_,  _prvt_ito_s_lado_direito_,  _prvt_nKs_lado_direito_,  _prvt_aur_lado_direito_,  _prvt_iur_lado_direito_,  _prvt_aKss_lado_direito_,  _prvt_iKss_lado_direito_,  _prvt_C_K2_lado_direito_,  _prvt_C_K1_lado_direito_,  _prvt_O_K_lado_direito_,  _prvt_I_K_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;
			__NEW_[0] = __OLD_[0] = -8.2420200000e+01;
			__NEW_[1] = __OLD_[1] = 1.1500100000e-01;
			__NEW_[2] = __OLD_[2] = 1.1500100000e-01;
			__NEW_[3] = __OLD_[3] = 1.2995000000e+03;
			__NEW_[4] = __OLD_[4] = 1.2995000000e+03;
			__NEW_[5] = __OLD_[5] = 0.0000000000e+00;
			__NEW_[6] = __OLD_[6] = 1.1268400000e+01;
			__NEW_[7] = __OLD_[7] = 1.2529000000e+02;
			__NEW_[8] = __OLD_[8] = 1.4910200000e-05;
			__NEW_[9] = __OLD_[9] = 9.5172600000e-11;
			__NEW_[10] = __OLD_[10] = 1.6774000000e-04;
			__NEW_[11] = __OLD_[11] = 9.3030800000e-19;
			__NEW_[12] = __OLD_[12] = 1.2421600000e-04;
			__NEW_[13] = __OLD_[13] = 5.7867900000e-09;
			__NEW_[14] = __OLD_[14] = 1.1981600000e-13;
			__NEW_[15] = __OLD_[15] = 4.9792300000e-19;
			__NEW_[16] = __OLD_[16] = 3.4584700000e-14;
			__NEW_[17] = __OLD_[17] = 1.8510600000e-14;
			__NEW_[18] = __OLD_[18] = 1.4237100000e+04;
			__NEW_[19] = __OLD_[19] = 2.0752000000e-02;
			__NEW_[20] = __OLD_[20] = 2.7913200000e-04;
			__NEW_[21] = __OLD_[21] = 7.1348300000e-07;
			__NEW_[22] = __OLD_[22] = 1.5317600000e-04;
			__NEW_[23] = __OLD_[23] = 6.7334500000e-07;
			__NEW_[24] = __OLD_[24] = 1.5578700000e-09;
			__NEW_[25] = __OLD_[25] = 1.1387900000e-02;
			__NEW_[26] = __OLD_[26] = 3.4278000000e-01;
			__NEW_[27] = __OLD_[27] = 1.4372000000e+05;
			__NEW_[28] = __OLD_[28] = 2.6556300000e-03;
			__NEW_[29] = __OLD_[29] = 9.9997700000e-01;
			__NEW_[30] = __OLD_[30] = 4.1706900000e-04;
			__NEW_[31] = __OLD_[31] = 9.9854300000e-01;
			__NEW_[32] = __OLD_[32] = 2.6275300000e-04;
			__NEW_[33] = __OLD_[33] = 4.1706900000e-04;
			__NEW_[34] = __OLD_[34] = 9.9854300000e-01;
			__NEW_[35] = __OLD_[35] = 4.1706900000e-04;
			__NEW_[36] = __OLD_[36] = 1.0000000000e+00;
			__NEW_[37] = __OLD_[37] = 6.4122900000e-04;
			__NEW_[38] = __OLD_[38] = 9.9251300000e-04;
			__NEW_[39] = __OLD_[39] = 1.7529800000e-04;
			__NEW_[40] = __OLD_[40] = 3.1912900000e-05;
			_prvt_time_new += _prvt_dtime;
			if(omp_get_thread_num()==_prvt_tree_thread[0])
			{
				_prvt_calc_i_stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
				_prvt_calc_Bi = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[1]),2.0000000000e+00))),(-1.0000000000e+00));
				_prvt_calc_Bss = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[2]),2.0000000000e+00))),(-1.0000000000e+00));
				_prvt_calc_BJSR = pow((1.0000000000e+00+((_prvt_CSQN_tot*_prvt_Km_CSQN)/pow((_prvt_Km_CSQN+__OLD_[3]),2.0000000000e+00))),(-1.0000000000e+00));
				_prvt_calc_J_rel = (_prvt_v1*(__OLD_[8]+__OLD_[9])*(__OLD_[3]-__OLD_[2])*__OLD_[5]);
				_prvt_calc_J_tr = ((__OLD_[4]-__OLD_[3])/_prvt_tau_tr);
				_prvt_calc_J_xfer = ((__OLD_[2]-__OLD_[1])/_prvt_tau_xfer);
				_prvt_calc_J_leak = (_prvt_v2*(__OLD_[4]-__OLD_[1]));
				_prvt_calc_J_up = ((_prvt_v3*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_up,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
				_prvt_calc_J_trpn = (((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))+(_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6])))-((_prvt_k_minus_htrpn*__OLD_[7])+(_prvt_k_minus_ltrpn*__OLD_[6])));
				_prvt_calc_i_CaL = (_prvt_g_CaL*__OLD_[11]*(__OLD_[0]-_prvt_E_CaL));
				_prvt_calc_i_pCa = ((_prvt_i_pCa_max*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_pCa,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
				_prvt_calc_i_NaCa = (((((((_prvt_k_NaCa*1.0000000000e+00)/(pow(_prvt_K_mNa,3.0000000000e+00)+pow(_prvt_Nao,3.0000000000e+00)))*1.0000000000e+00)/(_prvt_K_mCa+_prvt_Cao))*1.0000000000e+00)/(1.0000000000e+00+(_prvt_k_sat*exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))))*((exp(((_prvt_eta*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[18],3.0000000000e+00)*_prvt_Cao)-(exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Nao,3.0000000000e+00)*__OLD_[1])));
				_prvt_calc_E_CaN = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Cao/__OLD_[1])));
				_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((((9.0000000000e-01*_prvt_Nao)+(1.0000000000e-01*_prvt_Ko))/((9.0000000000e-01*__OLD_[18])+(1.0000000000e-01*__OLD_[27])))));
				_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ko/__OLD_[27])));
				_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[39]*(__OLD_[0]-(((_prvt_R*_prvt_T)/_prvt_F)*log((((9.8000000000e-01*_prvt_Ko)+(2.0000000000e-02*_prvt_Nao))/((9.8000000000e-01*__OLD_[27])+(2.0000000000e-02*__OLD_[18])))))));
				_prvt_calc_sigma = ((1.0000000000e+00/7.0000000000e+00)*(exp((_prvt_Nao/6.7300000000e+04))-1.0000000000e+00));
				_prvt_calc_O_ClCa = (2.0000000000e-01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.6700000000e+01))/7.8000000000e+00))));
				_prvt_calc_i_Nab = (_prvt_g_Nab*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_Kto_s = (_prvt_g_Kto_s*__OLD_[30]*__OLD_[31]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_K1 = ((((2.9380000000e-01*_prvt_Ko)/(_prvt_Ko+2.1000000000e+02))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp((8.9600000000e-02*(__OLD_[0]-_prvt_calc_E_K)))));
				_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[32],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_Kur = (_prvt_g_Kur*__OLD_[33]*__OLD_[34]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_Kss = (_prvt_g_Kss*__OLD_[35]*__OLD_[36]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_Cab = (_prvt_g_Cab*(__OLD_[0]-_prvt_calc_E_CaN));
				_prvt_calc_i_Na = (_prvt_g_Na*__OLD_[21]*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_Kto_f = (_prvt_g_Kto_f*pow(__OLD_[28],3.0000000000e+00)*__OLD_[29]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_f_NaK = (1.0000000000e+00/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.6500000000e-02*_prvt_calc_sigma*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_ClCa = (((_prvt_g_ClCa*_prvt_calc_O_ClCa*__OLD_[1])/(__OLD_[1]+_prvt_Km_Cl))*(__OLD_[0]-_prvt_E_Cl));
				_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_calc_f_NaK*1.0000000000e+00)/(1.0000000000e+00+pow((_prvt_Km_Nai/__OLD_[18]),1.5000000000e+00)))*_prvt_Ko)/(_prvt_Ko+_prvt_Km_Ko));
				__K1_[0]= (-(_prvt_calc_i_CaL+_prvt_calc_i_pCa+_prvt_calc_i_NaCa+_prvt_calc_i_Cab+_prvt_calc_i_Na+_prvt_calc_i_Nab+_prvt_calc_i_NaK+_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kur+_prvt_calc_i_Kss+_prvt_calc_i_Kr+_prvt_calc_i_ClCa+_prvt_calc_i_stim));
				__NEW_[0]= __K1_[0] * _prvt_dtime + __OLD_[0];
				__K1_[1]= (_prvt_calc_Bi*((_prvt_calc_J_leak+_prvt_calc_J_xfer)-(_prvt_calc_J_up+_prvt_calc_J_trpn+((((_prvt_calc_i_Cab+_prvt_calc_i_pCa)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vmyo*_prvt_F)))));
				__NEW_[1]= __K1_[1] * _prvt_dtime + __OLD_[1];
				__K1_[2]= (_prvt_calc_Bss*(((_prvt_calc_J_rel*_prvt_VJSR)/_prvt_Vss)-(((_prvt_calc_J_xfer*_prvt_Vmyo)/_prvt_Vss)+((_prvt_calc_i_CaL*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vss*_prvt_F)))));
				__NEW_[2]= __K1_[2] * _prvt_dtime + __OLD_[2];
				__K1_[3]= (_prvt_calc_BJSR*(_prvt_calc_J_tr-_prvt_calc_J_rel));
				__NEW_[3]= __K1_[3] * _prvt_dtime + __OLD_[3];
				__K1_[4]= ((((_prvt_calc_J_up-_prvt_calc_J_leak)*_prvt_Vmyo)/_prvt_VNSR)-((_prvt_calc_J_tr*_prvt_VJSR)/_prvt_VNSR));
				__NEW_[4]= __K1_[4] * _prvt_dtime + __OLD_[4];
				__K1_[5]= (((-4.0000000000e-02)*__OLD_[5])-(((1.0000000000e-01*_prvt_calc_i_CaL)/_prvt_i_CaL_max)*exp(((-pow((__OLD_[0]-5.0000000000e+00),2.0000000000e+00))/6.4800000000e+02))));
				__NEW_[5]= __K1_[5] * _prvt_dtime + __OLD_[5];
				__K1_[18]= (((-(_prvt_calc_i_Na+_prvt_calc_i_Nab+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
				__NEW_[18]= __K1_[18] * _prvt_dtime + __OLD_[18];
				__K1_[27]= (((-((_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kss+_prvt_calc_i_Kur+_prvt_calc_i_Kr)-(2.0000000000e+00*_prvt_calc_i_NaK)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
				__NEW_[27]= __K1_[27] * _prvt_dtime + __OLD_[27];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[1])
			{
				_prvt_calc_P_C1 = (1.0000000000e+00-(__OLD_[10]+__OLD_[8]+__OLD_[9]));
				__K1_[8]= (((_prvt_k_plus_a*pow(__OLD_[2],_prvt_n)*_prvt_calc_P_C1)+(_prvt_k_minus_b*__OLD_[9])+(_prvt_k_minus_c*__OLD_[10]))-((_prvt_k_minus_a*__OLD_[8])+(_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])+(_prvt_k_plus_c*__OLD_[8])));
				__NEW_[8]= __K1_[8] * _prvt_dtime + __OLD_[8];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[2])
			{
				_prvt_calc_C1 = (1.0000000000e+00-(__OLD_[11]+__OLD_[12]+__OLD_[12]+__OLD_[13]+__OLD_[14]+__OLD_[15]+__OLD_[16]+__OLD_[17]));
				_prvt_calc_alpha = ((4.0000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01))*((1.0000000000e+00+(7.0000000000e-01*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.0000000000e+01))))-(7.5000000000e-01*exp(((-pow((__OLD_[0]+2.0000000000e+01),2.0000000000e+00))/4.0000000000e+02)))))/(1.0000000000e+00+(1.2000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01)))));
				_prvt_calc_beta = (5.0000000000e-02*exp(((-(__OLD_[0]+1.2000000000e+01))/1.3000000000e+01)));
				_prvt_calc_gamma = ((_prvt_Kpc_max*__OLD_[2])/(_prvt_Kpc_half+__OLD_[2]));
				_prvt_calc_Kpcf = (1.3000000000e+01*(1.0000000000e+00-exp(((-pow((__OLD_[0]+1.4500000000e+01),2.0000000000e+00))/1.0000000000e+02))));
				__K1_[11]= (((_prvt_calc_alpha*__OLD_[14])+(_prvt_Kpcb*__OLD_[15])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[11]))))-((4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(_prvt_calc_gamma*__OLD_[11])));
				__NEW_[11]= __K1_[11] * _prvt_dtime + __OLD_[11];
				__K1_[12]= (((4.0000000000e+00*_prvt_calc_alpha*_prvt_calc_C1)+(2.0000000000e+00*_prvt_calc_beta*__OLD_[13]))-((_prvt_calc_beta*__OLD_[12])+(3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])));
				__NEW_[12]= __K1_[12] * _prvt_dtime + __OLD_[12];
				__K1_[13]= (((3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])+(3.0000000000e+00*_prvt_calc_beta*__OLD_[14]))-((2.0000000000e+00*_prvt_calc_beta*__OLD_[13])+(2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])));
				__NEW_[13]= __K1_[13] * _prvt_dtime + __OLD_[13];
				__K1_[14]= (((2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])+(4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(1.0000000000e-02*((4.0000000000e+00*_prvt_Kpcb*_prvt_calc_beta*__OLD_[15])-(_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])))+(2.0000000000e-03*((4.0000000000e+00*_prvt_calc_beta*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[14])))+(4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17]))-((3.0000000000e+00*_prvt_calc_beta*__OLD_[14])+(_prvt_calc_alpha*__OLD_[14])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14])));
				__NEW_[14]= __K1_[14] * _prvt_dtime + __OLD_[14];
				__K1_[15]= (((_prvt_calc_gamma*__OLD_[11])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[17])-(_prvt_calc_Kpcf*__OLD_[15])))+(1.0000000000e-02*((_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*_prvt_calc_Kpcf*__OLD_[15]))))-(_prvt_Kpcb*__OLD_[15]));
				__NEW_[15]= __K1_[15] * _prvt_dtime + __OLD_[15];
				__K1_[16]= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[11])-(_prvt_calc_alpha*__OLD_[16])))+(_prvt_Kpcb*__OLD_[17])+(2.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*__OLD_[16]))))-(_prvt_calc_gamma*__OLD_[16]));
				__NEW_[16]= __K1_[16] * _prvt_dtime + __OLD_[16];
				__K1_[17]= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[15])-(_prvt_calc_alpha*__OLD_[17])))+(_prvt_calc_gamma*__OLD_[16])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14]))-((4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17])+(_prvt_Kpcb*__OLD_[17])));
				__NEW_[17]= __K1_[17] * _prvt_dtime + __OLD_[17];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[3])
			{
				_prvt_calc_C_Na3 = (1.0000000000e+00-(__OLD_[21]+__OLD_[20]+__OLD_[19]+__OLD_[22]+__OLD_[23]+__OLD_[24]+__OLD_[25]+__OLD_[26]));
				_prvt_calc_alpha_Na11 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.7000000000e+01)))+(2.0000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
				_prvt_calc_alpha_Na12 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+01)))+(2.3000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
				_prvt_calc_alpha_Na13 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.2000000000e+01)))+(2.5000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
				_prvt_calc_beta_Na11 = (1.9170000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/2.0300000000e+01)));
				_prvt_calc_beta_Na12 = (2.0000000000e-01*exp(((-(__OLD_[0]-2.5000000000e+00))/2.0300000000e+01)));
				_prvt_calc_beta_Na13 = (2.2000000000e-01*exp(((-(__OLD_[0]-7.5000000000e+00))/2.0300000000e+01)));
				_prvt_calc_alpha_Na3 = (7.0000000000e-07*exp(((-(__OLD_[0]+7.0000000000e+00))/7.7000000000e+00)));
				_prvt_calc_beta_Na3 = (8.4000000000e-03+(2.0000000000e-05*(__OLD_[0]+7.0000000000e+00)));
				_prvt_calc_alpha_Na2 = (1.0000000000e+00/((1.8849500000e-01*exp(((-(__OLD_[0]+7.0000000000e+00))/1.6600000000e+01)))+3.9395600000e-01));
				_prvt_calc_beta_Na2 = ((_prvt_calc_alpha_Na13*_prvt_calc_alpha_Na2*_prvt_calc_alpha_Na3)/(_prvt_calc_beta_Na13*_prvt_calc_beta_Na3));
				_prvt_calc_alpha_Na4 = (_prvt_calc_alpha_Na2/1.0000000000e+03);
				_prvt_calc_beta_Na4 = _prvt_calc_alpha_Na3;
				_prvt_calc_alpha_Na5 = (_prvt_calc_alpha_Na2/9.5000000000e+04);
				_prvt_calc_beta_Na5 = (_prvt_calc_alpha_Na3/5.0000000000e+01);
				__K1_[19]= (((_prvt_calc_alpha_Na11*_prvt_calc_C_Na3)+(_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[19])+(_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na3*__OLD_[19])));
				__NEW_[19]= __K1_[19] * _prvt_dtime + __OLD_[19];
				__K1_[20]= (((_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na3*__OLD_[22]))-((_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na3*__OLD_[20])));
				__NEW_[20]= __K1_[20] * _prvt_dtime + __OLD_[20];
				__K1_[21]= (((_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na2*__OLD_[22]))-((_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na2*__OLD_[21])));
				__NEW_[21]= __K1_[21] * _prvt_dtime + __OLD_[21];
				__K1_[22]= (((_prvt_calc_alpha_Na2*__OLD_[21])+(_prvt_calc_beta_Na3*__OLD_[20])+(_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na12*__OLD_[25]))-((_prvt_calc_beta_Na2*__OLD_[22])+(_prvt_calc_alpha_Na3*__OLD_[22])+(_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na12*__OLD_[22])));
				__NEW_[22]= __K1_[22] * _prvt_dtime + __OLD_[22];
				__K1_[23]= (((_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na5*__OLD_[24]))-((_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na5*__OLD_[23])));
				__NEW_[23]= __K1_[23] * _prvt_dtime + __OLD_[23];
				__K1_[24]= ((_prvt_calc_alpha_Na5*__OLD_[23])-(_prvt_calc_beta_Na5*__OLD_[24]));
				__NEW_[24]= __K1_[24] * _prvt_dtime + __OLD_[24];
				__K1_[25]= (((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_beta_Na12*__OLD_[22])+(_prvt_calc_beta_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_alpha_Na12*__OLD_[25])+(_prvt_calc_alpha_Na3*__OLD_[25])));
				__NEW_[25]= __K1_[25] * _prvt_dtime + __OLD_[25];
				__K1_[26]= (((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_beta_Na3*_prvt_calc_C_Na3))-((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_alpha_Na3*__OLD_[26])));
				__NEW_[26]= __K1_[26] * _prvt_dtime + __OLD_[26];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[4])
			{
				_prvt_calc_alpha_a = (1.8064000000e-01*exp((3.5770000000e-02*(__OLD_[0]+3.0000000000e+01))));
				_prvt_calc_beta_a = (3.9560000000e-01*exp(((-6.2370000000e-02)*(__OLD_[0]+3.0000000000e+01))));
				__K1_[28]= ((_prvt_calc_alpha_a*(1.0000000000e+00-__OLD_[28]))-(_prvt_calc_beta_a*__OLD_[28]));
				__NEW_[28]= __K1_[28] * _prvt_dtime + __OLD_[28];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[5])
			{
				_prvt_calc_alpha_i = ((1.5200000000e-04*exp(((-(__OLD_[0]+1.3500000000e+01))/7.0000000000e+00)))/((6.7083000000e-03*exp(((-(__OLD_[0]+3.3500000000e+01))/7.0000000000e+00)))+1.0000000000e+00));
				_prvt_calc_beta_i = ((9.5000000000e-04*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))/((5.1335000000e-02*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))+1.0000000000e+00));
				__K1_[29]= ((_prvt_calc_alpha_i*(1.0000000000e+00-__OLD_[29]))-(_prvt_calc_beta_i*__OLD_[29]));
				__NEW_[29]= __K1_[29] * _prvt_dtime + __OLD_[29];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[6])
			{
				_prvt_calc_ass = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2500000000e+01))/7.7000000000e+00))));
				_prvt_calc_tau_ta_s = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
				_prvt_calc_tau_aur = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
				_prvt_calc_tau_Kss = ((3.9300000000e+01*exp(((-8.6200000000e-02)*__OLD_[0])))+1.3170000000e+01);
				__K1_[30]= ((_prvt_calc_ass-__OLD_[30])/_prvt_calc_tau_ta_s);
				__NEW_[30]= __K1_[30] * _prvt_dtime + __OLD_[30];
				__K1_[33]= ((_prvt_calc_ass-__OLD_[33])/_prvt_calc_tau_aur);
				__NEW_[33]= __K1_[33] * _prvt_dtime + __OLD_[33];
				__K1_[35]= ((_prvt_calc_ass-__OLD_[35])/_prvt_calc_tau_Kss);
				__NEW_[35]= __K1_[35] * _prvt_dtime + __OLD_[35];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[7])
			{
				_prvt_calc_iss = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00))));
				_prvt_calc_tau_ti_s = (2.7000000000e+02+(1.0500000000e+03/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
				_prvt_calc_tau_iur = (1.2000000000e+03-(1.7000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
				__K1_[31]= ((_prvt_calc_iss-__OLD_[31])/_prvt_calc_tau_ti_s);
				__NEW_[31]= __K1_[31] * _prvt_dtime + __OLD_[31];
				__K1_[34]= ((_prvt_calc_iss-__OLD_[34])/_prvt_calc_tau_iur);
				__NEW_[34]= __K1_[34] * _prvt_dtime + __OLD_[34];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[8])
			{
				_prvt_calc_alpha_n = ((4.8133300000e-06*(__OLD_[0]+2.6500000000e+01))/(1.0000000000e+00-exp(((-1.2800000000e-01)*(__OLD_[0]+2.6500000000e+01)))));
				_prvt_calc_beta_n = (9.5333300000e-05*exp(((-3.8000000000e-02)*(__OLD_[0]+2.6500000000e+01))));
				__K1_[32]= ((_prvt_calc_alpha_n*(1.0000000000e+00-__OLD_[32]))-(_prvt_calc_beta_n*__OLD_[32]));
				__NEW_[32]= __K1_[32] * _prvt_dtime + __OLD_[32];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[9])
			{
				_prvt_calc_C_K0 = (1.0000000000e+00-(__OLD_[38]+__OLD_[37]+__OLD_[39]+__OLD_[40]));
				_prvt_calc_alpha_a0 = (2.2348000000e-02*exp((1.1760000000e-02*__OLD_[0])));
				_prvt_calc_beta_a0 = (4.7002000000e-02*exp(((-6.3100000000e-02)*__OLD_[0])));
				__K1_[38]= (((_prvt_calc_alpha_a0*_prvt_calc_C_K0)+(_prvt_kb*__OLD_[37]))-((_prvt_calc_beta_a0*__OLD_[38])+(_prvt_kf*__OLD_[38])));
				__NEW_[38]= __K1_[38] * _prvt_dtime + __OLD_[38];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[10])
			{
				_prvt_calc_alpha_a1 = (1.3733000000e-02*exp((3.8198000000e-02*__OLD_[0])));
				_prvt_calc_beta_a1 = (6.8900000000e-05*exp(((-4.1780000000e-02)*__OLD_[0])));
				_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (9.0821000000e-02*exp((2.3391000000e-02*(__OLD_[0]+5.0000000000e+00))));
				_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (6.4970000000e-03*exp(((-3.2680000000e-02)*(__OLD_[0]+5.0000000000e+00))));
				__K1_[37]= (((_prvt_kf*__OLD_[38])+(_prvt_calc_beta_a1*__OLD_[39]))-((_prvt_kb*__OLD_[37])+(_prvt_calc_alpha_a1*__OLD_[37])));
				__NEW_[37]= __K1_[37] * _prvt_dtime + __OLD_[37];
				__K1_[39]= (((_prvt_calc_alpha_a1*__OLD_[37])+(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]))-((_prvt_calc_beta_a1*__OLD_[39])+(_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])));
				__NEW_[39]= __K1_[39] * _prvt_dtime + __OLD_[39];
				__K1_[40]= ((_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])-(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]));
				__NEW_[40]= __K1_[40] * _prvt_dtime + __OLD_[40];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[11])
			{
				__K1_[6]= ((_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6]))-(_prvt_k_minus_ltrpn*__OLD_[6]));
				__NEW_[6]= __K1_[6] * _prvt_dtime + __OLD_[6];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[12])
			{
				__K1_[7]= ((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))-(_prvt_k_minus_htrpn*__OLD_[7]));
				__NEW_[7]= __K1_[7] * _prvt_dtime + __OLD_[7];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[13])
			{
				__K1_[9]= ((_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])-(_prvt_k_minus_b*__OLD_[9]));
				__NEW_[9]= __K1_[9] * _prvt_dtime + __OLD_[9];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[14])
			{
				__K1_[10]= ((_prvt_k_plus_c*__OLD_[8])-(_prvt_k_minus_c*__OLD_[10]));
				__NEW_[10]= __K1_[10] * _prvt_dtime + __OLD_[10];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[15])
			{
				__K1_[36]= 0.0000000000e+00;
				__NEW_[36]= __K1_[36] * _prvt_dtime + __OLD_[36];
			}
			//store the old iteration in a aux 
			__TEMP_ = __OLD_;
			__OLD_ = __OLD_AUX_;
			__OLD_AUX_ = __TEMP_;
			//steps ahead with euler
			__TEMP_ = __NEW_;
			__NEW_ = __OLD_;
			__OLD_ = __TEMP_;
			//as threads devem começar o  laço ao mesmo tempo
// 			timeval t1, t2;
// 			gettimeofday(&t1, NULL);
			
			#pragma omp barrier
// 			gettimeofday(&t2, NULL);
		/*	#pragma omp critical
			{
			    int aux = (int)((int)t2.tv_usec - (int)t1.tv_usec);
			    tempoespera += (aux<0)?0:aux;
			}
		*/	
			while(_prvt_time_new<=_prvt_finalTime){
				_prvt_time_new += _prvt_dtime;
				if(omp_get_thread_num()==_prvt_tree_thread[0])
				{
					_prvt_calc_i_stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
					_prvt_calc_Bi = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[1]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_Bss = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[2]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_BJSR = pow((1.0000000000e+00+((_prvt_CSQN_tot*_prvt_Km_CSQN)/pow((_prvt_Km_CSQN+__OLD_[3]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_J_rel = (_prvt_v1*(__OLD_[8]+__OLD_[9])*(__OLD_[3]-__OLD_[2])*__OLD_[5]);
					_prvt_calc_J_tr = ((__OLD_[4]-__OLD_[3])/_prvt_tau_tr);
					_prvt_calc_J_xfer = ((__OLD_[2]-__OLD_[1])/_prvt_tau_xfer);
					_prvt_calc_J_leak = (_prvt_v2*(__OLD_[4]-__OLD_[1]));
					_prvt_calc_J_up = ((_prvt_v3*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_up,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
					_prvt_calc_J_trpn = (((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))+(_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6])))-((_prvt_k_minus_htrpn*__OLD_[7])+(_prvt_k_minus_ltrpn*__OLD_[6])));
					_prvt_calc_i_CaL = (_prvt_g_CaL*__OLD_[11]*(__OLD_[0]-_prvt_E_CaL));
					_prvt_calc_i_pCa = ((_prvt_i_pCa_max*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_pCa,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
					_prvt_calc_i_NaCa = (((((((_prvt_k_NaCa*1.0000000000e+00)/(pow(_prvt_K_mNa,3.0000000000e+00)+pow(_prvt_Nao,3.0000000000e+00)))*1.0000000000e+00)/(_prvt_K_mCa+_prvt_Cao))*1.0000000000e+00)/(1.0000000000e+00+(_prvt_k_sat*exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))))*((exp(((_prvt_eta*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[18],3.0000000000e+00)*_prvt_Cao)-(exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Nao,3.0000000000e+00)*__OLD_[1])));
					_prvt_calc_E_CaN = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Cao/__OLD_[1])));
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((((9.0000000000e-01*_prvt_Nao)+(1.0000000000e-01*_prvt_Ko))/((9.0000000000e-01*__OLD_[18])+(1.0000000000e-01*__OLD_[27])))));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ko/__OLD_[27])));
					_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[39]*(__OLD_[0]-(((_prvt_R*_prvt_T)/_prvt_F)*log((((9.8000000000e-01*_prvt_Ko)+(2.0000000000e-02*_prvt_Nao))/((9.8000000000e-01*__OLD_[27])+(2.0000000000e-02*__OLD_[18])))))));
					_prvt_calc_sigma = ((1.0000000000e+00/7.0000000000e+00)*(exp((_prvt_Nao/6.7300000000e+04))-1.0000000000e+00));
					_prvt_calc_O_ClCa = (2.0000000000e-01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.6700000000e+01))/7.8000000000e+00))));
					_prvt_calc_i_Nab = (_prvt_g_Nab*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Kto_s = (_prvt_g_Kto_s*__OLD_[30]*__OLD_[31]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K1 = ((((2.9380000000e-01*_prvt_Ko)/(_prvt_Ko+2.1000000000e+02))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp((8.9600000000e-02*(__OLD_[0]-_prvt_calc_E_K)))));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[32],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Kur = (_prvt_g_Kur*__OLD_[33]*__OLD_[34]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Kss = (_prvt_g_Kss*__OLD_[35]*__OLD_[36]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Cab = (_prvt_g_Cab*(__OLD_[0]-_prvt_calc_E_CaN));
					_prvt_calc_i_Na = (_prvt_g_Na*__OLD_[21]*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Kto_f = (_prvt_g_Kto_f*pow(__OLD_[28],3.0000000000e+00)*__OLD_[29]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_f_NaK = (1.0000000000e+00/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.6500000000e-02*_prvt_calc_sigma*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_ClCa = (((_prvt_g_ClCa*_prvt_calc_O_ClCa*__OLD_[1])/(__OLD_[1]+_prvt_Km_Cl))*(__OLD_[0]-_prvt_E_Cl));
					_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_calc_f_NaK*1.0000000000e+00)/(1.0000000000e+00+pow((_prvt_Km_Nai/__OLD_[18]),1.5000000000e+00)))*_prvt_Ko)/(_prvt_Ko+_prvt_Km_Ko));
					__K2_[0]= (-(_prvt_calc_i_CaL+_prvt_calc_i_pCa+_prvt_calc_i_NaCa+_prvt_calc_i_Cab+_prvt_calc_i_Na+_prvt_calc_i_Nab+_prvt_calc_i_NaK+_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kur+_prvt_calc_i_Kss+_prvt_calc_i_Kr+_prvt_calc_i_ClCa+_prvt_calc_i_stim));
					_prvt_aux_tol = fabs(__OLD_[0])*_prvt_rel_tol_;
					__TOL_[0] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[0] = fabs((_prvt_dtime/2) * (__K1_[0] - __K2_[0])/__TOL_[0]);
					__K2_[1]= (_prvt_calc_Bi*((_prvt_calc_J_leak+_prvt_calc_J_xfer)-(_prvt_calc_J_up+_prvt_calc_J_trpn+((((_prvt_calc_i_Cab+_prvt_calc_i_pCa)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vmyo*_prvt_F)))));
					_prvt_aux_tol = fabs(__OLD_[1])*_prvt_rel_tol_;
					__TOL_[1] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[1] = fabs((_prvt_dtime/2) * (__K1_[1] - __K2_[1])/__TOL_[1]);
					__K2_[2]= (_prvt_calc_Bss*(((_prvt_calc_J_rel*_prvt_VJSR)/_prvt_Vss)-(((_prvt_calc_J_xfer*_prvt_Vmyo)/_prvt_Vss)+((_prvt_calc_i_CaL*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vss*_prvt_F)))));
					_prvt_aux_tol = fabs(__OLD_[2])*_prvt_rel_tol_;
					__TOL_[2] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[2] = fabs((_prvt_dtime/2) * (__K1_[2] - __K2_[2])/__TOL_[2]);
					__K2_[3]= (_prvt_calc_BJSR*(_prvt_calc_J_tr-_prvt_calc_J_rel));
					_prvt_aux_tol = fabs(__OLD_[3])*_prvt_rel_tol_;
					__TOL_[3] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[3] = fabs((_prvt_dtime/2) * (__K1_[3] - __K2_[3])/__TOL_[3]);
					__K2_[4]= ((((_prvt_calc_J_up-_prvt_calc_J_leak)*_prvt_Vmyo)/_prvt_VNSR)-((_prvt_calc_J_tr*_prvt_VJSR)/_prvt_VNSR));
					_prvt_aux_tol = fabs(__OLD_[4])*_prvt_rel_tol_;
					__TOL_[4] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[4] = fabs((_prvt_dtime/2) * (__K1_[4] - __K2_[4])/__TOL_[4]);
					__K2_[5]= (((-4.0000000000e-02)*__OLD_[5])-(((1.0000000000e-01*_prvt_calc_i_CaL)/_prvt_i_CaL_max)*exp(((-pow((__OLD_[0]-5.0000000000e+00),2.0000000000e+00))/6.4800000000e+02))));
					_prvt_aux_tol = fabs(__OLD_[5])*_prvt_rel_tol_;
					__TOL_[5] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[5] = fabs((_prvt_dtime/2) * (__K1_[5] - __K2_[5])/__TOL_[5]);
					__K2_[18]= (((-(_prvt_calc_i_Na+_prvt_calc_i_Nab+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
					_prvt_aux_tol = fabs(__OLD_[18])*_prvt_rel_tol_;
					__TOL_[18] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[18] = fabs((_prvt_dtime/2) * (__K1_[18] - __K2_[18])/__TOL_[18]);
					__K2_[27]= (((-((_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kss+_prvt_calc_i_Kur+_prvt_calc_i_Kr)-(2.0000000000e+00*_prvt_calc_i_NaK)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
					_prvt_aux_tol = fabs(__OLD_[27])*_prvt_rel_tol_;
					__TOL_[27] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[27] = fabs((_prvt_dtime/2) * (__K1_[27] - __K2_[27])/__TOL_[27]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[1])
				{
					_prvt_calc_P_C1 = (1.0000000000e+00-(__OLD_[10]+__OLD_[8]+__OLD_[9]));
					__K2_[8]= (((_prvt_k_plus_a*pow(__OLD_[2],_prvt_n)*_prvt_calc_P_C1)+(_prvt_k_minus_b*__OLD_[9])+(_prvt_k_minus_c*__OLD_[10]))-((_prvt_k_minus_a*__OLD_[8])+(_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])+(_prvt_k_plus_c*__OLD_[8])));
					_prvt_aux_tol = fabs(__OLD_[8])*_prvt_rel_tol_;
					__TOL_[8] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[8] = fabs((_prvt_dtime/2) * (__K1_[8] - __K2_[8])/__TOL_[8]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[2])
				{
					_prvt_calc_C1 = (1.0000000000e+00-(__OLD_[11]+__OLD_[12]+__OLD_[12]+__OLD_[13]+__OLD_[14]+__OLD_[15]+__OLD_[16]+__OLD_[17]));
					_prvt_calc_alpha = ((4.0000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01))*((1.0000000000e+00+(7.0000000000e-01*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.0000000000e+01))))-(7.5000000000e-01*exp(((-pow((__OLD_[0]+2.0000000000e+01),2.0000000000e+00))/4.0000000000e+02)))))/(1.0000000000e+00+(1.2000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01)))));
					_prvt_calc_beta = (5.0000000000e-02*exp(((-(__OLD_[0]+1.2000000000e+01))/1.3000000000e+01)));
					_prvt_calc_gamma = ((_prvt_Kpc_max*__OLD_[2])/(_prvt_Kpc_half+__OLD_[2]));
					_prvt_calc_Kpcf = (1.3000000000e+01*(1.0000000000e+00-exp(((-pow((__OLD_[0]+1.4500000000e+01),2.0000000000e+00))/1.0000000000e+02))));
					__K2_[11]= (((_prvt_calc_alpha*__OLD_[14])+(_prvt_Kpcb*__OLD_[15])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[11]))))-((4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(_prvt_calc_gamma*__OLD_[11])));
					_prvt_aux_tol = fabs(__OLD_[11])*_prvt_rel_tol_;
					__TOL_[11] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[11] = fabs((_prvt_dtime/2) * (__K1_[11] - __K2_[11])/__TOL_[11]);
					__K2_[12]= (((4.0000000000e+00*_prvt_calc_alpha*_prvt_calc_C1)+(2.0000000000e+00*_prvt_calc_beta*__OLD_[13]))-((_prvt_calc_beta*__OLD_[12])+(3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])));
					_prvt_aux_tol = fabs(__OLD_[12])*_prvt_rel_tol_;
					__TOL_[12] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[12] = fabs((_prvt_dtime/2) * (__K1_[12] - __K2_[12])/__TOL_[12]);
					__K2_[13]= (((3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])+(3.0000000000e+00*_prvt_calc_beta*__OLD_[14]))-((2.0000000000e+00*_prvt_calc_beta*__OLD_[13])+(2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])));
					_prvt_aux_tol = fabs(__OLD_[13])*_prvt_rel_tol_;
					__TOL_[13] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[13] = fabs((_prvt_dtime/2) * (__K1_[13] - __K2_[13])/__TOL_[13]);
					__K2_[14]= (((2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])+(4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(1.0000000000e-02*((4.0000000000e+00*_prvt_Kpcb*_prvt_calc_beta*__OLD_[15])-(_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])))+(2.0000000000e-03*((4.0000000000e+00*_prvt_calc_beta*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[14])))+(4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17]))-((3.0000000000e+00*_prvt_calc_beta*__OLD_[14])+(_prvt_calc_alpha*__OLD_[14])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14])));
					_prvt_aux_tol = fabs(__OLD_[14])*_prvt_rel_tol_;
					__TOL_[14] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[14] = fabs((_prvt_dtime/2) * (__K1_[14] - __K2_[14])/__TOL_[14]);
					__K2_[15]= (((_prvt_calc_gamma*__OLD_[11])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[17])-(_prvt_calc_Kpcf*__OLD_[15])))+(1.0000000000e-02*((_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*_prvt_calc_Kpcf*__OLD_[15]))))-(_prvt_Kpcb*__OLD_[15]));
					_prvt_aux_tol = fabs(__OLD_[15])*_prvt_rel_tol_;
					__TOL_[15] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[15] = fabs((_prvt_dtime/2) * (__K1_[15] - __K2_[15])/__TOL_[15]);
					__K2_[16]= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[11])-(_prvt_calc_alpha*__OLD_[16])))+(_prvt_Kpcb*__OLD_[17])+(2.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*__OLD_[16]))))-(_prvt_calc_gamma*__OLD_[16]));
					_prvt_aux_tol = fabs(__OLD_[16])*_prvt_rel_tol_;
					__TOL_[16] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[16] = fabs((_prvt_dtime/2) * (__K1_[16] - __K2_[16])/__TOL_[16]);
					__K2_[17]= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[15])-(_prvt_calc_alpha*__OLD_[17])))+(_prvt_calc_gamma*__OLD_[16])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14]))-((4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17])+(_prvt_Kpcb*__OLD_[17])));
					_prvt_aux_tol = fabs(__OLD_[17])*_prvt_rel_tol_;
					__TOL_[17] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[17] = fabs((_prvt_dtime/2) * (__K1_[17] - __K2_[17])/__TOL_[17]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[3])
				{
					_prvt_calc_C_Na3 = (1.0000000000e+00-(__OLD_[21]+__OLD_[20]+__OLD_[19]+__OLD_[22]+__OLD_[23]+__OLD_[24]+__OLD_[25]+__OLD_[26]));
					_prvt_calc_alpha_Na11 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.7000000000e+01)))+(2.0000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_alpha_Na12 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+01)))+(2.3000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_alpha_Na13 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.2000000000e+01)))+(2.5000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_beta_Na11 = (1.9170000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_beta_Na12 = (2.0000000000e-01*exp(((-(__OLD_[0]-2.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_beta_Na13 = (2.2000000000e-01*exp(((-(__OLD_[0]-7.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_alpha_Na3 = (7.0000000000e-07*exp(((-(__OLD_[0]+7.0000000000e+00))/7.7000000000e+00)));
					_prvt_calc_beta_Na3 = (8.4000000000e-03+(2.0000000000e-05*(__OLD_[0]+7.0000000000e+00)));
					_prvt_calc_alpha_Na2 = (1.0000000000e+00/((1.8849500000e-01*exp(((-(__OLD_[0]+7.0000000000e+00))/1.6600000000e+01)))+3.9395600000e-01));
					_prvt_calc_beta_Na2 = ((_prvt_calc_alpha_Na13*_prvt_calc_alpha_Na2*_prvt_calc_alpha_Na3)/(_prvt_calc_beta_Na13*_prvt_calc_beta_Na3));
					_prvt_calc_alpha_Na4 = (_prvt_calc_alpha_Na2/1.0000000000e+03);
					_prvt_calc_beta_Na4 = _prvt_calc_alpha_Na3;
					_prvt_calc_alpha_Na5 = (_prvt_calc_alpha_Na2/9.5000000000e+04);
					_prvt_calc_beta_Na5 = (_prvt_calc_alpha_Na3/5.0000000000e+01);
					__K2_[19]= (((_prvt_calc_alpha_Na11*_prvt_calc_C_Na3)+(_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[19])+(_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na3*__OLD_[19])));
					_prvt_aux_tol = fabs(__OLD_[19])*_prvt_rel_tol_;
					__TOL_[19] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[19] = fabs((_prvt_dtime/2) * (__K1_[19] - __K2_[19])/__TOL_[19]);
					__K2_[20]= (((_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na3*__OLD_[22]))-((_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na3*__OLD_[20])));
					_prvt_aux_tol = fabs(__OLD_[20])*_prvt_rel_tol_;
					__TOL_[20] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[20] = fabs((_prvt_dtime/2) * (__K1_[20] - __K2_[20])/__TOL_[20]);
					__K2_[21]= (((_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na2*__OLD_[22]))-((_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na2*__OLD_[21])));
					_prvt_aux_tol = fabs(__OLD_[21])*_prvt_rel_tol_;
					__TOL_[21] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[21] = fabs((_prvt_dtime/2) * (__K1_[21] - __K2_[21])/__TOL_[21]);
					__K2_[22]= (((_prvt_calc_alpha_Na2*__OLD_[21])+(_prvt_calc_beta_Na3*__OLD_[20])+(_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na12*__OLD_[25]))-((_prvt_calc_beta_Na2*__OLD_[22])+(_prvt_calc_alpha_Na3*__OLD_[22])+(_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na12*__OLD_[22])));
					_prvt_aux_tol = fabs(__OLD_[22])*_prvt_rel_tol_;
					__TOL_[22] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[22] = fabs((_prvt_dtime/2) * (__K1_[22] - __K2_[22])/__TOL_[22]);
					__K2_[23]= (((_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na5*__OLD_[24]))-((_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na5*__OLD_[23])));
					_prvt_aux_tol = fabs(__OLD_[23])*_prvt_rel_tol_;
					__TOL_[23] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[23] = fabs((_prvt_dtime/2) * (__K1_[23] - __K2_[23])/__TOL_[23]);
					__K2_[24]= ((_prvt_calc_alpha_Na5*__OLD_[23])-(_prvt_calc_beta_Na5*__OLD_[24]));
					_prvt_aux_tol = fabs(__OLD_[24])*_prvt_rel_tol_;
					__TOL_[24] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[24] = fabs((_prvt_dtime/2) * (__K1_[24] - __K2_[24])/__TOL_[24]);
					__K2_[25]= (((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_beta_Na12*__OLD_[22])+(_prvt_calc_beta_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_alpha_Na12*__OLD_[25])+(_prvt_calc_alpha_Na3*__OLD_[25])));
					_prvt_aux_tol = fabs(__OLD_[25])*_prvt_rel_tol_;
					__TOL_[25] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[25] = fabs((_prvt_dtime/2) * (__K1_[25] - __K2_[25])/__TOL_[25]);
					__K2_[26]= (((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_beta_Na3*_prvt_calc_C_Na3))-((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_alpha_Na3*__OLD_[26])));
					_prvt_aux_tol = fabs(__OLD_[26])*_prvt_rel_tol_;
					__TOL_[26] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[26] = fabs((_prvt_dtime/2) * (__K1_[26] - __K2_[26])/__TOL_[26]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[4])
				{
					_prvt_calc_alpha_a = (1.8064000000e-01*exp((3.5770000000e-02*(__OLD_[0]+3.0000000000e+01))));
					_prvt_calc_beta_a = (3.9560000000e-01*exp(((-6.2370000000e-02)*(__OLD_[0]+3.0000000000e+01))));
					__K2_[28]= ((_prvt_calc_alpha_a*(1.0000000000e+00-__OLD_[28]))-(_prvt_calc_beta_a*__OLD_[28]));
					_prvt_aux_tol = fabs(__OLD_[28])*_prvt_rel_tol_;
					__TOL_[28] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[28] = fabs((_prvt_dtime/2) * (__K1_[28] - __K2_[28])/__TOL_[28]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[5])
				{
					_prvt_calc_alpha_i = ((1.5200000000e-04*exp(((-(__OLD_[0]+1.3500000000e+01))/7.0000000000e+00)))/((6.7083000000e-03*exp(((-(__OLD_[0]+3.3500000000e+01))/7.0000000000e+00)))+1.0000000000e+00));
					_prvt_calc_beta_i = ((9.5000000000e-04*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))/((5.1335000000e-02*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))+1.0000000000e+00));
					__K2_[29]= ((_prvt_calc_alpha_i*(1.0000000000e+00-__OLD_[29]))-(_prvt_calc_beta_i*__OLD_[29]));
					_prvt_aux_tol = fabs(__OLD_[29])*_prvt_rel_tol_;
					__TOL_[29] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[29] = fabs((_prvt_dtime/2) * (__K1_[29] - __K2_[29])/__TOL_[29]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[6])
				{
					_prvt_calc_ass = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2500000000e+01))/7.7000000000e+00))));
					_prvt_calc_tau_ta_s = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
					_prvt_calc_tau_aur = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
					_prvt_calc_tau_Kss = ((3.9300000000e+01*exp(((-8.6200000000e-02)*__OLD_[0])))+1.3170000000e+01);
					__K2_[30]= ((_prvt_calc_ass-__OLD_[30])/_prvt_calc_tau_ta_s);
					_prvt_aux_tol = fabs(__OLD_[30])*_prvt_rel_tol_;
					__TOL_[30] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[30] = fabs((_prvt_dtime/2) * (__K1_[30] - __K2_[30])/__TOL_[30]);
					__K2_[33]= ((_prvt_calc_ass-__OLD_[33])/_prvt_calc_tau_aur);
					_prvt_aux_tol = fabs(__OLD_[33])*_prvt_rel_tol_;
					__TOL_[33] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[33] = fabs((_prvt_dtime/2) * (__K1_[33] - __K2_[33])/__TOL_[33]);
					__K2_[35]= ((_prvt_calc_ass-__OLD_[35])/_prvt_calc_tau_Kss);
					_prvt_aux_tol = fabs(__OLD_[35])*_prvt_rel_tol_;
					__TOL_[35] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[35] = fabs((_prvt_dtime/2) * (__K1_[35] - __K2_[35])/__TOL_[35]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[7])
				{
					_prvt_calc_iss = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00))));
					_prvt_calc_tau_ti_s = (2.7000000000e+02+(1.0500000000e+03/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
					_prvt_calc_tau_iur = (1.2000000000e+03-(1.7000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
					__K2_[31]= ((_prvt_calc_iss-__OLD_[31])/_prvt_calc_tau_ti_s);
					_prvt_aux_tol = fabs(__OLD_[31])*_prvt_rel_tol_;
					__TOL_[31] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[31] = fabs((_prvt_dtime/2) * (__K1_[31] - __K2_[31])/__TOL_[31]);
					__K2_[34]= ((_prvt_calc_iss-__OLD_[34])/_prvt_calc_tau_iur);
					_prvt_aux_tol = fabs(__OLD_[34])*_prvt_rel_tol_;
					__TOL_[34] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[34] = fabs((_prvt_dtime/2) * (__K1_[34] - __K2_[34])/__TOL_[34]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[8])
				{
					_prvt_calc_alpha_n = ((4.8133300000e-06*(__OLD_[0]+2.6500000000e+01))/(1.0000000000e+00-exp(((-1.2800000000e-01)*(__OLD_[0]+2.6500000000e+01)))));
					_prvt_calc_beta_n = (9.5333300000e-05*exp(((-3.8000000000e-02)*(__OLD_[0]+2.6500000000e+01))));
					__K2_[32]= ((_prvt_calc_alpha_n*(1.0000000000e+00-__OLD_[32]))-(_prvt_calc_beta_n*__OLD_[32]));
					_prvt_aux_tol = fabs(__OLD_[32])*_prvt_rel_tol_;
					__TOL_[32] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[32] = fabs((_prvt_dtime/2) * (__K1_[32] - __K2_[32])/__TOL_[32]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[9])
				{
					_prvt_calc_C_K0 = (1.0000000000e+00-(__OLD_[38]+__OLD_[37]+__OLD_[39]+__OLD_[40]));
					_prvt_calc_alpha_a0 = (2.2348000000e-02*exp((1.1760000000e-02*__OLD_[0])));
					_prvt_calc_beta_a0 = (4.7002000000e-02*exp(((-6.3100000000e-02)*__OLD_[0])));
					__K2_[38]= (((_prvt_calc_alpha_a0*_prvt_calc_C_K0)+(_prvt_kb*__OLD_[37]))-((_prvt_calc_beta_a0*__OLD_[38])+(_prvt_kf*__OLD_[38])));
					_prvt_aux_tol = fabs(__OLD_[38])*_prvt_rel_tol_;
					__TOL_[38] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[38] = fabs((_prvt_dtime/2) * (__K1_[38] - __K2_[38])/__TOL_[38]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[10])
				{
					_prvt_calc_alpha_a1 = (1.3733000000e-02*exp((3.8198000000e-02*__OLD_[0])));
					_prvt_calc_beta_a1 = (6.8900000000e-05*exp(((-4.1780000000e-02)*__OLD_[0])));
					_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (9.0821000000e-02*exp((2.3391000000e-02*(__OLD_[0]+5.0000000000e+00))));
					_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (6.4970000000e-03*exp(((-3.2680000000e-02)*(__OLD_[0]+5.0000000000e+00))));
					__K2_[37]= (((_prvt_kf*__OLD_[38])+(_prvt_calc_beta_a1*__OLD_[39]))-((_prvt_kb*__OLD_[37])+(_prvt_calc_alpha_a1*__OLD_[37])));
					_prvt_aux_tol = fabs(__OLD_[37])*_prvt_rel_tol_;
					__TOL_[37] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[37] = fabs((_prvt_dtime/2) * (__K1_[37] - __K2_[37])/__TOL_[37]);
					__K2_[39]= (((_prvt_calc_alpha_a1*__OLD_[37])+(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]))-((_prvt_calc_beta_a1*__OLD_[39])+(_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])));
					_prvt_aux_tol = fabs(__OLD_[39])*_prvt_rel_tol_;
					__TOL_[39] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[39] = fabs((_prvt_dtime/2) * (__K1_[39] - __K2_[39])/__TOL_[39]);
					__K2_[40]= ((_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])-(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]));
					_prvt_aux_tol = fabs(__OLD_[40])*_prvt_rel_tol_;
					__TOL_[40] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[40] = fabs((_prvt_dtime/2) * (__K1_[40] - __K2_[40])/__TOL_[40]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[11])
				{
					__K2_[6]= ((_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6]))-(_prvt_k_minus_ltrpn*__OLD_[6]));
					_prvt_aux_tol = fabs(__OLD_[6])*_prvt_rel_tol_;
					__TOL_[6] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[6] = fabs((_prvt_dtime/2) * (__K1_[6] - __K2_[6])/__TOL_[6]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[12])
				{
					__K2_[7]= ((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))-(_prvt_k_minus_htrpn*__OLD_[7]));
					_prvt_aux_tol = fabs(__OLD_[7])*_prvt_rel_tol_;
					__TOL_[7] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[7] = fabs((_prvt_dtime/2) * (__K1_[7] - __K2_[7])/__TOL_[7]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[13])
				{
					__K2_[9]= ((_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])-(_prvt_k_minus_b*__OLD_[9]));
					_prvt_aux_tol = fabs(__OLD_[9])*_prvt_rel_tol_;
					__TOL_[9] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[9] = fabs((_prvt_dtime/2) * (__K1_[9] - __K2_[9])/__TOL_[9]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[14])
				{
					__K2_[10]= ((_prvt_k_plus_c*__OLD_[8])-(_prvt_k_minus_c*__OLD_[10]));
					_prvt_aux_tol = fabs(__OLD_[10])*_prvt_rel_tol_;
					__TOL_[10] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[10] = fabs((_prvt_dtime/2) * (__K1_[10] - __K2_[10])/__TOL_[10]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[15])
				{
					__K2_[36]= 0.0000000000e+00;
					_prvt_aux_tol = fabs(__OLD_[36])*_prvt_rel_tol_;
					__TOL_[36] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[36] = fabs((_prvt_dtime/2) * (__K1_[36] - __K2_[36])/__TOL_[36]);
				}
				_prvt_time_new -= _prvt_dtime;
				
// 				timeval t3, t4;
// 				gettimeofday(&t3, NULL);
				
				#pragma omp barrier
				
// 				gettimeofday(&t4, NULL);
				
		/*		#pragma omp critical
				{
				    int aux = (int)((int)t4.tv_usec - (int)t3.tv_usec);
				    tempoespera += (aux<0)?0:aux;
				    
				}
		*/		
				
				///adapt the time step
				double __greatestError_=0.0;
				for(int k=0;k<numEDO;k++){
						if(__ERROR_[k] > __greatestError_)
							__greatestError_ = __ERROR_[k];
				}
				__greatestError_ += __tiny_;
				int flag = this->getErrorCode(__greatestError_, 1.0);
				_prvt_previous_dt = _prvt_dtime;
				//it doesn't accept the solution and cut h in a half
				if(flag==-1){
					//throw the results away and compute again
					_prvt_dtime = _prvt_dtime/_DECREASE_DT_2_; //cut time step in a half
					if(omp_get_thread_num()==_prvt_tree_thread[0])
					{
						__NEW_[0] = __K1_[0] * _prvt_dtime + __OLD_AUX_[0];
						__NEW_[1] = __K1_[1] * _prvt_dtime + __OLD_AUX_[1];
						__NEW_[2] = __K1_[2] * _prvt_dtime + __OLD_AUX_[2];
						__NEW_[3] = __K1_[3] * _prvt_dtime + __OLD_AUX_[3];
						__NEW_[4] = __K1_[4] * _prvt_dtime + __OLD_AUX_[4];
						__NEW_[5] = __K1_[5] * _prvt_dtime + __OLD_AUX_[5];
						__NEW_[18] = __K1_[18] * _prvt_dtime + __OLD_AUX_[18];
						__NEW_[27] = __K1_[27] * _prvt_dtime + __OLD_AUX_[27];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[8] = __K1_[8] * _prvt_dtime + __OLD_AUX_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[11] = __K1_[11] * _prvt_dtime + __OLD_AUX_[11];
						__NEW_[12] = __K1_[12] * _prvt_dtime + __OLD_AUX_[12];
						__NEW_[13] = __K1_[13] * _prvt_dtime + __OLD_AUX_[13];
						__NEW_[14] = __K1_[14] * _prvt_dtime + __OLD_AUX_[14];
						__NEW_[15] = __K1_[15] * _prvt_dtime + __OLD_AUX_[15];
						__NEW_[16] = __K1_[16] * _prvt_dtime + __OLD_AUX_[16];
						__NEW_[17] = __K1_[17] * _prvt_dtime + __OLD_AUX_[17];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[19] = __K1_[19] * _prvt_dtime + __OLD_AUX_[19];
						__NEW_[20] = __K1_[20] * _prvt_dtime + __OLD_AUX_[20];
						__NEW_[21] = __K1_[21] * _prvt_dtime + __OLD_AUX_[21];
						__NEW_[22] = __K1_[22] * _prvt_dtime + __OLD_AUX_[22];
						__NEW_[23] = __K1_[23] * _prvt_dtime + __OLD_AUX_[23];
						__NEW_[24] = __K1_[24] * _prvt_dtime + __OLD_AUX_[24];
						__NEW_[25] = __K1_[25] * _prvt_dtime + __OLD_AUX_[25];
						__NEW_[26] = __K1_[26] * _prvt_dtime + __OLD_AUX_[26];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[28] = __K1_[28] * _prvt_dtime + __OLD_AUX_[28];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[29] = __K1_[29] * _prvt_dtime + __OLD_AUX_[29];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[30] = __K1_[30] * _prvt_dtime + __OLD_AUX_[30];
						__NEW_[33] = __K1_[33] * _prvt_dtime + __OLD_AUX_[33];
						__NEW_[35] = __K1_[35] * _prvt_dtime + __OLD_AUX_[35];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[31] = __K1_[31] * _prvt_dtime + __OLD_AUX_[31];
						__NEW_[34] = __K1_[34] * _prvt_dtime + __OLD_AUX_[34];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[32] = __K1_[32] * _prvt_dtime + __OLD_AUX_[32];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[38] = __K1_[38] * _prvt_dtime + __OLD_AUX_[38];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[37] = __K1_[37] * _prvt_dtime + __OLD_AUX_[37];
						__NEW_[39] = __K1_[39] * _prvt_dtime + __OLD_AUX_[39];
						__NEW_[40] = __K1_[40] * _prvt_dtime + __OLD_AUX_[40];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[6] = __K1_[6] * _prvt_dtime + __OLD_AUX_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[7] = __K1_[7] * _prvt_dtime + __OLD_AUX_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[13])
					{
						__NEW_[9] = __K1_[9] * _prvt_dtime + __OLD_AUX_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[14])
					{
						__NEW_[10] = __K1_[10] * _prvt_dtime + __OLD_AUX_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[15])
					{
						__NEW_[36] = __K1_[36] * _prvt_dtime + __OLD_AUX_[36];
					}
				__TEMP_ = __NEW_;
				__NEW_ = __OLD_;
				__OLD_ = __TEMP_;
				}else{//it accepts the solutions
					/*if(_prvt_savingRate!=0){
						#pragma omp single
						{
							this->V_old_ = __OLD_AUX_[0];
							this->V_new_ = __OLD_[0];
							this->Cai_old_ = __OLD_AUX_[1];
							this->Cai_new_ = __OLD_[1];
							this->Cass_old_ = __OLD_AUX_[2];
							this->Cass_new_ = __OLD_[2];
							this->CaJSR_old_ = __OLD_AUX_[3];
							this->CaJSR_new_ = __OLD_[3];
							this->CaNSR_old_ = __OLD_AUX_[4];
							this->CaNSR_new_ = __OLD_[4];
							this->P_RyR_old_ = __OLD_AUX_[5];
							this->P_RyR_new_ = __OLD_[5];
							this->LTRPN_Ca_old_ = __OLD_AUX_[6];
							this->LTRPN_Ca_new_ = __OLD_[6];
							this->HTRPN_Ca_old_ = __OLD_AUX_[7];
							this->HTRPN_Ca_new_ = __OLD_[7];
							this->P_O1_old_ = __OLD_AUX_[8];
							this->P_O1_new_ = __OLD_[8];
							this->P_O2_old_ = __OLD_AUX_[9];
							this->P_O2_new_ = __OLD_[9];
							this->P_C2_old_ = __OLD_AUX_[10];
							this->P_C2_new_ = __OLD_[10];
							this->O_old_ = __OLD_AUX_[11];
							this->O_new_ = __OLD_[11];
							this->C2_old_ = __OLD_AUX_[12];
							this->C2_new_ = __OLD_[12];
							this->C3_old_ = __OLD_AUX_[13];
							this->C3_new_ = __OLD_[13];
							this->C4_old_ = __OLD_AUX_[14];
							this->C4_new_ = __OLD_[14];
							this->I1_old_ = __OLD_AUX_[15];
							this->I1_new_ = __OLD_[15];
							this->I2_old_ = __OLD_AUX_[16];
							this->I2_new_ = __OLD_[16];
							this->I3_old_ = __OLD_AUX_[17];
							this->I3_new_ = __OLD_[17];
							this->Nai_old_ = __OLD_AUX_[18];
							this->Nai_new_ = __OLD_[18];
							this->C_Na2_old_ = __OLD_AUX_[19];
							this->C_Na2_new_ = __OLD_[19];
							this->C_Na1_old_ = __OLD_AUX_[20];
							this->C_Na1_new_ = __OLD_[20];
							this->O_Na_old_ = __OLD_AUX_[21];
							this->O_Na_new_ = __OLD_[21];
							this->IF_Na_old_ = __OLD_AUX_[22];
							this->IF_Na_new_ = __OLD_[22];
							this->I1_Na_old_ = __OLD_AUX_[23];
							this->I1_Na_new_ = __OLD_[23];
							this->I2_Na_old_ = __OLD_AUX_[24];
							this->I2_Na_new_ = __OLD_[24];
							this->IC_Na2_old_ = __OLD_AUX_[25];
							this->IC_Na2_new_ = __OLD_[25];
							this->IC_Na3_old_ = __OLD_AUX_[26];
							this->IC_Na3_new_ = __OLD_[26];
							this->Ki_old_ = __OLD_AUX_[27];
							this->Ki_new_ = __OLD_[27];
							this->ato_f_old_ = __OLD_AUX_[28];
							this->ato_f_new_ = __OLD_[28];
							this->ito_f_old_ = __OLD_AUX_[29];
							this->ito_f_new_ = __OLD_[29];
							this->ato_s_old_ = __OLD_AUX_[30];
							this->ato_s_new_ = __OLD_[30];
							this->ito_s_old_ = __OLD_AUX_[31];
							this->ito_s_new_ = __OLD_[31];
							this->nKs_old_ = __OLD_AUX_[32];
							this->nKs_new_ = __OLD_[32];
							this->aur_old_ = __OLD_AUX_[33];
							this->aur_new_ = __OLD_[33];
							this->iur_old_ = __OLD_AUX_[34];
							this->iur_new_ = __OLD_[34];
							this->aKss_old_ = __OLD_AUX_[35];
							this->aKss_new_ = __OLD_[35];
							this->iKss_old_ = __OLD_AUX_[36];
							this->iKss_new_ = __OLD_[36];
							this->C_K2_old_ = __OLD_AUX_[37];
							this->C_K2_new_ = __OLD_[37];
							this->C_K1_old_ = __OLD_AUX_[38];
							this->C_K1_new_ = __OLD_[38];
							this->O_K_old_ = __OLD_AUX_[39];
							this->O_K_new_ = __OLD_[39];
							this->I_K_old_ = __OLD_AUX_[40];
							this->I_K_new_ = __OLD_[40];
							this->previous_dt = _prvt_previous_dt;
							this->dtime = _prvt_dtime;
							this->time_new = _prvt_time_new;
							save_step(fileptr, _ADAP_DT_);
						}
					}*/
					if(flag==3){
						_prvt_dtime = _prvt_dtime*_DECREASE_DT_;
					}else if(flag==4){
						_prvt_dtime = _prvt_dtime*_INCREASE_DT_;
					}else if(flag==0){
						//it just doesnt do anything
					}else{
						printf("flag: %d\n", flag);
					}
					if(_prvt_dtime > _prvt_maxStep && _prvt_maxStep!=0){
						_prvt_dtime = _prvt_maxStep;
					}else if(_prvt_dtime==0){
						printf("Error: Time step is zero.\n");
						break;
					}
					//it steps the method ahead, with euler solution
					if(omp_get_thread_num()==_prvt_tree_thread[0])
					{
						__NEW_[0] = __K2_[0] * _prvt_dtime + __OLD_[0];
						__NEW_[1] = __K2_[1] * _prvt_dtime + __OLD_[1];
						__NEW_[2] = __K2_[2] * _prvt_dtime + __OLD_[2];
						__NEW_[3] = __K2_[3] * _prvt_dtime + __OLD_[3];
						__NEW_[4] = __K2_[4] * _prvt_dtime + __OLD_[4];
						__NEW_[5] = __K2_[5] * _prvt_dtime + __OLD_[5];
						__NEW_[18] = __K2_[18] * _prvt_dtime + __OLD_[18];
						__NEW_[27] = __K2_[27] * _prvt_dtime + __OLD_[27];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[8] = __K2_[8] * _prvt_dtime + __OLD_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[11] = __K2_[11] * _prvt_dtime + __OLD_[11];
						__NEW_[12] = __K2_[12] * _prvt_dtime + __OLD_[12];
						__NEW_[13] = __K2_[13] * _prvt_dtime + __OLD_[13];
						__NEW_[14] = __K2_[14] * _prvt_dtime + __OLD_[14];
						__NEW_[15] = __K2_[15] * _prvt_dtime + __OLD_[15];
						__NEW_[16] = __K2_[16] * _prvt_dtime + __OLD_[16];
						__NEW_[17] = __K2_[17] * _prvt_dtime + __OLD_[17];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[19] = __K2_[19] * _prvt_dtime + __OLD_[19];
						__NEW_[20] = __K2_[20] * _prvt_dtime + __OLD_[20];
						__NEW_[21] = __K2_[21] * _prvt_dtime + __OLD_[21];
						__NEW_[22] = __K2_[22] * _prvt_dtime + __OLD_[22];
						__NEW_[23] = __K2_[23] * _prvt_dtime + __OLD_[23];
						__NEW_[24] = __K2_[24] * _prvt_dtime + __OLD_[24];
						__NEW_[25] = __K2_[25] * _prvt_dtime + __OLD_[25];
						__NEW_[26] = __K2_[26] * _prvt_dtime + __OLD_[26];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[28] = __K2_[28] * _prvt_dtime + __OLD_[28];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[29] = __K2_[29] * _prvt_dtime + __OLD_[29];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[30] = __K2_[30] * _prvt_dtime + __OLD_[30];
						__NEW_[33] = __K2_[33] * _prvt_dtime + __OLD_[33];
						__NEW_[35] = __K2_[35] * _prvt_dtime + __OLD_[35];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[31] = __K2_[31] * _prvt_dtime + __OLD_[31];
						__NEW_[34] = __K2_[34] * _prvt_dtime + __OLD_[34];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[32] = __K2_[32] * _prvt_dtime + __OLD_[32];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[38] = __K2_[38] * _prvt_dtime + __OLD_[38];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[37] = __K2_[37] * _prvt_dtime + __OLD_[37];
						__NEW_[39] = __K2_[39] * _prvt_dtime + __OLD_[39];
						__NEW_[40] = __K2_[40] * _prvt_dtime + __OLD_[40];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[6] = __K2_[6] * _prvt_dtime + __OLD_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[7] = __K2_[7] * _prvt_dtime + __OLD_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[13])
					{
						__NEW_[9] = __K2_[9] * _prvt_dtime + __OLD_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[14])
					{
						__NEW_[10] = __K2_[10] * _prvt_dtime + __OLD_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[15])
					{
						__NEW_[36] = __K2_[36] * _prvt_dtime + __OLD_[36];
					}
					//store the old iteration in a aux 
					__TEMP_ = __OLD_;
					__OLD_ = __OLD_AUX_;
					__OLD_AUX_ = __TEMP_;
					//steps ahead with euler
					__TEMP_ = __NEW_;
					__NEW_ = __OLD_;
					__OLD_ = __TEMP_;
					//change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração
					__TEMP_	= __K2_;
					__K2_	= __K1_;
					__K1_	= __TEMP_;
					//sums the old dtime - the variable dtime is alreaady updated
					_prvt_time_new += _prvt_previous_dt;
				}//FIM ELSE
				/*timeval t5, t6;
				gettimeofday(&t5, NULL);*/
				#pragma omp barrier
				/*gettimeofday(&t6, NULL);
				#pragma omp critical
				{
				    int aux = (int)((int)t6.tv_usec - (int)t5.tv_usec);
				    tempoespera += (aux<0)?0:aux;
				}*/
			}
// 			printf("%.1f\n", (float)tempoespera/(float)(nThreads));
		}
		
	}
	void Solveode::addt2_OMP(double finalTime, FILE *fileptr, int nThreads){
		omp_set_num_threads(nThreads);
		double *__NEW_ = (double*)malloc(sizeof(double)*numEDO);
		double *__OLD_ = (double*)malloc(sizeof(double)*numEDO);
		double *__OLD_AUX_ = (double*)malloc(sizeof(double)*numEDO);
		double *__TOL_ = (double*)malloc(sizeof(double)*numEDO);
		double *__ERROR_ = (double*)malloc(sizeof(double)*numEDO);
		double *__K1_  = (double*)malloc(sizeof(double)*numEDO);
		double *__K2_  = (double*)malloc(sizeof(double)*numEDO);
		double *__TEMP_;
		this->time_new = this->time;
		this->timeSaving = this->time;
		this->previous_dt = this->dtime;
		#pragma omp parallel firstprivate(__NEW_, __OLD_, __TEMP_, __TOL_, __K1_,__K2_, __ERROR_,__OLD_AUX_)
		{
			const double _beta_safety_ = 0.8;
			int *_prvt_tree_thread = tree_thread;
			double _prvt_rel_tol_=reltol__, _prvt_abs_tol_ =abstol__, _prvt_aux_tol;
			const double __tiny_ = pow(_prvt_abs_tol_, 2.0);
			//private parameters
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -8.0000000000e+01,  _prvt_stim_start = 2.0000000000e+01,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 7.1430000000e+01,  _prvt_stim_duration = 1.5000000000e+00,  _prvt_Acap = 1.5340000000e-04,  _prvt_Cm = 1.0000000000e+00,  _prvt_Vmyo = 2.5840000000e-05,  _prvt_F = 9.6500000000e+01,  _prvt_VJSR = 1.2000000000e-07,  _prvt_Vss = 1.4850000000e-09,  _prvt_VNSR = 2.0980000000e-06,  _prvt_CMDN_tot = 5.0000000000e+01,  _prvt_Km_CMDN = 2.3800000000e-01,  _prvt_CSQN_tot = 1.5000000000e+04,  _prvt_Km_CSQN = 8.0000000000e+02,  _prvt_v1 = 4.5000000000e+00,  _prvt_tau_tr = 2.0000000000e+01,  _prvt_tau_xfer = 8.0000000000e+00,  _prvt_v2 = 1.7400000000e-05,  _prvt_v3 = 4.5000000000e-01,  _prvt_Km_up = 5.0000000000e-01,  _prvt_k_plus_htrpn = 2.3700000000e-03,  _prvt_HTRPN_tot = 1.4000000000e+02,  _prvt_k_plus_ltrpn = 3.2700000000e-02,  _prvt_LTRPN_tot = 7.0000000000e+01,  _prvt_k_minus_htrpn = 3.2000000000e-05,  _prvt_k_minus_ltrpn = 1.9600000000e-02,  _prvt_i_CaL_max = 7.0000000000e+00,  _prvt_k_plus_a = 6.0750000000e-03,  _prvt_n = 4.0000000000e+00,  _prvt_k_minus_b = 9.6500000000e-01,  _prvt_k_minus_c = 8.0000000000e-04,  _prvt_k_minus_a = 7.1250000000e-02,  _prvt_k_plus_b = 4.0500000000e-03,  _prvt_m = 3.0000000000e+00,  _prvt_k_plus_c = 9.0000000000e-03,  _prvt_g_CaL = 1.7290000000e-01,  _prvt_E_CaL = 6.3000000000e+01,  _prvt_Kpcb = 5.0000000000e-04,  _prvt_Kpc_max = 2.3324000000e-01,  _prvt_Kpc_half = 2.0000000000e+01,  _prvt_i_pCa_max = 1.0000000000e+00,  _prvt_Km_pCa = 5.0000000000e-01,  _prvt_k_NaCa = 2.9280000000e+02,  _prvt_K_mNa = 8.7500000000e+04,  _prvt_Nao = 1.4000000000e+05,  _prvt_K_mCa = 1.3800000000e+03,  _prvt_Cao = 1.8000000000e+03,  _prvt_k_sat = 1.0000000000e-01,  _prvt_eta = 3.5000000000e-01,  _prvt_R = 8.3140000000e+00,  _prvt_T = 2.9800000000e+02,  _prvt_g_Cab = 3.6700000000e-04,  _prvt_g_Na = 1.3000000000e+01,  _prvt_Ko = 5.4000000000e+03,  _prvt_g_Nab = 2.6000000000e-03,  _prvt_g_Kto_f = 4.0670000000e-01,  _prvt_g_Kto_s = 0.0000000000e+00,  _prvt_g_Ks = 5.7500000000e-03,  _prvt_g_Kur = 1.6000000000e-01,  _prvt_g_Kss = 5.0000000000e-02,  _prvt_g_Kr = 7.8000000000e-02,  _prvt_kf = 2.3761000000e-02,  _prvt_kb = 3.6778000000e-02,  _prvt_i_NaK_max = 8.8000000000e-01,  _prvt_Km_Nai = 2.1000000000e+04,  _prvt_Km_Ko = 1.5000000000e+03,  _prvt_g_ClCa = 1.0000000000e+01,  _prvt_Km_Cl = 1.0000000000e+01,  _prvt_E_Cl = -4.0000000000e+01, 
			//private aux variables
			 _prvt_calc_i_stim=0.0,  _prvt_calc_Bi=0.0,  _prvt_calc_Bss=0.0,  _prvt_calc_BJSR=0.0,  _prvt_calc_J_rel=0.0,  _prvt_calc_J_tr=0.0,  _prvt_calc_J_xfer=0.0,  _prvt_calc_J_leak=0.0,  _prvt_calc_J_up=0.0,  _prvt_calc_J_trpn=0.0,  _prvt_calc_P_C1=0.0,  _prvt_calc_i_CaL=0.0,  _prvt_calc_C1=0.0,  _prvt_calc_alpha=0.0,  _prvt_calc_beta=0.0,  _prvt_calc_gamma=0.0,  _prvt_calc_Kpcf=0.0,  _prvt_calc_i_pCa=0.0,  _prvt_calc_i_NaCa=0.0,  _prvt_calc_i_Cab=0.0,  _prvt_calc_E_CaN=0.0,  _prvt_calc_i_Na=0.0,  _prvt_calc_E_Na=0.0,  _prvt_calc_C_Na3=0.0,  _prvt_calc_alpha_Na11=0.0,  _prvt_calc_alpha_Na12=0.0,  _prvt_calc_alpha_Na13=0.0,  _prvt_calc_beta_Na11=0.0,  _prvt_calc_beta_Na12=0.0,  _prvt_calc_beta_Na13=0.0,  _prvt_calc_alpha_Na3=0.0,  _prvt_calc_beta_Na3=0.0,  _prvt_calc_alpha_Na2=0.0,  _prvt_calc_beta_Na2=0.0,  _prvt_calc_alpha_Na4=0.0,  _prvt_calc_beta_Na4=0.0,  _prvt_calc_alpha_Na5=0.0,  _prvt_calc_beta_Na5=0.0,  _prvt_calc_i_Nab=0.0,  _prvt_calc_i_Kto_f=0.0,  _prvt_calc_E_K=0.0,  _prvt_calc_alpha_a=0.0,  _prvt_calc_beta_a=0.0,  _prvt_calc_alpha_i=0.0,  _prvt_calc_beta_i=0.0,  _prvt_calc_i_Kto_s=0.0,  _prvt_calc_ass=0.0,  _prvt_calc_iss=0.0,  _prvt_calc_tau_ta_s=0.0,  _prvt_calc_tau_ti_s=0.0,  _prvt_calc_i_K1=0.0,  _prvt_calc_i_Ks=0.0,  _prvt_calc_alpha_n=0.0,  _prvt_calc_beta_n=0.0,  _prvt_calc_i_Kur=0.0,  _prvt_calc_tau_aur=0.0,  _prvt_calc_tau_iur=0.0,  _prvt_calc_i_Kss=0.0,  _prvt_calc_tau_Kss=0.0,  _prvt_calc_i_Kr=0.0,  _prvt_calc_C_K0=0.0,  _prvt_calc_alpha_a0=0.0,  _prvt_calc_beta_a0=0.0,  _prvt_calc_alpha_a1=0.0,  _prvt_calc_beta_a1=0.0,  _prvt_calc_i_NaK=0.0,  _prvt_calc_f_NaK=0.0,  _prvt_calc_sigma=0.0,  _prvt_calc_i_ClCa=0.0,  _prvt_calc_O_ClCa=0.0,  _prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current=0.0,  _prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current=0.0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_Cai_lado_direito_,  _prvt_Cass_lado_direito_,  _prvt_CaJSR_lado_direito_,  _prvt_CaNSR_lado_direito_,  _prvt_P_RyR_lado_direito_,  _prvt_LTRPN_Ca_lado_direito_,  _prvt_HTRPN_Ca_lado_direito_,  _prvt_P_O1_lado_direito_,  _prvt_P_O2_lado_direito_,  _prvt_P_C2_lado_direito_,  _prvt_O_lado_direito_,  _prvt_C2_lado_direito_,  _prvt_C3_lado_direito_,  _prvt_C4_lado_direito_,  _prvt_I1_lado_direito_,  _prvt_I2_lado_direito_,  _prvt_I3_lado_direito_,  _prvt_Nai_lado_direito_,  _prvt_C_Na2_lado_direito_,  _prvt_C_Na1_lado_direito_,  _prvt_O_Na_lado_direito_,  _prvt_IF_Na_lado_direito_,  _prvt_I1_Na_lado_direito_,  _prvt_I2_Na_lado_direito_,  _prvt_IC_Na2_lado_direito_,  _prvt_IC_Na3_lado_direito_,  _prvt_Ki_lado_direito_,  _prvt_ato_f_lado_direito_,  _prvt_ito_f_lado_direito_,  _prvt_ato_s_lado_direito_,  _prvt_ito_s_lado_direito_,  _prvt_nKs_lado_direito_,  _prvt_aur_lado_direito_,  _prvt_iur_lado_direito_,  _prvt_aKss_lado_direito_,  _prvt_iKss_lado_direito_,  _prvt_C_K2_lado_direito_,  _prvt_C_K1_lado_direito_,  _prvt_O_K_lado_direito_,  _prvt_I_K_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;
			__NEW_[0] = __OLD_[0] = -8.2420200000e+01;
			__NEW_[1] = __OLD_[1] = 1.1500100000e-01;
			__NEW_[2] = __OLD_[2] = 1.1500100000e-01;
			__NEW_[3] = __OLD_[3] = 1.2995000000e+03;
			__NEW_[4] = __OLD_[4] = 1.2995000000e+03;
			__NEW_[5] = __OLD_[5] = 0.0000000000e+00;
			__NEW_[6] = __OLD_[6] = 1.1268400000e+01;
			__NEW_[7] = __OLD_[7] = 1.2529000000e+02;
			__NEW_[8] = __OLD_[8] = 1.4910200000e-05;
			__NEW_[9] = __OLD_[9] = 9.5172600000e-11;
			__NEW_[10] = __OLD_[10] = 1.6774000000e-04;
			__NEW_[11] = __OLD_[11] = 9.3030800000e-19;
			__NEW_[12] = __OLD_[12] = 1.2421600000e-04;
			__NEW_[13] = __OLD_[13] = 5.7867900000e-09;
			__NEW_[14] = __OLD_[14] = 1.1981600000e-13;
			__NEW_[15] = __OLD_[15] = 4.9792300000e-19;
			__NEW_[16] = __OLD_[16] = 3.4584700000e-14;
			__NEW_[17] = __OLD_[17] = 1.8510600000e-14;
			__NEW_[18] = __OLD_[18] = 1.4237100000e+04;
			__NEW_[19] = __OLD_[19] = 2.0752000000e-02;
			__NEW_[20] = __OLD_[20] = 2.7913200000e-04;
			__NEW_[21] = __OLD_[21] = 7.1348300000e-07;
			__NEW_[22] = __OLD_[22] = 1.5317600000e-04;
			__NEW_[23] = __OLD_[23] = 6.7334500000e-07;
			__NEW_[24] = __OLD_[24] = 1.5578700000e-09;
			__NEW_[25] = __OLD_[25] = 1.1387900000e-02;
			__NEW_[26] = __OLD_[26] = 3.4278000000e-01;
			__NEW_[27] = __OLD_[27] = 1.4372000000e+05;
			__NEW_[28] = __OLD_[28] = 2.6556300000e-03;
			__NEW_[29] = __OLD_[29] = 9.9997700000e-01;
			__NEW_[30] = __OLD_[30] = 4.1706900000e-04;
			__NEW_[31] = __OLD_[31] = 9.9854300000e-01;
			__NEW_[32] = __OLD_[32] = 2.6275300000e-04;
			__NEW_[33] = __OLD_[33] = 4.1706900000e-04;
			__NEW_[34] = __OLD_[34] = 9.9854300000e-01;
			__NEW_[35] = __OLD_[35] = 4.1706900000e-04;
			__NEW_[36] = __OLD_[36] = 1.0000000000e+00;
			__NEW_[37] = __OLD_[37] = 6.4122900000e-04;
			__NEW_[38] = __OLD_[38] = 9.9251300000e-04;
			__NEW_[39] = __OLD_[39] = 1.7529800000e-04;
			__NEW_[40] = __OLD_[40] = 3.1912900000e-05;
			_prvt_time_new += _prvt_dtime;
			if(omp_get_thread_num()==_prvt_tree_thread[0])
			{
				_prvt_calc_i_stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
				_prvt_calc_Bi = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[1]),2.0000000000e+00))),(-1.0000000000e+00));
				_prvt_calc_Bss = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[2]),2.0000000000e+00))),(-1.0000000000e+00));
				_prvt_calc_BJSR = pow((1.0000000000e+00+((_prvt_CSQN_tot*_prvt_Km_CSQN)/pow((_prvt_Km_CSQN+__OLD_[3]),2.0000000000e+00))),(-1.0000000000e+00));
				_prvt_calc_J_rel = (_prvt_v1*(__OLD_[8]+__OLD_[9])*(__OLD_[3]-__OLD_[2])*__OLD_[5]);
				_prvt_calc_J_tr = ((__OLD_[4]-__OLD_[3])/_prvt_tau_tr);
				_prvt_calc_J_xfer = ((__OLD_[2]-__OLD_[1])/_prvt_tau_xfer);
				_prvt_calc_J_leak = (_prvt_v2*(__OLD_[4]-__OLD_[1]));
				_prvt_calc_J_up = ((_prvt_v3*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_up,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
				_prvt_calc_J_trpn = (((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))+(_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6])))-((_prvt_k_minus_htrpn*__OLD_[7])+(_prvt_k_minus_ltrpn*__OLD_[6])));
				_prvt_calc_i_CaL = (_prvt_g_CaL*__OLD_[11]*(__OLD_[0]-_prvt_E_CaL));
				_prvt_calc_i_pCa = ((_prvt_i_pCa_max*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_pCa,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
				_prvt_calc_i_NaCa = (((((((_prvt_k_NaCa*1.0000000000e+00)/(pow(_prvt_K_mNa,3.0000000000e+00)+pow(_prvt_Nao,3.0000000000e+00)))*1.0000000000e+00)/(_prvt_K_mCa+_prvt_Cao))*1.0000000000e+00)/(1.0000000000e+00+(_prvt_k_sat*exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))))*((exp(((_prvt_eta*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[18],3.0000000000e+00)*_prvt_Cao)-(exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Nao,3.0000000000e+00)*__OLD_[1])));
				_prvt_calc_E_CaN = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Cao/__OLD_[1])));
				_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((((9.0000000000e-01*_prvt_Nao)+(1.0000000000e-01*_prvt_Ko))/((9.0000000000e-01*__OLD_[18])+(1.0000000000e-01*__OLD_[27])))));
				_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ko/__OLD_[27])));
				_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[39]*(__OLD_[0]-(((_prvt_R*_prvt_T)/_prvt_F)*log((((9.8000000000e-01*_prvt_Ko)+(2.0000000000e-02*_prvt_Nao))/((9.8000000000e-01*__OLD_[27])+(2.0000000000e-02*__OLD_[18])))))));
				_prvt_calc_sigma = ((1.0000000000e+00/7.0000000000e+00)*(exp((_prvt_Nao/6.7300000000e+04))-1.0000000000e+00));
				_prvt_calc_O_ClCa = (2.0000000000e-01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.6700000000e+01))/7.8000000000e+00))));
				_prvt_calc_i_Nab = (_prvt_g_Nab*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_Kto_s = (_prvt_g_Kto_s*__OLD_[30]*__OLD_[31]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_K1 = ((((2.9380000000e-01*_prvt_Ko)/(_prvt_Ko+2.1000000000e+02))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp((8.9600000000e-02*(__OLD_[0]-_prvt_calc_E_K)))));
				_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[32],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_Kur = (_prvt_g_Kur*__OLD_[33]*__OLD_[34]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_Kss = (_prvt_g_Kss*__OLD_[35]*__OLD_[36]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_Cab = (_prvt_g_Cab*(__OLD_[0]-_prvt_calc_E_CaN));
				_prvt_calc_i_Na = (_prvt_g_Na*__OLD_[21]*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_Kto_f = (_prvt_g_Kto_f*pow(__OLD_[28],3.0000000000e+00)*__OLD_[29]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_f_NaK = (1.0000000000e+00/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.6500000000e-02*_prvt_calc_sigma*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_ClCa = (((_prvt_g_ClCa*_prvt_calc_O_ClCa*__OLD_[1])/(__OLD_[1]+_prvt_Km_Cl))*(__OLD_[0]-_prvt_E_Cl));
				_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_calc_f_NaK*1.0000000000e+00)/(1.0000000000e+00+pow((_prvt_Km_Nai/__OLD_[18]),1.5000000000e+00)))*_prvt_Ko)/(_prvt_Ko+_prvt_Km_Ko));
				__K1_[0]= (-(_prvt_calc_i_CaL+_prvt_calc_i_pCa+_prvt_calc_i_NaCa+_prvt_calc_i_Cab+_prvt_calc_i_Na+_prvt_calc_i_Nab+_prvt_calc_i_NaK+_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kur+_prvt_calc_i_Kss+_prvt_calc_i_Kr+_prvt_calc_i_ClCa+_prvt_calc_i_stim));
				__NEW_[0]= __K1_[0] * _prvt_dtime + __OLD_[0];
				__K1_[1]= (_prvt_calc_Bi*((_prvt_calc_J_leak+_prvt_calc_J_xfer)-(_prvt_calc_J_up+_prvt_calc_J_trpn+((((_prvt_calc_i_Cab+_prvt_calc_i_pCa)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vmyo*_prvt_F)))));
				__NEW_[1]= __K1_[1] * _prvt_dtime + __OLD_[1];
				__K1_[2]= (_prvt_calc_Bss*(((_prvt_calc_J_rel*_prvt_VJSR)/_prvt_Vss)-(((_prvt_calc_J_xfer*_prvt_Vmyo)/_prvt_Vss)+((_prvt_calc_i_CaL*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vss*_prvt_F)))));
				__NEW_[2]= __K1_[2] * _prvt_dtime + __OLD_[2];
				__K1_[3]= (_prvt_calc_BJSR*(_prvt_calc_J_tr-_prvt_calc_J_rel));
				__NEW_[3]= __K1_[3] * _prvt_dtime + __OLD_[3];
				__K1_[4]= ((((_prvt_calc_J_up-_prvt_calc_J_leak)*_prvt_Vmyo)/_prvt_VNSR)-((_prvt_calc_J_tr*_prvt_VJSR)/_prvt_VNSR));
				__NEW_[4]= __K1_[4] * _prvt_dtime + __OLD_[4];
				__K1_[5]= (((-4.0000000000e-02)*__OLD_[5])-(((1.0000000000e-01*_prvt_calc_i_CaL)/_prvt_i_CaL_max)*exp(((-pow((__OLD_[0]-5.0000000000e+00),2.0000000000e+00))/6.4800000000e+02))));
				__NEW_[5]= __K1_[5] * _prvt_dtime + __OLD_[5];
				__K1_[18]= (((-(_prvt_calc_i_Na+_prvt_calc_i_Nab+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
				__NEW_[18]= __K1_[18] * _prvt_dtime + __OLD_[18];
				__K1_[27]= (((-((_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kss+_prvt_calc_i_Kur+_prvt_calc_i_Kr)-(2.0000000000e+00*_prvt_calc_i_NaK)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
				__NEW_[27]= __K1_[27] * _prvt_dtime + __OLD_[27];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[1])
			{
				_prvt_calc_P_C1 = (1.0000000000e+00-(__OLD_[10]+__OLD_[8]+__OLD_[9]));
				__K1_[8]= (((_prvt_k_plus_a*pow(__OLD_[2],_prvt_n)*_prvt_calc_P_C1)+(_prvt_k_minus_b*__OLD_[9])+(_prvt_k_minus_c*__OLD_[10]))-((_prvt_k_minus_a*__OLD_[8])+(_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])+(_prvt_k_plus_c*__OLD_[8])));
				__NEW_[8]= __K1_[8] * _prvt_dtime + __OLD_[8];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[2])
			{
				_prvt_calc_C1 = (1.0000000000e+00-(__OLD_[11]+__OLD_[12]+__OLD_[12]+__OLD_[13]+__OLD_[14]+__OLD_[15]+__OLD_[16]+__OLD_[17]));
				_prvt_calc_alpha = ((4.0000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01))*((1.0000000000e+00+(7.0000000000e-01*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.0000000000e+01))))-(7.5000000000e-01*exp(((-pow((__OLD_[0]+2.0000000000e+01),2.0000000000e+00))/4.0000000000e+02)))))/(1.0000000000e+00+(1.2000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01)))));
				_prvt_calc_beta = (5.0000000000e-02*exp(((-(__OLD_[0]+1.2000000000e+01))/1.3000000000e+01)));
				_prvt_calc_gamma = ((_prvt_Kpc_max*__OLD_[2])/(_prvt_Kpc_half+__OLD_[2]));
				_prvt_calc_Kpcf = (1.3000000000e+01*(1.0000000000e+00-exp(((-pow((__OLD_[0]+1.4500000000e+01),2.0000000000e+00))/1.0000000000e+02))));
				__K1_[11]= (((_prvt_calc_alpha*__OLD_[14])+(_prvt_Kpcb*__OLD_[15])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[11]))))-((4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(_prvt_calc_gamma*__OLD_[11])));
				__NEW_[11]= __K1_[11] * _prvt_dtime + __OLD_[11];
				__K1_[12]= (((4.0000000000e+00*_prvt_calc_alpha*_prvt_calc_C1)+(2.0000000000e+00*_prvt_calc_beta*__OLD_[13]))-((_prvt_calc_beta*__OLD_[12])+(3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])));
				__NEW_[12]= __K1_[12] * _prvt_dtime + __OLD_[12];
				__K1_[13]= (((3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])+(3.0000000000e+00*_prvt_calc_beta*__OLD_[14]))-((2.0000000000e+00*_prvt_calc_beta*__OLD_[13])+(2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])));
				__NEW_[13]= __K1_[13] * _prvt_dtime + __OLD_[13];
				__K1_[14]= (((2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])+(4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(1.0000000000e-02*((4.0000000000e+00*_prvt_Kpcb*_prvt_calc_beta*__OLD_[15])-(_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])))+(2.0000000000e-03*((4.0000000000e+00*_prvt_calc_beta*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[14])))+(4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17]))-((3.0000000000e+00*_prvt_calc_beta*__OLD_[14])+(_prvt_calc_alpha*__OLD_[14])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14])));
				__NEW_[14]= __K1_[14] * _prvt_dtime + __OLD_[14];
				__K1_[15]= (((_prvt_calc_gamma*__OLD_[11])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[17])-(_prvt_calc_Kpcf*__OLD_[15])))+(1.0000000000e-02*((_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*_prvt_calc_Kpcf*__OLD_[15]))))-(_prvt_Kpcb*__OLD_[15]));
				__NEW_[15]= __K1_[15] * _prvt_dtime + __OLD_[15];
				__K1_[16]= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[11])-(_prvt_calc_alpha*__OLD_[16])))+(_prvt_Kpcb*__OLD_[17])+(2.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*__OLD_[16]))))-(_prvt_calc_gamma*__OLD_[16]));
				__NEW_[16]= __K1_[16] * _prvt_dtime + __OLD_[16];
				__K1_[17]= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[15])-(_prvt_calc_alpha*__OLD_[17])))+(_prvt_calc_gamma*__OLD_[16])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14]))-((4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17])+(_prvt_Kpcb*__OLD_[17])));
				__NEW_[17]= __K1_[17] * _prvt_dtime + __OLD_[17];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[3])
			{
				_prvt_calc_C_Na3 = (1.0000000000e+00-(__OLD_[21]+__OLD_[20]+__OLD_[19]+__OLD_[22]+__OLD_[23]+__OLD_[24]+__OLD_[25]+__OLD_[26]));
				_prvt_calc_alpha_Na11 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.7000000000e+01)))+(2.0000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
				_prvt_calc_alpha_Na12 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+01)))+(2.3000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
				_prvt_calc_alpha_Na13 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.2000000000e+01)))+(2.5000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
				_prvt_calc_beta_Na11 = (1.9170000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/2.0300000000e+01)));
				_prvt_calc_beta_Na12 = (2.0000000000e-01*exp(((-(__OLD_[0]-2.5000000000e+00))/2.0300000000e+01)));
				_prvt_calc_beta_Na13 = (2.2000000000e-01*exp(((-(__OLD_[0]-7.5000000000e+00))/2.0300000000e+01)));
				_prvt_calc_alpha_Na3 = (7.0000000000e-07*exp(((-(__OLD_[0]+7.0000000000e+00))/7.7000000000e+00)));
				_prvt_calc_beta_Na3 = (8.4000000000e-03+(2.0000000000e-05*(__OLD_[0]+7.0000000000e+00)));
				_prvt_calc_alpha_Na2 = (1.0000000000e+00/((1.8849500000e-01*exp(((-(__OLD_[0]+7.0000000000e+00))/1.6600000000e+01)))+3.9395600000e-01));
				_prvt_calc_beta_Na2 = ((_prvt_calc_alpha_Na13*_prvt_calc_alpha_Na2*_prvt_calc_alpha_Na3)/(_prvt_calc_beta_Na13*_prvt_calc_beta_Na3));
				_prvt_calc_alpha_Na4 = (_prvt_calc_alpha_Na2/1.0000000000e+03);
				_prvt_calc_beta_Na4 = _prvt_calc_alpha_Na3;
				_prvt_calc_alpha_Na5 = (_prvt_calc_alpha_Na2/9.5000000000e+04);
				_prvt_calc_beta_Na5 = (_prvt_calc_alpha_Na3/5.0000000000e+01);
				__K1_[19]= (((_prvt_calc_alpha_Na11*_prvt_calc_C_Na3)+(_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[19])+(_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na3*__OLD_[19])));
				__NEW_[19]= __K1_[19] * _prvt_dtime + __OLD_[19];
				__K1_[20]= (((_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na3*__OLD_[22]))-((_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na3*__OLD_[20])));
				__NEW_[20]= __K1_[20] * _prvt_dtime + __OLD_[20];
				__K1_[21]= (((_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na2*__OLD_[22]))-((_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na2*__OLD_[21])));
				__NEW_[21]= __K1_[21] * _prvt_dtime + __OLD_[21];
				__K1_[22]= (((_prvt_calc_alpha_Na2*__OLD_[21])+(_prvt_calc_beta_Na3*__OLD_[20])+(_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na12*__OLD_[25]))-((_prvt_calc_beta_Na2*__OLD_[22])+(_prvt_calc_alpha_Na3*__OLD_[22])+(_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na12*__OLD_[22])));
				__NEW_[22]= __K1_[22] * _prvt_dtime + __OLD_[22];
				__K1_[23]= (((_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na5*__OLD_[24]))-((_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na5*__OLD_[23])));
				__NEW_[23]= __K1_[23] * _prvt_dtime + __OLD_[23];
				__K1_[24]= ((_prvt_calc_alpha_Na5*__OLD_[23])-(_prvt_calc_beta_Na5*__OLD_[24]));
				__NEW_[24]= __K1_[24] * _prvt_dtime + __OLD_[24];
				__K1_[25]= (((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_beta_Na12*__OLD_[22])+(_prvt_calc_beta_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_alpha_Na12*__OLD_[25])+(_prvt_calc_alpha_Na3*__OLD_[25])));
				__NEW_[25]= __K1_[25] * _prvt_dtime + __OLD_[25];
				__K1_[26]= (((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_beta_Na3*_prvt_calc_C_Na3))-((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_alpha_Na3*__OLD_[26])));
				__NEW_[26]= __K1_[26] * _prvt_dtime + __OLD_[26];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[4])
			{
				_prvt_calc_alpha_a = (1.8064000000e-01*exp((3.5770000000e-02*(__OLD_[0]+3.0000000000e+01))));
				_prvt_calc_beta_a = (3.9560000000e-01*exp(((-6.2370000000e-02)*(__OLD_[0]+3.0000000000e+01))));
				__K1_[28]= ((_prvt_calc_alpha_a*(1.0000000000e+00-__OLD_[28]))-(_prvt_calc_beta_a*__OLD_[28]));
				__NEW_[28]= __K1_[28] * _prvt_dtime + __OLD_[28];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[5])
			{
				_prvt_calc_alpha_i = ((1.5200000000e-04*exp(((-(__OLD_[0]+1.3500000000e+01))/7.0000000000e+00)))/((6.7083000000e-03*exp(((-(__OLD_[0]+3.3500000000e+01))/7.0000000000e+00)))+1.0000000000e+00));
				_prvt_calc_beta_i = ((9.5000000000e-04*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))/((5.1335000000e-02*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))+1.0000000000e+00));
				__K1_[29]= ((_prvt_calc_alpha_i*(1.0000000000e+00-__OLD_[29]))-(_prvt_calc_beta_i*__OLD_[29]));
				__NEW_[29]= __K1_[29] * _prvt_dtime + __OLD_[29];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[6])
			{
				_prvt_calc_ass = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2500000000e+01))/7.7000000000e+00))));
				_prvt_calc_tau_ta_s = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
				_prvt_calc_tau_aur = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
				_prvt_calc_tau_Kss = ((3.9300000000e+01*exp(((-8.6200000000e-02)*__OLD_[0])))+1.3170000000e+01);
				__K1_[30]= ((_prvt_calc_ass-__OLD_[30])/_prvt_calc_tau_ta_s);
				__NEW_[30]= __K1_[30] * _prvt_dtime + __OLD_[30];
				__K1_[33]= ((_prvt_calc_ass-__OLD_[33])/_prvt_calc_tau_aur);
				__NEW_[33]= __K1_[33] * _prvt_dtime + __OLD_[33];
				__K1_[35]= ((_prvt_calc_ass-__OLD_[35])/_prvt_calc_tau_Kss);
				__NEW_[35]= __K1_[35] * _prvt_dtime + __OLD_[35];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[7])
			{
				_prvt_calc_iss = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00))));
				_prvt_calc_tau_ti_s = (2.7000000000e+02+(1.0500000000e+03/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
				_prvt_calc_tau_iur = (1.2000000000e+03-(1.7000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
				__K1_[31]= ((_prvt_calc_iss-__OLD_[31])/_prvt_calc_tau_ti_s);
				__NEW_[31]= __K1_[31] * _prvt_dtime + __OLD_[31];
				__K1_[34]= ((_prvt_calc_iss-__OLD_[34])/_prvt_calc_tau_iur);
				__NEW_[34]= __K1_[34] * _prvt_dtime + __OLD_[34];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[8])
			{
				_prvt_calc_alpha_n = ((4.8133300000e-06*(__OLD_[0]+2.6500000000e+01))/(1.0000000000e+00-exp(((-1.2800000000e-01)*(__OLD_[0]+2.6500000000e+01)))));
				_prvt_calc_beta_n = (9.5333300000e-05*exp(((-3.8000000000e-02)*(__OLD_[0]+2.6500000000e+01))));
				__K1_[32]= ((_prvt_calc_alpha_n*(1.0000000000e+00-__OLD_[32]))-(_prvt_calc_beta_n*__OLD_[32]));
				__NEW_[32]= __K1_[32] * _prvt_dtime + __OLD_[32];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[9])
			{
				_prvt_calc_C_K0 = (1.0000000000e+00-(__OLD_[38]+__OLD_[37]+__OLD_[39]+__OLD_[40]));
				_prvt_calc_alpha_a0 = (2.2348000000e-02*exp((1.1760000000e-02*__OLD_[0])));
				_prvt_calc_beta_a0 = (4.7002000000e-02*exp(((-6.3100000000e-02)*__OLD_[0])));
				__K1_[38]= (((_prvt_calc_alpha_a0*_prvt_calc_C_K0)+(_prvt_kb*__OLD_[37]))-((_prvt_calc_beta_a0*__OLD_[38])+(_prvt_kf*__OLD_[38])));
				__NEW_[38]= __K1_[38] * _prvt_dtime + __OLD_[38];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[10])
			{
				_prvt_calc_alpha_a1 = (1.3733000000e-02*exp((3.8198000000e-02*__OLD_[0])));
				_prvt_calc_beta_a1 = (6.8900000000e-05*exp(((-4.1780000000e-02)*__OLD_[0])));
				_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (9.0821000000e-02*exp((2.3391000000e-02*(__OLD_[0]+5.0000000000e+00))));
				_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (6.4970000000e-03*exp(((-3.2680000000e-02)*(__OLD_[0]+5.0000000000e+00))));
				__K1_[37]= (((_prvt_kf*__OLD_[38])+(_prvt_calc_beta_a1*__OLD_[39]))-((_prvt_kb*__OLD_[37])+(_prvt_calc_alpha_a1*__OLD_[37])));
				__NEW_[37]= __K1_[37] * _prvt_dtime + __OLD_[37];
				__K1_[39]= (((_prvt_calc_alpha_a1*__OLD_[37])+(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]))-((_prvt_calc_beta_a1*__OLD_[39])+(_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])));
				__NEW_[39]= __K1_[39] * _prvt_dtime + __OLD_[39];
				__K1_[40]= ((_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])-(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]));
				__NEW_[40]= __K1_[40] * _prvt_dtime + __OLD_[40];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[11])
			{
				__K1_[6]= ((_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6]))-(_prvt_k_minus_ltrpn*__OLD_[6]));
				__NEW_[6]= __K1_[6] * _prvt_dtime + __OLD_[6];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[12])
			{
				__K1_[7]= ((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))-(_prvt_k_minus_htrpn*__OLD_[7]));
				__NEW_[7]= __K1_[7] * _prvt_dtime + __OLD_[7];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[13])
			{
				__K1_[9]= ((_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])-(_prvt_k_minus_b*__OLD_[9]));
				__NEW_[9]= __K1_[9] * _prvt_dtime + __OLD_[9];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[14])
			{
				__K1_[10]= ((_prvt_k_plus_c*__OLD_[8])-(_prvt_k_minus_c*__OLD_[10]));
				__NEW_[10]= __K1_[10] * _prvt_dtime + __OLD_[10];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[15])
			{
				__K1_[36]= 0.0000000000e+00;
				__NEW_[36]= __K1_[36] * _prvt_dtime + __OLD_[36];
			}
			//store the old iteration in a aux 
			__TEMP_ = __OLD_;
			__OLD_ = __OLD_AUX_;
			__OLD_AUX_ = __TEMP_;
			//steps ahead with euler
			__TEMP_ = __NEW_;
			__NEW_ = __OLD_;
			__OLD_ = __TEMP_;
			//as threads devem começar o  laço ao mesmo tempo
			#pragma omp barrier
			while(_prvt_time_new<=_prvt_finalTime){
				_prvt_time_new += _prvt_dtime;
				if(omp_get_thread_num()==_prvt_tree_thread[0])
				{
					_prvt_calc_i_stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
					_prvt_calc_Bi = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[1]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_Bss = pow((1.0000000000e+00+((_prvt_CMDN_tot*_prvt_Km_CMDN)/pow((_prvt_Km_CMDN+__OLD_[2]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_BJSR = pow((1.0000000000e+00+((_prvt_CSQN_tot*_prvt_Km_CSQN)/pow((_prvt_Km_CSQN+__OLD_[3]),2.0000000000e+00))),(-1.0000000000e+00));
					_prvt_calc_J_rel = (_prvt_v1*(__OLD_[8]+__OLD_[9])*(__OLD_[3]-__OLD_[2])*__OLD_[5]);
					_prvt_calc_J_tr = ((__OLD_[4]-__OLD_[3])/_prvt_tau_tr);
					_prvt_calc_J_xfer = ((__OLD_[2]-__OLD_[1])/_prvt_tau_xfer);
					_prvt_calc_J_leak = (_prvt_v2*(__OLD_[4]-__OLD_[1]));
					_prvt_calc_J_up = ((_prvt_v3*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_up,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
					_prvt_calc_J_trpn = (((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))+(_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6])))-((_prvt_k_minus_htrpn*__OLD_[7])+(_prvt_k_minus_ltrpn*__OLD_[6])));
					_prvt_calc_i_CaL = (_prvt_g_CaL*__OLD_[11]*(__OLD_[0]-_prvt_E_CaL));
					_prvt_calc_i_pCa = ((_prvt_i_pCa_max*pow(__OLD_[1],2.0000000000e+00))/(pow(_prvt_Km_pCa,2.0000000000e+00)+pow(__OLD_[1],2.0000000000e+00)));
					_prvt_calc_i_NaCa = (((((((_prvt_k_NaCa*1.0000000000e+00)/(pow(_prvt_K_mNa,3.0000000000e+00)+pow(_prvt_Nao,3.0000000000e+00)))*1.0000000000e+00)/(_prvt_K_mCa+_prvt_Cao))*1.0000000000e+00)/(1.0000000000e+00+(_prvt_k_sat*exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))))*((exp(((_prvt_eta*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[18],3.0000000000e+00)*_prvt_Cao)-(exp((((_prvt_eta-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Nao,3.0000000000e+00)*__OLD_[1])));
					_prvt_calc_E_CaN = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Cao/__OLD_[1])));
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((((9.0000000000e-01*_prvt_Nao)+(1.0000000000e-01*_prvt_Ko))/((9.0000000000e-01*__OLD_[18])+(1.0000000000e-01*__OLD_[27])))));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ko/__OLD_[27])));
					_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[39]*(__OLD_[0]-(((_prvt_R*_prvt_T)/_prvt_F)*log((((9.8000000000e-01*_prvt_Ko)+(2.0000000000e-02*_prvt_Nao))/((9.8000000000e-01*__OLD_[27])+(2.0000000000e-02*__OLD_[18])))))));
					_prvt_calc_sigma = ((1.0000000000e+00/7.0000000000e+00)*(exp((_prvt_Nao/6.7300000000e+04))-1.0000000000e+00));
					_prvt_calc_O_ClCa = (2.0000000000e-01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.6700000000e+01))/7.8000000000e+00))));
					_prvt_calc_i_Nab = (_prvt_g_Nab*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Kto_s = (_prvt_g_Kto_s*__OLD_[30]*__OLD_[31]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K1 = ((((2.9380000000e-01*_prvt_Ko)/(_prvt_Ko+2.1000000000e+02))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp((8.9600000000e-02*(__OLD_[0]-_prvt_calc_E_K)))));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[32],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Kur = (_prvt_g_Kur*__OLD_[33]*__OLD_[34]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Kss = (_prvt_g_Kss*__OLD_[35]*__OLD_[36]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Cab = (_prvt_g_Cab*(__OLD_[0]-_prvt_calc_E_CaN));
					_prvt_calc_i_Na = (_prvt_g_Na*__OLD_[21]*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Kto_f = (_prvt_g_Kto_f*pow(__OLD_[28],3.0000000000e+00)*__OLD_[29]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_f_NaK = (1.0000000000e+00/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.6500000000e-02*_prvt_calc_sigma*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_ClCa = (((_prvt_g_ClCa*_prvt_calc_O_ClCa*__OLD_[1])/(__OLD_[1]+_prvt_Km_Cl))*(__OLD_[0]-_prvt_E_Cl));
					_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_calc_f_NaK*1.0000000000e+00)/(1.0000000000e+00+pow((_prvt_Km_Nai/__OLD_[18]),1.5000000000e+00)))*_prvt_Ko)/(_prvt_Ko+_prvt_Km_Ko));
					__K2_[0]= (-(_prvt_calc_i_CaL+_prvt_calc_i_pCa+_prvt_calc_i_NaCa+_prvt_calc_i_Cab+_prvt_calc_i_Na+_prvt_calc_i_Nab+_prvt_calc_i_NaK+_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kur+_prvt_calc_i_Kss+_prvt_calc_i_Kr+_prvt_calc_i_ClCa+_prvt_calc_i_stim));
					_prvt_aux_tol = fabs(__OLD_[0])*_prvt_rel_tol_;
					__TOL_[0] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[0] = fabs((_prvt_dtime/2) * (__K1_[0] - __K2_[0])/__TOL_[0]);
					__K2_[1]= (_prvt_calc_Bi*((_prvt_calc_J_leak+_prvt_calc_J_xfer)-(_prvt_calc_J_up+_prvt_calc_J_trpn+((((_prvt_calc_i_Cab+_prvt_calc_i_pCa)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vmyo*_prvt_F)))));
					_prvt_aux_tol = fabs(__OLD_[1])*_prvt_rel_tol_;
					__TOL_[1] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[1] = fabs((_prvt_dtime/2) * (__K1_[1] - __K2_[1])/__TOL_[1]);
					__K2_[2]= (_prvt_calc_Bss*(((_prvt_calc_J_rel*_prvt_VJSR)/_prvt_Vss)-(((_prvt_calc_J_xfer*_prvt_Vmyo)/_prvt_Vss)+((_prvt_calc_i_CaL*_prvt_Acap*_prvt_Cm)/(2.0000000000e+00*_prvt_Vss*_prvt_F)))));
					_prvt_aux_tol = fabs(__OLD_[2])*_prvt_rel_tol_;
					__TOL_[2] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[2] = fabs((_prvt_dtime/2) * (__K1_[2] - __K2_[2])/__TOL_[2]);
					__K2_[3]= (_prvt_calc_BJSR*(_prvt_calc_J_tr-_prvt_calc_J_rel));
					_prvt_aux_tol = fabs(__OLD_[3])*_prvt_rel_tol_;
					__TOL_[3] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[3] = fabs((_prvt_dtime/2) * (__K1_[3] - __K2_[3])/__TOL_[3]);
					__K2_[4]= ((((_prvt_calc_J_up-_prvt_calc_J_leak)*_prvt_Vmyo)/_prvt_VNSR)-((_prvt_calc_J_tr*_prvt_VJSR)/_prvt_VNSR));
					_prvt_aux_tol = fabs(__OLD_[4])*_prvt_rel_tol_;
					__TOL_[4] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[4] = fabs((_prvt_dtime/2) * (__K1_[4] - __K2_[4])/__TOL_[4]);
					__K2_[5]= (((-4.0000000000e-02)*__OLD_[5])-(((1.0000000000e-01*_prvt_calc_i_CaL)/_prvt_i_CaL_max)*exp(((-pow((__OLD_[0]-5.0000000000e+00),2.0000000000e+00))/6.4800000000e+02))));
					_prvt_aux_tol = fabs(__OLD_[5])*_prvt_rel_tol_;
					__TOL_[5] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[5] = fabs((_prvt_dtime/2) * (__K1_[5] - __K2_[5])/__TOL_[5]);
					__K2_[18]= (((-(_prvt_calc_i_Na+_prvt_calc_i_Nab+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
					_prvt_aux_tol = fabs(__OLD_[18])*_prvt_rel_tol_;
					__TOL_[18] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[18] = fabs((_prvt_dtime/2) * (__K1_[18] - __K2_[18])/__TOL_[18]);
					__K2_[27]= (((-((_prvt_calc_i_Kto_f+_prvt_calc_i_Kto_s+_prvt_calc_i_K1+_prvt_calc_i_Ks+_prvt_calc_i_Kss+_prvt_calc_i_Kur+_prvt_calc_i_Kr)-(2.0000000000e+00*_prvt_calc_i_NaK)))*_prvt_Acap*_prvt_Cm)/(_prvt_Vmyo*_prvt_F));
					_prvt_aux_tol = fabs(__OLD_[27])*_prvt_rel_tol_;
					__TOL_[27] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[27] = fabs((_prvt_dtime/2) * (__K1_[27] - __K2_[27])/__TOL_[27]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[1])
				{
					_prvt_calc_P_C1 = (1.0000000000e+00-(__OLD_[10]+__OLD_[8]+__OLD_[9]));
					__K2_[8]= (((_prvt_k_plus_a*pow(__OLD_[2],_prvt_n)*_prvt_calc_P_C1)+(_prvt_k_minus_b*__OLD_[9])+(_prvt_k_minus_c*__OLD_[10]))-((_prvt_k_minus_a*__OLD_[8])+(_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])+(_prvt_k_plus_c*__OLD_[8])));
					_prvt_aux_tol = fabs(__OLD_[8])*_prvt_rel_tol_;
					__TOL_[8] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[8] = fabs((_prvt_dtime/2) * (__K1_[8] - __K2_[8])/__TOL_[8]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[2])
				{
					_prvt_calc_C1 = (1.0000000000e+00-(__OLD_[11]+__OLD_[12]+__OLD_[12]+__OLD_[13]+__OLD_[14]+__OLD_[15]+__OLD_[16]+__OLD_[17]));
					_prvt_calc_alpha = ((4.0000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01))*((1.0000000000e+00+(7.0000000000e-01*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.0000000000e+01))))-(7.5000000000e-01*exp(((-pow((__OLD_[0]+2.0000000000e+01),2.0000000000e+00))/4.0000000000e+02)))))/(1.0000000000e+00+(1.2000000000e-01*exp(((__OLD_[0]+1.2000000000e+01)/1.0000000000e+01)))));
					_prvt_calc_beta = (5.0000000000e-02*exp(((-(__OLD_[0]+1.2000000000e+01))/1.3000000000e+01)));
					_prvt_calc_gamma = ((_prvt_Kpc_max*__OLD_[2])/(_prvt_Kpc_half+__OLD_[2]));
					_prvt_calc_Kpcf = (1.3000000000e+01*(1.0000000000e+00-exp(((-pow((__OLD_[0]+1.4500000000e+01),2.0000000000e+00))/1.0000000000e+02))));
					__K2_[11]= (((_prvt_calc_alpha*__OLD_[14])+(_prvt_Kpcb*__OLD_[15])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[11]))))-((4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(_prvt_calc_gamma*__OLD_[11])));
					_prvt_aux_tol = fabs(__OLD_[11])*_prvt_rel_tol_;
					__TOL_[11] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[11] = fabs((_prvt_dtime/2) * (__K1_[11] - __K2_[11])/__TOL_[11]);
					__K2_[12]= (((4.0000000000e+00*_prvt_calc_alpha*_prvt_calc_C1)+(2.0000000000e+00*_prvt_calc_beta*__OLD_[13]))-((_prvt_calc_beta*__OLD_[12])+(3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])));
					_prvt_aux_tol = fabs(__OLD_[12])*_prvt_rel_tol_;
					__TOL_[12] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[12] = fabs((_prvt_dtime/2) * (__K1_[12] - __K2_[12])/__TOL_[12]);
					__K2_[13]= (((3.0000000000e+00*_prvt_calc_alpha*__OLD_[12])+(3.0000000000e+00*_prvt_calc_beta*__OLD_[14]))-((2.0000000000e+00*_prvt_calc_beta*__OLD_[13])+(2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])));
					_prvt_aux_tol = fabs(__OLD_[13])*_prvt_rel_tol_;
					__TOL_[13] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[13] = fabs((_prvt_dtime/2) * (__K1_[13] - __K2_[13])/__TOL_[13]);
					__K2_[14]= (((2.0000000000e+00*_prvt_calc_alpha*__OLD_[13])+(4.0000000000e+00*_prvt_calc_beta*__OLD_[11])+(1.0000000000e-02*((4.0000000000e+00*_prvt_Kpcb*_prvt_calc_beta*__OLD_[15])-(_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])))+(2.0000000000e-03*((4.0000000000e+00*_prvt_calc_beta*__OLD_[16])-(_prvt_calc_Kpcf*__OLD_[14])))+(4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17]))-((3.0000000000e+00*_prvt_calc_beta*__OLD_[14])+(_prvt_calc_alpha*__OLD_[14])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14])));
					_prvt_aux_tol = fabs(__OLD_[14])*_prvt_rel_tol_;
					__TOL_[14] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[14] = fabs((_prvt_dtime/2) * (__K1_[14] - __K2_[14])/__TOL_[14]);
					__K2_[15]= (((_prvt_calc_gamma*__OLD_[11])+(1.0000000000e-03*((_prvt_calc_alpha*__OLD_[17])-(_prvt_calc_Kpcf*__OLD_[15])))+(1.0000000000e-02*((_prvt_calc_alpha*_prvt_calc_gamma*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*_prvt_calc_Kpcf*__OLD_[15]))))-(_prvt_Kpcb*__OLD_[15]));
					_prvt_aux_tol = fabs(__OLD_[15])*_prvt_rel_tol_;
					__TOL_[15] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[15] = fabs((_prvt_dtime/2) * (__K1_[15] - __K2_[15])/__TOL_[15]);
					__K2_[16]= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[11])-(_prvt_calc_alpha*__OLD_[16])))+(_prvt_Kpcb*__OLD_[17])+(2.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[14])-(4.0000000000e+00*_prvt_calc_beta*__OLD_[16]))))-(_prvt_calc_gamma*__OLD_[16]));
					_prvt_aux_tol = fabs(__OLD_[16])*_prvt_rel_tol_;
					__TOL_[16] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[16] = fabs((_prvt_dtime/2) * (__K1_[16] - __K2_[16])/__TOL_[16]);
					__K2_[17]= (((1.0000000000e-03*((_prvt_calc_Kpcf*__OLD_[15])-(_prvt_calc_alpha*__OLD_[17])))+(_prvt_calc_gamma*__OLD_[16])+(1.0000000000e+00*_prvt_calc_gamma*_prvt_calc_Kpcf*__OLD_[14]))-((4.0000000000e+00*_prvt_calc_beta*_prvt_Kpcb*__OLD_[17])+(_prvt_Kpcb*__OLD_[17])));
					_prvt_aux_tol = fabs(__OLD_[17])*_prvt_rel_tol_;
					__TOL_[17] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[17] = fabs((_prvt_dtime/2) * (__K1_[17] - __K2_[17])/__TOL_[17]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[3])
				{
					_prvt_calc_C_Na3 = (1.0000000000e+00-(__OLD_[21]+__OLD_[20]+__OLD_[19]+__OLD_[22]+__OLD_[23]+__OLD_[24]+__OLD_[25]+__OLD_[26]));
					_prvt_calc_alpha_Na11 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.7000000000e+01)))+(2.0000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_alpha_Na12 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+01)))+(2.3000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_alpha_Na13 = (3.8020000000e+00/((1.0270000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.2000000000e+01)))+(2.5000000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/1.5000000000e+02)))));
					_prvt_calc_beta_Na11 = (1.9170000000e-01*exp(((-(__OLD_[0]+2.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_beta_Na12 = (2.0000000000e-01*exp(((-(__OLD_[0]-2.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_beta_Na13 = (2.2000000000e-01*exp(((-(__OLD_[0]-7.5000000000e+00))/2.0300000000e+01)));
					_prvt_calc_alpha_Na3 = (7.0000000000e-07*exp(((-(__OLD_[0]+7.0000000000e+00))/7.7000000000e+00)));
					_prvt_calc_beta_Na3 = (8.4000000000e-03+(2.0000000000e-05*(__OLD_[0]+7.0000000000e+00)));
					_prvt_calc_alpha_Na2 = (1.0000000000e+00/((1.8849500000e-01*exp(((-(__OLD_[0]+7.0000000000e+00))/1.6600000000e+01)))+3.9395600000e-01));
					_prvt_calc_beta_Na2 = ((_prvt_calc_alpha_Na13*_prvt_calc_alpha_Na2*_prvt_calc_alpha_Na3)/(_prvt_calc_beta_Na13*_prvt_calc_beta_Na3));
					_prvt_calc_alpha_Na4 = (_prvt_calc_alpha_Na2/1.0000000000e+03);
					_prvt_calc_beta_Na4 = _prvt_calc_alpha_Na3;
					_prvt_calc_alpha_Na5 = (_prvt_calc_alpha_Na2/9.5000000000e+04);
					_prvt_calc_beta_Na5 = (_prvt_calc_alpha_Na3/5.0000000000e+01);
					__K2_[19]= (((_prvt_calc_alpha_Na11*_prvt_calc_C_Na3)+(_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[19])+(_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na3*__OLD_[19])));
					_prvt_aux_tol = fabs(__OLD_[19])*_prvt_rel_tol_;
					__TOL_[19] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[19] = fabs((_prvt_dtime/2) * (__K1_[19] - __K2_[19])/__TOL_[19]);
					__K2_[20]= (((_prvt_calc_alpha_Na12*__OLD_[19])+(_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na3*__OLD_[22]))-((_prvt_calc_beta_Na12*__OLD_[20])+(_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na3*__OLD_[20])));
					_prvt_aux_tol = fabs(__OLD_[20])*_prvt_rel_tol_;
					__TOL_[20] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[20] = fabs((_prvt_dtime/2) * (__K1_[20] - __K2_[20])/__TOL_[20]);
					__K2_[21]= (((_prvt_calc_alpha_Na13*__OLD_[20])+(_prvt_calc_beta_Na2*__OLD_[22]))-((_prvt_calc_beta_Na13*__OLD_[21])+(_prvt_calc_alpha_Na2*__OLD_[21])));
					_prvt_aux_tol = fabs(__OLD_[21])*_prvt_rel_tol_;
					__TOL_[21] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[21] = fabs((_prvt_dtime/2) * (__K1_[21] - __K2_[21])/__TOL_[21]);
					__K2_[22]= (((_prvt_calc_alpha_Na2*__OLD_[21])+(_prvt_calc_beta_Na3*__OLD_[20])+(_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na12*__OLD_[25]))-((_prvt_calc_beta_Na2*__OLD_[22])+(_prvt_calc_alpha_Na3*__OLD_[22])+(_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na12*__OLD_[22])));
					_prvt_aux_tol = fabs(__OLD_[22])*_prvt_rel_tol_;
					__TOL_[22] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[22] = fabs((_prvt_dtime/2) * (__K1_[22] - __K2_[22])/__TOL_[22]);
					__K2_[23]= (((_prvt_calc_alpha_Na4*__OLD_[22])+(_prvt_calc_beta_Na5*__OLD_[24]))-((_prvt_calc_beta_Na4*__OLD_[23])+(_prvt_calc_alpha_Na5*__OLD_[23])));
					_prvt_aux_tol = fabs(__OLD_[23])*_prvt_rel_tol_;
					__TOL_[23] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[23] = fabs((_prvt_dtime/2) * (__K1_[23] - __K2_[23])/__TOL_[23]);
					__K2_[24]= ((_prvt_calc_alpha_Na5*__OLD_[23])-(_prvt_calc_beta_Na5*__OLD_[24]));
					_prvt_aux_tol = fabs(__OLD_[24])*_prvt_rel_tol_;
					__TOL_[24] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[24] = fabs((_prvt_dtime/2) * (__K1_[24] - __K2_[24])/__TOL_[24]);
					__K2_[25]= (((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_beta_Na12*__OLD_[22])+(_prvt_calc_beta_Na3*__OLD_[25]))-((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_alpha_Na12*__OLD_[25])+(_prvt_calc_alpha_Na3*__OLD_[25])));
					_prvt_aux_tol = fabs(__OLD_[25])*_prvt_rel_tol_;
					__TOL_[25] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[25] = fabs((_prvt_dtime/2) * (__K1_[25] - __K2_[25])/__TOL_[25]);
					__K2_[26]= (((_prvt_calc_beta_Na11*__OLD_[25])+(_prvt_calc_beta_Na3*_prvt_calc_C_Na3))-((_prvt_calc_alpha_Na11*__OLD_[26])+(_prvt_calc_alpha_Na3*__OLD_[26])));
					_prvt_aux_tol = fabs(__OLD_[26])*_prvt_rel_tol_;
					__TOL_[26] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[26] = fabs((_prvt_dtime/2) * (__K1_[26] - __K2_[26])/__TOL_[26]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[4])
				{
					_prvt_calc_alpha_a = (1.8064000000e-01*exp((3.5770000000e-02*(__OLD_[0]+3.0000000000e+01))));
					_prvt_calc_beta_a = (3.9560000000e-01*exp(((-6.2370000000e-02)*(__OLD_[0]+3.0000000000e+01))));
					__K2_[28]= ((_prvt_calc_alpha_a*(1.0000000000e+00-__OLD_[28]))-(_prvt_calc_beta_a*__OLD_[28]));
					_prvt_aux_tol = fabs(__OLD_[28])*_prvt_rel_tol_;
					__TOL_[28] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[28] = fabs((_prvt_dtime/2) * (__K1_[28] - __K2_[28])/__TOL_[28]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[5])
				{
					_prvt_calc_alpha_i = ((1.5200000000e-04*exp(((-(__OLD_[0]+1.3500000000e+01))/7.0000000000e+00)))/((6.7083000000e-03*exp(((-(__OLD_[0]+3.3500000000e+01))/7.0000000000e+00)))+1.0000000000e+00));
					_prvt_calc_beta_i = ((9.5000000000e-04*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))/((5.1335000000e-02*exp(((__OLD_[0]+3.3500000000e+01)/7.0000000000e+00)))+1.0000000000e+00));
					__K2_[29]= ((_prvt_calc_alpha_i*(1.0000000000e+00-__OLD_[29]))-(_prvt_calc_beta_i*__OLD_[29]));
					_prvt_aux_tol = fabs(__OLD_[29])*_prvt_rel_tol_;
					__TOL_[29] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[29] = fabs((_prvt_dtime/2) * (__K1_[29] - __K2_[29])/__TOL_[29]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[6])
				{
					_prvt_calc_ass = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2500000000e+01))/7.7000000000e+00))));
					_prvt_calc_tau_ta_s = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
					_prvt_calc_tau_aur = ((4.9300000000e-01*exp(((-6.2900000000e-02)*__OLD_[0])))+2.0580000000e+00);
					_prvt_calc_tau_Kss = ((3.9300000000e+01*exp(((-8.6200000000e-02)*__OLD_[0])))+1.3170000000e+01);
					__K2_[30]= ((_prvt_calc_ass-__OLD_[30])/_prvt_calc_tau_ta_s);
					_prvt_aux_tol = fabs(__OLD_[30])*_prvt_rel_tol_;
					__TOL_[30] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[30] = fabs((_prvt_dtime/2) * (__K1_[30] - __K2_[30])/__TOL_[30]);
					__K2_[33]= ((_prvt_calc_ass-__OLD_[33])/_prvt_calc_tau_aur);
					_prvt_aux_tol = fabs(__OLD_[33])*_prvt_rel_tol_;
					__TOL_[33] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[33] = fabs((_prvt_dtime/2) * (__K1_[33] - __K2_[33])/__TOL_[33]);
					__K2_[35]= ((_prvt_calc_ass-__OLD_[35])/_prvt_calc_tau_Kss);
					_prvt_aux_tol = fabs(__OLD_[35])*_prvt_rel_tol_;
					__TOL_[35] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[35] = fabs((_prvt_dtime/2) * (__K1_[35] - __K2_[35])/__TOL_[35]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[7])
				{
					_prvt_calc_iss = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00))));
					_prvt_calc_tau_ti_s = (2.7000000000e+02+(1.0500000000e+03/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
					_prvt_calc_tau_iur = (1.2000000000e+03-(1.7000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+4.5200000000e+01)/5.7000000000e+00)))));
					__K2_[31]= ((_prvt_calc_iss-__OLD_[31])/_prvt_calc_tau_ti_s);
					_prvt_aux_tol = fabs(__OLD_[31])*_prvt_rel_tol_;
					__TOL_[31] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[31] = fabs((_prvt_dtime/2) * (__K1_[31] - __K2_[31])/__TOL_[31]);
					__K2_[34]= ((_prvt_calc_iss-__OLD_[34])/_prvt_calc_tau_iur);
					_prvt_aux_tol = fabs(__OLD_[34])*_prvt_rel_tol_;
					__TOL_[34] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[34] = fabs((_prvt_dtime/2) * (__K1_[34] - __K2_[34])/__TOL_[34]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[8])
				{
					_prvt_calc_alpha_n = ((4.8133300000e-06*(__OLD_[0]+2.6500000000e+01))/(1.0000000000e+00-exp(((-1.2800000000e-01)*(__OLD_[0]+2.6500000000e+01)))));
					_prvt_calc_beta_n = (9.5333300000e-05*exp(((-3.8000000000e-02)*(__OLD_[0]+2.6500000000e+01))));
					__K2_[32]= ((_prvt_calc_alpha_n*(1.0000000000e+00-__OLD_[32]))-(_prvt_calc_beta_n*__OLD_[32]));
					_prvt_aux_tol = fabs(__OLD_[32])*_prvt_rel_tol_;
					__TOL_[32] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[32] = fabs((_prvt_dtime/2) * (__K1_[32] - __K2_[32])/__TOL_[32]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[9])
				{
					_prvt_calc_C_K0 = (1.0000000000e+00-(__OLD_[38]+__OLD_[37]+__OLD_[39]+__OLD_[40]));
					_prvt_calc_alpha_a0 = (2.2348000000e-02*exp((1.1760000000e-02*__OLD_[0])));
					_prvt_calc_beta_a0 = (4.7002000000e-02*exp(((-6.3100000000e-02)*__OLD_[0])));
					__K2_[38]= (((_prvt_calc_alpha_a0*_prvt_calc_C_K0)+(_prvt_kb*__OLD_[37]))-((_prvt_calc_beta_a0*__OLD_[38])+(_prvt_kf*__OLD_[38])));
					_prvt_aux_tol = fabs(__OLD_[38])*_prvt_rel_tol_;
					__TOL_[38] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[38] = fabs((_prvt_dtime/2) * (__K1_[38] - __K2_[38])/__TOL_[38]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[10])
				{
					_prvt_calc_alpha_a1 = (1.3733000000e-02*exp((3.8198000000e-02*__OLD_[0])));
					_prvt_calc_beta_a1 = (6.8900000000e-05*exp(((-4.1780000000e-02)*__OLD_[0])));
					_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current = (9.0821000000e-02*exp((2.3391000000e-02*(__OLD_[0]+5.0000000000e+00))));
					_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current = (6.4970000000e-03*exp(((-3.2680000000e-02)*(__OLD_[0]+5.0000000000e+00))));
					__K2_[37]= (((_prvt_kf*__OLD_[38])+(_prvt_calc_beta_a1*__OLD_[39]))-((_prvt_kb*__OLD_[37])+(_prvt_calc_alpha_a1*__OLD_[37])));
					_prvt_aux_tol = fabs(__OLD_[37])*_prvt_rel_tol_;
					__TOL_[37] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[37] = fabs((_prvt_dtime/2) * (__K1_[37] - __K2_[37])/__TOL_[37]);
					__K2_[39]= (((_prvt_calc_alpha_a1*__OLD_[37])+(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]))-((_prvt_calc_beta_a1*__OLD_[39])+(_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])));
					_prvt_aux_tol = fabs(__OLD_[39])*_prvt_rel_tol_;
					__TOL_[39] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[39] = fabs((_prvt_dtime/2) * (__K1_[39] - __K2_[39])/__TOL_[39]);
					__K2_[40]= ((_prvt_calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[39])-(_prvt_calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current*__OLD_[40]));
					_prvt_aux_tol = fabs(__OLD_[40])*_prvt_rel_tol_;
					__TOL_[40] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[40] = fabs((_prvt_dtime/2) * (__K1_[40] - __K2_[40])/__TOL_[40]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[11])
				{
					__K2_[6]= ((_prvt_k_plus_ltrpn*__OLD_[1]*(_prvt_LTRPN_tot-__OLD_[6]))-(_prvt_k_minus_ltrpn*__OLD_[6]));
					_prvt_aux_tol = fabs(__OLD_[6])*_prvt_rel_tol_;
					__TOL_[6] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[6] = fabs((_prvt_dtime/2) * (__K1_[6] - __K2_[6])/__TOL_[6]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[12])
				{
					__K2_[7]= ((_prvt_k_plus_htrpn*__OLD_[1]*(_prvt_HTRPN_tot-__OLD_[7]))-(_prvt_k_minus_htrpn*__OLD_[7]));
					_prvt_aux_tol = fabs(__OLD_[7])*_prvt_rel_tol_;
					__TOL_[7] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[7] = fabs((_prvt_dtime/2) * (__K1_[7] - __K2_[7])/__TOL_[7]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[13])
				{
					__K2_[9]= ((_prvt_k_plus_b*pow(__OLD_[2],_prvt_m)*__OLD_[8])-(_prvt_k_minus_b*__OLD_[9]));
					_prvt_aux_tol = fabs(__OLD_[9])*_prvt_rel_tol_;
					__TOL_[9] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[9] = fabs((_prvt_dtime/2) * (__K1_[9] - __K2_[9])/__TOL_[9]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[14])
				{
					__K2_[10]= ((_prvt_k_plus_c*__OLD_[8])-(_prvt_k_minus_c*__OLD_[10]));
					_prvt_aux_tol = fabs(__OLD_[10])*_prvt_rel_tol_;
					__TOL_[10] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[10] = fabs((_prvt_dtime/2) * (__K1_[10] - __K2_[10])/__TOL_[10]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[15])
				{
					__K2_[36]= 0.0000000000e+00;
					_prvt_aux_tol = fabs(__OLD_[36])*_prvt_rel_tol_;
					__TOL_[36] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[36] = fabs((_prvt_dtime/2) * (__K1_[36] - __K2_[36])/__TOL_[36]);
				}
				_prvt_time_new -= _prvt_dtime;
				#pragma omp barrier
				///adapt the time step
				double __greatestError_=0.0;
				for(int k=0;k<numEDO;k++){
						if(__ERROR_[k] > __greatestError_)
							__greatestError_ = __ERROR_[k];
				}
				__greatestError_ += __tiny_;
				_prvt_previous_dt = _prvt_dtime;
				_prvt_dtime = _beta_safety_ * _prvt_dtime * sqrt(1.0/__greatestError_);
				if(__greatestError_>=1){
					//throw the results away and compute again
					if(omp_get_thread_num()==_prvt_tree_thread[0])
					{
						__NEW_[0] = __K1_[0] * _prvt_dtime + __OLD_AUX_[0];
						__NEW_[1] = __K1_[1] * _prvt_dtime + __OLD_AUX_[1];
						__NEW_[2] = __K1_[2] * _prvt_dtime + __OLD_AUX_[2];
						__NEW_[3] = __K1_[3] * _prvt_dtime + __OLD_AUX_[3];
						__NEW_[4] = __K1_[4] * _prvt_dtime + __OLD_AUX_[4];
						__NEW_[5] = __K1_[5] * _prvt_dtime + __OLD_AUX_[5];
						__NEW_[18] = __K1_[18] * _prvt_dtime + __OLD_AUX_[18];
						__NEW_[27] = __K1_[27] * _prvt_dtime + __OLD_AUX_[27];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[8] = __K1_[8] * _prvt_dtime + __OLD_AUX_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[11] = __K1_[11] * _prvt_dtime + __OLD_AUX_[11];
						__NEW_[12] = __K1_[12] * _prvt_dtime + __OLD_AUX_[12];
						__NEW_[13] = __K1_[13] * _prvt_dtime + __OLD_AUX_[13];
						__NEW_[14] = __K1_[14] * _prvt_dtime + __OLD_AUX_[14];
						__NEW_[15] = __K1_[15] * _prvt_dtime + __OLD_AUX_[15];
						__NEW_[16] = __K1_[16] * _prvt_dtime + __OLD_AUX_[16];
						__NEW_[17] = __K1_[17] * _prvt_dtime + __OLD_AUX_[17];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[19] = __K1_[19] * _prvt_dtime + __OLD_AUX_[19];
						__NEW_[20] = __K1_[20] * _prvt_dtime + __OLD_AUX_[20];
						__NEW_[21] = __K1_[21] * _prvt_dtime + __OLD_AUX_[21];
						__NEW_[22] = __K1_[22] * _prvt_dtime + __OLD_AUX_[22];
						__NEW_[23] = __K1_[23] * _prvt_dtime + __OLD_AUX_[23];
						__NEW_[24] = __K1_[24] * _prvt_dtime + __OLD_AUX_[24];
						__NEW_[25] = __K1_[25] * _prvt_dtime + __OLD_AUX_[25];
						__NEW_[26] = __K1_[26] * _prvt_dtime + __OLD_AUX_[26];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[28] = __K1_[28] * _prvt_dtime + __OLD_AUX_[28];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[29] = __K1_[29] * _prvt_dtime + __OLD_AUX_[29];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[30] = __K1_[30] * _prvt_dtime + __OLD_AUX_[30];
						__NEW_[33] = __K1_[33] * _prvt_dtime + __OLD_AUX_[33];
						__NEW_[35] = __K1_[35] * _prvt_dtime + __OLD_AUX_[35];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[31] = __K1_[31] * _prvt_dtime + __OLD_AUX_[31];
						__NEW_[34] = __K1_[34] * _prvt_dtime + __OLD_AUX_[34];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[32] = __K1_[32] * _prvt_dtime + __OLD_AUX_[32];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[38] = __K1_[38] * _prvt_dtime + __OLD_AUX_[38];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[37] = __K1_[37] * _prvt_dtime + __OLD_AUX_[37];
						__NEW_[39] = __K1_[39] * _prvt_dtime + __OLD_AUX_[39];
						__NEW_[40] = __K1_[40] * _prvt_dtime + __OLD_AUX_[40];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[6] = __K1_[6] * _prvt_dtime + __OLD_AUX_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[7] = __K1_[7] * _prvt_dtime + __OLD_AUX_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[13])
					{
						__NEW_[9] = __K1_[9] * _prvt_dtime + __OLD_AUX_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[14])
					{
						__NEW_[10] = __K1_[10] * _prvt_dtime + __OLD_AUX_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[15])
					{
						__NEW_[36] = __K1_[36] * _prvt_dtime + __OLD_AUX_[36];
					}
				__TEMP_ = __NEW_;
				__NEW_ = __OLD_;
				__OLD_ = __TEMP_;
				}else{//it accepts the solutions
					if(_prvt_savingRate!=0){
						#pragma omp single
						{
							this->V_old_ = __OLD_AUX_[0];
							this->V_new_ = __OLD_[0];
							this->Cai_old_ = __OLD_AUX_[1];
							this->Cai_new_ = __OLD_[1];
							this->Cass_old_ = __OLD_AUX_[2];
							this->Cass_new_ = __OLD_[2];
							this->CaJSR_old_ = __OLD_AUX_[3];
							this->CaJSR_new_ = __OLD_[3];
							this->CaNSR_old_ = __OLD_AUX_[4];
							this->CaNSR_new_ = __OLD_[4];
							this->P_RyR_old_ = __OLD_AUX_[5];
							this->P_RyR_new_ = __OLD_[5];
							this->LTRPN_Ca_old_ = __OLD_AUX_[6];
							this->LTRPN_Ca_new_ = __OLD_[6];
							this->HTRPN_Ca_old_ = __OLD_AUX_[7];
							this->HTRPN_Ca_new_ = __OLD_[7];
							this->P_O1_old_ = __OLD_AUX_[8];
							this->P_O1_new_ = __OLD_[8];
							this->P_O2_old_ = __OLD_AUX_[9];
							this->P_O2_new_ = __OLD_[9];
							this->P_C2_old_ = __OLD_AUX_[10];
							this->P_C2_new_ = __OLD_[10];
							this->O_old_ = __OLD_AUX_[11];
							this->O_new_ = __OLD_[11];
							this->C2_old_ = __OLD_AUX_[12];
							this->C2_new_ = __OLD_[12];
							this->C3_old_ = __OLD_AUX_[13];
							this->C3_new_ = __OLD_[13];
							this->C4_old_ = __OLD_AUX_[14];
							this->C4_new_ = __OLD_[14];
							this->I1_old_ = __OLD_AUX_[15];
							this->I1_new_ = __OLD_[15];
							this->I2_old_ = __OLD_AUX_[16];
							this->I2_new_ = __OLD_[16];
							this->I3_old_ = __OLD_AUX_[17];
							this->I3_new_ = __OLD_[17];
							this->Nai_old_ = __OLD_AUX_[18];
							this->Nai_new_ = __OLD_[18];
							this->C_Na2_old_ = __OLD_AUX_[19];
							this->C_Na2_new_ = __OLD_[19];
							this->C_Na1_old_ = __OLD_AUX_[20];
							this->C_Na1_new_ = __OLD_[20];
							this->O_Na_old_ = __OLD_AUX_[21];
							this->O_Na_new_ = __OLD_[21];
							this->IF_Na_old_ = __OLD_AUX_[22];
							this->IF_Na_new_ = __OLD_[22];
							this->I1_Na_old_ = __OLD_AUX_[23];
							this->I1_Na_new_ = __OLD_[23];
							this->I2_Na_old_ = __OLD_AUX_[24];
							this->I2_Na_new_ = __OLD_[24];
							this->IC_Na2_old_ = __OLD_AUX_[25];
							this->IC_Na2_new_ = __OLD_[25];
							this->IC_Na3_old_ = __OLD_AUX_[26];
							this->IC_Na3_new_ = __OLD_[26];
							this->Ki_old_ = __OLD_AUX_[27];
							this->Ki_new_ = __OLD_[27];
							this->ato_f_old_ = __OLD_AUX_[28];
							this->ato_f_new_ = __OLD_[28];
							this->ito_f_old_ = __OLD_AUX_[29];
							this->ito_f_new_ = __OLD_[29];
							this->ato_s_old_ = __OLD_AUX_[30];
							this->ato_s_new_ = __OLD_[30];
							this->ito_s_old_ = __OLD_AUX_[31];
							this->ito_s_new_ = __OLD_[31];
							this->nKs_old_ = __OLD_AUX_[32];
							this->nKs_new_ = __OLD_[32];
							this->aur_old_ = __OLD_AUX_[33];
							this->aur_new_ = __OLD_[33];
							this->iur_old_ = __OLD_AUX_[34];
							this->iur_new_ = __OLD_[34];
							this->aKss_old_ = __OLD_AUX_[35];
							this->aKss_new_ = __OLD_[35];
							this->iKss_old_ = __OLD_AUX_[36];
							this->iKss_new_ = __OLD_[36];
							this->C_K2_old_ = __OLD_AUX_[37];
							this->C_K2_new_ = __OLD_[37];
							this->C_K1_old_ = __OLD_AUX_[38];
							this->C_K1_new_ = __OLD_[38];
							this->O_K_old_ = __OLD_AUX_[39];
							this->O_K_new_ = __OLD_[39];
							this->I_K_old_ = __OLD_AUX_[40];
							this->I_K_new_ = __OLD_[40];
							this->previous_dt = _prvt_previous_dt;
							this->dtime = _prvt_dtime;
							this->time_new = _prvt_time_new;
							save_step(fileptr, _ADAP_DT_);
						}
					}
					if(_prvt_dtime > _prvt_maxStep && _prvt_maxStep!=0){
						_prvt_dtime = _prvt_maxStep;
					}else if(_prvt_dtime==0){
						printf("Error: Time step is zero.\n");
						break;
					}
					//it steps the method ahead, with euler solution
					if(omp_get_thread_num()==_prvt_tree_thread[0])
					{
						__NEW_[0] = __K2_[0] * _prvt_dtime + __OLD_[0];
						__NEW_[1] = __K2_[1] * _prvt_dtime + __OLD_[1];
						__NEW_[2] = __K2_[2] * _prvt_dtime + __OLD_[2];
						__NEW_[3] = __K2_[3] * _prvt_dtime + __OLD_[3];
						__NEW_[4] = __K2_[4] * _prvt_dtime + __OLD_[4];
						__NEW_[5] = __K2_[5] * _prvt_dtime + __OLD_[5];
						__NEW_[18] = __K2_[18] * _prvt_dtime + __OLD_[18];
						__NEW_[27] = __K2_[27] * _prvt_dtime + __OLD_[27];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[8] = __K2_[8] * _prvt_dtime + __OLD_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[11] = __K2_[11] * _prvt_dtime + __OLD_[11];
						__NEW_[12] = __K2_[12] * _prvt_dtime + __OLD_[12];
						__NEW_[13] = __K2_[13] * _prvt_dtime + __OLD_[13];
						__NEW_[14] = __K2_[14] * _prvt_dtime + __OLD_[14];
						__NEW_[15] = __K2_[15] * _prvt_dtime + __OLD_[15];
						__NEW_[16] = __K2_[16] * _prvt_dtime + __OLD_[16];
						__NEW_[17] = __K2_[17] * _prvt_dtime + __OLD_[17];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[19] = __K2_[19] * _prvt_dtime + __OLD_[19];
						__NEW_[20] = __K2_[20] * _prvt_dtime + __OLD_[20];
						__NEW_[21] = __K2_[21] * _prvt_dtime + __OLD_[21];
						__NEW_[22] = __K2_[22] * _prvt_dtime + __OLD_[22];
						__NEW_[23] = __K2_[23] * _prvt_dtime + __OLD_[23];
						__NEW_[24] = __K2_[24] * _prvt_dtime + __OLD_[24];
						__NEW_[25] = __K2_[25] * _prvt_dtime + __OLD_[25];
						__NEW_[26] = __K2_[26] * _prvt_dtime + __OLD_[26];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[28] = __K2_[28] * _prvt_dtime + __OLD_[28];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[29] = __K2_[29] * _prvt_dtime + __OLD_[29];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[30] = __K2_[30] * _prvt_dtime + __OLD_[30];
						__NEW_[33] = __K2_[33] * _prvt_dtime + __OLD_[33];
						__NEW_[35] = __K2_[35] * _prvt_dtime + __OLD_[35];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[31] = __K2_[31] * _prvt_dtime + __OLD_[31];
						__NEW_[34] = __K2_[34] * _prvt_dtime + __OLD_[34];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[32] = __K2_[32] * _prvt_dtime + __OLD_[32];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[38] = __K2_[38] * _prvt_dtime + __OLD_[38];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[37] = __K2_[37] * _prvt_dtime + __OLD_[37];
						__NEW_[39] = __K2_[39] * _prvt_dtime + __OLD_[39];
						__NEW_[40] = __K2_[40] * _prvt_dtime + __OLD_[40];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[6] = __K2_[6] * _prvt_dtime + __OLD_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[7] = __K2_[7] * _prvt_dtime + __OLD_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[13])
					{
						__NEW_[9] = __K2_[9] * _prvt_dtime + __OLD_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[14])
					{
						__NEW_[10] = __K2_[10] * _prvt_dtime + __OLD_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[15])
					{
						__NEW_[36] = __K2_[36] * _prvt_dtime + __OLD_[36];
					}
					//store the old iteration in a aux 
					__TEMP_ = __OLD_;
					__OLD_ = __OLD_AUX_;
					__OLD_AUX_ = __TEMP_;
					//steps ahead with euler
					__TEMP_ = __NEW_;
					__NEW_ = __OLD_;
					__OLD_ = __TEMP_;
					//change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração
					__TEMP_	= __K2_;
					__K2_	= __K1_;
					__K1_	= __TEMP_;
					//sums the old dtime - the variable dtime is alreaady updated
					_prvt_time_new += _prvt_previous_dt;
				}//FIM ELSE
				#pragma omp barrier
			}
		}
	}
	int Solveode::getErrorCode(double error,double tolerance){
		if (error <0.5*tolerance){
			//Accept current solution, and increase the size of the next mechanics step.
			return 4;
		}if (error>=0.5*tolerance && error<tolerance){
			//Accept current solution, and hold the step size fixed for the next mechanics step.
			return 0;
		}if (error>=tolerance && error<2*tolerance){
			//Accept current solution, but decrease the size of the next mechanics step.
			return 3;
		}else{
			//Throw current results away, cut the step size in half, and retry.
			return -1;
		}
	}

	//Get Methods

	double Solveode::getVariables(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return V_old_;    break;
		case 1:		return Cai_old_;    break;
		case 2:		return Cass_old_;    break;
		case 3:		return CaJSR_old_;    break;
		case 4:		return CaNSR_old_;    break;
		case 5:		return P_RyR_old_;    break;
		case 6:		return LTRPN_Ca_old_;    break;
		case 7:		return HTRPN_Ca_old_;    break;
		case 8:		return P_O1_old_;    break;
		case 9:		return P_O2_old_;    break;
		case 10:		return P_C2_old_;    break;
		case 11:		return O_old_;    break;
		case 12:		return C2_old_;    break;
		case 13:		return C3_old_;    break;
		case 14:		return C4_old_;    break;
		case 15:		return I1_old_;    break;
		case 16:		return I2_old_;    break;
		case 17:		return I3_old_;    break;
		case 18:		return Nai_old_;    break;
		case 19:		return C_Na2_old_;    break;
		case 20:		return C_Na1_old_;    break;
		case 21:		return O_Na_old_;    break;
		case 22:		return IF_Na_old_;    break;
		case 23:		return I1_Na_old_;    break;
		case 24:		return I2_Na_old_;    break;
		case 25:		return IC_Na2_old_;    break;
		case 26:		return IC_Na3_old_;    break;
		case 27:		return Ki_old_;    break;
		case 28:		return ato_f_old_;    break;
		case 29:		return ito_f_old_;    break;
		case 30:		return ato_s_old_;    break;
		case 31:		return ito_s_old_;    break;
		case 32:		return nKs_old_;    break;
		case 33:		return aur_old_;    break;
		case 34:		return iur_old_;    break;
		case 35:		return aKss_old_;    break;
		case 36:		return iKss_old_;    break;
		case 37:		return C_K2_old_;    break;
		case 38:		return C_K1_old_;    break;
		case 39:		return O_K_old_;    break;
		case 40:		return I_K_old_;    break;
		default:	return 1;    break;
		}
	}

	double Solveode::getLadoDireito(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return V_lado_direito_;    break;
		case 1:		return Cai_lado_direito_;    break;
		case 2:		return Cass_lado_direito_;    break;
		case 3:		return CaJSR_lado_direito_;    break;
		case 4:		return CaNSR_lado_direito_;    break;
		case 5:		return P_RyR_lado_direito_;    break;
		case 6:		return LTRPN_Ca_lado_direito_;    break;
		case 7:		return HTRPN_Ca_lado_direito_;    break;
		case 8:		return P_O1_lado_direito_;    break;
		case 9:		return P_O2_lado_direito_;    break;
		case 10:		return P_C2_lado_direito_;    break;
		case 11:		return O_lado_direito_;    break;
		case 12:		return C2_lado_direito_;    break;
		case 13:		return C3_lado_direito_;    break;
		case 14:		return C4_lado_direito_;    break;
		case 15:		return I1_lado_direito_;    break;
		case 16:		return I2_lado_direito_;    break;
		case 17:		return I3_lado_direito_;    break;
		case 18:		return Nai_lado_direito_;    break;
		case 19:		return C_Na2_lado_direito_;    break;
		case 20:		return C_Na1_lado_direito_;    break;
		case 21:		return O_Na_lado_direito_;    break;
		case 22:		return IF_Na_lado_direito_;    break;
		case 23:		return I1_Na_lado_direito_;    break;
		case 24:		return I2_Na_lado_direito_;    break;
		case 25:		return IC_Na2_lado_direito_;    break;
		case 26:		return IC_Na3_lado_direito_;    break;
		case 27:		return Ki_lado_direito_;    break;
		case 28:		return ato_f_lado_direito_;    break;
		case 29:		return ito_f_lado_direito_;    break;
		case 30:		return ato_s_lado_direito_;    break;
		case 31:		return ito_s_lado_direito_;    break;
		case 32:		return nKs_lado_direito_;    break;
		case 33:		return aur_lado_direito_;    break;
		case 34:		return iur_lado_direito_;    break;
		case 35:		return aKss_lado_direito_;    break;
		case 36:		return iKss_lado_direito_;    break;
		case 37:		return C_K2_lado_direito_;    break;
		case 38:		return C_K1_lado_direito_;    break;
		case 39:		return O_K_lado_direito_;    break;
		case 40:		return I_K_lado_direito_;    break;
		default:	return 1;    break;
		}
	}

	double Solveode::getParameters(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return time;    break;
		case 1:		return stim_amplitude;    break;
		case 2:		return stim_start;    break;
		case 3:		return stim_end;    break;
		case 4:		return stim_period;    break;
		case 5:		return stim_duration;    break;
		case 6:		return Acap;    break;
		case 7:		return Cm;    break;
		case 8:		return Vmyo;    break;
		case 9:		return F;    break;
		case 10:		return VJSR;    break;
		case 11:		return Vss;    break;
		case 12:		return VNSR;    break;
		case 13:		return CMDN_tot;    break;
		case 14:		return Km_CMDN;    break;
		case 15:		return CSQN_tot;    break;
		case 16:		return Km_CSQN;    break;
		case 17:		return v1;    break;
		case 18:		return tau_tr;    break;
		case 19:		return tau_xfer;    break;
		case 20:		return v2;    break;
		case 21:		return v3;    break;
		case 22:		return Km_up;    break;
		case 23:		return k_plus_htrpn;    break;
		case 24:		return HTRPN_tot;    break;
		case 25:		return k_plus_ltrpn;    break;
		case 26:		return LTRPN_tot;    break;
		case 27:		return k_minus_htrpn;    break;
		case 28:		return k_minus_ltrpn;    break;
		case 29:		return i_CaL_max;    break;
		case 30:		return k_plus_a;    break;
		case 31:		return n;    break;
		case 32:		return k_minus_b;    break;
		case 33:		return k_minus_c;    break;
		case 34:		return k_minus_a;    break;
		case 35:		return k_plus_b;    break;
		case 36:		return m;    break;
		case 37:		return k_plus_c;    break;
		case 38:		return g_CaL;    break;
		case 39:		return E_CaL;    break;
		case 40:		return Kpcb;    break;
		case 41:		return Kpc_max;    break;
		case 42:		return Kpc_half;    break;
		case 43:		return i_pCa_max;    break;
		case 44:		return Km_pCa;    break;
		case 45:		return k_NaCa;    break;
		case 46:		return K_mNa;    break;
		case 47:		return Nao;    break;
		case 48:		return K_mCa;    break;
		case 49:		return Cao;    break;
		case 50:		return k_sat;    break;
		case 51:		return eta;    break;
		case 52:		return R;    break;
		case 53:		return T;    break;
		case 54:		return g_Cab;    break;
		case 55:		return g_Na;    break;
		case 56:		return Ko;    break;
		case 57:		return g_Nab;    break;
		case 58:		return g_Kto_f;    break;
		case 59:		return g_Kto_s;    break;
		case 60:		return g_Ks;    break;
		case 61:		return g_Kur;    break;
		case 62:		return g_Kss;    break;
		case 63:		return g_Kr;    break;
		case 64:		return kf;    break;
		case 65:		return kb;    break;
		case 66:		return i_NaK_max;    break;
		case 67:		return Km_Nai;    break;
		case 68:		return Km_Ko;    break;
		case 69:		return g_ClCa;    break;
		case 70:		return Km_Cl;    break;
		case 71:		return E_Cl;    break;
		default:	break;
		}
	}

	double Solveode::getFreeVariable()
	{
		return dtime;
	}

	//Get Methods - Variables

	Variables Solveode::get_Variables()
	{
		Variables v("|V#|Cai#|Cass#|CaJSR#|CaNSR#|P_RyR#|LTRPN_Ca#|HTRPN_Ca#|P_O1#|P_O2#|P_C2#|O#|C2#|C3#|C4#|I1#|I2#|I3#|Nai#|C_Na2#|C_Na1#|O_Na#|IF_Na#|I1_Na#|I2_Na#|IC_Na2#|IC_Na3#|Ki#|ato_f#|ito_f#|ato_s#|ito_s#|nKs#|aur#|iur#|aKss#|iKss#|C_K2#|C_K1#|O_K#|I_K#");
		return v;
	}
	Variables Solveode::get_Parameters()
	{
		Variables v("|time#|stim_amplitude#|stim_start#|stim_end#|stim_period#|stim_duration#|Acap#|Cm#|Vmyo#|F#|VJSR#|Vss#|VNSR#|CMDN_tot#|Km_CMDN#|CSQN_tot#|Km_CSQN#|v1#|tau_tr#|tau_xfer#|v2#|v3#|Km_up#|k_plus_htrpn#|HTRPN_tot#|k_plus_ltrpn#|LTRPN_tot#|k_minus_htrpn#|k_minus_ltrpn#|i_CaL_max#|k_plus_a#|n#|k_minus_b#|k_minus_c#|k_minus_a#|k_plus_b#|m#|k_plus_c#|g_CaL#|E_CaL#|Kpcb#|Kpc_max#|Kpc_half#|i_pCa_max#|Km_pCa#|k_NaCa#|K_mNa#|Nao#|K_mCa#|Cao#|k_sat#|eta#|R#|T#|g_Cab#|g_Na#|Ko#|g_Nab#|g_Kto_f#|g_Kto_s#|g_Ks#|g_Kur#|g_Kss#|g_Kr#|kf#|kb#|i_NaK_max#|Km_Nai#|Km_Ko#|g_ClCa#|Km_Cl#|E_Cl#");
		return v;
	}
	Variables Solveode::get_FreeVariable()
	{
		Variables v("|time#");
		return v;
	}

	void Solveode::setParametersFromFile(char *filename)
	{
		FILE *file;
		if((file = fopen(filename, "r")) == NULL)
		{
			fprintf(stderr,"ERROR - setParametersFromFile - Unable to open file %s\n", filename);
			exit(1);
		}
		double value;
		int k = 0;
		Variables v = get_Parameters();
		int s = v.getQuantity();
		for(;k<s;k++)
		{
			fscanf(file,"%lf", &value);
			setParameters(k, value);
		}
		fclose(file);
	}

	void Solveode::setVariablesFromFile(char *filename)
	{
		FILE *file;
		if((file = fopen(filename, "r")) == NULL)
		{
			fprintf(stderr,"ERROR - setVariablesFromFile - Unable to open file %s\n", filename);
			exit(1);
		}
		double value;
		int k = 0;
		Variables v = get_Variables();
		int s = v.getQuantity();
		for(;k<s;k++)
		{
			fscanf(file,"%lf", &value);
			setVariables(k, value);
		}
		fclose(file);
	}

	void Solveode::setFreeVariableFromFile(char *filename)
	{
		FILE *file;
		if((file = fopen(filename, "r")) == NULL)
		{
			fprintf(stderr,"ERROR - setFreeVariableFromFile - Unable to open file %s\n", filename);
			exit(1);
		}
		double value;
		fscanf(file,"%lf", &value);
			setFreeVariable(value);
		fclose(file);
	}
	double Solveode::jacobian(double final, double svRate)
	{
		time_new = time;

		int num_results__		= (int)(final / svRate);
		int iterations	= (int)(final / this->dtime);
		if(time_vec__ != NULL)free( time_vec__);
			time_vec__ = (double *)malloc(sizeof(double)*num_results__);
		V_old_ = V_ini_;
		if(V != NULL)free( V);
			V = (double *)malloc(sizeof(double)*num_results__);
		Cai_old_ = Cai_ini_;
		if(Cai != NULL)free( Cai);
			Cai = (double *)malloc(sizeof(double)*num_results__);
		Cass_old_ = Cass_ini_;
		if(Cass != NULL)free( Cass);
			Cass = (double *)malloc(sizeof(double)*num_results__);
		CaJSR_old_ = CaJSR_ini_;
		if(CaJSR != NULL)free( CaJSR);
			CaJSR = (double *)malloc(sizeof(double)*num_results__);
		CaNSR_old_ = CaNSR_ini_;
		if(CaNSR != NULL)free( CaNSR);
			CaNSR = (double *)malloc(sizeof(double)*num_results__);
		P_RyR_old_ = P_RyR_ini_;
		if(P_RyR != NULL)free( P_RyR);
			P_RyR = (double *)malloc(sizeof(double)*num_results__);
		LTRPN_Ca_old_ = LTRPN_Ca_ini_;
		if(LTRPN_Ca != NULL)free( LTRPN_Ca);
			LTRPN_Ca = (double *)malloc(sizeof(double)*num_results__);
		HTRPN_Ca_old_ = HTRPN_Ca_ini_;
		if(HTRPN_Ca != NULL)free( HTRPN_Ca);
			HTRPN_Ca = (double *)malloc(sizeof(double)*num_results__);
		P_O1_old_ = P_O1_ini_;
		if(P_O1 != NULL)free( P_O1);
			P_O1 = (double *)malloc(sizeof(double)*num_results__);
		P_O2_old_ = P_O2_ini_;
		if(P_O2 != NULL)free( P_O2);
			P_O2 = (double *)malloc(sizeof(double)*num_results__);
		P_C2_old_ = P_C2_ini_;
		if(P_C2 != NULL)free( P_C2);
			P_C2 = (double *)malloc(sizeof(double)*num_results__);
		O_old_ = O_ini_;
		if(O != NULL)free( O);
			O = (double *)malloc(sizeof(double)*num_results__);
		C2_old_ = C2_ini_;
		if(C2 != NULL)free( C2);
			C2 = (double *)malloc(sizeof(double)*num_results__);
		C3_old_ = C3_ini_;
		if(C3 != NULL)free( C3);
			C3 = (double *)malloc(sizeof(double)*num_results__);
		C4_old_ = C4_ini_;
		if(C4 != NULL)free( C4);
			C4 = (double *)malloc(sizeof(double)*num_results__);
		I1_old_ = I1_ini_;
		if(I1 != NULL)free( I1);
			I1 = (double *)malloc(sizeof(double)*num_results__);
		I2_old_ = I2_ini_;
		if(I2 != NULL)free( I2);
			I2 = (double *)malloc(sizeof(double)*num_results__);
		I3_old_ = I3_ini_;
		if(I3 != NULL)free( I3);
			I3 = (double *)malloc(sizeof(double)*num_results__);
		Nai_old_ = Nai_ini_;
		if(Nai != NULL)free( Nai);
			Nai = (double *)malloc(sizeof(double)*num_results__);
		C_Na2_old_ = C_Na2_ini_;
		if(C_Na2 != NULL)free( C_Na2);
			C_Na2 = (double *)malloc(sizeof(double)*num_results__);
		C_Na1_old_ = C_Na1_ini_;
		if(C_Na1 != NULL)free( C_Na1);
			C_Na1 = (double *)malloc(sizeof(double)*num_results__);
		O_Na_old_ = O_Na_ini_;
		if(O_Na != NULL)free( O_Na);
			O_Na = (double *)malloc(sizeof(double)*num_results__);
		IF_Na_old_ = IF_Na_ini_;
		if(IF_Na != NULL)free( IF_Na);
			IF_Na = (double *)malloc(sizeof(double)*num_results__);
		I1_Na_old_ = I1_Na_ini_;
		if(I1_Na != NULL)free( I1_Na);
			I1_Na = (double *)malloc(sizeof(double)*num_results__);
		I2_Na_old_ = I2_Na_ini_;
		if(I2_Na != NULL)free( I2_Na);
			I2_Na = (double *)malloc(sizeof(double)*num_results__);
		IC_Na2_old_ = IC_Na2_ini_;
		if(IC_Na2 != NULL)free( IC_Na2);
			IC_Na2 = (double *)malloc(sizeof(double)*num_results__);
		IC_Na3_old_ = IC_Na3_ini_;
		if(IC_Na3 != NULL)free( IC_Na3);
			IC_Na3 = (double *)malloc(sizeof(double)*num_results__);
		Ki_old_ = Ki_ini_;
		if(Ki != NULL)free( Ki);
			Ki = (double *)malloc(sizeof(double)*num_results__);
		ato_f_old_ = ato_f_ini_;
		if(ato_f != NULL)free( ato_f);
			ato_f = (double *)malloc(sizeof(double)*num_results__);
		ito_f_old_ = ito_f_ini_;
		if(ito_f != NULL)free( ito_f);
			ito_f = (double *)malloc(sizeof(double)*num_results__);
		ato_s_old_ = ato_s_ini_;
		if(ato_s != NULL)free( ato_s);
			ato_s = (double *)malloc(sizeof(double)*num_results__);
		ito_s_old_ = ito_s_ini_;
		if(ito_s != NULL)free( ito_s);
			ito_s = (double *)malloc(sizeof(double)*num_results__);
		nKs_old_ = nKs_ini_;
		if(nKs != NULL)free( nKs);
			nKs = (double *)malloc(sizeof(double)*num_results__);
		aur_old_ = aur_ini_;
		if(aur != NULL)free( aur);
			aur = (double *)malloc(sizeof(double)*num_results__);
		iur_old_ = iur_ini_;
		if(iur != NULL)free( iur);
			iur = (double *)malloc(sizeof(double)*num_results__);
		aKss_old_ = aKss_ini_;
		if(aKss != NULL)free( aKss);
			aKss = (double *)malloc(sizeof(double)*num_results__);
		iKss_old_ = iKss_ini_;
		if(iKss != NULL)free( iKss);
			iKss = (double *)malloc(sizeof(double)*num_results__);
		C_K2_old_ = C_K2_ini_;
		if(C_K2 != NULL)free( C_K2);
			C_K2 = (double *)malloc(sizeof(double)*num_results__);
		C_K1_old_ = C_K1_ini_;
		if(C_K1 != NULL)free( C_K1);
			C_K1 = (double *)malloc(sizeof(double)*num_results__);
		O_K_old_ = O_K_ini_;
		if(O_K != NULL)free( O_K);
			O_K = (double *)malloc(sizeof(double)*num_results__);
		I_K_old_ = I_K_ini_;
		if(I_K != NULL)free( I_K);
			I_K = (double *)malloc(sizeof(double)*num_results__);
		this->timeSaving = dtime;

		double diff=0;
		int counter=0;
		for (int i = 0; i< iterations;i++ )
		{
			this->time_new += dtime;

			rightHandSideFunction.function(this);
			this->V_new_ = this->V_old_ + this->V_lado_direito_ * this->dtime;
			this->Cai_new_ = this->Cai_old_ + this->Cai_lado_direito_ * this->dtime;
			this->Cass_new_ = this->Cass_old_ + this->Cass_lado_direito_ * this->dtime;
			this->CaJSR_new_ = this->CaJSR_old_ + this->CaJSR_lado_direito_ * this->dtime;
			this->CaNSR_new_ = this->CaNSR_old_ + this->CaNSR_lado_direito_ * this->dtime;
			this->P_RyR_new_ = this->P_RyR_old_ + this->P_RyR_lado_direito_ * this->dtime;
			this->LTRPN_Ca_new_ = this->LTRPN_Ca_old_ + this->LTRPN_Ca_lado_direito_ * this->dtime;
			this->HTRPN_Ca_new_ = this->HTRPN_Ca_old_ + this->HTRPN_Ca_lado_direito_ * this->dtime;
			this->P_O1_new_ = this->P_O1_old_ + this->P_O1_lado_direito_ * this->dtime;
			this->P_O2_new_ = this->P_O2_old_ + this->P_O2_lado_direito_ * this->dtime;
			this->P_C2_new_ = this->P_C2_old_ + this->P_C2_lado_direito_ * this->dtime;
			this->O_new_ = this->O_old_ + this->O_lado_direito_ * this->dtime;
			this->C2_new_ = this->C2_old_ + this->C2_lado_direito_ * this->dtime;
			this->C3_new_ = this->C3_old_ + this->C3_lado_direito_ * this->dtime;
			this->C4_new_ = this->C4_old_ + this->C4_lado_direito_ * this->dtime;
			this->I1_new_ = this->I1_old_ + this->I1_lado_direito_ * this->dtime;
			this->I2_new_ = this->I2_old_ + this->I2_lado_direito_ * this->dtime;
			this->I3_new_ = this->I3_old_ + this->I3_lado_direito_ * this->dtime;
			this->Nai_new_ = this->Nai_old_ + this->Nai_lado_direito_ * this->dtime;
			this->C_Na2_new_ = this->C_Na2_old_ + this->C_Na2_lado_direito_ * this->dtime;
			this->C_Na1_new_ = this->C_Na1_old_ + this->C_Na1_lado_direito_ * this->dtime;
			this->O_Na_new_ = this->O_Na_old_ + this->O_Na_lado_direito_ * this->dtime;
			this->IF_Na_new_ = this->IF_Na_old_ + this->IF_Na_lado_direito_ * this->dtime;
			this->I1_Na_new_ = this->I1_Na_old_ + this->I1_Na_lado_direito_ * this->dtime;
			this->I2_Na_new_ = this->I2_Na_old_ + this->I2_Na_lado_direito_ * this->dtime;
			this->IC_Na2_new_ = this->IC_Na2_old_ + this->IC_Na2_lado_direito_ * this->dtime;
			this->IC_Na3_new_ = this->IC_Na3_old_ + this->IC_Na3_lado_direito_ * this->dtime;
			this->Ki_new_ = this->Ki_old_ + this->Ki_lado_direito_ * this->dtime;
			this->ato_f_new_ = this->ato_f_old_ + this->ato_f_lado_direito_ * this->dtime;
			this->ito_f_new_ = this->ito_f_old_ + this->ito_f_lado_direito_ * this->dtime;
			this->ato_s_new_ = this->ato_s_old_ + this->ato_s_lado_direito_ * this->dtime;
			this->ito_s_new_ = this->ito_s_old_ + this->ito_s_lado_direito_ * this->dtime;
			this->nKs_new_ = this->nKs_old_ + this->nKs_lado_direito_ * this->dtime;
			this->aur_new_ = this->aur_old_ + this->aur_lado_direito_ * this->dtime;
			this->iur_new_ = this->iur_old_ + this->iur_lado_direito_ * this->dtime;
			this->aKss_new_ = this->aKss_old_ + this->aKss_lado_direito_ * this->dtime;
			this->iKss_new_ = this->iKss_old_ + this->iKss_lado_direito_ * this->dtime;
			this->C_K2_new_ = this->C_K2_old_ + this->C_K2_lado_direito_ * this->dtime;
			this->C_K1_new_ = this->C_K1_old_ + this->C_K1_lado_direito_ * this->dtime;
			this->O_K_new_ = this->O_K_old_ + this->O_K_lado_direito_ * this->dtime;
			this->I_K_new_ = this->I_K_old_ + this->I_K_lado_direito_ * this->dtime;
			diff =  _agos_round(this->time_new - timeSaving, 6);
			if(diff==0){
				this->timeSaving += svRate;
				time_vec__[counter] = this->time_new;
				V[counter] = this->V_new_;
				Cai[counter] = this->Cai_new_;
				Cass[counter] = this->Cass_new_;
				CaJSR[counter] = this->CaJSR_new_;
				CaNSR[counter] = this->CaNSR_new_;
				P_RyR[counter] = this->P_RyR_new_;
				LTRPN_Ca[counter] = this->LTRPN_Ca_new_;
				HTRPN_Ca[counter] = this->HTRPN_Ca_new_;
				P_O1[counter] = this->P_O1_new_;
				P_O2[counter] = this->P_O2_new_;
				P_C2[counter] = this->P_C2_new_;
				O[counter] = this->O_new_;
				C2[counter] = this->C2_new_;
				C3[counter] = this->C3_new_;
				C4[counter] = this->C4_new_;
				I1[counter] = this->I1_new_;
				I2[counter] = this->I2_new_;
				I3[counter] = this->I3_new_;
				Nai[counter] = this->Nai_new_;
				C_Na2[counter] = this->C_Na2_new_;
				C_Na1[counter] = this->C_Na1_new_;
				O_Na[counter] = this->O_Na_new_;
				IF_Na[counter] = this->IF_Na_new_;
				I1_Na[counter] = this->I1_Na_new_;
				I2_Na[counter] = this->I2_Na_new_;
				IC_Na2[counter] = this->IC_Na2_new_;
				IC_Na3[counter] = this->IC_Na3_new_;
				Ki[counter] = this->Ki_new_;
				ato_f[counter] = this->ato_f_new_;
				ito_f[counter] = this->ito_f_new_;
				ato_s[counter] = this->ato_s_new_;
				ito_s[counter] = this->ito_s_new_;
				nKs[counter] = this->nKs_new_;
				aur[counter] = this->aur_new_;
				iur[counter] = this->iur_new_;
				aKss[counter] = this->aKss_new_;
				iKss[counter] = this->iKss_new_;
				C_K2[counter] = this->C_K2_new_;
				C_K1[counter] = this->C_K1_new_;
				O_K[counter] = this->O_K_new_;
				I_K[counter] = this->I_K_new_;
				counter++;
			}
			this->V_old_ = this->V_new_;
			this->Cai_old_ = this->Cai_new_;
			this->Cass_old_ = this->Cass_new_;
			this->CaJSR_old_ = this->CaJSR_new_;
			this->CaNSR_old_ = this->CaNSR_new_;
			this->P_RyR_old_ = this->P_RyR_new_;
			this->LTRPN_Ca_old_ = this->LTRPN_Ca_new_;
			this->HTRPN_Ca_old_ = this->HTRPN_Ca_new_;
			this->P_O1_old_ = this->P_O1_new_;
			this->P_O2_old_ = this->P_O2_new_;
			this->P_C2_old_ = this->P_C2_new_;
			this->O_old_ = this->O_new_;
			this->C2_old_ = this->C2_new_;
			this->C3_old_ = this->C3_new_;
			this->C4_old_ = this->C4_new_;
			this->I1_old_ = this->I1_new_;
			this->I2_old_ = this->I2_new_;
			this->I3_old_ = this->I3_new_;
			this->Nai_old_ = this->Nai_new_;
			this->C_Na2_old_ = this->C_Na2_new_;
			this->C_Na1_old_ = this->C_Na1_new_;
			this->O_Na_old_ = this->O_Na_new_;
			this->IF_Na_old_ = this->IF_Na_new_;
			this->I1_Na_old_ = this->I1_Na_new_;
			this->I2_Na_old_ = this->I2_Na_new_;
			this->IC_Na2_old_ = this->IC_Na2_new_;
			this->IC_Na3_old_ = this->IC_Na3_new_;
			this->Ki_old_ = this->Ki_new_;
			this->ato_f_old_ = this->ato_f_new_;
			this->ito_f_old_ = this->ito_f_new_;
			this->ato_s_old_ = this->ato_s_new_;
			this->ito_s_old_ = this->ito_s_new_;
			this->nKs_old_ = this->nKs_new_;
			this->aur_old_ = this->aur_new_;
			this->iur_old_ = this->iur_new_;
			this->aKss_old_ = this->aKss_new_;
			this->iKss_old_ = this->iKss_new_;
			this->C_K2_old_ = this->C_K2_new_;
			this->C_K1_old_ = this->C_K1_new_;
			this->O_K_old_ = this->O_K_new_;
			this->I_K_old_ = this->I_K_new_;
		}
		double h_jac_[numEDO];
		double quociente = 1000.0;
		h_jac_[0] = fabs(_agos_min(V, num_results__) / _agos_max(V, num_results__) );
		h_jac_[1] = fabs(_agos_min(Cai, num_results__) / _agos_max(Cai, num_results__) );
		h_jac_[2] = fabs(_agos_min(Cass, num_results__) / _agos_max(Cass, num_results__) );
		h_jac_[3] = fabs(_agos_min(CaJSR, num_results__) / _agos_max(CaJSR, num_results__) );
		h_jac_[4] = fabs(_agos_min(CaNSR, num_results__) / _agos_max(CaNSR, num_results__) );
		h_jac_[5] = fabs(_agos_min(P_RyR, num_results__) / _agos_max(P_RyR, num_results__) );
		h_jac_[6] = fabs(_agos_min(LTRPN_Ca, num_results__) / _agos_max(LTRPN_Ca, num_results__) );
		h_jac_[7] = fabs(_agos_min(HTRPN_Ca, num_results__) / _agos_max(HTRPN_Ca, num_results__) );
		h_jac_[8] = fabs(_agos_min(P_O1, num_results__) / _agos_max(P_O1, num_results__) );
		h_jac_[9] = fabs(_agos_min(P_O2, num_results__) / _agos_max(P_O2, num_results__) );
		h_jac_[10] = fabs(_agos_min(P_C2, num_results__) / _agos_max(P_C2, num_results__) );
		h_jac_[11] = fabs(_agos_min(O, num_results__) / _agos_max(O, num_results__) );
		h_jac_[12] = fabs(_agos_min(C2, num_results__) / _agos_max(C2, num_results__) );
		h_jac_[13] = fabs(_agos_min(C3, num_results__) / _agos_max(C3, num_results__) );
		h_jac_[14] = fabs(_agos_min(C4, num_results__) / _agos_max(C4, num_results__) );
		h_jac_[15] = fabs(_agos_min(I1, num_results__) / _agos_max(I1, num_results__) );
		h_jac_[16] = fabs(_agos_min(I2, num_results__) / _agos_max(I2, num_results__) );
		h_jac_[17] = fabs(_agos_min(I3, num_results__) / _agos_max(I3, num_results__) );
		h_jac_[18] = fabs(_agos_min(Nai, num_results__) / _agos_max(Nai, num_results__) );
		h_jac_[19] = fabs(_agos_min(C_Na2, num_results__) / _agos_max(C_Na2, num_results__) );
		h_jac_[20] = fabs(_agos_min(C_Na1, num_results__) / _agos_max(C_Na1, num_results__) );
		h_jac_[21] = fabs(_agos_min(O_Na, num_results__) / _agos_max(O_Na, num_results__) );
		h_jac_[22] = fabs(_agos_min(IF_Na, num_results__) / _agos_max(IF_Na, num_results__) );
		h_jac_[23] = fabs(_agos_min(I1_Na, num_results__) / _agos_max(I1_Na, num_results__) );
		h_jac_[24] = fabs(_agos_min(I2_Na, num_results__) / _agos_max(I2_Na, num_results__) );
		h_jac_[25] = fabs(_agos_min(IC_Na2, num_results__) / _agos_max(IC_Na2, num_results__) );
		h_jac_[26] = fabs(_agos_min(IC_Na3, num_results__) / _agos_max(IC_Na3, num_results__) );
		h_jac_[27] = fabs(_agos_min(Ki, num_results__) / _agos_max(Ki, num_results__) );
		h_jac_[28] = fabs(_agos_min(ato_f, num_results__) / _agos_max(ato_f, num_results__) );
		h_jac_[29] = fabs(_agos_min(ito_f, num_results__) / _agos_max(ito_f, num_results__) );
		h_jac_[30] = fabs(_agos_min(ato_s, num_results__) / _agos_max(ato_s, num_results__) );
		h_jac_[31] = fabs(_agos_min(ito_s, num_results__) / _agos_max(ito_s, num_results__) );
		h_jac_[32] = fabs(_agos_min(nKs, num_results__) / _agos_max(nKs, num_results__) );
		h_jac_[33] = fabs(_agos_min(aur, num_results__) / _agos_max(aur, num_results__) );
		h_jac_[34] = fabs(_agos_min(iur, num_results__) / _agos_max(iur, num_results__) );
		h_jac_[35] = fabs(_agos_min(aKss, num_results__) / _agos_max(aKss, num_results__) );
		h_jac_[36] = fabs(_agos_min(iKss, num_results__) / _agos_max(iKss, num_results__) );
		h_jac_[37] = fabs(_agos_min(C_K2, num_results__) / _agos_max(C_K2, num_results__) );
		h_jac_[38] = fabs(_agos_min(C_K1, num_results__) / _agos_max(C_K1, num_results__) );
		h_jac_[39] = fabs(_agos_min(O_K, num_results__) / _agos_max(O_K, num_results__) );
		h_jac_[40] = fabs(_agos_min(I_K, num_results__) / _agos_max(I_K, num_results__) );
		for(int l=0;l<numEDO;l++){
			h_jac_[l] = (h_jac_[l]==0 || h_jac_[l]==AGOS_NAN || h_jac_[l]==AGOS_INF)?this->dtime:h_jac_[l];
		}
		this->timeSaving = this->dtime;

		this->time_new = this->time;

		V_old_ = V_ini_;
		Cai_old_ = Cai_ini_;
		Cass_old_ = Cass_ini_;
		CaJSR_old_ = CaJSR_ini_;
		CaNSR_old_ = CaNSR_ini_;
		P_RyR_old_ = P_RyR_ini_;
		LTRPN_Ca_old_ = LTRPN_Ca_ini_;
		HTRPN_Ca_old_ = HTRPN_Ca_ini_;
		P_O1_old_ = P_O1_ini_;
		P_O2_old_ = P_O2_ini_;
		P_C2_old_ = P_C2_ini_;
		O_old_ = O_ini_;
		C2_old_ = C2_ini_;
		C3_old_ = C3_ini_;
		C4_old_ = C4_ini_;
		I1_old_ = I1_ini_;
		I2_old_ = I2_ini_;
		I3_old_ = I3_ini_;
		Nai_old_ = Nai_ini_;
		C_Na2_old_ = C_Na2_ini_;
		C_Na1_old_ = C_Na1_ini_;
		O_Na_old_ = O_Na_ini_;
		IF_Na_old_ = IF_Na_ini_;
		I1_Na_old_ = I1_Na_ini_;
		I2_Na_old_ = I2_Na_ini_;
		IC_Na2_old_ = IC_Na2_ini_;
		IC_Na3_old_ = IC_Na3_ini_;
		Ki_old_ = Ki_ini_;
		ato_f_old_ = ato_f_ini_;
		ito_f_old_ = ito_f_ini_;
		ato_s_old_ = ato_s_ini_;
		ito_s_old_ = ito_s_ini_;
		nKs_old_ = nKs_ini_;
		aur_old_ = aur_ini_;
		iur_old_ = iur_ini_;
		aKss_old_ = aKss_ini_;
		iKss_old_ = iKss_ini_;
		C_K2_old_ = C_K2_ini_;
		C_K1_old_ = C_K1_ini_;
		O_K_old_ = O_K_ini_;
		I_K_old_ = I_K_ini_;
		double jacobian[numEDO][numEDO];
		double edos_new_aux_[numEDO];
		double edos_aux_[numEDO];
		FILE *filejac;
		filejac = fopen("dat/jacobian.dat", "wb");
		int counter2=0;
		for (int i = 0; i< iterations;i++ ){
			time_new += dtime;
			rightHandSideFunction.function(this);
			this->V_new_ = this->dtime*(this->V_lado_direito_) + this->V_old_;
			this->Cai_new_ = this->dtime*(this->Cai_lado_direito_) + this->Cai_old_;
			this->Cass_new_ = this->dtime*(this->Cass_lado_direito_) + this->Cass_old_;
			this->CaJSR_new_ = this->dtime*(this->CaJSR_lado_direito_) + this->CaJSR_old_;
			this->CaNSR_new_ = this->dtime*(this->CaNSR_lado_direito_) + this->CaNSR_old_;
			this->P_RyR_new_ = this->dtime*(this->P_RyR_lado_direito_) + this->P_RyR_old_;
			this->LTRPN_Ca_new_ = this->dtime*(this->LTRPN_Ca_lado_direito_) + this->LTRPN_Ca_old_;
			this->HTRPN_Ca_new_ = this->dtime*(this->HTRPN_Ca_lado_direito_) + this->HTRPN_Ca_old_;
			this->P_O1_new_ = this->dtime*(this->P_O1_lado_direito_) + this->P_O1_old_;
			this->P_O2_new_ = this->dtime*(this->P_O2_lado_direito_) + this->P_O2_old_;
			this->P_C2_new_ = this->dtime*(this->P_C2_lado_direito_) + this->P_C2_old_;
			this->O_new_ = this->dtime*(this->O_lado_direito_) + this->O_old_;
			this->C2_new_ = this->dtime*(this->C2_lado_direito_) + this->C2_old_;
			this->C3_new_ = this->dtime*(this->C3_lado_direito_) + this->C3_old_;
			this->C4_new_ = this->dtime*(this->C4_lado_direito_) + this->C4_old_;
			this->I1_new_ = this->dtime*(this->I1_lado_direito_) + this->I1_old_;
			this->I2_new_ = this->dtime*(this->I2_lado_direito_) + this->I2_old_;
			this->I3_new_ = this->dtime*(this->I3_lado_direito_) + this->I3_old_;
			this->Nai_new_ = this->dtime*(this->Nai_lado_direito_) + this->Nai_old_;
			this->C_Na2_new_ = this->dtime*(this->C_Na2_lado_direito_) + this->C_Na2_old_;
			this->C_Na1_new_ = this->dtime*(this->C_Na1_lado_direito_) + this->C_Na1_old_;
			this->O_Na_new_ = this->dtime*(this->O_Na_lado_direito_) + this->O_Na_old_;
			this->IF_Na_new_ = this->dtime*(this->IF_Na_lado_direito_) + this->IF_Na_old_;
			this->I1_Na_new_ = this->dtime*(this->I1_Na_lado_direito_) + this->I1_Na_old_;
			this->I2_Na_new_ = this->dtime*(this->I2_Na_lado_direito_) + this->I2_Na_old_;
			this->IC_Na2_new_ = this->dtime*(this->IC_Na2_lado_direito_) + this->IC_Na2_old_;
			this->IC_Na3_new_ = this->dtime*(this->IC_Na3_lado_direito_) + this->IC_Na3_old_;
			this->Ki_new_ = this->dtime*(this->Ki_lado_direito_) + this->Ki_old_;
			this->ato_f_new_ = this->dtime*(this->ato_f_lado_direito_) + this->ato_f_old_;
			this->ito_f_new_ = this->dtime*(this->ito_f_lado_direito_) + this->ito_f_old_;
			this->ato_s_new_ = this->dtime*(this->ato_s_lado_direito_) + this->ato_s_old_;
			this->ito_s_new_ = this->dtime*(this->ito_s_lado_direito_) + this->ito_s_old_;
			this->nKs_new_ = this->dtime*(this->nKs_lado_direito_) + this->nKs_old_;
			this->aur_new_ = this->dtime*(this->aur_lado_direito_) + this->aur_old_;
			this->iur_new_ = this->dtime*(this->iur_lado_direito_) + this->iur_old_;
			this->aKss_new_ = this->dtime*(this->aKss_lado_direito_) + this->aKss_old_;
			this->iKss_new_ = this->dtime*(this->iKss_lado_direito_) + this->iKss_old_;
			this->C_K2_new_ = this->dtime*(this->C_K2_lado_direito_) + this->C_K2_old_;
			this->C_K1_new_ = this->dtime*(this->C_K1_lado_direito_) + this->C_K1_old_;
			this->O_K_new_ = this->dtime*(this->O_K_lado_direito_) + this->O_K_old_;
			this->I_K_new_ = this->dtime*(this->I_K_lado_direito_) + this->I_K_old_;
			diff =  _agos_round(this->time_new - timeSaving, 6);
			if(diff==0){
				this->timeSaving += svRate;
				V[counter2] = V_new_;
				Cai[counter2] = Cai_new_;
				Cass[counter2] = Cass_new_;
				CaJSR[counter2] = CaJSR_new_;
				CaNSR[counter2] = CaNSR_new_;
				P_RyR[counter2] = P_RyR_new_;
				LTRPN_Ca[counter2] = LTRPN_Ca_new_;
				HTRPN_Ca[counter2] = HTRPN_Ca_new_;
				P_O1[counter2] = P_O1_new_;
				P_O2[counter2] = P_O2_new_;
				P_C2[counter2] = P_C2_new_;
				O[counter2] = O_new_;
				C2[counter2] = C2_new_;
				C3[counter2] = C3_new_;
				C4[counter2] = C4_new_;
				I1[counter2] = I1_new_;
				I2[counter2] = I2_new_;
				I3[counter2] = I3_new_;
				Nai[counter2] = Nai_new_;
				C_Na2[counter2] = C_Na2_new_;
				C_Na1[counter2] = C_Na1_new_;
				O_Na[counter2] = O_Na_new_;
				IF_Na[counter2] = IF_Na_new_;
				I1_Na[counter2] = I1_Na_new_;
				I2_Na[counter2] = I2_Na_new_;
				IC_Na2[counter2] = IC_Na2_new_;
				IC_Na3[counter2] = IC_Na3_new_;
				Ki[counter2] = Ki_new_;
				ato_f[counter2] = ato_f_new_;
				ito_f[counter2] = ito_f_new_;
				ato_s[counter2] = ato_s_new_;
				ito_s[counter2] = ito_s_new_;
				nKs[counter2] = nKs_new_;
				aur[counter2] = aur_new_;
				iur[counter2] = iur_new_;
				aKss[counter2] = aKss_new_;
				iKss[counter2] = iKss_new_;
				C_K2[counter2] = C_K2_new_;
				C_K1[counter2] = C_K1_new_;
				O_K[counter2] = O_K_new_;
				I_K[counter2] = I_K_new_;
				time_vec__[counter2] = time_new;
				counter2++;
				//salva os valores das variaveis diferenciaveis em auxiliares
				for(int l=0;l<numEDO;l++){
					edos_aux_[l] = this->getLadoDireito(l);
					edos_new_aux_[l] = this->getVariables(l);
				}
				//para cada coluna
				for(int k=0;k<numEDO;k++){
					//escolhe uma variavel k, e calcula a derivada de todas as equacoes em relacao a k
					this->setVariables(k, edos_new_aux_[k]+h_jac_[k]);
					//calcula tudo de novo com o novo valor de k
					rightHandSideFunction.function(this);
					for(int j=0;j<numEDO;j++){//para cada linha 
						jacobian[j][k] = (this->getLadoDireito(j) - edos_aux_[j])/h_jac_[k];
						fprintf(filejac,"%f\t", jacobian[j][k]);
					}
					fprintf(filejac,"\n");
					//agora tem que voltar para o que estava antes, sem somar dtime a variavel k
					for(int l=0;l<numEDO;l++){
						this->setVariables(l, edos_new_aux_[l]);
					}
					rightHandSideFunction.function(this);
				}
				fprintf(filejac,"\n");
			}
			V_old_ = V_new_;
			Cai_old_ = Cai_new_;
			Cass_old_ = Cass_new_;
			CaJSR_old_ = CaJSR_new_;
			CaNSR_old_ = CaNSR_new_;
			P_RyR_old_ = P_RyR_new_;
			LTRPN_Ca_old_ = LTRPN_Ca_new_;
			HTRPN_Ca_old_ = HTRPN_Ca_new_;
			P_O1_old_ = P_O1_new_;
			P_O2_old_ = P_O2_new_;
			P_C2_old_ = P_C2_new_;
			O_old_ = O_new_;
			C2_old_ = C2_new_;
			C3_old_ = C3_new_;
			C4_old_ = C4_new_;
			I1_old_ = I1_new_;
			I2_old_ = I2_new_;
			I3_old_ = I3_new_;
			Nai_old_ = Nai_new_;
			C_Na2_old_ = C_Na2_new_;
			C_Na1_old_ = C_Na1_new_;
			O_Na_old_ = O_Na_new_;
			IF_Na_old_ = IF_Na_new_;
			I1_Na_old_ = I1_Na_new_;
			I2_Na_old_ = I2_Na_new_;
			IC_Na2_old_ = IC_Na2_new_;
			IC_Na3_old_ = IC_Na3_new_;
			Ki_old_ = Ki_new_;
			ato_f_old_ = ato_f_new_;
			ito_f_old_ = ito_f_new_;
			ato_s_old_ = ato_s_new_;
			ito_s_old_ = ito_s_new_;
			nKs_old_ = nKs_new_;
			aur_old_ = aur_new_;
			iur_old_ = iur_new_;
			aKss_old_ = aKss_new_;
			iKss_old_ = iKss_new_;
			C_K2_old_ = C_K2_new_;
			C_K1_old_ = C_K1_new_;
			O_K_old_ = O_K_new_;
			I_K_old_ = I_K_new_;
		}
		fclose(filejac);
			//FILE *fileptr;
			//char filename[12];
			//sprintf(filename,"%dthread.dat",numThreads);
			//fileptr = fopen(filename, "wb");
			//for(int i =0;i<num_results__;i++){
			//    fprintf(fileptr,"%f %f\n", time_vec__[i], V[i]);
			//}
			//fclose(fileptr);
		return 0;
	}

	double* Solveode::getIndependentVar()
	{
		return time_vec__;
	}

	double* Solveode::getSolution(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return V;    break;
		case 1:		return Cai;    break;
		case 2:		return Cass;    break;
		case 3:		return CaJSR;    break;
		case 4:		return CaNSR;    break;
		case 5:		return P_RyR;    break;
		case 6:		return LTRPN_Ca;    break;
		case 7:		return HTRPN_Ca;    break;
		case 8:		return P_O1;    break;
		case 9:		return P_O2;    break;
		case 10:		return P_C2;    break;
		case 11:		return O;    break;
		case 12:		return C2;    break;
		case 13:		return C3;    break;
		case 14:		return C4;    break;
		case 15:		return I1;    break;
		case 16:		return I2;    break;
		case 17:		return I3;    break;
		case 18:		return Nai;    break;
		case 19:		return C_Na2;    break;
		case 20:		return C_Na1;    break;
		case 21:		return O_Na;    break;
		case 22:		return IF_Na;    break;
		case 23:		return I1_Na;    break;
		case 24:		return I2_Na;    break;
		case 25:		return IC_Na2;    break;
		case 26:		return IC_Na3;    break;
		case 27:		return Ki;    break;
		case 28:		return ato_f;    break;
		case 29:		return ito_f;    break;
		case 30:		return ato_s;    break;
		case 31:		return ito_s;    break;
		case 32:		return nKs;    break;
		case 33:		return aur;    break;
		case 34:		return iur;    break;
		case 35:		return aKss;    break;
		case 36:		return iKss;    break;
		case 37:		return C_K2;    break;
		case 38:		return C_K1;    break;
		case 39:		return O_K;    break;
		case 40:		return I_K;    break;
		default:	return NULL;    break;
		}
	}
	void Solveode::reInitCVODE()
	{

		flag__ = CVodeReInit(cvode_mem_cvode__, f__, time, dependent_variable__, CV_SV, reltol__, &abstol__);
		if (check_flag(&flag__, "CVodeReInit", 1))
			exit(1);
	}

	void Solveode::setCVODEMaxStep(double maxstep)
	{

			CVodeSetMaxNumSteps(cvode_mem_cvode__, 1000000);

			realtype initDT = this->getFreeVariable();
			CVodeSetInitStep(cvode_mem_cvode__, initDT);
			flag__ = CVodeSetMaxStep(cvode_mem_cvode__, maxstep);
			if (check_flag(&flag__, "CVodeSetMaxStep", 1))
				exit(1);
	}


	double Solveode::ifnumber_0(){
		if(((time_new>=stim_start)&&(time_new<=stim_end)&&(((time_new-stim_start)-(floor(((time_new-stim_start)/stim_period))*stim_period))<=stim_duration))){
			return (stim_amplitude);
		}else{
			return (0.0000000000e+00);
		}
	}

static int check_flag(void *flagvalue, char *funcname, int opt){
	int *errflag;
	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
		return(1);}
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",funcname, *errflag);
			return(1); }}
	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
		return(1); }
	return 0;
}


	void Solveode::solveCVODE(int firstcall__, int steps__, int num_results__, int method__, char *fileName,int cv_method__)
	{
		FILE *f = fopen(fileName, "wb");
		if(!f){
			fprintf(stderr,"ERROR - solveToFile - Unable to open file %s\n",fileName);
			exit(1);
		}
		double dtL, dtM, dtMax=0.0,  dtMin=9990.0 ;
		dtM = 0.0;
		int iout = 1;
		time_new = time;
		realtype tout = time_new+dtime;

		realtype cvodeDT;
		int counter_it__ = 0;
		int offset_step = steps__ / num_results__;
		int aux = steps__%num_results__;
		int num_iterations_bak = steps__ + 1;
		if(firstcall__==1){
			if(steps__ <= 0)
				steps__ = 1;
			dependent_variable__ = N_VNew_Serial(41);
			if(check_flag((void *)dependent_variable__, "N_VNew_Serial", 0))
			exit(1);
			depvar__ = (double*)malloc(sizeof(double)*41);
			if(depvar__ == NULL){
			fprintf(stderr, "ERROR Cannot allocate memory for depvar__\n");
			exit(0);
			}
			NV_Ith_S(dependent_variable__, 0) = V_ini_;
			NV_Ith_S(dependent_variable__, 1) = Cai_ini_;
			NV_Ith_S(dependent_variable__, 2) = Cass_ini_;
			NV_Ith_S(dependent_variable__, 3) = CaJSR_ini_;
			NV_Ith_S(dependent_variable__, 4) = CaNSR_ini_;
			NV_Ith_S(dependent_variable__, 5) = P_RyR_ini_;
			NV_Ith_S(dependent_variable__, 6) = LTRPN_Ca_ini_;
			NV_Ith_S(dependent_variable__, 7) = HTRPN_Ca_ini_;
			NV_Ith_S(dependent_variable__, 8) = P_O1_ini_;
			NV_Ith_S(dependent_variable__, 9) = P_O2_ini_;
			NV_Ith_S(dependent_variable__, 10) = P_C2_ini_;
			NV_Ith_S(dependent_variable__, 11) = O_ini_;
			NV_Ith_S(dependent_variable__, 12) = C2_ini_;
			NV_Ith_S(dependent_variable__, 13) = C3_ini_;
			NV_Ith_S(dependent_variable__, 14) = C4_ini_;
			NV_Ith_S(dependent_variable__, 15) = I1_ini_;
			NV_Ith_S(dependent_variable__, 16) = I2_ini_;
			NV_Ith_S(dependent_variable__, 17) = I3_ini_;
			NV_Ith_S(dependent_variable__, 18) = Nai_ini_;
			NV_Ith_S(dependent_variable__, 19) = C_Na2_ini_;
			NV_Ith_S(dependent_variable__, 20) = C_Na1_ini_;
			NV_Ith_S(dependent_variable__, 21) = O_Na_ini_;
			NV_Ith_S(dependent_variable__, 22) = IF_Na_ini_;
			NV_Ith_S(dependent_variable__, 23) = I1_Na_ini_;
			NV_Ith_S(dependent_variable__, 24) = I2_Na_ini_;
			NV_Ith_S(dependent_variable__, 25) = IC_Na2_ini_;
			NV_Ith_S(dependent_variable__, 26) = IC_Na3_ini_;
			NV_Ith_S(dependent_variable__, 27) = Ki_ini_;
			NV_Ith_S(dependent_variable__, 28) = ato_f_ini_;
			NV_Ith_S(dependent_variable__, 29) = ito_f_ini_;
			NV_Ith_S(dependent_variable__, 30) = ato_s_ini_;
			NV_Ith_S(dependent_variable__, 31) = ito_s_ini_;
			NV_Ith_S(dependent_variable__, 32) = nKs_ini_;
			NV_Ith_S(dependent_variable__, 33) = aur_ini_;
			NV_Ith_S(dependent_variable__, 34) = iur_ini_;
			NV_Ith_S(dependent_variable__, 35) = aKss_ini_;
			NV_Ith_S(dependent_variable__, 36) = iKss_ini_;
			NV_Ith_S(dependent_variable__, 37) = C_K2_ini_;
			NV_Ith_S(dependent_variable__, 38) = C_K1_ini_;
			NV_Ith_S(dependent_variable__, 39) = O_K_ini_;
			NV_Ith_S(dependent_variable__, 40) = I_K_ini_;
			it_countx = 0;
			int nonlineariteration__;
			if(method__ == CV_BDF) nonlineariteration__ = CV_NEWTON;
			else nonlineariteration__ = CV_FUNCTIONAL;
			cvode_mem_cvode__ = CVodeCreate(method__, nonlineariteration__);
			if (check_flag((void *)cvode_mem_cvode__, "CVodeCreate", 0))
			exit(1);
			flag__ = CVodeMalloc(cvode_mem_cvode__, f__, time, dependent_variable__, CV_SS, reltol__, &abstol__);
			if (check_flag(&flag__, "CVodeMalloc", 1))
			exit(1);
			this->setCVODEMaxStep(this->maxStep);
			switch(cv_method__){
			case 1 :
			flag__ = CVDense(cvode_mem_cvode__, 41);
			if (check_flag(&flag__, "CVDense", 1))	exit(1);
			break;
			case 2:
			flag__ = CVDiag(cvode_mem_cvode__);
			if (check_flag(&flag__, "CVDiag", 1))	exit(1);
			break;
			case 3:
			flag__ = CVBand(cvode_mem_cvode__, 41, NULL, NULL);
			if (check_flag(&flag__, "CVBand", 1))	exit(1);
			break;
			}
			CVodeSetFdata(cvode_mem_cvode__, (void*)this);
			V_old_ = V_ini_;
			Cai_old_ = Cai_ini_;
			Cass_old_ = Cass_ini_;
			CaJSR_old_ = CaJSR_ini_;
			CaNSR_old_ = CaNSR_ini_;
			P_RyR_old_ = P_RyR_ini_;
			LTRPN_Ca_old_ = LTRPN_Ca_ini_;
			HTRPN_Ca_old_ = HTRPN_Ca_ini_;
			P_O1_old_ = P_O1_ini_;
			P_O2_old_ = P_O2_ini_;
			P_C2_old_ = P_C2_ini_;
			O_old_ = O_ini_;
			C2_old_ = C2_ini_;
			C3_old_ = C3_ini_;
			C4_old_ = C4_ini_;
			I1_old_ = I1_ini_;
			I2_old_ = I2_ini_;
			I3_old_ = I3_ini_;
			Nai_old_ = Nai_ini_;
			C_Na2_old_ = C_Na2_ini_;
			C_Na1_old_ = C_Na1_ini_;
			O_Na_old_ = O_Na_ini_;
			IF_Na_old_ = IF_Na_ini_;
			I1_Na_old_ = I1_Na_ini_;
			I2_Na_old_ = I2_Na_ini_;
			IC_Na2_old_ = IC_Na2_ini_;
			IC_Na3_old_ = IC_Na3_ini_;
			Ki_old_ = Ki_ini_;
			ato_f_old_ = ato_f_ini_;
			ito_f_old_ = ito_f_ini_;
			ato_s_old_ = ato_s_ini_;
			ito_s_old_ = ito_s_ini_;
			nKs_old_ = nKs_ini_;
			aur_old_ = aur_ini_;
			iur_old_ = iur_ini_;
			aKss_old_ = aKss_ini_;
			iKss_old_ = iKss_ini_;
			C_K2_old_ = C_K2_ini_;
			C_K1_old_ = C_K1_ini_;
			O_K_old_ = O_K_ini_;
			I_K_old_ = I_K_ini_;
		}
		while(1){
			flag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &time_new, CV_NORMAL);
			V_old_ = NV_Ith_S(dependent_variable__, 0);
			Cai_old_ = NV_Ith_S(dependent_variable__, 1);
			Cass_old_ = NV_Ith_S(dependent_variable__, 2);
			CaJSR_old_ = NV_Ith_S(dependent_variable__, 3);
			CaNSR_old_ = NV_Ith_S(dependent_variable__, 4);
			P_RyR_old_ = NV_Ith_S(dependent_variable__, 5);
			LTRPN_Ca_old_ = NV_Ith_S(dependent_variable__, 6);
			HTRPN_Ca_old_ = NV_Ith_S(dependent_variable__, 7);
			P_O1_old_ = NV_Ith_S(dependent_variable__, 8);
			P_O2_old_ = NV_Ith_S(dependent_variable__, 9);
			P_C2_old_ = NV_Ith_S(dependent_variable__, 10);
			O_old_ = NV_Ith_S(dependent_variable__, 11);
			C2_old_ = NV_Ith_S(dependent_variable__, 12);
			C3_old_ = NV_Ith_S(dependent_variable__, 13);
			C4_old_ = NV_Ith_S(dependent_variable__, 14);
			I1_old_ = NV_Ith_S(dependent_variable__, 15);
			I2_old_ = NV_Ith_S(dependent_variable__, 16);
			I3_old_ = NV_Ith_S(dependent_variable__, 17);
			Nai_old_ = NV_Ith_S(dependent_variable__, 18);
			C_Na2_old_ = NV_Ith_S(dependent_variable__, 19);
			C_Na1_old_ = NV_Ith_S(dependent_variable__, 20);
			O_Na_old_ = NV_Ith_S(dependent_variable__, 21);
			IF_Na_old_ = NV_Ith_S(dependent_variable__, 22);
			I1_Na_old_ = NV_Ith_S(dependent_variable__, 23);
			I2_Na_old_ = NV_Ith_S(dependent_variable__, 24);
			IC_Na2_old_ = NV_Ith_S(dependent_variable__, 25);
			IC_Na3_old_ = NV_Ith_S(dependent_variable__, 26);
			Ki_old_ = NV_Ith_S(dependent_variable__, 27);
			ato_f_old_ = NV_Ith_S(dependent_variable__, 28);
			ito_f_old_ = NV_Ith_S(dependent_variable__, 29);
			ato_s_old_ = NV_Ith_S(dependent_variable__, 30);
			ito_s_old_ = NV_Ith_S(dependent_variable__, 31);
			nKs_old_ = NV_Ith_S(dependent_variable__, 32);
			aur_old_ = NV_Ith_S(dependent_variable__, 33);
			iur_old_ = NV_Ith_S(dependent_variable__, 34);
			aKss_old_ = NV_Ith_S(dependent_variable__, 35);
			iKss_old_ = NV_Ith_S(dependent_variable__, 36);
			C_K2_old_ = NV_Ith_S(dependent_variable__, 37);
			C_K1_old_ = NV_Ith_S(dependent_variable__, 38);
			O_K_old_ = NV_Ith_S(dependent_variable__, 39);
			I_K_old_ = NV_Ith_S(dependent_variable__, 40);
			if (check_flag(&flag__, "CVode", 1)){printf("pau check_flag na %d iteracao\n", iout); break;}
			if (flag__ == CV_SUCCESS){
				if((iout != aux) && ((iout - aux)%offset_step == 0)) {
				CVodeGetLastStep(cvode_mem_cvode__, &cvodeDT);
					fprintf(f,"%.8e %.8e %.8e\n", time_new, V_old_,cvodeDT );
					counter_it__++;
				}
				iout++;
				tout += dtime; // timestep
				flag__ = CVodeGetLastStep(cvode_mem_cvode__, &dtL);
				dtM += dtL;
				if(dtL>dtMax) dtMax = dtL;
				if(dtL<dtMin) dtMin = dtL;
			}
			if (iout == num_iterations_bak ){ break;}
		}
		fclose(f);
		long int nsteps;
		long int fails;
		int flag = CVodeGetNumSteps(cvode_mem_cvode__ , &nsteps);
		flag = CVodeGetNumErrTestFails(cvode_mem_cvode__, &fails);
		
		
		printf("max: %e min: %e medias: %e %e %d fails %d\n", dtMax, dtMin, dtM/nsteps, dtM/iout , nsteps, fails);
	}
static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__){
	Solveode *ode = (Solveode *) f_data__;
	for(int i = 0; i<41; i++)
		ode->setVariables( i ,NV_Ith_S(dependent_variable__, i));
	ode->setParameters(0,time);
	rightHandSideFunction.function(ode);
	for(int i = 0; i<41; i++)
		NV_Ith_S(dep_var_dot__, i) = ode->getLadoDireito(i);
	return 0;
}



float __agos_factorial(int f){
	if(f>=0 & f<2)
		return 1.0;
	else if(f < 0)
		return 0.0/0.0;
	for(int i=f-1; i>=2; i--)
		f *= i;
	return (float)f;
}
double _agos_max(double* vector, int size){
	double max =vector[0];
	int i;
	for(i=1;i<size;i++){
		if(vector[i]>max) max = vector[i];
	}
	return max;
}
double _agos_min(double* vector, int size){
	double min = vector[0];
	int i;
	for(i=1;i<size;i++){
		if(vector[i]<min) min = vector[i];
	}
	return min;
}
double _agos_round( double x, int places )
{
	double const shift = powf( 10.0f, places );
	x *= shift;
	x = floorf( x + 0.5f );
	x /= shift;
	return x;
}
