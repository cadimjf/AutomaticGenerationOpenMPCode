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
#define forestSize 13
#include <sys/time.h>
#include "greedy.cpp"
#define _INCREASE_DT_ 1.5
#define _DECREASE_DT_ 0.65
#define _DECREASE_DT_2_ 2
#define numEDO 19
#define numAux 70

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
	double R; 	 // joule_per_mole_kelvin
	double T; 	 // kelvin
	double F; 	 // coulomb_per_millimole
	double Na_o; 	 // millimolar
	double K_o; 	 // millimolar
	double P_kna; 	 // dimensionless
	double Ca_o; 	 // millimolar
	double g_K1; 	 // nanoS_per_picoF
	double g_Kr; 	 // nanoS_per_picoF
	double g_Ks; 	 // nanoS_per_picoF
	double g_Na; 	 // nanoS_per_picoF
	double g_bna; 	 // nanoS_per_picoF
	double g_CaL; 	 // litre_per_farad_second
	double g_bca; 	 // nanoS_per_picoF
	double g_to; 	 // nanoS_per_picoF
	double P_NaK; 	 // picoA_per_picoF
	double K_mk; 	 // millimolar
	double K_mNa; 	 // millimolar
	double K_NaCa; 	 // picoA_per_picoF
	double gamma; 	 // dimensionless
	double alpha; 	 // dimensionless
	double Km_Nai; 	 // millimolar
	double Km_Ca; 	 // millimolar
	double K_sat; 	 // dimensionless
	double g_pCa; 	 // picoA_per_picoF
	double K_pCa; 	 // millimolar
	double g_pK; 	 // nanoS_per_picoF
	double V_rel; 	 // per_millisecond
	double Vmax_up; 	 // millimolar_per_millisecond
	double K_up; 	 // millimolar
	double V_leak; 	 // per_millisecond
	double V_xfer; 	 // per_millisecond
	double k3; 	 // per_millisecond
	double k4; 	 // per_millisecond
	double k1_prime; 	 // per_millimolar2_per_millisecond
	double k2_prime; 	 // per_millimolar_per_millisecond
	double max_sr; 	 // dimensionless
	double min_sr; 	 // dimensionless
	double EC; 	 // millimolar
	double Buf_c; 	 // millimolar
	double K_buf_c; 	 // millimolar
	double Buf_sr; 	 // millimolar
	double K_buf_sr; 	 // millimolar
	double Buf_ss; 	 // millimolar
	double K_buf_ss; 	 // millimolar
	double V_sr; 	 // micrometre3
	double V_c; 	 // micrometre3
	double Cm; 	 // microF
	double V_ss; 	 // micrometre3
	double calc_i_Stim; 	 // picoA_per_picoF
	double calc_E_Na; 	 // millivolt
	double calc_E_K; 	 // millivolt
	double calc_E_Ks; 	 // millivolt
	double calc_E_Ca; 	 // millivolt
	double calc_alpha_K1; 	 // dimensionless
	double calc_beta_K1; 	 // dimensionless
	double calc_xK1_inf; 	 // dimensionless
	double calc_i_K1; 	 // picoA_per_picoF
	double calc_i_Kr; 	 // picoA_per_picoF
	double calc_xr1_inf; 	 // dimensionless
	double calc_alpha_xr1; 	 // dimensionless
	double calc_beta_xr1; 	 // dimensionless
	double calc_tau_xr1; 	 // millisecond
	double calc_xr2_inf; 	 // dimensionless
	double calc_alpha_xr2; 	 // dimensionless
	double calc_beta_xr2; 	 // dimensionless
	double calc_tau_xr2; 	 // millisecond
	double calc_i_Ks; 	 // picoA_per_picoF
	double calc_xs_inf; 	 // dimensionless
	double calc_alpha_xs; 	 // dimensionless
	double calc_beta_xs; 	 // dimensionless
	double calc_tau_xs; 	 // millisecond
	double calc_i_Na; 	 // picoA_per_picoF
	double calc_m_inf; 	 // dimensionless
	double calc_alpha_m; 	 // dimensionless
	double calc_beta_m; 	 // dimensionless
	double calc_tau_m; 	 // millisecond
	double calc_h_inf; 	 // dimensionless
	double calc_alpha_h; 	 // per_millisecond
	double calc_beta_h; 	 // per_millisecond
	double calc_tau_h; 	 // millisecond
	double calc_j_inf; 	 // dimensionless
	double calc_alpha_j; 	 // per_millisecond
	double calc_beta_j; 	 // per_millisecond
	double calc_tau_j; 	 // millisecond
	double calc_i_b_Na; 	 // picoA_per_picoF
	double calc_i_CaL; 	 // picoA_per_picoF
	double calc_d_inf; 	 // dimensionless
	double calc_alpha_d; 	 // dimensionless
	double calc_beta_d; 	 // dimensionless
	double calc_gamma_d; 	 // millisecond
	double calc_tau_d; 	 // millisecond
	double calc_f_inf; 	 // dimensionless
	double calc_tau_f; 	 // millisecond
	double calc_f2_inf; 	 // dimensionless
	double calc_tau_f2; 	 // millisecond
	double calc_fCass_inf; 	 // dimensionless
	double calc_tau_fCass; 	 // millisecond
	double calc_i_b_Ca; 	 // picoA_per_picoF
	double calc_i_to; 	 // picoA_per_picoF
	double calc_s_inf; 	 // dimensionless
	double calc_tau_s; 	 // millisecond
	double calc_r_inf; 	 // dimensionless
	double calc_tau_r; 	 // millisecond
	double calc_i_NaK; 	 // picoA_per_picoF
	double calc_i_NaCa; 	 // picoA_per_picoF
	double calc_i_p_Ca; 	 // picoA_per_picoF
	double calc_i_p_K; 	 // picoA_per_picoF
	double calc_i_rel; 	 // millimolar_per_millisecond
	double calc_i_up; 	 // millimolar_per_millisecond
	double calc_i_leak; 	 // millimolar_per_millisecond
	double calc_i_xfer; 	 // millimolar_per_millisecond
	double calc_O; 	 // dimensionless
	double calc_k1; 	 // per_millimolar2_per_millisecond
	double calc_k2; 	 // per_millimolar_per_millisecond
	double calc_kcasr; 	 // dimensionless
	double calc_Ca_i_bufc; 	 // dimensionless
	double calc_Ca_sr_bufsr; 	 // dimensionless
	double calc_Ca_ss_bufss; 	 // dimensionless
	double dtime, *time_vec__;
	double time_new;

	//functions variables
	double *V;
	double V_new_, V_old_, V_ini_, V_lado_direito_;
	double *Xr1;
	double Xr1_new_, Xr1_old_, Xr1_ini_, Xr1_lado_direito_;
	double *Xr2;
	double Xr2_new_, Xr2_old_, Xr2_ini_, Xr2_lado_direito_;
	double *Xs;
	double Xs_new_, Xs_old_, Xs_ini_, Xs_lado_direito_;
	double *m;
	double m_new_, m_old_, m_ini_, m_lado_direito_;
	double *h;
	double h_new_, h_old_, h_ini_, h_lado_direito_;
	double *j;
	double j_new_, j_old_, j_ini_, j_lado_direito_;
	double *d;
	double d_new_, d_old_, d_ini_, d_lado_direito_;
	double *f;
	double f_new_, f_old_, f_ini_, f_lado_direito_;
	double *f2;
	double f2_new_, f2_old_, f2_ini_, f2_lado_direito_;
	double *fCass;
	double fCass_new_, fCass_old_, fCass_ini_, fCass_lado_direito_;
	double *s;
	double s_new_, s_old_, s_ini_, s_lado_direito_;
	double *r;
	double r_new_, r_old_, r_ini_, r_lado_direito_;
	double *R_prime;
	double R_prime_new_, R_prime_old_, R_prime_ini_, R_prime_lado_direito_;
	double *Ca_i;
	double Ca_i_new_, Ca_i_old_, Ca_i_ini_, Ca_i_lado_direito_;
	double *Ca_SR;
	double Ca_SR_new_, Ca_SR_old_, Ca_SR_ini_, Ca_SR_lado_direito_;
	double *Ca_ss;
	double Ca_ss_new_, Ca_ss_old_, Ca_ss_ini_, Ca_ss_lado_direito_;
	double *Na_i;
	double Na_i_new_, Na_i_old_, Na_i_ini_, Na_i_lado_direito_;
	double *K_i;
	double K_i_new_, K_i_old_, K_i_ini_, K_i_lado_direito_;

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
	inline double ifnumber_1();
	inline double ifnumber_2();
	inline double ifnumber_3();
	inline double ifnumber_4();
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
	__AGOS->calc_i_Stim = (((__AGOS->time_new>=__AGOS->stim_start)&&(__AGOS->time_new<=__AGOS->stim_end)&&(((__AGOS->time_new-__AGOS->stim_start)-(floor(((__AGOS->time_new-__AGOS->stim_start)/__AGOS->stim_period))*__AGOS->stim_period))<=__AGOS->stim_duration)))
?(__AGOS->stim_amplitude)
:(0.0000000000e+00);
	__AGOS->calc_E_Na = (((__AGOS->R*__AGOS->T)/__AGOS->F)*log((__AGOS->Na_o/__AGOS->Na_i_old_)));
	__AGOS->calc_E_K = (((__AGOS->R*__AGOS->T)/__AGOS->F)*log((__AGOS->K_o/__AGOS->K_i_old_)));
	__AGOS->calc_E_Ks = (((__AGOS->R*__AGOS->T)/__AGOS->F)*log(((__AGOS->K_o+(__AGOS->P_kna*__AGOS->Na_o))/(__AGOS->K_i_old_+(__AGOS->P_kna*__AGOS->Na_i_old_)))));
	__AGOS->calc_E_Ca = (((5.0000000000e-01*__AGOS->R*__AGOS->T)/__AGOS->F)*log((__AGOS->Ca_o/__AGOS->Ca_i_old_)));
	__AGOS->calc_i_CaL = ((((__AGOS->g_CaL*__AGOS->d_old_*__AGOS->f_old_*__AGOS->f2_old_*__AGOS->fCass_old_*4.0000000000e+00*(__AGOS->V_old_-1.5000000000e+01)*pow(__AGOS->F,2.0000000000e+00))/(__AGOS->R*__AGOS->T))*((2.5000000000e-01*__AGOS->Ca_ss_old_*exp(((2.0000000000e+00*(__AGOS->V_old_-1.5000000000e+01)*__AGOS->F)/(__AGOS->R*__AGOS->T))))-__AGOS->Ca_o))/(exp(((2.0000000000e+00*(__AGOS->V_old_-1.5000000000e+01)*__AGOS->F)/(__AGOS->R*__AGOS->T)))-1.0000000000e+00));
	__AGOS->calc_i_NaK = (((((__AGOS->P_NaK*__AGOS->K_o)/(__AGOS->K_o+__AGOS->K_mk))*__AGOS->Na_i_old_)/(__AGOS->Na_i_old_+__AGOS->K_mNa))/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T))))+(3.5300000000e-02*exp((((-__AGOS->V_old_)*__AGOS->F)/(__AGOS->R*__AGOS->T))))));
	__AGOS->calc_i_NaCa = ((__AGOS->K_NaCa*((exp(((__AGOS->gamma*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))*pow(__AGOS->Na_i_old_,3.0000000000e+00)*__AGOS->Ca_o)-(exp((((__AGOS->gamma-1.0000000000e+00)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))*pow(__AGOS->Na_o,3.0000000000e+00)*__AGOS->Ca_i_old_*__AGOS->alpha)))/((pow(__AGOS->Km_Nai,3.0000000000e+00)+pow(__AGOS->Na_o,3.0000000000e+00))*(__AGOS->Km_Ca+__AGOS->Ca_o)*(1.0000000000e+00+(__AGOS->K_sat*exp((((__AGOS->gamma-1.0000000000e+00)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))))));
	__AGOS->calc_i_p_Ca = ((__AGOS->g_pCa*__AGOS->Ca_i_old_)/(__AGOS->Ca_i_old_+__AGOS->K_pCa));
	__AGOS->calc_i_up = (__AGOS->Vmax_up/(1.0000000000e+00+(pow(__AGOS->K_up,2.0000000000e+00)/pow(__AGOS->Ca_i_old_,2.0000000000e+00))));
	__AGOS->calc_i_leak = (__AGOS->V_leak*(__AGOS->Ca_SR_old_-__AGOS->Ca_i_old_));
	__AGOS->calc_i_xfer = (__AGOS->V_xfer*(__AGOS->Ca_ss_old_-__AGOS->Ca_i_old_));
	__AGOS->calc_kcasr = (__AGOS->max_sr-((__AGOS->max_sr-__AGOS->min_sr)/(1.0000000000e+00+pow((__AGOS->EC/__AGOS->Ca_SR_old_),2.0000000000e+00))));
	__AGOS->calc_Ca_i_bufc = (1.0000000000e+00/(1.0000000000e+00+((__AGOS->Buf_c*__AGOS->K_buf_c)/pow((__AGOS->Ca_i_old_+__AGOS->K_buf_c),2.0000000000e+00))));
	__AGOS->calc_Ca_sr_bufsr = (1.0000000000e+00/(1.0000000000e+00+((__AGOS->Buf_sr*__AGOS->K_buf_sr)/pow((__AGOS->Ca_SR_old_+__AGOS->K_buf_sr),2.0000000000e+00))));
	__AGOS->calc_Ca_ss_bufss = (1.0000000000e+00/(1.0000000000e+00+((__AGOS->Buf_ss*__AGOS->K_buf_ss)/pow((__AGOS->Ca_ss_old_+__AGOS->K_buf_ss),2.0000000000e+00))));
	__AGOS->calc_alpha_K1 = (1.0000000000e-01/(1.0000000000e+00+exp((6.0000000000e-02*((__AGOS->V_old_-__AGOS->calc_E_K)-2.0000000000e+02)))));
	__AGOS->calc_beta_K1 = (((3.0000000000e+00*exp((2.0000000000e-04*((__AGOS->V_old_-__AGOS->calc_E_K)+1.0000000000e+02))))+exp((1.0000000000e-01*((__AGOS->V_old_-__AGOS->calc_E_K)-1.0000000000e+01))))/(1.0000000000e+00+exp(((-5.0000000000e-01)*(__AGOS->V_old_-__AGOS->calc_E_K)))));
	__AGOS->calc_i_Kr = (__AGOS->g_Kr*__AGOS->Xr1_old_*__AGOS->Xr2_old_*(__AGOS->V_old_-__AGOS->calc_E_K)*pow((__AGOS->K_o/5.4000000000e+00),1.0/2.0));
	__AGOS->calc_i_Ks = (__AGOS->g_Ks*pow(__AGOS->Xs_old_,2.0000000000e+00)*(__AGOS->V_old_-__AGOS->calc_E_Ks));
	__AGOS->calc_i_Na = (__AGOS->g_Na*pow(__AGOS->m_old_,3.0000000000e+00)*__AGOS->h_old_*__AGOS->j_old_*(__AGOS->V_old_-__AGOS->calc_E_Na));
	__AGOS->calc_i_b_Na = (__AGOS->g_bna*(__AGOS->V_old_-__AGOS->calc_E_Na));
	__AGOS->calc_i_b_Ca = (__AGOS->g_bca*(__AGOS->V_old_-__AGOS->calc_E_Ca));
	__AGOS->calc_i_to = (__AGOS->g_to*__AGOS->r_old_*__AGOS->s_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_p_K = ((__AGOS->g_pK*(__AGOS->V_old_-__AGOS->calc_E_K))/(1.0000000000e+00+exp(((2.5000000000e+01-__AGOS->V_old_)/5.9800000000e+00))));
	__AGOS->calc_xK1_inf = (__AGOS->calc_alpha_K1/(__AGOS->calc_alpha_K1+__AGOS->calc_beta_K1));
	__AGOS->calc_k1 = (__AGOS->k1_prime/__AGOS->calc_kcasr);
	__AGOS->calc_k2 = (__AGOS->k2_prime*__AGOS->calc_kcasr);
	__AGOS->calc_i_K1 = (__AGOS->g_K1*__AGOS->calc_xK1_inf*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_O = ((__AGOS->calc_k1*pow(__AGOS->Ca_ss_old_,2.0000000000e+00)*__AGOS->R_prime_old_)/(__AGOS->k3+(__AGOS->calc_k1*pow(__AGOS->Ca_ss_old_,2.0000000000e+00))));
	__AGOS->calc_i_rel = (__AGOS->V_rel*__AGOS->calc_O*(__AGOS->Ca_SR_old_-__AGOS->Ca_ss_old_));
	__AGOS->V_lado_direito_= (-(__AGOS->calc_i_K1+__AGOS->calc_i_to+__AGOS->calc_i_Kr+__AGOS->calc_i_Ks+__AGOS->calc_i_CaL+__AGOS->calc_i_NaK+__AGOS->calc_i_Na+__AGOS->calc_i_b_Na+__AGOS->calc_i_NaCa+__AGOS->calc_i_b_Ca+__AGOS->calc_i_p_K+__AGOS->calc_i_p_Ca+__AGOS->calc_i_Stim));
	__AGOS->R_prime_lado_direito_= (((-__AGOS->calc_k2)*__AGOS->Ca_ss_old_*__AGOS->R_prime_old_)+(__AGOS->k4*(1.0000000000e+00-__AGOS->R_prime_old_)));
	__AGOS->Ca_i_lado_direito_= (__AGOS->calc_Ca_i_bufc*(((((__AGOS->calc_i_leak-__AGOS->calc_i_up)*__AGOS->V_sr)/__AGOS->V_c)+__AGOS->calc_i_xfer)-((((__AGOS->calc_i_b_Ca+__AGOS->calc_i_p_Ca)-(2.0000000000e+00*__AGOS->calc_i_NaCa))*__AGOS->Cm)/(2.0000000000e+00*__AGOS->V_c*__AGOS->F))));
	__AGOS->Ca_SR_lado_direito_= (__AGOS->calc_Ca_sr_bufsr*(__AGOS->calc_i_up-(__AGOS->calc_i_rel+__AGOS->calc_i_leak)));
	__AGOS->Ca_ss_lado_direito_= (__AGOS->calc_Ca_ss_bufss*(((((-__AGOS->calc_i_CaL)*__AGOS->Cm)/(2.0000000000e+00*__AGOS->V_ss*__AGOS->F))+((__AGOS->calc_i_rel*__AGOS->V_sr)/__AGOS->V_ss))-((__AGOS->calc_i_xfer*__AGOS->V_c)/__AGOS->V_ss)));
	__AGOS->Na_i_lado_direito_= (((-(__AGOS->calc_i_Na+__AGOS->calc_i_b_Na+(3.0000000000e+00*__AGOS->calc_i_NaK)+(3.0000000000e+00*__AGOS->calc_i_NaCa)))/(__AGOS->V_c*__AGOS->F))*__AGOS->Cm);
	__AGOS->K_i_lado_direito_= (((-((__AGOS->calc_i_K1+__AGOS->calc_i_to+__AGOS->calc_i_Kr+__AGOS->calc_i_Ks+__AGOS->calc_i_p_K+__AGOS->calc_i_Stim)-(2.0000000000e+00*__AGOS->calc_i_NaK)))/(__AGOS->V_c*__AGOS->F))*__AGOS->Cm);
} //fim

void __tree2__( Solveode *__AGOS){
	__AGOS->calc_xr1_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-2.6000000000e+01)-__AGOS->V_old_)/7.0000000000e+00))));
	__AGOS->calc_alpha_xr1 = (4.5000000000e+02/(1.0000000000e+00+exp((((-4.5000000000e+01)-__AGOS->V_old_)/1.0000000000e+01))));
	__AGOS->calc_beta_xr1 = (6.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+3.0000000000e+01)/1.1500000000e+01))));
	__AGOS->calc_tau_xr1 = (1.0000000000e+00*__AGOS->calc_alpha_xr1*__AGOS->calc_beta_xr1);
	__AGOS->Xr1_lado_direito_= ((__AGOS->calc_xr1_inf-__AGOS->Xr1_old_)/__AGOS->calc_tau_xr1);
} //fim

void __tree3__( Solveode *__AGOS){
	__AGOS->calc_xr2_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+8.8000000000e+01)/2.4000000000e+01))));
	__AGOS->calc_alpha_xr2 = (3.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__AGOS->V_old_)/2.0000000000e+01))));
	__AGOS->calc_beta_xr2 = (1.1200000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_-6.0000000000e+01)/2.0000000000e+01))));
	__AGOS->calc_tau_xr2 = (1.0000000000e+00*__AGOS->calc_alpha_xr2*__AGOS->calc_beta_xr2);
	__AGOS->Xr2_lado_direito_= ((__AGOS->calc_xr2_inf-__AGOS->Xr2_old_)/__AGOS->calc_tau_xr2);
} //fim

void __tree4__( Solveode *__AGOS){
	__AGOS->calc_xs_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-5.0000000000e+00)-__AGOS->V_old_)/1.4000000000e+01))));
	__AGOS->calc_alpha_xs = (1.4000000000e+03/pow((1.0000000000e+00+exp(((5.0000000000e+00-__AGOS->V_old_)/6.0000000000e+00))),1.0/2.0));
	__AGOS->calc_beta_xs = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_-3.5000000000e+01)/1.5000000000e+01))));
	__AGOS->calc_tau_xs = ((1.0000000000e+00*__AGOS->calc_alpha_xs*__AGOS->calc_beta_xs)+8.0000000000e+01);
	__AGOS->Xs_lado_direito_= ((__AGOS->calc_xs_inf-__AGOS->Xs_old_)/__AGOS->calc_tau_xs);
} //fim

void __tree5__( Solveode *__AGOS){
	__AGOS->calc_m_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp((((-5.6860000000e+01)-__AGOS->V_old_)/9.0300000000e+00))),2.0000000000e+00));
	__AGOS->calc_alpha_m = (1.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__AGOS->V_old_)/5.0000000000e+00))));
	__AGOS->calc_beta_m = ((1.0000000000e-01/(1.0000000000e+00+exp(((__AGOS->V_old_+3.5000000000e+01)/5.0000000000e+00))))+(1.0000000000e-01/(1.0000000000e+00+exp(((__AGOS->V_old_-5.0000000000e+01)/2.0000000000e+02)))));
	__AGOS->calc_tau_m = (1.0000000000e+00*__AGOS->calc_alpha_m*__AGOS->calc_beta_m);
	__AGOS->m_lado_direito_= ((__AGOS->calc_m_inf-__AGOS->m_old_)/__AGOS->calc_tau_m);
} //fim

void __tree6__( Solveode *__AGOS){
	__AGOS->calc_h_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__AGOS->V_old_+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
	__AGOS->calc_alpha_h = ((__AGOS->V_old_<(-4.0000000000e+01)))
?((5.7000000000e-02*exp(((-(__AGOS->V_old_+8.0000000000e+01))/6.8000000000e+00))))
:(0.0000000000e+00);
	__AGOS->calc_beta_h = ((__AGOS->V_old_<(-4.0000000000e+01)))
?(((2.7000000000e+00*exp((7.9000000000e-02*__AGOS->V_old_)))+(3.1000000000e+05*exp((3.4850000000e-01*__AGOS->V_old_)))))
:((7.7000000000e-01/(1.3000000000e-01*(1.0000000000e+00+exp(((__AGOS->V_old_+1.0660000000e+01)/(-1.1100000000e+01)))))));
	__AGOS->calc_tau_h = (1.0000000000e+00/(__AGOS->calc_alpha_h+__AGOS->calc_beta_h));
	__AGOS->h_lado_direito_= ((__AGOS->calc_h_inf-__AGOS->h_old_)/__AGOS->calc_tau_h);
} //fim

void __tree7__( Solveode *__AGOS){
	__AGOS->calc_j_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__AGOS->V_old_+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
	__AGOS->calc_alpha_j = ((__AGOS->V_old_<(-4.0000000000e+01)))
?(((((((-2.5428000000e+04)*exp((2.4440000000e-01*__AGOS->V_old_)))-(6.9480000000e-06*exp(((-4.3910000000e-02)*__AGOS->V_old_))))*(__AGOS->V_old_+3.7780000000e+01))/1.0000000000e+00)/(1.0000000000e+00+exp((3.1100000000e-01*(__AGOS->V_old_+7.9230000000e+01))))))
:(0.0000000000e+00);
	__AGOS->calc_beta_j = ((__AGOS->V_old_<(-4.0000000000e+01)))
?(((2.4240000000e-02*exp(((-1.0520000000e-02)*__AGOS->V_old_)))/(1.0000000000e+00+exp(((-1.3780000000e-01)*(__AGOS->V_old_+4.0140000000e+01))))))
:(((6.0000000000e-01*exp((5.7000000000e-02*__AGOS->V_old_)))/(1.0000000000e+00+exp(((-1.0000000000e-01)*(__AGOS->V_old_+3.2000000000e+01))))));
	__AGOS->calc_tau_j = (1.0000000000e+00/(__AGOS->calc_alpha_j+__AGOS->calc_beta_j));
	__AGOS->j_lado_direito_= ((__AGOS->calc_j_inf-__AGOS->j_old_)/__AGOS->calc_tau_j);
} //fim

void __tree8__( Solveode *__AGOS){
	__AGOS->calc_d_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-8.0000000000e+00)-__AGOS->V_old_)/7.5000000000e+00))));
	__AGOS->calc_alpha_d = ((1.4000000000e+00/(1.0000000000e+00+exp((((-3.5000000000e+01)-__AGOS->V_old_)/1.3000000000e+01))))+2.5000000000e-01);
	__AGOS->calc_beta_d = (1.4000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+5.0000000000e+00)/5.0000000000e+00))));
	__AGOS->calc_gamma_d = (1.0000000000e+00/(1.0000000000e+00+exp(((5.0000000000e+01-__AGOS->V_old_)/2.0000000000e+01))));
	__AGOS->calc_tau_d = ((1.0000000000e+00*__AGOS->calc_alpha_d*__AGOS->calc_beta_d)+__AGOS->calc_gamma_d);
	__AGOS->d_lado_direito_= ((__AGOS->calc_d_inf-__AGOS->d_old_)/__AGOS->calc_tau_d);
} //fim

void __tree9__( Solveode *__AGOS){
	__AGOS->calc_f_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+2.0000000000e+01)/7.0000000000e+00))));
	__AGOS->calc_tau_f = ((1.1025000000e+03*exp(((-pow((__AGOS->V_old_+2.7000000000e+01),2.0000000000e+00))/2.2500000000e+02)))+(2.0000000000e+02/(1.0000000000e+00+exp(((1.3000000000e+01-__AGOS->V_old_)/1.0000000000e+01))))+(1.8000000000e+02/(1.0000000000e+00+exp(((__AGOS->V_old_+3.0000000000e+01)/1.0000000000e+01))))+2.0000000000e+01);
	__AGOS->f_lado_direito_= ((__AGOS->calc_f_inf-__AGOS->f_old_)/__AGOS->calc_tau_f);
} //fim

void __tree10__( Solveode *__AGOS){
	__AGOS->calc_f2_inf = ((6.7000000000e-01/(1.0000000000e+00+exp(((__AGOS->V_old_+3.5000000000e+01)/7.0000000000e+00))))+3.3000000000e-01);
	__AGOS->calc_tau_f2 = ((5.6200000000e+02*exp(((-pow((__AGOS->V_old_+2.7000000000e+01),2.0000000000e+00))/2.4000000000e+02)))+(3.1000000000e+01/(1.0000000000e+00+exp(((2.5000000000e+01-__AGOS->V_old_)/1.0000000000e+01))))+(8.0000000000e+01/(1.0000000000e+00+exp(((__AGOS->V_old_+3.0000000000e+01)/1.0000000000e+01)))));
	__AGOS->f2_lado_direito_= ((__AGOS->calc_f2_inf-__AGOS->f2_old_)/__AGOS->calc_tau_f2);
} //fim

void __tree11__( Solveode *__AGOS){
	__AGOS->calc_fCass_inf = ((6.0000000000e-01/(1.0000000000e+00+pow((__AGOS->Ca_ss_old_/5.0000000000e-02),2.0000000000e+00)))+4.0000000000e-01);
	__AGOS->calc_tau_fCass = ((8.0000000000e+01/(1.0000000000e+00+pow((__AGOS->Ca_ss_old_/5.0000000000e-02),2.0000000000e+00)))+2.0000000000e+00);
	__AGOS->fCass_lado_direito_= ((__AGOS->calc_fCass_inf-__AGOS->fCass_old_)/__AGOS->calc_tau_fCass);
} //fim

void __tree12__( Solveode *__AGOS){
	__AGOS->calc_s_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+2.0000000000e+01)/5.0000000000e+00))));
	__AGOS->calc_tau_s = ((8.5000000000e+01*exp(((-pow((__AGOS->V_old_+4.5000000000e+01),2.0000000000e+00))/3.2000000000e+02)))+(5.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_-2.0000000000e+01)/5.0000000000e+00))))+3.0000000000e+00);
	__AGOS->s_lado_direito_= ((__AGOS->calc_s_inf-__AGOS->s_old_)/__AGOS->calc_tau_s);
} //fim

void __tree13__( Solveode *__AGOS){
	__AGOS->calc_r_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((2.0000000000e+01-__AGOS->V_old_)/6.0000000000e+00))));
	__AGOS->calc_tau_r = ((9.5000000000e+00*exp(((-pow((__AGOS->V_old_+4.0000000000e+01),2.0000000000e+00))/1.8000000000e+03)))+8.0000000000e-01);
	__AGOS->r_lado_direito_= ((__AGOS->calc_r_inf-__AGOS->r_old_)/__AGOS->calc_tau_r);
} //fim

void __AGOS_EQUATIONS__( Solveode *__AGOS){
		const double time_new = __AGOS->time_new;
		const double time = 0.0000000000e+00;
		const double stim_amplitude = -5.2000000000e+01;
		const double stim_start = 1.0000000000e+02;
		const double stim_end = 1.0000000000e+05;
		const double stim_period = 5.0000000000e+02;
		const double stim_duration = 1.5000000000e+00;
		const double R = 8.3144720000e+03;
		const double T = 3.1000000000e+02;
		const double F = 9.6485341500e+04;
		const double Na_o = 1.4000000000e+02;
		const double K_o = 5.4000000000e+00;
		const double P_kna = 3.0000000000e-02;
		const double Ca_o = 2.0000000000e+00;
		const double g_K1 = 5.4050000000e+00;
		const double g_Kr = 1.5300000000e-01;
		const double g_Ks = 9.8000000000e-02;
		const double g_Na = 1.4838000000e+01;
		const double g_bna = 2.9000000000e-04;
		const double g_CaL = 3.9800000000e-05;
		const double g_bca = 5.9200000000e-04;
		const double g_to = 2.9400000000e-01;
		const double P_NaK = 2.7240000000e+00;
		const double K_mk = 1.0000000000e+00;
		const double K_mNa = 4.0000000000e+01;
		const double K_NaCa = 1.0000000000e+03;
		const double gamma = 3.5000000000e-01;
		const double alpha = 2.5000000000e+00;
		const double Km_Nai = 8.7500000000e+01;
		const double Km_Ca = 1.3800000000e+00;
		const double K_sat = 1.0000000000e-01;
		const double g_pCa = 1.2380000000e-01;
		const double K_pCa = 5.0000000000e-04;
		const double g_pK = 1.4600000000e-02;
		const double V_rel = 1.0200000000e-01;
		const double Vmax_up = 6.3750000000e-03;
		const double K_up = 2.5000000000e-04;
		const double V_leak = 3.6000000000e-04;
		const double V_xfer = 3.8000000000e-03;
		const double k3 = 6.0000000000e-02;
		const double k4 = 5.0000000000e-03;
		const double k1_prime = 1.5000000000e-01;
		const double k2_prime = 4.5000000000e-02;
		const double max_sr = 2.5000000000e+00;
		const double min_sr = 1.0000000000e+00;
		const double EC = 1.5000000000e+00;
		const double Buf_c = 2.0000000000e-01;
		const double K_buf_c = 1.0000000000e-03;
		const double Buf_sr = 1.0000000000e+01;
		const double K_buf_sr = 3.0000000000e-01;
		const double Buf_ss = 4.0000000000e-01;
		const double K_buf_ss = 2.5000000000e-04;
		const double V_sr = 1.0940000000e-03;
		const double V_c = 1.6404000000e-02;
		const double Cm = 1.8500000000e-01;
		const double V_ss = 5.4680000000e-05;
		const double V_old_= __AGOS->V_old_;
		const double Xr1_old_= __AGOS->Xr1_old_;
		const double Xr2_old_= __AGOS->Xr2_old_;
		const double Xs_old_= __AGOS->Xs_old_;
		const double m_old_= __AGOS->m_old_;
		const double h_old_= __AGOS->h_old_;
		const double j_old_= __AGOS->j_old_;
		const double d_old_= __AGOS->d_old_;
		const double f_old_= __AGOS->f_old_;
		const double f2_old_= __AGOS->f2_old_;
		const double fCass_old_= __AGOS->fCass_old_;
		const double s_old_= __AGOS->s_old_;
		const double r_old_= __AGOS->r_old_;
		const double R_prime_old_= __AGOS->R_prime_old_;
		const double Ca_i_old_= __AGOS->Ca_i_old_;
		const double Ca_SR_old_= __AGOS->Ca_SR_old_;
		const double Ca_ss_old_= __AGOS->Ca_ss_old_;
		const double Na_i_old_= __AGOS->Na_i_old_;
		const double K_i_old_= __AGOS->K_i_old_;
	const double calc_i_Stim = (((time_new>=stim_start)&&(time_new<=stim_end)&&(((time_new-stim_start)-(floor(((time_new-stim_start)/stim_period))*stim_period))<=stim_duration)))
?(stim_amplitude)
:(0.0000000000e+00);
	const double calc_E_Na = (((R*T)/F)*log((Na_o/Na_i_old_)));
	const double calc_E_K = (((R*T)/F)*log((K_o/K_i_old_)));
	const double calc_E_Ks = (((R*T)/F)*log(((K_o+(P_kna*Na_o))/(K_i_old_+(P_kna*Na_i_old_)))));
	const double calc_E_Ca = (((5.0000000000e-01*R*T)/F)*log((Ca_o/Ca_i_old_)));
	const double calc_i_CaL = ((((g_CaL*d_old_*f_old_*f2_old_*fCass_old_*4.0000000000e+00*(V_old_-1.5000000000e+01)*pow(F,2.0000000000e+00))/(R*T))*((2.5000000000e-01*Ca_ss_old_*exp(((2.0000000000e+00*(V_old_-1.5000000000e+01)*F)/(R*T))))-Ca_o))/(exp(((2.0000000000e+00*(V_old_-1.5000000000e+01)*F)/(R*T)))-1.0000000000e+00));
	const double calc_i_NaK = (((((P_NaK*K_o)/(K_o+K_mk))*Na_i_old_)/(Na_i_old_+K_mNa))/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*V_old_*F)/(R*T))))+(3.5300000000e-02*exp((((-V_old_)*F)/(R*T))))));
	const double calc_i_NaCa = ((K_NaCa*((exp(((gamma*V_old_*F)/(R*T)))*pow(Na_i_old_,3.0000000000e+00)*Ca_o)-(exp((((gamma-1.0000000000e+00)*V_old_*F)/(R*T)))*pow(Na_o,3.0000000000e+00)*Ca_i_old_*alpha)))/((pow(Km_Nai,3.0000000000e+00)+pow(Na_o,3.0000000000e+00))*(Km_Ca+Ca_o)*(1.0000000000e+00+(K_sat*exp((((gamma-1.0000000000e+00)*V_old_*F)/(R*T)))))));
	const double calc_i_p_Ca = ((g_pCa*Ca_i_old_)/(Ca_i_old_+K_pCa));
	const double calc_i_up = (Vmax_up/(1.0000000000e+00+(pow(K_up,2.0000000000e+00)/pow(Ca_i_old_,2.0000000000e+00))));
	const double calc_i_leak = (V_leak*(Ca_SR_old_-Ca_i_old_));
	const double calc_i_xfer = (V_xfer*(Ca_ss_old_-Ca_i_old_));
	const double calc_kcasr = (max_sr-((max_sr-min_sr)/(1.0000000000e+00+pow((EC/Ca_SR_old_),2.0000000000e+00))));
	const double calc_Ca_i_bufc = (1.0000000000e+00/(1.0000000000e+00+((Buf_c*K_buf_c)/pow((Ca_i_old_+K_buf_c),2.0000000000e+00))));
	const double calc_Ca_sr_bufsr = (1.0000000000e+00/(1.0000000000e+00+((Buf_sr*K_buf_sr)/pow((Ca_SR_old_+K_buf_sr),2.0000000000e+00))));
	const double calc_Ca_ss_bufss = (1.0000000000e+00/(1.0000000000e+00+((Buf_ss*K_buf_ss)/pow((Ca_ss_old_+K_buf_ss),2.0000000000e+00))));
	const double calc_alpha_K1 = (1.0000000000e-01/(1.0000000000e+00+exp((6.0000000000e-02*((V_old_-calc_E_K)-2.0000000000e+02)))));
	const double calc_beta_K1 = (((3.0000000000e+00*exp((2.0000000000e-04*((V_old_-calc_E_K)+1.0000000000e+02))))+exp((1.0000000000e-01*((V_old_-calc_E_K)-1.0000000000e+01))))/(1.0000000000e+00+exp(((-5.0000000000e-01)*(V_old_-calc_E_K)))));
	const double calc_i_Kr = (g_Kr*Xr1_old_*Xr2_old_*(V_old_-calc_E_K)*pow((K_o/5.4000000000e+00),1.0/2.0));
	const double calc_i_Ks = (g_Ks*pow(Xs_old_,2.0000000000e+00)*(V_old_-calc_E_Ks));
	const double calc_i_Na = (g_Na*pow(m_old_,3.0000000000e+00)*h_old_*j_old_*(V_old_-calc_E_Na));
	const double calc_i_b_Na = (g_bna*(V_old_-calc_E_Na));
	const double calc_i_b_Ca = (g_bca*(V_old_-calc_E_Ca));
	const double calc_i_to = (g_to*r_old_*s_old_*(V_old_-calc_E_K));
	const double calc_i_p_K = ((g_pK*(V_old_-calc_E_K))/(1.0000000000e+00+exp(((2.5000000000e+01-V_old_)/5.9800000000e+00))));
	const double calc_xK1_inf = (calc_alpha_K1/(calc_alpha_K1+calc_beta_K1));
	const double calc_k1 = (k1_prime/calc_kcasr);
	const double calc_k2 = (k2_prime*calc_kcasr);
	const double calc_i_K1 = (g_K1*calc_xK1_inf*(V_old_-calc_E_K));
	const double calc_O = ((calc_k1*pow(Ca_ss_old_,2.0000000000e+00)*R_prime_old_)/(k3+(calc_k1*pow(Ca_ss_old_,2.0000000000e+00))));
	const double calc_i_rel = (V_rel*calc_O*(Ca_SR_old_-Ca_ss_old_));
	__AGOS->V_lado_direito_= (-(calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_CaL+calc_i_NaK+calc_i_Na+calc_i_b_Na+calc_i_NaCa+calc_i_b_Ca+calc_i_p_K+calc_i_p_Ca+calc_i_Stim));
	__AGOS->R_prime_lado_direito_= (((-calc_k2)*Ca_ss_old_*R_prime_old_)+(k4*(1.0000000000e+00-R_prime_old_)));
	__AGOS->Ca_i_lado_direito_= (calc_Ca_i_bufc*(((((calc_i_leak-calc_i_up)*V_sr)/V_c)+calc_i_xfer)-((((calc_i_b_Ca+calc_i_p_Ca)-(2.0000000000e+00*calc_i_NaCa))*Cm)/(2.0000000000e+00*V_c*F))));
	__AGOS->Ca_SR_lado_direito_= (calc_Ca_sr_bufsr*(calc_i_up-(calc_i_rel+calc_i_leak)));
	__AGOS->Ca_ss_lado_direito_= (calc_Ca_ss_bufss*(((((-calc_i_CaL)*Cm)/(2.0000000000e+00*V_ss*F))+((calc_i_rel*V_sr)/V_ss))-((calc_i_xfer*V_c)/V_ss)));
	__AGOS->Na_i_lado_direito_= (((-(calc_i_Na+calc_i_b_Na+(3.0000000000e+00*calc_i_NaK)+(3.0000000000e+00*calc_i_NaCa)))/(V_c*F))*Cm);
	__AGOS->K_i_lado_direito_= (((-((calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_p_K+calc_i_Stim)-(2.0000000000e+00*calc_i_NaK)))/(V_c*F))*Cm);
	const double calc_xr1_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-2.6000000000e+01)-V_old_)/7.0000000000e+00))));
	const double calc_alpha_xr1 = (4.5000000000e+02/(1.0000000000e+00+exp((((-4.5000000000e+01)-V_old_)/1.0000000000e+01))));
	const double calc_beta_xr1 = (6.0000000000e+00/(1.0000000000e+00+exp(((V_old_+3.0000000000e+01)/1.1500000000e+01))));
	const double calc_tau_xr1 = (1.0000000000e+00*calc_alpha_xr1*calc_beta_xr1);
	__AGOS->Xr1_lado_direito_= ((calc_xr1_inf-Xr1_old_)/calc_tau_xr1);
	const double calc_xr2_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+8.8000000000e+01)/2.4000000000e+01))));
	const double calc_alpha_xr2 = (3.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-V_old_)/2.0000000000e+01))));
	const double calc_beta_xr2 = (1.1200000000e+00/(1.0000000000e+00+exp(((V_old_-6.0000000000e+01)/2.0000000000e+01))));
	const double calc_tau_xr2 = (1.0000000000e+00*calc_alpha_xr2*calc_beta_xr2);
	__AGOS->Xr2_lado_direito_= ((calc_xr2_inf-Xr2_old_)/calc_tau_xr2);
	const double calc_xs_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-5.0000000000e+00)-V_old_)/1.4000000000e+01))));
	const double calc_alpha_xs = (1.4000000000e+03/pow((1.0000000000e+00+exp(((5.0000000000e+00-V_old_)/6.0000000000e+00))),1.0/2.0));
	const double calc_beta_xs = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_-3.5000000000e+01)/1.5000000000e+01))));
	const double calc_tau_xs = ((1.0000000000e+00*calc_alpha_xs*calc_beta_xs)+8.0000000000e+01);
	__AGOS->Xs_lado_direito_= ((calc_xs_inf-Xs_old_)/calc_tau_xs);
	const double calc_m_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp((((-5.6860000000e+01)-V_old_)/9.0300000000e+00))),2.0000000000e+00));
	const double calc_alpha_m = (1.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-V_old_)/5.0000000000e+00))));
	const double calc_beta_m = ((1.0000000000e-01/(1.0000000000e+00+exp(((V_old_+3.5000000000e+01)/5.0000000000e+00))))+(1.0000000000e-01/(1.0000000000e+00+exp(((V_old_-5.0000000000e+01)/2.0000000000e+02)))));
	const double calc_tau_m = (1.0000000000e+00*calc_alpha_m*calc_beta_m);
	__AGOS->m_lado_direito_= ((calc_m_inf-m_old_)/calc_tau_m);
	const double calc_h_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((V_old_+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
	const double calc_alpha_h = ((V_old_<(-4.0000000000e+01)))
?((5.7000000000e-02*exp(((-(V_old_+8.0000000000e+01))/6.8000000000e+00))))
:(0.0000000000e+00);
	const double calc_beta_h = ((V_old_<(-4.0000000000e+01)))
?(((2.7000000000e+00*exp((7.9000000000e-02*V_old_)))+(3.1000000000e+05*exp((3.4850000000e-01*V_old_)))))
:((7.7000000000e-01/(1.3000000000e-01*(1.0000000000e+00+exp(((V_old_+1.0660000000e+01)/(-1.1100000000e+01)))))));
	const double calc_tau_h = (1.0000000000e+00/(calc_alpha_h+calc_beta_h));
	__AGOS->h_lado_direito_= ((calc_h_inf-h_old_)/calc_tau_h);
	const double calc_j_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((V_old_+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
	const double calc_alpha_j = ((V_old_<(-4.0000000000e+01)))
?(((((((-2.5428000000e+04)*exp((2.4440000000e-01*V_old_)))-(6.9480000000e-06*exp(((-4.3910000000e-02)*V_old_))))*(V_old_+3.7780000000e+01))/1.0000000000e+00)/(1.0000000000e+00+exp((3.1100000000e-01*(V_old_+7.9230000000e+01))))))
:(0.0000000000e+00);
	const double calc_beta_j = ((V_old_<(-4.0000000000e+01)))
?(((2.4240000000e-02*exp(((-1.0520000000e-02)*V_old_)))/(1.0000000000e+00+exp(((-1.3780000000e-01)*(V_old_+4.0140000000e+01))))))
:(((6.0000000000e-01*exp((5.7000000000e-02*V_old_)))/(1.0000000000e+00+exp(((-1.0000000000e-01)*(V_old_+3.2000000000e+01))))));
	const double calc_tau_j = (1.0000000000e+00/(calc_alpha_j+calc_beta_j));
	__AGOS->j_lado_direito_= ((calc_j_inf-j_old_)/calc_tau_j);
	const double calc_d_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-8.0000000000e+00)-V_old_)/7.5000000000e+00))));
	const double calc_alpha_d = ((1.4000000000e+00/(1.0000000000e+00+exp((((-3.5000000000e+01)-V_old_)/1.3000000000e+01))))+2.5000000000e-01);
	const double calc_beta_d = (1.4000000000e+00/(1.0000000000e+00+exp(((V_old_+5.0000000000e+00)/5.0000000000e+00))));
	const double calc_gamma_d = (1.0000000000e+00/(1.0000000000e+00+exp(((5.0000000000e+01-V_old_)/2.0000000000e+01))));
	const double calc_tau_d = ((1.0000000000e+00*calc_alpha_d*calc_beta_d)+calc_gamma_d);
	__AGOS->d_lado_direito_= ((calc_d_inf-d_old_)/calc_tau_d);
	const double calc_f_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+2.0000000000e+01)/7.0000000000e+00))));
	const double calc_tau_f = ((1.1025000000e+03*exp(((-pow((V_old_+2.7000000000e+01),2.0000000000e+00))/2.2500000000e+02)))+(2.0000000000e+02/(1.0000000000e+00+exp(((1.3000000000e+01-V_old_)/1.0000000000e+01))))+(1.8000000000e+02/(1.0000000000e+00+exp(((V_old_+3.0000000000e+01)/1.0000000000e+01))))+2.0000000000e+01);
	__AGOS->f_lado_direito_= ((calc_f_inf-f_old_)/calc_tau_f);
	const double calc_f2_inf = ((6.7000000000e-01/(1.0000000000e+00+exp(((V_old_+3.5000000000e+01)/7.0000000000e+00))))+3.3000000000e-01);
	const double calc_tau_f2 = ((5.6200000000e+02*exp(((-pow((V_old_+2.7000000000e+01),2.0000000000e+00))/2.4000000000e+02)))+(3.1000000000e+01/(1.0000000000e+00+exp(((2.5000000000e+01-V_old_)/1.0000000000e+01))))+(8.0000000000e+01/(1.0000000000e+00+exp(((V_old_+3.0000000000e+01)/1.0000000000e+01)))));
	__AGOS->f2_lado_direito_= ((calc_f2_inf-f2_old_)/calc_tau_f2);
	const double calc_fCass_inf = ((6.0000000000e-01/(1.0000000000e+00+pow((Ca_ss_old_/5.0000000000e-02),2.0000000000e+00)))+4.0000000000e-01);
	const double calc_tau_fCass = ((8.0000000000e+01/(1.0000000000e+00+pow((Ca_ss_old_/5.0000000000e-02),2.0000000000e+00)))+2.0000000000e+00);
	__AGOS->fCass_lado_direito_= ((calc_fCass_inf-fCass_old_)/calc_tau_fCass);
	const double calc_s_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+2.0000000000e+01)/5.0000000000e+00))));
	const double calc_tau_s = ((8.5000000000e+01*exp(((-pow((V_old_+4.5000000000e+01),2.0000000000e+00))/3.2000000000e+02)))+(5.0000000000e+00/(1.0000000000e+00+exp(((V_old_-2.0000000000e+01)/5.0000000000e+00))))+3.0000000000e+00);
	__AGOS->s_lado_direito_= ((calc_s_inf-s_old_)/calc_tau_s);
	const double calc_r_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((2.0000000000e+01-V_old_)/6.0000000000e+00))));
	const double calc_tau_r = ((9.5000000000e+00*exp(((-pow((V_old_+4.0000000000e+01),2.0000000000e+00))/1.8000000000e+03)))+8.0000000000e-01);
	__AGOS->r_lado_direito_= ((calc_r_inf-r_old_)/calc_tau_r);
} //fim
typedef struct str__rightHandSideFunction{
	void (*function)(Solveode*);
}typ_rightHandSideFunction;
typ_rightHandSideFunction rightHandSideFunction;
typ_rightHandSideFunction forest[13];
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
		time = 0.00000000e+00;
		stim_amplitude = -5.20000000e+01;
		stim_start = 1.00000000e+02;
		stim_end = 1.00000000e+05;
		stim_period = 5.00000000e+02;
		stim_duration = 1.50000000e+00;
		R = 8.31447200e+03;
		T = 3.10000000e+02;
		F = 9.64853415e+04;
		Na_o = 1.40000000e+02;
		K_o = 5.40000000e+00;
		P_kna = 3.00000000e-02;
		Ca_o = 2.00000000e+00;
		g_K1 = 5.40500000e+00;
		g_Kr = 1.53000000e-01;
		g_Ks = 9.80000000e-02;
		g_Na = 1.48380000e+01;
		g_bna = 2.90000000e-04;
		g_CaL = 3.98000000e-05;
		g_bca = 5.92000000e-04;
		g_to = 2.94000000e-01;
		P_NaK = 2.72400000e+00;
		K_mk = 1.00000000e+00;
		K_mNa = 4.00000000e+01;
		K_NaCa = 1.00000000e+03;
		gamma = 3.50000000e-01;
		alpha = 2.50000000e+00;
		Km_Nai = 8.75000000e+01;
		Km_Ca = 1.38000000e+00;
		K_sat = 1.00000000e-01;
		g_pCa = 1.23800000e-01;
		K_pCa = 5.00000000e-04;
		g_pK = 1.46000000e-02;
		V_rel = 1.02000000e-01;
		Vmax_up = 6.37500000e-03;
		K_up = 2.50000000e-04;
		V_leak = 3.60000000e-04;
		V_xfer = 3.80000000e-03;
		k3 = 6.00000000e-02;
		k4 = 5.00000000e-03;
		k1_prime = 1.50000000e-01;
		k2_prime = 4.50000000e-02;
		max_sr = 2.50000000e+00;
		min_sr = 1.00000000e+00;
		EC = 1.50000000e+00;
		Buf_c = 2.00000000e-01;
		K_buf_c = 1.00000000e-03;
		Buf_sr = 1.00000000e+01;
		K_buf_sr = 3.00000000e-01;
		Buf_ss = 4.00000000e-01;
		K_buf_ss = 2.50000000e-04;
		V_sr = 1.09400000e-03;
		V_c = 1.64040000e-02;
		Cm = 1.85000000e-01;
		V_ss = 5.46800000e-05;
		dtime = 0.0; time_vec__ = NULL;
		V = NULL;
		V_ini_ = -8.54230000e+01;
		Xr1 = NULL;
		Xr1_ini_ = 1.65000000e-02;
		Xr2 = NULL;
		Xr2_ini_ = 4.73000000e-01;
		Xs = NULL;
		Xs_ini_ = 1.74000000e-02;
		m = NULL;
		m_ini_ = 1.65000000e-03;
		h = NULL;
		h_ini_ = 7.49000000e-01;
		j = NULL;
		j_ini_ = 6.78800000e-01;
		d = NULL;
		d_ini_ = 3.28800000e-05;
		f = NULL;
		f_ini_ = 7.02600000e-01;
		f2 = NULL;
		f2_ini_ = 9.52600000e-01;
		fCass = NULL;
		fCass_ini_ = 9.94200000e-01;
		s = NULL;
		s_ini_ = 9.99998000e-01;
		r = NULL;
		r_ini_ = 2.34700000e-08;
		R_prime = NULL;
		R_prime_ini_ = 8.97800000e-01;
		Ca_i = NULL;
		Ca_i_ini_ = 1.53000000e-04;
		Ca_SR = NULL;
		Ca_SR_ini_ = 4.27200000e+00;
		Ca_ss = NULL;
		Ca_ss_ini_ = 4.20000000e-04;
		Na_i = NULL;
		Na_i_ini_ = 1.01320000e+01;
		K_i = NULL;
		K_i_ini_ = 1.38520000e+02;
		abstol__ = abs;
		reltol__ = rel;
		it_countx = 0;
	}
	Solveode::~Solveode()
	{
		if(V != NULL) free(V);
		if(Xr1 != NULL) free(Xr1);
		if(Xr2 != NULL) free(Xr2);
		if(Xs != NULL) free(Xs);
		if(m != NULL) free(m);
		if(h != NULL) free(h);
		if(j != NULL) free(j);
		if(d != NULL) free(d);
		if(f != NULL) free(f);
		if(f2 != NULL) free(f2);
		if(fCass != NULL) free(fCass);
		if(s != NULL) free(s);
		if(r != NULL) free(r);
		if(R_prime != NULL) free(R_prime);
		if(Ca_i != NULL) free(Ca_i);
		if(Ca_SR != NULL) free(Ca_SR);
		if(Ca_ss != NULL) free(Ca_ss);
		if(Na_i != NULL) free(Na_i);
		if(K_i != NULL) free(K_i);
	}

	int Solveode::setVariables(int indVariable, double value_new)
	{
		switch(indVariable)
		{
		case 0:		V_ini_ = V_old_= value_new;    break;
		case 1:		Xr1_ini_ = Xr1_old_= value_new;    break;
		case 2:		Xr2_ini_ = Xr2_old_= value_new;    break;
		case 3:		Xs_ini_ = Xs_old_= value_new;    break;
		case 4:		m_ini_ = m_old_= value_new;    break;
		case 5:		h_ini_ = h_old_= value_new;    break;
		case 6:		j_ini_ = j_old_= value_new;    break;
		case 7:		d_ini_ = d_old_= value_new;    break;
		case 8:		f_ini_ = f_old_= value_new;    break;
		case 9:		f2_ini_ = f2_old_= value_new;    break;
		case 10:		fCass_ini_ = fCass_old_= value_new;    break;
		case 11:		s_ini_ = s_old_= value_new;    break;
		case 12:		r_ini_ = r_old_= value_new;    break;
		case 13:		R_prime_ini_ = R_prime_old_= value_new;    break;
		case 14:		Ca_i_ini_ = Ca_i_old_= value_new;    break;
		case 15:		Ca_SR_ini_ = Ca_SR_old_= value_new;    break;
		case 16:		Ca_ss_ini_ = Ca_ss_old_= value_new;    break;
		case 17:		Na_i_ini_ = Na_i_old_= value_new;    break;
		case 18:		K_i_ini_ = K_i_old_= value_new;    break;
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
		case 6:		R = value_new;   break;
		case 7:		T = value_new;   break;
		case 8:		F = value_new;   break;
		case 9:		Na_o = value_new;   break;
		case 10:		K_o = value_new;   break;
		case 11:		P_kna = value_new;   break;
		case 12:		Ca_o = value_new;   break;
		case 13:		g_K1 = value_new;   break;
		case 14:		g_Kr = value_new;   break;
		case 15:		g_Ks = value_new;   break;
		case 16:		g_Na = value_new;   break;
		case 17:		g_bna = value_new;   break;
		case 18:		g_CaL = value_new;   break;
		case 19:		g_bca = value_new;   break;
		case 20:		g_to = value_new;   break;
		case 21:		P_NaK = value_new;   break;
		case 22:		K_mk = value_new;   break;
		case 23:		K_mNa = value_new;   break;
		case 24:		K_NaCa = value_new;   break;
		case 25:		gamma = value_new;   break;
		case 26:		alpha = value_new;   break;
		case 27:		Km_Nai = value_new;   break;
		case 28:		Km_Ca = value_new;   break;
		case 29:		K_sat = value_new;   break;
		case 30:		g_pCa = value_new;   break;
		case 31:		K_pCa = value_new;   break;
		case 32:		g_pK = value_new;   break;
		case 33:		V_rel = value_new;   break;
		case 34:		Vmax_up = value_new;   break;
		case 35:		K_up = value_new;   break;
		case 36:		V_leak = value_new;   break;
		case 37:		V_xfer = value_new;   break;
		case 38:		k3 = value_new;   break;
		case 39:		k4 = value_new;   break;
		case 40:		k1_prime = value_new;   break;
		case 41:		k2_prime = value_new;   break;
		case 42:		max_sr = value_new;   break;
		case 43:		min_sr = value_new;   break;
		case 44:		EC = value_new;   break;
		case 45:		Buf_c = value_new;   break;
		case 46:		K_buf_c = value_new;   break;
		case 47:		Buf_sr = value_new;   break;
		case 48:		K_buf_sr = value_new;   break;
		case 49:		Buf_ss = value_new;   break;
		case 50:		K_buf_ss = value_new;   break;
		case 51:		V_sr = value_new;   break;
		case 52:		V_c = value_new;   break;
		case 53:		Cm = value_new;   break;
		case 54:		V_ss = value_new;   break;
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
		
		system("ps aux | grep ./main ");
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
		system("ps aux | grep ./main ");
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
		this->Xr1_new_ = this->Xr1_old_ = this->Xr1_ini_;
		this->Xr2_new_ = this->Xr2_old_ = this->Xr2_ini_;
		this->Xs_new_ = this->Xs_old_ = this->Xs_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h_new_ = this->h_old_ = this->h_ini_;
		this->j_new_ = this->j_old_ = this->j_ini_;
		this->d_new_ = this->d_old_ = this->d_ini_;
		this->f_new_ = this->f_old_ = this->f_ini_;
		this->f2_new_ = this->f2_old_ = this->f2_ini_;
		this->fCass_new_ = this->fCass_old_ = this->fCass_ini_;
		this->s_new_ = this->s_old_ = this->s_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->R_prime_new_ = this->R_prime_old_ = this->R_prime_ini_;
		this->Ca_i_new_ = this->Ca_i_old_ = this->Ca_i_ini_;
		this->Ca_SR_new_ = this->Ca_SR_old_ = this->Ca_SR_ini_;
		this->Ca_ss_new_ = this->Ca_ss_old_ = this->Ca_ss_ini_;
		this->Na_i_new_ = this->Na_i_old_ = this->Na_i_ini_;
		this->K_i_new_ = this->K_i_old_ = this->K_i_ini_;
		this->time_new = this->time;
		this->timeSaving = this->time;
		if(savingRate!=0.0)
			this->save_step(fileptr, _EULER_);//save the initial conditions
		while(this->time_new<=finalTime){
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			this->V_new_ = this->dtime*(this->V_lado_direito_) + this->V_old_;
			this->Xr1_new_ = this->dtime*(this->Xr1_lado_direito_) + this->Xr1_old_;
			this->Xr2_new_ = this->dtime*(this->Xr2_lado_direito_) + this->Xr2_old_;
			this->Xs_new_ = this->dtime*(this->Xs_lado_direito_) + this->Xs_old_;
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h_new_ = this->dtime*(this->h_lado_direito_) + this->h_old_;
			this->j_new_ = this->dtime*(this->j_lado_direito_) + this->j_old_;
			this->d_new_ = this->dtime*(this->d_lado_direito_) + this->d_old_;
			this->f_new_ = this->dtime*(this->f_lado_direito_) + this->f_old_;
			this->f2_new_ = this->dtime*(this->f2_lado_direito_) + this->f2_old_;
			this->fCass_new_ = this->dtime*(this->fCass_lado_direito_) + this->fCass_old_;
			this->s_new_ = this->dtime*(this->s_lado_direito_) + this->s_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->R_prime_new_ = this->dtime*(this->R_prime_lado_direito_) + this->R_prime_old_;
			this->Ca_i_new_ = this->dtime*(this->Ca_i_lado_direito_) + this->Ca_i_old_;
			this->Ca_SR_new_ = this->dtime*(this->Ca_SR_lado_direito_) + this->Ca_SR_old_;
			this->Ca_ss_new_ = this->dtime*(this->Ca_ss_lado_direito_) + this->Ca_ss_old_;
			this->Na_i_new_ = this->dtime*(this->Na_i_lado_direito_) + this->Na_i_old_;
			this->K_i_new_ = this->dtime*(this->K_i_lado_direito_) + this->K_i_old_;
			//save results on a file
			if(savingRate!=0.0)
			{
				this->save_step(fileptr, _EULER_);
			}
		this->V_old_ = this->V_new_;
		this->Xr1_old_ = this->Xr1_new_;
		this->Xr2_old_ = this->Xr2_new_;
		this->Xs_old_ = this->Xs_new_;
		this->m_old_ = this->m_new_;
		this->h_old_ = this->h_new_;
		this->j_old_ = this->j_new_;
		this->d_old_ = this->d_new_;
		this->f_old_ = this->f_new_;
		this->f2_old_ = this->f2_new_;
		this->fCass_old_ = this->fCass_new_;
		this->s_old_ = this->s_new_;
		this->r_old_ = this->r_new_;
		this->R_prime_old_ = this->R_prime_new_;
		this->Ca_i_old_ = this->Ca_i_new_;
		this->Ca_SR_old_ = this->Ca_SR_new_;
		this->Ca_ss_old_ = this->Ca_ss_new_;
		this->Na_i_old_ = this->Na_i_new_;
		this->K_i_old_ = this->K_i_new_;
		}
	}
	void Solveode::rungeKutta2ndOrder(double finalTime, FILE *fileptr){
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->Xr1_new_ = this->Xr1_old_ = this->Xr1_ini_;
		this->Xr2_new_ = this->Xr2_old_ = this->Xr2_ini_;
		this->Xs_new_ = this->Xs_old_ = this->Xs_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h_new_ = this->h_old_ = this->h_ini_;
		this->j_new_ = this->j_old_ = this->j_ini_;
		this->d_new_ = this->d_old_ = this->d_ini_;
		this->f_new_ = this->f_old_ = this->f_ini_;
		this->f2_new_ = this->f2_old_ = this->f2_ini_;
		this->fCass_new_ = this->fCass_old_ = this->fCass_ini_;
		this->s_new_ = this->s_old_ = this->s_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->R_prime_new_ = this->R_prime_old_ = this->R_prime_ini_;
		this->Ca_i_new_ = this->Ca_i_old_ = this->Ca_i_ini_;
		this->Ca_SR_new_ = this->Ca_SR_old_ = this->Ca_SR_ini_;
		this->Ca_ss_new_ = this->Ca_ss_old_ = this->Ca_ss_ini_;
		this->Na_i_new_ = this->Na_i_old_ = this->Na_i_ini_;
		this->K_i_new_ = this->K_i_old_ = this->K_i_ini_;
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
			this->Xr1_new_ = this->dtime*(this->Xr1_lado_direito_) + this->Xr1_old_;
			this->Xr2_new_ = this->dtime*(this->Xr2_lado_direito_) + this->Xr2_old_;
			this->Xs_new_ = this->dtime*(this->Xs_lado_direito_) + this->Xs_old_;
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h_new_ = this->dtime*(this->h_lado_direito_) + this->h_old_;
			this->j_new_ = this->dtime*(this->j_lado_direito_) + this->j_old_;
			this->d_new_ = this->dtime*(this->d_lado_direito_) + this->d_old_;
			this->f_new_ = this->dtime*(this->f_lado_direito_) + this->f_old_;
			this->f2_new_ = this->dtime*(this->f2_lado_direito_) + this->f2_old_;
			this->fCass_new_ = this->dtime*(this->fCass_lado_direito_) + this->fCass_old_;
			this->s_new_ = this->dtime*(this->s_lado_direito_) + this->s_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->R_prime_new_ = this->dtime*(this->R_prime_lado_direito_) + this->R_prime_old_;
			this->Ca_i_new_ = this->dtime*(this->Ca_i_lado_direito_) + this->Ca_i_old_;
			this->Ca_SR_new_ = this->dtime*(this->Ca_SR_lado_direito_) + this->Ca_SR_old_;
			this->Ca_ss_new_ = this->dtime*(this->Ca_ss_lado_direito_) + this->Ca_ss_old_;
			this->Na_i_new_ = this->dtime*(this->Na_i_lado_direito_) + this->Na_i_old_;
			this->K_i_new_ = this->dtime*(this->K_i_lado_direito_) + this->K_i_old_;
			//stores the old variables in a vector
			for(int i=0;i<numEDO;i++){
				edos_old_aux_[i] = this->getVariables(i);
				edos_rightside_aux_[i] = this->getLadoDireito(i);
			}
			//steps one iteration ahead;
			this->V_old_ = this->V_new_;
			this->Xr1_old_ = this->Xr1_new_;
			this->Xr2_old_ = this->Xr2_new_;
			this->Xs_old_ = this->Xs_new_;
			this->m_old_ = this->m_new_;
			this->h_old_ = this->h_new_;
			this->j_old_ = this->j_new_;
			this->d_old_ = this->d_new_;
			this->f_old_ = this->f_new_;
			this->f2_old_ = this->f2_new_;
			this->fCass_old_ = this->fCass_new_;
			this->s_old_ = this->s_new_;
			this->r_old_ = this->r_new_;
			this->R_prime_old_ = this->R_prime_new_;
			this->Ca_i_old_ = this->Ca_i_new_;
			this->Ca_SR_old_ = this->Ca_SR_new_;
			this->Ca_ss_old_ = this->Ca_ss_new_;
			this->Na_i_old_ = this->Na_i_new_;
			this->K_i_old_ = this->K_i_new_;
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			//computes the runge kutta second order method
			this->V_new_ = ( this->V_lado_direito_ + edos_rightside_aux_[0] ) * this->dtime/2 + edos_old_aux_[0];
			this->Xr1_new_ = ( this->Xr1_lado_direito_ + edos_rightside_aux_[1] ) * this->dtime/2 + edos_old_aux_[1];
			this->Xr2_new_ = ( this->Xr2_lado_direito_ + edos_rightside_aux_[2] ) * this->dtime/2 + edos_old_aux_[2];
			this->Xs_new_ = ( this->Xs_lado_direito_ + edos_rightside_aux_[3] ) * this->dtime/2 + edos_old_aux_[3];
			this->m_new_ = ( this->m_lado_direito_ + edos_rightside_aux_[4] ) * this->dtime/2 + edos_old_aux_[4];
			this->h_new_ = ( this->h_lado_direito_ + edos_rightside_aux_[5] ) * this->dtime/2 + edos_old_aux_[5];
			this->j_new_ = ( this->j_lado_direito_ + edos_rightside_aux_[6] ) * this->dtime/2 + edos_old_aux_[6];
			this->d_new_ = ( this->d_lado_direito_ + edos_rightside_aux_[7] ) * this->dtime/2 + edos_old_aux_[7];
			this->f_new_ = ( this->f_lado_direito_ + edos_rightside_aux_[8] ) * this->dtime/2 + edos_old_aux_[8];
			this->f2_new_ = ( this->f2_lado_direito_ + edos_rightside_aux_[9] ) * this->dtime/2 + edos_old_aux_[9];
			this->fCass_new_ = ( this->fCass_lado_direito_ + edos_rightside_aux_[10] ) * this->dtime/2 + edos_old_aux_[10];
			this->s_new_ = ( this->s_lado_direito_ + edos_rightside_aux_[11] ) * this->dtime/2 + edos_old_aux_[11];
			this->r_new_ = ( this->r_lado_direito_ + edos_rightside_aux_[12] ) * this->dtime/2 + edos_old_aux_[12];
			this->R_prime_new_ = ( this->R_prime_lado_direito_ + edos_rightside_aux_[13] ) * this->dtime/2 + edos_old_aux_[13];
			this->Ca_i_new_ = ( this->Ca_i_lado_direito_ + edos_rightside_aux_[14] ) * this->dtime/2 + edos_old_aux_[14];
			this->Ca_SR_new_ = ( this->Ca_SR_lado_direito_ + edos_rightside_aux_[15] ) * this->dtime/2 + edos_old_aux_[15];
			this->Ca_ss_new_ = ( this->Ca_ss_lado_direito_ + edos_rightside_aux_[16] ) * this->dtime/2 + edos_old_aux_[16];
			this->Na_i_new_ = ( this->Na_i_lado_direito_ + edos_rightside_aux_[17] ) * this->dtime/2 + edos_old_aux_[17];
			this->K_i_new_ = ( this->K_i_lado_direito_ + edos_rightside_aux_[18] ) * this->dtime/2 + edos_old_aux_[18];
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
			this->Xr1_old_ = this->Xr1_new_;
			this->Xr2_old_ = this->Xr2_new_;
			this->Xs_old_ = this->Xs_new_;
			this->m_old_ = this->m_new_;
			this->h_old_ = this->h_new_;
			this->j_old_ = this->j_new_;
			this->d_old_ = this->d_new_;
			this->f_old_ = this->f_new_;
			this->f2_old_ = this->f2_new_;
			this->fCass_old_ = this->fCass_new_;
			this->s_old_ = this->s_new_;
			this->r_old_ = this->r_new_;
			this->R_prime_old_ = this->R_prime_new_;
			this->Ca_i_old_ = this->Ca_i_new_;
			this->Ca_SR_old_ = this->Ca_SR_new_;
			this->Ca_ss_old_ = this->Ca_ss_new_;
			this->Na_i_old_ = this->Na_i_new_;
			this->K_i_old_ = this->K_i_new_;
		}
	}
	void Solveode::addt(double finalTime, FILE *fileptr){
		double _tolerances_[numEDO];
		double _aux_tol=0.0;
		int desc;
		double maxDt = this->dtime, minDt = this->dtime;
		//initializes the variables
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->Xr1_new_ = this->Xr1_old_ = this->Xr1_ini_;
		this->Xr2_new_ = this->Xr2_old_ = this->Xr2_ini_;
		this->Xs_new_ = this->Xs_old_ = this->Xs_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h_new_ = this->h_old_ = this->h_ini_;
		this->j_new_ = this->j_old_ = this->j_ini_;
		this->d_new_ = this->d_old_ = this->d_ini_;
		this->f_new_ = this->f_old_ = this->f_ini_;
		this->f2_new_ = this->f2_old_ = this->f2_ini_;
		this->fCass_new_ = this->fCass_old_ = this->fCass_ini_;
		this->s_new_ = this->s_old_ = this->s_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->R_prime_new_ = this->R_prime_old_ = this->R_prime_ini_;
		this->Ca_i_new_ = this->Ca_i_old_ = this->Ca_i_ini_;
		this->Ca_SR_new_ = this->Ca_SR_old_ = this->Ca_SR_ini_;
		this->Ca_ss_new_ = this->Ca_ss_old_ = this->Ca_ss_ini_;
		this->Na_i_new_ = this->Na_i_old_ = this->Na_i_ini_;
		this->K_i_new_ = this->K_i_old_ = this->K_i_ini_;
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
		int _cont_=0;
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
						this->Xr1_new_ = edos_new_euler_[1];
						this->Xr2_new_ = edos_new_euler_[2];
						this->Xs_new_ = edos_new_euler_[3];
						this->m_new_ = edos_new_euler_[4];
						this->h_new_ = edos_new_euler_[5];
						this->j_new_ = edos_new_euler_[6];
						this->d_new_ = edos_new_euler_[7];
						this->f_new_ = edos_new_euler_[8];
						this->f2_new_ = edos_new_euler_[9];
						this->fCass_new_ = edos_new_euler_[10];
						this->s_new_ = edos_new_euler_[11];
						this->r_new_ = edos_new_euler_[12];
						this->R_prime_new_ = edos_new_euler_[13];
						this->Ca_i_new_ = edos_new_euler_[14];
						this->Ca_SR_new_ = edos_new_euler_[15];
						this->Ca_ss_new_ = edos_new_euler_[16];
						this->Na_i_new_ = edos_new_euler_[17];
						this->K_i_new_ = edos_new_euler_[18];
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
		printf("Dt max: %e dt min %e, %e %d\n", maxDt, minDt, _soma_/_cont_, _cont_);
		}
		printf("::desc %d\n", desc);
	}
	void Solveode::addt2(double finalTime, FILE *fileptr){
		const double _beta_safety_ = 0.8;
		double _tolerances_[numEDO];
		double _aux_tol=0.0;
		int desc=0;
		double maxDt = this->dtime, minDt = this->dtime;
		//initializes the variables
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->Xr1_new_ = this->Xr1_old_ = this->Xr1_ini_;
		this->Xr2_new_ = this->Xr2_old_ = this->Xr2_ini_;
		this->Xs_new_ = this->Xs_old_ = this->Xs_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h_new_ = this->h_old_ = this->h_ini_;
		this->j_new_ = this->j_old_ = this->j_ini_;
		this->d_new_ = this->d_old_ = this->d_ini_;
		this->f_new_ = this->f_old_ = this->f_ini_;
		this->f2_new_ = this->f2_old_ = this->f2_ini_;
		this->fCass_new_ = this->fCass_old_ = this->fCass_ini_;
		this->s_new_ = this->s_old_ = this->s_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->R_prime_new_ = this->R_prime_old_ = this->R_prime_ini_;
		this->Ca_i_new_ = this->Ca_i_old_ = this->Ca_i_ini_;
		this->Ca_SR_new_ = this->Ca_SR_old_ = this->Ca_SR_ini_;
		this->Ca_ss_new_ = this->Ca_ss_old_ = this->Ca_ss_ini_;
		this->Na_i_new_ = this->Na_i_old_ = this->Na_i_ini_;
		this->K_i_new_ = this->K_i_old_ = this->K_i_ini_;
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
		int _cont_=0;
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
				//restore the old values to do it again
				desc++;
				for(int i=0; i<numEDO; i++)
					this->setVariables(i, edos_old_aux_[i]);
				//throw the results away and compute again
				}else{//it accepts the solutions
					if(savingRate!=0.0){
					//restore the previous value of old variables
					for(int i=0; i<numEDO; i++)
						this->setVariables(i, edos_old_aux_[i]);
					this->V_new_ = edos_new_euler_[0];
					this->Xr1_new_ = edos_new_euler_[1];
					this->Xr2_new_ = edos_new_euler_[2];
					this->Xs_new_ = edos_new_euler_[3];
					this->m_new_ = edos_new_euler_[4];
					this->h_new_ = edos_new_euler_[5];
					this->j_new_ = edos_new_euler_[6];
					this->d_new_ = edos_new_euler_[7];
					this->f_new_ = edos_new_euler_[8];
					this->f2_new_ = edos_new_euler_[9];
					this->fCass_new_ = edos_new_euler_[10];
					this->s_new_ = edos_new_euler_[11];
					this->r_new_ = edos_new_euler_[12];
					this->R_prime_new_ = edos_new_euler_[13];
					this->Ca_i_new_ = edos_new_euler_[14];
					this->Ca_SR_new_ = edos_new_euler_[15];
					this->Ca_ss_new_ = edos_new_euler_[16];
					this->Na_i_new_ = edos_new_euler_[17];
					this->K_i_new_ = edos_new_euler_[18];
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
		if(savingRate!=0){
		printf("Dt max: %e dt min %e, %e  %d\n", maxDt, minDt, _soma_/_cont_ , _cont_);
		}
		printf("::desc %d\n", desc);
	}
	void Solveode::euler_OMP(double finalTime, FILE *fileptr, int nThreads){
		omp_set_num_threads(nThreads);
		double *__NEW_ = (double*)malloc(sizeof(double)*numEDO);
		double *__OLD_ = (double*)malloc(sizeof(double)*numEDO);
		double *temp;
		this->time_new = this->time;
		this->timeSaving = this->time;
		this->previous_dt = this->dtime;
		#pragma omp parallel firstprivate(__OLD_, __NEW_, temp)
		{
			//private parameters
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -5.2000000000e+01,  _prvt_stim_start = 1.0000000000e+02,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 5.0000000000e+02,  _prvt_stim_duration = 1.5000000000e+00,  _prvt_R = 8.3144720000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_F = 9.6485341500e+04,  _prvt_Na_o = 1.4000000000e+02,  _prvt_K_o = 5.4000000000e+00,  _prvt_P_kna = 3.0000000000e-02,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_g_K1 = 5.4050000000e+00,  _prvt_g_Kr = 1.5300000000e-01,  _prvt_g_Ks = 9.8000000000e-02,  _prvt_g_Na = 1.4838000000e+01,  _prvt_g_bna = 2.9000000000e-04,  _prvt_g_CaL = 3.9800000000e-05,  _prvt_g_bca = 5.9200000000e-04,  _prvt_g_to = 2.9400000000e-01,  _prvt_P_NaK = 2.7240000000e+00,  _prvt_K_mk = 1.0000000000e+00,  _prvt_K_mNa = 4.0000000000e+01,  _prvt_K_NaCa = 1.0000000000e+03,  _prvt_gamma = 3.5000000000e-01,  _prvt_alpha = 2.5000000000e+00,  _prvt_Km_Nai = 8.7500000000e+01,  _prvt_Km_Ca = 1.3800000000e+00,  _prvt_K_sat = 1.0000000000e-01,  _prvt_g_pCa = 1.2380000000e-01,  _prvt_K_pCa = 5.0000000000e-04,  _prvt_g_pK = 1.4600000000e-02,  _prvt_V_rel = 1.0200000000e-01,  _prvt_Vmax_up = 6.3750000000e-03,  _prvt_K_up = 2.5000000000e-04,  _prvt_V_leak = 3.6000000000e-04,  _prvt_V_xfer = 3.8000000000e-03,  _prvt_k3 = 6.0000000000e-02,  _prvt_k4 = 5.0000000000e-03,  _prvt_k1_prime = 1.5000000000e-01,  _prvt_k2_prime = 4.5000000000e-02,  _prvt_max_sr = 2.5000000000e+00,  _prvt_min_sr = 1.0000000000e+00,  _prvt_EC = 1.5000000000e+00,  _prvt_Buf_c = 2.0000000000e-01,  _prvt_K_buf_c = 1.0000000000e-03,  _prvt_Buf_sr = 1.0000000000e+01,  _prvt_K_buf_sr = 3.0000000000e-01,  _prvt_Buf_ss = 4.0000000000e-01,  _prvt_K_buf_ss = 2.5000000000e-04,  _prvt_V_sr = 1.0940000000e-03,  _prvt_V_c = 1.6404000000e-02,  _prvt_Cm = 1.8500000000e-01,  _prvt_V_ss = 5.4680000000e-05, 
			//private aux variables
			 _prvt_calc_i_Stim=0,  _prvt_calc_E_Na=0,  _prvt_calc_E_K=0,  _prvt_calc_E_Ks=0,  _prvt_calc_E_Ca=0,  _prvt_calc_alpha_K1=0,  _prvt_calc_beta_K1=0,  _prvt_calc_xK1_inf=0,  _prvt_calc_i_K1=0,  _prvt_calc_i_Kr=0,  _prvt_calc_xr1_inf=0,  _prvt_calc_alpha_xr1=0,  _prvt_calc_beta_xr1=0,  _prvt_calc_tau_xr1=0,  _prvt_calc_xr2_inf=0,  _prvt_calc_alpha_xr2=0,  _prvt_calc_beta_xr2=0,  _prvt_calc_tau_xr2=0,  _prvt_calc_i_Ks=0,  _prvt_calc_xs_inf=0,  _prvt_calc_alpha_xs=0,  _prvt_calc_beta_xs=0,  _prvt_calc_tau_xs=0,  _prvt_calc_i_Na=0,  _prvt_calc_m_inf=0,  _prvt_calc_alpha_m=0,  _prvt_calc_beta_m=0,  _prvt_calc_tau_m=0,  _prvt_calc_h_inf=0,  _prvt_calc_alpha_h=0,  _prvt_calc_beta_h=0,  _prvt_calc_tau_h=0,  _prvt_calc_j_inf=0,  _prvt_calc_alpha_j=0,  _prvt_calc_beta_j=0,  _prvt_calc_tau_j=0,  _prvt_calc_i_b_Na=0,  _prvt_calc_i_CaL=0,  _prvt_calc_d_inf=0,  _prvt_calc_alpha_d=0,  _prvt_calc_beta_d=0,  _prvt_calc_gamma_d=0,  _prvt_calc_tau_d=0,  _prvt_calc_f_inf=0,  _prvt_calc_tau_f=0,  _prvt_calc_f2_inf=0,  _prvt_calc_tau_f2=0,  _prvt_calc_fCass_inf=0,  _prvt_calc_tau_fCass=0,  _prvt_calc_i_b_Ca=0,  _prvt_calc_i_to=0,  _prvt_calc_s_inf=0,  _prvt_calc_tau_s=0,  _prvt_calc_r_inf=0,  _prvt_calc_tau_r=0,  _prvt_calc_i_NaK=0,  _prvt_calc_i_NaCa=0,  _prvt_calc_i_p_Ca=0,  _prvt_calc_i_p_K=0,  _prvt_calc_i_rel=0,  _prvt_calc_i_up=0,  _prvt_calc_i_leak=0,  _prvt_calc_i_xfer=0,  _prvt_calc_O=0,  _prvt_calc_k1=0,  _prvt_calc_k2=0,  _prvt_calc_kcasr=0,  _prvt_calc_Ca_i_bufc=0,  _prvt_calc_Ca_sr_bufsr=0,  _prvt_calc_Ca_ss_bufss=0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_Xr1_lado_direito_,  _prvt_Xr2_lado_direito_,  _prvt_Xs_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h_lado_direito_,  _prvt_j_lado_direito_,  _prvt_d_lado_direito_,  _prvt_f_lado_direito_,  _prvt_f2_lado_direito_,  _prvt_fCass_lado_direito_,  _prvt_s_lado_direito_,  _prvt_r_lado_direito_,  _prvt_R_prime_lado_direito_,  _prvt_Ca_i_lado_direito_,  _prvt_Ca_SR_lado_direito_,  _prvt_Ca_ss_lado_direito_,  _prvt_Na_i_lado_direito_,  _prvt_K_i_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate;
			__NEW_[0] = __OLD_[0] = -8.5423000000e+01;
			__NEW_[1] = __OLD_[1] = 1.6500000000e-02;
			__NEW_[2] = __OLD_[2] = 4.7300000000e-01;
			__NEW_[3] = __OLD_[3] = 1.7400000000e-02;
			__NEW_[4] = __OLD_[4] = 1.6500000000e-03;
			__NEW_[5] = __OLD_[5] = 7.4900000000e-01;
			__NEW_[6] = __OLD_[6] = 6.7880000000e-01;
			__NEW_[7] = __OLD_[7] = 3.2880000000e-05;
			__NEW_[8] = __OLD_[8] = 7.0260000000e-01;
			__NEW_[9] = __OLD_[9] = 9.5260000000e-01;
			__NEW_[10] = __OLD_[10] = 9.9420000000e-01;
			__NEW_[11] = __OLD_[11] = 9.9999800000e-01;
			__NEW_[12] = __OLD_[12] = 2.3470000000e-08;
			__NEW_[13] = __OLD_[13] = 8.9780000000e-01;
			__NEW_[14] = __OLD_[14] = 1.5300000000e-04;
			__NEW_[15] = __OLD_[15] = 4.2720000000e+00;
			__NEW_[16] = __OLD_[16] = 4.2000000000e-04;
			__NEW_[17] = __OLD_[17] = 1.0132000000e+01;
			__NEW_[18] = __OLD_[18] = 1.3852000000e+02;
			int *_prvt_tree_thread = tree_thread;
			while(_prvt_time_new<=_prvt_finalTime){
				_prvt_time_new += _prvt_dtime;
				if(omp_get_thread_num()==tree_thread[0])
				{
					_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[17])));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[18])));
					_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[18]+(_prvt_P_kna*__OLD_[17])))));
					_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[14])));
					_prvt_calc_i_CaL = ((((_prvt_g_CaL*__OLD_[7]*__OLD_[8]*__OLD_[9]*__OLD_[10]*4.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*((2.5000000000e-01*__OLD_[16]*exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))))-_prvt_Ca_o))/(exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00));
					_prvt_calc_i_NaK = (((((_prvt_P_NaK*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk))*__OLD_[17])/(__OLD_[17]+_prvt_K_mNa))/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.5300000000e-02*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_NaCa = ((_prvt_K_NaCa*((exp(((_prvt_gamma*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[17],3.0000000000e+00)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,3.0000000000e+00)*__OLD_[14]*_prvt_alpha)))/((pow(_prvt_Km_Nai,3.0000000000e+00)+pow(_prvt_Na_o,3.0000000000e+00))*(_prvt_Km_Ca+_prvt_Ca_o)*(1.0000000000e+00+(_prvt_K_sat*exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))))));
					_prvt_calc_i_p_Ca = ((_prvt_g_pCa*__OLD_[14])/(__OLD_[14]+_prvt_K_pCa));
					_prvt_calc_i_up = (_prvt_Vmax_up/(1.0000000000e+00+(pow(_prvt_K_up,2.0000000000e+00)/pow(__OLD_[14],2.0000000000e+00))));
					_prvt_calc_i_leak = (_prvt_V_leak*(__OLD_[15]-__OLD_[14]));
					_prvt_calc_i_xfer = (_prvt_V_xfer*(__OLD_[16]-__OLD_[14]));
					_prvt_calc_kcasr = (_prvt_max_sr-((_prvt_max_sr-_prvt_min_sr)/(1.0000000000e+00+pow((_prvt_EC/__OLD_[15]),2.0000000000e+00))));
					_prvt_calc_Ca_i_bufc = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_c*_prvt_K_buf_c)/pow((__OLD_[14]+_prvt_K_buf_c),2.0000000000e+00))));
					_prvt_calc_Ca_sr_bufsr = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_sr*_prvt_K_buf_sr)/pow((__OLD_[15]+_prvt_K_buf_sr),2.0000000000e+00))));
					_prvt_calc_Ca_ss_bufss = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_ss*_prvt_K_buf_ss)/pow((__OLD_[16]+_prvt_K_buf_ss),2.0000000000e+00))));
					_prvt_calc_alpha_K1 = (1.0000000000e-01/(1.0000000000e+00+exp((6.0000000000e-02*((__OLD_[0]-_prvt_calc_E_K)-2.0000000000e+02)))));
					_prvt_calc_beta_K1 = (((3.0000000000e+00*exp((2.0000000000e-04*((__OLD_[0]-_prvt_calc_E_K)+1.0000000000e+02))))+exp((1.0000000000e-01*((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01))))/(1.0000000000e+00+exp(((-5.0000000000e-01)*(__OLD_[0]-_prvt_calc_E_K)))));
					_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[1]*__OLD_[2]*(__OLD_[0]-_prvt_calc_E_K)*pow((_prvt_K_o/5.4000000000e+00),1.0/2.0));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
					_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*__OLD_[6]*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_to = (_prvt_g_to*__OLD_[12]*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_p_K = ((_prvt_g_pK*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/5.9800000000e+00))));
					_prvt_calc_xK1_inf = (_prvt_calc_alpha_K1/(_prvt_calc_alpha_K1+_prvt_calc_beta_K1));
					_prvt_calc_k1 = (_prvt_k1_prime/_prvt_calc_kcasr);
					_prvt_calc_k2 = (_prvt_k2_prime*_prvt_calc_kcasr);
					_prvt_calc_i_K1 = (_prvt_g_K1*_prvt_calc_xK1_inf*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_O = ((_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00)*__OLD_[13])/(_prvt_k3+(_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00))));
					_prvt_calc_i_rel = (_prvt_V_rel*_prvt_calc_O*(__OLD_[15]-__OLD_[16]));
					_prvt_V_lado_direito_= (-(_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_CaL+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_NaCa+_prvt_calc_i_b_Ca+_prvt_calc_i_p_K+_prvt_calc_i_p_Ca+_prvt_calc_i_Stim));
					__NEW_[0]= _prvt_V_lado_direito_ * _prvt_dtime + __OLD_[0];
					_prvt_R_prime_lado_direito_= (((-_prvt_calc_k2)*__OLD_[16]*__OLD_[13])+(_prvt_k4*(1.0000000000e+00-__OLD_[13])));
					__NEW_[13]= _prvt_R_prime_lado_direito_ * _prvt_dtime + __OLD_[13];
					_prvt_Ca_i_lado_direito_= (_prvt_calc_Ca_i_bufc*(((((_prvt_calc_i_leak-_prvt_calc_i_up)*_prvt_V_sr)/_prvt_V_c)+_prvt_calc_i_xfer)-((((_prvt_calc_i_b_Ca+_prvt_calc_i_p_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Cm)/(2.0000000000e+00*_prvt_V_c*_prvt_F))));
					__NEW_[14]= _prvt_Ca_i_lado_direito_ * _prvt_dtime + __OLD_[14];
					_prvt_Ca_SR_lado_direito_= (_prvt_calc_Ca_sr_bufsr*(_prvt_calc_i_up-(_prvt_calc_i_rel+_prvt_calc_i_leak)));
					__NEW_[15]= _prvt_Ca_SR_lado_direito_ * _prvt_dtime + __OLD_[15];
					_prvt_Ca_ss_lado_direito_= (_prvt_calc_Ca_ss_bufss*(((((-_prvt_calc_i_CaL)*_prvt_Cm)/(2.0000000000e+00*_prvt_V_ss*_prvt_F))+((_prvt_calc_i_rel*_prvt_V_sr)/_prvt_V_ss))-((_prvt_calc_i_xfer*_prvt_V_c)/_prvt_V_ss)));
					__NEW_[16]= _prvt_Ca_ss_lado_direito_ * _prvt_dtime + __OLD_[16];
					_prvt_Na_i_lado_direito_= (((-(_prvt_calc_i_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
					__NEW_[17]= _prvt_Na_i_lado_direito_ * _prvt_dtime + __OLD_[17];
					_prvt_K_i_lado_direito_= (((-((_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_p_K+_prvt_calc_i_Stim)-(2.0000000000e+00*_prvt_calc_i_NaK)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
					__NEW_[18]= _prvt_K_i_lado_direito_ * _prvt_dtime + __OLD_[18];
				}
				if(omp_get_thread_num()==tree_thread[1])
				{
					_prvt_calc_xr1_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-2.6000000000e+01)-__OLD_[0])/7.0000000000e+00))));
					_prvt_calc_alpha_xr1 = (4.5000000000e+02/(1.0000000000e+00+exp((((-4.5000000000e+01)-__OLD_[0])/1.0000000000e+01))));
					_prvt_calc_beta_xr1 = (6.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.1500000000e+01))));
					_prvt_calc_tau_xr1 = (1.0000000000e+00*_prvt_calc_alpha_xr1*_prvt_calc_beta_xr1);
					_prvt_Xr1_lado_direito_= ((_prvt_calc_xr1_inf-__OLD_[1])/_prvt_calc_tau_xr1);
					__NEW_[1]= _prvt_Xr1_lado_direito_ * _prvt_dtime + __OLD_[1];
				}
				if(omp_get_thread_num()==tree_thread[2])
				{
					_prvt_calc_xr2_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+8.8000000000e+01)/2.4000000000e+01))));
					_prvt_calc_alpha_xr2 = (3.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/2.0000000000e+01))));
					_prvt_calc_beta_xr2 = (1.1200000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-6.0000000000e+01)/2.0000000000e+01))));
					_prvt_calc_tau_xr2 = (1.0000000000e+00*_prvt_calc_alpha_xr2*_prvt_calc_beta_xr2);
					_prvt_Xr2_lado_direito_= ((_prvt_calc_xr2_inf-__OLD_[2])/_prvt_calc_tau_xr2);
					__NEW_[2]= _prvt_Xr2_lado_direito_ * _prvt_dtime + __OLD_[2];
				}
				if(omp_get_thread_num()==tree_thread[3])
				{
					_prvt_calc_xs_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-5.0000000000e+00)-__OLD_[0])/1.4000000000e+01))));
					_prvt_calc_alpha_xs = (1.4000000000e+03/pow((1.0000000000e+00+exp(((5.0000000000e+00-__OLD_[0])/6.0000000000e+00))),1.0/2.0));
					_prvt_calc_beta_xs = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-3.5000000000e+01)/1.5000000000e+01))));
					_prvt_calc_tau_xs = ((1.0000000000e+00*_prvt_calc_alpha_xs*_prvt_calc_beta_xs)+8.0000000000e+01);
					_prvt_Xs_lado_direito_= ((_prvt_calc_xs_inf-__OLD_[3])/_prvt_calc_tau_xs);
					__NEW_[3]= _prvt_Xs_lado_direito_ * _prvt_dtime + __OLD_[3];
				}
				if(omp_get_thread_num()==tree_thread[4])
				{
					_prvt_calc_m_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp((((-5.6860000000e+01)-__OLD_[0])/9.0300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_m = (1.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/5.0000000000e+00))));
					_prvt_calc_beta_m = ((1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/5.0000000000e+00))))+(1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]-5.0000000000e+01)/2.0000000000e+02)))));
					_prvt_calc_tau_m = (1.0000000000e+00*_prvt_calc_alpha_m*_prvt_calc_beta_m);
					_prvt_m_lado_direito_= ((_prvt_calc_m_inf-__OLD_[4])/_prvt_calc_tau_m);
					__NEW_[4]= _prvt_m_lado_direito_ * _prvt_dtime + __OLD_[4];
				}
				if(omp_get_thread_num()==tree_thread[5])
				{
					_prvt_calc_h_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_h = ((__OLD_[0]<(-4.0000000000e+01)))
?((5.7000000000e-02*exp(((-(__OLD_[0]+8.0000000000e+01))/6.8000000000e+00))))
:(0.0000000000e+00);
					_prvt_calc_beta_h = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.7000000000e+00*exp((7.9000000000e-02*__OLD_[0])))+(3.1000000000e+05*exp((3.4850000000e-01*__OLD_[0])))))
:((7.7000000000e-01/(1.3000000000e-01*(1.0000000000e+00+exp(((__OLD_[0]+1.0660000000e+01)/(-1.1100000000e+01)))))));
					_prvt_calc_tau_h = (1.0000000000e+00/(_prvt_calc_alpha_h+_prvt_calc_beta_h));
					_prvt_h_lado_direito_= ((_prvt_calc_h_inf-__OLD_[5])/_prvt_calc_tau_h);
					__NEW_[5]= _prvt_h_lado_direito_ * _prvt_dtime + __OLD_[5];
				}
				if(omp_get_thread_num()==tree_thread[6])
				{
					_prvt_calc_j_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((((((-2.5428000000e+04)*exp((2.4440000000e-01*__OLD_[0])))-(6.9480000000e-06*exp(((-4.3910000000e-02)*__OLD_[0]))))*(__OLD_[0]+3.7780000000e+01))/1.0000000000e+00)/(1.0000000000e+00+exp((3.1100000000e-01*(__OLD_[0]+7.9230000000e+01))))))
:(0.0000000000e+00);
					_prvt_calc_beta_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.4240000000e-02*exp(((-1.0520000000e-02)*__OLD_[0])))/(1.0000000000e+00+exp(((-1.3780000000e-01)*(__OLD_[0]+4.0140000000e+01))))))
:(((6.0000000000e-01*exp((5.7000000000e-02*__OLD_[0])))/(1.0000000000e+00+exp(((-1.0000000000e-01)*(__OLD_[0]+3.2000000000e+01))))));
					_prvt_calc_tau_j = (1.0000000000e+00/(_prvt_calc_alpha_j+_prvt_calc_beta_j));
					_prvt_j_lado_direito_= ((_prvt_calc_j_inf-__OLD_[6])/_prvt_calc_tau_j);
					__NEW_[6]= _prvt_j_lado_direito_ * _prvt_dtime + __OLD_[6];
				}
				if(omp_get_thread_num()==tree_thread[7])
				{
					_prvt_calc_d_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-8.0000000000e+00)-__OLD_[0])/7.5000000000e+00))));
					_prvt_calc_alpha_d = ((1.4000000000e+00/(1.0000000000e+00+exp((((-3.5000000000e+01)-__OLD_[0])/1.3000000000e+01))))+2.5000000000e-01);
					_prvt_calc_beta_d = (1.4000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.0000000000e+00)/5.0000000000e+00))));
					_prvt_calc_gamma_d = (1.0000000000e+00/(1.0000000000e+00+exp(((5.0000000000e+01-__OLD_[0])/2.0000000000e+01))));
					_prvt_calc_tau_d = ((1.0000000000e+00*_prvt_calc_alpha_d*_prvt_calc_beta_d)+_prvt_calc_gamma_d);
					_prvt_d_lado_direito_= ((_prvt_calc_d_inf-__OLD_[7])/_prvt_calc_tau_d);
					__NEW_[7]= _prvt_d_lado_direito_ * _prvt_dtime + __OLD_[7];
				}
				if(omp_get_thread_num()==tree_thread[8])
				{
					_prvt_calc_f_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/7.0000000000e+00))));
					_prvt_calc_tau_f = ((1.1025000000e+03*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.2500000000e+02)))+(2.0000000000e+02/(1.0000000000e+00+exp(((1.3000000000e+01-__OLD_[0])/1.0000000000e+01))))+(1.8000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01))))+2.0000000000e+01);
					_prvt_f_lado_direito_= ((_prvt_calc_f_inf-__OLD_[8])/_prvt_calc_tau_f);
					__NEW_[8]= _prvt_f_lado_direito_ * _prvt_dtime + __OLD_[8];
				}
				if(omp_get_thread_num()==tree_thread[9])
				{
					_prvt_calc_f2_inf = ((6.7000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/7.0000000000e+00))))+3.3000000000e-01);
					_prvt_calc_tau_f2 = ((5.6200000000e+02*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.4000000000e+02)))+(3.1000000000e+01/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/1.0000000000e+01))))+(8.0000000000e+01/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01)))));
					_prvt_f2_lado_direito_= ((_prvt_calc_f2_inf-__OLD_[9])/_prvt_calc_tau_f2);
					__NEW_[9]= _prvt_f2_lado_direito_ * _prvt_dtime + __OLD_[9];
				}
				if(omp_get_thread_num()==tree_thread[10])
				{
					_prvt_calc_fCass_inf = ((6.0000000000e-01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+4.0000000000e-01);
					_prvt_calc_tau_fCass = ((8.0000000000e+01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+2.0000000000e+00);
					_prvt_fCass_lado_direito_= ((_prvt_calc_fCass_inf-__OLD_[10])/_prvt_calc_tau_fCass);
					__NEW_[10]= _prvt_fCass_lado_direito_ * _prvt_dtime + __OLD_[10];
				}
				if(omp_get_thread_num()==tree_thread[11])
				{
					_prvt_calc_s_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/5.0000000000e+00))));
					_prvt_calc_tau_s = ((8.5000000000e+01*exp(((-pow((__OLD_[0]+4.5000000000e+01),2.0000000000e+00))/3.2000000000e+02)))+(5.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-2.0000000000e+01)/5.0000000000e+00))))+3.0000000000e+00);
					_prvt_s_lado_direito_= ((_prvt_calc_s_inf-__OLD_[11])/_prvt_calc_tau_s);
					__NEW_[11]= _prvt_s_lado_direito_ * _prvt_dtime + __OLD_[11];
				}
				if(omp_get_thread_num()==tree_thread[12])
				{
					_prvt_calc_r_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((2.0000000000e+01-__OLD_[0])/6.0000000000e+00))));
					_prvt_calc_tau_r = ((9.5000000000e+00*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.8000000000e+03)))+8.0000000000e-01);
					_prvt_r_lado_direito_= ((_prvt_calc_r_inf-__OLD_[12])/_prvt_calc_tau_r);
					__NEW_[12]= _prvt_r_lado_direito_ * _prvt_dtime + __OLD_[12];
				}
				//synchronizing all threads
				#pragma omp barrier
				/*if(_prvt_savingRate!=0){
					#pragma omp single
					{
						this->V_old_ = __OLD_[0];
						this->V_new_ = __NEW_[0];
						this->Xr1_old_ = __OLD_[1];
						this->Xr1_new_ = __NEW_[1];
						this->Xr2_old_ = __OLD_[2];
						this->Xr2_new_ = __NEW_[2];
						this->Xs_old_ = __OLD_[3];
						this->Xs_new_ = __NEW_[3];
						this->m_old_ = __OLD_[4];
						this->m_new_ = __NEW_[4];
						this->h_old_ = __OLD_[5];
						this->h_new_ = __NEW_[5];
						this->j_old_ = __OLD_[6];
						this->j_new_ = __NEW_[6];
						this->d_old_ = __OLD_[7];
						this->d_new_ = __NEW_[7];
						this->f_old_ = __OLD_[8];
						this->f_new_ = __NEW_[8];
						this->f2_old_ = __OLD_[9];
						this->f2_new_ = __NEW_[9];
						this->fCass_old_ = __OLD_[10];
						this->fCass_new_ = __NEW_[10];
						this->s_old_ = __OLD_[11];
						this->s_new_ = __NEW_[11];
						this->r_old_ = __OLD_[12];
						this->r_new_ = __NEW_[12];
						this->R_prime_old_ = __OLD_[13];
						this->R_prime_new_ = __NEW_[13];
						this->Ca_i_old_ = __OLD_[14];
						this->Ca_i_new_ = __NEW_[14];
						this->Ca_SR_old_ = __OLD_[15];
						this->Ca_SR_new_ = __NEW_[15];
						this->Ca_ss_old_ = __OLD_[16];
						this->Ca_ss_new_ = __NEW_[16];
						this->Na_i_old_ = __OLD_[17];
						this->Na_i_new_ = __NEW_[17];
						this->K_i_old_ = __OLD_[18];
						this->K_i_new_ = __NEW_[18];
						this->time_new = _prvt_time_new;
						save_step(fileptr, _EULER_);
					}
				}*/
				temp = __OLD_;
				__OLD_ = __NEW_;
				__NEW_= temp;
			}
		}
	}
	void Solveode::addt_OMP(double finalTime, FILE *fileptr, int nThreads){
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
			int *_prvt_tree_thread = tree_thread;
			double _prvt_rel_tol_=reltol__, _prvt_abs_tol_ =abstol__, _prvt_aux_tol;
			const double __tiny_ = pow(_prvt_abs_tol_, 2.0);
			//private parameters
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -5.2000000000e+01,  _prvt_stim_start = 1.0000000000e+02,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 5.0000000000e+02,  _prvt_stim_duration = 1.5000000000e+00,  _prvt_R = 8.3144720000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_F = 9.6485341500e+04,  _prvt_Na_o = 1.4000000000e+02,  _prvt_K_o = 5.4000000000e+00,  _prvt_P_kna = 3.0000000000e-02,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_g_K1 = 5.4050000000e+00,  _prvt_g_Kr = 1.5300000000e-01,  _prvt_g_Ks = 9.8000000000e-02,  _prvt_g_Na = 1.4838000000e+01,  _prvt_g_bna = 2.9000000000e-04,  _prvt_g_CaL = 3.9800000000e-05,  _prvt_g_bca = 5.9200000000e-04,  _prvt_g_to = 2.9400000000e-01,  _prvt_P_NaK = 2.7240000000e+00,  _prvt_K_mk = 1.0000000000e+00,  _prvt_K_mNa = 4.0000000000e+01,  _prvt_K_NaCa = 1.0000000000e+03,  _prvt_gamma = 3.5000000000e-01,  _prvt_alpha = 2.5000000000e+00,  _prvt_Km_Nai = 8.7500000000e+01,  _prvt_Km_Ca = 1.3800000000e+00,  _prvt_K_sat = 1.0000000000e-01,  _prvt_g_pCa = 1.2380000000e-01,  _prvt_K_pCa = 5.0000000000e-04,  _prvt_g_pK = 1.4600000000e-02,  _prvt_V_rel = 1.0200000000e-01,  _prvt_Vmax_up = 6.3750000000e-03,  _prvt_K_up = 2.5000000000e-04,  _prvt_V_leak = 3.6000000000e-04,  _prvt_V_xfer = 3.8000000000e-03,  _prvt_k3 = 6.0000000000e-02,  _prvt_k4 = 5.0000000000e-03,  _prvt_k1_prime = 1.5000000000e-01,  _prvt_k2_prime = 4.5000000000e-02,  _prvt_max_sr = 2.5000000000e+00,  _prvt_min_sr = 1.0000000000e+00,  _prvt_EC = 1.5000000000e+00,  _prvt_Buf_c = 2.0000000000e-01,  _prvt_K_buf_c = 1.0000000000e-03,  _prvt_Buf_sr = 1.0000000000e+01,  _prvt_K_buf_sr = 3.0000000000e-01,  _prvt_Buf_ss = 4.0000000000e-01,  _prvt_K_buf_ss = 2.5000000000e-04,  _prvt_V_sr = 1.0940000000e-03,  _prvt_V_c = 1.6404000000e-02,  _prvt_Cm = 1.8500000000e-01,  _prvt_V_ss = 5.4680000000e-05, 
			//private aux variables
			 _prvt_calc_i_Stim=0.0,  _prvt_calc_E_Na=0.0,  _prvt_calc_E_K=0.0,  _prvt_calc_E_Ks=0.0,  _prvt_calc_E_Ca=0.0,  _prvt_calc_alpha_K1=0.0,  _prvt_calc_beta_K1=0.0,  _prvt_calc_xK1_inf=0.0,  _prvt_calc_i_K1=0.0,  _prvt_calc_i_Kr=0.0,  _prvt_calc_xr1_inf=0.0,  _prvt_calc_alpha_xr1=0.0,  _prvt_calc_beta_xr1=0.0,  _prvt_calc_tau_xr1=0.0,  _prvt_calc_xr2_inf=0.0,  _prvt_calc_alpha_xr2=0.0,  _prvt_calc_beta_xr2=0.0,  _prvt_calc_tau_xr2=0.0,  _prvt_calc_i_Ks=0.0,  _prvt_calc_xs_inf=0.0,  _prvt_calc_alpha_xs=0.0,  _prvt_calc_beta_xs=0.0,  _prvt_calc_tau_xs=0.0,  _prvt_calc_i_Na=0.0,  _prvt_calc_m_inf=0.0,  _prvt_calc_alpha_m=0.0,  _prvt_calc_beta_m=0.0,  _prvt_calc_tau_m=0.0,  _prvt_calc_h_inf=0.0,  _prvt_calc_alpha_h=0.0,  _prvt_calc_beta_h=0.0,  _prvt_calc_tau_h=0.0,  _prvt_calc_j_inf=0.0,  _prvt_calc_alpha_j=0.0,  _prvt_calc_beta_j=0.0,  _prvt_calc_tau_j=0.0,  _prvt_calc_i_b_Na=0.0,  _prvt_calc_i_CaL=0.0,  _prvt_calc_d_inf=0.0,  _prvt_calc_alpha_d=0.0,  _prvt_calc_beta_d=0.0,  _prvt_calc_gamma_d=0.0,  _prvt_calc_tau_d=0.0,  _prvt_calc_f_inf=0.0,  _prvt_calc_tau_f=0.0,  _prvt_calc_f2_inf=0.0,  _prvt_calc_tau_f2=0.0,  _prvt_calc_fCass_inf=0.0,  _prvt_calc_tau_fCass=0.0,  _prvt_calc_i_b_Ca=0.0,  _prvt_calc_i_to=0.0,  _prvt_calc_s_inf=0.0,  _prvt_calc_tau_s=0.0,  _prvt_calc_r_inf=0.0,  _prvt_calc_tau_r=0.0,  _prvt_calc_i_NaK=0.0,  _prvt_calc_i_NaCa=0.0,  _prvt_calc_i_p_Ca=0.0,  _prvt_calc_i_p_K=0.0,  _prvt_calc_i_rel=0.0,  _prvt_calc_i_up=0.0,  _prvt_calc_i_leak=0.0,  _prvt_calc_i_xfer=0.0,  _prvt_calc_O=0.0,  _prvt_calc_k1=0.0,  _prvt_calc_k2=0.0,  _prvt_calc_kcasr=0.0,  _prvt_calc_Ca_i_bufc=0.0,  _prvt_calc_Ca_sr_bufsr=0.0,  _prvt_calc_Ca_ss_bufss=0.0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_Xr1_lado_direito_,  _prvt_Xr2_lado_direito_,  _prvt_Xs_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h_lado_direito_,  _prvt_j_lado_direito_,  _prvt_d_lado_direito_,  _prvt_f_lado_direito_,  _prvt_f2_lado_direito_,  _prvt_fCass_lado_direito_,  _prvt_s_lado_direito_,  _prvt_r_lado_direito_,  _prvt_R_prime_lado_direito_,  _prvt_Ca_i_lado_direito_,  _prvt_Ca_SR_lado_direito_,  _prvt_Ca_ss_lado_direito_,  _prvt_Na_i_lado_direito_,  _prvt_K_i_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;
			__NEW_[0] = __OLD_[0] = -8.5423000000e+01;
			__NEW_[1] = __OLD_[1] = 1.6500000000e-02;
			__NEW_[2] = __OLD_[2] = 4.7300000000e-01;
			__NEW_[3] = __OLD_[3] = 1.7400000000e-02;
			__NEW_[4] = __OLD_[4] = 1.6500000000e-03;
			__NEW_[5] = __OLD_[5] = 7.4900000000e-01;
			__NEW_[6] = __OLD_[6] = 6.7880000000e-01;
			__NEW_[7] = __OLD_[7] = 3.2880000000e-05;
			__NEW_[8] = __OLD_[8] = 7.0260000000e-01;
			__NEW_[9] = __OLD_[9] = 9.5260000000e-01;
			__NEW_[10] = __OLD_[10] = 9.9420000000e-01;
			__NEW_[11] = __OLD_[11] = 9.9999800000e-01;
			__NEW_[12] = __OLD_[12] = 2.3470000000e-08;
			__NEW_[13] = __OLD_[13] = 8.9780000000e-01;
			__NEW_[14] = __OLD_[14] = 1.5300000000e-04;
			__NEW_[15] = __OLD_[15] = 4.2720000000e+00;
			__NEW_[16] = __OLD_[16] = 4.2000000000e-04;
			__NEW_[17] = __OLD_[17] = 1.0132000000e+01;
			__NEW_[18] = __OLD_[18] = 1.3852000000e+02;
			_prvt_time_new += _prvt_dtime;
			if(omp_get_thread_num()==_prvt_tree_thread[0])
			{
				_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
				_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[17])));
				_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[18])));
				_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[18]+(_prvt_P_kna*__OLD_[17])))));
				_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[14])));
				_prvt_calc_i_CaL = ((((_prvt_g_CaL*__OLD_[7]*__OLD_[8]*__OLD_[9]*__OLD_[10]*4.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*((2.5000000000e-01*__OLD_[16]*exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))))-_prvt_Ca_o))/(exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00));
				_prvt_calc_i_NaK = (((((_prvt_P_NaK*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk))*__OLD_[17])/(__OLD_[17]+_prvt_K_mNa))/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.5300000000e-02*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_NaCa = ((_prvt_K_NaCa*((exp(((_prvt_gamma*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[17],3.0000000000e+00)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,3.0000000000e+00)*__OLD_[14]*_prvt_alpha)))/((pow(_prvt_Km_Nai,3.0000000000e+00)+pow(_prvt_Na_o,3.0000000000e+00))*(_prvt_Km_Ca+_prvt_Ca_o)*(1.0000000000e+00+(_prvt_K_sat*exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))))));
				_prvt_calc_i_p_Ca = ((_prvt_g_pCa*__OLD_[14])/(__OLD_[14]+_prvt_K_pCa));
				_prvt_calc_i_up = (_prvt_Vmax_up/(1.0000000000e+00+(pow(_prvt_K_up,2.0000000000e+00)/pow(__OLD_[14],2.0000000000e+00))));
				_prvt_calc_i_leak = (_prvt_V_leak*(__OLD_[15]-__OLD_[14]));
				_prvt_calc_i_xfer = (_prvt_V_xfer*(__OLD_[16]-__OLD_[14]));
				_prvt_calc_kcasr = (_prvt_max_sr-((_prvt_max_sr-_prvt_min_sr)/(1.0000000000e+00+pow((_prvt_EC/__OLD_[15]),2.0000000000e+00))));
				_prvt_calc_Ca_i_bufc = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_c*_prvt_K_buf_c)/pow((__OLD_[14]+_prvt_K_buf_c),2.0000000000e+00))));
				_prvt_calc_Ca_sr_bufsr = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_sr*_prvt_K_buf_sr)/pow((__OLD_[15]+_prvt_K_buf_sr),2.0000000000e+00))));
				_prvt_calc_Ca_ss_bufss = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_ss*_prvt_K_buf_ss)/pow((__OLD_[16]+_prvt_K_buf_ss),2.0000000000e+00))));
				_prvt_calc_alpha_K1 = (1.0000000000e-01/(1.0000000000e+00+exp((6.0000000000e-02*((__OLD_[0]-_prvt_calc_E_K)-2.0000000000e+02)))));
				_prvt_calc_beta_K1 = (((3.0000000000e+00*exp((2.0000000000e-04*((__OLD_[0]-_prvt_calc_E_K)+1.0000000000e+02))))+exp((1.0000000000e-01*((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01))))/(1.0000000000e+00+exp(((-5.0000000000e-01)*(__OLD_[0]-_prvt_calc_E_K)))));
				_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[1]*__OLD_[2]*(__OLD_[0]-_prvt_calc_E_K)*pow((_prvt_K_o/5.4000000000e+00),1.0/2.0));
				_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
				_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*__OLD_[6]*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
				_prvt_calc_i_to = (_prvt_g_to*__OLD_[12]*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_p_K = ((_prvt_g_pK*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/5.9800000000e+00))));
				_prvt_calc_xK1_inf = (_prvt_calc_alpha_K1/(_prvt_calc_alpha_K1+_prvt_calc_beta_K1));
				_prvt_calc_k1 = (_prvt_k1_prime/_prvt_calc_kcasr);
				_prvt_calc_k2 = (_prvt_k2_prime*_prvt_calc_kcasr);
				_prvt_calc_i_K1 = (_prvt_g_K1*_prvt_calc_xK1_inf*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_O = ((_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00)*__OLD_[13])/(_prvt_k3+(_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00))));
				_prvt_calc_i_rel = (_prvt_V_rel*_prvt_calc_O*(__OLD_[15]-__OLD_[16]));
				__K1_[0]= (-(_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_CaL+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_NaCa+_prvt_calc_i_b_Ca+_prvt_calc_i_p_K+_prvt_calc_i_p_Ca+_prvt_calc_i_Stim));
				__NEW_[0]= __K1_[0] * _prvt_dtime + __OLD_[0];
				__K1_[13]= (((-_prvt_calc_k2)*__OLD_[16]*__OLD_[13])+(_prvt_k4*(1.0000000000e+00-__OLD_[13])));
				__NEW_[13]= __K1_[13] * _prvt_dtime + __OLD_[13];
				__K1_[14]= (_prvt_calc_Ca_i_bufc*(((((_prvt_calc_i_leak-_prvt_calc_i_up)*_prvt_V_sr)/_prvt_V_c)+_prvt_calc_i_xfer)-((((_prvt_calc_i_b_Ca+_prvt_calc_i_p_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Cm)/(2.0000000000e+00*_prvt_V_c*_prvt_F))));
				__NEW_[14]= __K1_[14] * _prvt_dtime + __OLD_[14];
				__K1_[15]= (_prvt_calc_Ca_sr_bufsr*(_prvt_calc_i_up-(_prvt_calc_i_rel+_prvt_calc_i_leak)));
				__NEW_[15]= __K1_[15] * _prvt_dtime + __OLD_[15];
				__K1_[16]= (_prvt_calc_Ca_ss_bufss*(((((-_prvt_calc_i_CaL)*_prvt_Cm)/(2.0000000000e+00*_prvt_V_ss*_prvt_F))+((_prvt_calc_i_rel*_prvt_V_sr)/_prvt_V_ss))-((_prvt_calc_i_xfer*_prvt_V_c)/_prvt_V_ss)));
				__NEW_[16]= __K1_[16] * _prvt_dtime + __OLD_[16];
				__K1_[17]= (((-(_prvt_calc_i_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
				__NEW_[17]= __K1_[17] * _prvt_dtime + __OLD_[17];
				__K1_[18]= (((-((_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_p_K+_prvt_calc_i_Stim)-(2.0000000000e+00*_prvt_calc_i_NaK)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
				__NEW_[18]= __K1_[18] * _prvt_dtime + __OLD_[18];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[1])
			{
				_prvt_calc_xr1_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-2.6000000000e+01)-__OLD_[0])/7.0000000000e+00))));
				_prvt_calc_alpha_xr1 = (4.5000000000e+02/(1.0000000000e+00+exp((((-4.5000000000e+01)-__OLD_[0])/1.0000000000e+01))));
				_prvt_calc_beta_xr1 = (6.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.1500000000e+01))));
				_prvt_calc_tau_xr1 = (1.0000000000e+00*_prvt_calc_alpha_xr1*_prvt_calc_beta_xr1);
				__K1_[1]= ((_prvt_calc_xr1_inf-__OLD_[1])/_prvt_calc_tau_xr1);
				__NEW_[1]= __K1_[1] * _prvt_dtime + __OLD_[1];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[2])
			{
				_prvt_calc_xr2_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+8.8000000000e+01)/2.4000000000e+01))));
				_prvt_calc_alpha_xr2 = (3.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/2.0000000000e+01))));
				_prvt_calc_beta_xr2 = (1.1200000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-6.0000000000e+01)/2.0000000000e+01))));
				_prvt_calc_tau_xr2 = (1.0000000000e+00*_prvt_calc_alpha_xr2*_prvt_calc_beta_xr2);
				__K1_[2]= ((_prvt_calc_xr2_inf-__OLD_[2])/_prvt_calc_tau_xr2);
				__NEW_[2]= __K1_[2] * _prvt_dtime + __OLD_[2];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[3])
			{
				_prvt_calc_xs_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-5.0000000000e+00)-__OLD_[0])/1.4000000000e+01))));
				_prvt_calc_alpha_xs = (1.4000000000e+03/pow((1.0000000000e+00+exp(((5.0000000000e+00-__OLD_[0])/6.0000000000e+00))),1.0/2.0));
				_prvt_calc_beta_xs = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-3.5000000000e+01)/1.5000000000e+01))));
				_prvt_calc_tau_xs = ((1.0000000000e+00*_prvt_calc_alpha_xs*_prvt_calc_beta_xs)+8.0000000000e+01);
				__K1_[3]= ((_prvt_calc_xs_inf-__OLD_[3])/_prvt_calc_tau_xs);
				__NEW_[3]= __K1_[3] * _prvt_dtime + __OLD_[3];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[4])
			{
				_prvt_calc_m_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp((((-5.6860000000e+01)-__OLD_[0])/9.0300000000e+00))),2.0000000000e+00));
				_prvt_calc_alpha_m = (1.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/5.0000000000e+00))));
				_prvt_calc_beta_m = ((1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/5.0000000000e+00))))+(1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]-5.0000000000e+01)/2.0000000000e+02)))));
				_prvt_calc_tau_m = (1.0000000000e+00*_prvt_calc_alpha_m*_prvt_calc_beta_m);
				__K1_[4]= ((_prvt_calc_m_inf-__OLD_[4])/_prvt_calc_tau_m);
				__NEW_[4]= __K1_[4] * _prvt_dtime + __OLD_[4];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[5])
			{
				_prvt_calc_h_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
				_prvt_calc_alpha_h = ((__OLD_[0]<(-4.0000000000e+01)))
?((5.7000000000e-02*exp(((-(__OLD_[0]+8.0000000000e+01))/6.8000000000e+00))))
:(0.0000000000e+00);
				_prvt_calc_beta_h = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.7000000000e+00*exp((7.9000000000e-02*__OLD_[0])))+(3.1000000000e+05*exp((3.4850000000e-01*__OLD_[0])))))
:((7.7000000000e-01/(1.3000000000e-01*(1.0000000000e+00+exp(((__OLD_[0]+1.0660000000e+01)/(-1.1100000000e+01)))))));
				_prvt_calc_tau_h = (1.0000000000e+00/(_prvt_calc_alpha_h+_prvt_calc_beta_h));
				__K1_[5]= ((_prvt_calc_h_inf-__OLD_[5])/_prvt_calc_tau_h);
				__NEW_[5]= __K1_[5] * _prvt_dtime + __OLD_[5];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[6])
			{
				_prvt_calc_j_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
				_prvt_calc_alpha_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((((((-2.5428000000e+04)*exp((2.4440000000e-01*__OLD_[0])))-(6.9480000000e-06*exp(((-4.3910000000e-02)*__OLD_[0]))))*(__OLD_[0]+3.7780000000e+01))/1.0000000000e+00)/(1.0000000000e+00+exp((3.1100000000e-01*(__OLD_[0]+7.9230000000e+01))))))
:(0.0000000000e+00);
				_prvt_calc_beta_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.4240000000e-02*exp(((-1.0520000000e-02)*__OLD_[0])))/(1.0000000000e+00+exp(((-1.3780000000e-01)*(__OLD_[0]+4.0140000000e+01))))))
:(((6.0000000000e-01*exp((5.7000000000e-02*__OLD_[0])))/(1.0000000000e+00+exp(((-1.0000000000e-01)*(__OLD_[0]+3.2000000000e+01))))));
				_prvt_calc_tau_j = (1.0000000000e+00/(_prvt_calc_alpha_j+_prvt_calc_beta_j));
				__K1_[6]= ((_prvt_calc_j_inf-__OLD_[6])/_prvt_calc_tau_j);
				__NEW_[6]= __K1_[6] * _prvt_dtime + __OLD_[6];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[7])
			{
				_prvt_calc_d_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-8.0000000000e+00)-__OLD_[0])/7.5000000000e+00))));
				_prvt_calc_alpha_d = ((1.4000000000e+00/(1.0000000000e+00+exp((((-3.5000000000e+01)-__OLD_[0])/1.3000000000e+01))))+2.5000000000e-01);
				_prvt_calc_beta_d = (1.4000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.0000000000e+00)/5.0000000000e+00))));
				_prvt_calc_gamma_d = (1.0000000000e+00/(1.0000000000e+00+exp(((5.0000000000e+01-__OLD_[0])/2.0000000000e+01))));
				_prvt_calc_tau_d = ((1.0000000000e+00*_prvt_calc_alpha_d*_prvt_calc_beta_d)+_prvt_calc_gamma_d);
				__K1_[7]= ((_prvt_calc_d_inf-__OLD_[7])/_prvt_calc_tau_d);
				__NEW_[7]= __K1_[7] * _prvt_dtime + __OLD_[7];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[8])
			{
				_prvt_calc_f_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/7.0000000000e+00))));
				_prvt_calc_tau_f = ((1.1025000000e+03*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.2500000000e+02)))+(2.0000000000e+02/(1.0000000000e+00+exp(((1.3000000000e+01-__OLD_[0])/1.0000000000e+01))))+(1.8000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01))))+2.0000000000e+01);
				__K1_[8]= ((_prvt_calc_f_inf-__OLD_[8])/_prvt_calc_tau_f);
				__NEW_[8]= __K1_[8] * _prvt_dtime + __OLD_[8];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[9])
			{
				_prvt_calc_f2_inf = ((6.7000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/7.0000000000e+00))))+3.3000000000e-01);
				_prvt_calc_tau_f2 = ((5.6200000000e+02*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.4000000000e+02)))+(3.1000000000e+01/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/1.0000000000e+01))))+(8.0000000000e+01/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01)))));
				__K1_[9]= ((_prvt_calc_f2_inf-__OLD_[9])/_prvt_calc_tau_f2);
				__NEW_[9]= __K1_[9] * _prvt_dtime + __OLD_[9];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[10])
			{
				_prvt_calc_fCass_inf = ((6.0000000000e-01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+4.0000000000e-01);
				_prvt_calc_tau_fCass = ((8.0000000000e+01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+2.0000000000e+00);
				__K1_[10]= ((_prvt_calc_fCass_inf-__OLD_[10])/_prvt_calc_tau_fCass);
				__NEW_[10]= __K1_[10] * _prvt_dtime + __OLD_[10];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[11])
			{
				_prvt_calc_s_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/5.0000000000e+00))));
				_prvt_calc_tau_s = ((8.5000000000e+01*exp(((-pow((__OLD_[0]+4.5000000000e+01),2.0000000000e+00))/3.2000000000e+02)))+(5.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-2.0000000000e+01)/5.0000000000e+00))))+3.0000000000e+00);
				__K1_[11]= ((_prvt_calc_s_inf-__OLD_[11])/_prvt_calc_tau_s);
				__NEW_[11]= __K1_[11] * _prvt_dtime + __OLD_[11];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[12])
			{
				_prvt_calc_r_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((2.0000000000e+01-__OLD_[0])/6.0000000000e+00))));
				_prvt_calc_tau_r = ((9.5000000000e+00*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.8000000000e+03)))+8.0000000000e-01);
				__K1_[12]= ((_prvt_calc_r_inf-__OLD_[12])/_prvt_calc_tau_r);
				__NEW_[12]= __K1_[12] * _prvt_dtime + __OLD_[12];
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
					_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[17])));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[18])));
					_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[18]+(_prvt_P_kna*__OLD_[17])))));
					_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[14])));
					_prvt_calc_i_CaL = ((((_prvt_g_CaL*__OLD_[7]*__OLD_[8]*__OLD_[9]*__OLD_[10]*4.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*((2.5000000000e-01*__OLD_[16]*exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))))-_prvt_Ca_o))/(exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00));
					_prvt_calc_i_NaK = (((((_prvt_P_NaK*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk))*__OLD_[17])/(__OLD_[17]+_prvt_K_mNa))/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.5300000000e-02*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_NaCa = ((_prvt_K_NaCa*((exp(((_prvt_gamma*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[17],3.0000000000e+00)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,3.0000000000e+00)*__OLD_[14]*_prvt_alpha)))/((pow(_prvt_Km_Nai,3.0000000000e+00)+pow(_prvt_Na_o,3.0000000000e+00))*(_prvt_Km_Ca+_prvt_Ca_o)*(1.0000000000e+00+(_prvt_K_sat*exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))))));
					_prvt_calc_i_p_Ca = ((_prvt_g_pCa*__OLD_[14])/(__OLD_[14]+_prvt_K_pCa));
					_prvt_calc_i_up = (_prvt_Vmax_up/(1.0000000000e+00+(pow(_prvt_K_up,2.0000000000e+00)/pow(__OLD_[14],2.0000000000e+00))));
					_prvt_calc_i_leak = (_prvt_V_leak*(__OLD_[15]-__OLD_[14]));
					_prvt_calc_i_xfer = (_prvt_V_xfer*(__OLD_[16]-__OLD_[14]));
					_prvt_calc_kcasr = (_prvt_max_sr-((_prvt_max_sr-_prvt_min_sr)/(1.0000000000e+00+pow((_prvt_EC/__OLD_[15]),2.0000000000e+00))));
					_prvt_calc_Ca_i_bufc = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_c*_prvt_K_buf_c)/pow((__OLD_[14]+_prvt_K_buf_c),2.0000000000e+00))));
					_prvt_calc_Ca_sr_bufsr = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_sr*_prvt_K_buf_sr)/pow((__OLD_[15]+_prvt_K_buf_sr),2.0000000000e+00))));
					_prvt_calc_Ca_ss_bufss = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_ss*_prvt_K_buf_ss)/pow((__OLD_[16]+_prvt_K_buf_ss),2.0000000000e+00))));
					_prvt_calc_alpha_K1 = (1.0000000000e-01/(1.0000000000e+00+exp((6.0000000000e-02*((__OLD_[0]-_prvt_calc_E_K)-2.0000000000e+02)))));
					_prvt_calc_beta_K1 = (((3.0000000000e+00*exp((2.0000000000e-04*((__OLD_[0]-_prvt_calc_E_K)+1.0000000000e+02))))+exp((1.0000000000e-01*((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01))))/(1.0000000000e+00+exp(((-5.0000000000e-01)*(__OLD_[0]-_prvt_calc_E_K)))));
					_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[1]*__OLD_[2]*(__OLD_[0]-_prvt_calc_E_K)*pow((_prvt_K_o/5.4000000000e+00),1.0/2.0));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
					_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*__OLD_[6]*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_to = (_prvt_g_to*__OLD_[12]*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_p_K = ((_prvt_g_pK*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/5.9800000000e+00))));
					_prvt_calc_xK1_inf = (_prvt_calc_alpha_K1/(_prvt_calc_alpha_K1+_prvt_calc_beta_K1));
					_prvt_calc_k1 = (_prvt_k1_prime/_prvt_calc_kcasr);
					_prvt_calc_k2 = (_prvt_k2_prime*_prvt_calc_kcasr);
					_prvt_calc_i_K1 = (_prvt_g_K1*_prvt_calc_xK1_inf*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_O = ((_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00)*__OLD_[13])/(_prvt_k3+(_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00))));
					_prvt_calc_i_rel = (_prvt_V_rel*_prvt_calc_O*(__OLD_[15]-__OLD_[16]));
					__K2_[0]= (-(_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_CaL+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_NaCa+_prvt_calc_i_b_Ca+_prvt_calc_i_p_K+_prvt_calc_i_p_Ca+_prvt_calc_i_Stim));
					_prvt_aux_tol = fabs(__OLD_[0])*_prvt_rel_tol_;
					__TOL_[0] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[0] = fabs((_prvt_dtime/2) * (__K1_[0] - __K2_[0])/__TOL_[0]);
					__K2_[13]= (((-_prvt_calc_k2)*__OLD_[16]*__OLD_[13])+(_prvt_k4*(1.0000000000e+00-__OLD_[13])));
					_prvt_aux_tol = fabs(__OLD_[13])*_prvt_rel_tol_;
					__TOL_[13] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[13] = fabs((_prvt_dtime/2) * (__K1_[13] - __K2_[13])/__TOL_[13]);
					__K2_[14]= (_prvt_calc_Ca_i_bufc*(((((_prvt_calc_i_leak-_prvt_calc_i_up)*_prvt_V_sr)/_prvt_V_c)+_prvt_calc_i_xfer)-((((_prvt_calc_i_b_Ca+_prvt_calc_i_p_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Cm)/(2.0000000000e+00*_prvt_V_c*_prvt_F))));
					_prvt_aux_tol = fabs(__OLD_[14])*_prvt_rel_tol_;
					__TOL_[14] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[14] = fabs((_prvt_dtime/2) * (__K1_[14] - __K2_[14])/__TOL_[14]);
					__K2_[15]= (_prvt_calc_Ca_sr_bufsr*(_prvt_calc_i_up-(_prvt_calc_i_rel+_prvt_calc_i_leak)));
					_prvt_aux_tol = fabs(__OLD_[15])*_prvt_rel_tol_;
					__TOL_[15] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[15] = fabs((_prvt_dtime/2) * (__K1_[15] - __K2_[15])/__TOL_[15]);
					__K2_[16]= (_prvt_calc_Ca_ss_bufss*(((((-_prvt_calc_i_CaL)*_prvt_Cm)/(2.0000000000e+00*_prvt_V_ss*_prvt_F))+((_prvt_calc_i_rel*_prvt_V_sr)/_prvt_V_ss))-((_prvt_calc_i_xfer*_prvt_V_c)/_prvt_V_ss)));
					_prvt_aux_tol = fabs(__OLD_[16])*_prvt_rel_tol_;
					__TOL_[16] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[16] = fabs((_prvt_dtime/2) * (__K1_[16] - __K2_[16])/__TOL_[16]);
					__K2_[17]= (((-(_prvt_calc_i_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
					_prvt_aux_tol = fabs(__OLD_[17])*_prvt_rel_tol_;
					__TOL_[17] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[17] = fabs((_prvt_dtime/2) * (__K1_[17] - __K2_[17])/__TOL_[17]);
					__K2_[18]= (((-((_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_p_K+_prvt_calc_i_Stim)-(2.0000000000e+00*_prvt_calc_i_NaK)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
					_prvt_aux_tol = fabs(__OLD_[18])*_prvt_rel_tol_;
					__TOL_[18] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[18] = fabs((_prvt_dtime/2) * (__K1_[18] - __K2_[18])/__TOL_[18]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[1])
				{
					_prvt_calc_xr1_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-2.6000000000e+01)-__OLD_[0])/7.0000000000e+00))));
					_prvt_calc_alpha_xr1 = (4.5000000000e+02/(1.0000000000e+00+exp((((-4.5000000000e+01)-__OLD_[0])/1.0000000000e+01))));
					_prvt_calc_beta_xr1 = (6.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.1500000000e+01))));
					_prvt_calc_tau_xr1 = (1.0000000000e+00*_prvt_calc_alpha_xr1*_prvt_calc_beta_xr1);
					__K2_[1]= ((_prvt_calc_xr1_inf-__OLD_[1])/_prvt_calc_tau_xr1);
					_prvt_aux_tol = fabs(__OLD_[1])*_prvt_rel_tol_;
					__TOL_[1] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[1] = fabs((_prvt_dtime/2) * (__K1_[1] - __K2_[1])/__TOL_[1]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[2])
				{
					_prvt_calc_xr2_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+8.8000000000e+01)/2.4000000000e+01))));
					_prvt_calc_alpha_xr2 = (3.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/2.0000000000e+01))));
					_prvt_calc_beta_xr2 = (1.1200000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-6.0000000000e+01)/2.0000000000e+01))));
					_prvt_calc_tau_xr2 = (1.0000000000e+00*_prvt_calc_alpha_xr2*_prvt_calc_beta_xr2);
					__K2_[2]= ((_prvt_calc_xr2_inf-__OLD_[2])/_prvt_calc_tau_xr2);
					_prvt_aux_tol = fabs(__OLD_[2])*_prvt_rel_tol_;
					__TOL_[2] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[2] = fabs((_prvt_dtime/2) * (__K1_[2] - __K2_[2])/__TOL_[2]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[3])
				{
					_prvt_calc_xs_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-5.0000000000e+00)-__OLD_[0])/1.4000000000e+01))));
					_prvt_calc_alpha_xs = (1.4000000000e+03/pow((1.0000000000e+00+exp(((5.0000000000e+00-__OLD_[0])/6.0000000000e+00))),1.0/2.0));
					_prvt_calc_beta_xs = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-3.5000000000e+01)/1.5000000000e+01))));
					_prvt_calc_tau_xs = ((1.0000000000e+00*_prvt_calc_alpha_xs*_prvt_calc_beta_xs)+8.0000000000e+01);
					__K2_[3]= ((_prvt_calc_xs_inf-__OLD_[3])/_prvt_calc_tau_xs);
					_prvt_aux_tol = fabs(__OLD_[3])*_prvt_rel_tol_;
					__TOL_[3] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[3] = fabs((_prvt_dtime/2) * (__K1_[3] - __K2_[3])/__TOL_[3]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[4])
				{
					_prvt_calc_m_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp((((-5.6860000000e+01)-__OLD_[0])/9.0300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_m = (1.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/5.0000000000e+00))));
					_prvt_calc_beta_m = ((1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/5.0000000000e+00))))+(1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]-5.0000000000e+01)/2.0000000000e+02)))));
					_prvt_calc_tau_m = (1.0000000000e+00*_prvt_calc_alpha_m*_prvt_calc_beta_m);
					__K2_[4]= ((_prvt_calc_m_inf-__OLD_[4])/_prvt_calc_tau_m);
					_prvt_aux_tol = fabs(__OLD_[4])*_prvt_rel_tol_;
					__TOL_[4] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[4] = fabs((_prvt_dtime/2) * (__K1_[4] - __K2_[4])/__TOL_[4]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[5])
				{
					_prvt_calc_h_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_h = ((__OLD_[0]<(-4.0000000000e+01)))
?((5.7000000000e-02*exp(((-(__OLD_[0]+8.0000000000e+01))/6.8000000000e+00))))
:(0.0000000000e+00);
					_prvt_calc_beta_h = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.7000000000e+00*exp((7.9000000000e-02*__OLD_[0])))+(3.1000000000e+05*exp((3.4850000000e-01*__OLD_[0])))))
:((7.7000000000e-01/(1.3000000000e-01*(1.0000000000e+00+exp(((__OLD_[0]+1.0660000000e+01)/(-1.1100000000e+01)))))));
					_prvt_calc_tau_h = (1.0000000000e+00/(_prvt_calc_alpha_h+_prvt_calc_beta_h));
					__K2_[5]= ((_prvt_calc_h_inf-__OLD_[5])/_prvt_calc_tau_h);
					_prvt_aux_tol = fabs(__OLD_[5])*_prvt_rel_tol_;
					__TOL_[5] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[5] = fabs((_prvt_dtime/2) * (__K1_[5] - __K2_[5])/__TOL_[5]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[6])
				{
					_prvt_calc_j_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((((((-2.5428000000e+04)*exp((2.4440000000e-01*__OLD_[0])))-(6.9480000000e-06*exp(((-4.3910000000e-02)*__OLD_[0]))))*(__OLD_[0]+3.7780000000e+01))/1.0000000000e+00)/(1.0000000000e+00+exp((3.1100000000e-01*(__OLD_[0]+7.9230000000e+01))))))
:(0.0000000000e+00);
					_prvt_calc_beta_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.4240000000e-02*exp(((-1.0520000000e-02)*__OLD_[0])))/(1.0000000000e+00+exp(((-1.3780000000e-01)*(__OLD_[0]+4.0140000000e+01))))))
:(((6.0000000000e-01*exp((5.7000000000e-02*__OLD_[0])))/(1.0000000000e+00+exp(((-1.0000000000e-01)*(__OLD_[0]+3.2000000000e+01))))));
					_prvt_calc_tau_j = (1.0000000000e+00/(_prvt_calc_alpha_j+_prvt_calc_beta_j));
					__K2_[6]= ((_prvt_calc_j_inf-__OLD_[6])/_prvt_calc_tau_j);
					_prvt_aux_tol = fabs(__OLD_[6])*_prvt_rel_tol_;
					__TOL_[6] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[6] = fabs((_prvt_dtime/2) * (__K1_[6] - __K2_[6])/__TOL_[6]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[7])
				{
					_prvt_calc_d_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-8.0000000000e+00)-__OLD_[0])/7.5000000000e+00))));
					_prvt_calc_alpha_d = ((1.4000000000e+00/(1.0000000000e+00+exp((((-3.5000000000e+01)-__OLD_[0])/1.3000000000e+01))))+2.5000000000e-01);
					_prvt_calc_beta_d = (1.4000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.0000000000e+00)/5.0000000000e+00))));
					_prvt_calc_gamma_d = (1.0000000000e+00/(1.0000000000e+00+exp(((5.0000000000e+01-__OLD_[0])/2.0000000000e+01))));
					_prvt_calc_tau_d = ((1.0000000000e+00*_prvt_calc_alpha_d*_prvt_calc_beta_d)+_prvt_calc_gamma_d);
					__K2_[7]= ((_prvt_calc_d_inf-__OLD_[7])/_prvt_calc_tau_d);
					_prvt_aux_tol = fabs(__OLD_[7])*_prvt_rel_tol_;
					__TOL_[7] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[7] = fabs((_prvt_dtime/2) * (__K1_[7] - __K2_[7])/__TOL_[7]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[8])
				{
					_prvt_calc_f_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/7.0000000000e+00))));
					_prvt_calc_tau_f = ((1.1025000000e+03*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.2500000000e+02)))+(2.0000000000e+02/(1.0000000000e+00+exp(((1.3000000000e+01-__OLD_[0])/1.0000000000e+01))))+(1.8000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01))))+2.0000000000e+01);
					__K2_[8]= ((_prvt_calc_f_inf-__OLD_[8])/_prvt_calc_tau_f);
					_prvt_aux_tol = fabs(__OLD_[8])*_prvt_rel_tol_;
					__TOL_[8] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[8] = fabs((_prvt_dtime/2) * (__K1_[8] - __K2_[8])/__TOL_[8]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[9])
				{
					_prvt_calc_f2_inf = ((6.7000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/7.0000000000e+00))))+3.3000000000e-01);
					_prvt_calc_tau_f2 = ((5.6200000000e+02*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.4000000000e+02)))+(3.1000000000e+01/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/1.0000000000e+01))))+(8.0000000000e+01/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01)))));
					__K2_[9]= ((_prvt_calc_f2_inf-__OLD_[9])/_prvt_calc_tau_f2);
					_prvt_aux_tol = fabs(__OLD_[9])*_prvt_rel_tol_;
					__TOL_[9] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[9] = fabs((_prvt_dtime/2) * (__K1_[9] - __K2_[9])/__TOL_[9]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[10])
				{
					_prvt_calc_fCass_inf = ((6.0000000000e-01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+4.0000000000e-01);
					_prvt_calc_tau_fCass = ((8.0000000000e+01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+2.0000000000e+00);
					__K2_[10]= ((_prvt_calc_fCass_inf-__OLD_[10])/_prvt_calc_tau_fCass);
					_prvt_aux_tol = fabs(__OLD_[10])*_prvt_rel_tol_;
					__TOL_[10] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[10] = fabs((_prvt_dtime/2) * (__K1_[10] - __K2_[10])/__TOL_[10]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[11])
				{
					_prvt_calc_s_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/5.0000000000e+00))));
					_prvt_calc_tau_s = ((8.5000000000e+01*exp(((-pow((__OLD_[0]+4.5000000000e+01),2.0000000000e+00))/3.2000000000e+02)))+(5.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-2.0000000000e+01)/5.0000000000e+00))))+3.0000000000e+00);
					__K2_[11]= ((_prvt_calc_s_inf-__OLD_[11])/_prvt_calc_tau_s);
					_prvt_aux_tol = fabs(__OLD_[11])*_prvt_rel_tol_;
					__TOL_[11] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[11] = fabs((_prvt_dtime/2) * (__K1_[11] - __K2_[11])/__TOL_[11]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[12])
				{
					_prvt_calc_r_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((2.0000000000e+01-__OLD_[0])/6.0000000000e+00))));
					_prvt_calc_tau_r = ((9.5000000000e+00*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.8000000000e+03)))+8.0000000000e-01);
					__K2_[12]= ((_prvt_calc_r_inf-__OLD_[12])/_prvt_calc_tau_r);
					_prvt_aux_tol = fabs(__OLD_[12])*_prvt_rel_tol_;
					__TOL_[12] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[12] = fabs((_prvt_dtime/2) * (__K1_[12] - __K2_[12])/__TOL_[12]);
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
				int flag = this->getErrorCode(__greatestError_, 1.0);
				_prvt_previous_dt = _prvt_dtime;
				//it doesn't accept the solution and cut h in a half
				if(flag==-1){
					//throw the results away and compute again
					_prvt_dtime = _prvt_dtime/_DECREASE_DT_2_; //cut time step in a half
					if(omp_get_thread_num()==_prvt_tree_thread[0])
					{
						__NEW_[0] = __K1_[0] * _prvt_dtime + __OLD_AUX_[0];
						__NEW_[13] = __K1_[13] * _prvt_dtime + __OLD_AUX_[13];
						__NEW_[14] = __K1_[14] * _prvt_dtime + __OLD_AUX_[14];
						__NEW_[15] = __K1_[15] * _prvt_dtime + __OLD_AUX_[15];
						__NEW_[16] = __K1_[16] * _prvt_dtime + __OLD_AUX_[16];
						__NEW_[17] = __K1_[17] * _prvt_dtime + __OLD_AUX_[17];
						__NEW_[18] = __K1_[18] * _prvt_dtime + __OLD_AUX_[18];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[1] = __K1_[1] * _prvt_dtime + __OLD_AUX_[1];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[2] = __K1_[2] * _prvt_dtime + __OLD_AUX_[2];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[3] = __K1_[3] * _prvt_dtime + __OLD_AUX_[3];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[4] = __K1_[4] * _prvt_dtime + __OLD_AUX_[4];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[5] = __K1_[5] * _prvt_dtime + __OLD_AUX_[5];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[6] = __K1_[6] * _prvt_dtime + __OLD_AUX_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[7] = __K1_[7] * _prvt_dtime + __OLD_AUX_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[8] = __K1_[8] * _prvt_dtime + __OLD_AUX_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[9] = __K1_[9] * _prvt_dtime + __OLD_AUX_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[10] = __K1_[10] * _prvt_dtime + __OLD_AUX_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[11] = __K1_[11] * _prvt_dtime + __OLD_AUX_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[12] = __K1_[12] * _prvt_dtime + __OLD_AUX_[12];
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
							this->Xr1_old_ = __OLD_AUX_[1];
							this->Xr1_new_ = __OLD_[1];
							this->Xr2_old_ = __OLD_AUX_[2];
							this->Xr2_new_ = __OLD_[2];
							this->Xs_old_ = __OLD_AUX_[3];
							this->Xs_new_ = __OLD_[3];
							this->m_old_ = __OLD_AUX_[4];
							this->m_new_ = __OLD_[4];
							this->h_old_ = __OLD_AUX_[5];
							this->h_new_ = __OLD_[5];
							this->j_old_ = __OLD_AUX_[6];
							this->j_new_ = __OLD_[6];
							this->d_old_ = __OLD_AUX_[7];
							this->d_new_ = __OLD_[7];
							this->f_old_ = __OLD_AUX_[8];
							this->f_new_ = __OLD_[8];
							this->f2_old_ = __OLD_AUX_[9];
							this->f2_new_ = __OLD_[9];
							this->fCass_old_ = __OLD_AUX_[10];
							this->fCass_new_ = __OLD_[10];
							this->s_old_ = __OLD_AUX_[11];
							this->s_new_ = __OLD_[11];
							this->r_old_ = __OLD_AUX_[12];
							this->r_new_ = __OLD_[12];
							this->R_prime_old_ = __OLD_AUX_[13];
							this->R_prime_new_ = __OLD_[13];
							this->Ca_i_old_ = __OLD_AUX_[14];
							this->Ca_i_new_ = __OLD_[14];
							this->Ca_SR_old_ = __OLD_AUX_[15];
							this->Ca_SR_new_ = __OLD_[15];
							this->Ca_ss_old_ = __OLD_AUX_[16];
							this->Ca_ss_new_ = __OLD_[16];
							this->Na_i_old_ = __OLD_AUX_[17];
							this->Na_i_new_ = __OLD_[17];
							this->K_i_old_ = __OLD_AUX_[18];
							this->K_i_new_ = __OLD_[18];
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
						__NEW_[13] = __K2_[13] * _prvt_dtime + __OLD_[13];
						__NEW_[14] = __K2_[14] * _prvt_dtime + __OLD_[14];
						__NEW_[15] = __K2_[15] * _prvt_dtime + __OLD_[15];
						__NEW_[16] = __K2_[16] * _prvt_dtime + __OLD_[16];
						__NEW_[17] = __K2_[17] * _prvt_dtime + __OLD_[17];
						__NEW_[18] = __K2_[18] * _prvt_dtime + __OLD_[18];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[1] = __K2_[1] * _prvt_dtime + __OLD_[1];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[2] = __K2_[2] * _prvt_dtime + __OLD_[2];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[3] = __K2_[3] * _prvt_dtime + __OLD_[3];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[4] = __K2_[4] * _prvt_dtime + __OLD_[4];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[5] = __K2_[5] * _prvt_dtime + __OLD_[5];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[6] = __K2_[6] * _prvt_dtime + __OLD_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[7] = __K2_[7] * _prvt_dtime + __OLD_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[8] = __K2_[8] * _prvt_dtime + __OLD_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[9] = __K2_[9] * _prvt_dtime + __OLD_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[10] = __K2_[10] * _prvt_dtime + __OLD_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[11] = __K2_[11] * _prvt_dtime + __OLD_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[12] = __K2_[12] * _prvt_dtime + __OLD_[12];
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
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -5.2000000000e+01,  _prvt_stim_start = 1.0000000000e+02,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 5.0000000000e+02,  _prvt_stim_duration = 1.5000000000e+00,  _prvt_R = 8.3144720000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_F = 9.6485341500e+04,  _prvt_Na_o = 1.4000000000e+02,  _prvt_K_o = 5.4000000000e+00,  _prvt_P_kna = 3.0000000000e-02,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_g_K1 = 5.4050000000e+00,  _prvt_g_Kr = 1.5300000000e-01,  _prvt_g_Ks = 9.8000000000e-02,  _prvt_g_Na = 1.4838000000e+01,  _prvt_g_bna = 2.9000000000e-04,  _prvt_g_CaL = 3.9800000000e-05,  _prvt_g_bca = 5.9200000000e-04,  _prvt_g_to = 2.9400000000e-01,  _prvt_P_NaK = 2.7240000000e+00,  _prvt_K_mk = 1.0000000000e+00,  _prvt_K_mNa = 4.0000000000e+01,  _prvt_K_NaCa = 1.0000000000e+03,  _prvt_gamma = 3.5000000000e-01,  _prvt_alpha = 2.5000000000e+00,  _prvt_Km_Nai = 8.7500000000e+01,  _prvt_Km_Ca = 1.3800000000e+00,  _prvt_K_sat = 1.0000000000e-01,  _prvt_g_pCa = 1.2380000000e-01,  _prvt_K_pCa = 5.0000000000e-04,  _prvt_g_pK = 1.4600000000e-02,  _prvt_V_rel = 1.0200000000e-01,  _prvt_Vmax_up = 6.3750000000e-03,  _prvt_K_up = 2.5000000000e-04,  _prvt_V_leak = 3.6000000000e-04,  _prvt_V_xfer = 3.8000000000e-03,  _prvt_k3 = 6.0000000000e-02,  _prvt_k4 = 5.0000000000e-03,  _prvt_k1_prime = 1.5000000000e-01,  _prvt_k2_prime = 4.5000000000e-02,  _prvt_max_sr = 2.5000000000e+00,  _prvt_min_sr = 1.0000000000e+00,  _prvt_EC = 1.5000000000e+00,  _prvt_Buf_c = 2.0000000000e-01,  _prvt_K_buf_c = 1.0000000000e-03,  _prvt_Buf_sr = 1.0000000000e+01,  _prvt_K_buf_sr = 3.0000000000e-01,  _prvt_Buf_ss = 4.0000000000e-01,  _prvt_K_buf_ss = 2.5000000000e-04,  _prvt_V_sr = 1.0940000000e-03,  _prvt_V_c = 1.6404000000e-02,  _prvt_Cm = 1.8500000000e-01,  _prvt_V_ss = 5.4680000000e-05, 
			//private aux variables
			 _prvt_calc_i_Stim=0.0,  _prvt_calc_E_Na=0.0,  _prvt_calc_E_K=0.0,  _prvt_calc_E_Ks=0.0,  _prvt_calc_E_Ca=0.0,  _prvt_calc_alpha_K1=0.0,  _prvt_calc_beta_K1=0.0,  _prvt_calc_xK1_inf=0.0,  _prvt_calc_i_K1=0.0,  _prvt_calc_i_Kr=0.0,  _prvt_calc_xr1_inf=0.0,  _prvt_calc_alpha_xr1=0.0,  _prvt_calc_beta_xr1=0.0,  _prvt_calc_tau_xr1=0.0,  _prvt_calc_xr2_inf=0.0,  _prvt_calc_alpha_xr2=0.0,  _prvt_calc_beta_xr2=0.0,  _prvt_calc_tau_xr2=0.0,  _prvt_calc_i_Ks=0.0,  _prvt_calc_xs_inf=0.0,  _prvt_calc_alpha_xs=0.0,  _prvt_calc_beta_xs=0.0,  _prvt_calc_tau_xs=0.0,  _prvt_calc_i_Na=0.0,  _prvt_calc_m_inf=0.0,  _prvt_calc_alpha_m=0.0,  _prvt_calc_beta_m=0.0,  _prvt_calc_tau_m=0.0,  _prvt_calc_h_inf=0.0,  _prvt_calc_alpha_h=0.0,  _prvt_calc_beta_h=0.0,  _prvt_calc_tau_h=0.0,  _prvt_calc_j_inf=0.0,  _prvt_calc_alpha_j=0.0,  _prvt_calc_beta_j=0.0,  _prvt_calc_tau_j=0.0,  _prvt_calc_i_b_Na=0.0,  _prvt_calc_i_CaL=0.0,  _prvt_calc_d_inf=0.0,  _prvt_calc_alpha_d=0.0,  _prvt_calc_beta_d=0.0,  _prvt_calc_gamma_d=0.0,  _prvt_calc_tau_d=0.0,  _prvt_calc_f_inf=0.0,  _prvt_calc_tau_f=0.0,  _prvt_calc_f2_inf=0.0,  _prvt_calc_tau_f2=0.0,  _prvt_calc_fCass_inf=0.0,  _prvt_calc_tau_fCass=0.0,  _prvt_calc_i_b_Ca=0.0,  _prvt_calc_i_to=0.0,  _prvt_calc_s_inf=0.0,  _prvt_calc_tau_s=0.0,  _prvt_calc_r_inf=0.0,  _prvt_calc_tau_r=0.0,  _prvt_calc_i_NaK=0.0,  _prvt_calc_i_NaCa=0.0,  _prvt_calc_i_p_Ca=0.0,  _prvt_calc_i_p_K=0.0,  _prvt_calc_i_rel=0.0,  _prvt_calc_i_up=0.0,  _prvt_calc_i_leak=0.0,  _prvt_calc_i_xfer=0.0,  _prvt_calc_O=0.0,  _prvt_calc_k1=0.0,  _prvt_calc_k2=0.0,  _prvt_calc_kcasr=0.0,  _prvt_calc_Ca_i_bufc=0.0,  _prvt_calc_Ca_sr_bufsr=0.0,  _prvt_calc_Ca_ss_bufss=0.0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_Xr1_lado_direito_,  _prvt_Xr2_lado_direito_,  _prvt_Xs_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h_lado_direito_,  _prvt_j_lado_direito_,  _prvt_d_lado_direito_,  _prvt_f_lado_direito_,  _prvt_f2_lado_direito_,  _prvt_fCass_lado_direito_,  _prvt_s_lado_direito_,  _prvt_r_lado_direito_,  _prvt_R_prime_lado_direito_,  _prvt_Ca_i_lado_direito_,  _prvt_Ca_SR_lado_direito_,  _prvt_Ca_ss_lado_direito_,  _prvt_Na_i_lado_direito_,  _prvt_K_i_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;
			__NEW_[0] = __OLD_[0] = -8.5423000000e+01;
			__NEW_[1] = __OLD_[1] = 1.6500000000e-02;
			__NEW_[2] = __OLD_[2] = 4.7300000000e-01;
			__NEW_[3] = __OLD_[3] = 1.7400000000e-02;
			__NEW_[4] = __OLD_[4] = 1.6500000000e-03;
			__NEW_[5] = __OLD_[5] = 7.4900000000e-01;
			__NEW_[6] = __OLD_[6] = 6.7880000000e-01;
			__NEW_[7] = __OLD_[7] = 3.2880000000e-05;
			__NEW_[8] = __OLD_[8] = 7.0260000000e-01;
			__NEW_[9] = __OLD_[9] = 9.5260000000e-01;
			__NEW_[10] = __OLD_[10] = 9.9420000000e-01;
			__NEW_[11] = __OLD_[11] = 9.9999800000e-01;
			__NEW_[12] = __OLD_[12] = 2.3470000000e-08;
			__NEW_[13] = __OLD_[13] = 8.9780000000e-01;
			__NEW_[14] = __OLD_[14] = 1.5300000000e-04;
			__NEW_[15] = __OLD_[15] = 4.2720000000e+00;
			__NEW_[16] = __OLD_[16] = 4.2000000000e-04;
			__NEW_[17] = __OLD_[17] = 1.0132000000e+01;
			__NEW_[18] = __OLD_[18] = 1.3852000000e+02;
			_prvt_time_new += _prvt_dtime;
			if(omp_get_thread_num()==_prvt_tree_thread[0])
			{
				_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
				_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[17])));
				_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[18])));
				_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[18]+(_prvt_P_kna*__OLD_[17])))));
				_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[14])));
				_prvt_calc_i_CaL = ((((_prvt_g_CaL*__OLD_[7]*__OLD_[8]*__OLD_[9]*__OLD_[10]*4.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*((2.5000000000e-01*__OLD_[16]*exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))))-_prvt_Ca_o))/(exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00));
				_prvt_calc_i_NaK = (((((_prvt_P_NaK*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk))*__OLD_[17])/(__OLD_[17]+_prvt_K_mNa))/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.5300000000e-02*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_NaCa = ((_prvt_K_NaCa*((exp(((_prvt_gamma*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[17],3.0000000000e+00)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,3.0000000000e+00)*__OLD_[14]*_prvt_alpha)))/((pow(_prvt_Km_Nai,3.0000000000e+00)+pow(_prvt_Na_o,3.0000000000e+00))*(_prvt_Km_Ca+_prvt_Ca_o)*(1.0000000000e+00+(_prvt_K_sat*exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))))));
				_prvt_calc_i_p_Ca = ((_prvt_g_pCa*__OLD_[14])/(__OLD_[14]+_prvt_K_pCa));
				_prvt_calc_i_up = (_prvt_Vmax_up/(1.0000000000e+00+(pow(_prvt_K_up,2.0000000000e+00)/pow(__OLD_[14],2.0000000000e+00))));
				_prvt_calc_i_leak = (_prvt_V_leak*(__OLD_[15]-__OLD_[14]));
				_prvt_calc_i_xfer = (_prvt_V_xfer*(__OLD_[16]-__OLD_[14]));
				_prvt_calc_kcasr = (_prvt_max_sr-((_prvt_max_sr-_prvt_min_sr)/(1.0000000000e+00+pow((_prvt_EC/__OLD_[15]),2.0000000000e+00))));
				_prvt_calc_Ca_i_bufc = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_c*_prvt_K_buf_c)/pow((__OLD_[14]+_prvt_K_buf_c),2.0000000000e+00))));
				_prvt_calc_Ca_sr_bufsr = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_sr*_prvt_K_buf_sr)/pow((__OLD_[15]+_prvt_K_buf_sr),2.0000000000e+00))));
				_prvt_calc_Ca_ss_bufss = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_ss*_prvt_K_buf_ss)/pow((__OLD_[16]+_prvt_K_buf_ss),2.0000000000e+00))));
				_prvt_calc_alpha_K1 = (1.0000000000e-01/(1.0000000000e+00+exp((6.0000000000e-02*((__OLD_[0]-_prvt_calc_E_K)-2.0000000000e+02)))));
				_prvt_calc_beta_K1 = (((3.0000000000e+00*exp((2.0000000000e-04*((__OLD_[0]-_prvt_calc_E_K)+1.0000000000e+02))))+exp((1.0000000000e-01*((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01))))/(1.0000000000e+00+exp(((-5.0000000000e-01)*(__OLD_[0]-_prvt_calc_E_K)))));
				_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[1]*__OLD_[2]*(__OLD_[0]-_prvt_calc_E_K)*pow((_prvt_K_o/5.4000000000e+00),1.0/2.0));
				_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
				_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*__OLD_[6]*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
				_prvt_calc_i_to = (_prvt_g_to*__OLD_[12]*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_p_K = ((_prvt_g_pK*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/5.9800000000e+00))));
				_prvt_calc_xK1_inf = (_prvt_calc_alpha_K1/(_prvt_calc_alpha_K1+_prvt_calc_beta_K1));
				_prvt_calc_k1 = (_prvt_k1_prime/_prvt_calc_kcasr);
				_prvt_calc_k2 = (_prvt_k2_prime*_prvt_calc_kcasr);
				_prvt_calc_i_K1 = (_prvt_g_K1*_prvt_calc_xK1_inf*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_O = ((_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00)*__OLD_[13])/(_prvt_k3+(_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00))));
				_prvt_calc_i_rel = (_prvt_V_rel*_prvt_calc_O*(__OLD_[15]-__OLD_[16]));
				__K1_[0]= (-(_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_CaL+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_NaCa+_prvt_calc_i_b_Ca+_prvt_calc_i_p_K+_prvt_calc_i_p_Ca+_prvt_calc_i_Stim));
				__NEW_[0]= __K1_[0] * _prvt_dtime + __OLD_[0];
				__K1_[13]= (((-_prvt_calc_k2)*__OLD_[16]*__OLD_[13])+(_prvt_k4*(1.0000000000e+00-__OLD_[13])));
				__NEW_[13]= __K1_[13] * _prvt_dtime + __OLD_[13];
				__K1_[14]= (_prvt_calc_Ca_i_bufc*(((((_prvt_calc_i_leak-_prvt_calc_i_up)*_prvt_V_sr)/_prvt_V_c)+_prvt_calc_i_xfer)-((((_prvt_calc_i_b_Ca+_prvt_calc_i_p_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Cm)/(2.0000000000e+00*_prvt_V_c*_prvt_F))));
				__NEW_[14]= __K1_[14] * _prvt_dtime + __OLD_[14];
				__K1_[15]= (_prvt_calc_Ca_sr_bufsr*(_prvt_calc_i_up-(_prvt_calc_i_rel+_prvt_calc_i_leak)));
				__NEW_[15]= __K1_[15] * _prvt_dtime + __OLD_[15];
				__K1_[16]= (_prvt_calc_Ca_ss_bufss*(((((-_prvt_calc_i_CaL)*_prvt_Cm)/(2.0000000000e+00*_prvt_V_ss*_prvt_F))+((_prvt_calc_i_rel*_prvt_V_sr)/_prvt_V_ss))-((_prvt_calc_i_xfer*_prvt_V_c)/_prvt_V_ss)));
				__NEW_[16]= __K1_[16] * _prvt_dtime + __OLD_[16];
				__K1_[17]= (((-(_prvt_calc_i_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
				__NEW_[17]= __K1_[17] * _prvt_dtime + __OLD_[17];
				__K1_[18]= (((-((_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_p_K+_prvt_calc_i_Stim)-(2.0000000000e+00*_prvt_calc_i_NaK)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
				__NEW_[18]= __K1_[18] * _prvt_dtime + __OLD_[18];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[1])
			{
				_prvt_calc_xr1_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-2.6000000000e+01)-__OLD_[0])/7.0000000000e+00))));
				_prvt_calc_alpha_xr1 = (4.5000000000e+02/(1.0000000000e+00+exp((((-4.5000000000e+01)-__OLD_[0])/1.0000000000e+01))));
				_prvt_calc_beta_xr1 = (6.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.1500000000e+01))));
				_prvt_calc_tau_xr1 = (1.0000000000e+00*_prvt_calc_alpha_xr1*_prvt_calc_beta_xr1);
				__K1_[1]= ((_prvt_calc_xr1_inf-__OLD_[1])/_prvt_calc_tau_xr1);
				__NEW_[1]= __K1_[1] * _prvt_dtime + __OLD_[1];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[2])
			{
				_prvt_calc_xr2_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+8.8000000000e+01)/2.4000000000e+01))));
				_prvt_calc_alpha_xr2 = (3.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/2.0000000000e+01))));
				_prvt_calc_beta_xr2 = (1.1200000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-6.0000000000e+01)/2.0000000000e+01))));
				_prvt_calc_tau_xr2 = (1.0000000000e+00*_prvt_calc_alpha_xr2*_prvt_calc_beta_xr2);
				__K1_[2]= ((_prvt_calc_xr2_inf-__OLD_[2])/_prvt_calc_tau_xr2);
				__NEW_[2]= __K1_[2] * _prvt_dtime + __OLD_[2];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[3])
			{
				_prvt_calc_xs_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-5.0000000000e+00)-__OLD_[0])/1.4000000000e+01))));
				_prvt_calc_alpha_xs = (1.4000000000e+03/pow((1.0000000000e+00+exp(((5.0000000000e+00-__OLD_[0])/6.0000000000e+00))),1.0/2.0));
				_prvt_calc_beta_xs = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-3.5000000000e+01)/1.5000000000e+01))));
				_prvt_calc_tau_xs = ((1.0000000000e+00*_prvt_calc_alpha_xs*_prvt_calc_beta_xs)+8.0000000000e+01);
				__K1_[3]= ((_prvt_calc_xs_inf-__OLD_[3])/_prvt_calc_tau_xs);
				__NEW_[3]= __K1_[3] * _prvt_dtime + __OLD_[3];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[4])
			{
				_prvt_calc_m_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp((((-5.6860000000e+01)-__OLD_[0])/9.0300000000e+00))),2.0000000000e+00));
				_prvt_calc_alpha_m = (1.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/5.0000000000e+00))));
				_prvt_calc_beta_m = ((1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/5.0000000000e+00))))+(1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]-5.0000000000e+01)/2.0000000000e+02)))));
				_prvt_calc_tau_m = (1.0000000000e+00*_prvt_calc_alpha_m*_prvt_calc_beta_m);
				__K1_[4]= ((_prvt_calc_m_inf-__OLD_[4])/_prvt_calc_tau_m);
				__NEW_[4]= __K1_[4] * _prvt_dtime + __OLD_[4];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[5])
			{
				_prvt_calc_h_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
				_prvt_calc_alpha_h = ((__OLD_[0]<(-4.0000000000e+01)))
?((5.7000000000e-02*exp(((-(__OLD_[0]+8.0000000000e+01))/6.8000000000e+00))))
:(0.0000000000e+00);
				_prvt_calc_beta_h = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.7000000000e+00*exp((7.9000000000e-02*__OLD_[0])))+(3.1000000000e+05*exp((3.4850000000e-01*__OLD_[0])))))
:((7.7000000000e-01/(1.3000000000e-01*(1.0000000000e+00+exp(((__OLD_[0]+1.0660000000e+01)/(-1.1100000000e+01)))))));
				_prvt_calc_tau_h = (1.0000000000e+00/(_prvt_calc_alpha_h+_prvt_calc_beta_h));
				__K1_[5]= ((_prvt_calc_h_inf-__OLD_[5])/_prvt_calc_tau_h);
				__NEW_[5]= __K1_[5] * _prvt_dtime + __OLD_[5];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[6])
			{
				_prvt_calc_j_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
				_prvt_calc_alpha_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((((((-2.5428000000e+04)*exp((2.4440000000e-01*__OLD_[0])))-(6.9480000000e-06*exp(((-4.3910000000e-02)*__OLD_[0]))))*(__OLD_[0]+3.7780000000e+01))/1.0000000000e+00)/(1.0000000000e+00+exp((3.1100000000e-01*(__OLD_[0]+7.9230000000e+01))))))
:(0.0000000000e+00);
				_prvt_calc_beta_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.4240000000e-02*exp(((-1.0520000000e-02)*__OLD_[0])))/(1.0000000000e+00+exp(((-1.3780000000e-01)*(__OLD_[0]+4.0140000000e+01))))))
:(((6.0000000000e-01*exp((5.7000000000e-02*__OLD_[0])))/(1.0000000000e+00+exp(((-1.0000000000e-01)*(__OLD_[0]+3.2000000000e+01))))));
				_prvt_calc_tau_j = (1.0000000000e+00/(_prvt_calc_alpha_j+_prvt_calc_beta_j));
				__K1_[6]= ((_prvt_calc_j_inf-__OLD_[6])/_prvt_calc_tau_j);
				__NEW_[6]= __K1_[6] * _prvt_dtime + __OLD_[6];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[7])
			{
				_prvt_calc_d_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-8.0000000000e+00)-__OLD_[0])/7.5000000000e+00))));
				_prvt_calc_alpha_d = ((1.4000000000e+00/(1.0000000000e+00+exp((((-3.5000000000e+01)-__OLD_[0])/1.3000000000e+01))))+2.5000000000e-01);
				_prvt_calc_beta_d = (1.4000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.0000000000e+00)/5.0000000000e+00))));
				_prvt_calc_gamma_d = (1.0000000000e+00/(1.0000000000e+00+exp(((5.0000000000e+01-__OLD_[0])/2.0000000000e+01))));
				_prvt_calc_tau_d = ((1.0000000000e+00*_prvt_calc_alpha_d*_prvt_calc_beta_d)+_prvt_calc_gamma_d);
				__K1_[7]= ((_prvt_calc_d_inf-__OLD_[7])/_prvt_calc_tau_d);
				__NEW_[7]= __K1_[7] * _prvt_dtime + __OLD_[7];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[8])
			{
				_prvt_calc_f_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/7.0000000000e+00))));
				_prvt_calc_tau_f = ((1.1025000000e+03*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.2500000000e+02)))+(2.0000000000e+02/(1.0000000000e+00+exp(((1.3000000000e+01-__OLD_[0])/1.0000000000e+01))))+(1.8000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01))))+2.0000000000e+01);
				__K1_[8]= ((_prvt_calc_f_inf-__OLD_[8])/_prvt_calc_tau_f);
				__NEW_[8]= __K1_[8] * _prvt_dtime + __OLD_[8];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[9])
			{
				_prvt_calc_f2_inf = ((6.7000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/7.0000000000e+00))))+3.3000000000e-01);
				_prvt_calc_tau_f2 = ((5.6200000000e+02*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.4000000000e+02)))+(3.1000000000e+01/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/1.0000000000e+01))))+(8.0000000000e+01/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01)))));
				__K1_[9]= ((_prvt_calc_f2_inf-__OLD_[9])/_prvt_calc_tau_f2);
				__NEW_[9]= __K1_[9] * _prvt_dtime + __OLD_[9];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[10])
			{
				_prvt_calc_fCass_inf = ((6.0000000000e-01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+4.0000000000e-01);
				_prvt_calc_tau_fCass = ((8.0000000000e+01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+2.0000000000e+00);
				__K1_[10]= ((_prvt_calc_fCass_inf-__OLD_[10])/_prvt_calc_tau_fCass);
				__NEW_[10]= __K1_[10] * _prvt_dtime + __OLD_[10];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[11])
			{
				_prvt_calc_s_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/5.0000000000e+00))));
				_prvt_calc_tau_s = ((8.5000000000e+01*exp(((-pow((__OLD_[0]+4.5000000000e+01),2.0000000000e+00))/3.2000000000e+02)))+(5.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-2.0000000000e+01)/5.0000000000e+00))))+3.0000000000e+00);
				__K1_[11]= ((_prvt_calc_s_inf-__OLD_[11])/_prvt_calc_tau_s);
				__NEW_[11]= __K1_[11] * _prvt_dtime + __OLD_[11];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[12])
			{
				_prvt_calc_r_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((2.0000000000e+01-__OLD_[0])/6.0000000000e+00))));
				_prvt_calc_tau_r = ((9.5000000000e+00*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.8000000000e+03)))+8.0000000000e-01);
				__K1_[12]= ((_prvt_calc_r_inf-__OLD_[12])/_prvt_calc_tau_r);
				__NEW_[12]= __K1_[12] * _prvt_dtime + __OLD_[12];
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
					_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[17])));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[18])));
					_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[18]+(_prvt_P_kna*__OLD_[17])))));
					_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[14])));
					_prvt_calc_i_CaL = ((((_prvt_g_CaL*__OLD_[7]*__OLD_[8]*__OLD_[9]*__OLD_[10]*4.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*((2.5000000000e-01*__OLD_[16]*exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))))-_prvt_Ca_o))/(exp(((2.0000000000e+00*(__OLD_[0]-1.5000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00));
					_prvt_calc_i_NaK = (((((_prvt_P_NaK*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk))*__OLD_[17])/(__OLD_[17]+_prvt_K_mNa))/(1.0000000000e+00+(1.2450000000e-01*exp((((-1.0000000000e-01)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T))))+(3.5300000000e-02*exp((((-__OLD_[0])*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_NaCa = ((_prvt_K_NaCa*((exp(((_prvt_gamma*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[17],3.0000000000e+00)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,3.0000000000e+00)*__OLD_[14]*_prvt_alpha)))/((pow(_prvt_Km_Nai,3.0000000000e+00)+pow(_prvt_Na_o,3.0000000000e+00))*(_prvt_Km_Ca+_prvt_Ca_o)*(1.0000000000e+00+(_prvt_K_sat*exp((((_prvt_gamma-1.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))))));
					_prvt_calc_i_p_Ca = ((_prvt_g_pCa*__OLD_[14])/(__OLD_[14]+_prvt_K_pCa));
					_prvt_calc_i_up = (_prvt_Vmax_up/(1.0000000000e+00+(pow(_prvt_K_up,2.0000000000e+00)/pow(__OLD_[14],2.0000000000e+00))));
					_prvt_calc_i_leak = (_prvt_V_leak*(__OLD_[15]-__OLD_[14]));
					_prvt_calc_i_xfer = (_prvt_V_xfer*(__OLD_[16]-__OLD_[14]));
					_prvt_calc_kcasr = (_prvt_max_sr-((_prvt_max_sr-_prvt_min_sr)/(1.0000000000e+00+pow((_prvt_EC/__OLD_[15]),2.0000000000e+00))));
					_prvt_calc_Ca_i_bufc = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_c*_prvt_K_buf_c)/pow((__OLD_[14]+_prvt_K_buf_c),2.0000000000e+00))));
					_prvt_calc_Ca_sr_bufsr = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_sr*_prvt_K_buf_sr)/pow((__OLD_[15]+_prvt_K_buf_sr),2.0000000000e+00))));
					_prvt_calc_Ca_ss_bufss = (1.0000000000e+00/(1.0000000000e+00+((_prvt_Buf_ss*_prvt_K_buf_ss)/pow((__OLD_[16]+_prvt_K_buf_ss),2.0000000000e+00))));
					_prvt_calc_alpha_K1 = (1.0000000000e-01/(1.0000000000e+00+exp((6.0000000000e-02*((__OLD_[0]-_prvt_calc_E_K)-2.0000000000e+02)))));
					_prvt_calc_beta_K1 = (((3.0000000000e+00*exp((2.0000000000e-04*((__OLD_[0]-_prvt_calc_E_K)+1.0000000000e+02))))+exp((1.0000000000e-01*((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01))))/(1.0000000000e+00+exp(((-5.0000000000e-01)*(__OLD_[0]-_prvt_calc_E_K)))));
					_prvt_calc_i_Kr = (_prvt_g_Kr*__OLD_[1]*__OLD_[2]*(__OLD_[0]-_prvt_calc_E_K)*pow((_prvt_K_o/5.4000000000e+00),1.0/2.0));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
					_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*__OLD_[6]*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_to = (_prvt_g_to*__OLD_[12]*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_p_K = ((_prvt_g_pK*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/5.9800000000e+00))));
					_prvt_calc_xK1_inf = (_prvt_calc_alpha_K1/(_prvt_calc_alpha_K1+_prvt_calc_beta_K1));
					_prvt_calc_k1 = (_prvt_k1_prime/_prvt_calc_kcasr);
					_prvt_calc_k2 = (_prvt_k2_prime*_prvt_calc_kcasr);
					_prvt_calc_i_K1 = (_prvt_g_K1*_prvt_calc_xK1_inf*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_O = ((_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00)*__OLD_[13])/(_prvt_k3+(_prvt_calc_k1*pow(__OLD_[16],2.0000000000e+00))));
					_prvt_calc_i_rel = (_prvt_V_rel*_prvt_calc_O*(__OLD_[15]-__OLD_[16]));
					__K2_[0]= (-(_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_CaL+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_NaCa+_prvt_calc_i_b_Ca+_prvt_calc_i_p_K+_prvt_calc_i_p_Ca+_prvt_calc_i_Stim));
					_prvt_aux_tol = fabs(__OLD_[0])*_prvt_rel_tol_;
					__TOL_[0] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[0] = fabs((_prvt_dtime/2) * (__K1_[0] - __K2_[0])/__TOL_[0]);
					__K2_[13]= (((-_prvt_calc_k2)*__OLD_[16]*__OLD_[13])+(_prvt_k4*(1.0000000000e+00-__OLD_[13])));
					_prvt_aux_tol = fabs(__OLD_[13])*_prvt_rel_tol_;
					__TOL_[13] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[13] = fabs((_prvt_dtime/2) * (__K1_[13] - __K2_[13])/__TOL_[13]);
					__K2_[14]= (_prvt_calc_Ca_i_bufc*(((((_prvt_calc_i_leak-_prvt_calc_i_up)*_prvt_V_sr)/_prvt_V_c)+_prvt_calc_i_xfer)-((((_prvt_calc_i_b_Ca+_prvt_calc_i_p_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa))*_prvt_Cm)/(2.0000000000e+00*_prvt_V_c*_prvt_F))));
					_prvt_aux_tol = fabs(__OLD_[14])*_prvt_rel_tol_;
					__TOL_[14] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[14] = fabs((_prvt_dtime/2) * (__K1_[14] - __K2_[14])/__TOL_[14]);
					__K2_[15]= (_prvt_calc_Ca_sr_bufsr*(_prvt_calc_i_up-(_prvt_calc_i_rel+_prvt_calc_i_leak)));
					_prvt_aux_tol = fabs(__OLD_[15])*_prvt_rel_tol_;
					__TOL_[15] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[15] = fabs((_prvt_dtime/2) * (__K1_[15] - __K2_[15])/__TOL_[15]);
					__K2_[16]= (_prvt_calc_Ca_ss_bufss*(((((-_prvt_calc_i_CaL)*_prvt_Cm)/(2.0000000000e+00*_prvt_V_ss*_prvt_F))+((_prvt_calc_i_rel*_prvt_V_sr)/_prvt_V_ss))-((_prvt_calc_i_xfer*_prvt_V_c)/_prvt_V_ss)));
					_prvt_aux_tol = fabs(__OLD_[16])*_prvt_rel_tol_;
					__TOL_[16] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[16] = fabs((_prvt_dtime/2) * (__K1_[16] - __K2_[16])/__TOL_[16]);
					__K2_[17]= (((-(_prvt_calc_i_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
					_prvt_aux_tol = fabs(__OLD_[17])*_prvt_rel_tol_;
					__TOL_[17] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[17] = fabs((_prvt_dtime/2) * (__K1_[17] - __K2_[17])/__TOL_[17]);
					__K2_[18]= (((-((_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_p_K+_prvt_calc_i_Stim)-(2.0000000000e+00*_prvt_calc_i_NaK)))/(_prvt_V_c*_prvt_F))*_prvt_Cm);
					_prvt_aux_tol = fabs(__OLD_[18])*_prvt_rel_tol_;
					__TOL_[18] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[18] = fabs((_prvt_dtime/2) * (__K1_[18] - __K2_[18])/__TOL_[18]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[1])
				{
					_prvt_calc_xr1_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-2.6000000000e+01)-__OLD_[0])/7.0000000000e+00))));
					_prvt_calc_alpha_xr1 = (4.5000000000e+02/(1.0000000000e+00+exp((((-4.5000000000e+01)-__OLD_[0])/1.0000000000e+01))));
					_prvt_calc_beta_xr1 = (6.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.1500000000e+01))));
					_prvt_calc_tau_xr1 = (1.0000000000e+00*_prvt_calc_alpha_xr1*_prvt_calc_beta_xr1);
					__K2_[1]= ((_prvt_calc_xr1_inf-__OLD_[1])/_prvt_calc_tau_xr1);
					_prvt_aux_tol = fabs(__OLD_[1])*_prvt_rel_tol_;
					__TOL_[1] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[1] = fabs((_prvt_dtime/2) * (__K1_[1] - __K2_[1])/__TOL_[1]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[2])
				{
					_prvt_calc_xr2_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+8.8000000000e+01)/2.4000000000e+01))));
					_prvt_calc_alpha_xr2 = (3.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/2.0000000000e+01))));
					_prvt_calc_beta_xr2 = (1.1200000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-6.0000000000e+01)/2.0000000000e+01))));
					_prvt_calc_tau_xr2 = (1.0000000000e+00*_prvt_calc_alpha_xr2*_prvt_calc_beta_xr2);
					__K2_[2]= ((_prvt_calc_xr2_inf-__OLD_[2])/_prvt_calc_tau_xr2);
					_prvt_aux_tol = fabs(__OLD_[2])*_prvt_rel_tol_;
					__TOL_[2] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[2] = fabs((_prvt_dtime/2) * (__K1_[2] - __K2_[2])/__TOL_[2]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[3])
				{
					_prvt_calc_xs_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-5.0000000000e+00)-__OLD_[0])/1.4000000000e+01))));
					_prvt_calc_alpha_xs = (1.4000000000e+03/pow((1.0000000000e+00+exp(((5.0000000000e+00-__OLD_[0])/6.0000000000e+00))),1.0/2.0));
					_prvt_calc_beta_xs = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-3.5000000000e+01)/1.5000000000e+01))));
					_prvt_calc_tau_xs = ((1.0000000000e+00*_prvt_calc_alpha_xs*_prvt_calc_beta_xs)+8.0000000000e+01);
					__K2_[3]= ((_prvt_calc_xs_inf-__OLD_[3])/_prvt_calc_tau_xs);
					_prvt_aux_tol = fabs(__OLD_[3])*_prvt_rel_tol_;
					__TOL_[3] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[3] = fabs((_prvt_dtime/2) * (__K1_[3] - __K2_[3])/__TOL_[3]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[4])
				{
					_prvt_calc_m_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp((((-5.6860000000e+01)-__OLD_[0])/9.0300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_m = (1.0000000000e+00/(1.0000000000e+00+exp((((-6.0000000000e+01)-__OLD_[0])/5.0000000000e+00))));
					_prvt_calc_beta_m = ((1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/5.0000000000e+00))))+(1.0000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]-5.0000000000e+01)/2.0000000000e+02)))));
					_prvt_calc_tau_m = (1.0000000000e+00*_prvt_calc_alpha_m*_prvt_calc_beta_m);
					__K2_[4]= ((_prvt_calc_m_inf-__OLD_[4])/_prvt_calc_tau_m);
					_prvt_aux_tol = fabs(__OLD_[4])*_prvt_rel_tol_;
					__TOL_[4] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[4] = fabs((_prvt_dtime/2) * (__K1_[4] - __K2_[4])/__TOL_[4]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[5])
				{
					_prvt_calc_h_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_h = ((__OLD_[0]<(-4.0000000000e+01)))
?((5.7000000000e-02*exp(((-(__OLD_[0]+8.0000000000e+01))/6.8000000000e+00))))
:(0.0000000000e+00);
					_prvt_calc_beta_h = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.7000000000e+00*exp((7.9000000000e-02*__OLD_[0])))+(3.1000000000e+05*exp((3.4850000000e-01*__OLD_[0])))))
:((7.7000000000e-01/(1.3000000000e-01*(1.0000000000e+00+exp(((__OLD_[0]+1.0660000000e+01)/(-1.1100000000e+01)))))));
					_prvt_calc_tau_h = (1.0000000000e+00/(_prvt_calc_alpha_h+_prvt_calc_beta_h));
					__K2_[5]= ((_prvt_calc_h_inf-__OLD_[5])/_prvt_calc_tau_h);
					_prvt_aux_tol = fabs(__OLD_[5])*_prvt_rel_tol_;
					__TOL_[5] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[5] = fabs((_prvt_dtime/2) * (__K1_[5] - __K2_[5])/__TOL_[5]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[6])
				{
					_prvt_calc_j_inf = (1.0000000000e+00/pow((1.0000000000e+00+exp(((__OLD_[0]+7.1550000000e+01)/7.4300000000e+00))),2.0000000000e+00));
					_prvt_calc_alpha_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((((((-2.5428000000e+04)*exp((2.4440000000e-01*__OLD_[0])))-(6.9480000000e-06*exp(((-4.3910000000e-02)*__OLD_[0]))))*(__OLD_[0]+3.7780000000e+01))/1.0000000000e+00)/(1.0000000000e+00+exp((3.1100000000e-01*(__OLD_[0]+7.9230000000e+01))))))
:(0.0000000000e+00);
					_prvt_calc_beta_j = ((__OLD_[0]<(-4.0000000000e+01)))
?(((2.4240000000e-02*exp(((-1.0520000000e-02)*__OLD_[0])))/(1.0000000000e+00+exp(((-1.3780000000e-01)*(__OLD_[0]+4.0140000000e+01))))))
:(((6.0000000000e-01*exp((5.7000000000e-02*__OLD_[0])))/(1.0000000000e+00+exp(((-1.0000000000e-01)*(__OLD_[0]+3.2000000000e+01))))));
					_prvt_calc_tau_j = (1.0000000000e+00/(_prvt_calc_alpha_j+_prvt_calc_beta_j));
					__K2_[6]= ((_prvt_calc_j_inf-__OLD_[6])/_prvt_calc_tau_j);
					_prvt_aux_tol = fabs(__OLD_[6])*_prvt_rel_tol_;
					__TOL_[6] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[6] = fabs((_prvt_dtime/2) * (__K1_[6] - __K2_[6])/__TOL_[6]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[7])
				{
					_prvt_calc_d_inf = (1.0000000000e+00/(1.0000000000e+00+exp((((-8.0000000000e+00)-__OLD_[0])/7.5000000000e+00))));
					_prvt_calc_alpha_d = ((1.4000000000e+00/(1.0000000000e+00+exp((((-3.5000000000e+01)-__OLD_[0])/1.3000000000e+01))))+2.5000000000e-01);
					_prvt_calc_beta_d = (1.4000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.0000000000e+00)/5.0000000000e+00))));
					_prvt_calc_gamma_d = (1.0000000000e+00/(1.0000000000e+00+exp(((5.0000000000e+01-__OLD_[0])/2.0000000000e+01))));
					_prvt_calc_tau_d = ((1.0000000000e+00*_prvt_calc_alpha_d*_prvt_calc_beta_d)+_prvt_calc_gamma_d);
					__K2_[7]= ((_prvt_calc_d_inf-__OLD_[7])/_prvt_calc_tau_d);
					_prvt_aux_tol = fabs(__OLD_[7])*_prvt_rel_tol_;
					__TOL_[7] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[7] = fabs((_prvt_dtime/2) * (__K1_[7] - __K2_[7])/__TOL_[7]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[8])
				{
					_prvt_calc_f_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/7.0000000000e+00))));
					_prvt_calc_tau_f = ((1.1025000000e+03*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.2500000000e+02)))+(2.0000000000e+02/(1.0000000000e+00+exp(((1.3000000000e+01-__OLD_[0])/1.0000000000e+01))))+(1.8000000000e+02/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01))))+2.0000000000e+01);
					__K2_[8]= ((_prvt_calc_f_inf-__OLD_[8])/_prvt_calc_tau_f);
					_prvt_aux_tol = fabs(__OLD_[8])*_prvt_rel_tol_;
					__TOL_[8] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[8] = fabs((_prvt_dtime/2) * (__K1_[8] - __K2_[8])/__TOL_[8]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[9])
				{
					_prvt_calc_f2_inf = ((6.7000000000e-01/(1.0000000000e+00+exp(((__OLD_[0]+3.5000000000e+01)/7.0000000000e+00))))+3.3000000000e-01);
					_prvt_calc_tau_f2 = ((5.6200000000e+02*exp(((-pow((__OLD_[0]+2.7000000000e+01),2.0000000000e+00))/2.4000000000e+02)))+(3.1000000000e+01/(1.0000000000e+00+exp(((2.5000000000e+01-__OLD_[0])/1.0000000000e+01))))+(8.0000000000e+01/(1.0000000000e+00+exp(((__OLD_[0]+3.0000000000e+01)/1.0000000000e+01)))));
					__K2_[9]= ((_prvt_calc_f2_inf-__OLD_[9])/_prvt_calc_tau_f2);
					_prvt_aux_tol = fabs(__OLD_[9])*_prvt_rel_tol_;
					__TOL_[9] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[9] = fabs((_prvt_dtime/2) * (__K1_[9] - __K2_[9])/__TOL_[9]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[10])
				{
					_prvt_calc_fCass_inf = ((6.0000000000e-01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+4.0000000000e-01);
					_prvt_calc_tau_fCass = ((8.0000000000e+01/(1.0000000000e+00+pow((__OLD_[16]/5.0000000000e-02),2.0000000000e+00)))+2.0000000000e+00);
					__K2_[10]= ((_prvt_calc_fCass_inf-__OLD_[10])/_prvt_calc_tau_fCass);
					_prvt_aux_tol = fabs(__OLD_[10])*_prvt_rel_tol_;
					__TOL_[10] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[10] = fabs((_prvt_dtime/2) * (__K1_[10] - __K2_[10])/__TOL_[10]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[11])
				{
					_prvt_calc_s_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+2.0000000000e+01)/5.0000000000e+00))));
					_prvt_calc_tau_s = ((8.5000000000e+01*exp(((-pow((__OLD_[0]+4.5000000000e+01),2.0000000000e+00))/3.2000000000e+02)))+(5.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]-2.0000000000e+01)/5.0000000000e+00))))+3.0000000000e+00);
					__K2_[11]= ((_prvt_calc_s_inf-__OLD_[11])/_prvt_calc_tau_s);
					_prvt_aux_tol = fabs(__OLD_[11])*_prvt_rel_tol_;
					__TOL_[11] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[11] = fabs((_prvt_dtime/2) * (__K1_[11] - __K2_[11])/__TOL_[11]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[12])
				{
					_prvt_calc_r_inf = (1.0000000000e+00/(1.0000000000e+00+exp(((2.0000000000e+01-__OLD_[0])/6.0000000000e+00))));
					_prvt_calc_tau_r = ((9.5000000000e+00*exp(((-pow((__OLD_[0]+4.0000000000e+01),2.0000000000e+00))/1.8000000000e+03)))+8.0000000000e-01);
					__K2_[12]= ((_prvt_calc_r_inf-__OLD_[12])/_prvt_calc_tau_r);
					_prvt_aux_tol = fabs(__OLD_[12])*_prvt_rel_tol_;
					__TOL_[12] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[12] = fabs((_prvt_dtime/2) * (__K1_[12] - __K2_[12])/__TOL_[12]);
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
						__NEW_[13] = __K1_[13] * _prvt_dtime + __OLD_AUX_[13];
						__NEW_[14] = __K1_[14] * _prvt_dtime + __OLD_AUX_[14];
						__NEW_[15] = __K1_[15] * _prvt_dtime + __OLD_AUX_[15];
						__NEW_[16] = __K1_[16] * _prvt_dtime + __OLD_AUX_[16];
						__NEW_[17] = __K1_[17] * _prvt_dtime + __OLD_AUX_[17];
						__NEW_[18] = __K1_[18] * _prvt_dtime + __OLD_AUX_[18];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[1] = __K1_[1] * _prvt_dtime + __OLD_AUX_[1];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[2] = __K1_[2] * _prvt_dtime + __OLD_AUX_[2];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[3] = __K1_[3] * _prvt_dtime + __OLD_AUX_[3];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[4] = __K1_[4] * _prvt_dtime + __OLD_AUX_[4];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[5] = __K1_[5] * _prvt_dtime + __OLD_AUX_[5];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[6] = __K1_[6] * _prvt_dtime + __OLD_AUX_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[7] = __K1_[7] * _prvt_dtime + __OLD_AUX_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[8] = __K1_[8] * _prvt_dtime + __OLD_AUX_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[9] = __K1_[9] * _prvt_dtime + __OLD_AUX_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[10] = __K1_[10] * _prvt_dtime + __OLD_AUX_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[11] = __K1_[11] * _prvt_dtime + __OLD_AUX_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[12] = __K1_[12] * _prvt_dtime + __OLD_AUX_[12];
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
							this->Xr1_old_ = __OLD_AUX_[1];
							this->Xr1_new_ = __OLD_[1];
							this->Xr2_old_ = __OLD_AUX_[2];
							this->Xr2_new_ = __OLD_[2];
							this->Xs_old_ = __OLD_AUX_[3];
							this->Xs_new_ = __OLD_[3];
							this->m_old_ = __OLD_AUX_[4];
							this->m_new_ = __OLD_[4];
							this->h_old_ = __OLD_AUX_[5];
							this->h_new_ = __OLD_[5];
							this->j_old_ = __OLD_AUX_[6];
							this->j_new_ = __OLD_[6];
							this->d_old_ = __OLD_AUX_[7];
							this->d_new_ = __OLD_[7];
							this->f_old_ = __OLD_AUX_[8];
							this->f_new_ = __OLD_[8];
							this->f2_old_ = __OLD_AUX_[9];
							this->f2_new_ = __OLD_[9];
							this->fCass_old_ = __OLD_AUX_[10];
							this->fCass_new_ = __OLD_[10];
							this->s_old_ = __OLD_AUX_[11];
							this->s_new_ = __OLD_[11];
							this->r_old_ = __OLD_AUX_[12];
							this->r_new_ = __OLD_[12];
							this->R_prime_old_ = __OLD_AUX_[13];
							this->R_prime_new_ = __OLD_[13];
							this->Ca_i_old_ = __OLD_AUX_[14];
							this->Ca_i_new_ = __OLD_[14];
							this->Ca_SR_old_ = __OLD_AUX_[15];
							this->Ca_SR_new_ = __OLD_[15];
							this->Ca_ss_old_ = __OLD_AUX_[16];
							this->Ca_ss_new_ = __OLD_[16];
							this->Na_i_old_ = __OLD_AUX_[17];
							this->Na_i_new_ = __OLD_[17];
							this->K_i_old_ = __OLD_AUX_[18];
							this->K_i_new_ = __OLD_[18];
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
						__NEW_[13] = __K2_[13] * _prvt_dtime + __OLD_[13];
						__NEW_[14] = __K2_[14] * _prvt_dtime + __OLD_[14];
						__NEW_[15] = __K2_[15] * _prvt_dtime + __OLD_[15];
						__NEW_[16] = __K2_[16] * _prvt_dtime + __OLD_[16];
						__NEW_[17] = __K2_[17] * _prvt_dtime + __OLD_[17];
						__NEW_[18] = __K2_[18] * _prvt_dtime + __OLD_[18];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[1] = __K2_[1] * _prvt_dtime + __OLD_[1];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[2] = __K2_[2] * _prvt_dtime + __OLD_[2];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[3] = __K2_[3] * _prvt_dtime + __OLD_[3];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[4] = __K2_[4] * _prvt_dtime + __OLD_[4];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[5] = __K2_[5] * _prvt_dtime + __OLD_[5];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[6] = __K2_[6] * _prvt_dtime + __OLD_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[7] = __K2_[7] * _prvt_dtime + __OLD_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[8] = __K2_[8] * _prvt_dtime + __OLD_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[9] = __K2_[9] * _prvt_dtime + __OLD_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[10] = __K2_[10] * _prvt_dtime + __OLD_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[11] = __K2_[11] * _prvt_dtime + __OLD_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[12] = __K2_[12] * _prvt_dtime + __OLD_[12];
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
		case 1:		return Xr1_old_;    break;
		case 2:		return Xr2_old_;    break;
		case 3:		return Xs_old_;    break;
		case 4:		return m_old_;    break;
		case 5:		return h_old_;    break;
		case 6:		return j_old_;    break;
		case 7:		return d_old_;    break;
		case 8:		return f_old_;    break;
		case 9:		return f2_old_;    break;
		case 10:		return fCass_old_;    break;
		case 11:		return s_old_;    break;
		case 12:		return r_old_;    break;
		case 13:		return R_prime_old_;    break;
		case 14:		return Ca_i_old_;    break;
		case 15:		return Ca_SR_old_;    break;
		case 16:		return Ca_ss_old_;    break;
		case 17:		return Na_i_old_;    break;
		case 18:		return K_i_old_;    break;
		default:	return 1;    break;
		}
	}

	double Solveode::getLadoDireito(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return V_lado_direito_;    break;
		case 1:		return Xr1_lado_direito_;    break;
		case 2:		return Xr2_lado_direito_;    break;
		case 3:		return Xs_lado_direito_;    break;
		case 4:		return m_lado_direito_;    break;
		case 5:		return h_lado_direito_;    break;
		case 6:		return j_lado_direito_;    break;
		case 7:		return d_lado_direito_;    break;
		case 8:		return f_lado_direito_;    break;
		case 9:		return f2_lado_direito_;    break;
		case 10:		return fCass_lado_direito_;    break;
		case 11:		return s_lado_direito_;    break;
		case 12:		return r_lado_direito_;    break;
		case 13:		return R_prime_lado_direito_;    break;
		case 14:		return Ca_i_lado_direito_;    break;
		case 15:		return Ca_SR_lado_direito_;    break;
		case 16:		return Ca_ss_lado_direito_;    break;
		case 17:		return Na_i_lado_direito_;    break;
		case 18:		return K_i_lado_direito_;    break;
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
		case 6:		return R;    break;
		case 7:		return T;    break;
		case 8:		return F;    break;
		case 9:		return Na_o;    break;
		case 10:		return K_o;    break;
		case 11:		return P_kna;    break;
		case 12:		return Ca_o;    break;
		case 13:		return g_K1;    break;
		case 14:		return g_Kr;    break;
		case 15:		return g_Ks;    break;
		case 16:		return g_Na;    break;
		case 17:		return g_bna;    break;
		case 18:		return g_CaL;    break;
		case 19:		return g_bca;    break;
		case 20:		return g_to;    break;
		case 21:		return P_NaK;    break;
		case 22:		return K_mk;    break;
		case 23:		return K_mNa;    break;
		case 24:		return K_NaCa;    break;
		case 25:		return gamma;    break;
		case 26:		return alpha;    break;
		case 27:		return Km_Nai;    break;
		case 28:		return Km_Ca;    break;
		case 29:		return K_sat;    break;
		case 30:		return g_pCa;    break;
		case 31:		return K_pCa;    break;
		case 32:		return g_pK;    break;
		case 33:		return V_rel;    break;
		case 34:		return Vmax_up;    break;
		case 35:		return K_up;    break;
		case 36:		return V_leak;    break;
		case 37:		return V_xfer;    break;
		case 38:		return k3;    break;
		case 39:		return k4;    break;
		case 40:		return k1_prime;    break;
		case 41:		return k2_prime;    break;
		case 42:		return max_sr;    break;
		case 43:		return min_sr;    break;
		case 44:		return EC;    break;
		case 45:		return Buf_c;    break;
		case 46:		return K_buf_c;    break;
		case 47:		return Buf_sr;    break;
		case 48:		return K_buf_sr;    break;
		case 49:		return Buf_ss;    break;
		case 50:		return K_buf_ss;    break;
		case 51:		return V_sr;    break;
		case 52:		return V_c;    break;
		case 53:		return Cm;    break;
		case 54:		return V_ss;    break;
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
		Variables v("|V#|Xr1#|Xr2#|Xs#|m#|h#|j#|d#|f#|f2#|fCass#|s#|r#|R_prime#|Ca_i#|Ca_SR#|Ca_ss#|Na_i#|K_i#");
		return v;
	}
	Variables Solveode::get_Parameters()
	{
		Variables v("|time#|stim_amplitude#|stim_start#|stim_end#|stim_period#|stim_duration#|R#|T#|F#|Na_o#|K_o#|P_kna#|Ca_o#|g_K1#|g_Kr#|g_Ks#|g_Na#|g_bna#|g_CaL#|g_bca#|g_to#|P_NaK#|K_mk#|K_mNa#|K_NaCa#|gamma#|alpha#|Km_Nai#|Km_Ca#|K_sat#|g_pCa#|K_pCa#|g_pK#|V_rel#|Vmax_up#|K_up#|V_leak#|V_xfer#|k3#|k4#|k1_prime#|k2_prime#|max_sr#|min_sr#|EC#|Buf_c#|K_buf_c#|Buf_sr#|K_buf_sr#|Buf_ss#|K_buf_ss#|V_sr#|V_c#|Cm#|V_ss#");
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
		Xr1_old_ = Xr1_ini_;
		if(Xr1 != NULL)free( Xr1);
			Xr1 = (double *)malloc(sizeof(double)*num_results__);
		Xr2_old_ = Xr2_ini_;
		if(Xr2 != NULL)free( Xr2);
			Xr2 = (double *)malloc(sizeof(double)*num_results__);
		Xs_old_ = Xs_ini_;
		if(Xs != NULL)free( Xs);
			Xs = (double *)malloc(sizeof(double)*num_results__);
		m_old_ = m_ini_;
		if(m != NULL)free( m);
			m = (double *)malloc(sizeof(double)*num_results__);
		h_old_ = h_ini_;
		if(h != NULL)free( h);
			h = (double *)malloc(sizeof(double)*num_results__);
		j_old_ = j_ini_;
		if(j != NULL)free( j);
			j = (double *)malloc(sizeof(double)*num_results__);
		d_old_ = d_ini_;
		if(d != NULL)free( d);
			d = (double *)malloc(sizeof(double)*num_results__);
		f_old_ = f_ini_;
		if(f != NULL)free( f);
			f = (double *)malloc(sizeof(double)*num_results__);
		f2_old_ = f2_ini_;
		if(f2 != NULL)free( f2);
			f2 = (double *)malloc(sizeof(double)*num_results__);
		fCass_old_ = fCass_ini_;
		if(fCass != NULL)free( fCass);
			fCass = (double *)malloc(sizeof(double)*num_results__);
		s_old_ = s_ini_;
		if(s != NULL)free( s);
			s = (double *)malloc(sizeof(double)*num_results__);
		r_old_ = r_ini_;
		if(r != NULL)free( r);
			r = (double *)malloc(sizeof(double)*num_results__);
		R_prime_old_ = R_prime_ini_;
		if(R_prime != NULL)free( R_prime);
			R_prime = (double *)malloc(sizeof(double)*num_results__);
		Ca_i_old_ = Ca_i_ini_;
		if(Ca_i != NULL)free( Ca_i);
			Ca_i = (double *)malloc(sizeof(double)*num_results__);
		Ca_SR_old_ = Ca_SR_ini_;
		if(Ca_SR != NULL)free( Ca_SR);
			Ca_SR = (double *)malloc(sizeof(double)*num_results__);
		Ca_ss_old_ = Ca_ss_ini_;
		if(Ca_ss != NULL)free( Ca_ss);
			Ca_ss = (double *)malloc(sizeof(double)*num_results__);
		Na_i_old_ = Na_i_ini_;
		if(Na_i != NULL)free( Na_i);
			Na_i = (double *)malloc(sizeof(double)*num_results__);
		K_i_old_ = K_i_ini_;
		if(K_i != NULL)free( K_i);
			K_i = (double *)malloc(sizeof(double)*num_results__);
		this->timeSaving = dtime;

		double diff=0;
		int counter=0;
		for (int i = 0; i< iterations;i++ )
		{
			this->time_new += dtime;

			rightHandSideFunction.function(this);
			this->V_new_ = this->V_old_ + this->V_lado_direito_ * this->dtime;
			this->Xr1_new_ = this->Xr1_old_ + this->Xr1_lado_direito_ * this->dtime;
			this->Xr2_new_ = this->Xr2_old_ + this->Xr2_lado_direito_ * this->dtime;
			this->Xs_new_ = this->Xs_old_ + this->Xs_lado_direito_ * this->dtime;
			this->m_new_ = this->m_old_ + this->m_lado_direito_ * this->dtime;
			this->h_new_ = this->h_old_ + this->h_lado_direito_ * this->dtime;
			this->j_new_ = this->j_old_ + this->j_lado_direito_ * this->dtime;
			this->d_new_ = this->d_old_ + this->d_lado_direito_ * this->dtime;
			this->f_new_ = this->f_old_ + this->f_lado_direito_ * this->dtime;
			this->f2_new_ = this->f2_old_ + this->f2_lado_direito_ * this->dtime;
			this->fCass_new_ = this->fCass_old_ + this->fCass_lado_direito_ * this->dtime;
			this->s_new_ = this->s_old_ + this->s_lado_direito_ * this->dtime;
			this->r_new_ = this->r_old_ + this->r_lado_direito_ * this->dtime;
			this->R_prime_new_ = this->R_prime_old_ + this->R_prime_lado_direito_ * this->dtime;
			this->Ca_i_new_ = this->Ca_i_old_ + this->Ca_i_lado_direito_ * this->dtime;
			this->Ca_SR_new_ = this->Ca_SR_old_ + this->Ca_SR_lado_direito_ * this->dtime;
			this->Ca_ss_new_ = this->Ca_ss_old_ + this->Ca_ss_lado_direito_ * this->dtime;
			this->Na_i_new_ = this->Na_i_old_ + this->Na_i_lado_direito_ * this->dtime;
			this->K_i_new_ = this->K_i_old_ + this->K_i_lado_direito_ * this->dtime;
			diff =  _agos_round(this->time_new - timeSaving, 5);
			if(diff==0){
				this->timeSaving += svRate;
				time_vec__[counter] = this->time_new;
				V[counter] = this->V_new_;
				Xr1[counter] = this->Xr1_new_;
				Xr2[counter] = this->Xr2_new_;
				Xs[counter] = this->Xs_new_;
				m[counter] = this->m_new_;
				h[counter] = this->h_new_;
				j[counter] = this->j_new_;
				d[counter] = this->d_new_;
				f[counter] = this->f_new_;
				f2[counter] = this->f2_new_;
				fCass[counter] = this->fCass_new_;
				s[counter] = this->s_new_;
				r[counter] = this->r_new_;
				R_prime[counter] = this->R_prime_new_;
				Ca_i[counter] = this->Ca_i_new_;
				Ca_SR[counter] = this->Ca_SR_new_;
				Ca_ss[counter] = this->Ca_ss_new_;
				Na_i[counter] = this->Na_i_new_;
				K_i[counter] = this->K_i_new_;
				counter++;
			}
			this->V_old_ = this->V_new_;
			this->Xr1_old_ = this->Xr1_new_;
			this->Xr2_old_ = this->Xr2_new_;
			this->Xs_old_ = this->Xs_new_;
			this->m_old_ = this->m_new_;
			this->h_old_ = this->h_new_;
			this->j_old_ = this->j_new_;
			this->d_old_ = this->d_new_;
			this->f_old_ = this->f_new_;
			this->f2_old_ = this->f2_new_;
			this->fCass_old_ = this->fCass_new_;
			this->s_old_ = this->s_new_;
			this->r_old_ = this->r_new_;
			this->R_prime_old_ = this->R_prime_new_;
			this->Ca_i_old_ = this->Ca_i_new_;
			this->Ca_SR_old_ = this->Ca_SR_new_;
			this->Ca_ss_old_ = this->Ca_ss_new_;
			this->Na_i_old_ = this->Na_i_new_;
			this->K_i_old_ = this->K_i_new_;
		}
		double h_jac_[numEDO];
		double quociente = 1000.0;
		h_jac_[0] = fabs(_agos_min(V, num_results__) / _agos_max(V, num_results__) );
		h_jac_[1] = fabs(_agos_min(Xr1, num_results__) / _agos_max(Xr1, num_results__) );
		h_jac_[2] = fabs(_agos_min(Xr2, num_results__) / _agos_max(Xr2, num_results__) );
		h_jac_[3] = fabs(_agos_min(Xs, num_results__) / _agos_max(Xs, num_results__) );
		h_jac_[4] = fabs(_agos_min(m, num_results__) / _agos_max(m, num_results__) );
		h_jac_[5] = fabs(_agos_min(h, num_results__) / _agos_max(h, num_results__) );
		h_jac_[6] = fabs(_agos_min(j, num_results__) / _agos_max(j, num_results__) );
		h_jac_[7] = fabs(_agos_min(d, num_results__) / _agos_max(d, num_results__) );
		h_jac_[8] = fabs(_agos_min(f, num_results__) / _agos_max(f, num_results__) );
		h_jac_[9] = fabs(_agos_min(f2, num_results__) / _agos_max(f2, num_results__) );
		h_jac_[10] = fabs(_agos_min(fCass, num_results__) / _agos_max(fCass, num_results__) );
		h_jac_[11] = fabs(_agos_min(s, num_results__) / _agos_max(s, num_results__) );
		h_jac_[12] = fabs(_agos_min(r, num_results__) / _agos_max(r, num_results__) );
		h_jac_[13] = fabs(_agos_min(R_prime, num_results__) / _agos_max(R_prime, num_results__) );
		h_jac_[14] = fabs(_agos_min(Ca_i, num_results__) / _agos_max(Ca_i, num_results__) );
		h_jac_[15] = fabs(_agos_min(Ca_SR, num_results__) / _agos_max(Ca_SR, num_results__) );
		h_jac_[16] = fabs(_agos_min(Ca_ss, num_results__) / _agos_max(Ca_ss, num_results__) );
		h_jac_[17] = fabs(_agos_min(Na_i, num_results__) / _agos_max(Na_i, num_results__) );
		h_jac_[18] = fabs(_agos_min(K_i, num_results__) / _agos_max(K_i, num_results__) );
		for(int l=0;l<numEDO;l++){
			h_jac_[l] = (h_jac_[l]==0 || h_jac_[l]==AGOS_NAN || h_jac_[l]==AGOS_INF)?this->dtime:h_jac_[l];
		}
		this->timeSaving = this->dtime;

		this->time_new = this->time;

		V_old_ = V_ini_;
		Xr1_old_ = Xr1_ini_;
		Xr2_old_ = Xr2_ini_;
		Xs_old_ = Xs_ini_;
		m_old_ = m_ini_;
		h_old_ = h_ini_;
		j_old_ = j_ini_;
		d_old_ = d_ini_;
		f_old_ = f_ini_;
		f2_old_ = f2_ini_;
		fCass_old_ = fCass_ini_;
		s_old_ = s_ini_;
		r_old_ = r_ini_;
		R_prime_old_ = R_prime_ini_;
		Ca_i_old_ = Ca_i_ini_;
		Ca_SR_old_ = Ca_SR_ini_;
		Ca_ss_old_ = Ca_ss_ini_;
		Na_i_old_ = Na_i_ini_;
		K_i_old_ = K_i_ini_;
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
			this->Xr1_new_ = this->dtime*(this->Xr1_lado_direito_) + this->Xr1_old_;
			this->Xr2_new_ = this->dtime*(this->Xr2_lado_direito_) + this->Xr2_old_;
			this->Xs_new_ = this->dtime*(this->Xs_lado_direito_) + this->Xs_old_;
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h_new_ = this->dtime*(this->h_lado_direito_) + this->h_old_;
			this->j_new_ = this->dtime*(this->j_lado_direito_) + this->j_old_;
			this->d_new_ = this->dtime*(this->d_lado_direito_) + this->d_old_;
			this->f_new_ = this->dtime*(this->f_lado_direito_) + this->f_old_;
			this->f2_new_ = this->dtime*(this->f2_lado_direito_) + this->f2_old_;
			this->fCass_new_ = this->dtime*(this->fCass_lado_direito_) + this->fCass_old_;
			this->s_new_ = this->dtime*(this->s_lado_direito_) + this->s_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->R_prime_new_ = this->dtime*(this->R_prime_lado_direito_) + this->R_prime_old_;
			this->Ca_i_new_ = this->dtime*(this->Ca_i_lado_direito_) + this->Ca_i_old_;
			this->Ca_SR_new_ = this->dtime*(this->Ca_SR_lado_direito_) + this->Ca_SR_old_;
			this->Ca_ss_new_ = this->dtime*(this->Ca_ss_lado_direito_) + this->Ca_ss_old_;
			this->Na_i_new_ = this->dtime*(this->Na_i_lado_direito_) + this->Na_i_old_;
			this->K_i_new_ = this->dtime*(this->K_i_lado_direito_) + this->K_i_old_;
			diff =  _agos_round(this->time_new - timeSaving, 5);
			if(diff==0){
				this->timeSaving += svRate;
				V[counter2] = V_new_;
				Xr1[counter2] = Xr1_new_;
				Xr2[counter2] = Xr2_new_;
				Xs[counter2] = Xs_new_;
				m[counter2] = m_new_;
				h[counter2] = h_new_;
				j[counter2] = j_new_;
				d[counter2] = d_new_;
				f[counter2] = f_new_;
				f2[counter2] = f2_new_;
				fCass[counter2] = fCass_new_;
				s[counter2] = s_new_;
				r[counter2] = r_new_;
				R_prime[counter2] = R_prime_new_;
				Ca_i[counter2] = Ca_i_new_;
				Ca_SR[counter2] = Ca_SR_new_;
				Ca_ss[counter2] = Ca_ss_new_;
				Na_i[counter2] = Na_i_new_;
				K_i[counter2] = K_i_new_;
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
			Xr1_old_ = Xr1_new_;
			Xr2_old_ = Xr2_new_;
			Xs_old_ = Xs_new_;
			m_old_ = m_new_;
			h_old_ = h_new_;
			j_old_ = j_new_;
			d_old_ = d_new_;
			f_old_ = f_new_;
			f2_old_ = f2_new_;
			fCass_old_ = fCass_new_;
			s_old_ = s_new_;
			r_old_ = r_new_;
			R_prime_old_ = R_prime_new_;
			Ca_i_old_ = Ca_i_new_;
			Ca_SR_old_ = Ca_SR_new_;
			Ca_ss_old_ = Ca_ss_new_;
			Na_i_old_ = Na_i_new_;
			K_i_old_ = K_i_new_;
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
		case 1:		return Xr1;    break;
		case 2:		return Xr2;    break;
		case 3:		return Xs;    break;
		case 4:		return m;    break;
		case 5:		return h;    break;
		case 6:		return j;    break;
		case 7:		return d;    break;
		case 8:		return f;    break;
		case 9:		return f2;    break;
		case 10:		return fCass;    break;
		case 11:		return s;    break;
		case 12:		return r;    break;
		case 13:		return R_prime;    break;
		case 14:		return Ca_i;    break;
		case 15:		return Ca_SR;    break;
		case 16:		return Ca_ss;    break;
		case 17:		return Na_i;    break;
		case 18:		return K_i;    break;
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
	double Solveode::ifnumber_1(){
		if((V_old_<(-4.0000000000e+01))){
			return ((5.7000000000e-02*exp(((-(V_old_+8.0000000000e+01))/6.8000000000e+00))));
		}else{
			return (0.0000000000e+00);
		}
	}
	double Solveode::ifnumber_2(){
		if((V_old_<(-4.0000000000e+01))){
			return (((2.7000000000e+00*exp((7.9000000000e-02*V_old_)))+(3.1000000000e+05*exp((3.4850000000e-01*V_old_)))));
		}else{
			return ((7.7000000000e-01/(1.3000000000e-01*(1.0000000000e+00+exp(((V_old_+1.0660000000e+01)/(-1.1100000000e+01)))))));
		}
	}
	double Solveode::ifnumber_3(){
		if((V_old_<(-4.0000000000e+01))){
			return (((((((-2.5428000000e+04)*exp((2.4440000000e-01*V_old_)))-(6.9480000000e-06*exp(((-4.3910000000e-02)*V_old_))))*(V_old_+3.7780000000e+01))/1.0000000000e+00)/(1.0000000000e+00+exp((3.1100000000e-01*(V_old_+7.9230000000e+01))))));
		}else{
			return (0.0000000000e+00);
		}
	}
	double Solveode::ifnumber_4(){
		if((V_old_<(-4.0000000000e+01))){
			return (((2.4240000000e-02*exp(((-1.0520000000e-02)*V_old_)))/(1.0000000000e+00+exp(((-1.3780000000e-01)*(V_old_+4.0140000000e+01))))));
		}else{
			return (((6.0000000000e-01*exp((5.7000000000e-02*V_old_)))/(1.0000000000e+00+exp(((-1.0000000000e-01)*(V_old_+3.2000000000e+01))))));
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
		double dtL, dtM, dtMax=0.0,  dtMin=999 ;
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
			dependent_variable__ = N_VNew_Serial(19);
			if(check_flag((void *)dependent_variable__, "N_VNew_Serial", 0))
			exit(1);
			depvar__ = (double*)malloc(sizeof(double)*19);
			if(depvar__ == NULL){
			fprintf(stderr, "ERROR Cannot allocate memory for depvar__\n");
			exit(0);
			}
			NV_Ith_S(dependent_variable__, 0) = V_ini_;
			NV_Ith_S(dependent_variable__, 1) = Xr1_ini_;
			NV_Ith_S(dependent_variable__, 2) = Xr2_ini_;
			NV_Ith_S(dependent_variable__, 3) = Xs_ini_;
			NV_Ith_S(dependent_variable__, 4) = m_ini_;
			NV_Ith_S(dependent_variable__, 5) = h_ini_;
			NV_Ith_S(dependent_variable__, 6) = j_ini_;
			NV_Ith_S(dependent_variable__, 7) = d_ini_;
			NV_Ith_S(dependent_variable__, 8) = f_ini_;
			NV_Ith_S(dependent_variable__, 9) = f2_ini_;
			NV_Ith_S(dependent_variable__, 10) = fCass_ini_;
			NV_Ith_S(dependent_variable__, 11) = s_ini_;
			NV_Ith_S(dependent_variable__, 12) = r_ini_;
			NV_Ith_S(dependent_variable__, 13) = R_prime_ini_;
			NV_Ith_S(dependent_variable__, 14) = Ca_i_ini_;
			NV_Ith_S(dependent_variable__, 15) = Ca_SR_ini_;
			NV_Ith_S(dependent_variable__, 16) = Ca_ss_ini_;
			NV_Ith_S(dependent_variable__, 17) = Na_i_ini_;
			NV_Ith_S(dependent_variable__, 18) = K_i_ini_;
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
			flag__ = CVDense(cvode_mem_cvode__, 19);
			if (check_flag(&flag__, "CVDense", 1))	exit(1);
			break;
			case 2:
			flag__ = CVDiag(cvode_mem_cvode__);
			if (check_flag(&flag__, "CVDiag", 1))	exit(1);
			break;
			case 3:
			flag__ = CVBand(cvode_mem_cvode__, 19, NULL, NULL);
			if (check_flag(&flag__, "CVBand", 1))	exit(1);
			break;
			}
			CVodeSetFdata(cvode_mem_cvode__, (void*)this);
			V_old_ = V_ini_;
			Xr1_old_ = Xr1_ini_;
			Xr2_old_ = Xr2_ini_;
			Xs_old_ = Xs_ini_;
			m_old_ = m_ini_;
			h_old_ = h_ini_;
			j_old_ = j_ini_;
			d_old_ = d_ini_;
			f_old_ = f_ini_;
			f2_old_ = f2_ini_;
			fCass_old_ = fCass_ini_;
			s_old_ = s_ini_;
			r_old_ = r_ini_;
			R_prime_old_ = R_prime_ini_;
			Ca_i_old_ = Ca_i_ini_;
			Ca_SR_old_ = Ca_SR_ini_;
			Ca_ss_old_ = Ca_ss_ini_;
			Na_i_old_ = Na_i_ini_;
			K_i_old_ = K_i_ini_;
		}
		while(1){
			flag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &time_new, CV_NORMAL);
			V_old_ = NV_Ith_S(dependent_variable__, 0);
			Xr1_old_ = NV_Ith_S(dependent_variable__, 1);
			Xr2_old_ = NV_Ith_S(dependent_variable__, 2);
			Xs_old_ = NV_Ith_S(dependent_variable__, 3);
			m_old_ = NV_Ith_S(dependent_variable__, 4);
			h_old_ = NV_Ith_S(dependent_variable__, 5);
			j_old_ = NV_Ith_S(dependent_variable__, 6);
			d_old_ = NV_Ith_S(dependent_variable__, 7);
			f_old_ = NV_Ith_S(dependent_variable__, 8);
			f2_old_ = NV_Ith_S(dependent_variable__, 9);
			fCass_old_ = NV_Ith_S(dependent_variable__, 10);
			s_old_ = NV_Ith_S(dependent_variable__, 11);
			r_old_ = NV_Ith_S(dependent_variable__, 12);
			R_prime_old_ = NV_Ith_S(dependent_variable__, 13);
			Ca_i_old_ = NV_Ith_S(dependent_variable__, 14);
			Ca_SR_old_ = NV_Ith_S(dependent_variable__, 15);
			Ca_ss_old_ = NV_Ith_S(dependent_variable__, 16);
			Na_i_old_ = NV_Ith_S(dependent_variable__, 17);
			K_i_old_ = NV_Ith_S(dependent_variable__, 18);
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
		long int nsteps, fails;
		int flag = CVodeGetNumSteps(cvode_mem_cvode__ , &nsteps);
		flag = CVodeGetNumErrTestFails(cvode_mem_cvode__, &fails);
		
		
		printf("max: %e min: %e medias: %e %e %d fails %d\n", dtMax, dtMin, dtM/nsteps, dtM/iout , nsteps, fails);
	}
static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__){
	Solveode *ode = (Solveode *) f_data__;
	for(int i = 0; i<19; i++)
		ode->setVariables( i ,NV_Ith_S(dependent_variable__, i));
	ode->setParameters(0,time);
	rightHandSideFunction.function(ode);
	for(int i = 0; i<19; i++)
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
