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
#define forestSize 15
#include <sys/time.h>
#include "greedy.cpp"
#define _INCREASE_DT_ 1.5
#define _DECREASE_DT_ 0.65
#define _DECREASE_DT_2_ 2
#define numEDO 22
#define numAux 60

class Solveode
{
public:
	int it_countx;
	int *tree_thread;
	double time; 	 // second
	double stim_amplitude; 	 // nanoA
	double stim_start; 	 // second
	double stim_end; 	 // second
	double stim_period; 	 // second
	double stim_duration; 	 // second
	double Cm; 	 // microF
	double R; 	 // joule_per_kilomole_kelvin
	double T; 	 // kelvin
	double F; 	 // coulomb_per_mole
	double Na_o; 	 // millimolar
	double K_o; 	 // millimolar
	double P_kna; 	 // dimensionless
	double Ca_o; 	 // millimolar
	double g_K1; 	 // microS
	double K_mk1; 	 // millimolar
	double g_Kr1; 	 // microS
	double g_Kr2; 	 // microS
	double g_Ks; 	 // microS
	double g_Na; 	 // microS
	double delta_m; 	 // millivolt
	double shift_h; 	 // millivolt
	double g_pna; 	 // microS
	double g_bna; 	 // microS
	double FrICa; 	 // dimensionless
	double P_Ca_L; 	 // nanoA_per_millimolar
	double P_CaK; 	 // dimensionless
	double P_CaNa; 	 // dimensionless
	double speed_d; 	 // dimensionless
	double delta_f; 	 // millivolt
	double speed_f; 	 // dimensionless
	double Km_f2; 	 // millimolar
	double R_decay; 	 // per_second
	double Km_f2ds; 	 // millimolar
	double g_bca; 	 // microS
	double g_to; 	 // microS
	double g_tos; 	 // dimensionless
	double i_NaK_max; 	 // nanoA
	double K_mK; 	 // millimolar
	double K_mNa; 	 // millimolar
	double FRiNaCa; 	 // dimensionless
	double k_NaCa; 	 // nanoA
	double gamma; 	 // dimensionless
	double n_NaCa; 	 // dimensionless
	double d_NaCa; 	 // dimensionless
	double K_cyca; 	 // millimolar
	double K_xcs; 	 // dimensionless
	double K_srca; 	 // millimolar
	double alpha_up; 	 // millimolar_per_second
	double beta_up; 	 // millimolar_per_second
	double K_m_Ca_cyt; 	 // millimolar
	double K_m_Ca_ds; 	 // millimolar
	double K_m_rel; 	 // per_second
	double K_leak_rate; 	 // per_second
	double radius; 	 // micrometre
	double length; 	 // micrometre
	double V_e_ratio; 	 // dimensionless
	double V_up_ratio; 	 // dimensionless
	double V_rel_ratio; 	 // dimensionless
	double V_ds_ratio; 	 // dimensionless
	double Kdecay; 	 // per_second
	double alpha_Calmod; 	 // per_millimolar_second
	double Calmod; 	 // millimolar
	double beta_Calmod; 	 // per_second
	double alpha_Trop; 	 // per_millimolar_second
	double Trop; 	 // millimolar
	double beta_Trop; 	 // per_second
	double calc_i_Stim; 	 // nanoA
	double calc_E_Na; 	 // millivolt
	double calc_E_K; 	 // millivolt
	double calc_E_Ks; 	 // millivolt
	double calc_E_Ca; 	 // millivolt
	double calc_E_mh; 	 // millivolt
	double calc_i_K1; 	 // nanoA
	double calc_i_Kr; 	 // nanoA
	double calc_alpha_xr1; 	 // per_second
	double calc_beta_xr1; 	 // per_second
	double calc_alpha_xr2; 	 // per_second
	double calc_beta_xr2; 	 // per_second
	double calc_i_Ks; 	 // nanoA
	double calc_alpha_xs; 	 // per_second
	double calc_beta_xs; 	 // per_second
	double calc_i_Na; 	 // nanoA
	double calc_E0_m; 	 // millivolt
	double calc_alpha_m; 	 // per_second
	double calc_beta_m; 	 // per_second
	double calc_alpha_h; 	 // per_second
	double calc_beta_h; 	 // per_second
	double calc_i_p_Na; 	 // nanoA
	double calc_i_b_Na; 	 // nanoA
	double calc_i_Ca_L_Ca_cyt; 	 // nanoA
	double calc_i_Ca_L_K_cyt; 	 // nanoA
	double calc_i_Ca_L_Na_cyt; 	 // nanoA
	double calc_i_Ca_L_Ca_ds; 	 // nanoA
	double calc_i_Ca_L_K_ds; 	 // nanoA
	double calc_i_Ca_L_Na_ds; 	 // nanoA
	double calc_i_Ca_L; 	 // nanoA
	double calc_E0_d; 	 // millivolt
	double calc_alpha_d; 	 // per_second
	double calc_beta_d; 	 // per_second
	double calc_E0_f; 	 // millivolt
	double calc_alpha_f; 	 // per_second
	double calc_beta_f; 	 // per_second
	double calc_i_b_Ca; 	 // nanoA
	double calc_i_to; 	 // nanoA
	double calc_alpha_s; 	 // per_second
	double calc_beta_s; 	 // per_second
	double calc_i_NaK; 	 // nanoA
	double calc_i_NaCa_cyt; 	 // nanoA
	double calc_i_NaCa_ds; 	 // nanoA
	double calc_i_NaCa; 	 // nanoA
	double calc_K_1; 	 // dimensionless
	double calc_K_2; 	 // millimolar
	double calc_i_up; 	 // millimolar_per_second
	double calc_i_trans; 	 // millimolar_per_second
	double calc_VoltDep; 	 // dimensionless
	double calc_CaiReg; 	 // dimensionless
	double calc_CadsReg; 	 // dimensionless
	double calc_RegBindSite; 	 // dimensionless
	double calc_ActRate; 	 // per_second
	double calc_InactRate; 	 // per_second
	double calc_SpeedRel; 	 // dimensionless
	double calc_PrecFrac; 	 // dimensionless
	double calc_i_rel; 	 // millimolar_per_second
	double calc_V_Cell; 	 // micrometre3
	double calc_V_i_ratio; 	 // dimensionless
	double calc_V_i; 	 // micrometre3
	double dtime, *time_vec__;
	double time_new;

	//functions variables
	double *V;
	double V_new_, V_old_, V_ini_, V_lado_direito_;
	double *xr1;
	double xr1_new_, xr1_old_, xr1_ini_, xr1_lado_direito_;
	double *xr2;
	double xr2_new_, xr2_old_, xr2_ini_, xr2_lado_direito_;
	double *xs;
	double xs_new_, xs_old_, xs_ini_, xs_lado_direito_;
	double *m;
	double m_new_, m_old_, m_ini_, m_lado_direito_;
	double *h;
	double h_new_, h_old_, h_ini_, h_lado_direito_;
	double *d;
	double d_new_, d_old_, d_ini_, d_lado_direito_;
	double *f;
	double f_new_, f_old_, f_ini_, f_lado_direito_;
	double *f2;
	double f2_new_, f2_old_, f2_ini_, f2_lado_direito_;
	double *f2ds;
	double f2ds_new_, f2ds_old_, f2ds_ini_, f2ds_lado_direito_;
	double *s;
	double s_new_, s_old_, s_ini_, s_lado_direito_;
	double *r;
	double r_new_, r_old_, r_ini_, r_lado_direito_;
	double *ActFrac;
	double ActFrac_new_, ActFrac_old_, ActFrac_ini_, ActFrac_lado_direito_;
	double *ProdFrac;
	double ProdFrac_new_, ProdFrac_old_, ProdFrac_ini_, ProdFrac_lado_direito_;
	double *Na_i;
	double Na_i_new_, Na_i_old_, Na_i_ini_, Na_i_lado_direito_;
	double *K_i;
	double K_i_new_, K_i_old_, K_i_ini_, K_i_lado_direito_;
	double *Ca_i;
	double Ca_i_new_, Ca_i_old_, Ca_i_ini_, Ca_i_lado_direito_;
	double *Ca_ds;
	double Ca_ds_new_, Ca_ds_old_, Ca_ds_ini_, Ca_ds_lado_direito_;
	double *Ca_up;
	double Ca_up_new_, Ca_up_old_, Ca_up_ini_, Ca_up_lado_direito_;
	double *Ca_rel;
	double Ca_rel_new_, Ca_rel_old_, Ca_rel_ini_, Ca_rel_lado_direito_;
	double *Ca_Calmod;
	double Ca_Calmod_new_, Ca_Calmod_old_, Ca_Calmod_ini_, Ca_Calmod_lado_direito_;
	double *Ca_Trop;
	double Ca_Trop_new_, Ca_Trop_old_, Ca_Trop_ini_, Ca_Trop_lado_direito_;

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
	inline double ifnumber_5();
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
	__AGOS->calc_E_mh = (((__AGOS->R*__AGOS->T)/__AGOS->F)*log(((__AGOS->Na_o+(1.2000000000e-01*__AGOS->K_o))/(__AGOS->Na_i_old_+(1.2000000000e-01*__AGOS->K_i_old_)))));
	__AGOS->calc_i_Ca_L_Ca_cyt = (((((1.0000000000e+00-__AGOS->FrICa)*4.0000000000e+00*__AGOS->P_Ca_L*__AGOS->d_old_*__AGOS->f_old_*__AGOS->f2_old_*(__AGOS->V_old_-5.0000000000e+01)*__AGOS->F)/(__AGOS->R*__AGOS->T))/(1.0000000000e+00-exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F*2.0000000000e+00)/(__AGOS->R*__AGOS->T)))))*((__AGOS->Ca_i_old_*exp(((1.0000000000e+02*__AGOS->F)/(__AGOS->R*__AGOS->T))))-(__AGOS->Ca_o*exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F*2.0000000000e+00)/(__AGOS->R*__AGOS->T))))));
	__AGOS->calc_i_Ca_L_K_cyt = (((((1.0000000000e+00-__AGOS->FrICa)*__AGOS->P_CaK*__AGOS->P_Ca_L*__AGOS->d_old_*__AGOS->f_old_*__AGOS->f2_old_*(__AGOS->V_old_-5.0000000000e+01)*__AGOS->F)/(__AGOS->R*__AGOS->T))/(1.0000000000e+00-exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F)/(__AGOS->R*__AGOS->T)))))*((__AGOS->K_i_old_*exp(((5.0000000000e+01*__AGOS->F)/(__AGOS->R*__AGOS->T))))-(__AGOS->K_o*exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F)/(__AGOS->R*__AGOS->T))))));
	__AGOS->calc_i_Ca_L_Na_cyt = (((((1.0000000000e+00-__AGOS->FrICa)*__AGOS->P_CaNa*__AGOS->P_Ca_L*__AGOS->d_old_*__AGOS->f_old_*__AGOS->f2_old_*(__AGOS->V_old_-5.0000000000e+01)*__AGOS->F)/(__AGOS->R*__AGOS->T))/(1.0000000000e+00-exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F)/(__AGOS->R*__AGOS->T)))))*((__AGOS->Na_i_old_*exp(((5.0000000000e+01*__AGOS->F)/(__AGOS->R*__AGOS->T))))-(__AGOS->Na_o*exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F)/(__AGOS->R*__AGOS->T))))));
	__AGOS->calc_i_Ca_L_Ca_ds = ((((__AGOS->FrICa*4.0000000000e+00*__AGOS->P_Ca_L*__AGOS->d_old_*__AGOS->f_old_*__AGOS->f2ds_old_*(__AGOS->V_old_-5.0000000000e+01)*__AGOS->F)/(__AGOS->R*__AGOS->T))/(1.0000000000e+00-exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F*2.0000000000e+00)/(__AGOS->R*__AGOS->T)))))*((__AGOS->Ca_i_old_*exp(((1.0000000000e+02*__AGOS->F)/(__AGOS->R*__AGOS->T))))-(__AGOS->Ca_o*exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F*2.0000000000e+00)/(__AGOS->R*__AGOS->T))))));
	__AGOS->calc_i_Ca_L_K_ds = ((((__AGOS->FrICa*__AGOS->P_CaK*__AGOS->P_Ca_L*__AGOS->d_old_*__AGOS->f_old_*__AGOS->f2ds_old_*(__AGOS->V_old_-5.0000000000e+01)*__AGOS->F)/(__AGOS->R*__AGOS->T))/(1.0000000000e+00-exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F)/(__AGOS->R*__AGOS->T)))))*((__AGOS->K_i_old_*exp(((5.0000000000e+01*__AGOS->F)/(__AGOS->R*__AGOS->T))))-(__AGOS->K_o*exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F)/(__AGOS->R*__AGOS->T))))));
	__AGOS->calc_i_Ca_L_Na_ds = ((((__AGOS->FrICa*__AGOS->P_CaNa*__AGOS->P_Ca_L*__AGOS->d_old_*__AGOS->f_old_*__AGOS->f2ds_old_*(__AGOS->V_old_-5.0000000000e+01)*__AGOS->F)/(__AGOS->R*__AGOS->T))/(1.0000000000e+00-exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F)/(__AGOS->R*__AGOS->T)))))*((__AGOS->Na_i_old_*exp(((5.0000000000e+01*__AGOS->F)/(__AGOS->R*__AGOS->T))))-(__AGOS->Na_o*exp((((-(__AGOS->V_old_-5.0000000000e+01))*__AGOS->F)/(__AGOS->R*__AGOS->T))))));
	__AGOS->calc_i_NaK = ((((__AGOS->i_NaK_max*__AGOS->K_o)/(__AGOS->K_mK+__AGOS->K_o))*__AGOS->Na_i_old_)/(__AGOS->K_mNa+__AGOS->Na_i_old_));
	__AGOS->calc_i_NaCa_cyt = (((1.0000000000e+00-__AGOS->FRiNaCa)*__AGOS->k_NaCa*((exp(((__AGOS->gamma*(__AGOS->n_NaCa-2.0000000000e+00)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))*pow(__AGOS->Na_i_old_,__AGOS->n_NaCa)*__AGOS->Ca_o)-(exp((((__AGOS->gamma-1.0000000000e+00)*(__AGOS->n_NaCa-2.0000000000e+00)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))*pow(__AGOS->Na_o,__AGOS->n_NaCa)*__AGOS->Ca_i_old_)))/((1.0000000000e+00+(__AGOS->d_NaCa*((__AGOS->Ca_i_old_*pow(__AGOS->Na_o,__AGOS->n_NaCa))+(__AGOS->Ca_o*pow(__AGOS->Na_i_old_,__AGOS->n_NaCa)))))*(1.0000000000e+00+(__AGOS->Ca_i_old_/6.9000000000e-03))));
	__AGOS->calc_i_NaCa_ds = ((__AGOS->FRiNaCa*__AGOS->k_NaCa*((exp(((__AGOS->gamma*(__AGOS->n_NaCa-2.0000000000e+00)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))*pow(__AGOS->Na_i_old_,__AGOS->n_NaCa)*__AGOS->Ca_o)-(exp((((__AGOS->gamma-1.0000000000e+00)*(__AGOS->n_NaCa-2.0000000000e+00)*__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))*pow(__AGOS->Na_o,__AGOS->n_NaCa)*__AGOS->Ca_ds_old_)))/((1.0000000000e+00+(__AGOS->d_NaCa*((__AGOS->Ca_ds_old_*pow(__AGOS->Na_o,__AGOS->n_NaCa))+(__AGOS->Ca_o*pow(__AGOS->Na_i_old_,__AGOS->n_NaCa)))))*(1.0000000000e+00+(__AGOS->Ca_ds_old_/6.9000000000e-03))));
	__AGOS->calc_K_1 = ((__AGOS->K_cyca*__AGOS->K_xcs)/__AGOS->K_srca);
	__AGOS->calc_i_trans = (5.0000000000e+01*(__AGOS->Ca_up_old_-__AGOS->Ca_rel_old_));
	__AGOS->calc_i_rel = (((pow((__AGOS->ActFrac_old_/(__AGOS->ActFrac_old_+2.5000000000e-01)),2.0000000000e+00)*__AGOS->K_m_rel)+__AGOS->K_leak_rate)*__AGOS->Ca_rel_old_);
	__AGOS->calc_V_Cell = (3.1415926540e+00*pow(__AGOS->radius,2.0000000000e+00)*__AGOS->length);
	__AGOS->calc_V_i_ratio = (((1.0000000000e+00-__AGOS->V_e_ratio)-__AGOS->V_up_ratio)-__AGOS->V_rel_ratio);
	__AGOS->calc_i_K1 = ((((__AGOS->g_K1*__AGOS->K_o)/(__AGOS->K_o+__AGOS->K_mk1))*(__AGOS->V_old_-__AGOS->calc_E_K))/(1.0000000000e+00+exp(((((__AGOS->V_old_-__AGOS->calc_E_K)-1.0000000000e+01)*__AGOS->F*1.2500000000e+00)/(__AGOS->R*__AGOS->T)))));
	__AGOS->calc_i_Kr = (((((__AGOS->g_Kr1*__AGOS->xr1_old_)+(__AGOS->g_Kr2*__AGOS->xr2_old_))*1.0000000000e+00)/(1.0000000000e+00+exp(((__AGOS->V_old_+9.0000000000e+00)/2.2400000000e+01))))*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_Ks = (__AGOS->g_Ks*pow(__AGOS->xs_old_,2.0000000000e+00)*(__AGOS->V_old_-__AGOS->calc_E_Ks));
	__AGOS->calc_i_Na = (__AGOS->g_Na*pow(__AGOS->m_old_,3.0000000000e+00)*__AGOS->h_old_*(__AGOS->V_old_-__AGOS->calc_E_mh));
	__AGOS->calc_i_p_Na = (((__AGOS->g_pna*1.0000000000e+00)/(1.0000000000e+00+exp(((-(__AGOS->V_old_+5.2000000000e+01))/8.0000000000e+00))))*(__AGOS->V_old_-__AGOS->calc_E_Na));
	__AGOS->calc_i_b_Na = (__AGOS->g_bna*(__AGOS->V_old_-__AGOS->calc_E_Na));
	__AGOS->calc_i_Ca_L = (__AGOS->calc_i_Ca_L_Ca_cyt+__AGOS->calc_i_Ca_L_K_cyt+__AGOS->calc_i_Ca_L_Na_cyt+__AGOS->calc_i_Ca_L_Ca_ds+__AGOS->calc_i_Ca_L_K_ds+__AGOS->calc_i_Ca_L_Na_ds);
	__AGOS->calc_i_b_Ca = (__AGOS->g_bca*(__AGOS->V_old_-__AGOS->calc_E_Ca));
	__AGOS->calc_i_to = (__AGOS->g_to*(__AGOS->g_tos+(__AGOS->s_old_*(1.0000000000e+00-__AGOS->g_tos)))*__AGOS->r_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_NaCa = (__AGOS->calc_i_NaCa_cyt+__AGOS->calc_i_NaCa_ds);
	__AGOS->calc_K_2 = (__AGOS->Ca_i_old_+(__AGOS->Ca_up_old_*__AGOS->calc_K_1)+(__AGOS->K_cyca*__AGOS->K_xcs)+__AGOS->K_cyca);
	__AGOS->calc_V_i = (__AGOS->calc_V_Cell*__AGOS->calc_V_i_ratio);
	__AGOS->calc_i_up = (((__AGOS->Ca_i_old_/__AGOS->calc_K_2)*__AGOS->alpha_up)-(((__AGOS->Ca_up_old_*__AGOS->calc_K_1)/__AGOS->calc_K_2)*__AGOS->beta_up));
	__AGOS->V_lado_direito_= (((-1.0000000000e+00)/__AGOS->Cm)*(__AGOS->calc_i_Stim+__AGOS->calc_i_K1+__AGOS->calc_i_to+__AGOS->calc_i_Kr+__AGOS->calc_i_Ks+__AGOS->calc_i_NaK+__AGOS->calc_i_Na+__AGOS->calc_i_b_Na+__AGOS->calc_i_p_Na+__AGOS->calc_i_Ca_L_Na_cyt+__AGOS->calc_i_Ca_L_Na_ds+__AGOS->calc_i_NaCa_cyt+__AGOS->calc_i_NaCa_ds+__AGOS->calc_i_Ca_L_Ca_cyt+__AGOS->calc_i_Ca_L_Ca_ds+__AGOS->calc_i_Ca_L_K_cyt+__AGOS->calc_i_Ca_L_K_ds+__AGOS->calc_i_b_Ca));
	__AGOS->Na_i_lado_direito_= (((-1.0000000000e+00)/(1.0000000000e+00*__AGOS->calc_V_i*__AGOS->F))*(__AGOS->calc_i_Na+__AGOS->calc_i_p_Na+__AGOS->calc_i_b_Na+(3.0000000000e+00*__AGOS->calc_i_NaK)+(3.0000000000e+00*__AGOS->calc_i_NaCa_cyt)+__AGOS->calc_i_Ca_L_Na_cyt+__AGOS->calc_i_Ca_L_Na_ds));
	__AGOS->K_i_lado_direito_= (((-1.0000000000e+00)/(1.0000000000e+00*__AGOS->calc_V_i*__AGOS->F))*((__AGOS->calc_i_K1+__AGOS->calc_i_Kr+__AGOS->calc_i_Ks+__AGOS->calc_i_Ca_L_K_cyt+__AGOS->calc_i_Ca_L_K_ds+__AGOS->calc_i_to)-(2.0000000000e+00*__AGOS->calc_i_NaK)));
	__AGOS->Ca_i_lado_direito_= (((((((-1.0000000000e+00)/(2.0000000000e+00*1.0000000000e+00*__AGOS->calc_V_i*__AGOS->F))*(((__AGOS->calc_i_Ca_L_Ca_cyt+__AGOS->calc_i_b_Ca)-(2.0000000000e+00*__AGOS->calc_i_NaCa_cyt))-(2.0000000000e+00*__AGOS->calc_i_NaCa_ds)))+(__AGOS->Ca_ds_old_*__AGOS->V_ds_ratio*__AGOS->Kdecay)+((__AGOS->calc_i_rel*__AGOS->V_rel_ratio)/__AGOS->calc_V_i_ratio))-((__AGOS->alpha_Calmod*__AGOS->Ca_i_old_*(__AGOS->Calmod-__AGOS->Ca_Calmod_old_))-(__AGOS->beta_Calmod*__AGOS->Ca_Calmod_old_)))-((__AGOS->alpha_Trop*__AGOS->Ca_i_old_*(__AGOS->Trop-__AGOS->Ca_Trop_old_))-(__AGOS->beta_Trop*__AGOS->Ca_Trop_old_)))-__AGOS->calc_i_up);
	__AGOS->Ca_ds_lado_direito_= ((((-1.0000000000e+00)*__AGOS->calc_i_Ca_L_Ca_ds)/(2.0000000000e+00*1.0000000000e+00*__AGOS->V_ds_ratio*__AGOS->calc_V_i*__AGOS->F))-(__AGOS->Ca_ds_old_*__AGOS->Kdecay));
	__AGOS->Ca_up_lado_direito_= (((__AGOS->calc_V_i_ratio/__AGOS->V_up_ratio)*__AGOS->calc_i_up)-__AGOS->calc_i_trans);
	__AGOS->Ca_rel_lado_direito_= (((__AGOS->V_up_ratio/__AGOS->V_rel_ratio)*__AGOS->calc_i_trans)-__AGOS->calc_i_rel);
} //fim

void __tree2__( Solveode *__AGOS){
	__AGOS->calc_alpha_xr1 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__AGOS->V_old_-5.0000000000e+00))/9.0000000000e+00))));
	__AGOS->calc_beta_xr1 = (5.0000000000e-02*exp(((-(__AGOS->V_old_-2.0000000000e+01))/1.5000000000e+01)));
	__AGOS->xr1_lado_direito_= ((__AGOS->calc_alpha_xr1*(1.0000000000e+00-__AGOS->xr1_old_))-(__AGOS->calc_beta_xr1*__AGOS->xr1_old_));
} //fim

void __tree3__( Solveode *__AGOS){
	__AGOS->calc_alpha_xr2 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__AGOS->V_old_-5.0000000000e+00))/9.0000000000e+00))));
	__AGOS->calc_beta_xr2 = (4.0000000000e-01*exp((-pow(((__AGOS->V_old_+3.0000000000e+01)/3.0000000000e+01),3.0000000000e+00))));
	__AGOS->xr2_lado_direito_= ((__AGOS->calc_alpha_xr2*(1.0000000000e+00-__AGOS->xr2_old_))-(__AGOS->calc_beta_xr2*__AGOS->xr2_old_));
} //fim

void __tree4__( Solveode *__AGOS){
	__AGOS->calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__AGOS->V_old_-4.0000000000e+01))/9.0000000000e+00))));
	__AGOS->calc_beta_xs = (1.0000000000e+00*exp(((-__AGOS->V_old_)/4.5000000000e+01)));
	__AGOS->xs_lado_direito_= ((__AGOS->calc_alpha_xs*(1.0000000000e+00-__AGOS->xs_old_))-(__AGOS->calc_beta_xs*__AGOS->xs_old_));
} //fim

void __tree5__( Solveode *__AGOS){
	__AGOS->calc_E0_m = (__AGOS->V_old_+4.1000000000e+01);
	__AGOS->calc_beta_m = (8.0000000000e+03*exp(((-5.6000000000e-02)*(__AGOS->V_old_+6.6000000000e+01))));
	__AGOS->calc_alpha_m = ((fabs(__AGOS->calc_E0_m)<__AGOS->delta_m))
?(2.0000000000e+03)
:(((2.0000000000e+02*__AGOS->calc_E0_m)/(1.0000000000e+00-exp(((-1.0000000000e-01)*__AGOS->calc_E0_m)))));
	__AGOS->m_lado_direito_= ((__AGOS->calc_alpha_m*(1.0000000000e+00-__AGOS->m_old_))-(__AGOS->calc_beta_m*__AGOS->m_old_));
} //fim

void __tree6__( Solveode *__AGOS){
	__AGOS->calc_alpha_h = (2.0000000000e+01*exp(((-1.2500000000e-01)*((__AGOS->V_old_+7.5000000000e+01)-__AGOS->shift_h))));
	__AGOS->calc_beta_h = (2.0000000000e+03/(1.0000000000e+00+(3.2000000000e+02*exp(((-1.0000000000e-01)*((__AGOS->V_old_+7.5000000000e+01)-__AGOS->shift_h))))));
	__AGOS->h_lado_direito_= ((__AGOS->calc_alpha_h*(1.0000000000e+00-__AGOS->h_old_))-(__AGOS->calc_beta_h*__AGOS->h_old_));
} //fim

void __tree7__( Solveode *__AGOS){
	__AGOS->calc_E0_d = ((__AGOS->V_old_+2.4000000000e+01)-5.0000000000e+00);
	__AGOS->calc_alpha_d = ((fabs(__AGOS->calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((3.0000000000e+01*__AGOS->calc_E0_d)/(1.0000000000e+00-exp(((-__AGOS->calc_E0_d)/4.0000000000e+00)))));
	__AGOS->calc_beta_d = ((fabs(__AGOS->calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((1.2000000000e+01*__AGOS->calc_E0_d)/(exp((__AGOS->calc_E0_d/1.0000000000e+01))-1.0000000000e+00)));
	__AGOS->d_lado_direito_= (__AGOS->speed_d*((__AGOS->calc_alpha_d*(1.0000000000e+00-__AGOS->d_old_))-(__AGOS->calc_beta_d*__AGOS->d_old_)));
} //fim

void __tree8__( Solveode *__AGOS){
	__AGOS->calc_E0_f = (__AGOS->V_old_+3.4000000000e+01);
	__AGOS->calc_beta_f = (1.2000000000e+01/(1.0000000000e+00+exp((((-1.0000000000e+00)*(__AGOS->V_old_+3.4000000000e+01))/4.0000000000e+00))));
	__AGOS->calc_alpha_f = ((fabs(__AGOS->calc_E0_f)<__AGOS->delta_f))
?(2.5000000000e+01)
:(((6.2500000000e+00*__AGOS->calc_E0_f)/(exp((__AGOS->calc_E0_f/4.0000000000e+00))-1.0000000000e+00)));
	__AGOS->f_lado_direito_= (__AGOS->speed_f*((__AGOS->calc_alpha_f*(1.0000000000e+00-__AGOS->f_old_))-(__AGOS->calc_beta_f*__AGOS->f_old_)));
} //fim

void __tree9__( Solveode *__AGOS){
	__AGOS->calc_alpha_s = (3.3000000000e-02*exp(((-__AGOS->V_old_)/1.7000000000e+01)));
	__AGOS->calc_beta_s = (3.3000000000e+01/(1.0000000000e+00+exp(((-1.2500000000e-01)*(__AGOS->V_old_+1.0000000000e+01)))));
	__AGOS->s_lado_direito_= ((__AGOS->calc_alpha_s*(1.0000000000e+00-__AGOS->s_old_))-(__AGOS->calc_beta_s*__AGOS->s_old_));
} //fim

void __tree10__( Solveode *__AGOS){
	__AGOS->calc_VoltDep = exp((8.0000000000e-02*(__AGOS->V_old_-4.0000000000e+01)));
	__AGOS->calc_CaiReg = (__AGOS->Ca_i_old_/(__AGOS->Ca_i_old_+__AGOS->K_m_Ca_cyt));
	__AGOS->calc_CadsReg = (__AGOS->Ca_ds_old_/(__AGOS->Ca_ds_old_+__AGOS->K_m_Ca_ds));
	__AGOS->calc_SpeedRel = ((__AGOS->V_old_<(-5.0000000000e+01)))
?(5.0000000000e+00)
:(1.0000000000e+00);
	__AGOS->calc_PrecFrac = ((1.0000000000e+00-__AGOS->ActFrac_old_)-__AGOS->ProdFrac_old_);
	__AGOS->calc_RegBindSite = (__AGOS->calc_CaiReg+((1.0000000000e+00-__AGOS->calc_CaiReg)*__AGOS->calc_CadsReg));
	__AGOS->calc_ActRate = ((0.0000000000e+00*__AGOS->calc_VoltDep)+(5.0000000000e+02*pow(__AGOS->calc_RegBindSite,2.0000000000e+00)));
	__AGOS->calc_InactRate = (6.0000000000e+01+(5.0000000000e+02*pow(__AGOS->calc_RegBindSite,2.0000000000e+00)));
	__AGOS->ActFrac_lado_direito_= ((__AGOS->calc_PrecFrac*__AGOS->calc_SpeedRel*__AGOS->calc_ActRate)-(__AGOS->ActFrac_old_*__AGOS->calc_SpeedRel*__AGOS->calc_InactRate));
	__AGOS->ProdFrac_lado_direito_= ((__AGOS->ActFrac_old_*__AGOS->calc_SpeedRel*__AGOS->calc_InactRate)-(__AGOS->calc_SpeedRel*1.0000000000e+00*__AGOS->ProdFrac_old_));
} //fim

void __tree11__( Solveode *__AGOS){
	__AGOS->f2_lado_direito_= (1.0000000000e+00-(1.0000000000e+00*((__AGOS->Ca_i_old_/(__AGOS->Km_f2+__AGOS->Ca_i_old_))+__AGOS->f2_old_)));
} //fim

void __tree12__( Solveode *__AGOS){
	__AGOS->f2ds_lado_direito_= (__AGOS->R_decay*(1.0000000000e+00-((__AGOS->Ca_ds_old_/(__AGOS->Km_f2ds+__AGOS->Ca_ds_old_))+__AGOS->f2ds_old_)));
} //fim

void __tree13__( Solveode *__AGOS){
	__AGOS->r_lado_direito_= (3.3300000000e+02*((1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+4.0000000000e+00))/5.0000000000e+00))))-__AGOS->r_old_));
} //fim

void __tree14__( Solveode *__AGOS){
	__AGOS->Ca_Calmod_lado_direito_= ((__AGOS->alpha_Calmod*__AGOS->Ca_i_old_*(__AGOS->Calmod-__AGOS->Ca_Calmod_old_))-(__AGOS->beta_Calmod*__AGOS->Ca_Calmod_old_));
} //fim

void __tree15__( Solveode *__AGOS){
	__AGOS->Ca_Trop_lado_direito_= ((__AGOS->alpha_Trop*__AGOS->Ca_i_old_*(__AGOS->Trop-__AGOS->Ca_Trop_old_))-(__AGOS->beta_Trop*__AGOS->Ca_Trop_old_));
} //fim

void __AGOS_EQUATIONS__( Solveode *__AGOS){
		const double time_new = __AGOS->time_new;
		const double time = 0.0000000000e+00;
		const double stim_amplitude = -6.0000000000e+00;
		const double stim_start = 1.0000000000e-01;
		const double stim_end = 1.0000000000e+05;
		const double stim_period = 1.0000000000e+00;
		const double stim_duration = 1.5000000000e-03;
		const double Cm = 9.5000000000e-05;
		const double R = 8.3144720000e+03;
		const double T = 3.1000000000e+02;
		const double F = 9.6485341500e+04;
		const double Na_o = 1.4000000000e+02;
		const double K_o = 4.0000000000e+00;
		const double P_kna = 3.0000000000e-02;
		const double Ca_o = 2.0000000000e+00;
		const double g_K1 = 5.0000000000e-01;
		const double K_mk1 = 1.0000000000e+01;
		const double g_Kr1 = 2.1000000000e-03;
		const double g_Kr2 = 1.3000000000e-03;
		const double g_Ks = 2.6000000000e-03;
		const double g_Na = 2.5000000000e+00;
		const double delta_m = 1.0000000000e-05;
		const double shift_h = 0.0000000000e+00;
		const double g_pna = 4.0000000000e-03;
		const double g_bna = 6.0000000000e-04;
		const double FrICa = 1.0000000000e+00;
		const double P_Ca_L = 1.0000000000e-01;
		const double P_CaK = 2.0000000000e-03;
		const double P_CaNa = 1.0000000000e-02;
		const double speed_d = 3.0000000000e+00;
		const double delta_f = 1.0000000000e-04;
		const double speed_f = 3.0000000000e-01;
		const double Km_f2 = 1.0000000000e+05;
		const double R_decay = 2.0000000000e+01;
		const double Km_f2ds = 1.0000000000e-03;
		const double g_bca = 2.5000000000e-04;
		const double g_to = 5.0000000000e-03;
		const double g_tos = 0.0000000000e+00;
		const double i_NaK_max = 7.0000000000e-01;
		const double K_mK = 1.0000000000e+00;
		const double K_mNa = 4.0000000000e+01;
		const double FRiNaCa = 1.0000000000e-03;
		const double k_NaCa = 5.0000000000e-04;
		const double gamma = 5.0000000000e-01;
		const double n_NaCa = 3.0000000000e+00;
		const double d_NaCa = 0.0000000000e+00;
		const double K_cyca = 3.0000000000e-04;
		const double K_xcs = 4.0000000000e-01;
		const double K_srca = 5.0000000000e-01;
		const double alpha_up = 4.0000000000e-01;
		const double beta_up = 3.0000000000e-02;
		const double K_m_Ca_cyt = 5.0000000000e-04;
		const double K_m_Ca_ds = 1.0000000000e-02;
		const double K_m_rel = 2.5000000000e+02;
		const double K_leak_rate = 5.0000000000e-02;
		const double radius = 1.2000000000e-02;
		const double length = 7.4000000000e-02;
		const double V_e_ratio = 4.0000000000e-01;
		const double V_up_ratio = 1.0000000000e-02;
		const double V_rel_ratio = 1.0000000000e-01;
		const double V_ds_ratio = 1.0000000000e-01;
		const double Kdecay = 1.0000000000e+01;
		const double alpha_Calmod = 1.0000000000e+05;
		const double Calmod = 2.0000000000e-02;
		const double beta_Calmod = 5.0000000000e+01;
		const double alpha_Trop = 1.0000000000e+05;
		const double Trop = 5.0000000000e-02;
		const double beta_Trop = 2.0000000000e+02;
		const double V_old_= __AGOS->V_old_;
		const double xr1_old_= __AGOS->xr1_old_;
		const double xr2_old_= __AGOS->xr2_old_;
		const double xs_old_= __AGOS->xs_old_;
		const double m_old_= __AGOS->m_old_;
		const double h_old_= __AGOS->h_old_;
		const double d_old_= __AGOS->d_old_;
		const double f_old_= __AGOS->f_old_;
		const double f2_old_= __AGOS->f2_old_;
		const double f2ds_old_= __AGOS->f2ds_old_;
		const double s_old_= __AGOS->s_old_;
		const double r_old_= __AGOS->r_old_;
		const double ActFrac_old_= __AGOS->ActFrac_old_;
		const double ProdFrac_old_= __AGOS->ProdFrac_old_;
		const double Na_i_old_= __AGOS->Na_i_old_;
		const double K_i_old_= __AGOS->K_i_old_;
		const double Ca_i_old_= __AGOS->Ca_i_old_;
		const double Ca_ds_old_= __AGOS->Ca_ds_old_;
		const double Ca_up_old_= __AGOS->Ca_up_old_;
		const double Ca_rel_old_= __AGOS->Ca_rel_old_;
		const double Ca_Calmod_old_= __AGOS->Ca_Calmod_old_;
		const double Ca_Trop_old_= __AGOS->Ca_Trop_old_;
	const double calc_i_Stim = (((time_new>=stim_start)&&(time_new<=stim_end)&&(((time_new-stim_start)-(floor(((time_new-stim_start)/stim_period))*stim_period))<=stim_duration)))
?(stim_amplitude)
:(0.0000000000e+00);
	const double calc_E_Na = (((R*T)/F)*log((Na_o/Na_i_old_)));
	const double calc_E_K = (((R*T)/F)*log((K_o/K_i_old_)));
	const double calc_E_Ks = (((R*T)/F)*log(((K_o+(P_kna*Na_o))/(K_i_old_+(P_kna*Na_i_old_)))));
	const double calc_E_Ca = (((5.0000000000e-01*R*T)/F)*log((Ca_o/Ca_i_old_)));
	const double calc_E_mh = (((R*T)/F)*log(((Na_o+(1.2000000000e-01*K_o))/(Na_i_old_+(1.2000000000e-01*K_i_old_)))));
	const double calc_i_Ca_L_Ca_cyt = (((((1.0000000000e+00-FrICa)*4.0000000000e+00*P_Ca_L*d_old_*f_old_*f2_old_*(V_old_-5.0000000000e+01)*F)/(R*T))/(1.0000000000e+00-exp((((-(V_old_-5.0000000000e+01))*F*2.0000000000e+00)/(R*T)))))*((Ca_i_old_*exp(((1.0000000000e+02*F)/(R*T))))-(Ca_o*exp((((-(V_old_-5.0000000000e+01))*F*2.0000000000e+00)/(R*T))))));
	const double calc_i_Ca_L_K_cyt = (((((1.0000000000e+00-FrICa)*P_CaK*P_Ca_L*d_old_*f_old_*f2_old_*(V_old_-5.0000000000e+01)*F)/(R*T))/(1.0000000000e+00-exp((((-(V_old_-5.0000000000e+01))*F)/(R*T)))))*((K_i_old_*exp(((5.0000000000e+01*F)/(R*T))))-(K_o*exp((((-(V_old_-5.0000000000e+01))*F)/(R*T))))));
	const double calc_i_Ca_L_Na_cyt = (((((1.0000000000e+00-FrICa)*P_CaNa*P_Ca_L*d_old_*f_old_*f2_old_*(V_old_-5.0000000000e+01)*F)/(R*T))/(1.0000000000e+00-exp((((-(V_old_-5.0000000000e+01))*F)/(R*T)))))*((Na_i_old_*exp(((5.0000000000e+01*F)/(R*T))))-(Na_o*exp((((-(V_old_-5.0000000000e+01))*F)/(R*T))))));
	const double calc_i_Ca_L_Ca_ds = ((((FrICa*4.0000000000e+00*P_Ca_L*d_old_*f_old_*f2ds_old_*(V_old_-5.0000000000e+01)*F)/(R*T))/(1.0000000000e+00-exp((((-(V_old_-5.0000000000e+01))*F*2.0000000000e+00)/(R*T)))))*((Ca_i_old_*exp(((1.0000000000e+02*F)/(R*T))))-(Ca_o*exp((((-(V_old_-5.0000000000e+01))*F*2.0000000000e+00)/(R*T))))));
	const double calc_i_Ca_L_K_ds = ((((FrICa*P_CaK*P_Ca_L*d_old_*f_old_*f2ds_old_*(V_old_-5.0000000000e+01)*F)/(R*T))/(1.0000000000e+00-exp((((-(V_old_-5.0000000000e+01))*F)/(R*T)))))*((K_i_old_*exp(((5.0000000000e+01*F)/(R*T))))-(K_o*exp((((-(V_old_-5.0000000000e+01))*F)/(R*T))))));
	const double calc_i_Ca_L_Na_ds = ((((FrICa*P_CaNa*P_Ca_L*d_old_*f_old_*f2ds_old_*(V_old_-5.0000000000e+01)*F)/(R*T))/(1.0000000000e+00-exp((((-(V_old_-5.0000000000e+01))*F)/(R*T)))))*((Na_i_old_*exp(((5.0000000000e+01*F)/(R*T))))-(Na_o*exp((((-(V_old_-5.0000000000e+01))*F)/(R*T))))));
	const double calc_i_NaK = ((((i_NaK_max*K_o)/(K_mK+K_o))*Na_i_old_)/(K_mNa+Na_i_old_));
	const double calc_i_NaCa_cyt = (((1.0000000000e+00-FRiNaCa)*k_NaCa*((exp(((gamma*(n_NaCa-2.0000000000e+00)*V_old_*F)/(R*T)))*pow(Na_i_old_,n_NaCa)*Ca_o)-(exp((((gamma-1.0000000000e+00)*(n_NaCa-2.0000000000e+00)*V_old_*F)/(R*T)))*pow(Na_o,n_NaCa)*Ca_i_old_)))/((1.0000000000e+00+(d_NaCa*((Ca_i_old_*pow(Na_o,n_NaCa))+(Ca_o*pow(Na_i_old_,n_NaCa)))))*(1.0000000000e+00+(Ca_i_old_/6.9000000000e-03))));
	const double calc_i_NaCa_ds = ((FRiNaCa*k_NaCa*((exp(((gamma*(n_NaCa-2.0000000000e+00)*V_old_*F)/(R*T)))*pow(Na_i_old_,n_NaCa)*Ca_o)-(exp((((gamma-1.0000000000e+00)*(n_NaCa-2.0000000000e+00)*V_old_*F)/(R*T)))*pow(Na_o,n_NaCa)*Ca_ds_old_)))/((1.0000000000e+00+(d_NaCa*((Ca_ds_old_*pow(Na_o,n_NaCa))+(Ca_o*pow(Na_i_old_,n_NaCa)))))*(1.0000000000e+00+(Ca_ds_old_/6.9000000000e-03))));
	const double calc_K_1 = ((K_cyca*K_xcs)/K_srca);
	const double calc_i_trans = (5.0000000000e+01*(Ca_up_old_-Ca_rel_old_));
	const double calc_i_rel = (((pow((ActFrac_old_/(ActFrac_old_+2.5000000000e-01)),2.0000000000e+00)*K_m_rel)+K_leak_rate)*Ca_rel_old_);
	const double calc_V_Cell = (3.1415926540e+00*pow(radius,2.0000000000e+00)*length);
	const double calc_V_i_ratio = (((1.0000000000e+00-V_e_ratio)-V_up_ratio)-V_rel_ratio);
	const double calc_i_K1 = ((((g_K1*K_o)/(K_o+K_mk1))*(V_old_-calc_E_K))/(1.0000000000e+00+exp(((((V_old_-calc_E_K)-1.0000000000e+01)*F*1.2500000000e+00)/(R*T)))));
	const double calc_i_Kr = (((((g_Kr1*xr1_old_)+(g_Kr2*xr2_old_))*1.0000000000e+00)/(1.0000000000e+00+exp(((V_old_+9.0000000000e+00)/2.2400000000e+01))))*(V_old_-calc_E_K));
	const double calc_i_Ks = (g_Ks*pow(xs_old_,2.0000000000e+00)*(V_old_-calc_E_Ks));
	const double calc_i_Na = (g_Na*pow(m_old_,3.0000000000e+00)*h_old_*(V_old_-calc_E_mh));
	const double calc_i_p_Na = (((g_pna*1.0000000000e+00)/(1.0000000000e+00+exp(((-(V_old_+5.2000000000e+01))/8.0000000000e+00))))*(V_old_-calc_E_Na));
	const double calc_i_b_Na = (g_bna*(V_old_-calc_E_Na));
	const double calc_i_Ca_L = (calc_i_Ca_L_Ca_cyt+calc_i_Ca_L_K_cyt+calc_i_Ca_L_Na_cyt+calc_i_Ca_L_Ca_ds+calc_i_Ca_L_K_ds+calc_i_Ca_L_Na_ds);
	const double calc_i_b_Ca = (g_bca*(V_old_-calc_E_Ca));
	const double calc_i_to = (g_to*(g_tos+(s_old_*(1.0000000000e+00-g_tos)))*r_old_*(V_old_-calc_E_K));
	const double calc_i_NaCa = (calc_i_NaCa_cyt+calc_i_NaCa_ds);
	const double calc_K_2 = (Ca_i_old_+(Ca_up_old_*calc_K_1)+(K_cyca*K_xcs)+K_cyca);
	const double calc_V_i = (calc_V_Cell*calc_V_i_ratio);
	const double calc_i_up = (((Ca_i_old_/calc_K_2)*alpha_up)-(((Ca_up_old_*calc_K_1)/calc_K_2)*beta_up));
	__AGOS->V_lado_direito_= (((-1.0000000000e+00)/Cm)*(calc_i_Stim+calc_i_K1+calc_i_to+calc_i_Kr+calc_i_Ks+calc_i_NaK+calc_i_Na+calc_i_b_Na+calc_i_p_Na+calc_i_Ca_L_Na_cyt+calc_i_Ca_L_Na_ds+calc_i_NaCa_cyt+calc_i_NaCa_ds+calc_i_Ca_L_Ca_cyt+calc_i_Ca_L_Ca_ds+calc_i_Ca_L_K_cyt+calc_i_Ca_L_K_ds+calc_i_b_Ca));
	__AGOS->Na_i_lado_direito_= (((-1.0000000000e+00)/(1.0000000000e+00*calc_V_i*F))*(calc_i_Na+calc_i_p_Na+calc_i_b_Na+(3.0000000000e+00*calc_i_NaK)+(3.0000000000e+00*calc_i_NaCa_cyt)+calc_i_Ca_L_Na_cyt+calc_i_Ca_L_Na_ds));
	__AGOS->K_i_lado_direito_= (((-1.0000000000e+00)/(1.0000000000e+00*calc_V_i*F))*((calc_i_K1+calc_i_Kr+calc_i_Ks+calc_i_Ca_L_K_cyt+calc_i_Ca_L_K_ds+calc_i_to)-(2.0000000000e+00*calc_i_NaK)));
	__AGOS->Ca_i_lado_direito_= (((((((-1.0000000000e+00)/(2.0000000000e+00*1.0000000000e+00*calc_V_i*F))*(((calc_i_Ca_L_Ca_cyt+calc_i_b_Ca)-(2.0000000000e+00*calc_i_NaCa_cyt))-(2.0000000000e+00*calc_i_NaCa_ds)))+(Ca_ds_old_*V_ds_ratio*Kdecay)+((calc_i_rel*V_rel_ratio)/calc_V_i_ratio))-((alpha_Calmod*Ca_i_old_*(Calmod-Ca_Calmod_old_))-(beta_Calmod*Ca_Calmod_old_)))-((alpha_Trop*Ca_i_old_*(Trop-Ca_Trop_old_))-(beta_Trop*Ca_Trop_old_)))-calc_i_up);
	__AGOS->Ca_ds_lado_direito_= ((((-1.0000000000e+00)*calc_i_Ca_L_Ca_ds)/(2.0000000000e+00*1.0000000000e+00*V_ds_ratio*calc_V_i*F))-(Ca_ds_old_*Kdecay));
	__AGOS->Ca_up_lado_direito_= (((calc_V_i_ratio/V_up_ratio)*calc_i_up)-calc_i_trans);
	__AGOS->Ca_rel_lado_direito_= (((V_up_ratio/V_rel_ratio)*calc_i_trans)-calc_i_rel);
	const double calc_alpha_xr1 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(V_old_-5.0000000000e+00))/9.0000000000e+00))));
	const double calc_beta_xr1 = (5.0000000000e-02*exp(((-(V_old_-2.0000000000e+01))/1.5000000000e+01)));
	__AGOS->xr1_lado_direito_= ((calc_alpha_xr1*(1.0000000000e+00-xr1_old_))-(calc_beta_xr1*xr1_old_));
	const double calc_alpha_xr2 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(V_old_-5.0000000000e+00))/9.0000000000e+00))));
	const double calc_beta_xr2 = (4.0000000000e-01*exp((-pow(((V_old_+3.0000000000e+01)/3.0000000000e+01),3.0000000000e+00))));
	__AGOS->xr2_lado_direito_= ((calc_alpha_xr2*(1.0000000000e+00-xr2_old_))-(calc_beta_xr2*xr2_old_));
	const double calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(V_old_-4.0000000000e+01))/9.0000000000e+00))));
	const double calc_beta_xs = (1.0000000000e+00*exp(((-V_old_)/4.5000000000e+01)));
	__AGOS->xs_lado_direito_= ((calc_alpha_xs*(1.0000000000e+00-xs_old_))-(calc_beta_xs*xs_old_));
	const double calc_E0_m = (V_old_+4.1000000000e+01);
	const double calc_beta_m = (8.0000000000e+03*exp(((-5.6000000000e-02)*(V_old_+6.6000000000e+01))));
	const double calc_alpha_m = ((fabs(calc_E0_m)<delta_m))
?(2.0000000000e+03)
:(((2.0000000000e+02*calc_E0_m)/(1.0000000000e+00-exp(((-1.0000000000e-01)*calc_E0_m)))));
	__AGOS->m_lado_direito_= ((calc_alpha_m*(1.0000000000e+00-m_old_))-(calc_beta_m*m_old_));
	const double calc_alpha_h = (2.0000000000e+01*exp(((-1.2500000000e-01)*((V_old_+7.5000000000e+01)-shift_h))));
	const double calc_beta_h = (2.0000000000e+03/(1.0000000000e+00+(3.2000000000e+02*exp(((-1.0000000000e-01)*((V_old_+7.5000000000e+01)-shift_h))))));
	__AGOS->h_lado_direito_= ((calc_alpha_h*(1.0000000000e+00-h_old_))-(calc_beta_h*h_old_));
	const double calc_E0_d = ((V_old_+2.4000000000e+01)-5.0000000000e+00);
	const double calc_alpha_d = ((fabs(calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((3.0000000000e+01*calc_E0_d)/(1.0000000000e+00-exp(((-calc_E0_d)/4.0000000000e+00)))));
	const double calc_beta_d = ((fabs(calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((1.2000000000e+01*calc_E0_d)/(exp((calc_E0_d/1.0000000000e+01))-1.0000000000e+00)));
	__AGOS->d_lado_direito_= (speed_d*((calc_alpha_d*(1.0000000000e+00-d_old_))-(calc_beta_d*d_old_)));
	const double calc_E0_f = (V_old_+3.4000000000e+01);
	const double calc_beta_f = (1.2000000000e+01/(1.0000000000e+00+exp((((-1.0000000000e+00)*(V_old_+3.4000000000e+01))/4.0000000000e+00))));
	const double calc_alpha_f = ((fabs(calc_E0_f)<delta_f))
?(2.5000000000e+01)
:(((6.2500000000e+00*calc_E0_f)/(exp((calc_E0_f/4.0000000000e+00))-1.0000000000e+00)));
	__AGOS->f_lado_direito_= (speed_f*((calc_alpha_f*(1.0000000000e+00-f_old_))-(calc_beta_f*f_old_)));
	const double calc_alpha_s = (3.3000000000e-02*exp(((-V_old_)/1.7000000000e+01)));
	const double calc_beta_s = (3.3000000000e+01/(1.0000000000e+00+exp(((-1.2500000000e-01)*(V_old_+1.0000000000e+01)))));
	__AGOS->s_lado_direito_= ((calc_alpha_s*(1.0000000000e+00-s_old_))-(calc_beta_s*s_old_));
	const double calc_VoltDep = exp((8.0000000000e-02*(V_old_-4.0000000000e+01)));
	const double calc_CaiReg = (Ca_i_old_/(Ca_i_old_+K_m_Ca_cyt));
	const double calc_CadsReg = (Ca_ds_old_/(Ca_ds_old_+K_m_Ca_ds));
	const double calc_SpeedRel = ((V_old_<(-5.0000000000e+01)))
?(5.0000000000e+00)
:(1.0000000000e+00);
	const double calc_PrecFrac = ((1.0000000000e+00-ActFrac_old_)-ProdFrac_old_);
	const double calc_RegBindSite = (calc_CaiReg+((1.0000000000e+00-calc_CaiReg)*calc_CadsReg));
	const double calc_ActRate = ((0.0000000000e+00*calc_VoltDep)+(5.0000000000e+02*pow(calc_RegBindSite,2.0000000000e+00)));
	const double calc_InactRate = (6.0000000000e+01+(5.0000000000e+02*pow(calc_RegBindSite,2.0000000000e+00)));
	__AGOS->ActFrac_lado_direito_= ((calc_PrecFrac*calc_SpeedRel*calc_ActRate)-(ActFrac_old_*calc_SpeedRel*calc_InactRate));
	__AGOS->ProdFrac_lado_direito_= ((ActFrac_old_*calc_SpeedRel*calc_InactRate)-(calc_SpeedRel*1.0000000000e+00*ProdFrac_old_));
	__AGOS->f2_lado_direito_= (1.0000000000e+00-(1.0000000000e+00*((Ca_i_old_/(Km_f2+Ca_i_old_))+f2_old_)));
	__AGOS->f2ds_lado_direito_= (R_decay*(1.0000000000e+00-((Ca_ds_old_/(Km_f2ds+Ca_ds_old_))+f2ds_old_)));
	__AGOS->r_lado_direito_= (3.3300000000e+02*((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+4.0000000000e+00))/5.0000000000e+00))))-r_old_));
	__AGOS->Ca_Calmod_lado_direito_= ((alpha_Calmod*Ca_i_old_*(Calmod-Ca_Calmod_old_))-(beta_Calmod*Ca_Calmod_old_));
	__AGOS->Ca_Trop_lado_direito_= ((alpha_Trop*Ca_i_old_*(Trop-Ca_Trop_old_))-(beta_Trop*Ca_Trop_old_));
} //fim
typedef struct str__rightHandSideFunction{
	void (*function)(Solveode*);
}typ_rightHandSideFunction;
typ_rightHandSideFunction rightHandSideFunction;
typ_rightHandSideFunction forest[15];
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
		time = 0.00000000e+00;
		stim_amplitude = -6.00000000e+00;
		stim_start = 1.00000000e-01;
		stim_end = 1.00000000e+05;
		stim_period = 1.00000000e+00;
		stim_duration = 1.50000000e-03;
		Cm = 9.50000000e-05;
		R = 8.31447200e+03;
		T = 3.10000000e+02;
		F = 9.64853415e+04;
		Na_o = 1.40000000e+02;
		K_o = 4.00000000e+00;
		P_kna = 3.00000000e-02;
		Ca_o = 2.00000000e+00;
		g_K1 = 5.00000000e-01;
		K_mk1 = 1.00000000e+01;
		g_Kr1 = 2.10000000e-03;
		g_Kr2 = 1.30000000e-03;
		g_Ks = 2.60000000e-03;
		g_Na = 2.50000000e+00;
		delta_m = 1.00000000e-05;
		shift_h = 0.00000000e+00;
		g_pna = 4.00000000e-03;
		g_bna = 6.00000000e-04;
		FrICa = 1.00000000e+00;
		P_Ca_L = 1.00000000e-01;
		P_CaK = 2.00000000e-03;
		P_CaNa = 1.00000000e-02;
		speed_d = 3.00000000e+00;
		delta_f = 1.00000000e-04;
		speed_f = 3.00000000e-01;
		Km_f2 = 1.00000000e+05;
		R_decay = 2.00000000e+01;
		Km_f2ds = 1.00000000e-03;
		g_bca = 2.50000000e-04;
		g_to = 5.00000000e-03;
		g_tos = 0.00000000e+00;
		i_NaK_max = 7.00000000e-01;
		K_mK = 1.00000000e+00;
		K_mNa = 4.00000000e+01;
		FRiNaCa = 1.00000000e-03;
		k_NaCa = 5.00000000e-04;
		gamma = 5.00000000e-01;
		n_NaCa = 3.00000000e+00;
		d_NaCa = 0.00000000e+00;
		K_cyca = 3.00000000e-04;
		K_xcs = 4.00000000e-01;
		K_srca = 5.00000000e-01;
		alpha_up = 4.00000000e-01;
		beta_up = 3.00000000e-02;
		K_m_Ca_cyt = 5.00000000e-04;
		K_m_Ca_ds = 1.00000000e-02;
		K_m_rel = 2.50000000e+02;
		K_leak_rate = 5.00000000e-02;
		radius = 1.20000000e-02;
		length = 7.40000000e-02;
		V_e_ratio = 4.00000000e-01;
		V_up_ratio = 1.00000000e-02;
		V_rel_ratio = 1.00000000e-01;
		V_ds_ratio = 1.00000000e-01;
		Kdecay = 1.00000000e+01;
		alpha_Calmod = 1.00000000e+05;
		Calmod = 2.00000000e-02;
		beta_Calmod = 5.00000000e+01;
		alpha_Trop = 1.00000000e+05;
		Trop = 5.00000000e-02;
		beta_Trop = 2.00000000e+02;
		dtime = 0.0; time_vec__ = NULL;
		V = NULL;
		V_ini_ = -9.28493330e+01;
		xr1 = NULL;
		xr1_ini_ = 1.03000000e-05;
		xr2 = NULL;
		xr2_ini_ = 2.00000000e-07;
		xs = NULL;
		xs_ini_ = 1.30200000e-03;
		m = NULL;
		m_ini_ = 1.62030000e-03;
		h = NULL;
		h_ini_ = 9.94403600e-01;
		d = NULL;
		d_ini_ = 0.00000000e+00;
		f = NULL;
		f_ini_ = 1.00000000e+00;
		f2 = NULL;
		f2_ini_ = 9.34919700e-01;
		f2ds = NULL;
		f2ds_ini_ = 9.65195800e-01;
		s = NULL;
		s_ini_ = 9.94864500e-01;
		r = NULL;
		r_ini_ = 0.00000000e+00;
		ActFrac = NULL;
		ActFrac_ini_ = 4.26140000e-03;
		ProdFrac = NULL;
		ProdFrac_ini_ = 4.06815400e-01;
		Na_i = NULL;
		Na_i_ini_ = 7.33212230e+00;
		K_i = NULL;
		K_i_ini_ = 1.36564428e+02;
		Ca_i = NULL;
		Ca_i_ini_ = 1.40000000e-05;
		Ca_ds = NULL;
		Ca_ds_ini_ = 1.88000000e-05;
		Ca_up = NULL;
		Ca_up_ini_ = 4.53188900e-01;
		Ca_rel = NULL;
		Ca_rel_ini_ = 4.48192700e-01;
		Ca_Calmod = NULL;
		Ca_Calmod_ini_ = 5.55500000e-04;
		Ca_Trop = NULL;
		Ca_Trop_ini_ = 3.54200000e-04;
		abstol__ = abs;
		reltol__ = rel;
		it_countx = 0;
	}
	Solveode::~Solveode()
	{
		if(V != NULL) free(V);
		if(xr1 != NULL) free(xr1);
		if(xr2 != NULL) free(xr2);
		if(xs != NULL) free(xs);
		if(m != NULL) free(m);
		if(h != NULL) free(h);
		if(d != NULL) free(d);
		if(f != NULL) free(f);
		if(f2 != NULL) free(f2);
		if(f2ds != NULL) free(f2ds);
		if(s != NULL) free(s);
		if(r != NULL) free(r);
		if(ActFrac != NULL) free(ActFrac);
		if(ProdFrac != NULL) free(ProdFrac);
		if(Na_i != NULL) free(Na_i);
		if(K_i != NULL) free(K_i);
		if(Ca_i != NULL) free(Ca_i);
		if(Ca_ds != NULL) free(Ca_ds);
		if(Ca_up != NULL) free(Ca_up);
		if(Ca_rel != NULL) free(Ca_rel);
		if(Ca_Calmod != NULL) free(Ca_Calmod);
		if(Ca_Trop != NULL) free(Ca_Trop);
	}

	int Solveode::setVariables(int indVariable, double value_new)
	{
		switch(indVariable)
		{
		case 0:		V_ini_ = V_old_= value_new;    break;
		case 1:		xr1_ini_ = xr1_old_= value_new;    break;
		case 2:		xr2_ini_ = xr2_old_= value_new;    break;
		case 3:		xs_ini_ = xs_old_= value_new;    break;
		case 4:		m_ini_ = m_old_= value_new;    break;
		case 5:		h_ini_ = h_old_= value_new;    break;
		case 6:		d_ini_ = d_old_= value_new;    break;
		case 7:		f_ini_ = f_old_= value_new;    break;
		case 8:		f2_ini_ = f2_old_= value_new;    break;
		case 9:		f2ds_ini_ = f2ds_old_= value_new;    break;
		case 10:		s_ini_ = s_old_= value_new;    break;
		case 11:		r_ini_ = r_old_= value_new;    break;
		case 12:		ActFrac_ini_ = ActFrac_old_= value_new;    break;
		case 13:		ProdFrac_ini_ = ProdFrac_old_= value_new;    break;
		case 14:		Na_i_ini_ = Na_i_old_= value_new;    break;
		case 15:		K_i_ini_ = K_i_old_= value_new;    break;
		case 16:		Ca_i_ini_ = Ca_i_old_= value_new;    break;
		case 17:		Ca_ds_ini_ = Ca_ds_old_= value_new;    break;
		case 18:		Ca_up_ini_ = Ca_up_old_= value_new;    break;
		case 19:		Ca_rel_ini_ = Ca_rel_old_= value_new;    break;
		case 20:		Ca_Calmod_ini_ = Ca_Calmod_old_= value_new;    break;
		case 21:		Ca_Trop_ini_ = Ca_Trop_old_= value_new;    break;
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
		case 6:		Cm = value_new;   break;
		case 7:		R = value_new;   break;
		case 8:		T = value_new;   break;
		case 9:		F = value_new;   break;
		case 10:		Na_o = value_new;   break;
		case 11:		K_o = value_new;   break;
		case 12:		P_kna = value_new;   break;
		case 13:		Ca_o = value_new;   break;
		case 14:		g_K1 = value_new;   break;
		case 15:		K_mk1 = value_new;   break;
		case 16:		g_Kr1 = value_new;   break;
		case 17:		g_Kr2 = value_new;   break;
		case 18:		g_Ks = value_new;   break;
		case 19:		g_Na = value_new;   break;
		case 20:		delta_m = value_new;   break;
		case 21:		shift_h = value_new;   break;
		case 22:		g_pna = value_new;   break;
		case 23:		g_bna = value_new;   break;
		case 24:		FrICa = value_new;   break;
		case 25:		P_Ca_L = value_new;   break;
		case 26:		P_CaK = value_new;   break;
		case 27:		P_CaNa = value_new;   break;
		case 28:		speed_d = value_new;   break;
		case 29:		delta_f = value_new;   break;
		case 30:		speed_f = value_new;   break;
		case 31:		Km_f2 = value_new;   break;
		case 32:		R_decay = value_new;   break;
		case 33:		Km_f2ds = value_new;   break;
		case 34:		g_bca = value_new;   break;
		case 35:		g_to = value_new;   break;
		case 36:		g_tos = value_new;   break;
		case 37:		i_NaK_max = value_new;   break;
		case 38:		K_mK = value_new;   break;
		case 39:		K_mNa = value_new;   break;
		case 40:		FRiNaCa = value_new;   break;
		case 41:		k_NaCa = value_new;   break;
		case 42:		gamma = value_new;   break;
		case 43:		n_NaCa = value_new;   break;
		case 44:		d_NaCa = value_new;   break;
		case 45:		K_cyca = value_new;   break;
		case 46:		K_xcs = value_new;   break;
		case 47:		K_srca = value_new;   break;
		case 48:		alpha_up = value_new;   break;
		case 49:		beta_up = value_new;   break;
		case 50:		K_m_Ca_cyt = value_new;   break;
		case 51:		K_m_Ca_ds = value_new;   break;
		case 52:		K_m_rel = value_new;   break;
		case 53:		K_leak_rate = value_new;   break;
		case 54:		radius = value_new;   break;
		case 55:		length = value_new;   break;
		case 56:		V_e_ratio = value_new;   break;
		case 57:		V_up_ratio = value_new;   break;
		case 58:		V_rel_ratio = value_new;   break;
		case 59:		V_ds_ratio = value_new;   break;
		case 60:		Kdecay = value_new;   break;
		case 61:		alpha_Calmod = value_new;   break;
		case 62:		Calmod = value_new;   break;
		case 63:		beta_Calmod = value_new;   break;
		case 64:		alpha_Trop = value_new;   break;
		case 65:		Trop = value_new;   break;
		case 66:		beta_Trop = value_new;   break;
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
		this->xr1_new_ = this->xr1_old_ = this->xr1_ini_;
		this->xr2_new_ = this->xr2_old_ = this->xr2_ini_;
		this->xs_new_ = this->xs_old_ = this->xs_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h_new_ = this->h_old_ = this->h_ini_;
		this->d_new_ = this->d_old_ = this->d_ini_;
		this->f_new_ = this->f_old_ = this->f_ini_;
		this->f2_new_ = this->f2_old_ = this->f2_ini_;
		this->f2ds_new_ = this->f2ds_old_ = this->f2ds_ini_;
		this->s_new_ = this->s_old_ = this->s_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->ActFrac_new_ = this->ActFrac_old_ = this->ActFrac_ini_;
		this->ProdFrac_new_ = this->ProdFrac_old_ = this->ProdFrac_ini_;
		this->Na_i_new_ = this->Na_i_old_ = this->Na_i_ini_;
		this->K_i_new_ = this->K_i_old_ = this->K_i_ini_;
		this->Ca_i_new_ = this->Ca_i_old_ = this->Ca_i_ini_;
		this->Ca_ds_new_ = this->Ca_ds_old_ = this->Ca_ds_ini_;
		this->Ca_up_new_ = this->Ca_up_old_ = this->Ca_up_ini_;
		this->Ca_rel_new_ = this->Ca_rel_old_ = this->Ca_rel_ini_;
		this->Ca_Calmod_new_ = this->Ca_Calmod_old_ = this->Ca_Calmod_ini_;
		this->Ca_Trop_new_ = this->Ca_Trop_old_ = this->Ca_Trop_ini_;
		this->time_new = this->time;
		this->timeSaving = this->time;
		if(savingRate!=0.0)
			this->save_step(fileptr, _EULER_);//save the initial conditions
		while(this->time_new<=finalTime){
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			this->V_new_ = this->dtime*(this->V_lado_direito_) + this->V_old_;
			this->xr1_new_ = this->dtime*(this->xr1_lado_direito_) + this->xr1_old_;
			this->xr2_new_ = this->dtime*(this->xr2_lado_direito_) + this->xr2_old_;
			this->xs_new_ = this->dtime*(this->xs_lado_direito_) + this->xs_old_;
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h_new_ = this->dtime*(this->h_lado_direito_) + this->h_old_;
			this->d_new_ = this->dtime*(this->d_lado_direito_) + this->d_old_;
			this->f_new_ = this->dtime*(this->f_lado_direito_) + this->f_old_;
			this->f2_new_ = this->dtime*(this->f2_lado_direito_) + this->f2_old_;
			this->f2ds_new_ = this->dtime*(this->f2ds_lado_direito_) + this->f2ds_old_;
			this->s_new_ = this->dtime*(this->s_lado_direito_) + this->s_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->ActFrac_new_ = this->dtime*(this->ActFrac_lado_direito_) + this->ActFrac_old_;
			this->ProdFrac_new_ = this->dtime*(this->ProdFrac_lado_direito_) + this->ProdFrac_old_;
			this->Na_i_new_ = this->dtime*(this->Na_i_lado_direito_) + this->Na_i_old_;
			this->K_i_new_ = this->dtime*(this->K_i_lado_direito_) + this->K_i_old_;
			this->Ca_i_new_ = this->dtime*(this->Ca_i_lado_direito_) + this->Ca_i_old_;
			this->Ca_ds_new_ = this->dtime*(this->Ca_ds_lado_direito_) + this->Ca_ds_old_;
			this->Ca_up_new_ = this->dtime*(this->Ca_up_lado_direito_) + this->Ca_up_old_;
			this->Ca_rel_new_ = this->dtime*(this->Ca_rel_lado_direito_) + this->Ca_rel_old_;
			this->Ca_Calmod_new_ = this->dtime*(this->Ca_Calmod_lado_direito_) + this->Ca_Calmod_old_;
			this->Ca_Trop_new_ = this->dtime*(this->Ca_Trop_lado_direito_) + this->Ca_Trop_old_;
			//save results on a file
			if(savingRate!=0.0)
			{
				this->save_step(fileptr, _EULER_);
			}
		this->V_old_ = this->V_new_;
		this->xr1_old_ = this->xr1_new_;
		this->xr2_old_ = this->xr2_new_;
		this->xs_old_ = this->xs_new_;
		this->m_old_ = this->m_new_;
		this->h_old_ = this->h_new_;
		this->d_old_ = this->d_new_;
		this->f_old_ = this->f_new_;
		this->f2_old_ = this->f2_new_;
		this->f2ds_old_ = this->f2ds_new_;
		this->s_old_ = this->s_new_;
		this->r_old_ = this->r_new_;
		this->ActFrac_old_ = this->ActFrac_new_;
		this->ProdFrac_old_ = this->ProdFrac_new_;
		this->Na_i_old_ = this->Na_i_new_;
		this->K_i_old_ = this->K_i_new_;
		this->Ca_i_old_ = this->Ca_i_new_;
		this->Ca_ds_old_ = this->Ca_ds_new_;
		this->Ca_up_old_ = this->Ca_up_new_;
		this->Ca_rel_old_ = this->Ca_rel_new_;
		this->Ca_Calmod_old_ = this->Ca_Calmod_new_;
		this->Ca_Trop_old_ = this->Ca_Trop_new_;
		}
	}
	void Solveode::rungeKutta2ndOrder(double finalTime, FILE *fileptr){
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->xr1_new_ = this->xr1_old_ = this->xr1_ini_;
		this->xr2_new_ = this->xr2_old_ = this->xr2_ini_;
		this->xs_new_ = this->xs_old_ = this->xs_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h_new_ = this->h_old_ = this->h_ini_;
		this->d_new_ = this->d_old_ = this->d_ini_;
		this->f_new_ = this->f_old_ = this->f_ini_;
		this->f2_new_ = this->f2_old_ = this->f2_ini_;
		this->f2ds_new_ = this->f2ds_old_ = this->f2ds_ini_;
		this->s_new_ = this->s_old_ = this->s_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->ActFrac_new_ = this->ActFrac_old_ = this->ActFrac_ini_;
		this->ProdFrac_new_ = this->ProdFrac_old_ = this->ProdFrac_ini_;
		this->Na_i_new_ = this->Na_i_old_ = this->Na_i_ini_;
		this->K_i_new_ = this->K_i_old_ = this->K_i_ini_;
		this->Ca_i_new_ = this->Ca_i_old_ = this->Ca_i_ini_;
		this->Ca_ds_new_ = this->Ca_ds_old_ = this->Ca_ds_ini_;
		this->Ca_up_new_ = this->Ca_up_old_ = this->Ca_up_ini_;
		this->Ca_rel_new_ = this->Ca_rel_old_ = this->Ca_rel_ini_;
		this->Ca_Calmod_new_ = this->Ca_Calmod_old_ = this->Ca_Calmod_ini_;
		this->Ca_Trop_new_ = this->Ca_Trop_old_ = this->Ca_Trop_ini_;
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
			this->xr1_new_ = this->dtime*(this->xr1_lado_direito_) + this->xr1_old_;
			this->xr2_new_ = this->dtime*(this->xr2_lado_direito_) + this->xr2_old_;
			this->xs_new_ = this->dtime*(this->xs_lado_direito_) + this->xs_old_;
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h_new_ = this->dtime*(this->h_lado_direito_) + this->h_old_;
			this->d_new_ = this->dtime*(this->d_lado_direito_) + this->d_old_;
			this->f_new_ = this->dtime*(this->f_lado_direito_) + this->f_old_;
			this->f2_new_ = this->dtime*(this->f2_lado_direito_) + this->f2_old_;
			this->f2ds_new_ = this->dtime*(this->f2ds_lado_direito_) + this->f2ds_old_;
			this->s_new_ = this->dtime*(this->s_lado_direito_) + this->s_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->ActFrac_new_ = this->dtime*(this->ActFrac_lado_direito_) + this->ActFrac_old_;
			this->ProdFrac_new_ = this->dtime*(this->ProdFrac_lado_direito_) + this->ProdFrac_old_;
			this->Na_i_new_ = this->dtime*(this->Na_i_lado_direito_) + this->Na_i_old_;
			this->K_i_new_ = this->dtime*(this->K_i_lado_direito_) + this->K_i_old_;
			this->Ca_i_new_ = this->dtime*(this->Ca_i_lado_direito_) + this->Ca_i_old_;
			this->Ca_ds_new_ = this->dtime*(this->Ca_ds_lado_direito_) + this->Ca_ds_old_;
			this->Ca_up_new_ = this->dtime*(this->Ca_up_lado_direito_) + this->Ca_up_old_;
			this->Ca_rel_new_ = this->dtime*(this->Ca_rel_lado_direito_) + this->Ca_rel_old_;
			this->Ca_Calmod_new_ = this->dtime*(this->Ca_Calmod_lado_direito_) + this->Ca_Calmod_old_;
			this->Ca_Trop_new_ = this->dtime*(this->Ca_Trop_lado_direito_) + this->Ca_Trop_old_;
			//stores the old variables in a vector
			for(int i=0;i<numEDO;i++){
				edos_old_aux_[i] = this->getVariables(i);
				edos_rightside_aux_[i] = this->getLadoDireito(i);
			}
			//steps one iteration ahead;
			this->V_old_ = this->V_new_;
			this->xr1_old_ = this->xr1_new_;
			this->xr2_old_ = this->xr2_new_;
			this->xs_old_ = this->xs_new_;
			this->m_old_ = this->m_new_;
			this->h_old_ = this->h_new_;
			this->d_old_ = this->d_new_;
			this->f_old_ = this->f_new_;
			this->f2_old_ = this->f2_new_;
			this->f2ds_old_ = this->f2ds_new_;
			this->s_old_ = this->s_new_;
			this->r_old_ = this->r_new_;
			this->ActFrac_old_ = this->ActFrac_new_;
			this->ProdFrac_old_ = this->ProdFrac_new_;
			this->Na_i_old_ = this->Na_i_new_;
			this->K_i_old_ = this->K_i_new_;
			this->Ca_i_old_ = this->Ca_i_new_;
			this->Ca_ds_old_ = this->Ca_ds_new_;
			this->Ca_up_old_ = this->Ca_up_new_;
			this->Ca_rel_old_ = this->Ca_rel_new_;
			this->Ca_Calmod_old_ = this->Ca_Calmod_new_;
			this->Ca_Trop_old_ = this->Ca_Trop_new_;
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			//computes the runge kutta second order method
			this->V_new_ = ( this->V_lado_direito_ + edos_rightside_aux_[0] ) * this->dtime/2 + edos_old_aux_[0];
			this->xr1_new_ = ( this->xr1_lado_direito_ + edos_rightside_aux_[1] ) * this->dtime/2 + edos_old_aux_[1];
			this->xr2_new_ = ( this->xr2_lado_direito_ + edos_rightside_aux_[2] ) * this->dtime/2 + edos_old_aux_[2];
			this->xs_new_ = ( this->xs_lado_direito_ + edos_rightside_aux_[3] ) * this->dtime/2 + edos_old_aux_[3];
			this->m_new_ = ( this->m_lado_direito_ + edos_rightside_aux_[4] ) * this->dtime/2 + edos_old_aux_[4];
			this->h_new_ = ( this->h_lado_direito_ + edos_rightside_aux_[5] ) * this->dtime/2 + edos_old_aux_[5];
			this->d_new_ = ( this->d_lado_direito_ + edos_rightside_aux_[6] ) * this->dtime/2 + edos_old_aux_[6];
			this->f_new_ = ( this->f_lado_direito_ + edos_rightside_aux_[7] ) * this->dtime/2 + edos_old_aux_[7];
			this->f2_new_ = ( this->f2_lado_direito_ + edos_rightside_aux_[8] ) * this->dtime/2 + edos_old_aux_[8];
			this->f2ds_new_ = ( this->f2ds_lado_direito_ + edos_rightside_aux_[9] ) * this->dtime/2 + edos_old_aux_[9];
			this->s_new_ = ( this->s_lado_direito_ + edos_rightside_aux_[10] ) * this->dtime/2 + edos_old_aux_[10];
			this->r_new_ = ( this->r_lado_direito_ + edos_rightside_aux_[11] ) * this->dtime/2 + edos_old_aux_[11];
			this->ActFrac_new_ = ( this->ActFrac_lado_direito_ + edos_rightside_aux_[12] ) * this->dtime/2 + edos_old_aux_[12];
			this->ProdFrac_new_ = ( this->ProdFrac_lado_direito_ + edos_rightside_aux_[13] ) * this->dtime/2 + edos_old_aux_[13];
			this->Na_i_new_ = ( this->Na_i_lado_direito_ + edos_rightside_aux_[14] ) * this->dtime/2 + edos_old_aux_[14];
			this->K_i_new_ = ( this->K_i_lado_direito_ + edos_rightside_aux_[15] ) * this->dtime/2 + edos_old_aux_[15];
			this->Ca_i_new_ = ( this->Ca_i_lado_direito_ + edos_rightside_aux_[16] ) * this->dtime/2 + edos_old_aux_[16];
			this->Ca_ds_new_ = ( this->Ca_ds_lado_direito_ + edos_rightside_aux_[17] ) * this->dtime/2 + edos_old_aux_[17];
			this->Ca_up_new_ = ( this->Ca_up_lado_direito_ + edos_rightside_aux_[18] ) * this->dtime/2 + edos_old_aux_[18];
			this->Ca_rel_new_ = ( this->Ca_rel_lado_direito_ + edos_rightside_aux_[19] ) * this->dtime/2 + edos_old_aux_[19];
			this->Ca_Calmod_new_ = ( this->Ca_Calmod_lado_direito_ + edos_rightside_aux_[20] ) * this->dtime/2 + edos_old_aux_[20];
			this->Ca_Trop_new_ = ( this->Ca_Trop_lado_direito_ + edos_rightside_aux_[21] ) * this->dtime/2 + edos_old_aux_[21];
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
			this->xr1_old_ = this->xr1_new_;
			this->xr2_old_ = this->xr2_new_;
			this->xs_old_ = this->xs_new_;
			this->m_old_ = this->m_new_;
			this->h_old_ = this->h_new_;
			this->d_old_ = this->d_new_;
			this->f_old_ = this->f_new_;
			this->f2_old_ = this->f2_new_;
			this->f2ds_old_ = this->f2ds_new_;
			this->s_old_ = this->s_new_;
			this->r_old_ = this->r_new_;
			this->ActFrac_old_ = this->ActFrac_new_;
			this->ProdFrac_old_ = this->ProdFrac_new_;
			this->Na_i_old_ = this->Na_i_new_;
			this->K_i_old_ = this->K_i_new_;
			this->Ca_i_old_ = this->Ca_i_new_;
			this->Ca_ds_old_ = this->Ca_ds_new_;
			this->Ca_up_old_ = this->Ca_up_new_;
			this->Ca_rel_old_ = this->Ca_rel_new_;
			this->Ca_Calmod_old_ = this->Ca_Calmod_new_;
			this->Ca_Trop_old_ = this->Ca_Trop_new_;
		}
	}
	void Solveode::addt(double finalTime, FILE *fileptr){
		double _tolerances_[numEDO];
		double _aux_tol=0.0;
		double maxDt = this->dtime, minDt = this->dtime;
		//initializes the variables
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->xr1_new_ = this->xr1_old_ = this->xr1_ini_;
		this->xr2_new_ = this->xr2_old_ = this->xr2_ini_;
		this->xs_new_ = this->xs_old_ = this->xs_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h_new_ = this->h_old_ = this->h_ini_;
		this->d_new_ = this->d_old_ = this->d_ini_;
		this->f_new_ = this->f_old_ = this->f_ini_;
		this->f2_new_ = this->f2_old_ = this->f2_ini_;
		this->f2ds_new_ = this->f2ds_old_ = this->f2ds_ini_;
		this->s_new_ = this->s_old_ = this->s_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->ActFrac_new_ = this->ActFrac_old_ = this->ActFrac_ini_;
		this->ProdFrac_new_ = this->ProdFrac_old_ = this->ProdFrac_ini_;
		this->Na_i_new_ = this->Na_i_old_ = this->Na_i_ini_;
		this->K_i_new_ = this->K_i_old_ = this->K_i_ini_;
		this->Ca_i_new_ = this->Ca_i_old_ = this->Ca_i_ini_;
		this->Ca_ds_new_ = this->Ca_ds_old_ = this->Ca_ds_ini_;
		this->Ca_up_new_ = this->Ca_up_old_ = this->Ca_up_ini_;
		this->Ca_rel_new_ = this->Ca_rel_old_ = this->Ca_rel_ini_;
		this->Ca_Calmod_new_ = this->Ca_Calmod_old_ = this->Ca_Calmod_ini_;
		this->Ca_Trop_new_ = this->Ca_Trop_old_ = this->Ca_Trop_ini_;
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
		int desc=0;
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
						this->xr1_new_ = edos_new_euler_[1];
						this->xr2_new_ = edos_new_euler_[2];
						this->xs_new_ = edos_new_euler_[3];
						this->m_new_ = edos_new_euler_[4];
						this->h_new_ = edos_new_euler_[5];
						this->d_new_ = edos_new_euler_[6];
						this->f_new_ = edos_new_euler_[7];
						this->f2_new_ = edos_new_euler_[8];
						this->f2ds_new_ = edos_new_euler_[9];
						this->s_new_ = edos_new_euler_[10];
						this->r_new_ = edos_new_euler_[11];
						this->ActFrac_new_ = edos_new_euler_[12];
						this->ProdFrac_new_ = edos_new_euler_[13];
						this->Na_i_new_ = edos_new_euler_[14];
						this->K_i_new_ = edos_new_euler_[15];
						this->Ca_i_new_ = edos_new_euler_[16];
						this->Ca_ds_new_ = edos_new_euler_[17];
						this->Ca_up_new_ = edos_new_euler_[18];
						this->Ca_rel_new_ = edos_new_euler_[19];
						this->Ca_Calmod_new_ = edos_new_euler_[20];
						this->Ca_Trop_new_ = edos_new_euler_[21];
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
		printf("%d\n", desc);
		if(savingRate!=0){
		printf("Dt max: %e dt min %e, %e %d\n", maxDt, minDt, _soma_/_cont_, _cont_);
		}
	}
	void Solveode::addt2(double finalTime, FILE *fileptr){
		const double _beta_safety_ = 0.8;
		double _tolerances_[numEDO];
		double _aux_tol=0.0;
		double maxDt = this->dtime, minDt = this->dtime;
		//initializes the variables
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->xr1_new_ = this->xr1_old_ = this->xr1_ini_;
		this->xr2_new_ = this->xr2_old_ = this->xr2_ini_;
		this->xs_new_ = this->xs_old_ = this->xs_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h_new_ = this->h_old_ = this->h_ini_;
		this->d_new_ = this->d_old_ = this->d_ini_;
		this->f_new_ = this->f_old_ = this->f_ini_;
		this->f2_new_ = this->f2_old_ = this->f2_ini_;
		this->f2ds_new_ = this->f2ds_old_ = this->f2ds_ini_;
		this->s_new_ = this->s_old_ = this->s_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->ActFrac_new_ = this->ActFrac_old_ = this->ActFrac_ini_;
		this->ProdFrac_new_ = this->ProdFrac_old_ = this->ProdFrac_ini_;
		this->Na_i_new_ = this->Na_i_old_ = this->Na_i_ini_;
		this->K_i_new_ = this->K_i_old_ = this->K_i_ini_;
		this->Ca_i_new_ = this->Ca_i_old_ = this->Ca_i_ini_;
		this->Ca_ds_new_ = this->Ca_ds_old_ = this->Ca_ds_ini_;
		this->Ca_up_new_ = this->Ca_up_old_ = this->Ca_up_ini_;
		this->Ca_rel_new_ = this->Ca_rel_old_ = this->Ca_rel_ini_;
		this->Ca_Calmod_new_ = this->Ca_Calmod_old_ = this->Ca_Calmod_ini_;
		this->Ca_Trop_new_ = this->Ca_Trop_old_ = this->Ca_Trop_ini_;
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
		int desc =0;
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
				for(int i=0; i<numEDO; i++)
					this->setVariables(i, edos_old_aux_[i]);
				desc++;
				//throw the results away and compute again
				}else{//it accepts the solutions
					if(savingRate!=0.0){
					//restore the previous value of old variables
					for(int i=0; i<numEDO; i++)
						this->setVariables(i, edos_old_aux_[i]);
					this->V_new_ = edos_new_euler_[0];
					this->xr1_new_ = edos_new_euler_[1];
					this->xr2_new_ = edos_new_euler_[2];
					this->xs_new_ = edos_new_euler_[3];
					this->m_new_ = edos_new_euler_[4];
					this->h_new_ = edos_new_euler_[5];
					this->d_new_ = edos_new_euler_[6];
					this->f_new_ = edos_new_euler_[7];
					this->f2_new_ = edos_new_euler_[8];
					this->f2ds_new_ = edos_new_euler_[9];
					this->s_new_ = edos_new_euler_[10];
					this->r_new_ = edos_new_euler_[11];
					this->ActFrac_new_ = edos_new_euler_[12];
					this->ProdFrac_new_ = edos_new_euler_[13];
					this->Na_i_new_ = edos_new_euler_[14];
					this->K_i_new_ = edos_new_euler_[15];
					this->Ca_i_new_ = edos_new_euler_[16];
					this->Ca_ds_new_ = edos_new_euler_[17];
					this->Ca_up_new_ = edos_new_euler_[18];
					this->Ca_rel_new_ = edos_new_euler_[19];
					this->Ca_Calmod_new_ = edos_new_euler_[20];
					this->Ca_Trop_new_ = edos_new_euler_[21];
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
		printf("%d\n", desc);
		if(savingRate!=0){
		printf("Dt max: %e dt min %e, %e %d\n", maxDt, minDt, _soma_/_cont_ , _cont_);
		}
	}
	void Solveode::euler_OMP(double finalTime, FILE *fileptr, int nThreads){
		omp_set_num_threads(nThreads);
		double *__NEW_ = (double*)malloc(sizeof(double)*numEDO);
		double *__OLD_ = (double*)malloc(sizeof(double)*numEDO);
		double *temp;
// 		int tempoespera=0;
		this->time_new = this->time;
		this->timeSaving = this->time;
		this->previous_dt = this->dtime;
		#pragma omp parallel firstprivate(__OLD_, __NEW_, temp)
// 		shared(tempoespera)
		{
		    
// 			tempoespera=0;
			//private parameters
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -6.0000000000e+00,  _prvt_stim_start = 1.0000000000e-01,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 1.0000000000e+00,  _prvt_stim_duration = 1.5000000000e-03,  _prvt_Cm = 9.5000000000e-05,  _prvt_R = 8.3144720000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_F = 9.6485341500e+04,  _prvt_Na_o = 1.4000000000e+02,  _prvt_K_o = 4.0000000000e+00,  _prvt_P_kna = 3.0000000000e-02,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_g_K1 = 5.0000000000e-01,  _prvt_K_mk1 = 1.0000000000e+01,  _prvt_g_Kr1 = 2.1000000000e-03,  _prvt_g_Kr2 = 1.3000000000e-03,  _prvt_g_Ks = 2.6000000000e-03,  _prvt_g_Na = 2.5000000000e+00,  _prvt_delta_m = 1.0000000000e-05,  _prvt_shift_h = 0.0000000000e+00,  _prvt_g_pna = 4.0000000000e-03,  _prvt_g_bna = 6.0000000000e-04,  _prvt_FrICa = 1.0000000000e+00,  _prvt_P_Ca_L = 1.0000000000e-01,  _prvt_P_CaK = 2.0000000000e-03,  _prvt_P_CaNa = 1.0000000000e-02,  _prvt_speed_d = 3.0000000000e+00,  _prvt_delta_f = 1.0000000000e-04,  _prvt_speed_f = 3.0000000000e-01,  _prvt_Km_f2 = 1.0000000000e+05,  _prvt_R_decay = 2.0000000000e+01,  _prvt_Km_f2ds = 1.0000000000e-03,  _prvt_g_bca = 2.5000000000e-04,  _prvt_g_to = 5.0000000000e-03,  _prvt_g_tos = 0.0000000000e+00,  _prvt_i_NaK_max = 7.0000000000e-01,  _prvt_K_mK = 1.0000000000e+00,  _prvt_K_mNa = 4.0000000000e+01,  _prvt_FRiNaCa = 1.0000000000e-03,  _prvt_k_NaCa = 5.0000000000e-04,  _prvt_gamma = 5.0000000000e-01,  _prvt_n_NaCa = 3.0000000000e+00,  _prvt_d_NaCa = 0.0000000000e+00,  _prvt_K_cyca = 3.0000000000e-04,  _prvt_K_xcs = 4.0000000000e-01,  _prvt_K_srca = 5.0000000000e-01,  _prvt_alpha_up = 4.0000000000e-01,  _prvt_beta_up = 3.0000000000e-02,  _prvt_K_m_Ca_cyt = 5.0000000000e-04,  _prvt_K_m_Ca_ds = 1.0000000000e-02,  _prvt_K_m_rel = 2.5000000000e+02,  _prvt_K_leak_rate = 5.0000000000e-02,  _prvt_radius = 1.2000000000e-02,  _prvt_length = 7.4000000000e-02,  _prvt_V_e_ratio = 4.0000000000e-01,  _prvt_V_up_ratio = 1.0000000000e-02,  _prvt_V_rel_ratio = 1.0000000000e-01,  _prvt_V_ds_ratio = 1.0000000000e-01,  _prvt_Kdecay = 1.0000000000e+01,  _prvt_alpha_Calmod = 1.0000000000e+05,  _prvt_Calmod = 2.0000000000e-02,  _prvt_beta_Calmod = 5.0000000000e+01,  _prvt_alpha_Trop = 1.0000000000e+05,  _prvt_Trop = 5.0000000000e-02,  _prvt_beta_Trop = 2.0000000000e+02, 
			//private aux variables
			 _prvt_calc_i_Stim=0,  _prvt_calc_E_Na=0,  _prvt_calc_E_K=0,  _prvt_calc_E_Ks=0,  _prvt_calc_E_Ca=0,  _prvt_calc_E_mh=0,  _prvt_calc_i_K1=0,  _prvt_calc_i_Kr=0,  _prvt_calc_alpha_xr1=0,  _prvt_calc_beta_xr1=0,  _prvt_calc_alpha_xr2=0,  _prvt_calc_beta_xr2=0,  _prvt_calc_i_Ks=0,  _prvt_calc_alpha_xs=0,  _prvt_calc_beta_xs=0,  _prvt_calc_i_Na=0,  _prvt_calc_E0_m=0,  _prvt_calc_alpha_m=0,  _prvt_calc_beta_m=0,  _prvt_calc_alpha_h=0,  _prvt_calc_beta_h=0,  _prvt_calc_i_p_Na=0,  _prvt_calc_i_b_Na=0,  _prvt_calc_i_Ca_L_Ca_cyt=0,  _prvt_calc_i_Ca_L_K_cyt=0,  _prvt_calc_i_Ca_L_Na_cyt=0,  _prvt_calc_i_Ca_L_Ca_ds=0,  _prvt_calc_i_Ca_L_K_ds=0,  _prvt_calc_i_Ca_L_Na_ds=0,  _prvt_calc_i_Ca_L=0,  _prvt_calc_E0_d=0,  _prvt_calc_alpha_d=0,  _prvt_calc_beta_d=0,  _prvt_calc_E0_f=0,  _prvt_calc_alpha_f=0,  _prvt_calc_beta_f=0,  _prvt_calc_i_b_Ca=0,  _prvt_calc_i_to=0,  _prvt_calc_alpha_s=0,  _prvt_calc_beta_s=0,  _prvt_calc_i_NaK=0,  _prvt_calc_i_NaCa_cyt=0,  _prvt_calc_i_NaCa_ds=0,  _prvt_calc_i_NaCa=0,  _prvt_calc_K_1=0,  _prvt_calc_K_2=0,  _prvt_calc_i_up=0,  _prvt_calc_i_trans=0,  _prvt_calc_VoltDep=0,  _prvt_calc_CaiReg=0,  _prvt_calc_CadsReg=0,  _prvt_calc_RegBindSite=0,  _prvt_calc_ActRate=0,  _prvt_calc_InactRate=0,  _prvt_calc_SpeedRel=0,  _prvt_calc_PrecFrac=0,  _prvt_calc_i_rel=0,  _prvt_calc_V_Cell=0,  _prvt_calc_V_i_ratio=0,  _prvt_calc_V_i=0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_xr1_lado_direito_,  _prvt_xr2_lado_direito_,  _prvt_xs_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h_lado_direito_,  _prvt_d_lado_direito_,  _prvt_f_lado_direito_,  _prvt_f2_lado_direito_,  _prvt_f2ds_lado_direito_,  _prvt_s_lado_direito_,  _prvt_r_lado_direito_,  _prvt_ActFrac_lado_direito_,  _prvt_ProdFrac_lado_direito_,  _prvt_Na_i_lado_direito_,  _prvt_K_i_lado_direito_,  _prvt_Ca_i_lado_direito_,  _prvt_Ca_ds_lado_direito_,  _prvt_Ca_up_lado_direito_,  _prvt_Ca_rel_lado_direito_,  _prvt_Ca_Calmod_lado_direito_,  _prvt_Ca_Trop_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate;
			__NEW_[0] = __OLD_[0] = -9.2849333000e+01;
			__NEW_[1] = __OLD_[1] = 1.0300000000e-05;
			__NEW_[2] = __OLD_[2] = 2.0000000000e-07;
			__NEW_[3] = __OLD_[3] = 1.3020000000e-03;
			__NEW_[4] = __OLD_[4] = 1.6203000000e-03;
			__NEW_[5] = __OLD_[5] = 9.9440360000e-01;
			__NEW_[6] = __OLD_[6] = 0.0000000000e+00;
			__NEW_[7] = __OLD_[7] = 1.0000000000e+00;
			__NEW_[8] = __OLD_[8] = 9.3491970000e-01;
			__NEW_[9] = __OLD_[9] = 9.6519580000e-01;
			__NEW_[10] = __OLD_[10] = 9.9486450000e-01;
			__NEW_[11] = __OLD_[11] = 0.0000000000e+00;
			__NEW_[12] = __OLD_[12] = 4.2614000000e-03;
			__NEW_[13] = __OLD_[13] = 4.0681540000e-01;
			__NEW_[14] = __OLD_[14] = 7.3321223000e+00;
			__NEW_[15] = __OLD_[15] = 1.3656442810e+02;
			__NEW_[16] = __OLD_[16] = 1.4000000000e-05;
			__NEW_[17] = __OLD_[17] = 1.8800000000e-05;
			__NEW_[18] = __OLD_[18] = 4.5318890000e-01;
			__NEW_[19] = __OLD_[19] = 4.4819270000e-01;
			__NEW_[20] = __OLD_[20] = 5.5550000000e-04;
			__NEW_[21] = __OLD_[21] = 3.5420000000e-04;
			int *_prvt_tree_thread = tree_thread;
			while(_prvt_time_new<=_prvt_finalTime){
				_prvt_time_new += _prvt_dtime;
				if(omp_get_thread_num()==tree_thread[0])
				{
					_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[14])));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[15])));
					_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[15]+(_prvt_P_kna*__OLD_[14])))));
					_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[16])));
					_prvt_calc_E_mh = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_Na_o+(1.2000000000e-01*_prvt_K_o))/(__OLD_[14]+(1.2000000000e-01*__OLD_[15])))));
					_prvt_calc_i_Ca_L_Ca_cyt = (((((1.0000000000e+00-_prvt_FrICa)*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_K_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Na_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Ca_ds = ((((_prvt_FrICa*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_K_ds = ((((_prvt_FrICa*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Na_ds = ((((_prvt_FrICa*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_K_o)/(_prvt_K_mK+_prvt_K_o))*__OLD_[14])/(_prvt_K_mNa+__OLD_[14]));
					_prvt_calc_i_NaCa_cyt = (((1.0000000000e+00-_prvt_FRiNaCa)*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[16])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[16]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[16]/6.9000000000e-03))));
					_prvt_calc_i_NaCa_ds = ((_prvt_FRiNaCa*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[17])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[17]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[17]/6.9000000000e-03))));
					_prvt_calc_K_1 = ((_prvt_K_cyca*_prvt_K_xcs)/_prvt_K_srca);
					_prvt_calc_i_trans = (5.0000000000e+01*(__OLD_[18]-__OLD_[19]));
					_prvt_calc_i_rel = (((pow((__OLD_[12]/(__OLD_[12]+2.5000000000e-01)),2.0000000000e+00)*_prvt_K_m_rel)+_prvt_K_leak_rate)*__OLD_[19]);
					_prvt_calc_V_Cell = (3.1415926540e+00*pow(_prvt_radius,2.0000000000e+00)*_prvt_length);
					_prvt_calc_V_i_ratio = (((1.0000000000e+00-_prvt_V_e_ratio)-_prvt_V_up_ratio)-_prvt_V_rel_ratio);
					_prvt_calc_i_K1 = ((((_prvt_g_K1*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk1))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01)*_prvt_F*1.2500000000e+00)/(_prvt_R*_prvt_T)))));
					_prvt_calc_i_Kr = (((((_prvt_g_Kr1*__OLD_[1])+(_prvt_g_Kr2*__OLD_[2]))*1.0000000000e+00)/(1.0000000000e+00+exp(((__OLD_[0]+9.0000000000e+00)/2.2400000000e+01))))*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
					_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*(__OLD_[0]-_prvt_calc_E_mh));
					_prvt_calc_i_p_Na = (((_prvt_g_pna*1.0000000000e+00)/(1.0000000000e+00+exp(((-(__OLD_[0]+5.2000000000e+01))/8.0000000000e+00))))*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Ca_L = (_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_Ca_L_Na_ds);
					_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_to = (_prvt_g_to*(_prvt_g_tos+(__OLD_[10]*(1.0000000000e+00-_prvt_g_tos)))*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_NaCa = (_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds);
					_prvt_calc_K_2 = (__OLD_[16]+(__OLD_[18]*_prvt_calc_K_1)+(_prvt_K_cyca*_prvt_K_xcs)+_prvt_K_cyca);
					_prvt_calc_V_i = (_prvt_calc_V_Cell*_prvt_calc_V_i_ratio);
					_prvt_calc_i_up = (((__OLD_[16]/_prvt_calc_K_2)*_prvt_alpha_up)-(((__OLD_[18]*_prvt_calc_K_1)/_prvt_calc_K_2)*_prvt_beta_up));
					_prvt_V_lado_direito_= (((-1.0000000000e+00)/_prvt_Cm)*(_prvt_calc_i_Stim+_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_p_Na+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds+_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds+_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_b_Ca));
					__NEW_[0]= _prvt_V_lado_direito_ * _prvt_dtime + __OLD_[0];
					_prvt_Na_i_lado_direito_= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(_prvt_calc_i_Na+_prvt_calc_i_p_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa_cyt)+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds));
					__NEW_[14]= _prvt_Na_i_lado_direito_ * _prvt_dtime + __OLD_[14];
					_prvt_K_i_lado_direito_= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*((_prvt_calc_i_K1+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_to)-(2.0000000000e+00*_prvt_calc_i_NaK)));
					__NEW_[15]= _prvt_K_i_lado_direito_ * _prvt_dtime + __OLD_[15];
					_prvt_Ca_i_lado_direito_= (((((((-1.0000000000e+00)/(2.0000000000e+00*1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(((_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_b_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa_cyt))-(2.0000000000e+00*_prvt_calc_i_NaCa_ds)))+(__OLD_[17]*_prvt_V_ds_ratio*_prvt_Kdecay)+((_prvt_calc_i_rel*_prvt_V_rel_ratio)/_prvt_calc_V_i_ratio))-((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20])))-((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21])))-_prvt_calc_i_up);
					__NEW_[16]= _prvt_Ca_i_lado_direito_ * _prvt_dtime + __OLD_[16];
					_prvt_Ca_ds_lado_direito_= ((((-1.0000000000e+00)*_prvt_calc_i_Ca_L_Ca_ds)/(2.0000000000e+00*1.0000000000e+00*_prvt_V_ds_ratio*_prvt_calc_V_i*_prvt_F))-(__OLD_[17]*_prvt_Kdecay));
					__NEW_[17]= _prvt_Ca_ds_lado_direito_ * _prvt_dtime + __OLD_[17];
					_prvt_Ca_up_lado_direito_= (((_prvt_calc_V_i_ratio/_prvt_V_up_ratio)*_prvt_calc_i_up)-_prvt_calc_i_trans);
					__NEW_[18]= _prvt_Ca_up_lado_direito_ * _prvt_dtime + __OLD_[18];
					_prvt_Ca_rel_lado_direito_= (((_prvt_V_up_ratio/_prvt_V_rel_ratio)*_prvt_calc_i_trans)-_prvt_calc_i_rel);
					__NEW_[19]= _prvt_Ca_rel_lado_direito_ * _prvt_dtime + __OLD_[19];
				}
				if(omp_get_thread_num()==tree_thread[1])
				{
					_prvt_calc_alpha_xr1 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
					_prvt_calc_beta_xr1 = (5.0000000000e-02*exp(((-(__OLD_[0]-2.0000000000e+01))/1.5000000000e+01)));
					_prvt_xr1_lado_direito_= ((_prvt_calc_alpha_xr1*(1.0000000000e+00-__OLD_[1]))-(_prvt_calc_beta_xr1*__OLD_[1]));
					__NEW_[1]= _prvt_xr1_lado_direito_ * _prvt_dtime + __OLD_[1];
				}
				if(omp_get_thread_num()==tree_thread[2])
				{
					_prvt_calc_alpha_xr2 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
					_prvt_calc_beta_xr2 = (4.0000000000e-01*exp((-pow(((__OLD_[0]+3.0000000000e+01)/3.0000000000e+01),3.0000000000e+00))));
					_prvt_xr2_lado_direito_= ((_prvt_calc_alpha_xr2*(1.0000000000e+00-__OLD_[2]))-(_prvt_calc_beta_xr2*__OLD_[2]));
					__NEW_[2]= _prvt_xr2_lado_direito_ * _prvt_dtime + __OLD_[2];
				}
				if(omp_get_thread_num()==tree_thread[3])
				{
					_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
					_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
					_prvt_xs_lado_direito_= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[3]))-(_prvt_calc_beta_xs*__OLD_[3]));
					__NEW_[3]= _prvt_xs_lado_direito_ * _prvt_dtime + __OLD_[3];
				}
				if(omp_get_thread_num()==tree_thread[4])
				{
					_prvt_calc_E0_m = (__OLD_[0]+4.1000000000e+01);
					_prvt_calc_beta_m = (8.0000000000e+03*exp(((-5.6000000000e-02)*(__OLD_[0]+6.6000000000e+01))));
					_prvt_calc_alpha_m = ((fabs(_prvt_calc_E0_m)<_prvt_delta_m))
?(2.0000000000e+03)
:(((2.0000000000e+02*_prvt_calc_E0_m)/(1.0000000000e+00-exp(((-1.0000000000e-01)*_prvt_calc_E0_m)))));
					_prvt_m_lado_direito_= ((_prvt_calc_alpha_m*(1.0000000000e+00-__OLD_[4]))-(_prvt_calc_beta_m*__OLD_[4]));
					__NEW_[4]= _prvt_m_lado_direito_ * _prvt_dtime + __OLD_[4];
				}
				if(omp_get_thread_num()==tree_thread[5])
				{
					_prvt_calc_alpha_h = (2.0000000000e+01*exp(((-1.2500000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))));
					_prvt_calc_beta_h = (2.0000000000e+03/(1.0000000000e+00+(3.2000000000e+02*exp(((-1.0000000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))))));
					_prvt_h_lado_direito_= ((_prvt_calc_alpha_h*(1.0000000000e+00-__OLD_[5]))-(_prvt_calc_beta_h*__OLD_[5]));
					__NEW_[5]= _prvt_h_lado_direito_ * _prvt_dtime + __OLD_[5];
				}
				if(omp_get_thread_num()==tree_thread[6])
				{
					_prvt_calc_E0_d = ((__OLD_[0]+2.4000000000e+01)-5.0000000000e+00);
					_prvt_calc_alpha_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((3.0000000000e+01*_prvt_calc_E0_d)/(1.0000000000e+00-exp(((-_prvt_calc_E0_d)/4.0000000000e+00)))));
					_prvt_calc_beta_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((1.2000000000e+01*_prvt_calc_E0_d)/(exp((_prvt_calc_E0_d/1.0000000000e+01))-1.0000000000e+00)));
					_prvt_d_lado_direito_= (_prvt_speed_d*((_prvt_calc_alpha_d*(1.0000000000e+00-__OLD_[6]))-(_prvt_calc_beta_d*__OLD_[6])));
					__NEW_[6]= _prvt_d_lado_direito_ * _prvt_dtime + __OLD_[6];
				}
				if(omp_get_thread_num()==tree_thread[7])
				{
					_prvt_calc_E0_f = (__OLD_[0]+3.4000000000e+01);
					_prvt_calc_beta_f = (1.2000000000e+01/(1.0000000000e+00+exp((((-1.0000000000e+00)*(__OLD_[0]+3.4000000000e+01))/4.0000000000e+00))));
					_prvt_calc_alpha_f = ((fabs(_prvt_calc_E0_f)<_prvt_delta_f))
?(2.5000000000e+01)
:(((6.2500000000e+00*_prvt_calc_E0_f)/(exp((_prvt_calc_E0_f/4.0000000000e+00))-1.0000000000e+00)));
					_prvt_f_lado_direito_= (_prvt_speed_f*((_prvt_calc_alpha_f*(1.0000000000e+00-__OLD_[7]))-(_prvt_calc_beta_f*__OLD_[7])));
					__NEW_[7]= _prvt_f_lado_direito_ * _prvt_dtime + __OLD_[7];
				}
				if(omp_get_thread_num()==tree_thread[8])
				{
					_prvt_calc_alpha_s = (3.3000000000e-02*exp(((-__OLD_[0])/1.7000000000e+01)));
					_prvt_calc_beta_s = (3.3000000000e+01/(1.0000000000e+00+exp(((-1.2500000000e-01)*(__OLD_[0]+1.0000000000e+01)))));
					_prvt_s_lado_direito_= ((_prvt_calc_alpha_s*(1.0000000000e+00-__OLD_[10]))-(_prvt_calc_beta_s*__OLD_[10]));
					__NEW_[10]= _prvt_s_lado_direito_ * _prvt_dtime + __OLD_[10];
				}
				if(omp_get_thread_num()==tree_thread[9])
				{
					_prvt_calc_VoltDep = exp((8.0000000000e-02*(__OLD_[0]-4.0000000000e+01)));
					_prvt_calc_CaiReg = (__OLD_[16]/(__OLD_[16]+_prvt_K_m_Ca_cyt));
					_prvt_calc_CadsReg = (__OLD_[17]/(__OLD_[17]+_prvt_K_m_Ca_ds));
					_prvt_calc_SpeedRel = ((__OLD_[0]<(-5.0000000000e+01)))
?(5.0000000000e+00)
:(1.0000000000e+00);
					_prvt_calc_PrecFrac = ((1.0000000000e+00-__OLD_[12])-__OLD_[13]);
					_prvt_calc_RegBindSite = (_prvt_calc_CaiReg+((1.0000000000e+00-_prvt_calc_CaiReg)*_prvt_calc_CadsReg));
					_prvt_calc_ActRate = ((0.0000000000e+00*_prvt_calc_VoltDep)+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
					_prvt_calc_InactRate = (6.0000000000e+01+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
					_prvt_ActFrac_lado_direito_= ((_prvt_calc_PrecFrac*_prvt_calc_SpeedRel*_prvt_calc_ActRate)-(__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate));
					__NEW_[12]= _prvt_ActFrac_lado_direito_ * _prvt_dtime + __OLD_[12];
					_prvt_ProdFrac_lado_direito_= ((__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate)-(_prvt_calc_SpeedRel*1.0000000000e+00*__OLD_[13]));
					__NEW_[13]= _prvt_ProdFrac_lado_direito_ * _prvt_dtime + __OLD_[13];
				}
				if(omp_get_thread_num()==tree_thread[10])
				{
					_prvt_f2_lado_direito_= (1.0000000000e+00-(1.0000000000e+00*((__OLD_[16]/(_prvt_Km_f2+__OLD_[16]))+__OLD_[8])));
					__NEW_[8]= _prvt_f2_lado_direito_ * _prvt_dtime + __OLD_[8];
				}
				if(omp_get_thread_num()==tree_thread[11])
				{
					_prvt_f2ds_lado_direito_= (_prvt_R_decay*(1.0000000000e+00-((__OLD_[17]/(_prvt_Km_f2ds+__OLD_[17]))+__OLD_[9])));
					__NEW_[9]= _prvt_f2ds_lado_direito_ * _prvt_dtime + __OLD_[9];
				}
				if(omp_get_thread_num()==tree_thread[12])
				{
					_prvt_r_lado_direito_= (3.3300000000e+02*((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+4.0000000000e+00))/5.0000000000e+00))))-__OLD_[11]));
					__NEW_[11]= _prvt_r_lado_direito_ * _prvt_dtime + __OLD_[11];
				}
				if(omp_get_thread_num()==tree_thread[13])
				{
					_prvt_Ca_Calmod_lado_direito_= ((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20]));
					__NEW_[20]= _prvt_Ca_Calmod_lado_direito_ * _prvt_dtime + __OLD_[20];
				}
				if(omp_get_thread_num()==tree_thread[14])
				{
					_prvt_Ca_Trop_lado_direito_= ((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21]));
					__NEW_[21]= _prvt_Ca_Trop_lado_direito_ * _prvt_dtime + __OLD_[21];
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
				    
// 				}
				/*if(_prvt_savingRate!=0){
					#pragma omp single
					{
						this->V_old_ = __OLD_[0];
						this->V_new_ = __NEW_[0];
						this->xr1_old_ = __OLD_[1];
						this->xr1_new_ = __NEW_[1];
						this->xr2_old_ = __OLD_[2];
						this->xr2_new_ = __NEW_[2];
						this->xs_old_ = __OLD_[3];
						this->xs_new_ = __NEW_[3];
						this->m_old_ = __OLD_[4];
						this->m_new_ = __NEW_[4];
						this->h_old_ = __OLD_[5];
						this->h_new_ = __NEW_[5];
						this->d_old_ = __OLD_[6];
						this->d_new_ = __NEW_[6];
						this->f_old_ = __OLD_[7];
						this->f_new_ = __NEW_[7];
						this->f2_old_ = __OLD_[8];
						this->f2_new_ = __NEW_[8];
						this->f2ds_old_ = __OLD_[9];
						this->f2ds_new_ = __NEW_[9];
						this->s_old_ = __OLD_[10];
						this->s_new_ = __NEW_[10];
						this->r_old_ = __OLD_[11];
						this->r_new_ = __NEW_[11];
						this->ActFrac_old_ = __OLD_[12];
						this->ActFrac_new_ = __NEW_[12];
						this->ProdFrac_old_ = __OLD_[13];
						this->ProdFrac_new_ = __NEW_[13];
						this->Na_i_old_ = __OLD_[14];
						this->Na_i_new_ = __NEW_[14];
						this->K_i_old_ = __OLD_[15];
						this->K_i_new_ = __NEW_[15];
						this->Ca_i_old_ = __OLD_[16];
						this->Ca_i_new_ = __NEW_[16];
						this->Ca_ds_old_ = __OLD_[17];
						this->Ca_ds_new_ = __NEW_[17];
						this->Ca_up_old_ = __OLD_[18];
						this->Ca_up_new_ = __NEW_[18];
						this->Ca_rel_old_ = __OLD_[19];
						this->Ca_rel_new_ = __NEW_[19];
						this->Ca_Calmod_old_ = __OLD_[20];
						this->Ca_Calmod_new_ = __NEW_[20];
						this->Ca_Trop_old_ = __OLD_[21];
						this->Ca_Trop_new_ = __NEW_[21];
						this->time_new = _prvt_time_new;
						save_step(fileptr, _EULER_);
					}
				}*/
				temp = __OLD_;
				__OLD_ = __NEW_;
				__NEW_= temp;
			}
// 			printf("%.1f\n", (float)tempoespera/(float)(nThreads));
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
// 		int tempoespera =0;
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
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -6.0000000000e+00,  _prvt_stim_start = 1.0000000000e-01,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 1.0000000000e+00,  _prvt_stim_duration = 1.5000000000e-03,  _prvt_Cm = 9.5000000000e-05,  _prvt_R = 8.3144720000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_F = 9.6485341500e+04,  _prvt_Na_o = 1.4000000000e+02,  _prvt_K_o = 4.0000000000e+00,  _prvt_P_kna = 3.0000000000e-02,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_g_K1 = 5.0000000000e-01,  _prvt_K_mk1 = 1.0000000000e+01,  _prvt_g_Kr1 = 2.1000000000e-03,  _prvt_g_Kr2 = 1.3000000000e-03,  _prvt_g_Ks = 2.6000000000e-03,  _prvt_g_Na = 2.5000000000e+00,  _prvt_delta_m = 1.0000000000e-05,  _prvt_shift_h = 0.0000000000e+00,  _prvt_g_pna = 4.0000000000e-03,  _prvt_g_bna = 6.0000000000e-04,  _prvt_FrICa = 1.0000000000e+00,  _prvt_P_Ca_L = 1.0000000000e-01,  _prvt_P_CaK = 2.0000000000e-03,  _prvt_P_CaNa = 1.0000000000e-02,  _prvt_speed_d = 3.0000000000e+00,  _prvt_delta_f = 1.0000000000e-04,  _prvt_speed_f = 3.0000000000e-01,  _prvt_Km_f2 = 1.0000000000e+05,  _prvt_R_decay = 2.0000000000e+01,  _prvt_Km_f2ds = 1.0000000000e-03,  _prvt_g_bca = 2.5000000000e-04,  _prvt_g_to = 5.0000000000e-03,  _prvt_g_tos = 0.0000000000e+00,  _prvt_i_NaK_max = 7.0000000000e-01,  _prvt_K_mK = 1.0000000000e+00,  _prvt_K_mNa = 4.0000000000e+01,  _prvt_FRiNaCa = 1.0000000000e-03,  _prvt_k_NaCa = 5.0000000000e-04,  _prvt_gamma = 5.0000000000e-01,  _prvt_n_NaCa = 3.0000000000e+00,  _prvt_d_NaCa = 0.0000000000e+00,  _prvt_K_cyca = 3.0000000000e-04,  _prvt_K_xcs = 4.0000000000e-01,  _prvt_K_srca = 5.0000000000e-01,  _prvt_alpha_up = 4.0000000000e-01,  _prvt_beta_up = 3.0000000000e-02,  _prvt_K_m_Ca_cyt = 5.0000000000e-04,  _prvt_K_m_Ca_ds = 1.0000000000e-02,  _prvt_K_m_rel = 2.5000000000e+02,  _prvt_K_leak_rate = 5.0000000000e-02,  _prvt_radius = 1.2000000000e-02,  _prvt_length = 7.4000000000e-02,  _prvt_V_e_ratio = 4.0000000000e-01,  _prvt_V_up_ratio = 1.0000000000e-02,  _prvt_V_rel_ratio = 1.0000000000e-01,  _prvt_V_ds_ratio = 1.0000000000e-01,  _prvt_Kdecay = 1.0000000000e+01,  _prvt_alpha_Calmod = 1.0000000000e+05,  _prvt_Calmod = 2.0000000000e-02,  _prvt_beta_Calmod = 5.0000000000e+01,  _prvt_alpha_Trop = 1.0000000000e+05,  _prvt_Trop = 5.0000000000e-02,  _prvt_beta_Trop = 2.0000000000e+02, 
			//private aux variables
			 _prvt_calc_i_Stim=0.0,  _prvt_calc_E_Na=0.0,  _prvt_calc_E_K=0.0,  _prvt_calc_E_Ks=0.0,  _prvt_calc_E_Ca=0.0,  _prvt_calc_E_mh=0.0,  _prvt_calc_i_K1=0.0,  _prvt_calc_i_Kr=0.0,  _prvt_calc_alpha_xr1=0.0,  _prvt_calc_beta_xr1=0.0,  _prvt_calc_alpha_xr2=0.0,  _prvt_calc_beta_xr2=0.0,  _prvt_calc_i_Ks=0.0,  _prvt_calc_alpha_xs=0.0,  _prvt_calc_beta_xs=0.0,  _prvt_calc_i_Na=0.0,  _prvt_calc_E0_m=0.0,  _prvt_calc_alpha_m=0.0,  _prvt_calc_beta_m=0.0,  _prvt_calc_alpha_h=0.0,  _prvt_calc_beta_h=0.0,  _prvt_calc_i_p_Na=0.0,  _prvt_calc_i_b_Na=0.0,  _prvt_calc_i_Ca_L_Ca_cyt=0.0,  _prvt_calc_i_Ca_L_K_cyt=0.0,  _prvt_calc_i_Ca_L_Na_cyt=0.0,  _prvt_calc_i_Ca_L_Ca_ds=0.0,  _prvt_calc_i_Ca_L_K_ds=0.0,  _prvt_calc_i_Ca_L_Na_ds=0.0,  _prvt_calc_i_Ca_L=0.0,  _prvt_calc_E0_d=0.0,  _prvt_calc_alpha_d=0.0,  _prvt_calc_beta_d=0.0,  _prvt_calc_E0_f=0.0,  _prvt_calc_alpha_f=0.0,  _prvt_calc_beta_f=0.0,  _prvt_calc_i_b_Ca=0.0,  _prvt_calc_i_to=0.0,  _prvt_calc_alpha_s=0.0,  _prvt_calc_beta_s=0.0,  _prvt_calc_i_NaK=0.0,  _prvt_calc_i_NaCa_cyt=0.0,  _prvt_calc_i_NaCa_ds=0.0,  _prvt_calc_i_NaCa=0.0,  _prvt_calc_K_1=0.0,  _prvt_calc_K_2=0.0,  _prvt_calc_i_up=0.0,  _prvt_calc_i_trans=0.0,  _prvt_calc_VoltDep=0.0,  _prvt_calc_CaiReg=0.0,  _prvt_calc_CadsReg=0.0,  _prvt_calc_RegBindSite=0.0,  _prvt_calc_ActRate=0.0,  _prvt_calc_InactRate=0.0,  _prvt_calc_SpeedRel=0.0,  _prvt_calc_PrecFrac=0.0,  _prvt_calc_i_rel=0.0,  _prvt_calc_V_Cell=0.0,  _prvt_calc_V_i_ratio=0.0,  _prvt_calc_V_i=0.0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_xr1_lado_direito_,  _prvt_xr2_lado_direito_,  _prvt_xs_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h_lado_direito_,  _prvt_d_lado_direito_,  _prvt_f_lado_direito_,  _prvt_f2_lado_direito_,  _prvt_f2ds_lado_direito_,  _prvt_s_lado_direito_,  _prvt_r_lado_direito_,  _prvt_ActFrac_lado_direito_,  _prvt_ProdFrac_lado_direito_,  _prvt_Na_i_lado_direito_,  _prvt_K_i_lado_direito_,  _prvt_Ca_i_lado_direito_,  _prvt_Ca_ds_lado_direito_,  _prvt_Ca_up_lado_direito_,  _prvt_Ca_rel_lado_direito_,  _prvt_Ca_Calmod_lado_direito_,  _prvt_Ca_Trop_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;
			__NEW_[0] = __OLD_[0] = -9.2849333000e+01;
			__NEW_[1] = __OLD_[1] = 1.0300000000e-05;
			__NEW_[2] = __OLD_[2] = 2.0000000000e-07;
			__NEW_[3] = __OLD_[3] = 1.3020000000e-03;
			__NEW_[4] = __OLD_[4] = 1.6203000000e-03;
			__NEW_[5] = __OLD_[5] = 9.9440360000e-01;
			__NEW_[6] = __OLD_[6] = 0.0000000000e+00;
			__NEW_[7] = __OLD_[7] = 1.0000000000e+00;
			__NEW_[8] = __OLD_[8] = 9.3491970000e-01;
			__NEW_[9] = __OLD_[9] = 9.6519580000e-01;
			__NEW_[10] = __OLD_[10] = 9.9486450000e-01;
			__NEW_[11] = __OLD_[11] = 0.0000000000e+00;
			__NEW_[12] = __OLD_[12] = 4.2614000000e-03;
			__NEW_[13] = __OLD_[13] = 4.0681540000e-01;
			__NEW_[14] = __OLD_[14] = 7.3321223000e+00;
			__NEW_[15] = __OLD_[15] = 1.3656442810e+02;
			__NEW_[16] = __OLD_[16] = 1.4000000000e-05;
			__NEW_[17] = __OLD_[17] = 1.8800000000e-05;
			__NEW_[18] = __OLD_[18] = 4.5318890000e-01;
			__NEW_[19] = __OLD_[19] = 4.4819270000e-01;
			__NEW_[20] = __OLD_[20] = 5.5550000000e-04;
			__NEW_[21] = __OLD_[21] = 3.5420000000e-04;
			_prvt_time_new += _prvt_dtime;
			if(omp_get_thread_num()==_prvt_tree_thread[0])
			{
				_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
				_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[14])));
				_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[15])));
				_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[15]+(_prvt_P_kna*__OLD_[14])))));
				_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[16])));
				_prvt_calc_E_mh = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_Na_o+(1.2000000000e-01*_prvt_K_o))/(__OLD_[14]+(1.2000000000e-01*__OLD_[15])))));
				_prvt_calc_i_Ca_L_Ca_cyt = (((((1.0000000000e+00-_prvt_FrICa)*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_K_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_Na_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_Ca_ds = ((((_prvt_FrICa*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_K_ds = ((((_prvt_FrICa*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_Na_ds = ((((_prvt_FrICa*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_K_o)/(_prvt_K_mK+_prvt_K_o))*__OLD_[14])/(_prvt_K_mNa+__OLD_[14]));
				_prvt_calc_i_NaCa_cyt = (((1.0000000000e+00-_prvt_FRiNaCa)*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[16])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[16]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[16]/6.9000000000e-03))));
				_prvt_calc_i_NaCa_ds = ((_prvt_FRiNaCa*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[17])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[17]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[17]/6.9000000000e-03))));
				_prvt_calc_K_1 = ((_prvt_K_cyca*_prvt_K_xcs)/_prvt_K_srca);
				_prvt_calc_i_trans = (5.0000000000e+01*(__OLD_[18]-__OLD_[19]));
				_prvt_calc_i_rel = (((pow((__OLD_[12]/(__OLD_[12]+2.5000000000e-01)),2.0000000000e+00)*_prvt_K_m_rel)+_prvt_K_leak_rate)*__OLD_[19]);
				_prvt_calc_V_Cell = (3.1415926540e+00*pow(_prvt_radius,2.0000000000e+00)*_prvt_length);
				_prvt_calc_V_i_ratio = (((1.0000000000e+00-_prvt_V_e_ratio)-_prvt_V_up_ratio)-_prvt_V_rel_ratio);
				_prvt_calc_i_K1 = ((((_prvt_g_K1*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk1))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01)*_prvt_F*1.2500000000e+00)/(_prvt_R*_prvt_T)))));
				_prvt_calc_i_Kr = (((((_prvt_g_Kr1*__OLD_[1])+(_prvt_g_Kr2*__OLD_[2]))*1.0000000000e+00)/(1.0000000000e+00+exp(((__OLD_[0]+9.0000000000e+00)/2.2400000000e+01))))*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
				_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*(__OLD_[0]-_prvt_calc_E_mh));
				_prvt_calc_i_p_Na = (((_prvt_g_pna*1.0000000000e+00)/(1.0000000000e+00+exp(((-(__OLD_[0]+5.2000000000e+01))/8.0000000000e+00))))*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_Ca_L = (_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_Ca_L_Na_ds);
				_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
				_prvt_calc_i_to = (_prvt_g_to*(_prvt_g_tos+(__OLD_[10]*(1.0000000000e+00-_prvt_g_tos)))*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_NaCa = (_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds);
				_prvt_calc_K_2 = (__OLD_[16]+(__OLD_[18]*_prvt_calc_K_1)+(_prvt_K_cyca*_prvt_K_xcs)+_prvt_K_cyca);
				_prvt_calc_V_i = (_prvt_calc_V_Cell*_prvt_calc_V_i_ratio);
				_prvt_calc_i_up = (((__OLD_[16]/_prvt_calc_K_2)*_prvt_alpha_up)-(((__OLD_[18]*_prvt_calc_K_1)/_prvt_calc_K_2)*_prvt_beta_up));
				__K1_[0]= (((-1.0000000000e+00)/_prvt_Cm)*(_prvt_calc_i_Stim+_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_p_Na+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds+_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds+_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_b_Ca));
				__NEW_[0]= __K1_[0] * _prvt_dtime + __OLD_[0];
				__K1_[14]= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(_prvt_calc_i_Na+_prvt_calc_i_p_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa_cyt)+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds));
				__NEW_[14]= __K1_[14] * _prvt_dtime + __OLD_[14];
				__K1_[15]= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*((_prvt_calc_i_K1+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_to)-(2.0000000000e+00*_prvt_calc_i_NaK)));
				__NEW_[15]= __K1_[15] * _prvt_dtime + __OLD_[15];
				__K1_[16]= (((((((-1.0000000000e+00)/(2.0000000000e+00*1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(((_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_b_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa_cyt))-(2.0000000000e+00*_prvt_calc_i_NaCa_ds)))+(__OLD_[17]*_prvt_V_ds_ratio*_prvt_Kdecay)+((_prvt_calc_i_rel*_prvt_V_rel_ratio)/_prvt_calc_V_i_ratio))-((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20])))-((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21])))-_prvt_calc_i_up);
				__NEW_[16]= __K1_[16] * _prvt_dtime + __OLD_[16];
				__K1_[17]= ((((-1.0000000000e+00)*_prvt_calc_i_Ca_L_Ca_ds)/(2.0000000000e+00*1.0000000000e+00*_prvt_V_ds_ratio*_prvt_calc_V_i*_prvt_F))-(__OLD_[17]*_prvt_Kdecay));
				__NEW_[17]= __K1_[17] * _prvt_dtime + __OLD_[17];
				__K1_[18]= (((_prvt_calc_V_i_ratio/_prvt_V_up_ratio)*_prvt_calc_i_up)-_prvt_calc_i_trans);
				__NEW_[18]= __K1_[18] * _prvt_dtime + __OLD_[18];
				__K1_[19]= (((_prvt_V_up_ratio/_prvt_V_rel_ratio)*_prvt_calc_i_trans)-_prvt_calc_i_rel);
				__NEW_[19]= __K1_[19] * _prvt_dtime + __OLD_[19];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[1])
			{
				_prvt_calc_alpha_xr1 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
				_prvt_calc_beta_xr1 = (5.0000000000e-02*exp(((-(__OLD_[0]-2.0000000000e+01))/1.5000000000e+01)));
				__K1_[1]= ((_prvt_calc_alpha_xr1*(1.0000000000e+00-__OLD_[1]))-(_prvt_calc_beta_xr1*__OLD_[1]));
				__NEW_[1]= __K1_[1] * _prvt_dtime + __OLD_[1];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[2])
			{
				_prvt_calc_alpha_xr2 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
				_prvt_calc_beta_xr2 = (4.0000000000e-01*exp((-pow(((__OLD_[0]+3.0000000000e+01)/3.0000000000e+01),3.0000000000e+00))));
				__K1_[2]= ((_prvt_calc_alpha_xr2*(1.0000000000e+00-__OLD_[2]))-(_prvt_calc_beta_xr2*__OLD_[2]));
				__NEW_[2]= __K1_[2] * _prvt_dtime + __OLD_[2];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[3])
			{
				_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
				_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
				__K1_[3]= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[3]))-(_prvt_calc_beta_xs*__OLD_[3]));
				__NEW_[3]= __K1_[3] * _prvt_dtime + __OLD_[3];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[4])
			{
				_prvt_calc_E0_m = (__OLD_[0]+4.1000000000e+01);
				_prvt_calc_beta_m = (8.0000000000e+03*exp(((-5.6000000000e-02)*(__OLD_[0]+6.6000000000e+01))));
				_prvt_calc_alpha_m = ((fabs(_prvt_calc_E0_m)<_prvt_delta_m))
?(2.0000000000e+03)
:(((2.0000000000e+02*_prvt_calc_E0_m)/(1.0000000000e+00-exp(((-1.0000000000e-01)*_prvt_calc_E0_m)))));
				__K1_[4]= ((_prvt_calc_alpha_m*(1.0000000000e+00-__OLD_[4]))-(_prvt_calc_beta_m*__OLD_[4]));
				__NEW_[4]= __K1_[4] * _prvt_dtime + __OLD_[4];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[5])
			{
				_prvt_calc_alpha_h = (2.0000000000e+01*exp(((-1.2500000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))));
				_prvt_calc_beta_h = (2.0000000000e+03/(1.0000000000e+00+(3.2000000000e+02*exp(((-1.0000000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))))));
				__K1_[5]= ((_prvt_calc_alpha_h*(1.0000000000e+00-__OLD_[5]))-(_prvt_calc_beta_h*__OLD_[5]));
				__NEW_[5]= __K1_[5] * _prvt_dtime + __OLD_[5];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[6])
			{
				_prvt_calc_E0_d = ((__OLD_[0]+2.4000000000e+01)-5.0000000000e+00);
				_prvt_calc_alpha_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((3.0000000000e+01*_prvt_calc_E0_d)/(1.0000000000e+00-exp(((-_prvt_calc_E0_d)/4.0000000000e+00)))));
				_prvt_calc_beta_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((1.2000000000e+01*_prvt_calc_E0_d)/(exp((_prvt_calc_E0_d/1.0000000000e+01))-1.0000000000e+00)));
				__K1_[6]= (_prvt_speed_d*((_prvt_calc_alpha_d*(1.0000000000e+00-__OLD_[6]))-(_prvt_calc_beta_d*__OLD_[6])));
				__NEW_[6]= __K1_[6] * _prvt_dtime + __OLD_[6];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[7])
			{
				_prvt_calc_E0_f = (__OLD_[0]+3.4000000000e+01);
				_prvt_calc_beta_f = (1.2000000000e+01/(1.0000000000e+00+exp((((-1.0000000000e+00)*(__OLD_[0]+3.4000000000e+01))/4.0000000000e+00))));
				_prvt_calc_alpha_f = ((fabs(_prvt_calc_E0_f)<_prvt_delta_f))
?(2.5000000000e+01)
:(((6.2500000000e+00*_prvt_calc_E0_f)/(exp((_prvt_calc_E0_f/4.0000000000e+00))-1.0000000000e+00)));
				__K1_[7]= (_prvt_speed_f*((_prvt_calc_alpha_f*(1.0000000000e+00-__OLD_[7]))-(_prvt_calc_beta_f*__OLD_[7])));
				__NEW_[7]= __K1_[7] * _prvt_dtime + __OLD_[7];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[8])
			{
				_prvt_calc_alpha_s = (3.3000000000e-02*exp(((-__OLD_[0])/1.7000000000e+01)));
				_prvt_calc_beta_s = (3.3000000000e+01/(1.0000000000e+00+exp(((-1.2500000000e-01)*(__OLD_[0]+1.0000000000e+01)))));
				__K1_[10]= ((_prvt_calc_alpha_s*(1.0000000000e+00-__OLD_[10]))-(_prvt_calc_beta_s*__OLD_[10]));
				__NEW_[10]= __K1_[10] * _prvt_dtime + __OLD_[10];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[9])
			{
				_prvt_calc_VoltDep = exp((8.0000000000e-02*(__OLD_[0]-4.0000000000e+01)));
				_prvt_calc_CaiReg = (__OLD_[16]/(__OLD_[16]+_prvt_K_m_Ca_cyt));
				_prvt_calc_CadsReg = (__OLD_[17]/(__OLD_[17]+_prvt_K_m_Ca_ds));
				_prvt_calc_SpeedRel = ((__OLD_[0]<(-5.0000000000e+01)))
?(5.0000000000e+00)
:(1.0000000000e+00);
				_prvt_calc_PrecFrac = ((1.0000000000e+00-__OLD_[12])-__OLD_[13]);
				_prvt_calc_RegBindSite = (_prvt_calc_CaiReg+((1.0000000000e+00-_prvt_calc_CaiReg)*_prvt_calc_CadsReg));
				_prvt_calc_ActRate = ((0.0000000000e+00*_prvt_calc_VoltDep)+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
				_prvt_calc_InactRate = (6.0000000000e+01+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
				__K1_[12]= ((_prvt_calc_PrecFrac*_prvt_calc_SpeedRel*_prvt_calc_ActRate)-(__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate));
				__NEW_[12]= __K1_[12] * _prvt_dtime + __OLD_[12];
				__K1_[13]= ((__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate)-(_prvt_calc_SpeedRel*1.0000000000e+00*__OLD_[13]));
				__NEW_[13]= __K1_[13] * _prvt_dtime + __OLD_[13];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[10])
			{
				__K1_[8]= (1.0000000000e+00-(1.0000000000e+00*((__OLD_[16]/(_prvt_Km_f2+__OLD_[16]))+__OLD_[8])));
				__NEW_[8]= __K1_[8] * _prvt_dtime + __OLD_[8];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[11])
			{
				__K1_[9]= (_prvt_R_decay*(1.0000000000e+00-((__OLD_[17]/(_prvt_Km_f2ds+__OLD_[17]))+__OLD_[9])));
				__NEW_[9]= __K1_[9] * _prvt_dtime + __OLD_[9];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[12])
			{
				__K1_[11]= (3.3300000000e+02*((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+4.0000000000e+00))/5.0000000000e+00))))-__OLD_[11]));
				__NEW_[11]= __K1_[11] * _prvt_dtime + __OLD_[11];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[13])
			{
				__K1_[20]= ((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20]));
				__NEW_[20]= __K1_[20] * _prvt_dtime + __OLD_[20];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[14])
			{
				__K1_[21]= ((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21]));
				__NEW_[21]= __K1_[21] * _prvt_dtime + __OLD_[21];
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
// // 			gettimeofday(&t2, NULL);
// 			#pragma omp critical
// 			{
// 			    int aux = (int)((int)t2.tv_usec - (int)t1.tv_usec);
// // 			    tempoespera += (aux<0)?0:aux;
// 			}
			
			while(_prvt_time_new<=_prvt_finalTime){
				_prvt_time_new += _prvt_dtime;
				if(omp_get_thread_num()==_prvt_tree_thread[0])
				{
					_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[14])));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[15])));
					_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[15]+(_prvt_P_kna*__OLD_[14])))));
					_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[16])));
					_prvt_calc_E_mh = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_Na_o+(1.2000000000e-01*_prvt_K_o))/(__OLD_[14]+(1.2000000000e-01*__OLD_[15])))));
					_prvt_calc_i_Ca_L_Ca_cyt = (((((1.0000000000e+00-_prvt_FrICa)*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_K_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Na_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Ca_ds = ((((_prvt_FrICa*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_K_ds = ((((_prvt_FrICa*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Na_ds = ((((_prvt_FrICa*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_K_o)/(_prvt_K_mK+_prvt_K_o))*__OLD_[14])/(_prvt_K_mNa+__OLD_[14]));
					_prvt_calc_i_NaCa_cyt = (((1.0000000000e+00-_prvt_FRiNaCa)*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[16])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[16]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[16]/6.9000000000e-03))));
					_prvt_calc_i_NaCa_ds = ((_prvt_FRiNaCa*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[17])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[17]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[17]/6.9000000000e-03))));
					_prvt_calc_K_1 = ((_prvt_K_cyca*_prvt_K_xcs)/_prvt_K_srca);
					_prvt_calc_i_trans = (5.0000000000e+01*(__OLD_[18]-__OLD_[19]));
					_prvt_calc_i_rel = (((pow((__OLD_[12]/(__OLD_[12]+2.5000000000e-01)),2.0000000000e+00)*_prvt_K_m_rel)+_prvt_K_leak_rate)*__OLD_[19]);
					_prvt_calc_V_Cell = (3.1415926540e+00*pow(_prvt_radius,2.0000000000e+00)*_prvt_length);
					_prvt_calc_V_i_ratio = (((1.0000000000e+00-_prvt_V_e_ratio)-_prvt_V_up_ratio)-_prvt_V_rel_ratio);
					_prvt_calc_i_K1 = ((((_prvt_g_K1*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk1))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01)*_prvt_F*1.2500000000e+00)/(_prvt_R*_prvt_T)))));
					_prvt_calc_i_Kr = (((((_prvt_g_Kr1*__OLD_[1])+(_prvt_g_Kr2*__OLD_[2]))*1.0000000000e+00)/(1.0000000000e+00+exp(((__OLD_[0]+9.0000000000e+00)/2.2400000000e+01))))*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
					_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*(__OLD_[0]-_prvt_calc_E_mh));
					_prvt_calc_i_p_Na = (((_prvt_g_pna*1.0000000000e+00)/(1.0000000000e+00+exp(((-(__OLD_[0]+5.2000000000e+01))/8.0000000000e+00))))*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Ca_L = (_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_Ca_L_Na_ds);
					_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_to = (_prvt_g_to*(_prvt_g_tos+(__OLD_[10]*(1.0000000000e+00-_prvt_g_tos)))*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_NaCa = (_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds);
					_prvt_calc_K_2 = (__OLD_[16]+(__OLD_[18]*_prvt_calc_K_1)+(_prvt_K_cyca*_prvt_K_xcs)+_prvt_K_cyca);
					_prvt_calc_V_i = (_prvt_calc_V_Cell*_prvt_calc_V_i_ratio);
					_prvt_calc_i_up = (((__OLD_[16]/_prvt_calc_K_2)*_prvt_alpha_up)-(((__OLD_[18]*_prvt_calc_K_1)/_prvt_calc_K_2)*_prvt_beta_up));
					__K2_[0]= (((-1.0000000000e+00)/_prvt_Cm)*(_prvt_calc_i_Stim+_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_p_Na+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds+_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds+_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_b_Ca));
					_prvt_aux_tol = fabs(__OLD_[0])*_prvt_rel_tol_;
					__TOL_[0] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[0] = fabs((_prvt_dtime/2) * (__K1_[0] - __K2_[0])/__TOL_[0]);
					__K2_[14]= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(_prvt_calc_i_Na+_prvt_calc_i_p_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa_cyt)+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds));
					_prvt_aux_tol = fabs(__OLD_[14])*_prvt_rel_tol_;
					__TOL_[14] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[14] = fabs((_prvt_dtime/2) * (__K1_[14] - __K2_[14])/__TOL_[14]);
					__K2_[15]= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*((_prvt_calc_i_K1+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_to)-(2.0000000000e+00*_prvt_calc_i_NaK)));
					_prvt_aux_tol = fabs(__OLD_[15])*_prvt_rel_tol_;
					__TOL_[15] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[15] = fabs((_prvt_dtime/2) * (__K1_[15] - __K2_[15])/__TOL_[15]);
					__K2_[16]= (((((((-1.0000000000e+00)/(2.0000000000e+00*1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(((_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_b_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa_cyt))-(2.0000000000e+00*_prvt_calc_i_NaCa_ds)))+(__OLD_[17]*_prvt_V_ds_ratio*_prvt_Kdecay)+((_prvt_calc_i_rel*_prvt_V_rel_ratio)/_prvt_calc_V_i_ratio))-((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20])))-((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21])))-_prvt_calc_i_up);
					_prvt_aux_tol = fabs(__OLD_[16])*_prvt_rel_tol_;
					__TOL_[16] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[16] = fabs((_prvt_dtime/2) * (__K1_[16] - __K2_[16])/__TOL_[16]);
					__K2_[17]= ((((-1.0000000000e+00)*_prvt_calc_i_Ca_L_Ca_ds)/(2.0000000000e+00*1.0000000000e+00*_prvt_V_ds_ratio*_prvt_calc_V_i*_prvt_F))-(__OLD_[17]*_prvt_Kdecay));
					_prvt_aux_tol = fabs(__OLD_[17])*_prvt_rel_tol_;
					__TOL_[17] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[17] = fabs((_prvt_dtime/2) * (__K1_[17] - __K2_[17])/__TOL_[17]);
					__K2_[18]= (((_prvt_calc_V_i_ratio/_prvt_V_up_ratio)*_prvt_calc_i_up)-_prvt_calc_i_trans);
					_prvt_aux_tol = fabs(__OLD_[18])*_prvt_rel_tol_;
					__TOL_[18] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[18] = fabs((_prvt_dtime/2) * (__K1_[18] - __K2_[18])/__TOL_[18]);
					__K2_[19]= (((_prvt_V_up_ratio/_prvt_V_rel_ratio)*_prvt_calc_i_trans)-_prvt_calc_i_rel);
					_prvt_aux_tol = fabs(__OLD_[19])*_prvt_rel_tol_;
					__TOL_[19] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[19] = fabs((_prvt_dtime/2) * (__K1_[19] - __K2_[19])/__TOL_[19]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[1])
				{
					_prvt_calc_alpha_xr1 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
					_prvt_calc_beta_xr1 = (5.0000000000e-02*exp(((-(__OLD_[0]-2.0000000000e+01))/1.5000000000e+01)));
					__K2_[1]= ((_prvt_calc_alpha_xr1*(1.0000000000e+00-__OLD_[1]))-(_prvt_calc_beta_xr1*__OLD_[1]));
					_prvt_aux_tol = fabs(__OLD_[1])*_prvt_rel_tol_;
					__TOL_[1] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[1] = fabs((_prvt_dtime/2) * (__K1_[1] - __K2_[1])/__TOL_[1]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[2])
				{
					_prvt_calc_alpha_xr2 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
					_prvt_calc_beta_xr2 = (4.0000000000e-01*exp((-pow(((__OLD_[0]+3.0000000000e+01)/3.0000000000e+01),3.0000000000e+00))));
					__K2_[2]= ((_prvt_calc_alpha_xr2*(1.0000000000e+00-__OLD_[2]))-(_prvt_calc_beta_xr2*__OLD_[2]));
					_prvt_aux_tol = fabs(__OLD_[2])*_prvt_rel_tol_;
					__TOL_[2] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[2] = fabs((_prvt_dtime/2) * (__K1_[2] - __K2_[2])/__TOL_[2]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[3])
				{
					_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
					_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
					__K2_[3]= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[3]))-(_prvt_calc_beta_xs*__OLD_[3]));
					_prvt_aux_tol = fabs(__OLD_[3])*_prvt_rel_tol_;
					__TOL_[3] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[3] = fabs((_prvt_dtime/2) * (__K1_[3] - __K2_[3])/__TOL_[3]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[4])
				{
					_prvt_calc_E0_m = (__OLD_[0]+4.1000000000e+01);
					_prvt_calc_beta_m = (8.0000000000e+03*exp(((-5.6000000000e-02)*(__OLD_[0]+6.6000000000e+01))));
					_prvt_calc_alpha_m = ((fabs(_prvt_calc_E0_m)<_prvt_delta_m))
?(2.0000000000e+03)
:(((2.0000000000e+02*_prvt_calc_E0_m)/(1.0000000000e+00-exp(((-1.0000000000e-01)*_prvt_calc_E0_m)))));
					__K2_[4]= ((_prvt_calc_alpha_m*(1.0000000000e+00-__OLD_[4]))-(_prvt_calc_beta_m*__OLD_[4]));
					_prvt_aux_tol = fabs(__OLD_[4])*_prvt_rel_tol_;
					__TOL_[4] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[4] = fabs((_prvt_dtime/2) * (__K1_[4] - __K2_[4])/__TOL_[4]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[5])
				{
					_prvt_calc_alpha_h = (2.0000000000e+01*exp(((-1.2500000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))));
					_prvt_calc_beta_h = (2.0000000000e+03/(1.0000000000e+00+(3.2000000000e+02*exp(((-1.0000000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))))));
					__K2_[5]= ((_prvt_calc_alpha_h*(1.0000000000e+00-__OLD_[5]))-(_prvt_calc_beta_h*__OLD_[5]));
					_prvt_aux_tol = fabs(__OLD_[5])*_prvt_rel_tol_;
					__TOL_[5] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[5] = fabs((_prvt_dtime/2) * (__K1_[5] - __K2_[5])/__TOL_[5]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[6])
				{
					_prvt_calc_E0_d = ((__OLD_[0]+2.4000000000e+01)-5.0000000000e+00);
					_prvt_calc_alpha_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((3.0000000000e+01*_prvt_calc_E0_d)/(1.0000000000e+00-exp(((-_prvt_calc_E0_d)/4.0000000000e+00)))));
					_prvt_calc_beta_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((1.2000000000e+01*_prvt_calc_E0_d)/(exp((_prvt_calc_E0_d/1.0000000000e+01))-1.0000000000e+00)));
					__K2_[6]= (_prvt_speed_d*((_prvt_calc_alpha_d*(1.0000000000e+00-__OLD_[6]))-(_prvt_calc_beta_d*__OLD_[6])));
					_prvt_aux_tol = fabs(__OLD_[6])*_prvt_rel_tol_;
					__TOL_[6] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[6] = fabs((_prvt_dtime/2) * (__K1_[6] - __K2_[6])/__TOL_[6]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[7])
				{
					_prvt_calc_E0_f = (__OLD_[0]+3.4000000000e+01);
					_prvt_calc_beta_f = (1.2000000000e+01/(1.0000000000e+00+exp((((-1.0000000000e+00)*(__OLD_[0]+3.4000000000e+01))/4.0000000000e+00))));
					_prvt_calc_alpha_f = ((fabs(_prvt_calc_E0_f)<_prvt_delta_f))
?(2.5000000000e+01)
:(((6.2500000000e+00*_prvt_calc_E0_f)/(exp((_prvt_calc_E0_f/4.0000000000e+00))-1.0000000000e+00)));
					__K2_[7]= (_prvt_speed_f*((_prvt_calc_alpha_f*(1.0000000000e+00-__OLD_[7]))-(_prvt_calc_beta_f*__OLD_[7])));
					_prvt_aux_tol = fabs(__OLD_[7])*_prvt_rel_tol_;
					__TOL_[7] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[7] = fabs((_prvt_dtime/2) * (__K1_[7] - __K2_[7])/__TOL_[7]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[8])
				{
					_prvt_calc_alpha_s = (3.3000000000e-02*exp(((-__OLD_[0])/1.7000000000e+01)));
					_prvt_calc_beta_s = (3.3000000000e+01/(1.0000000000e+00+exp(((-1.2500000000e-01)*(__OLD_[0]+1.0000000000e+01)))));
					__K2_[10]= ((_prvt_calc_alpha_s*(1.0000000000e+00-__OLD_[10]))-(_prvt_calc_beta_s*__OLD_[10]));
					_prvt_aux_tol = fabs(__OLD_[10])*_prvt_rel_tol_;
					__TOL_[10] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[10] = fabs((_prvt_dtime/2) * (__K1_[10] - __K2_[10])/__TOL_[10]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[9])
				{
					_prvt_calc_VoltDep = exp((8.0000000000e-02*(__OLD_[0]-4.0000000000e+01)));
					_prvt_calc_CaiReg = (__OLD_[16]/(__OLD_[16]+_prvt_K_m_Ca_cyt));
					_prvt_calc_CadsReg = (__OLD_[17]/(__OLD_[17]+_prvt_K_m_Ca_ds));
					_prvt_calc_SpeedRel = ((__OLD_[0]<(-5.0000000000e+01)))
?(5.0000000000e+00)
:(1.0000000000e+00);
					_prvt_calc_PrecFrac = ((1.0000000000e+00-__OLD_[12])-__OLD_[13]);
					_prvt_calc_RegBindSite = (_prvt_calc_CaiReg+((1.0000000000e+00-_prvt_calc_CaiReg)*_prvt_calc_CadsReg));
					_prvt_calc_ActRate = ((0.0000000000e+00*_prvt_calc_VoltDep)+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
					_prvt_calc_InactRate = (6.0000000000e+01+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
					__K2_[12]= ((_prvt_calc_PrecFrac*_prvt_calc_SpeedRel*_prvt_calc_ActRate)-(__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate));
					_prvt_aux_tol = fabs(__OLD_[12])*_prvt_rel_tol_;
					__TOL_[12] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[12] = fabs((_prvt_dtime/2) * (__K1_[12] - __K2_[12])/__TOL_[12]);
					__K2_[13]= ((__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate)-(_prvt_calc_SpeedRel*1.0000000000e+00*__OLD_[13]));
					_prvt_aux_tol = fabs(__OLD_[13])*_prvt_rel_tol_;
					__TOL_[13] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[13] = fabs((_prvt_dtime/2) * (__K1_[13] - __K2_[13])/__TOL_[13]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[10])
				{
					__K2_[8]= (1.0000000000e+00-(1.0000000000e+00*((__OLD_[16]/(_prvt_Km_f2+__OLD_[16]))+__OLD_[8])));
					_prvt_aux_tol = fabs(__OLD_[8])*_prvt_rel_tol_;
					__TOL_[8] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[8] = fabs((_prvt_dtime/2) * (__K1_[8] - __K2_[8])/__TOL_[8]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[11])
				{
					__K2_[9]= (_prvt_R_decay*(1.0000000000e+00-((__OLD_[17]/(_prvt_Km_f2ds+__OLD_[17]))+__OLD_[9])));
					_prvt_aux_tol = fabs(__OLD_[9])*_prvt_rel_tol_;
					__TOL_[9] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[9] = fabs((_prvt_dtime/2) * (__K1_[9] - __K2_[9])/__TOL_[9]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[12])
				{
					__K2_[11]= (3.3300000000e+02*((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+4.0000000000e+00))/5.0000000000e+00))))-__OLD_[11]));
					_prvt_aux_tol = fabs(__OLD_[11])*_prvt_rel_tol_;
					__TOL_[11] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[11] = fabs((_prvt_dtime/2) * (__K1_[11] - __K2_[11])/__TOL_[11]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[13])
				{
					__K2_[20]= ((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20]));
					_prvt_aux_tol = fabs(__OLD_[20])*_prvt_rel_tol_;
					__TOL_[20] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[20] = fabs((_prvt_dtime/2) * (__K1_[20] - __K2_[20])/__TOL_[20]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[14])
				{
					__K2_[21]= ((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21]));
					_prvt_aux_tol = fabs(__OLD_[21])*_prvt_rel_tol_;
					__TOL_[21] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[21] = fabs((_prvt_dtime/2) * (__K1_[21] - __K2_[21])/__TOL_[21]);
				}
				_prvt_time_new -= _prvt_dtime;
// 				timeval t3, t4;
// 				gettimeofday(&t3, NULL);
				#pragma omp barrier
// 				gettimeofday(&t4, NULL);
// 				#pragma omp critical
// 				{
// 				    int aux = (int)((int)t4.tv_usec - (int)t3.tv_usec);
// 				    tempoespera += (aux<0)?0:aux;
				    
// 				}
				
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
						__NEW_[14] = __K1_[14] * _prvt_dtime + __OLD_AUX_[14];
						__NEW_[15] = __K1_[15] * _prvt_dtime + __OLD_AUX_[15];
						__NEW_[16] = __K1_[16] * _prvt_dtime + __OLD_AUX_[16];
						__NEW_[17] = __K1_[17] * _prvt_dtime + __OLD_AUX_[17];
						__NEW_[18] = __K1_[18] * _prvt_dtime + __OLD_AUX_[18];
						__NEW_[19] = __K1_[19] * _prvt_dtime + __OLD_AUX_[19];
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
						__NEW_[10] = __K1_[10] * _prvt_dtime + __OLD_AUX_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[12] = __K1_[12] * _prvt_dtime + __OLD_AUX_[12];
						__NEW_[13] = __K1_[13] * _prvt_dtime + __OLD_AUX_[13];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[8] = __K1_[8] * _prvt_dtime + __OLD_AUX_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[9] = __K1_[9] * _prvt_dtime + __OLD_AUX_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[11] = __K1_[11] * _prvt_dtime + __OLD_AUX_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[13])
					{
						__NEW_[20] = __K1_[20] * _prvt_dtime + __OLD_AUX_[20];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[14])
					{
						__NEW_[21] = __K1_[21] * _prvt_dtime + __OLD_AUX_[21];
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
							this->xr1_old_ = __OLD_AUX_[1];
							this->xr1_new_ = __OLD_[1];
							this->xr2_old_ = __OLD_AUX_[2];
							this->xr2_new_ = __OLD_[2];
							this->xs_old_ = __OLD_AUX_[3];
							this->xs_new_ = __OLD_[3];
							this->m_old_ = __OLD_AUX_[4];
							this->m_new_ = __OLD_[4];
							this->h_old_ = __OLD_AUX_[5];
							this->h_new_ = __OLD_[5];
							this->d_old_ = __OLD_AUX_[6];
							this->d_new_ = __OLD_[6];
							this->f_old_ = __OLD_AUX_[7];
							this->f_new_ = __OLD_[7];
							this->f2_old_ = __OLD_AUX_[8];
							this->f2_new_ = __OLD_[8];
							this->f2ds_old_ = __OLD_AUX_[9];
							this->f2ds_new_ = __OLD_[9];
							this->s_old_ = __OLD_AUX_[10];
							this->s_new_ = __OLD_[10];
							this->r_old_ = __OLD_AUX_[11];
							this->r_new_ = __OLD_[11];
							this->ActFrac_old_ = __OLD_AUX_[12];
							this->ActFrac_new_ = __OLD_[12];
							this->ProdFrac_old_ = __OLD_AUX_[13];
							this->ProdFrac_new_ = __OLD_[13];
							this->Na_i_old_ = __OLD_AUX_[14];
							this->Na_i_new_ = __OLD_[14];
							this->K_i_old_ = __OLD_AUX_[15];
							this->K_i_new_ = __OLD_[15];
							this->Ca_i_old_ = __OLD_AUX_[16];
							this->Ca_i_new_ = __OLD_[16];
							this->Ca_ds_old_ = __OLD_AUX_[17];
							this->Ca_ds_new_ = __OLD_[17];
							this->Ca_up_old_ = __OLD_AUX_[18];
							this->Ca_up_new_ = __OLD_[18];
							this->Ca_rel_old_ = __OLD_AUX_[19];
							this->Ca_rel_new_ = __OLD_[19];
							this->Ca_Calmod_old_ = __OLD_AUX_[20];
							this->Ca_Calmod_new_ = __OLD_[20];
							this->Ca_Trop_old_ = __OLD_AUX_[21];
							this->Ca_Trop_new_ = __OLD_[21];
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
						__NEW_[14] = __K2_[14] * _prvt_dtime + __OLD_[14];
						__NEW_[15] = __K2_[15] * _prvt_dtime + __OLD_[15];
						__NEW_[16] = __K2_[16] * _prvt_dtime + __OLD_[16];
						__NEW_[17] = __K2_[17] * _prvt_dtime + __OLD_[17];
						__NEW_[18] = __K2_[18] * _prvt_dtime + __OLD_[18];
						__NEW_[19] = __K2_[19] * _prvt_dtime + __OLD_[19];
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
						__NEW_[10] = __K2_[10] * _prvt_dtime + __OLD_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[12] = __K2_[12] * _prvt_dtime + __OLD_[12];
						__NEW_[13] = __K2_[13] * _prvt_dtime + __OLD_[13];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[8] = __K2_[8] * _prvt_dtime + __OLD_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[9] = __K2_[9] * _prvt_dtime + __OLD_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[11] = __K2_[11] * _prvt_dtime + __OLD_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[13])
					{
						__NEW_[20] = __K2_[20] * _prvt_dtime + __OLD_[20];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[14])
					{
						__NEW_[21] = __K2_[21] * _prvt_dtime + __OLD_[21];
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
// 				timeval t5, t6;
// 				gettimeofday(&t5, NULL);
				#pragma omp barrier
// 				gettimeofday(&t6, NULL);
// 				#pragma omp critical
// 				{
// 				    int aux = (int)((int)t6.tv_usec - (int)t5.tv_usec);
// 				    tempoespera += (aux<0)?0:aux;
// 				}
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
			double  _prvt_time = 0.0000000000e+00,  _prvt_stim_amplitude = -6.0000000000e+00,  _prvt_stim_start = 1.0000000000e-01,  _prvt_stim_end = 1.0000000000e+05,  _prvt_stim_period = 1.0000000000e+00,  _prvt_stim_duration = 1.5000000000e-03,  _prvt_Cm = 9.5000000000e-05,  _prvt_R = 8.3144720000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_F = 9.6485341500e+04,  _prvt_Na_o = 1.4000000000e+02,  _prvt_K_o = 4.0000000000e+00,  _prvt_P_kna = 3.0000000000e-02,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_g_K1 = 5.0000000000e-01,  _prvt_K_mk1 = 1.0000000000e+01,  _prvt_g_Kr1 = 2.1000000000e-03,  _prvt_g_Kr2 = 1.3000000000e-03,  _prvt_g_Ks = 2.6000000000e-03,  _prvt_g_Na = 2.5000000000e+00,  _prvt_delta_m = 1.0000000000e-05,  _prvt_shift_h = 0.0000000000e+00,  _prvt_g_pna = 4.0000000000e-03,  _prvt_g_bna = 6.0000000000e-04,  _prvt_FrICa = 1.0000000000e+00,  _prvt_P_Ca_L = 1.0000000000e-01,  _prvt_P_CaK = 2.0000000000e-03,  _prvt_P_CaNa = 1.0000000000e-02,  _prvt_speed_d = 3.0000000000e+00,  _prvt_delta_f = 1.0000000000e-04,  _prvt_speed_f = 3.0000000000e-01,  _prvt_Km_f2 = 1.0000000000e+05,  _prvt_R_decay = 2.0000000000e+01,  _prvt_Km_f2ds = 1.0000000000e-03,  _prvt_g_bca = 2.5000000000e-04,  _prvt_g_to = 5.0000000000e-03,  _prvt_g_tos = 0.0000000000e+00,  _prvt_i_NaK_max = 7.0000000000e-01,  _prvt_K_mK = 1.0000000000e+00,  _prvt_K_mNa = 4.0000000000e+01,  _prvt_FRiNaCa = 1.0000000000e-03,  _prvt_k_NaCa = 5.0000000000e-04,  _prvt_gamma = 5.0000000000e-01,  _prvt_n_NaCa = 3.0000000000e+00,  _prvt_d_NaCa = 0.0000000000e+00,  _prvt_K_cyca = 3.0000000000e-04,  _prvt_K_xcs = 4.0000000000e-01,  _prvt_K_srca = 5.0000000000e-01,  _prvt_alpha_up = 4.0000000000e-01,  _prvt_beta_up = 3.0000000000e-02,  _prvt_K_m_Ca_cyt = 5.0000000000e-04,  _prvt_K_m_Ca_ds = 1.0000000000e-02,  _prvt_K_m_rel = 2.5000000000e+02,  _prvt_K_leak_rate = 5.0000000000e-02,  _prvt_radius = 1.2000000000e-02,  _prvt_length = 7.4000000000e-02,  _prvt_V_e_ratio = 4.0000000000e-01,  _prvt_V_up_ratio = 1.0000000000e-02,  _prvt_V_rel_ratio = 1.0000000000e-01,  _prvt_V_ds_ratio = 1.0000000000e-01,  _prvt_Kdecay = 1.0000000000e+01,  _prvt_alpha_Calmod = 1.0000000000e+05,  _prvt_Calmod = 2.0000000000e-02,  _prvt_beta_Calmod = 5.0000000000e+01,  _prvt_alpha_Trop = 1.0000000000e+05,  _prvt_Trop = 5.0000000000e-02,  _prvt_beta_Trop = 2.0000000000e+02, 
			//private aux variables
			 _prvt_calc_i_Stim=0.0,  _prvt_calc_E_Na=0.0,  _prvt_calc_E_K=0.0,  _prvt_calc_E_Ks=0.0,  _prvt_calc_E_Ca=0.0,  _prvt_calc_E_mh=0.0,  _prvt_calc_i_K1=0.0,  _prvt_calc_i_Kr=0.0,  _prvt_calc_alpha_xr1=0.0,  _prvt_calc_beta_xr1=0.0,  _prvt_calc_alpha_xr2=0.0,  _prvt_calc_beta_xr2=0.0,  _prvt_calc_i_Ks=0.0,  _prvt_calc_alpha_xs=0.0,  _prvt_calc_beta_xs=0.0,  _prvt_calc_i_Na=0.0,  _prvt_calc_E0_m=0.0,  _prvt_calc_alpha_m=0.0,  _prvt_calc_beta_m=0.0,  _prvt_calc_alpha_h=0.0,  _prvt_calc_beta_h=0.0,  _prvt_calc_i_p_Na=0.0,  _prvt_calc_i_b_Na=0.0,  _prvt_calc_i_Ca_L_Ca_cyt=0.0,  _prvt_calc_i_Ca_L_K_cyt=0.0,  _prvt_calc_i_Ca_L_Na_cyt=0.0,  _prvt_calc_i_Ca_L_Ca_ds=0.0,  _prvt_calc_i_Ca_L_K_ds=0.0,  _prvt_calc_i_Ca_L_Na_ds=0.0,  _prvt_calc_i_Ca_L=0.0,  _prvt_calc_E0_d=0.0,  _prvt_calc_alpha_d=0.0,  _prvt_calc_beta_d=0.0,  _prvt_calc_E0_f=0.0,  _prvt_calc_alpha_f=0.0,  _prvt_calc_beta_f=0.0,  _prvt_calc_i_b_Ca=0.0,  _prvt_calc_i_to=0.0,  _prvt_calc_alpha_s=0.0,  _prvt_calc_beta_s=0.0,  _prvt_calc_i_NaK=0.0,  _prvt_calc_i_NaCa_cyt=0.0,  _prvt_calc_i_NaCa_ds=0.0,  _prvt_calc_i_NaCa=0.0,  _prvt_calc_K_1=0.0,  _prvt_calc_K_2=0.0,  _prvt_calc_i_up=0.0,  _prvt_calc_i_trans=0.0,  _prvt_calc_VoltDep=0.0,  _prvt_calc_CaiReg=0.0,  _prvt_calc_CadsReg=0.0,  _prvt_calc_RegBindSite=0.0,  _prvt_calc_ActRate=0.0,  _prvt_calc_InactRate=0.0,  _prvt_calc_SpeedRel=0.0,  _prvt_calc_PrecFrac=0.0,  _prvt_calc_i_rel=0.0,  _prvt_calc_V_Cell=0.0,  _prvt_calc_V_i_ratio=0.0,  _prvt_calc_V_i=0.0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_xr1_lado_direito_,  _prvt_xr2_lado_direito_,  _prvt_xs_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h_lado_direito_,  _prvt_d_lado_direito_,  _prvt_f_lado_direito_,  _prvt_f2_lado_direito_,  _prvt_f2ds_lado_direito_,  _prvt_s_lado_direito_,  _prvt_r_lado_direito_,  _prvt_ActFrac_lado_direito_,  _prvt_ProdFrac_lado_direito_,  _prvt_Na_i_lado_direito_,  _prvt_K_i_lado_direito_,  _prvt_Ca_i_lado_direito_,  _prvt_Ca_ds_lado_direito_,  _prvt_Ca_up_lado_direito_,  _prvt_Ca_rel_lado_direito_,  _prvt_Ca_Calmod_lado_direito_,  _prvt_Ca_Trop_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;
			__NEW_[0] = __OLD_[0] = -9.2849333000e+01;
			__NEW_[1] = __OLD_[1] = 1.0300000000e-05;
			__NEW_[2] = __OLD_[2] = 2.0000000000e-07;
			__NEW_[3] = __OLD_[3] = 1.3020000000e-03;
			__NEW_[4] = __OLD_[4] = 1.6203000000e-03;
			__NEW_[5] = __OLD_[5] = 9.9440360000e-01;
			__NEW_[6] = __OLD_[6] = 0.0000000000e+00;
			__NEW_[7] = __OLD_[7] = 1.0000000000e+00;
			__NEW_[8] = __OLD_[8] = 9.3491970000e-01;
			__NEW_[9] = __OLD_[9] = 9.6519580000e-01;
			__NEW_[10] = __OLD_[10] = 9.9486450000e-01;
			__NEW_[11] = __OLD_[11] = 0.0000000000e+00;
			__NEW_[12] = __OLD_[12] = 4.2614000000e-03;
			__NEW_[13] = __OLD_[13] = 4.0681540000e-01;
			__NEW_[14] = __OLD_[14] = 7.3321223000e+00;
			__NEW_[15] = __OLD_[15] = 1.3656442810e+02;
			__NEW_[16] = __OLD_[16] = 1.4000000000e-05;
			__NEW_[17] = __OLD_[17] = 1.8800000000e-05;
			__NEW_[18] = __OLD_[18] = 4.5318890000e-01;
			__NEW_[19] = __OLD_[19] = 4.4819270000e-01;
			__NEW_[20] = __OLD_[20] = 5.5550000000e-04;
			__NEW_[21] = __OLD_[21] = 3.5420000000e-04;
			_prvt_time_new += _prvt_dtime;
			if(omp_get_thread_num()==_prvt_tree_thread[0])
			{
				_prvt_calc_i_Stim = (((_prvt_time_new>=_prvt_stim_start)&&(_prvt_time_new<=_prvt_stim_end)&&(((_prvt_time_new-_prvt_stim_start)-(floor(((_prvt_time_new-_prvt_stim_start)/_prvt_stim_period))*_prvt_stim_period))<=_prvt_stim_duration)))
?(_prvt_stim_amplitude)
:(0.0000000000e+00);
				_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[14])));
				_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[15])));
				_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[15]+(_prvt_P_kna*__OLD_[14])))));
				_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[16])));
				_prvt_calc_E_mh = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_Na_o+(1.2000000000e-01*_prvt_K_o))/(__OLD_[14]+(1.2000000000e-01*__OLD_[15])))));
				_prvt_calc_i_Ca_L_Ca_cyt = (((((1.0000000000e+00-_prvt_FrICa)*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_K_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_Na_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_Ca_ds = ((((_prvt_FrICa*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_K_ds = ((((_prvt_FrICa*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_Ca_L_Na_ds = ((((_prvt_FrICa*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
				_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_K_o)/(_prvt_K_mK+_prvt_K_o))*__OLD_[14])/(_prvt_K_mNa+__OLD_[14]));
				_prvt_calc_i_NaCa_cyt = (((1.0000000000e+00-_prvt_FRiNaCa)*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[16])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[16]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[16]/6.9000000000e-03))));
				_prvt_calc_i_NaCa_ds = ((_prvt_FRiNaCa*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[17])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[17]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[17]/6.9000000000e-03))));
				_prvt_calc_K_1 = ((_prvt_K_cyca*_prvt_K_xcs)/_prvt_K_srca);
				_prvt_calc_i_trans = (5.0000000000e+01*(__OLD_[18]-__OLD_[19]));
				_prvt_calc_i_rel = (((pow((__OLD_[12]/(__OLD_[12]+2.5000000000e-01)),2.0000000000e+00)*_prvt_K_m_rel)+_prvt_K_leak_rate)*__OLD_[19]);
				_prvt_calc_V_Cell = (3.1415926540e+00*pow(_prvt_radius,2.0000000000e+00)*_prvt_length);
				_prvt_calc_V_i_ratio = (((1.0000000000e+00-_prvt_V_e_ratio)-_prvt_V_up_ratio)-_prvt_V_rel_ratio);
				_prvt_calc_i_K1 = ((((_prvt_g_K1*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk1))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01)*_prvt_F*1.2500000000e+00)/(_prvt_R*_prvt_T)))));
				_prvt_calc_i_Kr = (((((_prvt_g_Kr1*__OLD_[1])+(_prvt_g_Kr2*__OLD_[2]))*1.0000000000e+00)/(1.0000000000e+00+exp(((__OLD_[0]+9.0000000000e+00)/2.2400000000e+01))))*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
				_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*(__OLD_[0]-_prvt_calc_E_mh));
				_prvt_calc_i_p_Na = (((_prvt_g_pna*1.0000000000e+00)/(1.0000000000e+00+exp(((-(__OLD_[0]+5.2000000000e+01))/8.0000000000e+00))))*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_Ca_L = (_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_Ca_L_Na_ds);
				_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
				_prvt_calc_i_to = (_prvt_g_to*(_prvt_g_tos+(__OLD_[10]*(1.0000000000e+00-_prvt_g_tos)))*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_NaCa = (_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds);
				_prvt_calc_K_2 = (__OLD_[16]+(__OLD_[18]*_prvt_calc_K_1)+(_prvt_K_cyca*_prvt_K_xcs)+_prvt_K_cyca);
				_prvt_calc_V_i = (_prvt_calc_V_Cell*_prvt_calc_V_i_ratio);
				_prvt_calc_i_up = (((__OLD_[16]/_prvt_calc_K_2)*_prvt_alpha_up)-(((__OLD_[18]*_prvt_calc_K_1)/_prvt_calc_K_2)*_prvt_beta_up));
				__K1_[0]= (((-1.0000000000e+00)/_prvt_Cm)*(_prvt_calc_i_Stim+_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_p_Na+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds+_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds+_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_b_Ca));
				__NEW_[0]= __K1_[0] * _prvt_dtime + __OLD_[0];
				__K1_[14]= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(_prvt_calc_i_Na+_prvt_calc_i_p_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa_cyt)+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds));
				__NEW_[14]= __K1_[14] * _prvt_dtime + __OLD_[14];
				__K1_[15]= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*((_prvt_calc_i_K1+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_to)-(2.0000000000e+00*_prvt_calc_i_NaK)));
				__NEW_[15]= __K1_[15] * _prvt_dtime + __OLD_[15];
				__K1_[16]= (((((((-1.0000000000e+00)/(2.0000000000e+00*1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(((_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_b_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa_cyt))-(2.0000000000e+00*_prvt_calc_i_NaCa_ds)))+(__OLD_[17]*_prvt_V_ds_ratio*_prvt_Kdecay)+((_prvt_calc_i_rel*_prvt_V_rel_ratio)/_prvt_calc_V_i_ratio))-((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20])))-((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21])))-_prvt_calc_i_up);
				__NEW_[16]= __K1_[16] * _prvt_dtime + __OLD_[16];
				__K1_[17]= ((((-1.0000000000e+00)*_prvt_calc_i_Ca_L_Ca_ds)/(2.0000000000e+00*1.0000000000e+00*_prvt_V_ds_ratio*_prvt_calc_V_i*_prvt_F))-(__OLD_[17]*_prvt_Kdecay));
				__NEW_[17]= __K1_[17] * _prvt_dtime + __OLD_[17];
				__K1_[18]= (((_prvt_calc_V_i_ratio/_prvt_V_up_ratio)*_prvt_calc_i_up)-_prvt_calc_i_trans);
				__NEW_[18]= __K1_[18] * _prvt_dtime + __OLD_[18];
				__K1_[19]= (((_prvt_V_up_ratio/_prvt_V_rel_ratio)*_prvt_calc_i_trans)-_prvt_calc_i_rel);
				__NEW_[19]= __K1_[19] * _prvt_dtime + __OLD_[19];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[1])
			{
				_prvt_calc_alpha_xr1 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
				_prvt_calc_beta_xr1 = (5.0000000000e-02*exp(((-(__OLD_[0]-2.0000000000e+01))/1.5000000000e+01)));
				__K1_[1]= ((_prvt_calc_alpha_xr1*(1.0000000000e+00-__OLD_[1]))-(_prvt_calc_beta_xr1*__OLD_[1]));
				__NEW_[1]= __K1_[1] * _prvt_dtime + __OLD_[1];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[2])
			{
				_prvt_calc_alpha_xr2 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
				_prvt_calc_beta_xr2 = (4.0000000000e-01*exp((-pow(((__OLD_[0]+3.0000000000e+01)/3.0000000000e+01),3.0000000000e+00))));
				__K1_[2]= ((_prvt_calc_alpha_xr2*(1.0000000000e+00-__OLD_[2]))-(_prvt_calc_beta_xr2*__OLD_[2]));
				__NEW_[2]= __K1_[2] * _prvt_dtime + __OLD_[2];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[3])
			{
				_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
				_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
				__K1_[3]= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[3]))-(_prvt_calc_beta_xs*__OLD_[3]));
				__NEW_[3]= __K1_[3] * _prvt_dtime + __OLD_[3];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[4])
			{
				_prvt_calc_E0_m = (__OLD_[0]+4.1000000000e+01);
				_prvt_calc_beta_m = (8.0000000000e+03*exp(((-5.6000000000e-02)*(__OLD_[0]+6.6000000000e+01))));
				_prvt_calc_alpha_m = ((fabs(_prvt_calc_E0_m)<_prvt_delta_m))
?(2.0000000000e+03)
:(((2.0000000000e+02*_prvt_calc_E0_m)/(1.0000000000e+00-exp(((-1.0000000000e-01)*_prvt_calc_E0_m)))));
				__K1_[4]= ((_prvt_calc_alpha_m*(1.0000000000e+00-__OLD_[4]))-(_prvt_calc_beta_m*__OLD_[4]));
				__NEW_[4]= __K1_[4] * _prvt_dtime + __OLD_[4];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[5])
			{
				_prvt_calc_alpha_h = (2.0000000000e+01*exp(((-1.2500000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))));
				_prvt_calc_beta_h = (2.0000000000e+03/(1.0000000000e+00+(3.2000000000e+02*exp(((-1.0000000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))))));
				__K1_[5]= ((_prvt_calc_alpha_h*(1.0000000000e+00-__OLD_[5]))-(_prvt_calc_beta_h*__OLD_[5]));
				__NEW_[5]= __K1_[5] * _prvt_dtime + __OLD_[5];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[6])
			{
				_prvt_calc_E0_d = ((__OLD_[0]+2.4000000000e+01)-5.0000000000e+00);
				_prvt_calc_alpha_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((3.0000000000e+01*_prvt_calc_E0_d)/(1.0000000000e+00-exp(((-_prvt_calc_E0_d)/4.0000000000e+00)))));
				_prvt_calc_beta_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((1.2000000000e+01*_prvt_calc_E0_d)/(exp((_prvt_calc_E0_d/1.0000000000e+01))-1.0000000000e+00)));
				__K1_[6]= (_prvt_speed_d*((_prvt_calc_alpha_d*(1.0000000000e+00-__OLD_[6]))-(_prvt_calc_beta_d*__OLD_[6])));
				__NEW_[6]= __K1_[6] * _prvt_dtime + __OLD_[6];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[7])
			{
				_prvt_calc_E0_f = (__OLD_[0]+3.4000000000e+01);
				_prvt_calc_beta_f = (1.2000000000e+01/(1.0000000000e+00+exp((((-1.0000000000e+00)*(__OLD_[0]+3.4000000000e+01))/4.0000000000e+00))));
				_prvt_calc_alpha_f = ((fabs(_prvt_calc_E0_f)<_prvt_delta_f))
?(2.5000000000e+01)
:(((6.2500000000e+00*_prvt_calc_E0_f)/(exp((_prvt_calc_E0_f/4.0000000000e+00))-1.0000000000e+00)));
				__K1_[7]= (_prvt_speed_f*((_prvt_calc_alpha_f*(1.0000000000e+00-__OLD_[7]))-(_prvt_calc_beta_f*__OLD_[7])));
				__NEW_[7]= __K1_[7] * _prvt_dtime + __OLD_[7];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[8])
			{
				_prvt_calc_alpha_s = (3.3000000000e-02*exp(((-__OLD_[0])/1.7000000000e+01)));
				_prvt_calc_beta_s = (3.3000000000e+01/(1.0000000000e+00+exp(((-1.2500000000e-01)*(__OLD_[0]+1.0000000000e+01)))));
				__K1_[10]= ((_prvt_calc_alpha_s*(1.0000000000e+00-__OLD_[10]))-(_prvt_calc_beta_s*__OLD_[10]));
				__NEW_[10]= __K1_[10] * _prvt_dtime + __OLD_[10];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[9])
			{
				_prvt_calc_VoltDep = exp((8.0000000000e-02*(__OLD_[0]-4.0000000000e+01)));
				_prvt_calc_CaiReg = (__OLD_[16]/(__OLD_[16]+_prvt_K_m_Ca_cyt));
				_prvt_calc_CadsReg = (__OLD_[17]/(__OLD_[17]+_prvt_K_m_Ca_ds));
				_prvt_calc_SpeedRel = ((__OLD_[0]<(-5.0000000000e+01)))
?(5.0000000000e+00)
:(1.0000000000e+00);
				_prvt_calc_PrecFrac = ((1.0000000000e+00-__OLD_[12])-__OLD_[13]);
				_prvt_calc_RegBindSite = (_prvt_calc_CaiReg+((1.0000000000e+00-_prvt_calc_CaiReg)*_prvt_calc_CadsReg));
				_prvt_calc_ActRate = ((0.0000000000e+00*_prvt_calc_VoltDep)+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
				_prvt_calc_InactRate = (6.0000000000e+01+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
				__K1_[12]= ((_prvt_calc_PrecFrac*_prvt_calc_SpeedRel*_prvt_calc_ActRate)-(__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate));
				__NEW_[12]= __K1_[12] * _prvt_dtime + __OLD_[12];
				__K1_[13]= ((__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate)-(_prvt_calc_SpeedRel*1.0000000000e+00*__OLD_[13]));
				__NEW_[13]= __K1_[13] * _prvt_dtime + __OLD_[13];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[10])
			{
				__K1_[8]= (1.0000000000e+00-(1.0000000000e+00*((__OLD_[16]/(_prvt_Km_f2+__OLD_[16]))+__OLD_[8])));
				__NEW_[8]= __K1_[8] * _prvt_dtime + __OLD_[8];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[11])
			{
				__K1_[9]= (_prvt_R_decay*(1.0000000000e+00-((__OLD_[17]/(_prvt_Km_f2ds+__OLD_[17]))+__OLD_[9])));
				__NEW_[9]= __K1_[9] * _prvt_dtime + __OLD_[9];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[12])
			{
				__K1_[11]= (3.3300000000e+02*((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+4.0000000000e+00))/5.0000000000e+00))))-__OLD_[11]));
				__NEW_[11]= __K1_[11] * _prvt_dtime + __OLD_[11];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[13])
			{
				__K1_[20]= ((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20]));
				__NEW_[20]= __K1_[20] * _prvt_dtime + __OLD_[20];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[14])
			{
				__K1_[21]= ((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21]));
				__NEW_[21]= __K1_[21] * _prvt_dtime + __OLD_[21];
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
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/__OLD_[14])));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/__OLD_[15])));
					_prvt_calc_E_Ks = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(_prvt_P_kna*_prvt_Na_o))/(__OLD_[15]+(_prvt_P_kna*__OLD_[14])))));
					_prvt_calc_E_Ca = (((5.0000000000e-01*_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Ca_o/__OLD_[16])));
					_prvt_calc_E_mh = (((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_Na_o+(1.2000000000e-01*_prvt_K_o))/(__OLD_[14]+(1.2000000000e-01*__OLD_[15])))));
					_prvt_calc_i_Ca_L_Ca_cyt = (((((1.0000000000e+00-_prvt_FrICa)*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_K_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Na_cyt = (((((1.0000000000e+00-_prvt_FrICa)*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[8]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Ca_ds = ((((_prvt_FrICa*4.0000000000e+00*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T)))))*((__OLD_[16]*exp(((1.0000000000e+02*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Ca_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F*2.0000000000e+00)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_K_ds = ((((_prvt_FrICa*_prvt_P_CaK*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[15]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_K_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_Ca_L_Na_ds = ((((_prvt_FrICa*_prvt_P_CaNa*_prvt_P_Ca_L*__OLD_[6]*__OLD_[7]*__OLD_[9]*(__OLD_[0]-5.0000000000e+01)*_prvt_F)/(_prvt_R*_prvt_T))/(1.0000000000e+00-exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T)))))*((__OLD_[14]*exp(((5.0000000000e+01*_prvt_F)/(_prvt_R*_prvt_T))))-(_prvt_Na_o*exp((((-(__OLD_[0]-5.0000000000e+01))*_prvt_F)/(_prvt_R*_prvt_T))))));
					_prvt_calc_i_NaK = ((((_prvt_i_NaK_max*_prvt_K_o)/(_prvt_K_mK+_prvt_K_o))*__OLD_[14])/(_prvt_K_mNa+__OLD_[14]));
					_prvt_calc_i_NaCa_cyt = (((1.0000000000e+00-_prvt_FRiNaCa)*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[16])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[16]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[16]/6.9000000000e-03))));
					_prvt_calc_i_NaCa_ds = ((_prvt_FRiNaCa*_prvt_k_NaCa*((exp(((_prvt_gamma*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(__OLD_[14],_prvt_n_NaCa)*_prvt_Ca_o)-(exp((((_prvt_gamma-1.0000000000e+00)*(_prvt_n_NaCa-2.0000000000e+00)*__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))*pow(_prvt_Na_o,_prvt_n_NaCa)*__OLD_[17])))/((1.0000000000e+00+(_prvt_d_NaCa*((__OLD_[17]*pow(_prvt_Na_o,_prvt_n_NaCa))+(_prvt_Ca_o*pow(__OLD_[14],_prvt_n_NaCa)))))*(1.0000000000e+00+(__OLD_[17]/6.9000000000e-03))));
					_prvt_calc_K_1 = ((_prvt_K_cyca*_prvt_K_xcs)/_prvt_K_srca);
					_prvt_calc_i_trans = (5.0000000000e+01*(__OLD_[18]-__OLD_[19]));
					_prvt_calc_i_rel = (((pow((__OLD_[12]/(__OLD_[12]+2.5000000000e-01)),2.0000000000e+00)*_prvt_K_m_rel)+_prvt_K_leak_rate)*__OLD_[19]);
					_prvt_calc_V_Cell = (3.1415926540e+00*pow(_prvt_radius,2.0000000000e+00)*_prvt_length);
					_prvt_calc_V_i_ratio = (((1.0000000000e+00-_prvt_V_e_ratio)-_prvt_V_up_ratio)-_prvt_V_rel_ratio);
					_prvt_calc_i_K1 = ((((_prvt_g_K1*_prvt_K_o)/(_prvt_K_o+_prvt_K_mk1))*(__OLD_[0]-_prvt_calc_E_K))/(1.0000000000e+00+exp(((((__OLD_[0]-_prvt_calc_E_K)-1.0000000000e+01)*_prvt_F*1.2500000000e+00)/(_prvt_R*_prvt_T)))));
					_prvt_calc_i_Kr = (((((_prvt_g_Kr1*__OLD_[1])+(_prvt_g_Kr2*__OLD_[2]))*1.0000000000e+00)/(1.0000000000e+00+exp(((__OLD_[0]+9.0000000000e+00)/2.2400000000e+01))))*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_Ks = (_prvt_g_Ks*pow(__OLD_[3],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_Ks));
					_prvt_calc_i_Na = (_prvt_g_Na*pow(__OLD_[4],3.0000000000e+00)*__OLD_[5]*(__OLD_[0]-_prvt_calc_E_mh));
					_prvt_calc_i_p_Na = (((_prvt_g_pna*1.0000000000e+00)/(1.0000000000e+00+exp(((-(__OLD_[0]+5.2000000000e+01))/8.0000000000e+00))))*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_Na = (_prvt_g_bna*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_Ca_L = (_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_Ca_L_Na_ds);
					_prvt_calc_i_b_Ca = (_prvt_g_bca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_to = (_prvt_g_to*(_prvt_g_tos+(__OLD_[10]*(1.0000000000e+00-_prvt_g_tos)))*__OLD_[11]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_NaCa = (_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds);
					_prvt_calc_K_2 = (__OLD_[16]+(__OLD_[18]*_prvt_calc_K_1)+(_prvt_K_cyca*_prvt_K_xcs)+_prvt_K_cyca);
					_prvt_calc_V_i = (_prvt_calc_V_Cell*_prvt_calc_V_i_ratio);
					_prvt_calc_i_up = (((__OLD_[16]/_prvt_calc_K_2)*_prvt_alpha_up)-(((__OLD_[18]*_prvt_calc_K_1)/_prvt_calc_K_2)*_prvt_beta_up));
					__K2_[0]= (((-1.0000000000e+00)/_prvt_Cm)*(_prvt_calc_i_Stim+_prvt_calc_i_K1+_prvt_calc_i_to+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_NaK+_prvt_calc_i_Na+_prvt_calc_i_b_Na+_prvt_calc_i_p_Na+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds+_prvt_calc_i_NaCa_cyt+_prvt_calc_i_NaCa_ds+_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_Ca_L_Ca_ds+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_b_Ca));
					_prvt_aux_tol = fabs(__OLD_[0])*_prvt_rel_tol_;
					__TOL_[0] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[0] = fabs((_prvt_dtime/2) * (__K1_[0] - __K2_[0])/__TOL_[0]);
					__K2_[14]= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(_prvt_calc_i_Na+_prvt_calc_i_p_Na+_prvt_calc_i_b_Na+(3.0000000000e+00*_prvt_calc_i_NaK)+(3.0000000000e+00*_prvt_calc_i_NaCa_cyt)+_prvt_calc_i_Ca_L_Na_cyt+_prvt_calc_i_Ca_L_Na_ds));
					_prvt_aux_tol = fabs(__OLD_[14])*_prvt_rel_tol_;
					__TOL_[14] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[14] = fabs((_prvt_dtime/2) * (__K1_[14] - __K2_[14])/__TOL_[14]);
					__K2_[15]= (((-1.0000000000e+00)/(1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*((_prvt_calc_i_K1+_prvt_calc_i_Kr+_prvt_calc_i_Ks+_prvt_calc_i_Ca_L_K_cyt+_prvt_calc_i_Ca_L_K_ds+_prvt_calc_i_to)-(2.0000000000e+00*_prvt_calc_i_NaK)));
					_prvt_aux_tol = fabs(__OLD_[15])*_prvt_rel_tol_;
					__TOL_[15] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[15] = fabs((_prvt_dtime/2) * (__K1_[15] - __K2_[15])/__TOL_[15]);
					__K2_[16]= (((((((-1.0000000000e+00)/(2.0000000000e+00*1.0000000000e+00*_prvt_calc_V_i*_prvt_F))*(((_prvt_calc_i_Ca_L_Ca_cyt+_prvt_calc_i_b_Ca)-(2.0000000000e+00*_prvt_calc_i_NaCa_cyt))-(2.0000000000e+00*_prvt_calc_i_NaCa_ds)))+(__OLD_[17]*_prvt_V_ds_ratio*_prvt_Kdecay)+((_prvt_calc_i_rel*_prvt_V_rel_ratio)/_prvt_calc_V_i_ratio))-((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20])))-((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21])))-_prvt_calc_i_up);
					_prvt_aux_tol = fabs(__OLD_[16])*_prvt_rel_tol_;
					__TOL_[16] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[16] = fabs((_prvt_dtime/2) * (__K1_[16] - __K2_[16])/__TOL_[16]);
					__K2_[17]= ((((-1.0000000000e+00)*_prvt_calc_i_Ca_L_Ca_ds)/(2.0000000000e+00*1.0000000000e+00*_prvt_V_ds_ratio*_prvt_calc_V_i*_prvt_F))-(__OLD_[17]*_prvt_Kdecay));
					_prvt_aux_tol = fabs(__OLD_[17])*_prvt_rel_tol_;
					__TOL_[17] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[17] = fabs((_prvt_dtime/2) * (__K1_[17] - __K2_[17])/__TOL_[17]);
					__K2_[18]= (((_prvt_calc_V_i_ratio/_prvt_V_up_ratio)*_prvt_calc_i_up)-_prvt_calc_i_trans);
					_prvt_aux_tol = fabs(__OLD_[18])*_prvt_rel_tol_;
					__TOL_[18] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[18] = fabs((_prvt_dtime/2) * (__K1_[18] - __K2_[18])/__TOL_[18]);
					__K2_[19]= (((_prvt_V_up_ratio/_prvt_V_rel_ratio)*_prvt_calc_i_trans)-_prvt_calc_i_rel);
					_prvt_aux_tol = fabs(__OLD_[19])*_prvt_rel_tol_;
					__TOL_[19] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[19] = fabs((_prvt_dtime/2) * (__K1_[19] - __K2_[19])/__TOL_[19]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[1])
				{
					_prvt_calc_alpha_xr1 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
					_prvt_calc_beta_xr1 = (5.0000000000e-02*exp(((-(__OLD_[0]-2.0000000000e+01))/1.5000000000e+01)));
					__K2_[1]= ((_prvt_calc_alpha_xr1*(1.0000000000e+00-__OLD_[1]))-(_prvt_calc_beta_xr1*__OLD_[1]));
					_prvt_aux_tol = fabs(__OLD_[1])*_prvt_rel_tol_;
					__TOL_[1] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[1] = fabs((_prvt_dtime/2) * (__K1_[1] - __K2_[1])/__TOL_[1]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[2])
				{
					_prvt_calc_alpha_xr2 = (5.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-5.0000000000e+00))/9.0000000000e+00))));
					_prvt_calc_beta_xr2 = (4.0000000000e-01*exp((-pow(((__OLD_[0]+3.0000000000e+01)/3.0000000000e+01),3.0000000000e+00))));
					__K2_[2]= ((_prvt_calc_alpha_xr2*(1.0000000000e+00-__OLD_[2]))-(_prvt_calc_beta_xr2*__OLD_[2]));
					_prvt_aux_tol = fabs(__OLD_[2])*_prvt_rel_tol_;
					__TOL_[2] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[2] = fabs((_prvt_dtime/2) * (__K1_[2] - __K2_[2])/__TOL_[2]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[3])
				{
					_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
					_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
					__K2_[3]= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[3]))-(_prvt_calc_beta_xs*__OLD_[3]));
					_prvt_aux_tol = fabs(__OLD_[3])*_prvt_rel_tol_;
					__TOL_[3] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[3] = fabs((_prvt_dtime/2) * (__K1_[3] - __K2_[3])/__TOL_[3]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[4])
				{
					_prvt_calc_E0_m = (__OLD_[0]+4.1000000000e+01);
					_prvt_calc_beta_m = (8.0000000000e+03*exp(((-5.6000000000e-02)*(__OLD_[0]+6.6000000000e+01))));
					_prvt_calc_alpha_m = ((fabs(_prvt_calc_E0_m)<_prvt_delta_m))
?(2.0000000000e+03)
:(((2.0000000000e+02*_prvt_calc_E0_m)/(1.0000000000e+00-exp(((-1.0000000000e-01)*_prvt_calc_E0_m)))));
					__K2_[4]= ((_prvt_calc_alpha_m*(1.0000000000e+00-__OLD_[4]))-(_prvt_calc_beta_m*__OLD_[4]));
					_prvt_aux_tol = fabs(__OLD_[4])*_prvt_rel_tol_;
					__TOL_[4] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[4] = fabs((_prvt_dtime/2) * (__K1_[4] - __K2_[4])/__TOL_[4]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[5])
				{
					_prvt_calc_alpha_h = (2.0000000000e+01*exp(((-1.2500000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))));
					_prvt_calc_beta_h = (2.0000000000e+03/(1.0000000000e+00+(3.2000000000e+02*exp(((-1.0000000000e-01)*((__OLD_[0]+7.5000000000e+01)-_prvt_shift_h))))));
					__K2_[5]= ((_prvt_calc_alpha_h*(1.0000000000e+00-__OLD_[5]))-(_prvt_calc_beta_h*__OLD_[5]));
					_prvt_aux_tol = fabs(__OLD_[5])*_prvt_rel_tol_;
					__TOL_[5] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[5] = fabs((_prvt_dtime/2) * (__K1_[5] - __K2_[5])/__TOL_[5]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[6])
				{
					_prvt_calc_E0_d = ((__OLD_[0]+2.4000000000e+01)-5.0000000000e+00);
					_prvt_calc_alpha_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((3.0000000000e+01*_prvt_calc_E0_d)/(1.0000000000e+00-exp(((-_prvt_calc_E0_d)/4.0000000000e+00)))));
					_prvt_calc_beta_d = ((fabs(_prvt_calc_E0_d)<1.0000000000e-04))
?(1.2000000000e+02)
:(((1.2000000000e+01*_prvt_calc_E0_d)/(exp((_prvt_calc_E0_d/1.0000000000e+01))-1.0000000000e+00)));
					__K2_[6]= (_prvt_speed_d*((_prvt_calc_alpha_d*(1.0000000000e+00-__OLD_[6]))-(_prvt_calc_beta_d*__OLD_[6])));
					_prvt_aux_tol = fabs(__OLD_[6])*_prvt_rel_tol_;
					__TOL_[6] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[6] = fabs((_prvt_dtime/2) * (__K1_[6] - __K2_[6])/__TOL_[6]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[7])
				{
					_prvt_calc_E0_f = (__OLD_[0]+3.4000000000e+01);
					_prvt_calc_beta_f = (1.2000000000e+01/(1.0000000000e+00+exp((((-1.0000000000e+00)*(__OLD_[0]+3.4000000000e+01))/4.0000000000e+00))));
					_prvt_calc_alpha_f = ((fabs(_prvt_calc_E0_f)<_prvt_delta_f))
?(2.5000000000e+01)
:(((6.2500000000e+00*_prvt_calc_E0_f)/(exp((_prvt_calc_E0_f/4.0000000000e+00))-1.0000000000e+00)));
					__K2_[7]= (_prvt_speed_f*((_prvt_calc_alpha_f*(1.0000000000e+00-__OLD_[7]))-(_prvt_calc_beta_f*__OLD_[7])));
					_prvt_aux_tol = fabs(__OLD_[7])*_prvt_rel_tol_;
					__TOL_[7] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[7] = fabs((_prvt_dtime/2) * (__K1_[7] - __K2_[7])/__TOL_[7]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[8])
				{
					_prvt_calc_alpha_s = (3.3000000000e-02*exp(((-__OLD_[0])/1.7000000000e+01)));
					_prvt_calc_beta_s = (3.3000000000e+01/(1.0000000000e+00+exp(((-1.2500000000e-01)*(__OLD_[0]+1.0000000000e+01)))));
					__K2_[10]= ((_prvt_calc_alpha_s*(1.0000000000e+00-__OLD_[10]))-(_prvt_calc_beta_s*__OLD_[10]));
					_prvt_aux_tol = fabs(__OLD_[10])*_prvt_rel_tol_;
					__TOL_[10] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[10] = fabs((_prvt_dtime/2) * (__K1_[10] - __K2_[10])/__TOL_[10]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[9])
				{
					_prvt_calc_VoltDep = exp((8.0000000000e-02*(__OLD_[0]-4.0000000000e+01)));
					_prvt_calc_CaiReg = (__OLD_[16]/(__OLD_[16]+_prvt_K_m_Ca_cyt));
					_prvt_calc_CadsReg = (__OLD_[17]/(__OLD_[17]+_prvt_K_m_Ca_ds));
					_prvt_calc_SpeedRel = ((__OLD_[0]<(-5.0000000000e+01)))
?(5.0000000000e+00)
:(1.0000000000e+00);
					_prvt_calc_PrecFrac = ((1.0000000000e+00-__OLD_[12])-__OLD_[13]);
					_prvt_calc_RegBindSite = (_prvt_calc_CaiReg+((1.0000000000e+00-_prvt_calc_CaiReg)*_prvt_calc_CadsReg));
					_prvt_calc_ActRate = ((0.0000000000e+00*_prvt_calc_VoltDep)+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
					_prvt_calc_InactRate = (6.0000000000e+01+(5.0000000000e+02*pow(_prvt_calc_RegBindSite,2.0000000000e+00)));
					__K2_[12]= ((_prvt_calc_PrecFrac*_prvt_calc_SpeedRel*_prvt_calc_ActRate)-(__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate));
					_prvt_aux_tol = fabs(__OLD_[12])*_prvt_rel_tol_;
					__TOL_[12] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[12] = fabs((_prvt_dtime/2) * (__K1_[12] - __K2_[12])/__TOL_[12]);
					__K2_[13]= ((__OLD_[12]*_prvt_calc_SpeedRel*_prvt_calc_InactRate)-(_prvt_calc_SpeedRel*1.0000000000e+00*__OLD_[13]));
					_prvt_aux_tol = fabs(__OLD_[13])*_prvt_rel_tol_;
					__TOL_[13] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[13] = fabs((_prvt_dtime/2) * (__K1_[13] - __K2_[13])/__TOL_[13]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[10])
				{
					__K2_[8]= (1.0000000000e+00-(1.0000000000e+00*((__OLD_[16]/(_prvt_Km_f2+__OLD_[16]))+__OLD_[8])));
					_prvt_aux_tol = fabs(__OLD_[8])*_prvt_rel_tol_;
					__TOL_[8] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[8] = fabs((_prvt_dtime/2) * (__K1_[8] - __K2_[8])/__TOL_[8]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[11])
				{
					__K2_[9]= (_prvt_R_decay*(1.0000000000e+00-((__OLD_[17]/(_prvt_Km_f2ds+__OLD_[17]))+__OLD_[9])));
					_prvt_aux_tol = fabs(__OLD_[9])*_prvt_rel_tol_;
					__TOL_[9] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[9] = fabs((_prvt_dtime/2) * (__K1_[9] - __K2_[9])/__TOL_[9]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[12])
				{
					__K2_[11]= (3.3300000000e+02*((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+4.0000000000e+00))/5.0000000000e+00))))-__OLD_[11]));
					_prvt_aux_tol = fabs(__OLD_[11])*_prvt_rel_tol_;
					__TOL_[11] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[11] = fabs((_prvt_dtime/2) * (__K1_[11] - __K2_[11])/__TOL_[11]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[13])
				{
					__K2_[20]= ((_prvt_alpha_Calmod*__OLD_[16]*(_prvt_Calmod-__OLD_[20]))-(_prvt_beta_Calmod*__OLD_[20]));
					_prvt_aux_tol = fabs(__OLD_[20])*_prvt_rel_tol_;
					__TOL_[20] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[20] = fabs((_prvt_dtime/2) * (__K1_[20] - __K2_[20])/__TOL_[20]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[14])
				{
					__K2_[21]= ((_prvt_alpha_Trop*__OLD_[16]*(_prvt_Trop-__OLD_[21]))-(_prvt_beta_Trop*__OLD_[21]));
					_prvt_aux_tol = fabs(__OLD_[21])*_prvt_rel_tol_;
					__TOL_[21] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[21] = fabs((_prvt_dtime/2) * (__K1_[21] - __K2_[21])/__TOL_[21]);
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
						__NEW_[14] = __K1_[14] * _prvt_dtime + __OLD_AUX_[14];
						__NEW_[15] = __K1_[15] * _prvt_dtime + __OLD_AUX_[15];
						__NEW_[16] = __K1_[16] * _prvt_dtime + __OLD_AUX_[16];
						__NEW_[17] = __K1_[17] * _prvt_dtime + __OLD_AUX_[17];
						__NEW_[18] = __K1_[18] * _prvt_dtime + __OLD_AUX_[18];
						__NEW_[19] = __K1_[19] * _prvt_dtime + __OLD_AUX_[19];
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
						__NEW_[10] = __K1_[10] * _prvt_dtime + __OLD_AUX_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[12] = __K1_[12] * _prvt_dtime + __OLD_AUX_[12];
						__NEW_[13] = __K1_[13] * _prvt_dtime + __OLD_AUX_[13];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[8] = __K1_[8] * _prvt_dtime + __OLD_AUX_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[9] = __K1_[9] * _prvt_dtime + __OLD_AUX_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[11] = __K1_[11] * _prvt_dtime + __OLD_AUX_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[13])
					{
						__NEW_[20] = __K1_[20] * _prvt_dtime + __OLD_AUX_[20];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[14])
					{
						__NEW_[21] = __K1_[21] * _prvt_dtime + __OLD_AUX_[21];
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
							this->xr1_old_ = __OLD_AUX_[1];
							this->xr1_new_ = __OLD_[1];
							this->xr2_old_ = __OLD_AUX_[2];
							this->xr2_new_ = __OLD_[2];
							this->xs_old_ = __OLD_AUX_[3];
							this->xs_new_ = __OLD_[3];
							this->m_old_ = __OLD_AUX_[4];
							this->m_new_ = __OLD_[4];
							this->h_old_ = __OLD_AUX_[5];
							this->h_new_ = __OLD_[5];
							this->d_old_ = __OLD_AUX_[6];
							this->d_new_ = __OLD_[6];
							this->f_old_ = __OLD_AUX_[7];
							this->f_new_ = __OLD_[7];
							this->f2_old_ = __OLD_AUX_[8];
							this->f2_new_ = __OLD_[8];
							this->f2ds_old_ = __OLD_AUX_[9];
							this->f2ds_new_ = __OLD_[9];
							this->s_old_ = __OLD_AUX_[10];
							this->s_new_ = __OLD_[10];
							this->r_old_ = __OLD_AUX_[11];
							this->r_new_ = __OLD_[11];
							this->ActFrac_old_ = __OLD_AUX_[12];
							this->ActFrac_new_ = __OLD_[12];
							this->ProdFrac_old_ = __OLD_AUX_[13];
							this->ProdFrac_new_ = __OLD_[13];
							this->Na_i_old_ = __OLD_AUX_[14];
							this->Na_i_new_ = __OLD_[14];
							this->K_i_old_ = __OLD_AUX_[15];
							this->K_i_new_ = __OLD_[15];
							this->Ca_i_old_ = __OLD_AUX_[16];
							this->Ca_i_new_ = __OLD_[16];
							this->Ca_ds_old_ = __OLD_AUX_[17];
							this->Ca_ds_new_ = __OLD_[17];
							this->Ca_up_old_ = __OLD_AUX_[18];
							this->Ca_up_new_ = __OLD_[18];
							this->Ca_rel_old_ = __OLD_AUX_[19];
							this->Ca_rel_new_ = __OLD_[19];
							this->Ca_Calmod_old_ = __OLD_AUX_[20];
							this->Ca_Calmod_new_ = __OLD_[20];
							this->Ca_Trop_old_ = __OLD_AUX_[21];
							this->Ca_Trop_new_ = __OLD_[21];
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
						__NEW_[14] = __K2_[14] * _prvt_dtime + __OLD_[14];
						__NEW_[15] = __K2_[15] * _prvt_dtime + __OLD_[15];
						__NEW_[16] = __K2_[16] * _prvt_dtime + __OLD_[16];
						__NEW_[17] = __K2_[17] * _prvt_dtime + __OLD_[17];
						__NEW_[18] = __K2_[18] * _prvt_dtime + __OLD_[18];
						__NEW_[19] = __K2_[19] * _prvt_dtime + __OLD_[19];
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
						__NEW_[10] = __K2_[10] * _prvt_dtime + __OLD_[10];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[12] = __K2_[12] * _prvt_dtime + __OLD_[12];
						__NEW_[13] = __K2_[13] * _prvt_dtime + __OLD_[13];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[8] = __K2_[8] * _prvt_dtime + __OLD_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[11])
					{
						__NEW_[9] = __K2_[9] * _prvt_dtime + __OLD_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[12])
					{
						__NEW_[11] = __K2_[11] * _prvt_dtime + __OLD_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[13])
					{
						__NEW_[20] = __K2_[20] * _prvt_dtime + __OLD_[20];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[14])
					{
						__NEW_[21] = __K2_[21] * _prvt_dtime + __OLD_[21];
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
		case 1:		return xr1_old_;    break;
		case 2:		return xr2_old_;    break;
		case 3:		return xs_old_;    break;
		case 4:		return m_old_;    break;
		case 5:		return h_old_;    break;
		case 6:		return d_old_;    break;
		case 7:		return f_old_;    break;
		case 8:		return f2_old_;    break;
		case 9:		return f2ds_old_;    break;
		case 10:		return s_old_;    break;
		case 11:		return r_old_;    break;
		case 12:		return ActFrac_old_;    break;
		case 13:		return ProdFrac_old_;    break;
		case 14:		return Na_i_old_;    break;
		case 15:		return K_i_old_;    break;
		case 16:		return Ca_i_old_;    break;
		case 17:		return Ca_ds_old_;    break;
		case 18:		return Ca_up_old_;    break;
		case 19:		return Ca_rel_old_;    break;
		case 20:		return Ca_Calmod_old_;    break;
		case 21:		return Ca_Trop_old_;    break;
		default:	return 1;    break;
		}
	}

	double Solveode::getLadoDireito(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return V_lado_direito_;    break;
		case 1:		return xr1_lado_direito_;    break;
		case 2:		return xr2_lado_direito_;    break;
		case 3:		return xs_lado_direito_;    break;
		case 4:		return m_lado_direito_;    break;
		case 5:		return h_lado_direito_;    break;
		case 6:		return d_lado_direito_;    break;
		case 7:		return f_lado_direito_;    break;
		case 8:		return f2_lado_direito_;    break;
		case 9:		return f2ds_lado_direito_;    break;
		case 10:		return s_lado_direito_;    break;
		case 11:		return r_lado_direito_;    break;
		case 12:		return ActFrac_lado_direito_;    break;
		case 13:		return ProdFrac_lado_direito_;    break;
		case 14:		return Na_i_lado_direito_;    break;
		case 15:		return K_i_lado_direito_;    break;
		case 16:		return Ca_i_lado_direito_;    break;
		case 17:		return Ca_ds_lado_direito_;    break;
		case 18:		return Ca_up_lado_direito_;    break;
		case 19:		return Ca_rel_lado_direito_;    break;
		case 20:		return Ca_Calmod_lado_direito_;    break;
		case 21:		return Ca_Trop_lado_direito_;    break;
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
		case 6:		return Cm;    break;
		case 7:		return R;    break;
		case 8:		return T;    break;
		case 9:		return F;    break;
		case 10:		return Na_o;    break;
		case 11:		return K_o;    break;
		case 12:		return P_kna;    break;
		case 13:		return Ca_o;    break;
		case 14:		return g_K1;    break;
		case 15:		return K_mk1;    break;
		case 16:		return g_Kr1;    break;
		case 17:		return g_Kr2;    break;
		case 18:		return g_Ks;    break;
		case 19:		return g_Na;    break;
		case 20:		return delta_m;    break;
		case 21:		return shift_h;    break;
		case 22:		return g_pna;    break;
		case 23:		return g_bna;    break;
		case 24:		return FrICa;    break;
		case 25:		return P_Ca_L;    break;
		case 26:		return P_CaK;    break;
		case 27:		return P_CaNa;    break;
		case 28:		return speed_d;    break;
		case 29:		return delta_f;    break;
		case 30:		return speed_f;    break;
		case 31:		return Km_f2;    break;
		case 32:		return R_decay;    break;
		case 33:		return Km_f2ds;    break;
		case 34:		return g_bca;    break;
		case 35:		return g_to;    break;
		case 36:		return g_tos;    break;
		case 37:		return i_NaK_max;    break;
		case 38:		return K_mK;    break;
		case 39:		return K_mNa;    break;
		case 40:		return FRiNaCa;    break;
		case 41:		return k_NaCa;    break;
		case 42:		return gamma;    break;
		case 43:		return n_NaCa;    break;
		case 44:		return d_NaCa;    break;
		case 45:		return K_cyca;    break;
		case 46:		return K_xcs;    break;
		case 47:		return K_srca;    break;
		case 48:		return alpha_up;    break;
		case 49:		return beta_up;    break;
		case 50:		return K_m_Ca_cyt;    break;
		case 51:		return K_m_Ca_ds;    break;
		case 52:		return K_m_rel;    break;
		case 53:		return K_leak_rate;    break;
		case 54:		return radius;    break;
		case 55:		return length;    break;
		case 56:		return V_e_ratio;    break;
		case 57:		return V_up_ratio;    break;
		case 58:		return V_rel_ratio;    break;
		case 59:		return V_ds_ratio;    break;
		case 60:		return Kdecay;    break;
		case 61:		return alpha_Calmod;    break;
		case 62:		return Calmod;    break;
		case 63:		return beta_Calmod;    break;
		case 64:		return alpha_Trop;    break;
		case 65:		return Trop;    break;
		case 66:		return beta_Trop;    break;
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
		Variables v("|V#|xr1#|xr2#|xs#|m#|h#|d#|f#|f2#|f2ds#|s#|r#|ActFrac#|ProdFrac#|Na_i#|K_i#|Ca_i#|Ca_ds#|Ca_up#|Ca_rel#|Ca_Calmod#|Ca_Trop#");
		return v;
	}
	Variables Solveode::get_Parameters()
	{
		Variables v("|time#|stim_amplitude#|stim_start#|stim_end#|stim_period#|stim_duration#|Cm#|R#|T#|F#|Na_o#|K_o#|P_kna#|Ca_o#|g_K1#|K_mk1#|g_Kr1#|g_Kr2#|g_Ks#|g_Na#|delta_m#|shift_h#|g_pna#|g_bna#|FrICa#|P_Ca_L#|P_CaK#|P_CaNa#|speed_d#|delta_f#|speed_f#|Km_f2#|R_decay#|Km_f2ds#|g_bca#|g_to#|g_tos#|i_NaK_max#|K_mK#|K_mNa#|FRiNaCa#|k_NaCa#|gamma#|n_NaCa#|d_NaCa#|K_cyca#|K_xcs#|K_srca#|alpha_up#|beta_up#|K_m_Ca_cyt#|K_m_Ca_ds#|K_m_rel#|K_leak_rate#|radius#|length#|V_e_ratio#|V_up_ratio#|V_rel_ratio#|V_ds_ratio#|Kdecay#|alpha_Calmod#|Calmod#|beta_Calmod#|alpha_Trop#|Trop#|beta_Trop#");
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
		xr1_old_ = xr1_ini_;
		if(xr1 != NULL)free( xr1);
			xr1 = (double *)malloc(sizeof(double)*num_results__);
		xr2_old_ = xr2_ini_;
		if(xr2 != NULL)free( xr2);
			xr2 = (double *)malloc(sizeof(double)*num_results__);
		xs_old_ = xs_ini_;
		if(xs != NULL)free( xs);
			xs = (double *)malloc(sizeof(double)*num_results__);
		m_old_ = m_ini_;
		if(m != NULL)free( m);
			m = (double *)malloc(sizeof(double)*num_results__);
		h_old_ = h_ini_;
		if(h != NULL)free( h);
			h = (double *)malloc(sizeof(double)*num_results__);
		d_old_ = d_ini_;
		if(d != NULL)free( d);
			d = (double *)malloc(sizeof(double)*num_results__);
		f_old_ = f_ini_;
		if(f != NULL)free( f);
			f = (double *)malloc(sizeof(double)*num_results__);
		f2_old_ = f2_ini_;
		if(f2 != NULL)free( f2);
			f2 = (double *)malloc(sizeof(double)*num_results__);
		f2ds_old_ = f2ds_ini_;
		if(f2ds != NULL)free( f2ds);
			f2ds = (double *)malloc(sizeof(double)*num_results__);
		s_old_ = s_ini_;
		if(s != NULL)free( s);
			s = (double *)malloc(sizeof(double)*num_results__);
		r_old_ = r_ini_;
		if(r != NULL)free( r);
			r = (double *)malloc(sizeof(double)*num_results__);
		ActFrac_old_ = ActFrac_ini_;
		if(ActFrac != NULL)free( ActFrac);
			ActFrac = (double *)malloc(sizeof(double)*num_results__);
		ProdFrac_old_ = ProdFrac_ini_;
		if(ProdFrac != NULL)free( ProdFrac);
			ProdFrac = (double *)malloc(sizeof(double)*num_results__);
		Na_i_old_ = Na_i_ini_;
		if(Na_i != NULL)free( Na_i);
			Na_i = (double *)malloc(sizeof(double)*num_results__);
		K_i_old_ = K_i_ini_;
		if(K_i != NULL)free( K_i);
			K_i = (double *)malloc(sizeof(double)*num_results__);
		Ca_i_old_ = Ca_i_ini_;
		if(Ca_i != NULL)free( Ca_i);
			Ca_i = (double *)malloc(sizeof(double)*num_results__);
		Ca_ds_old_ = Ca_ds_ini_;
		if(Ca_ds != NULL)free( Ca_ds);
			Ca_ds = (double *)malloc(sizeof(double)*num_results__);
		Ca_up_old_ = Ca_up_ini_;
		if(Ca_up != NULL)free( Ca_up);
			Ca_up = (double *)malloc(sizeof(double)*num_results__);
		Ca_rel_old_ = Ca_rel_ini_;
		if(Ca_rel != NULL)free( Ca_rel);
			Ca_rel = (double *)malloc(sizeof(double)*num_results__);
		Ca_Calmod_old_ = Ca_Calmod_ini_;
		if(Ca_Calmod != NULL)free( Ca_Calmod);
			Ca_Calmod = (double *)malloc(sizeof(double)*num_results__);
		Ca_Trop_old_ = Ca_Trop_ini_;
		if(Ca_Trop != NULL)free( Ca_Trop);
			Ca_Trop = (double *)malloc(sizeof(double)*num_results__);
		this->timeSaving = dtime;

		double diff=0;
		int counter=0;
		for (int i = 0; i< iterations;i++ )
		{
			this->time_new += dtime;

			rightHandSideFunction.function(this);
			this->V_new_ = this->V_old_ + this->V_lado_direito_ * this->dtime;
			this->xr1_new_ = this->xr1_old_ + this->xr1_lado_direito_ * this->dtime;
			this->xr2_new_ = this->xr2_old_ + this->xr2_lado_direito_ * this->dtime;
			this->xs_new_ = this->xs_old_ + this->xs_lado_direito_ * this->dtime;
			this->m_new_ = this->m_old_ + this->m_lado_direito_ * this->dtime;
			this->h_new_ = this->h_old_ + this->h_lado_direito_ * this->dtime;
			this->d_new_ = this->d_old_ + this->d_lado_direito_ * this->dtime;
			this->f_new_ = this->f_old_ + this->f_lado_direito_ * this->dtime;
			this->f2_new_ = this->f2_old_ + this->f2_lado_direito_ * this->dtime;
			this->f2ds_new_ = this->f2ds_old_ + this->f2ds_lado_direito_ * this->dtime;
			this->s_new_ = this->s_old_ + this->s_lado_direito_ * this->dtime;
			this->r_new_ = this->r_old_ + this->r_lado_direito_ * this->dtime;
			this->ActFrac_new_ = this->ActFrac_old_ + this->ActFrac_lado_direito_ * this->dtime;
			this->ProdFrac_new_ = this->ProdFrac_old_ + this->ProdFrac_lado_direito_ * this->dtime;
			this->Na_i_new_ = this->Na_i_old_ + this->Na_i_lado_direito_ * this->dtime;
			this->K_i_new_ = this->K_i_old_ + this->K_i_lado_direito_ * this->dtime;
			this->Ca_i_new_ = this->Ca_i_old_ + this->Ca_i_lado_direito_ * this->dtime;
			this->Ca_ds_new_ = this->Ca_ds_old_ + this->Ca_ds_lado_direito_ * this->dtime;
			this->Ca_up_new_ = this->Ca_up_old_ + this->Ca_up_lado_direito_ * this->dtime;
			this->Ca_rel_new_ = this->Ca_rel_old_ + this->Ca_rel_lado_direito_ * this->dtime;
			this->Ca_Calmod_new_ = this->Ca_Calmod_old_ + this->Ca_Calmod_lado_direito_ * this->dtime;
			this->Ca_Trop_new_ = this->Ca_Trop_old_ + this->Ca_Trop_lado_direito_ * this->dtime;
			diff =  _agos_round(this->time_new - timeSaving, 5);
			if(diff==0){
				this->timeSaving += svRate;
				time_vec__[counter] = this->time_new;
				V[counter] = this->V_new_;
				xr1[counter] = this->xr1_new_;
				xr2[counter] = this->xr2_new_;
				xs[counter] = this->xs_new_;
				m[counter] = this->m_new_;
				h[counter] = this->h_new_;
				d[counter] = this->d_new_;
				f[counter] = this->f_new_;
				f2[counter] = this->f2_new_;
				f2ds[counter] = this->f2ds_new_;
				s[counter] = this->s_new_;
				r[counter] = this->r_new_;
				ActFrac[counter] = this->ActFrac_new_;
				ProdFrac[counter] = this->ProdFrac_new_;
				Na_i[counter] = this->Na_i_new_;
				K_i[counter] = this->K_i_new_;
				Ca_i[counter] = this->Ca_i_new_;
				Ca_ds[counter] = this->Ca_ds_new_;
				Ca_up[counter] = this->Ca_up_new_;
				Ca_rel[counter] = this->Ca_rel_new_;
				Ca_Calmod[counter] = this->Ca_Calmod_new_;
				Ca_Trop[counter] = this->Ca_Trop_new_;
				counter++;
			}
			this->V_old_ = this->V_new_;
			this->xr1_old_ = this->xr1_new_;
			this->xr2_old_ = this->xr2_new_;
			this->xs_old_ = this->xs_new_;
			this->m_old_ = this->m_new_;
			this->h_old_ = this->h_new_;
			this->d_old_ = this->d_new_;
			this->f_old_ = this->f_new_;
			this->f2_old_ = this->f2_new_;
			this->f2ds_old_ = this->f2ds_new_;
			this->s_old_ = this->s_new_;
			this->r_old_ = this->r_new_;
			this->ActFrac_old_ = this->ActFrac_new_;
			this->ProdFrac_old_ = this->ProdFrac_new_;
			this->Na_i_old_ = this->Na_i_new_;
			this->K_i_old_ = this->K_i_new_;
			this->Ca_i_old_ = this->Ca_i_new_;
			this->Ca_ds_old_ = this->Ca_ds_new_;
			this->Ca_up_old_ = this->Ca_up_new_;
			this->Ca_rel_old_ = this->Ca_rel_new_;
			this->Ca_Calmod_old_ = this->Ca_Calmod_new_;
			this->Ca_Trop_old_ = this->Ca_Trop_new_;
		}
		double h_jac_[numEDO];
		double quociente = 1000.0;
		h_jac_[0] = fabs(_agos_min(V, num_results__) / _agos_max(V, num_results__) );
		h_jac_[1] = fabs(_agos_min(xr1, num_results__) / _agos_max(xr1, num_results__) );
		h_jac_[2] = fabs(_agos_min(xr2, num_results__) / _agos_max(xr2, num_results__) );
		h_jac_[3] = fabs(_agos_min(xs, num_results__) / _agos_max(xs, num_results__) );
		h_jac_[4] = fabs(_agos_min(m, num_results__) / _agos_max(m, num_results__) );
		h_jac_[5] = fabs(_agos_min(h, num_results__) / _agos_max(h, num_results__) );
		h_jac_[6] = fabs(_agos_min(d, num_results__) / _agos_max(d, num_results__) );
		h_jac_[7] = fabs(_agos_min(f, num_results__) / _agos_max(f, num_results__) );
		h_jac_[8] = fabs(_agos_min(f2, num_results__) / _agos_max(f2, num_results__) );
		h_jac_[9] = fabs(_agos_min(f2ds, num_results__) / _agos_max(f2ds, num_results__) );
		h_jac_[10] = fabs(_agos_min(s, num_results__) / _agos_max(s, num_results__) );
		h_jac_[11] = fabs(_agos_min(r, num_results__) / _agos_max(r, num_results__) );
		h_jac_[12] = fabs(_agos_min(ActFrac, num_results__) / _agos_max(ActFrac, num_results__) );
		h_jac_[13] = fabs(_agos_min(ProdFrac, num_results__) / _agos_max(ProdFrac, num_results__) );
		h_jac_[14] = fabs(_agos_min(Na_i, num_results__) / _agos_max(Na_i, num_results__) );
		h_jac_[15] = fabs(_agos_min(K_i, num_results__) / _agos_max(K_i, num_results__) );
		h_jac_[16] = fabs(_agos_min(Ca_i, num_results__) / _agos_max(Ca_i, num_results__) );
		h_jac_[17] = fabs(_agos_min(Ca_ds, num_results__) / _agos_max(Ca_ds, num_results__) );
		h_jac_[18] = fabs(_agos_min(Ca_up, num_results__) / _agos_max(Ca_up, num_results__) );
		h_jac_[19] = fabs(_agos_min(Ca_rel, num_results__) / _agos_max(Ca_rel, num_results__) );
		h_jac_[20] = fabs(_agos_min(Ca_Calmod, num_results__) / _agos_max(Ca_Calmod, num_results__) );
		h_jac_[21] = fabs(_agos_min(Ca_Trop, num_results__) / _agos_max(Ca_Trop, num_results__) );
		for(int l=0;l<numEDO;l++){
			h_jac_[l] = (h_jac_[l]==0 || h_jac_[l]==AGOS_NAN || h_jac_[l]==AGOS_INF)?this->dtime:h_jac_[l];
		}
		this->timeSaving = this->dtime;

		this->time_new = this->time;

		V_old_ = V_ini_;
		xr1_old_ = xr1_ini_;
		xr2_old_ = xr2_ini_;
		xs_old_ = xs_ini_;
		m_old_ = m_ini_;
		h_old_ = h_ini_;
		d_old_ = d_ini_;
		f_old_ = f_ini_;
		f2_old_ = f2_ini_;
		f2ds_old_ = f2ds_ini_;
		s_old_ = s_ini_;
		r_old_ = r_ini_;
		ActFrac_old_ = ActFrac_ini_;
		ProdFrac_old_ = ProdFrac_ini_;
		Na_i_old_ = Na_i_ini_;
		K_i_old_ = K_i_ini_;
		Ca_i_old_ = Ca_i_ini_;
		Ca_ds_old_ = Ca_ds_ini_;
		Ca_up_old_ = Ca_up_ini_;
		Ca_rel_old_ = Ca_rel_ini_;
		Ca_Calmod_old_ = Ca_Calmod_ini_;
		Ca_Trop_old_ = Ca_Trop_ini_;
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
			this->xr1_new_ = this->dtime*(this->xr1_lado_direito_) + this->xr1_old_;
			this->xr2_new_ = this->dtime*(this->xr2_lado_direito_) + this->xr2_old_;
			this->xs_new_ = this->dtime*(this->xs_lado_direito_) + this->xs_old_;
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h_new_ = this->dtime*(this->h_lado_direito_) + this->h_old_;
			this->d_new_ = this->dtime*(this->d_lado_direito_) + this->d_old_;
			this->f_new_ = this->dtime*(this->f_lado_direito_) + this->f_old_;
			this->f2_new_ = this->dtime*(this->f2_lado_direito_) + this->f2_old_;
			this->f2ds_new_ = this->dtime*(this->f2ds_lado_direito_) + this->f2ds_old_;
			this->s_new_ = this->dtime*(this->s_lado_direito_) + this->s_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->ActFrac_new_ = this->dtime*(this->ActFrac_lado_direito_) + this->ActFrac_old_;
			this->ProdFrac_new_ = this->dtime*(this->ProdFrac_lado_direito_) + this->ProdFrac_old_;
			this->Na_i_new_ = this->dtime*(this->Na_i_lado_direito_) + this->Na_i_old_;
			this->K_i_new_ = this->dtime*(this->K_i_lado_direito_) + this->K_i_old_;
			this->Ca_i_new_ = this->dtime*(this->Ca_i_lado_direito_) + this->Ca_i_old_;
			this->Ca_ds_new_ = this->dtime*(this->Ca_ds_lado_direito_) + this->Ca_ds_old_;
			this->Ca_up_new_ = this->dtime*(this->Ca_up_lado_direito_) + this->Ca_up_old_;
			this->Ca_rel_new_ = this->dtime*(this->Ca_rel_lado_direito_) + this->Ca_rel_old_;
			this->Ca_Calmod_new_ = this->dtime*(this->Ca_Calmod_lado_direito_) + this->Ca_Calmod_old_;
			this->Ca_Trop_new_ = this->dtime*(this->Ca_Trop_lado_direito_) + this->Ca_Trop_old_;
			diff =  _agos_round(this->time_new - timeSaving, 5);
			if(diff==0){
				this->timeSaving += svRate;
				V[counter2] = V_new_;
				xr1[counter2] = xr1_new_;
				xr2[counter2] = xr2_new_;
				xs[counter2] = xs_new_;
				m[counter2] = m_new_;
				h[counter2] = h_new_;
				d[counter2] = d_new_;
				f[counter2] = f_new_;
				f2[counter2] = f2_new_;
				f2ds[counter2] = f2ds_new_;
				s[counter2] = s_new_;
				r[counter2] = r_new_;
				ActFrac[counter2] = ActFrac_new_;
				ProdFrac[counter2] = ProdFrac_new_;
				Na_i[counter2] = Na_i_new_;
				K_i[counter2] = K_i_new_;
				Ca_i[counter2] = Ca_i_new_;
				Ca_ds[counter2] = Ca_ds_new_;
				Ca_up[counter2] = Ca_up_new_;
				Ca_rel[counter2] = Ca_rel_new_;
				Ca_Calmod[counter2] = Ca_Calmod_new_;
				Ca_Trop[counter2] = Ca_Trop_new_;
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
			xr1_old_ = xr1_new_;
			xr2_old_ = xr2_new_;
			xs_old_ = xs_new_;
			m_old_ = m_new_;
			h_old_ = h_new_;
			d_old_ = d_new_;
			f_old_ = f_new_;
			f2_old_ = f2_new_;
			f2ds_old_ = f2ds_new_;
			s_old_ = s_new_;
			r_old_ = r_new_;
			ActFrac_old_ = ActFrac_new_;
			ProdFrac_old_ = ProdFrac_new_;
			Na_i_old_ = Na_i_new_;
			K_i_old_ = K_i_new_;
			Ca_i_old_ = Ca_i_new_;
			Ca_ds_old_ = Ca_ds_new_;
			Ca_up_old_ = Ca_up_new_;
			Ca_rel_old_ = Ca_rel_new_;
			Ca_Calmod_old_ = Ca_Calmod_new_;
			Ca_Trop_old_ = Ca_Trop_new_;
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
		case 1:		return xr1;    break;
		case 2:		return xr2;    break;
		case 3:		return xs;    break;
		case 4:		return m;    break;
		case 5:		return h;    break;
		case 6:		return d;    break;
		case 7:		return f;    break;
		case 8:		return f2;    break;
		case 9:		return f2ds;    break;
		case 10:		return s;    break;
		case 11:		return r;    break;
		case 12:		return ActFrac;    break;
		case 13:		return ProdFrac;    break;
		case 14:		return Na_i;    break;
		case 15:		return K_i;    break;
		case 16:		return Ca_i;    break;
		case 17:		return Ca_ds;    break;
		case 18:		return Ca_up;    break;
		case 19:		return Ca_rel;    break;
		case 20:		return Ca_Calmod;    break;
		case 21:		return Ca_Trop;    break;
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
		if((fabs(calc_E0_m)<delta_m)){
			return (2.0000000000e+03);
		}else{
			return (((2.0000000000e+02*calc_E0_m)/(1.0000000000e+00-exp(((-1.0000000000e-01)*calc_E0_m)))));
		}
	}
	double Solveode::ifnumber_2(){
		if((fabs(calc_E0_d)<1.0000000000e-04)){
			return (1.2000000000e+02);
		}else{
			return (((3.0000000000e+01*calc_E0_d)/(1.0000000000e+00-exp(((-calc_E0_d)/4.0000000000e+00)))));
		}
	}
	double Solveode::ifnumber_3(){
		if((fabs(calc_E0_d)<1.0000000000e-04)){
			return (1.2000000000e+02);
		}else{
			return (((1.2000000000e+01*calc_E0_d)/(exp((calc_E0_d/1.0000000000e+01))-1.0000000000e+00)));
		}
	}
	double Solveode::ifnumber_4(){
		if((fabs(calc_E0_f)<delta_f)){
			return (2.5000000000e+01);
		}else{
			return (((6.2500000000e+00*calc_E0_f)/(exp((calc_E0_f/4.0000000000e+00))-1.0000000000e+00)));
		}
	}
	double Solveode::ifnumber_5(){
		if((V_old_<(-5.0000000000e+01))){
			return (5.0000000000e+00);
		}else{
			return (1.0000000000e+00);
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
		double dtL, dtM, dtMax=0.0,  dtMin=990.0 ;
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
			dependent_variable__ = N_VNew_Serial(22);
			if(check_flag((void *)dependent_variable__, "N_VNew_Serial", 0))
			exit(1);
			depvar__ = (double*)malloc(sizeof(double)*22);
			if(depvar__ == NULL){
			fprintf(stderr, "ERROR Cannot allocate memory for depvar__\n");
			exit(0);
			}
			NV_Ith_S(dependent_variable__, 0) = V_ini_;
			NV_Ith_S(dependent_variable__, 1) = xr1_ini_;
			NV_Ith_S(dependent_variable__, 2) = xr2_ini_;
			NV_Ith_S(dependent_variable__, 3) = xs_ini_;
			NV_Ith_S(dependent_variable__, 4) = m_ini_;
			NV_Ith_S(dependent_variable__, 5) = h_ini_;
			NV_Ith_S(dependent_variable__, 6) = d_ini_;
			NV_Ith_S(dependent_variable__, 7) = f_ini_;
			NV_Ith_S(dependent_variable__, 8) = f2_ini_;
			NV_Ith_S(dependent_variable__, 9) = f2ds_ini_;
			NV_Ith_S(dependent_variable__, 10) = s_ini_;
			NV_Ith_S(dependent_variable__, 11) = r_ini_;
			NV_Ith_S(dependent_variable__, 12) = ActFrac_ini_;
			NV_Ith_S(dependent_variable__, 13) = ProdFrac_ini_;
			NV_Ith_S(dependent_variable__, 14) = Na_i_ini_;
			NV_Ith_S(dependent_variable__, 15) = K_i_ini_;
			NV_Ith_S(dependent_variable__, 16) = Ca_i_ini_;
			NV_Ith_S(dependent_variable__, 17) = Ca_ds_ini_;
			NV_Ith_S(dependent_variable__, 18) = Ca_up_ini_;
			NV_Ith_S(dependent_variable__, 19) = Ca_rel_ini_;
			NV_Ith_S(dependent_variable__, 20) = Ca_Calmod_ini_;
			NV_Ith_S(dependent_variable__, 21) = Ca_Trop_ini_;
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
			flag__ = CVDense(cvode_mem_cvode__, 22);
			if (check_flag(&flag__, "CVDense", 1))	exit(1);
			break;
			case 2:
			flag__ = CVDiag(cvode_mem_cvode__);
			if (check_flag(&flag__, "CVDiag", 1))	exit(1);
			break;
			case 3:
			flag__ = CVBand(cvode_mem_cvode__, 22, NULL, NULL);
			if (check_flag(&flag__, "CVBand", 1))	exit(1);
			break;
			}
			CVodeSetFdata(cvode_mem_cvode__, (void*)this);
			V_old_ = V_ini_;
			xr1_old_ = xr1_ini_;
			xr2_old_ = xr2_ini_;
			xs_old_ = xs_ini_;
			m_old_ = m_ini_;
			h_old_ = h_ini_;
			d_old_ = d_ini_;
			f_old_ = f_ini_;
			f2_old_ = f2_ini_;
			f2ds_old_ = f2ds_ini_;
			s_old_ = s_ini_;
			r_old_ = r_ini_;
			ActFrac_old_ = ActFrac_ini_;
			ProdFrac_old_ = ProdFrac_ini_;
			Na_i_old_ = Na_i_ini_;
			K_i_old_ = K_i_ini_;
			Ca_i_old_ = Ca_i_ini_;
			Ca_ds_old_ = Ca_ds_ini_;
			Ca_up_old_ = Ca_up_ini_;
			Ca_rel_old_ = Ca_rel_ini_;
			Ca_Calmod_old_ = Ca_Calmod_ini_;
			Ca_Trop_old_ = Ca_Trop_ini_;
		}
		while(1){
			flag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &time_new, CV_NORMAL);
			V_old_ = NV_Ith_S(dependent_variable__, 0);
			xr1_old_ = NV_Ith_S(dependent_variable__, 1);
			xr2_old_ = NV_Ith_S(dependent_variable__, 2);
			xs_old_ = NV_Ith_S(dependent_variable__, 3);
			m_old_ = NV_Ith_S(dependent_variable__, 4);
			h_old_ = NV_Ith_S(dependent_variable__, 5);
			d_old_ = NV_Ith_S(dependent_variable__, 6);
			f_old_ = NV_Ith_S(dependent_variable__, 7);
			f2_old_ = NV_Ith_S(dependent_variable__, 8);
			f2ds_old_ = NV_Ith_S(dependent_variable__, 9);
			s_old_ = NV_Ith_S(dependent_variable__, 10);
			r_old_ = NV_Ith_S(dependent_variable__, 11);
			ActFrac_old_ = NV_Ith_S(dependent_variable__, 12);
			ProdFrac_old_ = NV_Ith_S(dependent_variable__, 13);
			Na_i_old_ = NV_Ith_S(dependent_variable__, 14);
			K_i_old_ = NV_Ith_S(dependent_variable__, 15);
			Ca_i_old_ = NV_Ith_S(dependent_variable__, 16);
			Ca_ds_old_ = NV_Ith_S(dependent_variable__, 17);
			Ca_up_old_ = NV_Ith_S(dependent_variable__, 18);
			Ca_rel_old_ = NV_Ith_S(dependent_variable__, 19);
			Ca_Calmod_old_ = NV_Ith_S(dependent_variable__, 20);
			Ca_Trop_old_ = NV_Ith_S(dependent_variable__, 21);
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
	for(int i = 0; i<22; i++)
		ode->setVariables( i ,NV_Ith_S(dependent_variable__, i));
	ode->setParameters(0,time);
	rightHandSideFunction.function(ode);
	for(int i = 0; i<22; i++)
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
