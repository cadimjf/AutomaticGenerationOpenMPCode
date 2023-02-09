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
#define forestSize 11
#include <sys/time.h>
#include "greedy.cpp"
#define _INCREASE_DT_ 1.5
#define _DECREASE_DT_ 0.65
#define _DECREASE_DT_2_ 2
#define numEDO 15
#define numAux 75

class Solveode
{
public:
	int it_countx;
	int *tree_thread;
	double time; 	 // second
	double dCell; 	 // dimensionless
	double Version; 	 // dimensionless
	double FCellConstant; 	 // dimensionless
	double CmCentre; 	 // microF
	double CmPeriphery; 	 // microF
	double g_Na_Centre_Published; 	 // microlitre_per_second
	double g_Na_Periphery_Published; 	 // microlitre_per_second
	double g_Na_Centre_0DCapable; 	 // microlitre_per_second
	double g_Na_Periphery_0DCapable; 	 // microlitre_per_second
	double g_Na_Centre_1DCapable; 	 // microlitre_per_second
	double g_Na_Periphery_1DCapable; 	 // microlitre_per_second
	double Na_o; 	 // millimolar
	double F; 	 // coulomb_per_mole
	double R; 	 // millijoule_per_mole_kelvin
	double T; 	 // kelvin
	double g_Ca_L_Centre_Published; 	 // microS
	double g_Ca_L_Periphery_Published; 	 // microS
	double g_Ca_L_Centre_0DCapable; 	 // microS
	double g_Ca_L_Periphery_0DCapable; 	 // microS
	double g_Ca_L_Centre_1DCapable; 	 // microS
	double g_Ca_L_Periphery_1DCapable; 	 // microS
	double E_Ca_L; 	 // millivolt
	double g_Ca_T_Centre_Published; 	 // microS
	double g_Ca_T_Periphery_Published; 	 // microS
	double g_Ca_T_Centre_0DCapable; 	 // microS
	double g_Ca_T_Periphery_0DCapable; 	 // microS
	double g_Ca_T_Centre_1DCapable; 	 // microS
	double g_Ca_T_Periphery_1DCapable; 	 // microS
	double E_Ca_T; 	 // millivolt
	double g_to_Centre_Published; 	 // microS
	double g_to_Periphery_Published; 	 // microS
	double g_to_Centre_0DCapable; 	 // microS
	double g_to_Periphery_0DCapable; 	 // microS
	double g_to_Centre_1DCapable; 	 // microS
	double g_to_Periphery_1DCapable; 	 // microS
	double g_sus_Centre_Published; 	 // microS
	double g_sus_Periphery_Published; 	 // microS
	double g_sus_Centre_0DCapable; 	 // microS
	double g_sus_Periphery_0DCapable; 	 // microS
	double g_sus_Centre_1DCapable; 	 // microS
	double g_sus_Periphery_1DCapable; 	 // microS
	double g_K_r_Centre_Published; 	 // microS
	double g_K_r_Periphery_Published; 	 // microS
	double g_K_r_Centre_0DCapable; 	 // microS
	double g_K_r_Periphery_0DCapable; 	 // microS
	double g_K_r_Centre_1DCapable; 	 // microS
	double g_K_r_Periphery_1DCapable; 	 // microS
	double g_K_s_Centre_Published; 	 // microS
	double g_K_s_Periphery_Published; 	 // microS
	double g_K_s_Centre_0DCapable; 	 // microS
	double g_K_s_Periphery_0DCapable; 	 // microS
	double g_K_s_Centre_1DCapable; 	 // microS
	double g_K_s_Periphery_1DCapable; 	 // microS
	double g_f_Na_Centre_Published; 	 // microS
	double g_f_Na_Periphery_Published; 	 // microS
	double g_f_Na_Centre_0DCapable; 	 // microS
	double g_f_Na_Periphery_0DCapable; 	 // microS
	double g_f_Na_Centre_1DCapable; 	 // microS
	double g_f_Na_Periphery_1DCapable; 	 // microS
	double g_f_K_Centre_Published; 	 // microS
	double g_f_K_Periphery_Published; 	 // microS
	double g_f_K_Centre_0DCapable; 	 // microS
	double g_f_K_Periphery_0DCapable; 	 // microS
	double g_f_K_Centre_1DCapable; 	 // microS
	double g_f_K_Periphery_1DCapable; 	 // microS
	double g_b_Na_Centre_Published; 	 // microS
	double g_b_Na_Periphery_Published; 	 // microS
	double g_b_Na_Centre_0DCapable; 	 // microS
	double g_b_Na_Periphery_0DCapable; 	 // microS
	double g_b_Na_Centre_1DCapable; 	 // microS
	double g_b_Na_Periphery_1DCapable; 	 // microS
	double g_b_K_Centre_Published; 	 // microS
	double g_b_K_Periphery_Published; 	 // microS
	double g_b_K_Centre_0DCapable; 	 // microS
	double g_b_K_Periphery_0DCapable; 	 // microS
	double g_b_K_Centre_1DCapable; 	 // microS
	double g_b_K_Periphery_1DCapable; 	 // microS
	double g_b_Ca_Centre_Published; 	 // microS
	double g_b_Ca_Periphery_Published; 	 // microS
	double g_b_Ca_Centre_0DCapable; 	 // microS
	double g_b_Ca_Periphery_0DCapable; 	 // microS
	double g_b_Ca_Centre_1DCapable; 	 // microS
	double g_b_Ca_Periphery_1DCapable; 	 // microS
	double k_NaCa_Centre_Published; 	 // nanoA
	double k_NaCa_Periphery_Published; 	 // nanoA
	double k_NaCa_Centre_0DCapable; 	 // nanoA
	double k_NaCa_Periphery_0DCapable; 	 // nanoA
	double k_NaCa_Centre_1DCapable; 	 // nanoA
	double k_NaCa_Periphery_1DCapable; 	 // nanoA
	double Na_i; 	 // millimolar
	double Ca_o; 	 // millimolar
	double gamma_NaCa; 	 // dimensionless
	double Ca_i; 	 // millimolar
	double d_NaCa; 	 // dimensionless
	double i_p_max_Centre_Published; 	 // nanoA
	double i_p_max_Periphery_Published; 	 // nanoA
	double i_p_max_Centre_0DCapable; 	 // nanoA
	double i_p_max_Periphery_0DCapable; 	 // nanoA
	double i_p_max_Centre_1DCapable; 	 // nanoA
	double i_p_max_Periphery_1DCapable; 	 // nanoA
	double K_m_Na; 	 // millimolar
	double K_o; 	 // millimolar
	double K_m_K; 	 // millimolar
	double i_Ca_p_max_Centre_Published; 	 // nanoA
	double i_Ca_p_max_Periphery_Published; 	 // nanoA
	double i_Ca_p_max_Centre_0DCapable; 	 // nanoA
	double i_Ca_p_max_Periphery_0DCapable; 	 // nanoA
	double i_Ca_p_max_Centre_1DCapable; 	 // nanoA
	double i_Ca_p_max_Periphery_1DCapable; 	 // nanoA
	double K_i; 	 // millimolar
	double calc_FCell; 	 // dimensionless
	double calc_Cm; 	 // microF
	double calc_g_Na; 	 // microlitre_per_second
	double calc_i_Na; 	 // nanoA
	double calc_m_infinity; 	 // dimensionless
	double calc_tau_m; 	 // second
	double calc_F_Na; 	 // dimensionless
	double calc_h; 	 // dimensionless
	double calc_h1_infinity; 	 // dimensionless
	double calc_h2_infinity; 	 // dimensionless
	double calc_tau_h1; 	 // second
	double calc_tau_h2; 	 // second
	double calc_g_Ca_L; 	 // microS
	double calc_i_Ca_L; 	 // nanoA
	double calc_alpha_d_L; 	 // per_second
	double calc_beta_d_L; 	 // per_second
	double calc_tau_d_L; 	 // second
	double calc_d_L_infinity; 	 // dimensionless
	double calc_alpha_f_L; 	 // per_second
	double calc_beta_f_L; 	 // per_second
	double calc_tau_f_L; 	 // second
	double calc_f_L_infinity; 	 // dimensionless
	double calc_g_Ca_T; 	 // microS
	double calc_i_Ca_T; 	 // nanoA
	double calc_alpha_d_T; 	 // per_second
	double calc_beta_d_T; 	 // per_second
	double calc_tau_d_T; 	 // second
	double calc_d_T_infinity; 	 // dimensionless
	double calc_alpha_f_T; 	 // per_second
	double calc_beta_f_T; 	 // per_second
	double calc_tau_f_T; 	 // second
	double calc_f_T_infinity; 	 // dimensionless
	double calc_g_to; 	 // microS
	double calc_g_sus; 	 // microS
	double calc_i_to; 	 // nanoA
	double calc_i_sus; 	 // nanoA
	double calc_q_infinity; 	 // dimensionless
	double calc_tau_q; 	 // second
	double calc_r_infinity; 	 // dimensionless
	double calc_tau_r; 	 // second
	double calc_g_K_r; 	 // microS
	double calc_i_K_r; 	 // nanoA
	double calc_P_a; 	 // dimensionless
	double calc_P_af_infinity; 	 // dimensionless
	double calc_tau_P_af; 	 // second
	double calc_P_as_infinity; 	 // dimensionless
	double calc_tau_P_as; 	 // second
	double calc_tau_P_i; 	 // second
	double calc_P_i_infinity; 	 // dimensionless
	double calc_g_K_s; 	 // microS
	double calc_i_K_s; 	 // nanoA
	double calc_alpha_xs; 	 // per_second
	double calc_beta_xs; 	 // per_second
	double calc_g_f_Na; 	 // microS
	double calc_i_f_Na; 	 // nanoA
	double calc_g_f_K; 	 // microS
	double calc_i_f_K; 	 // nanoA
	double calc_alpha_y; 	 // per_second
	double calc_beta_y; 	 // per_second
	double calc_g_b_Na; 	 // microS
	double calc_i_b_Na; 	 // nanoA
	double calc_g_b_K; 	 // microS
	double calc_i_b_K; 	 // nanoA
	double calc_g_b_Ca; 	 // microS
	double calc_i_b_Ca; 	 // nanoA
	double calc_k_NaCa; 	 // nanoA
	double calc_i_NaCa; 	 // nanoA
	double calc_i_p_max; 	 // nanoA
	double calc_i_p; 	 // nanoA
	double calc_i_Ca_p_max; 	 // nanoA
	double calc_i_Ca_p; 	 // nanoA
	double calc_E_Na; 	 // millivolt
	double calc_E_K; 	 // millivolt
	double calc_E_Ca; 	 // millivolt
	double calc_E_K_s; 	 // millivolt
	double dtime, *time_vec__;
	double time_new;

	//functions variables
	double *V;
	double V_new_, V_old_, V_ini_, V_lado_direito_;
	double *m;
	double m_new_, m_old_, m_ini_, m_lado_direito_;
	double *h1;
	double h1_new_, h1_old_, h1_ini_, h1_lado_direito_;
	double *h2;
	double h2_new_, h2_old_, h2_ini_, h2_lado_direito_;
	double *d_L;
	double d_L_new_, d_L_old_, d_L_ini_, d_L_lado_direito_;
	double *f_L;
	double f_L_new_, f_L_old_, f_L_ini_, f_L_lado_direito_;
	double *d_T;
	double d_T_new_, d_T_old_, d_T_ini_, d_T_lado_direito_;
	double *f_T;
	double f_T_new_, f_T_old_, f_T_ini_, f_T_lado_direito_;
	double *q;
	double q_new_, q_old_, q_ini_, q_lado_direito_;
	double *r;
	double r_new_, r_old_, r_ini_, r_lado_direito_;
	double *P_af;
	double P_af_new_, P_af_old_, P_af_ini_, P_af_lado_direito_;
	double *P_as;
	double P_as_new_, P_as_old_, P_as_ini_, P_as_lado_direito_;
	double *P_i;
	double P_i_new_, P_i_old_, P_i_ini_, P_i_lado_direito_;
	double *xs;
	double xs_new_, xs_old_, xs_ini_, xs_lado_direito_;
	double *y;
	double y_new_, y_old_, y_ini_, y_lado_direito_;

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
	inline double ifnumber_6();
	inline double ifnumber_7();
	inline double ifnumber_8();
	inline double ifnumber_9();
	inline double ifnumber_10();
	inline double ifnumber_11();
	inline double ifnumber_12();
	inline double ifnumber_13();
	inline double ifnumber_14();
	inline double ifnumber_15();
	inline double ifnumber_16();
	inline double ifnumber_17();
	inline double ifnumber_18();
	inline double ifnumber_19();
	inline double ifnumber_20();
	inline double ifnumber_21();
	inline double ifnumber_22();
	inline double ifnumber_23();
	inline double ifnumber_24();
	inline double ifnumber_25();
	inline double ifnumber_26();
	inline double ifnumber_27();
	inline double ifnumber_28();
	inline double ifnumber_29();
	inline double ifnumber_30();
	inline double ifnumber_31();
	inline double ifnumber_32();
	inline double ifnumber_33();
	inline double ifnumber_34();
	inline double ifnumber_35();
	inline double ifnumber_36();
	inline double ifnumber_37();
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
	__AGOS->calc_FCell = ((__AGOS->Version==0.0000000000e+00))
?(((1.0700000000e+00*((3.0000000000e+00*__AGOS->dCell)-1.0000000000e-01))/(3.0000000000e+00*(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*__AGOS->dCell)-2.0500000000e+00))/2.9500000000e-01)))))))
:(((__AGOS->Version==1.0000000000e+00))
?(((__AGOS->FCellConstant*__AGOS->dCell)/(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*__AGOS->dCell)-2.0500000000e+00))/2.9500000000e-01))))))
:(((1.0700000000e+00*2.9000000000e+01*__AGOS->dCell)/(3.0000000000e+01*(1.0000000000e+00+(7.7450000000e-01*exp(((-((2.9000000000e+01*__AGOS->dCell)-2.4500000000e+01))/1.9500000000e+00))))))))
;
	__AGOS->calc_F_Na = ((__AGOS->Version==0.0000000000e+00))
?((((9.5200000000e-02*exp(((-6.3000000000e-02)*(__AGOS->V_old_+3.4400000000e+01))))/(1.0000000000e+00+(1.6600000000e+00*exp(((-2.2500000000e-01)*(__AGOS->V_old_+6.3700000000e+01))))))+8.6900000000e-02))
:((((9.5180000000e-02*exp(((-6.3060000000e-02)*(__AGOS->V_old_+3.4400000000e+01))))/(1.0000000000e+00+(1.6620000000e+00*exp(((-2.2510000000e-01)*(__AGOS->V_old_+6.3700000000e+01))))))+8.6930000000e-02));
	__AGOS->calc_alpha_f_L = ((__AGOS->Version==1.0000000000e+00))
?(((3.7500000000e+00*(__AGOS->V_old_+2.8000000000e+01))/(exp(((__AGOS->V_old_+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)))
:(((3.1200000000e+00*(__AGOS->V_old_+2.8000000000e+01))/(exp(((__AGOS->V_old_+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
	__AGOS->calc_beta_f_L = ((__AGOS->Version==1.0000000000e+00))
?((3.0000000000e+01/(1.0000000000e+00+exp(((-(__AGOS->V_old_+2.8000000000e+01))/4.0000000000e+00)))))
:((2.5000000000e+01/(1.0000000000e+00+exp(((-(__AGOS->V_old_+2.8000000000e+01))/4.0000000000e+00)))));
	__AGOS->calc_f_L_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+4.5000000000e+01)/5.0000000000e+00))));
	__AGOS->calc_beta_f_T = ((__AGOS->Version==1.0000000000e+00))
?((1.5000000000e+01*exp(((__AGOS->V_old_+7.1000000000e+01)/1.5380000000e+01))))
:((1.5000000000e+01*exp(((__AGOS->V_old_+7.1700000000e+01)/1.5380000000e+01))));
	__AGOS->calc_f_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+7.1000000000e+01)/9.0000000000e+00))));
	__AGOS->calc_P_a = ((6.0000000000e-01*__AGOS->P_af_old_)+(4.0000000000e-01*__AGOS->P_as_old_));
	__AGOS->calc_E_Na = (((__AGOS->R*__AGOS->T)/__AGOS->F)*log((__AGOS->Na_o/__AGOS->Na_i)));
	__AGOS->calc_E_K = (((__AGOS->R*__AGOS->T)/__AGOS->F)*log((__AGOS->K_o/__AGOS->K_i)));
	__AGOS->calc_E_Ca = (((__AGOS->R*__AGOS->T)/(2.0000000000e+00*__AGOS->F))*log((__AGOS->Ca_o/__AGOS->Ca_i)));
	__AGOS->calc_E_K_s = ((__AGOS->Version==0.0000000000e+00))
?((((__AGOS->R*__AGOS->T)/__AGOS->F)*log(((__AGOS->K_o+(1.2000000000e-01*__AGOS->Na_o))/(__AGOS->K_i+(1.2000000000e-01*__AGOS->Na_i))))))
:((((__AGOS->R*__AGOS->T)/__AGOS->F)*log(((__AGOS->K_o+(3.0000000000e-02*__AGOS->Na_o))/(__AGOS->K_i+(3.0000000000e-02*__AGOS->Na_i))))));
	__AGOS->calc_Cm = (__AGOS->CmCentre+(__AGOS->calc_FCell*(__AGOS->CmPeriphery-__AGOS->CmCentre)));
	__AGOS->calc_g_Na = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_Na_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_Na_Periphery_Published-__AGOS->g_Na_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_Na_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_Na_Periphery_0DCapable-__AGOS->g_Na_Centre_0DCapable))))
:((__AGOS->g_Na_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_Na_Periphery_1DCapable-__AGOS->g_Na_Centre_1DCapable)))))
;
	__AGOS->calc_h = (((1.0000000000e+00-__AGOS->calc_F_Na)*__AGOS->h1_old_)+(__AGOS->calc_F_Na*__AGOS->h2_old_));
	__AGOS->calc_g_Ca_L = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_Ca_L_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_Ca_L_Periphery_Published-__AGOS->g_Ca_L_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_Ca_L_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_Ca_L_Periphery_0DCapable-__AGOS->g_Ca_L_Centre_0DCapable))))
:((__AGOS->g_Ca_L_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_Ca_L_Periphery_1DCapable-__AGOS->g_Ca_L_Centre_1DCapable)))))
;
	__AGOS->calc_tau_f_L = ((__AGOS->Version==1.0000000000e+00))
?(((1.2000000000e+00-(2.0000000000e-01*__AGOS->calc_FCell))/(__AGOS->calc_alpha_f_L+__AGOS->calc_beta_f_L)))
:((1.0000000000e+00/(__AGOS->calc_alpha_f_L+__AGOS->calc_beta_f_L)));
	__AGOS->calc_g_Ca_T = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_Ca_T_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_Ca_T_Periphery_Published-__AGOS->g_Ca_T_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_Ca_T_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_Ca_T_Periphery_0DCapable-__AGOS->g_Ca_T_Centre_0DCapable))))
:((__AGOS->g_Ca_T_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_Ca_T_Periphery_1DCapable-__AGOS->g_Ca_T_Centre_1DCapable)))))
;
	__AGOS->calc_alpha_f_T = ((__AGOS->Version==1.0000000000e+00))
?((1.5300000000e+01*exp(((-(__AGOS->V_old_+7.1000000000e+01+(7.0000000000e-01*__AGOS->calc_FCell)))/8.3300000000e+01))))
:((1.5300000000e+01*exp(((-(__AGOS->V_old_+7.1700000000e+01))/8.3300000000e+01))));
	__AGOS->calc_g_to = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_to_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_to_Periphery_Published-__AGOS->g_to_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_to_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_to_Periphery_0DCapable-__AGOS->g_to_Centre_0DCapable))))
:((__AGOS->g_to_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_to_Periphery_1DCapable-__AGOS->g_to_Centre_1DCapable)))))
;
	__AGOS->calc_g_sus = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_sus_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_sus_Periphery_Published-__AGOS->g_sus_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_sus_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_sus_Periphery_0DCapable-__AGOS->g_sus_Centre_0DCapable))))
:((__AGOS->g_sus_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_sus_Periphery_1DCapable-__AGOS->g_sus_Centre_1DCapable)))))
;
	__AGOS->calc_g_K_r = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_K_r_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_K_r_Periphery_Published-__AGOS->g_K_r_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_K_r_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_K_r_Periphery_0DCapable-__AGOS->g_K_r_Centre_0DCapable))))
:((__AGOS->g_K_r_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_K_r_Periphery_1DCapable-__AGOS->g_K_r_Centre_1DCapable)))))
;
	__AGOS->calc_g_K_s = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_K_s_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_K_s_Periphery_Published-__AGOS->g_K_s_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_K_s_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_K_s_Periphery_0DCapable-__AGOS->g_K_s_Centre_0DCapable))))
:((__AGOS->g_K_s_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_K_s_Periphery_1DCapable-__AGOS->g_K_s_Centre_1DCapable)))))
;
	__AGOS->calc_g_f_Na = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_f_Na_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_f_Na_Periphery_Published-__AGOS->g_f_Na_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_f_Na_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_f_Na_Periphery_0DCapable-__AGOS->g_f_Na_Centre_0DCapable))))
:((__AGOS->g_f_Na_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_f_Na_Periphery_1DCapable-__AGOS->g_f_Na_Centre_1DCapable)))))
;
	__AGOS->calc_g_f_K = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_f_K_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_f_K_Periphery_Published-__AGOS->g_f_K_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_f_K_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_f_K_Periphery_0DCapable-__AGOS->g_f_K_Centre_0DCapable))))
:((__AGOS->g_f_K_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_f_K_Periphery_1DCapable-__AGOS->g_f_K_Centre_1DCapable)))))
;
	__AGOS->calc_g_b_Na = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_b_Na_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_b_Na_Periphery_Published-__AGOS->g_b_Na_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_b_Na_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_b_Na_Periphery_0DCapable-__AGOS->g_b_Na_Centre_0DCapable))))
:((__AGOS->g_b_Na_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_b_Na_Periphery_1DCapable-__AGOS->g_b_Na_Centre_1DCapable)))))
;
	__AGOS->calc_g_b_K = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_b_K_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_b_K_Periphery_Published-__AGOS->g_b_K_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_b_K_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_b_K_Periphery_0DCapable-__AGOS->g_b_K_Centre_0DCapable))))
:((__AGOS->g_b_K_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_b_K_Periphery_1DCapable-__AGOS->g_b_K_Centre_1DCapable)))))
;
	__AGOS->calc_g_b_Ca = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->g_b_Ca_Centre_Published+(__AGOS->calc_FCell*(__AGOS->g_b_Ca_Periphery_Published-__AGOS->g_b_Ca_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->g_b_Ca_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->g_b_Ca_Periphery_0DCapable-__AGOS->g_b_Ca_Centre_0DCapable))))
:((__AGOS->g_b_Ca_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->g_b_Ca_Periphery_1DCapable-__AGOS->g_b_Ca_Centre_1DCapable)))))
;
	__AGOS->calc_k_NaCa = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->k_NaCa_Centre_Published+(__AGOS->calc_FCell*(__AGOS->k_NaCa_Periphery_Published-__AGOS->k_NaCa_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->k_NaCa_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->k_NaCa_Periphery_0DCapable-__AGOS->k_NaCa_Centre_0DCapable))))
:((__AGOS->k_NaCa_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->k_NaCa_Periphery_1DCapable-__AGOS->k_NaCa_Centre_1DCapable)))))
;
	__AGOS->calc_i_p_max = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->i_p_max_Centre_Published+(__AGOS->calc_FCell*(__AGOS->i_p_max_Periphery_Published-__AGOS->i_p_max_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->i_p_max_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->i_p_max_Periphery_0DCapable-__AGOS->i_p_max_Centre_0DCapable))))
:((__AGOS->i_p_max_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->i_p_max_Periphery_1DCapable-__AGOS->i_p_max_Centre_1DCapable)))))
;
	__AGOS->calc_i_Ca_p_max = ((__AGOS->Version==0.0000000000e+00))
?((__AGOS->i_Ca_p_max_Centre_Published+(__AGOS->calc_FCell*(__AGOS->i_Ca_p_max_Periphery_Published-__AGOS->i_Ca_p_max_Centre_Published))))
:(((__AGOS->Version==1.0000000000e+00))
?((__AGOS->i_Ca_p_max_Centre_0DCapable+(__AGOS->calc_FCell*(__AGOS->i_Ca_p_max_Periphery_0DCapable-__AGOS->i_Ca_p_max_Centre_0DCapable))))
:((__AGOS->i_Ca_p_max_Centre_1DCapable+(__AGOS->calc_FCell*(__AGOS->i_Ca_p_max_Periphery_1DCapable-__AGOS->i_Ca_p_max_Centre_1DCapable)))))
;
	__AGOS->calc_i_Ca_L = (__AGOS->calc_g_Ca_L*((__AGOS->f_L_old_*__AGOS->d_L_old_)+(6.0000000000e-03/(1.0000000000e+00+exp(((-(__AGOS->V_old_+1.4100000000e+01))/6.0000000000e+00)))))*(__AGOS->V_old_-__AGOS->E_Ca_L));
	__AGOS->calc_i_Ca_T = (__AGOS->calc_g_Ca_T*__AGOS->d_T_old_*__AGOS->f_T_old_*(__AGOS->V_old_-__AGOS->E_Ca_T));
	__AGOS->calc_tau_f_T = (1.0000000000e+00/(__AGOS->calc_alpha_f_T+__AGOS->calc_beta_f_T));
	__AGOS->calc_i_to = (__AGOS->calc_g_to*__AGOS->q_old_*__AGOS->r_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_sus = (__AGOS->calc_g_sus*__AGOS->r_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_K_r = (__AGOS->calc_g_K_r*__AGOS->calc_P_a*__AGOS->P_i_old_*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_K_s = (__AGOS->calc_g_K_s*pow(__AGOS->xs_old_,2.0000000000e+00)*(__AGOS->V_old_-__AGOS->calc_E_K_s));
	__AGOS->calc_i_f_Na = ((__AGOS->Version!=2.0000000000e+00))
?((__AGOS->calc_g_f_Na*__AGOS->y_old_*(__AGOS->V_old_-__AGOS->calc_E_Na)))
:((__AGOS->calc_g_f_Na*__AGOS->y_old_*(__AGOS->V_old_-7.7600000000e+01)));
	__AGOS->calc_i_f_K = ((__AGOS->Version!=2.0000000000e+00))
?((__AGOS->calc_g_f_K*__AGOS->y_old_*(__AGOS->V_old_-__AGOS->calc_E_K)))
:((__AGOS->calc_g_f_K*__AGOS->y_old_*(__AGOS->V_old_+1.0200000000e+02)));
	__AGOS->calc_i_b_Na = (__AGOS->calc_g_b_Na*(__AGOS->V_old_-__AGOS->calc_E_Na));
	__AGOS->calc_i_b_K = (__AGOS->calc_g_b_K*(__AGOS->V_old_-__AGOS->calc_E_K));
	__AGOS->calc_i_b_Ca = (__AGOS->calc_g_b_Ca*(__AGOS->V_old_-__AGOS->calc_E_Ca));
	__AGOS->calc_i_NaCa = ((__AGOS->Version==0.0000000000e+00))
?(((__AGOS->calc_k_NaCa*((pow(__AGOS->Na_i,3.0000000000e+00)*__AGOS->Ca_o*exp((3.7430000000e-02*__AGOS->V_old_*__AGOS->gamma_NaCa)))-(pow(__AGOS->Na_o,3.0000000000e+00)*__AGOS->Ca_i*exp((3.7400000000e-02*__AGOS->V_old_*(__AGOS->gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(__AGOS->d_NaCa*((__AGOS->Ca_i*pow(__AGOS->Na_o,3.0000000000e+00))+(__AGOS->Ca_o*pow(__AGOS->Na_i,3.0000000000e+00)))))))
:(((__AGOS->calc_k_NaCa*((pow(__AGOS->Na_i,3.0000000000e+00)*__AGOS->Ca_o*exp((3.7430000000e-02*__AGOS->V_old_*__AGOS->gamma_NaCa)))-(pow(__AGOS->Na_o,3.0000000000e+00)*__AGOS->Ca_i*exp((3.7430000000e-02*__AGOS->V_old_*(__AGOS->gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(__AGOS->d_NaCa*((__AGOS->Ca_i*pow(__AGOS->Na_o,3.0000000000e+00))+(__AGOS->Ca_o*pow(__AGOS->Na_i,3.0000000000e+00)))))));
	__AGOS->calc_i_p = ((__AGOS->calc_i_p_max*pow((__AGOS->Na_i/(__AGOS->K_m_Na+__AGOS->Na_i)),3.0000000000e+00)*pow((__AGOS->K_o/(__AGOS->K_m_K+__AGOS->K_o)),2.0000000000e+00)*1.6000000000e+00)/(1.5000000000e+00+exp(((-(__AGOS->V_old_+6.0000000000e+01))/4.0000000000e+01))));
	__AGOS->calc_i_Ca_p = ((__AGOS->calc_i_Ca_p_max*__AGOS->Ca_i)/(__AGOS->Ca_i+4.0000000000e-04));
	__AGOS->calc_i_Na = (((((__AGOS->calc_g_Na*pow(__AGOS->m_old_,3.0000000000e+00)*__AGOS->calc_h*__AGOS->Na_o*pow(__AGOS->F,2.0000000000e+00))/(__AGOS->R*__AGOS->T))*(exp((((__AGOS->V_old_-__AGOS->calc_E_Na)*__AGOS->F)/(__AGOS->R*__AGOS->T)))-1.0000000000e+00))/(exp(((__AGOS->V_old_*__AGOS->F)/(__AGOS->R*__AGOS->T)))-1.0000000000e+00))*__AGOS->V_old_);
	__AGOS->V_lado_direito_= (((-1.0000000000e+00)/__AGOS->calc_Cm)*(__AGOS->calc_i_Na+__AGOS->calc_i_Ca_L+__AGOS->calc_i_Ca_T+__AGOS->calc_i_to+__AGOS->calc_i_sus+__AGOS->calc_i_K_r+__AGOS->calc_i_K_s+__AGOS->calc_i_f_Na+__AGOS->calc_i_f_K+__AGOS->calc_i_b_Na+__AGOS->calc_i_b_Ca+__AGOS->calc_i_b_K+__AGOS->calc_i_NaCa+__AGOS->calc_i_p+__AGOS->calc_i_Ca_p));
	__AGOS->f_L_lado_direito_= ((__AGOS->calc_f_L_infinity-__AGOS->f_L_old_)/__AGOS->calc_tau_f_L);
	__AGOS->f_T_lado_direito_= ((__AGOS->calc_f_T_infinity-__AGOS->f_T_old_)/__AGOS->calc_tau_f_T);
} //fim

void __tree2__( Solveode *__AGOS){
	__AGOS->calc_m_infinity = ((__AGOS->Version==0.0000000000e+00))
?(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-__AGOS->V_old_)/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)))
:(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+3.0320000000e+01))/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
	__AGOS->calc_tau_m = ((__AGOS->Version==0.0000000000e+00))
?(((6.2470000000e-04/((8.3200000000e-01*exp(((-3.3500000000e-01)*(__AGOS->V_old_+5.6700000000e+01))))+(6.2700000000e-01*exp((8.2000000000e-02*(__AGOS->V_old_+6.5010000000e+01))))))+4.0000000000e-05))
:(((6.2470000000e-04/((8.3221660000e-01*exp(((-3.3566000000e-01)*(__AGOS->V_old_+5.6706200000e+01))))+(6.2740000000e-01*exp((8.2300000000e-02*(__AGOS->V_old_+6.5013100000e+01))))))+4.5690000000e-05));
	__AGOS->m_lado_direito_= ((__AGOS->calc_m_infinity-__AGOS->m_old_)/__AGOS->calc_tau_m);
} //fim

void __tree3__( Solveode *__AGOS){
	__AGOS->calc_h1_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+6.6100000000e+01)/6.4000000000e+00))));
	__AGOS->calc_tau_h1 = (((3.7170000000e-06*exp(((-2.8150000000e-01)*(__AGOS->V_old_+1.7110000000e+01))))/(1.0000000000e+00+(3.7320000000e-03*exp(((-3.4260000000e-01)*(__AGOS->V_old_+3.7760000000e+01))))))+5.9770000000e-04);
	__AGOS->calc_tau_h2 = (((3.1860000000e-08*exp(((-6.2190000000e-01)*(__AGOS->V_old_+1.8800000000e+01))))/(1.0000000000e+00+(7.1890000000e-05*exp(((-6.6830000000e-01)*(__AGOS->V_old_+3.4070000000e+01))))))+3.5560000000e-03);
	__AGOS->calc_h2_infinity = __AGOS->calc_h1_infinity;
	__AGOS->h1_lado_direito_= ((__AGOS->calc_h1_infinity-__AGOS->h1_old_)/__AGOS->calc_tau_h1);
	__AGOS->h2_lado_direito_= ((__AGOS->calc_h2_infinity-__AGOS->h2_old_)/__AGOS->calc_tau_h2);
} //fim

void __tree4__( Solveode *__AGOS){
	__AGOS->calc_alpha_d_L = ((__AGOS->Version==0.0000000000e+00))
?(((((-2.8380000000e+01)*(__AGOS->V_old_+3.5000000000e+01))/(exp(((-(__AGOS->V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__AGOS->V_old_)/(exp(((-2.0800000000e-01)*__AGOS->V_old_))-1.0000000000e+00))))
:(((__AGOS->Version==1.0000000000e+00))
?(((((-2.8390000000e+01)*(__AGOS->V_old_+3.5000000000e+01))/(exp(((-(__AGOS->V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__AGOS->V_old_)/(exp(((-2.0800000000e-01)*__AGOS->V_old_))-1.0000000000e+00))))
:(((((-2.8400000000e+01)*(__AGOS->V_old_+3.5000000000e+01))/(exp(((-(__AGOS->V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__AGOS->V_old_)/(exp(((-2.0800000000e-01)*__AGOS->V_old_))-1.0000000000e+00)))))
;
	__AGOS->calc_beta_d_L = ((__AGOS->Version==1.0000000000e+00))
?(((1.1430000000e+01*(__AGOS->V_old_-5.0000000000e+00))/(exp((4.0000000000e-01*(__AGOS->V_old_-5.0000000000e+00)))-1.0000000000e+00)))
:(((1.1420000000e+01*(__AGOS->V_old_-5.0000000000e+00))/(exp((4.0000000000e-01*(__AGOS->V_old_-5.0000000000e+00)))-1.0000000000e+00)));
	__AGOS->calc_d_L_infinity = ((__AGOS->Version==0.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+2.3100000000e+01))/6.0000000000e+00)))))
:(((__AGOS->Version==1.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+2.2300000000e+01+(8.0000000000e-01*__AGOS->calc_FCell)))/6.0000000000e+00)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+2.2200000000e+01))/6.0000000000e+00))))))
;
	__AGOS->calc_tau_d_L = (2.0000000000e+00/(__AGOS->calc_alpha_d_L+__AGOS->calc_beta_d_L));
	__AGOS->d_L_lado_direito_= ((__AGOS->calc_d_L_infinity-__AGOS->d_L_old_)/__AGOS->calc_tau_d_L);
} //fim

void __tree5__( Solveode *__AGOS){
	__AGOS->calc_alpha_d_T = (1.0680000000e+03*exp(((__AGOS->V_old_+2.6300000000e+01)/3.0000000000e+01)));
	__AGOS->calc_beta_d_T = (1.0680000000e+03*exp(((-(__AGOS->V_old_+2.6300000000e+01))/3.0000000000e+01)));
	__AGOS->calc_d_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+3.7000000000e+01))/6.8000000000e+00))));
	__AGOS->calc_tau_d_T = (1.0000000000e+00/(__AGOS->calc_alpha_d_T+__AGOS->calc_beta_d_T));
	__AGOS->d_T_lado_direito_= ((__AGOS->calc_d_T_infinity-__AGOS->d_T_old_)/__AGOS->calc_tau_d_T);
} //fim

void __tree6__( Solveode *__AGOS){
	__AGOS->calc_q_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+5.9370000000e+01)/1.3100000000e+01))));
	__AGOS->calc_tau_q = ((__AGOS->Version==0.0000000000e+00))
?((1.0100000000e-02+(6.5170000000e-02/(5.7000000000e-01*exp(((-8.0000000000e-02)*(__AGOS->V_old_+4.9000000000e+01)))))+(2.4000000000e-05*exp((1.0000000000e-01*(__AGOS->V_old_+5.0930000000e+01))))))
:(((__AGOS->Version==1.0000000000e+00))
?(((1.0000000000e-03/3.0000000000e+00)*(3.0310000000e+01+(1.9550000000e+02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__AGOS->V_old_+3.9000000000e+01+(1.0000000000e+01*__AGOS->calc_FCell)))))+(7.1740000000e-01*exp(((2.7190000000e-01-(1.7190000000e-01*__AGOS->calc_FCell))*1.0000000000e+00*(__AGOS->V_old_+4.0930000000e+01+(1.0000000000e+01*__AGOS->calc_FCell))))))))))
:((1.0100000000e-02+(6.5170000000e-02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__AGOS->V_old_+3.9000000000e+01))))+(7.1740000000e-01*exp((2.7190000000e-01*(__AGOS->V_old_+4.0930000000e+01)))))))))
;
	__AGOS->q_lado_direito_= ((__AGOS->calc_q_infinity-__AGOS->q_old_)/__AGOS->calc_tau_q);
} //fim

void __tree7__( Solveode *__AGOS){
	__AGOS->calc_r_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_-1.0930000000e+01))/1.9700000000e+01))));
	__AGOS->calc_tau_r = ((__AGOS->Version==0.0000000000e+00))
?((1.0000000000e-03*(2.9800000000e+00+(1.5590000000e+01/((1.0370000000e+00*exp((9.0000000000e-02*(__AGOS->V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.2000000000e-01)*(__AGOS->V_old_+2.3840000000e+01)))))))))
:(((__AGOS->Version==1.0000000000e+00))
?((2.5000000000e-03*(1.1910000000e+00+(7.8380000000e+00/((1.0370000000e+00*exp((9.0120000000e-02*(__AGOS->V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__AGOS->V_old_+2.3840000000e+01)))))))))
:((1.0000000000e-03*(2.9800000000e+00+(1.9590000000e+01/((1.0370000000e+00*exp((9.0120000000e-02*(__AGOS->V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__AGOS->V_old_+2.3840000000e+01))))))))))
;
	__AGOS->r_lado_direito_= ((__AGOS->calc_r_infinity-__AGOS->r_old_)/__AGOS->calc_tau_r);
} //fim

void __tree8__( Solveode *__AGOS){
	__AGOS->calc_P_af_infinity = ((__AGOS->Version!=2.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+1.4200000000e+01))/1.0600000000e+01)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__AGOS->V_old_+1.3200000000e+01))/1.0600000000e+01)))));
	__AGOS->calc_tau_P_af = ((__AGOS->Version!=2.0000000000e+00))
?((1.0000000000e+00/((3.7200000000e+01*exp(((__AGOS->V_old_-9.0000000000e+00)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__AGOS->V_old_-9.0000000000e+00))/2.2500000000e+01))))))
:((1.0000000000e+00/((3.7200000000e+01*exp(((__AGOS->V_old_-1.0000000000e+01)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__AGOS->V_old_-1.0000000000e+01))/2.2500000000e+01))))));
	__AGOS->calc_tau_P_as = ((__AGOS->Version!=2.0000000000e+00))
?((1.0000000000e+00/((4.2000000000e+00*exp(((__AGOS->V_old_-9.0000000000e+00)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__AGOS->V_old_-9.0000000000e+00))/2.1600000000e+01))))))
:((1.0000000000e+00/((4.2000000000e+00*exp(((__AGOS->V_old_-1.0000000000e+01)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__AGOS->V_old_-1.0000000000e+01))/2.1600000000e+01))))));
	__AGOS->calc_P_as_infinity = __AGOS->calc_P_af_infinity;
	__AGOS->P_af_lado_direito_= ((__AGOS->calc_P_af_infinity-__AGOS->P_af_old_)/__AGOS->calc_tau_P_af);
	__AGOS->P_as_lado_direito_= ((__AGOS->calc_P_as_infinity-__AGOS->P_as_old_)/__AGOS->calc_tau_P_as);
} //fim

void __tree9__( Solveode *__AGOS){
	__AGOS->calc_tau_P_i = ((__AGOS->Version==0.0000000000e+00))
?(2.0000000000e-03)
:(((__AGOS->Version==1.0000000000e+00))
?(2.0000000000e-03)
:(6.0000000000e-03))
;
	__AGOS->calc_P_i_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__AGOS->V_old_+1.8600000000e+01)/1.0100000000e+01))));
	__AGOS->P_i_lado_direito_= ((__AGOS->calc_P_i_infinity-__AGOS->P_i_old_)/__AGOS->calc_tau_P_i);
} //fim

void __tree10__( Solveode *__AGOS){
	__AGOS->calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__AGOS->V_old_-4.0000000000e+01))/9.0000000000e+00))));
	__AGOS->calc_beta_xs = (1.0000000000e+00*exp(((-__AGOS->V_old_)/4.5000000000e+01)));
	__AGOS->xs_lado_direito_= ((__AGOS->calc_alpha_xs*(1.0000000000e+00-__AGOS->xs_old_))-(__AGOS->calc_beta_xs*__AGOS->xs_old_));
} //fim

void __tree11__( Solveode *__AGOS){
	__AGOS->calc_alpha_y = ((__AGOS->Version==0.0000000000e+00))
?((1.0000000000e+00*exp(((-(__AGOS->V_old_+7.8910000000e+01))/2.6620000000e+01))))
:((1.0000000000e+00*exp(((-(__AGOS->V_old_+7.8910000000e+01))/2.6630000000e+01))));
	__AGOS->calc_beta_y = (1.0000000000e+00*exp(((__AGOS->V_old_+7.5130000000e+01)/2.1250000000e+01)));
	__AGOS->y_lado_direito_= ((__AGOS->calc_alpha_y*(1.0000000000e+00-__AGOS->y_old_))-(__AGOS->calc_beta_y*__AGOS->y_old_));
} //fim

void __AGOS_EQUATIONS__( Solveode *__AGOS){
		const double time_new = __AGOS->time_new;
		const double time = 0.0000000000e+00;
		const double dCell = 0.0000000000e+00;
		const double Version = 1.0000000000e+00;
		const double FCellConstant = 1.0309347000e+00;
		const double CmCentre = 2.0000000000e-05;
		const double CmPeriphery = 6.5000000000e-05;
		const double g_Na_Centre_Published = 0.0000000000e+00;
		const double g_Na_Periphery_Published = 1.2000000000e-06;
		const double g_Na_Centre_0DCapable = 0.0000000000e+00;
		const double g_Na_Periphery_0DCapable = 1.2040000000e-06;
		const double g_Na_Centre_1DCapable = 0.0000000000e+00;
		const double g_Na_Periphery_1DCapable = 3.7000000000e-07;
		const double Na_o = 1.4000000000e+02;
		const double F = 9.6845000000e+04;
		const double R = 8.3140000000e+03;
		const double T = 3.1000000000e+02;
		const double g_Ca_L_Centre_Published = 5.8000000000e-03;
		const double g_Ca_L_Periphery_Published = 6.5900000000e-02;
		const double g_Ca_L_Centre_0DCapable = 5.7938000000e-03;
		const double g_Ca_L_Periphery_0DCapable = 6.5886480000e-02;
		const double g_Ca_L_Centre_1DCapable = 8.2000000000e-03;
		const double g_Ca_L_Periphery_1DCapable = 6.5900000000e-02;
		const double E_Ca_L = 4.6400000000e+01;
		const double g_Ca_T_Centre_Published = 4.3000000000e-03;
		const double g_Ca_T_Periphery_Published = 1.3900000000e-02;
		const double g_Ca_T_Centre_0DCapable = 4.2780600000e-03;
		const double g_Ca_T_Periphery_0DCapable = 1.3882300000e-02;
		const double g_Ca_T_Centre_1DCapable = 2.1000000000e-03;
		const double g_Ca_T_Periphery_1DCapable = 6.9400000000e-03;
		const double E_Ca_T = 4.5000000000e+01;
		const double g_to_Centre_Published = 4.9100000000e-03;
		const double g_to_Periphery_Published = 3.6490000000e-02;
		const double g_to_Centre_0DCapable = 4.9050000000e-03;
		const double g_to_Periphery_0DCapable = 3.6495000000e-02;
		const double g_to_Centre_1DCapable = 4.9050000000e-03;
		const double g_to_Periphery_1DCapable = 3.6500000000e-02;
		const double g_sus_Centre_Published = 6.6500000000e-05;
		const double g_sus_Periphery_Published = 1.1400000000e-02;
		const double g_sus_Centre_0DCapable = 6.6455040000e-05;
		const double g_sus_Periphery_0DCapable = 1.1383760000e-02;
		const double g_sus_Centre_1DCapable = 2.6600000000e-04;
		const double g_sus_Periphery_1DCapable = 1.1400000000e-02;
		const double g_K_r_Centre_Published = 7.9700000000e-04;
		const double g_K_r_Periphery_Published = 1.6000000000e-02;
		const double g_K_r_Centre_0DCapable = 7.9704000000e-04;
		const double g_K_r_Periphery_0DCapable = 1.6000000000e-02;
		const double g_K_r_Centre_1DCapable = 7.3800000000e-04;
		const double g_K_r_Periphery_1DCapable = 2.0800000000e-02;
		const double g_K_s_Centre_Published = 5.1800000000e-04;
		const double g_K_s_Periphery_Published = 1.0400000000e-02;
		const double g_K_s_Centre_0DCapable = 3.4450000000e-04;
		const double g_K_s_Periphery_0DCapable = 1.0400000000e-02;
		const double g_K_s_Centre_1DCapable = 3.4500000000e-04;
		const double g_K_s_Periphery_1DCapable = 1.0400000000e-02;
		const double g_f_Na_Centre_Published = 5.4800000000e-04;
		const double g_f_Na_Periphery_Published = 6.9000000000e-03;
		const double g_f_Na_Centre_0DCapable = 5.4650000000e-04;
		const double g_f_Na_Periphery_0DCapable = 6.8750000000e-03;
		const double g_f_Na_Centre_1DCapable = 4.3700000000e-04;
		const double g_f_Na_Periphery_1DCapable = 5.5000000000e-03;
		const double g_f_K_Centre_Published = 5.4800000000e-04;
		const double g_f_K_Periphery_Published = 6.9000000000e-03;
		const double g_f_K_Centre_0DCapable = 5.4650000000e-04;
		const double g_f_K_Periphery_0DCapable = 6.8750000000e-03;
		const double g_f_K_Centre_1DCapable = 4.3700000000e-04;
		const double g_f_K_Periphery_1DCapable = 5.5000000000e-03;
		const double g_b_Na_Centre_Published = 5.8000000000e-05;
		const double g_b_Na_Periphery_Published = 1.8900000000e-04;
		const double g_b_Na_Centre_0DCapable = 5.8181800000e-05;
		const double g_b_Na_Periphery_0DCapable = 1.8880000000e-04;
		const double g_b_Na_Centre_1DCapable = 5.8000000000e-05;
		const double g_b_Na_Periphery_1DCapable = 1.8900000000e-04;
		const double g_b_K_Centre_Published = 2.5200000000e-05;
		const double g_b_K_Periphery_Published = 8.1900000000e-05;
		const double g_b_K_Centre_0DCapable = 2.5236360000e-05;
		const double g_b_K_Periphery_0DCapable = 8.1892000000e-05;
		const double g_b_K_Centre_1DCapable = 2.5200000000e-05;
		const double g_b_K_Periphery_1DCapable = 8.1900000000e-05;
		const double g_b_Ca_Centre_Published = 1.3200000000e-05;
		const double g_b_Ca_Periphery_Published = 4.3000000000e-05;
		const double g_b_Ca_Centre_0DCapable = 1.3236000000e-05;
		const double g_b_Ca_Periphery_0DCapable = 4.2952000000e-05;
		const double g_b_Ca_Centre_1DCapable = 1.3230000000e-05;
		const double g_b_Ca_Periphery_1DCapable = 4.2900000000e-05;
		const double k_NaCa_Centre_Published = 2.7000000000e-06;
		const double k_NaCa_Periphery_Published = 8.8000000000e-06;
		const double k_NaCa_Centre_0DCapable = 2.7229000000e-06;
		const double k_NaCa_Periphery_0DCapable = 8.8358400000e-06;
		const double k_NaCa_Centre_1DCapable = 2.8000000000e-06;
		const double k_NaCa_Periphery_1DCapable = 8.8000000000e-06;
		const double Na_i = 8.0000000000e+00;
		const double Ca_o = 2.0000000000e+00;
		const double gamma_NaCa = 5.0000000000e-01;
		const double Ca_i = 1.0000000000e-04;
		const double d_NaCa = 1.0000000000e-04;
		const double i_p_max_Centre_Published = 4.7800000000e-02;
		const double i_p_max_Periphery_Published = 1.6000000000e-01;
		const double i_p_max_Centre_0DCapable = 4.7825450000e-02;
		const double i_p_max_Periphery_0DCapable = 1.5519360000e-01;
		const double i_p_max_Centre_1DCapable = 4.7800000000e-02;
		const double i_p_max_Periphery_1DCapable = 1.6000000000e-01;
		const double K_m_Na = 5.6400000000e+00;
		const double K_o = 5.4000000000e+00;
		const double K_m_K = 6.2100000000e-01;
		const double i_Ca_p_max_Centre_Published = 0.0000000000e+00;
		const double i_Ca_p_max_Periphery_Published = 0.0000000000e+00;
		const double i_Ca_p_max_Centre_0DCapable = 0.0000000000e+00;
		const double i_Ca_p_max_Periphery_0DCapable = 0.0000000000e+00;
		const double i_Ca_p_max_Centre_1DCapable = 4.2000000000e-03;
		const double i_Ca_p_max_Periphery_1DCapable = 3.3390000000e-02;
		const double K_i = 1.4000000000e+02;
		const double V_old_= __AGOS->V_old_;
		const double m_old_= __AGOS->m_old_;
		const double h1_old_= __AGOS->h1_old_;
		const double h2_old_= __AGOS->h2_old_;
		const double d_L_old_= __AGOS->d_L_old_;
		const double f_L_old_= __AGOS->f_L_old_;
		const double d_T_old_= __AGOS->d_T_old_;
		const double f_T_old_= __AGOS->f_T_old_;
		const double q_old_= __AGOS->q_old_;
		const double r_old_= __AGOS->r_old_;
		const double P_af_old_= __AGOS->P_af_old_;
		const double P_as_old_= __AGOS->P_as_old_;
		const double P_i_old_= __AGOS->P_i_old_;
		const double xs_old_= __AGOS->xs_old_;
		const double y_old_= __AGOS->y_old_;
	const double calc_FCell = ((Version==0.0000000000e+00))
?(((1.0700000000e+00*((3.0000000000e+00*dCell)-1.0000000000e-01))/(3.0000000000e+00*(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*dCell)-2.0500000000e+00))/2.9500000000e-01)))))))
:(((Version==1.0000000000e+00))
?(((FCellConstant*dCell)/(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*dCell)-2.0500000000e+00))/2.9500000000e-01))))))
:(((1.0700000000e+00*2.9000000000e+01*dCell)/(3.0000000000e+01*(1.0000000000e+00+(7.7450000000e-01*exp(((-((2.9000000000e+01*dCell)-2.4500000000e+01))/1.9500000000e+00))))))))
;
	const double calc_F_Na = ((Version==0.0000000000e+00))
?((((9.5200000000e-02*exp(((-6.3000000000e-02)*(V_old_+3.4400000000e+01))))/(1.0000000000e+00+(1.6600000000e+00*exp(((-2.2500000000e-01)*(V_old_+6.3700000000e+01))))))+8.6900000000e-02))
:((((9.5180000000e-02*exp(((-6.3060000000e-02)*(V_old_+3.4400000000e+01))))/(1.0000000000e+00+(1.6620000000e+00*exp(((-2.2510000000e-01)*(V_old_+6.3700000000e+01))))))+8.6930000000e-02));
	const double calc_alpha_f_L = ((Version==1.0000000000e+00))
?(((3.7500000000e+00*(V_old_+2.8000000000e+01))/(exp(((V_old_+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)))
:(((3.1200000000e+00*(V_old_+2.8000000000e+01))/(exp(((V_old_+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
	const double calc_beta_f_L = ((Version==1.0000000000e+00))
?((3.0000000000e+01/(1.0000000000e+00+exp(((-(V_old_+2.8000000000e+01))/4.0000000000e+00)))))
:((2.5000000000e+01/(1.0000000000e+00+exp(((-(V_old_+2.8000000000e+01))/4.0000000000e+00)))));
	const double calc_f_L_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+4.5000000000e+01)/5.0000000000e+00))));
	const double calc_beta_f_T = ((Version==1.0000000000e+00))
?((1.5000000000e+01*exp(((V_old_+7.1000000000e+01)/1.5380000000e+01))))
:((1.5000000000e+01*exp(((V_old_+7.1700000000e+01)/1.5380000000e+01))));
	const double calc_f_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+7.1000000000e+01)/9.0000000000e+00))));
	const double calc_P_a = ((6.0000000000e-01*P_af_old_)+(4.0000000000e-01*P_as_old_));
	const double calc_E_Na = (((R*T)/F)*log((Na_o/Na_i)));
	const double calc_E_K = (((R*T)/F)*log((K_o/K_i)));
	const double calc_E_Ca = (((R*T)/(2.0000000000e+00*F))*log((Ca_o/Ca_i)));
	const double calc_E_K_s = ((Version==0.0000000000e+00))
?((((R*T)/F)*log(((K_o+(1.2000000000e-01*Na_o))/(K_i+(1.2000000000e-01*Na_i))))))
:((((R*T)/F)*log(((K_o+(3.0000000000e-02*Na_o))/(K_i+(3.0000000000e-02*Na_i))))));
	const double calc_Cm = (CmCentre+(calc_FCell*(CmPeriphery-CmCentre)));
	const double calc_g_Na = ((Version==0.0000000000e+00))
?((g_Na_Centre_Published+(calc_FCell*(g_Na_Periphery_Published-g_Na_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_Na_Centre_0DCapable+(calc_FCell*(g_Na_Periphery_0DCapable-g_Na_Centre_0DCapable))))
:((g_Na_Centre_1DCapable+(calc_FCell*(g_Na_Periphery_1DCapable-g_Na_Centre_1DCapable)))))
;
	const double calc_h = (((1.0000000000e+00-calc_F_Na)*h1_old_)+(calc_F_Na*h2_old_));
	const double calc_g_Ca_L = ((Version==0.0000000000e+00))
?((g_Ca_L_Centre_Published+(calc_FCell*(g_Ca_L_Periphery_Published-g_Ca_L_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_Ca_L_Centre_0DCapable+(calc_FCell*(g_Ca_L_Periphery_0DCapable-g_Ca_L_Centre_0DCapable))))
:((g_Ca_L_Centre_1DCapable+(calc_FCell*(g_Ca_L_Periphery_1DCapable-g_Ca_L_Centre_1DCapable)))))
;
	const double calc_tau_f_L = ((Version==1.0000000000e+00))
?(((1.2000000000e+00-(2.0000000000e-01*calc_FCell))/(calc_alpha_f_L+calc_beta_f_L)))
:((1.0000000000e+00/(calc_alpha_f_L+calc_beta_f_L)));
	const double calc_g_Ca_T = ((Version==0.0000000000e+00))
?((g_Ca_T_Centre_Published+(calc_FCell*(g_Ca_T_Periphery_Published-g_Ca_T_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_Ca_T_Centre_0DCapable+(calc_FCell*(g_Ca_T_Periphery_0DCapable-g_Ca_T_Centre_0DCapable))))
:((g_Ca_T_Centre_1DCapable+(calc_FCell*(g_Ca_T_Periphery_1DCapable-g_Ca_T_Centre_1DCapable)))))
;
	const double calc_alpha_f_T = ((Version==1.0000000000e+00))
?((1.5300000000e+01*exp(((-(V_old_+7.1000000000e+01+(7.0000000000e-01*calc_FCell)))/8.3300000000e+01))))
:((1.5300000000e+01*exp(((-(V_old_+7.1700000000e+01))/8.3300000000e+01))));
	const double calc_g_to = ((Version==0.0000000000e+00))
?((g_to_Centre_Published+(calc_FCell*(g_to_Periphery_Published-g_to_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_to_Centre_0DCapable+(calc_FCell*(g_to_Periphery_0DCapable-g_to_Centre_0DCapable))))
:((g_to_Centre_1DCapable+(calc_FCell*(g_to_Periphery_1DCapable-g_to_Centre_1DCapable)))))
;
	const double calc_g_sus = ((Version==0.0000000000e+00))
?((g_sus_Centre_Published+(calc_FCell*(g_sus_Periphery_Published-g_sus_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_sus_Centre_0DCapable+(calc_FCell*(g_sus_Periphery_0DCapable-g_sus_Centre_0DCapable))))
:((g_sus_Centre_1DCapable+(calc_FCell*(g_sus_Periphery_1DCapable-g_sus_Centre_1DCapable)))))
;
	const double calc_g_K_r = ((Version==0.0000000000e+00))
?((g_K_r_Centre_Published+(calc_FCell*(g_K_r_Periphery_Published-g_K_r_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_K_r_Centre_0DCapable+(calc_FCell*(g_K_r_Periphery_0DCapable-g_K_r_Centre_0DCapable))))
:((g_K_r_Centre_1DCapable+(calc_FCell*(g_K_r_Periphery_1DCapable-g_K_r_Centre_1DCapable)))))
;
	const double calc_g_K_s = ((Version==0.0000000000e+00))
?((g_K_s_Centre_Published+(calc_FCell*(g_K_s_Periphery_Published-g_K_s_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_K_s_Centre_0DCapable+(calc_FCell*(g_K_s_Periphery_0DCapable-g_K_s_Centre_0DCapable))))
:((g_K_s_Centre_1DCapable+(calc_FCell*(g_K_s_Periphery_1DCapable-g_K_s_Centre_1DCapable)))))
;
	const double calc_g_f_Na = ((Version==0.0000000000e+00))
?((g_f_Na_Centre_Published+(calc_FCell*(g_f_Na_Periphery_Published-g_f_Na_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_f_Na_Centre_0DCapable+(calc_FCell*(g_f_Na_Periphery_0DCapable-g_f_Na_Centre_0DCapable))))
:((g_f_Na_Centre_1DCapable+(calc_FCell*(g_f_Na_Periphery_1DCapable-g_f_Na_Centre_1DCapable)))))
;
	const double calc_g_f_K = ((Version==0.0000000000e+00))
?((g_f_K_Centre_Published+(calc_FCell*(g_f_K_Periphery_Published-g_f_K_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_f_K_Centre_0DCapable+(calc_FCell*(g_f_K_Periphery_0DCapable-g_f_K_Centre_0DCapable))))
:((g_f_K_Centre_1DCapable+(calc_FCell*(g_f_K_Periphery_1DCapable-g_f_K_Centre_1DCapable)))))
;
	const double calc_g_b_Na = ((Version==0.0000000000e+00))
?((g_b_Na_Centre_Published+(calc_FCell*(g_b_Na_Periphery_Published-g_b_Na_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_b_Na_Centre_0DCapable+(calc_FCell*(g_b_Na_Periphery_0DCapable-g_b_Na_Centre_0DCapable))))
:((g_b_Na_Centre_1DCapable+(calc_FCell*(g_b_Na_Periphery_1DCapable-g_b_Na_Centre_1DCapable)))))
;
	const double calc_g_b_K = ((Version==0.0000000000e+00))
?((g_b_K_Centre_Published+(calc_FCell*(g_b_K_Periphery_Published-g_b_K_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_b_K_Centre_0DCapable+(calc_FCell*(g_b_K_Periphery_0DCapable-g_b_K_Centre_0DCapable))))
:((g_b_K_Centre_1DCapable+(calc_FCell*(g_b_K_Periphery_1DCapable-g_b_K_Centre_1DCapable)))))
;
	const double calc_g_b_Ca = ((Version==0.0000000000e+00))
?((g_b_Ca_Centre_Published+(calc_FCell*(g_b_Ca_Periphery_Published-g_b_Ca_Centre_Published))))
:(((Version==1.0000000000e+00))
?((g_b_Ca_Centre_0DCapable+(calc_FCell*(g_b_Ca_Periphery_0DCapable-g_b_Ca_Centre_0DCapable))))
:((g_b_Ca_Centre_1DCapable+(calc_FCell*(g_b_Ca_Periphery_1DCapable-g_b_Ca_Centre_1DCapable)))))
;
	const double calc_k_NaCa = ((Version==0.0000000000e+00))
?((k_NaCa_Centre_Published+(calc_FCell*(k_NaCa_Periphery_Published-k_NaCa_Centre_Published))))
:(((Version==1.0000000000e+00))
?((k_NaCa_Centre_0DCapable+(calc_FCell*(k_NaCa_Periphery_0DCapable-k_NaCa_Centre_0DCapable))))
:((k_NaCa_Centre_1DCapable+(calc_FCell*(k_NaCa_Periphery_1DCapable-k_NaCa_Centre_1DCapable)))))
;
	const double calc_i_p_max = ((Version==0.0000000000e+00))
?((i_p_max_Centre_Published+(calc_FCell*(i_p_max_Periphery_Published-i_p_max_Centre_Published))))
:(((Version==1.0000000000e+00))
?((i_p_max_Centre_0DCapable+(calc_FCell*(i_p_max_Periphery_0DCapable-i_p_max_Centre_0DCapable))))
:((i_p_max_Centre_1DCapable+(calc_FCell*(i_p_max_Periphery_1DCapable-i_p_max_Centre_1DCapable)))))
;
	const double calc_i_Ca_p_max = ((Version==0.0000000000e+00))
?((i_Ca_p_max_Centre_Published+(calc_FCell*(i_Ca_p_max_Periphery_Published-i_Ca_p_max_Centre_Published))))
:(((Version==1.0000000000e+00))
?((i_Ca_p_max_Centre_0DCapable+(calc_FCell*(i_Ca_p_max_Periphery_0DCapable-i_Ca_p_max_Centre_0DCapable))))
:((i_Ca_p_max_Centre_1DCapable+(calc_FCell*(i_Ca_p_max_Periphery_1DCapable-i_Ca_p_max_Centre_1DCapable)))))
;
	const double calc_i_Ca_L = (calc_g_Ca_L*((f_L_old_*d_L_old_)+(6.0000000000e-03/(1.0000000000e+00+exp(((-(V_old_+1.4100000000e+01))/6.0000000000e+00)))))*(V_old_-E_Ca_L));
	const double calc_i_Ca_T = (calc_g_Ca_T*d_T_old_*f_T_old_*(V_old_-E_Ca_T));
	const double calc_tau_f_T = (1.0000000000e+00/(calc_alpha_f_T+calc_beta_f_T));
	const double calc_i_to = (calc_g_to*q_old_*r_old_*(V_old_-calc_E_K));
	const double calc_i_sus = (calc_g_sus*r_old_*(V_old_-calc_E_K));
	const double calc_i_K_r = (calc_g_K_r*calc_P_a*P_i_old_*(V_old_-calc_E_K));
	const double calc_i_K_s = (calc_g_K_s*pow(xs_old_,2.0000000000e+00)*(V_old_-calc_E_K_s));
	const double calc_i_f_Na = ((Version!=2.0000000000e+00))
?((calc_g_f_Na*y_old_*(V_old_-calc_E_Na)))
:((calc_g_f_Na*y_old_*(V_old_-7.7600000000e+01)));
	const double calc_i_f_K = ((Version!=2.0000000000e+00))
?((calc_g_f_K*y_old_*(V_old_-calc_E_K)))
:((calc_g_f_K*y_old_*(V_old_+1.0200000000e+02)));
	const double calc_i_b_Na = (calc_g_b_Na*(V_old_-calc_E_Na));
	const double calc_i_b_K = (calc_g_b_K*(V_old_-calc_E_K));
	const double calc_i_b_Ca = (calc_g_b_Ca*(V_old_-calc_E_Ca));
	const double calc_i_NaCa = ((Version==0.0000000000e+00))
?(((calc_k_NaCa*((pow(Na_i,3.0000000000e+00)*Ca_o*exp((3.7430000000e-02*V_old_*gamma_NaCa)))-(pow(Na_o,3.0000000000e+00)*Ca_i*exp((3.7400000000e-02*V_old_*(gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(d_NaCa*((Ca_i*pow(Na_o,3.0000000000e+00))+(Ca_o*pow(Na_i,3.0000000000e+00)))))))
:(((calc_k_NaCa*((pow(Na_i,3.0000000000e+00)*Ca_o*exp((3.7430000000e-02*V_old_*gamma_NaCa)))-(pow(Na_o,3.0000000000e+00)*Ca_i*exp((3.7430000000e-02*V_old_*(gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(d_NaCa*((Ca_i*pow(Na_o,3.0000000000e+00))+(Ca_o*pow(Na_i,3.0000000000e+00)))))));
	const double calc_i_p = ((calc_i_p_max*pow((Na_i/(K_m_Na+Na_i)),3.0000000000e+00)*pow((K_o/(K_m_K+K_o)),2.0000000000e+00)*1.6000000000e+00)/(1.5000000000e+00+exp(((-(V_old_+6.0000000000e+01))/4.0000000000e+01))));
	const double calc_i_Ca_p = ((calc_i_Ca_p_max*Ca_i)/(Ca_i+4.0000000000e-04));
	const double calc_i_Na = (((((calc_g_Na*pow(m_old_,3.0000000000e+00)*calc_h*Na_o*pow(F,2.0000000000e+00))/(R*T))*(exp((((V_old_-calc_E_Na)*F)/(R*T)))-1.0000000000e+00))/(exp(((V_old_*F)/(R*T)))-1.0000000000e+00))*V_old_);
	__AGOS->V_lado_direito_= (((-1.0000000000e+00)/calc_Cm)*(calc_i_Na+calc_i_Ca_L+calc_i_Ca_T+calc_i_to+calc_i_sus+calc_i_K_r+calc_i_K_s+calc_i_f_Na+calc_i_f_K+calc_i_b_Na+calc_i_b_Ca+calc_i_b_K+calc_i_NaCa+calc_i_p+calc_i_Ca_p));
	__AGOS->f_L_lado_direito_= ((calc_f_L_infinity-f_L_old_)/calc_tau_f_L);
	__AGOS->f_T_lado_direito_= ((calc_f_T_infinity-f_T_old_)/calc_tau_f_T);
	const double calc_m_infinity = ((Version==0.0000000000e+00))
?(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-V_old_)/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)))
:(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+3.0320000000e+01))/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
	const double calc_tau_m = ((Version==0.0000000000e+00))
?(((6.2470000000e-04/((8.3200000000e-01*exp(((-3.3500000000e-01)*(V_old_+5.6700000000e+01))))+(6.2700000000e-01*exp((8.2000000000e-02*(V_old_+6.5010000000e+01))))))+4.0000000000e-05))
:(((6.2470000000e-04/((8.3221660000e-01*exp(((-3.3566000000e-01)*(V_old_+5.6706200000e+01))))+(6.2740000000e-01*exp((8.2300000000e-02*(V_old_+6.5013100000e+01))))))+4.5690000000e-05));
	__AGOS->m_lado_direito_= ((calc_m_infinity-m_old_)/calc_tau_m);
	const double calc_h1_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+6.6100000000e+01)/6.4000000000e+00))));
	const double calc_tau_h1 = (((3.7170000000e-06*exp(((-2.8150000000e-01)*(V_old_+1.7110000000e+01))))/(1.0000000000e+00+(3.7320000000e-03*exp(((-3.4260000000e-01)*(V_old_+3.7760000000e+01))))))+5.9770000000e-04);
	const double calc_tau_h2 = (((3.1860000000e-08*exp(((-6.2190000000e-01)*(V_old_+1.8800000000e+01))))/(1.0000000000e+00+(7.1890000000e-05*exp(((-6.6830000000e-01)*(V_old_+3.4070000000e+01))))))+3.5560000000e-03);
	const double calc_h2_infinity = calc_h1_infinity;
	__AGOS->h1_lado_direito_= ((calc_h1_infinity-h1_old_)/calc_tau_h1);
	__AGOS->h2_lado_direito_= ((calc_h2_infinity-h2_old_)/calc_tau_h2);
	const double calc_alpha_d_L = ((Version==0.0000000000e+00))
?(((((-2.8380000000e+01)*(V_old_+3.5000000000e+01))/(exp(((-(V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*V_old_)/(exp(((-2.0800000000e-01)*V_old_))-1.0000000000e+00))))
:(((Version==1.0000000000e+00))
?(((((-2.8390000000e+01)*(V_old_+3.5000000000e+01))/(exp(((-(V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*V_old_)/(exp(((-2.0800000000e-01)*V_old_))-1.0000000000e+00))))
:(((((-2.8400000000e+01)*(V_old_+3.5000000000e+01))/(exp(((-(V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*V_old_)/(exp(((-2.0800000000e-01)*V_old_))-1.0000000000e+00)))))
;
	const double calc_beta_d_L = ((Version==1.0000000000e+00))
?(((1.1430000000e+01*(V_old_-5.0000000000e+00))/(exp((4.0000000000e-01*(V_old_-5.0000000000e+00)))-1.0000000000e+00)))
:(((1.1420000000e+01*(V_old_-5.0000000000e+00))/(exp((4.0000000000e-01*(V_old_-5.0000000000e+00)))-1.0000000000e+00)));
	const double calc_d_L_infinity = ((Version==0.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+2.3100000000e+01))/6.0000000000e+00)))))
:(((Version==1.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+2.2300000000e+01+(8.0000000000e-01*calc_FCell)))/6.0000000000e+00)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+2.2200000000e+01))/6.0000000000e+00))))))
;
	const double calc_tau_d_L = (2.0000000000e+00/(calc_alpha_d_L+calc_beta_d_L));
	__AGOS->d_L_lado_direito_= ((calc_d_L_infinity-d_L_old_)/calc_tau_d_L);
	const double calc_alpha_d_T = (1.0680000000e+03*exp(((V_old_+2.6300000000e+01)/3.0000000000e+01)));
	const double calc_beta_d_T = (1.0680000000e+03*exp(((-(V_old_+2.6300000000e+01))/3.0000000000e+01)));
	const double calc_d_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+3.7000000000e+01))/6.8000000000e+00))));
	const double calc_tau_d_T = (1.0000000000e+00/(calc_alpha_d_T+calc_beta_d_T));
	__AGOS->d_T_lado_direito_= ((calc_d_T_infinity-d_T_old_)/calc_tau_d_T);
	const double calc_q_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+5.9370000000e+01)/1.3100000000e+01))));
	const double calc_tau_q = ((Version==0.0000000000e+00))
?((1.0100000000e-02+(6.5170000000e-02/(5.7000000000e-01*exp(((-8.0000000000e-02)*(V_old_+4.9000000000e+01)))))+(2.4000000000e-05*exp((1.0000000000e-01*(V_old_+5.0930000000e+01))))))
:(((Version==1.0000000000e+00))
?(((1.0000000000e-03/3.0000000000e+00)*(3.0310000000e+01+(1.9550000000e+02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(V_old_+3.9000000000e+01+(1.0000000000e+01*calc_FCell)))))+(7.1740000000e-01*exp(((2.7190000000e-01-(1.7190000000e-01*calc_FCell))*1.0000000000e+00*(V_old_+4.0930000000e+01+(1.0000000000e+01*calc_FCell))))))))))
:((1.0100000000e-02+(6.5170000000e-02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(V_old_+3.9000000000e+01))))+(7.1740000000e-01*exp((2.7190000000e-01*(V_old_+4.0930000000e+01)))))))))
;
	__AGOS->q_lado_direito_= ((calc_q_infinity-q_old_)/calc_tau_q);
	const double calc_r_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_-1.0930000000e+01))/1.9700000000e+01))));
	const double calc_tau_r = ((Version==0.0000000000e+00))
?((1.0000000000e-03*(2.9800000000e+00+(1.5590000000e+01/((1.0370000000e+00*exp((9.0000000000e-02*(V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.2000000000e-01)*(V_old_+2.3840000000e+01)))))))))
:(((Version==1.0000000000e+00))
?((2.5000000000e-03*(1.1910000000e+00+(7.8380000000e+00/((1.0370000000e+00*exp((9.0120000000e-02*(V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(V_old_+2.3840000000e+01)))))))))
:((1.0000000000e-03*(2.9800000000e+00+(1.9590000000e+01/((1.0370000000e+00*exp((9.0120000000e-02*(V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(V_old_+2.3840000000e+01))))))))))
;
	__AGOS->r_lado_direito_= ((calc_r_infinity-r_old_)/calc_tau_r);
	const double calc_P_af_infinity = ((Version!=2.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+1.4200000000e+01))/1.0600000000e+01)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+1.3200000000e+01))/1.0600000000e+01)))));
	const double calc_tau_P_af = ((Version!=2.0000000000e+00))
?((1.0000000000e+00/((3.7200000000e+01*exp(((V_old_-9.0000000000e+00)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(V_old_-9.0000000000e+00))/2.2500000000e+01))))))
:((1.0000000000e+00/((3.7200000000e+01*exp(((V_old_-1.0000000000e+01)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(V_old_-1.0000000000e+01))/2.2500000000e+01))))));
	const double calc_tau_P_as = ((Version!=2.0000000000e+00))
?((1.0000000000e+00/((4.2000000000e+00*exp(((V_old_-9.0000000000e+00)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(V_old_-9.0000000000e+00))/2.1600000000e+01))))))
:((1.0000000000e+00/((4.2000000000e+00*exp(((V_old_-1.0000000000e+01)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(V_old_-1.0000000000e+01))/2.1600000000e+01))))));
	const double calc_P_as_infinity = calc_P_af_infinity;
	__AGOS->P_af_lado_direito_= ((calc_P_af_infinity-P_af_old_)/calc_tau_P_af);
	__AGOS->P_as_lado_direito_= ((calc_P_as_infinity-P_as_old_)/calc_tau_P_as);
	const double calc_tau_P_i = ((Version==0.0000000000e+00))
?(2.0000000000e-03)
:(((Version==1.0000000000e+00))
?(2.0000000000e-03)
:(6.0000000000e-03))
;
	const double calc_P_i_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((V_old_+1.8600000000e+01)/1.0100000000e+01))));
	__AGOS->P_i_lado_direito_= ((calc_P_i_infinity-P_i_old_)/calc_tau_P_i);
	const double calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(V_old_-4.0000000000e+01))/9.0000000000e+00))));
	const double calc_beta_xs = (1.0000000000e+00*exp(((-V_old_)/4.5000000000e+01)));
	__AGOS->xs_lado_direito_= ((calc_alpha_xs*(1.0000000000e+00-xs_old_))-(calc_beta_xs*xs_old_));
	const double calc_alpha_y = ((Version==0.0000000000e+00))
?((1.0000000000e+00*exp(((-(V_old_+7.8910000000e+01))/2.6620000000e+01))))
:((1.0000000000e+00*exp(((-(V_old_+7.8910000000e+01))/2.6630000000e+01))));
	const double calc_beta_y = (1.0000000000e+00*exp(((V_old_+7.5130000000e+01)/2.1250000000e+01)));
	__AGOS->y_lado_direito_= ((calc_alpha_y*(1.0000000000e+00-y_old_))-(calc_beta_y*y_old_));
} //fim
typedef struct str__rightHandSideFunction{
	void (*function)(Solveode*);
}typ_rightHandSideFunction;
typ_rightHandSideFunction rightHandSideFunction;
typ_rightHandSideFunction forest[11];
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
		time = 0.00000000e+00;
		dCell = 0.00000000e+00;
		Version = 1.00000000e+00;
		FCellConstant = 1.03093470e+00;
		CmCentre = 2.00000000e-05;
		CmPeriphery = 6.50000000e-05;
		g_Na_Centre_Published = 0.00000000e+00;
		g_Na_Periphery_Published = 1.20000000e-06;
		g_Na_Centre_0DCapable = 0.00000000e+00;
		g_Na_Periphery_0DCapable = 1.20400000e-06;
		g_Na_Centre_1DCapable = 0.00000000e+00;
		g_Na_Periphery_1DCapable = 3.70000000e-07;
		Na_o = 1.40000000e+02;
		F = 9.68450000e+04;
		R = 8.31400000e+03;
		T = 3.10000000e+02;
		g_Ca_L_Centre_Published = 5.80000000e-03;
		g_Ca_L_Periphery_Published = 6.59000000e-02;
		g_Ca_L_Centre_0DCapable = 5.79380000e-03;
		g_Ca_L_Periphery_0DCapable = 6.58864800e-02;
		g_Ca_L_Centre_1DCapable = 8.20000000e-03;
		g_Ca_L_Periphery_1DCapable = 6.59000000e-02;
		E_Ca_L = 4.64000000e+01;
		g_Ca_T_Centre_Published = 4.30000000e-03;
		g_Ca_T_Periphery_Published = 1.39000000e-02;
		g_Ca_T_Centre_0DCapable = 4.27806000e-03;
		g_Ca_T_Periphery_0DCapable = 1.38823000e-02;
		g_Ca_T_Centre_1DCapable = 2.10000000e-03;
		g_Ca_T_Periphery_1DCapable = 6.94000000e-03;
		E_Ca_T = 4.50000000e+01;
		g_to_Centre_Published = 4.91000000e-03;
		g_to_Periphery_Published = 3.64900000e-02;
		g_to_Centre_0DCapable = 4.90500000e-03;
		g_to_Periphery_0DCapable = 3.64950000e-02;
		g_to_Centre_1DCapable = 4.90500000e-03;
		g_to_Periphery_1DCapable = 3.65000000e-02;
		g_sus_Centre_Published = 6.65000000e-05;
		g_sus_Periphery_Published = 1.14000000e-02;
		g_sus_Centre_0DCapable = 6.64550400e-05;
		g_sus_Periphery_0DCapable = 1.13837600e-02;
		g_sus_Centre_1DCapable = 2.66000000e-04;
		g_sus_Periphery_1DCapable = 1.14000000e-02;
		g_K_r_Centre_Published = 7.97000000e-04;
		g_K_r_Periphery_Published = 1.60000000e-02;
		g_K_r_Centre_0DCapable = 7.97040000e-04;
		g_K_r_Periphery_0DCapable = 1.60000000e-02;
		g_K_r_Centre_1DCapable = 7.38000000e-04;
		g_K_r_Periphery_1DCapable = 2.08000000e-02;
		g_K_s_Centre_Published = 5.18000000e-04;
		g_K_s_Periphery_Published = 1.04000000e-02;
		g_K_s_Centre_0DCapable = 3.44500000e-04;
		g_K_s_Periphery_0DCapable = 1.04000000e-02;
		g_K_s_Centre_1DCapable = 3.45000000e-04;
		g_K_s_Periphery_1DCapable = 1.04000000e-02;
		g_f_Na_Centre_Published = 5.48000000e-04;
		g_f_Na_Periphery_Published = 6.90000000e-03;
		g_f_Na_Centre_0DCapable = 5.46500000e-04;
		g_f_Na_Periphery_0DCapable = 6.87500000e-03;
		g_f_Na_Centre_1DCapable = 4.37000000e-04;
		g_f_Na_Periphery_1DCapable = 5.50000000e-03;
		g_f_K_Centre_Published = 5.48000000e-04;
		g_f_K_Periphery_Published = 6.90000000e-03;
		g_f_K_Centre_0DCapable = 5.46500000e-04;
		g_f_K_Periphery_0DCapable = 6.87500000e-03;
		g_f_K_Centre_1DCapable = 4.37000000e-04;
		g_f_K_Periphery_1DCapable = 5.50000000e-03;
		g_b_Na_Centre_Published = 5.80000000e-05;
		g_b_Na_Periphery_Published = 1.89000000e-04;
		g_b_Na_Centre_0DCapable = 5.81818000e-05;
		g_b_Na_Periphery_0DCapable = 1.88800000e-04;
		g_b_Na_Centre_1DCapable = 5.80000000e-05;
		g_b_Na_Periphery_1DCapable = 1.89000000e-04;
		g_b_K_Centre_Published = 2.52000000e-05;
		g_b_K_Periphery_Published = 8.19000000e-05;
		g_b_K_Centre_0DCapable = 2.52363600e-05;
		g_b_K_Periphery_0DCapable = 8.18920000e-05;
		g_b_K_Centre_1DCapable = 2.52000000e-05;
		g_b_K_Periphery_1DCapable = 8.19000000e-05;
		g_b_Ca_Centre_Published = 1.32000000e-05;
		g_b_Ca_Periphery_Published = 4.30000000e-05;
		g_b_Ca_Centre_0DCapable = 1.32360000e-05;
		g_b_Ca_Periphery_0DCapable = 4.29520000e-05;
		g_b_Ca_Centre_1DCapable = 1.32300000e-05;
		g_b_Ca_Periphery_1DCapable = 4.29000000e-05;
		k_NaCa_Centre_Published = 2.70000000e-06;
		k_NaCa_Periphery_Published = 8.80000000e-06;
		k_NaCa_Centre_0DCapable = 2.72290000e-06;
		k_NaCa_Periphery_0DCapable = 8.83584000e-06;
		k_NaCa_Centre_1DCapable = 2.80000000e-06;
		k_NaCa_Periphery_1DCapable = 8.80000000e-06;
		Na_i = 8.00000000e+00;
		Ca_o = 2.00000000e+00;
		gamma_NaCa = 5.00000000e-01;
		Ca_i = 1.00000000e-04;
		d_NaCa = 1.00000000e-04;
		i_p_max_Centre_Published = 4.78000000e-02;
		i_p_max_Periphery_Published = 1.60000000e-01;
		i_p_max_Centre_0DCapable = 4.78254500e-02;
		i_p_max_Periphery_0DCapable = 1.55193600e-01;
		i_p_max_Centre_1DCapable = 4.78000000e-02;
		i_p_max_Periphery_1DCapable = 1.60000000e-01;
		K_m_Na = 5.64000000e+00;
		K_o = 5.40000000e+00;
		K_m_K = 6.21000000e-01;
		i_Ca_p_max_Centre_Published = 0.00000000e+00;
		i_Ca_p_max_Periphery_Published = 0.00000000e+00;
		i_Ca_p_max_Centre_0DCapable = 0.00000000e+00;
		i_Ca_p_max_Periphery_0DCapable = 0.00000000e+00;
		i_Ca_p_max_Centre_1DCapable = 4.20000000e-03;
		i_Ca_p_max_Periphery_1DCapable = 3.33900000e-02;
		K_i = 1.40000000e+02;
		dtime = 0.0; time_vec__ = NULL;
		V = NULL;
		V_ini_ = -3.90135585e+01;
		m = NULL;
		m_ini_ = 9.23617017e-02;
		h1 = NULL;
		h1_ini_ = 1.59053803e-02;
		h2 = NULL;
		h2_ini_ = 1.44521611e-02;
		d_L = NULL;
		d_L_ini_ = 4.80490089e-02;
		f_L = NULL;
		f_L_ini_ = 4.87798452e-01;
		d_T = NULL;
		d_T_ini_ = 4.20740474e-01;
		f_T = NULL;
		f_T_ini_ = 3.89684206e-02;
		q = NULL;
		q_ini_ = 2.97605397e-01;
		r = NULL;
		r_ini_ = 6.44029503e-02;
		P_af = NULL;
		P_af_ini_ = 1.30342012e-01;
		P_as = NULL;
		P_as_ini_ = 4.69609560e-01;
		P_i = NULL;
		P_i_ini_ = 8.79933753e-01;
		xs = NULL;
		xs_ini_ = 8.22938272e-02;
		y = NULL;
		y_ini_ = 3.88929176e-02;
		abstol__ = abs;
		reltol__ = rel;
		it_countx = 0;
	}
	Solveode::~Solveode()
	{
		if(V != NULL) free(V);
		if(m != NULL) free(m);
		if(h1 != NULL) free(h1);
		if(h2 != NULL) free(h2);
		if(d_L != NULL) free(d_L);
		if(f_L != NULL) free(f_L);
		if(d_T != NULL) free(d_T);
		if(f_T != NULL) free(f_T);
		if(q != NULL) free(q);
		if(r != NULL) free(r);
		if(P_af != NULL) free(P_af);
		if(P_as != NULL) free(P_as);
		if(P_i != NULL) free(P_i);
		if(xs != NULL) free(xs);
		if(y != NULL) free(y);
	}

	int Solveode::setVariables(int indVariable, double value_new)
	{
		switch(indVariable)
		{
		case 0:		V_ini_ = V_old_= value_new;    break;
		case 1:		m_ini_ = m_old_= value_new;    break;
		case 2:		h1_ini_ = h1_old_= value_new;    break;
		case 3:		h2_ini_ = h2_old_= value_new;    break;
		case 4:		d_L_ini_ = d_L_old_= value_new;    break;
		case 5:		f_L_ini_ = f_L_old_= value_new;    break;
		case 6:		d_T_ini_ = d_T_old_= value_new;    break;
		case 7:		f_T_ini_ = f_T_old_= value_new;    break;
		case 8:		q_ini_ = q_old_= value_new;    break;
		case 9:		r_ini_ = r_old_= value_new;    break;
		case 10:		P_af_ini_ = P_af_old_= value_new;    break;
		case 11:		P_as_ini_ = P_as_old_= value_new;    break;
		case 12:		P_i_ini_ = P_i_old_= value_new;    break;
		case 13:		xs_ini_ = xs_old_= value_new;    break;
		case 14:		y_ini_ = y_old_= value_new;    break;
		default:	return 1;    break;
		}
		return 0;
	}

	int Solveode::setParameters(int indVariable, double value_new)
	{
		switch(indVariable)
		{
		case 0:		time = value_new;   break;
		case 1:		dCell = value_new;   break;
		case 2:		Version = value_new;   break;
		case 3:		FCellConstant = value_new;   break;
		case 4:		CmCentre = value_new;   break;
		case 5:		CmPeriphery = value_new;   break;
		case 6:		g_Na_Centre_Published = value_new;   break;
		case 7:		g_Na_Periphery_Published = value_new;   break;
		case 8:		g_Na_Centre_0DCapable = value_new;   break;
		case 9:		g_Na_Periphery_0DCapable = value_new;   break;
		case 10:		g_Na_Centre_1DCapable = value_new;   break;
		case 11:		g_Na_Periphery_1DCapable = value_new;   break;
		case 12:		Na_o = value_new;   break;
		case 13:		F = value_new;   break;
		case 14:		R = value_new;   break;
		case 15:		T = value_new;   break;
		case 16:		g_Ca_L_Centre_Published = value_new;   break;
		case 17:		g_Ca_L_Periphery_Published = value_new;   break;
		case 18:		g_Ca_L_Centre_0DCapable = value_new;   break;
		case 19:		g_Ca_L_Periphery_0DCapable = value_new;   break;
		case 20:		g_Ca_L_Centre_1DCapable = value_new;   break;
		case 21:		g_Ca_L_Periphery_1DCapable = value_new;   break;
		case 22:		E_Ca_L = value_new;   break;
		case 23:		g_Ca_T_Centre_Published = value_new;   break;
		case 24:		g_Ca_T_Periphery_Published = value_new;   break;
		case 25:		g_Ca_T_Centre_0DCapable = value_new;   break;
		case 26:		g_Ca_T_Periphery_0DCapable = value_new;   break;
		case 27:		g_Ca_T_Centre_1DCapable = value_new;   break;
		case 28:		g_Ca_T_Periphery_1DCapable = value_new;   break;
		case 29:		E_Ca_T = value_new;   break;
		case 30:		g_to_Centre_Published = value_new;   break;
		case 31:		g_to_Periphery_Published = value_new;   break;
		case 32:		g_to_Centre_0DCapable = value_new;   break;
		case 33:		g_to_Periphery_0DCapable = value_new;   break;
		case 34:		g_to_Centre_1DCapable = value_new;   break;
		case 35:		g_to_Periphery_1DCapable = value_new;   break;
		case 36:		g_sus_Centre_Published = value_new;   break;
		case 37:		g_sus_Periphery_Published = value_new;   break;
		case 38:		g_sus_Centre_0DCapable = value_new;   break;
		case 39:		g_sus_Periphery_0DCapable = value_new;   break;
		case 40:		g_sus_Centre_1DCapable = value_new;   break;
		case 41:		g_sus_Periphery_1DCapable = value_new;   break;
		case 42:		g_K_r_Centre_Published = value_new;   break;
		case 43:		g_K_r_Periphery_Published = value_new;   break;
		case 44:		g_K_r_Centre_0DCapable = value_new;   break;
		case 45:		g_K_r_Periphery_0DCapable = value_new;   break;
		case 46:		g_K_r_Centre_1DCapable = value_new;   break;
		case 47:		g_K_r_Periphery_1DCapable = value_new;   break;
		case 48:		g_K_s_Centre_Published = value_new;   break;
		case 49:		g_K_s_Periphery_Published = value_new;   break;
		case 50:		g_K_s_Centre_0DCapable = value_new;   break;
		case 51:		g_K_s_Periphery_0DCapable = value_new;   break;
		case 52:		g_K_s_Centre_1DCapable = value_new;   break;
		case 53:		g_K_s_Periphery_1DCapable = value_new;   break;
		case 54:		g_f_Na_Centre_Published = value_new;   break;
		case 55:		g_f_Na_Periphery_Published = value_new;   break;
		case 56:		g_f_Na_Centre_0DCapable = value_new;   break;
		case 57:		g_f_Na_Periphery_0DCapable = value_new;   break;
		case 58:		g_f_Na_Centre_1DCapable = value_new;   break;
		case 59:		g_f_Na_Periphery_1DCapable = value_new;   break;
		case 60:		g_f_K_Centre_Published = value_new;   break;
		case 61:		g_f_K_Periphery_Published = value_new;   break;
		case 62:		g_f_K_Centre_0DCapable = value_new;   break;
		case 63:		g_f_K_Periphery_0DCapable = value_new;   break;
		case 64:		g_f_K_Centre_1DCapable = value_new;   break;
		case 65:		g_f_K_Periphery_1DCapable = value_new;   break;
		case 66:		g_b_Na_Centre_Published = value_new;   break;
		case 67:		g_b_Na_Periphery_Published = value_new;   break;
		case 68:		g_b_Na_Centre_0DCapable = value_new;   break;
		case 69:		g_b_Na_Periphery_0DCapable = value_new;   break;
		case 70:		g_b_Na_Centre_1DCapable = value_new;   break;
		case 71:		g_b_Na_Periphery_1DCapable = value_new;   break;
		case 72:		g_b_K_Centre_Published = value_new;   break;
		case 73:		g_b_K_Periphery_Published = value_new;   break;
		case 74:		g_b_K_Centre_0DCapable = value_new;   break;
		case 75:		g_b_K_Periphery_0DCapable = value_new;   break;
		case 76:		g_b_K_Centre_1DCapable = value_new;   break;
		case 77:		g_b_K_Periphery_1DCapable = value_new;   break;
		case 78:		g_b_Ca_Centre_Published = value_new;   break;
		case 79:		g_b_Ca_Periphery_Published = value_new;   break;
		case 80:		g_b_Ca_Centre_0DCapable = value_new;   break;
		case 81:		g_b_Ca_Periphery_0DCapable = value_new;   break;
		case 82:		g_b_Ca_Centre_1DCapable = value_new;   break;
		case 83:		g_b_Ca_Periphery_1DCapable = value_new;   break;
		case 84:		k_NaCa_Centre_Published = value_new;   break;
		case 85:		k_NaCa_Periphery_Published = value_new;   break;
		case 86:		k_NaCa_Centre_0DCapable = value_new;   break;
		case 87:		k_NaCa_Periphery_0DCapable = value_new;   break;
		case 88:		k_NaCa_Centre_1DCapable = value_new;   break;
		case 89:		k_NaCa_Periphery_1DCapable = value_new;   break;
		case 90:		Na_i = value_new;   break;
		case 91:		Ca_o = value_new;   break;
		case 92:		gamma_NaCa = value_new;   break;
		case 93:		Ca_i = value_new;   break;
		case 94:		d_NaCa = value_new;   break;
		case 95:		i_p_max_Centre_Published = value_new;   break;
		case 96:		i_p_max_Periphery_Published = value_new;   break;
		case 97:		i_p_max_Centre_0DCapable = value_new;   break;
		case 98:		i_p_max_Periphery_0DCapable = value_new;   break;
		case 99:		i_p_max_Centre_1DCapable = value_new;   break;
		case 100:		i_p_max_Periphery_1DCapable = value_new;   break;
		case 101:		K_m_Na = value_new;   break;
		case 102:		K_o = value_new;   break;
		case 103:		K_m_K = value_new;   break;
		case 104:		i_Ca_p_max_Centre_Published = value_new;   break;
		case 105:		i_Ca_p_max_Periphery_Published = value_new;   break;
		case 106:		i_Ca_p_max_Centre_0DCapable = value_new;   break;
		case 107:		i_Ca_p_max_Periphery_0DCapable = value_new;   break;
		case 108:		i_Ca_p_max_Centre_1DCapable = value_new;   break;
		case 109:		i_Ca_p_max_Periphery_1DCapable = value_new;   break;
		case 110:		K_i = value_new;   break;
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
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h1_new_ = this->h1_old_ = this->h1_ini_;
		this->h2_new_ = this->h2_old_ = this->h2_ini_;
		this->d_L_new_ = this->d_L_old_ = this->d_L_ini_;
		this->f_L_new_ = this->f_L_old_ = this->f_L_ini_;
		this->d_T_new_ = this->d_T_old_ = this->d_T_ini_;
		this->f_T_new_ = this->f_T_old_ = this->f_T_ini_;
		this->q_new_ = this->q_old_ = this->q_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->P_af_new_ = this->P_af_old_ = this->P_af_ini_;
		this->P_as_new_ = this->P_as_old_ = this->P_as_ini_;
		this->P_i_new_ = this->P_i_old_ = this->P_i_ini_;
		this->xs_new_ = this->xs_old_ = this->xs_ini_;
		this->y_new_ = this->y_old_ = this->y_ini_;
		this->time_new = this->time;
		this->timeSaving = this->time;
		if(savingRate!=0.0)
			this->save_step(fileptr, _EULER_);//save the initial conditions
		while(this->time_new<=finalTime){
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			this->V_new_ = this->dtime*(this->V_lado_direito_) + this->V_old_;
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h1_new_ = this->dtime*(this->h1_lado_direito_) + this->h1_old_;
			this->h2_new_ = this->dtime*(this->h2_lado_direito_) + this->h2_old_;
			this->d_L_new_ = this->dtime*(this->d_L_lado_direito_) + this->d_L_old_;
			this->f_L_new_ = this->dtime*(this->f_L_lado_direito_) + this->f_L_old_;
			this->d_T_new_ = this->dtime*(this->d_T_lado_direito_) + this->d_T_old_;
			this->f_T_new_ = this->dtime*(this->f_T_lado_direito_) + this->f_T_old_;
			this->q_new_ = this->dtime*(this->q_lado_direito_) + this->q_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->P_af_new_ = this->dtime*(this->P_af_lado_direito_) + this->P_af_old_;
			this->P_as_new_ = this->dtime*(this->P_as_lado_direito_) + this->P_as_old_;
			this->P_i_new_ = this->dtime*(this->P_i_lado_direito_) + this->P_i_old_;
			this->xs_new_ = this->dtime*(this->xs_lado_direito_) + this->xs_old_;
			this->y_new_ = this->dtime*(this->y_lado_direito_) + this->y_old_;
			//save results on a file
			if(savingRate!=0.0)
			{
				this->save_step(fileptr, _EULER_);
			}
		this->V_old_ = this->V_new_;
		this->m_old_ = this->m_new_;
		this->h1_old_ = this->h1_new_;
		this->h2_old_ = this->h2_new_;
		this->d_L_old_ = this->d_L_new_;
		this->f_L_old_ = this->f_L_new_;
		this->d_T_old_ = this->d_T_new_;
		this->f_T_old_ = this->f_T_new_;
		this->q_old_ = this->q_new_;
		this->r_old_ = this->r_new_;
		this->P_af_old_ = this->P_af_new_;
		this->P_as_old_ = this->P_as_new_;
		this->P_i_old_ = this->P_i_new_;
		this->xs_old_ = this->xs_new_;
		this->y_old_ = this->y_new_;
		}
	}
	void Solveode::rungeKutta2ndOrder(double finalTime, FILE *fileptr){
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h1_new_ = this->h1_old_ = this->h1_ini_;
		this->h2_new_ = this->h2_old_ = this->h2_ini_;
		this->d_L_new_ = this->d_L_old_ = this->d_L_ini_;
		this->f_L_new_ = this->f_L_old_ = this->f_L_ini_;
		this->d_T_new_ = this->d_T_old_ = this->d_T_ini_;
		this->f_T_new_ = this->f_T_old_ = this->f_T_ini_;
		this->q_new_ = this->q_old_ = this->q_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->P_af_new_ = this->P_af_old_ = this->P_af_ini_;
		this->P_as_new_ = this->P_as_old_ = this->P_as_ini_;
		this->P_i_new_ = this->P_i_old_ = this->P_i_ini_;
		this->xs_new_ = this->xs_old_ = this->xs_ini_;
		this->y_new_ = this->y_old_ = this->y_ini_;
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
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h1_new_ = this->dtime*(this->h1_lado_direito_) + this->h1_old_;
			this->h2_new_ = this->dtime*(this->h2_lado_direito_) + this->h2_old_;
			this->d_L_new_ = this->dtime*(this->d_L_lado_direito_) + this->d_L_old_;
			this->f_L_new_ = this->dtime*(this->f_L_lado_direito_) + this->f_L_old_;
			this->d_T_new_ = this->dtime*(this->d_T_lado_direito_) + this->d_T_old_;
			this->f_T_new_ = this->dtime*(this->f_T_lado_direito_) + this->f_T_old_;
			this->q_new_ = this->dtime*(this->q_lado_direito_) + this->q_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->P_af_new_ = this->dtime*(this->P_af_lado_direito_) + this->P_af_old_;
			this->P_as_new_ = this->dtime*(this->P_as_lado_direito_) + this->P_as_old_;
			this->P_i_new_ = this->dtime*(this->P_i_lado_direito_) + this->P_i_old_;
			this->xs_new_ = this->dtime*(this->xs_lado_direito_) + this->xs_old_;
			this->y_new_ = this->dtime*(this->y_lado_direito_) + this->y_old_;
			//stores the old variables in a vector
			for(int i=0;i<numEDO;i++){
				edos_old_aux_[i] = this->getVariables(i);
				edos_rightside_aux_[i] = this->getLadoDireito(i);
			}
			//steps one iteration ahead;
			this->V_old_ = this->V_new_;
			this->m_old_ = this->m_new_;
			this->h1_old_ = this->h1_new_;
			this->h2_old_ = this->h2_new_;
			this->d_L_old_ = this->d_L_new_;
			this->f_L_old_ = this->f_L_new_;
			this->d_T_old_ = this->d_T_new_;
			this->f_T_old_ = this->f_T_new_;
			this->q_old_ = this->q_new_;
			this->r_old_ = this->r_new_;
			this->P_af_old_ = this->P_af_new_;
			this->P_as_old_ = this->P_as_new_;
			this->P_i_old_ = this->P_i_new_;
			this->xs_old_ = this->xs_new_;
			this->y_old_ = this->y_new_;
			this->time_new	+= this->dtime;
			rightHandSideFunction.function(this);
			//computes the runge kutta second order method
			this->V_new_ = ( this->V_lado_direito_ + edos_rightside_aux_[0] ) * this->dtime/2 + edos_old_aux_[0];
			this->m_new_ = ( this->m_lado_direito_ + edos_rightside_aux_[1] ) * this->dtime/2 + edos_old_aux_[1];
			this->h1_new_ = ( this->h1_lado_direito_ + edos_rightside_aux_[2] ) * this->dtime/2 + edos_old_aux_[2];
			this->h2_new_ = ( this->h2_lado_direito_ + edos_rightside_aux_[3] ) * this->dtime/2 + edos_old_aux_[3];
			this->d_L_new_ = ( this->d_L_lado_direito_ + edos_rightside_aux_[4] ) * this->dtime/2 + edos_old_aux_[4];
			this->f_L_new_ = ( this->f_L_lado_direito_ + edos_rightside_aux_[5] ) * this->dtime/2 + edos_old_aux_[5];
			this->d_T_new_ = ( this->d_T_lado_direito_ + edos_rightside_aux_[6] ) * this->dtime/2 + edos_old_aux_[6];
			this->f_T_new_ = ( this->f_T_lado_direito_ + edos_rightside_aux_[7] ) * this->dtime/2 + edos_old_aux_[7];
			this->q_new_ = ( this->q_lado_direito_ + edos_rightside_aux_[8] ) * this->dtime/2 + edos_old_aux_[8];
			this->r_new_ = ( this->r_lado_direito_ + edos_rightside_aux_[9] ) * this->dtime/2 + edos_old_aux_[9];
			this->P_af_new_ = ( this->P_af_lado_direito_ + edos_rightside_aux_[10] ) * this->dtime/2 + edos_old_aux_[10];
			this->P_as_new_ = ( this->P_as_lado_direito_ + edos_rightside_aux_[11] ) * this->dtime/2 + edos_old_aux_[11];
			this->P_i_new_ = ( this->P_i_lado_direito_ + edos_rightside_aux_[12] ) * this->dtime/2 + edos_old_aux_[12];
			this->xs_new_ = ( this->xs_lado_direito_ + edos_rightside_aux_[13] ) * this->dtime/2 + edos_old_aux_[13];
			this->y_new_ = ( this->y_lado_direito_ + edos_rightside_aux_[14] ) * this->dtime/2 + edos_old_aux_[14];
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
			this->m_old_ = this->m_new_;
			this->h1_old_ = this->h1_new_;
			this->h2_old_ = this->h2_new_;
			this->d_L_old_ = this->d_L_new_;
			this->f_L_old_ = this->f_L_new_;
			this->d_T_old_ = this->d_T_new_;
			this->f_T_old_ = this->f_T_new_;
			this->q_old_ = this->q_new_;
			this->r_old_ = this->r_new_;
			this->P_af_old_ = this->P_af_new_;
			this->P_as_old_ = this->P_as_new_;
			this->P_i_old_ = this->P_i_new_;
			this->xs_old_ = this->xs_new_;
			this->y_old_ = this->y_new_;
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
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h1_new_ = this->h1_old_ = this->h1_ini_;
		this->h2_new_ = this->h2_old_ = this->h2_ini_;
		this->d_L_new_ = this->d_L_old_ = this->d_L_ini_;
		this->f_L_new_ = this->f_L_old_ = this->f_L_ini_;
		this->d_T_new_ = this->d_T_old_ = this->d_T_ini_;
		this->f_T_new_ = this->f_T_old_ = this->f_T_ini_;
		this->q_new_ = this->q_old_ = this->q_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->P_af_new_ = this->P_af_old_ = this->P_af_ini_;
		this->P_as_new_ = this->P_as_old_ = this->P_as_ini_;
		this->P_i_new_ = this->P_i_old_ = this->P_i_ini_;
		this->xs_new_ = this->xs_old_ = this->xs_ini_;
		this->y_new_ = this->y_old_ = this->y_ini_;
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
		int desc=0;
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
						this->m_new_ = edos_new_euler_[1];
						this->h1_new_ = edos_new_euler_[2];
						this->h2_new_ = edos_new_euler_[3];
						this->d_L_new_ = edos_new_euler_[4];
						this->f_L_new_ = edos_new_euler_[5];
						this->d_T_new_ = edos_new_euler_[6];
						this->f_T_new_ = edos_new_euler_[7];
						this->q_new_ = edos_new_euler_[8];
						this->r_new_ = edos_new_euler_[9];
						this->P_af_new_ = edos_new_euler_[10];
						this->P_as_new_ = edos_new_euler_[11];
						this->P_i_new_ = edos_new_euler_[12];
						this->xs_new_ = edos_new_euler_[13];
						this->y_new_ = edos_new_euler_[14];
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
		double maxDt = this->dtime, minDt = this->dtime;
		int desc=0;
		//initializes the variables
		this->previous_dt = this->dtime;
		//initializes the variables
		this->V_new_ = this->V_old_ = this->V_ini_;
		this->m_new_ = this->m_old_ = this->m_ini_;
		this->h1_new_ = this->h1_old_ = this->h1_ini_;
		this->h2_new_ = this->h2_old_ = this->h2_ini_;
		this->d_L_new_ = this->d_L_old_ = this->d_L_ini_;
		this->f_L_new_ = this->f_L_old_ = this->f_L_ini_;
		this->d_T_new_ = this->d_T_old_ = this->d_T_ini_;
		this->f_T_new_ = this->f_T_old_ = this->f_T_ini_;
		this->q_new_ = this->q_old_ = this->q_ini_;
		this->r_new_ = this->r_old_ = this->r_ini_;
		this->P_af_new_ = this->P_af_old_ = this->P_af_ini_;
		this->P_as_new_ = this->P_as_old_ = this->P_as_ini_;
		this->P_i_new_ = this->P_i_old_ = this->P_i_ini_;
		this->xs_new_ = this->xs_old_ = this->xs_ini_;
		this->y_new_ = this->y_old_ = this->y_ini_;
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
			    desc++;
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
					this->m_new_ = edos_new_euler_[1];
					this->h1_new_ = edos_new_euler_[2];
					this->h2_new_ = edos_new_euler_[3];
					this->d_L_new_ = edos_new_euler_[4];
					this->f_L_new_ = edos_new_euler_[5];
					this->d_T_new_ = edos_new_euler_[6];
					this->f_T_new_ = edos_new_euler_[7];
					this->q_new_ = edos_new_euler_[8];
					this->r_new_ = edos_new_euler_[9];
					this->P_af_new_ = edos_new_euler_[10];
					this->P_as_new_ = edos_new_euler_[11];
					this->P_i_new_ = edos_new_euler_[12];
					this->xs_new_ = edos_new_euler_[13];
					this->y_new_ = edos_new_euler_[14];
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
			double  _prvt_time = 0.0000000000e+00,  _prvt_dCell = 0.0000000000e+00,  _prvt_Version = 1.0000000000e+00,  _prvt_FCellConstant = 1.0309347000e+00,  _prvt_CmCentre = 2.0000000000e-05,  _prvt_CmPeriphery = 6.5000000000e-05,  _prvt_g_Na_Centre_Published = 0.0000000000e+00,  _prvt_g_Na_Periphery_Published = 1.2000000000e-06,  _prvt_g_Na_Centre_0DCapable = 0.0000000000e+00,  _prvt_g_Na_Periphery_0DCapable = 1.2040000000e-06,  _prvt_g_Na_Centre_1DCapable = 0.0000000000e+00,  _prvt_g_Na_Periphery_1DCapable = 3.7000000000e-07,  _prvt_Na_o = 1.4000000000e+02,  _prvt_F = 9.6845000000e+04,  _prvt_R = 8.3140000000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_g_Ca_L_Centre_Published = 5.8000000000e-03,  _prvt_g_Ca_L_Periphery_Published = 6.5900000000e-02,  _prvt_g_Ca_L_Centre_0DCapable = 5.7938000000e-03,  _prvt_g_Ca_L_Periphery_0DCapable = 6.5886480000e-02,  _prvt_g_Ca_L_Centre_1DCapable = 8.2000000000e-03,  _prvt_g_Ca_L_Periphery_1DCapable = 6.5900000000e-02,  _prvt_E_Ca_L = 4.6400000000e+01,  _prvt_g_Ca_T_Centre_Published = 4.3000000000e-03,  _prvt_g_Ca_T_Periphery_Published = 1.3900000000e-02,  _prvt_g_Ca_T_Centre_0DCapable = 4.2780600000e-03,  _prvt_g_Ca_T_Periphery_0DCapable = 1.3882300000e-02,  _prvt_g_Ca_T_Centre_1DCapable = 2.1000000000e-03,  _prvt_g_Ca_T_Periphery_1DCapable = 6.9400000000e-03,  _prvt_E_Ca_T = 4.5000000000e+01,  _prvt_g_to_Centre_Published = 4.9100000000e-03,  _prvt_g_to_Periphery_Published = 3.6490000000e-02,  _prvt_g_to_Centre_0DCapable = 4.9050000000e-03,  _prvt_g_to_Periphery_0DCapable = 3.6495000000e-02,  _prvt_g_to_Centre_1DCapable = 4.9050000000e-03,  _prvt_g_to_Periphery_1DCapable = 3.6500000000e-02,  _prvt_g_sus_Centre_Published = 6.6500000000e-05,  _prvt_g_sus_Periphery_Published = 1.1400000000e-02,  _prvt_g_sus_Centre_0DCapable = 6.6455040000e-05,  _prvt_g_sus_Periphery_0DCapable = 1.1383760000e-02,  _prvt_g_sus_Centre_1DCapable = 2.6600000000e-04,  _prvt_g_sus_Periphery_1DCapable = 1.1400000000e-02,  _prvt_g_K_r_Centre_Published = 7.9700000000e-04,  _prvt_g_K_r_Periphery_Published = 1.6000000000e-02,  _prvt_g_K_r_Centre_0DCapable = 7.9704000000e-04,  _prvt_g_K_r_Periphery_0DCapable = 1.6000000000e-02,  _prvt_g_K_r_Centre_1DCapable = 7.3800000000e-04,  _prvt_g_K_r_Periphery_1DCapable = 2.0800000000e-02,  _prvt_g_K_s_Centre_Published = 5.1800000000e-04,  _prvt_g_K_s_Periphery_Published = 1.0400000000e-02,  _prvt_g_K_s_Centre_0DCapable = 3.4450000000e-04,  _prvt_g_K_s_Periphery_0DCapable = 1.0400000000e-02,  _prvt_g_K_s_Centre_1DCapable = 3.4500000000e-04,  _prvt_g_K_s_Periphery_1DCapable = 1.0400000000e-02,  _prvt_g_f_Na_Centre_Published = 5.4800000000e-04,  _prvt_g_f_Na_Periphery_Published = 6.9000000000e-03,  _prvt_g_f_Na_Centre_0DCapable = 5.4650000000e-04,  _prvt_g_f_Na_Periphery_0DCapable = 6.8750000000e-03,  _prvt_g_f_Na_Centre_1DCapable = 4.3700000000e-04,  _prvt_g_f_Na_Periphery_1DCapable = 5.5000000000e-03,  _prvt_g_f_K_Centre_Published = 5.4800000000e-04,  _prvt_g_f_K_Periphery_Published = 6.9000000000e-03,  _prvt_g_f_K_Centre_0DCapable = 5.4650000000e-04,  _prvt_g_f_K_Periphery_0DCapable = 6.8750000000e-03,  _prvt_g_f_K_Centre_1DCapable = 4.3700000000e-04,  _prvt_g_f_K_Periphery_1DCapable = 5.5000000000e-03,  _prvt_g_b_Na_Centre_Published = 5.8000000000e-05,  _prvt_g_b_Na_Periphery_Published = 1.8900000000e-04,  _prvt_g_b_Na_Centre_0DCapable = 5.8181800000e-05,  _prvt_g_b_Na_Periphery_0DCapable = 1.8880000000e-04,  _prvt_g_b_Na_Centre_1DCapable = 5.8000000000e-05,  _prvt_g_b_Na_Periphery_1DCapable = 1.8900000000e-04,  _prvt_g_b_K_Centre_Published = 2.5200000000e-05,  _prvt_g_b_K_Periphery_Published = 8.1900000000e-05,  _prvt_g_b_K_Centre_0DCapable = 2.5236360000e-05,  _prvt_g_b_K_Periphery_0DCapable = 8.1892000000e-05,  _prvt_g_b_K_Centre_1DCapable = 2.5200000000e-05,  _prvt_g_b_K_Periphery_1DCapable = 8.1900000000e-05,  _prvt_g_b_Ca_Centre_Published = 1.3200000000e-05,  _prvt_g_b_Ca_Periphery_Published = 4.3000000000e-05,  _prvt_g_b_Ca_Centre_0DCapable = 1.3236000000e-05,  _prvt_g_b_Ca_Periphery_0DCapable = 4.2952000000e-05,  _prvt_g_b_Ca_Centre_1DCapable = 1.3230000000e-05,  _prvt_g_b_Ca_Periphery_1DCapable = 4.2900000000e-05,  _prvt_k_NaCa_Centre_Published = 2.7000000000e-06,  _prvt_k_NaCa_Periphery_Published = 8.8000000000e-06,  _prvt_k_NaCa_Centre_0DCapable = 2.7229000000e-06,  _prvt_k_NaCa_Periphery_0DCapable = 8.8358400000e-06,  _prvt_k_NaCa_Centre_1DCapable = 2.8000000000e-06,  _prvt_k_NaCa_Periphery_1DCapable = 8.8000000000e-06,  _prvt_Na_i = 8.0000000000e+00,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_gamma_NaCa = 5.0000000000e-01,  _prvt_Ca_i = 1.0000000000e-04,  _prvt_d_NaCa = 1.0000000000e-04,  _prvt_i_p_max_Centre_Published = 4.7800000000e-02,  _prvt_i_p_max_Periphery_Published = 1.6000000000e-01,  _prvt_i_p_max_Centre_0DCapable = 4.7825450000e-02,  _prvt_i_p_max_Periphery_0DCapable = 1.5519360000e-01,  _prvt_i_p_max_Centre_1DCapable = 4.7800000000e-02,  _prvt_i_p_max_Periphery_1DCapable = 1.6000000000e-01,  _prvt_K_m_Na = 5.6400000000e+00,  _prvt_K_o = 5.4000000000e+00,  _prvt_K_m_K = 6.2100000000e-01,  _prvt_i_Ca_p_max_Centre_Published = 0.0000000000e+00,  _prvt_i_Ca_p_max_Periphery_Published = 0.0000000000e+00,  _prvt_i_Ca_p_max_Centre_0DCapable = 0.0000000000e+00,  _prvt_i_Ca_p_max_Periphery_0DCapable = 0.0000000000e+00,  _prvt_i_Ca_p_max_Centre_1DCapable = 4.2000000000e-03,  _prvt_i_Ca_p_max_Periphery_1DCapable = 3.3390000000e-02,  _prvt_K_i = 1.4000000000e+02, 
			//private aux variables
			 _prvt_calc_FCell=0,  _prvt_calc_Cm=0,  _prvt_calc_g_Na=0,  _prvt_calc_i_Na=0,  _prvt_calc_m_infinity=0,  _prvt_calc_tau_m=0,  _prvt_calc_F_Na=0,  _prvt_calc_h=0,  _prvt_calc_h1_infinity=0,  _prvt_calc_h2_infinity=0,  _prvt_calc_tau_h1=0,  _prvt_calc_tau_h2=0,  _prvt_calc_g_Ca_L=0,  _prvt_calc_i_Ca_L=0,  _prvt_calc_alpha_d_L=0,  _prvt_calc_beta_d_L=0,  _prvt_calc_tau_d_L=0,  _prvt_calc_d_L_infinity=0,  _prvt_calc_alpha_f_L=0,  _prvt_calc_beta_f_L=0,  _prvt_calc_tau_f_L=0,  _prvt_calc_f_L_infinity=0,  _prvt_calc_g_Ca_T=0,  _prvt_calc_i_Ca_T=0,  _prvt_calc_alpha_d_T=0,  _prvt_calc_beta_d_T=0,  _prvt_calc_tau_d_T=0,  _prvt_calc_d_T_infinity=0,  _prvt_calc_alpha_f_T=0,  _prvt_calc_beta_f_T=0,  _prvt_calc_tau_f_T=0,  _prvt_calc_f_T_infinity=0,  _prvt_calc_g_to=0,  _prvt_calc_g_sus=0,  _prvt_calc_i_to=0,  _prvt_calc_i_sus=0,  _prvt_calc_q_infinity=0,  _prvt_calc_tau_q=0,  _prvt_calc_r_infinity=0,  _prvt_calc_tau_r=0,  _prvt_calc_g_K_r=0,  _prvt_calc_i_K_r=0,  _prvt_calc_P_a=0,  _prvt_calc_P_af_infinity=0,  _prvt_calc_tau_P_af=0,  _prvt_calc_P_as_infinity=0,  _prvt_calc_tau_P_as=0,  _prvt_calc_tau_P_i=0,  _prvt_calc_P_i_infinity=0,  _prvt_calc_g_K_s=0,  _prvt_calc_i_K_s=0,  _prvt_calc_alpha_xs=0,  _prvt_calc_beta_xs=0,  _prvt_calc_g_f_Na=0,  _prvt_calc_i_f_Na=0,  _prvt_calc_g_f_K=0,  _prvt_calc_i_f_K=0,  _prvt_calc_alpha_y=0,  _prvt_calc_beta_y=0,  _prvt_calc_g_b_Na=0,  _prvt_calc_i_b_Na=0,  _prvt_calc_g_b_K=0,  _prvt_calc_i_b_K=0,  _prvt_calc_g_b_Ca=0,  _prvt_calc_i_b_Ca=0,  _prvt_calc_k_NaCa=0,  _prvt_calc_i_NaCa=0,  _prvt_calc_i_p_max=0,  _prvt_calc_i_p=0,  _prvt_calc_i_Ca_p_max=0,  _prvt_calc_i_Ca_p=0,  _prvt_calc_E_Na=0,  _prvt_calc_E_K=0,  _prvt_calc_E_Ca=0,  _prvt_calc_E_K_s=0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h1_lado_direito_,  _prvt_h2_lado_direito_,  _prvt_d_L_lado_direito_,  _prvt_f_L_lado_direito_,  _prvt_d_T_lado_direito_,  _prvt_f_T_lado_direito_,  _prvt_q_lado_direito_,  _prvt_r_lado_direito_,  _prvt_P_af_lado_direito_,  _prvt_P_as_lado_direito_,  _prvt_P_i_lado_direito_,  _prvt_xs_lado_direito_,  _prvt_y_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate;
			__NEW_[0] = __OLD_[0] = -3.9013558536e+01;
			__NEW_[1] = __OLD_[1] = 9.2361701692e-02;
			__NEW_[2] = __OLD_[2] = 1.5905380261e-02;
			__NEW_[3] = __OLD_[3] = 1.4452161090e-02;
			__NEW_[4] = __OLD_[4] = 4.8049008950e-02;
			__NEW_[5] = __OLD_[5] = 4.8779845203e-01;
			__NEW_[6] = __OLD_[6] = 4.2074047435e-01;
			__NEW_[7] = __OLD_[7] = 3.8968420558e-02;
			__NEW_[8] = __OLD_[8] = 2.9760539675e-01;
			__NEW_[9] = __OLD_[9] = 6.4402950262e-02;
			__NEW_[10] = __OLD_[10] = 1.3034201158e-01;
			__NEW_[11] = __OLD_[11] = 4.6960956028e-01;
			__NEW_[12] = __OLD_[12] = 8.7993375273e-01;
			__NEW_[13] = __OLD_[13] = 8.2293827208e-02;
			__NEW_[14] = __OLD_[14] = 3.8892917590e-02;
			int *_prvt_tree_thread = tree_thread;
			while(_prvt_time_new<=_prvt_finalTime){
				_prvt_time_new += _prvt_dtime;
				if(omp_get_thread_num()==tree_thread[0])
				{
					_prvt_calc_FCell = ((_prvt_Version==0.0000000000e+00))
?(((1.0700000000e+00*((3.0000000000e+00*_prvt_dCell)-1.0000000000e-01))/(3.0000000000e+00*(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01)))))))
:(((_prvt_Version==1.0000000000e+00))
?(((_prvt_FCellConstant*_prvt_dCell)/(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01))))))
:(((1.0700000000e+00*2.9000000000e+01*_prvt_dCell)/(3.0000000000e+01*(1.0000000000e+00+(7.7450000000e-01*exp(((-((2.9000000000e+01*_prvt_dCell)-2.4500000000e+01))/1.9500000000e+00))))))))
;
					_prvt_calc_F_Na = ((_prvt_Version==0.0000000000e+00))
?((((9.5200000000e-02*exp(((-6.3000000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6600000000e+00*exp(((-2.2500000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6900000000e-02))
:((((9.5180000000e-02*exp(((-6.3060000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6620000000e+00*exp(((-2.2510000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6930000000e-02));
					_prvt_calc_alpha_f_L = ((_prvt_Version==1.0000000000e+00))
?(((3.7500000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)))
:(((3.1200000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
					_prvt_calc_beta_f_L = ((_prvt_Version==1.0000000000e+00))
?((3.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))))
:((2.5000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))));
					_prvt_calc_f_L_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5000000000e+01)/5.0000000000e+00))));
					_prvt_calc_beta_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5000000000e+01*exp(((__OLD_[0]+7.1000000000e+01)/1.5380000000e+01))))
:((1.5000000000e+01*exp(((__OLD_[0]+7.1700000000e+01)/1.5380000000e+01))));
					_prvt_calc_f_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+7.1000000000e+01)/9.0000000000e+00))));
					_prvt_calc_P_a = ((6.0000000000e-01*__OLD_[10])+(4.0000000000e-01*__OLD_[11]));
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/_prvt_Na_i)));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/_prvt_K_i)));
					_prvt_calc_E_Ca = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Ca_o/_prvt_Ca_i)));
					_prvt_calc_E_K_s = ((_prvt_Version==0.0000000000e+00))
?((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(1.2000000000e-01*_prvt_Na_o))/(_prvt_K_i+(1.2000000000e-01*_prvt_Na_i))))))
:((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(3.0000000000e-02*_prvt_Na_o))/(_prvt_K_i+(3.0000000000e-02*_prvt_Na_i))))));
					_prvt_calc_Cm = (_prvt_CmCentre+(_prvt_calc_FCell*(_prvt_CmPeriphery-_prvt_CmCentre)));
					_prvt_calc_g_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_Published-_prvt_g_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_0DCapable-_prvt_g_Na_Centre_0DCapable))))
:((_prvt_g_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_1DCapable-_prvt_g_Na_Centre_1DCapable)))))
;
					_prvt_calc_h = (((1.0000000000e+00-_prvt_calc_F_Na)*__OLD_[2])+(_prvt_calc_F_Na*__OLD_[3]));
					_prvt_calc_g_Ca_L = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_L_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_Published-_prvt_g_Ca_L_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_L_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_0DCapable-_prvt_g_Ca_L_Centre_0DCapable))))
:((_prvt_g_Ca_L_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_1DCapable-_prvt_g_Ca_L_Centre_1DCapable)))))
;
					_prvt_calc_tau_f_L = ((_prvt_Version==1.0000000000e+00))
?(((1.2000000000e+00-(2.0000000000e-01*_prvt_calc_FCell))/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)))
:((1.0000000000e+00/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)));
					_prvt_calc_g_Ca_T = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_T_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_Published-_prvt_g_Ca_T_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_T_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_0DCapable-_prvt_g_Ca_T_Centre_0DCapable))))
:((_prvt_g_Ca_T_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_1DCapable-_prvt_g_Ca_T_Centre_1DCapable)))))
;
					_prvt_calc_alpha_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5300000000e+01*exp(((-(__OLD_[0]+7.1000000000e+01+(7.0000000000e-01*_prvt_calc_FCell)))/8.3300000000e+01))))
:((1.5300000000e+01*exp(((-(__OLD_[0]+7.1700000000e+01))/8.3300000000e+01))));
					_prvt_calc_g_to = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_to_Centre_Published+(_prvt_calc_FCell*(_prvt_g_to_Periphery_Published-_prvt_g_to_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_to_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_0DCapable-_prvt_g_to_Centre_0DCapable))))
:((_prvt_g_to_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_1DCapable-_prvt_g_to_Centre_1DCapable)))))
;
					_prvt_calc_g_sus = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_sus_Centre_Published+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_Published-_prvt_g_sus_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_sus_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_0DCapable-_prvt_g_sus_Centre_0DCapable))))
:((_prvt_g_sus_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_1DCapable-_prvt_g_sus_Centre_1DCapable)))))
;
					_prvt_calc_g_K_r = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_r_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_Published-_prvt_g_K_r_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_r_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_0DCapable-_prvt_g_K_r_Centre_0DCapable))))
:((_prvt_g_K_r_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_1DCapable-_prvt_g_K_r_Centre_1DCapable)))))
;
					_prvt_calc_g_K_s = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_s_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_Published-_prvt_g_K_s_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_s_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_0DCapable-_prvt_g_K_s_Centre_0DCapable))))
:((_prvt_g_K_s_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_1DCapable-_prvt_g_K_s_Centre_1DCapable)))))
;
					_prvt_calc_g_f_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_Published-_prvt_g_f_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_0DCapable-_prvt_g_f_Na_Centre_0DCapable))))
:((_prvt_g_f_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_1DCapable-_prvt_g_f_Na_Centre_1DCapable)))))
;
					_prvt_calc_g_f_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_Published-_prvt_g_f_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_0DCapable-_prvt_g_f_K_Centre_0DCapable))))
:((_prvt_g_f_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_1DCapable-_prvt_g_f_K_Centre_1DCapable)))))
;
					_prvt_calc_g_b_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_Published-_prvt_g_b_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_0DCapable-_prvt_g_b_Na_Centre_0DCapable))))
:((_prvt_g_b_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_1DCapable-_prvt_g_b_Na_Centre_1DCapable)))))
;
					_prvt_calc_g_b_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_Published-_prvt_g_b_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_0DCapable-_prvt_g_b_K_Centre_0DCapable))))
:((_prvt_g_b_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_1DCapable-_prvt_g_b_K_Centre_1DCapable)))))
;
					_prvt_calc_g_b_Ca = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Ca_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_Published-_prvt_g_b_Ca_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Ca_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_0DCapable-_prvt_g_b_Ca_Centre_0DCapable))))
:((_prvt_g_b_Ca_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_1DCapable-_prvt_g_b_Ca_Centre_1DCapable)))))
;
					_prvt_calc_k_NaCa = ((_prvt_Version==0.0000000000e+00))
?((_prvt_k_NaCa_Centre_Published+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_Published-_prvt_k_NaCa_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_k_NaCa_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_0DCapable-_prvt_k_NaCa_Centre_0DCapable))))
:((_prvt_k_NaCa_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_1DCapable-_prvt_k_NaCa_Centre_1DCapable)))))
;
					_prvt_calc_i_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_Published-_prvt_i_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_0DCapable-_prvt_i_p_max_Centre_0DCapable))))
:((_prvt_i_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_1DCapable-_prvt_i_p_max_Centre_1DCapable)))))
;
					_prvt_calc_i_Ca_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_Published-_prvt_i_Ca_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_0DCapable-_prvt_i_Ca_p_max_Centre_0DCapable))))
:((_prvt_i_Ca_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_1DCapable-_prvt_i_Ca_p_max_Centre_1DCapable)))))
;
					_prvt_calc_i_Ca_L = (_prvt_calc_g_Ca_L*((__OLD_[5]*__OLD_[4])+(6.0000000000e-03/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4100000000e+01))/6.0000000000e+00)))))*(__OLD_[0]-_prvt_E_Ca_L));
					_prvt_calc_i_Ca_T = (_prvt_calc_g_Ca_T*__OLD_[6]*__OLD_[7]*(__OLD_[0]-_prvt_E_Ca_T));
					_prvt_calc_tau_f_T = (1.0000000000e+00/(_prvt_calc_alpha_f_T+_prvt_calc_beta_f_T));
					_prvt_calc_i_to = (_prvt_calc_g_to*__OLD_[8]*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_sus = (_prvt_calc_g_sus*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K_r = (_prvt_calc_g_K_r*_prvt_calc_P_a*__OLD_[12]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K_s = (_prvt_calc_g_K_s*pow(__OLD_[13],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K_s));
					_prvt_calc_i_f_Na = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_Na)))
:((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-7.7600000000e+01)));
					_prvt_calc_i_f_K = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_K)))
:((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]+1.0200000000e+02)));
					_prvt_calc_i_b_Na = (_prvt_calc_g_b_Na*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_K = (_prvt_calc_g_b_K*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_b_Ca = (_prvt_calc_g_b_Ca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_NaCa = ((_prvt_Version==0.0000000000e+00))
?(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7400000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))))
:(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7430000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))));
					_prvt_calc_i_p = ((_prvt_calc_i_p_max*pow((_prvt_Na_i/(_prvt_K_m_Na+_prvt_Na_i)),3.0000000000e+00)*pow((_prvt_K_o/(_prvt_K_m_K+_prvt_K_o)),2.0000000000e+00)*1.6000000000e+00)/(1.5000000000e+00+exp(((-(__OLD_[0]+6.0000000000e+01))/4.0000000000e+01))));
					_prvt_calc_i_Ca_p = ((_prvt_calc_i_Ca_p_max*_prvt_Ca_i)/(_prvt_Ca_i+4.0000000000e-04));
					_prvt_calc_i_Na = (((((_prvt_calc_g_Na*pow(__OLD_[1],3.0000000000e+00)*_prvt_calc_h*_prvt_Na_o*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*(exp((((__OLD_[0]-_prvt_calc_E_Na)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))/(exp(((__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))*__OLD_[0]);
					_prvt_V_lado_direito_= (((-1.0000000000e+00)/_prvt_calc_Cm)*(_prvt_calc_i_Na+_prvt_calc_i_Ca_L+_prvt_calc_i_Ca_T+_prvt_calc_i_to+_prvt_calc_i_sus+_prvt_calc_i_K_r+_prvt_calc_i_K_s+_prvt_calc_i_f_Na+_prvt_calc_i_f_K+_prvt_calc_i_b_Na+_prvt_calc_i_b_Ca+_prvt_calc_i_b_K+_prvt_calc_i_NaCa+_prvt_calc_i_p+_prvt_calc_i_Ca_p));
					__NEW_[0]= _prvt_V_lado_direito_ * _prvt_dtime + __OLD_[0];
					_prvt_f_L_lado_direito_= ((_prvt_calc_f_L_infinity-__OLD_[5])/_prvt_calc_tau_f_L);
					__NEW_[5]= _prvt_f_L_lado_direito_ * _prvt_dtime + __OLD_[5];
					_prvt_f_T_lado_direito_= ((_prvt_calc_f_T_infinity-__OLD_[7])/_prvt_calc_tau_f_T);
					__NEW_[7]= _prvt_f_T_lado_direito_ * _prvt_dtime + __OLD_[7];
				}
				if(omp_get_thread_num()==tree_thread[1])
				{
					_prvt_calc_m_infinity = ((_prvt_Version==0.0000000000e+00))
?(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-__OLD_[0])/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)))
:(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.0320000000e+01))/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
					_prvt_calc_tau_m = ((_prvt_Version==0.0000000000e+00))
?(((6.2470000000e-04/((8.3200000000e-01*exp(((-3.3500000000e-01)*(__OLD_[0]+5.6700000000e+01))))+(6.2700000000e-01*exp((8.2000000000e-02*(__OLD_[0]+6.5010000000e+01))))))+4.0000000000e-05))
:(((6.2470000000e-04/((8.3221660000e-01*exp(((-3.3566000000e-01)*(__OLD_[0]+5.6706200000e+01))))+(6.2740000000e-01*exp((8.2300000000e-02*(__OLD_[0]+6.5013100000e+01))))))+4.5690000000e-05));
					_prvt_m_lado_direito_= ((_prvt_calc_m_infinity-__OLD_[1])/_prvt_calc_tau_m);
					__NEW_[1]= _prvt_m_lado_direito_ * _prvt_dtime + __OLD_[1];
				}
				if(omp_get_thread_num()==tree_thread[2])
				{
					_prvt_calc_h1_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+6.6100000000e+01)/6.4000000000e+00))));
					_prvt_calc_tau_h1 = (((3.7170000000e-06*exp(((-2.8150000000e-01)*(__OLD_[0]+1.7110000000e+01))))/(1.0000000000e+00+(3.7320000000e-03*exp(((-3.4260000000e-01)*(__OLD_[0]+3.7760000000e+01))))))+5.9770000000e-04);
					_prvt_calc_tau_h2 = (((3.1860000000e-08*exp(((-6.2190000000e-01)*(__OLD_[0]+1.8800000000e+01))))/(1.0000000000e+00+(7.1890000000e-05*exp(((-6.6830000000e-01)*(__OLD_[0]+3.4070000000e+01))))))+3.5560000000e-03);
					_prvt_calc_h2_infinity = _prvt_calc_h1_infinity;
					_prvt_h1_lado_direito_= ((_prvt_calc_h1_infinity-__OLD_[2])/_prvt_calc_tau_h1);
					__NEW_[2]= _prvt_h1_lado_direito_ * _prvt_dtime + __OLD_[2];
					_prvt_h2_lado_direito_= ((_prvt_calc_h2_infinity-__OLD_[3])/_prvt_calc_tau_h2);
					__NEW_[3]= _prvt_h2_lado_direito_ * _prvt_dtime + __OLD_[3];
				}
				if(omp_get_thread_num()==tree_thread[3])
				{
					_prvt_calc_alpha_d_L = ((_prvt_Version==0.0000000000e+00))
?(((((-2.8380000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((_prvt_Version==1.0000000000e+00))
?(((((-2.8390000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((((-2.8400000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00)))))
;
					_prvt_calc_beta_d_L = ((_prvt_Version==1.0000000000e+00))
?(((1.1430000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)))
:(((1.1420000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)));
					_prvt_calc_d_L_infinity = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.3100000000e+01))/6.0000000000e+00)))))
:(((_prvt_Version==1.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2300000000e+01+(8.0000000000e-01*_prvt_calc_FCell)))/6.0000000000e+00)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2200000000e+01))/6.0000000000e+00))))))
;
					_prvt_calc_tau_d_L = (2.0000000000e+00/(_prvt_calc_alpha_d_L+_prvt_calc_beta_d_L));
					_prvt_d_L_lado_direito_= ((_prvt_calc_d_L_infinity-__OLD_[4])/_prvt_calc_tau_d_L);
					__NEW_[4]= _prvt_d_L_lado_direito_ * _prvt_dtime + __OLD_[4];
				}
				if(omp_get_thread_num()==tree_thread[4])
				{
					_prvt_calc_alpha_d_T = (1.0680000000e+03*exp(((__OLD_[0]+2.6300000000e+01)/3.0000000000e+01)));
					_prvt_calc_beta_d_T = (1.0680000000e+03*exp(((-(__OLD_[0]+2.6300000000e+01))/3.0000000000e+01)));
					_prvt_calc_d_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.7000000000e+01))/6.8000000000e+00))));
					_prvt_calc_tau_d_T = (1.0000000000e+00/(_prvt_calc_alpha_d_T+_prvt_calc_beta_d_T));
					_prvt_d_T_lado_direito_= ((_prvt_calc_d_T_infinity-__OLD_[6])/_prvt_calc_tau_d_T);
					__NEW_[6]= _prvt_d_T_lado_direito_ * _prvt_dtime + __OLD_[6];
				}
				if(omp_get_thread_num()==tree_thread[5])
				{
					_prvt_calc_q_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.9370000000e+01)/1.3100000000e+01))));
					_prvt_calc_tau_q = ((_prvt_Version==0.0000000000e+00))
?((1.0100000000e-02+(6.5170000000e-02/(5.7000000000e-01*exp(((-8.0000000000e-02)*(__OLD_[0]+4.9000000000e+01)))))+(2.4000000000e-05*exp((1.0000000000e-01*(__OLD_[0]+5.0930000000e+01))))))
:(((_prvt_Version==1.0000000000e+00))
?(((1.0000000000e-03/3.0000000000e+00)*(3.0310000000e+01+(1.9550000000e+02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01+(1.0000000000e+01*_prvt_calc_FCell)))))+(7.1740000000e-01*exp(((2.7190000000e-01-(1.7190000000e-01*_prvt_calc_FCell))*1.0000000000e+00*(__OLD_[0]+4.0930000000e+01+(1.0000000000e+01*_prvt_calc_FCell))))))))))
:((1.0100000000e-02+(6.5170000000e-02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01))))+(7.1740000000e-01*exp((2.7190000000e-01*(__OLD_[0]+4.0930000000e+01)))))))))
;
					_prvt_q_lado_direito_= ((_prvt_calc_q_infinity-__OLD_[8])/_prvt_calc_tau_q);
					__NEW_[8]= _prvt_q_lado_direito_ * _prvt_dtime + __OLD_[8];
				}
				if(omp_get_thread_num()==tree_thread[6])
				{
					_prvt_calc_r_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]-1.0930000000e+01))/1.9700000000e+01))));
					_prvt_calc_tau_r = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e-03*(2.9800000000e+00+(1.5590000000e+01/((1.0370000000e+00*exp((9.0000000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.2000000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:(((_prvt_Version==1.0000000000e+00))
?((2.5000000000e-03*(1.1910000000e+00+(7.8380000000e+00/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:((1.0000000000e-03*(2.9800000000e+00+(1.9590000000e+01/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01))))))))))
;
					_prvt_r_lado_direito_= ((_prvt_calc_r_infinity-__OLD_[9])/_prvt_calc_tau_r);
					__NEW_[9]= _prvt_r_lado_direito_ * _prvt_dtime + __OLD_[9];
				}
				if(omp_get_thread_num()==tree_thread[7])
				{
					_prvt_calc_P_af_infinity = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4200000000e+01))/1.0600000000e+01)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.3200000000e+01))/1.0600000000e+01)))));
					_prvt_calc_tau_P_af = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-9.0000000000e+00)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.2500000000e+01))))))
:((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-1.0000000000e+01)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.2500000000e+01))))));
					_prvt_calc_tau_P_as = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-9.0000000000e+00)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.1600000000e+01))))))
:((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-1.0000000000e+01)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.1600000000e+01))))));
					_prvt_calc_P_as_infinity = _prvt_calc_P_af_infinity;
					_prvt_P_af_lado_direito_= ((_prvt_calc_P_af_infinity-__OLD_[10])/_prvt_calc_tau_P_af);
					__NEW_[10]= _prvt_P_af_lado_direito_ * _prvt_dtime + __OLD_[10];
					_prvt_P_as_lado_direito_= ((_prvt_calc_P_as_infinity-__OLD_[11])/_prvt_calc_tau_P_as);
					__NEW_[11]= _prvt_P_as_lado_direito_ * _prvt_dtime + __OLD_[11];
				}
				if(omp_get_thread_num()==tree_thread[8])
				{
					_prvt_calc_tau_P_i = ((_prvt_Version==0.0000000000e+00))
?(2.0000000000e-03)
:(((_prvt_Version==1.0000000000e+00))
?(2.0000000000e-03)
:(6.0000000000e-03))
;
					_prvt_calc_P_i_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+1.8600000000e+01)/1.0100000000e+01))));
					_prvt_P_i_lado_direito_= ((_prvt_calc_P_i_infinity-__OLD_[12])/_prvt_calc_tau_P_i);
					__NEW_[12]= _prvt_P_i_lado_direito_ * _prvt_dtime + __OLD_[12];
				}
				if(omp_get_thread_num()==tree_thread[9])
				{
					_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
					_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
					_prvt_xs_lado_direito_= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[13]))-(_prvt_calc_beta_xs*__OLD_[13]));
					__NEW_[13]= _prvt_xs_lado_direito_ * _prvt_dtime + __OLD_[13];
				}
				if(omp_get_thread_num()==tree_thread[10])
				{
					_prvt_calc_alpha_y = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6620000000e+01))))
:((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6630000000e+01))));
					_prvt_calc_beta_y = (1.0000000000e+00*exp(((__OLD_[0]+7.5130000000e+01)/2.1250000000e+01)));
					_prvt_y_lado_direito_= ((_prvt_calc_alpha_y*(1.0000000000e+00-__OLD_[14]))-(_prvt_calc_beta_y*__OLD_[14]));
					__NEW_[14]= _prvt_y_lado_direito_ * _prvt_dtime + __OLD_[14];
				}
				//synchronizing all threads
				#pragma omp barrier
				/*if(_prvt_savingRate!=0){
					#pragma omp single
					{
						this->V_old_ = __OLD_[0];
						this->V_new_ = __NEW_[0];
						this->m_old_ = __OLD_[1];
						this->m_new_ = __NEW_[1];
						this->h1_old_ = __OLD_[2];
						this->h1_new_ = __NEW_[2];
						this->h2_old_ = __OLD_[3];
						this->h2_new_ = __NEW_[3];
						this->d_L_old_ = __OLD_[4];
						this->d_L_new_ = __NEW_[4];
						this->f_L_old_ = __OLD_[5];
						this->f_L_new_ = __NEW_[5];
						this->d_T_old_ = __OLD_[6];
						this->d_T_new_ = __NEW_[6];
						this->f_T_old_ = __OLD_[7];
						this->f_T_new_ = __NEW_[7];
						this->q_old_ = __OLD_[8];
						this->q_new_ = __NEW_[8];
						this->r_old_ = __OLD_[9];
						this->r_new_ = __NEW_[9];
						this->P_af_old_ = __OLD_[10];
						this->P_af_new_ = __NEW_[10];
						this->P_as_old_ = __OLD_[11];
						this->P_as_new_ = __NEW_[11];
						this->P_i_old_ = __OLD_[12];
						this->P_i_new_ = __NEW_[12];
						this->xs_old_ = __OLD_[13];
						this->xs_new_ = __NEW_[13];
						this->y_old_ = __OLD_[14];
						this->y_new_ = __NEW_[14];
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
			double  _prvt_time = 0.0000000000e+00,  _prvt_dCell = 0.0000000000e+00,  _prvt_Version = 1.0000000000e+00,  _prvt_FCellConstant = 1.0309347000e+00,  _prvt_CmCentre = 2.0000000000e-05,  _prvt_CmPeriphery = 6.5000000000e-05,  _prvt_g_Na_Centre_Published = 0.0000000000e+00,  _prvt_g_Na_Periphery_Published = 1.2000000000e-06,  _prvt_g_Na_Centre_0DCapable = 0.0000000000e+00,  _prvt_g_Na_Periphery_0DCapable = 1.2040000000e-06,  _prvt_g_Na_Centre_1DCapable = 0.0000000000e+00,  _prvt_g_Na_Periphery_1DCapable = 3.7000000000e-07,  _prvt_Na_o = 1.4000000000e+02,  _prvt_F = 9.6845000000e+04,  _prvt_R = 8.3140000000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_g_Ca_L_Centre_Published = 5.8000000000e-03,  _prvt_g_Ca_L_Periphery_Published = 6.5900000000e-02,  _prvt_g_Ca_L_Centre_0DCapable = 5.7938000000e-03,  _prvt_g_Ca_L_Periphery_0DCapable = 6.5886480000e-02,  _prvt_g_Ca_L_Centre_1DCapable = 8.2000000000e-03,  _prvt_g_Ca_L_Periphery_1DCapable = 6.5900000000e-02,  _prvt_E_Ca_L = 4.6400000000e+01,  _prvt_g_Ca_T_Centre_Published = 4.3000000000e-03,  _prvt_g_Ca_T_Periphery_Published = 1.3900000000e-02,  _prvt_g_Ca_T_Centre_0DCapable = 4.2780600000e-03,  _prvt_g_Ca_T_Periphery_0DCapable = 1.3882300000e-02,  _prvt_g_Ca_T_Centre_1DCapable = 2.1000000000e-03,  _prvt_g_Ca_T_Periphery_1DCapable = 6.9400000000e-03,  _prvt_E_Ca_T = 4.5000000000e+01,  _prvt_g_to_Centre_Published = 4.9100000000e-03,  _prvt_g_to_Periphery_Published = 3.6490000000e-02,  _prvt_g_to_Centre_0DCapable = 4.9050000000e-03,  _prvt_g_to_Periphery_0DCapable = 3.6495000000e-02,  _prvt_g_to_Centre_1DCapable = 4.9050000000e-03,  _prvt_g_to_Periphery_1DCapable = 3.6500000000e-02,  _prvt_g_sus_Centre_Published = 6.6500000000e-05,  _prvt_g_sus_Periphery_Published = 1.1400000000e-02,  _prvt_g_sus_Centre_0DCapable = 6.6455040000e-05,  _prvt_g_sus_Periphery_0DCapable = 1.1383760000e-02,  _prvt_g_sus_Centre_1DCapable = 2.6600000000e-04,  _prvt_g_sus_Periphery_1DCapable = 1.1400000000e-02,  _prvt_g_K_r_Centre_Published = 7.9700000000e-04,  _prvt_g_K_r_Periphery_Published = 1.6000000000e-02,  _prvt_g_K_r_Centre_0DCapable = 7.9704000000e-04,  _prvt_g_K_r_Periphery_0DCapable = 1.6000000000e-02,  _prvt_g_K_r_Centre_1DCapable = 7.3800000000e-04,  _prvt_g_K_r_Periphery_1DCapable = 2.0800000000e-02,  _prvt_g_K_s_Centre_Published = 5.1800000000e-04,  _prvt_g_K_s_Periphery_Published = 1.0400000000e-02,  _prvt_g_K_s_Centre_0DCapable = 3.4450000000e-04,  _prvt_g_K_s_Periphery_0DCapable = 1.0400000000e-02,  _prvt_g_K_s_Centre_1DCapable = 3.4500000000e-04,  _prvt_g_K_s_Periphery_1DCapable = 1.0400000000e-02,  _prvt_g_f_Na_Centre_Published = 5.4800000000e-04,  _prvt_g_f_Na_Periphery_Published = 6.9000000000e-03,  _prvt_g_f_Na_Centre_0DCapable = 5.4650000000e-04,  _prvt_g_f_Na_Periphery_0DCapable = 6.8750000000e-03,  _prvt_g_f_Na_Centre_1DCapable = 4.3700000000e-04,  _prvt_g_f_Na_Periphery_1DCapable = 5.5000000000e-03,  _prvt_g_f_K_Centre_Published = 5.4800000000e-04,  _prvt_g_f_K_Periphery_Published = 6.9000000000e-03,  _prvt_g_f_K_Centre_0DCapable = 5.4650000000e-04,  _prvt_g_f_K_Periphery_0DCapable = 6.8750000000e-03,  _prvt_g_f_K_Centre_1DCapable = 4.3700000000e-04,  _prvt_g_f_K_Periphery_1DCapable = 5.5000000000e-03,  _prvt_g_b_Na_Centre_Published = 5.8000000000e-05,  _prvt_g_b_Na_Periphery_Published = 1.8900000000e-04,  _prvt_g_b_Na_Centre_0DCapable = 5.8181800000e-05,  _prvt_g_b_Na_Periphery_0DCapable = 1.8880000000e-04,  _prvt_g_b_Na_Centre_1DCapable = 5.8000000000e-05,  _prvt_g_b_Na_Periphery_1DCapable = 1.8900000000e-04,  _prvt_g_b_K_Centre_Published = 2.5200000000e-05,  _prvt_g_b_K_Periphery_Published = 8.1900000000e-05,  _prvt_g_b_K_Centre_0DCapable = 2.5236360000e-05,  _prvt_g_b_K_Periphery_0DCapable = 8.1892000000e-05,  _prvt_g_b_K_Centre_1DCapable = 2.5200000000e-05,  _prvt_g_b_K_Periphery_1DCapable = 8.1900000000e-05,  _prvt_g_b_Ca_Centre_Published = 1.3200000000e-05,  _prvt_g_b_Ca_Periphery_Published = 4.3000000000e-05,  _prvt_g_b_Ca_Centre_0DCapable = 1.3236000000e-05,  _prvt_g_b_Ca_Periphery_0DCapable = 4.2952000000e-05,  _prvt_g_b_Ca_Centre_1DCapable = 1.3230000000e-05,  _prvt_g_b_Ca_Periphery_1DCapable = 4.2900000000e-05,  _prvt_k_NaCa_Centre_Published = 2.7000000000e-06,  _prvt_k_NaCa_Periphery_Published = 8.8000000000e-06,  _prvt_k_NaCa_Centre_0DCapable = 2.7229000000e-06,  _prvt_k_NaCa_Periphery_0DCapable = 8.8358400000e-06,  _prvt_k_NaCa_Centre_1DCapable = 2.8000000000e-06,  _prvt_k_NaCa_Periphery_1DCapable = 8.8000000000e-06,  _prvt_Na_i = 8.0000000000e+00,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_gamma_NaCa = 5.0000000000e-01,  _prvt_Ca_i = 1.0000000000e-04,  _prvt_d_NaCa = 1.0000000000e-04,  _prvt_i_p_max_Centre_Published = 4.7800000000e-02,  _prvt_i_p_max_Periphery_Published = 1.6000000000e-01,  _prvt_i_p_max_Centre_0DCapable = 4.7825450000e-02,  _prvt_i_p_max_Periphery_0DCapable = 1.5519360000e-01,  _prvt_i_p_max_Centre_1DCapable = 4.7800000000e-02,  _prvt_i_p_max_Periphery_1DCapable = 1.6000000000e-01,  _prvt_K_m_Na = 5.6400000000e+00,  _prvt_K_o = 5.4000000000e+00,  _prvt_K_m_K = 6.2100000000e-01,  _prvt_i_Ca_p_max_Centre_Published = 0.0000000000e+00,  _prvt_i_Ca_p_max_Periphery_Published = 0.0000000000e+00,  _prvt_i_Ca_p_max_Centre_0DCapable = 0.0000000000e+00,  _prvt_i_Ca_p_max_Periphery_0DCapable = 0.0000000000e+00,  _prvt_i_Ca_p_max_Centre_1DCapable = 4.2000000000e-03,  _prvt_i_Ca_p_max_Periphery_1DCapable = 3.3390000000e-02,  _prvt_K_i = 1.4000000000e+02, 
			//private aux variables
			 _prvt_calc_FCell=0.0,  _prvt_calc_Cm=0.0,  _prvt_calc_g_Na=0.0,  _prvt_calc_i_Na=0.0,  _prvt_calc_m_infinity=0.0,  _prvt_calc_tau_m=0.0,  _prvt_calc_F_Na=0.0,  _prvt_calc_h=0.0,  _prvt_calc_h1_infinity=0.0,  _prvt_calc_h2_infinity=0.0,  _prvt_calc_tau_h1=0.0,  _prvt_calc_tau_h2=0.0,  _prvt_calc_g_Ca_L=0.0,  _prvt_calc_i_Ca_L=0.0,  _prvt_calc_alpha_d_L=0.0,  _prvt_calc_beta_d_L=0.0,  _prvt_calc_tau_d_L=0.0,  _prvt_calc_d_L_infinity=0.0,  _prvt_calc_alpha_f_L=0.0,  _prvt_calc_beta_f_L=0.0,  _prvt_calc_tau_f_L=0.0,  _prvt_calc_f_L_infinity=0.0,  _prvt_calc_g_Ca_T=0.0,  _prvt_calc_i_Ca_T=0.0,  _prvt_calc_alpha_d_T=0.0,  _prvt_calc_beta_d_T=0.0,  _prvt_calc_tau_d_T=0.0,  _prvt_calc_d_T_infinity=0.0,  _prvt_calc_alpha_f_T=0.0,  _prvt_calc_beta_f_T=0.0,  _prvt_calc_tau_f_T=0.0,  _prvt_calc_f_T_infinity=0.0,  _prvt_calc_g_to=0.0,  _prvt_calc_g_sus=0.0,  _prvt_calc_i_to=0.0,  _prvt_calc_i_sus=0.0,  _prvt_calc_q_infinity=0.0,  _prvt_calc_tau_q=0.0,  _prvt_calc_r_infinity=0.0,  _prvt_calc_tau_r=0.0,  _prvt_calc_g_K_r=0.0,  _prvt_calc_i_K_r=0.0,  _prvt_calc_P_a=0.0,  _prvt_calc_P_af_infinity=0.0,  _prvt_calc_tau_P_af=0.0,  _prvt_calc_P_as_infinity=0.0,  _prvt_calc_tau_P_as=0.0,  _prvt_calc_tau_P_i=0.0,  _prvt_calc_P_i_infinity=0.0,  _prvt_calc_g_K_s=0.0,  _prvt_calc_i_K_s=0.0,  _prvt_calc_alpha_xs=0.0,  _prvt_calc_beta_xs=0.0,  _prvt_calc_g_f_Na=0.0,  _prvt_calc_i_f_Na=0.0,  _prvt_calc_g_f_K=0.0,  _prvt_calc_i_f_K=0.0,  _prvt_calc_alpha_y=0.0,  _prvt_calc_beta_y=0.0,  _prvt_calc_g_b_Na=0.0,  _prvt_calc_i_b_Na=0.0,  _prvt_calc_g_b_K=0.0,  _prvt_calc_i_b_K=0.0,  _prvt_calc_g_b_Ca=0.0,  _prvt_calc_i_b_Ca=0.0,  _prvt_calc_k_NaCa=0.0,  _prvt_calc_i_NaCa=0.0,  _prvt_calc_i_p_max=0.0,  _prvt_calc_i_p=0.0,  _prvt_calc_i_Ca_p_max=0.0,  _prvt_calc_i_Ca_p=0.0,  _prvt_calc_E_Na=0.0,  _prvt_calc_E_K=0.0,  _prvt_calc_E_Ca=0.0,  _prvt_calc_E_K_s=0.0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h1_lado_direito_,  _prvt_h2_lado_direito_,  _prvt_d_L_lado_direito_,  _prvt_f_L_lado_direito_,  _prvt_d_T_lado_direito_,  _prvt_f_T_lado_direito_,  _prvt_q_lado_direito_,  _prvt_r_lado_direito_,  _prvt_P_af_lado_direito_,  _prvt_P_as_lado_direito_,  _prvt_P_i_lado_direito_,  _prvt_xs_lado_direito_,  _prvt_y_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;
			__NEW_[0] = __OLD_[0] = -3.9013558536e+01;
			__NEW_[1] = __OLD_[1] = 9.2361701692e-02;
			__NEW_[2] = __OLD_[2] = 1.5905380261e-02;
			__NEW_[3] = __OLD_[3] = 1.4452161090e-02;
			__NEW_[4] = __OLD_[4] = 4.8049008950e-02;
			__NEW_[5] = __OLD_[5] = 4.8779845203e-01;
			__NEW_[6] = __OLD_[6] = 4.2074047435e-01;
			__NEW_[7] = __OLD_[7] = 3.8968420558e-02;
			__NEW_[8] = __OLD_[8] = 2.9760539675e-01;
			__NEW_[9] = __OLD_[9] = 6.4402950262e-02;
			__NEW_[10] = __OLD_[10] = 1.3034201158e-01;
			__NEW_[11] = __OLD_[11] = 4.6960956028e-01;
			__NEW_[12] = __OLD_[12] = 8.7993375273e-01;
			__NEW_[13] = __OLD_[13] = 8.2293827208e-02;
			__NEW_[14] = __OLD_[14] = 3.8892917590e-02;
			_prvt_time_new += _prvt_dtime;
			if(omp_get_thread_num()==_prvt_tree_thread[0])
			{
				_prvt_calc_FCell = ((_prvt_Version==0.0000000000e+00))
?(((1.0700000000e+00*((3.0000000000e+00*_prvt_dCell)-1.0000000000e-01))/(3.0000000000e+00*(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01)))))))
:(((_prvt_Version==1.0000000000e+00))
?(((_prvt_FCellConstant*_prvt_dCell)/(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01))))))
:(((1.0700000000e+00*2.9000000000e+01*_prvt_dCell)/(3.0000000000e+01*(1.0000000000e+00+(7.7450000000e-01*exp(((-((2.9000000000e+01*_prvt_dCell)-2.4500000000e+01))/1.9500000000e+00))))))))
;
				_prvt_calc_F_Na = ((_prvt_Version==0.0000000000e+00))
?((((9.5200000000e-02*exp(((-6.3000000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6600000000e+00*exp(((-2.2500000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6900000000e-02))
:((((9.5180000000e-02*exp(((-6.3060000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6620000000e+00*exp(((-2.2510000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6930000000e-02));
				_prvt_calc_alpha_f_L = ((_prvt_Version==1.0000000000e+00))
?(((3.7500000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)))
:(((3.1200000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
				_prvt_calc_beta_f_L = ((_prvt_Version==1.0000000000e+00))
?((3.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))))
:((2.5000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))));
				_prvt_calc_f_L_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5000000000e+01)/5.0000000000e+00))));
				_prvt_calc_beta_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5000000000e+01*exp(((__OLD_[0]+7.1000000000e+01)/1.5380000000e+01))))
:((1.5000000000e+01*exp(((__OLD_[0]+7.1700000000e+01)/1.5380000000e+01))));
				_prvt_calc_f_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+7.1000000000e+01)/9.0000000000e+00))));
				_prvt_calc_P_a = ((6.0000000000e-01*__OLD_[10])+(4.0000000000e-01*__OLD_[11]));
				_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/_prvt_Na_i)));
				_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/_prvt_K_i)));
				_prvt_calc_E_Ca = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Ca_o/_prvt_Ca_i)));
				_prvt_calc_E_K_s = ((_prvt_Version==0.0000000000e+00))
?((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(1.2000000000e-01*_prvt_Na_o))/(_prvt_K_i+(1.2000000000e-01*_prvt_Na_i))))))
:((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(3.0000000000e-02*_prvt_Na_o))/(_prvt_K_i+(3.0000000000e-02*_prvt_Na_i))))));
				_prvt_calc_Cm = (_prvt_CmCentre+(_prvt_calc_FCell*(_prvt_CmPeriphery-_prvt_CmCentre)));
				_prvt_calc_g_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_Published-_prvt_g_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_0DCapable-_prvt_g_Na_Centre_0DCapable))))
:((_prvt_g_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_1DCapable-_prvt_g_Na_Centre_1DCapable)))))
;
				_prvt_calc_h = (((1.0000000000e+00-_prvt_calc_F_Na)*__OLD_[2])+(_prvt_calc_F_Na*__OLD_[3]));
				_prvt_calc_g_Ca_L = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_L_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_Published-_prvt_g_Ca_L_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_L_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_0DCapable-_prvt_g_Ca_L_Centre_0DCapable))))
:((_prvt_g_Ca_L_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_1DCapable-_prvt_g_Ca_L_Centre_1DCapable)))))
;
				_prvt_calc_tau_f_L = ((_prvt_Version==1.0000000000e+00))
?(((1.2000000000e+00-(2.0000000000e-01*_prvt_calc_FCell))/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)))
:((1.0000000000e+00/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)));
				_prvt_calc_g_Ca_T = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_T_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_Published-_prvt_g_Ca_T_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_T_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_0DCapable-_prvt_g_Ca_T_Centre_0DCapable))))
:((_prvt_g_Ca_T_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_1DCapable-_prvt_g_Ca_T_Centre_1DCapable)))))
;
				_prvt_calc_alpha_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5300000000e+01*exp(((-(__OLD_[0]+7.1000000000e+01+(7.0000000000e-01*_prvt_calc_FCell)))/8.3300000000e+01))))
:((1.5300000000e+01*exp(((-(__OLD_[0]+7.1700000000e+01))/8.3300000000e+01))));
				_prvt_calc_g_to = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_to_Centre_Published+(_prvt_calc_FCell*(_prvt_g_to_Periphery_Published-_prvt_g_to_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_to_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_0DCapable-_prvt_g_to_Centre_0DCapable))))
:((_prvt_g_to_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_1DCapable-_prvt_g_to_Centre_1DCapable)))))
;
				_prvt_calc_g_sus = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_sus_Centre_Published+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_Published-_prvt_g_sus_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_sus_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_0DCapable-_prvt_g_sus_Centre_0DCapable))))
:((_prvt_g_sus_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_1DCapable-_prvt_g_sus_Centre_1DCapable)))))
;
				_prvt_calc_g_K_r = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_r_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_Published-_prvt_g_K_r_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_r_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_0DCapable-_prvt_g_K_r_Centre_0DCapable))))
:((_prvt_g_K_r_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_1DCapable-_prvt_g_K_r_Centre_1DCapable)))))
;
				_prvt_calc_g_K_s = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_s_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_Published-_prvt_g_K_s_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_s_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_0DCapable-_prvt_g_K_s_Centre_0DCapable))))
:((_prvt_g_K_s_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_1DCapable-_prvt_g_K_s_Centre_1DCapable)))))
;
				_prvt_calc_g_f_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_Published-_prvt_g_f_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_0DCapable-_prvt_g_f_Na_Centre_0DCapable))))
:((_prvt_g_f_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_1DCapable-_prvt_g_f_Na_Centre_1DCapable)))))
;
				_prvt_calc_g_f_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_Published-_prvt_g_f_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_0DCapable-_prvt_g_f_K_Centre_0DCapable))))
:((_prvt_g_f_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_1DCapable-_prvt_g_f_K_Centre_1DCapable)))))
;
				_prvt_calc_g_b_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_Published-_prvt_g_b_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_0DCapable-_prvt_g_b_Na_Centre_0DCapable))))
:((_prvt_g_b_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_1DCapable-_prvt_g_b_Na_Centre_1DCapable)))))
;
				_prvt_calc_g_b_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_Published-_prvt_g_b_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_0DCapable-_prvt_g_b_K_Centre_0DCapable))))
:((_prvt_g_b_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_1DCapable-_prvt_g_b_K_Centre_1DCapable)))))
;
				_prvt_calc_g_b_Ca = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Ca_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_Published-_prvt_g_b_Ca_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Ca_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_0DCapable-_prvt_g_b_Ca_Centre_0DCapable))))
:((_prvt_g_b_Ca_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_1DCapable-_prvt_g_b_Ca_Centre_1DCapable)))))
;
				_prvt_calc_k_NaCa = ((_prvt_Version==0.0000000000e+00))
?((_prvt_k_NaCa_Centre_Published+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_Published-_prvt_k_NaCa_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_k_NaCa_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_0DCapable-_prvt_k_NaCa_Centre_0DCapable))))
:((_prvt_k_NaCa_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_1DCapable-_prvt_k_NaCa_Centre_1DCapable)))))
;
				_prvt_calc_i_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_Published-_prvt_i_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_0DCapable-_prvt_i_p_max_Centre_0DCapable))))
:((_prvt_i_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_1DCapable-_prvt_i_p_max_Centre_1DCapable)))))
;
				_prvt_calc_i_Ca_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_Published-_prvt_i_Ca_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_0DCapable-_prvt_i_Ca_p_max_Centre_0DCapable))))
:((_prvt_i_Ca_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_1DCapable-_prvt_i_Ca_p_max_Centre_1DCapable)))))
;
				_prvt_calc_i_Ca_L = (_prvt_calc_g_Ca_L*((__OLD_[5]*__OLD_[4])+(6.0000000000e-03/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4100000000e+01))/6.0000000000e+00)))))*(__OLD_[0]-_prvt_E_Ca_L));
				_prvt_calc_i_Ca_T = (_prvt_calc_g_Ca_T*__OLD_[6]*__OLD_[7]*(__OLD_[0]-_prvt_E_Ca_T));
				_prvt_calc_tau_f_T = (1.0000000000e+00/(_prvt_calc_alpha_f_T+_prvt_calc_beta_f_T));
				_prvt_calc_i_to = (_prvt_calc_g_to*__OLD_[8]*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_sus = (_prvt_calc_g_sus*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_K_r = (_prvt_calc_g_K_r*_prvt_calc_P_a*__OLD_[12]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_K_s = (_prvt_calc_g_K_s*pow(__OLD_[13],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K_s));
				_prvt_calc_i_f_Na = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_Na)))
:((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-7.7600000000e+01)));
				_prvt_calc_i_f_K = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_K)))
:((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]+1.0200000000e+02)));
				_prvt_calc_i_b_Na = (_prvt_calc_g_b_Na*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_b_K = (_prvt_calc_g_b_K*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_b_Ca = (_prvt_calc_g_b_Ca*(__OLD_[0]-_prvt_calc_E_Ca));
				_prvt_calc_i_NaCa = ((_prvt_Version==0.0000000000e+00))
?(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7400000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))))
:(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7430000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))));
				_prvt_calc_i_p = ((_prvt_calc_i_p_max*pow((_prvt_Na_i/(_prvt_K_m_Na+_prvt_Na_i)),3.0000000000e+00)*pow((_prvt_K_o/(_prvt_K_m_K+_prvt_K_o)),2.0000000000e+00)*1.6000000000e+00)/(1.5000000000e+00+exp(((-(__OLD_[0]+6.0000000000e+01))/4.0000000000e+01))));
				_prvt_calc_i_Ca_p = ((_prvt_calc_i_Ca_p_max*_prvt_Ca_i)/(_prvt_Ca_i+4.0000000000e-04));
				_prvt_calc_i_Na = (((((_prvt_calc_g_Na*pow(__OLD_[1],3.0000000000e+00)*_prvt_calc_h*_prvt_Na_o*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*(exp((((__OLD_[0]-_prvt_calc_E_Na)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))/(exp(((__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))*__OLD_[0]);
				__K1_[0]= (((-1.0000000000e+00)/_prvt_calc_Cm)*(_prvt_calc_i_Na+_prvt_calc_i_Ca_L+_prvt_calc_i_Ca_T+_prvt_calc_i_to+_prvt_calc_i_sus+_prvt_calc_i_K_r+_prvt_calc_i_K_s+_prvt_calc_i_f_Na+_prvt_calc_i_f_K+_prvt_calc_i_b_Na+_prvt_calc_i_b_Ca+_prvt_calc_i_b_K+_prvt_calc_i_NaCa+_prvt_calc_i_p+_prvt_calc_i_Ca_p));
				__NEW_[0]= __K1_[0] * _prvt_dtime + __OLD_[0];
				__K1_[5]= ((_prvt_calc_f_L_infinity-__OLD_[5])/_prvt_calc_tau_f_L);
				__NEW_[5]= __K1_[5] * _prvt_dtime + __OLD_[5];
				__K1_[7]= ((_prvt_calc_f_T_infinity-__OLD_[7])/_prvt_calc_tau_f_T);
				__NEW_[7]= __K1_[7] * _prvt_dtime + __OLD_[7];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[1])
			{
				_prvt_calc_m_infinity = ((_prvt_Version==0.0000000000e+00))
?(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-__OLD_[0])/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)))
:(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.0320000000e+01))/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
				_prvt_calc_tau_m = ((_prvt_Version==0.0000000000e+00))
?(((6.2470000000e-04/((8.3200000000e-01*exp(((-3.3500000000e-01)*(__OLD_[0]+5.6700000000e+01))))+(6.2700000000e-01*exp((8.2000000000e-02*(__OLD_[0]+6.5010000000e+01))))))+4.0000000000e-05))
:(((6.2470000000e-04/((8.3221660000e-01*exp(((-3.3566000000e-01)*(__OLD_[0]+5.6706200000e+01))))+(6.2740000000e-01*exp((8.2300000000e-02*(__OLD_[0]+6.5013100000e+01))))))+4.5690000000e-05));
				__K1_[1]= ((_prvt_calc_m_infinity-__OLD_[1])/_prvt_calc_tau_m);
				__NEW_[1]= __K1_[1] * _prvt_dtime + __OLD_[1];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[2])
			{
				_prvt_calc_h1_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+6.6100000000e+01)/6.4000000000e+00))));
				_prvt_calc_tau_h1 = (((3.7170000000e-06*exp(((-2.8150000000e-01)*(__OLD_[0]+1.7110000000e+01))))/(1.0000000000e+00+(3.7320000000e-03*exp(((-3.4260000000e-01)*(__OLD_[0]+3.7760000000e+01))))))+5.9770000000e-04);
				_prvt_calc_tau_h2 = (((3.1860000000e-08*exp(((-6.2190000000e-01)*(__OLD_[0]+1.8800000000e+01))))/(1.0000000000e+00+(7.1890000000e-05*exp(((-6.6830000000e-01)*(__OLD_[0]+3.4070000000e+01))))))+3.5560000000e-03);
				_prvt_calc_h2_infinity = _prvt_calc_h1_infinity;
				__K1_[2]= ((_prvt_calc_h1_infinity-__OLD_[2])/_prvt_calc_tau_h1);
				__NEW_[2]= __K1_[2] * _prvt_dtime + __OLD_[2];
				__K1_[3]= ((_prvt_calc_h2_infinity-__OLD_[3])/_prvt_calc_tau_h2);
				__NEW_[3]= __K1_[3] * _prvt_dtime + __OLD_[3];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[3])
			{
				_prvt_calc_alpha_d_L = ((_prvt_Version==0.0000000000e+00))
?(((((-2.8380000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((_prvt_Version==1.0000000000e+00))
?(((((-2.8390000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((((-2.8400000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00)))))
;
				_prvt_calc_beta_d_L = ((_prvt_Version==1.0000000000e+00))
?(((1.1430000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)))
:(((1.1420000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)));
				_prvt_calc_d_L_infinity = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.3100000000e+01))/6.0000000000e+00)))))
:(((_prvt_Version==1.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2300000000e+01+(8.0000000000e-01*_prvt_calc_FCell)))/6.0000000000e+00)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2200000000e+01))/6.0000000000e+00))))))
;
				_prvt_calc_tau_d_L = (2.0000000000e+00/(_prvt_calc_alpha_d_L+_prvt_calc_beta_d_L));
				__K1_[4]= ((_prvt_calc_d_L_infinity-__OLD_[4])/_prvt_calc_tau_d_L);
				__NEW_[4]= __K1_[4] * _prvt_dtime + __OLD_[4];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[4])
			{
				_prvt_calc_alpha_d_T = (1.0680000000e+03*exp(((__OLD_[0]+2.6300000000e+01)/3.0000000000e+01)));
				_prvt_calc_beta_d_T = (1.0680000000e+03*exp(((-(__OLD_[0]+2.6300000000e+01))/3.0000000000e+01)));
				_prvt_calc_d_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.7000000000e+01))/6.8000000000e+00))));
				_prvt_calc_tau_d_T = (1.0000000000e+00/(_prvt_calc_alpha_d_T+_prvt_calc_beta_d_T));
				__K1_[6]= ((_prvt_calc_d_T_infinity-__OLD_[6])/_prvt_calc_tau_d_T);
				__NEW_[6]= __K1_[6] * _prvt_dtime + __OLD_[6];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[5])
			{
				_prvt_calc_q_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.9370000000e+01)/1.3100000000e+01))));
				_prvt_calc_tau_q = ((_prvt_Version==0.0000000000e+00))
?((1.0100000000e-02+(6.5170000000e-02/(5.7000000000e-01*exp(((-8.0000000000e-02)*(__OLD_[0]+4.9000000000e+01)))))+(2.4000000000e-05*exp((1.0000000000e-01*(__OLD_[0]+5.0930000000e+01))))))
:(((_prvt_Version==1.0000000000e+00))
?(((1.0000000000e-03/3.0000000000e+00)*(3.0310000000e+01+(1.9550000000e+02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01+(1.0000000000e+01*_prvt_calc_FCell)))))+(7.1740000000e-01*exp(((2.7190000000e-01-(1.7190000000e-01*_prvt_calc_FCell))*1.0000000000e+00*(__OLD_[0]+4.0930000000e+01+(1.0000000000e+01*_prvt_calc_FCell))))))))))
:((1.0100000000e-02+(6.5170000000e-02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01))))+(7.1740000000e-01*exp((2.7190000000e-01*(__OLD_[0]+4.0930000000e+01)))))))))
;
				__K1_[8]= ((_prvt_calc_q_infinity-__OLD_[8])/_prvt_calc_tau_q);
				__NEW_[8]= __K1_[8] * _prvt_dtime + __OLD_[8];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[6])
			{
				_prvt_calc_r_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]-1.0930000000e+01))/1.9700000000e+01))));
				_prvt_calc_tau_r = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e-03*(2.9800000000e+00+(1.5590000000e+01/((1.0370000000e+00*exp((9.0000000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.2000000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:(((_prvt_Version==1.0000000000e+00))
?((2.5000000000e-03*(1.1910000000e+00+(7.8380000000e+00/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:((1.0000000000e-03*(2.9800000000e+00+(1.9590000000e+01/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01))))))))))
;
				__K1_[9]= ((_prvt_calc_r_infinity-__OLD_[9])/_prvt_calc_tau_r);
				__NEW_[9]= __K1_[9] * _prvt_dtime + __OLD_[9];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[7])
			{
				_prvt_calc_P_af_infinity = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4200000000e+01))/1.0600000000e+01)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.3200000000e+01))/1.0600000000e+01)))));
				_prvt_calc_tau_P_af = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-9.0000000000e+00)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.2500000000e+01))))))
:((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-1.0000000000e+01)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.2500000000e+01))))));
				_prvt_calc_tau_P_as = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-9.0000000000e+00)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.1600000000e+01))))))
:((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-1.0000000000e+01)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.1600000000e+01))))));
				_prvt_calc_P_as_infinity = _prvt_calc_P_af_infinity;
				__K1_[10]= ((_prvt_calc_P_af_infinity-__OLD_[10])/_prvt_calc_tau_P_af);
				__NEW_[10]= __K1_[10] * _prvt_dtime + __OLD_[10];
				__K1_[11]= ((_prvt_calc_P_as_infinity-__OLD_[11])/_prvt_calc_tau_P_as);
				__NEW_[11]= __K1_[11] * _prvt_dtime + __OLD_[11];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[8])
			{
				_prvt_calc_tau_P_i = ((_prvt_Version==0.0000000000e+00))
?(2.0000000000e-03)
:(((_prvt_Version==1.0000000000e+00))
?(2.0000000000e-03)
:(6.0000000000e-03))
;
				_prvt_calc_P_i_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+1.8600000000e+01)/1.0100000000e+01))));
				__K1_[12]= ((_prvt_calc_P_i_infinity-__OLD_[12])/_prvt_calc_tau_P_i);
				__NEW_[12]= __K1_[12] * _prvt_dtime + __OLD_[12];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[9])
			{
				_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
				_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
				__K1_[13]= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[13]))-(_prvt_calc_beta_xs*__OLD_[13]));
				__NEW_[13]= __K1_[13] * _prvt_dtime + __OLD_[13];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[10])
			{
				_prvt_calc_alpha_y = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6620000000e+01))))
:((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6630000000e+01))));
				_prvt_calc_beta_y = (1.0000000000e+00*exp(((__OLD_[0]+7.5130000000e+01)/2.1250000000e+01)));
				__K1_[14]= ((_prvt_calc_alpha_y*(1.0000000000e+00-__OLD_[14]))-(_prvt_calc_beta_y*__OLD_[14]));
				__NEW_[14]= __K1_[14] * _prvt_dtime + __OLD_[14];
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
					_prvt_calc_FCell = ((_prvt_Version==0.0000000000e+00))
?(((1.0700000000e+00*((3.0000000000e+00*_prvt_dCell)-1.0000000000e-01))/(3.0000000000e+00*(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01)))))))
:(((_prvt_Version==1.0000000000e+00))
?(((_prvt_FCellConstant*_prvt_dCell)/(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01))))))
:(((1.0700000000e+00*2.9000000000e+01*_prvt_dCell)/(3.0000000000e+01*(1.0000000000e+00+(7.7450000000e-01*exp(((-((2.9000000000e+01*_prvt_dCell)-2.4500000000e+01))/1.9500000000e+00))))))))
;
					_prvt_calc_F_Na = ((_prvt_Version==0.0000000000e+00))
?((((9.5200000000e-02*exp(((-6.3000000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6600000000e+00*exp(((-2.2500000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6900000000e-02))
:((((9.5180000000e-02*exp(((-6.3060000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6620000000e+00*exp(((-2.2510000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6930000000e-02));
					_prvt_calc_alpha_f_L = ((_prvt_Version==1.0000000000e+00))
?(((3.7500000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)))
:(((3.1200000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
					_prvt_calc_beta_f_L = ((_prvt_Version==1.0000000000e+00))
?((3.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))))
:((2.5000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))));
					_prvt_calc_f_L_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5000000000e+01)/5.0000000000e+00))));
					_prvt_calc_beta_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5000000000e+01*exp(((__OLD_[0]+7.1000000000e+01)/1.5380000000e+01))))
:((1.5000000000e+01*exp(((__OLD_[0]+7.1700000000e+01)/1.5380000000e+01))));
					_prvt_calc_f_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+7.1000000000e+01)/9.0000000000e+00))));
					_prvt_calc_P_a = ((6.0000000000e-01*__OLD_[10])+(4.0000000000e-01*__OLD_[11]));
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/_prvt_Na_i)));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/_prvt_K_i)));
					_prvt_calc_E_Ca = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Ca_o/_prvt_Ca_i)));
					_prvt_calc_E_K_s = ((_prvt_Version==0.0000000000e+00))
?((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(1.2000000000e-01*_prvt_Na_o))/(_prvt_K_i+(1.2000000000e-01*_prvt_Na_i))))))
:((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(3.0000000000e-02*_prvt_Na_o))/(_prvt_K_i+(3.0000000000e-02*_prvt_Na_i))))));
					_prvt_calc_Cm = (_prvt_CmCentre+(_prvt_calc_FCell*(_prvt_CmPeriphery-_prvt_CmCentre)));
					_prvt_calc_g_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_Published-_prvt_g_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_0DCapable-_prvt_g_Na_Centre_0DCapable))))
:((_prvt_g_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_1DCapable-_prvt_g_Na_Centre_1DCapable)))))
;
					_prvt_calc_h = (((1.0000000000e+00-_prvt_calc_F_Na)*__OLD_[2])+(_prvt_calc_F_Na*__OLD_[3]));
					_prvt_calc_g_Ca_L = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_L_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_Published-_prvt_g_Ca_L_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_L_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_0DCapable-_prvt_g_Ca_L_Centre_0DCapable))))
:((_prvt_g_Ca_L_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_1DCapable-_prvt_g_Ca_L_Centre_1DCapable)))))
;
					_prvt_calc_tau_f_L = ((_prvt_Version==1.0000000000e+00))
?(((1.2000000000e+00-(2.0000000000e-01*_prvt_calc_FCell))/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)))
:((1.0000000000e+00/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)));
					_prvt_calc_g_Ca_T = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_T_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_Published-_prvt_g_Ca_T_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_T_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_0DCapable-_prvt_g_Ca_T_Centre_0DCapable))))
:((_prvt_g_Ca_T_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_1DCapable-_prvt_g_Ca_T_Centre_1DCapable)))))
;
					_prvt_calc_alpha_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5300000000e+01*exp(((-(__OLD_[0]+7.1000000000e+01+(7.0000000000e-01*_prvt_calc_FCell)))/8.3300000000e+01))))
:((1.5300000000e+01*exp(((-(__OLD_[0]+7.1700000000e+01))/8.3300000000e+01))));
					_prvt_calc_g_to = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_to_Centre_Published+(_prvt_calc_FCell*(_prvt_g_to_Periphery_Published-_prvt_g_to_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_to_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_0DCapable-_prvt_g_to_Centre_0DCapable))))
:((_prvt_g_to_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_1DCapable-_prvt_g_to_Centre_1DCapable)))))
;
					_prvt_calc_g_sus = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_sus_Centre_Published+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_Published-_prvt_g_sus_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_sus_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_0DCapable-_prvt_g_sus_Centre_0DCapable))))
:((_prvt_g_sus_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_1DCapable-_prvt_g_sus_Centre_1DCapable)))))
;
					_prvt_calc_g_K_r = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_r_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_Published-_prvt_g_K_r_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_r_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_0DCapable-_prvt_g_K_r_Centre_0DCapable))))
:((_prvt_g_K_r_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_1DCapable-_prvt_g_K_r_Centre_1DCapable)))))
;
					_prvt_calc_g_K_s = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_s_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_Published-_prvt_g_K_s_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_s_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_0DCapable-_prvt_g_K_s_Centre_0DCapable))))
:((_prvt_g_K_s_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_1DCapable-_prvt_g_K_s_Centre_1DCapable)))))
;
					_prvt_calc_g_f_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_Published-_prvt_g_f_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_0DCapable-_prvt_g_f_Na_Centre_0DCapable))))
:((_prvt_g_f_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_1DCapable-_prvt_g_f_Na_Centre_1DCapable)))))
;
					_prvt_calc_g_f_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_Published-_prvt_g_f_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_0DCapable-_prvt_g_f_K_Centre_0DCapable))))
:((_prvt_g_f_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_1DCapable-_prvt_g_f_K_Centre_1DCapable)))))
;
					_prvt_calc_g_b_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_Published-_prvt_g_b_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_0DCapable-_prvt_g_b_Na_Centre_0DCapable))))
:((_prvt_g_b_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_1DCapable-_prvt_g_b_Na_Centre_1DCapable)))))
;
					_prvt_calc_g_b_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_Published-_prvt_g_b_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_0DCapable-_prvt_g_b_K_Centre_0DCapable))))
:((_prvt_g_b_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_1DCapable-_prvt_g_b_K_Centre_1DCapable)))))
;
					_prvt_calc_g_b_Ca = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Ca_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_Published-_prvt_g_b_Ca_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Ca_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_0DCapable-_prvt_g_b_Ca_Centre_0DCapable))))
:((_prvt_g_b_Ca_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_1DCapable-_prvt_g_b_Ca_Centre_1DCapable)))))
;
					_prvt_calc_k_NaCa = ((_prvt_Version==0.0000000000e+00))
?((_prvt_k_NaCa_Centre_Published+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_Published-_prvt_k_NaCa_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_k_NaCa_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_0DCapable-_prvt_k_NaCa_Centre_0DCapable))))
:((_prvt_k_NaCa_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_1DCapable-_prvt_k_NaCa_Centre_1DCapable)))))
;
					_prvt_calc_i_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_Published-_prvt_i_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_0DCapable-_prvt_i_p_max_Centre_0DCapable))))
:((_prvt_i_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_1DCapable-_prvt_i_p_max_Centre_1DCapable)))))
;
					_prvt_calc_i_Ca_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_Published-_prvt_i_Ca_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_0DCapable-_prvt_i_Ca_p_max_Centre_0DCapable))))
:((_prvt_i_Ca_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_1DCapable-_prvt_i_Ca_p_max_Centre_1DCapable)))))
;
					_prvt_calc_i_Ca_L = (_prvt_calc_g_Ca_L*((__OLD_[5]*__OLD_[4])+(6.0000000000e-03/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4100000000e+01))/6.0000000000e+00)))))*(__OLD_[0]-_prvt_E_Ca_L));
					_prvt_calc_i_Ca_T = (_prvt_calc_g_Ca_T*__OLD_[6]*__OLD_[7]*(__OLD_[0]-_prvt_E_Ca_T));
					_prvt_calc_tau_f_T = (1.0000000000e+00/(_prvt_calc_alpha_f_T+_prvt_calc_beta_f_T));
					_prvt_calc_i_to = (_prvt_calc_g_to*__OLD_[8]*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_sus = (_prvt_calc_g_sus*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K_r = (_prvt_calc_g_K_r*_prvt_calc_P_a*__OLD_[12]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K_s = (_prvt_calc_g_K_s*pow(__OLD_[13],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K_s));
					_prvt_calc_i_f_Na = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_Na)))
:((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-7.7600000000e+01)));
					_prvt_calc_i_f_K = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_K)))
:((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]+1.0200000000e+02)));
					_prvt_calc_i_b_Na = (_prvt_calc_g_b_Na*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_K = (_prvt_calc_g_b_K*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_b_Ca = (_prvt_calc_g_b_Ca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_NaCa = ((_prvt_Version==0.0000000000e+00))
?(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7400000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))))
:(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7430000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))));
					_prvt_calc_i_p = ((_prvt_calc_i_p_max*pow((_prvt_Na_i/(_prvt_K_m_Na+_prvt_Na_i)),3.0000000000e+00)*pow((_prvt_K_o/(_prvt_K_m_K+_prvt_K_o)),2.0000000000e+00)*1.6000000000e+00)/(1.5000000000e+00+exp(((-(__OLD_[0]+6.0000000000e+01))/4.0000000000e+01))));
					_prvt_calc_i_Ca_p = ((_prvt_calc_i_Ca_p_max*_prvt_Ca_i)/(_prvt_Ca_i+4.0000000000e-04));
					_prvt_calc_i_Na = (((((_prvt_calc_g_Na*pow(__OLD_[1],3.0000000000e+00)*_prvt_calc_h*_prvt_Na_o*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*(exp((((__OLD_[0]-_prvt_calc_E_Na)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))/(exp(((__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))*__OLD_[0]);
					__K2_[0]= (((-1.0000000000e+00)/_prvt_calc_Cm)*(_prvt_calc_i_Na+_prvt_calc_i_Ca_L+_prvt_calc_i_Ca_T+_prvt_calc_i_to+_prvt_calc_i_sus+_prvt_calc_i_K_r+_prvt_calc_i_K_s+_prvt_calc_i_f_Na+_prvt_calc_i_f_K+_prvt_calc_i_b_Na+_prvt_calc_i_b_Ca+_prvt_calc_i_b_K+_prvt_calc_i_NaCa+_prvt_calc_i_p+_prvt_calc_i_Ca_p));
					_prvt_aux_tol = fabs(__OLD_[0])*_prvt_rel_tol_;
					__TOL_[0] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[0] = fabs((_prvt_dtime/2) * (__K1_[0] - __K2_[0])/__TOL_[0]);
					__K2_[5]= ((_prvt_calc_f_L_infinity-__OLD_[5])/_prvt_calc_tau_f_L);
					_prvt_aux_tol = fabs(__OLD_[5])*_prvt_rel_tol_;
					__TOL_[5] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[5] = fabs((_prvt_dtime/2) * (__K1_[5] - __K2_[5])/__TOL_[5]);
					__K2_[7]= ((_prvt_calc_f_T_infinity-__OLD_[7])/_prvt_calc_tau_f_T);
					_prvt_aux_tol = fabs(__OLD_[7])*_prvt_rel_tol_;
					__TOL_[7] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[7] = fabs((_prvt_dtime/2) * (__K1_[7] - __K2_[7])/__TOL_[7]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[1])
				{
					_prvt_calc_m_infinity = ((_prvt_Version==0.0000000000e+00))
?(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-__OLD_[0])/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)))
:(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.0320000000e+01))/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
					_prvt_calc_tau_m = ((_prvt_Version==0.0000000000e+00))
?(((6.2470000000e-04/((8.3200000000e-01*exp(((-3.3500000000e-01)*(__OLD_[0]+5.6700000000e+01))))+(6.2700000000e-01*exp((8.2000000000e-02*(__OLD_[0]+6.5010000000e+01))))))+4.0000000000e-05))
:(((6.2470000000e-04/((8.3221660000e-01*exp(((-3.3566000000e-01)*(__OLD_[0]+5.6706200000e+01))))+(6.2740000000e-01*exp((8.2300000000e-02*(__OLD_[0]+6.5013100000e+01))))))+4.5690000000e-05));
					__K2_[1]= ((_prvt_calc_m_infinity-__OLD_[1])/_prvt_calc_tau_m);
					_prvt_aux_tol = fabs(__OLD_[1])*_prvt_rel_tol_;
					__TOL_[1] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[1] = fabs((_prvt_dtime/2) * (__K1_[1] - __K2_[1])/__TOL_[1]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[2])
				{
					_prvt_calc_h1_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+6.6100000000e+01)/6.4000000000e+00))));
					_prvt_calc_tau_h1 = (((3.7170000000e-06*exp(((-2.8150000000e-01)*(__OLD_[0]+1.7110000000e+01))))/(1.0000000000e+00+(3.7320000000e-03*exp(((-3.4260000000e-01)*(__OLD_[0]+3.7760000000e+01))))))+5.9770000000e-04);
					_prvt_calc_tau_h2 = (((3.1860000000e-08*exp(((-6.2190000000e-01)*(__OLD_[0]+1.8800000000e+01))))/(1.0000000000e+00+(7.1890000000e-05*exp(((-6.6830000000e-01)*(__OLD_[0]+3.4070000000e+01))))))+3.5560000000e-03);
					_prvt_calc_h2_infinity = _prvt_calc_h1_infinity;
					__K2_[2]= ((_prvt_calc_h1_infinity-__OLD_[2])/_prvt_calc_tau_h1);
					_prvt_aux_tol = fabs(__OLD_[2])*_prvt_rel_tol_;
					__TOL_[2] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[2] = fabs((_prvt_dtime/2) * (__K1_[2] - __K2_[2])/__TOL_[2]);
					__K2_[3]= ((_prvt_calc_h2_infinity-__OLD_[3])/_prvt_calc_tau_h2);
					_prvt_aux_tol = fabs(__OLD_[3])*_prvt_rel_tol_;
					__TOL_[3] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[3] = fabs((_prvt_dtime/2) * (__K1_[3] - __K2_[3])/__TOL_[3]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[3])
				{
					_prvt_calc_alpha_d_L = ((_prvt_Version==0.0000000000e+00))
?(((((-2.8380000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((_prvt_Version==1.0000000000e+00))
?(((((-2.8390000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((((-2.8400000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00)))))
;
					_prvt_calc_beta_d_L = ((_prvt_Version==1.0000000000e+00))
?(((1.1430000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)))
:(((1.1420000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)));
					_prvt_calc_d_L_infinity = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.3100000000e+01))/6.0000000000e+00)))))
:(((_prvt_Version==1.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2300000000e+01+(8.0000000000e-01*_prvt_calc_FCell)))/6.0000000000e+00)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2200000000e+01))/6.0000000000e+00))))))
;
					_prvt_calc_tau_d_L = (2.0000000000e+00/(_prvt_calc_alpha_d_L+_prvt_calc_beta_d_L));
					__K2_[4]= ((_prvt_calc_d_L_infinity-__OLD_[4])/_prvt_calc_tau_d_L);
					_prvt_aux_tol = fabs(__OLD_[4])*_prvt_rel_tol_;
					__TOL_[4] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[4] = fabs((_prvt_dtime/2) * (__K1_[4] - __K2_[4])/__TOL_[4]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[4])
				{
					_prvt_calc_alpha_d_T = (1.0680000000e+03*exp(((__OLD_[0]+2.6300000000e+01)/3.0000000000e+01)));
					_prvt_calc_beta_d_T = (1.0680000000e+03*exp(((-(__OLD_[0]+2.6300000000e+01))/3.0000000000e+01)));
					_prvt_calc_d_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.7000000000e+01))/6.8000000000e+00))));
					_prvt_calc_tau_d_T = (1.0000000000e+00/(_prvt_calc_alpha_d_T+_prvt_calc_beta_d_T));
					__K2_[6]= ((_prvt_calc_d_T_infinity-__OLD_[6])/_prvt_calc_tau_d_T);
					_prvt_aux_tol = fabs(__OLD_[6])*_prvt_rel_tol_;
					__TOL_[6] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[6] = fabs((_prvt_dtime/2) * (__K1_[6] - __K2_[6])/__TOL_[6]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[5])
				{
					_prvt_calc_q_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.9370000000e+01)/1.3100000000e+01))));
					_prvt_calc_tau_q = ((_prvt_Version==0.0000000000e+00))
?((1.0100000000e-02+(6.5170000000e-02/(5.7000000000e-01*exp(((-8.0000000000e-02)*(__OLD_[0]+4.9000000000e+01)))))+(2.4000000000e-05*exp((1.0000000000e-01*(__OLD_[0]+5.0930000000e+01))))))
:(((_prvt_Version==1.0000000000e+00))
?(((1.0000000000e-03/3.0000000000e+00)*(3.0310000000e+01+(1.9550000000e+02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01+(1.0000000000e+01*_prvt_calc_FCell)))))+(7.1740000000e-01*exp(((2.7190000000e-01-(1.7190000000e-01*_prvt_calc_FCell))*1.0000000000e+00*(__OLD_[0]+4.0930000000e+01+(1.0000000000e+01*_prvt_calc_FCell))))))))))
:((1.0100000000e-02+(6.5170000000e-02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01))))+(7.1740000000e-01*exp((2.7190000000e-01*(__OLD_[0]+4.0930000000e+01)))))))))
;
					__K2_[8]= ((_prvt_calc_q_infinity-__OLD_[8])/_prvt_calc_tau_q);
					_prvt_aux_tol = fabs(__OLD_[8])*_prvt_rel_tol_;
					__TOL_[8] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[8] = fabs((_prvt_dtime/2) * (__K1_[8] - __K2_[8])/__TOL_[8]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[6])
				{
					_prvt_calc_r_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]-1.0930000000e+01))/1.9700000000e+01))));
					_prvt_calc_tau_r = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e-03*(2.9800000000e+00+(1.5590000000e+01/((1.0370000000e+00*exp((9.0000000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.2000000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:(((_prvt_Version==1.0000000000e+00))
?((2.5000000000e-03*(1.1910000000e+00+(7.8380000000e+00/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:((1.0000000000e-03*(2.9800000000e+00+(1.9590000000e+01/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01))))))))))
;
					__K2_[9]= ((_prvt_calc_r_infinity-__OLD_[9])/_prvt_calc_tau_r);
					_prvt_aux_tol = fabs(__OLD_[9])*_prvt_rel_tol_;
					__TOL_[9] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[9] = fabs((_prvt_dtime/2) * (__K1_[9] - __K2_[9])/__TOL_[9]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[7])
				{
					_prvt_calc_P_af_infinity = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4200000000e+01))/1.0600000000e+01)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.3200000000e+01))/1.0600000000e+01)))));
					_prvt_calc_tau_P_af = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-9.0000000000e+00)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.2500000000e+01))))))
:((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-1.0000000000e+01)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.2500000000e+01))))));
					_prvt_calc_tau_P_as = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-9.0000000000e+00)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.1600000000e+01))))))
:((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-1.0000000000e+01)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.1600000000e+01))))));
					_prvt_calc_P_as_infinity = _prvt_calc_P_af_infinity;
					__K2_[10]= ((_prvt_calc_P_af_infinity-__OLD_[10])/_prvt_calc_tau_P_af);
					_prvt_aux_tol = fabs(__OLD_[10])*_prvt_rel_tol_;
					__TOL_[10] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[10] = fabs((_prvt_dtime/2) * (__K1_[10] - __K2_[10])/__TOL_[10]);
					__K2_[11]= ((_prvt_calc_P_as_infinity-__OLD_[11])/_prvt_calc_tau_P_as);
					_prvt_aux_tol = fabs(__OLD_[11])*_prvt_rel_tol_;
					__TOL_[11] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[11] = fabs((_prvt_dtime/2) * (__K1_[11] - __K2_[11])/__TOL_[11]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[8])
				{
					_prvt_calc_tau_P_i = ((_prvt_Version==0.0000000000e+00))
?(2.0000000000e-03)
:(((_prvt_Version==1.0000000000e+00))
?(2.0000000000e-03)
:(6.0000000000e-03))
;
					_prvt_calc_P_i_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+1.8600000000e+01)/1.0100000000e+01))));
					__K2_[12]= ((_prvt_calc_P_i_infinity-__OLD_[12])/_prvt_calc_tau_P_i);
					_prvt_aux_tol = fabs(__OLD_[12])*_prvt_rel_tol_;
					__TOL_[12] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[12] = fabs((_prvt_dtime/2) * (__K1_[12] - __K2_[12])/__TOL_[12]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[9])
				{
					_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
					_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
					__K2_[13]= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[13]))-(_prvt_calc_beta_xs*__OLD_[13]));
					_prvt_aux_tol = fabs(__OLD_[13])*_prvt_rel_tol_;
					__TOL_[13] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[13] = fabs((_prvt_dtime/2) * (__K1_[13] - __K2_[13])/__TOL_[13]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[10])
				{
					_prvt_calc_alpha_y = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6620000000e+01))))
:((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6630000000e+01))));
					_prvt_calc_beta_y = (1.0000000000e+00*exp(((__OLD_[0]+7.5130000000e+01)/2.1250000000e+01)));
					__K2_[14]= ((_prvt_calc_alpha_y*(1.0000000000e+00-__OLD_[14]))-(_prvt_calc_beta_y*__OLD_[14]));
					_prvt_aux_tol = fabs(__OLD_[14])*_prvt_rel_tol_;
					__TOL_[14] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[14] = fabs((_prvt_dtime/2) * (__K1_[14] - __K2_[14])/__TOL_[14]);
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
						__NEW_[5] = __K1_[5] * _prvt_dtime + __OLD_AUX_[5];
						__NEW_[7] = __K1_[7] * _prvt_dtime + __OLD_AUX_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[1] = __K1_[1] * _prvt_dtime + __OLD_AUX_[1];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[2] = __K1_[2] * _prvt_dtime + __OLD_AUX_[2];
						__NEW_[3] = __K1_[3] * _prvt_dtime + __OLD_AUX_[3];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[4] = __K1_[4] * _prvt_dtime + __OLD_AUX_[4];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[6] = __K1_[6] * _prvt_dtime + __OLD_AUX_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[8] = __K1_[8] * _prvt_dtime + __OLD_AUX_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[9] = __K1_[9] * _prvt_dtime + __OLD_AUX_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[10] = __K1_[10] * _prvt_dtime + __OLD_AUX_[10];
						__NEW_[11] = __K1_[11] * _prvt_dtime + __OLD_AUX_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[12] = __K1_[12] * _prvt_dtime + __OLD_AUX_[12];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[13] = __K1_[13] * _prvt_dtime + __OLD_AUX_[13];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[14] = __K1_[14] * _prvt_dtime + __OLD_AUX_[14];
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
							this->m_old_ = __OLD_AUX_[1];
							this->m_new_ = __OLD_[1];
							this->h1_old_ = __OLD_AUX_[2];
							this->h1_new_ = __OLD_[2];
							this->h2_old_ = __OLD_AUX_[3];
							this->h2_new_ = __OLD_[3];
							this->d_L_old_ = __OLD_AUX_[4];
							this->d_L_new_ = __OLD_[4];
							this->f_L_old_ = __OLD_AUX_[5];
							this->f_L_new_ = __OLD_[5];
							this->d_T_old_ = __OLD_AUX_[6];
							this->d_T_new_ = __OLD_[6];
							this->f_T_old_ = __OLD_AUX_[7];
							this->f_T_new_ = __OLD_[7];
							this->q_old_ = __OLD_AUX_[8];
							this->q_new_ = __OLD_[8];
							this->r_old_ = __OLD_AUX_[9];
							this->r_new_ = __OLD_[9];
							this->P_af_old_ = __OLD_AUX_[10];
							this->P_af_new_ = __OLD_[10];
							this->P_as_old_ = __OLD_AUX_[11];
							this->P_as_new_ = __OLD_[11];
							this->P_i_old_ = __OLD_AUX_[12];
							this->P_i_new_ = __OLD_[12];
							this->xs_old_ = __OLD_AUX_[13];
							this->xs_new_ = __OLD_[13];
							this->y_old_ = __OLD_AUX_[14];
							this->y_new_ = __OLD_[14];
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
						__NEW_[5] = __K2_[5] * _prvt_dtime + __OLD_[5];
						__NEW_[7] = __K2_[7] * _prvt_dtime + __OLD_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[1] = __K2_[1] * _prvt_dtime + __OLD_[1];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[2] = __K2_[2] * _prvt_dtime + __OLD_[2];
						__NEW_[3] = __K2_[3] * _prvt_dtime + __OLD_[3];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[4] = __K2_[4] * _prvt_dtime + __OLD_[4];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[6] = __K2_[6] * _prvt_dtime + __OLD_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[8] = __K2_[8] * _prvt_dtime + __OLD_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[9] = __K2_[9] * _prvt_dtime + __OLD_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[10] = __K2_[10] * _prvt_dtime + __OLD_[10];
						__NEW_[11] = __K2_[11] * _prvt_dtime + __OLD_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[12] = __K2_[12] * _prvt_dtime + __OLD_[12];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[13] = __K2_[13] * _prvt_dtime + __OLD_[13];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[14] = __K2_[14] * _prvt_dtime + __OLD_[14];
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
			double  _prvt_time = 0.0000000000e+00,  _prvt_dCell = 0.0000000000e+00,  _prvt_Version = 1.0000000000e+00,  _prvt_FCellConstant = 1.0309347000e+00,  _prvt_CmCentre = 2.0000000000e-05,  _prvt_CmPeriphery = 6.5000000000e-05,  _prvt_g_Na_Centre_Published = 0.0000000000e+00,  _prvt_g_Na_Periphery_Published = 1.2000000000e-06,  _prvt_g_Na_Centre_0DCapable = 0.0000000000e+00,  _prvt_g_Na_Periphery_0DCapable = 1.2040000000e-06,  _prvt_g_Na_Centre_1DCapable = 0.0000000000e+00,  _prvt_g_Na_Periphery_1DCapable = 3.7000000000e-07,  _prvt_Na_o = 1.4000000000e+02,  _prvt_F = 9.6845000000e+04,  _prvt_R = 8.3140000000e+03,  _prvt_T = 3.1000000000e+02,  _prvt_g_Ca_L_Centre_Published = 5.8000000000e-03,  _prvt_g_Ca_L_Periphery_Published = 6.5900000000e-02,  _prvt_g_Ca_L_Centre_0DCapable = 5.7938000000e-03,  _prvt_g_Ca_L_Periphery_0DCapable = 6.5886480000e-02,  _prvt_g_Ca_L_Centre_1DCapable = 8.2000000000e-03,  _prvt_g_Ca_L_Periphery_1DCapable = 6.5900000000e-02,  _prvt_E_Ca_L = 4.6400000000e+01,  _prvt_g_Ca_T_Centre_Published = 4.3000000000e-03,  _prvt_g_Ca_T_Periphery_Published = 1.3900000000e-02,  _prvt_g_Ca_T_Centre_0DCapable = 4.2780600000e-03,  _prvt_g_Ca_T_Periphery_0DCapable = 1.3882300000e-02,  _prvt_g_Ca_T_Centre_1DCapable = 2.1000000000e-03,  _prvt_g_Ca_T_Periphery_1DCapable = 6.9400000000e-03,  _prvt_E_Ca_T = 4.5000000000e+01,  _prvt_g_to_Centre_Published = 4.9100000000e-03,  _prvt_g_to_Periphery_Published = 3.6490000000e-02,  _prvt_g_to_Centre_0DCapable = 4.9050000000e-03,  _prvt_g_to_Periphery_0DCapable = 3.6495000000e-02,  _prvt_g_to_Centre_1DCapable = 4.9050000000e-03,  _prvt_g_to_Periphery_1DCapable = 3.6500000000e-02,  _prvt_g_sus_Centre_Published = 6.6500000000e-05,  _prvt_g_sus_Periphery_Published = 1.1400000000e-02,  _prvt_g_sus_Centre_0DCapable = 6.6455040000e-05,  _prvt_g_sus_Periphery_0DCapable = 1.1383760000e-02,  _prvt_g_sus_Centre_1DCapable = 2.6600000000e-04,  _prvt_g_sus_Periphery_1DCapable = 1.1400000000e-02,  _prvt_g_K_r_Centre_Published = 7.9700000000e-04,  _prvt_g_K_r_Periphery_Published = 1.6000000000e-02,  _prvt_g_K_r_Centre_0DCapable = 7.9704000000e-04,  _prvt_g_K_r_Periphery_0DCapable = 1.6000000000e-02,  _prvt_g_K_r_Centre_1DCapable = 7.3800000000e-04,  _prvt_g_K_r_Periphery_1DCapable = 2.0800000000e-02,  _prvt_g_K_s_Centre_Published = 5.1800000000e-04,  _prvt_g_K_s_Periphery_Published = 1.0400000000e-02,  _prvt_g_K_s_Centre_0DCapable = 3.4450000000e-04,  _prvt_g_K_s_Periphery_0DCapable = 1.0400000000e-02,  _prvt_g_K_s_Centre_1DCapable = 3.4500000000e-04,  _prvt_g_K_s_Periphery_1DCapable = 1.0400000000e-02,  _prvt_g_f_Na_Centre_Published = 5.4800000000e-04,  _prvt_g_f_Na_Periphery_Published = 6.9000000000e-03,  _prvt_g_f_Na_Centre_0DCapable = 5.4650000000e-04,  _prvt_g_f_Na_Periphery_0DCapable = 6.8750000000e-03,  _prvt_g_f_Na_Centre_1DCapable = 4.3700000000e-04,  _prvt_g_f_Na_Periphery_1DCapable = 5.5000000000e-03,  _prvt_g_f_K_Centre_Published = 5.4800000000e-04,  _prvt_g_f_K_Periphery_Published = 6.9000000000e-03,  _prvt_g_f_K_Centre_0DCapable = 5.4650000000e-04,  _prvt_g_f_K_Periphery_0DCapable = 6.8750000000e-03,  _prvt_g_f_K_Centre_1DCapable = 4.3700000000e-04,  _prvt_g_f_K_Periphery_1DCapable = 5.5000000000e-03,  _prvt_g_b_Na_Centre_Published = 5.8000000000e-05,  _prvt_g_b_Na_Periphery_Published = 1.8900000000e-04,  _prvt_g_b_Na_Centre_0DCapable = 5.8181800000e-05,  _prvt_g_b_Na_Periphery_0DCapable = 1.8880000000e-04,  _prvt_g_b_Na_Centre_1DCapable = 5.8000000000e-05,  _prvt_g_b_Na_Periphery_1DCapable = 1.8900000000e-04,  _prvt_g_b_K_Centre_Published = 2.5200000000e-05,  _prvt_g_b_K_Periphery_Published = 8.1900000000e-05,  _prvt_g_b_K_Centre_0DCapable = 2.5236360000e-05,  _prvt_g_b_K_Periphery_0DCapable = 8.1892000000e-05,  _prvt_g_b_K_Centre_1DCapable = 2.5200000000e-05,  _prvt_g_b_K_Periphery_1DCapable = 8.1900000000e-05,  _prvt_g_b_Ca_Centre_Published = 1.3200000000e-05,  _prvt_g_b_Ca_Periphery_Published = 4.3000000000e-05,  _prvt_g_b_Ca_Centre_0DCapable = 1.3236000000e-05,  _prvt_g_b_Ca_Periphery_0DCapable = 4.2952000000e-05,  _prvt_g_b_Ca_Centre_1DCapable = 1.3230000000e-05,  _prvt_g_b_Ca_Periphery_1DCapable = 4.2900000000e-05,  _prvt_k_NaCa_Centre_Published = 2.7000000000e-06,  _prvt_k_NaCa_Periphery_Published = 8.8000000000e-06,  _prvt_k_NaCa_Centre_0DCapable = 2.7229000000e-06,  _prvt_k_NaCa_Periphery_0DCapable = 8.8358400000e-06,  _prvt_k_NaCa_Centre_1DCapable = 2.8000000000e-06,  _prvt_k_NaCa_Periphery_1DCapable = 8.8000000000e-06,  _prvt_Na_i = 8.0000000000e+00,  _prvt_Ca_o = 2.0000000000e+00,  _prvt_gamma_NaCa = 5.0000000000e-01,  _prvt_Ca_i = 1.0000000000e-04,  _prvt_d_NaCa = 1.0000000000e-04,  _prvt_i_p_max_Centre_Published = 4.7800000000e-02,  _prvt_i_p_max_Periphery_Published = 1.6000000000e-01,  _prvt_i_p_max_Centre_0DCapable = 4.7825450000e-02,  _prvt_i_p_max_Periphery_0DCapable = 1.5519360000e-01,  _prvt_i_p_max_Centre_1DCapable = 4.7800000000e-02,  _prvt_i_p_max_Periphery_1DCapable = 1.6000000000e-01,  _prvt_K_m_Na = 5.6400000000e+00,  _prvt_K_o = 5.4000000000e+00,  _prvt_K_m_K = 6.2100000000e-01,  _prvt_i_Ca_p_max_Centre_Published = 0.0000000000e+00,  _prvt_i_Ca_p_max_Periphery_Published = 0.0000000000e+00,  _prvt_i_Ca_p_max_Centre_0DCapable = 0.0000000000e+00,  _prvt_i_Ca_p_max_Periphery_0DCapable = 0.0000000000e+00,  _prvt_i_Ca_p_max_Centre_1DCapable = 4.2000000000e-03,  _prvt_i_Ca_p_max_Periphery_1DCapable = 3.3390000000e-02,  _prvt_K_i = 1.4000000000e+02, 
			//private aux variables
			 _prvt_calc_FCell=0.0,  _prvt_calc_Cm=0.0,  _prvt_calc_g_Na=0.0,  _prvt_calc_i_Na=0.0,  _prvt_calc_m_infinity=0.0,  _prvt_calc_tau_m=0.0,  _prvt_calc_F_Na=0.0,  _prvt_calc_h=0.0,  _prvt_calc_h1_infinity=0.0,  _prvt_calc_h2_infinity=0.0,  _prvt_calc_tau_h1=0.0,  _prvt_calc_tau_h2=0.0,  _prvt_calc_g_Ca_L=0.0,  _prvt_calc_i_Ca_L=0.0,  _prvt_calc_alpha_d_L=0.0,  _prvt_calc_beta_d_L=0.0,  _prvt_calc_tau_d_L=0.0,  _prvt_calc_d_L_infinity=0.0,  _prvt_calc_alpha_f_L=0.0,  _prvt_calc_beta_f_L=0.0,  _prvt_calc_tau_f_L=0.0,  _prvt_calc_f_L_infinity=0.0,  _prvt_calc_g_Ca_T=0.0,  _prvt_calc_i_Ca_T=0.0,  _prvt_calc_alpha_d_T=0.0,  _prvt_calc_beta_d_T=0.0,  _prvt_calc_tau_d_T=0.0,  _prvt_calc_d_T_infinity=0.0,  _prvt_calc_alpha_f_T=0.0,  _prvt_calc_beta_f_T=0.0,  _prvt_calc_tau_f_T=0.0,  _prvt_calc_f_T_infinity=0.0,  _prvt_calc_g_to=0.0,  _prvt_calc_g_sus=0.0,  _prvt_calc_i_to=0.0,  _prvt_calc_i_sus=0.0,  _prvt_calc_q_infinity=0.0,  _prvt_calc_tau_q=0.0,  _prvt_calc_r_infinity=0.0,  _prvt_calc_tau_r=0.0,  _prvt_calc_g_K_r=0.0,  _prvt_calc_i_K_r=0.0,  _prvt_calc_P_a=0.0,  _prvt_calc_P_af_infinity=0.0,  _prvt_calc_tau_P_af=0.0,  _prvt_calc_P_as_infinity=0.0,  _prvt_calc_tau_P_as=0.0,  _prvt_calc_tau_P_i=0.0,  _prvt_calc_P_i_infinity=0.0,  _prvt_calc_g_K_s=0.0,  _prvt_calc_i_K_s=0.0,  _prvt_calc_alpha_xs=0.0,  _prvt_calc_beta_xs=0.0,  _prvt_calc_g_f_Na=0.0,  _prvt_calc_i_f_Na=0.0,  _prvt_calc_g_f_K=0.0,  _prvt_calc_i_f_K=0.0,  _prvt_calc_alpha_y=0.0,  _prvt_calc_beta_y=0.0,  _prvt_calc_g_b_Na=0.0,  _prvt_calc_i_b_Na=0.0,  _prvt_calc_g_b_K=0.0,  _prvt_calc_i_b_K=0.0,  _prvt_calc_g_b_Ca=0.0,  _prvt_calc_i_b_Ca=0.0,  _prvt_calc_k_NaCa=0.0,  _prvt_calc_i_NaCa=0.0,  _prvt_calc_i_p_max=0.0,  _prvt_calc_i_p=0.0,  _prvt_calc_i_Ca_p_max=0.0,  _prvt_calc_i_Ca_p=0.0,  _prvt_calc_E_Na=0.0,  _prvt_calc_E_K=0.0,  _prvt_calc_E_Ca=0.0,  _prvt_calc_E_K_s=0.0, 
			//private right hand side variables
			 _prvt_V_lado_direito_,  _prvt_m_lado_direito_,  _prvt_h1_lado_direito_,  _prvt_h2_lado_direito_,  _prvt_d_L_lado_direito_,  _prvt_f_L_lado_direito_,  _prvt_d_T_lado_direito_,  _prvt_f_T_lado_direito_,  _prvt_q_lado_direito_,  _prvt_r_lado_direito_,  _prvt_P_af_lado_direito_,  _prvt_P_as_lado_direito_,  _prvt_P_i_lado_direito_,  _prvt_xs_lado_direito_,  _prvt_y_lado_direito_, 
			//private time variables
			_prvt_time_new = this->time_new,
			_prvt_dtime = this->dtime,
			_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;
			__NEW_[0] = __OLD_[0] = -3.9013558536e+01;
			__NEW_[1] = __OLD_[1] = 9.2361701692e-02;
			__NEW_[2] = __OLD_[2] = 1.5905380261e-02;
			__NEW_[3] = __OLD_[3] = 1.4452161090e-02;
			__NEW_[4] = __OLD_[4] = 4.8049008950e-02;
			__NEW_[5] = __OLD_[5] = 4.8779845203e-01;
			__NEW_[6] = __OLD_[6] = 4.2074047435e-01;
			__NEW_[7] = __OLD_[7] = 3.8968420558e-02;
			__NEW_[8] = __OLD_[8] = 2.9760539675e-01;
			__NEW_[9] = __OLD_[9] = 6.4402950262e-02;
			__NEW_[10] = __OLD_[10] = 1.3034201158e-01;
			__NEW_[11] = __OLD_[11] = 4.6960956028e-01;
			__NEW_[12] = __OLD_[12] = 8.7993375273e-01;
			__NEW_[13] = __OLD_[13] = 8.2293827208e-02;
			__NEW_[14] = __OLD_[14] = 3.8892917590e-02;
			_prvt_time_new += _prvt_dtime;
			if(omp_get_thread_num()==_prvt_tree_thread[0])
			{
				_prvt_calc_FCell = ((_prvt_Version==0.0000000000e+00))
?(((1.0700000000e+00*((3.0000000000e+00*_prvt_dCell)-1.0000000000e-01))/(3.0000000000e+00*(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01)))))))
:(((_prvt_Version==1.0000000000e+00))
?(((_prvt_FCellConstant*_prvt_dCell)/(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01))))))
:(((1.0700000000e+00*2.9000000000e+01*_prvt_dCell)/(3.0000000000e+01*(1.0000000000e+00+(7.7450000000e-01*exp(((-((2.9000000000e+01*_prvt_dCell)-2.4500000000e+01))/1.9500000000e+00))))))))
;
				_prvt_calc_F_Na = ((_prvt_Version==0.0000000000e+00))
?((((9.5200000000e-02*exp(((-6.3000000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6600000000e+00*exp(((-2.2500000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6900000000e-02))
:((((9.5180000000e-02*exp(((-6.3060000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6620000000e+00*exp(((-2.2510000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6930000000e-02));
				_prvt_calc_alpha_f_L = ((_prvt_Version==1.0000000000e+00))
?(((3.7500000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)))
:(((3.1200000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
				_prvt_calc_beta_f_L = ((_prvt_Version==1.0000000000e+00))
?((3.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))))
:((2.5000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))));
				_prvt_calc_f_L_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5000000000e+01)/5.0000000000e+00))));
				_prvt_calc_beta_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5000000000e+01*exp(((__OLD_[0]+7.1000000000e+01)/1.5380000000e+01))))
:((1.5000000000e+01*exp(((__OLD_[0]+7.1700000000e+01)/1.5380000000e+01))));
				_prvt_calc_f_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+7.1000000000e+01)/9.0000000000e+00))));
				_prvt_calc_P_a = ((6.0000000000e-01*__OLD_[10])+(4.0000000000e-01*__OLD_[11]));
				_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/_prvt_Na_i)));
				_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/_prvt_K_i)));
				_prvt_calc_E_Ca = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Ca_o/_prvt_Ca_i)));
				_prvt_calc_E_K_s = ((_prvt_Version==0.0000000000e+00))
?((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(1.2000000000e-01*_prvt_Na_o))/(_prvt_K_i+(1.2000000000e-01*_prvt_Na_i))))))
:((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(3.0000000000e-02*_prvt_Na_o))/(_prvt_K_i+(3.0000000000e-02*_prvt_Na_i))))));
				_prvt_calc_Cm = (_prvt_CmCentre+(_prvt_calc_FCell*(_prvt_CmPeriphery-_prvt_CmCentre)));
				_prvt_calc_g_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_Published-_prvt_g_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_0DCapable-_prvt_g_Na_Centre_0DCapable))))
:((_prvt_g_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_1DCapable-_prvt_g_Na_Centre_1DCapable)))))
;
				_prvt_calc_h = (((1.0000000000e+00-_prvt_calc_F_Na)*__OLD_[2])+(_prvt_calc_F_Na*__OLD_[3]));
				_prvt_calc_g_Ca_L = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_L_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_Published-_prvt_g_Ca_L_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_L_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_0DCapable-_prvt_g_Ca_L_Centre_0DCapable))))
:((_prvt_g_Ca_L_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_1DCapable-_prvt_g_Ca_L_Centre_1DCapable)))))
;
				_prvt_calc_tau_f_L = ((_prvt_Version==1.0000000000e+00))
?(((1.2000000000e+00-(2.0000000000e-01*_prvt_calc_FCell))/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)))
:((1.0000000000e+00/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)));
				_prvt_calc_g_Ca_T = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_T_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_Published-_prvt_g_Ca_T_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_T_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_0DCapable-_prvt_g_Ca_T_Centre_0DCapable))))
:((_prvt_g_Ca_T_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_1DCapable-_prvt_g_Ca_T_Centre_1DCapable)))))
;
				_prvt_calc_alpha_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5300000000e+01*exp(((-(__OLD_[0]+7.1000000000e+01+(7.0000000000e-01*_prvt_calc_FCell)))/8.3300000000e+01))))
:((1.5300000000e+01*exp(((-(__OLD_[0]+7.1700000000e+01))/8.3300000000e+01))));
				_prvt_calc_g_to = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_to_Centre_Published+(_prvt_calc_FCell*(_prvt_g_to_Periphery_Published-_prvt_g_to_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_to_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_0DCapable-_prvt_g_to_Centre_0DCapable))))
:((_prvt_g_to_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_1DCapable-_prvt_g_to_Centre_1DCapable)))))
;
				_prvt_calc_g_sus = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_sus_Centre_Published+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_Published-_prvt_g_sus_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_sus_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_0DCapable-_prvt_g_sus_Centre_0DCapable))))
:((_prvt_g_sus_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_1DCapable-_prvt_g_sus_Centre_1DCapable)))))
;
				_prvt_calc_g_K_r = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_r_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_Published-_prvt_g_K_r_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_r_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_0DCapable-_prvt_g_K_r_Centre_0DCapable))))
:((_prvt_g_K_r_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_1DCapable-_prvt_g_K_r_Centre_1DCapable)))))
;
				_prvt_calc_g_K_s = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_s_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_Published-_prvt_g_K_s_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_s_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_0DCapable-_prvt_g_K_s_Centre_0DCapable))))
:((_prvt_g_K_s_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_1DCapable-_prvt_g_K_s_Centre_1DCapable)))))
;
				_prvt_calc_g_f_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_Published-_prvt_g_f_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_0DCapable-_prvt_g_f_Na_Centre_0DCapable))))
:((_prvt_g_f_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_1DCapable-_prvt_g_f_Na_Centre_1DCapable)))))
;
				_prvt_calc_g_f_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_Published-_prvt_g_f_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_0DCapable-_prvt_g_f_K_Centre_0DCapable))))
:((_prvt_g_f_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_1DCapable-_prvt_g_f_K_Centre_1DCapable)))))
;
				_prvt_calc_g_b_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_Published-_prvt_g_b_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_0DCapable-_prvt_g_b_Na_Centre_0DCapable))))
:((_prvt_g_b_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_1DCapable-_prvt_g_b_Na_Centre_1DCapable)))))
;
				_prvt_calc_g_b_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_Published-_prvt_g_b_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_0DCapable-_prvt_g_b_K_Centre_0DCapable))))
:((_prvt_g_b_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_1DCapable-_prvt_g_b_K_Centre_1DCapable)))))
;
				_prvt_calc_g_b_Ca = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Ca_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_Published-_prvt_g_b_Ca_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Ca_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_0DCapable-_prvt_g_b_Ca_Centre_0DCapable))))
:((_prvt_g_b_Ca_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_1DCapable-_prvt_g_b_Ca_Centre_1DCapable)))))
;
				_prvt_calc_k_NaCa = ((_prvt_Version==0.0000000000e+00))
?((_prvt_k_NaCa_Centre_Published+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_Published-_prvt_k_NaCa_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_k_NaCa_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_0DCapable-_prvt_k_NaCa_Centre_0DCapable))))
:((_prvt_k_NaCa_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_1DCapable-_prvt_k_NaCa_Centre_1DCapable)))))
;
				_prvt_calc_i_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_Published-_prvt_i_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_0DCapable-_prvt_i_p_max_Centre_0DCapable))))
:((_prvt_i_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_1DCapable-_prvt_i_p_max_Centre_1DCapable)))))
;
				_prvt_calc_i_Ca_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_Published-_prvt_i_Ca_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_0DCapable-_prvt_i_Ca_p_max_Centre_0DCapable))))
:((_prvt_i_Ca_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_1DCapable-_prvt_i_Ca_p_max_Centre_1DCapable)))))
;
				_prvt_calc_i_Ca_L = (_prvt_calc_g_Ca_L*((__OLD_[5]*__OLD_[4])+(6.0000000000e-03/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4100000000e+01))/6.0000000000e+00)))))*(__OLD_[0]-_prvt_E_Ca_L));
				_prvt_calc_i_Ca_T = (_prvt_calc_g_Ca_T*__OLD_[6]*__OLD_[7]*(__OLD_[0]-_prvt_E_Ca_T));
				_prvt_calc_tau_f_T = (1.0000000000e+00/(_prvt_calc_alpha_f_T+_prvt_calc_beta_f_T));
				_prvt_calc_i_to = (_prvt_calc_g_to*__OLD_[8]*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_sus = (_prvt_calc_g_sus*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_K_r = (_prvt_calc_g_K_r*_prvt_calc_P_a*__OLD_[12]*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_K_s = (_prvt_calc_g_K_s*pow(__OLD_[13],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K_s));
				_prvt_calc_i_f_Na = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_Na)))
:((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-7.7600000000e+01)));
				_prvt_calc_i_f_K = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_K)))
:((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]+1.0200000000e+02)));
				_prvt_calc_i_b_Na = (_prvt_calc_g_b_Na*(__OLD_[0]-_prvt_calc_E_Na));
				_prvt_calc_i_b_K = (_prvt_calc_g_b_K*(__OLD_[0]-_prvt_calc_E_K));
				_prvt_calc_i_b_Ca = (_prvt_calc_g_b_Ca*(__OLD_[0]-_prvt_calc_E_Ca));
				_prvt_calc_i_NaCa = ((_prvt_Version==0.0000000000e+00))
?(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7400000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))))
:(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7430000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))));
				_prvt_calc_i_p = ((_prvt_calc_i_p_max*pow((_prvt_Na_i/(_prvt_K_m_Na+_prvt_Na_i)),3.0000000000e+00)*pow((_prvt_K_o/(_prvt_K_m_K+_prvt_K_o)),2.0000000000e+00)*1.6000000000e+00)/(1.5000000000e+00+exp(((-(__OLD_[0]+6.0000000000e+01))/4.0000000000e+01))));
				_prvt_calc_i_Ca_p = ((_prvt_calc_i_Ca_p_max*_prvt_Ca_i)/(_prvt_Ca_i+4.0000000000e-04));
				_prvt_calc_i_Na = (((((_prvt_calc_g_Na*pow(__OLD_[1],3.0000000000e+00)*_prvt_calc_h*_prvt_Na_o*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*(exp((((__OLD_[0]-_prvt_calc_E_Na)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))/(exp(((__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))*__OLD_[0]);
				__K1_[0]= (((-1.0000000000e+00)/_prvt_calc_Cm)*(_prvt_calc_i_Na+_prvt_calc_i_Ca_L+_prvt_calc_i_Ca_T+_prvt_calc_i_to+_prvt_calc_i_sus+_prvt_calc_i_K_r+_prvt_calc_i_K_s+_prvt_calc_i_f_Na+_prvt_calc_i_f_K+_prvt_calc_i_b_Na+_prvt_calc_i_b_Ca+_prvt_calc_i_b_K+_prvt_calc_i_NaCa+_prvt_calc_i_p+_prvt_calc_i_Ca_p));
				__NEW_[0]= __K1_[0] * _prvt_dtime + __OLD_[0];
				__K1_[5]= ((_prvt_calc_f_L_infinity-__OLD_[5])/_prvt_calc_tau_f_L);
				__NEW_[5]= __K1_[5] * _prvt_dtime + __OLD_[5];
				__K1_[7]= ((_prvt_calc_f_T_infinity-__OLD_[7])/_prvt_calc_tau_f_T);
				__NEW_[7]= __K1_[7] * _prvt_dtime + __OLD_[7];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[1])
			{
				_prvt_calc_m_infinity = ((_prvt_Version==0.0000000000e+00))
?(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-__OLD_[0])/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)))
:(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.0320000000e+01))/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
				_prvt_calc_tau_m = ((_prvt_Version==0.0000000000e+00))
?(((6.2470000000e-04/((8.3200000000e-01*exp(((-3.3500000000e-01)*(__OLD_[0]+5.6700000000e+01))))+(6.2700000000e-01*exp((8.2000000000e-02*(__OLD_[0]+6.5010000000e+01))))))+4.0000000000e-05))
:(((6.2470000000e-04/((8.3221660000e-01*exp(((-3.3566000000e-01)*(__OLD_[0]+5.6706200000e+01))))+(6.2740000000e-01*exp((8.2300000000e-02*(__OLD_[0]+6.5013100000e+01))))))+4.5690000000e-05));
				__K1_[1]= ((_prvt_calc_m_infinity-__OLD_[1])/_prvt_calc_tau_m);
				__NEW_[1]= __K1_[1] * _prvt_dtime + __OLD_[1];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[2])
			{
				_prvt_calc_h1_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+6.6100000000e+01)/6.4000000000e+00))));
				_prvt_calc_tau_h1 = (((3.7170000000e-06*exp(((-2.8150000000e-01)*(__OLD_[0]+1.7110000000e+01))))/(1.0000000000e+00+(3.7320000000e-03*exp(((-3.4260000000e-01)*(__OLD_[0]+3.7760000000e+01))))))+5.9770000000e-04);
				_prvt_calc_tau_h2 = (((3.1860000000e-08*exp(((-6.2190000000e-01)*(__OLD_[0]+1.8800000000e+01))))/(1.0000000000e+00+(7.1890000000e-05*exp(((-6.6830000000e-01)*(__OLD_[0]+3.4070000000e+01))))))+3.5560000000e-03);
				_prvt_calc_h2_infinity = _prvt_calc_h1_infinity;
				__K1_[2]= ((_prvt_calc_h1_infinity-__OLD_[2])/_prvt_calc_tau_h1);
				__NEW_[2]= __K1_[2] * _prvt_dtime + __OLD_[2];
				__K1_[3]= ((_prvt_calc_h2_infinity-__OLD_[3])/_prvt_calc_tau_h2);
				__NEW_[3]= __K1_[3] * _prvt_dtime + __OLD_[3];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[3])
			{
				_prvt_calc_alpha_d_L = ((_prvt_Version==0.0000000000e+00))
?(((((-2.8380000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((_prvt_Version==1.0000000000e+00))
?(((((-2.8390000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((((-2.8400000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00)))))
;
				_prvt_calc_beta_d_L = ((_prvt_Version==1.0000000000e+00))
?(((1.1430000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)))
:(((1.1420000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)));
				_prvt_calc_d_L_infinity = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.3100000000e+01))/6.0000000000e+00)))))
:(((_prvt_Version==1.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2300000000e+01+(8.0000000000e-01*_prvt_calc_FCell)))/6.0000000000e+00)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2200000000e+01))/6.0000000000e+00))))))
;
				_prvt_calc_tau_d_L = (2.0000000000e+00/(_prvt_calc_alpha_d_L+_prvt_calc_beta_d_L));
				__K1_[4]= ((_prvt_calc_d_L_infinity-__OLD_[4])/_prvt_calc_tau_d_L);
				__NEW_[4]= __K1_[4] * _prvt_dtime + __OLD_[4];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[4])
			{
				_prvt_calc_alpha_d_T = (1.0680000000e+03*exp(((__OLD_[0]+2.6300000000e+01)/3.0000000000e+01)));
				_prvt_calc_beta_d_T = (1.0680000000e+03*exp(((-(__OLD_[0]+2.6300000000e+01))/3.0000000000e+01)));
				_prvt_calc_d_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.7000000000e+01))/6.8000000000e+00))));
				_prvt_calc_tau_d_T = (1.0000000000e+00/(_prvt_calc_alpha_d_T+_prvt_calc_beta_d_T));
				__K1_[6]= ((_prvt_calc_d_T_infinity-__OLD_[6])/_prvt_calc_tau_d_T);
				__NEW_[6]= __K1_[6] * _prvt_dtime + __OLD_[6];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[5])
			{
				_prvt_calc_q_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.9370000000e+01)/1.3100000000e+01))));
				_prvt_calc_tau_q = ((_prvt_Version==0.0000000000e+00))
?((1.0100000000e-02+(6.5170000000e-02/(5.7000000000e-01*exp(((-8.0000000000e-02)*(__OLD_[0]+4.9000000000e+01)))))+(2.4000000000e-05*exp((1.0000000000e-01*(__OLD_[0]+5.0930000000e+01))))))
:(((_prvt_Version==1.0000000000e+00))
?(((1.0000000000e-03/3.0000000000e+00)*(3.0310000000e+01+(1.9550000000e+02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01+(1.0000000000e+01*_prvt_calc_FCell)))))+(7.1740000000e-01*exp(((2.7190000000e-01-(1.7190000000e-01*_prvt_calc_FCell))*1.0000000000e+00*(__OLD_[0]+4.0930000000e+01+(1.0000000000e+01*_prvt_calc_FCell))))))))))
:((1.0100000000e-02+(6.5170000000e-02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01))))+(7.1740000000e-01*exp((2.7190000000e-01*(__OLD_[0]+4.0930000000e+01)))))))))
;
				__K1_[8]= ((_prvt_calc_q_infinity-__OLD_[8])/_prvt_calc_tau_q);
				__NEW_[8]= __K1_[8] * _prvt_dtime + __OLD_[8];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[6])
			{
				_prvt_calc_r_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]-1.0930000000e+01))/1.9700000000e+01))));
				_prvt_calc_tau_r = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e-03*(2.9800000000e+00+(1.5590000000e+01/((1.0370000000e+00*exp((9.0000000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.2000000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:(((_prvt_Version==1.0000000000e+00))
?((2.5000000000e-03*(1.1910000000e+00+(7.8380000000e+00/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:((1.0000000000e-03*(2.9800000000e+00+(1.9590000000e+01/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01))))))))))
;
				__K1_[9]= ((_prvt_calc_r_infinity-__OLD_[9])/_prvt_calc_tau_r);
				__NEW_[9]= __K1_[9] * _prvt_dtime + __OLD_[9];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[7])
			{
				_prvt_calc_P_af_infinity = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4200000000e+01))/1.0600000000e+01)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.3200000000e+01))/1.0600000000e+01)))));
				_prvt_calc_tau_P_af = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-9.0000000000e+00)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.2500000000e+01))))))
:((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-1.0000000000e+01)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.2500000000e+01))))));
				_prvt_calc_tau_P_as = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-9.0000000000e+00)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.1600000000e+01))))))
:((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-1.0000000000e+01)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.1600000000e+01))))));
				_prvt_calc_P_as_infinity = _prvt_calc_P_af_infinity;
				__K1_[10]= ((_prvt_calc_P_af_infinity-__OLD_[10])/_prvt_calc_tau_P_af);
				__NEW_[10]= __K1_[10] * _prvt_dtime + __OLD_[10];
				__K1_[11]= ((_prvt_calc_P_as_infinity-__OLD_[11])/_prvt_calc_tau_P_as);
				__NEW_[11]= __K1_[11] * _prvt_dtime + __OLD_[11];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[8])
			{
				_prvt_calc_tau_P_i = ((_prvt_Version==0.0000000000e+00))
?(2.0000000000e-03)
:(((_prvt_Version==1.0000000000e+00))
?(2.0000000000e-03)
:(6.0000000000e-03))
;
				_prvt_calc_P_i_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+1.8600000000e+01)/1.0100000000e+01))));
				__K1_[12]= ((_prvt_calc_P_i_infinity-__OLD_[12])/_prvt_calc_tau_P_i);
				__NEW_[12]= __K1_[12] * _prvt_dtime + __OLD_[12];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[9])
			{
				_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
				_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
				__K1_[13]= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[13]))-(_prvt_calc_beta_xs*__OLD_[13]));
				__NEW_[13]= __K1_[13] * _prvt_dtime + __OLD_[13];
			}
			if(omp_get_thread_num()==_prvt_tree_thread[10])
			{
				_prvt_calc_alpha_y = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6620000000e+01))))
:((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6630000000e+01))));
				_prvt_calc_beta_y = (1.0000000000e+00*exp(((__OLD_[0]+7.5130000000e+01)/2.1250000000e+01)));
				__K1_[14]= ((_prvt_calc_alpha_y*(1.0000000000e+00-__OLD_[14]))-(_prvt_calc_beta_y*__OLD_[14]));
				__NEW_[14]= __K1_[14] * _prvt_dtime + __OLD_[14];
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
					_prvt_calc_FCell = ((_prvt_Version==0.0000000000e+00))
?(((1.0700000000e+00*((3.0000000000e+00*_prvt_dCell)-1.0000000000e-01))/(3.0000000000e+00*(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01)))))))
:(((_prvt_Version==1.0000000000e+00))
?(((_prvt_FCellConstant*_prvt_dCell)/(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*_prvt_dCell)-2.0500000000e+00))/2.9500000000e-01))))))
:(((1.0700000000e+00*2.9000000000e+01*_prvt_dCell)/(3.0000000000e+01*(1.0000000000e+00+(7.7450000000e-01*exp(((-((2.9000000000e+01*_prvt_dCell)-2.4500000000e+01))/1.9500000000e+00))))))))
;
					_prvt_calc_F_Na = ((_prvt_Version==0.0000000000e+00))
?((((9.5200000000e-02*exp(((-6.3000000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6600000000e+00*exp(((-2.2500000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6900000000e-02))
:((((9.5180000000e-02*exp(((-6.3060000000e-02)*(__OLD_[0]+3.4400000000e+01))))/(1.0000000000e+00+(1.6620000000e+00*exp(((-2.2510000000e-01)*(__OLD_[0]+6.3700000000e+01))))))+8.6930000000e-02));
					_prvt_calc_alpha_f_L = ((_prvt_Version==1.0000000000e+00))
?(((3.7500000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)))
:(((3.1200000000e+00*(__OLD_[0]+2.8000000000e+01))/(exp(((__OLD_[0]+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
					_prvt_calc_beta_f_L = ((_prvt_Version==1.0000000000e+00))
?((3.0000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))))
:((2.5000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]+2.8000000000e+01))/4.0000000000e+00)))));
					_prvt_calc_f_L_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+4.5000000000e+01)/5.0000000000e+00))));
					_prvt_calc_beta_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5000000000e+01*exp(((__OLD_[0]+7.1000000000e+01)/1.5380000000e+01))))
:((1.5000000000e+01*exp(((__OLD_[0]+7.1700000000e+01)/1.5380000000e+01))));
					_prvt_calc_f_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+7.1000000000e+01)/9.0000000000e+00))));
					_prvt_calc_P_a = ((6.0000000000e-01*__OLD_[10])+(4.0000000000e-01*__OLD_[11]));
					_prvt_calc_E_Na = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_Na_o/_prvt_Na_i)));
					_prvt_calc_E_K = (((_prvt_R*_prvt_T)/_prvt_F)*log((_prvt_K_o/_prvt_K_i)));
					_prvt_calc_E_Ca = (((_prvt_R*_prvt_T)/(2.0000000000e+00*_prvt_F))*log((_prvt_Ca_o/_prvt_Ca_i)));
					_prvt_calc_E_K_s = ((_prvt_Version==0.0000000000e+00))
?((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(1.2000000000e-01*_prvt_Na_o))/(_prvt_K_i+(1.2000000000e-01*_prvt_Na_i))))))
:((((_prvt_R*_prvt_T)/_prvt_F)*log(((_prvt_K_o+(3.0000000000e-02*_prvt_Na_o))/(_prvt_K_i+(3.0000000000e-02*_prvt_Na_i))))));
					_prvt_calc_Cm = (_prvt_CmCentre+(_prvt_calc_FCell*(_prvt_CmPeriphery-_prvt_CmCentre)));
					_prvt_calc_g_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_Published-_prvt_g_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_0DCapable-_prvt_g_Na_Centre_0DCapable))))
:((_prvt_g_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Na_Periphery_1DCapable-_prvt_g_Na_Centre_1DCapable)))))
;
					_prvt_calc_h = (((1.0000000000e+00-_prvt_calc_F_Na)*__OLD_[2])+(_prvt_calc_F_Na*__OLD_[3]));
					_prvt_calc_g_Ca_L = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_L_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_Published-_prvt_g_Ca_L_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_L_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_0DCapable-_prvt_g_Ca_L_Centre_0DCapable))))
:((_prvt_g_Ca_L_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_L_Periphery_1DCapable-_prvt_g_Ca_L_Centre_1DCapable)))))
;
					_prvt_calc_tau_f_L = ((_prvt_Version==1.0000000000e+00))
?(((1.2000000000e+00-(2.0000000000e-01*_prvt_calc_FCell))/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)))
:((1.0000000000e+00/(_prvt_calc_alpha_f_L+_prvt_calc_beta_f_L)));
					_prvt_calc_g_Ca_T = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_Ca_T_Centre_Published+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_Published-_prvt_g_Ca_T_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_Ca_T_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_0DCapable-_prvt_g_Ca_T_Centre_0DCapable))))
:((_prvt_g_Ca_T_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_Ca_T_Periphery_1DCapable-_prvt_g_Ca_T_Centre_1DCapable)))))
;
					_prvt_calc_alpha_f_T = ((_prvt_Version==1.0000000000e+00))
?((1.5300000000e+01*exp(((-(__OLD_[0]+7.1000000000e+01+(7.0000000000e-01*_prvt_calc_FCell)))/8.3300000000e+01))))
:((1.5300000000e+01*exp(((-(__OLD_[0]+7.1700000000e+01))/8.3300000000e+01))));
					_prvt_calc_g_to = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_to_Centre_Published+(_prvt_calc_FCell*(_prvt_g_to_Periphery_Published-_prvt_g_to_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_to_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_0DCapable-_prvt_g_to_Centre_0DCapable))))
:((_prvt_g_to_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_to_Periphery_1DCapable-_prvt_g_to_Centre_1DCapable)))))
;
					_prvt_calc_g_sus = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_sus_Centre_Published+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_Published-_prvt_g_sus_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_sus_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_0DCapable-_prvt_g_sus_Centre_0DCapable))))
:((_prvt_g_sus_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_sus_Periphery_1DCapable-_prvt_g_sus_Centre_1DCapable)))))
;
					_prvt_calc_g_K_r = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_r_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_Published-_prvt_g_K_r_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_r_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_0DCapable-_prvt_g_K_r_Centre_0DCapable))))
:((_prvt_g_K_r_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_r_Periphery_1DCapable-_prvt_g_K_r_Centre_1DCapable)))))
;
					_prvt_calc_g_K_s = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_K_s_Centre_Published+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_Published-_prvt_g_K_s_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_K_s_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_0DCapable-_prvt_g_K_s_Centre_0DCapable))))
:((_prvt_g_K_s_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_K_s_Periphery_1DCapable-_prvt_g_K_s_Centre_1DCapable)))))
;
					_prvt_calc_g_f_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_Published-_prvt_g_f_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_0DCapable-_prvt_g_f_Na_Centre_0DCapable))))
:((_prvt_g_f_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_Na_Periphery_1DCapable-_prvt_g_f_Na_Centre_1DCapable)))))
;
					_prvt_calc_g_f_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_f_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_Published-_prvt_g_f_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_f_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_0DCapable-_prvt_g_f_K_Centre_0DCapable))))
:((_prvt_g_f_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_f_K_Periphery_1DCapable-_prvt_g_f_K_Centre_1DCapable)))))
;
					_prvt_calc_g_b_Na = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Na_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_Published-_prvt_g_b_Na_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Na_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_0DCapable-_prvt_g_b_Na_Centre_0DCapable))))
:((_prvt_g_b_Na_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Na_Periphery_1DCapable-_prvt_g_b_Na_Centre_1DCapable)))))
;
					_prvt_calc_g_b_K = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_K_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_Published-_prvt_g_b_K_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_K_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_0DCapable-_prvt_g_b_K_Centre_0DCapable))))
:((_prvt_g_b_K_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_K_Periphery_1DCapable-_prvt_g_b_K_Centre_1DCapable)))))
;
					_prvt_calc_g_b_Ca = ((_prvt_Version==0.0000000000e+00))
?((_prvt_g_b_Ca_Centre_Published+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_Published-_prvt_g_b_Ca_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_g_b_Ca_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_0DCapable-_prvt_g_b_Ca_Centre_0DCapable))))
:((_prvt_g_b_Ca_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_g_b_Ca_Periphery_1DCapable-_prvt_g_b_Ca_Centre_1DCapable)))))
;
					_prvt_calc_k_NaCa = ((_prvt_Version==0.0000000000e+00))
?((_prvt_k_NaCa_Centre_Published+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_Published-_prvt_k_NaCa_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_k_NaCa_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_0DCapable-_prvt_k_NaCa_Centre_0DCapable))))
:((_prvt_k_NaCa_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_k_NaCa_Periphery_1DCapable-_prvt_k_NaCa_Centre_1DCapable)))))
;
					_prvt_calc_i_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_Published-_prvt_i_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_0DCapable-_prvt_i_p_max_Centre_0DCapable))))
:((_prvt_i_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_p_max_Periphery_1DCapable-_prvt_i_p_max_Centre_1DCapable)))))
;
					_prvt_calc_i_Ca_p_max = ((_prvt_Version==0.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_Published+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_Published-_prvt_i_Ca_p_max_Centre_Published))))
:(((_prvt_Version==1.0000000000e+00))
?((_prvt_i_Ca_p_max_Centre_0DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_0DCapable-_prvt_i_Ca_p_max_Centre_0DCapable))))
:((_prvt_i_Ca_p_max_Centre_1DCapable+(_prvt_calc_FCell*(_prvt_i_Ca_p_max_Periphery_1DCapable-_prvt_i_Ca_p_max_Centre_1DCapable)))))
;
					_prvt_calc_i_Ca_L = (_prvt_calc_g_Ca_L*((__OLD_[5]*__OLD_[4])+(6.0000000000e-03/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4100000000e+01))/6.0000000000e+00)))))*(__OLD_[0]-_prvt_E_Ca_L));
					_prvt_calc_i_Ca_T = (_prvt_calc_g_Ca_T*__OLD_[6]*__OLD_[7]*(__OLD_[0]-_prvt_E_Ca_T));
					_prvt_calc_tau_f_T = (1.0000000000e+00/(_prvt_calc_alpha_f_T+_prvt_calc_beta_f_T));
					_prvt_calc_i_to = (_prvt_calc_g_to*__OLD_[8]*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_sus = (_prvt_calc_g_sus*__OLD_[9]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K_r = (_prvt_calc_g_K_r*_prvt_calc_P_a*__OLD_[12]*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_K_s = (_prvt_calc_g_K_s*pow(__OLD_[13],2.0000000000e+00)*(__OLD_[0]-_prvt_calc_E_K_s));
					_prvt_calc_i_f_Na = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_Na)))
:((_prvt_calc_g_f_Na*__OLD_[14]*(__OLD_[0]-7.7600000000e+01)));
					_prvt_calc_i_f_K = ((_prvt_Version!=2.0000000000e+00))
?((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]-_prvt_calc_E_K)))
:((_prvt_calc_g_f_K*__OLD_[14]*(__OLD_[0]+1.0200000000e+02)));
					_prvt_calc_i_b_Na = (_prvt_calc_g_b_Na*(__OLD_[0]-_prvt_calc_E_Na));
					_prvt_calc_i_b_K = (_prvt_calc_g_b_K*(__OLD_[0]-_prvt_calc_E_K));
					_prvt_calc_i_b_Ca = (_prvt_calc_g_b_Ca*(__OLD_[0]-_prvt_calc_E_Ca));
					_prvt_calc_i_NaCa = ((_prvt_Version==0.0000000000e+00))
?(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7400000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))))
:(((_prvt_calc_k_NaCa*((pow(_prvt_Na_i,3.0000000000e+00)*_prvt_Ca_o*exp((3.7430000000e-02*__OLD_[0]*_prvt_gamma_NaCa)))-(pow(_prvt_Na_o,3.0000000000e+00)*_prvt_Ca_i*exp((3.7430000000e-02*__OLD_[0]*(_prvt_gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(_prvt_d_NaCa*((_prvt_Ca_i*pow(_prvt_Na_o,3.0000000000e+00))+(_prvt_Ca_o*pow(_prvt_Na_i,3.0000000000e+00)))))));
					_prvt_calc_i_p = ((_prvt_calc_i_p_max*pow((_prvt_Na_i/(_prvt_K_m_Na+_prvt_Na_i)),3.0000000000e+00)*pow((_prvt_K_o/(_prvt_K_m_K+_prvt_K_o)),2.0000000000e+00)*1.6000000000e+00)/(1.5000000000e+00+exp(((-(__OLD_[0]+6.0000000000e+01))/4.0000000000e+01))));
					_prvt_calc_i_Ca_p = ((_prvt_calc_i_Ca_p_max*_prvt_Ca_i)/(_prvt_Ca_i+4.0000000000e-04));
					_prvt_calc_i_Na = (((((_prvt_calc_g_Na*pow(__OLD_[1],3.0000000000e+00)*_prvt_calc_h*_prvt_Na_o*pow(_prvt_F,2.0000000000e+00))/(_prvt_R*_prvt_T))*(exp((((__OLD_[0]-_prvt_calc_E_Na)*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))/(exp(((__OLD_[0]*_prvt_F)/(_prvt_R*_prvt_T)))-1.0000000000e+00))*__OLD_[0]);
					__K2_[0]= (((-1.0000000000e+00)/_prvt_calc_Cm)*(_prvt_calc_i_Na+_prvt_calc_i_Ca_L+_prvt_calc_i_Ca_T+_prvt_calc_i_to+_prvt_calc_i_sus+_prvt_calc_i_K_r+_prvt_calc_i_K_s+_prvt_calc_i_f_Na+_prvt_calc_i_f_K+_prvt_calc_i_b_Na+_prvt_calc_i_b_Ca+_prvt_calc_i_b_K+_prvt_calc_i_NaCa+_prvt_calc_i_p+_prvt_calc_i_Ca_p));
					_prvt_aux_tol = fabs(__OLD_[0])*_prvt_rel_tol_;
					__TOL_[0] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[0] = fabs((_prvt_dtime/2) * (__K1_[0] - __K2_[0])/__TOL_[0]);
					__K2_[5]= ((_prvt_calc_f_L_infinity-__OLD_[5])/_prvt_calc_tau_f_L);
					_prvt_aux_tol = fabs(__OLD_[5])*_prvt_rel_tol_;
					__TOL_[5] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[5] = fabs((_prvt_dtime/2) * (__K1_[5] - __K2_[5])/__TOL_[5]);
					__K2_[7]= ((_prvt_calc_f_T_infinity-__OLD_[7])/_prvt_calc_tau_f_T);
					_prvt_aux_tol = fabs(__OLD_[7])*_prvt_rel_tol_;
					__TOL_[7] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[7] = fabs((_prvt_dtime/2) * (__K1_[7] - __K2_[7])/__TOL_[7]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[1])
				{
					_prvt_calc_m_infinity = ((_prvt_Version==0.0000000000e+00))
?(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-__OLD_[0])/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)))
:(pow((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.0320000000e+01))/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
					_prvt_calc_tau_m = ((_prvt_Version==0.0000000000e+00))
?(((6.2470000000e-04/((8.3200000000e-01*exp(((-3.3500000000e-01)*(__OLD_[0]+5.6700000000e+01))))+(6.2700000000e-01*exp((8.2000000000e-02*(__OLD_[0]+6.5010000000e+01))))))+4.0000000000e-05))
:(((6.2470000000e-04/((8.3221660000e-01*exp(((-3.3566000000e-01)*(__OLD_[0]+5.6706200000e+01))))+(6.2740000000e-01*exp((8.2300000000e-02*(__OLD_[0]+6.5013100000e+01))))))+4.5690000000e-05));
					__K2_[1]= ((_prvt_calc_m_infinity-__OLD_[1])/_prvt_calc_tau_m);
					_prvt_aux_tol = fabs(__OLD_[1])*_prvt_rel_tol_;
					__TOL_[1] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[1] = fabs((_prvt_dtime/2) * (__K1_[1] - __K2_[1])/__TOL_[1]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[2])
				{
					_prvt_calc_h1_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+6.6100000000e+01)/6.4000000000e+00))));
					_prvt_calc_tau_h1 = (((3.7170000000e-06*exp(((-2.8150000000e-01)*(__OLD_[0]+1.7110000000e+01))))/(1.0000000000e+00+(3.7320000000e-03*exp(((-3.4260000000e-01)*(__OLD_[0]+3.7760000000e+01))))))+5.9770000000e-04);
					_prvt_calc_tau_h2 = (((3.1860000000e-08*exp(((-6.2190000000e-01)*(__OLD_[0]+1.8800000000e+01))))/(1.0000000000e+00+(7.1890000000e-05*exp(((-6.6830000000e-01)*(__OLD_[0]+3.4070000000e+01))))))+3.5560000000e-03);
					_prvt_calc_h2_infinity = _prvt_calc_h1_infinity;
					__K2_[2]= ((_prvt_calc_h1_infinity-__OLD_[2])/_prvt_calc_tau_h1);
					_prvt_aux_tol = fabs(__OLD_[2])*_prvt_rel_tol_;
					__TOL_[2] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[2] = fabs((_prvt_dtime/2) * (__K1_[2] - __K2_[2])/__TOL_[2]);
					__K2_[3]= ((_prvt_calc_h2_infinity-__OLD_[3])/_prvt_calc_tau_h2);
					_prvt_aux_tol = fabs(__OLD_[3])*_prvt_rel_tol_;
					__TOL_[3] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[3] = fabs((_prvt_dtime/2) * (__K1_[3] - __K2_[3])/__TOL_[3]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[3])
				{
					_prvt_calc_alpha_d_L = ((_prvt_Version==0.0000000000e+00))
?(((((-2.8380000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((_prvt_Version==1.0000000000e+00))
?(((((-2.8390000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00))))
:(((((-2.8400000000e+01)*(__OLD_[0]+3.5000000000e+01))/(exp(((-(__OLD_[0]+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*__OLD_[0])/(exp(((-2.0800000000e-01)*__OLD_[0]))-1.0000000000e+00)))))
;
					_prvt_calc_beta_d_L = ((_prvt_Version==1.0000000000e+00))
?(((1.1430000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)))
:(((1.1420000000e+01*(__OLD_[0]-5.0000000000e+00))/(exp((4.0000000000e-01*(__OLD_[0]-5.0000000000e+00)))-1.0000000000e+00)));
					_prvt_calc_d_L_infinity = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.3100000000e+01))/6.0000000000e+00)))))
:(((_prvt_Version==1.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2300000000e+01+(8.0000000000e-01*_prvt_calc_FCell)))/6.0000000000e+00)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+2.2200000000e+01))/6.0000000000e+00))))))
;
					_prvt_calc_tau_d_L = (2.0000000000e+00/(_prvt_calc_alpha_d_L+_prvt_calc_beta_d_L));
					__K2_[4]= ((_prvt_calc_d_L_infinity-__OLD_[4])/_prvt_calc_tau_d_L);
					_prvt_aux_tol = fabs(__OLD_[4])*_prvt_rel_tol_;
					__TOL_[4] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[4] = fabs((_prvt_dtime/2) * (__K1_[4] - __K2_[4])/__TOL_[4]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[4])
				{
					_prvt_calc_alpha_d_T = (1.0680000000e+03*exp(((__OLD_[0]+2.6300000000e+01)/3.0000000000e+01)));
					_prvt_calc_beta_d_T = (1.0680000000e+03*exp(((-(__OLD_[0]+2.6300000000e+01))/3.0000000000e+01)));
					_prvt_calc_d_T_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+3.7000000000e+01))/6.8000000000e+00))));
					_prvt_calc_tau_d_T = (1.0000000000e+00/(_prvt_calc_alpha_d_T+_prvt_calc_beta_d_T));
					__K2_[6]= ((_prvt_calc_d_T_infinity-__OLD_[6])/_prvt_calc_tau_d_T);
					_prvt_aux_tol = fabs(__OLD_[6])*_prvt_rel_tol_;
					__TOL_[6] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[6] = fabs((_prvt_dtime/2) * (__K1_[6] - __K2_[6])/__TOL_[6]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[5])
				{
					_prvt_calc_q_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+5.9370000000e+01)/1.3100000000e+01))));
					_prvt_calc_tau_q = ((_prvt_Version==0.0000000000e+00))
?((1.0100000000e-02+(6.5170000000e-02/(5.7000000000e-01*exp(((-8.0000000000e-02)*(__OLD_[0]+4.9000000000e+01)))))+(2.4000000000e-05*exp((1.0000000000e-01*(__OLD_[0]+5.0930000000e+01))))))
:(((_prvt_Version==1.0000000000e+00))
?(((1.0000000000e-03/3.0000000000e+00)*(3.0310000000e+01+(1.9550000000e+02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01+(1.0000000000e+01*_prvt_calc_FCell)))))+(7.1740000000e-01*exp(((2.7190000000e-01-(1.7190000000e-01*_prvt_calc_FCell))*1.0000000000e+00*(__OLD_[0]+4.0930000000e+01+(1.0000000000e+01*_prvt_calc_FCell))))))))))
:((1.0100000000e-02+(6.5170000000e-02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(__OLD_[0]+3.9000000000e+01))))+(7.1740000000e-01*exp((2.7190000000e-01*(__OLD_[0]+4.0930000000e+01)))))))))
;
					__K2_[8]= ((_prvt_calc_q_infinity-__OLD_[8])/_prvt_calc_tau_q);
					_prvt_aux_tol = fabs(__OLD_[8])*_prvt_rel_tol_;
					__TOL_[8] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[8] = fabs((_prvt_dtime/2) * (__K1_[8] - __K2_[8])/__TOL_[8]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[6])
				{
					_prvt_calc_r_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]-1.0930000000e+01))/1.9700000000e+01))));
					_prvt_calc_tau_r = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e-03*(2.9800000000e+00+(1.5590000000e+01/((1.0370000000e+00*exp((9.0000000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.2000000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:(((_prvt_Version==1.0000000000e+00))
?((2.5000000000e-03*(1.1910000000e+00+(7.8380000000e+00/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01)))))))))
:((1.0000000000e-03*(2.9800000000e+00+(1.9590000000e+01/((1.0370000000e+00*exp((9.0120000000e-02*(__OLD_[0]+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(__OLD_[0]+2.3840000000e+01))))))))))
;
					__K2_[9]= ((_prvt_calc_r_infinity-__OLD_[9])/_prvt_calc_tau_r);
					_prvt_aux_tol = fabs(__OLD_[9])*_prvt_rel_tol_;
					__TOL_[9] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[9] = fabs((_prvt_dtime/2) * (__K1_[9] - __K2_[9])/__TOL_[9]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[7])
				{
					_prvt_calc_P_af_infinity = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.4200000000e+01))/1.0600000000e+01)))))
:((1.0000000000e+00/(1.0000000000e+00+exp(((-(__OLD_[0]+1.3200000000e+01))/1.0600000000e+01)))));
					_prvt_calc_tau_P_af = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-9.0000000000e+00)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.2500000000e+01))))))
:((1.0000000000e+00/((3.7200000000e+01*exp(((__OLD_[0]-1.0000000000e+01)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.2500000000e+01))))));
					_prvt_calc_tau_P_as = ((_prvt_Version!=2.0000000000e+00))
?((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-9.0000000000e+00)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-9.0000000000e+00))/2.1600000000e+01))))))
:((1.0000000000e+00/((4.2000000000e+00*exp(((__OLD_[0]-1.0000000000e+01)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(__OLD_[0]-1.0000000000e+01))/2.1600000000e+01))))));
					_prvt_calc_P_as_infinity = _prvt_calc_P_af_infinity;
					__K2_[10]= ((_prvt_calc_P_af_infinity-__OLD_[10])/_prvt_calc_tau_P_af);
					_prvt_aux_tol = fabs(__OLD_[10])*_prvt_rel_tol_;
					__TOL_[10] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[10] = fabs((_prvt_dtime/2) * (__K1_[10] - __K2_[10])/__TOL_[10]);
					__K2_[11]= ((_prvt_calc_P_as_infinity-__OLD_[11])/_prvt_calc_tau_P_as);
					_prvt_aux_tol = fabs(__OLD_[11])*_prvt_rel_tol_;
					__TOL_[11] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[11] = fabs((_prvt_dtime/2) * (__K1_[11] - __K2_[11])/__TOL_[11]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[8])
				{
					_prvt_calc_tau_P_i = ((_prvt_Version==0.0000000000e+00))
?(2.0000000000e-03)
:(((_prvt_Version==1.0000000000e+00))
?(2.0000000000e-03)
:(6.0000000000e-03))
;
					_prvt_calc_P_i_infinity = (1.0000000000e+00/(1.0000000000e+00+exp(((__OLD_[0]+1.8600000000e+01)/1.0100000000e+01))));
					__K2_[12]= ((_prvt_calc_P_i_infinity-__OLD_[12])/_prvt_calc_tau_P_i);
					_prvt_aux_tol = fabs(__OLD_[12])*_prvt_rel_tol_;
					__TOL_[12] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[12] = fabs((_prvt_dtime/2) * (__K1_[12] - __K2_[12])/__TOL_[12]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[9])
				{
					_prvt_calc_alpha_xs = (1.4000000000e+01/(1.0000000000e+00+exp(((-(__OLD_[0]-4.0000000000e+01))/9.0000000000e+00))));
					_prvt_calc_beta_xs = (1.0000000000e+00*exp(((-__OLD_[0])/4.5000000000e+01)));
					__K2_[13]= ((_prvt_calc_alpha_xs*(1.0000000000e+00-__OLD_[13]))-(_prvt_calc_beta_xs*__OLD_[13]));
					_prvt_aux_tol = fabs(__OLD_[13])*_prvt_rel_tol_;
					__TOL_[13] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[13] = fabs((_prvt_dtime/2) * (__K1_[13] - __K2_[13])/__TOL_[13]);
				}
				if(omp_get_thread_num()==_prvt_tree_thread[10])
				{
					_prvt_calc_alpha_y = ((_prvt_Version==0.0000000000e+00))
?((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6620000000e+01))))
:((1.0000000000e+00*exp(((-(__OLD_[0]+7.8910000000e+01))/2.6630000000e+01))));
					_prvt_calc_beta_y = (1.0000000000e+00*exp(((__OLD_[0]+7.5130000000e+01)/2.1250000000e+01)));
					__K2_[14]= ((_prvt_calc_alpha_y*(1.0000000000e+00-__OLD_[14]))-(_prvt_calc_beta_y*__OLD_[14]));
					_prvt_aux_tol = fabs(__OLD_[14])*_prvt_rel_tol_;
					__TOL_[14] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;
					__ERROR_[14] = fabs((_prvt_dtime/2) * (__K1_[14] - __K2_[14])/__TOL_[14]);
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
						__NEW_[5] = __K1_[5] * _prvt_dtime + __OLD_AUX_[5];
						__NEW_[7] = __K1_[7] * _prvt_dtime + __OLD_AUX_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[1] = __K1_[1] * _prvt_dtime + __OLD_AUX_[1];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[2] = __K1_[2] * _prvt_dtime + __OLD_AUX_[2];
						__NEW_[3] = __K1_[3] * _prvt_dtime + __OLD_AUX_[3];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[4] = __K1_[4] * _prvt_dtime + __OLD_AUX_[4];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[6] = __K1_[6] * _prvt_dtime + __OLD_AUX_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[8] = __K1_[8] * _prvt_dtime + __OLD_AUX_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[9] = __K1_[9] * _prvt_dtime + __OLD_AUX_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[10] = __K1_[10] * _prvt_dtime + __OLD_AUX_[10];
						__NEW_[11] = __K1_[11] * _prvt_dtime + __OLD_AUX_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[12] = __K1_[12] * _prvt_dtime + __OLD_AUX_[12];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[13] = __K1_[13] * _prvt_dtime + __OLD_AUX_[13];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[14] = __K1_[14] * _prvt_dtime + __OLD_AUX_[14];
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
							this->m_old_ = __OLD_AUX_[1];
							this->m_new_ = __OLD_[1];
							this->h1_old_ = __OLD_AUX_[2];
							this->h1_new_ = __OLD_[2];
							this->h2_old_ = __OLD_AUX_[3];
							this->h2_new_ = __OLD_[3];
							this->d_L_old_ = __OLD_AUX_[4];
							this->d_L_new_ = __OLD_[4];
							this->f_L_old_ = __OLD_AUX_[5];
							this->f_L_new_ = __OLD_[5];
							this->d_T_old_ = __OLD_AUX_[6];
							this->d_T_new_ = __OLD_[6];
							this->f_T_old_ = __OLD_AUX_[7];
							this->f_T_new_ = __OLD_[7];
							this->q_old_ = __OLD_AUX_[8];
							this->q_new_ = __OLD_[8];
							this->r_old_ = __OLD_AUX_[9];
							this->r_new_ = __OLD_[9];
							this->P_af_old_ = __OLD_AUX_[10];
							this->P_af_new_ = __OLD_[10];
							this->P_as_old_ = __OLD_AUX_[11];
							this->P_as_new_ = __OLD_[11];
							this->P_i_old_ = __OLD_AUX_[12];
							this->P_i_new_ = __OLD_[12];
							this->xs_old_ = __OLD_AUX_[13];
							this->xs_new_ = __OLD_[13];
							this->y_old_ = __OLD_AUX_[14];
							this->y_new_ = __OLD_[14];
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
						__NEW_[5] = __K2_[5] * _prvt_dtime + __OLD_[5];
						__NEW_[7] = __K2_[7] * _prvt_dtime + __OLD_[7];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[1])
					{
						__NEW_[1] = __K2_[1] * _prvt_dtime + __OLD_[1];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[2])
					{
						__NEW_[2] = __K2_[2] * _prvt_dtime + __OLD_[2];
						__NEW_[3] = __K2_[3] * _prvt_dtime + __OLD_[3];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[3])
					{
						__NEW_[4] = __K2_[4] * _prvt_dtime + __OLD_[4];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[4])
					{
						__NEW_[6] = __K2_[6] * _prvt_dtime + __OLD_[6];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[5])
					{
						__NEW_[8] = __K2_[8] * _prvt_dtime + __OLD_[8];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[6])
					{
						__NEW_[9] = __K2_[9] * _prvt_dtime + __OLD_[9];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[7])
					{
						__NEW_[10] = __K2_[10] * _prvt_dtime + __OLD_[10];
						__NEW_[11] = __K2_[11] * _prvt_dtime + __OLD_[11];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[8])
					{
						__NEW_[12] = __K2_[12] * _prvt_dtime + __OLD_[12];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[9])
					{
						__NEW_[13] = __K2_[13] * _prvt_dtime + __OLD_[13];
					}
					if(omp_get_thread_num()==_prvt_tree_thread[10])
					{
						__NEW_[14] = __K2_[14] * _prvt_dtime + __OLD_[14];
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
		case 1:		return m_old_;    break;
		case 2:		return h1_old_;    break;
		case 3:		return h2_old_;    break;
		case 4:		return d_L_old_;    break;
		case 5:		return f_L_old_;    break;
		case 6:		return d_T_old_;    break;
		case 7:		return f_T_old_;    break;
		case 8:		return q_old_;    break;
		case 9:		return r_old_;    break;
		case 10:		return P_af_old_;    break;
		case 11:		return P_as_old_;    break;
		case 12:		return P_i_old_;    break;
		case 13:		return xs_old_;    break;
		case 14:		return y_old_;    break;
		default:	return 1;    break;
		}
	}

	double Solveode::getLadoDireito(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return V_lado_direito_;    break;
		case 1:		return m_lado_direito_;    break;
		case 2:		return h1_lado_direito_;    break;
		case 3:		return h2_lado_direito_;    break;
		case 4:		return d_L_lado_direito_;    break;
		case 5:		return f_L_lado_direito_;    break;
		case 6:		return d_T_lado_direito_;    break;
		case 7:		return f_T_lado_direito_;    break;
		case 8:		return q_lado_direito_;    break;
		case 9:		return r_lado_direito_;    break;
		case 10:		return P_af_lado_direito_;    break;
		case 11:		return P_as_lado_direito_;    break;
		case 12:		return P_i_lado_direito_;    break;
		case 13:		return xs_lado_direito_;    break;
		case 14:		return y_lado_direito_;    break;
		default:	return 1;    break;
		}
	}

	double Solveode::getParameters(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return time;    break;
		case 1:		return dCell;    break;
		case 2:		return Version;    break;
		case 3:		return FCellConstant;    break;
		case 4:		return CmCentre;    break;
		case 5:		return CmPeriphery;    break;
		case 6:		return g_Na_Centre_Published;    break;
		case 7:		return g_Na_Periphery_Published;    break;
		case 8:		return g_Na_Centre_0DCapable;    break;
		case 9:		return g_Na_Periphery_0DCapable;    break;
		case 10:		return g_Na_Centre_1DCapable;    break;
		case 11:		return g_Na_Periphery_1DCapable;    break;
		case 12:		return Na_o;    break;
		case 13:		return F;    break;
		case 14:		return R;    break;
		case 15:		return T;    break;
		case 16:		return g_Ca_L_Centre_Published;    break;
		case 17:		return g_Ca_L_Periphery_Published;    break;
		case 18:		return g_Ca_L_Centre_0DCapable;    break;
		case 19:		return g_Ca_L_Periphery_0DCapable;    break;
		case 20:		return g_Ca_L_Centre_1DCapable;    break;
		case 21:		return g_Ca_L_Periphery_1DCapable;    break;
		case 22:		return E_Ca_L;    break;
		case 23:		return g_Ca_T_Centre_Published;    break;
		case 24:		return g_Ca_T_Periphery_Published;    break;
		case 25:		return g_Ca_T_Centre_0DCapable;    break;
		case 26:		return g_Ca_T_Periphery_0DCapable;    break;
		case 27:		return g_Ca_T_Centre_1DCapable;    break;
		case 28:		return g_Ca_T_Periphery_1DCapable;    break;
		case 29:		return E_Ca_T;    break;
		case 30:		return g_to_Centre_Published;    break;
		case 31:		return g_to_Periphery_Published;    break;
		case 32:		return g_to_Centre_0DCapable;    break;
		case 33:		return g_to_Periphery_0DCapable;    break;
		case 34:		return g_to_Centre_1DCapable;    break;
		case 35:		return g_to_Periphery_1DCapable;    break;
		case 36:		return g_sus_Centre_Published;    break;
		case 37:		return g_sus_Periphery_Published;    break;
		case 38:		return g_sus_Centre_0DCapable;    break;
		case 39:		return g_sus_Periphery_0DCapable;    break;
		case 40:		return g_sus_Centre_1DCapable;    break;
		case 41:		return g_sus_Periphery_1DCapable;    break;
		case 42:		return g_K_r_Centre_Published;    break;
		case 43:		return g_K_r_Periphery_Published;    break;
		case 44:		return g_K_r_Centre_0DCapable;    break;
		case 45:		return g_K_r_Periphery_0DCapable;    break;
		case 46:		return g_K_r_Centre_1DCapable;    break;
		case 47:		return g_K_r_Periphery_1DCapable;    break;
		case 48:		return g_K_s_Centre_Published;    break;
		case 49:		return g_K_s_Periphery_Published;    break;
		case 50:		return g_K_s_Centre_0DCapable;    break;
		case 51:		return g_K_s_Periphery_0DCapable;    break;
		case 52:		return g_K_s_Centre_1DCapable;    break;
		case 53:		return g_K_s_Periphery_1DCapable;    break;
		case 54:		return g_f_Na_Centre_Published;    break;
		case 55:		return g_f_Na_Periphery_Published;    break;
		case 56:		return g_f_Na_Centre_0DCapable;    break;
		case 57:		return g_f_Na_Periphery_0DCapable;    break;
		case 58:		return g_f_Na_Centre_1DCapable;    break;
		case 59:		return g_f_Na_Periphery_1DCapable;    break;
		case 60:		return g_f_K_Centre_Published;    break;
		case 61:		return g_f_K_Periphery_Published;    break;
		case 62:		return g_f_K_Centre_0DCapable;    break;
		case 63:		return g_f_K_Periphery_0DCapable;    break;
		case 64:		return g_f_K_Centre_1DCapable;    break;
		case 65:		return g_f_K_Periphery_1DCapable;    break;
		case 66:		return g_b_Na_Centre_Published;    break;
		case 67:		return g_b_Na_Periphery_Published;    break;
		case 68:		return g_b_Na_Centre_0DCapable;    break;
		case 69:		return g_b_Na_Periphery_0DCapable;    break;
		case 70:		return g_b_Na_Centre_1DCapable;    break;
		case 71:		return g_b_Na_Periphery_1DCapable;    break;
		case 72:		return g_b_K_Centre_Published;    break;
		case 73:		return g_b_K_Periphery_Published;    break;
		case 74:		return g_b_K_Centre_0DCapable;    break;
		case 75:		return g_b_K_Periphery_0DCapable;    break;
		case 76:		return g_b_K_Centre_1DCapable;    break;
		case 77:		return g_b_K_Periphery_1DCapable;    break;
		case 78:		return g_b_Ca_Centre_Published;    break;
		case 79:		return g_b_Ca_Periphery_Published;    break;
		case 80:		return g_b_Ca_Centre_0DCapable;    break;
		case 81:		return g_b_Ca_Periphery_0DCapable;    break;
		case 82:		return g_b_Ca_Centre_1DCapable;    break;
		case 83:		return g_b_Ca_Periphery_1DCapable;    break;
		case 84:		return k_NaCa_Centre_Published;    break;
		case 85:		return k_NaCa_Periphery_Published;    break;
		case 86:		return k_NaCa_Centre_0DCapable;    break;
		case 87:		return k_NaCa_Periphery_0DCapable;    break;
		case 88:		return k_NaCa_Centre_1DCapable;    break;
		case 89:		return k_NaCa_Periphery_1DCapable;    break;
		case 90:		return Na_i;    break;
		case 91:		return Ca_o;    break;
		case 92:		return gamma_NaCa;    break;
		case 93:		return Ca_i;    break;
		case 94:		return d_NaCa;    break;
		case 95:		return i_p_max_Centre_Published;    break;
		case 96:		return i_p_max_Periphery_Published;    break;
		case 97:		return i_p_max_Centre_0DCapable;    break;
		case 98:		return i_p_max_Periphery_0DCapable;    break;
		case 99:		return i_p_max_Centre_1DCapable;    break;
		case 100:		return i_p_max_Periphery_1DCapable;    break;
		case 101:		return K_m_Na;    break;
		case 102:		return K_o;    break;
		case 103:		return K_m_K;    break;
		case 104:		return i_Ca_p_max_Centre_Published;    break;
		case 105:		return i_Ca_p_max_Periphery_Published;    break;
		case 106:		return i_Ca_p_max_Centre_0DCapable;    break;
		case 107:		return i_Ca_p_max_Periphery_0DCapable;    break;
		case 108:		return i_Ca_p_max_Centre_1DCapable;    break;
		case 109:		return i_Ca_p_max_Periphery_1DCapable;    break;
		case 110:		return K_i;    break;
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
		Variables v("|V#|m#|h1#|h2#|d_L#|f_L#|d_T#|f_T#|q#|r#|P_af#|P_as#|P_i#|xs#|y#");
		return v;
	}
	Variables Solveode::get_Parameters()
	{
		Variables v("|time#|dCell#|Version#|FCellConstant#|CmCentre#|CmPeriphery#|g_Na_Centre_Published#|g_Na_Periphery_Published#|g_Na_Centre_0DCapable#|g_Na_Periphery_0DCapable#|g_Na_Centre_1DCapable#|g_Na_Periphery_1DCapable#|Na_o#|F#|R#|T#|g_Ca_L_Centre_Published#|g_Ca_L_Periphery_Published#|g_Ca_L_Centre_0DCapable#|g_Ca_L_Periphery_0DCapable#|g_Ca_L_Centre_1DCapable#|g_Ca_L_Periphery_1DCapable#|E_Ca_L#|g_Ca_T_Centre_Published#|g_Ca_T_Periphery_Published#|g_Ca_T_Centre_0DCapable#|g_Ca_T_Periphery_0DCapable#|g_Ca_T_Centre_1DCapable#|g_Ca_T_Periphery_1DCapable#|E_Ca_T#|g_to_Centre_Published#|g_to_Periphery_Published#|g_to_Centre_0DCapable#|g_to_Periphery_0DCapable#|g_to_Centre_1DCapable#|g_to_Periphery_1DCapable#|g_sus_Centre_Published#|g_sus_Periphery_Published#|g_sus_Centre_0DCapable#|g_sus_Periphery_0DCapable#|g_sus_Centre_1DCapable#|g_sus_Periphery_1DCapable#|g_K_r_Centre_Published#|g_K_r_Periphery_Published#|g_K_r_Centre_0DCapable#|g_K_r_Periphery_0DCapable#|g_K_r_Centre_1DCapable#|g_K_r_Periphery_1DCapable#|g_K_s_Centre_Published#|g_K_s_Periphery_Published#|g_K_s_Centre_0DCapable#|g_K_s_Periphery_0DCapable#|g_K_s_Centre_1DCapable#|g_K_s_Periphery_1DCapable#|g_f_Na_Centre_Published#|g_f_Na_Periphery_Published#|g_f_Na_Centre_0DCapable#|g_f_Na_Periphery_0DCapable#|g_f_Na_Centre_1DCapable#|g_f_Na_Periphery_1DCapable#|g_f_K_Centre_Published#|g_f_K_Periphery_Published#|g_f_K_Centre_0DCapable#|g_f_K_Periphery_0DCapable#|g_f_K_Centre_1DCapable#|g_f_K_Periphery_1DCapable#|g_b_Na_Centre_Published#|g_b_Na_Periphery_Published#|g_b_Na_Centre_0DCapable#|g_b_Na_Periphery_0DCapable#|g_b_Na_Centre_1DCapable#|g_b_Na_Periphery_1DCapable#|g_b_K_Centre_Published#|g_b_K_Periphery_Published#|g_b_K_Centre_0DCapable#|g_b_K_Periphery_0DCapable#|g_b_K_Centre_1DCapable#|g_b_K_Periphery_1DCapable#|g_b_Ca_Centre_Published#|g_b_Ca_Periphery_Published#|g_b_Ca_Centre_0DCapable#|g_b_Ca_Periphery_0DCapable#|g_b_Ca_Centre_1DCapable#|g_b_Ca_Periphery_1DCapable#|k_NaCa_Centre_Published#|k_NaCa_Periphery_Published#|k_NaCa_Centre_0DCapable#|k_NaCa_Periphery_0DCapable#|k_NaCa_Centre_1DCapable#|k_NaCa_Periphery_1DCapable#|Na_i#|Ca_o#|gamma_NaCa#|Ca_i#|d_NaCa#|i_p_max_Centre_Published#|i_p_max_Periphery_Published#|i_p_max_Centre_0DCapable#|i_p_max_Periphery_0DCapable#|i_p_max_Centre_1DCapable#|i_p_max_Periphery_1DCapable#|K_m_Na#|K_o#|K_m_K#|i_Ca_p_max_Centre_Published#|i_Ca_p_max_Periphery_Published#|i_Ca_p_max_Centre_0DCapable#|i_Ca_p_max_Periphery_0DCapable#|i_Ca_p_max_Centre_1DCapable#|i_Ca_p_max_Periphery_1DCapable#|K_i#");
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
		m_old_ = m_ini_;
		if(m != NULL)free( m);
			m = (double *)malloc(sizeof(double)*num_results__);
		h1_old_ = h1_ini_;
		if(h1 != NULL)free( h1);
			h1 = (double *)malloc(sizeof(double)*num_results__);
		h2_old_ = h2_ini_;
		if(h2 != NULL)free( h2);
			h2 = (double *)malloc(sizeof(double)*num_results__);
		d_L_old_ = d_L_ini_;
		if(d_L != NULL)free( d_L);
			d_L = (double *)malloc(sizeof(double)*num_results__);
		f_L_old_ = f_L_ini_;
		if(f_L != NULL)free( f_L);
			f_L = (double *)malloc(sizeof(double)*num_results__);
		d_T_old_ = d_T_ini_;
		if(d_T != NULL)free( d_T);
			d_T = (double *)malloc(sizeof(double)*num_results__);
		f_T_old_ = f_T_ini_;
		if(f_T != NULL)free( f_T);
			f_T = (double *)malloc(sizeof(double)*num_results__);
		q_old_ = q_ini_;
		if(q != NULL)free( q);
			q = (double *)malloc(sizeof(double)*num_results__);
		r_old_ = r_ini_;
		if(r != NULL)free( r);
			r = (double *)malloc(sizeof(double)*num_results__);
		P_af_old_ = P_af_ini_;
		if(P_af != NULL)free( P_af);
			P_af = (double *)malloc(sizeof(double)*num_results__);
		P_as_old_ = P_as_ini_;
		if(P_as != NULL)free( P_as);
			P_as = (double *)malloc(sizeof(double)*num_results__);
		P_i_old_ = P_i_ini_;
		if(P_i != NULL)free( P_i);
			P_i = (double *)malloc(sizeof(double)*num_results__);
		xs_old_ = xs_ini_;
		if(xs != NULL)free( xs);
			xs = (double *)malloc(sizeof(double)*num_results__);
		y_old_ = y_ini_;
		if(y != NULL)free( y);
			y = (double *)malloc(sizeof(double)*num_results__);
		this->timeSaving = dtime;

		double diff=0;
		int counter=0;
		for (int i = 0; i< iterations;i++ )
		{
			this->time_new += dtime;

			rightHandSideFunction.function(this);
			this->V_new_ = this->V_old_ + this->V_lado_direito_ * this->dtime;
			this->m_new_ = this->m_old_ + this->m_lado_direito_ * this->dtime;
			this->h1_new_ = this->h1_old_ + this->h1_lado_direito_ * this->dtime;
			this->h2_new_ = this->h2_old_ + this->h2_lado_direito_ * this->dtime;
			this->d_L_new_ = this->d_L_old_ + this->d_L_lado_direito_ * this->dtime;
			this->f_L_new_ = this->f_L_old_ + this->f_L_lado_direito_ * this->dtime;
			this->d_T_new_ = this->d_T_old_ + this->d_T_lado_direito_ * this->dtime;
			this->f_T_new_ = this->f_T_old_ + this->f_T_lado_direito_ * this->dtime;
			this->q_new_ = this->q_old_ + this->q_lado_direito_ * this->dtime;
			this->r_new_ = this->r_old_ + this->r_lado_direito_ * this->dtime;
			this->P_af_new_ = this->P_af_old_ + this->P_af_lado_direito_ * this->dtime;
			this->P_as_new_ = this->P_as_old_ + this->P_as_lado_direito_ * this->dtime;
			this->P_i_new_ = this->P_i_old_ + this->P_i_lado_direito_ * this->dtime;
			this->xs_new_ = this->xs_old_ + this->xs_lado_direito_ * this->dtime;
			this->y_new_ = this->y_old_ + this->y_lado_direito_ * this->dtime;
			diff =  _agos_round(this->time_new - timeSaving, 10);
			if(diff==0){
				this->timeSaving += svRate;
				time_vec__[counter] = this->time_new;
				V[counter] = this->V_new_;
				m[counter] = this->m_new_;
				h1[counter] = this->h1_new_;
				h2[counter] = this->h2_new_;
				d_L[counter] = this->d_L_new_;
				f_L[counter] = this->f_L_new_;
				d_T[counter] = this->d_T_new_;
				f_T[counter] = this->f_T_new_;
				q[counter] = this->q_new_;
				r[counter] = this->r_new_;
				P_af[counter] = this->P_af_new_;
				P_as[counter] = this->P_as_new_;
				P_i[counter] = this->P_i_new_;
				xs[counter] = this->xs_new_;
				y[counter] = this->y_new_;
				counter++;
			}
			this->V_old_ = this->V_new_;
			this->m_old_ = this->m_new_;
			this->h1_old_ = this->h1_new_;
			this->h2_old_ = this->h2_new_;
			this->d_L_old_ = this->d_L_new_;
			this->f_L_old_ = this->f_L_new_;
			this->d_T_old_ = this->d_T_new_;
			this->f_T_old_ = this->f_T_new_;
			this->q_old_ = this->q_new_;
			this->r_old_ = this->r_new_;
			this->P_af_old_ = this->P_af_new_;
			this->P_as_old_ = this->P_as_new_;
			this->P_i_old_ = this->P_i_new_;
			this->xs_old_ = this->xs_new_;
			this->y_old_ = this->y_new_;
		}
		double h_jac_[numEDO];
		double quociente = 1000.0;
		h_jac_[0] = fabs(_agos_min(V, num_results__) / _agos_max(V, num_results__) );
		h_jac_[1] = fabs(_agos_min(m, num_results__) / _agos_max(m, num_results__) );
		h_jac_[2] = fabs(_agos_min(h1, num_results__) / _agos_max(h1, num_results__) );
		h_jac_[3] = fabs(_agos_min(h2, num_results__) / _agos_max(h2, num_results__) );
		h_jac_[4] = fabs(_agos_min(d_L, num_results__) / _agos_max(d_L, num_results__) );
		h_jac_[5] = fabs(_agos_min(f_L, num_results__) / _agos_max(f_L, num_results__) );
		h_jac_[6] = fabs(_agos_min(d_T, num_results__) / _agos_max(d_T, num_results__) );
		h_jac_[7] = fabs(_agos_min(f_T, num_results__) / _agos_max(f_T, num_results__) );
		h_jac_[8] = fabs(_agos_min(q, num_results__) / _agos_max(q, num_results__) );
		h_jac_[9] = fabs(_agos_min(r, num_results__) / _agos_max(r, num_results__) );
		h_jac_[10] = fabs(_agos_min(P_af, num_results__) / _agos_max(P_af, num_results__) );
		h_jac_[11] = fabs(_agos_min(P_as, num_results__) / _agos_max(P_as, num_results__) );
		h_jac_[12] = fabs(_agos_min(P_i, num_results__) / _agos_max(P_i, num_results__) );
		h_jac_[13] = fabs(_agos_min(xs, num_results__) / _agos_max(xs, num_results__) );
		h_jac_[14] = fabs(_agos_min(y, num_results__) / _agos_max(y, num_results__) );
		for(int l=0;l<numEDO;l++){
			h_jac_[l] = (h_jac_[l]==0 || h_jac_[l]==AGOS_NAN || h_jac_[l]==AGOS_INF)?this->dtime:h_jac_[l];
		}
		this->timeSaving = this->dtime;

		this->time_new = this->time;

		V_old_ = V_ini_;
		m_old_ = m_ini_;
		h1_old_ = h1_ini_;
		h2_old_ = h2_ini_;
		d_L_old_ = d_L_ini_;
		f_L_old_ = f_L_ini_;
		d_T_old_ = d_T_ini_;
		f_T_old_ = f_T_ini_;
		q_old_ = q_ini_;
		r_old_ = r_ini_;
		P_af_old_ = P_af_ini_;
		P_as_old_ = P_as_ini_;
		P_i_old_ = P_i_ini_;
		xs_old_ = xs_ini_;
		y_old_ = y_ini_;
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
			this->m_new_ = this->dtime*(this->m_lado_direito_) + this->m_old_;
			this->h1_new_ = this->dtime*(this->h1_lado_direito_) + this->h1_old_;
			this->h2_new_ = this->dtime*(this->h2_lado_direito_) + this->h2_old_;
			this->d_L_new_ = this->dtime*(this->d_L_lado_direito_) + this->d_L_old_;
			this->f_L_new_ = this->dtime*(this->f_L_lado_direito_) + this->f_L_old_;
			this->d_T_new_ = this->dtime*(this->d_T_lado_direito_) + this->d_T_old_;
			this->f_T_new_ = this->dtime*(this->f_T_lado_direito_) + this->f_T_old_;
			this->q_new_ = this->dtime*(this->q_lado_direito_) + this->q_old_;
			this->r_new_ = this->dtime*(this->r_lado_direito_) + this->r_old_;
			this->P_af_new_ = this->dtime*(this->P_af_lado_direito_) + this->P_af_old_;
			this->P_as_new_ = this->dtime*(this->P_as_lado_direito_) + this->P_as_old_;
			this->P_i_new_ = this->dtime*(this->P_i_lado_direito_) + this->P_i_old_;
			this->xs_new_ = this->dtime*(this->xs_lado_direito_) + this->xs_old_;
			this->y_new_ = this->dtime*(this->y_lado_direito_) + this->y_old_;
			diff =  _agos_round(this->time_new - timeSaving, 10);
			if(diff==0){
				this->timeSaving += svRate;
				V[counter2] = V_new_;
				m[counter2] = m_new_;
				h1[counter2] = h1_new_;
				h2[counter2] = h2_new_;
				d_L[counter2] = d_L_new_;
				f_L[counter2] = f_L_new_;
				d_T[counter2] = d_T_new_;
				f_T[counter2] = f_T_new_;
				q[counter2] = q_new_;
				r[counter2] = r_new_;
				P_af[counter2] = P_af_new_;
				P_as[counter2] = P_as_new_;
				P_i[counter2] = P_i_new_;
				xs[counter2] = xs_new_;
				y[counter2] = y_new_;
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
			m_old_ = m_new_;
			h1_old_ = h1_new_;
			h2_old_ = h2_new_;
			d_L_old_ = d_L_new_;
			f_L_old_ = f_L_new_;
			d_T_old_ = d_T_new_;
			f_T_old_ = f_T_new_;
			q_old_ = q_new_;
			r_old_ = r_new_;
			P_af_old_ = P_af_new_;
			P_as_old_ = P_as_new_;
			P_i_old_ = P_i_new_;
			xs_old_ = xs_new_;
			y_old_ = y_new_;
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
		case 1:		return m;    break;
		case 2:		return h1;    break;
		case 3:		return h2;    break;
		case 4:		return d_L;    break;
		case 5:		return f_L;    break;
		case 6:		return d_T;    break;
		case 7:		return f_T;    break;
		case 8:		return q;    break;
		case 9:		return r;    break;
		case 10:		return P_af;    break;
		case 11:		return P_as;    break;
		case 12:		return P_i;    break;
		case 13:		return xs;    break;
		case 14:		return y;    break;
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
		if((Version==0.0000000000e+00)){
			return (((1.0700000000e+00*((3.0000000000e+00*dCell)-1.0000000000e-01))/(3.0000000000e+00*(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*dCell)-2.0500000000e+00))/2.9500000000e-01)))))));
		}else if((Version==1.0000000000e+00)){
			return (((FCellConstant*dCell)/(1.0000000000e+00+(7.7450000000e-01*exp(((-((3.0000000000e+00*dCell)-2.0500000000e+00))/2.9500000000e-01))))));
		} else{
			return (((1.0700000000e+00*2.9000000000e+01*dCell)/(3.0000000000e+01*(1.0000000000e+00+(7.7450000000e-01*exp(((-((2.9000000000e+01*dCell)-2.4500000000e+01))/1.9500000000e+00)))))));
		}
	}
	double Solveode::ifnumber_1(){
		if((Version==0.0000000000e+00)){
			return ((g_Na_Centre_Published+(calc_FCell*(g_Na_Periphery_Published-g_Na_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_Na_Centre_0DCapable+(calc_FCell*(g_Na_Periphery_0DCapable-g_Na_Centre_0DCapable))));
		} else{
			return ((g_Na_Centre_1DCapable+(calc_FCell*(g_Na_Periphery_1DCapable-g_Na_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_2(){
		if((Version==0.0000000000e+00)){
			return (pow((1.0000000000e+00/(1.0000000000e+00+exp(((-V_old_)/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
		}else{
			return (pow((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+3.0320000000e+01))/5.4600000000e+00)))),(1.0000000000e+00/3.0000000000e+00)));
		}
	}
	double Solveode::ifnumber_3(){
		if((Version==0.0000000000e+00)){
			return (((6.2470000000e-04/((8.3200000000e-01*exp(((-3.3500000000e-01)*(V_old_+5.6700000000e+01))))+(6.2700000000e-01*exp((8.2000000000e-02*(V_old_+6.5010000000e+01))))))+4.0000000000e-05));
		}else{
			return (((6.2470000000e-04/((8.3221660000e-01*exp(((-3.3566000000e-01)*(V_old_+5.6706200000e+01))))+(6.2740000000e-01*exp((8.2300000000e-02*(V_old_+6.5013100000e+01))))))+4.5690000000e-05));
		}
	}
	double Solveode::ifnumber_4(){
		if((Version==0.0000000000e+00)){
			return ((((9.5200000000e-02*exp(((-6.3000000000e-02)*(V_old_+3.4400000000e+01))))/(1.0000000000e+00+(1.6600000000e+00*exp(((-2.2500000000e-01)*(V_old_+6.3700000000e+01))))))+8.6900000000e-02));
		}else{
			return ((((9.5180000000e-02*exp(((-6.3060000000e-02)*(V_old_+3.4400000000e+01))))/(1.0000000000e+00+(1.6620000000e+00*exp(((-2.2510000000e-01)*(V_old_+6.3700000000e+01))))))+8.6930000000e-02));
		}
	}
	double Solveode::ifnumber_5(){
		if((Version==0.0000000000e+00)){
			return ((g_Ca_L_Centre_Published+(calc_FCell*(g_Ca_L_Periphery_Published-g_Ca_L_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_Ca_L_Centre_0DCapable+(calc_FCell*(g_Ca_L_Periphery_0DCapable-g_Ca_L_Centre_0DCapable))));
		} else{
			return ((g_Ca_L_Centre_1DCapable+(calc_FCell*(g_Ca_L_Periphery_1DCapable-g_Ca_L_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_6(){
		if((Version==0.0000000000e+00)){
			return (((((-2.8380000000e+01)*(V_old_+3.5000000000e+01))/(exp(((-(V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*V_old_)/(exp(((-2.0800000000e-01)*V_old_))-1.0000000000e+00))));
		}else if((Version==1.0000000000e+00)){
			return (((((-2.8390000000e+01)*(V_old_+3.5000000000e+01))/(exp(((-(V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*V_old_)/(exp(((-2.0800000000e-01)*V_old_))-1.0000000000e+00))));
		} else{
			return (((((-2.8400000000e+01)*(V_old_+3.5000000000e+01))/(exp(((-(V_old_+3.5000000000e+01))/2.5000000000e+00))-1.0000000000e+00))-((8.4900000000e+01*V_old_)/(exp(((-2.0800000000e-01)*V_old_))-1.0000000000e+00))));
		}
	}
	double Solveode::ifnumber_7(){
		if((Version==1.0000000000e+00)){
			return (((1.1430000000e+01*(V_old_-5.0000000000e+00))/(exp((4.0000000000e-01*(V_old_-5.0000000000e+00)))-1.0000000000e+00)));
		}else{
			return (((1.1420000000e+01*(V_old_-5.0000000000e+00))/(exp((4.0000000000e-01*(V_old_-5.0000000000e+00)))-1.0000000000e+00)));
		}
	}
	double Solveode::ifnumber_8(){
		if((Version==0.0000000000e+00)){
			return ((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+2.3100000000e+01))/6.0000000000e+00)))));
		}else if((Version==1.0000000000e+00)){
			return ((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+2.2300000000e+01+(8.0000000000e-01*calc_FCell)))/6.0000000000e+00)))));
		} else{
			return ((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+2.2200000000e+01))/6.0000000000e+00)))));
		}
	}
	double Solveode::ifnumber_9(){
		if((Version==1.0000000000e+00)){
			return (((3.7500000000e+00*(V_old_+2.8000000000e+01))/(exp(((V_old_+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
		}else{
			return (((3.1200000000e+00*(V_old_+2.8000000000e+01))/(exp(((V_old_+2.8000000000e+01)/4.0000000000e+00))-1.0000000000e+00)));
		}
	}
	double Solveode::ifnumber_10(){
		if((Version==1.0000000000e+00)){
			return ((3.0000000000e+01/(1.0000000000e+00+exp(((-(V_old_+2.8000000000e+01))/4.0000000000e+00)))));
		}else{
			return ((2.5000000000e+01/(1.0000000000e+00+exp(((-(V_old_+2.8000000000e+01))/4.0000000000e+00)))));
		}
	}
	double Solveode::ifnumber_11(){
		if((Version==1.0000000000e+00)){
			return (((1.2000000000e+00-(2.0000000000e-01*calc_FCell))/(calc_alpha_f_L+calc_beta_f_L)));
		}else{
			return ((1.0000000000e+00/(calc_alpha_f_L+calc_beta_f_L)));
		}
	}
	double Solveode::ifnumber_12(){
		if((Version==0.0000000000e+00)){
			return ((g_Ca_T_Centre_Published+(calc_FCell*(g_Ca_T_Periphery_Published-g_Ca_T_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_Ca_T_Centre_0DCapable+(calc_FCell*(g_Ca_T_Periphery_0DCapable-g_Ca_T_Centre_0DCapable))));
		} else{
			return ((g_Ca_T_Centre_1DCapable+(calc_FCell*(g_Ca_T_Periphery_1DCapable-g_Ca_T_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_13(){
		if((Version==1.0000000000e+00)){
			return ((1.5300000000e+01*exp(((-(V_old_+7.1000000000e+01+(7.0000000000e-01*calc_FCell)))/8.3300000000e+01))));
		}else{
			return ((1.5300000000e+01*exp(((-(V_old_+7.1700000000e+01))/8.3300000000e+01))));
		}
	}
	double Solveode::ifnumber_14(){
		if((Version==1.0000000000e+00)){
			return ((1.5000000000e+01*exp(((V_old_+7.1000000000e+01)/1.5380000000e+01))));
		}else{
			return ((1.5000000000e+01*exp(((V_old_+7.1700000000e+01)/1.5380000000e+01))));
		}
	}
	double Solveode::ifnumber_15(){
		if((Version==0.0000000000e+00)){
			return ((g_to_Centre_Published+(calc_FCell*(g_to_Periphery_Published-g_to_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_to_Centre_0DCapable+(calc_FCell*(g_to_Periphery_0DCapable-g_to_Centre_0DCapable))));
		} else{
			return ((g_to_Centre_1DCapable+(calc_FCell*(g_to_Periphery_1DCapable-g_to_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_16(){
		if((Version==0.0000000000e+00)){
			return ((g_sus_Centre_Published+(calc_FCell*(g_sus_Periphery_Published-g_sus_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_sus_Centre_0DCapable+(calc_FCell*(g_sus_Periphery_0DCapable-g_sus_Centre_0DCapable))));
		} else{
			return ((g_sus_Centre_1DCapable+(calc_FCell*(g_sus_Periphery_1DCapable-g_sus_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_17(){
		if((Version==0.0000000000e+00)){
			return ((1.0100000000e-02+(6.5170000000e-02/(5.7000000000e-01*exp(((-8.0000000000e-02)*(V_old_+4.9000000000e+01)))))+(2.4000000000e-05*exp((1.0000000000e-01*(V_old_+5.0930000000e+01))))));
		}else if((Version==1.0000000000e+00)){
			return (((1.0000000000e-03/3.0000000000e+00)*(3.0310000000e+01+(1.9550000000e+02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(V_old_+3.9000000000e+01+(1.0000000000e+01*calc_FCell)))))+(7.1740000000e-01*exp(((2.7190000000e-01-(1.7190000000e-01*calc_FCell))*1.0000000000e+00*(V_old_+4.0930000000e+01+(1.0000000000e+01*calc_FCell))))))))));
		} else{
			return ((1.0100000000e-02+(6.5170000000e-02/((5.6860000000e-01*exp(((-8.1610000000e-02)*(V_old_+3.9000000000e+01))))+(7.1740000000e-01*exp((2.7190000000e-01*(V_old_+4.0930000000e+01))))))));
		}
	}
	double Solveode::ifnumber_18(){
		if((Version==0.0000000000e+00)){
			return ((1.0000000000e-03*(2.9800000000e+00+(1.5590000000e+01/((1.0370000000e+00*exp((9.0000000000e-02*(V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.2000000000e-01)*(V_old_+2.3840000000e+01)))))))));
		}else if((Version==1.0000000000e+00)){
			return ((2.5000000000e-03*(1.1910000000e+00+(7.8380000000e+00/((1.0370000000e+00*exp((9.0120000000e-02*(V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(V_old_+2.3840000000e+01)))))))));
		} else{
			return ((1.0000000000e-03*(2.9800000000e+00+(1.9590000000e+01/((1.0370000000e+00*exp((9.0120000000e-02*(V_old_+3.0610000000e+01))))+(3.6900000000e-01*exp(((-1.1900000000e-01)*(V_old_+2.3840000000e+01)))))))));
		}
	}
	double Solveode::ifnumber_19(){
		if((Version==0.0000000000e+00)){
			return ((g_K_r_Centre_Published+(calc_FCell*(g_K_r_Periphery_Published-g_K_r_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_K_r_Centre_0DCapable+(calc_FCell*(g_K_r_Periphery_0DCapable-g_K_r_Centre_0DCapable))));
		} else{
			return ((g_K_r_Centre_1DCapable+(calc_FCell*(g_K_r_Periphery_1DCapable-g_K_r_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_20(){
		if((Version!=2.0000000000e+00)){
			return ((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+1.4200000000e+01))/1.0600000000e+01)))));
		}else{
			return ((1.0000000000e+00/(1.0000000000e+00+exp(((-(V_old_+1.3200000000e+01))/1.0600000000e+01)))));
		}
	}
	double Solveode::ifnumber_21(){
		if((Version!=2.0000000000e+00)){
			return ((1.0000000000e+00/((3.7200000000e+01*exp(((V_old_-9.0000000000e+00)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(V_old_-9.0000000000e+00))/2.2500000000e+01))))));
		}else{
			return ((1.0000000000e+00/((3.7200000000e+01*exp(((V_old_-1.0000000000e+01)/1.5900000000e+01)))+(9.6000000000e-01*exp(((-(V_old_-1.0000000000e+01))/2.2500000000e+01))))));
		}
	}
	double Solveode::ifnumber_22(){
		if((Version!=2.0000000000e+00)){
			return ((1.0000000000e+00/((4.2000000000e+00*exp(((V_old_-9.0000000000e+00)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(V_old_-9.0000000000e+00))/2.1600000000e+01))))));
		}else{
			return ((1.0000000000e+00/((4.2000000000e+00*exp(((V_old_-1.0000000000e+01)/1.7000000000e+01)))+(1.5000000000e-01*exp(((-(V_old_-1.0000000000e+01))/2.1600000000e+01))))));
		}
	}
	double Solveode::ifnumber_23(){
		if((Version==0.0000000000e+00)){
			return (2.0000000000e-03);
		}else if((Version==1.0000000000e+00)){
			return (2.0000000000e-03);
		} else{
			return (6.0000000000e-03);
		}
	}
	double Solveode::ifnumber_24(){
		if((Version==0.0000000000e+00)){
			return ((g_K_s_Centre_Published+(calc_FCell*(g_K_s_Periphery_Published-g_K_s_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_K_s_Centre_0DCapable+(calc_FCell*(g_K_s_Periphery_0DCapable-g_K_s_Centre_0DCapable))));
		} else{
			return ((g_K_s_Centre_1DCapable+(calc_FCell*(g_K_s_Periphery_1DCapable-g_K_s_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_25(){
		if((Version==0.0000000000e+00)){
			return ((g_f_Na_Centre_Published+(calc_FCell*(g_f_Na_Periphery_Published-g_f_Na_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_f_Na_Centre_0DCapable+(calc_FCell*(g_f_Na_Periphery_0DCapable-g_f_Na_Centre_0DCapable))));
		} else{
			return ((g_f_Na_Centre_1DCapable+(calc_FCell*(g_f_Na_Periphery_1DCapable-g_f_Na_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_26(){
		if((Version!=2.0000000000e+00)){
			return ((calc_g_f_Na*y_old_*(V_old_-calc_E_Na)));
		}else{
			return ((calc_g_f_Na*y_old_*(V_old_-7.7600000000e+01)));
		}
	}
	double Solveode::ifnumber_27(){
		if((Version==0.0000000000e+00)){
			return ((g_f_K_Centre_Published+(calc_FCell*(g_f_K_Periphery_Published-g_f_K_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_f_K_Centre_0DCapable+(calc_FCell*(g_f_K_Periphery_0DCapable-g_f_K_Centre_0DCapable))));
		} else{
			return ((g_f_K_Centre_1DCapable+(calc_FCell*(g_f_K_Periphery_1DCapable-g_f_K_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_28(){
		if((Version!=2.0000000000e+00)){
			return ((calc_g_f_K*y_old_*(V_old_-calc_E_K)));
		}else{
			return ((calc_g_f_K*y_old_*(V_old_+1.0200000000e+02)));
		}
	}
	double Solveode::ifnumber_29(){
		if((Version==0.0000000000e+00)){
			return ((1.0000000000e+00*exp(((-(V_old_+7.8910000000e+01))/2.6620000000e+01))));
		}else{
			return ((1.0000000000e+00*exp(((-(V_old_+7.8910000000e+01))/2.6630000000e+01))));
		}
	}
	double Solveode::ifnumber_30(){
		if((Version==0.0000000000e+00)){
			return ((g_b_Na_Centre_Published+(calc_FCell*(g_b_Na_Periphery_Published-g_b_Na_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_b_Na_Centre_0DCapable+(calc_FCell*(g_b_Na_Periphery_0DCapable-g_b_Na_Centre_0DCapable))));
		} else{
			return ((g_b_Na_Centre_1DCapable+(calc_FCell*(g_b_Na_Periphery_1DCapable-g_b_Na_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_31(){
		if((Version==0.0000000000e+00)){
			return ((g_b_K_Centre_Published+(calc_FCell*(g_b_K_Periphery_Published-g_b_K_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_b_K_Centre_0DCapable+(calc_FCell*(g_b_K_Periphery_0DCapable-g_b_K_Centre_0DCapable))));
		} else{
			return ((g_b_K_Centre_1DCapable+(calc_FCell*(g_b_K_Periphery_1DCapable-g_b_K_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_32(){
		if((Version==0.0000000000e+00)){
			return ((g_b_Ca_Centre_Published+(calc_FCell*(g_b_Ca_Periphery_Published-g_b_Ca_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((g_b_Ca_Centre_0DCapable+(calc_FCell*(g_b_Ca_Periphery_0DCapable-g_b_Ca_Centre_0DCapable))));
		} else{
			return ((g_b_Ca_Centre_1DCapable+(calc_FCell*(g_b_Ca_Periphery_1DCapable-g_b_Ca_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_33(){
		if((Version==0.0000000000e+00)){
			return ((k_NaCa_Centre_Published+(calc_FCell*(k_NaCa_Periphery_Published-k_NaCa_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((k_NaCa_Centre_0DCapable+(calc_FCell*(k_NaCa_Periphery_0DCapable-k_NaCa_Centre_0DCapable))));
		} else{
			return ((k_NaCa_Centre_1DCapable+(calc_FCell*(k_NaCa_Periphery_1DCapable-k_NaCa_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_34(){
		if((Version==0.0000000000e+00)){
			return (((calc_k_NaCa*((pow(Na_i,3.0000000000e+00)*Ca_o*exp((3.7430000000e-02*V_old_*gamma_NaCa)))-(pow(Na_o,3.0000000000e+00)*Ca_i*exp((3.7400000000e-02*V_old_*(gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(d_NaCa*((Ca_i*pow(Na_o,3.0000000000e+00))+(Ca_o*pow(Na_i,3.0000000000e+00)))))));
		}else{
			return (((calc_k_NaCa*((pow(Na_i,3.0000000000e+00)*Ca_o*exp((3.7430000000e-02*V_old_*gamma_NaCa)))-(pow(Na_o,3.0000000000e+00)*Ca_i*exp((3.7430000000e-02*V_old_*(gamma_NaCa-1.0000000000e+00))))))/(1.0000000000e+00+(d_NaCa*((Ca_i*pow(Na_o,3.0000000000e+00))+(Ca_o*pow(Na_i,3.0000000000e+00)))))));
		}
	}
	double Solveode::ifnumber_35(){
		if((Version==0.0000000000e+00)){
			return ((i_p_max_Centre_Published+(calc_FCell*(i_p_max_Periphery_Published-i_p_max_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((i_p_max_Centre_0DCapable+(calc_FCell*(i_p_max_Periphery_0DCapable-i_p_max_Centre_0DCapable))));
		} else{
			return ((i_p_max_Centre_1DCapable+(calc_FCell*(i_p_max_Periphery_1DCapable-i_p_max_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_36(){
		if((Version==0.0000000000e+00)){
			return ((i_Ca_p_max_Centre_Published+(calc_FCell*(i_Ca_p_max_Periphery_Published-i_Ca_p_max_Centre_Published))));
		}else if((Version==1.0000000000e+00)){
			return ((i_Ca_p_max_Centre_0DCapable+(calc_FCell*(i_Ca_p_max_Periphery_0DCapable-i_Ca_p_max_Centre_0DCapable))));
		} else{
			return ((i_Ca_p_max_Centre_1DCapable+(calc_FCell*(i_Ca_p_max_Periphery_1DCapable-i_Ca_p_max_Centre_1DCapable))));
		}
	}
	double Solveode::ifnumber_37(){
		if((Version==0.0000000000e+00)){
			return ((((R*T)/F)*log(((K_o+(1.2000000000e-01*Na_o))/(K_i+(1.2000000000e-01*Na_i))))));
		}else{
			return ((((R*T)/F)*log(((K_o+(3.0000000000e-02*Na_o))/(K_i+(3.0000000000e-02*Na_i))))));
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
		double dtL, dtM, dtMax=0.0,  dtMin=9990 ;
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
			dependent_variable__ = N_VNew_Serial(15);
			if(check_flag((void *)dependent_variable__, "N_VNew_Serial", 0))
			exit(1);
			depvar__ = (double*)malloc(sizeof(double)*15);
			if(depvar__ == NULL){
			fprintf(stderr, "ERROR Cannot allocate memory for depvar__\n");
			exit(0);
			}
			NV_Ith_S(dependent_variable__, 0) = V_ini_;
			NV_Ith_S(dependent_variable__, 1) = m_ini_;
			NV_Ith_S(dependent_variable__, 2) = h1_ini_;
			NV_Ith_S(dependent_variable__, 3) = h2_ini_;
			NV_Ith_S(dependent_variable__, 4) = d_L_ini_;
			NV_Ith_S(dependent_variable__, 5) = f_L_ini_;
			NV_Ith_S(dependent_variable__, 6) = d_T_ini_;
			NV_Ith_S(dependent_variable__, 7) = f_T_ini_;
			NV_Ith_S(dependent_variable__, 8) = q_ini_;
			NV_Ith_S(dependent_variable__, 9) = r_ini_;
			NV_Ith_S(dependent_variable__, 10) = P_af_ini_;
			NV_Ith_S(dependent_variable__, 11) = P_as_ini_;
			NV_Ith_S(dependent_variable__, 12) = P_i_ini_;
			NV_Ith_S(dependent_variable__, 13) = xs_ini_;
			NV_Ith_S(dependent_variable__, 14) = y_ini_;
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
			flag__ = CVDense(cvode_mem_cvode__, 15);
			if (check_flag(&flag__, "CVDense", 1))	exit(1);
			break;
			case 2:
			flag__ = CVDiag(cvode_mem_cvode__);
			if (check_flag(&flag__, "CVDiag", 1))	exit(1);
			break;
			case 3:
			flag__ = CVBand(cvode_mem_cvode__, 15, NULL, NULL);
			if (check_flag(&flag__, "CVBand", 1))	exit(1);
			break;
			}
			CVodeSetFdata(cvode_mem_cvode__, (void*)this);
			V_old_ = V_ini_;
			m_old_ = m_ini_;
			h1_old_ = h1_ini_;
			h2_old_ = h2_ini_;
			d_L_old_ = d_L_ini_;
			f_L_old_ = f_L_ini_;
			d_T_old_ = d_T_ini_;
			f_T_old_ = f_T_ini_;
			q_old_ = q_ini_;
			r_old_ = r_ini_;
			P_af_old_ = P_af_ini_;
			P_as_old_ = P_as_ini_;
			P_i_old_ = P_i_ini_;
			xs_old_ = xs_ini_;
			y_old_ = y_ini_;
		}
		while(1){
			flag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &time_new, CV_NORMAL);
			V_old_ = NV_Ith_S(dependent_variable__, 0);
			m_old_ = NV_Ith_S(dependent_variable__, 1);
			h1_old_ = NV_Ith_S(dependent_variable__, 2);
			h2_old_ = NV_Ith_S(dependent_variable__, 3);
			d_L_old_ = NV_Ith_S(dependent_variable__, 4);
			f_L_old_ = NV_Ith_S(dependent_variable__, 5);
			d_T_old_ = NV_Ith_S(dependent_variable__, 6);
			f_T_old_ = NV_Ith_S(dependent_variable__, 7);
			q_old_ = NV_Ith_S(dependent_variable__, 8);
			r_old_ = NV_Ith_S(dependent_variable__, 9);
			P_af_old_ = NV_Ith_S(dependent_variable__, 10);
			P_as_old_ = NV_Ith_S(dependent_variable__, 11);
			P_i_old_ = NV_Ith_S(dependent_variable__, 12);
			xs_old_ = NV_Ith_S(dependent_variable__, 13);
			y_old_ = NV_Ith_S(dependent_variable__, 14);
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
	for(int i = 0; i<15; i++)
		ode->setVariables( i ,NV_Ith_S(dependent_variable__, i));
	ode->setParameters(0,time);
	rightHandSideFunction.function(ode);
	for(int i = 0; i<15; i++)
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
