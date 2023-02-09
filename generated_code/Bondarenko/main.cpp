#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <malloc.h>
#include <malloc.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "Solveode.hpp"
// #include "Solveode.Funciona.Overhead.hpp"
#include "Stopwatch.h"

#include<sys/time.h>
#include<sys/resource.h>


#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

double getMemoryUsage(){
    double t;
    struct timeval tim;

    struct rusage ru, ru2;
    getrusage(RUSAGE_SELF, &ru);
    getrusage(RUSAGE_CHILDREN, &ru2);
    
//     tim=ru.ru_stime;
//     t=(double)tim.tv_sec * 1000000.0 + (double)tim.tv_usec;
    printf("Memoria: %u %u\n",ru.ru_maxrss, ru2.ru_maxrss);
    using std::cout;
    using std::endl;
    
    double vm, rss;
//     process_mem_usage(vm, rss);
    cout << "VM: " << vm << "; RSS: " << rss << endl;
    
}
using namespace std;

double
dt = 2e-4,
finalTime= 10000,
maxStep = 1.5,
abstol =0.01, reltol=0.01,
abstol_cvode=1e-6, reltol_cvode=1e-6;

void runMethods(double savingRate, int nThreads, int eqType){
    char fileNameADDT[120];
  
    Stopwatch watch;
    Solveode *addt = new Solveode(abstol, reltol, eqType);
    addt->setFreeVariable(dt);
    addt->setMaxStep(maxStep);
    
    if(savingRate==0)	watch.start();
    
//     sprintf(fileNameADDT, "dat/addt_%d_%d.dat", eqType,nThreads);
//     struct mallinfo inf;
//     struct mallinfo inf2;
//     inf = mallinfo();
//     double mem = ((double)((double)inf.arena + (double)inf.hblkhd)/*/(1024.*1024.)*/);
//     
    
//     int *teste = (int*)malloc(sizeof(int)*1e+5);
    addt->solver(4, finalTime, savingRate/10,fileNameADDT,nThreads);
    
    
//     inf2 = mallinfo();
//     double mem2 = ((double)((double)inf2.arena + (double)inf2.hblkhd)/*/(1024.*1024.)*/);
//     printf("\n%f %f %f\n", mem, mem2, mem2-mem);
    
    if(savingRate==0){
	watch.stop();
	printf("ADDT %d ms\n",watch.timeMS());
    }
    addt->~Solveode();
    free(addt);
    
   
    char fileNameADDT2[120];
    Stopwatch watch2;
    
    Solveode *addt2 = new Solveode(abstol, reltol, eqType);
    addt2->setFreeVariable(dt);
    addt2->setMaxStep(maxStep);
    
    
    if(savingRate==0.0)	watch2.start();
	
	sprintf(fileNameADDT2, "dat/addt2_%d_%d.dat", eqType,nThreads);
//     system("ps aux | grep ./main");
    addt2->solver(5, finalTime, savingRate/10, fileNameADDT2,nThreads);
// 	system("ps aux | grep ./main");
	if(savingRate==0){
	    watch2.stop();
	    printf("ADDT2 tempo %d ms\n",watch2.timeMS());
	}
	
	addt2->~Solveode();
	free(addt2);
//     }



/* 

if(savingRate!=0){
	char folderName[100];
	sprintf(folderName, "dat/img_%d_%d",eqType,nThreads);
	char cmdMkdir[100];
	sprintf(cmdMkdir, "mkdir %s",folderName);
// 	system(cmdMkdir);
	
	char cmdCP2[100];
	sprintf(cmdCP2, "cp dat/*.dat %s",folderName);
// 	system(cmdCP2);
	
	char cmdCP[100];
	sprintf(cmdCP, "mv dat/*.png %s",folderName);
	
	char commandADDT[120];
	sprintf(commandADDT, "./error %s dat/rk2_%d_%d.dat", fileNameADDT, 1,1);
	system(commandADDT);
	
	char commandADDT2[120];
	sprintf(commandADDT2, "./error %s dat/rk2_%d_%d.dat", fileNameADDT2, 1, 1);
	system(commandADDT2);
	
	

    }
    */
}


void run(double savingRate,int nThreads, int eqType){
   
   if(savingRate!=0 && eqType==1 && nThreads==1){
	Stopwatch watch2;
	Solveode *rk2 = new Solveode(abstol, reltol, eqType);
	rk2->setFreeVariable(dt);
	rk2->setMaxStep(maxStep);
	if(savingRate==0)	watch2.start();
	char fileNameRK2[120];
	sprintf(fileNameRK2, "dat/rk2_%d_%d.dat",eqType,nThreads);
	printf("%f\n",savingRate);
// 	getchar();
	rk2->solver(3, finalTime, savingRate/10,fileNameRK2,nThreads);

	if(savingRate==0){
	    watch2.stop();
	    printf("RK2 %d ms\n",watch2.timeMS());
	}
	rk2->~Solveode();
	free(rk2);
	
    }
    
/*
    Stopwatch watch;

	Solveode *EULER = new Solveode(abstol, reltol, eqType);
	EULER->setFreeVariable(dt);
	EULER->setMaxStep(maxStep);
	if(savingRate==0)	watch.start();
	char fileNameEuler[120];
	sprintf(fileNameEuler, "dat/euler_%d_%d.dat",eqType,nThreads);
	EULER->solver(0, finalTime, savingRate,fileNameEuler,nThreads);
	if(savingRate==0){
	watch.stop();
	printf("EULER %d ms\n",watch.timeMS());
	}
	EULER->~Solveode();
	free(EULER);
    char commandEuler[120];

    //sprintf(commandEuler, "./error %s dat/rk2_%d_%d.dat", fileNameEuler, eqType,nThreads);
    sprintf(commandEuler, "./error %s dat/rk2_%d_%d.dat", fileNameEuler, 1,1);
    if(savingRate!=0){
	system(commandEuler);
    }
*/
/*
   if(nThreads==1)
    {
	///cvode   
	int cvodemethod = CV_BDF;
	//cvodemethod = CV_ADAMS;
	double dtCvode =  0.01;
	double mxdtCvode = 0.05;
	int iteracoes = (int)(finalTime/dtCvode);
	int iteracoesSalvas;
	if(savingRate==0)
	    iteracoesSalvas = 1;
	else{
// 	    iteracoesSalvas =  iteracoes;
	    iteracoesSalvas =  finalTime*2;
	}
	int cv_method = 1; //1 CVdense
		//	2 CVDiag
		//	3 CvBand 
	Stopwatch watch3;
	Solveode *cvode = new Solveode(abstol_cvode, reltol_cvode, eqType);
	cvode->setFreeVariable(dtCvode);
	
	cvode->setMaxStep(mxdtCvode);
	if(savingRate==0)	watch3.start();
	char fileName[100];
	sprintf(fileName, "dat/cvode_%d_%d.dat",eqType, nThreads);
// 	system("ps aux | grep ./main");
	
	
// 	struct mallinfo inf;
// 	inf = mallinfo();
// 	double mem = ((double)((double)inf.arena + (double)inf.hblkhd)/(1024.*1024.));
// 	
// 	getMemoryUsage();
// system("ps aux | grep ./main");
	cvode->solveCVODE(1, iteracoes, iteracoesSalvas, cvodemethod, fileName, cv_method);
// 	system("ps aux | grep ./main");
// 	getMemoryUsage();
// 	
// 	struct mallinfo inf2;
// 	inf2 = mallinfo();
// 	double mem2 = ((double)((double)inf2.arena + (double)inf2.hblkhd)/(1024.*1024.));

// 	printf("\n%f %f %f\n", mem, mem2, mem2-mem);
// 	system("ps aux | grep ./main");
	
	
	if(savingRate==0){
	    watch3.stop();
	    printf("CVODE %d ms\n",watch3.timeMS());
	}
	cvode->~Solveode();
	free(cvode);
	if(savingRate!=0){
	    char command[120];
// 	    sprintf(command, "./error %s dat/rk2_%d_%d.dat", fileName, 1,nThreads);
// 	    system(command);

	}
    }*/

//     runMethods(savingRate,nThreads, eqType);
     
}


int main(int argc, char** argv) {
 
    int eqType = atoi(argv[1]);
    //for(int eqType=1; eqType<=5;eqType++)
    {
	if(eqType==1) printf("::Agos:::::::::\n");
	if(eqType==2)  printf("::PycmlSimples:::::::::\n");
	if(eqType==3) printf("::PE:::::::::\n");
	if(eqType==4) printf("::Lut:::::::::\n");
	if(eqType==5) printf("::Pe+Lut:::::::::\n");
	
	
	
	//135168.000000
    }
    if(eqType!=6){
	printf("threads: 1\n");
	run(1, 1, eqType);
// 	run(1, 2, eqType);
	
// 	run(0, 1, eqType);
// 	run(0, 2, eqType);
// 	run(0, 3, eqType);	
// 	run(0, 4, eqType);	
	
// 	
	int numThreads=1;
	if(eqType==1)
	    for(numThreads=1;numThreads<=4;numThreads++)
	    {
// 		printf("threads: %d\n", numThreads);
// 		run(0, numThreads, eqType);
// 		run(0, numThreads, eqType);
// 		run(0, numThreads, eqType);
// 		run(0, numThreads, eqType);
// 		run(0, numThreads, eqType);
// 		run(1, numThreads, eqType);
		
	    }
	    
	
    }else{
	    
	printf("\nJACOBIAN--------\n");
	Solveode *jac = new Solveode(abstol, reltol, 1);
	jac->setFreeVariable(2e-4);
	jac->setMaxStep(0);
	jac->jacobian(1000, 1);
	jac->~Solveode();
	free(jac);
	
    }
//     printf("\n--------\n");
       
}

