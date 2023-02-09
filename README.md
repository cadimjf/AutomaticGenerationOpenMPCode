# AutomaticGenerationOpenMPCode

Please cite this work as 

https://doi.org/10.1007/s00607-012-0268-y


Campos, R.S., Campos, F.O., Gomes, J.M. et al. 
Comparing high performance techniques for the automatic generation of efficient solvers of cardiac cell models. 
Computing 95 (Suppl 1), 639–660 (2013). 


In this work, the AGOS right hand side function was modified in order to use
OpenMP to compute the equations in parallel. 
Since the numerical methods implemented in AGOS are iterative and one step always depends on a previous one, it is not
possible to execute the main loop in parallel. Therefore some adaptations on AGOS
translator were done in order to allow this code to execute in parallel. In summary, the translator analyzes all equations of the system to detect the dependency among them
and automatically determines which equations can be computed in parallel, resulting
in performances up to twice faster than the sequential version of the code.
