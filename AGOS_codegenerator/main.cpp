#include "compile.h"

 
int main(int argc, char **argv) 
{
	xmlDoc *doc = NULL;
	xmlNode *root_element = NULL;
	bool dll = false, cvode = false;

	
	if(argc != 6 && argc != 7 && argc != 8)
	{
		fprintf(stdout, "Builder [src-file] [classname] [dest-data-Diff vars] [dest-data-Parameters] [dest-data-Free var] {-cvode} {-dll}\n");
		exit(1);
	}
	
	
	LIBXML_TEST_VERSION
	doc = xmlReadFile (argv[1], NULL, XML_PARSE_NOBLANKS);
	
	if (doc == NULL)
	{
		printf ("error: could not parse file %s\n", "teste.xml");
		return (2);
	}	
	root_element = xmlDocGetRootElement (doc);
	if(argc == 7)
	if(!strcmp(argv[6],"-dll")){
		dll = true;			
	}else if(!strcmp(argv[6],"-cvode")){
		cvode = true;
	} else{
		printf("Invalid Option %s \n", argv[6]);
		return 1;
	}
	
	if(argc == 8)
	if((!strcmp(argv[6],"-dll") && !strcmp(argv[7],"-cvode")) || (!strcmp(argv[7],"-dll") && !strcmp(argv[6],"-cvode")))
	{
		dll = true;
		cvode = true;
	}else{
		printf("Invalid Options %s and %s\n", argv[6], argv[7]);
		return 1;
	}
		
	
	// inicializando a variavel global
	pointer = root_element;
	if(princ() != SUCCESS){
		fprintf(stderr,"Unknow parsing error\n");
		return 1;
	}
	
	if(rewind_list(difvarlist) == NULL)
	{
		printf("ERROR - Differential Equation not found\n");
		exit(1);
	}
	// building param list and associating initial values
	initial_values(root_element);
	fix_duplicate_equations(root_element, alglist);
	if(!dll && !cvode)
		compile_API(argv[2]);
	else if(dll && !cvode){
		compile_DAPI(argv[2]);
	}else if(!dll && cvode){
		compile_Sundials_API(argv[2]);
	}else {
		compile_Sundials_DAPI(argv[2]);
	}
	saveLists(argv[3],argv[4],argv[5]);		
		
	return 0;
}
