#ifndef _COMPILE_H_
#define _COMPILE_H_
#include "fsm.h"
#include <string.h>

int print_only_diff(FILE *file, DiffList *list);
int print_diff(FILE *file, DiffList *list);
int print_alg(FILE *file, AlgList *list);
int print_right_alg(FILE *file, char *classname, AlgList *list, TokenNode *orderedlist);
int print_if(FILE *file,char* classname, IfList *list);
int print_if2(FILE *file, IfList *list, char *ifNumber, int printVector=0, int printEquationName=1);
int print_eq(FILE *file, TokenNode *t, int printVector=0, int printIfFunction=1, int printEquationName=1);
int	compile_Sundials_API(char *classname);	
int	compile_Sundials_DAPI(char *classname);
int	compile_API(char *classname);
int	compile_DAPI(char *classname);	// dynamic library
void initial_values(xmlNode *root);
int set_initial_value(char *name, double initvalue, char *units);
void saveLists(char* par, char* diff, char *freev);

// Fix duplicated algebraic variables
void fix_duplicate_equations(xmlNode *root, AlgList *alg);
void fix(xmlNode *root, AlgList *list);
void change_duplicated_equations(int eqnumber, int bottom, int top, char *cname, xmlNode *comp);
void fix_node(TokenNode *node, char *oldname, char *newname);

int ifcounter = 0;

char *prefixo="";

typedef struct typ_node_{
    char *equationName;
    void *root;//list of equations that this node needs
    int treeId;
    int flag;
    int height;
}_node_;

typedef struct typ_graph_node_{
    _node_ *equation;
    void *next;
}_list_item_;




int print_only_diff(FILE *file, DiffList *list)
{
    if(file == NULL)
    {
	printf("ERROR - Can't write in file, print_only_diff");
	exit(1);
    }
    DiffList *curl = rewind_list(list);
    TokenNode *cur = NULL;
    int i =0;
    while(curl != NULL)
    {
	fprintf(file,"\n\t\tdepvar__[%d]= ", i);
	cur = curl->diffheader->eq->next;
	print_eq(file, cur);
	fprintf(file,";\t// %d", curl->number);			
	curl = curl->next;
	i++;
    }
    return 0;
}

int print_diff(FILE *file, DiffList *list)
{
    if(file == NULL)
    {
	printf("ERROR - Can't write in file, print_diff");
	exit(1);
    }
    
    DiffList *curl = rewind_list(list);
    TokenNode *cur = NULL;
    while(curl != NULL)
    {
	fprintf(file,"\n\t\t\t%s_new_=d%s*(", curl->diffheader->diffvar.content,curl->diffheader->freevar.content);
	//char teste[255];
	//strcpy (teste, curl->diffheader->diffvar.content);
	cur = curl->diffheader->eq->next;
	print_eq(file, cur);
	fprintf(file,")+%s_old_;\t// %d", curl->diffheader->diffvar.content, curl->number);			
	//fprintf(file,"\n\t\t\tfprintf(stdout, \"%s_new_:  %%f \", %s_new_);\n", teste, teste);
	//fprintf(file,"\n\t\t\tfprintf(stdout, \"TESTE\\n\");\n");
	curl = curl->next;
    }
    //fprintf(file,"\n\t\t\tfprintf(stdout, \"\\n\");");
    return 0;
}

int print_alg(FILE *file, char *classname, AlgList *list)
{
    if(file == NULL)
    {
	printf("ERROR - Can't write in file, print_alg");
	exit(1);
    }
    AlgList *curl = rewind_list(list);
    TokenNode *cur = NULL;
    while(curl != NULL)
    {
	fprintf(file,"\n\tdouble %s::calc_%s(){\n\t\treturn (",classname, curl->eq->token.content);
	cur = curl->eq;
	cur = cur->next->next;
	print_eq(file, cur);
	fprintf(file,");\t//%d\n\t}",curl->number);
	curl = curl->next;
    }
    return 0;
}

int print_if(FILE *file, char *classname, IfList *list)
{
    if(file == NULL)
    {
	printf("ERROR - Can't write in file, print_alg");
	exit(1);
    }
    IfList *curl = rewind_list(list);
    TokenNode *cur = NULL;
    int i = 0;
    while(curl != NULL)
    {
	fprintf(file,"\n\tdouble %s::ifnumber_%d(){\n\t\tif(",classname,curl->ifheader->if_counter);
	// printing the condition
	cur = curl->ifheader->cond;
	print_eq(file, cur);	// printing if condition
	fprintf(file,"){\n\t\t\treturn (");
	cur = curl->ifheader->piece;
	print_eq(file, cur);	// printing the return of if
	fprintf(file,");\n\t\t}");
	while(curl->next != NULL && curl->next->ifheader->if_counter == curl->ifheader->if_counter){
	    curl = curl->next;					
	    fprintf(file,"else if(");
	    
	    cur = curl->ifheader->cond;
	    print_eq(file, cur);	// printing if condition
	    fprintf(file,"){\n\t\t\treturn (");
	    cur = curl->ifheader->piece;
	    print_eq(file, cur);	// printing the return of if
	    fprintf(file,");\n\t\t} ");
	}
	fprintf(file,"else{\n\t\t\t");
	fprintf(file,"return (");
	cur = curl->ifheader->other;
	print_eq(file, cur);	// printing the return of else
	fprintf(file,");\n\t\t}\n\t}");	
	
	i++;
	curl = curl->next;
    }
    return 0;
}


// int print_if2(FILE *file, IfList *list, char *ifNumber, int printVector)
int print_if2(FILE *file, IfList *list, char *ifNumber, int printVector, int printEquationName)
{
    if(file == NULL)
    {
	printf("ERROR - Can't write in file, print_alg");
	exit(1);
    }
    IfList *curl = rewind_list(list);
    TokenNode *cur = NULL;
    int i = 0;
    while(curl != NULL)
    {
	char if_name[120];
	sprintf(if_name, "%d", curl->ifheader->if_counter);
	if( strcmp( if_name, ifNumber) == 0){
	    fprintf(file,"\n");
	    
	    fprintf(file,"(");
	    // printing the condition
	    cur = curl->ifheader->cond;
	    print_eq(file, cur,printVector, 0, printEquationName);	// printing if condition
	    fprintf(file,")\n?(");
	    cur = curl->ifheader->piece;
	    print_eq(file, cur,printVector, 0, printEquationName );	// printing the return of if
	    fprintf(file,")\n:");
	    int cont=0;
	    while(curl->next != NULL && curl->next->ifheader->if_counter == curl->ifheader->if_counter){
		curl = curl->next;					
		fprintf(file,"((");
		cur = curl->ifheader->cond;
		print_eq(file, cur,printVector, 0, printEquationName);	// printing if condition
		fprintf(file,")\n?(");
		cur = curl->ifheader->piece;
		print_eq(file, cur, printVector, 0, printEquationName);	// printing the return of if
		fprintf(file,")\n:");
		fprintf(file,"(");
		cur = curl->ifheader->other;
		print_eq(file, cur, printVector, 0, printEquationName);	// printing the return of else
		fprintf(file,")");	
		fprintf(file,")\n");
		cont++;
	    }
	    if(cont==0){
		fprintf(file,"(");
		cur = curl->ifheader->other;
		print_eq(file, cur, printVector, 0, printEquationName);	// printing the return of else
		fprintf(file,")");	
	    }
	    
	    i++;
	}
	curl = curl->next;
	
    }
    return 0;
}

//funções do grafo de precedencia
int searchEquation(_node_ *graph, int totalSize, char* search){
    for(int i=0;i<totalSize;i++){
	if(strcmp(graph[i].equationName,search)==0){//se as strings sao iguais, strcmp retorna 0
	    return i;
	}
    }
    return -1;
}

//looks up for a node's child
void labelChildren(_node_ *graph, _node_ *node, int totalSize, int idTree);
void labelParents(_node_ *graph, int totalSize, char* search, int idTree);

void labelChildren(_node_ *graph, _node_ *node,int totalSize, int idTree){
    if(node->treeId==0){
	node->treeId = idTree;
    }
    if(node->flag==0){
	node->flag = 1;
    }
    _list_item_ *auxAchaUltimo = (_list_item_*)node->root;
    while(auxAchaUltimo!=NULL){
	if(auxAchaUltimo->equation->flag == 0 ){
	    labelChildren(graph, auxAchaUltimo->equation, totalSize, idTree);
	}
	auxAchaUltimo = (_list_item_*)auxAchaUltimo->next;
    }
    labelParents(graph, totalSize, node->equationName, idTree);
    
}

void labelParents(_node_ *graph, int totalSize, char* search, int idTree){
    for(int i=0;i<totalSize;i++){
	_list_item_ *auxAchaUltimo = (_list_item_*)graph[i].root;
	while(auxAchaUltimo!=NULL){
	    if(strcmp(auxAchaUltimo->equation->equationName, search)==0){
		labelChildren(graph, &graph[i], totalSize, idTree);
	    }
	    auxAchaUltimo = (_list_item_*)auxAchaUltimo->next;
	}
    }
}



int classifyTrees(_node_ *graph, int totalSize){
    int lastId = 0;
    for(int i=0;i<totalSize;i++){
	
	if ( graph[i].treeId==0 ){
	    lastId++;
	    labelChildren(graph, &graph[i], totalSize, lastId);
	}
	
    }
    return lastId;
}

void insereDependencia(TokenNode *t, _node_* graph, int nodeId, int totalSize){
    _node_ *equationItem = &graph[nodeId];
    //encontra o indice da equação no grafo
    int index = searchEquation(graph, totalSize, t->token.content);
    _list_item_ *aux;
    if(index!=-1){
	//monta um item da lista de depencia, com um ponteiro para a equaçao
	aux = new _list_item_;
	aux->equation = &graph[index];
	aux->next = NULL;
    }else{ //caso nao tenha nenhuma dependencia, retorna -1 e o root será NULL;
    aux = NULL;
    }
    //encontra o ultimo item da lista de depencias
    _list_item_ *auxAchaUltimo = (_list_item_*) equationItem->root;
    if(auxAchaUltimo==NULL){ //a raiz é nula
	equationItem->root = aux;
    }else{//a raiz nao eh nula
	while(auxAchaUltimo->next!=NULL){
	    auxAchaUltimo = (_list_item_*)auxAchaUltimo->next;
	}
	//adiciona na ultima posicao
	auxAchaUltimo->next = (void*)aux;
    }
}


void criaListaPrecedencia(TokenNode *t, _node_* graph, int nodeId, int totalSize)
{
    TokenNode *curalg = rewind_list(algvarlist);
    //searching in alglist
    while(curalg != NULL)
    {
	printf("%s\n", curalg->token.content);
	
	curalg = curalg->next;
    }
    printf("--------****************---\n");
    
    
    
    
    
    
    _node_ *equationItem = &graph[nodeId];
    TokenNode *cur = t;
    while(cur != NULL)
    {	
	printf("%s ::::: %s %d\n",  cur->token.content, equationItem->equationName,cur->token.type);
	// a variavel está na lista de variaveis algebricas
	if(list_has_var(cur->token, algvarlist))
	{
	    printf("%s ::::: %s\n",  cur->token.content, equationItem->equationName);
	    insereDependencia(cur, graph, nodeId, totalSize);
	}
// 	a variavel eh do tipo piecewise
	else if(cur->token.type == PI_W){
	    //cur->token.content contem o numero do If, que está na lista de dependencia da equação
	    IfList *curl  = rewind_list(iflist);
	    TokenNode *curIF = NULL;
	    int i = 0;
	    while(curl != NULL)
	    {
		char _if_number_[3];
		sprintf(_if_number_, "%d", curl->ifheader->if_counter);
		//cur->token.content contem o numero do If, que está na lista de dependencia da equação
		//curl->ifheader->if_counter contem o numero do if, que esta na lista que esta sendo percorrida nesta função
		if( strcmp( _if_number_, cur->token.content)==0){
		    curIF = curl->ifheader->cond;
		    while(curIF != NULL)
		    {	
			if(list_has_var(curIF->token, algvarlist))
			{
			    insereDependencia(curIF, graph, nodeId, totalSize);
			}
			curIF = curIF->next;
		    }
		    
		    curIF = curl->ifheader->piece;
		    while(curIF != NULL)
		    {	
			if(list_has_var(curIF->token, algvarlist))
			{
			    insereDependencia(curIF, graph, nodeId, totalSize);
			}
			curIF = curIF->next;
		    }
		    
		    curIF = curl->ifheader->other;
		    while(curIF != NULL)
		    {	
			if(list_has_var(curIF->token, algvarlist))
			{
			    insereDependencia(curIF, graph, nodeId, totalSize);
			}
			curIF = curIF->next;
		    }
		    
		}
		// printing the condition
		curIF = curl->ifheader->cond;
		curIF = curl->ifheader->piece;
		while(curl->next != NULL && curl->next->ifheader->if_counter == curl->ifheader->if_counter){
		    curl = curl->next;					
		    curIF = curl->ifheader->cond;
		    curIF = curl->ifheader->piece;
		}
		curIF = curl->ifheader->other;
		curl = curl->next;
	    }
	// a variavel é a
	}else {
	}
	cur = cur->next;
    }
}

void printGrafo(_node_ *graph, int totalSize)
{
    for(int i=0;i<totalSize;i++){
	printf("[%d] %s (%d):", graph[i].height, graph[i].equationName, graph[i].treeId);
	_list_item_ *auxAchaUltimo = (_list_item_*)graph[i].root;
	while(auxAchaUltimo!=NULL){
	    printf(" -> %s:", auxAchaUltimo->equation->equationName);	
	    auxAchaUltimo = (_list_item_*)auxAchaUltimo->next;
	}
	printf("\n");
    }
}


int counTreeItems(_node_ *graph, int totalSize, int idTree)
{
    int count=0;
    for(int i=0;i<totalSize;i++){
	if(graph[i].treeId == idTree)
	    count++;
    }
    return count;
}

/**
**********************
* FILE *file:	arquivo onde as equações serão escritas
**********************
* TokenNode *t:	estrutura de dados do AGOS
**********************
* int printVector:	se igual a 1, imprime o vetor OLD[i] para as variáveis diferenciaveis.
*			se igual a 0, imprime o nome da variavel: V_old_
**********************
* int printIfFunction: se igual a 1, imprime uma chamada de função, para as equações condicionais.
*			se igual a 0, imprime as equações condicionais com o formato 
*			(condição)?(equação se verdadeiro):(equação se falso);
**********************						
* int printEquationName: se igual a 1, imprime o nome da equação, quando está está sendo usada no lado direito.
*			 se igual a 0, imprime a equação em si, e não o nome. Isto é utilizado para gerar equações sem avaliação parcial
**********************
*/
int print_eq(FILE *file, TokenNode *t, int printVector, int printIfFunction, int printEquationName)
{
    if(file == NULL)
    {
	printf("ERROR - Can't write in file, print_alg");
	exit(1);
    }	
    TokenNode *cur = t;
    while(cur != NULL)
    {	
	if(list_has_var(cur->token, algvarlist))
	{
	    //imprimir a equação ao inves do nome dela
	    if(printEquationName==0){
		TokenNode *curlCalcALG = rewind_list(algvarlist);
		TokenNode *curCalcALG  = NULL;
		AlgList *curalgCalcALG = NULL;
		while(curlCalcALG != NULL)
		{
		    if(strcmp(cur->token.content, curlCalcALG->token.content)==0){
			fprintf(file, "(");//abre parenteses, para garantir que a equação será computada com a precedencia correta
			curalgCalcALG = rewind_list(alglist);
			while(strcmp(curalgCalcALG->eq->token.content, curlCalcALG->token.content)){
			    curalgCalcALG = curalgCalcALG->next;
			}
			curCalcALG = curalgCalcALG->eq;
			curCalcALG = curCalcALG->next->next;
			//print_eq(file, curCalcALG);
			print_eq(file, curCalcALG, printVector, printIfFunction, printEquationName);
			//final da funcao
			curlCalcALG = curlCalcALG->next;
			fprintf(file, ")");//fecha parenteses
			break;
		    }else{
			    curCalcALG = curCalcALG->next->next;
			    curlCalcALG = curlCalcALG->next;
		    }
		}
	    }
	    else//imprime o nome da equação
	    {
		fprintf(file, "%scalc_",prefixo);
		fprintf(file, "%s", cur->token.content);
	    }
	}else{
	    if(cur->token.type == PI_W){
		if(printIfFunction==1){
		    fprintf(file, "%sifnumber_%s()", prefixo, cur->token.content);
		}else{
		    print_if2(file, iflist, cur->token.content, printVector, printEquationName);
		}
		//fprintf(file, "ifnumber_%d()", p_if_counter);
	    }else{
		if(!strcmp(cur->token.content, difflist->diffheader->freevar.content)){
		    fprintf(file,"%s%s_new",prefixo,cur->token.content);		
		}else{
		    if(!cur->token.isdiff || cur->token.type == CPAR || cur->token.type == OPAR){
			if(printVector==0){
			    if(cur->token.type == VARI){
				fprintf(file, "%s%s",prefixo, cur->token.content);
			    }else
				fprintf(file, "%s",cur->token.content);//nao por prefixo = numeros e sinais
			}else{
			    // TODO 22/jun/11 -> imprimir o vetor com indice ao inves de VARDIFD_OLD_
			    if(!list_has_var(cur->token, difvarlist)){
				if( cur->token.type == VARI ){
				    fprintf(file, "%s%s",prefixo, cur->token.content);
				}else{
				    fprintf(file, "%s",cur->token.content);//nao por prefixo = numeros e sinais
				}
			    }
			}
		    }
		    if(list_has_var(cur->token, difvarlist))
		    {
			
			if(!cur->token.isdiff){
			    if(printVector==0){
				fprintf(file, "_old_");
			    }else{//TODO 22/jun/11 -> vetor con indice
				int indice =0;
				DiffList *curdiff = rewind_list(difflist);
				while(curdiff!=NULL)
				{
				    if(!strcmp((const char*)curdiff->diffheader->diffvar.content, cur->token.content)){
					break;
				    }
				    curdiff = curdiff->next;
				    indice++;
				}
				
				fprintf(file, "__OLD_[%d]",indice);//nao por prefixo = numeros e sinais
			    }
			}else{
			    DiffList *curdiff = rewind_list(difflist);
			    int found = 0;
			    while(curdiff!=NULL && !found)
			    {
				if(!strcmp((const char*)curdiff->diffheader->diffvar.content, cur->token.content))
				    found = 1;
				else curdiff = curdiff->next;
			    }
			    if(found)
			    {
				printf("2:%s\n",cur->token.content);/// //////////////////////
				
				print_eq(file, curdiff->diffheader->eq->next,printVector,printIfFunction, printEquationName);
			    }else{
				printf("ERROR - printfDiff Differential equation referenced, but not defined\n");
				exit(1);								
			    }
			}
		    }
		}
	    }
	}
	cur = cur->next;
    }
    
    return 0;
}

void initial_values(xmlNode *root)
{
    parvarlist = create_param_list(varlist, difvarlist, algvarlist);
    xmlNode *curmodel = root;
    xmlNode *curvar;
    
    while(curmodel != NULL)
    {
	// if the file is just a MathML not a CellML
	if(!strcmp((char*)curmodel, "math"))
	    break;
	curmodel = get_component(curmodel);
	if(curmodel != NULL)
	{
	    curvar = get_variable(curmodel);
	    while(curvar != NULL)
	    {
		char *name = "";
		double initvalue = 0.0;
		char *units = "";
		// looking for initial_values in the variable
		
		if(curvar!= NULL){
		    xmlAttr *attr = curvar->properties;
		    while(attr!= NULL)
		    {
			if(!strcmp("initial_value", (char*)attr->name))
			    break;
			attr = attr->next;
		    }
		    // if a initial value was found
		    if(attr != NULL)
		    {
			initvalue = atof((char*)xmlNodeGetContent((xmlNode*)attr));
		    }	
		    // take the variable name
		    attr = curvar->properties;
		    while(attr!= NULL)
		    {
			if(!strcmp("name", (char*)attr->name))
			    break;
			attr = attr->next;
		    }
		    if(attr != NULL)
		    {
			name= (char*)xmlNodeGetContent((xmlNode*)attr);
		    }	
		    attr = curvar->properties;
		    while(attr!= NULL)
		    {
			if(!strcmp("units", (char*)attr->name))
			    break;
			attr = attr->next;
		    }
		    if(attr != NULL)
		    {
			units = (char*)xmlNodeGetContent((xmlNode*)attr);
		    }
		    set_initial_value(name, initvalue, units);
		    
		}
		curvar = get_variable(curvar);
	    }
	}
    }
    
}

int set_initial_value(char *name, double initvalue, char *units)
{
    // creating the param list
    
    TokenNode *curdif = rewind_list(difvarlist);
    
    
    int found = 0;
    // searching in difflist
    while(curdif != NULL & !found)
    {
	if(!strcmp(strNoSpace(name), strNoSpace(curdif->token.content)))
	{
	    curdif->units = units;
	    if(initvalue != 0.0)
		curdif->initialvalue = initvalue;
	    found = 1;
	}
	
	curdif = curdif->next;
    }
    TokenNode *curalg = rewind_list(algvarlist);
    //searching in alglist
    while(curalg != NULL & !found)
    {
	if(!strcmp(name, strNoSpace(curalg->token.content)))
	{
	    curalg->units = units;
	    if(initvalue != 0.0)
		curalg->initialvalue = initvalue;
	    found = 1;
	}
	
	curalg = curalg->next;
    }
    TokenNode *curpar = rewind_list(parvarlist);
    //searching in parlist
    while(curpar != NULL & !found)
    {
	if(!strcmp(name, strNoSpace(curpar->token.content)))
	{
	    curpar->units = units;
	    if(initvalue != 0.0)
		curpar->initialvalue = initvalue;
	    found = 1;
	}
	curpar = curpar->next;
    }			
    if(!found) return 0;
    else return 1;
}

void saveLists(char* diff, char* par, char *freev)
{
    TokenNode *curpar;
    
    char fname[100], funits[100],finit[100];
    memset(fname,100,sizeof(char));
    sprintf(fname, "%s", par);
    sprintf(funits, "%su", par);
    sprintf(finit, "%sv", par);
    
    FILE *fn = fopen(fname, "wb");
    FILE *fu = fopen(funits, "wb");
    FILE *fv = fopen(finit, "wb");
    
    if(!(fn && fu && fv))
    {
	printf("ERROR - Can't create files - saveLists 1\n");
	exit(1);
    }
    curpar = rewind_list(parvarlist);
    while(curpar != NULL)
    {
	if(curpar->units == NULL)
	    curpar->units = "";
	fprintf(fn, "%s ",curpar->token.content);
	fprintf(fu, "%s ",curpar->units);
	fprintf(fv, "%.8e ",curpar->initialvalue);
	curpar = curpar->next;
    }
    fclose(fn);
    fclose(fu);
    fclose(fv);
    sprintf(fname, "%s", diff);
    sprintf(funits, "%su", diff);
    sprintf(finit, "%sv", diff);
    fn = fopen(fname, "wb");
    fu = fopen(funits, "wb");
    fv = fopen(finit, "wb");
    if(!(fn && fu && fv))
    {
	printf("ERROR - Can't create files - saveLists 2\n");
	exit(1);
    }
    curpar = rewind_list(difvarlist);
    while(curpar != NULL)
    {
	if(curpar->units == NULL)
	    curpar->units = "";
	fprintf(fn, "%s ",curpar->token.content);
	fprintf(fu, "%s ",curpar->units);
	fprintf(fv, "%.8e ",curpar->initialvalue);		
	curpar = curpar->next;
    }
    fclose(fn);
    fclose(fu);
    fclose(fv);
    sprintf(fname, "%s", freev);
    sprintf(funits, "%su", freev);
    sprintf(finit, "%sv", freev);
    fn = fopen(fname, "wb");
    fu = fopen(funits, "wb");
    fv = fopen(finit, "wb");
    
    if(!(fn && fu && fv))
    {
	printf("ERROR - Can't create files - saveLists 3\n");
	exit(1);
    }
    curpar = rewind_list(parvarlist);
    
    fprintf(fn, "%s ",curpar->token.content);
    fprintf(fu, "%s ",curpar->units);
    fprintf(fv, "%.8e ",curpar->initialvalue);		
    
    fclose(fn);
    fclose(fu);
    fclose(fv);	
}

void fix_duplicate_equations(xmlNode *root, AlgList *alg)
{
    xmlNode *curroot = root;
    TokenNode *list = NULL;
    AlgList *dpllist = NULL;
    AlgList *curalg = rewind_list(alg);
    
    // creating a list of duplicated equations
    while(curalg != NULL)
    {
	Token t = curalg->eq->token;
	if(!list_has_var(t, list))
	    list = attach_token(t, list);
	else{
	    dpllist = attach_alg_eq(NULL, dpllist);	// PODE DAR PAU AKI
	    dpllist->number = curalg->number;
	}
	curalg = curalg->next;
    }
    
    dpllist = rewind_list(dpllist);		
    fix(curroot, dpllist);
}

void fix(xmlNode *root, AlgList *list)
{
    xmlNode *curroot = root;
    AlgList *dpllist = list;
    int counter = -1;	// to count the equations
    //searching for the duplicated equation
    while((curroot = get_component(curroot)) != NULL)
    {
	xmlNode *math = get_math(curroot);		
	xmlNode *apply = NULL;
	if(math !=NULL)
	    apply = math->children;
	int component = counter;
	while(apply!=NULL)
	{
	    if(!strcmp((char*)apply->name, "apply"))
		counter++;
	    apply=apply->next;
	}
	// getting the component name*/
	    xmlAttr *name = curroot->properties;
	    char *cname = NULL;
	    while(name != NULL)
	    {
		if(!strcmp((char*)name->name, "name"))
		    break;
		name=name->next;
	    }
	    if(name != NULL)
	    {
		//printf("--------%s\n",name->children->content);
		cname = (char*)xmlNodeGetContent((xmlNode *)name);
	    }
	    if(dpllist != NULL)
		
		while(counter >= dpllist->number && component < dpllist->number)
		{
		    printf("A duplicated algebraic equation was found at component %s\nTrying to fix the problem\nPlease verify if the equations generated are correct\n", cname);
		    
		    change_duplicated_equations(dpllist->number,component+1, counter, cname, curroot);
		    dpllist = dpllist->next;
		    if(dpllist == NULL)
			break;
		}
		
    }		
}

void change_duplicated_equations(int eqnumber, int bottom, int top, char *cname, xmlNode *comp)
{
    AlgList *alg = rewind_list(alglist);
    DiffList *diff = rewind_list(difflist);
    char *oldname = NULL;
    
    while(alg != NULL)
    {
	if(alg->number == eqnumber)
	    break;
	alg = alg->next;
    }
    if(alg != NULL)
	oldname = alg->eq->token.content;
    
    char *newname = (char*)calloc(strlen(oldname)+strlen(cname)+13, sizeof(char));
    sprintf(newname,"%s_duplicated_%s",oldname,cname);
    
    Token t;
    t.content = newname;
    algvarlist = attach_token(t,algvarlist);
    /*
    xmlNode *var = comp;
    while((var = get_variable(var)) != NULL)
		{
		    xmlAttr *attr = var->properties;
		    while(attr != NULL)
		{
		    if(!strcmp((char*)attr->name, "name"))
			break;
		    attr = attr->next;
		}
		if(attr != NULL)
		{
		    attr->children->content = (xmlChar*)newname;
		    break;
		}
		}
		*/
    alg = rewind_list(alglist);
    
    while(alg != NULL)
    {
	if(alg->number >= bottom && alg->number <= top)
	{
	    TokenNode *cur = alg->eq;
	    fix_node(cur, oldname, newname);
	}		
	if(alg->number > top)
	    break;
	alg = alg->next;
    }
    
    while(diff!=NULL)
    {
	if(diff->number >= bottom && diff->number <= top)
	{
	    TokenNode *cur = diff->diffheader->eq;
	    fix_node(cur, oldname, newname);
	}		
	if(diff->number > top)
	    break;
	diff = diff->next;
    }
}

void fix_node(TokenNode *node, char *oldname, char *newname)
{
    TokenNode *cur = node;
    while(cur != NULL)
    {
	if(cur->token.type == PI_W)
	{
	    int num = atoi(cur->token.content);
	    //printf("IFNUMBER %d\n", num);
	    IfList *ilist = rewind_list(iflist);
	    int i;
	    
	    for(i=0;i<num;i++)
		ilist = ilist->next;
	    fix_node(ilist->ifheader->cond, oldname, newname);
	    fix_node(ilist->ifheader->piece, oldname, newname);
	    fix_node(ilist->ifheader->other, oldname, newname);
	}
	
	if(!strcmp(cur->token.content, oldname))
	    cur->token.content = newname;
	cur = cur->next;
    }
}

int	compile_API(char *classname)
{
    
    if(eq_counter <= 0)
    {
	printf("ERROR - Equation not found\n");
	exit(1);
    }
    
    
    // teste da insercao de precedencia 08/03/2008
    preced_alg_list = create_preced_alg_list(rewind_list(alglist));
    AlgList *cural = rewind_list(preced_alg_list);
    TokenNode *resolved_dep_list = NULL;
    
    TokenNode *cur = rewind_list(algvarlist);
    
    
    
    while(preced_alg_list != NULL)
    {
	cural = rewind_list(preced_alg_list);
	while(cural != NULL){
	    //	printf("%s $$$$$\n", cural->eq->token.content);
	    cural = cural->next;
	}
	cural = rewind_list(preced_alg_list);
	
	//	printf("ASDFASDFASDFAS\n");
	while(cural != NULL){
	    
	    TokenNode *cureq = cural->eq->next;
	    //	printf("%s \n", cural->eq->token.content);
	    //				getchar();
	    if(cureq == NULL){
		resolved_dep_list = add_list(cural->eq->token, resolved_dep_list);
		preced_alg_list = delete_from_list(preced_alg_list, cural->eq->token);
		cural = cural->next;
		//					
		//					printf("Lista %s\n", resolved_dep_list->token.content);
		//					getchar();
	    }else{
		
		while(cureq !=NULL){
		    if(list_has_var(cureq->token, resolved_dep_list)){
			//printf("%s \n", cureq->token.content);
			//getchar();
			if(cureq->prev != NULL){
			    cureq->prev->next = cureq->next;
			}
			if(cureq->next != NULL){
			    cureq->next->prev = cureq->prev;
			}
		    }
		    cureq = cureq->next;
		}
		cural = cural->next;					
	    }	
	    
	}
	
    }
    //		TokenNode *currl = rewind_list(resolved_dep_list);
    //					while(currl != NULL){
	//						printf("LRD %s\n", currl->token.content);
	//						currl = currl->next;
	//		}
	// fim - teste da insercao de precedencia 08/03/2008
	FILE *file;
	
	char *filename = (char *)calloc(strlen(classname)+5, sizeof(char*));
	
	sprintf(filename,"%s.hpp",(const char*)classname);
	
	file = fopen(filename, "wb");
	
	fprintf(file, "#include <stdio.h>\n#include <math.h>\n#include \"MCutil.hpp\"\n");
	//fprintf(file, "#include <omp.h>\n#include \"Stopwatch.h\"\n");
	fprintf(file, "#define _H_FACTOR_ 2\n");
	fprintf(file, "#define RECURSION_LIMIT 10\n");
	
	fprintf(file, "#define _EULER_ 0\n");
	fprintf(file, "#define _ADAP_DT_ 4\n");
	fprintf(file, "#define _RK2_ 3\n");
	
	fprintf(file, "#include <omp.h>\n");
	fprintf(file, "#define _LOWER_LIMIT_   0.2\n");
	fprintf(file, "#define _UPPER_LIMIT_   1.0\n");
	
	
	int numCalcs =0;
	while(cur != NULL)
	{
	    cur = cur->next;
	    numCalcs++;
	}
	
	cur = rewind_list(difvarlist);
	int numEdos=0;
	while(cur != NULL)
	{
	    cur = cur->next;
	    numEdos++;
	}
	
	int totalSize = numCalcs + numEdos;
	
	//grafo que armazena as dependencias de cada equação
	_node_ *grafoDependencias = (_node_*) malloc( (totalSize) * sizeof (_node_) );
	
	cur = rewind_list(algvarlist);
	numCalcs =0;
	while(cur != NULL)
	{
	    TokenNode *curlCalc = rewind_list(resolved_dep_list);
	    TokenNode *curCalc = NULL;
	    AlgList *curalgCalc = NULL;
	    int cont=0;
	    while(curlCalc != NULL)
	    {
		if(cont==numCalcs){
		    curalgCalc = rewind_list(alglist);
		    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
			curalgCalc = curalgCalc->next;
		    }
		    // 		printf("%s ->", curalgCalc->eq->token.content);
		    grafoDependencias[numCalcs].equationName = curalgCalc->eq->token.content; 
		    grafoDependencias[numCalcs].root = NULL;
		    grafoDependencias[numCalcs].treeId = 0;
		    grafoDependencias[numCalcs].flag = 0;
		    grafoDependencias[numCalcs].height = -1;
		    curCalc = curalgCalc->eq;
		    curCalc = curCalc->next->next;
		    
		    criaListaPrecedencia(curCalc,grafoDependencias, numCalcs, totalSize);
		    
		    curlCalc = curlCalc->next;
		    break;
		}else{
		    curCalc = curCalc->next->next;
		    curlCalc = curlCalc->next;
		}
		cont++;
	    }
	    cur = cur->next;
	    numCalcs++;
	    // 	printf("\n");
	}
	
	
	cur = rewind_list(difvarlist);
	numEdos=0;
	while(cur != NULL)
	{
	    DiffList *curlDiff = rewind_list(difflist);
	    TokenNode *curDiff = NULL;
	    int cont=0;
	    while(curlDiff != NULL)
	    {
		if(cont==numEdos){
		    curDiff = curlDiff->diffheader->eq->next;
		    // 		printf("::::\t\t%s \n", curlDiff->diffheader->diffvar.content);
		    grafoDependencias[numCalcs + numEdos].equationName = curlDiff->diffheader->diffvar.content;
		    grafoDependencias[numCalcs + numEdos].root = NULL;
		    grafoDependencias[numCalcs + numEdos].treeId = 0;
		    grafoDependencias[numCalcs + numEdos].flag = 0;
		    grafoDependencias[numCalcs+ numEdos].height = -1;
		    criaListaPrecedencia(curDiff, grafoDependencias, numCalcs + numEdos, totalSize);
		    curlDiff = curlDiff->next;
		    break;
		}else{
		    curDiff = curlDiff->diffheader->eq->next;
		    curlDiff = curlDiff->next;
		}
		cont++;
	    }
	    cur = cur->next;
	    numEdos++;
	    // 	printf("\n");
	}
	//     	printf("%d\n", totalSize);
	int treeNumber = classifyTrees(grafoDependencias, totalSize);
	//     printf("treeNumber: %d\n", treeNumber);
	    printGrafo(grafoDependencias, totalSize);

	
	
	
	
	
	fprintf(file,"#define numEDO %d\n",numEdos);
	fprintf(file,"#define numAux %d\n",numCalcs);
	
	fprintf(file, "\nclass %s\n{\npublic:\n", classname);
	fputs("\tint it_countx;\n",file);
	
	
	/******* DECLARACAO DAS VARIAVEIS *******/
	// building parvaslist
	
	// 	    TokenNode *cur = rewind_list(parvarlist);
	cur = rewind_list(parvarlist);
	while(cur != NULL)
	{
	    fprintf(file,"\tdouble %s; \t // %s\n",cur->token.content, cur->units);
	    //fprintf(file,"\tdouble %s_new;\n",cur->name);
	    cur = cur->next;
	}
	
	cur = rewind_list(algvarlist);
	while(cur != NULL)
	{
	    fprintf(file,"\tdouble calc_%s; \t // %s\n",cur->token.content, cur->units);
	    //fprintf(file,"\tdouble %s_new;\n",cur->name);
	    cur = cur->next;
	}
	
	
	fprintf(file,"\tdouble d%s, *%s_vec__;\n",difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
	fprintf(file,"\tdouble %s_new;\n",difflist->diffheader->freevar.content);
	
	fprintf(file,"\n\t//functions variables\n");
	cur = rewind_list(difvarlist);
	while(cur != NULL)
	{
	    fprintf(file,"\tdouble *%s;\n", cur->token.content);
	    fprintf(file,"\tdouble %s_new_, %s_old_, %s_ini_, %s_lado_direito_;\n", cur->token.content,cur->token.content,cur->token.content, cur->token.content);
	    cur = cur->next;
	}
	/******* FIM * DECLARACAO DAS VARIAVEIS *******/
	
	fprintf(file,"\npublic:\n");
	fprintf(file, "\t%s(double);\n", classname);
	fprintf(file, "\t~%s();\n", classname);
	
	
	///---------------ADDT----------------------
	fprintf(file, "\tdouble STEP_TOLERANCE_;\n");
	fprintf(file, "\tvoid getRightHandSide();\n");
	fprintf(file, "\tvoid getRightHandSideParallel();\n");    
	fprintf(file, "\tvoid setMaxStep(double);\n");
	fprintf(file, "\tint getErrorCode(double);\n");
	fprintf(file, "\tdouble solveToFileOMP(double finalTime = 0, double savingRate = 0, int method=1, char *fileName=\"output.dat\", int nThreads=1);\n");
	fprintf(file, "\tvoid explicitEulerStep();\n");
	fprintf(file, "\tvoid rungekutta2_Parallel( );\n");
	fprintf(file, "\tvoid rungekutta2_Single( );\n");
	fprintf(file, "\tvoid adaptiveDt_Parallel();\n");
	fprintf(file, "\tint adaptiveDt_Single();\n");
	
	fprintf(file, "\tint adaptiveDt();\n");
	//fprintf(file, "\tint adaptiveTrapezoidalStep();\n");
	fprintf(file, "\tdouble errorAux;\n");
	
	fprintf(file, "\tint count1; int count0;\n");
	fprintf(file, "\tint count3; int count4;\n");
	fprintf(file, "\tint recursion_counter;\n");
	fprintf(file, "\tdouble maxStep;\n");
	///-------------------------------------
	
	fputs("\tint setVariables(int, double);\n",file);
	fputs("\tint setParameters(int, double);\n",file);
	fputs("\tint setFreeVariable(double);\n",file);
	
	fputs("\tdouble getVariables(int);\n",file);
	fputs("\tdouble getLadoDireito(int);\n",file);
	fputs("\tvoid executeTree(void* tree, int treeSize);\n",file);
	fputs("\tdouble getParameters(int);\n",file);
	fputs("\tdouble getFreeVariable();\n",file);
	
	fputs("\tVariables get_Parameters();\n",file);
	fputs("\tVariables get_Variables();\n",file);
	fputs("\tVariables get_FreeVariable();\n",file);
	
	fputs("\tvoid setParametersFromFile(char*);\n", file);
	fputs("\tvoid setVariablesFromFile(char*);\n", file);
	fputs("\tvoid setFreeVariableFromFile(char*);\n", file);
	
	fputs("\tdouble solve(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0, int numThreads=1);\n",file);
	fputs("\tdouble jacobian(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0, int numThreads=1);\n",file);
	fputs("\tdouble* getSolution(int indVariable);\n",file);
	fputs("\tdouble* getIndependentVar();\n",file);
	fputs("\tdouble solveToFile(char *filename, char *fileaccess = \"\", int firstcall__ = 0, int num_iterations__ = 0, int num_results__ = 0);\n",file);
	fputs("public:\n",file);
	
	int iif;
	for(iif =0; iif< ifs; iif++)
	{
	    fprintf(file, "\tinline double ifnumber_%d();\n", iif);
	}
	
	fprintf(file,"private:\n");
	fprintf(file,"\tdouble errorAD_;\n");
	fprintf(file,"\tdouble edos_euler_[numEDO];\n");
	fprintf(file,"\tdouble edos_rk2_[numEDO];\n");
	fprintf(file,"\tdouble edos_fn_[numEDO];\n");
	fprintf(file,"\tdouble ladodir_fn_[numEDO];\n");
	
	fprintf(file,"};\n\n");
	
	fprintf(file,"#define AGOS_NAN (0.0/0.0)\n#define AGOS_INF (1.0/0.0)\n");
	fprintf(file,"#define __agos_xor(a,b) (!(a && b) && (a || b))\nfloat __agos_factorial(int);\n");
	fprintf(file,"double _agos_max(double*,int);\n");
	fprintf(file,"double _agos_min(double*,int);\n");
	fprintf(file,"double _agos_round( double, int);\n");
	
	
	//alterado por cadim, 11/09/2009
	//colocando o openmp
	fprintf(file," typedef struct str_funcao{\n\tvoid (*funcao)(Solveode*);\n}typ_funcao;\n\n");
	
	
	
	fprintf(file,"typedef struct str_forest{\n");
	fprintf(file,"typ_funcao* tree;\n");
	fprintf(file,"int treeSize;\n");
	fprintf(file,"}typ_forest;\n");
	
	prefixo="classe->";
	//imprime as funções auxiliares
	
	cur = rewind_list(algvarlist);
	numCalcs =0;
	while(cur != NULL)
	{
	    TokenNode *curlCalc = rewind_list(resolved_dep_list);
	    TokenNode *curCalc = NULL;
	    AlgList *curalgCalc = NULL;
	    int cont=0;
	    
	    while(curlCalc != NULL)
	    {
		if(cont==numCalcs){
		    curalgCalc = rewind_list(alglist);
		    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
			curalgCalc = curalgCalc->next;
		    }
		    fprintf(file,"\nvoid __%s__(%s *classe){\n", curalgCalc->eq->token.content, classname);
		    fprintf(file,"\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
		    curCalc = curalgCalc->eq;
		    curCalc = curCalc->next->next;
		    print_eq(file, curCalc);
		    fprintf(file,";\t//%d",curalgCalc->number);
		    //final da funcao
		    fprintf(file,"\n}\n");
		    curlCalc = curlCalc->next;
		    break;
		}else{
		    curCalc = curCalc->next->next;
		    curlCalc = curlCalc->next;
		}
		cont++;
	    }
	    
	    cur = cur->next;
	    numCalcs++;
	}
	
	fprintf(file,"\n\t//DEPENDENT VARIABLES\n");
	cur = rewind_list(difvarlist);
	numEdos=0;
	while(cur != NULL)
	{
	    fprintf(file,"\nvoid __%s__(%s *classe){\n",cur->token.content,classname);
	    DiffList *curlDiff = rewind_list(difflist);
	    TokenNode *curDiff = NULL;
	    int cont=0;
	    while(curlDiff != NULL)
	    {
		if(cont==numEdos){
		    curDiff = curlDiff->diffheader->eq->next;
		    fprintf(file,"\t%s%s_lado_direito_= ",prefixo, curlDiff->diffheader->diffvar.content);
		    print_eq(file, curDiff);
		    fprintf(file,";\n");
		    //fprintf(file,"\t%s%s_new_= %sd%s*(%s%s_lado_direito_)+%s%s_old_;\t// %d\n",prefixo, curlDiff->diffheader->diffvar.content, prefixo, curlDiff->diffheader->freevar.content,prefixo, curlDiff->diffheader->diffvar.content,prefixo, curlDiff->diffheader->diffvar.content, curlDiff->number);
		    //fprintf(file,"\t//%s%s_old_ = %s%s_new_;\n",prefixo, curlDiff->diffheader->diffvar.content, prefixo, curlDiff->diffheader->diffvar.content);
		    curlDiff = curlDiff->next;
		    break;
		}else{
		    curDiff = curlDiff->diffheader->eq->next;
		    curlDiff = curlDiff->next;
		}
		cont++;
		
	    }
	    fprintf(file,"\n}\n");
	    cur = cur->next;
	    numEdos++;
	}
	prefixo="";
	
	//     fprintf(file,"typ_funcao edos[numEDO];\n");
	//     fprintf(file,"typ_funcao auxiliares[numAux];\n");
	
	//     fprintf(file,"typ_funcao equacoes[numAux];\n");
	fprintf(file,"\t#define forestSize %d\n", treeNumber);
	fprintf(file,"\ttyp_forest forest[forestSize];\n");
	for(int i=1;i<=treeNumber;i++){
	    int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
	    fprintf(file,"\t#define tree%dSize %d\n", i, treeSize);
	    fprintf(file,"\ttyp_funcao tree%d[tree%dSize];\n", i, i);
	}
	
	
	
	
	fprintf(file,"\n\t//Constructor, initializes all variables with 0.0 or default CellML initial values\n");
	fprintf(file, "\t%s::%s(double tol)\n\t{\n", classname,classname);
	
	for(int i=1;i<=treeNumber;i++){
	    int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
	    int countTreeItems = 0;
	    for(int j=0; j<totalSize; j++){
		if(grafoDependencias[j].treeId == i){
		    fprintf(file,"\t\ttree%d[%d].function = __%s__;\n", i, countTreeItems, grafoDependencias[j].equationName);
		    countTreeItems++;
		}
	    }
	    fprintf(file,"\t\tforest[%d].tree = tree%d;\n\n", i-1, i);
	    fprintf(file,"\t\tforest[%d].treeSize  = tree%dSize;\n\n", i-1, i);
	}
	
	cur = rewind_list(parvarlist); 		
		
	while(cur != NULL)
	{
	    fprintf(file,"\t\t%s = %.8e;\n", cur->token.content, cur->initialvalue);
	    cur = cur->next;
	}
	
	fprintf(file,"\t\td%s = 0.0; %s_vec__ = NULL;\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
	
	cur = rewind_list(difvarlist);
	while(cur != NULL)
	{
	    fprintf(file,"\t\t%s = NULL;\n",cur->token.content);
	    fprintf(file,"\t\t%s_ini_ = %.8e;\n",cur->token.content, cur->initialvalue);
	   
	    cur = cur->next;
	}
	fprintf(file,"\t\tit_countx = 0;\n");
	fprintf(file,"\t}\n");
	
	/****** FIM CONSTRUTOR *************************/
	
	fprintf(file, "\t%s::~%s()\n\t{\n", classname,classname);
	
	cur = rewind_list(difvarlist);
	while(cur != NULL){
	    fprintf(file, "\t\tif(%s != NULL) free(%s);\n", cur->token.content, cur->token.content);
	    cur = cur->next;
	    
	}
	fprintf(file, "\t}\n");		
	
	/****** METODOS SET ****************************/
	fprintf(file,"\n\tint %s::setVariables(int indVariable, double value_new)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
	
	cur = rewind_list(difvarlist);
	int numVar = 0;
	while(cur != NULL)
	{
	    fprintf(file,"\t\tcase %d:\t\t%s_ini_ = %s_old_= value_new;    break;\n", numVar, cur->token.content, cur->token.content);
	    numVar++;
	    cur = cur->next;
	}
	fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t\treturn 0;\n\t}\n");
	
	fprintf(file,"\n\tint %s::setParameters(int indVariable, double value_new)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
	
	cur = rewind_list(parvarlist);
	numVar = 0;
	while(cur != NULL)
	{
	    fprintf(file,"\t\tcase %d:\t\t%s = value_new;   break;\n", numVar, cur->token.content);
	    numVar++;
	    cur = cur->next;
	}
	fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t\treturn 0;\n\t}\n");
	
	
	fprintf(file,"\n\tint %s::setFreeVariable(double value_new)\n\t{\n\t\t",classname);
	fprintf(file,"d%s = value_new;\n\t}\n", difflist->diffheader->freevar.content);
	
	/// addt
	fprintf(file,"\n\tvoid %s::setMaxStep(double maxSt)\n\t{\n\t\t",classname);
	fprintf(file,"this->maxStep=maxSt;\n\t}\n");	
	
	fprintf(file,"\n\tvoid %s::executeTree(void* p, int treeSize){\n",classname);
	fprintf(file,"\t\ttyp_funcao* tree = (typ_funcao*)p;\n");
	fprintf(file,"\t\tfor(int k=0;k<treeSize;k++)\n");
	fprintf(file,"\t\t\ttree[k].funcao(this);\n");
	fprintf(file,"\t}\n");
	
	fprintf(file,"\n\tvoid %s::getRightHandSide(){\n",classname);
	fprintf(file,"\t\tfor(int k=0;k<forestSize;k++)\n");
	fprintf(file,"\t\t\tthis->executeTree((void*)forest[k].tree, forest[k].treeSize);\n");
	fprintf(file,"\t}\n");
	
	fprintf(file,"\n\tvoid %s::getRightHandSideParallel(){\n",classname);
	fprintf(file,"\t\t#pragma omp for schedule (dynamic)\n");
	fprintf(file,"\t\tfor(int k=0;k<forestSize;k++)\n");
	fprintf(file,"\t\t\tthis->executeTree((void*)forest[k].tree, forest[k].treeSize);\n");
	fprintf(file,"\t}\n");
	
	
	fprintf(file,"\n\tvoid %s::explicitEulerStep(){\n",classname);
	cur = rewind_list(difvarlist);
	while(cur != NULL)
	{
	    fprintf(file,"\t\tthis->%s_new_ = this->d%s*(this->%s_lado_direito_) + this->%s_old_;\n", 
		    cur->token.content, difflist->diffheader->freevar.content,cur->token.content,cur->token.content);
		    cur = cur->next;
	}
	fprintf(file,"\t}\n");
	
	fprintf(file,"\tint Solveode::getErrorCode(double error){\n");
	fprintf(file,"\t\tthis->errorAux = error;\n");
	fprintf(file,"\t\tif (error <0.5*_UPPER_LIMIT_*STEP_TOLERANCE_){\n");
	fprintf(file,"\t\t\t//Accept current solution, and increase the size of the next mechanics step.\n");
	fprintf(file,"\t\t\tcount4++;\n");
	fprintf(file,"\t\t\treturn 4;\n");
	fprintf(file,"\t\t}if (error>=0.5*_UPPER_LIMIT_*STEP_TOLERANCE_ && error<_UPPER_LIMIT_*STEP_TOLERANCE_){\n");
	fprintf(file,"\t\t\t//Accept current solution, and hold the step size fixed for the next mechanics step.\n");
	fprintf(file,"\t\t\tcount0++;\n");
	fprintf(file,"\t\t\treturn 0;\n");
	fprintf(file,"\t\t}if (error>=_UPPER_LIMIT_*STEP_TOLERANCE_ && error<2*_UPPER_LIMIT_*STEP_TOLERANCE_){\n");
	fprintf(file,"\t\t\t//Accept current solution, but decrease the size of the next mechanics step.\n");
	fprintf(file,"\t\t\tcount3++;\n");
	fprintf(file,"\t\t\treturn 3;\n");
	fprintf(file,"\t\t}else{//Throw current results away, cut the step size in half, and retry.\n");
	fprintf(file,"\t\t\tcount1++;\n");
	fprintf(file,"\t\t\treturn -1;\n");
	fprintf(file,"\t\t}\n");
	fprintf(file,"\t}\n");
	
	fprintf(file,"\t///runge kutta 2nd order's method\n");
	fprintf(file,"\t//Yn+1 = Yn +0.5h(Fn + F(Tn+H, Yn + FnH))\n");
	fprintf(file,"\t// double fn = equacao(t,y);\n");
	fprintf(file,"\t// double fnh = equacao(t+h, y+fn*h);\n");
	fprintf(file,"\t// return y + (fn + fnh)*h/2;\n");
	
	//Metodo runge kutta sem openmp
	fprintf(file,"\tvoid Solveode::rungekutta2_Single(){\n");
	fprintf(file,"\t\t//salva os valores das variaveis diferenciaveis em auxiliares\n");
	fprintf(file,"\t\tfor(int l=0;l<numEDO;l++){\n");
	fprintf(file,"\t\t\tladodir_fn_[l] = this->getLadoDireito(l);\n");
	fprintf(file,"\t\t\tedos_fn_[l] = this->getVariables(l);\n");
	fprintf(file,"\t\t\t//soma h às variaveis dt e todas as diferenciaveis, para encontrar Fn+1\n");
	fprintf(file,"\t\t\tthis->setVariables(l, edos_fn_[l] + this->d%s*ladodir_fn_[l]);//y+fn*h\n", 
		difflist->diffheader->freevar.content);
		fprintf(file,"\t\t}\n");
		fprintf(file,"\t\tthis->%s_new+= this->d%s;\n",difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
		fprintf(file,"\t\t//calcula os lados direitos de Fn+h\n");
		fprintf(file,"\t\tthis->getRightHandSide();\n");
		fprintf(file,"\t\t//volta o dt\n");
		fprintf(file,"\t\tthis->%s_new-= this->d%s;\n",difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
		fprintf(file,"\t\t// return y + (fn + fnh)*h/2;\n");
		cur = rewind_list(difvarlist);
		int countEDOS=0;
		while(cur != NULL)
		{
		    fprintf(file,"\t\tthis->%s_new_ = edos_fn_[%d] + (ladodir_fn_[%d]+this->%s_lado_direito_) * this->d%s/2;\n", 
			    cur->token.content,countEDOS,countEDOS, cur->token.content,difflist->diffheader->freevar.content);
			    cur = cur->next;
			    countEDOS++;
		}
		
		fprintf(file,"\t\tfor(int k=0;k<numEDO;k++)\tthis->setVariables(k, edos_fn_[k]);//volta com os verdadeiros _old_\n");
		fprintf(file,"\t}\n");
		
		
		//metodo runge kutta com openmp
		fprintf(file,"\tvoid Solveode::rungekutta2_Parallel(){\n");
		fprintf(file,"\t\t#pragma omp single\n");
		fprintf(file,"\t\t{\n");
		fprintf(file,"\t\t\t//salva os valores das variaveis diferenciaveis em auxiliares\n");
		fprintf(file,"\t\t\tfor(int l=0;l<numEDO;l++){\n");
		fprintf(file,"\t\t\t\tladodir_fn_[l] = this->getLadoDireito(l);\n");
		fprintf(file,"\t\t\t\tedos_fn_[l] = this->getVariables(l);\n");
		fprintf(file,"\t\t\t\t//soma h às variaveis dt e todas as diferenciaveis, para encontrar Fn+1\n");
		fprintf(file,"\t\t\t\tthis->setVariables(l, edos_fn_[l] + this->d%s*ladodir_fn_[l]);//y+fn*h\n", 
			difflist->diffheader->freevar.content);
			fprintf(file,"\t\t\t}\n");
			fprintf(file,"\t\t\tthis->%s_new+= this->d%s;\n",difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
			fprintf(file,"\t\t}\n");
			
			fprintf(file,"\t\t//calcula os lados direitos de Fn+h\n");
			fprintf(file,"\t\tthis->getRightHandSideParallel();\n");
			
			fprintf(file,"\t\t#pragma omp single\n");
			fprintf(file,"\t\t{\n");
			fprintf(file,"\t\t\t//volta o dt\n");
			fprintf(file,"\t\t\tthis->%s_new-= this->d%s;\n",difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
			fprintf(file,"\t\t\t// return y + (fn + fnh)*h/2;\n");
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\tthis->%s_new_ = edos_fn_[%d] + (ladodir_fn_[%d]+this->%s_lado_direito_) * this->d%s/2;\n", 
				cur->token.content,countEDOS,countEDOS, cur->token.content,difflist->diffheader->freevar.content);
				cur = cur->next;
				countEDOS++;
			}
			
			fprintf(file,"\t\t\tfor(int k=0;k<numEDO;k++)\tthis->setVariables(k, edos_fn_[k]);//volta com os verdadeiros _old_\n");
			fprintf(file,"\t\t}\n");
			fprintf(file,"\t}\n");
			
			//método adaptativo 
			fprintf(file,"\tint Solveode::adaptiveDt_Single(){\n");
			fprintf(file,"\t\t//calcula os lados direitos\n");
			fprintf(file,"\t\tthis->getRightHandSide();\n");
			fprintf(file,"\t\t//calcula uma iteração pelo metodo de euler\n");
			fprintf(file,"\t\tthis->explicitEulerStep();\n");
			fprintf(file,"\t\t//armazena o resultado obtido por euler em vetores.\n");
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\tedos_euler_[%d] = this->%s_new_;\n",
				countEDOS, cur->token.content);
				cur = cur->next;
				countEDOS++;
			}
			
			fprintf(file,"\t\t//calcula as edos pelo metodo trapezoidal, para o mesmo passo\n");
			fprintf(file,"\t\tthis->rungekutta2_Single();\n");
			
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\tedos_rk2_[%d] = this->%s_new_;\n", 
				countEDOS, cur->token.content);
				cur = cur->next;
				countEDOS++;
			}
			fprintf(file,"\t\t//calcula as diferencas entre cada metodo\n");
			fprintf(file,"\t\tif(edos_rk2_[0]!=0)\n");
			fprintf(file,"\t\t\terrorAD_ = fabs((edos_rk2_[0]-edos_euler_[0])/edos_rk2_[0]) ;//calcula o erro relativo para iniciar\n");
			fprintf(file,"\t\tdouble auxError=0;\n");
			fprintf(file,"\t\tfor(int k=0;k<numEDO;k++){\n");
			fprintf(file,"\t\t\tif(edos_rk2_[k]!=0){\n");
			fprintf(file,"\t\t\t\tauxError = fabs((edos_rk2_[k]-edos_euler_[k])/edos_rk2_[k]) ;//calcula o erro relativo\n");
			fprintf(file,"\t\t\t\tif(auxError>errorAD_)\n");
			fprintf(file,"\t\t\t\t\terrorAD_ = auxError;\n");
			fprintf(file,"\t\t\t}\n");
			fprintf(file,"\t\t}\n");
			fprintf(file,"\t\tif(this->dtime==0.0){\n");
			fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - adaptiveDt - dtime reach the limit (0.0)\\n\");\n");
			fprintf(file,"\t\t\texit(1);\n");
			fprintf(file,"\t\t}\n");    
			fprintf(file,"\t\treturn this->getErrorCode(errorAD_);\n");
			fprintf(file,"\t}\n");
			
			//metodo adaptativo paralelo
			fprintf(file,"\tvoid Solveode::adaptiveDt_Parallel(){\n");
			
			fprintf(file,"\t\t//calcula os lados direitos\n");
			fprintf(file,"\t\tthis->getRightHandSideParallel();\n");
			fprintf(file,"\t\t#pragma omp for schedule(static)\n");
			fprintf(file,"\t\tfor(int l=0;l<numEDO;l++){\n");
			fprintf(file,"\t\t\tladodir_fn_[l] = this->getLadoDireito(l);\n");
			fprintf(file,"\t\t\tedos_fn_[l] = this->getVariables(l);//salva os valores das variaveis diferenciaveis em auxiliares\n");
			fprintf(file,"\t\t\tedos_euler_[l] = edos_fn_[l] + this->dtime*ladodir_fn_[l];\n");
			fprintf(file,"\t\t\tthis->setVariables(l, edos_euler_[l]);//y+fn*h\n");
			fprintf(file,"\t\t}\n");
			
			fprintf(file,"\t\t#pragma omp single\n");
			fprintf(file,"\t\t{\n");
			fprintf(file,"\t\t\tthis->time_new+= this->dtime;\n");
			fprintf(file,"\t\t}\n");
			fprintf(file,"\t\t//calcula os lados direitos de Fn+h\n");
			fprintf(file,"\t\tthis->getRightHandSideParallel();\n");
			fprintf(file,"\t\t#pragma omp single\n");
			fprintf(file,"\t\t{\n");
			fprintf(file,"\t\t\t//volta o dt\n");
			fprintf(file,"\t\t\tthis->time_new-= this->dtime;\n");
			fprintf(file,"\t\t\t// return y + (fn + fnh)*h/2;\n");
			fprintf(file,"\t\t\t// 	    this->V_new_ = edos_fn_[0] + (ladodir_fn_[0]+this->V_lado_direito_) * this->dtime/2;\n");
			fprintf(file,"\t\t}\n");
			
			fprintf(file,"\t\t#pragma omp for schedule(static)\n");
			fprintf(file,"\t\tfor(int k=0;k<numEDO;k++){\n");
			fprintf(file,"\t\t\tedos_rk2_[k] = edos_fn_[k] + (ladodir_fn_[k]+this->getLadoDireito(k)) * this->dtime/2;\n");
			fprintf(file,"\t\t\tthis->setVariables(k, edos_fn_[k]);//volta com os verdadeiros _old_\n");
			fprintf(file,"\t\t}\n");
			fprintf(file,"\t}\n");
			//fim
			
			
			fprintf(file,"\tdouble Solveode::solveToFileOMP(double finalTime, double savingRate, int method, char *fileName, int nThreads, double tol)\n");
			fprintf(file,"\t{\n");
			fprintf(file,"\t\tthis->STEP_TOLERANCE_ = tol;\n");
			fprintf(file,"\t\tFILE *fileptr;\n");
			fprintf(file,"\t\tif(savingRate!=0.0)\n");
			fprintf(file,"\t\t\tfileptr = fopen(fileName, \"w+\");\n");
			
			fprintf(file,"\t\t%s_new = %s;\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
			fprintf(file,"\t\tdouble timeSaving = %s_new;\n",difflist->diffheader->freevar.content);
			fprintf(file,"\t\tdouble previousDtime = this->d%s;\n",difflist->diffheader->freevar.content);
			
			fprintf(file,"\t\tif(finalTime <= 0){\n");
			fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - solveToFile - negative finalTime is not allowed\\n\");\n");
			fprintf(file,"\t\t\texit(1);\n");
			fprintf(file,"\t\t}\n");
			fprintf(file,"\t\tif(savingRate < 0){\n");
			fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - solveToFile - negative saving rate is not allowed\\n\");\n");
			fprintf(file,"\t\t\texit(1);\n");
			fprintf(file,"\t\t}\n");
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\t%s_old_ = %s_ini_;\n", 
				cur->token.content,cur->token.content);
				cur = cur->next;
				countEDOS++;
			}
			fprintf(file,"\t\tint counter_it__ = 1, it_countx = 1;\n");
			fprintf(file,"\t\tdouble countSaving=0;\n");
			fprintf(file,"\t\tdouble sumError=0;\n");
			fprintf(file,"\t\tdouble auxDtime = dtime;\n");
			fprintf(file,"\t\tint flag=0;\n");
			fprintf(file,"\t\tint countRe = 0, cont2=0, cont0=0, cont3=0, cont4=0;\n");
			
			fprintf(file,"\t\tdouble maxDt=dtime, minDt=dtime, maxError=0, minError=9999999;\n");
			fprintf(file,"\t\tomp_set_num_threads(nThreads);\n");
			fprintf(file,"\t\t#pragma omp parallel\n");
			fprintf(file,"\t\twhile(time_new<finalTime){\n");
			fprintf(file,"\t\t\tthis->recursion_counter = 0 ;\n");
			fprintf(file,"\t\t\tswitch(method){\n");
			fprintf(file,"\t\t\t\tcase _EULER_:\n");
			
			fprintf(file,"\t\t\t\t\tthis->getRightHandSideParallel();\n");
			fprintf(file,"\t\t\t\t\t#pragma omp for schedule(static)\n");
			fprintf(file,"\t\t\t\t\tfor(int l=0;l<numEDO;l++){\n");
			fprintf(file,"\t\t\t\t\t\tedos_euler_[l] = this->getVariables(l) + this->dtime*this->getLadoDireito(l);\n");
			fprintf(file,"\t\t\t\t\t}\n");
			fprintf(file,"\t\t\t\t\t#pragma omp single\n");
			fprintf(file,"\t\t\t\t\t{\n");
			    
			//this->V_new_ = edos_euler_[0];
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t\t\t\tthis->%s_new_ = edos_euler_[%d];\n", 
				cur->token.content, countEDOS);
				cur = cur->next;
			    countEDOS++;
			}
			fprintf(file,"\t\t\t\t}\n");
			fprintf(file,"\t\t\t\tbreak;\n");
			
			
			fprintf(file,"\t\t\t\tcase _ADAP_DT_:\n");
			fprintf(file,"\t\t\t\t\tdo{\n");
			
			fprintf(file,"\t\t\t\t\t\tthis->adaptiveDt_Parallel();\n");
			fprintf(file,"\t\t\t\t\t\t#pragma omp single\n");
			fprintf(file,"\t\t\t\t\t\t{\n");
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
// 			    this->V_new_ = edos_rk2_[0];
			    fprintf(file,"\t\t\t\t\t\t\tthis->%s_new_ = edos_rk2_[%d];\n", 
				    cur->token.content, countEDOS);
				    cur = cur->next;
			    countEDOS++;
			}
			fprintf(file,"\t\t\t\t\t\t\t//calcula as diferencas entre cada metodo\n");
			fprintf(file,"\t\t\t\t\t\t\tif(edos_rk2_[0]!=0)\n");
			fprintf(file,"\t\t\t\t\t\t\t\terrorAD_ = fabs((edos_rk2_[0]-edos_euler_[0])/edos_rk2_[0]) ;//calcula o erro relativo para iniciar\n");
			fprintf(file,"\t\t\t\t\t\t\tdouble auxError=0;\n");
			fprintf(file,"\t\t\t\t\t\t\tfor(int k=0;k<numEDO;k++){\n");
			fprintf(file,"\t\t\t\t\t\t\t\tif(edos_rk2_[k]!=0){\n");
			fprintf(file,"\t\t\t\t\t\t\t\t\tauxError = fabs((edos_rk2_[k]-edos_euler_[k])/edos_rk2_[k]) ;//calcula o erro relativo\n");
			fprintf(file,"\t\t\t\t\t\t\t\t\tif(auxError>errorAD_)\n");
			fprintf(file,"\t\t\t\t\t\t\t\t\terrorAD_ = auxError;\n");
			fprintf(file,"\t\t\t\t\t\t\t\t}\n");
			fprintf(file,"\t\t\t\t\t\t\t}\n");
			fprintf(file,"\t\t\t\t\t\t\tif(this->dtime==0.0){\n");
			fprintf(file,"\t\t\t\t\t\t\t\tfprintf(stderr,\"ERROR - adaptiveDt - dtime reach the limit (0.0)\\n\");\n");
			fprintf(file,"\t\t\t\t\t\t\t\texit(1);\n");
			fprintf(file,"\t\t\t\t\t\t\t}\n");
				
			fprintf(file,"\t\t\t\t\t\t\tflag = this->getErrorCode(errorAD_);\n");
			fprintf(file,"\t\t\t\t\t\t\tif(flag==-1){\n");
			fprintf(file,"\t\t\t\t\t\t\t\tthis->dtime = this->dtime/_H_FACTOR_;\n");
			fprintf(file,"\t\t\t\t\t\t\t\tcountRe++;\n");
			fprintf(file,"\t\t\t\t\t\t\t}\n");
			fprintf(file,"\t\t\t\t\t\t}\n");
			fprintf(file,"\t\t\t\t\t}while(flag==-1);\n");
			fprintf(file,"\t\t\t\t\tbreak;\n");
			
			fprintf(file,"\t\t\t\tcase _RK2_:\n");
			fprintf(file,"\t\t\t\t\tthis->getRightHandSideParallel();\n");
			fprintf(file,"\t\t\t\t\tthis->rungekutta2_Parallel();\n");
			fprintf(file,"\t\t\t\tbreak;\n");
			fprintf(file,"\t\t\t\tdefault:\n");
			fprintf(file,"\t\t\t\t\tprintf(\"Invalid option. Choose 1 for Explicit Euler and 3 for Adaptive step size method.\\n\");\n");
			fprintf(file,"\t\t\t\t\texit(1);\n");
			fprintf(file,"\t\t\t}\n");
			
			fprintf(file,"\t\t\t#pragma omp single\n");
			fprintf(file,"\t\t\t{\n");    
			fprintf(file,"\t\t\t\tif(savingRate!=0.0){\n");
			fprintf(file,"\t\t\t\t\tdouble diff = _agos_round(this->time_new,32) - _agos_round(timeSaving,32);\n");
			fprintf(file,"\t\t\t\t\tif(diff==0){\n");
			fprintf(file,"\t\t\t\t\t\tfprintf(fileptr,\"%%.5e %%.10e %%.4e %%.10e %%.10e %%.10e %%.10e\\n\",\n");
			fprintf(file,"\t\t\t\t\t\tthis->time_new, this->V_new_, auxDtime*45000,this->errorAux,STEP_TOLERANCE_/2,  STEP_TOLERANCE_, 2*STEP_TOLERANCE_);\n");
			fprintf(file,"\t\t\t\t\t\tcounter_it__++;\n");
			fprintf(file,"\t\t\t\t\t\ttimeSaving = timeSaving + savingRate;\n");
			fprintf(file,"\t\t\t\t\t}else if(diff > 0){\n");
			fprintf(file,"\t\t\t\t\t\t//salva resultados que passaram do tempo de salvar\n");
			fprintf(file,"\t\t\t\t\t\tauxDtime = this->dtime;\n");
			fprintf(file,"\t\t\t\t\t\tdouble _temp_[numEDO];\n");
			fprintf(file,"\t\t\t\t\t\tdouble auxTimenew = this->time_new;\n");
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t\t\t\t_temp_[%d] = this->%s_new_;\n", 
				countEDOS, cur->token.content);
				cur = cur->next;
				countEDOS++;
			}
			
			fprintf(file,"\t\t\t\t\t\tdouble told = this->time_new - previousDtime;\n");
			fprintf(file,"\t\t\t\t\t\t//encontra-se o dt adequado pra salvar na hora certa\n");
			fprintf(file,"\t\t\t\t\t\tdouble dtTemp =   timeSaving - (told);\n");
			
			fprintf(file,"\t\t\t\t\t\tthis->dtime = dtTemp;\n");
			fprintf(file,"\t\t\t\t\t\tthis->time_new = timeSaving - dtTemp;\n");
			
			fprintf(file,"\t\t\t\t\t\t//verifica se a iteração acima alcançou o tempo real\n");
			fprintf(file,"\t\t\t\t\t\twhile(auxTimenew >= timeSaving){\n");
			fprintf(file,"\t\t\t\t\t\t\tthis->time_new += this->dtime;\n");
			fprintf(file,"\t\t\t\t\t\t\t//calcula mais uma iteração\n");
			fprintf(file,"\t\t\t\t\t\t\tif( method ==_ADAP_DT_ ){\n");
			fprintf(file,"\t\t\t\t\t\t\t\tthis->adaptiveDt_Single();\n");
			fprintf(file,"\t\t\t\t\t\t\t}else if( method == _EULER_ ){\n");
			fprintf(file,"\t\t\t\t\t\t\t\tthis->getRightHandSide();\n");
			fprintf(file,"\t\t\t\t\t\t\t\tthis->explicitEulerStep();\n");
			fprintf(file,"\t\t\t\t\t\t\t}else if( method = _RK2_ ){\n");
			fprintf(file,"\t\t\t\t\t\t\t\tthis->getRightHandSide();\n");
			fprintf(file,"\t\t\t\t\t\t\t\tthis->rungekutta2_Single();\n");
			fprintf(file,"\t\t\t\t\t\t\t}\n");
			fprintf(file,"\t\t\t\t\t\t\tfprintf(fileptr,\"%%.5e %%.10e %%.4e %%.10e %%.10e %%.10e %%.10e\\n\",\n");
			fprintf(file,"\t\t\t\t\t\t\tthis->time_new, this->V_new_, auxDtime*45000,this->errorAux,STEP_TOLERANCE_/2,  STEP_TOLERANCE_, 2*STEP_TOLERANCE_);\n");
			fprintf(file,"\t\t\t\t\t\t\ttimeSaving += savingRate;\n");
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t\t\t\t\tthis->%s_old_ = this->%s_new_;\n", 
				cur->token.content, cur->token.content);
				cur = cur->next;
				countEDOS++;
			}
			fprintf(file,"\t\t\t\t\t\t\t//coloca o dtime igual ao savingRate para bater exatamente com a proxima iteração de salvar\n");
			fprintf(file,"\t\t\t\t\t\t\tthis->dtime = savingRate;\n");
			fprintf(file,"\t\t\t\t\t\t}//fimwhile\n");
			fprintf(file,"\t\t\t\t\t\t//recupera as variaveis antigas, para nao ter que recalcular o passo\n");
			
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t\t\t\tthis->%s_new_ = _temp_[%d];\n", cur->token.content, countEDOS);
			    cur = cur->next;
			    countEDOS++;
			}
			
			
			fprintf(file,"\t\t\t\t\t\tthis->time_new = auxTimenew;\n");
			fprintf(file,"\t\t\t\t\t\tthis->dtime  = auxDtime;\n");
			fprintf(file,"\t\t\t\t\t\tcounter_it__++;\n");
			
			fprintf(file,"\t\t\t\t\t}//fim else if\n");
			fprintf(file,"\t\t\t\t}//fim savingrate=0\n");
			
			fprintf(file,"\t\t\t\tthis->%s_new += this->d%s;\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			fprintf(file,"\t\t\t\tpreviousDtime = this->d%s;\n",difflist->diffheader->freevar.content);
			
			cur = rewind_list(difvarlist);
			countEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t\tthis->%s_old_ = this->%s_new_;\n", 
				cur->token.content, cur->token.content);
				cur = cur->next;
				countEDOS++;
			}
			fprintf(file,"\t\t\t\tit_countx++;\n");
			fprintf(file,"\t\t\t\tif(flag==2){\n");
			fprintf(file,"\t\t\t\t\tthis->dtime = _H_FACTOR_*this->dtime;\n");
			fprintf(file,"\t\t\t\t\tcont2++;\n");
			fprintf(file,"\t\t\t\t}else if(flag==3){\n");
			fprintf(file,"\t\t\t\t\tthis->dtime = this->dtime*0.65;\n");
			fprintf(file,"\t\t\t\t\tcont3++;\n");
			fprintf(file,"\t\t\t\t}else if(flag==4){\n");
			fprintf(file,"\t\t\t\t\tthis->dtime = this->dtime*1.5;\n");
			fprintf(file,"\t\t\t\t\tcont4++;\n");
			fprintf(file,"\t\t\t\t}else if(flag==0){\n");
			fprintf(file,"\t\t\t\t\tcont0++;\n");
			fprintf(file,"\t\t\t\t}else{printf(\"flag: %%d\\n\", flag);}\n");
			
			fprintf(file,"\t\t\t\tif(this->dtime > this->maxStep && this->maxStep!=0){\n");
			fprintf(file,"\t\t\t\t\tthis->dtime = this->maxStep;\n");
			fprintf(file,"\t\t\t\t}\n");
			fprintf(file,"\t\t\t}//fim pragma omp single\n");
			fprintf(file,"\t\t} //fim while\n");
			fprintf(file,"\t\tif(savingRate!=0.0)\n");
			fprintf(file,"\t\t\tfclose(fileptr);\n");
			fprintf(file,"\t\tprintf(\"%%d\\n\",it_countx);\n");
			fprintf(file,"\t\t////printf(\"%%d %%d %%d %%d\\n\", count1, count0, count3, count4);\n");
			fprintf(file,"\t\t////printf(\"%%d %%d %%d %%d\\n\", countRe,  cont0,  cont3,  cont4);\n");
			
			fprintf(file,"\t\treturn it_countx;\n");
			fprintf(file,"\t}\n");
			/******* FIM SET ***********************/
			
			/******* METODOS GET ********************/
			fprintf(file,"\n\t//Get Methods\n");
			
			fprintf(file,"\n\tdouble %s::getVariables(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
			
			cur = rewind_list(difvarlist);
			numVar = 0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\tcase %d:\t\treturn %s_old_;    break;\n", numVar, cur->token.content);
			    numVar++;
			    cur = cur->next;
			}
			fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t}\n");
			
			
			//////
			fprintf(file,"\n\tdouble %s::getLadoDireito(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
			
			cur = rewind_list(difvarlist);
			numVar = 0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\tcase %d:\t\treturn %s_lado_direito_;    break;\n", numVar, cur->token.content);
			    numVar++;
			    cur = cur->next;
			}
			fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t}\n");
			
			fprintf(file,"\n\tdouble %s::getParameters(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
			cur = rewind_list(parvarlist);
			numVar = 0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\tcase %d:\t\treturn %s;    break;\n", numVar, cur->token.content);
			    numVar++;
			    cur = cur->next;
			}
			fprintf(file,"\t\tdefault:\tbreak;\n\t\t}\n\t}\n");
			
			fprintf(file,"\n\tdouble %s::getFreeVariable()\n\t{\n",classname);
			
			fprintf(file,"\t\treturn d%s;\n\t}\n", difflist->diffheader->freevar.content);
			
			
			fprintf(file,"\n\t//Get Methods - Variables\n\n");
			fprintf(file, "\tVariables %s::get_Variables()\n\t{\n\t\tVariables v(\"",classname);
			
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"|%s#",cur->token.content);
			    cur = cur->next;
			}
			fprintf(file, "\");\n\t\treturn v;\n\t}\n");
			
			fprintf(file, "\tVariables %s::get_Parameters()\n\t{\n\t\tVariables v(\"",classname);
			
			cur = rewind_list(parvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"|%s#",cur->token.content);
			    cur = cur->next;
			}
			fprintf(file, "\");\n\t\treturn v;\n\t}\n");	
			
			fprintf(file, "\tVariables %s::get_FreeVariable()\n\t{\n\t\tVariables v(\"",classname);
			
			fprintf(file,"|%s#",difflist->diffheader->freevar.content); 		
			fprintf(file, "\");\n\t\treturn v;\n\t}\n");	
			
			/********SET DO ARQUIVO **********************/
			fprintf(file,"\n\tvoid %s::setParametersFromFile(char *filename)\n\t{\n",classname);
			fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
			fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setParametersFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
			fprintf(file,"\tdouble value;\n\t\tint k = 0;\n\t\tVariables v = get_Parameters();\n\t\tint s = v.getQuantity();\n\t");
			fprintf(file,"\tfor(;k<s;k++)\n\t\t{\n\t\t");
			fprintf(file,"\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetParameters(k, value);\n\t\t}\n\t\tfclose(file);\n\t}\n");
			
			fprintf(file,"\n\tvoid %s::setVariablesFromFile(char *filename)\n\t{\n", classname);
			fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
			fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setVariablesFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
			fprintf(file,"\tdouble value;\n\t\tint k = 0;\n\t\tVariables v = get_Variables();\n\t\tint s = v.getQuantity();\n\t");
			fprintf(file,"\tfor(;k<s;k++)\n\t\t{\n\t\t");
			fprintf(file,"\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetVariables(k, value);\n\t\t}\n\t\tfclose(file);\n\t}\n");
			
			fprintf(file,"\n\tvoid %s::setFreeVariableFromFile(char *filename)\n\t{\n", classname);
			fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
			fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setFreeVariableFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
			fprintf(file,"\tdouble value;\n\t\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetFreeVariable(value);\n\t\tfclose(file);\n\t}\n");
			/********FIM SET DO ARQUIVO ******************/
			
			/******* METODO SOLVE ******************/
			fprintf(file,"\tdouble %s::solve(int firstcall__, int num_iterations__, int num_results__, int numThreads)\n\t{\n\n",classname);		
			fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
			fprintf(file,"\t\tstatic int num_results_bak = 0;\n");
			fprintf(file,"\t\tstatic int offset_step = 1;\n");
			fprintf(file,"\t\tif(firstcall__){\n");
			fprintf(file,"\t\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			
			fprintf(file,"\t\t\tif(num_results__ <= 0)\n\t\t\t\tnum_results__ = 1;\n\t\t\tif(num_iterations__ <= 0)\n\t\t\t\tnum_iterations__ = 1;\n");
			fprintf(file,"\t\t\toffset_step = num_iterations__ / num_results__;\n");
			
			fprintf(file,"\t\t\tif(%s_vec__ != NULL)free( %s_vec__);\n\t\t\t%s_vec__ = (double *)malloc(sizeof(double)*num_results__);\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
			    fprintf(file,"\t\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*num_results__);\n", cur->token.content,cur->token.content,cur->token.content);
			    cur = cur->next;
			}
			fprintf(file,"\t\t\tnum_results_bak = num_results__;\n");
			fprintf(file,"\t\t\tnum_iterations_bak = num_iterations__;\n");
			
			fprintf(file,"\t\t}\n");
			
			fprintf(file,"\t\tint counter_it__ = 1, it_countx = 1;\n");
			fprintf(file,"\t\tint aux = num_iterations_bak%%num_results_bak;\n");
			//fprintf(file,"\t\tomp_set_num_threads(numThreads);\n");
			//fprintf(file,"\t\t#pragma omp parallel for ordered schedule(static)\n");
			
			fprintf(file,"\t\tfor(it_countx = 1; it_countx<=num_iterations_bak; it_countx++){\n");
			fprintf(file,"\t\t\t%s_new += d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			
			fprintf(file,"\t\t\t////for(int k=0;k<numAux;k++)\t");
			fprintf(file,"////auxiliares[k].funcao(this);\n");
			
			fprintf(file,"\t\t\t////for(int k=0;k<numEDO;k++)\t");
			fprintf(file,"////edos[k].funcao(this);\n");
			fprintf(file,"\t\t\tthis->getRightHandSide();\n");
			
			//fprintf(file,"\n\t\t\t#pragma omp critical\n");
			fprintf(file,"\t\t\tif(it_countx != aux && (it_countx-aux)%%offset_step == 0){\n");
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t\t%s[counter_it__] = %s_new_;\n", cur->token.content,cur->token.content);
			    cur = cur->next;
			}
			fprintf(file, "\t\t\t\t%s_vec__[counter_it__] = %s_new;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			fprintf(file,"\t\t\t\tcounter_it__++;\n");
			fprintf(file,"\t\t\t}\n");
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t%s_old_ = %s_new_;\n", cur->token.content,cur->token.content);
			    cur = cur->next;
			}
			fprintf(file,"\t\t}\n");
			
			fprintf(file,"\t\t\t//FILE *fileptr;\n");
			fprintf(file,"\t\t\t//char filename[12];\n");
			fprintf(file,"\t\t\t//sprintf(filename,\"%%dthread.dat\",numThreads);\n");
			fprintf(file,"\t\t\t//fileptr = fopen(filename, \"wb\");\n");
			fprintf(file,"\t\t\t//for(int i =0;i<num_results__;i++){\n");
			fprintf(file,"\t\t\t//    fprintf(fileptr,\"%%f %%f\\n\", time_vec__[i], V[i]);\n");
			fprintf(file,"\t\t\t//}\n");
			fprintf(file,"\t\t\t//fclose(fileptr);\n");
			fprintf(file, "\t\treturn (num_iterations_bak%%offset_step)*d%s;\n", difflist->diffheader->freevar.content);
			fprintf(file,"\t}\n");
			/****** FIM SOLUCAO ******************/
			
			
			
			
			/******* METODO JACOBIAN ******************/
			fprintf(file,"\tdouble %s::jacobian(int firstcall__, int num_iterations__, int num_results__, int numThreads)\n\t{\n\n",classname);
			fprintf(file,"\t\tthis->solve(firstcall__, num_iterations__, num_results__, numThreads);\n");
			fprintf(file,"\t\tdouble h_jac_[numEDO];\n");
			fprintf(file,"\t\tdouble quociente = 1000.0;\n");
			cur = rewind_list(difvarlist);
			int contaEDOS=0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\th_jac_[%d] = abs(_agos_max(%s, num_results__) - _agos_min(%s, num_results__)) / quociente;\n", contaEDOS, cur->token.content,cur->token.content, cur->token.content,cur->token.content);
			    cur = cur->next;
			    contaEDOS++;
			}
			
			
			fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
			fprintf(file,"\t\tstatic int num_results_bak = 0;\n");
			fprintf(file,"\t\tstatic int offset_step = 1;\n\n");
			
			fprintf(file,"\t\tdouble jacobian[numEDO][numEDO];\n");
			fprintf(file,"\t\tdouble edos_new_aux_[numEDO];\n");
			fprintf(file,"\t\tdouble edos_aux_[numEDO];\n");
			fprintf(file,"\t\tFILE *filejac;\n");
			fprintf(file,"\t\tfilejac = fopen(\"dat/jacobian.dat\", \"wb\");\n");
			
			fprintf(file,"\t\tif(firstcall__){\n");
			fprintf(file,"\t\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			
			fprintf(file,"\t\t\tif(num_results__ <= 0)\n\t\t\t\tnum_results__ = 1;\n\t\t\tif(num_iterations__ <= 0)\n\t\t\t\tnum_iterations__ = 1;\n");
			fprintf(file,"\t\t\toffset_step = num_iterations__ / num_results__;\n");
			
			fprintf(file,"\t\t\tif(%s_vec__ != NULL)free( %s_vec__);\n\t\t\t%s_vec__ = (double *)malloc(sizeof(double)*num_results__);\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
			    fprintf(file,"\t\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*num_results__);\n", cur->token.content,cur->token.content,cur->token.content);
			    cur = cur->next;
			}
			fprintf(file,"\t\t\tnum_results_bak = num_results__;\n");
			fprintf(file,"\t\t\tnum_iterations_bak = num_iterations__;\n");
			
			fprintf(file,"\t\t}\n");
			
			fprintf(file,"\t\tint counter_it__ = 1, it_countx = 1;\n");
			fprintf(file,"\t\tint aux = num_iterations_bak%%num_results_bak;\n");
			//fprintf(file,"\t\tomp_set_num_threads(numThreads);\n");
			fprintf(file,"\t\tfor(it_countx = 1; it_countx<=num_iterations_bak; it_countx++){\n");
			fprintf(file,"\t\t\t%s_new += d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			
			fprintf(file,"\t\t\t////for(int k=0;k<numAux;k++)\t");
			fprintf(file,"////auxiliares[k].funcao(this);\n");
			
			fprintf(file,"\t\t\t////for(int k=0;k<numEDO;k++)\t");
			fprintf(file,"////edos[k].funcao(this);\n");
			fprintf(file,"\t\t\tthis->getRightHandSide();\n");
			fprintf(file,"\t\t\tif(it_countx != aux && (it_countx-aux)%%offset_step == 0){\n");
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t\t%s[counter_it__] = %s_new_;\n", cur->token.content,cur->token.content);
			    cur = cur->next;
			}
			fprintf(file, "\t\t\t\t%s_vec__[counter_it__] = %s_new;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			fprintf(file,"\t\t\t\tcounter_it__++;\n");
			fprintf(file,"\t\t\t\t//salva os valores das variaveis diferenciaveis em auxiliares\n");
			fprintf(file,"\t\t\t\tfor(int l=0;l<numEDO;l++){\n");
			fprintf(file,"\t\t\t\t\tedos_aux_[l] = this->getLadoDireito(l);\n");
			fprintf(file,"\t\t\t\t\tedos_new_aux_[l] = this->getVariables(l);\n");
			fprintf(file,"\t\t\t\t}\n");
			
			fprintf(file,"\t\t\t\t//para cada coluna\n");
			fprintf(file,"\t\t\t\tfor(int k=0;k<numEDO;k++){\n");
			fprintf(file,"\t\t\t\t\t//escolhe uma variavel k, e calcula a derivada de todas as equacoes em relacao a k\n");
			fprintf(file,"\t\t\t\t\tthis->setVariables(k, edos_new_aux_[k]+h_jac_[k]);\n");
			fprintf(file,"\t\t\t\t\t//calcula tudo de novo com o novo valor de k\n");
			fprintf(file,"\t\t\t\t\t////for(int i=0;i<numAux;i++)\tauxiliares[i].funcao(this);\n");
			fprintf(file,"\t\t\t\t\t////for(int i=0;i<numEDO;i++)\tedos[i].funcao(this);\n");
			fprintf(file,"\t\t\t\t\tthis->getRightHandSide();\n");
			fprintf(file,"\t\t\t\t\tfor(int j=0;j<numEDO;j++){//para cada linha \n");
			fprintf(file,"\t\t\t\t\t\tjacobian[j][k] = (this->getLadoDireito(j) - edos_aux_[j])/h_jac_[k];\n");
			fprintf(file,"\t\t\t\t\t\tfprintf(filejac,\"%%f\\t\", jacobian[j][k]);\n");
			fprintf(file,"\t\t\t\t\t}\n");
			fprintf(file,"\t\t\t\t\tfprintf(filejac,\"\\n\");\n");
			fprintf(file,"\t\t\t\t\t//agora tem que voltar para o que estava antes, sem somar dtime a variavel k\n");
			fprintf(file,"\t\t\t\t\tfor(int l=0;l<numEDO;l++){\n");
			fprintf(file,"\t\t\t\t\t\tthis->setVariables(l, edos_new_aux_[l]);\n");
			fprintf(file,"\t\t\t\t\t}\n");
			fprintf(file,"\t\t\t\t\t////for(int i=0;i<numAux;i++)\tauxiliares[i].funcao(this);\n");
			fprintf(file,"\t\t\t\t\t////for(int i=0;i<numEDO;i++)\tedos[i].funcao(this);\n");
			fprintf(file,"\t\t\t\t\tthis->getRightHandSide();\n");
			fprintf(file,"\t\t\t\t}\n");
			fprintf(file,"\t\t\t\tfprintf(filejac,\"\\n\");\n");
			fprintf(file,"\t\t\t}\n");
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t%s_old_ = %s_new_;\n", cur->token.content,cur->token.content);
			    cur = cur->next;
			}
			fprintf(file,"\t\t}\n");
			fprintf(file,"\t\tfclose(filejac);\n");
			fprintf(file,"\t\t\t//FILE *fileptr;\n");
			fprintf(file,"\t\t\t//char filename[12];\n");
			fprintf(file,"\t\t\t//sprintf(filename,\"%%dthread.dat\",numThreads);\n");
			fprintf(file,"\t\t\t//fileptr = fopen(filename, \"wb\");\n");
			fprintf(file,"\t\t\t//for(int i =0;i<num_results__;i++){\n");
			fprintf(file,"\t\t\t//    fprintf(fileptr,\"%%f %%f\\n\", time_vec__[i], V[i]);\n");
			fprintf(file,"\t\t\t//}\n");
			fprintf(file,"\t\t\t//fclose(fileptr);\n");
			fprintf(file, "\t\treturn (num_iterations_bak%%offset_step)*d%s;\n", difflist->diffheader->freevar.content);
			fprintf(file,"\t}\n");
			/****** FIM JACOBIAN ******************/
			
			fprintf(file, "\n\tdouble* %s::getIndependentVar()\n\t{\n", classname);
			fprintf(file, "\t\treturn %s_vec__;\n\t}\n", difflist->diffheader->freevar.content);
			
			fprintf(file,"\n\tdouble* %s::getSolution(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
			/*** GET SOLUTION ********************/
			
			cur = rewind_list(difvarlist);
			numVar = 0;
			while(cur != NULL)
			{
			    fprintf(file,"\t\tcase %d:\t\treturn %s;    break;\n", numVar, cur->token.content);
			    numVar++;		
			    cur = cur->next;
			}
			
			fprintf(file,"\t\tdefault:\treturn NULL;    break;\n\t\t}\n\t}\n");
			/*** FIM GET SOLUTION ****************/
			
			/*** SOLUCAO EM DISCO ****************/
			fprintf(file,"\tdouble %s::solveToFile(char *filename, char *fileaccess, int firstcall__, int num_iterations__, int num_results__)\n\t{\n\n",classname);
			
			fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
			fprintf(file,"\t\tstatic int num_results_bak = 0;\n");
			fprintf(file,"\t\tstatic int offset_step = 1;\n");
			fprintf(file,"\t\tstatic char *fileaccess_bak = \"\";\n");
			fprintf(file,"\t\tif(firstcall__){\n");
			fprintf(file,"\t\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			
			fprintf(file,"\t\t\tif(num_results__ <= 0)\n\t\t\t\tnum_results__ = 1;\n\t\t\tif(num_iterations__ <= 0)\n\t\t\t\tnum_iterations__ = 1;\n");
			fprintf(file,"\t\t\toffset_step = num_iterations__ / num_results__;\n");
			
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
			    fprintf(file,"\t\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*num_results__);\n", cur->token.content,cur->token.content,cur->token.content);
			    cur = cur->next;
			}
			fprintf(file,"\t\t\tnum_results_bak = num_results__;\n");
			fprintf(file,"\t\t\tnum_iterations_bak = num_iterations__;\n");
			fprintf(file,"\t\t\tfileaccess_bak = fileaccess;\n");
			fprintf(file,"\t\t}\n");
			
			fprintf(file,"\t\tFILE *file = fopen(filename, fileaccess_bak);\n");
			fprintf(file,"\t\tif(!file){\n\t\t\tfprintf(stderr,\"ERROR - solveToFile - Unable to open file %%s\\n\",filename);\n\t\t\texit(1);\n\t\t}\n");
			
			fprintf(file,"\t\tint counter_it__ = 1, it_countx = 1;\n");
			fprintf(file,"\t\tint aux = num_iterations_bak%%num_results_bak;\n");
			fprintf(file,"\t\tfor(it_countx = 1; it_countx<=num_iterations_bak; it_countx++){\n");
			fprintf(file,"\t\t\t%s_new += d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			
			// 08/03/2008//
			print_right_alg(file, classname, alglist, resolved_dep_list);
			///////////////
			print_diff(file, difflist);
			
			
			fprintf(file,"\n\n\t\t\tif(it_countx != aux && (it_countx-aux)%%offset_step == 0){\n");
			fprintf(file,"\t\t\t\tfprintf(file,\"%%.8e");
			int list_size = get_list_size(difvarlist);
			for(int i = 0; i < list_size; i++)
			    fprintf(file," %%.8e");
			fprintf(file,"\\n\",%s_new", difflist->diffheader->freevar.content);
			
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,",%s_new_", cur->token.content);
			    cur = cur->next;
			}
			fprintf(file,");\n");
			fprintf(file,"\t\t\t\tcounter_it__++;\n");
			fprintf(file,"\t\t\t}\n");
			cur = rewind_list(difvarlist);
			while(cur != NULL)
			{
			    fprintf(file,"\t\t\t%s_old_ = %s_new_;\n", cur->token.content,cur->token.content);
			    cur = cur->next;
			}
			fprintf(file,"\t\t}\n");
			fprintf(file,"\t\tfclose(file);\n");
			fprintf(file, "\t\treturn (num_iterations_bak%%offset_step)*d%s;\n", difflist->diffheader->freevar.content);
			fprintf(file,"\t}\n");
			
			/****** FIM SOLUCAO ******************/
			/******* FIM SOLUCAO EM DISCO ****************/
			//print_alg(file, classname, alglist);
			print_if(file, classname, iflist);
			
			fprintf(file,"\n\nfloat __agos_factorial(int f){\n\tif(f>=0 & f<2)\n\t\treturn 1.0;\n\telse if(f < 0)\n\t\treturn 0.0/0.0;\n\tfor(int i=f-1; i>=2; i--)\n\t\tf *= i;\n\treturn (float)f;\n}\n");
			fprintf(file,"double _agos_max(double* vector, int size){\n");
			fprintf(file,"\tdouble max =vector[0];\n");
			fprintf(file,"\tint i;\n");
			fprintf(file,"\tfor(i=1;i<size;i++){\n");
			fprintf(file,"\t\tif(vector[i]>max) max = vector[i];\n");
			fprintf(file,"\t}\n");
			fprintf(file,"\treturn max;\n}\n");
			
			fprintf(file,"double _agos_min(double* vector, int size){\n");
			fprintf(file,"\tdouble min = vector[0];\n");
			fprintf(file,"\tint i;\n");
			
			fprintf(file,"\tfor(i=1;i<size;i++){\n");
			fprintf(file,"\t\tif(vector[i]<min) min = vector[i];\n");
			fprintf(file,"\t}\n");
			fprintf(file,"\treturn min;\n}\n");
			
			fprintf(file,"double _agos_round( double x, int places )\n");
			fprintf(file,"{\n");
			fprintf(file,"\tdouble const shift = powf( 10.0f, places );\n");
			fprintf(file,"\tx *= shift;\n");
			fprintf(file,"\tx = floorf( x + 0.5f );\n");
			fprintf(file,"\tx /= shift;\n");
			fprintf(file,"\treturn x;\n");
			fprintf(file,"}\n");
			
			free(grafoDependencias); //TODO
			fclose(file);
			return 0;
			}
			
			
			//this function prints the cvode functions
			
			
			    
			    int	compile_Sundials_DAPI(char *classname)
			    {
				if(eq_counter <= 0)
				{
				    printf("ERROR - Equation not found\n");
				    exit(1);
				}
				
				// teste da insercao de precedencia 08/03/2008
				preced_alg_list = create_preced_alg_list(rewind_list(alglist));
				AlgList *cural = rewind_list(preced_alg_list);
				TokenNode *resolved_dep_list = NULL;
				
				
				while(preced_alg_list != NULL)
				{
				    cural = rewind_list(preced_alg_list);
				    while(cural != NULL){
					//	printf("%s $$$$$\n", cural->eq->token.content);
					cural = cural->next;
				    }
				    cural = rewind_list(preced_alg_list);
				    
				    //	printf("ASDFASDFASDFAS\n");
				    while(cural != NULL){
					
					TokenNode *cureq = cural->eq->next;
					//	printf("%s \n", cural->eq->token.content);
					//				getchar();
					if(cureq == NULL){
					    resolved_dep_list = add_list(cural->eq->token, resolved_dep_list);
					    preced_alg_list = delete_from_list(preced_alg_list, cural->eq->token);
					    cural = cural->next;
					    //					
					    //					printf("Lista %s\n", resolved_dep_list->token.content);
					    //					getchar();
					}else{
					    
					    while(cureq !=NULL){
						if(list_has_var(cureq->token, resolved_dep_list)){
						    //printf("%s \n", cureq->token.content);
						    //getchar();
						    if(cureq->prev != NULL){
							cureq->prev->next = cureq->next;
						    }
						    if(cureq->next != NULL){
							cureq->next->prev = cureq->prev;
						    }
						}
						cureq = cureq->next;
					    }
					    cural = cural->next;					
					}	
					
				    }
				    
				}
				//		TokenNode *currl = rewind_list(resolved_dep_list);
				//					while(currl != NULL){
			    //						printf("LRD %s\n", currl->token.content);
			    //						currl = currl->next;
			    //		}
			    
			    
			    // fim - teste da insercao de precedencia 08/03/2008
			    
			    
			    FILE *file;
			    char *filename = (char *)calloc(strlen(classname)+5, sizeof(char*));
			    
			    sprintf(filename,"%s.hpp",(const char*)classname);
			    
			    file = fopen(filename, "wb");
			    fprintf(file, "#include \"solver.h\"\n");
			    fprintf(file, "#include <stdio.h>\n#include <math.h>\n\n");
			    
			    fprintf(file, "#include \"sundials_types.h\"\n#include \"cvode.h\"\n#include \"cvode_dense.h\"\n#include \"nvector_serial.h\"\n#include \"sundials_dense.h\"\n");
			    
			    fprintf(file, "#define ABSTOL 1.0E-06\n");
			    fprintf(file, "#define RELTOL 1.0E-04\n");
			    
			    fprintf(file, "static int check_flag(void *flagvalue, char *funcname, int opt);\n");
			    fprintf(file, "static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__);\n");
			    
			    fprintf(file, "class %s: public Solver\n{\n\t//PARAMETERS\npublic:\n", classname);
			    fputs("\tint it_countx;\n",file);
			    
			    
			    /******* DECLARACAO DAS VARIAVEIS *******/
			    // building parvaslist
			    
			    
			    
			    TokenNode *cur = rewind_list(parvarlist);
			    while(cur != NULL)
			    {
				fprintf(file,"\tdouble %s; \t // %s\n",cur->token.content, cur->units);
				cur = cur->next;
			    }
			    
			    //08/03/2008
			    cur = rewind_list(algvarlist);
			    while(cur != NULL)
			    {
				fprintf(file,"\tdouble calc_%s; \t // %s\n",cur->token.content, cur->units);
				//fprintf(file,"\tdouble %s_new;\n",cur->name);
				cur = cur->next;
			    }
			    
			    
			    
			    fprintf(file,"\tdouble d%s, *%s_vec__;\n",difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
			    fprintf(file,"\tdouble %s_new;\n",difflist->diffheader->freevar.content);
			    
			    fprintf(file,"\n\t//DEPENDENT VARIABLES\n");
			    cur = rewind_list(difvarlist);
			    while(cur != NULL)
			    {
				fprintf(file,"\tdouble *%s;\n", cur->token.content);
				fprintf(file,"\tdouble %s_new_, %s_old_, %s_ini_, %s_lado_direito_;\n", cur->token.content,cur->token.content,cur->token.content,cur->token.content);
				cur = cur->next;
			    }
			    // VARIAVEIS PARA O CVODE 
			    fprintf(file, "\t// CVODE VARIABLES\n");
			    fprintf(file,"\trealtype reltol__;\n");
			    fprintf(file,"\tvoid *cvode_mem_cvode__;\n");
			    fprintf(file, "\tN_Vector dependent_variable__, abstol__;\n");
			    fprintf(file, "\tint flag__, flagr__;\n");
			    
			    fputs("\tdouble *depvar__;\n",file);
			    /******* FIM * DECLARACAO DAS VARIAVEIS *******/
			    
			    fprintf(file,"\npublic:\n");
			    fprintf(file, "\t%s(double val_abstol__ = ABSTOL, double val_reltol__= RELTOL, int method__ = CV_BDF);\n", classname);
			    fprintf(file, "\t~%s();\n", classname);
			    
			    fputs("\tvirtual int setVariables(int, double);\n",file);
			    fputs("\tvirtual int setParameters(int, double);\n",file);
			    fputs("\tvirtual int setFreeVariable(double);\n",file);
			    
			    fputs("\tvirtual void setReltol(double value);\n",file);
			    fputs("\tvirtual void setAbstol(int index, double value);\n",file);
			    fputs("\tvirtual void reInitCVODE();\n",file);
			    
			    fputs("\tvirtual void setCVODEMaxStep(double maxstep);\n", file);
			    
			    
			    fputs("\tvirtual double getVariables(int);\n",file);
			    fputs("\tvirtual double getParameters(int);\n",file);
			    fputs("\tvirtual double getFreeVariable();\n",file);
			    
			    fputs("\tvirtual Variables get_Parameters();\n",file);
			    fputs("\tvirtual Variables get_Variables();\n",file);
			    fputs("\tvirtual Variables get_FreeVariable();\n",file);
			    
			    fputs("\tvirtual void setParametersFromFile(char*);\n", file);
			    fputs("\tvirtual void setVariablesFromFile(char*);\n", file);
			    fputs("\tvirtual void setFreeVariableFromFile(char*);\n", file);
			    
			    fputs("\tvirtual double* solveDiff();\n",file);
			    fputs("\tvirtual void solveCVODE(int firstcall__ = 0, int steps__ = 0);\n", file);
			    fputs("\tvirtual double solve(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0,int numThreads=1);\n",file);
			    fputs("\tdouble jacobian(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0, int numThreads=1);\n",file);
			    fputs("\tvirtual double* getSolution(int indVariable);\n",file);
			    fputs("\tvirtual double* getIndependentVar();\n",file);
			    fputs("\tvirtual double solveToFile(char *filename, char *fileaccess = \"\", int firstcall__ = 0, int num_iterations__ = 0, int num_results__ = 0);\n",file);
			    fputs("\tvirtual void solveCVODEToFile(char *filename, char *fileaccess = \"\", int firstcall__ = 0, int steps__ = 0);\n",file);
			    //fputs("\tint solveToFile(char *filename, int numero_iteracoes__, int buffer_size);\n",file);
			    fputs("public:\n",file);
			    
			    
			    // 08/03/2008
			    //		cur = rewind_list(algvarlist);
			    //		while(cur != NULL)
			    //		{
				//			fprintf(file, "\tinline double calc_%s();\n", cur->token.content);
				//			cur = cur->next;
				//		}
				int iif;
				for(iif =0; iif< ifs; iif++)
				{
				    fprintf(file, "\tinline double ifnumber_%d();\n", iif);
				}
				fprintf(file,"};\n\n");
				fprintf(file,"extern \"C\" {\n");
				fprintf(file, "\tSolver *CreateObject(int w);\n");
				fprintf(file, "\tSolver **CreateObject2D(int w, int h);\n");
				fprintf(file, "\tint ode_sizeof();\n");	
				fprintf(file,"}\n");
				
				
				
				// AQUI SEPARA EM 2 ARQUIVOS
				fclose(file);
				
				
				char *cppfilename = (char *)calloc(strlen(classname)+5, sizeof(char*));
				sprintf(cppfilename,"%s.cpp",(const char*)classname);
				
				file = fopen(cppfilename, "wb");
				
				
				fprintf(file, "#include \"%s.hpp\"\n", classname);
				fprintf(file,"#define AGOS_NAN NAN\n#define AGOS_INF \n");
				fprintf(file,"#define __agos_xor(a,b) (!(a && b) && (a || b))\nfloat __agos_factorial(int);\n");
				fprintf(file,"double _agos_max(double*,int);\n");
				fprintf(file,"double _agos_min(double*,int);\n");
				
				//alterado por cadim, 11/09/2009
				//colocando o openmp
				fprintf(file," typedef struct str_funcao{\n\tvoid (*funcao)(Solveode*);\n}typ_funcao;\n\n");
				
				prefixo="classe->";
				//imprime as funções auxiliares
				
				cur = rewind_list(algvarlist);
				int numCalcs =0;
				while(cur != NULL)
				{
				    fprintf(file,"\nvoid _f_calc_%s_(%s *classe){\n",cur->token.content,classname);
				    
				    // 		    file, classname, alglist, resolved_dep_list);
				    // 		   (FILE *file, char *classname, AlgList *list, TokenNode *orderedlist)
				    TokenNode *curlCalc = rewind_list(resolved_dep_list);
				    TokenNode *curCalc = NULL;
				    AlgList *curalgCalc = NULL;
				    int cont=0;
				    while(curlCalc != NULL)
				    {
					if(cont==numCalcs){
					    curalgCalc = rewind_list(alglist);
					    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content))
						curalgCalc = curalgCalc->next;			    
					    fprintf(file,"\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
					    curCalc = curalgCalc->eq;
					    curCalc = curCalc->next->next;
					    print_eq(file, curCalc);
					    fprintf(file,";\t//%d",curalgCalc->number);
					    curlCalc = curlCalc->next;
					    break;
					}else{
					    curCalc = curCalc->next->next;
					    curlCalc = curlCalc->next;
					}
					cont++;
				    }
				    
				    
				    //final da funcao
				    fprintf(file,"\n}\n");
				    cur = cur->next;
				    numCalcs++;
				}
				
				
				
				fprintf(file,"\n\t//DEPENDENT VARIABLES\n");
				cur = rewind_list(difvarlist);
				int numEdos=0;
				while(cur != NULL)
				{
				    fprintf(file,"\nvoid _f_%s_(%s *classe){\n",cur->token.content,classname);
				    DiffList *curlDiff = rewind_list(difflist);
				    TokenNode *curDiff = NULL;
				    int cont=0;
				    while(curlDiff != NULL)
				    {
					if(cont==numEdos){
					    curDiff = curlDiff->diffheader->eq->next;
					    fprintf(file,"\t%s%s_lado_direito_= ",prefixo, curlDiff->diffheader->diffvar.content);
					    print_eq(file, curDiff);
					    fprintf(file,";\n");
					    fprintf(file,"\t%s%s_new_= %sd%s*(%s%s_lado_direito_)+%s%s_old_;\t// %d\n",prefixo, curlDiff->diffheader->diffvar.content, prefixo, curlDiff->diffheader->freevar.content,prefixo, curlDiff->diffheader->diffvar.content,prefixo, curlDiff->diffheader->diffvar.content, curlDiff->number);
					    fprintf(file,"\t//%s%s_old_ = %s%s_new_;\n",prefixo, curlDiff->diffheader->diffvar.content, prefixo, curlDiff->diffheader->diffvar.content);
					    curlDiff = curlDiff->next;
					    break;
					}else{
					    curDiff = curlDiff->diffheader->eq->next;
					    curlDiff = curlDiff->next;
					}
					cont++;
					
				    }
				    fprintf(file,"\n}\n");
				    cur = cur->next;
				    numEdos++;
				}
				prefixo="";
				fprintf(file,"#define numEDO %d\n",numEdos);
				fprintf(file,"#define numAux %d\n",numCalcs);
				
				fprintf(file,"typ_funcao edos[numEDO];\n");
				fprintf(file,"typ_funcao auxiliares[numAux];\n");
				
				//fim -trecho
				
				fprintf(file,"\n\t//CONSTRUCTOR\n");
				fprintf(file, "\t%s::%s(double val_abstol__, double val_reltol__, int method__)\n\t{\n", classname,classname);
				//alterado por cadim
				//11/09/2009
				//openmp
				
				cur = rewind_list(algvarlist);
				numCalcs=0;
				while(cur != NULL){
				    fprintf(file,"\t\tauxiliares[%d].funcao = _f_calc_%s_;\n", numCalcs, cur->token.content);
				    numCalcs++;
				    cur = cur->next;
				}
				
				cur = rewind_list(difvarlist);
				numEdos=0;
				while(cur != NULL){
				    fprintf(file,"\t\tedos[%d].funcao = _f_%s_;\n", numEdos, cur->token.content);
				    // 			    fprintf(file,"auxiliares[%d].funcao = _f_%s_;\n", numCalcs, cur->token.content);
				    numEdos++;
				    // 			    numCalcs++;
				    cur = cur->next;
				}
				
				
				
				
				cur = rewind_list(parvarlist); 		
				
				while(cur != NULL)
				{
				    fprintf(file,"\t\t%s = %.8e;\n", cur->token.content, cur->initialvalue);
				    //fprintf(file,"\t\t%s_new = %e;\n", cur->name, cur->value);
				    cur = cur->next;
				}
				
				fprintf(file,"\t\td%s = 0.0; %s_vec__ = NULL;\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
				
				// INICIALIZACAO CVODE
				fprintf(file, "\n");
				fprintf(file, "\t\tdependent_variable__ = N_VNew_Serial(%d);\n", get_list_size(difvarlist));
				fprintf(file, "\t\tif(check_flag((void *)dependent_variable__, \"N_VNew_Serial\", 0))\n");
				fprintf(file, "\t\t\texit(1);\n");
				fprintf(file, "\t\tabstol__ = N_VNew_Serial(%d);\n", get_list_size(difvarlist));
				fprintf(file, "\t\tif (check_flag((void *)abstol__, \"N_VNew_Serial\", 0))\n");
				fprintf(file, "\t\t\texit(1);\n\n");
				
				
				fprintf(file,"\t\tdepvar__ = (double*)malloc(sizeof(double)*%d);\n", get_list_size(difvarlist));
				fprintf(file,"\t\tif(depvar__ == NULL){\n\t\t\tfprintf(stderr, \"ERROR Cannot allocate memory for depvar__\\n\");\n\t\t\texit(0);\n\t\t}\n");
				
				
				cur = rewind_list(difvarlist);
				int i =0;
				while(cur != NULL)
				{
				    fprintf(file,"\t\t%s = NULL;\n",cur->token.content);
				    fprintf(file,"\t\tNV_Ith_S(dependent_variable__, %d) = %s_ini_ = %.8e;", i, cur->token.content, cur->initialvalue);
				    fprintf(file,"\tNV_Ith_S(abstol__, %d) = val_abstol__;\n", i);
				    cur = cur->next;
				    i++;
				}
				fprintf(file, "\t\treltol__ = val_reltol__;\n");
				fprintf(file,"\t\tit_countx = 0;\n");
				fprintf(file,"\t\tint nonlineariteration__;\n");
				fprintf(file, "\t\tif(method__ == CV_BDF) nonlineariteration__ = CV_NEWTON;\n");
				fprintf(file, "\t\telse nonlineariteration__ = CV_FUNCTIONAL;\n");
				fprintf(file, "\t\tcvode_mem_cvode__ = CVodeCreate(method__, nonlineariteration__);\n");
				fprintf(file, "\t\tif (check_flag((void *)cvode_mem_cvode__, \"CVodeCreate\", 0))\n");
				fprintf(file, "\t\t\texit(1);\n\n");
				fprintf(file, "\t\tflag__ = CVodeMalloc(cvode_mem_cvode__, f__, %s, dependent_variable__, CV_SV, reltol__, abstol__);\n", rewind_list(parvarlist)->token.content);
				fprintf(file, "\t\tif (check_flag(&flag__, \"CVodeMalloc\", 1))\n");
				fprintf(file, "\t\t\texit(1);\n\n");
				fprintf(file, "\t\tflag__ = CVDense(cvode_mem_cvode__, %d);\n", get_list_size(difvarlist));
				fprintf(file, "\t\tif (check_flag(&flag__, \"CVDense\", 1))\n");
				fprintf(file, "\t\t\texit(1);\n");
				
				
				//Alterado por Ricardo 
				//fprintf(file, "\t\tflag__ = CVodeSetMaxStep(cvode_mem_cvode__, %s);\n", rewind_list(parvarlist)->token.content);
				//fprintf(file, "\t\tfprintf(stderr, \"%%p\t\", maxStep);\n");
				// 		fprintf(file, "\t\tflag__ = CVodeSetMaxStep(cvode_mem_cvode__, maxStep);\n" );
				// 		fprintf(file, "\t\tif (check_flag(&flag__, \"CVodeSetMaxStep\", 1))\n");
				// 		fprintf(file, "\t\t\texit(1);\n\n");
				
				
				//fim-trecho	
				
				fprintf(file, "\t\tCVodeSetFdata(cvode_mem_cvode__, (void*)this);\n");
				
				fprintf(file,"\t}\n");
				
				/****** FIM CONSTRUTOR *************************/
				
				/****** DESTRUTOR ******************************/
				fprintf(file, "\t%s::~%s()\n\t{\n", classname,classname);
				fprintf(file, "\t\tif(depvar__ != NULL) free(depvar__);\n");
				
				cur = rewind_list(difvarlist);
				while(cur != NULL){
				    fprintf(file, "\t\tif(%s != NULL) free(%s);\n", cur->token.content, cur->token.content);
				    cur = cur->next;
				    
				}
				
				
				fprintf(file, "\t\tN_VDestroy_Serial(dependent_variable__);\n");
				fprintf(file, "\t\tN_VDestroy_Serial(abstol__);\n");
				fprintf(file, "\t\tCVodeFree(&cvode_mem_cvode__);\n");
				fprintf(file, "\t}\n");
				
				
				
				/****** FIM DESTRUTOR **************************/
				
				
				/****** METODOS SET ****************************/
				fprintf(file,"\n\tint %s::setVariables(int indVariable, double value_new)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
				
				cur = rewind_list(difvarlist);
				int numVar = 0;
				while(cur != NULL)
				{
				    fprintf(file,"\t\tcase %d:\t\tNV_Ith_S(dependent_variable__, %d) = %s_old_ = %s_ini_ = value_new;    break;\n", numVar, numVar, cur->token.content, cur->token.content);
				    numVar++;
				    cur = cur->next;
				}
				fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t\treturn 0;\n\t}\n");
				
				fprintf(file,"\n\tint %s::setParameters(int indVariable, double value_new)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
				
				cur = rewind_list(parvarlist);
				numVar = 0;
				while(cur != NULL)
				{
				    fprintf(file,"\t\tcase %d:\t\t%s = value_new;   break;\n", numVar, cur->token.content);
				    numVar++;
				    cur = cur->next;
				}
				fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t\treturn 0;\n\t}\n");
				
				
				fprintf(file,"\n\tint %s::setFreeVariable(double value_new)\n\t{\n\t\t",classname);
				fprintf(file,"d%s = value_new;\n\t}\n", difflist->diffheader->freevar.content);	
				
				/******* FIM SET ***********************/
				
				/******* METODOS GET ********************/
				fprintf(file,"\n\t//Get Methods\n");
				
				fprintf(file,"\n\tdouble %s::getVariables(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
				
				cur = rewind_list(difvarlist);
				numVar = 0;
				while(cur != NULL)
				{
				    fprintf(file,"\t\tcase %d:\t\treturn %s_old_;    break;\n", numVar, cur->token.content);
				    numVar++;
				    cur = cur->next;
				}
				fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t}\n");
				
				/////////////
				fprintf(file,"\n\tdouble %s::getLadoDireito(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
				
				cur = rewind_list(difvarlist);
				numVar = 0;
				while(cur != NULL)
				{
				    fprintf(file,"\t\tcase %d:\t\treturn %s_lado_direito_;    break;\n", numVar, cur->token.content);
				    numVar++;
				    cur = cur->next;
				}
				fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t}\n");
				
				
				
				fprintf(file,"\n\tdouble %s::getParameters(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
				cur = rewind_list(parvarlist);
				numVar = 0;
				while(cur != NULL)
				{
				    fprintf(file,"\t\tcase %d:\t\treturn %s;    break;\n", numVar, cur->token.content);
				    numVar++;
				    cur = cur->next;
				}
				fprintf(file,"\t\tdefault:\tbreak;\n\t\t}\n\t}\n");
				
				fprintf(file,"\n\tdouble %s::getFreeVariable()\n\t{\n",classname);
				
				fprintf(file,"\t\treturn d%s;\n\t}\n", difflist->diffheader->freevar.content);
				
				
				fprintf(file,"\n\t//Get Methods - Variables\n\n");
				
				fprintf(file, "\tVariables %s::get_Variables()\n\t{\n\t\tVariables v(\"",classname);
				
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,"|%s#",cur->token.content);
				    cur = cur->next;
				}
				fprintf(file, "\");\n\t\treturn v;\n\t}\n");
				
				fprintf(file, "\tVariables %s::get_Parameters()\n\t{\n\t\tVariables v(\"",classname);
				
				cur = rewind_list(parvarlist);
				while(cur != NULL)
				{
				    fprintf(file,"|%s#",cur->token.content);
				    cur = cur->next;
				}
				fprintf(file, "\");\n\t\treturn v;\n\t}\n");	
				
				fprintf(file, "\tVariables %s::get_FreeVariable()\n\t{\n\t\tVariables v(\"",classname);
				
				fprintf(file,"|%s#",difflist->diffheader->freevar.content); 		
				fprintf(file, "\");\n\t\treturn v;\n\t}\n");	
				
				
				/********SET DO ARQUIVO **********************/
				fprintf(file,"\n\tvoid %s::setParametersFromFile(char *filename)\n\t{\n",classname);
				fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
				fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setParametersFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
				fprintf(file,"\tdouble value;\n\t\tint k = 0;\n\t\tVariables v = get_Parameters();\n\t\tint s = v.getQuantity();\n\t");
				fprintf(file,"\tfor(;k<s;k++)\n\t\t{\n\t\t");
				fprintf(file,"\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetParameters(k, value);\n\t\t}\n\t\tfclose(file);\n\t}\n");
				
				fprintf(file,"\n\tvoid %s::setVariablesFromFile(char *filename)\n\t{\n", classname);
				fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
				fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setVariablesFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
				fprintf(file,"\tdouble value;\n\t\tint k = 0;\n\t\tVariables v = get_Variables();\n\t\tint s = v.getQuantity();\n\t");
				fprintf(file,"\tfor(;k<s;k++)\n\t\t{\n\t\t");
				fprintf(file,"\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetVariables(k, value);\n\t\t}\n\t\tfclose(file);\n\t}\n");
				
				fprintf(file,"\n\tvoid %s::setFreeVariableFromFile(char *filename)\n\t{\n", classname);
				fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
				fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setFreeVariableFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
				fprintf(file,"\tdouble value;\n\t\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetFreeVariable(value);\n\t\tfclose(file);\n\t}\n");
				/********FIM SET DO ARQUIVO ******************/
				
				/** METODO SOLVEDIFF *************************************************/
				fprintf(file,"\n\tdouble* %s::solveDiff()\n\t{\n",classname);
				// 08/03/2008//
				print_right_alg(file, classname, alglist, resolved_dep_list);
				///////////////
				print_only_diff(file, difflist);
				fprintf(file, "\n\t\treturn depvar__;\n\t}\n\n");
				/** FIM SOLVEDIFF ****************************************************/
				
				
				/** SOLVECVODE *******************************************************/
				fprintf(file,"\n\tvoid %s::solveCVODE(int firstcall__, int steps__)\n\t{\n",classname);
				fprintf(file,"\t\tint iout = 0;\n");
				fprintf(file,"\t\trealtype tout = %s_new+d%s;\n\n",rewind_list(parvarlist)->token.content,rewind_list(parvarlist)->token.content );
				
				fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
				fprintf(file,"\t\tif(firstcall__){\n");
				fprintf(file,"\t\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
				fprintf(file,"\t\t\ttout = %s + d%s;\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
				
				fprintf(file,"\t\t\tif(steps__ <= 0)\n\t\t\t\tsteps__ = 1;\n");
				
				fprintf(file,"\t\t\tif(%s_vec__ != NULL)free( %s_vec__);\n\t\t\t%s_vec__ = (double *)malloc(sizeof(double)*steps__);\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
				    fprintf(file,"\t\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*steps__);\n", cur->token.content,cur->token.content,cur->token.content);
				    cur = cur->next;
				}
				fprintf(file,"\t\t\tnum_iterations_bak = steps__;\n");
				fprintf(file,"\t\t}\n");
				
				fprintf(file,"\t\twhile(1){\n");
				fprintf(file,"\t\t\tflag__ = CVodeSetStopTime(cvode_mem_cvode__, tout);\n");
				fprintf(file,"\t\t\tflag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &%s_new, CV_NORMAL_TSTOP);\n",rewind_list(parvarlist)->token.content);
				//Trecho alterado por Ricardo Silva Campos. Eu comentei a linha acima e troquei pelas abaixo
				//fprintf(file, "\t\tfprintf(stderr, \"%%p\t\", maxStep);\n");
				//fprintf(file, "\t\tprintf(\"Maxstep: %s\n\", maxStep);\n");
				//fprintf(file,"\t\t\t flag__ = CVodeSetMaxStep(cvode_mem_cvode__, maxStep);\n" );
				//fprintf(file,"\t\t\tflag__ = CVodeSetMaxStep(cvode_mem_cvode__, dtime);\n");
				//fprintf(file,"\t\t\tflag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &time_new, CV_NORMAL);\n");
				//fim-trecho
				
				fprintf(file,"\t\t\t%s_vec__[iout] = %s_new;\n",rewind_list(parvarlist)->token.content, rewind_list(parvarlist)->token.content); 
				cur = rewind_list(difvarlist);
				i =0;
				while(cur != NULL){
				    fprintf(file, "\t\t\t%s_old_ = %s[iout] = NV_Ith_S(dependent_variable__, %d);\n", cur->token.content, cur->token.content, i);
				    cur = cur->next;
				    i++;
				}
				
				fprintf(file,"\t\t\tif (check_flag(&flag__, \"CVode\", 1)) break;\n");
				fprintf(file,"\t\t\tif (flag__ == CV_SUCCESS){\n");
				fprintf(file,"\t\t\t\tiout++;\n\t\t\t\ttout += d%s; // timestep\n", rewind_list(parvarlist)->token.content);
				fprintf(file,"\t\t\t}\n\t\t\tif (iout == num_iterations_bak ) break;\n\t\t}\n"); 	
				
				
				fprintf(file, "\t}\n");
				
				/** FIM SOLVECVODE ***************************************************/ 				
				
				/******* METODO SOLVE ******************/
				fprintf(file,"\tdouble %s::solve(int firstcall__, int num_iterations__, int num_results__, int numThreads)\n\t{\n\n",classname);
				
				//fprintf(file,"// 3 testeeeeeeeeeeeeee");
				
				fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
				fprintf(file,"\t\tstatic int num_results_bak = 0;\n");
				fprintf(file,"\t\tstatic int offset_step = 1;\n");
				fprintf(file,"\t\tif(firstcall__){\n");
				fprintf(file,"\t\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
				
				fprintf(file,"\t\t\tif(num_results__ <= 0)\n\t\t\t\tnum_results__ = 1;\n\t\t\tif(num_iterations__ <= 0)\n\t\t\t\tnum_iterations__ = 1;\n");
				fprintf(file,"\t\t\toffset_step = num_iterations__ / num_results__;\n");
				
				fprintf(file,"\t\t\tif(%s_vec__ != NULL)free( %s_vec__);\n\t\t\t%s_vec__ = (double *)malloc(sizeof(double)*num_results__);\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
				    fprintf(file,"\t\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*num_results__);\n", cur->token.content,cur->token.content,cur->token.content);
				    cur = cur->next;
				}
				fprintf(file,"\t\t\tnum_results_bak = num_results__;\n");
				fprintf(file,"\t\t\tnum_iterations_bak = num_iterations__;\n");
				fprintf(file,"\t\t}\n");
				
				fprintf(file,"\t\tint counter_it__ = 1, it_countx = 1;\n");
				fprintf(file,"\t\tint aux = num_iterations_bak%%num_results_bak;\n");
				
				//fprintf(file,"\t\tomp_set_num_threads(numThreads);\n");
				//fprintf(file,"\t\t#pragma omp parallel for ordered schedule(static)\n");
				
				fprintf(file,"\t\tfor(it_countx = 1; it_countx<=num_iterations_bak; it_countx++){\n");
				fprintf(file,"\t\t\t%s_new += d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
				
				// 08/03/2008//
				//        print_right_alg(file, classname, alglist, resolved_dep_list);
				///////////////
				//        print_diff(file, difflist);
				
				
				fprintf(file,"\t\t\tfor(int k=0;k<numAux;k++)\t");
				fprintf(file,"auxiliares[k].funcao(this);\n");
				
				fprintf(file,"\t\t\tfor(int k=0;k<numEDO;k++)\t");
				fprintf(file,"edos[k].funcao(this);\n");
				
				//fprintf(file,"\t\t\t#pragma omp critical\n");
				fprintf(file,"\t\t\tif(it_countx != aux && (it_countx-aux)%%offset_step == 0){\n");
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,"\t\t\t\t%s[counter_it__] = %s_new_;\n", cur->token.content,cur->token.content);
				    cur = cur->next;
				}
				fprintf(file, "\t\t\t\t%s_vec__[counter_it__] = %s_new;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
				fprintf(file,"\t\t\t\tcounter_it__++;\n");
				fprintf(file,"\t\t\t}\n");
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,"\t\t\t%s_old_ = %s_new_;\n", cur->token.content,cur->token.content);
				    cur = cur->next;
				}
				fprintf(file,"\t\t}\n");
				
				
				fprintf(file,"\t\t\t//FILE *fileptr;\n");
				fprintf(file,"\t\t\t//char filename[12];\n");
				fprintf(file,"\t\t\t//sprintf(filename,\"%%dthread.dat\",numThreads);\n");
				fprintf(file,"\t\t\t//fileptr = fopen(filename, \"wb\");\n");
				fprintf(file,"\t\t\t//for(int i =0;i<num_results__;i++){\n");
				fprintf(file,"\t\t\t//    fprintf(fileptr,\"%%f %%f\\n\", time_vec__[i], V[i]);\n");
				fprintf(file,"\t\t\t//}\n");
				fprintf(file,"\t\t\t//fclose(fileptr);\n");
				
				
				
				
				// 		fprintf(file,"\t\t\t \n",);
				fprintf(file, "\t\treturn (num_iterations_bak%%offset_step)*d%s;\n", difflist->diffheader->freevar.content);
				fprintf(file,"\t}\n");
				/****** FIM SOLUCAO ******************/
				
				fprintf(file, "\n\tdouble* %s::getIndependentVar()\n\t{\n", classname);
				fprintf(file, "\t\treturn %s_vec__;\n\t}\n", difflist->diffheader->freevar.content);
				
				
				fprintf(file,"\n\tdouble* %s::getSolution(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
				/*** GET SOLUTION ********************/
				
				cur = rewind_list(difvarlist);
				numVar = 0;
				while(cur != NULL)
				{
				    fprintf(file,"\t\tcase %d:\t\treturn %s;    break;\n", numVar, cur->token.content);
				    numVar++;		
				    cur = cur->next;
				}
				
				fprintf(file,"\t\tdefault:\treturn NULL;    break;\n\t\t}\n\t}\n");
				/*** FIM GET SOLUTION ****************/
				
				/*** SOLUCAO EM DISCO ****************/
				fprintf(file,"\tdouble %s::solveToFile(char *filename, char *fileaccess, int firstcall__, int num_iterations__, int num_results__)\n\t{\n\n",classname);
				
				fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
				fprintf(file,"\t\tstatic int num_results_bak = 0;\n");
				fprintf(file,"\t\tstatic int offset_step = 1;\n");
				fprintf(file,"\t\tstatic char *fileaccess_bak = \"\";\n");
				fprintf(file,"\t\tif(firstcall__){\n");
				fprintf(file,"\t\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
				
				fprintf(file,"\t\t\tif(num_results__ <= 0)\n\t\t\t\tnum_results__ = 1;\n\t\t\tif(num_iterations__ <= 0)\n\t\t\t\tnum_iterations__ = 1;\n");
				fprintf(file,"\t\t\toffset_step = num_iterations__ / num_results__;\n");
				
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
				    fprintf(file,"\t\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*num_results__);\n", cur->token.content,cur->token.content,cur->token.content);
				    cur = cur->next;
				}
				fprintf(file,"\t\t\tnum_results_bak = num_results__;\n");
				fprintf(file,"\t\t\tnum_iterations_bak = num_iterations__;\n");
				fprintf(file,"\t\t\tfileaccess_bak = fileaccess;\n");
				fprintf(file,"\t\t}\n");
				
				fprintf(file,"\t\tFILE *file = fopen(filename, fileaccess_bak);\n");
				fprintf(file,"\t\tif(!file){\n\t\t\tfprintf(stderr,\"ERROR - solveToFile - Unable to open file %%s\\n\",filename);\n\t\t\texit(1);\n\t\t}\n");
				
				fprintf(file,"\t\tint counter_it__ = 1, it_countx = 1;\n");
				fprintf(file,"\t\tint aux = num_iterations_bak%%num_results_bak;\n");
				fprintf(file,"\t\tfor(it_countx = 1; it_countx<=num_iterations_bak; it_countx++){\n");
				fprintf(file,"\t\t\t%s_new += d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
				
				// 08/03/2008//
				print_right_alg(file, classname, alglist, resolved_dep_list);
				///////////////
				print_diff(file, difflist);
				
				fprintf(file,"\n\n\t\t\tif(it_countx != aux && (it_countx-aux)%%offset_step == 0){\n");
				fprintf(file,"\t\t\t\tfprintf(file,\"%%.8e");
				int list_size = get_list_size(difvarlist);
				for(int i = 0; i < list_size; i++)
				    fprintf(file," %%.8e");
				fprintf(file,"\\n\",%s_new", difflist->diffheader->freevar.content);
				
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,",%s_new_", cur->token.content);
				    cur = cur->next;
				}
				fprintf(file,");\n");
				fprintf(file,"\t\t\t\tcounter_it__++;\n");
				fprintf(file,"\t\t\t}\n");
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,"\t\t\t%s_old_ = %s_new_;\n", cur->token.content,cur->token.content);
				    cur = cur->next;
				}
				fprintf(file,"\t\t}\n");
				fprintf(file,"\t\tfclose(file);\n");
				fprintf(file, "\t\treturn (num_iterations_bak%%offset_step)*d%s;\n", difflist->diffheader->freevar.content);
				fprintf(file,"\t}\n");
				
				/****** FIM SOLUCAO ******************/
				
				fprintf(file,"\tvoid %s::solveCVODEToFile(char *filename, char *fileaccess, int firstcall__, int steps__)\n\t{\n\n",classname);
				
				fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
				fprintf(file,"\t\tstatic char *fileaccess_bak = \"\";\n");
				fprintf(file,"\t\tif(firstcall__){\n");
				
				
				fprintf(file,"\t\t\tif(steps__ <= 0)\n\t\t\t\tsteps__ = 1;\n");
				fprintf(file,"\t\t\tnum_iterations_bak = steps__;\n");
				fprintf(file,"\t\t\tfileaccess_bak = fileaccess;\n");
				fprintf(file,"\t\t}\n");
				
				fprintf(file,"\t\tFILE *file = fopen(filename, fileaccess_bak);\n");
				fprintf(file,"\t\tif(!file){\n\t\t\tfprintf(stderr,\"ERROR - solveCVODEToFile - Unable to open file %%s\\n\",filename);\n\t\t\texit(1);\n\t\t}\n");
				fprintf(file,"\t\tsolveCVODE(firstcall__, num_iterations_bak);\n");
				fprintf(file,"\t\tfor(int i=0;i<num_iterations_bak; i++){\n");
				fprintf(file,"\t\t\tfprintf(file,\"%%.8e");
				list_size = get_list_size(difvarlist);
				for(int i = 0; i < list_size; i++)
				    fprintf(file," %%.8e");
				fprintf(file,"\\n\",%s_vec__[i]", difflist->diffheader->freevar.content);
				cur = rewind_list(difvarlist);
				while(cur != NULL)
				{
				    fprintf(file,",%s[i]", cur->token.content);
				    cur = cur->next;
				}
				fprintf(file,");\n\t\t}\n");
				fprintf(file,"\t\tfclose(file);\n");
				fprintf(file,"\t}\n\n");
				
				
				fprintf(file,"\tvoid %s::reInitCVODE()\n\t{\n\n",classname);
				fprintf(file,"\t\tflag__ = CVodeReInit(cvode_mem_cvode__, f__, %s, dependent_variable__, CV_SV, reltol__, abstol__);\n", difflist->diffheader->freevar.content);
				fprintf(file,"\t\tif (check_flag(&flag__, \"CVodeReInit\", 1))\n\t\t\texit(1);\n");
				fprintf(file, "\t}\n\n");
				
				fprintf(file, "\tvoid %s::setCVODEMaxStep(double maxstep)\n\t{\n\n",classname, file);
				//Alterado por Ricardo 
				fprintf(file, "\t\tflag__ = CVodeSetMaxStep(cvode_mem_cvode__, maxstep);\n", rewind_list(parvarlist)->token.content);
				fprintf(file, "\t\tif (check_flag(&flag__, \"CVodeSetMaxStep\", 1))\n");
				fprintf(file, "\t\t\texit(1);\n\n");
				fprintf(file, "\t}\n\n");
				
				
				fprintf(file, "\tvoid %s::setReltol(double value_new)\n\t{\n",classname);
				fprintf(file, "\t\treltol__ = value_new;\n\t}\n\n");
				
				fprintf(file, "\tvoid %s::setAbstol(int index, double value_new)\n\t{\n",classname);
				fprintf(file, "\t\tswitch(index){\n");
				cur = rewind_list(difvarlist);
				numVar = 0;
				while(cur != NULL)
				{
				    fprintf(file,"\t\tcase %d:\t\tNV_Ith_S(abstol__, %d) = value_new;    break;\n", numVar, numVar);
				    numVar++;
				    cur = cur->next;
				}
				fprintf(file,"\t\tdefault: fprintf(stderr,\"ERROR - setAbstol - index = %%d out of bounds\\n\",index);    break;\n\t\t}\n\t}\n");
				
				
				
				/******* FIM SOLUCAO EM DISCO ****************/
				//print_alg(file, classname, alglist);
				print_if(file, classname, iflist);
				
				
				fprintf(file, "\n\n");
				fprintf(file, "static int check_flag(void *flagvalue, char *funcname, int opt){\n");
				fprintf(file, "\tint *errflag;\n");
				fprintf(file, "\tif (opt == 0 && flagvalue == NULL) {\n");
				fprintf(file, "\t\tfprintf(stderr, \"\\nSUNDIALS_ERROR: %%s() failed - returned NULL pointer\\n\\n\",funcname);\n");
				fprintf(file, "\t\treturn(1);}\n");
				fprintf(file, "\telse if (opt == 1) {\n");
				fprintf(file, "\t\terrflag = (int *) flagvalue;\n");
				fprintf(file, "\t\tif (*errflag < 0) {\n");
				fprintf(file, "\t\t\tfprintf(stderr, \"\\nSUNDIALS_ERROR: %%s() failed with flag = %%d\\n\\n\",funcname, *errflag);\n");
				fprintf(file, "\t\t\treturn(1); }}\n");
				fprintf(file, "\telse if (opt == 2 && flagvalue == NULL) {\n");
				fprintf(file, "\t\tfprintf(stderr, \"\\nMEMORY_ERROR: %%s() failed - returned NULL pointer\\n\\n\",funcname);\n");
				fprintf(file, "\t\treturn(1); }\n");
				fprintf(file, "\treturn 0;\n}\n\n");
				
				fprintf(file, "static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__){\n");
				fprintf(file, "\t%s *ode = (%s *) f_data__;\n", classname, classname);
				
				
				int listsize = get_list_size(difvarlist);
				for(int i=0; i<listsize; i++)
				    fprintf(file, "\tode->setVariables( %d ,NV_Ith_S(dependent_variable__, %d));\n", i, i);
				
				fprintf(file, "\tode->setParameters(0,time);\n");
				fprintf(file, "\tdouble *t = ode->solveDiff();\n");
				
				for(int i=0; i<listsize; i++)
				    fprintf(file, "\tNV_Ith_S(dep_var_dot__, %d) = t[%d];\n", i, i);
				fprintf(file, "\treturn 0;\n}\n\n");
				
				fprintf(file,"Solver *CreateObject(int w){\n");
				fprintf(file,"\t%s *solve_ode = new %s[w];\n", classname, classname);
				fprintf(file,"\treturn ((Solver*)solve_ode);\n");
				fprintf(file,"}\n");
				
				fprintf(file,"\n\nSolver **CreateObject2D(int w, int h){\n");
				fprintf(file,"\t%s **solve_ode = new %s*[h];\n", classname, classname);
				fprintf(file,"\tfor(int i=0;i<h; i++)\n");
				fprintf(file,"\t\tsolve_ode[i] = new %s[w];\n",classname);
				fprintf(file,"\treturn ((Solver**)solve_ode);\n");
				fprintf(file,"}\n");
				
				fprintf(file,"int ode_sizeof(){\n");
				fprintf(file,"\treturn sizeof(%s);\n", classname);
				fprintf(file,"}\n");
				fprintf(file,"\n\nfloat __agos_factorial(int f){\n\tif(f>=0 & f<2)\n\t\treturn 1.0;\n\telse if(f < 0)\n\t\treturn 0.0/0.0;\n\tfor(int i=f-1; i>=2; i--)\n\t\tf *= i;\n\treturn (float)f;\n}\n");
				fprintf(file,"double _agos_max(double* vector, int size){\n");
				fprintf(file,"\tdouble max =vector[0];\n");
				fprintf(file,"\tint i;\n");
				fprintf(file,"\tfor(i=1;i<size;i++){\n");
				fprintf(file,"\t\tif(vector[i]>max) max = vector[i];\n");
				fprintf(file,"\t}\n");
				fprintf(file,"\treturn max;\n}\n");
				
				fprintf(file,"double _agos_min(double* vector, int size){\n");
				fprintf(file,"\tdouble min = vector[0];\n");
				fprintf(file,"\tint i;\n");
				
				fprintf(file,"\tfor(i=1;i<size;i++){\n");
				fprintf(file,"\t\tif(vector[i]<min) min = vector[i];\n");
				fprintf(file,"\t}\n");
				fprintf(file,"\treturn min;\n}\n");
				fclose(file);
				return 0;	
				}
				
				int	compile_DAPI(char *classname)
				{
				    if(eq_counter <= 0)
				    {
					printf("ERROR - Equation not found\n");
					exit(1);
				    }
				    
				    
				    // teste da insercao de precedencia 08/03/2008
				    preced_alg_list = create_preced_alg_list(rewind_list(alglist));
				    AlgList *cural = rewind_list(preced_alg_list);
				    TokenNode *resolved_dep_list = NULL;
				    
				    
				    while(preced_alg_list != NULL)
				    {
					cural = rewind_list(preced_alg_list);
					while(cural != NULL){
					    //	printf("%s $$$$$\n", cural->eq->token.content);
					    cural = cural->next;
					}
					cural = rewind_list(preced_alg_list);
					
					//	printf("ASDFASDFASDFAS\n");
					while(cural != NULL){
					    
					    TokenNode *cureq = cural->eq->next;
					    //	printf("%s \n", cural->eq->token.content);
					    //				getchar();
					    if(cureq == NULL){
						resolved_dep_list = add_list(cural->eq->token, resolved_dep_list);
						preced_alg_list = delete_from_list(preced_alg_list, cural->eq->token);
						cural = cural->next;
						//					
						//					printf("Lista %s\n", resolved_dep_list->token.content);
						//					getchar();
					    }else{
						
						while(cureq !=NULL){
						    if(list_has_var(cureq->token, resolved_dep_list)){
							//printf("%s \n", cureq->token.content);
							//getchar();
							if(cureq->prev != NULL){
							    cureq->prev->next = cureq->next;
							}
							if(cureq->next != NULL){
							    cureq->next->prev = cureq->prev;
							}
						    }
						    cureq = cureq->next;
						}
						cural = cural->next;					
					    }	
					    
					}
					
				    }
				    //		TokenNode *currl = rewind_list(resolved_dep_list);
				    //					while(currl != NULL){
				    //						printf("LRD %s\n", currl->token.content);
				    //						currl = currl->next;
				    //		}
				    
				    
				    // fim - teste da insercao de precedencia 08/03/2008
				    
				    
				    FILE *file;
				    char *filename = (char *)calloc(strlen(classname)+5, sizeof(char*));
				    
				    sprintf(filename,"%s.hpp",(const char*)classname);
				    
				    file = fopen(filename, "wb");
				    fprintf(file, "#include \"solver.h\"\n");
				    fprintf(file, "#include <stdio.h>\n#include <math.h>\n#include \"MCutil.hpp\"\n");
				    fprintf(file, "#include <omp.h>\n#include \"Stopwatch.h\"\n");
				    fprintf(file, "class %s: public Solver\n{\n\t//PARAMETERS\npublic:\n", classname);
				    fputs("\tint it_countx;\n",file);
				    
				    
				    /******* DECLARACAO DAS VARIAVEIS *******/
				    // building parvaslist
				    
				    TokenNode *cur = rewind_list(parvarlist);
				    while(cur != NULL)
				    {
					fprintf(file,"\tdouble %s; \t // %s\n",cur->token.content, cur->units);
					//fprintf(file,"\tdouble %s_new;\n",cur->name);
					cur = cur->next;
				    }
				    
				    //08/03/2008
				    cur = rewind_list(algvarlist);
				    while(cur != NULL)
				    {
					fprintf(file,"\tdouble calc_%s; \t // %s\n",cur->token.content, cur->units);
					//fprintf(file,"\tdouble %s_new;\n",cur->name);
					cur = cur->next;
				    }
				    
				    
				    fprintf(file,"\tdouble d%s, *%s_vec__;\n",difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
				    fprintf(file,"\tdouble %s_new;\n",difflist->diffheader->freevar.content);
				    
				    fprintf(file,"\n\t//function variables\n");
				    cur = rewind_list(difvarlist);
				    while(cur != NULL)
				    {
					fprintf(file,"\tdouble *%s;\n", cur->token.content);
					fprintf(file,"\tdouble %s_new_, %s_old_, %s_ini_,%s_lado_direito_;\n", cur->token.content,cur->token.content,cur->token.content, cur->token.content);
					cur = cur->next;
				    }
				    /******* FIM * DECLARACAO DAS VARIAVEIS *******/
				    
				    fprintf(file,"\npublic:\n");
				    fprintf(file, "\t%s();\n", classname);
				    fprintf(file, "\t~%s();\n", classname);
				    
				    fputs("\tvirtual int setVariables(int, double);\n",file);
				    fputs("\tvirtual int setParameters(int, double);\n",file);
				    fputs("\tvirtual int setFreeVariable(double);\n",file);
				    
				    fputs("\tvirtual double getVariables(int);\n",file);
				    fputs("\tvirtual double getParameters(int);\n",file);
				    fputs("\tvirtual double getFreeVariable();\n",file);
				    
				    fputs("\tvirtual Variables get_Parameters();\n",file);
				    fputs("\tvirtual Variables get_Variables();\n",file);
				    fputs("\tvirtual Variables get_FreeVariable();\n",file);
				    
				    fputs("\tvirtual void setParametersFromFile(char*);\n", file);
				    fputs("\tvirtual void setVariablesFromFile(char*);\n", file);
				    fputs("\tvirtual void setFreeVariableFromFile(char*);\n", file);
				    
				    fputs("\tvirtual double solve(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0, int numThreads=1);\n",file);
				    fputs("\tvirtual double* getSolution(int indVariable);\n",file);
				    fputs("\tvirtual double* getIndependentVar();\n",file);
				    fputs("\tvirtual double solveToFile(char *filename, char *fileaccess = \"\", int firstcall__ = 0, int num_iterations__ = 0, int num_results__ = 0);\n",file);
				    //fputs("\tint solveToFile(char *filename, int numero_iteracoes__, int buffer_size);\n",file);
				    fputs("public:\n",file);
				    
				    // 08/03/2008
				    //		cur = rewind_list(algvarlist);
				    //		while(cur != NULL)
				    //		{
					//			fprintf(file, "\tinline double calc_%s();\n", cur->token.content);
					//			cur = cur->next;
					//		}
					int iif;
					for(iif =0; iif< ifs; iif++)
					{
					    fprintf(file, "\tinline double ifnumber_%d();\n", iif);
					}
					fprintf(file,"};\n\n");
					
					fprintf(file,"extern \"C\" {\n");
					fprintf(file, "\tSolver *CreateObject(int w);\n");
					fprintf(file, "\tSolver **CreateObject2D(int w, int h);\n");
					fprintf(file, "\tint ode_sizeof();\n");				
					fprintf(file,"}\n");
					
					
					
					// AQUI SEPARA EM 2 ARQUIVOS
					fclose(file);
					
					
					char *cppfilename = (char *)calloc(strlen(classname)+5, sizeof(char*));
					sprintf(cppfilename,"%s.cpp",(const char*)classname);
					
					file = fopen(cppfilename, "wb");
					
					fprintf(file, "#include \"%s.hpp\"\n", classname);
					
					fprintf(file,"#define AGOS_NAN NAN\n#define AGOS_INF INFINITY\n");
					fprintf(file,"#define __agos_xor(a,b) (!(a && b) && (a || b))\nfloat __agos_factorial(int);\n");
					fprintf(file,"double _agos_max(double*,int);\n");
					fprintf(file,"double _agos_min(double*,int);\n");
					fprintf(file, "#include <omp.h>\n#include \"Stopwatch.h\"\n");
					
					//alterado por cadim, 11/09/2009
					//colocando o openmp
					fprintf(file," typedef struct str_funcao{\n\tvoid (*funcao)(Solveode*);\n}typ_funcao;\n\n");
					
					prefixo="classe->";
					//imprime as funções auxiliares
					
					cur = rewind_list(algvarlist);
					int numCalcs =0;
					while(cur != NULL)
					{
					    fprintf(file,"\nvoid _f_calc_%s_(%s *classe){\n",cur->token.content,classname);
					    
					    // 		    file, classname, alglist, resolved_dep_list);
					    // 		   (FILE *file, char *classname, AlgList *list, TokenNode *orderedlist)
					    TokenNode *curlCalc = rewind_list(resolved_dep_list);
					    TokenNode *curCalc = NULL;
					    AlgList *curalgCalc = NULL;
					    int cont=0;
					    while(curlCalc != NULL)
					    {
						if(cont==numCalcs){
						    curalgCalc = rewind_list(alglist);
						    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content))
							curalgCalc = curalgCalc->next;			    
						    fprintf(file,"\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
						    curCalc = curalgCalc->eq;
						    curCalc = curCalc->next->next;
						    print_eq(file, curCalc);
						    fprintf(file,";\t//%d",curalgCalc->number);
						    curlCalc = curlCalc->next;
						    break;
						}else{
						    curCalc = curCalc->next->next;
						    curlCalc = curlCalc->next;
						}
						cont++;
					    }
					    
					    
					    //final da funcao
					    fprintf(file,"\n}\n");
					    cur = cur->next;
					    numCalcs++;
					}
					
					
					
					fprintf(file,"\n\t//DEPENDENT VARIABLES\n");
					cur = rewind_list(difvarlist);
					int numEdos=0;
					while(cur != NULL)
					{
					    fprintf(file,"\nvoid _f_%s_(%s *classe){\n",cur->token.content,classname);
					    DiffList *curlDiff = rewind_list(difflist);
					    TokenNode *curDiff = NULL;
					    int cont=0;
					    while(curlDiff != NULL)
					    {
						if(cont==numEdos){
						    curDiff = curlDiff->diffheader->eq->next;
						    fprintf(file,"\t%s%s_lado_direito_= ",prefixo, curlDiff->diffheader->diffvar.content);
						    print_eq(file, curDiff);
						    fprintf(file,";\n");
						    fprintf(file,"\t%s%s_new_= %sd%s*(%s%s_lado_direito_)+%s%s_old_;\t// %d\n",prefixo, curlDiff->diffheader->diffvar.content, prefixo, curlDiff->diffheader->freevar.content,prefixo, curlDiff->diffheader->diffvar.content,prefixo, curlDiff->diffheader->diffvar.content, curlDiff->number);
						    fprintf(file,"\t//%s%s_old_ = %s%s_new_;\n",prefixo, curlDiff->diffheader->diffvar.content, prefixo, curlDiff->diffheader->diffvar.content);
						    curlDiff = curlDiff->next;
						    break;
						}else{
						    curDiff = curlDiff->diffheader->eq->next;
						    curlDiff = curlDiff->next;
						}
						cont++;
						
					    }
					    fprintf(file,"\n}\n");
					    cur = cur->next;
					    numEdos++;
					}
					prefixo="";
					fprintf(file,"#define numEDO %d\n",numEdos);
					fprintf(file,"#define numAux %d\n",numCalcs);
					
					fprintf(file,"typ_funcao edos[numEDO];\n");
					fprintf(file,"typ_funcao auxiliares[numAux];\n");
					
					//fim -trecho
					fprintf(file,"\n\t//Constructor, initializes all variables with 0.0 or default CellML initial values\n");
					fprintf(file, "\t%s::%s()\n\t{\n", classname,classname);
					//alterado por cadim
					//11/09/2009
					//openmp
					
					cur = rewind_list(algvarlist);
					numCalcs=0;
					while(cur != NULL){
					    fprintf(file,"\t\tauxiliares[%d].funcao = _f_calc_%s_;\n", numCalcs, cur->token.content);
					    numCalcs++;
					    cur = cur->next;
					}
					
					cur = rewind_list(difvarlist);
					numEdos=0;
					while(cur != NULL){
					    fprintf(file,"\t\tedos[%d].funcao = _f_%s_;\n", numEdos, cur->token.content);
					    // 			    fprintf(file,"auxiliares[%d].funcao = _f_%s_;\n", numCalcs, cur->token.content);
					    numEdos++;
					    // 			    numCalcs++;
					    cur = cur->next;
					}
					
					cur = rewind_list(parvarlist); 		
					
					while(cur != NULL)
					{
					    fprintf(file,"\t\t%s = %.8e;\n", cur->token.content, cur->initialvalue);
					    //fprintf(file,"\t\t%s_new = %e;\n", cur->name, cur->value);
					    cur = cur->next;
					}
					
					fprintf(file,"\t\td%s = 0.0; %s_vec__ = NULL;\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
					
					cur = rewind_list(difvarlist);
					while(cur != NULL)
					{
					    fprintf(file,"\t\t%s = NULL;\n",cur->token.content);
					    fprintf(file,"\t\t%s_ini_ = %.8e;\n",cur->token.content, cur->initialvalue);
					    cur = cur->next;
					}
					fprintf(file,"\t\tit_countx = 0;\n");
					fprintf(file,"\t}\n");
					
					/****** FIM CONSTRUTOR *************************/
					
					
					fprintf(file, "\t%s::~%s()\n\t{\n", classname,classname);
					
					cur = rewind_list(difvarlist);
					while(cur != NULL){
					    fprintf(file, "\t\tif(%s != NULL) free(%s);\n", cur->token.content, cur->token.content);
					    cur = cur->next;
					    
					}
					fprintf(file, "\t}\n");
					
					
					/****** METODOS SET ****************************/
					fprintf(file,"\n\tint %s::setVariables(int indVariable, double value_new)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
					
					cur = rewind_list(difvarlist);
					int numVar = 0;
					while(cur != NULL)
					{
					    fprintf(file,"\t\tcase %d:\t\t%s_ini_ %s_old_= value_new;    break;\n", numVar, cur->token.content,cur->token.content);
					    numVar++;
					    cur = cur->next;
					}
					fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t\treturn 0;\n\t}\n");
					
					fprintf(file,"\n\tint %s::setParameters(int indVariable, double value_new)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
					
					cur = rewind_list(parvarlist);
					numVar = 0;
					while(cur != NULL)
					{
					    fprintf(file,"\t\tcase %d:\t\t%s = value_new;   break;\n", numVar, cur->token.content);
					    numVar++;
					    cur = cur->next;
					}
					fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t\treturn 0;\n\t}\n");
					
					
					fprintf(file,"\n\tint %s::setFreeVariable(double value_new)\n\t{\n\t\t",classname);
					fprintf(file,"d%s = value_new;\n\t}\n", difflist->diffheader->freevar.content);	
					
					/******* FIM SET ***********************/
					
					/******* METODOS GET ********************/
					fprintf(file,"\n\t//Get Methods\n");
					
					fprintf(file,"\n\tdouble %s::getVariables(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
					
					cur = rewind_list(difvarlist);
					numVar = 0;
					while(cur != NULL)
					{
					    fprintf(file,"\t\tcase %d:\t\treturn %s_old_;    break;\n", numVar, cur->token.content);
					    numVar++;
					    cur = cur->next;
					}
					fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t}\n");
					//////////////////
					fprintf(file,"\n\tdouble %s::getLadoDireito(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
					
					cur = rewind_list(difvarlist);
					numVar = 0;
					while(cur != NULL)
					{
					    fprintf(file,"\t\tcase %d:\t\treturn %s_lado_direito_;    break;\n", numVar, cur->token.content);
					    numVar++;
					    cur = cur->next;
					}
					fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t}\n");
					
					fprintf(file,"\n\tdouble %s::getParameters(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
					cur = rewind_list(parvarlist);
					numVar = 0;
					while(cur != NULL)
					{
					    fprintf(file,"\t\tcase %d:\t\treturn %s;    break;\n", numVar, cur->token.content);
					    numVar++;
					    cur = cur->next;
					}
					fprintf(file,"\t\tdefault:\tbreak;\n\t\t}\n\t}\n");
					
					fprintf(file,"\n\tdouble %s::getFreeVariable()\n\t{\n",classname);
					
					fprintf(file,"\t\treturn d%s;\n\t}\n", difflist->diffheader->freevar.content);
					
					
					fprintf(file,"\n\t//Get Methods - Variables\n\n");
					fprintf(file, "\tVariables %s::get_Variables()\n\t{\n\t\tVariables v(\"",classname);
					
					cur = rewind_list(difvarlist);
					while(cur != NULL)
					{
					    fprintf(file,"|%s#",cur->token.content);
					    cur = cur->next;
					}
					fprintf(file, "\");\n\t\treturn v;\n\t}\n");
					
					fprintf(file, "\tVariables %s::get_Parameters()\n\t{\n\t\tVariables v(\"",classname);
					
					cur = rewind_list(parvarlist);
					while(cur != NULL)
					{
					    fprintf(file,"|%s#",cur->token.content);
					    cur = cur->next;
					}
					fprintf(file, "\");\n\t\treturn v;\n\t}\n");	
					
					fprintf(file, "\tVariables %s::get_FreeVariable()\n\t{\n\t\tVariables v(\"",classname);
					
					fprintf(file,"|%s#",difflist->diffheader->freevar.content); 		
					fprintf(file, "\");\n\t\treturn v;\n\t}\n");	
					
					
					/********SET DO ARQUIVO **********************/
					fprintf(file,"\n\tvoid %s::setParametersFromFile(char *filename)\n\t{\n",classname);
					fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
					fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setParametersFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
					fprintf(file,"\tdouble value;\n\t\tint k = 0;\n\t\tVariables v = get_Parameters();\n\t\tint s = v.getQuantity();\n\t");
					fprintf(file,"\tfor(;k<s;k++)\n\t\t{\n\t\t");
					fprintf(file,"\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetParameters(k, value);\n\t\t}\n\t\tfclose(file);\n\t}\n");
					
					fprintf(file,"\n\tvoid %s::setVariablesFromFile(char *filename)\n\t{\n", classname);
					fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
					fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setVariablesFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
					fprintf(file,"\tdouble value;\n\t\tint k = 0;\n\t\tVariables v = get_Variables();\n\t\tint s = v.getQuantity();\n\t");
					fprintf(file,"\tfor(;k<s;k++)\n\t\t{\n\t\t");
					fprintf(file,"\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetVariables(k, value);\n\t\t}\n\t\tfclose(file);\n\t}\n");
					
					fprintf(file,"\n\tvoid %s::setFreeVariableFromFile(char *filename)\n\t{\n", classname);
					fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
					fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setFreeVariableFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
					fprintf(file,"\tdouble value;\n\t\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetFreeVariable(value);\n\t\tfclose(file);\n\t}\n");
					/********FIM SET DO ARQUIVO ******************/
					
					/******* METODO SOLVE ******************/
					fprintf(file,"\tdouble %s::solve(int firstcall__, int num_iterations__, int num_results__, int numThreads)\n\t{\n\n",classname);
					//fprintf(file,"//4 testeeeeeeeeeeeeee");
					
					fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
					fprintf(file,"\t\tstatic int num_results_bak = 0;\n");
					fprintf(file,"\t\tstatic int offset_step = 1;\n");
					fprintf(file,"\t\tif(firstcall__){\n");
					fprintf(file,"\t\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
					
					fprintf(file,"\t\t\tif(num_results__ <= 0)\n\t\t\t\tnum_results__ = 1;\n\t\t\tif(num_iterations__ <= 0)\n\t\t\t\tnum_iterations__ = 1;\n");
					fprintf(file,"\t\t\toffset_step = num_iterations__ / num_results__;\n");
					
					fprintf(file,"\t\t\tif(%s_vec__ != NULL)free( %s_vec__);\n\t\t\t%s_vec__ = (double *)malloc(sizeof(double)*num_results__);\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
					cur = rewind_list(difvarlist);
					while(cur != NULL)
					{
					    fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
					    fprintf(file,"\t\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*num_results__);\n", cur->token.content,cur->token.content,cur->token.content);
					    cur = cur->next;
					}
					fprintf(file,"\t\t\tnum_results_bak = num_results__;\n");
					fprintf(file,"\t\t\tnum_iterations_bak = num_iterations__;\n");
					fprintf(file,"\t\t}\n");
					
					fprintf(file,"\t\tint counter_it__ = 1, it_countx = 1;\n");
					fprintf(file,"\t\tint aux = num_iterations_bak%%num_results_bak;\n");
					//fprintf(file,"\t\tomp_set_num_threads(numThreads);\n");
					//fprintf(file,"\t\t#pragma omp parallel for ordered schedule(static)\n");
					
					fprintf(file,"\t\tfor(it_countx = 1; it_countx<=num_iterations_bak; it_countx++){\n");
					fprintf(file,"\t\t\t%s_new += d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
					
					// 08/03/2008//
					//        print_right_alg(file, classname, alglist, resolved_dep_list);
					///////////////
					//        print_diff(file, difflist);
					
					fprintf(file,"\t\t\tfor(int k=0;k<numAux;k++)\t");
					fprintf(file,"auxiliares[k].funcao(this);\n");
					
					fprintf(file,"\t\t\tfor(int k=0;k<numEDO;k++)\t");
					fprintf(file,"edos[k].funcao(this);\n");
					
					
					//fprintf(file,"\n\t\t\t#pragma omp critical\n");
					fprintf(file,"\t\t\tif(it_countx != aux && (it_countx-aux)%%offset_step == 0){\n");
					cur = rewind_list(difvarlist);
					while(cur != NULL)
					{
					    fprintf(file,"\t\t\t\t%s[counter_it__] = %s_new_;\n", cur->token.content,cur->token.content);
					    cur = cur->next;
					}
					fprintf(file, "\t\t\t\t%s_vec__[counter_it__] = %s_new;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
					fprintf(file,"\t\t\t\tcounter_it__++;\n");
					fprintf(file,"\t\t\t}\n");
					cur = rewind_list(difvarlist);
					while(cur != NULL)
					{
					    fprintf(file,"\t\t\t%s_old_ = %s_new_;\n", cur->token.content,cur->token.content);
					    cur = cur->next;
					}
					fprintf(file,"\t\t}\n");
					
					fprintf(file,"\t\t\t//FILE *fileptr;\n");
					fprintf(file,"\t\t\t//char filename[12];\n");
					fprintf(file,"\t\t\t//sprintf(filename,\"%%dthread.dat\",numThreads);\n");
					fprintf(file,"\t\t\t//fileptr = fopen(filename, \"wb\");\n");
					fprintf(file,"\t\t\t//for(int i =0;i<num_results__;i++){\n");
					fprintf(file,"\t\t\t//    fprintf(fileptr,\"%%f %%f\\n\", time_vec__[i], V[i]);\n");
					fprintf(file,"\t\t\t//}\n");
					fprintf(file,"\t\t\t//fclose(fileptr);\n");
					// 		fprintf(file,"\t\t\t \n",);
					fprintf(file, "\t\treturn (num_iterations_bak%%offset_step)*d%s;\n", difflist->diffheader->freevar.content);
					fprintf(file,"\t}\n");
					/****** FIM SOLUCAO ******************/
					
					fprintf(file, "\n\tdouble* %s::getIndependentVar()\n\t{\n", classname);
					fprintf(file, "\t\treturn %s_vec__;\n\t}\n", difflist->diffheader->freevar.content);
					
					
					fprintf(file,"\n\tdouble* %s::getSolution(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
					/*** GET SOLUTION ********************/
					
					cur = rewind_list(difvarlist);
					numVar = 0;
					while(cur != NULL)
					{
					    fprintf(file,"\t\tcase %d:\t\treturn %s;    break;\n", numVar, cur->token.content);
					    numVar++;		
					    cur = cur->next;
					}
					
					fprintf(file,"\t\tdefault:\treturn NULL;    break;\n\t\t}\n\t}\n");
					/*** FIM GET SOLUTION ****************/
					
					/*** SOLUCAO EM DISCO ****************/
					fprintf(file,"\tdouble %s::solveToFile(char *filename, char *fileaccess, int firstcall__, int num_iterations__, int num_results__)\n\t{\n\n",classname);
					
					fprintf(file,"\t\tstatic int num_iterations_bak = 0;\n");
					fprintf(file,"\t\tstatic int num_results_bak = 0;\n");
					fprintf(file,"\t\tstatic int offset_step = 1;\n");
					fprintf(file,"\t\tstatic char *fileaccess_bak = \"\";\n");
					fprintf(file,"\t\tif(firstcall__){\n");
					fprintf(file,"\t\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
					
					fprintf(file,"\t\t\tif(num_results__ <= 0)\n\t\t\t\tnum_results__ = 1;\n\t\t\tif(num_iterations__ <= 0)\n\t\t\t\tnum_iterations__ = 1;\n");
					fprintf(file,"\t\t\toffset_step = num_iterations__ / num_results__;\n");
					
					cur = rewind_list(difvarlist);
					while(cur != NULL)
					{
					    fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
					    fprintf(file,"\t\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*num_results__);\n", cur->token.content,cur->token.content,cur->token.content);
					    cur = cur->next;
					}
					fprintf(file,"\t\t\tnum_results_bak = num_results__;\n");
					fprintf(file,"\t\t\tnum_iterations_bak = num_iterations__;\n");
					fprintf(file,"\t\t\tfileaccess_bak = fileaccess;\n");
					fprintf(file,"\t\t}\n");
					
					fprintf(file,"\t\tFILE *file = fopen(filename, fileaccess_bak);\n");
					fprintf(file,"\t\tif(!file){\n\t\t\tfprintf(stderr,\"ERROR - solveToFile - Unable to open file %%s\\n\",filename);\n\t\t\texit(1);\n\t\t}\n");
					
					fprintf(file,"\t\tint counter_it__ = 1, it_countx = 1;\n");
					fprintf(file,"\t\tint aux = num_iterations_bak%%num_results_bak;\n");
					fprintf(file,"\t\tfor(it_countx = 1; it_countx<=num_iterations_bak; it_countx++){\n");
					fprintf(file,"\t\t\t%s_new += d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
					
					// 08/03/2008//
					print_right_alg(file, classname, alglist, resolved_dep_list);
					///////////////
					print_diff(file, difflist);
					
					fprintf(file,"\n\n\t\t\tif(it_countx != aux && (it_countx-aux)%%offset_step == 0){\n");
					fprintf(file,"\t\t\t\tfprintf(file,\"%%.8e");
					int list_size = get_list_size(difvarlist);
					for(int i = 0; i < list_size; i++)
					    fprintf(file," %%.8e");
					fprintf(file,"\\n\",%s_new", difflist->diffheader->freevar.content);
					
					cur = rewind_list(difvarlist);
					while(cur != NULL)
					{
					    fprintf(file,",%s_new_", cur->token.content);
					    cur = cur->next;
					}
					fprintf(file,");\n");
					fprintf(file,"\t\t\t\tcounter_it__++;\n");
					fprintf(file,"\t\t\t}\n");
					cur = rewind_list(difvarlist);
					while(cur != NULL)
					{
					    fprintf(file,"\t\t\t%s_old_ = %s_new_;\n", cur->token.content,cur->token.content);
					    cur = cur->next;
					}
					fprintf(file,"\t\t}\n");
					fprintf(file,"\t\tfclose(file);\n");
					fprintf(file, "\t\treturn (num_iterations_bak%%offset_step)*d%s;\n", difflist->diffheader->freevar.content);
					fprintf(file,"\t}\n");
					
					/****** FIM SOLUCAO ******************/
					
					
					
					
					/******* FIM SOLUCAO EM DISCO ****************/
					// 08/03/2008
					//print_alg(file, classname, alglist);
					print_if(file, classname, iflist);
					
					
					fprintf(file,"\n\nSolver *CreateObject(int w){\n");
					fprintf(file,"\t%s *solve_ode = new %s[w];\n", classname, classname);
					fprintf(file,"\treturn ((Solver*)solve_ode);\n");
					fprintf(file,"}\n");
					
					fprintf(file,"\n\nSolver **CreateObject2D(int w, int h){\n");
					fprintf(file,"\t%s **solve_ode = new %s*[h];\n", classname, classname);
					fprintf(file,"\tfor(int i=0;i<h; i++)\n");
					fprintf(file,"\t\tsolve_ode[i] = new %s[w];\n",classname);
					fprintf(file,"\treturn ((Solver**)solve_ode);\n");
					fprintf(file,"}\n");
					
					fprintf(file,"int ode_sizeof(){\n");
					fprintf(file,"\treturn sizeof(%s);\n", classname);
					fprintf(file,"}\n");
					
					fprintf(file,"\n\nfloat __agos_factorial(int f){\n\tif(f>=0 & f<2)\n\t\treturn 1.0;\n\telse if(f < 0)\n\t\treturn 0.0/0.0;\n\tfor(int i=f-1; i>=2; i--)\n\t\tf *= i;\n\treturn (float)f;\n}\n");
					fprintf(file,"double _agos_max(double* vector, int size){\n");
					fprintf(file,"\tdouble max =vector[0];\n");
					fprintf(file,"\tint i;\n");
					fprintf(file,"\tfor(i=1;i<size;i++){\n");
					fprintf(file,"\t\tif(vector[i]>max) max = vector[i];\n");
					fprintf(file,"\t}\n");
					fprintf(file,"\treturn max;\n}\n");
					
					fprintf(file,"double _agos_min(double* vector, int size){\n");
					fprintf(file,"\tdouble min = vector[0];\n");
					fprintf(file,"\tint i;\n");
					
					fprintf(file,"\tfor(i=1;i<size;i++){\n");
					fprintf(file,"\t\tif(vector[i]<min) min = vector[i];\n");
					fprintf(file,"\t}\n");
					fprintf(file,"\treturn min;\n}\n");
					fclose(file);
					return 0;
					}
					
					
					// 08/03/2008 ///////
					int print_right_alg(FILE *file, char *classname, AlgList *list, TokenNode *orderedlist)
					{
					    if(file == NULL)
					    {
						printf("ERROR - Can't write in file, print_alg");
						exit(1);
					    }
					    TokenNode *curl = rewind_list(orderedlist);
					    TokenNode *cur = NULL;
					    AlgList *curalg = NULL;
					    
					    while(curl != NULL)
					    {
						curalg = rewind_list(list);
						while(strcmp(curalg->eq->token.content, curl->token.content))
						    curalg = curalg->next;
						
						fprintf(file,"\n\t\t\tcalc_%s = ",curalg->eq->token.content);
						cur = curalg->eq;
						cur = cur->next->next;
						print_eq(file, cur);
						fprintf(file,";\t//%d",curalg->number);
						curl = curl->next;
					    }
					    fprintf(file,"\n");
					    return 0;
					}
					
					
					
					
					
					
					
					
					
					
					
					
					
					
	///TODO novo agos sundials pycml omp adaptativo
	int compile_Sundials_API(char *classname)
	{
	    
	    if(eq_counter <= 0)
	    {
		printf("ERROR - Equation not found\n");
		exit(1);
	    }
	    
	    
	    // teste da insercao de precedencia 08/03/2008
	    preced_alg_list = create_preced_alg_list(rewind_list(alglist));
	    AlgList *cural = rewind_list(preced_alg_list);
	    TokenNode *resolved_dep_list = NULL;
	    
	    TokenNode *cur = rewind_list(algvarlist);
	    
	    
	    
	    while(preced_alg_list != NULL)
	    {
		cural = rewind_list(preced_alg_list);
		while(cural != NULL){
		    //	printf("%s $$$$$\n", cural->eq->token.content);
		    cural = cural->next;
		}
		cural = rewind_list(preced_alg_list);
		
		//	printf("ASDFASDFASDFAS\n");
		while(cural != NULL){
		    
		    TokenNode *cureq = cural->eq->next;
		    //	printf("%s \n", cural->eq->token.content);
		    //				getchar();
		    if(cureq == NULL){
			resolved_dep_list = add_list(cural->eq->token, resolved_dep_list);
			preced_alg_list = delete_from_list(preced_alg_list, cural->eq->token);
			cural = cural->next;
			//					
			//					printf("Lista %s\n", resolved_dep_list->token.content);
			//					getchar();
		    }else{
			
			while(cureq !=NULL){
			    if(list_has_var(cureq->token, resolved_dep_list)){
				//printf("%s \n", cureq->token.content);
				//getchar();
				if(cureq->prev != NULL){
				    cureq->prev->next = cureq->next;
				}
				if(cureq->next != NULL){
				    cureq->next->prev = cureq->prev;
				}
			    }
			    cureq = cureq->next;
			}
			cural = cural->next;					
		    }	
		    
		}
		
	    }
	    //		TokenNode *currl = rewind_list(resolved_dep_list);
	    //					while(currl != NULL){
	    //						printf("LRD %s\n", currl->token.content);
	    //						currl = currl->next;
	    //		}
	    // fim - teste da insercao de precedencia 08/03/2008
	    FILE *file;
	    
	    char *filename = (char *)calloc(strlen(classname)+5, sizeof(char*));
	    
	    sprintf(filename,"%s.hpp",(const char*)classname);
	    
	    file = fopen(filename, "wb");
	    
	    fprintf(file, "#include \"sundials_types.h\"\n#include \"cvode.h\"\n#include \"cvode_dense.h\"\n#include \"nvector_serial.h\"\n#include \"sundials_dense.h\"\n");
	    fprintf(file, "#include \"cvode_diag.h\"\n");
	    fprintf(file, "#include \"cvode_band.h\"\n");
	    
	    fprintf(file, "#define ABSTOL 1.0E-06\n");
	    fprintf(file, "#define RELTOL 1.0E-04\n");
	    
	    fprintf(file, "static int check_flag(void *flagvalue, char *funcname, int opt);\n");
	    fprintf(file, "static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__);\n");
		
		
	    fprintf(file, "#include <stdio.h>\n#include <math.h>\n#include \"MCutil.hpp\"\n");
	    fprintf(file, "#define _H_FACTOR_ 2\n");
	    fprintf(file, "#define RECURSION_LIMIT 10\n");
	    
	    fprintf(file, "#define _EULER_ 0\n");
	    fprintf(file, " #define _ADAP_DT_ 4\n");
	    fprintf(file, "#define _RK2_ 3\n");
	    fprintf(file, "#define _ADAP_DT2_ 5\n");
	    
	    
	    fprintf(file, "#include <omp.h>\n");
	    fprintf(file, "#define _LOWER_LIMIT_   0.2\n");
	    fprintf(file, "#define _UPPER_LIMIT_   1.0\n");
	    
	    
	    int numCalcs =0;
	    while(cur != NULL)
	    {
		cur = cur->next;
		numCalcs++;
	    }
		
	    cur = rewind_list(difvarlist);
	    int numEdos=0;
	    while(cur != NULL)
	    {
		cur = cur->next;
		numEdos++;
	    }
	    
	    int totalSize = numCalcs + numEdos;
	    
	    //grafo que armazena as dependencias de cada equação
	    _node_ *grafoDependencias = (_node_*) malloc( (totalSize) * sizeof (_node_) );
	    
	    cur = rewind_list(algvarlist);
	    numCalcs =0;
	    while(cur != NULL)
	    {
		TokenNode *curlCalc = rewind_list(resolved_dep_list);
		TokenNode *curCalc = NULL;
		AlgList *curalgCalc = NULL;
		int cont=0;
		while(curlCalc != NULL)
		{
		    if(cont==numCalcs){
			curalgCalc = rewind_list(alglist);
			while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
			    curalgCalc = curalgCalc->next;
			}
			// 		printf("%s ->", curalgCalc->eq->token.content);
			grafoDependencias[numCalcs].equationName = curalgCalc->eq->token.content; 
			grafoDependencias[numCalcs].root = NULL;
			grafoDependencias[numCalcs].treeId = 0;
			grafoDependencias[numCalcs].flag = 0;
			grafoDependencias[numCalcs].height = -1;
			curCalc = curalgCalc->eq;
			curCalc = curCalc->next->next;
			
			criaListaPrecedencia(curCalc,grafoDependencias, numCalcs, totalSize);
			
			curlCalc = curlCalc->next;
			break;
		    }else{
			curCalc = curCalc->next->next;
			curlCalc = curlCalc->next;
		    }
		    cont++;
		}
		cur = cur->next;
		numCalcs++;
		// 	printf("\n");
	    }

	    cur = rewind_list(difvarlist);
	    numEdos=0;
	    while(cur != NULL)
	    {
		DiffList *curlDiff = rewind_list(difflist);
		TokenNode *curDiff = NULL;
		int cont=0;
		while(curlDiff != NULL)
		{
		    if(cont==numEdos){
			curDiff = curlDiff->diffheader->eq->next;
			// 		printf("::::\t\t%s \n", curlDiff->diffheader->diffvar.content);
			grafoDependencias[numCalcs + numEdos].equationName = curlDiff->diffheader->diffvar.content;
			grafoDependencias[numCalcs + numEdos].root = NULL;
			grafoDependencias[numCalcs + numEdos].treeId = 0;
			grafoDependencias[numCalcs + numEdos].flag = 0;
			grafoDependencias[numCalcs + numEdos].height = -1;
			criaListaPrecedencia(curDiff, grafoDependencias, numCalcs + numEdos, totalSize);
			curlDiff = curlDiff->next;
			break;
		    }else{
			curDiff = curlDiff->diffheader->eq->next;
			curlDiff = curlDiff->next;
		    }
		    cont++;
		}
		cur = cur->next;
		numEdos++;
		// 	printf("\n");
	    }
	    //     	printf("%d\n", totalSize);
	    int treeNumber = classifyTrees(grafoDependencias, totalSize);
	    //     printf("treeNumber: %d\n", treeNumber);
	    	    printGrafo(grafoDependencias, totalSize);
	    
	    
	    
	    fprintf(file,"/// equationType\n");
	    fprintf(file,"#define eqAgos 1\n");
	    fprintf(file,"#define eqPyCmlSimples 2\n");
	    fprintf(file,"#define eqPyCmlPE 3\n");
	    fprintf(file,"#define eqPyCmlLUT 4\n");
	    fprintf(file,"#define eqPyCmlPE_LUT 5\n");
	    fprintf(file,"#define forestSize %d\n", treeNumber);
	    
	    fprintf(file,"#include <sys/time.h>\n", treeNumber);
	    fprintf(file,"#include \"greedy.cpp\"\n", treeNumber);
	    
	    fprintf(file,"#define _INCREASE_DT_ 1.5\n");
	    fprintf(file,"#define _DECREASE_DT_ 0.65\n");
	    fprintf(file,"#define _DECREASE_DT_2_ 2\n");
	    
	    fprintf(file,"#define numEDO %d\n",numEdos);
	    fprintf(file,"#define numAux %d\n",numCalcs);
	    
	    fprintf(file, "\nclass %s\n{\npublic:\n", classname);
	    fputs("\tint it_countx;\n",file);
	    fputs("\tint *tree_thread;\n",file);
	    
	    /******* DECLARACAO DAS VARIAVEIS *******/
	    // building parvaslist
	    
	    // 	    TokenNode *cur = rewind_list(parvarlist);
	    cur = rewind_list(parvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\tdouble %s; \t // %s\n",cur->token.content, cur->units);
		//fprintf(file,"\tdouble %s_new;\n",cur->name);
		cur = cur->next;
	    }
	    
	    cur = rewind_list(algvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\tdouble calc_%s; \t // %s\n",cur->token.content, cur->units);
		//fprintf(file,"\tdouble %s_new;\n",cur->name);
		cur = cur->next;
	    }
		
		
	    fprintf(file,"\tdouble d%s, *%s_vec__;\n",difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
	    fprintf(file,"\tdouble %s_new;\n",difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\n\t//functions variables\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\tdouble *%s;\n", cur->token.content);
		fprintf(file,"\tdouble %s_new_, %s_old_, %s_ini_, %s_lado_direito_;\n", cur->token.content,cur->token.content,cur->token.content, cur->token.content);
		cur = cur->next;
	    }
	    /******* FIM * DECLARACAO DAS VARIAVEIS *******/
	    
	    fprintf(file,"\npublic:\n");
	    fprintf(file, "\t%s(double,double,int);\n", classname);
	    fprintf(file, "\t~%s();\n", classname);
	    
		
	    ///---------------ADDT----------------------
	    fprintf(file, "\tdouble STEP_TOLERANCE_;\n");
	    fprintf(file, "\tvoid getRightHandSide();\n");
	    fprintf(file, "\tvoid getRightHandSideParallel();\n");    
	    fprintf(file, "\tvoid setMaxStep(double);\n");
	    fprintf(file, "\tint getErrorCode(double, double);\n");
	    fprintf(file, "\tdouble solveToFileOMP(double finalTime = 0, double savingRate = 0, int method=1, char *fileName=\"output.dat\", int nThreads=1, double tol=0.2);\n");
	    fprintf(file, "\tvoid explicitEulerStep();\n");
	    fprintf(file, "\tvoid rungekutta2_Parallel( );\n");
	    fprintf(file, "\tvoid rungekutta2_Single( );\n");
	    fprintf(file, "\tvoid adaptiveDt_Parallel();\n");
	    fprintf(file, "\tint adaptiveDt_Single();\n");
	    
	    fprintf(file, "\tint adaptiveDt();\n");
	    //fprintf(file, "\tint adaptiveTrapezoidalStep();\n");
	    fprintf(file, "\tdouble errorAux;\n");
	    
	    fprintf(file, "\tint count1; int count0;\n");
	    fprintf(file, "\tint count3; int count4;\n");
	    fprintf(file, "\tint recursion_counter;\n");
	    fprintf(file, "\tdouble maxStep;\n");
	    ///-------------------------------------
	    
	    // VARIAVEIS PARA O CVODE 
	    fprintf(file, "\t// CVODE VARIABLES\n");
	    fprintf(file,"\trealtype reltol__, abstol__;\n");
	    fprintf(file,"\tvoid *cvode_mem_cvode__;\n");
	    fprintf(file, "\tN_Vector dependent_variable__;\n");
	    fprintf(file, "\tint flag__, flagr__;\n");
	    fputs("\tdouble *depvar__;\n",file);
	    
	    fputs("\tint setVariables(int, double);\n",file);
	    fputs("\tint setParameters(int, double);\n",file);
	    fputs("\tint setFreeVariable(double);\n",file);
	    
	    
	    fputs("\tdouble getVariables(int);\n",file);
	    fputs("\tdouble getLadoDireito(int);\n",file);
	    fputs("\tvoid executeTree(void* tree, int treeSize);\n",file);
	    fputs("\tdouble getParameters(int);\n",file);
	    fputs("\tdouble getFreeVariable();\n",file);
	    
	    fputs("\tVariables get_Parameters();\n",file);
	    fputs("\tVariables get_Variables();\n",file);
	    fputs("\tVariables get_FreeVariable();\n",file);
	    
	    fputs("\tvoid setParametersFromFile(char*);\n", file);
	    fputs("\tvoid setVariablesFromFile(char*);\n", file);
	    fputs("\tvoid setFreeVariableFromFile(char*);\n", file);
	    
	    fputs("\tdouble solve(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0, int numThreads=1);\n",file);
	    fputs("\tdouble jacobian(double final, double svRate);\n",file);
	    fputs("\tdouble* getSolution(int indVariable);\n",file);
	    fputs("\tdouble* getIndependentVar();\n",file);
	    fputs("\tdouble solveToFile(char *filename, char *fileaccess = \"\", int firstcall__ = 0, int num_iterations__ = 0, int num_results__ = 0);\n",file);
	    fputs("public:\n",file);
	    
	    int iif;
	    for(iif =0; iif< ifs; iif++)
	    {
		fprintf(file, "\tinline double ifnumber_%d();\n", iif);
	    }
	    fprintf(file, "\tvoid reInitCVODE();\n");
	    fprintf(file, "\tvoid setCVODEMaxStep(double maxstep);\n");
	    fprintf(file, "\tvoid solveCVODE(int firstcall__ = 0, int steps__ = 0, int num_results__ =0, int method__=0,char *fileName=\"saida.out\", int cv_method__=1);\n");
	    
	    fprintf(file, "\tvoid solver(int method=0,double finalTime=0, double svRate=0, char* fileName=\"\", int threads=1);\n");
	    
	    
	    fprintf(file,"private:\n");
	    fprintf(file,"\tdouble errorAD_;\n");
	    fprintf(file,"\tdouble edos_euler_[numEDO];\n");
	    fprintf(file,"\tdouble edos_rk2_[numEDO];\n");
	    fprintf(file,"\tdouble edos_fn_[numEDO];\n");
	    fprintf(file,"\tdouble ladodir_fn_[numEDO];\n");
	    fprintf(file,"\tvoid save_step(FILE *fileptr, int method);\n");
	    
	    
	    fprintf(file,"\tvoid euler(double finalTime, FILE *file);\n");
	    fprintf(file,"\tvoid euler_OMP(double finalTime, FILE *fileptr, int nThreads);\n");
	    fprintf(file,"\tvoid rungeKutta2ndOrder(double finalTime, FILE *fileptr);\n");
	    fprintf(file,"\tvoid addt(double finalTime, FILE *fileptr);\n");
	    fprintf(file,"\tvoid addt2(double finalTime, FILE *fileptr);\n");
	    fprintf(file,"\tvoid addt_OMP(double finalTime, FILE *fileptr, int nThreads);\n");
	    fprintf(file,"\tvoid addt2_OMP(double finalTime, FILE *fileptr, int nThreads);\n");
	    
	    
	    
	    fprintf(file,"\tdouble timeSaving;\n");
	    fprintf(file,"\tdouble previous_dt;\n");
	    fprintf(file,"\tdouble savingRate;\n");
	    
	    fprintf(file,"};\n\n");
	    fprintf(file,"#define AGOS_NAN (0.0/0.0)\n#define AGOS_INF (1.0/0.0)\n");
	    fprintf(file,"#define __agos_xor(a,b) (!(a && b) && (a || b))\nfloat __agos_factorial(int);\n");
	    fprintf(file,"double _agos_max(double*,int);\n");
	    fprintf(file,"double _agos_min(double*,int);\n");
	    fprintf(file,"double _agos_round( double, int);\n");
	    fprintf(file," typedef struct str_funcao{\n\tvoid (*funcao)(Solveode*);\n}typ_funcao;\n\n");
	    fprintf(file,"typedef struct str_forest{\n");
	    fprintf(file,"typ_funcao* tree;\n");
	    fprintf(file,"int treeSize;\n");
	    fprintf(file,"}typ_forest;\n");
	    
	    
	    
	    /// /////////////////////
	    
	    /// ////////////////
	    prefixo="__AGOS->"; 
	    
	    
	    /*//TODO TESTE AGOS
	    prefixo="";
	    
	    fprintf(file,"\t\tconst double %s_new = __AGOS->%s_new;\n", difflist->diffheader->freevar.content, 
	    difflist->diffheader->freevar.content);
	    
	    cur = rewind_list( parvarlist); 		
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tconst double %s = %.10e;\n", cur->token.content, cur->initialvalue);
		// 		fprintf(file,"\t\tconst double %s = __AGOS->%s;\n", cur->token.content, cur->token.content);
		cur = cur->next;
	}
	cur = rewind_list(difvarlist);
	while(cur != NULL)
	{
	    fprintf(file,"\t\tconst double %s_old_= __AGOS->%s_old_;\n", cur->token.content, cur->token.content);
	    cur = cur->next;
	}
	//TODO FIM TESTE AGOS
	*/    
	    
	    /** ************** TODO
	    ********* imprime as equações, sem funções, separadas por arvores
	    *****/
	    for(int i=1;i<=treeNumber;i++){
		fprintf(file,"\nvoid __tree%d__( Solveode *__AGOS){\n", i);
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		
		
		int countTreeItems = 0;
		
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			//procurar por esta equação calc e depois pelo diff:  grafoDependencias[j].equationName);
			/// imprime os calc de uma arvore
			cur = rewind_list(algvarlist);
			numCalcs =0;
			while(cur != NULL)
			{
			    TokenNode *curlCalc = rewind_list(resolved_dep_list);
			    TokenNode *curCalc = NULL;
			    AlgList *curalgCalc = NULL;
			    int cont=0;
			    while(curlCalc != NULL)
			    {
				if(cont==numCalcs){
				    curalgCalc = rewind_list(alglist);
				    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
					curalgCalc = curalgCalc->next;
				    }
				    if(strcmp(grafoDependencias[j].equationName,curalgCalc->eq->token.content) ==0){
					fprintf(file,"\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
					curCalc = curalgCalc->eq;
					curCalc = curCalc->next->next;
					print_eq(file, curCalc, 0, 0);
					curlCalc = curlCalc->next;
					fprintf(file,";\n");
					break;
				    }
				    
				}else{
				    curCalc = curCalc->next->next;
				    curlCalc = curlCalc->next;
				}
				cont++;
			    }
			    
			    cur = cur->next;
			    numCalcs++;
			}
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t%s%s_lado_direito_= ",prefixo, curlDiff->diffheader->diffvar.content);
					print_eq(file, curDiff,0,0);
					fprintf(file,";\n");    
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"} //fim\n");
	    }
	    
	    
	    /// ///////////////
	    
	    /// /////////////////////
	    
	    
	    prefixo="__AGOS->"; 
	    fprintf(file,"\nvoid __AGOS_EQUATIONS__( Solveode *__AGOS){\n");
	    
	    //TODO TESTE AGOS
	    prefixo="";
	    
	    fprintf(file,"\t\tconst double %s_new = __AGOS->%s_new;\n", difflist->diffheader->freevar.content, 
		    difflist->diffheader->freevar.content);
	    
	    cur = rewind_list( parvarlist); 		
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tconst double %s = %.10e;\n", cur->token.content, cur->initialvalue);
// 		fprintf(file,"\t\tconst double %s = __AGOS->%s;\n", cur->token.content, cur->token.content);
		cur = cur->next;
	    }
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tconst double %s_old_= __AGOS->%s_old_;\n", cur->token.content, cur->token.content);
		cur = cur->next;
	    }
	    //TODO FIM TESTE AGOS
	    
	    
	    /** ************** TODO
	    ********* imprime as equações, sem funções, separadas por arvores
	    *****/
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		
		
		int countTreeItems = 0;
		
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			//procurar por esta equação calc e depois pelo diff:  grafoDependencias[j].equationName);
			    /// imprime os calc de uma arvore
			    cur = rewind_list(algvarlist);
			    numCalcs =0;
			    while(cur != NULL)
			    {
				TokenNode *curlCalc = rewind_list(resolved_dep_list);
				TokenNode *curCalc = NULL;
				AlgList *curalgCalc = NULL;
				int cont=0;
				while(curlCalc != NULL)
				{
				    if(cont==numCalcs){
					curalgCalc = rewind_list(alglist);
					while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
					    curalgCalc = curalgCalc->next;
					}
					if(strcmp(grafoDependencias[j].equationName,curalgCalc->eq->token.content) ==0){
// 					    fprintf(file,"\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
					    //TODO TESTE AGOS
					    fprintf(file,"\tconst double calc_%s = ",curalgCalc->eq->token.content);
					    curCalc = curalgCalc->eq;
					    curCalc = curCalc->next->next;
					    print_eq(file, curCalc, 0, 0);
					    curlCalc = curlCalc->next;
					    fprintf(file,";\n");
					    break;
					}
					
				    }else{
					curCalc = curCalc->next->next;
					curlCalc = curlCalc->next;
				    }
				    cont++;
				}
				
				cur = cur->next;
				numCalcs++;
			    }
		    
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
// 					fprintf(file,"\t%s%s_lado_direito_= ",prefixo, curlDiff->diffheader->diffvar.content);
					//TODO TESTE AGOS
					fprintf(file,"\t%s%s_lado_direito_= ","__AGOS->", curlDiff->diffheader->diffvar.content);
					
					print_eq(file, curDiff,0,0);
					fprintf(file,";\n");    
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
		    
		    }
		}
		
	    }
	    fprintf(file,"} //fim\n");
	    prefixo = "";
	    
	    
	    //imprime as equações diferenciais, sem nenhuma avaliação parcial. 
	    //Todas as variáveis algébricas são escritas no código da equação diferencial;
	    ///imprime as equações diferenciais
	    cur = rewind_list(difvarlist);
	    numEdos=0;
	    prefixo="__AGOS->"; 
	    
	    while(cur != NULL)
	    {
		DiffList *curlDiff = rewind_list(difflist);
		TokenNode *curDiff = NULL;
		int cont=0;
		while(curlDiff != NULL)
		{
		    if(cont==numEdos){
			curDiff = curlDiff->diffheader->eq->next;
			fprintf(file,"\tdouble %s_completa(Solveode *__AGOS){\n", curlDiff->diffheader->diffvar.content);
// 			fprintf(file,"\t\t%s%s_lado_direito_= ",prefixo, curlDiff->diffheader->diffvar.content);
			fprintf(file,"\t\treturn  ");
			print_eq(file, curDiff,0,0, 0);
			fprintf(file,";\n");    
			fprintf(file,"\t}\n");
			curlDiff = curlDiff->next;
							
		    }else{
			curDiff = curlDiff->diffheader->eq->next;
			curlDiff = curlDiff->next;
		    }
		    cont++;
		    
		}
		cur = cur->next;
		numEdos++;
	    }
	    prefixo = "";
	    
	    /** **********Fim***********/
	    
	    fprintf(file,"typedef struct str__rightHandSideFunction{\n");
	    fprintf(file,"\tvoid (*function)(Solveode*);\n");
	    fprintf(file,"}typ_rightHandSideFunction;\n");
	    fprintf(file,"typ_rightHandSideFunction rightHandSideFunction;\n");
	    fprintf(file,"typ_rightHandSideFunction forest[%d];\n", treeNumber);
	    
	    prefixo="";
	    /*fprintf(file,"\tint forestSize ");
	    for(int i=1;i<=treeNumber;i++){
		fprintf(file,", tree%dSize", i);
	    }
	    fprintf(file,";\n");
	    */
	   /* fprintf(file,"\ttyp_forest* forest;\n");
	    for(int i=1;i<=treeNumber;i++){
		fprintf(file,"\ttyp_funcao* tree%d;\n", i);
	    }*/
	    
	    fprintf(file,"\t#include \"simple.hpp\"\n");
	    fprintf(file,"\t#include \"pe.hpp\"\n");
	    fprintf(file,"\t#include \"lut.hpp\"\n");
	    fprintf(file,"\t#include \"pe_lut.hpp\"\n");
	    
	    fprintf(file,"\n\t//Constructor, initializes all variables with 0.0 or default CellML initial values\n");
	    fprintf(file, "\t%s::%s(double abs, double rel, int equationType)\n\t{\n", classname,classname);

	    fprintf(file, "\t\tif(equationType==eqAgos){\n");
	    fprintf(file, "\t\t\trightHandSideFunction.function = __AGOS_EQUATIONS__;\n");
	    fprintf(file, "\t\t}else if(equationType==eqPyCmlSimples){\n");
	    fprintf(file, "\t\t\trightHandSideFunction.function = __equation__;\n");
	    fprintf(file, "\t\t}else if(equationType==eqPyCmlPE){\n");
	    fprintf(file, "\t\t\trightHandSideFunction.function = __equation__pe_;\n");
	    fprintf(file, "\t\t}else if(equationType==eqPyCmlLUT){\n");
	    fprintf(file, "\t\t\trightHandSideFunction.function = __equation__lut_;\n");
	    fprintf(file, "\t\t}else if(equationType==eqPyCmlPE_LUT){\n");
	    fprintf(file, "\t\t\trightHandSideFunction.function = __equation__pe_lut_;\n");
	    fprintf(file, "\t\t}\n");
	    
	    for(int i=0;i<treeNumber;i++){
		fprintf(file,"\t\tforest[%d].function = __tree%d__;\n", i, i+1);
	    }
    
	    cur = rewind_list( parvarlist); 		
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t%s = %.8e;\n", cur->token.content, cur->initialvalue);
		cur = cur->next;
	    }
	    
	    fprintf(file,"\t\td%s = 0.0; %s_vec__ = NULL;\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
	    
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t%s = NULL;\n",cur->token.content);
		fprintf(file,"\t\t%s_ini_ = %.8e;\n",cur->token.content, cur->initialvalue);
		
		cur = cur->next;
	    }
	    fprintf(file,"\t\tabstol__ = abs;\n");
	    fprintf(file,"\t\treltol__ = rel;\n");
	    
	    fprintf(file,"\t\tit_countx = 0;\n");
	    fprintf(file,"\t}\n");
	    
	    /****** FIM CONSTRUTOR *************************/
	    
	    fprintf(file, "\t%s::~%s()\n\t{\n", classname,classname);
	    
	    cur = rewind_list(difvarlist);
	    while(cur != NULL){
		fprintf(file, "\t\tif(%s != NULL) free(%s);\n", cur->token.content, cur->token.content);
		cur = cur->next;
		
	    }
	    fprintf(file, "\t}\n");		
	    
	    /****** METODOS SET ****************************/
	    fprintf(file,"\n\tint %s::setVariables(int indVariable, double value_new)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
	    
	    cur = rewind_list(difvarlist);
	    int numVar = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tcase %d:\t\t%s_ini_ = %s_old_= value_new;    break;\n", numVar, cur->token.content, cur->token.content);
		numVar++;
		cur = cur->next;
	    }
	    fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t\treturn 0;\n\t}\n");
	    
	    fprintf(file,"\n\tint %s::setParameters(int indVariable, double value_new)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
	    
	    cur = rewind_list(parvarlist);
	    numVar = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tcase %d:\t\t%s = value_new;   break;\n", numVar, cur->token.content);
		numVar++;
		cur = cur->next;
	    }
	    fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t\treturn 0;\n\t}\n");
	    
	    
	    fprintf(file,"\n\tint %s::setFreeVariable(double value_new)\n\t{\n\t\t",classname);
	    fprintf(file,"d%s = value_new;\n\t}\n", difflist->diffheader->freevar.content);
	    
	    /// addt
	    fprintf(file,"\n\tvoid %s::setMaxStep(double maxSt)\n\t{\n\t\t",classname);
	    fprintf(file,"this->maxStep=maxSt;\n\t}\n");	
	    
	    
	    
	    fprintf(file,"\tvoid %s::solver(int method, double finalTime, double svRate, char* fileName, int threads)\n",classname);
	    fprintf(file,"\t{\n");
	    fprintf(file,"\t\tthis->savingRate = svRate;\n");
	    fprintf(file,"\t\tif(finalTime <= 0){\n");
	    fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - solveToFile - negative finalTime is not allowed\\n\");\n");
	    fprintf(file,"\t\t\texit(1);\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tif(savingRate < 0){\n");
	    fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - solveToFile - negative saving rate is not allowed\\n\");\n");
	    fprintf(file,"\t\t\texit(1);\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tFILE *fileptr;\n");
	    fprintf(file,"\t\tif(savingRate!=0.0)\n");
	    fprintf(file,"\t\t\tfileptr = fopen(fileName, \"w+\");\n");
	    fprintf(file,"\t\tif(threads==1)\n");
	    fprintf(file,"\t\t{\n");
	    fprintf(file,"\t\t\tif(method==_EULER_){\n");
	    fprintf(file,"\t\t\t\tthis->euler(finalTime, fileptr);\n");
	    fprintf(file,"\t\t\t}else if(method==_RK2_){\n");
	    fprintf(file,"\t\t\t\tthis->rungeKutta2ndOrder(finalTime, fileptr);\n");
	    fprintf(file,"\t\t\t}else if(method == _ADAP_DT_){\n");
	    fprintf(file,"\t\t\t\tthis->addt(finalTime, fileptr);\n");
	    fprintf(file,"\t\t\t}else if(method == _ADAP_DT2_){\n");
	    fprintf(file,"\t\t\t\tthis->addt2(finalTime, fileptr);\n");
	    fprintf(file,"\t\t\t}else{\n");
	    fprintf(file,"\t\t\t\tprintf(\"Invalid option!\\n\");\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t}else//OMP\n");
	    fprintf(file,"\t\t{\n");
	    
	    fprintf(file,"\t\t\ttree_thread = (int*)malloc(sizeof(double)*forestSize);\n");
	    fprintf(file,"\t\t\tint *_jobs = (int*)malloc(sizeof(double)*forestSize);\n");
	    fprintf(file,"\t\t\tfor(int i=0; i< forestSize; i++){\n");
	    fprintf(file,"\t\t\t\ttimeval t1, t2;\n");
	    fprintf(file,"\t\t\t\tgettimeofday(&t1, NULL);\n");
	    fprintf(file,"\t\t\t\tfor (int j=0; j< 100; j++)\n");
	    fprintf(file,"\t\t\t\t\tforest[i].function(this); //invoca cada uma das arvores da floresta\n");
	    fprintf(file,"\t\t\t\t//stop timer\n");
	    fprintf(file,"\t\t\t\tgettimeofday(&t2, NULL);\n");
	    fprintf(file,"\t\t\t\t//compute and print the elapsed time in millisec\n");
	    fprintf(file,"\t\t\t\tdouble elapsedTime = (t2.tv_sec - t1.tv_sec);      // sec to ms\n");
	    fprintf(file,"\t\t\t\t//elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms\n");
	    fprintf(file,"\t\t\t\t_jobs[i] = (int)(t2.tv_usec - t1.tv_usec) + 1;\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t\tgreedyPartitionProblem(_jobs, forestSize, threads, tree_thread);\n");
	    
    	    
	    fprintf(file,"\t\t\tif(method==_EULER_){\n");
	    fprintf(file,"\t\t\t\tthis->euler_OMP(finalTime, fileptr, threads);\n");
	    fprintf(file,"\t\t\t}else if(method == _ADAP_DT_){\n");
	    fprintf(file,"\t\t\t\tthis->addt_OMP(finalTime, fileptr, threads);\n");
	    fprintf(file,"\t\t\t}else if(method == _ADAP_DT2_){\n");
	    fprintf(file,"\t\t\t\tthis->addt2_OMP(finalTime, fileptr, threads);\n");
	    fprintf(file,"\t\t\t}else{\n");
	    fprintf(file,"\t\t\t\tprintf(\"Invalid option!\\n\");\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tif(savingRate!=0.0)\n");
	    fprintf(file,"\t\t\tfclose(fileptr);\n");
	    fprintf(file,"\t}\n");
	    ////////
	    

	    fprintf(file,"\tvoid %s::save_step(FILE *file, int method){\n",classname);
	    fprintf(file,"\t\tdouble time_aux = this->time_new;\n");
	    fprintf(file,"\t\tdouble diff =  _agos_round(this->time_new - timeSaving, 10);\n");

	    fprintf(file,"\t\tif(diff==0){\n");
	    fprintf(file,"\t\t\tfprintf(file,  \"%%.8f %%.8f %%.8f\\n\", this->time_new, this->V_new_,previous_dt);\n");
	    fprintf(file,"\t\t\tthis->timeSaving += savingRate;\n");
	    fprintf(file,"\t\t}else if(diff > 0){\n");
	    fprintf(file,"\t\t\t//salva resultados que passaram do tempo de salvar\n");
	    fprintf(file,"\t\t\tdouble *edos_new_aux_;\n");
	    fprintf(file,"\t\t\tdouble *edos_old_aux_;\n");
	    fprintf(file,"\t\t\tedos_old_aux_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\t\tedos_new_aux_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\t\tfor(int i = 0; i < numEDO; i++){\n");
	    fprintf(file,"\t\t\t\tedos_old_aux_[i] = this->getVariables(i);\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t\tdouble _righthandside_[numEDO];\n");

	    fprintf(file,"\t\t\tdouble told = this->time_new - this->previous_dt;\n");
	    fprintf(file,"\t\t\t//encontra-se o dt adequado pra salvar na hora certa\n");
	    fprintf(file,"\t\t\tdouble dtTemp =   this->timeSaving - (told);\n");
	    fprintf(file,"\t\t\tthis->time_new = told;\n");

	    fprintf(file,"\t\t\t//verifica se a iteração acima alcançou o tempo real\n");
	    fprintf(file,"\t\t\tdouble euler_ladodireito[numEDO], rk2_res[numEDO], old_aux[numEDO];\n");
	    fprintf(file,"\t\t\twhile(time_aux >= timeSaving){\n");
	    fprintf(file,"\t\t\t\tthis->time_new += dtTemp;\n");
	    fprintf(file,"\t\t\t\t//calcula mais uma iteração\n");
	    fprintf(file,"\t\t\t\trightHandSideFunction.function(this);\n");
	    fprintf(file,"\t\t\t\t//euler\n");	    
	    fprintf(file,"\t\t\t\tfor(int i = 0; i < numEDO; i++){\n");
	    fprintf(file,"\t\t\t\t\told_aux[i] = this->getVariables(i);\n");
	    fprintf(file,"\t\t\t\t\tedos_new_aux_[i] = this->getVariables(i) + this->getLadoDireito(i)*dtTemp;\n");
	    fprintf(file,"\t\t\t\t\teuler_ladodireito[i] = this->getLadoDireito(i);\n");
	    fprintf(file,"\t\t\t\t\tthis->setVariables(i, edos_new_aux_[i]);\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t//se RK2\n");
	    fprintf(file,"\t\t\t\tif(method==_RK2_)\n");
	    fprintf(file,"\t\t\t\t{\n");
	    fprintf(file,"\t\t\t\t\tthis->time_new += dtTemp;\n");
	    fprintf(file,"\t\t\t\t\trightHandSideFunction.function(this);\n");
	    fprintf(file,"\t\t\t\t\t//rk2\n");
	    fprintf(file,"\t\t\t\t\tfor(int i = 0; i < numEDO; i++){\n");
	    fprintf(file,"\t\t\t\t\t\trk2_res[i] = old_aux[i] +(this->getLadoDireito(i) + euler_ladodireito[i])*(dtTemp/2);\n");
	    fprintf(file,"\t\t\t\t\t}		\n");
	    fprintf(file,"\t\t\t\t\tthis->time_new -= dtTemp;\n");
	    fprintf(file,"\t\t\t\t}\n");
	    
	    fprintf(file,"\t\t\t\tif(method==_RK2_ ){\n");
	    fprintf(file,"\t\t\t\t\tfor(int i = 0; i < numEDO; i++){\n");
	    fprintf(file,"\t\t\t\t\t\tthis->setVariables(i, rk2_res[i] );\n");
	    fprintf(file,"\t\t\t\t\t}   \n");
	    fprintf(file,"\t\t\t\t\tfprintf(file,\"%%.8e %%.8e %%.8e\\n\", time_new, rk2_res[0],previous_dt);\n");
	    fprintf(file,"\t\t\t\t}else{\n");
	    fprintf(file,"\t\t\t\t\tfor(int i = 0; i < numEDO; i++){\n");
	    fprintf(file,"\t\t\t\t\t\tthis->setVariables(i, edos_new_aux_[i] );\n");
	    fprintf(file,"\t\t\t\t\t}   \n");
	    fprintf(file,"\t\t\t\t\tfprintf(file,\"%%.8e %%.8e %%.8e\\n\", time_new, edos_new_aux_[0],previous_dt);\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\tthis->timeSaving += savingRate;\n");
	    fprintf(file,"\t\t\t\t//coloca o dtime igual ao savingRate para bater exatamente com a proxima iteração de salvar\n");
	    fprintf(file,"\t\t\t\tdtTemp = savingRate;\n");
	    fprintf(file,"\t\t\t}//fimwhile\n");
	    fprintf(file,"\t\t\t//volta com os valores old originais\n");
	    fprintf(file,"\t\t\tfor(int i = 0; i < numEDO; i++){\n");
	    fprintf(file,"\t\t\t\tthis->setVariables(i, edos_old_aux_[i] );\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t\tthis->time_new = time_aux;\n");
	    fprintf(file,"\t\t}//fim else if\n");
    	    fprintf(file,"\t}\n");

/// ***************	    


	fprintf(file,"\tvoid %s::euler(double finalTime, FILE *fileptr){\n",classname);
	    fprintf(file,"\t\tthis->previous_dt = this->dtime;\n");
	    fprintf(file,"\t\t//initializes the variables\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tthis->%s_new_ = this->%s_old_ = this->%s_ini_;\n", 
			cur->token.content, cur->token.content, cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file,"\t\tthis->%s_new = this->%s;\n", difflist->diffheader->freevar.content , difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->timeSaving = this->%s;\n", difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\t\tif(savingRate!=0.0)\n");
	    fprintf(file,"\t\t\tthis->save_step(fileptr, _EULER_);//save the initial conditions\n");
		
	    fprintf(file,"\t\twhile(this->%s_new<=finalTime){\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\tthis->%s_new	+= this->d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\trightHandSideFunction.function(this);\n");
	    
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\tthis->%s_new_ = this->d%s*(this->%s_lado_direito_) + this->%s_old_;\n", 
			cur->token.content, difflist->diffheader->freevar.content,cur->token.content,cur->token.content);
		cur = cur->next;
	    }
	    
	    fprintf(file,"\t\t\t//save results on a file\n");
	    fprintf(file,"\t\t\tif(savingRate!=0.0)\n");
	    fprintf(file,"\t\t\t{\n");
	    fprintf(file,"\t\t\t\tthis->save_step(fileptr, _EULER_);\n");
	    fprintf(file,"\t\t\t}\n");
	    
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tthis->%s_old_ = this->%s_new_;\n", 
			cur->token.content, cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t}\n");
	    
	    /// ************
	    fprintf(file,"\tvoid %s::rungeKutta2ndOrder(double finalTime, FILE *fileptr){\n",classname);
	    fprintf(file,"\t\tthis->previous_dt = this->dtime;\n");
	    fprintf(file,"\t\t//initializes the variables\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tthis->%s_new_ = this->%s_old_ = this->%s_ini_;\n", 
			cur->token.content, cur->token.content, cur->token.content);
			cur = cur->next;
	    }
	    fprintf(file,"\t\tthis->%s_new = this->%s;\n", difflist->diffheader->freevar.content ,
		    difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->timeSaving = this->%s;\n", difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\t\tif(savingRate!=0.0)\n");
	    fprintf(file,"\t\tthis->save_step(fileptr, _RK2_);//save the initial conditions\n");
	    fprintf(file,"\t\tdouble edos_old_aux_[numEDO];\n");
	    fprintf(file,"\t\tdouble edos_rightside_aux_[numEDO];\n");
	    fprintf(file,"\t\twhile(this->%s_new<=finalTime){\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\tthis->%s_new	+= this->d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\trightHandSideFunction.function(this);\n");
	    
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\tthis->%s_new_ = this->d%s*(this->%s_lado_direito_) + this->%s_old_;\n", 
			cur->token.content, difflist->diffheader->freevar.content,cur->token.content,cur->token.content);
			cur = cur->next;
	    }
	    fprintf(file,"\t\t\t//stores the old variables in a vector\n");
	    fprintf(file,"\t\t\tfor(int i=0;i<numEDO;i++){\n");
	    fprintf(file,"\t\t\t\tedos_old_aux_[i] = this->getVariables(i);\n");
	    fprintf(file,"\t\t\t\tedos_rightside_aux_[i] = this->getLadoDireito(i);\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t\t//steps one iteration ahead;\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\tthis->%s_old_ = this->%s_new_;\n", 
			cur->token.content, cur->token.content);
			cur = cur->next;
	    }
	    fprintf(file,"\t\t\tthis->%s_new	+= this->d%s;\n", 
		    difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\t\t\trightHandSideFunction.function(this);\n");
		    
	    fprintf(file,"\t\t\t//computes the runge kutta second order method\n");
		    
	    cur = rewind_list(difvarlist);
	    int count=0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\tthis->%s_new_ = ( this->%s_lado_direito_ + edos_rightside_aux_[%d] ) * this->d%s/2 + edos_old_aux_[%d];\n", 
			cur->token.content, cur->token.content, count, difflist->diffheader->freevar.content, count);
		cur = cur->next;
		count++;    
	    }
	    
	    fprintf(file,"\t\t\tthis->%s_new	-= this->d%s;//step back, to save in the right time\n", 
		    difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\tfor(int i=0;i<numEDO;i++){\n");
	    fprintf(file,"\t\t\t\tthis->setVariables(i, edos_old_aux_[i]);\n");
	    fprintf(file,"\t\t\t}\n");
	    
	    fprintf(file,"\t\t\t//save results on a file\n");
	    fprintf(file,"\t\t\tif(savingRate!=0.0)\n");
	    fprintf(file,"\t\t\t{\n");
	    fprintf(file,"\t\t\t\tthis->save_step(fileptr,_RK2_);\n");
	    fprintf(file,"\t\t\t}\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\tthis->%s_old_ = this->%s_new_;\n", 
			cur->token.content, cur->token.content);
			cur = cur->next;
	    }
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t}\n");
	    
	    /// ///////********
	   
	   fprintf(file,"\tvoid %s::addt(double finalTime, FILE *fileptr){\n",classname);
	    fprintf(file,"\t\tdouble _tolerances_[numEDO];\n");
	    fprintf(file,"\t\tdouble _aux_tol=0.0;\n");
	    fprintf(file,"\t\tdouble maxDt = this->dtime, minDt = this->dtime;\n");
	    fprintf(file,"\t\t//initializes the variables\n");
	    fprintf(file,"\t\tthis->previous_dt = this->dtime;\n");
	    fprintf(file,"\t\t//initializes the variables\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tthis->%s_new_ = this->%s_old_ = this->%s_ini_;\n", 
			cur->token.content, cur->token.content, cur->token.content);
			cur = cur->next;
	    }
	    fprintf(file,"\t\tthis->%s_new = this->%s;\n", difflist->diffheader->freevar.content ,
		    difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->timeSaving = this->%s;\n", difflist->diffheader->freevar.content);
	    
	    
	    
	    fprintf(file,"\t\tif(savingRate!=0)\n");
	    fprintf(file,"\t\t\tsave_step(fileptr, _ADAP_DT_);//save the initial conditions\n");
		    
	    fprintf(file,"\t\tdouble edos_old_aux_[numEDO];\n");
	    fprintf(file,"\t\tdouble edos_new_rk2_[numEDO],edos_new_euler_[numEDO];\n");
	    
	    fprintf(file,"\t\tdouble *_k1__ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *_k2__ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *_k_aux__;\n");
	    fprintf(file,"\t\tint iMaiorErro=0;\n");
	    fprintf(file,"\t\tint _cont_=0;\n");
	    fprintf(file,"\t\tdouble _soma_=0.0;\n");
		
	    fprintf(file,"\t\tthis->%s_new += this->d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content); 
	    fprintf(file,"\t\trightHandSideFunction.function(this);\n");
	    
	    fprintf(file,"\t\tfor(int i=0;i<numEDO;i++){\n");
	    fprintf(file,"\t\t\t_k1__[i] = this->getLadoDireito(i);\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tconst double __tiny_ = pow(this->abstol__, 2.0);\n");
	    
	    fprintf(file,"\t\twhile(this->%s_new<=finalTime){\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\tfor(int i=0; i<numEDO; i++){\n");
	    fprintf(file,"\t\t\t\t//stores the old variables in a vector\n");
	    fprintf(file,"\t\t\t\tedos_old_aux_[i] = this->getVariables(i);\n");
	    fprintf(file,"\t\t\t\t//computes euler method\n");
	    fprintf(file,"\t\t\t\tedos_new_euler_[i] = _k1__[i] * this->dtime + edos_old_aux_[i];\n");
	    fprintf(file,"\t\t\t\t//steps ahead to compute the rk2 method\n");
	    fprintf(file,"\t\t\t\tthis->setVariables(i, edos_new_euler_[i]);\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t\t//steps time ahead\n");
	    fprintf(file,"\t\t\tthis->time_new += this->dtime;\n");
	    fprintf(file,"\t\t\t//computes the right-hand side one step ahead\n");
	    fprintf(file,"\t\t\trightHandSideFunction.function(this);\n");
	    fprintf(file,"\t\t\t//restore the original old value\n");
	    fprintf(file,"\t\t\tthis->time_new -= this->dtime;//step back\n");
	    
	    fprintf(file,"\t\t\tdouble greatestError=0.0, auxError=0.0;\n");
	    fprintf(file,"\t\t\tfor(int i=0;i<numEDO;i++){\n");
	    fprintf(file,"\t\t\t\t//stores the new evaluation \n");
	    fprintf(file,"\t\t\t\t_k2__[i] = this->getLadoDireito(i);\n");
	    fprintf(file,"\t\t\t\t_aux_tol = fabs(edos_new_euler_[i])*reltol__;\n");
	    fprintf(file,"\t\t\t\t_tolerances_[i] = (abstol__ > _aux_tol )?abstol__:_aux_tol;\n");
	    fprintf(file,"\t\t\t\t//finds the greatest error between  the steps\n");
	    fprintf(file,"\t\t\t\tauxError = fabs(( (this->dtime/2)*(_k1__[i] - _k2__[i])) / _tolerances_[i]);\n");
	    fprintf(file,"\t\t\t\tif(auxError > greatestError){\n");
	    fprintf(file,"\t\t\t\t\tgreatestError = auxError;\n");
	    fprintf(file,"\t\t\t\t\tiMaiorErro = i;\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t}\n");
	    
	    fprintf(file,"\t\t\t///adapt the time step\n");
	    fprintf(file,"\t\t\tgreatestError+=__tiny_;\n");
	    fprintf(file,"\t\t\tint flag = this->getErrorCode(greatestError, 1.0);\n");
	    fprintf(file,"\t\t\tthis->previous_dt = this->dtime;\n");
	    fprintf(file,"\t\t\t//it doesn't accept the solution and cut h in a half\n");
	    fprintf(file,"\t\t\tif(flag==-1){\n");
	    fprintf(file,"\t\t\t\t//restore the old values to do it again\n");
	    fprintf(file,"\t\t\t\tfor(int i=0; i<numEDO; i++)\n");
	    fprintf(file,"\t\t\t\t\tthis->setVariables(i, edos_old_aux_[i]);\n");
	    fprintf(file,"\t\t\t\t//throw the results away and compute again\n");
	    fprintf(file,"\t\t\t\tthis->d%s = this->d%s/_DECREASE_DT_2_;//cut time step in a half\n", 
		    difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
			
	    fprintf(file,"\t\t\t}else{//it accepts the solutions\n");
	    fprintf(file,"\t\t\t\tif(savingRate!=0){\n");
	    fprintf(file,"\t\t\t\t\t//restore the previous value of old variables\n");
	    fprintf(file,"\t\t\t\t\tfor(int i=0; i<numEDO; i++)\n");
	    fprintf(file,"\t\t\t\t\t\tthis->setVariables(i, edos_old_aux_[i]);\n");
	    cur = rewind_list(difvarlist);
	    count =0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t\t\t\tthis->%s_new_ = edos_new_euler_[%d];\n", cur->token.content, 
			count);
		cur = cur->next;
		count++;
	    }
	    
	    
	    fprintf(file,"\t\t\t\t\tsave_step(fileptr, _ADAP_DT_);\n");
	    fprintf(file,"\t\t\t\t\tif(this->dtime>maxDt) maxDt = this->dtime;\n");
	    fprintf(file,"\t\t\t\t\tif(this->dtime<minDt) minDt = this->dtime;\n");
	    fprintf(file,"\t\t\t\t\t_soma_+=this->d%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t\t_cont_++;\n");
	    fprintf(file,"\t\t\t\t}\n");
	    
	    fprintf(file,"\t\t\t\tif(flag==3){\n");
	    fprintf(file,"\t\t\t\t\tthis->d%s = this->d%s*_DECREASE_DT_;\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t}else if(flag==4){\n");
	    fprintf(file,"\t\t\t\t\tthis->d%s = this->d%s*_INCREASE_DT_;\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t}else if(flag==0){\n");
	    fprintf(file,"\t\t\t\t\t//it just doesnt do anything\n");
	    fprintf(file,"\t\t\t\t}else{\n");
	    fprintf(file,"\t\t\t\t\tprintf(\"flag: %%d\\n\", flag);\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\tif(this->d%s > maxStep && maxStep!=0){\n",difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t\tthis->d%s = maxStep;\n",difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t}else if(this->d%s==0){\n",difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t\tprintf(\"Error: Time step is zero.\\n\");\n");
	    fprintf(file,"\t\t\t\t\treturn;\n");
	    fprintf(file,"\t\t\t\t}\n");
			
	    fprintf(file,"\t\t\t\t//change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração\n");
	    fprintf(file,"\t\t\t\t_k_aux__	= _k2__;\n");
	    fprintf(file,"\t\t\t\t_k2__		= _k1__;\n");
	    fprintf(file,"\t\t\t\t_k1__		= _k_aux__;\n");
	    
	    fprintf(file,"\t\t\t\t//sums the old dtime - the variable dtime is alreaady updated\n");
	    fprintf(file,"\t\t\t\tthis->time_new += this->previous_dt;\n");
	    
	    fprintf(file,"\t\t\t\t//it steps the method ahead, with euler solution\n");
	    fprintf(file,"\t\t\t\tfor(int i=0;i<numEDO;i++){\n");
	    fprintf(file,"\t\t\t\t\tthis->setVariables(i, edos_new_euler_[i]);\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t}\n");
		
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tif(savingRate!=0){\n");
	    fprintf(file,"\t\tprintf(\"Dt max: %%e dt min %%e, %%e %%d\\n\", maxDt, minDt, _soma_/_cont_, _cont_);\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t}\n");
	
	    /** *************** */
	    /// ///////********
	    
	    fprintf(file,"\tvoid %s::addt2(double finalTime, FILE *fileptr){\n",classname);
	    fprintf(file,"\t\tconst double _beta_safety_ = 0.8;\n");
	    fprintf(file,"\t\tdouble _tolerances_[numEDO];\n");
	    fprintf(file,"\t\tdouble _aux_tol=0.0;\n");
	    fprintf(file,"\t\tdouble maxDt = this->dtime, minDt = this->dtime;\n");
	    fprintf(file,"\t\t//initializes the variables\n");
	    fprintf(file,"\t\tthis->previous_dt = this->dtime;\n");
	    fprintf(file,"\t\t//initializes the variables\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tthis->%s_new_ = this->%s_old_ = this->%s_ini_;\n", 
			cur->token.content, cur->token.content, cur->token.content);
			cur = cur->next;
	    }
	    fprintf(file,"\t\tthis->%s_new = this->%s;\n", difflist->diffheader->freevar.content ,
		difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->timeSaving = this->%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tif(savingRate!=0)\n");
	    fprintf(file,"\t\t\tsave_step(fileptr, _ADAP_DT_);//save the initial conditions\n");
	    
	    fprintf(file,"\t\tdouble edos_old_aux_[numEDO];\n");
	    fprintf(file,"\t\tdouble edos_new_rk2_[numEDO],edos_new_euler_[numEDO];\n");
	    
	    fprintf(file,"\t\tdouble *_k1__ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *_k2__ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *_k_aux__;\n");
	    fprintf(file,"\t\tint iMaiorErro=0;\n");
	    fprintf(file,"\t\tint _cont_=0;\n");
	    fprintf(file,"\t\tdouble _soma_=0.0;\n");
		
	    fprintf(file,"\t\tthis->%s_new += this->d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\trightHandSideFunction.function(this);\n");
			
	    fprintf(file,"\t\tfor(int i=0;i<numEDO;i++){\n");
	    fprintf(file,"\t\t\t_k1__[i] = this->getLadoDireito(i);\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tconst double __tiny_ = pow(abstol__, 2.0);\n");
	    
	    fprintf(file,"\t\twhile(this->%s_new<=finalTime){\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\tfor(int i=0; i<numEDO; i++){\n");
	    fprintf(file,"\t\t\t\t//stores the old variables in a vector\n");
	    fprintf(file,"\t\t\t\tedos_old_aux_[i] = this->getVariables(i);\n");
	    fprintf(file,"\t\t\t\t//computes euler method\n");
	    fprintf(file,"\t\t\t\tedos_new_euler_[i] = _k1__[i] * this->dtime + edos_old_aux_[i];\n");
	    fprintf(file,"\t\t\t\t//steps ahead to compute the rk2 method\n");
	    fprintf(file,"\t\t\t\tthis->setVariables(i, edos_new_euler_[i]);\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t\t//steps time ahead\n");
	    fprintf(file,"\t\t\tthis->time_new += this->dtime;\n");
	    fprintf(file,"\t\t\t//computes the right-hand side one step ahead\n");
	    fprintf(file,"\t\t\trightHandSideFunction.function(this);\n");
	    fprintf(file,"\t\t\t//restore the original old value\n");
	    fprintf(file,"\t\t\tthis->time_new -= this->dtime;//step back\n");
	    
	    fprintf(file,"\t\t\tdouble greatestError=0.0, auxError=0.0;\n");
	    fprintf(file,"\t\t\tfor(int i=0;i<numEDO;i++){\n");
	    fprintf(file,"\t\t\t\t//stores the new evaluation \n");
	    fprintf(file,"\t\t\t\t_k2__[i] = this->getLadoDireito(i);\n");
	    fprintf(file,"\t\t\t\t_aux_tol = fabs(edos_new_euler_[i])*reltol__;\n");
	    fprintf(file,"\t\t\t\t_tolerances_[i] = (abstol__ > _aux_tol )?abstol__:_aux_tol;\n");
	    fprintf(file,"\t\t\t\t//finds the greatest error between  the steps\n");
	    fprintf(file,"\t\t\t\tauxError = fabs(( (this->dtime/2)*(_k1__[i] - _k2__[i])) / _tolerances_[i]);\n");
	    fprintf(file,"\t\t\t\tif(auxError > greatestError){\n");
	    fprintf(file,"\t\t\t\t\tgreatestError = auxError;\n");
	    fprintf(file,"\t\t\t\t\tiMaiorErro = i;\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t}\n");
	    
	    fprintf(file,"\t\t\t///adapt the time step\n");
	    fprintf(file,"\t\t\tgreatestError += __tiny_;\n");
	    fprintf(file,"\t\t\tthis->previous_dt = this->d%s;\n",difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\t\t\t///adapt the time step\n");
	    fprintf(file,"\t\t\tthis->d%s = _beta_safety_ * this->d%s * sqrt(1.0/greatestError);\n", 
		    difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\t\t\t//it doesn't accept the solution\n");
	    fprintf(file,"\t\t\tif(greatestError>=1.0){\n");		    
	    fprintf(file,"\t\t\t\t//restore the old values to do it again\n");
	    fprintf(file,"\t\t\t\tfor(int i=0; i<numEDO; i++)\n");
	    fprintf(file,"\t\t\t\t\tthis->setVariables(i, edos_old_aux_[i]);\n");
	    fprintf(file,"\t\t\t\t//throw the results away and compute again\n");		    
	    fprintf(file,"\t\t\t\t}else{//it accepts the solutions\n");
	    fprintf(file,"\t\t\t\t\tif(savingRate!=0.0){\n");
	    fprintf(file,"\t\t\t\t\t//restore the previous value of old variables\n");
	    fprintf(file,"\t\t\t\t\tfor(int i=0; i<numEDO; i++)\n");
	    fprintf(file,"\t\t\t\t\t\tthis->setVariables(i, edos_old_aux_[i]);\n");
	    cur = rewind_list(difvarlist);
	    count =0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t\t\tthis->%s_new_ = edos_new_euler_[%d];\n", cur->token.content, 
			count);
		cur = cur->next;
		count++;
	    }
		    
	    fprintf(file,"\t\t\t\t\tsave_step(fileptr, _ADAP_DT_);\n");
	    fprintf(file,"\t\t\t\t\tif(this->dtime>maxDt) maxDt = this->dtime;\n");
	    fprintf(file,"\t\t\t\t\tif(this->dtime<minDt) minDt = this->dtime;\n");
	    fprintf(file,"\t\t\t\t\t_soma_+=this->d%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t\t_cont_++;\n");
	    fprintf(file,"\t\t\t\t}\n");
		
	    fprintf(file,"\t\t\t\tif(this->dtime > maxStep && maxStep!=0){\n");
	    fprintf(file,"\t\t\t\t\tthis->dtime = maxStep;\n");
	    fprintf(file,"\t\t\t\t}else if(this->dtime==0){\n");
	    fprintf(file,"\t\t\t\t\tprintf(\"Error: Time step is zero.\\n\");\n");
	    fprintf(file,"\t\t\t\t\treturn;\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t//change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração\n");
	    fprintf(file,"\t\t\t\t_k_aux__	= _k2__;\n");
	    fprintf(file,"\t\t\t\t_k2__	= _k1__;\n");
	    fprintf(file,"\t\t\t\t_k1__	= _k_aux__;\n");
		
	    fprintf(file,"\t\t\t\t//sums the old dtime - the variable dtime is alreaady updated\n");
	    fprintf(file,"\t\t\t\tthis->time_new += this->previous_dt;\n");
	    fprintf(file,"\t\t\t\t//it steps the method ahead, with euler solution\n");
	    fprintf(file,"\t\t\t\tfor(int i=0;i<numEDO;i++){\n");
	    fprintf(file,"\t\t\t\t\tthis->setVariables(i, edos_new_euler_[i]);\n");
	    fprintf(file,"\t\t\t\t}\n");
		
	    fprintf(file,"\t\t\t}\n");
	   
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tif(savingRate!=0){\n");
	    fprintf(file,"\t\tprintf(\"Dt max: %%e dt min %%e, %%e  %%d\\n\", maxDt, minDt, _soma_/_cont_ , _cont_);\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t}\n");
	    /** *************** */
	
	    fprintf(file,"\tvoid %s::euler_OMP(double finalTime, FILE *fileptr, int nThreads){\n",classname);
	    fprintf(file,"\t\tomp_set_num_threads(nThreads);\n");
	    fprintf(file,"\t\tdouble *__NEW_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__OLD_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *temp;\n");
	    fprintf(file,"\t\tthis->%s_new = this->%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->timeSaving = this->%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->previous_dt = this->d%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t#pragma omp parallel firstprivate(__OLD_, __NEW_, temp)\n");
	    fprintf(file,"\t\t{\n");
	    //private parameter
	    cur = rewind_list(parvarlist);
	    fprintf(file,"\t\t\t//private parameters\n");
	    fprintf(file,"\t\t\tdouble ");
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_%s = %.10e, ",cur->token.content, cur->initialvalue);
		cur = cur->next;
	    }
	    
	    fprintf(file,"\n\t\t\t//private aux variables\n\t\t\t");
	    cur = rewind_list(algvarlist);
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_calc_%s=0, ",cur->token.content);
		//fprintf(file,"\tdouble %s_new;\n",cur->name);
		cur = cur->next;
	    }
	    fprintf(file,"\n\t\t\t//private right hand side variables\n\t\t\t");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_%s_lado_direito_, ", cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file,"\n\t\t\t//private time variables\n");
	    fprintf(file,"\t\t\t_prvt_%s_new = this->%s_new,\n",
		    difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t_prvt_d%s = this->d%s,\n", 
		    difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t_prvt_finalTime = finalTime, _prvt_savingRate = savingRate;\n");
		
	    cur = rewind_list(difvarlist);
	    count = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t__NEW_[%d] = __OLD_[%d] = %.10e;\n",
			count, count, cur->initialvalue);
		cur = cur->next;
		count++;
	    }
	    fprintf(file,"\t\t\tint *_prvt_tree_thread = tree_thread;\n");
	    fprintf(file,"\t\t\twhile(_prvt_%s_new<=_prvt_finalTime){\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t_prvt_%s_new += _prvt_d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		
		
		int countTreeItems = 0;
		fprintf(file,"\t\t\t\tif(omp_get_thread_num()==tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			//procurar por esta equação calc e depois pelo diff:  grafoDependencias[j].equationName);
			/// imprime os calc de uma arvore
			cur = rewind_list(algvarlist);
			numCalcs =0;
			while(cur != NULL)
			{
			    TokenNode *curlCalc = rewind_list(resolved_dep_list);
			    TokenNode *curCalc = NULL;
			    AlgList *curalgCalc = NULL;
			    int cont=0;
			    while(curlCalc != NULL)
			    {
				if(cont==numCalcs){
				    curalgCalc = rewind_list(alglist);
				    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
					curalgCalc = curalgCalc->next;
				    }
				    if(strcmp(grafoDependencias[j].equationName,curalgCalc->eq->token.content) ==0){
					fprintf(file,"\t\t\t\t\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
					curCalc = curalgCalc->eq;
					curCalc = curCalc->next->next;
					print_eq(file, curCalc, 1, 0);
					curlCalc = curlCalc->next;
					fprintf(file,";\n");
					break;
				    }
				    
				}else{
				    curCalc = curCalc->next->next;
				    curlCalc = curlCalc->next;
				}
				cont++;
			    }
			    
			    cur = cur->next;
			    numCalcs++;
			}
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t\t%s%s_lado_direito_= ",prefixo, curlDiff->diffheader->diffvar.content);
					print_eq(file, curDiff, 1,0);
					fprintf(file,";\n");    
					fprintf(file,"\t\t\t\t\t__NEW_[%d]= %s%s_lado_direito_ * %sd%s + __OLD_[%d];\n",cont, prefixo, curlDiff->diffheader->diffvar.content, prefixo,difflist->diffheader->freevar.content, cont);
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t\t}\n");
	    }
	    fprintf(file,"\t\t\t\t//synchronizing all threads\n");
	    fprintf(file,"\t\t\t\t#pragma omp barrier\n");
	    
	    fprintf(file,"\t\t\t\tif(_prvt_savingRate!=0){\n");
	    fprintf(file,"\t\t\t\t\t#pragma omp single\n");
	    fprintf(file,"\t\t\t\t\t{\n");
	    count = 0;
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t\t\t\tthis->%s_old_ = __OLD_[%d];\n",
			cur->token.content, count);
		fprintf(file,"\t\t\t\t\t\tthis->%s_new_ = __NEW_[%d];\n",
			cur->token.content,count);
		cur = cur->next;
		count++;
	    }
	    
	    fprintf(file,"\t\t\t\t\t\tthis->%s_new = _prvt_time_new;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t\t\tsave_step(fileptr, _EULER_);\n");
	    fprintf(file,"\t\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t}\n");
	    
	    fprintf(file,"\t\t\t\ttemp = __OLD_;\n");
	    fprintf(file,"\t\t\t\t__OLD_ = __NEW_;\n");
	    fprintf(file,"\t\t\t\t__NEW_= temp;\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t}\n");
	    /// ///////////////////
	    /// //////////////////////////
	    
	    
	    fprintf(file,"\tvoid %s::addt_OMP(double finalTime, FILE *fileptr, int nThreads){\n",classname);
	    fprintf(file,"\t\tomp_set_num_threads(nThreads);\n");
	    fprintf(file,"\t\tdouble *__NEW_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__OLD_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__OLD_AUX_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__TOL_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__ERROR_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__K1_  = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__K2_  = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__TEMP_;\n");
	    
	    
	    fprintf(file,"\t\tthis->%s_new = this->%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->timeSaving = this->%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->previous_dt = this->d%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t#pragma omp parallel firstprivate(__NEW_, __OLD_, __TEMP_, __TOL_, __K1_,__K2_, __ERROR_,__OLD_AUX_)\n");
	    fprintf(file,"\t\t{\n");
	    fprintf(file,"\t\t\tint *_prvt_tree_thread = tree_thread;\n");	    
	    fprintf(file,"\t\t\tdouble _prvt_rel_tol_=reltol__, _prvt_abs_tol_ =abstol__, _prvt_aux_tol;\n");	    
	    fprintf(file,"\t\t\tconst double __tiny_ = pow(_prvt_abs_tol_, 2.0);\n");
	    //private parameter
	    cur = rewind_list(parvarlist);
	    fprintf(file,"\t\t\t//private parameters\n");
	    fprintf(file,"\t\t\tdouble ");
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_%s = %.10e, ",cur->token.content, cur->initialvalue);
		cur = cur->next;
	    }
	    fprintf(file,"\n\t\t\t//private aux variables\n\t\t\t");
	    cur = rewind_list(algvarlist);
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_calc_%s=0.0, ",cur->token.content);
		//fprintf(file,"\tdouble %s_new;\n",cur->name);
		cur = cur->next;
	    }
	    fprintf(file,"\n\t\t\t//private right hand side variables\n\t\t\t");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_%s_lado_direito_, ",
			cur->token.content);
			cur = cur->next;
	    }
	    fprintf(file,"\n\t\t\t//private time variables\n");
	    fprintf(file,"\t\t\t_prvt_%s_new = this->%s_new,\n",
		difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t_prvt_d%s = this->d%s,\n", 
		difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
		fprintf(file,"\t\t\t_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;\n");
	    
	    cur = rewind_list(difvarlist);
	    count = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t__NEW_[%d] = __OLD_[%d] = %.10e;\n",
			count, count, cur->initialvalue);
			cur = cur->next;
			count++;
	    }
	    
	    fprintf(file,"\t\t\t_prvt_%s_new += _prvt_d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		int countTreeItems = 0;
		fprintf(file,"\t\t\tif(omp_get_thread_num()==_prvt_tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			//procurar por esta equação calc e depois pelo diff:  grafoDependencias[j].equationName);
			/// imprime os calc de uma arvore
			cur = rewind_list(algvarlist);
			numCalcs =0;
			while(cur != NULL)
			{
			    TokenNode *curlCalc = rewind_list(resolved_dep_list);
			    TokenNode *curCalc = NULL;
			    AlgList *curalgCalc = NULL;
			    int cont=0;
			    while(curlCalc != NULL)
			    {
				if(cont==numCalcs){
				    curalgCalc = rewind_list(alglist);
				    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
					curalgCalc = curalgCalc->next;
				    }
				    if(strcmp(grafoDependencias[j].equationName,curalgCalc->eq->token.content) ==0){
					fprintf(file,"\t\t\t\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
					curCalc = curalgCalc->eq;
					curCalc = curCalc->next->next;
					print_eq(file, curCalc,1,0);
					curlCalc = curlCalc->next;
					fprintf(file,";\n");
					break;
				    }
				}else{
				    curCalc = curCalc->next->next;
				    curlCalc = curlCalc->next;
				}
				cont++;
			    }
			    
			    cur = cur->next;
			    numCalcs++;
			}
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t__K1_[%d]= ",cont);
					print_eq(file, curDiff,1,0);
					fprintf(file,";\n");    
					fprintf(file,"\t\t\t\t__NEW_[%d]= __K1_[%d] * _prvt_dtime + __OLD_[%d];\n",cont,cont,cont);
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t}\n");
		
	    }//fim for das arvoresa
	    
	    fprintf(file,"\t\t\t//store the old iteration in a aux \n");
	    fprintf(file,"\t\t\t__TEMP_ = __OLD_;\n");
	    fprintf(file,"\t\t\t__OLD_ = __OLD_AUX_;\n");
	    fprintf(file,"\t\t\t__OLD_AUX_ = __TEMP_;\n");
	    
	    
	    fprintf(file,"\t\t\t//steps ahead with euler\n");
	    fprintf(file,"\t\t\t__TEMP_ = __NEW_;\n");
	    fprintf(file,"\t\t\t__NEW_ = __OLD_;\n");
	    fprintf(file,"\t\t\t__OLD_ = __TEMP_;\n");
	    
	    
	    fprintf(file,"\t\t\t//as threads devem começar o  laço ao mesmo tempo\n");
	    fprintf(file,"\t\t\t#pragma omp barrier\n");
	    
	    fprintf(file,"\t\t\twhile(_prvt_%s_new<=_prvt_finalTime){\n", difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\t\t\t\t_prvt_%s_new += _prvt_d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		int countTreeItems = 0;
		fprintf(file,"\t\t\t\tif(omp_get_thread_num()==_prvt_tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			//procurar por esta equação calc e depois pelo diff:  grafoDependencias[j].equationName);
			/// imprime os calc de uma arvore
			cur = rewind_list(algvarlist);
			numCalcs =0;
			while(cur != NULL)
			{
			    TokenNode *curlCalc = rewind_list(resolved_dep_list);
			    TokenNode *curCalc = NULL;
			    AlgList *curalgCalc = NULL;
			    int cont=0;
			    while(curlCalc != NULL)
			    {
				if(cont==numCalcs){
				    curalgCalc = rewind_list(alglist);
				    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
					curalgCalc = curalgCalc->next;
				    }
				    if(strcmp(grafoDependencias[j].equationName,curalgCalc->eq->token.content) ==0){
					fprintf(file,"\t\t\t\t\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
					curCalc = curalgCalc->eq;
					curCalc = curCalc->next->next;
					print_eq(file, curCalc,1,0);
					curlCalc = curlCalc->next;
					fprintf(file,";\n");
					break;
				    }
				}else{
				    curCalc = curCalc->next->next;
				    curlCalc = curlCalc->next;
				}
				cont++;
			    }
			    
			    cur = cur->next;
			    numCalcs++;
			}
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t\t__K2_[%d]= ",cont);
					print_eq(file, curDiff,1,0);
					fprintf(file,";\n");    
					fprintf(file,"\t\t\t\t\t_prvt_aux_tol = fabs(__OLD_[%d])*_prvt_rel_tol_;\n",cont);
					fprintf(file,"\t\t\t\t\t__TOL_[%d] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;\n",cont);
					fprintf(file,"\t\t\t\t\t__ERROR_[%d] = fabs((_prvt_dtime/2) * (__K1_[%d] - __K2_[%d])/__TOL_[%d]);\n",cont,cont,cont,cont); 
					
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t\t}\n");
		
	    }//fim for das arvoresa
	    
	    fprintf(file,"\t\t\t\t_prvt_%s_new -= _prvt_d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t#pragma omp barrier\n");
	    
	    fprintf(file,"\t\t\t\t///adapt the time step\n");
	    fprintf(file,"\t\t\t\tdouble __greatestError_=0.0;\n");
	    fprintf(file,"\t\t\t\tfor(int k=0;k<numEDO;k++){\n");
	    fprintf(file,"\t\t\t\t\t\tif(__ERROR_[k] > __greatestError_)\n");
	    fprintf(file,"\t\t\t\t\t\t\t__greatestError_ = __ERROR_[k];\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t__greatestError_ += __tiny_;\n");
	    fprintf(file,"\t\t\t\tint flag = this->getErrorCode(__greatestError_, 1.0);\n");
	    fprintf(file,"\t\t\t\t_prvt_previous_dt = _prvt_dtime;\n");
	    fprintf(file,"\t\t\t\t//it doesn't accept the solution and cut h in a half\n");
	    fprintf(file,"\t\t\t\tif(flag==-1){\n");
	    fprintf(file,"\t\t\t\t\t//throw the results away and compute again\n");
	    fprintf(file,"\t\t\t\t\t_prvt_dtime = _prvt_dtime/_DECREASE_DT_2_; //cut time step in a half\n");
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		int countTreeItems = 0;
		fprintf(file,"\t\t\t\t\tif(omp_get_thread_num()==_prvt_tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t\t\t__NEW_[%d] = __K1_[%d] * _prvt_dtime + __OLD_AUX_[%d];\n",cont,cont,cont); 
					
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t\t\t}\n");
		
	    }//fim for das arvoresa
	    
	    fprintf(file,"\t\t\t\t__TEMP_ = __NEW_;\n");
	    fprintf(file,"\t\t\t\t__NEW_ = __OLD_;\n");
	    fprintf(file,"\t\t\t\t__OLD_ = __TEMP_;\n");
	    
	    fprintf(file,"\t\t\t\t}else{//it accepts the solutions\n");
	    fprintf(file,"\t\t\t\t\tif(_prvt_savingRate!=0){\n");
	    fprintf(file,"\t\t\t\t\t\t#pragma omp single\n");
	    fprintf(file,"\t\t\t\t\t\t{\n");
	    count = 0;
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t\t\t\t\tthis->%s_old_ = __OLD_AUX_[%d];\n",
		    cur->token.content, count);
		fprintf(file,"\t\t\t\t\t\t\tthis->%s_new_ = __OLD_[%d];\n",
		    cur->token.content,count);
		cur = cur->next;
		count++;
	    }
	    fprintf(file,"\t\t\t\t\t\t\tthis->previous_dt = _prvt_previous_dt;\n");
	    fprintf(file,"\t\t\t\t\t\t\tthis->dtime = _prvt_dtime;\n");
	    
	    fprintf(file,"\t\t\t\t\t\t\tthis->%s_new = _prvt_time_new;\n",difflist->diffheader->freevar.content );
	    fprintf(file,"\t\t\t\t\t\t\tsave_step(fileptr, _ADAP_DT_);\n");
	    fprintf(file,"\t\t\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t\t}\n");

	    fprintf(file,"\t\t\t\t\tif(flag==3){\n");
	    fprintf(file,"\t\t\t\t\t\t_prvt_dtime = _prvt_dtime*_DECREASE_DT_;\n");
	    fprintf(file,"\t\t\t\t\t}else if(flag==4){\n");
	    fprintf(file,"\t\t\t\t\t\t_prvt_dtime = _prvt_dtime*_INCREASE_DT_;\n");
	    fprintf(file,"\t\t\t\t\t}else if(flag==0){\n");
	    fprintf(file,"\t\t\t\t\t\t//it just doesnt do anything\n");
	    fprintf(file,"\t\t\t\t\t}else{\n");
	    fprintf(file,"\t\t\t\t\t\tprintf(\"flag: %%d\\n\", flag);\n");
	    fprintf(file,"\t\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t\tif(_prvt_dtime > _prvt_maxStep && _prvt_maxStep!=0){\n");
	    fprintf(file,"\t\t\t\t\t\t_prvt_dtime = _prvt_maxStep;\n");
	    fprintf(file,"\t\t\t\t\t}else if(_prvt_dtime==0){\n");
	    fprintf(file,"\t\t\t\t\t\tprintf(\"Error: Time step is zero.\\n\");\n");
	    fprintf(file,"\t\t\t\t\t\tbreak;\n");
	    fprintf(file,"\t\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t\t//it steps the method ahead, with euler solution\n");
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		int countTreeItems = 0;
		fprintf(file,"\t\t\t\t\tif(omp_get_thread_num()==_prvt_tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t\t\t__NEW_[%d] = __K2_[%d] * _prvt_dtime + __OLD_[%d];\n",cont,cont,cont); 
					
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t\t\t}\n");
		
	    }//fim for das arvoresa
	    fprintf(file,"\t\t\t\t\t//store the old iteration in a aux \n");
	    fprintf(file,"\t\t\t\t\t__TEMP_ = __OLD_;\n");
	    fprintf(file,"\t\t\t\t\t__OLD_ = __OLD_AUX_;\n");
	    fprintf(file,"\t\t\t\t\t__OLD_AUX_ = __TEMP_;\n");
	    
	    fprintf(file,"\t\t\t\t\t//steps ahead with euler\n");
	    fprintf(file,"\t\t\t\t\t__TEMP_ = __NEW_;\n");
	    fprintf(file,"\t\t\t\t\t__NEW_ = __OLD_;\n");
	    fprintf(file,"\t\t\t\t\t__OLD_ = __TEMP_;\n");
	    
	    fprintf(file,"\t\t\t\t\t//change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração\n");
	    fprintf(file,"\t\t\t\t\t__TEMP_	= __K2_;\n");
	    fprintf(file,"\t\t\t\t\t__K2_	= __K1_;\n");
	    fprintf(file,"\t\t\t\t\t__K1_	= __TEMP_;\n");
	    
	    fprintf(file,"\t\t\t\t\t//sums the old dtime - the variable dtime is alreaady updated\n");
	    fprintf(file,"\t\t\t\t\t_prvt_time_new += _prvt_previous_dt;\n");
	    fprintf(file,"\t\t\t\t}//FIM ELSE\n");
	    fprintf(file,"\t\t\t\t#pragma omp barrier\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t}\n");

	    /// ///////////////////
	    
	    fprintf(file,"\tvoid %s::addt2_OMP(double finalTime, FILE *fileptr, int nThreads){\n",classname);
	    fprintf(file,"\t\tomp_set_num_threads(nThreads);\n");
	    fprintf(file,"\t\tdouble *__NEW_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__OLD_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__OLD_AUX_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__TOL_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__ERROR_ = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__K1_  = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__K2_  = (double*)malloc(sizeof(double)*numEDO);\n");
	    fprintf(file,"\t\tdouble *__TEMP_;\n");
	    
	    
	    fprintf(file,"\t\tthis->%s_new = this->%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->timeSaving = this->%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->previous_dt = this->d%s;\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t#pragma omp parallel firstprivate(__NEW_, __OLD_, __TEMP_, __TOL_, __K1_,__K2_, __ERROR_,__OLD_AUX_)\n");
	    fprintf(file,"\t\t{\n");
	    fprintf(file,"\t\t\tconst double _beta_safety_ = 0.8;\n");
	    fprintf(file,"\t\t\tint *_prvt_tree_thread = tree_thread;\n");	    
	    fprintf(file,"\t\t\tdouble _prvt_rel_tol_=reltol__, _prvt_abs_tol_ =abstol__, _prvt_aux_tol;\n");	    
	    fprintf(file,"\t\t\tconst double __tiny_ = pow(_prvt_abs_tol_, 2.0);\n");	    
	    //private parameter
	    cur = rewind_list(parvarlist);
	    fprintf(file,"\t\t\t//private parameters\n");
	    fprintf(file,"\t\t\tdouble ");
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_%s = %.10e, ",cur->token.content, cur->initialvalue);
		cur = cur->next;
	    }
	    fprintf(file,"\n\t\t\t//private aux variables\n\t\t\t");
	    cur = rewind_list(algvarlist);
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_calc_%s=0.0, ",cur->token.content);
		//fprintf(file,"\tdouble %s_new;\n",cur->name);
		cur = cur->next;
	    }
	    fprintf(file,"\n\t\t\t//private right hand side variables\n\t\t\t");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file," _prvt_%s_lado_direito_, ",
			cur->token.content);
			cur = cur->next;
	    }
	    fprintf(file,"\n\t\t\t//private time variables\n");
	    fprintf(file,"\t\t\t_prvt_%s_new = this->%s_new,\n",
		    difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t_prvt_d%s = this->d%s,\n", 
		    difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t_prvt_finalTime = finalTime, _prvt_savingRate = savingRate, _prvt_previous_dt = this->previous_dt, _prvt_maxStep = this->maxStep;\n");
	    
	    cur = rewind_list(difvarlist);
	    count = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t__NEW_[%d] = __OLD_[%d] = %.10e;\n",
			count, count, cur->initialvalue);
			cur = cur->next;
			count++;
	    }
	    
	    fprintf(file,"\t\t\t_prvt_%s_new += _prvt_d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		int countTreeItems = 0;
		fprintf(file,"\t\t\tif(omp_get_thread_num()==_prvt_tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			//procurar por esta equação calc e depois pelo diff:  grafoDependencias[j].equationName);
			/// imprime os calc de uma arvore
			cur = rewind_list(algvarlist);
			numCalcs =0;
			while(cur != NULL)
			{
			    TokenNode *curlCalc = rewind_list(resolved_dep_list);
			    TokenNode *curCalc = NULL;
			    AlgList *curalgCalc = NULL;
			    int cont=0;
			    while(curlCalc != NULL)
			    {
				if(cont==numCalcs){
				    curalgCalc = rewind_list(alglist);
				    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
					curalgCalc = curalgCalc->next;
				    }
				    if(strcmp(grafoDependencias[j].equationName,curalgCalc->eq->token.content) ==0){
					fprintf(file,"\t\t\t\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
					curCalc = curalgCalc->eq;
					curCalc = curCalc->next->next;
					print_eq(file, curCalc,1,0);
					curlCalc = curlCalc->next;
					fprintf(file,";\n");
					break;
				    }
				}else{
				    curCalc = curCalc->next->next;
				    curlCalc = curlCalc->next;
				}
				cont++;
			    }
			    
			    cur = cur->next;
			    numCalcs++;
			}
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t__K1_[%d]= ",cont);
					print_eq(file, curDiff,1,0);
					fprintf(file,";\n");    
					fprintf(file,"\t\t\t\t__NEW_[%d]= __K1_[%d] * _prvt_dtime + __OLD_[%d];\n",cont,cont,cont);
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t}\n");
		
	    }//fim for das arvoresa
	    
	    fprintf(file,"\t\t\t//store the old iteration in a aux \n");
	    fprintf(file,"\t\t\t__TEMP_ = __OLD_;\n");
	    fprintf(file,"\t\t\t__OLD_ = __OLD_AUX_;\n");
	    fprintf(file,"\t\t\t__OLD_AUX_ = __TEMP_;\n");
	    
	    
	    fprintf(file,"\t\t\t//steps ahead with euler\n");
	    fprintf(file,"\t\t\t__TEMP_ = __NEW_;\n");
	    fprintf(file,"\t\t\t__NEW_ = __OLD_;\n");
	    fprintf(file,"\t\t\t__OLD_ = __TEMP_;\n");
	    
	    
	    fprintf(file,"\t\t\t//as threads devem começar o  laço ao mesmo tempo\n");
	    fprintf(file,"\t\t\t#pragma omp barrier\n");
	    
	    fprintf(file,"\t\t\twhile(_prvt_%s_new<=_prvt_finalTime){\n", difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\t\t\t\t_prvt_%s_new += _prvt_d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		int countTreeItems = 0;
		fprintf(file,"\t\t\t\tif(omp_get_thread_num()==_prvt_tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			//procurar por esta equação calc e depois pelo diff:  grafoDependencias[j].equationName);
			/// imprime os calc de uma arvore
			cur = rewind_list(algvarlist);
			numCalcs =0;
			while(cur != NULL)
			{
			    TokenNode *curlCalc = rewind_list(resolved_dep_list);
			    TokenNode *curCalc = NULL;
			    AlgList *curalgCalc = NULL;
			    int cont=0;
			    while(curlCalc != NULL)
			    {
				if(cont==numCalcs){
				    curalgCalc = rewind_list(alglist);
				    while(strcmp(curalgCalc->eq->token.content, curlCalc->token.content)){
					curalgCalc = curalgCalc->next;
				    }
				    if(strcmp(grafoDependencias[j].equationName,curalgCalc->eq->token.content) ==0){
					fprintf(file,"\t\t\t\t\t%scalc_%s = ",prefixo, curalgCalc->eq->token.content);
					curCalc = curalgCalc->eq;
					curCalc = curCalc->next->next;
					print_eq(file, curCalc,1,0);
					curlCalc = curlCalc->next;
					fprintf(file,";\n");
					break;
				    }
				}else{
				    curCalc = curCalc->next->next;
				    curlCalc = curlCalc->next;
				}
				cont++;
			    }
			    
			    cur = cur->next;
			    numCalcs++;
			}
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t\t__K2_[%d]= ",cont);
					print_eq(file, curDiff,1,0);
					fprintf(file,";\n");    
					fprintf(file,"\t\t\t\t\t_prvt_aux_tol = fabs(__OLD_[%d])*_prvt_rel_tol_;\n",cont);
					fprintf(file,"\t\t\t\t\t__TOL_[%d] = (_prvt_abs_tol_ > _prvt_aux_tol )?_prvt_abs_tol_:_prvt_aux_tol;\n",cont);
					fprintf(file,"\t\t\t\t\t__ERROR_[%d] = fabs((_prvt_dtime/2) * (__K1_[%d] - __K2_[%d])/__TOL_[%d]);\n",cont,cont,cont,cont); 
					
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t\t}\n");
		
	    }//fim for das arvoresa
	    
	    fprintf(file,"\t\t\t\t_prvt_%s_new -= _prvt_d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\t#pragma omp barrier\n");
	    
	    fprintf(file,"\t\t\t\t///adapt the time step\n");
	    fprintf(file,"\t\t\t\tdouble __greatestError_=0.0;\n");
	    fprintf(file,"\t\t\t\tfor(int k=0;k<numEDO;k++){\n");
	    fprintf(file,"\t\t\t\t\t\tif(__ERROR_[k] > __greatestError_)\n");
	    fprintf(file,"\t\t\t\t\t\t\t__greatestError_ = __ERROR_[k];\n");
	    fprintf(file,"\t\t\t\t}\n");
	    
	    fprintf(file,"\t\t\t\t__greatestError_ += __tiny_;\n");
	    fprintf(file,"\t\t\t\t_prvt_previous_dt = _prvt_dtime;\n");
	    fprintf(file,"\t\t\t\t_prvt_dtime = _beta_safety_ * _prvt_dtime * sqrt(1.0/__greatestError_);\n");
	    fprintf(file,"\t\t\t\tif(__greatestError_>=1){\n");
	    fprintf(file,"\t\t\t\t\t//throw the results away and compute again\n");
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		int countTreeItems = 0;
		fprintf(file,"\t\t\t\t\tif(omp_get_thread_num()==_prvt_tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t\t\t__NEW_[%d] = __K1_[%d] * _prvt_dtime + __OLD_AUX_[%d];\n",cont,cont,cont); 
					
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t\t\t}\n");
		
	    }//fim for das arvoresa
	    
	    fprintf(file,"\t\t\t\t__TEMP_ = __NEW_;\n");
	    fprintf(file,"\t\t\t\t__NEW_ = __OLD_;\n");
	    fprintf(file,"\t\t\t\t__OLD_ = __TEMP_;\n");
	    
	    fprintf(file,"\t\t\t\t}else{//it accepts the solutions\n");
	    fprintf(file,"\t\t\t\t\tif(_prvt_savingRate!=0){\n");
	    fprintf(file,"\t\t\t\t\t\t#pragma omp single\n");
	    fprintf(file,"\t\t\t\t\t\t{\n");
	    count = 0;
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t\t\t\t\tthis->%s_old_ = __OLD_AUX_[%d];\n",
			cur->token.content, count);
			fprintf(file,"\t\t\t\t\t\t\tthis->%s_new_ = __OLD_[%d];\n",
				cur->token.content,count);
				cur = cur->next;
				count++;
	    }
	    fprintf(file,"\t\t\t\t\t\t\tthis->previous_dt = _prvt_previous_dt;\n");
	    fprintf(file,"\t\t\t\t\t\t\tthis->dtime = _prvt_dtime;\n");
	    
	    fprintf(file,"\t\t\t\t\t\t\tthis->%s_new = _prvt_time_new;\n",difflist->diffheader->freevar.content );
	    fprintf(file,"\t\t\t\t\t\t\tsave_step(fileptr, _ADAP_DT_);\n");
	    fprintf(file,"\t\t\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t\t}\n");
	    
	    
	    fprintf(file,"\t\t\t\t\tif(_prvt_dtime > _prvt_maxStep && _prvt_maxStep!=0){\n");
	    fprintf(file,"\t\t\t\t\t\t_prvt_dtime = _prvt_maxStep;\n");
	    fprintf(file,"\t\t\t\t\t}else if(_prvt_dtime==0){\n");
	    fprintf(file,"\t\t\t\t\t\tprintf(\"Error: Time step is zero.\\n\");\n");
	    fprintf(file,"\t\t\t\t\t\tbreak;\n");
	    fprintf(file,"\t\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t\t//it steps the method ahead, with euler solution\n");
	    prefixo = "_prvt_";
	    for(int i=1;i<=treeNumber;i++){
		int treeSize = counTreeItems(grafoDependencias, totalSize, i);//
		int countTreeItems = 0;
		fprintf(file,"\t\t\t\t\tif(omp_get_thread_num()==_prvt_tree_thread[%d])\n", i-1);
		fprintf(file,"\t\t\t\t\t{\n");
		
		for(int j=0; j<totalSize; j++){
		    if(grafoDependencias[j].treeId == i){
			
			///imprime as equações diferenciais
			cur = rewind_list(difvarlist);
			numEdos=0;
			while(cur != NULL)
			{
			    DiffList *curlDiff = rewind_list(difflist);
			    TokenNode *curDiff = NULL;
			    int cont=0;
			    while(curlDiff != NULL)
			    {
				if(cont==numEdos){
				    if(strcmp(grafoDependencias[j].equationName, curlDiff->diffheader->diffvar.content) ==0){
					curDiff = curlDiff->diffheader->eq->next;
					fprintf(file,"\t\t\t\t\t\t__NEW_[%d] = __K2_[%d] * _prvt_dtime + __OLD_[%d];\n",cont,cont,cont); 
					
					curlDiff = curlDiff->next;
					break;
				    }				
				}else{
				    curDiff = curlDiff->diffheader->eq->next;
				    curlDiff = curlDiff->next;
				}
				cont++;
				
			    }
			    cur = cur->next;
			    numEdos++;
			}
			
		    }
		}
		fprintf(file,"\t\t\t\t\t}\n");
		
	    }//fim for das arvoresa
	    fprintf(file,"\t\t\t\t\t//store the old iteration in a aux \n");
	    fprintf(file,"\t\t\t\t\t__TEMP_ = __OLD_;\n");
	    fprintf(file,"\t\t\t\t\t__OLD_ = __OLD_AUX_;\n");
	    fprintf(file,"\t\t\t\t\t__OLD_AUX_ = __TEMP_;\n");
	    
	    fprintf(file,"\t\t\t\t\t//steps ahead with euler\n");
	    fprintf(file,"\t\t\t\t\t__TEMP_ = __NEW_;\n");
	    fprintf(file,"\t\t\t\t\t__NEW_ = __OLD_;\n");
	    fprintf(file,"\t\t\t\t\t__OLD_ = __TEMP_;\n");
	    
	    fprintf(file,"\t\t\t\t\t//change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração\n");
	    fprintf(file,"\t\t\t\t\t__TEMP_	= __K2_;\n");
	    fprintf(file,"\t\t\t\t\t__K2_	= __K1_;\n");
	    fprintf(file,"\t\t\t\t\t__K1_	= __TEMP_;\n");
	    
	    fprintf(file,"\t\t\t\t\t//sums the old dtime - the variable dtime is alreaady updated\n");
	    fprintf(file,"\t\t\t\t\t_prvt_time_new += _prvt_previous_dt;\n");
	    fprintf(file,"\t\t\t\t}//FIM ELSE\n");
	    fprintf(file,"\t\t\t\t#pragma omp barrier\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t}\n");
			    
	    /// ///////////////////
	/// ///////////////////
	    fprintf(file,"\tint Solveode::getErrorCode(double error,double tolerance){\n");
	    fprintf(file,"\t\tif (error <0.5*tolerance){\n");
	    fprintf(file,"\t\t\t//Accept current solution, and increase the size of the next mechanics step.\n");
	    fprintf(file,"\t\t\treturn 4;\n");
	    fprintf(file,"\t\t}if (error>=0.5*tolerance && error<tolerance){\n");
	    fprintf(file,"\t\t\t//Accept current solution, and hold the step size fixed for the next mechanics step.\n");
	    fprintf(file,"\t\t\treturn 0;\n");
	    fprintf(file,"\t\t}if (error>=tolerance && error<2*tolerance){\n");
	    fprintf(file,"\t\t\t//Accept current solution, but decrease the size of the next mechanics step.\n");
	    fprintf(file,"\t\t\treturn 3;\n");
	    fprintf(file,"\t\t}else{\n");
	    fprintf(file,"\t\t\t//Throw current results away, cut the step size in half, and retry.\n");
	    fprintf(file,"\t\t\treturn -1;\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t}\n");
	/// ///////////////////
	
	    /******* FIM SET ***********************/
	    
	    /******* METODOS GET ********************/
	    fprintf(file,"\n\t//Get Methods\n");
	    
	    fprintf(file,"\n\tdouble %s::getVariables(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
	    
	    cur = rewind_list(difvarlist);
	    numVar = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tcase %d:\t\treturn %s_old_;    break;\n", numVar, cur->token.content);
		numVar++;
		cur = cur->next;
	    }
	    fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t}\n");
	    
	    
	    //////
	    fprintf(file,"\n\tdouble %s::getLadoDireito(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
				
	    cur = rewind_list(difvarlist);
	    numVar = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tcase %d:\t\treturn %s_lado_direito_;    break;\n", numVar, cur->token.content);
		numVar++;
		cur = cur->next;
	    }
	    fprintf(file,"\t\tdefault:\treturn 1;    break;\n\t\t}\n\t}\n");
	    
	    fprintf(file,"\n\tdouble %s::getParameters(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
	    cur = rewind_list(parvarlist);
	    numVar = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tcase %d:\t\treturn %s;    break;\n", numVar, cur->token.content);
		numVar++;
		cur = cur->next;
	    }
	    fprintf(file,"\t\tdefault:\tbreak;\n\t\t}\n\t}\n");
	    
	    fprintf(file,"\n\tdouble %s::getFreeVariable()\n\t{\n",classname);
	    
	    fprintf(file,"\t\treturn d%s;\n\t}\n", difflist->diffheader->freevar.content);
	    
	    
	    fprintf(file,"\n\t//Get Methods - Variables\n\n");
	    fprintf(file, "\tVariables %s::get_Variables()\n\t{\n\t\tVariables v(\"",classname);
	    
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"|%s#",cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file, "\");\n\t\treturn v;\n\t}\n");
	    
	    fprintf(file, "\tVariables %s::get_Parameters()\n\t{\n\t\tVariables v(\"",classname);
	    
	    cur = rewind_list(parvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"|%s#",cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file, "\");\n\t\treturn v;\n\t}\n");	
	    
	    fprintf(file, "\tVariables %s::get_FreeVariable()\n\t{\n\t\tVariables v(\"",classname);
	    
	    fprintf(file,"|%s#",difflist->diffheader->freevar.content); 		
	    fprintf(file, "\");\n\t\treturn v;\n\t}\n");	
	    
	    /********SET DO ARQUIVO **********************/
	    fprintf(file,"\n\tvoid %s::setParametersFromFile(char *filename)\n\t{\n",classname);
	    fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
	    fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setParametersFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
	    fprintf(file,"\tdouble value;\n\t\tint k = 0;\n\t\tVariables v = get_Parameters();\n\t\tint s = v.getQuantity();\n\t");
	    fprintf(file,"\tfor(;k<s;k++)\n\t\t{\n\t\t");
	    fprintf(file,"\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetParameters(k, value);\n\t\t}\n\t\tfclose(file);\n\t}\n");
	    
	    fprintf(file,"\n\tvoid %s::setVariablesFromFile(char *filename)\n\t{\n", classname);
	    fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
	    fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setVariablesFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
	    fprintf(file,"\tdouble value;\n\t\tint k = 0;\n\t\tVariables v = get_Variables();\n\t\tint s = v.getQuantity();\n\t");
	    fprintf(file,"\tfor(;k<s;k++)\n\t\t{\n\t\t");
	    fprintf(file,"\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetVariables(k, value);\n\t\t}\n\t\tfclose(file);\n\t}\n");
	    
	    fprintf(file,"\n\tvoid %s::setFreeVariableFromFile(char *filename)\n\t{\n", classname);
	    fprintf(file,"\t\tFILE *file;\n\t\tif((file = fopen(filename, \"r\")) == NULL)\n\t\t{\n");
	    fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - setFreeVariableFromFile - Unable to open file %%s\\n\", filename);\n\t\t\texit(1);\n\t\t}\n\t");
	    fprintf(file,"\tdouble value;\n\t\tfscanf(file,\"%%lf\", &value);\n\t\t\tsetFreeVariable(value);\n\t\tfclose(file);\n\t}\n");
	    /********FIM SET DO ARQUIVO ******************/
	    
	   
	    
	    
	    
	    
	    /******* METODO JACOBIAN ******************/
	    
	    fprintf(file,"\tdouble %s::jacobian(double final, double svRate)\n",classname);
	    fprintf(file,"\t{\n");		
	    
	    fprintf(file,"\t\t%s_new = %s;\n\n",difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tint num_results__		= (int)(final / svRate);\n");
	    fprintf(file,"\t\tint iterations	= (int)(final / this->dtime);\n");
	    
	    fprintf(file,"\t\tif(%s_vec__ != NULL)free( %s_vec__);\n\t\t\t%s_vec__ = (double *)malloc(sizeof(double)*num_results__);\n", difflist->diffheader->freevar.content,difflist->diffheader->freevar.content,difflist->diffheader->freevar.content);
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
		fprintf(file,"\t\tif(%s != NULL)free( %s);\n\t\t\t%s = (double *)malloc(sizeof(double)*num_results__);\n", cur->token.content,cur->token.content,cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file,"\t\tthis->timeSaving = d%s;\n\n", difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\t\tdouble diff=0;\n");
	    fprintf(file,"\t\tint counter=0;\n");
	    fprintf(file,"\t\tfor (int i = 0; i< iterations;i++ )\n");
	    fprintf(file,"\t\t{\n");
	    fprintf(file,"\t\t\tthis->%s_new += d%s;\n\n", difflist->diffheader->freevar.content, 
		    difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\trightHandSideFunction.function(this);\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\tthis->%s_new_ = this->%s_old_ + this->%s_lado_direito_ * this->d%s;\n", 
			cur->token.content, cur->token.content, cur->token.content, difflist->diffheader->freevar.content);
		cur = cur->next;
	    }
		
	    fprintf(file,"\t\t\tdiff =  _agos_round(this->%s_new - timeSaving, 10);\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\tif(diff==0){\n");
	    fprintf(file,"\t\t\t\tthis->timeSaving += svRate;\n");
	    fprintf(file,"\t\t\t\t%s_vec__[counter] = this->%s_new;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t\t%s[counter] = this->%s_new_;\n", 
			cur->token.content, cur->token.content);
		cur = cur->next;
	    }
	    
	   
	   fprintf(file,"\t\t\t\tcounter++;\n");
	   fprintf(file,"\t\t\t}\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\tthis->%s_old_ = this->%s_new_;\n", 
		    cur->token.content, cur->token.content);
		cur = cur->next;
	    }

	    fprintf(file,"\t\t}\n");
	    
	    fprintf(file,"\t\tdouble h_jac_[numEDO];\n");
	    fprintf(file,"\t\tdouble quociente = 1000.0;\n");
	    
	    cur = rewind_list(difvarlist);
	    int contaEDOS=0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\th_jac_[%d] = fabs(_agos_min(%s, num_results__) / _agos_max(%s, num_results__) );\n", contaEDOS, cur->token.content, cur->token.content);
		cur = cur->next;
		contaEDOS++;
	    }
	    fprintf(file,"\t\tfor(int l=0;l<numEDO;l++){\n");
	    fprintf(file,"\t\t\th_jac_[l] = (h_jac_[l]==0 || h_jac_[l]==AGOS_NAN || h_jac_[l]==AGOS_INF)?this->dtime:h_jac_[l];\n");
	    fprintf(file,"\t\t}\n");
	    
	    fprintf(file,"\t\tthis->timeSaving = this->d%s;\n\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tthis->%s_new = this->%s;\n\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t%s_old_ = %s_ini_;\n", cur->token.content ,cur->token.content);
		cur = cur->next;
	    }
	    
	    
	    fprintf(file,"\t\tdouble jacobian[numEDO][numEDO];\n");
	    fprintf(file,"\t\tdouble edos_new_aux_[numEDO];\n");
	    fprintf(file,"\t\tdouble edos_aux_[numEDO];\n");
	    fprintf(file,"\t\tFILE *filejac;\n");
	    fprintf(file,"\t\tfilejac = fopen(\"dat/jacobian.dat\", \"wb\");\n");
	    
	    
	    fprintf(file,"\t\tint counter2=0;\n");
	    fprintf(file,"\t\tfor (int i = 0; i< iterations;i++ ){\n");
	    fprintf(file,"\t\t\t%s_new += d%s;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\trightHandSideFunction.function(this);\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\tthis->%s_new_ = this->d%s*(this->%s_lado_direito_) + this->%s_old_;\n", 
			cur->token.content, difflist->diffheader->freevar.content,cur->token.content,cur->token.content);
			cur = cur->next;
	    }
	    
	    
	    fprintf(file,"\t\t\tdiff =  _agos_round(this->time_new - timeSaving, 10);\n");
	    fprintf(file,"\t\t\tif(diff==0){\n");
	    fprintf(file,"\t\t\t\tthis->timeSaving += svRate;\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t\t%s[counter2] = %s_new_;\n", cur->token.content,cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file, "\t\t\t\t%s_vec__[counter2] = %s_new;\n", difflist->diffheader->freevar.content, difflist->diffheader->freevar.content);
	    fprintf(file,"\t\t\t\tcounter2++;\n");
	    fprintf(file,"\t\t\t\t//salva os valores das variaveis diferenciaveis em auxiliares\n");
	    fprintf(file,"\t\t\t\tfor(int l=0;l<numEDO;l++){\n");
	    fprintf(file,"\t\t\t\t\tedos_aux_[l] = this->getLadoDireito(l);\n");
	    fprintf(file,"\t\t\t\t\tedos_new_aux_[l] = this->getVariables(l);\n");
	    fprintf(file,"\t\t\t\t}\n");
	    
	    fprintf(file,"\t\t\t\t//para cada coluna\n");
	    fprintf(file,"\t\t\t\tfor(int k=0;k<numEDO;k++){\n");
	    fprintf(file,"\t\t\t\t\t//escolhe uma variavel k, e calcula a derivada de todas as equacoes em relacao a k\n");
	    fprintf(file,"\t\t\t\t\tthis->setVariables(k, edos_new_aux_[k]+h_jac_[k]);\n");
	    fprintf(file,"\t\t\t\t\t//calcula tudo de novo com o novo valor de k\n");
	    fprintf(file,"\t\t\t\t\trightHandSideFunction.function(this);\n");
	    fprintf(file,"\t\t\t\t\tfor(int j=0;j<numEDO;j++){//para cada linha \n");
	    fprintf(file,"\t\t\t\t\t\tjacobian[j][k] = (this->getLadoDireito(j) - edos_aux_[j])/h_jac_[k];\n");
	    fprintf(file,"\t\t\t\t\t\tfprintf(filejac,\"%%f\\t\", jacobian[j][k]);\n");
	    fprintf(file,"\t\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t\tfprintf(filejac,\"\\n\");\n");
	    fprintf(file,"\t\t\t\t\t//agora tem que voltar para o que estava antes, sem somar dtime a variavel k\n");
	    fprintf(file,"\t\t\t\t\tfor(int l=0;l<numEDO;l++){\n");
	    fprintf(file,"\t\t\t\t\t\tthis->setVariables(l, edos_new_aux_[l]);\n");
	    fprintf(file,"\t\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\t\trightHandSideFunction.function(this);\n");
	    fprintf(file,"\t\t\t\t}\n");
	    fprintf(file,"\t\t\t\tfprintf(filejac,\"\\n\");\n");
	    fprintf(file,"\t\t\t}\n");
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,"\t\t\t%s_old_ = %s_new_;\n", cur->token.content,cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tfclose(filejac);\n");
	    fprintf(file,"\t\t\t//FILE *fileptr;\n");
	    fprintf(file,"\t\t\t//char filename[12];\n");
	    fprintf(file,"\t\t\t//sprintf(filename,\"%%dthread.dat\",numThreads);\n");
	    fprintf(file,"\t\t\t//fileptr = fopen(filename, \"wb\");\n");
	    fprintf(file,"\t\t\t//for(int i =0;i<num_results__;i++){\n");
	    fprintf(file,"\t\t\t//    fprintf(fileptr,\"%%f %%f\\n\", time_vec__[i], V[i]);\n");
	    fprintf(file,"\t\t\t//}\n");
	    fprintf(file,"\t\t\t//fclose(fileptr);\n");
	    fprintf(file, "\t\treturn 0;\n");
	    fprintf(file,"\t}\n");
	    /****** FIM JACOBIAN ******************/
	    prefixo="";
	    fprintf(file, "\n\tdouble* %s::getIndependentVar()\n\t{\n", classname);
	    fprintf(file, "\t\treturn %s_vec__;\n\t}\n", difflist->diffheader->freevar.content);
	    
	    fprintf(file,"\n\tdouble* %s::getSolution(int indVariable)\n\t{\n\t\tswitch(indVariable)\n\t\t{\n",classname);
	    /*** GET SOLUTION ********************/
	    
	    cur = rewind_list(difvarlist);
	    numVar = 0;
	    while(cur != NULL)
	    {
		fprintf(file,"\t\tcase %d:\t\treturn %s;    break;\n", numVar, cur->token.content);
		numVar++;		
		cur = cur->next;
	    }
	    
	    fprintf(file,"\t\tdefault:\treturn NULL;    break;\n\t\t}\n\t}\n");

	    //print_alg(file, classname, alglist);
// 				print_if(file, classname, iflist);
	    
	    free(grafoDependencias); 
	    
	    fprintf(file,"\tvoid %s::reInitCVODE()\n\t{\n\n",classname);
	    fprintf(file,"\t\tflag__ = CVodeReInit(cvode_mem_cvode__, f__, %s, dependent_variable__, CV_SV, reltol__, &abstol__);\n", difflist->diffheader->freevar.content);
	    fprintf(file,"\t\tif (check_flag(&flag__, \"CVodeReInit\", 1))\n\t\t\texit(1);\n");
	    fprintf(file, "\t}\n\n");
	    
	    fprintf(file, "\tvoid %s::setCVODEMaxStep(double maxstep)\n\t{\n\n",classname, file);
	    //Alterado por Ricardo 
	    fprintf(file, "\t\t\tCVodeSetMaxNumSteps(cvode_mem_cvode__, 1000000);\n\n");
	    
	    fprintf(file, "\t\t\trealtype initDT = this->getFreeVariable();\n");
	    fprintf(file, "\t\t\tCVodeSetInitStep(cvode_mem_cvode__, initDT);\n");
	    
	    fprintf(file, "\t\t\tflag__ = CVodeSetMaxStep(cvode_mem_cvode__, maxstep);\n");
	    fprintf(file, "\t\t\tif (check_flag(&flag__, \"CVodeSetMaxStep\", 1))\n");
	    fprintf(file, "\t\t\t\texit(1);\n");
	    
	    fprintf(file, "\t}\n\n"); 		
	    /******* FIM SOLUCAO EM DISCO ****************/
	    // 08/03/2008 print_alg(file, classname, alglist);
	    print_if(file, classname, iflist);
	    
	    
	    fprintf(file, "\n\n");
	    fprintf(file, "static int check_flag(void *flagvalue, char *funcname, int opt){\n");
	    fprintf(file, "\tint *errflag;\n");
	    fprintf(file, "\tif (opt == 0 && flagvalue == NULL) {\n");
	    fprintf(file, "\t\tfprintf(stderr, \"\\nSUNDIALS_ERROR: %%s() failed - returned NULL pointer\\n\\n\",funcname);\n");
	    fprintf(file, "\t\treturn(1);}\n");
	    fprintf(file, "\telse if (opt == 1) {\n");
	    fprintf(file, "\t\terrflag = (int *) flagvalue;\n");
	    fprintf(file, "\t\tif (*errflag < 0) {\n");
	    fprintf(file, "\t\t\tfprintf(stderr, \"\\nSUNDIALS_ERROR: %%s() failed with flag = %%d\\n\\n\",funcname, *errflag);\n");
	    fprintf(file, "\t\t\treturn(1); }}\n");
	    fprintf(file, "\telse if (opt == 2 && flagvalue == NULL) {\n");
	    fprintf(file, "\t\tfprintf(stderr, \"\\nMEMORY_ERROR: %%s() failed - returned NULL pointer\\n\\n\",funcname);\n");
	    fprintf(file, "\t\treturn(1); }\n");
	    fprintf(file, "\treturn 0;\n}\n\n");

	    /** SOLVECVODE *******************************************************/
	    fprintf(file,"\n\tvoid %s::solveCVODE(int firstcall__, int steps__, int num_results__, int method__, char *fileName,int cv_method__)\n\t{\n",classname);
	    fprintf(file,"\t\tFILE *f = fopen(fileName, \"wb\");\n");
	    fprintf(file,"\t\tif(!f){\n");
	    fprintf(file,"\t\t\tfprintf(stderr,\"ERROR - solveToFile - Unable to open file %%s\\n\",fileName);\n");
	    fprintf(file,"\t\t\texit(1);\n");
	    fprintf(file,"\t\t}\n");
	    fprintf(file,"\t\tdouble dtL, dtM, dtMax=0.0,  dtMin=0.0 ;\n");
	    
	    fprintf(file,"\t\tdtM = 0.0;\n");
	    fprintf(file,"\t\tint iout = 1;\n");
	    fprintf(file,"\t\t%s_new = %s;\n", rewind_list(parvarlist)->token.content, rewind_list(parvarlist)->token.content);
	    fprintf(file,"\t\trealtype tout = %s_new+d%s;\n\n",rewind_list(parvarlist)->token.content,rewind_list(parvarlist)->token.content );
	    fprintf(file,"\t\trealtype cvodeDT;\n");
	    
	    fprintf(file,"\t\tint counter_it__ = 0;\n");
	    fprintf(file,"\t\tint offset_step = steps__ / num_results__;\n");
	    fprintf(file,"\t\tint aux = steps__%%num_results__;\n");
	    fprintf(file,"\t\tint num_iterations_bak = steps__ + 1;\n");
	    
	    fprintf(file,"\t\tif(firstcall__==1){\n");
	    
	    fprintf(file,"\t\t\tif(steps__ <= 0)\n\t\t\t\tsteps__ = 1;\n");
	    fprintf(file,"\t\t\tdependent_variable__ = N_VNew_Serial(%d);\n",numEdos);
	    fprintf(file,"\t\t\tif(check_flag((void *)dependent_variable__, \"N_VNew_Serial\", 0))\n");
	    fprintf(file,"\t\t\texit(1);\n");
	    fprintf(file,"\t\t\tdepvar__ = (double*)malloc(sizeof(double)*%d);\n",numEdos);
	    fprintf(file,"\t\t\tif(depvar__ == NULL){\n");
	    fprintf(file,"\t\t\tfprintf(stderr, \"ERROR Cannot allocate memory for depvar__\\n\");\n");
	    fprintf(file,"\t\t\texit(0);\n");
	    fprintf(file,"\t\t\t}\n");
	    
	    int i=0;
	    cur = rewind_list(difvarlist);
	    while(cur != NULL){
		fprintf(file,"\t\t\tNV_Ith_S(dependent_variable__, %d) = %s_ini_;\n", i, cur->token.content);
		cur = cur->next;
		i++;
	    }
	    fprintf(file,"\t\t\tit_countx = 0;\n");
	    fprintf(file,"\t\t\tint nonlineariteration__;\n");
	    fprintf(file,"\t\t\tif(method__ == CV_BDF) nonlineariteration__ = CV_NEWTON;\n");
	    fprintf(file,"\t\t\telse nonlineariteration__ = CV_FUNCTIONAL;\n");
	    fprintf(file,"\t\t\tcvode_mem_cvode__ = CVodeCreate(method__, nonlineariteration__);\n");
	    fprintf(file,"\t\t\tif (check_flag((void *)cvode_mem_cvode__, \"CVodeCreate\", 0))\n");
	    fprintf(file,"\t\t\texit(1);\n");
	    fprintf(file,"\t\t\tflag__ = CVodeMalloc(cvode_mem_cvode__, f__, time, dependent_variable__, CV_SS, reltol__, &abstol__);\n");
	    fprintf(file,"\t\t\tif (check_flag(&flag__, \"CVodeMalloc\", 1))\n");
	    fprintf(file,"\t\t\texit(1);\n");
	    
	    fprintf(file,"\t\t\tthis->setCVODEMaxStep(this->maxStep);\n");
// 	    fprintf(file,"\t\t\tthis->reInitCVODE();\n");
	    fprintf(file,"\t\t\tswitch(cv_method__){\n");
	    fprintf(file,"\t\t\tcase 1 :\n");
	    fprintf(file,"\t\t\tflag__ = CVDense(cvode_mem_cvode__, %d);\n", numEdos);
	    fprintf(file,"\t\t\tif (check_flag(&flag__, \"CVDense\", 1))	exit(1);\n");
	    fprintf(file,"\t\t\tbreak;\n");
	    fprintf(file,"\t\t\tcase 2:\n");
	    fprintf(file,"\t\t\tflag__ = CVDiag(cvode_mem_cvode__);\n");
	    fprintf(file,"\t\t\tif (check_flag(&flag__, \"CVDiag\", 1))	exit(1);\n");
	    fprintf(file,"\t\t\tbreak;\n");
	    fprintf(file,"\t\t\tcase 3:\n");
	    fprintf(file,"\t\t\tflag__ = CVBand(cvode_mem_cvode__, %d, NULL, NULL);\n", numEdos);
	    fprintf(file,"\t\t\tif (check_flag(&flag__, \"CVBand\", 1))	exit(1);\n");
	    fprintf(file,"\t\t\tbreak;\n");
	    fprintf(file,"\t\t\t}\n");
	    fprintf(file,"\t\t\tCVodeSetFdata(cvode_mem_cvode__, (void*)this);\n");
				
	    cur = rewind_list(difvarlist);
	    while(cur != NULL){
		fprintf(file,"\t\t\t%s_old_ = %s_ini_;\n", cur->token.content,cur->token.content);
		cur = cur->next;
	    }
	    
	    fprintf(file,"\t\t}\n");
	    
	    fprintf(file,"\t\twhile(1){\n");
	    fprintf(file,"\t\t\tflag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &%s_new, CV_NORMAL);\n",rewind_list(parvarlist)->token.content);
	    
	    cur = rewind_list(difvarlist);
	    i =0;
	    while(cur != NULL){
		fprintf(file, "\t\t\t%s_old_ = NV_Ith_S(dependent_variable__, %d);\n", cur->token.content, i);
		cur = cur->next;
		i++;
	    }
	    
	    fprintf(file,"\t\t\tif (check_flag(&flag__, \"CVode\", 1)){printf(\"pau check_flag na %%d iteracao\\n\", iout); break;}\n");
	    
	    fprintf(file,"\t\t\tif (flag__ == CV_SUCCESS){\n");
	    fprintf(file,"\t\t\t\tif((iout != aux) && ((iout - aux)%%offset_step == 0)) {\n");
// 	    fprintf(file,"\t\t\t\t\tfprintf(f,\"%%f %%f\\n\", time_new, V_old_);\n");
	    fprintf(file,"\t\t\t\tCVodeGetLastStep(cvode_mem_cvode__, &cvodeDT);\n");
	
	    /**fprintf(file,"\t\t\t\tfprintf(f,\"%%.8e");
	    int list_size = get_list_size(difvarlist);
	    for(int i = 0; i < list_size; i++)
		fprintf(file," %%.8e");
	    fprintf(file,"\\n\",%s_new", difflist->diffheader->freevar.content);
	    
	    cur = rewind_list(difvarlist);
	    while(cur != NULL)
	    {
		fprintf(file,",%s_old_", cur->token.content);
		cur = cur->next;
	    }
	    fprintf(file,");\n");*/
	    fprintf(file,"\t\t\t\t\tfprintf(f,\"%%.8e %%.8e %%.8e\\n\", time_new, V_old_,cvodeDT );\n");
	    
	    fprintf(file,"\t\t\t\t\tcounter_it__++;\n");
	    fprintf(file,"\t\t\t\t}\n");

	    fprintf(file,"\t\t\t\tiout++;\n");
	    fprintf(file,"\t\t\t\ttout += dtime; // timestep\n");
	    fprintf(file,"\t\t\t\tflag__ = CVodeGetLastStep(cvode_mem_cvode__, &dtL);\n");
	    fprintf(file,"\t\t\t\tdtM += dtL;\n");
	    fprintf(file,"\t\t\t\tif(dtL>dtMax) dtMax = dtL;\n");
	    fprintf(file,"\t\t\t\tif(dtL<dtMin) dtMin = dtL;\n");
	    
	    fprintf(file,"\t\t\t}\n");
	
	    fprintf(file,"\t\t\tif (iout == num_iterations_bak ){ break;}\n");
	    
	    fprintf(file, "\t\t}\n");
	    fprintf(file, "\t\tfclose(f);\n");
	    fprintf(file, "\t\tlong int nsteps;\n");
	    fprintf(file, "\t\tint flag = CVodeGetNumSteps(cvode_mem_cvode__ , &nsteps);\n");
	    fprintf(file, "\t\tprintf(\"max: %%e min: %%e media: %%e %%d\\n\", dtMax, dtMin, dtM/nsteps, nsteps);\n");
	    fprintf(file, "\t}\n");
	    
	    /** FIM SOLVECVODE ***************************************************/ 				
	    
	    fprintf(file, "static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__){\n");
	    fprintf(file, "\t%s *ode = (%s *) f_data__;\n", classname, classname);
	    int listsize = get_list_size(difvarlist);
	    fprintf(file, "\tfor(int i = 0; i<%d; i++)\n", listsize);
	    fprintf(file, "\t\tode->setVariables( i ,NV_Ith_S(dependent_variable__, i));\n");
	    fprintf(file, "\tode->setParameters(0,time);\n");
	    fprintf(file, "\trightHandSideFunction.function(ode);\n");
	    fprintf(file, "\tfor(int i = 0; i<%d; i++)\n",  listsize);
	    fprintf(file, "\t\tNV_Ith_S(dep_var_dot__, i) = ode->getLadoDireito(i);\n");
	    fprintf(file, "\treturn 0;\n}\n\n");
	    
	    fprintf(file,"\n\nfloat __agos_factorial(int f){\n\tif(f>=0 & f<2)\n\t\treturn 1.0;\n\telse if(f < 0)\n\t\treturn 0.0/0.0;\n\tfor(int i=f-1; i>=2; i--)\n\t\tf *= i;\n\treturn (float)f;\n}\n");
	    
	    fprintf(file,"double _agos_max(double* vector, int size){\n");
	    fprintf(file,"\tdouble max =vector[0];\n");
	    fprintf(file,"\tint i;\n");
	    fprintf(file,"\tfor(i=1;i<size;i++){\n");
	    fprintf(file,"\t\tif(vector[i]>max) max = vector[i];\n");
	    fprintf(file,"\t}\n");
	    fprintf(file,"\treturn max;\n}\n");
	    
	    fprintf(file,"double _agos_min(double* vector, int size){\n");
	    fprintf(file,"\tdouble min = vector[0];\n");
	    fprintf(file,"\tint i;\n");
	    
	    fprintf(file,"\tfor(i=1;i<size;i++){\n");
	    fprintf(file,"\t\tif(vector[i]<min) min = vector[i];\n");
	    fprintf(file,"\t}\n");
	    fprintf(file,"\treturn min;\n}\n");

	    fprintf(file,"double _agos_round( double x, int places )\n");
	    fprintf(file,"{\n");
	    fprintf(file,"\tdouble const shift = powf( 10.0f, places );\n");
	    fprintf(file,"\tx *= shift;\n");
	    fprintf(file,"\tx = floorf( x + 0.5f );\n");
	    fprintf(file,"\tx /= shift;\n");
	    fprintf(file,"\treturn x;\n");
	    fprintf(file,"}\n");
	    
	    
	    fclose(file);
	    return 0;
	}
	

#endif _COMPILE_H_