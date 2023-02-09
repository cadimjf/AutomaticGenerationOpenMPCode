
#ifndef _LEX_H_
#define _LEX_H_

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <string.h>
#include "fsm_util.h"



// TOKEN type definitions /////////////////////////////////////////////////////

#define	END		0	// end of file
#define	I_OP	1	// infix operator
#define	P_OP	2	// prefix operator
#define R_OP    12	// relational operator < > = <= >=
#define L_OP	16	// logical operator and or
#define PIOP	3	// mix operator only for '-' it may be unary
#define	NUMB	4	// number
#define	VARI	5	// variable
#define	DIFF	6	// diff operator
#define	BVAR	7	// bvar specification
#define OPAR	8	// open parenthesis (apply)
#define CPAR	9	// close parenthesis (#)
#define MATH	11	// math token
#define PI_W	13	// piecewise token
#define PIEC	14	// piece token
#define	OT_W	15	// otherwise token
#define ROOT	16	// root operator
#define DGRE	17	// used to especify root degree

#define RDF		18
#define MODL	19
#define COMP	20
#define VARB	21
#define UNIT	22



#define STACK_SIZE 50

// GLOBAL VARS ////////////////////////////////////////////////////////////////
xmlNode *pointer;
xmlNode *prev_pointer;

int		lex_index = -1;
int		lex_index_bound = STACK_SIZE;
xmlNode **lex_stack = NULL;



// FUNCTIONS PROTOTYPES ///////////////////////////////////////////////////////
void		init_lex(xmlNode *node);
Token		lex();
xmlNode* 	next(xmlNode *ptr);
Token		get_token(xmlNode *ptr);
void		lex_push(xmlNode *ptr);
xmlNode*	lex_pop();
int 		lex_empty();
char*		strNoSpace(char* str);

// FUNCTIONS IMPLEMENTATION ///////////////////////////////////////////////////

/* init_lex *******************************************************************
 * Initializes the stack and the xmlNode pointer with root node of the        *
 * document to be parsed                                                      *
 * ***************************************************************************/
void init_lex(xmlNode *node)
{
	if(lex_stack != NULL)
		free(lex_stack);
	lex_stack = (xmlNode **)calloc(STACK_SIZE,sizeof(xmlNode*));
	pointer = node;
}


/* lex ************************************************************************
 * Returns the next token in the DOM tree                                     *
 * this function sets the token values                                        *
 *****************************************************************************/
 
Token lex()
{
	Token token;
	token.tag = NULL;
	pointer = next(pointer);
	// if the pointer is NULL try to pop the stack
	if(pointer == NULL)
	{
		if(!lex_empty())
		{
			pointer = lex_pop()->next;
			token.tag = "#";
			token.type = CPAR;
			token.content = ")";
			}else{ token.tag = "END"; token.type = END;}
	}		
	else
	{	
		prev_pointer = pointer;
		// ignoring text nodes	
		while(pointer->type == XML_TEXT_NODE)
		{
			pointer = pointer->next;
			if(pointer == NULL)
				break;
		}
		if(pointer != NULL) // we may have a problem with the else that don't exists
		{
			token = token = get_token(pointer);//.tag = (char*)pointer->name;
			if(pointer->children != NULL)
			{ 
				// only push to stack these tags
				if(!strcmp((char*)pointer->name, "apply") ||
				!strcmp((char*)pointer->name, "math") ||
				!strcmp((char*)pointer->name, "bvar") ||
				!strcmp((char*)pointer->name, "piecewise") ||
				!strcmp((char*)pointer->name, "piece") ||
				!strcmp((char*)pointer->name, "degree") ||
				!strcmp((char*)pointer->name, "model") ||
				!strcmp((char*)pointer->name, "component") ||
				!strcmp((char*)pointer->name, "otherwise"))
				{
					lex_push(pointer);
					pointer = pointer->children;
				}else pointer = pointer->next;
			}
			else pointer = pointer->next;
		}
	}		
	
	return token;
}


/* get_token ******************************************************************
 * Translate node into token information                                      *
 * ***************************************************************************/
 
Token get_token(xmlNode *ptr)
{
	Token token;
	memset(&token,0,sizeof(Token));
	
	if(!strcmp((char*)ptr->name,"RDF"))
	{
		token.tag = (char*)ptr->name;
		token.content = token.tag;
		token.type = RDF;
		
	}else if(!strcmp((char*)ptr->name,"model"))
	{
		token.tag = (char*)ptr->name;
		token.content = token.tag;
		token.type = MODL;
		
		
	}else if(!strcmp((char*)ptr->name,"units"))
	{
		token.tag = (char*)ptr->name;
		token.content = token.tag;
		token.type = UNIT;
		
		
	}else if(!strcmp((char*)ptr->name,"component"))
	{
		token.tag = (char*)ptr->name;
		token.content = token.tag;
		token.type = COMP;
		
	}else if(!strcmp((char*)ptr->name,"variable"))
	{
		token.tag = (char*)ptr->name;
		token.content = token.tag;
		token.type = VARB;
		
	}else if(!strcmp((char*)ptr->name,"diff"))
	{
		token.tag = (char*)ptr->name;
		token.type = DIFF;
		
	}else if(!strcmp((char*)ptr->name,"apply"))
	{
		token.content = "(";
		token.tag = (char*)ptr->name;
		token.type = OPAR;
		
	}else if(!strcmp((char*)ptr->name,"plus"))
	{
		token.content = "+";
		token.tag = (char*)ptr->name;
		token.type = I_OP;
	}else if(!strcmp((char*)ptr->name,"divide"))
	{
		token.content = "/";
		token.tag = (char*)ptr->name;
		token.type = I_OP;
	}else if(!strcmp((char*)ptr->name,"times"))
	{
		token.content = "*";
		token.tag = (char*)ptr->name;
		token.type = I_OP;
	}else if(!strcmp((char*)ptr->name,"minus"))
	{
		token.content = "-";
		token.tag = (char*)ptr->name;
		token.type = PIOP;
	}else if(!strcmp((char*)ptr->name,"quotient"))
	{
		token.content = "/";
		token.tag = (char*)ptr->name;
		token.type = I_OP;
	}else if(!strcmp((char*)ptr->name,"rem"))
	{
		token.content = "%";
		token.tag = (char*)ptr->name;
		token.type = I_OP;
	}else if(!strcmp((char*)ptr->name,"math"))
	{
		token.tag = (char*)ptr->name;
		token.type = MATH;
	}else if(!strcmp((char*)ptr->name,"power"))
	{
		token.content = "pow(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"exp"))
	{
		token.content = "exp(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"floor"))
	{
		token.content = "floor(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"xor"))
	{
		token.content = "__agos_xor(";
		token.tag = (char*)ptr->name;
		token.type = L_OP;
	}else if(!strcmp((char*)ptr->name,"ln"))
	{
		token.content = "log(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"log"))
	{
		token.content = "log10(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"sin"))
	{
		token.content = "sin(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"cos"))
	{
		token.content = "cos(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"tan"))
	{
		token.content = "tan(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"csc"))
	{
		token.content = "1/sin(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"sec"))
	{
		token.content = "1/cos(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"cot"))
	{
		token.content = "1/tan(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"sinh"))
	{
		token.content = "sinh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"cosh"))
	{
		token.content = "cosh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"tanh"))
	{
		token.content = "tanh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"csch"))
	{
		token.content = "1/sinh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"sech"))
	{
		token.content = "1/cosh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"coth"))
	{
		token.content = "1/tanh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arcsin"))
	{
		token.content = "asin(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arccos"))
	{
		token.content = "acos(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arctan"))
	{
		token.content = "atan(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arccsc"))
	{
		token.content = "1/asin(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arcsec"))
	{
		token.content = "1/acos(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arccot"))
	{
		token.content = "1/atan(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arcsinh"))
	{
		token.content = "asinh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arccosh"))
	{
		token.content = "acosh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arctanh"))
	{
		token.content = "atanh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arccsch"))
	{
		token.content = "1/asinh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arcsech"))
	{
		token.content = "1/acosh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"arccoth"))
	{
		token.content = "1/atanh(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"ceiling"))
	{
		token.content = "ceil(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"abs"))
	{
		token.content = "fabs(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}if(!strcmp((char*)ptr->name,"eq"))
	{
		token.tag = (char*)ptr->name;
		token.content = "=";
		token.type = R_OP;
	}else if(!strcmp((char*)ptr->name,"ci"))
	{
		token.tag = (char*)ptr->name;
		token.content = strNoSpace((char*)xmlNodeGetContent(ptr));
		token.type = VARI;
		//printf("%s#",token.content);
		// including a var in the var list
		varlist = add_list(token, varlist);
	}
	else if(!strcmp((char*)ptr->name,"cn"))
	{
		token.tag = (char*)ptr->name;
		//token.content = strNoSpace((char*)xmlNodeGetContent(ptr));
		//Alterado por Ricardo
		// mais algumas modificacoes por Ronan atof substituido por strtod e escrevendo %.15e ao inves de %lf
		double teste= strtod(strNoSpace((char*)xmlNodeGetContent(ptr)),NULL);
		char strteste[255];
// 		sprintf(strteste,"%.15e", teste);
		sprintf(strteste,"%.10e", teste);
		token.content=  strNoSpace((char*)strteste);
		//fim-trecho
		token.type = NUMB;
	}
	// CONSTANTS //////////////////////////////////////////////////////////
	else if(!strcmp((char*)ptr->name,"pi"))
	{
		token.tag = (char*)ptr->name;
		token.content = "3.141592653589793116";
		token.type = NUMB;
	}
	else if(!strcmp((char*)ptr->name,"exponentiale"))
	{
		token.tag = (char*)ptr->name;
		token.content = "2.718281828459045091"; 
		token.type = NUMB;
	}
	else if(!strcmp((char*)ptr->name,"true"))
	{
		token.tag = (char*)ptr->name;
		token.content = "1"; 
		token.type = NUMB;
	}
	else if(!strcmp((char*)ptr->name,"false"))
	{
		token.tag = (char*)ptr->name;
		token.content = "0"; 
		token.type = NUMB;
	}
	else if(!strcmp((char*)ptr->name,"infinity"))
	{
		token.tag = (char*)ptr->name;
		token.content = "AGOS_INF"; 
		token.type = NUMB;
	}
	else if(!strcmp((char*)ptr->name,"notanumber"))
	{
		token.tag = (char*)ptr->name;
		token.content = "AGOS_NAN"; 
		token.type = NUMB;
	}
	///////////////////////////////////////////////////////////////////////
	else if(!strcmp((char*)ptr->name,"bvar"))
	{
		token.tag = (char*)ptr->name;
		token.content = NULL;
		token.type = BVAR;
	}else if(!strcmp((char*)ptr->name,"piecewise"))
	{
		token.tag = (char*)ptr->name;
		token.content = NULL;
		token.type = PI_W;
	}else if(!strcmp((char*)ptr->name,"piece"))
	{
		token.tag = (char*)ptr->name;
		token.content = NULL;
		token.type = PIEC;
	}else if(!strcmp((char*)ptr->name,"otherwise"))
	{
		token.tag = (char*)ptr->name;
		token.content = NULL;
		token.type = OT_W;
	}
	else if(!strcmp((char*)ptr->name,"gt"))
	{
		token.tag = (char*)ptr->name;
		token.content = ">";
		token.type = R_OP;
	}
	else if(!strcmp((char*)ptr->name,"lt"))
	{
		token.tag = (char*)ptr->name;
		token.content = "<";
		token.type = R_OP;
	}else if(!strcmp((char*)ptr->name,"geq"))
	{
		token.content = ">=";
		token.tag = (char*)ptr->name;
		token.type = R_OP;
	}else if(!strcmp((char*)ptr->name,"leq"))
	{
		token.content = "<=";
		token.tag = (char*)ptr->name;
		token.type = R_OP;
	}else if(!strcmp((char*)ptr->name,"neq"))
	{
		token.content = "!=";
		token.tag = (char*)ptr->name;
		token.type = R_OP;
	}else if(!strcmp((char*)ptr->name,"and"))
	{
		token.tag = (char*)ptr->name;
		token.content = "&&";
		token.type = L_OP;
	}else if(!strcmp((char*)ptr->name,"or"))
	{
		token.tag = (char*)ptr->name;
		token.content = "||";
		token.type = L_OP;
	}else if(!strcmp((char*)ptr->name,"root"))
	{
		token.content = "pow(";
		token.tag = (char*)ptr->name;
		token.type = ROOT;
	}else if(!strcmp((char*)ptr->name,"not"))
	{
		token.content = "!("; 
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}else if(!strcmp((char*)ptr->name,"degree"))
	{
		token.content = "";
		token.tag = (char*)ptr->name;
		token.type = DGRE;
	}else if(!strcmp((char*)ptr->name,"factorial"))
	{
		token.content = "__agos_factorial(";
		token.tag = (char*)ptr->name;
		token.type = P_OP;
	}
	
	return token;
}
/** Arithmetic_algebra *****************************************************
*   quotient, divide, plus, minus, times, rem. and, or, not                *
***************************************************************************/
/*
else if(!xmlStrcmp(name,(const xmlChar *)"not"))
{
	mmlnode->name = "!";     mmlnode->type = prefix;    
}
*/
/** Relations **************************************************************
*   eq, neq, gt, lt, geq, leq                                              *
***************************************************************************/    

/** Functions **************************************************************
*  exp, ln, log, sin, cos, tan, sec, csc, cot, sinh, cosh, tanh, sech,     *
*  csch, coth, arcsin, arccos, arctan, arcsec, arccsc, arccot, arcsinh,    *
*  arccosh, arctanh, arcsech, arccsch, arccoth,                            *
***************************************************************************/

void lex_push(xmlNode *ptr)
{	
	if(lex_index >= lex_index_bound)
	{
		lex_index_bound += STACK_SIZE/2;
		lex_stack = (xmlNode **)realloc(lex_stack, sizeof(xmlNode*)*lex_index_bound);
		if(lex_stack == NULL)
		{
			printf("ERROR - lex stack reallocating memory failed\n");
			exit(1);
		}		
	}
	lex_index++;
	lex_stack[lex_index] = ptr;	
}

xmlNode* lex_pop() 
{
	xmlNode *ptr = NULL;
	if(lex_index >= 0)
	{
		ptr = lex_stack[lex_index];
		lex_index--;
	}
	return ptr;
}

int lex_empty()
{
	if(lex_index >= 0)
		return 0;
	else return 1;
}

xmlNode* next(xmlNode *ptr)
{
	xmlNode *cur=ptr;
	if(cur != NULL)
	while(!strcmp((char*)cur->name,"text")|| !strcmp((char*)cur->name,"comment"))
	{
		cur = cur->next;
		if(cur == NULL)
			break;
	}
	return cur;		
}

/* strNoSpace *****************************************************************
 * Returns a new string containing no spaces                                  *
 *****************************************************************************/  

char* strNoSpace(char* str)
{
	int size = strlen(str);
	int space = 0;
	int i;
	for(i=0;i<size;i++)
		if(str[i] == ' ')
			space++;
	
	char *strn = (char*)calloc(size-space +1, sizeof(char));

	space=0;
	for(i=0;i<size;i++)
		if(str[i] != ' ')
		{
			strn[space] = str[i];
			space++;
		}
	return strn;
}

#endif
