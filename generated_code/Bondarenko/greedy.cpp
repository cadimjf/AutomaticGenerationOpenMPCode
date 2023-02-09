#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <malloc.h>
typedef struct str_elemento{
    int indice;
    int value;
    str_elemento *next;
}typ_elemento;

typedef struct str_raiz{
    typ_elemento *root;
    int totalSum;
}typ_raiz;


bool quickSort(typ_elemento **trees, int elements);
double index_min(typ_raiz* vector, int size);
int greedyPartitionProblem(int *vector,int elements, int subSets, int *threads);


//  quickSort
//
//  This public-domain C implementation by Darel Rex Finley.
//
//  * Returns YES if sort was successful, or NO if the nested
//    pivots went too deep, in which case your array will have
//    been re-ordered, but probably not sorted correctly.
//
//  * This function assumes it is called with valid parameters.
//
//  * Example calls:
//    quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
//    quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7
bool quickSort(typ_elemento **trees, int elements) {
    
#define  MAX_LEVELS  1000

int  piv, beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R ;

typ_elemento *piv_i = (typ_elemento*)malloc(sizeof(typ_elemento));

    beg[0]=0; end[0]=elements;
    while (i>=0) 
    {
	L=beg[i]; R=end[i]-1;
	if (L<R) {
	    piv_i = trees[L];
	    piv=trees[L]->value;
	    
	    if (i==MAX_LEVELS-1) return 0;
	    while (L<R) {
		while (trees[R]->value>=piv && L<R) R--;
		if (L<R){
		    int aux = L;
		    trees[aux]=trees[R];		    
		    L++;
		}
		while (trees[L]->value<=piv && L<R) L++; 
		if (L<R){
		    int aux = R;
		    trees[aux]=trees[L]; 
		    R--;
		}
	    }
	    trees[L]=piv_i;
	    beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L; 
	    
	}
	    else {
		i--; 
	    }
    }
    return 1; 
}

// finds the mininum
double index_min(typ_raiz* vector, int size){
    int min = vector[0].totalSum;
    int i_min = 0;
    int i;
    for(i=1;i<size;i++){
	if(vector[i].totalSum<min){
	    min = vector[i].totalSum;
	    i_min = i;
	}
    }
    return i_min;
}

//receives a vector and make a partition in n sub sets
//vector -> vector with the complete set;
// elements -> number of elements in vector;
//subsets -> number of subsets
//threads -> output vector containing the jobs distributed among threads. 
//		its size is the number of elements, and it stores the number of the thread that is 
//		going to solve the job i

int greedyPartitionProblem(int *vector,int elements, int subSets, int *threads_tree)
{
    int teste = 4;
//     printf("%d\n", --teste);printf("%d", teste);getchar();getchar();
    typ_elemento **trees = (typ_elemento**)malloc(sizeof(typ_elemento)*elements);
    for(int i=0; i<elements;i++){
	trees[i] = (typ_elemento*)malloc(sizeof(typ_elemento));
	trees[i]->indice	= i;
	trees[i]->value 	= vector[i];
	trees[i]->next 	= NULL;
    }
    
    
//     
//     printf("\n");
//     for(int i=0; i<elements; i++) printf( "[%d]->%d ", trees[i]->indice, trees[i]->value); printf("\n");
    
    //sorts the vector
    quickSort(trees, elements);
//     for(int i=0; i<elements; i++) printf( "[%d]->%d ", trees[i]->indice, trees[i]->value); printf("\n");    
    
    //creates the subset struct 
    typ_raiz *subSets_list;
    subSets_list = (typ_raiz*)malloc(sizeof(typ_raiz)*subSets);
    //init
    for(int i=0;i<subSets; i++){
	subSets_list[i].root = NULL;
	subSets_list[i].totalSum = 0;
    }
//     printf("%d %d:\n", elements, subSets);
    
    //começa do maior
    for(int i=elements-1; i>=0; i--)
    {
// 	printf("%d de %d:\n", i, elements);
	//finds the subset with minimal total sum
	int  min_index = index_min(subSets_list, subSets);
	//creates a new node struct
	typ_elemento *novoElemento = (typ_elemento*)malloc(sizeof(typ_elemento));
	novoElemento->indice	= trees[i]->indice;
	novoElemento->value 	= trees[i]->value;
	novoElemento->next 	= NULL;
// 	printf("min %d acrescenta %d\n",  min_index, i);
	//inserts in the proper place
	subSets_list[min_index].totalSum += trees[i]->value;
	typ_elemento *element = subSets_list[min_index].root ;
	typ_elemento *temp;
	
	if(element == NULL){    
	    subSets_list[min_index].root = novoElemento;
	}else
	{
	    while(element->next != NULL){
		element = element->next;
	    }
	    element->next = novoElemento;
	}
	
    }
    //prints
    for(int i=0;i<subSets; i++){
// 	printf("Total: %d:\t", subSets_list[i].totalSum);
	typ_elemento *element = subSets_list[i].root;
	while(element != NULL){
// 	    printf("[%d]->%d\t", element->indice, element->value);
	    threads_tree[element->indice] = i;
	    element = element->next;
	}
// 	printf( "\n");    
    }
    
}


/*int main(int argc, char** argv) {
    int vector[16] = {14, 37, 60, 9, 120, 12, 16,   9,    2,   17,    1,    1,    3,    1,    1,    10};
    //sorts the vector
   
    
    for(int i=0; i<16; i++)
    {
	printf( "%d\t", vector[i]);
    }
    printf( "\n");
       
    greedyPartitionProblem(vector,16, 10);
    
}*/

