#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{


    // inicializa gerador de nr aleatoreos
    srand(20192);
    
    SistLinear_t *SL = lerSistLinear ();
    /*alocaSistLinear(5);
  	inicializaSistLinear (SL, 5, 12);*/
  	prnSistLinear (SL);

  	real_t *x = (real_t*)malloc(SL->n * sizeof(double));

 	eliminacaoGauss(SL, x, 1);
  	printf("ELIMINACAO GAUSS\n");
  	prnVetor(x, SL->n);

	for (int i = 0; i < SL->n; i++) {
  	x[i] = 0;
 	}

	gaussJacobi(SL, x, EPS);
	printf("GAUSS JACOBI\n");
	prnVetor(x, SL->n);
	
	for (int i = 0; i < SL->n; i++) {
  	x[i] = 0;
 	}

	gaussSeidel(SL, x, EPS);
	printf("GAUSS SEIDEL\n");
	prnVetor(x, SL->n);
 // printf("agora tem que ser o mesmo de antes!\n");
  //prnSistLinear (SL);
}

