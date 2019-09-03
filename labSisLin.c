#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{


  // inicializa gerador de nr aleatoreos
  srand(20192);
  
  //lendo sistema 
  SistLinear_t *SL = lerSistLinear ();

  // gerando sistema
 // SistLinear_t *SL = alocaSistLinear(10);
	//inicializaSistLinear (SL, 5, 12); // o diagonalmente dominante não está correto

  printf("----------------- SISTEMA LINENAR -----------------\n\n");

	prnSistLinear (SL);

  // alocando vetor solução
	real_t *x = (real_t*)malloc(SL->n * sizeof(double));

  printf("\n\n----------------- ELIMINAÇÃO GAUSS -----------------\n\n");
  for (int i = 0; i < SL->n; i++) x[i] = 0.0;
  double tempo = timestamp();
  eliminacaoGauss(SL, x, 0);
  printf("\nTempo: %lf\n", timestamp() - tempo);
  printf("Norma L2: %f\n", normaL2Residuo(SL, x));
  prnVetor(x, SL->n);

  printf("\n\n----------------- GAUSS - JACOBI -----------------\n\n");  
	for (int i = 0; i < SL->n; i++) x[i] = 0.0;
  tempo = timestamp();
	gaussJacobi(SL, x, EPS);
  printf("\nTempo: %lf\n", timestamp() - tempo);
  printf("Norma L2: %f\n", normaL2Residuo(SL, x));
  prnVetor(x, SL->n);
	

  printf("\n\n----------------- GAUSS - SEIDEL -----------------\n\n");
	for (int i = 0; i < SL->n; i++) x[i] = 0.0;
  tempo = timestamp();
	gaussSeidel(SL, x, EPS);
  printf("\nTempo: %lf\n", timestamp() - tempo);
  printf("Norma L2: %f\n", normaL2Residuo(SL, x));
  prnVetor(x, SL->n);
}

