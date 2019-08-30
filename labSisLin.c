#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
    // inicializa gerador de nr aleatoreos
    srand(20192);
    
  SistLinear_t *SL = alocaSistLinear(5);
  inicializaSistLinear (SL, 0, 12);
  //prnSistLinear (SL);

  real_t *x = malloc(SL->n * sizeof(double));

  eliminacaoGauss(SL, x, 1);
}

