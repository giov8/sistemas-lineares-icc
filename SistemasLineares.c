#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "SistemasLineares.h"

#define POS(n, i, j) ((n*i)+j)
#define N SL->n

SistLinear_t *copiaSL(SistLinear_t *SL_origem){
  SistLinear_t *SL_dest = alocaSistLinear(SL_origem->n);

  memcpy(SL_dest->A, SL_origem->A, SL_origem->n*SL_origem->n*sizeof(real_t));
  memcpy(SL_dest->b, SL_origem->b, SL_origem->n*sizeof(real_t));

  return SL_dest;
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
{

}


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param pivotamento flag para indicar se o pivotamento parcial deve ser feito (!=0)

  \return código de erro. 0 em caso de sucesso.
*/


/*int posicao(int n, int i, int j) {
  return n * i + j;
}*/

/*Função que encontra o MÁXIMO DA COLUNA e devolve a coluna que tem o maior valor*/
int encontraMax(SistLinear_t *SL, int i) {
  int LinhaMaior = i; 

  for (int k = i+1; k < SL->n; ++k) {
    if (fabs(SL->A[POS(SL->n,k,i)]) > fabs(SL->A[POS(SL->n,LinhaMaior,i)]))
      LinhaMaior = k;
  }
  return LinhaMaior;
}

void trocaLinha(SistLinear_t *SL, int i, int iPivo) {
  // troca linha i pela linha iPivo

  real_t *aux = (real_t*)malloc (SL->n*sizeof(real_t));

  memcpy(aux, &SL->A[POS(SL->n, i,0)], SL->n*sizeof(real_t)); // Copia a linha i para o vetor auxiliar
  memcpy(&SL->A[POS(SL->n, i, 0)], &SL->A[POS(SL->n,iPivo, 0)], SL->n*sizeof(real_t)); // Copia a linha iPivo para a linha i
  memcpy(&SL->A[POS(SL->n, iPivo, 0)], aux, SL->n*sizeof(real_t)); // Copia o vetor auxiliar para a linha iPivo

  real_t aux2 = SL->b[i];
  SL->b[i] = SL->b[iPivo];
  SL->b[iPivo] = aux2;

  free(aux);
}

int eliminacaoGauss (SistLinear_t *SL, real_t *x, int pivotamento)
{

 SistLinear_t *SLcp = copiaSL(SL);

 for (int i = 0; i < SLcp->n; ++i) {
    if (pivotamento) {
      int iPivo = encontraMax(SLcp, i);
      if ( i != iPivo )
        trocaLinha(SLcp, i, iPivo);
    }
    for(int k=i+1; k < SL->n; ++k) {
       double m = SLcp->A[POS(SLcp->n,k,i)] / SLcp->A[POS(SLcp->n,i,i)];
       SLcp->A[POS(SLcp->n,k,i)] = 0.0;
       for(int j=i+1; j < SLcp->n; ++j)
          SLcp->A[POS(SLcp->n,k,j)] -= SLcp->A[POS(SLcp->n,i,j)] * m;
       SLcp->b[k] -= SL->b[i] * m;
    } 
 }
 //RETROSUBSTITUIÇÃO
  real_t soma;
  int i = N - 1;
  x[i] = SL->b[i] / SL->A[POS(N,i,i)];
  for (i = N+2; i >= 0; i--) {
    soma = SL->b[i];
    for (int j = i+1; j < N; j++) {
      soma -= SL->A[POS(N,i,j)];
    }
    x[i] = soma/SL->A[POS(N,i,i)];
  }

   liberaSistLinear(SLcp);
  return 0;
}


/*!
  \brief Método de Gauss-Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, real_t erro)
{
  real_t *xk = (real_t*)malloc(N*sizeof(real_t));
  real_t soma = 0;
  real_t norma = 0.0;
  int i, k = 0;

  for (k = 0; k < MAXIT; ++k) {
    for (int i = 0; i < N; ++i) {

      for (int j = 0; j < N; j++) if (j != i) soma = soma + SL->A[POS(N, i, j)] * x[j];

      xk[i] = (1/SL->A[POS(N, i, i)]) * (SL->b[i] - soma);
    }
  
    for (i = 0; i < N; i++) {
      real_t diff = fabs(x[i]-xk[i]);
      if (norma < diff) norma = diff;
    }

    if (norma < erro) {
      memcpy(x,xk,N*sizeof(real_t));
    }
  }

  memcpy(x,xk,N*sizeof(real_t));
  free(xk);
  if (k == MAXIT)
    return -1; // não houve convergencia!
  return k; // retorna nmro de iterações 
}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro)
{
  real_t norma = 1.0 + erro; // valor só para entrar no FOR
  real_t *xk = (real_t*)malloc(N*sizeof(real_t));
  int k, i, j;

  for (k = 0; norma > erro; k++) {
    for (i = 0; i < N; i++) {
      real_t soma = 0.0;
      for (j = 0; j < i; j++) soma += SL->A[POS(N,i,j)] * xk[j];
      for (j = i+1; j < N; j++) soma += SL->A[POS(N, i, j)] * x[j];
      xk[i] = (SL->b[i] - soma)/SL->A[POS(N,i,i)];
    }
    // calculando norma
    norma = 0.0;
    for (i = 0; i < N; i++) {
      real_t diff = fabs(x[i] - xk[i]);
      if (norma < diff) norma = diff;
    }
    memcpy(x,xk,N*sizeof(real_t));
  }
}


// Alocaçao de memória
SistLinear_t* alocaSistLinear (unsigned int tam)
{
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  if ( SL ) {
    SL->A = (real_t *) malloc(tam * tam * sizeof(real_t));
    SL->b = (real_t *) malloc(tam * sizeof(real_t));

    if (!(SL->A) || !(SL->b))
      liberaSistLinear(SL);
  }
  
  SL->n = tam;

  return (SL);
}

// Liberacao de memória
void liberaSistLinear (SistLinear_t *SL)
{
  free(SL->A);
  free(SL->b);
  free(SL);
}

/*!
  \brief Cria coeficientes e termos independentes do SL
  *
  \param SL Ponteiro para o sistema linear
  \param tipo Tipo de sistema linear a ser criado. Pode ser: comSolucao,
  eqNula, eqProporcional, eqCombLinear, hilbert 
  \param coef_max Maior valor para coeficientes e termos independentes
*/
void inicializaSistLinear (SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max)
{
  unsigned int tam = SL->n;
  // para gerar valores no intervalo [0,coef_max[
  real_t invRandMax = ((real_t)coef_max / (real_t)RAND_MAX);

  // inicializa vetor b
  for (unsigned int i=0; i<tam; ++i) {
    SL->b[i] = (real_t)rand() * invRandMax;
  }
    
  if (tipo == hilbert) {
    for (unsigned int i=0; i<tam; ++i) {
      for (unsigned int j=0; j<tam; ++j)  {
	      SL->A[i*tam+j] = 1.0 / (real_t)(i+j+1);
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<tam; ++i) {
      for (unsigned int j=0; j<tam; ++j)  {
	      SL->A[i*tam+j] = (real_t)rand() * invRandMax;
      }
    }
    if (tipo == eqNula) {
      // sorteia eq a ser "nula"
      unsigned int nula = rand() % tam;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[nula*tam+j] = 0.0;
      }
      SL->b[nula] = 0.0;
    } 
    else if (tipo == eqProporcional) {
      // sorteia eq a ser "proporcional" e valor
      unsigned int propDst = rand() % tam;
      unsigned int propSrc = (propDst + 1) % tam;
      real_t mult = (real_t)rand() * invRandMax;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[propDst*tam+j] = SL->A[propSrc*tam+j] * mult;
      }
      SL->b[propDst] = SL->b[propSrc] * mult;
    } 
    else if (tipo == eqCombLinear) {
      // sorteia eq a ser "combLinear"
      unsigned int combDst = rand() % tam;
      unsigned int combSrc1 = (combDst + 1) % tam;
      unsigned int combSrc2 = (combDst + 2) % tam;
      for (unsigned int j=0; j<tam; ++j) {
	      SL->A[combDst*tam+j] = SL->A[combSrc1*tam+j] + SL->A[combSrc2*tam+j];
      }
      SL->b[combDst] = SL->b[combSrc1] + SL->b[combSrc2];
    }
    else if (tipo == diagDominante) {
      // aumenta o expoente dos termos da diagonal principal
      for (unsigned int i=0; i<tam; ++i) {
        SL->A[i*tam+i] *= (real_t)tam;
      }
    }

  }
}


SistLinear_t *lerSistLinear ()
{
  unsigned int n;
  SistLinear_t *SL;
  
  scanf("%d",&n);

  SL = alocaSistLinear (n);
  
  for(int i=0; i < n; ++i)
    for(int j=0; j < n; ++j)
      scanf ("%g", &SL->A[i*n+j]);

  for(int i=0; i < n; ++i)
    scanf ("%g", &SL->b[i]);
  
  return SL;
}


void prnSistLinear (SistLinear_t *SL)
{
  int n=SL->n;

  for(int i=0; i < n; ++i) {
    printf("\n");
    for(int j=0; j < n; ++j)
      printf ("%10g", SL->A[i*n+j]);
    printf ("   |   %g", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor (real_t *v, unsigned int n)
{
  int i;

  printf ("\n");
  for(i=0; i < n; ++i)
      printf ("%10g ", v[i]);
  printf ("\n\n");

}

