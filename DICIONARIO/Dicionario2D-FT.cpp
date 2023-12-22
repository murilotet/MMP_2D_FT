

#include <stdio.h>
#include <stdlib.h>
#include "Dicionario2D-FT.h"

/********************************************************************************************/
/*					Metodos da classe Bloco				    */
/********************************************************************************************/

Bloco :: Bloco(void) {
	pixel = NULL;
	linhas = colunas = 0;
}
void Bloco :: Dimensiona(int lin, int col) {
//dimensiona o bloco como uma matriz 2D
int n;
/*	if(pixel != NULL) free(pixel);
	if(lin*col == 0) {
		linhas = colunas = 0;
		return;
	}
	if((pixel = (int **)malloc(lin*sizeof(int *))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	for(n = 0; n < lin; n++) {
		if((pixel[n] = (int *)malloc(col*sizeof(int))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}
	}*/
	if(lin*col == 0) {
		for(n = 0; n < linhas; n++) delete [] pixel[n];
		delete [] pixel;
	}
	else {
		pixel = new int* [lin];
		for(n = 0; n < lin; n++) pixel[n] = new int [col];
	}
	linhas = lin;
	colunas = col;
}

float Bloco :: Distancia(Bloco A) {
//Calcula a distancia ate o bloco A (erro quadratico)
int n, m;
float D=0;
float Dp;
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			Dp = pixel[n][m] - A.pixel[n][m];
			D = D + Dp*Dp;
		}
	}
	return(D);
}

float Bloco :: DistanciaF(Bloco A, float Da) {
//Calcula a distancia ate o bloco A (erro quadratico). Contudo, se D ultrapassar Da, aborta o
// calculo e retorna 2*Da
int n, m;
float D=0;
float Dp;
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			Dp = pixel[n][m] - A.pixel[n][m];
			D = D + Dp*Dp;
			if(D > Da) return(2*Da);
		}
	}
	return(D);
}

void Bloco :: Constante(int a) {
//Inicializa todos os pixels do bloco com o valor a
int n, m;
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			pixel[n][m] = a;
		}
	}

}

void Bloco :: Copia(Bloco *A, int n0, int m0) {
//Copia um bloco de tamanho linhas x colunas de um bloco A a partir da posicao (n0, m0)
int n, m;
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			pixel[n][m] = A->pixel[n+n0][m+m0];
		}
	}

}


void Bloco :: ConcatenaLinhas(Bloco *A1, Bloco *A2) {
//Copia o conteudo de A1 e embaixo o conteudo de A2
int n, m;
	for(n = 0; n < linhas/2; n++) {
		for(m = 0; m < colunas; m++) {
			pixel[n][m] = A1->pixel[n][m];
			pixel[n+linhas/2][m] = A2->pixel[n][m];
		}
	}

}

void Bloco :: ConcatenaColunas(Bloco *A1, Bloco *A2) {
//Copia o conteudo de A1 e ao lado o conteudo de A2
int n, m;
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas/2; m++) {
			pixel[n][m] = A1->pixel[n][m];
			pixel[n][m+colunas/2] = A2->pixel[n][m];
		}
	}

}

void Bloco :: ExpandeLinhas(void) {
//Upsampling por 2 das linhas
int n, m;
Bloco *A;
	
	A = new Bloco;
	A->Dimensiona(linhas, colunas);
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			A->pixel[n][m] = pixel[n][m];
		}
	}

	for(n = 0; n < linhas; n++) delete [] pixel[n];	
	delete [] pixel;

	Dimensiona(2*A->linhas, A->colunas);
	for(n = 0; n < A->linhas; n++) {
		for(m = 0; m < A->colunas; m++) {
			pixel[2*n][m] = A->pixel[n][m];
			if(n < A->linhas-1) pixel[2*n+1][m] = (A->pixel[n][m]+A->pixel[n+1][m])/2;
			//else pixel[2*n+1][m] = 2*A->pixel[n][m] - pixel[2*n-1][m];
			//if(pixel[2*n+1][m] > 255) pixel[2*n+1][m] = 255;
			//if(pixel[2*n+1][m] < 0) pixel[2*n+1][m] = 0;
			else pixel[2*n+1][m] = A->pixel[n][m];
		}
	}

	delete A;
}

void Bloco :: ExpandeColunas(void) {
//Upsampling por 2 das colunas
int n, m;
	
Bloco *A;
	
	A = new Bloco;
	A->Dimensiona(linhas, colunas);
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			A->pixel[n][m] = pixel[n][m];
		}
	}

	for(n = 0; n < linhas; n++) delete [] pixel[n];	
	delete [] pixel;

	Dimensiona(A->linhas, 2*A->colunas);
	for(n = 0; n < A->linhas; n++) {
		for(m = 0; m < A->colunas; m++) {
			pixel[n][2*m] = A->pixel[n][m];
			if(m < A->colunas-1) pixel[n][2*m+1] = (A->pixel[n][m]+A->pixel[n][m+1])/2;
			//else pixel[n][2*m+1] = 2*A->pixel[n][m] - pixel[n][2*m-1];
			//if(pixel[n][2*m+1] > 255) pixel[n][2*m+1] = 255;
			//if(pixel[n][2*m+1] < 0) pixel[n][2*m+1] = 0;
			else pixel[n][2*m+1] = A->pixel[n][m];
		}
	}
	delete A;
}

void Bloco :: ContraiLinhas(void) {
//Downsampling por 2 das linhas
int n, m;
	
Bloco *A;
	
	A = new Bloco;
	A->Dimensiona(linhas, colunas);
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			A->pixel[n][m] = pixel[n][m];
		}
	}

	for(n = 0; n < linhas; n++) delete [] pixel[n];	
	delete [] pixel;

	Dimensiona(A->linhas/2, A->colunas);
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			pixel[n][m] = (A->pixel[2*n][m]+A->pixel[2*n+1][m])/2;
		}
	}
	delete A;
}

void Bloco :: ContraiColunas(void) {
//Downsampling por 2 das linhas
int n, m;
	
Bloco *A;
	
	A = new Bloco;
	A->Dimensiona(linhas, colunas);
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			A->pixel[n][m] = pixel[n][m];
		}
	}

	for(n = 0; n < linhas; n++) delete [] pixel[n];	
	delete [] pixel;

	Dimensiona(A->linhas, A->colunas/2);
	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			pixel[n][m] = (A->pixel[n][2*m]+A->pixel[n][2*m+1])/2;
		}
	}
	delete A;
}

char Bloco :: Compara(Bloco *A) {
//Compara com o bloco A. Se igual, retorna 0
int n, m;

	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			if(pixel[n][m] != A->pixel[n][m]) return(1);
		}
	}
	return(0);
}

void Bloco :: Escala(Bloco *A) {
//Copia os elementos de A porém aplicando uma transformação de escala de modo a adaptar a dimensão
 //short int exps1[BUFF];
 Bloco *EXP1;
 long int aux, q, n, m, aux1, aux2, es1, k;


 	if((linhas == A->linhas)&&(colunas == A->colunas)) {
   		for(n = 0; n < linhas; n++) {
     			for(m = 0; m < colunas; m++) {
       				pixel[n][m] = A->pixel[n][m];
    			 }
  		 }
  		 return;
	 }

	EXP1 = new Bloco;
	EXP1->Dimensiona(linhas, A->colunas);

	/* Processa as colunas */

    	aux = A->linhas-1;

    	for(q = 0; q < A->colunas; q++) {

		if((aux >= 0)&&(aux < linhas)) {

			for(n = 0; n < linhas; n++) {

		   		m = (aux*n)/(linhas);
				if(m+1 <= aux) aux1 = A->pixel[m+1][q];
				else aux1 = A->pixel[m][q];
				aux2 = A->pixel[m][q];
				es1 = ((aux1-aux2)*((aux*n)%(linhas)))/(linhas) + aux2;
				//exps1[n*A->colunas+q] = es1;
				EXP1->pixel[n][q] = es1;
			}
		}

		if(aux >= linhas) {

			for(n = 0; n < linhas; n++) {


				es1 = 0;
				for(k = -aux/2; k <= aux/2; k++) {
					m = (aux*n+k)/(linhas);
					if(m < 0) m = -m;
					if(m+1 <= aux) aux1 = A->pixel[m+1][q];
					else aux1 = A->pixel[m][q];
					aux2 = A->pixel[m][q];
					es1 = es1+((aux1-aux2)*((aux*n+k)%(linhas)))/(linhas) + aux2;
				}
				es1 = es1/(aux);
				//exps1[n*A->colunas+q] = es1;
				EXP1->pixel[n][q] = es1;
			}

		}

    	}
	/* Processa as linhas */
	
	aux = A->colunas-1;

   	for(q = 0; q < linhas; q++) {

		if((aux >= 0)&&(aux < colunas)) {

			for(n = 0; n < colunas; n++) {

		   		m = (aux*n)/(colunas);
				//if(m+1 <= aux) aux1 = exps1[q*Mi+m+1];
				//else aux1 = exps1[q*Mi+m];
				//aux2 = exps1[q*Mi+m];
				if(m+1 <= aux) aux1 = EXP1->pixel[q][m+1];
				else aux1 = EXP1->pixel[q][m];
				aux2 = EXP1->pixel[q][m];
				es1 = ((aux1-aux2)*((aux*n)%(colunas)))/(colunas) + aux2;
				pixel[q][n] = es1;
			}
		}

		if(aux >= colunas) {

			for(n = 0; n < colunas; n++) {


				es1 = 0;
				for(k = -aux/2; k <= aux/2; k++) {
					m = (aux*n+k)/(colunas);
					if(m < 0) m = -m;
					//if(m+1 <= aux) aux1 = exps1[q*Mi+m+1];
					//else aux1 = exps1[q*Mi+m];
					//aux2 = exps1[q*Mi+m];
					if(m+1 <= aux) aux1 = EXP1->pixel[q][m+1];
					else aux1 = EXP1->pixel[q][m];
					aux2 = EXP1->pixel[q][m];
					es1 = es1+((aux1-aux2)*((aux*n+k)%(colunas)))/(colunas) + aux2;
				}
				es1 = es1/(aux);
				pixel[q][n] = es1;
			}

		}

    	}	
	delete EXP1;
}

/********************************************************************************************/
/*				Metodos da classe Dicionario				    */
/********************************************************************************************/
Dicionario :: Dicionario(void) {

	TamSubDic = NULL;
	Elemento = NULL;
	linhas = colunas = NULL;
	NumEscalasV = NumEscalasH = 0;
}

void Dicionario :: Dimensiona(int NescV, int NescH, int Lmax, int Cmax) {
//Dimensiona o dicionario com Nescalas^2 subdicionarios de escalas diferentes e com dimensao maxima (na escala 0) Lmax x Cmax. Inicializa o tamanho de cada subdicionario com 0
int n, m, lin, col;
	//Aloca memoria para os elementos
	if(Elemento != NULL) {
	  for(n = 0; n < NescV; n++) {
	      for(m = 0; m < NescH; m++) free(Elemento[n][m]);
	      free(Elemento[n]);
	  }
	  free(Elemento);
	}
	if((Elemento = (Bloco ****)malloc(NescV*sizeof(Bloco ***))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	for(n = 0; n < NescV; n++) {
		if((Elemento[n] = (Bloco ***)malloc(NescH*sizeof(Bloco **))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}
	}
	//Neste ponto, Elemento[escalaV][escalaH] aponta o subdicionario correspondente a escala 

	NumEscalasV = NescV;
	NumEscalasH = NescH;

	//Aloca memoria para os tamanhos
	if(TamSubDic != NULL) {
		for(n = 0; n < NescV; n++) free(TamSubDic[n]);
		free(TamSubDic);
	}
	if((TamSubDic = (int **)malloc(NescV*sizeof(int *))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	for(n = 0; n < NescV; n++) {
		if((TamSubDic[n] = (int *)malloc(NescH*sizeof(int))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}
	}
	//neste ponto, TamSubDic[escalaV][escalaH] indica o tamanho do subdicionario correspondente

	//Faz todos os subdicionarios vazios
	for(n = 0; n < NescV; n++) {
		for(m = 0; m < NescH; m++) {
			TamSubDic[n][m] = 0;
		}
	}


	//Aloca memoria para a dimensoes
	if(linhas != NULL)  {
		for(n = 0; n < NescV; n++) free(linhas[n]);
		free(linhas);
	}
	if((linhas = (int **)malloc(NescV*sizeof(int *))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	for(n = 0; n < NescV; n++) {
		if((linhas[n] = (int *)malloc(NescH*sizeof(int))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}
	}
	if(colunas != NULL)  {
		for(n = 0; n < NescV; n++) free(colunas[n]);
		free(colunas);
	}
	if((colunas = (int **)malloc(NescV*sizeof(int *))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	for(n = 0; n < NescV; n++) {
		if((colunas[n] = (int *)malloc(NescH*sizeof(int))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}
	}

	//Calcula a dimensao correspondente a cada escala

	lin = Lmax;
	for(n = 0; n < NescV; n++) {
		col = Cmax;	
		for(m = 0; m < NescH; m++) {
			linhas[n][m] = lin;
			colunas[n][m] = col;
			col = col/2;
		}
		lin = lin/2;
	}
	//neste ponto, linhas[escala] e colunas[escala] indicam a dimensao dos elementos que compoe o subdicionario correspondente
  	for(n = 0; n < NumEscalasV; n++) {
		for(m = 0; m < NumEscalasH; m++) {
			printf("linhas[%d][%d] = %d, colunas[%d][%d] = %d\n", n, m, linhas[n][m], n, m, colunas[n][m]);
		}
	}
}

void Dicionario :: Inicializa(int inicio, int fim, int passo) {
//Inicializa dicionario com blocos constantes..
int n, m, k, Nelementos=0;

	//Calcula o numero de elementos a inserir
	for(n = inicio; n <= fim; n = n + passo) Nelementos++;
	
	for(n = 0; n < NumEscalasV; n++) {
		for(m = 0; m < NumEscalasH; m++) {
			//Aloca memoria para os elementos de cada subsubdicionario
			if((Elemento[n][m] = (Bloco **)malloc(Nelementos*sizeof(Bloco *))) == NULL) {
				printf("Erro na alocaçâo dinâmica de memória\n");
				exit(0);
			}
			//inicializa os elementos de cada subsubdicionario
			for(k = 0; k < Nelementos; k++) {
				Elemento[n][m][k] = new Bloco;
				Elemento[n][m][k]->Dimensiona(linhas[n][m], colunas[n][m]);
				Elemento[n][m][k]->Constante(inicio+passo*k);
			}
			TamSubDic[n][m] = Nelementos;
		}
	}	

}

char Dicionario :: IncluiSubDic(Bloco *A, int escalaV, int escalaH) {
//Inclui a matriz A no Subdicionario [escala] e retorna 0. Se o elemento ja estiver em algum dos subdicionarios da mesma escala, nao inclui e retorna 1
	int n, m;

	//Verifica se as dimensoes sao compativeis
	if((A->linhas != linhas[escalaV][escalaH])||(A->colunas != colunas[escalaV][escalaH])) {
		printf("Erro: dimensâo incompatível. Nâo foi possível incluir no subdicionaário.\n");
		exit(0);
	} 

	//Verifica se A ja pertence ao subdicionario
	for(n = 0; n < TamSubDic[escalaV][escalaH]; n++) {
		//if(Elemento[escala][n]->Distancia(A) == 0) return;
		if(Elemento[escalaV][escalaH][n]->Compara(A) == 0) return(1);
	}
	
	//aloca memoria para o novo elemento
	TamSubDic[escalaV][escalaH]++;
	if(TamSubDic[escalaV][escalaH] == 1) {
		if((Elemento[escalaV][escalaH] = (Bloco **)malloc(sizeof(Bloco *))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}

	}
	else {
		if((Elemento[escalaV][escalaH] = (Bloco **)realloc(Elemento[escalaV][escalaH], TamSubDic[escalaV][escalaH]*sizeof(Bloco *))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}
	}
	Elemento[escalaV][escalaH][TamSubDic[escalaV][escalaH]-1] = new Bloco;
	//Inclui elemento
	Elemento[escalaV][escalaH][TamSubDic[escalaV][escalaH]-1]->Dimensiona(A->linhas, A->colunas);
	for(n = 0; n < A->linhas; n++) {
		for(m = 0; m < A->colunas; m++) {
			Elemento[escalaV][escalaH][TamSubDic[escalaV][escalaH]-1]->pixel[n][m] = A->pixel[n][m];
		}
	}
	//Elemento[escala][TamSubDic[escala]-1]->Copia(A, 0, 0);
	//*Elemento[escala][TamSubDic[escala]-1] = A;

	return(0);
}

char Dicionario :: AdicionaAoSubDic(Bloco *A, int escalaV, int escalaH) {
//Inclui a matriz A no Subdicionario [escala] e retorna 0. Não testa se já está presente
	int n, m;

	//Verifica se as dimensoes sao compativeis
	if((A->linhas != linhas[escalaV][escalaH])||(A->colunas != colunas[escalaV][escalaH])) {
		printf("Erro: dimensâo incompatível. Nâo foi possível incluir no subdicionaário.\n");
		exit(0);
	} 
	
	//aloca memoria para o novo elemento
	TamSubDic[escalaV][escalaH]++;
	if(TamSubDic[escalaV][escalaH] == 1) {
		if((Elemento[escalaV][escalaH] = (Bloco **)malloc(sizeof(Bloco *))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}

	}
	else {
		if((Elemento[escalaV][escalaH] = (Bloco **)realloc(Elemento[escalaV][escalaH], TamSubDic[escalaV][escalaH]*sizeof(Bloco *))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}
	}
	Elemento[escalaV][escalaH][TamSubDic[escalaV][escalaH]-1] = new Bloco;
	//Inclui elemento
	Elemento[escalaV][escalaH][TamSubDic[escalaV][escalaH]-1]->Dimensiona(A->linhas, A->colunas);
	for(n = 0; n < A->linhas; n++) {
		for(m = 0; m < A->colunas; m++) {
			Elemento[escalaV][escalaH][TamSubDic[escalaV][escalaH]-1]->pixel[n][m] = A->pixel[n][m];
		}
	}
	//Elemento[escala][TamSubDic[escala]-1]->Copia(A, 0, 0);
	//*Elemento[escala][TamSubDic[escala]-1] = A;

	return(0);
}

/********************************************************************************************/
/*				Metodos da classe Imagem256				    */
/********************************************************************************************/
Imagem256 :: Imagem256(void) {

	linhas = colunas = 0;
	pixel = NULL;
}
void Imagem256 :: Dimensiona(int lin, int col) {
//dimensiona a imagem como uma matriz 2D
int n;
	if(pixel != NULL) {
	    for(n = 0; n < lin; n++) free(pixel[n]);
	    free(pixel);
	}
	if((pixel = (unsigned char **)malloc(lin*sizeof(unsigned char *))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	for(n = 0; n < lin; n++) {
		if((pixel[n] = (unsigned char *)malloc(col*sizeof(unsigned char))) == NULL) {
			printf("Erro na alocaçâo dinâmica de memória\n");
			exit(0);
		}
	}
	linhas = lin;
	colunas = col;
}

void Imagem256 :: CarregaDoArquivo(char *nome) {
	FILE *ifp;
	int n,m;

	if((ifp = fopen(nome, "rb")) == NULL) {
		printf("Erro na abertura do arquivo de entrada\n");
		exit(0);
	}

	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			pixel[n][m] = fgetc(ifp);
		}
	}
	fclose(ifp);
}
void Imagem256 :: GravaNoArquivo(char *nome) {
	FILE *ofp;
	int n,m;

	if((ofp = fopen(nome, "wb")) == NULL) {
		printf("Erro na abertura do arquivo de saida\n");
		exit(0);
	}

	for(n = 0; n < linhas; n++) {
		for(m = 0; m < colunas; m++) {
			fputc(pixel[n][m], ofp);
		}
	}
	fclose(ofp);
}

Bloco Imagem256 :: ExtraiBloco(int n, int m, int N, int M) {
	Bloco A;
	
	int k, l;
	
	A.Dimensiona(N,M);
	for(k = 0; k < N; k++) {
		for(l = 0; l < M; l++) {
			A.pixel[k][l] = pixel[k+n][l+m];
		}
	}
	return(A);
}

void Imagem256 :: EscreveBloco(Bloco A, int n, int m) {
	
	int k, l;
	
	for(k = 0; k < A.linhas; k++) {
		for(l = 0; l < A.colunas; l++) {
			pixel[k+n][l+m] = A.pixel[k][l];
		}
	}

}
