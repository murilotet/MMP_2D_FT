#include <stdio.h>
#include <stdlib.h>
#include "../DICIONARIO/Dicionario2D-FT.h"
#include "../CODARITMETICO/LONGLONG/Aritmetico.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>


class DecodificadorMMPRD2D {
public:
	Dicionario Dic;
	Imagem256 Ima;
	ArithmeticDecoder Decod;
	int Nmax, Mmax, ProfV, ProfH;
	void Decodifica(int N, int M, int linhas, int colunas, char *Entrada, char *Saida);
	void DesegmentaBloco(int escalaV, int escalaH, Bloco *AC);
};
void DecodificadorMMPRD2D :: Decodifica(int N, int M, int linhas, int colunas, char *Entrada, char *Saida) {

	int n, m;
	Bloco A;

	Nmax = N;
	Mmax = M;
	//Calcula o numero de escalas
	ProfV = 1;
	n = Nmax;
	while(n > 1) {
		ProfV++;
		n = n/2;
	}
	ProfH = 1;
	n = Mmax;
	while(n > 1) {
		ProfH++;
		n = n/2;
	}
	printf("Tamanho do bloco: %d x %d, número de escalas: %d x %d\n", Nmax, Mmax, ProfV,ProfH);


	//Inicializa o dicionario 
	Dic.Dimensiona(ProfV, ProfH, Nmax, Mmax);
	//Dic.Inicializa(0, 255, 1);
	Dic.Inicializa(0, 255, 1);

	//Prepara imagem de saida
	Ima.Dimensiona(linhas, colunas);

	//Inicializa contextos do codificador aritmetico
	for(n = 0; n < ProfV; n++) { 
		for(m = 0; m < ProfH; m++) { //2 contextos para cada escala: flags, e indices

			Decod.add_model();   //adiciona contexto para flags 2*escala)
			Decod.start_model(3, 2*(ProfH*n+m)); //Inicializa contexto com 3 simbolos

			Decod.add_model();   //adiciona contexto para indice (2*escala+1)
			Decod.start_model(Dic.TamSubDic[n][m], 2*(ProfH*n+m)+1); //Inicializa contexto com TamSubDic[n] simbolos
		}
	}

	//Abre arquivo de entrada
	Decod.set_input_file(Entrada);
	//Varre a imagem de entrada
	for(n = 0; n < Ima.linhas; n = n + Nmax) {
		for(m = 0; m < Ima.colunas; m = m + Mmax) {
			DesegmentaBloco(0, 0, &A);
			Ima.EscreveBloco(A, n, m);
			printf("Processou o bloco (%d %d)\n", n, m);
		}	
	}

	//Termina
	Decod.done_decoding();
	//Grava imagem no arquivo de saida
	Ima.GravaNoArquivo(Saida);

	printf("Dicionários:\n");
	for(n = 0; n < Dic.NumEscalasV; n++) {
		for(m = 0; m < Dic.NumEscalasH; m++) {
			printf("  TamSubDic[%d][%d] = %d:\n", n, m, Dic.TamSubDic[n][m]);
		}	
	}
}

void DecodificadorMMPRD2D :: DesegmentaBloco(int escalaV, int escalaH, Bloco *AC) {
	int k, i, n, m, n_thread;
	Bloco *AS, *AS1, *AR1, *AR2;
char fla;
	if((escalaV < ProfV-1)||(escalaH < ProfH-1)) {	
		k = Decod.decode_symbol(2*(escalaV*ProfH+escalaH)); //decodifica flag 
		Decod.update_model(k, 2*(escalaV*ProfH+escalaH));  //Atualiza tabelas do decodificador 
	}
	else k = 0;

	if(k == 0) {

		i = Decod.decode_symbol(2*(escalaV*ProfH+escalaH)+1); //Decodifica indice
		Decod.update_model(i, 2*(escalaV*ProfH+escalaH)+1);  //Atualiza tabelas
		AC->Dimensiona(Dic.linhas[escalaV][escalaH], Dic.colunas[escalaV][escalaH]);
		AC->Copia(Dic.Elemento[escalaV][escalaH][i], 0, 0);
		return;
	}

	if(k == 1) { 

		//Desegmenta na vertical
		AR1 = new Bloco;
		AR2 = new Bloco;
	
		DesegmentaBloco(escalaV+1, escalaH, AR1);
		DesegmentaBloco(escalaV+1, escalaH, AR2);

	}
	else { 

		//Desegmenta na horizontal
		AR1 = new Bloco;
		AR2 = new Bloco;
	
		DesegmentaBloco(escalaV, escalaH+1, AR1);
		DesegmentaBloco(escalaV, escalaH+1, AR2);

	}

	AC->Dimensiona(Dic.linhas[escalaV][escalaH], Dic.colunas[escalaV][escalaH]);
	if(AR1->linhas != Dic.linhas[escalaV][escalaH]) AC->ConcatenaLinhas(AR1, AR2);
	else AC->ConcatenaColunas(AR1, AR2);

	delete AR1;
	delete AR2;
	AS = new Bloco [32];
	#pragma omp parallel for private(n,m,n_thread)
	for(n = 0; n < Dic.NumEscalasV; n++) {
		for(m = 0; m < Dic.NumEscalasH; m++) {
			n_thread = omp_get_thread_num();
			AS[n_thread].Dimensiona(Dic.linhas[n][m], Dic.colunas[n][m]);
			AS[n_thread].Escala(AC);
			if(Dic.IncluiSubDic(&AS[n_thread], n, m) == 0) {
				Decod.add_new_char(2*(n*ProfH+m)+1);
			}
			AS[n_thread].Dimensiona(0,0);
		}
	}
	delete [] AS;	
	//Dic.Inclui(AC);
/*	if(Dic.IncluiSubDic(AC, escalaV, escalaH) == 0) {
		Decod.add_new_char(2*(escalaV*ProfH+escalaH)+1);
	}

	AS1 = new Bloco;
	AS1->Dimensiona(AC->linhas, AC->colunas);
	AS1->Copia(AC, 0, 0);
	for(m = 1; m <= escalaH; m++) {		
		AS1->ExpandeColunas();
		if(Dic.IncluiSubDic(AS1, escalaV, escalaH-m) == 0) {
			Decod.add_new_char(2*(escalaV*ProfH+escalaH-m)+1);
		}
	}
	delete AS1;	
	AS1 = new Bloco;
	AS1->Dimensiona(AC->linhas, AC->colunas);
	AS1->Copia(AC, 0, 0);
	for(m = escalaH+1; m < Dic.NumEscalasH; m++) {
		AS1->ContraiColunas();
		if(Dic.IncluiSubDic(AS1, escalaV, m) == 0) {
			Decod.add_new_char(2*(escalaV*ProfH+m)+1);
		}
	}	
	delete AS1;			

	AS = new Bloco;
	AS->Dimensiona(AC->linhas, AC->colunas);
	AS->Copia(AC, 0, 0);
	for(n = 1; n <= escalaV; n++) {  //Aqui temos escalaV-n
		AS->ExpandeLinhas();
		if(Dic.IncluiSubDic(AS, escalaV-n, escalaH) == 0) {
			Decod.add_new_char(2*((escalaV-n)*ProfH+escalaH)+1);
		}
		AS1 = new Bloco;
		AS1->Dimensiona(AS->linhas, AS->colunas);
		AS1->Copia(AS, 0, 0);
		for(m = 1; m <= escalaH; m++) {		
			AS1->ExpandeColunas();
			if(Dic.IncluiSubDic(AS1, escalaV-n, escalaH-m) == 0) {
				Decod.add_new_char(2*((escalaV-n)*ProfH+escalaH-m)+1);
			}
		}
		delete AS1;	
		AS1 = new Bloco;
		AS1->Dimensiona(AS->linhas, AS->colunas);
		AS1->Copia(AS, 0, 0);
		for(m = escalaH+1; m < Dic.NumEscalasH; m++) {
			AS1->ContraiColunas();
			if(Dic.IncluiSubDic(AS1, escalaV-n, m) == 0) {
				Decod.add_new_char(2*((escalaV-n)*ProfH+m)+1);
			}
		}	
		delete AS1;			
	}
	delete AS;	
	AS = new Bloco;
	AS->Dimensiona(AC->linhas, AC->colunas);
	AS->Copia(AC, 0, 0);
	for(n = escalaV+1; n < Dic.NumEscalasV; n++) {
		AS->ContraiLinhas();
		if(Dic.IncluiSubDic(AS, n, escalaH) == 0) {
			Decod.add_new_char(2*(n*ProfH+escalaH)+1);
		}
		AS1 = new Bloco;
		AS1->Dimensiona(AS->linhas, AS->colunas);
		AS1->Copia(AS, 0, 0);
		for(m = 1; m <= escalaH; m++) {		
			AS1->ExpandeColunas();
			if(Dic.IncluiSubDic(AS1, n, escalaH-m) == 0) {
				Decod.add_new_char(2*(n*ProfH+escalaH-m)+1);
			}
		}
		delete AS1;	
		AS1 = new Bloco;
		AS1->Dimensiona(AS->linhas, AS->colunas);
		AS1->Copia(AS, 0, 0);
		for(m = escalaH+1; m < Dic.NumEscalasH; m++) {
			AS1->ContraiColunas();
			if(Dic.IncluiSubDic(AS1, n, m) == 0) {
				Decod.add_new_char(2*(n*ProfH+m)+1);
			}
		}	
		delete AS1;			
	}	
	delete AS;*/	
	
}

int main(int argc, char **argv) {
	DecodificadorMMPRD2D MMP;

	int NMAX = 8;
	int MMAX = 8;
	int LINHAS=512;
	int COLUNAS=512;
	char ArquivoDeEntrada[256];
	char ArquivoDeSaida[256];
	struct timespec T0, T1;
	double accum;

	int n, m;

	for(n = 0; n < 256; n++) ArquivoDeEntrada[n] = 0;
	//Le parametros de entrada
	for(n = 0; n < argc; n++) {
		if(strcmp(argv[n], "-lin") == 0) LINHAS = atoi(argv[n+1]); 
		if(strcmp(argv[n], "-col") == 0) COLUNAS = atoi(argv[n+1]); 
		if(strcmp(argv[n], "-nmax") == 0) NMAX = atoi(argv[n+1]); 
		if(strcmp(argv[n], "-mmax") == 0) MMAX = atoi(argv[n+1]); 
		if(strcmp(argv[n], "-f") == 0) strcpy(ArquivoDeEntrada, argv[n+1]); 
	}
	if(strcmp(ArquivoDeEntrada, "") == 0) {
		printf("Favor especificar o arquivo de entrada (-f nomearq)\n");
		exit(0);
	}
	strcpy(ArquivoDeSaida, ArquivoDeEntrada);
	strcat(ArquivoDeSaida, ".rec");

	clock_gettime(CLOCK_REALTIME, &T0);

	MMP.Decodifica(NMAX, MMAX, LINHAS, COLUNAS, ArquivoDeEntrada, ArquivoDeSaida);
	clock_gettime(CLOCK_REALTIME, &T1);

	accum = ( T1.tv_sec - T0.tv_sec ) + ( T1.tv_nsec - T0.tv_nsec ) / 1000000000.0;
	
	printf( "Tempo de execuçâo: %lf segundos\n", accum );


}

