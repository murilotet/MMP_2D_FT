
#include <stdio.h>
#include <stdlib.h>
#include "../DICIONARIO/Dicionario2D-FT.h"
#include "../CODARITMETICO/LONGLONG/Aritmetico.h"
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>

class Custo {
public:
	float Jmin;
	int Imin;
};

class CodificadorMMPRD2D {
public:
	Dicionario Dic;		//Armazena as matrizes do dicionario em multiplas escalas
	Imagem256 Ima;		//Armazena a imagem a ser comprimida
	ArithmeticCoder Cod;	//Interface com o arquivo de saida. Faz codificacao de entropia
	int Nmax, Mmax;		//Dimensao maxima dos blocos da imagem a serem processados 
	int ProfV, ProfH;	//Profundidade maxima da arvore de segmentacao
	float Lambda;		//Parametro de controle do ponto de operacao do algoritmo	
	char *CodSeg;		//String com o codigo de segmentacao (flags)
	int Pt_seg;		//Pointer para percorrer CodSeg
	Custo ****Cotimo;	//Mapa de índices para acelerar a execução
	int NumThreads;		//Numero de threads para execuçaõ paralela
	void Codifica(int N, int M, float Lam, int linhas, int colunas, char *Entrada, char *Saida, int n_threads);
	void SegmentaBloco(Bloco A, int escalaV, int escalaH, Bloco *AC);
	float ArvoreOtima(Bloco A, int escalaV, int escalaH, char *CodSeg0, int PosL, int PosC);
	Custo MenorJ_MP(Bloco A, int escalaV, int escalaH);
	Custo MenorJ(Bloco A, int escalaV, int escalaH);
	void CalculaMapa(Bloco A);
};
void CodificadorMMPRD2D :: Codifica(int N, int M, float Lam, int linhas, int colunas, char *Entrada, char *Saida, int n_threads) {
//Divide a imagem Ima em blocos de tamanho Nmax x Mmax e os processa sequencialmente
	int n, m, k;
	Bloco A, AC;

	NumThreads = omp_get_num_procs();
	if(n_threads < NumThreads) NumThreads = n_threads;

	Nmax = N;
	Mmax = M;
	Lambda = Lam;
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

	//aloca memoria para a string de segmentacao
	if((CodSeg = (char *)malloc(2*Nmax*Mmax*sizeof(char))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}

	//Inicializa o dicionario 
	Dic.Dimensiona(ProfV, ProfH, Nmax, Mmax);
	//Dic.Inicializa(0, 255, 1);
	Dic.Inicializa(0, 255, 1);

	//Aloca memoria para o mapa de indices
	Cotimo = new Custo*** [ProfV];
	for(n = 0; n < ProfV; n++) {
		Cotimo[n] = new Custo** [ProfH];
		for(m = 0; m < ProfH; m++) {
			Cotimo[n][m] = new Custo* [Nmax/Dic.linhas[n][m]];
			for(k = 0; k < Nmax/Dic.linhas[n][m]; k++) {
				Cotimo[n][m][k] = new Custo [Mmax/Dic.colunas[n][m]];
			}
		}
	}

	//Abre imagem de entrada
	Ima.Dimensiona(linhas, colunas);
	Ima.CarregaDoArquivo(Entrada);

	//Inicializa contextos do codificador aritmetico
	for(n = 0; n < ProfV; n++) { 
		for(m = 0; m < ProfH; m++) { //2 contextos para cada escala: flags, e indices

			Cod.add_model();   //adiciona contexto para flags (2*(NumEscalasH*escalaV+escalaH))
			Cod.start_model(3, 2*(ProfH*n+m)); //Inicializa contexto com 3 simbolos

			Cod.add_model();   //adiciona contexto para indices (2*escala+1)
			Cod.start_model(Dic.TamSubDic[n][m], 2*(ProfH*n+m)+1); //Inicializa contexto com  TamSubDic[n][m] simbolos
		}
	}

	//Abre arquivo de saida
	Cod.set_output_file(Saida);
	//Varre a imagem de entrada
	for(n = 0; n < Ima.linhas; n = n + Nmax) {
		for(m = 0; m < Ima.colunas; m = m + Mmax) {
			A = Ima.ExtraiBloco(n, m, Nmax, Mmax);
			CalculaMapa(A);
			ArvoreOtima(A, 0, 0, CodSeg, 0, 0);
			printf("CodSeg = %s\n", CodSeg);
			SegmentaBloco(A, 0, 0, &AC);
			A.Dimensiona(0, 0); //Libera memoria usada
			printf("Processou o bloco (%d %d)\n", n, m);
		}	
	}

	//Termina
	Cod.done_encoding();
	printf("Dicionários:\n");
	for(n = 0; n < Dic.NumEscalasV; n++) {
		for(m = 0; m < Dic.NumEscalasH; m++) {
			printf("  TamSubDic[%d][%d] = %d:\n", n, m, Dic.TamSubDic[n][m]);
		}	
	}


}
float CodificadorMMPRD2D :: ArvoreOtima(Bloco A,  int escalaV, int escalaH, char *CodSeg0, int PosL, int PosC) {
//Determina a melhor arvore de segmentacao para o bloco A, segundo um criterio R-D
	float J0, J1, J2, J3, J4, JV, JH;
	char *CodSeg1, *CodSeg2, *CodSeg3, *CodSeg4, V;
	Custo Cmin;
	Bloco A1, A2, A3, A4;

	if((escalaV == 0)&&(escalaH == 0)) {
		CodSeg0[0] = 0;  // CodSeg0 = ""
		Pt_seg = -1;
	}

	//Cmin = MenorJ(A, escalaV, escalaH);
	Cmin = Cotimo[escalaV][escalaH][PosL/Dic.linhas[escalaV][escalaH]][PosC/Dic.colunas[escalaV][escalaH]];
	J0 = Cmin.Jmin;		//Custo lagrangeano do nodo

	//if(Dic.linhas[escala]*Dic.colunas[escala] == 1) {
	//if((escalaV == ProfV-1)&&(escalaH == ProfH-1)) {
	if(A.linhas*A.colunas == 1) {
		strcat(CodSeg0, "M");
		return(J0);
	}


	if((CodSeg1 = (char *)malloc(2*A.linhas*A.colunas*sizeof(char))) == NULL) {
	//if((CodSeg1 = (char *)malloc(2*Nmax*Mmax*sizeof(char))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	if((CodSeg2 = (char *)malloc(2*A.linhas*A.colunas*sizeof(char))) == NULL) {
	//if((CodSeg2 = (char *)malloc(2*Nmax*Mmax*sizeof(char))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	if((CodSeg3 = (char *)malloc(2*A.linhas*A.colunas*sizeof(char))) == NULL) {
	//if((CodSeg1 = (char *)malloc(2*Nmax*Mmax*sizeof(char))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	if((CodSeg4 = (char *)malloc(2*A.linhas*A.colunas*sizeof(char))) == NULL) {
	//if((CodSeg2 = (char *)malloc(2*Nmax*Mmax*sizeof(char))) == NULL) {
		printf("Erro na alocaçâo dinâmica de memória\n");
		exit(0);
	}
	CodSeg1[0] = 0;		//CodSeg1 = ""
	CodSeg2[0] = 0;		//CodSeg2 = ""
	CodSeg3[0] = 0;		//CodSeg3 = ""
	CodSeg4[0] = 0;		//CodSeg4 = ""

	//Segmenta A em dois subblocos, A1 e A2 na direção vertical
	if(escalaV < ProfV-1) {
		A1.Dimensiona(Dic.linhas[escalaV+1][escalaH], Dic.colunas[escalaV+1][escalaH]);
		A2.Dimensiona(Dic.linhas[escalaV+1][escalaH], Dic.colunas[escalaV+1][escalaH]);

		A1.Copia(&A, 0, 0);
		A2.Copia(&A, A.linhas/2, 0);

		J1 = ArvoreOtima(A1, escalaV+1, escalaH, CodSeg1, PosL, PosC);
		J2 = ArvoreOtima(A2, escalaV+1, escalaH, CodSeg2, PosL+Dic.linhas[escalaV+1][escalaH], PosC);

		JV = J1+J2+Lambda*Cod.rate(1,2*(escalaV*ProfH+escalaH));

		A1.Dimensiona(0, 0);  //Libera a memoria alocada para A1
		A2.Dimensiona(0, 0);  //Libera a memoria alocada para A2
	}
	//Segmenta A em dois subblocos, A1 e A2 na direção horizontal
	if(escalaH < ProfH-1) {
		A3.Dimensiona(Dic.linhas[escalaV][escalaH+1], Dic.colunas[escalaV][escalaH+1]);
		A4.Dimensiona(Dic.linhas[escalaV][escalaH+1], Dic.colunas[escalaV][escalaH+1]);

		A3.Copia(&A, 0, 0);
		A4.Copia(&A, 0, A.colunas/2);

		J3 = ArvoreOtima(A3, escalaV, escalaH+1, CodSeg3, PosL, PosC);
		J4 = ArvoreOtima(A4, escalaV, escalaH+1, CodSeg4, PosL, PosC+Dic.colunas[escalaV][escalaH+1]);

		JH = J3+J4+Lambda*Cod.rate(2,2*(escalaV*ProfH+escalaH));

		A3.Dimensiona(0, 0);  //Libera a memoria alocada para A1
		A4.Dimensiona(0, 0);  //Libera a memoria alocada para A2
	}
	//Decide se é melhor segmentar na vertical ou na horizontal
	if(escalaV < ProfV-1) { //se tem JV
		if(escalaH < ProfH-1) { //se tem JH
	
			if(JV < JH) V = 1; //se JV < JH
			else V = 0; //Se JV > JH
		}
		else V = 1; //Se nao tem JH
	}
	else V = 0; //Se nao tem JV	
		
	J0 = J0+Lambda*Cod.rate(0,2*(escalaV*ProfH+escalaH));
	if(V == 1) J1 = JV;
	else J1 = JH;
	if(J0 <= J1) {
		strcat(CodSeg0, "M");
		free(CodSeg1);
		free(CodSeg2);
		free(CodSeg3);
		free(CodSeg4);
		return(J0);
	}
	else {
		if(V == 1) {
			strcat(CodSeg0, "V");
			strcat(CodSeg0, CodSeg1);
			strcat(CodSeg0, CodSeg2);
		}
		else {		
			strcat(CodSeg0, "H");
			strcat(CodSeg0, CodSeg3);
			strcat(CodSeg0, CodSeg4);
		}
		free(CodSeg1);
		free(CodSeg2);
		free(CodSeg3);
		free(CodSeg4);
		return(J1);
	}	
}

Custo CodificadorMMPRD2D :: MenorJ_MP(Bloco A,  int escalaV, int escalaH) {
//Escolhe o melhor elemento no dicionario para representar o bloco A, segundo um criterio R-D
	Custo Otimo;
	int indice, Nt, *Imin, n;
	float J, LF, *Jmin;

	Nt = NumThreads;
	Jmin = new float [Nt];
	Imin = new int [Nt];

	Imin[0] = 0;
	Jmin[0] = Dic.Elemento[escalaV][escalaH][0]->Distancia(A)+Lambda*(Cod.rate(0,2*(escalaV*ProfH+escalaH)+1));

	for(n = 1; n < Nt; n++) {
	    Imin[n] = Imin[0];
	    Jmin[n] = Jmin[0];
	}
	omp_set_num_threads(Nt);
	#pragma omp parallel for private(indice,LF,J,n)
	for(indice = 1; indice < Dic.TamSubDic[escalaV][escalaH]; indice++) { 
		LF = Lambda*(Cod.rate(indice,2*(escalaV*ProfH+escalaH)+1));
		n = omp_get_thread_num();
		if(LF < Jmin[n]) {
			J = Dic.Elemento[escalaV][escalaH][indice]->DistanciaF(A, Jmin[n]-LF)+LF;
			if(J < Jmin[n]) {
				Jmin[n] = J;
				Imin[n] = indice;
			}
		}
	}

	Otimo.Imin = Imin[0];
	Otimo.Jmin = Jmin[0];
	for(n = 1; n < Nt; n++) {
	      if(Jmin[n] < Otimo.Jmin) {
		  Otimo.Jmin = Jmin[n];
		  Otimo.Imin = Imin[n];
	      }
	}

	delete [] Jmin;
	delete [] Imin; 

/*	Otimo.Imin = 0;
	Otimo.Jmin = Dic.Elemento[escalaV][escalaH][0]->Distancia(A)+Lambda*(Cod.rate(0,2*(escalaV*ProfH+escalaH)+1));

	for(indice = 1; indice < Dic.TamSubDic[escalaV][escalaH]; indice++) { 
		LF = Lambda*(Cod.rate(indice,2*(escalaV*ProfH+escalaH)+1));
		if(LF < Otimo.Jmin) {
			J = Dic.Elemento[escalaV][escalaH][indice]->DistanciaF(A, Otimo.Jmin-LF)+LF;
			if(J < Otimo.Jmin) {
				Otimo.Jmin = J;
				Otimo.Imin = indice;
			}
		}
	}*/

	return(Otimo);
}

Custo CodificadorMMPRD2D :: MenorJ(Bloco A,  int escalaV, int escalaH) {
//Escolhe o melhor elemento no dicionario para representar o bloco A, segundo um criterio R-D
	Custo Otimo;
	int indice;
	float J, LF;

	Otimo.Imin = 0;
	Otimo.Jmin = Dic.Elemento[escalaV][escalaH][0]->Distancia(A)+Lambda*(Cod.rate(0,2*(escalaV*ProfH+escalaH)+1));

	for(indice = 1; indice < Dic.TamSubDic[escalaV][escalaH]; indice++) { 
		LF = Lambda*(Cod.rate(indice,2*(escalaV*ProfH+escalaH)+1));
		if(LF < Otimo.Jmin) {
			J = Dic.Elemento[escalaV][escalaH][indice]->DistanciaF(A, Otimo.Jmin-LF)+LF;
			if(J < Otimo.Jmin) {
				Otimo.Jmin = J;
				Otimo.Imin = indice;
			}
		}
	}
	return(Otimo);
}

void CodificadorMMPRD2D :: SegmentaBloco(Bloco A, int escalaV, int escalaH, Bloco *AC) {
//Faz a segmentacao, codificacao e atualizacao do dicionario
	Custo Cmin;
	Bloco *A1, *A2, *AR1, *AR2, *AS, *AS1;
	int n,m, n_thread;

	Pt_seg++;
	if(CodSeg[Pt_seg] == 'M') {

		Cmin = MenorJ_MP(A, escalaV, escalaH);
		if((escalaV < ProfV-1)||(escalaH < ProfH-1)) {
			Cod.encode_symbol(0, 2*(escalaV*ProfH+escalaH)); //Codifica flag 1 (se bloco > 1 x 1)
			Cod.update_model(0, 2*(escalaV*ProfH+escalaH));  //Atualiza tabelas do codificador 
		}
		Cod.encode_symbol(Cmin.Imin, 2*(escalaV*ProfH+escalaH)+1); //Codifica indice
		Cod.update_model(Cmin.Imin, 2*(escalaV*ProfH+escalaH)+1);  //Atualiza tabelas

		//return(*Dic.Elemento[escala][Cmin.Imin]);
		AC->Dimensiona(Dic.linhas[escalaV][escalaH], Dic.colunas[escalaV][escalaH]);
		AC->Copia(Dic.Elemento[escalaV][escalaH][Cmin.Imin], 0, 0);
		return;
	}

	if(CodSeg[Pt_seg] == 'V') {	
		Cod.encode_symbol(1, 2*(escalaV*ProfH+escalaH)); //Codifica flag 1
		Cod.update_model(1, 2*(escalaV*ProfH+escalaH));  //Atualiza tabelas do codificador 

		//segmenta na vertical
		A1 = new Bloco;
		A2 = new Bloco;

		A1->Dimensiona(Dic.linhas[escalaV+1][escalaH], Dic.colunas[escalaV+1][escalaH]);
		A2->Dimensiona(Dic.linhas[escalaV+1][escalaH], Dic.colunas[escalaV+1][escalaH]);
		A1->Copia(&A, 0, 0);
		A2->Copia(&A, A.linhas/2, 0);

		AR1 = new Bloco;
		AR2 = new Bloco;

		SegmentaBloco(*A1, escalaV+1, escalaH, AR1);
		SegmentaBloco(*A2, escalaV+1, escalaH, AR2);

		delete A1;
		delete A2;
	}
	else {
		Cod.encode_symbol(2, 2*(escalaV*ProfH+escalaH)); //Codifica flag 2
		Cod.update_model(2, 2*(escalaV*ProfH+escalaH));  //Atualiza tabelas do codificador 

		//segmenta na horizontal
		A1 = new Bloco;
		A2 = new Bloco;

		A1->Dimensiona(Dic.linhas[escalaV][escalaH+1], Dic.colunas[escalaV][escalaH+1]);
		A2->Dimensiona(Dic.linhas[escalaV][escalaH+1], Dic.colunas[escalaV][escalaH+1]);
		A1->Copia(&A, 0, 0);
		A2->Copia(&A, 0, A.colunas/2);

		AR1 = new Bloco;
		AR2 = new Bloco;

		SegmentaBloco(*A1, escalaV, escalaH+1, AR1);
		SegmentaBloco(*A2, escalaV, escalaH+1, AR2);

		delete A1;
		delete A2;
	}


	AC->Dimensiona(Dic.linhas[escalaV][escalaH], Dic.colunas[escalaV][escalaH]);
	if(AR1->linhas != Dic.linhas[escalaV][escalaH]) AC->ConcatenaLinhas(AR1, AR2);
	else AC->ConcatenaColunas(AR1, AR2);

	delete AR1;
	delete AR2;

	//Dic.Inclui(AC);
	AS = new Bloco [NumThreads];
	omp_set_num_threads(NumThreads);
	#pragma omp parallel for private(n,m,n_thread)
	for(n = 0; n < Dic.NumEscalasV; n++) {
		for(m = 0; m < Dic.NumEscalasH; m++) {
			n_thread = omp_get_thread_num();
			AS[n_thread].Dimensiona(Dic.linhas[n][m], Dic.colunas[n][m]);
			AS[n_thread].Escala(AC);
			if(Dic.IncluiSubDic(&AS[n_thread], n, m) == 0) {
				Cod.add_new_char(2*(n*ProfH+m)+1);
			}
			AS[n_thread].Dimensiona(0,0);
		}
	}
	delete [] AS;	
/*	AS = new Bloco;

	for(n = 0; n < Dic.NumEscalasV; n++) {
		for(m = 0; m < Dic.NumEscalasH; m++) {
			
			AS->Dimensiona(Dic.linhas[n][m], Dic.colunas[n][m]);
			AS->Escala(AC);
			if(Dic.IncluiSubDic(AS, n, m) == 0) {
				Cod.add_new_char(2*(n*ProfH+m)+1);
			}
			AS->Dimensiona(0,0);
		}
	}
	delete AS;*/	
	
}

void CodificadorMMPRD2D :: CalculaMapa(Bloco A) {

	int n,m,k,l,p, n_thread;
	Custo J;
	Bloco *B;

	B = new Bloco [NumThreads];
	//for(n = 0; n < Dic.NumEscalasV; n++) {
		//for(m = 0; m < Dic.NumEscalasH; m++) {
	omp_set_num_threads(NumThreads);
	#pragma omp parallel for private(p,k,l,n,m,n_thread) schedule(dynamic)
	for(p = 0; p < Dic.NumEscalasV*Dic.NumEscalasH; p++) {
			n = p/Dic.NumEscalasH;
			m = p%Dic.NumEscalasH;

			n_thread = omp_get_thread_num();
			B[n_thread].Dimensiona(Dic.linhas[n][m], Dic.colunas[n][m]);
			for(k = 0; k < Nmax/Dic.linhas[n][m]; k++) {
				for(l = 0; l < Mmax/Dic.colunas[n][m]; l++) {
					n_thread = omp_get_thread_num();
					B[n_thread].Copia(&A, k*Dic.linhas[n][m], l*Dic.colunas[n][m]);
					Cotimo[n][m][k][l] = MenorJ(B[n_thread], n, m);
				}
			}
			B[n_thread].Dimensiona(0,0);
		//}
	}
	delete [] B;
}


int main(int argc, char **argv) {
	CodificadorMMPRD2D MMP;

	int NMAX = 8;
	int MMAX = 8;
	int LINHAS=512;
	int COLUNAS=512;
	char ArquivoDeEntrada[256];
	char ArquivoDeSaida[256];
	char ArquivoDeStat[256];
	float LAMBDA = 30;
	FILE *ofp;
	char Estatisticas=0;
	int NUMTHREADS;
	struct timespec T0, T1;
	double accum;

	int n, m, k;

	NUMTHREADS = omp_get_num_procs();

	for(n = 0; n < 256; n++) ArquivoDeEntrada[n] = 0;
	//Le parametros de entrada
	for(n = 0; n < argc; n++) {
		if(strcmp(argv[n], "-lin") == 0) LINHAS = atoi(argv[n+1]); 
		if(strcmp(argv[n], "-col") == 0) COLUNAS = atoi(argv[n+1]); 
		if(strcmp(argv[n], "-nmax") == 0) NMAX = atoi(argv[n+1]); 
		if(strcmp(argv[n], "-mmax") == 0) MMAX = atoi(argv[n+1]); 
		if(strcmp(argv[n], "-lambda") == 0) LAMBDA = atof(argv[n+1]); 
		if(strcmp(argv[n], "-f") == 0) strcpy(ArquivoDeEntrada, argv[n+1]); 
		if(strcmp(argv[n], "-stat") == 0) Estatisticas = 1; 
		if(strcmp(argv[n], "-nt") == 0) NUMTHREADS = atoi(argv[n+1]); 
	}
	if(strcmp(ArquivoDeEntrada, "") == 0) {
		printf("Favor especificar o arquivo de entrada (-f nomearq)\n");
		exit(0);
	}
	strcpy(ArquivoDeSaida, ArquivoDeEntrada);
	//strcpy(strchr(ArquivoDeSaida, '.')  , ".mmp");
	strcat(ArquivoDeSaida, ".mmp");

	clock_gettime(CLOCK_REALTIME, &T0);

	MMP.Codifica(NMAX, MMAX, LAMBDA, LINHAS, COLUNAS, ArquivoDeEntrada, ArquivoDeSaida, NUMTHREADS);

	clock_gettime(CLOCK_REALTIME, &T1);

	accum = ( T1.tv_sec - T0.tv_sec ) + ( T1.tv_nsec - T0.tv_nsec ) / 1000000000.0;
	
	printf( "Tempo de execuçâo: %lf segundos\n", accum );

	if(Estatisticas == 1) {
		for(n = 0; n < MMP.Dic.NumEscalasV; n++) {
			for(m = 0; m < MMP.Dic.NumEscalasH; m++) {
				sprintf(ArquivoDeStat, "Estatisticas%dx%d.txt", NMAX/(1<<n), MMAX/(1<<m));
				ofp = fopen(ArquivoDeStat, "w");
				for(k = 0; k < MMP.Dic.TamSubDic[n][m]; k++) {
					fprintf(ofp, "%d %f\n",k,pow(2,-MMP.Cod.rate(k,2*(n*MMP.ProfH+m)+1)));
				}
				fclose(ofp);
			}
		}
	}

}

