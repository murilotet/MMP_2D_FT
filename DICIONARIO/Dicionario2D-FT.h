/********************************************************************************************/
/*			Classes para o dicionario de matrizes				    */
/********************************************************************************************/

class Bloco {  //Define um bloco como uma matriz
public:
	int linhas, colunas;  //dimensoes da matriz
	int **pixel;          //elementos da matriz
	Bloco(void);
	void Dimensiona(int lin, int col); //Aloca memoria
	float Distancia(Bloco A); //Calcula distancia (erro quadratico) para o bloco A
	float DistanciaF(Bloco A, float Da); //Calcula distancia para A, abortando se D > Da
	void Constante(int a); //Inicializa todos os elementos com a
	void Copia(Bloco *A, int n0, int m0); //Copia um bloco de A a partir de (n0, m0)
	void ConcatenaLinhas(Bloco *A1, Bloco *A2);  //
	void ConcatenaColunas(Bloco *A1, Bloco *A2);  //
	void ExpandeLinhas(void);
	void ExpandeColunas(void);
	void ContraiLinhas(void);
	void ContraiColunas(void);
	char Compara(Bloco *A);
	void Escala(Bloco *A); //Copia A transformando para dimensao corrente
};


class Dicionario { //Dicionario de matrizes composto por subdicionarios particionados
public:
	int NumEscalasV, NumEscalasH;	          //Numero de escalas distintas em cada direcao
	int **TamSubDic;           //Numero de elementos de cada subdicionario
	Bloco ****Elemento;        //ponteiro para o bloco (Elemento[escalaV][escalaH][indice])
 	int **linhas, **colunas;     //Dimensao dos elementos dos subdicionarios a cada escala
	Dicionario(void);
	void Dimensiona(int NescV, int NescH, int Lmax, int Cmax);
	void Inicializa(int inicio, int fim, int passo);
	void Inclui(Bloco A);
	char IncluiSubDic(Bloco *A, int escalaV, int escalaH);
	char AdicionaAoSubDic(Bloco *A, int escalaV, int escalaH);
};

class Imagem256 {
public:
	int linhas, colunas;
	unsigned char **pixel;
	Imagem256(void);
	void Dimensiona(int lin, int col);
	void CarregaDoArquivo(char *nome);
	void GravaNoArquivo(char *nome);
	Bloco ExtraiBloco(int n, int m, int N, int M);
	void EscreveBloco(Bloco A, int n, int m);
};
