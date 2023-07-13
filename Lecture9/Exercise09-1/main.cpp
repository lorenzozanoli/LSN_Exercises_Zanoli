#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "../../random/random.h"
#include <cmath>
#include "punto.h"
#include <cassert>
#include <string>

using namespace std;


const int G = 1000; //generazioni
const int M = 5000; //numero cromosomi
const int n = 34; //cities

int seed[4];
Random rnd;

void Initialize(bool mode, Punto* cities, int** pop); 
void Permuta(int* v); 
double Lenght(Punto* cities, int* v);
void Fitness(Punto* cities, int** pop, double* fit, double* sum_fit);
void Order(double* fit, double* fit_ordered);
double Ave (double* fit_ordered);
int Selection(int** pop, double* fit, double* sum_fit, double* fit_ordered);
void Crossover(int lavusca, int dizanoli, int label, int** pop, int** pop_new);
void Mutation(int label, int** pop_new);
void Check(int** pop);
void Report(int g, double* fit_ordered);
int BestPoint(double* fit);



int main (int argc, char *argv[]){

	if(argc!=2){
		cerr << "Usage: ./main.exe mode" << endl;
		cerr << "mode=0: Circonferenza R=1" << endl;
		cerr << "mode=1: Quadrato l=1" << endl;
		return -1;
	}

	Punto* cities = new Punto[n]; 
	bool mode = atoi(argv[1]);
	
	int** pop = new int*[M]; //Popolazione = vettore di vettore di interi
	for(int i=0; i<M; i++) pop[i] = new int[n]; //Cromosoma = vettore di interi
	double* fit = new double[M]; //Vettore in cui metterò le lunghezze di percorso
	double* sum_fit = new double;

	int** pop_new = new int*[M]; // Nuova popolazione: generazione successiva
	for(int i=0; i<M; i++) pop_new[i] = new int[n];

	double* fit_ordered = new double[M];
	int lavusca, dizanoli;

	ofstream output("best_individuals.txt");
	ofstream output2("path.txt");
	
	 
    // INIZIALIZZO E OPERO SULLA GENERAZIONE ZERO

	Initialize(mode, cities, pop);
	Check(pop);
	Fitness(cities, pop, fit, sum_fit);
	Order(fit, fit_ordered);
	output << "0" <<  "\t" << fit_ordered[0] << "\t" << Ave(fit_ordered) << endl;
	Report(0, fit_ordered);


	// CICLO SU TUTTE LE GENERAZIONI
	for(int g=1; g<G; g++){

		// CICLO SU TUTTI I CROMOSOMI
		for(int i=0; i<M/2; i++){ 
		    lavusca = Selection(pop, fit, sum_fit, fit_ordered);    //restituisce l'indice
		    dizanoli = Selection(pop, fit, sum_fit, fit_ordered);  
            Crossover(lavusca, dizanoli, i, pop, pop_new);
			Mutation(i, pop_new);
        }

		// CONTROLLO POPOLAZIONE BUONA
		Check(pop_new);
		// SOSTITUISCO LA POP NUOVA ALLA VECCHIA
		for(int i=0; i<M; i++){
			delete[] pop[i];
			pop[i] = pop_new[i];
		}
		delete [] pop;
		pop = pop_new;
		// RI-CREO SPAZIO PER LA POP NUOVA
		pop_new = new int*[M]; // per registare la nuova pop
	    for(int i=0; i<M; i++) pop_new[i] = new int[n];

	  // NUOVO CALCOLO FITNESS
		Fitness(cities, pop, fit, sum_fit);
	    Order(fit, fit_ordered);
	    output << "\t" << g << "\t" << fit_ordered[0] << "\t" << Ave(fit_ordered) << endl;
		if(g<10 or g%50==0) Report(g, fit_ordered);
	}
	
	// REGISTRO IL MIGLIOR PERCORSO PER L'ULTIMA GENERZIONE 
    int best = BestPoint(fit);
	int label = 0;
    for(int i=0; i<n; i++){
		label = pop[best][i]-1;
		output2 << "\t" << cities[label].GetX();
		output2 << "\t" << cities[label].GetY();
		output2 << endl;
	}
	output2 << "\t" << cities[pop[best][0]-1].GetX();
	output2 << "\t" << cities[pop[best][0]-1].GetY();
	output2 << endl;
	
    // PULIZIA FINALE
	for(int i=0; i<M; i++) delete[] pop[i];
	delete[] pop;
	delete[] fit;
	delete[] fit_ordered;
	delete[] cities;
    output.close();
	output2.close();

	
   return 0;

}

//===========================FUNZIONI===============================

// prepara randomgen genera vettore di città e prima generazione casuale
void Initialize(bool mode, Punto* cities, int** pop){
	// Preparo random gen
	 int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
 	 Primes.close();
	 ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
	// Genero vettore di città
	 double R=1;
	 double theta;
	 for(int i=0; i<n; i++){
	   if(mode){ 
		   cities[i].SetPoint(rnd.Rannyu(0,R), rnd.Rannyu(0,R));
	   }
	   else{
		   theta = rnd.Rannyu(0,2*M_PI);
		   cities[i].SetPoint(R*cos(theta), R*sin(theta));
	   }
   }

	// Preparo randomly la generazione zero
   for(int i=0; i<M; i++){
		 for(int j=0; j<n; j++) pop[i][j] = j+1;
	   Permuta(pop[i]);
	 }
}

// permuta vettore di n interi lasciando invariato il primo elemento
void Permuta(int* v){
	 for(int i=1; i<n; i++){
	 int label = rnd.Dice(n-1); // numero da 1 a n-1
	 int appo = v[i];
	 v[i] = v[label];
	 v[label] = appo;
	}
}

// dà la lunghezza della strada per un cromosoma
double Lenght(Punto* cities, int* v){
	double l=0;
	for(int i=0; i<n-1; i++)  l += (cities[v[i]-1]-cities[v[i+1]-1]).Modulo2();
	l += (cities[v[0]-1]-cities[v[n-1]-1]).Modulo2();

	return l;
}

// calcola il fitness di una popolazione e registra su fit e sum_fit
void Fitness(Punto* cities, int** pop, double* fit, double* sum_fit){
	 *sum_fit = 0;
	 for(int i=0; i<M; i++){
		 fit[i]=Lenght(cities, pop[i]);
		 *sum_fit += pow(1./fit[i],3);
	 }
}

//ordine fit e lo mette in fit_ordered
void Order(double* fit, double* fit_ordered){
	 for(int i=0; i<M; i++) fit_ordered[i]=fit[i]; //copio fit in fit_ordered
	 int pos=0;
	 double appo=0;
   for(int i=0; i<M-1; i++){
		 pos=i;
		 for(int j=i+1; j<M; j++)
		 	 if(fit_ordered[j]<fit_ordered[pos]) pos=j;
	   appo = fit_ordered[i];
		 fit_ordered[i] = fit_ordered[pos];
		 fit_ordered[pos] = appo;
	 }
}

// dà la media dei primi M/2 elementi di fit_ordered
double Ave (double* fit_ordered){
	 double accu = 0;
	 for(int i=0; i<M/2; i++) accu+=fit_ordered[i];
	 return accu/(M/2);
}

// seleziona la label di un genitore
int Selection(int** pop, double* fit, double* sum_fit, double* fit_ordered){

  double r=rnd.Rannyu();
	int j = floor(M*pow(r,5));
	double appo=fit_ordered[j];
	int conta=0;
	while(appo!=fit[conta]){
		conta++;
	} 
	return conta;	
}

// attua crossing over e posiziona in pop_new
void Crossover(int lavusca, int dizanoli, int label, int** pop, int** pop_new){
	 if(rnd.Rannyu()<0.65){ // avviene crossover
		 int pos = rnd.Dice(n); //num da 1 a n
		 //cout << "taglio in posizione " << pos << endl;
		 //cout << pos << " numeri buoni a sx" << endl;
		 //cout << n-pos << " numeri mancanti a dx" << endl;
		 int pos_copy = pos;
		 bool check = 0;
		 for(int i=0; i<pos; i++){ //ricopio a sx del taglio
			 pop_new[2*label][i] = pop[lavusca][i];
			 pop_new[2*label+1][i] = pop[dizanoli][i];
		 }

		 //cout << "Inizio setaccio: " << endl;
		 for(int i=1; i<n; i++){ //riempio a dx il primo figlio
			 check=0;
			 //cout << i << ") ";
			 for(int j=1; j<pos; j++){
				 if(pop[dizanoli][i]==pop_new[2*label][j]){
					 check=1;
					 //cout << "Il  " << pop[dizanoli][i] << " c'è" << endl;
				 }
			 }	 
			 if(!check){
				 pop_new[2*label][pos] = pop[dizanoli][i];
				// cout << "manca il " << pop[dizanoli][i] << " in pos " << pos << endl;
				 pos++;
				 //cout <<"incremento: pos=" << pos << endl;	 
			 }
		 }
		 //cerr << "fin qui OK";
		 //cerr << "Inizio setaccio: " << endl;
		 for(int i=1; i<n; i++){ // riempio a dx il secondo figlio
			 check=0;
			 //cerr << i << ") ";
			 for(int j=1; j<pos_copy; j++){
				 if(pop[lavusca][i]==pop_new[2*label+1][j]){
					 check=1;
					 //cerr << "Il  " << pop[lavusca][i] << " c'è" << endl;
				 }
			 }
			 if(!check){
				 pop_new[2*label+1][pos_copy] = pop[lavusca][i];
				 //cout << "manca il " << pop[lavusca][i] << " in pos " << pos_copy << endl;
				 pos_copy++;
				 //cout <<"incremento: pos_copy=" << pos_copy << endl;	
			 }
		 }
		 if(pos!=n or pos_copy!=n) cout << "Errore in crossover" << endl;
	 }
	 else{ // non avviene crossover
		 for(int i=0; i<n; i++){
			 pop_new[2*label][i] = pop[2*label][i];
			 pop_new[2*label+1][i] = pop[2*label+1][i];
		}
	}
}

// muta la popolazione e la mette in pop_new
void Mutation(int label, int** pop_new){
	
    //MUTAZIONE 1: Scambio singolo (ne tento 5)
	if(rnd.Rannyu()<0.1){
		for(int i=0; i<5; i++){
			int a = rnd.Dice(n-1); //pos 1
		    int b = rnd.Dice(n-1); // pos 2
		   swap( pop_new[2*label][a], pop_new[2*label][b] );
		}
	}

	//MUTAZIONE 2: Shift di celle 
	if(rnd.Rannyu()<0.1){
		int m = rnd.Dice(n-1); // numero celle contigue
		int* v = new int[n];
		v[0]=1;
		for(int i=1; i<n; i++){
			if(i+m<n) v[i]=pop_new[2*label][i+m];
			else v[i]=pop_new[2*label][i+m+1-n];
		}
		for(int i=0; i<n; i++) pop_new[2*label][i]=v[i];
		delete[] v;
	}

	//MUTAZIONE 3: Permutazione di celle
	if(rnd.Rannyu()<0.1){
		 int start = rnd.Dice(n-11);
		 int m = rnd.Dice(5);
		 int start2 = start + m + rnd.Dice(n-m-start-1-6);
		 int* v = new int[n];
		 v[0]=1;
		 for(int i=1; i<start; i++) v[i]=pop_new[2*label][i];
		 for(int i=start; i<start+m; i++) v[i-start+start2]=pop_new[2*label][i];
		 for(int i=start+m; i<start2; i++) v[i]=pop_new[2*label][i];
		 for(int i=start2; i<start2+m; i++) v[i-start2+start]=pop_new[2*label][i];
		 for(int i=start2+m; i<n; i++) v[i]=pop_new[2*label][i];
		 for(int i=0; i<n; i++) pop_new[2*label][i]=v[i];
		 delete[] v;
	}

	//MUTAZIONE 4: Inversione
	if(rnd.Rannyu()<0.1){
		 int m = rnd.Dice(n-1);
		 int start = rnd.Dice(n-1);
		 int* v = new int[m];
		 v[0]=1;
		 for(int i=0; i<m; i++){
			 if(start+m-1-i<n) v[i]=pop_new[2*label][start+m-i-1];
			 else v[i]=pop_new[2*label][start+m-1-i-n+1];
		 }
		 for(int i=0; i<m; i++){
			 if(start+i<n) pop_new[2*label][start+i]=v[i];
			 else pop_new[2*label][start+i-n+1]=v[i];
		 } 
		 delete[]v;
	}
}

// controlla se la popolazione è ben definita
void Check(int** pop){
	 bool check=0;
	 for(int i=0; i<M; i++){
		 if(pop[i][0]!=1){
			 cout << "Error in creating new gen: doesn't start with 1" << endl;
			 exit(1);
		 }
		 for(int num=2; num<n+1; num++){
			 check=0;
			 for(int j=1; j<n; j++){
				 if(pop[i][j]==num) check=1;
			 }
			 if(!check){
				 cout << "Error in creating new gen: absent " << num << endl;
				 exit(1);
			 }
		 }
	 }
} 

// stampa il resoconto della generazione
void Report(int g, double* fit_ordered){
	cout << endl;
	cout << "GENEREATION " << g << " OF "<< G-1 << endl;
	cout << "Best fit: " << fit_ordered[0] << endl;
	cout << "______________________________________" << endl;
}

// dà la label del miglior percorso in una generazione
int BestPoint(double* fit){
	 int best = 0;
	 for(int i=1; i<M; i++){
		 if(fit[i]<fit[best]) best = i;
	 }
	 return best;
}