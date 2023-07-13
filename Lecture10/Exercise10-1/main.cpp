#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "../../random/random.h"
#include <cmath>
#include "punto.h"
#include <cassert>
#include <string>
#include <mpi.h>
//mpiexec -n 8 ./main.exe
using namespace std;



// DICHIARAZIONE DELLE VARIABILI GLOBALI

 const int G = 500; //generazioni
 const int M = 2000; //numero cromosomi
 const int n = 50; //cities
 const int G_migr = 50; //ogni quante generazioni migrare
 const int M_migr = 10; //quanti cromosomi scambiare
 const int Proc = 8; //numero processi

 int seed[4]; // per random gen
 Random rnd;


// DICHIARAZIONE DI FUNZIONI

// prepara randomgen genera vettore di città e prima generazione casuale
 void Initialize(bool mode, int rank, Punto* city, int** pop); 

// permuta vettore di n interi lasciando invariato il primo elemento
 void Permuta(int* v); //tot=n
 void Permuta(int* v, int tot);

// dà la lunghezza della strada per un cromosoma
 double Lenght(Punto* city, int* v);

// calcola il fitness di una popolazione e registra su fit e sum_fit
 void Fitness(Punto* city, int** pop, double* fit, double* sum_fit);

// ordina vettore fit in fit_ordered
 void Order(double* fit, double* fit_ordered);

// dà la media dei primi M/2 elementi di fit_ordered
 double Ave (double* fit_ordered);

// seleziona la label di un genitore
 int Selection(int** pop, double* fit, double* sum_fit, double* fit_ordered);

// attua crossing over e posiziona in pop_new
 void Crossover(int lavusca, int dizanoli, int label, int** pop, int** pop_new);

// muta in pop_new
 void Mutation(int label, int** pop_new);

// controlla se la popolazione è ben definita
 void Check(int** pop);

// stampa il resoconto della generazione
 void Report(int g, double* fit_ordered);

// dà la label del miglior percorso in una generazione
 int BestPoint(double* fit);

// dice se un intero è già contenuto in vettore di interi 
 bool is_there(int* v, int size, int new_int);



int main (int argc, char *argv[]){

	// DIVISIONE IN NODI

	 int size;
	 int rank;
	 MPI_Init(&argc, &argv);
	 MPI_Comm_size(MPI_COMM_WORLD, &size);
	 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	 cout << "Sono il nodo " << rank << " di " << size << endl;

	 int* exc_order = new int[Proc]; // ordine degli scambi
	 int* node_send = new int[Proc]; // vettore con nodo a cui mandare 
	 int* node_recv = new int[Proc]; // vettore con nodo da cui ricevere
	 // OSS: fungeranno da flag
   for(int i=0; i<Proc; i++) exc_order[i]=i;
	 int* to_send = new int[n]; // cromosoma da mandare
	 int* to_recv = new int[n]; // cromosoma da ricevere
	 int* pos_sent = new int[M_migr]; // vettore posizioni cromosomi mandati
	 MPI_Status* stat = new MPI_Status[Proc];
	 MPI_Request* req = new MPI_Request[Proc];

	
  // ALLOCAZIONE DINAMICA DELLE VARIABILI

	 Punto* city = new Punto[n]; // vettore di città
	 bool mode = 1; // quadrato
	
	 int** pop = new int*[M]; // popolazione
	 for(int i=0; i<M; i++) pop[i] = new int[n];
	 double* fit = new double[M]; 
	 double* sum_fit = new double;

	 int** pop_new = new int*[M]; // per registare la nuova pop
	 for(int i=0; i<M; i++) pop_new[i] = new int[n];

	 double* fit_ordered = new double[M];
	 int lavusca, dizanoli;
	 double lowest; 

	 ofstream output("best_individuals" + to_string(rank) + ".txt");
	 ofstream output2("path" + to_string(rank) + ".txt");
	
	 
  // INIZIALIZZO E OPERO SULLA GENERAZIONE ZERO

	 Initialize(mode, rank, city, pop);
	 Check(pop);
	 Fitness(city, pop, fit, sum_fit);
	 Order(fit, fit_ordered);
	 lowest = fit_ordered[0];
	 output << "0" <<  "\t" << fit_ordered[0] << "\t" << Ave(fit_ordered) << endl;
	 Report(0, fit_ordered);

	for(int i=0; i<n; i++){
		cout << "rank " << rank << ". City " << i << ") " << city[i].GetX() << endl;
		}

	// CICLO SU TUTTE LE GENERAZIONI
	
	 for(int g=1; g<G; g++){

		// CICLO SU TUTTI I CROMOSOMI
		 for(int i=0; i<M/2; i++){ 
		   lavusca = Selection(pop, fit, sum_fit, fit_ordered);
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
		 Fitness(city, pop, fit, sum_fit);
	   Order(fit, fit_ordered);
	   output << g << "\t" << fit_ordered[0] << "\t" << Ave(fit_ordered) << endl;
		 if(fit_ordered[0]<lowest) lowest = fit_ordered[0];
		 if(g<30 or g%50==0){
			 cout << "Rank " << rank << ": ";
			 Report(g, fit_ordered);
			 }

		// SCAMBIO FRA CONTINENTI: 
		// 1) IL NODO 0 DECIDE GLI ACCOPPIAMENTI
		// 2) LI COMUNICA A TUTTI
		// 3) AVVENGA LO SCAMBIO
		 if(g%G_migr==0){
			// 1) 
			 if(rank==0) Permuta(exc_order, Proc);
			// 2)
			 MPI_Bcast(&exc_order[0], Proc, MPI_INTEGER, 0, MPI_COMM_WORLD);
			// 3)
			 for(int i=0; i<M_migr; i++){
				// 3A: trovo la posiz del i-esimo migliore (senza ripetermi!)
				 double appo=fit_ordered[i];
	       pos_sent[i]=0;
	       while( fit[pos_sent[i]]!=appo  or  is_there(pos_sent,i, pos_sent[i]) )  pos_sent[i]++;
		    // 3B: copio il i-esimo cromosoma in to_send
         for(int j=0; j<n; j++) to_send[j] = pop[pos_sent[i]][j];
				// 3C: scambi (assegno flags, uso le MPI functions)
				 for(int j=0; j<Proc; j++){
					 if(rank==exc_order[j]){
						 if(j==0) node_recv[rank] = exc_order[Proc-1];
						 else node_recv[rank] = exc_order[j-1];
						 if(j==Proc-1) node_send[rank] = exc_order[0];
						 else node_send[rank] = exc_order[j+1];
					 }
				 } 
				 MPI_Isend(&to_send[0], n, MPI_INTEGER, node_send[rank], node_send[rank], MPI_COMM_WORLD, &req[rank]);
				 MPI_Recv(&to_recv[0], n, MPI_INTEGER, node_recv[rank], rank, MPI_COMM_WORLD, &stat[rank]);
				 
				// 3D: copio il cromosoma ricevuto al posto giusto
				 for(int j=0; j<n; j++) pop[pos_sent[i]][j] = to_recv[j];
			 }
		 }
		 
	 }
	
	// REGISTRO IL MIGLIOR PERCORSO PER L'ULTIMA GENERZIONE 
   int best = BestPoint(fit);
	 int label = 0;
   for(int i=0; i<n; i++){
		 label = pop[best][i]-1;
		 output2 << "\t" << city[label].GetX();
		 output2 << "\t" << city[label].GetY();
		 output2 << endl;
	 }
	 output2 << "\t" << city[pop[best][0]-1].GetX();
	 output2 << "\t" << city[pop[best][0]-1].GetY();
	 output2 << endl;
	


  // PULIZIA FINALE
	 for(int i=0; i<M; i++) delete[] pop[i];
	 delete[] pop;
	 delete[] fit;
	 delete[] fit_ordered;
	 delete[] city;
   output.close();
	 output2.close();

	 MPI_Finalize();

	
   return 0;

}

//_______________________________________________________________
// DEFINIZIONE DELLE FUNZIONI UTILIZZATE
//_______________________________________________________________


void Initialize(bool mode, int rank, Punto* city, int** pop){
	// Preparo random gen
	 int p1, p2;
   ifstream Primes("Primes");
	 for(int i=0; i<2*rank; i++) Primes >> p1 >> p2;
   Primes >> p1 >> p2 ;
 	 Primes.close();
	 ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
	// Genero vettore di città
	 /*double R=1;
	 double theta;
	 for(int i=0; i<n; i++){
	   if(mode){ 
		   city[i].SetPoint(rnd.Rannyu(0,R), rnd.Rannyu(0,R));
	   }
	   else{
		   theta = rnd.Rannyu(0,2*M_PI);
		   city[i].SetPoint(R*cos(theta), R*sin(theta));
	   }
   }*/
	 ifstream in;
	 in.open("American_capitals.dat", ios::in);
	 double x, y;
	 for(int i=0; i<n; i++){
		 in >> x >> y;
		 city[i].SetPoint(x, y);
	 }
	 in.close();
	
	// Preparo randomly la generazione zero
   for(int i=0; i<M; i++){
		 for(int j=0; j<n; j++) pop[i][j] = j+1;
	   Permuta(pop[i]);
	 }
}

void Permuta(int* v){
	 for(int i=1; i<n; i++){
	 int label = rnd.Dice(n-1); // numero da 1 a n-1
	 int appo = v[i];
	 v[i] = v[label];
	 v[label] = appo;
	}
}

void Permuta(int* v, int tot){
	 for(int i=1; i<tot; i++){
	 int label = rnd.Dice(tot-1); // numero da 1 a n-1
	 int appo = v[i];
	 v[i] = v[label];
	 v[label] = appo;
	}
}

double Lenght(Punto* city, int* v){
	 double l=0;
	 for(int i=0; i<n-1; i++)  l += (city[v[i]-1]-city[v[i+1]-1]).Modulo();
	 l += (city[v[0]-1]-city[v[n-1]-1]).Modulo();
	 return l;
}

void Fitness(Punto* city, int** pop, double* fit, double* sum_fit){
	 *sum_fit = 0;
	 for(int i=0; i<M; i++){
		 fit[i]=Lenght(city, pop[i]);
		 *sum_fit += pow(1./fit[i],3);
	 }
}

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

double Ave (double* fit_ordered){
	 double accu = 0;
	 for(int i=0; i<M/2; i++) accu+=fit_ordered[i];
	 return accu/(M/2);
}

int Selection(int** pop, double* fit, double* sum_fit, double* fit_ordered){

  double r=rnd.Rannyu(); //IMPLEMENTAZIONE ALTERNATIVA
	int j = floor(M*pow(r,5));
	double appo=fit_ordered[j];
	int conta=0;
	while(appo!=fit[conta]){
		conta++;
	} 
	return conta;	
}

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

void Mutation(int label, int** pop_new){
	double prob=0.1;
	//prima mutazione: single swap
	 if(rnd.Rannyu()<prob){
		 for(int i=0; i<10; i++){
			 int a = rnd.Dice(n-1); //pos 1
		   int b = rnd.Dice(n-1); // pos 2
		   swap( pop_new[2*label][a], pop_new[2*label][b] );
		 }
	 }
	 if(rnd.Rannyu()<prob){
		 for(int i=0; i<prob; i++){
			 int a = rnd.Dice(n-1); //pos 1
		   int b = rnd.Dice(n-1); // pos 2
		   swap( pop_new[2*label+1][a], pop_new[2*label+1][b] );
		 }
	 }
	//seconda mutazione: shift di m celle contigue
	 if(rnd.Rannyu()<prob){
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
	 if(rnd.Rannyu()<prob){
		 int m = rnd.Dice(n-1); // numero celle contigue
		 int* v = new int[n];
		 v[0]=1;
		 for(int i=1; i<n; i++){
			 if(i+m<n) v[i]=pop_new[2*label+1][i+m];
			 else v[i]=pop_new[2*label+1][i+m+1-n];
		 }
		 for(int i=0; i<n; i++) pop_new[2*label+1][i]=v[i];
		 delete[] v;
		 }
		// terza mutazione: permutazione celle contigue con altre
	 if(rnd.Rannyu()<prob){
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
   if(rnd.Rannyu()<prob){
		 int start = rnd.Dice(n-11);
		 int m = rnd.Dice(5);
		 int start2 = start + m + rnd.Dice(n-m-start-1-6);
		 int* v = new int[n];
		 v[0]=1;
		 for(int i=1; i<start; i++) v[i]=pop_new[2*label+1][i];
		 for(int i=start; i<start+m; i++) v[i-start+start2]=pop_new[2*label+1][i];
		 for(int i=start+m; i<start2; i++) v[i]=pop_new[2*label+1][i];
		 for(int i=start2; i<start2+m; i++) v[i-start2+start]=pop_new[2*label+1][i];
		 for(int i=start2+m; i<n; i++) v[i]=pop_new[2*label+1][i];
		 for(int i=0; i<n; i++) pop_new[2*label+1][i]=v[i];
		 delete[] v;
	 }
	// quarta mutazione: inversione
	 if(rnd.Rannyu()<prob){
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
	if(rnd.Rannyu()<prob){
		 int start = rnd.Dice(n-3);
		 int m = rnd.Dice(n-start-1);
		 int* v = new int[m];
		 v[0]=1;
		 for(int i=0; i<m; i++){
			 if(start+m-1-i<n) v[i]=pop_new[2*label+1][start+m-i-1];
			 else v[i]=pop_new[2*label+1][start+m-1-i-n+1];
		 }
		 for(int i=0; i<m; i++){
			 if(start+i<n) pop_new[2*label+1][start+i]=v[i];
			 else pop_new[2*label+1][start+i-n+1]=v[i];
		 }
	 }
}
	
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

void Report(int g, double* fit_ordered){
	cout << endl;
	cout << "GENEREATION " << g << " OF "<< G-1 << endl;
	cout << "Best fit: " << fit_ordered[0] << endl;
	cout << "______________________________________" << endl;
}

int BestPoint(double* fit){
	 int best = 0;
	 for(int i=1; i<M; i++){
		 if(fit[i]<fit[best]) best = i;
	 }
	 return best;
}


bool is_there(int* v, int size, int new_int){
	 bool check=0;
	 for(int i=0; i<size; i++){
		 if(v[i]==new_int) check=1;
	 }
	 return check;	
}