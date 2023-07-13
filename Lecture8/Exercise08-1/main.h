#include "../../random/random.h"
#include "../../mylibraries/symbols.h"

Random rnd;  

int nstep = 1e6;                      // numero totale di ripetizioni
int nblk = 100;                   // numero di blocchi
int steps = (int)nstep/nblk;          // numero di step in ogni blocco
double L;                   // lato del cubo visitabile / sigma per caso gaussiano
double span;
bool gauss;
int np_rend = 10000;
int mod = (int)nstep/np_rend;
int contastep = 0;

ofstream WriteResult;
ofstream WritePos;
ofstream WriteTraj;
ofstream HistoParam


double sgm = 1;
double mu = 1;

bool SA = 0;
int histofill_blk = 100;

// variabili di appoggio per SA
//double T = 10;          // regolabile da input.dat
double beta = 1.0;    // 1/T, avendo dichiarato una T iniziale prima
double db = 1;          // delta beta per ogni passaggio successivo di annealing. Regolabile da input.dat
int step = 0;           // contatore di step
int step_in_beta = 10;  // step per ogni temperatura

double deltaH = 1;
double errH = 0; // basta che sia < deltaH, per entrare nel ciclo
double oldH = 0;
double newH = 0;
double d_mu = 0;
double d_sgm = 0;
double Lmu = 1; // ampiezza iniziale dei salti di mu. Da impostare da input.dat
double Lsgm = 1; // idem per sgm

double H = 0;           
double var_prog_H = 0;  
double mean_prog_H = 0;

double x0; 
double y;  // posizione prima del passo
double x;  // posizione dopo il passo
double dx;

double q, A;
int attempted = 0, accepted = 0;


//=======================================================================Functions

double EvalH(double x) {
    double exp1 = exp(-(x-mu)*(x-mu)/(2*sgm*sgm)); 
    double exp2 = exp(-(x+mu)*(x+mu)/(2*sgm*sgm));

    double psi = exp1+exp2;
    double V = pow(x,4) - 2.5*pow(x,2);

    double Hpsi_kin = -0.5*( exp1*(pow(x-mu,2)/pow(sgm,4)-(pow(sgm,-2))) + exp2*(pow(x+mu,2)/pow(sgm,4)-(pow(sgm,-2))));

    return Hpsi_kin/psi + V;
}

double Psi2(double x) {
    return pow(exp(-(x-mu)*(x-mu)/(2*sgm*sgm))+exp(-(x+mu)*(x+mu)/(2*sgm*sgm)),2); 
}

void Input(void){

  ifstream ReadInput;

  cout << BOLDCYAN << "EXERCISE 8 - SIMULATED ANNEALING" << RESET << endl;

  ReadInput.open("input.dat");
    
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> L;
  ReadInput >> x0;
  ReadInput >> mu;
  ReadInput >> sgm;
  ReadInput >> SA;
  ReadInput >> histofill_blk;

  ReadInput >> beta;       
  ReadInput >> db;      
  ReadInput >> step_in_beta; 

  ReadInput >> Lmu;     // larghezza iniziale passi mu
  ReadInput >> Lsgm;    // idem per sigma

  rnd.SetGenerator();
  
  x = x0;
  steps = nstep/nblk;

  cout << "Total number of steps = " << nstep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << steps << endl << endl;
  cout << "Step lenght = " << L << endl;
  cout << "Initial position = " << x0 << endl;
  cout << "Num punti per riempire histo = " << histofill_blk*nblk << endl;

  ReadInput.close();

}

void MetroMove(){

    dx = {rnd.Rannyu(-L,L)};
    x = y + dx;
    attempted++;

    q = Psi2(x)/Psi2(y);
    A = min(1,q);

    if (A==1) // accetto direttamente
    {
        y = x;
        accepted++;
    }
    else {
        if(rnd.Rannyu() < A){
            y = x;
            accepted++;
        }
    }
}

void MeasureAverageH(){

    WriteResult.open("output/result.out");
    WritePos.open("output/pos.out");       //Queste servono per l'istogramma della funzione d'onda

    mean_prog_H = 0;
    var_prog_H = 0;

    for(int i = 0; i < nblk; i++)
    {     
        H = 0;
        attempted = 0;
        accepted = 0;

        for (int j = 0; j < steps; j++){

            MetroMove();
            H += EvalH(y);
            
            if(!SA && j%(steps/histofill_blk) == 0) WritePos << y << endl;
            
        }

        if(!SA){
            cout << "Block # " << i+1 << endl;
            cout << "Acceptance rate:   " << (double)accepted/attempted << endl;
            cout << "-----------------------------------" << endl;
        }

        H /= steps;   
        mean_prog_H += H;
        var_prog_H += H*H;

        if(!SA){
            if (WriteResult.is_open()){            
                if(i == 0) WriteResult << H << " " << mean_prog_H/(i+1) << " " << 0 << endl;
                else WriteResult << H << " " << mean_prog_H/(i+1) << " " << error(mean_prog_H/(i+1), var_prog_H/(i+1), i) << " " << endl;
            } else {
                if (!WriteResult.is_open()) cerr << "PROBLEM: Unable to open result.out" << endl;
            }  
        }
    }

    WriteResult.close();
    WritePos.close();
}
