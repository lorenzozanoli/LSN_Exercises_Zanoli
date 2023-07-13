/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int i = 0; i < 2000; i++) Move(metro);
  
  ofstream spin;
  spin.open("spin.dat");

  for(int i = 0; i < 3000; i++) Move(metro);
  
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    for (int i=0; i<nspin; ++i) spin << s[i] << "  ";
    spin << endl;
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  spin.close();
  return 0;
}

void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;

  ReadInput >> config;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  cout << "Initial configuration:" << endl;
  if(!restart) {
    for (int i=0; i<nspin; ++i){
      if(config == 0) {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
      } else s[i] = -1;
    }
  } else {
    ifstream spin;
    spin.open("config.final");
    for (int i=0; i<nspin; ++i) {
      spin >> s[i];
    }
    spin.close();

  }
  //for(int i = 0; i < nspin; i++) cout << s[i] << endl;
  cout << endl << endl;

  
  //Evaluate energy etc. of the initial configuration
    Measure();

  //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Move(int metro)
{
  int o;  //Seleziona la particella da flipapre
  //double p, energy_old, energy_new, sm;
  //double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //Metropolis
      int new_spin = - s[o];     //Inverto lo spin presentw
      attempted++;
      double DeltaEnergy = 2*Boltzmann(new_spin, o);    //Calcola la differenza di energia 
      //Ora si tratta di accettare o no la mossa
      double A = min(1, exp(- beta * DeltaEnergy));
      if(A == 1) {              //Accetto automaticamente se DeltaEnergy < 0
        s[o] = new_spin;
        accepted++;
      } else if (rnd.Rannyu()  < A) {         //Se DeltaEnergy > 0 non rifiuto subito, dò un'altra chance
        s[o] = new_spin;
        accepted++;
      }
      //Qui se non si verifica nessuno dei due if allora la mossa non viene accettata
      //In Gibbs questa possibilità non è contemplata
    }
    else //Gibbs sampling
    {
      attempted++;
      accepted++;       //So già che verrà accettata prima ancora di farla
      int new_spin = pow(1, rnd.Rannyu(0, 2));    //Riscelgo un nuovo spin a caso
      double DeltaEnergy = 2*Boltzmann(new_spin, o); 
      double prob = 1./(1 + exp(beta * DeltaEnergy));   //Probabilità dello scambio
      if(rnd.Rannyu() < prob) s[o] = new_spin;
      else s[o] = - new_spin;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  // int bin;
  double H = 0.0, u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i) 
  {
    H = -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    u += H; 
    m += s[i];
  } 
  // energy
  walker[iu] = u; // La media si fa in Averages
  walker[ic] = u*u; // 
  walker[im] = m;                 
  walker[ix] = beta*(m*m); // non beta*(x - m*m), perché h = 0  
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
   ofstream Ene, Heat, Mag, Chi;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    //Adesso i vettori blk_av contengono le medie a blocchi

    if(h == 0) {
      //Energia Interna
      Ene.open("output/output_ene_0.dat",ios::app);
      stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
      glob_av[iu]  += stima_u;
      glob_av2[iu] += stima_u*stima_u;
      err_u=Error(glob_av[iu],glob_av2[iu],iblk);
      Ene << iblk <<  "\t"<< stima_u << "\t"<< glob_av[iu]/(double)iblk << "\t"<< err_u << endl;
      Ene.close();

      //Capacità Termica
      Heat.open("output/output_heat_0.dat",ios::app);
      stima_c = beta*beta * (blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2))/(double)nspin; 
      glob_av[ic]  += stima_c;
      glob_av2[ic] += stima_c*stima_c;
      err_c=Error(glob_av[ic],glob_av2[ic],iblk);
      Heat << iblk <<  "\t"<< stima_c << "\t"<< glob_av[ic]/(double)iblk << "\t"<< err_c << endl;
      Heat.close();

      //Suscettività Magnetica
      Chi.open("output/output_chi_0.dat",ios::app);
      stima_x = blk_av[ix]/blk_norm/(double)nspin;
      glob_av[ix]  += stima_x;
      glob_av2[ix] += stima_x*stima_x;
      err_x=Error(glob_av[ix],glob_av2[ix],iblk);
      Chi << iblk <<  "\t"<< stima_x << "\t"<< glob_av[ix]/(double)iblk << "\t"<< err_x << endl;
      Chi.close();
    } else {
      //Magnetizzazione
      Mag.open("output/output_mag_0.dat",ios::app);
      stima_m = blk_av[im]/blk_norm/(double)nspin;
      glob_av[im]  += stima_m;
      glob_av2[im] += stima_m*stima_m;
      err_m=Error(glob_av[im],glob_av2[im],iblk);
      Mag << iblk <<  "\t"<< stima_m << "\t"<< glob_av[im]/(double)iblk << "\t"<< err_m << endl;
      Mag.close();
    }

    cout << "----------------------------" << endl << endl;
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double min(double a, double b){
  if(a < b) return a;
  else return b;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
