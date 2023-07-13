#include "../../random/random.h"
#include "walker.h"
using namespace std;

double GenerateAngle(Random* myRand);
double* GenerateAngle3D(Random* myRand);

int main(){

    intro();
    Random* myRand = new Random();
    myRand -> SetGenerator();

    //int M = 1000;     numero totale esperimenti
    int N = 100;        //Numero di blocchi
    int steps = 100;    //passi random walk
    double a = 1.;

    Walker myWalker;

    vector<double> r2(N);                                
    vector<double> ave(N,0);
    vector<double> av2(N,0);
    vector<double> err_prog(N,0);

    const char* dirName = "results";
    CreateDirectory(dirName);
    vector<string> path = {"results/results_discrete.txt", "results/results_continue.txt", "results/discrete_walk.txt", "results/continue_walk.txt", "results/angle3d_distribution.txt"};
    ofstream file_1(path[0]);
    ofstream file_2(path[1]);

    //L'idea è salvare nel vettore r2 le distanze per ogni passo, per 100 walks alla volta
    //Una volta fatti 100 walks, queste informazioni vengono inserite nel vettore ave
    //Da questo vettore ave si procede a calcolare errore 

    for(int k = 0; k < N; k++) {
		r2 = vector<double>(steps,0.);
		for(int j = 0; j < N; j++) {                        //Ciclo sul singolo blocco, 100 walks
			myWalker.Reset();		
			for(int i = 0; i < steps; i++) {                //Ciclo sui passi del singolo walk
				r2[i]+=pow(myWalker.GetR(), 2);             //Salvo r² 
                myWalker.DiscreteStep(myRand, a);
			}
		}
		for(int i=0; i<steps; i++) r2[i]/=N;                            //Normalizzazione
		for(int i=0; i<steps; ++i) {
			ave[i]=( ave[i] * k + r2[i] )/((double)(k+1.) );   
			av2[i]=( av2[i] * k + r2[i]*r2[i])/((double)(k+1.) );               
		}
	}

    //Evaluating final error on √⟨r²〉
	for(int i=0; i<steps; ++i) err_prog[i] = sqrt((av2[i]-ave[i]*ave[i])/((double)N));      //Error on ⟨r²〉
	for(int i=0; i<steps; ++i){
		ave[i] = sqrt(ave[i]);                                                              
		err_prog[i] = 1./(2*ave[i])*err_prog[i];                                            //Error propagation
	    file_1 << i << "\t" << ave[i] << "\t" << err_prog[i] << endl;
	}

    //Evaluating a single walk
    myWalker.Reset();
    for(int i = 0; i < steps/2.; i++) {
        myWalker.PrintFile(path[2]);
        myWalker.DiscreteStep(myRand, a);
    }

    myWalker.Reset();
    std::fill(r2.begin(), r2.end(), 0);
    std::fill(ave.begin(), ave.end(), 0);
    std::fill(av2.begin(), av2.end(), 0);
    std::fill(err_prog.begin(), err_prog.end(), 0);        

    //Verifico riempimento uniforme della sfera
    ofstream file_5(path[4]);
    for(int i = 0; i < 1000; i++) {
        file_5 << GenerateAngle3D(myRand)[0] << "\t" << GenerateAngle3D(myRand)[1] << "\t";
        //Will also print a fallacious extraction in order to confront
        file_5 << myRand -> Rannyu(0, M_PI) << "\t" << myRand -> Rannyu(0, 2*M_PI) << endl;
    }

    //Ora ripeto il procedimento per il cammino continuo
    for(int k = 0; k < N; k++) {
		r2 = vector<double>(steps,0.);
		for(int j = 0; j < N; j++) {                        //Ciclo sul singolo blocco, 100 walks
			myWalker.Reset();		
			for(int i = 0; i < steps; i++) {                //Ciclo sui passi del singolo walk
				r2[i]+=pow(myWalker.GetR(), 2);             //Salvo r² 
                myWalker.ContinueStep(myRand, a);
			}
		}
		for(int i=0; i<steps; i++) r2[i]/=N;                            //Normalizzazione
		for(int i=0; i<steps; ++i) {
			ave[i]=( ave[i] * k + r2[i] )/((double)(k+1.) );   
			av2[i]=( av2[i] * k + r2[i]*r2[i])/((double)(k+1.) );               
		}
	}

    //Evaluating final error on √⟨r²〉
	for(int i=0; i<steps; ++i) err_prog[i] = sqrt((av2[i]-ave[i]*ave[i])/((double)N));      //Error on ⟨r²〉
	for(int i=0; i<steps; ++i){
		ave[i] = sqrt(ave[i]);                                                              
		err_prog[i] = 1./(2*ave[i])*err_prog[i];                                            //Error propagation
	    file_2 << i << "\t" << ave[i] << "\t" << err_prog[i] << endl;
	}

    //Evaluating a single walk
    myWalker.Reset();
    for(int i = 0; i < steps/2.; i++) {
        myWalker.PrintFile(path[3]);
        myWalker.ContinueStep(myRand, a);
    }
    


    cout << endl;
    return 0;
}

double GenerateAngle(Random* myRand) {

    double theta = 0;

    //In order to preserve uniform sampling, the angle will be generated in a square [-1,1]×[-1,1]
    double x = myRand -> Rannyu(-1,1);
    double y = myRand -> Rannyu(-1,1);
    while(x*x + y*y > 1) {                      //if the point lies outside the circle it needs to be recalculated
        x = myRand -> Rannyu(-1,1);
        y = myRand -> Rannyu(-1,1);
    }
    //In order to know the effective angle a distintion between quadrants has to be introduced
    //arccos(x) is the inverse function of cos(x) in the interval [0,π]
    if(y > 0) theta = acos(x/sqrt(x*x+y*y));
    else if  (y < 0) theta = M_PI + acos(x/sqrt(x*x+y*y));
    
    return theta;
}

double* GenerateAngle3D(Random* myRand) {

    double* angle = new double[2];
    double* position = new double[3];
    //Extract a point on the cube
    for(int i = 0; i < 3; i++) position[i] = myRand -> Rannyu(-1, 1);
    do {                      //if the point lies outside the circle it needs to be recalculated
        position[0] = myRand -> Rannyu(-1,1);
        position[1] = myRand -> Rannyu(-1,1);
        position[2] = myRand -> Rannyu(-1,1);
    }while(position[0]*position[0] + position[1]*position[1] + position[2]*position[2] > 1);

    //θ = arccos(z), considering that θ∈[0,π]
    angle[0] = acos(position[2]);

    if (position[0] != 0) angle[1] = atan2(position[1],position[0]);
    else if (position[1] > 0) angle[1] = M_PI / 2.0;
    else if (position[1] < 0) angle[1] = -M_PI / 2.0;
    else angle[1] = 0.0;

    return angle;
}