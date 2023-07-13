#include "../../random/random.h"
#include "../../mylibraries/symbols.h"
using namespace std;

int main(){

    intro();
    Random* myRand = new Random();
    myRand -> SetGenerator();

    //Generate exponentially distributed numbers
    int dim = 100000;
    double gamma_exp = 1;
    double gamma = 1;
    double mu = 0;

    //Print results
    const char* dirName = "results";
    CreateDirectory(dirName);
    ofstream file1("results/distributions.txt");
    if (!file1.is_open()) {
        cerr << "Error creating distributions.txt" << endl;
        return 1;
    }
    cout << "Printing vectors on file distributions.txt" << endl;
    for(int i = 0; i < dim; i++) {
        double exp = myRand -> Exponential(gamma_exp);
        double lor = myRand -> Lorentz(mu, gamma);
        file1 << exp << "\t" << lor << endl;
    }
    file1.close();

    //Testing Central Limit Theorem
    vector<string> files = {"results/uniform.txt", "results/exponential.txt","results/lorentzian.txt"};
    ofstream uniform(files[0]);
    ofstream exponential(files[1]);
    ofstream lorentzian(files[2]);
    vector<int> N = {1, 2, 10, 100};


    for(int j = 0; j < 10E4; j++) {
        double unif_acc = 0;
        double expo_acc = 0;
        double lor_acc = 0;
        for(int i = 0; i < 4; i++) {
            for(int k = 0; k < N[i]; k++) {
                unif_acc += myRand -> Rannyu();
                expo_acc += myRand -> Exponential(gamma_exp);
                lor_acc += myRand -> Lorentz(mu, gamma);
            }
            unif_acc /= N[i];
            expo_acc /= N[i];
            lor_acc /= N[i];
            uniform << unif_acc;
            exponential << expo_acc;
            lorentzian << lor_acc;

            if(i < 3) {
                uniform << "\t";
                exponential << "\t";
                lorentzian << "\t";
            } else {
                uniform << endl;
                exponential << endl;
                lorentzian << endl;
            }

            unif_acc = 0;
            expo_acc = 0;
            lor_acc = 0;

        }
    }

    uniform.close();
    exponential.close();
    lorentzian.close();

    cout << endl;
    return 0;
}