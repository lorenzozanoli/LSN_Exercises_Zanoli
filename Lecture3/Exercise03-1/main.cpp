#include "../../random/random.h"
#include "../../mylibraries/symbols.h"
using namespace std;

struct parameter {
    double S_0 = 100;
    double T = 1;
    double K = 100;
    double r = 0.1;
    double sigma = 0.25;
};

double max(double a, double b);
double Call(Random* myRand, parameter param, bool flag);
double Pull(Random* myRand, parameter param, bool flag);

int main(){

    intro();
    Random* myRand = new Random();
    myRand -> SetGenerator();

    parameter param;
    param.S_0 = 100;
    param.T = 1;
    param.K = 100;
    param.r = 0.1;
    param.sigma = 0.25;

    int M = 100000;                                     //Total number of throws
    int N = 100;                                        //Number of blocks
    int L = int(M/N);                                   //Number of throws in each block

    vector<double> x(N);                                //Block counter
    vector<double> ave(N,0);
    vector<double> av2(N,0);
    vector<double> sum_prog(N,0);
    vector<double> sum_2_prog(N,0);
    vector<double> err_prog(N,0);
    std::iota(x.begin(), x.end(), 0);                   //Fills with 0,1,2,3...

    const char* dirName = "results";
    CreateDirectory(dirName);
    vector<string> path = {"results/results_call_discrete.txt", "results/results_pull_discrete.txt", "results/results_call_onestep.txt", "results/results_pull_onestep.txt"};

    //Evaluating call option with one step
    for(int i = 0; i < N; i++) {                        //Cycle on all blocks
        double sum1 = 0;
        for(int j = 0; j < L; j++) {                    //Cycle on the single block
            sum1 += Call(myRand, param, true);
        }
        //Normalising mean and variance
        ave[i] = sum1/L;
        av2[i] = (ave[i])*(ave[i]);
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){       
            sum_prog[i] += ave[j];
            sum_2_prog[i] += av2[j];
        }
        sum_prog[i]/=(i+1);
        sum_2_prog[i]/=(i+1);
        if(i == 0) {
            err_prog[i] = 0;
        } else err_prog[i] = sqrt((sum_2_prog[i] - sum_prog[i]*sum_prog[i])/double(i));
    }
    PrintResults(path[2], x, sum_prog, err_prog);

    std::fill(ave.begin(), ave.end(), 0);
    std::fill(av2.begin(), av2.end(), 0);
    std::fill(sum_prog.begin(), sum_prog.end(), 0);
    std::fill(err_prog.begin(), err_prog.end(), 0);
    std::fill(sum_2_prog.begin(), sum_2_prog.end(), 0);

    //Evaluating pull option with one step
    for(int i = 0; i < N; i++) {                            //Cycle on all blocks
        double sum1 = 0;
        for(int j = 0; j < L; j++) {                        //Cycle on the single block
            sum1 += Pull(myRand, param, true);
        }
        //Normalising mean and variance
        ave[i] = sum1/L;
        av2[i] = (ave[i])*(ave[i]);
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){       
            sum_prog[i] += ave[j];
            sum_2_prog[i] += av2[j];
        }
        sum_prog[i]/=(i+1);
        sum_2_prog[i]/=(i+1);
        if(i == 0) {
            err_prog[i] = 0;
        } else err_prog[i] = sqrt((sum_2_prog[i] - sum_prog[i]*sum_prog[i])/double(i));
    }
    PrintResults(path[3], x, sum_prog, err_prog);

    std::fill(ave.begin(), ave.end(), 0);
    std::fill(av2.begin(), av2.end(), 0);
    std::fill(sum_prog.begin(), sum_prog.end(), 0);
    std::fill(err_prog.begin(), err_prog.end(), 0);
    std::fill(sum_2_prog.begin(), sum_2_prog.end(), 0);

    //Evaluating call option with discrete step
    for(int i = 0; i < N; i++) {                            //Cycle on all blocks
        double sum1 = 0;
        for(int j = 0; j < L; j++) {                        //Cycle on the single block
            sum1 += Call(myRand, param, false);
        }
        //Normalising mean and variance
        ave[i] = sum1/L;
        av2[i] = (ave[i])*(ave[i]);
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){       
            sum_prog[i] += ave[j];
            sum_2_prog[i] += av2[j];
        }
        sum_prog[i]/=(i+1);
        sum_2_prog[i]/=(i+1);
        if(i == 0) {
            err_prog[i] = 0;
        } else err_prog[i] = sqrt((sum_2_prog[i] - sum_prog[i]*sum_prog[i])/double(i));
    }
    PrintResults(path[0], x, sum_prog, err_prog);

    std::fill(ave.begin(), ave.end(), 0);
    std::fill(av2.begin(), av2.end(), 0);
    std::fill(sum_prog.begin(), sum_prog.end(), 0);
    std::fill(err_prog.begin(), err_prog.end(), 0);
    std::fill(sum_2_prog.begin(), sum_2_prog.end(), 0);

    //Evaluating call option with discrete step
    for(int i = 0; i < N; i++) {                            //Cycle on all blocks
        double sum1 = 0;
        for(int j = 0; j < L; j++) {                        //Cycle on the single block
            sum1 += Pull(myRand, param, false);
        }
        //Normalising mean and variance
        ave[i] = sum1/L;
        av2[i] = (ave[i])*(ave[i]);
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i+1; j++){       
            sum_prog[i] += ave[j];
            sum_2_prog[i] += av2[j];
        }
        sum_prog[i]/=(i+1);
        sum_2_prog[i]/=(i+1);
        if(i == 0) {
            err_prog[i] = 0;
        } else err_prog[i] = sqrt((sum_2_prog[i] - sum_prog[i]*sum_prog[i])/double(i));
    }
    PrintResults(path[1], x, sum_prog, err_prog);

    cout << endl;
    return 0;
}



double max(double a, double b) {
    if (a > b) return a;
    else return b;
}

double Call(Random* myRand, parameter param, bool flag) {
    
    //Performs one step is flag is true
    if(flag == true) {
        double C = 0;
        for(int i = 0; i < 1; i++) {
            double Z = myRand -> Gauss(0, 1);
            double rad = param.sigma*Z*sqrt(param.T);
            double e = (param.r-param.sigma*param.sigma/2.)*param.T;
            double S_i = param.S_0 * exp(e + rad);
            C += exp(-param.r*param.T) * max(0, S_i - param.K);
        }
        return C;
    }
    else { //Proceeds in discreet steps if flag is false
        double C = 0;
        for(int i = 0; i < 1; i++) {
            double S_i = param.S_0;
            double told = 0;
            for(double t = 0; t < 1; t += 0.01) {
                double Z = myRand -> Gauss(0, 1);
                double rad = param.sigma*Z*sqrt((t-told));
                double e = (param.r-param.sigma*param.sigma/2.)*(t-told);
                S_i *= exp( e + rad );
                told = t;
            }
        C += exp(-param.r*param.T) * max(0, S_i - param.K);
        }
        return C;
    }   
}

double Pull(Random* myRand, parameter param, bool flag){

    if(flag == true) {
        double P = 0;
        for(int i = 0; i < 1; i++) {
            double Z = myRand -> Gauss(0, 1);
            double rad = param.sigma*Z*sqrt(param.T);
            double e = (param.r-param.sigma*param.sigma/2.)*param.T;
            double S_i = param.S_0 * exp(e + rad);
            P += exp(-param.r*param.T) * max(0, param.K - S_i);
        }
        return P;
    } else {
        double P = 0;
        for(int i = 0; i < 1; i++) {
            double S_i = param.S_0;
            double told = 0;
            for(double t = 0; t < 1; t += 0.01) {
                double Z = myRand -> Gauss(0, 1);
                double rad = param.sigma*Z*sqrt((t-told));
                double e = (param.r-param.sigma*param.sigma/2.)*(t-told);
                S_i *= exp( e + rad );
                told = t;
            }
        P += exp(-param.r*param.T) * max(0, param.K - S_i);
        }
        return P;        
    }
}