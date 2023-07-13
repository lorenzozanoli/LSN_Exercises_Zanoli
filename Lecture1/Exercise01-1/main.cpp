#include "../../random/random.h"
#include "../../mylibraries/symbols.h"
using namespace std;

int main(){

    intro();
    Random* myRand = new Random();
    myRand -> SetGenerator();

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
    string path_1 = "results/results_r.txt";
    string path_2 = "results/results_r2.txt";
    string path_3 ="results/results_chi2.txt";

    //Evaluating integral of r in [0,1]
    for(int i = 0; i < N; i++) {                        //Cycle on all blocks
        double sum1 = 0;
        for(int j = 0; j < L; j++) {                    //Cycle on the single block
            sum1 += myRand -> Rannyu();                 //THIS IS THE ONLY THING TO CHANGE
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
    PrintResults(path_1, x, sum_prog, err_prog);

    //Evaluating integral of r² in [0,1]
    //Empyting all vectors;
    std::fill(ave.begin(), ave.end(), 0);
    std::fill(av2.begin(), av2.end(), 0);
    std::fill(sum_prog.begin(), sum_prog.end(), 0);
    std::fill(err_prog.begin(), err_prog.end(), 0);
    std::fill(sum_2_prog.begin(), sum_2_prog.end(), 0);

    for(int i = 0; i < N; i++) {                            //Cycle on all blocks
        double sum1 = 0;
        for(int j = 0; j < L; j++) {                        //Cycle on the single block
            sum1 += pow(myRand -> Rannyu()-0.5, 2);         //THIS IS THE ONLY THING TO CHANGE
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
    PrintResults(path_2, x, sum_prog, err_prog);

    //Chi2 Evaluation of Uniform extraction
    double r = 0;
    int intervals = 100;                                 //Number of intervals in [0,1]
    int n = 10000;                                         //Number of extractions per experiment
    int exp = 10000;                                     //Number of experiments
    vector<int> counter(intervals, 0);                   //Counts elements in each interval
    double chi2 = 0;

    ofstream file3(path_3);
        if (!file3.is_open()) {
        cerr << "Error creating results_chi2.txt" << endl;
        return -1;
    }

    for(int i = 0; i < exp; i++) {                      //Repetition of the experiments
        chi2 = 0;
        std::fill(counter.begin(), counter.end(), 0);
        for(int j = 0; j < n; j++) {                    //Extracting random numbers and counting appearance in each bin
            r = myRand -> Rannyu();
            counter[(int)floor(r * intervals)]++;
            
        }
        for(int k = 0; k < intervals; k++) {            //Evaluating chi2 per each bin and summing
            double expected = n/intervals;
            double observed = counter[k];
            chi2 += pow(observed - expected, 2)/expected;
        }
        file3 << chi2 << endl;
    }
    cout << "File with " << LCHI SUP2 << " results has been printed" << endl;
    
	file3.close();
    cout << endl;
    return 0;
}