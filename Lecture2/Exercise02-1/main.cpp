#include "../../random/random.h"
#include "../../mylibraries/symbols.h"
using namespace std;

int main(){

    intro();
    Random* myRand = new Random();
    myRand -> SetGenerator();

    int M = 1E5;                //Total throws in each block
    int N = 100;                //Number of blocks

    auto f = [](double x) -> double{return (M_PI/2.)*cos(M_PI*x/2.);};
    //Sampling with y = 2(1-x), the cumulative is the following
    auto sampling_1 = [](double y) -> double{return 1 - sqrt(1-y);};
    //Sampling with y = π/2(1-x²) by accept reject method
    auto sampling_2 = [](double x) -> double{return (M_PI/2.)*(1-x*x);};

    const char* dirName = "results";
    CreateDirectory(dirName);
    vector<string> path = {"results/results_uniform.txt", "results/results_linear.txt", "results/results_parabolic.txt", "results/results_antithetic.txt"};

    //Will not be using vectors in order to create less variables
    double sum_prog_uniform = 0;
    double sum_prog_linear = 0;
    double sum_prog_parabolic = 0;
    double sum_prog_antithetic = 0;
    double sum_prog_2_uniform = 0;
    double sum_prog_2_linear = 0;
    double sum_prog_2_parabolic = 0;
    double sum_prog_2_antithetic= 0;
    double err_prog_uniform = 0;
    double err_prog_linear = 0;
    double err_prog_parabolic = 0;
    double err_prog_antithetic = 0;

    ofstream file_1(path[0]);
    ofstream file_2(path[1]);
    ofstream file_3(path[2]);
    ofstream file_4(path[3]);

    for(int i = 0; i < N; i++) {
        double mean_uniform = 0;
        double mean_linear = 0;
        double mean_parabolic = 0;
        double mean_antithetic = 0;
        for(int j = 0; j < M; j++) {
            //Integral with uniform sampling
            double x = myRand -> Rannyu();
            double y = f(x);
            mean_uniform += y;
            //Integral with linear importance sampling (cumulative)
            x = myRand -> SamplingCumulative(sampling_1, 0, 1);
            y = f(x) / (2-2*x);
            mean_linear += y;
            //Integral with parabolic importance sampling (accept-reject)
            x = myRand -> SamplingAR(sampling_2, 0, 1, 3./2.);
            y = f(x) / (1.5*(1-x*x));
            mean_parabolic += y;
            //Antithetic variates      ------> È QUESTO QUI CHE DEVI GUARDARE + LEZ 2 SLIDE 44
            x = myRand -> Rannyu();
            y = 0.5*(f(x) - f(1-x));
            mean_antithetic += y;
            
        }
        //Printing results for uniform sampling
        mean_uniform /= M;
        sum_prog_uniform += mean_uniform;
        sum_prog_2_uniform += mean_uniform*mean_uniform;
        err_prog_uniform = sqrt((sum_prog_2_uniform/(i+1) - pow(sum_prog_uniform/(i+1),2))/i);
        file_1 << i << "\t" << sum_prog_uniform/(i+1) << "\t";
        if(i==0) file_1 << 0 << endl;
        else file_1 << err_prog_uniform << endl;

        //Printing results for linear sampling
        mean_linear /= M;
        sum_prog_linear += mean_linear;
        sum_prog_2_linear += mean_linear*mean_linear;
        err_prog_linear = sqrt((sum_prog_2_linear/(i+1) - pow(sum_prog_linear/(i+1),2))/i);
        file_2 << i+1 << "\t" << sum_prog_linear/(i+1) << "\t";
        if(i==0) file_2 << 0 << endl;
        else file_2 << err_prog_linear << endl;

        //Printing results for parabolic sampling
        mean_parabolic /= M;
        sum_prog_parabolic += mean_parabolic;
        sum_prog_2_parabolic += mean_parabolic*mean_parabolic;
        err_prog_parabolic = sqrt((sum_prog_2_parabolic/(i+1) - pow(sum_prog_parabolic/(i+1),2))/i);
        file_3 << i+1 << "\t" << sum_prog_linear/(i+1) << "\t";
        if(i==0) file_3 << 0 << endl;
        else file_3 << err_prog_parabolic << endl;

        //Printing results for antithetic variates
        mean_antithetic /= M;
        sum_prog_antithetic += mean_antithetic;
        sum_prog_2_antithetic += mean_antithetic*mean_antithetic;
        err_prog_antithetic = sqrt((sum_prog_2_antithetic/(i+1) - pow(sum_prog_antithetic/(i+1),2))/i);
        file_4 << i+1 << "\t" << sum_prog_linear/(i+1) << "\t";
        if(i==0) file_4 << 0 << endl;
        else file_4 << err_prog_antithetic << endl;

        
    }

    cout << endl;
    return 0;
}
