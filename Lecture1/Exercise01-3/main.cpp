#include "../../random/random.h"
#include "../../mylibraries/symbols.h"
using namespace std;

double GenerateAngle(Random* myRand);

int main(){

    intro();
    Random* myRand = new Random();
    myRand -> SetGenerator();

    double length = 1.;              //Lenght of the stick   
    double d = 5.;              //Distance between lines
    int M=1E5;                  //Total throws
    int N=100;                  //Number of blocks
    int L = int(M/N);           //Throws per block
    int hit = 0;
    int tot = 0;

    //Vectors for data blocking
    vector<double> x(N);                                //Block counter
    vector<double> ave(N,0);
    vector<double> av2(N,0);
    vector<double> sum_prog(N,0);
    vector<double> sum_2_prog(N,0);
    vector<double> err_prog(N,0);
    std::iota(x.begin(), x.end(), 0);                   //Fills with 0,1,2,3...

    const char* dirName = "results";
    CreateDirectory(dirName);
    vector<string> path = {"results/results_pi.txt", "results/results_angle.txt"};
    ofstream file1(path[0]);
    ofstream file2(path[1]);
    int counter = 0;

    for(int i = 0; i < N; i++) {                        //Cycle on all blocks
        double pi = 0;
        for(int j = 0; j < L; j++) {                    //Cycle on the single block

            //Repeat the experiment 100 times for a single π value
            for(int k = 0; k < 100; k++) {
                //Here is the effective Buffon Experiment
                //Generate first end of the stick in a quadrant [0,1]×[0,3]
                double y_start = myRand -> Rannyu(0, d);                //Only the vertical coordinate is needed
                double theta = GenerateAngle(myRand);
                if(counter < 500) file2 << theta << endl;               //the angle is printed in order to verify uniform sampling
                //generate second end of the stick
                double y_end = y_start + sin(theta);
                if(y_end < 0 || y_end > d) hit++;
                tot++;
                counter++;                                              //Only needed to print angles
            }
            pi += (2 * length * tot)/(hit * d);

        }
        //Normalising mean and variance
        ave[i] = pi/L;
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
    if(y > 0) {
        theta = acos(x/sqrt(x*x+y*y));
    } else if  (y < 0) {
        theta = M_PI + acos(x/sqrt(x*x+y*y));
    }
    return theta;
}