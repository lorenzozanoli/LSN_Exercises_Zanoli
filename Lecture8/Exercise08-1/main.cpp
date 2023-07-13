#include "main.h"
using namespace std;
//==============================================

int main (){

    Input();

    if(!SA){
        MeasureAverageH();
    }

    if(SA){
        WriteTraj.open("output/traj.out");

        H = 0; // riazzero tutto prima di iniziare SA
        attempted = 0;
        accepted = 0;

        while (beta < 500){            
            for(int i = 0; i < step_in_beta; i++){

                oldH = mean_prog_H/nblk;
                // equivalente di passo tentativo: esploro mu e sgm in funzione di beta
                //più vado avanti meno esploro, mi avvicino all'equilibrio quindi 
                d_mu = rnd.Rannyu(-Lmu,Lmu)/beta; 
                d_sgm = rnd.Rannyu(-Lsgm,Lsgm)/beta;

                mu += d_mu;
                sgm += d_sgm;

                MeasureAverageH();
                attempted++;
                
                // calcolo probabilità
                newH = mean_prog_H/nblk;
                errH = error(mean_prog_H/nblk,var_prog_H/nblk,(nblk-1));
                deltaH = oldH - newH;

                q = exp(beta*(deltaH)); // con +, perchè devo andare verso il minimo di H
                A = min(1,q);

                // giudico con metropolis, confrontantdo newH con oldH
                if (A==1){
                    oldH = newH;
                    accepted++;
                }
                else{
                    if(rnd.Rannyu() < A){
                        oldH = newH;
                        accepted++;
                    } else {mu -= d_mu; sgm -= d_sgm;} // se non accetto la mossa, ripristino
                }
                step++;

            }

            cout << "step:    " << step << endl 
                 << "beta:    " << beta << endl 
                 << "H:       " << oldH << endl
                 << "errH:    " << errH << endl
                 << "mu:      " << mu << endl
                 << "sgm:     " << sgm << endl
                 << endl; 
            WriteTraj << beta << " " << mu << " " << sgm << " " << oldH  << " " << errH << endl;

            // abbasso temperatura (beta = 1/T)
            beta += db;
            
        }
        cout << endl
            << ",=======================================" << endl
            << "| Minimo di H: " << CYAN << oldH << " +/- " << errH << RESET << endl
            << "| mu:          " << CYAN << mu << RESET << endl
            << "| sigma:       " << CYAN << sgm << RESET << endl
            << "'=======================================" << endl << endl;

        //Ora voglio calcolare l'errore

        for(int i = 0; i < 10000; i++) {
            oldH = mean_prog_H/nblk;
            //Ora beta è fissato, voglio solo plottarli per calcolare l'errore
            d_mu = rnd.Rannyu(-Lmu,Lmu)/beta; 
            d_sgm = rnd.Rannyu(-Lsgm,Lsgm)/beta;

            mu += d_mu;
            sgm += d_sgm;

            MeasureAverageH();
            attempted++;
                
            // calcolo probabilità
            newH = mean_prog_H/nblk;
            errH = error(mean_prog_H/nblk,var_prog_H/nblk,(nblk-1));
            deltaH = oldH - newH;

            q = exp(beta*(deltaH)); // con +, perchè devo andare verso il minimo di H
            A = min(1,q);

            // giudico con metropolis, confrontantdo newH con oldH
            if (A==1){
                oldH = newH;
                accepted++;
            }
            else{
                if(rnd.Rannyu() < A){
                    oldH = newH;
                    accepted++;
                } else {mu -= d_mu; sgm -= d_sgm;} // se non accetto la mossa, ripristino
            }
            HistoParameters << mu << "\t" << sgm << endl;
            step++;
        }
        HistoParameters.close();

    }
    WriteTraj.close();



    return 0;
}

// ================================================================================
