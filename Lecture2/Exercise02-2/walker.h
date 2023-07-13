#ifndef __FunzioneBase_h__
#define __FunzioneBase_h__
#include "../../mylibraries/symbols.h"

class Walker {

    public:
    Walker() {m_x = 0, m_y = 0, m_z = 0;};
    Walker(double x, double y, double z) {m_x = x, m_y = y, m_z = z;};

    double GetX() const {return m_x;};
    double GetY() const {return m_y;};
    double GetZ() const {return m_z;};
    double GetR() const {return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);};
    void Print() { std::cout << "[" << m_x << ", " << m_y << ", " << m_z << "]" << endl; }; 
    void PrintFile(string path) {
        ofstream file(path, ios::app);
        file << m_x << "\t" << m_y << "\t" << m_z << endl;
        file.close();
    }

    Walker Reset() {
        m_x = 0;
        m_y = 0;
        m_z = 0;
        return *this;
    }

    Walker DiscreteStep(Random* myRand, double a) {

        // Choosing direction: 0 = x, 1 = y, 2 = z
        int d = int(myRand->Rannyu(0, 3));
        // Choosing verse: 0 = forward, 1 = backwards
        int sign = int(myRand->Rannyu(0, 2));

        double step = pow(-1, sign)  * a;
        if (d == 0) m_x += step;
        else if (d == 1) m_y += step;
        else if (d == 2) m_z += step;
            
        return *this;
 
    }

    Walker ContinueStep(Random* myRand, double a) {
        // Choosing direction: 0 = θ, 1 = φ
        double* angle = GenerateAngle3D(myRand);

        m_x += a * cos(angle[0]) * sin(angle[1]);
        m_y += a * sin(angle[0]) * sin(angle[1]);
        m_z += a * cos(angle[0]);

        return *this;
    } 

    

    private:
    double m_x, m_y, m_z;

};

#endif