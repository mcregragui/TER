#ifndef _TIME_SCHEME_H
#include "Eigen/Dense"
#include "FiniteVolume.h"






class TimeScheme
{
    protected:
        //les deux vecteirs solution
        Eigen:: VectorXd _rho;
        Eigen:: VectorXd _T  ;
        Eigen:: VectorXd _xi;
        Eigen:: VectorXd _rcp;
        Eigen:: VectorXd _temp; //sert à stocké _rho(n-1) 
        Eigen:: VectorXd _temp1; //sert à stocké _rho(n-1)
        //Time
        double _t; 
        // Pointeur vers la classe FiniteVolume
        FiniteVolume* _fin_vol;

    public:
        // Constructeur par défaut
        TimeScheme();
        // Une etape du shémat de temps
        void Advance_rho();
        void Advance_T();
        void initialize();
        Eigen:: VectorXd Get_rho(){return _rho;};
        Eigen:: VectorXd Get_T(){return _T;};
};


#define _TIME_SCHEME_H
#endif