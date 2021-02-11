#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme() :  _t(0.)
{
    cout<<"TSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"<<endl;
    Eigen:: VectorXd rho(100);
    _rho=rho;


}


void TimeScheme::initialize()
{
    int n=10;
    cout<<"ok"<<endl;
    Eigen:: VectorXd rho(100);
    cout<<rho<<endl;
    _rho=rho;
    cout<<"TSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"<<endl;
    _T.resize(n*n);
    _xi.resize(n*n);
    cout<<"TSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"<<endl;
    
    for (int i = 0; i < pow(n,2); i++)
    {
        cout<<_rho.size()<<endl;
        _rho(i)=1500;
        cout<<"TSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"<<endl;

        _T(i)=100;

        cout<<"TSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"<<endl;

        _xi(i)=0;
    }  
        cout<<"size de T  : \n "<< _T.size()<<endl;
}



// Euler Explicite
void TimeScheme::Advance_rho()
{
    double rhov=1500.;
    double rhop=1000.;
    double rcpv=15000.;
    double rcpp=10000.;
    double dt  =10E-5;
    double Ta  =6000. ;
    _temp=_rho;
    _temp1=_rcp;
    //remlissage de rho 
    for (int i=0; i<_rho.size(); i++)
    {
        _rho(i)=_rho(i)-rhov*dt*1000.*((_rho(i)-rhop)/(rhov-rhop))*exp(-Ta/_T(i));
        _xi(i) =(rhov-_rho(i))/(rhov-rhop);
        _rcp(i)=(1.-_xi(i))*rcpv+_xi(i)*rcpp;
         
    }
}

void TimeScheme::Advance_T()
{
    double rhov=1500.;
    double rhop=1000.;
    double dt  =10E-2;
    double Ta  =6000.;
    double Lm  =130.*10E3;
    _temp=_rho;
    _fin_vol->Build_flux_mat_and_rhs(_t, _xi);
    SparseMatrix<double> A(_fin_vol->Get_flux_matrix());
    VectorXd b(_fin_vol->Get_BC_RHS());
    VectorXd Flux;
    Flux=A*_T+b;
    Advance_rho();
    for (int i=0; i<_xi.size(); i++)
    {
       _T(i)=200.-Lm*_rho(i)/_rcp(i)+Flux(i)*dt/_rcp(i)+_temp1(i)*(_T(i)-200.)/_rcp(i)+Lm*_temp(i)/_rcp(i);
    }

    _t=_t+dt;
    
}








#define _TIME_SCHEME_CPP
#endif