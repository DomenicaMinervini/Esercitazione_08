/*************************************************************
**************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
**************************************************************
*************************************************************/

#include <iostream>
#include <cmath> 

#ifndef __Funzioni__
#define __Funzioni__

double Psi(double, double, double); //Funzione d'onda.
//Rapporto p tra le probabilità nuova e vecchia.
double Rapporto_probabilita(double, double); 
double Minimum(double, double); //Minimo tra p e 1.
double T(double, double, double); //Contributo dell'energia cinetica.
double V(double); //Energia potenziale.
double L(double, double, double); //Funzione costo.
double Beta(double); //Inverso della temperatura.
//Funzione che campiona il peso di Boltzmann variando mu e sigma.
double Boltzmann(double, double, double); 

#endif //__Funzioni__

/*************************************************************
**************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
**************************************************************
*************************************************************/
