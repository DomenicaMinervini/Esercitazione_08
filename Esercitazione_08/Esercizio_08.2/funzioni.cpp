/*************************************************************
**************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
**************************************************************
*************************************************************/

#include "funzioni.h"

using namespace std;

double Psi(double x, double mu, double sigma){ 
   double sigma_2=pow(sigma,2);
   double xm=pow((x-mu),2);
   double xp=pow((x+mu),2);
   return exp(-0.5*xm/sigma_2)+exp(-0.5*xp/sigma_2);
}

//Rapporto p tra le probabilità nuova e vecchia con le funzioni d'onda.
double Rapporto_probabilita(double psi_new, double psi_old){ 
   return (psi_new*psi_new)/(psi_old*psi_old);
}

double Minimum(double one, double two){ //Minimo tra p e 1.
   if(one < two) return one;
   else return two;
}

double T(double x, double mu, double sigma){ //Contributo dell'energia cinetica.
   double sigma_2=pow(sigma,2);
   double sigma_4=pow(sigma,4);
   double xm=pow((x-mu),2);
   double xp=pow((x+mu),2);
   return (-0.5/Psi(x, mu, sigma))*((-1.0/sigma_2+xm/sigma_4)*exp(-0.5*xm/sigma_2)+(-1.0/sigma_2+xp/sigma_4)*exp(-0.5*xp/sigma_2));  
}

double V(double x){ //Energia potenziale.
   return pow(x,4)-(5.0/2.0)*pow(x,2);
}

double L(double x, double mu, double sigma){ //Funzione costo L, l'energia totale del sistema, che deve diminuire.
   return T(x, mu, sigma)+V(x);
}

double Beta(double temp){ //Inverso della temperatura.
   return 1.0/temp;
}

//Calcola il rapporto delle probabilità con il peso di Boltzmann.
double Boltzmann(double temp, double en_old, double en_new){
   //return exp(-Beta*(en_new-en_old);
   return exp((1.0/temp)*(en_old-en_new));

}


/*************************************************************
**************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
**************************************************************
*************************************************************/
