/*************************************************************
**************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
**************************************************************
*************************************************************/

#include <string>
#include "random.h"
#include "funzioni.h"

#ifndef __VAR_GS__
#define __VAR_GS__

//VARIABILI

//Numeri random.
int seed[4];
Random rnd;

//Parametri variazionali. 
double mu_old, mu_new;
double sigma_old, sigma_new; 
double delta_mu = 0.2; //Ampiezze degli intervalli in cui pesco mu e sigma.
double delta_sigma = 0.1; //Va bene anche la coppia 1.5 e 0.75, con accettazione del 68%.

//Definizione dei blocchi.
const int blocchi = 20; //Numero di blocchi, 20.
const int lanci = 1000; //Numero di elementi per ogni blocco, 1000.
int step = 10; //Numero di passi, cicli di Metropolis per mu e sigma per una temperatura.
int iblk; //Indice di blocco per la media a blocchi (Esercizio 08.1).
int iste; //Indice di Metropolis per mu e sigma.

//Percentuale di accettazione.
int accettati, totale; //Numero di accettati per Metropolis mu e sigma e su tutta la simulazione.
double percentuale, per_tot; //Relative percentuali di accettazione.

//Variabili del problema legate alla temperatura.
double temp = 1.0; //Temperatura iniziale che decresce.
int annealing = 300; //Numero di cicli sulla temperatura. 
int item; //Indice che conta il decrescere della temperatura.
double delta_t; //Passo con cui decresce la temperatura.
double beta; //Inverso della temperatura.

//Variabili del problema legate alla posizione.
double x_old, x_new; //posizione
double psi_old, psi_new; //funzione d'onda

//Variabili del problema legate all'energia.
double en_old, en_new; //Energie vecchia e nuova legate a x_old e a x_new.

//Algoritmo di Metropolis.
double delta = 2.43; //passo, scelto in modo che % accettazione sia del 50%. 
//percentuali ottenute nel test: 50.01%, 50.30%, 50.75%, 50.60%, 50.18%, per mu = 1.0 e sigma = 0.5. 
double p, alpha, value; //rapporto probabilita', probabilita' scelta, valore tra 0 e 1 per spostamenti x.
double pro, art, val; //rapporto probabilita', probabilita' scelta, valore tra 0 e 1 per mu e sigma.

//Media a blocchi.
double r, sum;
double mean[blocchi], square[blocchi];
double sum_mean, sum_square;
double mean_mean, mean_square;

//Varianza, deviazione standard, incertezza statistica.
double varianza, stdev, incertezza;
double incer_old, incer_new;

//FUNZIONI

//Funzione che da il set al seme del generatore di numeri casuali.
void Input();
//Funzione che azzera i valori all'inizio di ogni ciclo della singola media a blocchi.
void Reset_block();
//Funzione che campiona la funzione d'onda, calcola il rapporto tra le probabilità 
//delle funzioni d'onda e sceglie la x più probabile per effettuare lo spostamento.
void Move(double, double);
//Funzione che muove i parametri mu e sigma con un altro algoritmo di Metropolis.
void Move_parameters(double);
//Funzione che fa l'equilibrazione del sistema prima della vera simulazione.
void Equilibrate(double, double);
//Funzione che calcola la funzione da integrare.
void Accumulate(double, double);
//Funzione che alloca i vettori delle medie di blocco e dei quadrati delle medie.
void Average(int);
//Funzione che calcola l'incertezza statistica: la barra di errore.
void Error(int);
//Funzione che azzera i valori somma prima di ogni media a blocchi.
void Reset_cycle();
//Funzione che calcola l'energia con media a blocchi progressiva.
//E' la funzione main dell'esercizio 08.1. 
//La eseguo ogni volta che calcolo en_old.
double Energy(double, double);

#endif //__VAR_GS__

/*************************************************************
**************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
**************************************************************
*************************************************************/
