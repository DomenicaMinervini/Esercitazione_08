/*************************************************************
**************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
**************************************************************
*************************************************************/

#include "VAR_GS.h"

using namespace std;

int main (int argc, char *argv[]){

   //Valori inizializzati.
   mu_old = atof(argv[1]); //Passo mu e sigma da terminale.
   sigma_old = atof(argv[2]);
   accettati=0;
   x_old=0;  

   int i; //Contatori.

   Input(); //Genera i valori casuali.  
   Equilibrate(mu_old, sigma_old); //Equilibrazione prima della simulazione.

   accettati=0; //Azzero ancora prima del ciclo: e' aumentato con Equilibrate().

   for(i=0; i<1; i++){
      en_old = Energy(mu_old, sigma_old); //Muovo la posizione con i vecchi valori di mu e sigma: trovo x_old.
      cout << en_old << endl;

      //Vario i parametri mu e sigma: diventano mu_new e sigma_new. Cerco i parametri nuovi
      //in due intervalli centrati nei valori vecchi e che si estendono di [-0.5,+0.5].
      //mu_new    =    mu_old+rnd.Rannyu()-0.5;
      //sigma_new = sigma_old+rnd.Rannyu()-0.5;

      //Con i valori nuovi di mu e sigma muovo di nuovo la posizione e calcolo 
   }
 
   rnd.SaveSeed();  
   return 0;
}


void Input(){

   //GENERATORE DI NUMERI CASUALI
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;	
	
}

void Reset(){	
   //Azzera i valori all'inizio di ogni ciclo.
   accettati=0;
   percentuale=0;
   r=0; 
   sum=0;  
}

void Move(double mu, double sigma){ //mu e sigma assegnati ad ogni uso di Move.

   //GENERATORE DELLE COORDINATE E FUNZIONE D'ONDA VECCHIA E NUOVA  
   //CALCOLO DELLA P(X): i valori selezionati con l'algoritmo di Metropolis fanno la P(x).
   //Li seleziona in base al rapporto delle probablita' nuova e vecchia della funzione d'onda.
    x_new = x_old+delta*(rnd.Rannyu()-0.5);
    psi_old = Psi(x_old, mu, sigma);
    psi_new = Psi(x_new, mu, sigma);

    //ALGORITMO DI METROPOLIS
    //Rapporto tra le funzioni d'onda nuova e vecchia al quadrato, ossia tra le probabilita'.
    p = Rapporto_probabilita(psi_new,psi_old);
	
    //Scelgo il minimo tra 1 e il rapporto delle probabilita'.
    alpha = Minimum(1.0, p);
    //Estraggo un valore in una distribuzione uniforme tra 0 e 1:
    //se alpha e' maggiore del valore estratto, accetto la mossa nuova,
    //altrimenti tengo la posizione in cui mi trovo.
    value = rnd.Rannyu();

    if(value<=alpha){ //Accetto la mossa: la posizione nuova diventa quella attuale.
       x_old = x_new;
       accettati++;
    }
	
}

void Equilibrate(double mu, double sigma){
   //Equilibrazione prima della simulazione. Il primo blocco serve per equilibrare il sistema.
   for(int i=0; i<lanci; i++){ //Questo blocco non interverrà nella media a blocchi.
      Move(mu, sigma);   
   }
}

void Accumulate(double mu, double sigma){	
    //CALCOLO DELLA G(X): le x_old selezionate prima servono a calcolare la funzione da integrare.
    r = L(x_old, mu, sigma);
    sum += r;	
}

void Average(int iblk){
   //Alloco i vettori delle medie di blocco e dei quadrati delle medie.
   mean[iblk] = sum/double(lanci); //Media dei valori generati in un blocco: la metto in un vettore, per ogni blocco.
   square[iblk] = (mean[iblk])*(mean[iblk]); //Quadrato della media dei valori: la metto in un vettore, per ogni blocco.

   sum_mean += mean[iblk]; //Somma degli elementi del vettore medie, dal valore 0 al valore numero di blocchi.
   sum_square += square[iblk]; //Somma degli elementi del vettore dei quadrati, dal valore 0 al valore numero di blocchi.

   mean_mean = sum_mean/double(iblk+1); //Calcolo la media dei valori medi di ogni blocco, aggiungendo ogni volta un blocco.
   mean_square = sum_square/double(iblk+1); //Faccio lo stesso per i quadrati delle medie.
}   

void Error(int iblk){
   if(iblk == 0){ //Nel caso in cui si abbia un solo blocco, l'incertezza statistica è posta a 0.
      varianza = 0;
      stdev = 0;
      incertezza = 0; 
   }
   else{
      varianza = mean_square - mean_mean*mean_mean; //Varianza per la media di 1 blocco, 2 blocchi, ... , N blocchi,
      stdev = sqrt(varianza); //lo stesso per la deviazione standard.
      incertezza = stdev/(sqrt(iblk)); //Incertezza stastistica del valor medio, se sono stati generati iblk blocchi.
      //Nella formula dell'incertezza statistica si divide per un fattore radice di N-1. Dato che il primo elemento del vettore
      //è il numero zero, N è già ridotto di 1. Quindi nella formula divido per la radice di N, ossia per iblK blocchi.
   }
}

//FUNZIONE MAIN DELL'ESERCIZIO 08.1.
double Energy(double mu, double sigma){

   //Scrivo su dei file esterni i valori della media progressiva dell'energia con l'incertezza e di %accettazione. 
   /*fstream Pos, Per; 
   Pos.open("Energia_GS.dat",ios::out);
   Per.open("percentuale.dat",ios::out);*/

   //Scrivo su un file esterno i valori della posizione x_old.  
   //Li ottengo campionando la funzione d'onda al quadrato.
   //Con la funzione Move() sposto la posizione nello spazio.
   //Ripeto per un certo numero di volte.
   //Pos << "Numero di blocchi, media progressiva per i blocchi e incertezza." << endl;

   for(iblk=0; iblk<blocchi; iblk++){ //Per ogni blocco,
      Reset(); //azzera i valori del blocco,

      for(int j=0; j<lanci; j++){ //all'interno di ciascun blocco.
     	  Move(mu, sigma); //Seleziona i valori di x,
 	  Accumulate(mu, sigma); //calcola i potenziali per ogni mossa accettata e li somma.
      }

      Average(iblk); //Calcola la media,
      Error(iblk); //calcola l'incertezza associata a questa media.
	  
      //step = (iblk+1)*lanci; //Numero di lanci a questo punto = numero blocchi * elementi in un blocco.
      //Pos << iblk+1 << "  " << mean_mean << "  " << incertezza << endl; //Scrivo i valori su un file esterno.
	  
   }

   //Calcolo la percentuale di accettazione e la scrivo su un file esterno.
   /*percentuale = 100.0*accettati/double(lanci);
   Per << "Percentuale di accettazione del blocco delle posizioni." << endl;
   Per << percentuale << "%" << endl;*/

   //Chiudo i file. 
   /*cout << "E' stato prodotto il file: 'posizione_energia.dat'." << endl; 
   cout << "E' stato prodotto il file: 'percentuale.dat'." << endl; 
   Pos.close();
   Per.close();*/

   return mean_mean;
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
