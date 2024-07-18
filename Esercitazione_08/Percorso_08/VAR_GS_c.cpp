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

   /* Test per trovare la combinazione delta_mu e delta_sigma tale che la % di accettazione del Metropolis sia del 50%.
   Le pongo = 0 una alla volta, e l'altra la passo da terminale. Trovo quando sia una sia l'altra danno circa il 50%.
   Con valori dell'ordine di 0.1-0.5 ho buoni valori dell'energia ma accettazione del 98%: troppo alta.
   Con la coppia 2.0 e 1.0 si ottiene accettazione del 55%, la migliore, ma le energie sono un po' sballate.
   Va bene la coppia delta_mu = 1.5/1.6 e delta_sigma = 0.75/0.8, con accettazione del 68% ed energie nella norma.
   Questi risultati si ottengono per temp = 1, blocchi = 20 e step = 50. Con step = 10 si ottiene accettazione del 50%.*/
   delta_mu = atof(argv[1]); 
   delta_sigma = atof(argv[2]);
   
   //Valori inizializzati.
   mu_old = 1.0;
   sigma_old = 0.5;
   x_old = 0;  
   totale = 0;

   Input(); //Genera i valori casuali.  
   Equilibrate(mu_old, sigma_old); //Equilibrazione prima della simulazione.
   
   fstream Nrg, Ems, Per; //Scrivo su tre file esterni: energia-incertezza, energia-mu-sigma e le percentuali di accetazione.
   Nrg.open("Energia_GS.dat",ios::app);
   Ems.open("Energia_mu_sigma.dat",ios::app);
   Per.open("Percentuali_acc.dat",ios::app);

   //Assegno i valori legati alla temperatura.
   delta_t = temp/annealing; 
   cout << "Temperatura, cicli, delta t: " << temp << "  " << annealing << "  " << delta_t << endl;
   
   for(item=0; item<annealing; item++){ //Cicli della temperatura che decresce.
      accettati = 0; //Azzero la % di accettazione all inizio del ciclo in cui cerco mu e sigma.
	   
      for(iste=0; iste<step; iste++){ //Cicli di Metropolis con cui vario mu e sigma.	  
         //Con i valori iniziali di mu e sigma muovo la posizione e calcolo il valore dell'energia corrispondente: en_old.  
         en_old = Energy(mu_old, sigma_old); //Muovo la posizione con i vecchi valori di mu e sigma: trovo x_old.
         Nrg << iste << "  " << mean_mean << "  " << incertezza << "  " << mu_old << "  " << sigma_old << endl; 

         //Vario i parametri mu e sigma: diventano mu_new e sigma_new. Cerco i parametri nuovi
         //in due intervalli centrati nei valori vecchi e che si estendono di [-0.5,+0.5].
         mu_new = mu_old + delta_mu*(rnd.Rannyu()-0.5);
         sigma_new = sigma_old + delta_sigma*(rnd.Rannyu()-0.5);

         //Con i valori nuovi di mu e sigma muovo di nuovo la posizione e calcolo il nuovo valore dell'energia: en_new.
         en_new = Energy(mu_new, sigma_new); //Muovo la posizione con i nuovi valori di mu e sigma: trovo x_old.
         Nrg << iste << "  " << mean_mean << "  " << incertezza << "  " << mu_new << "  " << sigma_new << endl; 

         //ALGORITMO DI METROPOLIS: confronto le energie ottenute con i parametri variazionali nuovi e vecchi.
         /*1) Rapporto tra le probabilità degli stati con energia nuova e con energia vecchia.
           Questo valore ci rivela le proprietà fisiche del sistema. La probabilità dipende dalla temperatura 
           e dall'energia in funzione della posizione.*/   
         pro = Boltzmann(temp, en_old, en_new); 
   
         /*2) Scelgo un solo valore tra 1 e il rapporto delle probabilita': il minimo. 
           Se art = 1 è un fatto gradito, perché vuol dire che la mossa sarà sicuramente accettata.*/
         art = Minimum(1.0, pro);
   
        /*3) Estraggo un valore in una distribuzione uniforme tra 0 e 1 e lo confronto con min, che sarà o 1 o pro.
          Se art e' maggiore del valore estratto, accetto la mossa nuova, altrimenti tengo i valori di mu e sigma che ho.
          Se art = 1 la mossa sarà accettata sicuramente, perché questo valore determinante è compreso in [0,1).
          L'elemento val completa l'algoritmo ed è l'ostacolo che la mossa proposta deve superare per essere accettata.*/	 	 
         val = rnd.Rannyu();

         if(val <= art){ 
	        mu_old = mu_new; //Accetto la mossa. I nuovi parametri variazionali diventano quelli attuali.
	        sigma_old = sigma_new; //Questi sono i valori dello stato scelto e ne calcolo l'energia.
            accettati++; //Numero di accettati per ogni blocco. Azzerato ad ogni ciclo di blocco.
	        totale++; //Numero di tentativi accettati di Metropolis dei parametri in tutta la simulazione. Non si azzera mai.
         }	
  	  
	 //en_old = Energy(mu_old, sigma_old); //Calcolo l'energia con i valori di mu e sigma selezionati.
	 //Ems << temp << "  " << iste << "  " << en_old << "  " << mu_old << "  " << sigma_old << endl;
      }

      percentuale = 100.0*accettati/step; //Il divisore è il numero totale di lanci nel singolo ciclo mu-sigma.
      Per << temp << "  " << percentuale << "%" << endl;  

      temp -= delta_t; //Abbasso la temperatura per il prossimo ciclo di un valore delta_t.	*=  
   }

   per_tot = (100.0*totale)/(annealing*step); //Il divisore è il numero totale di lanci di tutta la simulazione.
   cout << "La percentuale di accettazione totale della simulazione e':" << per_tot << "%." << endl;   
	  
   cout << "E' stato prodotto il file: 'Energia_GS.dat'." << endl; 
   cout << "E' stato prodotto il file: 'Energia_mu_sigma.dat'." << endl; 
   cout << "E' stato prodotto il file: 'Percentuali_acc.dat'." << endl; 

   Nrg.close();
   Ems.close();
   Per.close();
   
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

void Reset_block(){	
   //Azzera i valori all'inizio di ogni ciclo della media a blocchi.
   r = 0; 
   sum = 0;  
}

void Reset_cycle(){	
   //Inizializza questi valori all'inizio di ogni ciclo della media a blocchi
   sum_mean = 0;  
   sum_square = 0;     
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
      //è il numero zero, N è già ridotto di 1. Quindi nella formula divido per la radice di N, ossia per iblk blocchi.
   }
}

//FUNZIONE MAIN DELL'ESERCIZIO 08.1.
double Energy(double mu, double sigma){
   Reset_cycle(); //Azzera i valori somma usati nella media a blocchi. Faccio più medie a blocchi nel corso della simulazione.
   
   //Scrivo su dei file esterni i valori della media progressiva dell'energia con l'incertezza e di %accettazione. 
   //fstream Nrg; 
   //Nrg.open("Energia_GS.dat",ios::app);

   //Scrivo su un file esterno i valori della posizione x_old. Li ottengo campionando la funzione d'onda al quadrato.
   //Con la funzione Move() sposto la posizione nello spazio. Ripeto per un certo numero di volte.

   for(iblk=0; iblk<blocchi; iblk++){ //Per ogni blocco,
      Reset_block(); //azzera i valori del blocco,
	  
      for(int j=0; j<lanci; j++){ //all'interno di ciascun blocco.
     	  Move(mu, sigma); //Seleziona i valori di x,
 	  Accumulate(mu, sigma); //calcola i potenziali per ogni mossa accettata e li somma.
      }

      Average(iblk); //Calcola la media,
	  
	  if(iblk == blocchi-1){ //Quando hai raggiunto l'ultimo elemento del blocco, ossia la media più accurata,
	     Error(iblk); //calcola l'incertezza associata a questa media, solo per l'ultimo blocco.
         //Nrg << mean_mean << "  " << incertezza << endl; //scrivi i valori su un file esterno.
      }	 
   }
   
   //Nrg.close();
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
