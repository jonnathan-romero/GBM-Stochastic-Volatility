// main.cpp. Jonnathan Romero, March 2017
// compile with...
// c++ -o hw3 homework3.cpp
// execute with...
// ./hw3

#include "Declarations.h"
#include <iostream>
#include <cmath>

const double s0    = 100;          // starting stock price
const double rr    = 0.05;         // risk free rate
const double T     = 0.5;          // time to espiration
const double V2    = 0.30 * 0.30;  // long term voltility
const double V2_0  = 0.35 * 0.35;  // volatility at time 0

const double alpha = 0.080;        // weight for long terms volatility
const double beta  = 0.885;        // weight for previous volatility
const double gmma  = 0.035;        // weight for current return * return

const int steps    = 126;          // steps for gbm (252 days in a year)

double BlackScholes(const double&, double, const double&, const double&, const double&);  // black scholes equation
double ImpliedVol (const double&, const double&, const double&, const double&, const double&); // implied volatility

struct Op{
	double K;                     // strike for this option
	double Pr;                    // price for this option under BS
	double EX;
	double EX2;
	double STD;
	double error;
	double impVol;                // implied vol
};

int main(){

	int seed = 1;                 // chose a seed
	int i;                        // looper, great movie btw
	int n = 0;

	MTUniform (seed);             // seed the RNG 

	double st;                    // current stock price
	double _st;                   // antithetic stock price
	double dt;                    // change in time for each step
	double mu;                    // drift parameter for gbm
	double  B;                    // brownian motion scaled to fit vol

	dt = T / steps;               // change in time is total time divided by steps

	Op calls[11];
    double sqdt   = sqrt(dt);     // sqrt of dt
    double U = 0;                 // for standard norm

    bool done = 0;
    double epsilon = 0.01;


	for(i=0; i<11; ++i){
		calls[i].K = 50 + (i*10); // strikes K = 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150
		calls[i].EX = 0;
		calls[i].EX2 = 0;
	}

	std::cout<<"--Errors-- (+-.01 to finish)"<<std::endl;

    while(!done){

		double sigma2 = V2_0;         // volatility used for gbm and option pricing
		double _sigma2 = V2_0;        // antithetic volatility
		st = s0;                      // current stock price is the starting stock price
		_st = s0;
		B  = 0;
			
		for (i = 1; i <= steps; ++i){

			U = PsiInv(MTUniform(0));

			mu += (rr - 0.5 * sigma2) * dt; // recalibrate mu for sigma2
			B = sqdt * sqrt(sigma2) * U; // brownian jump and scale with sigma2

			st *= exp(mu + B); // update stock price

			sigma2 = gmma * V2  +  beta * sigma2  +  alpha * B * B / dt; //garch

			////////////////////////////////////antithetic/////////////////////////////////////

			mu = (rr - 0.5 * _sigma2) * dt; // recalibrate mu for sigma2
			B = sqdt * sqrt(_sigma2) * -U; // brownian jump and scale with sigma2

			_st *= exp(mu + B);

			_sigma2 = gmma * V2  +  beta * _sigma2  +  alpha * B * B / dt; //garch

		}

		++n;

		for(i=0; i<11; ++i){
			calls[i].Pr = exp(-T*rr)*((std::max((st-calls[i].K),0.0)+std::max((_st-calls[i].K),0.0))/2);
			calls[i].EX = (calls[i].EX*(n-1) +calls[i].Pr)/(n);
			calls[i].EX2 = (calls[i].EX2*(n-1) +calls[i].Pr*calls[i].Pr)/(n);
		}



		if(n%100000==0){

			done = 1;

			for(i=0; i<11; ++i){
				calls[i].STD = sqrt(calls[i].EX2-calls[i].EX*calls[i].EX);
				calls[i].error = 1.96 * calls[i].STD/sqrt(n);
				printf("%2.3f, ",calls[i].error);
				if(calls[i].error>epsilon){
					done=0;
				}
			}
			std::cout<<std::endl;
		}
	}

	std::cout<<"\n\n--Prices of Options (95% Confidence Interval .01) && (K=50,60,...,150)"<<std::endl;

	for(i=0; i<11; ++i)
		printf("%8.4f, ",calls[i].EX);
		std::cout<<std::endl;

	std::cout<<"\n\n--Implied Volatility "<<std::endl;

	for(i=0; i<11; ++i){
		calls[i].impVol=ImpliedVol(T,s0,calls[i].K,rr,calls[i].EX);
		printf("%2.5f, ", calls[i].impVol);
	}
	
	std::cout<<std::endl;
	//

	return 0;
}



#include "Definitions.h"

double BlackScholes (const double& tau, double s, const double& k, const double& sigma, const double& r) {

   double value, d_plus, d_minus, v, d0;

   v = sigma * sqrt(tau);

   if (v < 0.0000001) {

      s *= exp(r*tau);
      value = (s > k ? s-k : 0);
      value *= exp(-r*tau);

   }else{
      d0 = (log (s/k) + r*tau) / v;
      d_plus  = d0 + 0.5 * v;
      d_minus = d0 - 0.5 * v;

      value = s * Psi(d_plus) - k * exp(-r*tau) * Psi(d_minus);
   }

   return (value>0?value:0);
}

double ImpliedVol (const double& tau, const double& s, const double& k, const double& r, const double& c) {

   int i;
   double step = 0.10, sigma = 0.0, value, value0;

   value0 = BlackScholes (tau, s, k, 0.0, r);

   if (value0 > c) {
   	//std::cout<<"Ouch! That really hurts. "<<std::endl;
   	return -1;}
   
   if (c < value0 + 0.0000001)
      return 0.0;

   value = value0;
   for (i = 1; i <= 5; i++) {

      while (value <= c) {

         sigma += step;
         value = BlackScholes (tau, s, k, sigma, r);
      }

      sigma -= step;
      value = BlackScholes (tau, s, k, sigma, r);

      step /= 10;
   }

   return sigma;

}