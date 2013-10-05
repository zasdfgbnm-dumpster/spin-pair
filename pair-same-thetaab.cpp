#include "pair.hpp" 

int main(){
	#pragma omp parallel for
	for(int n=1;n<=8;n++) {
		auto pairs = random_pairs(n);
		double w1 = pairs[0].omegar();
		stringstream fnstream;
		fnstream << "n=" << n << ",w1=" << w1/1_MHz << "MHz.txt";
		string suffix = fnstream.str();
		if(n<=4) {
			ofstream pairs_out(string("pairs-")+suffix);
			pairs_simulation(w1,pairs,pairs_out);
		}
		ofstream pes_out(string("pes-")+suffix);
		pes_simulation(w1,pairs,pes_out);
	}
}