#include "pair.hpp" 

int main(){
	for(int n=1;n<=4;n++) {
		auto pairs = random_pairs(n);
		double cor = pairs[0].omegar();
		double ratios[] = { 0.9,0.95,0.98,1,1.02,1.05,1.1 };
		for(auto &ratio:ratios) {
			double w1 = cor*ratio;
			stringstream fnstream;
			fnstream << "n=" << n << ",w1=" << w1/1_MHz << "MHz.txt";
			string suffix = fnstream.str();
			ofstream pairs_out(string("pairs-")+suffix);
			pairs_simulation(w1,pairs,pairs_out);
		}
	}
}