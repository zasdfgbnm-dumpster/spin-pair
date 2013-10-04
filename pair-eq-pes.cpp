#include "pair.hpp" 

int main(){
	for(int n=1;n<=4;n++) {
		auto pairs = random_pairs(n);
		std::vector<double> omegars(n);
		transform(pairs.begin(),pairs.end(),omegars.begin(),[](spin_pair p){return p.omegar();});
		for(auto &cor:omegars) {
			double ratios[3] = { 0.9,1,1.1 };
			for(auto &ratio:ratios) {
				double w1 = cor*ratio;
				stringstream fnstream;
				fnstream << "n=" << n << ",w1=" << w1/1_MHz << "MHz.txt";
				string suffix = fnstream.str();
				ofstream pairs_out(string("pairs-")+suffix);
				ofstream pes_out(string("pes-")+suffix);
				pairs_simulation(w1,pairs,pairs_out);
				pes_simulation(w1,pairs,pes_out);
			}
		}
	}
}