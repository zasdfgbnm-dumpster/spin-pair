#include <limits>
#include "pair.hpp"

constexpr int threads = 4;
constexpr int splits = 200;
constexpr int per_thread_n = 100000;
constexpr int output_per = 10000;

const double unit[2] = { 1_MHz,1_kHz };
const int minimum[2] = { -16,-25 };
const int maximum[2] = { 31,40 };

constexpr int N = threads*per_thread_n;
const double bs[2] = { static_cast<double>(maximum[0]-minimum[0])/splits, static_cast<double>(maximum[1]-minimum[1])/splits };
const double bs2 = bs[0]*bs[1];

inline int get_index(double value,int dim){
	value = value/unit[dim];
	value -= minimum[dim];
	value /= bs[dim];
	return static_cast<int>(value);
}

inline double get_coord(int i,int dim) {
	return bs[dim]*(0.5+i)+minimum[dim];
}

int dist_count_per_thread[threads][splits][splits] = { 0 };
int dist_count[splits][splits] = { 0 };

inline void mc1(int thread) {
	auto pair = random_pair();
	double value[2] = { pair.coeff_Sxpes(),pair.coeff_Sz0_Szpes() };
	int index[2] = { get_index(value[0],0),get_index(value[1],1) };
	dist_count_per_thread[thread][index[0]][index[1]]++;
}

int main(){
	#pragma omp parallel for
	for(int i=0;i<threads;i++){
		for(int j=0;j<per_thread_n;j++) {
			mc1(i);
			if( j%output_per == 0 )
				cout << "[" << i << "] " << j/output_per << "/" << per_thread_n/output_per << endl;
		}
	}
	
	for(int i=0;i<threads;i++)
		for(int j=0;j<splits;j++)
			for(int k=0;k<splits;k++)
				dist_count[j][k] += dist_count_per_thread[i][j][k];
	
	ofstream out("result.txt");
	for(int j=0;j<splits;j++) {
		for(int k=0;k<splits;k++) {
			double value = static_cast<double>(dist_count[j][k])/N/bs2;
			double x = get_coord(j,0);
			double y = get_coord(k,1);
			out << x << "\t" << y << "\t" << value << endl;
		}
	}
	out.close();
}