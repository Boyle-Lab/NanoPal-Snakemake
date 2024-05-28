#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main() {
	for (string cigar, input; getline(cin, input);) {
		cigar=input+"0E";

		int number=0;
		int sc1, sc2, I, M, D, F, h1, h2, eq;
		sc1=sc2=I=M=D=F=h1=h2=eq=0;
		char cig;

		stringstream ss_cigar1;
		ss_cigar1.clear();
		ss_cigar1.str(cigar);

		ss_cigar1>>number;
		ss_cigar1>>cig;

		for(;!(cig=='E'&&number==0);){
			if(cig=='S'&&F==0) {sc1=sc1+number;F=1;}
			else if(cig=='H'&&F==0) h1=h1+number;
			else if (cig=='M') {M=M+number;F=1;}
			else if (cig=='='||cig=='X') {eq=eq+number;F=1;}
			else if (cig=='I') {I=I+number;F=1;}
			else if (cig=='D'||cig=='N') {D=D+number;F=1;}
			else if (cig=='S'&&F==1) sc2=sc2+number;
			else if (cig=='H'&&F==1) h2=h2+number;
			ss_cigar1>>number;
			ss_cigar1>>cig;
		}

		cout
			<< "sc1" << '\t' << sc1 << '\t'
			<< "M" << '\t' << M << '\t'
			<< "eq" << '\t' << eq << '\t'
			<< "I" << '\t' << I << '\t'
			<< "D" << '\t' << D << '\t'
			<< "sc2" << '\t' << sc2 << '\t'
			<< "h1" << '\t' << h1 << '\t'
			<< "h2" << '\t' << h2 << endl;
	}

	return 0;
}
