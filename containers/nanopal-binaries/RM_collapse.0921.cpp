
#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cstdlib>

using namespace std;

int main(){
	ifstream f_in;
	ofstream f_out;


	f_in.open("input_rm_cluster.txt");
	f_out.open("output_rm_cluster.txt");

	// Count input lines.
	string input;
	int line1=0;
	for(int i=0; !f_in.eof(); i++){
		getline(f_in, input);
		line1=i;
	}
	f_in.close();
	f_in.clear();
	f_in.open("input_rm_cluster.txt");

	// This is *stack allocated* here, not heap allocated, so it explodes
	// when the line count is huge.
	string rm[line1][5];
	string cluster[line1][2];
	int flag[line1];

	/* int rm_flag[line1][6]; */
	for(int i=0;i!=line1;i++){
		f_in>>rm[i][0]; // chr10
		f_in>>rm[i][1]; // start
		f_in>>rm[i][2]; // end
		f_in>>rm[i][3]; // uuid
		f_in>>rm[i][4]; // NON/L1HS/etc
		flag[i]=0;      // always sets flag to 0
	}

	// Loop through every line of the file.
	for(int i=0;i!=line1;i++){
		// If this hasn't been marked already.
		if(flag[i]==0){
			// Mark it.
			flag[i]=1;

			// Loop through every line.
			for(int j=0;j!=line1;j++){

				// If we find any other lines that aren't already marked.
				if(i!=j&&flag[j]==0){

					// If the UUIDs match?
					if(rm[i][3]==rm[j][3]){
						// Mark the other line
						flag[j]=1;

						// And update this line to concat to i_type/j_type.
						rm[i][4]=rm[i][4]+"/"+rm[j][4];
					}
				}
			}
			// Print the line.
			f_out<<rm[i][0]<<'\t'<<rm[i][1]<<'\t'<<rm[i][2]<<'\t'<<rm[i][3]<<'\t'<<rm[i][4]<<endl;
		}
	}


}

