//copyright by ArthurZhou @ UMich&FUDAN&HUST
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
#include <numeric>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

using namespace std;

int main(){

	/* int bin=50; */
	ifstream file1;

	// Count lines in input file.
	file1.open("input_cluster.txt");

	int input_int;
	string input;
	int line=0;

	for(int i=0;!file1.eof();i++){
		getline(file1,input);
		line=i;
	}

	file1.close();
	file1.clear();


	// Parse input file.  Lines look like:
	//                                                    |
	// chr1    10014616        10014717         NON 1 3 + | id1
	// chr1    100784028       100784129        NON 1 3 + | id2
	// chr1    101123651       101123752        NON 1 5 + | id3
	// chr1    101418520       101418621        NON 1 5 + | id4
	// chr1    101476959       101477060        NON 1 5 + | id5
	file1.open("input_cluster.txt");

	// clus - [(chr1, tail of line…, +), …]
	// loc  - [(start, end), …]
	// flag - [0, …]
	// num  - [1, …]
	// end  - [(0, 1), …]    (0, 1) for 3, or (1, 0) for 5
	string clus[line][3];
	int    loc[line][2];
	int    flag[line];
	int    num[line];
	int    end[line][2];
	//getchar();
	for(int i=0;i!=line;i++){
		file1>> clus[i][0];
		file1>> loc[i][0];
		// Already a 100bp range from inter step, so it doesn't extend
		// another 100 here I think.
		//loc[i][0]=loc[i][0]-50;
		file1>> loc[i][1];
		//loc[i][1]=loc[i][1]+50;
		file1>> input;
		file1>> input;
		end[i][0]=0;
		end[i][1]=0;
		file1>> input_int;
		if(input_int == 3){
			end[i][1]=1;
		}
		else if(input_int == 5){
			end[i][0]=1;
		}
		file1>> clus[i][2];
		getline(file1, clus[i][1]);
		flag[i]=0;
		num[i]=1;
	}

	//getchar();
	// Loop through each line.
	for(int i =0;i!=line;i++) {
		// If this line hasn't already been clustered by an earlier line
		if(flag[i] == 0) {
			// Loop through all the other lines.
			for(int j=0; j !=line;j++) {
				if(flag[j] == 0 // If a line hasn't already been clustered.
				   && i!=j      // don't cluster with itself
				   && clus[i][0] == clus[j][0]    // if they're on the same chromosome
				   && clus[i][2] == clus[j][2]) { // and the same strand
					if(!(
						// i_start > j_end || i_end < j_start
						// i.e. they do NOT overlap,
						// but this is NEGATED above, so
						// this fires when they DO overlap.

						// i_start > j_end
						(loc[i][0] > loc[j][1])

						// i_end   < j_start
					     || (loc[i][1] < loc[j][0])
					)){

						// j_start < i_start
						if(loc[j][0] < loc[i][0]) loc[i][0] = loc[j][0]; // extend i_start

						// j_end > i_end
						if(loc[j][1] > loc[i][1]) loc[i][1] = loc[j][1]; // extend i_end

						// mark j as processed, but NOT i
						flag[j] = 1;

						// incf count of things we've merged into i
						num[i] = num[i]+1;

						// I *think* this is counting the total number of reads on
						// the 5' and 3' ends.
						end[i][0] = end[i][0] + end[j][0];
						end[i][1] = end[i][1] + end[j][1];

						// Concat the tails
						clus[i][1] = clus[i][1] + "/" + clus[j][1];
					}
				}
			}
		}
	}

	ofstream file2;
	file2.open("clustered.txt");

	for(int i=0;i!=line;i++){
		// for all the merged-into lines
		if(flag[i]==0){
			// print fields
			file2 << clus[i][0] << '\t'  // chr
			      << loc[i][0]  << '\t'  // (merged/widened) start
			      << loc[i][1]  << '\t'  // (merged/widened) end
			      << num[i]     << '\t'  // number of things we've merged into this line
			      << end[i][0]  << '\t'  // number of reads on 5' end
			      << end[i][1]  << '\t'  // number of reads on 3' end
			      << clus[i][2] << '\t'  // strand
			      << clus[i][1] << endl; // concatenated tails
		}
	}


}
