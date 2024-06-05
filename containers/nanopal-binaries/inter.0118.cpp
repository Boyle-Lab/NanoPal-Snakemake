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
	ifstream f_q;
	ifstream f_s;
	ofstream f_inter;

	f_q.open("Q.txt");
	f_s.open("S.txt");
	f_inter.open("inter.txt");

	string input;
	int q_lines = 0;
	int s_lines = 0;

	// Count lines in Q and S.
	for(int i=0;!f_q.eof();i++){
		getline(f_q,input);
		q_lines=i;
	}

	for(int i=0;!f_s.eof();i++){
		getline(f_s,input);
		s_lines=i;
	}

	f_q.close();
	f_q.clear();
	f_s.close();
	f_s.clear();
	f_q.open("Q.txt");
	f_s.open("S.txt");

	string **Q_info;
	Q_info=new string*[q_lines];
	for(int i=0;i!=q_lines;i++) Q_info[i]=new string[2];

	string **S_info;
	S_info=new string*[s_lines];
	for(int i=0;i!=s_lines;i++) S_info[i]=new string[2];

	int **Q_loc;
	Q_loc=new int*[q_lines];
	for(int i=0;i!=q_lines;i++) Q_loc[i]=new int[2];

	int **S_loc;
	S_loc=new int*[s_lines];
	for(int i=0;i!=s_lines;i++) S_loc[i]=new int[2];

	// These appear to parse the input lines like this:
	//
	//     chr1 42018289 42018289 cluster_L1_1_PAV 22.12.0.314156.0/1.9.9.2.2.16.0.0.cluster0_chr1_42018290_42018428_42018290_42018428_NA12878_PALMER + L1HS 4568 6023 L1HS
	//
	// By doing the equivalent of line.split(' ', 3) and storing the (±50)
	// start/end in _loc and the chromosome and tail of the string in _info:
	//
	//	Q_info[0] = ["chr1", "cluster_L1_1_PAV 22.12.0.314156.0/1.9.9.2.2.16.0.0.cluster0_chr1_42018290_42018428_42018290_42018428_NA12878_PALMER + L1HS 4568 6023 L1HS"]
	//	Q_loc[0]  = [42018239, 42018339]
	for(int i=0; i!=q_lines; i++){
		f_q >> Q_info[i][0];

		f_q >> Q_loc[i][0];
		Q_loc[i][0] = Q_loc[i][0]-50;

		f_q >> Q_loc[i][1];
		Q_loc[i][1] = Q_loc[i][1]+50;

		getline(f_q, Q_info[i][1]);
	}

	for(int i=0;i!=s_lines;i++){
		f_s >> S_info[i][0];

		f_s >> S_loc[i][0];
		S_loc[i][0] = S_loc[i][0]-50;

		f_s >> S_loc[i][1];
		S_loc[i][1] = S_loc[i][1]+50;

		getline(f_s, S_info[i][1]);
	}

	cout << q_lines << " " << s_lines << " ";

	if (q_lines) {
		cout << Q_loc[q_lines-1][1];
	} else {
		cout << "NA";
	}

	cout << " ";

	if (s_lines) {
		cout << S_loc[s_lines-1][1];
	} else {
		cout << "NA";
	}

	cout << endl;

	for(int i=0; i!=q_lines; i++){
		int flag=0;
		int qstart = Q_loc[i][0];
		int qend   = Q_loc[i][1];

		// For each line in Q, loop through every line in S
		for(int j=0; j != s_lines; j++) {
			int sstart = S_loc[j][0];
			int send   = S_loc[j][1];

			if(Q_info[i][0] == S_info[j][0]) { // if they're on the same chromosome
				if(!((qstart > send) || (qend < sstart))) { // if the ranges overlap
					// concat both lines (WITH modified start/end!!) and print
					f_inter << Q_info[i][0]
						<< '\t' << Q_loc[i][0]
						<< '\t' << Q_loc[i][1]
						<< '\t' << Q_info[i][1]
						<< '\t' << S_info[j][0]
						<< '\t' << S_loc[j][0]
						<< '\t' << S_loc[j][1]
						<< '\t' << S_info[j][1]
						<< endl;
					// sets the flag to avoid dumping the original line, but doesn't break or anything, so we might get multiple if multiple things overlap.
					flag=1;
				}
			}
		}

		// otherwise just dump the original line (BUT with the modified
		// start/end!!) from Q
		if(flag == 0) {
			f_inter << Q_info[i][0]
				<< '\t' << Q_loc[i][0]
				<< '\t' << Q_loc[i][1]
				<< '\t' << Q_info[i][1]
				<< endl;
		}
	}
}
