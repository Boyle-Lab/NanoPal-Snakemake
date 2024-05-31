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

	cout << q_lines
		<< " " << s_lines
		<< " " << Q_loc[q_lines-1][1]
		<< " " << S_loc[s_lines-1][1]
		<< endl;

	for(int i=0; i!=q_lines; i++){
		int flag=0;
		for(int j=0; j != s_lines; j++){
			if(Q_info[i][0] == S_info[j][0]) {
				if(!((Q_loc[i][0] > S_loc[j][1]) || (Q_loc[i][1] < S_loc[j][0]))) {
					f_inter << Q_info[i][0]
						<< '\t' << Q_loc[i][0]
						<< '\t' << Q_loc[i][1]
						<< '\t' << Q_info[i][1]
						<< '\t' << S_info[j][0]
						<< '\t' << S_loc[j][0]
						<< '\t' << S_loc[j][1]
						<< '\t' << S_info[j][1]
						<< endl;
					flag=1;
				}
			}
		}

		if(flag == 0) {
			f_inter << Q_info[i][0]
				<< '\t' << Q_loc[i][0]
				<< '\t' << Q_loc[i][1]
				<< '\t' << Q_info[i][1]
				<< endl;
		}
	}
}
