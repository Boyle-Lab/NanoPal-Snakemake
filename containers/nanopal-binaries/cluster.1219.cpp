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
	file1.open("input_cluster.txt");

	string clus[line][3];
	int loc[line][2];
	int flag[line];
	int num[line];
	int end[line][2];
	//getchar();
	for(int i=0;i!=line;i++){

		file1>> clus[i][0];
		file1>> loc[i][0];
		//loc[i][0]=loc[i][0]-50;
		file1>> loc[i][1];
		//loc[i][1]=loc[i][1]+50;
		file1>> input;
		file1>> input;
		end[i][0]=0;
		end[i][1]=0;
		file1>> input_int;
		if(input_int==3){
			end[i][1]=1;
		}
		else if(input_int==5){
			end[i][0]=1;
		}
		file1>>clus[i][2];
		getline(file1, clus[i][1]);
		flag[i]=0;
		num[i]=1;
	}

	//getchar();
	for(int i=0;i!=line;i++){
		if(flag[i]==0){
			for(int j=0;j!=line;j++){
				if(flag[j]==0&&i!=j&&clus[i][0]==clus[j][0]&&clus[i][2]==clus[j][2]){
					if(!((loc[i][0]>loc[j][1])||(loc[i][1]<loc[j][0]))){

						if(loc[j][0]<loc[i][0]) loc[i][0]=loc[j][0];
						if(loc[j][1]>loc[i][1]) loc[i][1]=loc[j][1];
						flag[j]=1;
						num[i]=num[i]+1;
						end[i][0]=end[i][0]+end[j][0];
						end[i][1]=end[i][1]+end[j][1];
					}
				}
			}
		}
	}

	ofstream file2;
	file2.open("clustered.txt");

	for(int i=0;i!=line;i++){
		if(flag[i]==0){
			file2<<clus[i][0]<<'\t'<<loc[i][0]<<'\t'<<loc[i][1]<<'\t'<<num[i]<<'\t'<<end[i][0]<<'\t'<<end[i][1]<<'\t'<<clus[i][2]<<endl;
		}
	}


}
