
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
	ifstream file1;
    ifstream file2;
    ofstream file5;
    
    
    file1.open("input_rm_cluster.txt");
file5.open("output_rm_cluster.txt");

    string input;
    int line1=0;
    for(int i=0;!file1.eof();i++){
        getline(file1,input);
        line1=i;
    }
    file1.close();
    file1.clear();
    file1.open("input_rm_cluster.txt");
    
    string rm[line1][5];
    string cluster[line1][2];
    int flag[line1];
    int rm_flag[line1][6];
    for(int i=0;i!=line1;i++){
        file1>>rm[i][0];
        file1>>rm[i][1];
	file1>>rm[i][2];
	file1>>rm[i][3];
	file1>>rm[i][4];
	flag[i]=0;
    }
    
    for(int i=0;i!=line1;i++){
        if(flag[i]==0){
            flag[i]=1;
            for(int j=0;j!=line1;j++){
                if(i!=j&&flag[j]==0){
                    if(rm[i][3]==rm[j][3]){
                        flag[j]=1;
			rm[i][4]=rm[i][4]+"/"+rm[j][4];
                    }
                }
            }
            file5<<rm[i][0]<<'\t'<<rm[i][1]<<'\t'<<rm[i][2]<<'\t'<<rm[i][3]<<'\t'<<rm[i][4]<<endl;
        }
    }
    
    
}

