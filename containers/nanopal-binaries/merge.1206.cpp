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
#include <vector>

using namespace std;
using std::vector;

int main(){
    ifstream file1;
    
    ofstream file5;
    
    file1.open("summary.final.2.txt.input");
    file5.open("summary.final.2.txt.output");
    
    
    string input;
    int line;
    
    for(int i=0;!file1.eof();i++){
        getline(file1,input);
        line=i;
    }
    cout<<line<<endl;
    
    string info[line][8];
    int num[line][10];
    int flag[line];
    file1.close();   
    file1.clear();
    file1.open("summary.final.2.txt.input"); 
    for(int i=0;i!=line;i++){
        file1>>info[i][0];
        file1>>num[i][0];
        file1>>num[i][1];
        file1>>num[i][2];
        file1>>num[i][3];
        file1>>num[i][4];
        file1>>num[i][5];
        file1>>num[i][6];
        file1>>num[i][7];
        file1>>num[i][8];
        file1>>info[i][1];
        file1>>info[i][2];
        file1>>info[i][3];
        file1>>info[i][4];
        file1>>info[i][5];
        file1>>info[i][6];
        file1>>info[i][7];
        file1>>num[i][9];
        flag[i]=1;
    }
  // cout<<info[1][7]<<endl; 
    
    for(int i=0;i!=line;i++){
  //      cout<<i<<endl;
	if(flag[i]==1){
            for(int j=0;j!=line;j++){
//cout<<j<<endl;
                if(flag[j]==1&&i!=j){
                    if(info[i][0]==info[j][0]){
                        flag[j]=0;
//cout<<"yes "<<i<<" "<<j<<endl;
                        num[i][1]=num[i][1]+num[j][1];
                        num[i][2]=num[i][2]+num[j][2];
//cout<<"num[i][2]"<<num[i][2]<<endl;
//cout<<"num[j][2]"<<num[j][2]<<endl;
                        num[i][3]=num[i][3]+num[j][3];
                        num[i][4]=num[i][4]+num[j][4];
                        num[i][5]=num[i][5]+num[j][5];
                        num[i][6]=num[i][6]+num[j][6];
                        num[i][7]=num[i][7]+num[j][7];
                        num[i][8]=num[i][8]+num[j][8];
                    }
                }
            }
            file5<<info[i][0]<<'\t'<<num[i][0]<<'\t'<<num[i][1]<<'\t'<<num[i][2]<<'\t'<<num[i][3]<<'\t'<<num[i][4]<<'\t'<<num[i][5]<<'\t'<<num[i][6]<<'\t'<<num[i][7]<<'\t'<<num[i][8]<<'\t'<<info[i][1]<<'\t'<<info[i][2]<<'\t'<<info[i][3]<<'\t'<<info[i][4]<<'\t'<<info[i][5]<<'\t'<<info[i][6]<<'\t'<<info[i][7]<<'\t'<<num[i][9]<<endl;
        }
    }
    
    
    
    file1.close();
    file5.close();

}
