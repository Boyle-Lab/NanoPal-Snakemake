
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
#include <iomanip>

using namespace std;

int main(){
    ifstream file1;
    ifstream file2;
    ofstream file5;
    
    
    file1.open("Query.txt");
    file2.open("Index.txt");
    file5.open("IF.output.txt");
    
    string input;
    int line1,line2;
    
    //column
    int col2=0;
    
    vector<string> item;
    
    getline(file2,input);
    item.push_back(input);
    for(auto it=item.begin();it!=item.end();it++){
        istringstream istr(*it);
        string str;
        while (istr >> str){
            col2++;
            
        }
    }
    cout << col2 <<endl;
    
    
    for(int i=0;!file1.eof();i++){
        getline(file1,input);
        line1=i;
    }
    
    
    for(int i=0;!file2.eof();i++){
        getline(file2,input);
        line2=i;
    }
    line2=line2+1;
    
    cout<< line1<< " "<<line2<<endl;
    
    file1.close();
    file1.clear();
    file2.close();
    file2.clear();
    file1.open("Query.txt");
    file2.open("Index.txt");
    
    
    
    
    string **Q_info;
    Q_info=new string*[line1];
    for(int i=0;i!=line1;i++) Q_info[i]=new string[2];
    string **S_info;
    S_info=new string*[line2];
    for(int i=0;i!=line2;i++) S_info[i]=new string[col2];
    /*
     int **Q_loc;
     Q_loc=new int*[line1];
     for(int i=0;i!=line1;i++) Q_loc[i]=new int[2];
     int **S_loc;
     S_loc=new int*[line2];
     for(int i=0;i!=line2;i++) S_loc[i]=new int[2];
     */
    for(int i=0;i!=line1;i++){
        file1>>Q_info[i][0];
        
        getline(file1,Q_info[i][1]);
    }
    
    
    for(int i=0;i!=line2;i++){
        for(int j=0;j!=col2;j++){
            file2>>S_info[i][j];
        }
    }
    
    cout<<line1<<" "<<line2<<" "<<Q_info[line1-1][1]<<" "<<S_info[line2-1][1]<<endl;
    
    
    
    for(int i=0;i!=line1;i++){
        int flag=0;
        for(int j=0;j!=line2;j++){
            if(Q_info[i][0]==S_info[j][0]){
                file5<<Q_info[i][0]<<'\t'<<Q_info[i][1];
                for(int k=1;k!=col2;k++){
                    file5<<'\t'<<S_info[j][k];
                }
                file5<<endl;
                flag=1;
                break;
            }
        }
        if(flag==0){
            file5<<Q_info[i][0]<<'\t'<<Q_info[i][1];
            for(int k=1;k!=col2;k++){
                file5<<'\t'<<"EMPTY";
            }
            file5<<endl;
        }
    }
    
    
    
}
