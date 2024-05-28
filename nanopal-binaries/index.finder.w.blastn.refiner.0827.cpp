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

using namespace std;
using std::vector;


int main(int argc, char *argv[]){
    
    
    cout<<"Example: ./index.finder.w.blastn.refiner.0827.o --mei 5900 --primer 30 --read 100"<<endl;
    
    
    int bin_mei, bin_prim, bin_read;
    int flag_mei=0;
    int flag_prim=0;
    int flag_read=0;
    string s_mei, s_prim, s_read;
    
    for(int i=1;i!=argc;++i){
        if(strncmp(argv[i],"--mei",5)==0){
            s_mei = argv[i+1];
            flag_mei = 1;
            bin_mei = std::stoi(s_mei);
        }
        
        if(strncmp(argv[i],"--primer",8)==0){
            s_prim = argv[i+1];
            flag_prim = 1;
            bin_prim = std::stoi(s_prim);
        }
        
        if(strncmp(argv[i],"--read",6)==0){
            s_read = argv[i+1];
            flag_read = 1;
            bin_read = std::stoi(s_read);
        }
    }
    
    if(flag_mei==0||flag_prim==0||flag_read==0){
        cout<<"Please have the correct input."<<endl;
        exit(1);
    }
    
    
    ifstream file1;
    ifstream file2;
    ifstream file3;
    
    ofstream file5;
    //ofstream file11;
    //ofstream file12;
    //ofstream file13;
    //ofstream file14;
    
    file1.open("Query.txt");
    file2.open("Index.txt");
    file3.open("Blastn.txt");
    
    file5.open("IF.output.txt");
    //file11.open("output.5.plus.txt");
    //file12.open("output.5.minus.txt");
    //file13.open("output.3.plus.txt");
    //file14.open("output.3.minus.txt");
    
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
    
    int line3;
    
    for(int i=0;!file3.eof();i++){
        getline(file3,input);
        line3=i;
    }
    
    cout<< line1<<" "<<line2<<" "<<line3<<endl;
    
    file1.close();
    file1.clear();
    file2.close();
    file2.clear();
    file3.close();
    file3.clear();
    file1.open("Query.txt");
    file2.open("Index.txt");
    file3.open("Blastn.txt");
    
    
    
    string **Q_info;
    Q_info=new string*[line1];
    for(int i=0;i!=line1;i++) Q_info[i]=new string[1];
    int **Q_loc;
    Q_loc=new int*[line1];
    for(int i=0;i!=line1;i++) Q_loc[i]=new int[5];
    
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
    
    string **B_info;
    B_info=new string*[line3];
    for(int i=0;i!=line3;i++) B_info[i]=new string[3];
    int **B_loc;
    B_loc=new int*[line3];
    for(int i=0;i!=line3;i++) B_loc[i]=new int[4];
    
    
    
    for(int i=0;i!=line1;i++){
        file1>>Q_info[i][0];
        file1>>Q_loc[i][0];
        file1>>Q_loc[i][1];
        file1>>Q_loc[i][2];
        file1>>Q_loc[i][3];
        file1>>Q_loc[i][4];
        //getline(file1,Q_info[i][1]);
    }
    

    for(int i=0;i!=line2;i++){
        for(int j=0;j!=col2;j++){
            file2>>S_info[i][j];
        }
    }


    for(int i=0;i!=line3;i++){
	    
        file3>>B_info[i][0];
        file3>>B_info[i][1];
        file3>>B_info[i][2];
        
        file3>>B_loc[i][0];
        file3>>B_loc[i][1];
        file3>>B_loc[i][2];
        file3>>B_loc[i][3];

	//getline(file1,Q_info[i][1]);
    }
    
    cout<<line1<<" "<<line2<<" "<<line3<<" "<<Q_info[line1-1][0]<<" "<<S_info[line2-1][1]<<" "<<B_info[line3-1][1]<<endl;
    
    
    
    for(int i=0;i!=line1;i++){
        
        file5<<Q_info[i][0]<<'\t'<<Q_loc[i][0]<<'\t'<<Q_loc[i][1]<<'\t'<<Q_loc[i][2]<<'\t'<<Q_loc[i][3]<<'\t'<<Q_loc[i][4];
        
        
        
        int flag_5p=0;
        int flag_5m=0;
        int flag_3p=0;
        int flag_3m=0;
        for(int j=0;j!=line3;j++){
            string::size_type idx;
            idx=B_info[j][1].find(Q_info[i][0]);
            if(idx != string::npos && B_loc[j][1]>bin_mei){
                if(B_loc[j][2]<(bin_prim+bin_read)||(B_loc[j][3]<(bin_prim+bin_read))){
                    if(B_loc[j][2]<B_loc[j][3]){
                        flag_5p=1;
                    }
                    else if(B_loc[j][2]>B_loc[j][3]){
                        flag_5m=1;
                    }
                }
                if(B_loc[j][2]>(Q_loc[i][0]-bin_read-bin_prim)||(B_loc[j][3]>(Q_loc[i][0]-bin_read-bin_prim))){
                    if(B_loc[j][2]<B_loc[j][3]){
                        flag_3p=1;
                    }
                    else if(B_loc[j][2]>B_loc[j][3]){
                        flag_3m=1;
                    }
                }
            }
        }
        
        file5<<'\t'<<flag_5p<<'\t'<<flag_5m<<'\t'<<flag_3p<<'\t'<<flag_3m;
        
        int flag=0;
        for(int j=0;j!=line2;j++){
            if(Q_info[i][0]==S_info[j][0]){
                
                for(int k=1;k!=col2;k++){
                    file5<<'\t'<<S_info[j][k];
                }
                file5<<endl;
                flag=1;
                break;
            }
        }
        if(flag==0){
            //file5<<Q_info[i][0]<<'\t'<<Q_info[i][1];
            for(int k=1;k!=col2;k++){
                file5<<'\t'<<"EMPTY";
            }
            file5<<endl;
        }
        
        
        
    }
    
    
    
}
