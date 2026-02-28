#include "function.hpp"

void function(int A[20][20], int B[20][20], int C[20][20]){
    for (int i=0; i<20; i++){
#pragma HLS PIPELINE II=1
        for (int j=0; j<20; j++){
#pragma HLS UNROLL factor = 4
            C[i][j]=0;
        }
    }
    for (int i=0; i<20; i++){
#pragma HLS PIPELINE II=1
        for (int j=0; j<20; j++){
#pragma HLS UNROLL factor = 4
            C[i][j]+=A[i][j]*B[j][i];
        }
    }
}
