#include "reactivity.hpp"
#include "HX1.hpp"
#include "HX2.hpp"
#include "thermal_hydraulics.hpp"
#include "neutronics.hpp"
#include <hls_stream.h>

void simulation_test(float rho_insertion, float &avg_rho, float &phi_mid) {
    // Define internal arrays
    float y_n[length_neutr];
    float y_th[length_th];
    float y_hx1[length_hx];
    float y_hx2[length_hx];
    float rho[N];

    // Initialize all arrays to zero
    for (int i = 0; i < length_neutr; i++) y_n[i] = 0;
    for (int i = 0; i < length_th; i++) y_th[i] = 0;
    for (int i = 0; i < length_hx; i++) y_hx1[i] = 0;
    for (int i = 0; i < length_hx; i++) y_hx2[i] = 0;
    for (int i = 0; i < N; i++) rho[i] = 0;

    // Initialize parameters
    Param_Neutronics params_neutr;

    // Perform simulation loop
    for (int step = 0; step < time_span; step++) {
        neutronics(y_n, rho, step, params_neutr);
        thermal_hydraulics(y_th, y_n, y_hx1[0], step);
        HX1(y_hx1, y_th[N - 1], step);
        HX2(y_hx2, y_hx1[Nx - 1], step);
        reactivity(y_th, y_th + N, step, rho, rho_insertion);
    }

    // Calculate the average rho
    avg_rho = 0.0;
    for (int i = 0; i < 200; i++) {
        avg_rho += rho[i];
    }
    avg_rho /= 200.0;

    // Get the mid-value of phi
    phi_mid = y_n[100];
}


//void simulation_test(
//    float y_n[length_neutr],
//    float y_th[length_th],
//    float y_hx1[length_hx],
//    float y_hx2[length_hx],
//    float rho[N],
//    float &rho_mid,         // Output for rho[N/2]
//    float &phi1_mid        // Output for y_n[100]
////    float &phi2_mid,        // Output for y_n[300]
////    float &fuel_mid,        // Output for y_th[100]
////    float &graphite_mid     // Output for y_th[300]
//) {
//    #pragma HLS INTERFACE s_axilite port=return bundle=CONTROL_BUS
//    #pragma HLS INTERFACE s_axilite port=rho_mid bundle=CONTROL_BUS
//    #pragma HLS INTERFACE s_axilite port=phi1_mid bundle=CONTROL_BUS
////    #pragma HLS INTERFACE s_axilite port=phi2_mid bundle=CONTROL_BUS
////    #pragma HLS INTERFACE s_axilite port=fuel_mid bundle=CONTROL_BUS
////    #pragma HLS INTERFACE s_axilite port=graphite_mid bundle=CONTROL_BUS
//    #pragma HLS INTERFACE m_axi depth=length_neutr port=y_n offset=slave
////    #pragma HLS INTERFACE m_axi depth=length_th port=y_th offset=slave
////    #pragma HLS INTERFACE m_axi depth=length_hx port=y_hx1 offset=slave
////    #pragma HLS INTERFACE m_axi depth=length_hx port=y_hx2 offset=slave
//    #pragma HLS INTERFACE m_axi depth=N port=rho offset=slave
//
//    #pragma HLS array_partition variable=y_n complete
////    #pragma HLS array_partition variable=y_th complete
////    #pragma HLS array_partition variable=y_hx1 complete
////    #pragma HLS array_partition variable=y_hx2 complete
//    #pragma HLS array_partition variable=rho complete
//
//    // Initialize all arrays to zero
//    for (int i = 0; i < length_neutr; i++) y_n[i] = 0;
//    for (int i = 0; i < length_th; i++) y_th[i] = 0;
//    for (int i = 0; i < length_hx; i++) y_hx1[i] = 0;
//    for (int i = 0; i < length_hx; i++) y_hx2[i] = 0;
//    for (int i = 0; i < N; i++) rho[i] = 0;
//
//    Param_Neutronics params_neutr;
//    for (int step = 0; step < time_span; step++) {
//        // Call functions (assuming they modify y_n, y_th, y_hx1, y_hx2, and rho arrays)
//        neutronics(y_n, rho, step, params_neutr);
//        thermal_hydraulics(y_th, y_n, y_hx1[0], step);
//        HX1(y_hx1, y_th[N - 1], step);
//        HX2(y_hx2, y_hx1[Nx - 1], step);
//        reactivity(y_th, y_th + N, step, rho);
//    }
//
//    // Set the specific outputs
//    rho_mid = 0.0;
//    for (int i=0; i<200; i++){
//    	rho_mid+=rho[i];
//    }
//    rho_mid/=200;
//    phi1_mid = y_n[100];
////    phi2_mid = y_n[300];
////    fuel_mid = y_th[100];
////    graphite_mid = y_th[300];
//}



//void simulation_test (float y_n[length_neutr], float y_th[length_th], float y_hx1[length_hx], float y_hx2[length_hx], float rho[N]){
//#pragma HLS array_partition variable=y_n complete
//    #pragma HLS array_partition variable=y_th complete
//    #pragma HLS array_partition variable=y_hx1 complete
//    #pragma HLS array_partition variable=y_hx2 complete
//    #pragma HLS array_partition variable=rho complete
//
//    // Initialize all arrays to zero
//    for (int i = 0; i < length_neutr; i++) y_n[i] = 0;
//    for (int i = 0; i < length_th; i++) y_th[i] = 0;
//    for (int i = 0; i < length_hx; i++) y_hx1[i] = 0;
//    for (int i = 0; i < length_hx; i++) y_hx2[i] = 0;
//    for (int i = 0; i < N; i++) rho[i] = 0;
//
//	Param_Neutronics params_neutr;
//	for (int step = 0; step< time_span; step++){
//		neutronics(y_n, rho, step, params_neutr);
//		thermal_hydraulics(y_th, y_n, y_hx1[0], step);
//		HX1(y_hx1, y_th[N - 1], step);
//		HX2(y_hx2, y_hx1[Nx - 1], step);
//		reactivity(y_th, y_th + N, step, rho);
//	}
//}

