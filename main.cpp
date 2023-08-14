// the last dance.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include<vector>
#include "polymer_class.h"
#include "polymer_generation.h"
#include "math_functions.h"
#include "helix_generation.h"
#include "simulations.h"
//#include "acceptance.h"
#include <fstream>
#include <chrono>




//TO DO.  
////////////////////////////////////////////////////////////////////////////////////////

// stop sampling alpha and beta. too much work not worth
// 
// SAMPLE UNLINK REGION NEEDS TO BE CHANGED AGAIN IF ZIP IS INCORPORATED. when doing the simulation where helix length
// was varied makes sim think that standard structures of length 3 are zipped which is not the intention for that simulation.
// 
// how does swivel work with zipped structures? before we assumed that the centre of mass was in the plane of 
// one of the base pairs however if we have an even number of base pairs e.g 4 this will not be the case, so is 
// it working correctly? probably not.
// 
// 
// can make metropolis structure more efficient. e,g only calculate search results, extensible structures when 
// a certain move has been accepted. after a series of rejections it should not be re calcualted.
// 
// 
// change zip so that we always sample the side of zip from which we can actually extend the helix. need to modify
// the acceptance rule potentially.
// 
// 
// structure id needs to be cleaned up i think. 
// // 
// 
// limiting case for zip where we zip the two tail ends of the chain. am i protected against that?
// 
// 
////////////////////////////////////////////////////////////////////////////////////////



int main()
{
    srand(time(NULL));
    std::int64_t init{ rand()};
    init = -1*init;
    ran2(&init);

    int NMC{ 100 }, i{ 0 }, N_monomers(10), NMC2{ 10000 }, init_update_rate{10};
    std::cout << "NMC is " << NMC << std::endl;
    //polymer p0(N_monomers, true);
    int big_NMC{ 10000 };

    std::vector<int> KH;
    bool overstrecth;
    std::vector<int> s1{ 0,5,14, 19};
    polymer short0(20, true);
    short0.set_search_results(s1);
    short0.link_move(overstrecth);
    short0.ovito_bonds("long_helix_8link_0");

    //std::vector<int> s2{ 0,2,17,19 };
    //polymer long0(20, true);
    //long0.set_search_results(s2);
    //long0.link_move(overstrecth);
    //long0.ovito_bonds("long_helix");



    //polymer p0(30, true);
    //p0.set_search_results(always_linked_structures);
    //p0.set_link_acceptance_prefactor(1000000);
    bool overstr;
    //for (int i = 0; i < 100; i++)
    //{
    //    std::cout << "ITERATION " << i << std::endl;
    //    int c{ 0 };
    //    polymer p0(30, true);
    //    p0.set_search_results(always_linked_structures);
    //    p0.set_link_acceptance_prefactor(1000000);

    //    p0.ovito_bonds("s0");
    //    p0.link_move(overstr, 3);
    //    p0.ovito_bonds("s1");
    //    p0.link_move(overstr, 3);
    //    if (overstr != true) { c++; }
    //    p0.ovito_bonds("s2");

    //    int k{ 0 };
    //    while (c != 2) {

    //        p0.link_move(overstr, 3);
    //        p0.ovito_bonds("s" + std::to_string(3 + 2*k));
    //        if (overstr != true) { c++; }
    //        p0.swivel_move();
    //        p0.ovito_bonds("s" + std::to_string(4 + 2*k));
    //        k++;
    //    }


    //    if (c == 2){
    //        std::cout << "stop";
    //    }
    //}

    std::vector<std::vector<int>> always_linked_structures{ { 0,3,8,11 }, { 13,15,20,22 }, { 4,6,16,18 }, {25,35,45,55}, {37,40,60,63}, {100,120,140,160},{75,85,170,180},{90,95,130,135} };

    std::vector<std::vector<int>> structurez{ {4,6,10,12}, {2,6,10,14} };
    std::vector<int> zip_r(4);
    polymer zippy(200, true);
    zippy.set_search_results(always_linked_structures);
    zippy.set_link_acceptance_prefactor(1000000);

    int c{ 0 }, unl{ 0 };
    zippy.ovito_bonds("s0");
    int k{ 1 };
    while (zippy.get_num_helices() < 4) {
        zippy.link_move(overstr, 3);
        zippy.ovito_bonds("s"+std::to_string(k));
        if (zippy.get_num_helices() > 1 && unl==0) {
            zippy.unlink_move(overstr);
            unl++;
        }
        zippy.swivel_move();
        zippy.ovito_bonds("s" + std::to_string(k+1));
        k = k + 2;
    }
    zippy.ovito_bonds("s1");
    zippy.link_move(overstr, 3);
    zippy.sample_unzip_region(zip_r, 0);
    bool success;
    int s_index, side;

    //std::vector<int> helix{ 0,2,7,9 };
    //polymer p0(10, true);
    //p0.set_search_results(helix);
    //p0.ovito_bonds("visUNBOUND");
    //bool yamum{ false };
    //p0.link_move(yamum);
    //p0.ovito_bonds("visBOUND");
    //p0.set_link_acceptance_prefactor(1.0);
    // reference region to calculate prefactor
    std::vector<int> regions{ 9,11,N_monomers - 3,N_monomers - 1 };
    //////////////////////////////////////////////////////////////////////////////////////////
    //std::string path0;

    //for (int i{ 0 }; i < big_NMC; i++) {
    //    path0 = "Results/swivel_sims/rep";
    //    path0 = path0 + std::to_string(i) + "_";
        //swivel_sim(N_monomers, NMC, init_update_rate);

    //}


    //////////////////////////////////////////////////////////////////////////////////////////
    // biased equilibration sim. starting from several different prefactors, show how we reach equilibrium. first run is arbitrary then
    // find equilibrium.
    //std::string path0;
    //for (int i = 0; i < big_NMC; i++)
    //{
    //    path0 = "Results/Equilibrium_Sim/rep";
    //    path0 = path0 + std::to_string(i) + "_";
    //    //std::cout << path0 << std::endl;
    //    equilibration_sim(10, 100, 10, path0);

    //}
    //////////////////////////////////////////////////////////////////////////////////////////
    //std::string path0;

    //for (int i{ 0 }; i < big_NMC; i++) {
    //    path0 = "Results/varying_helix_length_sims/rep";
    //    path0 = path0 + std::to_string(i) + "_";
    //    long_helix_simulation(10, NMC, init_update_rate, path0);

    //}

    ////////////////////////////////////////////////////////////////////////////////////////////
    //std::string path0;
    //for (int i{ 0 }; i < 1; i++) {
    //    path0 = "Results/varying_dangle_sims/rep";
    //    path0 = path0 + std::to_string(i) + "_";
    //    long_dangle_simulation(N_monomers, 10000, 100, path0);

    //}
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    //print_1d_doub_vec(generate_linear_array(0, 2, 9));
    //std::string path0;
    //for (int i{ 0 }; i < 1; i++) {
    //    path0 = "Results/varying_dangle_sims/rep";
    //    path0 = path0 + std::to_string(i) + "_";
    //    WCA_sim(N_monomers, 300, 300, path0);

    //}
    //////////////////////////////////////////////////////////////////////////////////////////

    //std::string path0;

    //for (int i{ 0 }; i < big_NMC; i++) {
    //    path0 = "Results/varying_linker_sims/rep";
    //    path0 = path0 + std::to_string(i) + "_";
    //    long_linker_simulation(N_monomers, 1000, 100, path0);

    //}

    //long_linker_simulation(N_monomers, NMC, init_update_rate, "Results/varying_linker_sims/");
    //////////////////////////////////////////////////////////////////////////////////////////

    // determine appropriate weight factor
    //p0.set_link_acceptance_prefactor(0.01);
    //p0.set_search_results(regions);
    //determine_g0(&p0, NMC, 2);
    //correlation_length(p0.get_link_acceptance_prefactor(), N_monomers, regions, NMC);
    ////// varying dangle
    //double acc_prefac{ p0.get_link_acceptance_prefactor() };
    //std::cout << "varying dangle length" << std::endl;
    //vary_dangle_sim(acc_prefac, N_monomers, NMC2, regions,"Results/varying dangle sims/");
    //////////////////////////////////////////////////////////////////////////////////////
    //std::vector<int> regions{ 0,2,N_monomers - 3,N_monomers - 2 };
    //regions = { 0,2,N_monomers - 3,N_monomers - 1 };
    //polymer p1(N_monomers, true);
    //p1.set_link_acceptance_prefactor(0.01);
    //p1.set_search_results(regions);
    //determine_g0(&p1, NMC, 100);
    //correlation_length(p1.get_link_acceptance_prefactor(), N_monomers, regions, NMC);
    //vary_WCA_potential(p1.get_link_acceptance_prefactor(), N_monomers, NMC2, regions, "Results/WCA_potential/");
    ////////////////////////////////////////////////////////////////////////////////////////
    //std::vector<std::vector<int>> structures{ { 2,4,23,25 }, { 15,17,49,51 }, { 40,42,57,59 } };
    //int chain_length{ 60 };

    //polymer p0(N_monomers, true);

    //polymer ref_p(chain_length, true);
    //ref_p.set_link_acceptance_prefactor(0.01);
    //p0.set_search_results(structures[1]);
    //determine_prefactor(&ref_p, NMC, 100);
    //correlation_length(ref_p.get_link_acceptance_prefactor(), chain_length, structures[1], NMC);
    //double reference_prefactor{ ref_p.get_link_acceptance_prefactor() };

    //polymer KH_OFF(chain_length, true);
    //KH_OFF.set_link_acceptance_prefactor(reference_prefactor);
    //KH_OFF.set_search_results(structures);
    ////determine_g0(&KH, NMC,100);
    ////std::cout << "past g0" << std::endl;
    //kissing_hairpin_swivel(&KH_OFF, NMC, false);

    //polymer KH_ON(chain_length, true);
    //KH_ON.set_search_results(structures);
    //KH_ON.set_link_acceptance_prefactor(reference_prefactor);
    //std::cout << "past KH with swivel off" << std::endl;
    //kissing_hairpin_swivel(&KH_ON, NMC, true);
    //////////////////////////////////////////////////////////////////////////////////////////
    //varying_helix_length(N_monomers, NMC, "Results/varying_helix_length_sims/");
    //////////////////////////////////////////////////////////////////////////////////////////
    //vary_dangle_sim(N_monomers, NMC, "Results/varying_dangle_sims/");
    

    //////////////////////////////////////////////////////////////////////////////////////////
    //std::cout<<partition_ni_unlinked(3);
    //int max_helix_length(20);
    //std::vector<double> helix_NS_Z_functions(max_helix_length-3);
    //std::string path0{ "Results/rep" };
    //for (int m = 0; m < big_NMC; m++)
    //{
    //    path0 = "Results/rep";
    //    path0 += std::to_string(m) + "_NS_part_functions.txt";
    //    for (int i = 0; i < max_helix_length - 3; i++)
    //    {
    //        helix_NS_Z_functions[i] = partition_ni_unlinked(i + 3, NMC);
    //    }
    //    double_vector_to_txt(path0, helix_NS_Z_functions);

    //    //helix_NS_Z_functions.clear();


    //}
    //////////////////////////////////////////////////////////////////////////////////////////

    //long_dangle_simulation(10, 100, 10, "Results/varying_dangle_sims/");
    //polymer p0(N_monomers, true);
    //std::vector<int> dh{ 5,7,27,29 };

    //}


    //std::vector<double> ideal_chain_pdf_vals(25);
    //double delta{ helix_separation() };
    //ideal_chain_pdf_vals[0] = 1;
    //ideal_chain_pdf_vals[1] = 1 / (2 * pi);
    //for (int i = 0; i < 23; i++) {
    //    ideal_chain_pdf_vals[i + 2] = ideal_chain_pdf(delta, i + 3);
    //}
    //double_vector_to_txt("ideal_chain_pdf_0_to_15.txt", ideal_chain_pdf_vals);

    //code to investigate e2e distance of middle monomer for polymer of length 71. 
    // results outputted to text file.

    //auto start = std::chrono::high_resolution_clock::now();
    //int N_sims{ 50000 };
    ////std::vector<double> sim_results(N_sims);
    //std::vector<double> sim_results(2*N_sims), temp(3);
    //std::vector<double> N_position{ 0,0,30 }, initial_position{ 0,0,0 };
    //int N_seg{ 50 };
    //int N_mid = N_seg / 2;
    //for (int i{ 0 }; i < N_sims; i++) {
    //    std::cout << i << std::endl;
    //    std::vector<std::vector<double>> p{ grow_chain(N_position,initial_position,N_seg) };
    //    temp = vector_subtraction(p[N_mid],initial_position);
    //    sim_results[2*i]=temp[2];
    //    temp = vector_subtraction(N_position, p[N_mid]);
    //    sim_results[2 * i + 1] = temp[2];

    //}
    //double_vector_to_txt("Results/yamakawa_results.txt", sim_results);
    //const char* path = "C:/Users/32471/Documents/Mémoire/testing yamakawa function/sim_results_x.txt";
    //std::ofstream file(path); //open in constructor

    //for (int k{ 0 }; k<sim_results.size(); k++) {
    //    file << sim_results[k] << std::endl;
    //}
    //std::ofstream vector("sim_results2.txt");
    //for (int k{ 0 }; k < N_sims; k++) {
    //    file << sim_results[k] << std::endl;
    //}


}
