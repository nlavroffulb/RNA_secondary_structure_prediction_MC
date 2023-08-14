#include "simulations.h"
#include "math_functions.h"
#include <fstream>
#include <iostream>
// this is the simulation where only one link move is allowed at a time. corresponds to the formation of a hairpin loop.
void hairpin_sim(polymer* p0, int NMC)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    
    while (i < NMC) {
        u_branch = rand2(0, 1);
        //if (p0->get_num_helices() == 2) {
        //    std::cout << "stop " << std::endl;
        //}
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0->link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0->unlink_move(m);
        }

        if (p0->linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }
        i++;
    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;

}

// might need a new strategy for convergence. it does seem correct, just the convergence rate isn't good.
void determine_prefactor(polymer* p0, int NMC, int update_rate, std::string filename)
{
    //polymer p0(N_monomers, true);

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{0}, uc_0{0};

    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };
    p0->set_link_acceptance_prefactor(init_prefactor);
    double r{0};
    std::vector<double> prefactors;
    prefactors.push_back(init_prefactor);
    int cc{ 0 };
    while (i < NMC) {
        // convergence condition.
        if ((i + 1) % update_rate == 0 && i != 0) {
            if (std::abs(std::log(r)) < epsilon) {
                cc++;
                if (cc >= 3) {
                    std::cout << "it converged" << std::endl;
                    prefactors.push_back(0.000001);
                    break;
                }
            }
            else if (i != 0) {
                cc = 0;
            }

        }

        std::cout << "iteration " << i << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0->link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0->unlink_move(m);
        }

        if (p0->linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }

        if ((i+1) % update_rate  == 0 && i!=0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;
          
            p0->set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0->get_link_acceptance_prefactor());
            bc_0 = bound_counts, uc_0 = unbound_counts;
            std::cout << "r = " << r << std::endl;
            std::cout << p0->get_link_acceptance_prefactor() << std::endl;
        }
        i++;
    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;

    double_vector_to_txt("Results/determine_g0.txt",prefactors);



}

double get_prefactor(std::vector<int> structure, int N_monomers, int NMC, int update_rate, std::string filename)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };

    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };

    polymer p0(N_monomers, true);
    p0.set_search_results(structure);
    p0.set_link_acceptance_prefactor(init_prefactor);
    double r{ 0 };
    std::vector<double> prefactors;
    prefactors.push_back(init_prefactor);
    int cc{ 0 };
    while (i < NMC) {
        // convergence condition.
        if ((i + 1) % update_rate == 0 && i != 0) {
            if (std::abs(std::log(r)) < epsilon) {
                cc++;
                if (cc >= 3) {
                    std::cout << "it converged" << std::endl;
                    //prefactors.push_back(0.000001);
                    return prefactors.back();
                    break;
                }
            }
            else if (i != 0) {
                cc = 0;
            }

        }

        std::cout << "iteration " << i << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0.link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0.unlink_move(m);
        }

        if (p0.linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }

        if ((i + 1) % update_rate == 0 && i != 0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;

            p0.set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0.get_link_acceptance_prefactor());
            bc_0 = bound_counts, uc_0 = unbound_counts;
            std::cout << "r = " << r << std::endl;
            std::cout << p0.get_link_acceptance_prefactor() << std::endl;
        }
        i++;
    }

    return prefactors.back();
}

double get_prefactor(double init_prefactor, std::vector<int> structure, int N_monomers, int NMC, int update_rate, std::string filename)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };

    double epsilon{ 2 };// for convergence

    polymer p0(N_monomers, true);
    p0.set_search_results(structure);
    p0.set_link_acceptance_prefactor(init_prefactor);
    double r{ 0 };
    std::vector<double> prefactors, r_values;
    prefactors.push_back(init_prefactor);
    int cc{ 0 };

    while (i < NMC) {
        // convergence condition.
        if ((i + 1) % update_rate == 0 && i != 0) {
            if (std::abs(std::log(r)) < epsilon) {
                cc++;
                if (cc >= 3) {
                    std::cout << "it converged" << std::endl;
                    break;
                }
            }
            else if (i != 0) {
                cc = 0;
            }

        }

        std::cout << "iteration " << i << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0.link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0.unlink_move(m);
        }

        if (p0.linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }

        if ((i + 1) % update_rate == 0 && i != 0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;

            p0.set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0.get_link_acceptance_prefactor());
            bc_0 = bound_counts, uc_0 = unbound_counts;
            std::cout << "r = " << r << std::endl;
            r_values.push_back(r);
            std::cout << p0.get_link_acceptance_prefactor() << std::endl;
        }
        i++;
    }

    if (filename != "") {
        double_vector_to_txt(filename + "prefactors0.txt", prefactors);
        double_vector_to_txt(filename + "r_values0.txt", r_values);
    }

    return prefactors.back();
}

double correlation_length(double acceptance_prefactor, int N_monomers, std::vector<int> structure, int NMC, std::string filename)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };    


    std::vector<int> correlation_lengths;
    bool link_state0{ false }, link_state{false};
    int last_switch{0};


    // a run to calculate the correlation time for a given prefactor. this is the average time that
    // the system takes to switch from the linked to the unlinked configuration.
    polymer p0(N_monomers, true);
    p0.set_link_acceptance_prefactor(acceptance_prefactor);
    p0.set_search_results(structure);
    while (i < NMC) {
        std::cout << "iteration " << i << " out of " << NMC << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0.link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0.unlink_move(m);
        }

        if (p0.linked() == true) {
            if (link_state0 == false) {
                link_state0 = true;
                correlation_lengths.push_back(i - last_switch+1);
                last_switch = i+1;
            }
            bound_counts++;
        }
        else {
            if (link_state0 == true) {
                link_state0 = false;
                correlation_lengths.push_back(i - last_switch+1);
                last_switch = i+1;
            }
            unbound_counts++;
        }
        i++;
    }

    // calculate correlation length
    //int_vector_to_txt(filename + "correlation_lengths" + ".txt", correlation_lengths);

    double avg_correlation{ static_cast<double>(sum_of_elements(correlation_lengths)) / static_cast<double>(correlation_lengths.size()) };
    int update_rate{ static_cast<int>(std::round(avg_correlation)) };

    // we want to make sure that we do a number of trials that is at least 10 times the correlation length 
    // so that there will be at least 10 updates of the prefactor.
    int NMC2{ NMC };
    if (avg_correlation / static_cast<int>(NMC) < 10) {
        NMC2 = 10 * NMC2;
    }
    int bc_0{ 0 }, uc_0{ 0 };
    double epsilon{ 2 };// for convergence
    int cc{ 0 };
    double r{ 0 };
    std::vector<double> prefactors{ acceptance_prefactor };

    i = 0;


    // another run to get a more precise prefactor. compared to the previous run, the re calculation of the 
    // prefactor happens at intervals equal to the calculated correlation time.
    while (i < NMC2) {
        // convergence condition.
        if ((i + 1) % update_rate == 0 && i != 0) {
            if (std::abs(std::log(r)) < epsilon) {
                cc++;
                if (cc >= 3) {
                    std::cout << "it converged" << std::endl;
                    //prefactors.push_back(0.000001);
                    break;
                }
            }
            else if (i != 0) {
                cc = 0;
            }

        }

        std::cout << "iteration " << i << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0.link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;
            bool m;
            p0.unlink_move(m);
        }

        if (p0.linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }

        if ((i + 1) % update_rate == 0 && i != 0) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;

            p0.set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0.get_link_acceptance_prefactor());
            bc_0 = bound_counts, uc_0 = unbound_counts;
            std::cout << "r = " << r << std::endl;
            std::cout << p0.get_link_acceptance_prefactor() << std::endl;
        }
        i++;
    }
    double_vector_to_txt(filename + "prefactors_with_corr_length" + ".txt", prefactors);

    return prefactors.back();

}


void varying_helix_length(int N_monomers, int NMC, std::string path)
{
    int i{ 0 },l{0};
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    int helix_max_length{ 20 };

    int N_0{N_monomers};


    bool skip_trial{ false };
    std::vector<double> old_weights, new_weights;
    std::vector<double> linker_weights(helix_max_length-3), helix_weights(helix_max_length - 3), 
        linker_weights_l, helix_weights_l, physical_rejections(helix_max_length - 3),
            old_linker_weights, new_linker_weights, link_probability(helix_max_length-3);


    std::vector<double> link_probabilities(helix_max_length - 3), prefactors(helix_max_length-3);
    double link_prob;
    //p0->set_link_acceptance_prefactor(1);

    std::vector<int> helix;
    int N;

    while (l < helix_max_length-3) {
        if (l == 14-3) {
            std::cout << "stop" << std::endl;
        }


        std::cout << "considering helices of length" << 3 + l << std::endl;
        helix = { 0,2 + l,N_0 + 2 * l - 3 - l,N_0 + 2 * l - 1 };
        N = N_0 + 2 * l;

        polymer p0(N, true);
        p0.set_search_results(helix);

        prefactors[l] = get_prefactor(helix, N, NMC, 100);
        prefactors[l] = correlation_length(prefactors[l], N, helix, NMC, std::to_string(l));
        p0.set_link_acceptance_prefactor(prefactors[l]);

        i = 0;
        bound_counts = 0;
        unbound_counts = 0;


        double overstretch_rejections{ 0 };
        std::vector<double> link_acceptances, unlink_acceptances;


        helix_weights_l.clear();
        old_linker_weights.clear();
        new_linker_weights.clear();
        old_weights.clear();
        new_weights.clear();
        link_acceptances.clear();
        unlink_acceptances.clear();
        overstretch_rejections = 0;

        double lacc, uacc;
        while (i < NMC) {

            skip_trial = false;
            std::cout << i << "out of " << NMC << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0.link_move(skip_trial);
                if (skip_trial == true) {
                    if (!p0.linked()) {
                        overstretch_rejections = overstretch_rejections + 1.0;
                    }
                }
                else {
                    lacc = p0.link_acceptance(1);
                    if (lacc != 0.0) {
                        link_acceptances.push_back(lacc);
                    }
                }
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0.unlink_move(skip_trial);
                uacc = p0.link_acceptance(0);
                if (skip_trial == false && uacc != 0.0) {
                    unlink_acceptances.push_back(p0.link_acceptance(0));
                }
            }

            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            if (skip_trial == false) {
                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

                if (old_w != 0 && new_w != 0) {

                    old_weights.push_back(p0.get_old_config_weight());
                    new_weights.push_back(p0.get_new_config_weight());

                    helix_weights_l.push_back(p0.get_helix_weight());

                    old_linker_weights.push_back(p0.get_hairpin_weight(false));
                    new_linker_weights.push_back(p0.get_hairpin_weight(true));

                }
            }

            i++;
        }
        helix_weights[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        link_probability[l] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
        physical_rejections[l] = overstretch_rejections;

        double_vector_to_txt(path + "link_acc_extra_bp=" + std::to_string(l) + ".txt", link_acceptances);
        double_vector_to_txt(path + "unlink_acc_extra_bp=" + std::to_string(l) + ".txt", unlink_acceptances);
        double_vector_to_txt(path + "old_linker_weights_extra_bp=" + std::to_string(l) + ".txt", old_linker_weights);
        double_vector_to_txt(path + "new_linker_weights_extra_bp=" + std::to_string(l) + ".txt", new_linker_weights);
        double_vector_to_txt(path + "old_weights_extra_bp=" + std::to_string(l) + ".txt", old_weights);
        double_vector_to_txt(path + "new_weights_extra_bp=" + std::to_string(l) + ".txt", new_weights);

        l++;


    }

    std::cout << "Unlink state: " << unbound_counts << std::endl;
    std::cout << "Link state: " << bound_counts << std::endl;

    double_vector_to_txt(path + "prefactors.txt", prefactors);
    double_vector_to_txt(path + "link_probabilities.txt", link_probability);
    double_vector_to_txt(path + "overstretch_rejections.txt", physical_rejections);
    double_vector_to_txt(path + "helix_weights.txt", helix_weights);



    print_1d_doub_vec(link_probabilities);
}



void vary_dangle_sim(int N_monomers, int NMC, std::string path)
{
    int max_dangle_length{ 20 };

    //data 
    std::vector<double> link_probability(max_dangle_length), link_acceptances, unlink_acceptances, prefactors(max_dangle_length);

    std::vector<double> linker_weights(max_dangle_length), helix_weights(max_dangle_length), 
        linker_weights_l, helix_weights_l, physical_rejections(max_dangle_length),
        old_linker_weights, new_linker_weights, old_dangle_weights, new_dangle_weights;
    double overstretch_rejections{ 0 };

    bool skip_trial{ false };


    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    std::vector<double> old_weights, new_weights;
    //p0->force_unbound_state();

    int N0{ N_monomers }, N{N0};
    std::vector<int> helix(4);

    for (int l{ 0 }; l < max_dangle_length; l++) {

        // the tail is going to get longer such that the total number of monomers in the polymer increases each time.
        N = N0 + 2*l;
        polymer p0(N, true);
        helix = { 0 + l, 2 + l, N - 3-l, N - 1-l };
        p0.set_search_results(helix);



        prefactors[l] = get_prefactor(helix, N, NMC, 100);
        prefactors[l] = correlation_length(prefactors[l], N, helix, NMC, std::to_string(l));
        p0.set_link_acceptance_prefactor(prefactors[l]);

        std::cout << "dangle length " << l << std::endl;

        i = 0;
        bound_counts = 0;
        unbound_counts = 0;

        helix_weights_l.clear();
        linker_weights_l.clear();
        old_weights.clear();
        new_weights.clear();
        link_acceptances.clear();
        unlink_acceptances.clear();
        overstretch_rejections = 0;
        double lacc, uacc;
        while (i < NMC) {
            skip_trial = false;
            std::cout << i << "out of " << NMC << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0.link_move(skip_trial);
                if (skip_trial == true) {
                    if (!p0.linked()) {
                        overstretch_rejections = overstretch_rejections + 1.0;
                    }
                }
                else {
                    lacc = p0.link_acceptance(1);
                    if (lacc != 0.0) {
                        link_acceptances.push_back(lacc);
                    }
                }
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0.unlink_move(skip_trial);
                uacc = p0.link_acceptance(0);
                if (skip_trial == false && uacc != 0.0) {
                    unlink_acceptances.push_back(p0.link_acceptance(0));
                }
            }

            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            if (skip_trial == false) {
                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

                if (old_w != 0 && new_w != 0) {

                    old_weights.push_back(p0.get_old_config_weight());
                    new_weights.push_back(p0.get_new_config_weight());

                    linker_weights_l.push_back(p0.get_hairpin_weight(true));
                    helix_weights_l.push_back(p0.get_helix_weight());

                    old_linker_weights.push_back(p0.get_hairpin_weight(false));
                    new_linker_weights.push_back(p0.get_hairpin_weight(true));

                    if (l != 0) {
                        old_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, false));
                        new_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, true));
                    }


                }
            }

            i++;
        }

        linker_weights[l] = sum_of_elements(linker_weights_l) / static_cast<double>(linker_weights_l.size());
        helix_weights[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        link_probability[l] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
        physical_rejections[l] = overstretch_rejections;


        double_vector_to_txt(path + "link_acc" + std::to_string(l) + ".txt", link_acceptances);
        double_vector_to_txt(path + "unlink_acc" + std::to_string(l) + ".txt", unlink_acceptances);


        double_vector_to_txt(path + "old_linker_weights_dangle=" + std::to_string(l) + ".txt", old_linker_weights);
        double_vector_to_txt(path + "new_linker_weights_dangle=" + std::to_string(l) + ".txt", new_linker_weights);

        double_vector_to_txt(path + "old_dangle_weights_dangle=" + std::to_string(l) + ".txt", old_dangle_weights);
        double_vector_to_txt(path + "new_dangle_weights_dangle=" + std::to_string(l) + ".txt", new_dangle_weights);

        double_vector_to_txt(path + "old_weights_dangle=" + std::to_string(l) + ".txt", old_weights);
        double_vector_to_txt(path + "new_weights_dangle=" + std::to_string(l) + ".txt", new_weights);

    }
    double_vector_to_txt(path + "prefactors.txt", prefactors);
    double_vector_to_txt(path + "helix_weights.txt", helix_weights);
    double_vector_to_txt(path + "linker_weights.txt", linker_weights);
    double_vector_to_txt(path + "overstretch_rejections.txt", physical_rejections);
    double_vector_to_txt(path + "link_probabilities.txt", link_probability);

}


void kissing_hairpin_swivel(polymer* p0, int NMC, bool swivel_allowed, std::string filename)
{
    double swivel_switch;
    swivel_allowed == true ? swivel_switch = 1 : swivel_switch = 2.0 / 3.0;
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };

    int kissing_hairpin_count{ 0 }, two_link_state{0}, one_link_state{0}, free_state{0};

    p0->force_unbound_state();

    while (i < NMC) {
        std::cout << i << std::endl;
        u_branch = rand2(0, swivel_switch);
        if (u_branch < 1.0/3.0) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0->link_move(m);
        }
        else if(u_branch <2.0/3.0) {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0->unlink_move(m);
        }
        else if (u_branch < 1) {
            p0->swivel_move();
        }

        if (p0->get_num_helices() == 0) {
            free_state++;
        }
        else if (p0->get_num_helices() == 1) {
            one_link_state++;
        }
        else if (p0->get_num_helices() == 2) {
            two_link_state++;
        }
        else if (p0->get_num_helices() == 3) {
            kissing_hairpin_count++;
        }
        i++;
    }
    std::cout << "out of main loop" << std::endl; 

    double P_kissing_hairpin{ static_cast<double>(kissing_hairpin_count) / static_cast<double>(NMC) }, 
    P_two_links{ static_cast<double>(two_link_state) / static_cast<double>(NMC) },
        P_one_link{ static_cast<double>(one_link_state) / static_cast<double>(NMC) },
        P_free{ static_cast<double>(free_state) / static_cast<double>(NMC) };
    std::string file_name;
    if (swivel_allowed == true) {
        file_name = "P_kissing_hairpin_ON.txt";
    }
    else {
        file_name = "P_kissing_hairpin_OFF.txt";
    }
    double_vector_to_txt("Results/" + file_name, {P_free, P_one_link, P_two_links,P_kissing_hairpin});

}

void old_new_weights_sim(polymer* p0, int NMC, std::string filename)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };


    std::vector<double> old_config_weights, new_config_weights;
    bool skip_trial{false};

    p0->force_unbound_state();

    while (i < NMC) {
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            p0->link_move(skip_trial);
        }
        else {

            p0->unlink_move(skip_trial);
        }

        if (p0->linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }
        if (skip_trial == false) {
            double old_w{ p0->get_old_config_weight() }, new_w{ p0->get_new_config_weight() };
            
            if (old_w != 0 && new_w != 0) {
                old_config_weights.push_back(p0->get_old_config_weight());
                new_config_weights.push_back(p0->get_new_config_weight());
            }
        }
        i++;
    }

    double_vector_to_txt("Results/old_vs_new_old_weights.txt", old_config_weights);
    double_vector_to_txt("Results/old_vs_new_new_weights.txt", new_config_weights);



}

void vary_WCA_potential(double prefactor, int N_monomers, int NMC, const std::vector<int>& search_result, std::string path)
{
    int i{ 0 }, k{ 0 };
    double n_runs{ 10.0 };
    std::vector<double> epsilons{ generate_linear_array(0,2,10) };

    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    bool skip_trial{ false };
    std::vector<double> physical_rejections(n_runs);
    int overstretch_rejections{ 0 };

    std::vector<double> old_weights, new_weights, helix_weights_l, linker_weights_l;
    std::vector<double> helix_weights(n_runs), linker_weights(n_runs), link_probability(n_runs), avg_linker_weight(n_runs);

    double current_epsilon;
    while (k < n_runs) {
        polymer p0(N_monomers, true);
        p0.set_link_acceptance_prefactor(prefactor);
        p0.set_search_results(search_result);

        current_epsilon = epsilons[k];
        set_WCA_parameter(current_epsilon);


        // reset temp data for new round
        i = 0;
        bound_counts = 0;
        unbound_counts = 0;
        linker_weights_l.clear();
        helix_weights_l.clear();
        old_weights.clear();
        new_weights.clear();
        overstretch_rejections = 0;

        while (i < NMC) {
            std::cout << i << "out of " << NMC << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                std::cout << "link branch" << std::endl;
                p0.link_move(skip_trial);
                if (skip_trial == true) {
                    overstretch_rejections = overstretch_rejections + 1.0;
                }
            }
            else {
                std::cout << "unlink branch" << std::endl;

                p0.unlink_move(skip_trial);
            }

            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            if (skip_trial == false) {
                double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

                if (old_w != 0 && new_w != 0) {
                    old_weights.push_back(p0.get_old_config_weight());
                    new_weights.push_back(p0.get_new_config_weight());

                    linker_weights_l.push_back(p0.get_hairpin_weight(true));
                    helix_weights_l.push_back(p0.get_helix_weight());
                }
            }

            i++;
        }
        linker_weights[k] = sum_of_elements(linker_weights_l) / static_cast<double>(linker_weights_l.size());
        helix_weights[k] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        link_probability[k] = static_cast<double>(bound_counts) / (static_cast<double>(NMC));
        physical_rejections[k] = overstretch_rejections;

        double_vector_to_txt(path + "old_weightsEPS=" + std::to_string(current_epsilon) + ".txt", old_weights);
        double_vector_to_txt(path + "new_weightsEPS=" + std::to_string(current_epsilon) + ".txt", new_weights);


        k++;
    }
    double_vector_to_txt(path + "helix_weights.txt", helix_weights);
    double_vector_to_txt(path + "linker_weights.txt", linker_weights);
    double_vector_to_txt(path + "overstretch_rejections.txt", physical_rejections);
    double_vector_to_txt(path + "link_probabilities.txt", link_probability);

}

void hairpin_biased_equilibrium(std::vector<int> structure, int N_monomers, int NMC, int update_rate, std::string filename)
{


    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };
    //int update_rate{ 500 };
    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };

    polymer p0(N_monomers, true);
    p0.set_search_results(structure);
    p0.set_link_acceptance_prefactor(init_prefactor);
    double r{ 0 };

    // convergence data
    std::vector<double> prefactors, r_values;
    prefactors.push_back(init_prefactor);

    // other data
    std::vector<double> old_weights, new_weights, linker_weights_l, helix_weights_l, old_linker_weights, new_linker_weights;

    int cc{ 0 };

    bool update_r{ true };

    while (i < NMC && update_r == true) {
        // convergence condition.
        if ((i + 1) % update_rate == 0 && i != 0) {
            if (std::abs(std::log(r)) < epsilon) {
                cc++;
                if (cc >= 3) {
                    std::cout << "it converged" << std::endl;
                    update_r == false;
                }
            }
            else if (i != 0) {
                cc = 0;
            }

        }

        std::cout << "iteration " << i << std::endl;
        u_branch = rand2(0, 1);
        if (u_branch < 0.5) {// link_branch
            //std::cout << "link branch" << std::endl;
            bool m;
            p0.link_move(m);
        }
        else {
            //std::cout << "unlink branch" << std::endl;

            bool m;
            p0.unlink_move(m);
        }

        // get data
        double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

        if (old_w != 0 && new_w != 0) {

            old_weights.push_back(p0.get_old_config_weight());
            new_weights.push_back(p0.get_new_config_weight());

            linker_weights_l.push_back(p0.get_hairpin_weight(true));
            helix_weights_l.push_back(p0.get_helix_weight());

            old_linker_weights.push_back(p0.get_hairpin_weight(false));
            new_linker_weights.push_back(p0.get_hairpin_weight(true));

        }
        // histogram of bound and unbound counts
        if (p0.linked() == true) {
            bound_counts++;
        }
        else {
            unbound_counts++;
        }
        // update the prefactor in the acceptance rule
        if ((i + 1) % update_rate == 0 && i != 0 && update_r == true) {

            r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
            r = 1 / r;
            p0.set_link_acceptance_prefactor(0, r);
            prefactors.push_back(p0.get_link_acceptance_prefactor());
            r_values.push_back(r);
            bc_0 = bound_counts, uc_0 = unbound_counts;
            std::cout << "r = " << r << std::endl;
            std::cout << p0.get_link_acceptance_prefactor() << std::endl;
        }

        i++;
    }
    double_vector_to_txt(filename + "prefactors.txt", prefactors);
    double_vector_to_txt(filename + "r_values.txt", r_values);

    double_vector_to_txt(filename + "helix_weights.txt", helix_weights_l);
    double_vector_to_txt(filename + "old_linker_weights.txt", old_linker_weights);
    double_vector_to_txt(filename + "new_linker_weights.txt", new_linker_weights);
    double_vector_to_txt(filename + "old_weights.txt", old_weights);
    double_vector_to_txt(filename + "new_weights.txt", new_weights);



}

void WCA_sim(int N_monomers, int NMC, int init_update_rate, std::string path)
{

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };
    //int update_rate{ 500 };
    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };

    double r{ 0 };

    int n_values{ 9 };
    std::vector<double> epsilons{ generate_linear_array(0,2,n_values) };
    

    // convergence data
    std::vector<double> prefactors, r_values;
    prefactors.push_back(init_prefactor);

    // other data
    std::vector<double> old_weights, new_weights, linker_weights_l, helix_weights_l, old_linker_weights, new_linker_weights,
        old_dangle_weights, new_dangle_weights;
    // averages data
    std::vector<double> avg_helix_weight(n_values), avg_prefactors(n_values), avg_old_weights(n_values), avg_new_weights(n_values);

    int N0{ N_monomers }, N{ N0 };
    std::vector<int> helix(4);
    int cc{ 0 };
    bool update_r{ true };
    double current_epsilon;

    for (int l{ 0 }; l < epsilons.size(); l++) {

        polymer p0(N, true);
        helix = { 0, 2 , N - 3, N - 1};
        p0.set_search_results(helix);
        p0.set_link_acceptance_prefactor(init_prefactor);
        current_epsilon = epsilons[l];
        set_WCA_parameter(current_epsilon);

        // resetting data 
        old_weights.clear();
        new_weights.clear();
        old_linker_weights.clear(), new_linker_weights.clear();
        helix_weights_l.clear();
        prefactors.clear();
        r_values.clear();
        cc = 0;
        i = 0;
        update_r = true;
        bound_counts = 0, unbound_counts = 0;
        bc_0 = 0, uc_0 = 0;
        r = 0;

        while (update_r == true) {

            //while (update_r == true) {
                // convergence condition.
            if (i % init_update_rate == 0 && i != 0) {
                if (std::abs(std::log(r)) < epsilon) {
                    cc++;
                    if (cc >= 3) {
                        std::cout << "it converged" << std::endl;
                        update_r = false;
                    }
                }
                else if (i != 0) {
                    cc = 0;
                }

            }

            std::cout << "iteration " << i << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                //std::cout << "link branch" << std::endl;
                bool m;
                p0.link_move(m);
            }
            else {
                //std::cout << "unlink branch" << std::endl;

                bool m;
                p0.unlink_move(m);
            }

            // get data
            double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

            //if (old_w != 0 && new_w != 0) {

            //    old_weights.push_back(p0.get_old_config_weight());
            //    new_weights.push_back(p0.get_new_config_weight());

            //    linker_weights_l.push_back(p0.get_hairpin_weight(true));
            //    helix_weights_l.push_back(p0.get_helix_weight());

            //    old_linker_weights.push_back(p0.get_hairpin_weight(false));
            //    new_linker_weights.push_back(p0.get_hairpin_weight(true));

            //    if (l != 0) {
            //        old_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, false));
            //        new_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, true));
            //    }


            //}
            // histogram of bound and unbound counts
            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            // update the prefactor in the acceptance rule
            if ((i + 1) % init_update_rate == 0 && i != 0 && update_r == true) {

                r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
                r = 1 / r;
                p0.set_link_acceptance_prefactor(0, r);
                prefactors.push_back(p0.get_link_acceptance_prefactor());
                r_values.push_back(r);
                bc_0 = bound_counts, uc_0 = unbound_counts;
                std::cout << "r = " << r << std::endl;
                std::cout << p0.get_link_acceptance_prefactor() << std::endl;
            }

            i++;
        }
        //double_vector_to_txt(path + std::to_string(current_epsilon) + "_old_weights.txt", old_weights);
        //double_vector_to_txt(path + std::to_string(current_epsilon) + "_new_weights.txt", new_weights);

        //avg_helix_weight[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        avg_prefactors[l] = sum_of_elements(prefactors) / static_cast<double>(prefactors.size());


    }
    double_vector_to_txt(path + "WCA_avg_prefactors.txt", avg_prefactors);
    //double_vector_to_txt(path + "WCA_avg_old_weights.txt", avg_old_weights);
    //double_vector_to_txt(path + "WCA_avg_new_weights.txt", avg_new_weights);
    //double_vector_to_txt(path + "dangle_avg_helix_weights.txt", avg_helix_weight);


}

//void swivel_sim(int N_monomers, int NMC, int init_update_rate, std::string path)
//{
//    //int repeats{ 1000 };
//    int success_rate{0}, phys_rejections{0};
//    std::vector<double> epsilons{ generate_log_array(0.00000001,10,10) };
//    double current_epsilon;
//
//    std::vector<std::vector<int>> always_linked_structures{ { 0,2,8,10 }, { 12,14,20,22 } };
//    std::vector<int> KH{ 4,6,16,18 };
//    int N = 24;
//
//    bool overstr{ false };
//    int i{ 0 };
//
//
//    double u_branch;
//    int bound_counts{ 0 }, unbound_counts{ 0 };
//    int bc_0{ 0 }, uc_0{ 0 };
//    double init_prefactor{ 1 };
//    double r{ 0 };
//    int cc{ 0 };
//    bool update_r{ true };
//
//    std::vector<double> prefactors;
//    std::vector<double>  avg_prefactors(epsilons.size()), physical_rejec_rate(epsilons.size());
//
//    //int update_rate{ 500 };
//    double tau{ 2 };// for convergence
//
//    std::vector<double> ON_physical_rejections(epsilons.size()), OFF_physical_rejections(epsilons.size());
//    std::vector<double> ON_P_link(epsilons.size()), OFF_P_link(epsilons.size());
//
//    // SWIVEL OFF
//    for (int l = 0; l < epsilons.size(); l++)
//    {
//        current_epsilon = epsilons[l];
//        set_WCA_parameter(current_epsilon);
//
//        // INITIALIZE STATE WITH TWO LINKERS
//        bool overstr = false;
//        phys_rejections = 0;
//        polymer p0(N, true);
//        p0.set_search_results(always_linked_structures[0]);
//        p0.always_accept_link_move(overstr);
//
//        p0.set_search_results(always_linked_structures[1]);
//        overstr = false;
//        while (i == 0 || overstr == true) {
//            p0.always_accept_link_move(overstr);
//            i++;
//        }
//        overstr = false;
//        update_r = true;
//
//        p0.set_link_acceptance_prefactor(1);
//        p0.set_search_results(KH);
//        prefactors.clear();
//        cc = 0;
//        i = 0;
//        update_r = true;
//        bound_counts = 0, unbound_counts = 0;
//        bc_0 = 0, uc_0 = 0;
//        r = 0;
//
//        // INITIALIZATION COMPLETE
//        // ATTEMPT KH
//        // SIMULATION LOOP. GO TO CONVERGENCE
//        while (update_r == true) {
//            overstr = false;
//            // convergence condition.
//            if (i % init_update_rate == 0 && i != 0 && update_r == true) {
//                if (std::abs(std::log(r)) < tau) {
//                    cc++;
//                    if (cc >= 3) {
//                        std::cout << "it converged" << std::endl;
//                        update_r = false;
//                    }
//                }
//                else if (i != 0) {
//                    cc = 0;
//                }
//
//            }
//
//            std::cout << "iteration " << i << std::endl;
//            u_branch = rand2(0, 1);
//            if (u_branch < 0.5) {// link_branch
//                //std::cout << "link branch" << std::endl;
//                p0.always_accept_link_move(overstr);
//                if (p0.get_num_helices()!=3 && overstr == true) {
//                    phys_rejections++;
//                }
//            }
//            else {
//
//                p0.unlink_move(overstr);
//            }
//
//
//
//            // histogram of bound and unbound counts
//            if (p0.get_num_helices() == 3) {
//                bound_counts++;
//            }
//            else {
//                unbound_counts++;
//            }
//            // update the prefactor in the acceptance rule
//            if ((i + 1) % init_update_rate == 0 && i != 0 && update_r == true) {
//
//                r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
//                r = 1 / r;
//                p0.set_link_acceptance_prefactor(0, r);
//                prefactors.push_back(p0.get_link_acceptance_prefactor());
//                bc_0 = bound_counts, uc_0 = unbound_counts;
//                std::cout << "r = " << r << std::endl;
//                std::cout << p0.get_link_acceptance_prefactor() << std::endl;
//            }
//
//            i++;
//        }
//        avg_prefactors[l] = sum_of_elements(prefactors) / static_cast<double>(prefactors.size());
//        physical_rejec_rate[l] = static_cast<double>(phys_rejections) / static_cast<double>(i);
//
//    }
//    double_vector_to_txt(path + "_OFF_avg_prefactors.txt", avg_prefactors);
//    double_vector_to_txt(path + "_OFF_physical_rejec_rate.txt", physical_rejec_rate);
//
//    // SWIVEL ON
//    for (int l = 0; l < epsilons.size(); l++)
//    {
//        current_epsilon = epsilons[l];
//        set_WCA_parameter(current_epsilon);
//
//        // INITIALIZE STATE WITH TWO LINKERS
//        //success_rate = 0;
//        //fail_rate = 0;
//        polymer p0(N, true);
//        p0.set_search_results(always_linked_structures[0]);
//        p0.always_accept_link_move(overstr);
//
//        p0.set_search_results(always_linked_structures[1]);
//        overstr = false;
//        while (i == 0 || overstr == true) {
//            p0.always_accept_link_move(overstr);
//            i++;
//        }
//        overstr = false;
//        update_r = true;
//
//        p0.set_search_results(KH);
//
//        bool overstr = false;
//        phys_rejections = 0;
//        prefactors.clear();
//        cc = 0;
//        i = 0;
//        update_r = true;
//        bound_counts = 0, unbound_counts = 0;
//        bc_0 = 0, uc_0 = 0;
//        r = 0;
//
//
//        // INITIALIZATION COMPLETE
//        // ATTEMPT KH
//        // SIMULATION LOOP. GO TO CONVERGENCE
//        while (update_r == true) {
//            overstr = false;
//            // convergence condition.
//            if (i % init_update_rate == 0 && i != 0 && update_r == true) {
//                if (std::abs(std::log(r)) < tau) {
//                    cc++;
//                    if (cc >= 3) {
//                        std::cout << "it converged" << std::endl;
//                        update_r = false;
//                    }
//                }
//                else if (i != 0) {
//                    cc = 0;
//                }
//
//            }
//
//            std::cout << "iteration " << i << std::endl;
//            u_branch = rand2(0, 1);
//            if (u_branch < 1.0/3.0) {// link_branch
//                //std::cout << "link branch" << std::endl;
//                p0.always_accept_link_move(overstr);
//                if (!p0.linked() && overstr == true) {
//                    phys_rejections++;
//                }
//            }
//            else if(u_branch< 2.0/3.0) {
//
//                p0.unlink_move(overstr);
//            }
//            else {
//                p0.swivel_move();
//            }
//
//
//
//            // histogram of bound and unbound counts
//            if (p0.linked() == true) {
//                bound_counts++;
//            }
//            else {
//                unbound_counts++;
//            }
//            // update the prefactor in the acceptance rule
//            if ((i + 1) % init_update_rate == 0 && i != 0 && update_r == true) {
//
//                r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
//                r = 1 / r;
//                p0.set_link_acceptance_prefactor(0, r);
//                prefactors.push_back(p0.get_link_acceptance_prefactor());
//                bc_0 = bound_counts, uc_0 = unbound_counts;
//                std::cout << "r = " << r << std::endl;
//                std::cout << p0.get_link_acceptance_prefactor() << std::endl;
//            }
//
//            i++;
//        }
//        avg_prefactors[l] = sum_of_elements(prefactors) / static_cast<double>(prefactors.size());
//        physical_rejec_rate[l] = static_cast<double>(phys_rejections) / static_cast<double>(i);
//
//    }
//
//    double_vector_to_txt(path + "_ON_avg_prefactors.txt", avg_prefactors);
//    double_vector_to_txt(path + "_ON_physical_rejec_rate.txt", physical_rejec_rate);
//
//}
void swivel_sim(int N_monomers, int NMC, int init_update_rate, std::string path)
{
    //int repeats{ 1000 };
    int success_rate{ 0 }, phys_rejections{ 0 };
    std::vector<double> epsilons{ generate_log_array(0.00000001,10,10) };
    double current_epsilon;

    std::vector<std::vector<int>> always_linked_structures{ { 0,2,8,10 }, { 12,14,20,22 } };
    std::vector<int> KH{ 4,6,16,18 };
    int N = 24;

    bool overstr{ false };
    int i{ 0 };


    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };
    double init_prefactor{ 1 };
    double r{ 0 };
    int cc{ 0 };
    bool update_r{ true };

    std::vector<double> prefactors;
    std::vector<double>  avg_prefactors(epsilons.size()), physical_rejec_rate(epsilons.size());

    //int update_rate{ 500 };
    double tau{ 2 };// for convergence

    std::vector<double> ON_physical_rejections_time(epsilons.size()), OFF_physical_rejection_time(epsilons.size());
    std::vector<double> ON_P_link(epsilons.size()), OFF_P_link(epsilons.size());

    bool quit_loop{ false };

    // HAVE TO DECIDE WHERE TO PUT THE MAIN BIG LOOP. ALREADY IN THIS LOOP EVERYTHING IS INDEPENDENT. WE INITIALIZE 
    // INDEPENDENTLY EACH TIME. SO WE COULD ALREADY HAVE THE AVERAGE HERE. ITS JUST A BIT RISKIER IF WE WANT TO SEE 
    // WHATS GOING ON. 
    // 
    // SWIVEL OFF
    for (int l = 0; l < epsilons.size(); l++)
    {
        current_epsilon = epsilons[l];
        set_WCA_parameter(current_epsilon);
        quit_loop = false;
        for (size_t i = 0; i < NMC; i++)
        {

        }
        while (!quit_loop) {
            // INITIALIZE STATE WITH TWO LINKERS
            bool overstr = false;
            phys_rejections = 0;
            polymer p0(N, true);
            p0.set_search_results(always_linked_structures[0]);
            p0.always_accept_link_move(overstr);

            p0.set_search_results(always_linked_structures[1]);
            overstr = false;
            while (i == 0 || overstr == true) {
                p0.always_accept_link_move(overstr);
                i++;
            }
            overstr = false;

            p0.set_link_acceptance_prefactor(1);
            p0.set_search_results(KH);
            prefactors.clear();

            // INITIALIZATION COMPLETE
            // ATTEMPT KH
            // SIMULATION LOOP. GO TO CONVERGENCE

            overstr = false;
            std::cout << "iteration " << i << std::endl;
            p0.always_accept_link_move(overstr);
            if (p0.get_num_helices() != 3 && overstr == true) {
                phys_rejections++;
            }

            i++;
            if (p0.get_num_helices() == 3) {
                quit_loop = true;
            }
        }
        //avg_prefactors[l] = sum_of_elements(prefactors) / static_cast<double>(prefactors.size());
        //physical_rejec_rate[l] = static_cast<double>(phys_rejections) / static_cast<double>(i);

    }
    double_vector_to_txt(path + "_OFF_avg_prefactors.txt", avg_prefactors);
    double_vector_to_txt(path + "_OFF_physical_rejec_rate.txt", physical_rejec_rate);

    // SWIVEL ON
    for (int l = 0; l < epsilons.size(); l++)
    {
        current_epsilon = epsilons[l];
        set_WCA_parameter(current_epsilon);

        // INITIALIZE STATE WITH TWO LINKERS
        //success_rate = 0;
        //fail_rate = 0;
        polymer p0(N, true);
        p0.set_search_results(always_linked_structures[0]);
        p0.always_accept_link_move(overstr);

        p0.set_search_results(always_linked_structures[1]);
        overstr = false;
        while (i == 0 || overstr == true) {
            p0.always_accept_link_move(overstr);
            i++;
        }
        overstr = false;
        update_r = true;

        p0.set_search_results(KH);

        bool overstr = false;
        phys_rejections = 0;
        prefactors.clear();
        cc = 0;
        i = 0;
        update_r = true;
        bound_counts = 0, unbound_counts = 0;
        bc_0 = 0, uc_0 = 0;
        r = 0;


        // INITIALIZATION COMPLETE
        // ATTEMPT KH
        // SIMULATION LOOP. GO TO CONVERGENCE
        while (update_r == true) {
            overstr = false;
            // convergence condition.
            if (i % init_update_rate == 0 && i != 0 && update_r == true) {
                if (std::abs(std::log(r)) < tau) {
                    cc++;
                    if (cc >= 3) {
                        std::cout << "it converged" << std::endl;
                        update_r = false;
                    }
                }
                else if (i != 0) {
                    cc = 0;
                }

            }

            std::cout << "iteration " << i << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 1.0 / 3.0) {// link_branch
                //std::cout << "link branch" << std::endl;
                p0.always_accept_link_move(overstr);
                if (!p0.linked() && overstr == true) {
                    phys_rejections++;
                }
            }
            else if (u_branch < 2.0 / 3.0) {

                p0.unlink_move(overstr);
            }
            else {
                p0.swivel_move();
            }



            // histogram of bound and unbound counts
            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            // update the prefactor in the acceptance rule
            if ((i + 1) % init_update_rate == 0 && i != 0 && update_r == true) {

                r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
                r = 1 / r;
                p0.set_link_acceptance_prefactor(0, r);
                prefactors.push_back(p0.get_link_acceptance_prefactor());
                bc_0 = bound_counts, uc_0 = unbound_counts;
                std::cout << "r = " << r << std::endl;
                std::cout << p0.get_link_acceptance_prefactor() << std::endl;
            }

            i++;
        }
        avg_prefactors[l] = sum_of_elements(prefactors) / static_cast<double>(prefactors.size());
        physical_rejec_rate[l] = static_cast<double>(phys_rejections) / static_cast<double>(i);

    }

    double_vector_to_txt(path + "_ON_avg_prefactors.txt", avg_prefactors);
    double_vector_to_txt(path + "_ON_physical_rejec_rate.txt", physical_rejec_rate);

}

void long_dangle_simulation(int N_init, int NMC, int update_rate, std::string path)
{

    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };
    //int update_rate{ 500 };
    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };
    int max_dangle_length{ 1 };

    double r{ 0 };

    // convergence data
    std::vector<double> prefactors, r_values;
    prefactors.push_back(init_prefactor);

    // other data
    std::vector<double> old_weights, new_weights, linker_weights_l, helix_weights_l, old_linker_weights, new_linker_weights,
        old_dangle_weights,new_dangle_weights;
    // averages data
    std::vector<double> avg_helix_weight(max_dangle_length), avg_prefactors(max_dangle_length);

    int N0{ N_init }, N{ N0 };
    std::vector<int> helix(4);
    int cc{ 0 };
    bool update_r{ true };

    for (int l{ 0 }; l < max_dangle_length; l++) {

        // the tail is going to get longer such that the total number of monomers in the polymer increases each time.
        N = N0 + 2 * l;
        polymer p0(N, true);
        helix = { 0 + l, 2 + l, N - 3 - l, N - 1 - l };
        p0.set_search_results(helix);
        p0.set_link_acceptance_prefactor(init_prefactor);


        // resetting data 
        old_weights.clear();
        new_weights.clear();
        old_linker_weights.clear(), new_linker_weights.clear();
        helix_weights_l.clear();
        prefactors.clear();
        r_values.clear();
        cc = 0;
        i = 0;
        update_r = true;
        bound_counts = 0, unbound_counts = 0;
        bc_0 = 0, uc_0 = 0;
        r = 0;

        while (i < NMC) {

        //while (update_r == true) {
            // convergence condition.
            if (i % update_rate == 0 && i != 0) {
                if (std::abs(std::log(r)) < epsilon) {
                    cc++;
                    if (cc >= 3) {
                        std::cout << "it converged" << std::endl;
                        update_r = false;
                    }
                }
                else if (i != 0) {
                    cc = 0;
                }

            }

            std::cout << "iteration " << i << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                //std::cout << "link branch" << std::endl;
                bool m;
                p0.link_move(m);
            }
            else {
                //std::cout << "unlink branch" << std::endl;

                bool m;
                p0.unlink_move(m);
            }

            // get data
            double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

            if (old_w != 0 && new_w != 0) {

                old_weights.push_back(p0.get_old_config_weight());
                new_weights.push_back(p0.get_new_config_weight());

                linker_weights_l.push_back(p0.get_hairpin_weight(true));
                helix_weights_l.push_back(p0.get_helix_weight());

                old_linker_weights.push_back(p0.get_hairpin_weight(false));
                new_linker_weights.push_back(p0.get_hairpin_weight(true));

                if (l != 0) {
                    old_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, false));
                    new_dangle_weights.push_back(p0.get_subsection_weight({ 0,l }, true));
                }


            }
            // histogram of bound and unbound counts
            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            // update the prefactor in the acceptance rule
            if ((i + 1) % update_rate == 0 && i != 0 && update_r == true) {

                r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
                r = 1 / r;
                p0.set_link_acceptance_prefactor(0, r);
                prefactors.push_back(p0.get_link_acceptance_prefactor());
                r_values.push_back(r);
                bc_0 = bound_counts, uc_0 = unbound_counts;
                std::cout << "r = " << r << std::endl;
                std::cout << p0.get_link_acceptance_prefactor() << std::endl;
            }

            i++;
        }
        //double_vector_to_txt(path + std::to_string(l) + "dangle_prefactors.txt", prefactors);
        //double_vector_to_txt(path + std::to_string(l) + "dangle_r_values.txt", r_values);

        //double_vector_to_txt(path + std::to_string(l) + "dangle_helix_weights.txt", helix_weights_l);
        //double_vector_to_txt(path + std::to_string(l) + "dangle_old_linker_weights.txt", old_linker_weights);
        //double_vector_to_txt(path + std::to_string(l) + "dangle_new_linker_weights.txt", new_linker_weights);
        double_vector_to_txt(path + std::to_string(l) + "dangle_old_weights.txt", old_weights);
        double_vector_to_txt(path + std::to_string(l) + "dangle_new_weights.txt", new_weights);

        //double_vector_to_txt(path + std::to_string(l) + "dangle_old_dangle_weights.txt", old_dangle_weights);
        //double_vector_to_txt(path + std::to_string(l) + "dangle_new_dangle_weights.txt", new_dangle_weights);

        avg_helix_weight[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        avg_prefactors[l] = sum_of_elements(prefactors) / static_cast<double>(prefactors.size());


    }
    double_vector_to_txt(path + "dangle_avg_prefactors.txt", avg_prefactors);
    double_vector_to_txt(path + "dangle_avg_helix_weights.txt", avg_helix_weight);




}

void long_helix_simulation(int N_init, int NMC, int update_rate, std::string path)
{
    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };
    //int update_rate{ 500 };
    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };
    int max_helix_bps{ 10 };

    double r{ 0 };

    // convergence data
    std::vector<double> prefactors, r_values;
    prefactors.push_back(init_prefactor);

    // other data
    std::vector<double> old_weights, new_weights, helix_weights_l, old_linker_weights, new_linker_weights;
    // averages data
    std::vector<double> avg_helix_weight(max_helix_bps), avg_prefactors(max_helix_bps);

    int N0{ N_init }, N{ N0 };
    std::vector<int> helix(4);
    int cc{ 0 };
    bool update_r{ true };

    for (int l{ 0 }; l < max_helix_bps; l++) {

        // the helix is going to get longer
        N = N0 + 2 * l;
        polymer p0(N, true);
        helix = { 0 , 2 + l, N - 3 - l, N - 1};
        p0.set_search_results(helix);
        p0.set_link_acceptance_prefactor(init_prefactor);


        // resetting data 
        old_weights.clear();
        new_weights.clear();
        old_linker_weights.clear(), new_linker_weights.clear();
        helix_weights_l.clear();
        prefactors.clear();
        r_values.clear();
        cc = 0;
        i = 0;
        update_r = true;
        bound_counts = 0, unbound_counts = 0;
        bc_0 = 0, uc_0 = 0;
        r = 0;


        while (update_r == true) {
            // convergence condition.
            if (i % update_rate == 0 && i != 0 && update_r == true) {
                if (std::abs(std::log(r)) < epsilon) {
                    cc++;
                    if (cc >= 3) {
                        std::cout << "it converged" << std::endl;
                        update_r = false;
                    }
                }
                else if (i != 0) {
                    cc = 0;
                }

            }

            std::cout << "iteration " << i << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                //std::cout << "link branch" << std::endl;
                bool m;
                p0.link_move(m);
            }
            else {
                //std::cout << "unlink branch" << std::endl;

                bool m;
                p0.unlink_move(m);
            }

            // get data
            double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

            if (old_w != 0 && new_w != 0) {

                old_weights.push_back(p0.get_old_config_weight());
                new_weights.push_back(p0.get_new_config_weight());

                helix_weights_l.push_back(p0.get_helix_weight());

                old_linker_weights.push_back(p0.get_hairpin_weight(false));
                new_linker_weights.push_back(p0.get_hairpin_weight(true));


            }
            // histogram of bound and unbound counts
            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            // update the prefactor in the acceptance rule
            if ((i + 1) % update_rate == 0 && i != 0 && update_r == true) {

                r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
                r = 1 / r;
                p0.set_link_acceptance_prefactor(0, r);
                prefactors.push_back(p0.get_link_acceptance_prefactor());
                r_values.push_back(r);
                bc_0 = bound_counts, uc_0 = unbound_counts;
                std::cout << "r = " << r << std::endl;
                std::cout << p0.get_link_acceptance_prefactor() << std::endl;
            }

            i++;
        }
        //double_vector_to_txt(path + std::to_string(l) + "long_helix_prefactors.txt", prefactors);
        //double_vector_to_txt(path + std::to_string(l) + "long_helix_r_values.txt", r_values);

        //double_vector_to_txt(path + std::to_string(l) + "long_helix_helix_weights.txt", helix_weights_l);
        //double_vector_to_txt(path + std::to_string(l) + "long_helix_old_linker_weights.txt", old_linker_weights);
        //double_vector_to_txt(path + std::to_string(l) + "long_helix_new_linker_weights.txt", new_linker_weights);
        //double_vector_to_txt(path + std::to_string(l) + "long_helix_old_weights.txt", old_weights);
        //double_vector_to_txt(path + std::to_string(l) + "long_helix_new_weights.txt", new_weights);


        avg_helix_weight[l] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        avg_prefactors[l] = sum_of_elements(prefactors) / static_cast<double>(prefactors.size());


    }
    double_vector_to_txt(path + "long_helix_avg_prefactors.txt", avg_prefactors);
    double_vector_to_txt(path + "long_helix_avg_helix_weights.txt", avg_helix_weight);



}

void long_linker_simulation(int N_init, int NMC, int update_rate, std::string path)
{


    int i{ 0 };
    double u_branch;
    int bound_counts{ 0 }, unbound_counts{ 0 };
    int bc_0{ 0 }, uc_0{ 0 };
    //int update_rate{ 500 };
    double epsilon{ 2 };// for convergence
    double init_prefactor{ 1 };
    int max_linker_length{ 10 };

    double r{ 0 };

    // convergence data
    std::vector<double> prefactors, r_values;
    prefactors.push_back(init_prefactor);

    // other data
    std::vector<double> old_weights, new_weights, helix_weights_l, old_linker_weights, new_linker_weights;
    // averages data
    std::vector<double> avg_helix_weight(max_linker_length), avg_prefactors(max_linker_length);

    int N0{ N_init }, N{ N0 };
    std::vector<int> helix(4);
    int cc{ 0 };
    bool update_r{ true };
    int h{ 0 };
    for (int l{ 0 }; l < max_linker_length; l++) {
        h = 13 - l;
        // the helix is going to get longer
        N = N0 + l;
        polymer p0(N, true);
        helix = { 0 , 2, N - 3, N - 1 };
        p0.set_search_results(helix);
        p0.set_link_acceptance_prefactor(init_prefactor);


        // resetting data 
        old_weights.clear();
        new_weights.clear();
        old_linker_weights.clear(), new_linker_weights.clear();
        helix_weights_l.clear();
        prefactors.clear();
        r_values.clear();
        cc = 0;
        i = 0;
        update_r = true;
        bound_counts = 0, unbound_counts = 0;
        bc_0 = 0, uc_0 = 0;
        r = 0;


        while (update_r == true) {
            // convergence condition.
            if (i % update_rate == 0 && i != 0) {
                if (std::abs(std::log(r)) < epsilon) {
                    cc++;
                    if (cc >= 3) {
                        std::cout << "it converged" << std::endl;
                        update_r = false;
                    }
                }
                else if (i != 0) {
                    cc = 0;
                }

            }

            std::cout << "iteration " << i << std::endl;
            u_branch = rand2(0, 1);
            if (u_branch < 0.5) {// link_branch
                //std::cout << "link branch" << std::endl;
                bool m;
                p0.link_move(m);
            }
            else {
                //std::cout << "unlink branch" << std::endl;

                bool m;
                p0.unlink_move(m);
            }

            // get data
            double old_w{ p0.get_old_config_weight() }, new_w{ p0.get_new_config_weight() };

            if (old_w != 0 && new_w != 0) {

                old_weights.push_back(p0.get_old_config_weight());
                new_weights.push_back(p0.get_new_config_weight());

                helix_weights_l.push_back(p0.get_helix_weight());

                old_linker_weights.push_back(p0.get_hairpin_weight(false));
                new_linker_weights.push_back(p0.get_hairpin_weight(true));


            }
            // histogram of bound and unbound counts
            if (p0.linked() == true) {
                bound_counts++;
            }
            else {
                unbound_counts++;
            }
            // update the prefactor in the acceptance rule
            if ((i + 1) % update_rate == 0 && i != 0 && update_r == true) {

                r = (static_cast<double>(bound_counts - bc_0) + 1) / (static_cast<double>(unbound_counts - uc_0) + 1);
                r = 1 / r;
                p0.set_link_acceptance_prefactor(0, r);
                prefactors.push_back(p0.get_link_acceptance_prefactor());
                r_values.push_back(r);
                bc_0 = bound_counts, uc_0 = unbound_counts;
                std::cout << "r = " << r << std::endl;
                std::cout << p0.get_link_acceptance_prefactor() << std::endl;
            }

            i++;
        }
        //double_vector_to_txt(path + std::to_string(l) + "long_linker_prefactors.txt", prefactors);
        //double_vector_to_txt(path + std::to_string(l) + "long_linker_r_values.txt", r_values);

        //double_vector_to_txt(path + std::to_string(l) + "long_linker_helix_weights.txt", helix_weights_l);
        //double_vector_to_txt(path + std::to_string(l) + "long_linker_old_linker_weights.txt", old_linker_weights);
        //double_vector_to_txt(path + std::to_string(l) + "long_linker_new_linker_weights.txt", new_linker_weights);
        //double_vector_to_txt(path + std::to_string(l) + "long_linker_old_weights.txt", old_weights);
        //double_vector_to_txt(path + std::to_string(l) + "long_linker_new_weights.txt", new_weights);


        //avg_helix_weight[h] = sum_of_elements(helix_weights_l) / static_cast<double>(helix_weights_l.size());
        avg_prefactors[h] = sum_of_elements(prefactors) / static_cast<double>(prefactors.size());


    }
    double_vector_to_txt(path + "long_linker_avg_prefactors.txt", avg_prefactors);
    double_vector_to_txt(path + "long_linker_avg_helix_weights.txt", avg_helix_weight);


}

void equilibration_sim(int N_monomers, int NMC,  int init_update_rate, std::string path)
{
    
    int N_values = 10;
    std::vector<double> initial_prefactors(N_values+1);
    double pm{ static_cast<double>(N_values) / 2.0 };

    initial_prefactors[0] = 1;
    for (int i = 0; i < N_values/2; i++)
    {
        initial_prefactors[1 + 2*i] = std::pow(10.0, i + 1);
        initial_prefactors[2 + 2*i] = std::pow(10.0, -(i + 1));

    }
    std::vector<int> helix{ 0,2,N_monomers - 3,N_monomers - 1 };
    double first_run_prefactor;
    for (int p = 0; p < initial_prefactors.size(); p++)
    {
        first_run_prefactor = get_prefactor(initial_prefactors[p], helix, N_monomers, NMC , init_update_rate, path + std::to_string(initial_prefactors[p]));
        //correlation_length(first_run_prefactor, N_monomers, helix, NMC, path + std::to_string(initial_prefactors[p]));

    }
}
