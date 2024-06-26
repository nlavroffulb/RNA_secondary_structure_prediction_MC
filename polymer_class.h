#pragma once
#pragma once

#include "monomer_class.h"
#include "helix_struct.h"
#include "rosenbluth_growth_class.h"
#include "polymer_generation.h"
#include "helix_generation.h"
#include "math_functions.h"




class polymer {
protected:
    monomers chain;
    helices helix_list;
    std::vector<std::vector<int>> search_results;
    std::vector<int> extendable_structures;

    std::vector<int> zipped_structures;

    std::vector<std::vector<double>> excluded_volume;


    // data for acceptance probability calculation.
    std::vector<std::vector<int>> new_growth_lims, old_growth_lims;
    std::vector<std::vector < double >> new_config_positions;

    rosenbluth_growth *old_weights;
    rosenbluth_growth* new_weights;

    // important for link_unlink
    double helix_R_factor;
    helix_struct* proposed_link_helix;

    // important for swivel
    double old_helix_weight, new_helix_weight;
    std::vector<double> new_running_centre;
    std::vector<int> zip_unzip_structure;

    double helix_interaction_weight{1.0};
    

    bool rosenbluth_switch{true};// 1 equals on, 0 equals off.

    // data for analysis
    std::vector<double> linker_weights;
    std::vector<double> helix_weights;
    std::vector<double> dangling_weights;
    double old_config_weight, new_config_weight;

    double link_acceptance_prefactor;

private:
    int N{ 0 };

public:

    //constructors
    polymer() = default;
    polymer(int num);
    polymer(int num, bool rosenbluth, int struct_size=3);
    polymer(int num, bool rosenbluth, std::vector<std::string>& sequence, int struct_size = 3);


    //destructor
    ~polymer();

    //useful functions if we want to print or visualized things.
    std::vector<std::vector<double>> generate_vector_of_monomer_positions();
    void print_monomer_positions();
    void print_bases();
    void print_bases1();
    void output_for_ovito4(std::string filename="");
    void ovito_bonds(std::string filename = "");
    std::vector<std::string> get_sequence();
    std::vector<std::vector<double>> get_subpositions(std::vector<int> limit);
    //operator overload function. allows us to write chain[i] as a reference to a monomer object.
    monomer* operator[](int i);
    int get_chain_length();
    ////////////////////////////////////////////////////////////////////////////////////////
    // at different moments in the simulation we have to update different member variables or other quantities
    //usually based on whether a move is accepted or rejected. we also need to search for a compatible regions.
    std::vector<std::vector<int>> structure_search(int n=3);
    void update_large_struct_list();
    void update_extensible_structures();
    void update_positions();
    void reset_positions();
    void neighbouring_linkers(std::vector<int> linked_monomers);
    void neighbouring_linkers();
    void update_excluded_volume(std::vector<std::vector<int>>& growth_limits, std::vector<int> helix = {});
    void get_linked_monomers(std::vector<int>& links);
    bool overlapping_monomers(std::vector<int>& reaction_region);
    void structure_extension(std::vector<int>& s, int& index);// returns (by reference) the possible extension of a double helix
    bool zip_structure_overlap(std::vector<int> s_z, int side);// we shouldn't zip into another double helix structure.
    ////////////////////////////////////////////////////////////////////////////////////////
    // our double helix structures are represented by [a b c d] where a b c and d are not consecutive to each other usually, they
    //define limits. sometimes its useful to "unpack" that representation: [a,a+1,...,b, c,c+1,...,d] and vice versa.
    void unpack_region(std::vector<int>& reaction_region);
    void pack_region(std::vector<int>& reaction_region);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // link specific functions
    bool linked();
    void sample_link_region(std::vector<int>& link, int& alpha, int& beta, int& ss_index);
    bool reject_link(std::vector<int>& link_region, int alpha, int beta);
    void link(std::vector<int>& link_region, int alpha, int beta, int ss_index);
    void link_update(int ss_index);
    void link_growth_limits(std::vector<int> helix, int alpha, int beta, std::vector<std::vector<int>>& growth_limits);
    void link_growth_limits1(std::vector<int> helix, int alpha, int beta, std::vector<std::vector<int>>& growth_limits, bool forward_move = true);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // unlink specific functions
    void sample_unlink_region(int& helix_index);
    void unlink(int helix_index);
    void unlink_update(int s_index); // if move is accepted
    void unlink_growth_limits(std::vector<int> helix, int alpha, int beta, std::vector<std::vector<int>>& growth_limits);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // zip related functions
    void zip_growth_limits(std::vector<int>& s, int side, std::vector<std::vector<int>>& growth_limits);
    bool reject_zip(std::vector<int>& s, std::vector < std::vector<double>> new_positions, int side);
    void zip(bool &success, int &sigma, int &s_index);
    void sample_zip_move(bool& success, int& s_index);
    void sample_zip_region(std::vector<int>& region_to_zip, std::vector<int>& extension, int side);
    void zip_update(int s_index, int side);
    std::vector<std::vector<double>> add_bp_to_helix(helix_struct* double_helix, int side);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // unzip related functions
    void sample_unzip_region(std::vector<int>& region_to_unzip, int side);
    void unzip(bool& success, int& side, int& s_index);
    void unzip_growth_linits(std::vector<int>& s, int side, std::vector<std::vector<int>>& growth_limits);
    void unzip_update(int s_index, int side);

    ////////// acceptance rule calculation functions ////////////////////////
    double link_acceptance(bool link_or_unlink);
    double zip_acceptance(bool zip_or_unzip, int side);
    double swivel_acceptance();
    double weight_gen_config(bool forward_move);
    double p_gen_ratio(std::vector<int> limit);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    void helix_excluded_volume_interaction(std::vector<std::vector<double>>& helix_positions);
    void grow_limits(std::vector<std::vector<int>>& limits, int alpha, bool forward_move);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // swivel specific functions
    bool reject_spin(helix_struct* double_helix, double angle, std::vector<double> rot_axis);
    bool reject_corkscrew(helix_struct* double_helix, double angle);
    bool reject_translate(helix_struct* double_helix, std::vector<double> &translation);
    void swivel_growth_limits(helix_struct* dh, std::vector<std::vector<int>>& limits);
    std::vector<double> centre_of_mass(std::vector<int> structure);
    void corkscrew_update(helix_struct* dh);
    void spin_update(helix_struct* dh, std::vector<double>& rot_axis, double angle);
    void translate_update(helix_struct* dh, std::vector<double>&translation);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<std::vector<double>> get_helix_monomer_positions(std::vector<int>& momonomers);

    void translate_move(helix_struct* dh);
    void corkscrew_move(helix_struct* dh);
    void spin_move(helix_struct* dh);
    void swivel_move();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    double configuration_energy();
    void add_linker_weight(bool link_move);
    double get_hairpin_weight();

    double get_subsection_weight(std::vector<int> limit, bool new_config);
    double get_hairpin_weight(bool new_weight);
    double get_helix_weight();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //void link_move(int multiples_ds=0);

    void always_accept_link_move(bool& overstretch_reject);
    void KH_unlink_move(bool& overstretch_reject, std::vector<std::vector<int>> unlink_prohibited);

    void link_move(bool& overstretch_reject, int multiple_ds=0);

    //void unlink_move();
    void unlink_move(bool& overstretch_reject);
    void force_unbound_state();
    void set_link_acceptance_prefactor(double prefactor, double multiplier=-1);
    double get_link_acceptance_prefactor();
    
    void set_search_results(std::vector<int> helix);
    void set_search_results(std::vector<std::vector<int>>& helices);

    int get_num_monomers();
    int get_num_helices();
    double get_old_config_weight();
    double get_new_config_weight();
};

