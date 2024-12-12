import anvil.server
import anvil.media
import random
import numpy as np
import re
import math
from sys import exit
from collections import Counter
import matplotlib.pyplot as plt
import anvil.mpl_util
import pandas as pd

# This is a server module. It runs on the Anvil server,
# rather than in the user's browser.
#
# To allow anvil.server.call() to call functions here, we mark
# them with @anvil.server.callable.

#### GENERAL FUNCTIONS ####

def chitosan_generator(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks):
    """Generates one chitosan molecule of either...
    A) a specific DP and overall FA. If the given FA_pattern equals the given
    overall FA, the chitosan molecule features a random PA. Otherwise, every
    nth unit (n = pattern) features either an overrepresentation of A-units
    (FA_pattern > FA) or of B-units (FA_pattern < FA).
    B) a specific DP and exactly defined blocks of A- and B-units, e.g. always
    3 As followed by always 2 Bs. It is chosen at random whether the molecule
    starts with an A- or B-block.

    returns chitosan molecule in format: ["B", "A", "B", "B", "A", ...]"""

    # OPTION A: specific FA and pattern
    if A_blocks == 0 and B_blocks == 0:
        # calculate FA_non_pattern of all residues that are not part of pattern
        pattern_spaces = int(np.ceil(DP / pattern))
        non_pattern_spaces = DP - pattern_spaces
        proportion_pattern = pattern_spaces / DP
        proportion_non_pattern = non_pattern_spaces / DP
        FA_non_pattern = ((overall_FA - proportion_pattern * FA_pattern) /
                          proportion_non_pattern)
        if FA_non_pattern < 0:
            print("ERROR: combination of overall_FA, FA_pattern and pattern not \
    sensible!")
            exit()
        # generate A/B-sequence lists for residues that are (not) part of pattern
        pattern_list = ["A" if random.choices(["A", "B"],
                                              [FA_pattern,
                                              1 - FA_pattern])
                        == ["A"]
                        else "B" for i in range(int(pattern_spaces))]
        non_pattern_list = ["A" if random.choices(["A", "B"],
                                                  [FA_non_pattern,
                                                  1 - FA_non_pattern])
                            == ["A"]
                            else "B" for i in range(int(non_pattern_spaces))]
        # insert pattern sequence into non-pattern sequence at every nth product
        nth_index = list(range(0, DP, pattern))
        position_pattern_list = 0
        for position in nth_index:
            non_pattern_list.insert(position, pattern_list[position_pattern_list])
            position_pattern_list += 1
        final_molecule = non_pattern_list

    # OPTION B: specific A- and B-blocks
    else:
        block_polymer = []
        A_seq = ["A"] * A_blocks
        B_seq = ["B"] * B_blocks
        A_or_B = random.choices(["A", "B"], [0.5, 0.5])
        if A_or_B == ["A"]:
            block_polymer.extend(A_seq)
        elif A_or_B == ["B"]:
            block_polymer.extend(B_seq)
        while len(block_polymer) < DP:
            if block_polymer[-1] == "A":
                block_polymer.extend(B_seq)
            elif block_polymer[-1] == "B":
                block_polymer.extend(A_seq)
        final_molecule = block_polymer[0:DP]
    return final_molecule


def chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks,
    number_molecules):
    """Generates library containing a given number of chitosan molecules
    created via chitosan_generator.

    returns chitosan library: ["molecule1", "molecule2", ...] with molecule1
    e.g. = "BABBBAAB" """

    library = []
    for i in range(0, number_molecules):
        molecule = "".join(chitosan_generator(DP, overall_FA, pattern,
                                              FA_pattern, A_blocks, B_blocks))
        library.append(molecule)
    return library


def enzyme(substrate, minus_specificity, plus_specificity, efficiency):
    """Simulates cleavage of an individual chitosan molecule (given as string)
    by an enzyme of given specificity. The plus and minus specificities always
    consist of three characters, e.g. AAA, ABB, AXX, AB., ...
    with A = unit A, B = unit B, X = A or B but subsite needs to be
    occupied, . = A or B or unoccupied.

    returns a list of products: [product1, product2, ...]"""

    minus = [*minus_specificity]
    minus = ["[^-]" if pref == "X" else pref for pref in minus]
    minus = ["." if pref == "_" else pref for pref in minus]
    plus = [*plus_specificity]
    plus = ["[^-]" if pref == "X" else pref for pref in plus]
    plus = ["." if pref == "_" else pref for pref in plus]
    clvg_yn = (True, False)
    clvg_prob = [efficiency, 1-efficiency]
    product_list = []
    dig_substrate = [*substrate]
    clvg_idx_list = []
    # identify potential cleavage sites and their indices
    for idx, monomer in enumerate(dig_substrate[:-5]):
        if re.match(minus[0], dig_substrate[idx]) and \
           re.match(minus[1], dig_substrate[idx + 1]) and \
           re.match(minus[2], dig_substrate[idx + 2]) and \
           re.match(plus[0], dig_substrate[idx + 3]) and \
           re.match(plus[1], dig_substrate[idx + 4]) and \
           re.match(plus[2], dig_substrate[idx + 5]):
            clvg_idx_list.append(idx)
    # shuffle to simulate random endo-cleavage of enzyme
    random.shuffle(clvg_idx_list)
    # cleave at identified sites if former cut did not destroy cleavage site
    for clvg_no in range(len(clvg_idx_list)):
        if re.match(minus[0], dig_substrate[clvg_idx_list[clvg_no]]) and \
           re.match(minus[1], dig_substrate[clvg_idx_list[clvg_no] + 1]) and \
           re.match(minus[2], dig_substrate[clvg_idx_list[clvg_no] + 2]) and \
           re.match(plus[0], dig_substrate[clvg_idx_list[clvg_no] + 3]) and \
           re.match(plus[1], dig_substrate[clvg_idx_list[clvg_no] + 4]) and \
           re.match(plus[2], dig_substrate[clvg_idx_list[clvg_no] + 5]):
            dig_substrate.insert(clvg_idx_list[clvg_no] + 3, "-")
            clvg_idx_list = [clvg_idx + 1 if (
                clvg_idx > clvg_idx_list[clvg_no]) else
                clvg_idx for clvg_idx in clvg_idx_list]
    dig_substrate_str = "".join(dig_substrate)
    dig_substrate_list = dig_substrate_str.split("-")
    # implement efficiency < 1 by putting products back together
    partial_dig_sub_list = [dig_substrate_list[0]]
    for product in dig_substrate_list[1:]:
        if random.choices(clvg_yn, weights=clvg_prob) == [True]:
            partial_dig_sub_list.append(product)
        else:
            partial_dig_sub_list[-1] = partial_dig_sub_list[-1]+product
    for product in range(1, len(partial_dig_sub_list)-1):
        product_list.append(partial_dig_sub_list[product])
    return product_list


def nano2(substrate, cuts):
    """Simulates cleavage of an individual chitosan molecule (given as string)
    using sodium nitrite (NaNO3) that cleaves after each B-unit

    returns a list of products: [product1, product2, ...]"""

    product_list = []
    dig_substrate = [*substrate]
    clvg_idx_list = []
    # identify potential cleavage sites and their indices
    for idx, monomer in enumerate(dig_substrate[:-1]):
        if re.match("B", dig_substrate[idx]) and \
           re.match(".", dig_substrate[idx + 1]):
            clvg_idx_list.append(idx)
    # shuffle to simulate random endo-cleavage
    random.shuffle(clvg_idx_list)
    # cleave at identified sites if former cut did not destroy cleavage site
    for clvg_no in range(len(clvg_idx_list)):
        if cuts > 0:
            if re.match("B", dig_substrate[clvg_idx_list[clvg_no]]) and \
               re.match(".", dig_substrate[clvg_idx_list[clvg_no] + 1]):
                dig_substrate.insert(clvg_idx_list[clvg_no] + 1, "-")
                clvg_idx_list = [clvg_idx + 1 if (
                    clvg_idx > clvg_idx_list[clvg_no]) else
                    clvg_idx for clvg_idx in clvg_idx_list]
                cuts -= 1
        else:
            break
    product_str = "".join(dig_substrate)
    product_list = product_str.split("-")
    return product_list


def FA_pattern_calc(overall_FA, pattern, strength):
    """Calculates the FA_pattern for a given overall FA and pattern that
    corresponds to a certain strength of this pattern (value between 0 and 1):
    0: the overall FA and the FA of the pattern are the same, resulting
    in a random PA molecule.
    1: the strength of the pattern is maximal, meaning the overrepresentation
    of A-units at every nth position is as strong as possible for the given
    overall FA.
    ONLY WORKS FOR A-PATTERNS AKA OVERREPRESENTATION OF A-UNITS!"""

    if strength == 0:
        FA_pattern = overall_FA
    else:
        strength_factor = 1 + strength * (pattern - 1)
        FA_pattern_raw = 0.98 * strength_factor * overall_FA
        if FA_pattern_raw > 1:
            FA_pattern = 1
        else:
            FA_pattern = FA_pattern_raw
    return(FA_pattern)


def calc_DP(oligomer):
    """Takes the string of an oligomer as input (e.g. A1B12) and returns the
    corresponding DP (e.g. 13)"""
    temp = "0"
    sum = 0
    for ch in oligomer:
        if ch.isdigit():
            temp += ch
        else:
            sum += int(temp)
            temp = "0"
    return sum + int(temp)


def diads_PA(library):
    """Determines diad frequencies (AA, AB, BA and BB) and calculates the PA
    value for a given library.

    returns respective values in a list: [AA, AB, BA, BB, PA]"""

    # count diads
    AA = 0
    AB = 0
    BA = 0
    BB = 0
    for molecule in library:
        for idx, monomer in enumerate(molecule[:-1]):
            if monomer == "A":
                if molecule[idx + 1] == "A":
                    AA += 1
                elif molecule[idx + 1] == "B":
                    AB += 1
            elif monomer == "B":
                if molecule[idx + 1] == "A":
                    BA += 1
                elif molecule[idx + 1] == "B":
                    BB += 1
    sum_all_diads = AA + AB + BA + BB
    # calculate diad frequencies
    F_AA = AA / sum_all_diads
    F_BB = BB / sum_all_diads
    F_BA = BA / sum_all_diads
    F_AB = AB / sum_all_diads
    F_AB_BA = (AB + BA) / sum_all_diads
    # calculate PA value based on Kumirska et al. (2009)
    PA = F_AB_BA / (2 * F_AA + F_AB_BA) + F_AB_BA / (2 * F_BB + F_AB_BA)
    output = [F_AA, F_AB, F_BA, F_BB, PA]
    return output


def triads(library):
    """Determines triad frequencies (AAA, AAB, BAA, ABA, ABB, BBA, BAB, BBB)
    and calculates their fractions for a given chitosan library.

    returns values in list: [AAA, AAB, BAA, ABA, ABB, BBA, BAB, BBB]"""

    # count triads
    AAA = 0
    AAB = 0
    BAA = 0
    ABA = 0
    ABB = 0
    BBA = 0
    BAB = 0
    BBB = 0
    for molecule in library:
        for idx, monomer in enumerate(molecule[:-2]):
            if monomer == "A":
                if molecule[idx + 1] == "A":
                    if molecule[idx + 2] == "A":
                        AAA += 1
                    elif molecule[idx + 2] == "B":
                        AAB += 1
                elif molecule[idx + 1] == "B":
                    if molecule[idx + 2] == "A":
                        ABA += 1
                    elif molecule[idx + 2] == "B":
                        ABB += 1
            elif monomer == "B":
                if molecule[idx + 1] == "A":
                    if molecule[idx + 2] == "A":
                        BAA += 1
                    elif molecule[idx + 2] == "B":
                        BAB += 1
                elif molecule[idx + 1] == "B":
                    if molecule[idx + 2] == "A":
                        BBA += 1
                    elif molecule[idx + 2] == "B":
                        BBB += 1
    sum_all_triads = AAA + AAB + BAA + ABA + ABB + BBA + BAB + BBB
    # calculate triad frequencies
    F_AAA = AAA / sum_all_triads
    F_AAB = AAB / sum_all_triads
    F_BAA = BAA / sum_all_triads
    F_ABA = ABA / sum_all_triads
    F_ABB = ABB / sum_all_triads
    F_BBA = BBA / sum_all_triads
    F_BAB = BAB / sum_all_triads
    F_BBB = BBB / sum_all_triads
    output = [F_AAA, F_AAB, F_BAA, F_ABA, F_ABB, F_BBA, F_BAB, F_BBB]
    return output


def blocks(product_list, A_cutoff, DP_cutoff, ion_eff):
    """Takes list of chitinosanase products (["BA", "BBBA", "BA"]) and
    calculates the block sizes An, Bn, Aw and Bw.
    It is possible to:
    -set an A_cutoff: products with more A-units than this cutoff are excluded
    due to their potential insolublity
    -set a DP_cutoff: products of higher DP are excluded because they are
    hard to detect via MS and not included when analyzing experimental data
    -consider or not consider the different ionization_efficiencies of the
    different products

    returns block_sizes: (block_An, block_Bn, block_Aw, block_Bw)"""

    # remove products bigger than DP_cutoff and with more As than A_cutoff
    a_b_tuple_list = []
    if A_cutoff == False and DP_cutoff == False:
        for oligomer in product_list:
            a_b_tuple_list.append((oligomer.count("A"), oligomer.count("B")))
    elif A_cutoff != False and DP_cutoff == False:
        for oligomer in product_list:
            if oligomer.count("A") <= A_cutoff:
                a_b_tuple_list.append((oligomer.count("A"), oligomer.count("B")))
    elif A_cutoff == False and DP_cutoff != False:
        for oligomer in product_list:
            if len(oligomer) <= DP_cutoff:
                a_b_tuple_list.append((oligomer.count("A"), oligomer.count("B")))
    elif A_cutoff != False and DP_cutoff != False:
        for oligomer in product_list:
            if len(oligomer) <= DP_cutoff and oligomer.count("A") <= A_cutoff:
                a_b_tuple_list.append((oligomer.count("A"), oligomer.count("B")))
    # create dictionary with numbers of each chitinosanase product
    product_dict = dict(Counter(a_b_tuple_list))
    # if applicable: consider different ionization efficiencies per DP
    if ion_eff:
        for product in product_dict:
            if sum(product) == 2:
                product_dict[product] *= 0.3
            elif sum(product) == 3:
                product_dict[product] *= 1.0
            elif sum(product) == 4:
                product_dict[product] *= 0.65
            elif sum(product) > 4 and sum(product) <= 8:
                product_dict[product] *= 0.35
            elif sum(product) > 8:
                product_dict[product] *= 0.2
    # calculate block-sizes, product[0] = number As, product[1] = number Bs
    block_An = 0
    block_Bn = 0
    Aw_sum = 0
    Bw_sum = 0
    DP_sum = 0
    sum_all_products = sum(product_dict.values())
    for product in product_dict:
        rel_abundance_product = (product_dict[product] /
                                 sum_all_products)
        product_DP = sum(product)
        block_An += product[0] * rel_abundance_product
        block_Bn += product[1] * rel_abundance_product
        Aw_sum += product[0] * rel_abundance_product * product_DP
        Bw_sum += product[1] * rel_abundance_product * product_DP
        DP_sum += rel_abundance_product * product_DP
    block_Aw = Aw_sum/DP_sum
    block_Bw = Bw_sum/DP_sum
    block_sizes = [block_An, block_Bn, block_Aw, block_Bw]
    return block_sizes


#### DP HISTOGRAM ####

@anvil.server.callable
def histogram_DP_FAp_molar_enzyme(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules, minus_specificity, plus_specificity, efficiency):
    # generate library
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and count products per DP and track FA
    product_dict = {}
    all_products = 0
    for substrate in library:
        product_list = enzyme(substrate, minus_specificity, plus_specificity, efficiency)
        for product in product_list:
            all_products += 1
            if len(product) in product_dict:
                product_dict[len(product)][0] += 1
                product_dict[len(product)][1] += product.count("A")
            else:
                product_dict[len(product)] = [1, product.count("A")]
    # calculate average FA and prepare output list
    DP_list = []
    average_FA_list = []
    mole_fraction_list = []
    for DP in range(min(product_dict), max(product_dict) + 1):
        if DP in product_dict:
            DP_list.append(DP)
            mole_fraction = product_dict[DP][0] / all_products
            mole_fraction_list.append(mole_fraction)
            average_FA = product_dict[DP][1] / (DP * product_dict[DP][0])
            average_FA_list.append(average_FA)
    output_list = [DP_list, mole_fraction_list, average_FA_list]
    return(output_list)

    
@anvil.server.callable
def histogram_DP_FAp_molar_nano2(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules, cuts):
    # generate library
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and count products per DP and track FA
    product_dict = {}
    all_products = 0
    for substrate in library:
        product_list = nano2(substrate, cuts)
        for product in product_list:
            all_products += 1
            if len(product) in product_dict:
                product_dict[len(product)][0] += 1
                product_dict[len(product)][1] += product.count("A")
            else:
                product_dict[len(product)] = [1, product.count("A")]
    # calculate average FA and prepare output list
    DP_list = []
    average_FA_list = []
    mole_fraction_list = []
    for DP in range(min(product_dict), max(product_dict) + 1):
        if DP in product_dict:
            DP_list.append(DP)
            mole_fraction = product_dict[DP][0] / all_products
            mole_fraction_list.append(mole_fraction)
            average_FA = product_dict[DP][1] / (DP * product_dict[DP][0])
            average_FA_list.append(average_FA)
    output_list = [DP_list, mole_fraction_list, average_FA_list]
    return(output_list)

    
@anvil.server.callable
def histogram_DP_FAp_weight_enzyme(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules, minus_specificity, plus_specificity, efficiency, mass_A, mass_B):
    # generate library
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and count products per DP and track FA
    product_dict = {}
    all_products = 0
    for substrate in library:
        product_list = enzyme(substrate, minus_specificity, plus_specificity, efficiency)
        for product in product_list:
            all_products += 1
            if len(product) in product_dict:
                product_dict[len(product)][0] += 1
                product_dict[len(product)][1] += product.count("A")
            else:
                product_dict[len(product)] = [1, product.count("A")]
    # calculate average FA and prepare output list
    DP_list = []
    average_FA_list = []
    mass_mole_fraction_list = []
    mass_mole_fraction_sum = 0
    weight_fraction_list = []
    for DP in range(min(product_dict), max(product_dict) + 1):
        if DP in product_dict:
            DP_list.append(DP)
            mole_fraction = product_dict[DP][0] / all_products
            average_FA = product_dict[DP][1] / (DP * product_dict[DP][0])
            average_FA_list.append(average_FA)
            average_mass = DP * (average_FA * mass_A + (1-average_FA) * mass_B) - (DP-1)*18
            mass_mole_fraction = mole_fraction * average_mass
            mass_mole_fraction_list.append(mass_mole_fraction)
            mass_mole_fraction_sum += mass_mole_fraction
    for mass_mole_fraction in mass_mole_fraction_list:
        weight_fraction = mass_mole_fraction / mass_mole_fraction_sum
        weight_fraction_list.append(weight_fraction)
    output_list = [DP_list, weight_fraction_list, average_FA_list]
    return(output_list)

    
@anvil.server.callable
def histogram_DP_FAp_weight_nano2(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules, cuts, mass_A, mass_B):
    # generate library
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and count products per DP and track FA
    product_dict = {}
    all_products = 0
    for substrate in library:
        product_list = nano2(substrate, cuts)
        for product in product_list:
            all_products += 1
            if len(product) in product_dict:
                product_dict[len(product)][0] += 1
                product_dict[len(product)][1] += product.count("A")
            else:
                product_dict[len(product)] = [1, product.count("A")]
    # calculate average FA and prepare output list
    DP_list = []
    average_FA_list = []
    mass_mole_fraction_list = []
    mass_mole_fraction_sum = 0
    weight_fraction_list = []
    for DP in range(min(product_dict), max(product_dict) + 1):
        if DP in product_dict:
            DP_list.append(DP)
            mole_fraction = product_dict[DP][0] / all_products
            average_FA = product_dict[DP][1] / (DP * product_dict[DP][0])
            average_FA_list.append(average_FA)
            average_mass = DP * (average_FA * mass_A + (1-average_FA) * mass_B) - (DP-1)*18 - 17 # 17 = mass difference between GlcN and anhydromannose
            mass_mole_fraction = mole_fraction * average_mass
            mass_mole_fraction_list.append(mass_mole_fraction)
            mass_mole_fraction_sum += mass_mole_fraction
    for mass_mole_fraction in mass_mole_fraction_list:
        weight_fraction = mass_mole_fraction / mass_mole_fraction_sum
        weight_fraction_list.append(weight_fraction)
    output_list = [DP_list, weight_fraction_list, average_FA_list]
    return(output_list)

    
@anvil.server.callable
def histogram_DP_strength_molar_enzyme(DP, overall_FA, pattern, strength, A_blocks, B_blocks, molecules, minus_specificity, plus_specificity, efficiency):
    # generate library
    FA_pattern = FA_pattern_calc(overall_FA, pattern, strength)
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and count products per DP and track FA
    product_dict = {}
    all_products = 0
    for substrate in library:
        product_list = enzyme(substrate, minus_specificity, plus_specificity, efficiency)
        for product in product_list:
            all_products += 1
            if len(product) in product_dict:
                product_dict[len(product)][0] += 1
                product_dict[len(product)][1] += product.count("A")
            else:
                product_dict[len(product)] = [1, product.count("A")]
    # calculate average FA and prepare output list
    DP_list = []
    average_FA_list = []
    mole_fraction_list = []
    for DP in range(min(product_dict), max(product_dict) + 1):
        if DP in product_dict:
            DP_list.append(DP)
            mole_fraction = product_dict[DP][0] / all_products
            mole_fraction_list.append(mole_fraction)
            average_FA = product_dict[DP][1] / (DP * product_dict[DP][0])
            average_FA_list.append(average_FA)
    output_list = [DP_list, mole_fraction_list, average_FA_list]
    return(output_list)

    
@anvil.server.callable
def histogram_DP_strength_molar_nano2(DP, overall_FA, pattern, strength, A_blocks, B_blocks, molecules, cuts):
    # generate library
    FA_pattern = FA_pattern_calc(overall_FA, pattern, strength)
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and count products per DP and track FA
    product_dict = {}
    all_products = 0
    for substrate in library:
        product_list = nano2(substrate, cuts)
        for product in product_list:
            all_products += 1
            if len(product) in product_dict:
                product_dict[len(product)][0] += 1
                product_dict[len(product)][1] += product.count("A")
            else:
                product_dict[len(product)] = [1, product.count("A")]
    # calculate average FA and prepare output list
    DP_list = []
    average_FA_list = []
    mole_fraction_list = []
    for DP in range(min(product_dict), max(product_dict) + 1):
        if DP in product_dict:
            DP_list.append(DP)
            mole_fraction = product_dict[DP][0] / all_products
            mole_fraction_list.append(mole_fraction)
            average_FA = product_dict[DP][1] / (DP * product_dict[DP][0])
            average_FA_list.append(average_FA)
    output_list = [DP_list, mole_fraction_list, average_FA_list]
    return(output_list)

    
@anvil.server.callable
def histogram_DP_strength_weight_enzyme(DP, overall_FA, pattern, strength, A_blocks, B_blocks, molecules, minus_specificity, plus_specificity, efficiency, mass_A, mass_B):
    # generate library
    FA_pattern = FA_pattern_calc(overall_FA, pattern, strength)
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and count products per DP and track FA
    product_dict = {}
    all_products = 0
    for substrate in library:
        product_list = enzyme(substrate, minus_specificity, plus_specificity, efficiency)
        for product in product_list:
            all_products += 1
            if len(product) in product_dict:
                product_dict[len(product)][0] += 1
                product_dict[len(product)][1] += product.count("A")
            else:
                product_dict[len(product)] = [1, product.count("A")]
    # calculate average FA and prepare output list
    DP_list = []
    average_FA_list = []
    mass_mole_fraction_list = []
    mass_mole_fraction_sum = 0
    weight_fraction_list = []
    for DP in range(min(product_dict), max(product_dict) + 1):
        if DP in product_dict:
            DP_list.append(DP)
            mole_fraction = product_dict[DP][0] / all_products
            average_FA = product_dict[DP][1] / (DP * product_dict[DP][0])
            average_FA_list.append(average_FA)
            average_mass = DP * (average_FA * mass_A + (1-average_FA) * mass_B) - (DP-1)*18
            mass_mole_fraction = mole_fraction * average_mass
            mass_mole_fraction_list.append(mass_mole_fraction)
            mass_mole_fraction_sum += mass_mole_fraction
    for mass_mole_fraction in mass_mole_fraction_list:
        weight_fraction = mass_mole_fraction / mass_mole_fraction_sum
        weight_fraction_list.append(weight_fraction)
    output_list = [DP_list, weight_fraction_list, average_FA_list]
    return(output_list)

    
@anvil.server.callable
def histogram_DP_strength_weight_nano2(DP, overall_FA, pattern, strength, A_blocks, B_blocks, molecules, cuts, mass_A, mass_B):
    # generate library
    FA_pattern = FA_pattern_calc(overall_FA, pattern, strength)
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and count products per DP and track FA
    product_dict = {}
    all_products = 0
    for substrate in library:
        product_list = nano2(substrate, cuts)
        for product in product_list:
            all_products += 1
            if len(product) in product_dict:
                product_dict[len(product)][0] += 1
                product_dict[len(product)][1] += product.count("A")
            else:
                product_dict[len(product)] = [1, product.count("A")]
    # calculate average FA and prepare output list
    DP_list = []
    average_FA_list = []
    mass_mole_fraction_list = []
    mass_mole_fraction_sum = 0
    weight_fraction_list = []
    for DP in range(min(product_dict), max(product_dict) + 1):
        if DP in product_dict:
            DP_list.append(DP)
            mole_fraction = product_dict[DP][0] / all_products
            average_FA = product_dict[DP][1] / (DP * product_dict[DP][0])
            average_FA_list.append(average_FA)
            average_mass = DP * (average_FA * mass_A + (1-average_FA) * mass_B) - (DP-1)*18 - 17 # 17 = mass difference between GlcN and anhydromannose
            mass_mole_fraction = mole_fraction * average_mass
            mass_mole_fraction_list.append(mass_mole_fraction)
            mass_mole_fraction_sum += mass_mole_fraction
    for mass_mole_fraction in mass_mole_fraction_list:
        weight_fraction = mass_mole_fraction / mass_mole_fraction_sum
        weight_fraction_list.append(weight_fraction)
    output_list = [DP_list, weight_fraction_list, average_FA_list]
    return(output_list)

    
@anvil.server.callable
def make_histogram(data, y_label, file_name):
    datax = data[0]
    min_x = min(datax)
    max_x = max(datax)
    datay = data[1]
    plt.figure(1, figsize=(20,7), dpi=400)    
    ax = plt.axes()
    plt.bar(datax, height=datay, color="#009688")
    plt.xlabel("DP", fontsize=28)
    plt.ylabel(y_label, fontsize=28)
    if max_x - min_x > 300:
        plt.xticks(np.arange(min_x, max_x+1, 50.0), fontsize=22)
    elif max_x - min_x > 200:
        plt.xticks(np.arange(min_x, max_x+1, 20.0), fontsize=22)
        plt.xticks(np.arange(min_x, max_x+1, 20.0), fontsize=22)
    elif max_x - min_x > 100:
        plt.xticks(np.arange(min_x, max_x+1, 10.0), fontsize=22)
    elif max_x - min_x > 80:
        plt.xticks(np.arange(min_x, max_x+1, 5.0), fontsize=22)
    elif max_x - min_x > 60:
        plt.xticks(np.arange(min_x, max_x+1, 1.0), fontsize=22)
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 4 != 0]
    elif max_x - min_x > 40:
        plt.xticks(np.arange(min_x, max_x+1, 1.0), fontsize=22)
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 3 != 0]
    elif max_x - min_x > 25:
        plt.xticks(np.arange(min_x, max_x+1, 1.0), fontsize=22)
        [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 2 != 0]
    elif max_x - min_x <=25:
        plt.xticks(np.arange(min_x, max_x+1, 1.0), fontsize=22)        
    plt.yticks(fontsize=22)
    # Return this plot as a PNG image in a Media object
    return anvil.mpl_util.plot_image(filename=("%s.png" % file_name))

    
@anvil.server.callable
def histogram_csv(data, ylabel, filename):
    df = pd.DataFrame(data).transpose()
    df.columns = ["DP of product", ylabel, "average fraction of unit A"]
    df.to_csv('/tmp/data.csv', index=False, encoding='utf-8-sig')
    csv_file = anvil.media.from_file('/tmp/data.csv', 'csv', ('%s.csv' % filename))
    return csv_file

#### PRODUCT PROFILE ####

@anvil.server.callable
def profile_FAp_molar_enzyme(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules, minus_specificity, plus_specificity, efficiency):
    # generate library
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and save product compositions
    profile_list = []
    oligomer_list = []
    for substrate in library:
        product_list = enzyme(substrate, minus_specificity,
                            plus_specificity, efficiency)
        for product in product_list:
            oligomer = "A%sB%s" % (product.count("A"), product.count("B"))
            profile_list.append(oligomer)
            if oligomer not in oligomer_list:
                oligomer_list.append(oligomer)
    # count products per composition
    length_profile_list = len(profile_list)
    output_list = []
    oligomer_list.sort()
    oligomer_list.sort(key=calc_DP)
    proportion_list = []
    for i in oligomer_list:
        proportion = profile_list.count(i) / length_profile_list
        proportion_list.append(proportion)
    output_list = [oligomer_list, proportion_list]
    return(output_list)

    
@anvil.server.callable
def profile_FAp_molar_nano2(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules, cuts):
    # generate library
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and save product compositions
    profile_list = []
    oligomer_list = []
    for substrate in library:
        product_list = nano2(substrate, cuts)
        for product in product_list:
            oligomer = "A%sB%s" % (product.count("A"), product.count("B"))
            profile_list.append(oligomer)
            if oligomer not in oligomer_list:
                oligomer_list.append(oligomer)
    # count products per composition
    length_profile_list = len(profile_list)
    output_list = []
    oligomer_list.sort()
    oligomer_list.sort(key=calc_DP)
    proportion_list = []
    for i in oligomer_list:
        proportion = profile_list.count(i) / length_profile_list
        proportion_list.append(proportion)
    output_list = [oligomer_list, proportion_list]
    return(output_list)


@anvil.server.callable
def profile_FAp_weight_enzyme(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules, minus_specificity, plus_specificity, efficiency, mass_A, mass_B):
    # generate library
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and save product compositions
    profile_list = []
    oligomer_list = []
    mass_all_products = 0
    for substrate in library:
        product_list = enzyme(substrate, minus_specificity,
                            plus_specificity, efficiency)
        for product in product_list:
            oligomer = "A%sB%s" % (product.count("A"), product.count("B"))
            mass_oligomer = product.count("A") * mass_A + product.count("B") * mass_B
            mass_all_products += mass_oligomer
            profile_list.append(oligomer)
            if oligomer not in oligomer_list:
                oligomer_list.append(oligomer)
    # calculate masses
    output_list = []
    oligomer_list.sort()
    oligomer_list.sort(key=calc_DP)
    proportion_list = []
    for oligomer in oligomer_list:
        number_A = int(re.search('(A)(\d+)', oligomer).group(2))
        number_B = int(re.search('(B)(\d+)', oligomer).group(2))
        mass_oligomer = number_A * mass_A + number_B * mass_B
        proportion = profile_list.count(oligomer) * mass_oligomer / mass_all_products
        proportion_list.append(proportion)
    output_list = [oligomer_list, proportion_list]
    return(output_list)

    
@anvil.server.callable
def profile_FAp_weight_nano2(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules, cuts, mass_A, mass_B):
    # generate library
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and save product compositions
    profile_list = []
    oligomer_list = []
    mass_all_products = 0
    for substrate in library:
        product_list = nano2(substrate, cuts)
        for product in product_list:
            oligomer = "A%sB%s" % (product.count("A"), product.count("B"))
            mass_oligomer = product.count("A") * mass_A + product.count("B") * mass_B
            mass_all_products += mass_oligomer
            profile_list.append(oligomer)
            if oligomer not in oligomer_list:
                oligomer_list.append(oligomer)
    # calculate masses
    output_list = []
    oligomer_list.sort()
    oligomer_list.sort(key=calc_DP)
    proportion_list = []
    for oligomer in oligomer_list:
        number_A = int(re.search('(A)(\d+)', oligomer).group(2))
        number_B = int(re.search('(B)(\d+)', oligomer).group(2))
        mass_oligomer = number_A * mass_A + number_B * mass_B
        proportion = profile_list.count(oligomer) * mass_oligomer / mass_all_products
        proportion_list.append(proportion)
    output_list = [oligomer_list, proportion_list]
    return(output_list)

    
@anvil.server.callable
def profile_strength_molar_enzyme(DP, overall_FA, pattern, strength, A_blocks, B_blocks, molecules, minus_specificity, plus_specificity, efficiency):
    # generate library
    FA_pattern = FA_pattern_calc(overall_FA, pattern, strength)
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and save product compositions
    profile_list = []
    oligomer_list = []
    for substrate in library:
        product_list = enzyme(substrate, minus_specificity,
                            plus_specificity, efficiency)
        for product in product_list:
            oligomer = "A%sB%s" % (product.count("A"), product.count("B"))
            profile_list.append(oligomer)
            if oligomer not in oligomer_list:
                oligomer_list.append(oligomer)
    # count products per composition
    length_profile_list = len(profile_list)
    output_list = []
    oligomer_list.sort()
    oligomer_list.sort(key=calc_DP)
    proportion_list = []
    for i in oligomer_list:
        proportion = profile_list.count(i) / length_profile_list
        proportion_list.append(proportion)
    output_list = [oligomer_list, proportion_list]
    return(output_list)


@anvil.server.callable
def profile_strength_molar_nano2(DP, overall_FA, pattern, strength, A_blocks, B_blocks, molecules, cuts):
    # generate library
    FA_pattern = FA_pattern_calc(overall_FA, pattern, strength)
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and save product compositions
    profile_list = []
    oligomer_list = []
    for substrate in library:
        product_list = nano2(substrate, cuts)
        for product in product_list:
            oligomer = "A%sB%s" % (product.count("A"), product.count("B"))
            profile_list.append(oligomer)
            if oligomer not in oligomer_list:
                oligomer_list.append(oligomer)
    # count products per composition
    length_profile_list = len(profile_list)
    output_list = []
    oligomer_list.sort()
    oligomer_list.sort(key=calc_DP)
    proportion_list = []
    for i in oligomer_list:
        proportion = profile_list.count(i) / length_profile_list
        proportion_list.append(proportion)
    output_list = [oligomer_list, proportion_list]
    return(output_list)


@anvil.server.callable
def profile_strength_weight_enzyme(DP, overall_FA, pattern, strength, A_blocks, B_blocks, molecules, minus_specificity, plus_specificity, efficiency, mass_A, mass_B):
    # generate library
    FA_pattern = FA_pattern_calc(overall_FA, pattern, strength)
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and save product compositions
    profile_list = []
    oligomer_list = []
    mass_all_products = 0
    for substrate in library:
        product_list = enzyme(substrate, minus_specificity,
                            plus_specificity, efficiency)
        for product in product_list:
            oligomer = "A%sB%s" % (product.count("A"), product.count("B"))
            mass_oligomer = product.count("A") * mass_A + product.count("B") * mass_B
            mass_all_products += mass_oligomer
            profile_list.append(oligomer)
            if oligomer not in oligomer_list:
                oligomer_list.append(oligomer)
    # calculate masses
    output_list = []
    oligomer_list.sort()
    oligomer_list.sort(key=calc_DP)
    proportion_list = []
    for oligomer in oligomer_list:
        number_A = int(re.search('(A)(\d+)', oligomer).group(2))
        number_B = int(re.search('(B)(\d+)', oligomer).group(2))
        mass_oligomer = number_A * mass_A + number_B * mass_B
        proportion = profile_list.count(oligomer) * mass_oligomer / mass_all_products
        proportion_list.append(proportion)
    output_list = [oligomer_list, proportion_list]
    return(output_list)

    
@anvil.server.callable
def profile_strength_weight_nano2(DP, overall_FA, pattern, strength, A_blocks, B_blocks, molecules, cuts, mass_A, mass_B):
    # generate library
    FA_pattern = FA_pattern_calc(overall_FA, pattern, strength)
    library = chitosan_library(DP, overall_FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
    # cleave and save product compositions
    profile_list = []
    oligomer_list = []
    mass_all_products = 0
    for substrate in library:
        product_list = nano2(substrate, cuts)
        for product in product_list:
            oligomer = "A%sB%s" % (product.count("A"), product.count("B"))
            mass_oligomer = product.count("A") * mass_A + product.count("B") * mass_B
            mass_all_products += mass_oligomer
            profile_list.append(oligomer)
            if oligomer not in oligomer_list:
                oligomer_list.append(oligomer)
    # calculate masses
    output_list = []
    oligomer_list.sort()
    oligomer_list.sort(key=calc_DP)
    proportion_list = []
    for oligomer in oligomer_list:
        number_A = int(re.search('(A)(\d+)', oligomer).group(2))
        number_B = int(re.search('(B)(\d+)', oligomer).group(2))
        mass_oligomer = number_A * mass_A + number_B * mass_B
        proportion = profile_list.count(oligomer) * mass_oligomer / mass_all_products
        proportion_list.append(proportion)
    output_list = [oligomer_list, proportion_list]
    return(output_list)

    
@anvil.server.callable
def make_profile(data, xlabel, file_name, DP_cutoff):
    data_w_cutoff = []
    legend_w_cutoff = []
    rest_counter = 0
    color_counter = 0
    for idx, oligomer in enumerate(data[0]):
        if calc_DP(oligomer) <= DP_cutoff:
            legend_w_cutoff.append(data[0][idx])
            data_w_cutoff.append(data[1][idx])
            color_counter += 1
        else:
            rest_counter += data[1][idx]
    legend_w_cutoff.append("> DP %s" % DP_cutoff)
    data_w_cutoff.append(rest_counter)
    df = pd.DataFrame(data_w_cutoff)
    colormap = ['#264653', '#2A9D8F', '#6ED8CC', '#E9C46A', '#F4A261', '#E76F51', '#BC4749', '#A7C957',
              '#6A994E', '#386641', '#24422A', '#800F2F', '#C9184A', '#EE6D91', '#F4A4BB', '#4D386B', '#896CB2', '#B4A1CE']
    times_colormap = math.floor(color_counter / len(colormap))
    additional_colors = color_counter % len(colormap)
    colors = times_colormap * colormap + colormap[0:additional_colors]
    colors.append('grey')
    df.T.plot(stacked=True, kind='barh', figsize=(20,5), color=colors)
    plt.legend(legend_w_cutoff, loc="upper center", ncol=9, fontsize=16)
    ax = plt.axes()
    ax.get_yaxis().set_visible(False)
    plt.xlabel(xlabel, fontsize=28)
    plt.xticks(fontsize=22)
    plt.xlim(0,1)
    plt.ylim(-0.5,2)
    plt.subplots_adjust(bottom=0.2)
    # Return this plot as a PNG image in a Media object
    return anvil.mpl_util.plot_image(filename=("%s.png" % file_name))

    
@anvil.server.callable
def profile_csv(data, ylabel, filename):
    df = pd.DataFrame(data).transpose()
    df.columns = ["oligomer", ylabel]
    df.to_csv('/tmp/data.csv', index=False, encoding='utf-8-sig')
    csv_file = anvil.media.from_file('/tmp/data.csv', 'csv', ('%s.csv' % filename))
    return csv_file

#### NMR ####

@anvil.server.callable
def diads_triads(DP, overall_FA, A_blocks, B_blocks, molecules, strength, pattern):
    output_dict = {}
    FA_list = []
    FA_pattern_list = []
    F_AA_list = []
    F_AB_list = []
    F_BA_list = []
    F_BB_list = []
    PA_list = []
    F_AAA_list = []
    F_AAB_list = []
    F_BAA_list = []
    F_ABA_list = []
    F_ABB_list = []
    F_BBA_list = []
    F_BAB_list = []
    F_BBB_list = []
    step = 0.02
    if A_blocks == 0 and B_blocks == 0:
        if overall_FA == None:
            FA_range = np.arange(step, 1, step)
        else:
            FA_range = [overall_FA]
    else:
        FA_range = [A_blocks/(A_blocks + B_blocks)]
    for FA in FA_range:
        FA_list.append(FA)
        FA_pattern = FA_pattern_calc(FA, pattern, strength)
        FA_pattern_list.append(FA_pattern)
        library = chitosan_library(DP, FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
        diads_PA_values = diads_PA(library)
        F_AA_list.append(diads_PA_values[0])
        F_AB_list.append(diads_PA_values[1])
        F_BA_list.append(diads_PA_values[2])
        F_BB_list.append(diads_PA_values[3])
        PA_list.append(diads_PA_values[4])
        triads_values = triads(library)
        F_AAA_list.append(triads_values[0])
        F_AAB_list.append(triads_values[1])
        F_BAA_list.append(triads_values[2])
        F_ABA_list.append(triads_values[3])
        F_ABB_list.append(triads_values[4])
        F_BBA_list.append(triads_values[5])
        F_BAB_list.append(triads_values[6])
        F_BBB_list.append(triads_values[7])    
    output_dict["average fraction unit A"] = FA_list
    output_dict["average fraction unit A in pattern"] = FA_pattern_list
    output_dict["F_AA"] = F_AA_list
    output_dict["F_AB"] = F_AB_list
    output_dict["F_BA"] = F_BA_list
    output_dict["F_BB"] = F_BB_list
    output_dict["PA"] = PA_list
    output_dict["F_AAA"] = F_AAA_list
    output_dict["F_AAB"] = F_AAB_list
    output_dict["F_BAA"] = F_BAA_list
    output_dict["F_ABA"] = F_ABA_list
    output_dict["F_ABB"] = F_ABB_list
    output_dict["F_BBA"] = F_BBA_list
    output_dict["F_BAB"] = F_BAB_list
    output_dict["F_BBB"] = F_BBB_list
    return(output_dict)

    
@anvil.server.callable
def make_PA_plot(data_nmr, overall_FA, filename):
    datax = data_nmr["average fraction unit A"]
    datay = data_nmr["PA"]
    plt.figure(1, figsize=(20,10), dpi=400)    
    ax = plt.axes()
    if overall_FA != None:
        plt.scatter(datax, datay, marker="o", s=120, c="#009688")
    else:
        plt.plot(datax, datay, color="#009688", linewidth=5.0)
    plt.xlabel("average fraction unit A", fontsize=28)
    plt.ylabel("PA", fontsize=28)
    plt.yticks(np.arange(0, 2.2, 0.2), fontsize=22)
    plt.xticks(np.arange(0, 1.1, 0.1), fontsize=22)
    plt.xlim(0,1)
    plt.ylim(-0.1,2.1)
    # Return this plot as a PNG image in a Media object
    return anvil.mpl_util.plot_image(filename=("%s.png" % filename))


@anvil.server.callable
def nmr_csv(data_nmr, filename):
    df = pd.DataFrame(data_nmr)
    df.to_csv('/tmp/data.csv', index=False, encoding='utf-8-sig', columns=["average fraction unit A", "average fraction unit A in pattern", "PA", "F_AA", "F_AB", "F_BA", "F_BB",
                                                                           "F_AAA", "F_AAB", "F_BAA", "F_ABA", "F_ABB", "F_BBA", "F_BAB", "F_BBB"])
    csv_file = anvil.media.from_file('/tmp/data.csv', 'csv', ('%s.csv' % filename))
    return csv_file


#### BLOCK SIZES ####

@anvil.server.callable
def block_sizes(DP, overall_FA, A_blocks, B_blocks, molecules, strength, pattern,
                efficiency, A_cutoff, DP_cutoff, ion_eff):
    output_dict = {}
    FA_list = []
    FA_pattern_list = []
    An_list = []
    Bn_list = []
    Aw_list = []
    Bw_list = []
    step = 0.02
    if A_blocks == 0 and B_blocks == 0:
        if overall_FA == None:
            FA_range = np.arange(step, 1, step)
        else:
            FA_range = [overall_FA]
    else:
        FA_range = [A_blocks/(A_blocks + B_blocks)]
    for FA in FA_range:
        FA_list.append(FA)
        FA_pattern = FA_pattern_calc(FA, pattern, strength)
        FA_pattern_list.append(FA_pattern)
        library = chitosan_library(DP, FA, pattern, FA_pattern, A_blocks, B_blocks, molecules)
        all_products = []
        for substrate in library:
            product_list = enzyme(substrate, ".BA", "XX.", efficiency)
            for product in product_list:
                all_products.append(product)
                all_products.append(product)
                all_products.append(product)
        output = blocks(all_products, A_cutoff, DP_cutoff, ion_eff)
        An_list.append(output[0])
        Bn_list.append(output[1])
        Aw_list.append(output[2])
        Bw_list.append(output[3])
    output_dict["average fraction unit A"] = FA_list
    output_dict["average fraction unit A in pattern"] = FA_pattern_list
    output_dict["block(A)n"] = An_list
    output_dict["block(B)n"] = Bn_list
    output_dict["block(A)w"] = Aw_list
    output_dict["block(B)w"] = Bw_list
    return(output_dict)

    
@anvil.server.callable
def make_blocks_plot(data_blocks, overall_FA, filename):
    datax = data_blocks["average fraction unit A"]
    data_An = data_blocks["block(A)n"]
    data_Bn = data_blocks["block(B)n"]
    data_Aw = data_blocks["block(A)w"]
    data_Bw = data_blocks["block(B)w"]
    plt.figure(1, figsize=(20,10), dpi=400)    
    ax = plt.axes()
    if overall_FA != None:
        plt.scatter(datax, data_Bn, marker="o", s=120, c='#E76F51', label = "block(B)n")
        plt.scatter(datax, data_Bw, marker="o", s=120, c='#E9C46A', label = "block(B)w")
        plt.scatter(datax, data_An, marker="o", s=120, c='#264653', label = "block(A)n")
        plt.scatter(datax, data_Aw, marker="o", s=120, c='#2A9D8F', label = "block(A)w")
    else:
        plt.plot(datax, data_Bn, color='#E76F51', label = "block(B)n", linewidth=5.0)
        plt.plot(datax, data_Bw, color='#E9C46A', label = "block(B)w", linewidth=5.0)
        plt.plot(datax, data_An, color='#264653', label = "block(A)n", linewidth=5.0)
        plt.plot(datax, data_Aw, color='#2A9D8F', label = "block(A)w", linewidth=5.0)
    plt.legend(loc="upper center", ncol=2, fontsize=25)
    plt.xlabel("average fraction unit A", fontsize=28)
    plt.ylabel("average block size", fontsize=28)
    plt.yticks(np.arange(0, 22, 2), fontsize=22)
    plt.xticks(np.arange(0, 1.1, 0.1), fontsize=22)
    plt.xlim(0,1)
    plt.ylim(-1,21)
    # Return this plot as a PNG image in a Media object
    return anvil.mpl_util.plot_image(filename=("%s.png" % filename))

    
@anvil.server.callable
def blocks_csv(data_blocks, filename):
    df = pd.DataFrame(data_blocks)
    df.to_csv('/tmp/data.csv', index=False, encoding='utf-8-sig', columns=["average fraction unit A", "average fraction unit A in pattern", "block(A)n", "block(B)n",
                                                                          "block(A)w", "block(B)w"])
    csv_file = anvil.media.from_file('/tmp/data.csv', 'csv', ('%s.csv' % filename))
    return csv_file
