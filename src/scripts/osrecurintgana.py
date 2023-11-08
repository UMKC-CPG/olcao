#!/usr/bin/env python3

"""
Collection of subroutines that are used to (1) produce analytical solutions to
multicenter gaussian type orbital integrals; (2) print Fortran code that, when
executed, will evaluate the analytical solutions for use in *either* a testing
or a production environment.
"""

import osrecurintglib as lib


def triad_search(triad, triads, increment, xyz_index, min_idx, max_idx):

    # Increment or decrement the appropriate cartesian angular momentum index.
    if (increment):
        triad[xyz_index] += 1
    else:
        triad[xyz_index] -= 1

    # Search the triads in the given range for a match of the modified triad.
    for j in range(min_idx, max_idx):
        if (triad == triads[j]):
            # Return the index where the match was found.
            return j

    # In the event that no match is found, then return -1.
    return -1


# Used by multiple recursive formulas.
def add_recursive_minus_terms(idx, icode, a, b, i, a_minus_idx, b_minus_idx,
                              triads, preFactorTag, vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Initialize the temp string so that we can build on it. If either minus
    #   indices are != -1 then we prefix the whole expression with inv_2zeta.
    if ((a_minus_idx != -1) or (b_minus_idx != -1)):
        temp_string = f" + inv_2zeta{vec_tag}*("
    else:
        temp_string = ""

    if (a_minus_idx != -1):
        if (triads[a][i] > 1):
            temp_string += f"{triads[a][i]}*"
        if (a_minus_idx+1 == 1 and b+1 == 1):
            if (icode == 4):  # Nuclear with two auxiliary integrals.
                temp_string += \
                        f"(preFactor{preFactorTag}({vec_tag_comma}{idx+1}) - "
                temp_string += \
                        f"preFactor{preFactorTag}({vec_tag_comma}{idx+2}))"
            elif ((icode == 5) or (icode == 7)):  # Momentum, dipole mom.
                temp_string += \
                        f"preFactor{preFactorTag}({vec_tag_comma}{idx+1})"
            else:  # Overlap, KE, Elec. 3C, Mass Vel.
                temp_string += f"preFactor{preFactorTag}{vec_tag}"
        else:
            if (icode == 4):  # (Nuclear with two auxiliary integrals.)
                temp_string += \
                        f"(pc({vec_tag_comma}{a_minus_idx+1},{b+1},{idx+1}) - "
                temp_string += \
                        f"pc({vec_tag_comma}{a_minus_idx+1},{b+1},{idx+2}))"
            elif ((icode == 5) or (icode == 7)):  # Momentum, dipole mom.
                temp_string += \
                        f"pc({vec_tag_comma}{a_minus_idx+1},{b+1},{idx+1})"
            else:  # Overlap, KE, Elec. 3C, Mass Vel.
                temp_string += f"pc({vec_tag_comma}{a_minus_idx+1},{b+1})"

    if ((a_minus_idx != -1) and (b_minus_idx != -1)):
        temp_string += " + "

    if (b_minus_idx != -1):
        if (triads[b][i] > 1):
            temp_string += f"{triads[b][i]}*"
        if (a+1 == 1 and b_minus_idx+1 == 1):
            if (icode == 4):  # Nuclear with two auxiliary integrals.
                temp_string += \
                        f"(preFactor{preFactorTag}({vec_tag_comma}{idx+1}) - "
                temp_string += \
                        f"preFactor{preFactorTag}({vec_tag_comma}{idx+2}))"
            elif ((icode == 5) or (icode == 7)): # Momentum or dipole
                temp_string += \
                        f"preFactor{preFactorTag}({vec_tag_comma}{idx+1})"
            else: # Everyone else
                temp_string += f"preFactor{preFactorTag}{vec_tag}"
        else:
            if (icode == 4): # Nuclear with two auxiliary integrals
                temp_string += \
                        f"(pc({vec_tag_comma}{a+1},{b_minus_idx+1},{idx+1}) - "
                temp_string += \
                        f"pc({vec_tag_comma}{a+1},{b_minus_idx+1},{idx+2}))"
            elif ((icode == 5) or (icode == 7)): # Momentum, dipole mom.
                temp_string += \
                        f"pc({vec_tag_comma}{a+1},{b_minus_idx+1},{idx+1})"
            else: # Everyone else
                temp_string += \
                        f"pc({vec_tag_comma}{a+1},{b_minus_idx+1})"

#    if (a_minus_idx != -1):
#        if (triads[a][i] > 1):
#            temp_string += f"{triads[a][i]}*"
#        if (a_minus_idx+1 == 1 and b+1 == 1):
#            if (m >= 0): # Account for multiple auxiliary integrals. (Nuclear)
#                temp_string += f"(preFactorN({vec_tag_comma}{m+1}) - "
#                temp_string += f"preFactorN({vec_tag_comma}{m+2}))"
#            elif (xyz >= 0): # Momentum xyz
#                temp_string += f"preFactorMM({vec_tag_comma}{xyz+1})"
#            else: # Everyone else
#                temp_string += f"preFactor{preFactorTag}{vec_tag}"
#        else:
#            if (m >= 0): # Account for multiple auxiliary integrals. (Nuclear)
#                temp_string += \
#                        f"(pc({vec_tag_comma}{a_minus_idx+1},{b+1},{m+1}) - "
#                temp_string += \
#                        f"pc({vec_tag_comma}{a_minus_idx+1},{b+1},{m+2}))"
#            elif (xyz >= 0): # Momentum xyz
#                temp_string += \
#                        f"pc({vec_tag_comma}{a_minus_idx+1},{b+1},{xyz+1})"
#            else: # Everyone else
#                temp_string += f"pc({vec_tag_comma}{a_minus_idx+1},{b+1})"
#
#    if ((a_minus_idx != -1) and (b_minus_idx != -1)):
#        temp_string += " + "
#
#    if (b_minus_idx != -1):
#        if (triads[b][i] > 1):
#            temp_string += f"{triads[b][i]}*"
#        if (a+1 == 1 and b_minus_idx+1 == 1):
#            if (m >= 0): # Account for multiple auxiliary integrals. (Nuclear)
#                temp_string += f"(preFactorN({vec_tag_comma}{m+1}) - "
#                temp_string += f"preFactorN({vec_tag_comma}{m+2}))"
#            elif (xyz >= 0): # Momentum xyz
#                temp_string += f"preFactorMM({vec_tag_comma}{xyz+1})"
#            else: # Everyone else
#                temp_string += f"preFactor{preFactorTag}{vec_tag}"
#        else:
#            if (m >= 0): # Account for multiple auxiliary integrals. (Nuclear)
#                temp_string += \
#                        f"(pc({vec_tag_comma}{a+1},{b_minus_idx+1},{m+1}) - "
#                temp_string += \
#                        f"pc({vec_tag_comma}{a+1},{b_minus_idx+1},{m+2}))"
#            elif (xyz >= 0): # Momentum xyz
#                temp_string += \
#                        f"pc({vec_tag_comma}{a+1},{b_minus_idx+1},{xyz+1})"
#            else: # Everyone else
#                temp_string += \
#                        f"pc({vec_tag_comma}{a+1},{b_minus_idx+1})"

    # Similarly to the beginning, if either minus index is != -1 then we need
    #   to close the parentheses.
    if (a_minus_idx != -1 or b_minus_idx != -1):
        temp_string += ")"

    return temp_string


# Used specifically by the kinetic energy recursive formula.
def add_kinetic_overlap_terms(a, b, i, minus_idx, plus_idx, triads, do_a,
                              vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Always wrap the terms with this initial string.
    temp_string = f" + 2*xi{vec_tag}*("

    # Always add the a+1_i term.
    if (do_a):
        temp_string += f"pc_ol({vec_tag_comma}{plus_idx+1},{b+1})"
    else:
        temp_string += f"pc_ol({vec_tag_comma}{a+1},{plus_idx+1})"

    # Add the a-1_i term if the a-1_i index is not less than zero.
    if (minus_idx != -1):
        if (do_a):
            temp_string += f" - inv_2zeta_a{vec_tag}*"
            if (triads[a][i] > 1):
                temp_string += f"{triads[a][i]}*"
            if (minus_idx+1 == 1 and b+1 == 1):
                temp_string += f"preFactorOL{vec_tag}"
            else:
                temp_string += f"pc_ol({vec_tag_comma}{minus_idx+1},{b+1})"
        else:
            temp_string += f" - inv_2zeta_b{vec_tag}*"
            if (triads[b][i] > 1):
                temp_string += f"{triads[b][i]}*"
            if (a+1 == 1 and minus_idx+1 == 1):
                temp_string += f"preFactorOL{vec_tag}"
            else:
                temp_string += f"pc_ol({vec_tag_comma}{a+1},{minus_idx+1})"
    
    # Close the outer wrap.
    temp_string += ")"

    return temp_string


# Used specifically by the momentum recursive formula.
def add_momentum_overlap_terms(a, b, do_a, vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Add the appropriate prefix term.
    if (do_a):
        temp_string = f" - zeta_b_zeta{vec_tag}*"
    else:
        temp_string = f" + zeta_a_zeta{vec_tag}*"

    # Add the overlap term (which is the same regardless of a or b.
    if (a+1 == 1 and b+1 == 1):
        temp_string += f"preFactorOL{vec_tag}"
    else:
        temp_string += f"pc_ol({vec_tag_comma}{a+1},{b+1})"

    return temp_string


# Used specifically by the mass velocity recursive formula.
def add_massvel_overlap_terms(a, b, i, minus_idx, plus_idx, triads, do_a,
                              vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Always wrap the terms with this initial string.
    temp_string = f" + 16*xi{vec_tag}*inv_8m3c2*("

    # Always add the a+1_i term.
    if (do_a):
        temp_string += f"pc_ke({vec_tag_comma}{plus_idx+1},{b+1})"
    else:
        temp_string += f"pc_ke({vec_tag_comma}{a+1},{plus_idx+1})"

    # Add the a-1_i term if the a-1_i index is not less than zero.
    if (minus_idx != -1):
        if (do_a):
            temp_string += f" - inv_2zeta_a{vec_tag}*"
            if (triads[a][i] > 1):
                temp_string += f"{triads[a][i]}*"
            if (minus_idx+1 == 1 and b+1 == 1):
                temp_string += f"preFactorKE{vec_tag}"
            else:
                temp_string += f"pc_ke({vec_tag_comma}{minus_idx+1},{b+1})"
        else:
            temp_string += f" - inv_2zeta_b{vec_tag}*"
            if (triads[b][i] > 1):
                temp_string += f"{triads[b][i]}*"
            if (a+1 == 1 and minus_idx+1 == 1):
                temp_string += f"preFactorKE{vec_tag}"
            else:
                temp_string += f"pc_ke({vec_tag_comma}{a+1},{minus_idx+1})"
    
    # Close the outer wrap.
    temp_string += ")"

    return temp_string

# Used specifically by the dipole moment recursive formula.
def add_dipole_overlap_terms(a, b, i, xyz, triads, do_a, vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Dipole moment
    N_i_mu = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        
    # Always wrap the terms with this initial string.
    temp_string = f" + inv_2zeta{vec_tag}*("

    if (do_a):
        temp_string += f"{N_i_mu[xyz][i]}*"
        if (a+1 == 1 and b+1 == 1):
            temp_string += f"preFactorOL{vec_tag}"
        else:
            temp_string += f"pc_ol({vec_tag_comma}{a+1},{b+1})"
    else:
        temp_string += f"{N_i_mu[xyz][i]}*"
        if (a+1 == 1 and b+1 == 1):
            temp_string += f"preFactorOL{vec_tag}"
        else:
            temp_string += f"pc_ol({vec_tag_comma}{a+1},{b+1})"
    
    # Close the outer wrap.
    temp_string += ")"

    return temp_string


# This subroutine is used for the 2-center overlap, the kinetic energy, the
#   3-center overlap (used ultimately for electron repulsion after a charge
#   fitting process), and also the momentum matrix elements.
def overlap(triads, add_matrix_tag, relabel3C, vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Make a convenient shorthand for the length of the triad array.
    num_triads = len(triads)

    # Initialize the solution matrix to empty strings.
    soln_mtx_ol = [[""]*num_triads for i in range(num_triads)]

    # Initialize the 0,0 (python) 1,1 (fotran) element.
    soln_mtx_ol[0][0] = f"preFactorOL{vec_tag}"

    # Produce a string for each element of the solution matrix.
    for a in range(num_triads):
        for b in range(num_triads):

            # Every element that we visit in the a,b nested loops should
            #   already be filled in. But, we check here just to confirm.
            if (soln_mtx_ol[a][b] == ""):
                print("Missing element: ", a, ",", b, ".")
                print(soln_mtx_ol)
                exit()

            # Use the current element to build out higher elements. We will
            #   add 1 to each of the x, y, z components of the current 'a'
            #   and 'b' triads to build new element possibilities. We will
            #   then search the triad list for a match and if found we will
            #   assemble the element accordingly.

            # Start by adding 1 to the 'a' triad x, y, z terms and searching
            #   for matches in the other (higher) triads.
            for i in range(3):

                # Search for the index numbers in the triad list associated
                #   with a+1_i, a-1_i, and b-1_i.

                # Search for the a+1_i first.
                a_plus_idx = triad_search(triads[a].copy(), triads, True, i,
                                          a+1, num_triads)

                # If no match was found, then go to the next i (x,y,z).
                if (a_plus_idx == -1):
                    break

                # If this solution matrix element is already done, then
                #   skip to the next i.
                if (soln_mtx_ol[a_plus_idx][b] != ""):
                    continue

                # For the current a+1_i we need to find the associated a-1_i.
                a_minus_idx = \
                        triad_search(triads[a].copy(), triads, False, i, 0, a)

                # For the current a+1_i we need to find the associated b-1_i.
                b_minus_idx = \
                        triad_search(triads[b].copy(), triads, False, i, 0, b)

                # At this point we have an a+1_i match and have found the
                #   indices for a-1_i and b-1_i (or the indices are -1 if
                #   no match was found).

                # Now, we assemble the [a_plus_idx][b] element.
                if (a+1 == 1 and b+1 == 1):
                    temp_string = f"PA({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorOL{vec_tag}"
                else:
                    temp_string = f"PA({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1})"

                # Add the overlap terms from lower cart. ang. mom.
                temp_string += add_recursive_minus_terms(1, 1, a, b, i,
                        a_minus_idx, b_minus_idx, triads, "OL", vectorize)

                soln_mtx_ol[a_plus_idx][b] = temp_string

            # Now do the same for the b triad. Start by adding a 1 to the 'b'
            #   triad x, y, z terms and searching for matches in the other
            #   higher triads.
            for i in range(3):

                # Search for the index numbers in the triad list associated
                #   with b+1_i, b-1_i, and a-1_i.
                
                # Search for the b+1_i first.
                b_plus_idx = triad_search(triads[b].copy(), triads, True, i,
                                          b+1, num_triads)

                # If no match was found, then go to the next i (x,y,z).
                if (b_plus_idx == -1):
                    break

                # If this solution matrix element is already done, then
                #   skip to the next i.
                if (soln_mtx_ol[a][b_plus_idx] != ""):
                    continue

                # For the current b+1_i we need to find the associated b-1_i.
                b_minus_idx = \
                        triad_search(triads[b].copy(), triads, False, i, 0, b)

                # For the current b+1_i we need to find the associated a-1_i.
                a_minus_idx = \
                        triad_search(triads[a].copy(), triads, False, i, 0, a)

                # At this point we have a b+1_i match and have found the
                #   indices for b-1_i and a-1_i (or the indices are -1 if
                #   no match was found.

                # Now we assemble the [a][b_plus_idx] element.
                if (a+1 == 1 and b+1 == 1):
                    temp_string = f"PB({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorOL{vec_tag}"
                else:
                    temp_string = f"PB({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1})"

                # Add the overlap terms from lower cart. ang. mom.
                temp_string += add_recursive_minus_terms(1, 1, a, b, i,
                        a_minus_idx, b_minus_idx, triads, "OL", vectorize)

                soln_mtx_ol[a][b_plus_idx] = temp_string

    # The overlap procedure is customizable for use in other integrals. Here
    #   we make name changes to accomodate those other uses.

    # If the overlap is explicitly used (as for the kinetic energy) then we
    #   relabel the pc terms with pc_ol to clearly identify them as being
    #   specifically from the overlap.
    if (add_matrix_tag):
        # Replace all instances of "pc" from the overlap solution matrix with
        #   "pc_ol". Thus, the overlap solution elements can be integrated
        #   with the kinetic energy (or other) solution elements.
        soln_mtx_ol = \
                [[y.replace("pc", "pc_ol") for y in x] for x in soln_mtx_ol]

    # For the three center overlap we relabel the PA and PB variables to match
    #   with the notation used in Obara Saika (1986). Also relabel inv_2zeta
    #   with inv_2zeta3C.
    if (relabel3C):
        soln_mtx_ol = \
                [[y.replace("PA", "GA") for y in x] for x in soln_mtx_ol]
        soln_mtx_ol = \
                [[y.replace("PB", "GB") for y in x] for x in soln_mtx_ol]
        soln_mtx_ol = \
                [[y.replace("inv_2zeta", "inv_2zeta3C") for y in x] \
                for x in soln_mtx_ol]

    return soln_mtx_ol


# Relies on the overlap solutions.
def kinetic(triads, vectorize, self_reference):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Make a convenient shorthand for the length of the triad array.
    num_triads = len(triads)

    # Initialize the solution matrix to empty strings.
    soln_mtx_ke = [[""]*num_triads for i in range(num_triads)]

    # Initialize the 0,0 (python) 1,1 (fotran) elements.
    soln_mtx_ke[0][0] = f"preFactorKE{vec_tag}"

    # Produce a string for each element of the solution matrix.
    for a in range(num_triads):
        for b in range(num_triads):

            # Every element that we visit in the a,b nested loops should
            #   already be filled in. But, we check here just to confirm.
            if (soln_mtx_ke[a][b] == ""):
                print("Missing element: ", a, ",", b, ".")
                print(soln_mtx_ke)
                exit()

            # Use the current element to build out higher elements. We will
            #   add 1 to each of the x, y, z components of the current 'a'
            #   and 'b' triads to build new element possibilities. We will
            #   then search the triad list for a match and if found we will
            #   assemble the element accordingly.

            # Start by adding 1 to the 'a' triad x, y, z terms and searching
            #   for matches in the other (higher) triads.
            for i in range(3):

                # Search for the index numbers in the triad list associated
                #   with a+1_i, a-1_i, and b-1_i.

                # Search for the a+1_i first.
                a_plus_idx = triad_search(triads[a].copy(), triads, True, i,
                                          a+1, num_triads)

                # If no match was found, then go to the next i (x,y,z).
                if (a_plus_idx == -1):
                    break

                # If this solution matrix element is already done, then
                #   skip to the next i.
                if (soln_mtx_ke[a_plus_idx][b] != ""):
                    continue

                # For the current a+1_i we need to find the associated a-1_i.
                a_minus_idx = \
                        triad_search(triads[a].copy(), triads, False, i, 0, a)

                # For the current a+1_i we need to find the associated b-1_i.
                b_minus_idx = \
                        triad_search(triads[b].copy(), triads, False, i, 0, b)

                # At this point we have an a+1_i match and have found the
                #   indices for a-1_i and b-1_i (or the indices are -1 if
                #   no match was found).

                # Now, we assemble the [a_plus_idx][b] element.
                if (a+1 == 1 and b+1 == 1):
                    temp_string = f"PA({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorKE{vec_tag}"
                else:
                    temp_string = f"PA({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1})"

                # Add the first set of terms from lower cart. ang. momentum.
                #   This process is exactly the same as the overlap and so
                #   we just reuse this recursive code.
                temp_string += add_recursive_minus_terms(1, 2, a, b, i,
                        a_minus_idx, b_minus_idx, triads, "KE", vectorize)

                # Add the second set of terms (both higher and lower angular
                #   momentum) that explicitly contain overlap integrals.
                temp_string += add_kinetic_overlap_terms(a, b, i, a_minus_idx,
                        a_plus_idx, triads, True, vectorize)

                soln_mtx_ke[a_plus_idx][b] = temp_string

            # Now do the same for the b triad. Start by adding a 1 to the 'b'
            #   triad x, y, z terms and searching for matches in the other
            #   higher triads.
            for i in range(3):

                # Search for the index numbers in the triad list associated
                #   with b+1_i, b-1_i, and a-1_i.
                
                # Search for the b+1_i first.
                b_plus_idx = triad_search(triads[b].copy(), triads, True, i,
                                          b+1, num_triads)

                # If no match was found, then go to the next i (x,y,z).
                if (b_plus_idx == -1):
                    break

                # If this solution matrix element is already done, then
                #   skip to the next i.
                if (soln_mtx_ke[a][b_plus_idx] != ""):
                    continue

                # For the current b+1_i we need to find the associated b-1_i.
                b_minus_idx = \
                        triad_search(triads[b].copy(), triads, False, i, 0, b)

                # For the current b+1_i we need to find the associated a-1_i.
                a_minus_idx = \
                        triad_search(triads[a].copy(), triads, False, i, 0, a)

                # At this point we have a b+1_i match and have found the
                #   indices for b-1_i and a-1_i (or the indices are -1 if
                #   no match was found.

                # Now we assemble the [a][b_plus_idx] element.
                if (a+1 == 1 and b+1 == 1):
                    temp_string = f"PB({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorKE{vec_tag}"
                else:
                    temp_string = f"PB({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1})"

                # Add the first set of terms from lower cart. ang. momentum.
                #   This process is exactly the same as the overlap and so
                #   we just reuse this recursive code.
                temp_string += add_recursive_minus_terms(1, 2, a, b, i,
                        a_minus_idx, b_minus_idx, triads, "KE", vectorize)

                # Add the second set of terms (both higher and lower angular
                #   momentum) that explicitly contain overlap integrals.
                temp_string += add_kinetic_overlap_terms(a, b, i, b_minus_idx,
                        b_plus_idx, triads, False, vectorize)

                soln_mtx_ke[a][b_plus_idx] = temp_string

    # If the kinetic energy solutions is used to compute another matrix, then
    #   we need to change the references to the recursive pc elements.
    if (self_reference):
        soln_mtx_ke = \
                [[y.replace("pc(", "pc_ke(") for y in x] \
                for x in soln_mtx_ke]

    return soln_mtx_ke


def nuclear(triads, m, vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Make a convenient shorthand for the length of the triad array.
    num_triads = len(triads)

    # Initialize the nuclear potential solution matrix to empty strings.
    soln_mtx_np = [[""]*num_triads for i in range(num_triads)] 

    # Initialize the 0,0 (python) 1,1 (fotran) element.
    soln_mtx_np[0][0] = f"preFactorN({vec_tag_comma}{m+1})"

    # Produce a string for each element of the solution matrix.
    for a in range(num_triads):
        for b in range(num_triads):

            # Every element that we visit in the a,b nested loops should
            #   already be filled in. But, we check here just to confirm.
            if (soln_mtx_np[a][b] == ""):
                print("Missing element: ", a, ",", b, ".")
                print(soln_mtx_np)
                exit()

            # Use the current element to build out higher elements. We will
            #   add 1 to each of the x, y, z components of the current 'a'
            #   and 'b' triads to build new element possibilities. We will
            #   then search the triad list for a match and if found we will
            #   assemble the element accordingly.

            # Start by adding 1 to the 'a' triad x, y, z terms and searching
            #   for matches in the other (higher) triads.
            for i in range(3):

                # Search for the index numbers in the triad list associated
                #   with a+1_i, a-1_i, and b-1_i.

                # Search for the a+1_i first.
                a_plus_idx = triad_search(triads[a].copy(), triads, True, i,
                                          a+1, num_triads)

                # If no match was found, then go to the next i (x,y,z).
                if (a_plus_idx == -1):
                    break

                # If this solution matrix element is already done, then
                #   skip to the next i.
                if (soln_mtx_np[a_plus_idx][b] != ""):
                    continue

                # For the current a+1_i we need to find the associated a-1_i.
                a_minus_idx = \
                        triad_search(triads[a].copy(), triads, False, i, 0, a)

                # For the current a+1_i we need to find the associated b-1_i.
                b_minus_idx = \
                        triad_search(triads[b].copy(), triads, False, i, 0, b)

                # At this point we have an a+1_i match and have found the
                #   indices for a-1_i and b-1_i (or the indices are -1 if
                #   no match was found).

                # Now, we assemble the [m][a_plus_idx][b] element.
                if (a+1 == 1 and b+1 == 1):
                    temp_string = f"GA({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorN({vec_tag_comma}{m+1}) - "
                    temp_string += f"GC({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorN({vec_tag_comma}{m+2})"
                else:
                    temp_string = f"GA({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1},{m+1}) - "
                    temp_string += f"GC({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1},{m+2})"

                # Add the overlap terms from lower cart. ang. mom.
                temp_string += add_recursive_minus_terms(m, 4, a, b, i,
                        a_minus_idx, b_minus_idx, triads, "N", vectorize)

                soln_mtx_np[a_plus_idx][b] = temp_string

            # Now do the same for the b triad. Start by adding a 1 to the 'b'
            #   triad x, y, z terms and searching for matches in the other
            #   higher triads.
            for i in range(3):

                # Search for the index numbers in the triad list associated
                #   with b+1_i, b-1_i, and a-1_i.

                # Search for the b+1_i first.
                b_plus_idx = triad_search(triads[b].copy(), triads, True, i,
                                          b+1, num_triads)

                # If no match was found, then go to the next i (x,y,z).
                if (b_plus_idx == -1):
                    break

                # If this solution matrix element is already done, then
                #   skip to the next i.
                if (soln_mtx_np[a][b_plus_idx] != ""):
                    continue

                # For the current b+1_i we need to find the associated b-1_i.
                b_minus_idx = \
                        triad_search(triads[b].copy(), triads, False, i, 0, b)

                # For the current b+1_i we need to find the associated a-1_i.
                a_minus_idx = \
                        triad_search(triads[a].copy(), triads, False, i, 0, a)

                # At this point we have a b+1_i match and have found the
                #   indices for b-1_i and a-1_i (or the indices are -1 if
                #   no match was found.

                # Now we assemble the [m][a][b_plus_idx] element.
                if (a+1 == 1 and b+1 == 1):
                    temp_string = f"GB({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorN({vec_tag_comma}{m+1}) - "
                    temp_string += f"GC({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorN({vec_tag_comma}{m+2})"
                else:
                    temp_string = f"GB({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1},{m+1}) - "
                    temp_string += f"GC({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1},{m+2})"

                # Add the overlap terms from lower cart. ang. mom.
                temp_string += add_recursive_minus_terms(m, 4, a, b, i,
                        a_minus_idx, b_minus_idx, triads, "N", vectorize)

                soln_mtx_np[a][b_plus_idx] = temp_string

    # For the three center nuclear we relabel the inv_2zeta variable for use
    #   in a three center integral (as inv_2zeta3C). Note that this differs
    #   from the Obara-Saika approach where the nuclear integral there does
    #   not include a third Gaussian.
    soln_mtx_np = \
            [[y.replace("inv_2zeta", "inv_2zeta3C") for y in x] \
            for x in soln_mtx_np]

    return soln_mtx_np


def momentum(triads, vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Make a convenient shorthand for the length of the triad array.
    num_triads = len(triads)

    # Initialize the momentum solution matrix to empty strings.
    soln_mtx_mm = [[[""]*3 for i in range(num_triads)] \
                   for j in range(num_triads)]

    # Initialize the 0,0,0..2 (python) 1,1,1..3 (fotran) elements.
    soln_mtx_mm[0][0][0] = f"preFactorMM({vec_tag_comma}1)"
    soln_mtx_mm[0][0][1] = f"preFactorMM({vec_tag_comma}2)"
    soln_mtx_mm[0][0][2] = f"preFactorMM({vec_tag_comma}3)"

    # Produce a string for each element of the solution matrix.
    for xyz in range(3):
        for a in range(num_triads):
            for b in range(num_triads):

                # Every element that we visit in the a,b nested loops should
                #   already be filled in. But, we check here just to confirm.
                if (soln_mtx_mm[a][b][xyz] == ""):
                    print("Missing element: ", a, ",", b, ".")
                    print(soln_mtx_mm)
                    exit()
    
                # Use the current element to build out higher elements. We
                #   add 1 to each of the x, y, z components of the current 'a'
                #   and 'b' triads to build new element possibilities. We will
                #   then search the triad list for a match and if found we
                #   assemble the element accordingly.
    
                # Start by adding 1 to the 'a' triad x, y, z terms and
                #   searching for matches in the other (higher) triads.
                for i in range(3):
    
                    # Search for the index numbers in the triad list
                    #   associated with a+1_i, a-1_i, and b-1_i.
    
                    # Search for the a+1_i first.
                    a_plus_idx = triad_search(triads[a].copy(), triads, True,
                                              i, a+1, num_triads)
    
                    # If no match was found, then go to the next i (x,y,z).
                    if (a_plus_idx == -1):
                        break
    
                    # If this solution matrix element is already done, then
                    #   skip to the next i.
                    if (soln_mtx_mm[a_plus_idx][b][xyz] != ""):
                        continue
    
                    # For the current a+1_i, find the associated a-1_i.
                    a_minus_idx = triad_search(triads[a].copy(), triads,
                                               False, i, 0, a)
    
                    # For the current a+1_i, find the associated b-1_i.
                    b_minus_idx = triad_search(triads[b].copy(), triads,
                                               False, i, 0, b)
    
                    # At this point we have an a+1_i match and have found the
                    #   indices for a-1_i and b-1_i (or the indices are -1 if
                    #   no match was found).
    
                    # Now, we assemble the [xyz][a_plus_idx][b] element.
                    if (a+1 == 1 and b+1 == 1):
                        temp_string = f"PA({vec_tag_comma}{i+1})*"
                        temp_string += f"preFactorMM({vec_tag_comma}{xyz+1})"
                    else:
                        temp_string = f"PA({vec_tag_comma}{i+1})*"
                        temp_string += f"pc({vec_tag_comma}{a+1},{b+1},{xyz+1})"
    
                    # Add the first set of terms from lower cart. ang. mom.
                    #   This process is exactly the same as the overlap and
                    #   so we just reuse this recursive code.
                    temp_string += add_recursive_minus_terms(xyz, 5, a, b, i,
                            a_minus_idx, b_minus_idx, triads, "MM", vectorize)

                    # Add the second term (which is of equal ang. mom.)
                    #   and that explicity contains an overlap integral.
                    if (xyz == i):
                        temp_string += add_momentum_overlap_terms(a, b, True,
                                vectorize)
    
                    soln_mtx_mm[a_plus_idx][b][xyz] = temp_string
    
                # Now do the same for the b triad. Start by adding a 1 to the
                #   'b' triad x, y, z terms and searching for matches in the
                #   other higher triads.
                for i in range(3):
    
                    # Search for the index numbers in the triad list
                    #   associated with b+1_i, b-1_i, and a-1_i.
                    
                    # Search for the b+1_i first.
                    b_plus_idx = triad_search(triads[b].copy(), triads, True,
                                              i, b+1, num_triads)
    
                    # If no match was found, then go to the next i (x,y,z).
                    if (b_plus_idx == -1):
                        break
    
                    # If this solution matrix element is already done, then
                    #   skip to the next i.
                    if (soln_mtx_mm[a][b_plus_idx][xyz] != ""):
                        continue
    
                    # For the current b+1_i, find the associated b-1_i.
                    b_minus_idx = triad_search(triads[b].copy(), triads,
                                               False, i, 0, b)
    
                    # For the current b+1_i, find the associated a-1_i.
                    a_minus_idx = triad_search(triads[a].copy(), triads,
                                               False, i, 0, a)
    
                    # At this point we have a b+1_i match and have found the
                    #   indices for b-1_i and a-1_i (or the indices are -1 if
                    #   no match was found.
    
                    # Now we assemble the [xyz][a][b_plus_idx] element.
                    if (a+1 == 1 and b+1 == 1):
                        temp_string = f"PB({vec_tag_comma}{i+1})*"
                        temp_string += f"preFactorMM({vec_tag_comma}{xyz+1})"
                    else:
                        temp_string = f"PB({vec_tag_comma}{i+1})*"
                        temp_string += f"pc({vec_tag_comma}{a+1},{b+1},{xyz+1})"
    
                    # Add the first set of terms from lower cart. ang. mom.
                    #   this process is exactly the same as the overlap and
                    #   so we just reuse this recursive code.
                    temp_string += add_recursive_minus_terms(xyz, 5, a, b, i,
                            a_minus_idx, b_minus_idx, triads, "MM", vectorize)
                    
                    # Add the second term (which is of equal ang. mom.)
                    #   and that explicity contains an overlap integral.
                    if (xyz == i):
                        temp_string += add_momentum_overlap_terms(a, b, False,
                                vectorize)
    
                    soln_mtx_mm[a][b_plus_idx][xyz] = temp_string

    return soln_mtx_mm


# Relies on the overlap and kinetic energy solutions.
# See the discussion about the nature of the mass velocity integral in the
#   osrecurintgnum.py file.
# Use the Obara-Saika recurrence relations given in Molecular electronic
#   structure theory. Equations 9.3.26 - 9.3.28 on page 348. The form looks
#   like Q_ij^4
def massvel(triads, vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Make a convenient shorthand for the length of the triad array.
    num_triads = len(triads)

    # Initialize the solution matrix to empty strings.
    soln_mtx_mv = [[""]*num_triads for i in range(num_triads)]

    # Initialize the 0,0 (python) 1,1 (fotran) elements.
    soln_mtx_mv[0][0] = f"preFactorMV{vec_tag}"

    # Produce a string for each element of the solution matrix.
    for a in range(num_triads):
        for b in range(num_triads):

            # Every element that we visit in the a,b nested loops should
            #   already be filled in. But, we check here just to confirm.
            if (soln_mtx_mv[a][b] == ""):
                print("Missing element: ", a, ",", b, ".")
                print(soln_mtx_mv)
                exit()

            # Use the current element to build out higher elements. We will
            #   add 1 to each of the x, y, z components of the current 'a'
            #   and 'b' triads to build new element possibilities. We will
            #   then search the triad list for a match and if found we will
            #   assemble the element accordingly.

            # Start by adding 1 to the 'a' triad x, y, z terms and searching
            #   for matches in the other (higher) triads.
            for i in range(3):

                # Search for the index numbers in the triad list associated
                #   with a+1_i, a-1_i, and b-1_i.

                # Search for the a+1_i first.
                a_plus_idx = triad_search(triads[a].copy(), triads, True, i,
                                          a+1, num_triads)

                # If no match was found, then go to the next i (x,y,z).
                if (a_plus_idx == -1):
                    break

                # If this solution matrix element is already done, then
                #   skip to the next i.
                if (soln_mtx_mv[a_plus_idx][b] != ""):
                    continue

                # For the current a+1_i we need to find the associated a-1_i.
                a_minus_idx = \
                        triad_search(triads[a].copy(), triads, False, i, 0, a)

                # For the current a+1_i we need to find the associated b-1_i.
                b_minus_idx = \
                        triad_search(triads[b].copy(), triads, False, i, 0, b)

                # At this point we have an a+1_i match and have found the
                #   indices for a-1_i and b-1_i (or the indices are -1 if
                #   no match was found).

                # Now, we assemble the [a_plus_idx][b] element.
                if (a+1 == 1 and b+1 == 1):
                    temp_string = f"PA({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorMV{vec_tag}"
                else:
                    temp_string = f"PA({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1})"

                # Add the first set of terms from lower cart. ang. momentum.
                #   This process is exactly the same as the overlap and so
                #   we just reuse this recursive code.
                temp_string += add_recursive_minus_terms(1, 6, a, b, i,
                        a_minus_idx, b_minus_idx, triads, "MV", vectorize)

                # Add the second set of terms (both higher and lower angular
                #   momentum) that explicitly contain overlap integrals.
                temp_string += add_massvel_overlap_terms(a, b, i, a_minus_idx,
                        a_plus_idx, triads, True, vectorize)

                soln_mtx_mv[a_plus_idx][b] = temp_string

            # Now do the same for the b triad. Start by adding a 1 to the 'b'
            #   triad x, y, z terms and searching for matches in the other
            #   higher triads.
            for i in range(3):

                # Search for the index numbers in the triad list associated
                #   with b+1_i, b-1_i, and a-1_i.
                
                # Search for the b+1_i first.
                b_plus_idx = triad_search(triads[b].copy(), triads, True, i,
                                          b+1, num_triads)

                # If no match was found, then go to the next i (x,y,z).
                if (b_plus_idx == -1):
                    break

                # If this solution matrix element is already done, then
                #   skip to the next i.
                if (soln_mtx_mv[a][b_plus_idx] != ""):
                    continue

                # For the current b+1_i we need to find the associated b-1_i.
                b_minus_idx = \
                        triad_search(triads[b].copy(), triads, False, i, 0, b)

                # For the current b+1_i we need to find the associated a-1_i.
                a_minus_idx = \
                        triad_search(triads[a].copy(), triads, False, i, 0, a)

                # At this point we have a b+1_i match and have found the
                #   indices for b-1_i and a-1_i (or the indices are -1 if
                #   no match was found.

                # Now we assemble the [a][b_plus_idx] element.
                if (a+1 == 1 and b+1 == 1):
                    temp_string = f"PB({vec_tag_comma}{i+1})*"
                    temp_string += f"preFactorMV{vec_tag}"
                else:
                    temp_string = f"PB({vec_tag_comma}{i+1})*"
                    temp_string += f"pc({vec_tag_comma}{a+1},{b+1})"

                # Add the first set of terms from lower cart. ang. momentum.
                #   This process is exactly the same as the overlap and so
                #   we just reuse this recursive code.
                temp_string += add_recursive_minus_terms(1, 6, a, b, i,
                        a_minus_idx, b_minus_idx, triads, "MV", vectorize)

                # Add the second set of terms (both higher and lower angular
                #   momentum) that explicitly contain overlap integrals.
                temp_string += add_massvel_overlap_terms(a, b, i, b_minus_idx,
                        b_plus_idx, triads, False, vectorize)

                soln_mtx_mv[a][b_plus_idx] = temp_string

    return soln_mtx_mv


# Relies on the overlap solutions.
def dipole(triads, vectorize):

    # If we need to vectorize, then add (:) or :, in the appropriate places.
    if (vectorize):
        vec_tag = "(:)"
        vec_tag_comma = ":,"
    else:
        vec_tag = ""
        vec_tag_comma = ""

    # Make a convenient shorthand for the length of the triad array.
    num_triads = len(triads)

    # Initialize the solution matrix to empty strings.
    #   soln_mtx_dm[a][b][xyz]
    soln_mtx_dm = [[[""]*3 for i in range(num_triads)] \
                   for j in range(num_triads)]

    # Initialize the 0,0,0..2,0..2 (python) 1,1,,1.3,1..3 (fotran) elements.
    # These respresent the (s|M|s) in the x-, y-, and z-direction of mu and 
    # the x-, y-. z-direction of A,B, and C
    soln_mtx_dm[0][0][0] = f"preFactorDM({vec_tag_comma}1)"
    soln_mtx_dm[0][0][1] = f"preFactorDM({vec_tag_comma}2)"
    soln_mtx_dm[0][0][2] = f"preFactorDM({vec_tag_comma}3)"


    # Produce a string for each element of the solution matrix.
    for xyz in range(3): # Go through the xyz directions
        for a in range(num_triads):
            for b in range(num_triads):

                # Every element that we visit in the a,b nested loops should
                #   already be filled in. But, we check here just to confirm.
                if (soln_mtx_dm[a][b][xyz] == ""):
                    print("Missing element: ", a, ",", b, ".")
                    print(soln_mtx_dm)
                    exit()

                # Use the current element to build out higher elements. We
                #   add 1 to each of the x, y, z components of the current 'a'
                #   and 'b' triads to build new element possibilities. We will
                #   then search the triad list for a match and if found we
                #   assemble the element accordingly.

                # Start by adding 1 to the 'a' triad x, y, z terms and
                #   searching for matches in the other (higher) triads.
                for i in range(3):

                    # Search for the index numbers in the triad list
                    #   associated with a+1_i, a-1_i, and b-1_i.

                    # Search for the a+1_i first.
                    a_plus_idx = triad_search(triads[a].copy(), triads, True,
                                              i, a+1, num_triads)

                    # If no match was found, then go to the next i (x,y,z).
                    if (a_plus_idx == -1):
                        break

                    # If this solution matrix element is already done, then
                    #   skip to the next i.
                    if (soln_mtx_dm[a_plus_idx][b][xyz] != ""):
                        continue

                    # For the current a+1_i, find the associated a-1_i.
                    a_minus_idx = triad_search(triads[a].copy(), triads,
                                               False, i, 0, a)

                    # For the current a+1_i, find the associated b-1_i.
                    b_minus_idx = triad_search(triads[b].copy(), triads,
                                               False, i, 0, b)

                    # At this point we have an a+1_i match and have found the
                    #   indices for a-1_i and b-1_i (or the indices are -1 if
                    #   no match was found).

                    # Now, we assemble the [xyz][a_plus_idx][b] element.
                    if (a+1 == 1 and b+1 == 1):
                        temp_string = f"PA({vec_tag_comma}{i+1})*"
                        temp_string += f"preFactorDM({vec_tag_comma}{xyz+1})"
                    else:
                        temp_string = f"PA({vec_tag_comma}{i+1})*"
                        temp_string += f"pc({vec_tag_comma}{a+1},{b+1},{xyz+1})"

                    # Add the first set of terms from lower cart. ang. mom.
                    #   This process is exactly the same as the overlap and
                    #   so we just reuse this recursive code.
                    temp_string += add_recursive_minus_terms(xyz, 7, a, b, i,
                            a_minus_idx, b_minus_idx, triads, "DM", vectorize)

                    # Add the second term (which is of equal ang. mom.)
                    #   and that explicity contains an overlap integral.
                    if (xyz == i):
                        temp_string += add_dipole_overlap_terms(a, b, i, xyz,
                                triads, True, vectorize)

                    soln_mtx_dm[a_plus_idx][b][xyz] = temp_string

                # Now do the same for the b triad. Start by adding a 1 to
                #   the 'b' triad x, y, z terms and searching for matches
                #   in the other higher triads.
                for i in range(3):

                    # Search for the index numbers in the triad list
                    #   associated with b+1_i, b-1_i, and a-1_i.

                    # Search for the b+1_i first.
                    b_plus_idx = triad_search(triads[b].copy(), triads,
                                              True, i, b+1, num_triads)

                    # If no match was found, then go to the next i (x,y,z).
                    if (b_plus_idx == -1):
                        break

                    # If this solution matrix element is already done, then
                    #   skip to the next i.
                    if (soln_mtx_dm[a][b_plus_idx][xyz] != ""):
                        continue

                    # For the current b+1_i, find the associated b-1_i.
                    b_minus_idx = triad_search(triads[b].copy(), triads,
                                               False, i, 0, b)

                    # For the current b+1_i, find the associated a-1_i.
                    a_minus_idx = triad_search(triads[a].copy(), triads,
                                               False, i, 0, a)

                    # At this point we have a b+1_i match and have found
                    #   the indices for b-1_i and a-1_i (or the indices are
                    #   -1 if no match was found.

                    # Now we assemble the [xyz][a][b_plus_idx] element.
                    if (a+1 == 1 and b+1 == 1):
                        temp_string = f"PB({vec_tag_comma}{i+1})*"
                        temp_string += f"preFactorDM({vec_tag_comma}{xyz+1})"
                    else:
                        temp_string = f"PB({vec_tag_comma}{i+1})*"
                        temp_string += \
                            f"pc({vec_tag_comma}{a+1},{b+1},{xyz+1})"

                    # Add the first set of terms from lower cart. ang. mom.
                    #   this process is exactly the same as the overlap and
                    #   so we just reuse this recursive code.
                    temp_string += add_recursive_minus_terms(xyz, 7, a, b, i,
                            a_minus_idx, b_minus_idx, triads, "DM", vectorize)

                    # Add the second term (which is of equal ang. mom.)
                    #   and that explicity contains an overlap integral.
                    if (xyz == i):
                        temp_string += add_dipole_overlap_terms(a, b, i, xyz,
                                triads, False, vectorize)

                    soln_mtx_dm[a][b_plus_idx][xyz] = temp_string

    return soln_mtx_dm


def print_production_overlap_vec(conversion, triads, matrix, lam_sh_list,
                                 lam_pc_list, f):

    head = """
   subroutine overlap2CIntgVec(A,B,numAlphaPairs,orderedAlphaPairs,sh,pc,&
         & numSegs,segIndices)

   use O_Kinds
   use O_Constants, only: pi, lAngMomCount

   implicit none

   ! sh(:,16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(:,20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (3), intent (in) :: A, B
   integer, intent (in) :: numAlphaPairs
   real (kind=double), dimension(2,numAlphaPairs), intent (in) :: &
         & orderedAlphaPairs
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(conversion)},{len(conversion)}), intent (out) :: sh" + """
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)}), intent (inout) :: pc" + """
   integer, dimension(lAngMomCount,lAngMomCount) :: numSegs
   integer, dimension(2,lAngMomCount,lAngMomCount,lAngMomCount), &
         & intent (in) :: segIndices

   ! Define local variables
   integer :: i, j
   real (kind=double) :: sum_d_sqrd
   real (kind=double), dimension(numAlphaPairs,3) :: P, PA, PB
   real (kind=double), dimension(numAlphaPairs) :: inv_zeta, inv_2zeta
   real (kind=double), dimension(numAlphaPairs) :: xi, preFactorOL

   ! Initialize local variables.
   inv_zeta(:) = 1.0d0 / (orderedAlphaPairs(1,:) + orderedAlphaPairs(2,:))
   inv_2zeta(:) = 0.5d0 * inv_zeta
   xi(:) = orderedAlphaPairs(1,:) * orderedAlphaPairs(2,:) * inv_zeta(:)
   do i = 1, 3
      P(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i)) &
            & * inv_zeta(:)
      PA(:,i) = P(:,i) - A(i)
      PB(:,i) = P(:,i) - B(i)
   enddo
   sum_d_sqrd = sum((A(:)-B(:))**2)
   preFactorOL(:) = ((pi*inv_zeta(:))**1.5d0)*exp(-xi(:)*sum_d_sqrd)

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix], [[1, 1, ""]],
                               lam_sh_list, lam_pc_list, True)

    foot = """
   end subroutine overlap2CIntgVec
    """
    f.write(foot)


def print_production_kinetic_vec(conversion, triads, matrix_ke, matrix_ol,
                                 lam_sh_list, lam_pc_list, f):

    head = """
   subroutine kinetic2CIntgVec(A,B,numAlphaPairs,orderedAlphaPairs,sh,pc,&
         & numSegs,segIndices)

   use O_Kinds
   use O_Constants, only: pi, lAngMomCount

   implicit none

   ! sh(:,16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(:,20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (3), intent (in) :: A, B
   integer, intent (in) :: numAlphaPairs
   real (kind=double), dimension(2,numAlphaPairs), intent (in) :: &
         & orderedAlphaPairs
   real (kind=double), dimension (numAlphaPairs,""" \
           +f"{len(conversion)},{len(conversion)}), intent (out) :: sh" + """
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)}), intent (inout) :: pc" + """
   integer, dimension(lAngMomCount,lAngMomCount) :: numSegs
   integer, dimension(2,lAngMomCount,lAngMomCount,lAngMomCount), &
         & intent (in) :: segIndices

   ! Define local variables
   integer :: i, j
   real (kind=double) :: sum_d_sqrd
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (numAlphaPairs,3) :: P, PA, PB
   real (kind=double), dimension (numAlphaPairs) :: inv_zeta, inv_2zeta, xi
   real (kind=double), dimension (numAlphaPairs) :: inv_2zeta_a, inv_2zeta_b
   real (kind=double), dimension (numAlphaPairs) :: preFactorOL, preFactorKE

   ! Initialize local variables.
   inv_zeta(:) = 1.0d0 / (orderedAlphaPairs(1,:) + orderedAlphaPairs(2,:))
   inv_2zeta(:) = 0.5d0 * inv_zeta(:)
   inv_2zeta_a(:) = 0.5d0 * orderedAlphaPairs(1,:)
   inv_2zeta_b(:) = 0.5d0 * orderedAlphaPairs(2,:)
   xi(:) = orderedAlphaPairs(1,:) * orderedAlphaPairs(2,:) * inv_zeta(:)
   do i = 1, 3
      P(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i)) &
            & * inv_zeta(:)
      PA(:,i) = P(:,i) - A(i)
      PB(:,i) = P(:,i) - B(i)
   enddo
   sum_d_sqrd = sum((A(:)-B(:))**2)
   preFactorOL(:) = ((pi*inv_zeta(:))**1.5d0)*exp(-xi(:)*sum_d_sqrd)
   preFactorKE(:) = xi(:)*(3.0d0 - 2.0d0*xi(:)*sum_d_sqrd)*preFactorOL(:)

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix_ol, matrix_ke],
                               [[1, 1, "_ol"], [2, 1, ""]], lam_sh_list,
                               lam_pc_list, True)

    foot = """
   end subroutine kinetic2CIntgVec
    """
    f.write(foot)


def print_production_electron_vec(conversion, triads, matrix, lam_sh_list,
                                  lam_pc_list, f):

    head = """
   subroutine electron3CIntgVec(a3,A,B,C,numAlphaPairs,orderedAlphaPairs,sh,&
         & pc,numSegs,segIndices)

   use O_Kinds
   use O_Constants, only: pi, lAngMomCount

   implicit none

   ! sh(:,16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(:,20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   integer, intent (in) :: numAlphaPairs
   real (kind=double), dimension(2,numAlphaPairs), intent (in) :: &
         & orderedAlphaPairs
   real (kind=double), dimension (numAlphaPairs,""" \
           +f"{len(conversion)},{len(conversion)}), intent (out) :: sh" + """
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)}), intent (inout) :: pc" + """
   integer, dimension(lAngMomCount,lAngMomCount) :: numSegs
   integer, dimension(2,lAngMomCount,lAngMomCount,lAngMomCount), &
         & intent (in) :: segIndices

   ! Define local variables
   integer :: i, j
   real (kind=double) :: sum_d_sqrd
   real (kind=double), dimension(numAlphaPairs) :: sum_PC_3C_sqrd
   real (kind=double), dimension(numAlphaPairs,3) :: G, GA, GB, P, PC_3C
   real (kind=double), dimension(numAlphaPairs) :: zeta, inv_zeta, inv_2zeta3C
   real (kind=double), dimension(numAlphaPairs) :: inv_zeta3C, xi, preFactorOL

   ! Initialize local variables.
   zeta(:) = orderedAlphaPairs(1,:) + orderedAlphaPairs(2,:)
   inv_zeta(:) = 1.0d0 / zeta(:)
   inv_zeta3C(:) = 1.0d0 / (zeta(:) + a3)
   inv_2zeta3C(:) = 0.5d0 * inv_zeta3C
   xi(:) = orderedAlphaPairs(1,:) * orderedAlphaPairs(2,:) * inv_zeta(:)
   do i = 1, 3
      P(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i)) &
            & * inv_zeta(:)
      G(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i) &
            & + a3*C(i)) * inv_zeta3C(:)
      GA(:,i) = G(:,i) - A(i)
      GB(:,i) = G(:,i) - B(i)
      PC_3C(:,i) = P(:,i) - C(i)
   enddo
   sum_d_sqrd = sum((A(:)-B(:))**2)
   sum_PC_3C_sqrd(:) = 0.0d0
   do i = 1, 3
      sum_PC_3C_sqrd(:) = sum_PC_3C_sqrd(:) + PC_3C(:,i)**2
   enddo
   prefactorOL(:) = ((pi*inv_zeta3C(:))**1.5) &
         & * exp(-xi(:)*sum_d_sqrd &
         &       - (zeta(:)*a3*inv_zeta3C(:))*sum_PC_3C_sqrd(:))

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix], [[3, 1, ""]],
                               lam_sh_list, lam_pc_list, True)

    foot = """
   end subroutine electron3CIntgVec
    """
    f.write(foot)


def print_production_nuclear_vec(conversion, triads, matrix, max_lam,
                                 lam_sh_list, lam_pc_list, f):

    head = """
   subroutine nuclear3CIntgVec(a3,A,B,C,numAlphaPairs,orderedAlphaPairs,sh,&
         & pc,numSegs,segIndices)

   use O_Kinds
   use O_Constants, only: pi, lAngMomCount

   implicit none

   ! sh(:,16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(:,20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   integer, intent (in) :: numAlphaPairs
   real (kind=double), dimension(2,numAlphaPairs), intent (in) :: &
         & orderedAlphaPairs
   real (kind=double), dimension (numAlphaPairs,""" \
           +f"{len(conversion)},{len(conversion)}), intent (out) :: sh" + """
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)},{2*max_lam+1})," \
           + f"intent (inout) :: pc" + """
   integer, dimension(lAngMomCount,lAngMomCount) :: numSegs
   integer, dimension(2,lAngMomCount,lAngMomCount,lAngMomCount), &
         & intent (in) :: segIndices

   ! Define local variables
   integer :: i, j, m
   real (kind=double) :: sum_d1_sqrd, sum_d2_sqrd, sum_d3_sqrd
   real (kind=double), dimension (3) :: d1, d2, d3
   real (kind=double), dimension (numAlphaPairs,3) :: G, GA, GB, GC, P, PC_3C
   real (kind=double), dimension (numAlphaPairs) :: zeta, zeta3C, sum_GC_sqrd
   real (kind=double), dimension (numAlphaPairs) :: inv_zeta3C, inv_2zeta3C, U
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{2*max_lam+1}" + """) :: preFactorN
   real (kind=double), dimension (numAlphaPairs,7) :: F

   ! Initialize local variables.
   zeta(:) = orderedAlphaPairs(1,:) + orderedAlphaPairs(2,:)
   zeta3C(:) = zeta(:) + a3
   inv_zeta3C(:) = 1.0d0 / zeta3C(:)
   inv_2zeta3C(:) = 0.5d0 * inv_zeta3C(:)
   do i = 1, 3
      P(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i)) &
            & / zeta(:)
      G(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i) &
            & + a3*C(i)) * inv_zeta3C(:)
      GA(:,i) = G(:,i) - A(i)
      GB(:,i) = G(:,i) - B(i)
      GC(:,i) = G(:,i) - C(i)
      PC_3C(:,i) = P(:,i) - C(i)
   enddo
   sum_d1_sqrd = sum((A(:)-B(:))**2)
   sum_d2_sqrd = sum((A(:)-C(:))**2)
   sum_d3_sqrd = sum((B(:)-C(:))**2)
   sum_GC_sqrd(:) = 0.0d0
   do i = 1, 3
      sum_GC_sqrd(:) = sum_GC_sqrd(:) + GC(:,i)**2
   enddo
   U(:) = zeta3C(:) * sum_GC_sqrd
   call boys_vec(numAlphaPairs,U,F)

   do m = 1, """ + f"{2*max_lam + 1}" + """
      preFactorN(:,m) = F(:,m) * 2.0d0*(pi*inv_zeta3C(:)) * exp(- &
            & (orderedAlphaPairs(1,:)*orderedAlphaPairs(2,:)*sum_d1_sqrd &
            & + orderedAlphaPairs(1,:)*a3*sum_d2_sqrd &
            & + orderedAlphaPairs(2,:)*a3*sum_d3_sqrd) * inv_zeta3C(:))
   enddo
"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix], [[4, 1, ""]],
                               lam_sh_list, lam_pc_list, True)

    foot = """
   end subroutine nuclear3CIntgVec
    """
    f.write(foot)


def print_production_momentum_vec(conversion, triads, matrix_mm, matrix_ol, 
                                  lam_sh_list, lam_pc_list, f):

    head = """
   subroutine momentum2CIntgVec(A,B,numAlphaPairs,orderedAlphaPairs,sh,pc,&
         & numSegs,segIndices)

   use O_Kinds
   use O_Constants, only: pi, lAngMomCount

   implicit none

   ! sh(:,16,16,3): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(:,20,20,3): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (3), intent (in) :: A, B
   integer, intent (in) :: numAlphaPairs
   real (kind=double), dimension(2,numAlphaPairs), intent (in) :: &
         & orderedAlphaPairs
   real (kind=double), dimension (numAlphaPairs,""" \
           +f"{len(conversion)},{len(conversion)},3), intent (out) :: sh" + """
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)},3), intent (inout) :: pc" + """
   integer, dimension(lAngMomCount,lAngMomCount) :: numSegs
   integer, dimension(2,lAngMomCount,lAngMomCount,lAngMomCount), &
         & intent (in) :: segIndices

   ! Define local variables
   integer :: i, j
   real (kind=double) :: sum_d_sqrd
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (numAlphaPairs,3) :: P, PA, PB, preFactorMM
   real (kind=double), dimension (numAlphaPairs) :: inv_zeta, inv_2zeta
   real (kind=double), dimension (numAlphaPairs) :: xi, preFactorOL
   real (kind=double), dimension (numAlphaPairs) :: zeta_a_zeta, zeta_b_zeta

   ! Initialize local variables.
   inv_zeta(:) = 1.0d0 / (orderedAlphaPairs(1,:) + orderedAlphaPairs(2,:))
   zeta_a_zeta(:) = orderedAlphaPairs(1,:) * inv_zeta(:)
   zeta_b_zeta(:) = orderedAlphaPairs(2,:) * inv_zeta(:)
   inv_2zeta(:) = 0.5d0 * inv_zeta(:)
   xi(:) = orderedAlphaPairs(1,:) * orderedAlphaPairs(2,:) * inv_zeta(:)
   do i = 1, 3
      P(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i)) &
            & * inv_zeta(:)
      PA(:,i) = P(:,i) - A(i)
      PB(:,i) = P(:,i) - B(i)
   enddo
   sum_d_sqrd = sum((A(:)-B(:))**2)
   preFactorOL(:) = ((pi*inv_zeta(:))**1.5d0)*exp(-xi(:)*sum_d_sqrd)

   ! Slightly wonky initialization of MM prefactor to avoid repeated
   !   reading of orderedAlphaPairs(1,:) and repeated multiplications.
   preFactorMM(:,1) = preFactorOL(:)*2.0d0*orderedAlphaPairs(1,:)
   preFactorMM(:,2) = preFactorMM(:,1)*PA(:,2)
   preFactorMM(:,3) = preFactorMM(:,1)*PA(:,3)
   preFactorMM(:,1) = preFactorMM(:,1)*PA(:,1)

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix_ol, matrix_mm],
                               [[1, 1, "_ol"], [5, 3, ""]], lam_sh_list,
                               lam_pc_list, True)

    foot = """
   end subroutine momentum2CIntgVec
    """
    f.write(foot)


def print_production_massvel_vec(conversion, triads, matrix_mv, matrix_ol,
                                 lam_sh_list, lam_pc_list, f):

    head = """
   subroutine massvel2CIntgVec(A,B,numAlphaPairs,orderedAlphaPairs,sh,pc,&
         & numSegs,segIndices)

   use O_Kinds
   use O_Constants, only: pi, lAngMomCount

   implicit none

   ! sh(:,16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(:,20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (3), intent (in) :: A, B
   integer, intent (in) :: numAlphaPairs
   real (kind=double), dimension(2,numAlphaPairs), intent (in) :: &
         & orderedAlphaPairs
   real (kind=double), dimension (numAlphaPairs,""" \
           +f"{len(conversion)},{len(conversion)}), intent (out) :: sh" + """
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)}), intent (inout) :: pc" + """
   integer, dimension(lAngMomCount,lAngMomCount) :: numSegs
   integer, dimension(2,lAngMomCount,lAngMomCount,lAngMomCount), &
         & intent (in) :: segIndices

   ! Define local variables
   integer :: i, j
   real (kind=double) :: sum_d_sqrd
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (numAlphaPairs,3) :: P, PA, PB
   real (kind=double), dimension (numAlphaPairs) :: inv_zeta, inv_2zeta, xi
   real (kind=double), dimension (numAlphaPairs) :: inv_2zeta_a, inv_2zeta_b
   real (kind=double), dimension (numAlphaPairs) :: preFactorOL, preFactorMV

   ! Initialize local variables.
   inv_zeta(:) = 1.0d0 / (orderedAlphaPairs(1,:) + orderedAlphaPairs(2,:))
   inv_2zeta(:) = 0.5d0 * inv_zeta(:)
   inv_2zeta_a(:) = 0.5d0 * orderedAlphaPairs(1,:)
   inv_2zeta_b(:) = 0.5d0 * orderedAlphaPairs(2,:)
   xi(:) = orderedAlphaPairs(1,:) * orderedAlphaPairs(2,:) * inv_zeta(:)
   do i = 1, 3
      P(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i)) &
            & * inv_zeta(:)
      PA(:,i) = P(:,i) - A(i)
      PB(:,i) = P(:,i) - B(i)
   enddo
   sum_d_sqrd = sum((A(:)-B(:))**2)
   preFactorOL(:) = ((pi*inv_zeta(:))**1.5d0)*exp(-xi(:)*sum_d_sqrd)
   preFactorMV(:) = xi(:)*(3.0d0 - 2.0d0*xi(:)*sum_d_sqrd)*preFactorOL(:)

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f,
                               [matrix_ol, matrix_ke, matrix_mv],
                               [[1, 1, "_ol"], [2, 1, "_ke"], [6, 1, ""]],
                               lam_sh_list, lam_pc_list, True)

    foot = """
   end subroutine massvel2CIntgVec
    """
    f.write(foot)


def print_production_dipole_vec(conversion, triads, matrix_dm, matrix_ol, 
                                lam_sh_list, lam_pc_list, f):

    head = """
   subroutine dipole3CIntgVec(A,B,C,numAlphaPairs,orderedAlphaPairs,sh,pc,&
         & numSegs,segIndices)

   use O_Kinds
   use O_Constants, only: pi, lAngMomCount

   implicit none

   ! sh(:,16,16,3): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(:,20,20,3): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (3), intent (in) :: A, B, C
   integer, intent (in) :: numAlphaPairs
   real (kind=double), dimension(2,numAlphaPairs), intent (in) :: &
         & orderedAlphaPairs
   real (kind=double), dimension (numAlphaPairs,""" \
           +f"{len(conversion)},{len(conversion)},3), intent (out) :: sh" + """
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)},3), intent (inout) :: pc" + """
   integer, dimension(lAngMomCount,lAngMomCount) :: numSegs
   integer, dimension(2,lAngMomCount,lAngMomCount,lAngMomCount), &
         & intent (in) :: segIndices

   ! Define local variables
   integer :: i, j
   real (kind=double) :: sum_d_sqrd
   real (kind=double), dimension (numAlphaPairs,""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (numAlphaPairs,3) :: P, PA, PB, PC_3C
   real (kind=double), dimension (numAlphaPairs,3) :: preFactorDM
   real (kind=double), dimension (numAlphaPairs) :: inv_zeta, inv_2zeta
   real (kind=double), dimension (numAlphaPairs) :: xi, preFactorOL, mu
   real (kind=double), dimension (numAlphaPairs) :: zeta_a_zeta, zeta_b_zeta

   ! Initialize local variables.
   inv_zeta(:) = 1.0d0 / (orderedAlphaPairs(1,:) + orderedAlphaPairs(2,:))
   zeta_a_zeta(:) = orderedAlphaPairs(1,:) * inv_zeta(:)
   zeta_b_zeta(:) = orderedAlphaPairs(2,:) * inv_zeta(:)
   inv_2zeta(:) = 0.5d0 * inv_zeta(:)
   xi(:) = orderedAlphaPairs(1,:) * orderedAlphaPairs(2,:) * inv_zeta(:)
   do i = 1, 3
      P(:,i) = (orderedAlphaPairs(1,:)*A(i) + orderedAlphaPairs(2,:)*B(i)) &
            & * inv_zeta(:)
      PA(:,i) = P(:,i) - A(i)
      PB(:,i) = P(:,i) - B(i)
      PC_3C(:,i) = P(:,i) - C(i)
   enddo
   sum_d_sqrd = sum((A(:)-B(:))**2)
   preFactorOL(:) = ((pi*inv_zeta(:))**1.5d0)*exp(-xi(:)*sum_d_sqrd)

   ! Initialization of DM prefactor.
   preFactorDM(:,1) = preFactorOL(:)*PC_3C(:,1)
   preFactorDM(:,2) = preFactorOL(:)*PC_3C(:,2)
   preFactorDM(:,3) = preFactorOL(:)*PC_3C(:,3)

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix_ol, matrix_dm],
                               [[1, 1, "_ol"], [7, 3, ""]], lam_sh_list,
                               lam_pc_list, True)

    foot = """
   end subroutine dipole3CIntgVec
    """
    f.write(foot)


def print_production_overlap(conversion, triads, matrix, lam_sh_list,
                             lam_pc_list, f):

    head = """
   subroutine overlap2CIntg(a1,a2,A,B,l1l2switch,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   integer, intent (in) :: l1l2switch
   real (kind=double), dimension (""" \
           +f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc" + """
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d**2))

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix], [[1, 1, ""]],
                               lam_sh_list, lam_pc_list, False)

    foot = """
   end subroutine overlap2CIntg
    """
    f.write(foot)


def print_production_kinetic(conversion, triads, matrix_ke, matrix_ol,
                             lam_sh_list, lam_pc_list,  f):

    head = """
   subroutine kinetic2CIntg(a1,a2,A,B,l1l2switch,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   integer, intent (in) :: l1l2switch
   real (kind=double), dimension (""" \
           +f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc" + """
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL, preFactorKE
   real (kind=double) :: inv_2zeta_a, inv_2zeta_b

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   inv_2zeta_a = 1.0d0 / (2.0d0 * a1)
   inv_2zeta_b = 1.0d0 / (2.0d0 * a2)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d*d))
   preFactorKE = xi*(3 - 2*xi*sum(d*d))*preFactorOL

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix_ol, matrix_ke],
                               [[1, 1, "_ol"], [2, 1, ""]], lam_sh_list,
                               lam_pc_list, False)

    foot = """
   end subroutine kinetic2CIntg
    """
    f.write(foot)


def print_production_electron(conversion, triads, matrix, lam_sh_list,
                              lam_pc_list, f):

    head = """
   subroutine electron3CIntg(a1,a2,a3,A,B,C,l1l2switch,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   integer, intent (in) :: l1l2switch
   real (kind=double), dimension (""" \
           +f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc" + """
   real (kind=double), dimension (3) :: G, GA, GB, P, PC_3C, d
   real (kind=double) :: zeta, zeta3C, inv_2zeta3C, xi
   real (kind=double) :: preFactorOL

   ! Initialize local variables.
   zeta = a1 + a2
   zeta3C = zeta + a3
   inv_2zeta3C = 1.0d0 / (2.0d0 * zeta3C)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   G = (a1*A + a2*B + a3*C) / zeta3C
   GA = G - A
   GB = G - B
   PC_3C = P - C
   d = A - B
   prefactorOL = ((pi/zeta3C)**1.5) &
         & * exp(-xi*sum(d**2)-(zeta*a3/zeta3C)*sum(PC_3C**2))

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix], [[3, 1, ""]],
                               lam_sh_list, lam_pc_list, False)

    foot = """
   end subroutine electron3CIntg
    """
    f.write(foot)


def print_production_nuclear(conversion, triads, matrix, max_lam, lam_sh_list,
                             lam_pc_list, f):

    head = """
   subroutine nuclear3CIntg(a1,a2,a3,A,B,C,l1l2switch,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   integer, intent (in) :: l1l2switch
   real (kind=double), dimension (""" \
    + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables
   integer :: m
   real (kind=double), dimension (""" \
    + f"{len(triads)},{len(triads)},{2*max_lam+1}) :: pc" + """
   real (kind=double), dimension (3) :: G, GA, GB, GC, P, PC_3C, d1, d2, d3
   real (kind=double) :: zeta, zeta3C, inv_2zeta3C, U
   real (kind=double), dimension (""" + f"{2*max_lam+1}" + """) :: preFactorN
   real (kind=double), dimension (7) :: F

   ! Initialize local variables.
   zeta = a1 + a2 
   zeta3C = zeta + a3
   inv_2zeta3C = 1.0d0 / (2.0d0 * zeta3C)
   P = (a1*A + a2*B) / zeta
   G = (a1*A + a2*B + a3*C) / zeta3C
   GA = G - A
   GB = G - B
   GC = G - C
   PC_3C = P - C
   d1 = A - B
   d2 = A - C
   d3 = B - C
   U = zeta3C * sum(GC**2)
   call boys(U,F)

   do m = 1, """ + f"{2*max_lam + 1}" + """
      preFactorN(m) = F(m) * 2.0d0*(pi/zeta3C) * &
            & exp(-sum(a1*a2*(d1**2) + a1*a3*(d2**2) + a2*a3*(d3**2))/zeta3C)
   enddo

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix], [[4, 1, ""]],
                               lam_sh_list, lam_pc_list, False)

    foot = """
   end subroutine nuclear3CIntg
    """
    f.write(foot)


def print_production_momentum(conversion, triads, matrix_mm, matrix_ol,
                              lam_sh_list, lam_pc_list, f):

    head = """
   subroutine momentum2CIntg(a1,a2,A,B,l1l2switch,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   integer, intent (in) :: l1l2switch
   real (kind=double), dimension (""" \
           +f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """

   ! Define local variables
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)},3) :: pc" + """
   real (kind=double), dimension (3) :: P, PA, PB, d, preFactorMM
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL
   real (kind=double) :: zeta_a_zeta, zeta_b_zeta

   ! Initialize local variables.
   zeta = a1 + a2
   zeta_a_zeta = a1/zeta
   zeta_b_zeta = a2/zeta
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d*d))
   preFactorMM(1) = preFactorOL*2.0d0*a1*PA(1) ! a or b arbitrary
   preFactorMM(2) = preFactorOL*2.0d0*a1*PA(2) ! a or b arbitrary
   preFactorMM(3) = preFactorOL*2.0d0*a1*PA(3) ! a or b arbitrary

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix_ol, matrix_mm],
                               [[1, 1, "_ol"], [5, 3, ""]],
                               lam_sh_list, lam_pc_list, False)

    foot = """
   end subroutine momentum2CIntg
    """
    f.write(foot)


def print_production_massvel(conversion, triads, matrix_mv, matrix_ke,
                             matrix_ol, lam_sh_list, lam_pc_list, f):

    head = """
   subroutine massvel2CIntg(a1,a2,A,B,l1l2switch,sh)

   use O_Kinds
   use O_Constants, only: pi, fineStructure

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   integer, intent (in) :: l1l2switch
   real (kind=double), dimension (""" \
           +f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ke" + """
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc" + """
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi
   real (kind=double) :: preFactorOL, preFactorKE
   real (kind=double) :: preFactor02, preFactor04
   real (kind=double) :: preFactor22, preFactorMV
   real (kind=double) :: inv_2zeta_a, inv_2zeta_b, inv_8m3c2

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   inv_2zeta_a = 1.0d0 / (2.0d0 * a1)
   inv_2zeta_b = 1.0d0 / (2.0d0 * a2)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d*d))
   preFactorKE = xi*(3 - 2*xi*sum(d*d))*preFactorOL
   preFactor02 = sum(PB(:)**2) + 3.0d0/(2.0d0*zeta)
   preFactor04 = sum(PB(:)**4) + sum(PB(:)**2)*3.0d0/zeta &
         & + 9.0d0/(4.0d0*zeta**2)
   preFactor22 = (PB(1)*PB(2))**2 + (PB(1)*PB(3))**2 + (PB(2)*PB(3))**2 &
         & + sum(PB(:)**2)/zeta + 3.0d0/(4.0d0*zeta**2)
   preFactorMV = (fineStructure * 0.001d0)**2 / 8.0d0 &
         & * (16*a2**4*preFactor04 - 80*a2**3*preFactor02 + 60*a2**2 &
         & + 32*a2**4 * preFactor22) * preFactorOL
   inv_8m3c2 = (fineStructure * 0.001d0)**2 / 8.0d0

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f,
                               [matrix_ol, matrix_ke, matrix_mv],
                               [[1, 1, "_ol"], [2, 1, "_ke"], [6, 1, ""]],
                               lam_sh_list, lam_pc_list, False)

    foot = """
   end subroutine massvel2CIntg
    """
    f.write(foot)


def print_production_dipole(conversion, triads, matrix_dm, matrix_ol,
                              lam_sh_list, lam_pc_list, f):

    head = """
   subroutine dipole3CIntg(a1,a2,A,B,C,l1l2switch,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B, C
   integer, intent (in) :: l1l2switch
   real (kind=double), dimension (""" \
           +f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """

   ! Define local variables
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)},3) :: pc" + """
   
   real (kind=double), dimension (3) :: P, PA, PB, PC_3C, d, preFactorDM
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL, mu
   real (kind=double), dimension (3,3) :: dip_mom


   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   d = A - B
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   PC_3C = P - C
   
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d**2))
   
   ! These are the (s|M(mu)|s) integrals
   preFactorDM(1) = preFactorOL*PC_3C(1)
   preFactorDM(2) = preFactorOL*PC_3C(2)
   preFactorDM(3) = preFactorOL*PC_3C(3)


"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    lib.print_production_pc_sh(conversion, triads, f, [matrix_ol, matrix_dm],
                               [[1, 1, "_ol"], [7, 3, ""]],
                               lam_sh_list, lam_pc_list, False)

    foot = """
   end subroutine dipole3CIntg
    """
    f.write(foot)



def print_test_overlap_ana(conversion, triads, matrix, f):

    # Print the subroutine header for the analytical portion.
    head = """
   subroutine overlap2CIntgAna(a1,a2,A,B,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables.
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d**2))

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    num_triads = len(triads)
    lib.print_pc(conversion, f, [matrix], [[1, 1, ""]], num_triads,
                 num_triads)

    # Print the subroutine foot.
    foot = """
   end subroutine overlap2CIntgAna
"""
    f.write(foot)


def print_test_kinetic_ana(conversion, triads, matrix_ke, matrix_ol, f):

    # Print the subroutine header for the analytical portion.
    head = """
   subroutine kinetic2CIntgAna(a1,a2,A,B,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables.
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL, preFactorKE
   real (kind=double) :: altPreFactorKE
   real (kind=double) :: inv_2zeta_a, inv_2zeta_b

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   inv_2zeta_a = 1.0d0 / (2.0d0 * a1)
   inv_2zeta_b = 1.0d0 / (2.0d0 * a2)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d*d))
   preFactorKE = xi*(3 - 2*xi*sum(d*d))*preFactorOL
   altPreFactorKE = a2*(3.0d0 - 2.0d0*a2*(sum(PB*PB) + 3.0d0/(2.0d0*zeta))) &
         & * preFactorOL

   write (6,*) "preFactorKE = ", preFactorKE
   write (6,*) "altPreFactorKE = ", altPreFactorKE

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    num_triads = len(triads)
    lib.print_pc(conversion, f, [matrix_ol, matrix_ke],
                 [[1, 1, "_ol"], [2, 1, ""]], num_triads, num_triads)

    # Print the subroutine foot.
    foot = """
   end subroutine kinetic2CIntgAna
"""
    f.write(foot)


def print_test_electron_ana(conversion, triads, matrix, f):

    # Print the subroutine header for the analytical portion.
    head = """
   subroutine electron3CIntgAna(a1,a2,a3,A,B,C,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables.
   real (kind=double), dimension (3) :: G, GA, GB, P, PC_3C, d
   real (kind=double) :: zeta, zeta3C, inv_2zeta3C, xi
   real (kind=double) :: preFactorOL

   ! Initialize local variables.
   zeta = a1 + a2
   zeta3C = zeta + a3
   inv_2zeta3C = 1.0d0 / (2.0d0 * zeta3C)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   G = (a1*A + a2*B + a3*C) / zeta3C
   GA = G - A
   GB = G - B
   PC_3C = P - C
   d = A - B
   prefactorOL = ((pi/zeta3C)**1.5) &
         & * exp(-xi*sum(d**2)-(zeta*a3/zeta3C)*sum(PC_3C**2))

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    num_triads = len(triads)
    lib.print_pc(conversion, f, [matrix], [[3, 1, ""]], num_triads,
                 num_triads)

    # Print the subroutine foot.
    foot = """
   end subroutine electron3CIntgAna
"""
    f.write(foot)


def print_test_nuclear_ana(conversion, triads, matrix, max_lam, f):

    # Note the number of triads.
    num_triads = len(triads)

    # Print the subroutine header for the analytical portion.
    head = """
   subroutine nuclear3CIntgAna(a1,a2,a3,A,B,C,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2, a3
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (""" \
    + f"{num_triads},{num_triads},{2*max_lam+1}), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
    + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables
   integer :: m
   real (kind=double), dimension (3) :: G, GA, GB, GC, P, PC_3C, d1, d2, d3
   real (kind=double) :: zeta, zeta3C, inv_2zeta3C, U
   real (kind=double), dimension (""" + f"{2*max_lam+1}" + """) :: preFactorN
   real (kind=double), dimension (7) :: F

   ! Initialize local variables.
   zeta = a1 + a2
   zeta3C = zeta + a3
   inv_2zeta3C = 1.0d0 / (2.0d0 * zeta3C)
   P = (a1*A + a2*B) / zeta
   G = (a1*A + a2*B + a3*C) / zeta3C
   GA = G - A
   GB = G - B
   GC = G - C
   PC_3C = P - C
   d1 = A - B
   d2 = A - C
   d3 = B - C
   U = zeta3C * sum(GC**2)
   call boys(U,F)

   do m = 1, """ + f"{2*max_lam + 1}" + """
      preFactorN(m) = F(m) * 2.0d0*(pi/zeta3C) * &
            & exp(-sum(a1*a2*(d1**2) + a1*a3*(d2**2) + a2*a3*(d3**2))/zeta3C)
   enddo

"""
    f.write(head)

    num_conversion = len(conversion)
    num_triads = len(triads)
    lib.print_nuclear_pc(f, matrix, num_triads, num_triads)
    lib.print_pc_to_sh(conversion, f, [4, 1, ""], 0, num_conversion,
                       num_conversion)

    # Print the subroutine foot and the boys function subroutine.
    foot = """
   end subroutine nuclear3CIntgAna
"""
    f.write(foot)


def print_boys_vec(f):

    # Print the boys function subroutine.
    head = """
   ! This subroutine computes F_N(T) - which is the Incomplete
   ! Gamma Function, aka the Boys Function. This comes about as a 
   ! result of the integral transforms that are required to 
   ! account for the (1/r) term that is present in the Coulomb 
   ! integral, to make the integral separable in the Cartesian 
   ! directions (x,y,z).  The formalism is Eq. (29) on page 41 of the 
   ! following document:
   !
   ! ****************************************************************
   ! Petersson T., Hellsing B., "A Detailed Derivation of Gaussian
   ! Orbital-Based Matrix Elements In Electron Structure Calculations."
   ! European Journal of Physics, Vol. 31, No. 1, Page 37.  
   ! ****************************************************************
   !
   ! The Boys Function itself is defined as:
   ! F_N(T) = integral [from 0 to 1] {(t**(2*N))*exp(-T*(t**2))}dt
   !
   ! where T = (a1 + a2 + a3)*sum(P - C)**2
   !
   ! For our purposes here, N will be an integer >=0 that is a linear
   ! combination of index values that come from the summations that
   ! compute the main body of the Coulomb integral, be it Nuclear
   ! Attraction, or the more complicated Electron Repulsion Integral.
   !
   ! Generally, the solution to this integral is expressed numerically
   ! in terms of standard error functions (erf), but this closed form
   ! expression becomes numerically unstable if T = 0 or is very small 
   ! (~10E-2 or less).  To avoid the large roundoff errors that occur
   ! as a result of this, a check has to be performed on the value of T.  
   ! If T <= 10E-2, the following formalism shown by Cook (referenced below)
   ! needs to be applied (Cook: p. 260):
   !
   ! F_N(T) = (1/2)*exp(-T)*sum[from i=0 to infinity] of
   !                            {gamma(N+1/2)/gamma(N+i+3/2)}
   !
   ! This infinite series is truncated at i = 6 for our desired accuracy,
   ! including for the case where T = 0.
   !
   ! More detailed in formation on the Boys Function and numerical 
   ! approximations of it are detailed in the following sources:
   !
   ! ****************************************************************
   ! Cook, David B., "Handbook of Computational Quantum Chemistry."
   ! Dover Publications (August 22, 2005).
   ! ISBN-10: 0486443078
   ! ISBN-13: 978-0486443072
   ! ****************************************************************
   !
   ! as well as:
   !
   ! ****************************************************************
   ! Helgaker T., Taylor P., "Gaussian Basis Sets and Molecular 
   ! Integrals", Ch. 12 from "Modern Electronic Structure Theory."
   ! Pp. 809 - 815.
   ! ****************************************************************
   subroutine boys_vec(numAlphaPairs,T,F)
   use o_kinds
   use o_constants, only: pi
   implicit none
   
   ! Define passed parameters
   integer :: numAlphaPairs
   real (kind=double), dimension(numAlphaPairs), intent (in) :: T
   real (kind=double), dimension(numAlphaPairs,7), intent (out) :: F

   ! Define local variables
   integer :: N, i, j
   integer :: lessCount, moreCount, max_count
   integer, dimension(numAlphaPairs) :: record_T
   real (kind=double) :: sqrt_pi
   real (kind=double), dimension(numAlphaPairs) :: less_T, more_T
   real (kind=double), allocatable, dimension(:,:) :: less_F, more_F
   real (kind=double), allocatable, dimension(:) :: erf_T
   real (kind=double), allocatable, dimension(:,:) :: exp_T, sqrt_T
   real (kind=double), allocatable, dimension(:,:,:) :: pow_T
   real (kind=double), allocatable, dimension(:,:) :: half_pow_T
   real (kind=double), dimension(7) :: a0, a1, a2, a3, a4, a5, a6, a7

   ! Initialize local variables that are used regardless of the method used
   !   to evaluate the boys function (expansion or analytic).

   ! As stated earler, if T is small (<= 10E-2), the closed form solution 
   !   to the Boys Function can't be used because it becomes numerically
   !   unstable.  Instead, the following expansion is used to approximate it
   !   for small values (Cook: p. 260).  The loop defines the Nth degree Boys
   !   Function.
   lessCount = 0
   moreCount = 0
   do i = 1, numAlphaPairs
      if (T(i) <= 10E-2) then
         lessCount = lessCount + 1
         less_T(lessCount) = T(i)
         record_T(i) = 1
      else
         moreCount = moreCount + 1
         more_T(moreCount) = T(i)
         record_T(i) = 2
      endif
   enddo

   ! Allocate space for the temporary variables used for both less and more.
   max_count = max(moreCount,lessCount)
   allocate(exp_T(max_count,2))
   allocate(sqrt_T(max_count,2))
   allocate(pow_T(max_count,1,2))
   if (lessCount > 0) then
      allocate(less_F(lessCount,7))
   endif
   if (moreCount > 0) then
      allocate(erf_T(moreCount))
      allocate(more_F(moreCount,7))
      allocate(half_pow_T(moreCount,7))
   endif

   ! Precompute relevant quantities for both less (index 1, the expansion
   !   based method) and more (index 2, the analytic based method).
   sqrt_pi = sqrt(pi)

   if (lessCount > 0) then
      exp_T(1:lessCount,1) = exp(-less_T(1:lessCount))
      sqrt_T(1:lessCount,1) = sqrt(less_T(1:lessCount))
      pow_T(1:lessCount,1,1) = less_T(1:lessCount)
      do i = 2, 6
         pow_T(1:lessCount,i,1) = pow_T(1:lessCount,i,1) * less_T(1:lessCount)
      enddo
   endif

   if (moreCount > 0) then
      exp_T(1:moreCount,2) = exp(-more_T(1:moreCount))
      sqrt_T(1:moreCount,2) = sqrt(more_T(1:moreCount))
      erf_T(1:moreCount) = erf(sqrt_T(1:moreCount,2))
      pow_T(1:moreCount,1,2) = more_T(1:moreCount)
      do i = 2, 6
         pow_T(1:moreCount,i,2) = pow_T(1:moreCount,i,2) * more_T(1:moreCount)
      enddo
   endif

   if (lessCount > 0) then
      ! The values (a0 through a7) are computed gamma(N + 0.5) values
      ! to be plugged into the below approximation of the boys Function:
      a0(1) = 1.7724538509055161     
      a0(2) = 0.88622692545275805
      a0(3) = 1.3293403881791370
      a0(4) = 3.3233509704478426
      a0(5) = 11.631728396567450
      a0(6) = 52.342777784553533
      a0(7) = 287.88527781504416
      a1(1) = 0.88622692545275805     
      a1(2) = 1.3293403881791370
      a1(3) = 3.3233509704478426
      a1(4) = 11.631728396567450
      a1(5) = 52.342777784553533
      a1(6) = 287.88527781504416
      a1(7) = 1871.2543057977896
      a2(1) = 1.3293403881791370     
      a2(2) = 3.3233509704478426
      a2(3) = 11.631728396567450
      a2(4) = 52.342777784553533
      a2(5) = 287.88527781504416
      a2(6) = 1871.2543057977896
      a2(7) = 14034.407293483402
      a3(1) = 3.3233509704478426     
      a3(2) = 11.631728396567450
      a3(3) = 52.342777784553533
      a3(4) = 287.88527781504416
      a3(5) = 1871.2543057977896
      a3(6) = 14034.407293483402
      a3(7) = 119292.46199460901
      a4(1) = 11.631728396567450     
      a4(2) = 52.342777784553533
      a4(3) = 287.88527781504416
      a4(4) = 1871.2543057977896
      a4(5) = 14034.407293483402
      a4(6) = 119292.46199460901
      a4(7) = 1133278.3889487833
      a5(1) = 52.342777784553533     
      a5(2) = 287.88527781504416
      a5(3) = 1871.2543057977896
      a5(4) = 14034.407293483402
      a5(5) = 119292.46199460901
      a5(6) = 1133278.3889487833
      a5(7) = 11899423.083962219
      a6(1) = 287.88527781504416     
      a6(2) = 1871.2543057977896
      a6(3) = 14034.407293483402
      a6(4) = 119292.46199460901
      a6(5) = 1133278.3889487833
      a6(6) = 11899423.083962219
      a6(7) = 136843365.46556622
      a7(1) = 1871.2543057977896     
      a7(2) = 14034.407293483402
      a7(3) = 119292.46199460901
      a7(4) = 1133278.3889487833
      a7(5) = 11899423.083962219
      a7(6) = 136843365.46556622
      a7(7) = 1710542068.3195722

      do N = 1, 7
         !
         ! S(N) = 0.5d0*exp_T*(gamma(N-1 + 0.5d0)/gamma(N-1 + 1.5d0)&
         !   &+ (gamma(N-1 + 0.5d0)*T)/gamma(N-1 + 2.5d0)& 
         !   &+ (gamma(N-1 + 0.5d0)*T**2)/gamma(N-1 + 3.5d0)& 
         !   &+ (gamma(N-1 + 0.5d0)*T**3)/gamma(N-1 + 4.5d0)&
         !   &+ (gamma(N-1 + 0.5d0)*T**4)/gamma(N-1 + 5.5d0)& 
         !   &+ (gamma(N-1 + 0.5d0)*T**5)/gamma(N-1 + 6.5d0)& 
         !   &+ (gamma(N-1 + 0.5d0)*T**6)/gamma(N-1 + 7.5d0))

         less_F(1:lessCount,N) = 0.5d0*exp_T(1:lessCount,1)*a0(N)*((1/a1(N)) &
               & + (pow_T(1:lessCount,1,1) / a2(N)) &
               & + ((pow_T(1:lessCount,2,1)) / a3(N)) &
               & + ((pow_T(1:lessCount,3,1)) / a4(N)) &
               & + ((pow_T(1:lessCount,4,1)) / a5(N)) &
               & + ((pow_T(1:lessCount,5,1)) / a6(N)) &
               & + ((pow_T(1:lessCount,6,1)) / a7(N)))
      end do
   endif

   if (moreCount > 0) then

      allocate(half_pow_T(moreCount,7))
      half_pow_T(:,1) = sqrt_T(:,2)
      do i = 2, 7
         half_pow_T(:,i) = half_pow_T(:,i-1) * more_T(:)
      enddo

more_F(1:moreCount,1) = 0.5d0*((sqrt_pi/(half_pow_T(1:moreCount,1))) &
      & * erf_T(1:moreCount))

more_F(1:moreCount,2) = ((sqrt_pi/(4.0d0*half_pow_T(1:moreCount,2))) &
      & * erf_T(1:moreCount) - exp_T(1:moreCount,2) &
      & * (0.5d0/more_T(1:moreCount)))

more_F(1:moreCount,3) = 6.0d0 &
      & * ((sqrt_pi/(16*half_pow_T(1:moreCount,3))) &
      & * erf_T(1:moreCount) - exp_T(1:moreCount,2) &
      & * (2.0d0/(24d0*more_T(1:moreCount)) &
      & + 1.0d0/(4.0d0*2.0d0*pow_T(1:moreCount,2,2))))

more_F(1:moreCount,4) = 60.0d0 &
      & * ((sqrt_pi/(64*half_pow_T(1:moreCount,4))) &
      & * erf_T(1:moreCount) - exp_T(1:moreCount,2) &
      & * (6.0d0/(720d0*more_T(1:moreCount)) &
      & + 2/((4)*24d0*pow_T(1:moreCount,2,2)) &
      & + 1/((16.0d0)*2d0*pow_T(1:moreCount,3,2))))

more_F(1:moreCount,5) = 840.0d0 &
      & * ((sqrt_pi/(256*half_pow_T(1:moreCount,5))) &
      & * erf_T(1:moreCount) - exp_T(1:moreCount,2) &
      & * (24.0d0/(40320d0*more_T(1:moreCount)) &
      & + 6.0d0/((4)*720d0*pow_T(1:moreCount,2,2)) &
      & + 2.0d0/((16)*24d0*pow_T(1:moreCount,3,2)) &
      & + 1.0d0/((64)*2d0*pow_T(1:moreCount,4,2))))
 
more_F(1:moreCount,6) = 15120.0d0 &
      & * ((sqrt_pi/(1024*half_pow_T(1:moreCount,6))) &
      & * erf_T(1:moreCount) - exp_T(1:moreCount,2) &
      & * (120.0d0/(3628800d0*more_T(1:moreCount)) &
      & + 24.0d0/((4)*40320d0*pow_T(1:moreCount,2,2)) &
      & + 6.0d0/((16)*720d0*pow_T(1:moreCount,3,2)) &
      & + 2.0d0/((64)*24d0*pow_T(1:moreCount,4,2)) &
      & + 1.0d0/((256)*2d0*pow_T(1:moreCount,5,2))))
 
more_F(1:moreCount,7) = 332640.0d0 &
      & * ((sqrt_pi/(4096*half_pow_T(1:moreCount,7))) &
      & * erf_T(1:moreCount) - exp_T(1:moreCount,2) &
      & * (720.0d0/(479001600d0*more_T(1:moreCount)) &
      & + 120.0d0/((4)*3628800d0*pow_T(1:moreCount,2,2)) &
      & + 24.0d0/((16)*40320d0*pow_T(1:moreCount,3,2)) &
      & + 6.0d0/((64)*720d0*pow_T(1:moreCount,4,2)) &
      & + 2.0d0/((256)*24d0*pow_T(1:moreCount,5,2)) &
      & + 1.0d0/((1024)*2d0*pow_T(1:moreCount,6,2))))
   endif
   

   ! Reassemble the total F from less and more components according to the
   !   recorded indices.
   do i = 1, 7
      lessCount = 0
      moreCount = 0
      do j = 1, numAlphaPairs
         if (record_T(j) == 1) then
            lessCount = lessCount + 1
            F(j,i) = less_F(lessCount,i)
         else
            moreCount = moreCount + 1
            F(j,i) = more_F(moreCount,i)
         endif
      enddo
   enddo

   ! Deallocate temporary variables.
   deallocate(exp_T)
   deallocate(sqrt_T)
   deallocate(pow_T)
   if (lessCount > 0) then
      deallocate(less_F)
   endif
   if (moreCount > 0) then
      deallocate(erf_T)
      deallocate(half_pow_T)
      deallocate(more_F)
   endif

end subroutine boys_vec
"""
    f.write(head)


def print_boys(f):

    # Print the boys function subroutine.
    head = """
   ! This subroutine computes F_N(T) - which is the Incomplete
   ! Gamma Function, aka the Boys Function. This comes about as a 
   ! result of the integral transforms that are required to 
   ! account for the (1/r) term that is present in the Coulomb 
   ! integral, to make the integral separable in the Cartesian 
   ! directions (x,y,z).  The formalism is Eq. (29) on page 41 of the 
   ! following document:
   !
   ! ****************************************************************
   ! Petersson T., Hellsing B., "A Detailed Derivation of Gaussian
   ! Orbital-Based Matrix Elements In Electron Structure Calculations."
   ! European Journal of Physics, Vol. 31, No. 1, Page 37.  
   ! ****************************************************************
   !
   ! The Boys Function itself is defined as:
   ! F_N(T) = integral [from 0 to 1] {(t**(2*N))*exp(-T*(t**2))}dt
   !
   ! where T = (a1 + a2 + a3)*sum(P - C)**2
   !
   ! For our purposes here, N will be an integer >=0 that is a linear
   ! combination of index values that come from the summations that
   ! compute the main body of the Coulomb integral, be it Nuclear
   ! Attraction, or the more complicated Electron Repulsion Integral.
   !
   ! Generally, the solution to this integral is expressed numerically
   ! in terms of standard error functions (erf), but this closed form
   ! expression becomes numerically unstable if T = 0 or is very small 
   ! (~10E-2 or less).  To avoid the large roundoff errors that occur
   ! as a result of this, a check has to be performed on the value of T.  
   ! If T <= 10E-2, the following formalism shown by Cook (referenced below)
   ! needs to be applied (Cook: p. 260):
   !
   ! F_N(T) = (1/2)*exp(-T)*sum[from i=0 to infinity] of
   !                            {gamma(N+1/2)/gamma(N+i+3/2)}
   !
   ! This infinite series is truncated at i = 6 for our desired accuracy,
   ! including for the case where T = 0.
   !
   ! More detailed in formation on the Boys Function and numerical 
   ! approximations of it are detailed in the following sources:
   !
   ! ****************************************************************
   ! Cook, David B., "Handbook of Computational Quantum Chemistry."
   ! Dover Publications (August 22, 2005).
   ! ISBN-10: 0486443078
   ! ISBN-13: 978-0486443072
   ! ****************************************************************
   !
   ! as well as:
   !
   ! ****************************************************************
   ! Helgaker T., Taylor P., "Gaussian Basis Sets and Molecular 
   ! Integrals", Ch. 12 from "Modern Electronic Structure Theory."
   ! Pp. 809 - 815.
   ! ****************************************************************
   subroutine boys(T,F)
   use o_kinds
   use o_constants, only: pi
   implicit none
   
   ! Define passed parameters
   real (kind=double), intent (in) :: T
   real (kind=double), dimension(7), intent (out) :: F

   ! Define local variables
   integer :: N
   real (kind=double) :: sqrt_pi, erf_T, exp_T
   real (kind=double) :: a0, a1, a2, a3, a4, a5, a6, a7
   real (kind=double), dimension(8) :: S

   ! Initialize local variables that are used regardless of the method used
   !   to evaluate the boys function (expansion or analytic).
   exp_T = exp(-T)

   sqrt_pi = sqrt(pi) 
   erf_T = erf(sqrt(T))
   ! As stated earler, if T is small (<= 10E-2), the closed form solution 
   !   to the Boys Function can't be used because it becomes numerically
   !   unstable.  Instead, the following expansion is used to approximate it
   !   for small values (Cook: p. 260).  The loop defines the Nth degree Boys
   !   Function.

   if (T <= 10E-2) then
      ! The values (a0 through a7) are computed gamma(N + 0.5) values
      ! to be plugged into the below approximation of the boys Function:

      ! N = 0

      do N = 0, 6
         !
         ! S(N + 1) = 0.5d0*exp_T*(gamma(N + 0.5d0)/gamma(N + 1.5d0)&
         !   &+ (gamma(N + 0.5d0)*T)/gamma(N + 2.5d0)& 
         !   &+ (gamma(N + 0.5d0)*T**2)/gamma(N + 3.5d0)& 
         !   &+ (gamma(N + 0.5d0)*T**3)/gamma(N + 4.5d0)&
         !   &+ (gamma(N + 0.5d0)*T**4)/gamma(N + 5.5d0)& 
         !   &+ (gamma(N + 0.5d0)*T**5)/gamma(N + 6.5d0)& 
         !   &+ (gamma(N + 0.5d0)*T**6)/gamma(N + 7.5d0))

         ! This approximation is explained in the above documentation,
         ! and is required when T <= 10E-2..
    if (N == 0) then
      a0 = 1.7724538509055161     
      a1 = 0.88622692545275805     
      a2 = 1.3293403881791370     
      a3 = 3.3233509704478426     
      a4 = 11.631728396567450     
      a5 = 52.342777784553533     
      a6 = 287.88527781504416     
      a7 = 1871.2543057977896     
    else if (N == 1) then
      a0 = 0.88622692545275805     
      a1 = 1.3293403881791370     
      a2 = 3.3233509704478426     
      a3 = 11.631728396567450     
      a4 = 52.342777784553533     
      a5 = 287.88527781504416     
      a6 = 1871.2543057977896     
      a7 = 14034.407293483402     
    else if (N == 2) then
      a0 = 1.3293403881791370     
      a1 = 3.3233509704478426     
      a2 = 11.631728396567450     
      a3 = 52.342777784553533     
      a4 = 287.88527781504416     
      a5 = 1871.2543057977896     
      a6 = 14034.407293483402     
      a7 = 119292.46199460901     
    else if (N == 3) then
      a0 = 3.3233509704478426     
      a1 = 11.631728396567450     
      a2 = 52.342777784553533     
      a3 = 287.88527781504416     
      a4 = 1871.2543057977896     
      a5 = 14034.407293483402     
      a6 = 119292.46199460901     
      a7 = 1133278.3889487833     
    else if (N == 4) then
      a0 = 11.631728396567450     
      a1 = 52.342777784553533     
      a2 = 287.88527781504416     
      a3 = 1871.2543057977896     
      a4 = 14034.407293483402     
      a5 = 119292.46199460901     
      a6 = 1133278.3889487833     
      a7 = 11899423.083962219     
    else if (N == 5) then
      a0 = 52.342777784553533     
      a1 = 287.88527781504416     
      a2 = 1871.2543057977896     
      a3 = 14034.407293483402     
      a4 = 119292.46199460901     
      a5 = 1133278.3889487833     
      a6 = 11899423.083962219     
      a7 = 136843365.46556622     
    else if (N == 6) then  
      a0 = 287.88527781504416     
      a1 = 1871.2543057977896     
      a2 = 14034.407293483402     
      a3 = 119292.46199460901     
      a4 = 1133278.3889487833     
      a5 = 11899423.083962219     
      a6 = 136843365.46556622     
      a7 = 1710542068.3195722
    end if
    S(N + 1) = 0.5d0*exp_T*a0*((1/a1) + (T/a2) + ((T**2)/a3)&
      &+ ((T**3)/a4) + ((T**4)/a5) + ((T**5)/a6) + ((T**6)/a7))
  end do
  F(1) = S(1)
  F(2) = S(2)
  F(3) = S(3)
  F(4) = S(4)
  F(5) = S(5)
  F(6) = S(6)
  F(7) = S(7)
else
  F(1) = 0.5*((sqrt_pi/(T**(0.5d0)))*erf_T)

  F(2) = ((sqrt_pi/(4*T**(1.5d0)))*erf_T-exp_T*(1/(2d0*T)))
  
  F(3) = 6*((sqrt_pi/(16*T**(2.5d0)))*erf_T-exp_T*(2/(24d0*T)+1/((4)*2d0*T**(2&
  &))))
  
  F(4) = 60*((sqrt_pi/(64*T**(3.5d0)))*erf_T-exp_T*(6/(720d0*T)+2/((4)*24d0*T*&
  &*(2))+1/((16)*2d0*T**(3))))
  
  F(5) = 840*((sqrt_pi/(256*T**(4.5d0)))*erf_T-exp_T*(24/(40320d0*T)+6/((4)*720&
  &d0*T**(2))+2/((16)*24d0*T**(3))+1/((64)*2d0*T**(4))))
  
  F(6) = 15120*((sqrt_pi/(1024*T**(5.5d0)))*erf_T-exp_T*(120/(3628800d0*T)+24/(&
  &(4)*40320d0*T**(2))+6/((16)*720d0*T**(3))+2/((64)*24d0*T**(4))+1/((256)*2d0&
  &*T**(5))))
  
  F(7) = 332640*((sqrt_pi/(4096*T**(6.5d0)))*erf_T-exp_T*(720/(479001600d0*T)+1&
  &20/((4)*3628800d0*T**(2))+24/((16)*40320d0*T**(3))+6/((64)*720d0*T**(4))+2/&
  &((256)*24d0*T**(5))+1/((1024)*2d0*T**(6))))
end if
end subroutine boys
"""
    f.write(head)


def print_test_momentum_ana(conversion, triads, matrix_mm, matrix_ol, f):

    # Print the subroutine header for the analytical portion.
    head = """
   subroutine momentum2CIntgAna(a1,a2,A,B,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)},3), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" \
           + """

   ! Define local variables.
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (3) :: P, PA, PB, d, preFactorMM
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL
   real (kind=double) :: zeta_a_zeta, zeta_b_zeta

   ! Initialize local variables.
   zeta = a1 + a2
   zeta_a_zeta = a1/zeta
   zeta_b_zeta = a2/zeta
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d**2))
   preFactorMM(1) = preFactorOL*2.0d0*a1*PA(1) ! a or b arbitrary"
   preFactorMM(2) = preFactorOL*2.0d0*a1*PA(2) ! a or b arbitrary"
   preFactorMM(3) = preFactorOL*2.0d0*a1*PA(3) ! a or b arbitrary"

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    num_triads = len(triads)
    lib.print_pc(conversion, f, [matrix_ol, matrix_mm],
                 [[1, 1, "_ol"], [5, 3, ""]], num_triads, num_triads)

    # Print the subroutine foot.
    foot = """
   end subroutine momentum2CIntgAna
"""
    f.write(foot)


def print_test_massvel_ana(conversion, triads, matrix_mv, matrix_ke,
                           matrix_ol, f):

    # Print the subroutine header for the analytical portion.
    head = """
   subroutine massvel2CIntgAna(a1,a2,A,B,pc,sh)

   use O_Kinds
   use O_Constants, only: pi, fineStructure

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)}), intent(out) :: sh" + """

   ! Define local variables.
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ke" + """
   real (kind=double), dimension (3) :: P, PA, PB, d
   real (kind=double) :: zeta, inv_2zeta, xi
   real (kind=double) :: preFactorOL, preFactorKE
   real (kind=double) :: preFactor02, preFactor04
   real (kind=double) :: preFactor22, preFactorMV
   real (kind=double) :: inv_2zeta_a, inv_2zeta_b, inv_8m3c2

   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   inv_2zeta_a = 1.0d0 / (2.0d0 * a1)
   inv_2zeta_b = 1.0d0 / (2.0d0 * a2)
   xi = a1 * a2 / zeta
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   d = A - B
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d*d))
   preFactorKE = xi*(3 - 2*xi*sum(d*d))*preFactorOL
   preFactor02 = sum(PB(:)**2) + 3.0d0/(2.0d0*zeta)
   preFactor04 = sum(PB(:)**4) + sum(PB(:)**2)*3.0d0/zeta &
         & + 9.0d0/(4.0d0*zeta**2)
   preFactor22 = (PB(1)*PB(2))**2 + (PB(1)*PB(3))**2 + (PB(2)*PB(3))**2 &
         & + sum(PB(:)**2)/zeta + 3.0d0/(4.0d0*zeta**2)
   preFactorMV = (fineStructure * 0.001d0)**2 / 8.0d0 &
         & * (16*a2**4*preFactor04 - 80*a2**3*preFactor02 + 60*a2**2 &
         & + 32*a2**4 * preFactor22) * preFactorOL
   inv_8m3c2 = (fineStructure * 0.001d0)**2 / 8.0d0

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    num_triads = len(triads)
    lib.print_pc(conversion, f, [matrix_ol, matrix_ke, matrix_mv],
                 [[1, 1, "_ol"], [2, 1, "_ke"], [6, 1, ""]], num_triads,
                 num_triads)

    # Print the subroutine foot.
    foot = """
   end subroutine massvel2CIntgAna
"""
    f.write(foot)


def print_test_dipole_ana(conversion, triads, matrix_dm, matrix_ol, f):

    # Print the subroutine header for the analytical portion.
    head = """
   subroutine dipole3CIntgAna(a1,a2,A,B,C,pc,sh)

   use O_Kinds
   use O_Constants, only: pi

   implicit none

   ! sh(16,16): 1,s; 2,x; 3,y; 4,z; 5,xy; 6,xz; 7,yz; 8,xx-yy;
   ! 9,2zz-xx-yy; 10,xyz; 11,xxz-yyz; 12,xxx-3yyx; 13,3xxy-yyy; 
   ! 14,2zzz-3xxz-3yyz; 15,4zzx-xxx-yyx; 16,4zzy-xxy-yyy

   ! pc(20,20): 1,s; 2,x; 3,y; 4,z; 5,xx; 6,yy; 7,zz; 8,xy; 9,xz;
   ! 10,yz; 11,xyz; 12,xxy; 13,xxz; 14,yyx; 15,yyz; 16,zzx; 17,zzy
   ! 18,xxx; 19,yyy; 20,zzz

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), intent (in) :: a1, a2
   real (kind=double), dimension (3), intent (in) :: A, B, C
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)},3), intent(out) :: pc" + """
   real (kind=double), dimension (""" \
           + f"{len(conversion)},{len(conversion)},3), intent(out) :: sh" + """

   ! Define local variables.
   real (kind=double), dimension (""" \
           + f"{len(triads)},{len(triads)}) :: pc_ol" + """
           
           
           
   real (kind=double), dimension (3) :: P, PA, PB, PC_3C, d, preFactorDM
   real (kind=double) :: zeta, inv_2zeta, xi, preFactorOL, mu
   real (kind=double), dimension (3,3) :: dip_mom


   ! Initialize local variables.
   zeta = a1 + a2
   inv_2zeta = 1.0d0 / (2.0d0 * zeta)
   xi = a1 * a2 / zeta
   d = A - B
   P = (a1*A + a2*B) / zeta
   PA = P - A
   PB = P - B
   PC_3C = P - C
   
   preFactorOL = ((pi/zeta)**1.5)*exp(-xi*sum(d**2))
   
   ! These are the (s|M(mu)|s) integrals
   preFactorDM(1) = preFactorOL*PC_3C(1)
   preFactorDM(2) = preFactorOL*PC_3C(2)
   preFactorDM(3) = preFactorOL*PC_3C(3)

"""
    f.write(head)

    # Print the pc and sh Gaussian terms.
    num_triads = len(triads)
    lib.print_pc(conversion, f, [matrix_ol, matrix_dm],
                 [[1, 1, "_ol"], [7, 3, ""]], num_triads, num_triads)

    # Print the subroutine foot.
    foot = """
   end subroutine dipole3CIntgAna
"""
    f.write(foot)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. In this case, there is no program. This
    #   should only be imported as a module.
    main()
