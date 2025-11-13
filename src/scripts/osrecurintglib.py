#!/usr/bin/env python3

# Matrices are printed in the order in which they are provided in the
#   matrix_list array. Each matrix is handled according to whatever special
#   requirements that the matrix may have associated with its icode value.
#   The icode value is the first element in the meta data associated with
#   each matrix. As a note, the second meta data element is the number of
#   terms that the matrix has (i.e., one or three (xyz) usually). The third
#   meta data element is the tag that should be appended to variable names.
# The matrix_meta contains:
#   [0] = icode
#   [1] = num terms (usually 1 or 3 (for xyz))
#   [2] = tag string (usually "_ol" or "_ke")
# icode = 1: Overlap matrix
# icode = 2: Kinetic energy matrix
# icode = 3: Three-center overlap matrix (electronic potential)
# icode = 4: Nuclear potential
# icode = 5: Mass velocity
# icode = 6: Momentum matrix
# icode = 7: Dipole moment
# icode = 8: Derivative of the kinetic energy
# icode = 9: Derivative of the electronic potential
# icode = 10: Derivative of the nuclear potential.
def print_production_pc_sh(conversion, triads, f, matrix_list, matrix_meta,
                           lam_sh_list, lam_pc_list, vectorize):

    if (vectorize):

        for n in range(len(matrix_list)):

            # The nuclear potential matrix is dealt with distinctly.
            if (matrix_meta[n][0] == 4):
                print_nuclear_pc_vec(f, matrix_list[n], lam_pc_list,
                                     len(triads))
            else:

                # If the matrix is identified as any of the kind that requires
                #   multiple indices (i.e., one for each xyz axis)
                for xyz in range(matrix_meta[n][1]):
                    for i in range(len(triads)):
                        for j in range(len(triads)):
                            print_one_pc(f, matrix_list[n], matrix_meta[n],
                                         xyz, lam_pc_list, j, i)

        # Assemble the sh terms for the last element in the matrix list. This
        #   is the only one that matters as a total result when the matrix
        #   list has a length > one. In that case, the previous matrices are
        #   used to assemble the final matrix solution.
        for xyz in range(matrix_meta[n][1]):
            for i in range(len(triads)):
               for j in range(len(triads)):
                   print_one_sh(conversion, f, lam_sh_list, matrix_meta[n],
                                xyz, j, i)
        '''
        # Print the leading overlap matrix pc terms if necessary.
        if (matrix2 != ""):
            for i in range(len(triads)):
                for j in range(len(triads)):
                    print_one_pc(f, matrix2, lam_pc_list, 0, j, i, "_ol")

        # Print the target pc terms for the 2C and 3C overlap and KE.
        if (spec_idx == 0):
            for i in range(len(triads)):
                for j in range(len(triads)):
                    print_one_pc(f, matrix1, lam_pc_list, spec_idx, j, i, "")

        # Print the target pc terms for the nuclear matrix.
        if (spec_idx < 0):
            print_nuclear_pc_vec(f, matrix1, lam_pc_list, len(triads))

        # Print the target pc terms for the momentum matrix.
        elif (spec_idx > 0):
            for xyz in range(3):
                for i in range(len(triads)):
                    for j in range(len(triads)):
                        print_one_pc(f, matrix1, lam_pc_list, xyz+1, j, i, "")

        # Assemble the sh terms for the 2C & 3C overlap, the KE, and the NP.
        if (spec_idx <= 0):
            for i in range(len(conversion)):
                for j in range(len(conversion)):
                    print_one_sh(conversion, f, lam_sh_list, spec_idx, j, i)

        # Assemble the sh terms for the momentum matrix.
        elif (spec_idx > 0):
            for xyz in range(3):
                for i in range(len(conversion)):
                    for j in range(len(conversion)):
                        print_one_sh(conversion, f, lam_sh_list, xyz+1, j, i)
        '''

    else:

        f.write("! 0000 0000 = 128,64,32,16  8,4,2,1 for b, a\n\n")

        f.write("\nif (l1l2switch .eq. 34) then ! 0010 0010\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 4, 4)

        f.write("\nelse if (l1l2switch .eq. 17) then ! 0001 0001\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 1, 1)

        f.write("\nelse if (l1l2switch .eq. 33) then ! 0010 0001\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 1, 4)

        f.write("\nelse if (l1l2switch .eq. 65) then ! 0100 0001\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 1, 10)

        f.write("\nelse if (l1l2switch .eq. 129) then ! 1000 0001\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 1, 20)

        f.write("\nelse if (l1l2switch .eq. 18) then ! 0001 0010\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 4, 1)

        f.write("\nelse if (l1l2switch .eq. 66) then ! 0100 0010\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 4, 10)

        f.write("\nelse if (l1l2switch .eq. 130) then ! 1000 0010\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 4, 20)

        f.write("\nelse if (l1l2switch .eq. 20) then ! 0001 0100\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 10, 1)

        f.write("\nelse if (l1l2switch .eq. 36) then ! 0010 0100\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 10, 4)

        f.write("\nelse if (l1l2switch .eq. 68) then ! 0100 0100\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 10, 10)

        f.write("\nelse if (l1l2switch .eq. 132) then ! 1000 0100\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 10, 20)

        f.write("\nelse if (l1l2switch .eq. 24) then ! 0001 1000\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 20, 1)

        f.write("\nelse if (l1l2switch .eq. 40) then ! 0010 1000\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 20, 4)

        f.write("\nelse if (l1l2switch .eq. 72) then ! 0100 1000\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 20, 10)

        f.write("\nelse ! switch 136; 1000 1000\n\n")
        print_pc(conversion, f, matrix_list, matrix_meta, 20, 20)

        f.write("end if\n\n")


def print_one_pc(f, matrix, matrix_meta, xyz, lam_pc_list, pc_idx_a,
                 pc_idx_b):

    # The matrix_meta contains:
    #   [0] = icode
    #   [1] = num terms (usually 1 or 3 (for xyz)).
    #   [2] = tag string

    idx1 = lam_pc_list[pc_idx_a] + 1
    idx2 = lam_pc_list[pc_idx_b] + 1

    # Create a substring that defines the range of the current segment.
    substring = f"segIndices(1,i,{idx1},{idx2}):"
    substring += f"segIndices(2,i,{idx1},{idx2})"

    string = f"do i = 1, numSegs({idx1},{idx2})\n"
    string += f"   do j = segIndices(1,i,{idx1},{idx2}),"
    string += f"segIndices(2,i,{idx1},{idx2})\n"
    f.write(string)

    if (matrix_meta[1] > 1):
        string = f"   pc(j,{pc_idx_a+1},{pc_idx_b+1},{xyz+1})"
        element = matrix[pc_idx_a][pc_idx_b][xyz].replace(":", "j")
    else:
        string = f"   pc{matrix_meta[2]}(j,{pc_idx_a+1},{pc_idx_b+1})"
        element = matrix[pc_idx_a][pc_idx_b].replace(":", "j")
    string += f" = {element}"
    print_cont_string(string, 80, 3, f, False)

    f.write("\n   enddo\nenddo\n\n")


#def print_one_sh(conversion, f, lam_sh_list, spec_idx, sh_idx_a, sh_idx_b):
def print_one_sh(conversion, f, lam_sh_list, matrix_meta, xyz, sh_idx_a,
                 sh_idx_b):

    idx1 = lam_sh_list[sh_idx_a] + 1
    idx2 = lam_sh_list[sh_idx_b] + 1

    # Assign the appropriate array index completion string according to the
    #   value of spec_idx (1, 2, 3 for momentum matrics, 1 for analytic
    #   nuclear potential but not for numerical nuclear potential).
    if (matrix_meta[1] == 1): # Overlap, KE, Nuclear, 3C-electron, mass vel.
        sh_idx = ""
        pc_idx = ""
        if (matrix_meta[0] == 4): # Nuclear
            pc_idx = ",1"
    elif (matrix_meta[1] == 3): # Momentum, dipole moment, Koverlap, derivatives
        sh_idx = f",{xyz+1}"
        pc_idx = f",{xyz+1}"
    '''
    if (spec_idx == 0):
        sh_idx = ""
        pc_idx = ""
    elif (spec_idx == -1):
        sh_idx = ""
        pc_idx = ",1"
    else:
        sh_idx = f",{spec_idx}"
        pc_idx = f",{spec_idx}"
    '''

    # Grab the conversion coefficients and indices for pc terms.
    orb_a_coeff = conversion[sh_idx_a][0]
    orb_a_index = conversion[sh_idx_a][1]
    orb_b_coeff = conversion[sh_idx_b][0]
    orb_b_index = conversion[sh_idx_b][1]
    
    # Create a substring that defines the range of the current segment.
    substring = f"segIndices(1,i,{idx1},{idx2}):"
    substring += f"segIndices(2,i,{idx1},{idx2})"

    # Initialize the output fortran string to "sh(#,#) = " or
    #   "sh(#,#,#) = ".
    string = f"do i = 1, numSegs({idx1},{idx2})\n"
    string += f"   do j = segIndices(1,i,{idx1},{idx2}),"
    string += f"segIndices(2,i,{idx1},{idx2})\n"
    f.write(string)
    string = f"   sh(j,{sh_idx_a+1},{sh_idx_b+1}{sh_idx}) = "

    # Append pc components to the string for each non-zero coeff pair.
    #   A tricky part is that we need to add some unknown number of
    #   terms together. So, the first time that we add a term we do
    #   not have a " + " prefix. Then, all subsequent terms prefix
    #   the term with a " + ".
    first_term_done = False
    for c in range(len(orb_a_coeff)):
        for d in range(len(orb_b_coeff)):

            # Only include terms with non-zero coefficients.
            if (orb_a_coeff[c]*orb_b_coeff[d] != 0):

                # Manage the prefix "+" sign.
                if (first_term_done):
                    if (orb_a_coeff[c]*orb_b_coeff[d] > 0):
                        string += " + "
                    else:
                        string += " - "
                if (not first_term_done):
                    first_term_done = True

                # Append the term with a few small efficiency checks.

                # If the coefficient product is "1", then we don't
                #   need to multiply by one.
                if (abs(orb_a_coeff[c]*orb_b_coeff[d]) == 1):
                    coeff = ""
                else:
                    coeff = f"{abs(orb_a_coeff[c]*orb_b_coeff[d])}*"

                string += f"{coeff}"
                string += f"pc(j,{orb_a_index[c]}," \
                        + f"{orb_b_index[d]}{pc_idx})"

    print_cont_string(string, 80, 3, f, True)

    f.write("   enddo\nenddo\n\n")


# Print the terms needed to compute the primitive Cartesian Gaussian integrals
#   in a production environment. Note that the key result matrix is always
#   given in matrix1. However, sometimes the overlap matrix is a needed
#   ingredient. In such a case it is supplied as matrix2.
def print_pc(conversion, f, matrix_list, matrix_meta, pc_max_idx_a,
             pc_max_idx_b):
#def print_pc(conversion, f, matrix1, matrix2, spec_idx, pc_max_idx_a,
#             pc_max_idx_b):

    # Compute the max indices for the pc-to-sh conversion.
    sh_max_idx_a = map_pc_to_sh(pc_max_idx_a)
    sh_max_idx_b = map_pc_to_sh(pc_max_idx_b)

    for n in range(len(matrix_list)):
        if (matrix_meta[n][0] == 4):
            # Use the special print subroutine for the nuclear matrix.
            print_nuclear_pc(f, matrix_list, matrix_meta[n], pc_max_idx_a,
                             pc_max_idx_b)
        else:
            for xyz in range(matrix_meta[n][1]):
                for b in range(pc_max_idx_b):
                    for a in range(pc_max_idx_a):
                        if (matrix_meta[n][1] > 1):
                            string = f"pc{matrix_meta[n][2]}({a+1},"
                            string += f"{b+1},{xyz+1}) = "
                            string += f"{matrix_list[n][a][b][xyz]}"
                        else:
                            string = f"pc{matrix_meta[n][2]}({a+1},{b+1}) = "
                            string += f"{matrix_list[n][a][b]}"
                        print_cont_string(string, 80, 3, f, True)

        # Use the number of terms from the meta data of the last matrix in the
        #   list to control the printing format. Also, this should only need
        #   to be done for the last matrix in the list.
        if (n == len(matrix_list)-1):
            for xyz in range(matrix_meta[n][1]):
                print_pc_to_sh(conversion, f, matrix_meta[n], xyz, sh_max_idx_a,
                               sh_max_idx_b)

    '''
    # If matrix2 is required, then print those elements first.
    if (matrix2 != ""):
        for b in range(pc_max_idx_b):
            for a in range(pc_max_idx_a):
                string = f"pc_ol({a+1},{b+1}) = {matrix2[a][b]}"
                print_cont_string(string, 80, 3, f, True)

    # If we are doing a matrix without xyz index then just print the one
    #   matrix followed by the print_pc_to_sh terms.
    if (spec_idx == 0):
        for b in range(pc_max_idx_b):
            for a in range(pc_max_idx_a):
                string = f"pc({a+1},{b+1}) = {matrix1[a][b]}"
                print_cont_string(string, 80, 3, f, True)

        # Construct the sh matrix elements from the pc elements.
        print_pc_to_sh(conversion, f, spec_idx, sh_max_idx_a,
                       sh_max_idx_b)

    # If we have xyz components, then print all components. After each component,
    #   print the associated sh terms.
    elif (spec_idx > 0):
        for xyz in range(3):
            for b in range(pc_max_idx_b):
                for a in range(pc_max_idx_a):
                    string = f"pc({a+1},{b+1},{xyz+1}) = " \
                            + f"{matrix1[a][b][xyz]}"
                    print_cont_string(string, 80, 3, f, True)

        for xyz in range(3):
            # Construct the sh matrix elements from the pc elements.
            print_pc_to_sh(conversion, f, xyz+1, sh_max_idx_a,
                           sh_max_idx_b)

    # If we are doing the nuclear matrix (or perhaps in the future other
    #   matrices with multiple "levels") we call a special print subroutine.
    elif (spec_idx == -1):

        # Use the special print subroutine for the nuclear matrix.
        print_nuclear_pc(f, matrix1, pc_max_idx_a, pc_max_idx_b)

        # Construct the sh matrix elements from the pc elements.
        print_pc_to_sh(conversion, f, spec_idx, sh_max_idx_a,
                       sh_max_idx_b)

    '''

def print_nuclear_pc_vec(f, matrix, matrix_meta, lam_pc_list, pc_max_idx):

    # Construct a list of only the terms that are needed from each m level.
    #   This is a bit of a funny algorithm so I'll try to explain what is
    #   happening to reduce confusion. I'm sure that this could be written
    #   in a better way, but until I have a bit more time to create such
    #   an algorithm this one will have to do.
    # There are three stages to the algorithm.
    # (1) First, we initialize a list of the needed pc terms to include
    #   all the terms from the m=0 (first) level.
    # (2) Then, we search through that linear list and examine each of the
    #   matrix element expressions associated with the members in the list.
    #   The examination consists of creating possible additional terms
    #   and then attempting to find those terms in the text of the needed_pc
    #   matrix element expressions. If a match is found and the search term
    #   is not already in the list of needed_pcs then we add it to the list.
    #   Note: A trick is that we iterate in a while loop along a linear
    #   list. We are searching a list from some start to an end and when
    #   we reach the end we reset the start and end points of the search
    #   (because the list might have grown). Once the list does not grow
    #   any more (because we've found all the necessary terms) then we can
    #   quit the loop.
    # (3) Finally, we traverse the list of needed_pc in a funny order. In
    #   each m level we must traverse from lower to higher, but we need to
    #   start with the highest m level going down to the lowest m level.
    #   Hence, we will reverse the m=0 list and each other m level list as
    #   it is added. Then the full list will be reverse before printing.
    #   This is because the higher m terms need to be known to compute many of
    #   the lower m terms, but within one m level many of the lower terms
    #   need to be known to compute higher terms in the same m level.

    # Initialize the m=0 list of needed_pc.
    needed_pc = []
    m = 0
    for b in range(pc_max_idx):
        for a in range(pc_max_idx):
            needed_pc.append([a, b, m])

    # Add the reversed list of needed pcs for the m=0 level to the final list.
    final_pc = needed_pc.copy()[::-1]

    # Search the contents of the needed_pc list for more terms to include.
    current_term = 0
    last_term = len(needed_pc)
    m = 1
    next_m_pc = []
    while (current_term < last_term):
        for b in range(pc_max_idx):
            for a in range(pc_max_idx):
                substring = f"(:,{a+1},{b+1},{m+1})"
                a_idx = needed_pc[current_term][0]
                b_idx = needed_pc[current_term][1]
                m_idx = needed_pc[current_term][2]
                if (matrix[a_idx][b_idx][m_idx].find(substring) != -1 and
                        [a, b, m] not in needed_pc and
                        [a, b, m] not in next_m_pc):
                    next_m_pc.append([a, b, m])

        # Go to the next term in the needed_pc list.
        current_term += 1

        # If the current_term and last_term are equal that means we have
        #   examined all the terms from one particular m level.
        # Now, we reverse the set of next_m_pcs that we found, append them
        #   to the needed_pc list, increase the m level, and then extend
        #   the index of the last term to include the newly appended terms.
        # If we didn't find any terms then we will quit the while loop because
        #   last_term will not go any higher than the current_term.
        if (current_term == last_term):

            # Go to the next m level.
            m += 1

            # Add a reversed version of the new pcs to the final list.
            final_pc = final_pc + next_m_pc[::-1]

            # Add the not reversed new pcs to the needed_pc list to continue
            #   the search.
            needed_pc = needed_pc + next_m_pc

            # Reset the next_m_pc for the search of the next m level.
            next_m_pc = []

            # Update the list position index.
            last_term = len(needed_pc)

    # Print the terms from the list in reverse order.
    final_pc.reverse()
    for term in final_pc:
        a_idx = term[0]
        b_idx = term[1]
        m_idx = term[2]

        idx1 = lam_pc_list[a_idx] + 1
        idx2 = lam_pc_list[b_idx] + 1

        # Create a substring that defines the range of the current segment.
        substring = f"segIndices(1,i,{idx1},{idx2}):"
        substring += f"segIndices(2,i,{idx1},{idx2})"

        string = f"do i = 1, numSegs({idx1},{idx2})\n"
        f.write(string)
        string = f"   pc({substring},{a_idx+1},{b_idx+1},{m_idx+1})"

        element = matrix[a_idx][b_idx][m_idx].replace(":", substring)
        string += f" = {element}\n"
        print_cont_string(string, 80, 3, f, False)

        f.write("enddo\n\n")

# Note that matrix_list will contain either 1 or 2 matrices in it while
#   matrix_meta will contain only the meta data for the first matrix in
#   any case.
def print_nuclear_pc(f, matrix_list, matrix_meta, pc_max_idx_a, pc_max_idx_b):

    # Construct a list of only the terms that are needed from each m level.
    #   This is a bit of a funny algorithm so I'll try to explain what is
    #   happening to reduce confusion. I'm sure that this could be written
    #   in a better way, but until I have a bit more time to create such
    #   an algorithm this one will have to do.
    # (1) First, we initialize a list of the needed pc terms to include
    #   all the terms from the m=0 (first) level. I.e., every single
    #   integral at the m=0 level is required. Therefore, we add every
    #   single m=0 term to the needed_pc list.
    # (2) Then, we search through that linear list and examine each of the
    #   matrix element expressions associated with the members in the list.
    #   The examination consists of creating possible additional terms (of
    #   higher m value than each of the current needed_pc integral terms)
    #   and then attempting to find those terms in the text of the needed_pc
    #   matrix element expressions. If a match is found and the search term
    #   is not already in the list of needed_pcs then we add it to the list.
    #   Note: A trick is that we iterate in a while loop along a linear
    #   list. We are searching a list from some start to an end and when
    #   we reach the end we reset the start and end points of the search
    #   (because the list might have grown). Once the list does not grow
    #   any more (because we've found all the necessary terms) then we can
    #   quit the loop.
    # (3) Then, we traverse the list of needed_pc in a funny order. In
    #   each m level we must traverse from lower to higher, but we need to
    #   start with the highest m level going down to the lowest m level.
    #   Hence, we will reverse the m=0 list and each other m level list as
    #   it is added. Then, the full list will be reversed before printing.
    #   This is because the higher m terms need to be known to compute many of
    #   the lower m terms, but within one m level many of the lower terms
    #   need to be known to compute higher terms in the same m level.
    # (4) Additionally, there are cases where the regular nuclear calculation
    #   results are used to compute the solutions for another (related)
    #   integral (e.g., for force calculations). In that case, we also search
    #   the terms in the other matrix to see if any terms in the nuclear
    #   matrix are required.

    # Initialize the m=0 list of needed_pc.
    needed_pc = []
    m = 0
    for b in range(pc_max_idx_b):
        for a in range(pc_max_idx_a):
            needed_pc.append([a, b, m])
    #for m in range(matrix_meta[1]):
    #    for b in range(pc_max_idx_b):
    #        for a in range(pc_max_idx_a):
    #            needed_pc.append([a, b, m])

    # Add the reversed list of needed pcs for the m=0 level to the final list.
    final_pc = needed_pc.copy()[::-1]
    #print(final_pc)

    # Search the contents of the needed_pc list for more terms to include.
    current_term = 0
    last_term = len(needed_pc)
    m = 1
    next_m_pc = []
    while (current_term < last_term):
        for b in range(pc_max_idx_b):
            for a in range(pc_max_idx_a):
                a_idx = needed_pc[current_term][0]
                b_idx = needed_pc[current_term][1]
                m_idx = needed_pc[current_term][2]
                if ((a_idx == 0) and (b_idx == 0)):
                    substring = f"preFactorN({m+1})"
                else:
                    substring = f"({a+1},{b+1},{m+1})"
                # If the substring is found, not already in the needed_pc
                #   list, and not already in the next_m_pc list, then we
                #   append it to the next_m_pc list.
                #print (a_idx, b_idx, m_idx, substring,
                #       matrix_list[0][a_idx][b_idx][m_idx],
                #       matrix_list[1][a_idx][b_idx][m_idx])
                if (matrix_list[0][a_idx][b_idx][m_idx].find(substring)
                        != -1 and
                        [a, b, m] not in needed_pc and
                        [a, b, m] not in next_m_pc):
                    next_m_pc.append([a, b, m])
                if (len(matrix_list) > 1 and m_idx <= 1 and
                        matrix_list[1][a_idx][b_idx][m_idx].find(substring)
                        != -1 and
                        [a, b, m] not in needed_pc and
                        [a, b, m] not in next_m_pc):
                    next_m_pc.append([a, b, m])

        # Go to the next term in the needed_pc list.
        current_term += 1

        # If the current_term and last_term are equal that means we have
        #   examined all the terms from one particular m level.
        # Now, we reverse the set of next_m_pcs that we found, append them
        #   to the needed_pc list, increase the m level, and then extend
        #   the index of the last term to include the newly appended terms.
        # If we didn't find any terms then we will quit the while loop because
        #   last_term will not go any higher than the current_term.
        if (current_term == last_term):

            # Go to the next m level.
            m += 1

            # Add a reversed version of the new pcs to the final list.
            final_pc = final_pc + next_m_pc[::-1]

            # Add the not reversed new pcs to the needed_pc list to continue
            #   the search.
            needed_pc = needed_pc + next_m_pc

            # Reset the next_m_pc for the search of the next m level.
            next_m_pc = []

            # Update the list position index.
            last_term = len(needed_pc)

    # Print the terms from the list in reverse order.
    final_pc.reverse()
    for term in final_pc:
        a_idx = term[0]
        b_idx = term[1]
        m_idx = term[2]
        string = f"pc{matrix_meta[2]}({a_idx+1},{b_idx+1},{m_idx+1}) = "
        string += f"{matrix_list[0][a_idx][b_idx][m_idx]}"
        print_cont_string(string, 80, 3, f, True)
    

# Assemble the spherical harmonic (sh) Gaussian integrals from the primitive
#   Cartesian (pc) Gaussian integrals. The value of spec_idx will be either
#   0, 1, 2, 3, or -1 for the cases of "no third index required", x, y, z or
#   "no special sh index required, but use index 1 for pc" respectively.
#   Note: spec_idx = special index.
#def print_pc_to_sh(conversion, f, spec_idx, sh_max_idx_a, sh_max_idx_b):
def print_pc_to_sh(conversion, f, matrix_meta, xyz, sh_max_idx_a,
                   sh_max_idx_b):

    # Assign the appropriate array index completion string according to the
    #   value of spec_idx (1, 2, 3 for momentum matrics, 1 for analytic
    #   nuclear potential but not for numerical nuclear potential).
    if (matrix_meta[1] == 1): # Overlap, KE, Nuclear, 3C-electron, mass vel.
        sh_idx = ""
        pc_idx = ""
        if (matrix_meta[0] == 4): # Nuclear
            pc_idx = ",1"
    elif (matrix_meta[1] == 3): # Momentum, dipole moment, Koverlap, derivatives
        sh_idx = f",{xyz+1}"
        pc_idx = f",{xyz+1}"
    '''
    if (spec_idx == 0):
        sh_idx = ""
        pc_idx = ""
    elif (spec_idx == -1):
        sh_idx = ""
        pc_idx = ",1"
    else:
        sh_idx = f",{spec_idx}"
        pc_idx = f",{spec_idx}"
    '''

    # Construct the sh matrix from the pc matrix and the conversion helper
    #   for the range of a and b orbitals requested.
    for b in range(sh_max_idx_b):
        for a in range(sh_max_idx_a):

            # Grab the conversion coefficients and indices for pc terms.
            orb_a_coeff = conversion[a][0]
            orb_a_index = conversion[a][1]
            orb_b_coeff = conversion[b][0]
            orb_b_index = conversion[b][1]

            # Initialize the output fortran string to "sh(#,#) = " or
            #   "sh(#,#,#) = ".
            string = f"sh({a+1},{b+1}{sh_idx}) = "

            # Append pc components to the string for each non-zero coeff pair.
            #   A tricky part is that we need to add some unknown number of
            #   terms together. So, the first time that we add a term we do
            #   not have a " + " prefix. Then, all subsequent terms prefix
            #   the term with a " + ".
            first_term_done = False
            for c in range(len(orb_a_coeff)):
                for d in range(len(orb_b_coeff)):

                    # Only include terms with non-zero coefficients.
                    if (orb_a_coeff[c]*orb_b_coeff[d] != 0):

                        # Manage the prefix "+" sign.
                        if (first_term_done):
                            if (orb_a_coeff[c]*orb_b_coeff[d] > 0):
                                string += " + "
                            else:
                                string += " - "
                        if (not first_term_done):
                            first_term_done = True

                        # Append the term with a few small efficiency checks.

                        # If the coefficient product is "1", then we don't
                        #   need to multiply by one.
                        if (abs(orb_a_coeff[c]*orb_b_coeff[d]) == 1):
                            coeff = ""
                        else:
                            coeff = f"{abs(orb_a_coeff[c]*orb_b_coeff[d])}*"

                        string += f"{coeff}"
                        string += f"pc({orb_a_index[c]}," \
                                + f"{orb_b_index[d]}{pc_idx})"

            print_cont_string(string, 80, 3, f, True)


def print_cont_string(string, max_col, indent, f, tailing_newlines):
    new_string = ""
    seg_length = max_col - 2 # 2 for the beginning and ending &s
    initial = 0
    final = seg_length
    if (len(string) > 1):

        # Compute an initial guess for the number of times the line needs to
        #   be continued.
        num_cont = len(string) // seg_length

        # Apply a correction if the string is exactly divisible by the segment
        #   length. Subtract one in that case.
        if (len(string) % seg_length == 0):
            num_cont -= 1

        # Add each of the continuation portions of the string to new_string.
        for i in range(num_cont):
            new_string += f"{string[initial : final]}&\n&"
            initial += seg_length 
            final += seg_length

        # Add the final portion of the string to the new_string.
        new_string += string[initial : len(string)]
        if (tailing_newlines):
            new_string += "\n\n"

        f.write(new_string)


def map_pc_to_sh(pc_max_idx):

    if (pc_max_idx == 1):
        sh_max_idx = 1
    elif (pc_max_idx == 4):
        sh_max_idx = 4
    elif (pc_max_idx == 10):
        sh_max_idx = 9
    elif (pc_max_idx == 20):
        sh_max_idx = 16
    elif (pc_max_idx == 32):
        sh_max_idx = 25

    return sh_max_idx


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. In this case, there is no program. This
    #   should only be imported as a module.
    main()
