#!/usr/bin/env python3


def parameters_and_defaults():
    # Note that the "CB" and "BB" below are for the different cases where the
    #   C site and B site are different: "CB"; and when the C site and B site
    #   are the same: "BB".
    param_dict = {
            "prod" : False, # T = production output; F = test output
            "vect" : False, # T (F) = use (!use) vectorization; !prod. only!
            "angmo" : 3, # highest angular mom. (0=s; 1=p; 2=d; 3=f; 4=g; etc.)
            "2COL" : False, # T (F) = do (!do) two-center overlap
            "2CKE" : False, # T (F) = do (!do) two-center kinetic energy
            "3CNP" : False, # T (F) = do (!do) three-center nuclear potential
            "3COL" : False, # T (F) = do (!do) three-center overlap (c=s-type)
            "2CMM" : False, # T (F) = do (!do) two-center momentum matrix
            "2CMV" : False, # T (F) = do (!do) two-center mass velocity
            "2CDM" : False, # T (F) = do (!do) two-center dipole moment
            "2CDKE" : False, # T (F) = do (!do) two-center derivative KE
            "3CDNPCB" : False, # T (F) = do (!do) three-center deriv. NP-CB
            "3CDNPBB" : False, # T (F) = do (!do) three-center deriv. NP-BB
            "3CDNPBC" : False, # T (F) = do (!do) three-center deriv. NP-BC
            "3CDOLCB" : False, # T (F) = do (!do) three-center deriv. OL-CB
            "3CDOLBB" : False, # T (F) = do (!do) three-center deriv. OL-BB
            "3CDOLBC" : False  # T (F) = do (!do) three-center deriv. OL-BC
            }
    return param_dict


if __name__ == '__main__':
    print(parameters_and_defaults())
