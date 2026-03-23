# CMake file for OLCAO set of scripts.

# Document the current location
message("At the set of scripts.")

set (SCRIPTS
   ab2pf.py addAtoms appendCol atomPosDX
   batchChange batchClean batchSub bondAnalysis bond_analysis.py
   cif2skl collectSpectra condense condense.py
   dead-md.py dump2skl dump2skl.py
   element_data.py element_pct.py
   gaussOverlapRT getPlane graspElems.py
   hdf5Comp.py
   insert3cbo interFit.py
   lmp2skl.py
   make_veusz_graph.py makeBDB makeBDB.py makeBOND makeColorSYBD makeFittedRhoV
   makeinput makePDOS makePotDB makePotDB.py makeReactions makeReactions.py
   makeSGDB makeSYBD makeSYBD.py massGather massPrep massSubmit modStruct
   makeinput.py modStruct.py
   nanoTube.py normalize
   olcao olcaoCompare olcaoObject.py olcao.py olcaoStructures.py
   org_radwavefn.py osrecurintg.py osrecurintgana.py osrecurintgnum.py
   osrecurintglib.py osrecurintg_makenum.py
   pdb2skl pdb2skl.py plotgraph plotgraph.py polcao pot2plot.py processPOPTC
   processPOPTC.py
   rpdf runCheck runIsoAtoms runIsoAtoms.py
   scanSpectra selectSYBD skl2isaacs skl2isaacs.py skl2lmp skl2pdb skl2pdb.py
   skl2vasp skl2vasp.py spectra2ODX struct2skl.py structure_control.py
   unpackOLCAODB unpackOLCAODB.py uolcao uolcao.py upolcao
   vasp2skl vasp2skl.py viewCell.py
   xcMesh xyz2skl xyz2skl.py)

set (MODULES AngleData.pm BondData.pm ElementData.pm StructureControl.pm)

set (RC_FILES
   bond_analysisrc.py
   collectSYBDrc.py condenserc.py
   dump2sklrc.py
   esviewrc.py
   hdf5Comprc.py
   graspElemsrc.py
   insert3cborc.py interFitrc.py
   makeBDBrc.py makeBONDrc.py makeFittedRhoVrc.py
   makeinputrc.py makePDOSrc.py makePotDBrc.py makeReactionsrc.py
   makeSYBDrc.py modStructrc.py
   osrecurintg_makenumrc.py osrecurintgrc.py
   pdb2sklrc.py plotgraphrc.py porosityrc.py pot2plotrc.py processPOPTCrc.py
   runIsoAtomsrc.py
   skl2isaacsrc.py skl2pdbrc.py skl2vasprc.py struct2sklrc.py
   uolcaorc.py
   vasp2sklrc.py viewCellrc.py
   xyz2sklrc.py)

install (PROGRAMS ${SCRIPTS} DESTINATION bin)
install (PROGRAMS ${MODULES} DESTINATION "bin/perl5")

# SPECIAL PROCESS FOR RC FILES: #################################
# Install RC files only if they do not already exist, so that user
# edits made after the initial install are preserved across future
# "make install" runs.  To force an overwrite (e.g. when upstream
# defaults change), reconfigure with:
#   cmake -DOLCAO_FORCE_RC_INSTALL=ON ...
# and re-run "make install".
option(OLCAO_FORCE_RC_INSTALL
    "Overwrite existing RC files on install even if they already exist"
    OFF)

foreach(RC_FILE ${RC_FILES})
    if(OLCAO_FORCE_RC_INSTALL)
        install(PROGRAMS ${RC_FILE} DESTINATION ".olcao")
    else()
        install(CODE "
            set(_src \"${CMAKE_CURRENT_SOURCE_DIR}/${RC_FILE}\")
            set(_dest
                \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/.olcao\")
            set(_dst \"\${_dest}/${RC_FILE}\")
            if(NOT EXISTS \"\${_dst}\")
                file(INSTALL DESTINATION \"\${_dest}\"
                    TYPE PROGRAM
                    FILES \"\${_src}\")
            else()
                message(STATUS
                    \"Skipping existing RC file: \${_dst}\")
            endif()
        ")
    endif()
endforeach()
