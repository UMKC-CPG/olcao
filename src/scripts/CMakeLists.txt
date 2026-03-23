# CMake file for OLCAO set of scripts.

# Document the current location
message("At the set of scripts.")

set (SCRIPTS
ab2pf.py addAtoms appendCol atomPosDX
batchChange batchClean batchSub bondAnalysis bond_analysis.py
cif2skl collectSpectra condense condense.py
dead-md.py dump2skl dump2skl.py
element_data.py element_pct.py
gaussOverlapRT
getPlane graspElems.py
hdf5Comp.py
insert3cbo interFit.py
lmp2skl.py
make_veusz_graph.py makeBDB makeBDB.py makeBOND makeColorSYBD makeFittedRhoV
makeinput makePDOS makePotDB makePotDB.py makeReactions makeReactions.py
makeSGDB makeSYBD makeSYBD.py massGather massPrep massSubmit modStruct
makeinput.py modStruct.py
nanoTube.py normalize
olcao olcaoCompare olcaoObject.py olcao.py olcaoStructures.py org_radwavefn.py
osrecurintg.py osrecurintgana.py osrecurintgnum.py osrecurintglib.py
osrecurintg_makenum.py
pdb2skl pdb2skl.py plotgraph plotgraph.py polcao pot2plot.py processPOPTC
processPOPTC.py
rpdf runCheck runIsoAtoms runIsoAtoms.py
scanSpectra selectSYBD skl2isaacs skl2isaacs.py skl2lmp skl2pdb skl2pdb.py
skl2vasp skl2vasp.py spectra2ODX struct2skl.py structure_control.py
unpackOLCAODB unpackOLCAODB.py uolcao uolcao.py upolcao
vasp2skl vasp2skl.py viewCell.py
xcMesh xyz2skl xyz2skl.py)

set (MODULES AngleData.pm BondData.pm ElementData.pm StructureControl.pm)

set (RC_FILES bond_analysisrc.py condenserc.py dump2sklrc.py makeBDBrc.py makeinputrc.py makePotDBrc.py makeReactionsrc.py makeSYBDrc.py modStructrc.py graspElemsrc.py pdb2sklrc.py processPOPTCrc.py runIsoAtomsrc.py skl2isaacsrc.py skl2pdbrc.py skl2vasprc.py struct2sklrc.py uolcaorc.py vasp2sklrc.py xyz2sklrc.py)

install (PROGRAMS ${SCRIPTS} DESTINATION bin)
install (PROGRAMS ${MODULES} DESTINATION "bin/perl5")
install (PROGRAMS ${RC_FILES} DESTINATION ".olcao")
