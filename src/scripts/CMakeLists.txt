# CMake file for OLCAO set of scripts.

# Document the current location
message("At the set of scripts.")

set (SCRIPTS
ab2pf.py addAtoms appendCol atomPosDX
batchChange batchClean batchSub bondAnalysis
cif2skl collectSpectra condense
dump2skl
ensembleCondense
gaussOverlapRT
getPlane graspElems.py
insert3cbo interFit.py
lmp2skl.py
make_veusz_graph.py makeBDB makeBOND makeColorSYBD makeFittedRhoV makeinput
makePDOS makePotDB makeReactions makeSGDB makeSYBD massGather massPrep
massSubmit modStruct
nanoTube.py normalize
olcao olcaoCompare olcaoObject.py olcao.py olcaoStructures.py org_radwavefn.py
osrecurintg.py osrecurintgana.py osrecurintgnum.py osrecurintglib.py
osrecurintg_makenum.py
pdb2skl plotgraph plotgraph.py polcao pot2plot.py processPOPTC
rpdf runCheck runIsoAtoms
scanSpectra selectSYBD skl2isaacs skl2lmp skl2pdb skl2vasp spectra2ODX
unpackOLCAODB uolcao upolcao
vasp2skl
xcMesh xyz2skl)

set (MODULES AngleData.pm BondData.pm ElementData.pm StructureControl.pm)

install (PROGRAMS ${SCRIPTS} DESTINATION bin)
install (PROGRAMS ${MODULES} DESTINATION "bin/perl5")
