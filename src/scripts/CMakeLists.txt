# CMake file for OLCAO set of scripts.

# Document the current location
message("At the set of scripts.")

set (SCRIPTS ab2pf.py addAtoms appendCol atomPosDX
batchChange batchClean batchSub bondAnalysis
cif2olcao collectSpectra condense copyInputs
dump2skl
gaussInt2C gaussInt3C gaussIntKE gaussIntMom gaussIntNucV gaussOverlapRT
insert3cbo intGen
lmp2skl.py
makeBDB makeBOND makeFittedRhoV makeinput makePDOS makePotDB makeReactions
makeSGDB makeSYBD massGather massPrep massSubmit mkdirhier modStruct
nanoTube.py normalize
olcao olcaoCompare olcaoObject.py olcao.py olcaoStructures.py
pdb2skl plotgraph polcao
rpdf runCheck runIsoAtoms
scanSpectra selectSYBD skl2pdb skl2vasp spectra2ODX
unpackOLCAODB
vasp2skl
xcMesh xyz2skl)

set (MODULES AngleData.pm BondData.pm ElementData.pm StructureControl.pm)

install (PROGRAMS ${SCRIPTS} DESTINATION bin)
install (PROGRAMS ${MODULES} DESTINATION "bin/perl5")
