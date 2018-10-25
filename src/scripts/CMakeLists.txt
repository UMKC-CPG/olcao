# CMake file for OLCAO set of scripts.

# Document the current location
message("At the set of scripts.")

set (SCRIPTS ab2pf.py addAtoms appendCol atomPosDX
batchChange batchClean batchSub bondAnalysis
cif2olcao collectSpectra copyInputs
dump2skl
gaussOverlapRT
insert3cbo intGen
lmp2skl.py
makeBDB makeBOND makeFittedRhoV makeinput makePDOS makePotDB makeSGDB
makeSYBD massGather massPrep massSubmit mkdirhier modStruct
nanoTube.py normalize
olcao olcao.py olcaoObject.py olcaoStructures.py
pdb2skl plotgraph polcao
rpdf runIsoAtoms
scanSpectra selectSYBD skl2vasp spectra2ODX
unpackOLCAODB
vasp2skl
xcMesh xyz2skl)

set (MODULES ElementData.pm BondData.pm AngleData.pm StructureControl.pm)

install (PROGRAMS ${SCRIPTS} DESTINATION bin)
install (PROGRAMS ${MODULES} DESTINATION "bin/perl5")
