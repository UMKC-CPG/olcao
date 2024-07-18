package StructureControl;

our $VERSION = '0.01';

use 5.008002;
use strict;
use warnings;
use POSIX qw(ceil floor);
use Env;
use lib "$OLCAO_BIN/perl5";
use ElementData;
use Math::Complex;
use Math::Trig;
use Math::MatrixReal;
#use Inline (C => Config => cc => 'gcc',
#            ld => 'gcc',
#            inc => '-I/usr/include');
use Inline C => 'DATA',
           NAME => 'StructureControl';

require Exporter;

our @ISA = qw(Exporter);

# Items to export into caller's namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use StructureControl ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw();


# Preloaded methods go here.


# Define the variables used in this package.

# Constants
my $pi = 3.1415926535897932384626433832795;
#my $epsilon = 0.00000000001;
my $epsilon = 0.00001;
my $bigInt = 1000000000;
my $bigReal = 1000000000.0;
my $bohrRad = 0.5291772180;
my $hartree = 27.21138386;

# File names and other names.
my $olcaoskl = "olcao.skl"; # Default olcao skeleton input file name.
my $modelxyz = "model.xyz"; # Default xyz input file name.
my $modelpdb = "model.pdb"; # Default pdb input file name.
my $modelhin = "model.hin"; # Default hin input file name.
my @systemTitle = "empty";  # Name for the project under consideration.
my $spaceGroupDB = "$OLCAO_DATA/spaceDB"; # Space group data location.
my $atomicBDB    = "$OLCAO_DATA/atomicBDB"; # Basis set data location.
my $atomicPDB    = "$OLCAO_DATA/atomicPDB"; # Potential Fn data location.
my $title;

# Atomic positions and attributes
my $numElements;     # The number of unique elements in the system.
my $numAtoms;        # The number of atoms in the system.
my $numAtomsExt;     # The number of relevant atoms in the periodic system.
my $numAtomUniqSpecies; # The number of unique atomic species in the system.
my @atomUniqSpecies; # A unique index number for each atomic species.
my @central2ExtItemMap;  # Maps from central cell atom # to extended list.
my @ext2CentralItemMap;  # Maps from extended atom # to central cell atom #.
my @datSklMap;       # A Map, given dat atom#, gives skl atom#.
my @sklDatMap;       # A Map, given skl atom#, gives dat atom#.
my @moleculeName;    # Name of the molecule to which each atom belongs.
my @moleculeSeqNum;  # A sequence number for the molecule of each atom.
my @residueName;     # Name of the residue to which each atom belongs.
my @residueSeqNum;   # A sequence number for the residue of each atom.
my @atomElementName; # The lowercase element name for each atom in the system.
my @atomElementID;   # An ID number for each atom saying which element it is.
my @elementList;     # List of the unique elements in the system.
my @elementCount;    # List of the # of atoms of each element in the system.
my @atomSpeciesID;   # The species assignments of each atom.
my @numSpecies;      # The number of species for each element.
my @speciesList;     # List of the unique species for each element. (Tag not #)
my @atomicZ;         # Atomic number (Z) of each atom in the system.
my @atomTag;         # Name for the atom usually from an HIN or PDB file.
my @connectionTag;   # List of atoms connected in undetermined way to each atom.
my @coordination;    # Elemental coordination of each atom.
my @coordinationSummary; # Elemental coordination summary for each element.
my @fractABC;        # Fractional coords. of each atom in a,b,c.
my @directABC;       # Direct space coords. of each atom in a,b,c.
my @directXYZ;       # Direct space coords. of each atom in x,y,z.
my @extDirectXYZ;    # Direct space coords. of each atom in x,y,z
                     #   for each replicated cell, structured by cells.
my @extDirectXYZList;# Direct space coords. of each relevant atom in x,y,z
                     #   in the periodic system, in a flat list.
my @extFractABCList; # Fractional coords. of each atom in a,b,c of each
                     #   relevant atom in the periodic system in a flat list.
my @maxPos;          # Maximum x,y,z value of all atoms in the system.
my @minPos;          # Minimum x,y,z value of all atoms in the system.

# General item attributes
my @skipItem1;       # Skip this item in a pair search. [1..#1]
my @skipItem2;       # Skip this item in a pair search. [1..#2]
my @skipPair;        # Skip this item pair in a loop. [1..#1][1..#2]

# Lattice parameters.
my @realLattice;     # Double array of [a,b,c][x,y,z] lattice parameters (A).
my @realLatticeInv;  # Inverse matrix of the real lattice parameters.
my @recipLattice;    # Reciprocal space lattice [a',b',c'][x,y,z] (A).[2pi*Inv]
my @mag;             # Magnitude of each a,b,c vector in @realLattice.
my @angle;           # Alpha(b,c),beta(a,c),gamma(a,b) in @realLattice. Radians
my @angleDeg;        # Alpha(b,c),beta(a,c),gamma(a,b) in @realLattice. Degrees
my @fullCellMag;     # Given cell (full cell) lattice magnitudes.
my @fullCellAngle;   # Given cell (full cell) angles in radians.
my @fullCellAngleDeg;# Given cell (full cell) angles in degrees.
my @magRecip;        # Magnitude of each vector k_a, k_b, k_c in @recipLattice.
my @angleRecip;      # Radians k_alpha(k_b, k_c), etc. in @recipLattice.
my @angleDegRecip;   # Degrees k_alpha(k_b, k_c), etc. in @recipLattice.
my @sineAngle;       # Sine of each angle: alpha, beta, gamma.
my @abcOrder;        # Order in which abc lattice vectors relate to @xyzOrder.
my @xyzOrder;        # Order in which xyz axes are used by @abcOrder;
my $realCellVolume;  # Volume of the real space cell in angstroms.
my $recipCellVolume; # Volume of the reciprocal space cell in angstroms.
my @sortedLatIndices;  # Indices used to sort the angles and lattice magnitudes.


# Brillouin zone descriptors.
my @bZoneSurfaces;   # Array holding surface normal vectors of the BZs.
my @bZoneVertices;   # Array holding vertices of each BZ.


# Extra lattice data and replicated cell data.
my $doFullCell;      # Flag requesting to use the conventional (full) cell.
my $spaceGroup;      # Space group name or number for this system from skel.
my $spaceGroupNum;   # Space group root number for this system.
my $spaceGroupSubNum;# Space group root sub number for this system (a=1, etc.).
my $spaceGroupName;  # Space group name for this system.
my $latticeType;     # A letter designating type face=F, body=I, prim=P, etc.
my @supercell;       # Definition of a supercell request.
my @scMirror;        # Flags to request mirroring of cells in a supercell.
my @negBit;          # Num. of replicated cells needed in the -abc directions.
my @posBit;          # Num. of replicated cells needed in the +abc directions.
my $buffer=10;       # Buffer space between the system and the simulation box
                     #   for the case of non-periodic systems (bio-molecules).

# References to elemental data from the database and implicit data about the
#   current system under study that is extracted from the elemental database.
my $elementNames_ref;
my $covalRadii_ref;
my $valenceCharge_ref;
my $minPotAlpha;
my $maxNumPotAlphas;
my $maxNumAtomAlphas;
my $maxNumValeStates; # Largest number of valence states from any given atom.

# Manipulation data.
my $rotMatrix;    # Three dimensional rotation matrix.

# Control parameters
my $limitDist;    # Distance threshold for taking some requested action.
my $limitDistSqrd;# Distance threshold for taking some requested action squared.

# Variables defining a selecting box.
my $borderZone=0; # Constrain certain calculations to atoms (items) within a
                  #   defined border.  0 = no constraints applied.  1 = atoms
                  #   must be within border, 2 = atoms must be outside border.
my $borderCoords_ref;
                  # Pointer to the coordinate system that has been used to
                  #   define the upper and lower boundaries of the border.
my $borderCoordsExt_ref;
                  # Exactly as above except that this one points to the data
                  #   set containing the extend atom coordinates as well.
my @borderLow;    # Lower border constraints.
my @borderHigh;   # Higher border constraints.

# Radial Pair Distribution Function constraints and results.
my $rpdfSelect=0; # Constrain certain calculations for the rpdf to only a
                  #   specific element pair.
my @rpdfSelectAtoms; # Element names for atom pair to use for rpdf calculation.
my @rpdf;         # Array result of the RPDF calculation.
#my @rpdfDistances; # Array holding all RPDF distances to be plotted (does not
#                  # include any information for skipped pairs)

# Minimal distance matrix data structures and results (makeinput).
my @minDist;      # Minimal distance between atom pairs accounting for PBC.
my @selfMinDist;  # Minimal distance between an atom and its replicated self.

# Bonding data structures and results.
my @covalRadii;   # The covalent radius of each atom.
my @numBonds;     # Number of bonds for each atom.
my @bonded;       # List of central cell atom numbers each atom is bonded to.
my @bondedExt;    # List of extended cell atom numbers each atom is bonded to.
my @bondLengthExt;# Bond lengths for each atom to its bonded pairs.
my @bondAnglesExt;# List of bond angle values and extended cell bonded pair
                  # atom numbers for each atom in the central cell.
my @numBondAngles;# Number of bond angles for each central cell atom.
my $numBondsTotal;# Number of bonds in the entire system.
my @bondTagID;    # Defines which unique ID group each bond belongs to.
my $numUniqueBondTags; # Num unique bonds (sorted lower to higher alphanumeric)
my @uniqueBondTags; # The list of unique tags for all bonds.
my $numAnglesTotal; # Number of bond angles in the entire system.
my $numUniqueAngleTags; # Num unique bond angles (sorted set of tags + angle)
my @uniqueAngleTags; # The list of unique tags for all bond angles.
my @angleTagID;   # Defines which unique ID group each angle belongs to.
my @angleBonded;  # For each atom and angle, the pair of bonded atoms.

# Bond orientational order data.
my @qBar;         # The complex multicomponent orientational value.
my @qOrder;       # The real scalar bond orientational correlation value.
my @Ylm_Coeff;    # Coefficients to the order 6 spherical harmonics.
my $Ylm_l;        # The order of qBar.  Num qBar components=2*$Ylm_l+1

# Spectral scan data structures and results.
my $scanElement="";
my $numScanPoints;
my @scanPoints;
my @extScanDist;
my @extPoreMap;

# Define the functions that return needed variables or, in the case of arrays,
#   references to needed variables.
sub getTitle
   {return $title;}

sub getPi
   {return $pi;}

sub getSystemTitleRef
   {return \@systemTitle;}

sub getSpaceGroupDB
   {return $spaceGroupDB;}

sub getAtomicBDB
   {return $atomicBDB;}

sub getAtomicPDB
   {return $atomicPDB;}

sub getBohrRad
   {return $bohrRad;}

sub getHartree
   {return $hartree;}

sub getNumElements
   {return $numElements;}

sub getNumSpeciesRef
   {return \@numSpecies;}

sub getNumAtoms
   {return $numAtoms;}

sub getNumAtomsExt
   {return $numAtomsExt;}

sub getAtomUniqSpecies
   {return \@atomUniqSpecies;}

sub getDatSklMapRef
   {return \@datSklMap;}

sub getSklDatMapRef
   {return \@sklDatMap;}

sub getAtomicZRef
   {return \@atomicZ;}

sub getCovalRadiiRef
   {return \@covalRadii;}

sub getAtomElementNameRef
   {return \@atomElementName;}

sub getAtomElementIDRef
   {return \@atomElementID;}

sub getElementListRef
   {return \@elementList;}

sub getElementCountRef
   {return \@elementCount;}

sub getAtomSpeciesIDRef
   {return \@atomSpeciesID;}

sub getSpeciesListRef
   {return \@speciesList;}

sub getDoFullCell
   {return $doFullCell;}

sub getRealLatticeRef
   {return \@realLattice;}

sub getRealLatticeInvRef
   {return \@realLatticeInv;}

sub getRecipLatticeRef
   {return \@recipLattice;}

sub getFractABCRef
   {return \@fractABC;}

sub getDirectABCRef
   {return \@directABC;}

sub getDirectXYZRef
   {return \@directXYZ;}

sub getDirectXYZExtListRef
   {return \@extDirectXYZList;}

sub getCoordinationRef
   {return \@coordination;}

sub getCoordinationSummaryRef
   {return \@coordinationSummary;}

sub getMaxXYZRef
   {return \@maxPos;}

sub getMinXYZRef
   {return \@minPos;}

sub getMagRef
   {return \@mag;}

sub getAngleRef
   {return \@angle;}

sub getFullCellMagRef
   {return \@fullCellMag;}

sub getFullCellAngleRef
   {return \@fullCellAngle;}

sub getMagRecipRef
   {return \@magRecip;}

sub getAngleRecipRef
   {return \@angleRecip;}

sub getSineAngleRef
   {return \@sineAngle;}

sub getSpaceGroupNum
   {return $spaceGroupNum;}

sub getSpaceGroupSubNum
   {return $spaceGroupSubNum;}

sub getSpaceGroup
   {return $spaceGroup;}

sub getSpaceGroupName
   {return $spaceGroupName;}

sub getLatticeType
   {return $latticeType;}

sub getNegBitRef
   {return \@negBit;}

sub getPosBitRef
   {return \@posBit;}

sub getRPDFRef
   {return \@rpdf;}

sub getMinDistRef
   {return \@minDist;}

sub getSelfMinDistRef
   {return \@selfMinDist;}

sub getNumBondsRef
   {return \@numBonds;}

sub getNumBondAnglesRef
   {return \@numBondAngles;}

sub getBondingListRef
   {return \@bonded;}

sub getBondingExtListRef
   {return \@bondedExt;}

sub getBondLengthExtRef
   {return \@bondLengthExt;}

sub getBondTagIDRef
   {return \@bondTagID;}

sub getNumBondsTotal
   {return $numBondsTotal;}

sub getUniqueBondTagsRef
   {return \@uniqueBondTags;}

sub getBondAngleExtRef
   {return \@bondAnglesExt;}

sub getAngleBondedRef
   {return \@angleBonded;}

sub getAngleTagIDRef
   {return \@angleTagID;}

sub getNumAnglesTotal
   {return $numAnglesTotal;}

sub getNumUniqueBondTags
   {return $numUniqueBondTags;}

sub getNumUniqueAngleTags
   {return $numUniqueAngleTags;}

sub getUniqueAngleTagsRef
   {return \@uniqueAngleTags;}

sub getExtScanDistRef
   {return \@extScanDist;}

sub getExtPoreMapRef
   {return \@extPoreMap;}

sub getQOrderRef
   {return \@qOrder;}

sub getExt2CentralItemMapRef
   {return \@ext2CentralItemMap;}

sub getCentral2ExtItemMapRef
   {return \@central2ExtItemMap;}

sub getMinPotAlpha
   {return $minPotAlpha;}

sub getMaxNumPotAlphas
   {return $maxNumPotAlphas;}

sub getMaxNumAtomAlphas
   {return $maxNumAtomAlphas;}

sub getMaxNumValeStates
   {return $maxNumValeStates;}

# Define functions to set these variables.

sub setNumAtoms
   {$numAtoms=$_[0];}

sub setLimitDist
{
   $limitDist = $_[0];
   $limitDistSqrd = $limitDist*$limitDist;
}

sub setABCXYZAssignmentOrder
{
   # Define local variables for passed parameters.
   my $abc_ref = $_[0];
   my $xyz_ref = $_[1];

   # Define other local variables.
   my $axis;

   foreach $axis (1..3)
   {
      $abcOrder[$axis] = $abc_ref->[$axis];
      $xyzOrder[$axis] = $xyz_ref->[$axis];
   }
}

sub setBuffer
   {$buffer=$_[0];}


sub setBorder
{
   # Input is [0] = flag, [1] = coordinate type, [2]..[7] = a1,b1,c1,a2,b2,c2
   #   with a1<a2 etc.  The borders can be given in any units (direct space xyz,
   #   direct space abc, or fractional abc). It is up to the user to be
   #   internally consistent with their programs.  (I.e. if the border is set
   #   in fractional coordinates, then don't ask to check whether an item is
   #   is inside it with direct space xyz coordinates.  The coordinate type
   #   variable is intended to make dealing with that issue easier.

   # Declare local variables.
   my $border;

   $borderZone = $_[0];  # 1 = Only atoms within; 2 = Only atoms outside.
   # 1 = directXYZ; 2 = directABC; 3 = fractABC.
   if ($_[1]==1)
   {
      $borderCoords_ref = \@directXYZ;
      $borderCoordsExt_ref = \@extDirectXYZList;
   }
   elsif ($_[1]==2)
      {$borderCoords_ref = \@directABC;}  # At present there is no extDirectABC.
   elsif ($_[1]==3)
   {
      $borderCoords_ref = \@fractABC;
      $borderCoordsExt_ref = \@extFractABCList;
   }
   foreach $border (1..3)
      {$borderLow[$border] = $_[$border+1];}
   foreach $border (1..3)
      {$borderHigh[$border] = $_[$border+4];}
}

sub setRPDFSelectAtoms
{
   $rpdfSelect = 1;
   $rpdfSelectAtoms[1] = $_[0];
   $rpdfSelectAtoms[2] = $_[1];
}

sub setYlm_l
   {$Ylm_l = $_[0];}

sub setXYZMeshPoints
{
   # Declare passed parameters.
   my @numMeshPoints;

   $numMeshPoints[1] = $_[0];
   $numMeshPoints[2] = $_[1];
   $numMeshPoints[3] = $_[2];
   
   # Declare local variables.
   my $xPoint;
   my $yPoint;
   my $zPoint;

   # Initialize all points of the pore map to zero
   foreach $xPoint (1..$numMeshPoints[1])
   {
      foreach $yPoint (1..$numMeshPoints[2])
      {
         foreach $zPoint (1..$numMeshPoints[3])
            {$extPoreMap[$xPoint][$yPoint][$zPoint] = 0;}
      }
   }
}

sub setScanPoints
{
   # Declare local variables.
   my $doMesh;
   my $axis;
   my $point;
   my $aPoint;
   my $bPoint;
   my $cPoint;
   my @numMeshPoints;
   my @delta;
   my @positionDiff;
   my @scanStartXYZ;
   my @scanEndXYZ;

   # Determine if the points list is defined by a line or a mesh.
   $doMesh = $_[0];
   if ($doMesh == 0)
   {
      # Obtain the passed variables describing the end points and the number of
      #   points to use between the end points.
      $numScanPoints = $_[4];
      foreach $axis (1..3)
         {$scanStartXYZ[$axis] = $_[$axis+4];}
      foreach $axis (1..3)
         {$scanEndXYZ[$axis] = $_[$axis+3+4];}


      # Compute the delta necessary in each direction to move from start to end
      #   in exactly $numScanPoints steps.
      foreach $axis (1..3)
      {
         $positionDiff[$axis] = $scanEndXYZ[$axis] - $scanStartXYZ[$axis];
         $delta[$axis] = $positionDiff[$axis] / ($numScanPoints-1);
      }

      # Compute the position of each point along the path.
      foreach $point (1..$numScanPoints)
      {
         foreach $axis (1..3)
            {$scanPoints[$point][$axis] = $scanStartXYZ[$axis] + 
                   ($point-1)*$delta[$axis];}
      }
   }
   else
   {

      # Obtain the passed variables describing the number of mesh points to
      #   use along each axis.
      $numMeshPoints[1] = $_[1];
      $numMeshPoints[2] = $_[2];
      $numMeshPoints[3] = $_[3];

      # Compute the total number of points in this scan.
      $numScanPoints = $numMeshPoints[1]*$numMeshPoints[2]*$numMeshPoints[3];

      # Compute the x,y,z delta necessary along each a,b,c axis.
      #   $delta[a,b,c][x,y,z]
      foreach $axis (1..3) # Loop over x,y,z
      {
         # Compute a axis contribution to x,y,z delta.
         $delta[1][$axis] = $realLattice[1][$axis] / ($numMeshPoints[1]-1);

         # Compute b axis contribution to x,y,z delta.
         $delta[2][$axis] = $realLattice[2][$axis] / ($numMeshPoints[2]-1);

         # Compute c axis contribution to x,y,z delta.
         $delta[3][$axis] = $realLattice[3][$axis] / ($numMeshPoints[3]-1);
      }

      # Compute the position of each point along the path.
      $point = 0;

      foreach $aPoint (0..$numMeshPoints[1]-1)
      {
      foreach $bPoint (0..$numMeshPoints[2]-1)
      {
      foreach $cPoint (0..$numMeshPoints[3]-1)
      {
         $point++;
         foreach $axis (1..3) # Loop over x,y,z
         {
            $scanPoints[$point][$axis] = $delta[1][$axis]*$aPoint +
                                         $delta[2][$axis]*$bPoint +
                                         $delta[3][$axis]*$cPoint;
         }
      }
      }
      }
   }
}

sub setScanElement
   {$scanElement = lc($_[0]);}

# Reset all data structures in the StructureControl so that another file can
#   be read in and manipulated.
sub reset
{
   $pi = 3.1415926535897932384626433832795;
   $epsilon = 0.00001;
   $bigInt = 1000000000;
   $bigReal = 1000000000.0;
   $bohrRad = 0.5291772180;
   $hartree = 27.21138386;
   $olcaoskl = "olcao.skl";
   $modelxyz = "model.xyz";
   $modelpdb = "model.pdb";
   $modelhin = "model.hin";
   undef (@systemTitle);
   @systemTitle = "empty";
   $spaceGroupDB = "$OLCAO_DATA/spaceDB";
   $atomicBDB    = "$OLCAO_DATA/atomicBDB";
   $atomicPDB    = "$OLCAO_DATA/atomicPDB";
   $title = "";
   $numElements = 0;
   $numAtoms = 0;
   $numAtomsExt = 0;
   undef (@central2ExtItemMap);
   undef @ext2CentralItemMap;
   undef @datSklMap;
   undef @sklDatMap;
   undef @moleculeName;
   undef @moleculeSeqNum;
   undef @residueName;
   undef @residueSeqNum;
   undef @atomElementName;
   undef @atomElementID;
   undef @elementList;
   undef @elementCount;
   undef @atomSpeciesID;
   undef @numSpecies;
   undef @atomUniqSpecies;
   undef @speciesList;
   undef @atomicZ;
   undef @atomTag;
   undef @connectionTag;
   undef @coordination;
   undef @coordinationSummary;
   undef @fractABC;
   undef @directABC;
   undef @directXYZ;
   undef @extDirectXYZ;
   undef @extDirectXYZList;
   undef @extFractABCList;
   undef @maxPos;
   undef @minPos;
   undef @skipItem1;
   undef @skipItem2;
   undef @skipPair;
   undef @realLattice;
   undef @realLatticeInv;
   undef @recipLattice;
   undef @mag;
   undef @angle;
   undef @angleDeg;
   undef @magRecip;
   undef @angleRecip;
   undef @angleDegRecip;
   undef @abcOrder;
   undef @xyzOrder;
   undef $realCellVolume;
   undef $recipCellVolume;
   undef $doFullCell;
   undef $spaceGroup;
   undef $spaceGroupNum;
   undef $spaceGroupSubNum;
   undef $spaceGroupName;
   undef $latticeType;
   undef @supercell;
   undef @scMirror;
   undef @sineAngle;
   undef @negBit;
   undef @posBit;
   $buffer=10;
   undef $elementNames_ref;
   undef $covalRadii_ref;
   undef $valenceCharge_ref;
   undef $minPotAlpha;
   undef $maxNumPotAlphas;
   undef $maxNumAtomAlphas;
   undef $maxNumValeStates;
   undef $rotMatrix;
   undef $limitDist;
   undef $limitDistSqrd;
   $borderZone=0;
   undef $borderCoords_ref;
   undef $borderCoordsExt_ref;
   undef @borderLow;
   undef @borderHigh;
   $rpdfSelect=0;
   undef @rpdfSelectAtoms;
   undef @rpdf;
   undef @minDist;
   undef @selfMinDist;
   undef @covalRadii;
   undef @numBonds;
   undef @bonded;
   undef @bondedExt;
   undef @bondLengthExt;
   undef @bondAnglesExt;
   undef @numBondAngles;
   undef $numBondsTotal;
   undef $numUniqueBondTags;
   undef @bondTagID;
   undef @uniqueBondTags;
   undef $numAnglesTotal;
   undef $numUniqueAngleTags;
   undef @uniqueAngleTags;
   undef @angleTagID;
   undef @angleBonded;
   undef @qBar;
   undef @qOrder;
   undef @Ylm_Coeff;
   undef $Ylm_l;
   $scanElement="";
   undef $numScanPoints;
   undef @scanPoints;
   undef @extScanDist;
   undef @extPoreMap;
}


# This subroutine will read the species used in the OLCAO program as recorded
#   in the olcao.fract-mi file.  This file was created from the makeinput
#   script.
sub readOLCAOSpecies
{
   # Declare local variables.
   my $structFile = "inputs/olcao.fract-mi";
   my $line="";
   my @values;
   my $tempName;
   my $atom;

   # This subroutine requires the dat to skl map.
   &readDatSklMap;

   open (OLCAO,"<$structFile") || die "Cannot open $structFile for reading\n";

   # Read past the lattice parameters and other unnecessary lines.
   $line = <OLCAO>;
   $line = <OLCAO>;
   $line = <OLCAO>;
   $line = <OLCAO>;
   $line = <OLCAO>;
   $line = <OLCAO>;
   $line = <OLCAO>;

   # Read through the atom positions and stop when the index number resets.
   #   This is the symbol that we are at the potential positions.
   foreach $atom (1..$numAtoms)
   {
      # Prepare the current line.
      @values = &prepLine(\*OLCAO,"",'\s+');

      # Extract the species number.
      $tempName = $values[0];
      @values = split(/[a-zA-Z]+/,$tempName);

      # Record the species ID.
      $atomSpeciesID[$datSklMap[$atom]] = $values[1];
   }
   close (OLCAO);
}

sub readDatSklMap
{
   # Declare local variables.
   my @values;
   my $line;
   my $atom;

   # Open the file containing the map;
   open (MAP,"<inputs/datSkl.map") || die "Cannot open datSkl.map.\n";

   # Read past the header.
   $line = <MAP>;

   # Loop to obtain the mapping.  The index number of the datSklMap array
   #   is the dat atom number, and the value stored is the original skel atom
   #   number.  Conversely, the index number of the sklDatMap array is the
   #   skeleton atom number, and the value stored is the dat atom number.
   foreach $atom (1..$numAtoms)
   {
      @values = &prepLine(\*MAP,"",'\s+');
      $datSklMap[$atom] = $values[1];
      $sklDatMap[$values[1]] = $atom;
   }

   # Close the map file.
   close (MAP);
}

sub readInputFile
{
   # Define the needed passed parameters.
   my $inputFile = $_[0];

   # Link to the atomic database.
   &setupDataBaseRef;

   # Decide which type of file to get the structural data from.
   if ($inputFile =~ /olcao/)
      {&readOLCAOSkl(@_);}
   elsif ($inputFile =~ /skl/)
      {&readOLCAOSkl(@_);}
   elsif ($inputFile =~ /pdb/)
      {&readPDB(@_);}
   elsif ($inputFile =~ /xyz/)
      {&readXYZ(@_);}
   elsif ($inputFile =~ /dat/) # Not ideal naming convention. FIX
      {&readStruct(@_);}
   elsif ($inputFile =~ /abc/)
      {&readABC(@_);}
   elsif ($inputFile =~ /hin/)
      {&readHIN(@_);}
   elsif ($inputFile =~ /cif/)
      {&readCIF(@_);}
   elsif ($inputFile =~ /lmp/)
      {&readLMP(@_);}
   else
   {
      print STDOUT "Unknown file type $inputFile.  Aborting\n";
      exit;
   }

   # Map element names to atomic Z number.
   &mapElementNumber;

   # Compute basic information about the system such as the maximum number of
   #   potential terms for any one potential site in the system, the maximum
   #   number of Gaussians in an atomic basis function of any one atom in the
   #   system.
   &computeImplicitInfo;
}

sub setupDataBaseRef
{
   # Initialize the element data from the database.
   ElementData::initElementData;

   # Obtain references to the database.
   $elementNames_ref  = ElementData::getElementNamesRef;
   $covalRadii_ref    = ElementData::getCovalRadiiRef;
#   $valenceCharge_ref = ElementData::getValenceChargeRef;
}


sub readOLCAOSkl
{
   # Declare passed parameters.
   my $inFile = $_[0];
   my $useFileSpecies = $_[1];

   # Declare local parameters.
   my $axisABC;
   my $axisXYZ;
   my $atom;
   my @values;
   my @tempValues;
   my $coordType;
   my $coordIndex;
   my $fixedLattice;
   my $lineNum; # Which current line number is being processed.
   my $numLines; # Total number of lines in the skeleton file.
   my @skeleton;
   my $gettingTitle;
   my $magAngleForm;

   # This is going to be a bit crude to make it easier to program.  Sorry.

   # The olcao.skl file MUST include a 'title' line. The title ends
   #   when the word 'end' is found on a line by itself. Between the two
   #   the user may include any information (e.g. citations etc.).
   # The file MUST include a line with 'cell' or 'cellxyz'. Optionally,
   #   the user may add 'fixed" after 'cell' or 'cellxyz'. The impact of that
   #   option is discussed below. First, however, if 'cell' is used, then
   #   the line after cell contains the a, b, c, alpha, beta, gamma lattice
   #   parameters in Angstroms and degrees.  If 'cellxyz' is used then the
   #   next three lines contain the lattice vectors in xyz format. The first
   #   line is the x,y,z components of the 'a' vector, the second line is the
   #   x,y,z components of the 'b' vector, and the third line is for 'c'. The
   #   'fixed' option will demand that the cell be used as given. If the
   #   'fixed' option is *not* given then the lattice parameters (a, b, c
   #   magnitudes and alpha, beta, gamma angles) will be automatically
   #   arranged such that a <= b <= c and alpha <= beta <= gamma. Note also
   #   that the if 'cell' is given then the a-axis will be colinear with the
   #   Cartesian x-axis, the b-axis will be in the Cartesian x-y plane, and
   #   the c-axis will point as needed out of the x-y plane. On the other
   #   hand, if the celxyz format is used, then the a, b, c axes will be
   #   oriented as given.
   # The file MUST include a line with either 'frac' or 'cart' that
   #   specifies the form of the atomic coordinates.  Also on the line there 
   #   MUST be a single integer specifying how many atoms are listed in the
   #   file.
   # The file MUST include a list of the atoms in the model equal to the
   #   number provided in the 'frac' or 'cart' line.  Each atom must have a
   #   tag followed by the three spatial coordinates in either fractional or
   #   cartesian coordinates as specified above.  After the coordinates for a
   #   given atom, an optional tag may exist which is used to add a set of
   #   covalently bonded atoms to the atom in that line.  The tags are defined
   #   as colon separated element lists (e.g. "H:H").  If such a tag is found
   #   on a line defining an O atom it will create water.
   # The file MUST include a line with 'space' and the space group # on the
   #   same line.  The space group can be given as either a number (which will
   #   use the standard default origin and unique axis), or you can specify
   #   the exact space group name using the uniform ASCII form for specifying
   #   space groups which was derived from the list given here:  "Uniform ASCII
   #   symbols for space groups, point groups, and crystal systems.  by
   #   P. SUSSE, Poster presented at the 17th General Meeting of the
   #   International Mineralogical Association, Toronto, Canada, Aug. 9 - 14,
   #   1998. s. Abstracts  page A62.  It is also available at the following
   #   web site:  http://www.saint-hilaire.ca/main/space.htm.  The space groups
   #   and alternate origins actually available are derived from the book:
   #   A Hypertext Book of Crystallographic Space Group Diagrams and Tables
   #   Copyright 1997-1999.  Birkbeck College, University of London.  Available
   #   at the web site:  http://img.chem.ucl.ac.uk/sgp/mainmenu.htm.  All
   #   names are from the book.  Notational scheme is from first reference.
   # The file MUST include a line with 'supercell' with three integers on that
   #   line indicating how many times each of the a,b,c axes should be
   #   duplicated. (It is common to use "1 1 1" to generate a single cell.)
   # The file MUST include a line with either 'full' or 'prim'.  This
   #   indicates a request for either a primitive cell or a full
   #   (non-primitive) cell.  The full cell is only attainable if it was given
   #   in the input file.  (i.e. This can convert from full to prim, but not
   #   the other way around.)

   # Open the file for reading.
   open (OLCAO,"<$inFile") || die "Could not open $inFile for reading\n";

   # We are going to read all the information into a temporary variable so that
   #   it can be processed in whatever way that we want afterward.
   @skeleton = <OLCAO>;
   $numLines = $#skeleton + 1;
   close (OLCAO);

   # Extract the title. The title is found between two special lines that must
   #   start with barewords "title" and "end". The extraction is a little bit
   #   tricky because we read a line and increment the line number counter.
   #   So, for a given line, we have to check if we found the "end" before we
   #   check for "title". This way the words "title" and "end" will not be
   #   included as part of the actual title.
   $lineNum = 0;
   $gettingTitle = 0;
   while (@values = &prepLine("",$skeleton[$lineNum++],'\s+'))
   {
      if (($values[0] eq "end") && ($#values == 0))
         {$gettingTitle = 0; last;}

      if ($gettingTitle == 1)
         {push (@systemTitle,$skeleton[$lineNum-1]);}

      if (($values[0] eq "title") && ($#values == 0))
         {$gettingTitle = 1;}

      # If the end of the file is reached and we fail to read the desired data
      #   then we need to quit and report the issue to the user.
      if ($lineNum > $numLines)
      {
         print <<ENDHELP;
Finished reading the skeleton file but did not successfully read the title.
The title must have the word 'title' on a line by itself followed by any
number of lines describing the system (e.g., citations, previous
calculations, other pertanent information, etc.). The title is closed by
giving the word 'end' on a line by itself. Please adjust your input
accordingly and resubmit. Exiting.
ENDHELP
         exit;
      }
   }

   # Extract the space group identifier. The value given here must exactly
   #   match the name of some file (or soft link to a file) in the space
   #   group database. (The database is usually found in
   #   $OLCAO_DIR/share/spaceDB.)
   $lineNum = 0;
   while (@values = &prepLine("",$skeleton[$lineNum++],'\s+'))
   {
      if (($values[0] eq "space") && ($#values == 1))
      {
         # Get the space group indicator (which might be either a full
         #   name in ASCII form or a number with optional center after an
         #   underscore.
         $spaceGroup = $values[1];

         # Open the file and extract a specific space group number and type.
         #   Use the space group ID (name or number) to open the space group
         #   operations file from the space group data base.
         open (SPACE,"<$spaceGroupDB/$spaceGroup") ||
               die "Cannot open $spaceGroup for reading in $spaceGroupDB.\n";

         # Read the first line to get the cell lattice type (primitive or non-
         #   primitive) and the space group name.
         @values = &prepLine(\*SPACE,"",'\s+');
         $latticeType = $values[0];
         $spaceGroupName = $values[1];

         # If there is an alternative definition for this space group it is
         #   labeled in the space group file with an index starting at 'a'.  In
         #   such a case, the file name is extended by "_a" for the complete
         #   name of the space group (i.e. not the number)
         #   (e.g. 227 = Fd3~m_a).
         if (($#values > 1) && ($values[2] =~ /[a-z]/))
            {$spaceGroupName = $spaceGroupName . "_$values[2]";}

         # Read the second line to get the root space group number and sub
         #   number.
         @values = &prepLine(\*SPACE,"",'\s+');
         $spaceGroupNum = $values[0];
         $spaceGroupSubNum = $values[1];

         close (SPACE);

         last;
      }

      # If the end of the file is reached and we fail to read the desired data
      #   then we need to quit and report the issue to the user.
      if ($lineNum > $numLines)
      {
         print <<ENDHELP;
Finished reading the skeleton file but did not successfully read the space
group information. The space group information is given by using the word
'space' on a line followed (on the same line) by the space group. The space
group must be given by using one of the available options from the
$spaceGroupDB directory.  Please adjust your input accordingly and resubmit.
Exiting.
ENDHELP
         exit;
      }
   }

   # Extract the supercell and any optionally included mirror plane numbers.
   $lineNum = 0;
   while (@values = &prepLine("",$skeleton[$lineNum++],'\s+'))
   {
      if (($values[0] eq "supercell") && (($#values == 3) || ($#values == 6)))
      {
         @supercell[1..3] = @values[1..3];

         if ($#values > 3)
            {@scMirror[1..3] = @values[4..6];}
         else
            {@scMirror[1..3] = (0,0,0);}

         # We've read this data so we can quite the while loop.
         last;
      }

      # If the end of the file is reached and we fail to read the desired data
      #   then we need to quit and report the issue to the user.
      if ($lineNum > $numLines)
      {
         print <<ENDHELP;
Finished reading the skeleton file but did not successfully read the supercell
information. The supercell information is given by using the word 'space' on a
line followed (on the same line) by at least three integers. The integers
define the number of duplicate cells to be created along each of the a, b, c
axes. If a single unit cell calculation is desired, then the integers should
be 1 1 1. Optionally, an additional three binary values (0 or 1) may be
appended to turn on mirroring of the cell along the axes designated with a 1.
If these numbers are not given then they are assumed to be 0 by default.
Please adjust your input accordingly and resubmit. Exiting.
ENDHELP
         exit;
      }
   }

   # It is expected that the lattice parameters will be for the conventional
   #   cell of the associated space group. Therefore, we need to extract a
   #   flag that indicates whether or not the conventional cell should be
   #   reduced to a primitive cell.
   $lineNum = 0;
   while (@values = &prepLine("",$skeleton[$lineNum++],'\s+'))
   {
      if (($values[0] eq "full") && ($#values == 0))
         {$doFullCell = 1; last;}
      elsif (($values[0] eq "prim") && ($#values == 0))
         {$doFullCell = 0; last;}

      # If the end of the file is reached and we fail to read the desired data
      #   then we need to quit and report the issue to the user.
      if ($lineNum > $numLines)
      {
         print <<ENDHELP;
Finished reading the skeleton file but did not successfully read the flag
information indicating whether to do a full cell or a primitive cell
calculation. The cell type information is given by using either the word
'full' or 'prim' on a line by itself. Please adjust your input accordingly
and resubmit. Exiting.
ENDHELP
         exit;
      }
   }

   # Extract the cell parameters in either a,b,c,alpha,beta,gamma format or
   #   in abc,xyz vector format. Initially assume that the lattice will be
   #   given in abc,xyz vector format.
   $magAngleForm = 0;
   $lineNum = 0;
   while (@values = &prepLine("",$skeleton[$lineNum++],'\s+'))
   {
      if ((($values[0] eq "cell") || ($values[0] eq "cellxyz")) &&
          ($#values <= 1))
      {
         # It is possible that for some reason the structure as defined in the
         #   skeleton file needs to be reoriented. For example, the structure
         #   may be triclinic and may have been given with a > b and c > b.
         #   Then, perhaps the user would like to have the structure reordered
         #   such that a < b < c. Perhaps there are some other issues with the
         #   angles too. On the line with "cell" or "cellxyz" the user can
         #   also specify information to direct the reordering of the lattice
         #   information.
         # This is the place where such information would be provided and
         #   read in. For now, we will not implement this capability.

         # Read in the lattice parameters depending on whether "cell" or
         #   "cellxyz" is given.
         if($values[0] eq "cell")
         {
            # Record that the magnitudes and angles were given.
            $magAngleForm = 1;

            # Read a,b,c magnitudes and alpha, beta, gamma angles.
            @values   = &prepLine("",$skeleton[$lineNum],'\s+');

            $mag[1]   = $values[0];
            $mag[2]   = $values[1];
            $mag[3]   = $values[2];
            $angleDeg[1] = $values[3];  # Degrees.
            $angleDeg[2] = $values[4];  # Degrees.
            $angleDeg[3] = $values[5];  # Degrees.
            $angle[1] = $pi/180.0 * $values[3];  # Radians.
            $angle[2] = $pi/180.0 * $values[4];  # Radians.
            $angle[3] = $pi/180.0 * $values[5];  # Radians.

            # The above cell parameters are the "operational" cell paramters.
            #   Somewhat problematically, they act as both the full cell and
            #   the primitive cell paramters. I.e., the same variable is used
            #   for both purposes. The values are simply overwritten when the
            #   primitive cell is created. Doing so simplifies some other work,
            #   but obviously it can cause other points of confusion. So, to
            #   help reduce confusion we will also store here the cell
            #   parameters as originally given and we will always consider
            #   these to be "full" cell parameters.
            $fullCellMag[1] = $mag[1];
            $fullCellMag[2] = $mag[2];
            $fullCellMag[3] = $mag[3];
            $fullCellAngleDeg[1] = $angleDeg[1];
            $fullCellAngleDeg[2] = $angleDeg[2];
            $fullCellAngleDeg[3] = $angleDeg[3];
            $fullCellAngle[1] = $angle[1];
            $fullCellAngle[2] = $angle[2];
            $fullCellAngle[3] = $angle[3];
         }
         else # $values[0] eq "cellxyz"
         {
            # Read a,b,c vectors in x,y,z format.
            foreach $axisABC (1..3)
            {
               @values   = &prepLine("",$skeleton[$lineNum++],'\s+');
               unshift @values, "";
               @{$realLattice[$axisABC]} = @values;
            }
         }
         last;
      }

      # If the end of the file is reached and we fail to read the desired data
      #   then we need to quit and report the issue to the user.
      if ($lineNum > $numLines)
      {
         print <<ENDHELP;
Finished reading the skeleton file but did not successfully read the cell
parameter information indicating the lattice vector magnitudes and angles or
the explicit lattice vectors. The cell parameters are given by using either
the word 'cell' or 'cellxyz' on a line by itself. If the word 'cell' is used
then the next line must contain 6 real numbers that define (in order) the
magnitudes of the lattice vectors: a, b, c and the angles between vectors:
alpha, beta, gamma. If the word 'cellxyz' is used then the next three lines
must contain (one after the other) the x, y, z vector representation of each
of the a, b, c vectors. (I.e., the next line contains the x,y,z values for the
'a' vector, then the next line contains the x,y,z values for the 'b' vector,
and then the final line describes the x,y,z values for the 'c' vector. In both
cases the units are Angstroms. Please adjust your input accordingly and
resubmit. Exiting.
ENDHELP
         exit;
      }
   }

   # If the lattice was given using a,b,c magnitudes and alpha,beta,gamma
   #   angles, then we need to compute the a,b,c x,y,z vectors. Otherwise,
   #   we were given the lattice directly in a,b,c x,y,z vectors and now
   #   we need to obtain the a,b,c magnitudes and alpha,beta,gamma angles.
   if ($magAngleForm == 1)
      {&getABCVectors;} # Following W. Setyawan, S. Curtarolo, Comp. Mat. Sci.
            #   vol. 49, pp. 299-312, (2010) as much as possible and extending
            #   from there according to "A Hypertext Book of Crystallographic
            #   Space Group Diagrams and Tables", Copyright 1997-1999.
            #   Birkbeck College, University of London.  Available at the web
            #   site:  http://img.chem.ucl.ac.uk/sgp/mainmenu.htm.
   else
      {&abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);}

   # Obtain the inverse of the real lattice.  This must be done now so that
   #   we can get the abc coordinates of atoms if we are given xyz.  It must
   #   be done again later after applying the spacegroup, supercell, and
   #   any reordering of the lattice parameters.
   &makeLatticeInv(\@realLattice,\@realLatticeInv,0);

   # Obtain the sine function of each angle.
   &getAngleSine(\@angle);

   # Initialize the atomic counter.
   $numAtoms = 0;

   # Extract the atomic positions and apply any sorted lattice order.
   $lineNum = 0;
   while (@values = &prepLine("",$skeleton[$lineNum++],'\s+'))
   {
      if (($values[0] =~ /^frac/) || ($values[0] =~ /^cart/))
      {
         # Save the form that the coordinates are given in.
         $coordType = $values[0];

         # Save the number of atoms in the system
         $numAtoms = $values[1];

         # Read the information for each atom. Note that we have applied here
         #   Thoreau's fourth theory of adaptation which states: That's not a
         #   "bug", that's a feature! If the user specifies a number atoms on
         #   the "frac" or "cart" line (say x) that is less than the number
         #   of atoms actually listed in the skeleton file (say q), then only
         #   the first x atoms will be read in and the q-x atoms at the end
         #   of the list will be ignored.
         # To help the poor user that misses this "feature" we print a notice
         #   here indicating the number of atoms being read in. Hopefully, the
         #   user will notice the discrepancy (if there is one) and will make
         #   a correction.
         print STDOUT "Reading in $numAtoms atoms.\n";
         foreach $atom (1..$numAtoms)
         {
            @values = &prepLine("",$skeleton[$lineNum++],'\s+');

            # Save the atom tag element name and species ID (1 if undefined).
            $atomTag[$atom] = lc($values[0]);
            if ($atomTag[$atom] !~ /[0-9]/)
               {$atomTag[$atom] = $atomTag[$atom] . "1";}

            # If a connection tag exists, then save it. Otherwise empty it.
            if ($#values > 3)
               {$connectionTag[$atom] = lc($values[4]);}
            else
               {$connectionTag[$atom] = "";}

            # Save the atomic coordinates in the form that they are presented
            #   (cart or fract) and then convert them to the other possible
            #   forms.
            if ($coordType =~ /frac/)
            {
               # Initialize the $fractABC[$atom][0] to "";
               $fractABC[$atom][0] = "";

               # Get the fractional a,b,c coordinates.
               foreach $axisABC (1..3)
                  {$fractABC[$atom][$axisABC] = $values[$axisABC];}

               # Compute the direct space coordinates in a,b,c.
               foreach $axisABC (1..3)
               {
                  $directABC[$atom][$axisABC] = 
                        $fractABC[$atom][$axisABC]*$mag[$axisABC];
               }
               # Calculate the atomic coordinates in x,y,z direct space.
               &getDirectXYZ($atom);
            }
            elsif ($coordType =~ /cart/)
            {
               # Get the x,y,z direct space coordinates.
               foreach $axisXYZ (1..3)
                  {$directXYZ[$atom][$axisXYZ] = $values[$axisXYZ];}

               # Convert to direct a,b,c and then to fractional a,b,c.
               &getDirectABC($atom);
               &getFractABC($atom);
            }
         }

         # We have read in all the atom information.
         last;
      }

      # If the end of the file is reached and we fail to read the desired data
      #   then we need to quit and report the issue to the user.
      if ($lineNum > $numLines)
      {
         print <<ENDHELP;
Finished reading the skeleton file but did not successfully read the atomic
coordinate information indicating the types and positions of the symmetry
equivalent atoms in the model. The atomic coordinates are given by using
either the word 'fract' or 'cart' on a line followed (on the same line) by an
integer X equal to the number of symmetry equivalent atoms listed in the
skeleton file. On each of the X subsequent lines one symmetry equivalent atom
type and coordinates are given. The type consists of the letters from the
periodic table of the elements designating the element followed by an integer
number indicating the type number. (E.g., for silicon type 1 you would write
"si1" or "Si1".) The coordinates are given in either a,b,c fractional units
if 'fract' was used or x,y,z Cartesian units (Angstroms) if 'cart' was used.
Take note that for all coordinates they should be specified to at least eight
decimal places to ensure correct application of symmetry operations.
Optionally, after the coordinates (and on the same line) a special tag may be
added for the purpose of adding additional covalently bonded atoms to the
atom given on that line. The tags are defined as colon separated element lists
(e.g. "H:H").  If a tag such as that example is found on a line defining an
oxygen atom then it will create a water molecule at that site. Please adjust
your input accordingly and resubmit. Exiting.
ENDHELP
         exit;
      }
   }


   # Add atoms listed in connection tags.
   &addConnectionAtoms;

   # Create a set of data structures that contain information about the
   #   relationship between atom numbers, elements, and species.  Also make
   #   lists of the elements in the system and species for each element.
   &createElementList;
   &createSpeciesData($useFileSpecies);

   # Apply the space group operations to the defined atoms.
   &applySpaceGroup;

   # Apply the supercell request to the cell.
   &applySupercell(@supercell[1..3],@scMirror[1..3]);

   # Under certain constraints we will now reorder the magnitudes and angles.
   #   To be in agreement with W. Setyawan, S. Curtarolo, Comp. Mat. Sci.,
   #   v. 49, pp. 299-312, (2010) we would try to obtain the ordering of
   #   alpha <= beta <= gamma and a <= b <= c where possible. However, we will
   #   realign with whatever the user specifies. Take note that due to
   #   symmetry, reordering is not needed or is limited for certain crystal
   #   lattices as follows:
   #   Cubic:     No reordering needed because a=b=c and alpha=beta=gamma=90.
   #   Hexagonal: No reordering needed. The basal plane axes will always be
   #              the a and b axes by convention while the c axis will always
   #              be perpindicular to the basal axes.
   #   Trigonal:  No reordering needed because for the non-hexagonal cells we
   #              have a=b=c and alpha=beta=gamma /= 90.
   #   Tetragonal: Reordering is probably unlikely. Usually a=b /= c is the
   #               conventional presentation. However, someone may want an
   #               unconventional presentation at some time so it is still
   #               permitted.
   #   Orthorhombic: Reordering of a,b,c may be requested.
   #   Monoclinic: No reordering required. The a and b axes are, by
   #               convention, always the same. The c axis may be greater or
   #               less than a,b but it is always perpendicular to the a,b
   #               plane. Similar to tetragonal.
   #   Triclinic: Sort angles first and then sort axes by magnitude. In the
   #              event that this is not possible the lattice vectors will be
   #              ordered first by their angles and then by their magnitudes.
#   if ($fixedLattice == 0)
#   {
#      &reorderLatticeParameters;
#
#      # Reobtain the inverse lattice and the reciprocal lattice vectors.
#      &makeLatticeInv(\@realLattice,\@realLatticeInv,0);
#      &makeLatticeInv(\@realLattice,\@recipLattice,1);
#      &abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);
#      &abcAlphaBetaGamma(\@recipLattice,\@magRecip,\@angleRecip,
#                         \@angleDegRecip,0);
#
#      # If the coordinates are not 'fixed' then they should be reordered as
#      #   the lattice parameters were. The tricky issue is that it is
#      #   possible (although probably stupid and therefor unlikely) that
#      #   someone could give the lattice parameters in non-orthogonal abc
#      #   axes and then provide the coordinates in cartesian xyz format.
#      #   The sorting is done according to abc,alpha,beta,gamma and so it
#      #   should be carried out on the abc form first and then propogated.
#      foreach $atom (1..$numAtoms)
#         {&reorderAtomCoordinates($atom);}
#   }

   # Create an interatomic distance list.
#   &makeInteratomicDist;
}


sub readPDB
{
   # Define passed parameters.
   my $inFile = $_[0];
   my $useFileSpecies = $_[1];

   # Declare local parameters.
   my $line;
   my $recordName;
   my $found;
   my @values;
   my @values2;
   my $explicitABC;
   my $axis;
   my $atom;
   my $species;
   my $element;
   my $currElementID;
   my $tempTag;
   my @pdbAtomTag; # An analog to the @atomTag, but specific for PDB files.
   my @pdbSpeciesList;

   # Assume the file will not contain an explicit set of lattice parameters.
   $explicitABC = 0;

   # Initialize a count of the number of atoms in the system.
   $numAtoms = 0;

   # Open the file for reading.
   open (INFILE,"<$inFile") || die "Could not open $inFile for reading\n";

   while ($line = <INFILE>)
   {
      # Prepare the line for usage.
      chomp $line;

      # Extract the record name and replace with a placeholder.
      $recordName = substr($line,0,6,"123456");

      if ($recordName =~ /TER/)
         {last;}
      elsif ($recordName =~ /CRYST1/)
      {
         # Mark that the crystal parameters are given explicitly.
         $explicitABC = 1;

         # Extract the a,b,c and alpha, beta, gamma lattice parameters.  Make
         #   sure that the angles are stored in radians.
         $mag[1]   = substr($line, 6,9,"123456789");
         $mag[2]   = substr($line,15,9,"123456789");
         $mag[3]   = substr($line,24,9,"123456789");
         $angle[1] = $pi/180.0 * substr($line,33,7,"1234567");
         $angle[2] = $pi/180.0 * substr($line,40,7,"1234567");
         $angle[3] = $pi/180.0 * substr($line,47,7,"1234567");

         # Convert the values to a,b,c into x,y,z vector lattice format.
print STDOUT "Should not be here\n";
         &getABCVectors;

         # Compute the sine of the angles.
         &getAngleSine(\@angle);
      }
      elsif ($recordName =~ /ATOM/)
      {
         # Increment the number of atoms in the system.
         $numAtoms++;

         # Extract the pdb atom tag for this atom. This gets one more character
         #   than defined in the PDB standard because some files rudely use this
         #   extra character when they are not supposed to. This grabs from
         #   columns 12 through 16.  Note that column 12 has no specified value
         #   and so could be used by either the atom name or the serial number
         #   field which is before the atom name field.  Thus, for poorly
         #   generated PDB files or those that do not follow the standard this
         #   particular read my get part of a number (if someone used column
         #   12 for the serial number when they shouldn't) and will break.
         # We should also note that the pdbAtomTag is very different than the
         #   atomTag used in other input files. It will never be an element +
         #   species number combination.
         $tempTag = lc(substr($line,11,5,"12345"));

         # We are now going to try to create an atom tag that is an element
         #   name + species number combination. First we need to get rid of
         #   the leading spaces and trailing spaces from the tempTag before
         #   we store it in the pdbAtomTag.
         @values = &prepLine("",$tempTag,'\s+');
         $pdbAtomTag[$numAtoms] = $values[0];
print STDOUT "PDB Atom Tag: $pdbAtomTag[$numAtoms]\n";

         # Now, if the element symbol is provided on this ATOM line of the PDB
         #   file, then we get that name and use it to construct the atomTag.
         #   This is the preferred way to specify the element for this atom.
         #   If the element symbol is *not* provided, then we make a guess as
         #   to the element name on the basis of the pdbAtomTag. Essentially,
         #   if we have to guess, then we will only be certain to guess
         #   correctly if the element symbol is one character. We will fail if
         #   it is two characters. I'm not sure at this point about how to fix
         #   this issue.
         if (length($line) >= 77)
            {$element = lc(substr($line,76,2,"12"));}
         else
            {$element = "";}
         chomp ($element);
print STDOUT "element = $element\n";
         if ($element ne "")
         {
            @values = &prepLine("",$element,'\s+'); # Get rid of spaces.
            @values2 = &prepLine("",$values[0],'[0-9]'); # Get rid of numbers.
            $atomElementName[$numAtoms] = $values2[0];
         }
         else # Or get the element name from the pdbAtomTag.
         {
            # Get the leading character and call that the element name.  NOTE
            #   that this does not work for elements with two letters in their
            #   name because the tag could be something like "CA" and we can't
            #   necessarily tell if this is Calcium or Carbon Type 'A'.
            #   Presently we just assume that there are no atoms like calcium
            #   in the model.
            @values = split(/[0-9]/,$pdbAtomTag[$numAtoms]);
            if ($values[0] eq "") # Get rid of leading spaces (if any).
               {shift @values;}
            $atomElementName[$numAtoms] = lc($values[0]);
         }

         # Extract the residue name for this atom.
         $residueName[$numAtoms] = substr($line,17,3,"123");

         # Extract the residue sequence number for this atom.
         $residueSeqNum[$numAtoms] = substr($line,22,4,"1234");

         # Assume only 1 molecule is present.
         $moleculeSeqNum[$numAtoms] = 1;

         # Extract the direct space x,y,z atomic coordinates.  This is odd to
         #   me since we define the crystal in terms of a,b,c coordinates, but
         #   the format specification says that the atom coordinates are to be
         #   given in orthogonal x,y,z.
         $directXYZ[$numAtoms][1] = substr($line,30,8,"12345678");
         $directXYZ[$numAtoms][2] = substr($line,38,8,"12345678");
         $directXYZ[$numAtoms][3] = substr($line,46,8,"12345678");

         # Assume that the species for this atom is 1 and create a temporary
         #   atomTag based on that assumption. This will be used to construct
         #   other element information later.
         $atomSpeciesID[$numAtoms] = 1;
         $atomTag[$numAtoms] = "$atomElementName[$numAtoms]" . "1";
      }
   }

   # Compute the crystal parameters if they are not given explicitly.
   if ($explicitABC == 0)
      {&computeCrystalParameters;}

   # Convert the directXYZ into directABC and then to fractABC.
   foreach $atom (1..$numAtoms)
   {
      &getDirectABC($atom);
      &getFractABC($atom);
   }

   # Now, if we want to use the species as defined in the pdb file we need to
   #   extract the species number from the uniqueness of the pdbAtomTag values
   #   for each element. This is very much like the createSpeciesData
   #   subroutine except that we don't have valid atomTag values. Instead we
   #   have pdbAtomTag values that don't represent element+IDnumber pairs.
   #   Thus, we are faced with a problem. Do we modify the createSpeciesData
   #   subroutine to accept a pdbAtomTag array as an edge case or do we just
   #   repeat essentially the same subroutine here. I am going to replicate the
   #   subroutine here to create the atomTag array and then just call the
   #   createSpeciesData subroutine as it normally would be called.
   &createElementList;
   if ($useFileSpecies == 1)
   {
      # Initialize the count of the number of unique species for each element.
      foreach $element (1..$numElements)
         {$numSpecies[$element] = 0;}

      foreach $atom (1..$numAtoms)
      {
         # Define the current element ID.
         $currElementID = $atomElementID[$atom];

         # Determine if the pdbAtomTag associated with the element for this
         #   atom already exists.
         $found = 0;
         foreach $species (1..$numSpecies[$currElementID])
         {
            if ($pdbAtomTag[$atom] eq $pdbSpeciesList[$currElementID][$species])
               {$found = $species; last;}
         }

         if ($found == 0)
         {
            $numSpecies[$currElementID]++;
            $pdbSpeciesList[$currElementID][$numSpecies[$currElementID]] =
                  $pdbAtomTag[$atom];
            $atomSpeciesID[$atom] = $numSpecies[$currElementID];
            $atomTag[$atom] = "$atomElementName[$atom]"
                  . "$numSpecies[$currElementID]";
         }
         else
         {
            $atomSpeciesID[$atom] = $found;
            $atomTag[$atom] = "$atomElementName[$atom]" . "$found";
         }
      }
   }
   &createSpeciesData($useFileSpecies);
}


sub readXYZ
{
   # Declare passed parameters.
   my $inFile = $_[0];
   my $useFileSpecies = $_[1];

   # Declare local parameters.
   my $atom;
   my $line;
   my @values;
   my @tempValues;
   my $explicitABC;
   my $axis;

   # Assume the file will not contain an explicit set of lattice parameters.
   $explicitABC = 0;

   # Open the file for reading.
   open (INFILE,"<$inFile") || die "Could not open $inFile for reading\n";

   # Read and prepare the first line that contains the number of atoms.
   @values = &prepLine(\*INFILE,"",'\s+');
   $numAtoms = $values[0];

   # Read the second line that contains information for the title.
   @values = &prepLine(\*INFILE,"",'\s+');
   $title = $values[0];

   # Read each of the atom element names and direct space xyz coordinates.
   foreach $atom (1..$numAtoms)
   {
      @values = &prepLine(\*INFILE,"",'\s+');
      $atomTag[$atom] = lc($values[0]); # Assume: 1 or 2 letter abbreviation
      $atomElementName[$atom] = lc($values[0]);
      $directXYZ[$atom][1] = $values[1];
      $directXYZ[$atom][2] = $values[2];
      $directXYZ[$atom][3] = $values[3];

      # Assume that the species for this atom is 1.
      $atomSpeciesID[$atom] = 1;

      # Assume an unassigned residue name.
      $residueName[$atom] = "-";

      # Assume an unassigned residue sequence number.
      $residueSeqNum[$atom] = 0;

      # Assume only 1 molecule is present.
      $moleculeSeqNum[$atom] = 1;
   }

   # Compute the crystal parameters if they are not given explicitly.
   if ($explicitABC == 0)
      {&computeCrystalParameters;}

   # Convert the directXYZ into directABC and then to fractABC.
   foreach $atom (1..$numAtoms)
   {
      &getDirectABC($atom);
      &getFractABC($atom);
   }

   # Create a set of data structures that contain information about the
   #   relationship between atom numbers, elements, and species.  Also make
   #   lists of the elements in the system and species for each element.
   &createElementList;
   &createSpeciesData($useFileSpecies);
}


sub readStruct
{
   # Declare passed parameters.
   my $inFile = $_[0];
   my $useFileSpecies = $_[1];

   # Declare local parameters.
   my $atom;
   my $line;
   my @values;
   my @tempValues;
   my $axis;

   # Assume an empty title.
   $title = "";

   # Open the file for reading.
   open (INFILE,"<$inFile") || die "Could not open $inFile for reading\n";

   # Read past the CELL_VECTORS header.
   <INFILE>;

   # Read the abc lattice parameters in vector xyz format.
   foreach $axis (1..3) # Loop over abc.
   {
      @values = &prepLine(\*INFILE,"",'\s+');
      $realLattice[$axis][1] = $values[0]*$bohrRad;
      $realLattice[$axis][2] = $values[1]*$bohrRad;
      $realLattice[$axis][3] = $values[2]*$bohrRad;
   }

   # Read past the NUM_ATOM_SITES header.
   <INFILE>;

   # Read and extract the number of atoms.
   @values = &prepLine(\*INFILE,"",'\s+');
   $numAtoms = $values[0];

   # Read past the NUM_TYPE_X_Y_Z_ELEM header.
   <INFILE>;

   # Read all the type identification of each atom, the atomic coordiantes, and
   #    the elemental designation of each atom.

   # Read each of the atom element names and direct space xyz coordinates.
   foreach $atom (1..$numAtoms)
   {
      @values = &prepLine(\*INFILE,"",'\s+');
      $atomElementName[$atom] = lc($values[5]);
      $directXYZ[$atom][1] = $values[2]*$bohrRad;
      $directXYZ[$atom][2] = $values[3]*$bohrRad;
      $directXYZ[$atom][3] = $values[4]*$bohrRad;

      # Store the type number for this atom.
      if ($useFileSpecies == 1)
         {$atomSpeciesID[$atom] = $values[1];}
      else # Assume that the species for this atom is 1.
         {$atomSpeciesID[$atom] = 1;}

      # Construct the atomTag.
      $atomTag[$atom] = "$atomElementName[$atom]"."$atomSpeciesID[$atom]\n";

      # Assume an unassigned residue name.
      $residueName[$atom] = "-";

      # Assume an unassigned residue sequence number.
      $residueSeqNum[$atom] = 0;

      # Assume only 1 molecule is present.
      $moleculeSeqNum[$atom] = 1;
   }

   # Obtain the lattice vectors in a,b,c,alpha,beta,gamma form.
   &abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);

   # Obtain the inverse of the real lattice.  This must be done now so that
   #   we can get the abc coordinates of atoms if we are given xyz.
   &makeLatticeInv(\@realLattice,\@realLatticeInv,0);

   # Obtain the sine function of each angle.
   &getAngleSine(\@angle);

   # Convert the directXYZ into directABC and then to fractABC.
   foreach $atom (1..$numAtoms)
   {
      &getDirectABC($atom);
      &getFractABC($atom);
   }

   # Create a set of data structures that contain information about the
   #   relationship between atom numbers, elements, and species.  Also make
   #   lists of the elements in the system and species for each element.
   &createElementList;
   &createSpeciesData($useFileSpecies);
}


sub readABC
{
   # Declare passed parameters.
   my $inFile = $_[0];
   my $useFileSpecies = $_[1];

   # Declare local parameters.
   my $line;
   my @values;
   my @tempValues;
   my $axis;

   # Open the file for reading.
   open (INFILE,"<$inFile") || die "Could not open $inFile for reading\n";

   # Read and prepare the first line that contains the lattice parameters.
   @values = &prepLine(\*INFILE,"",'\s+');

   # Fill in the lattice parameters from the first line.
   $realLattice[1][1]=$values[0]; #ax
   $realLattice[1][2]=$values[1]; #ay
   $realLattice[1][3]=$values[2]; #az
   $realLattice[2][1]=$values[3]; #bx
   $realLattice[2][2]=$values[4]; #by
   $realLattice[2][3]=$values[5]; #bz
   $realLattice[3][1]=$values[6]; #cx
   $realLattice[3][2]=$values[7]; #cy
   $realLattice[3][3]=$values[8]; #cz

   # Compute the magnitude of the a,b,c vectors, the angles between then, and
   #   the sine of each angle.
   &abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);

   # Initialize a count of the number of atoms.
   $numAtoms = 0;

   # Read the remaining lines to get the atom names and fractional a,b,c
   #   coordinates.
   while ($line = <INFILE>)
   {
      $numAtoms++;

      # Prepare the line for processing.
      @values = &prepLine("",$line,'\s+');

      # Store atom tag (ID given in olcao.skl file) and atom element.
      $atomTag[$numAtoms] = lc($values[0]);
      @tempValues = split(/[0-9]/,"$atomTag[$numAtoms]");
      $atomElementName[$numAtoms] = lc($tempValues[0]);

      # Assume that the species is 1.
      $atomSpeciesID[$numAtoms] = 1;

      # Initialize the $fractABC[$atom][0] to "";
      $fractABC[$numAtoms][0] = "";

      # Get the fractional a,b,c coordinates.
      foreach $axis (1..3)
         {$fractABC[$numAtoms][$axis] = $values[$axis];}

      # Compute the direct space coordinates in a,b,c.
      foreach $axis (1..3)
         {$directABC[$numAtoms][$axis]=$fractABC[$numAtoms][$axis]*$mag[$axis];}

      # Calculate the atomic coordinates in x,y,z direct space.
      &getDirectXYZ($numAtoms);
   }
   close (INFILE);

   # Create a set of data structures that contain information about the
   #   relationship between atom numbers, elements, and species.  Also make
   #   lists of the elements in the system and species for each element.
   &createElementList;
   &createSpeciesData($useFileSpecies);
}

# At present this will only read cif files that are symmetry P1 and in the
#   format created by the makeinput program.  Sorry, perhaps in the future
#   this can be made more general.
sub readCIF
{
   # Declare passed parameters.
   my $inFile = $_[0];
   my $useFileSpecies = $_[1];

   # Declare local parameters.
   my $line;
   my @values;
   my @tempTags;
   my $axisABC;

   # Open the file for reading.
   open (INFILE,"<$inFile") || die "Could not open $inFile for reading\n";

   # Assume space group 1.
   $spaceGroup = "1_a";

   # Get the title.
   @values = &prepLine(\*INFILE,"",'\s+');
   @values = split(/data_/,$values[0]);
   $systemTitle[0] = "$values[1]";

   # Read up to the cell parameters.
   <INFILE>;
   <INFILE>;
   <INFILE>;
   <INFILE>;
   <INFILE>;

   # Read the cell parameters.
   @values = &prepLine(\*INFILE,"",'\s+');
   $mag[1] = $values[1];
   @values = &prepLine(\*INFILE,"",'\s+');
   $mag[2] = $values[1];
   @values = &prepLine(\*INFILE,"",'\s+');
   $mag[3] = $values[1];
   @values = &prepLine(\*INFILE,"",'\s+');
   $angle[1] = $pi/180.0 * $values[1];
   @values = &prepLine(\*INFILE,"",'\s+');
   $angle[2] = $pi/180.0 * $values[1];
   @values = &prepLine(\*INFILE,"",'\s+');
   $angle[3] = $pi/180.0 * $values[1];

   # Convert the values to a,b,c in x,y,z vector format.
   &getABCVectors;

   # Read up to the atomic positions.
   <INFILE>;
   <INFILE>;
   <INFILE>;
   <INFILE>;
   <INFILE>;
   <INFILE>;

   # Initialize the number of atoms.
   $numAtoms = 0;

   # Read the atomic species, elements, and fractional positions.
   while ($line = <INFILE>)
   {
      # Prepare the line for processing.
      @values = &prepLine("",$line,'\s+');

      # Increment the number of atoms in the system.
      $numAtoms++;

      # Save the atom tag element name and species ID (1 if undefined).
      @tempTags = split(/_/,$values[0]);
      $atomTag[$numAtoms] = lc($tempTags[0]);

      # Initialize the $fractABC[$atom][0] to "";
      $fractABC[$numAtoms][0] = "";

      # Get the fractional a,b,c coordinates.
      foreach $axisABC (1..3)
         {$fractABC[$numAtoms][$axisABC] = $values[$axisABC+1];}

      # Compute the direct space coordinates in a,b,c.
      foreach $axisABC (1..3)
      {
         $directABC[$numAtoms][$axisABC] = 
               $fractABC[$numAtoms][$axisABC]*$mag[$axisABC];
      }
      # Calculate the atomic coordinates in x,y,z direct space.
      &getDirectXYZ($numAtoms);
   }

   # Create a set of data structures that contain information about the
   #   relationship between atom numbers, elements, and species.  Also make
   #   lists of the elements in the system and species for each element.
   &createElementList;
   &createSpeciesData($useFileSpecies);
}


sub readHIN
{
   # Declare passed parameters.
   my $inFile = $_[0];
   my $useFileSpecies = $_[1];

   # Declare local parameters.
   my $line;
   my $lineID;
   my $axis;
   my $atom;
   my @values;
   my @tempValues;
   my $explicitABC;
   my $numMolecules;
   my $numResidues;
   my $currentResName;

   # Assume that the file will not contain an explicit set of ABC parameters.
   $explicitABC = 0;

   # Initialize a count of the number of atoms in the system.
   $numAtoms = 0;

   # Initialize a count of the number of molecules in the system.
   $numMolecules = 0;

   # Open the file for reading.
   open (INFILE,"<$inFile") || die "Could not open $inFile for reading\n";

   while ($line = <INFILE>)
   {
      # Prepare the line for processing.
      @values = &prepLine("",$line,'\s+');

      # Extract the line type identifier.
      $lineID = lc($values[0]);

      if (($lineID =~ /mol/) && ($lineID !~ /end/))
      {
         $numMolecules++;
         $numResidues = 0;
      }
      elsif (($lineID =~ /res/) && ($lineID !~ /end/))
      {
         $numResidues++;
         $currentResName = lc($values[2]);
      }
      elsif ($lineID =~ /atom/)
      {
         # Note that some atoms are actually lone pairs:  "Lp".  These are not
         #   really atoms but they occupy a place in the count of atoms in the
         #   molecule in the file.  For this reason the numAtoms will not be
         #   incremented.

         #Go to the next line if the atom is designated as a lone pair electron.
         if (lc("$values[3]") eq "lp")
            {next;}

         $numAtoms++;

         # Save the atom tag.
         $atomTag[$numAtoms] = $values[2];

         # Record the residue name for this atom.
         $residueName[$numAtoms] = $currentResName;

         # Record the residue sequence number for this atom.
         $residueSeqNum[$numAtoms] = $numResidues;

         # Record the molecule sequence number for this atom.
         $moleculeSeqNum[$numAtoms] = $numMolecules;

         # Record the direct space x,y,z atomic coordinates.
         $directXYZ[$numAtoms][1] = $values[7];
         $directXYZ[$numAtoms][2] = $values[8];
         $directXYZ[$numAtoms][3] = $values[9];

         # Record the element name for this atom.
         $atomElementName[$numAtoms] = lc($values[3]);

         # Assume that the species for this atom is 1.
         $atomSpeciesID[$numAtoms] = 1;
      }
   }

   # Compute the crystal parameters if they are not given explicitly.
   if ($explicitABC == 0)
      {&computeCrystalParameters;}

   # Convert the directXYZ into directABC and then to fractABC.
   foreach $atom (1..$numAtoms)
   {
      &getDirectABC($atom);
      &getFractABC($atom);
   }

   # Create a set of data structures that contain information about the
   #   relationship between atom numbers, elements, and species.  Also make
   #   lists of the elements in the system and species for each element.
   &createElementList;
   &createSpeciesData($useFileSpecies);
}


# Read a lammps dump file while also pulling information from an associated
#   lammps data file.
sub readLMP
{
   # Declare passed parameters.
   my $inFile = $_[0];
   my $useFileSpecies = $_[1];

   # Declare local parameters.
   my $line;
   my @values;
   my @tempTags;
   my $axisABC;

   # Open the file for reading.
   open (INFILE,"<$inFile") || die "Could not open $inFile for reading\n";

   # Assume space group 1.
   $spaceGroup = "1_a";

   # Get the title.
   @values = &prepLine(\*INFILE,"",'\s+');
   @values = split(/data_/,$values[0]);
   $systemTitle[0] = "$values[1]";

   # Read up to the cell parameters.
   <INFILE>;
   <INFILE>;
   <INFILE>;
}

sub reorderLatticeParameters
{
   # Define local variables.
   my $abcAxis;
   my $xyzAxis;
   my @tempLattice;

   # Find the order of the angles from smallest to largest via a stable sort.
   use sort 'stable';
   @sortedLatIndices = sort{$angle[$a] <=> $angle[$b]} 1..3;
   @angle = @angle[@sortedLatIndices];
   @angleDeg = @angleDeg[@sortedLatIndices];

   # The cell vector magnitudes must be assigned in agreement with the angles.
   #   Alpha is the angle between b and c axes, beta is the angle between a
   #   and c axes, and gamma is the angle between a and b axes. Therefore the
   #   old axis a must follow the same index as alpha, b must follow beta, and
   #   c must follow gamma.
   @mag = @mag[@sortedLatIndices];

   # From that initial ordering we need to further refine the ordering in the
   #   event that more that one angle is the same but where the magnitudes
   #   differ. For example an orthorhombic cell has alpha=beta=gamma=90 but
   #   the lattice vector magnitudes will be different. We want to sort the
   #   vectors by their magnitudes from least to greatest.
   # First we check if all three angles are the same. If so, then we can just
   #   directly sort the indices for the magnitudes. Otherwise we need to
   #   check the angles in pairs.
   if ((abs($angle[1] - $angle[2]) < $epsilon) &&
       (abs($angle[2] - $angle[3]) < $epsilon))
   {
      @sortedLatIndices = sort{$mag[$a] <=> $mag[$b]} @sortedLatIndices[1..3];
      @mag = @mag[@sortedLatIndices];
   }
   else
   {
      if ((abs($angle[1] - $angle[3]) < $epsilon) && ($mag[1] > $mag[3]))
      {
         ($sortedLatIndices[1],$sortedLatIndices[3]) =
               ($sortedLatIndices[3],$sortedLatIndices[1]);
         ($mag[1],$mag[3]) = ($mag[3],$mag[1]);
      }
      if ((abs($angle[1] - $angle[2]) < $epsilon) && ($mag[1] > $mag[2]))
      {
         ($sortedLatIndices[1],$sortedLatIndices[2]) =
               ($sortedLatIndices[2],$sortedLatIndices[1]);
         ($mag[1],$mag[2]) = ($mag[2],$mag[1]);
      }
      if ((abs($angle[2] - $angle[3]) < $epsilon) && ($mag[2] > $mag[3]))
      {
         ($sortedLatIndices[2],$sortedLatIndices[3]) =
               ($sortedLatIndices[3],$sortedLatIndices[2]);
         ($mag[2],$mag[3]) = ($mag[3],$mag[2]);
      }
   }

   # Now, we have to apply the sorted order to the vector (abc,xyz)
   #   representation of the lattice parameters.
   foreach $xyzAxis (1..3)
   {
      # Store a temporary copy of the current xyz for each abc.
      foreach $abcAxis (1..3)
         {$tempLattice[$abcAxis][$xyzAxis] = $realLattice[$abcAxis][$xyzAxis];}
      
      # For the current xyz, copy the sorted abc out of the tempLattice.
      foreach $abcAxis (1..3)
      {
         $realLattice[$abcAxis][$xyzAxis] =
               $tempLattice[$sortedLatIndices[$abcAxis]][$xyzAxis];
      }
   }

   # At this point, the angles are certain to be in sorted order from least to
   #   greatest. Then, within the constraint of the ordered angles, the
   #   magnitudes are also in sorted order from least to greatest. After that,
   #   we have sorted the abc,xyz representation of the lattice vectors.
   #   Finally, the indices that reflect the final sort are obtained.
}


sub reorderAtomCoordinates
{
   # Give nice names to the passed parameters.
   my $atom = $_[0];

   # Define local variables.
   my @tempABC;
   my $axisABC;

   # Copy the current abc coordinates into a temporary variable.
   foreach $axisABC (1..3)
      {$tempABC[$axisABC] = $fractABC[$atom][$axisABC];}

   # Apply the reordering learned from the reordering of the lattice params.
   foreach $axisABC (1..3)
      {$fractABC[$atom][$axisABC] = $tempABC[$sortedLatIndices[$axisABC]];}

   # Propogate the reordering to all other representations of the coordinates.

   # Compute the direct space coordinates in a,b,c.
   foreach $axisABC (1..3)
   {
      $directABC[$atom][$axisABC] = 
            $fractABC[$atom][$axisABC]*$mag[$axisABC];
   }

   # Calculate the atomic coordinates in x,y,z direct space.
   &getDirectXYZ($atom);
}


sub addConnectionAtoms
{
   # Define local variables.
   my $axis;
   my $atom;
   my @values;
   my @connAtom;
   my $numConnAtoms;
   my @connOneUnit;   # Unit vector for the first connection atom.
   my $initAtomRadius;
   my $connAtomRadius;
   my $initElementName;
   my @connElementNames;
   my $bondLength;
   my $theta;
   my $phi;
   my $phiPrime;
   my @tempDirectXYZ;
   my $beta = 1.91113553; # radians = 109.5 degrees

   foreach $atom (1..$numAtoms)
   {
      # If the list is empty, then skip to the next atom.
      if ($connectionTag[$atom] eq "")
         {next;}

      # Get the list of atoms connected to this atom.
      @connElementNames = &prepLine("",$connectionTag[$atom],':');
      unshift(@connElementNames,"");

      # Get the number of atoms to be added.
      $numConnAtoms = $#connElementNames;

      # If the number of connected atoms is greater than 4 then we abort
      #   because we don't know how to add such a pattern of atoms.
      if ($numConnAtoms > 4)
         {die "Cannot connect more than four atoms for atom number $atom.\n";}

      # Get the element name of the initial atom (the one being connected to).
      @values = &prepLine("","$atomTag[$atom]",'[0-9]+');
      $initElementName = lc($values[0]);

      # Get the covalent radius of the initial atom in angstroms.
      $initAtomRadius = $covalRadii_ref->[&getElementNumber($initElementName)];

      # Treat the first connection atom.  If we reached this point, then there
      #   is always at least one connection atom.  The first connection atom is
      #   always added a distance away equal to the sum of the covalent radii,
      #   and in a direction that is randomly determined in spherical coords.

      # Increment the number of atoms in the system.
      $numAtoms++;

      # Record the atom tag for this atom.
      $atomTag[$numAtoms] = lc($connElementNames[1]) . "1";

      # Get the covalent radius of this connection atom in angstroms.
      $connAtomRadius =
            $covalRadii_ref->[&getElementNumber(lc($connElementNames[1]))];

      # Find the bond distance as the sum of covalent radii of the connection
      #   atom and the initial atom (the one being connected to).  (In a.u.)
      $bondLength = $initAtomRadius + $connAtomRadius;

      # Choose a random theta (0..pi) and phi (0..2pi).
      $theta = rand($pi);
      $phi   = rand(2*$pi);

      # Compute the x,y,z coordinates of the connected atom (including the
      #   translation from the cell origin).
      $directXYZ[$numAtoms][1] = sin($theta)*cos($phi)*$bondLength +
                                 $directXYZ[$atom][1];
      $directXYZ[$numAtoms][2] = sin($theta)*sin($phi)*$bondLength +
                                 $directXYZ[$atom][2];
      $directXYZ[$numAtoms][3] = cos($theta)*$bondLength +
                                 $directXYZ[$atom][3];

      # Convert to direct a,b,c and then to fractional a,b,c.
      &getDirectABC($numAtoms);
      &getFractABC($numAtoms);

      # If we don't have a second connection atom then we skip to the next atom.
      if ($numConnAtoms < 2)
         {next;}

      # Compute information about the first atom that is useful for all the
      #   remaining atoms to be added.

      # Compute a unit vector based on the first atom and translated to the
      #   cell origin.
      $connOneUnit[1] = $directXYZ[$numAtoms-1][1] - $directXYZ[$atom][1];
      $connOneUnit[2] = $directXYZ[$numAtoms-1][2] - $directXYZ[$atom][2];
      $connOneUnit[3] = $directXYZ[$numAtoms-1][3] - $directXYZ[$atom][3];
      foreach $axis (1..3)
      {
         $connOneUnit[$axis] = ($directXYZ[$numAtoms][$axis] -
                                $directXYZ[$atom][$axis]) / $bondLength;
      }

      # Treat the second connection atom.  This atom's position is constrained
      #   by the position of the first.  So, we will find it by:  (1)  Compute
      #   the second atom position just as for the first except with the new
      #   bond length, an angle $theta that is greater by $beta radians, $phi
      #   kept the same, and at the cell origin.  (2)  Rotate this new position
      #   by a random angle $phiPrime about the axis created by a unit vector
      #   in the direction of the first connected atom.  (As if the first atom
      #   where at the cell origin.)

      # Increment the number of atoms in the system.
      $numAtoms++;

      # Record the atom tag for this atom.
      $atomTag[$numAtoms] = lc($connElementNames[2]) . "1";

      # Get the covalent radius of this connection atom.
      $connAtomRadius =
            $covalRadii_ref->[&getElementNumber(lc($connElementNames[2]))];

      # Find the bond distance as the sum of covalent radii of the connection
      #   atom and the initial atom (the one being connected to).
      $bondLength = $initAtomRadius + $connAtomRadius;

      # Compute the second atom position as described in (1) above.
      $tempDirectXYZ[1] = sin($theta+$beta)*cos($phi)*$bondLength;
      $tempDirectXYZ[2] = sin($theta+$beta)*sin($phi)*$bondLength;
      $tempDirectXYZ[3] = cos($theta+$beta)*$bondLength;

      # Determine the random angle 0..2pi that this point should be rotated
      #   about the unit vector.
      $phiPrime = rand(2*$pi);

      # Compute the rotation of this point about the unit vector by $phiPrime.
      &rotateArbAxis(\@tempDirectXYZ,\@connOneUnit,$phiPrime);

      # Shift the point according to the initial atom's position in the cell.
      foreach $axis (1..3)
         {$tempDirectXYZ[$axis] += $directXYZ[$atom][$axis];}

      # Copy the results to the official record.
      foreach $axis (1..3)
         {$directXYZ[$numAtoms][$axis] = $tempDirectXYZ[$axis];}

      # Convert to direct a,b,c and then to fractional a,b,c.
      &getDirectABC($numAtoms);
      &getFractABC($numAtoms);

      # If we don't have a third connection atom then we skip to the next atom.
      if ($numConnAtoms < 3)
         {next;}

      # Treat the third connection atom.  This will be done exactly as for the
      #   second except that the angle of rotation will be $phiPrime + $beta
      #   radians (109.5 degrees).

      # Increment the number of atoms in the system.
      $numAtoms++;

      # Record the atom tag for this atom.
      $atomTag[$numAtoms] = lc($connElementNames[3]) . "1";

      # Get the covalent radius of this connection atom.
      $connAtomRadius =
            $covalRadii_ref->[&getElementNumber(lc($connElementNames[2]))];

      # Find the bond distance as the sum of covalent radii of the connection
      #   atom and the initial atom (the one being connected to).
      $bondLength = $initAtomRadius + $connAtomRadius;

      # Compute the third atom position as described in (1) above for the
      #   second atom.  (It starts the same way.)
      $tempDirectXYZ[1] = sin($theta+$beta)*cos($phi)*$bondLength;
      $tempDirectXYZ[2] = sin($theta+$beta)*sin($phi)*$bondLength;
      $tempDirectXYZ[3] = cos($theta+$beta)*$bondLength;

      # Determine the angle that this point should be rotated about the unit
      #   vector.
      $phiPrime += 2.0944;  # 120 degrees.

      # Compute the rotation of this point about the unit vector by $phiPrime.
      &rotateArbAxis(\@tempDirectXYZ,\@connOneUnit,$phiPrime);

      # Shift the point according to the initial atom's position in the cell.
      foreach $axis (1..3)
         {$tempDirectXYZ[$axis] += $directXYZ[$atom][$axis];}

      # Copy the results to the official record.
      foreach $axis (1..3)
         {$directXYZ[$numAtoms][$axis] = $tempDirectXYZ[$axis];}

      # Convert to direct a,b,c and then to fractional a,b,c.
      &getDirectABC($numAtoms);
      &getFractABC($numAtoms);

      # If we don't have a fourth connection atom then we skip to the next atom.
      if ($numConnAtoms < 4)
         {next;}

      # Treat the fourth connection atom.  This will be done exactly as for the
      #   second except that the angle of rotation will be $phiPrime - $beta
      #   radians.

      # Increment the number of atoms in the system.
      $numAtoms++;

      # Record the atom tag for this atom.
      $atomTag[$numAtoms] = lc($connElementNames[4]) . "1";

      # Get the covalent radius of this connection atom.
      $connAtomRadius =
            $covalRadii_ref->[&getElementNumber(lc($connElementNames[2]))];

      # Find the bond distance as the sum of covalent radii of the connection
      #   atom and the initial atom (the one being connected to).
      $bondLength = $initAtomRadius + $connAtomRadius;

      # Compute the fourth atom position as described in (1) above for the
      #   second atom.  (It starts the same way.)
      $tempDirectXYZ[1] = sin($theta+$beta)*cos($phi)*$bondLength;
      $tempDirectXYZ[2] = sin($theta+$beta)*sin($phi)*$bondLength;
      $tempDirectXYZ[3] = cos($theta+$beta)*$bondLength;

      # Determine the angle that this point should be rotated about the unit
      #   vector.
      $phiPrime += 2.0944;  # 120 degrees

      # Compute the rotation of this point about the unit vector by $phiPrime.
      &rotateArbAxis(\@tempDirectXYZ,\@connOneUnit,$phiPrime);

      # Shift the point according to the initial atom's position in the cell.
      foreach $axis (1..3)
         {$tempDirectXYZ[$axis] += $directXYZ[$atom][$axis];}

      # Copy the results to the official record.
      foreach $axis (1..3)
         {$directXYZ[$numAtoms][$axis] = $tempDirectXYZ[$axis];}

      # Convert to direct a,b,c and then to fractional a,b,c.
      &getDirectABC($numAtoms);
      &getFractABC($numAtoms);
   }
}


# This subroutine will read a given bondAnalysis.bl file and will populate the
#   relevant data structures. For this subroutine, we will not assume that we
#   have access to any other information about the system under study. I.e.,
#   we may not have previously read in any type of associated olcao.skl file.
# This reads BOND LENGTHS.
sub readBondAnalysisBL
{
   # Define passed parameters.
   my $inFile = $_[0];
   my $numAtoms = $_[1];

   # Define local variables.
   my $atom;
   my $bond;
   my $line;
   my $row;
   my $tag;
   my $found;
   my @values;
   my @values2;
   my $atomTag;
   my $bondTag;
   my $atomNumber;
   my $bondIndex;
   my $numBondRows;
   my $numBondsInRow;

   
   # Open and read the bonds into an array of arrays.
   open (BONDS,"<$inFile") || die "Cannot open $inFile for reading.\n";

   # Assume that all atoms have zero bonds.
   foreach $atom (1..$numAtoms)
      {$numBonds[$atom] = 0;}

   # Initialize a count of the total number of bonds in the system.
   $numBondsTotal = 0;

   # Initialize a count of the number of unique bond tags.
   $numUniqueBondTags = 0;

   # Read each set of bond information in the given bondAnalysis file.
   while ($line = <BONDS>)
   {
      # Get the number of bonds, the number of bond rows, and the atom number
      #   of the current atom.
      @values = &prepLine("",$line,'\s+');
      @values2 = split(/_/,"$values[0]");
      $atomTag = $values2[0]; # Something like b1 or c2 or fe2, etc.
      $atomNumber = $values2[1];
      $numBonds[$atomNumber] = $values[2];
      $numBondRows = ceil($values[2]/4.0);
      $bondIndex = 0; # Initialize an indexer for the bonds for this atom.
      foreach $row (1..$numBondRows)
      {
         # Compute the number of bonds in this row (may be 4 or less).
         @values = &prepLine(\*BONDS,"",'\s+');
         $numBondsInRow = ($#values+1)/2;

         # Read the atom number of the bonded pair.
         foreach $bond (1..$numBondsInRow)
         {
            $bondIndex++;
            $numBondsTotal++;
            @values2 = split(/_/,"@values[2*($bond-1)]"); # Like b1, c2, fe2,...
            $bonded[$atomNumber][$bondIndex] = $values2[1];
            $bondLengthExt[$atomNumber][$bondIndex] = @values[2*($bond-1)+1];

            if ("$atomTag" gt "$values2[0]")
               {$bondTag = "$values2[0] $atomTag";}
            else
               {$bondTag = "$atomTag $values2[0]";}

            $found = 0;
            foreach $tag (1..$numUniqueBondTags)
            {
              if ($bondTag eq $uniqueBondTags[$tag])
                 {$found = $tag;}
            }

            # If we didn't find the tag in the list of known tags, then
            #   increase the number of known tags and store the tag.
            if ($found == 0)
            {
               $numUniqueBondTags++;
               $found = $numUniqueBondTags;
               $uniqueBondTags[$found] = $bondTag;
            }

            # For the current bond, store the unique bond ID tag number.
            $bondTagID[$atomNumber][$bondIndex] = $found;
         }
      }
   }

   # Cut the bond count in half because all bonds are double-listed in the
   #   bondAnalysis.bl files.
   $numBondsTotal /= 2;

   close (BONDS);
}


# This subroutine will read a given bondAnalysis.ba file and will populate the
#   relevant data structures. For this subroutine, we will not assume that we
#   have access to any other information about the system under study. I.e.,
#   we may not have previously read in any type of associated olcao.skl file.
# This reads BOND ANGLES.
sub readBondAnalysisBA
{
   # Define passed parameters.
   my $inFile = $_[0];
   my $numAtoms = $_[1];

   # Define local variables.
   my $tag;
   my $atom;
   my $line;
   my $angle;
   my $found;
   my @values;
   my @values2;
   my $atomNumber;
   my $angleTag;
   my $currentBondAngle;
   my $decimal;
   
   # Open and read the bond angles into an array of arrays.
   open (BONDS,"<$inFile") || die "Cannot open $inFile for reading.\n";

   # Assume that all atoms have zero bond angles.
   foreach $atom (1..$numAtoms)
      {$numBondAngles[$atom] = 0;}

   # Initialize a count of the total number of bond angles in the system.
   $numAnglesTotal = 0;

   # Initialize a count of the total number of unique bond angle configurations
   #   defined by the combination of element-species tags and bond angles.
   $numUniqueAngleTags = 0;

   # Read each set of bond angle information in the given bondAnalysis file.
   while ($line = <BONDS>)
   {
      # Get the number of bond angles where this atom serves as a vertex. Also,
      #   increment the number of bond angles in total.
      @values = &prepLine("",$line,'\s+');
      $atomNumber = $values[0];
      $numBondAngles[$atomNumber] = $values[4];
      $numAnglesTotal += $numBondAngles[$atomNumber];

      # Read information about each bond angle.
      foreach $angle (1..$numBondAngles[$atomNumber])
      {
         # Read the next bond angle data set.
         @values = &prepLine(\*BONDS,"",'\s+');

         # Round the bond angle to the nearest integer or half-integer. Take
         #   special note of the different arrangement of parenthesis in the
         #   rounding "int" call.
         $decimal = ($values[$#values] - int($values[$#values]));
         if ($decimal > 0.75)
            {$currentBondAngle = int($values[$#values])+1.0;}
         elsif ($decimal > 0.5)
            {$currentBondAngle = int($values[$#values])+0.5;}
         elsif ($decimal > 0.25)
            {$currentBondAngle = int($values[$#values])+0.5;}
         else
            {$currentBondAngle = int($values[$#values]);}
         $bondAnglesExt[$atomNumber][$angle] = $currentBondAngle;

         # Classify the bond angle type as a unique element species combination
         #   where the 1st and 3rd tags can be in any order (and so we will
         #   sort them to make identifying unique sets easy). Physically, the
         #   order is unimportant though. However, the second tag must remain
         #   as the second because that signifies the vertex atom of the angle.
         if ($values[0] le $values[2])
            {$angleTag="@values[0..2] $currentBondAngle";}
         else
            {$angleTag="$values[2] $values[1] $values[0] $currentBondAngle";}

         # Search the list of known tags looking for the one we just made.
         #   Note that when $uniqueAngleTags is empty, then $numUniqueAngleTags
         #   is equal to 0 and so the inner part of the loop will not be run.
         $found = 0;
         foreach $tag (1..$numUniqueAngleTags)
         {
            if ($angleTag eq $uniqueAngleTags[$tag])
               {$found = $tag;}
         }

         # If we don't find the tag in the list of known tags, then increase
         #   the number of known tags and store the tag.
         if ($found == 0)
         {
            $numUniqueAngleTags++;
            $found = $numUniqueAngleTags;
            $uniqueAngleTags[$found] = $angleTag;
         }

         # Identify the current atom.
         $atom = $values[4];

         # For the current angle, store the associated unique tag number.
         $angleTagID[$atom][$angle] = $found;

         # Store the two atoms associated with this bond angle vertex.
         $angleBonded[$atom][$angle][1] = $values[3];
         $angleBonded[$atom][$angle][2] = $values[5];
      }
   }

   close (BONDS);
}


sub rotateArbAxis
{
   # Define passed parameters.
   my $point_ref    = $_[0];
   my $unitVect_ref = $_[1];
   my $rotAngle     = $_[2];

   # Define local variables.
   my @orig;

   # Establish the origin of the rotation axis.
   @orig = (0,0,0,0);

   # Create a rotation matrix for the requested rotation.
   &defineRotMatrix($rotAngle,$unitVect_ref);

   # Apply the rotation matrix.
   &rotateOnePoint($point_ref,\@orig);
}


sub defineRotMatrix
{
   # Define passed variables.
   my $rotAngle = $_[0];   # Angle of rotation.
   my $rot_ref  = $_[1];   # Axis of rotation.

   # Define local variables.
   my $c; # cos angle
   my $s; # sin angle
   my $t; # 1-cos angle

   # Create the rotation matrix shown
   #   here:  http://mathworld.wolfram.com/RotationFormula.html
   #   here:  http://www.gamedev.net/reference/articles/article1199.asp
   #      There is an error in the final matrix of the above reference, but it
   #      has the complete derivation without errors so I kept it here.  The
   #      error is the row 1 col 3 term.  It should be txz-sy NOT txz-xy.
   #   here:  http://www.fastgraph.com/makegames/3drotation/
   #   here:  http://en.wikipedia.org/wiki/Rotation_matrix
   $c = cos($rotAngle);
   $s = sin($rotAngle);
   $t = 1.0-cos($rotAngle);
   $rotMatrix = new Math::MatrixReal(3,3);
   $rotMatrix->assign(1,1,$t*$rot_ref->[1]**2 + $c);
   $rotMatrix->assign(1,2,$t*$rot_ref->[1]*$rot_ref->[2] + $s*$rot_ref->[3]);
   $rotMatrix->assign(1,3,$t*$rot_ref->[1]*$rot_ref->[3] - $s*$rot_ref->[2]);
   $rotMatrix->assign(2,1,$t*$rot_ref->[1]*$rot_ref->[2] - $s*$rot_ref->[3]);
   $rotMatrix->assign(2,2,$t*$rot_ref->[2]**2 + $c);
   $rotMatrix->assign(2,3,$t*$rot_ref->[2]*$rot_ref->[3] + $s*$rot_ref->[1]);
   $rotMatrix->assign(3,1,$t*$rot_ref->[1]*$rot_ref->[3] + $s*$rot_ref->[2]);
   $rotMatrix->assign(3,2,$t*$rot_ref->[2]*$rot_ref->[3] - $s*$rot_ref->[1]);
   $rotMatrix->assign(3,3,$t*$rot_ref->[3]**2 + $c);
}


sub rotateOnePoint
{
   # Define passed parameters.
   my $point_ref = $_[0];
   my $orig_ref  = $_[1];

   # Define local variables.
   my $pointVector;
   my $rotPointVector;

   # Create a vector from the point position.
   $pointVector = Math::MatrixReal->new_from_rows(
         [[$point_ref->[1] - $orig_ref->[1],
           $point_ref->[2] - $orig_ref->[2],
           $point_ref->[3] - $orig_ref->[3]]]);

   # Apply the rotation matrix.
   $rotPointVector = $pointVector->multiply($rotMatrix);

   # Save into the given point position.
   $point_ref->[1] = $rotPointVector->element(1,1) + $orig_ref->[1];
   $point_ref->[2] = $rotPointVector->element(1,2) + $orig_ref->[2];
   $point_ref->[3] = $rotPointVector->element(1,3) + $orig_ref->[3];
}

sub rotateAllAtoms
{
   # Define passed parameters.
   my $orig_ref = $_[0];   # Origin of the rotation axis.
   my $rot_ref  = $_[1];   # Axis of rotation.
   my $rotAngle = $_[2];   # Angle of rotation.

   # Define local variables.
   my $atom;
   my $axis;
   my @directXYZCopy;

   # Establish the rotation matrix.
   &defineRotMatrix($rotAngle,$rot_ref);

   # Make a copy of the directXYZ locations in case any atoms need to be
   #   re-rotated because they went outside the simulation box.
   foreach $atom (1..$numAtoms)
   {
      foreach $axis (1..3)
         {$directXYZCopy[$atom][$axis] = $directXYZ[$atom][$axis];}
   }

   # For each atom in the model rotate it in the requested plane by the
   #   requested number of degrees.
   foreach $atom (1..$numAtoms)
   {
      # Apply the rotation matrix to the current atom.
      &rotateOnePoint(\@{$directXYZ[$atom]},$orig_ref);

      # Propogate to other representations.
      &getDirectABC($atom);
      &getFractABC($atom);
   }

   # It is almost certain that some atoms will be outside the simulation box
   #   after the rotation.  These need to be properly addressed since we cannot
   #   simply shift them to the other side of the box and claim the application
   #   of periodic boundary conditions.  We must instead do the following:
   # (1) Restore the atom to its original position.
   # (2) Shift the atom by one lattice length in the direction of the offense.
   # (3) Reapply the rotation.

   # (Essentially the trick is this:  Save a copy of the original XYZ
   #   positions. (This has already been done.) Apply the requested rotation to
   #   all coordinate systems (directABC, fractABC, directXYZ, but not the
   #   copy of the original). (This has also already been done.)  Then restore
   #   ONLY the directXYZ to its original value (before rotation).  Then check
   #   to see which atoms are outside the box ***according to the now rotated
   #   fractABC***.  If an atom is outside the box, then use the restored
   #   directXYZ to get the original fractABC.  This original fractABC is then
   #   shifted so that when the rotation is reapplied it will remain in the
   #   box.  (All of that bit starting with the check is accomplished inside
   #   of checkBoundingBox.) Once the original fractABC is shifted, the
   #   directXYZ and directABC are derived from it.)  Ugly.

   # Restore the original positions of the directXYZ positions for all atoms.
   foreach $atom (1..$numAtoms)
   {
      foreach $axis (1..3)
         {$directXYZ[$atom][$axis] = $directXYZCopy[$atom][$axis];}
   }

   # Make sure that all atoms are inside the simulation box and shifted
   #   properly. This will address (1) and (2) above.
   &checkBoundingBox(1);

   # Repeat the rotation.  (3) from above.
   foreach $atom (1..$numAtoms)
   {
      # Apply the rotation matrix to the current atom.
      &rotateOnePoint(\@{$directXYZ[$atom]},$orig_ref);

      # Propogate to other representations.
      &getDirectABC($atom);
      &getFractABC($atom);
   }

   # At this point there should be no atoms outside the simulation box at all.
}

sub prepSurface
{
   # Define passed parameters.
   my $hkl_ref = $_[0];

   # Define the local variables.
   my $i;
   my $j;
   my $k;
   my $l;
   my $m;
   my $found;
   my $product;
   my $outside;
   my $same;
   my $onFace;
   my @tempFractABC;
   my $rotAngle;
   my @rotAxis;
   my $normalizer;
   my @diagonal;
   my $diagonalMag;
   my $atom;
   my $newAtom;
   my $face;
   my @faceNormal;
   my @faceCenter;
   my @face2Atom;
   my @face2CellCenter;
   my $hklIndex;
   my $leastCommonMult;
   my @maxRep;
   my @currentRep;
   my $cell;
   my $numCells;
   my @cellDims;
   my @cellDist;
   my @lattPointTemp;
   my $xyzAxis;            # An axis iterator over xyz axes.
   my $abcAxis;            # An axis iterator over abc axes.
   my $uvwMag;
   my @uvw;                # Real space vector associated with hkl.
   my @uvwNormalLattice;
   my @uvwRealLattice;     # Lattice of a real cell defining @uvw vector.
   my @newRealLattice;     # Actual new lattice for simulation.
   my @newRealLatticeCenter;
   my @maxXYZ;
   my @origin;             # Vector for the origin (0,0,0).
   my @xAxis;              # Vector for the x-axis (1,0,0).
   my @yAxis;              # Vector for the y-axis (0,1,0).
   my @yzProj;             # Vector for the yz projection of a real lat. vect.
   my $newNumAtoms;        # Local version to replace old one.
   my @newMag;             # Local version to replace old one.
   my @newAtomElementName; # Local version to replace old one.
   my @newAtomElementID;   # Local version to replace old one.
   my @newAtomSpeciesID;   # Local version to replace old one.
   my @newDirectXYZ;       # Local version to replace old one.
   my @newAtomTag;         # Local version to replace old one.
   my @hklRecipLattice;
   my @tempVector;
   my @sortedIndices;  # Index=current sorted atom; value=where it came from.

   # Pre-compute the diagonal of the current cell. It is used later.
   @diagonal = (0,0,0,0); #Uses 1..3
   foreach $abcAxis (1..3)
   {
      foreach $xyzAxis (1..3)
         {$diagonal[$xyzAxis] += $realLattice[$abcAxis][$xyzAxis];}
   }
   $diagonalMag = sqrt($diagonal[1]**2 + $diagonal[2]**2 + $diagonal[3]**2);

   # We are given a reciprocal space hkl. This defines a plane in real space
   #   that intercepts the primitive real space cell at a1/h a2/k a3/l.
   #   However, we need to find the real space lattice vector uvw that defines
   #   a plane with the same orientation but for which the normal vector is
   #   longer. In fact, the length must be such that it points to a specific
   #   real space lattice site. (I.e. the vector components are integral
   #   multiples of the real space lattice vector components.)
   # Therefore, we first compute the least common multiple of the hkl values
   #   with the caveat that an hkl value of zero will not be used to compute
   #   the least common multiple.
   $leastCommonMult = 0;
   $product = 1;
   foreach $abcAxis (1..3)
   {
      if ($hkl_ref->[$abcAxis] != 0)
         {$product *= $hkl_ref->[$abcAxis];}
   }
   foreach $i (1..$product)
   {
      # Assume that this will be *the* $i that is the least common multiple.
      $found = 1;

      # Check each hkl to see if it is equal to zero. Only then do we check if
      #   this $i can be divided by the hkl value.
      foreach $hklIndex (1..3)
      {
         if ($hkl_ref->[$hklIndex] != 0)
         {
            if ($i % $hkl_ref->[$hklIndex] != 0)
               {$found = 0; last;}
         }
      }

      # If all of the hkl values evenly divide into $i, then $i is the least
      #   common multiple we are looking for. (Obviously ignoring the hkl=0
      #   values.)
      if ($found == 1)
         {$leastCommonMult = $i; last;}
   }

   # Create two similar expressions. The first is called the uvwNormalLattice.
   #   The uvwNormalLattice takes the form of three abc vectors, each expressed
   #   in Cartesian xyz coordinates. The sum of the three vectors is the vector
   #   that is normal to the hkl defined plane. Note that the uvwNormalLattice
   #   makes use of integer multiples of the original realLattice vectors. The
   #   uvwNormalLattice is created by multiplying by 1/hkl of the appropriate
   #   hkl vector *and* the least common multiple of the hkl values that was
   #   determined above. The exception to that process is the case where one
   #   of the hkl values is zero. In that case the associated uvwNormalLattice
   #   vector is set to zero as well. Essentially, the resultant uvw vector
   #   (R = ua_1 + va_2 + wa_3) will be the normal to the plane defined by the
   #   hkl vector: K = hb_1 + kb_2 + lb_3 where b_1,2,3 are reciprocal lattice
   #   primitive vectors (expressed in xyz coordinates).
   # Second, create a uvwRealLattice (abc vectors expressed in xyz coordinates)
   #   that is the exact same as the uvwNormalLattice except for the cases
   #   where the hkl value is zero. In this case the uvwRealLattice will have
   #   a value equal to the original realLattice vector.
   # So, in effect the uvwNormalLattice defines a normal vector while the
   #   uvwRealLattice defines an actual lattice.
   foreach $abcAxis (1..3)
   {
      foreach $xyzAxis (1..3)
      {
         if ($hkl_ref->[$abcAxis] == 0)
         {
            $uvwNormalLattice[$abcAxis][$xyzAxis] = 0.0;
            $uvwRealLattice[$abcAxis][$xyzAxis] =
                  $realLattice[$abcAxis][$xyzAxis];
         }
         else
         {
            $uvwNormalLattice[$abcAxis][$xyzAxis] = $leastCommonMult * 
                  $realLattice[$abcAxis][$xyzAxis] / $hkl_ref->[$abcAxis];
            $uvwRealLattice[$abcAxis][$xyzAxis] = $leastCommonMult * 
                  $realLattice[$abcAxis][$xyzAxis] / $hkl_ref->[$abcAxis];
         }
      }
   }

   # The accumulated sum of the uvwNormal vectors (in xyz coordinates) is
   #   the normal vector for the plane defined by the hkl reciprocal lattice
   #   vector.
   @uvw = (0,0,0,0);
   foreach $xyzAxis (1..3)
   {
      $uvw[$xyzAxis] = 0.0;
      foreach $abcAxis (1..3)
         {$uvw[$xyzAxis] += $uvwNormalLattice[$abcAxis][$xyzAxis];}
   }
   $uvwMag = sqrt($uvw[1]**2 + $uvw[2]**2 + $uvw[3]**2);

   # However, the new a and b axes are a bit trickier. The essential question
   #   is the following: Given a lattice point defined by the uvwRealLattice,
   #   what is the minimal distance that one needs to "travel" in the plane
   #   perpendicular to the new c-axis before encounting that lattice point
   #   again? This direction and distance will be the a-axis for our new real
   #   lattice. This question must be answered a second time in some different
   #   direction to determine the b-axis for the new real lattice.
   # In one case the question is very easy to answer. That is the case where
   #   the one of the hkl values is 0 because the lattice in that direction
   #   will be unchanged (except that it may be rotated). As an example,
   #   consider (1 1 0) as the hkl values. The original cell will be rotated
   #   such that the vector from 0,0,0 to 1,1,0 becomes the new (0 0 |c|) axis.
   #   Thus, if we make the old z-axis become the new a-axis, then the a-axis
   #   will have a value of (|old z mag| 0 0). After that, we need to find
   #   the magnitude of the b-axis in a similar way as will be described below.
   # If we don't have the easy case, then we need to look for two other lattice
   #   points that are not all co-linear with the uvw vector point but which
   #   are positioned on the plane defined by the uvw vector. (In the easy
   #   case we just need one and can ignore searching in the direction of the
   #   easy lattice vector.)
   # We will find the lattice points by executing a spiral search out from the
   #   origin in steps of uvwRealLattice. This follows the same algorithm as
   #   the makeLattice subroutine in lattice.f90 of the olcao code. Note that
   #   uvwRealMatrix has dimensions [abc][xyz] while in the f90 code the
   #   matrix has dimensions [xyz][abc]. To begin the spiral search we first
   #   need to define the limits of the spiral.
   foreach $i (1..3)
   {
      if ($hkl_ref->[$i] == 0)
         {$maxRep[$i] = 5;}
      else
         {$maxRep[$i] = 100;}  # Assuming for now that this will be big enough.
   }

   # 
   $numCells = 0;
   foreach $i (0..$maxRep[1]+1)
   {
      $currentRep[1] = $i - int((1+$maxRep[1])/2);
      foreach $j (0..$maxRep[2]+1)
      {
         $currentRep[2] = $j - int((1+$maxRep[2])/2);
         foreach $k (0..$maxRep[3]+1)
         {
            $currentRep[3] = $k - int((1+$maxRep[3])/2);

            # Create a lattice point for the current set of lattice repititions
            #   by summing the xyz contributions of the uvwRealLattice times
            #   the appropriate repitiiion number.
            @lattPointTemp = (0,0,0,0);
            foreach $l (1..3) # $abcAxis
            {
               foreach $m (1..3) # xyzAxis
                  {$lattPointTemp[$m] += ($uvwRealLattice[$l][$m] * 
                                          $currentRep[$l]);}
            }

            # Store the result for later sorting and analysis.
            $numCells += 1;
            $cellDist[$numCells] = sqrt($lattPointTemp[1]**2 +
                                        $lattPointTemp[2]**2 +
                                        $lattPointTemp[3]**2);
            $cellDims[1][$numCells] = $lattPointTemp[1];
            $cellDims[2][$numCells] = $lattPointTemp[2];
            $cellDims[3][$numCells] = $lattPointTemp[3];
         }
      }
   }

   # Sort the cells according to the magnitude of their distance from the
   #   origin. (First shift the cellDist to make sorting correct.)
   # After applying the sort, each index of the @sortedIndices array will
   #   contain an index value that references the original cellDist array. So,
   #   for example, if the cellDist array holds (37.3, 83.9, 11.5, 403.9, 2.1)
   #   then the sortedIndices array will hold (4, 2, 0, 1, 3). You should
   #   interpret that as saying "what was in index 4 is now in index 0", "what
   #   was in index 2 is now in index 1", "what was in index 0 is now in index
   #   2", etc.
   # Note importantly that the cellDist array will not actually be sorted here.
   #   We only have the indices that each cell distance should go in to so as
   #   to sort it (or anything that had the same original order and that needs
   #   to be sorted along with it (e.g., the cellDims)).
   shift(@cellDist);
   @sortedIndices = sort {$cellDist[$a] <=> $cellDist[$b]}
         0..$#cellDist;

   # Adjust the indices to match the indices that start from 1.
   unshift (@cellDist,"");
   unshift (@sortedIndices,"");
   foreach $cell (1..$numCells)
      {$sortedIndices[$cell]++;}

   # Apply this set of indices to the ordered list of cellDims.
   &applySort(\@{$cellDims[1]},\@sortedIndices);
   &applySort(\@{$cellDims[2]},\@sortedIndices);
   &applySort(\@{$cellDims[3]},\@sortedIndices);
   &applySort(\@cellDist,\@sortedIndices);

   # Now we can finally establish the newRealLattice and newMag. We do this by
   #   checking each of the lattice points in cellDims to determine if the
   #   vector defined by that lattice point is perpendicular to the uvw vector.
   #   (Obviously we skip the trivial (0,0,0) case.)  When we find the first
   #   cellDims vector that is perpendicular to uvw we can identify the
   #   distance to this lattice point as newMag[2] and the vector as the new
   #   b-axis vector for newRealLattice. Any future lattice point that is found
   #   must be perpendicular to uvw *and* non-colinear with the first
   #   discovered lattice point. The second vector that is found that satisfies
   #   the criteria will be the new c-axis and the distance to it will be the
   #   newMag[3] value.
   @newMag = (0,0,0,0);
   foreach $cell (2..$numCells) # Start at 2 to skip (0,0,0).
   {
      @tempVector = (0,$cellDims[1][$cell],$cellDims[2][$cell],
                       $cellDims[3][$cell]);
      if (abs(&dotProduct(\@tempVector,\@uvw)) < $epsilon)
      {
         if ($newMag[2] == 0)
         {
            $newMag[2] = $cellDist[$cell];
            $newRealLattice[2][1] = $cellDims[1][$cell];
            $newRealLattice[2][2] = $cellDims[2][$cell];
            $newRealLattice[2][3] = $cellDims[3][$cell];
         }
         elsif (abs(&crossProductMag(\@tempVector,\@{$newRealLattice[2]})) >
                    $epsilon)
         {
            $newMag[3] = $cellDist[$cell];
            $newRealLattice[3][1] = $cellDims[1][$cell];
            $newRealLattice[3][2] = $cellDims[2][$cell];
            $newRealLattice[3][3] = $cellDims[3][$cell];
            last;
         }
      }
   }

   # We reserve the a-axis to be the axis perpendicular to the plane of the
   #   surface. Later, the cell will be rotated so that this axis is co-linear
   #   with the Cartesean x-axis.
   $newMag[1] = $uvwMag;
   $newRealLattice[1][1] = $uvw[1];
   $newRealLattice[1][2] = $uvw[2];
   $newRealLattice[1][3] = $uvw[3];


   # At this point we can put the newRealLattice parameters into realLattice,
   #   and similarly for the lattice vector magnitudes. We do this here becaues
   #   we need to get some directABC coordinates of atoms to eleiminate
   #   duplicates.  Many part of the remainder of this subroutine do not
   #   refer yet to the realLattice, but still refer to the newRealLattice.
   #   That can be changed in the future.
   foreach $abcAxis (1..3)
   {
      foreach $xyzAxis (1..3)
      {
         $realLattice[$abcAxis][$xyzAxis] =
               $newRealLattice[$abcAxis][$xyzAxis];
      }
      $mag[$abcAxis] = $newMag[$abcAxis];
   }

   # Obtain the inverse lattice, the reciprocal lattice vectors, and the
   #   angles between the lattice vectors. Note that this is only temporary
   #   and that these will need to be re-obtained later after the lattice
   #   goes through a series of rotations later. We do this now though to make
   #   the process of identifying duplicate atoms in the new cell easier.
   &makeLatticeInv(\@realLattice,\@realLatticeInv,0);
   &makeLatticeInv(\@realLattice,\@recipLattice,1);
   &abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);

   # At this point our new real lattice is defined with its periodic boundary
   #   conditions matching up with lattice points of the original real lattice.
   #   What we want to do now is fill this new lattice with atoms from a super
   #   lattice (made from the original real lattice). The approach for the time
   #   being will be to iterate through each of the lattices defined by the
   #   cellDims points and iterate through each atom in each of those cells
   #   checking to see if the selected atom is inside the new real lattice. A
   #   future (better) approach would figure out a way to reduce the number of
   #   cells to have to check. Note that this is tricky because one could
   #   conceivably have an oddly shaped cell that passes through the new cell
   #   but which doesn't have any lattice points within the new cell. At any
   #   rate...
   # The procedure for doing the check is the following:
   # First: compute the center point of the new real lattice.
   # Second: compute the center point of each face of the new real lattice.
   # Third: make a set of normal vectors for each face of the new real lattice
   #   with the condition that the normals point outward.
   # Fourth: iterate through each atom of each replicated cell and:
   # Fourth': compute a difference vector that points from the center of each
   #   face to the selected atomic point.
   # Fourth'': compute the dot product of the difference vector and the normal
   #   vector of each face. If any of the dot products has a negative sign,
   #   then the atom is outside of the cell.

   # Compute position of the center point of the new real lattice.
   foreach $xyzAxis (1..3)
   {
      $newRealLatticeCenter[$xyzAxis] = ($newRealLattice[1][$xyzAxis] + 
                                         $newRealLattice[2][$xyzAxis] +
                                         $newRealLattice[3][$xyzAxis])/2.0
   }


   # Compute the position of the center of each face of the new real lattice.
   foreach $xyzAxis (1..3)
   {
      $faceCenter[1][$xyzAxis] = ($newRealLattice[1][$xyzAxis] +
                                  $newRealLattice[2][$xyzAxis]) / 2.0;
      $faceCenter[2][$xyzAxis] = $faceCenter[1][$xyzAxis] +
                                 $newRealLattice[3][$xyzAxis];
      $faceCenter[3][$xyzAxis] = ($newRealLattice[1][$xyzAxis] +
                                  $newRealLattice[3][$xyzAxis]) / 2.0;
      $faceCenter[4][$xyzAxis] = $faceCenter[3][$xyzAxis] +
                                 $newRealLattice[2][$xyzAxis];
      $faceCenter[5][$xyzAxis] = ($newRealLattice[2][$xyzAxis] +
                                  $newRealLattice[3][$xyzAxis]) / 2.0;
      $faceCenter[6][$xyzAxis] = $faceCenter[5][$xyzAxis] +
                                 $newRealLattice[1][$xyzAxis];
   }

   # Compute normal vectors for each face.
   @origin = (0,0,0,0); # Uses indices 1..3.
   @{$faceNormal[1]} = &getPlaneNormal(\@origin,\@{$newRealLattice[1]},
                                    \@{$newRealLattice[2]});
   @{$faceNormal[3]} = &getPlaneNormal(\@origin,\@{$newRealLattice[1]},
                                    \@{$newRealLattice[3]});
   @{$faceNormal[5]} = &getPlaneNormal(\@origin,\@{$newRealLattice[2]},
                                    \@{$newRealLattice[3]});
   foreach $xyzAxis (1..3)
   {
      $faceNormal[2][$xyzAxis] = -1.0*$faceNormal[1][$xyzAxis];
      $faceNormal[4][$xyzAxis] = -1.0*$faceNormal[3][$xyzAxis];
      $faceNormal[6][$xyzAxis] = -1.0*$faceNormal[5][$xyzAxis];
   }

   # Make sure that the normal vectors point out. Because d = |a||b|*cos(theta)
   #   we have that the sign of d is determined by whether or not the two
   #   vectors face the same direction. Therefore, we create a vector that
   #   points from the face center to the cell center and compare it to the
   #   face normals computed above. If the sign is positive, then the normal
   #   should be reversed in sign (because it too is pointing "in" like the
   #   vector from the face to the center). If the sign is negative, then the
   #   normal can stay as it is.
   foreach $face (1..6)
   {
      # Compute a vector pointing from the face center to the cell center.
      foreach $xyzAxis (1..3)
         {$face2CellCenter[$xyzAxis] = $newRealLatticeCenter[$xyzAxis] - 
                                       $faceCenter[$face][$xyzAxis];}

      # Compute the dot product, compare, and flip if needed.
      $product = &dotProduct(\@{$faceNormal[$face]},\@face2CellCenter);
      if ($product > 0.0)
      {
         foreach $xyzAxis (1..3)
            {$faceNormal[$face][$xyzAxis] = -1.0*$faceNormal[$face][$xyzAxis];}
      }

      # Make the face normal a unit vector.
      $normalizer = sqrt($faceNormal[$face][1]**2 + $faceNormal[$face][2]**2 +
                         $faceNormal[$face][3]**2);
      foreach $xyzAxis (1..3)
         {$faceNormal[$face][$xyzAxis] /= $normalizer;}
   }


   # Now we have the tools to check if each atom in a replicated original cell
   #   is (or is not) inside of the new cell for the surface model. Note that
   #   we can avoid checking the atoms of a given cell if the magnitude of the
   #   distance from the origin to cell lattice point is greater than the
   #   magnitude of the diagonal of the original cell. Additionally, we will
   #   pre-compute the maximum xyz values to help eliminate duplicate atoms.
   $newNumAtoms = 0;
   foreach $cell (1..$numCells)
   {
      if ($cellDist[$cell] > $diagonalMag)
         {next;}

      foreach $atom (1..$numAtoms)
      {
         # Assume that this atom will be included.
         $newNumAtoms++;

         # Compute the position of this atom of this replicated cell.
         foreach $xyzAxis (1..3)
            {$newDirectXYZ[$newNumAtoms][$xyzAxis] =
                  $directXYZ[$atom][$xyzAxis] + $cellDims[$xyzAxis][$cell];}

         # Assume that the atom is inside the new cell and the check. If it is
         #   actually outside, then we abort on this atom and go to the next.
         $outside = 0;
         foreach $face (1..6)
         {
            # Compute a vector from each new cell face to the position of this
            #   replicated atom. Then apply the same dot product test as above
            #   to determine if the atom is inside or outside of the cell.
            foreach $xyzAxis (1..3)
               {$face2Atom[$xyzAxis] = $newDirectXYZ[$newNumAtoms][$xyzAxis] - 
                                       $faceCenter[$face][$xyzAxis];}

            $product = &dotProduct(\@{$faceNormal[$face]},\@face2Atom);
            if ($product > $epsilon)
               {$outside = 1; last;}
         }
         if ($outside == 1)
            {$newNumAtoms--; next;}

         # Now we need to prevent the inclusion of atoms that copy another atom
         #   by being on a face. That is, consider an atom at 0,0,0 and another
         #   atom at 0,0,c. These are the same atom and one of them must be
         #   eliminated. We will achieve this goal by setting the appropriate
         #   x, y, or z coordinate equal to zero for any atomic coordinate that
         #   is equal to the maximum allowed xyz coordinate. Then, the next bit
         #   of code will remove any duplicate.
         $onFace = 0;
         @tempFractABC = &directXYZ2fractABC(\@{$newDirectXYZ[$newNumAtoms]},
               \@realLatticeInv);
         unshift(@tempFractABC,"");
         foreach $abcAxis (1..3)
         {
            if (abs(abs($tempFractABC[$abcAxis]) - 1.0) < $epsilon)
               {$onFace = 1; last;}
         }
         if ($onFace == 1)
            {$newNumAtoms--; next;}


         # Check if this atom is an exact copy of an already included newAtom.
         $same = 0;
         foreach $newAtom (1..$newNumAtoms-1)
         {
            if ((abs($newDirectXYZ[$newAtom][1] -
                     $newDirectXYZ[$newNumAtoms][1]) < $epsilon) &&
                (abs($newDirectXYZ[$newAtom][2] -
                     $newDirectXYZ[$newNumAtoms][2]) < $epsilon) && 
                (abs($newDirectXYZ[$newAtom][3] -
                     $newDirectXYZ[$newNumAtoms][3]) < $epsilon))
               {$same = 1; last;}
         }
         if ($same == 1)
            {$newNumAtoms--; next;}

         # If we make it this far, then we will keep the atom.
         $newAtomElementName[$newNumAtoms] = $atomElementName[$atom];
         $newAtomElementID[$newNumAtoms]   = $atomElementID[$atom];
         $newAtomSpeciesID[$newNumAtoms]   = $atomSpeciesID[$atom];
         $newAtomTag[$newNumAtoms]         = $atomTag[$atom];
      }
   }

   # Update the system information.
   undef (@directXYZ);
   undef (@atomElementName);
   undef (@atomElementID);
   undef (@atomSpeciesID);
   undef (@atomTag);
   $numAtoms        = $newNumAtoms;
   @directXYZ       = @newDirectXYZ;
   @atomElementName = @newAtomElementName;
   @atomElementID   = @newAtomElementID;
   @atomSpeciesID   = @newAtomSpeciesID;
   @atomTag         = @newAtomTag;
   undef (@newDirectXYZ);
   undef (@newAtomElementName);
   undef (@newAtomElementID);
   undef (@newAtomSpeciesID);
   undef (@newAtomTag);
   undef (@cellDims);
   undef (@cellDist);


   # Now we are going to rotate the lattice so that the a-axis is perpendicular
   #   to the plane of the surface and co-linear with the Cartesean x-axis.

   # Compute a vector that is perpendicular to the plane formed by the new real
   #   space vector (@uvw) and the cartesian x-axis. This will be the axis of
   #   rotation.
   @xAxis  = (0,1,0,0); # Uses indices 1..3.
   @rotAxis = &getPlaneNormal(\@origin,\@xAxis,\@uvw);

   # Make the rotation vector a unit vector.
   $normalizer = sqrt($rotAxis[1]**2 + $rotAxis[2]**2 + $rotAxis[3]**2);
   foreach $xyzAxis (1..3)
      {$rotAxis[$xyzAxis] /= $normalizer;}

   # Compute the angle of rotation between the real space vector and the x-axis.
   $rotAngle = &getVectorAngle(\@uvw,\@xAxis);

   # In the event that the uvw vector is pointing into any positive y quadrant
   #   then we need to rotate *back* to the x-axis. In that case we invert the
   #   rotation axis.
   if ($uvw[2] > 0.0)
   {
      foreach $xyzAxis (1..3)
         {$rotAxis[$xyzAxis] = -1.0*$rotAxis[$xyzAxis];}
   }

   # Establish the rotation matrix.
   &defineRotMatrix($rotAngle,\@rotAxis);

   # Apply the rotation matrix to the current lattice vectors.
   &rotateOnePoint(\@{$realLattice[1]},\@origin);
   &rotateOnePoint(\@{$realLattice[2]},\@origin);
   &rotateOnePoint(\@{$realLattice[3]},\@origin);
   foreach $abcAxis (1..3)
   {
      foreach $xyzAxis (1..3)
      {
         if (abs($realLattice[$abcAxis][$xyzAxis]) < $epsilon)
            {$realLattice[$abcAxis][$xyzAxis] = 0.0;}
      }
   }

   # For each atom in the model rotate it in the requested plane by the
   #   requested number of degrees (associated with the defined rotation
   #   matrix).
   foreach $atom (1..$numAtoms)
      {&rotateOnePoint(\@{$directXYZ[$atom]},\@origin);}

   # Now we need to rotate the b-axis so that it lies in the Cartesean x-y
   #   plane with zero component in the z-axis.

   # The rotation axis will be the current a-axis (which is also the current
   #   x-axis).
   @rotAxis  = (0,1,0,0); # Uses indices 1..3.

   # Compute the angle of rotation between the current y-z component of the
   #   b-axis and the Cartesean x-y plane. This is a bit tricky. Just recall
   #   by this point we have the lattice a-axis and x-axis aligned, but that
   #   the b-axis is pointing in some xyz direction. We want to rotate it into
   #   the x-y plane so that it has zero z-component. Thus, the angle that it
   #   needs to rotate through is the same as the angle between the y-axis
   #   and the vector that is a projection of the b-axis onto the y-z plane.
   # If we just tried to compute the angle between the b-axis and the y-axis
   #   then we would get a totally bogus angle if the b-axis was already in
   #   the x-y plane but not aligned with the y-axis. We would not need to
   #   rotate the b-axis in that case, but we would have computed a non-zero
   #   angle for the rotation.
   @yzProj = (0,0,$realLattice[2][2],$realLattice[2][3]); # Uses indices 1..3
   @yAxis  = (0,0,1,0); # Uses indices 1..3
   $rotAngle = &getVectorAngle(\@yzProj,\@yAxis);

   # In the event that the yzProj vector is pointing into any positive z
   #   quadrant then we need to rotate *back* to the y-axis. In that case we
   #   invert the rotation axis.
   if ($yzProj[2] > 0.0)
      {@rotAxis = (0,-1,0,0);}

   # Establish the rotation matrix.
   &defineRotMatrix($rotAngle,\@rotAxis);

   # Apply the rotation matrix to the current lattice vectors.
   &rotateOnePoint(\@{$realLattice[1]},\@origin);
   &rotateOnePoint(\@{$realLattice[2]},\@origin);
   &rotateOnePoint(\@{$realLattice[3]},\@origin);
   foreach $abcAxis (1..3)
   {
      foreach $xyzAxis (1..3)
      {
         if (abs($realLattice[$abcAxis][$xyzAxis]) < $epsilon)
            {$realLattice[$abcAxis][$xyzAxis] = 0.0;}
      }
   }

   # For each atom in the model rotate it in the requested plane by the
   #   requested number of degrees (associated with the defined rotation
   #   matrix).
   foreach $atom (1..$numAtoms)
      {&rotateOnePoint(\@{$directXYZ[$atom]},\@origin);}

   # Reobtain the inverse lattice and the reciprocal lattice vectors.
   &makeLatticeInv(\@realLattice,\@realLatticeInv,0);
   &makeLatticeInv(\@realLattice,\@recipLattice,1);
   &abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);

   # Obtain the directABC and fractABC coordinates of the atoms in the model
   #   such that they match the directXYZ.
   foreach $atom (1..$numAtoms)
   {
      &getDirectABC($atom);
      &getFractABC($atom);
   }
}

sub applySpaceGroup
{
   # Define local variables.
   my $axis;
   my $atom;
   my @values;
   my @applySpaceGroupIn;

   # Use the space group ID (name or number) to open the space group operations
   #   file from the space group data base.
   open (SPACE,"<$spaceGroupDB/$spaceGroup") ||
         die "Cannot open $spaceGroup for reading in $spaceGroupDB.\n";

   # Read the first line to get the cell lattice type (primitive or non-
   #   primitive) and the space group name.
   @values = &prepLine(\*SPACE,"",'\s+');
   $latticeType = $values[0];
   $spaceGroupName = $values[1];

   # If there is an alternative definition for this space group it is labeled
   #   in the space group file with an index starting at 'a'.  In such a case,
   #   the file name is extended by "_a" for the complete name of the space
   #   group (i.e. not the number)  (e.g. 227 = Fd3~m_a).
   if (($#values > 1) && ($values[2] =~ /[a-z]/))
      {$spaceGroupName = $spaceGroupName . "_$values[2]";}

   # Read the second line to get the root space group number.
   @values = &prepLine(\*SPACE,"",'\s+');
   $spaceGroupNum = $values[0];
   $spaceGroupSubNum = $values[1];

   # Read remaining lines and prepare input for applySpaceGroup.f90 program.
   @applySpaceGroupIn = <SPACE>;
   unshift (@applySpaceGroupIn,
         "$latticeType\n$spaceGroupNum  $spaceGroupSubNum\n");
   close (SPACE);

   # Add the flag indicating whether to do a full (possibly non-primitive)
   #   or primitive cell.
   if ($doFullCell == 1)
      {push (@applySpaceGroupIn,"1\n");}
   else
      {push (@applySpaceGroupIn,"0\n");}

   # Add lattice parameters to the applySpaceGroupIn array.
   foreach $axis (1..3)
      {push (@applySpaceGroupIn,"@{$realLattice[$axis]}\n");}

   # Add atomic data to the applySpaceGroupIn array.
   push (@applySpaceGroupIn,"$numAtoms\n");
   foreach $atom (1..$numAtoms)
   {
      push (@applySpaceGroupIn,"$atomElementID[$atom] $atomSpeciesID[$atom]" .
            " @{$fractABC[$atom]}\n");
   }

   # Create the hidden input file for the applySpaceGroup program.
   open (SGINPUT,">sginput") || die "Cannot open sginput for writing\n";
   print SGINPUT @applySpaceGroupIn;
   close (SGINPUT);

   # Apply the space group to this structure.
   system ("applySpaceGroup");

   # Open the output file for reading and processing.
   open (SGOUTPUT,"<sgoutput") || die "Cannot open sgoutput for reading\n";

   # Reassign the lattice parameters and propagate to other representations.
   foreach $axis (1..3)
   {
      @{$realLattice[$axis]} = &prepLine(\*SGOUTPUT,"",'\s+');
      unshift (@{$realLattice[$axis]},"");
   }
   &abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);

   # Read the number of atoms in the newly defined system.
   @values = &prepLine(\*SGOUTPUT,"",'\s+');
   $numAtoms = $values[0];

   # Undefine the old data sets based on atom number = array index number.
   #   These will have to be redefined.
   undef (@atomTag);
   undef (@atomElementName);
   undef (@atomElementID);
   undef (@atomSpeciesID);
   undef (@fractABC);

   # Read each of the atoms and their associated element ID and species ID
   #   numbers.  Then rebuild the above undefined arrays.  Note that the
   #   arrays @elementList and @speciesList should be unaffected by the
   #   application of space group operations because the number of different
   #   elements and number of different species for each element should still
   #   be the same.
   foreach $atom (1..$numAtoms)
   {
      @values = &prepLine(\*SGOUTPUT,"",'\s+');
      $atomElementID[$atom]   = $values[0];
      $atomSpeciesID[$atom]   = $values[1];
      $fractABC[$atom][1]     = $values[2];
      $fractABC[$atom][2]     = $values[3];
      $fractABC[$atom][3]     = $values[4];
      $atomElementName[$atom] = $elementList[$values[0]];
      $atomTag[$atom]         = $speciesList[$values[0]][$values[1]];
   }

   # Propagate the coordinate information to other forms.
   foreach $atom (1..$numAtoms)
   {
      &getDirectXYZ($atom);
      &getDirectABC($atom);
   }

   close (SGOUTPUT);
}



sub applySupercell
{
   # Define passed parameters.
   my @scRequest;
   $scRequest[1] = $_[0];
   $scRequest[2] = $_[1];
   $scRequest[3] = $_[2];
   my @scMirror;
   $scMirror[1]  = $_[3];
   $scMirror[2]  = $_[4];
   $scMirror[3]  = $_[5];

   # Define the local variables.
   my $axisXYZ;            # The xyz axes of the cell's vector description.
   my $axisABC;            # The abc axis of the cell.
   my $atom;
   my $cell;
   my $atomAxis;           # The abc axis of each atom.
   my $newNumAtoms;        # Local version to replace old one.
   my @newAtomElementName; # Local version to replace old one.
   my @newAtomElementID;   # Local version to replace old one.
   my @newAtomSpeciesID;   # Local version to replace old one.
   my @newFractABC;        # Local version to replace old one.
   my @newAtomTag;         # Local version to replace old one.

   # Compute the new lattice parameters in both a,b,c,alpha,beta,gamma and
   #   [a,b,c][x,y,z] form.
   foreach $axisABC (1..3)
   {
      # Increase the dimensions of the model along the current axis by a
      #   factor of the number of cells in the direction of the current axis.
      $mag[$axisABC] *= $scRequest[$axisABC];

      # Increase each x,y,z component of each a,b,c vector in the same way.
      foreach $axisXYZ (1..3)
         {$realLattice[$axisABC][$axisXYZ] *= $scRequest[$axisABC];}
   }

   # Apply the supercell requests along each axis to each atom.
   foreach $axisABC (1..3)
   {
      # There is no need to replicate atoms if the supercell request in this
      #   direction is equal to one.
      if ($scRequest[$axisABC] == 1)
         {next;}

      # Loop over each atom and replicate them to fill the new cell.
      $newNumAtoms = 0;
      foreach $atom (1..$numAtoms)
      {
         foreach $cell (1..$scRequest[$axisABC])
         {
            $newNumAtoms++;
            foreach $atomAxis (1..3)
            {
               if ($atomAxis == $axisABC)
               {
                  if (($scMirror[$axisABC] == 0) || (($cell % 2) == 1))
                  {
                     $newFractABC[$newNumAtoms][$atomAxis] = 
                           $fractABC[$atom][$atomAxis] *
                           1.0 / $scRequest[$atomAxis] +
                           ($cell-1.0) / $scRequest[$atomAxis];
                  }
                  else
                  {
                     $newFractABC[$newNumAtoms][$atomAxis] =
                           (1.0-$fractABC[$atom][$atomAxis]) *
                           1.0 / $scRequest[$atomAxis] +
                           ($cell-1.0) / $scRequest[$atomAxis];
                  }
               }
               else
               {
                  $newFractABC[$newNumAtoms][$atomAxis] =
                        $fractABC[$atom][$atomAxis];
               }
            }
            $newAtomElementName[$newNumAtoms] = $atomElementName[$atom];
            $newAtomElementID[$newNumAtoms]   = $atomElementID[$atom];
            $newAtomSpeciesID[$newNumAtoms]   = $atomSpeciesID[$atom];
            $newAtomTag[$newNumAtoms]         = $atomTag[$atom];
         }
      }

      # Update the system information necessary for the next axisABC.
      undef (@fractABC);
      undef (@atomElementName);
      undef (@atomElementID);
      undef (@atomSpeciesID);
      undef (@atomTag);
      $numAtoms        = $newNumAtoms;
      @fractABC        = @newFractABC;
      @atomElementName = @newAtomElementName;
      @atomElementID   = @newAtomElementID;
      @atomSpeciesID   = @newAtomSpeciesID;
      @atomTag         = @newAtomTag;
      undef (@newFractABC);
      undef (@newAtomElementName);
      undef (@newAtomElementID);
      undef (@newAtomSpeciesID);
      undef (@newAtomTag);
   }

   # Reobtain the inverse lattice and the reciprocal lattice vectors.
   &makeLatticeInv(\@realLattice,\@realLatticeInv,0);
   &makeLatticeInv(\@realLattice,\@recipLattice,1);
   &abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);
   &abcAlphaBetaGamma(\@recipLattice,\@magRecip,\@angleRecip,
                      \@angleDegRecip,0);

   # Obtain the directABC and directXYZ coordinates of the atoms in the model
   #   such that they match the fractionalABC.
   foreach $atom (1..$numAtoms)
   {
      &getDirectXYZ($atom);
      &getDirectABC($atom);
   }
}


sub createSpeciesData
{

   # Define passed parameters.
   my $useFileSpecies = $_[0];

   # Declare local variables.
   my $atom;
   my $element;
   my $species;
   my $currElementID;
   my $found;

   # If we ask to not use the species designated in the file then we must
   #   set the number of species of each element to 1, set all the species IDs
   #   to 1, set the speciesList for each species to the element name catted
   #   with a 1, and then return.
   if ($useFileSpecies == 0)
   {
      foreach $element (1..$numElements)
         {$numSpecies[$element] = 1;}
      foreach $atom (1..$numAtoms)
      {
         $speciesList[$atomElementID[$atom]][1] = $atomElementName[$atom] ."1";
         $atomSpeciesID[$atom] = 1;
      }

      return;
   }

   # Proceed only if we will use the species designated in the input file and
   #   stored in the @atomTag data set.

   # Initialize the count of the number of unique species for each element.
   foreach $element (1..$numElements)
      {$numSpecies[$element] = 0;}

   foreach $atom (1..$numAtoms)
   {
      # Define the current element ID.
      $currElementID = $atomElementID[$atom];

      # Determine if the tag associated with the element for this atom
      #   already exists.
      $found = 0;
      foreach $species (1..$numSpecies[$currElementID])
      {
         if ($atomTag[$atom] eq "$speciesList[$currElementID][$species]")
         {
            $found = $species;
            last;
         }
      }

      if ($found == 0)
      {
         $numSpecies[$currElementID]++;
         $speciesList[$currElementID][$numSpecies[$currElementID]] =
               $atomTag[$atom];
         $atomSpeciesID[$atom] = $numSpecies[$currElementID];
      }
      else
         {$atomSpeciesID[$atom] = $found;}
   }
}


sub createElementList
{
   # Declare local variables.
   my $atom;
   my $element;
   my $found;
   my @values;

   # Initialize the element count and element list.
   $numElements = 0;
   $elementList[0] = "";

   foreach $atom (1..$numAtoms)
   {
      # Define the element name for this atom.
      @values = &prepLine("","$atomTag[$atom]",'[0-9]+');
      $atomElementName[$atom] = $values[0];

      # Find the element name of this atom in the list of unique elements.
      $found = 0;
      foreach $element (1..$numElements)
      {
         if ($atomElementName[$atom] eq "$elementList[$element]")
         {
            $found = $element;
            last;
         }
      }

      # Add a new element to the list of unique elements if it was not found.
      if ($found == 0)
      {
         $numElements++;
         $elementList[$numElements] = "$atomElementName[$atom]";
         $atomElementID[$atom] = $numElements;
      }
      else
         {$atomElementID[$atom] = $found;}
   }
}


sub computeImplicitInfo
{
   # Define local variables
   my $atom;
   my $element;
   my $basis;
   my $orbitalType;
   my $currentZ;
   my $tempValue;
   my @tempArray;

   $minPotAlpha = 1000.0;
   $maxNumPotAlphas = 0;
   $maxNumAtomAlphas = 0;
   $maxNumValeStates = 0;
   foreach $atom (1..$numAtoms)
   {
      $currentZ = $atomicZ[$atom];

      $tempValue = ElementData::getMinTermPot($currentZ);
      if ($tempValue < $minPotAlpha)
         {$minPotAlpha = $tempValue;}

      $tempValue = ElementData::getNumTermsPot($currentZ);
      if ($tempValue > $maxNumPotAlphas)
         {$maxNumPotAlphas = $tempValue;}

      @tempArray = ElementData::getNumTermsWF($currentZ);
      if ($tempArray[0] > $maxNumAtomAlphas)
         {$maxNumAtomAlphas = $tempArray[0];}

      $tempValue = 0;
      foreach $basis (1..3) # MB+FB+EB
      {
         @tempArray = ElementData::getValeOrbitals($basis,$currentZ);
         $tempValue += 1*$tempArray[0]; #s
         $tempValue += 3*$tempArray[1]; #p
         $tempValue += 5*$tempArray[2]; #d
         $tempValue += 7*$tempArray[3]; #f
      }
      if ($tempValue > $maxNumValeStates)
         {$maxNumValeStates = $tempValue;}
   }


}


sub computeCrystalParameters
{
   # Declare local variables.
   my $axis;

   # Find the maximum and minimum values of all atoms for each x,y,z axis.
   &getMinMaxXYZ;

   # We will assume an orthorhomic box that contains all the atoms
   #   with a buffer on all sides of the system. Even though the box is
   #   orthorhombic, we will make it space group 1_a.
   $spaceGroup = "1_a";
   $spaceGroupName = "P1";
   $spaceGroupNum = 1;
   $spaceGroupSubNum = 1;

   # Define the angles.
   $angle[1] = $pi/2.0;
   $angle[2] = $pi/2.0;
   $angle[3] = $pi/2.0;

   # Compute the sine of each angle.
   &getAngleSine(\@angle);

   # Define the cell dimensions.
   foreach $axis (1..3)
      {$mag[$axis] = $maxPos[$axis] - $minPos[$axis] + $buffer;}

   # Define the real lattice parameters (a,b,c) in x,y,z vector form.
   &getABCVectors;

   # Generate the inverse of the real space lattice.
   &makeLatticeInv(\@realLattice,\@realLatticeInv,0);

   # Demand that there be no atoms with negative positions and that the
   #   system be centered.
   &shiftXYZCenter;
}


# Compute the magnitudes and angles between vectors in some lattice (real or
#   reciprocal). Originially this was written just for real lattices, but it
#   was adapted to serve reciprocal lattices too later on. Hence some of the
#   notes are written more in terms of real lattices. However, for reciprocal
#   lattices we use k_a, k_b, k_c, k_alpha, k_beta, and k_gamma in place of
#   a, b, c, alpha, beta, gamma. Otherwise everything is the same.
sub abcAlphaBetaGamma
{
   # Define local names for passed variables.
   my $lattice_ref = $_[0];
   my $mag_ref = $_[1];
   my $angle_ref = $_[2];
   my $angleDeg_ref = $_[3];
   my $getSineAngle = $_[4];

   # Declare local variables.
   my $axis;

   foreach $axis (1..3) # Loop over a,b,c
   {
      # Compute the magnitude of each vector.
      $mag_ref->[$axis]=sqrt($lattice_ref->[$axis][1]*$lattice_ref->[$axis][1]+
                             $lattice_ref->[$axis][2]*$lattice_ref->[$axis][2]+
                             $lattice_ref->[$axis][3]*$lattice_ref->[$axis][3]);
   }

   # Compute the angles (alpha,beta,gamma) between lattice vectors.
   # Angle = arccos(v1 dot v2) where v1,v2 are normalize vectors.
   # Alpha = angle between b and c.
   # Beta  = angle between a and c.
   # Gamma = angle between a and b.
   # All in radians.
   $angle_ref->[1] = &acos(($lattice_ref->[2][1]*$lattice_ref->[3][1] +
                           $lattice_ref->[2][2]*$lattice_ref->[3][2] +
                           $lattice_ref->[2][3]*$lattice_ref->[3][3])/
                           ($mag_ref->[2]*$mag_ref->[3]));
   $angle_ref->[2] = &acos(($lattice_ref->[1][1]*$lattice_ref->[3][1] +
                           $lattice_ref->[1][2]*$lattice_ref->[3][2] +
                           $lattice_ref->[1][3]*$lattice_ref->[3][3])/
                           ($mag_ref->[1]*$mag_ref->[3]));
   $angle_ref->[3] = &acos(($lattice_ref->[1][1]*$lattice_ref->[2][1] +
                           $lattice_ref->[1][2]*$lattice_ref->[2][2] +
                           $lattice_ref->[1][3]*$lattice_ref->[2][3])/
                           ($mag_ref->[1]*$mag_ref->[2]));
   $angleDeg_ref->[1] = $angle_ref->[1] * 180.0/$pi;
   $angleDeg_ref->[2] = $angle_ref->[2] * 180.0/$pi;
   $angleDeg_ref->[3] = $angle_ref->[3] * 180.0/$pi;

   # Compute the sine of each angle. Presently this only works for real
   #   lattices.
   if ($getSineAngle == 1)
      {&getAngleSine($angle_ref);}
}

sub getAngleSine
{
   # Define passed parameters.
   my $angle_ref = $_[0];

   # Define local variables.
   my $axis;

   # Compute the sine of each angle.
   foreach $axis (1..3)
      {$sineAngle[$axis] = sin($angle_ref->[$axis]);}
}

sub getABCVectors
{
   # If the space group has not been defined yet, then we will assume a space
   #   group of 1.
   if (length $spaceGroup == 0)
      {$spaceGroup = 1;}

   # Define local variables.
   my $abcAxis;
   my $xyzAxis;
   my $uniqueAxis;

   # Given the space group, define the lattice vectors according to W. Setyawan,
   #   S. Curtarolo, Comp. Mat. Sci. v. 49, pp. 299-312, (2010).

   # Initialize the realLattice[0][0] components to "".
   foreach $abcAxis (1..3)
      {$realLattice[$abcAxis][0] = "";}

   # Initially assume that the lattice will use the conventional cell.
   if ($spaceGroupNum <= 2) # Triclinic: alpha<beta<gamma then a<b<c if possible
   {

      # We first assume that a and x are co-axial.
      $realLattice[1][1] = $mag[1];
      $realLattice[1][2] = 0.0;
      $realLattice[1][3] = 0.0;
   
      # Then, b is in the x,y plane.
      $realLattice[2][1] = $mag[2] * cos($angle[3]);  # (b*cos(gamma))
      $realLattice[2][2] = $mag[2] * sin($angle[3]);  # (b*sin(gamma))
      $realLattice[2][3] = 0.0;
   
      # Then, c is a mix of all x,y,z directions.
      $realLattice[3][1] = $mag[3] * cos($angle[2]);  # (c*cos(beta))
      $realLattice[3][2] = $mag[3] * 
                           (cos($angle[1]) - cos($angle[3])*cos($angle[2])) /
                            sin($angle[3]);
      $realLattice[3][3] = $mag[3] * sqrt(1.0 - cos($angle[2])**2 -
                           ($realLattice[3][2]/$mag[3])**2);
   }
   elsif ($spaceGroupNum <= 15) # Monoclinic: a<=b<=c; alpha<90; beta,gamma=90
   {
      # Define the conditions in which a y-axis unique setting is used.
      if ((($spaceGroupNum == 3) && ($spaceGroupSubNum <= 2)) ||
          (($spaceGroupNum == 4) && ($spaceGroupSubNum <= 2)) ||
          (($spaceGroupNum == 5) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 6) && ($spaceGroupSubNum <= 2)) ||
          (($spaceGroupNum == 7) && ($spaceGroupSubNum <= 5)) ||
          (($spaceGroupNum == 8) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 9) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 10) && ($spaceGroupSubNum <= 2)) ||
          (($spaceGroupNum == 11) && ($spaceGroupSubNum <= 2)) ||
          (($spaceGroupNum == 12) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 13) && ($spaceGroupSubNum <= 5)) ||
          (($spaceGroupNum == 14) && ($spaceGroupSubNum <= 5)) ||
          (($spaceGroupNum == 15) && ($spaceGroupSubNum <= 7)))
         {$uniqueAxis = 2;} # y(b)-axis unique
      elsif ((($spaceGroupNum == 3) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 4) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 5) && ($spaceGroupSubNum <= 8)) ||
          (($spaceGroupNum == 6) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 7) && ($spaceGroupSubNum <= 10)) ||
          (($spaceGroupNum == 8) && ($spaceGroupSubNum <= 8)) ||
          (($spaceGroupNum == 9) && ($spaceGroupSubNum <= 8)) ||
          (($spaceGroupNum == 10) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 11) && ($spaceGroupSubNum <= 4)) ||
          (($spaceGroupNum == 12) && ($spaceGroupSubNum <= 8)) ||
          (($spaceGroupNum == 13) && ($spaceGroupSubNum <= 10)) ||
          (($spaceGroupNum == 14) && ($spaceGroupSubNum <= 10)) ||
          (($spaceGroupNum == 15) && ($spaceGroupSubNum <= 14)))
         {$uniqueAxis = 3;} # z(c)-axis unique
      else
         {$uniqueAxis = 1;} # x(a)-axis unique

      if ($uniqueAxis == 1) # y(b)-axis unique
      {
         # Define the a axis as (a*sin(beta), 0, a*cos(beta))
         $realLattice[1][1] = $mag[1] * sin($angle[2]);
         $realLattice[1][2] = 0.0;
         $realLattice[1][3] = $mag[1] * cos($angle[2]);

         # Define the b axis as (0, b, 0)
         $realLattice[2][1] = 0.0;
         $realLattice[2][2] = $mag[2];
         $realLattice[2][3] = 0.0;
   
         # Define the c axis as (0, 0, c)
         $realLattice[3][1] = 0.0;
         $realLattice[3][2] = 0.0;
         $realLattice[3][3] = $mag[3];
      }
      elsif ($uniqueAxis == 2) # z(c)-axis unique
      {
         # Define the a axis as (a, 0, 0)
         $realLattice[1][1] = $mag[1];
         $realLattice[1][2] = 0.0;
         $realLattice[1][3] = 0.0;
   
         # Define the b axis as (b*cos(gamma), b*sin(gamma), 0)
         $realLattice[2][1] = $mag[2] * cos($angle[3]);
         $realLattice[2][2] = $mag[2] * sin($angle[3]);
         $realLattice[2][3] = 0.0;
   
         # Define the c axis as (0, 0, c)
         $realLattice[3][1] = 0.0;
         $realLattice[3][2] = 0.0;
         $realLattice[3][3] = $mag[3];
      }
      else # ($uniqueAxis == 3) # x(a)-axis unique
      {
         # Define the a axis as (a, 0, 0)
         $realLattice[1][1] = $mag[1];
         $realLattice[1][2] = 0.0;
         $realLattice[1][3] = 0.0;
   
         # Define the b axis as (0, b, 0)
         $realLattice[2][1] = 0.0;
         $realLattice[2][2] = $mag[2];
         $realLattice[2][3] = 0.0;
   
         # Define the c axis as (0, c*cos(alpha), c*sin(alpha))
         $realLattice[3][1] = 0.0;
         $realLattice[3][2] = $mag[3] * cos($angle[1]);
         $realLattice[3][3] = $mag[3] * sin($angle[1]);
      }
   }
   elsif ($spaceGroupNum <= 74) # Orthorhombic: a<b<c; alpha=beta=gamma=90
   {
      # Define the a axis as (a, 0, 0)
      $realLattice[1][1] = $mag[1];
      $realLattice[1][2] = 0.0;
      $realLattice[1][3] = 0.0;

      # Define the b axis as (0, b, 0)
      $realLattice[2][1] = 0.0;
      $realLattice[2][2] = $mag[2];
      $realLattice[2][3] = 0.0;

      # Define the c axis as (0, 0, c)
      $realLattice[3][1] = 0.0;
      $realLattice[3][2] = 0.0;
      $realLattice[3][3] = $mag[3];
   }
   elsif ($spaceGroupNum <= 142) # Tetragonal: a=b/=c; alpha=beta=gamma=90
   {
      # Define the a axis as (a, 0, 0)
      $realLattice[1][1] = $mag[1];
      $realLattice[1][2] = 0.0;
      $realLattice[1][3] = 0.0;

      # Define the b axis as (0, b, 0)
      $realLattice[2][1] = 0.0;
      $realLattice[2][2] = $mag[2];
      $realLattice[2][3] = 0.0;

      # Define the c axis as (0, 0, c)
      $realLattice[3][1] = 0.0;
      $realLattice[3][2] = 0.0;
      $realLattice[3][3] = $mag[3];
   }
   elsif ($spaceGroupNum <= 167) # Trigonal: Two Options:
   {
      # Hexagonal cell: a=b/=c; alpha=beta=90, gamma=120 occurs for P,
      #   H-centered (made from three hexagonal cells) and R with hexagonal
      #   axes (made from three primitive rhombohedral cells).
      if (($latticeType eq "P") || ($latticeType eq "H") ||
          (($latticeType eq "R") && ($spaceGroupSubNum == 2))) # Hexagonal cell
      {
         # Define the a axis as (a/2, -(a*sqrt(3))/2, 0)
         $realLattice[1][1] = $mag[1]/2.0;
         $realLattice[1][2] = -($mag[1]*sqrt(3.0)/2.0);
         $realLattice[1][3] = 0.0;

         # Define the b axis as (b/2,  (b*sqrt(3))/2, 0)
         $realLattice[2][1] = $mag[2]/2.0;
         $realLattice[2][2] =  ($mag[2]*sqrt(3.0)/2.0);
         $realLattice[2][3] = 0.0;

         # Define the c axis as (0, 0, c)
         $realLattice[3][1] = 0.0;
         $realLattice[3][2] = 0.0;
         $realLattice[3][3] = $mag[3];
      }
      # Rhombohedral cell: a=b=c; alpha=beta=gamma/=90
      else # Primitive Rhombohedral cell
      {
         # Define the a axis as (a*cos(alpha/2), -a*sin(alpha/2), 0)
         $realLattice[1][1] =  $mag[1]*cos($angle[1]/2.0);
         $realLattice[1][2] = -$mag[1]*sin($angle[1]/2.0);
         $realLattice[1][3] =  0.0;

         # Define the b axis as (a*cos(alpha/2), a*sin(alpha/2), 0)
         $realLattice[2][1] = $mag[1]*cos($angle[1]/2.0);
         $realLattice[2][2] = $mag[1]*sin($angle[1]/2.0);
         $realLattice[2][3] = 0.0;

         # Define the c axis as: (a*cos(alpha)/cos(alpha/2), 0,
         #   a*sqrt(1 - cos(alpha)**2 / cos(alpha/2)**2))
         $realLattice[3][1] = $mag[1]*cos($angle[1])/cos($angle[1]/2.0);
         $realLattice[3][2] = 0.0;
         $realLattice[3][3] = $mag[1]*sqrt(1.0 - cos($angle[1])**2 / 
               cos($angle[1]/2.0)**2);
      }
   }
   elsif ($spaceGroupNum <= 194) # Hexagonal: a=b/=c; alpha=beta=90; gamma=120
   {
      # Define the a axis as (a/2, -(a*sqrt(3))/2, 0)
      $realLattice[1][1] = $mag[1]/2.0;
      $realLattice[1][2] = -($mag[1]*sqrt(3.0)/2.0);
      $realLattice[1][3] = 0.0;

      # Define the b axis as (b/2,  (b*sqrt(3))/2, 0)
      $realLattice[2][1] = $mag[2]/2.0;
      $realLattice[2][2] =  ($mag[2]*sqrt(3.0)/2.0);
      $realLattice[2][3] = 0.0;

      # Define the c axis as (0, 0, c)
      $realLattice[3][1] = 0.0;
      $realLattice[3][2] = 0.0;
      $realLattice[3][3] = $mag[3];
   }
   elsif ($spaceGroupNum <= 230) # Cubic: a=b=c; alpha=beta=gamma=90
   {
      # Define the a axis as (a, 0, 0)
      $realLattice[1][1] = $mag[1];
      $realLattice[1][2] = 0.0;
      $realLattice[1][3] = 0.0;

      # Define the b axis as (0, b, 0)
      $realLattice[2][1] = 0.0;
      $realLattice[2][2] = $mag[2];
      $realLattice[2][3] = 0.0;

      # Define the c axis as (0, 0, c)
      $realLattice[3][1] = 0.0;
      $realLattice[3][2] = 0.0;
      $realLattice[3][3] = $mag[3];
   }

   
   # Correct numerical errors.
   foreach $abcAxis (1..3)
   {
      foreach $xyzAxis (1..3)
      {
         if (abs($realLattice[$abcAxis][$xyzAxis]) < $epsilon)
            {$realLattice[$abcAxis][$xyzAxis] = 0.0;}
      }
   }
}


sub makeLatticeInv
{
   # Define passed parameters.
   my $inLattice_ref = $_[0];
   my $outLattice_ref = $_[1];
   my $piFactor = $_[2];

   # Define local variables.
   my $string;
   my $matrix;
   my $matrixInv;
   my $axisABC;
   my $axisXYZ;

   # Initialize the matrix string.
   $string = "";

   # Copy the in lattice matrix to $matrix.
   foreach $axisABC (1..3)
      {$string = $string . "\[ $inLattice_ref->[$axisABC][1] ".
                              "$inLattice_ref->[$axisABC][2] ".
                              "$inLattice_ref->[$axisABC][3] \]\n";}
   $matrix = Math::MatrixReal->new_from_string($string);

   # Invert it.
   $matrixInv = $matrix->inverse();

   # Copy it to the standard form used in the rest of this program and include
   #   the factor of 2Pi if requested.
   foreach $axisXYZ (1..3)
   {
      foreach $axisABC (1..3)
      {
         $outLattice_ref->[$axisXYZ][$axisABC] =
               $matrixInv->element($axisXYZ,$axisABC);
         if ($piFactor == 1)
            {$outLattice_ref->[$axisXYZ][$axisABC] *= (2.0 * $pi);}
         elsif ($piFactor == -1)
            {$outLattice_ref->[$axisXYZ][$axisABC] /= (2.0 * $pi);}
      }
   }
}


sub computeBZ_WS
{
   # Define passed parameters.
   my $inLattice_ref = $_[0];
   my $zoneLevel = $_[1];
   my $zone_ref = $_[2];

   # Define local variables.

   # 
}

# Given some (real/reciprocal) set of lattice vectors, this subroutine will
#   produce the lattice volume.
sub makeLatticeVolume
{
   # Define passed paramenters.
   my $inLattice_ref = $_[0];

   # Define local variables.
   my $i;
   my $iCycle1;
   my $iCycle2;
   my $j;
   my $jCycle1;
   my $jCycle2;
   my @tempLattice;
   my $cellVolume;

   for ($i=1;$i<=3;$i++)
   {
      $cellVolume = 0;
      $iCycle1 = ($i % 3) + 1;
      $iCycle2 = (($i+1) % 3) + 1;
      for ($j=1;$j<=3;$j++)
      {
         $jCycle1 = ($j % 3) + 1;
         $jCycle2 = (($j+1) % 3) + 1;
         $tempLattice[$j][$i] = $inLattice_ref->[$jCycle1][$iCycle1] *
                                $inLattice_ref->[$jCycle2][$iCycle2] -
                                $inLattice_ref->[$jCycle2][$iCycle1] *
                                $inLattice_ref->[$jCycle1][$iCycle2];
         $cellVolume = $cellVolume + $inLattice_ref->[$j][$i] *
                                     $tempLattice[$j][$i];
      }
      for ($j=1;$j<=3;$j++)
         {$tempLattice[$j][$i] = 2.0 * $pi * $tempLattice[$j][$i] / 
                                 $cellVolume;}  # In Angstroms!
   }

   return $cellVolume;
}


sub getMinMaxXYZ
{
   # Define local variables.
   my $atom;
   my $axis;

   $maxPos[1] = -$bigReal;  $minPos[1] = $bigReal;  # x
   $maxPos[2] = -$bigReal;  $minPos[2] = $bigReal;  # y
   $maxPos[3] = -$bigReal;  $minPos[3] = $bigReal;  # z
   foreach $atom (1..$numAtoms)
   {
      foreach $axis (1..3)
      {
         if ($directXYZ[$atom][$axis] > $maxPos[$axis])
            {$maxPos[$axis] = $directXYZ[$atom][$axis];}
         if ($directXYZ[$atom][$axis] < $minPos[$axis])
            {$minPos[$axis] = $directXYZ[$atom][$axis];}
      }
   }
}
# Get the direct XYZ coordinates of an atom given the fractional ABC ones.
sub getDirectXYZ
{
   # Define passed parameters.
   my $currentAtom = $_[0];

   # Declare local variables.
   my @fractCoordsABC;
   my @directCoordsXYZ;
   my @values;

   $fractCoordsABC[0] = $fractABC[$currentAtom][1];
   $fractCoordsABC[1] = $fractABC[$currentAtom][2];
   $fractCoordsABC[2] = $fractABC[$currentAtom][3];

   @values = &getDirectXYZPoint(@fractCoordsABC);

   $directXYZ[$currentAtom][1] = $values[1];
   $directXYZ[$currentAtom][2] = $values[2];
   $directXYZ[$currentAtom][3] = $values[3];
}

sub getDirectXYZPoint
{
   # Define passed parameter.
   my @Pabc;

   $Pabc[1] = $_[0];
   $Pabc[2] = $_[1];
   $Pabc[3] = $_[2];

   # Declare local variables.
   my @Pxyz;
   my $xyzAxis;
   my $abcAxis;

   # Pxyz = atom position in xyz orthogonal coordinates.
   # Pabc = atom position in abc fractional coordinates.
   # ax ay az = xyz component coefficients of a lattice vector.
   # bx by bz = xyz component coefficients of b lattice vector.
   # cx cy cz = xyz component coefficients of c lattice vector.

   # Px = (Pa*ax + Pb*bx + Pc*cx) x
   # Py = (Pa*ay + Pb*by + Pc*cy) y
   # Pz = (Pa*az + Pb*bz + Pc*cz) z

   # Loop over x,y,z orthogonal axes.
   foreach $xyzAxis (1..3)
   {
      # Initialize the x,y,z coordinate of this atom.
      $Pxyz[$xyzAxis] = 0.0;

      # Loop over a,b,c fractional directions.
      foreach $abcAxis (1..3)
      {
         # Example:  point[x] = sum fractABC[a,b,c] * latt[a,b,c][x]
         $Pxyz[$xyzAxis] += $Pabc[$abcAxis] * $realLattice[$abcAxis][$xyzAxis];
      }
   }

   return (@Pxyz);
}

sub getDirectABC
{
   # Define passed parameters.
   my $currentAtom = $_[0];

   # Declare local variables.
   my $abcAxis;
   my $xyzAxis;


   # Pxyz = atom position in xyz orthogonal coordinates.
   # Pabc = atom position in abc direct space coordinates.
   # L-1  = Real lattice matrix inverse.
   # Pabc = Pxyz L-1

   foreach $abcAxis (1..3)
   {
      # Initialize the directABC coordinate for this axis.
      $directABC[$currentAtom][$abcAxis] = 0;

      foreach $xyzAxis (1..3)
      {
         $directABC[$currentAtom][$abcAxis] +=
               $directXYZ[$currentAtom][$xyzAxis] *
               $realLatticeInv[$xyzAxis][$abcAxis];
      }

      # Note that we now have the fractional ABC coordinates but we are not
      #   going to do anything with them now because I don't want to mess stuff
      #   up by making things confusing (any more than they already are).

      # Multiply by the lattice magnitude to get the direct space ABC coords.
      $directABC[$currentAtom][$abcAxis] *= $mag[$abcAxis];
   }
}

sub directXYZ2fractABC
{
   # Define passed parameters.
   my $posXYZ_ref = $_[0];
   my $realLatticeInv_ref = $_[1];

   # Define local variables.
   my @posABC;
   my $abcAxis;
   my $xyzAxis;

   # Pxyz = atom position in xyz orthogonal coordinates.
   # Pabc = atom position in abc direct space coordinates.
   # L-1  = Real lattice matrix inverse.
   # Pabc = Pxyz L-1

   foreach $abcAxis (1..3)
   {
      # Initialize the direct space ABC coordinate for this axis.
      $posABC[$abcAxis] = 0;

      foreach $xyzAxis (1..3)
      {
         $posABC[$abcAxis] += $posXYZ_ref->[$xyzAxis] *
               $realLatticeInv_ref->[$xyzAxis][$abcAxis];
      }
   }

   return (@posABC[1..3]);
}

sub directXYZ2directABC
{
   # Define passed parameters.
   my $posXYZ_ref = $_[0];
   my $realLatticeInv_ref = $_[1];

   # Define local variables.
   my @posABC;
   my $abcAxis;
   my $xyzAxis;

   # Pxyz = atom position in xyz orthogonal coordinates.
   # Pabc = atom position in abc direct space coordinates.
   # L-1  = Real lattice matrix inverse.
   # Pabc = Pxyz L-1

   foreach $abcAxis (1..3)
   {
      # Initialize the direct space ABC coordinate for this axis.
      $posABC[$abcAxis] = 0;

      foreach $xyzAxis (1..3)
      {
         $posABC[$abcAxis] += $posXYZ_ref->[$xyzAxis] *
               $realLatticeInv_ref->[$xyzAxis][$abcAxis];
      }

      # Multiply by the lattice magnitude to get the direct space ABC coords.
      $posABC[$abcAxis] *= $mag[$abcAxis];
   }

   return (@posABC[1..3]);
}

sub fractABC2directXYZ
{
   # Define passed parameters.
   my $posABC_ref = $_[0];  # Fractional ABC coordinates.
   my $realLattice_ref = $_[1];

   # Declare local variables.
   my @posXYZ;
   my $xyzAxis;
   my $abcAxis;

   # Pxyz = atom position in xyz orthogonal coordinates.
   # Pabc = atom position in abc fractional coordinates.
   # ax ay az = xyz component coefficients of a lattice vector.
   # bx by bz = xyz component coefficients of b lattice vector.
   # cx cy cz = xyz component coefficients of c lattice vector.

   # Px = (Pa*ax + Pb*bx + Pc*cx) x
   # Py = (Pa*ay + Pb*by + Pc*cy) y
   # Pz = (Pa*az + Pb*bz + Pc*cz) z

   # Loop over x,y,z orthogonal axes.
   foreach $xyzAxis (1..3)
   {
      # Initialize the x,y,z coordinate of this atom.
      $posXYZ[$xyzAxis] = 0.0;

      # Loop over a,b,c fractional directions.
      foreach $abcAxis (1..3)
      {
         # Example:  atom[x][1] = sum fractABC[a,b,c][1] * latt[a,b,c][x]
         $posXYZ[$xyzAxis] += $posABC_ref->[$abcAxis] *
               $realLattice_ref->[$abcAxis][$xyzAxis];
      }
   }
}

sub getFractABC
{
   # Declare local variables.
   my $currentAtom = $_[0];
   my $axis;

   # Initialize the fractABC[$currentAtom][0] to "" for easy printing later.
   $fractABC[$currentAtom][0] = "";

   # Use the existing direct space ABC coordinates and the vector magnitudes
   #   to get the fractional ABC coordinates.
   foreach $axis (1..3)
   {
      $fractABC[$currentAtom][$axis] = $directABC[$currentAtom][$axis] /
            $mag[$axis];
   }
}

sub min
{
   if ($_[0] < $_[1])
      {return $_[0];}
   else
      {return $_[1];}
}


# This subroutine will translate linearly all atoms along each orthogonal axis
#   to make sure that the system (as a whole) is centered in the simulation
#   box.  This subroutine is only usable once the minimum position for each
#   axis has been found and is really only applicable to molecular systems.
sub shiftXYZCenter
{
   # Declare local variables.
   my $atom;
   my $axis;
   my @centerPos; # This is the Atom System Center
   my @centerCell; # This is the center of the simulation cell.
   my @centerDiff; # The difference between the atom system and cell centers.

   # Compute the midpoint of the set of atoms on the basis of the max and min
   #   positions of atoms. AtomSysCenter = Min + (Max-Min)/2 = (Max+Min)/2
   foreach $axis (1..3)
      {$centerPos[$axis] = ($maxPos[$axis] + $minPos[$axis])/2.0;}

   # Compute the midpoint of the simulation cell.
   foreach $axis (1..3)
      {$centerCell[$axis] = $mag[$axis]/2.0;}

   # Compute the difference between the cell center and the atom system center.
   #   Diff = AtomSysCenter - CellCenter
   foreach $axis (1..3)
      {$centerDiff[$axis] = ($centerPos[$axis] - $centerCell[$axis]);}

   # Recompute the direct XYZ coordiantes of each atom.
   foreach $atom (1..$numAtoms)
   {
      foreach $axis (1..3)
      {
         # Shift the current atom's coordiante by the difference between the
         #   atom system center and the cell center. This will center the
         #   system of atoms about the center of the cell.
         $directXYZ[$atom][$axis] -= $centerDiff[$axis];
      }
   }
}

sub translateAtoms
{
   # Define passed parameters.
   my $translation_ref = $_[0];

   # Declare local variables.
   my $atom;
   my $axis;

   # If the translation request uses the special 0 0 0 designation, then shift
   #   the atoms to the center of the cell. Otherwise, shift that atoms
   #   according to the requested translation.
   if (($translation_ref->[1] == 0) && ($translation_ref->[2] == 0) &&
       ($translation_ref->[3] == 0))
   {
      # Find the center and then shift the atoms to the center.
      &getMinMaxXYZ;
      &shiftXYZCenter;

      # Convert to direct a,b,c and then to fractional a,b,c.
      foreach $atom (1..$numAtoms)
      {
         &getDirectABC($atom);
         &getFractABC($atom);
      }
   }
   else
   {
      foreach $atom (1..$numAtoms)
      {
         foreach $axis (1..3)
            {$directXYZ[$atom][$axis] += $translation_ref->[$axis];}

         # Convert to direct a,b,c and then to fractional a,b,c.
         &getDirectABC($atom);
         &getFractABC($atom);
      }
   }

   # Make sure that all atoms are inside the simulation box.
   &checkBoundingBox(0);
}

sub insertVacuum
{
   # Define passed parameters.
   my $vacAxis = $_[0];
   my $vacAmt  = $_[1];

   # Define local variables.
   my $atom;

   # Increase the magnitude of the requested lattice direction.
   $mag[$vacAxis] += $vacAmt;

   # Get new ABC vectors (in XYZ form) based on the new lattice vectors.
   &getABCVectors;

   # Reobtain the inverse lattice and the reciprocal lattice vectors.
   &makeLatticeInv(\@realLattice,\@realLatticeInv,0);
   &makeLatticeInv(\@realLattice,\@recipLattice,1);
   &abcAlphaBetaGamma(\@realLattice,\@mag,\@angle,\@angleDeg,1);

   # Get new fractional ABC atomic positions based on the old direct space
   #   ABC atomic positions which are still the same.  Note that the direct
   #   space XYZ atomic coordinates are also still the same.
   foreach $atom (1..$numAtoms)
      {&getFractABC($atom);}
}

# The purpose of this subroutine is to remove a designated set of atoms from
#   the currently stored model. The basic algorithm is to make a duplicate
#   list of all the atoms that contains only those that should not be cut.
#   Then, any data structures that need to be updated will be.
sub cutBlock
{
   # Define local variables for passed parameters.
   my $zone = $_[0];
   my $abcxyzFlag = $_[1];
   my $blockBorder_ref = $_[2];

   # Define local variables.
   my $atom;
   my $axis;
   my $cut;
   my $foundInside;
   my $newNumAtoms;        # Local version to replace old one.
   my @newAtomElementName; # Local version to replace old one.
   my @newAtomElementID;   # Local version to replace old one.
   my @newAtomSpeciesID;   # Local version to replace old one.
   my @newFractABC;        # Local version to replace old one.
   my @newDirectABC;       # Local version to replace old one.
   my @newDirectXYZ;       # Local version to replace old one.
   my @newAtomTag;         # Local version to replace old one.
   my @currentCoords;

   # Note that the block borders could have been defined to have a letter as
   #   the "to" value. This indicates a request for the maximum. Check and
   #   modify accordingly.
   if ($blockBorder_ref->[1][2] eq "a")
      {$blockBorder_ref->[1][2] = $mag[1];}
   if ($blockBorder_ref->[2][2] eq "b")
      {$blockBorder_ref->[2][2] = $mag[2];}
   if ($blockBorder_ref->[3][2] eq "c")
      {$blockBorder_ref->[3][2] = $mag[3];}

   # Initialize the new count of the number of atoms.
   $newNumAtoms = 0;

   # Loop over each atom and determine if it should be kept or not.
   foreach $atom (1..$numAtoms)
   {
      # Set a variable that defines whether or not this atom is found inside
      #   the selected region.  The default is to assume it is inside.
      $foundInside = 1;

      # We can check versus abc or xyz axes. Make that determination here and
      #   gather the coordinates we will use for comparison.
      if ($abcxyzFlag == 0) # abc
      {
         foreach $axis (1..3)
            {$currentCoords[$axis] = $directABC[$atom][$axis];}
      }
      else
      {
         foreach $axis (1..3)
            {$currentCoords[$axis] = $directXYZ[$atom][$axis];}
      }

      # Determine if the current atom is "inside" the block.
      foreach $axis (1..3)
      {
         if (($currentCoords[$axis] < $blockBorder_ref->[$axis][1]) ||
             ($currentCoords[$axis] > $blockBorder_ref->[$axis][2]))
            {$foundInside = 0;}
      }

      # If the inside zone (0) is the zone to cut and the atom is found inside
      #   then we turn on $cut. If the outside zone is the zone to cut, and
      #   the atom is found outside, the we also turn on $cut. Otherwise, we
      #   leave $cut off.
      if (($zone == 0) && ($foundInside == 1)) # Cutting inside, found inside
         {$cut = 1;}
      elsif (($zone == 1) && ($foundInside == 0)) # Cutting outside, found out
         {$cut = 1;}
      else
         {$cut = 0;}

      # Act, if we have decided to keep this atom.
      if ($cut == 0)
      {
         $newNumAtoms++;

         foreach $axis (1..3)
         {
            $newFractABC[$newNumAtoms][$axis]  = $fractABC[$atom][$axis];
            $newDirectABC[$newNumAtoms][$axis] = $directABC[$atom][$axis];
            $newDirectXYZ[$newNumAtoms][$axis] = $directXYZ[$atom][$axis];
         }
         $newAtomElementName[$newNumAtoms] = $atomElementName[$atom];
         $newAtomElementID[$newNumAtoms]   = $atomElementID[$atom];
         $newAtomSpeciesID[$newNumAtoms]   = $atomSpeciesID[$atom];
         $newAtomTag[$newNumAtoms]         = $atomTag[$atom];
      }
   }

   # Update the system information.
   undef (@fractABC);
   undef (@directABC);
   undef (@directXYZ);
   undef (@atomElementName);
   undef (@atomElementID);
   undef (@atomSpeciesID);
   undef (@atomTag);
   $numAtoms        = $newNumAtoms;
   @fractABC        = @newFractABC;
   @directABC       = @newDirectABC;
   @directXYZ       = @newDirectXYZ;
   @atomElementName = @newAtomElementName;
   @atomElementID   = @newAtomElementID;
   @atomSpeciesID   = @newAtomSpeciesID;
   @atomTag         = @newAtomTag;
   undef (@newFractABC);
   undef (@newDirectABC);
   undef (@newDirectXYZ);
   undef (@newAtomElementName);
   undef (@newAtomElementID);
   undef (@newAtomSpeciesID);
   undef (@newAtomTag);
}

# The purpose of this subroutine is to remove a designated set of atoms from
#   the currently stored model. The basic algorithm is to make a duplicate
#   list of all the atoms that contains only those that should not be cut.
#   Then, any data structures that need to be updated will be.
sub cutSphere
{
   # Define local variables for passed parameters.
   my $zone          = $_[0];
   my $abcxyzFlag    = $_[1];
   my $sphereRad     = $_[2];
   my $targetAtom    = $_[3];
   my $sphereLoc_ref = $_[4];

   # Define local variables.
   my $atom;
   my $axis;
   my $cut;
   my $foundInside;
   my $distance;
   my $newNumAtoms;        # Local version to replace old one.
   my @newAtomElementName; # Local version to replace old one.
   my @newAtomElementID;   # Local version to replace old one.
   my @newAtomSpeciesID;   # Local version to replace old one.
   my @newFractABC;        # Local version to replace old one.
   my @newDirectABC;       # Local version to replace old one.
   my @newDirectXYZ;       # Local version to replace old one.
   my @newAtomTag;         # Local version to replace old one.
   my @currentCoords;

   # In the case that a particular atom site was give, we will define the
   #   sphere location in XYZ coordinates according to the atom position.
   if ($targetAtom != 0)
   {
      # Copy the XYZ coordinates of the targeted atom into the sphere location.
      foreach $axis (1..3)
         {$sphereLoc_ref->[$axis] = $directXYZ[$targetAtom][$axis];}

      # Make sure that the program knows we want to use XYZ coordinates.
      $abcxyzFlag = 1;
   }

   # Initialize the new count of the number of atoms.
   $newNumAtoms = 0;

   # Loop over each atom and determine if it should be kept or not.
   foreach $atom (1..$numAtoms)
   {
      # Set a variable that defines whether or not this atom is found inside
      #   the selected region.  The default is to assume it is inside.
      $foundInside = 1;

      # We can check versus abc or xyz axes. Make that determination here and
      #   gather the coordinates we will use for comparison.
      if ($abcxyzFlag == 0) # abc
      {
         foreach $axis (1..3)
            {$currentCoords[$axis] = $directABC[$atom][$axis];}
      }
      else
      {
         foreach $axis (1..3)
            {$currentCoords[$axis] = $directXYZ[$atom][$axis];}
      }

      # Determine if the current atom is "inside" the sphere. First compute
      #   the distance and then compare to the sphere radius.
      $distance = sqrt(($currentCoords[1]-$sphereLoc_ref->[1])**2 +
                       ($currentCoords[2]-$sphereLoc_ref->[2])**2 +
                       ($currentCoords[3]-$sphereLoc_ref->[3])**2);
      if ($distance > $sphereRad)
         {$foundInside = 0;}

      # If the inside zone (0) is the zone to cut and the atom is found inside
      #   then we turn on $cut. If the outside zone is the zone to cut, and
      #   the atom is found outside, the we also turn on $cut. Otherwise, we
      #   leave $cut off.
      if (($zone == 0) && ($foundInside == 1)) # Cutting inside, found inside
         {$cut = 1;}
      elsif (($zone == 1) && ($foundInside == 0)) # Cutting outside, found out
         {$cut = 1;}
      else
         {$cut = 0;}

      # Act, if we have decided to keep this atom.
      if ($cut == 0)
      {
         $newNumAtoms++;

         foreach $axis (1..3)
         {
            $newFractABC[$newNumAtoms][$axis]  = $fractABC[$atom][$axis];
            $newDirectABC[$newNumAtoms][$axis] = $directABC[$atom][$axis];
            $newDirectXYZ[$newNumAtoms][$axis] = $directXYZ[$atom][$axis];
         }
         $newAtomElementName[$newNumAtoms] = $atomElementName[$atom];
         $newAtomElementID[$newNumAtoms]   = $atomElementID[$atom];
         $newAtomSpeciesID[$newNumAtoms]   = $atomSpeciesID[$atom];
         $newAtomTag[$newNumAtoms]         = $atomTag[$atom];
      }
   }

   # Update the system information.
   undef (@fractABC);
   undef (@directABC);
   undef (@directXYZ);
   undef (@atomElementName);
   undef (@atomElementID);
   undef (@atomSpeciesID);
   undef (@atomTag);
   $numAtoms        = $newNumAtoms;
   @fractABC        = @newFractABC;
   @directABC       = @newDirectABC;
   @directXYZ       = @newDirectXYZ;
   @atomElementName = @newAtomElementName;
   @atomElementID   = @newAtomElementID;
   @atomSpeciesID   = @newAtomSpeciesID;
   @atomTag         = @newAtomTag;
   undef (@newFractABC);
   undef (@newDirectABC);
   undef (@newDirectXYZ);
   undef (@newAtomElementName);
   undef (@newAtomElementID);
   undef (@newAtomSpeciesID);
   undef (@newAtomTag);
}

# The goal of this subroutine is to take a crystal in its given lattice and
#   deduce from it a set of orthorhombic lattice parameters that retains
#   periodic boundary conditions.  Essentially we will make a=x, b=y, and c=z
#   while retaining the magnitudes of the ax, by, and cz components for a, b,
#   and c.  This can only be done for hexagonal lattices now that have already
#   been doubled (via the -sc option) in the direction of the lattice vector
#   that will be modified.  (i.e. you have to plan ahead to use this, it won't
#   do everything automatically.)
sub makeOrtho
{
   # Define local variables.
   my $axisABC;
   my $axisXYZ;
   my $atom;
   my $axis;

   # Clearly, all the angles will have to be 90 degrees (= pi/2 radians).
   foreach $axis (1..3)
   {
      $angle[$axis]    = $pi/2.0;
      $angleDeg[$axis] = 90.0;
   }

   # Further, the magnitudes of the ABC vectors will be the values currently
   #   in the ax, by, and cz indices.
   foreach $axisABC (1..3)
      {$mag[$axisABC] = $realLattice[$axisABC][$axisABC];}

   # Eliminate components that are off axis (i.e. not ax, by, or cz.)
   foreach $axisABC (1..3)
   {
      foreach $axisXYZ (1..3)
      {
         # The ABC axis and XYZ axis *should* be colinear so we don't need to
         #   adjust this vector.
         if ($axisABC==$axisXYZ)
            {next;}

         # It is necessary to make this ABC,XYZ axis pair orthogonal.
         $realLattice[$axisABC][$axisXYZ] = 0.0;
      }
   }

   # Obtain the inverse of the real lattice.  This must be done now so that
   #   we can get the abc coordinates of atoms if we are given xyz.  It must
   #   be done again later after applying the spacegroup and supercell.
   &makeLatticeInv(\@realLattice,\@realLatticeInv,0);

   # Finally, it is necessary to make sure that all atoms in the model are
   #   within the simulation box.  This is done by first recomputing the
   #   directABC and fractABC atomic coordinates from the directXYZ coordinates
   #   and the newly defined lattice.  Then we simply call the checkBoundingBox
   #   subroutine above to do the rest.
   foreach $atom (1..$numAtoms)
   {
      &getDirectABC($atom);
      &getFractABC($atom);
   }

   &checkBoundingBox(0);
}

sub applyPerturbation
{
   # Define passed parameters
    my $maxPertMag = $_[0];

   # Define local variables.
   my $numAtoms;
   my $atom;
   my $axis;
   my @perturbRThetaPhi; # [1]=R; [2]=Theta; [3]=Phi
   my @perturbXYZ;  # [1] = X; [2] = Y; [3] = Z;

   # For each atom in the model, we will apply the necessary translation.  Note
   #   that these values will be adjusted for periodicity and will be
   #   propogated to the other representations in the StructureControl module.
   foreach $atom (1..$numAtoms)
   {
      # Compute the XYZ perturbation for this atom based on random R, Theta,
      #   and Phi values where the maximum R value was defined on the command
      #   line in angstroms.
      $perturbRThetaPhi[1] = rand($maxPertMag);
      $perturbRThetaPhi[2] = rand($pi);
      $perturbRThetaPhi[3] = rand(2.0*$pi);
      $perturbXYZ[1] = $perturbRThetaPhi[1] * sin($perturbRThetaPhi[2]) *
                       cos($perturbRThetaPhi[3]);
      $perturbXYZ[2] = $perturbRThetaPhi[1] * sin($perturbRThetaPhi[2]) *
                       sin($perturbRThetaPhi[3]);
      $perturbXYZ[3] = $perturbRThetaPhi[1] * cos($perturbRThetaPhi[2]);

      foreach $axis (1..3)
         {$directXYZ[$atom][$axis] += $perturbXYZ[$axis];}

      # Propogate to other representations.
      &getDirectABC($atom);
      &getFractABC($atom);
   }

   # Make sure that all atoms are inside the simulation box.
   &checkBoundingBox(0);
}

sub applyFilter
{
   # Define passed parameters.
   my $minDistFilter = $_[0];

   # Define local variables.
   my $axis;
   my $atom;
   my $neighborAtom;
   my $rejectAtom;
   my $newNumAtoms;
   my $numRejected;
   my @rejectedAtoms;
   # Local collection of data for forming the new atom list.
   my @fractABCLocal;
   my @atomElementNameLocal;
   my @atomElementIDLocal;
   my @atomSpeciesIDLocal;
   # This subroutine does not preserve all the information in the
   #   StructureControl module.  In particular, the maps and some extended
   #   positions will not be preserved and may need to be recalculated if they
   #   are going to be used for anything immediately after this.

   # Set the limitDist value for periodic cell interaction using the
   #   minDistFilter interaction factor for this operation.
   &setLimitDist($minDistFilter);

   # Create the minimal distance matrix.
   &createMinDistMatrix;

   # Initialize the new number of atoms in the model.
   $newNumAtoms = 0;

   # Initialize the number of atoms that have been rejected from the model.
   $numRejected = 0;

   # Consider each atom in turn and determine if it has any neighbors that
   #   are too close.  Any neighbor that is too close is recorded and
   #   considered to be removed from the system.
ATOM: foreach $atom (1..$numAtoms)
   {
      # Print a running count of progress on the screen so the user does not
      #   get too bored waiting.  This process can take a while for large
      #   systems.
      if ($atom%10 == 0)
         {print STDOUT "|";}
      else
         {print STDOUT ".";}

      if ($atom%50 == 0)
         {print STDOUT " $atom\n";}

      # Determine if this atom has already been rejected.  If so, then we skip
      #   to the next atom.
      foreach $rejectAtom (1..$numRejected)
      {
         if ($rejectedAtoms[$rejectAtom] == $atom)
            {next ATOM;}
      }

      # Increment the number of atoms in the new model and record the pertanent
      #   information to the local copies of data.
      $newNumAtoms++;
      foreach $axis (1..3)
         {$fractABCLocal[$newNumAtoms][$axis] = $fractABC[$atom][$axis];}
      $atomElementNameLocal[$newNumAtoms] = $atomElementName[$atom];
      $atomElementIDLocal[$newNumAtoms] = $atomElementID[$atom];
      $atomSpeciesIDLocal[$newNumAtoms] = $atomSpeciesID[$atom];

      # Determine if there are any atoms that are too close to the current one
      #   and which must therefore be rejected from the model.
      foreach $neighborAtom ($atom+1..$numAtomsExt)
      {
         if ($minDist[$atom][$neighborAtom] <= $minDistFilter)
         {
print STDOUT "Rejecting $ext2CentralItemMap[$neighborAtom]\n";
            $rejectedAtoms[++$numRejected] = $ext2CentralItemMap[$neighborAtom];
         }
      }
   }

   if ($numAtoms%50 != 0)
      {print STDOUT "\n";}

   # Update the number of atoms in the system.
   &setNumAtoms($newNumAtoms);

   # Save these new atom positions in the referenced fractABC and propogate
   #   them to the other representations.  Also, record the necessary naming
   #   data.
   foreach $atom (1..$numAtoms)
   {
      foreach $axis (1..3)
         {$fractABC[$atom][$axis] = $fractABCLocal[$atom][$axis];}
      &getDirectXYZ($atom);
      &getDirectABC($atom);
      $atomElementName[$atom] = $atomElementNameLocal[$atom];
      $atomElementID[$atom] = $atomElementIDLocal[$atom];
      $atomSpeciesID[$atom] = $atomSpeciesIDLocal[$atom];
   }
}

sub checkBoundingBox
{
   # Define the passed parameters.
   my $shiftStyle = $_[0];  # 0=translation; 1=rotation.

   # Define local variables.
   my $atom;
   my $axis;
   my $move;
   my @doMove;

   # Check that each direction of the abc fractional representation is inside
   #   the simulation box.   If it is not, then we will have to move the atom
   #   and re-propogate the position.  The rotation situation is very ugly
   #   because the directXYZ has already been reset to the original non-rotated
   #   position because it must be re-rotated after certain atoms are shifted
   #   to NOT go out of the simulation box when the rotation is applied.  The
   #   effect is that the fractABC used to determine if the atom is outside
   #   the box is NOT the coordinate position that must be shifted by +,-1.
   #   The coordinate position to shift is in the directXYZ now and must be
   #   retrieved through the directABC to fractABC first.  Ugly.

   # Check all the atoms.
   foreach $atom (1..$numAtoms)
   {
      # Assume that the atom will not be moved along any axis.
      @doMove = (0,0,0,0);  # The first index is never used.

      # Identify which axes the atom should be moved along and in which
      #   direction (+ or -).
      foreach $axis (1..3)
      {
         if($fractABC[$atom][$axis] > 1.0)
            {$doMove[$axis] = 1;}
         elsif($fractABC[$atom][$axis] < 0.0)
            {$doMove[$axis] = -1;}
      }

      # For the shiftStyle of 0 we can simply move the atom along each axis
      #   directly.  However, for shiftStyle 1 we are not going to shift this
      #   atom to a periodic cell because the problem was that the atom was
      #   *rotated* outside the box, not shifted outside the box.  Therefor, we
      #   will restore all original fractABCs from a saved directXYZ copy and
      #   **then** shift that original fractABC so that when the rotation is
      #   reapplied it **will not** become shifted outside the box.  Annoying.
      if ($shiftStyle == 1)
      {
         &getDirectABC($atom);
         &getFractABC($atom);
      }

      # Now we can apply the shifts to each axis.
      foreach $axis (1..3)
      {
         if($doMove[$axis] == 1)
            {$fractABC[$atom][$axis] -= 1.0;}
         elsif($doMove[$axis] == -1)
            {$fractABC[$atom][$axis] += 1.0;}
      }

      foreach $move (@doMove)
      {
         if ($move != 0)
         {
            &getDirectXYZ($atom);
            &getDirectABC($atom);
            &getFractABC($atom);
            last;
         }
      }
   }
}


# This subroutine will sort all the atoms according to element name.
#   NOTICE:  It is important that if this subroutine is used that it be
#   called before any calculations be performed such as the creation of the
#   minimal distance matrix or the extended atom positions including the
#   effects of periodic boundary conditions.  This subroutine should only be
#   used immediatly after reading a data set.  Otherwise some data structure
#   could be created from the unsorted data, then the data is sorted and the
#   maping between the data structure and the atom data is lost.
sub sortAtoms
{
   # Define local variables.
   my $atom;
   my $axis;
   my @elementNames;
   my @sortedIndices;  # Index=current sorted atom; value=where it came from.

   # Initialize the list of atomic element names.
   foreach $atom (1..$numAtoms)
      {push (@elementNames,lc($atomElementName[$atom]));}

   # Sort the element names to obtain the sorted indices. For a more detailed
   #   explaination of what is being done here please search for the sorting
   #   work done on cellSizes and cellDims.
   @sortedIndices = sort {$elementNames[$a] cmp $elementNames[$b]}
         0..$#elementNames;

   # Adjust the indices to match the indices that start from 1.
   unshift (@sortedIndices,"");
   foreach $atom (1..$numAtoms)
      {$sortedIndices[$atom]++;}

   # Apply this set of indices to each of the ordered atomic lists that exists.
   &applySort(\@atomElementName,\@sortedIndices);
   &applySort(\@atomElementID,\@sortedIndices);
   &applySort(\@atomSpeciesID,\@sortedIndices);
   &applySort(\@atomicZ,\@sortedIndices);
   &applySort(\@atomTag,\@sortedIndices);
   &applySort(\@datSklMap,\@sortedIndices);
   &applySort(\@sklDatMap,\@sortedIndices);
   &applySort(\@moleculeName,\@sortedIndices);
   &applySort(\@moleculeSeqNum,\@sortedIndices);
   &applySort(\@residueName,\@sortedIndices);
   &applySort(\@residueSeqNum,\@sortedIndices);
   foreach $axis (1..3)
   {
      &applySort(\@{$fractABC[$axis]},\@sortedIndices);
      &applySort(\@{$directABC[$axis]},\@sortedIndices);
      &applySort(\@{$directXYZ[$axis]},\@sortedIndices);
   }
}

sub applySort
{
   # Define passed parameters.
   my $listToSort_ref;
   my $listIndices_ref;

   # Define local variables.
   my @tempList;
   my $numListElements;
   my $listElement;

   # Assigen passed parameters;
   $listToSort_ref = $_[0];
   $listIndices_ref = $_[1];

   # Compute the number of elements in the list to sort.  Keep in mind that the
   #   data to be sorted starts at array index 1 and so array index 0 is empty.
   #   Hence the need for the -1.
   $numListElements = scalar(@$listToSort_ref)-1;

   if ($numListElements > 0)
   {
      # Copy the list to sort to a temporary array.
      foreach $listElement (1..$numListElements)
         {$tempList[$listElement] = $listToSort_ref->[$listElement];}

      # Apply the sorted indices to sort the temporary array back into the
      #   original array.
      foreach $listElement (1..$numListElements)
         {$listToSort_ref->[$listElement] =
                $tempList[$listIndices_ref->[$listElement]];}
   }
}


sub printOLCAO
{
   # Define passed paramters.
   my $fileHandle = $_[0];
   my $title = $_[1];
   my $style = $_[2];

   # Declare local variables.
   my $axis;
   my $atom;
   my @angleDegrees;
   my @atoms;

   # Store the name for each atom.
   foreach $atom (1..$numAtoms)
      {$atoms[$atom] = "$atomElementName[$atom]" . "$atomSpeciesID[$atom]";}

   # Convert the angles into degrees for printing.
   foreach $axis (1..3)
      {$angleDegrees[$axis] = $angle[$axis] * 180/$pi;}

   print $fileHandle "title\n";
   print $fileHandle "$title\n";
   print $fileHandle "@systemTitle\n";
   print $fileHandle "end\n";
   print $fileHandle "cell\n";
   foreach $axis (1..3)
      {printf $fileHandle "%13.8f",$mag[$axis];}
   foreach $axis (1..3)
      {printf $fileHandle "%13.8f",$angleDegrees[$axis];}
   print $fileHandle "\n$style $numAtoms\n";

   # Print the atomic positions in the appropriate style.
   if ($style eq "fractional")
   {
      foreach $atom (1..$numAtoms)
      {
         printf $fileHandle "%6s %12.8f %12.8f %12.8f\n",$atoms[$atom],
               $fractABC[$atom][1],$fractABC[$atom][2],$fractABC[$atom][3];
      }
   }
   else # Cartesian
   {
      foreach $atom (1..$numAtoms)
      {
         printf $fileHandle "%6s %12.8f %12.8f %12.8f\n",$atoms[$atom],
               $directXYZ[$atom][1],$directXYZ[$atom][2],$directXYZ[$atom][3];
      }
   }

   print $fileHandle "space 1_a\n";
   print $fileHandle "supercell 1 1 1\n";
   print $fileHandle "full\n";

   close($fileHandle);
}


sub printOLCAOMap
{
   # Define passed parameters.
   my $fileHandle = $_[0];

   # Define local variables.
   my $atom;

   # Print the header.
   print $fileHandle "OATOM#   MOL#   RES#   RESNAME   ELEMENT   TYPE\n";

   # Print the information for each atom.
   foreach $atom (1..$numAtoms)
   {
      printf $fileHandle "%6d%7d%7d%10s%10s%7s\n",$atom,$moleculeSeqNum[$atom],
            $residueSeqNum[$atom],$residueName[$atom],$atomElementName[$atom],
            $atomTag[$atom];
   }

   close ($fileHandle);
}


sub printVASP
{
   # Define passed parameters.
   my $cellType   = $_[0];
   my $potType    = $_[1];
   my $subPotType = $_[2];
   my $doGamma    = $_[3];
   my $jobType    = $_[4];

   # Define local variables.
   my $element;
   my $currElement;
   my $atom;
   my $axisABC;
   my $axisXYZ;
   my $potDir;
   my $subPotDir;
   my @VASPkp;
   my $ismear;
   my $isif;
   my @magAngstroms;  # Cell lattice parameter magnitudes in angstroms.

   # Compute the magnitude of the lattice parameters in angstroms.
   $magAngstroms[1] = $mag[1]*$bohrRad;
   $magAngstroms[2] = $mag[2]*$bohrRad;
   $magAngstroms[3] = $mag[3]*$bohrRad;

   # Compute the number of atoms of each element in the model.
   &createElementList;
   &countElementAtoms;

   # Open the POSCAR for printing.
   open (POSCAR,">POSCAR") || die "Cannot open POSCAR for writing.\n";


   # BEGIN PRINTING POSCAR


   # Print the POSCAR header using the first line of the systemTitle as the
   #   name for the system.
   print POSCAR " System $systemTitle[0]\n";

   # Print the scaling factor.
   print POSCAR "1.000\n";

   # Print the lattice vectors.
   foreach $axisABC (1..3)
   {
      foreach $axisXYZ (1..3)
         {printf POSCAR "%16.8e",$realLattice[$axisABC][$axisXYZ];}
      print POSCAR "\n";
   }

   # Print the number of atoms of each element.
   print POSCAR "$elementCount[1]";
   foreach $element (2..$numElements)
      {print POSCAR " $elementCount[$element]";}
   print POSCAR "\n";

   # Print the coordinates in fractional units.  NOTE in VASP that "direct" is
   #   used in comparison to Cartesian.  The mapping is as follows:
   #   VASP(direct) = OLCAO(fract); VASP(Cartesian) = OLCAO(direct).  Try not
   #   to get confused.  ><;
   print POSCAR "direct\n";
   foreach $atom (1..$numAtoms)
   {
      foreach $axisABC (1..3)
         {printf POSCAR "%16.8e",$fractABC[$atom][$axisABC];}
      print POSCAR "\n";
   }

   # Close the POSCAR file.
   close (POSCAR);



   # BEGIN PRINTING POTCAR


   # Make sure that we have access to the VASP potential database.
   if (! -d $VASPPOT_DIR)
   {
      print STDOUT "The VASP potential database variable is not defined.\n\n";
      print STDOUT "Please define \$VASPPOT_DIR in the environment and have\n";
      print STDOUT "   it point to the directory that contains directories\n";
      print STDOUT "   for the other potentials.\n";
   }

   # Clear previous data.
   system("rm -f POTCAR");

   # Decompress each potential file in element order and append to the POTCAR.
   foreach $element (1..$numElements)
   {
      # Get the element name and for the case of paw type potentials make the
      #   first character a capital letter.
      $currElement = lc($elementList[$element]);
      if ($potType =~ /paw/)
         {$currElement = ucfirst($currElement);}

      # Start looking for the requested potential file's directory.
      if ($potType eq "potLDA")
         {$potDir = $VASPPOT_USPP_LDA;}
      elsif ($potType eq "potGGA")
         {$potDir = $VASPPOT_USPP_GGA;}
      elsif ($potType eq "potpawLDA")
         {$potDir = $VASPPOT_PAW_LDA;}
      elsif ($potType eq "potpawGGA")
         {$potDir = $VASPPOT_PAW_GGA;}
      elsif ($potType eq "potpawPBE")
         {$potDir = $VASPPOT_PAW_PBE;}
      elsif ($potType eq "potpawLDA5x")
         {$potDir = $VASPPOT_PAW_LDA5x;}
      elsif ($potType eq "potpawPBE5x")
         {$potDir = $VASPPOT_PAW_PBE5x;}

      # Check for the existence of the requested sub potential type directory.
      $subPotDir = $potDir . "/" . $currElement . $subPotType;
      if (-d $subPotDir)
      {
         if (-f "$subPotDir/POTCAR.Z")
            {system("zcat $subPotDir/POTCAR.Z >> POTCAR");}
         elsif (-f "$subPotDir/POTCAR")
            {system("cat $subPotDir/POTCAR >> POTCAR");}
         else
            {die "Cannot find the POTCAR in $subPotDir\n";}
      }
      else
      {
         # First check for the existence of the potential without the sub type.
         #   If it is found, then use it.  If it isn't found, then die.
         $subPotDir = $potDir . "/" . $currElement;
         print STDOUT "Attempting a potential of the $subPotDir type.\n";

         if (-d $subPotDir)
         {
            if (-f "$subPotDir/POTCAR.Z")
               {system("zcat $subPotDir/POTCAR.Z >> POTCAR");}
            elsif (-f "$subPotDir/POTCAR")
               {system("cat $subPotDir/POTCAR >> POTCAR");}
            else
               {die "Cannot find the POTCAR in $subPotDir\n";}
         }
      }
   }


   # Print the KPOINTS file.

   # Assume a default value for the KPOINTS based on the size of the cell.
   if ($doGamma == 1)
      {$VASPkp[1]=1; $VASPkp[2]=1; $VASPkp[3]=1;}
   else
   {
      foreach $axisABC (1..3)
      {
         if ($magAngstroms[$axisABC] > 10)
            {$VASPkp[$axisABC] = 1;}
         elsif ($magAngstroms[$axisABC] > 5)
            {$VASPkp[$axisABC] = 2;}
         else
            {$VASPkp[$axisABC] = 4;}
      }
   }

   # Open the KPOINTS file.
   open (KPOINT,">KPOINTS") || die "Cannot open KPOINTS for writing.\n";

   # Print the KPOINTS file data.
   if ($doGamma == 1)
   {
   print KPOINT <<ENDKPOINTS;
Gamma KPoint only
0
G
1 1 1
0 0 0
ENDKPOINTS
   }
   else
   {
   print KPOINT <<ENDKPOINTS;
Automatic Mesh Generation
0
Monkhorst-Pack
$VASPkp[1] $VASPkp[2] $VASPkp[3]
0 0 0
ENDKPOINTS
   }

   # Close the KPOINTS file.
   close (KPOINT);



   # Print the INCAR file.

   # Determine the value of the $ismear variable.  If there are more than four
   #   kpoints then use the linear analytic tetrahedron method.  Otherwise
   #   just use the Gaussian method.
   if ($VASPkp[1] * $VASPkp[2] * $VASPkp[3] > 4)
      {$ismear = -5;}
   else
      {$ismear = 0;}

   # Modify parameters based on the job type.
   if ($jobType eq "relaxfull")
      {$isif = 3;}
   elsif ($jobType eq "relaxion1")
      {$isif = 2;}
   elsif ($jobType eq "relaxion2")
      {$isif = 4;}
   elsif ($jobType eq "relaxvol")
      {$isif = 7;}
   else
   {
      print STDOUT "This job type has not been implemented yet.\n";
      exit;
   }

   # Open the INCAR file.
   open (INCAR,">INCAR") || die "Cannot open INCAR for writing.\n";

   # Print the INCAR file data.
   print INCAR <<ENDINCAR;
System = $systemTitle[0]

ISMEAR = $ismear

PREC  = Accurate  ! low, medium, normal are other options. Use suitable one.
ENCUT = 600 eV    ! Decide considering the crystal size and accuracy you want.
EDIFF = 1.0E-5    ! Enegy difference coverg. limit for electronic optimization.
EDIFFG = -1.0E-3  ! Enegy difference covergence limit for ionic optimization.

IBRION = 1        ! 0 for MD, 1 best, 2 for diff relax.
NSW    = 100      ! Total number of ionic steps.
ISIF   = $isif        ! 2 and 4 ionic, 7 volume and 3 both.

LREAL  = Auto     ! Proj. on real space. use FALSE (default) for recip space.
NPAR   = 12       ! Best sqrt of NCPUs used. should be >= NCPUs/32.
ALGO   = Fast     ! default is Normal.
LCHARGE = .FALSE. ! Do not print the CHGCAR.
LWAVE   = .FALSE. ! Do not print the WAVCAR.
ENDINCAR

   # Close the INCAR file.
   close (INCAR);

}

sub printPDB
{
   # Define passed parameters.
   my $pdbFile = $_[0];
   my $pdbFormat = $_[1];

   # Define local variables.
   my $atom;

   # Open the PDB file for writing.
   open (PDB,">$pdbFile") || die "Cannot open $pdbFile for writing.\n";

   # Print the HEADER and header-like REMARKS.
   print PDB "HEADER\n";
   print PDB "REMARK   1\n";
   print PDB "REMARK   1 Created by printPDB.\n";
   print PDB "REMARK   1\n";

   # Only print crystal cell information if pdbForm indicates that the format
   #   should include it.
   if ($pdbFormat eq "crystal")
   {
      printf PDB "%6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%10s\n","CRYST1",$mag[1],
         $mag[2],$mag[3],$angleDeg[1],$angleDeg[2],$angleDeg[3],$spaceGroupName;
   }

   # Print each of the atoms. The format is specified here:
   # www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM  
   foreach $atom (1..$numAtoms)
   {
      printf PDB "%6s%5d%3s%-2d%6s%4d%2s%8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s\n",
         "ATOM  ",$atom,uc($atomElementName[$atom]),$atomSpeciesID[$atom],
         " MOL A",1,"  ",$directXYZ[$atom][1],$directXYZ[$atom][2],
         $directXYZ[$atom][3],1.0,0.0,"          ",uc($atomElementName[$atom]);
   }

   # Terminate the ATOM section.
   print PDB "TER\n";

   # Print connectivity information.

   close (PDB);
}

sub printCIF
{
   # Define passed parameters.
   my $currentName             = $_[0];
   my $cellName                = $_[1];
   my $sortedAtomElementID_ref = $_[2];
   my $sortedAtomSpeciesID_ref = $_[3];
   my $sortedAtomTypeID_ref    = $_[4];
   my $sortedFractABC_ref      = $_[5];

   # Define local variables.
   my $i;
   my $element;
   my $species;
   my $type;
   my $fileName;
   my @magAngstroms;  # Cell lattice parameter magnitudes in angstroms.

   # Compute the magnitude of the lattice parameters in angstroms.
   $magAngstroms[1] = $mag[1]*$bohrRad;
   $magAngstroms[2] = $mag[2]*$bohrRad;
   $magAngstroms[3] = $mag[3]*$bohrRad;

   # Open the CIF file for writing.
   $fileName = $currentName . ".cif";
   open (CIF,">$fileName") || die "Cannot open $fileName for writing\n";

   #Prepare the header for the file.
   print CIF "data_$currentName\n";
   print CIF "_symmetry_cell_setting $cellName\n";
   print CIF "_symmetry_space_group_name_H-M 'P 1'\n";
   print CIF "loop_\n";
   print CIF "   _symmetry_equiv_pos_as_xyz\n";
   print CIF "              X,Y,Z\n";

   #Print the crystal lattice information.
   print CIF "_cell_length_a $magAngstroms[1]\n";
   print CIF "_cell_length_b $magAngstroms[2]\n";
   print CIF "_cell_length_c $magAngstroms[3]\n";
   print CIF "_cell_angle_alpha $angleDeg[1]\n";
   print CIF "_cell_angle_beta $angleDeg[2]\n";
   print CIF "_cell_angle_gamma $angleDeg[3]\n";


   #Prepare to print the atoms
   print CIF "loop_\n";
   print CIF "   _atom_site_label\n";
   print CIF "   _atom_site_type_symbol\n";
   print CIF "   _atom_site_fract_x\n";
   print CIF "   _atom_site_fract_y\n";
   print CIF "   _atom_site_fract_z\n";

   #Loop to print all the atoms in the cell.
   for ($i=1;$i<=$numAtoms;$i++)
   {
      # Determine the element name and species ID number.
      $element = $elementList[$sortedAtomElementID_ref->[$i]];
      $species = $sortedAtomSpeciesID_ref->[$i];
      $type    = $sortedAtomTypeID_ref->[$i];

      printf CIF "   %2s%-s%1s%-3s %2s %12.6f %12.6f %12.6f\n",$element,
            $species,"_",$type,ucfirst($element),
            $sortedFractABC_ref->[$i][1],
            $sortedFractABC_ref->[$i][2],$sortedFractABC_ref->[$i][3];
   }

   #Close it.
   close (CIF);
}


sub printLMP
{
   # Define short names for passed parameters.
   my $lmpFile = $_[0];

   # Define local variables.
   my $atom;
   my $element;
   my $species;
   my $xy;
   my $xz;
   my $yz;
   my $numUniqueSpecies;
   my $currentElement;
   my $currentSpecies;
   my $atomicMasses_ref;
   $atomicMasses_ref = ElementData::getAtomicMassesRef;

   # Compute the number of unique species.
   $numUniqueSpecies = 0;
   foreach $element (1..$numElements)
   {
      foreach $species (1..$numSpecies[$element])
         {$numUniqueSpecies++;}
   }

   # Open the lammps file for writing.
   open (LMP, ">$lmpFile") || die "Cannot open $lmpFile for writing.\n";

   # Print the header.
   print LMP "$lmpFile\n\n";

   # Print the number of items and number of types.
   print LMP "$numAtoms atoms\n";
   print LMP "$numUniqueSpecies atom types\n";

   # Print the orthogonal cell information.
   print LMP "0.000 $mag[1] xlo xhi\n";
   print LMP "0.000 $mag[2] ylo yhi\n";
   print LMP "0.000 $mag[3] zlo zhi\n";

   # If the cell is non-orthogonal, then compute and print the tilt parameters.
   #   See: https://lammps.sandia.gov/doc/Howto_triclinic.html
   if (($angleDeg[1] != 90) || ($angleDeg[2] != 90) || ($angleDeg[3] != 90))
   {
      $xy = $mag[2] * cos($angle[3]);
      $xz = $mag[3] * cos($angle[2]);
      $yz = ($mag[2]*$mag[3]*cos($angle[1]) - $xy*$xz)
            / sqrt($mag[2]*$mag[2] - $xy*$xy);
      print LMP "$xy $xz $yz xy xz yz\n\n";
   }

   # Print information about the masses.
   print LMP "Masses\n\n";
   $currentElement = 0;
   $currentSpecies = 0;
   $numUniqueSpecies = 0;
   foreach $atom (1..$numAtoms)
   {
      if ($atomElementID[$atom] != $currentElement)
      {
         $currentElement = $atomElementID[$atom];
         $currentSpecies = $atomSpeciesID[$atom];
         $numUniqueSpecies++;
         print LMP "$numUniqueSpecies $atomicMasses_ref->[$atomicZ[$atom]] ";
         print LMP "# $atomElementName[$atom] $currentSpecies\n";
      }
      elsif ($atomSpeciesID[$atom] != $currentSpecies)
      {
         $currentSpecies = $atomSpeciesID[$atom];
         $numUniqueSpecies++;
         print LMP "$numUniqueSpecies $atomicMasses_ref->[$atomicZ[$atom]] ";
         print LMP "# $atomElementName[$atom] $currentSpecies\n";
      }

      $atomUniqSpecies[$atom] = $currentSpecies;
   }

   # Print information about the atoms.
   print LMP "\nAtoms\n\n";
   foreach $atom (1..$numAtoms)
   {
      print LMP "$atom $atomUniqSpecies[$atom] ";
      print LMP "$directXYZ[$atom][1] $directXYZ[$atom][2] ";
      print LMP "$directXYZ[$atom][3] ";
      print LMP "# $atomElementName[$atom] $currentSpecies\n";
   }

   # Close the lammps file.
   close (LMP);
}



# Get the number of atoms of each element in the system.
sub countElementAtoms
{
   # Define local variables.
   my $atom;
   my $element;

   # Initialize the element count.
   foreach $element (1..$numElements)
      {$elementCount[$element] = 0;}

   foreach $atom (1..$numAtoms)
      {$elementCount[$atomElementID[$atom]]++;}
}


# Map the element name to the element number (Z) for each atom in the system.
#   This subroutine assums that the link to the element data base has been made.
sub mapElementNumber
{
   # Declare local variables.
   my $atom;
   my $element;

   foreach $atom (1..$numAtoms)
      {$atomicZ[$atom] = &getElementNumber($atomElementName[$atom]);}
}


# Return the element number (Z) of a particular atom.
#   This subroutine assums that the link to the element data base has been made.
sub getElementNumber
{
   # Define passed parameters.
   my $givenElement = $_[0];

   # Declare local variables.
   my $element;

   foreach $element (1..$#$elementNames_ref)
   {
      if ($givenElement eq $elementNames_ref->[$element])
         {return $element;last;}
   }
}

sub assignCovalRadii
{
   # Declare local variables.
   my $atom;

   foreach $atom (1..$numAtoms)
      {$covalRadii[$atom] = $covalRadii_ref->[$atomicZ[$atom]];}
}


sub createMinDistMatrix
{
   # Define passed parameters.
   my $numItems1     = $_[0];
   my $numItems2     = $_[1];
   my $items1ABC_ref = $_[2];
   my $items2ABC_ref = $_[3];
   my $items1XYZ_ref = $_[4];
   my $items2XYZ_ref = $_[5];

   # Define local variables.
   my $item;
   my $numItemsTotal;
   my $numItemsTotalExt;
   my @itemsTotalABC;
   my @itemsTotalXYZ;
   my @itemsTotalExt;
   my @itemsTotalExtXYZ;
   my @itemsTotalExtXYZList;
   my @itemsTotalExtABCList;

   # Join the two lists.
   $numItemsTotal = 0;
   foreach $item (1..$numItems1)
   {
      $itemsTotalABC[++$numItemsTotal] = $items1ABC_ref->[$item];
      $itemsTotalXYZ[$numItemsTotal]   = $items1XYZ_ref->[$item];
   }
   foreach $item (1..$numItems2)
   {
      $itemsTotalABC[++$numItemsTotal] = $items2ABC_ref->[$item];
      $itemsTotalXYZ[$numItemsTotal]   = $items2XYZ_ref->[$item];
   }

   # Determine the number of cells in each direction needed to account for all
   #   intercell interactions.
   &computeNumCells($numItemsTotal,\@itemsTotalABC,\@negBit,\@posBit);

   # Compute the position of each item in each of the non-origin cells where
   #   it is needed because it interacts with an item in the origin cell.
   $numItemsTotalExt = &computeExtendedPos(\@negBit,\@posBit,$numItemsTotal,
         $numItemsTotalExt,\@itemsTotalABC,\@central2ExtItemMap,
         \@ext2CentralItemMap,\@itemsTotalExtXYZ,\@itemsTotalExtXYZList,
         \@itemsTotalExtABCList);

   # Compute and record the interatomic distances between all items.  Flag #1
   #   is used to record the minDist matrix from the general purpose
   #   "obtainAtomicInteraction" subroutine.
   &obtainAtomicInteraction(1,$numItemsTotal,$numItemsTotalExt,
         \@itemsTotalXYZ,\@itemsTotalExtXYZList);
}


sub createBondingList
{
   # Determine the number of cells in each direction needed to account for all
   #   intercell interactions.
   &computeNumCells($numAtoms,\@fractABC,\@negBit,\@posBit);

   # Compute the position of each atom in each of the non-origin cells where
   #   it is needed because it interacts with an atom in the origin cell.
   $numAtomsExt = &computeExtendedPos(\@negBit,\@posBit,$numAtoms,
         $numAtomsExt,\@fractABC,\@central2ExtItemMap,\@ext2CentralItemMap,
         \@extDirectXYZ,\@extDirectXYZList,\@extFractABCList);

   # Assign the covalent radius of each atom from the database.
   &assignCovalRadii;

   # Compute and record which atoms are bonded to which other atoms.  Flag #2
   #   is used to record extended bonding information from the general purpose
   #   "obtainAtomicInteraction" subroutine.
   &obtainAtomicInteraction(2,$numAtoms,$numAtomsExt,\@directXYZ,
         \@extDirectXYZList);

   # Map the extended atom positions back to the original central cell.
   &mapExtToCentral;
}

sub computePoreMap
{
   # Define passed parameters.
   my $resolution_ref = $_[0]; # Step size in angstroms x,y,z (not abc yet).
   my $numMeshPoints_ref = $_[1]; # In each x,y,z direction (not abc yet).

   # Declare local variables.
   my @indexNum;
   my $axis;
   my $atom;
   my $xPoint;
   my $yPoint;
   my $zPoint;
   my $currCovalRadiusSqrd;
   my @miniCubeSize; # Circumscribes a sphere of size equal to the coval rad.
#   my @coordDiff; # Used by the Perl implementation of CdistanceSqrd
   my $distanceSqrd;

   # Determine the number of cells in each direction needed to account for all
   #   intercell interactions.
   &computeNumCells($numAtoms,\@fractABC,\@negBit,\@posBit);

   # Compute the position of each atom in each of the non-origin cells where
   #   it is needed because it interacts with an atom in the origin cell.
   $numAtomsExt = &computeExtendedPos(\@negBit,\@posBit,$numAtoms,
         $numAtomsExt,\@fractABC,\@central2ExtItemMap,\@ext2CentralItemMap,
         \@extDirectXYZ,\@extDirectXYZList,\@extFractABCList);

   # For each atom, visit its coordinate and fill in the pore map that
   #   surrounds it within a scaled distance of the colavent radius with "1"s.
   foreach $atom (1..$numAtomsExt)
   {
      # Get the already scaled covalent radius for this atom.
      $currCovalRadiusSqrd =
         ($covalRadii_ref->[$atomicZ[$ext2CentralItemMap[$atom]]])**2;

      # Find the closest positional index numbers in the pore map to the atom
      #   site.
      foreach $axis (1..3) 
      {
         # Round the ratio of the position in angstroms over the resolution.
         $indexNum[$axis]
            = int($extDirectXYZList[$atom][$axis]/$resolution_ref->[$axis]);

         # Compute the exact size of the mini cube for this atom.
         $miniCubeSize[$axis] = sqrt($currCovalRadiusSqrd)
               / $resolution_ref->[$axis];

         # Round the minicube size to find grid positions.
         if (($miniCubeSize[$axis] - int($miniCubeSize[$axis])) > 0.5)
            {$miniCubeSize[$axis] = int($miniCubeSize[$axis]) + 1;}
         else
            {$miniCubeSize[$axis] = int($miniCubeSize[$axis]);}
      }

      # Triple loop over x,y,z of minicube dimensions without indentation.
      foreach $xPoint (($indexNum[1]-$miniCubeSize[1])
                     ..($indexNum[1]+$miniCubeSize[1]))
      {
         # Only consider mesh points inside the central system cell.
         if (($xPoint < 1) or ($xPoint > $numMeshPoints_ref->[1]))
            {next;}
      foreach $yPoint (($indexNum[2]-$miniCubeSize[2])
                     ..($indexNum[2]+$miniCubeSize[2]))
      {
         # Only consider mesh points inside the central system cell.
         if (($yPoint < 1) or ($yPoint > $numMeshPoints_ref->[2]))
            {next;}
      foreach $zPoint (($indexNum[3]-$miniCubeSize[3])
                     ..($indexNum[3]+$miniCubeSize[3]))
      {
         # Only consider mesh points inside the central system cell.
         if (($zPoint < 1) or ($zPoint > $numMeshPoints_ref->[3]))
            {next;}

         # Compute the square of the exact distance between this mesh
         #   point and the atom.
# This is the Inline C implementation.
         $distanceSqrd = CdistanceSqrd($limitDistSqrd,$limitDist,
            $extDirectXYZList[$atom][1],$extDirectXYZList[$atom][2],
            $extDirectXYZList[$atom][3],$xPoint*$resolution_ref->[1],
            $yPoint*$resolution_ref->[2],$zPoint*$resolution_ref->[3]);

# Perl implementation of CdistanceSqrd         
#         foreach $axis (1..3)
#         {
#            $coordDiff[$axis] = $xPoint*$resolution_ref->[$axis]
#                              - $extDirectXYZList[$atom][$axis];
#            if ($coordDiff[$axis] > $limitDist)
#               {$distanceSqrd = $bigReal; last;}
#         }
#
#         # At this point the limit dist criterion has been met so we can
#         #   compute the square of the distance.
#         $distanceSqrd = $coordDiff[1]*$coordDiff[1]
#                       + $coordDiff[2]*$coordDiff[2]
#                       + $coordDiff[3]*$coordDiff[3];

         if ($distanceSqrd < $currCovalRadiusSqrd)
            {$extPoreMap[$xPoint][$yPoint][$zPoint] = 1;}
      }
      }
      }
   }
}

sub createExtBondingList
{
   # Determine the number of cells in each direction needed to account for all
   #   intercell interactions.
   &computeNumCells($numAtoms,\@fractABC,\@negBit,\@posBit);

   # Compute the position of each atom in each of the non-origin cells where
   #   it is needed because it interacts with an atom in the origin cell.
   $numAtomsExt = &computeExtendedPos(\@negBit,\@posBit,$numAtoms,
         $numAtomsExt,\@fractABC,\@central2ExtItemMap,\@ext2CentralItemMap,
         \@extDirectXYZ,\@extDirectXYZList,\@extFractABCList);

   # Assign the covalent radius of each atom from the database.
   &assignCovalRadii;

   # Compute and record which atoms are bonded to which other atoms.  Flag #2
   #   is used to record extended bonding information from the general purpose
   #   "obtainAtomicInteraction" subroutine.
   &obtainAtomicInteraction(2,$numAtoms,$numAtomsExt,\@directXYZ,
         \@extDirectXYZList);
}


sub computeRPDF
{
   # Determine the number of cells in each direction needed to account for all
   #   intercell interactions.
   &computeNumCells($numAtoms,\@fractABC,\@negBit,\@posBit);

   # Compute the position of each atom in each of the non-origin cells where
   #   it is needed because it interacts with an atom in the origin cell.
   $numAtomsExt = &computeExtendedPos(\@negBit,\@posBit,$numAtoms,
         $numAtomsExt,\@fractABC,\@central2ExtItemMap,\@ext2CentralItemMap,
         \@extDirectXYZ,\@extDirectXYZList,\@extFractABCList);

   # Compute and record the RPDF.  Flag #3 is used to specify RPDF for this
   #   general subroutine.
   &obtainAtomicInteraction(3,$numAtoms,$numAtomsExt,\@directXYZ,
         \@extDirectXYZList);
}

sub createQList
{
   # Determine the number of cells in each direction needed to account for all
   #   intercell interactions.
   &computeNumCells($numAtoms,\@fractABC,\@negBit,\@posBit);

   # Compute the position of each atom in each of the non-origin cells where
   #   it is needed because it interacts with an atom in the origin cell.
   $numAtomsExt = &computeExtendedPos(\@negBit,\@posBit,$numAtoms,
         $numAtomsExt,\@fractABC,\@central2ExtItemMap,\@ext2CentralItemMap,
         \@extDirectXYZ,\@extDirectXYZList,\@extFractABCList);

   # Assign the covalent radius of each atom from the database.
   &assignCovalRadii;

   # Compute and record the BOO.  Flag #4 is used to record BOO specific 
   #   information for this general purpose subroutine.
   &obtainAtomicInteraction(4,$numAtoms,$numAtomsExt,\@directXYZ,
         \@extDirectXYZList);
}

sub computeAtomMeshDist
{
   # Determine the number of cells in each direction needed to account for all
   #   intercell interactions.
   &computeNumCells($numAtoms,\@fractABC,\@negBit,\@posBit);

   # Compute the position of each atom in each of the non-origin cells where
   #   it is needed because it interacts with an atom in the origin cell.
   $numAtomsExt = &computeExtendedPos(\@negBit,\@posBit,$numAtoms,
         $numAtomsExt,\@fractABC,\@central2ExtItemMap,\@ext2CentralItemMap,
         \@extDirectXYZ,\@extDirectXYZList,\@extFractABCList);

   # Compute and record the distance between each atom (of the requested
   #   element) and the points of the defined scan line or mesh.  The #5 is the
   #   flag requesting this type of action from the general purpose subroutine.
   &obtainAtomicInteraction(5,$numScanPoints,$numAtomsExt,\@scanPoints,
         \@extDirectXYZList);
}

sub computeBondMeshDist
{
   # Determine the number of cells in each direction needed to account for all
   #   intercell interactions.
   &computeNumCells($numBondsTotal,\@fractABC,\@negBit,\@posBit);

   # Compute the position of each atom in each of the non-origin cells where
   #   it is needed because it interacts with an atom in the origin cell.
   $numAtomsExt = &computeExtendedPos(\@negBit,\@posBit,$numAtoms,
         $numAtomsExt,\@fractABC,\@central2ExtItemMap,\@ext2CentralItemMap,
         \@extDirectXYZ,\@extDirectXYZList,\@extFractABCList);

   # Compute and record to distance between each atom (of the requested
   #   element) and the points of the defined scan line or mesh.  The #5 is the
   #   flag requesting this type of action from the general purpose subroutine.
   &obtainAtomicInteraction(5,$numScanPoints,$numAtomsExt,\@scanPoints,
         \@extDirectXYZList);
}


sub computeNumCells
{
   # Define the passed paramters.
   my $numItems      = $_[0];  # Could be atoms or bonds or something else.
   my $fractItemABC  = $_[1];  # Could be atoms or bonds or something else.
   my $negBit        = $_[2];
   my $posBit        = $_[3];

   # Declare local variables
   my @currentSineAngle;
   my @extLimitDist;
   my $item;  # Could be atoms or bonds or something else.
   my $axis;

   # Get the minimum of the sine function for angles involving a, then b,
   #   and then c lattice vectors.  This is needed in the next step.
   $currentSineAngle[1] = &min($sineAngle[2],$sineAngle[3]); # a,c  a,b
   $currentSineAngle[2] = &min($sineAngle[1],$sineAngle[3]); # b,c  a,b
   $currentSineAngle[3] = &min($sineAngle[2],$sineAngle[1]); # a,c  b,c

   # Compute the maximum modification to the limit distance for each axis.
   foreach $axis (1..3)
      {$extLimitDist[$axis]=$limitDist/$mag[$axis]/$currentSineAngle[$axis];}

   # To make a significant efficiency gain we will now determine which periodic
   #   cells each item needs to be searched for in.  (Consider an item at the
   #   origin of a large cell.  Other items in the origin cell will seach for
   #   this item in the origin cell (x0), at the +x1 cell, but not at the -x1
   #   cell.  Items near the origin item are still too far from the replication
   #   of the origin item in the -x1 cell to need to count it.  However, items
   #   in the origin cell, but near the +x1 border will want to include the
   #   origin item in the +x1 replicated cell when they search for interacting
   #   items.  This all means that we need to know the position of the origin
   #   item in the origin cell, and in the +x1 cell, but not in the -x1 cell.)
   foreach $item (1..$numItems)
   {
      # We want to eliminate as much work as possible.  Given the position of
      #   an item and the value of the limit distance we can determine the
      #   number of periodic cells that should be considered in each direction.
      #   The difficult aspect is to make this work for arbitrary lattice
      #   parameters.  Consider a hexagonal cell and an item near one border.
      #   We cannot just add the limit distance in the a direction to see if
      #   this atom needs to be searched for in the -a1 cell.  (Recall that
      #   this test is counter intuitive.  If the original position + the limit
      #   crosses into the +a1 cell then this item must be searched for in the
      #   -a1 cell.)  However, the a,b,c position of the item + limitDist in
      #   the a direction may not penetrate into the next cell, *but*, the
      #   a,b,c position of the item + limitDist in the direction perpendicular
      #   to b might penetrate into the next cell.  Back to the hexagonal cell
      #   example:  b-----    (c is out of the screen).
      #              \   .\
      #               \    \
      #                c----a
      # The distance to the cell wall in the a direction is greater than the
      #   distance to the cell wall in the direction perpendicular to the wall.
      #   The difference is small,  but if the limit distance is just the right
      #   length, then this item would not be searched for when it should be.
      #   To correct for this problem we increase the limit distance by a
      #   factor of 1/sin(angle) where the angle is alpha, beta, or gamma as
      #   needed when we add it to the a,b,c item positions.
      #
      # For each bit direction we have to consider the effects of the two
      #   angles.  Example:  The a direction must consider alpha and gamma
      #   because those are the angles between the a,c vectors and a,b vectors
      #   respectively.  The maximum of either case is chosen.

      # Check if this position is outside the defined bounding box.
      if (&itemOutOfBounds($item,$borderCoords_ref))
      {
         foreach $axis (1..3)
         {
            $negBit->[$item][$axis] = 0;
            $posBit->[$item][$axis] = 0;
         }

         # Go to the next item in the main loop.
         next;
      }

      # Check if this position is one of the requested tags.  Always true if
      #   no specific tags were requested.
      if (&itemNotRequested($item,\@atomElementName,\@rpdfSelectAtoms))
      {
         foreach $axis (1..3)
         {
            $negBit->[$item][$axis] = 0;
            $posBit->[$item][$axis] = 0;
         }

         # Go to the next item in the main loop.
         next;
      }


      # Compute the number of cells in each direction (a,b,c).
      foreach $axis (1..3)
      {
         $negBit->[$item][$axis] =
               floor($fractItemABC->[$item][$axis] + $extLimitDist[$axis]);

         $posBit->[$item][$axis] =
               floor(1.0-$fractItemABC->[$item][$axis] + $extLimitDist[$axis]);
      }

      # Each bit indicates the number of periodic cells in that direction that
      #   need to be searched in *by* other items *for* item $item.
   }
}



sub computeExtendedPos
{
   # Define passed parameters.
   my $negBit               = $_[0];
   my $posBit               = $_[1];
   my $numItems             = $_[2];
   my $numItemsExt          = $_[3];  # Explicit return value.
   my $fractItemABC         = $_[4];
   my $central2ExtItemMap   = $_[5];
   my $ext2CentralItemMap   = $_[6];
   my $extDirectXYZItem     = $_[7];
   my $extDirectXYZListItem = $_[8];
   my $extFractABCListItem  = $_[9];

   # Define local variables.
   my @tempFractABC;
   my @tempXYZ;
   my $cellBitA; # Perl will not let me use @cellBit[1..3] here.
   my $cellBitB;
   my $cellBitC;
   my $item;
   my $axis;

   # Initialize the number of extended item positions.
   $numItemsExt = 0;


   # Loop through all the necessary replicated cells for each item and store
   #   the direct space x,y,z coordinates of the item in that cell.
   foreach $item (1..$numItems)
   {
      foreach $cellBitA (-$negBit->[$item][1]..$posBit->[$item][1])
      {
      foreach $cellBitB (-$negBit->[$item][2]..$posBit->[$item][2])
      {
      foreach $cellBitC (-$negBit->[$item][3]..$posBit->[$item][3])
      {
         # Obtain the fractional coordinates as affected by the cell shifting.
         $tempFractABC[1] = $fractItemABC->[$item][1] + $cellBitA;
         $tempFractABC[2] = $fractItemABC->[$item][2] + $cellBitB;
         $tempFractABC[3] = $fractItemABC->[$item][3] + $cellBitC;

         # Obtain the x,y,z coordinates of this item.
         $tempXYZ[1] = $tempFractABC[1]*$realLattice[1][1] +
                       $tempFractABC[2]*$realLattice[2][1] +
                       $tempFractABC[3]*$realLattice[3][1];
         $tempXYZ[2] = $tempFractABC[1]*$realLattice[1][2] +
                       $tempFractABC[2]*$realLattice[2][2] +
                       $tempFractABC[3]*$realLattice[3][2];
         $tempXYZ[3] = $tempFractABC[1]*$realLattice[1][3] +
                       $tempFractABC[2]*$realLattice[2][3] +
                       $tempFractABC[3]*$realLattice[3][3];

         # Compute the x, y, z coordinates of the tempFract a, b, c position
         #   coordinates.
         $extDirectXYZItem->[1][$cellBitA+$negBit->[$item][1]]
                               [$cellBitB+$negBit->[$item][2]]
                               [$cellBitC+$negBit->[$item][3]][$item] =
                               $tempXYZ[1];
         $extDirectXYZItem->[2][$cellBitA+$negBit->[$item][1]]
                               [$cellBitB+$negBit->[$item][2]]
                               [$cellBitC+$negBit->[$item][3]][$item] =
                               $tempXYZ[2];
         $extDirectXYZItem->[3][$cellBitA+$negBit->[$item][1]]
                               [$cellBitB+$negBit->[$item][2]]
                               [$cellBitC+$negBit->[$item][3]][$item] =
                               $tempXYZ[3];

         # Add this item to the extended (flat) list.
         $numItemsExt++;
         $extDirectXYZListItem->[$numItemsExt][1] = $tempXYZ[1];
         $extDirectXYZListItem->[$numItemsExt][2] = $tempXYZ[2];
         $extDirectXYZListItem->[$numItemsExt][3] = $tempXYZ[3];

         $extFractABCListItem->[$numItemsExt][1] = $tempFractABC[1];
         $extFractABCListItem->[$numItemsExt][2] = $tempFractABC[2];
         $extFractABCListItem->[$numItemsExt][3] = $tempFractABC[3];

         # Create the mapping between the item numbers of the items in the
         #   original central cell + item numbers of the SAME CENTRAL CELL ATOMS
         #   as listed in the extended periodic system list.  (E.g. item#1 in
         #   the original central cell is the first item of the central cell
         #   duh.  However, item#1 of the extended list is actually a periodic
         #   image of the first central cell item in the negBit direction.  We
         #   have to keep track of this.  The index number is the original item
         #   number, and the value is the item number in the periodic system.)
         if (($cellBitA == 0) && ($cellBitB == 0) && ($cellBitC == 0))
            {$central2ExtItemMap->[$item] = $numItemsExt;}

         # Create a mapping between the item numbers of the items in the
         #   extended cell to their original item number in the central cell
         #   as listed in the non-extended system list.  (E.g. item#1 in the
         #   extended system is some replication of a central cell atom, this
         #   records which central cell atom.  The index number is the
         #   replicated extended item number and the value is the original
         #   central cell item number.
         $ext2CentralItemMap->[$numItemsExt] = $item;
      }
      }
      }
   }

   return $numItemsExt;
}


sub mapExtToCentral
{
   # Define local variables.
   my $atom;
   my $bond;

   foreach $atom (1..$numAtoms)
   {
      foreach $bond (1..$numBonds[$atom])
         {$bonded[$atom][$bond]=$ext2CentralItemMap[$bondedExt[$atom][$bond]];}
   }
}



# The purpose of this subroutine is to compute one of a number of possible
#   types of interactions between a set of points within the original cell
#   (possibly the atom positions) and another set of points both within
#   the original cell and in neighboring replicated cells.
# Case #1 computes a minimal distance matrix.  This is the minimum distance
#   between item pairs including the effects of periodic boundary conditions
#   At least one item from each pair must be from the central cell.
# Case #2 computes a list of which items are "linked" to which other items, and
#   how many links each item has.  The index numbers for items that are in the
#   central cell are used to index the links between items.  (For example, if
#   bondedExt[1][4] == 77 it means that the fourth link (bond) for item #1 in
#   the central cell list is linked (bonded) to item #77 from the extended cell
#   list of items.)  This is often used to identify which pairs of atoms are
#   bonded.  This data can be reduced such that all the indices refer to the
#   central cell.  Sometimes one listing is better than another for sharing
#   data with other programs or people.
# Case #3 computes the radial pair distribution function for the given set of
#   items with the accompanying tag and space restrictions.
# Case #4 is used for the bond orientational order calculation.  This is not
#   usable at the moment.
# Case #5 computes distances between points in a scan line and atomic positions.
sub obtainAtomicInteraction
{
   # Define passed parameters.
   my $interactionType    = $_[0];
   my $numItems1          = $_[1];
   my $numItems2          = $_[2];
   my $item1XYZ_ref       = $_[3];
   my $item2XYZ_ref       = $_[4];

   # Declare local variables.
   my $cellBitA;
   my $cellBitB;
   my $cellBitC;
   my @diff;
   my $distance;
   my $item1;
   my $item2;
   my $axis;
   my $extAtom2;
   my $item2Start;
   my $theta;
   my $phi;

   # Initialize values for the needed interaction
   &initializeInteractionData($interactionType,$numItems1,$numItems2);

   foreach $item1 (1..$numItems1)
   {
      # Print a running count of progress on the screen so the user does not
      #   get too bored waiting.  This process can take a while for large
      #   systems.
      if ($item1%10 == 0)
         {print STDOUT "|";}
      else
         {print STDOUT ".";}

      if ($item1%50 == 0)
         {print STDOUT " $item1\n";}


      # Check if this item should be skipped.
      if ($skipItem1[$item1])
         {next;}


      foreach $item2 (1..$numItems2)
      {
         # Check if this item should be skipped.
         if ($skipItem2[$item2])
            {next;}

         # Check if this item pair should be skipped.
         if ($skipPair[$item1][$item2])
            {next;}

         # Compute the difference in position between these two items along
         #   each of the x,y,z axes.
         foreach $axis (1..3)
            {$diff[$axis]=$item1XYZ_ref->[$item1][$axis] -
                          $item2XYZ_ref->[$item2][$axis];}

         # Calculate the direct distance between the two items.
# This approach uses C Inline.
         $distance = Cdistance($limitDistSqrd,$limitDist,
               $item1XYZ_ref->[$item1][1],$item1XYZ_ref->[$item1][2],
               $item1XYZ_ref->[$item1][3],$item2XYZ_ref->[$item2][1],
               $item2XYZ_ref->[$item2][2],$item2XYZ_ref->[$item2][3]);

# Perl code version of the Inline C.
#         # Check that no axis is beyond the limit distance. If one is, then
#         #   abandon this distance calculation.
#         foreach $axis (1..3)
#         {
#            if ($diff[$axis] > $limitDist)
#               {$distance = $bigReal; next;}
#         }
#
#         # If we are within the limit distance, then compute the square root of
#         #    the distance. Otherwise, abort.
#         $distance = $diff[1]*$diff[1] + $diff[2]*$diff[2] + $diff[3]*$diff[3];
#         if ($distance < $limitDistSqrd)
#            {$distance = sqrt($distance);}
#         else
#            {$distance = $bigReal; next;}

# This check is used when Cdistance is called.
         # In the event that the distance between the two items is judged to
         #   be beyond a set "limitDist" then the Cdistance subroutine will
         #   return a quantity equal to $bigReal. Thus, this interaction datum
         #   should not be recorded.
         if ($distance == $bigReal)
            {next;}

         # Record the requested information based on these relative
         #   spherical coordinates, item numbers, and interaction type.
         &recordInteractionData($interactionType,$distance,\@diff,
               $item1,$item2);
      }
   }

   # Finalize the collected data according to the interaction type.
   &finalizeInteractionData($interactionType,$numItems1,$numItems2);

   if ($numItems1%50 != 0)
      {print STDOUT "\n";}
}


sub initializeInteractionData
{
   # Define passed parameters.
   my $interactionType = $_[0];
   my $numItems1       = $_[1];
   my $numItems2       = $_[2];

   # Declare local variables
   my $axis;
   my $item;
   my $item1;
   my $item2;
   my $point;
   my $Ylm_m;      # Used for bond orientational order (boo).


   if ($interactionType == 1)  # Minimal distance matrix.
   {
      # This will find the minimal distances between central cell item pairs
      #   including the effects of periodicity.  Watch carefully for the use of
      #   ext2CentralItemMap[$item] and just $item!  We ust search all ext
      #   items but record the results for central cell items only.

      # Reset the minimal distance matrix and self min distance matrix.
      undef @selfMinDist;
      undef @minDist;

      # Don't skip any atoms but skip pairs where item2 is < item1.
      foreach $item1 (1..$numItems1)
      {
         $skipItem1[$item1] = 0;
         foreach $item2 (1..$numItems2)
         {
            if ($ext2CentralItemMap[$item2] < $item1)
               {$skipPair[$item1][$item2] = 1;}
            else
               {$skipPair[$item1][$item2] = 0;}
         }
      }
      foreach $item2 (1..$numItems2)
         {$skipItem2[$item2] = 0;}

      # Initialize the self minimal distance to a very large number.
      foreach $item (1..$numItems1)
         {$selfMinDist[$item] = $bigReal;}

      # Initialize the minimal distances between central cell items to very
      #   large numbers.
      foreach $item1 (1..$numItems1)
      {
         foreach $item2 ($item1..$numItems1)
            {$minDist[$item1][$item2] = $bigReal;}
      }
   }
   elsif ($interactionType == 2) # check for extended bonding
   {
      # Don't skip any atoms or pairs of atoms.
      foreach $item1 (1..$numItems1)
      {
         $skipItem1[$item1] = 0;
         foreach $item2 (1..$numItems2)
            {$skipPair[$item1][$item2] = 0;}
      }
      foreach $item2 (1..$numItems2)
         {$skipItem2[$item2] = 0;}

      # Initialize the count of the number of bonds for each atom in the
      #   central cell.
      foreach $item (1..$numItems1)
         {$numBonds[$item] = 0;}
   }
   elsif ($interactionType == 3) # rpdf
   {
      # Skip certain atoms if they are outside a defined bounding box.
      foreach $item1 (1..$numItems1)
         {$skipItem1[$item1] = &itemOutOfBounds($item1,$borderCoords_ref);}
      foreach $item2 (1..$numItems2)
         {$skipItem2[$item2] = &itemOutOfBounds($item2,$borderCoordsExt_ref);}

      # Skip atom pairs if they are not of the requested elemental types.
      foreach $item1 (1..$numItems1)
      {
         foreach $item2 (1..$numItems2)
         {
            $skipPair[$item1][$item2] =
                  &pairNotOfElement($item1,$ext2CentralItemMap[$item2]);
         }
      }

      # Initialize the values of the rpdf to zero so that it can be accumulated.
      foreach $point (1..1000)
         {$rpdf[$point] = 0.0;}
   }
   elsif ($interactionType == 4) # bond orientational order
   {
      # Don't skip any atoms or pairs.
      foreach $item1 (1..$numItems1)
      {
         $skipItem1[$item1] = 0;
         foreach $item2 (1..$numItems2)
            {$skipPair[$item1][$item2] = 0;}
      }
      foreach $item2 (1..$numItems2)
         {$skipItem2[$item2] = 0;}

      # Initialize the count of the number of bonds for each atom.
      foreach $item (1..$numItems1)
         {$numBonds[$item] = 0;}

      # Initialize the qBar for each atom.  (qBar = complex multi component
      #   vector which starts as the sum of the Ylm (spherical harmonics) for
      #   each bonded atom.)  I know that there is a complex module for perl
      #   but it gave me a lot of trouble and so I am using this lame method.
      foreach $item (1..$numItems1)
      {
         foreach $Ylm_m (1..2*$Ylm_l+1)
         {
            $qBar[$item][$Ylm_m][1] = 0.0;  # Real
            $qBar[$item][$Ylm_m][2] = 0.0;  # Imag
         }
      }

      # Initialize the coefficients for the Ylm.
      &initYlmCoeff;
   }
   elsif ($interactionType == 5) # Min dist of ext items to grid positions.
   {
      # Do not skip any scan points.
      foreach $item1 (1..$numItems1)
         {$skipItem1[$item1] = 0;}

      # Skip items that are not of the requested scan element.
      foreach $item2 (1..$numItems2)
      {
         if ($atomElementName[$ext2CentralItemMap[$item2]] ne $scanElement)
            {$skipItem2[$item2] = 1;}
         else
            {$skipItem2[$item2] = 0;}
      }

      # Do not skip any pairs.
      foreach $item1 (1..$numItems1)
      {
         foreach $item2 (1..$numItems2)
            {$skipPair[$item1][$item2]=0;}
      }

      # Initialize the distance between the mesh point and the extended item.
      foreach $item1 (1..$numItems1)
      {
         foreach $item2 (1..$numItems2)
            {$extScanDist[$item1][$item2] = $bigReal;}
      }
   }
   else
      {die "Unknown interaction type $interactionType.  Aborting.\n";}
}

sub recordInteractionData
{
   # Rename passed parameters for easier reference.
   my $interactionType  = $_[0];
   my $distance         = $_[1];
   my $diff_ref         = $_[2];
   my $item1            = $_[3];
   my $item2            = $_[4];

   # Define local variables
   my $theta;
   my $phi;

   if ($interactionType == 1)  # Minimal distance matrix.
   {
      if ($item1 == $ext2CentralItemMap[$item2])
      {
         # Record the self minimal distance excluding 0 distance.
         if (($distance != 0) && ($distance < $selfMinDist[$item1]))
            {$selfMinDist[$item1] = $distance;}
      }
      else
      {
         # Record the value in the minimal distance matrix.
         if ($distance < $minDist[$item1][$ext2CentralItemMap[$item2]])
            {$minDist[$item1][$ext2CentralItemMap[$item2]] = $distance;}
      }
   }
   elsif ($interactionType == 2) # check for bonding
   {
      # If the distance is <= the sum of the covalent radii of
      #   the two atoms then we say these atoms are bonded.  Note that
      #   the bonding factor was already included when we obtained the
      #   covalent radii.
      if (($distance <= $covalRadii[$item1] +
                        $covalRadii[$ext2CentralItemMap[$item2]]) &&
          ($distance != 0))
      {
         # Increment the number of bonds for the central cell atom.
         $numBonds[$item1]++;

         # Store the periodic item number for this bond.  NOTE:  VERY IMPORTANT
         #   AND VERY TRICKY (and probably very bad).  The value of $item1 is
         #   in reference to the atom number list associated with the central
         #   cell and the value of $item2 is in reference to the atom number
         #   list associated with the periodic system.  Recall that atom 1 in
         #   the central cell list does not correspond with atom 1 in the
         #   periodic list.  Don't get confused!
         $bondedExt[$item1][$numBonds[$item1]] = $item2;

         # Record the distance between these two items.  (Bond Length)
         $bondLengthExt[$item1][$numBonds[$item1]] = $distance;
      }
   }
   elsif ($interactionType == 3) # rpdf
   {
      if (int($distance*100) != 0)
      {
          $rpdf[int($distance*100)] += 1;

#         # Add the distance to the rpdfDistances array. Broadening will occur
#         #   after the array has been filled. Note that the distances are
#         #   are converted to the graph scale (graph point = 100 * distance).
#         push(@rpdfDistances,$distance*100);
      }
   }
   elsif ($interactionType == 4) # Bond orientational order
   {
      # If the distance is <= the sum of the covalent radii of
      #   the two atoms then we say these atoms are bonded.  Note that
      #   the bonding factor was already included when we obtained the
      #   covalent radii.
      if (($distance <= $covalRadii[$item1] +
                        $covalRadii[$ext2CentralItemMap[$item2]]) &&
          ($distance != 0))
      {
         # Increment the number of bonds for the first central cell atom.
         $numBonds[$item1]++;

         # Store the periodic item number for this bond.  NOTE:  VERY IMPORTANT
         #   AND VERY TRICKY (and probably very bad).  The value of $item1 is
         #   in reference to the atom number list associated with the central
         #   cell and the value of $item2 is in reference to the atom number
         #   list associated with the periodic system.  Recall that atom 1 in
         #   the central cell list does not correspond with atom 1 in the
         #   periodic list.  Don't get confused!  This may not be needed for a
         #   bond orientational order calculation.
         $bondedExt[$item1][$numBonds[$item1]] = $item2;

         # Get the angles theta and phi.
         ($theta,$phi) = &sphericalAngles($diff_ref);

         # Accumulate the values of qBar for the atom in the central cell.
         &accumulateQBar($item1,$theta,$phi);
      }
   }
   elsif ($interactionType == 5) # Min dist of ext items to grid positions.
   {
      # Record the distance between the mesh point and the extended item.
      $extScanDist[$item1][$item2] = $distance;
   }
   else
      {die "Unknown interaction type $interactionType.  Aborting.\n";}
}

sub finalizeInteractionData
{
   # Define passed parameters.
   my $interactionType = $_[0];
   my $numItems1       = $_[1];
   my $numItems2       = $_[2];

   if ($interactionType == 1)  # makeinput minimal distance matrix.
   {
      my $item1;
      my $item2;

      # Complete the minimal distance matrix for the self reference.
      foreach my $item (1..$numItems1)
         {$minDist[$item][$item] = 0;}

      # Complete the symmetric portion of the minimal distance matrix.
      foreach $item1 (1..$numItems1-1)
      {
         foreach $item2 ($item1+1..$numItems1)
            {$minDist[$item2][$item1] = $minDist[$item1][$item2];}
      }
   }
   elsif ($interactionType == 2) # check for bonding
      {} # Do nothing
   elsif ($interactionType == 3) # rpdf
   {
      # Create a local point.
      my $point;
      my @rpdfBroadened;

      # Apply Gaussian broadening to the points
      my $rpdfSigma = 0.01;
      &gaussianBroaden($rpdfSigma,1,1000,\@rpdf,\@rpdfBroadened);

      # Divide by r^2 to account for the effect of more atoms being
      #   present for each greater distance.  To divide by r^2 we divide
      #   by index^2 / 100^2 because the index numbers are 100 times greater
      #   than the distances.
      foreach $point (1..1000)
         {$rpdf[$point] = $rpdfBroadened[$point]/($point*$point/100.0/100.0);}
   }
   elsif ($interactionType == 4) # Bond orientational order
   {
      # Divide the qBar for each atom by the number of bonds for that atom.
      &qBarPerBond;

      # Normalize each qBar vector.
      &qBarNormalize;

      # Accumulate the qBar correlation factor.
      &qBarCorrelation;
   }
   elsif ($interactionType == 5) # Min dist of ext items to grid positions.
      {} # Do nothing.
   else
      {die "Unknown interaction type $interactionType.  Aborting.\n";}
}


sub createCoordinationList
{
   # Define local variables.
   my $atom;
   my $bond;
   my $element;
   my $found;
   my @bondedElem;
   my $uniqueBondedElem_ref;
   my $uniqueBondedElemCount_ref;

   # If the length of the numBonds array equals zero, then we have not yet
   #   computed the bonds between atoms and so we should do that.  If however,
   #   the array length is finite, then we have to check to be sure that the
   #   extended calculation is mapped back to the central cell.
   if (scalar(@numBonds) == 0)
      {&createBondingList;}
   elsif (scalar(@bonded) == 0)
      {&mapExtToCentral;}


   foreach $atom (1..$numAtoms)
   {
      # Initialize the coordination for this atom to be empty.
      $coordination[$atom] = "";

      # Compare the position of this atom with the box the user defined (or the
      #   default).  This also considers if the user wanted only atoms inside
      #   the box or only outside.  The return value is 1 if the atom is out
      #   of the desired region and it is 0 if the atom is inside the desired
      #   region.
      if (&itemOutOfBounds($atom,$borderCoords_ref))
         {next;}

      # Construct an array containing the element names of the atoms bonded
      #   to the current $atom.  (Make sure that the array is empty for each
      #   iteration to the next.)
      @bondedElem = ("");
      foreach $bond (1..$numBonds[$atom])
         {$bondedElem[$bond] = $atomElementName[$bonded[$atom][$bond]];}

      # Get the unique elements of the array and a count of how many times each
      #   element appeared in the data.  Note that the flags request
      #   comparison by ASCII characters and for the unique elements to be
      #   sorted by those ASCII characters.
      ($uniqueBondedElem_ref,$uniqueBondedElemCount_ref) =
            &getUniqueArray(\@bondedElem,0,1);

      # Construct the coordination string for this atom.
      foreach $element (1..$#{$uniqueBondedElem_ref})
      {
         $coordination[$atom] = $coordination[$atom] .
                                "$uniqueBondedElemCount_ref->[$element]" .
                                "$uniqueBondedElem_ref->[$element]";
      }
   }
}

# Create a coordination summary for each element in the system.
sub createCoordinationSummary
{
   # Define local variables.
   my $atom;
   my $element;
   my $coordType;
   my @coordinationList;  # Set of coordinations for all atoms of each element.
   my $uniqueCoord_ref;
   my $uniqueCoordCount_ref;

   # Check to be sure that the coordination data has been obtained.
   if (scalar(@coordination) == 0)
      {&createCoordinationList;}

   # Initialize the coordination summary for each element so that the zero
   #   index is empty.
   foreach $element (1..$numElements)
      {push (@{$coordinationList[$element]},"");}

   # Append the coordination of each atom onto the coordination summary for
   #   the associated element. 
   foreach $atom (1..$numAtoms)
   {
      if ($coordination[$atom] eq "")
         {next;}
      push (@{$coordinationList[$atomElementID[$atom]]},$coordination[$atom]);
   }

   # Get the unique coordinations for each element along with the number of
   #   times that each coordination occurs.
   foreach $element (1..$numElements)
   {
      # Initialize the coordination summary for this element to be empty.
      $coordinationSummary[$element] = "";

      # Get an array containing the unique coordinations that occur and another
      #   array with the number of times that each coordination occurs.
      ($uniqueCoord_ref,$uniqueCoordCount_ref) = 
            &getUniqueArray(\@{$coordinationList[$element]},0,1);

      # Construct a string for this element that describes the different
      #   coordinations that this element has.
      foreach $coordType (1..$#{$uniqueCoord_ref})
      {
         $coordinationSummary[$element] = $coordinationSummary[$element] .
               "[$uniqueCoordCount_ref->[$coordType]]" .
               "$uniqueCoord_ref->[$coordType]  ";
      }
   }
}


sub computeBondAnglesExt
{
   # Define local variables.
   my $atom;
   my $bond1;
   my $bond2;
   my $a;
   my $b;
   my $c;
   my $gamma;
   my $cosGamma;
   my $xyzAxis;
   my @abDiff;
   my $sqrDiff;

   foreach $atom (1..$numAtoms)
   {
      # Compare the position of this atom with the box the user defined (or the
      #   default).  This also considers if the user wanted only atoms inside
      #   the box or only outside.  The return value is 1 if the atom is out
      #   of the desired region and it is 0 if the atom is inside the desired
      #   region.
      if (&itemOutOfBounds($atom,$borderCoords_ref))
         {next;}

      # Initialize the count of the number of bond angles for this atom.
      $numBondAngles[$atom] = 0;

      # Use the law of cosines to compute the angle between this vertex atom
      #   and its bonded atom pairs:  c^2 = a^2 + b^2 - 2ab*cos(gamma).
      #   Gamma = the angle we wish to discover, a = distance from vertex atom
      #   to bonded atom #1, b = distance from vertex atom to bonded atom #2
      #   and c = distance between two bonded atoms.

      foreach $bond1 (1..$numBonds[$atom]-1)
      {
         # Short name for the distance between vertex and this bonded atom.
         $a = $bondLengthExt[$atom][$bond1];

         foreach $bond2 ($bond1+1..$numBonds[$atom])
         {
            # Increment the count of the number of bond angles.
            $numBondAngles[$atom]++;

            # Short name for the distance between vertex and this bonded atom.
            $b = $bondLengthExt[$atom][$bond2];

            # Now we must compute the distance between bonded atom #1 and
            #   bonded atom #2.  First we get the difference vector between the
            #   a and b positions, then we find the magnitude of that distance.
            foreach $xyzAxis (1..3)
            {
               $abDiff[$xyzAxis] =
                     $extDirectXYZList[$bondedExt[$atom][$bond1]][$xyzAxis] -
                     $extDirectXYZList[$bondedExt[$atom][$bond2]][$xyzAxis];
            }
            $c = sqrt($abDiff[1]**2 + $abDiff[2]**2 + $abDiff[3]**2);

            # Apply the law of cosines:  c^2 = a^2 + b^2 - 2ab cos(gamma) to
            #   solve for gamma (the angle between two bonds) where c is the
            #   distance between the two bonded atoms, a is the distance from
            #   the central (vertex) atom to one bonded atom and b is the
            #   distance from the vertex to the other.
            $cosGamma = ($a*$a + $b*$b - $c*$c)/(2*$a*$b);
            if (abs($cosGamma + 1.0) < $epsilon)
               {$gamma = $pi;}
            else
               {$gamma = &acos(($a*$a + $b*$b - $c*$c)/(2*$a*$b));}

            # Record for general use in degrees.
            $bondAnglesExt[$atom][$numBondAngles[$atom]] = $gamma*180.0/$pi;
         }
      }
   }
}


sub sphericalAngles
{
   # Define passed parameters.
   my $diff_ref = $_[0];

   # Define local variables.
   my $axis;
   my $theta;
   my $phi;

   # Correct for numerical error.
   foreach $axis (1..3)
   {
      if (abs($diff_ref->[$axis]) < $epsilon)
         {$diff_ref->[$axis] = 0.0;}
   }

   # Compute theta.
   $theta = atan2(sqrt($diff_ref->[1]**2 + $diff_ref->[2]**2),
                  $diff_ref->[3]);

   # Compute Phi.  In the case that phi is undefined we set it to zero.
   if (($diff_ref->[1] == 0.0) && ($diff_ref->[2] == 0.0))
      {$phi = 0.0;}
   else
      {$phi = atan2($diff_ref->[2],$diff_ref->[1]);}

   return ($theta,$phi);
}

# Determine if this item should be excluded because it is not labeled with the
#   requested tag (if such a request was made).
sub itemNotRequested
{
   # Define passed parameters.
   my $item        = $_[0];
   my $tagList     = $_[1];
   my $requestList = $_[2];

   # Define return variable.
   my $skipItem = 1;  # Assume that we will skip this item.

   # Define local variables.
   my $tag;

   # If the request list has negative length then every item is not skipped.
   if ($#$requestList < 0)
      {$skipItem=0;}
   else
   {
      # If the tag of the current item in the tag list is present in the list
      #   of requested tags then we will not skip this item.
      foreach $tag (1..$#$requestList)
      {
         if ($tagList->[$item] eq $requestList->[$tag])
            {$skipItem=0;last;}
      }
   }

   return $skipItem;

#   # In the case of a calculation with a constraint on the types of items to
#   #   consider, we here ignore items other than the requested ones.
#   if ($rpdfSelect == 1)
#   {
#      # Skip to the next atom if the current atom is not one of the selected
#      #   elements.
#      if (($atomElementName[$item] ne $rpdfSelectAtoms[1]) &&
#          ($atomElementName[$item] ne $rpdfSelectAtoms[2]))
#         {$skipItem = 1;}
#   }
#   elsif ($scanElement ne "")
#   {
#      # Skip to the next atom if the current atom is not the element we are
#      #   scanning for.
#      if ($atomElementName[$item] ne $scanElement)
#         {$skipItem = 1;}
#   }
}


# Determine if this atom pair should be excluded because they are not of the
#   correct element pair (if such a request was made).
sub pairNotOfElement
{
   # Define passed variables.
   my $item1    = $_[0];
   my $item2    = $_[1];

   # Define local variables.
   my $skipItem = 0;  # Assume that we will not skip the item2.

   # In the case of an RPDF calculation with a constraint on the types of
   #   elements to consider, we here ignore atom pairs that do not exactly
   #   match the requested pair.
   if ($rpdfSelect == 1)
   {
      # We will skip to the next $j if we do not have a match between
      #   the ($i name and $atom1) + ($j name and $atom2) or
      #   ($i name and $atom2) + ($j name and $atom1).
      if (!((($atomElementName[$item1] eq $rpdfSelectAtoms[1]) &&
             ($atomElementName[$item2] eq $rpdfSelectAtoms[2])) ||
            (($atomElementName[$item1] eq $rpdfSelectAtoms[2]) &&
             ($atomElementName[$item2] eq $rpdfSelectAtoms[1]))))
         {$skipItem = 1;}
   }

   return $skipItem;
}


# Determine if this item should be excluded because it is outside the defined
#   boundaries (if any were defeined).
sub itemOutOfBounds
{
   # Define passed paramters.
   my $item           = $_[0];
   my $itemCoords_ref = $_[1];

   # Define the return variable.
   my $skipItem;

   # Declare local variables.
   my $axis;

   # In the case of a calculation with a constraint on the positions of the
   #   items (e.g. RPDF), we here ignore items other than those either inside
   #   or outside the defined borders as requested.
   if ($borderZone == 1)  # Exclude atoms outside defined borders.
   {
      # Assume that we will not skip this item.
      $skipItem = 0;

      foreach $axis (1..3)
      {
         if (($itemCoords_ref->[$item][$axis] < $borderLow[$axis]-$epsilon) ||
             ($itemCoords_ref->[$item][$axis] > $borderHigh[$axis]+$epsilon))
            {$skipItem = 1;last;}
      }
   }
   elsif ($borderZone == 2) # Exclude atoms inside defined borders.
   {
      # Assume that we will skip this item.
      $skipItem = 1;

      foreach $axis (1..3)
      {
         if (($itemCoords_ref->[$item][$axis] < $borderLow[$axis]-$epsilon) ||
             ($itemCoords_ref->[$item][$axis] > $borderHigh[$axis]+$epsilon))
            {$skipItem = 0;last;}
      }
   }

   return $skipItem;
}

# Given three points, this subroutine will find a normal vector to the plane
#   defined by those points.
sub getPlaneNormal
{
   # Define passed parameters.
   my $point1_ref = $_[0];
   my $point2_ref = $_[1];
   my $point3_ref = $_[2];

   # Define local variables.
   my @normal;
   my @midPoint;
   my @centroid;
   my @vectorDiff1;
   my @vectorDiff2;

   # Determine the center point of the three given points.  This is done by
   #   first finding the midpoint between points 1 and 2.  Then find the point
   #   1/3 of the way between point 3 and the midpoint, closer to the midpoint.
   #   This is the triangle's centroid.
   $midPoint[1] = ($point1_ref->[1] + $point2_ref->[1])/2.0;
   $midPoint[2] = ($point1_ref->[2] + $point2_ref->[2])/2.0;
   $midPoint[3] = ($point1_ref->[3] + $point2_ref->[3])/2.0;

   $centroid[1] = $midPoint[1] + ($point3_ref->[1] - $midPoint[1])/3.0;
   $centroid[2] = $midPoint[2] + ($point3_ref->[2] - $midPoint[2])/3.0;
   $centroid[3] = $midPoint[3] + ($point3_ref->[3] - $midPoint[3])/3.0;

   # Determine the normal to this plane.  This is obtained from the cross
   #   product of two vectors.  The first vector is the difference between
   #   point2 and point1 while the second vector is the difference between
   #   point3 and point1.

   $vectorDiff1[1] = $point2_ref->[1] - $point1_ref->[1];
   $vectorDiff1[2] = $point2_ref->[2] - $point1_ref->[2];
   $vectorDiff1[3] = $point2_ref->[3] - $point1_ref->[3];

   $vectorDiff2[1] = $point3_ref->[1] - $point1_ref->[1];
   $vectorDiff2[2] = $point3_ref->[2] - $point1_ref->[2];
   $vectorDiff2[3] = $point3_ref->[3] - $point1_ref->[3];

   $normal[1] = $vectorDiff1[2]*$vectorDiff2[3]-$vectorDiff1[3]*$vectorDiff2[2];
   $normal[2] =-$vectorDiff1[1]*$vectorDiff2[3]+$vectorDiff1[3]*$vectorDiff2[1];
   $normal[3] = $vectorDiff1[1]*$vectorDiff2[2]-$vectorDiff1[2]*$vectorDiff2[1];

   return (@normal);
}


sub dotProduct
{
   # Define passed parameters.
   my $vector1_ref = $_[0];
   my $vector2_ref = $_[1];

   # Define local variables.
   my $product;

   $product = $vector1_ref->[1]*$vector2_ref->[1] +
              $vector1_ref->[2]*$vector2_ref->[2] +
              $vector1_ref->[3]*$vector2_ref->[3];

   return $product;
}


sub crossProduct
{
   # Define passed parameters.
   my $vector1_ref = $_[0];
   my $vector2_ref = $_[1];

   # Define local variables.
   my @product;

   $product[0] = 0.0;
   $product[1] = $vector1_ref->[2]*$vector2_ref->[3] - 
                 $vector1_ref->[3]*$vector2_ref->[2];
   $product[2] = $vector1_ref->[3]*$vector2_ref->[1] - 
                 $vector1_ref->[1]*$vector2_ref->[3];
   $product[3] = $vector1_ref->[1]*$vector2_ref->[2] - 
                 $vector1_ref->[2]*$vector2_ref->[1];

   return @product;
}


sub crossProductMag
{
   # Define passed parameters.
   my $vector1_ref = $_[0];
   my $vector2_ref = $_[1];

   # Define local variables.
   my @crossProductVector;
   my $productMag;

   @crossProductVector = &crossProduct($vector1_ref,$vector2_ref);

   $productMag = sqrt($crossProductVector[1]**2 +
                      $crossProductVector[2]**2 +
                      $crossProductVector[3]**2);

   return $productMag;
}


sub normalizedCrossProduct
{
   # Define passed parameters.
   my $vector1_ref = $_[0];
   my $vector2_ref = $_[1];

   # Define local variables.
   my $axis;
   my @crossProductVector;
   my $productMag;

   @crossProductVector = &crossProduct($vector1_ref,$vector2_ref);

   $productMag = sqrt($crossProductVector[1]**2 +
                      $crossProductVector[2]**2 +
                      $crossProductVector[3]**2);

   foreach $axis (1..3)
      {$crossProductVector[$axis] /= $productMag;}

   return @crossProductVector;
}

sub getVectorAngle
{
   # Define passed parameters.
   my $vector1_ref = $_[0];
   my $vector2_ref = $_[1];

   # Define local variables.
   my $angle;
   my $mag1;
   my $mag2;

   # Compute the magnitudes of the two vectors.
   $mag1 = sqrt($vector1_ref->[1]*$vector1_ref->[1] +
                $vector1_ref->[2]*$vector1_ref->[2] +
                $vector1_ref->[3]*$vector1_ref->[3]);
   $mag2 = sqrt($vector2_ref->[1]*$vector2_ref->[1] +
                $vector2_ref->[2]*$vector2_ref->[2] +
                $vector2_ref->[3]*$vector2_ref->[3]);

   # Compute the angle between the vectors from:
   #   theta = arccos(v1 . v2) / (|v1| |v2|).

   $angle = &acos(($vector1_ref->[1]*$vector2_ref->[1] +
                   $vector1_ref->[2]*$vector2_ref->[2] +
                   $vector1_ref->[3]*$vector2_ref->[3]) / ($mag1*$mag2));

   # Correct for slight deviations from 90 degrees which may occur because
   #   we do not express some vector elements as infinity. Instead they are
   #   expressed as large numbers.
   if (abs($angle - $pi/2.0) < $epsilon)
      {$angle = $pi/2.0;}

   return $angle;
}

# This subroutine will take any array reference and will create and return a
#   new array that contains only the unique elements of the original array.
#   Uniqueness is determined either according to numerical or character
#   equality.  Also returned is a reference to an array that contains a count
#   of the number of times that each element appeared in the list.
sub getUniqueArray
{
   # Note that we assume the given array to be indexed from 1 to the the number
   #   of elements in the array.  The returned array is also indexed from 1 to
   #   the number of unique elements found.

   # Declare the passed variables.
   my $array_ref = $_[0];
   my $dataForm = $_[1];   # 0 = character; 1 = numerical
   my $doSort = $_[2];     # 0 = no sorting; 1 = sort

   # Declare local variables.
   my @uniqueData;
   my @uniqueCount;
   my $element;
   my $uniqueElement;
   my $foundIndex;
   my $numUnique;
   my @indices;

   # Initialize the number of unique data elements found.
   $numUnique = 0;

   # Consider each element in the given array.  If it is the first element
   #   it is certainly unique so far.  Other elements are added to the list of
   #   unique elements only if they are not found to be in the unique list
   #   already.
   foreach $element (1..$#{$array_ref})
   {
      if ($element == 1)
      {
         $numUnique = 1;
         $uniqueData[1]  = $array_ref->[1];
         $uniqueCount[1] = 1;
      }
      else
      {
         $foundIndex = 0;
         foreach $uniqueElement (1..$numUnique)
         {
            if ($dataForm == 0)
            {
               if ($uniqueData[$uniqueElement] eq $array_ref->[$element])
                  {$foundIndex = $uniqueElement; last;}
            }
            else
            {
               if ($uniqueData[$uniqueElement][1] == $array_ref->[$element])
                  {$foundIndex = $uniqueElement; last;}
            }
         }

         if ($foundIndex == 0)
         {
            $numUnique++;
            $uniqueData[$numUnique]  = $array_ref->[$element];
            $uniqueCount[$numUnique] = 1;
         }
         else
            {$uniqueCount[$foundIndex]++;}
      }
   }

   if ($doSort == 1)
   {
      # Now we will sort the two arrays.  The element names will be sorted in
      #   the requested order and the number of times each element appears will
      #   follow the indices of the names array.  This will be accomplished by
      #   first creating an array of the sorted indices and then assigning
      #   both arrays according to those indices.
      if ($dataForm == 0)
      {
         @indices = sort{$uniqueData[$a] cmp $uniqueData[$b]} 1..$numUnique;
         @uniqueData  = @uniqueData[@indices];
         @uniqueCount = @uniqueCount[@indices];
      }
      else
      {
         @indices = sort{$uniqueData[$a] <=> $uniqueData[$b]} 1..$numUnique;
         @uniqueData  = @uniqueData[@indices];
         @uniqueCount = @uniqueCount[@indices];
      }

      # Make the zero index element empty and have the returned data start from
      #   index one.
      unshift (@uniqueData,"");
      unshift (@uniqueCount,"");
   }

   return (\@uniqueData,\@uniqueCount);
}


sub initYlmCoeff
{
   if ($Ylm_l == 0)
      {$Ylm_Coeff[1] = 0.5 * sqrt(1.0/$pi);}
   elsif ($Ylm_l == 1)
   {
      $Ylm_Coeff[1] = 0.5 * sqrt(3.0/2.0/$pi);
      $Ylm_Coeff[2] = 0.5 * sqrt(3.0/$pi);
      $Ylm_Coeff[3] =-0.5 * sqrt(3.0/2.0/$pi);
   }
   elsif ($Ylm_l == 6)
   {
      # Initialize the values of the Ylm for order 6.
      $Ylm_Coeff[1]  = 1.0/64.0 * sqrt(3003.0/$pi);
      $Ylm_Coeff[2]  = 3.0/32.0 * sqrt(1001.0/$pi);
      $Ylm_Coeff[3]  = 3.0/32.0 * sqrt(91.0/2.0/$pi);
      $Ylm_Coeff[4]  = 1.0/32.0 * sqrt(1365.0/$pi);
      $Ylm_Coeff[5]  = 1.0/64.0 * sqrt(1365.0/$pi);
      $Ylm_Coeff[6]  = 1.0/16.0 * sqrt(273.0/2.0/$pi);
      $Ylm_Coeff[7]  = 1.0/32.0 * sqrt(13.0/$pi);
      $Ylm_Coeff[8]  =-1.0/16.0 * sqrt(273.0/2.0/$pi);
      $Ylm_Coeff[9]  = 1.0/64.0 * sqrt(1365.0/$pi);
      $Ylm_Coeff[10] =-1.0/32.0 * sqrt(1365.0/$pi);
      $Ylm_Coeff[11] = 3.0/32.0 * sqrt(91.0/2.0/$pi);
      $Ylm_Coeff[12] =-3.0/32.0 * sqrt(1001.0/$pi);
      $Ylm_Coeff[13] = 1.0/64.0 * sqrt(3003.0/$pi);
   }
}


sub accumulateQBar
{

   # Define passed parameters.
   my $atom  = $_[0];
   my $theta = $_[1];
   my $phi   = $_[2];

#print STDOUT "atom=$atom theta=$theta   phi=$phi\n";

   if ($Ylm_l == 0)
   {
      $qBar[$atom][1][1] += $Ylm_Coeff[1];
      $qBar[$atom][1][2] += 0.0;
      print STDOUT "$qBar[$atom][1][1]\n";
   }
   elsif ($Ylm_l == 1)
   {
      $qBar[$atom][1][1] += $Ylm_Coeff[1] * cos(-1.0*$phi) * sin($theta);
      $qBar[$atom][2][1] += $Ylm_Coeff[2] * cos($theta);
      $qBar[$atom][3][1] += $Ylm_Coeff[3] * cos(1.0*$phi) * sin($theta);
      $qBar[$atom][1][2] += $Ylm_Coeff[1] * sin(-1.0*$phi) * sin($theta);
      $qBar[$atom][2][2] += 0.0;
      $qBar[$atom][3][2] += $Ylm_Coeff[3] * sin(1.0*$phi) * sin($theta);
#print STDOUT "$qBar[$atom][1][1] $qBar[$atom][2][1] $qBar[$atom][3][1]\n";
#print STDOUT "$qBar[$atom][1][2] $qBar[$atom][2][2] $qBar[$atom][3][2]\n";
   }
   elsif ($Ylm_l == 6)
   {
# REAL
      $qBar[$atom][1][1] += $Ylm_Coeff[1]*cos(-6.0*$phi) *
            sin($theta)**6;
      $qBar[$atom][2][1] += $Ylm_Coeff[2]*cos(-5.0*$phi) *
            sin($theta)**5*cos($theta);
      $qBar[$atom][3][1] += $Ylm_Coeff[3]*cos(-4.0*$phi) *
            sin($theta)**4*(11.0*cos($theta)**2 - 1.0);
      $qBar[$atom][4][1] += $Ylm_Coeff[4]*cos(-3.0*$phi) *
            sin($theta)**3*(11.0*cos($theta)**3 - 3.0*cos($theta));
      $qBar[$atom][5][1] += $Ylm_Coeff[5]*cos(-2.0*$phi) *
            sin($theta)**2*(33.0*cos($theta)**4 - 18.0*cos($theta)**2 + 1.0);
      $qBar[$atom][6][1] += $Ylm_Coeff[6]*cos(-1.0*$phi) * sin($theta) *
            (33.0*cos($theta)**5 - 30.0*cos($theta)**3 + 5.0*cos($theta));
      $qBar[$atom][7][1] += (231.0*cos($theta)**6 - 315.0*cos($theta)**4 +
             105.0*cos($theta)**2 - 5.0) * $Ylm_Coeff[7];
      $qBar[$atom][8][1] += $Ylm_Coeff[8]*cos(1.0*$phi) * sin($theta) *
            (33.0*cos($theta)**5 - 30.0*cos($theta)**3 + 5.0*cos($theta));
      $qBar[$atom][9][1] += $Ylm_Coeff[9]*cos(2.0*$phi) *
            sin($theta)**2*(33.0*cos($theta)**4 - 18.0*cos($theta)**2 + 1.0);
      $qBar[$atom][10][1] += $Ylm_Coeff[10]*cos(3.0*$phi) *
            sin($theta)**3*(11.0*cos($theta)**3 - 3.0*cos($theta));
      $qBar[$atom][11][1] += $Ylm_Coeff[11]*cos(4.0*$phi) *
            sin($theta)**4*(11.0*cos($theta)**2 - 1.0);
      $qBar[$atom][12][1] += $Ylm_Coeff[12]*cos(5.0*$phi) *
            sin($theta)**5*cos($theta);
      $qBar[$atom][13][1] += $Ylm_Coeff[13]*cos(6.0*$phi) *
            sin($theta)**6;

# IMAGINARY
      $qBar[$atom][1][2] += $Ylm_Coeff[1]*sin(-6.0*$phi) *
            sin($theta)**6;
      $qBar[$atom][2][2] += $Ylm_Coeff[2]*sin(-5.0*$phi) *
            sin($theta)**5*cos($theta);
      $qBar[$atom][3][2] += $Ylm_Coeff[3]*sin(-4.0*$phi) *
            sin($theta)**4*(11.0*cos($theta)**2 - 1.0);
      $qBar[$atom][4][2] += $Ylm_Coeff[4]*sin(-3.0*$phi) *
            sin($theta)**3*(11.0*cos($theta)**3 - 3.0*cos($theta));
      $qBar[$atom][5][2] += $Ylm_Coeff[5]*sin(-2.0*$phi) *
            sin($theta)**2*(33.0*cos($theta)**4 - 18.0*cos($theta)**2 + 1.0);
      $qBar[$atom][6][2] += $Ylm_Coeff[6]*sin(-1.0*$phi) * sin($theta) *
            (33.0*cos($theta)**5 - 30.0*cos($theta)**3 + 5.0*cos($theta));
      $qBar[$atom][7][2] += (231.0*cos($theta)**6 - 315.0*cos($theta)**4 +
             105.0*cos($theta)**2 - 5.0) * $Ylm_Coeff[7];
      $qBar[$atom][8][2] += $Ylm_Coeff[8]*sin(1.0*$phi) * sin($theta) *
            (33.0*cos($theta)**5 - 30.0*cos($theta)**3 + 5.0*cos($theta));
      $qBar[$atom][9][2] += $Ylm_Coeff[9]*sin(2.0*$phi) *
            sin($theta)**2*(33.0*cos($theta)**4 - 18.0*cos($theta)**2 + 1.0);
      $qBar[$atom][10][2] += $Ylm_Coeff[10]*sin(3.0*$phi) *
            sin($theta)**3*(11.0*cos($theta)**3 - 3.0*cos($theta));
      $qBar[$atom][11][2] += $Ylm_Coeff[11]*sin(4.0*$phi) *
            sin($theta)**4*(11.0*cos($theta)**2 - 1.0);
      $qBar[$atom][12][2] += $Ylm_Coeff[12]*sin(5.0*$phi) *
            sin($theta)**5*cos($theta);
      $qBar[$atom][13][2] += $Ylm_Coeff[13]*sin(6.0*$phi) *
            sin($theta)**6;
   }
}


sub qBarPerBond
{
   # Define local variables.
   my $atom;
   my $Ylm_m;

   foreach $atom (1..$numAtoms)
   {
#print STDOUT "numBonds[$atom] = $numBonds[$atom]\n";
#print STDOUT "$qBar[$atom][1][1] $qBar[$atom][2][1] $qBar[$atom][3][1]\n";
#print STDOUT "$qBar[$atom][1][2] $qBar[$atom][2][2] $qBar[$atom][3][2]\n";
      foreach $Ylm_m (1..2*$Ylm_l+1)
      {
         $qBar[$atom][$Ylm_m][1] /= $numBonds[$atom];
         $qBar[$atom][$Ylm_m][2] /= $numBonds[$atom];
      }
#print STDOUT "$qBar[$atom][1][1] $qBar[$atom][2][1] $qBar[$atom][3][1]\n";
#print STDOUT "$qBar[$atom][1][2] $qBar[$atom][2][2] $qBar[$atom][3][2]\n";
   }
}


sub qBarNormalize
{
   # Define local variables.
   my $atom;
   my $Ylm_m;
   my $magnitude;

   # Compute the total magnitude of this complex multi component vector.
   foreach $atom (1..$numAtoms)
   {
      # Initialize the magnitude for accumulation.
      $magnitude = 0.0;

      # Accumulate the sum of squares of the vector components.
      foreach $Ylm_m (1..2*$Ylm_l+1)
         {$magnitude+=($qBar[$atom][$Ylm_m][1]**2+$qBar[$atom][$Ylm_m][2]**2);}

      # Obtain the magnitude of the vector in complex space.
      $magnitude = sqrt($magnitude);

#print STDOUT "atom=$atom magnitude=$magnitude\n";
      # Divide each component by the magnitude to make it a normalized vector.
#print STDOUT "$qBar[$atom][1][1] $qBar[$atom][2][1] $qBar[$atom][3][1]\n";
#print STDOUT "$qBar[$atom][1][2] $qBar[$atom][2][2] $qBar[$atom][3][2]\n";
      foreach $Ylm_m (1..2*$Ylm_l+1)
      {
         $qBar[$atom][$Ylm_m][1] /= $magnitude;
         $qBar[$atom][$Ylm_m][2] /= $magnitude;
#         $qBar[$atom][$Ylm_m][1] = abs($qBar[$atom][$Ylm_m][1]) / $magnitude;
#         $qBar[$atom][$Ylm_m][2] = abs($qBar[$atom][$Ylm_m][2]) / $magnitude;
      }
#print STDOUT "$qBar[$atom][1][1] $qBar[$atom][2][1] $qBar[$atom][3][1]\n";
#print STDOUT "$qBar[$atom][1][2] $qBar[$atom][2][2] $qBar[$atom][3][2]\n";
#print STDOUT "$qBar[$atom][1]\n";
   }
}

sub qBarCorrelation
{
   # Define local variables.
   my $atom;
   my $bondedAtom;
   my $Ylm_m;

   &mapExtToCentral;

   foreach $atom (1..$numAtoms)
   {
      # Initialize the bond orientational order parameter for this atom.
      $qOrder[$atom] = 0.0;

      # Accumulate the dot product of $atom's normalized qBar with the
      #   normalized qBar of each bonded atom.
      foreach $bondedAtom (1..$numBonds[$atom])
      {
         foreach $Ylm_m (1..2*$Ylm_l+1)
         {
            $qOrder[$atom] += $qBar[$atom][$Ylm_m][1] *
                              $qBar[$bonded[$atom][$bondedAtom]][$Ylm_m][1] +
                              $qBar[$atom][$Ylm_m][1] *
                             -$qBar[$bonded[$atom][$bondedAtom]][$Ylm_m][2] +
                              $qBar[$atom][$Ylm_m][2] *
                              $qBar[$bonded[$atom][$bondedAtom]][$Ylm_m][1] +
                              $qBar[$atom][$Ylm_m][2] *
                              $qBar[$bonded[$atom][$bondedAtom]][$Ylm_m][2];
#print STDOUT "qOrder[$atom] = $qOrder[$atom]\n";
         }
      }

      # Divide by the total number of bonded atoms.
      $qOrder[$atom] /= $numBonds[$atom];
   }
}


sub stableSort
{
   # Define passed parameters.
   my $value_ref = $_[0];
   my $index_ref = $_[1];
   my $numItems  = $_[2];

   # Define local variables.
   my $i;
   my $j;
   my $tempIndex;
   my $tempValue;

   # Sort by using a simple insertion sort.  Will improve the algorithm if it
   #   becomes a performance issue.  The important thing is that the sort must
   #   be stable, and the indexMap must be created via swapping since it may
   #   not be in sequential order when it arrives in this subroutine.

   for ($i=2;$i<=$numAtoms;++$i)
   {
      $tempIndex = $index_ref->[$i];
      $tempValue = $value_ref->[$i];
      for ($j=$i-1;($j>=1) && ($value_ref->[$j] > $tempValue);--$j)
      {
         $index_ref->[$j+1] = $index_ref->[$j];
         $value_ref->[$j+1] = $value_ref->[$j];
      }
      $index_ref->[$j+1] = $tempIndex;
      $value_ref->[$j+1] = $tempValue;
   }
}


sub prepLine
{
   # Define passed parameters.
   my $fileHandle = $_[0];
   my $line = $_[1];
   my $splitter = $_[2];

   # Declare local variables.
   my @values;

   # Read the line if necessary.
   if ("$fileHandle" ne "")
      {$line = <$fileHandle>;}

   # Chomp, split, and shift it.
   chomp $line;
   @values = split(/$splitter/,$line);
   if ($#values == -1)
      {$values[0] = "";}
   if (($values[0] eq "") && ($#values != 0))
      {shift @values;}

   return @values;
}

sub findNumber
{
   # Define passed parameters.
   my $itemToFind = $_[0];
   my $array = $_[1];
   my $startIndex = $_[2];
   my $endIndex = $_[3];

   # Define local variables.
   my $arrayIndex;
   my $found;

   # Assume that the item will not be found in the array.
   $found = 0;

   # Search the array for the item.
   foreach $arrayIndex ($startIndex..$endIndex)
   {
      if ($array->[$arrayIndex] == $itemToFind)
         {$found = $arrayIndex; last;}
   }

   return $found;
}

sub findString
{
   # Define passed parameters.
   my $itemToFind = $_[0];
   my $array = $_[1];
   my $startIndex = $_[2];
   my $endIndex = $_[3];

   # Define local variables.
   my $arrayIndex;
   my $found;

   # Assume that the item will not be found in the array.
   $found = 0;

   # Search the array for the item.
   foreach $arrayIndex ($startIndex..$endIndex)
   {
      if ($array->[$arrayIndex] eq $itemToFind)
         {$found = $arrayIndex; last;}
   }

   return $found;
}

sub intRemainder
{
   # Define passed parameters.
   my $num1 = $_[0];
   my $num2 = $_[1];

   # Compute and return the remainder.
#   return (($num1/$num2) - int($num1/$num2));
   return ($num2 - int($num2));
}


sub gaussianBroaden
{
   # Define passed parameters.
   my $sigma = $_[0];
   my $start = $_[1]; #Start loop over graph_ref at this value
   my $end = $_[2]; #End loop over graph_ref at this value
   my $points_ref = $_[3]; #Points at which to center broadening Gaussians
                    # Broadening will be applied at each of these points, so
                    # make sure there aren't any that need to be skipped
   my $graph_ref = $_[4]; #Graph scale on which broadening will be accumulated

   # Define local variables
   my $sigmaSqrtPi = $sigma * sqrt($pi);
   my $point;
   my $graphValue;
   my $broadenDiff;
   my $expTerm;
   my @graph;
   my $currPoint;

   foreach $point ($start..$end)
   {

      $currPoint = $points_ref->[$point];

      foreach $graphValue ($start..$end)
      {
         # Determine the difference between the point in question and each
         #   place on the graph scale
         $broadenDiff = ($point - $graphValue)/100.0;

         # Determine the exponential term for the broadening factor. This is
         # the "blah" part of gaussian = normalConst * e^(-blah).
         $expTerm = ($broadenDiff*$broadenDiff) / ($sigma*$sigma);

         # If graphValue is far enough from point that the contribution at
         #   graphValue of a Gaussian centered at point would be negligible,
         #   don't bother completing the broadening.
         # Note that 1/sigmaSqrtPi is the normalization constant.
         if ($expTerm < 20.0)
            {$graph[$graphValue] += $currPoint * exp(-$expTerm) / $sigmaSqrtPi;}
      }
   }

   foreach $graphValue ($start..$end)
      {$graph_ref->[$graphValue] = $graph[$graphValue];}
}


1;
__DATA__
# Below is stub documentation for this module.

=head1 NAME

StructureControl - Perl extension for reading and computing atomic structure
information.

=head1 SYNOPSIS


=head1 DESCRIPTION



=head2 EXPORT

None by default.



=head1 SEE ALSO

No other related documents yet.

There is presently no mailing list for this module.

There is presently no website for this module.

=head1 AUTHOR

Paul Rulis, E<lt>rulisp@umkc.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2007 by Paul Rulis

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.


=cut

__C__

double Cdistance(double limitDistSqrd, double limitDist,
                 double p1x, double p1y, double p1z,
                 double p2x, double p2y, double p2z)
{
   double diff[] = {p2x-p1x, p2y-p1y, p2z-p1z};
   double r;
   int i;
   int tooFar;

   tooFar = 0;
   for (i=0;i<=2;i++)
   {
      if (diff[i] > limitDist)
      {
         tooFar = 1;
         break;
      }
   }

   if (tooFar)
      {return 1000000000.0;}
   else
   {
      r = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
      if (r<limitDistSqrd)
         {return sqrt(r);}
      else
         {return 1000000000.0;}
   }
}

double CdistanceSqrd(double limitDistSqrd, double limitDist,
                     double p1x, double p1y, double p1z,
                     double p2x, double p2y, double p2z)
{
   double diff[] = {p2x-p1x, p2y-p1y, p2z-p1z};
   double r;
   int i;
   int tooFar;

   tooFar = 0;
   for (i=0;i<=2;i++)
   {
      if (diff[i] > limitDist)
      {
         tooFar = 1;
         break;
      }
   }

   if (tooFar)
      {return 1000000000.0;}
   else
   {
      r = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
      if (r<limitDistSqrd)
         {return r;}
      else
         {return 1000000000.0;}
   }
}
