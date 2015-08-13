package ElementData;

use 5.008002;
use strict;
use warnings;
use Env;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use ElementData ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw();

our $VERSION = '0.01';


# Preloaded methods go here.


# Define the arrays that hold elemental data.
my @covalRadii;      # Covalent radii of each element in Angstroms.
my @coreCharge;      # Number of core electrons in each s,p,d,f.
my @valeCharge;      # Number of valence electrons in each s,p,d,f.
my @elementNames;    # Abbreviation of element names from the periodic table.
my @numTermsWF;      # Num gaussians to represent orbital wavefns.
my @minTermWF;       # The minimum exponential gaussian coefficient (wavefn).
my @maxTermWF;       # The maximal exponential gaussian coefficient (wavefn).
my @numTermsPot;     # Num gaussians to represent the potential function.
my @minTermPot;      # The minimum exponential gaussian coefficient (pot).
my @maxTermPot;      # The maximal exponential gaussian coefficient (pot).
my @coreOrbitals;    # Num of s,p,d,f core orbitals for each element.
my @valeOrbitals;    # Num of s,p,d,f valence orbitals [1..3][1..numElements]
                     #   1=MB; 2=FB beyond MB only; 3=EB beyond FB only.
my @orbitalTerms;    # Strings identifying which terms to use for each orbital
                     #   l QN of each atom.
my @colorDX;         # Color in openDX on scale of 1-100.
my @greyDX;          # Greyscale value in openDX on scale of 1-100.

# Define other important variables.
my $maxQN_l;
my $numElements = 103;
my $initData = 0; # Set to 1 when the initialization portion is run.

sub initElementData
{
   # Define local variables.
   my @values;
   my $line;
   my $element;
   my $basis;
   my $QN_l;

   # Check to see if the initialization was run already and skip the rest if so.
   if ($initData == 0)
      {$initData=1;}
   else
      {return;}

   # Open the element data file.
   open (EDATA,"<$OLCAO_DATA/elements.dat") ||
         die "Cannot open elements.dat, aborting.\n";

   # Read the number of elements.
   <EDATA>; # Read past the header.
   @values = &prepLine(\*EDATA,$line,'\s+');
   $numElements = $values[0];

   # Read the max l quantum number of all atoms.
   <EDATA>; # Read past the header.
   @values = &prepLine(\*EDATA,$line,'\s+');
   $maxQN_l = $values[0];

   # Read the element names.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $elementNames[$element] = $values[0];
   }

   # Read the covalent radii.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $covalRadii[$element] = $values[0];
   }

   # Read the designation of core orbitals for each atom.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      foreach $QN_l (0..$maxQN_l)
         {$coreOrbitals[$element][$QN_l] = $values[$QN_l];}
   }

   # Read the designation of minimal basis orbitals beyond the core, followed
   #   by the full beyond the minimal and then the extended beyond the full.
   foreach $basis (1..3)
   {
      <EDATA>; # Read past the header.
      foreach $element (1..$numElements)
      {
         @values = &prepLine(\*EDATA,$line,'\s+');
         foreach $QN_l (0..$maxQN_l)
            {$valeOrbitals[$basis][$element][$QN_l] = $values[$QN_l];}
      }
   }

   # Read the core charge for each spdf orbital of each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      foreach $QN_l (0..$maxQN_l)
         {$coreCharge[$element][$QN_l] = $values[$QN_l];}
   }

   # Read the valence charge for each spdf orbital of each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      foreach $QN_l (0..$maxQN_l)
         {$valeCharge[$element][$QN_l] = $values[$QN_l];}
   }

   # Read the number of radial wave function Gaussians for each spdf orbital
   #   of each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      foreach $QN_l (0..$maxQN_l)
         {$numTermsWF[$element][$QN_l] = $values[$QN_l];}
   }

   # Read the minimum exponential alpha for the radial wave functions for
   #   each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $minTermWF[$element] = $values[0];
   }

   # Read the maximum exponential alpha for the radial wave functions for
   #   each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $maxTermWF[$element] = $values[0];
   }

   # Read the number of potential terms for each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $numTermsPot[$element] = $values[0];
   }

   # Read the minimum exponential alpha for the potential function for
   #   each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $minTermPot[$element] = $values[0];
   }

   # Read the maximum exponential alpha for the potential function for
   #   each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $maxTermPot[$element] = $values[0];
   }

   # Read data strings for each spdf orbital type of each element. The data
   #   strings are boolean identifiers for which Gaussians in the orbital
   #   expansion should be used.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      foreach $QN_l (0..$maxQN_l)
      {
         @values = &prepLine(\*EDATA,$line,'a');# Need to get whole string.
#         if ($values[0] == -1)
#            {$orbitalTerms[$QN_l][$element] = $values[$QN_l];}
          $orbitalTerms[$QN_l][$element] = $values[0];
      }
   }

   # Read the openDX color assignment for each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $colorDX[$element] = $values[0];
   }

   # Read the openDX grey scale assignment for each element.
   <EDATA>; # Read past the header.
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $greyDX[$element] = $values[0];
   }

   close (EDATA);
}

# Define a local version of the prepLine subroutine.  Eventually, this should
#   be redesigned as a module unto itself, independent from StructureControl.pm.
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
   if ($values[0] eq "")
      {shift @values;}

   return @values;
}


#$elementNames[1]   = "h";
#$elementNames[2]   = "he";
#$elementNames[3]   = "li";
#$elementNames[4]   = "be";
#$elementNames[5]   = "b";
#$elementNames[6]   = "c";
#$elementNames[7]   = "n";
#$elementNames[8]   = "o";
#$elementNames[9]   = "f";
#$elementNames[10]  = "ne";
#$elementNames[11]  = "na";
#$elementNames[12]  = "mg";
#$elementNames[13]  = "al";
#$elementNames[14]  = "si";
#$elementNames[15]  = "p";
#$elementNames[16]  = "s";
#$elementNames[17]  = "cl";
#$elementNames[18]  = "ar";
#$elementNames[19]  = "k";
#$elementNames[20]  = "ca";
#$elementNames[21]  = "sc";
#$elementNames[22]  = "ti";
#$elementNames[23]  = "v";
#$elementNames[24]  = "cr";
#$elementNames[25]  = "mn";
#$elementNames[26]  = "fe";
#$elementNames[27]  = "co";
#$elementNames[28]  = "ni";
#$elementNames[29]  = "cu";
#$elementNames[30]  = "zn";
#$elementNames[31]  = "ga";
#$elementNames[32]  = "ge";
#$elementNames[33]  = "as";
#$elementNames[34]  = "se";
#$elementNames[35]  = "br";
#$elementNames[36]  = "kr";
#$elementNames[37]  = "rb";
#$elementNames[38]  = "sr";
#$elementNames[39]  = "y";
#$elementNames[40]  = "zr";
#$elementNames[41]  = "nb";
#$elementNames[42]  = "mo";
#$elementNames[43]  = "tc";
#$elementNames[44]  = "ru";
#$elementNames[45]  = "rh";
#$elementNames[46]  = "pd";
#$elementNames[47]  = "ag";
#$elementNames[48]  = "cd";
#$elementNames[49]  = "in";
#$elementNames[50]  = "sn";
#$elementNames[51]  = "sb";
#$elementNames[52]  = "te";
#$elementNames[53]  = "i";
#$elementNames[54]  = "xe";
#$elementNames[55]  = "cs";
#$elementNames[56]  = "ba";
#$elementNames[57]  = "la";
#$elementNames[58]  = "ce";
#$elementNames[59]  = "pr";
#$elementNames[60]  = "nd";
#$elementNames[61]  = "pm";
#$elementNames[62]  = "sm";
#$elementNames[63]  = "eu";
#$elementNames[64]  = "gd";
#$elementNames[65]  = "tb";
#$elementNames[66]  = "dy";
#$elementNames[67]  = "ho";
#$elementNames[68]  = "er";
#$elementNames[69]  = "tm";
#$elementNames[70]  = "yb";
#$elementNames[71]  = "lu";
#$elementNames[72]  = "hf";
#$elementNames[73]  = "ta";
#$elementNames[74]  = "w";
#$elementNames[75]  = "re";
#$elementNames[76]  = "os";
#$elementNames[77]  = "ir";
#$elementNames[78]  = "pt";
#$elementNames[79]  = "au";
#$elementNames[80]  = "hg";
#$elementNames[81]  = "tl";
#$elementNames[82]  = "pb";
#$elementNames[83]  = "bi";
#$elementNames[84]  = "po";
#$elementNames[85]  = "at";
#$elementNames[86]  = "rn";
#$elementNames[87]  = "fr";
#$elementNames[88]  = "ra";
#$elementNames[89]  = "ac";
#$elementNames[90]  = "th";
#$elementNames[91]  = "pa";
#$elementNames[92]  = "u";
#$elementNames[93]  = "np";
#$elementNames[94]  = "pu";
#$elementNames[95]  = "am";
#$elementNames[96]  = "cm";
#$elementNames[97]  = "bk";
#$elementNames[98]  = "cf";
#$elementNames[99]  = "es";
#$elementNames[100] = "fm";
#$elementNames[101] = "md";
#$elementNames[102] = "no";
#$elementNames[103] = "lr";
#
## Covalent Radii of each element.
#$covalRadii[1]   = 0.32;
#$covalRadii[2]   = 0.93;
#$covalRadii[3]   = 1.23;
#$covalRadii[4]   = 0.90;
#$covalRadii[5]   = 0.82;
#$covalRadii[6]   = 0.77;
#$covalRadii[7]   = 0.75;
#$covalRadii[8]   = 0.73;
#$covalRadii[9]   = 0.72;
#$covalRadii[10]  = 0.71;
#$covalRadii[11]  = 1.54;
#$covalRadii[12]  = 1.36;
#$covalRadii[13]  = 1.18;
#$covalRadii[14]  = 1.11;
#$covalRadii[15]  = 1.06;
#$covalRadii[16]  = 1.02;
#$covalRadii[17]  = 0.99;
#$covalRadii[18]  = 0.98;
#$covalRadii[19]  = 2.03;
#$covalRadii[20]  = 1.74;
#$covalRadii[21]  = 1.44;
#$covalRadii[22]  = 1.32;
#$covalRadii[23]  = 1.22;
#$covalRadii[24]  = 1.18;
#$covalRadii[25]  = 1.17;
#$covalRadii[26]  = 1.17;
#$covalRadii[27]  = 1.16;
#$covalRadii[28]  = 1.15;
#$covalRadii[29]  = 1.17;
#$covalRadii[30]  = 1.25;
#$covalRadii[31]  = 1.26;
#$covalRadii[32]  = 1.22;
#$covalRadii[33]  = 1.20;
#$covalRadii[34]  = 1.16;
#$covalRadii[35]  = 1.14;
#$covalRadii[36]  = 1.12;
#$covalRadii[37]  = 2.16;
#$covalRadii[38]  = 1.91;
#$covalRadii[39]  = 1.62;
#$covalRadii[40]  = 1.45;
#$covalRadii[41]  = 1.34;
#$covalRadii[42]  = 1.30;
#$covalRadii[43]  = 1.27;
#$covalRadii[44]  = 1.25;
#$covalRadii[45]  = 1.25;
#$covalRadii[46]  = 1.28;
#$covalRadii[47]  = 1.34;
#$covalRadii[48]  = 1.48;
#$covalRadii[49]  = 1.44;
#$covalRadii[50]  = 1.41;
#$covalRadii[51]  = 1.40;
#$covalRadii[52]  = 1.36;
#$covalRadii[53]  = 1.33;
#$covalRadii[54]  = 1.31;
#$covalRadii[55]  = 2.35;
#$covalRadii[56]  = 1.98;
#$covalRadii[57]  = 1.69;
#$covalRadii[58]  = 1.65;
#$covalRadii[59]  = 1.65;
#$covalRadii[60]  = 1.64;
#$covalRadii[61]  = 1.63;
#$covalRadii[62]  = 1.62;
#$covalRadii[63]  = 1.85;
#$covalRadii[64]  = 1.61;
#$covalRadii[65]  = 1.59;
#$covalRadii[66]  = 1.59;
#$covalRadii[67]  = 1.58;
#$covalRadii[68]  = 1.57;
#$covalRadii[69]  = 1.56;
#$covalRadii[70]  = 1.74;
#$covalRadii[71]  = 1.56;
#$covalRadii[72]  = 1.44;
#$covalRadii[73]  = 1.34;
#$covalRadii[74]  = 1.30;
#$covalRadii[75]  = 1.28;
#$covalRadii[76]  = 1.26;
#$covalRadii[77]  = 1.27;
#$covalRadii[78]  = 1.30;
#$covalRadii[79]  = 1.34;
#$covalRadii[80]  = 1.49;
#$covalRadii[81]  = 1.48;
#$covalRadii[82]  = 1.47;
#$covalRadii[83]  = 1.46;
#$covalRadii[84]  = 1.46;
#$covalRadii[85]  = 1.45;
#$covalRadii[86]  = 0.00;
#$covalRadii[87]  = 0.00;
#$covalRadii[88]  = 0.00;
#$covalRadii[89]  = 0.00;
#$covalRadii[90]  = 1.65;
#$covalRadii[91]  = 0.00;
#$covalRadii[92]  = 1.42;
#$covalRadii[93]  = 0.00;
#$covalRadii[94]  = 0.00;
#$covalRadii[95]  = 0.00;
#$covalRadii[96]  = 0.00;
#$covalRadii[97]  = 0.00;
#$covalRadii[98]  = 0.00;
#$covalRadii[99]  = 0.00;
#$covalRadii[100] = 0.00;
#$covalRadii[101] = 0.00;
#$covalRadii[102] = 0.00;
#$covalRadii[103] = 0.00;
#
#
#
#@{$coreOrbitals[1]}   = (0,0,0,0);  # No core
#@{$coreOrbitals[2]}   = (0,0,0,0);  # No core
#@{$coreOrbitals[3]}   = (1,0,0,0);  # 1s
#@{$coreOrbitals[4]}   = (1,0,0,0);  # 1s
#@{$coreOrbitals[5]}   = (1,0,0,0);  # 1s
#@{$coreOrbitals[6]}   = (1,0,0,0);  # 1s
#@{$coreOrbitals[7]}   = (1,0,0,0);  # 1s
#@{$coreOrbitals[8]}   = (1,0,0,0);  # 1s
#@{$coreOrbitals[9]}   = (1,0,0,0);  # 1s
#@{$coreOrbitals[10]}  = (1,0,0,0);  # 1s
#@{$coreOrbitals[11]}  = (2,1,0,0);  # 1s,2s 2p
#@{$coreOrbitals[12]}  = (2,1,0,0);  # 1s,2s 2p
#@{$coreOrbitals[13]}  = (2,1,0,0);  # 1s,2s 2p
#@{$coreOrbitals[14]}  = (2,1,0,0);  # 1s,2s 2p
#@{$coreOrbitals[15]}  = (2,1,0,0);  # 1s,2s 2p
#@{$coreOrbitals[16]}  = (2,1,0,0);  # 1s,2s 2p
#@{$coreOrbitals[17]}  = (2,1,0,0);  # 1s,2s 2p
#@{$coreOrbitals[18]}  = (2,1,0,0);  # 1s,2s 2p
#@{$coreOrbitals[19]}  = (3,1,0,0);  # 1s-3s 2p
#@{$coreOrbitals[20]}  = (3,1,0,0);  # 1s-3s 2p
#@{$coreOrbitals[21]}  = (3,1,0,0);  # 1s-3s 2p
#@{$coreOrbitals[22]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[23]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[24]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[25]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[26]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[27]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[28]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[29]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[30]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[31]}  = (3,2,0,0);  # 1s-3s 2p,3p
#@{$coreOrbitals[32]}  = (3,2,1,0);  # 1s-3s 2p,3p 3d
#@{$coreOrbitals[33]}  = (3,2,1,0);  # 1s-3s 2p,3p 3d
#@{$coreOrbitals[34]}  = (3,2,1,0);  # 1s-3s 2p,3p 3d
#@{$coreOrbitals[35]}  = (3,2,1,0);  # 1s-3s 2p,3p 3d
#@{$coreOrbitals[36]}  = (3,2,1,0);  # 1s-3s 2p,3p 3d
#@{$coreOrbitals[37]}  = (4,2,1,0);  # 1s-4s 2p,3p 3d
#@{$coreOrbitals[38]}  = (4,2,1,0);  # 1s-4s 2p,3p 3d
#@{$coreOrbitals[39]}  = (4,2,1,0);  # 1s-4s 2p,3p 3d
#@{$coreOrbitals[40]}  = (4,2,1,0);  # 1s-4s 2p,3p 3d
#@{$coreOrbitals[41]}  = (4,2,1,0);  # 1s-4s 2p,3p 3d
#@{$coreOrbitals[42]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[43]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[44]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[45]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[46]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[47]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[48]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[49]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[50]}  = (4,3,1,0);  # 1s-4s 2p-4p 3d
#@{$coreOrbitals[51]}  = (4,3,2,0);  # 1s-4s 2p-4p 3d,4d
#@{$coreOrbitals[52]}  = (4,3,2,0);  # 1s-4s 2p-4p 3d,4d
#@{$coreOrbitals[53]}  = (4,3,2,0);  # 1s-4s 2p-4p 3d,4d
#@{$coreOrbitals[54]}  = (4,3,2,0);  # 1s-4s 2p-4p 3d,4d
#@{$coreOrbitals[55]}  = (4,3,2,0);  # 1s-4s 2p-4p 3d,4d
#@{$coreOrbitals[56]}  = (4,3,2,0);  # 1s-4s 2p-4p 3d,4d
#@{$coreOrbitals[57]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[58]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[59]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[60]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[61]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[62]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[63]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[64]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[65]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[66]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[67]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[68]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[69]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[70]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[71]}  = (5,3,2,0);  # 1s-5s 2p-4p 3d,4d
#@{$coreOrbitals[72]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[73]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[74]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[75]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[76]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[77]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[78]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[79]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[80]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[81]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[82]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[83]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[84]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[85]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[86]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[87]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[88]}  = (5,4,2,1);  # 1s-5s 2p-5p 3d,4d 4f
#@{$coreOrbitals[89]}  = (6,4,3,1);  # 1s-6s 2p-5p 3d-5d 4f
#@{$coreOrbitals[90]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[91]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[92]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[93]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[94]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[95]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[96]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[97]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[98]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[99]}  = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[100]} = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[101]} = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[102]} = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#@{$coreOrbitals[103]} = (6,5,3,1);  # 1s-6s 2p-6p 3d-5d 4f
#
#
## Minimal Basis beyond the core orbitals.
#@{$valeOrbitals[1][1]}   = (1,0,0,0);  # 1s
#@{$valeOrbitals[1][2]}   = (1,0,0,0);  # 1s
#@{$valeOrbitals[1][3]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[1][4]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[1][5]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[1][6]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[1][7]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[1][8]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[1][9]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[1][10]}  = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[1][11]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[1][12]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[1][13]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[1][14]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[1][15]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[1][16]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[1][17]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[1][18]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[1][19]}  = (1,1,1,0);  # 4s 3p 3d
#@{$valeOrbitals[1][20]}  = (1,1,1,0);  # 4s 3p 3d
#@{$valeOrbitals[1][21]}  = (1,1,1,0);  # 4s 3p 3d
#@{$valeOrbitals[1][22]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][23]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][24]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][25]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][26]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][27]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][28]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][29]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][30]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][31]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[1][32]}  = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[1][33]}  = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[1][34]}  = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[1][35]}  = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[1][36]}  = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[1][37]}  = (1,1,1,0);  # 5s 4p 4d
#@{$valeOrbitals[1][38]}  = (1,1,1,0);  # 5s 4p 4d
#@{$valeOrbitals[1][39]}  = (1,1,1,0);  # 5s 4p 4d
#@{$valeOrbitals[1][40]}  = (1,1,1,0);  # 5s 4p 4d
#@{$valeOrbitals[1][41]}  = (1,1,1,0);  # 5s 4p 4d
#@{$valeOrbitals[1][42]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][43]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][44]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][45]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][46]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][47]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][48]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][49]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][50]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[1][51]}  = (1,1,0,0);  # 5s 5p
#@{$valeOrbitals[1][52]}  = (1,1,0,0);  # 5s 5p
#@{$valeOrbitals[1][53]}  = (1,1,0,0);  # 5s 5p
#@{$valeOrbitals[1][54]}  = (1,1,0,0);  # 5s 5p
#@{$valeOrbitals[1][55]}  = (2,1,1,0);  # 5s,6s 5p 5d
#@{$valeOrbitals[1][56]}  = (2,1,1,0);  # 5s,6s 5p 5d
#@{$valeOrbitals[1][57]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[1][58]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][59]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][60]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][61]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][62]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][63]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][64]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][65]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][66]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][67]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][68]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][69]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][70]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][71]}  = (1,1,1,1);  # 6s 5p 5d 4f
#@{$valeOrbitals[1][72]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][73]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][74]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][75]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][76]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][77]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][78]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][79]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][80]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][81]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][82]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][83]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][84]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][85]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][86]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[1][87]}  = (1,1,1,0);  # 7s 6p 5d
#@{$valeOrbitals[1][88]}  = (1,1,1,0);  # 7s 6p 5d
#@{$valeOrbitals[1][89]}  = (1,1,1,0);  # 7s 6p 5d
#@{$valeOrbitals[1][90]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][91]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][92]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][93]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][94]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][95]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][96]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][97]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][98]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][99]}  = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][100]} = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][101]} = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][102]} = (1,1,1,1);  # 7s 7p 6d 5f
#@{$valeOrbitals[1][103]} = (1,1,1,1);  # 7s 7p 6d 5f
#
#
#
## Full Basis Orbitals (beyond the minimal basis).
#@{$valeOrbitals[2][1]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[2][2]}   = (1,1,0,0);  # 2s 2p
#@{$valeOrbitals[2][3]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[2][4]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[2][5]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[2][6]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[2][7]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[2][8]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[2][9]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[2][10]}  = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[2][11]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[2][12]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[2][13]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[2][14]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[2][15]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[2][16]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[2][17]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[2][18]}  = (1,1,1,0);  # 4s 4p 3d
#@{$valeOrbitals[2][19]}  = (1,1,1,0);  # 5s 4p 4d
#@{$valeOrbitals[2][20]}  = (1,1,1,0);  # 5s 4p 4d
#@{$valeOrbitals[2][21]}  = (1,1,1,0);  # 5s 4p 4d
#@{$valeOrbitals[2][22]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][23]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][24]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][25]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][26]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][27]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][28]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][29]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][30]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][31]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][32]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][33]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][34]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][35]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][36]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[2][37]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[2][38]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[2][39]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[2][40]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[2][41]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[2][42]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][43]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][44]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][45]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][46]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][47]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][48]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][49]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][50]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][51]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][52]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][53]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][54]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[2][55]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][56]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][57]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][58]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][59]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][60]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][61]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][62]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][63]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][64]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][65]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][66]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][67]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][68]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][69]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][70]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][71]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[2][72]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][73]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][74]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][75]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][76]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][77]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][78]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][79]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][80]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][81]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][82]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][83]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][84]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][85]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][86]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[2][87]}  = (1,1,1,0);  # 8s 7p 6d
#@{$valeOrbitals[2][88]}  = (1,1,1,0);  # 8s 7p 6d
#@{$valeOrbitals[2][89]}  = (1,1,1,0);  # 8s 7p 6d
#@{$valeOrbitals[2][90]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][91]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][92]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][93]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][94]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][95]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][96]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][97]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][98]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][99]}  = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][100]} = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][101]} = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][102]} = (1,1,1,1);  # 8s 8p 7d 6f
#@{$valeOrbitals[2][103]} = (1,1,1,1);  # 8s 8p 7d 6f
#
#
#
## Extended Basis Orbitals (beyond the full basis).
#@{$valeOrbitals[3][1]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[3][2]}   = (1,1,0,0);  # 3s 3p
#@{$valeOrbitals[3][3]}   = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[3][4]}   = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[3][5]}   = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[3][6]}   = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[3][7]}   = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[3][8]}   = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[3][9]}   = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[3][10]}  = (1,1,0,0);  # 4s 4p
#@{$valeOrbitals[3][11]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[3][12]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[3][13]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[3][14]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[3][15]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[3][16]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[3][17]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[3][18]}  = (1,1,1,0);  # 5s 5p 4d
#@{$valeOrbitals[3][19]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[3][20]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[3][21]}  = (1,1,1,0);  # 6s 5p 5d
#@{$valeOrbitals[3][22]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][23]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][24]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][25]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][26]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][27]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][28]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][29]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][30]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][31]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][32]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][33]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][34]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][35]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][36]}  = (1,1,1,0);  # 6s 6p 5d
#@{$valeOrbitals[3][37]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[3][38]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[3][39]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[3][40]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[3][41]}  = (1,1,1,0);  # 7s 6p 6d
#@{$valeOrbitals[3][42]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][43]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][44]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][45]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][46]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][47]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][48]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][49]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][50]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][51]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][52]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][53]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][54]}  = (1,1,1,0);  # 7s 7p 6d
#@{$valeOrbitals[3][55]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][56]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][57]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][58]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][59]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][60]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][61]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][62]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][63]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][64]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][65]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][66]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][67]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][68]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][69]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][70]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][71]}  = (1,1,1,0);  # 8s 7p 7d
#@{$valeOrbitals[3][72]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][73]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][74]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][75]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][76]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][77]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][78]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][79]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][80]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][81]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][82]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][83]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][84]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][85]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][86]}  = (1,1,1,0);  # 8s 8p 7d
#@{$valeOrbitals[3][87]}  = (1,1,1,0);  # 9s 8p 7d
#@{$valeOrbitals[3][88]}  = (1,1,1,0);  # 9s 8p 7d
#@{$valeOrbitals[3][89]}  = (1,1,1,0);  # 9s 8p 7d
#@{$valeOrbitals[3][90]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][91]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][92]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][93]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][94]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][95]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][96]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][97]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][98]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][99]}  = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][100]} = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][101]} = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][102]} = (1,1,1,1);  # 9s 9p 8d 7f
#@{$valeOrbitals[3][103]} = (1,1,1,1);  # 9s 9p 8d 7f
#
## Define the number of electrons in the occupied orbitals above the core.
#@{$valeCharge[1]}   = (1,0,0,0);
#@{$valeCharge[2]}   = (2,0,0,0);
#@{$valeCharge[3]}   = (1,0,0,0);
#@{$valeCharge[4]}   = (2,0,0,0);
#@{$valeCharge[5]}   = (2,1,0,0);
#@{$valeCharge[6]}   = (2,2,0,0);
#@{$valeCharge[7]}   = (2,3,0,0);
#@{$valeCharge[8]}   = (2,4,0,0);
#@{$valeCharge[9]}   = (2,5,0,0);
#@{$valeCharge[10]}  = (2,6,0,0);
#@{$valeCharge[11]}  = (1,0,0,0);
#@{$valeCharge[12]}  = (2,0,0,0);
#@{$valeCharge[13]}  = (2,1,0,0);
#@{$valeCharge[14]}  = (2,2,0,0);
#@{$valeCharge[15]}  = (2,3,0,0);
#@{$valeCharge[16]}  = (2,4,0,0);
#@{$valeCharge[17]}  = (2,5,0,0);
#@{$valeCharge[18]}  = (2,6,0,0);
#@{$valeCharge[19]}  = (1,6,0,0);
#@{$valeCharge[20]}  = (2,6,0,0);
#@{$valeCharge[21]}  = (2,6,1,0);
#@{$valeCharge[22]}  = (2,0,2,0);
#@{$valeCharge[23]}  = (2,0,3,0);
#@{$valeCharge[24]}  = (1,0,5,0);
#@{$valeCharge[25]}  = (2,0,5,0);
#@{$valeCharge[26]}  = (2,0,6,0);
#@{$valeCharge[27]}  = (2,0,7,0);
#@{$valeCharge[28]}  = (2,0,8,0);
#@{$valeCharge[29]}  = (1,0,10,0);
#@{$valeCharge[30]}  = (2,0,10,0);
#@{$valeCharge[31]}  = (2,1,10,0);
#@{$valeCharge[32]}  = (2,2,0,0);
#@{$valeCharge[33]}  = (2,3,0,0);
#@{$valeCharge[34]}  = (2,4,0,0);
#@{$valeCharge[35]}  = (2,5,0,0);
#@{$valeCharge[36]}  = (2,6,0,0);
#@{$valeCharge[37]}  = (1,6,0,0);
#@{$valeCharge[38]}  = (2,6,0,0);
#@{$valeCharge[39]}  = (2,6,1,0);
#@{$valeCharge[40]}  = (2,6,2,0);
#@{$valeCharge[41]}  = (1,6,4,0);
#@{$valeCharge[42]}  = (1,0,5,0);
#@{$valeCharge[43]}  = (2,0,5,0);
#@{$valeCharge[44]}  = (1,0,7,0);
#@{$valeCharge[45]}  = (1,0,8,0);
#@{$valeCharge[46]}  = (0,0,10,0);
#@{$valeCharge[47]}  = (1,0,10,0);
#@{$valeCharge[48]}  = (2,0,10,0);
#@{$valeCharge[49]}  = (2,1,10,0);
#@{$valeCharge[50]}  = (2,2,10,0);
#@{$valeCharge[51]}  = (2,3,0,0);
#@{$valeCharge[52]}  = (2,4,0,0);
#@{$valeCharge[53]}  = (2,5,0,0);
#@{$valeCharge[54]}  = (2,6,0,0);
#@{$valeCharge[55]}  = (3,6,0,0); #5s2 6s1
#@{$valeCharge[56]}  = (4,6,0,0); #5s2 6s2
#@{$valeCharge[57]}  = (2,6,1,0);
#@{$valeCharge[58]}  = (2,6,1,1);
#@{$valeCharge[59]}  = (2,6,0,3);
#@{$valeCharge[60]}  = (2,6,0,4);
#@{$valeCharge[61]}  = (2,6,0,5);
#@{$valeCharge[62]}  = (2,6,0,6);
#@{$valeCharge[63]}  = (2,6,0,7);
#@{$valeCharge[64]}  = (2,6,1,7);
#@{$valeCharge[65]}  = (2,6,0,9);
#@{$valeCharge[66]}  = (2,6,0,10);
#@{$valeCharge[67]}  = (2,6,0,11);
#@{$valeCharge[68]}  = (2,6,0,12);
#@{$valeCharge[69]}  = (2,6,0,13);
#@{$valeCharge[70]}  = (2,6,0,14);
#@{$valeCharge[71]}  = (2,6,1,14);
#@{$valeCharge[72]}  = (2,0,2,0);
#@{$valeCharge[73]}  = (2,0,3,0);
#@{$valeCharge[74]}  = (2,0,4,0);
#@{$valeCharge[75]}  = (2,0,5,0);
#@{$valeCharge[76]}  = (2,0,6,0);
#@{$valeCharge[77]}  = (2,0,7,0);
#@{$valeCharge[78]}  = (1,0,9,0);
#@{$valeCharge[79]}  = (1,0,10,0);
#@{$valeCharge[80]}  = (2,0,10,0);
#@{$valeCharge[81]}  = (2,1,10,0);
#@{$valeCharge[82]}  = (2,2,10,0);
#@{$valeCharge[83]}  = (2,3,10,0);
#@{$valeCharge[84]}  = (2,4,10,0);
#@{$valeCharge[85]}  = (2,5,10,0);
#@{$valeCharge[86]}  = (2,6,10,0);
#@{$valeCharge[87]}  = (3,6,10,0); #6s2 7s1
#@{$valeCharge[88]}  = (4,6,10,0); #6s2 7s2
#@{$valeCharge[89]}  = (2,6,1,0);
#@{$valeCharge[90]}  = (2,0,2,0);
#@{$valeCharge[91]}  = (2,0,1,2);
#@{$valeCharge[92]}  = (2,0,1,3);
#@{$valeCharge[93]}  = (2,0,1,4);
#@{$valeCharge[94]}  = (2,0,0,6);
#@{$valeCharge[95]}  = (2,0,0,7);
#@{$valeCharge[96]}  = (2,0,1,7);
#@{$valeCharge[97]}  = (2,0,0,9);
#@{$valeCharge[98]}  = (2,0,0,10);
#@{$valeCharge[99]}  = (2,0,0,11);
#@{$valeCharge[100]} = (2,0,0,12);
#@{$valeCharge[101]} = (2,0,0,13);
#@{$valeCharge[102]} = (2,0,0,14);
#@{$valeCharge[103]} = (2,0,1,14);
#
#@{$numTermsWF[1]}   = (16,16,0,0);
#@{$numTermsWF[2]}   = (16,16,0,0);
#@{$numTermsWF[3]}   = (20,20,0,0);
#@{$numTermsWF[4]}   = (20,20,0,0);
#@{$numTermsWF[5]}   = (20,20,0,0);
#@{$numTermsWF[6]}   = (20,20,0,0);
#@{$numTermsWF[7]}   = (20,20,0,0);
#@{$numTermsWF[8]}   = (20,20,0,0);
#@{$numTermsWF[9]}   = (20,20,0,0);
#@{$numTermsWF[10]}  = (20,20,0,0);
#@{$numTermsWF[11]}  = (21,21,12,0);
#@{$numTermsWF[12]}  = (20,20,12,0);
#@{$numTermsWF[13]}  = (20,20,12,0);
#@{$numTermsWF[14]}  = (21,21,12,0);
#@{$numTermsWF[15]}  = (21,21,12,0);
#@{$numTermsWF[16]}  = (21,21,12,0);
#@{$numTermsWF[17]}  = (21,21,12,0);
#@{$numTermsWF[18]}  = (21,21,12,0);
#@{$numTermsWF[19]}  = (22,22,16,0);
#@{$numTermsWF[20]}  = (22,22,16,0);
#@{$numTermsWF[21]}  = (22,22,16,0);
#@{$numTermsWF[22]}  = (22,22,16,0);
#@{$numTermsWF[23]}  = (22,22,16,0);
#@{$numTermsWF[24]}  = (22,22,16,0);
#@{$numTermsWF[25]}  = (22,22,16,0);
#@{$numTermsWF[26]}  = (22,22,16,0);
#@{$numTermsWF[27]}  = (22,22,16,0);
#@{$numTermsWF[28]}  = (22,22,16,0);
#@{$numTermsWF[29]}  = (22,22,16,0);
#@{$numTermsWF[30]}  = (23,23,17,0);
#@{$numTermsWF[31]}  = (23,23,17,0);
#@{$numTermsWF[32]}  = (23,23,17,0);
#@{$numTermsWF[33]}  = (23,23,17,0);
#@{$numTermsWF[34]}  = (23,23,17,0);
#@{$numTermsWF[35]}  = (23,23,17,0);
#@{$numTermsWF[36]}  = (23,23,17,0);
#@{$numTermsWF[37]}  = (24,24,18,10);
#@{$numTermsWF[38]}  = (24,24,18,10);
#@{$numTermsWF[39]}  = (24,24,18,10);
#@{$numTermsWF[40]}  = (24,24,18,10);
#@{$numTermsWF[41]}  = (24,24,18,10);
#@{$numTermsWF[42]}  = (24,24,18,10);
#@{$numTermsWF[43]}  = (24,24,18,10);
#@{$numTermsWF[44]}  = (24,24,18,10);
#@{$numTermsWF[45]}  = (24,24,18,10);
#@{$numTermsWF[46]}  = (24,24,18,10);
#@{$numTermsWF[47]}  = (24,24,18,10);
#@{$numTermsWF[48]}  = (24,24,18,10);
#@{$numTermsWF[49]}  = (24,24,18,10);
#@{$numTermsWF[50]}  = (24,24,18,10);
#@{$numTermsWF[51]}  = (24,24,18,10);
#@{$numTermsWF[52]}  = (24,24,18,10);
#@{$numTermsWF[53]}  = (24,24,18,10);
#@{$numTermsWF[54]}  = (24,24,18,10);
#@{$numTermsWF[55]}  = (26,26,20,10);
#@{$numTermsWF[56]}  = (26,26,20,10);
#@{$numTermsWF[57]}  = (26,26,20,10);
#@{$numTermsWF[58]}  = (26,26,20,10);
#@{$numTermsWF[59]}  = (26,26,20,10);
#@{$numTermsWF[60]}  = (26,26,20,10);
#@{$numTermsWF[61]}  = (26,26,20,10);
#@{$numTermsWF[62]}  = (26,26,20,10);
#@{$numTermsWF[63]}  = (26,26,20,10);
#@{$numTermsWF[64]}  = (26,26,20,10);
#@{$numTermsWF[65]}  = (26,26,20,10);
#@{$numTermsWF[66]}  = (26,26,20,10);
#@{$numTermsWF[67]}  = (26,26,20,10);
#@{$numTermsWF[68]}  = (26,26,20,10);
#@{$numTermsWF[69]}  = (26,26,20,10);
#@{$numTermsWF[70]}  = (26,26,20,10);
#@{$numTermsWF[71]}  = (26,26,20,10);
#@{$numTermsWF[72]}  = (28,28,22,12);
#@{$numTermsWF[73]}  = (28,28,22,12);
#@{$numTermsWF[74]}  = (28,28,22,12);
#@{$numTermsWF[75]}  = (28,28,22,12);
#@{$numTermsWF[76]}  = (28,28,22,12);
#@{$numTermsWF[77]}  = (28,28,22,12);
#@{$numTermsWF[78]}  = (28,28,22,12);
#@{$numTermsWF[79]}  = (28,28,22,12);
#@{$numTermsWF[80]}  = (28,28,22,12);
#@{$numTermsWF[81]}  = (28,28,22,12);
#@{$numTermsWF[82]}  = (28,28,22,12);
#@{$numTermsWF[83]}  = (28,28,22,12);
#@{$numTermsWF[84]}  = (28,28,22,12);
#@{$numTermsWF[85]}  = (28,28,22,12);
#@{$numTermsWF[86]}  = (28,28,22,12);
#@{$numTermsWF[87]}  = (30,30,24,14);
#@{$numTermsWF[88]}  = (30,30,24,14);
#@{$numTermsWF[89]}  = (30,30,24,14);
#@{$numTermsWF[90]}  = (30,30,24,14);
#@{$numTermsWF[91]}  = (30,30,24,14);
#@{$numTermsWF[92]}  = (30,30,24,14);
#@{$numTermsWF[93]}  = (30,30,24,14);
#@{$numTermsWF[94]}  = (30,30,24,14);
#@{$numTermsWF[95]}  = (30,30,24,14);
#@{$numTermsWF[96]}  = (30,30,24,14);
#@{$numTermsWF[97]}  = (30,30,24,14);
#@{$numTermsWF[98]}  = (30,30,24,14);
#@{$numTermsWF[99]}  = (30,30,24,14);
#@{$numTermsWF[100]} = (30,30,24,14);
#@{$numTermsWF[101]} = (30,30,24,14);
#@{$numTermsWF[102]} = (30,30,24,14);
#@{$numTermsWF[103]} = (30,30,24,14);
#
#
#$minTermWF[1]   = "0.12";
#$minTermWF[2]   = "0.12";
#$minTermWF[3]   = "0.12";
#$minTermWF[4]   = "0.12";
#$minTermWF[5]   = "0.12";
#$minTermWF[6]   = "0.12";
#$minTermWF[7]   = "0.12";
#$minTermWF[8]   = "0.12";
#$minTermWF[9]   = "0.12";
#$minTermWF[10]  = "0.12";
#$minTermWF[11]  = "0.12";
#$minTermWF[12]  = "0.12";
#$minTermWF[13]  = "0.12";
#$minTermWF[14]  = "0.12";
#$minTermWF[15]  = "0.12";
#$minTermWF[16]  = "0.12";
#$minTermWF[17]  = "0.12";
#$minTermWF[18]  = "0.12";
#$minTermWF[19]  = "0.12";
#$minTermWF[20]  = "0.12";
#$minTermWF[21]  = "0.12";
#$minTermWF[22]  = "0.12";
#$minTermWF[23]  = "0.12";
#$minTermWF[24]  = "0.12";
#$minTermWF[25]  = "0.12";
#$minTermWF[26]  = "0.12";
#$minTermWF[27]  = "0.12";
#$minTermWF[28]  = "0.12";
#$minTermWF[29]  = "0.12";
#$minTermWF[30]  = "0.12";
#$minTermWF[31]  = "0.12";
#$minTermWF[32]  = "0.12";
#$minTermWF[33]  = "0.12";
#$minTermWF[34]  = "0.12";
#$minTermWF[35]  = "0.12";
#$minTermWF[36]  = "0.12";
#$minTermWF[37]  = "0.12";
#$minTermWF[38]  = "0.12";
#$minTermWF[39]  = "0.12";
#$minTermWF[40]  = "0.12";
#$minTermWF[41]  = "0.12";
#$minTermWF[42]  = "0.12";
#$minTermWF[43]  = "0.12";
#$minTermWF[44]  = "0.12";
#$minTermWF[45]  = "0.12";
#$minTermWF[46]  = "0.12";
#$minTermWF[47]  = "0.12";
#$minTermWF[48]  = "0.12";
#$minTermWF[49]  = "0.12";
#$minTermWF[50]  = "0.12";
#$minTermWF[51]  = "0.12";
#$minTermWF[52]  = "0.12";
#$minTermWF[53]  = "0.12";
#$minTermWF[54]  = "0.12";
#$minTermWF[55]  = "0.12";
#$minTermWF[56]  = "0.12";
#$minTermWF[57]  = "0.12";
#$minTermWF[58]  = "0.12";
#$minTermWF[59]  = "0.12";
#$minTermWF[60]  = "0.12";
#$minTermWF[61]  = "0.12";
#$minTermWF[62]  = "0.12";
#$minTermWF[63]  = "0.12";
#$minTermWF[64]  = "0.12";
#$minTermWF[65]  = "0.12";
#$minTermWF[66]  = "0.12";
#$minTermWF[67]  = "0.12";
#$minTermWF[68]  = "0.12";
#$minTermWF[69]  = "0.12";
#$minTermWF[70]  = "0.12";
#$minTermWF[71]  = "0.12";
#$minTermWF[72]  = "0.12";
#$minTermWF[73]  = "0.12";
#$minTermWF[74]  = "0.12";
#$minTermWF[75]  = "0.12";
#$minTermWF[76]  = "0.12";
#$minTermWF[77]  = "0.12";
#$minTermWF[78]  = "0.12";
#$minTermWF[79]  = "0.12";
#$minTermWF[80]  = "0.12";
#$minTermWF[81]  = "0.12";
#$minTermWF[82]  = "0.12";
#$minTermWF[83]  = "0.12";
#$minTermWF[84]  = "0.12";
#$minTermWF[85]  = "0.12";
#$minTermWF[86]  = "0.12";
#$minTermWF[87]  = "0.12";
#$minTermWF[88]  = "0.12";
#$minTermWF[89]  = "0.12";
#$minTermWF[90]  = "0.12";
#$minTermWF[91]  = "0.12";
#$minTermWF[92]  = "0.12";
#$minTermWF[93]  = "0.12";
#$minTermWF[94]  = "0.12";
#$minTermWF[95]  = "0.12";
#$minTermWF[96]  = "0.12";
#$minTermWF[97]  = "0.12";
#$minTermWF[98]  = "0.12";
#$minTermWF[99]  = "0.12";
#$minTermWF[100] = "0.12";
#$minTermWF[101] = "0.12";
#$minTermWF[102] = "0.12";
#$minTermWF[103] = "0.12";
#
#
#$maxTermWF[1]   = "10000";
#$maxTermWF[2]   = "10000";
#$maxTermWF[3]   = "50000";
#$maxTermWF[4]   = "50000";
#$maxTermWF[5]   = "50000";
#$maxTermWF[6]   = "50000";
#$maxTermWF[7]   = "50000";
#$maxTermWF[8]   = "50000";
#$maxTermWF[9]   = "50000";
#$maxTermWF[10]  = "50000";
#$maxTermWF[11]  = "100000";
#$maxTermWF[12]  = "100000";
#$maxTermWF[13]  = "100000";
#$maxTermWF[14]  = "200000";
#$maxTermWF[15]  = "100000";
#$maxTermWF[16]  = "100000";
#$maxTermWF[17]  = "100000";
#$maxTermWF[18]  = "100000";
#$maxTermWF[19]  = "500000";
#$maxTermWF[20]  = "500000";
#$maxTermWF[21]  = "500000";
#$maxTermWF[22]  = "500000";
#$maxTermWF[23]  = "500000";
#$maxTermWF[24]  = "500000";
#$maxTermWF[25]  = "500000";
#$maxTermWF[26]  = "500000";
#$maxTermWF[27]  = "500000";
#$maxTermWF[28]  = "500000";
#$maxTermWF[29]  = "500000";
#$maxTermWF[30]  = "1000000";
#$maxTermWF[31]  = "1000000";
#$maxTermWF[32]  = "1000000";
#$maxTermWF[33]  = "1000000";
#$maxTermWF[34]  = "1000000";
#$maxTermWF[35]  = "1000000";
#$maxTermWF[36]  = "1000000";
#$maxTermWF[37]  = "5000000";
#$maxTermWF[38]  = "5000000";
#$maxTermWF[39]  = "5000000";
#$maxTermWF[40]  = "5000000";
#$maxTermWF[41]  = "5000000";
#$maxTermWF[42]  = "5000000";
#$maxTermWF[43]  = "5000000";
#$maxTermWF[44]  = "5000000";
#$maxTermWF[45]  = "5000000";
#$maxTermWF[46]  = "5000000";
#$maxTermWF[47]  = "5000000";
#$maxTermWF[48]  = "5000000";
#$maxTermWF[49]  = "5000000";
#$maxTermWF[50]  = "5000000";
#$maxTermWF[51]  = "5000000";
#$maxTermWF[52]  = "5000000";
#$maxTermWF[53]  = "5000000";
#$maxTermWF[54]  = "5000000";
#$maxTermWF[55]  = "10000000";
#$maxTermWF[56]  = "10000000";
#$maxTermWF[57]  = "10000000";
#$maxTermWF[58]  = "10000000";
#$maxTermWF[59]  = "10000000";
#$maxTermWF[60]  = "10000000";
#$maxTermWF[61]  = "10000000";
#$maxTermWF[62]  = "10000000";
#$maxTermWF[63]  = "10000000";
#$maxTermWF[64]  = "10000000";
#$maxTermWF[65]  = "10000000";
#$maxTermWF[66]  = "10000000";
#$maxTermWF[67]  = "10000000";
#$maxTermWF[68]  = "10000000";
#$maxTermWF[69]  = "10000000";
#$maxTermWF[70]  = "10000000";
#$maxTermWF[71]  = "50000000";
#$maxTermWF[72]  = "50000000";
#$maxTermWF[73]  = "50000000";
#$maxTermWF[74]  = "50000000";
#$maxTermWF[75]  = "50000000";
#$maxTermWF[76]  = "50000000";
#$maxTermWF[77]  = "50000000";
#$maxTermWF[78]  = "50000000";
#$maxTermWF[79]  = "50000000";
#$maxTermWF[80]  = "50000000";
#$maxTermWF[81]  = "50000000";
#$maxTermWF[82]  = "50000000";
#$maxTermWF[83]  = "50000000";
#$maxTermWF[84]  = "50000000";
#$maxTermWF[85]  = "50000000";
#$maxTermWF[86]  = "50000000";
#$maxTermWF[87]  = "50000000";
#$maxTermWF[88]  = "50000000";
#$maxTermWF[89]  = "50000000";
#$maxTermWF[90]  = "50000000";
#$maxTermWF[91]  = "50000000";
#$maxTermWF[92]  = "50000000";
#$maxTermWF[93]  = "50000000";
#$maxTermWF[94]  = "50000000";
#$maxTermWF[95]  = "50000000";
#$maxTermWF[96]  = "50000000";
#$maxTermWF[97]  = "50000000";
#$maxTermWF[98]  = "50000000";
#$maxTermWF[99]  = "50000000";
#$maxTermWF[100] = "50000000";
#$maxTermWF[101] = "50000000";
#$maxTermWF[102] = "50000000";
#$maxTermWF[103] = "50000000";
#
#
## Number of gaussian functions to represent each atomic potential.
#$numTermsPot[1]   = 6;
#$numTermsPot[2]   = 10;
#$numTermsPot[3]   = 16;
#$numTermsPot[4]   = 16;
#$numTermsPot[5]   = 16;
#$numTermsPot[6]   = 16;
#$numTermsPot[7]   = 16;
#$numTermsPot[8]   = 16;
#$numTermsPot[9]   = 16;
#$numTermsPot[10]  = 16;
#$numTermsPot[11]  = 16;
#$numTermsPot[12]  = 16;
#$numTermsPot[13]  = 16;
#$numTermsPot[14]  = 16;
#$numTermsPot[15]  = 20;
#$numTermsPot[16]  = 20;
#$numTermsPot[17]  = 20;
#$numTermsPot[18]  = 20;
#$numTermsPot[19]  = 20;
#$numTermsPot[20]  = 22;
#$numTermsPot[21]  = 22;
#$numTermsPot[22]  = 22;
#$numTermsPot[23]  = 22;
#$numTermsPot[24]  = 22;
#$numTermsPot[25]  = 23;
#$numTermsPot[26]  = 23;
#$numTermsPot[27]  = 22;
#$numTermsPot[28]  = 27;
#$numTermsPot[29]  = 27;
#$numTermsPot[30]  = 26;
#$numTermsPot[31]  = 25;
#$numTermsPot[32]  = 24;
#$numTermsPot[33]  = 21;
#$numTermsPot[34]  = 23;
#$numTermsPot[35]  = 24;
#$numTermsPot[36]  = 26;
#$numTermsPot[37]  = 26;
#$numTermsPot[38]  = 24;
#$numTermsPot[39]  = 24;
#$numTermsPot[40]  = 26;
#$numTermsPot[41]  = 26;
#$numTermsPot[42]  = 26;
#$numTermsPot[43]  = 26;
#$numTermsPot[44]  = 26;
#$numTermsPot[45]  = 26;
#$numTermsPot[46]  = 26;
#$numTermsPot[47]  = 26;
#$numTermsPot[48]  = 28;
#$numTermsPot[49]  = 32;
#$numTermsPot[50]  = 34;
#$numTermsPot[51]  = 34;
#$numTermsPot[52]  = 34;
#$numTermsPot[53]  = 34;
#$numTermsPot[54]  = 34;
#$numTermsPot[55]  = 32;
#$numTermsPot[56]  = 31;
#$numTermsPot[57]  = 32;
#$numTermsPot[58]  = 32;
#$numTermsPot[59]  = 32;
#$numTermsPot[60]  = 32;
#$numTermsPot[61]  = 32;
#$numTermsPot[62]  = 32;
#$numTermsPot[63]  = 32;
#$numTermsPot[64]  = 32;
#$numTermsPot[65]  = 32;
#$numTermsPot[66]  = 32;
#$numTermsPot[67]  = 32;
#$numTermsPot[68]  = 32;
#$numTermsPot[69]  = 32;
#$numTermsPot[70]  = 31;
#$numTermsPot[71]  = 35;
#$numTermsPot[72]  = 38;
#$numTermsPot[73]  = 42;
#$numTermsPot[74]  = 50;
#$numTermsPot[75]  = 42;
#$numTermsPot[76]  = 42;
#$numTermsPot[77]  = 42;
#$numTermsPot[78]  = 42;
#$numTermsPot[79]  = 40;
#$numTermsPot[80]  = 42;
#$numTermsPot[81]  = 42;
#$numTermsPot[82]  = 42;
#$numTermsPot[83]  = 40;
#$numTermsPot[84]  = 42;
#$numTermsPot[85]  = 42;
#$numTermsPot[86]  = 42;
#$numTermsPot[87]  = 42;
#$numTermsPot[88]  = 42;
#$numTermsPot[89]  = 42;
#$numTermsPot[90]  = 42;
#$numTermsPot[91]  = 42;
#$numTermsPot[92]  = 42;
#$numTermsPot[93]  = 42;
#$numTermsPot[94]  = 42;
#$numTermsPot[95]  = 42;
#$numTermsPot[96]  = 42;
#$numTermsPot[97]  = 42;
#$numTermsPot[98]  = 42;
#$numTermsPot[99]  = 42;
#$numTermsPot[100] = 42;
#$numTermsPot[101] = 42;
#$numTermsPot[102] = 42;
#$numTermsPot[103] = 42;
#
#
## Define the minimum Gaussian alpha for the potential fns of each element.
#$minTermPot[1]   = 0.2;
#$minTermPot[2]   = 0.3;
#$minTermPot[3]   = 0.1;
#$minTermPot[4]   = 0.2;
#$minTermPot[5]   = 0.2;
#$minTermPot[6]   = 0.2;
#$minTermPot[7]   = 0.25;
#$minTermPot[8]   = 0.3;
#$minTermPot[9]   = 0.3;
#$minTermPot[10]  = 0.3;
#$minTermPot[11]  = 0.1;
#$minTermPot[12]  = 0.15;
#$minTermPot[13]  = 0.15;
#$minTermPot[14]  = 0.15;
#$minTermPot[15]  = 0.15;
#$minTermPot[16]  = 0.2;
#$minTermPot[17]  = 0.25;
#$minTermPot[18]  = 0.3;
#$minTermPot[19]  = 0.08;
#$minTermPot[20]  = 0.08;
#$minTermPot[21]  = 0.135;
#$minTermPot[22]  = 0.15;
#$minTermPot[23]  = 0.15;
#$minTermPot[24]  = 0.15;
#$minTermPot[25]  = 0.15;
#$minTermPot[26]  = 0.15;
#$minTermPot[27]  = 0.15;
#$minTermPot[28]  = 0.15;
#$minTermPot[29]  = 0.15;
#$minTermPot[30]  = 0.12;
#$minTermPot[31]  = 0.15;
#$minTermPot[32]  = 0.10;
#$minTermPot[33]  = 0.20;
#$minTermPot[34]  = 0.22;
#$minTermPot[35]  = 0.25;
#$minTermPot[36]  = 0.3;
#$minTermPot[37]  = 0.1;
#$minTermPot[38]  = 0.1;
#$minTermPot[39]  = 0.15;
#$minTermPot[40]  = 0.1;
#$minTermPot[41]  = 0.15;
#$minTermPot[42]  = 0.15;
#$minTermPot[43]  = 0.15;
#$minTermPot[44]  = 0.15;
#$minTermPot[45]  = 0.15;
#$minTermPot[46]  = 0.15;
#$minTermPot[47]  = 0.15;
#$minTermPot[48]  = 0.15;
#$minTermPot[49]  = 0.12;
#$minTermPot[50]  = 0.15;
#$minTermPot[51]  = 0.1;
#$minTermPot[52]  = 0.1;
#$minTermPot[53]  = 0.25;
#$minTermPot[54]  = 0.3;
#$minTermPot[55]  = 0.1;
#$minTermPot[56]  = 0.1;
#$minTermPot[57]  = 0.15;
#$minTermPot[58]  = 0.15;
#$minTermPot[59]  = 0.15;
#$minTermPot[60]  = 0.15;
#$minTermPot[61]  = 0.15;
#$minTermPot[62]  = 0.15;
#$minTermPot[63]  = 0.15;
#$minTermPot[64]  = 0.15;
#$minTermPot[65]  = 0.15;
#$minTermPot[66]  = 0.15;
#$minTermPot[67]  = 0.15;
#$minTermPot[68]  = 0.15;
#$minTermPot[69]  = 0.15;
#$minTermPot[70]  = 0.15;
#$minTermPot[71]  = 0.15;
#$minTermPot[72]  = 0.15;
#$minTermPot[73]  = 0.15;
#$minTermPot[74]  = 0.15;
#$minTermPot[75]  = 0.15;
#$minTermPot[76]  = 0.15;
#$minTermPot[77]  = 0.15;
#$minTermPot[78]  = 0.15;
#$minTermPot[79]  = 0.15;
#$minTermPot[80]  = 0.15;
#$minTermPot[81]  = 0.15;
#$minTermPot[82]  = 0.15;
#$minTermPot[83]  = 0.15;
#$minTermPot[84]  = 0.15;
#$minTermPot[85]  = 0.15;
#$minTermPot[86]  = 0.2;
#$minTermPot[87]  = 0.1;
#$minTermPot[88]  = 0.1;
#$minTermPot[89]  = 0.1;
#$minTermPot[90]  = 0.15;
#$minTermPot[91]  = 0.15;
#$minTermPot[92]  = 0.15;
#$minTermPot[93]  = 0.15;
#$minTermPot[94]  = 0.15;
#$minTermPot[95]  = 0.15;
#$minTermPot[96]  = 0.15;
#$minTermPot[97]  = 0.15;
#$minTermPot[98]  = 0.15;
#$minTermPot[99]  = 0.15;
#$minTermPot[100] = 0.15;
#$minTermPot[101] = 0.15;
#$minTermPot[102] = 0.15;
#$minTermPot[103] = 0.15;
#
#
## Define the minimum Gaussian alpha for the potential fns of each element.
#$maxTermPot[1]   = 1000;
#$maxTermPot[2]   = 1000;
#$maxTermPot[3]   = 1000000;
#$maxTermPot[4]   = 1000000;
#$maxTermPot[5]   = 1000000;
#$maxTermPot[6]   = 1000000;
#$maxTermPot[7]   = 1000000;
#$maxTermPot[8]   = 1000000;
#$maxTermPot[9]   = 1000000;
#$maxTermPot[10]  = 1000000;
#$maxTermPot[11]  = 1000000;
#$maxTermPot[12]  = 1000000;
#$maxTermPot[13]  = 1000000;
#$maxTermPot[14]  = 1000000;
#$maxTermPot[15]  = 1000000;
#$maxTermPot[16]  = 1000000;
#$maxTermPot[17]  = 1000000;
#$maxTermPot[18]  = 1000000;
#$maxTermPot[19]  = 1000000;
#$maxTermPot[20]  = 1000000;
#$maxTermPot[21]  = 1000000;
#$maxTermPot[22]  = 10000000;
#$maxTermPot[23]  = 10000000;
#$maxTermPot[24]  = 10000000;
#$maxTermPot[25]  = 10000000;
#$maxTermPot[26]  = 10000000;
#$maxTermPot[27]  = 10000000;
#$maxTermPot[28]  = 10000000;
#$maxTermPot[29]  = 10000000;
#$maxTermPot[30]  = 10000000;
#$maxTermPot[31]  = 10000000;
#$maxTermPot[32]  = 10000000;
#$maxTermPot[33]  = 10000000;
#$maxTermPot[34]  = 10000000;
#$maxTermPot[35]  = 10000000;
#$maxTermPot[36]  = 10000000;
#$maxTermPot[37]  = 10000000;
#$maxTermPot[38]  = 10000000;
#$maxTermPot[39]  = 10000000;
#$maxTermPot[40]  = 10000000;
#$maxTermPot[41]  = 10000000;
#$maxTermPot[42]  = 10000000;
#$maxTermPot[43]  = 10000000;
#$maxTermPot[44]  = 10000000;
#$maxTermPot[45]  = 10000000;
#$maxTermPot[46]  = 10000000;
#$maxTermPot[47]  = 10000000;
#$maxTermPot[48]  = 10000000;
#$maxTermPot[49]  = 100000000;
#$maxTermPot[50]  = 100000000;
#$maxTermPot[51]  = 100000000;
#$maxTermPot[52]  = 100000000;
#$maxTermPot[53]  = 100000000;
#$maxTermPot[54]  = 100000000;
#$maxTermPot[55]  = 100000000;
#$maxTermPot[56]  = 100000000;
#$maxTermPot[57]  = 100000000;
#$maxTermPot[58]  = 100000000;
#$maxTermPot[59]  = 100000000;
#$maxTermPot[60]  = 100000000;
#$maxTermPot[61]  = 100000000;
#$maxTermPot[62]  = 100000000;
#$maxTermPot[63]  = 100000000;
#$maxTermPot[64]  = 100000000;
#$maxTermPot[65]  = 100000000;
#$maxTermPot[66]  = 100000000;
#$maxTermPot[67]  = 100000000;
#$maxTermPot[68]  = 100000000;
#$maxTermPot[69]  = 100000000;
#$maxTermPot[70]  = 100000000;
#$maxTermPot[71]  = 100000000;
#$maxTermPot[72]  = 100000000;
#$maxTermPot[73]  = 100000000;
#$maxTermPot[74]  = 100000000;
#$maxTermPot[75]  = 100000000;
#$maxTermPot[76]  = 100000000;
#$maxTermPot[77]  = 100000000;
#$maxTermPot[78]  = 100000000;
#$maxTermPot[79]  = 100000000;
#$maxTermPot[80]  = 100000000;
#$maxTermPot[81]  = 100000000;
#$maxTermPot[82]  = 100000000;
#$maxTermPot[83]  = 1000000000;
#$maxTermPot[84]  = 1000000000;
#$maxTermPot[85]  = 1000000000;
#$maxTermPot[86]  = 1000000000;
#$maxTermPot[87]  = 1000000000;
#$maxTermPot[88]  = 1000000000;
#$maxTermPot[89]  = 1000000000;
#$maxTermPot[90]  = 1000000000;
#$maxTermPot[91]  = 1000000000;
#$maxTermPot[92]  = 1000000000;
#$maxTermPot[93]  = 1000000000;
#$maxTermPot[94]  = 1000000000;
#$maxTermPot[95]  = 1000000000;
#$maxTermPot[96]  = 1000000000;
#$maxTermPot[97]  = 1000000000;
#$maxTermPot[98]  = 1000000000;
#$maxTermPot[99]  = 1000000000;
#$maxTermPot[100] = 1000000000;
#$maxTermPot[101] = 1000000000;
#$maxTermPot[102] = 1000000000;
#$maxTermPot[103] = 1000000000;
#
#
## Which gaussian terms to use for each QN_l orbital type.  Index 1=orbital
##   type (0=s, 1=p, 2=d, 3=f), index 2=element ID number.
#$orbitalTerms[0][1]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][1]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][2]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][2]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][3]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][3]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][4]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][4]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][5]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][5]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][6]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][6]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][7]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][7]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][8]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][8]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][9]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][9]  = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][10] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][10] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[0][11] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][11] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][11] = "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][12] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][12] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][12] = "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][13] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][13] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][13] = "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][14] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][14] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][14] = "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][15] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][15] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][15] = "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][16] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][16] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][16] = "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][17] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][17] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][17] = "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][18] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][18] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][18] = "1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][19] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][19] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][19] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][20] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][20] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][20] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][21] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][21] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][21] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][22] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][22] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][22] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][23] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][23] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][23] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][24] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][24] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][24] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][25] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][25] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][25] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][26] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][26] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][26] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][27] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][27] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][27] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][28] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][28] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][28] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][29] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][29] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][29] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][30] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][30] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][30] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][31] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][31] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][31] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][32] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][32] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][32] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][33] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][33] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][33] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][34] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][34] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][34] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][35] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][35] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][35] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][36] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][36] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][36] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[0][37] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][37] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][37] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][37] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][38] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][38] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][38] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][38] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][39] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][39] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][39] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][39] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][40] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][40] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][40] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][40] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][41] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][41] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][41] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][41] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][42] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][42] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][42] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][42] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][43] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][43] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][43] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][43] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][44] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][44] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][44] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][44] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][45] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][45] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][45] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][45] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][46] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][46] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][46] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][46] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][47] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][47] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][47] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][47] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][48] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][48] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][48] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][48] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][49] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][49] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][49] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][49] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][50] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][50] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][50] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][50] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][51] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][51] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][51] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][51] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][52] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][52] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][52] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][52] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][53] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][53] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][53] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][53] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][54] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][54] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][54] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][54] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][55] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][55] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][55] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][55] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][56] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][56] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][56] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][56] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][57] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][57] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][57] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][57] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][58] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][58] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][58] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][58] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][59] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][59] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][59] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][59] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][60] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][60] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][60] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][60] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][61] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][61] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][61] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][61] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][62] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][62] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][62] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][62] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][63] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][63] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][63] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][63] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][64] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][64] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][64] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][64] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][65] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][65] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][65] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][65] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][66] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][66] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][66] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][66] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][67] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][67] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][67] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][67] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][68] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][68] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][68] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][68] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][69] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][69] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][69] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][69] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][70] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][70] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][70] = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][70] = "1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][71] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][71] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][71] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][71] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][72] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][72] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][72] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][72] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][73] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][73] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][73] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][73] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][74] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][74] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][74] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][74] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][75] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][75] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][75] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][75] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][76] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][76] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][76] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][76] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][77] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][77] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][77] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][77] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][78] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][78] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][78] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][78] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][79] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][79] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][79] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][79] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][80] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][80] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][80] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][80] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][81] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][81] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][81] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][81] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][82] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][82] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][82] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][82] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][83] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][83] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][83] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][83] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][84] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][84] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][84] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][84] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][85] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][85] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][85] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][85] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][86] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][86] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][86] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0";
#$orbitalTerms[3][86] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][87] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][87] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][87] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][87] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][88] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][88] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][88] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][88] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][89] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][89] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][89] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][89] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][90] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][90] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][90] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][90] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][91] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][91] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][91] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][91] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][92] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][92] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][92] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][92] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][93] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][93] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][93] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][93] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][94] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][94] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][94] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][94] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][95] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][95] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][95] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][95] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][96] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][96] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][96] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][96] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][97] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][97] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][97] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][97] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][98] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][98] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][98] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][98] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][99] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][99] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][99] ="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][99] ="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][100]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][100]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][100]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][100]="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][101]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][101]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][101]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][101]="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][102]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][102]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][102]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][102]="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#$orbitalTerms[0][103]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[1][103]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1";
#$orbitalTerms[2][103]="1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0";
#$orbitalTerms[3][103]="1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
#
#
## Assign color values to each atom.
#$colorDX[1]   = 10;
#$colorDX[2]   = 10;
#$colorDX[3]   = 75;
#$colorDX[4]   = 10;
#$colorDX[5]   = 10;
#$colorDX[6]   = 66;
#$colorDX[7]   = 1 ;
#$colorDX[8]   = 100;
#$colorDX[9]   = 10;
#$colorDX[10]  = 10;
#$colorDX[11]  = 10;
#$colorDX[12]  = 10;
#$colorDX[13]  = 10;
#$colorDX[14]  = 10;
#$colorDX[15]  = 10;
#$colorDX[16]  = 10;
#$colorDX[17]  = 10;
#$colorDX[18]  = 10;
#$colorDX[19]  = 10;
#$colorDX[20]  = 10;
#$colorDX[21]  = 10;
#$colorDX[22]  = 10;
#$colorDX[23]  = 10;
#$colorDX[24]  = 10;
#$colorDX[25]  = 10;
#$colorDX[26]  = 10;
#$colorDX[27]  = 10;
#$colorDX[28]  = 10;
#$colorDX[29]  = 10;
#$colorDX[30]  = 10;
#$colorDX[31]  = 10;
#$colorDX[32]  = 10;
#$colorDX[33]  = 10;
#$colorDX[34]  = 10;
#$colorDX[35]  = 10;
#$colorDX[36]  = 10;
#$colorDX[37]  = 10;
#$colorDX[38]  = 10;
#$colorDX[39]  = 10;
#$colorDX[40]  = 10;
#$colorDX[41]  = 10;
#$colorDX[42]  = 10;
#$colorDX[43]  = 10;
#$colorDX[44]  = 10;
#$colorDX[45]  = 10;
#$colorDX[46]  = 10;
#$colorDX[47]  = 10;
#$colorDX[48]  = 10;
#$colorDX[49]  = 10;
#$colorDX[50]  = 10;
#$colorDX[51]  = 10;
#$colorDX[52]  = 10;
#$colorDX[53]  = 10;
#$colorDX[54]  = 10;
#$colorDX[55]  = 10;
#$colorDX[56]  = 10;
#$colorDX[57]  = 10;
#$colorDX[58]  = 10;
#$colorDX[59]  = 10;
#$colorDX[60]  = 10;
#$colorDX[61]  = 10;
#$colorDX[62]  = 10;
#$colorDX[63]  = 10;
#$colorDX[64]  = 10;
#$colorDX[65]  = 10;
#$colorDX[66]  = 10;
#$colorDX[67]  = 10;
#$colorDX[68]  = 10;
#$colorDX[69]  = 10;
#$colorDX[70]  = 10;
#$colorDX[71]  = 10;
#$colorDX[72]  = 10;
#$colorDX[73]  = 10;
#$colorDX[74]  = 10;
#$colorDX[75]  = 10;
#$colorDX[76]  = 10;
#$colorDX[77]  = 10;
#$colorDX[78]  = 10;
#$colorDX[79]  = 10;
#$colorDX[80]  = 10;
#$colorDX[81]  = 10;
#$colorDX[82]  = 10;
#$colorDX[83]  = 10;
#$colorDX[84]  = 10;
#$colorDX[85]  = 10;
#$colorDX[86]  = 10;
#$colorDX[87]  = 10;
#$colorDX[88]  = 10;
#$colorDX[89]  = 10;
#$colorDX[90]  = 10;
#$colorDX[91]  = 10;
#$colorDX[92]  = 10;
#$colorDX[93]  = 10;
#$colorDX[94]  = 10;
#$colorDX[95]  = 10;
#$colorDX[96]  = 10;
#$colorDX[97]  = 10;
#$colorDX[98]  = 10;
#$colorDX[99]  = 10;
#$colorDX[100] = 10;
#$colorDX[101] = 10;
#$colorDX[102] = 10;
#$colorDX[103] = 10;
#;
#;
#;
## Assign grey scale values to each atom.;
#$greyDX[1]   = 10;
#$greyDX[2]   = 10;
#$greyDX[3]   = 75;
#$greyDX[4]   = 10;
#$greyDX[5]   = 10;
#$greyDX[6]   = 66;
#$greyDX[7]   = 1 ;
#$greyDX[8]   = 100;
#$greyDX[9]   = 10;
#$greyDX[10]  = 10;
#$greyDX[11]  = 10;
#$greyDX[12]  = 10;
#$greyDX[13]  = 10;
#$greyDX[14]  = 10;
#$greyDX[15]  = 10;
#$greyDX[16]  = 10;
#$greyDX[17]  = 10;
#$greyDX[18]  = 10;
#$greyDX[19]  = 10;
#$greyDX[20]  = 10;
#$greyDX[21]  = 10;
#$greyDX[22]  = 10;
#$greyDX[23]  = 10;
#$greyDX[24]  = 10;
#$greyDX[25]  = 10;
#$greyDX[26]  = 10;
#$greyDX[27]  = 10;
#$greyDX[28]  = 10;
#$greyDX[29]  = 10;
#$greyDX[30]  = 10;
#$greyDX[31]  = 10;
#$greyDX[32]  = 10;
#$greyDX[33]  = 10;
#$greyDX[34]  = 10;
#$greyDX[35]  = 10;
#$greyDX[36]  = 10;
#$greyDX[37]  = 10;
#$greyDX[38]  = 10;
#$greyDX[39]  = 10;
#$greyDX[40]  = 10;
#$greyDX[41]  = 10;
#$greyDX[42]  = 10;
#$greyDX[43]  = 10;
#$greyDX[44]  = 10;
#$greyDX[45]  = 10;
#$greyDX[46]  = 10;
#$greyDX[47]  = 10;
#$greyDX[48]  = 10;
#$greyDX[49]  = 10;
#$greyDX[50]  = 10;
#$greyDX[51]  = 10;
#$greyDX[52]  = 10;
#$greyDX[53]  = 10;
#$greyDX[54]  = 10;
#$greyDX[55]  = 10;
#$greyDX[56]  = 10;
#$greyDX[57]  = 10;
#$greyDX[58]  = 10;
#$greyDX[59]  = 10;
#$greyDX[60]  = 10;
#$greyDX[61]  = 10;
#$greyDX[62]  = 10;
#$greyDX[63]  = 10;
#$greyDX[64]  = 10;
#$greyDX[65]  = 10;
#$greyDX[66]  = 10;
#$greyDX[67]  = 10;
#$greyDX[68]  = 10;
#$greyDX[69]  = 10;
#$greyDX[70]  = 10;
#$greyDX[71]  = 10;
#$greyDX[72]  = 10;
#$greyDX[73]  = 10;
#$greyDX[74]  = 10;
#$greyDX[75]  = 10;
#$greyDX[76]  = 10;
#$greyDX[77]  = 10;
#$greyDX[78]  = 10;
#$greyDX[79]  = 10;
#$greyDX[80]  = 10;
#$greyDX[81]  = 10;
#$greyDX[82]  = 10;
#$greyDX[83]  = 10;
#$greyDX[84]  = 10;
#$greyDX[85]  = 10;
#$greyDX[86]  = 10;
#$greyDX[87]  = 10;
#$greyDX[88]  = 10;
#$greyDX[89]  = 10;
#$greyDX[90]  = 10;
#$greyDX[91]  = 10;
#$greyDX[92]  = 10;
#$greyDX[93]  = 10;
#$greyDX[94]  = 10;
#$greyDX[95]  = 10;
#$greyDX[96]  = 10;
#$greyDX[97]  = 10;
#$greyDX[98]  = 10;
#$greyDX[99]  = 10;
#$greyDX[100] = 10;
#$greyDX[101] = 10;
#$greyDX[102] = 10;
#$greyDX[103] = 10;


# Define subroutines that modify the default values of the database.

sub applyBondFactor
{
   my $bondingFactor = $_[0];

   # Multiply all the covalent radii by the bonding factor.
   if ($bondingFactor != 1.0)
   {
      foreach my $atom (1..$#covalRadii)
         {$covalRadii[$atom] *= $bondingFactor;}
   }
}


# Define functions that return references to the requested data structures or
#   particular values.

sub getNumElements
   {return $numElements;}

sub getMaxQN_l
   {return $maxQN_l;}

sub getCoreOrbitals
   {return @{$coreOrbitals[$_[0]]};}

sub getValeOrbitals
   {return @{$valeOrbitals[$_[0]][$_[1]]};}  # $_[0]=basis; $_[1]=element

sub getOrbitalTerms
   {return $orbitalTerms[$_[0]][$_[1]];} #$_[0]=orbital type; $_[1]=element;

sub getNumTermsWF
   {return @{$numTermsWF[$_[0]]};}

sub getMinTermWF
   {return $minTermWF[$_[0]];}

sub getMaxTermWF
   {return $maxTermWF[$_[0]];}

sub getNumTermsPot
   {return $numTermsPot[$_[0]];}

sub getMinTermPot
   {return $minTermPot[$_[0]];}

sub getMaxTermPot
   {return $maxTermPot[$_[0]];}

sub getCovalRadiiRef
   {return \@covalRadii;}

sub getElementNamesRef
   {return \@elementNames;}

sub getValeCharge
   {return @{$valeCharge[$_[0]]};}

sub getValeChargeRef
   {return \@valeCharge;}

sub getMinTermRef
   {return \@minTermWF;}

sub getMaxTermRef
   {return \@maxTermWF;}

sub getColorDXRef
   {return \@colorDX;}

sub getGreyDXRef
   {return \@greyDX;}

sub getElementZ
{
   # Define passed variables.
   my $elementName = $_[0];

   # Declare local variables.
   my $element;

   foreach $element (1..$numElements)
   {
      if (lc($elementName) eq $elementNames[$element])
         {return $element;}
   }
}

1;
__END__
# Below is stub documentation for this module.

=head1 NAME

ElementData - Perl extension for storing information from the periodic table of
the elements and other data.

=head1 SYNOPSIS

  use ElementData;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for ElementData, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Paul Rulis/Ph/Ching, E<lt>rulis@activestate.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2007 by Paul Rulis/Ph/Ching

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
