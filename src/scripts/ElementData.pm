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
my @atomicMasses;    # Atomic mass of each element in atomic mass units.
my @covalRadii;      # Covalent radii of each element in Angstroms.
my @numUJElectrons;  # Number of electrons in the highest d or f orbital.
my @ljPairCoeffs;    # LJ pair coefficients for LAMMPS MD simulations.
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
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "NUM_ELEMENTS")
      {die "Expecting NUM_ELEMENTS tag in $OLCAO_DATA/elements.dat";}
   @values = &prepLine(\*EDATA,$line,'\s+');
   $numElements = $values[0];

   # Read the max l quantum number of all atoms.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "MAX_QN_L")
      {die "Expecting MAX_QN_L tag in $OLCAO_DATA/elements.dat";}
   @values = &prepLine(\*EDATA,$line,'\s+');
   $maxQN_l = $values[0];

   # Read the element names.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "ELEMENT_NAMES")
      {die "Expecting ELEMENT_NAMES tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $elementNames[$element] = $values[0];
   }

   # Read the atomic masses.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "MASS")
      {die "Expecting MASS tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $atomicMasses[$element] = $values[0];
   }

   # Read the covalent radii.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "COVALENT_RADII")
      {die "Expecting COVALENT_RADII tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $covalRadii[$element] = $values[0];
   }

   # Read the number of electrons in the highest d or f orbital in the ground
   #   state. Used for UJ calculations.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "NUM_UJ_ELECTRONS")
      {die "Expecting NUM_UJ_ELECTRONS tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $numUJElectrons[$element] = $values[0];
   }

   # Read the LJ potential coefficients for each atom.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "LJ_PAIR_COEFFS")
      {die "Expecting LJ_PAIR_COEFFS tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $ljPairCoeffs[$element][1] = $values[0];
      $ljPairCoeffs[$element][2] = $values[1];
   }

   # Read the designation of core orbitals for each atom.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "CORE_ORBITALS")
      {die "Expecting CORE_ORBITALS tag in $OLCAO_DATA/elements.dat";}
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
      @values = &prepLine(\*EDATA,$line,'\s+');
      if (($basis == 1) and ($values[0] ne "MB_BEYOND_CORE"))
         {die "Expecting MB_BEYOND_CORE tag in $OLCAO_DATA/elements.dat";}
      elsif (($basis == 2) and ($values[0] ne "FB_BEYOND_MB"))
         {die "Expecting FB_BEYOND_MB tag in $OLCAO_DATA/elements.dat";}
      elsif (($basis == 3) and ($values[0] ne "EB_BEYOND_FB"))
         {die "Expecting EB_BEYOND_FB tag in $OLCAO_DATA/elements.dat";}
      foreach $element (1..$numElements)
      {
         @values = &prepLine(\*EDATA,$line,'\s+');
         foreach $QN_l (0..$maxQN_l)
            {$valeOrbitals[$basis][$element][$QN_l] = $values[$QN_l];}
      }
   }

   # Read the core charge for each spdf orbital of each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "CORE_CHARGE")
      {die "Expecting CORE_CHARGE tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      foreach $QN_l (0..$maxQN_l)
         {$coreCharge[$element][$QN_l] = $values[$QN_l];}
   }

   # Read the valence charge for each spdf orbital of each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "VALE_CHARGE")
      {die "Expecting VALE_CHARGE tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      foreach $QN_l (0..$maxQN_l)
         {$valeCharge[$element][$QN_l] = $values[$QN_l];}
   }

   # Read the number of radial wave function Gaussians for each spdf orbital
   #   of each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "NUM_RWF_TERMS_SPDF")
      {die "Expecting NUM_RWF_TERMS_SPDF tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      foreach $QN_l (0..$maxQN_l)
         {$numTermsWF[$element][$QN_l] = $values[$QN_l];}
   }

   # Read the minimum exponential alpha for the radial wave functions for
   #   each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "MIN_RWF_ALPHA")
      {die "Expecting MIN_RWF_ALPHA tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $minTermWF[$element] = $values[0];
   }

   # Read the maximum exponential alpha for the radial wave functions for
   #   each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "MAX_RWF_ALPHA")
      {die "Expecting MAX_RWF_ALPHA tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $maxTermWF[$element] = $values[0];
   }

   # Read the number of potential terms for each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "NUM_POT_TERMS")
      {die "Expecting NUM_POT_TERMS tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $numTermsPot[$element] = $values[0];
   }

   # Read the minimum exponential alpha for the potential function for
   #   each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "MIN_POT_ALPHA")
      {die "Expecting MIN_POT_ALPHA tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $minTermPot[$element] = $values[0];
   }

   # Read the maximum exponential alpha for the potential function for
   #   each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "MAX_POT_ALPHA")
      {die "Expecting MAX_POT_ALPHA tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $maxTermPot[$element] = $values[0];
   }

   # Read data strings for each spdf orbital type of each element. The data
   #   strings are boolean identifiers for which Gaussians in the orbital
   #   expansion should be used.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "GAUSS_TERMS_SPDF")
      {die "Expecting GAUSS_TERMS_SPDF tag in $OLCAO_DATA/elements.dat";}
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
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "ODX_COLOR")
      {die "Expecting ODX_COLOR tag in $OLCAO_DATA/elements.dat";}
   foreach $element (1..$numElements)
   {
      @values = &prepLine(\*EDATA,$line,'\s+');
      $colorDX[$element] = $values[0];
   }

   # Read the openDX grey scale assignment for each element.
   @values = &prepLine(\*EDATA,$line,'\s+');
   if ($values[0] ne "ODX_GREY")
      {die "Expecting ODX_GREY tag in $OLCAO_DATA/elements.dat";}
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

sub getAtomicMassesRef
   {return \@atomicMasses;}

sub getCovalRadiiRef
   {return \@covalRadii;}

sub getNumUJElectrons
   {return \@numUJElectrons;}

sub getLJPairCoeffs
   {return \@ljPairCoeffs;}

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
