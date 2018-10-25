package BondData;

use 5.008002;
use strict;
use warnings;
use Env;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration use BondData ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw();

our $VERSION = '0.01';


# Preloaded methods go here.


# Define the arrays that hold bond data.
my @hookeBondCoeffs;

# Define other important variables.
my $numHookeBonds;
my $initData = 0; # Set to 1 when the initialization portion is run.

sub initBondData
{
   # Define local variables.
   my @values;
   my $line;
   my $bond;

   # Check to see if the initialization was run already and skip the rest if so.
   if ($initData == 0)
      {$initData=1;}
   else
      {return;}

   # Open the bond data file.
   open (BDATA,"<$OLCAO_DATA/bonds.dat") ||
         die "Cannot open bonds.dat, aborting.\n";

   # Read the number of bonds.
   @values = &prepLine(\*BDATA,$line,'\s+');
   if ($values[0] ne "NUM_HOOKE_BONDS")
      {die "Expecting NUM_HOOKE_BONDS tag in $OLCAO_DATA/bonds.dat";}
   @values = &prepLine(\*BDATA,$line,'\s+');
   $numHookeBonds = $values[0];

   # Read the Hooke bond coefficients.
   @values = &prepLine(\*BDATA,$line,'\s+');
   if ($values[0] ne "HOOKE_BOND_COEFFS")
      {die "Expecting HOOKE_BOND_COEFFS tag in $OLCAO_DATA/bonds.dat";}
   foreach $bond (1..$numHookeBonds)
   {
      @values = &prepLine(\*BDATA,$line,'\s+');
      $hookeBondCoeffs[$bond][1] = $values[0]; # Z number of first element
      $hookeBondCoeffs[$bond][2] = $values[1]; # Z number of second element
      $hookeBondCoeffs[$bond][3] = $values[2]; # k spring constant
      $hookeBondCoeffs[$bond][4] = $values[3]; # Rest bond length in Angstroms.
   }

   close (BDATA);
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


# Define functions that return references to the requested data structures or
#   particular values.

sub getNumHookeBonds
   {return $numHookeBonds;}

sub getHookeBondCoeffsRef
   {return \@hookeBondCoeffs;}

1;
__END__
# Below is stub documentation for this module.

=head1 NAME

BondData - Perl extension for storing information about bonds.

=head1 SYNOPSIS

  use BondData;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for BondData, created by h2xs. It looks like the
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
