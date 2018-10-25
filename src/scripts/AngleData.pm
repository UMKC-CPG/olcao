package AngleData;

use 5.008002;
use strict;
use warnings;
use Env;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration use AngleData ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw();

our $VERSION = '0.01';


# Preloaded methods go here.


# Define the arrays that hold angle data.
my @hookeAngleCoeffs;

# Define other important variables.
my $numHookeAngles;
my $initData = 0; # Set to 1 when the initialization portion is run.

sub initAngleData
{
   # Define local variables.
   my @values;
   my $line;
   my $angle;

   # Check to see if the initialization was run already and skip the rest if so.
   if ($initData == 0)
      {$initData=1;}
   else
      {return;}

   # Open the angle data file.
   open (ADATA,"<$OLCAO_DATA/angles.dat") ||
         die "Cannot open angles.dat, aborting.\n";

   # Read the number of angles.
   @values = &prepLine(\*ADATA,$line,'\s+');
   if ($values[0] ne "NUM_HOOKE_ANGLES")
      {die "Expecting NUM_HOOKE_ANGLES tag in $OLCAO_DATA/angles.dat";}
   @values = &prepLine(\*ADATA,$line,'\s+');
   $numHookeAngles = $values[0];

   # Read the Hooke angle coefficients.
   @values = &prepLine(\*ADATA,$line,'\s+');
   if ($values[0] ne "HOOKE_ANGLE_COEFFS")
      {die "Expecting HOOKE_ANGLE_COEFFS tag in $OLCAO_DATA/angles.dat";}
   foreach $angle (1..$numHookeAngles)
   {
      @values = &prepLine(\*ADATA,$line,'\s+');
      $hookeAngleCoeffs[$angle][1] = $values[0]; # Z number of first element
      $hookeAngleCoeffs[$angle][2] = $values[1]; # Z number of vertex element
      $hookeAngleCoeffs[$angle][3] = $values[2]; # Z number of second element
      $hookeAngleCoeffs[$angle][4] = $values[3]; # k angle spring constant
      $hookeAngleCoeffs[$angle][5] = $values[4]; # Rest angle in degrees.
      $hookeAngleCoeffs[$angle][6] = $values[5]; # Tolerance for angle relax.
   }

   close (ADATA);
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

sub getNumHookeAngles
   {return $numHookeAngles;}

sub getHookeAngleCoeffs
   {return \@hookeAngleCoeffs;}

1;
__END__
# Below is stub documentation for this module.

=head1 NAME

AngleData - Perl extension for storing information about angles.

=head1 SYNOPSIS

  use AngleData;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for AngleData, created by h2xs. It looks like the
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
