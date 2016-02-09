#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//ORIGINAL PROGRAM
////////////////////////////////////////////////////////////////////////////////
//    fillatoms -- C code for atomistic polycrystalline structure generation. //
//    Copyright (C) 2008  Sirirat Traiviratana                                //
//    http://code.google.com/p/fillatoms/                                     //
//                                                                            //
//    This program is free software: you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program.  If not, see http://www.gnu.org/licenses/.     //
////////////////////////////////////////////////////////////////////////////////
// MODIFICATIONS and most documentation by Paul Rulis: rulisp@umkc.edu
// I renamed it to polycryst because I thought that was a bit more informative.
// Copyright (C) 2015 Paul Rulis
// This program is included as part of the OLCAO distribution available at
//    https://github.com/UMKC-CPG.
// I'm not sure how the licensing should work out, but this is a distinct bit
//   of code that is not linked with any OLCAO code. So OLCAO is licensed under
//   the Educational Community License v2 and this will remain licensed under
//   the GNU Public License version 3 as described above.
#include "polycryst.h"
#include "global.h"

/*
USAGE:  ./polycryst inputFile


The structure of the input file is as follows:
Line 1: 6 real numbers.  The first three are the minimum x, y, z coordinates of
   the polycrystalline cell while the last three are the maximum x, y, z
   coordinates of the polycrystalline cell. The cell does not have to be
   centered on any specific point, the minimum does not need to be 0,0,0, and
   the cell does not have to be cubic. (Although, by definition it will be
   orthorhombic.) The coordinates may be selected for any reason of
   convenience. NOTE: All of these numbers are assumed to be in units of
   Angstroms regardless of the units for the crystalline cell.
Line 2: 1 integer, 1 real number, 2 integers, 1 real number. The first number
   is the integer number of seed crystal points that you would like to have in
   the polycrystalline cell. The second is a penalty (between 0.0 and 1.0) that
   is used to scale the covalent radius value of each atom. The covalent radius
   is used to control the elimination of atoms that are in the grain boundary
   but which are too close to each other. If two atoms are closer than the sum
   of the scaled covalent radii then one will be subject to elimination. The
   next integer is a number for initializing a pseudorandom number generator.
   If 0 is given, then the time function will be used. The next integer is a
   1 to indicate that the information for the crystalline cell is in units of
   Angstroms and it is a 0 to indicate that it is in units of Bohr radii. The
   final real number is a proximity control. It is always in units of Angstroms
   and it will be used to determine whether or not an atom is in a grain
   boundary or not.
Line 3: The x, y, z components of the a lattice vector. Must be x, 0, 0 in
   units given by the angst variable in Line 2.
Line 4: The x, y, z components of the b lattice vector. Must be 0, y, 0 in
   units given by the angst variable in Line 2.
Line 5: The x, y, z components of the c lattice vector. Must be 0, 0, z in
   units given by the angst variable in Line 2.
Line 6: 1 integer. This is the number of atoms in the crystalline cell with no
   further symmetry operations to be applied. (I.e. all atoms are specified as
   if the crystal is in spacegroup 1.)
Remaining Lines: 2 integers, 3 real numbers, 1 string. The first integer is
   simply an index number that counts the atom number. The second integer is
   the type designation for the atom on this line. The next three real numbers
   are the Cartesian coordinates of the atom in units specific by the angst
   variable given on Line 2. The final string is the abbreviation for the atom
   as from the Periodic Table of the Elements. Note that presently the name
   needs to be in all lower-case for correct matching to the data in the
   elements.dat file.
*/

/* Here are some examples:
Example 1: BCC Fe:
0.0 0.0 0.0 70.0 70.0 70.0
14 0.8 0 0 0.3
      5.41689986      0.00000000      0.00000000
      0.00000000      5.41689986      0.00000000
      0.00000000      0.00000000      5.41689986
2
    1     1         0.00000000         0.00000000         0.00000000 fe
    2     1         2.70844993         2.70844993         2.70844993 fe

Example 2: FCC Fe
0.0 0.0 0.0 70.0 70.0 70.0
14 0.8 0 0 0.3
      6.78978587      0.00000000      0.00000000
      0.00000000      6.78978587      0.00000000
      0.00000000      0.00000000      6.78978587
4
    1     1         0.00000000         0.00000000         0.00000000 fe
    2     1         0.00000000         3.39489294         3.39489294 fe
    3     1         3.39489294         0.00000000         3.39489294 fe
    4     1         3.39489294         3.39489294         0.00000000 fe

Example 3: HCP Ni
0.0 0.0 0.0 70.0 70.0 70.0
10 0.8 0 0 0.3
      5.00777416      0.00000000      0.00000000
      0.00000000      8.67371996      0.00000000
      0.00000000      0.00000000      8.18251401
2
    1     1         0.00000000         2.89123979         2.04562850 ni
    2     1         2.50388708         7.22809942         2.04562850 ni

Example 4: Fluorapatite
0.0 0.0 0.0 100.0 100.0 100.0
10 0.8 0 0 0.3
     17.66421364      0.00000000      0.00000000
      0.00000000     30.59531376      0.00000000
      0.00000000      0.00000000     12.97221300
84
    1     1        17.66421364        10.19842772         0.01426943 ca
    2     1         8.83210682        25.49608460         0.01426943 ca
    3     1         8.83210682         5.09922916         6.50037593 ca
    4     1         0.00000000        20.39688604         6.50037593 ca
    5     1         8.83210682         5.09922916        12.95794356 ca
    6     1         0.00000000        20.39688604        12.95794356 ca
    7     1        17.66421364        10.19842772         6.47183706 ca
    8     1         8.83210682        25.49608460         6.47183706 ca
    9     2        13.15453990        15.12326359         3.24305325 ca
   10     2         4.32243308        30.42092047         3.24305325 ca
   11     2        13.34178057         0.17439329         9.72915975 ca
   12     2         4.50967374        15.47205017         9.72915975 ca
   13     2        15.65402613         3.83053328         3.24305325 ca
   14     2         6.82191931        19.12819016         3.24305325 ca
   15     2         6.51986126        11.64151689         3.24305325 ca
   16     2        15.35196808        26.93917377         3.24305325 ca
   17     2         2.31224557         3.65613999         9.72915975 ca
   18     2        11.14435239        18.95379688         9.72915975 ca
   19     2        10.84229434        11.46712360         9.72915975 ca
   20     2         2.01018751        26.76478048         9.72915975 ca
   21     3         3.73421476         5.66013305         3.24305325  p
   22     3        12.56632159        20.95778993         3.24305325  p
   23     3         5.09789206         9.63752384         9.72915975  p
   24     3        13.92999888        24.93518072         9.72915975  p
   25     3        10.89528698         0.40385814         3.24305325  p
   26     3         2.06318015        15.70151502         3.24305325  p
   27     3        11.86681873         9.23366569         3.24305325  p
   28     3         3.03471190        24.53132258         3.24305325  p
   29     3        14.62950174         6.06399119         9.72915975  p
   30     3         5.79739492        21.36164807         9.72915975  p
   31     3        15.60103349        14.89379874         9.72915975  p
   32     3         6.76892667        30.19145562         9.72915975  p
   33     4         1.47142900         7.42242312         3.24305325  o
   34     4        10.30353582        22.72008000         3.24305325  o
   35     4         7.36067783         7.87523376         9.72915975  o
   36     4        16.19278465        23.17289064         9.72915975  o
   37     4         1.66838498        12.86074014         3.24305325  o
   38     4        10.50049180        28.15839702         3.24305325  o
   39     4        14.52439967        10.31215050         3.24305325  o
   40     4         5.69229285        25.60980739         3.24305325  o
   41     4        11.97192080         4.98550638         9.72915975  o
   42     4         3.13981398        20.28316326         9.72915975  o
   43     4         7.16372184         2.43691674         9.72915975  o
   44     4        15.99582867        17.73457362         9.72915975  o
   45     5         6.29110969         7.18224991         3.24305325  o
   46     5        15.12321651        22.47990679         3.24305325  o
   47     5         2.54099713         8.11540698         9.72915975  o
   48     5        11.37310396        23.41306386         9.72915975  o
   49     5         8.29864757         1.85713555         3.24305325  o
   50     5        17.13075439        17.15479243         3.24305325  o
   51     5        11.90656321         6.25827143         3.24305325  o
   52     5         3.07445638        21.55592831         3.24305325  o
   53     5        14.58975726         9.03938545         9.72915975  o
   54     5         5.75765044        24.33704233         9.72915975  o
   55     5         0.53345925        13.44052134         9.72915975  o
   56     5         9.36556607        28.73817822         9.72915975  o
   57     6         3.71655055         3.90396204         0.91713546  o
   58     6        12.54865737        19.20161892         0.91713546  o
   59     6         5.11555627        11.39369485         7.40324196  o
   60     6        13.94766309        26.69135173         7.40324196  o
   61     6        12.42500788         1.26664599         0.91713546  o
   62     6         3.59290106        16.56430287         0.91713546  o
   63     6        10.35476204        10.12704886         0.91713546  o
   64     6         1.52265522        25.42470574         0.91713546  o
   65     6        16.14155843         5.17060803         7.40324196  o
   66     6         7.30945161        20.46826491         7.40324196  o
   67     6        14.07131259        14.03101089         7.40324196  o
   68     6         5.23920577        29.32866777         7.40324196  o
   69     6         5.11555627        11.39369485        12.05507754  o
   70     6        13.94766309        26.69135173        12.05507754  o
   71     6         3.71655055         3.90396204         5.56897104  o
   72     6        12.54865737        19.20161892         5.56897104  o
   73     6        14.07131259        14.03101089        12.05507754  o
   74     6         5.23920577        29.32866777        12.05507754  o
   75     6        16.14155843         5.17060803        12.05507754  o
   76     6         7.30945161        20.46826491        12.05507754  o
   77     6        10.35476204        10.12704886         5.56897104  o
   78     6         1.52265522        25.42470574         5.56897104  o
   79     6        12.42500788         1.26664599         5.56897104  o
   80     6         3.59290106        16.56430287         5.56897104  o
   81     7         0.00000000         0.00000000         3.24305325  f
   82     7         8.83210682        15.29765688         3.24305325  f
   83     7         0.00000000         0.00000000         9.72915975  f
   84     7         8.83210682        15.29765688         9.72915975  f

Example 5: Silicon
0.0 0.0 0.0 70.0 70.0 70.0
13 0.8 0 0 0.3
     10.26253553      0.00000000      0.00000000
      0.00000000     10.26253553      0.00000000
      0.00000000      0.00000000     10.26253553
8
    1     1       0.00000000         0.00000000         0.00000000 si
    2     1       2.56563388         2.56563388         2.56563388 si
    3     1       0.00000000         5.13126776         5.13126776 si
    4     1       2.56563388         7.69690165         7.69690165 si
    5     1       5.13126776         0.00000000         5.13126776 si
    6     1       7.69690165         2.56563388         7.69690165 si
    7     1       5.13126776         5.13126776         0.00000000 si
    8     1       7.69690165         7.69690165         2.56563388 si

Example 6: ZnS
0.0 0.0 0.0 70.0 70.0 70.0
10 0.8 0 0 0.3
     10.23192196      0.00000000      0.00000000
      0.00000000     10.23192196      0.00000000
      0.00000000      0.00000000     10.23192196
8
    1     1         0.00000000         0.00000000         0.00000000 zn
    2     1         0.00000000         5.11596098         5.11596098 zn
    3     1         5.11596098         0.00000000         5.11596098 zn
    4     1         5.11596098         5.11596098         0.00000000 zn
    5     2         2.55798049         2.55798049         2.55798049  s
    6     2         2.55798049         7.67394147         7.67394147  s
    7     2         7.67394147         2.55798049         7.67394147  s
    8     2         7.67394147         7.67394147         2.55798049  s
*/

// Some notes.
// I use Polygon, Polyhedron, and Voronoi region mostly interchangably. The
//   second term is more likely to be used if I am trying to refer specifically
//   to the data structure that holds Polyhedron information.  (A struct).

int main(int argc, char *argv[])
{
  int i,j,k,n; // Loop index variables
  int a; // This is a bit of a dummy variable that gets used for a variety of
         //   different things including: Dummy read variable and counter for
         //   the number of atoms (I think).
  int nos; // Number of seeds in the central cell. (AKA #sitecenters)
  int noj; // Number of Voronoi vertices in the central cell. (AKA #junctions)
  int noppp; // A dummy variable for reading in the number of Voronoi vertices
             //   for a given Polyhedron.
  int nopoly; // The number of Voronoi polygons in all 27 cells.
  int nop; // The number of Voronoi vertices in all 27 cells.
  int nof; // The number of Faces in all Voronoi Polygons.
  int noppf; // The number of vertices for a given Face.
  int pfcatom; // The number of atoms from a perfect crystal that would fit in
               //   the simulation cell given in the "box" input file.
  int extra; // A strang bit used to turn on some "debugging" print statements.
             //   For some reason it is hard coded in the middle of the program
             //   to be turned off at the moment.
  int totalAtomsPrinted;
  double lowerLimit = -1e100; // A max value can't be any lower than this.
  double upperLimit =  1e100; // A min value can't be any higher than this.
  double x,y,z,boxvolume,unitcellvolume;
  double penalty;  // Used to control the minimum distance that atoms from
        // two different grains are allowed to get.
  double pi=4.0*atan2(1.0,1.0);  // Definition of Pi.
  double offset[3] = {-1.0,0.0,1.0};  // Tool for making a 27 cell supercell.
  Vector *p; // An array holding Vectors to all the Voronoi vertices.
  Vector *site; // First holds the site center coordinates for each seed
                //   crystal in the central polycrystalline cell. Later holds
                //   the site center coordinates for all seed crystal centers
                //   including the periodic replicated ones.
  int *sitecenter; // Holds index numbers of sites that are inside the central
                   // cell.
  int *junction; // Similar to sitecenter except it holds index numbers for
                 // polygon vertices that are inside the central cell.
  Orientation Or; // An orientation object.
  Polyhedron *ph; // Array of all Polyhedra including periodic ones.
  Face *f;  // Array of all faces including those in periodic cells.
  FILE *fi,*fj;  // Poorly named file pointers.
  int *filledlist; // An array to track which polygons have been filled.
  int fill=0;
  time_t start_time;
  unsigned int seed;
  ElementList eList;
  Crystal cryst;
  int angst;
  int dummy;
  char tempName[80];
  char line[80];
  Vector temp; // A temporary vector used for intermediate calculations.

  // Assign constants.
  bohr = 0.5291772180; // From NIST (m) (e-10).

  // Initialize the list of element names and covalent radii from a data file.
  fillElementList(&eList);

  printf("Start reading input files.\n");

  // Open the box input file for reading.
  fj = fopen(argv[1],"r");

  // Read the cartesian min and max values for the x, y, z axes.
  fgets(line,80,fj);
  sscanf(line,"%lf %lf %lf %lf %lf %lf",&gmx,&gmy,&gmz,&px,&py,&pz);
  printf("  X: %8.2f %8.2f\n  Y: %8.2f %8.2f\n  Z: %8.2f %8.2f\n",
        gmx,px,gmy,py,gmz,pz);

  // Read the number of seeds, the "penalty", a number that controls the random
  //   seeding of the cell, and a flag indicating whether or not the lattice
  //   parameters and atomic coordinates are in angstroms (1) or atomic units
  //   (2).
  fgets(line,80,fj);
  sscanf(line,"%d %lf %d %d %lf",&nos,&penalty,&seed,&angst,&proximity);
  printf("%d %lf %d %d %lf\n",nos,penalty,seed,angst,proximity);

  // Read the lattice parameters for the crystal
  fgets(line,80,fj);
  sscanf(line,"%lf %lf %lf",&cryst.a.x, &cryst.a.y, &cryst.a.z);
  printf ("%lf %lf %lf\n",cryst.a.x, cryst.a.y, cryst.a.z);

  fgets(line,80,fj);
  sscanf(line,"%lf %lf %lf",&cryst.b.x, &cryst.b.y, &cryst.b.z);
  printf ("%lf %lf %lf\n",cryst.b.x, cryst.b.y, cryst.b.z);

  fgets(line,80,fj);
  sscanf(line,"%lf %lf %lf",&cryst.c.x, &cryst.c.y, &cryst.c.z);
  printf ("%lf %lf %lf\n",cryst.c.x, cryst.c.y, cryst.c.z);

  // Convert lattice parameters to Angstroms if necessary.
  if (angst == 0)
  {
     cryst.a.x *= bohr; cryst.a.y *= bohr; cryst.a.z *= bohr;
     cryst.b.x *= bohr; cryst.b.y *= bohr; cryst.b.z *= bohr;
     cryst.c.x *= bohr; cryst.c.y *= bohr; cryst.c.z *= bohr;
  }

  // Read the number of atoms.
  fgets(line,80,fj);
  sscanf(line,"%d",&cryst.numAtoms);
  printf ("%d\n",cryst.numAtoms);

  // Allocate space to hold the atomic coordinates, Z number, covalent radius
  //   and type ID number (index to the typeZ array) of each atom. Additionally,
  //   make a list of the Z number of unique types in the same order that they
  //   are encountered. More space is allocated than will typically be used for
  //   typeZ, but we assume a possible worst case scenario where every atom is
  //   a unique element.
  cryst.atomCoords = (Vector *)malloc(cryst.numAtoms * sizeof(Vector));
  cryst.atomZ      = (int *)   malloc(cryst.numAtoms * sizeof(int));
  cryst.covRad     = (double *)malloc(cryst.numAtoms * sizeof(double));
  cryst.atomTypeID = (int *)   malloc(cryst.numAtoms * sizeof(int));
  cryst.typeZ      = (int *)   malloc(cryst.numAtoms * sizeof(int)); // This
        // array will only use numTypes array slots, but we allocate more than
        // typically used to deal with the worst case scenario.

  // Initialize the number of types to 0.
  cryst.numTypes = 0;

  // Read the Cartesian coordinates of each atom plus the element name. Then
  //   convert the coordinates to Angstroms if necessary and convert the name
  //   to an atomic Z number and store it.
  for (i=0;i<=cryst.numAtoms-1;i++)
  {
     // Read the data.
     fgets(line,80,fj);
     printf ("%s",line);
     sscanf(line,"%d %d %lf %lf %lf %s",&dummy,&cryst.atomTypeID[i],
           &cryst.atomCoords[i].x,&cryst.atomCoords[i].y,
           &cryst.atomCoords[i].z,tempName);

     // Convert the coordinates to Angstroms if necessary.
     if (angst == 0)
     {
        cryst.atomCoords[i].x *= bohr;
        cryst.atomCoords[i].y *= bohr;
        cryst.atomCoords[i].z *= bohr;
     }

     // Map the element name to an atomic Z number and the covalent radius.
     getElementZandCovRad(tempName,&eList,&cryst.atomZ[i],&cryst.covRad[i],
           &penalty);

     // If this atom is the first atom of this element, then we store this
     //   element as an addition to a list of unique elements in the system.
     addToTypeList(&i,&cryst);
  }

  // Close the input file.
  fclose(fj);

  // This is also a bit odd. The seed variable is used to select a seed for the
  //   pseudorandom number generator. If the seed is zero then we use the time
  //   to get a seed, if the seed is 1 we use the number 1 for the seed. If the
  //   seed is some other number then we use that number for the seed. Thus,
  //   the use of seed==1 is just a subset of the last option.  (Note that the
  //   word "seed" here is not in reference to the concept of a seed crystal
  //   for grain growth. It is a seed number for pseudorandom number
  //   generation.)
  if (seed==0) {seed=time(0);srand(seed);printf("Time seed\n");}
  else if (seed==1) {srand(1);printf("Default seed\n");}
  else {srand(seed);printf("User seed\n");}

  printf("Seed for srand call was %u.\n",seed);

  // Create an array of Vectors to store the positions of the seed points for
  //   crystalline grain growth.
  site = (Vector *)malloc(nos*sizeof(Vector));

  // Create two files that will contain the locations of the seed crystals in
  //   the cell and the neighboring 26 replicated cells surrounding it. The two
  //   files are identical except that siteq is used as input for qhull and
  //   thus needs to have two additional header lines indicating the dimension
  //   of the problem and the number of points.
  fi = fopen("siteq","w");
  fj = fopen("site","w");
  fprintf(fi,"3\n");  // Indicate the dimension.
  fprintf(fi,"%d\n",nos*27);  // Indicate the number of seed crystal locations
                              // in the central cell and 26 neighboring cells.

  // Create the locations of the seed crystals in the central cell.
  for (i=0;i<nos;i++) {
    site[i].x = RandomWithIn(gmx,px);
    site[i].y = RandomWithIn(gmy,py);
    site[i].z = RandomWithIn(gmz,pz);
  }

  // Sorry, but this is a bit annoying.  The original variable that contained
  //   the upper bound for the cell size is converted into a variable that
  //   contains the length of the cell in each x,y,z dimension.
  px = px - gmx; py = py - gmy; pz = pz - gmz;

  // Compute the volume of the cell.
  boxvolume = px*py*pz;

  // Compute the volume of the crystal lattice.
  temp = cross(&cryst.b,&cryst.c);
  unitcellvolume = dot(&cryst.a,&temp);
  if (unitcellvolume < 0.0)
     {unitcellvolume *= -1.0;}
  printf("Unit cell volume = %lf\n",unitcellvolume);

  // This output should be consolidated.
  printf("%d grains with ",nos);
  printf("grain size = approx. %f Angstrom.\n",
      2.0*pow((3.0*boxvolume/4.0/(double)nos/pi),(1.0/3.0)));

  // Record the positions of the seed crystal sites in all 27 cells (in two
  //   almost identical files). The reason for doing this appears to be to
  //   allow for direct reading of the seed site locations in the future. I'm
  //   not sure why the information in the second file (site) couldn't just be
  //   stored in an array of Vectors.
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
        for (n=0;n<nos;n++) {
          fprintf(fi,"%15.10f %15.10f %15.10f\n",
              site[n].x+offset[i]*px,
              site[n].y+offset[j]*py,
              site[n].z+offset[k]*pz);
          fprintf(fj,"%15.10f %15.10f %15.10f\n",
              site[n].x+offset[i]*px,
              site[n].y+offset[j]*py,
              site[n].z+offset[k]*pz);
        }
      }
    }
  }
  fclose(fi);
  fclose(fj);

  // Free the site array.  It will be reallocated later to hold 


  // Call the qhull program, feeding it the siteq file as input. The siteq file
  //   indicates that this is a three dimensional problem and the number of
  //   seed points in the block of 27.
  // The command line parameters have the following meaning and impact:
  //   v    - Voronoi diagram as the dual of the Delaunay triangulation
  //   o    - OFF format (if 'v', outputs Voronoi regions)
  //   s    - summary of results (default)
  //   Fv   - count plus vertices for each facet
  //          for 'v', Voronoi diagram as Voronoi vertices for pairs of sites
  //   TO <file> - output results to <file>, may be enclosed in single quotes
  system("cat siteq | qhull v o s Fv TO input");

  // This call to qhull produces a file called "input" that contains:
  //   The first line is the dimension.
  //   The second line is the number of Voronoi vertices, the number of Voronoi
  //      regions, and (the qhull documentation says) the number of ridges, but
  //      I think that this number is something else that isn't used by this
  //      code.
  //   The first vertex listed is the infinity vertex [-10.101 -10.101 -10.101]
  //   All other vertices are simply 3D vectors from the origin to the points
  //      in the 27 cell space where the Voronoi vertices are to be found.
  //   Then, the Voronoi region associated with each seed crystal point in the
  //      27 cell space is defined by listing the index numbers of the vertices
  //      associated with that region. On each line the first number is the
  //      number of vertices and the numbers that follow are the indices of the
  //      Voronoi vertices.
  //   The next line contains a single integer for the number of facets.
  //   The next set of lines are hard to fathom from the Qhull documentation.
  //      It says essentially the following:
  //   Then Qhull prints the index numbers for the vertices of each facet on a
  //      separate line. The first number is the number of vertices for the
  //      current facet. The remaining numbers are the indices of the
  //      corresponding points.
  //   However, I believe the truth to be closer to:
  //   Then Qhull prints a line with a series of numbers. The first number is
  //      the number of other numbers on the line. The next two numbers are
  //      index numbers indicating which two Voronoi regions the face
  //      separates.  The remaining numbers are index numbers for the vertices
  //      associated with this face. Thus the first number is equal to the
  //      number of vertices for this face + 2.  Go figure.

  // Open the output of qhull for reading.  Additionally, open the previously
  //   written data that contains all of the seed crystal points.
  fi = fopen("input","r");
  fj = fopen("site","r");

  // Read in the number of dimensions (into a dummy variable), the number of
  //   Voronoi vertices (nop), the number of Voronoi regions (polygons), and
  //   some random bit which basically gets ignored but then treated as a
  //   "verbose" flag (which makes almost no sense).
  fscanf(fi,"%d %d %d %d",&a,&nop,&nopoly,&extra);extra=0; //FIXME extra=verbose

  // This is a bit of verbosity.
  if (extra) printf("Nop     :%d\n",nop);

  // Create an array of Vectors to store the Voronoi vertices.  There is a
  //   really annoying issue here.  The name of this array of Vectors "p" is
  //   the same as the name of an array of integers inside the Face and
  //   Polyhedron structs. This "p" holds Vectors to the Voronoi vertices.
  //   The other "p"s hold integers that are index numbers for vertices
  //   associated with the Voronoi Polyhedron and Faces.
  p = (Vector *)malloc(nop*sizeof(Vector));

  // Create an array of Polyhedron structs to store the Voronoi polyhedra.
  ph = (Polyhedron *)malloc(nopoly*sizeof(Polyhedron));

  // Create an array of Vectors to store the centers of the Polyhedron structs.
  //   Note that this will just hold the same sites that were written out to
  //   the "site" file earlier.  (Also annoying, this variable "site" was used
  //   earlier to hold specifically just the positions of the original randomly
  //   determined seed points inside the central cell.  This array is capable
  //   of holding those and the other replicated seed points in the neighboring
  //   cells. Later on, an array "sitecenter" will hold the index numbers of
  //   the "site" Vectors that are inside the central cell.
  site = (Vector *)malloc(nopoly*sizeof(Vector));

  // Create an array of integers that will act as indices. The index number of
  //   the sitecenter array is an index for the number of seeds. The value in
  //   the sitecenter array at a given index is itself an index number that
  //   points to the particular seed crystal sites in the site array that are
  //   inside the central cell.  
  sitecenter = (int *)malloc(nopoly*sizeof(int));

  // This array serves much the same purpose as the sitecenter array except
  //   that it indicates which Voronoi vertices are inside the central cell
  //   out of all the vertices.
  junction = (int *)malloc(nop*sizeof(int));

  // For each crystal seed location in the 27 cells we will determine if the
  //   seed location is actually inside the central cell.  (Note that the
  //   loop limit of nopoly is misleading and could just as easily (and more
  //   clearly) have been a different variable defined much earlier on such as
  //   numSeedCrystals.) The value stored in sitecenter is just an index number
  //   for the site. Note that points outside the central cell will not have an
  //   index value assigned to them and so the array indices may be undefined
  //   depending on the compiler policies. This is dangerous and it should be
  //   corrected. Final note: The data that is read in here is just the site
  //   data that was written before.
  nos = 0;
  for (i=0;i<nopoly;i++) {
    fscanf(fj,"%lf %lf %lf",&x,&y,&z);
    site[i].x = x;
    site[i].y = y;
    site[i].z = z;
    if (inthebox(&site[i])) { sitecenter[nos] = i; nos++; }
  }
  fclose(fj);
  
  // Compute the maximum possible number of atoms that could conceivably ever
  //   appear in the cell.  The leading 4 comes from an FCC full unit cell. The
  //   trailing 0.01 appears to attempt to account for the slightly higher
  //   density that is possible in a grain boundary vs. a bulk grain. This is a
  //   slightly dangerous calculation and it should be improved.
//  maxatom = (int)(boxvolume / unitcellvolume * 4.01);
  maxatom = (int)(boxvolume / unitcellvolume * ((double)cryst.numAtoms+0.05));
        // Changed to (hopefully) reduce the liklihood of an attepted access
        // outside of array bounds.

  // Compute the number of atoms that would be present in a perfect crystal.
  pfcatom = (int)(boxvolume / unitcellvolume * cryst.numAtoms);

  // Allocate space to hold the an array of "atom" particles. Each atom is a
  //   struct that holds the associated polyhedron index number, the type # of
  //   the atom, a tag indicating whether the atom is in the bulk of the grain
  //   or in the boundary region of the grain, an ID # indicating the element
  //   that the atom is, the orientation of the grain, and a Vector for
  //   the location of the atom.
  atom = (Particle *)malloc(maxatom*sizeof(Particle));
  printf("Allocated %ld particles (perfect crystal %d particles).\n",
        maxatom,pfcatom);

  // For each Voronoi vertex, determin if it is inside the central cell. If it
  //   is, then record the index number of the vertex in the junction array.
  //   Note that as with the sitecenter array above, the junction array
  //   elements will only have a value assigned if the vertex is inside the
  //   central cell. This will also store all of the Voronoi vertex points
  //   from all 27 cells. (Only the central cell vertices will be referenced
  //   by the junction array.) Note that the first (i-0) vertex is the
  //   so-called infinity vertex at -10.101, -10.101, -10.101.
  noj = 0;  // Counter for the number of junctions.
  for (i=0;i<nop;i++) {
    fscanf(fi,"%lf %lf %lf",&x,&y,&z); // Read in the Voronoi vertex.
    p[i].x = x;
    p[i].y = y;
    p[i].z = z;
    if (inthebox(&p[i])) { junction[noj] = i; noj++; }
    if (extra) 
      printf("P:%d %14.10f %14.10f %14.10f\n",i,p[i].x,p[i].y,p[i].z);
  }

  if (extra) printf("Nopoly  :%d\n",nopoly);


  // For each polyhedron, read in the index number of the vertices associated
  //   with it.  Additionally, mark the polyhedron as being empty.  (I assume
  //   that we will only fill the polyhedra that have a component in the
  //   central cell.)
  for (j=0;j<nopoly;j++) {
    ph[j].blank = 1; // Initialize to the assumption of an empty Voronoi region.
    ph[j].p0 = &p[0]; // Every Polyhedron includes a ref to the infinity vertex.
    fscanf(fi,"%d",&noppp); // Read in the number of vertices for this Polyhed.
    if (extra) printf("Noppp   :%d\n",noppp);
    ph[j].nop = noppp;  // Store the number of vertices for this Polyhedron.

    // Read in the index numbers of the Voronoi vertices for this Polyhedron.
    //   Note that the p[] array has a hard limit for its dimension of 100.
    //   That seems reasonable, but it isn't the best way to do things.
    for (i=0;i<noppp;i++) {
      fscanf(fi,"%d",&a);
      ph[j].p[i] = a; if (extra) printf("%d ",a);
    }
  }

  // Read in the number of faces (facets) and allocate space to hold each Face.
  fscanf(fi,"%d",&nof);
  f = (Face *)malloc(nof*sizeof(Face));

  // For each face, read the number of vertices for this face and then read
  //   the index numbers for the vertices.  Take special note of the
  //   documentation above regarding the meaning of the numbers on each line.
  for (i=0;i<nof;i++) {
    if (extra) printf("Face:%d\n",i);
    fscanf(fi,"%d",&noppf);  // Get the number of vertices for this face. (+2)
    if (extra) printf("%d ",noppf);
    f[i].noppf = noppf;  // Store the number of vertices for this face. (+2)

    // Read in the index number for each vertex associated with this Face. Note
    //   that this p[] array has a hard limit for its dimension of 50. Same
    //   comment as above.  Also note that there is a good chance that the
    //   first two numbers are not vertex indices but are instead Voronoi
    //   region indices.
    for (k=0;k<noppf;k++) {
      fscanf(fi,"%d",&a);
      if (extra) printf("%d ",a);
      f[i].p[k] = a;
    }
    if (extra) printf("\n");
  }
  if (extra) printf("\n");
  fclose(fi);
  
  printf("\nFinish reading, counting neighbors of each polyhedron.\n");

  // Now, for each polygon in the 27 cell system we loop over the list of all
  //   Faces in the system and determine which Faces are associated with this
  //   Polyhedron.
  for (j=0;j<nopoly;j++) {
    ph[j].nof = 0;  // Initialize the number of faces for this Polyhedron to 0.
    ph[j].com = site[j];  // Record the seed point for this Polyhedron.

    // Establish initial upper and lower limits for the extent of the
    //   Polyhedron. These will be computed exactly later, but the initial
    //   values need to be larger than the "largest possible" min value for the
    //   min Vector and smaller than the "smallest possible" max value for the
    //   max Vector.
    ph[j].min.x = upperLimit;
    ph[j].min.y = upperLimit;
    ph[j].min.z = upperLimit;
    ph[j].max.x = lowerLimit;
    ph[j].max.y = lowerLimit;
    ph[j].max.z = lowerLimit;

    // At the end of this loop the current j Polyhedron will have a list of
    //   Faces associated with it and it will have stored the index number of
    //   the Polyhedron on the other side of each Face as well as the unit
    //   normal vector and the center point for the Face.
    for (i=0;i<nof;i++) {

      // Check if the first or second number in the list of vertices for this
      //   Face are the same as the current Polygon index number and thus that
      //   the current Face is associated with the current Polygon. (Recall
      //   that the first two numbers in this list are Polyhedron index numbers
      //   and not vertex index numbers.)
      if ((f[i].p[0]==j)||(f[i].p[1]==j)) {
        if (extra) {
          printf("Polyhedron:%d has Face:%d -> ",j,i);
          for (k=0;k<f[i].noppf;k++) printf("%d ",f[i].p[k]);
          printf("\n");
        }

        // If the first index matches the index of the current Polyhedron.
        if (f[i].p[0]==j) {
          // Compute a vector pointing from the seed location at index 0 to
          //   the seed location at index 1.
          temp = va2b(&site[f[i].p[0]],&site[f[i].p[1]]);

          // Compute a unit vector from the seed at 0 toward the seed at 1
          //   which will become the unit normal for this Face within this
          //   Polygon.
          ph[j].n[ph[j].nof] = unit(&temp);

          // Similarly, obtain the identification of the Polyhedron neighbor
          //   number on the other side of this Face for this j Polyhedron.
          ph[j].ngh[ph[j].nof] = f[i].p[1]; 
        }
        // If the second index matches the index of the current Polyhedron then
        //   do all the same stuff just the other way around.
        else {
          temp = va2b(&site[f[i].p[1]],&site[f[i].p[0]]);
          ph[j].n[ph[j].nof] = unit(&temp);
          ph[j].ngh[ph[j].nof] = f[i].p[0]; 
        }

        // The center point of the face can be computed as the midpoint of the
        //   vector obtained from the sum of the two seed location vectors.
        ph[j].c[ph[j].nof] = center2(&site[f[i].p[1]],&site[f[i].p[0]]);

        // Because a Face was found we increment the number of Faces for this
        //   Polyhedron.
        ph[j].nof++;
      }
    }
    if (extra) printf("Nof:%d \n",ph[j].nof);
  }

  printf("\nSearch Box generating.\n");

  // Allocate two link-list heads for *each* Polyhedron. The first will hold
  //   a list of all the atoms for this Polyhedron, and the second will hold a
  //   subset list of all the atoms. This subset is the set of atoms that are
  //   "close" to the edge of the Polyhedron and which may be considered as
  //   being in the grain boundary.
  searchbox = (struct node**)malloc(nopoly*sizeof(struct node*));
  gb        = (struct node**)malloc(nopoly*sizeof(struct node*));
  printf("Allocated %d searchboxes(polyhedron).\n",nopoly);

  // Make the linked list heads point to NULL.
  for (j=0;j<nopoly;j++) searchbox[j] = gb[j] = NULL;


  printf("\nStart filling up polyhedrons.\n");

  // Initialize a running count of the number of atoms.
  ano=0;

  // Allocate an array to track which Polygons have been filled. (Add more
  //   detail later.)
  filledlist = (int *)malloc(nopoly*sizeof(int));

  time(&start_time);

  // Loop over the all of the Polyhedra in the 27 cell system and for each
  //   Polyhedron that is in the central cell and which has not yet been
  //   filled, determine the orientation of the grain and fill it and all
  //   periodic adjacent Polyhedra.
  for (i=0;i<nopoly;i++) {

    // Check if the current Polyhedron is both in the central cell and not yet
    //   filled. If it is either outside the central cell or already filled,
    //   then we ignore it. However, this does not mean that a Polyhedron that
    //   outside the central cell will never be filled.  Instead, during the
    //   filling of a Polyhedron a check is made to identify which other
    //   Polyhedra that are periodic copies of the central one need to be
    //   filled.
    if ((inthebox(&site[i])) && (ph[i].blank)) {

      // Establish the orientation of the grain for the current Polyderon.
      Or = MakeOrientation(0);

      //FILL POLYGONS
      i = FillAt(&ph[0],nopoly,i,&Or,&ph[i].com,&cryst);
      for (j=0;j<ph[i].pr_count;j++) {
        k = FillAt(&ph[0],nopoly,0,&Or,&ph[i].periodic[j],&cryst);
      }
      printpercent(start_time,100.*(double)ano/(double)maxatom);
    } 
  }

  // Get the number of filled Polyhedra and also create an array (of length
  //   equal to the number of filled Polyhedra) where each element in the array
  //   holds the index number of a filled Polyhedron.
  fill = UpdateFilledList(&ph[0],nopoly,&filledlist[0]);
  printf("Filled %ld particles (%ld different)\n",ano,ano-pfcatom);

  
  printf("\nFilter atoms, too close within neighbors\n");

  // Go through the list of filled Polyhedra. For each one, consider all of its
  //   Faces and initiate a search of the pairs of atoms from Polyhedra on
  //   opposite sides of the Face. Obtain a running count of the number of
  //   atoms that are identified as being too close (n).
  time(&start_time);n = 0;
  for (i=0;i<fill;i++) {
    for (j=0;j<ph[filledlist[i]].nof;j++) {
      a = Search2Boxes(filledlist[i],ph[filledlist[i]].ngh[j]);
      n = n + a;
    }
  }
  printf("Total %d to be deleted.\n",n);
  printpercent(start_time,100.0); 
  
  
  // Obtain the total number of atoms in the box.
  for (i=0;i<fill;i++)
    {a = countAtomsInBox(filledlist[i],a);}
  totalAtomsPrinted = a-1;

  printf("\nPrint out all boxes (and their atoms) into one LAMMPS file.\n");

  printLammpsIn(&cryst,&eList);

  printCellHeader(totalAtomsPrinted,nos,penalty,seed,&cryst,&eList);

  time(&start_time);a = 1;
  for (i=0;i<fill;i++) {
    a = printCell(filledlist[i],a);
    if (((i+1)*100/fill) % 10 == 0) 
      printpercent(start_time,100.*(double)(i+1)/(double)fill); 
  }

  //Probably need to free a lot more things.
  free(atom);       atom = NULL;
  free(searchbox);  searchbox = NULL;
  free(gb);         gb = NULL;
  
  return 0;
} 

// END MAIN ################################################### END MAIN
