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
//MODIFICATIONS and most documentation by Paul Rulis: rulisp@umkc.edu

typedef struct {
  double  x;        // Vx
  double  y;        // Vy
  double  z;        // Vz
} Vector;

typedef struct {
  int     noppf;    // # of points per face
  int     p[50];    // points number
  Vector  center;   // center (cg) of face
  Vector  normal;   // normal of face
} Face;

typedef struct {
  int     nof;      // # of faces
  int     noa;      // # of atoms/particles in this polyhedron
  int     nop;      // # of vertices associated with this Polyhedron
  int     start;    // first particle # in
  int     pr_count; // # of periodic adjacent polyhedrons
  int     blank;    // 1 = blank, 0 = filled
  int     ngh[25];  // neighbor polyhedron number
  int     p[100];   // point number
  Vector  *p0;      // pointer pointing to the first point
  Vector  c[30];    // centers of faces
  Vector  n[30];    // normals of faces
  Vector  periodic[8];  // reference point for filling adjacent polyhedra
  Vector  com;      // center of mass // NOT REALLY A CENTER OF MASS.  This is
                    //   really just the seed location.
  Vector  min;      // absolute max location, change with orientation
  Vector  max;      // absolute min location, change with orientation
} Polyhedron;

typedef struct {
  int     ph;       // Polyhedron No
  int     or;       // Orientation
  int     nature;   // Atom Nature (1=in grain boundary; 2=in grain bulk)
  int     atomTypeID; // This is the index number for this atom in the crystal
                      // typeZ array.
  double  covRad;   // Covalent radius
  Vector  r;        // Location vector
} Particle;

typedef struct {
  double rotM[3][3];  // rotation matrix
  int    or;          // orientation number
} Orientation;

struct node {         // search box link-list
  int           no;    
  struct node*  next;
};

typedef struct {
  char **elementName;
  double *covRad; // Covalent radius.
} ElementList;

typedef struct {

  // A given crystal will have some number of atoms.  Each atom will have a
  //   type ID number and an atomic Z number. The atomic Z number corresponds
  //   to the index number in the ElementList data structure that identifies
  //   the name of the element from the Periodic Table and also the covalent
  //   radius of the atom. The type ID number corresponds to the index number
  //   of a short list of the atomic Z number *for each type*.
  // As an example of the last part, consider a crystal with 2 types of Ca,
  //   4 types of O, 1 type of H and 1 type of P.  (Hydroxyapatite). After
  //   reading all of the data, we need an array of the format:
  //   {Ca, Ca, O, O, O, O, H, P}. All of the atoms of type Ca2 will have the
  //   same index number and all of the atoms of type O3 will have the same
  //   index number while all the atoms of Ca1 will have some other index
  //   number and all of the atoms of type O2 will have yet another index
  //   number.
  // The purpose of this is to coordinate between the LAMMPS "in" file and the
  //   LAMMPS "data" file. The example system would have 8 atom types and a
  //   eam potential file would need to have 8 element names listed after it
  //   to indicate what element each type was. Similarly, the LAMMPS data file
  //   would have
  // Additionally, it is necessary to have a short list of the 
  //   
  // Lattice vectors
  Vector a;
  Vector b;
  Vector c;

  int numAtoms;  // The number of atoms in the crystal cell.
  int numTypes;  // The number of atomic types of all elements in the crystal.
  int *typeZ;    // Atomic Z number of each type in the crystal cell. This
                 // array has indices from 0 to numTypes-1.
  int *atomZ;    // The atomic Z number of each atom in the crystal cell. This
                 // array has indices from 0 to numAtoms-1.
  int *atomTypeID;  // The sequential type ID number of each atom in the cell.
                    // This array has indices from 0 to numAtoms-1.
  double *covRad; // The covalent radius of each atom in the crystal cell. This
                  // array has indices from 0 to numAtoms-1.
  Vector *atomCoords;
  
} Crystal;



#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
