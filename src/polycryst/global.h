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
//#define penalty 0.5 // allow 50% of 2*r0
//#define lattice 3.615 // FCC (Cu) Copper Lattice
//#define lattice 3.52  // FCC (Ni) Nickel Lattice

Particle *atom;                 // pointer to atom
double px,py,pz,gmx,gmy,gmz;    // global domain definition
//double bx,by,bz;                // search box size
//int nx,ny,nz;                   // number of search box in each direction
long ano,maxatom;               // running atom # and max #
struct node** searchbox;        // pointer to search box
struct node** gb;               // pointer to GB of each grain
int numElements;                // # of elements from the periodic table for
                                // which the name and covalent radius is given.
double bohr;                    // Conversion factor for a.u. to Angstroms.
double proximity;               // This is a distance (in angstroms) that is
      // used to identify whether or not an atom is in a grain boundary or not.
      // Only atoms in a grain boundary are checked for positional overlap with
      // other atoms.

extern int UpdateFilledList(Polyhedron *ph,int nopoly,int *filledlist);
extern int NotInThisList(int *list,int size,int check);
extern int PolyhedronHasThis(Vector *r,Polyhedron *ph,int nopoly);
extern Vector latsize(Vector *a,Vector *r,int mm, Crystal *cryst);
extern int SaveAtom(Vector *rt,int j,int or,int noa,int tno,double covRad,
   int atomTypeID);
extern void bound(Vector *min,Vector *max, Vector *p);
extern int FindNeighbor(Polyhedron *ph,int nopoly,int j,int *tobefilledlist,int tobefilled);
extern int FillAt(Polyhedron *ph,int nopoly,int j,Orientation *or,Vector *r,
Crystal *cryst);
extern Polyhedron MakeBox(double *rotM);
extern double dcos(double v);
extern Vector rotation(Vector p,double *r);
extern Vector rotationT(Vector p,double *r);
extern void printvector(int i,Vector *r);
extern void printatom(Vector *a, int value, int no);
extern void printlammps(Vector *a, int atomno, int tno, int fno);
//extern int inSearchBox(Vector *v);
extern int inpolyhedron(Vector *x, Polyhedron *p);
extern int inthebox(Vector *r);
extern Vector center2(Vector *r, Vector *s);
extern Vector shift(Vector *center, Vector *normal);
extern Vector center(Vector *r, Vector *s, Vector *t);
extern Vector center4(Vector *r, Vector *s, Vector *t, Vector *u);
extern Vector Vnormal(Vector *r, Vector *s, Vector *t, Polyhedron *ph, Face *fc);
extern Vector va2b(Vector *r, Vector *s);
extern double len(Vector *a);
extern double lenA2B(Vector *r, Vector *s);
extern Vector unit(Vector *a);
extern Vector unit3(double x,double y,double z);
extern Vector vscale(Vector *a,double value);
extern Vector vmodify(Vector *a,double value,int idx);
extern Vector cross(Vector *r, Vector *s );
extern double dot(Vector *a, Vector *b);
extern Vector periodic(Vector *r);
extern Orientation MakeOrientation(int k);
extern void EulerAngles(double *rotM,double phi,double theta,double psi);
extern void EulerParameters(double *rotM,double a,double b,double c,double d);
extern double RandomWithIn(double a,double b);
extern double degree(double radian);
extern double radian(double degree);
extern void Push(struct node** headref,int no);
//extern void makeVoid(Vector *r,double vradius);
extern int countAtomsInBox(int boxno,int count);
extern int printbox(int boxno,int count);
extern int Search2Boxes(int boxA,int boxB);
extern void SearchDelete(struct node* lista,struct node* listb,double tol);
//extern int i2n(int i,int j,int k);
//extern void n2i(int n,int *i,int *j,int *k);
//extern void SearchPack(int i,int j,int k,double tol);
extern void printpercent(time_t start,double perc);
extern double  runline(Vector *r,double x,double y,double z,double dist);
extern int Vreorder(Vector *v,int i);
extern Vector spherical(Vector *v);
extern int intriangle(Vector *v);
extern Vector vpositive(Vector *a);
extern int printCellHeader(int totalAtomsPrinted, int nos, double penalty,
      unsigned int seed, Crystal *cryst, ElementList *eList);
extern int printCell(int boxno,int count);
extern int printLammpsIn(Crystal *cryst, ElementList *eList);
extern int fillElementList(ElementList *eList);
extern int getElementZandCovRad (char *elementName, ElementList *eList,
      int *atomZ, double *covRad, double *penalty);
extern int addToTypeList(int *currentAtomNum, Crystal *cryst);

// Define element names and corresponding.
