#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
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
#include "polycryst.h"
#include "global.h"

int UpdateFilledList(Polyhedron *ph,int nopoly,int *filledlist)
{
  int i,fill=0;

  // For each Polyhedron, if it is not blank then we record its index number in
  //   an array and then increment the number of filled Polyhedra.
  for (i=0;i<nopoly;i++) {
    if (!(ph->blank)) {
      *(filledlist+fill) = i;
      fill++;
    }
    ph++;
  }
  return fill;
}

int NotInThisList(int *list,int size,int check)
{
  int i,go=0;

  for (i=0;i<size;i++) {
    if (check == *list) go++;
    list++;
  }
  if (go == 0) return 1;
  else return 0;
}

int PolyhedronHasThis(Vector *r,Polyhedron *ph,int nopoly)
{
  int i,j;

  for (i=0;i<nopoly;i++) {
    if (inpolyhedron(r,ph) < 1) {
      j = i;
      break;
    } 
    ph++;
  }
  return j;
}

// The vision to have when reading this is that the original Polyhedron center
//   Vector has been rotated to point in a random direction (but it still has
//   the same magnitude). The Polyhedron is surrounded by a set of bounding
//   limits, but the distances between the min and max x,y,z limits can be any
//   real numbers. The goal here is to slightly extend the bounding limits in
//   such a way as to make the distances be an integral number of units defined
//   by the crystal cell lattice parameter. This is done by starting from the
//   Polyhedron center Vector and marching toward the min or max bounding
//   limits in units of lattice parameter lengths until we *just* surpass the
//   current limit. This will define the new bounding limit. The purpose is to
//   have an integer number of crystal lattice steps that will completely
//   cover the Polyhedron bounding limits.
// Vector a is either the min(mm==0) or the max(mm==1) extent of the Polyhedron
//   that we are attempting to fill.
// Vector r is the Vector for the Polyhedron center rotated according to the
//   rotation matrix.
Vector latsize(Vector *a, Vector *r, int mm, Crystal *cryst)
{
  Vector v = *r;

  //  March torward the current min in units of crystal lattice lengths.
  if (mm == 0) {
    while (v.x > a->x) {
      v.x -= cryst->a.x;
    }
    while (v.y > a->y) {
      v.y -= cryst->b.y;
    }
    while (v.z > a->z) {
      v.z -= cryst->c.z;
    }
  }
  //  March torward the current max in units of crystal lattice lengths.
  else if (mm == 1) {
    while (v.x < a->x) {
      v.x += cryst->a.x;
    }
    while (v.y < a->y) {
      v.y += cryst->b.y;
    }
    while (v.z < a->z) {
      v.z += cryst->c.z;
    }
  }
  return v;
}
              
int SaveAtom(Vector *rt,int j,int or,int noa,int tno,double covRad,
      int atomTypeID)
{
  // Save the location of the atom.
  atom[ano].r = *rt;

  // Identify the "nature" of the atom.  Recall that the variable tno comes
  //   from inpolyhedron and if that is 0 then the atom was in the bulk of the
  //   Polyhedron, but if it was -1 that meant it was "near" the edge where
  //   "near" is/was originally a hard coded value of 0.3.
  atom[ano].nature = tno+2; //GB = 1, internal = 2

  // Save the index number of the Polyhedron associated with this atom.
  atom[ano].ph = j;

  // Save the orientation number (which has no realy value for us).
  atom[ano].or = or;

  // Save the covalent radius.
  atom[ano].covRad = covRad;

  // Save the atom type ID number. This is the index number of the typeZ array.
  atom[ano].atomTypeID = atomTypeID;

  // Save the atom number into the searchbox and if it is close to the boundry
  //   between Polyhedra then also save the atom number into the grain boundry
  //   (gb) list for this Polyhedron.
  Push(&searchbox[j],ano);
  if (tno==-1) Push(&gb[j],ano);
  
  ano++; // Increment the overall count of atoms in the box.
  noa++; // Increment the number of atoms in this Polyhedron.

  return noa;
}

void bound(Vector *min,Vector *max, Vector *p)
{
  min->x = min(min->x,p->x);
  min->y = min(min->y,p->y);
  min->z = min(min->z,p->z);
  max->x = max(max->x,p->x);
  max->y = max(max->y,p->y);
  max->z = max(max->z,p->z);
  return;
}

int FindNeighbor(Polyhedron *ph,int nopoly,int j,int *tobefilledlist,int tobefilled)
{
  Polyhedron *phj, *phk;
  int i,k,count;
  Vector *r;

  phj=ph+j;

  //printf("Polyhedron %d of %d with %d faces %d points\n",j,nopoly,phj->nof,phj->nop);
  //printf("Neighbor:\n");
  for (i=0;i<phj->nof;i++) {
    phk=ph + phj->ngh[i];
    //printf("%d contain %d points ",phj->ngh[i],phk->nop);
    count=0;
    for (k=0;k<phk->nop;k++) {
      r=phk->p0 + phk->p[k]; // r is the point that bound this polyhedron
      //printvector(k,r);
      if (inthebox(r)) {
        count++;
      }
    }
    if (count>0) {
      //printf("In The Box %d / %d\n",count,phk->nop);
      if (phk->blank) {
        *(tobefilledlist + tobefilled) = phj->ngh[i];
        tobefilled++;
      }
    } else {
      
      //printf("\n");
    }
  }
  

  return tobefilled;
}

// Note that the Polyhedron that is passed to this subroutine is always the
//   first Polyhedron in the array list. Pointer arithmetic using the j integer
//   is used to ID the acutal targeted Polyhedron that we will play with.
// Polyhedron *ph that is targeted by this subroutine may be either inside or
//   outside of the central cell. Also, the index from the FILL POLYGON
//   comment section in fillatoms.c is "i" and here it is received as "j". This
//   may cause some confusion so take note.
// The vector *r passed in is the center point of the *targeted* Polyhedron.
// The basic idea of this subroutine is to fill a Polyhedron with atoms
//   according to the crystal cell definition. The orientation of the atomic
//   lattice is random and governed by the rotation matrix (or->rotM).
// Note that whenever we discuss rotation here it should be understood that the
//   axis of rotation is a randomly oriented line that passes through the
//   point (0,0,0).  (One helpful way to think through the operations in this
//   subroutine is to consider the case of an Identity rotation matrix so that,
//   in effect, no Vectors are rotated.)
// The algorithm that is followed is to:
//   (1) Compute the bounding limits of the rotated passed Polyhedron.
//   (2) Create a central cell copy, rotate it the same way, and then compute
//       its bounding limits.
//   (3) Redefine the Polyhedron bounding limits such that it contains only the
//       portion of the Polyhedron that is *inside* the central cell.
//   (4) Further redefine the Polyhedron bounding limits such that the lengths
//       of its sides are an integral number of crystal cell lattice lengths.
//       This will allow for easy filling and also recognition of which
//       periodic Polyhedra need to be used.
//   (5) Use three nested loops to fill in the Polyhedra and identify the need
//       for periodic Polyhedra.
int FillAt(Polyhedron *ph,int nopoly,int j,Orientation *or,Vector *r,
      Crystal *cryst)
{
  Polyhedron box,*phh,*phk,*phi; // The naming convention for these variables
        // needs to be improved. phk is always used to point to the head of the
        // Polyhedron array list. phh will point to the targeted Polyhedron.
  Vector xyz,tmp,rt,rr,pr,dummy;
  Vector *p0; // This points to the first vertex Vector of a Polyhedron which
              // is always the infinity vertex.
  Vector p1;  // This is a temp Vector that holds a rotated vertex Vector.
  int noa=0;
  double *rotM = &(or->rotM[0][0]);
  int i,k,m,extend[8],ext=0,inph;

  // Initialize the dummy vector so that min and max values can be computed.
  dummy.x = 0.0;
  dummy.y = 0.0;
  dummy.z = 0.0;

  // Make a pointer to the head of the Polyhedron array list.
  phk = ph;

  // Compute a rotated center Vector for the passed Polyhedron. Recall that the
  //   axis of rotation is a randomly oriented line passing though (0,0,0).
  rr = rotation(*r,rotM);

  // If j==0 then one of two possible conditions is true. Either (1) this is
  //   the first iteration of the loop to fill Polygons that are *inside* the
  //   the central cell (see FILL POLYGONS section in fillatoms.c), or (2) we
  //   are in the loop after the FILL POLYGONS comment so that we are going to
  //   fill the periodic Polyhedra associated with the Polyhedron that was
  //   found inside the central cell.
  // Somehow, I am hoping that the first condition can never happen because of
  //   the ordering of the Polyhedron list. (I.e. the first Polyhedron in the
  //   array of Polyhedra can never be found *inside* the central cell. However,
  //   this needs to be confirmed. Based on a cursory analysis though, I
  //   believe that this is true.) It might be better for the call to FillAt
  //   for the periodic Polyhedra to send the integer -1 instead. Then there is
  //   no chance for confusion.
  // Therefore, the most likely scenario if j==0 is that we have a periodic
  //   point. The goal of this if and for loop is to identify which polyhedron
  //   is the "owner" of the periodic point.
  if (j==0) { //Search for polyhedron[j] to fill
    for (i=0;i<nopoly;i++) {

      // Check if this periodic point is in the current i indexed Polyhedron.
      //   Note that pointer arithmetic is used below (ph++) to increment the
      //   choice of Polyhedron that is sent to the inpolyhedron subroutine.
      if (inpolyhedron(r,ph) < 1) {

        // Set index j to the index of the periodic Polyhedron that holds the
        //   center point indicated by Vector r.
        j = i;

        // Establish a second pointer to that Polyhedron and exit the loop.
        //   This is done because phh is set to ph+j below in the condition
        //   that j!=0.
        phh = ph; 
        break;
      }
      ph++; // Iterate to the next Polyhedron.
    }
  } else { // This is the case where we already know which Polyhedron to fill.
    // Establish phh as *the* pointer to the Polyhedron of interest.
    phh = ph + j;
  }


  // If it's blank then fill it, but take note that phh may be pointing to a
  //   Polyhedron that is inside or outside of the central cell.
  if (phh->blank) {

    // Record the atom # of the first atom in this Polyhedron. Recall that ano
    //   is a running total of the number of atoms.
    phh->start = ano;
    p0 = phh->p0;  // Set p0 to point to the infinity point vertex.

    // Loop over all vertices associated with this Polyhedron to set the min
    //   and max values of its x,y,z Cartesian coordinates *in the rotated
    //   frame*.
    for (i=0;i<phh->nop;i++) {

      // Rotate the vertex point according to the crystal grain orientation.
      p1 = rotation(*(p0+phh->p[i]),rotM);

      // Set the min and max Vectors if the point p1 has any value that pushes
      //   beyond the already defined limits.  Note that the phh-min and
      //   phh-max values are initialized to very large magnitude numbers so
      //   that the first point p1 will automatically be associated with both
      //   the min and the max.
      bound(&(phh->min),&(phh->max),&p1);
    }

    // Make a cubic Polyhedron that is the same size as the central cell but
    //   which is centered just off-center from the actual central cell and
    //   which is rotated just as the current Polygon and vertices of interest.
    // Some thing to be aware of is that the Vectors for the face centers and
    //   normals are all expressed in reference to the *non-rotated* box.
    box = MakeBox(rotM);

    // There are now a few different scenarios to consider:
    // (1) The target Polyhedron, whose center Vector has been rotated, is
    //   wholly inside the rotated box. All of the Polyhedron vertices are
    //   inside the rotated box. (This is the same as the unrotated Polyhedron
    //   being wholly inside the unrotated box.)
    // (2) The target Polyhedron, whose center Vector has been rotated, is
    //   partially inside the rotated box. The center Vector and some of the
    //   Polyhedron vertices are inside the rotated box, but some of the
    //   vertices are outside the rotated box.
    // (3) The target Polyhedron, whose center Vector has been rotated, is
    //   partially outside the rotated box. This is distinguished from case (2)
    //   by the center Vector being *outside* the box. As with case (2) though,
    //   some vertices are inside and some are outside of the box. There is no
    //   requirement in either case (2) or (3) that a majority of vertices be
    //   inside or outside. That allotment is completely decided by the random
    //   distribution of seed points.
    // (4) It is never possible for the target Polyhedron, whose center Vector
    //   has been rotated, to be wholly outside the rotated box. This is
    //   because the action that calls this subroutine will (first) never send
    //   a Polyhedron whos center is outside the box except for the case where
    //   a periodic Polyhedron is sent. However, so-called periodic Polyhedra
    //   are *only* sent if they were added to a specific list. To be added to
    //   that list, the periodic Polyhedron (with the center Vector outside the
    //   box) needed to have an "original" Polyhedron with its center Vector
    //   inside the box and part of its volume outside the box (i.e. case (2)).
    //   Such an event is called a case (3).

    // In all four cases we have a phh->min and phh->max that define a bounding
    //   limit for the *rotated* Polyhedron. Additionally, we have a box.min
    //   and a box.max that define the bounds of the rotated central cell.

    // Replace phh->min with the *max* of either phh->min or box.min.
    bound(&dummy,&(phh->min),&box.min);

    // Replace phh-max with the *min* of either phh->max or box.max. The
    //   rationale for these two calls to bound is to create a bounding limit
    //   that contains only the portion of the Polyhedron that is inside the
    //   the box. If the Polyhedron (including all vertices) is completely
    //   outside the box then the phh->min will necessariily become > the value
    //   of phh->max. Consequently, the while-loop set will not execute and
    //   this Polyhedron will be labeled as "blank". (It may already be blank.)
    bound(&(phh->max),&dummy,&box.max);

    // This adjustment will modify the bounding limits of the Polyhedron by
    //   creating a new set of limits that are (1) in integral units of the
    //   crystal cell lattice parameter and (2) are *just* beyond the current
    //   limits.
    // Because the bounding limits were (until here) exactly confined to the
    //   intersection of the Polyhedron and the box, this action will in effect
    //   extend that limit to include some small region that is both inside
    //   the Polyhedron but outside the box. Because we only want to record
    //   atoms that are actually inside the box, this will be the cue to
    //   make a record of the need for using a periodic copy of the current
    //   Polyhedron to fill in the region on the "other side" of the box.
    phh->min = latsize(&(phh->min),&rr,0,cryst);
    phh->max = latsize(&(phh->max),&rr,1,cryst);

    // Start in the left, lower, near corner and determine if an atom at this
    //   point would be inside the Polyhedron and the box. If it is inside
    //   both, then save it. If it is not inside the Polyhedron then don't save
    //   it and don't bother checking if it is inside the box. If it is inside
    //   the Polyhedron, but not inside the box, then this is a cue that a
    //   periodic Polyhedron needs to be used to fill in an empty space on the
    //   "other side" of the box.
    // The other atoms of the crystal cell are also considered here. If one
    //   is inside both the Polyhedron and the box, then save it. There is no
    //   need to make a note of a need for a periodic Polyhedron because the
    //   first atom (for now, this works only because the FCC first atom is at
    //   (0,0,0) of its crystal cell) will always be a sufficient trigger.
    xyz.z = phh->min.z; // Near corner
    while (xyz.z <= phh->max.z) {
      xyz.y = phh->min.y; // Lower corner
      while (xyz.y <= phh->max.y) {
        xyz.x = phh->min.x; // Left corner
        while (xyz.x <= phh->max.x) {

          // Because only the min and max Vectors of the Polyhedron struct
          //   represent the rotated Polyhedron, we need to apply a reverse
          //   rotation to the current xyz point (which is based off of the min
          //   and max limit Vectors) to determine if it is inside the
          //   Polyhedron or not (because the Faces and Normals are all still
          //   in reference to the non-rotated Polyhedron).
          // The reason that we have to ask this at all is because the bounding
          //   limits define a cube while the Polyhedron is just that: a
          //   Polyhedron which will (by definition of the limits) not fill the
          //   bounding limits.  The region outside the Polyhedron but inside
          //   the bounding limits should not be filled with atoms at this
          //   point in time. That region belongs to another Polyhedron.
          rt = rotationT(xyz,rotM);

          // Check if we are inside this Polyhedron.
          inph = inpolyhedron(&rt,phh);
          if (inph < 1) {
            // Check if we are inside the box.
            if (inpolyhedron(&rt,&box) < 1) {
              noa = SaveAtom(&rt,j,or->or,noa,inph,cryst->covRad[0],
                             cryst->atomTypeID[0]);
            }
            else {
              // We are in the Polyhedron but not in domain (periodic) so we
              //   need to make use of a periodic copy.

              // Find the Vector that is a periodic copy of the current atom.
              pr = periodic(&rt);

              // Determine the index number of the Polyhedron that holds the
              //   periodic atom.
              k = PolyhedronHasThis(&pr,phk,nopoly);

              // Make a pointer to that k indexed Polyhedron.
              phi = phk + k;

              // In the case that the k indexed Polyhedron has not yet been
              //   identified as one that is needed for "periodic" filling
              //   purposes (and which also has not yet been filled).
              if ((NotInThisList(&extend[0],ext,k)) && (phi->blank)) {
                // Not already stored, not filled, then save it.
                extend[ext] = k; phh->periodic[ext] = pr; ext++;
              }
            }
          }

          // Iterate through the other three atoms of the crystal FCC cell.
          for (m=1; m<=cryst->numAtoms-1; m++)
          {
             tmp.x = xyz.x + cryst->atomCoords[m].x;
             tmp.y = xyz.y + cryst->atomCoords[m].y;
             tmp.z = xyz.z + cryst->atomCoords[m].z;
             rt = rotationT(tmp,rotM);
             inph = inpolyhedron(&rt,phh);
             if ((inph<1)&&(inpolyhedron(&rt,&box) < 1)) {
               noa = SaveAtom(&rt,j,or->or,noa,inph,cryst->covRad[m],
                     cryst->atomTypeID[m]);
             }
          }
          xyz.x += cryst->a.x;
        }
        xyz.y += cryst->b.y;
      }
      xyz.z += cryst->c.z;
    } // END while

    // Check if atoms have been added to the list of atoms for this Polyhedron.
    //   If no atoms have been added, then mark the Polyhedron as blank. If
    //   atoms have been added then record the number for this Polyhedron,
    //   make note of the number of periodic Polyhedra
    if (noa==0) {
      phh->blank = 1;
    } else {
      phh->noa = noa;
      phh->pr_count = ext;
      phh->blank = 0;
      // printf("Ph:%5d Ano = %10d\n",j,ano);
    }
  }
  return j;
}

Polyhedron MakeBox(double *rotM)
{
  // Establish the cube size as being equal to the whole central cell dimension.
  double minx=gmx, miny=gmy, minz=gmz; 
  double maxx=px+gmx, maxy=py+gmy, maxz=pz+gmz;

  int i;
  Vector pb[8],r;
  Face fb[6];
  Polyhedron box;

  // Create boundaries for the box that retain the size, but which shift the
  //   values of the boundaries on all sides by -0.01. These boundaries are
  //   just opposite corners of the box.
  // Left (-x), Lower (-y), Near (-z)
  pb[0].x = minx-0.01;  pb[0].y = miny-0.01;  pb[0].z = minz-0.01;
  // Right (+x), Upper (+y), Far (+z)
  pb[6].x = maxx-0.01;  pb[6].y = maxy-0.01;  pb[6].z = maxz-0.01;

  // Compute the center point of the box. Note that it will be just slightly
  //   off from the center point of the given cell.
  box.com.x = (minx+maxx)/2.0;
  box.com.y = (miny+maxy)/2.0;
  box.com.z = (minz+maxz)/2.0;

  // Establish values for the other 6 corner points of the box.
  pb[1].x = pb[6].x;  pb[1].y = pb[0].y;  pb[1].z = pb[0].z; //Right,Lower,Near
  pb[2].x = pb[6].x;  pb[2].y = pb[6].y;  pb[2].z = pb[0].z; //Right,Upper,Near
  pb[3].x = pb[0].x;  pb[3].y = pb[6].y;  pb[3].z = pb[0].z; //Left, Upper,near
  pb[4].x = pb[0].x;  pb[4].y = pb[0].y;  pb[4].z = pb[6].z; //Left, Lower,Far
  pb[5].x = pb[6].x;  pb[5].y = pb[0].y;  pb[5].z = pb[6].z; //Right,Lower,Far
  pb[7].x = pb[0].x;  pb[7].y = pb[6].y;  pb[7].z = pb[6].z; //Left, Upper,Far

  // Set initial values before a search for the min and max Cartesian
  //   coordinates of the rotated box. The original way (all 0.0) is dangerous
  //   because it assumes that the box will be centered somewhat near zero so
  //   that all min values will eventually be negative and all max values will
  //   eventually be positive. If someone does not specify the central cell on
  //   the input line in this way though, then this will not work. Thus, the
  //   following is done instead.
  minx = miny = minz =  1e100; // A min value can't be any higher than this.
  maxx = maxy = maxz = -1e100; // A max value can't be any lower than this.

  // Rotate each corner point of the box and then use its new coordinates to
  //   establish the overall min and max coordinates.
  for (i=0;i<8;i++) {
    r = rotation(pb[i],rotM);   
    minx = min(minx,r.x);
    maxx = max(maxx,r.x);
    miny = min(miny,r.y);
    maxy = max(maxy,r.y);
    minz = min(minz,r.z);
    maxz = max(maxz,r.z);
  }
  box.min.x = minx;
  box.max.x = maxx;
  box.min.y = miny;
  box.max.y = maxy;
  box.min.z = minz;
  box.max.z = maxz;

  // Now, the center point of each face is determined. Note that this is in
  //   reference to the un-rotated box.
  fb[0].center = center4(&pb[1],&pb[2],&pb[5],&pb[6]);
  fb[1].center = center4(&pb[0],&pb[3],&pb[4],&pb[7]);
  fb[2].center = center4(&pb[2],&pb[3],&pb[6],&pb[7]);
  fb[3].center = center4(&pb[0],&pb[1],&pb[4],&pb[5]);
  fb[4].center = center4(&pb[4],&pb[5],&pb[6],&pb[7]);
  fb[5].center = center4(&pb[0],&pb[1],&pb[2],&pb[3]);
  fb[0].normal = Vnormal(&pb[1],&pb[2],&pb[5],&box,&fb[0]);
  fb[1].normal = Vnormal(&pb[0],&pb[3],&pb[4],&box,&fb[1]);
  fb[2].normal = Vnormal(&pb[2],&pb[3],&pb[6],&box,&fb[2]);
  fb[3].normal = Vnormal(&pb[0],&pb[1],&pb[4],&box,&fb[3]);
  fb[4].normal = Vnormal(&pb[4],&pb[5],&pb[6],&box,&fb[4]);
  fb[5].normal = Vnormal(&pb[0],&pb[1],&pb[2],&box,&fb[5]);

  //  Record the center and normal vectors for each face of the box.
  for (i=0;i<6;i++) {
    box.c[i] = fb[i].center;
    box.n[i] = fb[i].normal;
  }
  box.nof = 6; // Set the number of faces to 6

  return box;
}

double dcos(double v)
{
  double pi=4.0*atan2(1.0,1.0);
  double rad=(v/180.)*pi;

  return cos(rad);
}

Vector rotation(Vector p,double *r)
{
/*
a = *(r+0); b = *(r+1); c = *(r+2)
d = *(r+3); e = *(r+4); f = *(r+5)
g = *(r+6); h = *(r+7); i = *(r+8)

   --     --   -- --   --              --   --   --
   | a b c |   | x |   | (ax + by + cz) |   | u.x |
   | d e f | x | y | = | (dx + ey + fz) | = | u.y |
   | g h i |   | z |   | (gx + hy + iz) |   | u.z |
   --     --   -- --   --              --   --   --
*/

  Vector u;

  u.x = *(r+0)*p.x + *(r+1)*p.y + *(r+2)*p.z;
  u.y = *(r+3)*p.x + *(r+4)*p.y + *(r+5)*p.z;
  u.z = *(r+6)*p.x + *(r+7)*p.y + *(r+8)*p.z;
  
  return u;
}

Vector rotationT(Vector p,double *r)
{
/*
a = *(r+0); b = *(r+1); c = *(r+2)
d = *(r+3); e = *(r+4); f = *(r+5)
g = *(r+6); h = *(r+7); i = *(r+8)

Use the Transpose of the rotation matrix to reverse the rotation operation.
   --     --   -- --   --              --   --   --
   | a d g |   | x |   | (ax + dy + gz) |   | u.x |
   | b e h | x | y | = | (bx + ey + hz) | = | u.y |
   | c f i |   | z |   | (cx + fy + iz) |   | u.z |
   --     --   -- --   --              --   --   --
*/

  Vector u;

  u.x = *(r+0)*p.x + *(r+3)*p.y + *(r+6)*p.z;
  u.y = *(r+1)*p.x + *(r+4)*p.y + *(r+7)*p.z;
  u.z = *(r+2)*p.x + *(r+5)*p.y + *(r+8)*p.z;
  
  return u;
}

/* These subroutines are not used anymore, if they ever were.
void printvector(int i,Vector *r)
{
  printf("%3d %14.7f %14.7f %14.7f\n",i,r->x,r->y,r->z);
  return;
}

void printatom(Vector *a, int value, int no)
{
  FILE *fj;
  char fname[20];

  if (no < 10) sprintf(fname,"polygon-000%d",no);
  else if (no < 100) sprintf(fname,"polygon-00%d",no);
  else if (no < 1000) sprintf(fname,"polygon-0%d",no);
  else sprintf(fname,"polygon-%d",no);
  fj = fopen(fname,"a");

  fprintf(fj,"%14.10f %14.10f %14.10f %d\n",a->x,a->y,a->z,value);
  fclose(fj);
  return;
}
*/

/*int inSearchBox(Vector *v)
{
  int i,j,k,r=0;
  double minx,miny,minz,maxx,maxy,maxz;

  for (i=0;i<nx;i++) {
    for (j=0;j<ny;j++) {
      for (k=0;k<nz;k++) {
        minx = gmx+i*bx;
        miny = gmy+j*by;
        minz = gmz+k*bz;
        maxx = minx + bx;
        maxy = miny + by;
        maxz = minz + bz;
        if ((v->x >= minx)&&(v->x <= maxx)&&
            (v->y >= miny)&&(v->y <= maxy)&&
            (v->z >= minz)&&(v->z <= maxz)) {
          r = i*ny*nz+j*nz+k;
          i = nx;
          j = ny;
          k = nz;
        }
      }
    }
  }
  return r;
}*/

int inpolyhedron(Vector *x, Polyhedron *p) //FIXME
{
  Vector a2x;
  int i,sum=0;
  double val;
  double v=-proximity;

  // Loop over each face of this Polyhedron.
  for (i=0;i<p->nof;i++) {

    // Compute a vector from the center of the current face to the point
    //   indicated by the passed vector x. For example, this vector could be
    //   the center of a Polyhedron.
    a2x = va2b(&p->c[i],x);

    // Compute the dot product between the difference vector computed above and
    //   the normal vector for the current face (which points outward from the
    //   cell center).
    val = dot(&a2x,&p->n[i]);

    // Because d = |a|*|b|*cos(theta) we have that the sign of d is determined
    //   by whether or not the position pointed to by the x vector is beyond(+)
    //   the face {with theta in [0,90)} or behind(-) the face {with theta in
    //   (90,180]}.

    // In the event that the given point is beyond the face, we record the
    //   face number for which this occured.
    if (val > 0.0) sum=i+1;

    // Otherwise val < 0.0 and we want to know how close the given vector is
    //   to pointing right on the face. Thus, if val is e.g. -0.2 or -0.1
    //   then it might be "too close" to a face.
    else v=max(v,val);
  }
  // In the condition that we found a face for which the given vector is
  //   "really really" close, then we set sum to -1.  The point is still inside
  //   the Polyhedron, but it is close to a face and thus may require further
  //   consideration.
  if ((v > -proximity)&&(sum==0)) {sum=-1;};

  // We can return three basic types of values here.  If we return 0 then the
  //   point given by vector x is inside the current Polyhedron and it is not
  //   near any face.  If we return -1, then it is still inside the current
  //   Polyhedron, but it is also potentially "too close" to a face.  (It may
  //   be good to allow some control over the measure of "too close" by the
  //   program user.) Finally, we may return some positive value that indicates
  //   the sequential number (the array index number + 1) of the last face that
  //   the point was "beyond". This defines the "nature" of the atom.
  return sum;
}

int inthebox(Vector *r)
{
  if ((r->x > gmx)&&(r->x < gmx+px))
    if ((r->y > gmy)&&(r->y < gmy+py))
      if ((r->z > gmz)&&(r->z < gmz+pz)) return 1;
      else return 0;
    else return 0;
  else return 0;
}

Vector center2(Vector *r, Vector *s)
{
  Vector c;
  
  c.x = (r->x + s->x)/2.0;
  c.y = (r->y + s->y)/2.0;
  c.z = (r->z + s->z)/2.0;
  return c;
}

Vector shift(Vector *center, Vector *normal) //FIXME
{
  Vector c;
  
  c = unit(normal);
  c = vscale(&c,1.2781*0.5);
  c = va2b(&c,center);
  return c;
} //FIXME

Vector center(Vector *r, Vector *s, Vector *t)
{
  Vector c;
  
  c.x = (r->x + s->x + t->x)/3.0;
  c.y = (r->y + s->y + t->y)/3.0;
  c.z = (r->z + s->z + t->z)/3.0;
  return c;
}

Vector center4(Vector *r, Vector *s, Vector *t, Vector *u)
{
  Vector c;
  
  c.x = (r->x + s->x + t->x + u->x)/4.0;
  c.y = (r->y + s->y + t->y + u->y)/4.0;
  c.z = (r->z + s->z + t->z + u->z)/4.0;
  return c;
}

Vector Vnormal(Vector *r, Vector *s, Vector *t, Polyhedron *ph, Face *fc)
{
  Vector u,v,n,o;
  double d;
  
  u = va2b(r,s);
  v = va2b(s,t);
  n = cross(&u,&v);
  o = va2b(&ph->com,&fc->center);
  d = dot(&o,&n);
  if (d<0.0) n = cross(&v,&u);
  return n;
}

// I have renamed the variables in this subroutine to more clearly indicate
//   that the resultant vector points from a to b.
Vector va2b(Vector *a, Vector *b)
{
  Vector u;

  u.x = b->x - a->x;
  u.y = b->y - a->y;
  u.z = b->z - a->z;
  return u;
}
/*Vector va2b(Vector *r, Vector *s)
{
  Vector u;

  u.x = s->x - r->x;
  u.y = s->y - r->y;
  u.z = s->z - r->z;
  return u;
}*/

double len(Vector *a)
{
  double l;

  l = a->x * a->x + a->y * a->y + a->z * a->z;
  l = sqrt(l);
  return l;
}

double lenA2B(Vector *r, Vector *s)
{
  Vector u;
  double l;

  u.x = s->x - r->x;
  u.y = s->y - r->y;
  u.z = s->z - r->z;
  l = u.x * u.x + u.y * u.y + u.z * u.z;
  l = sqrt(l);

  return l;
}

Vector unit(Vector *a)
{
  Vector u;
  double l=len(a);
  
  u.x = a->x / l;
  u.y = a->y / l;
  u.z = a->z / l;
  return u;
}

Vector unit3(double x,double y,double z)
{
  Vector u;
  double l=sqrt(x*x + y*y + z*z);
  
  u.x = x / l;
  u.y = y / l;
  u.z = z / l;
  return u;
}

Vector vscale(Vector *a,double value)
{
  Vector s;

  s.x = a->x * value;
  s.y = a->y * value;
  s.z = a->z * value;
  return s;
}

Vector vmodify(Vector *a,double value,int idx)
{
  Vector u;
  
  if ((idx>=0)&&(idx<=2)) {
    if (idx==0) {
      u.x = a->x + value;
      u.y = a->y;
      u.z = a->z;
    } else if (idx==1) {
      u.x = a->x;
      u.y = a->y + value;
      u.z = a->z;
    } else if (idx==2) {
      u.x = a->x;
      u.y = a->y;
      u.z = a->z + value;
    }
  }

  return u;
}

Vector cross(Vector *r, Vector *s )
{
  Vector u;
  
  u.x = r->y * s->z - r->z * s->y;
  u.y = r->z * s->x - r->x * s->z;
  u.z = r->x * s->y - r->y * s->x;
  return u;
}

double dot(Vector *a, Vector *b)
{
  double d;

  d = a->x * b->x + a->y * b->y + a->z * b->z;
  return d;
}

Vector periodic(Vector *r)
{
  Vector u = *r;

  if (r->x < gmx) u.x = r->x + px;
  else if (r->x > (gmx+px)) u.x = r->x - px;
 
  if (r->y < gmy) u.y = r->y + py;
  else if (r->y > (gmy+py)) u.y = r->y - py;
  
  if (r->z < gmz) u.z = r->z + pz;
  else if (r->z > (gmz+pz)) u.z = r->z - pz;

  return u;
}

int Triangle2Color(Vector *r)
{
  int color;
  Vector d100,d110,d111;
  int l1,l2,l3;
  
  d100 = unit3(1.0,0.0,0.0);
  d110 = unit3(1.0,1.0,0.0);
  d111 = unit3(1.0,1.0,1.0);
  l1 = (int)(15.0*(1.0-lenA2B(r,&d100)));
  l2 = (int)(15.0*(1.0-lenA2B(r,&d110)));
  l3 = (int)(15.0*(1.0-lenA2B(r,&d111)));
  color=l1*16*16 + l2*16 + l3;
  
  return(color);
}

// Two different methods are available here for making a rotation matrix.
// (1) If the incoming variable k is 1,2,3 then a rotation matrix for a very
//    specific angle is created. The method is based on Euler Angles.
// (2) If the incoming variable k is 0, then a random rotation matrix is
//    created.
Orientation MakeOrientation(int k)
{
  Orientation or;
  double theta; // Used only for the EulerAngle method.
  double phiOver2,a,b,c,d; // Used only for the EulerParameter method.
  FILE *fi;
  Vector Z,u,p;
  int i;

  if (k==1) { // <100> || Z-axis
    //or.or = k;
    EulerAngles(&or.rotM[0][0],radian(90.0),radian(90.0),radian(0.0));
  } else if (k==2) { // <110> || Z-axis
    //or.or = k;
    EulerAngles(&or.rotM[0][0],radian(135.0),radian(90.0),radian(270.0));
  } else if (k==3) { // <111> || Z-axis
    //or.or = k;
    theta = radian(54.735610317245);
    EulerAngles(&or.rotM[0][0],radian(135.0),theta,radian(270.0));
  } else {
    a = RandomWithIn(-1.0,1.0);
    phiOver2  = acos(a);

    // Create a random unit vector.
    u.x = RandomWithIn(-1.0,1.0);
    u.y = RandomWithIn(-1.0,1.0);
    u.z = RandomWithIn(-1.0,1.0);
//    b = RandomWithIn(-1.0,1.0); u.x = b;
//    c = RandomWithIn(-1.0,1.0); u.y = c;
//    d = RandomWithIn(-1.0,1.0); u.z = d;
    u = unit(&u);

    // Compute the other Euler Parameters.
    b = u.x * sin(phiOver2);
    c = u.y * sin(phiOver2);
    d = u.z * sin(phiOver2);
    EulerParameters(&or.rotM[0][0],a,b,c,d); //Produce a rotation matrix.
    //or.or = (int)((a+1.)*5000 + (b+1.)*500 + (c+1.)*50 + (d+1.)*5);
  }

  Z.x = Z.y = 0.0; Z.z = 1.0;
  Z = rotationT(Z,&or.rotM[0][0]);  // <001> in crystal coordinate
  
  p = unit(&Z); 
  p = vpositive(&p);
  i=0; while ( ( !(intriangle(&p)) ) && (i<5) ) i = Vreorder(&p,i);
  or.or = Triangle2Color(&p);

  fi = fopen("euler_parameters","a");
  fprintf(fi,"%5d%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f\n",
    or.or,a,b,c,d,Z.x,Z.y,Z.z,p.x,p.y,p.z);
  fclose(fi);

  return or;
}

// See (for example) Equation 48.4 in Corben and Stehle.
void EulerAngles(double *rotM,double phi,double theta,double psi)
{ // phi = [0,2*pi] : theta = [0,pi] : psi = [0,2*pi]
  *(rotM+0) = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi);
  *(rotM+1) = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi);
  *(rotM+2) = sin(psi)*sin(theta);
  *(rotM+3) = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi);
  *(rotM+4) = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi);
  *(rotM+5) = cos(psi)*sin(theta);
  *(rotM+6) = sin(theta)*sin(phi);
  *(rotM+7) = -sin(theta)*cos(phi);
  *(rotM+8) = cos(theta);
  return;
}

// See (for example) https://en.wikipedia.org/wiki/Euler-Rodrigues_formula
void EulerParameters(double *rotM,double a,double b,double c,double d)
{ // a = e0 = cos(phi/2) : bi+cj+dk = e1i+e2j+e3k => normalize * sin(phi/2)
  *(rotM+0) = a*a + b*b - c*c - d*d;
  *(rotM+1) = 2*(b*c + a*d);
  *(rotM+2) = 2*(b*d - a*c);
  *(rotM+3) = 2*(b*c - a*d);
  *(rotM+4) = a*a - b*b + c*c - d*d;
  *(rotM+5) = 2*(c*d + a*b);
  *(rotM+6) = 2*(b*d + a*c);
  *(rotM+7) = 2*(c*d - a*b);
  *(rotM+8) = a*a - b*b - c*c + d*d;
  return;
}

double RandomWithIn(double a,double b)
{
  double r;
  
  r = (double)rand() / ((double)(RAND_MAX) + (double)(1));
  return r*( max(a,b) - min(a,b) ) + min(a,b);  
}

double degree(double radian)
{
  double pi=4.0*atan2(1.0,1.0);
  return (radian/pi)*180.0;
}

double radian(double degree)
{
  double pi=4.0*atan2(1.0,1.0);
  return (degree/180.)*pi;
}

void Push(struct node** headref,int no)
{
  struct node* newnode = malloc(sizeof(struct node));
  newnode->no = no;
  newnode->next = *headref;
  *headref = newnode;
}

/*void makeVoid(Vector *r,double vradius)
{
  int boxno = inSearchBox(r);
  int a,b,c,i,j,k;
  struct node* ca;
  struct node* lista;
  double dist;

  n2i(boxno,&a,&b,&c);
  for (i=a-1;i<a+2;i++)
    for (j=b-1;j<b+2;j++)
      for (k=c-1;k<c+2;k++) {
        lista = searchbox[i2n(i,j,k)];
        for (ca=lista;ca!=NULL;ca=ca->next) {
          dist = lenA2B(&atom[ca->no].r,r);
          if (dist < vradius) {
            atom[ca->no].ph = -2;
          }
        }
      }
  return;
}*/

int countAtomsInBox(int boxno,int count)
{
  struct node* ca;
//  Vector r;  // Set but ununsed.

  for (ca=searchbox[boxno];ca!=NULL;ca=ca->next) {
//    r = atom[ca->no].r;  // Set but unused.
    if (atom[ca->no].ph > 0)
      {count++;}
  }
  return(count);
}

// This subroutine is not used any more. The printing of the .nature value is
//   not used in the new subroutine printCell. Instead, the atom type is
//   printed.
int printbox(int boxno,int count)
{
  struct node* ca;
  Vector r;
  FILE *fi,*fj;
  char fname[20],lname[20];

  sprintf(fname,"polygon-%04d",boxno);
  sprintf(lname,"Lpolygon-%04d",boxno);
  fi = fopen(fname,"w");
  fj = fopen(lname,"w");

  for (ca=searchbox[boxno];ca!=NULL;ca=ca->next) {
    r = atom[ca->no].r;
    fprintf(fi,"%14.10f %14.10f %14.10f %d\n",r.x,r.y,r.z,atom[ca->no].or);
    if (atom[ca->no].ph > 0) {
      fprintf(fj,"%7d%7d%14.8f%14.8f%14.8f\n",count,atom[ca->no].nature,
            r.x,r.y,r.z);
      count++;
    }
  }
  fclose(fi);
  fclose(fj);

  return(count);
}

//Added by Rulis June 2015.
int printCellHeader(int totalAtomsPrinted, int nos, double penalty,
      unsigned int seed, Crystal *cryst, ElementList *eList)
{
  int i;
  FILE *fj;
  char fname[20];

  sprintf(fname,"lammps-polycryst.dat");
  fj = fopen(fname,"w");

  fprintf(fj,"#Generated by the fillatoms program\n\n");
  fprintf(fj,"#Input file to regenerate this configuration:\n");
  fprintf(fj,"#%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
        gmx,gmy,gmz,px+gmx,py+gmy,pz+gmz);
  fprintf(fj,"#   %d %8.3f %u %d %8.3f\n",nos,penalty,seed,1,proximity);
  fprintf(fj,"#%15.10f   %15.10f   %15.10f\n",
        cryst->a.x, cryst->a.y, cryst->a.z);
  fprintf(fj,"#%15.10f   %15.10f   %15.10f\n",
        cryst->b.x, cryst->b.y, cryst->b.z);
  fprintf(fj,"#%15.10f   %15.10f   %15.10f\n",
        cryst->c.x, cryst->c.y, cryst->c.z);
  fprintf(fj,"#%d\n",cryst->numAtoms);
  for (i=0;i<cryst->numAtoms;i++)
     {fprintf(fj,"# %d %d %15.10f %15.10f %15.10f %s\n",i,cryst->atomTypeID[i],
            cryst->atomCoords[i].x,cryst->atomCoords[i].y,
            cryst->atomCoords[i].z,eList->elementName[cryst->atomZ[i]]);}
  fprintf(fj,"\n");
  fprintf(fj,"%d atoms\n",totalAtomsPrinted);
  fprintf(fj,"%d atom types\n",cryst->numTypes);
  fprintf(fj,"%14.10f %14.10f xlo xhi\n",gmx,px+gmx);
  fprintf(fj,"%14.10f %14.10f ylo yhi\n",gmy,py+gmy);
  fprintf(fj,"%14.10f %14.10f zlo zhi\n\n",gmz,pz+gmz);
  fprintf(fj,"Atoms # atomic\n\n");

  fclose(fj);

  return(1);
}


int printCell(int boxno,int count)
{
  struct node* ca;
  Vector r;
  FILE *fj;
  char fname[20];

  sprintf(fname,"lammps-polycryst.dat");
  fj = fopen(fname,"a");

  for (ca=searchbox[boxno];ca!=NULL;ca=ca->next) {
    r = atom[ca->no].r;
    if (atom[ca->no].ph > 0) {
      fprintf(fj,"%7d%7d%14.8f%14.8f%14.8f\n",count,atom[ca->no].atomTypeID,
            r.x,r.y,r.z);
      count++;
    }
  }
  fclose(fj);

  return(count);
}

int printLammpsIn(Crystal *cryst, ElementList *eList)
{
  int i;
  FILE *fj;
  char fname[20];

  sprintf(fname,"lammps-polycryst.in");
  fj = fopen(fname,"w");

  fprintf (fj,"#A lammps command file for minimizing and obtaining cohesive\n");
  fprintf (fj,"#energy.  Modeled after the work of Mark Tschopp, 2010.\n");
  fprintf (fj,"\n");
  fprintf (fj,"#Initialize simulation\n");
  fprintf (fj,"\n");
  fprintf (fj,"clear\n");
  fprintf (fj,"units metal\n");
  fprintf (fj,"dimension 3\n");
  fprintf (fj,"boundary p p p\n");
  fprintf (fj,"atom_style atomic\n");
  fprintf (fj,"#atom_modify map array sort 10000 0.0 #May improve perform.\n");
  fprintf (fj,"\n");
  fprintf (fj,"#Read in the system configuration\n");
  fprintf (fj,"\n");
  fprintf (fj,"read_data polycryst.dat\n");
  fprintf (fj,"\n");
  fprintf (fj,"#Define the interatomic potential\n");
  fprintf (fj,"\n");
  fprintf (fj,"pair_style eam/alloy\n");
  fprintf (fj,"pair_coeff * * XXX.lammps.eam");
  for (i=0;i<=cryst->numTypes-1;i++)
     {fprintf (fj," %s",eList->elementName[cryst->typeZ[i]]);}
  fprintf (fj,"\n");
  fprintf (fj,"neighbor 2.0 bin\n");
  fprintf (fj,"neigh_modify delay 10 check yes\n");
  fprintf (fj,"\n");
  fprintf (fj,"#Define compute settings\n");
  fprintf (fj,"\n");
  fprintf (fj,"compute eng all pe/atom\n");
  fprintf (fj,"compute eatoms all reduce sum c_eng\n");
  fprintf (fj,"dump coords all atom 100 polycryst.dump\n");
  fprintf (fj,"\n");
  fprintf (fj,"#Execution code\n");
  fprintf (fj,"\n");
  fprintf (fj,"reset_timestep 0\n");
  fprintf (fj,"fix 1 all box/relax iso 0.0 vmax 0.001\n");
  fprintf (fj,"thermo 10\n");
  fprintf (fj,"thermo_style custom step pe lx ly lz press ");
  fprintf (fj,      "pxx pyy pzz c_eatoms\n");
  fprintf (fj,"min_style cg\n");
  fprintf (fj,"minimize 1e-25 1e-25 5000 10000\n");
  fprintf (fj,"\n");
  fprintf (fj,"variable natoms equal \"count(all)\"\n");
  fprintf (fj,"variable teng equal \"c_eatoms\"\n");
  fprintf (fj,"variable length equal \"lx\"\n");
  fprintf (fj,"variable ecoh equal \"v_teng/v_natoms\"\n");
  fprintf (fj,"\n");
  fprintf (fj,"print \"Total energy (eV) = ${teng};\"\n");
  fprintf (fj,"print \"Number of atoms = ${natoms};\"\n");
  fprintf (fj,"print \"Lattice constant (Angstroms) = ${length};\"\n");
  fprintf (fj,"print \"Cohesive energy (eV) = ${ecoh};\"\n");

  fclose(fj);

  return(1);
}

int Search2Boxes(int boxA,int boxB)
{
  struct node* ca;
  struct node* cb;
  double dist;
  int count=0; //FIXME

  // We know that the Polyhedron indexed at boxA was filled and that it must
  //   therefore have some atoms in a grain boundary (gb) region. The first
  //   for-loop will consider each atom in the gb list in turn.
  for (ca=gb[boxA];ca!=NULL;ca=ca->next) {

    // Now, we do the same for the Polyhedron indexed at boxB (which is a
    //   Polyhedron on the other side of a Face from Polyhedron boxA).
    for (cb=gb[boxB];cb!=NULL;cb=cb->next) {

      // Compute the distance between the two identified gb atoms.
      dist = lenA2B(&atom[ca->no].r,&atom[cb->no].r);
      //val=min(val,dist);

      // If the sum of the covalent radii (which were stored with a
      //   multiplicative penalty) of the two atoms is greater than the
      //   distance between the atoms then the atoms are considered to be "too
      //   close" and atom "b" is labeled to be ignored (unless atom "a" is
      //   already labeled to be ignored).
      if ((dist < (atom[ca->no].covRad+atom[cb->no].covRad)) &&
          (atom[ca->no].ph != 0))
      {
        count++;
        atom[cb->no].ph = 0;
        printf("%lf  %lf dist=%lf\n",atom[ca->no].covRad,atom[cb->no].covRad,
              dist);
      }
    }
  }

  return(count);
}

void SearchDelete(struct node* lista,struct node* listb,double tol)
{
  struct node* ca;
  struct node* cb;
  double dist;

  for (ca=lista;ca!=NULL;ca=ca->next) {
    for (cb=listb;cb!=NULL;cb=cb->next) {
      if ((atom[ca->no].ph >= 0)&&(atom[ca->no].ph < atom[cb->no].ph)) {
        dist = lenA2B(&atom[ca->no].r,&atom[cb->no].r);
        if (( dist > 1e-5 )&&( dist < tol )) {
          atom[cb->no].ph = -1;
        }
      }
    }
  }
  return;
}

/*int i2n(int i,int j,int k)
{
  if (i==-1) i = nx-1;
  if (j==-1) j = ny-1;
  if (k==-1) k = nz-1;
  if (i==nx) i = 0;
  if (j==ny) j = 0;
  if (k==nz) k = 0;
  return i*ny*nz+j*nz+k;
}

void n2i(int n,int *i,int *j,int *k)
{
  *i = n / (ny*nz);
  n -= *i * ny * nz;
  *j = n / nz;
  *k = n - *j * nz;
  return;
}

void SearchPack(int i,int j,int k,double tol)
{
  SearchDelete(searchbox[i2n(i,j,k)],searchbox[i2n(i,j,k)],tol);
  SearchDelete(searchbox[i2n(i,j,k)],searchbox[i2n(i,j,k+1)],tol);
  SearchDelete(searchbox[i2n(i,j,k)],searchbox[i2n(i,j+1,k+1)],tol);
  SearchDelete(searchbox[i2n(i,j,k)],searchbox[i2n(i,j+1,k)],tol);
  SearchDelete(searchbox[i2n(i,j,k)],searchbox[i2n(i+1,j,k)],tol);
  SearchDelete(searchbox[i2n(i,j,k)],searchbox[i2n(i+1,j,k+1)],tol);
  SearchDelete(searchbox[i2n(i,j,k)],searchbox[i2n(i+1,j+1,k+1)],tol);
  SearchDelete(searchbox[i2n(i,j,k)],searchbox[i2n(i+1,j+1,k)],tol);
  return;
}*/

void printpercent(time_t start,double perc)
{
  time_t end;
  double runtime;
  int hours,minutes,seconds;

  time(&end);
  runtime = difftime(end,start);
  hours   = runtime / (60*60);
  minutes = (runtime - hours*60*60)/60;
  seconds = runtime - hours*60*60 - minutes*60;

  printf("%6.2f%% done in %2d hours, %2d minutes, %2d seconds.\n",perc,hours,minutes,seconds);

  return;
}

double  runline(Vector *r,double x,double y,double z,double dist) //FIXME
{
  Vector p,d;
  double t;

  p.x = x - r->x; d.x = x; 
  p.y = y - r->y; d.y = y; 
  p.z = z - r->z; d.z = z;
  t = lenA2B(r,&d);
  p = unit(&p);

  r->x = r->x + dist*p.x;
  r->y = r->y + dist*p.y;
  r->z = r->z + dist*p.z;

  return t;
} //FIXME

int Vreorder(Vector *v,int i)
{
  double a,b,c;
  
  a = v->x; b = v->y; c = v->z;
  if (i==0) {
    v->x = b; v->y = c; v->z = a; i = 1;
  } else if (i==1) {
    v->x = b; v->y = c; v->z = a; i = 2;
  } else if (i==2) {
    v->x = b; v->y = a; v->z = c; i = 3;
  } else if (i==3) {
    v->x = b; v->y = c; v->z = a; i = 4;
  } else if (i==4) {
    v->x = b; v->y = c; v->z = a; i = 5;
  } else {
    printf("Something might be wrong!!\n");
  }
  return i;
}

Vector spherical(Vector *v)
{
  double rho,phi,theta;
  Vector r;

  r.x = rho   = len(v);             // radial distance
  r.y = phi   = acos(v->z / rho);   // zenith angle from pos-z (north pole)
  r.z = theta = atan(v->y / v->x);  // azimuth angle from pos-x
  return r;
}

int intriangle(Vector *v)           
{ // that 100, 110, 111 triangle
  Vector r = spherical(v),n,o,a2x;
  double pi=4.0*atan2(1.0,1.0);
  
  o.x = o.y = o.z = 0.0;            // origin
  n.x = 0.0; n.y = -1.0; n.z = 1.0; // normal = cross(100,111)
  a2x = va2b(&o,v);
  if ((r.z < pi/4)&&( dot(&a2x,&n)<=0 )) return 1;
  else return 0;
}

Vector vpositive(Vector *a)
{
  Vector s = *a;

  if (s.x < 0.0) s.x *= -1.0;
  if (s.y < 0.0) s.y *= -1.0;
  if (s.z < 0.0) s.z *= -1.0;
  return s;
}

int fillElementList(ElementList *eList)
{
   int i;
   FILE *eData;
   char line[80];

   eData = fopen("elements.dat","r");

   // Read the number of elements.
   fgets(line,80,eData);
   sscanf(line,"%d",&numElements);


   // Allocate space to hold the element name list and associated covalent
   //   radii. The trick is that we want the index number of the array to
   //   match the atomic Z number.  Thus we make index 0 empty and allocate
   //   numElements+1 amount of space.
   eList->elementName = (char**)malloc((numElements+1)*sizeof(char*));
   for (i=0;i<numElements+1;i++)
      {eList->elementName[i] = (char *)malloc(3*sizeof(char));}
   eList->covRad      = (double *)malloc((numElements+1)*sizeof(double));

   // There are no elements with Z=0.
   eList->elementName[0] = "empty";
   eList->covRad[0]      = 0.0;

   // Read the element names so that the index number matches the Z number.
   for (i=1;i<=numElements;i++)
   {
      fgets(line,80,eData);
      sscanf(line,"%s",eList->elementName[i]);
   }

   // Read the element covalent radii to correspond with the element names.
   for (i=1;i<=numElements;i++)
   {
      fgets(line,80,eData);
      sscanf(line,"%lf",&eList->covRad[i]);
   }

   fclose(eData);

   return 1;
}

// We are passed in a string holding the element name for this current atom and
//   we will return the atomic Z number and covalent radius for this atom.
int getElementZandCovRad (char *elementName, ElementList *eList, int *atomZ,
      double *covRad, double *penalty)
{
   int i;

   // First we iterate through all the names in the list of element names from
   //   the periodic table looking for a match.  (This should probably deal
   //   with case appropriately at some point in the future.)
   for (i = 1; i<=numElements; i++)
   {
      if (strcmp(elementName,eList->elementName[i])==0)
      {
         *atomZ = i;
         *covRad = eList->covRad[i] * *penalty;
         break;
      }
   }

   return 1;
}

int addToTypeList (int *currentAtomNum, Crystal *cryst)
{
   // If the type number of the currently read in atom is greater than the
   //   number of types found so far, then store the Z number of this type in
   //   the typeZ array in the numTypes slot and then increment numTypes. The
   //   assumption built in here is that the atomic types are in increasing
   //   order in the input file.  It is okay to have types such as 1,2,3,2 but
   //   it is not okay to have 1,3,2 because the 2 will not be recognized as a
   //   new type.
   if (cryst->atomTypeID[*currentAtomNum] > cryst->numTypes)
   {
      cryst->typeZ[cryst->numTypes] = cryst->atomZ[*currentAtomNum];
      cryst->numTypes++;
   }

   return 1;
}
