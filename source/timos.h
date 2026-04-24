//-------------------------------------------------
// Copyright: 
// Biomedical Imaging Division 
// School of Biomedical Engineering and Sciences
// Virginia Tech
// 2009
//----------------------------------------------------

// Main head file for TIMOS.

#ifndef TIMOS_H
#define TIMOS_H

// Structure for node or vector
struct TNode{
  double X;
  double Y;
  double Z;
};

// Structure for triangle
// 1> Each triangle has three nodes. 
//    The information of the three nodes and the area of the triangle 
//    is only needed in pre-processor and post-processor.

struct TTriNode{
  // N is the node index list of a triangle
  // We only need the first three spaces (N[0] to N[2]).
  // The reason we declare N[4] is that this will align the data to 8 bytes.
  int N[4];

  // The area of the triangle
  double Area;
};

// Structure for triangle
// 2> A triangle could connect at most two tetrahedra (Elements).
//    If a triangle only connects with one tetrahedron, 
//    this triangle is a boundary (surface) trianlgle.
// 3> TTrianlge stores the normal vector of the triangle.
//    It also stores the plane function of that triangle.

struct TTriangle{
     // The index of this triangle in BTrig list.
  int BoundaryIdx;

  // If this triangle is on the boundary, it is 1, otherwise, 2.
  int Num_Elem; 

  // Tetrahedron elements contain this triangle
  int ElemIdx[2]; // the initial value should be -1 (no such Elem). 

  // The normal vector of the triangle.
  // X*X + Y*Y + Z*Z == 1
  double X;
  double Y;
  double Z;

  // The triangle is on a plane of   a x + b y + c z + d = 0
  // Since we already know the normal vector of this triangle,
  // we have: a = X,  b = Y, and c = Z;
  double d;
};

// This is a temp structure.
struct TSimpleTri{
  int N[3];
  int NumConElem;
  int ElemIdx[2];
};

// Structures for elements (tetrahedron)
// 1> Each element contains four nodes and four triangles.
//    For reason of efficiency, we move the node and triangle information into TElemNode

struct TElemNode{
  // N[0] to N[3] are the node index list of the tetrahedron
  // The indexs of nodes are from small to large
  int N[4];
  // T[0] to T[3] are the triangle index list of the tetrahedron
  /* The order in the triangle index is not true in this new version anymore and there is no need to keep this order. 
  // T[0] should be nodes 0, 1, 2
  // T[1] should be nodes 0, 1, 3
  // T[2] should be nodes 0, 2, 3
  // T[3] should be nodes 1, 2, 3
  */
  int T[4]; 

  // Volumn of the tetrahedron
  double Vol;

};

struct TElem{ 
  // MedIdx
  int MedIdx;
  int Blank;

  int AdjElemIdx[4]; 
  // If AdjElemIdx[i] > 0, then it is the AdjElemIdx
  // If AdjElemIdx[i] < 0, then it's negation is the boundary triangle index

  double TriNorm [16];
  // Norm vectors of the four triangles. 
  // All the normal vectors should point in the tetrahedron 
  // The normal vector of the triangle which is point into the internal of the tetrahedron
  // X*X + Y*Y + Z*Z == 1
  // The triangle is on a plane of   a x + b y + c z + d = 0, with:
  //  double a = X;
  //  double b = Y;
  //  double c = Z;

  //  double X;
  //  double Y;
  //  double Z;
  //  double d;

  //  double Area    [4];
};

// Temp structure for elements
struct TSElem{
  int Idx;
  int NCElem;
  int ElemList[4];
};

// Structure for optical parameters
struct TMedOptic{
  double mua; // \mu_a
  double mus; // \mu_s
  double g;   // g
  double OneMinsGG; // 1-g*g
  double OneMinsG;  // 1-g
  double OneAddGG;  // 1+g*g
  double OneAddG;   // 1+g
  double TwoG;      // 2*g
  double pdwa;      // mua/(mua+mus)
  double RefIdx;    // the refraction index of the medium
  double MUAMUS;
  double IMUAMUS;
};

// Structure for source
struct TSource{
  // Type for the source
  // 1:  isotropic point source inside the specimen
  //     The source file will provide:
  //     X, Y, Z: the position of the point source
  //     The program will compute the Elem_Idx

  // 2:  isotropic region source inside a tetrahedron
  //     The source file will provide:
  //     Elem_Idx

  // 11: pencil beam source on the surface of the specimen 
  //     In this version:
  //     The source file will provide 
  //     (a) a point on a boundary triangle.
  //     (b) a vector of the source direct, which should be normal to the boundary triangle.
  //     (c) the element index, which contains the boundary triangle

  // 12: surface triangle region source
  //     In this type of source, we assume the light is normal to 
  //     the surface triangle and point in the object

  enum {ISOPOINT = 1, ISOREGION = 2, PENBEAM = 11, SURFTRIANGLE = 12};


  int           SourceType; 

  // The Element index 
  int           ElemIdx;

  // The surface triangle nodes index
  int           SurfTriNodes[3];
  int           SurfTriIdx;

  // The postion for isotropic point source
  // The point on the surface for pencil laser beam
  TNode         Position;

  // The incident angle for the pencil laser
  TNode         IncAngle;  

  // 
  long long int NumPhoton;

};

// Structure defined for absorption quantification
struct TPhotonInfo{
  int Time_Type; 
  // The last bit:  0 : internal element, 1: surface triangle
  // Other bits are time index
  int Idx;
  // Elem Index of Surface Triangle Index
  double LostWeight;
};
  
struct TPhoton{
  int    Cur_Elem;       // variable to store the current element index of the photon
  int    Cur_Med;        // variables to store the current med index of the photon
  double X, Y, Z, Path;  // 
  double UX, UY, UZ;     // direction
  double Weight;         // weight of the photon
  double Step;           // variables to store the step size of the photon
};

#endif
