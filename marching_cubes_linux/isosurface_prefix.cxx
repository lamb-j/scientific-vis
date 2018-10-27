#include <chrono>
#include <cstdio>
#include <cstdlib>
typedef std::chrono::high_resolution_clock Clock;

extern int triCase[256][16];

void AddTriangle(float *tl, int tIndex, float X1, float Y1, float Z1, float X2, float Y2, float Z2, float X3, float Y3, float Z3)
{
  tl[9*tIndex+0] = X1;
  tl[9*tIndex+1] = Y1;
  tl[9*tIndex+2] = Z1;
  tl[9*tIndex+3] = X2;
  tl[9*tIndex+4] = Y2;
  tl[9*tIndex+5] = Z2;
  tl[9*tIndex+6] = X3;
  tl[9*tIndex+7] = Y3;
  tl[9*tIndex+8] = Z3;
}

float InterpLin (float A, float B, float X, float fA, float fB) 
{
  return fA + ( (X - A) / (B - A) )*(fB - fA);
}

int GetNumberOfCells(const int *dims) 
{
  return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
}

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
  idx[0] = cellId%(dims[0]-1);
  idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
  idx[2] = cellId/((dims[0]-1)*(dims[1]-1));
}

void BoundingCubeForCell(const float *X, const float *Y, const float *Z, const int *dims,
    int cellId, float *bcube)
{
  if (cellId >= GetNumberOfCells(dims)) {
    printf("Bad cellID in BondingCubeForCell\n");
    return;
  }

  int idx[3];
  GetLogicalCellIndex(idx, cellId, dims);

  bcube[0] = X[idx[0]];
  bcube[1] = X[idx[0] + 1];

  bcube[2] = Y[idx[1]];
  bcube[3] = Y[idx[1] + 1];

  bcube[4] = Z[idx[2]];
  bcube[5] = Z[idx[2]+ 1];
}

int IdentifyCase(float *f_val, float iso_value) 
{
  int rv = 0;

  if (f_val[0] > iso_value) rv += 1;
  if (f_val[1] > iso_value) rv += 2;
  if (f_val[2] > iso_value) rv += 4;
  if (f_val[3] > iso_value) rv += 8;
  if (f_val[4] > iso_value) rv += 16;
  if (f_val[5] > iso_value) rv += 32;
  if (f_val[6] > iso_value) rv += 64;
  if (f_val[7] > iso_value) rv += 128;

  return rv;
};

float *GenerateIsosurfacePrefix(float *X, float *Y, float *Z, float *F, const int *dims, float iso_value, int *num_tris_total) 
{
  static int edge_to_vertex[12][2]  = 
  { {0,1},  {1,3},  {2,3},  {0,2},
    {4,5},  {5,7},  {6,7},  {4,6},
    {0,4},  {1,5},  {2,6},  {3,7} };

  static int edge_to_bcube[12][6]  = 
  { {0,1,2,2,4,4}, {1,1,2,3,4,4}, {0,1,3,3,4,4}, {0,0,2,3,4,4},  
    {0,1,2,2,5,5}, {1,1,2,3,5,5}, {0,1,3,3,5,5}, {0,0,2,3,5,5},  
    {0,0,2,2,4,5}, {1,1,2,2,4,5}, {0,0,3,3,4,5}, {1,1,3,3,4,5}};

  int nCells = GetNumberOfCells(dims);

  // Determine which cells contain triangles
  int *tri_occupied = (int *) malloc(sizeof(int) * nCells); 
  int *tri_prefix = (int *) malloc(sizeof(int) * nCells); 
  int *tri_icase = (int *) malloc(sizeof(int) * nCells); 
  tri_prefix[0] = 0;

  for (int cell_id = 0; cell_id < nCells; ++cell_id) {

    // get Field values
    int idx[3];
    GetLogicalCellIndex(idx, cell_id, dims);

    float f_val[8];
    int vx = idx[0], vy = idx[1], vz = idx[2];
    int dz = dims[0]*dims[1], dy = dims[0];

    f_val[0] = F[vz*dz + vy*dy + vx];               // v0
    f_val[1] = F[vz*dz + vy*dy + (vx+1)];           // v1
    f_val[2] = F[vz*dz + (vy+1)*dy + vx];           // v2
    f_val[3] = F[vz*dz + (vy+1)*dy + (vx+1)];       // v3
    f_val[4] = F[(vz+1)*dz + vy*dy + vx];           // v4
    f_val[5] = F[(vz+1)*dz + vy*dy + (vx+1)];       // v5
    f_val[6] = F[(vz+1)*dz + (vy+1)*dy + vx];       // v6
    f_val[7] = F[(vz+1)*dz + (vy+1)*dy + (vx+1)];   // v7

    int icase = IdentifyCase(f_val, iso_value);
    tri_icase[cell_id] = icase;

    int num_tris = 0;
    for (int i = 0; triCase[icase][i] != -1; i+=3) num_tris++;

    tri_occupied[cell_id] = num_tris;
    if (cell_id > 0) 
      tri_prefix[cell_id] = tri_prefix[cell_id-1] + tri_occupied[cell_id-1];
  }

  *num_tris_total = tri_occupied[nCells - 1] + tri_prefix[nCells - 1];

  // list of triangles
  float *tl = (float *) malloc(sizeof(float) * (*num_tris_total) * 9);

  // Generate Triangles
  for (int cell_id = 0; cell_id < nCells; ++cell_id) {

    if (tri_occupied[cell_id] == 0) continue;

    //int num_tris = tri_occupied[cell_id];
    int tri_index = tri_prefix[cell_id];
    int icase = tri_icase[cell_id];

    // get X, Y, Z values
    float bcube[6];   // X1, X2, Y1, Y2, Z1, Z2
    BoundingCubeForCell(X, Y, Z, dims, cell_id, bcube);

    // get Field values
    int idx[3];
    GetLogicalCellIndex(idx, cell_id, dims);

    float f_val[8];
    int vx = idx[0], vy = idx[1], vz = idx[2];
    int dz = dims[0]*dims[1], dy = dims[0];

    f_val[0] = F[vz*dz + vy*dy + vx];               // v0
    f_val[1] = F[vz*dz + vy*dy + (vx+1)];           // v1
    f_val[2] = F[vz*dz + (vy+1)*dy + vx];           // v2
    f_val[3] = F[vz*dz + (vy+1)*dy + (vx+1)];       // v3
    f_val[4] = F[(vz+1)*dz + vy*dy + vx];           // v4
    f_val[5] = F[(vz+1)*dz + vy*dy + (vx+1)];       // v5
    f_val[6] = F[(vz+1)*dz + (vy+1)*dy + vx];       // v6
    f_val[7] = F[(vz+1)*dz + (vy+1)*dy + (vx+1)];   // v7

    // draw triangles 
    for (int i = 0, t = 0; triCase[icase][i] != -1; i+=3, t++) {
      float tri[3][3]; // (x1, y1, z1), (x2, y2, z2), (x3, y3, z3)

      for (int e = 0; e < 3; ++e) {
        int edge = triCase[icase][i + e];

        // get field values for this edge
        float f1 = f_val[edge_to_vertex[edge][0]];
        float f2 = f_val[edge_to_vertex[edge][1]];

        // determine X, Y, Z values of edge 
        float x1 = bcube[ edge_to_bcube[edge][0] ];
        float x2 = bcube[ edge_to_bcube[edge][1] ];
        float y1 = bcube[ edge_to_bcube[edge][2] ];
        float y2 = bcube[ edge_to_bcube[edge][3] ];
        float z1 = bcube[ edge_to_bcube[edge][4] ];
        float z2 = bcube[ edge_to_bcube[edge][5] ];

        // interpolate X, Y values of new point
        tri[e][0] = InterpLin(f1, f2, iso_value, x1, x2);
        tri[e][1] = InterpLin(f1, f2, iso_value, y1, y2);
        tri[e][2] = InterpLin(f1, f2, iso_value, z1, z2);
      }

      AddTriangle(tl, tri_index + t, 
          tri[0][0], tri[0][1], tri[0][2], 
          tri[1][0], tri[1][1], tri[1][2],
          tri[2][0], tri[2][1], tri[2][2]); 
    }

  }

  return tl;
}
