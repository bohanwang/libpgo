#define MATRIX_SCALE3X4(a,scalar)\
  (a)[0] *= (scalar);\
  (a)[1] *= (scalar);\
  (a)[2] *= (scalar);\
  (a)[3] *= (scalar);\
  (a)[4] *= (scalar);\
  (a)[5] *= (scalar);\
  (a)[6] *= (scalar);\
  (a)[7] *= (scalar);\
  (a)[8] *= (scalar);\
  (a)[9] *= (scalar);\
  (a)[10] *= (scalar);\
  (a)[11] *= (scalar);

// C = A * B
#define MATRIX_MULTIPLY3X3X4(a,b,c)\
(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[4]+(a)[2]*(b)[8];\
(c)[1]=(a)[0]*(b)[1]+(a)[1]*(b)[5]+(a)[2]*(b)[9];\
(c)[2]=(a)[0]*(b)[2]+(a)[1]*(b)[6]+(a)[2]*(b)[10];\
(c)[3]=(a)[0]*(b)[3]+(a)[3]*(b)[7]+(a)[2]*(b)[11];\
(c)[4]=(a)[3]*(b)[0]+(a)[4]*(b)[4]+(a)[5]*(b)[8];\
(c)[5]=(a)[3]*(b)[1]+(a)[4]*(b)[5]+(a)[5]*(b)[9];\
(c)[6]=(a)[3]*(b)[2]+(a)[4]*(b)[6]+(a)[5]*(b)[10];\
(c)[7]=(a)[3]*(b)[3]+(a)[4]*(b)[7]+(a)[5]*(b)[11];\
(c)[8]=(a)[6]*(b)[0]+(a)[7]*(b)[4]+(a)[8]*(b)[8];\
(c)[9]=(a)[6]*(b)[1]+(a)[7]*(b)[5]+(a)[8]*(b)[9];\
(c)[10]=(a)[6]*(b)[2]+(a)[7]*(b)[6]+(a)[8]*(b)[10];\
(c)[11]=(a)[6]*(b)[3]+(a)[7]*(b)[7]+(a)[8]*(b)[11];

#define MATRIX_MULTIPLY3X4X3(a,b,c)\
(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[3]+(a)[2]*(b)[6]+(a)[3]*(b)[9];\
(c)[1]=(a)[0]*(b)[1]+(a)[1]*(b)[4]+(a)[2]*(b)[7]+(a)[3]*(b)[10];\
(c)[2]=(a)[0]*(b)[2]+(a)[1]*(b)[5]+(a)[2]*(b)[8]+(a)[3]*(b)[11];\
(c)[3]=(a)[4]*(b)[0]+(a)[5]*(b)[3]+(a)[6]*(b)[6]+(a)[7]*(b)[9];\
(c)[4]=(a)[4]*(b)[1]+(a)[5]*(b)[4]+(a)[6]*(b)[7]+(a)[7]*(b)[10];\
(c)[5]=(a)[4]*(b)[2]+(a)[5]*(b)[5]+(a)[6]*(b)[8]+(a)[7]*(b)[11];\
(c)[6]=(a)[8]*(b)[0]+(a)[9]*(b)[3]+(a)[10]*(b)[6]+(a)[11]*(b)[9];\
(c)[7]=(a)[8]*(b)[1]+(a)[9]*(b)[4]+(a)[10]*(b)[7]+(a)[11]*(b)[10];\
(c)[8]=(a)[8]*(b)[2]+(a)[9]*(b)[5]+(a)[10]*(b)[8]+(a)[11]*(b)[11];

#define MATRIX_MULTIPLY4X4X3(a,b,c)\
(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[3]+(a)[2]*(b)[6]+(a)[3]*(b)[9];\
(c)[1]=(a)[0]*(b)[1]+(a)[1]*(b)[4]+(a)[2]*(b)[7]+(a)[3]*(b)[10];\
(c)[2]=(a)[0]*(b)[2]+(a)[1]*(b)[5]+(a)[2]*(b)[8]+(a)[3]*(b)[11];\
(c)[3]=(a)[4]*(b)[0]+(a)[5]*(b)[3]+(a)[6]*(b)[6]+(a)[7]*(b)[9];\
(c)[4]=(a)[4]*(b)[1]+(a)[5]*(b)[4]+(a)[6]*(b)[7]+(a)[7]*(b)[10];\
(c)[5]=(a)[4]*(b)[2]+(a)[5]*(b)[5]+(a)[6]*(b)[8]+(a)[7]*(b)[11];\
(c)[6]=(a)[8]*(b)[0]+(a)[9]*(b)[3]+(a)[10]*(b)[6]+(a)[11]*(b)[9];\
(c)[7]=(a)[8]*(b)[1]+(a)[9]*(b)[4]+(a)[10]*(b)[7]+(a)[11]*(b)[10];\
(c)[8]=(a)[8]*(b)[2]+(a)[9]*(b)[5]+(a)[10]*(b)[8]+(a)[11]*(b)[11];\
(c)[9]=(a)[12]*(b)[0]+(a)[13]*(b)[3]+(a)[14]*(b)[6]+(a)[15]*(b)[9];\
(c)[10]=(a)[12]*(b)[1]+(a)[13]*(b)[4]+(a)[14]*(b)[7]+(a)[15]*(b)[10];\
(c)[11]=(a)[12]*(b)[2]+(a)[13]*(b)[5]+(a)[14]*(b)[8]+(a)[15]*(b)[11];

#define MATRIX_MULTIPLY4X3X3(a,b,c)\
(c)[0]=(a)[0]*(b)[0]+(a)[1]*(b)[3]+(a)[2]*(b)[6];\
(c)[1]=(a)[0]*(b)[1]+(a)[1]*(b)[4]+(a)[2]*(b)[7];\
(c)[2]=(a)[0]*(b)[2]+(a)[1]*(b)[5]+(a)[2]*(b)[8];\
(c)[3]=(a)[3]*(b)[0]+(a)[4]*(b)[3]+(a)[5]*(b)[6];\
(c)[4]=(a)[3]*(b)[1]+(a)[4]*(b)[4]+(a)[5]*(b)[7];\
(c)[5]=(a)[3]*(b)[2]+(a)[4]*(b)[5]+(a)[5]*(b)[8];\
(c)[6]=(a)[6]*(b)[0]+(a)[7]*(b)[3]+(a)[8]*(b)[6];\
(c)[7]=(a)[6]*(b)[1]+(a)[7]*(b)[4]+(a)[8]*(b)[7];\
(c)[8]=(a)[6]*(b)[2]+(a)[7]*(b)[5]+(a)[8]*(b)[8];\
(c)[9]=(a)[9]*(b)[0]+(a)[10]*(b)[3]+(a)[11]*(b)[6];\
(c)[10]=(a)[9]*(b)[1]+(a)[10]*(b)[4]+(a)[11]*(b)[7];\
(c)[11]=(a)[9]*(b)[2]+(a)[10]*(b)[5]+(a)[11]*(b)[8];

// C = A+B
#define MATRIX_ADD3X4(a,b)\
  (a)[0] += (b)[0];\
  (a)[1] += (b)[1];\
  (a)[2] += (b)[2];\
  (a)[3] += (b)[3];\
  (a)[4] += (b)[4];\
  (a)[5] += (b)[5];\
  (a)[6] += (b)[6];\
  (a)[7] += (b)[7];\
  (a)[8] += (b)[8];\
  (a)[9] += (b)[9];\
  (a)[10] += (b)[10];\
  (a)[11] += (b)[11];

#define MATRIX_TRANSPOSE4X3(A,B)\
(B)[0] = (A)[0]; (B)[1] = (A)[3]; (B)[2] = (A)[6]; (B)[3]=(A)[9];\
(B)[4] = (A)[1]; (B)[5] = (A)[4]; (B)[6] = (A)[7]; (B)[7]=(A)[10];\
(B)[8] = (A)[2]; (B)[9] = (A)[5]; (B)[10] = (A)[8]; (B)[11]=(A)[11];

#define MATRIX_TRANSPOSE3X4(A,B)\
(B)[0] = (A)[0]; (B)[1] = (A)[4]; (B)[2] = (A)[8]; \
(B)[3] = (A)[1]; (B)[4] = (A)[5]; (B)[5] = (A)[9];\
(B)[6] = (A)[2]; (B)[7] = (A)[6]; (B)[8] = (A)[10];\
(B)[9] = (A)[3]; (B)[10] = (A)[7]; (B)[11] = (A)[11];

#define MATRIX_SET4X3(a,b)\
  (a)[0] = (b)[0];\
  (a)[1] = (b)[1];\
  (a)[2] = (b)[2];\
  (a)[3] = (b)[3];\
  (a)[4] = (b)[4];\
  (a)[5] = (b)[5];\
  (a)[6] = (b)[6];\
  (a)[7] = (b)[7];\
  (a)[8] = (b)[8];\
  (a)[9] = (b)[9];\
  (a)[10] = (b)[10];\
  (a)[11] = (b)[11];

// c = a - b
#define MATRIX_SUBTRACT4X4(a,b,c)\
  (c)[0] = (a)[0] - (b)[0];\
  (c)[1] = (a)[1] - (b)[1];\
  (c)[2] = (a)[2] - (b)[2];\
  (c)[3] = (a)[3] - (b)[3];\
  (c)[4] = (a)[4] - (b)[4];\
  (c)[5] = (a)[5] - (b)[5];\
  (c)[6] = (a)[6] - (b)[6];\
  (c)[7] = (a)[7] - (b)[7];\
  (c)[8] = (a)[8] - (b)[8];\
  (c)[9] = (a)[9] - (b)[9];\
  (c)[10] = (a)[10] - (b)[10];\
  (c)[11] = (a)[11] - (b)[11];\
  (c)[12] = (a)[12] - (b)[12];\
  (c)[13] = (a)[13] - (b)[13];\
  (c)[14] = (a)[14] - (b)[14];\
  (c)[15] = (a)[15] - (b)[15];