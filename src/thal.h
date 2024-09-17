/*
 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2009,2010,
               2011,2012
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
 All rights reserved.

       This file is part of primer3 software suite.

       This software suite is is free software;
       you can redistribute it and/or modify it under the terms
       of the GNU General Public License as published by the Free
       Software Foundation; either version 2 of the License, or (at
       your option) any later version.

       This software is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this software (file gpl-2.0.txt in the source
       distribution); if not, write to the Free Software
       Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _THAL_H
#define _THAL_H


#include "limits.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "setjmp.h"
#include "ctype.h"
#include "math.h"
#include "unistd.h"
#include "float.h" /* ! mul ei ole float.h-d includes DBL_MAX */

namespace primer3thal
{

/* The maximum length of _one_ of the two sequences being aligned in a
   thermodynamic alignment. In other words, the length of one sequence
   must be <= THAL_MAX_ALIGN, but the other sequence can be longer.
   The rationale behind this value (60) is that this is the maxium
   reasonable length for nearest neighbor models. It is the maxium
   length at which we can restrict our model to only two states of
   melting: fully intact duplex or completely dissociated single
   strands. */
#ifndef THAL_MAX_ALIGN
#define THAL_MAX_ALIGN 60
#endif

/* The maxium length of the other sequence in a thermodynamic
   alignment. This value can be increased, though alignments against
   very long sequences will be quite slow. As of 2012-05-18, we only
   potentially see sequences longer this when checking for mispriming
   in the template ('max_template_mispriming') in libprimer3.c, which
   is really designed to find sites of ectopic primer very close (a
   few kilobases) from the location of the cadidate primer. */
#ifndef THAL_MAX_SEQ
#define THAL_MAX_SEQ   10000
#endif

#ifndef MIN_HRPN_LOOP
#define MIN_HRPN_LOOP 3 /*  minimum size of hairpin loop */
#endif

/* table where bp-s enthalpies, that retrieve to the most stable Tm, are saved */
#ifdef EnthalpyDPT
# undef EnthalpyDPT
#endif
#define EnthalpyDPT(i, j) enthalpyDPT[(j) + ((i-1)*len3) - (1)]

/* table where bp-s entropies, that retrieve to the most stable Tm, are saved */
#ifdef EntropyDPT
# undef EntropyDPT
#endif
#define EntropyDPT(i, j) entropyDPT[(j) + ((i-1)*len3) - (1)]

/* entropies of most stable hairpin terminal bp */
#ifndef SEND5
# define SEND5(i) send5[i]
#endif

/* enthalpies of most stable hairpin terminal bp */
#ifndef HEND5
# define HEND5(i) hend5[i]
#endif

#define bpIndx(a, b) BPI[a][b] /* for traceing matrix BPI */
#define atPenaltyS(a, b) atpS[a][b]
#define atPenaltyH(a, b) atpH[a][b]

#define SMALL_NON_ZERO 0.000001
#define DBL_EQ(X,Y) (((X) - (Y)) < (SMALL_NON_ZERO) ? (1) : (2)) /* 1 when numbers are equal */

#define isFinite(x) (x < _INFINITY / 2)

#define isPositive(x) ((x) > 0 ? (1) : (0))

/*** BEGIN CONSTANTS ***/
const double _INFINITY = 999999.0;
const double THAL_ERROR_SCORE = -_INFINITY;
static const double R = 1.9872; /* cal/Kmol */
static const double ILAS = (-300 / 310.15); /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
static const double ILAH = 0.0; /* Internal Loop EntHalpy Asymmetry correction */
static const double AT_H = 2200.0; /* AT penalty */
static const double AT_S = 6.9; /* AT penalty */
static const double MinEntropyCutoff = -2500.0; /* to filter out non-existing entropies */
static const double MinEntropy = -3224.0; /* initiation */
static const double G2 = 0.0; /* structures w higher G are considered to be unstabile */
const double ABSOLUTE_ZERO = 273.15;
const double TEMP_KELVIN = 310.15;
const int MAX_LOOP = 30; /* the maximum size of loop that can be calculated; for larger loops formula must be implemented */
const int MIN_LOOP = 0;

/* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
static const int BPI[5][5] =  {
     {0, 0, 0, 1, 0}, /* A, C, G, T, N; */
     {0, 0, 1, 0, 0},
     {0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0},
     {0, 0, 0, 0, 0}};

/* BEGIN TYPEDEFs */

typedef enum thal_alignment_type {
  thal_any = 1,
  thal_end1 = 2,
  thal_end2 = 3,
  thal_hairpin = 4,
} thal_alignment_type;

/* Structure for passing arguments to THermodynamic ALignment calculation */
typedef struct {
   int debug; /* if non zero, print debugging info to stderr */
   thal_alignment_type type; /* one of the
	      1 THAL_ANY, (by default)
	      2 THAL_END1,
	      3 THAL_END2,
	      4 THAL_HAIRPIN */
   int maxLoop;  /* maximum size of loop to consider; longer than 30 bp are not allowed */
   double mv; /* concentration of monovalent cations */
   double dv; /* concentration of divalent cations */
   double dntp; /* concentration of dNTP-s */
   double dna_conc; /* concentration of oligonucleotides */
   double temp; /* temperature from which hairpin structures will be calculated */
   int temponly; /* if non zero, print only temperature to stderr */
   int dimer; /* if non zero, dimer structure is calculated */
} thal_args;

/* Structure for receiving results from the thermodynamic alignment calculation */
typedef struct {
   double temp;
   int align_end_1;
   int align_end_2;
} thal_results;

/*** END OF TYPEDEFS ***/

/*** BEGIN STRUCTs ***/

struct triloop {
  char loop[5];
  double value; };

struct tetraloop {
  char loop[6];
  double value; };

struct tracer /* structure for tracebacku - unimolecular str */ {
  int i;
  int j;
  int mtrx; /* [0 1] EntropyDPT/EnthalpyDPT*/
  struct tracer* next;
};

static char* parampath = NULL; /* path to parameter files */
static int numTriloops; /* hairpin triloop penalties */
static int numTetraloops; /* hairpin tetraloop penalties */
static double atpS[5][5]; /* AT penalty */
static double atpH[5][5]; /* AT penalty */
static double *send5, *hend5; /* calc 5'  */
/* w/o init not constant anymore, cause for unimolecular and bimolecular foldings there are different values */
static double dplx_init_H; /* initiation enthalpy; for duplex 200, for unimolecular structure 0 */
static double dplx_init_S; /* initiation entropy; for duplex -5.7, for unimoleculat structure 0 */
static double saltCorrection; /* value calculated by saltCorrectS, includes correction for monovalent and divalent cations */
static double RC; /* universal gas constant multiplied w DNA conc - for melting temperature */
static double SHleft; /* var that helps to find str w highest melting temperature */
static int bestI, bestJ; /* starting position of most stable str */
static double* enthalpyDPT; /* matrix for values of enthalpy */
static double* entropyDPT; /* matrix for values of entropy */
static unsigned char *oligo1, *oligo2; /* inserted oligo sequenced */
static unsigned char *numSeq1, *numSeq2; /* same as oligo1 and oligo2 but converted to numbers */
static int len1, len2, len3; /* length of sequense 1 and 2 *//* 17.02.2009 int temponly;*/ /* print only temperature of the predicted structure */
static double dangleEntropies3[5][5][5]; /* thermodynamic paramteres for 3' dangling ends */
static double dangleEnthalpies3[5][5][5]; /* ther params for 3' dangling ends */
static double dangleEntropies5[5][5][5];  /* ther params for 5' dangling ends */
static double dangleEnthalpies5[5][5][5]; /* ther params for 5' dangling ends */
static double stackEntropies[5][5][5][5]; /* ther params for perfect match pairs */
static double stackEnthalpies[5][5][5][5]; /* ther params for perfect match pairs */
static double stackint2Entropies[5][5][5][5]; /*ther params for perfect match and internal mm */
static double stackint2Enthalpies[5][5][5][5]; /* ther params for perfect match and internal mm*/
static double interiorLoopEntropies[30]; /* interior loop params according to length of the loop */
static double bulgeLoopEntropies[30]; /* bulge loop params according to length of the loop */
static double hairpinLoopEntropies[30]; /* hairpin loop params accordint to length of the loop */
static double interiorLoopEnthalpies[30]; /* same as interiorLoopEntropies but values of entropy */
static double bulgeLoopEnthalpies[30]; /* same as bulgeLoopEntropies but values of entropy */
static double hairpinLoopEnthalpies[30]; /* same as hairpinLoopEntropies but values of entropy */
static double tstackEntropies[5][5][5][5]; /* ther params for terminal mismatches */
static double tstackEnthalpies[5][5][5][5]; /* ther params for terminal mismatches */
static double tstack2Entropies[5][5][5][5]; /* ther params for internal terminal mismatches */
static double tstack2Enthalpies[5][5][5][5]; /* ther params for internal terminal mismatches */
static struct triloop* triloopEntropies = NULL; /* ther penalties for given triloop seq-s */
static struct triloop* triloopEnthalpies = NULL; /* ther penalties for given triloop seq-s */
static struct tetraloop* tetraloopEntropies = NULL; /* ther penalties for given tetraloop seq-s */
static struct tetraloop* tetraloopEnthalpies = NULL; /* ther penalties for given tetraloop seq-s */


static int 
equal(double a, double b)
{
   return a == b;
}

static void 
strcatc(char* str, char c)
{
   str[strlen(str) + 1] = 0;
   str[strlen(str)] = c;
}

static int 
length_unsig_char(const unsigned char * str)
{
   int i = 0;
   while(*(str++)) {
      i++;
      if(i == INT_MAX)
	return -1;
   }
   return i;
}

static unsigned char 
str2int(char c)
{
   switch (c) {
    case 'A': case '0':
      return 0;
    case 'C': case '1':
      return 1;
    case 'G': case '2':
      return 2;
    case 'T': case '3':
      return 3;
   }
   return 4;
}


/* memory stuff */
static void* 
safe_realloc(void* ptr, size_t n)
{
   ptr = realloc(ptr, n);
   if (ptr == NULL) {
#ifdef DEBUG
      fputs("Error in realloc()\n", stderr);
#endif
   }
   return ptr;
}

static double* 
safe_recalloc(double* ptr, int m, int n)
{
  return (double*) safe_realloc(ptr, m * n * sizeof(double));
}

static void* 
safe_calloc(size_t m, size_t n)
{
   void* ptr;
   if (!(ptr = calloc(m, n))) {
#ifdef DEBUG
      fputs("Error in calloc()\n", stderr);
#endif
   }
   return ptr;
}

static void* 
safe_malloc(size_t n)
{
   void* ptr;
   if (!(ptr = malloc(n))) {
#ifdef DEBUG
      fputs("Error in malloc()\n", stderr);
#endif
   }
   return ptr;
}

static int 
max5(double a, double b, double c, double d, double e)
{
   if(a > b && a > c && a > d && a > e) return 1;
   else if(b > c && b > d && b > e) return 2;
   else if(c > d && c > e) return 3;
   else if(d > e) return 4;
   else return 5;
}

static void 
push(struct tracer** stack, int i, int j, int mtrx)
{
   struct tracer* new_top;
   new_top = (struct tracer*) safe_malloc(sizeof(struct tracer));
   new_top->i = i;
   new_top->j = j;
   new_top->mtrx = mtrx;
   new_top->next = *stack;
   *stack = new_top;
}

static void 
reverse(unsigned char *s)
{
   int i,j;
   char c;
   for (i = 0, j = length_unsig_char(s)-1; i < j; i++, j--) {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}

static double 
saltCorrectS (double mv, double dv, double dntp)
{
   if(dv<=0) dntp=dv;
   return 0.368*((log((mv+120*(sqrt(fmax(0.0, dv-dntp))))/1000)));
}


#define INIT_BUF_SIZE 1024

static char*
p3_read_line(FILE *file)
{
  static size_t ssz;
  static char *s = NULL;

  size_t remaining_size;
  char *p, *n;

  if (NULL == s) {
    ssz = INIT_BUF_SIZE;
    s = (char *) safe_malloc(ssz);
  }
  p = s;
  remaining_size = ssz;
  while (1) {
    if (fgets(p, remaining_size, file) == NULL) /* End of file. */
      return p == s ? NULL : s;

    if ((n = strchr(p, '\n')) != NULL) {
      *n = '\0';
      return s;
    }

    /* We did not get the whole line. */

    if (ssz >= INT_MAX / 2)
      ssz = INT_MAX;
    else {
      ssz *= 2;
    }
    s = (char *) safe_realloc(s, ssz);
    p = strchr(s, '\0');
    remaining_size = ssz - (p - s);
  }
}


/* These functions are needed as "inf" cannot be read on Windows directly */
static double 
readDouble(FILE *file)
{
  double result;
  char *line = p3_read_line(file);
  /* skip any spaces at beginning of the line */
  while (isspace(*line)) line++;
  if (!strncmp(line, "inf", 3))
    return _INFINITY;
  sscanf(line, "%lf", &result);
  return result;
}

/* Reads a line containing 4 doubles, which can be specified as "inf". */
static void
readLoop(FILE *file, double *v1, double *v2, double *v3)
{
  char *line = p3_read_line(file);
  char *p = line, *q;
  /* skip first number on the line */
  while (isspace(*p)) p++;
  while (isdigit(*p)) p++;
  while (isspace(*p)) p++;
  /* read second number */
  q = p;
  while (!isspace(*q)) q++;
  *q = '\0'; q++;
  if (!strcmp(p, "inf")) *v1 = _INFINITY;
  else sscanf(p, "%lf", v1);
  while (isspace(*q)) q++;
  /* read third number */
  p = q;
  while (!isspace(*p)) p++;
  *p = '\0'; p++;
  if (!strcmp(q, "inf")) *v2 = _INFINITY;
  else sscanf(q, "%lf", v2);
  while (isspace(*p)) p++;
  /* read last number */
  q = p;
  while (!isspace(*q) && (*q != '\0')) q++;
  *q = '\0';
  if (!strcmp(p, "inf")) *v3 = _INFINITY;
  else sscanf(p, "%lf", v3);
}

/* Reads a line containing a short string and a double, used for reading a triloop or tetraloop. */
static int
readTLoop(FILE *file, char *s, double *v, int triloop)
{
  char *line = p3_read_line(file);
  if (!line) return -1;
  char *p = line, *q;
  /* skip first spaces */
  while (isspace(*p)) p++;
  /* read the string */
  q = p;
  while (isalpha(*q)) q++;
  *q = '\0'; q++;
  if (triloop) {
    strncpy(s, p, 5);   /*triloop string has 5 characters*/
    s[5] = '\0';
  } else {
    strncpy(s, p, 6);   /*tetraloop string has 6 characters*/
    s[6] = '\0';
  }
  /* skip all spaces */
  while (isspace(*q)) q++;
  p = q;
  while (!isspace(*p) && (*p != '\0')) p++;
  *p = '\0';
  if (!strcmp(q, "inf")) *v = _INFINITY;
  else sscanf(q, "%lg", v);
  return 0;
} 


static FILE* 
openParamFile(const char* fname)
{
   FILE* file;
   char* paramdir;
   file = fopen(fname, "rt");
   if (!file) {
      paramdir = (char*) safe_malloc(strlen(parampath) + strlen(fname) + 1);
      strcpy(paramdir, parampath);
      strcat(paramdir, fname);
      if (!(file = fopen(paramdir, "rt"))) {
	fputs("Error opening params file\n", stderr);
      }
      free(paramdir);
   }
   return file;
}

static void 
getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5])
{
   int i, j, ii, jj;
   FILE *sFile, *hFile;
   sFile = openParamFile("stack.ds");
   hFile = openParamFile("stack.dh");
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
	 for (j = 0; j < 5; ++j) {
	    for (jj = 0; jj < 5; ++jj) {
	       if (i == 4 || j == 4 || ii == 4 || jj == 4) {
		  stackEntropies[i][ii][j][jj] = -1.0;
		  stackEnthalpies[i][ii][j][jj] = _INFINITY;
	       } else {
                  stackEntropies[i][ii][j][jj] = readDouble(sFile);
                  stackEnthalpies[i][ii][j][jj] = readDouble(hFile);
                  if (!isFinite(stackEntropies[i][ii][j][jj]) || !isFinite(stackEnthalpies[i][ii][j][jj])) {
		     stackEntropies[i][ii][j][jj] = -1.0;
		     stackEnthalpies[i][ii][j][jj] = _INFINITY;
		  }
	       }
	    }
	 }
      }
   }
   fclose(sFile);
   fclose(hFile);
}

static void 
getStackint2(double stackint2Entropies[5][5][5][5], double stackint2Enthalpies[5][5][5][5])
{
   int i, j, ii, jj;
   FILE *sFile, *hFile;
   sFile = openParamFile("stackmm.ds");
   hFile = openParamFile("stackmm.dh");
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
	 for (j = 0; j < 5; ++j) {
	    for (jj = 0; jj < 5; ++jj) {
	       if (i == 4 || j == 4 || ii == 4 || jj == 4) {
		  stackint2Entropies[i][ii][j][jj] = -1.0;
		  stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
	       } else {
                  stackint2Entropies[i][ii][j][jj] = readDouble(sFile);
                  stackint2Enthalpies[i][ii][j][jj] = readDouble(hFile);
                  if (!isFinite(stackint2Entropies[i][ii][j][jj]) || !isFinite(stackint2Enthalpies[i][ii][j][jj])) {
		     stackint2Entropies[i][ii][j][jj] = -1.0;
		     stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
		  }
	       }
	    }
	 }
      }
   }
   fclose(sFile);
   fclose(hFile);
}


static void 
getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5], double dangleEnthalpies5[5][5][5])
{
   int i, j, k;
   FILE *sFile, *hFile;
   sFile = openParamFile("dangle.ds");
   hFile = openParamFile("dangle.dh");
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
	  if (i == 4 || j == 4) {
	     dangleEntropies3[i][k][j] = -1.0;
	     dangleEnthalpies3[i][k][j] = _INFINITY;
	  } else if (k == 4) {
	     dangleEntropies3[i][k][j] = -1.0;
	     dangleEnthalpies3[i][k][j] = _INFINITY;
	  } else {
	     dangleEntropies3[i][k][j] = readDouble(sFile);
             dangleEnthalpies3[i][k][j] = readDouble(hFile);
             if(!isFinite(dangleEntropies3[i][k][j]) || !isFinite(dangleEnthalpies3[i][k][j])) {
		dangleEntropies3[i][k][j] = -1.0;
		dangleEnthalpies3[i][k][j] = _INFINITY;	     
             }
	  }
       }

   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
	  if (i == 4 || j == 4) {
	     dangleEntropies5[i][j][k] = -1.0;
	     dangleEnthalpies5[i][j][k] = _INFINITY;
	  } else if (k == 4) {
	     dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
	  } else {
	     dangleEntropies5[i][j][k] = readDouble(sFile);
             dangleEnthalpies5[i][j][k] = readDouble(hFile);
             if(!isFinite(dangleEntropies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k])) {
		dangleEntropies5[i][j][k] = -1.0;
		dangleEnthalpies5[i][j][k] = _INFINITY;
	     }
	  }
       }
   fclose(sFile);
   fclose(hFile);
}

static void 
getLoop(double hairpinLoopEntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropies[30],
	double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30])
{
   int k;
   FILE *sFile, *hFile;
   sFile = openParamFile("loops.ds");
   hFile = openParamFile("loops.dh");
   for (k = 0; k < 30; ++k) {
      readLoop(sFile, &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k]);
      readLoop(hFile, &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k]);
   }
   fclose(sFile);
   fclose(hFile);
}

static void 
getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5])
{
   int i1, j1, i2, j2;
   FILE *sFile, *hFile;
   sFile = openParamFile("tstack_tm_inf.ds");
   hFile = openParamFile("tstack.dh");
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
	 for (j2 = 0; j2 < 5; ++j2)
	   if (i1 == 4 || j1 == 4) {
	      tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
	      tstackEntropies[i1][i2][j1][j2] = -1.0;
	   } else if (i2 == 4 || j2 == 4) {
	      tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
	      tstackEnthalpies[i1][i2][j1][j2] = 0.0;
	   } else {
              tstackEntropies[i1][i2][j1][j2] = readDouble(sFile);
              tstackEnthalpies[i1][i2][j1][j2] = readDouble(hFile);
              if (!isFinite(tstackEntropies[i1][i2][j1][j2]) || !isFinite(tstackEnthalpies[i1][i2][j1][j2])) {
		 tstackEntropies[i1][i2][j1][j2] = -1.0;
		 tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
	      }
	   }
   fclose(sFile);
   fclose(hFile);
}

static void 
getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5])
{

   int i1, j1, i2, j2;
   FILE *sFile, *hFile;
   sFile = openParamFile("tstack2.ds");
   hFile = openParamFile("tstack2.dh");
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
	 for (j2 = 0; j2 < 5; ++j2)
	   if (i1 == 4 || j1 == 4)  {
	      tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
	      tstack2Entropies[i1][i2][j1][j2] = -1.0;
	   } else if (i2 == 4 || j2 == 4) {
	      tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
	      tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
	   } else {
              tstack2Entropies[i1][i2][j1][j2] = readDouble(sFile);
              tstack2Enthalpies[i1][i2][j1][j2] = readDouble(hFile);
              if (!isFinite(tstack2Entropies[i1][i2][j1][j2]) || !isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
		 tstack2Entropies[i1][i2][j1][j2] = -1.0;
		 tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
	      }
	   }
   fclose(sFile);
   fclose(hFile);
}

static void 
getTriloop(struct triloop** triloopEntropies, struct triloop** triloopEnthalpies, int* num)
{
   FILE *sFile, *hFile;
   int i, size;
   double value;
   sFile = openParamFile("triloop.ds");
   *num = 0;
   size = 16;
   *triloopEntropies = (struct triloop*) safe_calloc(16, sizeof(struct triloop));
   while (readTLoop(sFile, (*triloopEntropies)[*num].loop, &value, 1) != -1) {
      for (i = 0; i < 5; ++i)
	(*triloopEntropies)[*num].loop[i] = str2int((*triloopEntropies)[*num].loop[i]);
      (*triloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size)	{
	 size *= 2;
	 *triloopEntropies = (struct triloop*) safe_realloc(*triloopEntropies, size * sizeof(struct triloop));
      }
   }
   *triloopEntropies = (struct triloop*) safe_realloc(*triloopEntropies, *num * sizeof(struct triloop));

   fclose(sFile);

   hFile = openParamFile("triloop.dh");
   *num = 0;
   size = 16;
   *triloopEnthalpies = (struct triloop*) safe_calloc(16, sizeof(struct triloop));
   while (readTLoop(hFile, (*triloopEnthalpies)[*num].loop, &value, 1) != -1) {
      for (i = 0; i < 5; ++i)
	(*triloopEnthalpies)[*num].loop[i] = str2int((*triloopEnthalpies)[*num].loop[i]);
      (*triloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *triloopEnthalpies = (struct triloop*) safe_realloc(*triloopEnthalpies, size * sizeof(struct triloop));
      }
   }
   *triloopEnthalpies = (struct triloop*) safe_realloc(*triloopEnthalpies, *num * sizeof(struct triloop));

   fclose(hFile);
}

static void 
getTetraloop(struct tetraloop** tetraloopEntropies, struct tetraloop** tetraloopEnthalpies, int* num)
{

   FILE *sFile, *hFile;
   int i, size;
   double value;
   sFile = openParamFile("tetraloop.ds");
   *num = 0;
   size = 16;
   *tetraloopEntropies = (struct tetraloop*) safe_calloc(16, sizeof(struct tetraloop));
   while (readTLoop(sFile, (*tetraloopEntropies)[*num].loop, &value, 0) != -1) {
      for (i = 0; i < 6; ++i)
	(*tetraloopEntropies)[*num].loop[i] = str2int((*tetraloopEntropies)[*num].loop[i]);
      (*tetraloopEntropies)[*num].value = value;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *tetraloopEntropies = (struct tetraloop*) safe_realloc(*tetraloopEntropies, size * sizeof(struct tetraloop));
      }
   }
   *tetraloopEntropies = (struct tetraloop*) safe_realloc(*tetraloopEntropies, *num * sizeof(struct tetraloop));
   fclose(sFile);

   hFile = openParamFile("tetraloop.dh");
   *num = 0;
   size = 16;
   *tetraloopEnthalpies = (struct tetraloop*) safe_calloc(16, sizeof(struct tetraloop));
   while (readTLoop(hFile, (*tetraloopEnthalpies)[*num].loop, &value, 0) != -1) {
      for (i = 0; i < 6; ++i)
	(*tetraloopEnthalpies)[*num].loop[i] = str2int((*tetraloopEnthalpies)[*num].loop[i]);
      (*tetraloopEnthalpies)[*num].value = value;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *tetraloopEnthalpies = (struct tetraloop*) safe_realloc(*tetraloopEnthalpies, size * sizeof(struct tetraloop));
      }
   }
   *tetraloopEnthalpies = (struct tetraloop*) safe_realloc(*tetraloopEnthalpies, *num * sizeof(struct tetraloop));

   fclose(hFile);
}

static void 
tableStartATS(double atp_value, double atpS[5][5]) {
  int i, j;
  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      atpS[i][j] = 0.00000000001;
  atpS[0][3] = atpS[3][0] = atp_value;
}


static void 
tableStartATH(double atp_value, double atpH[5][5]) {
  int i, j;
  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      atpH[i][j] = 0.0;
  atpH[0][3] = atpH[3][0] = atp_value;
}

static int 
comp3loop(const void* loop1, const void* loop2)
{

     int i;
     const unsigned char* h1 = (const unsigned char*) loop1;
     const struct triloop *h2 = (const struct triloop*) loop2;

     for (i = 0; i < 5; ++i)
         if (h1[i] < h2->loop[i])
             return -1;
       else if (h1[i] > h2->loop[i])
           return 1;

     return 0;
}

static int 
comp4loop(const void* loop1, const void* loop2)
{
   int i;
   const unsigned char* h1 = (const unsigned char*) loop1;
   const struct tetraloop *h2 = (const struct tetraloop*) loop2;

   for (i = 0; i < 6; ++i)
     if (h1[i] < h2->loop[i])
       return -1;
   else if (h1[i] > h2->loop[i])
     return 1;

   return 0;
}


static void 
initMatrix()
{
   int i, j;
   for (i = 1; i <= len1; ++i) {
      for (j = 1; j <= len2; ++j) {
	 if (bpIndx(numSeq1[i], numSeq2[j]) == 0)  {
	    EnthalpyDPT(i, j) = _INFINITY;
	    EntropyDPT(i, j) = -1.0;
	 } else {
	    EnthalpyDPT(i, j) = 0.0;
	    EntropyDPT(i, j) = MinEntropy;
	 }
      }
   }
}

static void 
initMatrix2()
{
   int i, j;
   for (i = 1; i <= len1; ++i)
     for (j = i; j <= len2; ++j)
       if (j - i < MIN_HRPN_LOOP + 1 || (bpIndx(numSeq1[i], numSeq1[j]) == 0)) {
	  EnthalpyDPT(i, j) = _INFINITY;
	  EntropyDPT(i, j) = -1.0;
       } else {
	  EnthalpyDPT(i, j) = 0.0;
	  EntropyDPT(i, j) = MinEntropy;
       }

}

static void 
LSH(int i, int j, double* EntropyEnthalpy)
{
   double S1, H1, T1, G1;
   double S2, H2, T2, G2;
   S1 = S2 = -1.0;
   H1 = H2 = -_INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {
      EntropyDPT(i, j) = -1.0;
      EnthalpyDPT(i, j) = _INFINITY;
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(!isFinite(H1) || G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }
   /** If there is two dangling ends at the same end of duplex **/
   if((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1 ) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
	dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
	dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
	 H2 = _INFINITY;
	 S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1<T2 && G2<0) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else if(G2<0) {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   } else if ((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
	 H2 = _INFINITY;
	 S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1<T2 && G2<0) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else if (G2<0){
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   } else if ((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
	 H2 = _INFINITY;
	 S2 = -1.0;
     G2 = 1.0;     
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1 < T2  && G2<0) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else if(G2<0) {
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   if(isFinite(H1)) {
      if(T1 < T2) {
	 EntropyEnthalpy[0] = S2;
	 EntropyEnthalpy[1] = H2;
      } else {
	 EntropyEnthalpy[0] = S1;
	 EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}

static void 
RSH(int i, int j, double* EntropyEnthalpy)
{
   double G1, G2;
   double S1, S2;
   double H1, H2;
   double T1, T2;
   S1 = S2 = -1.0;
   H1 = H2 = _INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(!isFinite(H1) || G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }
   
   if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
	dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
	dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
    G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
	 H2 = _INFINITY;
	 S2 = -1.0;
     G2 = 1.0;
      }
      
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1 < T2 && G2<0) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else if(G2<0){
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }

   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2 >0) {
	 H2 = _INFINITY;
	 S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1 < T2 && G2<0) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else if(G2<0){
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }

   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
	 H2 = _INFINITY;
	 S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
	 T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
	 if(T1 < T2 && G2<0) {
	    S1 = S2;
	    H1 = H2;
	    T1 = T2;
	 }
      } else if (G2<0){
	 S1 = S2;
	 H1 = H2;
	 T1 = T2;
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   if(isFinite(H1)) {
      if(T1 < T2) {
	 EntropyEnthalpy[0] = S2;
	 EntropyEnthalpy[1] = H2;
      } else {
	 EntropyEnthalpy[0] = S1;
	 EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}

static double 
Ss(int i, int j, int k)
{
   if(k==2) {
      if (i >= j)
	return -1.0;
      if (i == len1 || j == len2 + 1)
	return -1.0;

      if (i > len1)
	i -= len1;
      if (j > len2)
	j -= len2;
      return stackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
   } else {
      return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}


static double 
Hs(int i, int j, int k)
{
   if(k==2) {
      if (i >= j)
	return _INFINITY;
      if (i == len1 || j == len2 + 1)
	return _INFINITY;

      if (i > len1)
	i -= len1;
      if (j > len2)
	j -= len2;
      if(isFinite(stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]])) {
	 return stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
      } else {
	 return _INFINITY;
      }
   } else {
      return stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}



static void 
maxTM(int i, int j)
{
   double S0, S1;
   double H0, H1;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double));
   double T1 = -_INFINITY;
   S0 = EntropyDPT(i, j);
   H0 = EnthalpyDPT(i, j);
   RSH(i,j,SH);
   double T0 = (H0 + dplx_init_H + SH[1]) /(S0 + dplx_init_S + SH[0] + RC); /* at current position */
   if(isFinite(EnthalpyDPT(i - 1, j - 1)) && isFinite(Hs(i - 1, j - 1, 1))) {
      S1 = (EntropyDPT(i - 1, j - 1) + Ss(i - 1, j - 1, 1));
      H1 = (EnthalpyDPT(i - 1, j - 1) + Hs(i - 1, j - 1, 1));
      T1 = (H1 + dplx_init_H + SH[1]) /(S1 + dplx_init_S + SH[0] + RC);
   } else {
      S1 = -1.0;
      H1 = _INFINITY;
      T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
   }
   
   if(S1 < MinEntropyCutoff) {
      /* to not give dH any value if dS is unreasonable */
      S1 = MinEntropy;
      H1 = 0.0;
   }
   if(S0 < MinEntropyCutoff) {
      /* to not give dH any value if dS is unreasonable */
      S0 = MinEntropy;
      H0 = 0.0;
   }
   if(T1 > T0) { 
      EntropyDPT(i, j) = S1;
      EnthalpyDPT(i, j) = H1;
   } else if(T0 >= T1) {
      EntropyDPT(i, j) = S0;
      EnthalpyDPT(i, j) = H0;
   }
   free(SH);
}

static void 
maxTM2(int i, int j)
{
   double S0, S1;
   double H0, H1;
   double T1 = -_INFINITY;
   S0 = EntropyDPT(i, j);
   H0 = EnthalpyDPT(i, j);
   double T0 = (H0 + dplx_init_H) /(S0 + dplx_init_S + RC);
   if(isFinite(EnthalpyDPT(i, j))) {
      S1 = (EntropyDPT(i + 1, j - 1) + Ss(i, j, 2));
      H1 = (EnthalpyDPT(i + 1, j - 1) + Hs(i, j, 2));
   } else {
      S1 = -1.0;
      H1 = _INFINITY;
   }
   T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
   if(S1 < MinEntropyCutoff) {
      S1 = MinEntropy;
      H1 = 0.0;
   }
   if(S0 < MinEntropyCutoff) {
      S0 = MinEntropy;
      H0 = 0.0;
   }

   if(T1 > T0) {
      EntropyDPT(i, j) = S1;
      EnthalpyDPT(i, j) = H1;
   } else {
      EntropyDPT(i, j) = S0;
      EnthalpyDPT(i, j) = H0;
   }
}

static void 
calc_bulge_internal(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback)
{
   int loopSize1, loopSize2, loopSize;
   double S,H,G1,G2;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double));
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   S = -1.0;
   H = _INFINITY;
   loopSize1 = ii - i - 1;
   loopSize2 = jj - j - 1;
#ifdef DEBUG
   if (ii <= i){
      fputs("Error in calc_bulge_internal(): ii is not greater than i\n", stderr);
   }
   if (jj <= j)
     fputs("Error in calc_bulge_internal(): jj is not greater than j\n", stderr);
#endif

#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      free(SH);
      return;
   }
#endif
   loopSize = loopSize1 + loopSize2-1;
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
					      the intervening nn-pair must be added */

	 if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
	    H = bulgeLoopEnthalpies[loopSize] +
	      stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
	    S = bulgeLoopEntropies[loopSize] +
	      stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
	 }
	 if(isPositive(H) || isPositive(S)){
	    H = _INFINITY;
	    S = -1.0;
	 }
	 H += EnthalpyDPT(i, j);
	 S += EntropyDPT(i, j);
	 if(!isFinite(H)) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 RSH(ii,jj,SH);
	 G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
	 G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*((EntropyDPT(ii, jj)+SH[0]));
	 if((G1< G2) || (traceback==1)) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

	 H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
	 H += EnthalpyDPT(i, j);

	 S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
	 S += EntropyDPT(i, j);
	 if(!isFinite(H)) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 if(isPositive(H) && isPositive(S)){ 
	    H = _INFINITY;
	    S = -1.0;
	 }
	
     RSH(ii,jj,SH);
	 G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
	 G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
	 if(G1< G2 || (traceback==1)){
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
	 
      }
   } else if (loopSize1 == 1 && loopSize2 == 1) {
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
	stackint2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      S += EntropyDPT(i, j);

      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
	stackint2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]];
      H += EnthalpyDPT(i, j);
      if(!isFinite(H)) {
	 H = _INFINITY;
	 S = -1.0;
      }
      if(isPositive(H) && isPositive(S)){
	 H = _INFINITY;
	 S = -1.0;
      }    
     RSH(ii,jj,SH);
	 G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
      G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
	   if((G1< G2) || traceback==1) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      free(SH);
      return;
   } else { /* only internal loops */
      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
	tstackEnthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]
	+ (ILAH * abs(loopSize1 - loopSize2));
      H += EnthalpyDPT(i, j);

      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
	tstackEntropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      S += EntropyDPT(i, j);
      if(!isFinite(H)) {
	 H = _INFINITY;
	 S = -1.0;
      }
   if(isPositive(H) && isPositive(S)){ 
	 H = _INFINITY;
	 S = -1.0;
      }
     RSH(ii,jj,SH);
     G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
     G2 = EnthalpyDPT(ii, jj)+SH[1]-TEMP_KELVIN*(EntropyDPT(ii, jj)+SH[0]);
     if((G1< G2) || (traceback==1)){
	     EntropyEnthalpy[0] = S;
	     EntropyEnthalpy[1] = H;
      }
   }
   free(SH);
   return;
}

static void 
calc_bulge_internal2(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback, int maxLoop)
{
   int loopSize1, loopSize2, loopSize;
   double T1, T2;
   double S,H;
   /* int N, N_loop; Triinu, please review */
   T1 = T2 = -_INFINITY;
   S = MinEntropy;
   H = 0.0;
   loopSize1 = ii - i - 1;
   loopSize2 = j - jj - 1;
   if (loopSize1 + loopSize2 > maxLoop) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   /* Triinu, please review the statements below. */
   /* if(i < (len1 -j)) { */
     /* N  = i; */
      /* N_loop = (i - 1); */
   /* } else { */
     /* N = len1-j;  */
      /* N_loop = len1 - j - 1; */
   /* } */
#ifdef DEBUG
   if (ii <= i)
     fputs("Error in calc_bulge_internal(): ii isn't greater than i\n", stderr);
   if (jj >= j)
     fputs("Error in calc_bulge_internal(): jj isn't less than j\n", stderr);
   if (ii >= jj)
     fputs("Error in calc_bulge_internal(): jj isn't greater than ii\n", stderr);

   if ((i <= len1 && len1 < ii) || (jj <= len2 && len2 < j))  {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
#endif

#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      return;
   }
#endif
#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      return;
   }
#endif

#ifdef DEBUG
   if (i > len1)
     i -= len1;
   if (ii > len1)
     ii -= len1;
   if (j > len2)
     j -= len2;
   if (jj > len2)
     jj -= len2;
#endif
   loopSize = loopSize1 + loopSize2 -1; /* for indx only */
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
					      the intervening nn-pair must be added */
	 if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
	    H = bulgeLoopEnthalpies[loopSize] +
	      stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
	    S = bulgeLoopEntropies[loopSize] +
	      stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
	 }
	 if(traceback!=1) {
	    H += EnthalpyDPT(ii, jj); /* bulge koos otsaga, st bulge i,j-ni */
	    S += EntropyDPT(ii, jj);
	 }

	 if(!isFinite(H)) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	
	 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
	 T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);

	 if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }

      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

	 H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii], numSeq2[jj]);
	 if(traceback!=1)
	   H += EnthalpyDPT(ii, jj);

	 S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii], numSeq2[jj]);
	 if(traceback!=1)
	   S += EntropyDPT(ii, jj);
	 if(!isFinite(H)) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	  
	 T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
	 T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);

	 if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      }
   } /* end of calculating bulges */
   else if (loopSize1 == 1 && loopSize2 == 1) {
      /* mismatch nearest neighbor parameters */

      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
	stackint2Entropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
	S += EntropyDPT(ii, jj);

      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
	stackint2Enthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]];
      if(traceback!=1)
	H += EnthalpyDPT(ii, jj);
      if(!isFinite(H)) {
	 H = _INFINITY;
	 S = -1.0;
      }
    
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);
      if((DBL_EQ(T1,T2) == 2) || traceback) {
	 if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
	    EntropyEnthalpy[0] = S;
	    EntropyEnthalpy[1] = H;
	 }
      }
      return;
   } else { /* only internal loops */

      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
	tstackEnthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]]
	+ (ILAH * abs(loopSize1 - loopSize2));
      if(traceback!=1)
	H += EnthalpyDPT(ii, jj);

      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
	tstackEntropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]] + (ILAS * abs(loopSize1 - loopSize2));
      if(traceback!=1)
	S += EntropyDPT(ii, jj);
      if(!isFinite(H)) {
	 H = _INFINITY;
	 S = -1.0;
      }
    
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);
      if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
	 EntropyEnthalpy[0] = S;
	 EntropyEnthalpy[1] = H;
      }
   }
   return;
}


static void 
fillMatrix(int maxLoop)
{
   int d, i, j, ii, jj;
   double* SH;

   SH = (double*) safe_malloc(2 * sizeof(double));
   for (i = 1; i <= len1; ++i) {
      for (j = 1; j <= len2; ++j) {
	 if(isFinite(EnthalpyDPT(i, j))) { /* if finite */
	    SH[0] = -1.0;
	    SH[1] = _INFINITY;
	    LSH(i,j,SH);
	    if(isFinite(SH[1])) {
	       EntropyDPT(i,j) = SH[0];
	       EnthalpyDPT(i,j) = SH[1];
	    }
	    if (i > 1 && j > 1) {
	       maxTM(i, j); /* stack: sets EntropyDPT(i, j) and EnthalpyDPT(i, j) */
	       for(d = 3; d <= maxLoop + 2; d++) { /* max=30, length over 30 is not allowed */
		  ii = i - 1;
		  jj = - ii - d + (j + i);
		  if (jj < 1) {
		     ii -= abs(jj-1);
		     jj = 1;
		  }
		  for (; ii > 0 && jj < j; --ii, ++jj) {
		     if (isFinite(EnthalpyDPT(ii, jj))) {
			SH[0] = -1.0;
			SH[1] = _INFINITY;
			calc_bulge_internal(ii, jj, i, j, SH,0);
			if(SH[0] < MinEntropyCutoff) {
			   /* to not give dH any value if dS is unreasonable */
			   SH[0] = MinEntropy;
			   SH[1] = 0.0;
			}
			if(isFinite(SH[1])) {
			   EnthalpyDPT(i, j) = SH[1];
			   EntropyDPT(i, j) = SH[0];
			}
		     }
		  }
	       }
	    } /* if */
	 }
      } /* for */
   } /* for */
   free(SH);
}

static void 
CBI(int i, int j, double* EntropyEnthalpy, int traceback, int maxLoop)
{
   int d, ii, jj;
   for (d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop; --d)
     for (ii = i + 1; ii < j - d && ii <= len1; ++ii) {
	jj = d + ii;
	if(traceback==0) {
	   EntropyEnthalpy[0] = -1.0;
	   EntropyEnthalpy[1] = _INFINITY;
	}
	if (isFinite(EnthalpyDPT(ii, jj)) && isFinite(EnthalpyDPT(i, j))) {
	   calc_bulge_internal2(i, j, ii, jj, EntropyEnthalpy, traceback,maxLoop);
	   if(isFinite(EntropyEnthalpy[1])) {
	      if(EntropyEnthalpy[0] < MinEntropyCutoff) {
		 EntropyEnthalpy[0] = MinEntropy;
		 EntropyEnthalpy[1] = 0.0;
	      }
	      if(traceback==0) {
		 EnthalpyDPT(i, j) = EntropyEnthalpy[1];
		 EntropyDPT(i, j) = EntropyEnthalpy[0];
	      }
	   }
	}
     }
   return;
}


static void 
calc_hairpin(int i, int j, double* EntropyEnthalpy, int traceback)
{
   int loopSize = j - i - 1;
   double G1, G2;
   G1 = G2 = -_INFINITY;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double));
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   if(loopSize < MIN_HRPN_LOOP) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      free(SH);
      return;
   }
   if (i <= len1 && len2 < j) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      free(SH);
      return;
   } else if (i > len2) {
      i -= len1;
      j -= len2;
   }
   if(loopSize <= 30) {
      EntropyEnthalpy[1] = hairpinLoopEnthalpies[loopSize - 1];
      EntropyEnthalpy[0] = hairpinLoopEntropies[loopSize - 1];
   } else {
      EntropyEnthalpy[1] = hairpinLoopEnthalpies[29];
      EntropyEnthalpy[0] = hairpinLoopEntropies[29];
   }

   if (loopSize > 3) { /* for loops 4 bp and more in length, terminal mm are accounted */
      EntropyEnthalpy[1] += tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
      EntropyEnthalpy[0] += tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
   } else if(loopSize == 3){ /* for loops 3 bp in length at-penalty is considered */
      EntropyEnthalpy[1] += atPenaltyH(numSeq1[i], numSeq1[j]);
      EntropyEnthalpy[0] += atPenaltyS(numSeq1[i], numSeq1[j]);
   }

   if (loopSize == 3) {	 /* closing AT-penalty (+), triloop bonus, hairpin of 3 (+) */
      struct triloop* loop;
      if (numTriloops) {
	 if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEnthalpies, numTriloops, sizeof(struct triloop), comp3loop)))
	   EntropyEnthalpy[1] += loop->value;
	 if ((loop = (struct triloop*) bsearch(numSeq1 + i, triloopEntropies, numTriloops, sizeof(struct triloop), comp3loop)))
	   EntropyEnthalpy[0] += loop->value;
      }
   } else if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
      struct tetraloop* loop;
      if (numTetraloops) {
	 if ((loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEnthalpies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
	    EntropyEnthalpy[1] += loop->value;
	 }
	 if ((loop = (struct tetraloop*) bsearch(numSeq1 + i, tetraloopEntropies, numTetraloops, sizeof(struct tetraloop), comp4loop))) {
	    EntropyEnthalpy[0] += loop->value;
	 }
      }
   }
   if(!isFinite(EntropyEnthalpy[1])) {
      EntropyEnthalpy[1] = _INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   if(isPositive(EntropyEnthalpy[1]) && isPositive(EntropyEnthalpy[0]) && (!isPositive(EnthalpyDPT(i, j)) || !isPositive(EntropyDPT(i, j)))) { /* if both, S and H are positive */
      EntropyEnthalpy[1] = _INFINITY;
      EntropyEnthalpy[0] = -1.0;
   }
   RSH(i,j,SH);
   G1 = EntropyEnthalpy[1]+SH[1] -TEMP_KELVIN*(EntropyEnthalpy[0]+SH[0]);
   G2 = EnthalpyDPT(i, j)+SH[1] -TEMP_KELVIN*(EntropyDPT(i, j)+SH[0]);
     if(G2 < G1 && traceback == 0) {
      EntropyEnthalpy[0] = EntropyDPT(i, j);
      EntropyEnthalpy[1] = EnthalpyDPT(i, j);
   }
   free(SH);
   return;
}

static double 
Sd5(int i, int j)
{
   return dangleEntropies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
}

static double 
Hd5(int i, int j)
{
   return dangleEnthalpies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]];
}

static double 
Sd3(int i, int j)
{
   return dangleEntropies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
}

static double 
Hd3(int i, int j)
{
   return dangleEnthalpies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]];
}

static double 
Ststack(int i, int j)
{
   return tstack2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
}

static double 
Htstack(int i, int j)
{ /* e.g AG_TC 210 */
   return tstack2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]];
}


static double 
END5_1(int i,int hs)
{
   int k;
   double max_tm; /* energy min */
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   S_max = S = -1.0;
   T1 = T2 = -_INFINITY;
   max_tm = -_INFINITY;
   for(k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
	 H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
	 S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
	 if(!isFinite(H) || H > 0 || S > 0) { /* H and S must be greater than 0 to avoid BS */
	    H = _INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
	 H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
	 S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
	 if(S > MinEntropyCutoff) {
	    H_max = H;
	    S_max = S;
	    max_tm = T1;
	 }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double 
END5_2(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
	 H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
	 S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
	 H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i);
	 S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
	 if(S > MinEntropyCutoff) {
	    H_max = H;
	    S_max = S;
	    max_tm = T1;
	 }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double 
END5_3(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
	 H = HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
	 S = SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
	 H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1);
	 S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
	 if(S > MinEntropyCutoff) {
	    H_max = H;
	    S_max = S;
	    max_tm = T1;
	 }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}

static double 
END5_4(int i,int hs)
{
   int k;
   double max_tm;
   double T1, T2;
   double H, S;
   double H_max, S_max;
   H_max = H = _INFINITY;
   T1 = T2 = max_tm = -_INFINITY;
   S_max = S = -1.0;
   for(k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k) {
      T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
      T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
      if(T1 >= T2) {
	 H = HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
	 S = SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
      } else {
	 H = 0 + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1);
	 S = 0 + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1);
	 if(!isFinite(H) || H > 0 || S > 0) {
	    H = _INFINITY;
	    S = -1.0;
	 }
	 T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
      }
      if(max_tm < T1) {
	 if(S > MinEntropyCutoff) {
	    H_max = H;
	    S_max = S;
	    max_tm = T1;
	 }
      }
   }
   if (hs == 1) return H_max;
   return S_max;
}


static void 
calc_terminal_bp(double temp) { /* compute exterior loop */
   int i;
   int max;
   double G;
   SEND5(0) = SEND5(1) = -1.0;
   HEND5(0) = HEND5(1) = _INFINITY;
   for(i = 2; i<=(len1); i++) {
      SEND5(i) = MinEntropy;
      HEND5(i) = 0;
   }

   /* adding terminal penalties to 3' end and to 5' end */
   for(i = 2; i <= len1; ++i) {
      double T1 = (HEND5(i - 1) + dplx_init_H) / (SEND5(i - 1) + dplx_init_S + RC);
      double T2 = (END5_1(i,1) + dplx_init_H) / (END5_1(i,2) + dplx_init_S + RC);
      double T3 = (END5_2(i,1) + dplx_init_H) / (END5_2(i,2) + dplx_init_S + RC);
      double T4 = (END5_3(i,1) + dplx_init_H) / (END5_3(i,2) + dplx_init_S + RC);
      double T5 = (END5_4(i,1) + dplx_init_H) / (END5_4(i,2) + dplx_init_S + RC);
      max = max5(T1,T2,T3,T4,T5);
      switch (max) {
       case 1:
	 SEND5(i) = SEND5(i - 1);
	 HEND5(i) = HEND5(i - 1);
	 break;
       case 2:
	 G = END5_1(i,1) - (temp * (END5_1(i,2)));
	 if(G < G2) {
	    SEND5(i) = END5_1(i,2);
	    HEND5(i) = END5_1(i,1);
	 } else {
	    SEND5(i) = SEND5(i - 1);
	    HEND5(i) = HEND5(i - 1);
	 }
	 break;
       case 3:
	 G = END5_2(i,1) - (temp * (END5_2(i,2)));
	 if(G < G2) {
	    SEND5(i) = END5_2(i,2);
	    HEND5(i) = END5_2(i,1);
	 } else {
	    SEND5(i) = SEND5(i - 1);
	    HEND5(i) = HEND5(i - 1);
	 }
	 break;
       case 4:
	 G = END5_3(i,1) - (temp * (END5_3(i,2)));
	 if(G < G2) {
	    SEND5(i) = END5_3(i,2);
	    HEND5(i) = END5_3(i,1);
	 } else {
	    SEND5(i) = SEND5(i - 1);
	    HEND5(i) = HEND5(i - 1);
	 }
	 break;
       case 5:
	 G = END5_4(i,1) - (temp * (END5_4(i,2)));
	 if(G < G2) {
	    SEND5(i) = END5_4(i,2);
	    HEND5(i) = END5_4(i,1);
	 } else {
	    SEND5(i) = SEND5(i - 1);
	    HEND5(i) = HEND5(i - 1);
	 }
	 break;
       default:
#ifdef DEBUG
	 printf ("WARNING: max5 returned character code %d ??\n", max);
#endif
	 break;
      }
   }
}

static void 
fillMatrix2(int maxLoop)
{
   int i, j;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double));
   for (j = 2; j <= len2; ++j)
     for (i = j - MIN_HRPN_LOOP - 1; i >= 1; --i) {
	if (isFinite(EnthalpyDPT(i, j))) {
	   SH[0] = -1.0;
	   SH[1] = _INFINITY;
	   maxTM2(i,j); /* calculate stack */
	   CBI(i, j, SH, 0,maxLoop); /* calculate Bulge and Internal loop and stack */
	   SH[0] = -1.0;
	   SH[1] = _INFINITY;
	   calc_hairpin(i, j, SH, 0);
	   if(isFinite(SH[1])) {
	      if(SH[0] < MinEntropyCutoff){ /* to not give dH any value if dS is unreasonable */
		 SH[0] = MinEntropy;
		 SH[1] = 0.0;
	      }
	      EntropyDPT(i,j) = SH[0];
	      EnthalpyDPT(i,j) = SH[1];
	   }
	}
     }
   free(SH);
}


static int 
symmetry_thermo(const unsigned char* seq)
{
   char s;
   char e;
   const unsigned char *seq_end=seq;
   int i = 0;
   int seq_len=length_unsig_char(seq);
   int mp = seq_len/2;
   if(seq_len%2==1) {
      return 0;
   }
   seq_end+=seq_len;
   seq_end--;
   while(i<mp) {
      i++;
      s=toupper(*seq);
      e=toupper(*seq_end);
      if ((s=='A' && e!='T')
	  || (s=='T' && e!='A')
	  || (e=='A' && s!='T')
	  || (e=='T' && s!='A')) {
	 return 0;
      }
      if ((s=='C' && e!='G')
	  || (s=='G' && e!='C')
	  || (e=='C' && s!='G')
	  || (e=='G' && s!='C')) {
	 return 0;
      }
      seq++;
      seq_end--;
   }
   return 1;
}

static void 
tracebacku(int* bp, int maxLoop) /* traceback for unimolecular structure */
{
   int i, j;
   i = j = 0;
   int ii, jj, k;
   struct tracer *top, *stack = NULL;
   double* SH1;
   double* SH2;
   double* EntropyEnthalpy;
   SH1 = (double*) safe_malloc(2 * sizeof(double));
   SH2 = (double*) safe_malloc(2 * sizeof(double));
   EntropyEnthalpy = (double*) safe_malloc(2 * sizeof(double));
   push(&stack,len1, 0, 1);
   while(stack) {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;
      if(top->mtrx==1) {
	 while (equal(SEND5(i), SEND5(i - 1)) && equal(HEND5(i), HEND5(i - 1))) /* if previous structure is the same as this one */
	   --i;
	 if (i == 0)
	   continue;
	 if (equal(SEND5(i), END5_1(i,2)) && equal(HEND5(i), END5_1(i,1))) {
	    for (k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k)
	      if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i)) &&
		  equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i))) {
		 push(&stack, k + 1, i,0);
		 break;
	      }
	    else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i)) &&
		     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i))) {
	       push(&stack, k + 1, i, 0);
	       push(&stack, k, 0, 1);
	       break;
	    }
	 }
	 else if (equal(SEND5(i), END5_2(i,2)) && equal(HEND5(i), END5_2(i,1))) {
	    for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
	      if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i)) &&
		  equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i))) {
		 push(&stack, k + 2, i, 0);
		 break;
	      }
	    else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) + Sd5(i, k + 2) + EntropyDPT(k + 2, i)) &&
		     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) + Hd5(i, k + 2) + EnthalpyDPT(k + 2, i))) {
	       push(&stack, k + 2, i, 0);
	       push(&stack, k, 0, 1);
	       break;
	    }
	 }
	 else if (equal(SEND5(i), END5_3(i,2)) && equal(HEND5(i), END5_3(i,1))) {
	    for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k)
	      if (equal(SEND5(i), atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1))
		  && equal(HEND5(i), atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
		 push(&stack, k + 1, i - 1, 0);
		 break;
	      }
	    else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) + Sd3(i - 1, k + 1) + EntropyDPT(k + 1, i - 1)) &&
		     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) + Hd3(i - 1, k + 1) + EnthalpyDPT(k + 1, i - 1))) {
	       push(&stack, k + 1, i - 1, 0); /* matrix 0  */
	       push(&stack, k, 0, 1); /* matrix 3 */
	       break;
	    }
	 }
	 else if(equal(SEND5(i), END5_4(i,2)) && equal(HEND5(i), END5_4(i,1))) {
	    for (k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k)
	      if (equal(SEND5(i), atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
		  equal(HEND5(i), atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1))) {
		 push(&stack, k + 2, i - 1, 0);
		 break;
	      }
	    else if (equal(SEND5(i), SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) + Ststack(i - 1, k + 2) + EntropyDPT(k + 2, i - 1)) &&
		     equal(HEND5(i), HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) + Htstack(i - 1, k + 2) + EnthalpyDPT(k + 2, i - 1)) ) {
	       push(&stack, k + 2, i - 1, 0);
	       push(&stack, k, 0, 1);
	       break;
	    }
	 }
      }
      else if(top->mtrx==0) {
	 bp[i - 1] = j;
	 bp[j - 1] = i;
	 SH1[0] = -1.0;
	 SH1[1] = _INFINITY;
	 calc_hairpin(i, j, SH1, 1); /* 1 means that we use this method in traceback */
	 SH2[0] = -1.0;
	 SH2[1] = _INFINITY;
	 CBI(i,j,SH2,2,maxLoop);
	 if (equal(EntropyDPT(i, j), Ss(i, j, 2) + EntropyDPT(i + 1, j - 1)) &&
	     equal(EnthalpyDPT(i, j), Hs(i, j, 2) + EnthalpyDPT(i + 1, j - 1))) {
	    push(&stack, i + 1, j - 1, 0);
	 }
	 else if (equal(EntropyDPT(i, j), SH1[0]) && equal(EnthalpyDPT(i,j), SH1[1]));
	 else if (equal(EntropyDPT(i, j), SH2[0]) && equal(EnthalpyDPT(i, j), SH2[1])) {
	    int d, done;
	    for (done = 0, d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop && !done; --d)
	      for (ii = i + 1; ii < j - d; ++ii) {
		 jj = d + ii;
		 EntropyEnthalpy[0] = -1.0;
		 EntropyEnthalpy[1] = _INFINITY;
		 calc_bulge_internal2(i, j, ii, jj,EntropyEnthalpy,1,maxLoop);
		 if (equal(EntropyDPT(i, j), EntropyEnthalpy[0] + EntropyDPT(ii, jj)) &&
		     equal(EnthalpyDPT(i, j), EntropyEnthalpy[1] + EnthalpyDPT(ii, jj))) {
		    push(&stack, ii, jj, 0);
		    ++done;
		    break;
		 }
	      }
	 } else {
	 }
      }
      free(top);
   }
   free(SH1);
   free(SH2);
   free(EntropyEnthalpy);
}


static void 
traceback(int i, int j, int* ps1, int* ps2, int maxLoop)
{
   int d, ii, jj, done;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double));
   ps1[i - 1] = j;
   ps2[j - 1] = i;
   while(1) {
      SH[0] = -1.0;
      SH[1] = _INFINITY;
      LSH(i,j,SH);
      if(equal(EntropyDPT(i,j),SH[0]) && equal(EnthalpyDPT(i,j),SH[1])) {
	 break;
      }
      done = 0;
      if (i > 1 && j > 1 && equal(EntropyDPT(i,j), Ss(i - 1, j - 1, 1) + EntropyDPT(i - 1, j - 1)) && equal(EnthalpyDPT(i,j), Hs(i - 1, j - 1, 1) + EnthalpyDPT(i - 1, j - 1))) {
	 i = i - 1;
	 j = j - 1;
	 ps1[i - 1] = j;
	 ps2[j - 1] = i;
	 done = 1;
      }
      for (d = 3; !done && d <= maxLoop + 2; ++d) {
	 ii = i - 1;
	 jj = -ii - d + (j + i);
	 if (jj < 1) {
	    ii -= abs(jj-1);
	    jj = 1;
	 }
	 for (; !done && ii > 0 && jj < j; --ii, ++jj) {
	    SH[0] = -1.0;
	    SH[1] = _INFINITY;
	    calc_bulge_internal(ii, jj, i, j, SH,1);
	    if (equal(EntropyDPT(i, j), SH[0]) && equal(EnthalpyDPT(i, j), SH[1])) {
	       i = ii;
	       j = jj;
	       ps1[i - 1] = j;
	       ps2[j - 1] = i;
	       done = 1;
	       break;
	    }
	 }
      }
   }
   free(SH);
}

inline void
drawDimer(int* ps1, int* ps2, double temp, double H, double S, int temponly, double t37, thal_results *o)
{
   int i, j, k, numSS1, numSS2, N;
   char* duplex[4];
   double G, t;
   t = G = 0;
   if (!isFinite(temp)){
      if(temponly==0) {
	 printf("No predicted secondary structures for given sequences\n");
      }
      o->temp = 0.0; /* lets use generalization here; this should rather be very negative value */
      std::cerr << "No predicted sec struc for given seq" << std::endl;
      return;
   } else {
      N=0;
      for(i=0;i<len1;i++){
	 if(ps1[i]>0) ++N;
      }
      for(i=0;i<len2;i++) {
	 if(ps2[i]>0) ++N;
      }
      N = (N/2) -1;
      t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
      if(temponly==0) {
	 G = (H) - (t37 * (S + (N * saltCorrection)));
	 S = S + (N * saltCorrection);
	 o->temp = (double) t;
	 /* maybe user does not need as precise as that */
	 /* printf("Thermodynamical values:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\tN = %d, SaltC=%f, RC=%f\n",
		len1, (double) S, (double) H, (double) G, (double) t, (int) N, saltCorrection, RC); */
	 printf("Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
		(double) S, (double) H, (double) G, (double) t);
      } else {
	 o->temp = (double) t;
	 return;
      }
   }

   duplex[0] = (char*) safe_malloc(len1 + len2 + 1);
   duplex[1] = (char*) safe_malloc(len1 + len2 + 1);
   duplex[2] = (char*) safe_malloc(len1 + len2 + 1);
   duplex[3] = (char*) safe_malloc(len1 + len2 + 1);
   duplex[0][0] = duplex[1][0] = duplex[2][0] = duplex[3][0] = 0;

   i = 0;
   numSS1 = 0;
   while (ps1[i++] == 0) ++numSS1;
   j = 0;
   numSS2 = 0;
   while (ps2[j++] == 0) ++numSS2;

   if (numSS1 >= numSS2){
      for (i = 0; i < numSS1; ++i) {
	 strcatc(duplex[0], oligo1[i]);
	 strcatc(duplex[1], ' ');
	 strcatc(duplex[2], ' ');
      }
      for (j = 0; j < numSS1 - numSS2; ++j) strcatc(duplex[3], ' ');
      for (j = 0; j < numSS2; ++j) strcatc(duplex[3], oligo2[j]);
   } else {
      for (j = 0; j < numSS2; ++j) {
	 strcatc(duplex[3], oligo2[j]);
	 strcatc(duplex[1], ' ');
	 strcatc(duplex[2], ' ');
      }
      for (i = 0; i < numSS2 - numSS1; ++i)
	strcatc(duplex[0], ' ');
      for (i = 0; i < numSS1; ++i)
	strcatc(duplex[0], oligo1[i]);
   }
   i = numSS1 + 1;
   j = numSS2 + 1;

   while (i <= len1) {
      while (i <= len1 && ps1[i - 1] != 0 && j <= len2 && ps2[j - 1] != 0) {
	 strcatc(duplex[0], ' ');
	 strcatc(duplex[1], oligo1[i - 1]);
	 strcatc(duplex[2], oligo2[j - 1]);
	 strcatc(duplex[3], ' ');
	 ++i;
	 ++j;
      }
      numSS1 = 0;
      while (i <= len1 && ps1[i - 1] == 0) {
	 strcatc(duplex[0], oligo1[i - 1]);
	 strcatc(duplex[1], ' ');
	 ++numSS1;
	 ++i;
      }
      numSS2 = 0;
      while (j <= len2 && ps2[j - 1] == 0) {
	 strcatc(duplex[2], ' ');
	 strcatc(duplex[3], oligo2[j - 1]);
	 ++numSS2;
	 ++j;
      }
      if (numSS1 < numSS2)
	for (k = 0; k < numSS2 - numSS1; ++k) {
	   strcatc(duplex[0], '-');
	   strcatc(duplex[1], ' ');
	}
      else if (numSS1 > numSS2)
	for (k = 0; k < numSS1 - numSS2; ++k) {
	   strcatc(duplex[2], ' ');
	   strcatc(duplex[3], '-');
	}
   }
   printf("SEQ\t");
   printf("%s\n", duplex[0]);
   printf("SEQ\t");
   printf("%s\n", duplex[1]);
   printf("STR\t");
   printf("%s\n", duplex[2]);
   printf("STR\t");
   printf("%s\n", duplex[3]);

   free(duplex[0]);
   free(duplex[1]);
   free(duplex[2]);
   free(duplex[3]);

   return;
}

inline bool
drawHairpin(int* bp, double mh, double ms, int temponly, double temp, thal_results *o)
{
   /* Plain text */
   int i, N;
   N = 0;
   double mg, t;
   if (!isFinite(ms) || !isFinite(mh)) {
      if(temponly == 0) {
	 printf("0\tdS = %g\tdH = %g\tinf\tinf\n", (double) ms,(double) mh);
#ifdef DEBUG
	 fputs("No temperature could be calculated\n",stderr);
#endif
      } else {
	 o->temp = 0.0; /* lets use generalization here */
	 std::cerr << "No predicted sec struc for given seq" << std::endl;
	 return false;
      }
   } else {
      if(temponly == 0) {
	 for (i = 1; i < len1; ++i) {
	    if(bp[i-1] > 0) N++;
	 }
      } else {
	 for (i = 1; i < len1; ++i) {
	    if(bp[i-1] > 0) N++;
	 }
      }
      t = (mh / (ms + (((N/2)-1) * saltCorrection))) - ABSOLUTE_ZERO;
      if(temponly == 0) {
	 mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));
	 ms = ms + (((N/2)-1) * saltCorrection);
	 o->temp = (double) t;
	 printf("Calculated thermodynamical parameters for dimer:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
		len1, (double) ms, (double) mh, (double) mg, (double) t);
      } else {
	 o->temp = (double) t;
	 return true;
      }
   }
   /* plain-text output */
   char* asciiRow;
   asciiRow = (char*) safe_malloc(len1);
   for(i = 0; i < len1; ++i) asciiRow[i] = '0';
   for(i = 1; i < len1+1; ++i) {
      if(bp[i-1] == 0) {
	 asciiRow[(i-1)] = '-';
      } else {
	 if(bp[i-1] > (i-1)) {
	    asciiRow[(bp[i-1]-1)]='\\';
	 } else  {
	    asciiRow[(bp[i-1]-1)]='/';
	 }
      }
   }
   printf("SEQ\t");
   for(i = 0; i < len1; ++i) printf("%c",asciiRow[i]);
   printf("\nSTR\t%s\n", oligo1);
   free(asciiRow);
   return true;
}


/* Read the thermodynamic values (parameters) from the parameter files
   in the directory specified by 'path'.  Return 0 on success and -1
   on error. The thermodynamic values are stored in multiple static
   variables. */
int 
get_thermodynamic_values(const char* path)
{
  parampath = (char*) safe_malloc((strlen(path) + 1) * sizeof(char));
  strcpy(parampath, path);

  getStack(stackEntropies, stackEnthalpies);
  /* verifyStackTable(stackEntropies, "entropy");
     verifyStackTable(stackEnthalpies, "enthalpy"); */ /* this is for code debugging */
  getStackint2(stackint2Entropies, stackint2Enthalpies);
  getDangle(dangleEntropies3, dangleEnthalpies3, dangleEntropies5, dangleEnthalpies5);
  getLoop(hairpinLoopEntropies, interiorLoopEntropies, bulgeLoopEntropies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies);
  getTstack(tstackEntropies, tstackEnthalpies);
  getTstack2(tstack2Entropies, tstack2Enthalpies);
  getTriloop(&triloopEntropies, &triloopEnthalpies, &numTriloops);
  getTetraloop(&tetraloopEntropies, &tetraloopEnthalpies, &numTetraloops);
  /* getting the AT-penalties */
  tableStartATS(AT_S, atpS);
  tableStartATH(AT_H, atpH);

  return 0;
}

void 
destroy_thal_structures()
{
  free(parampath);
  free(triloopEntropies);
  free(triloopEnthalpies);
  free(tetraloopEntropies);
  free(tetraloopEnthalpies);
}

/* central method: execute all sub-methods for calculating secondary
   structure for dimer or for monomer */
inline bool
thal(const unsigned char *oligo_f, const unsigned char *oligo_r, const thal_args *a, thal_results *o) {
  if (NULL == o) {
    std::cerr << "No result object!" << std::endl;
    return false;
  }

   double* SH;
   int i, j;
   int len_f, len_r;
   int k;
   int *bp;
   double mh, ms;
   double G1, bestG;

   send5 = hend5 = NULL;
   enthalpyDPT = entropyDPT = NULL;
   numSeq1 = numSeq2 = NULL;
   oligo1 = oligo2 = NULL;
   o->temp = THAL_ERROR_SCORE;
   errno = 0; 

   if (NULL == oligo_f) {
     std::cerr << "NULL first sequence" << std::endl;
     return false;
   }
   if (NULL == oligo_r) {
     std::cerr << "NULL second sequence" << std::endl;
     return false;
   }
   len_f = length_unsig_char(oligo_f);
   len_r = length_unsig_char(oligo_r);
   if ((len_f > THAL_MAX_ALIGN) && (len_r > THAL_MAX_ALIGN)) {
     std::cerr << "Both sequences longer than " << THAL_MAX_ALIGN << " for thermodynamic alignment" << std::endl;
     return false;
   }
   if (len_f > THAL_MAX_SEQ) {
     std::cerr << "Target sequence length > maximum allowed!" << std::endl;
     return false;
   }
   if (len_r > THAL_MAX_SEQ) {
     std::cerr << "Target sequence length > maximum allowed!" << std::endl;
     return false;
   }

   if (NULL == a) {
     std::cerr << "NULL 'in' pointer" << std::endl;
     return false;
   }
   if (a->type != thal_any && a->type != thal_end1 && a->type != thal_end2 && a->type != thal_hairpin) {
     std::cerr << "Illegal type" << std::endl;
     return false;
   }
   o->align_end_1 = -1;
   o->align_end_2 = -1;
   if ('\0' == oligo_f[0]) {
     std::cerr << "Empty first sequence" << std::endl;
     o->temp = 0.0;
     return false;
   }
   if ('\0' == oligo_r[0]) {
     std::cerr << "Empty second sequence" << std::endl;
     o->temp = 0.0;
     return false;
   }
   if (0 == len_f) {
      o->temp = 0.0;
      return false;
   }
   if (0 == len_r) {
      o->temp = 0.0;
      return false;
   }
   if(a->type!=3) {
      oligo1 = (unsigned char*) safe_malloc((len_f + 1) * sizeof(unsigned char));
      oligo2 = (unsigned char*) safe_malloc((len_r + 1) * sizeof(unsigned char));
      strcpy((char*)oligo1,(const char*)oligo_f);
      strcpy((char*)oligo2,(const char*)oligo_r);
   } else  {
      oligo1 = (unsigned char*) safe_malloc((len_r + 1) * sizeof(unsigned char));
      oligo2 = (unsigned char*) safe_malloc((len_f + 1) * sizeof(unsigned char));
      strcpy((char*)oligo1,(const char*)oligo_r);
      strcpy((char*)oligo2,(const char*)oligo_f);
   }
   /*** INIT values for unimolecular and bimolecular structures ***/
   unsigned char *oligo2_rev = NULL;
   if (a->type==4) { /* unimolecular folding */
      len2 = length_unsig_char(oligo2);
      len3 = len2 -1;
      dplx_init_H = 0.0;
      dplx_init_S = -0.00000000001;
      RC=0;
   } else if(a->type!=4) {
      /* hybridization of two oligos */
      dplx_init_H = 200;
      dplx_init_S = -5.7;
      if(symmetry_thermo(oligo1) && symmetry_thermo(oligo2)) {
	 RC = R  * log(a->dna_conc/1000000000.0);
      } else {
	 RC = R  * log(a->dna_conc/4000000000.0);
      }
      if(a->type!=3) {
	 oligo2_rev = (unsigned char*) safe_malloc((length_unsig_char(oligo_r) + 1) * sizeof(unsigned char));
	 strcpy((char*)oligo2_rev,(const char*)oligo_r);
      } else {
	 oligo2_rev = (unsigned char*) safe_malloc((length_unsig_char(oligo_f) + 1) * sizeof(unsigned char));
	 strcpy((char*)oligo2_rev,(const char*)oligo_f);
      }
      reverse(oligo2_rev); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
      free(oligo2);
      oligo2=NULL;
      oligo2=&oligo2_rev[0];
   } else {
     std::cerr << "Wrong alignment type!" << std::endl;
     o->temp = THAL_ERROR_SCORE;
     errno=0;
     return false;
   }
   len1 = length_unsig_char(oligo1);
   len2 = length_unsig_char(oligo2);
   /* convert nucleotides to numbers */
   numSeq1 = (unsigned char*) safe_realloc(numSeq1, len1 + 2);
   numSeq2 = (unsigned char*) safe_realloc(numSeq2, len2 + 2);

   /*** Calc part of the salt correction ***/
   saltCorrection=saltCorrectS(a->mv,a->dv,a->dntp); /* salt correction for entropy, must be multiplied with N, which is
						   the total number of phosphates in the duplex divided by 2; 8bp dplx N=7 */

   if(a->type == 4){ /* monomer */
      /* terminal basepairs */
      send5 = (double*) safe_realloc(send5, (len1 + 1) * sizeof(double));
      hend5 = (double*) safe_realloc(hend5, (len1 + 1) * sizeof(double));
   }
   for(i = 0; i < len1; i++) oligo1[i] = toupper(oligo1[i]);
   for(i = 0; i < len2; i++) oligo2[i] = toupper(oligo2[i]);
   for(i = 1; i <= len1; ++i) numSeq1[i] = str2int(oligo1[i - 1]);
   for(i = 1; i <= len2; ++i) numSeq2[i] = str2int(oligo2[i - 1]);
   numSeq1[0] = numSeq1[len1 + 1] = numSeq2[0] = numSeq2[len2 + 1] = 4; /* mark as N-s */
   if (a->type==4) { /* calculate structure of monomer */
      enthalpyDPT = safe_recalloc(enthalpyDPT, len1, len2);
      entropyDPT = safe_recalloc(entropyDPT, len1, len2);
      initMatrix2();
      fillMatrix2(a->maxLoop);
      calc_terminal_bp(a->temp);
      mh = HEND5(len1);
      ms = SEND5(len1);
      o->align_end_1 = (int) mh;
      o->align_end_2 = (int) ms;
      bp = (int*) safe_calloc(len1, sizeof(int));
      for (k = 0; k < len1; ++k) bp[k] = 0;
      if(isFinite(mh)) {
	 tracebacku(bp, a->maxLoop);
	 /* traceback for unimolecular structure */
	 if (!drawHairpin(bp, mh, ms, a->temponly,a->temp, o)) return false; /* if temponly=1 then return after printing basic therm data */
      } else if(a->temponly==0) {
	    fputs("No secondary structure could be calculated\n",stderr);
      }

      if (o->temp==-_INFINITY) o->temp=0.0;
      free(bp);
      free(enthalpyDPT);
      free(entropyDPT);
      free(numSeq1);
      free(numSeq2);
      free(send5);
      free(hend5);
      free(oligo1);
      free(oligo2);
   } else if(a->type!=4) { /* Hybridization of two moleculs */
      len3 = len2;
      enthalpyDPT = safe_recalloc(enthalpyDPT, len1, len2); /* dyn. programming table for dS and dH */
      entropyDPT = safe_recalloc(entropyDPT, len1, len2); /* enthalpyDPT is 3D array represented as 1D array */
      initMatrix();
      fillMatrix(a->maxLoop);
      SH = (double*) safe_malloc(2 * sizeof(double));
      /* calculate terminal basepairs */
      bestI = bestJ = 0; 
      G1 = bestG = _INFINITY;
      if(a->type==1)
	for (i = 1; i <= len1; i++) {
	   for (j = 1; j <= len2; j++) {
	      RSH(i, j, SH);
	      SH[0] = SH[0]+SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0 */
	      SH[1] = SH[1]+SMALL_NON_ZERO;
	      G1 = (EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) - TEMP_KELVIN*(EntropyDPT(i, j) + SH[0] + dplx_init_S);  
	      if(G1<bestG){
		 bestG = G1;
		 bestI = i;
		 bestJ = j;
	      }
	   }
	}
      int *ps1, *ps2;
      ps1 = (int*) safe_calloc(len1, sizeof(int));
      ps2 = (int*) safe_calloc(len2, sizeof(int));
      for (i = 0; i < len1; ++i)
	ps1[i] = 0;
      for (j = 0; j < len2; ++j)
	ps2[j] = 0;
      if(a->type == 2 || a->type == 3)	{
	 /* THAL_END1 */
	 bestI = bestJ = 0;
	 bestI = len1;
	 i = len1;
	 G1 = bestG = _INFINITY;
	 for (j = 1; j <= len2; ++j) {
	    RSH(i, j, SH);
	    SH[0] = SH[0]+SMALL_NON_ZERO; /* this adding is done for compiler, optimization -O2 vs -O0,
					   that compiler could understand that SH is changed in this cycle */
	    SH[1] = SH[1]+SMALL_NON_ZERO;
	    G1 = (EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) - TEMP_KELVIN*(EntropyDPT(i, j) + SH[0] + dplx_init_S);  
		if(G1<bestG){
		   bestG = G1;
      	           bestJ = j;
	    }
	 }
      }
      if (!isFinite(bestG)) bestI = bestJ = 1;
      double dH, dS;
      RSH(bestI, bestJ, SH);
      dH = EnthalpyDPT(bestI, bestJ)+ SH[1] + dplx_init_H;
      dS = (EntropyDPT(bestI, bestJ) + SH[0] + dplx_init_S);
      /* tracebacking */
      for (i = 0; i < len1; ++i)
	ps1[i] = 0;
      for (j = 0; j < len2; ++j)
	ps2[j] = 0;
      if(isFinite(EnthalpyDPT(bestI, bestJ))){
	 traceback(bestI, bestJ, ps1, ps2, a->maxLoop);
	 drawDimer(ps1, ps2, SHleft, dH, dS, a->temponly,a->temp, o);
	 o->align_end_1=bestI;
	 o->align_end_2=bestJ;
      } else  {
	 o->temp = 0.0;
	 /* fputs("No secondary structure could be calculated\n",stderr); */
      }
      free(ps1);
      free(ps2);
      free(SH);
      free(oligo2_rev);
      free(enthalpyDPT);
      free(entropyDPT);
      free(numSeq1);
      free(numSeq2);
      free(oligo1);
   }
   return true;
}


/* Set default args */
void 
set_thal_default_args(thal_args *a)
{
   a->debug = 0;
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.8; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = TEMP_KELVIN; /* Kelvin */
   a->temponly = 1; /* return only melting temperature of predicted structure */
   a->dimer = 1; /* by default dimer structure is calculated */
}

/* Set default args for oligo */
void
set_thal_oligo_default_args(thal_args *a)    
{
   a->debug = 0;
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.0; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = TEMP_KELVIN; /* Kelvin */
   a->temponly = 1; /* return only melting temperature of predicted structure */
   a->dimer = 1; /* by default dimer structure is calculated */
}

}

#endif
