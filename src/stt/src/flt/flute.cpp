// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#include "stt/flute.h"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "utl/decode.h"

// Use flute LUT file reader.
#define LUT_FILE 1
// Init LUTs from base64 encoded string variables.
#define LUT_VAR 2
// Init LUTs from base64 encoded string variables
// and check against LUTs from file reader.
#define LUT_VAR_CHECK 3

// Set this to LUT_FILE, LUT_VAR, or LUT_VAR_CHECK.
// #define LUT_SOURCE LUT_FILE
// #define LUT_SOURCE LUT_VAR_CHECK
#define LUT_SOURCE LUT_VAR

#if FLUTE_D <= 7
#define MGROUP 5040 / 4  // Max. # of groups, 7! = 5040
#define MPOWV 15         // Max. # of POWVs per group
#elif FLUTE_D == 8
#define MGROUP 40320 / 4  // Max. # of groups, 8! = 40320
#define MPOWV 33          // Max. # of POWVs per group
#elif FLUTE_D == 9
#define MGROUP 362880 / 4  // Max. # of groups, 9! = 362880
#define MPOWV 79           // Max. # of POWVs per group
#endif

namespace stt {

namespace flt {

struct Flute::csoln
{
  unsigned char parent;
  unsigned char seg[11];  // Add: 0..i, Sub: j..10; seg[i+1]=seg[j-1]=0
  unsigned char rowcol[FLUTE_D - 2];  // row = rowcol[]/16, col = rowcol[]%16,
  unsigned char neighbor[2 * FLUTE_D - 2];
};

struct point
{
  int x, y;
  int o;
};

template <class T>
inline T ADIFF(T x, T y)
{
  if (x > y) {
    return (x - y);
  }
  return (y - x);
}

////////////////////////////////////////////////////////////////

#if LUT_SOURCE == LUT_FILE || LUT_SOURCE == LUT_VAR_CHECK
static void readLUTfiles(LUT_TYPE LUT, NUMSOLN_TYPE numsoln)
{
  unsigned char charnum[256], line[32], *linep, c;
  FILE *fpwv, *fprt;
  int d, i, j, k, kk, ns, nn;

  for (i = 0; i <= 255; i++) {
    if ('0' <= i && i <= '9') {
      charnum[i] = i - '0';
    } else if (i >= 'A') {
      charnum[i] = i - 'A' + 10;
    } else {  // if (i=='$' || i=='\n' || ... )
      charnum[i] = 0;
    }
  }

  fpwv = fopen(FLUTE_POWVFILE, "r");
  if (fpwv == NULL) {
    printf("Error in opening %s\n", FLUTE_POWVFILE);
    exit(1);
  }

#if FLUTE_ROUTING == 1
  fprt = fopen(FLUTE_POSTFILE, "r");
  if (fprt == NULL) {
    printf("Error in opening %s\n", FLUTE_POSTFILE);
    exit(1);
  }
#endif

  for (d = 4; d <= FLUTE_D; d++) {
    fscanf(fpwv, "d=%d", &d);
    fgetc(fpwv);  // '/n'
#if FLUTE_ROUTING == 1
    fscanf(fprt, "d=%d", &d);
    fgetc(fprt);  // '/n'
#endif
    for (k = 0; k < numgrp[d]; k++) {
      ns = (int) charnum[fgetc(fpwv)];

      if (ns == 0) {  // same as some previous group
        fscanf(fpwv, "%d", &kk);
        fgetc(fpwv);  // '/n'
        numsoln[d][k] = numsoln[d][kk];
        (*LUT)[d][k] = (*LUT)[d][kk];
      } else {
        fgetc(fpwv);  // '\n'
        numsoln[d][k] = ns;
        auto p = std::make_shared<struct csoln[]>(ns);
        for (i = 1; i <= ns; i++) {
          linep = (unsigned char*) fgets((char*) line, 32, fpwv);
          p->parent = charnum[*(linep++)];
          j = 0;
          while ((p->seg[j++] = charnum[*(linep++)]) != 0)
            ;
          j = 10;
          while ((p->seg[j--] = charnum[*(linep++)]) != 0)
            ;
#if FLUTE_ROUTING == 1
          nn = 2 * d - 2;
          fread(line, 1, d - 2, fprt);
          linep = line;
          for (j = d; j < nn; j++) {
            c = charnum[*(linep++)];
            p->rowcol[j - d] = c;
          }
          fread(line, 1, nn / 2 + 1, fprt);
          linep = line;  // last char \n
          for (j = 0; j < nn;) {
            c = *(linep++);
            p->neighbor[j++] = c / 16;
            p->neighbor[j++] = c % 16;
          }
#endif
          p++;
        }
        (*LUT)[d][k] = std::move(p);
      }
    }
  }
  fclose(fpwv);
#if FLUTE_ROUTING == 1
  fclose(fprt);
#endif
}
#endif

////////////////////////////////////////////////////////////////

#if LUT_SOURCE == LUT_VAR_CHECK
static void checkLUT(LUT_TYPE LUT1,
                     NUMSOLN_TYPE numsoln1,
                     LUT_TYPE LUT2,
                     NUMSOLN_TYPE numsoln2);
#endif

extern const char* post9[];
extern const char* powv9[];

void Flute::readLUT()
{
  makeLUT(LUT_, numsoln_);

#if LUT_SOURCE == LUT_FILE
  readLUTfiles(LUT, numsoln);
  lut_valid_d = FLUTE_D;

#elif LUT_SOURCE == LUT_VAR
  // Only init to d=8 on startup because d=9 is big and slow.
  initLUT(lut_initial_d, LUT_, numsoln_);

#elif LUT_SOURCE == LUT_VAR_CHECK
  readLUTfiles(LUT, numsoln);
  // Temporaries to compare to file results.
  LUT_TYPE LUT_;
  NUMSOLN_TYPE numsoln_;
  makeLUT(LUT_, numsoln_);
  initLUT(FLUTE_D, LUT_, numsoln_);
  checkLUT(LUT, numsoln, LUT_, numsoln_);
  deleteLUT(LUT_, numsoln_);
#endif
}

void Flute::makeLUT(LUT_TYPE& LUT, NUMSOLN_TYPE& numsoln)
{
  LUT = new boost::multi_array<std::shared_ptr<struct csoln[]>, 2>(
      boost::extents[FLUTE_D + 1][MGROUP]);
  numsoln = new int*[FLUTE_D + 1];
  for (int d = 4; d <= FLUTE_D; d++) {
    numsoln[d] = new int[MGROUP];
  }
}

void Flute::deleteLUT()
{
  deleteLUT(LUT_, numsoln_);
}

void Flute::deleteLUT(LUT_TYPE& LUT, NUMSOLN_TYPE& numsoln)
{
  if (LUT) {
    delete LUT;
    for (int d = 4; d <= FLUTE_D; d++) {
      delete[] numsoln[d];
    }
    delete[] numsoln;
  }
}

static unsigned char charNum(unsigned char c)
{
  if (isdigit(c)) {
    return c - '0';
  }
  if (c >= 'A') {
    return c - 'A' + 10;
  }
  return 0;
}

inline const char* readDecimalInt(const char* s, int& value)
{
  value = 0;
  bool negative = (*s == '-');
  if (negative || *s == '+') {
    ++s;
  }
  constexpr int zero_code = int('0');
  while (*s >= '0' && *s <= '9') {
    value = 10 * value + (int(*s) - zero_code);
    ++s;
  }
  if (negative) {
    value = -value;
  }
  return s;
}

// Init LUTs from base64 encoded string variables.
void Flute::initLUT(int to_d, LUT_TYPE LUT, NUMSOLN_TYPE numsoln)
{
  std::string pwv_string = utl::base64_decode(powv9);
  const char* pwv = pwv_string.c_str();

#if FLUTE_ROUTING == 1
  std::string prt_string = utl::base64_decode(post9);
  const char* prt = prt_string.c_str();
#endif

  for (int d = 4; d <= to_d; d++) {
    if (pwv[0] == 'd' && pwv[1] == '=') {
      pwv = readDecimalInt(pwv + 2, d);
    }
    ++pwv;
#if FLUTE_ROUTING == 1
    if (prt[0] == 'd' && prt[1] == '=') {
      prt = readDecimalInt(prt + 2, d);
    }
    ++prt;
#endif
    for (int k = 0; k < numgrp[d]; k++) {
      int ns = charNum(*pwv++);
      if (ns == 0) {  // same as some previous group
        int kk;
        pwv = readDecimalInt(pwv, kk) + 1;
        numsoln[d][k] = numsoln[d][kk];
        (*LUT)[d][k] = (*LUT)[d][kk];
      } else {
        pwv++;  // '\n'
        numsoln[d][k] = ns;
        struct csoln* p = new struct csoln[ns];
        (*LUT)[d][k] = std::shared_ptr<struct csoln[]>(p);
        for (int i = 1; i <= ns; i++) {
          p->parent = charNum(*pwv++);

          int j = 0;
          unsigned char ch, seg;
          do {
            ch = *pwv++;
            seg = charNum(ch);
            p->seg[j++] = seg;
          } while (seg != 0);

          j = 10;
          if (ch == '\n') {
            p->seg[j] = 0;
          } else {
            do {
              ch = *pwv++;
              seg = charNum(ch);
              p->seg[j--] = seg;
            } while (seg != 0);
          }

#if FLUTE_ROUTING == 1
          int nn = 2 * d - 2;
          for (int j = d; j < nn; j++) {
            p->rowcol[j - d] = charNum(*prt++);
          }

          for (int j = 0; j < nn;) {
            unsigned char c = *prt++;
            p->neighbor[j++] = c / 16;
            p->neighbor[j++] = c % 16;
          }
          prt++;  // \n
#endif
          p++;
        }
      }
    }
  }
  lut_valid_d_ = to_d;
}

void Flute::ensureLUT(int d)
{
  if (LUT_ == nullptr) {
    readLUT();
  }
  if (d > lut_valid_d_ && d <= FLUTE_D) {
    initLUT(FLUTE_D, LUT_, numsoln_);
  }
}

#if LUT_SOURCE == LUT_VAR_CHECK
static void checkLUT(LUT_TYPE LUT1,
                     NUMSOLN_TYPE numsoln1,
                     LUT_TYPE LUT2,
                     NUMSOLN_TYPE numsoln2)
{
  for (int d = 4; d <= FLUTE_D; d++) {
    for (int k = 0; k < numgrp[d]; k++) {
      int ns1 = numsoln1[d][k];
      int ns2 = numsoln2[d][k];
      if (ns1 != ns2) {
        printf("numsoln[%d][%d] mismatch\n", d, k);
      }
      struct csoln* soln1 = LUT1[d][k].get();
      struct csoln* soln2 = LUT2[d][k].get();
      if (soln1->parent != soln2->parent) {
        printf("LUT[%d][%d]->parent mismatch\n", d, k);
      }
      for (int j = 0; soln1->seg[j] != 0; j++) {
        if (soln1->seg[j] != soln2->seg[j]) {
          printf("LUT[%d][%d]->seg[%d] mismatch\n", d, k, j);
        }
      }
      for (int j = 10; soln1->seg[j] != 0; j--) {
        if (soln1->seg[j] != soln2->seg[j]) {
          printf("LUT[%d][%d]->seg[%d] mismatch\n", d, k, j);
        }
      }
      int nn = 2 * d - 2;
      for (int j = d; j < nn; j++) {
        if (soln1->rowcol[j - d] != soln2->rowcol[j - d]) {
          printf("LUT[%d][%d]->rowcol[%d] mismatch\n", d, k, j);
        }
      }
      for (int j = 0; j < nn; j++) {
        if (soln1->neighbor[j] != soln2->neighbor[j]) {
          printf("LUT[%d][%d]->neighbor[%d] mismatch\n", d, k, j);
        }
      }
    }
  }
}
#endif

////////////////////////////////////////////////////////////////

int Flute::flute_wl(int d,
                    const std::vector<int>& x,
                    const std::vector<int>& y,
                    int acc)
{
  int minval, l, xu, xl, yu, yl;
  std::vector<int> xs, ys;
  int i, j, minidx;
  std::vector<int> s;
  struct point **ptp, *tmpp;
  struct point* pt;

  /* allocate the dynamic pieces on the heap rather than the stack */
  xs.resize(d);
  ys.resize(d);
  s.resize(d);
  pt = (struct point*) malloc(sizeof(struct point) * (d + 1));
  ptp = (struct point**) malloc(sizeof(struct point*) * (d + 1));

  if (d == 2) {
    l = ADIFF(x[0], x[1]) + ADIFF(y[0], y[1]);
  } else if (d == 3) {
    if (x[0] > x[1]) {
      xu = std::max(x[0], x[2]);
      xl = std::min(x[1], x[2]);
    } else {
      xu = std::max(x[1], x[2]);
      xl = std::min(x[0], x[2]);
    }
    if (y[0] > y[1]) {
      yu = std::max(y[0], y[2]);
      yl = std::min(y[1], y[2]);
    } else {
      yu = std::max(y[1], y[2]);
      yl = std::min(y[0], y[2]);
    }
    l = (xu - xl) + (yu - yl);
  } else {
    ensureLUT(d);

    for (i = 0; i < d; i++) {
      pt[i].x = x[i];
      pt[i].y = y[i];
      ptp[i] = &pt[i];
    }

    // sort x
    for (i = 0; i < d - 1; i++) {
      minval = ptp[i]->x;
      minidx = i;
      for (j = i + 1; j < d; j++) {
        if (minval > ptp[j]->x) {
          minval = ptp[j]->x;
          minidx = j;
        }
      }
      tmpp = ptp[i];
      ptp[i] = ptp[minidx];
      ptp[minidx] = tmpp;
    }

#if FLUTE_REMOVE_DUPLICATE_PIN == 1
    ptp[d] = &pt[d];
    ptp[d]->x = ptp[d]->y = -999999;
    j = 0;
    for (i = 0; i < d; i++) {
      for (k = i + 1; ptp[k]->x == ptp[i]->x; k++) {
        if (ptp[k]->y == ptp[i]->y) {  // pins k and i are the same
          break;
        }
      }
      if (ptp[k]->x != ptp[i]->x) {
        ptp[j++] = ptp[i];
      }
    }
    d = j;
#endif

    for (i = 0; i < d; i++) {
      xs[i] = ptp[i]->x;
      ptp[i]->o = i;
    }

    // sort y to find s[]
    for (i = 0; i < d - 1; i++) {
      minval = ptp[i]->y;
      minidx = i;
      for (j = i + 1; j < d; j++) {
        if (minval > ptp[j]->y) {
          minval = ptp[j]->y;
          minidx = j;
        }
      }
      ys[i] = ptp[minidx]->y;
      s[i] = ptp[minidx]->o;
      ptp[minidx] = ptp[i];
    }
    ys[d - 1] = ptp[d - 1]->y;
    s[d - 1] = ptp[d - 1]->o;

    l = flutes_wl(d, xs, ys, s, acc);
  }
  free(pt);
  free(ptp);

  return l;
}

// xs[] and ys[] are coords in x and y in sorted order
// s[] is a list of nodes in increasing y direction
//   if nodes are indexed in the order of increasing x coord
//   i.e., s[i] = s_i as defined in paper
// The points are (xs[s[i]], ys[i]) for i=0..d-1
//             or (xs[i], ys[si[i]]) for i=0..d-1

int Flute::flutes_wl_RDP(int d,
                         std::vector<int> xs,
                         std::vector<int> ys,
                         std::vector<int> s,
                         int acc)
{
  int i, j, ss;

  ensureLUT(d);

  for (i = 0; i < d - 1; i++) {
    if (xs[s[i]] == xs[s[i + 1]] && ys[i] == ys[i + 1]) {
      if (s[i] < s[i + 1]) {
        ss = s[i + 1];
      } else {
        ss = s[i];
        s[i] = s[i + 1];
      }
      for (j = i + 2; j < d; j++) {
        ys[j - 1] = ys[j];
        s[j - 1] = s[j];
      }
      for (j = ss + 1; j < d; j++) {
        xs[j - 1] = xs[j];
      }
      for (j = 0; j <= d - 2; j++) {
        if (s[j] > ss) {
          s[j]--;
        }
      }
      i--;
      d--;
    }
  }
  return flutes_wl_ALLD(d, xs, ys, s, acc);
}

// For low-degree, i.e., 2 <= d <= FLUTE_D
int Flute::flutes_wl_LD(int d,
                        const std::vector<int>& xs,
                        const std::vector<int>& ys,
                        const std::vector<int>& s)
{
  int k, pi, i, j;
  struct csoln* rlist;
  int dd[2 * FLUTE_D - 2];  // 0..FLUTE_D-2 for v, FLUTE_D-1..2*D-3 for h
  int minl, sum, l[MPOWV + 1];

  if (d <= 3) {
    minl = xs[d - 1] - xs[0] + ys[d - 1] - ys[0];
  } else {
    ensureLUT(d);

    k = 0;
    if (s[0] < s[2]) {
      k++;
    }
    if (s[1] < s[2]) {
      k++;
    }

    for (i = 3; i <= d - 1; i++) {  // p0=0 always, skip i=1 for symmetry
      pi = s[i];
      for (j = d - 1; j > i; j--) {
        if (s[j] < s[i]) {
          pi--;
        }
      }
      k = pi + (i + 1) * k;
    }

    if (k < numgrp[d]) {  // no horizontal flip
      for (i = 1; i <= d - 3; i++) {
        dd[i] = ys[i + 1] - ys[i];
        dd[d - 1 + i] = xs[i + 1] - xs[i];
      }
    } else {
      k = 2 * numgrp[d] - 1 - k;
      for (i = 1; i <= d - 3; i++) {
        dd[i] = ys[i + 1] - ys[i];
        dd[d - 1 + i] = xs[d - 1 - i] - xs[d - 2 - i];
      }
    }

    minl = l[0] = xs[d - 1] - xs[0] + ys[d - 1] - ys[0];
    rlist = (*LUT_)[d][k].get();
    for (i = 0; rlist->seg[i] > 0; i++) {
      minl += dd[rlist->seg[i]];
    }

    l[1] = minl;
    j = 2;
    while (j <= numsoln_[d][k]) {
      rlist++;
      sum = l[rlist->parent];
      for (i = 0; rlist->seg[i] > 0; i++) {
        sum += dd[rlist->seg[i]];
      }
      for (i = 10; rlist->seg[i] > 0; i--) {
        sum -= dd[rlist->seg[i]];
      }
      minl = std::min(minl, sum);
      l[j++] = sum;
    }
  }

  return minl;
}

// For medium-degree, i.e., FLUTE_D+1 <= d
int Flute::flutes_wl_MD(int d,
                        const std::vector<int>& xs,
                        const std::vector<int>& ys,
                        const std::vector<int>& s,
                        int acc)
{
  float pnlty, dx, dy;
  float *score, *penalty;
  int xydiff;
  int ll, minl;
  int extral = 0;
  std::vector<int> x1, x2, y1, y2;
  std::vector<int> distx, disty;
  int i, r, p, maxbp, nbp, bp, ub, lb, n1, n2, newacc;
  int ms, mins, maxs, minsi, maxsi, degree;
  int return_val;
  std::vector<int> si, s1, s2;

  degree = d + 1;
  score = (float*) malloc(sizeof(float) * (2 * degree));
  penalty = (float*) malloc(sizeof(float) * (degree));

  x1.resize(degree);
  x2.resize(degree);
  y1.resize(degree);
  y2.resize(degree);
  distx.resize(degree);
  disty.resize(degree);
  si.resize(degree);
  s1.resize(degree);
  s2.resize(degree);

  ensureLUT(d);

  if (s[0] < s[d - 1]) {
    ms = std::max(s[0], s[1]);
    for (i = 2; i <= ms; i++) {
      ms = std::max(ms, s[i]);
    }
    if (ms <= d - 3) {
      for (i = 0; i <= ms; i++) {
        x1[i] = xs[i];
        y1[i] = ys[i];
        s1[i] = s[i];
      }
      x1[ms + 1] = xs[ms];
      y1[ms + 1] = ys[ms];
      s1[ms + 1] = ms + 1;

      s2[0] = 0;
      for (i = 1; i <= d - 1 - ms; i++) {
        s2[i] = s[i + ms] - ms;
      }

      std::vector<int> tmp_xs(xs.begin() + ms, xs.end());
      std::vector<int> tmp_ys(ys.begin() + ms, ys.end());
      return_val = flutes_wl_LMD(ms + 2, x1, y1, s1, acc)
                   + flutes_wl_LMD(d - ms, tmp_xs, tmp_ys, s2, acc);
      free(score);
      free(penalty);

      return return_val;
    }
  } else {  // (s[0] > s[d-1])
    ms = std::min(s[0], s[1]);
    for (i = 2; i <= d - 1 - ms; i++) {
      ms = std::min(ms, s[i]);
    }
    if (ms >= 2) {
      x1[0] = xs[ms];
      y1[0] = ys[0];
      s1[0] = s[0] - ms + 1;
      for (i = 1; i <= d - 1 - ms; i++) {
        x1[i] = xs[i + ms - 1];
        y1[i] = ys[i];
        s1[i] = s[i] - ms + 1;
      }
      x1[d - ms] = xs[d - 1];
      y1[d - ms] = ys[d - 1 - ms];
      s1[d - ms] = 0;

      s2[0] = ms;
      for (i = 1; i <= ms; i++) {
        s2[i] = s[i + d - 1 - ms];
      }

      std::vector<int> tmp_ys(ys.begin() + d - 1 - ms, ys.end());
      return_val = flutes_wl_LMD(d + 1 - ms, x1, y1, s1, acc)
                   + flutes_wl_LMD(ms + 1, xs, tmp_ys, s2, acc);
      free(score);
      free(penalty);
      return return_val;
    }
  }

  // Find inverse si[] of s[]
  for (r = 0; r < d; r++) {
    si[s[r]] = r;
  }

  // Determine breaking directions and positions dp[]
  lb = (d - 2 * acc + 2) / 4;
  if (lb < 2) {
    lb = 2;
  }
  ub = d - 1 - lb;

  // Compute scores
#define AAWL 0.6
#define BBWL 0.3
  float CCWL = 7.4 / ((d + 10.) * (d - 3.));
  float DDWL = 4.8 / (d - 1);

  // Compute penalty[]
  dx = CCWL * (xs[d - 2] - xs[1]);
  dy = CCWL * (ys[d - 2] - ys[1]);
  for (r = d / 2, pnlty = 0; r >= 0; r--, pnlty += dx) {
    penalty[r] = pnlty, penalty[d - 1 - r] = pnlty;
  }
  for (r = d / 2 - 1, pnlty = dy; r >= 0; r--, pnlty += dy) {
    penalty[s[r]] += pnlty, penalty[s[d - 1 - r]] += pnlty;
  }
  // #define CCWL 0.16
  //     for (r=0; r<d; r++)
  //         penalty[r] = abs(d-1-r-r)*dx + abs(d-1-si[r]-si[r])*dy;

  // Compute distx[], disty[]
  xydiff = (xs[d - 1] - xs[0]) - (ys[d - 1] - ys[0]);
  if (s[0] < s[1]) {
    mins = s[0], maxs = s[1];
  } else {
    mins = s[1], maxs = s[0];
  }
  if (si[0] < si[1]) {
    minsi = si[0], maxsi = si[1];
  } else {
    minsi = si[1], maxsi = si[0];
  }
  for (r = 2; r <= ub; r++) {
    if (s[r] < mins) {
      mins = s[r];
    } else if (s[r] > maxs) {
      maxs = s[r];
    }
    distx[r] = xs[maxs] - xs[mins];
    if (si[r] < minsi) {
      minsi = si[r];
    } else if (si[r] > maxsi) {
      maxsi = si[r];
    }
    disty[r] = ys[maxsi] - ys[minsi] + xydiff;
  }

  if (s[d - 2] < s[d - 1]) {
    mins = s[d - 2], maxs = s[d - 1];
  } else {
    mins = s[d - 1], maxs = s[d - 2];
  }
  if (si[d - 2] < si[d - 1]) {
    minsi = si[d - 2], maxsi = si[d - 1];
  } else {
    minsi = si[d - 1], maxsi = si[d - 2];
  }
  for (r = d - 3; r >= lb; r--) {
    if (s[r] < mins) {
      mins = s[r];
    } else if (s[r] > maxs) {
      maxs = s[r];
    }
    distx[r] += xs[maxs] - xs[mins];
    if (si[r] < minsi) {
      minsi = si[r];
    } else if (si[r] > maxsi) {
      maxsi = si[r];
    }
    disty[r] += ys[maxsi] - ys[minsi];
  }

  nbp = 0;
  for (r = lb; r <= ub; r++) {
    if (si[r] == 0 || si[r] == d - 1) {
      score[nbp] = (xs[r + 1] - xs[r - 1]) - penalty[r]
                   - AAWL * (ys[d - 2] - ys[1]) - DDWL * disty[r];
    } else {
      score[nbp] = (xs[r + 1] - xs[r - 1]) - penalty[r]
                   - BBWL * (ys[si[r] + 1] - ys[si[r] - 1]) - DDWL * disty[r];
    }
    nbp++;

    if (s[r] == 0 || s[r] == d - 1) {
      score[nbp] = (ys[r + 1] - ys[r - 1]) - penalty[s[r]]
                   - AAWL * (xs[d - 2] - xs[1]) - DDWL * distx[r];
    } else {
      score[nbp] = (ys[r + 1] - ys[r - 1]) - penalty[s[r]]
                   - BBWL * (xs[s[r] + 1] - xs[s[r] - 1]) - DDWL * distx[r];
    }
    nbp++;
  }

  if (acc <= 3) {
    newacc = 1;
  } else {
    newacc = acc / 2;
    if (acc >= nbp) {
      acc = nbp - 1;
    }
  }

  minl = (int) INT_MAX;
  for (i = 0; i < acc; i++) {
    maxbp = 0;
    for (bp = 1; bp < nbp; bp++) {
      if (score[maxbp] < score[bp]) {
        maxbp = bp;
      }
    }
    score[maxbp] = -9e9;

#define BreakPt(bp) ((bp) / 2 + lb)
#define BreakInX(bp) ((bp) % 2 == 0)
    p = BreakPt(maxbp);
    // Breaking in p
    if (BreakInX(maxbp)) {  // break in x
      n1 = n2 = 0;
      for (r = 0; r < d; r++) {
        if (s[r] < p) {
          s1[n1] = s[r];
          y1[n1] = ys[r];
          n1++;
        } else if (s[r] > p) {
          s2[n2] = s[r] - p;
          y2[n2] = ys[r];
          n2++;
        } else {  // if (s[r] == p)  i.e.,  r = si[p]
          s1[n1] = p;
          s2[n2] = 0;
          if (r == d - 1 || r == d - 2) {
            y1[n1] = y2[n2] = ys[r - 1];
            extral = ys[r] - ys[r - 1];
          } else if (r == 0 || r == 1) {
            y1[n1] = y2[n2] = ys[r + 1];
            extral = ys[r + 1] - ys[r];
          } else {
            y1[n1] = y2[n2] = ys[r];
            extral = 0;
          }
          n1++;
          n2++;
        }
      }
      std::vector<int> tmp_xs(xs.begin() + p, xs.end());
      ll = extral + flutes_wl_LMD(p + 1, xs, y1, s1, newacc)
           + flutes_wl_LMD(d - p, tmp_xs, y2, s2, newacc);
    } else {  // if (!BreakInX(maxbp))
      n1 = n2 = 0;
      for (r = 0; r < d; r++) {
        if (si[r] < p) {
          s1[si[r]] = n1;
          x1[n1] = xs[r];
          n1++;
        } else if (si[r] > p) {
          s2[si[r] - p] = n2;
          x2[n2] = xs[r];
          n2++;
        } else {  // if (si[r] == p)  i.e.,  r = s[p]
          s1[p] = n1;
          s2[0] = n2;
          if (r == d - 1 || r == d - 2) {
            x1[n1] = x2[n2] = xs[r - 1];
            extral = xs[r] - xs[r - 1];
          } else if (r == 0 || r == 1) {
            x1[n1] = x2[n2] = xs[r + 1];
            extral = xs[r + 1] - xs[r];
          } else {
            x1[n1] = x2[n2] = xs[r];
            extral = 0;
          }
          n1++;
          n2++;
        }
      }
      std::vector<int> tmp_ys(ys.begin() + p, ys.end());
      ll = extral + flutes_wl_LMD(p + 1, x1, ys, s1, newacc)
           + flutes_wl_LMD(d - p, x2, tmp_ys, s2, newacc);
    }
    if (minl > ll) {
      minl = ll;
    }
  }
  return_val = minl;

  free(score);
  free(penalty);
  return return_val;
}

static bool orderx(const point* a, const point* b)
{
  return a->x < b->x;
}

static bool ordery(const point* a, const point* b)
{
  return a->y < b->y;
}

Tree Flute::flute(const std::vector<int>& x, const std::vector<int>& y, int acc)
{
  std::vector<int> xs, ys;
  int minval;
  std::vector<int> s;
  int i, j, minidx;
  struct point *pt, *tmpp;
  Tree t;
  int d = x.size();

  if (d < 2) {
    t.deg = 1;
    t.branch.resize(0);
    t.length = 0;
    return t;
  }

  if (d == 2) {
    t.deg = 2;
    t.length = ADIFF(x[0], x[1]) + ADIFF(y[0], y[1]);
    t.branch.resize(2);
    t.branch[0].x = x[0];
    t.branch[0].y = y[0];
    t.branch[0].n = 1;
    t.branch[1].x = x[1];
    t.branch[1].y = y[1];
    t.branch[1].n = 1;
  } else {
    ensureLUT(d);

    xs.resize(d);
    ys.resize(d);
    s.resize(d);
    pt = (struct point*) malloc(sizeof(struct point) * (d + 1));
    std::vector<point*> ptp(d + 1);

    for (i = 0; i < d; i++) {
      pt[i].x = x[i];
      pt[i].y = y[i];
      ptp[i] = &pt[i];
    }

    // sort x
    if (d < 200) {
      for (i = 0; i < d - 1; i++) {
        minval = ptp[i]->x;
        minidx = i;
        for (j = i + 1; j < d; j++) {
          if (minval > ptp[j]->x) {
            minval = ptp[j]->x;
            minidx = j;
          }
        }
        tmpp = ptp[i];
        ptp[i] = ptp[minidx];
        ptp[minidx] = tmpp;
      }
    } else {
      std::stable_sort(ptp.begin(), ptp.end() - 1, orderx);
    }

#if FLUTE_REMOVE_DUPLICATE_PIN == 1
    ptp[d] = &pt[d];
    ptp[d]->x = ptp[d]->y = -999999;
    j = 0;
    for (i = 0; i < d; i++) {
      for (k = i + 1; ptp[k]->x == ptp[i]->x; k++) {
        if (ptp[k]->y == ptp[i]->y) {  // pins k and i are the same
          break;
        }
      }
      if (ptp[k]->x != ptp[i]->x) {
        ptp[j++] = ptp[i];
      }
    }
    d = j;
#endif

    for (i = 0; i < d; i++) {
      xs[i] = ptp[i]->x;
      ptp[i]->o = i;
    }

    // sort y to find s[]
    if (d < 200) {
      for (i = 0; i < d - 1; i++) {
        minval = ptp[i]->y;
        minidx = i;
        for (j = i + 1; j < d; j++) {
          if (minval > ptp[j]->y) {
            minval = ptp[j]->y;
            minidx = j;
          }
        }
        ys[i] = ptp[minidx]->y;
        s[i] = ptp[minidx]->o;
        ptp[minidx] = ptp[i];
      }
      ys[d - 1] = ptp[d - 1]->y;
      s[d - 1] = ptp[d - 1]->o;
    } else {
      std::stable_sort(ptp.begin(), ptp.end() - 1, ordery);
      for (i = 0; i < d; i++) {
        ys[i] = ptp[i]->y;
        s[i] = ptp[i]->o;
      }
    }

    t = flutes(xs, ys, s, acc);

    free(pt);
  }

  return t;
}

// xs[] and ys[] are coords in x and y in sorted order
// s[] is a list of nodes in increasing y direction
//   if nodes are indexed in the order of increasing x coord
//   i.e., s[i] = s_i as defined in paper
// The points are (xs[s[i]], ys[i]) for i=0..d-1
//             or (xs[i], ys[si[i]]) for i=0..d-1

Tree Flute::flutes_RDP(int d,
                       std::vector<int> xs,
                       std::vector<int> ys,
                       std::vector<int> s,
                       int acc)
{
  int i, j, ss;

  ensureLUT(d);

  for (i = 0; i < d - 1; i++) {
    if (xs[s[i]] == xs[s[i + 1]] && ys[i] == ys[i + 1]) {
      if (s[i] < s[i + 1]) {
        ss = s[i + 1];
      } else {
        ss = s[i];
        s[i] = s[i + 1];
      }
      for (j = i + 2; j < d; j++) {
        ys[j - 1] = ys[j];
        s[j - 1] = s[j];
      }
      for (j = ss + 1; j < d; j++) {
        xs[j - 1] = xs[j];
      }
      for (j = 0; j <= d - 2; j++) {
        if (s[j] > ss) {
          s[j]--;
        }
      }
      i--;
      d--;
    }
  }
  return flutes_ALLD(d, xs, ys, s, acc);
}

// For low-degree, i.e., 2 <= d <= FLUTE_D
Tree Flute::flutes_LD(int d,
                      const std::vector<int>& xs,
                      const std::vector<int>& ys,
                      const std::vector<int>& s)
{
  int k, pi, i, j;
  struct csoln *rlist, *bestrlist;
  int dd[2 * FLUTE_D - 2];  // 0..D-2 for v, D-1..2*D-3 for h
  int minl, sum, l[MPOWV + 1];
  int hflip;
  Tree t;

  t.deg = d;
  t.branch.resize(2 * d - 2);
  if (d == 2) {
    minl = xs[1] - xs[0] + ys[1] - ys[0];
    t.branch[0].x = xs[s[0]];
    t.branch[0].y = ys[0];
    t.branch[0].n = 1;
    t.branch[1].x = xs[s[1]];
    t.branch[1].y = ys[1];
    t.branch[1].n = 1;
  } else if (d == 3) {
    minl = xs[2] - xs[0] + ys[2] - ys[0];
    t.branch[0].x = xs[s[0]];
    t.branch[0].y = ys[0];
    t.branch[0].n = 3;
    t.branch[1].x = xs[s[1]];
    t.branch[1].y = ys[1];
    t.branch[1].n = 3;
    t.branch[2].x = xs[s[2]];
    t.branch[2].y = ys[2];
    t.branch[2].n = 3;
    t.branch[3].x = xs[1];
    t.branch[3].y = ys[1];
    t.branch[3].n = 3;
  } else {
    ensureLUT(d);

    k = 0;
    if (s[0] < s[2]) {
      k++;
    }
    if (s[1] < s[2]) {
      k++;
    }

    for (i = 3; i <= d - 1; i++) {  // p0=0 always, skip i=1 for symmetry
      pi = s[i];
      for (j = d - 1; j > i; j--) {
        if (s[j] < s[i]) {
          pi--;
        }
      }
      k = pi + (i + 1) * k;
    }

    if (k < numgrp[d]) {  // no horizontal flip
      hflip = 0;
      for (i = 1; i <= d - 3; i++) {
        dd[i] = ys[i + 1] - ys[i];
        dd[d - 1 + i] = xs[i + 1] - xs[i];
      }
    } else {
      hflip = 1;
      k = 2 * numgrp[d] - 1 - k;
      for (i = 1; i <= d - 3; i++) {
        dd[i] = ys[i + 1] - ys[i];
        dd[d - 1 + i] = xs[d - 1 - i] - xs[d - 2 - i];
      }
    }

    minl = l[0] = xs[d - 1] - xs[0] + ys[d - 1] - ys[0];
    rlist = (*LUT_)[d][k].get();
    for (i = 0; rlist->seg[i] > 0; i++) {
      minl += dd[rlist->seg[i]];
    }
    bestrlist = rlist;
    l[1] = minl;
    j = 2;
    while (j <= numsoln_[d][k]) {
      rlist++;
      sum = l[rlist->parent];
      for (i = 0; rlist->seg[i] > 0; i++) {
        sum += dd[rlist->seg[i]];
      }
      for (i = 10; rlist->seg[i] > 0; i--) {
        sum -= dd[rlist->seg[i]];
      }
      if (sum < minl) {
        minl = sum;
        bestrlist = rlist;
      }
      l[j++] = sum;
    }

    t.branch[0].x = xs[s[0]];
    t.branch[0].y = ys[0];
    t.branch[1].x = xs[s[1]];
    t.branch[1].y = ys[1];
    for (i = 2; i < d - 2; i++) {
      t.branch[i].x = xs[s[i]];
      t.branch[i].y = ys[i];
      t.branch[i].n = bestrlist->neighbor[i];
    }
    t.branch[d - 2].x = xs[s[d - 2]];
    t.branch[d - 2].y = ys[d - 2];
    t.branch[d - 1].x = xs[s[d - 1]];
    t.branch[d - 1].y = ys[d - 1];
    if (hflip) {
      if (s[1] < s[0]) {
        t.branch[0].n = bestrlist->neighbor[1];
        t.branch[1].n = bestrlist->neighbor[0];
      } else {
        t.branch[0].n = bestrlist->neighbor[0];
        t.branch[1].n = bestrlist->neighbor[1];
      }
      if (s[d - 1] < s[d - 2]) {
        t.branch[d - 2].n = bestrlist->neighbor[d - 1];
        t.branch[d - 1].n = bestrlist->neighbor[d - 2];
      } else {
        t.branch[d - 2].n = bestrlist->neighbor[d - 2];
        t.branch[d - 1].n = bestrlist->neighbor[d - 1];
      }
      for (i = d; i < 2 * d - 2; i++) {
        t.branch[i].x = xs[d - 1 - bestrlist->rowcol[i - d] % 16];
        t.branch[i].y = ys[bestrlist->rowcol[i - d] / 16];
        t.branch[i].n = bestrlist->neighbor[i];
      }
    } else {  // !hflip
      if (s[0] < s[1]) {
        t.branch[0].n = bestrlist->neighbor[1];
        t.branch[1].n = bestrlist->neighbor[0];
      } else {
        t.branch[0].n = bestrlist->neighbor[0];
        t.branch[1].n = bestrlist->neighbor[1];
      }
      if (s[d - 2] < s[d - 1]) {
        t.branch[d - 2].n = bestrlist->neighbor[d - 1];
        t.branch[d - 1].n = bestrlist->neighbor[d - 2];
      } else {
        t.branch[d - 2].n = bestrlist->neighbor[d - 2];
        t.branch[d - 1].n = bestrlist->neighbor[d - 1];
      }
      for (i = d; i < 2 * d - 2; i++) {
        t.branch[i].x = xs[bestrlist->rowcol[i - d] % 16];
        t.branch[i].y = ys[bestrlist->rowcol[i - d] / 16];
        t.branch[i].n = bestrlist->neighbor[i];
      }
    }
  }
  t.length = minl;

  return t;
}

// For medium-degree, i.e., FLUTE_D+1 <= d
Tree Flute::flutes_MD(int d,
                      const std::vector<int>& xs,
                      const std::vector<int>& ys,
                      const std::vector<int>& s,
                      int acc)
{
  float *score, *penalty, pnlty, dx, dy;
  int ms, mins, maxs, minsi, maxsi;
  int i, r, p, maxbp, bestbp = 0, bp, nbp, ub, lb, n1, n2, newacc;
  int nn1 = 0;
  int nn2 = 0;
  std::vector<int> si, s1, s2;
  int degree;
  Tree t, t1, t2, bestt1, bestt2;
  int ll, minl, coord1, coord2;
  std::vector<int> distx, disty;
  int xydiff;
  std::vector<int> x1, x2, y1, y2;

  degree = d + 1;
  score = (float*) malloc(sizeof(float) * (2 * degree));
  penalty = (float*) malloc(sizeof(float) * (degree));

  x1.resize(degree);
  x2.resize(degree);
  y1.resize(degree);
  y2.resize(degree);
  distx.resize(degree);
  disty.resize(degree);
  si.resize(degree);
  s1.resize(degree);
  s2.resize(degree);

  if (s[0] < s[d - 1]) {
    ms = std::max(s[0], s[1]);
    for (i = 2; i <= ms; i++) {
      ms = std::max(ms, s[i]);
    }
    if (ms <= d - 3) {
      for (i = 0; i <= ms; i++) {
        x1[i] = xs[i];
        y1[i] = ys[i];
        s1[i] = s[i];
      }
      x1[ms + 1] = xs[ms];
      y1[ms + 1] = ys[ms];
      s1[ms + 1] = ms + 1;

      s2[0] = 0;
      for (i = 1; i <= d - 1 - ms; i++) {
        s2[i] = s[i + ms] - ms;
      }

      t1 = flutes_LMD(ms + 2, x1, y1, s1, acc);

      std::vector<int> tmp_xs(xs.begin() + ms, xs.end());
      std::vector<int> tmp_ys(ys.begin() + ms, ys.end());
      t2 = flutes_LMD(d - ms, tmp_xs, tmp_ys, s2, acc);
      t = dmergetree(t1, t2);

      free(score);
      free(penalty);

      return t;
    }
  } else {  // (s[0] > s[d-1])
    ms = std::min(s[0], s[1]);
    for (i = 2; i <= d - 1 - ms; i++) {
      ms = std::min(ms, s[i]);
    }
    if (ms >= 2) {
      x1[0] = xs[ms];
      y1[0] = ys[0];
      s1[0] = s[0] - ms + 1;
      for (i = 1; i <= d - 1 - ms; i++) {
        x1[i] = xs[i + ms - 1];
        y1[i] = ys[i];
        s1[i] = s[i] - ms + 1;
      }
      x1[d - ms] = xs[d - 1];
      y1[d - ms] = ys[d - 1 - ms];
      s1[d - ms] = 0;

      s2[0] = ms;
      for (i = 1; i <= ms; i++) {
        s2[i] = s[i + d - 1 - ms];
      }

      t1 = flutes_LMD(d + 1 - ms, x1, y1, s1, acc);

      std::vector<int> tmp_ys(ys.begin() + d - 1 - ms, ys.end());
      t2 = flutes_LMD(ms + 1, xs, tmp_ys, s2, acc);
      t = dmergetree(t1, t2);

      free(score);
      free(penalty);

      return t;
    }
  }

  // Find inverse si[] of s[]
  for (r = 0; r < d; r++) {
    si[s[r]] = r;
  }

  // Determine breaking directions and positions dp[]
  lb = (d - 2 * acc + 2) / 4;
  if (lb < 2) {
    lb = 2;
  }
  ub = d - 1 - lb;

  // Compute scores
#define AA 0.6  // 2.0*BB
#define BB 0.3
  float CC = 7.4 / ((d + 10.) * (d - 3.));
  float DD = 4.8 / (d - 1);

  // Compute penalty[]
  dx = CC * (xs[d - 2] - xs[1]);
  dy = CC * (ys[d - 2] - ys[1]);
  for (r = d / 2, pnlty = 0; r >= 2; r--, pnlty += dx) {
    penalty[r] = pnlty, penalty[d - 1 - r] = pnlty;
  }
  penalty[1] = pnlty, penalty[d - 2] = pnlty;
  penalty[0] = pnlty, penalty[d - 1] = pnlty;
  for (r = d / 2 - 1, pnlty = dy; r >= 2; r--, pnlty += dy) {
    penalty[s[r]] += pnlty, penalty[s[d - 1 - r]] += pnlty;
  }
  penalty[s[1]] += pnlty, penalty[s[d - 2]] += pnlty;
  penalty[s[0]] += pnlty, penalty[s[d - 1]] += pnlty;
  // #define CC 0.16
  // #define v(r) ((r==0||r==1||r==d-2||r==d-1) ? d-3 : abs(d-1-r-r))
  //     for (r=0; r<d; r++)
  //         penalty[r] = v(r)*dx + v(si[r])*dy;

  // Compute distx[], disty[]
  xydiff = (xs[d - 1] - xs[0]) - (ys[d - 1] - ys[0]);
  if (s[0] < s[1]) {
    mins = s[0], maxs = s[1];
  } else {
    mins = s[1], maxs = s[0];
  }
  if (si[0] < si[1]) {
    minsi = si[0], maxsi = si[1];
  } else {
    minsi = si[1], maxsi = si[0];
  }
  for (r = 2; r <= ub; r++) {
    if (s[r] < mins) {
      mins = s[r];
    } else if (s[r] > maxs) {
      maxs = s[r];
    }
    distx[r] = xs[maxs] - xs[mins];
    if (si[r] < minsi) {
      minsi = si[r];
    } else if (si[r] > maxsi) {
      maxsi = si[r];
    }
    disty[r] = ys[maxsi] - ys[minsi] + xydiff;
  }

  if (s[d - 2] < s[d - 1]) {
    mins = s[d - 2], maxs = s[d - 1];
  } else {
    mins = s[d - 1], maxs = s[d - 2];
  }
  if (si[d - 2] < si[d - 1]) {
    minsi = si[d - 2], maxsi = si[d - 1];
  } else {
    minsi = si[d - 1], maxsi = si[d - 2];
  }
  for (r = d - 3; r >= lb; r--) {
    if (s[r] < mins) {
      mins = s[r];
    } else if (s[r] > maxs) {
      maxs = s[r];
    }
    distx[r] += xs[maxs] - xs[mins];
    if (si[r] < minsi) {
      minsi = si[r];
    } else if (si[r] > maxsi) {
      maxsi = si[r];
    }
    disty[r] += ys[maxsi] - ys[minsi];
  }

  nbp = 0;
  for (r = lb; r <= ub; r++) {
    if (si[r] <= 1) {
      score[nbp] = (xs[r + 1] - xs[r - 1]) - penalty[r] - AA * (ys[2] - ys[1])
                   - DD * disty[r];
    } else if (si[r] >= d - 2) {
      score[nbp] = (xs[r + 1] - xs[r - 1]) - penalty[r]
                   - AA * (ys[d - 2] - ys[d - 3]) - DD * disty[r];
    } else {
      score[nbp] = (xs[r + 1] - xs[r - 1]) - penalty[r]
                   - BB * (ys[si[r] + 1] - ys[si[r] - 1]) - DD * disty[r];
    }
    nbp++;

    if (s[r] <= 1) {
      score[nbp] = (ys[r + 1] - ys[r - 1]) - penalty[s[r]]
                   - AA * (xs[2] - xs[1]) - DD * distx[r];
    } else if (s[r] >= d - 2) {
      score[nbp] = (ys[r + 1] - ys[r - 1]) - penalty[s[r]]
                   - AA * (xs[d - 2] - xs[d - 3]) - DD * distx[r];
    } else {
      score[nbp] = (ys[r + 1] - ys[r - 1]) - penalty[s[r]]
                   - BB * (xs[s[r] + 1] - xs[s[r] - 1]) - DD * distx[r];
    }
    nbp++;
  }

  if (acc <= 3) {
    newacc = 1;
  } else {
    newacc = acc / 2;
    if (acc >= nbp) {
      acc = nbp - 1;
    }
  }

  minl = (int) INT_MAX;
  bestt1.branch.clear();
  bestt2.branch.clear();
  for (i = 0; i < acc; i++) {
    maxbp = 0;
    for (bp = 1; bp < nbp; bp++) {
      if (score[maxbp] < score[bp]) {
        maxbp = bp;
      }
    }
    score[maxbp] = -9e9;

#define BreakPt(bp) ((bp) / 2 + lb)
#define BreakInX(bp) ((bp) % 2 == 0)
    p = BreakPt(maxbp);
    // Breaking in p
    if (BreakInX(maxbp)) {  // break in x
      n1 = n2 = 0;
      for (r = 0; r < d; r++) {
        if (s[r] < p) {
          s1[n1] = s[r];
          y1[n1] = ys[r];
          n1++;
        } else if (s[r] > p) {
          s2[n2] = s[r] - p;
          y2[n2] = ys[r];
          n2++;
        } else {  // if (s[r] == p)  i.e.,  r = si[p]
          s1[n1] = p;
          s2[n2] = 0;
          y1[n1] = y2[n2] = ys[r];
          nn1 = n1;
          nn2 = n2;
          n1++;
          n2++;
        }
      }

      t1 = flutes_LMD(p + 1, xs, y1, s1, newacc);

      std::vector<int> tmp_xs(xs.begin() + p, xs.end());
      t2 = flutes_LMD(d - p, tmp_xs, y2, s2, newacc);
      ll = t1.length + t2.length;
      coord1 = t1.branch[t1.branch[nn1].n].y;
      coord2 = t2.branch[t2.branch[nn2].n].y;
      if (t2.branch[nn2].y > std::max(coord1, coord2)) {
        ll -= t2.branch[nn2].y - std::max(coord1, coord2);
      } else if (t2.branch[nn2].y < std::min(coord1, coord2)) {
        ll -= std::min(coord1, coord2) - t2.branch[nn2].y;
      }
    } else {  // if (!BreakInX(maxbp))
      n1 = n2 = 0;
      for (r = 0; r < d; r++) {
        if (si[r] < p) {
          s1[si[r]] = n1;
          x1[n1] = xs[r];
          n1++;
        } else if (si[r] > p) {
          s2[si[r] - p] = n2;
          x2[n2] = xs[r];
          n2++;
        } else {  // if (si[r] == p)  i.e.,  r = s[p]
          s1[p] = n1;
          s2[0] = n2;
          x1[n1] = x2[n2] = xs[r];
          n1++;
          n2++;
        }
      }

      t1 = flutes_LMD(p + 1, x1, ys, s1, newacc);

      std::vector<int> tmp_ys(ys.begin() + p, ys.end());
      t2 = flutes_LMD(d - p, x2, tmp_ys, s2, newacc);
      ll = t1.length + t2.length;
      coord1 = t1.branch[t1.branch[p].n].x;
      coord2 = t2.branch[t2.branch[0].n].x;
      if (t2.branch[0].x > std::max(coord1, coord2)) {
        ll -= t2.branch[0].x - std::max(coord1, coord2);
      } else if (t2.branch[0].x < std::min(coord1, coord2)) {
        ll -= std::min(coord1, coord2) - t2.branch[0].x;
      }
    }
    if (minl > ll) {
      minl = ll;
      bestt1 = t1;
      bestt2 = t2;
      bestbp = maxbp;
    }
  }

#if FLUTE_LOCAL_REFINEMENT == 1
  if (BreakInX(bestbp)) {
    t = hmergetree(std::move(bestt1), std::move(bestt2), s);
    local_refinement(degree, &t, si[BreakPt(bestbp)]);
  } else {
    t = vmergetree(std::move(bestt1), std::move(bestt2));
    local_refinement(degree, &t, BreakPt(bestbp));
  }
#else
  if (BreakInX(bestbp)) {
    t = hmergetree(bestt1, bestt2, s);
  } else {
    t = vmergetree(bestt1, bestt2);
  }
#endif

  free(score);
  free(penalty);

  return t;
}

Tree Flute::dmergetree(Tree t1, Tree t2)
{
  int i, d, prev, curr, next, offset1, offset2;
  Tree t;

  t.deg = d = t1.deg + t2.deg - 2;
  t.length = t1.length + t2.length;
  t.branch.resize(2 * d - 2);
  offset1 = t2.deg - 2;
  offset2 = 2 * t1.deg - 4;

  for (i = 0; i <= t1.deg - 2; i++) {
    t.branch[i].x = t1.branch[i].x;
    t.branch[i].y = t1.branch[i].y;
    t.branch[i].n = t1.branch[i].n + offset1;
  }
  for (i = t1.deg - 1; i <= d - 1; i++) {
    t.branch[i].x = t2.branch[i - t1.deg + 2].x;
    t.branch[i].y = t2.branch[i - t1.deg + 2].y;
    t.branch[i].n = t2.branch[i - t1.deg + 2].n + offset2;
  }
  for (i = d; i <= d + t1.deg - 3; i++) {
    t.branch[i].x = t1.branch[i - offset1].x;
    t.branch[i].y = t1.branch[i - offset1].y;
    t.branch[i].n = t1.branch[i - offset1].n + offset1;
  }
  for (i = d + t1.deg - 2; i <= 2 * d - 3; i++) {
    t.branch[i].x = t2.branch[i - offset2].x;
    t.branch[i].y = t2.branch[i - offset2].y;
    t.branch[i].n = t2.branch[i - offset2].n + offset2;
  }

  prev = t2.branch[0].n + offset2;
  curr = t1.branch[t1.deg - 1].n + offset1;
  next = t.branch[curr].n;
  while (curr != next) {
    t.branch[curr].n = prev;
    prev = curr;
    curr = next;
    next = t.branch[curr].n;
  }
  t.branch[curr].n = prev;

  return t;
}

Tree Flute::hmergetree(Tree t1, Tree t2, const std::vector<int>& s)
{
  int i, prev, curr, next, extra, offset1, offset2;
  int p, n1, n2;
  int nn1 = 0;
  int nn2 = 0;
  int ii = 0;

  int coord1, coord2;
  Tree t;

  t.deg = t1.deg + t2.deg - 1;
  t.length = t1.length + t2.length;
  t.branch.resize(2 * t.deg - 2);
  offset1 = t2.deg - 1;
  offset2 = 2 * t1.deg - 3;

  p = t1.deg - 1;
  n1 = n2 = 0;
  for (i = 0; i < t.deg; i++) {
    if (s[i] < p) {
      t.branch[i].x = t1.branch[n1].x;
      t.branch[i].y = t1.branch[n1].y;
      t.branch[i].n = t1.branch[n1].n + offset1;
      n1++;
    } else if (s[i] > p) {
      t.branch[i].x = t2.branch[n2].x;
      t.branch[i].y = t2.branch[n2].y;
      t.branch[i].n = t2.branch[n2].n + offset2;
      n2++;
    } else {
      t.branch[i].x = t2.branch[n2].x;
      t.branch[i].y = t2.branch[n2].y;
      t.branch[i].n = t2.branch[n2].n + offset2;
      nn1 = n1;
      nn2 = n2;
      ii = i;
      n1++;
      n2++;
    }
  }
  for (i = t.deg; i <= t.deg + t1.deg - 3; i++) {
    t.branch[i].x = t1.branch[i - offset1].x;
    t.branch[i].y = t1.branch[i - offset1].y;
    t.branch[i].n = t1.branch[i - offset1].n + offset1;
  }
  for (i = t.deg + t1.deg - 2; i <= 2 * t.deg - 4; i++) {
    t.branch[i].x = t2.branch[i - offset2].x;
    t.branch[i].y = t2.branch[i - offset2].y;
    t.branch[i].n = t2.branch[i - offset2].n + offset2;
  }
  extra = 2 * t.deg - 3;
  coord1 = t1.branch[t1.branch[nn1].n].y;
  coord2 = t2.branch[t2.branch[nn2].n].y;
  if (t2.branch[nn2].y > std::max(coord1, coord2)) {
    t.branch[extra].y = std::max(coord1, coord2);
    t.length -= t2.branch[nn2].y - t.branch[extra].y;
  } else if (t2.branch[nn2].y < std::min(coord1, coord2)) {
    t.branch[extra].y = std::min(coord1, coord2);
    t.length -= t.branch[extra].y - t2.branch[nn2].y;
  } else {
    t.branch[extra].y = t2.branch[nn2].y;
  }
  t.branch[extra].x = t2.branch[nn2].x;
  t.branch[extra].n = t.branch[ii].n;
  t.branch[ii].n = extra;

  prev = extra;
  curr = t1.branch[nn1].n + offset1;
  next = t.branch[curr].n;
  while (curr != next) {
    t.branch[curr].n = prev;
    prev = curr;
    curr = next;
    next = t.branch[curr].n;
  }
  t.branch[curr].n = prev;

  return t;
}

Tree Flute::vmergetree(Tree t1, Tree t2)
{
  int i, prev, curr, next, extra, offset1, offset2;
  int coord1, coord2;
  Tree t;

  t.deg = t1.deg + t2.deg - 1;
  t.length = t1.length + t2.length;
  t.branch.resize(2 * t.deg - 2);
  offset1 = t2.deg - 1;
  offset2 = 2 * t1.deg - 3;

  for (i = 0; i <= t1.deg - 2; i++) {
    t.branch[i].x = t1.branch[i].x;
    t.branch[i].y = t1.branch[i].y;
    t.branch[i].n = t1.branch[i].n + offset1;
  }
  for (i = t1.deg - 1; i <= t.deg - 1; i++) {
    t.branch[i].x = t2.branch[i - t1.deg + 1].x;
    t.branch[i].y = t2.branch[i - t1.deg + 1].y;
    t.branch[i].n = t2.branch[i - t1.deg + 1].n + offset2;
  }
  for (i = t.deg; i <= t.deg + t1.deg - 3; i++) {
    t.branch[i].x = t1.branch[i - offset1].x;
    t.branch[i].y = t1.branch[i - offset1].y;
    t.branch[i].n = t1.branch[i - offset1].n + offset1;
  }
  for (i = t.deg + t1.deg - 2; i <= 2 * t.deg - 4; i++) {
    t.branch[i].x = t2.branch[i - offset2].x;
    t.branch[i].y = t2.branch[i - offset2].y;
    t.branch[i].n = t2.branch[i - offset2].n + offset2;
  }
  extra = 2 * t.deg - 3;
  coord1 = t1.branch[t1.branch[t1.deg - 1].n].x;
  coord2 = t2.branch[t2.branch[0].n].x;
  if (t2.branch[0].x > std::max(coord1, coord2)) {
    t.branch[extra].x = std::max(coord1, coord2);
    t.length -= t2.branch[0].x - t.branch[extra].x;
  } else if (t2.branch[0].x < std::min(coord1, coord2)) {
    t.branch[extra].x = std::min(coord1, coord2);
    t.length -= t.branch[extra].x - t2.branch[0].x;
  } else {
    t.branch[extra].x = t2.branch[0].x;
  }
  t.branch[extra].y = t2.branch[0].y;
  t.branch[extra].n = t.branch[t1.deg - 1].n;
  t.branch[t1.deg - 1].n = extra;

  prev = extra;
  curr = t1.branch[t1.deg - 1].n + offset1;
  next = t.branch[curr].n;
  while (curr != next) {
    t.branch[curr].n = prev;
    prev = curr;
    curr = next;
    next = t.branch[curr].n;
  }
  t.branch[curr].n = prev;

  return t;
}

void Flute::local_refinement(int deg, Tree* tp, int p)
{
  int d, dd, i, ii, j, prev, curr, next, root;
  std::vector<int> SteinerPin, index, ss;
  int degree;
  std::vector<int> x, xs, ys;
  Tree tt;

  degree = deg + 1;
  SteinerPin.resize(2 * degree);
  index.resize(2 * degree);
  x.resize(degree);
  xs.resize(degree);
  ys.resize(degree);
  ss.resize(degree);

  d = tp->deg;
  root = tp->branch[p].n;

  // Reverse edges to point to root
  prev = root;
  curr = tp->branch[prev].n;
  next = tp->branch[curr].n;
  while (curr != next) {
    tp->branch[curr].n = prev;
    prev = curr;
    curr = next;
    next = tp->branch[curr].n;
  }
  tp->branch[curr].n = prev;
  tp->branch[root].n = root;

  // Find Steiner nodes that are at pins
  for (i = d; i <= 2 * d - 3; i++) {
    SteinerPin[i] = -1;
  }
  for (i = 0; i < d; i++) {
    next = tp->branch[i].n;
    if (tp->branch[i].x == tp->branch[next].x
        && tp->branch[i].y == tp->branch[next].y) {
      SteinerPin[next] = i;  // Steiner 'next' at Pin 'i'
    }
  }
  SteinerPin[root] = p;

  // Find pins that are directly connected to root
  dd = 0;
  for (i = 0; i < d; i++) {
    curr = tp->branch[i].n;
    if (SteinerPin[curr] == i) {
      curr = tp->branch[curr].n;
    }
    while (SteinerPin[curr] < 0) {
      curr = tp->branch[curr].n;
    }
    if (curr == root) {
      x[dd] = tp->branch[i].x;
      if (SteinerPin[tp->branch[i].n] == i && tp->branch[i].n != root) {
        index[dd++] = tp->branch[i].n;  // Steiner node
      } else {
        index[dd++] = i;  // Pin
      }
    }
  }

  if (4 <= dd && dd <= FLUTE_D) {
    // Find Steiner nodes that are directly connected to root
    ii = dd;
    for (i = 0; i < dd; i++) {
      curr = tp->branch[index[i]].n;
      while (SteinerPin[curr] < 0) {
        index[ii++] = curr;
        SteinerPin[curr] = INT_MAX;
        curr = tp->branch[curr].n;
      }
    }
    index[ii] = root;

    for (ii = 0; ii < dd; ii++) {
      ss[ii] = 0;
      for (j = 0; j < ii; j++) {
        if (x[j] < x[ii]) {
          ss[ii]++;
        }
      }
      for (j = ii + 1; j < dd; j++) {
        if (x[j] <= x[ii]) {
          ss[ii]++;
        }
      }
      xs[ss[ii]] = x[ii];
      ys[ii] = tp->branch[index[ii]].y;
    }

    tt = flutes_LD(dd, xs, ys, ss);

    // Find new wirelength
    tp->length += tt.length;
    for (ii = 0; ii < 2 * dd - 3; ii++) {
      i = index[ii];
      j = tp->branch[i].n;
      tp->length -= ADIFF(tp->branch[i].x, tp->branch[j].x)
                    + ADIFF(tp->branch[i].y, tp->branch[j].y);
    }

    // Copy tt into t
    for (ii = 0; ii < dd; ii++) {
      tp->branch[index[ii]].n = index[tt.branch[ii].n];
    }
    for (; ii <= 2 * dd - 3; ii++) {
      tp->branch[index[ii]].x = tt.branch[ii].x;
      tp->branch[index[ii]].y = tt.branch[ii].y;
      tp->branch[index[ii]].n = index[tt.branch[ii].n];
    }
  }
}

int Flute::wirelength(Tree t)
{
  int i, j;
  int l = 0;

  for (i = 0; i < 2 * t.deg - 2; i++) {
    j = t.branch[i].n;
    l += ADIFF(t.branch[i].x, t.branch[j].x)
         + ADIFF(t.branch[i].y, t.branch[j].y);
  }

  return l;
}

// Output in a format that can be plotted by gnuplot
void Flute::plottree(Tree t)
{
  int i;

  for (i = 0; i < 2 * t.deg - 2; i++) {
    printf("%d %d\n", t.branch[i].x, t.branch[i].y);
    printf("%d %d\n\n", t.branch[t.branch[i].n].x, t.branch[t.branch[i].n].y);
  }
}

}  // namespace flt

}  // namespace stt
