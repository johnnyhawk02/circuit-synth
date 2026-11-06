// Per-sample hot path of the circuit solver (src/solver.js), compiled to WebAssembly.
// Build: sh kernel/build.sh  (clang --target=wasm32, no libc)
//
// The JS side builds the netlist, the matrices (P, Q, CP) and does everything that is not
// per-sample. Each sample it calls newton() (Newton on the nonlinear port voltages,
// v = P·e + Q·f(v)) and accept_caps() (capacitor states from the solution).
// Everything lives in the module's linear memory; Ctx holds the pointers.

typedef struct {
  int m, k, nc, cols, ncross, luOK, na, iters, factors;
  double *P, *Q, *CP, *qmax;
  double *e, *eLast, *vp, *fp, *pl, *dd, *dv, *r, *Ja, *ra, *ddF, *tolv, *xdF;
  int *act, *pos, *piv, *keep, *isAct;
  int *kind, *partner;
  double *is, *nvt, *vcrit, *vsat, *kk, *vos, *tanhStep, *beta, *vto, *soft, *xd;
  int *crossI, *crossX;
  double *capGeq, *capIeq, *capV, *capI; int *capBe;
} Ctx;

// no libc: the compiler may still emit calls to these
void *memset(void *d, int v, unsigned long n) { unsigned char *p = d; while (n--) *p++ = (unsigned char)v; return d; }
void *memcpy(void *d, const void *s, unsigned long n) { unsigned char *p = d; const unsigned char *q = s; while (n--) *p++ = *q++; return d; }

enum { EXP = 0, TANH = 1, JFET = 2, JFET3 = 3, ZERO = 4, CLAMP = 5 };
#define GMIN 1e-12

static inline double fabs_(double x) { return __builtin_fabs(x); }

// exp: range reduction to |r| <= ln2/2, Taylor to r^12 (rel. error ~2e-16)
static double fexp(double x) {
  if (x > 708) x = 708;
  if (x < -708) return 0;
  double n = __builtin_rint(x * 1.4426950408889634);
  double r = (x - n * 0.693145751953125) - n * 1.4286068203094173e-06;
  double p = 1 + r * (1 + r * (0.5 + r * (1.0 / 6 + r * (1.0 / 24 + r * (1.0 / 120 + r * (1.0 / 720 +
             r * (1.0 / 5040 + r * (1.0 / 40320 + r * (1.0 / 362880 + r * (1.0 / 3628800 +
             r * (1.0 / 39916800 + r * (1.0 / 479001600))))))))))));
  union { double d; long long i; } u; u.i = (long long)(n + 1023) << 52;
  return p * u.d;
}
// natural log of x > 0: x = m·2^e with m in [√½, √2), log m = 2·atanh((m-1)/(m+1))
static double flog(double x) {
  if (!(x > 0)) return -1e300;
  union { double d; long long i; } u; u.d = x;
  int e = (int)((u.i >> 52) & 0x7ff) - 1023;
  u.i = (u.i & 0x000fffffffffffffLL) | 0x3ff0000000000000LL;   // m in [1, 2)
  double m = u.d;
  if (m > 1.4142135623730951) { m *= 0.5; e++; }
  double s = (m - 1) / (m + 1), s2 = s * s;
  double t = s * (2 + s2 * (2.0 / 3 + s2 * (2.0 / 5 + s2 * (2.0 / 7 + s2 * (2.0 / 9 + s2 * (2.0 / 11 +
             s2 * (2.0 / 13 + s2 * (2.0 / 15 + s2 * (2.0 / 17 + s2 * (2.0 / 19 + s2 * (2.0 / 21)))))))))));
  return t + e * 0.6931471805599453;
}
static inline double flog1p(double y) {   // used only for y = exp(a) >= 0
  if (y < 1e-5) return y * (1 - y * (0.5 - y * (1.0 / 3)));
  return flog(1 + y);
}
static inline double ftanh(double x) {
  if (x > 20) return 1;
  if (x < -20) return -1;
  if (fabs_(x) < 1e-3) { double x2 = x * x; return x * (1 - x2 * (1.0 / 3 - x2 * (2.0 / 15))); }
  double e = fexp(2 * x);
  return (e - 1) / (e + 1);
}

static double pnjlim(double vnew, double vold, double vt, double vcrit) {
  if (vnew > vcrit && fabs_(vnew - vold) > 2 * vt) {
    if (vold > 0) {
      double arg = 1 + (vnew - vold) / vt;
      vnew = arg > 0 ? vold + vt * flog(arg) : vcrit;
    } else {
      vnew = vt * flog(vnew / vt);
    }
  }
  return vnew;
}

// Shichman-Hodges N-JFET, vds >= 0, softplus-smoothed overdrive: id, d/dvgs, d/dvds
static void jfetId(double beta, double vto, double soft, double vgs, double vds, double *o) {
  double x = (vgs - vto) / soft;
  double vov = x > 30 ? soft * x : soft * flog1p(fexp(x));
  double dvov = x > 30 ? 1 : 1 / (1 + fexp(-x));
  if (vds < vov) { o[0] = beta * (2 * vov - vds) * vds; o[1] = 2 * beta * vds * dvov; o[2] = 2 * beta * (vov - vds); }
  else { o[0] = beta * vov * vov; o[1] = 2 * beta * vov * dvov; o[2] = 0; }
}

static void evalPort(Ctx *c, int i, double v, double *f, double *d) {
  switch (c->kind[i]) {
    case EXP: {
      double nvt = c->nvt[i], e = fexp(v / nvt);
      *f = c->is[i] * (e - 1) + GMIN * v;
      *d = c->is[i] * e / nvt + GMIN;
      break;
    }
    case TANH: {
      double k = c->kk[i], vs = c->vsat[i], t = ftanh(k * (v + c->vos[i]));
      *f = vs * t;
      *d = vs * k * (1 - t * t) + 1e-9;
      break;
    }
    case JFET3: {
      double vgs = c->vp[c->partner[i]], vds = v, o[3];
      if (vds >= 0) jfetId(c->beta[i], c->vto[i], c->soft[i], vgs, vds, o);
      else {
        double q[3]; jfetId(c->beta[i], c->vto[i], c->soft[i], vgs - vds, -vds, q);
        o[0] = -q[0]; o[1] = -q[1]; o[2] = q[1] + q[2];
      }
      *f = o[0] + GMIN * vds; *d = o[2] + GMIN; c->xd[i] = o[1];
      break;
    }
    case ZERO: *f = 0; *d = 0; break;
    case CLAMP: {
      double s = c->soft[i], vs = c->vsat[i], a = (v - vs) / s, b = (-v - vs) / s;
      double spa = a > 30 ? a : flog1p(fexp(a)), spb = b > 30 ? b : flog1p(fexp(b));
      double sa = a > 30 ? 1 : 1 / (1 + fexp(-a)), sb = b > 30 ? 1 : 1 / (1 + fexp(-b));
      *f = c->kk[i] * s * (spa - spb); *d = c->kk[i] * (sa + sb) + GMIN;
      break;
    }
    case JFET: {
      double s = c->soft[i], x = (v - c->vto[i]) / s;
      double sp = x > 30 ? x : flog1p(fexp(x));
      double sig = x > 30 ? 1 : 1 / (1 + fexp(-x));
      double u = s * sp;
      *f = c->beta[i] * u * u + GMIN * v;
      *d = 2 * c->beta[i] * u * sig + GMIN;
      break;
    }
  }
}

static int smallLU(double *A, int m, int *piv) {
  for (int k = 0; k < m; k++) {
    int p = k; double max = fabs_(A[k * m + k]);
    for (int i = k + 1; i < m; i++) { double v = fabs_(A[i * m + k]); if (v > max) { max = v; p = i; } }
    if (!(max > 1e-300)) return 0;
    piv[k] = p;
    if (p != k) for (int j = 0; j < m; j++) { double t = A[k * m + j]; A[k * m + j] = A[p * m + j]; A[p * m + j] = t; }
    double inv = 1 / A[k * m + k];
    A[k * m + k] = inv;   // store the reciprocal: the solve multiplies
    for (int i = k + 1; i < m; i++) {
      double f = (A[i * m + k] *= inv);
      if (f != 0) { double *ri = A + i * m, *rk = A + k * m; for (int j = k + 1; j < m; j++) ri[j] -= f * rk[j]; }
    }
  }
  return 1;
}
static void smallLUSolve(const double *A, int m, const int *piv, double *b) {
  for (int k = 0; k < m; k++) { int p = piv[k]; if (p != k) { double t = b[k]; b[k] = b[p]; b[p] = t; } }
  for (int i = 1; i < m; i++) { double s = b[i]; const double *a = A + i * m; for (int j = 0; j < i; j++) s -= a[j] * b[j]; b[i] = s; }
  for (int i = m - 1; i >= 0; i--) { double s = b[i]; const double *a = A + i * m; for (int j = i + 1; j < m; j++) s -= a[j] * b[j]; b[i] = s * a[i]; }
}

// Jacobian I - Q·D (+ JFET cross terms). Ports with negligible slope leave an identity
// column: only the active block is factored, the rest is back-substituted.
static void fillJa(Ctx *c, int a) {
  int m = c->m; double *Ja = c->Ja, *Q = c->Q, *ddF = c->ddF; int *act = c->act;
  for (int ii = 0; ii < a; ii++) {
    const double *q = Q + act[ii] * m; double *row = Ja + ii * a;
    for (int jj = 0; jj < a; jj++) { int j = act[jj]; row[jj] = -q[j] * ddF[j]; }
    row[ii] += 1;
  }
  for (int x = 0; x < c->ncross; x++) {
    int ci = c->crossI[x], xa = c->pos[c->crossX[x]]; double xd = c->xdF[x];
    for (int ii = 0; ii < a; ii++) Ja[ii * a + xa] -= Q[act[ii] * m + ci] * xd;
  }
}
static int factor(Ctx *c) {
  int m = c->m, a = 0;
  for (int j = 0; j < m; j++) if (c->keep[j] || c->dd[j] * c->qmax[j] > 1e-6) { c->pos[j] = a; c->act[a++] = j; }
  c->na = a;
  for (int j = 0; j < m; j++) c->ddF[j] = c->dd[j];
  for (int x = 0; x < c->ncross; x++) c->xdF[x] = c->xd[c->crossI[x]];
  fillJa(c, a);
  return smallLU(c->Ja, a, c->piv);
}
static void solve(Ctx *c, double *r) {
  int m = c->m, a = c->na; double *ra = c->ra, *Q = c->Q, *ddF = c->ddF; int *act = c->act;
  for (int ii = 0; ii < a; ii++) ra[ii] = r[act[ii]];
  smallLUSolve(c->Ja, a, c->piv, ra);
  if (a == m) { for (int i = 0; i < m; i++) r[i] = ra[i]; return; }
  int *isAct = c->isAct;
  for (int i = 0; i < m; i++) isAct[i] = 0;
  for (int ii = 0; ii < a; ii++) isAct[act[ii]] = 1;
  for (int i = 0; i < m; i++) {
    if (isAct[i]) continue;
    double s = r[i]; const double *q = Q + i * m;
    for (int jj = 0; jj < a; jj++) { int j = act[jj]; s += q[j] * ddF[j] * ra[jj]; }
    r[i] = s;
  }
  for (int ii = 0; ii < a; ii++) r[act[ii]] = ra[ii];
}

__attribute__((export_name("newton")))
int newton(Ctx *c, int maxIter) {
  int m = c->m, k = c->k;
  double *P = c->P, *Q = c->Q, *e = c->e, *pl = c->pl, *vp = c->vp, *fp = c->fp, *dd = c->dd, *r = c->r, *dv = c->dv, *tolv = c->tolv;
  for (int i = 0; i < m; i++) { double s = 0; const double *p = P + i * k; for (int j = 0; j < k; j++) s += p[j] * e[j]; pl[i] = s; }
  if (m == 0) return 1;
  double prevErr = 1e300; int refactor = 1;
  for (int it = 0; it < maxIter; it++) {
    for (int i = 0; i < m; i++) evalPort(c, i, vp[i], fp + i, dd + i);
    for (int i = 0; i < m; i++) {   // residual F = v - pl - Q f
      double s = 0; const double *q = Q + i * m;
      for (int j = 0; j < m; j++) s += q[j] * fp[j];
      r[i] = -(vp[i] - pl[i] - s);
    }
    c->iters++;
    if (refactor) {
      if (!factor(c)) { c->luOK = 0; break; }
      c->luOK = 1; refactor = 0; c->factors++;
    }
    solve(c, r);
    int converged = 1; double err = 0;
    for (int i = 0; i < m; i++) {
      int kind = c->kind[i]; double step = r[i], old = vp[i], vn = old + step;
      if (kind == EXP) vn = pnjlim(vn, old, c->nvt[i], c->vcrit[i]);
      else if (kind == TANH) { double lim = c->tanhStep[i]; if (step > lim) vn = old + lim; else if (step < -lim) vn = old - lim; }
      else { double lim = kind == CLAMP ? 0.5 : 1; if (step > lim) vn = old + lim; else if (step < -lim) vn = old - lim; }
      double d = vn - old, ad = fabs_(d);
      if (!(ad <= tolv[i])) converged = 0;
      if (ad > err) err = ad;
      dv[i] = d; vp[i] = vn;
    }
    // ill-conditioned current-fed tail nodes: stagnation below 2 µV counts as converged
    if (!converged && err < 2e-6 && err > 0.25 * prevErr) converged = 1;
    if (!converged && (it > 0 ? err > 0.1 * prevErr : err > 1e-2)) refactor = 1;
    prevErr = err;
    if (converged) {
      // last update was tiny: move the port currents along their slopes
      for (int i = 0; i < m; i++) fp[i] += dd[i] * dv[i];
      for (int x = 0; x < c->ncross; x++) { int ci = c->crossI[x]; fp[ci] += c->xd[ci] * dv[c->crossX[x]]; }
      return 1;
    }
  }
  for (int i = 0; i < m; i++) evalPort(c, i, vp[i], fp + i, dd + i);
  c->luOK = 0;
  return 0;
}

// capacitor voltages/currents from the accepted solution; trapezoidal or backward Euler
__attribute__((export_name("accept_caps")))
void accept_caps(Ctx *c, int dc) {
  int k = c->k, m = c->m, cols = c->cols;
  double *e = c->e, *fp = c->fp;
  for (int j = 0; j < k; j++) c->eLast[j] = e[j];
  for (int i = 0; i < c->nc; i++) {
    const double *R = c->CP + i * cols; double v = 0;
    for (int j = 0; j < k; j++) v += R[j] * e[j];
    for (int j = 0; j < m; j++) v += R[k + j] * fp[j];
    double g = c->capGeq[i];
    double cur = dc ? 0 : g * v - c->capIeq[i];
    c->capI[i] = cur; c->capV[i] = v;
    c->capIeq[i] = dc ? 0 : (c->capBe[i] ? g * v : g * v + cur);
  }
}

// e = [sources, capacitor history]
__attribute__((export_name("load_caps")))
void load_caps(Ctx *c, int ns) {
  for (int j = 0; j < c->nc; j++) c->e[ns + j] = c->capIeq[j];
}

// memory for the JS allocator starts after the module's own data
extern unsigned char __heap_base;
__attribute__((export_name("heap_base")))
int heap_base(void) { return (int)&__heap_base; }
