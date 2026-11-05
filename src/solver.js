// Component-level circuit solver: Modified Nodal Analysis (MNA)
// + Newton-Raphson per sample + trapezoidal integration for capacitors.
// Works in Node (ES modules) and in the browser / AudioWorklet.

const VT = 0.025852; // thermal voltage at 27 °C
const GMIN = 1e-12;  // tiny conductance across every junction, helps convergence

// ---------- small dense LU solve (circuits here are < 50 unknowns) ----------
function solveInPlace(A, b, n) {
  for (let k = 0; k < n; k++) {
    // partial pivoting
    let p = k, max = Math.abs(A[k * n + k]);
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(A[i * n + k]);
      if (v > max) { max = v; p = i; }
    }
    if (max < 1e-30) throw new Error('Singular matrix (floating node?) at row ' + k);
    if (p !== k) {
      for (let j = 0; j < n; j++) { const t = A[k * n + j]; A[k * n + j] = A[p * n + j]; A[p * n + j] = t; }
      const t = b[k]; b[k] = b[p]; b[p] = t;
    }
    const piv = A[k * n + k];
    for (let i = k + 1; i < n; i++) {
      const f = A[i * n + k] / piv;
      if (f === 0) continue;
      for (let j = k; j < n; j++) A[i * n + j] -= f * A[k * n + j];
      b[i] -= f * b[k];
    }
  }
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let j = i + 1; j < n; j++) s -= A[i * n + j] * b[j];
    b[i] = s / A[i * n + i];
  }
  return b;
}

// SPICE-style junction voltage limiting: stops exp() from blowing up mid-iteration
function pnjlim(vnew, vold, vt, vcrit) {
  if (vnew > vcrit && Math.abs(vnew - vold) > 2 * vt) {
    if (vold > 0) {
      const arg = 1 + (vnew - vold) / vt;
      vnew = arg > 0 ? vold + vt * Math.log(arg) : vcrit;
    } else {
      vnew = vt * Math.log(vnew / vt);
    }
  }
  return vnew;
}

// ---------- the circuit ----------
export class Circuit {
  constructor() {
    this.nodeIndex = new Map([['0', -1], ['gnd', -1]]);
    this.nodeCount = 0;
    this.branchCount = 0;
    this.elements = [];
    this.inputs = new Map(); // name -> voltage source element, settable per sample
    this.maxStep = 3; // volts per Newton iteration
  }

  node(name) {
    if (!this.nodeIndex.has(name)) this.nodeIndex.set(name, this.nodeCount++);
    return this.nodeIndex.get(name);
  }
  branch() { return this.branchCount++; } // resolved to matrix row in finalize()

  add(el) { this.elements.push(el); return el; }

  // --- component helpers (nodes are names, '0' is ground) ---
  R(a, b, ohms) { return this.add(new Resistor(this.node(a), this.node(b), ohms)); }
  C(a, b, farads, v0 = 0) { return this.add(new Capacitor(this.node(a), this.node(b), farads, v0)); }
  V(p, n, volts, name) {
    const el = this.add(new VSource(this.node(p), this.node(n), volts, this.branch()));
    if (name) this.inputs.set(name, el);
    return el;
  }
  D(a, k, params) { return this.add(new Diode(this.node(a), this.node(k), params)); }
  NPN(c, b, e, params) { return this.add(new BJT(this.node(c), this.node(b), this.node(e), +1, params)); }
  PNP(c, b, e, params) { return this.add(new BJT(this.node(c), this.node(b), this.node(e), -1, params)); }
  // Op-amp macro-model: out = vsat * tanh(gain * (in+ - in-) / vsat), with rails ±vsat
  OPAMP(inp, inn, out, params) {
    return this.add(new OpAmp(this.node(inp), this.node(inn), this.node(out), this.branch(), params));
  }

  set(name, volts) { this.inputs.get(name).volts = volts; }
  v(name) { const i = this.nodeIndex.get(name); return i < 0 ? 0 : this.x[i]; }

  finalize(sampleRate) {
    this.n = this.nodeCount + this.branchCount;
    for (const el of this.elements) if (el.br !== undefined) el.row = this.nodeCount + el.br;
    this.A = new Float64Array(this.n * this.n);
    this.b = new Float64Array(this.n);
    this.x = new Float64Array(this.n);
    this.xOld = new Float64Array(this.n);
    this.Abase = new Float64Array(this.n * this.n);
    this.nonlinear = this.elements.filter(e => e.nonlinear);
    this.linear = this.elements.filter(e => !e.nonlinear);
    this.dynamic = this.elements.filter(e => e.dynamic);
    this.setSampleRate(sampleRate);
    this.stats = { steps: 0, iters: 0, failures: 0 };
    return this;
  }

  setSampleRate(sr) { this.h = 1 / sr; this.setStep(this.h); }
  setStep(h) { for (const el of this.dynamic || this.elements) el.setStep?.(h); this.baseDirty = true; }
  // anything that changes a linear part's value (e.g. a pot) must call this
  invalidate() { this.baseDirty = true; }

  // Solve the DC operating point (caps open) by ramping supplies in.
  dcop() {
    for (const el of this.dynamic) el.dcMode = true;
    this.baseDirty = true;
    const srcs = this.elements.filter(e => e instanceof VSource);
    const full = srcs.map(s => s.volts);
    for (let k = 1; k <= 20; k++) {
      srcs.forEach((s, i) => (s.volts = full[i] * k / 20));
      this.newton(200);
    }
    for (const el of this.dynamic) { el.dcMode = false; el.accept(this.x); }
    this.baseDirty = true;
  }

  newton(maxIter = 50, tol = 1e-6) {
    const { A, b, x, n } = this;
    for (let it = 0; it < maxIter; it++) {
      if (this.baseDirty) {
        this.Abase.fill(0);
        for (const el of this.linear) el.stampA(this.Abase, n);
        this.baseDirty = false;
      }
      A.set(this.Abase); b.fill(0);
      for (const el of this.linear) el.stampB(b);
      for (const el of this.nonlinear) el.stamp(A, b, x, n);
      this.xOld.set(x);
      solveInPlace(A, b, n);
      let converged = true;
      // damping: cap the largest node-voltage move per iteration (helps comparators flip)
      let maxDv = 0;
      for (let i = 0; i < this.nodeCount; i++) maxDv = Math.max(maxDv, Math.abs(b[i] - x[i]));
      const alpha = maxDv > this.maxStep ? this.maxStep / maxDv : 1;
      for (let i = 0; i < n; i++) {
        if (Math.abs(b[i] - x[i]) > tol + 1e-4 * Math.abs(b[i])) converged = false;
        x[i] += alpha * (b[i] - x[i]);
      }
      if (alpha < 1) converged = false;
      // junctions limit their own voltages inside stamp(); after this we just test
      this.stats.iters++;
      if (converged && it > 0) return true;
      if (this.nonlinear.length === 0) return true;
    }
    this.stats.retries = (this.stats.retries || 0) + 1;
    return false;
  }

  // Advance one sample. Returns true if Newton converged.
  step() {
    this.stats.steps++;
    return this._step(this.h, 0);
  }

  // If Newton fails (e.g. a comparator flipping), retry as two half steps.
  _step(h, depth) {
    const save = (this._save ||= []).length > depth ? this._save[depth] : (this._save[depth] = new Float64Array(this.n));
    save.set(this.x);
    if (this.newton()) {
      for (const el of this.dynamic) el.accept(this.x);
      return true;
    }
    if (depth >= 6) { this.stats.failures++; for (const el of this.dynamic) el.accept(this.x); return false; }
    this.x.set(save);
    this.setStep(h / 2);
    const ok = this._step(h / 2, depth + 1) && this._step(h / 2, depth + 1);
    this.setStep(h);
    return ok;
  }
}

// ---------- stamping helpers ----------
function stampG(A, n, a, b, g) {
  if (a >= 0) A[a * n + a] += g;
  if (b >= 0) A[b * n + b] += g;
  if (a >= 0 && b >= 0) { A[a * n + b] -= g; A[b * n + a] -= g; }
}
function stampI(bv, a, b, i) { // current i flowing from a to b through the element
  if (a >= 0) bv[a] -= i;
  if (b >= 0) bv[b] += i;
}
const V = (x, i) => (i < 0 ? 0 : x[i]);

class Resistor {
  constructor(a, b, r) { this.a = a; this.b = b; this.g = 1 / r; }
  stampA(A, n) { stampG(A, n, this.a, this.b, this.g); }
  stampB() {}
}

class Capacitor {
  constructor(a, b, c, v0) { Object.assign(this, { a, b, c, v: v0, i: 0, dynamic: true }); }
  setStep(h) { this.geq = 2 * this.c / h; }
  // trapezoidal companion: i_new = geq*v_new - (geq*v_old + i_old)
  stampA(A, n) { stampG(A, n, this.a, this.b, this.dcMode ? GMIN : this.geq); }
  stampB(b) { if (!this.dcMode) stampI(b, this.a, this.b, -(this.geq * this.v + this.i)); }
  accept(x) {
    const vn = V(x, this.a) - V(x, this.b);
    this.i = this.dcMode ? 0 : this.geq * vn - (this.geq * this.v + this.i);
    this.v = vn;
  }
}

class VSource {
  constructor(p, nn, volts, br) { Object.assign(this, { p, nn, volts, br }); }
  stampA(A, n) {
    const r = this.row;
    if (this.p >= 0) { A[this.p * n + r] += 1; A[r * n + this.p] += 1; }
    if (this.nn >= 0) { A[this.nn * n + r] -= 1; A[r * n + this.nn] -= 1; }
  }
  stampB(b) { b[this.row] += this.volts; }
}

class Diode {
  constructor(a, k, p = {}) {
    this.a = a; this.k = k; this.nonlinear = true;
    this.is = p.is ?? 2.52e-9; this.nvt = (p.n ?? 1.752) * VT; // 1N4148-ish
    this.vcrit = this.nvt * Math.log(this.nvt / (Math.SQRT2 * this.is));
    this.vd = 0;
  }
  stamp(A, b, x, n) {
    const vd = pnjlim(V(x, this.a) - V(x, this.k), this.vd, this.nvt, this.vcrit);
    this.vd = vd;
    const e = Math.exp(vd / this.nvt);
    const id = this.is * (e - 1), gd = this.is * e / this.nvt + GMIN;
    stampG(A, n, this.a, this.k, gd);
    stampI(b, this.a, this.k, id - gd * vd);
  }
}

function stampTerminal(A, bv, n, node, b, e, c, pol, i0, gbe, gbc, vbe, vbc) {
  if (node < 0) return;
  const row = node * n;
  if (b >= 0) A[row + b] += gbe + gbc;
  if (e >= 0) A[row + e] -= gbe;
  if (c >= 0) A[row + c] -= gbc;
  bv[node] -= pol * (i0 - gbe * vbe - gbc * vbc);
}

// Ebers-Moll transport model. pol = +1 NPN, -1 PNP. Defaults ~ 2N3904.
class BJT {
  constructor(c, b, e, pol, p = {}) {
    Object.assign(this, { c, b, e, pol, nonlinear: true });
    this.is = p.is ?? 6.734e-15; this.bf = p.bf ?? 416; this.br = p.br ?? 0.7371;
    this.vt = VT;
    this.vcrit = VT * Math.log(VT / (Math.SQRT2 * this.is));
    this.vbe = 0; this.vbc = 0;
  }
  stamp(A, bv, x, n) {
    const { c, b, e, pol, is, bf, br, vt } = this;
    let vbe = pol * (V(x, b) - V(x, e));
    let vbc = pol * (V(x, b) - V(x, c));
    vbe = this.vbe = pnjlim(vbe, this.vbe, vt, this.vcrit);
    vbc = this.vbc = pnjlim(vbc, this.vbc, vt, this.vcrit);
    const ef = Math.exp(vbe / vt), er = Math.exp(vbc / vt);
    const it = is * (ef - er);                    // transport current C->E
    const ibf = is / bf * (ef - 1), ibr = is / br * (er - 1);
    const ic = it - ibr, ib = ibf + ibr;
    const gf = is * ef / vt, gr = is * er / vt;   // d(exp terms)/dv
    // partial derivatives (in the pol-normalised frame)
    const dic_dvbe = gf, dic_dvbc = -gr - gr / br;
    const dib_dvbe = gf / bf + GMIN, dib_dvbc = gr / br + GMIN;
    // Currents into terminals (real frame): Ic into c, Ib into b, Ie = -(Ic+Ib) into e.
    // Linearise: I = I0 + dI/dvbe*(vbe-vbe0) + dI/dvbc*(vbc-vbc0), vbe = pol*(Vb-Ve), vbc = pol*(Vb-Vc)
    // Terminal currents (into the device, i.e. leaving the node) in the real frame:
    //   pol * (i0 + gbe*(vbe - vbe0) + gbc*(vbc - vbc0)),  vbe = pol*(Vb-Ve), vbc = pol*(Vb-Vc)
    // pol*pol = 1, so the Jacobian entries don't depend on polarity.
    stampTerminal(A, bv, n, c, b, e, c, pol, ic, dic_dvbe, dic_dvbc, vbe, vbc);
    stampTerminal(A, bv, n, b, b, e, c, pol, ib, dib_dvbe, dib_dvbc, vbe, vbc);
    stampTerminal(A, bv, n, e, b, e, c, pol, -(ic + ib), -(dic_dvbe + dib_dvbe), -(dic_dvbc + dib_dvbc), vbe, vbc);
  }
}

// Op-amp: ideal inputs, output is a voltage source = vsat*tanh(gain*vd/vsat)
// behind a small output resistance. Gain-bandwidth is not modelled (yet).
class OpAmp {
  constructor(inp, inn, out, br, p = {}) {
    Object.assign(this, { inp, inn, out, br, nonlinear: true });
    this.gain = p.gain ?? 2e4; this.vsat = p.vsat ?? 10.5; this.vos = p.vos ?? 0; // input offset voltage
  }
  stamp(A, b, x, n) {
    const r = this.row, { inp, inn, out, gain, vsat } = this;
    const vd = V(x, inp) - V(x, inn) + this.vos;
    const t = Math.tanh(gain * vd / vsat);
    const f = vsat * t, df = gain * (1 - t * t);
    // branch current flows out of the output node into the circuit
    if (out >= 0) { A[out * n + r] += 1; A[r * n + out] += 1; }
    // constraint: Vout - f(vd) = 0  ->  Vout - df*vd = f - df*vd0
    if (inp >= 0) A[r * n + inp] -= df;
    if (inn >= 0) A[r * n + inn] += df;
    b[r] += f - df * (vd - this.vos);
  }
}

export { VT };
