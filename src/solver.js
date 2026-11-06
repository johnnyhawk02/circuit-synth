// Component-level circuit solver.
//
// Netlist -> Modified Nodal Analysis (MNA). Capacitors use trapezoidal (or backward-Euler)
// companion models. Every nonlinear part is written as one or more "ports": a port voltage
// that is a linear combination of node voltages, and a port current f(v) that is injected
// back into the nodes through a fixed pattern. Because the rest of the circuit is linear,
// we precompute G^-1 once (the "DK method" / nodal K-method) and Newton only iterates on
// the small port vector each sample:
//
//   v = P·e + Q·f(v)       e = sources + capacitor history,  P, Q precomputed
//
// Each sample costs O(m^3) in the number of ports m (2 per full BJT, 1 per diode,
// forward-only BJT, JFET or op-amp) instead of O(n^3) in the node count. When a part
// value changes (a pot, a switch) the matrices are rebuilt once, on the next step.
//
// Works in Node (ES modules) and in the browser / AudioWorklet.

const VT = 0.025852; // thermal voltage at 27 °C
const GMIN = 1e-12;  // tiny conductance across every junction, helps convergence

// ---------- dense LU (used only when the linear part changes) ----------
function luFactor(A, n, piv) {
  for (let k = 0; k < n; k++) {
    let p = k, max = Math.abs(A[k * n + k]);
    for (let i = k + 1; i < n; i++) { const v = Math.abs(A[i * n + k]); if (v > max) { max = v; p = i; } }
    if (max < 1e-300) throw new Error('Singular matrix (floating node?) at row ' + k);
    piv[k] = p;
    if (p !== k) for (let j = 0; j < n; j++) { const t = A[k * n + j]; A[k * n + j] = A[p * n + j]; A[p * n + j] = t; }
    const d = A[k * n + k];
    for (let i = k + 1; i < n; i++) {
      const f = (A[i * n + k] /= d);
      if (f !== 0) for (let j = k + 1; j < n; j++) A[i * n + j] -= f * A[k * n + j];
    }
  }
}
function luSolve(A, n, piv, b) {
  for (let k = 0; k < n; k++) { const p = piv[k]; if (p !== k) { const t = b[k]; b[k] = b[p]; b[p] = t; } }
  for (let i = 1; i < n; i++) { let s = b[i]; for (let j = 0; j < i; j++) s -= A[i * n + j] * b[j]; b[i] = s; }
  for (let i = n - 1; i >= 0; i--) { let s = b[i]; for (let j = i + 1; j < n; j++) s -= A[i * n + j] * b[j]; b[i] = s / A[i * n + i]; }
}

// small in-place Gaussian elimination for the per-sample port Newton step
function smallSolve(J, r, m) {
  for (let k = 0; k < m; k++) {
    let p = k, max = Math.abs(J[k * m + k]);
    for (let i = k + 1; i < m; i++) { const v = Math.abs(J[i * m + k]); if (v > max) { max = v; p = i; } }
    if (max < 1e-300) return false;
    if (p !== k) {
      for (let j = k; j < m; j++) { const t = J[k * m + j]; J[k * m + j] = J[p * m + j]; J[p * m + j] = t; }
      const t = r[k]; r[k] = r[p]; r[p] = t;
    }
    const d = J[k * m + k];
    for (let i = k + 1; i < m; i++) {
      const f = J[i * m + k] / d;
      if (f === 0) continue;
      for (let j = k + 1; j < m; j++) J[i * m + j] -= f * J[k * m + j];
      r[i] -= f * r[k];
    }
  }
  for (let i = m - 1; i >= 0; i--) {
    let s = r[i];
    for (let j = i + 1; j < m; j++) s -= J[i * m + j] * r[j];
    r[i] = s / J[i * m + i];
  }
  return true;
}

function smallLU(A, m, piv) {
  for (let k = 0; k < m; k++) {
    let p = k, max = Math.abs(A[k * m + k]);
    for (let i = k + 1; i < m; i++) { const v = Math.abs(A[i * m + k]); if (v > max) { max = v; p = i; } }
    if (!(max > 1e-300)) return false;
    piv[k] = p;
    if (p !== k) for (let j = 0; j < m; j++) { const t = A[k * m + j]; A[k * m + j] = A[p * m + j]; A[p * m + j] = t; }
    const d = A[k * m + k];
    for (let i = k + 1; i < m; i++) {
      const f = (A[i * m + k] /= d);
      if (f !== 0) for (let j = k + 1; j < m; j++) A[i * m + j] -= f * A[k * m + j];
    }
  }
  return true;
}
function smallLUSolve(A, m, piv, b) {
  for (let k = 0; k < m; k++) { const p = piv[k]; if (p !== k) { const t = b[k]; b[k] = b[p]; b[p] = t; } }
  for (let i = 1; i < m; i++) { let s = b[i]; const o = i * m; for (let j = 0; j < i; j++) s -= A[o + j] * b[j]; b[i] = s; }
  for (let i = m - 1; i >= 0; i--) { let s = b[i]; const o = i * m; for (let j = i + 1; j < m; j++) s -= A[o + j] * b[j]; b[i] = s / A[o + i]; }
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

// ---------- port kinds ----------
// Each port: { kind, nv: [[unknownIndex, coef]...] (port voltage), inj: [[row, coef]...] (where f goes),
//              params }. eval(v) -> sets port.f and port.d (df/dv).
const EXP = 0, TANH = 1, JFET = 2, JFET3 = 3, ZERO = 4, CLAMP = 5;

// Shichman-Hodges N-JFET drain current for vds >= 0, softplus-smoothed overdrive.
// returns [id, d id/d vgs, d id/d vds]
function jfetId(beta, vto, soft, vgs, vds) {
  const x = (vgs - vto) / soft;
  const vov = x > 30 ? soft * x : soft * Math.log1p(Math.exp(x));
  const dvov = x > 30 ? 1 : 1 / (1 + Math.exp(-x));
  if (vds < vov) return [beta * (2 * vov - vds) * vds, 2 * beta * vds * dvov, 2 * beta * (vov - vds)];
  return [beta * vov * vov, 2 * beta * vov * dvov, 0];
}

function evalPort(p, v) {
  switch (p.kind) {
    case EXP: {
      const e = Math.exp(v / p.nvt);
      p.f = p.is * (e - 1) + GMIN * v;
      p.d = p.is * e / p.nvt + GMIN;
      break;
    }
    case TANH: {
      const t = Math.tanh(p.k * (v + p.vos));
      p.f = p.vsat * t;
      p.d = p.vsat * p.k * (1 - t * t) + 1e-9;
      break;
    }
    case JFET3: { // port voltage = vds; partner port carries vgs
      const vgs = p.partner.v, vds = v;
      let r;
      if (vds >= 0) r = jfetId(p.beta, p.vto, p.soft, vgs, vds);
      else { const q = jfetId(p.beta, p.vto, p.soft, vgs - vds, -vds); r = [-q[0], -q[1], q[1] + q[2]]; }
      p.f = r[0] + GMIN * vds; p.d = r[2] + GMIN; p.xd = r[1];
      break;
    }
    case ZERO: p.f = 0; p.d = 0; break;
    case CLAMP: { // soft rail clamp: conducts only beyond ±vsat (smooth hinge, slope k siemens)
      const s = p.soft, a = (v - p.vsat) / s, b = (-v - p.vsat) / s;
      const spa = a > 30 ? a : Math.log1p(Math.exp(a)), spb = b > 30 ? b : Math.log1p(Math.exp(b));
      const sa = a > 30 ? 1 : 1 / (1 + Math.exp(-a)), sb = b > 30 ? 1 : 1 / (1 + Math.exp(-b));
      p.f = p.k * s * (spa - spb); p.d = p.k * (sa + sb) + GMIN;
      break;
    }
    case JFET: { // square law in saturation, softplus-smoothed at pinch-off
      const x = (v - p.vto) / p.soft;
      const sp = x > 30 ? x : Math.log1p(Math.exp(x));
      const sig = x > 30 ? 1 : 1 / (1 + Math.exp(-x));
      const u = p.soft * sp;
      p.f = p.beta * u * u + GMIN * v;
      p.d = 2 * p.beta * u * sig + GMIN;
      break;
    }
  }
}

// all ports share one shape so the hot loop stays monomorphic
class Port {
  constructor(p) {
    this.kind = p.kind; this.is = p.is ?? 0; this.nvt = p.nvt ?? 1; this.vcrit = p.vcrit ?? 0;
    this.vsat = p.vsat ?? 0; this.k = p.k ?? 0; this.vos = p.vos ?? 0; this.tanhStep = p.tanhStep ?? 1;
    this.beta = p.beta ?? 0; this.vto = p.vto ?? 0; this.soft = p.soft ?? 1; this.partner = null;
    this.nv = p.nv; this.inj = p.inj; this.v = 0; this.f = 0; this.d = 0; this.xd = 0;
  }
}

// ---------- the circuit ----------
export class Circuit {
  constructor({ integrator = 'trap' } = {}) {
    this.nodeIndex = new Map([['0', -1], ['gnd', -1]]);
    this.nodeCount = 0;
    this.branchCount = 0;
    this.resistors = [];
    this.caps = [];
    this.sources = [];
    this.devices = [];   // nonlinear parts, each owning ports
    this.inputs = new Map();
    this.integrator = integrator;
    this.probes = [];
    // every field the hot path touches exists from the start, in one order, so all circuits
    // share a hidden class and property loads stay monomorphic
    this.opamps = []; this._oaN = 0; this.watches = null; this._st = null; this.dirty = true;
    this.n = 0; this.m = 0; this.ns = 0; this.nc = 0; this.k = 0; this.ports = null; this.cross = null;
    this.e = null; this.vp = null; this.fp = null; this.pl = null; this.J = null; this.r = null;
    this.vSave = null; this.piv = null; this.dd = null; this.dv = null; this.vPrev = null; this.tolv = null;
    this.predict = false; this.luOK = false; this.stats = null; this.eLast = null;
    this.levels = null; this.dcLevel = null; this.cur = null; this.level = 0;
    this.X = null; this.cols = 0; this.P = null; this.Q = null; this.CP = null;
    this.hold = 0; this.startLevel = 0; this.deepest = 0; this.chord = false; this.hasPrev = false;
    this.h = 0; this.maxLevel = 8;
  }

  node(name) {
    if (!this.nodeIndex.has(name)) this.nodeIndex.set(name, this.nodeCount++);
    return this.nodeIndex.get(name);
  }
  branch() { return this.branchCount++; }

  // --- component helpers (nodes are names, '0' is ground) ---
  R(a, b, ohms) {
    const r = { a: this.node(a), b: this.node(b), g: 1 / ohms, circuit: this };
    // r.set(ohms) changes the value later (pots, switches); the matrices rebuild on the next step
    r.set = ohms => { r.g = 1 / ohms; this.dirty = true; };
    this.resistors.push(r);
    return r;
  }
  // Potentiometer a—(wiper)—b, total ohms, pos 0 (at a) .. 1 (at b). Returns { set(pos) }.
  POT(a, b, wiper, ohms, pos = 0.5) {
    const ra = this.R(a, wiper, 1), rb = this.R(wiper, b, 1);
    const pot = { pos, set: p => { pot.pos = p; ra.set(Math.max(ohms * p, 0.01)); rb.set(Math.max(ohms * (1 - p), 0.01)); } };
    pot.set(pos);
    return pot;
  }
  // { be: true } integrates this capacitor with backward Euler (for very fast parasitic /
  // internal capacitors, where trapezoidal would ring at Nyquist)
  C(a, b, farads, v0 = 0, opts = {}) {
    const c = { a: this.node(a), b: this.node(b), c: farads, v: v0, i: 0, ieq: 0, geq: 0,
      be: opts.be ?? this.integrator === 'be' };
    // c.set(farads) swaps the part (e.g. a range switch); charge is kept, current reset
    c.set = f => { c.c = f; c.i = 0; this.dirty = true; };
    this.caps.push(c);
    return c;
  }
  V(p, n, volts, name) {
    const s = { p: this.node(p), n: this.node(n), volts, br: this.branch() };
    this.sources.push(s);
    if (name) this.inputs.set(name, s);
    return s;
  }
  D(a, k, p = {}) {
    const A = this.node(a), K = this.node(k);
    const nvt = (p.n ?? 1.752) * VT, is = p.is ?? 2.52e-9; // 1N4148-ish
    const port = { kind: EXP, is, nvt, vcrit: nvt * Math.log(nvt / (Math.SQRT2 * is)),
      nv: [[A, 1], [K, -1]], inj: [[A, -1], [K, 1]] };
    return this._dev([port]);
  }
  // Bipolar transistors. Defaults ~ 2N3904 / 2N3906 class small-signal parts.
  // { simple: true } drops the base-collector junction (forward-active only): half the ports,
  // fine for any transistor that never saturates (ladder stages, followers, diff pairs).
  NPN(c, b, e, p) { return this._bjt(c, b, e, +1, p); }
  PNP(c, b, e, p) { return this._bjt(c, b, e, -1, p); }
  _bjt(c, b, e, pol, p = {}) {
    const C = this.node(c), B = this.node(b), E = this.node(e);
    const is = p.is ?? 6.734e-15, bf = p.bf ?? 416, br = p.br ?? 0.7371;
    const vcrit = VT * Math.log(VT / (Math.SQRT2 * is));
    // forward junction: transport current flows C->E, base current is if/bf
    const fwd = { kind: EXP, is, nvt: VT, vcrit,
      nv: [[B, pol], [E, -pol]],
      inj: [[C, -pol], [B, -pol / bf], [E, pol * (1 + 1 / bf)]] };
    // optional junction capacitances (b-e, b-c): give regenerative switching circuits
    // (flip-flops, Schmitt triggers) real dynamics so each sample has a unique solution
    if (p.cj) { this.C(b, e, p.cj, 0, { be: true }); this.C(b, c, p.cj / 3, 0, { be: true }); }
    if (p.simple) return this._dev([fwd]);
    const rev = { kind: EXP, is, nvt: VT, vcrit,
      nv: [[B, pol], [C, -pol]],
      inj: [[C, pol * (1 + 1 / br)], [B, -pol / br], [E, -pol]] };
    return this._dev([fwd, rev]);
  }
  // N-channel JFET in saturation (drain current depends on Vgs only). ~2N4303 defaults.
  NJF(d, g, s, p = {}) {
    const D = this.node(d), G = this.node(g), S = this.node(s);
    const port = { kind: JFET, beta: p.beta ?? 0.6e-3, vto: p.vto ?? -3, soft: 0.02,
      nv: [[G, 1], [S, -1]], inj: [[D, -1], [S, 1]] };
    return this._dev([port]);
  }
  // N-channel JFET with triode region (use it as an analog switch). ~J111/2N4391 defaults.
  NJFS(d, g, s, p = {}) {
    const D = this.node(d), G = this.node(g), S = this.node(s);
    const pg = { kind: ZERO, nv: [[G, 1], [S, -1]], inj: [] };
    const pd = { kind: JFET3, beta: p.beta ?? 1.5e-3, vto: p.vto ?? -3, soft: 0.02, partner: pg,
      nv: [[D, 1], [S, -1]], inj: [[D, -1], [S, 1]] };
    return this._dev([pd, pg]);
  }
  // Op-amp macro-model: out = vsat * tanh(gain * (in+ - in- + vos) / vsat), ideal inputs,
  // zero output impedance. Stands in for both IC op-amps and Moog's discrete Q1-Q5 adders.
  //
  // { gbw: Hz, slew: V/s } adds a dominant pole and slew limit (an internal integrator node,
  // backward Euler): use it for op-amps used as comparators, whose regenerative flip has no
  // continuous solution with an infinitely fast op-amp.
  OPAMP(inp, inn, out, p = {}) {
    const P = this.node(inp), N = this.node(inn), O = this.node(out);
    const br = this.branch();
    const vsat = p.vsat ?? 10.5, gain = p.gain ?? 2e4;
    let ports;
    if (p.gbw) {
      // internal node Y integrates gm*vd (slew-limited), is clamped at the rails, and a
      // unity buffer drives the output
      const id = (this._oaN = (this._oaN || 0) + 1), y = '__oa' + id;
      const Y = this.node(y), wt = 2 * Math.PI * p.gbw, sr = p.slew ?? 0.5e6, Ci = 1e-9;
      this.C(y, '0', Ci, 0, { be: true });
      this.R(y, '0', gain / (wt * Ci));
      const gm = { kind: TANH, vsat: Ci * sr, k: wt / sr, vos: p.vos ?? 0,
        nv: [[P, 1], [N, -1]], inj: [[Y, 1]], tanhStep: 1e9 };
      // internal node clamped just beyond the rails (stands in for the output stage limits)
      const clamp = { kind: CLAMP, vsat: vsat, k: 1, soft: 0.05, nv: [[Y, 1]], inj: [[Y, -1]] };
      ports = [gm, clamp];
      const dev = this._dev(ports);
      dev.out = O; dev.br = br; dev.y = Y;
      this.opamps.push(dev);
      return dev;
    } else {
      ports = [{ kind: TANH, vsat, k: gain / vsat, vos: p.vos ?? 0,
        nv: [[P, 1], [N, -1]], inj: [[{ br }, 1]], tanhStep: p.step ?? 1 }];
    }
    const dev = this._dev(ports);
    dev.out = O; dev.br = br;
    this.opamps.push(dev);
    return dev;
  }
  _dev(ports) { const d = { ports }; this.devices.push(d); return d; }

  set(name, volts) { this.inputs.get(name).volts = volts; }
  probe(...names) { for (const n of names) if (!this.probes.includes(n)) this.probes.push(n); return this; }

  // ---------- compile ----------
  finalize(sampleRate) {
    this.n = this.nodeCount + this.branchCount;
    const row = x => (typeof x === 'object' ? this.nodeCount + x.br : x);
    for (const s of this.sources) s.row = this.nodeCount + s.br;
    this.ports = [];
    const map = new Map();
    for (const d of this.devices) for (const raw of d.ports) {
      const p = new Port(raw);
      p.inj = raw.inj.map(([r, c]) => [row(r), c]).filter(([r]) => r >= 0);
      p.nv = raw.nv.filter(([r]) => r >= 0);
      map.set(raw, p);
      this.ports.push(p);
    }
    for (const d of this.devices) for (const raw of d.ports) if (raw.partner) map.get(raw).partner = map.get(raw.partner);
    this.m = this.ports.length;
    this.cross = [];
    this.ports.forEach((p, i) => { if (p.partner) this.cross.push({ i, x: this.ports.indexOf(p.partner) }); });
    this.ns = this.sources.length;
    this.nc = this.caps.length;
    this.k = this.ns + this.nc;
    this.e = new Float64Array(this.k);
    this.vp = new Float64Array(this.m);   // port voltages
    this.fp = new Float64Array(this.m);   // port currents
    this.pl = new Float64Array(this.m);   // linear part of port voltages this step
    this.J = new Float64Array(this.m * this.m);
    this.r = new Float64Array(this.m);
    this.vSave = new Float64Array(this.m);
    this.piv = new Int32Array(this.m);
    this.dd = new Float64Array(this.m);
    this.dv = new Float64Array(this.m);
    this.vPrev = new Float64Array(this.m);
    // convergence tolerance per port, in port volts: 1 µV on junctions; for op-amp inputs
    // scaled by the gain so the output is good to ~1 µV
    this.tolv = Float64Array.from(this.ports, p => p.kind === TANH ? 1e-4 / Math.max(1, p.vsat * p.k) : 1e-5);
    this.predict = true;   // linear extrapolation of the port voltages (about halves the iterations)
    this.luOK = false;
    this.stats = { steps: 0, iters: 0, failures: 0, rebuilds: 0, substeps: 0, factors: 0, fast: 0 };
    this.eLast = new Float64Array(this.k);
    this.levels = []; this.dcLevel = null; this.cur = null; this.level = 0;
    this.X = null; this.cols = 0; this.P = null; this.Q = null; this.CP = null;
    this.hold = 0; this.startLevel = 0; this.deepest = 0; this.chord = false; this.hasPrev = false;
    this.h = 1 / sampleRate;
    this.dirty = true;
    this.maxLevel = 8;   // smallest sub-step h/256
    return this;
  }

  setSampleRate(sr) { this.h = 1 / sr; this.dirty = true; }
  invalidate() { this.dirty = true; } // kept for older callers that poke r.g directly

  // Build G for step size h (h = 0: DC, caps open), factor it, and precompute P/Q (ports),
  // capacitor and probe readouts. Returns a "level" object; levels are cached per step size.
  _build(h) {
    const { n, m, k, ns, nc } = this;
    const G = new Float64Array(n * n);
    const add = (i, j, v) => { if (i >= 0 && j >= 0) G[i * n + j] += v; };
    const stampG = (a, b, g) => { add(a, a, g); add(b, b, g); add(a, b, -g); add(b, a, -g); };
    for (const r of this.resistors) stampG(r.a, r.b, r.g);
    const geq = new Float64Array(nc);
    this.caps.forEach((c, i) => {
      geq[i] = h === 0 ? GMIN : (c.be ? c.c / h : 2 * c.c / h);
      stampG(c.a, c.b, geq[i]);
    });
    for (const s of this.sources) { add(s.p, s.row, 1); add(s.row, s.p, 1); add(s.n, s.row, -1); add(s.row, s.n, -1); }
    for (const d of this.opamps || []) {
      const r = this.nodeCount + d.br; add(d.out, r, 1); add(r, d.out, 1);
      if (d.y !== undefined) add(r, d.y, -1); // buffered op-amp: out - y = 0
    }
    // tiny leak on every node so a node that only touches caps/junctions is never singular
    for (let i = 0; i < this.nodeCount; i++) G[i * n + i] += GMIN;
    const piv = new Int32Array(n);
    luFactor(G, n, piv);
    // columns: excitation (k) then port injections (m); solve G X = [B M]
    const cols = k + m;
    const X = new Float64Array(n * cols);
    const col = new Float64Array(n);
    const solveCol = (j, entries) => {
      col.fill(0);
      for (const [r, v] of entries) if (r >= 0) col[r] += v;
      luSolve(G, n, piv, col);
      for (let i = 0; i < n; i++) X[i * cols + j] = col[i];
    };
    this.sources.forEach((s, j) => solveCol(j, [[s.row, 1]]));
    this.caps.forEach((c, j) => solveCol(ns + j, [[c.a, 1], [c.b, -1]]));
    this.ports.forEach((p, j) => solveCol(k + j, p.inj));
    const rowOf = (entries, out, from, len) => {
      for (let j = 0; j < len; j++) {
        let s = 0;
        for (const [r, c] of entries) if (r >= 0) s += c * X[r * cols + from + j];
        out[j] = s;
      }
    };
    const P = new Float64Array(m * k), Q = new Float64Array(m * m);
    this.ports.forEach((p, i) => {
      rowOf(p.nv, P.subarray(i * k, i * k + k), 0, k);
      rowOf(p.nv, Q.subarray(i * m, i * m + m), k, m);
    });
    const CP = new Float64Array(nc * cols);
    this.caps.forEach((c, i) => rowOf([[c.a, 1], [c.b, -1]], CP.subarray(i * cols, (i + 1) * cols), 0, cols));
    this.stats.rebuilds++;
    return { h, X, cols, P, Q, CP, geq };
  }

  // Switch to the matrices for step size h / 2^level (0 = full step, -1 = DC).
  _use(level) {
    if (this.dirty) { this.levels = []; this.dcLevel = null; this.dirty = false; this.cur = null; }
    let L;
    if (level < 0) L = this.dcLevel ||= this._build(0);
    else L = this.levels[level] ||= this._build(this.h / (1 << level));
    if (L === this.cur) return;
    this.X = L.X; this.cols = L.cols; this.P = L.P; this.Q = L.Q; this.CP = L.CP;
    // geq changed: re-derive each capacitor's history source from its state (v, i)
    for (let i = 0; i < this.nc; i++) {
      const c = this.caps[i];
      c.geq = L.geq[i];
      c.ieq = level < 0 ? 0 : (c.be ? c.geq * c.v : c.geq * c.v + c.i);
    }
    this.cur = L;
    this.level = level;
    this.luOK = false;
  }
  rebuild() { this.dirty = true; this._use(this.level); }

  _loadE() {
    const e = this.e;
    for (let j = 0; j < this.ns; j++) e[j] = this.sources[j].volts;
    for (let j = 0; j < this.nc; j++) e[this.ns + j] = this.caps[j].ieq;
  }

  // Newton on the port vector. Returns true if converged.
  // circuit.chord = true reuses the factored Jacobian across iterations and samples and only
  // refactors when convergence slows down (off by default: measured slower on these circuits).
  newton(maxIter = 60) {
    const { m, k, P, Q, J, r, vp, fp, pl, e, ports, cross, piv, dd, tolv, dv } = this;
    for (let i = 0; i < m; i++) {
      let s = 0; const o = i * k;
      for (let j = 0; j < k; j++) s += P[o + j] * e[j];
      pl[i] = s;
    }
    if (m === 0) return true;
    let prevErr = Infinity, refactor = !this.luOK || !this.chord;
    // (non-chord: refactor at the first iteration, then reuse for small corrections)
    for (let it = 0; it < maxIter; it++) {
      for (let i = 0; i < m; i++) ports[i].v = vp[i];
      for (let i = 0; i < m; i++) { const p = ports[i]; evalPort(p, vp[i]); fp[i] = p.f; dd[i] = p.d; }
      // residual F = v - pl - Q f
      for (let i = 0; i < m; i++) {
        let s = 0; const o = i * m;
        for (let j = 0; j < m; j++) s += Q[o + j] * fp[j];
        r[i] = -(vp[i] - pl[i] - s);
      }
      this.stats.iters++;
      let fresh = false;
      if (refactor) {
        // Jacobian I - Q D (D diagonal plus a few cross terms), LU-factored in place
        for (let i = 0; i < m; i++) {
          const o = i * m;
          for (let j = 0; j < m; j++) J[o + j] = -Q[o + j] * dd[j];
          for (let c = 0; c < cross.length; c++) J[o + cross[c].x] -= Q[o + cross[c].i] * ports[cross[c].i].xd;
          J[o + i] += 1;
        }
        if (!smallLU(J, m, piv)) { this.luOK = false; break; }
        this.luOK = true; fresh = true; refactor = false;
        this.stats.factors++;
      }
      smallLUSolve(J, m, piv, r);
      let converged = true, err = 0;
      for (let i = 0; i < m; i++) {
        const p = ports[i], kind = p.kind, step = r[i], old = vp[i];
        let vn = old + step;
        if (kind === EXP) vn = pnjlim(vn, old, p.nvt, p.vcrit);
        else if (kind === TANH) { const lim = p.tanhStep; if (step > lim) vn = old + lim; else if (step < -lim) vn = old - lim; }
        else { const lim = kind === CLAMP ? 0.5 : 1; if (step > lim) vn = old + lim; else if (step < -lim) vn = old - lim; }
        const d = vn - old, ad = Math.abs(d);
        if (!(ad <= tolv[i])) converged = false;
        if (ad > err) err = ad;
        dv[i] = d;
        vp[i] = vn;
      }
      // a node fed only by transistor currents (a diff-pair tail) makes the port equations
      // ill-conditioned: round-off then sets a noise floor ~1e-7 V. Stagnation below 2 µV
      // is accepted as converged.
      if (!converged && err < 2e-6 && err > 0.25 * prevErr) converged = true;
      // slow contraction with a stale Jacobian: refactor on the next iteration
      // the Jacobian is factored on the first iteration of every sample; later iterations
      // (usually just the confirming one) reuse it unless convergence is slow
      if (!converged && (it > 0 ? err > 0.1 * prevErr : (this.chord ? err > 1e-3 : err > 1e-2))) refactor = true;
      prevErr = err;
      if (converged) {
        // the last update was tiny: move the port currents along their slopes instead of
        // re-evaluating every exp/tanh
        for (let i = 0; i < m; i++) { const p = ports[i]; p.v = vp[i]; fp[i] += dd[i] * dv[i]; p.f = fp[i]; }
        for (let c = 0; c < cross.length; c++) { const ci = cross[c].i; fp[ci] += ports[ci].xd * dv[cross[c].x]; ports[ci].f = fp[ci]; }
        return true;
      }
    }
    for (let i = 0; i < m; i++) ports[i].v = vp[i];
    for (let i = 0; i < m; i++) { evalPort(ports[i], vp[i]); fp[i] = ports[i].f; }
    this.luOK = false;
    return false;
  }

  // value of row vector R (length k+m) against current excitation and port currents
  _read(R, o, e) {
    const k = this.k, m = this.m, fp = this.fp;
    let s = 0;
    for (let j = 0; j < k; j++) s += R[o + j] * e[j];
    o += k;
    for (let j = 0; j < m; j++) s += R[o + j] * fp[j];
    return s;
  }

  _acceptCaps(dc = false) {
    // keep the excitation this solve used, so v() after step() reads the same solution
    this.eLast.set(this.e);
    const cols = this.cols, CP = this.CP, e = this.e;
    for (let i = 0; i < this.nc; i++) {
      const c = this.caps[i];
      const v = this._read(CP, i * cols, e);
      c.i = dc ? 0 : c.geq * v - c.ieq;
      c.v = v;
      c.ieq = dc ? 0 : (c.be ? c.geq * c.v : c.geq * c.v + c.i);
    }
  }

  // Solve the DC operating point (caps open) by ramping supplies in.
  dcop() {
    this.dirty = true;
    this._use(-1);
    const full = this.sources.map(s => s.volts);
    for (let s = 1; s <= 20; s++) {
      this.sources.forEach((src, i) => (src.volts = full[i] * s / 20));
      this._loadE();
      this.newton(300);
    }
    this.sources.forEach((src, i) => (src.volts = full[i]));
    this._loadE();
    this.newton(300);
    this._acceptCaps(true);
    this._use(0);
    // re-express the solution with the transient matrices so v() reads it correctly
    this._loadE();
    this.eLast.set(this.e);
  }

  // Advance one sample. Returns true if Newton converged. A step that fails is retried as
  // two half steps, recursively (down to h/2^maxLevel); matrices for each level are cached.
  // Fast-transition watch: while node `name` is strictly between lo and hi (e.g. a comparator
  // output between its rails) the step is taken as 2^level sub-steps. Regenerative flips
  // then resolve at their real speed instead of being smeared (or frozen) by the coarse step.
  watch(name, lo, hi, level = 4) {
    (this.watches ||= []).push({ name, lo, hi, level, idx: -1 });
    return this;
  }
  _transition() {
    let lvl = 0;
    for (const w of this.watches) {
      const v = this.vNode(w.idx);
      if (v > w.lo && v < w.hi && w.level > lvl) lvl = w.level;
    }
    return lvl;
  }
  _save() {
    const st = this._st ||= { vp: new Float64Array(this.m), cv: new Float64Array(this.nc), ci: new Float64Array(this.nc), e: new Float64Array(this.k) };
    st.vp.set(this.vp); st.e.set(this.eLast);
    for (let i = 0; i < this.nc; i++) { st.cv[i] = this.caps[i].v; st.ci[i] = this.caps[i].i; }
    st.level = this.level;
  }
  _restore() {
    const st = this._st;
    this.vp.set(st.vp); this.eLast.set(st.e);
    for (let i = 0; i < this.nc; i++) {
      const c = this.caps[i]; c.v = st.cv[i]; c.i = st.ci[i];
      c.ieq = c.be ? c.geq * c.v : c.geq * c.v + c.i;
    }
  }

  step() {
    if (this.watches) {
      if (this.watches[0].idx === -1) for (const w of this.watches) w.idx = this.nodeIndex.get(w.name);
      this.stats.steps++;
      let lvl = this._transition();
      if (lvl > 0) { this.stats.fast++; return this._advanceN(lvl); }
      this._save();
      const ok = this._advanceN(0);
      lvl = this._transition();
      if (lvl === 0) return ok;
      this._use(0); this._restore();
      this.stats.fast++;
      return this._advanceN(lvl);
    }
    return this._step();
  }
  _step() {
    this.stats.steps++;
    // after a step needed sub-steps, start the next few at that depth (no wasted retries)
    const start = this.hold > 0 ? this.startLevel : 0;
    if (this.hold > 0) this.hold--;
    this.deepest = 0;
    const ok = this._advanceN(start);
    if (this.deepest > 0) { this.startLevel = this.deepest; this.hold = 2; }
    return ok;
  }
  _advanceN(level) { // advance a full step using 2^level sub-steps
    let ok = true;
    for (let i = 0; i < (1 << level); i++) ok = this._advance(level) && ok;
    return ok;
  }
  _advance(level) {
    const switched = this.level !== level;
    this._use(level);
    this._loadE();
    this.vSave.set(this.vp);
    // optional linear extrapolation of the port voltages from the last two samples
    if (this.predict) {
      const vp = this.vp, vPrev = this.vPrev, m = this.m;
      if (level === 0 && !switched && this.hasPrev) {
        for (let i = 0; i < m; i++) { const x = vp[i]; vp[i] = 2 * x - vPrev[i]; vPrev[i] = x; }
      } else { vPrev.set(vp); this.hasPrev = level === 0; }
    }
    const last = level >= this.maxLevel;
    if (this.newton(last ? 100 : (level === 0 ? 15 : 25))) { this._acceptCaps(); return true; }
    if (last) { this.stats.failures++; this._acceptCaps(); return false; }
    this.vp.set(this.vSave);
    // a full step that fails is almost always a switching event: go straight to h/16
    const next = Math.min(this.maxLevel, level === 0 ? 4 : level + 1), n = 1 << (next - level);
    if (next > this.deepest) this.deepest = next;
    this.stats.substeps++;
    let ok = true;
    for (let i = 0; i < n; i++) ok = this._advance(next) && ok;
    return ok;
  }

  // Node voltage by name, or vNode(index) with an index from nodeIndex (faster in a loop).
  v(name) {
    const idx = this.nodeIndex.get(name);
    if (idx === undefined) throw new Error('unknown node ' + name);
    return this.vNode(idx);
  }
  vNode(idx) {
    if (idx < 0) return 0;
    if (this.dirty) this._use(this.level);
    return this._read(this.X, idx * this.cols, this.eLast);
  }
}

export { VT };
