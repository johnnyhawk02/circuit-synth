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

import { K, alloc, buffer, CTX } from './kernel.js';

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

// ---------- port kinds ----------
// Each port: { kind, nv: [[unknownIndex, coef]...] (port voltage), inj: [[row, coef]...] (where f goes),
//              params }. eval(v) -> sets port.f and port.d (df/dv).
// (the port equations themselves are evaluated in kernel/kernel.c)
const EXP = 0, TANH = 1, JFET = 2, JFET3 = 3, ZERO = 4, CLAMP = 5;

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
    this.opamps = []; this._oaN = 0; this.vars = []; this.cache = new Map(); this.watches = null; this._st = null; this.dirty = true;
    this.n = 0; this.m = 0; this.ns = 0; this.nc = 0; this.k = 0; this.ports = null; this.cross = null;
    this.e = null; this.vp = null; this.fp = null; this.pl = null; this.J = null; this.r = null;
    this.vSave = null; this.piv = null; this.dd = null; this.dv = null; this.vPrev = null; this.tolv = null;
    this.predict = false; this.luOK = false; this.stats = null; this.eLast = null;
    this.levels = null; this.dcLevel = null; this.cur = null; this.level = 0;
    this.X = null; this.cols = 0; this.P = null; this.Q = null; this.CP = null;
    this.hold = 0; this.startLevel = 0; this.deepest = 0; this.chord = false; this.hasPrev = false;
    this.h = 0; this.maxLevel = 8;
    this.Ja = null; this.act = null; this.pos = null; this.na = 0; this.ra = null; this.ddF = null;
    this.ctx = 0; this.off = null; this.buf = null; this.ci32 = null; this.f64n = null; this.i32n = null;
    this.wP = null; this.wQ = null; this.wCP = null; this.wqmax = null; this.capGeq = null; this.capIeq = null; this.capV = null; this.capI = null;
    this.isAct = null; this.xdF = null; this.keep = null; this.qmax = null; this.cfg = null;
  }

  node(name) {
    if (!this.nodeIndex.has(name)) this.nodeIndex.set(name, this.nodeCount++);
    return this.nodeIndex.get(name);
  }
  branch() { return this.branchCount++; }
  // a part whose value can change (pot, switch, range cap): its value is part of the key
  // that the cached matrices are filed under, so toggling back to an earlier setting
  // (a gate switch, say) reuses the matrices instead of rebuilding them
  _variable(x) {
    if (!x.isVar) { x.isVar = true; this.vars.push(x); }
    this.dirty = true;
  }

  // --- component helpers (nodes are names, '0' is ground) ---
  R(a, b, ohms) {
    const r = { a: this.node(a), b: this.node(b), g: 1 / ohms, circuit: this };
    // r.set(ohms) changes the value later (pots, switches); the matrices rebuild on the next step
    r.set = ohms => { const g = 1 / ohms; if (g === r.g) return; r.g = g; this._variable(r); };
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
    c.set = f => { if (f === c.c) return; c.c = f; c.i = 0; if (this.capI) this.capI[c.idx] = 0; this._variable(c); };
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
    this.vSave = new Float64Array(this.m);
    this.vPrev = new Float64Array(this.m);
    this._layout();
    this.predict = true;   // linear extrapolation of the port voltages (about halves the iterations)
    const self = this;
    this.stats = { steps: 0, failures: 0, rebuilds: 0, substeps: 0, fast: 0,
      get iters() { return self.ci32[self.ctx / 4 + 7]; }, get factors() { return self.ci32[self.ctx / 4 + 8]; } };
    this.levels = []; this.dcLevel = null; this.cur = null; this.level = 0;
    this.X = null; this.cols = 0; this.P = null; this.Q = null; this.CP = null;
    this.hold = 0; this.startLevel = 0; this.deepest = 0; this.chord = false; this.hasPrev = false;
    this.h = 1 / sampleRate;
    this.dirty = true;
    this.maxLevel = 8;   // smallest sub-step h/256
    return this;
  }

  // Lay this circuit's arrays out in the kernel's memory and fill the constant parts.
  _layout() {
    const { m, k, nc } = this, cols = k + m, nx = this.cross.length;
    const f64 = { P: m * k, Q: m * m, CP: nc * cols, qmax: m, e: k, eLast: k, vp: m, fp: m, pl: m, dd: m, dv: m,
      r: m, Ja: m * m, ra: m, ddF: m, tolv: m, xdF: nx, is: m, nvt: m, vcrit: m, vsat: m, kk: m, vos: m,
      tanhStep: m, beta: m, vto: m, soft: m, xd: m, capGeq: nc, capIeq: nc, capV: nc, capI: nc };
    const i32 = { act: m, pos: m, piv: m, keep: m, isAct: m, kind: m, partner: m, crossI: nx, crossX: nx, capBe: nc };
    this.ctx = alloc(CTX.length * 4);
    this.off = {};
    for (const [name, n] of Object.entries(f64)) this.off[name] = alloc(Math.max(1, n) * 8);
    for (const [name, n] of Object.entries(i32)) this.off[name] = alloc(Math.max(1, n) * 4);
    this.f64n = f64; this.i32n = i32;
    this._views();
    const c = this.ci32, b = this.ctx / 4;
    c[b] = m; c[b + 1] = k; c[b + 2] = nc; c[b + 3] = cols; c[b + 4] = nx;
    CTX.forEach((name, i) => { if (name in this.off) c[b + i] = this.off[name]; });
    const P = this.ports;
    P.forEach((p, i) => {
      this.kind[i] = p.kind; this.partner[i] = p.partner ? P.indexOf(p.partner) : 0;
      this.is[i] = p.is; this.nvt[i] = p.nvt; this.vcrit[i] = p.vcrit; this.vsat[i] = p.vsat; this.kk[i] = p.k;
      this.vos[i] = p.vos; this.tanhStep[i] = p.tanhStep; this.beta[i] = p.beta; this.vto[i] = p.vto; this.soft[i] = p.soft;
      // convergence tolerance per port, in port volts: 10 µV on junctions; for op-amp inputs
      // scaled by the gain so the output is good to ~0.1 mV
      this.tolv[i] = p.kind === TANH ? 1e-4 / Math.max(1, p.vsat * p.k) : 1e-5;
    });
    this.cross.forEach((x, j) => { this.crossI[j] = x.i; this.crossX[j] = x.x; this.keep[x.i] = 1; this.keep[x.x] = 1; });
    this.caps.forEach((cp, i) => { cp.idx = i; this.capV[i] = cp.v; this.capI[i] = cp.i; this.capBe[i] = cp.be ? 1 : 0; });
  }
  // (re)make the typed-array views; needed again whenever the kernel memory grows
  _views() {
    const buf = this.buf = buffer();
    this.ci32 = new Int32Array(buf);
    const W = { P: 'wP', Q: 'wQ', CP: 'wCP', qmax: 'wqmax' };   // the current level's matrices
    for (const [name, n] of Object.entries(this.f64n)) this[W[name] || name] = new Float64Array(buf, this.off[name], n);
    for (const [name, n] of Object.entries(this.i32n)) this[name] = new Int32Array(buf, this.off[name], n);
    this.cur = null;   // level matrices must be copied in again
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
    // largest |Q| in each column: a port with slope d moves the others by at most qmax·d
    const qmax = new Float64Array(m);
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) qmax[j] = Math.max(qmax[j], Math.abs(Q[i * m + j]));
    return { h, X, cols, P, Q, CP, geq, qmax };
  }

  // Switch to the matrices for step size h / 2^level (0 = full step, -1 = DC).
  _use(level) {
    if (this.buf !== buffer()) this._views();
    if (this.dirty) this._lookup();
    let L;
    if (level < 0) L = this.cfg.dc ||= this._build(0);
    else L = this.levels[level] ||= this._build(this.h / (1 << level));
    if (L === this.cur) return;
    this.X = L.X; this.cols = L.cols; this.P = L.P; this.Q = L.Q; this.CP = L.CP;
    this.wP.set(L.P); this.wQ.set(L.Q); this.wCP.set(L.CP); this.wqmax.set(L.qmax);
    // geq changed: re-derive each capacitor's history source from its state (v, i)
    const { capGeq, capIeq, capV, capI } = this;
    for (let i = 0; i < this.nc; i++) {
      const g = capGeq[i] = L.geq[i];
      capIeq[i] = level < 0 ? 0 : (this.caps[i].be ? g * capV[i] : g * capV[i] + capI[i]);
    }
    this.cur = L;
    this.level = level;
  }
  rebuild() { this.cache.clear(); this.dirty = true; this._use(this.level); }
  // find (or start) the matrix set for the current part values; keeps the 16 most recent
  _lookup() {
    let key = String(this.h);
    for (const x of this.vars) key += ',' + (x.g ?? x.c);
    let set = this.cache.get(key);
    if (set) this.cache.delete(key);
    else {
      set = { levels: [], dc: null };
      if (this.cache.size >= 16) this.cache.delete(this.cache.keys().next().value);
    }
    this.cache.set(key, set);
    this.levels = set.levels; this.cfg = set;
    this.dirty = false; this.cur = null;
  }

  _loadE() {
    const e = this.e, src = this.sources;
    for (let j = 0; j < this.ns; j++) e[j] = src[j].volts;
    K.load_caps(this.ctx, this.ns);
  }

  // Newton on the port vector (in the WebAssembly kernel). Returns true if converged.
  newton(maxIter = 60) { return K.newton(this.ctx, maxIter) === 1; }

  // value of row vector R (length k+m) against current excitation and port currents
  _read(R, o, e) {
    const k = this.k, m = this.m, fp = this.fp;
    let s = 0;
    for (let j = 0; j < k; j++) s += R[o + j] * e[j];
    o += k;
    for (let j = 0; j < m; j++) s += R[o + j] * fp[j];
    return s;
  }

  // capacitor states from the accepted solution; also keeps the excitation this solve used,
  // so v() after step() reads the same solution
  _acceptCaps(dc = false) { K.accept_caps(this.ctx, dc ? 1 : 0); }

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
    st.cv.set(this.capV); st.ci.set(this.capI);
    st.level = this.level;
  }
  _restore() {
    const st = this._st;
    this.vp.set(st.vp); this.eLast.set(st.e);
    this.capV.set(st.cv); this.capI.set(st.ci);
    const { capGeq, capIeq, capV, capI } = this;
    for (let i = 0; i < this.nc; i++) capIeq[i] = this.caps[i].be ? capGeq[i] * capV[i] : capGeq[i] * capV[i] + capI[i];
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
    if (this.dirty || this.buf !== buffer()) this._use(this.level);
    return this._read(this.X, idx * this.cols, this.eLast);
  }
}

export { VT };
