// node test/run-tests.mjs  — checks each module against hand-calculated values
import { mixer, lfo, vca } from '../src/modules.js';
import { vco, Filter904A, vca902, env911 } from '../src/moog.js';

const SR = 48000;
let fails = 0;
const check = (name, ok, detail) => {
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}  ${detail}`);
  if (!ok) fails++;
};

// --- mixer: out = -(in1 + in2 + in3)
{
  const m = mixer(SR);
  const o = m.process({ in1: 1, in2: 2, in3: -0.5 }).out;
  check('mixer sum', Math.abs(o + 2.5) < 1e-3, `out=${o.toFixed(4)} V (expect -2.5)`);
  const clip = m.process({ in1: 8, in2: 8, in3: 0 }).out;
  check('mixer rail clip', clip < -10 && clip > -10.6, `out=${clip.toFixed(3)} V (expect ~ -10.5)`);
}

// --- LFO: 5 Hz, triangle ~ ±5.25 V
{
  const l = lfo(SR);
  const secs = 3, n = SR * secs;
  let crossings = 0, firstEdge = -1, lastEdge = -1, prev = 0, tmax = -99, tmin = 99;
  const t0 = performance.now();
  for (let i = 0; i < n; i++) {
    const { tri, sq } = l.process();
    if (i > SR / 2) { tmax = Math.max(tmax, tri); tmin = Math.min(tmin, tri); }
    if (prev <= 0 && sq > 0) { if (i > SR / 2) { if (firstEdge < 0) firstEdge = i; lastEdge = i; crossings++; } }
    prev = sq;
  }
  const ms = performance.now() - t0;
  const f = (crossings - 1) / ((lastEdge - firstEdge) / SR);
  check('lfo frequency', Math.abs(f - 5) <= 0.25, `${f.toFixed(3)} Hz (expect 5)`);
  check('lfo triangle amplitude', Math.abs(tmax - 5.25) < 0.4 && Math.abs(tmin + 5.25) < 0.4,
    `+${tmax.toFixed(2)} / ${tmin.toFixed(2)} V (expect ±5.25)`);
  const s = l.circuit.stats;
  console.log(`      lfo: ${(s.iters / s.steps).toFixed(2)} Newton iters/step, ${s.failures} non-converged, ${s.retries || 0} sub-stepped, ` +
    `${(ms / 1000 / secs * 100).toFixed(1)}% of real time (4x oversampled)`);
}

// --- VCA: gain rises with CV, zero CV is (near) silent, big input soft-clips
{
  const v = vca(SR);
  const amp = (cv, inAmp) => {
    let peak = 0;
    for (let i = 0; i < SR / 10; i++) {
      const x = inAmp * Math.sin(2 * Math.PI * 440 * i / SR);
      const { out } = v.process({ in: x, cv });
      if (i > SR / 20) peak = Math.max(peak, Math.abs(out));
    }
    return peak;
  };
  const g0 = amp(0, 1), g5 = amp(5, 1), g10 = amp(10, 1);
  check('vca off at cv=0', g0 < 0.01, `peak=${g0.toFixed(4)} V`);
  check('vca gain grows with cv', g10 > g5 * 1.6 && g5 > 0.1, `cv5=${g5.toFixed(3)} V, cv10=${g10.toFixed(3)} V`);
  const small = amp(10, 0.5), big = amp(10, 5);
  check('vca soft clip (10x input < 10x output)', big < small * 10 * 0.85,
    `0.5V in -> ${small.toFixed(3)} V, 5V in -> ${big.toFixed(3)} V`);
  const s = v.circuit.stats;
  console.log(`      vca: ${(s.iters / s.steps).toFixed(2)} Newton iters/step, ${s.failures} non-converged`);
}

// --- VCO: 1 V/oct tracking (period from rising zero crossings of the saw about its mean)
{
  const o = vco(SR);
  const freq = cv => {
    const n = SR / 2, buf = new Float64Array(n);
    for (let i = 0; i < n; i++) buf[i] = o.process({ cv }).saw;
    let mean = 0; for (let i = n / 2; i < n; i++) mean += buf[i]; mean /= n / 2;
    const xs = [];
    for (let i = n / 2 + 1; i < n; i++) if (buf[i - 1] < mean && buf[i] >= mean) xs.push(i);
    return (xs.length - 1) * SR / (xs[xs.length - 1] - xs[0]);
  };
  const f0 = freq(0), f1 = freq(1), f2 = freq(2);
  const c1 = 1200 * Math.log2(f1 / f0 / 2), c2 = 1200 * Math.log2(f2 / f1 / 2);
  check('vco 1 V/oct', Math.abs(c1) < 25 && Math.abs(c2) < 25,
    `${f0.toFixed(1)} / ${f1.toFixed(1)} / ${f2.toFixed(1)} Hz, octave error ${c1.toFixed(1)}, ${c2.toFixed(1)} cents`);
}

// --- 904A: passes low frequencies, attenuates far above cutoff (24 dB/oct)
{
  const amp = f => {
    const fl = new Filter904A(SR);
    fl.setCutoff(0.5);
    const n = SR / 2, ys = new Float64Array(n);
    for (let i = 0; i < n; i++) ys[i] = fl.process(0.1 * Math.sin(2 * Math.PI * f * i / SR));
    let mean = 0, peak = 0;
    for (let i = n / 2; i < n; i++) mean += ys[i]; mean /= n / 2;
    for (let i = n / 2; i < n; i++) peak = Math.max(peak, Math.abs(ys[i] - mean));
    return peak;
  };
  const lo = amp(50), hi = amp(1600);
  const db = 20 * Math.log10(hi / lo);
  check('904A low-pass rolloff', db < -40, `50 Hz ${lo.toFixed(4)} V, 1.6 kHz ${hi.toFixed(5)} V (${db.toFixed(1)} dB)`);
}

// --- 902: gain grows exponentially with control voltage
{
  const m = vca902(SR);
  const amp = cv => {
    let peak = 0, mean = 0;
    const n = SR / 5, ys = new Float64Array(n);
    for (let i = 0; i < n; i++) ys[i] = m.process({ in: 0.5 * Math.sin(2 * Math.PI * 440 * i / SR), cv1: cv }).o12;
    for (let i = n / 2; i < n; i++) mean += ys[i]; mean /= n / 2;
    for (let i = n / 2; i < n; i++) peak = Math.max(peak, Math.abs(ys[i] - mean));
    return peak;
  };
  const a1 = amp(1), a2 = amp(2), a3 = amp(3);
  check('902 exponential gain', a2 / a1 > 2 && a3 / a2 > 2, `cv 1/2/3 V -> ${a1.toFixed(4)} / ${a2.toFixed(4)} / ${a3.toFixed(4)} V`);
}

// --- 911: attack rises, settles to sustain, releases to ~0
{
  const e = env911(SR);
  e.setAttack(0.1); e.setDecay(0.3); e.setSustain(0.5); e.setRelease(0.1);
  const run = secs => { let y = 0, pk = 0; for (let i = 0; i < secs * SR / e.decimate; i++) { y = e.process({}).out; pk = Math.max(pk, y); } return [y, pk]; };
  const [idle] = run(0.1);
  e.gate(true); const [held, peak] = run(1.5);
  e.gate(false); const [rel] = run(1.5);
  check('911 envelope ADSR shape', peak > 4 && held < peak * 0.8 && held > 1 && Math.abs(rel) < 0.2 && Math.abs(idle) < 0.2,
    `idle ${idle.toFixed(3)}, peak ${peak.toFixed(2)}, sustain ${held.toFixed(2)}, released ${rel.toFixed(3)} V`);
}

console.log(fails ? `\n${fails} test(s) failed` : '\nall tests passed');
process.exit(fails ? 1 : 0);
