// node test/run-tests.mjs  — checks each module against hand-calculated values
import { mixer, lfo, vca } from '../src/modules.js';

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

console.log(fails ? `\n${fails} test(s) failed` : '\nall tests passed');
process.exit(fails ? 1 : 0);
