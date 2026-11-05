// A tiny patch: LFO triangle -> VCA CV, test sine -> VCA input -> mixer.
// Everything below the oscillator is solved at component level.
import { lfo, vca, mixer } from './modules.js';

export function renderPatch({ sampleRate = 48000, seconds = 3, toneHz = 220, lfoRateOhms = 100e3,
                              inputVolts = 5, onProgress } = {}) {
  const l = lfo(sampleRate, { rate: lfoRateOhms });
  const v = vca(sampleRate);
  const m = mixer(sampleRate);
  const n = Math.floor(sampleRate * seconds);
  const out = new Float32Array(n), cvTrace = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const tri = l.process().tri;                // ±5.25 V
    const cv = tri + 5.25;                      // shift to 0..10.5 V
    const x = inputVolts * Math.sin(2 * Math.PI * toneHz * i / sampleRate);
    const y = v.process({ in: x, cv }).out;
    out[i] = m.process({ in1: y }).out / 5;     // ±5 V synth level -> ±1 audio
    cvTrace[i] = cv;
    if (onProgress && i % 4800 === 0) onProgress(i / n);
  }
  return { out, cvTrace, stats: { lfo: l.circuit.stats, vca: v.circuit.stats } };
}
