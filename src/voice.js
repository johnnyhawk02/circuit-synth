// The patch: step sequencer -> VCO -> 904A -> 902, with a 911 envelope on the 902 and
// (through an attenuator) on the 904A. Runs sample by sample; used by the AudioWorklet and
// by the offline renderer.
import { vco, Filter904A, vca902, env911 } from './moog.js';

export const QUALITY = {
  eco: { vcoOS: 2, filterOS: 1 },
  hq: { vcoOS: 4, filterOS: 2 },
};

export class MoogPatch {
  constructor(sampleRate, { quality = 'eco' } = {}) {
    const q = QUALITY[quality] || QUALITY.eco;
    this.sr = sampleRate;
    this.vco = vco(sampleRate, { oversample: q.vcoOS });
    this.vcf = new Filter904A(sampleRate, { oversample: q.filterOS });
    this.vca = vca902(sampleRate);
    this.vca.trimBias();
    this.eg = env911(sampleRate);
    this.p = {
      tempo: 120, gateLen: 0.5, running: true,
      steps: [0, 12, 0, 10, 0, 7, 3, 12].map(n => ({ note: n, on: true })),
      envAmt: 0.5, kbd: 0.5, transpose: 0, level: 0.2,
    };
    this.set({ cutoff: 0.45, resonance: 0.4, attack: 0.05, decay: 0.35, sustain: 0.3, release: 0.25, range: 2 });
    this.t = 0; this.step = -1; this.gate = false; this.cv = 0;
    this.envPrev = 0; this.envNext = 0; this.envPhase = 0;
    this.hp = 0; this.hpPrev = 0;
  }
  set(o) {
    for (const [k, v] of Object.entries(o)) {
      switch (k) {
        case 'cutoff': this.vcf.setCutoff(v); break;
        case 'resonance': this.vcf.setResonance(v); break;
        case 'range': this.vcf.setRange(v); break;
        case 'attack': this.eg.setAttack(v); break;
        case 'decay': this.eg.setDecay(v); break;
        case 'sustain': this.eg.setSustain(v); break;
        case 'release': this.eg.setRelease(v); break;
        case 'steps': this.p.steps = v; break;
        case 'tune': this.vco.setTune(v); break;
        default: this.p[k] = v;
      }
    }
  }
  // fill `out` (Float32Array); optional `scope` gets {audio, env} per sample
  render(out, scope) {
    const p = this.p, sr = this.sr, stepLen = sr * 60 / p.tempo / 4; // 16th notes
    const dec = this.eg.decimate;
    for (let i = 0; i < out.length; i++) {
      // --- sequencer (digital, like a gate/CV source patched in) ---
      if (p.running) {
        const pos = this.t / stepLen, s = Math.floor(pos) % p.steps.length;
        if (s !== this.step) {
          this.step = s;
          const st = p.steps[s];
          if (st.on) { this.cv = (st.note + p.transpose) / 12; this.gate = true; this.eg.gate(true); }
        }
        if (this.gate && pos - Math.floor(pos) >= p.gateLen) { this.gate = false; this.eg.gate(false); }
        this.t++;
      } else if (this.gate) { this.gate = false; this.eg.gate(false); }
      // --- 911 at sr/decimate, linearly interpolated ---
      if (this.envPhase === 0) { this.envPrev = this.envNext; this.envNext = this.eg.process().out; }
      const env = this.envPrev + (this.envNext - this.envPrev) * (this.envPhase / dec);
      this.envPhase = (this.envPhase + 1) % dec;
      // --- audio path ---
      const saw = this.vco.process({ cv: this.cv }).saw;
      const f = this.vcf.process(saw, this.cv * p.kbd, env * p.envAmt, 0);
      const y = this.vca.process({ in: f, cv1: env }).o12;
      // output DC blocker (the audio interface's coupling cap), 5 Hz
      this.hp = 0.99935 * (this.hp + y - this.hpPrev); this.hpPrev = y;
      out[i] = this.hp * p.level;
      if (scope) scope(out[i], env, saw);
    }
    return out;
  }
}
