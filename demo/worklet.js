// AudioWorklet: runs the whole component-level patch in the audio thread.
import { MoogPatch } from '../src/voice.js';

class MoogProcessor extends AudioWorkletProcessor {
  constructor(options) {
    super();
    const t0 = Date.now();
    this.patch = new MoogPatch(sampleRate, { quality: options.processorOptions?.quality || 'eco' });
    if (options.processorOptions?.params) this.patch.set(options.processorOptions.params);
    this.scopeA = new Float32Array(2048); this.scopeE = new Float32Array(2048); this.si = 0;
    this.busyMs = 0; this.blocks = 0; this.lastStep = -1;
    this.alive = true;
    this.port.onmessage = e => { if (e.data.set) this.patch.set(e.data.set); if (e.data.stop) this.alive = false; };
    this.port.postMessage({ ready: true, initMs: Date.now() - t0 });
  }
  process(_inputs, outputs) {
    if (!this.alive) return false;
    const out = outputs[0][0];
    const t0 = Date.now(); // ms resolution, but unbiased when averaged over many blocks
    this.patch.render(out, (a, env) => {
      this.scopeA[this.si] = a; this.scopeE[this.si] = env; this.si = (this.si + 1) & 2047;
    });
    for (let c = 1; c < outputs[0].length; c++) outputs[0][c].set(out);
    this.busyMs += Date.now() - t0;
    this.blocks++;
    if (this.patch.step !== this.lastStep) { this.lastStep = this.patch.step; this.port.postMessage({ step: this.lastStep }); }
    if (this.blocks % 32 === 0) {
      const blockMs = 32 * out.length / sampleRate * 1000;
      const a = new Float32Array(2048), e = new Float32Array(2048);
      for (let i = 0; i < 2048; i++) { a[i] = this.scopeA[(this.si + i) & 2047]; e[i] = this.scopeE[(this.si + i) & 2047]; }
      this.port.postMessage({ load: this.busyMs / blockMs, scope: a, env: e }, [a.buffer, e.buffer]);
      this.busyMs = 0;
    }
    return true;
  }
}
registerProcessor('moog-patch', MoogProcessor);
