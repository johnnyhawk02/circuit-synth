// node test/render-wav.mjs  -> writes demo/patch.wav (16-bit mono)
import { writeFileSync } from 'node:fs';
import { renderPatch } from '../src/patch.js';
const sr = 48000;
const t0 = performance.now();
const { out, stats } = renderPatch({ sampleRate: sr, seconds: 3 });
const secs = (performance.now() - t0) / 1000;
const buf = Buffer.alloc(44 + out.length * 2);
buf.write('RIFF', 0); buf.writeUInt32LE(36 + out.length * 2, 4); buf.write('WAVE', 8);
buf.write('fmt ', 12); buf.writeUInt32LE(16, 16); buf.writeUInt16LE(1, 20); buf.writeUInt16LE(1, 22);
buf.writeUInt32LE(sr, 24); buf.writeUInt32LE(sr * 2, 28); buf.writeUInt16LE(2, 32); buf.writeUInt16LE(16, 34);
buf.write('data', 36); buf.writeUInt32LE(out.length * 2, 40);
let peak = 0;
out.forEach((s, i) => { peak = Math.max(peak, Math.abs(s)); buf.writeInt16LE(Math.max(-32767, Math.min(32767, Math.round(s * 32767))), 44 + i * 2); });
writeFileSync(new URL('../demo/patch.wav', import.meta.url), buf);
console.log(`rendered 3 s in ${secs.toFixed(2)} s (${(secs / 3 * 100).toFixed(0)}% of real time), peak ${peak.toFixed(3)}`);
