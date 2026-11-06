// node test/render-moog.mjs [seconds] [eco|hq]  -> writes demo/moog.wav (16-bit mono)
// Renders the live demo's patch (sequencer -> VCO -> 904A -> 902, 911 envelope) offline.
import { writeFileSync } from 'node:fs';
import { MoogPatch } from '../src/voice.js';
const sr = 48000, seconds = +(process.argv[2] || 4), quality = process.argv[3] || 'eco';
const t0 = performance.now();
const patch = new MoogPatch(sr, { quality });
const init = performance.now() - t0;
const out = new Float32Array(sr * seconds), block = new Float32Array(128);
const t1 = performance.now();
for (let i = 0; i < out.length; i += 128) { patch.render(block); out.set(block.subarray(0, Math.min(128, out.length - i)), i); }
const secs = (performance.now() - t1) / 1000;
const buf = Buffer.alloc(44 + out.length * 2);
buf.write('RIFF', 0); buf.writeUInt32LE(36 + out.length * 2, 4); buf.write('WAVE', 8);
buf.write('fmt ', 12); buf.writeUInt32LE(16, 16); buf.writeUInt16LE(1, 20); buf.writeUInt16LE(1, 22);
buf.writeUInt32LE(sr, 24); buf.writeUInt32LE(sr * 2, 28); buf.writeUInt16LE(2, 32); buf.writeUInt16LE(16, 34);
buf.write('data', 36); buf.writeUInt32LE(out.length * 2, 40);
let peak = 0;
out.forEach((s, i) => { peak = Math.max(peak, Math.abs(s)); buf.writeInt16LE(Math.max(-32767, Math.min(32767, Math.round(s * 32767))), 44 + i * 2); });
writeFileSync(new URL('../demo/moog.wav', import.meta.url), buf);
console.log(`${quality}: start-up ${init.toFixed(0)} ms, rendered ${seconds} s in ${secs.toFixed(2)} s (${(secs / seconds * 100).toFixed(0)}% of real time), peak ${peak.toFixed(3)}`);
