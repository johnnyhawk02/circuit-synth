// Loads the WebAssembly solver kernel (kernel/kernel.c) synchronously and hands out memory
// for each circuit's arrays. Works in Node, the browser and an AudioWorklet.
import { KERNEL_WASM } from './kernel-wasm.js';

function decode(b64) {
  if (typeof atob === 'function') {
    const s = atob(b64), out = new Uint8Array(s.length);
    for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
    return out;
  }
  // AudioWorkletGlobalScope has no atob: decode by hand
  const T = new Uint8Array(128);
  const A = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
  for (let i = 0; i < 64; i++) T[A.charCodeAt(i)] = i;
  const n = b64.length, pad = b64.endsWith('==') ? 2 : b64.endsWith('=') ? 1 : 0;
  const out = new Uint8Array(n / 4 * 3 - pad);
  for (let i = 0, o = 0; i < n; i += 4) {
    const v = T[b64.charCodeAt(i)] << 18 | T[b64.charCodeAt(i + 1)] << 12 | T[b64.charCodeAt(i + 2)] << 6 | T[b64.charCodeAt(i + 3)];
    if (o < out.length) out[o++] = v >> 16 & 255;
    if (o < out.length) out[o++] = v >> 8 & 255;
    if (o < out.length) out[o++] = v & 255;
  }
  return out;
}

const instance = new WebAssembly.Instance(new WebAssembly.Module(decode(KERNEL_WASM)), {});
export const K = instance.exports;
const memory = K.memory;
let top = (K.heap_base() + 15) & ~15;

// Bump allocator (circuits live as long as the page). Growing the memory replaces its
// ArrayBuffer, so circuits compare their cached buffer with buffer() and re-make views.
export function alloc(bytes) {
  const at = top;
  top = (top + bytes + 15) & ~15;
  const need = top - memory.buffer.byteLength;
  if (need > 0) memory.grow(Math.ceil(need / 65536) + 16);
  return at;
}
export const buffer = () => memory.buffer;

// Ctx field order, matching the struct in kernel.c (every field is 4 bytes)
export const CTX = [
  'm', 'k', 'nc', 'cols', 'ncross', 'luOK', 'na', 'iters', 'factors',
  'P', 'Q', 'CP', 'qmax',
  'e', 'eLast', 'vp', 'fp', 'pl', 'dd', 'dv', 'r', 'Ja', 'ra', 'ddF', 'tolv', 'xdF',
  'act', 'pos', 'piv', 'keep', 'isAct',
  'kind', 'partner',
  'is', 'nvt', 'vcrit', 'vsat', 'kk', 'vos', 'tanhStep', 'beta', 'vto', 'soft', 'xd',
  'crossI', 'crossX',
  'capGeq', 'capIeq', 'capV', 'capI', 'capBe',
];
