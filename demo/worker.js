import { renderPatch } from '../src/patch.js';
self.onmessage = e => {
  const t0 = performance.now();
  const res = renderPatch({ ...e.data, onProgress: p => self.postMessage({ progress: p }) });
  self.postMessage({ done: true, out: res.out, cv: res.cvTrace, ms: performance.now() - t0, stats: res.stats });
};
