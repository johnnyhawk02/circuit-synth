# Component-level synth prototype (JS)

Every module is a real schematic written as a netlist and solved sample by sample:
Modified Nodal Analysis + Newton-Raphson, trapezoidal capacitors, 4x oversampling.
No DSP shortcuts inside the modules.

## Run

```
node test/run-tests.mjs      # checks each module against hand-calculated values
node test/render-wav.mjs     # writes demo/patch.wav (LFO -> VCA patch, 3 s)
python3 -m http.server       # then open http://localhost:8000/demo/
```
The browser demo needs a local server (ES modules don't load from file://).

## Files

- `src/solver.js` – the circuit engine. Parts: R, C, V (settable inputs), diode (Shockley),
  NPN/PNP (Ebers-Moll), op-amp (tanh rail-clipping macro-model with input offset).
  Junction voltage limiting, Newton damping, automatic half-step retry when a step fails.
- `src/modules.js` – the schematics:
  - **mixer**: 3-input inverting op-amp summer, 100k resistors.
  - **lfo**: op-amp integrator + Schmitt trigger. f = R3 / (4·Rrate·C·R2) = 5 Hz at 100k / 1 µF.
  - **vca**: discrete BJT differential pair, op-amp CV-to-tail-current, op-amp difference amp.
- `src/patch.js` – patches them: LFO triangle -> VCA CV, sine -> VCA -> mixer.

## Adding a module

```js
const c = new Circuit();
c.V('vcc', '0', 12); c.V('vee', '0', -12);
c.V('in_src', '0', 0, 'in');          // named input you set each sample
c.R('in_src', 'n1', 10e3);
c.NPN('col', 'base', 'emit');        // c, b, e
c.OPAMP('plus', 'minus', 'out');
new Module(c, ['in'], ['out'], 48000, 4);
```

## Known limits / next steps

- Speed: ~1.1x real time in Node for LFO + VCA + mixer, dominated by the dense LU solve.
  Next: DK method (precompute the linear part, iterate only on nonlinear ports) or sparse LU.
- The Schmitt trigger's regenerative flip has no continuous solution with an ideal-bandwidth
  op-amp; ~5 steps per second don't fully converge but the waveform is correct (5.004 Hz measured).
  Next: add gain-bandwidth to the op-amp model.
- Modules are solved separately and joined by ideal voltage ports (no loading between them).
- Next modules from the schematic guide: Thomas Henry VCO-1 core, Moog transistor ladder, ADSR.
