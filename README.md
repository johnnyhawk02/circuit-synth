# Component-level synth prototype (JS)

Every module is a real schematic written as a netlist and solved sample by sample in real time.
No DSP shortcuts inside the modules: resistors, capacitors, transistors, JFETs, diodes and
op-amp macro-models, nothing else.

## Live demo

GitHub Pages serves the browser demo from the repo root (`index.html` redirects to `demo/`).

- `demo/index.html` – **live Moog patch**: 8-step sequencer → VCO → 904A ladder filter →
  902 VCA, with a 911 envelope on the VCA and (through an attenuator) on the filter.
  Runs in an AudioWorklet; every knob changes a real component value (pot position, cap value).
- `demo/lfo-vca.html` – the earlier LFO → diff-pair VCA demo.

Quality: **Eco** (VCO 2x oversampled, filter 1x) runs at roughly 50% of one core on a
2.1 GHz server Xeon (headless Chrome), so a recent laptop has plenty of headroom. **HQ** (4x / 2x)
runs at about 85%.

## Run

```
node test/run-tests.mjs        # checks each module against hand-calculated values
node test/render-moog.mjs 4    # writes demo/moog.wav (the live patch, 4 s, offline)
node test/render-wav.mjs       # writes demo/patch.wav (old LFO -> VCA patch)
python3 -m http.server         # then open http://localhost:8000/demo/
```
The browser demo needs a local server (ES modules and worklets don't load from file://).

## Files

- `src/solver.js` – the circuit engine, DK-method style: the linear network is folded into
  precomputed matrices per step size, and Newton runs only on the nonlinear port voltages.
  Parts: R, POT, C (trapezoidal or backward Euler per cap), V (settable inputs), diode, NPN/PNP
  (Ebers–Moll, or a cheaper forward-only model; optional junction capacitance), JFET (saturation
  and a triode switch), op-amp (tanh rail clip, or a gain-bandwidth pole + slew limit for
  comparators). Sub-stepping (h/2…h/256, cached matrices) on failed steps and while a watched
  node is mid-transition; linear predictor for the port voltages.
- `kernel/kernel.c` – the per-sample hot path (Newton on the port voltages, port models,
  capacitor update) in plain C, compiled to WebAssembly with SIMD; about 1.7x faster than the
  same code in JavaScript, and reusable as-is in a C++ plugin. `sh kernel/build.sh` rebuilds
  it (needs clang with the wasm32 target) and embeds it in `src/kernel-wasm.js`, which
  `src/kernel.js` instantiates synchronously, so nothing else changes for Node or the worklet.
  Ports whose slope is negligible (a cut-off transistor, an idle clamp) are left out of the
  factorisation, which makes the 911's 25-port flip-flop circuit about 35% cheaper.
- `src/modules.js` – mixer, LFO, diff-pair VCA (first prototypes).
- `src/moog.js` – the Moog modules, from the original service drawings:
  - **vco**: exponential converter (matched pair + reference op-amp), integrator, JFET reset,
    Schmitt comparator. Generic saw core (the 921 uses ICs we don't model yet). 1 V/oct within
    ~10 cents over three octaves, falls flat at the top (−97 cents at C7).
  - **Filter904A**: input buffer Q22/Q23, adder + Q6 exponential converter, the 4-pole ladder
    Q7–Q16 (bias chain 680/220/220/220/68 Ω, range caps 1.2 / 0.3 / 0.075 µF), output buffer
    Q17–Q21. DC points match the drawing's annotations within ~0.1 V.
  - **vca902**: adder, Q10 exponential converter, Q6–Q9 gain cell, Q11/Q12 PNP output pair,
    Q13 current source. `trimBias()` sets OUTPUT BIAS like the service procedure.
  - **env911**: 10 µF timing cap, 1 MΩ attack/decay/release pots, JFET buffer, discrete
    transistor flip-flop (Q9–Q14). Zero trim set automatically.
- `src/voice.js` – the patch and step sequencer, shared by the worklet and the offline renderer.

## Assumptions (things the drawings don't settle)

- 904A regeneration is fed back from the output pair (Q20 collector) through C8 into R22.
- 904A modules are solved as three circuits joined at transistor bases; the ladder sees the
  adder's output as an ideal source (the loading there is one base current).
- 911 flip-flop transistors get 1 nF junction capacitance so the latch flips in a few sub-steps.
- A 4.7 µF patch-cable coupling cap in front of the 904A (the VCO saw sits on a DC offset).
- The VCO is a generic expo/integrator/Schmitt core, not a 921 copy.

## Known limits / next steps

- Each module is solved separately and joined by ideal voltage ports (no loading between modules).
- Speed: the patch takes about half a core. The kernel is already C; a native plugin would
  also move the matrix building and the module glue out of JavaScript.
- VCO high-frequency tracking needs a bulk-resistance compensation like the real 921.
