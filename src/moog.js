// Moog modular modules at component level, from the original schematics:
//   904A Voltage Controlled Low Pass Filter  (drawing 1149, rev E)
//   902  Voltage Controlled Amplifier         (service manual fig. 9, 993-041813)
//   911  Envelope Generator                   (service manual fig. 17, 993-041815)
// plus a discrete saw-core VCO (not a Moog design, see vco() below).
//
// Part numbers in comments (R12, Q7 ...) are the schematic designators. NPNs are 2N3392,
// PNPs 2N4058 in all three Moog modules. Where the scan left something open, the choice made
// here is written next to it as ASSUMPTION.
import { Circuit } from './solver.js';
import { Module } from './modules.js';

const NPN_3392 = { bf: 250 };          // 2N3392: hFE 150-300
const SIMPLE = { ...NPN_3392, simple: true };

// ----------------------------------------------------------------------------------------------
// VCO: classic exponential saw core. Not a published schematic: an expo pair with an op-amp
// reference, an op-amp integrator, a Schmitt-trigger comparator and an NPN reset transistor.
//   f = Iref * 2^CV / (C * (Vupper - Vlower)),  C = 10 nF, thresholds 5.0 V / 0.5 V
//   Iref = 12 V / Rref; Rref = 4.08 M puts CV = 0 at 65.4 Hz (C2).
// ----------------------------------------------------------------------------------------------
export function vco(sampleRate, { oversample = 2 } = {}) {
  const c = new Circuit();
  c.V('vcc', '0', 12); c.V('vee', '0', -12);
  c.V('cv_src', '0', 0, 'cv');
  // 1 V/oct scaling: 100k / 1.823k divider gives VT*ln2 = 17.92 mV per volt at the base
  c.R('cv_src', 'b2', 100e3);
  c.R('b2', '0', 1.823e3);
  // expo pair: reference transistor Q1 held at Iref by op-amp A1
  const rref = c.R('vcc', 'ref_n', 4.08e6);
  c.NPN('ref_n', '0', 'ee', SIMPLE);          // Q1
  c.OPAMP('0', 'ref_n', 'a1o');
  c.R('a1o', 'ee', 4.7e3);
  c.NPN('sum', 'b2', 'ee', SIMPLE);            // Q2, sinks the core current from the integrator
  // integrator A2: saw ramps up at I / C
  c.C('sum', 'saw', 10e-9);
  c.OPAMP('0', 'sum', 'saw');
  // reset JFET (with 220Ω in series) across the integrating capacitor; no gate current reaches
  // the summing node. The gate sits at (cmp - 12 V)/2: off at -11 V, on (-0.75 V) once the
  // comparator is past +6 V, so the Schmitt has latched before the discharge starts.
  c.NJFS('sawr', 'jg', 'sum');
  c.R('saw', 'sawr', 220);
  c.R('cmp', 'jg', 47e3); c.R('jg', 'vee', 47e3);
  // Schmitt trigger A3: upper threshold 5.0 V, lower 0.5 V
  c.R('vcc', 'vref', 43e3); c.R('vref', '0', 10e3);
  c.R('saw', 'sp', 21.4e3); c.R('cmp', 'sp', 100e3);
  c.OPAMP('sp', 'vref', 'cmp', { vos: 1e-3, gbw: 1e6, slew: 13e6 });
  // The saw is taken straight from the integrator op-amp output (0.5 V .. 5 V ramp).
  // sub-step while the comparator is between its rails (h/16 at 4x: 0.33 µs)
  c.watch('cmp', -10, 10, oversample >= 4 ? 4 : 5);  // rarely fires: failed steps sub-step too
  const m = new Module(c, ['cv'], ['saw'], sampleRate, oversample, { init: 'zero' });
  m.setTune = semis => rref.set(4.08e6 / Math.pow(2, semis / 12));
  return m;
}

// ----------------------------------------------------------------------------------------------
// 904A low pass filter. Three circuits stepped in lockstep, joined at transistor bases
// (where the loading is a base current only):
//   input buffer  Q22/Q23 + regeneration pot R22
//   ladder        adder (op-amp macro for the discrete Q1-Q5 adder) + Q6 expo + Q7..Q16 ladder
//   output buffer Q17/Q21 followers, Q19/Q20 PNP pair, Q18 current source
// ----------------------------------------------------------------------------------------------
export class Filter904A {
  constructor(sampleRate, { oversample = 2 } = {}) {
    this.os = oversample;
    const fs = sampleRate * oversample;

    // --- input buffer ---
    const ib = new Circuit();
    ib.V('p12', '0', 12); ib.V('m6', '0', -6);
    ib.V('sig', '0', 0, 'in');
    ib.V('fb', '0', 0, 'fb');                    // from the output buffer (Q20 collector)
    // patch-cable coupling capacitor (not on the 904A drawing: the VCO saw has a DC offset)
    ib.C('sig', 'sigc', 4.7e-6);
    ib.R('sigc', 'b23', 8.2e3); ib.R('b23', '0', 1e3);         // R20, R21
    // ASSUMPTION: R22's far end is ground; its other end gets the output through C8
    ib.C('fb', 'c8n', 2.5e-6);                                // C8 2.5 µF
    this.regen = ib.POT('0', 'c8n', 'rw', 50e3, 0);          // R22 50k REGENERATION
    ib.R('rw', 'b22', 1.8e3); ib.R('b22', '0', 1e3);          // R11 1.8k select, R10 1k
    ib.NPN('c22', 'b22', 'e22', SIMPLE);                      // Q22
    ib.NPN('c23', 'b23', 'e23', SIMPLE);                      // Q23
    ib.R('e22', 'e23', 1e3);                                  // R16
    ib.R('e22', 'm6', 10e3); ib.R('e23', 'm6', 10e3);         // R15, R18
    ib.R('c22', 'r17a', 100); ib.R('c23', 'r17b', 100);       // R14, R19
    ib.POT('r17a', 'r17b', 's', 10, 0.5);                     // R17 10Ω balance
    ib.R('s', '0', 330); ib.C('s', '0', 10e-6);              // R12₂ 330, C3 10 µF
    ib.R('s', 'w', 330); ib.R('w', 'p12', 330);              // R12, R42
    this.ib = ib.probe('c22', 'c23', 'w');
    ib.finalize(fs); ib.dcop();

    // --- adder + expo bias (control path, its own small circuit; joins the ladder at Q6's
    // emitter, which the adder output drives at low impedance, and Q6's diode-biased base) ---
    const ct = new Circuit();
    ct.V('p12', '0', 12); ct.V('m6', '0', -6);
    for (const k of ['cv1', 'cv2', 'cv3']) { ct.V(k + '_s', '0', 0, k); ct.R(k + '_s', 'sum', 100e3); } // R1-R3
    this.freq = ct.POT('m6', '0', 'fcw', 5e3, 0.5);         // R11₂ 5k FIXED CONTROL VOLTAGE
    ct.R('fcw', 'sum', 51e3);                                 // R12₁ 51k
    // feedback R8 (1k trim) + R9 1.3k; trimmed to 0.49k so the loop gives VT*ln2 per volt
    ct.R('sum', 'aout', 1.3e3 + 0.49e3);
    ct.OPAMP('0', 'sum', 'aout');
    // Q6's base: one diode drop (CL-1 fed from +12 via R30 750, R24 2.5k trim)
    ct.R('p12', 'r30', 750); ct.R('r30', 'b6', 1.25e3); ct.D('b6', '0');
    this.ct = ct.probe('aout', 'b6');
    ct.finalize(fs); ct.dcop();

    // --- ladder ---
    const ld = new Circuit();
    ld.V('p12', '0', 12);
    ld.V('b8_src', '0', ib.v('c22'), 'b8');
    ld.V('b7_src', '0', ib.v('c23'), 'b7');
    ld.V('aout', '0', ct.v('aout'), 'aout');
    ld.V('b6', '0', ct.v('b6'), 'b6');
    ld.NPN('c6', 'b6', 'aout', SIMPLE);                       // Q6
    ld.R('c6', 'w25', 470);                                   // R23
    // R25 balance trim (value not legible: ASSUMPTION 100Ω, centred) with R45 47Ω across
    ld.POT('e8', 'e7', 'w25', 100, 0.5); ld.R('e8', 'e7', 47);
    ld.NPN('u1', 'b8_src', 'e8', SIMPLE);                     // Q8
    ld.NPN('l1', 'b7_src', 'e7', SIMPLE);                     // Q7
    // ladder rungs: bases on the bias chain R25₂ 680 / R26-28 220 / R29 68 from +12
    ld.R('0', 'bias1', 680); ld.R('bias1', 'bias2', 220); ld.R('bias2', 'bias3', 220);
    ld.R('bias3', 'bias4', 220); ld.R('bias4', 'p12', 68);
    const rung = [['u1', 'l1', 'u2', 'l2', 'bias1'], ['u2', 'l2', 'u3', 'l3', 'bias2'],
                  ['u3', 'l3', 'u4', 'l4', 'bias3'], ['u4', 'l4', 'p12', 'p12', 'bias4']];
    for (const [ue, le, uc, lc, b] of rung) {               // Q10/Q9, Q12/Q11, Q14/Q13, Q16/Q15
      ld.NPN(uc, b, ue, SIMPLE);
      ld.NPN(lc, b, le, SIMPLE);
    }
    this.caps = ['1', '2', '3', '4'].map(k => ld.C('u' + k, 'l' + k, 0.075e-6)); // S1 range 3
    this.ld = ld.probe('u4', 'l4', 'aout');
    ld.finalize(fs); ld.dcop();

    // --- output buffer ---
    const ob = new Circuit();
    ob.V('p12', '0', 12); ob.V('m6', '0', -6);
    ob.V('u4s', '0', ld.v('u4'), 'u4'); ob.V('l4s', '0', ld.v('l4'), 'l4');
    // ASSUMPTION: the R41/R40/R34 bias string is modelled as 33k from the R40 wiper (the
    // ~8.1 V node W shared with the input buffer) to each PNP base; 330k to ground each.
    ob.V('w', '0', ib.v('w'));
    ob.NPN('p12', 'u4s', 'e21', SIMPLE); ob.R('e21', 'm6', 47e3); // Q21, R43
    ob.NPN('p12', 'l4s', 'e17', SIMPLE); ob.R('e17', 'm6', 47e3); // Q17, R31
    ob.C('e21', 'b20', 10e-6); ob.C('e17', 'b19', 10e-6);          // C4, C9
    ob.R('w', 'b20', 33e3); ob.R('b20', '0', 330e3);
    ob.R('w', 'b19', 33e3); ob.R('b19', '0', 330e3);
    const PNP = { bf: 200, simple: true };
    ob.PNP('c20', 'b20', 'et', PNP);                         // Q20
    ob.PNP('c19', 'b19', 'et', PNP);                         // Q19
    ob.PNP('et', 'b18', 'e18', PNP);                         // Q18 tail current source
    ob.R('e18', 'p12', 68); ob.R('b18', 'p12', 470); ob.R('b18', '0', 2.2e3); // R38, R37, R36
    ob.R('c20', 'bn', 680); ob.R('c19', 'bn', 680);          // R39, R35
    ob.R('bn', 'm6', 33); ob.C('bn', '0', 400e-6);           // R46, C6
    ob.C('c19', 'out', 2.5e-6); ob.R('out', '0', 33e3);      // C7, R32 -> SIGNAL OUTPUT
    this.ob = ob.probe('out', 'c20');
    ob.finalize(fs); ob.dcop();
    ob.predict = false; // its current-fed tail node makes extrapolated guesses fail to converge

    // regeneration: the scan's C8 connection doesn't close a resonant (negative at fc) loop
    // with Q19's collector, so it is taken from the opposite phase, Q20's collector.
    this.fbDC = ob.v('c20');
    ib.set('fb', this.fbDC);
    ib.dcop();
    this.out = 0;
  }
  setCutoff(pos) { this.freq.set(pos); }           // fixed control pot, 0..1
  setResonance(pos) { this.regen.set(pos); }       // R22, 0..1
  setRange(i) { const f = [1.2e-6, 0.3e-6, 0.075e-6][i]; for (const c of this.caps) c.set(f); }
  // in: audio volts; cv1..cv3: control inputs (1 V/oct)
  process(input, cv1 = 0, cv2 = 0, cv3 = 0) {
    const { ib, ld, ob } = this;
    ib.set('in', input);
    const { ct } = this;
    ct.set('cv1', cv1); ct.set('cv2', cv2); ct.set('cv3', cv3);
    ct.step();
    ld.set('aout', ct.v('aout')); ld.set('b6', ct.v('b6'));
    let acc = 0;
    for (let s = 0; s < this.os; s++) {
      ib.step();
      ld.set('b8', ib.v('c22')); ld.set('b7', ib.v('c23'));
      ld.step();
      ob.set('u4', ld.v('u4')); ob.set('l4', ld.v('l4'));
      ob.step();
      ib.set('fb', ob.v('c20'));
      acc += ob.v('out');
    }
    return (this.out = acc / this.os);
  }
}

// ----------------------------------------------------------------------------------------------
// 902 VCA (EXP mode). Adder -> Q10 expo -> Q8/Q9 gain pair fed by the Q6/Q7 input pair,
// Q11/Q12 PNP output pair with Q13 current source.
// ----------------------------------------------------------------------------------------------
export function vca902(sampleRate, { oversample = 1 } = {}) {
  const c = new Circuit();
  c.V('p12', '0', 12); c.V('m6', '0', -6);
  c.V('sig', '0', 0, 'in');
  for (const k of ['cv1', 'cv2', 'cv3']) { c.V(k + '_s', '0', 0, k); c.R(k + '_s', 'sum', 100e3); } // R13-15
  const fixed = c.POT('m6', '0', 'fcw', 5e3, 0);   // FIXED CONTROL VOLTAGE (0 = fully off)
  c.R('fcw', 'sum', 100e3);                        // R39
  c.R('sum', 'aout', 3.9e3);                       // R17 (EXP position of SW1)
  c.OPAMP('0', 'sum', 'aout');
  c.R('p12', 'r36', 750); c.R('r36', 'b10', 1.25e3); c.D('b10', '0'); // R36, P4 INPUT BIAS, CL-1
  c.NPN('c10', 'b10', 'aout', SIMPLE);             // Q10 expo converter
  c.R('c10', 'w35', 470);                          // R35
  c.POT('e8', 'e9', 'w35', 100, 0.5); c.R('e8', 'e9', 100); // Q8-Q9 balance pot, 100Ω
  // input pair
  c.R('sig', 'b6', 8.2e3); c.R('b6', '0', 1e3);    // R28, R29 (lower input)
  c.R('b7', '0', 1e3);                              // R34 (upper input R33 unused)
  c.NPN('c6', 'b6', 'e6', SIMPLE); c.NPN('c7', 'b7', 'e7', SIMPLE); // Q6, Q7
  c.R('e6', 'e7', 1e3); c.R('e6', 'm6', 10e3); c.R('e7', 'm6', 10e3); // R30, R31, R32
  c.R('c6', 'p1a', 100); c.R('c7', 'p1b', 100);    // R18, R19
  c.POT('p1a', 'p1b', 'y', 25, 0.5);               // P1 INPUT BALANCE
  c.R('y', '0', 330); c.C('y', '0', 10e-6);        // R7, C2
  c.R('y', 'z', 330); c.R('z', 'p12', 470);        // R8, R10
  // gain pair
  c.NPN('c8', 'c6', 'e8', SIMPLE); c.NPN('c9', 'c7', 'e9', SIMPLE); // Q8, Q9
  c.R('z', 'c8', 150); c.R('z', 'c9', 150);        // R20, R21
  // output pair
  const PNP = { bf: 200 };
  c.PNP('o11', 'c8', 'e11', { ...PNP, simple: true }); c.PNP('o12', 'c9', 'e12', { ...PNP, simple: true });
  c.R('e11', 'p3a', 15); c.R('e12', 'p3b', 15);    // R25, R26
  c.POT('p3a', 'p3b', 'c13', 2.5, 0.5);            // P3 OUTPUT BALANCE
  c.PNP('c13', 'b13', 'e13', PNP);                 // Q13
  c.R('e13', 'p12', 120);                          // R11
  c.R('p12', 'p2a', 820);                          // R12
  const outBias = c.POT('p2a', 'p2b', 'b13', 250, 0.5); // P2 OUTPUT BIAS
  c.R('p2b', '0', 2.7e3);                          // R22
  c.R('o11', 'm6', 680); c.R('o12', 'm6', 680);    // R37, R38
  const m = new Module(c, ['in', 'cv1', 'cv2', 'cv3'], ['o11', 'o12'], sampleRate, oversample);
  // alignment step 6: OUTPUT BIAS trimmed so the signal output sits at 0 V
  m.trimBias = () => {
    let lo = 0, hi = 1;
    for (let i = 0; i < 30; i++) {
      const mid = (lo + hi) / 2; outBias.set(mid); c.dcop();
      if (c.v('o12') > 0) hi = mid; else lo = mid;  // wiper toward p2b -> lower base -> more current
    }
  };
  m.setFixed = pos => fixed.set(pos);
  return m;
}

// ----------------------------------------------------------------------------------------------
// 911 envelope generator (S-trigger input: a switch to ground while the gate is on).
// Q1 attack switch, Q3 release switch, Q2 decay switch into the Q8 sustain follower,
// Q11-Q14 flip-flop, Q9/Q10 Schmitt comparator that ends the attack via C3, Q4 attack kill.
// Output buffer: Q5 JFET (2N4303) follower, Q6 current sink, Q7 emitter follower.
// Solved with backward Euler: the flip-flop's 220 pF and C3 pulse are much faster than
// a sample and trapezoidal would ring on them.
// ----------------------------------------------------------------------------------------------
export function env911(sampleRate, { decimate = 8 } = {}) {
  const c = new Circuit({ integrator: 'be' });
  c.V('p12', '0', 12); c.V('m6', '0', -6);
  const F = { bf: 250 }, Pp = { bf: 200 };
  // trigger and switches
  c.R('p12', 't', 1e3);                                  // R2
  const sw = c.R('t', '0', 1e9);                         // S-trigger contact
  c.R('t', 'b3', 10e3); c.R('b3', '0', 1e3);             // R4, R5
  c.NPN('q3c', 'b3', '0', F);                            // Q3 release switch
  c.R('q3c', 'r8', 100);                                 // R8
  const rel = c.R('r8', 'x', 1e6);                       // R12 RELEASE (T3)
  c.R('t', 'bq', 4.7e3);                                 // R3
  c.PNP('bq', 'j', 'p12', Pp);                           // Q4 attack kill
  c.R('bq', 'b1', 4.7e3); c.R('p12', 'b1', 1e3);         // R7, R6
  c.PNP('q1c', 'b1', 'p12', Pp);                         // Q1 attack switch
  c.R('q1c', 'a1', 100);                                 // R9
  const att = c.R('a1', 'x', 1e6);                       // R10 ATTACK (T1)
  c.C('x', 'm6', 10e-6);                                 // C1 10 µF timing capacitor
  // flip-flop
  c.R('p12', 'j', 470); c.C('j', '0', 220e-12); c.R('j', 'k1', 4.7e3); // R19, 220p, R20
  c.R('p12', 'k2', 4.7e3);                               // R23
  // the flip-flop and comparator get 1 nF junction capacitance (real parts: ~10 pF) so their
  // regenerative flips are a few µs long and solvable by sub-stepping; at envelope time
  // scales this makes no audible difference
  const FF = { ...F, cj: 1e-9 };
  c.NPN('k1', 'b11', '0', FF); c.NPN('k1', 'b12', '0', FF);// Q11, Q12
  c.NPN('k2', 'b13', '0', FF); c.NPN('k2', 'b14', '0', FF);// Q13, Q14
  c.R('k1', 'b13', 47e3); c.R('k2', 'b12', 47e3);        // R21, R22
  c.R('t', 'f', 3.3e3); c.R('f', '0', 3.3e3);            // R1, R40
  c.R('f', 'b14', 3.3e3);                                // R25: trigger off resets to "attack"
  c.R('p', 'b11', 3.3e3);                                // R18: comparator pulse sets "decay"
  // decay switch and sustain level
  c.R('k2', 'b2', 3.3e3);                                // R24
  c.NPN('d2', 'b2', 'es', F);                            // Q2
  const dec = c.R('x', 'r13', 1e6);                      // R11 DECAY (T2)
  c.R('r13', 'd2', 100);                                 // R13
  c.R('p12', 'es', 12e3);                                // R26
  c.PNP('m6', 'wsus', 'es', { ...Pp, simple: true });    // Q8 sustain follower
  // R27 is "selected"; R27 20k + R28 9.5k (trim) puts full sustain at the attack peak
  c.R('p12', 'stop', 29.5e3);
  const sus = c.POT('stop', 'mm', 'wsus', 25e3, 1);      // R29 25k SUSTAIN (pos 0 = full)
  c.R('mm', 'r30', 150); c.R('r30', '0', 125);           // R31, R30 Esus zero trim (centred)
  c.R('mm', 'm6', 2.2e3);                                // R32
  // output buffer
  c.R('x', 'jg', 10e3);
  c.NJF('jd', 'jg', 'js');                               // Q5 2N4303
  c.R('p12', 'jd', 330);                                 // R41
  const zero = c.POT('0', 'm6', 'w14', 2.5e3, 0.5);      // R14 ZERO TRIM
  c.NPN('js', 'w14', 'e6', { ...F, simple: true });      // Q6 current sink
  c.R('e6', 'm6', 750);                                  // R15
  c.NPN('c7', 'js', 'out', { ...F, simple: true });      // Q7
  c.R('p12', 'c7', 100); c.R('out', 'm6', 10e3);         // R16, R17
  // comparator
  c.R('out', 'b9', 33e3); c.R('b9', 'm6', 33e3);         // R33, R34
  c.NPN('c9', 'b9', 'ec', FF); c.NPN('c10', 'c9', 'ec', FF); // Q9, Q10
  c.R('p12', 'c9', 22e3); c.R('p12', 'c10', 22e3);       // R35, R37
  c.R('ec', 'm6', 3.3e3); c.D('0', 'ec');                // R36, CL-1
  c.C('c10', 'pp', 0.01e-6);                             // C3
  c.R('pp', 'p', 3.3e3); c.R('p', '0', 3.3e3);           // R38, R39
  const m = new Module(c, [], ['out', 'x'], sampleRate / decimate, 1);
  // alignment step 1: zero output with no trigger
  let lo = 0, hi = 1;
  for (let i = 0; i < 25; i++) {
    const mid = (lo + hi) / 2; zero.set(mid); c.dcop();
    if (c.v('out') > 0) hi = mid; else lo = mid;          // toward -6 -> less sink current -> higher out
  }
  // map 0..1 knob to the 1M pot (+ the 100Ω end stops), log taper
  const tpot = r => x => r.set(Math.max(1e3, 1e6 * x * x));
  m.setAttack = tpot(att); m.setDecay = tpot(dec); m.setRelease = tpot(rel);
  m.setSustain = x => sus.set(1 - x);
  m.gate = on => sw.set(on ? 10 : 1e9);
  m.decimate = decimate;
  return m;
}
