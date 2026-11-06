// Synth modules, each a real schematic expressed as a netlist.
// Every module: build(circuit) wires parts; inputs are voltage sources you set
// per sample; outputs are node voltages you read back.
import { Circuit } from './solver.js';

const VCC = 12, VEE = -12;

function rails(c) {
  c.V('vcc', '0', VCC);
  c.V('vee', '0', VEE);
}

// Runs a circuit at sampleRate * oversample, holds inputs across sub-steps,
// and averages the outputs back down (a crude but safe decimator).
export class Module {
  // init: 'dc' solves the DC operating point first; 'zero' starts with every capacitor
  // discharged (needed for oscillators, whose DC point is an unstable balance).
  constructor(circuit, inputs, outputs, sampleRate, oversample = 4, { init = 'dc' } = {}) {
    Object.assign(this, { circuit, inputs, outputs, oversample });
    circuit.probe(...outputs);
    circuit.finalize(sampleRate * oversample);
    if (init === 'dc') circuit.dcop();
    this.out = Object.fromEntries(outputs.map(o => [o, 0]));
    this.idx = outputs.map(o => circuit.nodeIndex.get(o));
    this.acc = new Float64Array(outputs.length);
  }
  process(inValues = {}) {
    const c = this.circuit, n = this.idx.length, acc = this.acc;
    for (const k of this.inputs) if (k in inValues) c.set(k, inValues[k]);
    acc.fill(0);
    for (let s = 0; s < this.oversample; s++) {
      c.step();
      for (let j = 0; j < n; j++) acc[j] += c.vNode(this.idx[j]);
    }
    for (let j = 0; j < n; j++) this.out[this.outputs[j]] = acc[j] / this.oversample;
    return this.out;
  }
}

// 1) MIXER: three-input inverting summer, unity gain per channel.
//    in1/in2/in3 -100k-> sum node -> op-amp(-) ; 100k feedback ; (+) to ground.
export function mixer(sampleRate, oversample = 1) {
  const c = new Circuit();
  rails(c);
  for (const k of ['in1', 'in2', 'in3']) {
    c.V(k + '_src', '0', 0, k);
    c.R(k + '_src', 'sum', 100e3);
  }
  c.R('sum', 'out', 100e3);
  c.OPAMP('0', 'sum', 'out');
  return new Module(c, ['in1', 'in2', 'in3'], ['out'], sampleRate, oversample);
}

// 2) LFO: op-amp integrator + Schmitt trigger (the textbook triangle/square pair).
//    f = R3 / (4 * Rrate * C * R2); triangle peak = Vsat * R2 / R3.
//    With R2 = 10k, R3 = 20k, C = 1uF, Rrate = 100k -> 5 Hz, ~ ±5.25 V triangle.
export function lfo(sampleRate, { rate = 100e3, cap = 1e-6, oversample = 4 } = {}) {
  const c = new Circuit();
  rails(c);
  // integrator: square -> Rrate -> int(-), C from int(-) to tri, (+) grounded
  const rRate = c.R('sq', 'int_n', rate);
  c.C('int_n', 'tri', cap);
  c.OPAMP('0', 'int_n', 'tri');
  // Schmitt trigger: (+) sees tri through R2 and square through R3 (positive feedback)
  c.R('tri', 'sch_p', 10e3);
  c.R('sq', 'sch_p', 20e3);
  // 1 mV input offset (real op-amps have one) breaks the symmetric start-up state
  c.OPAMP('sch_p', '0', 'sq', { vos: 1e-3, gbw: 90e3, slew: 0.5e6 }); // slow pole: see OPAMP gbw note
  const m = new Module(c, [], ['tri', 'sq'], sampleRate, oversample, { init: 'zero' });
  m.setRate = ohms => rRate.set(ohms);
  return m;
}

// 3) VCA: discrete BJT differential pair (Moog 902 family idea).
//    - Signal is attenuated ~1/101 (5 V in -> ~50 mV, where tanh saturation starts) into Q1's base, Q2's base is grounded.
//    - Tail current: op-amp inverts CV (0..10 V -> 0..-10 V) and a 47k resistor
//      from the joined emitters to that node sinks I_tail ~ (CV - 0.65 V) / 47k.
//    - Collectors have 10k loads; an op-amp difference amp (gain 4.7) removes
//      the CV bleed-through (common mode) and gives the output.
//    Gain is roughly proportional to I_tail; large inputs soft-clip as tanh.
export function vca(sampleRate, oversample = 4) {
  const c = new Circuit();
  rails(c);
  c.V('in_src', '0', 0, 'in');
  c.V('cv_src', '0', 0, 'cv');
  // input attenuator
  c.R('in_src', 'b1', 100e3);
  c.R('b1', '0', 1e3);
  c.R('b2', '0', 1e3);
  // differential pair
  c.R('vcc', 'c1', 10e3);
  c.R('vcc', 'c2', 10e3);
  c.NPN('c1', 'b1', 'e');
  c.NPN('c2', 'b2', 'e');
  // CV inverter + tail resistor
  c.R('cv_src', 'inv_n', 100e3);
  c.R('inv_n', 'tail_ref', 100e3);
  c.OPAMP('0', 'inv_n', 'tail_ref');
  c.R('e', 'tail_ref', 47e3);
  // difference amplifier: out = 4.7 * (c2 - c1)
  c.R('c1', 'da_n', 10e3);
  c.R('da_n', 'out', 47e3);
  c.R('c2', 'da_p', 10e3);
  c.R('da_p', '0', 47e3);
  c.OPAMP('da_p', 'da_n', 'out');
  return new Module(c, ['in', 'cv'], ['out'], sampleRate, oversample);
}
