export class ModelMeta {
  uid: number;
  name: string;
  states: string[];
  globalParams: string[];
  globalParamsDefault: number[];
  forcings: string[];
  hasDam: boolean;
}

export const models: Array<ModelMeta> = [{
  uid: 190,
  name: "Constant runoff",
  states: ["q", "sp", "ss"],
  globalParams: ["v_r", "lambda_1", "lambda_2", "RC", "v_h", "v_g"],
  globalParamsDefault: [0.33, 0.20, -0.1, 0.33, 0.1, 2.2917e-5],
  forcings: ["Precipitation", "Evaporation"],
  hasDam: false
  }, {
  uid: 254,
  name: "Top Layer Hillsope",
  states: ["q", "sp", "st", "ss"],
  globalParams: ["v_0", "lambda_1", "lambda_2", "v_h", "k_3", "k_I_factor", "h_b", "S_L", "A", "B", "exponent", "v_B"],
  globalParamsDefault: [0.33, 0.20, -0.1, 0.02, 2.0425e-6, 0.02, 0.5, 0.10, 0.0, 99.0, 3.0, 0.75],
  forcings: ["Precipitation", "Evaporation", "Reservoirs"],
  hasDam: false
}]
