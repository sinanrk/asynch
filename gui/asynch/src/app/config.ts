import { ModelMeta } from 'app/models';

export class ErrorCtlConfig {
  facmin: number;
  facmax: number;
  fac: number;
}

export class BuffersConfig {
  num_step: number;
  num_transfer: number;
  num_discont: number;
  
  constructor() {
    this.num_step = 30;
    this.num_transfer = 10;
    this.num_discont = 30;
  }
}

export class SolverConfig {
  method: number;
  error_ctl: ErrorCtlConfig;
  buffers: BuffersConfig;
  tolerances: Array<number>;  
  
  constructor() {
    this.error_ctl = new ErrorCtlConfig();
    this.buffers = new BuffersConfig();
    this.tolerances = new Array<number>(1e-3, 1e-6, 1e-3, 1e-6);
  }
}

export class TimeserieConfig {
  filename: string;
}

export class DbTimeserieConfig extends TimeserieConfig {
  step: number;
}

export class ForcingStateConfig {
  filename: string;
  index: number;
}

export class ForcingsConfig {
  timeseries: Array<TimeserieConfig>;
  state: ForcingStateConfig;
}

export class AsynchConfig {
  model: number;
  topology: string;
  begin: Date;
  end: Date;
  globalParams: Array<number>;
  dams: string;
  forcings: ForcingsConfig;
  solver: SolverConfig;  
  scratch_folder: string;
  outputs: any;
  
  constructor() {
    this.scratch_folder = '/tmp';
    this.outputs = {
      timeseries: {},
      peaks: {},
      snaphosts: {}
    };
  } 
  
  init(meta: ModelMeta) {
    this.model = meta.uid;
    
    //Init forcings
    this.forcings = new ForcingsConfig();
    this.forcings.timeseries = new Array<TimeserieConfig>(meta.forcings.length);
    this.globalParams = meta.globalParamsDefault;
  }
  
  addSolver() {
    this.solver = new SolverConfig();
  }
  
  removeSolver() {
    delete this.solver;
  }
}
