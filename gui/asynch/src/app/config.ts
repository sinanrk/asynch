import { ModelMeta } from './models';

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

export class ForcingConfig {
  filename: string;
  step: number;
}

//export class DbTimeserieConfig extends TimeserieConfig {
//  filename: string;
//  step: number;
//}

export class ForcingStateConfig {
  filename: string;
  index: number;
}

export class ForcingsConfig {
  timeseries: Array<string | ForcingConfig>;
  state: ForcingStateConfig;
  dams: string;
}


export class TimeserieConfig {
  
}

interface PeakConfig {
  
}

interface SnapshotConfig {
  
}

interface OutputsConfig {
  timeseries?: TimeserieConfig;
  peaks?: PeakConfig;
  snapshots?: SnapshotConfig;
}


export class AsynchConfig {
  model: number;
  topology: string;
  begin: Date;
  end: Date;
  globalParams: Array<number>;
  scratch_folder: string;
  forcings?: ForcingsConfig;
  solver?: SolverConfig;  
  outputs?: OutputsConfig;
  
  constructor() {
    this.scratch_folder = '/tmp';
    this.outputs = {
      timeseries: {},
      peaks: {},
      snapshots: {}
    };
  } 
  
  init(meta: ModelMeta) {
    this.model = meta.uid;
    
    //Init forcings
    this.forcings = new ForcingsConfig();
    this.forcings.timeseries = new Array<ForcingConfig>(meta.forcings.length);
    this.globalParams = meta.globalParamsDefault;
  }
  
  addSolver() {
    this.solver = new SolverConfig();
  }
  
  removeSolver() {
    delete this.solver;
  }
}
