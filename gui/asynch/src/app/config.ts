export class ModelConfig {
  uid: number;
  name: string;
  states: string[];
  globalParams: string[];
  forcings: string[];  
}

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

export class AsynchConfig {
  model: number;
  topology: string;
  begin: Date;
  end: Date;
  forcings: Array<string>;
  solver: SolverConfig;  
  scratch_folder: string;
  
  constructor() {
    //this.solver = new SolverConfig();
    
    this.forcings = new Array<string>();
    
    this.scratch_folder = '/tmp';
  } 
  
  addSolver() {
    this.solver = new SolverConfig();
  }
  
  removeSolver() {
    delete this.solver;
  }
}
