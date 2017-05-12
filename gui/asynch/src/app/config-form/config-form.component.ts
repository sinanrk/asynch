import { Component, OnInit } from '@angular/core';
import { ModelConfig, AsynchConfig } from '../config';

@Component({
  selector: 'app-config-form',
  templateUrl: './config-form.component.html',
  styleUrls: ['./config-form.component.scss']
})
export class ConfigFormComponent implements OnInit {

  models: ModelConfig[] = [{
    uid: 190,
    name: 'Constant runoff',
    states: ['q', 'sp', 'ss'],
    globalParams: ['v_r', 'lambda_1', 'lambda_2', 'RC', 'v_h', 'v_g'],
    forcings: ['Precipitation', 'Evaporation']
  }, {
    uid: 254,
    name: 'Top Layer Hillsope',
    states: ['q', 'sp', 'st', 'ss'],
    globalParams: ['v_0', 'lambda_1', 'lambda_2', 'v_h', 'k_3', 'k_I_factor', 'h_b', 'S_L', 'A', 'B', 'exponent', 'v_B'],
    forcings: ['Precipitation', 'Evaporation', 'Reservoirs']
  }];
  
  config: AsynchConfig;
  selectedModel: ModelConfig;
  
  constructor() {
    this.config = new AsynchConfig();
  }

  ngOnInit() {
  }
  
  onModelChange(uid): void {
    this.selectedModel = this.models.find((model) => {return model.uid === uid;})
  }
  
  copyToClipboard(element) {
    var selection = window.getSelection();            
    var range = document.createRange();
    range.selectNodeContents(element);
    selection.removeAllRanges();
    selection.addRange(range);
    
    document.execCommand('copy');
  }
  
  useDefaultSolverConfig(useDefault: boolean) {
    if (useDefault) {
      this.config.removeSolver();
    } else  {
      this.config.addSolver();
    }
      
  }
}
