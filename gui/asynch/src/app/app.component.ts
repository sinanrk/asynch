import { Component } from '@angular/core';
import { MdDialog, MdDialogModule } from '@angular/material';
import * as FileSaver from 'file-saver';

import { AsynchConfig } from './config';
import { SaveAsDialogComponent } from './saveas-dialog/saveas-dialog.component';

@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.scss'],
})
export class AppComponent {
  config: AsynchConfig;
  //dialog: MdDialog;
  
  constructor(public dialog: MdDialog) {
    this.config = new AsynchConfig();
  }
  
  onNew() {
    this.config = new AsynchConfig();    
  }
  
  onLoad(file: File) {
    var that = this;
    var reader = new FileReader();

    // Closure to capture the file information.
    reader.onload = (function(theFile) {
      return function(e) {
        that.config = JSON.parse(e.target.result);
      };
    })(file);
    
    reader.readAsText(file);
  }
  
  onSave() {
    //Open a Save As Dialog
    let dialogRef = this.dialog.open(SaveAsDialogComponent);
    dialogRef.afterClosed().subscribe(filename => {
      if (filename) {
            var blob = new Blob([JSON.stringify(this.config)], {type: "text/plain;charset=utf-8"});
            FileSaver.saveAs(blob, filename);
        }
    });
  }
  
  copyToClipboard(element) {
    var selection = window.getSelection();            
    var range = document.createRange();
    range.selectNodeContents(element);
    selection.removeAllRanges();
    selection.addRange(range);
    
    document.execCommand('copy');
  }
}
