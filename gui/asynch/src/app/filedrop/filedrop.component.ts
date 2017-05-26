import { Component, Output, EventEmitter } from '@angular/core';

@Component({
  selector: 'app-filedrop',
  templateUrl: './filedrop.component.html',
  styleUrls: ['./filedrop.component.scss']
})
export class FileDropComponent {

  @Output() public dragFileAccepted: EventEmitter<Object> = new EventEmitter();
  public isHovering: boolean = false;

  private onDragFileOverStart(event) {
    if (!this.isHovering) {
      this.isHovering = true;
    }
    this.preventDefaultAndStopPropagation(event);
    return false;
  };

  private onDragFileOverEnd(event): any {
    this.preventDefaultAndStopPropagation(event);
    return false;
  }

  private onDragFileAccepted(acceptedFile: File): any {
      if (this.dragFileAccepted) {
        this.dragFileAccepted.emit({ file: acceptedFile });
      }
  }

  private onDragFileDrop(event: any): any {
    this.preventDefaultAndStopPropagation(event);
    this.FileSelectHandler(event);
  }

  private preventDefaultAndStopPropagation(event) {
    event.preventDefault();
    event.stopPropagation();
  }

  // file selection handler (can be called from drag, or from a file-requestor select box)
  public FileSelectHandler(e) {
    this.isHovering = false;      // cancel the hover
    var files = e.target.files || e.dataTransfer.files;     // fetch FileList object

    // process all File objects
    for (var i = 0, f; f = files[i]; i++) {
      this.onDragFileAccepted(f);
    }
  }

}
