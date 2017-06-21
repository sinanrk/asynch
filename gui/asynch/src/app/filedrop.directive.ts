import { Directive, ElementRef, Output, EventEmitter, HostListener } from '@angular/core';

@Directive({
  selector: '[filedrop]'
//  templateUrl: './filedrop.component.html',
//  styleUrls: ['./filedrop.component.scss']
})
export class FileDropDirective {

  @Output() public dragFileAccepted: EventEmitter<Object> = new EventEmitter();
  public isHovering: boolean = false;
  
//  constructor(el: ElementRef, renderer: Renderer) {
//    renderer.listen(el.nativeElement, 'dragover', this.onDragFileOverStart);
//    renderer.listen(el.nativeElement, 'dragleave', this.onDragFileOverEnd);
//    renderer.listen(el.nativeElement, 'drop', this.onDragFileDrop);
//  }

  @HostListener('dragover', ['$event']) onDragFileOverStart(event) {
    if (!this.isHovering) {
      this.isHovering = true;
    }
    this.preventDefaultAndStopPropagation(event);
    return false;
  };

  @HostListener('dragleave', ['$event']) onDragFileOverEnd(event): any {
    this.preventDefaultAndStopPropagation(event);
    return false;
  }

  @HostListener('drop', ['$event']) onDragFileDrop(event: any): any {
    this.preventDefaultAndStopPropagation(event);
    this.FileSelectHandler(event);
  }
  
  private onDragFileAccepted(acceptedFile: File): any {
      if (this.dragFileAccepted) {
        this.dragFileAccepted.emit({ file: acceptedFile });
      }
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
