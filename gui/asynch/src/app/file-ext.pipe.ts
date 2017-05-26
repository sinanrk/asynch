import { Pipe, PipeTransform } from '@angular/core';

@Pipe({
  name: 'fileExt'
})
export class FileExtPipe implements PipeTransform {

  transform(filename: string, ext: string): boolean {
    if (filename) {
      let cur = filename.substr((~-filename.lastIndexOf('.') >>> 0) + 2);
      return cur === ext;
    } else {
      return false;
    }
  }

}
