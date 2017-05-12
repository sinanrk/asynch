import { Directive, Input, OnChanges, SimpleChanges} from '@angular/core';
import { NG_VALIDATORS, AbstractControl, Validator, ValidatorFn, Validators } from '@angular/forms';

export function fileExtValidator(fileExt: string[]): ValidatorFn {
  return (control: AbstractControl): {[key: string]: any} => {
    const filename = control.value || '';
    const ext = filename.substr((~-filename.lastIndexOf('.') >>> 0) + 2);
    const no = fileExt.some((e) => {return e === ext;} );
    return (!no) ? {'fileExt': {filename}} : null;
  };
}

@Directive({
  selector: '[fileExt]',
  providers: [{
    provide: NG_VALIDATORS,
    useExisting: FileExtValidatorDirective,
    multi: true }]
})
export class FileExtValidatorDirective implements Validator, OnChanges  {
  @Input() fileExt: string;
  private valFn = Validators.nullValidator;

//  constructor(public model: NgModel) {
//    console.error('init')
//  }
  
  ngOnChanges(changes: SimpleChanges): void {
    const change = changes['fileExt'];
    if (change) {
      const val: string = change.currentValue;
      const ar = val.split(',');
      this.valFn = fileExtValidator(ar);
    } else {
      this.valFn = Validators.nullValidator;
    }
  }
   
  validate(control: AbstractControl): {[key: string]: any} {
    return this.valFn(control);
  }

  static getFileExtension(filename: string) {
    return filename.substr((~-filename.lastIndexOf('.') >>> 0) + 2);
  }
  
}
