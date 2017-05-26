import { Component } from '@angular/core';
import { MdDialogRef } from '@angular/material';

@Component({
  templateUrl: './saveas-dialog.component.html',
  styleUrls: ['./saveas-dialog.component.scss']
})
export class SaveAsDialogComponent {

  constructor(public dialogRef: MdDialogRef<SaveAsDialogComponent>) {}
}
