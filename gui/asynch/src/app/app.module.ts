import 'hammerjs';

import { BrowserModule } from '@angular/platform-browser';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { NgModule } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { HttpModule } from '@angular/http';
import { FlexLayoutModule } from '@angular/flex-layout';

import { MdToolbarModule, MdButtonModule, MdCheckboxModule, MdInputModule, MdSelectModule, MdCardModule, MdSidenavModule, MdListModule, MdDialogModule, MdDatepickerModule, MdNativeDateModule } from '@angular/material';

import { AppComponent } from './app.component';
import { ConfigFormComponent } from './config-form/config-form.component';
import { SaveAsDialogComponent } from './saveas-dialog/saveas-dialog.component';
import { FileExtValidatorDirective } from './file-ext-valid.directive';
import { FileDropDirective } from './filedrop.directive';
import { FileExtPipe } from './file-ext.pipe';


@NgModule({
  declarations: [
    AppComponent,
    ConfigFormComponent,
    FileExtValidatorDirective,
    FileExtPipe,
    SaveAsDialogComponent,
    FileDropDirective
  ],
  entryComponents: [
    SaveAsDialogComponent
  ],
  imports: [
    BrowserModule,
    BrowserAnimationsModule,
    FormsModule,
    HttpModule,
    FlexLayoutModule,
    MdToolbarModule,
    MdButtonModule,
    MdCheckboxModule,
    MdInputModule,
    MdSelectModule,
    MdCardModule,
    MdSidenavModule,
    MdListModule,
    MdDialogModule,
    MdDatepickerModule,
    MdNativeDateModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
