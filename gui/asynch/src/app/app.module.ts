import 'hammerjs';

import { BrowserModule } from '@angular/platform-browser';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { NgModule } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { HttpModule } from '@angular/http';
import { FlexLayoutModule } from '@angular/flex-layout';

import { MdToolbarModule, MdButtonModule, MdCheckboxModule, MdInputModule, MdSelectModule, MdCardModule, MdSidenavModule, MdListModule } from '@angular/material';

import { AppComponent } from './app.component';
import { ConfigFormComponent } from './config-form/config-form.component';
import { FileExtValidatorDirective } from './file-ext-valid.directive';

@NgModule({
  declarations: [
    AppComponent,
    ConfigFormComponent,
    FileExtValidatorDirective
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
    MdListModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
