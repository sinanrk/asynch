import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { SaveAsDialogComponent } from './saveas-dialog.component';

describe('SaveAsDialogComponent', () => {
  let component: SaveAsDialogComponent;
  let fixture: ComponentFixture<SaveAsDialogComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ SaveAsDialogComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SaveAsDialogComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
