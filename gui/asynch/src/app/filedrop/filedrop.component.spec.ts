import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { FiledropComponent } from './filedrop.component';

describe('FiledropComponent', () => {
  let component: FiledropComponent;
  let fixture: ComponentFixture<FiledropComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ FiledropComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(FiledropComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should be created', () => {
    expect(component).toBeTruthy();
  });
});
