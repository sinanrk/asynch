import { AsynchPage } from './app.po';

describe('asynch App', () => {
  let page: AsynchPage;

  beforeEach(() => {
    page = new AsynchPage();
  });

  it('should display message saying app works', () => {
    page.navigateTo();
    expect(page.getParagraphText()).toEqual('app works!');
  });
});
