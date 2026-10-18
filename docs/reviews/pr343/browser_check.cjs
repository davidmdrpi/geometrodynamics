// Real-browser check of PR #343 HTML: console errors, contact readout, screenshots. Needs a global Playwright.
//   node docs/reviews/pr343/browser_check.cjs docs/visualizations/finite_neck.html <output dir>
const {chromium} = require(require('child_process').execSync('npm root -g').toString().trim() + '/playwright');
const path = require('path');
(async () => {
  const page_path = 'file://' + path.resolve(process.argv[2]);
  const out = process.argv[3];
  const browser = await chromium.launch();
  const errors = [];
  for (const [name, w, motion] of [['desktop', 1200, 'no-preference'], ['phone', 390, 'reduce']]) {
    const page = await browser.newPage({viewport: {width: w, height: 900}, reducedMotion: motion});
    page.on('pageerror', e => errors.push(`${name}: ${e.message}`));
    page.on('console', m => { if (m.type() === 'error') errors.push(`${name}: ${m.text()}`); });
    await page.goto(page_path);
    await page.waitForTimeout(300);
    const status = await page.textContent('#status');
    const play = await page.textContent('#play');
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
    console.log(name, '| play button:', play, '| horizontal overflow px:', overflow, '|', status);
    if (name === 'desktop') {
      await page.click('#play').catch(() => {});
      for (const [tag, js] of [['before', 'tc-.1'], ['at', 'tc'], ['after02', 'tc+.2'], ['after1', 'tc+1']]) {
        await page.evaluate(js => { t = eval(js); playing = false; draw(); }, js);
        await page.locator('#scene').screenshot({path: `${out}/html_${tag}.png`});
      }
      await page.click('#off'); console.log('wave off:', await page.textContent('#status'), '| contact disabled:', await page.isDisabled('#contact'));
    }
    await page.close();
  }
  console.log('errors:', errors.length ? errors : 'none');
  await browser.close();
})();
