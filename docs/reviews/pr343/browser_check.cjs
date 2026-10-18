// Headless-Chromium check of an HTML page: console errors, horizontal overflow, screenshots.
//   node browser_check.cjs <page.html> <output dir> [js to run before the screenshot]...
const {chromium} = require(require('child_process').execSync('npm root -g').toString().trim() + '/playwright');
const path = require('path');
(async () => {
  const [page_path, out, ...scripts] = process.argv.slice(2);
  const browser = await chromium.launch(); const errors = [];
  for (const [name, w] of [['desktop', 1200], ['phone', 390]]) {
    const page = await browser.newPage({viewport: {width: w, height: 900}, reducedMotion: 'reduce'});
    page.on('pageerror', e => errors.push(`${name}: ${e.message}`));
    page.on('console', m => { if (m.type() === 'error') errors.push(`${name}: ${m.text()}`); });
    await page.goto('file://' + path.resolve(page_path)); await page.waitForTimeout(400);
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth - innerWidth);
    const status = await page.evaluate(() => (document.getElementById('status') || {}).textContent || '');
    console.log(`${name}: overflow ${overflow}px | ${status.slice(0, 160)}`);
    if (name === 'desktop') for (const [k, js] of scripts.entries()) {
      await page.evaluate(js); await page.waitForTimeout(200);
      await page.screenshot({path: `${out}/shot_${k}.png`, fullPage: true});
    }
    await page.close();
  }
  console.log('errors:', errors.length ? errors : 'none'); await browser.close();
})();
