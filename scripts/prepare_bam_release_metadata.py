"""Set the actual release date before committing/tagging; or validate tag metadata."""
import argparse
from datetime import date
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def prepare(root=ROOT, release_date=None):
    root=Path(root)
    version=json.loads((root/'docs/releases/bam-foundations-evidence.json').read_text())['release']
    cff=root/'CITATION.cff'
    changelog=root/'CHANGELOG.md'
    if release_date is not None:
        day=date.fromisoformat(release_date).isoformat()
        text=re.sub(r'^date-released:.*\n?', '', cff.read_text(), flags=re.M)
        citation=text.rstrip()+f'\ndate-released: {day}\n'
        text,n=re.subn(r'^## '+re.escape(version)+r' — .*$', f'## {version} — {day}',changelog.read_text(),count=1,flags=re.M)
        if n!=1:raise ValueError('Release changelog entry missing')
        cff.write_text(citation)
        changelog.write_text(text)
    matches=re.findall(r'^date-released: (\d{4}-\d{2}-\d{2})$',cff.read_text(),flags=re.M)
    if len(matches)!=1:raise ValueError('Set the actual release date before tagging with --date YYYY-MM-DD')
    day=date.fromisoformat(matches[0]).isoformat()
    if f'## {version} — {day}' not in changelog.read_text().splitlines():
        raise ValueError('Citation and changelog release dates differ')
    return day


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--date', help='Actual release date, set before committing and tagging')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    print(prepare(release_date=args.date))
