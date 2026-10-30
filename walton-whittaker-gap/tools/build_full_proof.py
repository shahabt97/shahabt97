"""Assemble FULL_PROOF.md from docs/FULL_PROOF.template.md.

Placeholders in the template:
    {{FILE:path}}   the source file `path` verbatim, in a fenced code block
    {{LOG:name}}    logs/name.txt verbatim, in a fenced code block
    {{SHA:path}}    the sha256 of the file `path`
usage: python -m tools.build_full_proof          (from the package root)
"""
import hashlib
import os
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def fence(text, lang):
    ticks = max([3] + [len(m) + 1 for m in re.findall(r'`+', text)])
    return '`' * ticks + lang + '\n' + text.rstrip('\n') + '\n' + '`' * ticks


def expand(m):
    kind, arg = m.group(1), m.group(2).strip()
    if kind == 'FILE':
        with open(os.path.join(ROOT, arg), encoding='utf-8') as fh:
            return fence(fh.read(), 'python' if arg.endswith('.py') else '')
    if kind == 'LOG':
        with open(os.path.join(ROOT, 'logs', arg + '.txt'), encoding='utf-8') as fh:
            return fence(fh.read(), 'text')
    with open(os.path.join(ROOT, arg), 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    with open(os.path.join(ROOT, 'docs', 'FULL_PROOF.template.md'), encoding='utf-8') as fh:
        template = fh.read()
    out = re.sub(r'\{\{(FILE|LOG|SHA):([^}]+)\}\}', expand, template)
    assert not re.search(r'\{\{(FILE|LOG|SHA):', out), "unresolved placeholder"
    path = os.path.join(ROOT, 'FULL_PROOF.md')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(out)
    print(f"wrote {os.path.normpath(path)}: {len(out)} characters, {out.count(chr(10)) + 1} lines")


if __name__ == '__main__':
    main()
