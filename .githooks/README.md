# Git hooks

This directory travels **inside** the repository, so the hooks are version-controlled
and shared with every collaborator. Git does not run hooks from a tracked directory
automatically, though — each person must point git at this folder once after cloning.

## What's here

- **`pre-commit`** — blocks any staged file **≥ 50 MB** from being committed. Large
  files bloat git history permanently (git keeps every blob forever, even after a
  later delete) and can hit GitHub's limits (50 MB warning / 100 MB hard reject).

## Install (run once after `git clone`)

From the repository root:

```sh
git config core.hooksPath .githooks
```

That tells git to look for hooks in `.githooks/` (this directory) instead of the
default, untracked `.git/hooks/`. The setting is stored in your **local**
`.git/config`, so every collaborator runs it once on their own machine.

If the hook ever doesn't fire, make sure it's executable:

```sh
chmod +x .githooks/pre-commit
```

## Usage

Once installed, the hook runs on every `git commit`. If a staged file is too big,
the commit is blocked with a message listing the offending file(s).

- **Raise the limit for one commit:** `MAX_MB=100 git commit ...`
- **Keep the file out of git:** add its path/pattern to `.gitignore`, then
  `git restore --staged <file>`
- **Bypass the hook once (use sparingly):** `git commit --no-verify`

## Verify it works

```sh
# create a 60 MB file, stage it, and try to commit
dd if=/dev/zero of=bigfile.bin bs=1m count=60
git add bigfile.bin
git commit -m "test"      # should be BLOCKED
git restore --staged bigfile.bin && rm bigfile.bin
```
