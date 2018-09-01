# Developer instructions

## Automated tests

### During development

[`checks.sh`](./checks.sh) executes tests which should be run when developing.
You should run these during development as frequently as suits you.

It is advised to be aware of the state of the tests when making a commit, which is best automated by setting up a pre-commit hook.
When you clone this repo, please make sure to set this script up as your pre-commit hook:

```sh
cd ../.git/hooks
ln -s pre-commit ../../dev/pre-commit.sh
```

You may choose to ignore these checks during development.

### When submitting PR

When you submit a PR, you should pay particular attention to the results of the tests.
To ensure you don't mistakenly push a branch with failing tests, you're advised to set up a pre-push hook:

```sh
cd ../.git/hooks
ln -s pre-push ../../dev/pre-push.sh
```

Unit tests are required to pass.

`mypy` typings tend to be flaky in my experience and thus are not all required to pass.
They are tested against a snapshot of the existing errors instead.
Note that line number changes are ignored in dev mode, so you may see changes which weren't flagged before.
You should update the snapshot with new line numbers before pushing.

If your changes diverge from the snapshot, you should update and commit it to the repo, so that the new errors can be discussed in the PR:

```sh
./update-mypy.sh
git add mypy_snapshot.out
```
