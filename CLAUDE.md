# Strict Development Workflow

To guarantee code quality and transparent change tracking, this repository enforces a strict workflow. Direct pushes to `master` are **strictly prohibited**.

## 1. Zero-Warning Policy (`-D warnings`)
Your local environment must treat all warnings as errors. Do not commit code with warnings.
* **Requirement:** Run your compiler/linter with the equivalent of `-D warnings` or `-Werror`.
* **Action:** Fix all warnings locally before pushing.

## 2. Mandatory Branching
All modifications must happen on an isolated branch to make changes against the `master` branch obvious during review.
```bash
git checkout -b feature/your-branch-name
```
## 3. GitHub CI/CD & Pull Requests
To integrate your code:

Open a PR: Push your branch and open a Pull Request targeting master.

Pass CI/CD: Our GitHub Actions pipeline will enforce the -D warnings rule. If the CI detects a single warning, the build will fail and the PR will be blocked.

Review & Merge: Code is only merged into master after passing CI and receiving maintainer approval.