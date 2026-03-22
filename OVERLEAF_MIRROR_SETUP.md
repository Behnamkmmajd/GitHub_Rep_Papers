# Overleaf Mirror Setup

This repository includes an automated GitHub Actions mirror for the LaTeX paper only.

Source folder mirrored:

- project1_linear_datadriven/Latex

The workflow file is:

- .github/workflows/mirror-latex-to-overleaf.yml

Setup required on GitHub:

1. Create a second GitHub repository for the Overleaf-facing project.
2. Keep that target repository empty, or let the first mirror push overwrite its default branch.
3. In this repository, open Settings -> Secrets and variables -> Actions.
4. Add a repository variable named OVERLEAF_MIRROR_REPO with the value owner/repo for the target repository.
5. Add a repository secret named OVERLEAF_MIRROR_TOKEN containing a fine-grained personal access token with Contents: Read and write access to the target repository.
6. Optionally add OVERLEAF_MIRROR_BRANCH if the target branch should be something other than main.

Behavior:

- Every push to main that changes project1_linear_datadriven/Latex triggers the workflow.
- The workflow extracts only that folder with git subtree split.
- It then force-pushes the extracted history to the target repository branch.

Notes:

- The default GitHub Actions token cannot push to a different repository, which is why a personal access token is required.
- Overleaf should be connected to the second repository, not to this monorepo.
