# Overleaf Mirror Setup

This repository includes an automated GitHub Actions mirror for the paper folders only.

Source folders mirrored into the Overleaf repository root:

- Project1/ from project1_linear_datadriven/Latex
- Project2/ from Project2_LCS/Latex

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

- Every push to main that changes either source LaTeX folder triggers the workflow.
- The workflow assembles a clean export repository containing Project1/ and Project2/ as top-level folders.
- It then force-pushes that assembled export to the target repository branch.

Notes:

- The default GitHub Actions token cannot push to a different repository, which is why a personal access token is required.
- Overleaf should be connected to the second repository, not to this monorepo.
