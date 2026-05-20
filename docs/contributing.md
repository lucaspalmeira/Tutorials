# Contributing Tutorials

Keep every tutorial in English and write it as a reproducible protocol, not only as a list of commands.

## Recommended Structure

Each tutorial should include:

- goal
- when to use the workflow
- required software
- expected input files
- command sequence
- expected outputs
- validation checks
- troubleshooting notes
- adaptation notes for other systems or clusters

## Add a New Tutorial Page

1. Create a Markdown file in `docs/tutorials/`.
2. Use a short, URL-friendly file name such as `new-workflow.md`.
3. Add the page to the `nav` section in `mkdocs.yml`.
4. Run `mkdocs build --strict`.
5. Commit the page together with any scripts or example files it requires.

## Using README Files as Tutorials

If a workflow is developed in its own folder, write that folder's `README.md` as the source tutorial in English. Then copy the finalized tutorial text into a page under `docs/tutorials/` or keep the docs page as the canonical version and link back to the workflow folder.

The published site should only expose finished English documentation. Draft notes, local paths, and cluster-specific secrets should stay out of the public pages.
