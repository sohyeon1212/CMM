# Clean Implementation Policy

CMM is implemented as a new codebase.

## Rules

- Do not copy source files, function bodies, UI forms, resource files, icons, docs, test
  fixtures, or project archives from a legacy application.
- Implement behavior from product requirements and public documentation for third-party
  libraries.
- Keep public algorithms behind small service interfaces so the desktop UI does not
  inherit old menu or dialog structure.
- Use new project data structures and file formats.
- Keep provenance notes outside the CMM repository when they are only needed for internal
  migration planning.

## Allowed References

- public package APIs and documentation
- mathematical descriptions from papers or official package examples
- user-authored requirements and expected workflow descriptions
- newly created toy models for tests

## Audit Checklist

- No legacy package names in source paths, entry points, config files, or docs
- No legacy project extensions as primary file formats
- No copied screenshots, icons, SVG maps, or bundled examples
- No ported tests that encode old file layout or UI labels
- New About dialog, settings directory, project archive manifest, and CLI command
