# Bundled third-party resources

## `e_coli_core.Core metabolism.json`

An Escher pathway map of the *Escherichia coli* core metabolic network, redistributed
unmodified and used to lay out CMM's flux maps.

| | |
|---|---|
| Retrieved from | <https://escher.github.io/1-0-0/6/maps/Escherichia%20coli/e_coli_core.Core%20metabolism.json> |
| Retrieved | 2026-08-20 |
| SHA-256 | `2caf3c036ff93f487a12c601adeb74353c174c96fc144730cc528b304b1141bc` |
| Schema | `https://escher.github.io/escher/jsonschema/1-0-0#` |
| Contents | 95 reactions, 462 nodes |

### Which licence applies, and how that was established

That URL is served by GitHub Pages from the repository
[`escher/escher.github.io`](https://github.com/escher/escher.github.io), where the file is
tracked at `1-0-0/6/maps/Escherichia coli/e_coli_core.Core metabolism.json`. The copy in that
repository is **byte-identical** to the copy here — the SHA-256 above matches it — so the file
CMM redistributes is that repository's file, and that repository's licence is the one that
governs it:

> The MIT License (MIT). This software is Copyright © 2019 The Regents of the University of
> California. All Rights Reserved.

MIT permits redistribution provided the copyright notice and permission notice travel with the
copy, which is what this file is for. (GitHub's licence detector reports "Other" for that
repository only because the extra `This software is Copyright …` line stops its automatic
matcher; the text below that line is verbatim MIT.)

The Escher *application* repository, `zakandrewking/escher`, is MIT as well but carries a
different copyright year (2015) and **does not contain these map files**, so it is not the
source of this one and its notice is not the notice to reproduce.

**Do not re-source this file from BiGG Models.** BiGG serves the same map under different
terms — free for educational, research and non-profit use, with commercial use requiring a
written agreement with UC San Diego — which are not compatible with redistributing it inside an
MIT-licensed package. The provenance recorded here is what makes CMM's redistribution sound.

### Citation

The map is a published, hand-drawn layout. Cite its paper in any work that uses it:

> King ZA, Dräger A, Ebrahim A, Sonnenschein N, Lewis NE, Palsson BØ. Escher: a web
> application for building, sharing, and embedding data-rich visualizations of biological
> pathways. *PLoS Comput Biol* 2015;11(8):e1004321. doi:10.1371/journal.pcbi.1004321

The file is redistributed **byte-for-byte as retrieved**; its `map_description` still carries
the upstream "Last Modified Fri Dec 05 2014" stamp. CMM neither edits nor regenerates it.
CMM's own licence (MIT) is in `LICENSE` at the repository root and covers CMM's code only.
