# Data Sources and Attribution

Frozen mapping release: **19 October 2026**. Package code version: **0.1.0**.
The bundled `manifest.json` records exact SHA256 hashes of upstream final tables
and compact projections. No participant data, cell barcodes or expression
matrices are distributed with the package.

## Primary Sources

- [Ensembl Compara orthology](https://www.ensembl.org/info/genome/compara/homology_method.html):
  original human-mouse stable-ID graph, archived human-fish fallback links,
  and release-116 exact-pair percentage identities. Current scores annotate
  frozen links; they do not replace curated link decisions.
- [ZFIN downloads](https://zfin.org/downloads): direct human-fish curation,
  direct mouse-fish corroboration and identifier cross-references. The frozen
  reports were retrieved on 19 October 2026, advertised 18 October 2026.
- [Human ortholog report](https://zfin.org/downloads/human_orthos.txt)
- [Mouse ortholog report](https://zfin.org/downloads/mouse_orthos.txt)
- [ZFIN-Ensembl associations](https://zfin.org/downloads/ensembl_1_to_1.txt)
- [Alliance of Genome Resources](https://www.alliancegenome.org/),
  [HGNC HCOP](https://www.genenames.org/tools/hcop/) and
  [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/) were supporting sources in
  the full annotated release. Their algorithms overlap; support counts are
  not independent votes or a calibrated posterior probability.

## Frozen ZFIN Checksums

| Report | SHA256 |
| --- | --- |
| human_orthos.txt | 96be5c7eb15c31c03debc08e030d33ef3017329947ca9fa86760281cbd9d98e0 |
| mouse_orthos.txt | 1512ecdce792b2e6e75a0ecf106f6ce8c16d15e280a216c6a41e3ed459265381 |
| Ensembl associations | d3b27ba9f7c43b80542c60ee5a1f78fe21a3a0e387adc6026c3441a5302a23ba |

The compact package data contain selected, resolved pairwise edges, final
human-perspective classes, complete-graph degrees and sequence scores. The
large annotated source tables additionally contain rejected candidates,
publications, identifiers and source-specific scores. A compact projection
must not be presented as the full annotated evidence database.

## Related Software

[mousipy](https://github.com/stefanpeidli/mousipy) is related mouse-to-human
AnnData translation software and inspired the request for a collaborator-
friendly package. heart-orthomap is a separate implementation. In particular,
its count-conservation and no-uppercase-fallback policies differ; consult
each project's documentation before treating outputs as interchangeable.

## Terms and Credits

The MIT license in this repository covers original code/documentation only.
It does not relicense upstream database content. Retain Ensembl, ZFIN,
Alliance, HGNC and NCBI attribution and consult their current data-use terms
when redistributing or combining source data. Database evidence and curator
decisions remain attributable to their original providers.

Terms checked on 5 November 2026:

- ZFIN data and annotations are **CC BY 4.0**. Attribute ZFIN, link the license,
  and identify changes. This package selects and reformats records, adds final
  mapping-policy annotations, and quarantines unresolved/colliding feature IDs;
  it does not distribute an unchanged complete ZFIN database. See the
  [ZFIN terms](https://zfin.atlassian.net/wiki/spaces/general/pages/1942160112/WARRANTY+AND+LIABILITY+DISCLAIMER+OWNERSHIP+AND+LIMITS+ON+USE)
  and [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- Ensembl-generated data are available without restriction; incorporated
  third-party material can retain separate constraints. See the
  [Ensembl data disclaimer](https://useast.ensembl.org/info/about/legal/disclaimer.html).
  Ensembl's software license is distinct from its data terms.
- HGNC data are CC0; attribution is recommended. See the
  [HGNC data license](https://www.genenames.org/about/license/).

The compact wheel includes a data notice with these attributions. Supporting
sources in the full annotated release are not relicensed under the code's MIT
license or presented as independent evidence votes.

Package and heart-focused resource: **Azad Alizada**, **Riley lab**,
**CellTalk HHD (Human Heart Development) consortium**. Biological corrections
should be submitted with stable IDs, exact source evidence, expected mapping
direction and a reproducible example; they must not silently change the frozen
human-mouse baseline.
