# Research contributions

This repository maintains work on the 547/727 and 571/599 theta prime races. Start at [the results index](RESULTS.md).

Extend the existing result when an argument, bound or certificate improves it. Give a distinct result in the same family a directory under `results/` and an index entry. Record its hypotheses, its verification command and which parts the proof checker or numerical certificate covers. Keep evidence and independent checks identifiable.

Use the following maintained homes when the mathematical question belongs to one of these families. Distinct questions can have their own project.

- [547/727 and 571/599 theta prime races](https://github.com/leecrossley/theta-host-split-inert-race)
- [Real quadratic Mertens crossings](https://github.com/leecrossley/quadratic-mertens-first-crossings)
- [Ordinary-prime windows](https://github.com/leecrossley/prime-window-all-widths)
- [Window criteria for the Riemann hypothesis](https://github.com/leecrossley/prime-window-rh-criterion)
- [Eisenstein prime-ideal race](https://github.com/leecrossley/eisenstein-prime-window-transition)
- [Integrated prime-counting errors in progressions](https://github.com/leecrossley/integrated-prime-errors-mod-11)

Keep original release tags, assets, commit links and Zenodo identifiers available. Cite a result through its original publication record. A new release in a family repository can have different contents from an earlier paper, so identify what was released and use the appropriate citation.

When a paper’s version changes, remove the previous edition’s DOI from its citation file before tagging the release. For a paper released from a collection, set its title, description and version in the root `.zenodo.json` on the release tag. After publication, check that the DOI resolves to that tag and add it to the paper’s `CITATION.cff` and the [edition index](RESULTS.md#edition-dois). If several papers share a release, give the shared archive DOI and the path to each paper. [Zenodo release metadata](https://help.zenodo.org/docs/github/describe-software/zenodo-json/)

## Release metadata

For a single-paper software release, work from the paper revision in a detached worktree. Set `paper` to the result directory name and `tag` to the intended tag name. With an abstract in the paper’s citation file, this Bash pipeline prints the Zenodo metadata, including the creators and licence. For an older citation file without an abstract, take the description from the paper. [Citation File Format converter](https://github.com/citation-file-format/cffconvert) [Zenodo metadata](https://help.zenodo.org/docs/github/describe-software/zenodo-json/)

```bash
set -o pipefail
cffconvert --infile "results/${paper:?}/CITATION.cff" --format zenodo |
  jq -e --arg tag "${tag:?}" '
    select((.description | length) > 0)
    | .version = $tag
    | .description = ("<p>" + (.description | @html) + "</p>")
    | .license = (.license.id | ascii_downcase)
    | .upload_type = "software"'
```

Save the output as `.zenodo.json` at the root of the release worktree and commit it there. Place the tag on that commit. Before publishing the GitHub release, `git show "$tag:.zenodo.json"` displays the title, description and version Zenodo will read. [Dyadic paper example](https://github.com/leecrossley/theta-host-split-inert-race/blob/consecutive-doubling-patterns-v1.4.0/.zenodo.json)
