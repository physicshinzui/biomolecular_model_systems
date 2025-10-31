# Biomolecular Model Systems

Curated protein and peptide model systems for benchmarking molecular dynamics (MD) or other structure-based workflows. Each entry is supplied as a PDB file with capped termini (where appropriate) and minimal preprocessing so you can drop the structures directly into your validation pipeline.

## What's Included

- Curated proteins and peptides frequently used for MD force-field validation
- Multiple conformational states (folded, extended, mutant variants) for several protein targets
- Top-level grouping into `proteins/` and `peptides/` for easier scripting and automation

## Directory Guide

### Proteins

| Path | Notes |
| --- | --- |
| `proteins/bcl-xl/` | Bcl-xL fragment derived from 1maz with wild-type, F143W, and F143L mutants. |
| `proteins/calmodulin/` | Capped calmodulin EF-hand fragment (`1cll_calmodulin_capped.pdb`). |
| `proteins/chignolin/` | CLN025 variants: 1uao model, pdbfixer-capped form, capped random coil. |
| `proteins/p53CTD/` | p53 C-terminal domain (1dt7) with raw, capped, and extended models. |
| `proteins/proteinG_B1/1pgb.pdb` | Protein G B1 domain (56 aa) for β-sheet benchmarking. |
| `proteins/trypcage/` | 1l2y model 1 and derived conformations for the Trp-cage mini protein. |
| `proteins/ubiquitin/1ubq.pdb` | Ubiquitin X-ray structure (chain A). |
| `proteins/villin_hp35/1yrf.pdb` | Villin headpiece HP35 NMR ensemble. |

### Peptides

| Path | Notes |
| --- | --- |
| `peptides/AAQAA_3/` | Tri-alanine-based helix former with flat, random, helical, and capped (`ac-`/`-nme`) conformations. |
| `peptides/angiotensinII/1n9v.pdb` | Angiotensin II octapeptide (model 1 of NMR ensemble). |
| `peptides/bradykinin/1brk.pdb` | Bradykinin nonapeptide in aqueous solution (model 1). |
| `peptides/dynorphin/` | Dynorphin A peptide from 2n2f (NMR model 1, extended, capped). |
| `peptides/oxytocin/2mgo.pdb` | Oxytocin cyclic nonapeptide model. |
| `peptides/small/` | Small reference peptides (`diala.pdb`, `met_enke.pdb`, `leu_enke.pdb`). |
| `peptides/substanceP/` | Substance P from 2ks9 (original and model 1). |
| `peptides/decaala.pdb` | Deca-alanine reference helix (capped). |
| `peptides/hivcap.pdb` | HIV capsid C-terminal peptide (capped). |
| `peptides/nuccoac.pdb` | Nucleocapsid C-terminal peptide (capped). |

## Usage

1. Pick the system that matches your validation target and copy the desired PDB file.
2. Inspect or visualise using PyMOL, VMD, or `mdtraj`:
   ```python
   import mdtraj as md
   traj = md.load_pdb("proteins/chignolin/2rvd_fixer_capped.pdb")
   print(traj.n_atoms, traj.topology)
   ```
3. Integrate into your MD workflow (e.g., convert to Amber, Gromacs, or OpenMM topologies) and run a short equilibration to confirm your pipeline behaves as expected.

## Contributing

Feel free to open a pull request or issue if you spot an error, want an additional conformation, or have a system that should live alongside the current set.

## Licence

MIT License

Copyright (c) 2025 Shinji Iida

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
