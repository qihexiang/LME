# LME: Layered Molecule Editor

LME is a molecule model editor designed with layers.

Plan to implement following functions:

- [x] Manual editable layer
  - [x] functions
  - [x] export
  - [x] load
  - [x] convert
- [x] Readonly base layer
  - [x] export
  - [x] load
  - [x] extractt
  - [x] convert
- [x] Substitute
  - [x] create
  - [x] add to editable layer
- [ ] Helper layers
  - [x] Symmetry layers
    - [x] Inverse layer
      - [x] functions
      - [x] export
      - [x] import
    - [x] Mirror layer
      - [x] functions
      - [x] export
      - [x] import
    - [x] Rotate layer 
      - [x] functions
      - [x] export
      - [x] import
  - [ ] Batch layers
- [ ] Multi-format read/write support
  - [x] mol2(Partial)

## RadiiTable

`libs.RadiiTable` uses file `assets/Radii.csv`, comes from <https://pubs.rsc.org/en/content/articlelanding/2008/DT/b801115j>