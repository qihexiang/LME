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
    - [ ] Diff Layer
      - [x] functions(partial)
      - [x] export
      - [ ] import
    - [ ] Bond Generation Layer
      - [x] functions(partial)
      - [x] export
      - [ ] import
  - [ ] ~~Batch layers~~
- [ ] Multi-format read/write support
  - [x] mol2(Partial)
- [ ] Workflow
  - [ ] Syntax
  - [ ] functions
    - [ ] load template
    - [ ] add substitutes and compound
    - [ ] import
    - [ ] bond modify
    - [ ] atom modify
    - [ ] output

## RadiiTable

`libs.RadiiTable` uses file `assets/Radii.csv`, comes from <https://pubs.rsc.org/en/content/articlelanding/2008/DT/b801115j>

## Workflows

Write workflows in a YAML file.

### Sections

#### meta

```yaml
meta:
  # Syntax version, default is "1" and now must be "1"
  version: "1"
  # Root directory for loading files
  # relative (from working directory, where command execute) or absolute
  # all relative paths below will be walk from this directory
  # default: "./"
  rootDirectory: "./"
  # list of directories stores substitutes files
  # Directories listed later have higher priority
  # Will implicitly use subsitutes in built-in library
  # default: "[]"
  substitutes: []
```

#### template

Load a `mol2` file as template

```yaml
# Specify the template mol2 file. Must provide.
template: "./template.mol2"
```

#### jobs

`jobs` is a sequence of workings to do based on the template.

Each item in `jobs` should contains follow properties:

- name(optional): The name of the item, must be unique.
- use: The name of the job
- with: Configuration for the job

This part is still under consideration, a possible syntax will look like this:

```yaml
jobs:
  - name: right_substitute
    use: add_substitutent
    with:
      entry: P11
      centers:
        - H33
        - H36
      substitutes:
        - H
        - CH3
        - CH2F
        - CHF2
        - CF3
    - name: right_substitute
      use: add_subsitutent
      with:
        entry: P12
        centers:
          - H34
          - H35
        substitutes: all
    - name: output_A1
      use: output
      with:
        pattern: "{right_substitute}.{left_subsititute}/A1.mol2"
    - name: add_hydrogen
      use: import
      with:
        targets: 
          - path: hydrogens.mol2
            rotation:
              axis: [0, 0, 1]
              center: [0, 0, 0]
              angle: 120
            translation: no
    - name: connect_Hs_to_cata
      use: modify_bond
      with:
        connect:
          - [HLeft, Mn1, 1]
          - [HRight, N1, 1]
    - name: output_A2
      use: output
      with:
        pattern: "{right_substitute}.{left_subsitute}/A2.mol2"
    - name: add_substrate
      use: import
      with:
        targets:
          - name: acetone
            path: substrate/acetone.mol2
            rotation: no
            translation: no
          - name: amine
            path: substrate/amine.mol2
            rotation: no
            translation: no
    - name: output_A3
      use: output
      with:
        pattern: "{right_substitute}.{left_substitute}/A3_{add_substrate}.mol2"
    - name: generate_A4
      use: modify_bond
      with:
        connect:
          - [HLeft, R1, 1.0]
          - [HRigth, R2, 1.0]
        disconnect:
          - [HLeft, Mn1]
          - [HRight, N2]
    - name: output_A4
      use: output
      with:
        pattern: "{right_substitute}.{left_substitute}/A4_{add_substrate}.mol2"
```
<!-- 
#### add_subsitute (job)

#### import (job)

#### add_atom (job)

#### modify_atom (job)

#### modify_bond (job)

#### output (job) -->