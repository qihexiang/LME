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
  - [x] Syntax
  - [ ] functions
    - [x] load template
    - [x] add substitute
    - [x] import
    - [x] bond modify
    - [ ] atom modify
    - [x] output

## RadiiTable

`libs.RadiiTable` uses file `assets/Radii.csv`, comes from <https://pubs.rsc.org/en/content/articlelanding/2008/DT/b801115j>

## Workflows

Write workflows in a YAML file.

### Conecpts

### Sections

#### version

```yaml
# Syntax version, default is "1" and now must be "1"
version: "1"
```

#### rootDirectory

```yaml
# Root directory for loading files
# relative (from working directory, where command execute) or absolute
# all relative paths below will be walk from this directory
# default: "./"
rootDirectory: "./"
```

#### substitutes

```yaml  
# list of directories stores substitutes files
# Directory in the left have higher priority
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

#### runners

jobs in workflow are processed with runners. Runners will process each item passthrough the workflow.

A runner function will be built when a job is start with given parameters. The first parameter is `options` comes from `with` section in a job, second parameter is `meta` which contains `version`, `rootDirectory` and `substitutes` specified in global of the yaml.

The runner shuold receive a `StaticLayer` and a `tags`(generated by jobs before) as parameters, and returns a processed `StaticLayer` and a string as tag. Return a list of `(StaticLayer, str)` is also accepted, while the runner should registered with `need_flat` as `True`.

To create a custom runner, you could use Python class like this:

```py3
class AddAHydrogen:
  def __init__(self, options, metas):
    self.options = options
    self.metas = metas
  
  def __call__(self, target, names):
    editable = target.to_editable_layer()
    editable.add_atoms([Atom("H", self.options["position"], "new_hydrogen")])
    return editable.to_static_layer(), "add_hydrogen_at_origin"
```

When use workflow, register your custom runner like this:

```py3
from LME.workflow import Workflow, default_runners
import yaml

runners = default_runners | {
  "add_a_hydrogen": (AddAHydrogen, False)
}

data = yaml.load("some_file", yaml.Loader)
Workflow(data, runners)
```

In the yaml file, call it with:

```yaml
jobs:
# some steps
- name: add_hydrogen_to_0_0_0
  use: add_a_hydrogen
  with:
    position: [0., 0., 0.]
# some steps
```
