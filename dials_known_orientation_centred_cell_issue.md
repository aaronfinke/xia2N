# `dials.index known_symmetry.A_matrix` corrupts a centred unit cell

## Summary

Supplying a known orientation to `dials.index` routes indexing through
`IndexerKnownOrientation`, which documents that it expects crystals in the
**primitive setting** and contains a loop to convert them. That loop assigns the
converted crystal to the loop variable and discards it, so the conversion never
happens. `Indexer._apply_symmetry_post_indexing` subsequently converts to the
centred setting a crystal that is already centred, and indexing returns a cell
roughly √2 larger in two directions with the monoclinic angle lost.

There is no error and no warning — indexing reports success, and the wrong cell
is written to `indexed.expt`.

Primitive space groups are unaffected, because there the change of basis is the
identity and discarding it is harmless. Any centred lattice (C, I, F, R) is
affected.

## Version

```
DIALS 3.dev.1508-gb78e3bffe   (dials-v3-30-0 build)
Python 3.13.14, macOS 15 (arm64)
```

## To reproduce

Index any dataset once to get a crystal model, then feed that model's `A`
matrix back in:

```bash
# 1. a normal indexing run, giving a C-centred cell
dials.index imported.expt strong.refl \
    indexing.known_symmetry.space_group="C 1 2 1" \
    indexing.known_symmetry.unit_cell=106.1,95.5,113.7,90,98.1,90
dials.refine indexed.expt indexed.refl

# 2. feed the refined orientation back in as a known orientation
A=$(dials.python -c "
from dxtbx.serialize import load
c = load.experiment_list('refined.expt', check_format=False).crystals()[0]
print(','.join(str(v) for v in c.get_A()))")

dials.index imported.expt strong.refl \
    indexing.known_symmetry.space_group="C 1 2 1" \
    indexing.known_symmetry.A_matrix=$A
```

**Expected** — the cell that was supplied:

```
(105.9, 95.686, 113.944, 90, 98.058, 90)
```

**Observed** — `indexed.expt` from step 2:

```
(183.909, 157.959, 101.751, 90, 89.997, 90)
```

Observed on neutron Laue-TOF data (simulated ESS NMX, inorganic pyrophosphatase,
`C 1 2 1`), but nothing about the failure is specific to the experiment type:
it is a setting conversion that does not happen.

## Cause

`src/dials/algorithms/indexing/known_orientation.py`:

```python
    def find_lattices(self):
        experiments = ExperimentList()
        # Note that the code below cannot be refactored into a single loop
        # as the known orientations may be a shared model between multiple experiments
        for cm in set(self.known_orientations):
            space_group = cm.get_space_group()
            cb_op_to_primitive = (
                space_group.info().change_of_basis_op_to_primitive_setting()
            )
            cm = cm.change_basis(cb_op_to_primitive)      # <-- discarded
        assert len(self.known_orientations) == len(self.experiments)
        for cm, expt in zip(self.known_orientations, self.experiments):
            # indexer expects crystals to be in primitive setting
            experiments.append(
                Experiment(
                    ...
                    crystal=cm,                           # <-- the original, unconverted
                    ...
                )
            )
        return experiments
```

`cm` in the first loop is rebound and then dropped; the second loop iterates
`self.known_orientations` afresh, so the crystals that reach `Experiment` are
the ones that were passed in. `Indexer._apply_symmetry_post_indexing`
(`indexer.py`) then applies the symmetry handler to a crystal that is already in
the centred setting, which is where the enlarged cell comes from.

### This is a regression

Before [#2534](https://github.com/dials/dials/pull/2534) (`a77bd2e77`,
2023-10-25, *"Handle shared models in dials.index with known orientations"*) the
conversion was inside the single loop that built the experiments, so the
converted crystal was the one used:

```python
        for cm, expt in zip(self.known_orientations, self.experiments):
            # indexer expects crystals to be in primitive setting
            space_group = cm.get_space_group()
            cb_op_to_primitive = (
                space_group.info().change_of_basis_op_to_primitive_setting()
            )
            cm = cm.change_basis(cb_op_to_primitive)
            experiments.append(Experiment(..., crystal=cm, ...))
```

Splitting it into two loops to handle shared models dropped the result.

## Suggested fix

Keep the deduplication that #2534 added, but carry the converted crystals
forward — `Crystal` is hashable, so a mapping preserves shared models:

```python
    def find_lattices(self):
        experiments = ExperimentList()
        # Convert once per distinct model, as the known orientations may be a
        # shared model between multiple experiments
        converted = {}
        for cm in set(self.known_orientations):
            cb_op_to_primitive = (
                cm.get_space_group().info().change_of_basis_op_to_primitive_setting()
            )
            converted[cm] = cm.change_basis(cb_op_to_primitive)
        assert len(self.known_orientations) == len(self.experiments)
        for cm, expt in zip(self.known_orientations, self.experiments):
            # indexer expects crystals to be in primitive setting
            experiments.append(
                Experiment(
                    imageset=expt.imageset,
                    beam=expt.beam,
                    detector=expt.detector,
                    goniometer=expt.goniometer,
                    scan=expt.scan,
                    crystal=converted[cm],
                    identifier=expt.identifier,
                )
            )
        return experiments
```

Verified by monkeypatching this `find_lattices` into the run above and passing
the same unconverted `A_matrix`:

```
cell WITH the proposed fix: (105.923, 95.698, 113.961, 90, 98.057, 90)
```

A regression test that seeds `dials.index` with a centred crystal and asserts
the returned cell matches the input would catch this.

## Workaround

Convert in the caller, so the skipped conversion has nothing left to do:

```python
cb_op = crystal.get_space_group().info().change_of_basis_op_to_primitive_setting()
A = crystal.change_basis(cb_op).get_A()
```

With that, the same seeded run returns `(105.846, 95.665, 113.918, 90, 98.039, 90)`
and matches the reference orientation to 0.02°.

Note that this workaround depends on the defect being present: the seed crystal
is built as `Crystal(A_matrix, known_symmetry.space_group)`, so the change of
basis computed in `find_lattices` is the non-identity centred→primitive op
whatever matrix is supplied. Once the conversion is restored, a caller passing a
primitive-setting matrix would have it converted a second time, so any such
workaround has to be removed with the fix.
