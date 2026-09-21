# Monkeypatch the proposed fix in-process, then index with a NON-primitive A matrix.
from dxtbx.model.experiment_list import Experiment, ExperimentList
from dials.algorithms.indexing import known_orientation


def fixed_find_lattices(self):
    experiments = ExperimentList()
    converted = {}
    for cm in set(self.known_orientations):
        cb_op = cm.get_space_group().info().change_of_basis_op_to_primitive_setting()
        converted[cm] = cm.change_basis(cb_op)
    assert len(self.known_orientations) == len(self.experiments)
    for cm, expt in zip(self.known_orientations, self.experiments):
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


known_orientation.IndexerKnownOrientation.find_lattices = fixed_find_lattices

import libtbx.phil
from dials.array_family import flex
from dials.command_line.index import index as run_index
from dials.command_line.index import phil_scope
from dxtbx.serialize import load

seed = load.experiment_list(
    "../orientation_1/refined.expt", check_format=False
).crystals()[0]
print("Crystal hashable:", isinstance(hash(seed), int))
A = ",".join(str(v) for v in seed.get_A())  # deliberately NOT converted by the caller
user = 'indexing.known_symmetry.space_group="C 1 2 1"\n'
user += f"indexing.known_symmetry.A_matrix={A}\n"
params = phil_scope.fetch(libtbx.phil.parse(user)).extract()
expts = load.experiment_list("imported.expt", check_format=True)
refl = [flex.reflection_table.from_file("strong.refl")]
out_expts, out_refl = run_index(expts, refl, params)
print("cell WITH the proposed fix:", out_expts.crystals()[0].get_unit_cell())
