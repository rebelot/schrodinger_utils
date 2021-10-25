import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
from schrodinger.application.desmond.packages import analysis, topo, traj, traj_util
from schrodinger.application.desmond.packages.energygroup import (
    EnergyGroupBase,
    # EnergyFacade,
    # SliceParams
    analyze,
)
# from schrodinger.job import launchapi
# from schrodinger.utils import cmdline


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Calculate energy properties from MD simulation, example:
            trj_ene.py md-out.cms md-out.cfg enegrp -e 'interprot:c.n A:c.n B:Coulomb,vdW' -e solv:solvent::Total -e solvprot:solvent:protein:Coulomb,vdW -e prot:protein::Total"""
    )
    parser.add_argument("cms", help="input cms file (-out.cfg)")
    parser.add_argument("cfg", help="simulation -out.cfg")
    parser.add_argument("out", help="output file(s) basename")
    parser.add_argument("-t", help="trajectory dir")
    parser.add_argument("-s", help="slice trajectory START:END:STEP")
    parser.add_argument(
        "-e",
        action="append",
        help="""Expression of type "NAME:ASL1:[ASL2]:TYPE(s)" where:
            - NAME is an arbitrary name to index the calculation;
            - ASL1 and (optional) ASL2 represent the selections used to calculate energies;
            - TYPE(s) may be a single or comma-separated list of keys representing the energy types to report;
            Avaialble energy types are: Coulomb, vdW, Bond, Angle, Torsion, Total.
            TYPE energies are reported within ASL1 or between ASL1 and ASL2 if the latter is specified.""",
    )
    return parser.parse_args()


# def get_job_spec_from_args(argv):
#     # first argument is this script
#     args = parse_args(argv[1:])
#     job_builder = launchapi.JobSpecificationArgsBuilder(argv, use_jobname_log=True)
#     job_builder.setInputDirectory(args.t)
#     job_builder.setInputFile(args.cms)
#     job_builder.setInputFile(args.cfg)
#     job_builder.setOutputFile(args.out + ".png")
#     for expr in args.e:
#         name = expr.split(":")[0]
#         job_builder.setOutputFile(args.out + f"_{name}.dat")
#     job_builder.setJobname("enecalc")
#     job_builder.setInputDirectory
#     return job_builder.getJobSpec()


class GroupEnergy(EnergyGroupBase):
    """
    Analyzer for energy group calculation with multiple molecule selection.
    """

    TYPE_MAP = {
        "Coulomb": "elec",
        "vdW": "vdw",
        "Bond": "stretch",
        "Angle": "angle",
        "Torsion": "dihedral",
        "Total": "Total",
    }

    def __init__(self, cms_model, group1, type, group2=None):
        """
        :type        cms_model: `Cms`
        :param          group1: index of the atoms to decompose energies
        :type  molecule_number: `List[int]`
        :param            type: pre-defined energy type(s), see TYPE_MAP
        :type             type: `List[str]`
        :param          group2: if provided, will return the `[type(s)]` energies between group1 and group2
        :type            type2: `List[int] or None`
        """
        # DUCK: this is used to cache the results, since all energy types are actually calculated
        self.key = (group1, group2)
        aids1 = set(cms_model.select_atom(group1))
        aids2 = set(cms_model.select_atom(group2)) if group2 else None
        groups = [aids1, aids2] if aids2 else [aids1]
        # DUCK: required by EnergyFacade
        self.kwargs = {
            "options": self.DEFAULT_OPTIONS,
            "groups": groups,
        }
        # dict key: (g1, g2) where gn is the nth group index:
        # results[time][_idx] -> EnergyComponent storing energies between indexed groups
        self._idx = frozenset([0, 1]) if group2 else frozenset([0])
        self._attr = [self.TYPE_MAP[t] for t in type]

    def getResult(self, result):
        """
        :return: frame-by-frame result of the `type(s)` energies (see TYPE_MAP)
                 of the atom group1 or between group1 and 2
        :rtype: `List[List(float)]`
        """
        return [[getattr(c[self._idx], attr) for attr in self._attr] for c in result]

def main():
    args = parse_args()

    if not args.t:
        _, cms, trj = traj_util.read_cms_and_traj(args.cms)
    else:
        _, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.t)

    slicer = (
        slice(*[int(v) if v else None for v in args.s.split(":")])
        if args.s
        else slice(None, None)
    )
    trj = [fr.copy() for fr in trj[slicer]]

    analyzers = []
    names = []
    for expr in args.e:
        name, asl1, asl2, types = expr.split(":")
        types = types.replace(" ", "").split(",")
        analyzers.append(GroupEnergy(cms, asl1, types, group2=asl2 or None))
        names.append(name)

    results = analyze(trj, cms, *analyzers, sim_cfg=args.cfg)
    results = [results]  # energygroup.py:1122 yeah... why?!

    # results dims = (n_ana, n_times , n_types)
    for name, ana, res in zip(names, analyzers, results):
        res = np.array(res)  # (n_times, n_types)
        with open(args.out + f"_{name}.dat", "w") as fh:
            fh.write("# " + " ".join(ana._attr) + "\n")
            fh.write("\n".join(" ".join(str(v) for v in r) for r in res))

        plt.plot([fr.time / 1000 for fr in trj], res.sum(axis=1))
        plt.xlabel("time (ns)")
        plt.ylabel(r"Energy ($kcal\ mol^\{{-1}}$")
    plt.legend(names)
    plt.savefig(args.out + ".png")


if __name__ == "__main__":
    # cmdline.main_wrapper(main, *sys.argv[1:])
    main()
