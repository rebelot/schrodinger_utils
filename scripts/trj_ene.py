import argparse
import sys
import pickle

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
    parser.add_argument("-p", help="plot data")
    parser.add_argument("-s", help="slice trajectory START:END:STEP")
    parser.add_argument(
        "-e",
        action="append",
        help="""Expression of type "NAME:ASL1:[ASL2]:TYPE(s)" where:
            - NAME is an arbitrary name to index the calculation;
            - ASL1 and (optional) ASL2 represent the selections used to calculate energies;
            - TYPE(s) may be a single or comma-separated list of keys representing the energy types to report;
            Avaialble energy types are: Coulomb, vdW, Bond, Angle, Torsion, Total.
            When referring to the Desmond User Manual (2.6.4): Coulomb = pair_elec + nonbonded_elec + far_exclusion, vdW = pair_vdw + nonbonded_vdw;
            although not directly inspectable, improper dihedral energy contribution is summed in the Total energy.
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
        # it could be improved, since also Self energies are calculate and could be retrieved
        # i.e. [GroupEnergy(cms, 'protein', 'Total', 'water'), GroupEnergy(cms, 'protein', 'Total')]
        # will have different keys but data could be retrived from 0 with key = frozenset({0})
        self.key = (group1, group2)
        aids1 = set(cms_model.select_atom(group1))
        aids2 = set(cms_model.select_atom(group2)) if group2 else None
        groups = [aids1, aids2] if aids2 else [aids1]
        # DUCK: required by EnergyFacade
        self.kwargs = {
            "options": self.DEFAULT_OPTIONS,
            "groups": groups,
        }
        # dict key: frozenset({g1, g2}) where gn is the nth group index:
        # results[time][_idx] -> EnergyComponent storing energies between indexed groups
        # when i == j the key is {i}
        self._idx = frozenset([0, group2 and 1 or 0])
        # energy types
        self._attr = [self.TYPE_MAP[t] for t in type]

    def getResult(self, result):
        """
        :type result: tuple[dict[Union(str,frozenset), EnergyComponent]]
        :return: frame-by-frame result of the `type(s)` energies (see TYPE_MAP)
                 of the atom group1 or between group1 and 2
        :rtype: `list[list[float]]`
        """
        return [[getattr(c[self._idx], etype) for etype in self._attr] for c in result]

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
        analyzers.append(GroupEnergy(cms, asl1, types, group2=asl2))
        names.append(name)

    results = analyze(trj, cms, *analyzers, sim_cfg=args.cfg)
    results = [results] if len(analyzers) == 1 else results  # energygroup.py:1122

    # results dims = (n_ana, n_times, n_types(ana))
    for name, ana, res in zip(names, analyzers, results):
        res = np.array(res)  # (n_times, n_types)
        with open(args.out + f"_{name}.dat", "w") as fh:
            fh.write("# " + " ".join(ana._attr) + "\n")
            fh.write("\n".join(" ".join(str(v) for v in r) for r in res))

        if args.p:
            plt.plot([fr.time / 1000 for fr in trj], res.sum(axis=1))
            plt.xlabel("time (ns)")
            plt.ylabel(r"Energy ($kcal\ mol^\{{-1}}$")
    if args.p:
        plt.legend(names)
        plt.savefig(args.out + ".png")


if __name__ == "__main__":
    # cmdline.main_wrapper(main, *sys.argv[1:])
    main()

#   energy_term = r"(angle|dihedral|far_exclusion|far_terms|nonbonded_elec|nonbonded_vdw|pair_elec|pair_vdw|stretch|Total)"
#   time=0.000000 en=7.37335063e+05 E_p=-1.82672213e+05 E_k=0.00000000e+00 E_x=1.03413102e+01 P=-1.30614149e+03 V=7.10988372e+05
#   Dispersion_Correction           (0.000000)      -2.11576729e+03
#   Self_Energy_Correction          (0.000000)      -9.17891509e+05
#   Net_Charge_Correction           (0.000000)      -5.15641662e-12
#   Global_Force_Sum                (0.000000)      0.00000000e+00                                  k = {0,0}           {0,1}           {1,1}       "all"
#    (group)                        (0.000000)            0               1               2         total
#   Kinetic                         (0.000000)      0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
#    (pair)                         (0.000000)         ( 0, 0)         ( 0, 1)         ( 0, 2)         ( 1, 1)         ( 1, 2)         ( 2, 2)      total
# B angle                  *        (0.000000)      9.90058603e+03  1.27680700e+03  1.56821997e+02  0.00000000e+00  1.56821997e+02  0.00000000e+00  1.14910370e+04
# B dihedral               *        (0.000000)      8.16555462e+03  6.89871612e+02  1.65259951e+01  1.44852127e+01  8.26299754e+00  0.00000000e+00  8.89470043e+03
# B improper                        (0.000000)      7.98710917e+01  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  7.98710917e+01
# N pair_elec  14       ,+ *        (0.000000)      -7.77588351e+02 -7.00995515e+02 -3.91073926e+02 3.99505213e+02  0.00000000e+00  0.00000000e+00  -1.47015258e+03
# N pair_vdw   14  elec :  * +,     (0.000000)      3.22434844e+03  3.86199112e+02  3.23889788e+02  -1.77210960e+01 0.00000000e+00  0.00000000e+00  3.91671625e+03
# B stretch             :  *  :     (0.000000)      5.52426398e+03  6.74776070e+02  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  6.19904005e+03
# N far_exclusion       :+ *  : vdw (0.000000)      8.79625472e+05  3.04464804e+03  7.56343060e+02  -6.98762576e+02 -4.11928217e+02 0.00000000e+00  8.82315772e+05
# N nonbonded_elec      '+ *  :     (0.000000)      -1.91022043e+05 2.34648035e+01  -1.85332078e+02 3.65781606e+01  9.06274987e+01  1.60126598e+01  -1.91040692e+05
# N nonbonded_vdw          * +'     (0.000000)      1.71651583e+04  -1.03316076e+03 -1.35890014e+02 -2.70047220e+01 -3.26228006e+01 -1.60229740e+00 1.59348777e+04
# N far_terms              *        (0.000000)      1.96384690e+04  -2.68010199e+04 -1.13082227e+04 9.65330892e+03  8.11512720e+03  1.71623000e+03  1.01389245e+03
#   Total                  *        (0.000000)      7.51524092e+05  -2.24394096e+04 -1.07669379e+04 9.36038911e+03  7.92628867e+03  1.73064036e+03  7.37335063e+05
#   Virial                          (0.000000)    -14055    -1899.2     1726.5    -753.19     -12985      -1200     178.48    -279.21     -13059
#   K.E.tensor                      (0.000000)         0          0          0          0          0          0          0          0          0
#   Pressure_Tensor                 (0.000000)   -1373.4    -185.59     168.71    -73.601    -1268.9    -117.26     17.441    -27.284    -1276.1
