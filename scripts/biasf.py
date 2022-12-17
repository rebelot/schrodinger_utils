import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "T", help="system temperature (default: 300)", default=300.0, type=float
    )
    parser.add_argument(
        "-b", help="convert biasfactor to kTemp", metavar="biasf", type=float
    )
    parser.add_argument(
        "-kT", help="convert kTemp to biasfactor", metavar="kTemp", type=float
    )
    parser.add_argument(
        "-u",
        help="units",
        choices=["kcal", "kj"],
        default="kcal",
        metavar="kTemp",
        type=float,
    )
    args = parser.parse_args()

    # Theory:
    # from https://www.plumed.org/doc-v2.5/user-doc/html/belfast-6.html
    # biasf = gamma = (T + dT) / T
    # kTemp = k * dT

    T = args.T
    if args.u == 'kcal':
        k = 0.001987204259  # kcal / (mol * K)
    elif args.u == 'kj':
        k = 0.008314462618 # kJ / (mol * K)
    else:
        raise ValueError('unknown units')

    if args.b:
        gamma = args.b
        dT = gamma * T - T
        print("kTemp =", k * dT)
        return

    if args.kT:
        kdT = args.kT
        gamma = (T + kdT / k) / T
        print("biasf =", gamma)
        return


if __name__ == "__main__":
    main()
