import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', help='convert biasfactor to kTemp', metavar='biasf', type=float)
    parser.add_argument('-kT', help='convert kTemp to biasfactor', metavar='kTemp', type=float)
    parser.add_argument('T', help='system temperature', type=float)
    args = parser.parse_args()

    # Theory:
    # from https://www.plumed.org/doc-v2.5/user-doc/html/belfast-6.html
    # biasf = gamma = (T + dT) / T
    # kTemp = k * dT

    T = args.T
    k = 0.001985875  # kcal / (mol * K)

    if args.b:
        gamma = args.b
        dT = gamma * T - T
        print('kTemp =', k * dT)
        return

    if args.kT:
        kdT = args.kT
        gamma = (T + kdT / k) / T
        print('biasf =', gamma)
        return

if __name__ == "__main__":
    main()

