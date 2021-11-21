import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Convert Desmond metadynamics kerseq file into PLUMED readable format\
                     Note that dihedral CVs are split into angle sin and cos and can be\
                     merged using atan2"
    )
    parser.add_argument("kerseq", help="<mtd>.kerseq file")
    parser.add_argument("-kT", help="kTemp value", type=float)
    parser.add_argument("-T", help="system temperature", type=float)
    parser.add_argument("out", help="output filename")
    parser.add_argument("header", help="comma-separated list of CV names")
    args = parser.parse_args()

    CV_names = args.header.split(",")
    OUT_HEADER = f'#! FIELDS time {" ".join(CV_names)} sigma_{" sigma_".join(CV_names)} height biasf\n#! SET multivariate false\n'
    biasf = (args.T + args.kT / 0.001985875) / args.T if args.kT > 0 else -1

    with open(args.kerseq, "r") as fh:
        with open(args.out, "w") as out:
            out.write(OUT_HEADER)

            # Example kerseq header
            #
            # Metadynamics kernel sequence file
            # Lines beginning with a '#' should be ignored
            # This file is suitable for input as an initial kernels file
            #
            #         time           height      center_0         width_0       center_1       width_1

            def skip_header(lines):
                for i, line in enumerate(lines):
                    if line[0] != "#":
                        return i

            start = 0 if not args.kT else 1
            stride = 1 if not args.kT else 2

            lines = fh.readlines()
            for line in lines[skip_header(lines) + start :len(lines):stride]:
                vals = line.strip().split()
                time = vals[0]
                height = vals[1]
                cvs = vals[2::2]
                sigmas = vals[3::2]
                out.write(
                    f"{time} {' '.join(cvs)} {' '.join(sigmas)} {height} {biasf}\n"
                )


if __name__ == "__main__":
    main()
