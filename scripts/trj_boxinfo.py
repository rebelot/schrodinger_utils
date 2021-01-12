from schrodinger.application.desmond.packages import traj_util
import numpy as np
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_cms', help='path to input cms file')
    parser.add_argument('out_csv', help='path to output csv file')
    args = parser.parse_args()

    in_cms = args.in_cms
    out_csv = args.out_csv
    _, _, trj = traj_util.read_cms_and_traj(in_cms)

    box = []
    for fr in trj:
        box.append([fr.box[0][0], fr.box[1][1], fr.box[2][2]])

    df = pd.DataFrame(np.matrix(box))
    df.to_csv(out_csv, index=False, header=['a','b','c'])

if __name__ == "__main__":
    main()
