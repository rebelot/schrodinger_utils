from schrodinger.application.desmond.packages import topo, traj, traj_util
from schrodinger.structure import Structure, PDBWriter
from schrodinger.structutils import analyze
import argparse

def main():
    parser = argparse.ArgumentParser(description='Fanculo!')
    parser.add_argument('cms', help='fanculo')
    parser.add_argument('-o', help='fanculo out base name', required=False)
    args = parser.parse_args()

    basename = args.o if args.o else args.cms.split['-out.cms'][0]

    msys, cms, trj = traj_util.read_cms_and_traj(args.cms)

    traj.write_traj(trj, basename + '.xtc', format=traj.Fmt.XTC)

    full_system_indices = analyze.evaluate_asl(cms, 'all')
    sys_st = cms.extract(full_system_indices)
    sys_st.write(basename + '.pdb', format='pdb')

if __name__ == "__main__":
    main()
