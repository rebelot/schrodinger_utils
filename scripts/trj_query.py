from schrodinger.application.desmond.packages import traj
import argparse
import math


def main():
    parser = argparse.ArgumentParser(description="Query various MD parameters")
    parser.add_argument('trj', help="Trajectory dir")
    parser.add_argument(
        '-fat', '--frame_at_time', help='Query nearest frame to time TIME (ps)',
        default=-1, type=float, metavar='TIME')
    parser.add_argument(
        '-tof', '--time_of_frame', help='Query the time (ps) of frame FRAME (0-based)',
        default=-1, type=int, metavar='FRAME')
    parser.add_argument(
        '-tt', help="Check timestep consistency", action='store_true')

    args = parser.parse_args()

    trj = traj.read_traj(args.trj)

    f = len(trj)
    t0 = trj[0].time
    te = trj[-1].time
    dt = (te - t0) / f

    print()
    print(f"frames:     {f}")
    print(f"start_time: {t0} ps")
    print(f"last_time:  {te} ps")
    print(f"interval:   {dt:.5f} ps/f", end=" ")

    if args.tt:
        test = sum((trj[i+1].time - trj[i].time) /
                   dt for i in range(f - 1))/f
        print(f"(consistency: {test:.2%})")
    else:
        print()

    # t = t0 + dt * f
    # frames = (te - t0) / dt
    # dt = (te - t0) / frames
    # f = (t - t0) / dt

    if args.frame_at_time >= 0:
        print('Nearest frame to time %.2f ps is %d'
              % (args.frame_at_time, math.floor((args.frame_at_time - t0) / dt)))

    if args.time_of_frame >= 0:
        try:
            print('Chemical time of frame %d is %.2f ps'
                  % (args.time_of_frame, trj[args.time_of_frame].time))
        except IndexError:
            print('IndexError: Index is out of bounds, remember that frame indexing is 0-based')



if __name__ == "__main__":
    main()
