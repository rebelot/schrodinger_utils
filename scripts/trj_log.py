#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import time
import re
import argparse


def loglogfile(LOGFILE, TIME):
    fh = open(LOGFILE, 'r')
    LOG = fh.readlines()
    last_line = LOG[-1]
    if re.match('Chemical time:', last_line):
        line = last_line.split()
        CT = float(line[2])/1000
        V = float(line[-1])
        COMP = CT/TIME
        SEC_LEFT = round(24 * 3600 * (TIME - CT) / V)
        DAY = SEC_LEFT / 86400
        HRS = SEC_LEFT / 2600 % 24
        MIN = SEC_LEFT % 3600 / 60
        SEC = SEC_LEFT % 60
        ETA = f'{DAY:02.0f}d:{HRS:02.0f}h:{MIN:02.0f}m:{SEC:02.0f}s'
        print(
            f'Completion: {CT:.2f} of {TIME} ns ({COMP:.2%})  @ {V} ns/day     ETA: {ETA:15s}\r', end='')
    else:
        m = re.findall('Chemical time:.*', ''.join(LOG))
        if len(m) == 0:
            print('Iterations not yet started')
        else:
            print(m[-1])
        return

    time.sleep(1)
    fh.close()


def main():
    parser = argparse.ArgumentParser('Log Desmond trajectories')
    parser.add_argument('logfile', help='input log file')
    parser.add_argument('-f', help='continuous logging', action='store_true')
    args = parser.parse_args()

    LOGFILE = args.logfile
    continuous = args.f

    try:
        fh = open(LOGFILE, 'r')
        LOG = fh.readlines()
    except FileNotFoundError:
        print(LOGFILE, 'does not exist')
        return

    m = re.findall(r'last_time.*=.*\d', ' '.join(LOG))
    if m:
        m = m[0]
    else:
        print('Cannot find last_time keyword')
        return

    TIME = float(m.split()[-1].replace('"', '')) / 1000
    fh.close()

    if continuous:
        while True:
            try:
                loglogfile(LOGFILE, TIME)
            except KeyboardInterrupt:
                break
    else:
        loglogfile(LOGFILE, TIME)

    print()



if __name__ == "__main__":
    main()
