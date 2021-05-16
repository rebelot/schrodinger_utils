import sys

from schrodinger.application.desmond.cms import Cms
from schrodinger.application.desmond.enhsamp import parse_mexpr

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write('usage: enhsamp structure.[mae|dms] mexpression.pot\n')
        exit(1)

    try:
        system = Cms(sys.argv[1])
    except Exception as e:
        sys.stderr.write('failed to read topology %s: %s\n' % (sys.argv[1], e))
        exit(1)

    try:
        mexp = open(sys.argv[2], 'r').read()
    except:
        sys.stderr.write('failed to read potential %s\n' % sys.argv[2])
        exit(1)
    cfg = 'backend.force.term = %s' % parse_mexpr(system, mexp)[10:]
    print(cfg)
