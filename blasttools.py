import sys


class BlastHit:
    def __init__(self, **kwargs):
        if 'blast7_line' in kwargs:
            self._parse_blast7_line(kwargs['blast7_line'])
        else:
            sys.exit('Unrecognized parameter combination to initialize BlastHit instance')

    def _parse_blast7_line(self, line):
        var = line.strip().split('\t')
        self.q_id, self.s_id = var[:2]
        self.identity = float(var[2])
        self.length, self.mismatches, self.gap_opens, \
            self.q_start, self.q_end, self.s_start, self.s_end = map(int, var[3:10])
        self.evalue, self.bit_score = map(float, var[10:])

    def set_q_length(self, thing):
        """
        set_q_length accepts either an int or a dict with a length given the q_id.
        """
        if isinstance(thing, int):
            self.q_length = thing
        elif isinstance(thing, dict):
            if self.q_id in thing:
                if isinstance(thing[self.q_id], int):
                    self.q_length = thing[self.q_id]
                else:
                    sys.exit('Dict entry for %s not an integer: %s' % (self.q_id, thing[self.q_id]))
            else:
                sys.exit('%s not in dict.' % self.q_id)
        else:
            sys.exit('Unexpected object type for length: (%s, %s)' % (thing, type(thing)))

        if self.q_length < 0:
            sys.exit('q_length must be larger than 0')

def parse_blast7_file(fname):
    for line in open(fname):
        if line[0] == '#':
            continue
        yield BlastHit(blast7_line=line)
