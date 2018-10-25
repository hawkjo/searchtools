import re
import itertools


class HMMERTextParser:
    query_line_re = re.compile('^Query:\s+(?P<gene_name>\S+)\s+\[M=(?P<model_len>\d+)\]\s*$')
    def __init__(self, fpath):
        self.sequence_scores = {}
        self.dom_annotations_given_seq_name = {}
        self._parse(fpath)

    def get_best_domain(self):
        if self.sequence_scores == {} and self.dom_annotations_given_seq_name == {}:
            return None
        return max(itertools.chain(*self.dom_annotations_given_seq_name.values()),
                   key=lambda x: x.score)

    def _parse(self, fpath):
        with open(fpath) as f:
            #------------------------------
            # Header
            #------------------------------
            line = next(f)
            while line.startswith('#'):
                line = next(f)
            assert line.strip() == '', (fpath, line)

            #------------------------------
            # Query name and length
            #------------------------------
            line = next(f)
            m = self.query_line_re.match(line)
            self.gene_name = m.group('gene_name')
            self.model_len = m.group('model_len')

            #------------------------------
            # Sequence Scores
            #------------------------------
            line = next(f)
            assert line.startswith('Scores for complete sequences (score includes all domains):'), fpath

            line = next(f)
            assert 'full sequence' in line and 'best 1 domain' in line and '#dom' in line, fpath

            line = next(f)
            assert line.strip().split() == ['E-value', 'score', 'bias', 'E-value', 'score', 'bias',
                                            'exp', 'N', 'Sequence', 'Description'], fpath

            line = next(f)
            assert set(line) == set([' ', '-', '\n']), fpath

            line = next(f)
            if line.strip() == '':
                # No hits found
                line = next(f)
                assert '[No hits detected that satisfy reporting thresholds]' in line, fpath
                return

            while line.strip() != '':
                if 'inclusion threshold' in line:
                    line = next(f)
                hmmer_seq_score = HMMERSequenceScore(line)
                self.sequence_scores[hmmer_seq_score.seq_name] = hmmer_seq_score
                self.dom_annotations_given_seq_name[hmmer_seq_score.seq_name] = []
                line = next(f)

            while line.strip() == '':
                line = next(f)

            #------------------------------
            # Domain Annotation
            #------------------------------
            assert line.startswith('Domain annotation for each sequence (and alignments):'), fpath
            line = next(f)
            assert line.startswith('>>'), fpath
            while line.startswith('>>'):
                # Parse each sequence's domain annotation
                var = line.strip().split()
                assert var[0] == '>>' and len(var) == 2, fpath
                seq_name = var[1]
                assert seq_name in self.dom_annotations_given_seq_name, fpath  # From sequence score parse
                
                line = next(f)
                if not line.strip().startswith('[No individual domains that satisfy reporting thresholds (although complete target did)]'):
                    assert line.strip().split() == ['#', 'score', 'bias', 'c-Evalue', 'i-Evalue',
                                                    'hmmfrom', 'hmm', 'to', 'alifrom', 'ali', 'to',
                                                    'envfrom', 'env', 'to', 'acc'], fpath

                    line = next(f)
                    assert set(line) == set([' ', '-', '\n']), fpath

                    line = next(f)
                    while line.strip() != '':
                        self.dom_annotations_given_seq_name[seq_name].append(
                            HMMERDomainAnnotation(seq_name, line)
                        )
                        line = next(f)

                while line.startswith(' ') or line.strip() == '':
                    line = next(f)

            # Test for correct number of domain annotations
            assert len(self.dom_annotations_given_seq_name[seq_name]) == self.sequence_scores[seq_name].N_dom, fpath

            #------------------------------
            # End of file
            #------------------------------
            assert line.startswith('Internal pipeline statistics summary:'), fpath
            line = next(f)
            assert set(line) <= set([' ', '-', '\n']), fpath

            for s in ['Query model(s):',
                      'Target sequences:',
                      'Passed MSV filter:',
                      'Passed bias filter:',
                      'Passed Vit filter:',
                      'Passed Fwd filter:',
                      'Initial search space (Z):',
                      'Domain search space  (domZ):',
                      '# CPU time:',
                      '# Mc/sec:',
                      '//',
                      '[ok]']:
                line = next(f)
                assert line.startswith(s), fpath

            try:
                line = next(f)
                assert False, 'Should have arrived at end of file %s.' % fpath
            except StopIteration:
                pass


class HMMERSequenceScore:
    def __init__(self, line):
        var = line.strip().split()
        self.seq_evalue, self.seq_score, self.seq_bias = map(float, var[:3])
        self.bestdom_evalue, self.bestdom_score, self.bestdom_bias = map(float, var[3:6])
        self.exp_dom = float(var[6])
        self.N_dom = int(var[7])
        self.seq_name = var[8]
        try:
            self.description = var[9]
        except:
            self.description = ''


class HMMERDomainAnnotation:
    def __init__(self, seq_name, line):
        self.seq_name = seq_name
        var = line.strip().split()
        assert len(var) == 16, seq_name
        self.dom_num = int(var[0])
        self.dom_significance = var[1]
        self.score, self.bias, self.c_evalue, self.i_evalue = map(float, var[2:6])
        self.hmmfrom, self.hmm_to = map(int, var[6:8])
        self.hmm_bounds = var[8]
        self.alifrom, self.ali_to = map(int, var[9:11])
        self.ali_bounds = var[11]
        self.envfrom, self.env_to = map(int, var[12:14])
        self.env_bounds = var[14]
        self.acc = float(var[15])

        # Store pythonic coordinates
        self.hmmstart = self.hmmfrom - 1
        self.hmmend = self.hmm_to
        self.alistart = self.alifrom - 1
        self.aliend = self.ali_to
        self.envstart = self.envfrom - 1
        self.envend = self.env_to

        self.hmm_len = self.hmmend - self.hmmstart
        self.ali_len = self.aliend - self.alistart
        self.env_len = self.envend - self.envstart

    def __str__(self):
        return '\t'.join(map(str, 
                             [self.seq_name,
                              self.dom_num, self.dom_significance,
                              self.score, self.bias, self.c_evalue, self.i_evalue,
                              self.hmmfrom, self.hmm_to, self.hmm_bounds,
                              self.alifrom, self.ali_to, self.ali_bounds,
                              self.envfrom, self.env_to, self.env_bounds,
                              self.acc]))
