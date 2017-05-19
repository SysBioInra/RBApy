
class FastaEntry:
    def __init__(self, id_, name, set_name, stoichiometry, sequence):
        self.id = id_
        self.name = name
        self.set_name = set_name
        self.stoichiometry = stoichiometry
        self.sequence = sequence

class FastaParser:
    def __init__(self, input_file):
        self._input_file = input_file
        self.entries = []
        try:
            with open(input_file, 'rU') as f:
                self._header = ''
                self._sequence = ''
                for line in f:
                    if line.strip() == '':
                        self._store_current_data()
                        continue
                    if self._read_header(line): continue
                    if self._read_sequence(line): continue
                    print(input_file + ': invalid line\n\t' + line)
                    raise UserWarning('Invalid input file')
                self._store_current_data()
        except IOError:
            print('Could not open file ' + input_file)
            raise UserWarning('Please provide file ' + input_file)

    def _read_header(self, line):
        if line[0] != '>': return False
        self._store_current_data()
        self._header = line[1:].rstrip()
        return True

    def _read_sequence(self, line):
        if self._header:
            self._sequence += line.strip()
        else:
            self._invalid_line(line)
        return True

    def _store_current_data(self):
        if self._header:
            try:
                [rba, id_, name, set_, sto] = self._header.split('|')
                sto = float(sto)
            except ValueError:
                self._invalid_header(self._header)
            if rba != 'rba': self._invalid_header(self._header)
            self.entries.append\
                (FastaEntry(id_, name, set_, sto, self._sequence))
        # reinitialize data
        self._header = ''
        self._sequence = ''

    def _invalid_header(self, line):
        print(self._input_file + ': invalid_header\n\t>' + line)
        print('Expected\n\t>rba|id|name|set_name|stoichiometry')
        raise UserWarning('Invalid input file')

    def _invalid_line(self, line):
        print(self._input_file + ': unexpected data\n\t' + line)
        raise UserWarning('Invalid input file')
