class ArcherVrnts:
    def __init__(self, archer):
        self.archer = archer
        self.header = self._get_header()
        self.vrnts = self._collect_vrnts()

    def _get_header(self):
        """
        Get the header, so we can rewrite it.
        :return:
        """
        header = []
        for entry in self.archer:
            if entry.startswith('#'):
                header.append(entry)
        return header

    def _collect_vrnts(self):
        """
        Get all of the variants together so we can determine what needs to be added.
        :return:
        """
        vrnts = {}
        for entry in self.archer.split('\n'):
            if not entry.startswith('#'):
                entry = entry.rstrip('\n').split()
                if entry:
                    uniq_key = (entry[0], entry[1], entry[3], entry[4])
                    vrnts[uniq_key] = entry
        return vrnts


class ArcherWriteVrnts:
    def __init__(self, filename):
        self.filename = filename

    def write_archer(self, archer):
        with open(self.filename, 'w') as outfile:
            outfile.write(archer)
