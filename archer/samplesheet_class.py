# Objects to read and write Illumina SampleSheets.
# SampleSheetReader(infilename)
# SampleSheetWriter(SampleSheetReader, outfilename, <column seperator (default ',')>, number of seperators)

import collections


class SampleSheetReader(object):
    """
    Read an Illumina SampleSheet and parse in to usable data structure.
    """

    def __init__(self, filename, dual_smoosh=True, sep=','):

        self.header_titles = ("Header", "Reads", "Settings", "Data")
        self.filename = filename

        # Define separator, usually comma.
        self._sep = sep
        # Check to make sure all rows have same number of columns.
        self.sep_count = self._val_count_cols()

        self.header = self._populate_sections(filename, self.header_titles[0], sep) 
        self.runid = self.header['Experiment Name'][0]
        self.date = self.header['Date'][0]
        self.assay = self.header['Assay'][0]
        
        self.reads = self._populate_sections(filename, self.header_titles[1], sep) 
        self.settings = self._populate_sections(filename, self.header_titles[2], sep) 

        # Special case of Data section.
        self.data = self._populate_sections(filename, self.header_titles[3], sep) 
        self._data_header = self._create_data_header()
        self.data = self._reformat_data_section(self.data, self._data_header)
        self.sample_cnt = len(self.data)
        # GET HOP-specific counts!!
        self.hop_cnt = self._get_hop_cnts(proj='QIAseq_V3_HOP')
        self.hop2_cnt = self._get_hop_cnts(proj='QIAseq_V3_HOP2')

        if dual_smoosh:
            self.dual_barcode_smoosh()

        # Exome trio specific
        self.trio = self._get_trio_parents()

        # For Archer workflows
        self.archer_barcodes = self._get_archer_barcodes(proj='QIAseq_V3_HEME2')

    def __iter__(self):
        """
        Define iteration over SampleSheet field dictionaries.
        """

        field_list = (self.header, self.reads, self.settings, self.data)

        for field in field_list:
            yield field

    def _get_trio_parents(self):
        """
        If there are entries in the SampleSheet that have trio info, get it.
        :return:
        """
        parentage = {}
        for samp, info in self.data.items():
            descrip = info['Description']
            if 'mother' in descrip or 'father' in descrip:
                relation = descrip.split(' ')[0]
                proband = descrip.split(' ')[1]
                if proband not in parentage:
                    parentage[proband] = {'proband': proband}
                parentage[proband][relation] = samp

        parentage_out = {}
        for samp, parents in parentage.items():
            parentage_out[parents['proband']] = parents
            parentage_out[parents['mother']] = parents
            parentage_out[parents['father']] = parents

        return parentage_out

    def _get_hop_cnts(self, proj):
        """
        If we send HOP samples alongside other samples, we don't want to count the other samples towards the
        HOP file count.
        :param proj:
        :return:
        """
        count = 0
        for specimen in self.data:
            this_proj = self.data[specimen]['Sample_Project']
            if this_proj == proj:
                count += 1
        return count

    def remove_specimen(self, specimen):
        """
        
        :return: 
        """
        self.data.pop(specimen, None)

    def dual_barcode_smoosh(self):
        """
        Take I7_Index and I5_Index and make in to one barcode.
        For CGD compatibility.
        :return: 
        """
        for specimen in self.data.values():
            if 'I5_Index_ID' in specimen:
                bc_1 = specimen['index']
                bc_1_id = specimen['I7_Index_ID']
                bc_2 = specimen['index2']
                bc_2_id = specimen['I5_Index_ID']
                if specimen['I5_Index_ID']:
                    new_bc = '_'.join([bc_1, bc_2])
                    new_bc_id = '_'.join([bc_1_id, bc_2_id])
                    specimen['I7_Index_ID'] = new_bc_id
                    specimen['index'] = new_bc
                specimen.pop('I5_Index_ID', None)
                specimen.pop('index2', None)

        if 'I5_Index_ID' in self._data_header:
            self._data_header.remove('I5_Index_ID')
            self._data_header.remove('index2')
            self.sep_count -= 2

    def _create_data_header(self):
        """
        Return a list that is the header of the Data section.
        """
        try:
            data_header = self.data['Sample_ID'][0]
        except KeyError:
            raise Exception("Key Sample_ID is not contained within the SampleSheet Data section.")

        return data_header

    def _reformat_data_section(self, data_dict, data_header):
        """
        Break the ordered dictionary up in to a new data structure that facilitates calling individual values from this
        section.
        data_header = [str1, str2, ...]
        data_dict = {'SampleID': [[data]]}
        """
        sample_dict = collections.OrderedDict()
        for key in self.data:
            if key != self._data_header[0]:
                sample_dict[key] = collections.OrderedDict({self._data_header[0]: self.data[key][0][0]})
                for i in range(1, len(data_header)):
                    sample_dict[key][data_header[i]] = data_dict[key][0][i]
        return sample_dict

    def _val_count_cols(self):
        """
        Check to make sure each row contains the same number of columns (defined by sep).
        """
        first = True
        handle = open(self.filename, 'rU')
       
        with handle as samplesheet:
            for line in samplesheet:

                if first:
                    first_count = line.count(self._sep)
                    first = False
                else:
                    if line.count(self._sep) != first_count:
                        raise Exception("SampleSheet malformed, column count differs between rows.")
                    
        return first_count

    @staticmethod
    def _populate_sections(filename, header_title, sep):
        """
        Fill each section of the SampleSheet in to a dictionary.
        """
        
        local_dict = collections.OrderedDict()
        filling = False

        reader = open(filename, 'rU')

        with reader as samplesheet:
            for line in samplesheet:
                if (filling and line.startswith('[') and line.replace(',', '').rstrip('\n')[1:-1]
                        != header_title.title()):
                    filling = False
                elif filling:
                    line = line.rstrip('\n').split(sep)
                    if line[0] not in local_dict:
                        local_dict[line[0]] = [line]
                    else:
                        local_dict[line[0]].append(line)
                elif not filling and line.startswith('[') \
                        and (line.replace(',', '').rstrip('\n')[1:-1] == header_title.title()):
                    filling = True

        return local_dict

    def add_samplesheet_fields(self, specimen, title, value):
        """
        These are the fields that will be added to this particular samplesheet
        record.  Note this currently is meant to add Data section columns.
        :param specimen:
        :param title:
        :param value:
        :return:
        """
        if title not in self.data[specimen]:
            self.data[specimen][title] = value
        if title not in self._data_header:
            self._data_header.append(title)
            self.sep_count += 1

    def _get_archer_barcodes(self, proj='ArcherV6'):
        """
        For Archer workflows, get a list of barcodes, so we can then get reportedvariants from CGD.
        :return:
        """
        barcodes = []
        for info in self.data.values():
            this_proj = info['Sample_Project']
            barcode = info['I7_Index_ID']
            if this_proj == proj:
                barcodes.append(barcode)
        return barcodes


class SampleSheetWriter(object):
    """
    Write an Illumina SampleSheet.
    """

    def __init__(self, samplesheet, filename, cgd, sep=','):

        self.samplesheet = samplesheet
        self.filename = filename
        self.cgd = cgd
        self.sep = sep
        self.write_sheet(filename, samplesheet, sep, samplesheet.sep_count)

    def write_sheet(self, filename, samplesheet, sep, sep_count):
        """
        Write the new SampleSheet from the samplesheet object.
        """

        handle_out = open(filename, 'w')

        i = 0
        no_data_header = True
        for field in samplesheet:
            handle_out.write(self._prepare_section_header(samplesheet.header_titles[i]))
            handle_out.write(sep*sep_count)
            handle_out.write('\n')
            if field != samplesheet.data:
                for entry in field:
                    if len(field[entry]) == 2:
                        for duplicate in field[entry]:
                            handle_out.write(sep.join(duplicate))
                            handle_out.write(self._fill_seps(field[entry][0],
                                                             sep_count))
                            handle_out.write('\n')
                    else:
                        handle_out.write(sep.join(field[entry][0]))
                        handle_out.write(self._fill_seps(field[entry][0],
                                                         sep_count))
                        handle_out.write('\n')
            else:
                for sample in field:

                    if no_data_header:
                        handle_out.write(sep.join([entry for entry in field[sample]]))
                        handle_out.write('\n')
                        no_data_header = False

                    if not no_data_header:
                        if sample in self.cgd:
                            handle_out.write(sep.join([field[sample][entry] for entry in field[sample]]))
                            handle_out.write('\n')

            i += 1

        handle_out.close()

    def _fill_seps(self, entry, sep_count):
        """
        :return:
        """

        if (len(entry)-1) < sep_count:
            num_sep = sep_count - (len(entry)-1)
            return num_sep*self.sep

        return ""

    @staticmethod
    def _prepare_section_header(title):
        """
        Prepare a section title to be written.
        """
        new_title = '[' + title.title() + ']'
        return new_title
