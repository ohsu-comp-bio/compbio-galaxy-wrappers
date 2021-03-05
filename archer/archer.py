import json
import requests


class ArcherOps:
    """
    This is how we get the data from Archer.

    [{'list_url': 'https://ohsukdldev.analysis.archerdx.com/rest_api/dependency-files/targeted_variant_file',
    'name': 'Targeted Mutations', 'format': 'vcf'},
    {'list_url': 'https://ohsukdldev.analysis.archerdx.com/rest_api/dependency-files/target_region_file',
    'name': 'Panel', 'format': 'gtf'},
    {'list_url': 'https://ohsukdldev.analysis.archerdx.com/rest_api/dependency-files/qc_target_region_file',
    'name': 'Regions of Interest', 'format': 'bed'}]
    """
    def __init__(self, token, url):
        self.token = token
        self.url_prefix = url
        self.header = {'authorization': 'Basic {}'.format(token)}

    def put_dep_file(self, filename, loc='dependency-files', dtype='targeted_region_file'):
        """
        Send a file of certain type to Archer.
        dependency_file_type *
        string
        (formData)
        The type of dependency file being uploaded
        Available values : target_region_file, targeted_variant_file, qc_target_region_file
        custom_file *
        file
        (formData)
        The file that is being added

        curl -X PUT "https://ohsukdldev.analysis.archerdx.com/rest_api/dependency-files/"
        -H "accept: application/json" -H "Content-Type: multipart/form-data" -u $user_email:$user_password
        -d {"dependency_file_type":"targeted_variant_file","ignore_warnings":"false","custom_file":{}}
        :return:
        """
        url = self._create_api_url(loc) + '/'
        payload = self._create_payload(dtype, filename)
        response = requests.put(url, data=payload, headers=self.header)
        print(response)
        return response.content.decode("utf-8")

    def _create_payload(self, dtype, filename, ignore=False):
        """
        {"dependency_file_type":"targeted_variant_file","ignore_warnings":"false","custom_file":{}}
        :return:
        """
        payload = {'dependency_file_type': dtype,
                   'ignore_warnings': str(ignore).lower(),
                   'custom_file': filename}
        return payload

    def _create_api_url(self, *args):
        """
        Take the url prefix and create a new endpoint url based on input.
        :return:
        """
        new_prefix = [self.url_prefix]
        for arg in args:
            new_prefix.append(arg)
        return '/'.join(new_prefix)

    def _get_dep_files(self, loc='dependency-files'):
        """
        Get the
        :return:
        """
        url = self._create_api_url(loc)
        response = requests.get(url, headers=self.header)
        return response.json()['results']

    def _get_targ_muts_loc(self, ftype='Targeted Mutations'):
        """
        Based on ftype, find all of the various targeted mutations files.
        :param ftype:
        :return:
        """
        for entry in self._get_dep_files():
            if entry['name'] == ftype:
                url = entry['list_url']
                response = requests.get(url, headers=self.header)
                return response.json()['results']
        raise FileNotFoundError(f"No {ftype} file found in Archer dependency-files.")

    @staticmethod
    def _get_latest_targ_vers(targ_matches):
        """
        Get the targeted mutations file with the most recent version.
        :return:
        """
        curr_match = 1.0
        best_name = None
        for entry in targ_matches:
            vers = entry.split(' ')[-1].lstrip('v')
            if float(vers) > curr_match:
                curr_match = float(vers)
                best_name = entry
        return best_name

    def _get_all_targ_muts_file(self, targets='Archer Comprehensive Targets'):
        """
        Get all "Targeted Mutations" files that match the $targets string.
        :return:
        """
        targ_list = []
        for entry in self._get_targ_muts_loc():
            if targets in entry['name']:
                targ_list.append(entry['name'])
        return targ_list

    def _create_targ_name(self):
        """
        Create the file name that we will be then grabbing.
        :return:
        """
        targ_list = self._get_all_targ_muts_file()
        return self._get_latest_targ_vers(targ_list)

    def get_targ_vcf(self):
        """
        Return the actual VCF that we want to do further work on.
        :return:
        """
        for entry in self._get_targ_muts_loc():
            if self._create_targ_name() == entry['name']:
                url = entry['file_url']
                response = requests.get(url, headers=self.header)
                return response.content.decode("utf-8"), entry['name']
