import requests
from requests_oauthlib import OAuth1


class CgdOps:
    """
    Manage the construction of CGD endpoints in order to grab reportedvariants.
    Typical reportedvariants endpoint url:
    https://kdlwebprod02.ohsu.edu/cgd/service/illumina/variants/samplerun/201020_NB551673_0011_AHM3W3BGXG/N702-A_S505-A
    Just need the runid/barcode ids.
    """
    def __init__(self, token, url='https://kdlwebprod02.ohsu.edu/cgd/service/illumina/variants/samplerun',
                 client_key='cgd-key'):
        self.url_prefix = url
        self.auth = OAuth1(client_key=client_key, client_secret=token, signature_type='auth_header')

    def get_url(self, terms):
        url = self._create_url(terms)
        r = requests.get(url, auth=self.auth)
        return r.json()

    def _create_url(self, terms):
        new_url = [self.url_prefix]
        new_url.extend(terms)
        return '/'.join(new_url)


class CgdPropertiesParser:
    def __init__(self, filename):
        self.filename = filename
        self.prop_vals = self._get_values()

    def _get_values(self, sep='='):
        """
        Get all the values in a simple config file.
        Comments will start with a #.
        Other values should just be key/val separated by =.
        :return:
        """
        vals = {}
        with open(self.filename, 'r') as myfile:
            for line in myfile:
                if not line.startswith('#'):
                    if sep in line:
                        line = line.rstrip('\n')
                        key = line.split(sep)[0]
                        val = line.split(sep)[1]
                        vals[key] = val
        return vals
