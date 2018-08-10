#!/usr/bin/env python

# DESCRIPTION: Reclassify variants when updated to variant annotations sources occur.
# INPUTS:      Updated HGMD DB
# USAGE: hgmd_ref_pull.py --hgmd <New HGMD DB>
# CODED BY: John Letaw

from __future__ import print_function
from jinja2 import Template
from mysql.connector import errorcode
import argparse
import mysql.connector


VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--chrom', help='Chromosome [1-22,X,Y]')
    parser.add_argument('--coord', help='Genomic Coordinate')
    parser.add_argument('--ref', help='Reference Allele')
    parser.add_argument('--alt', help='Alternate Allele')
    parser.add_argument('--creds', help='Credentials to access HGMD.')
    parser.add_argument('--outfile', help='Output HTML')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def file_type(filename):
    """
    Check for the file type and return the extension.
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type-and-uncompress
    """

    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip"
        }

    max_len = max(len(x) for x in magic_dict)

    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return None


class MysqlRetrieve():
    """
    NEEDS: dictionary {(chrom, pos, ref, alt): [HGMD pmid's]}
    CREATE: Updated version of input.
            {(chrom, pos, ref, alt): True or False whether new and old match}
    """

    def __init__(self, cgd_vars, creds):

        """
        Need from allmut and extrarefs tables.
        """

        self.cgd_vars = cgd_vars
        try:
            self.cnx = mysql.connector.connect(user=creds[0], password=creds[1], host=creds[2], port=creds[3], database=creds[4])
            print("Connection successful")
        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password")
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist")
            else:
                print(err)

        self.new_refs = self._new_refs()
        self.cnx.close()

    def _new_refs(self):
        """
        Get the new list of HGMD references.
        :return:
        """
        pmids = {}
        for variant in self.cgd_vars:
            pmid_list = []
            coord = variant[1]
            chrom = variant[0]
            ref = variant[2]
            alt = variant[3]

            query = "SELECT extrarefs.pmid FROM extrarefs INNER JOIN hgmd_hg19_vcf ON hgmd_hg19_vcf.id=extrarefs.acc_num WHERE hgmd_hg19_vcf.pos = %s AND hgmd_hg19_vcf.chrom = %s AND hgmd_hg19_vcf.ref = %s AND hgmd_hg19_vcf.alt = %s"
            cursor = self.cnx.cursor()
            cursor.execute(query, (coord, chrom, ref, alt))
            for entry in cursor:
                pmid_list.append(str(entry[0]))

            query = "SELECT allmut.pmid FROM allmut INNER JOIN hgmd_hg19_vcf ON hgmd_hg19_vcf.id=allmut.acc_num WHERE hgmd_hg19_vcf.pos = %s AND hgmd_hg19_vcf.chrom = %s AND hgmd_hg19_vcf.ref = %s AND hgmd_hg19_vcf.alt = %s"
            cursor = self.cnx.cursor()
            cursor.execute(query, (coord, chrom, ref, alt))
            for entry in cursor:
                pmid_list.append(str(entry[0]))

            try:
                pmid_list.remove('None')
            except ValueError:
                pass

            pmids[variant] = pmid_list

        cursor.close()
        return pmids


class TemplateHandler():
    def __init__(self, html_tmpl):
        self.j2_env = Template(html_tmpl)

    def render_template(self, context):
        return self.j2_env.render(context)

    def create_index_html(self, output, context):
        with open(output, 'w') as f:
            html = self.render_template(context)
            f.write(html)

def my_html():
    """
    This is the HTML blob that we will fill in.
    :return:
    """
    HTML = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<HTML>
    <HEAD>
        <TITLE>HGMD References Output</TITLE>	
    </HEAD>
    <BODY>
        <h1>HGMD References for {{chrom}}:{{coord}} {{ref}}&gt;{{alt}}</h1>
        <ul>
            {% for pmid in refs %}
            <li><h2><a href="https://www.ncbi.nlm.nih.gov/pubmed/{{pmid}}" target="_blank">pmid:{{pmid}}</a></h2></li>
            {% endfor %}
        </ul>
    </BODY>
</HTML>"""
    return HTML


def build_dict(variants, refs):
    """
    Put together the dictionary that will be passed to the template renderer.
    :return:
    """
    context = {}
    for variant in variants:
        context['chrom'] = variant[0]
        context['coord'] = variant[1]
        context['ref'] = variant[2]
        context['alt'] = variant[3]
        context['refs'] = refs[variant]
    return context


class CredParser():
    def __init__(self, fname, lines=(1,2,3,4,5)):
        self.lines = lines
        self.fname = open(fname, 'rU')
        self.result = []
        for entry in self.find_lines():
            self.result.append(entry)

    def find_lines(self):
        i = 0
        with self.fname as myfile:
            for line in myfile:
                if i in self.lines:
                    yield line
                i += 1


def create_creds(raw_creds):
    """
    Get the credentials packed up in a meaningful way.
    :return:
    """
    return [x.split()[0] for x in raw_creds]


def main():

    args = supply_args()
    variants = [(args.chrom, args.coord, args.ref, args.alt)]
    creds = create_creds(CredParser(args.creds).result)
    my_mysql = MysqlRetrieve(variants, creds)
    context = build_dict(variants, my_mysql.new_refs)
    TemplateHandler(my_html()).create_index_html(args.outfile, context)

if __name__ == "__main__":
    main()
