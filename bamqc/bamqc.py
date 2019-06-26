"""
Rewrite of rgFastQC.py for BamQC.

SYNOPSIS

    bamqc.py -i input_file -j input_file.name -o output_html_file [-d output_directory]
        [-f fastq|bam|sam] [-n job_name] [-c contaminant_file] [-e fastqc_executable]

"""

import re
import os
import shutil
import subprocess
import optparse
import tempfile
import glob
import gzip
import bz2
import zipfile

class BamQCRunner(object):
    
    def __init__(self,opts=None):
        '''
        Initializes an object to run BamQC in Galaxy. To start the process, use the function run_fastqc()
        '''
        
        # Check whether the options are specified and saves them into the object
        assert opts != None
        self.opts = opts

    def prepare_command_line(self):
        '''
        Develops the Commandline to run BamQC in Galaxy
        '''
        
        # Check whether a given file compression format is valid
        # This prevents uncompression of already uncompressed files
#        infname = self.opts.inputfilename
        infname = self.opts.input ### http://dev.list.galaxyproject.org/BamQC-wrapper-not-seeing-files-at-gzipped-td4666363.html
        linf = infname.lower()
        trimext = False
        # decompression at upload currently does NOT remove this now bogus ending - fastqc will barf
        # patched may 29 2013 until this is fixed properly
        if ( linf.endswith('.gz') or linf.endswith('.gzip') ): 
            f = gzip.open(self.opts.input)
            try:
                f.readline()
            except:
                trimext = True
            f.close()
        elif linf.endswith('bz2'):
            f = bz2.open(self.opts.input,'rb')
            try:
                f.readline()
            except:
                trimext = True
            f.close()
        elif linf.endswith('.zip'):
            if not zipfile.is_zipfile(self.opts.input):
                trimext = True
        if trimext:
	   f = open(self.opts.input)
	   try:
	       f.readline()
	   except:
	       raise Exception("Input file corruption, could not identify the filetype")
           infname = os.path.splitext(infname)[0]
        
        # Replace unwanted or problematic charaters in the input file name
        self.fastqinfilename = re.sub(ur'[^a-zA-Z0-9_\-\.]', '_', os.path.basename(infname))
        # check that the symbolic link gets a proper ending, fastqc seems to ignore the given format otherwise
        if 'fastq' in opts.informat:
            # with fastq the .ext is ignored, but when a format is actually passed it must comply with fastqc's 
            # accepted formats..
            opts.informat = 'fastq'
        elif not self.fastqinfilename.endswith(opts.informat):
            self.fastqinfilename += '.%s' % opts.informat

        # Build the Commandline from the given parameters
        command_line = [opts.executable, '--outdir %s' % opts.outputdir]
        if opts.contaminants != None:
            command_line.append('--contaminants %s' % opts.contaminants)
        if opts.limits != None:
	    command_line.append('--limits %s' % opts.limits)
        command_line.append('--quiet')
        command_line.append('--extract') # to access the output text file
        command_line.append(self.fastqinfilename)
        command_line.append('-f %s' % opts.informat)
        command_line.append('-t ${GALAXY_SLOTS:-4}')
        self.command_line = ' '.join(command_line)

    def copy_output_file_to_dataset(self):
        '''
        Retrieves the output html and text files from the output directory and copies them to the Galaxy output files
        '''
        
        # retrieve html file
        result_file = glob.glob(opts.outputdir + '/*html')
        with open(result_file[0], 'rb') as fsrc:
            with open(self.opts.htmloutput, 'wb') as fdest:
                shutil.copyfileobj(fsrc, fdest)
        
        # retrieve text file
        text_file = glob.glob(opts.outputdir + '/*/fastqc_data.txt')
        with open(text_file[0], 'rb') as fsrc:
            with open(self.opts.textoutput, 'wb') as fdest:
                shutil.copyfileobj(fsrc, fdest)

    def run_fastqc(self):
        '''
        Executes BamQC. Make sure the mandatory import parameters input, inputfilename, outputdir and htmloutput have been specified in the options (opts)
        '''
        
        # Create a log file
        dummy,tlog = tempfile.mkstemp(prefix='rgBamQC',suffix=".log",dir=self.opts.outputdir)
        sout = open(tlog, 'w')
        
        self.prepare_command_line()
        sout.write(self.command_line)
        sout.write('\n')
        sout.write("Creating symlink\n") # between the input (.dat) file and the given input file name
        os.symlink(self.opts.input, self.fastqinfilename)
        sout.write("check_call\n")
        subprocess.check_call(self.command_line, shell=True)
        sout.write("Copying working %s file to %s \n" % (self.fastqinfilename, self.opts.htmloutput))
        self.copy_output_file_to_dataset()
        sout.write("Finished")
        sout.close()

if __name__ == '__main__':
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-j', '--inputfilename', default=None)
    op.add_option('-o', '--htmloutput', default=None)
    op.add_option('-t', '--textoutput', default=None)
    op.add_option('-d', '--outputdir', default="/tmp/shortread")
    op.add_option('-f', '--informat', default='fastq')
    op.add_option('-n', '--namejob', default='rgBamQC')
    op.add_option('-c', '--contaminants', default=None)
    op.add_option('-l', '--limits', default=None)
    op.add_option('-e', '--executable', default='fastqc')
    opts, args = op.parse_args()
    
    assert opts.input != None
    assert opts.inputfilename != None
    assert opts.htmloutput != None
    assert os.path.isfile(opts.executable),'##rgBamQC.py error - cannot find executable %s' % opts.executable
    if not os.path.exists(opts.outputdir): 
        os.makedirs(opts.outputdir)
    
    fastqc_runner = BamQCRunner(opts)
    fastqc_runner.run_fastqc()
