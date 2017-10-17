# OHSU CompBio Tools :bowtie:

A suite of galaxy wrapped tools.

* add_af: 
Add the AF and DP annotations to the SAMPLE column in a VCF.  NOTE: This is an old piece of code that needs to be split up in to its component parts and rewritten.
* add_hotspots:
If you have forced variant calls at specific genomic loci, this will allow you to merge the hotspots VCF with a non-hotspots VCF.  Mainly used to handle collisions, and to apply FILTER column annotations.
* annotate_vcf_with_bed:
Given a basic BED file, annotate any instances of VCF variants overlapping regions defined in said BED file.  FILTER column annotation will be applied.
* breakdancer:
Detect structural variants. (https://github.com/genome/breakdancer) (http://breakdancer.sourceforge.net/)




### From Jaclyn:

This directory contains a pointer to UCSC's pcawg_tools repository. 

To set up this repository correctly, you will need to add a submodule to your repo.

If you have not already forked and cloned this repository, the following command will
allow you to pull the submodule gracefully:

`git clone --recursive git@github.com:OHSUCompBio/docker_tools.git`

If you have already forked and cloned this repository, but have recently updated, the pcawg_tools
directory will be empty. `cd` into pcawg_tools and type `git submodule init`. 

You will see this:

> Submodule 'pcwag_tools' (git@github.com:ucscCancer/pcawg_tools.git)
> registers for path '../pcawg_tools'

At this point, the directory will still be empty. Now perform a `git submodule update`.
You should see pcawg_tools being cloned into that folder, with a final note:

> Submodule path '../pcawg_tools': checked out...

