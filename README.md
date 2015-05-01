# OHSU CompBio Tools :bowtie:

A suite of galaxy wrapped tools.

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

You will most likely want to contribute to this module. To do so, treat the pcwag_tools
directory as a repository. Branch and submit pull requests as you wish.

