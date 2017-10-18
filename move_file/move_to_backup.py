#!/usr/bin/env python
# Move a source file to a destination file.


import sys
import shutil


def main():

    source_file = sys.argv[1]
    dest_file = sys.argv[2]

    try:
        shutil.move(source_file, dest_file)
        print("Moved " + source_file + " to " + dest_file)
    except:
        print("Move failed.")


if __name__ == "__main__":
    main()


