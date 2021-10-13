import gzip

def file_type(filename):
    """
    Check for the file type and return the extension.
    http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type-and-uncompress
    """

    magic_dict = {
        b"\x1f\x8b\x08": "gz",
        b"\x42\x5a\x68": "bz2",
        b"\x50\x4b\x03\x04": "zip"
        }

    max_len = max(len(x) for x in magic_dict)

    with open(filename, 'rb') as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return None

def handle_return(filename):
    """
    Based off of which file type we're dealing with, return the correct file handle object.
    """
    file_t = file_type(filename)
    print(file_t)
    if file_t == 'gz':
        return gzip.open(filename, 'rt')
    elif not file_t:
        return open(filename, 'rU')
    else:
        raise Exception("Compressed file type " + file_t + " not implemented.")
