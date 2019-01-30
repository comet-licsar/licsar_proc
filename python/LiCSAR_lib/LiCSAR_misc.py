import os
import sys
import subprocess as subp
################################################################################
#CD class
################################################################################
class cd:
    """Context manager for changing the current working directory, used with
    while to provide a temporary path change"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

################################################################################
# Usage exception class
################################################################################
class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg

################################################################################
# Grep function (grep wrapper?)
################################################################################
def grep(arg,file):
    res = subp.check_output(['grep',arg,file])
    return res
