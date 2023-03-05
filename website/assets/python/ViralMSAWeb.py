#! /usr/bin/env python3
import ViralMSA

def checkOutputOverride(args):
    return str.encode("'Usage: minimap2'")

if ('arguments' in globals()):
    ViralMSA.sys.argv = arguments.split()
    ViralMSA.subprocess.check_output = checkOutputOverride
    ViralMSA.main()