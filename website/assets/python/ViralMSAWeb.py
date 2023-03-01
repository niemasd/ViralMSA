#! /usr/bin/env python3
from sys import argv
import ViralMSA

if ('arguments' in globals()):
    ViralMSA.sys.argv = arguments.split()
    ViralMSA.main()