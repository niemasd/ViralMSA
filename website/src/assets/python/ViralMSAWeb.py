#! /usr/bin/env python3
import ViralMSA

# for retrieving REFS and REF_NAMES from ViralMSA for preloading indexes in web implementation
REFS = ViralMSA.REFS
REF_NAMES = ViralMSA.REF_NAMES
VERSION = ViralMSA.VERSION

old_subprocess_call = ViralMSA.subprocess.call

# override to ensure that ViralMSA sanity check for minimap2 passes for web implementation
def check_output_override(args):
    return str.encode("'Usage: minimap2'")

# override minimap2 call to run Javascript adaptation instead
def call_override(command, stderr):
    if (command[0] != 'minimap2'):
        # run old subprocess call if not minimap2 (won't work because Pyodide doesn't support subprocess.call)
        return old_subprocess_call(command, stderr=stderr)
    
    # run minimap2, from javascript
    minimap2Override(command)

if ('arguments' in globals()):
    # set command line arguments, from javascript
    ViralMSA.sys.argv = arguments.split()
    # set overrides
    ViralMSA.subprocess.check_output = check_output_override
    ViralMSA.subprocess.call = call_override
    # run ViralMSA
    ViralMSA.main()