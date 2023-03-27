#! /usr/bin/env python3
import ViralMSA
import asyncio

old_subprocess_call = ViralMSA.subprocess.call

def check_output_override(args):
    return str.encode("'Usage: minimap2'")

def call_override(command, stderr):
    if (command[0] != 'minimap2'):
        return old_subprocess_call(command, stderr=stderr)
    
    minimap2Override(command)

if ('arguments' in globals()):
    ViralMSA.sys.argv = arguments.split()
    ViralMSA.subprocess.check_output = check_output_override
    ViralMSA.subprocess.call = call_override
    ViralMSA.main()