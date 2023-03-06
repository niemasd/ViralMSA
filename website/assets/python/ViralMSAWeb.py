#! /usr/bin/env python3
import ViralMSA
import asyncio

oldSubprocessCall = ViralMSA.subprocess.call

def checkOutputOverride(args):
    return str.encode("'Usage: minimap2'")

async def callOverride(command, stderr):
    if (command[0] != 'minimap2'):
        return oldSubprocessCall(command, stderr=stderr)
    
    print('hit')
    await minimap2Override(command)
    print('done')

if ('arguments' in globals()):
    ViralMSA.sys.argv = arguments.split()
    ViralMSA.subprocess.check_output = checkOutputOverride
    ViralMSA.subprocess.call = callOverride
    asyncio.get_running_loop().create_task(ViralMSA.main())