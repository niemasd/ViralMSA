#! /usr/bin/env python3
import ViralMSA
import asyncio

old_subprocess_call = ViralMSA.subprocess.call

def check_output_override(args):
    return str.encode("'Usage: minimap2'")

async def call_override(command, stderr):
    if (command[0] != 'minimap2'):
        return old_subprocess_call(command, stderr=stderr)
    
    await minimap2Override(command)

async def add_success_callback(task, callback):
    await task
    callback()

if ('arguments' in globals()):
    ViralMSA.sys.argv = arguments.split()
    ViralMSA.subprocess.check_output = check_output_override
    ViralMSA.subprocess.call = call_override
    main_task = asyncio.create_task(ViralMSA.main())
    main_task.add_done_callback(ViralMSAFinish)
    asyncio.get_running_loop().run_until_complete(main_task)