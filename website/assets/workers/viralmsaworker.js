importScripts("https://cdn.jsdelivr.net/pyodide/v0.23.4/full/pyodide.js");

// constants & global variables
const PATH_TO_PYODIDE_ROOT = "/home/pyodide/";
let downloadResults = false;
let pyodide; 
// holds minimap2 output and is used to block thread until minimap2 is done
let mm2FinishedBuffer; 
// where (in Pyodide FS) to write minimap2 output
let mmiOutput; 
// holds python code of web wrapper for ViralMSA
let ViralMSAWeb;
let REFS; 
let REF_NAMES;

// webworker api
self.onmessage = async (event) => {
    // send data over to main thread to be downloaded to user's computer
    if (event.data.getResults) {
        if (downloadResults) {
            self.postMessage({
                'download': [
                    ['sequence.fas.aln', pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + "output/sequence.fas.aln", { encoding: "utf8" })],
                ]
            })
        } else {
            self.postMessage({
                'error': "No results to download."
            })
        }
    // run ViralMSA
    } else if (event.data.run) {    
        await runViralMSA(event.data.inputSeq, event.data.refSeq, event.data.refID, event.data.omitRef);
    // initialize minimap2 output buffer
    } else if (event.data.arraybuffer) {
        mm2FinishedBuffer = event.data.arraybuffer;
    }
} 

const init = async () => {
    // load pyodide
    pyodide = await loadPyodide({
        stdout: (text) => {
            self.postMessage({'pyodideConsole': "STDOUT: " + text + "\n"})
        },
        stderr: (text) => {
            self.postMessage({'pyodideConsole': "STDERR: " + text + "\n"})
        },
    });

    // load micropip, a package manager for Pyodide
    await pyodide.loadPackage("micropip");
    const micropip = pyodide.pyimport("micropip");
    
    // install biopython, a ViralMSA dependency
    await micropip.install('biopython');

    // create cache directory for ViralMSA sequences and indexes 
    pyodide.FS.mkdir(PATH_TO_PYODIDE_ROOT + 'cache');

    // load in ViralMSA.py
    pyodide.FS.writeFile(PATH_TO_PYODIDE_ROOT + 'ViralMSA.py', await (await fetch("https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py")).text(), { encoding: "utf8" });

    // load in ViralMSAWeb.py
    ViralMSAWeb = await (await fetch("../python/ViralMSAWeb.py")).text()

    // get REFS and REF_NAMES for preloaded reference sequences and indexes
    pyodide.runPython(ViralMSAWeb)
    REFS = pyodide.globals.get('REFS').toJs()
    REF_NAMES = pyodide.globals.get('REF_NAMES').toJs()
    self.postMessage({
        'init': 'done',
        'REFS': REFS, 
        'REF_NAMES': REF_NAMES,
        'VERSION': pyodide.globals.get('VERSION'),
    })
}

init();

// only needs referenceSequence or refID (providing refID means using a preloaded reference sequence and index, which is later fetched)
const runViralMSA = async (inputSequences, referenceSequence, refID, omitRef) => {
    // reset global variable
    downloadResults = false;

    // remove sequence.fas
    if (pyodide.FS.readdir(PATH_TO_PYODIDE_ROOT).includes('sequence.fas')) {
        pyodide.FS.unlink(PATH_TO_PYODIDE_ROOT + 'sequence.fas');
    }

    // remove reference.fas
    if (pyodide.FS.readdir(PATH_TO_PYODIDE_ROOT).includes('reference.fas')) {
        pyodide.FS.unlink(PATH_TO_PYODIDE_ROOT + 'reference.fas');
    }

    // remove output folder
    if (pyodide.FS.readdir(PATH_TO_PYODIDE_ROOT).includes('output')) {
        for (const file of pyodide.FS.readdir(PATH_TO_PYODIDE_ROOT + 'output')) {
            if (file === '.' || file === '..') continue;
            pyodide.FS.unlink(PATH_TO_PYODIDE_ROOT + 'output/' + file)
        }
        pyodide.FS.rmdir(PATH_TO_PYODIDE_ROOT + 'output', true);
    }
    
    // write provided files to Pyodide
    pyodide.FS.writeFile(PATH_TO_PYODIDE_ROOT + 'sequence.fas', inputSequences, { encoding: "utf8" });

	let args = undefined;
    // preloaded reference sequence and index  
    if (refID) {
        // get reference sequence virus name to use in command line args
        refVirus = [...REFS].find(([key, value]) => value === refID)[0];

        // only fetch reference sequence and index if not already in cache
        if (!pyodide.FS.readdir(`${PATH_TO_PYODIDE_ROOT}/cache/`).includes(refID)) {
            // get reference sequence and index
            referenceSequence = await (await fetch("https://raw.githubusercontent.com/niemasd/viralmsa/master/ref_genomes/" + refID + "/" + refID + ".fas")).text();
            refIndex = new Uint8Array(await (await fetch("https://raw.githubusercontent.com/niemasd/viralmsa/master/ref_genomes/" + refID + "/" + refID + ".fas.mmi")).arrayBuffer());
        
            // write reference sequence and index to Pyodide
            pyodide.FS.mkdir(`${PATH_TO_PYODIDE_ROOT}/cache/${refID}`)
            pyodide.FS.writeFile(`${PATH_TO_PYODIDE_ROOT}/cache/${refID}/${refID}.fas`, referenceSequence, { encoding: "utf8" });
            pyodide.FS.writeFile(`${PATH_TO_PYODIDE_ROOT}/cache/${refID}/${refID}.fas.mmi`, refIndex, { encoding: "binary" });
        }

        // set global args variable
        args = `./ViralMSA.py -e email@address.com -s sequence.fas -o output -r ${refVirus} --viralmsa_dir cache`;
    } else {
        pyodide.FS.writeFile(PATH_TO_PYODIDE_ROOT + 'reference.fas', referenceSequence, { encoding: "utf8" });

        // set global args variable
        args = "./ViralMSA.py -e email@address.com -s sequence.fas -o output -r reference.fas --viralmsa_dir cache";
    }

	if (omitRef) {
		args += " --omit_ref";
	}

	pyodide.globals.set("arguments", args);

    
    // set global minimapOverride variable to use BioWASM
    // handles both index building and alignment
    const minimap2Override = (PyProxy) => {
        const command = PyProxy.toJs();
        
        // reset minimap2FinishedBuffer
        Atomics.store(mm2FinishedBuffer, 0, 0);

        // build minimap2 index
        if (command.includes('-d')) {
            // remember mmi output path and change to target.mmi for minimap2 in BioWASM
            mmiOutput = command[command.length - 2];
            command[command.length - 2] =  "target.fas.mmi";
            
            // mount fasta file to BioWASM
            const fastaFile = command[command.length - 1];
            
            // change fasta file name to target.fas for minimap2 in BioWASM
            command[command.length - 1] = "target.fas";

            // post message to main thread to run minimap2
            self.postMessage({
                'runminimap2': 'buildIndex',
                'command': command, 
                'inputSeq': pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + fastaFile, { encoding: "utf8" })
            });

        // alignment
        } else if (command.includes('-a')) {
            // swap out third to last argument (output file)
            command[command.length - 3] = "sequence.fas.sam";

            // swap out second to last argument - use reference.fas instead of reference.fas.mmi and remove pyodide path
            command[command.length - 2] = "target.fas.mmi";

            // swap out last argument - use sequence instead of pyodide path
            command[command.length - 1] = "sequence.fas";

            // run minimap2 in BioWASM
            self.postMessage({
                'runminimap2': 'alignment',
                'command': command, 
                'refSeq': pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + "sequence.fas", { encoding: "utf8" })
            });
        }

        // 300 second timeout
        Atomics.wait(mm2FinishedBuffer, 0, 0, 300000)

        // get minimap2 output and strip trailing null bytes
        const mm2FinishedArray = new Uint8Array(mm2FinishedBuffer.buffer);
        let lastIndex = mm2FinishedArray.length - 1;
        while (lastIndex >= 0 && mm2FinishedArray[lastIndex] === 0) {
            lastIndex--;
        }
        // create new array with no trailing null bytes
        const strippedUint8Array = new Uint8Array(lastIndex + 1);
        strippedUint8Array.set(mm2FinishedArray.subarray(0, lastIndex + 1));

        if (command.includes('-d')) {
            // write mmi file to pyodide
            pyodide.FS.writeFile(mmiOutput, strippedUint8Array, { encoding: "binary" })
        } else if (command.includes('-a')) {
            // write sam file to pyodide
            pyodide.FS.writeFile(PATH_TO_PYODIDE_ROOT + "output/sequence.fas.sam", strippedUint8Array, { encoding: "binary" })
        }
    }
    
    pyodide.globals.set("minimap2Override", minimap2Override);

    // run ViralMSAWeb.py
    pyodide.runPython(ViralMSAWeb);

    // after finished
    downloadResults = true;
    self.postMessage({'finished': true})
}
