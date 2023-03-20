importScripts("https://cdn.jsdelivr.net/pyodide/v0.22.1/full/pyodide.js");
importScripts("http://localhost:8000/assets/biowasm/aioli.js")

// webworker api
self.onmessage = async (event) => {
    if (event.data.getResults) {
        if (downloadResults) {
            self.postMessage({
                'download': [
                    ['sequence.fas.aln', pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + "output/sequence.fas.aln", { encoding: "utf8" })],
                    ['sequence.fas.sam', pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + "output/sequence.fas.sam", { encoding: "utf8" })],
                ]
            })
        } else {
            self.postMessage({
                'error': "No results to download."
            })
        }
    } else if (event.data.run) {
        if (!pyodide) {
            return;
        }
    
        await runViralMSA(event.data.inputSeq, event.data.refSeq);
    }

} 

// constants & global variables
const PATH_TO_PYODIDE_ROOT = "/home/pyodide/";
let downloadResults = false;
let startTime = new Date().getTime();
let pyodide; 
let CLI;

const init = async () => {

    // load minimap2 from BioWASM
    CLI = await new Aioli([{
        tool: "minimap2",
        version: "2.22",
        urlPrefix: "http://localhost:8000/assets/biowasm",
    }]);
    
    // load pyodide
    pyodide = await loadPyodide({
        stdout: (text) => {
            self.postMessage({'pyodideConsole': "STDOUT: " + text + "\n"})
            console.log(text);
        },
        stderr: (text) => {
            self.postMessage({'pyodideConsole': "STDERR: " + text + "\n"})
            console.err(text);
        },
    });
        
    // load micropip, a package manager for Pyodide
    await pyodide.loadPackage("micropip");
    const micropip = pyodide.pyimport("micropip");
    
    // install biopython, a ViralMSA dependency
    await micropip.install('biopython');
    
    // create cache directory
    pyodide.FS.mkdir(PATH_TO_PYODIDE_ROOT + 'cache');
    
    self.postMessage({'init': 'done'})
}

init();

const runViralMSA = async (inputSequences, referenceSequence) => {
    startTime = new Date().getTime();
    downloadResults = false;

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
    pyodide.FS.writeFile(PATH_TO_PYODIDE_ROOT + 'reference.fas', referenceSequence, { encoding: "utf8" });
    
    // load in ViralMSA.py
    pyodide.FS.writeFile(PATH_TO_PYODIDE_ROOT + 'ViralMSA.py', await (await fetch("../python/ViralMSA.py")).text(), { encoding: "utf8" });

    // set global args variable
    let args = "./ViralMSA.py -e email@address.com -s sequence.fas -o output -r reference.fas --viralmsa_dir cache";
    pyodide.globals.set("arguments", args);
    
    // set global minimapOverride variable to use BioWASM
    const minimap2Override = async (PyProxy) => {
        const command = PyProxy.toJs();
        
        // build minimap2 index
        if (command.includes('-d')) {
            // remember mmi output path and change to target.mmi for minimap2 in BioWASM
            const mmiOutput = command[command.length - 2];
            command[command.length - 2] =  "target.fas.mmi";
            
            // mount fasta file to BioWASM
            const fastaFile = command[command.length - 1];
            await CLI.mount([{
                name: "target.fas",
                data: pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + fastaFile, { encoding: "utf8" }),
            }]);
            
            // change fasta file name to target.fas for minimap2 in BioWASM
            command[command.length - 1] = "target.fas";

            // run minimap2 in BioWASM
            const output = await CLI.exec(command.join(' '));
            
            // copy over mmi file to Pyodide
            const mmi = await CLI.fs.readFile("target.fas.mmi", { encoding: "binary" });
            pyodide.FS.writeFile(mmiOutput, mmi, { encoding: "binary" })

        // alignment
        } else if (command.includes('-a')) {
            // swap out third to last argument (output file)
            const oldOutput = command[command.length - 3];
            command[command.length - 3] = "sequence.fas.sam";

            // swap out second to last argument - use reference.fas instead of reference.fas.mmi and remove pyodide path
            command[command.length - 2] = "target.fas.mmi";

            // mount sequence.fas
            CLI.mount([{
                name: "sequence.fas",
                data: pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + "sequence.fas", { encoding: "utf8" }),
            }]);

            // swap out last argument - use sequence instead of pyodide path
            command[command.length - 1] = "sequence.fas";

            // run minimap2 in BioWASM
            const output = await CLI.exec(command.join(' '));

            // copy over sam file to Pyodide
            const sam = await CLI.fs.readFile("sequence.fas.sam", { encoding: "utf8" });
            pyodide.FS.writeFile(oldOutput, sam, { encoding: "utf8" });
        }
        
    }
    pyodide.globals.set("minimap2Override", minimap2Override);

    // callback for when ViralMSA.py finishes
    const ViralMSAFinish = () => {
        downloadResults = true;
        const outputPreview = pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + "output/sequence.fas.aln", { encoding: "utf8" }).substring(0, 100000)
        const samPreview = pyodide.FS.readFile(PATH_TO_PYODIDE_ROOT + "output/sequence.fas.sam", { encoding: "utf8" }).substring(0, 100000)
        self.postMessage({'outputPreview': outputPreview, 'samPreview': samPreview, 'duration': (new Date().getTime() - startTime) / 1000})
    }
    pyodide.globals.set("ViralMSAFinish", ViralMSAFinish);

    // run ViralMSAWeb.py
    // TODO: Change to actual ViralMSAWeb.py
    pyodide.runPython(await (await fetch("../python/ViralMSAWeb.py")).text());
}
