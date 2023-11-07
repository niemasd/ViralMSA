let CLI;
let mm2FinishedBuffer;

const init = async () => {
	// load minimap2 from BioWASM
	const Aioli = (await import("@biowasm/aioli/dist/aioli")).default;
	CLI = await new Aioli([{
		tool: "minimap2",
		version: "2.22",
		urlPrefix: `${self.location.origin}${import.meta.env.BASE_URL || ''}tools/minimap2`,
	}]);
}

init();

// listen for messages from main thread
self.onmessage = (event) => {
	if (event.data.arraybuffer) {
		mm2FinishedBuffer = event.data.arraybuffer;
	} else if (event.data.writeIndex) {
		CLI.fs.writeFile("target.fas.mmi", event.data.writeIndex, { encoding: "binary" });
	} else if (event.data.runminimap2) {
		runMinimap2(event.data.command, event.data.inputSeq, event.data.refSeq);
	}
}

// run minimap2 with provided command, sequences
const runMinimap2 = async (command, inputSeq, refSeq) => {
	// reset minimap2 output buffer
	mm2FinishedBuffer.fill(0);

	// build minimap2 index
	if (command.includes('-d')) {
		await CLI.mount([{
			name: "target.fas",
			data: inputSeq,
		}]);

		// run minimap2 in BioWASM
		await CLI.exec(command.join(' '));
		
		// send over output file data (minimap2 index)
		self.postMessage({ 'minimap2done': 'buildIndex', 'mmi': await CLI.fs.readFile("target.fas.mmi", { encoding: "binary" }) })
		
		// alignment
	} else if (command.includes('-a')) {
		// mount sequence.fas
		await CLI.mount([{
			name: "sequence.fas",
			data: refSeq,
		}]);
		
		// run minimap2 in BioWASM
		const output = await CLI.exec(command.join(' '));
		self.postMessage({ 'log': output })

		// send over output file data (sequence alignment / map file)
		self.postMessage({ 'minimap2done': 'alignment', 'sam': await CLI.fs.readFile("sequence.fas.sam", { encoding: "binary" }) })
	}
}
