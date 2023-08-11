import React, { Component, Fragment } from 'react'

import './App.scss'
import loadingCircle from './assets/loading.png'

import {
	LOG,
	CLEAR_LOG,
	MAX_SHARED_ARRAY_BUFFER_SIZE,
	VIRAL_MSA_LINK,
	VIRAL_MSA_REPO_STRUCTURE_LINK
} from './constants.js';

const viralMSAWorker = new Worker(new URL('./assets/workers/viralmsaworker.js', import.meta.url));
const minimap2Worker = new Worker(new URL('./assets/workers/minimap2worker.js', import.meta.url));

export class App extends Component {
	constructor(props) {
		super(props)

		this.state = {
			REFS: undefined,
			REF_NAMES: undefined,

			exampleInput: undefined,
			useExampleInput: false,
			inputFile: undefined,

			refGenomes: new Set(),
			preloadRefOptions: undefined,
			preloadedRef: undefined,

			refFile: undefined,

			omitRef: false,

			startTime: new Date().getTime(),
			timeElapsed: undefined,
			running: false,
			done: false,

			siteTitle: "ViralMSA",
			siteReady: false,
			expandedContainer: undefined,

			sharedArray: undefined,
		}
	}

	async componentDidMount() {
		// Setup WebWorkers
		viralMSAWorker.onmessage = this.handleViralMSAMessage;
		minimap2Worker.onmessage = this.handleMinimap2Message;
		// Setup shared array buffer for waiting for minimap2 to finish and transmitting file data
		const sharedArrayBuffer = new SharedArrayBuffer(MAX_SHARED_ARRAY_BUFFER_SIZE);
		const sharedArray = new Int32Array(sharedArrayBuffer);
		viralMSAWorker.postMessage({ 'arraybuffer': sharedArray })
		minimap2Worker.postMessage({ 'arraybuffer': sharedArray })
		this.setState({ sharedArray })

		// Other initialization
		this.getViralMSAVersion();
		this.fetchPreloadedRef();
		this.initPreloadedRefs();
		this.fetchExampleInput();
	}

	getViralMSAVersion = async () => {
		try {
			const VERSION = [...(await (await fetch(VIRAL_MSA_LINK)).text()).matchAll(/VERSION = '([0-9|.]*)'/gm)][0][1];
			const siteTitle = "ViralMSA v" + VERSION;
			this.setState({ siteTitle });
			document.title = siteTitle;
		} catch (error) {
			console.log(error);
		}
	}

	fetchPreloadedRef = async () => {
		const res = await fetch(VIRAL_MSA_REPO_STRUCTURE_LINK);
		const json = await res.json();
		const refGenomes = new Set();
		for (const file of json.tree) {
			if (file.path.startsWith("ref_genomes/")) {
				refGenomes.add(file.path.split("/")[1]);
			}
		}
		this.setState({ refGenomes });
	}

	initPreloadedRefs = () => {
		const preloadRefInterval = setInterval(() => {
			if (this.state.REFS && this.state.REF_NAMES && this.state.refGenomes.size > 0) {
				clearInterval(preloadRefInterval);
				const preloadRefOptions = [];
				for (const REF_NAME_MAP of this.state.REF_NAMES) {
					const REF_NAME_MAP_TYPE = [...REF_NAME_MAP[1]]
					for (const REF_NAME of REF_NAME_MAP_TYPE) {
						const virus = REF_NAME[0];
						const commonName = REF_NAME[1];
						preloadRefOptions.push(
							<option value={this.state.REFS.get(virus)} key={commonName}>{commonName}</option>
						)
					}
				}

				preloadRefOptions.sort((a, b) => a.key.localeCompare(b.key));
				this.setState({ siteReady: true, preloadRefOptions })
			}
		}, 250)
	}

	fetchExampleInput = async () => {
		this.setState({
			exampleInput: await (await fetch("https://raw.githubusercontent.com/niemasd/viralmsa/master/example/example_hiv.fas")).text()
		})
	}

	handleViralMSAMessage = (event) => {
		if (event.data.error) {
			// error handling
			this.setState({ running: false, done: false, timeElapsed: undefined })
			alert(event.data.error);
		} else if (event.data.init) {
			// done loading pyodide / ViralMSA 
			this.setState({ REFS: event.data.REFS, REF_NAMES: event.data.REF_NAMES })
			LOG("ViralMSA Loaded.\n")
		} else if (event.data.download) {
			// download results
			for (const download of event.data.download) {
				// first element of array is filename, second element is content
				LOG(`Downloading ${download[0]}`)
				this.downloadFile(download[0], download[1])
			}
		} else if (event.data.pyodideConsole) {
			// updating console
			LOG(event.data.pyodideConsole, false)
		} else if (event.data.finished) {
			// on ViralMSA finish
			this.setState({ done: true, timeElapsed: (new Date().getTime() - this.state.startTime) / 1000 })
		} else if (event.data.runminimap2) {
			// Pyodide call to run minimap2 
			if (event.data.runminimap2 === 'alignment') {
				minimap2Worker.postMessage({ 'runminimap2': 'alignment', 'command': event.data.command, 'refSeq': event.data.refSeq });
			} else if (event.data.runminimap2 === 'buildIndex') {
				minimap2Worker.postMessage({ 'runminimap2': 'buildIndex', 'command': event.data.command, 'inputSeq': event.data.inputSeq });
			}
		}
	}

	handleMinimap2Message = (event) => {
		// Minimap2 done running
		if (event.data.minimap2done) {
			const fileData = event.data.minimap2done === 'alignment' ? event.data.sam : event.data.mmi;

			// adjust array size to be divisible by 4
			const adjustedArray = new Uint8Array(Math.ceil(fileData.length / 4) * 4);
			adjustedArray.set(fileData);
			// update shared array buffer
			this.state.sharedArray.set(new Uint32Array(adjustedArray.buffer), 0)

			// notify ViralMSA WebWorker that Minimap2 is done
			Atomics.notify(this.state.sharedArray, 0);
		}
	}

	setInputFile = (event) => {
		this.setState({ useExampleInput: false, inputFile: event.target.files[0] })
	}

	setPreloadedRef = (event) => {
		this.setState({ preloadedRef: event.target.value === 'undefined' ? undefined : event.target.value, })
	}

	setRefFile = (event) => {
		this.setState({ refFile: event.target.files[0] })
	}

	clearRefFile = () => {
		this.setState({ refFile: undefined })
		document.getElementById('ref-sequence').value = null;
	}

	toggleOmitRef = () => {
		this.setState(prevState => ({ omitRef: !prevState.omitRef }))
	}

	toggleExampleData = () => {
		this.setState(prevState => ({ useExampleInput: !prevState.useExampleInput }))
	}

	runViralMSA = async () => {
		// validation
		if (!this.state.useExampleInput && this.state.inputFile === undefined) {
			alert("Please upload an input sequence file.");
			return;
		}

		if (this.state.preloadedRef === undefined && this.state.refFile === undefined) {
			alert("Please upload or select a reference sequence file.");
			return;
		}

		// validation passed
		// clear console and runtime record
		CLEAR_LOG();
		this.setState({ running: true, done: false, timeElapsed: undefined, startTime: new Date().getTime() })

		let inputSeq;
		let refSeq;
		let refIndex;
		let refID;

		// sending file data to webworker
		LOG("Reading input sequence file...")
		if (this.state.useExampleInput) {
			inputSeq = this.state.exampleInput;
		} else {
			const inputSeqReader = new FileReader();
			inputSeqReader.readAsText(this.state.inputFile, "UTF-8");
			inputSeqReader.onload = (e) => {
				inputSeq = e.target.result;
			}
		}

		if (this.state.refFile) {
			const refSeqReader = new FileReader();
			refSeqReader.readAsText(this.state.refFile, "UTF-8");
			refSeqReader.onload = (e) => {
				refSeq = e.target.result;
			}
		} else {
			// only need to provide refID when using a preloaded reference sequence and index
			refID = this.state.preloadedRef;
			// write reference index to minimap2 since it will never be built
			refIndex = new Uint8Array(await (await fetch("https://raw.githubusercontent.com/niemasd/viralmsa/master/ref_genomes/" + refID + "/" + refID + ".fas.mmi")).arrayBuffer());
			minimap2Worker.postMessage({ writeIndex: refIndex })
		}

		// wait until file data is read and then run ViralMSA
		const interval = setInterval(() => {
			if (inputSeq && (refSeq || refID)) {
				clearInterval(interval);
				viralMSAWorker.postMessage({ 'run': 'viralmsa', 'inputSeq': inputSeq, 'refSeq': refSeq, 'refID': refID, 'omitRef': this.state.omitRef });
			}
		}, 100);
	}

	downloadResults = () => {
		viralMSAWorker.postMessage({ 'getResults': 'all' });
	}

	downloadFile = (filename, text) => {
		var a = document.createElement('a');
		a.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
		a.setAttribute('download', filename);
		document.body.appendChild(a);
		a.click();
		document.body.removeChild(a);
	}

	toggleExpandContainer = (container) => {
		this.setState(prevState => {
			return { expandedContainer: prevState.expandedContainer === container ? undefined : container }
		});
	}

	render() {
		return (
			<div className="root">
				<h2 className="mt-5 mb-2 text-center" >{this.state.siteTitle}</h2>
				<p className="text-center my-3">
					WebAssembly implementation of <a href="https://www.github.com/niemasd/ViralMSA" target="_blank"
						rel="noreferrer">ViralMSA</a>.
				</p>
				<div id="loading" className={this.state.siteReady ? 'd-none' : 'mt-4'}>
					<h5 className="text-center me-2">Loading </h5>
					<img className="loading-circle mb-2" src={loadingCircle} alt="loading" />
				</div>
				<div id="content" className={`${this.state.siteReady ? '' : 'd-none'} mt-4 mb-4`}>
					<div id="input" className={`${this.state.expandedContainer === 'input' && 'full-width-container'} ${this.state.expandedContainer === 'output' && 'd-none'}`}>
						<div id="input-header" className="mb-3">
							<h5 className="my-0">Input</h5>
							<h4 className="my-0">
								<i className={`bi bi-${this.state.expandedContainer === 'input' ? 'arrows-angle-contract' : 'arrows-fullscreen'}`} onClick={() => this.toggleExpandContainer('input')}></i>
							</h4>
						</div>
						<div id="ref-seq-container">
							<div id="input-sequences-container" className="mb-3">
								<label htmlFor="input-sequences" className="form-label">Input Sequences (FASTA Format)</label>
								<input className="form-control" type="file" id="input-sequences" onChange={this.setInputFile} />
								{this.state.useExampleInput &&
									<p className="mt-2"><strong>Using Loaded Example Data: <a
										href="https://raw.githubusercontent.com/niemasd/viralmsa/master/example/example_hiv.fas"
										target="_blank" rel="noreferrer">example_hiv.fas</a></strong></p>
								}
							</div>

							<label htmlFor="common-sequences" className="form-label mt-2">
								Select Preloaded Reference Sequence
								{this.state.refFile !== undefined &&
									<span className='mt-2 text-warning'>
										<strong>&nbsp;(Warning: Using Uploaded Reference File)</strong>
									</span>
								}
							</label>
							<select className="form-select" aria-label="Default select example" id="common-sequences" value={this.state.preloadedRef} onChange={this.setPreloadedRef}>
								<option value="undefined">Select a Reference Sequence</option>
								{this.state.preloadRefOptions}
							</select>

							<h5 className="mt-2 text-center">&#8213; OR &#8213;</h5>

							<div>
								<label htmlFor="ref-sequence" className="form-label">Upload Reference Sequence</label>
								<div className="input-group">
									<input className="form-control" type="file" id="ref-sequence" onChange={this.setRefFile} aria-describedby="ref-sequence-addon" />
									<button className="btn btn-outline-danger" type="button" id="ref-sequence-addon" onClick={this.clearRefFile}><i className="bi bi-trash"></i></button>
								</div>
							</div>

							<div className="form-check mt-4">
								<input className="form-check-input" type="checkbox" value="" id="omit-ref" checked={this.state.omitRef} onChange={this.toggleOmitRef} />
								<label className="form-check-label" htmlFor="omit-ref">
									Omit Reference Sequence from Output
								</label>
							</div>
						</div>

						<button type="button" className={`mt-3 w-100 btn ${this.state.useExampleInput ? 'btn-success' : 'btn-warning'}`} onClick={this.toggleExampleData}>
							Load Example Data {this.state.useExampleInput ? '(Currently Using Example Data)' : ''}
						</button>
						<button type="button" className="mt-3 btn btn-primary w-100" onClick={this.runViralMSA}>Run ViralMSA</button>
					</div>
					<div id="output" className={`${this.state.expandedContainer === 'output' && 'full-width-container'} ${this.state.expandedContainer === 'input' && 'd-none'}`}>
						<div id="output-header" className="mb-3">
							<h5 className="my-0">Console</h5>
							<h4 className="my-0">
								<i className={`bi bi-${this.state.expandedContainer === 'output' ? 'arrows-angle-contract' : 'arrows-fullscreen'}`} onClick={() => this.toggleExpandContainer('output')}></i>
							</h4>
						</div>
						<textarea className="form-control" id="output-console" rows="3"></textarea>
						<button type="button" className="mt-4 btn btn-primary w-100" disabled={!this.state.done} onClick={this.downloadResults}>Download Results</button>
						<div id="duration">
							{this.state.timeElapsed &&
								<p id="duration-text" className="my-3">Total runtime: {this.state.timeElapsed} seconds</p>
							}
							{this.state.running && !this.state.done &&
								<Fragment>
									Running ... &nbsp;
									<img id="running-loading-circle" className="loading-circle ms-2" src={loadingCircle}
										alt="loading" />
								</Fragment>
							}
						</div>
					</div>
				</div>
			</div>
		)
	}
}

export default App