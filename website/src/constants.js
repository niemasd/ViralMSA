export const OUTPUT_ID = "output-console"; 
export const MAX_SHARED_ARRAY_BUFFER_SIZE = 209715200;
export const VIRAL_MSA_LINK = "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py";
export const EXAMPLE_INPUT_FILE = "https://raw.githubusercontent.com/niemasd/viralmsa/master/example/example_hiv.fas";
export const EXAMPLE_PRELOADED_REF = "NC_001802";

export const CLEAR_LOG = () => {
	const textArea = document.getElementById(OUTPUT_ID);
	textArea.value = "";
}

export const LOG = (output, extraFormat = true) => {
	const textArea = document.getElementById(OUTPUT_ID);
	const date = new Date();
	textArea.value += (extraFormat ? `${getTimeWithMilliseconds(date)}: ` : '') + output + (extraFormat ? '\n' : '');
}

export const getTimeWithMilliseconds = date => {
	const t = date.toLocaleTimeString();
	return `${t.substring(0, 7)}.${("00" + date.getMilliseconds()).slice(-3)}`;
}