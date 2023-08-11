export const OUTPUT_ID = "output-console"; 
export const MAX_SHARED_ARRAY_BUFFER_SIZE = 209715200;
export const VIRAL_MSA_LINK = "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py";
export const VIRAL_MSA_REPO_STRUCTURE_LINK = "https://api.github.com/repos/niemasd/viralmsa/git/trees/master?recursive=1";

export const CLEAR_LOG = () => {
	const textArea = document.getElementById(OUTPUT_ID);
	textArea.value = "";
}

export const LOG = (output, includeTime = true) => {
	const textArea = document.getElementById(OUTPUT_ID);
	const date = new Date();
	textArea.value += (includeTime ? `${getTimeWithMilliseconds(date)}: ` : '') + output;
}

export const getTimeWithMilliseconds = date => {
	const t = date.toLocaleTimeString();
	return `${t.substring(0, 7)}.${("00" + date.getMilliseconds()).slice(-3)}`;
}