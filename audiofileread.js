const fs = require('fs');
const csv = require('csv-parser');
const { promisify } = require('util');
const struct = require('python-struct');

const NUMBERFORMATS = {
    1: 'int',
    3: 'float',
    0xFFFE: 'extended',
};

const SUBFORMAT_FLOAT = [3, 0, 16, 128, 0, 0, 170, 0, 56, 155, 113];
const SUBFORMAT_INT = [1, 0, 16, 128, 0, 0, 170, 0, 56, 155, 113];

const TYPES_DIRECT = {
    'FLOAT64LE': '<d',
    'FLOAT32LE': '<f',
    'S16LE': '<h',
    'S32LE': '<i',
};

const TYPES_INDIRECT = {
    'S24LE': { 'pattern': 'sssx', 'endian': 'little' },
    'S24LE3': { 'pattern': 'sss', 'endian': 'little' },
};

const SCALEFACTOR = {
    'FLOAT64LE': 1.0,
    'FLOAT32LE': 1.0,
    'S16LE': 2 ** 15,
    'S24LE': 2 ** 23,
    'S24LE3': 2 ** 23,
    'S32LE': 2 ** 31,
};

const BYTESPERSAMPLE = {
    'FLOAT64LE': 8,
    'FLOAT32LE': 4,
    'S16LE': 2,
    'S24LE': 4,
    'S24LE3': 3,
    'S32LE': 4,
};

async function readCoeffs(conf) {
    if (conf.type === 'Raw') {
        const fname = conf.filename;
        const sampleformat = conf.format || 'TEXT';
        let readNbr = conf.read_bytes_lines !== undefined ? conf.read_bytes_lines : null;
        if (readNbr === 0) {
            readNbr = null;
        }
        const skipNbr = conf.skip_bytes_lines || 0;
        const values = await readRawCoeffs(fname, sampleformat, skipNbr, readNbr);
        return values;
    } else if (conf.type === 'Wav') {
        const channel = conf.channel || 0;
        const fname = conf.filename;
        const values = await readWavCoeffs(fname, channel);
        return values;
    }
}

async function readRawCoeffs(filename, sampleformat = null, skipNbr = 0, readNbr = null) {
    if (readNbr === 0) {
        readNbr = null;
    }
    let values;
    if (sampleformat === 'TEXT') {
        values = await readTextCoeffs(filename, skipNbr, readNbr);
    } else {
        if (TYPES_DIRECT.hasOwnProperty(sampleformat)) {
            values = await readBinaryDirectCoeffs(filename, sampleformat, skipNbr, readNbr);
        } else if (TYPES_INDIRECT.hasOwnProperty(sampleformat)) {
            values = await readBinaryIndirectCoeffs(filename, sampleformat, skipNbr, readNbr);
        } else {
            throw new Error(`Unsupported format ${sampleformat}`);
        }
    }
    return values;
}

async function readTextCoeffs(fname, skipLines, readLines) {
    if (readLines === 0) {
        readLines = null;
    }
    const results = [];
    return new Promise((resolve, reject) => {
        fs.createReadStream(fname)
            .pipe(csv())
            .on('data', (row) => {
                if (results.length >= skipLines) {
                    results.push(parseFloat(row[0]));
                }
                if (readLines !== null && results.length >= skipLines + readLines) {
                    this.end();
                }
            })
            .on('end', () => {
                resolve(results);
            })
            .on('error', (error) => {
                reject(error);
            });
    });
}

async function readBinaryDirectCoeffs(fname, sampleformat, skipBytes, readBytes) {
    if (readBytes === null) {
        readBytes = fs.statSync(fname).size;
    }
    const datatype = TYPES_DIRECT[sampleformat];
    const factor = SCALEFACTOR[sampleformat];
    const fileBuffer = await promisify(fs.readFile)(fname);
    const data = fileBuffer.slice(skipBytes, skipBytes + readBytes);
    const values = [];
    for (let i = 0; i < data.length; i += BYTESPERSAMPLE[sampleformat]) {
        values.push(struct.unpack(datatype, data.slice(i, i + BYTESPERSAMPLE[sampleformat]))[0] / factor);
    }
    return values;
}

async function readBinaryIndirectCoeffs(fname, sampleformat, skipBytes, readBytes) {
    if (readBytes === null) {
        readBytes = fs.statSync(fname).size;
    }
    const pattern = TYPES_INDIRECT[sampleformat].pattern;
    const endian = TYPES_INDIRECT[sampleformat].endian;
    const factor = SCALEFACTOR[sampleformat];
    const fileBuffer = await promisify(fs.readFile)(fname);
    const data = fileBuffer.slice(skipBytes, skipBytes + readBytes);
    const values = [];
    for (let i = 0; i < data.length; i += BYTESPERSAMPLE[sampleformat]) {
        const val = struct.unpack(pattern, data.slice(i, i + BYTESPERSAMPLE[sampleformat]));
        values.push(Buffer.concat(val).readIntLE(0, BYTESPERSAMPLE[sampleformat]) / factor);
    }
    return values;
}

async function readWavCoeffs(fname, channel) {
    const params = await readWavHeader(fname);
    if (!params) {
        throw new Error(`Invalid or unsupported wav file '${fname}'`);
    }
    if (channel >= params.channels) {
        throw new Error(`Can't read channel ${channel} from ${fname} which has ${params.channels} channels`);
    }
    const allValues = await readRawCoeffs(
        fname,
        params.sampleformat,
        params.dataoffset,
        params.datalength
    );
    const values = allValues.filter((_, index) => index % params.channels === channel);
    return values;
}

function analyzeWavChunk(type, start, length, file, wavInfo) {
    if (type === "fmt ") {
        const data = file.slice(start, start + length);
        wavInfo.sampleformat = NUMBERFORMATS[struct.unpack("<H", data.slice(0, 2))[0]] || "unknown"; // Uint16
        wavInfo.channels = struct.unpack("<H", data.slice(2, 4))[0];        // Uint16
        wavInfo.samplerate = struct.unpack("<L", data.slice(4, 8))[0];      // Uint32
        wavInfo.byterate = struct.unpack("<L", data.slice(8, 12))[0];       // Uint32
        wavInfo.bytesperframe = struct.unpack("<H", data.slice(12, 14))[0]; // Uint16
        wavInfo.bitspersample = struct.unpack("<H", data.slice(14, 16))[0]; // Uint16
        const bytesPerSample = wavInfo.bytesperframe / wavInfo.channels;

        // Handle extended fmt chunk
        if (wavInfo.sampleformat === "extended") {
            if (length !== 40) {
                console.log("Invalid extended wav header");
                return;
            }
            const cbSize = struct.unpack("<H", data.slice(16, 18))[0];             // Uint16
            const validBitsPerSample = struct.unpack("<H", data.slice(18, 20))[0]; // Uint16
            if (cbSize !== 22 || validBitsPerSample !== wavInfo.bitspersample) {
                console.log("Invalid extended wav header");
                return;
            }
            const _channelMask = struct.unpack("<L", data.slice(20, 24))[0];       // Uint32
            const subformat = struct.unpack("<LHHBBBBBBBB", data.slice(24, 40));
            if (arraysEqual(subformat, SUBFORMAT_FLOAT)) {
                wavInfo.sampleformat = "float";
            } else if (arraysEqual(subformat, SUBFORMAT_INT)) {
                wavInfo.sampleformat = "int";
            } else {
                wavInfo.sampleformat = "unknown";
            }
        }

        if (wavInfo.sampleformat === "int") {
            if (wavInfo.bitspersample === 16) {
                wavInfo.sampleformat = "S16LE";
            } else if (wavInfo.bitspersample === 24 && bytesPerSample === 3) {
                wavInfo.sampleformat = "S24LE3";
            } else if (wavInfo.bitspersample === 24 && bytesPerSample === 4) {
                wavInfo.sampleformat = "S24LE";
            } else if (wavInfo.bitspersample === 32) {
                wavInfo.sampleformat = "S32LE";
            }
        } else if (wavInfo.sampleformat === "float") {
            if (wavInfo.bitspersample === 32) {
                wavInfo.sampleformat = "FLOAT32LE";
            } else if (wavInfo.bitspersample === 64) {
                wavInfo.sampleformat = "FLOAT64LE";
            }
        } else {
            wavInfo.sampleformat = "unknown";
        }
    } else if (type === "data") {
        wavInfo.dataoffset = start + 8;
        wavInfo.datalength = length;
    }
}

function arraysEqual(a, b) {
    if (a.length !== b.length) return false;
    for (let i = 0; i < a.length; i++) {
        if (a[i] !== b[i]) return false;
    }
    return true;
}
async function readWavHeader(filename) {
    try {
        const fileBuffer = await fs.promises.readFile(filename);

        // Read fixed header
        const bufHeader = fileBuffer.slice(0, 12);
        // Verify that the correct identifiers are present
        if (bufHeader.slice(0, 4).toString() !== "RIFF" || bufHeader.slice(8, 12).toString() !== "WAVE") {
            console.log("Input file is not a standard WAV file");
            return;
        }

        const wavInfo = {
            dataoffset: null,
            datalength: null,
            sampleformat: null,
            bitspersample: null,
            channels: null,
            byterate: null,
            samplerate: null,
            bytesperframe: null,
        };

        // Get file length
        const inputFilesize = fileBuffer.length;

        let nextChunkLocation = 12; // skip the fixed header
        while (nextChunkLocation < inputFilesize) {
            const bufChunkHeader = fileBuffer.slice(nextChunkLocation, nextChunkLocation + 8);
            const chunkType = bufChunkHeader.slice(0, 4).toString('utf8');
            const chunkLength = struct.unpack("<L", bufChunkHeader.slice(4, 8))[0];
            analyzeWavChunk(
                chunkType, nextChunkLocation, chunkLength, fileBuffer, wavInfo
            );
            nextChunkLocation += 8 + chunkLength;
        }
        if (wavInfo.datalength !== null && wavInfo.sampleformat !== null && wavInfo.sampleformat !== "unknown") {
            return wavInfo;
        }
    } catch (err) {
        console.log(`Could not open input file: "${filename}", error: ${err}`);
        return;
    }
}
