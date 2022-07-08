package main

import (
	"bufio"
	"bytes"
	"errors"
	"io"
)

// A struct for one Fasta record
type FastaRecord struct {
	ID          string
	Description string
	Seq         []byte
	Count_A     int
	Count_T     int
	Count_G     int
	Count_C     int
	Score       int64 // this is for e.g., genome completeness
	Idx         int
	encoded     bool
}

// Encode a fasta record, panics if the record is already encoded or if there are
// invalid nucleotides
func (FR *FastaRecord) MustEncode() {
	if FR.encoded {
		panic("Fasta record is already encoded")
	}
	EA := MakeEncodingArray()
	for i, nuc := range FR.Seq {
		if EA[nuc] == 0 {
			panic("invalid nucleotide in file: \"" + string(nuc) + "\"")
		}
		FR.Seq[i] = EA[nuc]
	}
	FR.encoded = true
}

// Decode a fasta record, panics if the record is already decoded
func (FR *FastaRecord) MustDecode() {
	if !FR.encoded {
		panic("Fasta record is already decoded")
	}
	DA := MakeDecodingArray()
	for i, nuc := range FR.Seq {
		FR.Seq[i] = DA[nuc]
	}
	FR.encoded = false
}

var (
	errBadlyFormedFasta = errors.New("Badly formed Fasta")
	errDifferentWidths  = errors.New("Different width sequences in alignment")
)

type Reader struct {
	r *bufio.Reader
}

func NewReader(f io.Reader) *Reader {
	return &Reader{r: bufio.NewReader(f)}
}

// Read reads one fasta record from the underlying reader. The final record is returned with error = nil,
// and the next call to Read() returns an empty FastaRecord struct and error = io.EOF.
func (r *Reader) Read() (FastaRecord, error) {

	var (
		buffer, line, peek []byte
		fields             [][]byte
		err                error
		FR                 FastaRecord
	)

	first := true

	for {

		if first {

			// "ReadBytes reads until the first occurrence of delim in the input,
			// returning a slice containing the data up to and including the delimiter.
			// If ReadBytes encounters an error before finding a delimiter,
			// it returns the data read before the error and the error itself (often io.EOF).
			// ReadBytes returns err != nil if and only if the returned data does not end in delim.
			// For simple uses, a Scanner may be more convenient."
			line, err = r.r.ReadBytes('\n')

			// return even if err == io.EOF, because the file should never end on a fasta header line
			if err != nil {
				return FastaRecord{}, err

				// if the header doesn't start with a > then something is also wrong
			} else if line[0] != '>' {
				return FastaRecord{}, errBadlyFormedFasta
			}

			drop := 0
			// Strip unix or dos newline characters from the header before setting the description.
			if line[len(line)-1] == '\n' {
				drop = 1
				if len(line) > 1 && line[len(line)-2] == '\r' {
					drop = 2
				}
				line = line[:len(line)-drop]
			}

			// split the header on whitespace
			fields = bytes.Fields(line[1:])
			// fasta ID
			FR.ID = string(fields[0])
			// fasta description
			FR.Description = string(line[1:])

			// we are no longer on a header line
			first = false

		} else {

			// peek at the first next byte of the underlying reader, in order
			// to see if we've reached the end of this record (or the file)
			peek, err = r.r.Peek(1)

			// both these cases are fine if first = false, so we can exit the loop and return the fasta record
			if err == io.EOF || peek[0] == '>' {
				err = nil
				break

				// other errors are returned along with an empty fasta record
			} else if err != nil {
				return FastaRecord{}, err
			}

			// If we've got this far, this should be a sequence line.
			// The err from ReadBytes() may be io.EOF if the file ends before a newline character, but this is okay because it will
			// be caught when we peek in the next iteration of the while loop.
			line, err = r.r.ReadBytes('\n')
			if err != nil && err != io.EOF {
				return FastaRecord{}, err
			}

			drop := 0
			// Strip unix or dos newline characters from the sequence before appending it.
			if line[len(line)-1] == '\n' {
				drop = 1
				if len(line) > 1 && line[len(line)-2] == '\r' {
					drop = 2
				}
				line = line[:len(line)-drop]
			}

			buffer = append(buffer, line...)
		}
	}

	FR.Seq = buffer

	return FR, err
}

func LoadAlignment(r io.Reader) ([]FastaRecord, error) {

	records := make([]FastaRecord, 0)
	reader := NewReader(r)

	first := true
	var w int

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		} else if err != nil {
			return []FastaRecord{}, err
		}
		record.MustEncode()

		if first {
			w = len(record.Seq)
			first = false
		} else if len(record.Seq) != w {
			return []FastaRecord{}, errDifferentWidths
		}

		records = append(records, record)
	}

	return records, nil
}

func StreamAlignment(r io.Reader, chnl chan FastaRecord, chnlerr chan error, cdone chan bool) {

	reader := NewReader(r)
	counter := 0

	first := true
	var w int

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		} else if err != nil {
			chnlerr <- err
			return
		}
		record.MustEncode()

		if first {
			w = len(record.Seq)
			first = false
		} else if len(record.Seq) != w {
			chnlerr <- errDifferentWidths
			return
		}

		record.Idx = counter
		counter++

		chnl <- record
	}

	cdone <- true
}

func MakeEncodingArray() [256]byte {
	var byteArray [256]byte

	byteArray['A'] = 136
	byteArray['a'] = 136
	byteArray['G'] = 72
	byteArray['g'] = 72
	byteArray['C'] = 40
	byteArray['c'] = 40
	byteArray['T'] = 24
	byteArray['t'] = 24
	byteArray['R'] = 192
	byteArray['r'] = 192
	byteArray['M'] = 160
	byteArray['m'] = 160
	byteArray['W'] = 144
	byteArray['w'] = 144
	byteArray['S'] = 96
	byteArray['s'] = 96
	byteArray['K'] = 80
	byteArray['k'] = 80
	byteArray['Y'] = 48
	byteArray['y'] = 48
	byteArray['V'] = 224
	byteArray['v'] = 224
	byteArray['H'] = 176
	byteArray['h'] = 176
	byteArray['D'] = 208
	byteArray['d'] = 208
	byteArray['B'] = 112
	byteArray['b'] = 112
	byteArray['N'] = 240
	byteArray['n'] = 240
	byteArray['-'] = 244
	byteArray['?'] = 242

	return byteArray
}

func MakeDecodingArray() [256]byte {
	var byteArray [256]byte

	byteArray[136] = 'A'
	byteArray[72] = 'G'
	byteArray[40] = 'C'
	byteArray[24] = 'T'
	byteArray[192] = 'R'
	byteArray[160] = 'M'
	byteArray[144] = 'W'
	byteArray[96] = 'S'
	byteArray[80] = 'K'
	byteArray[48] = 'Y'
	byteArray[224] = 'V'
	byteArray[176] = 'H'
	byteArray[208] = 'D'
	byteArray[112] = 'B'
	byteArray[240] = 'N'
	byteArray[244] = '-'
	byteArray[242] = '?'

	return byteArray
}
