/**
 * @file compress.cpp
 * @brief GZIP compression utilities
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2013
 */

#include "utils/compress.h"
#include "utils/asserter.h"
#include <stdexcept>
#include <cstring>

namespace dominiqs {

static const size_t GZ_BUF_SIZE = 16384;

static void gz_compress(FILE* in, gzFile out)
{
	char buffer[GZ_BUF_SIZE];
	for (;;)
	{
		size_t readcnt = fread(buffer, 1, sizeof(buffer), in);
		if (ferror(in)) throw std::runtime_error("Error reading file!");
		if (readcnt == 0) break;
		int writtencnt = gzwrite(out, buffer, readcnt);
		if (writtencnt != (int)readcnt) throw std::runtime_error("Error writing compressed data to file!");
	}
}

static void gz_uncompress(gzFile in, FILE* out)
{
	char buffer[GZ_BUF_SIZE];
	for (;;)
	{
		int readcnt = gzread(in, buffer, sizeof(buffer));
		if (readcnt < 0) throw std::runtime_error("Error reading compressed file!");
		if (readcnt == 0) break;
		size_t writtencnt = fwrite(buffer, 1, (size_t)readcnt, out);
		if (readcnt != (int)writtencnt) throw std::runtime_error("Error writing to file!");
	}
}

void gzipCompress(const std::string& inPath, const std::string& outPath, const std::string& mode, bool removeOriginal)
{
	FILE* in = fopen(inPath.c_str(), "rb");
	if (!in) throw std::runtime_error("Cannot open file " + inPath);
	gzFile out = gzopen(outPath.c_str(), mode.c_str());
	if (!out) throw std::runtime_error("Cannot open file " + outPath);
	gz_compress(in, out);
	fclose(in);
	gzclose(out);
	if (removeOriginal) unlink(inPath.c_str());
}

void gzipUncompress(const std::string& inPath, const std::string& outPath, bool removeOriginal)
{
	gzFile in = gzopen(inPath.c_str(), "rb");
	if (!in) throw std::runtime_error("Cannot open file " + inPath);
	FILE* out = fopen(outPath.c_str(), "wb");
	if (!out) throw std::runtime_error("Cannot open file " + outPath);
	gz_uncompress(in, out);
	fclose(out);
	gzclose(in);
	if (removeOriginal) unlink(inPath.c_str());
}

/**
 * Internal Class gzstreambuf
 */

gzstreambuf* gzstreambuf::open(const char* name, int open_mode)
{
	if (is_open()) return nullptr;
	mode = open_mode;
	// no append nor read/write mode
	if ((mode & std::ios::ate) || (mode & std::ios::app)) return nullptr;
	if ((mode & std::ios::in) && (mode & std::ios::out)) return nullptr;
	std::string fmode;
	if (mode & std::ios::in) fmode = "rb";
	else if (mode & std::ios::out) fmode = "wb";
	file = gzopen(name, fmode.c_str());
	if (file == nullptr) return nullptr;
	opened = true;
	return this;
}

gzstreambuf* gzstreambuf::close()
{
	if (is_open())
	{
		sync();
		opened = false;
		if (gzclose(file) == Z_OK) return this;
	}
	return nullptr;
}

int gzstreambuf::underflow() // used for input buffer only
{
	if (gptr() && (gptr() < egptr())) return *reinterpret_cast<unsigned char *>(gptr());
	if (!(mode & std::ios::in) || !opened) return EOF;

	// Josuttis' implementation of inbuf
	int n_putback = gptr() - eback();
	if ( n_putback > 4) n_putback = 4;
	std::memcpy(buffer + (4 - n_putback), gptr() - n_putback, n_putback);
	int num = gzread(file, buffer + 4, GZ_BUFFER_SIZE - 4);
	if (num <= 0) return EOF;

	// reset buffer pointers
	setg(buffer + (4 - n_putback), // beginning of putback area
		buffer + 4, // read position
		buffer + 4 + num); // end of buffer

	// return next character
	return *reinterpret_cast<unsigned char *>(gptr());
}

int gzstreambuf::flush_buffer()
{
	// Separate the writing of the buffer from overflow() and sync() operation.
	int w = pptr() - pbase();
	if (gzwrite(file, pbase(), w) != w) return EOF;
	pbump( -w);
	return w;
}

int gzstreambuf::overflow(int c) // used for output buffer only
{
	if (!(mode & std::ios::out) || !opened) return EOF;
	if (c != EOF)
	{
		*pptr() = c;
		pbump(1);
	}
	if (flush_buffer() == EOF) return EOF;
	return c;
}

int gzstreambuf::sync()
{
	// Changed to use flush_buffer() instead of overflow( EOF)
	// which caused improper behavior with std::endl and flush(),
	// bug reported by Vincent Ricard.
	if (pptr() && (pptr() > pbase()))
	{
		if (flush_buffer() == EOF) return -1;
	}
	return 0;
}

} // namespace dominiqs
